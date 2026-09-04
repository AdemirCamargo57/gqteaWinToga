"""Tests for the Bond Angle Analysis tool (``bondAngle.py``).

Run from ``venv/src/``::

    python -m pytest tests/ -q

Mirrors tests/test_bond.py: the Toga GUI is never instantiated, and these drive
``BondAngleAnalyser`` (the logic class) rather than ``BondAngleUI``. Angles are
checked against values that can be worked out by hand.
"""
import os
import glob
import math
import asyncio
from statistics import variance

import numpy as np
import pytest

import bondAngle
from bondAngle import BondAngleAnalyser


class FakeDialog:
    def __init__(self, title, message):
        self.title, self.message = title, message


class FakeWindow:
    def __init__(self):
        self.dialogs = []

    async def dialog(self, dlg):
        self.dialogs.append((dlg.title, dlg.message))
        return None

    def messages(self):
        return " ".join(msg for _, msg in self.dialogs)


class FakeInput:
    def __init__(self, value=""):
        self.value = value


class FakeSwitch:
    def __init__(self, value=False):
        self.value = value


@pytest.fixture(autouse=True)
def _patch_toga(monkeypatch):
    monkeypatch.setattr(bondAngle.toga, "InfoDialog", FakeDialog)


@pytest.fixture(autouse=True)
def _isolate_plot_state():
    from displayPlots import DisplayPlots
    DisplayPlots.saved_plot_files = []
    DisplayPlots.saved_plot_data = []
    yield
    DisplayPlots.saved_plot_files = []
    DisplayPlots.saved_plot_data = []


# --------------------------------------------------------------------------- #
# Fixture builders                                                              #
# --------------------------------------------------------------------------- #
def write_xyz(path, frames, symbols, trailer=""):
    with open(path, "w") as f:
        for coords in frames:
            f.write(f"{len(symbols)}\ncomment\n")
            for s, (x, y, z) in zip(symbols, coords):
                f.write(f"{s} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write(trailer)


def triplet_traj(tmp_path, frames, symbols=("O", "C", "N"), trailer="", name="traj.xyz"):
    """Each entry of `frames` is a 3x3 array of i, j (vertex), k positions."""
    traj = tmp_path / name
    write_xyz(traj, [np.asarray(fr, float) for fr in frames], list(symbols),
              trailer=trailer)
    return traj


def right_angle_frame(scale=1.0):
    """i-j-k with a 90 degree vertex at j."""
    return [[scale, 0, 0], [0, 0, 0], [0, scale, 0]]


def angle_frame(degrees_wanted):
    """i-j-k whose vertex angle is exactly `degrees_wanted`."""
    rad = math.radians(degrees_wanted)
    return [[1, 0, 0], [0, 0, 0], [math.cos(rad), math.sin(rad), 0]]


def make_analyser(out_dir, traj, num_atoms, n_frames, *, max_angle="180",
                  time_step="5.0", sampling="10", temperature="300",
                  atom_labels="1 2 3", bin_width="1.0", cell="",
                  jacobian=False, smooth_fe=False, stub_plots=True):
    # No progress bar and no progress label: this tool reports status as text.
    a = BondAngleAnalyser()
    a.main_window = FakeWindow()
    a.multi_line_text = FakeInput()
    a.output_dir = str(out_dir)
    if stub_plots:
        a.save_plots = lambda *args, **kw: None
    a.display_plots = lambda *args, **kw: None
    a.textInput_max_angle = FakeInput(max_angle)
    a.textInput_time_step = FakeInput(time_step)
    a.textInput_sampling_interval = FakeInput(sampling)
    a.textInput_temperature = FakeInput(temperature)
    a.textInput_atom_labels = FakeInput(atom_labels)
    a.textInput_bin_width = FakeInput(bin_width)
    a.textInput_cell_lengths = FakeInput(cell)
    a.switch_show_plots = FakeSwitch(True)
    a.switch_save_csv = FakeSwitch(False)
    a.switch_use_jacobian = FakeSwitch(jacobian)
    a.switch_smooth_fe = FakeSwitch(smooth_fe)
    a.trajec = str(traj)
    a.num_atoms = num_atoms
    a.total_frame_number = n_frames
    return a


def run_analysis(a):
    a._reset_analysis_state()
    if not asyncio.run(a.read_params(None)):
        return False
    if not asyncio.run(a.bond_angle()):
        return False
    if not asyncio.run(a.distribution_function()):
        return False
    return asyncio.run(a.free_energy())


def angle_values(a):
    return [v for _, v in a.angles]


# --------------------------------------------------------------------------- #
# Core geometry                                                                 #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("wanted", [30.0, 60.0, 90.0, 120.0, 179.0])
def test_angle_matches_hand_computed_value(tmp_path, wanted):
    traj = triplet_traj(tmp_path, [angle_frame(wanted)])
    a = make_analyser(tmp_path, traj, 3, 1)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    assert angle_values(a) == pytest.approx([wanted])


def test_angle_is_independent_of_arm_length(tmp_path):
    """Only the direction of the two arms matters, not their length."""
    traj = triplet_traj(tmp_path, [right_angle_frame(1.0), right_angle_frame(7.5)])
    a = make_analyser(tmp_path, traj, 3, 2)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    assert angle_values(a) == pytest.approx([90.0, 90.0])


def test_time_axis_uses_timestep_times_sampling(tmp_path):
    traj = triplet_traj(tmp_path, [right_angle_frame(), right_angle_frame()])
    a = make_analyser(tmp_path, traj, 3, 2, time_step="5.0", sampling="10")
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    step_ps = 5.0 * 10 * BondAngleAnalyser.atufs / 1000.0
    assert [t for t, _ in a.angles] == pytest.approx([0.0, step_ps])


def test_degenerate_geometry_is_skipped(tmp_path):
    """An arm of zero length has no defined angle; that frame is dropped."""
    collapsed = [[0, 0, 0], [0, 0, 0], [0, 1, 0]]   # i sits on the vertex
    traj = triplet_traj(tmp_path, [right_angle_frame(), collapsed, right_angle_frame()])
    a = make_analyser(tmp_path, traj, 3, 3)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    assert angle_values(a) == pytest.approx([90.0, 90.0])


# --------------------------------------------------------------------------- #
# Parser robustness (shared reader)                                             #
# --------------------------------------------------------------------------- #
def test_trailing_whitespace_lines_do_not_discard_the_run(tmp_path):
    traj = triplet_traj(tmp_path, [angle_frame(90), angle_frame(120)],
                        trailer="   \n\n")
    a = make_analyser(tmp_path, traj, 3, 2)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    assert angle_values(a) == pytest.approx([90.0, 120.0])


def test_truncated_final_frame_keeps_the_complete_frames(tmp_path):
    traj = triplet_traj(tmp_path, [angle_frame(90), angle_frame(120)])
    with open(traj, "a") as f:
        f.write("3\ncomment\nO 1.0 0.0 0.0\n")  # two atoms missing
    a = make_analyser(tmp_path, traj, 3, 3)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    assert angle_values(a) == pytest.approx([90.0, 120.0])
    assert "truncated" in a.main_window.messages().lower()


def test_frame_header_atom_count_mismatch_is_reported(tmp_path):
    traj = tmp_path / "traj.xyz"
    with open(traj, "w") as f:
        f.write("3\ncomment\nO 1.0 0.0 0.0\nC 0.0 0.0 0.0\nN 0.0 1.0 0.0\n")
        f.write("4\ncomment\nO 1.0 0.0 0.0\nC 0.0 0.0 0.0\nN 0.0 1.0 0.0\nH 2.0 0.0 0.0\n")
    a = make_analyser(tmp_path, traj, 3, 2)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    assert angle_values(a) == pytest.approx([90.0])
    assert "atom count" in a.main_window.messages().lower()


def test_unparseable_coordinates_on_unselected_atoms_are_skipped(tmp_path):
    traj = tmp_path / "traj.xyz"
    with open(traj, "w") as f:
        for _ in range(2):
            f.write("4\ncomment\n")
            f.write("O 1.0 0.0 0.0\nC 0.0 0.0 0.0\nN 0.0 1.0 0.0\n")
            f.write("X ????? ????? ?????\n")
    a = make_analyser(tmp_path, traj, 4, 2, atom_labels="1 2 3")
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    assert angle_values(a) == pytest.approx([90.0, 90.0])


def test_empty_trajectory_is_reported(tmp_path):
    traj = tmp_path / "empty.xyz"
    traj.write_text("")
    a = make_analyser(tmp_path, traj, 3, 0)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is False
    assert "no frames" in a.main_window.messages().lower()


# --------------------------------------------------------------------------- #
# Periodic boundaries                                                           #
# --------------------------------------------------------------------------- #
def test_blank_cell_field_leaves_vectors_unwrapped(tmp_path):
    """i at x=9, vertex at origin, k at (1,1,0): 45 degrees without PBC."""
    traj = triplet_traj(tmp_path, [[[9, 0, 0], [0, 0, 0], [1, 1, 0]]])
    a = make_analyser(tmp_path, traj, 3, 1, cell="")
    assert asyncio.run(a.read_params(None)) is True
    assert a.cell_lengths is None
    assert asyncio.run(a.bond_angle()) is True
    assert angle_values(a) == pytest.approx([45.0])


def test_minimum_image_changes_the_angle_across_a_boundary(tmp_path):
    """The same triplet in a 10 A box: the i arm wraps to -1, giving 135."""
    traj = triplet_traj(tmp_path, [[[9, 0, 0], [0, 0, 0], [1, 1, 0]]])
    a = make_analyser(tmp_path, traj, 3, 1, cell="10 10 10")
    assert asyncio.run(a.read_params(None)) is True
    assert a.cell_lengths == pytest.approx([10.0, 10.0, 10.0])
    assert asyncio.run(a.bond_angle()) is True
    assert angle_values(a) == pytest.approx([135.0])


def test_minimum_image_is_invariant_to_shifting_the_box(tmp_path):
    base = triplet_traj(tmp_path, [[[1, 1, 1], [0, 0, 0], [0, 2, 0]]], name="base.xyz")
    shifted = triplet_traj(tmp_path, [[[11, 1, -9], [0, 0, 0], [0, 2, 0]]],
                           name="shifted.xyz")
    a = make_analyser(tmp_path, base, 3, 1, cell="10 10 10")
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    b = make_analyser(tmp_path, shifted, 3, 1, cell="10 10 10")
    assert asyncio.run(b.read_params(None)) is True
    assert asyncio.run(b.bond_angle()) is True
    assert angle_values(a) == pytest.approx(angle_values(b))


def test_arm_longer_than_half_the_box_is_reported(tmp_path):
    """A diagonal arm can survive minimum image longer than min(a,b,c)/2."""
    traj = triplet_traj(tmp_path, [[[4, 4, 4], [0, 0, 0], [0, 2, 0]]])
    a = make_analyser(tmp_path, traj, 3, 1, cell="10 10 10")
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    assert "half" in a.main_window.messages().lower()


@pytest.mark.parametrize("cell", ["10 10", "10 10 10 10", "10 0 10", "10 -1 10", "a b c"])
def test_invalid_cell_lengths_are_rejected(tmp_path, cell):
    traj = triplet_traj(tmp_path, [right_angle_frame()])
    a = make_analyser(tmp_path, traj, 3, 1, cell=cell)
    assert asyncio.run(a.read_params(None)) is False


# --------------------------------------------------------------------------- #
# Histogram and statistics                                                      #
# --------------------------------------------------------------------------- #
def test_angles_beyond_max_angle_are_counted_not_clamped(tmp_path):
    """Out-of-range angles are excluded and reported, never piled into the end bin."""
    frames = [angle_frame(30), angle_frame(30), angle_frame(30), angle_frame(120)]
    traj = triplet_traj(tmp_path, frames)
    a = make_analyser(tmp_path, traj, 3, 4, max_angle="90", bin_width="1.0")
    a._reset_analysis_state()
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    assert asyncio.run(a.distribution_function()) is True
    assert a.out_of_range_count == 1
    assert a.histogram_sum == 3
    assert "out of range" in a.main_window.messages().lower()
    # The last bin must not have absorbed the 120 degree sample.
    assert a.histogram[-1] == pytest.approx(0.0)


def test_histogram_bins_match_floor_division(tmp_path):
    frames = [angle_frame(30.5), angle_frame(31.5), angle_frame(31.9), angle_frame(45.5)]
    traj = triplet_traj(tmp_path, frames)
    a = make_analyser(tmp_path, traj, 3, 4, max_angle="180", bin_width="1.0")
    assert run_analysis(a) is True
    counts = [round(p * a.histogram_sum / 100.0) for p in a.histogram]
    assert counts[30] == 1     # bin [30, 31)
    assert counts[31] == 2     # bin [31, 32)
    assert counts[45] == 1     # bin [45, 46)
    assert sum(counts) == 4


def test_statistics_match_sample_variance(tmp_path):
    wanted = [30.0, 45.0, 60.0, 75.0, 90.0]
    traj = triplet_traj(tmp_path, [angle_frame(w) for w in wanted])
    a = make_analyser(tmp_path, traj, 3, len(wanted))
    assert run_analysis(a) is True
    assert a.stats["variance"] == pytest.approx(variance(wanted))
    assert a.stats["std_dev"] == pytest.approx(math.sqrt(variance(wanted)))
    assert a.stats["average_angle"] == pytest.approx(sum(wanted) / len(wanted))
    assert a.stats["largest_angle"] == pytest.approx(max(wanted))
    assert a.stats["smallest_angle"] == pytest.approx(min(wanted))


# --------------------------------------------------------------------------- #
# Free energy                                                                   #
# --------------------------------------------------------------------------- #
def test_free_energy_without_jacobian_is_minus_rt_ln_p(tmp_path):
    frames = [angle_frame(30.5), angle_frame(30.5), angle_frame(60.5)]
    traj = triplet_traj(tmp_path, frames)
    a = make_analyser(tmp_path, traj, 3, 3, bin_width="1.0", temperature="300",
                      jacobian=False, smooth_fe=False)
    assert run_analysis(a) is True
    R, T = 0.001987204, 300.0
    expected = {30.5: -R * T * math.log(2 / 3), 60.5: -R * T * math.log(1 / 3)}
    got = {round(x, 4): fe for x, fe in a.free_energy_pairs}
    assert got == pytest.approx(expected)


def test_jacobian_switch_applies_the_sine_correction(tmp_path):
    """With the switch on, W = -RT ln[P(theta)/sin(theta)]."""
    frames = [angle_frame(30.5), angle_frame(30.5), angle_frame(60.5)]
    traj = triplet_traj(tmp_path, frames)
    a = make_analyser(tmp_path, traj, 3, 3, bin_width="1.0", temperature="300",
                      jacobian=True, smooth_fe=False)
    assert run_analysis(a) is True
    R, T = 0.001987204, 300.0
    expected = {
        30.5: -R * T * math.log((2 / 3) / math.sin(math.radians(30.5))),
        60.5: -R * T * math.log((1 / 3) / math.sin(math.radians(60.5))),
    }
    got = {round(x, 4): fe for x, fe in a.free_energy_pairs}
    assert got == pytest.approx(expected)


def test_smoothing_switch_changes_the_free_energy(tmp_path):
    wanted = [30.0 + i * 0.7 for i in range(40)]
    traj = triplet_traj(tmp_path, [angle_frame(w) for w in wanted])
    plain = make_analyser(tmp_path / "p", traj, 3, len(wanted), smooth_fe=False)
    smoothed = make_analyser(tmp_path / "s", traj, 3, len(wanted), smooth_fe=True)
    assert run_analysis(plain) is True
    assert run_analysis(smoothed) is True
    assert [fe for _, fe in plain.free_energy_pairs] != \
        pytest.approx([fe for _, fe in smoothed.free_energy_pairs])


# --------------------------------------------------------------------------- #
# Text-only reporting and outputs                                               #
# --------------------------------------------------------------------------- #
def test_no_png_files_are_written(tmp_path):
    traj = triplet_traj(tmp_path, [angle_frame(30), angle_frame(60), angle_frame(90)])
    a = make_analyser(tmp_path, traj, 3, 3, stub_plots=False)
    assert run_analysis(a) is True
    assert glob.glob(os.path.join(str(tmp_path), "*.png")) == []
    assert len(glob.glob(os.path.join(str(tmp_path), "*.dat"))) == 3
    assert len(a.saved_plot_data) == 3
    assert a.saved_plot_files == []


def test_frames_counter_reports_status_as_text(tmp_path):
    traj = triplet_traj(tmp_path, [right_angle_frame()] * 4)
    a = make_analyser(tmp_path, traj, 3, 0)

    async def fake_open(widget):
        a.trajec = str(traj)
        a.num_atoms = 3
        a.output_dir = str(tmp_path)

    a.open_file_dialog = fake_open
    asyncio.run(a.frames_counter(None))
    assert a.total_frame_number == 4
    assert "Trajectory loaded" in a.multi_line_text.value
    assert "frames: 4" in a.multi_line_text.value
    assert not hasattr(a, "progress_label")


def test_calculation_status_is_reported_as_text(tmp_path):
    traj = triplet_traj(tmp_path, [angle_frame(90), angle_frame(100)])
    a = make_analyser(tmp_path, traj, 3, 2)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is True
    assert "Calculating bond angles" in a.multi_line_text.value


def test_output_files_are_written_with_the_triplet_tag(tmp_path):
    traj = triplet_traj(tmp_path, [angle_frame(90), angle_frame(100)])
    a = make_analyser(tmp_path, traj, 3, 2)
    assert run_analysis(a) is True
    assert asyncio.run(a.save_summary()) is True
    for name in ("angles_O1_C2_N3.dat", "angles_distribution_O1_C2_N3.dat",
                 "angle_free_energy_O1_C2_N3.dat", "summary_O1_C2_N3.txt"):
        assert os.path.exists(tmp_path / name), name


def test_csv_export_writes_three_files(tmp_path):
    traj = triplet_traj(tmp_path, [angle_frame(90), angle_frame(100)])
    a = make_analyser(tmp_path, traj, 3, 2)
    assert run_analysis(a) is True
    assert asyncio.run(a.export_csv()) is True
    for name in ("angles_O1_C2_N3.csv", "angles_distribution_O1_C2_N3.csv",
                 "angle_free_energy_O1_C2_N3.csv"):
        assert os.path.exists(tmp_path / name), name


def test_final_summary_reports_parameters_and_results(tmp_path):
    traj = triplet_traj(tmp_path, [angle_frame(90), angle_frame(100), angle_frame(110)])
    a = make_analyser(tmp_path, traj, 3, 3, cell="10 10 10", bin_width="1.0",
                      temperature="300")
    assert run_analysis(a) is True
    assert asyncio.run(a.save_summary()) is True
    text = a.multi_line_text.value
    assert text.startswith("BOND ANGLE ANALYSIS - COMPLETED")
    for heading in ("PARAMETERS", "RESULTS"):
        assert heading in text, heading
    for token in ("Maximum angle", "Bin width", "Temperature", "300",
                  "O1-C2-N3", "10 10 10"):
        assert token in text, token
    for token in ("Largest angle", "Smallest angle", "Average angle",
                  "Std deviation", "Lowest free energy", "Frames processed"):
        assert token in text, token
    assert "Calculating bond angles" not in text
    with open(tmp_path / "summary_O1_C2_N3.txt") as f:
        assert f.read() == text


def test_unwritable_output_reports_instead_of_raising(tmp_path):
    traj = triplet_traj(tmp_path, [angle_frame(90), angle_frame(100)])
    os.makedirs(tmp_path / "angles_O1_C2_N3.dat", exist_ok=True)
    a = make_analyser(tmp_path, traj, 3, 2)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_angle()) is False
    assert a.main_window.dialogs != []


# --------------------------------------------------------------------------- #
# Parameter validation matrix                                                   #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("field,value", [
    ("max_angle", ""), ("max_angle", "0"), ("max_angle", "-1"), ("max_angle", "abc"),
    ("time_step", ""), ("time_step", "0"),
    ("sampling", ""), ("sampling", "0"), ("sampling", "1.5"),
    ("temperature", ""), ("temperature", "0"),
    ("bin_width", ""), ("bin_width", "0"), ("bin_width", "200"),
    ("atom_labels", ""), ("atom_labels", "1 2"), ("atom_labels", "1 2 3 4"),
    ("atom_labels", "0 1 2"), ("atom_labels", "1 2 99"), ("atom_labels", "x y z"),
    ("atom_labels", "1 2 2"),
])
def test_invalid_parameters_are_rejected(tmp_path, field, value):
    traj = triplet_traj(tmp_path, [right_angle_frame()])
    a = make_analyser(tmp_path, traj, 3, 1, **{field: value})
    assert asyncio.run(a.read_params(None)) is False
    assert a.main_window.dialogs != []
