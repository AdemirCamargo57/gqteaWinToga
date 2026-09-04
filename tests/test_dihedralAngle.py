"""Tests for the Dihedral Angle Analysis tool (``dihedralAngle.py``).

Run from ``venv/src/``::

    python -m pytest tests/ -q

Mirrors tests/test_bond.py and tests/test_bondAngle.py: the Toga GUI is never
instantiated, and these drive ``DihedralAngleAnalyser`` rather than the UI class.

Test geometry: with i=(0,1,0), j=(0,0,0), k=(1,0,0) and
l = k + (0, cos phi, -sin phi), the i-j-k-l dihedral is exactly +phi, which
makes every expected value here hand-checkable.
"""
import os
import glob
import math
import asyncio
from statistics import variance

import numpy as np
import pytest

import dihedralAngle
from dihedralAngle import DihedralAngleAnalyser


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
    monkeypatch.setattr(dihedralAngle.toga, "InfoDialog", FakeDialog)


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
SYMBOLS = ("C", "N", "O", "S")


def write_xyz(path, frames, symbols, trailer=""):
    with open(path, "w") as f:
        for coords in frames:
            f.write(f"{len(symbols)}\ncomment\n")
            for s, (x, y, z) in zip(symbols, coords):
                f.write(f"{s} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write(trailer)


def dihedral_frame(phi_deg):
    """Four atoms whose i-j-k-l dihedral is exactly `phi_deg`."""
    r = math.radians(phi_deg)
    return [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, math.cos(r), -math.sin(r)]]


def quad_traj(tmp_path, frames, symbols=SYMBOLS, trailer="", name="traj.xyz"):
    traj = tmp_path / name
    write_xyz(traj, [np.asarray(fr, float) for fr in frames], list(symbols),
              trailer=trailer)
    return traj


def make_analyser(out_dir, traj, num_atoms, n_frames, *, max_angle="360",
                  time_step="5.0", sampling="10", temperature="300",
                  atom_labels="1 2 3 4", bin_width="1.0", cell="",
                  wrap_180=False, smooth_fe=False, stub_plots=True):
    # No progress bar and no progress label: this tool reports status as text.
    a = DihedralAngleAnalyser()
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
    a.switch_wrap_180 = FakeSwitch(wrap_180)
    a.switch_smooth_fe = FakeSwitch(smooth_fe)
    a.trajec = str(traj)
    a.num_atoms = num_atoms
    a.total_frame_number = n_frames
    return a


def run_analysis(a):
    a._reset_analysis_state()
    if not asyncio.run(a.read_params(None)):
        return False
    if not asyncio.run(a.dihedral_angle(None)):
        return False
    if not asyncio.run(a.distribution_function()):
        return False
    return asyncio.run(a.dihedral_free_energy())


def dihedral_values(a):
    return [v for _, v in a.dihedral]


# --------------------------------------------------------------------------- #
# Core geometry                                                                 #
# --------------------------------------------------------------------------- #
# Coordinates are written to the .xyz with six decimals, so a recovered angle
# carries about 1e-5 degrees of quantisation error. Assert to 1e-4.
ANGLE_TOL = 1e-4


@pytest.mark.parametrize("phi", [0.0, 30.0, 90.0, 150.0, 179.0])
def test_positive_dihedral_matches_hand_computed_value(tmp_path, phi):
    traj = quad_traj(tmp_path, [dihedral_frame(phi)])
    a = make_analyser(tmp_path, traj, 4, 1, wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    assert dihedral_values(a) == pytest.approx([phi], abs=ANGLE_TOL)


@pytest.mark.parametrize("phi", [-30.0, -90.0, -150.0])
def test_negative_dihedral_is_kept_when_wrapping(tmp_path, phi):
    traj = quad_traj(tmp_path, [dihedral_frame(phi)])
    a = make_analyser(tmp_path, traj, 4, 1, wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    assert dihedral_values(a) == pytest.approx([phi], abs=ANGLE_TOL)


@pytest.mark.parametrize("phi,expected", [(-30.0, 330.0), (-90.0, 270.0), (60.0, 60.0)])
def test_without_wrapping_angles_map_to_0_360(tmp_path, phi, expected):
    traj = quad_traj(tmp_path, [dihedral_frame(phi)])
    a = make_analyser(tmp_path, traj, 4, 1, wrap_180=False)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    assert dihedral_values(a) == pytest.approx([expected], abs=ANGLE_TOL)


def test_time_axis_uses_timestep_times_sampling(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60), dihedral_frame(60)])
    a = make_analyser(tmp_path, traj, 4, 2, time_step="5.0", sampling="10")
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    step_ps = 5.0 * 10 * DihedralAngleAnalyser.atufs / 1000.0
    assert [t for t, _ in a.dihedral] == pytest.approx([0.0, step_ps])


def test_collinear_geometry_is_skipped(tmp_path):
    """Collinear atoms give no dihedral plane; that frame is dropped."""
    collinear = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]
    traj = quad_traj(tmp_path, [dihedral_frame(60), collinear, dihedral_frame(60)])
    a = make_analyser(tmp_path, traj, 4, 3, wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    assert dihedral_values(a) == pytest.approx([60.0, 60.0])
    assert a.skipped_degenerate == 1


# --------------------------------------------------------------------------- #
# The wrap-mode histogram regression                                            #
# --------------------------------------------------------------------------- #
def test_wrapped_negative_angles_are_histogrammed(tmp_path):
    """Regression: negative angles used to get a negative bin index and vanish.

    With wrap on, the histogram must span [-180, 180], not [0, max_angle].
    """
    traj = quad_traj(tmp_path, [dihedral_frame(-90), dihedral_frame(-90),
                                dihedral_frame(90)])
    a = make_analyser(tmp_path, traj, 4, 3, wrap_180=True, bin_width="1.0")
    assert run_analysis(a) is True
    assert a.histogram_sum == 3          # all three counted, not just the positive one
    assert a.out_of_range_count == 0
    assert min(a.bin_centers) < 0        # the range covers negative angles
    assert a.bin_centers[0] == pytest.approx(-179.5)
    assert a.bin_centers[-1] == pytest.approx(179.5)


def test_wrapped_free_energy_covers_both_signs(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(-90), dihedral_frame(-90),
                                dihedral_frame(90)])
    a = make_analyser(tmp_path, traj, 4, 3, wrap_180=True, bin_width="1.0",
                      smooth_fe=False)
    assert run_analysis(a) is True
    populated = {round(x, 1) for x, _ in a.free_energy_pairs}
    assert populated == {-89.5, 90.5}


# --------------------------------------------------------------------------- #
# Parser robustness (shared reader)                                             #
# --------------------------------------------------------------------------- #
def test_trailing_whitespace_lines_do_not_discard_the_run(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60), dihedral_frame(90)],
                     trailer="   \n\n")
    a = make_analyser(tmp_path, traj, 4, 2, wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    assert dihedral_values(a) == pytest.approx([60.0, 90.0])


def test_truncated_final_frame_keeps_the_complete_frames(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60), dihedral_frame(90)])
    with open(traj, "a") as f:
        f.write("4\ncomment\nC 0.0 1.0 0.0\n")
    a = make_analyser(tmp_path, traj, 4, 3, wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    assert dihedral_values(a) == pytest.approx([60.0, 90.0])
    assert "truncated" in a.main_window.messages().lower()


def test_frame_header_atom_count_mismatch_is_reported(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60)])
    with open(traj, "a") as f:
        f.write("5\ncomment\n")
        for row in dihedral_frame(90) + [[9, 9, 9]]:
            f.write(f"C {row[0]:.6f} {row[1]:.6f} {row[2]:.6f}\n")
    a = make_analyser(tmp_path, traj, 4, 2, wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    assert dihedral_values(a) == pytest.approx([60.0])
    assert "atom count" in a.main_window.messages().lower()


def test_unparseable_coordinates_on_unselected_atoms_are_skipped(tmp_path):
    traj = tmp_path / "traj.xyz"
    with open(traj, "w") as f:
        for _ in range(2):
            f.write("5\ncomment\n")
            for s, row in zip(SYMBOLS, dihedral_frame(60)):
                f.write(f"{s} {row[0]:.6f} {row[1]:.6f} {row[2]:.6f}\n")
            f.write("X ????? ????? ?????\n")
    a = make_analyser(tmp_path, traj, 5, 2, wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    assert dihedral_values(a) == pytest.approx([60.0, 60.0])


def test_empty_trajectory_is_reported(tmp_path):
    traj = tmp_path / "empty.xyz"
    traj.write_text("")
    a = make_analyser(tmp_path, traj, 4, 0)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is False
    assert "no frames" in a.main_window.messages().lower()


# --------------------------------------------------------------------------- #
# Periodic boundaries                                                           #
# --------------------------------------------------------------------------- #
def test_minimum_image_changes_the_dihedral_across_a_boundary(tmp_path):
    """Shifting l by one box vector must not change the wrapped result."""
    base = quad_traj(tmp_path, [dihedral_frame(60)], name="base.xyz")
    shifted_frame = [list(p) for p in dihedral_frame(60)]
    shifted_frame[3] = [shifted_frame[3][0] + 10.0, shifted_frame[3][1],
                        shifted_frame[3][2] - 10.0]
    shifted = quad_traj(tmp_path, [shifted_frame], name="shifted.xyz")

    a = make_analyser(tmp_path, base, 4, 1, cell="10 10 10", wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True

    b = make_analyser(tmp_path, shifted, 4, 1, cell="10 10 10", wrap_180=True)
    assert asyncio.run(b.read_params(None)) is True
    assert asyncio.run(b.dihedral_angle(None)) is True

    assert dihedral_values(a) == pytest.approx(dihedral_values(b))


def test_blank_cell_field_disables_periodic_boundaries(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60)])
    a = make_analyser(tmp_path, traj, 4, 1, cell="", wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert a.cell_lengths is None


@pytest.mark.parametrize("cell", ["10 10", "10 10 10 10", "10 0 10", "a b c"])
def test_invalid_cell_lengths_are_rejected(tmp_path, cell):
    traj = quad_traj(tmp_path, [dihedral_frame(60)])
    a = make_analyser(tmp_path, traj, 4, 1, cell=cell)
    assert asyncio.run(a.read_params(None)) is False


# --------------------------------------------------------------------------- #
# Histogram and statistics                                                      #
# --------------------------------------------------------------------------- #
def test_angles_beyond_max_angle_are_counted_not_clamped(tmp_path):
    frames = [dihedral_frame(30), dihedral_frame(30), dihedral_frame(30),
              dihedral_frame(150)]
    traj = quad_traj(tmp_path, frames)
    a = make_analyser(tmp_path, traj, 4, 4, max_angle="90", bin_width="1.0",
                      wrap_180=False)
    a._reset_analysis_state()
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    assert asyncio.run(a.distribution_function()) is True
    assert a.out_of_range_count == 1
    assert a.histogram_sum == 3
    assert "out of range" in a.main_window.messages().lower()
    assert a.histogram[-1] == pytest.approx(0.0)


def test_statistics_match_sample_variance(tmp_path):
    wanted = [30.0, 45.0, 60.0, 75.0, 90.0]
    traj = quad_traj(tmp_path, [dihedral_frame(w) for w in wanted])
    a = make_analyser(tmp_path, traj, 4, len(wanted), wrap_180=True)
    assert run_analysis(a) is True
    assert a.stats["variance"] == pytest.approx(variance(wanted))
    assert a.stats["std_dev"] == pytest.approx(math.sqrt(variance(wanted)))
    assert a.stats["average_dihedral"] == pytest.approx(sum(wanted) / len(wanted))
    assert a.stats["largest_dihedral"] == pytest.approx(max(wanted))
    assert a.stats["smallest_dihedral"] == pytest.approx(min(wanted))


# --------------------------------------------------------------------------- #
# Free energy (uniform measure: no Jacobian for a dihedral)                      #
# --------------------------------------------------------------------------- #
def test_free_energy_is_minus_rt_ln_p(tmp_path):
    frames = [dihedral_frame(30.5), dihedral_frame(30.5), dihedral_frame(60.5)]
    traj = quad_traj(tmp_path, frames)
    a = make_analyser(tmp_path, traj, 4, 3, bin_width="1.0", temperature="300",
                      wrap_180=True, smooth_fe=False)
    assert run_analysis(a) is True
    R, T = 0.001987204, 300.0
    expected = {30.5: -R * T * math.log(2 / 3), 60.5: -R * T * math.log(1 / 3)}
    got = {round(x, 4): fe for x, fe in a.free_energy_pairs}
    assert got == pytest.approx(expected)


def test_no_jacobian_switch_is_exposed(tmp_path):
    """A dihedral has a uniform measure, so there is nothing to divide out."""
    traj = quad_traj(tmp_path, [dihedral_frame(60)])
    a = make_analyser(tmp_path, traj, 4, 1)
    assert not hasattr(a, "switch_use_jacobian")
    assert not hasattr(a, "switch_jacobian")
    assert asyncio.run(a.read_params(None)) is True
    assert getattr(a, "apply_jacobian", False) is False


def test_smoothing_switch_changes_the_free_energy(tmp_path):
    wanted = [30.0 + i * 0.7 for i in range(40)]
    traj = quad_traj(tmp_path, [dihedral_frame(w) for w in wanted])
    plain = make_analyser(tmp_path / "p", traj, 4, len(wanted), wrap_180=True,
                          smooth_fe=False)
    smoothed = make_analyser(tmp_path / "s", traj, 4, len(wanted), wrap_180=True,
                             smooth_fe=True)
    assert run_analysis(plain) is True
    assert run_analysis(smoothed) is True
    assert [fe for _, fe in plain.free_energy_pairs] != \
        pytest.approx([fe for _, fe in smoothed.free_energy_pairs])


# --------------------------------------------------------------------------- #
# Text-only reporting and outputs                                               #
# --------------------------------------------------------------------------- #
def test_no_png_files_are_written(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(30), dihedral_frame(60),
                                dihedral_frame(90)])
    a = make_analyser(tmp_path, traj, 4, 3, wrap_180=True, stub_plots=False)
    assert run_analysis(a) is True
    assert glob.glob(os.path.join(str(tmp_path), "*.png")) == []
    assert len(glob.glob(os.path.join(str(tmp_path), "*.dat"))) == 3
    assert len(a.saved_plot_data) == 3
    assert a.saved_plot_files == []


def test_frames_counter_reports_status_as_text(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60)] * 4)
    a = make_analyser(tmp_path, traj, 4, 0)

    async def fake_open(widget):
        a.trajec = str(traj)
        a.num_atoms = 4
        a.output_dir = str(tmp_path)

    a.open_file_dialog = fake_open
    asyncio.run(a.frames_counter(None))
    assert a.total_frame_number == 4
    assert "Trajectory loaded" in a.multi_line_text.value
    assert "frames: 4" in a.multi_line_text.value
    assert not hasattr(a, "progress_label")


def test_calculation_status_is_reported_as_text(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60), dihedral_frame(90)])
    a = make_analyser(tmp_path, traj, 4, 2, wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is True
    assert "Calculating dihedral angles" in a.multi_line_text.value


def test_output_files_are_written_with_the_quadruplet_tag(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60), dihedral_frame(90)])
    a = make_analyser(tmp_path, traj, 4, 2, wrap_180=True)
    assert run_analysis(a) is True
    assert asyncio.run(a.save_summary()) is True
    for name in ("dihedral_C1_N2_O3_S4.dat", "dihedral_distribution_C1_N2_O3_S4.dat",
                 "dihedral_free_energy_C1_N2_O3_S4.dat", "summary_C1_N2_O3_S4.txt"):
        assert os.path.exists(tmp_path / name), name


def test_csv_export_writes_three_files(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60), dihedral_frame(90)])
    a = make_analyser(tmp_path, traj, 4, 2, wrap_180=True)
    assert run_analysis(a) is True
    assert asyncio.run(a.export_csv()) is True
    for name in ("dihedral_C1_N2_O3_S4.csv", "dihedral_distribution_C1_N2_O3_S4.csv",
                 "dihedral_free_energy_C1_N2_O3_S4.csv"):
        assert os.path.exists(tmp_path / name), name


def test_final_summary_reports_parameters_and_results(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60), dihedral_frame(90),
                                dihedral_frame(120)])
    a = make_analyser(tmp_path, traj, 4, 3, cell="10 10 10", bin_width="1.0",
                      temperature="300", wrap_180=True)
    assert run_analysis(a) is True
    assert asyncio.run(a.save_summary()) is True
    text = a.multi_line_text.value
    assert text.startswith("DIHEDRAL ANGLE ANALYSIS - COMPLETED")
    for heading in ("PARAMETERS", "RESULTS"):
        assert heading in text, heading
    for token in ("Bin width", "Temperature", "300", "C1-N2-O3-S4", "10 10 10",
                  "Wrap to [-180,180]"):
        assert token in text, token
    for token in ("Largest dihedral", "Smallest dihedral", "Average dihedral",
                  "Std deviation", "Lowest free energy", "Frames processed"):
        assert token in text, token
    assert "Calculating dihedral angles" not in text
    with open(tmp_path / "summary_C1_N2_O3_S4.txt") as f:
        assert f.read() == text


def test_unwritable_output_reports_instead_of_raising(tmp_path):
    traj = quad_traj(tmp_path, [dihedral_frame(60), dihedral_frame(90)])
    os.makedirs(tmp_path / "dihedral_C1_N2_O3_S4.dat", exist_ok=True)
    a = make_analyser(tmp_path, traj, 4, 2, wrap_180=True)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.dihedral_angle(None)) is False
    assert a.main_window.dialogs != []


# --------------------------------------------------------------------------- #
# Parameter validation matrix                                                   #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("field,value", [
    ("max_angle", ""), ("max_angle", "0"), ("max_angle", "-1"), ("max_angle", "abc"),
    ("max_angle", "400"),
    ("time_step", ""), ("time_step", "0"),
    ("sampling", ""), ("sampling", "0"), ("sampling", "1.5"),
    ("temperature", ""), ("temperature", "0"),
    ("bin_width", ""), ("bin_width", "0"), ("bin_width", "400"),
    ("atom_labels", ""), ("atom_labels", "1 2 3"), ("atom_labels", "1 2 3 4 5"),
    ("atom_labels", "0 1 2 3"), ("atom_labels", "1 2 3 99"),
    ("atom_labels", "w x y z"), ("atom_labels", "1 2 3 3"),
])
def test_invalid_parameters_are_rejected(tmp_path, field, value):
    traj = quad_traj(tmp_path, [dihedral_frame(60)])
    a = make_analyser(tmp_path, traj, 4, 1, **{field: value})
    assert asyncio.run(a.read_params(None)) is False
    assert a.main_window.dialogs != []
