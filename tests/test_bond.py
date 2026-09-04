"""Tests for the Bond Length Analysis tool (``bond.py``).

Run from ``venv/src/``::

    python -m pytest tests/ -q

The Toga GUI is never instantiated: dialogs, text inputs, switches and the
progress bar are routed through stubs, so everything runs headless. Following
the house pattern, these drive ``BondAnalyser`` (the logic class) directly
rather than ``BondUI``.

Physics is checked against distances that can be worked out by hand, so the
tests stay meaningful if the implementation is refactored again.
"""
import os
import glob
import math
import asyncio
from statistics import variance

import numpy as np
import pytest

import bond
from bond import BondAnalyser


# --------------------------------------------------------------------------- #
# Toga stubs                                                                    #
# --------------------------------------------------------------------------- #
class FakeDialog:
    def __init__(self, title, message):
        self.title, self.message = title, message


class FakeWindow:
    def __init__(self):
        self.dialogs = []          # list of (title, message)

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
    """Route every dialog through a capturable stub instead of the GUI backend."""
    monkeypatch.setattr(bond.toga, "InfoDialog", FakeDialog)


@pytest.fixture(autouse=True)
def _isolate_plot_state():
    """``DisplayPlots`` keeps its plot buffers as class attributes, so reset
    them per test to stop one test's curves leaking into the next."""
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
    """Write a standard multi-frame .xyz; ``trailer`` is appended verbatim."""
    with open(path, "w") as f:
        for coords in frames:
            f.write(f"{len(symbols)}\ncomment\n")
            for s, (x, y, z) in zip(symbols, coords):
                f.write(f"{s} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write(trailer)


def two_atom_traj(tmp_path, separations, symbols=("C", "O"), trailer=""):
    """Trajectory of two atoms placed ``d`` apart along x, one frame per d."""
    traj = tmp_path / "traj.xyz"
    frames = [np.array([[0.0, 0.0, 0.0], [d, 0.0, 0.0]]) for d in separations]
    write_xyz(traj, frames, list(symbols), trailer=trailer)
    return traj


def make_analyser(out_dir, traj, num_atoms, n_frames, *, max_r="5.0",
                  time_step="5.0", sampling="10", temperature="300",
                  atom_labels="1 2", bin_width="0.1", cell="",
                  jacobian=False, smooth_fe=False, stub_plots=True):
    # No progress bar and no progress label are attached: this tool reports
    # status as text in multi_line_text, so touching either would AttributeError.
    a = BondAnalyser()
    a.main_window = FakeWindow()
    a.multi_line_text = FakeInput()
    a.output_dir = str(out_dir)
    if stub_plots:
        a.save_plots = lambda *args, **kw: None
    a.display_plots = lambda *args, **kw: None
    a.textInput_max_r = FakeInput(max_r)
    a.textInput_time_step = FakeInput(time_step)
    a.textInput_sampling_interval = FakeInput(sampling)
    a.textInput_temperature = FakeInput(temperature)
    a.textInput_atom_labels = FakeInput(atom_labels)
    a.textInput_bin_width = FakeInput(bin_width)
    a.textInput_cell_lengths = FakeInput(cell)
    a.switch_show_plots = FakeSwitch(True)
    a.switch_save_csv = FakeSwitch(False)
    a.switch_jacobian = FakeSwitch(jacobian)
    a.switch_smooth_fe = FakeSwitch(smooth_fe)
    a.trajec = str(traj)
    a.num_atoms = num_atoms
    a.total_frame_number = n_frames
    return a


def run_analysis(a):
    """Drive the stage chain the way ``BondUI.workflow`` does."""
    a._reset_analysis_state()
    if not asyncio.run(a.read_params(None)):
        return False
    if not asyncio.run(a.bond_length(None)):
        return False
    if not asyncio.run(a.distribution_function()):
        return False
    return asyncio.run(a.free_energy())


def bond_values(a):
    return [r for _, r in a.bond_lengths]


# --------------------------------------------------------------------------- #
# Core geometry                                                                 #
# --------------------------------------------------------------------------- #
def test_bond_length_matches_hand_computed_distance(tmp_path):
    traj = two_atom_traj(tmp_path, [1.0, 1.5, 2.0])
    a = make_analyser(tmp_path, traj, 2, 3)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is True
    assert bond_values(a) == pytest.approx([1.0, 1.5, 2.0])


def test_time_axis_uses_timestep_times_sampling(tmp_path):
    traj = two_atom_traj(tmp_path, [1.0, 1.0])
    a = make_analyser(tmp_path, traj, 2, 2, time_step="5.0", sampling="10")
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is True
    step_ps = 5.0 * 10 * BondAnalyser.atufs / 1000.0
    assert [t for t, _ in a.bond_lengths] == pytest.approx([0.0, step_ps])


# --------------------------------------------------------------------------- #
# Parser robustness                                                             #
# --------------------------------------------------------------------------- #
def test_trailing_blank_line_does_not_discard_the_run(tmp_path):
    """A trailing newline is clean EOF, not a malformed frame.

    Regression: the old loop treated "\\n" as a frame header, raised
    "Invalid atom line format", and threw away every frame already parsed.
    """
    traj = two_atom_traj(tmp_path, [1.0, 1.2], trailer="\n")
    a = make_analyser(tmp_path, traj, 2, 2)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is True
    assert bond_values(a) == pytest.approx([1.0, 1.2])
    assert a.main_window.dialogs == []


def test_trailing_whitespace_lines_do_not_discard_the_run(tmp_path):
    traj = two_atom_traj(tmp_path, [1.0, 1.2], trailer="   \n\n")
    a = make_analyser(tmp_path, traj, 2, 2)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is True
    assert bond_values(a) == pytest.approx([1.0, 1.2])


def test_truncated_final_frame_keeps_the_complete_frames(tmp_path):
    """A half-written last frame warns but preserves the good data."""
    traj = two_atom_traj(tmp_path, [1.0, 1.2])
    with open(traj, "a") as f:
        f.write("2\ncomment\nC 0.000000 0.000000 0.000000\n")  # second atom missing
    a = make_analyser(tmp_path, traj, 2, 3)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is True
    assert bond_values(a) == pytest.approx([1.0, 1.2])
    assert "truncated" in a.main_window.messages().lower()


def test_frame_header_atom_count_mismatch_is_reported(tmp_path):
    """A frame declaring a different atom count would desync framing."""
    traj = tmp_path / "traj.xyz"
    with open(traj, "w") as f:
        f.write("2\ncomment\nC 0.0 0.0 0.0\nO 1.0 0.0 0.0\n")
        f.write("3\ncomment\nC 0.0 0.0 0.0\nO 1.5 0.0 0.0\nH 2.0 0.0 0.0\n")
    a = make_analyser(tmp_path, traj, 2, 2)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is True
    assert bond_values(a) == pytest.approx([1.0])
    assert "atom count" in a.main_window.messages().lower()


def test_malformed_selected_atom_line_is_reported(tmp_path):
    traj = tmp_path / "traj.xyz"
    with open(traj, "w") as f:
        f.write("2\ncomment\nC nope nope nope\nO 1.0 0.0 0.0\n")
    a = make_analyser(tmp_path, traj, 2, 1)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is False
    assert a.main_window.dialogs != []


def test_unparseable_coordinates_on_unselected_atoms_are_skipped(tmp_path):
    """Only the two selected atoms are parsed, so junk elsewhere is inert.

    This is the observable half of the performance fix: the reader no longer
    splits and float-converts every atom line in every frame.
    """
    traj = tmp_path / "traj.xyz"
    with open(traj, "w") as f:
        for d in (1.0, 1.4):
            f.write("3\ncomment\n")
            f.write("C 0.0 0.0 0.0\n")
            f.write(f"O {d} 0.0 0.0\n")
            f.write("X ????? ????? ?????\n")
    a = make_analyser(tmp_path, traj, 3, 2, atom_labels="1 2")
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is True
    assert bond_values(a) == pytest.approx([1.0, 1.4])


# --------------------------------------------------------------------------- #
# Periodic boundaries (optional, orthorhombic minimum image)                     #
# --------------------------------------------------------------------------- #
def test_blank_cell_field_leaves_distances_unwrapped(tmp_path):
    traj = two_atom_traj(tmp_path, [9.0])
    a = make_analyser(tmp_path, traj, 2, 1, cell="", max_r="12")
    assert asyncio.run(a.read_params(None)) is True
    assert a.cell_lengths is None
    assert asyncio.run(a.bond_length(None)) is True
    assert bond_values(a) == pytest.approx([9.0])


def test_minimum_image_shortens_a_bond_across_the_boundary(tmp_path):
    """Atoms at x=0 and x=9 in a 10 A box are 1 A apart, not 9 A."""
    traj = two_atom_traj(tmp_path, [9.0])
    a = make_analyser(tmp_path, traj, 2, 1, cell="10 10 10", max_r="5")
    assert asyncio.run(a.read_params(None)) is True
    assert a.cell_lengths == pytest.approx([10.0, 10.0, 10.0])
    assert asyncio.run(a.bond_length(None)) is True
    assert bond_values(a) == pytest.approx([1.0])


def test_minimum_image_is_invariant_to_shifting_the_box(tmp_path):
    """Translating every atom by a whole box vector must not change the bond."""
    base = tmp_path / "base.xyz"
    write_xyz(base, [np.array([[1.0, 1.0, 1.0], [3.5, 1.0, 1.0]])], ["C", "O"])
    shifted = tmp_path / "shifted.xyz"
    write_xyz(shifted, [np.array([[11.0, 1.0, -9.0], [3.5, 1.0, 1.0]])], ["C", "O"])

    a = make_analyser(tmp_path, base, 2, 1, cell="10 10 10")
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is True

    b = make_analyser(tmp_path, shifted, 2, 1, cell="10 10 10")
    assert asyncio.run(b.read_params(None)) is True
    assert asyncio.run(b.bond_length(None)) is True

    assert bond_values(a) == pytest.approx(bond_values(b))


def test_max_r_beyond_half_the_box_is_rejected(tmp_path):
    """Minimum image cannot resolve distances past min(a,b,c)/2."""
    traj = two_atom_traj(tmp_path, [1.0])
    a = make_analyser(tmp_path, traj, 2, 1, cell="10 10 10", max_r="6")
    assert asyncio.run(a.read_params(None)) is False
    assert "half" in a.main_window.messages().lower()


@pytest.mark.parametrize("cell", ["10 10", "10 10 10 10", "10 0 10", "10 -1 10", "a b c"])
def test_invalid_cell_lengths_are_rejected(tmp_path, cell):
    traj = two_atom_traj(tmp_path, [1.0])
    a = make_analyser(tmp_path, traj, 2, 1, cell=cell)
    assert asyncio.run(a.read_params(None)) is False


# --------------------------------------------------------------------------- #
# Histogram and statistics                                                      #
# --------------------------------------------------------------------------- #
def test_bonds_beyond_max_r_are_counted_and_reported(tmp_path):
    """Out-of-range samples are dropped from the histogram, so say so."""
    traj = two_atom_traj(tmp_path, [1.0, 1.0, 1.0, 4.0])
    a = make_analyser(tmp_path, traj, 2, 4, max_r="2.0", bin_width="0.1")
    a._reset_analysis_state()
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is True
    assert asyncio.run(a.distribution_function()) is True
    assert a.out_of_range_count == 1
    assert a.histogram_sum == 3
    assert "out of range" in a.main_window.messages().lower()


def test_histogram_bins_match_floor_division(tmp_path):
    """Binning semantics are unchanged: bin index is floor(r / bin_width)."""
    traj = two_atom_traj(tmp_path, [0.05, 0.15, 0.15, 0.25])
    a = make_analyser(tmp_path, traj, 2, 4, max_r="1.0", bin_width="0.1")
    assert run_analysis(a) is True
    counts = [round(p * a.histogram_sum / 100.0) for p in a.histogram]
    assert counts[:3] == [1, 2, 1]
    assert sum(counts) == 4


def test_statistics_match_sample_variance(tmp_path):
    """Characterisation: moving to numpy must keep the ddof=1 convention."""
    seps = [1.0, 1.1, 1.25, 1.4, 1.05]
    traj = two_atom_traj(tmp_path, seps)
    a = make_analyser(tmp_path, traj, 2, len(seps), max_r="3.0")
    assert run_analysis(a) is True
    assert a.stats["variance"] == pytest.approx(variance(seps))
    assert a.stats["std_dev"] == pytest.approx(math.sqrt(variance(seps)))
    assert a.stats["average_bond"] == pytest.approx(sum(seps) / len(seps))
    assert a.stats["largest_bond"] == pytest.approx(max(seps))
    assert a.stats["smallest_bond"] == pytest.approx(min(seps))


# --------------------------------------------------------------------------- #
# Free energy                                                                   #
# --------------------------------------------------------------------------- #
def test_free_energy_without_jacobian_is_minus_rt_ln_p(tmp_path):
    traj = two_atom_traj(tmp_path, [1.05, 1.05, 1.15])
    a = make_analyser(tmp_path, traj, 2, 3, max_r="2.0", bin_width="0.1",
                      temperature="300", jacobian=False)
    assert run_analysis(a) is True
    R, T = 0.001987204, 300.0
    expected = {round(r, 4): -R * T * math.log(p)
                for r, p in [(1.05, 2 / 3), (1.15, 1 / 3)]}
    got = {round(r, 4): fe for r, fe in a.free_energy_pairs}
    assert got == pytest.approx(expected)


def test_jacobian_switch_applies_the_r_squared_correction(tmp_path):
    """With the switch on, W(r) = -RT ln[P(r)/r^2]."""
    traj = two_atom_traj(tmp_path, [1.05, 1.05, 1.15])
    a = make_analyser(tmp_path, traj, 2, 3, max_r="2.0", bin_width="0.1",
                      temperature="300", jacobian=True)
    assert run_analysis(a) is True
    R, T = 0.001987204, 300.0
    expected = {round(r, 4): -R * T * math.log(p / r ** 2)
                for r, p in [(1.05, 2 / 3), (1.15, 1 / 3)]}
    got = {round(r, 4): fe for r, fe in a.free_energy_pairs}
    assert got == pytest.approx(expected)


def test_jacobian_switch_changes_the_free_energy(tmp_path):
    """Guard against the switch being read but ignored."""
    seps = [1.05, 1.05, 1.15]
    traj = two_atom_traj(tmp_path, seps)
    plain = make_analyser(tmp_path / "p", traj, 2, 3, max_r="2.0", jacobian=False)
    corrected = make_analyser(tmp_path / "c", traj, 2, 3, max_r="2.0", jacobian=True)
    assert run_analysis(plain) is True
    assert run_analysis(corrected) is True
    plain_fe = [fe for _, fe in plain.free_energy_pairs]
    corrected_fe = [fe for _, fe in corrected.free_energy_pairs]
    assert plain_fe != pytest.approx(corrected_fe)


# --------------------------------------------------------------------------- #
# Output and error handling                                                     #
# --------------------------------------------------------------------------- #
def test_output_files_are_written_with_the_pair_tag(tmp_path):
    traj = two_atom_traj(tmp_path, [1.0, 1.1], symbols=("C", "O"))
    a = make_analyser(tmp_path, traj, 2, 2, max_r="3.0")
    assert run_analysis(a) is True
    assert asyncio.run(a.save_summary()) is True
    for name in ("bond_C1_O2.dat", "bond_distribution_C1_O2.dat",
                 "free_energy_C1_O2.dat", "summary_C1_O2.txt"):
        assert os.path.exists(tmp_path / name), name


def test_unwritable_output_reports_instead_of_raising(tmp_path):
    """A blocked output path must surface a dialog, not an unhandled exception."""
    traj = two_atom_traj(tmp_path, [1.0, 1.1], symbols=("C", "O"))
    # Occupy the .dat path with a directory so opening it for writing fails.
    os.makedirs(tmp_path / "bond_C1_O2.dat", exist_ok=True)
    a = make_analyser(tmp_path, traj, 2, 2, max_r="3.0")
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is False
    assert a.main_window.dialogs != []


def test_empty_trajectory_is_reported(tmp_path):
    traj = tmp_path / "empty.xyz"
    traj.write_text("")
    a = make_analyser(tmp_path, traj, 2, 0)
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is False
    assert "no frames" in a.main_window.messages().lower()


def test_csv_export_writes_three_files(tmp_path):
    traj = two_atom_traj(tmp_path, [1.0, 1.1], symbols=("C", "O"))
    a = make_analyser(tmp_path, traj, 2, 2, max_r="3.0")
    assert run_analysis(a) is True
    assert asyncio.run(a.export_csv()) is True
    for name in ("bond_C1_O2.csv", "bond_distribution_C1_O2.csv",
                 "free_energy_C1_O2.csv"):
        assert os.path.exists(tmp_path / name), name


# --------------------------------------------------------------------------- #
# Text-only reporting: no PNG clutter, no progress bar, no frame label           #
# --------------------------------------------------------------------------- #
def test_no_png_files_are_written(tmp_path):
    """Figures go to the interactive viewer only; data files are still written."""
    traj = two_atom_traj(tmp_path, [1.0, 1.1, 1.2], symbols=("C", "O"))
    a = make_analyser(tmp_path, traj, 2, 3, max_r="3.0", stub_plots=False)
    assert run_analysis(a) is True

    assert glob.glob(os.path.join(str(tmp_path), "*.png")) == []
    assert len(glob.glob(os.path.join(str(tmp_path), "*.dat"))) == 3
    # The curves are still recorded for the interactive viewer.
    assert len(a.saved_plot_data) == 3
    assert a.saved_plot_files == []


def test_frames_counter_reports_status_as_text(tmp_path):
    """Loading a trajectory reports into the text box, not a frame-count label."""
    traj = two_atom_traj(tmp_path, [1.0, 1.1, 1.2, 1.3])
    a = make_analyser(tmp_path, traj, 2, 0)  # frame count not known yet

    async def fake_open(widget):  # stand in for the GUI file-open dialog
        a.trajec = str(traj)
        a.num_atoms = 2
        a.output_dir = str(tmp_path)

    a.open_file_dialog = fake_open
    asyncio.run(a.frames_counter(None))

    assert a.total_frame_number == 4
    assert "Trajectory loaded" in a.multi_line_text.value
    assert "frames: 4" in a.multi_line_text.value
    assert not hasattr(a, "progress_label")  # no frame label is used


def test_calculation_status_is_reported_as_text(tmp_path):
    traj = two_atom_traj(tmp_path, [1.0, 1.1])
    a = make_analyser(tmp_path, traj, 2, 2, max_r="3.0")
    assert asyncio.run(a.read_params(None)) is True
    assert asyncio.run(a.bond_length(None)) is True
    assert "Calculating bond lengths" in a.multi_line_text.value


def test_final_summary_reports_parameters_and_results(tmp_path):
    """The run ends with one summary of the inputs used and what came out."""
    traj = two_atom_traj(tmp_path, [1.0, 1.1, 1.2], symbols=("C", "O"))
    a = make_analyser(tmp_path, traj, 2, 3, max_r="3.0", bin_width="0.1",
                      temperature="300", cell="10 10 10")
    assert run_analysis(a) is True
    assert asyncio.run(a.save_summary()) is True

    text = a.multi_line_text.value
    for heading in ("PARAMETERS", "RESULTS"):
        assert heading in text, heading
    # Parameters actually used
    for token in ("Maximum r", "3.0", "Bin width", "0.1", "Temperature", "300",
                  "C1-O2", "10 10 10"):
        assert token in text, token
    # Results obtained
    for token in ("Largest bond", "Smallest bond", "Average bond",
                  "Std deviation", "Lowest free energy", "Frames processed"):
        assert token in text, token
    # No leftover progress chatter from the calculation stages.
    assert "Calculating bond lengths" not in text
    # The same summary is what lands in the .txt file.
    with open(tmp_path / "summary_C1_O2.txt") as f:
        assert f.read() == text


def test_summary_reports_the_jacobian_and_pbc_settings(tmp_path):
    traj = two_atom_traj(tmp_path, [1.0, 1.1, 1.2], symbols=("C", "O"))
    a = make_analyser(tmp_path, traj, 2, 3, max_r="3.0", cell="", jacobian=True)
    assert run_analysis(a) is True
    assert asyncio.run(a.save_summary()) is True
    text = a.multi_line_text.value
    assert "none" in text                      # no periodic boundaries
    assert "Lowest PMF" in text                # Jacobian form was used
    assert "yes" in text


# --------------------------------------------------------------------------- #
# Parameter validation matrix                                                   #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("field,value", [
    ("max_r", ""), ("max_r", "0"), ("max_r", "-1"), ("max_r", "abc"),
    ("time_step", ""), ("time_step", "0"),
    ("sampling", ""), ("sampling", "0"), ("sampling", "1.5"),
    ("temperature", ""), ("temperature", "0"),
    ("bin_width", ""), ("bin_width", "0"), ("bin_width", "10"),
    ("atom_labels", ""), ("atom_labels", "1"), ("atom_labels", "1 2 3"),
    ("atom_labels", "0 1"), ("atom_labels", "1 99"), ("atom_labels", "x y"),
])
def test_invalid_parameters_are_rejected(tmp_path, field, value):
    traj = two_atom_traj(tmp_path, [1.0])
    a = make_analyser(tmp_path, traj, 2, 1, **{field: value})
    assert asyncio.run(a.read_params(None)) is False
    assert a.main_window.dialogs != []
