"""Classical MRT tool: tolerance in frames or fs, and exit counting.

The tolerance is the largest number of consecutive sampled frames an atom may
spend outside the shell without the excursion counting as an exit. These
tests pin the fs -> frames conversion (rounded down, exact for multiples of
the frame interval), that equivalent frames/fs tolerances give identical
output files, and that a zero tolerance counts every exit. The UI half is
driven headless through ``__new__`` and stub widgets; no window is created.
"""
import asyncio
import importlib.util
import math
import threading
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import meanResidenceTime_classical as mc
from meanResidenceTime_classical import ClassicalMeanResidenceTimeUI, MRTAnalyzer

# Occupancy of the moving O atom, frame by frame (1 = inside the 2.5 A shell):
# gaps of 1, 2 and 3 frames between in-shell frames.
PATTERN = [1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1]


# ------------------------------------------------------------------
# Conversion
# ------------------------------------------------------------------
class TestConversion:
    def test_frame_interval_from_au_and_fs(self):
        assert mc.frame_interval_fs(10.0, "a.u.", 4) == pytest.approx(40 * mc.ATU_FS)
        assert mc.frame_interval_fs(0.5, "fs", 20) == pytest.approx(10.0)
        with pytest.raises(ValueError, match="time step unit"):
            mc.frame_interval_fs(1.0, "ps", 1)

    def test_au_interval_matches_the_old_time_axis_exactly(self):
        """The ps time axis of an a.u. run is bit-identical to the pre-unit code."""
        atufs, step, sampling = 0.02418884326505, 7.0, 13
        assert mc.frame_interval_fs(step, "a.u.", sampling) / 1000.0 == (atufs * step * sampling) / 1000.0

    def test_frames_must_be_whole_and_non_negative(self):
        assert mc.tolerance_in_frames(3, "frames", 0.5) == 3
        assert mc.tolerance_in_frames(4.0, "frames", 0.5) == 4
        with pytest.raises(ValueError, match="whole number"):
            mc.tolerance_in_frames(2.5, "frames", 0.5)
        for bad in (-1, float("nan"), float("inf")):
            with pytest.raises(ValueError):
                mc.tolerance_in_frames(bad, "frames", 0.5)
        with pytest.raises(ValueError, match="unit"):
            mc.tolerance_in_frames(1, "ps", 0.5)

    def test_fs_is_rounded_down(self):
        dt = 0.9675537306
        assert mc.tolerance_in_frames(0.0, "fs", dt) == 0
        assert mc.tolerance_in_frames(0.99 * dt, "fs", dt) == 0   # shorter than one frame
        assert mc.tolerance_in_frames(2.5 * dt, "fs", dt) == 2
        assert mc.tolerance_in_frames(2.999 * dt, "fs", dt) == 2

    @pytest.mark.parametrize("dt", [0.1, 0.3, 0.9675537306, 40 * mc.ATU_FS, 2.5, 1e-3])
    def test_n_frames_and_n_dt_fs_agree(self, dt):
        """Exact multiples survive floating-point division (0.3 / 0.1 = 2.9999999999999996)."""
        for n in range(0, 200):
            assert mc.tolerance_in_frames(n * dt, "fs", dt) == mc.tolerance_in_frames(n, "frames", dt) == n

    def test_fs_needs_a_positive_interval(self):
        with pytest.raises(ValueError, match="greater than zero"):
            mc.tolerance_in_frames(5.0, "fs", 0.0)


# ------------------------------------------------------------------
# Exits
# ------------------------------------------------------------------
class TestExits:
    def test_zero_tolerance_counts_every_exit(self):
        assert MRTAnalyzer.count_exits(PATTERN, 0) == 3
        assert MRTAnalyzer.count_exits([1, 0, 1, 0, 1, 0], 0) == 3   # the last one, still outside, too

    @pytest.mark.parametrize("tolerance, exits", [(0, 3), (1, 2), (2, 1), (3, 0), (10, 0)])
    def test_bridged_gaps_are_not_exits(self, tolerance, exits):
        corrected = MRTAnalyzer.apply_tolerance_to_occupancy(np.array([PATTERN]), tolerance)[0]
        assert MRTAnalyzer.count_exits(corrected, tolerance) == exits

    def test_an_excursion_at_the_end_counts_only_beyond_the_tolerance(self):
        row = [1, 1, 0, 0]
        assert MRTAnalyzer.count_exits(row, 0) == 1
        assert MRTAnalyzer.count_exits(row, 1) == 1     # already 2 frames out
        assert MRTAnalyzer.count_exits(row, 2) == 0     # not (yet) more than the tolerance

    def test_entering_late_or_never_leaving_is_no_exit(self):
        assert MRTAnalyzer.count_exits([0, 0, 1, 1], 0) == 0
        assert MRTAnalyzer.count_exits([1, 1, 1], 0) == 0
        assert MRTAnalyzer.count_exits([], 0) == 0


# ------------------------------------------------------------------
# The tool, headless
# ------------------------------------------------------------------
def write_trajectory(path, pattern=PATTERN):
    """Na at the origin, O1 at 1.5 A (inside) or 4 A (outside), O2 always inside."""
    lines = []
    for inside in pattern:
        lines += ["3", "frame", "Na 0.0 0.0 0.0", f"O {1.5 if inside else 4.0} 0.0 0.0", "O 0.0 1.0 0.0"]
    path.write_text("\n".join(lines) + "\n")
    return path


class RecordedDialog:
    def __init__(self, title, message):
        self.title, self.message = title, message


@pytest.fixture(autouse=True)
def recorded_dialogs(monkeypatch):
    monkeypatch.setattr(mc.toga, "InfoDialog", RecordedDialog)


class StubWindow:
    def __init__(self):
        self.dialogs = []

    async def dialog(self, d):
        self.dialogs.append(d)


def make_tool(tmp_path, tolerance="", tolerance_unit="frames", time_step="10", time_step_unit="a.u.",
              sampling="4"):
    trajectory = write_trajectory(tmp_path / "TRAJEC.xyz")
    tool = ClassicalMeanResidenceTimeUI.__new__(ClassicalMeanResidenceTimeUI)
    tool.trajec, tool.output_dir = str(trajectory), str(tmp_path)
    tool.total_frame_number, tool.num_atoms = len(PATTERN), 3
    tool.main_window = StubWindow()
    w = lambda value: SimpleNamespace(value=value)   # noqa: E731
    tool.textInput_outer_radius, tool.textInput_center_label = w("2.5"), w("1")
    tool.textInput_atom_symbol, tool.textInput_atom_labels = w("O"), w("0")
    tool.textInput_time_step, tool.textInput_sampling = w(time_step), w(sampling)
    tool.selection_time_step_unit = w(time_step_unit)
    tool.textInput_tolerance, tool.selection_tolerance_unit = w(tolerance), w(tolerance_unit)
    tool.multi_line_text = w("")
    return tool


def run(tool):
    assert asyncio.run(tool.read_params(None)), tool.main_window.dialogs[-1].message
    tool.calc_mrt()
    tool.mrt_summary()
    out = Path(tool.output_dir)
    return {name: (out / name).read_text() for name in ("mrt.dat", "mrt_summary.dat")}


def summary_values(text, key):
    return [line.split()[-1] for line in text.splitlines() if f" {key}." in line]


class TestTool:
    def test_default_tolerance_is_zero_and_counts_every_exit(self, tmp_path):
        tool = make_tool(tmp_path)          # blank tolerance field
        files = run(tool)
        assert tool.tolerance_frames == 0
        assert summary_values(files["mrt_summary.dat"], "exits") == ["3", "0"]   # moving O1, static O2
        assert "every exit from the shell is counted" in tool.multi_line_text.value

    def test_mrt_divides_by_exits(self, tmp_path):
        tool = make_tool(tmp_path, tolerance="0")
        run(tool)
        dt_ps = 40 * mc.ATU_FS / 1000.0
        label, total, mrt, exits = tool.mrt_results[0]
        assert (label, exits) == (2, 3)
        assert total == pytest.approx(PATTERN.count(1) * dt_ps)
        assert mrt == pytest.approx(total / 3)
        o2 = tool.mrt_results[1]
        assert o2[3] == 0 and o2[2] == pytest.approx(o2[1])   # never left: mrt = total rt

    @pytest.mark.parametrize("frames", [0, 1, 2, 3])
    def test_frames_and_equivalent_fs_give_identical_files(self, tmp_path, frames):
        dt = mc.frame_interval_fs(10.0, "a.u.", 4)
        (tmp_path / "frames").mkdir()
        (tmp_path / "fs").mkdir()
        by_frames = run(make_tool(tmp_path / "frames", tolerance=str(frames)))
        by_fs = run(make_tool(tmp_path / "fs", tolerance=repr(frames * dt), tolerance_unit="fs"))
        assert by_fs == by_frames                              # mrt.dat and mrt_summary.dat, byte for byte

    def test_a_time_between_multiples_acts_as_the_lower_one(self, tmp_path):
        dt = mc.frame_interval_fs(10.0, "a.u.", 4)
        tool = make_tool(tmp_path, tolerance=repr(2.7 * dt), tolerance_unit="fs")
        files = run(tool)
        assert tool.tolerance_frames == 2
        assert summary_values(files["mrt_summary.dat"], "exits")[0] == "1"
        assert f"{2.7 * dt:g} fs -> 2 frames" in tool.multi_line_text.value

    def test_time_step_in_fs_matches_the_same_step_in_au(self, tmp_path):
        (tmp_path / "au").mkdir()
        (tmp_path / "fs").mkdir()
        au = make_tool(tmp_path / "au", tolerance="25", tolerance_unit="fs")
        fs = make_tool(tmp_path / "fs", tolerance="25", tolerance_unit="fs",
                       time_step=repr(10 * mc.ATU_FS), time_step_unit="fs")
        run(au)
        run(fs)
        assert au.tolerance_frames == fs.tolerance_frames == math.floor(25 / au.frame_interval_fs) == 25
        assert fs.frame_interval_fs == pytest.approx(au.frame_interval_fs)

    def test_summary_header_records_the_tolerance(self, tmp_path):
        files = run(make_tool(tmp_path, tolerance="2"))
        header = files["mrt_summary.dat"]
        assert "Tolerance (frames)" in header and "Tolerance (fs)" in header
        assert "Time between frames (fs)" in header

    @pytest.mark.parametrize("text, unit, message", [
        ("2.5", "frames", "whole number"),
        ("-1", "fs", "zero or a positive"),
        ("abc", "frames", "not a number"),
    ])
    def test_bad_tolerances_are_refused_with_a_message(self, tmp_path, text, unit, message):
        tool = make_tool(tmp_path, tolerance=text, tolerance_unit=unit)
        assert asyncio.run(tool.read_params(None)) is False
        assert message in tool.main_window.dialogs[-1].message

    def test_the_tolerance_label_follows_its_unit(self):
        tool = ClassicalMeanResidenceTimeUI.__new__(ClassicalMeanResidenceTimeUI)
        tool.tolerance_label = SimpleNamespace(text="")
        tool.textInput_tolerance = SimpleNamespace(placeholder="")
        tool.selection_tolerance_unit = SimpleNamespace(value="fs")
        tool.on_tolerance_unit_change(None)
        assert tool.tolerance_label.text == "Tolerance (fs):" and "0 = count every exit" in \
            tool.textInput_tolerance.placeholder
        tool.selection_tolerance_unit.value = "frames"
        tool.on_tolerance_unit_change(None)
        assert tool.tolerance_label.text == "Tolerance (frames):"

    def test_the_help_text_documents_the_rules(self):
        from help import HelpGqteaWin
        text = HelpGqteaWin.help_mrt
        for phrase in ("ROUNDED DOWN", "dt (fs) = time step (fs) x sampling interval", "The default is 0",
                       "exits"):
            assert phrase in text


# ------------------------------------------------------------------
# Golden outputs: the one-pass reader reproduces the two-pass loop
# ------------------------------------------------------------------
REFERENCE_DIR = Path(__file__).parent / "data" / "mrt_classical"
_spec = importlib.util.spec_from_file_location("mrt_make_reference", REFERENCE_DIR / "make_reference.py")
make_reference = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(make_reference)


class TestGoldenOutputs:
    """The reference files were written by the two-pass reader of commit c1e50e3."""

    @pytest.mark.parametrize("case", sorted(make_reference.CASES))
    def test_output_files_are_byte_identical(self, tmp_path, case):
        tool = make_reference.run_case(tmp_path, *make_reference.CASES[case])
        for name in make_reference.OUTPUTS:
            assert (tmp_path / name).read_bytes() == (REFERENCE_DIR / case / name).read_bytes(), name
        assert tool.center_atom_symbol == ("Na" if make_reference.CASES[case][0] == 1 else "O")

    def test_columns_follow_the_order_of_first_entry(self):
        header = (REFERENCE_DIR / "o_tol0" / "mrt.dat").read_text().splitlines()[0].split()[2:]
        labels = [int(h[1:]) for h in header]
        assert labels != sorted(labels)     # the fixture really exercises the ordering

    @pytest.mark.parametrize("atom, frame", [(2, 0), (5, 40), (8, 150), (14, 299), (29, 77)])
    def test_an_atom_exactly_at_the_cutoff_is_inside(self, atom, frame):
        """The distance is the float64 the per-pair np.linalg.norm loop computed, and the test is <=."""
        text = make_reference.TRAJECTORY.read_text().splitlines()
        stride = 1 + 3 * make_reference.N_WATERS + 2

        def xyz(label):
            return np.array([float(v) for v in text[frame * stride + 1 + label].split()[1:]])

        d = np.linalg.norm(xyz(atom) - xyz(1))
        path = str(make_reference.TRAJECTORY)
        at = mc.scan_trajectory(path, 1, "O", [], d)
        assert at.occupancy[at.shell_atom_list.index(atom)][frame] == 1
        below = mc.scan_trajectory(path, 1, "O", [], np.nextafter(d, 0))
        assert atom not in below.shell_atom_list or below.occupancy[below.shell_atom_list.index(atom)][frame] == 0


# ------------------------------------------------------------------
# Vectorized tolerance and exits == the original loops
# ------------------------------------------------------------------
def loop_fill_gaps(occupancy, tolerance_frames):
    """apply_tolerance_to_occupancy as it was before vectorization."""
    corrected = occupancy.copy()
    if tolerance_frames <= 0:
        return corrected
    n_objects, n_frames = corrected.shape
    for obj in range(n_objects):
        arr = corrected[obj]
        i = 0
        while i < n_frames:
            if arr[i] == 1:
                i += 1
                continue
            start_zero = i
            while i < n_frames and arr[i] == 0:
                i += 1
            left_is_one = start_zero - 1 >= 0 and arr[start_zero - 1] == 1
            right_is_one = i < n_frames and arr[i] == 1
            if left_is_one and right_is_one and i - start_zero <= tolerance_frames:
                arr[start_zero:i] = 1
    return corrected


def loop_count_exits(row, tolerance_frames):
    """count_exits as it was before vectorization."""
    row = [int(v) for v in row]
    n, exits, i = len(row), 0, 0
    while i < n - 1:
        if row[i] == 1 and row[i + 1] == 0:
            j = i + 1
            while j < n and row[j] == 0:
                j += 1
            if j < n or (j - i - 1) > tolerance_frames:
                exits += 1
            i = j
        else:
            i += 1
    return exits


class TestVectorizedMatchesLoops:
    @pytest.mark.parametrize("n_frames", [1, 2, 3, 7, 50, 400])
    @pytest.mark.parametrize("p_inside", [0.0, 0.2, 0.5, 0.9, 1.0])
    def test_random_occupancies(self, n_frames, p_inside):
        rng = np.random.default_rng(n_frames * 100 + int(p_inside * 10))
        occupancy = (rng.random((40, n_frames)) < p_inside).astype(np.int8)
        for tolerance in range(0, 8):
            filled = MRTAnalyzer.apply_tolerance_to_occupancy(occupancy, tolerance)
            assert np.array_equal(filled, loop_fill_gaps(occupancy, tolerance))
            assert filled.dtype == occupancy.dtype
            expected = [loop_count_exits(row, tolerance) for row in filled]
            assert mc.exits_per_atom(filled, tolerance).tolist() == expected
            assert [MRTAnalyzer.count_exits(row, tolerance) for row in filled] == expected

    def test_long_runs(self):
        rng = np.random.default_rng(3)
        occupancy = np.repeat(rng.random((12, 300)) < 0.5, rng.integers(1, 9, 300), axis=1).astype(int)
        for tolerance in (0, 3, 10, 40):
            filled = MRTAnalyzer.apply_tolerance_to_occupancy(occupancy, tolerance)
            assert np.array_equal(filled, loop_fill_gaps(occupancy, tolerance))
            assert mc.exits_per_atom(filled, tolerance).tolist() == [loop_count_exits(r, tolerance) for r in filled]

    def test_the_input_is_not_modified(self):
        occupancy = np.array([PATTERN])
        before = occupancy.copy()
        MRTAnalyzer.apply_tolerance_to_occupancy(occupancy, 5)
        assert np.array_equal(occupancy, before)

    def test_no_atoms(self):
        empty = np.zeros((0, 12), np.int8)
        assert MRTAnalyzer.apply_tolerance_to_occupancy(empty, 3).shape == (0, 12)
        assert mc.exits_per_atom(empty, 0).shape == (0,)


# ------------------------------------------------------------------
# Counting frames and reading malformed trajectories
# ------------------------------------------------------------------
def trajectory_lines(pattern=PATTERN):
    """The lines write_trajectory writes."""
    lines = []
    for inside in pattern:
        lines += ["3", "frame", "Na 0.0 0.0 0.0", f"O {1.5 if inside else 4.0} 0.0 0.0", "O 0.0 1.0 0.0"]
    return lines


def write_lines(path, lines, end="\n"):
    path.write_bytes((end.join(lines) + end).encode())
    return str(path)


def scan(path, n_frames=len(PATTERN), exclude=()):
    return mc.scan_trajectory(path, 1, "O", list(exclude), 2.5, n_frames=n_frames)


class TestCountFrames:
    def test_fixture(self):
        assert mc.count_frames(str(make_reference.TRAJECTORY)) == (31, 300)

    @pytest.mark.parametrize("end", ["\n", "\r\n"])
    @pytest.mark.parametrize("tail", ["", "\n", "\n\n  \n"])
    def test_line_endings_and_trailing_blank_lines(self, tmp_path, end, tail):
        path = tmp_path / "t.xyz"
        path.write_bytes((end.join(trajectory_lines()) + end + tail.replace("\n", end)).encode())
        assert mc.count_frames(str(path)) == (3, len(PATTERN))

    def test_no_final_newline(self, tmp_path):
        path = tmp_path / "t.xyz"
        path.write_bytes("\n".join(trajectory_lines()).encode())
        assert mc.count_frames(str(path)) == (3, len(PATTERN))

    def test_progress_reaches_one(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mc, "COUNT_BLOCK_BYTES", 64)
        seen = []
        mc.count_frames(write_lines(tmp_path / "t.xyz", trajectory_lines()), seen.append)
        assert len(seen) > 1 and seen == sorted(seen) and seen[-1] == 1.0

    @pytest.mark.parametrize("lines, message", [
        (trajectory_lines()[:-1], "not a whole number"),                                # truncated
        (trajectory_lines()[:5] + [""] + trajectory_lines()[5:], "not a whole number"),  # blank line inside
        ([""], "empty"),
        (["abc", "x"], "number of atoms"),
    ])
    def test_bad_files(self, tmp_path, lines, message):
        with pytest.raises(ValueError, match=message):
            mc.count_frames(write_lines(tmp_path / "t.xyz", lines))


class TestScanTrajectory:
    def test_trailing_blank_lines_are_ignored(self, tmp_path):
        result = scan(write_lines(tmp_path / "t.xyz", trajectory_lines() + ["", "  "]))
        assert result.shell_atom_list == [2, 3]
        assert result.occupancy[0].tolist() == PATTERN

    @pytest.mark.parametrize("chunk", [1, 2, 7, 100])
    def test_block_size_does_not_change_the_result(self, chunk):
        args = (str(make_reference.TRAJECTORY), 1, "H", [], 3.0)
        a = mc.scan_trajectory(*args, chunk_frames=chunk)
        b = mc.scan_trajectory(*args)
        assert a.shell_atom_list == b.shell_atom_list and np.array_equal(a.occupancy, b.occupancy)

    def test_the_frames_are_counted_when_not_given(self, tmp_path):
        assert scan(write_lines(tmp_path / "t.xyz", trajectory_lines()), n_frames=None).occupancy.shape == (2, 12)

    def test_a_changed_atom_count_is_refused(self, tmp_path):
        lines = trajectory_lines()
        lines[10] = "4"         # header of frame 3
        with pytest.raises(ValueError, match="Frame 3 should start with the atom count 3"):
            scan(write_lines(tmp_path / "t.xyz", lines))

    def test_a_short_file_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="Invalid atom record found in frame 12"):
            scan(write_lines(tmp_path / "t.xyz", trajectory_lines()[:-1]))

    def test_frames_beyond_the_count_are_refused(self, tmp_path):
        with pytest.raises(ValueError, match="load it again with Browse"):
            scan(write_lines(tmp_path / "t.xyz", trajectory_lines()), n_frames=10)

    @pytest.mark.parametrize("index", [3, 23])     # an analysed O in frame 1 / frame 5
    def test_a_malformed_line_of_an_analysed_atom_is_refused(self, tmp_path, index):
        lines = trajectory_lines()
        lines[index] = "O 1.5 0.0"
        with pytest.raises(ValueError, match=f"Invalid atom record found in frame {index // 5 + 1}"):
            scan(write_lines(tmp_path / "t.xyz", lines))

    def test_a_malformed_line_of_an_atom_not_analysed_is_not_read(self, tmp_path):
        lines = trajectory_lines()
        lines[24] = "O garbage"          # atom 3, frame 5; atom 3 is excluded
        result = scan(write_lines(tmp_path / "t.xyz", lines), exclude=[3])
        assert result.shell_atom_list == [2]

    def test_frame_1_is_checked_line_by_line(self, tmp_path):
        lines = trajectory_lines()
        lines[4] = "O garbage"          # atom 3, frame 1, even though it is excluded
        with pytest.raises(ValueError, match="frame 1"):
            scan(write_lines(tmp_path / "t.xyz", lines), exclude=[3])

    def test_a_changed_atom_order_is_refused(self, tmp_path):
        lines = trajectory_lines()
        lines[18] = "H 1.5 0.0 0.0"     # atom 2 of frame 4
        with pytest.raises(ValueError, match="frame 4 atom 2 is H"):
            scan(write_lines(tmp_path / "t.xyz", lines))

    def test_a_bad_number_is_refused(self, tmp_path):
        lines = trajectory_lines()
        lines[13] = "O 1.5 0.0 x"
        with pytest.raises(ValueError, match="could not convert"):
            scan(write_lines(tmp_path / "t.xyz", lines))


# ------------------------------------------------------------------
# The window: worker thread, spinner, status line
# ------------------------------------------------------------------
class RecordingLabel:
    def __init__(self):
        self.history = []

    @property
    def text(self):
        return self.history[-1] if self.history else ""

    @text.setter
    def text(self, value):
        self.history.append(value)


class StubSpinner:
    def __init__(self):
        self.running, self.events = False, []

    def start(self):
        self.running = True
        self.events.append("start")

    def stop(self):
        self.running = False
        self.events.append("stop")


def make_window_tool(tmp_path, lines=None):
    tool = make_tool(tmp_path)
    if lines is not None:
        write_lines(tmp_path / "TRAJEC.xyz", lines)
    tool.activity_indicator, tool.status_label = StubSpinner(), RecordingLabel()
    tool.browse_button, tool.btn_execute = SimpleNamespace(enabled=True), SimpleNamespace(enabled=True)
    tool.plots_shown = 0
    tool.save_plots = lambda *args, **kwargs: None

    def show():
        tool.plots_shown += 1
    tool.display_plots = show
    return tool


class TestWindow:
    def test_the_analysis_runs_in_a_worker_thread_with_the_spinner(self, tmp_path, monkeypatch):
        monkeypatch.setattr(mc, "STATUS_INTERVAL", 0.0)
        tool = make_window_tool(tmp_path)
        threads, original = [], tool.run_analysis

        def run_analysis(progress=None):
            threads.append(threading.get_ident())
            assert tool.activity_indicator.running and not tool.btn_execute.enabled
            original(progress)
        tool.run_analysis = run_analysis
        asyncio.run(tool.workflow(None))
        assert threads and threads[0] != threading.get_ident()
        assert tool.activity_indicator.events == ["start", "stop"]
        assert tool.btn_execute.enabled and tool.browse_button.enabled
        history = tool.status_label.history
        assert "Reading the trajectory... frame 12 of 12" in history
        assert "Applying the tolerance..." in history
        assert history[-1].startswith("Finished in") and "2 O atoms entered the shell" in history[-1]
        assert "FINISHED SUCCESSFULLY" in tool.multi_line_text.value
        assert tool.plots_shown == 1
        assert (tmp_path / "mrt_total.dat").exists()

    def test_an_error_stops_the_spinner_and_is_shown(self, tmp_path):
        lines = trajectory_lines()
        lines[13] = "O 1.5"
        tool = make_window_tool(tmp_path, lines)
        asyncio.run(tool.workflow(None))
        assert tool.activity_indicator.events == ["start", "stop"]
        assert tool.btn_execute.enabled and tool.browse_button.enabled
        assert tool.status_label.text == "The MRT calculation stopped with an error."
        assert "Invalid atom record found in frame 3" in tool.main_window.dialogs[-1].message
        assert tool.plots_shown == 0

    def test_a_failing_progress_callback_does_not_stop_the_analysis(self, tmp_path):
        tool = make_window_tool(tmp_path)
        assert asyncio.run(tool.read_params(None))

        def broken(message, force):
            raise RuntimeError("status line gone")
        tool.run_analysis(broken)
        assert (tmp_path / "mrt_summary.dat").exists()

    def test_a_closed_window_is_not_touched(self, tmp_path):
        tool = make_window_tool(tmp_path)
        tool._on_close(None)
        asyncio.run(tool.workflow(None))
        assert tool.activity_indicator.events == [] and tool.plots_shown == 0
        assert (tmp_path / "mrt.dat").exists()      # the job itself still finished

    @staticmethod
    def browse(tool, path):
        async def choose(widget):
            tool.trajec = str(path)
        tool.open_file_dialog = choose
        asyncio.run(tool.frames_counter(None))

    def test_browse_counts_frames_with_the_spinner(self, tmp_path):
        tool = make_window_tool(tmp_path)
        tool.total_frame_number = tool.num_atoms = None
        self.browse(tool, tmp_path / "TRAJEC.xyz")
        assert (tool.num_atoms, tool.total_frame_number) == (3, 12)
        assert tool.activity_indicator.events == ["start", "stop"]
        assert tool.browse_button.enabled
        assert "Number of frames -->  12" in tool.multi_line_text.value
        assert "Number of lines  -->  60" in tool.multi_line_text.value
        assert tool.status_label.text.startswith("12 frames of 3 atoms")

    def test_browse_of_a_bad_file_leaves_no_frame_count(self, tmp_path):
        tool = make_window_tool(tmp_path, trajectory_lines()[:-1])
        self.browse(tool, tmp_path / "TRAJEC.xyz")
        assert tool.total_frame_number == 0
        assert "not a whole number" in tool.main_window.dialogs[-1].message
        assert tool.activity_indicator.events == ["start", "stop"] and tool.browse_button.enabled
