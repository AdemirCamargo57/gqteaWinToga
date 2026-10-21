"""Classical MRT tool: tolerance in frames or fs, and exit counting.

The tolerance is the largest number of consecutive sampled frames an atom may
spend outside the shell without the excursion counting as an exit. These
tests pin the fs -> frames conversion (rounded down, exact for multiples of
the frame interval), that equivalent frames/fs tolerances give identical
output files, and that a zero tolerance counts every exit. The UI half is
driven headless through ``__new__`` and stub widgets; no window is created.
"""
import asyncio
import math
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


class StubProgress:
    value = 0

    def start(self):
        pass

    def stop(self):
        pass


def make_tool(tmp_path, tolerance="", tolerance_unit="frames", time_step="10", time_step_unit="a.u.",
              sampling="4"):
    trajectory = write_trajectory(tmp_path / "TRAJEC.xyz")
    tool = ClassicalMeanResidenceTimeUI.__new__(ClassicalMeanResidenceTimeUI)
    tool.trajec, tool.output_dir = str(trajectory), str(tmp_path)
    tool.total_frame_number, tool.num_atoms = len(PATTERN), 3
    tool.main_window = StubWindow()
    tool.progress_bar = StubProgress()
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
