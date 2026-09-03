"""Tests for the Radial Distribution Function tool (``radialDistribution.py``).

Run from ``venv/src/``::

    python -m pytest tests/ -q

The Toga GUI is never instantiated: dialogs are routed through a stub and the
interactive plot viewer is replaced with a no-op, so everything runs headless.
Physics is checked against an analytic reference (a simple-cubic lattice) rather
than a frozen golden file, so the tests stay meaningful if the implementation is
refactored again.
"""
import os
import glob
import asyncio

import numpy as np
import pytest

import radialDistribution as rd
from radialDistribution import RadialAnalyser
from displayPlots import DisplayPlots
import plotViewer


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


class FakeInput:
    def __init__(self, value=""):
        self.value = value


@pytest.fixture(autouse=True)
def _patch_toga(monkeypatch):
    """Route every dialog through a capturable stub instead of the GUI backend."""
    monkeypatch.setattr(rd.toga, "InfoDialog", FakeDialog)


# --------------------------------------------------------------------------- #
# Fixture builders                                                              #
# --------------------------------------------------------------------------- #
def write_xyz(path, frames, symbols):
    with open(path, "w") as f:
        for coords in frames:
            f.write(f"{len(symbols)}\ncomment\n")
            for s, (x, y, z) in zip(symbols, coords):
                f.write(f"{s} {x:.6f} {y:.6f} {z:.6f}\n")


def simple_cubic(spacing=3.0, reps=3):
    """Simple-cubic lattice of identical atoms -> (coords, symbols, box_length)."""
    pts = [(i * spacing, j * spacing, k * spacing)
           for i in range(reps) for j in range(reps) for k in range(reps)]
    coords = np.array(pts, dtype=float)
    symbols = ["O"] * len(pts)
    return coords, symbols, spacing * reps


def make_analyser(out_dir, traj, num_atoms, n_frames, *, box="10 10 10",
                  radius="4", bin_width="0.2", symbol="O", shell_center="1",
                  atom_list="0", axis="4 5"):
    a = RadialAnalyser()
    a.main_window = FakeWindow()
    a.multi_line_text = FakeInput()
    a.output_dir = str(out_dir)
    a.display_plots = lambda *args, **kw: None  # never spawn the GUI viewer
    a.textInput_radius = FakeInput(radius)
    a.textInput_bin_width = FakeInput(bin_width)
    a.textInput_atom_list = FakeInput(atom_list)
    a.textInput_shell_center = FakeInput(shell_center)
    a.textInput_atom_symbol = FakeInput(symbol)
    a.textInput_cell_lattices = FakeInput(box)
    a.textInput_axis = FakeInput(axis)
    a.trajec = str(traj)
    a.num_atoms = num_atoms
    a.total_frame_number = n_frames
    return a


def read_dat(out_dir):
    dat = glob.glob(os.path.join(str(out_dir), "*.dat"))[0]
    return np.loadtxt(dat, skiprows=2)  # columns: r, g(r), Integral


@pytest.fixture
def tiny_traj(tmp_path):
    """A 4-atom, 2-frame trajectory (atoms 1-2 are O, 3-4 are H)."""
    traj = tmp_path / "tiny.xyz"
    coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], float)
    write_xyz(traj, [coords, coords + 0.1], ["O", "O", "H", "H"])
    return traj


# --------------------------------------------------------------------------- #
# Physics: minimum-image RDF against an analytic reference                      #
# --------------------------------------------------------------------------- #
def test_simple_cubic_coordination(tmp_path):
    """Simple cubic (s=3, box=9): 6 neighbours at 3.0, 12 at 3*sqrt(2)=4.243."""
    coords, symbols, box = simple_cubic()
    traj = tmp_path / "sc.xyz"
    write_xyz(traj, [coords], symbols)
    a = make_analyser(tmp_path, traj, len(symbols), 1,
                      box=f"{box} {box} {box}", radius="4.4", bin_width="0.1")
    assert asyncio.run(a.read_params(None)) is True
    asyncio.run(a.calc_rdf())

    coordination = read_dat(tmp_path)[:, 2]
    assert coordination[-1] == pytest.approx(18.0)  # 6 + 12 within r < 4.4
    assert 6.0 in coordination                       # first-shell plateau
    assert np.all(np.diff(coordination) >= -1e-12)   # monotonic non-decreasing

    # The completion summary is written to the status box.
    status = a.multi_line_text.value
    assert "COMPLETED" in status
    assert "Output: RDF_O1_O.dat" in status
    assert "Coordination" in status


def test_pbc_shift_invariance(tmp_path):
    """Shifting every atom by one full box (unwrapped coords) must not change g(r).

    This is the regression guard for the minimum-image / periodic-boundary fix.
    """
    coords, symbols, box = simple_cubic()
    base, shifted = tmp_path / "base.xyz", tmp_path / "shifted.xyz"
    write_xyz(base, [coords], symbols)
    write_xyz(shifted, [coords + np.array([box, 0.0, 0.0])], symbols)

    out1, out2 = tmp_path / "o1", tmp_path / "o2"
    out1.mkdir()
    out2.mkdir()
    a1 = make_analyser(out1, base, len(symbols), 1, box=f"{box} {box} {box}", radius="4.4")
    a2 = make_analyser(out2, shifted, len(symbols), 1, box=f"{box} {box} {box}", radius="4.4")
    assert asyncio.run(a1.read_params(None)) is True
    assert asyncio.run(a2.read_params(None)) is True
    asyncio.run(a1.calc_rdf())
    asyncio.run(a2.calc_rdf())

    d1, d2 = read_dat(out1), read_dat(out2)
    np.testing.assert_allclose(d1[:, 1], d2[:, 1], atol=1e-9)  # g(r)
    np.testing.assert_allclose(d1[:, 2], d2[:, 2], atol=1e-9)  # Integral


def test_r_is_bin_centers(tmp_path):
    coords, symbols, box = simple_cubic()
    traj = tmp_path / "sc.xyz"
    write_xyz(traj, [coords], symbols)
    a = make_analyser(tmp_path, traj, len(symbols), 1,
                      box=f"{box} {box} {box}", radius="4.0", bin_width="0.2")
    assert asyncio.run(a.read_params(None)) is True
    asyncio.run(a.calc_rdf())

    r = read_dat(tmp_path)[:, 0]
    assert r[0] == pytest.approx(0.1)          # centre of the [0.0, 0.2] bin
    assert np.allclose(np.diff(r), 0.2)


# --------------------------------------------------------------------------- #
# Parameter validation (read_params -> bool)                                    #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("overrides,expected", [
    ({}, True),
    ({"radius": "-1"}, False),
    ({"radius": "abc"}, False),
    ({"bin_width": "0"}, False),
    ({"bin_width": "99"}, False),          # bin width > radius
    ({"atom_list": "9"}, False),           # label out of range
    ({"atom_list": "0 3"}, False),         # 0 mixed with real labels
    ({"shell_center": "0"}, False),
    ({"shell_center": "99"}, False),
    ({"box": "10 10"}, False),             # not exactly three lattices
    ({"box": "10 0 10"}, False),           # non-positive lattice
    ({"radius": "6"}, False),              # > min(a,b,c)/2 == 5
    ({"axis": "4"}, False),                # need two values
    ({"axis": "0 5"}, False),              # non-positive limit
    ({"symbol": "Xx"}, False),             # target symbol absent -> rho == 0
])
def test_read_params_validation(tmp_path, tiny_traj, overrides, expected):
    params = dict(box="10 10 10", radius="4", bin_width="0.2", axis="4 5")
    params.update(overrides)
    a = make_analyser(tmp_path, tiny_traj, 4, 2, **params)
    assert asyncio.run(a.read_params(None)) is expected


def test_read_params_requires_loaded_file(tmp_path, tiny_traj):
    a = make_analyser(tmp_path, tiny_traj, 4, 2)
    del a.trajec  # simulate "user never clicked Browse"
    assert asyncio.run(a.read_params(None)) is False


# --------------------------------------------------------------------------- #
# Frame handling: clean EOF vs malformed/truncated frames (C8)                  #
# --------------------------------------------------------------------------- #
def _run_calc(out_dir, traj, n_frames):
    a = make_analyser(out_dir, traj, 4, n_frames)
    assert asyncio.run(a.read_params(None)) is True
    a.main_window.dialogs.clear()
    asyncio.run(a.calc_rdf())
    warns = [m for t, m in a.main_window.dialogs if t == "Warning"]
    errs = [m for t, m in a.main_window.dialogs if t == "Error"]
    return warns, errs


def test_calc_rdf_clean_no_warning(tmp_path):
    traj = tmp_path / "clean.xyz"
    c = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], float)
    write_xyz(traj, [c, c, c], ["O", "O", "H", "H"])
    warns, errs = _run_calc(tmp_path, traj, 3)
    assert not warns and not errs


def test_calc_rdf_trailing_blank_no_warning(tmp_path):
    traj = tmp_path / "blank.xyz"
    with open(traj, "w") as f:
        f.write("4\nc\nO 0 0 0\nO 1 0 0\nH 2 0 0\nH 3 0 0\n")
        f.write("4\nc\nO 0 0 0\nO 1 0 0\nH 2 0 0\nH 3 0 0\n")
        f.write("\n")  # trailing blank line == clean end of data
    warns, errs = _run_calc(tmp_path, traj, 2)
    assert not warns and not errs


def test_calc_rdf_truncated_frame_warns(tmp_path):
    traj = tmp_path / "trunc.xyz"
    with open(traj, "w") as f:
        f.write("4\nc\nO 0 0 0\nO 1 0 0\nH 2 0 0\nH 3 0 0\n")
        f.write("4\nc\nO 0 0 0\nO 1 0 0\n")  # declares 4, provides 2
    warns, errs = _run_calc(tmp_path, traj, 2)
    assert warns and not errs


def test_calc_rdf_garbage_header_warns(tmp_path):
    traj = tmp_path / "garbage.xyz"
    with open(traj, "w") as f:
        f.write("4\nc\nO 0 0 0\nO 1 0 0\nH 2 0 0\nH 3 0 0\n")
        f.write("garbage\n")  # not an atom count
    warns, errs = _run_calc(tmp_path, traj, 2)
    assert warns and not errs


# --------------------------------------------------------------------------- #
# Shared plotting extension (optional xlim/ylim), backward compatibility        #
# --------------------------------------------------------------------------- #
def test_save_plots_limits_roundtrip(tmp_path):
    class _Plotter(DisplayPlots):
        pass

    t = _Plotter()
    t.output_dir = str(tmp_path)
    t.saved_plot_files = []
    t.saved_plot_data = []

    t.save_plots(1, [0, 1, 2], [0, 1, 4], "x", "y", "no limits")
    t.save_plots(2, [0, 1, 2], [0, 1, 4], "x", "y", "limited",
                 xlim=(0, 5), ylim=(0, 30))

    e_plain, e_limited = t.saved_plot_data
    assert "xlim" not in e_plain and "ylim" not in e_plain     # backward compatible
    assert e_limited["xlim"] == [0.0, 5.0]
    assert e_limited["ylim"] == [0.0, 30.0]

    # The standalone viewer must accept both manifest shapes without error.
    plotViewer.build_figures(t.saved_plot_data)


def test_save_plots_save_png_flag(tmp_path):
    """save_png=False records viewer data but writes no PNG; True still does."""
    class _Plotter(DisplayPlots):
        pass

    t = _Plotter()
    t.output_dir = str(tmp_path)
    t.saved_plot_files = []
    t.saved_plot_data = []

    # save_png=False: data recorded for the interactive viewer, no file on disk.
    t.save_plots(1, [0, 1, 2], [0, 1, 4], "x", "y", "no png", save_png=False)
    assert len(t.saved_plot_data) == 1
    assert t.saved_plot_files == []
    assert glob.glob(os.path.join(str(tmp_path), "*.png")) == []

    # Default (save_png=True) still writes a PNG (unchanged for other tools).
    t.save_plots(2, [0, 1, 2], [0, 1, 4], "x", "y", "with png")
    assert len(t.saved_plot_files) == 1
    assert len(glob.glob(os.path.join(str(tmp_path), "*.png"))) == 1


def test_save_multiseries_plot_save_png_flag(tmp_path):
    """save_multiseries_plot honours save_png the same way (overlay figures).

    Used by autocorrelationFunction.py for total + partial VDOS overlays.
    """
    class _Plotter(DisplayPlots):
        pass

    t = _Plotter()
    t.output_dir = str(tmp_path)
    t.saved_plot_files = []
    t.saved_plot_data = []

    series = [([0, 1, 2], [0, 1, 4], "Total"), ([0, 1, 2], [0, 2, 8], "O")]

    # save_png=False: series data recorded for the viewer, no file on disk.
    t.save_multiseries_plot(3, series, "f", "P", "no png", save_png=False)
    assert len(t.saved_plot_data) == 1
    assert "series" in t.saved_plot_data[0]
    assert t.saved_plot_files == []
    assert glob.glob(os.path.join(str(tmp_path), "*.png")) == []

    # Default (save_png=True) still writes a PNG.
    t.save_multiseries_plot(4, series, "f", "P", "with png")
    assert len(t.saved_plot_files) == 1
    assert len(glob.glob(os.path.join(str(tmp_path), "*.png"))) == 1


def test_display_plots_spawns_viewer_source_and_frozen(tmp_path, monkeypatch):
    """The interactive viewer must be spawned in BOTH source and frozen builds.

    Regression guard: a frozen build previously skipped the subprocess entirely
    and fell back to static PNG windows. It must now re-launch the app with the
    --plot-viewer flag instead.
    """
    from displayPlots import DisplayPlots
    import displayPlots as dp

    calls = []
    monkeypatch.setattr(dp.subprocess, "Popen", lambda cmd, *a, **k: calls.append(cmd))

    def fresh_plotter():
        t = DisplayPlots()
        t.output_dir = str(tmp_path)
        t.saved_plot_files = []
        t.saved_plot_data = [{"x": [0, 1], "y": [0, 1], "xlabel": "x",
                              "ylabel": "y", "title": "t"}]
        return t

    # Source build: hand a Python interpreter the plotViewer.py script.
    monkeypatch.setattr(dp.sys, "frozen", False, raising=False)
    fresh_plotter().display_plots()
    assert calls[-1][1].endswith("plotViewer.py")
    assert calls[-1][2].endswith(".json")

    # Frozen build: re-launch this same executable with the --plot-viewer flag.
    monkeypatch.setattr(dp.sys, "frozen", True, raising=False)
    fresh_plotter().display_plots()
    assert calls[-1][1] == "--plot-viewer"
    assert calls[-1][2].endswith(".json")


def test_calc_rdf_writes_no_png(tmp_path):
    """The RDF workflow must leave no left-over PNG image files next to the data."""
    coords, symbols, box = simple_cubic()
    traj = tmp_path / "sc.xyz"
    write_xyz(traj, [coords], symbols)
    a = make_analyser(tmp_path, traj, len(symbols), 1,
                      box=f"{box} {box} {box}", radius="4.4", bin_width="0.1")
    assert asyncio.run(a.read_params(None)) is True
    asyncio.run(a.calc_rdf())

    assert glob.glob(os.path.join(str(tmp_path), "*.png")) == []   # no PNG clutter
    assert len(glob.glob(os.path.join(str(tmp_path), "*.dat"))) == 1  # data still written


# --------------------------------------------------------------------------- #
# Trajectory loading reports status as text (no progress bar / frame label)     #
# --------------------------------------------------------------------------- #
def test_frames_counter_reports_status(tmp_path):
    traj = tmp_path / "load.xyz"
    c = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], float)
    write_xyz(traj, [c, c, c], ["O", "O", "H", "H"])

    a = make_analyser(tmp_path, traj, 4, 0)  # frame count not known yet

    async def fake_open(widget):  # stand in for the GUI file-open dialog
        a.trajec = str(traj)
        a.num_atoms = 4
        a.output_dir = str(tmp_path)

    a.open_file_dialog = fake_open
    asyncio.run(a.frames_counter(None))

    assert a.total_frame_number == 3
    assert "Trajectory loaded" in a.multi_line_text.value
    assert "frames: 3" in a.multi_line_text.value
    assert not hasattr(a, "progress_label")  # no frame label is used
