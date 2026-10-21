"""Regenerate the classical-MRT golden fixture and its reference outputs.

    python tests/data/mrt_classical/make_reference.py        (from venv/src)

``na_water.xyz`` is a synthetic Na+ in water trajectory (atom 1 = Na, then
O, H, H per water) whose waters random-walk in and out of a ~3 A shell.
The reference outputs in ``<case>/`` were first written by the
two-pass reader of commit c1e50e3, before the trajectory reader was
rewritten, and the tests require the current code to reproduce them byte
for byte. Regenerate them only when the output format is *meant* to change.
"""
import asyncio
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
TRAJECTORY = HERE / "na_water.xyz"

# name: (centre label, symbol, cutoff, tolerance in frames, excluded labels)
CASES = {
    "o_tol0": (1, "O", "3.2", "0", "0"),
    "o_tol3_excluded": (1, "O", "3.2", "3", "5 11"),
    "o_centre_o": (2, "O", "3.0", "2", "0"),
    "h_tol1": (1, "H", "3.0", "1", "0"),
}
N_WATERS, N_FRAMES, BOX = 10, 300, 9.0
OUTPUTS = ("mrt.dat", "mrt_summary.dat", "mrt_total.dat")


def write_trajectory(path=TRAJECTORY, seed=7):
    rng = np.random.default_rng(seed)
    oxygens = rng.uniform(0, BOX, (N_WATERS, 3))
    h1, h2 = np.array([0.757, 0.586, 0.0]), np.array([-0.757, 0.586, 0.0])
    n_atoms = 1 + 3 * N_WATERS
    symbols = ["Na"] + ["O", "H", "H"] * N_WATERS
    lines = []
    for frame in range(N_FRAMES):
        oxygens = (oxygens + rng.normal(0, 0.25, oxygens.shape)) % BOX
        pos = np.empty((n_atoms, 3))
        pos[0] = BOX / 2
        pos[1::3], pos[2::3], pos[3::3] = oxygens, oxygens + h1, oxygens + h2
        lines += [f"{n_atoms:5d}", f" STEP: {frame * 10}"]
        lines += [f"{s:<2s} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in zip(symbols, pos)]
    path.write_text("\n".join(lines) + "\n", newline="\n")


def run_case(out_dir, centre, symbol, cutoff, tolerance, excluded):
    """Run the tool headless (stub widgets, no window) and write its three files."""
    sys.path.insert(0, str(HERE.parents[2]))
    from meanResidenceTime_classical import ClassicalMeanResidenceTimeUI

    async def dialog(d):
        raise AssertionError(d.message)

    n_frames = TRAJECTORY.read_text().count(" STEP:")
    tool = ClassicalMeanResidenceTimeUI.__new__(ClassicalMeanResidenceTimeUI)
    tool.trajec, tool.output_dir = str(TRAJECTORY), str(out_dir)
    tool.total_frame_number, tool.num_atoms = n_frames, 1 + 3 * N_WATERS
    tool.main_window = SimpleNamespace(dialog=dialog)
    w =lambda value: SimpleNamespace(value=value)   # noqa: E731
    tool.textInput_outer_radius, tool.textInput_center_label = w(cutoff), w(str(centre))
    tool.textInput_atom_symbol, tool.textInput_atom_labels = w(symbol), w(excluded)
    tool.textInput_time_step, tool.textInput_sampling = w("10"), w("4")
    tool.selection_time_step_unit = w("a.u.")
    tool.textInput_tolerance, tool.selection_tolerance_unit = w(tolerance), w("frames")
    tool.multi_line_text = w("")
    assert asyncio.run(tool.read_params(None))
    tool.save_plots = lambda *args, **kwargs: None    # no PNG next to the references
    tool.run_analysis()
    return tool


if __name__ == "__main__":
    if not TRAJECTORY.exists():
        write_trajectory()
    for name, case in CASES.items():
        out = HERE / name
        out.mkdir(exist_ok=True)
        run_case(out, *case)
        print(name, ", ".join(f"{f} {(out / f).stat().st_size} B" for f in OUTPUTS))
