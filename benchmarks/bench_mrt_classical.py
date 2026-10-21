"""Time the classical MRT tool against its pre-optimization version and check identical output.

    python benchmarks/bench_mrt_classical.py                     (from venv/src)
    python benchmarks/bench_mrt_classical.py --frames 2000 15121 --symbol H --cutoff 2.6
    python benchmarks/bench_mrt_classical.py --trajectory path/to/TRAJEC.xyz --centre 1 --symbol O

Synthetic Na+ in water trajectories (atom 1 = Na, 143 waters = 430 atoms by
default) are generated once into the system temp dir. The legacy code is the
module as it was at ``--legacy`` (default c1e50e3, the last two-pass reader),
read with ``git show``. Every variant runs in a fresh process, so the first
run includes the cold start (imports, file not yet read by this process);
the second run is in the same process. Times are medians over ``--repeat``
processes; the operating system's file cache is warm after the first.
mrt.dat, mrt_summary.dat and mrt_total.dat must be byte-identical.
"""
import argparse
import asyncio
import ctypes
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace

SRC = Path(__file__).resolve().parents[1]
OUTPUTS = ("mrt.dat", "mrt_summary.dat", "mrt_total.dat")


def generate(path, n_waters, n_frames, seed=1):
    import numpy as np
    rng = np.random.default_rng(seed)
    box = (n_waters / 0.0334) ** (1 / 3)            # liquid-water density
    oxygens = rng.uniform(0, box, (n_waters, 3))
    h1, h2 = np.array([0.757, 0.586, 0.0]), np.array([-0.757, 0.586, 0.0])
    n_atoms = 1 + 3 * n_waters
    symbols = ["Na"] + ["O", "H", "H"] * n_waters
    with open(path, "w", newline="\n") as f:
        for frame in range(n_frames):
            oxygens = (oxygens + rng.normal(0, 0.08, oxygens.shape)) % box
            pos = np.empty((n_atoms, 3))
            pos[0] = box / 2 + rng.normal(0, 0.02, 3)
            pos[1::3], pos[2::3], pos[3::3] = oxygens, oxygens + h1, oxygens + h2
            f.write(f"{n_atoms:6d}\n STEP: {frame * 10}\n")
            f.write("".join(f"{s:<2s} {x:16.8f} {y:16.8f} {z:16.8f}\n" for s, (x, y, z) in zip(symbols, pos)))


def peak_memory_mb():
    """Peak working set of this process (Windows), else None."""
    if os.name != "nt":
        return None

    class Counters(ctypes.Structure):
        _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong)] + [
            (name, ctypes.c_size_t) for name in (
                "PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage", "QuotaPagedPoolUsage",
                "QuotaPeakNonPagedPoolUsage", "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]
    counters = Counters()
    counters.cb = ctypes.sizeof(counters)
    get_process = ctypes.windll.kernel32.GetCurrentProcess
    get_process.restype = ctypes.c_void_p            # a pointer-sized pseudo-handle
    get_info = ctypes.windll.psapi.GetProcessMemoryInfo
    get_info.argtypes = [ctypes.c_void_p, ctypes.POINTER(Counters), ctypes.c_ulong]
    if not get_info(get_process(), ctypes.byref(counters), counters.cb):
        return None
    return counters.PeakWorkingSetSize / 1e6


# ---------------------------------------------------------------- one worker process
def worker(variant, module_dir, trajectory, out_dir, centre, symbol, cutoff, tolerance):
    t0 = time.perf_counter()
    sys.path.insert(0, module_dir)
    sys.path.insert(1, str(SRC))
    import matplotlib
    matplotlib.use("Agg")
    import meanResidenceTime_classical as mc
    t_import = time.perf_counter() - t0

    class Spinner:
        def start(self): pass
        def stop(self): pass

    async def dialog(d):
        raise RuntimeError(d.message)

    tool = mc.ClassicalMeanResidenceTimeUI.__new__(mc.ClassicalMeanResidenceTimeUI)
    w = lambda value: SimpleNamespace(value=value, text="", enabled=True)   # noqa: E731
    tool.main_window = SimpleNamespace(dialog=dialog)
    tool.progress_bar = tool.activity_indicator = Spinner()
    tool.progress_label, tool.status_label = w(""), w("")
    tool.browse_button = tool.btn_execute = w("")
    tool.multi_line_text, tool.textInput_file = w(""), w("")
    tool.textInput_outer_radius, tool.textInput_center_label = w(str(cutoff)), w(str(centre))
    tool.textInput_atom_symbol, tool.textInput_atom_labels = w(symbol), w("0")
    tool.textInput_time_step, tool.textInput_sampling = w("10"), w("4")
    tool.selection_time_step_unit = w("a.u.")
    tool.textInput_tolerance, tool.selection_tolerance_unit = w(str(tolerance)), w("frames")
    tool.save_plots = lambda *args, **kwargs: None      # the figure (identical data) is not timed
    tool.display_plots = lambda: None

    async def choose(widget):
        tool.trajec, tool.output_dir = trajectory, out_dir
        with open(trajectory) as f:
            tool.num_atoms = int(f.readline().split()[0])
    tool.open_file_dialog = choose

    times = {"import": t_import}
    t = time.perf_counter()
    asyncio.run(tool.frames_counter(None))
    times["browse"] = time.perf_counter() - t
    for run in ("run 1", "run 2"):
        t = time.perf_counter()
        assert asyncio.run(tool.read_params(None))
        if variant == "legacy":
            tool.calc_mrt()
            tool.mrt_summary()
            asyncio.run(tool.write_mrt_total())
        else:
            tool.run_analysis()
        times[run] = time.perf_counter() - t
    times["peak_mb"] = peak_memory_mb()
    print(json.dumps(times))


def legacy_module_dir(revision):
    folder = Path(tempfile.gettempdir()) / f"mrt_classical_{revision}"
    folder.mkdir(exist_ok=True)
    source = subprocess.run(["git", "show", f"{revision}:meanResidenceTime_classical.py"],
                            cwd=SRC, capture_output=True, check=True).stdout
    (folder / "meanResidenceTime_classical.py").write_bytes(source)
    return str(folder)


def run_worker(variant, module_dir, trajectory, out_dir, args):
    command = [sys.executable, __file__, "--worker", variant, module_dir, trajectory, out_dir,
               str(args.centre), args.symbol, str(args.cutoff), str(args.tolerance)]
    result = subprocess.run(command, capture_output=True, text=True, cwd=SRC)
    if result.returncode:
        sys.exit(f"{variant} worker failed:\n{result.stderr}")
    return json.loads(result.stdout.strip().splitlines()[-1])


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        v, module_dir, traj, out, centre, symbol, cutoff, tol = sys.argv[2:]
        return worker(v, module_dir, traj, out, int(centre), symbol, float(cutoff), int(tol))

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--frames", type=int, nargs="+", default=[2000, 15121])
    parser.add_argument("--waters", type=int, default=143)
    parser.add_argument("--trajectory", help="time this file instead of synthetic ones")
    parser.add_argument("--centre", type=int, default=1)
    parser.add_argument("--symbol", default="O")
    parser.add_argument("--cutoff", type=float, default=3.2)
    parser.add_argument("--tolerance", type=int, default=5)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--legacy", default="c1e50e3")
    parser.add_argument("--no-legacy", action="store_true", help="time only the current code")
    args = parser.parse_args()

    if args.trajectory:
        trajectories = [args.trajectory]
    else:
        trajectories = []
        for n in args.frames:
            path = Path(tempfile.gettempdir()) / f"mrt_bench_{1 + 3 * args.waters}atoms_{n}frames.xyz"
            if not path.exists():
                print(f"generating {path} ...", flush=True)
                generate(path, args.waters, n)
            trajectories.append(str(path))

    variants = {"current": str(SRC)}
    if not args.no_legacy:
        variants = {"legacy": legacy_module_dir(args.legacy), **variants}
    keys = ("import", "browse", "run 1", "run 2", "peak_mb")
    for trajectory in trajectories:
        print(f"\n{trajectory}: {os.path.getsize(trajectory) / 1e6:.0f} MB, symbol {args.symbol}, "
              f"cutoff {args.cutoff}, tolerance {args.tolerance} frames, median of {args.repeat} processes")
        print(f"  {'':<8s}" + "".join(f"{k:>10s}" for k in keys))
        outputs = {}
        for variant, module_dir in variants.items():
            out_dir = tempfile.mkdtemp(prefix=f"mrt_{variant}_")
            runs = [run_worker(variant, module_dir, trajectory, out_dir, args) for _ in range(args.repeat)]
            row = [statistics.median(r[k] for r in runs) if runs[0][k] is not None else float("nan") for k in keys]
            print(f"  {variant:<8s}" + "".join(f"{v:>9.2f}{'M' if k == 'peak_mb' else 's'}" for k, v in zip(keys, row)))
            outputs[variant] = {name: (Path(out_dir) / name).read_bytes() for name in OUTPUTS}
        if "legacy" in outputs:
            same = all(outputs["legacy"][n] == outputs["current"][n] for n in OUTPUTS)
            print(f"  output files byte-identical to {args.legacy}: {same}")
            if not same:
                sys.exit(1)


if __name__ == "__main__":
    main()
