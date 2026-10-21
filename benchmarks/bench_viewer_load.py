"""Time the 3D viewer's trajectory loading against its previous version and check identical frames.

    python benchmarks/bench_viewer_load.py                       (from venv/src)
    python benchmarks/bench_viewer_load.py --trajectory path/to/TRAJEC.xyz

Times ``MolecularViewer.load_xyz_trajectory`` (reading + frame-0 bonds) in a
fresh process per run, for the module as it is and as it was at ``--legacy``
(default c1e50e3, the line-by-line reader), read with ``git show``. The
synthetic trajectories are the ones benchmarks/bench_mrt_classical.py
generates (it is imported for that). The frames, offsets and element lists
of both versions are compared through a hash (PYTHONHASHSEED fixed).
"""
import argparse
import asyncio
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SRC = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_mrt_classical import generate, legacy_module_dir, peak_memory_mb   # noqa: E402


def worker(module_dir, trajectory):
    sys.path.insert(0, module_dir)
    sys.path.insert(1, str(SRC))
    import molecularViewer as mv
    viewer = mv.MolecularViewer()
    t = time.perf_counter()
    asyncio.run(viewer.load_xyz_trajectory(trajectory))
    seconds = time.perf_counter() - t
    digest = hash((tuple(map(tuple, viewer.frames)), tuple(viewer.frame_offsets), tuple(viewer.frame_elements)))
    print(json.dumps({"load": seconds, "peak_mb": peak_memory_mb(), "frames": len(viewer.frames),
                      "digest": digest}))


def legacy_viewer_dir(revision):
    folder = Path(legacy_module_dir(revision))     # also holds that revision's meanResidenceTime_classical
    source = subprocess.run(["git", "show", f"{revision}:molecularViewer.py"],
                            cwd=SRC, capture_output=True, check=True).stdout
    (folder / "molecularViewer.py").write_bytes(source)
    return str(folder)


def run(module_dir, trajectory):
    env = dict(os.environ, PYTHONHASHSEED="0")
    result = subprocess.run([sys.executable, __file__, "--worker", module_dir, trajectory],
                            capture_output=True, text=True, cwd=SRC, env=env)
    if result.returncode:
        sys.exit(result.stderr)
    return json.loads(result.stdout.strip().splitlines()[-1])


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        return worker(sys.argv[2], sys.argv[3])
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--frames", type=int, nargs="+", default=[2000, 15121])
    parser.add_argument("--waters", type=int, default=143)
    parser.add_argument("--trajectory")
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--legacy", default="c1e50e3")
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

    variants = {"legacy": legacy_viewer_dir(args.legacy), "current": str(SRC)}
    for trajectory in trajectories:
        print(f"\n{trajectory}: {os.path.getsize(trajectory) / 1e6:.0f} MB, median of {args.repeat} processes")
        digests = {}
        for name, module_dir in variants.items():
            runs = [run(module_dir, trajectory) for _ in range(args.repeat)]
            load = statistics.median(r["load"] for r in runs)
            peak = statistics.median(r["peak_mb"] for r in runs) if runs[0]["peak_mb"] else float("nan")
            print(f"  {name:<8s} load {load:7.2f} s   peak {peak:7.0f} MB   {runs[0]['frames']} frames")
            digests[name] = {r["digest"] for r in runs}
        same = digests["legacy"] == digests["current"] and len(digests["current"]) == 1
        print(f"  frames, offsets and elements identical to {args.legacy}: {same}")
        if not same:
            sys.exit(1)


if __name__ == "__main__":
    main()
