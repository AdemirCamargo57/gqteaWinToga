"""Time runatom.exe on every library case and measure event-loop lag while it runs.

Usage (from venv/src):  python benchmarks/bench_runatom.py [--lag]
"""

import asyncio
import shutil
import statistics
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import runatomInput as ri        # noqa: E402
import runatomRunner as rr       # noqa: E402


def library_cases():
    for makefile in sorted(rr.library_dir().rglob("Makefile")):
        folder = makefile.parent
        if folder.parent == rr.library_dir() or (folder / "V_ps.adat").exists():
            continue
        yield folder


def time_library(exe: str):
    times, failed = [], []
    for folder in library_cases():
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copytree(folder, tmp, dirs_exist_ok=True)
            case = ri.read_case_folder(tmp)
            started = time.perf_counter()
            try:
                results = rr.run_job(rr.runatom_command(exe), rr.plan_stages(case), tmp)
            except rr.RunatomError as problem:      # a defective upstream input, e.g. 034-Se
                print(f"{folder.name:34s} failed: {problem}")
                failed.append(folder.name)
                continue
            total = time.perf_counter() - started
        times.append(total)
        stages = ", ".join(f"{r.stage.kind} {r.seconds:.2f}" for r in results)
        print(f"{folder.name:34s} {total:6.2f} s  ({stages})")
    print(f"\n{len(times)} cases: median {statistics.median(times):.2f} s, max {max(times):.2f} s")
    if failed:
        print(f"{len(failed)} case(s) rejected by runatom itself (the Linux build fails the same way): {', '.join(failed)}")


async def measure_lag(exe: str, folder: Path):
    loop = asyncio.get_running_loop()
    lags, done = [], asyncio.Event()

    async def ticker():
        while not done.is_set():
            before = loop.time()
            await asyncio.sleep(0.01)
            lags.append(loop.time() - before - 0.01)

    with tempfile.TemporaryDirectory() as tmp:
        shutil.copytree(folder, tmp, dirs_exist_ok=True)
        case = ri.read_case_folder(tmp)
        tick = asyncio.create_task(ticker())
        progress = lambda info: loop.call_soon_threadsafe(lambda: None)
        await asyncio.to_thread(rr.run_job, rr.runatom_command(exe), rr.plan_stages(case), tmp, progress)
        done.set()
        await tick
    lags.sort()
    print(f"event-loop lag during a run: p95 {1000 * lags[int(0.95 * len(lags))]:.1f} ms, max {1000 * lags[-1]:.1f} ms")


if __name__ == "__main__":
    exe = rr.find_runatom()
    if exe is None:
        sys.exit("runatom.exe not found")
    time_library(exe)
    if "--lag" in sys.argv:
        asyncio.run(measure_lag(exe, rr.library_dir() / "056-Ba" / "056-Ba-gpw-nsp-campos"))
