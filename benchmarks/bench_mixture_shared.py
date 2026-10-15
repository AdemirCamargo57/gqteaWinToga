"""Benchmark the mixture and shared-wall packing engines (legacy vs numpy vs numba).

    python benchmarks/bench_mixture_shared.py

Water/methanol 60:40, 50 000 attempts per box, fixed seeds; checks that every
engine produced the identical structure. Not run by the test suite.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import solventPacking  # noqa: E402
from mixtureSolventBox import MixtureSolvationBox  # noqa: E402
from sharedWallBox import SharedWallBox  # noqa: E402

WATER = [["O", 0.0, 0.0, 0.1173], ["H", 0.0, 0.7572, -0.4692], ["H", 0.0, -0.7572, -0.4692]]
MEOH = [["C", -0.0465, 0.6630, 0.0], ["O", -0.0465, -0.7570, 0.0], ["H", -1.0861, 0.9729, 0.0],
        ["H", 0.4366, 1.0816, 0.8904], ["H", 0.4366, 1.0816, -0.8904], ["H", 0.8624, -1.0857, 0.0]]


def mixture(length, engine):
    b = MixtureSolvationBox()
    b.lattice_vectors = [length] * 3
    b.target_density, b.max_attempts, b.vdw_scale, b.random_seed = 0.9, 50000, 0.8, 7
    b.composition_A, b.composition_B = 0.6, 0.4
    b.coords_A, b.coords_B = [a[:] for a in WATER], [a[:] for a in MEOH]
    b.packing_engine = engine
    b.prepare_molecules()
    b.estimate_target_molecule_count()
    b.build_box()
    return b, b.data, f"{b.count_solvent_A + b.count_solvent_B}/{b.num_target_solvent_molecules}"


def shared(length, engine):
    b = SharedWallBox()
    b.shared_z = length
    b.lattice_vectors = [length, length, 2 * length]
    for box, comp, seed in ((b.boxes[0], 1.0, 3), (b.boxes[1], 0.6, 4)):
        box.update(target_density=0.95, max_attempts=50000, vdw_scale=0.8, random_seed=seed,
                   composition_A=comp, composition_B=1.0 - comp,
                   coords_A=[a[:] for a in WATER], coords_B=[a[:] for a in MEOH])
    b.packing_engine = engine
    b.build_box()
    placed = sum(x["count_solvent_A"] + x["count_solvent_B"] for x in b.boxes)
    target = sum(x["num_target_solvent_molecules"] for x in b.boxes)
    return b, b.combined_data, f"{placed}/{target}"


def main():
    engines = ["legacy", "numpy"] + (["numba"] if solventPacking.numba_available() else [])
    if "numba" in engines:
        mixture(10.0, "numba")  # compile outside the timings
    print(f"{'builder':>8} {'box':>5} {'engine':>7} {'placed':>11} {'time/s':>8} {'speedup':>8} identical")
    for name, fn, sizes in (("mixture", mixture, (20.0, 30.0)), ("shared", shared, (20.0, 26.0))):
        for size in sizes:
            reference = legacy_time = None
            for engine in engines:
                started = time.perf_counter()
                _, data, placed = fn(size, engine)
                elapsed = time.perf_counter() - started
                coords = np.array([a[1:] for a in data], dtype=float)
                if reference is None:
                    reference, legacy_time, identical = coords, elapsed, "ref"
                else:
                    identical = "yes" if np.array_equal(coords, reference) else "NO"
                print(f"{name:>8} {size:5.1f} {engine:>7} {placed:>11} {elapsed:8.2f} "
                      f"{legacy_time / elapsed:7.1f}x {identical}", flush=True)


if __name__ == "__main__":
    main()
