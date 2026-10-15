"""Benchmark the single-solute solvent-box packing engines and the pre-optimizer.

    python benchmarks/bench_solvent_box.py [--sizes 15 20 30] [--skip-legacy-above 20]

For each cubic water box (target density 0.997 g/cm3, 100 000 attempts,
seed 7) every engine packs the box from the same seed; the script checks that
all engines produced the identical structure, then times a rigid and a
flexible pre-optimization of the packed box. Not run by the test suite.
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import solventPacking  # noqa: E402
from single_solute_solvent_box import SingleSoluteSolventBox  # noqa: E402

WATER = [["O", 0.0, 0.0, 0.1173], ["H", 0.0, 0.7572, -0.4692], ["H", 0.0, -0.7572, -0.4692]]


def builder(length, engine, attempts):
    b = SingleSoluteSolventBox()
    b.lattice_vectors = [length] * 3
    b.target_density = 0.997
    b.max_attempts = attempts
    b.vdw_scale = 0.8
    b.random_seed = 7
    b.spacing = 0.0
    b.min_distance = 0.0
    b.packing_engine = engine
    b.solvent_coords = [a[:] for a in WATER]
    b.center_of_mass()
    return b


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", type=float, nargs="+", default=[15.0, 20.0, 30.0])
    parser.add_argument("--attempts", type=int, default=100000)
    parser.add_argument("--skip-legacy-above", type=float, default=1e9)
    args = parser.parse_args()

    engines = ["legacy", "numpy"] + (["numba"] if solventPacking.numba_available() else [])
    if "numba" in engines:  # compile once, outside the timings
        builder(10.0, "numba", 200).build_box()

    print(f"{'box':>6} {'engine':>7} {'placed':>12} {'attempts':>9} {'time/s':>9} {'speedup':>8} identical")
    for size in args.sizes:
        reference = None
        legacy_time = None
        for engine in engines:
            if engine == "legacy" and size > args.skip_legacy_above:
                print(f"{size:6.1f} {engine:>7}      skipped")
                continue
            b = builder(size, engine, args.attempts)
            b.build_box()
            coords = np.array([a[1:] for a in b.data], dtype=float)
            if reference is None:
                reference = coords
                identical = "ref"
            else:
                identical = "yes" if np.array_equal(coords, reference) else "NO"
            if engine == "legacy":
                legacy_time = b.build_seconds
            speed = f"{legacy_time / b.build_seconds:7.1f}x" if legacy_time else "      -"
            print(f"{size:6.1f} {engine:>7} {b.solvent_insertion_count:5d}/{b.num_target_solvent_molecules:<6d}"
                  f" {b.solvent_attempts_used:9d} {b.build_seconds:9.2f} {speed} {identical}")
        for mode in ("rigid", "flexible"):
            b.preopt_mode = mode
            b.preopt_contact_scale = 0.9  # relax to a larger clearance than packing used (0.8)
            b.run_preoptimization()
            r = b.preopt_result
            print(f"{'':6} preopt {mode:>8}: {b.preopt_seconds:6.2f} s, {r.iterations} it, clashes "
                  f"{r.report_initial['clashes']} -> {r.report_final['clashes']}, "
                  f"max internal change {r.max_internal_deviation:.1e} A")



if __name__ == "__main__":
    main()
