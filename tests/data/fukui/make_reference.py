"""Regenerate the Multiwfn reference data used by the Fukui tests.

NOT run by the test suite (it needs Multiwfn). Run it by hand, from the source
directory, only when a phenol_*.fchk fixture changes:

    python tests/data/fukui/make_reference.py path\\to\\Multiwfn.exe

It copies phenol_N/Np1/Nm1.fchk to N.fch/N+1.fch/N-1.fch in a temporary folder
and drives Multiwfn main function 22 (verified with Multiwfn 2026.7.15):
  option 2            -> CDFT.txt                  (saved as phenol_CDFT.txt)
  option 3, grid 1    -> low-quality grid; options 5-8 export f+.cub, f-.cub,
                         f0.cub, DD.cub
From the cubes it keeps 300 reproducible grid nodes and writes
multiwfn_reference/phenol_grid_sample.txt with columns
    x_bohr y_bohr z_bohr f+ f- f0 DD
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "multiwfn_reference"
SEED = 20261016
N_POINTS = 300
FILES = {"N.fch": "phenol_N.fchk", "N+1.fch": "phenol_Np1.fchk", "N-1.fch": "phenol_Nm1.fchk"}
PATHS = "N.fch\nN+1.fch\nN-1.fch\n"
SCRIPT = "22\n2\n" + PATHS + "3\n" + PATHS + "1\n5\n6\n7\n8\n0\n0\nq\n"  # N.fch is argv[1]


def read_cube(path: Path):
    lines = path.read_text().splitlines()
    head = lines[2].split()
    natoms, origin = int(head[0]), np.array([float(x) for x in head[1:4]])
    shape = tuple(int(lines[3 + k].split()[0]) for k in range(3))
    spacing = float(lines[3].split()[1])
    data = np.array([float(v) for l in lines[6 + natoms:] for v in l.split()]).reshape(shape)
    return origin, spacing, data


def main(multiwfn: str) -> None:
    settings = Path(multiwfn).resolve().parent / "settings.ini"
    lead = "" if settings.is_file() else "\n"  # without settings.ini Multiwfn pauses for ENTER
    with tempfile.TemporaryDirectory() as tmp:
        for target, source in FILES.items():
            shutil.copyfile(HERE / source, Path(tmp) / target)
        if settings.is_file():
            shutil.copyfile(settings, Path(tmp) / "settings.ini")
        subprocess.run([multiwfn, "N.fch"], input=lead + SCRIPT, text=True, cwd=tmp, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.copyfile(Path(tmp) / "CDFT.txt", HERE / "phenol_CDFT.txt")
        cubes = {k: read_cube(Path(tmp) / f"{k}.cub") for k in ("f+", "f-", "f0", "DD")}
    origin, spacing, first = cubes["f+"]
    rng = np.random.default_rng(SEED)
    idx = np.column_stack([rng.integers(n, size=N_POINTS) for n in first.shape])
    xyz = origin + spacing * idx
    OUT.mkdir(exist_ok=True)
    with (OUT / "phenol_grid_sample.txt").open("w") as handle:
        handle.write("# Multiwfn 22/3 low grid, options 5-8, sampled at 300 grid nodes\n")
        handle.write("# x_bohr y_bohr z_bohr f+ f- f0 DD   (a.u.)\n")
        for (i, j, k), (x, y, z) in zip(idx, xyz):
            values = [cubes[key][2][i, j, k] for key in ("f+", "f-", "f0", "DD")]
            handle.write(f"{x:.6f} {y:.6f} {z:.6f} " + " ".join(f"{v:.5e}" for v in values) + "\n")
    print("wrote phenol_CDFT.txt and multiwfn_reference/phenol_grid_sample.txt")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else shutil.which("Multiwfn.exe") or "Multiwfn.exe")
