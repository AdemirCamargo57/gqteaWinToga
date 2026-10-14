"""Regenerate the cubegen reference values used by the ESP tests.

NOT run by the test suite (it needs Gaussian). Run it by hand, from the
source directory, only when a fixture in tests/data/fchk changes:

    python tests/data/esp_reference/make_reference.py C:\\G09W\\cubegen.exe

For every tests/data/fchk/<name>.fchk it draws 300 reproducible points,
0.5-8 bohr from a random atom, asks cubegen for the SCF density and the ESP
there, and writes tests/data/esp_reference/<name>.txt with columns
    x_A  y_A  z_A  rho_au  esp_au
cubegen reads free-format point coordinates from stdin **in Angstrom**
(verified 2026-10-14 with G09W: a point typed as 0.1173 lands on the water
oxygen at 0.1173 A, not 0.1173 bohr), and echoes them back unchanged.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

BOHR_RADIUS_ANGSTROM = 0.529177210544
HERE = Path(__file__).resolve().parent
FCHK_DIR = HERE.parent / "fchk"
N_POINTS = 300
SEED = 20261014


def fchk_coordinates_bohr(path: Path) -> np.ndarray:
    lines = path.read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("Current cartesian coordinates"))
    count = int(lines[start].split()[-1])
    values = []
    k = start + 1
    while len(values) < count:
        values += [float(x) for x in lines[k].split()]
        k += 1
    return np.array(values).reshape(-1, 3)


def sample_points_angstrom(coords_bohr: np.ndarray) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    atoms = coords_bohr[rng.integers(len(coords_bohr), size=N_POINTS)]
    direction = rng.normal(size=(N_POINTS, 3))
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    distance = rng.uniform(0.5, 8.0, size=N_POINTS)
    return (atoms + direction * distance[:, None]) * BOHR_RADIUS_ANGSTROM


def run_cubegen(cubegen: str, fchk: Path, kind: str, points: np.ndarray) -> np.ndarray:
    env = dict(os.environ, GAUSS_EXEDIR=str(Path(cubegen).parent))
    with tempfile.TemporaryDirectory() as tmp:
        cube = Path(tmp) / "out.cub"
        stdin = "".join(f"{x:.10f} {y:.10f} {z:.10f}\n" for x, y, z in points)
        subprocess.run([cubegen, "0", f"{kind}=SCF", str(fchk), str(cube), "-5", "h"],
                       input=stdin, text=True, env=env, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        lines = cube.read_text().splitlines()
    natoms = abs(int(lines[2].split()[0]))
    data = np.array([[float(v) for v in l.split()] for l in lines[6 + natoms:] if l.strip()])
    if data.shape[0] != len(points) or not np.allclose(data[:, :3], points, atol=1e-8):
        raise RuntimeError(f"cubegen did not echo the requested points for {fchk.name}")
    return data[:, 3]


def main(cubegen: str) -> None:
    for fchk in sorted(FCHK_DIR.glob("*.fchk")):
        points = sample_points_angstrom(fchk_coordinates_bohr(fchk))
        rho = run_cubegen(cubegen, fchk, "density", points)
        esp = run_cubegen(cubegen, fchk, "potential", points)
        out = HERE / f"{fchk.stem}.txt"
        with out.open("w") as handle:
            handle.write(f"# cubegen reference for {fchk.name}: density=SCF and potential=SCF\n")
            handle.write("# x_A y_A z_A rho_au esp_au   (points sent to cubegen in Angstrom)\n")
            for (x, y, z), r, v in zip(points, rho, esp):
                handle.write(f"{x:.10f} {y:.10f} {z:.10f} {r:.12f} {v:.12f}\n")
        print(f"wrote {out.name}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else r"C:\G09W\cubegen.exe")
