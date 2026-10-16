"""Regenerate proatomDensities.py: spherically averaged free-atom densities, H-Ar.

NOT run by the test suite or the application (it needs Gaussian). Run it by
hand from the source directory only to change the level or add elements:

    python tools/make_proatom_densities.py [C:\\G09W] [work-folder]

Level (the one Multiwfn's built-in Hirshfeld densities use for Z <= 18,
Multiwfn manual Appendix 3 / section 6.3): UB3LYP/cc-pVQZ, ground-state
multiplicity, neutral atom. Each atom's SCF density is evaluated with
gaussianWavefunction on a 110-point Lebedev shell (exact for the angular
dependence of an s-g density) at each of 600 radii log-spaced from 1e-5 to
40 bohr, and averaged. Existing <El>.fchk files in the work folder are reused.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import lebedev_rule

SRC = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SRC))
import gaussianWavefunction as gw  # noqa: E402

LEVEL = "UB3LYP/cc-pVQZ"
GROUND_STATE_MULTIPLICITY = {
    "H": 2, "He": 1, "Li": 2, "Be": 1, "B": 2, "C": 3, "N": 4, "O": 3, "F": 2, "Ne": 1,
    "Na": 2, "Mg": 1, "Al": 2, "Si": 3, "P": 4, "S": 3, "Cl": 2, "Ar": 1,
}
RADII_BOHR = np.geomspace(1e-5, 40.0, 600)
LEBEDEV_ORDER = 17  # 110 points


def run_atom(symbol: str, gaussian_dir: Path, work: Path) -> Path:
    fchk = work / f"{symbol}.fchk"
    if fchk.is_file():
        return fchk
    gjf = work / f"{symbol}.gjf"
    gjf.write_bytes(  # LF and float coordinates: G09W rejects CRLF here and integer coordinates
        (f"%chk={symbol}.chk\n%mem=1GB\n#p {LEVEL} scf=(tight,maxcycle=300)\n\n"
         f"{symbol} free atom\n\n0 {GROUND_STATE_MULTIPLICITY[symbol]}\n{symbol} 0.0 0.0 0.0\n\n").encode())
    env = dict(os.environ, GAUSS_EXEDIR=str(gaussian_dir), GAUSS_SCRDIR=str(work))
    subprocess.run([str(gaussian_dir / "g09.exe"), gjf.name, f"{symbol}.out"], cwd=work, env=env, check=True)
    if "Normal termination" not in (work / f"{symbol}.out").read_text(errors="replace"):
        raise RuntimeError(f"Gaussian failed for {symbol}; see {work / (symbol + '.out')}")
    subprocess.run([str(gaussian_dir / "formchk.exe"), f"{symbol}.chk", fchk.name], cwd=work, env=env, check=True)
    return fchk


def spherical_average(fchk: Path) -> np.ndarray:
    w = gw.read_fchk(fchk)
    xyz, weights = lebedev_rule(LEBEDEV_ORDER)
    weights = weights / weights.sum()
    points = (RADII_BOHR[:, None, None] * xyz.T[None, :, :]).reshape(-1, 3)
    # numpy path: unscreened, so the tail stays positive out to 40 bohr
    rho = gw.electron_density(w, points, use_numba=False).reshape(len(RADII_BOHR), -1)
    return rho @ weights


def main(gaussian_dir: str, work_dir: str) -> None:
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    tables = {}
    for symbol in GROUND_STATE_MULTIPLICITY:
        tables[symbol] = spherical_average(run_atom(symbol, Path(gaussian_dir), work))
        print(f"{symbol}: {np.trapezoid(4 * np.pi * RADII_BOHR ** 2 * tables[symbol], RADII_BOHR):.5f} electrons")

    def row(values):
        return "(" + ", ".join(f"{v:.9e}" for v in values) + ")"

    lines = [
        '"""Spherically averaged free-atom electron densities for Hirshfeld charges (GENERATED).',
        "",
        "Do not edit: regenerate with tools/make_proatom_densities.py (needs Gaussian).",
        f"Level: {LEVEL}, neutral ground-state atoms (Multiwfn's level for Z <= 18), each",
        "density averaged over a 110-point Lebedev shell at every radius. Units: bohr and",
        "electrons/bohr^3. Interpolate log(rho) linearly in r (hirshfeldCharges.proatom_density).",
        '"""',
        "",
        f'LEVEL = "{LEVEL}"',
        f"RADII_BOHR = {row(RADII_BOHR)}",
        "DENSITIES = {",
    ]
    lines += [f'    "{s}": {row(v)},' for s, v in tables.items()]
    lines.append("}")
    (SRC / "proatomDensities.py").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote proatomDensities.py")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else r"C:\G09W",
         sys.argv[2] if len(sys.argv) > 2 else str(Path.cwd() / "proatom_work"))
