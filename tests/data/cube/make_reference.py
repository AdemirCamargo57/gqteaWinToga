"""Regenerate the real cubegen cubes used by the Cube tab tests.

NOT run by the test suite (it needs Gaussian). Run it by hand, from the
source directory, only when these fixtures must change:

    python tests/data/cube/make_reference.py C:\\G09W\\cubegen.exe

From tests/data/fchk/eteno_b3lyp_631g.fchk (ethylene, B3LYP/6-31G, 16
electrons, HOMO = MO 8) it writes, on cubegen's coarse grid (-2, 3 points
per bohr) to keep the files small:
    ethylene_mo78.cube     orbital cube with MOs 7 and 8 together (negative atom count)
    ethylene_density.cube  total SCF density

Checked when generated (2026-10-19, G09W): MO 8 read from the two-MO file is
bit-identical to cubegen's MO=Homo cube; both orbitals integrate to 1 within
3e-4 and to each other to 4e-9. The density integrates to only ~12 of the 16
electrons: cubegen screens the density near nuclei and a 1/3-bohr grid
under-samples the 1s cores, so no test relies on that integral.
"""

import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FCHK = HERE.parent / "fchk" / "eteno_b3lyp_631g.fchk"
JOBS = {
    "ethylene_mo78.cube": "MO=7,8",
    "ethylene_density.cube": "Density=SCF",
}


def main(cubegen: str) -> None:
    env = dict(os.environ, GAUSS_EXEDIR=str(Path(cubegen).parent))
    for name, kind in JOBS.items():
        out = HERE / name
        subprocess.run([cubegen, "0", kind, str(FCHK), str(out), "-2", "h"], env=env, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"wrote {out.name} ({out.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
