"""Regenerate the cubegen MO cubes of tests/data/orbitals (<name>.cube.gz).

NOT run by the test suite (it needs Gaussian's cubegen). Run it by hand from
the source directory (verified with G09W, C:\\G09W):

    python tests/data/orbitals/make_reference.py [path\\to\\cubegen.exe]

Each entry of REFERENCES is one orbital on cubegen's coarse grid (-2: 3
points per bohr), written as a plain cube ("h" format) and gzip-compressed.
tests/test_molecularOrbitals.py evaluates the same orbital at the cube's grid
points and compares.
"""

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent

# (output name, fchk relative to tests/data, cubegen kind)
REFERENCES = (
    ("water_5d_mo5", "fchk/water_631gd_5d.fchk", "MO=5"),                 # restricted HOMO
    ("water_5d_mo6", "fchk/water_631gd_5d.fchk", "MO=6"),                 # restricted LUMO
    ("ch3_u_alpha5", "fchk/ch3_ub3lyp_631gd.fchk", "AMO=5"),              # unrestricted alpha HOMO
    ("ch3_u_beta5", "fchk/ch3_ub3lyp_631gd.fchk", "BMO=5"),               # unrestricted beta LUMO
    ("ch3_ro_mo5", "orbitals/ch3_rob3lyp_631gd.fchk", "MO=5"),            # restricted-open SOMO
)


def main(cubegen: str) -> None:
    cubegen = Path(cubegen)
    env = dict(os.environ, GAUSS_EXEDIR=str(cubegen.parent))
    for name, fchk, kind in REFERENCES:
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            shutil.copy(DATA / fchk, work / "in.fchk")
            subprocess.run([str(cubegen), "0", kind, "in.fchk", "out.cube", "-2", "h"], cwd=work, env=env,
                           check=True, capture_output=True, timeout=600)
            with open(work / "out.cube", "rb") as src, gzip.open(HERE / f"{name}.cube.gz", "wb") as dst:
                shutil.copyfileobj(src, dst)
            print(f"{name}: {kind} of {fchk}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else r"C:\G09W\cubegen.exe")
