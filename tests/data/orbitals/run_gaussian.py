"""Regenerate the Gaussian fixtures of tests/data/orbitals (the .log and .fchk files).

NOT run by the test suite (it needs Gaussian). Run it by hand from the source
directory (verified with G09W, C:\\G09W):

    python tests/data/orbitals/run_gaussian.py [path\\to\\g09.exe]

Each inputs/<name>.gjf is run in a temporary folder; <name>.log and the
formchk'ed <name>.fchk are copied here.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAMES = ("ch3_rob3lyp_631gd",)


def main(g09: str) -> None:
    g09 = Path(g09)
    formchk = g09.with_name("formchk.exe")
    env = dict(os.environ, GAUSS_EXEDIR=str(g09.parent))
    for name in NAMES:
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            shutil.copy(HERE / "inputs" / f"{name}.gjf", work / f"{name}.gjf")
            env["GAUSS_SCRDIR"] = str(work)
            subprocess.run([str(g09), f"{name}.gjf", f"{name}.log"], cwd=work, env=env, check=True, timeout=3600)
            log = (work / f"{name}.log").read_text(errors="replace")
            if "Normal termination" not in log:
                raise SystemExit(f"{name}: Gaussian did not terminate normally.")
            subprocess.run([str(formchk), f"{name}.chk", f"{name}.fchk"], cwd=work, env=env, check=True,
                           capture_output=True, timeout=600)
            shutil.copy(work / f"{name}.log", HERE / f"{name}.log")
            shutil.copy(work / f"{name}.fchk", HERE / f"{name}.fchk")
            print(f"{name}: done")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else r"C:\G09W\g09.exe")
