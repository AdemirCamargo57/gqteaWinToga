"""Regenerate the Gaussian 16 outputs of tests/data/gaussian_output (never run by the test suite).

Needs Gaussian 16 (Rev. C.01) in WSL with g16 and formchk on the login shell's PATH. Run it from the
source directory on Windows:

    python tests/data/gaussian_output/run_gaussian.py

Each inputs/<name>.gjf is copied to /tmp/gqout in WSL (WSL's own filesystem: /mnt/c is very slow),
run as  g16 < <name>.gjf > <name>.log, and the logs are copied back. h2o_opt_freq is stored as .out
(both suffixes are tested); h2co_opt_td_excited.chk is formchk'ed and stored gzipped. The
h2co_opt_td_gen job only proves that ChkBasis runs: its output is not kept.
"""

import gzip
import shutil
import subprocess
import sys
from pathlib import Path, PureWindowsPath

HERE = Path(__file__).resolve().parent
WORK = "/tmp/gqout"
NAMES = ("h2o_opt", "h2o_opt_freq", "h2co_opt_maxcycles", "h2o_opt_then_failed", "h2co_opt_td", "h2co_opt_td_gen")
FAILING = ("h2co_opt_maxcycles", "h2o_opt_then_failed")     # they end in an error termination on purpose
STEPS = {"h2o_opt_freq": 2, "h2o_opt_then_failed": 1, "h2co_opt_td": 2, "h2co_opt_td_gen": 2}


def wsl(path: Path) -> str:
    p = PureWindowsPath(str(path))
    return "/mnt/" + p.drive[0].lower() + "/" + "/".join(p.parts[1:])


def bash(command: str) -> str:
    return subprocess.run(["wsl", "-e", "bash", "-lc", command], check=True, capture_output=True, text=True,
                          timeout=3600).stdout


def main() -> None:
    bash(f"rm -rf {WORK} && mkdir -p {WORK} && cp {wsl(HERE / 'inputs')}/*.gjf {WORK}/")
    for name in NAMES:
        bash(f"cd {WORK} && g16 < {name}.gjf > {name}.log || true")
        count = int(bash(f"grep -c 'Normal termination' {WORK}/{name}.log || true").strip() or 0)
        if count != STEPS.get(name, 1 if name not in FAILING else 0):
            raise SystemExit(f"{name}: {count} Normal termination(s), not the expected number.")
        print(f"{name}: {count} Normal termination(s)")
    bash(f"cd {WORK} && formchk h2co_opt_td_excited.chk h2co_opt_td_excited.fchk")
    bash(f"cp {WORK}/*.log {WORK}/h2co_opt_td_excited.fchk {wsl(HERE)}/")
    (HERE / "h2o_opt_freq.log").replace(HERE / "h2o_opt_freq.out")
    (HERE / "h2co_opt_td_gen.log").unlink()
    fchk = HERE / "h2co_opt_td_excited.fchk"
    with open(fchk, "rb") as src, gzip.open(HERE / "h2co_opt_td_excited.fchk.gz", "wb", compresslevel=9) as dst:
        shutil.copyfileobj(src, dst)
    fchk.unlink()


if __name__ == "__main__":
    if sys.platform != "win32":
        raise SystemExit("Run it on Windows; it drives Gaussian 16 through WSL.")
    main()
