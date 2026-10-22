"""Regenerate the Multiwfn references for tests/test_ntoAnalysis.py.

NOT run by the test suite (it needs Multiwfn). Run it by hand from the source
directory with the path of Multiwfn.exe (verified with Multiwfn 2026.7.15):

    python tests/data/nto/make_reference.py path\\to\\Multiwfn.exe

For every (fixture, state) in nto_fixtures.REFERENCE_STATES, Multiwfn's NTO
function (main function 18, subfunction 6) reads the fixture's .fchk and
Gaussian output and exports the NTOs as an .fch (option 2). That file is stored
gzip-compressed as <fixture>_state<n>_multiwfn.fch.gz; its orbital energies are
the NTO eigenvalues (occupied ascending, virtual descending).
"""

import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))  # tests/, for nto_fixtures

from nto_fixtures import FIXTURES, REFERENCE_STATES  # noqa: E402

COMMANDS = "18\n6\n{log}\n{state}\n2\nref.fch\n0\nq\n"


def main(multiwfn: str) -> None:
    settings = Path(multiwfn).resolve().parent / "settings.ini"
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        lead = "\n"  # without settings.ini in the working folder Multiwfn pauses for ENTER first
        if settings.is_file():
            shutil.copy(settings, work / "settings.ini")
            lead = ""
        for name, states in REFERENCE_STATES.items():
            fchk, log = FIXTURES[name]
            shutil.copy(HERE / fchk, work / fchk)
            shutil.copy(HERE / log, work / log)
            for state in states:
                (work / "ref.fch").unlink(missing_ok=True)
                proc = subprocess.run([multiwfn, fchk], input=lead + COMMANDS.format(log=log, state=state),
                                      text=True, cwd=work, capture_output=True, timeout=600)
                if not (work / "ref.fch").is_file():
                    raise SystemExit(f"{name} state {state}: Multiwfn wrote no .fch\n{proc.stdout[-2000:]}")
                out = HERE / f"{name}_state{state}_multiwfn.fch.gz"
                with open(work / "ref.fch", "rb") as src, gzip.open(out, "wb", compresslevel=9) as dst:
                    shutil.copyfileobj(src, dst)
                lines = proc.stdout.splitlines()
                summary = [l.strip() for l in lines if "eigenvalues" in l.lower() or "Deviation" in l]
                print(f"{name} state {state}: {out.name}; " + " | ".join(summary))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
