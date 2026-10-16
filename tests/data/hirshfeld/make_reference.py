"""Regenerate multiwfn_hirshfeld.txt, the Multiwfn reference for tests/test_hirshfeldCharges.py.

NOT run by the test suite (it needs Multiwfn). Run it by hand from the source
directory only when a fixture here changes:

    python tests/data/hirshfeld/make_reference.py path\\to\\Multiwfn.exe

The fixtures are G09W B3LYP/6-31G(d) single points at rough, unoptimized
geometries (inputs/*.gjf): only the density matters for comparing charges.
For each <name>.fchk Multiwfn main function 7, option 1 (Hirshfeld), option 1
(built-in free-atom densities) prints "Final atomic charges"
(verified with Multiwfn 2026.7.15); they are written as "<name> <charge>",
one line per atom in file order.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIXTURES = ["nh3bf3", "ch3scl", "h3po", "lih", "nah", "mgh2", "alh3", "sih4"]


def final_charges(stdout: str):
    lines = stdout.splitlines()
    start = next(i for i, l in enumerate(lines) if "Final atomic charges" in l)
    out = []
    for line in lines[start + 1:]:
        if not line.strip().startswith("Atom"):
            break
        out.append(float(line.split()[-1]))
    return out


def main(multiwfn: str) -> None:
    settings = Path(multiwfn).resolve().parent / "settings.ini"
    lead = "" if settings.is_file() else "\n"  # without settings.ini Multiwfn pauses for ENTER
    rows = ["# Multiwfn 2026.7.15 Hirshfeld charges (main function 7, option 1, built-in densities)",
            "# fixture charge   (one line per atom, file order)"]
    with tempfile.TemporaryDirectory() as tmp:
        if settings.is_file():
            shutil.copyfile(settings, Path(tmp) / "settings.ini")
        for name in FIXTURES:
            shutil.copyfile(HERE / f"{name}.fchk", Path(tmp) / "in.fch")
            proc = subprocess.run([multiwfn, "in.fch"], input=lead + "7\n1\n1\nn\n0\nq\n", text=True,
                                  cwd=tmp, capture_output=True, check=True)
            rows += [f"{name} {q:.8f}" for q in final_charges(proc.stdout)]
    (HERE / "multiwfn_hirshfeld.txt").write_text("\n".join(rows) + "\n")
    print(f"wrote multiwfn_hirshfeld.txt ({len(rows) - 2} atoms)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else shutil.which("Multiwfn.exe") or "Multiwfn.exe")
