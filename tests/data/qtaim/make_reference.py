"""Regenerate the Multiwfn references for tests/test_qtaimAnalysis.py.

NOT run by the test suite (it needs Multiwfn). Run it by hand from the source
directory with the path of Multiwfn.exe (verified with Multiwfn 2026.7.15):

    python tests/data/qtaim/make_reference.py path\\to\\Multiwfn.exe

For each fixture, Multiwfn's topology module (main function 2) searches CPs
from the nuclei, atom-pair midpoints, triangle and tetrahedron centres
(modes 2, 3, 4, 5, default parameters), generates the bond paths (8) and the
ring paths (9), prints the path summary (-5, 1) and writes every CP's
properties to CPprop.txt (7, 0). The parsed values go to <name>_multiwfn.txt:

    cp <index> <type> <x> <y> <z>                      (bohr)
    prop <rho> <lap> <G> <K> <V> <H> <l1> <l2> <l3> <ellipticity> <eta>
    atoms <i> <j>                                      (BCP "Connected atoms", 1-based)
    path <from> <to> <length>                          (bohr, one line per half path)
"""

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
FCHK = HERE.parent / "fchk"
FIXTURES = ["water_631gd_6d.fchk", "eteno_b3lyp_631g.fchk", "benzene_631gd.fchk", "ch3_ub3lyp_631gd.fchk",
            "water_dimer_631gd.fchk", "bicyclopentane_631gd.fchk"]
# Option 9 is only sent when a first pass found both (3,+1) and (3,+3): without
# them Multiwfn prints an error and waits for ENTER, swallowing the next command.
COMMANDS = "2\n2\n3\n4\n5\n8\n{ring}-5\n1\n0\n7\n0\n-10\nq\n"
NUMBER = r"[-+]?\d+\.\d+(?:[EeDd][-+]?\d+)?"

PROPERTY_LABELS = [
    ("rho", "Density of all electrons:"),
    ("lap", "Laplacian of electron density:"),
    ("G", "Lagrangian kinetic energy G(r):"),
    ("K", "Hamiltonian kinetic energy K(r):"),
    ("V", "Potential energy density V(r):"),
    ("H", "Energy density E(r) or H(r):"),
]


def fixture_path(name: str) -> Path:
    return HERE / name if (HERE / name).is_file() else FCHK / name


def parse_cpprop(text: str):
    blocks = re.split(r"-{5,}\s+CP\s+(\d+),\s+Type\s+(\(3,[-+]\d\))\s+-{5,}", text)
    cps = []
    for k in range(1, len(blocks), 3):
        index, kind, body = int(blocks[k]), blocks[k + 1], blocks[k + 2]
        pos = [float(x) for x in re.search(r"Position \(Bohr\):\s+(" + NUMBER + r")\s+(" + NUMBER + r")\s+("
                                           + NUMBER + ")", body).groups()]
        props = []
        for _, label in PROPERTY_LABELS:
            props.append(float(re.search(re.escape(label) + r"\s+(" + NUMBER + ")", body).group(1)))
        eig = re.search(r"Eigenvalues of Hessian:\s+(" + NUMBER + r")\s+(" + NUMBER + r")\s+(" + NUMBER + ")", body)
        props += [float(x) for x in eig.groups()]
        for label in ("Ellipticity of electron density:", "eta index:"):
            m = re.search(re.escape(label) + r"\s+(" + NUMBER + ")", body)
            props.append(float(m.group(1)) if m else float("nan"))
        atoms = re.search(r"Connected atoms:\s+(\d+)\(\w+\s*\)\s+--\s+(\d+)\(", body)
        cps.append((index, kind, pos, props, tuple(int(a) for a in atoms.groups()) if atoms else None))
    return cps


def parse_paths(stdout: str):
    pattern = r"#\s+\d+\s+CP:\s+(\d+)\s+\(3,[-+]\d\)\s+--->\s+CP:\s+(\d+)\s+\(3,[-+]\d\)\s+Length:\s+(" + NUMBER + ")"
    return [(int(a), int(b), float(c)) for a, b, c in re.findall(pattern, stdout)]


def main(multiwfn: str) -> None:
    settings = Path(multiwfn).resolve().parent / "settings.ini"
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        # without settings.ini in the working folder Multiwfn pauses for ENTER first
        lead = "" if (work / "settings.ini").is_file() else "\n"
        if settings.is_file():
            shutil.copy(settings, work / "settings.ini")
            lead = ""
        for name in FIXTURES:
            shutil.copy(fixture_path(name), work / "in.fch")
            ring = ""
            for _ in range(2):
                (work / "CPprop.txt").unlink(missing_ok=True)
                proc = subprocess.run([multiwfn, "in.fch"], input=lead + COMMANDS.format(ring=ring), text=True,
                                      cwd=work, capture_output=True, timeout=600)
                cps = parse_cpprop((work / "CPprop.txt").read_text(errors="replace"))
                kinds = {kind for _, kind, _, _, _ in cps}
                if ring or not {"(3,+1)", "(3,+3)"} <= kinds:
                    break
                ring = "9\n"
            paths = parse_paths(proc.stdout)
            rows = [f"# Multiwfn 2026.7.15 topology of {name}: main function 2, modes 2 3 4 5, paths 8 9, "
                    "default parameters"]
            for index, kind, pos, props, atoms in cps:
                rows.append(f"cp {index} {kind} " + " ".join(f"{x:.12f}" for x in pos))
                rows.append("prop " + " ".join(f"{x:.10E}" for x in props))
                if atoms:
                    rows.append(f"atoms {atoms[0]} {atoms[1]}")
            for a, b, length in paths:
                rows.append(f"path {a} {b} {length:.5f}")
            out = HERE / (Path(name).stem + "_multiwfn.txt")
            out.write_text("\n".join(rows) + "\n")
            print(f"{out.name}: {len(cps)} CPs, {len(paths)} paths")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else shutil.which("Multiwfn.exe") or "Multiwfn.exe")
