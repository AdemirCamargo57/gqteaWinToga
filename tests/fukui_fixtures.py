"""Paths of the Fukui test fixtures (not a test module).

tests/data/fukui/phenol_{N,Np1,Nm1}.fchk: G09W B3LYP/6-31G(d) scf=tight single points
at Multiwfn's examples/phenol.xyz geometry, charge/multiplicity 0 1, -1 2, 1 2
(inputs/*.gjf). phenol_CDFT.txt: Multiwfn 2026.7.15 main function 22 option 2 on
those three files. multiwfn_reference/*.cub: Multiwfn option 3 (low grid), options
5-8, on the same files -- regenerate with tests/data/fukui/make_reference.py.
"""
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data" / "fukui"
N, PLUS, MINUS = (DATA / f"phenol_{s}.fchk" for s in ("N", "Np1", "Nm1"))
CDFT = DATA / "phenol_CDFT.txt"
REFERENCE_DIR = DATA / "multiwfn_reference"
