"""Suggested runatom inputs for an element (Tools > Vanderbilt USPP Generator).

No Toga here (a subprocess test pins it). Three sources, best first:

1. a tested case of the bundled Work/ library, read with runatomInput's own
   reader and suggested unmodified (PBE preferred);
2. a same-group template: a computed all-electron input plus the
   pseudopotential values of the nearest library element in the same group,
   copied unchanged, with a warning;
3. a computed starter: the all-electron input only, the radii left blank.

Every suggestion is a dict of the window's form fields, so it reaches a file
through case_from_fields and format_adat like anything the user types.
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import runatomInput as ri
from help import AtomicData

MAX_DEFAULTS_Z = 92
PBE = 5.0
FUNCTIONAL_PREFERENCE = (5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0)   # PBE, PW91, other GGA, LDA
ELEMENT_FOLDER = re.compile(r"^(\d{3})-([A-Z][a-z]?)$")
L_LETTERS = "spdf"


@dataclass(frozen=True)
class LibraryVariant:
    folder: Path
    z: int
    file_symbol: str             # the Makefile's ATOM ("ti", "Ti"): kept, so file names match the library
    exfact: float
    irel: int
    xion: float
    nlcc: bool                   # ps.ifpcor == 1
    core: Tuple[int, ...]        # nnlz of the ncores core shells
    valence: Tuple[int, ...]     # nnlz of the nvales valence shells, in file order
    n_tests: int
    label: str


def shell_name(nnlz: int) -> str:
    """320 -> '3d'."""
    return f"{nnlz // 100}{L_LETTERS[(nnlz // 10) % 10]}"


def library_elements(roots: Sequence[Path]) -> Dict[str, List[Path]]:
    """'Ti' -> its case folders (those holding a Makefile), from the folder names alone."""
    found: Dict[str, List[Path]] = {}
    for root in roots:
        root = Path(root)
        if not root.is_dir():
            continue
        for element_dir in sorted(root.iterdir()):
            match = ELEMENT_FOLDER.match(element_dir.name)
            if not match or int(match.group(1)) != AtomicData.atomic_numbers.get(match.group(2)):
                continue
            cases = sorted(p.parent for p in element_dir.glob("*/Makefile"))
            if cases:
                found.setdefault(match.group(2), []).extend(cases)
    return found


def _semicore(valence: Sequence[int]) -> List[int]:
    """Valence s and p shells below the outermost n."""
    top = max(nnlz // 100 for nnlz in valence)
    return [nnlz for nnlz in valence if nnlz // 100 < top and (nnlz // 10) % 10 <= 1]


def _describe(folder: Path, case: ri.CaseInputs) -> LibraryVariant:
    ae, ps = case.ae, case.ps
    shells = [s.nnlz for s in ae.shells]
    core, valence = tuple(shells[:ps.ncores]), tuple(shells[ps.ncores:])
    parts = [ri.xc_label(ae.exfact)]
    if ae.irel == 2:
        parts.append("scalar-relativistic")
    semicore = _semicore(valence)
    if semicore:
        parts.append("semicore " + " ".join(shell_name(s) for s in semicore))
    if ps.ifpcor == 1:
        parts.append("NLCC")
    if ae.xion:
        parts.append(f"{ae.xion:+g} ion")
    parts.append(f"{len(case.tests)} test" + ("" if len(case.tests) == 1 else "s"))
    return LibraryVariant(folder, round(ae.z), case.symbol, ae.exfact + 0.0, ae.irel, ae.xion,
                          ps.ifpcor == 1, core, valence, len(case.tests),
                          f"{folder.name}: " + ", ".join(parts))


def element_variants(symbol: str, roots: Sequence[Path]) -> Tuple[List[LibraryVariant], List[str]]:
    """The element's readable library cases, and one note per case that was skipped."""
    variants, notes = [], []
    for folder in library_elements(roots).get(symbol, []):
        try:
            case = ri.read_case_folder(folder)
        except (OSError, ValueError, KeyError) as exc:
            notes.append(f"{folder.name} skipped: {exc}")
            continue
        if round(case.ae.z) != AtomicData.atomic_numbers[symbol]:
            notes.append(f"{folder.name} skipped: z = {case.ae.z:g} is not {symbol}")
            continue
        variants.append(_describe(folder, case))
    return variants, notes


def rank_variants(variants: Sequence[LibraryVariant]) -> List[LibraryVariant]:
    """PBE first, then PW91, other GGA, LDA; then scalar-relativistic, most tests, folder name."""
    def rank(v: LibraryVariant):
        order = (FUNCTIONAL_PREFERENCE.index(v.exfact) if v.exfact in FUNCTIONAL_PREFERENCE
                 else len(FUNCTIONAL_PREFERENCE))
        return order, v.irel != 2, -v.n_tests, v.folder.name
    return sorted(variants, key=rank)
