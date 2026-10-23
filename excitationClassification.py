"""Descriptor-based classification of an electronic excitation, with its sources.

Pure Python + numpy, no Toga. The rules summarise the hole-electron descriptors
the way Multiwfn's manual does (section 3.21.1.1, table on p. 245: LE, single-
direction CT, centrosymmetric CT, Rydberg) with explicit thresholds. Only D/Delta_r
(Guido et al. 2013), Sr/Lambda (Peach et al. 2008) and the sign of t (Multiwfn
manual) have literature values; the Delta sigma and diffuse-electron thresholds are
this program's heuristics and are labelled as such wherever they are reported. When
the descriptors do not fit one pattern the verdict says so instead of forcing a type.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Optional

import numpy as np

CT_DISTANCE = 2.0             # Angstrom: Delta_r criterion of Guido et al., applied to D as well
LOW_OVERLAP = 0.3             # Lambda value of Peach et al., applied to Sr as well
DIFFUSE_SPREAD = 1.0          # Angstrom, Delta sigma: heuristic of this program
DIFFUSE_ELECTRON = 0.5        # fraction of the electron on diffuse basis functions: heuristic of this program
NORMALIZATION_TOLERANCE = 0.01
IFCT_CT_MAJORITY = 50.0       # percent
DIFFUSE_EXPONENT = 0.05       # bohr^-2: a shell whose smallest exponent is below this is diffuse (non-H)
DIFFUSE_EXPONENT_H = 0.03     # idem for hydrogen

VERDICT_INCOMPLETE = "Undetermined (incomplete data)"
VERDICT_CT = "Single-direction charge transfer (CT)"
VERDICT_RYDBERG = "Rydberg-like excitation"
VERDICT_MULTI_CT = "Possible centrosymmetric / multidirectional CT"
VERDICT_LE = "Local excitation (LE)"
VERDICT_INTERMEDIATE = "Intermediate / undetermined"

REFERENCES = {
    "multiwfn": "T. Lu, F. Chen, Multiwfn: A Multifunctional Wavefunction Analyzer, J. Comput. Chem. 33, 580-592 "
                "(2012), DOI 10.1002/jcc.22885; T. Lu, A comprehensive electron wavefunction analysis toolbox for "
                "chemists, Multiwfn, J. Chem. Phys. 161, 082503 (2024), DOI 10.1063/5.0216272.",
    "hole_electron": "Z. Liu, T. Lu, Q. Chen, Carbon 165, 461-467 (2020), DOI 10.1016/j.carbon.2020.05.023 "
                     "(hole-electron analysis, as Multiwfn asks it to be cited).",
    "le_bahers": "T. Le Bahers, C. Adamo, I. Ciofini, J. Chem. Theory Comput. 7, 2498-2506 (2011) (D, RMSD, H, t "
                 "and C functions, introduced for density differences).",
    "delta_r": "C. A. Guido, P. Cortona, B. Mennucci, C. Adamo, J. Chem. Theory Comput. 9, 3118-3126 (2013) "
               "(Delta_r index; 2.0 A as the LE/CT boundary).",
    "lambda": "M. J. G. Peach, P. Benfield, T. Helgaker, D. J. Tozer, J. Chem. Phys. 128, 044118 (2008) "
              "(Lambda index; low values flag excitations that GGA/B3LYP describe poorly).",
    "ghost": "M. Campetella, F. Maschietto, M. J. Frisch, G. Scalmani, I. Ciofini, C. Adamo, J. Comput. Chem. 38, "
             "2151-2156 (2017) (ghost-hunter index).",
    "ghost_fix": "J. Chem. Phys. 154, 204102 (2021) (corrected form of the ghost-hunter index used by Multiwfn).",
    "ifct": "Interfragment charge transfer (IFCT): Multiwfn manual, section 3.21.8 (no separate publication; "
            "cite Multiwfn).",
    "hirshfeld": "F. L. Hirshfeld, Theor. Chim. Acta 44, 129-138 (1977) (stockholder partition).",
    "becke": "A. D. Becke, J. Chem. Phys. 88, 2547-2553 (1988) (multicentre numerical integration).",
    "mulliken": "R. S. Mulliken, J. Chem. Phys. 23, 1833-1840 (1955) (population analysis).",
    "exciton": "Hole-electron Coulomb attraction (exciton binding): Multiwfn manual 3.21.1.1, Theory 6; see "
               "J. Chem. Phys. 143, 244905 (2015) for a discussion.",
}


def references_for(keys: Iterable[str]) -> List[str]:
    """The reference texts of ``keys``, in REFERENCES order, without duplicates."""
    wanted = set(keys)
    return [text for key, text in REFERENCES.items() if key in wanted]


# Functionals with full or high exact exchange at long range (Multiwfn's ghost-state comment names wB97XD,
# CAM-B3LYP and M06-2X); matched against the fchk method name (e.g. RCAM-B3LYP, UWB97XD).
LONG_RANGE_SAFE = ("CAM-", "WB97", "LC-", "LCW", "M11", "M062X", "M06HF", "BHANDHLYP", "HF", "CIS")


def long_range_safe(method: str) -> bool:
    """True when ``method`` (fchk line 2, e.g. 'RCAM-B3LYP') handles long-range CT: range-separated or
    high-exchange hybrids, Hartree-Fock/CIS."""
    name = method.upper().lstrip()
    if name[:1] in ("R", "U") and not name.startswith(("RO",)):
        name = name[1:]
    return any(name.startswith(tag) for tag in LONG_RANGE_SAFE)


@dataclass(frozen=True)
class Descriptors:
    """What the rules look at. Lengths in Angstrom."""
    d: float
    t: float
    sr: float
    delta_sigma: float
    delta_r: Optional[float]
    lambda_: Optional[float]
    normalization_deviation: float
    hole_integral: float
    electron_integral: float
    has_diffuse: bool
    diffuse_electron: float          # fraction of the electron on diffuse basis functions (Mulliken)
    energy_ev: float
    ghost_ev: Optional[float] = None
    method: str = ""                 # fchk method, e.g. 'RCAM-B3LYP' (warnings only)


@dataclass
class Classification:
    verdict: str
    rule: str
    evidence: List[str]
    checks: List[str]
    confidence: str
    fragment_assignment: Optional[str]
    orbital_note: str
    warnings: List[str] = field(default_factory=list)


ORBITAL_NOTE = ("The orbital character of the transition (pi->pi*, n->pi*, d-d, MLCT, LMCT, ...) is not assigned: "
                "it needs the hole and electron isosurfaces (or the NTOs) and the chemical context. The dominant MO "
                "pairs and the atoms carrying hole and electron are listed above as evidence; metal-to-ligand "
                "character can be quantified by defining the metal as a fragment (IFCT).")


def has_diffuse_functions(wfn) -> bool:
    return bool(diffuse_function_mask(wfn).any())


def diffuse_function_mask(wfn) -> np.ndarray:
    """True for every fchk basis function of a shell whose smallest exponent is diffuse."""
    mask = []
    for s in wfn.shells:
        limit = DIFFUSE_EXPONENT_H if int(wfn.atomic_numbers[s.atom]) == 1 else DIFFUSE_EXPONENT
        mask += [bool(s.exponents.min() < limit)] * s.n_functions
    return np.array(mask, dtype=bool)


def _f(value: float, digits: int = 2) -> str:
    return f"{value:.{digits}f}"


def classify(d: Descriptors, ifct=None) -> Classification:
    """Apply the rules in order (first match wins); ``ifct`` is an excitationAnalysis.IFCTResult or None."""
    evidence = [
        f"D = {_f(d.d)} A ({'>=' if d.d >= CT_DISTANCE else '<'} {CT_DISTANCE} A: Guido et al. 2013 Delta_r "
        "criterion, applied to D, which equals Delta_r for a single orbital pair)",
        f"t = {_f(d.t)} A ({'> 0: hole and electron separated' if d.t > 0 else '<= 0: hole and electron not separated'}"
        " along the CT direction; Multiwfn manual)",
        f"Sr = {_f(d.sr)} ({'>=' if d.sr >= LOW_OVERLAP else '<'} {LOW_OVERLAP}: Peach et al. 2008 Lambda value, "
        "applied to Sr, which equals Lambda for a single orbital pair)",
        f"Delta sigma = {_f(d.delta_sigma)} A ({'>=' if abs(d.delta_sigma) >= DIFFUSE_SPREAD else '<'} "
        f"{DIFFUSE_SPREAD} A in magnitude: heuristic threshold of this program)",
    ]
    if d.has_diffuse:
        evidence.append(f"electron on diffuse basis functions = {100 * d.diffuse_electron:.0f} % "
                        f"({'>=' if d.diffuse_electron >= DIFFUSE_ELECTRON else '<'} {100 * DIFFUSE_ELECTRON:.0f} %: "
                        "heuristic threshold of this program; Mulliken population)")
    warnings: List[str] = []
    separated = d.t > 0
    far = d.d >= CT_DISTANCE
    spread = d.delta_sigma >= DIFFUSE_SPREAD

    incomplete = [label for label, bad in (
        ("configuration coefficients", abs(d.normalization_deviation) > NORMALIZATION_TOLERANCE),
        ("hole integral", abs(d.hole_integral - 1.0) > NORMALIZATION_TOLERANCE),
        ("electron integral", abs(d.electron_integral - 1.0) > NORMALIZATION_TOLERANCE)) if bad]
    if incomplete:
        verdict, rule = VERDICT_INCOMPLETE, (
            "Rule 1: " + ", ".join(incomplete) + f" deviate from their expected values by more than "
            f"{NORMALIZATION_TOLERANCE}; the indices are not reliable (rerun Gaussian with IOp(9/40=4)).")
    elif far and separated:
        verdict, rule = VERDICT_CT, "Rule 2: D >= 2.0 A and t > 0."
    elif spread and not far and d.has_diffuse and d.diffuse_electron >= DIFFUSE_ELECTRON:
        verdict, rule = VERDICT_RYDBERG, ("Rule 3: Delta sigma >= 1.0 A, D < 2.0 A and most of the electron on "
                                          "diffuse basis functions (the manual: the electron consists of very "
                                          "diffuse MOs). Charge-transfer character is not excluded by this rule.")
    elif spread and not far and not separated and d.sr >= LOW_OVERLAP:
        verdict, rule = VERDICT_MULTI_CT, ("Rule 4: D < 2.0 A, t <= 0, Delta sigma >= 1.0 A and Sr >= 0.3 "
                                           "(manual: D vanishes for centrosymmetric CT while Delta sigma is large); "
                                           "needs fragment confirmation (IFCT).")
    elif not far and not separated and abs(d.delta_sigma) < DIFFUSE_SPREAD and d.sr >= LOW_OVERLAP:
        verdict, rule = VERDICT_LE, "Rule 5: D < 2.0 A, t <= 0, |Delta sigma| < 1.0 A and Sr >= 0.3."
    else:
        verdict, rule = VERDICT_INTERMEDIATE, ("Rule 6: the descriptors do not fit one pattern ("
                                               + _disagreement(d) + "); no type is assigned.")
    if spread and not far and not d.has_diffuse:
        warnings.append("The electron is much more spread out than the hole, but the basis set has no diffuse "
                        "functions: Rydberg character cannot be assigned in this basis.")

    checks, disagree = [], []
    if d.delta_r is not None:
        ct_like = d.delta_r >= CT_DISTANCE
        checks.append(f"Delta_r = {_f(d.delta_r)} A ({'CT-like' if ct_like else 'local'} by the 2.0 A criterion)")
        if verdict in (VERDICT_CT, VERDICT_LE):
            if ct_like != (verdict == VERDICT_CT):
                disagree.append("Delta_r")
    if d.lambda_ is not None:
        low = d.lambda_ < LOW_OVERLAP
        checks.append(f"Lambda = {_f(d.lambda_)} ({'low overlap' if low else 'substantial overlap'} by the 0.3 "
                      "value)")
        if low and verdict == VERDICT_LE:  # a low Lambda flags CT; a CT state need not have one
            disagree.append("Lambda")
    if verdict in (VERDICT_INCOMPLETE, VERDICT_INTERMEDIATE):
        confidence = "none: no type assigned"
    elif disagree:
        confidence = (f"low: the orbital-based {' and '.join(disagree)} disagree(s); Multiwfn's author considers D "
                      "and Sr more reliable because they include the coupling between configurations")
    elif verdict in (VERDICT_MULTI_CT, VERDICT_RYDBERG):
        confidence = "tentative: confirm with the hole/electron isosurfaces" + (
            " and fragments (IFCT)" if verdict == VERDICT_MULTI_CT else "")
    else:
        confidence = "consistent: the independent indices agree"

    safe = long_range_safe(d.method)
    if d.ghost_ev is not None and d.energy_ev < d.ghost_ev:
        text = (f"The excitation energy ({d.energy_ev:.3f} eV) is below the ghost-hunter index ({d.ghost_ev:.3f} eV): "
                "possibly a spurious (ghost) CT state. This is a necessary, not a sufficient condition (Multiwfn "
                "manual 3.21.7), meaningful only for one-dimensional CT.")
        if safe:
            text += (f" With {d.method}, which has full or high exact exchange at long range, ghost states are not "
                     "expected (Multiwfn's comment).")
        warnings.append(text)
    if (d.delta_r is not None and d.delta_r >= CT_DISTANCE) or (d.lambda_ is not None and d.lambda_ < LOW_OVERLAP):
        text = ("Long-range or low-overlap excitation (Delta_r >= 2.0 A or Lambda < 0.3): pure and low-exact-"
                "exchange functionals (BLYP, PBE, B3LYP, PBE0) are unreliable here; a range-separated hybrid such "
                "as CAM-B3LYP or wB97XD is needed (manual 3.21.4; Peach et al.).")
        text += (f" The method used ({d.method}) is of that kind." if safe else
                 f" The method used ({d.method or 'unknown'}) is not recognised as one." )
        warnings.append(text)
    return Classification(verdict, rule, evidence, checks, confidence, _fragment_assignment(ifct), ORBITAL_NOTE,
                          warnings)


def _disagreement(d: Descriptors) -> str:
    parts = []
    if d.d < CT_DISTANCE and d.t > 0:
        parts.append("t > 0 although D < 2.0 A: partial charge transfer")
    if d.d >= CT_DISTANCE and d.t <= 0:
        parts.append("D >= 2.0 A although t <= 0: hole and electron are broad compared with their separation")
    if d.sr < LOW_OVERLAP and d.delta_sigma < DIFFUSE_SPREAD:
        parts.append("small Sr without a diffuse electron")
    if d.delta_sigma >= DIFFUSE_SPREAD and d.sr < LOW_OVERLAP and not (d.has_diffuse and
                                                                       d.diffuse_electron >= DIFFUSE_ELECTRON):
        parts.append("spread-out electron with small Sr but not on diffuse functions")
    if d.delta_sigma <= -DIFFUSE_SPREAD:
        parts.append("the hole is much more spread out than the electron")
    return "; ".join(parts) or "mixed indicators"


def _fragment_assignment(ifct) -> Optional[str]:
    if ifct is None or len(ifct.fragments) < 2:
        return None
    names = [f.name for f in ifct.fragments]
    if ifct.intrinsic_ct >= IFCT_CT_MAJORITY:
        r, s, net = ifct.net_transfers()[0]
        return (f"Charge transfer between fragments: intrinsic CT {ifct.intrinsic_ct:.1f} % >= "
                f"{IFCT_CT_MAJORITY:.0f} %; largest net transfer {names[r]} -> {names[s]} ({net:.3f} e). "
                "This depends on the fragment definition.")
    k = int(np.argmax(np.diag(ifct.matrix)))
    return (f"Local excitation with respect to the fragments: intrinsic LE {ifct.intrinsic_le:.1f} %; mainly within "
            f"{names[k]} ({100 * ifct.matrix[k, k]:.1f} % redistributed there). This depends on the fragment "
            "definition.")
