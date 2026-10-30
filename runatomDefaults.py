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
from typing import Dict, List, Optional, Sequence, Tuple

import runatomInput as ri
from help import AtomicData
from periodicTableData import element_block, group_and_period
from runatomStarterData import STARTER_EIGENVALUES, STARTER_SETTINGS, STARTER_XION, UNCONVERGED

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


# --- computed starter: configuration, core/valence split, eigenvalue guesses ---
MADELUNG_ORDER = sorted(((n, l) for n in range(1, 8) for l in range(min(n, 4))),
                        key=lambda nl: (sum(nl), nl[0]))
SLATER_N_STAR = {1: 1.0, 2: 2.0, 3: 3.0, 4: 3.7, 5: 4.0, 6: 4.2, 7: 4.4}   # n = 7 extrapolated
SYMBOLS_BY_Z = {z: s for s, z in AtomicData.atomic_numbers.items()}


def element_symbol(z: int) -> str:
    return SYMBOLS_BY_Z[int(round(z))]


def _checked_z(symbol: str) -> int:
    z = AtomicData.atomic_numbers.get(symbol)
    if z is None:
        raise ValueError(f"{symbol!r} is not an element symbol")
    if z > MAX_DEFAULTS_Z:
        raise ValueError(f"{symbol}: Z = {z} > {MAX_DEFAULTS_Z}; suggested values are available up to uranium")
    return z


def ground_state_configuration(z: int) -> List[Tuple[int, int, float]]:
    """(n, l, electrons) of the neutral atom in filling order: Madelung's rule,
    then AtomicData.ground_state_exceptions; empty shells left out."""
    occupation, left = {}, z
    for n, l in MADELUNG_ORDER:
        occupation[(n, l)] = min(left, 2 * (2 * l + 1))
        left -= occupation[(n, l)]
    occupation.update(AtomicData.ground_state_exceptions.get(element_symbol(z), {}))
    shells = [(n, l, float(occupation[(n, l)])) for n, l in MADELUNG_ORDER if occupation[(n, l)] > 0]
    assert sum(e for _, _, e in shells) == z, f"Z = {z}: configuration does not add up"
    return shells


def _slater_group(n: int, l: int):
    return (n, "sp") if l <= 1 else (n, l)


# Slater's order: (1s)(2s,2p)(3s,3p)(3d)(4s,4p)(4d)(4f)(5s,5p)(5d)(5f)...
SLATER_GROUPS = sorted(set(_slater_group(n, l) for n, l in MADELUNG_ORDER),
                       key=lambda g: (g[0], 0 if g[1] == "sp" else g[1]))


def eigenvalue_guess(z: int, config, n: int, l: int) -> float:
    """Slater-rule hydrogenic estimate of the (n, l) eigenvalue in Ry: -(Z_eff / n*)**2."""
    group = _slater_group(n, l)
    position = SLATER_GROUPS.index(group)
    screening = 0.0
    for m, k, electrons in config:
        other = _slater_group(m, k)
        if other == group:
            screening += (electrons - 1 if (m, k) == (n, l) else electrons) * (0.30 if n == 1 else 0.35)
        elif SLATER_GROUPS.index(other) < position:
            if l >= 2 or m <= n - 2:
                screening += electrons
            elif m == n - 1:
                screening += 0.85 * electrons
    z_eff = max(z - screening, 1.0)
    return -(z_eff / SLATER_N_STAR[n]) ** 2


def _adds_semicore(symbol: str) -> bool:
    """Groups 1-2 from period 3 and groups 3-6 from period 4 carry (n-1)s and (n-1)p in the valence."""
    group, period = group_and_period(symbol)
    if element_block(symbol) == "f":
        return False
    return (group <= 2 and period >= 3) or (3 <= group <= 6 and period >= 4)


def split_core_valence(z: int, config) -> Tuple[list, list, List[str]]:
    """Core and valence shells (each in filling order) and notes.

    Valence: the outermost ns and np, the (n-1)d of the d block (group 12
    included; the p block keeps its filled d in the core), the (n-2)f of the f
    block, plus semicore (n-1)s and (n-1)p where _adds_semicore says so.
    """
    symbol = element_symbol(z)
    group, period = group_and_period(symbol)
    block = element_block(symbol)
    semicore = _adds_semicore(symbol)

    def is_valence(n, l):
        if n == period and l <= 1:
            return True
        if l == 2 and n == period - 1:
            return group <= 12          # d and f block (f counts as group 3); the p block keeps it in the core
        if l == 3 and n == period - 2:
            return block == "f"
        return semicore and n == period - 1 and l <= 1

    valence = [s for s in config if is_valence(s[0], s[1])]
    core = [s for s in config if s not in valence]
    notes = []
    while len(valence) > ri.MAX_VALENCE:
        moved = valence.pop(0)
        core = sorted(core + [moved], key=lambda s: MADELUNG_ORDER.index((s[0], s[1])))
        notes.append(f"{symbol}: {shell_name(moved[0] * 100 + moved[1] * 10)} moved to the core "
                     f"(runatom allows {ri.MAX_VALENCE} valence shells).")
    return core, valence, notes


# The 000-Ex template's printing and convergence controls for a generation run.
STARTER_PS_CONTROLS = {"ps_print": "1  1  3", "ps_logd": "2.001  -6  2  80", "ps_conv": "1e-10  1e-09  0.5  0"}
# No log derivatives (ilogd 0), but rinit still needs rlogd on a usable mesh point
# (runatomInput.log_derivative_index); thresh, tol, damp, maxit 0 = program defaults.
STARTER_AE_CONTROLS = ri.Controls(1, 0, 0, 0, 0, 2.0, -2.0, 2.0, 40, 0.0, 0.0, 0.0, 0, "")
# Mesh and SCF settings of the starter: a larger box than the 80-bohr default for
# diffuse d atoms and stronger mixing damping for open d shells. Atoms that need
# others have them in runatomStarterData.STARTER_SETTINGS (found by the generator,
# checked by tests/test_runatomStarterExe.py).
STARTER_RMAX = 150.0
STARTER_BBSF = 40.0
STARTER_DAMP = 0.2
STARTER_MAXIT = 1000


def _nnlz(shell) -> int:
    return shell[0] * 100 + shell[1] * 10


def _occupation_text(electrons: float) -> str:
    return str(int(electrons)) if float(electrons).is_integer() else f"{electrons:g}"


def starter_settings(z: int) -> Tuple[float, float, float, int]:
    """(rmax, bbsf, damp, maxit) runatom converges this starter atom with."""
    return STARTER_SETTINGS.get(z, (STARTER_RMAX, STARTER_BBSF, STARTER_DAMP, STARTER_MAXIT))


def ionised(config, xion: float) -> List[Tuple[int, int, float]]:
    """The configuration with xion electrons taken from the outermost s shell."""
    if not xion:
        return list(config)
    outer = max((s for s in config if s[1] == 0), key=lambda s: s[0])
    return [(n, l, e - xion) if (n, l, e) == outer else (n, l, e) for n, l, e in config]


def starter_ae_input(symbol: str, guesses: Optional[Dict[int, float]] = None,
                     rmax: Optional[float] = None, bbsf: Optional[float] = None,
                     damp: Optional[float] = None, maxit: Optional[int] = None,
                     xion: Optional[float] = None) -> ri.AEInput:
    """The starter atom, PBE, scalar-relativistic: neutral, or the +1 ion where
    runatom does not converge the neutral atom (runatomStarterData.STARTER_XION).

    Eigenvalue guesses (nnlz -> Ry) default to runatom's own converged values
    (runatomStarterData, generated by tools/make_runatom_starter_data.py);
    a shell missing there gets Slater's estimate. Settings left None come
    from starter_settings and STARTER_XION.
    """
    z = _checked_z(symbol)
    xion = STARTER_XION.get(z, 0.0) if xion is None else xion
    neutral = ground_state_configuration(z)
    core, valence, _ = split_core_valence(z, neutral)
    config = ionised(core + valence, xion)
    known = STARTER_EIGENVALUES.get(z, {}) if guesses is None else guesses
    shells = [ri.Shell(_nnlz(s), s[2], known.get(_nnlz(s), round(eigenvalue_guess(z, config, s[0], s[1]), 4)))
              for s in config]
    own = starter_settings(z)
    rmax, bbsf, damp, maxit = [own[k] if v is None else v for k, v in enumerate((rmax, bbsf, damp, maxit))]
    controls = ri.replace(STARTER_AE_CONTROLS, title=symbol, damp=damp, maxit=maxit)
    return ri.AEInput(controls, float(z), xion, PBE, rmax, 0.0, bbsf, 2, shells)


def generation_tag(symbol: str) -> str:
    """The outer valence shells in the library's style: 's2' (Mg), 's2d2' (Ti)."""
    z = _checked_z(symbol)
    core, valence, _ = split_core_valence(z, ground_state_configuration(z))
    semicore = _semicore_shells(symbol, valence)
    starter = ionised(core + valence, STARTER_XION.get(z, 0.0))[len(core):]
    outer = [s for s, neutral in zip(starter, valence) if neutral not in semicore]
    return "".join(f"{L_LETTERS[l]}{_occupation_text(e)}" for _, l, e in outer)


def _semicore_shells(symbol: str, valence) -> list:
    _, period = group_and_period(symbol)
    return [s for s in valence if s[0] == period - 1 and s[1] <= 1] if _adds_semicore(symbol) else []


def starter_fields(symbol: str) -> Dict[str, str]:
    """Every form field of the computed starter; the pseudization choices left blank."""
    z = _checked_z(symbol)
    ae = starter_ae_input(symbol)
    _, valence, _ = split_core_valence(z, ground_state_configuration(z))
    fields = {key: "" for key in ri.FORM_FIELDS}
    fields.update(ri.ae_fields(ae))
    fields.update(STARTER_PS_CONTROLS)
    nang = max(l for _, l, _ in valence) + 1
    fields.update(symbol=symbol, generation=generation_tag(symbol), title=symbol,
                  ps_partition=f"{len(ae.shells) - len(valence)}  {len(valence)}  {nang}")
    return fields


def rank_variants(variants: Sequence[LibraryVariant]) -> List[LibraryVariant]:
    """PBE first, then PW91, other GGA, LDA; then scalar-relativistic, most tests, folder name."""
    def rank(v: LibraryVariant):
        order = (FUNCTIONAL_PREFERENCE.index(v.exfact) if v.exfact in FUNCTIONAL_PREFERENCE
                 else len(FUNCTIONAL_PREFERENCE))
        return order, v.irel != 2, -v.n_tests, v.folder.name
    return sorted(variants, key=rank)


# --- same-group templates and the suggestions the window offers --------------
STARTER_LABEL = "Computed starter (radii blank)"
COPIED_PS_KEYS = ("ps_print", "ps_logd", "ps_conv", "ps_bessel", "ps_keyps", "ps_rinner", "ps_rcloc", "ps_rc")
PROJECTOR_PS_KEYS = ("ps_betas", "ps_npf", "ps_local", "ps_qfunc", "ps_fixes")


@dataclass(frozen=True)
class Suggestion:
    kind: str                    # "library", "template" or "starter"
    label: str                   # the window's Selection item
    fields: Dict[str, str]       # every form field
    folder: Optional[Path]       # the library case, or the template's donor case
    note: str                    # one line for the window
    warnings: Tuple[str, ...]    # shown in a dialog after the fields are filled


def template_donors(symbol: str, roots: Sequence[Path]) -> List[str]:
    """Library elements of the same group, nearest period first, the heavier on a tie."""
    group, period = group_and_period(symbol)
    candidates = [s for s in library_elements(roots)
                  if s != symbol and group_and_period(s)[0] == group]
    return sorted(candidates, key=lambda s: (abs(group_and_period(s)[1] - period),
                                             -AtomicData.atomic_numbers[s]))


def _shifted(nnlz: int, periods: int) -> int:
    """The same shell n periods away; a shell that cannot exist (l >= n) becomes the lowest one of that l."""
    n, l = nnlz // 100 + periods, (nnlz // 10) % 10
    return (max(n, l + 1) * 10 + l) * 10


@dataclass(frozen=True)
class _ValenceMap:
    ae: ri.AEInput                  # the starter's shells, core first, valence in the donor's order
    ncores: int
    index: Dict[int, int]           # donor valence position -> template position (1-based)
    dropped: Dict[int, Tuple[int, float]]   # donor position -> (shifted nnlz, donor's ee) of dropped shells


def _template_ae(symbol: str, donor: ri.CaseInputs) -> Tuple[Optional[_ValenceMap], str]:
    """The starter atom with the donor's valence list, n shifted by the period difference.

    The template keeps exactly the starter's shells (whose all-electron run is known
    to converge); only which of them are valence, and in what order, follows the
    donor. A donor valence shell the target holds no electron in (Ca's empty 3d for
    Mg) is dropped: it need not be bound in the target. Returns (map, "") or
    (None, the reason the donor's valence list does not fit).
    """
    starter = starter_ae_input(symbol)
    z = round(starter.z)
    _, valence, _ = split_core_valence(z, ground_state_configuration(z))
    periods = group_and_period(symbol)[1] - group_and_period(element_symbol(donor.ae.z))[1]
    donor_valence = donor.ae.shells[donor.ps.ncores:]
    mapped = [_shifted(s.nnlz, periods) for s in donor_valence]
    if len(set(mapped)) != len(mapped):
        return None, f"the donor's valence shells do not map one to one onto {symbol}'s"
    own = {s.nnlz: s for s in starter.shells}
    semicore = {_nnlz(s) for s in _semicore_shells(symbol, valence)}
    missing = [_nnlz(s) for s in valence if _nnlz(s) not in mapped and _nnlz(s) not in semicore
               and own[_nnlz(s)].wwnl > 0]
    if missing:
        names = " ".join(shell_name(n) for n in missing)
        return None, (f"{symbol}'s occupied {names} valence shell has no counterpart among the donor's "
                      f"valence shells ({' '.join(shell_name(s.nnlz) for s in donor_valence)})")
    index, dropped, kept = {}, {}, []
    for position, (nnlz, shell) in enumerate(zip(mapped, donor_valence), 1):
        if nnlz in own and own[nnlz].wwnl > 0:
            kept.append(own[nnlz])
            index[position] = len(kept)
        else:
            dropped[position] = (nnlz, shell.ee)
    core = [s for s in starter.shells if s not in kept]
    return _ValenceMap(ri.replace(starter, shells=core + kept), len(core), index, dropped), ""


def _remapped_ps(donor: ri.CaseInputs, valence: _ValenceMap) -> ri.PSInput:
    """The donor's generation input with its valence references moved to the template's positions;
    a reference to a dropped shell becomes keyee 0 with that shell's energy from the donor input."""
    def remap(keyee: int, energy: float) -> Tuple[int, float]:
        if keyee <= 0:
            return keyee, energy
        if keyee in valence.dropped:
            return 0, valence.dropped[keyee][1]
        return valence.index[keyee], energy
    betas = [ri.replace(b, keyee=k, eeread=e) for b in donor.ps.betas for k, e in [remap(b.keyee, b.eeread)]]
    keyee_loc, eloc = remap(donor.ps.keyee_loc, donor.ps.eloc)
    return ri.replace(donor.ps, ncores=valence.ncores, nvales=len(valence.ae.shells) - valence.ncores,
                      betas=betas, keyee_loc=keyee_loc, eloc=eloc)


def _psqf2_minimum(ps: ri.PSInput) -> int:
    """The least nqf runatom's psqf2 accepts for ifqopt 3 (see runatomInput.check_ps)."""
    constrained = any(1 <= f.ibeta <= len(ps.betas) and ps.betas[f.ibeta - 1].iptype != 3 for f in ps.fixes)
    return 7 if constrained else 6


def template_fields(symbol: str, donor_folder: Path, donor: ri.CaseInputs) -> Tuple[Dict[str, str], List[str]]:
    """The starter of `symbol` with the pseudopotential tab of a same-group library case."""
    donor_symbol = element_symbol(donor.ae.z)
    periods = group_and_period(symbol)[1] - group_and_period(donor_symbol)[1]
    where = ("the same period" if periods == 0 else
             f"{abs(periods)} period{'s' if abs(periods) > 1 else ''} {'below' if periods > 0 else 'above'}")
    fields = starter_fields(symbol)
    donor_ps = ri.ps_fields(donor.ps)
    fields.update({key: donor_ps[key] for key in COPIED_PS_KEYS})
    valence, reason = _template_ae(symbol, donor)
    copied = "rc, rcloc, rinner, rpcor" + (", the projectors and the local potential" if valence else "")
    warnings = [
        f"{symbol} has no case in the USPP library. The All-electron tab was computed from its "
        f"configuration (PBE, scalar-relativistic). The Pseudopotential tab ({copied}) was copied unchanged "
        f"from {donor_symbol} ({donor_folder.name}, {ri.xc_label(donor.ae.exfact)}, same group, {where}). "
        f"These values were chosen for {donor_symbol}, not {symbol}: they are a starting point, not a tested "
        "pseudopotential. Review them, add transferability tests (Tests tab) and compare with all-electron "
        "results before production use."]
    if valence is None:
        partition = fields["ps_partition"].split()
        fields["ps_partition"] = f"{partition[0]}  {partition[1]}  {donor.ps.nang}"
        warnings.append("Only the radii and controls were copied: " + reason + ", so the projectors, the "
                        "local potential and the Q-function settings were left blank. Enter them on the "
                        "Pseudopotential tab.")
    else:
        ps = _remapped_ps(donor, valence)
        fields.update(ri.ae_fields(valence.ae))
        for position, (nnlz, energy) in sorted(valence.dropped.items()):
            users = [k for k, b in enumerate(donor.ps.betas, 1) if b.keyee == position]
            local = donor.ps.keyee_loc == position
            if users or local:
                what = ", ".join([f"projector {k}" for k in users] + (["the local potential"] if local else []))
                source = shell_name(donor.ae.shells[donor.ps.ncores + position - 1].nnlz)
                warnings.append(f"{donor_symbol}'s {what} took its reference energy from {donor_symbol}'s "
                                f"{source} shell, which is {shell_name(nnlz)} in {symbol} and holds no electron "
                                f"there; it now uses eeread = {energy:g} Ry ({donor_symbol}'s eigenvalue guess "
                                "for that shell). Review it.")
        if ps.ifqopt != 3:
            least = _psqf2_minimum(ps)
            warnings.append(f"ifqopt was set to 3 ({donor_symbol}'s case uses {ps.ifqopt}): runatom accepts "
                            "only ifqopt 3 for GGA functionals such as PBE.")
            if ps.nqf < least:
                warnings.append(f"nqf was raised from {ps.nqf} to {least}, the least runatom accepts "
                                "for ifqopt 3.")
            ps = ri.replace(ps, ifqopt=3, nqf=max(ps.nqf, least))
        remapped = ri.ps_fields(ps)
        fields["ps_partition"] = remapped["ps_partition"]
        fields.update({key: remapped[key] for key in PROJECTOR_PS_KEYS})
    if donor.ae.exfact + 0.0 != PBE:
        warnings.append(f"{donor_symbol}'s case was made for {ri.xc_label(donor.ae.exfact)}; this template is "
                        "PBE, so its radii and reference energies were tuned for another functional.")
    try:
        problems = ri.check_case(ri.case_from_fields(fields))
    except ValueError as exc:
        problems = [str(exc)]
    warnings += [f"Check: {p}" for p in problems]
    return fields, warnings


def _starter_warnings(symbol: str) -> List[str]:
    z = AtomicData.atomic_numbers[symbol]
    if STARTER_XION.get(z):
        return [f"runatom does not converge the neutral {symbol} atom (its 4f state), so the starter is the "
                f"{STARTER_XION[z]:+g} ion with one 6s electron removed, as the library does for Ti. Change xion "
                "and the occupations on the All-electron tab for another configuration."]
    return []


def suggestions(symbol: str, roots: Sequence[Path]) -> Tuple[List[Suggestion], List[str]]:
    """What the window offers for `symbol`, best first, and notes on library cases that were skipped."""
    _checked_z(symbol)
    variants, notes = element_variants(symbol, roots)
    if variants:
        options = []
        for variant in rank_variants(variants):
            warnings = [] if variant.exfact == PBE else [
                f"No PBE case for {symbol} in the library; loaded {variant.folder.name} "
                f"({ri.xc_label(variant.exfact)}). Changing exfact to 5 alone is not enough: GGA needs "
                "ifqopt 3, and the radii were tuned for this functional."]
            options.append(Suggestion("library", variant.label, ri.case_to_fields(ri.read_case_folder(variant.folder)),
                                      variant.folder, f"From {variant.folder.name}. Every field stays editable.",
                                      tuple(warnings)))
        return options, notes

    skipped = [f"Library case {n}" for n in notes]
    starter_warnings = _starter_warnings(symbol)
    options = []
    for donor_symbol in template_donors(symbol, roots):
        donor_variants, _ = element_variants(donor_symbol, roots)
        if not donor_variants:
            continue
        donor = rank_variants(donor_variants)[0]
        fields, warnings = template_fields(symbol, donor.folder, ri.read_case_folder(donor.folder))
        options.append(Suggestion("template", f"Template from {donor.folder.name}", fields, donor.folder,
                                  f"Template from {donor.folder.name}: review before use. Every field stays "
                                  "editable.", tuple(warnings + starter_warnings + skipped)))
        break
    if not library_elements(roots):
        reason = ("The USPP library was not found, so no tested values could be suggested or copied. "
                  "Enter rc, rcloc, rinner, the projectors and the local potential, or import a case.")
    elif options:
        reason = None
    else:
        reason = (f"No element of group {group_and_period(symbol)[0]} is in the USPP library, so no radii "
                  "could be copied. Enter rc, rcloc, rinner, the projectors and the local potential on the "
                  "Pseudopotential tab, or import a case.")
    starter_note = "Computed starter: the Pseudopotential tab is still to be filled. Every field stays editable."
    options.append(Suggestion("starter", STARTER_LABEL, starter_fields(symbol), None, starter_note,
                              tuple(([reason] if reason else []) + starter_warnings + (skipped if not options else []))))
    return options, notes
