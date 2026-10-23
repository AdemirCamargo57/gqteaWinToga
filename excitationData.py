"""Excited-state configuration coefficients for the NTO tab.

Pure Python + numpy, no Toga. Reads the two kinds of file Multiwfn accepts
for its electron excitation analyses (manual section 3.21.A.3), with
Multiwfn's rules (excittrans.f90, loadallexcinfo / loadexccoeff):

* a Gaussian output (.log/.out) of a CIS, TDA or TD job. Only the *last*
  "Excitation energies and oscillator strengths:" block is read (an excited-
  state optimisation prints one per geometry). A state's configurations are
  the "i -> a" / "i <- a" lines right after its header. Gaussian prints only
  coefficients above 0.1 unless the route has IOp(9/40=4) (above 1e-4).
* Multiwfn's plain-text format: "Excited State <n> <multiplicity> <eV>"
  headers followed by the same configuration lines; blank lines separate states.

Orbital numbers are 1-based, as printed. In an unrestricted file each index
carries A (alpha) or B (beta). Closed-shell coefficients follow Gaussian's
convention: sum of squares 0.5, not 1.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

SPIN_ALPHA, SPIN_BETA = "alpha", "beta"
GEOMETRY_TOLERANCE = 1e-4  # Angstrom; Gaussian prints 6 decimals
BOHR_RADIUS_ANGSTROM = 0.529177210544  # coordinateConverter imports Toga, so the constant is repeated

GAUSSIAN_BLOCK = "Excitation energies and oscillator strengths:"
GAUSSIAN_MARKERS = ("Entering Gaussian System", "Gaussian, Inc.")
GAUSSIAN_STATE = re.compile(r"^\s*Excited State\s+(\d+):\s+(\S+)\s+([-+]?\d+\.\d+)\s+eV\s+([-+]?\d+\.\d+)\s+nm"
                            r"\s+f=\s*([-+]?\d+\.\d+)(?:\s+<S\*\*2>=\s*([-+]?\d+\.\d+))?")
TEXT_STATE = re.compile(r"^\s*Excited State\s+(\d+)\s+(-?\d+)\s+([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)\s*$")
# "10 ->103  0.00017": Gaussian drops the space before a three-digit index.
CONFIGURATION = re.compile(r"^\s*(\d+)([AB]?)\s*(->|<-)\s*(\d+)([AB]?)\s+([-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?)\s*$")
DIPOLE_TABLE = "Ground to excited state transition electric dipole moments (Au):"
DIPOLE_ROW = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
NBASIS = re.compile(r"NBasis=\s*(\d+)")
NBSUSE = re.compile(r"NBsUse=\s*(\d+)")
ELECTRONS = re.compile(r"(\d+)\s+alpha electrons\s+(\d+)\s+beta electrons")
MULTIPLICITY_NAMES = {1: "Singlet", 3: "Triplet"}


@dataclass(frozen=True)
class Configuration:
    """One orbital pair: occupied -> virtual (excitation) or occupied <- virtual (de-excitation)."""
    occupied: int        # 1-based MO number
    virtual: int
    spin: str            # SPIN_ALPHA (also for a restricted reference) or SPIN_BETA
    deexcitation: bool   # the "<-" lines of TDHF/TD-DFT
    coefficient: float
    line: int            # 1-based line in the file


@dataclass
class ExcitedState:
    number: int                    # as numbered in the file (1-based)
    label: str                     # 'Singlet-A"', '2.004-E"', 'Singlet' (plain text) ...
    multiplicity: int              # 1, 3, ... ; 0 = undefined (unrestricted reference)
    energy_ev: float
    wavelength_nm: Optional[float]
    oscillator: Optional[float]
    s_squared: Optional[float]
    configurations: Tuple[Configuration, ...]

    def norm_terms(self) -> Tuple[float, float]:
        """(sum of squared excitation coefficients, sum of squared de-excitation coefficients)."""
        x2 = sum(c.coefficient ** 2 for c in self.configurations if not c.deexcitation)
        y2 = sum(c.coefficient ** 2 for c in self.configurations if c.deexcitation)
        return x2, y2

    def title(self) -> str:
        text = f"{self.number}: {self.label} {self.energy_ev:.4f} eV"
        if self.wavelength_nm is not None:
            text += f" {self.wavelength_nm:.2f} nm"
        if self.oscillator is not None:
            text += f" f={self.oscillator:.4f}"
        return text


@dataclass
class ExcitationFile:
    path: str
    source: str                        # "gaussian" | "text"
    states: List[ExcitedState]
    unrestricted: bool                 # indices carry A/B
    blocks: int = 1                    # excitation blocks found (the last is read)
    n_basis: Optional[int] = None      # Gaussian only, for the consistency checks
    n_independent: Optional[int] = None
    n_alpha: Optional[int] = None
    n_beta: Optional[int] = None
    atomic_numbers: Optional[np.ndarray] = None
    coordinates_angstrom: Optional[np.ndarray] = None
    transition_dipoles: Dict[int, Tuple[float, float, float]] = field(default_factory=dict)  # a.u., Gaussian only

    def state(self, number: int) -> ExcitedState:
        for s in self.states:
            if s.number == number:
                return s
        raise ValueError(f"{Path(self.path).name} has no excited state {number}; "
                         f"its states are 1-{len(self.states)}.")


# ----------------------------------------------------------------------
# Parsing
# ----------------------------------------------------------------------
def _number(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def _configuration(match, line_number: int, name: str) -> Configuration:
    occ, occ_spin, arrow, vir, vir_spin, coef = match.groups()
    if occ_spin != vir_spin:
        if occ_spin and vir_spin:
            raise ValueError(f"{name}, line {line_number}: {occ}{occ_spin} {arrow} {vir}{vir_spin} is a "
                             "spin-flip pair; spin-flip excitations are not supported.")
        raise ValueError(f"{name}, line {line_number}: only one side of the pair carries an A/B spin label.")
    spin = SPIN_BETA if occ_spin == "B" else SPIN_ALPHA
    return Configuration(int(occ), int(vir), spin, arrow == "<-", _number(coef), line_number)


def _read_states(lines: List[str], start: int, header, make_state, name: str) -> Tuple[List[ExcitedState], bool]:
    """States from ``lines[start:]``: each header followed by its consecutive configuration lines,
    and whether the configuration lines carry A/B spin labels (an unrestricted reference).

    A header numbered other than the next expected number ends the list (in a
    Gaussian output, a later summary of state 1 must not start a new list).
    """
    states: List[ExcitedState] = []
    labelled = set()
    n = start
    while n < len(lines):
        match = header.match(lines[n])
        n += 1
        if not match:
            continue
        if int(match.group(1)) != len(states) + 1:
            break
        configs: List[Configuration] = []
        while n < len(lines):
            line_match = CONFIGURATION.match(lines[n])
            if not line_match:
                break
            configs.append(_configuration(line_match, n + 1, name))
            labelled.add(bool(line_match.group(2)))
            n += 1
        states.append(make_state(match, tuple(configs)))
    for s in states:
        if not s.configurations:
            raise ValueError(f"{name}: excited state {s.number} lists no configuration coefficients.")
        seen = set()
        for c in s.configurations:
            key = (c.occupied, c.virtual, c.spin, c.deexcitation)
            if key in seen:
                raise ValueError(f"{name}, line {c.line}: the pair {c.occupied} -> {c.virtual} appears twice "
                                 f"in excited state {s.number}.")
            seen.add(key)
    if len(labelled) > 1:
        raise ValueError(f"{name} mixes configuration lines with and without A/B spin labels.")
    return states, labelled == {True}


def _gaussian_state(match, configs) -> ExcitedState:
    number, label, ev, nm, f, s2 = match.groups()
    multiplicity = 1 if label.startswith("Singlet") else 3 if label.startswith("Triplet") else 0
    return ExcitedState(int(number), label, multiplicity, float(ev), float(nm), float(f),
                        float(s2) if s2 is not None else None, configs)


def _text_state(match, configs) -> ExcitedState:
    number, multiplicity, ev = match.groups()
    multiplicity = int(multiplicity)
    label = MULTIPLICITY_NAMES.get(multiplicity, f"Multiplicity {multiplicity}" if multiplicity else "State")
    return ExcitedState(int(number), label, multiplicity, _number(ev), None, None, None, configs)


def _last_orientation(lines: List[str]):
    """(atomic numbers, Angstrom coordinates) of the last Standard (else Input) orientation block."""
    for title in ("Standard orientation:", "Input orientation:"):
        starts = [k for k, line in enumerate(lines) if title in line]
        if not starts:
            continue
        k = starts[-1] + 5  # title, dashes, two heading lines, dashes
        z, xyz = [], []
        while k < len(lines) and not lines[k].strip().startswith("---"):
            parts = lines[k].split()
            z.append(int(parts[1]))
            xyz.append([float(v) for v in parts[3:6]])
            k += 1
        return np.array(z, dtype=int), np.array(xyz, dtype=float).reshape(-1, 3)
    return None, None


def _last_dipole_table(lines: List[str]) -> Dict[int, Tuple[float, float, float]]:
    """{state: (X, Y, Z)} of the last ground-to-excited transition electric dipole table (a.u.)."""
    starts = [k for k, line in enumerate(lines) if DIPOLE_TABLE in line]
    table = {}
    if starts:
        for line in lines[starts[-1] + 2:]:  # title, column headings
            match = DIPOLE_ROW.match(line)
            if not match:
                break
            table[int(match.group(1))] = tuple(float(v) for v in match.groups()[1:])
    return table


def _last_int(pattern, text: str, group: int = 1) -> Optional[int]:
    found = pattern.findall(text)
    if not found:
        return None
    last = found[-1]
    return int(last if isinstance(last, str) else last[group - 1])


def parse_excitation_text(text: str, name: str, path: str = "") -> ExcitationFile:
    """Parse a Gaussian output or a plain-text excitation file. ValueError names every problem."""
    lines = text.splitlines()
    if GAUSSIAN_BLOCK in text:
        starts = [k for k, line in enumerate(lines) if GAUSSIAN_BLOCK in line]
        states, unrestricted = _read_states(lines, starts[-1] + 1, GAUSSIAN_STATE, _gaussian_state, name)
        source = "gaussian"
    elif any(marker in text for marker in GAUSSIAN_MARKERS):
        raise ValueError(f"{name} is a Gaussian output without excited states. Run a CIS/TDA/TD job "
                         "(with IOp(9/40=4)) and use its output.")
    else:
        starts = [0]
        states, unrestricted = _read_states(lines, 0, TEXT_STATE, _text_state, name)
        source = "text"
    if not states:
        raise ValueError(f"{name} contains no excited states: it is neither a Gaussian CIS/TDA/TD output "
                         "nor a plain-text file of 'Excited State' blocks.")
    data = ExcitationFile(path=path or name, source=source, states=states,
                          unrestricted=unrestricted, blocks=len(starts))
    if source == "gaussian":
        data.n_basis = _last_int(NBASIS, text)
        data.n_independent = _last_int(NBSUSE, text)
        electrons = ELECTRONS.findall(text)
        if electrons:
            data.n_alpha, data.n_beta = (int(v) for v in electrons[-1])
        data.atomic_numbers, data.coordinates_angstrom = _last_orientation(lines)
        data.transition_dipoles = _last_dipole_table(lines)
    return data


def read_excitation_file(path) -> ExcitationFile:
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ValueError(f"Cannot read {path}: {exc}") from exc
    return parse_excitation_text(text, path.name, str(path))


# ----------------------------------------------------------------------
# Consistency with the formatted checkpoint
# ----------------------------------------------------------------------
def check_against_wavefunction(data: ExcitationFile, wfn) -> List[str]:
    """Refuse an excitation file that does not belong to ``wfn``; notes (never errors) otherwise."""
    name, fchk = Path(data.path).name, Path(wfn.path).name
    if wfn.alpha_mo is None:
        raise ValueError(f"{fchk} has no MO coefficients; use the .fchk of the CIS/TD job.")
    if wfn.method.upper().startswith("RO"):
        raise ValueError(f"{fchk} is a restricted open-shell ({wfn.method}) reference; excited-state analyses "
                         "support restricted closed-shell and unrestricted references only.")
    if data.unrestricted != wfn.is_unrestricted:
        kinds = ("unrestricted (A/B labels)", "restricted") if data.unrestricted else ("restricted", "unrestricted")
        raise ValueError(f"{name} is {kinds[0]} but {fchk} is {kinds[1]}; use files of the same job.")
    if data.atomic_numbers is not None:  # the molecule first: the most telling mismatch
        if data.atomic_numbers.tolist() != [int(z) for z in wfn.atomic_numbers]:
            raise ValueError(f"{name} and {fchk} describe different atoms "
                             f"({len(data.atomic_numbers)} vs {wfn.n_atoms}); use files of the same job.")
        shift = np.linalg.norm(data.coordinates_angstrom - wfn.coordinates * BOHR_RADIUS_ANGSTROM, axis=1)
        worst = int(np.argmax(shift))
        if shift[worst] > GEOMETRY_TOLERANCE:
            raise ValueError(f"{name} and {fchk} have different geometries (atom {worst + 1} differs by "
                             f"{shift[worst]:.2e} A). Use the .fchk of the same job.")
    for label, a, b in (("basis functions", data.n_basis, wfn.n_basis),
                        ("independent functions", data.n_independent, wfn.n_mo),
                        ("alpha electrons", data.n_alpha, wfn.n_alpha),
                        ("beta electrons", data.n_beta, wfn.n_beta)):
        if a is not None and a != b:
            raise ValueError(f"{name} has {a} {label} but {fchk} has {b}; use files of the same job.")
    occupied ={SPIN_ALPHA: wfn.n_alpha, SPIN_BETA: wfn.n_beta}
    for s in data.states:
        for c in s.configurations:
            n_occ = occupied[c.spin]
            if not (1 <= c.occupied <= n_occ < c.virtual <= wfn.n_mo):
                raise ValueError(f"{name}, excited state {s.number}, line {c.line}: {c.occupied} -> {c.virtual} is "
                                 f"not an occupied -> virtual pair of {fchk} ({n_occ} occupied {c.spin} orbitals "
                                 f"of {wfn.n_mo}).")
    notes = []
    if data.blocks > 1:
        notes.append(f"{name} holds {data.blocks} excitation blocks (e.g. an excited-state optimisation); "
                     "the last one is used, as Multiwfn does.")
    return notes


def sibling_excitation_file(fchk_path) -> Optional[Path]:
    """The Gaussian output beside an .fchk: <stem>.log, else <stem>.out (Multiwfn's ENTER default)."""
    fchk_path = Path(fchk_path)
    for suffix in (".log", ".out"):
        candidate = fchk_path.with_suffix(suffix)
        if candidate.is_file():
            return candidate
    return None
