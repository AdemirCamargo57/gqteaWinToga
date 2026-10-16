"""Tools > Fukui Function Analysis.

Fukui functions f+(r), f-(r), f0(r) and the dual descriptor Df(r) from three
Gaussian formatted checkpoints of the same molecule at the same geometry and
level of theory: the N-electron state and the N+p / N-q states (p = q = 1
normally; p, q > 1 for (quasi-)degenerate frontier orbitals, Multiwfn manual
section 3.25.4).

    f+  = [rho(N+p) - rho(N)] / p        nucleophilic attack
    f-  = [rho(N) - rho(N-q)] / q        electrophilic attack
    f0  = (f+ + f-) / 2                  radical attack
    Df  = f+ - f-                        dual descriptor (> 0 nucleophilic,
                                         < 0 electrophilic attack)

Real-space functions are computed here, from the exact densities
(gaussianWavefunction). Atom-condensed indices come from Multiwfn main
function 22 option 2 (Hirshfeld charges of the three states), either run here
as an external engine or imported from an existing CDFT.txt -- this module
does not compute Hirshfeld charges itself. Global indices follow Multiwfn's
conventions (hardness = I - A, without the 1/2).

These quantities indicate *likely* reactive regions; they are not
predictions of reaction outcome.
"""
from __future__ import annotations

import asyncio
import os
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import cubeFile as cf
import gaussianWavefunction as gw
import molecularSurface as ms
from coordinateConverter import BOHR_RADIUS_ANGSTROM
from espViewer import ESPSurfaceViewer, SurfaceMarker, UNIT_AU
from help import AtomicData, HelpGqteaWin
from molecularViewer import save_canvas_image

HARTREE_TO_EV = 27.211386245988   # CODATA 2018; reproduces Multiwfn's eV columns
TCE_HOMO_HARTREE = -0.335198      # E_HOMO(TCE), B3LYP/6-31G*, as used by Multiwfn (manual 3.25.2)
DEGENERACY_EV = 0.01              # Multiwfn's automatic FMO degeneracy threshold (manual 4.22.3.1)
DEFAULT_FUKUI_ISOVALUE = 0.005    # a.u.; the manual's examples use 0.005-0.01
DEFAULT_SURFACE_ISOVALUE = 0.01   # rho isosurface for mapping f (manual 4.12.4)
DEFAULT_SPACING_BOHR = 0.25
COORDINATE_TOLERANCE_BOHR = 1e-5
CONDENSED_TOLERANCE = 2e-4        # CDFT.txt prints 4 decimals
MULTIWFN_TIMEOUT_S = 1800

FUNCTIONS = ("f+", "f-", "f0", "DD")
FUNCTION_TITLES = {
    "f+": "f+  (nucleophilic attack)",
    "f-": "f-  (electrophilic attack)",
    "f0": "f0  (radical attack)",
    "DD": "dual descriptor  (+ nucleophilic / - electrophilic)",
}
CUBE_SUFFIXES = {"f+": "_fplus.cub", "f-": "_fminus.cub", "f0": "_fzero.cub", "DD": "_dual.cub"}
REPORT_SUFFIX = "_Fukui.txt"

CONDENSED_MULTIWFN = "Run Multiwfn (Hirshfeld)"
CONDENSED_IMPORT = "Import Multiwfn CDFT.txt"
CONDENSED_NONE = "Skip condensed indices"
CONDENSED_SOURCES = (CONDENSED_MULTIWFN, CONDENSED_IMPORT, CONDENSED_NONE)

ProgressCallback = Callable[[float, str], None]


# ----------------------------------------------------------------------
# Input parsing
# ----------------------------------------------------------------------
def parse_positive_float(text: Optional[str], default: float, name: str, upper: float) -> float:
    """Blank -> default; otherwise a number in (0, upper)."""
    text = (text or "").strip()
    if not text:
        return default
    try:
        value = float(text)
    except ValueError:
        raise ValueError(f"{name}: '{text}' is not a number.") from None
    if not 0.0 < value < upper:
        raise ValueError(f"{name} must be greater than 0 and smaller than {upper:g}.")
    return value


def resolve_output_dir(fchk_path: str, chosen: Optional[str]) -> str:
    chosen = (chosen or "").strip()
    return chosen if chosen else str(Path(fchk_path).resolve().parent)


# ----------------------------------------------------------------------
# Three states: validation
# ----------------------------------------------------------------------
def strip_reference_prefix(method: str) -> str:
    """'RB3LYP', 'UB3LYP', 'ROB3LYP' -> 'B3LYP' (the closed- and open-shell states differ only there)."""
    m = method.upper()
    for prefix in ("RO", "R", "U"):
        if m.startswith(prefix) and len(m) > len(prefix):
            return m[len(prefix):]
    return m


def same_basis(a: gw.FchkWavefunction, b: gw.FchkWavefunction) -> bool:
    if len(a.shells) != len(b.shells):
        return False
    return all(sa.l == sb.l and sa.pure == sb.pure and sa.atom == sb.atom
               and sa.exponents.shape == sb.exponents.shape
               and np.allclose(sa.exponents, sb.exponents, rtol=1e-10, atol=0.0)
               for sa, sb in zip(a.shells, b.shells))


def validate_states(n: gw.FchkWavefunction, plus: gw.FchkWavefunction,
                    minus: gw.FchkWavefunction) -> Tuple[int, int, List[str]]:
    """(p, q, notes). ValueError when the three files cannot form a finite difference."""
    for label, w in (("N+1", plus), ("N-1", minus)):
        if w.n_atoms != n.n_atoms or not np.array_equal(w.atomic_numbers, n.atomic_numbers):
            raise ValueError(f"The {label} file holds different atoms from the N file; "
                             "all three files must describe the same molecule.")
        shift = float(np.abs(w.coordinates - n.coordinates).max())
        if shift > COORDINATE_TOLERANCE_BOHR:
            raise ValueError(
                f"The {label} geometry differs from the N geometry by up to "
                f"{shift * BOHR_RADIUS_ANGSTROM:.2e} A. All three states must use the optimized "
                "N-state geometry (the derivative is taken at fixed nuclei): run single points "
                "for the ions, do not re-optimize them.")
        if not same_basis(w, n):
            raise ValueError(f"The {label} file uses a different basis set from the N file.")
        if strip_reference_prefix(w.method) != strip_reference_prefix(n.method):
            raise ValueError(f"The {label} file was computed with {w.method}, the N file with "
                             f"{n.method}; all three states need the same method.")
    p = plus.n_electrons - n.n_electrons
    q = n.n_electrons - minus.n_electrons
    if p <= 0 or q <= 0:
        raise ValueError(
            f"Electron counts N = {n.n_electrons}, 'N+1' = {plus.n_electrons}, "
            f"'N-1' = {minus.n_electrons}: the N+1 file must hold more electrons than N and the "
            "N-1 file fewer. Were two files swapped?")
    notes: List[str] = []
    for label, w, k in (("N+1" if p == 1 else f"N+{p}", plus, p), ("N-1" if q == 1 else f"N-{q}", minus, q)):
        if abs(w.multiplicity - n.multiplicity) != k:
            notes.append(f"The {label} state has multiplicity {w.multiplicity} (N state: {n.multiplicity}). "
                         "Usually adding/removing electrons into/from a non-degenerate orbital changes "
                         f"the multiplicity by {k}; confirm this is the lowest-energy spin state.")
    for label, w in (("N", n), ("N+1", plus), ("N-1", minus)):
        if w.other_densities:
            notes.append(f"The {label} file also contains {', '.join(w.other_densities)}; "
                         f"the {gw.DENSITY_SECTION} is used here.")
    return p, q, notes


def frontier_degeneracy(w: gw.FchkWavefunction) -> Optional[Tuple[int, int]]:
    """(LUMO degeneracy p, HOMO degeneracy q) of a closed-shell state, Multiwfn's 0.01 eV rule."""
    e = w.alpha_orbital_energies
    if w.multiplicity != 1 or e is None or not 0 < w.n_alpha < len(e):
        return None
    e = e * HARTREE_TO_EV
    homo, lumo = e[w.n_alpha - 1], e[w.n_alpha]
    q = int(np.sum(np.abs(e[:w.n_alpha] - homo) < DEGENERACY_EV))
    p = int(np.sum(np.abs(e[w.n_alpha:] - lumo) < DEGENERACY_EV))
    return p, q


def has_diffuse_functions(basis_set: str) -> bool:
    b = basis_set.upper()
    return "+" in b or b.startswith("AUG-") or "D-AUG" in b


# ----------------------------------------------------------------------
# Formulas
# ----------------------------------------------------------------------
def fukui_from_densities(rho_n, rho_plus, rho_minus, p: int = 1, q: int = 1) -> Dict[str, np.ndarray]:
    f_plus = (np.asarray(rho_plus) - rho_n) / p
    f_minus = (np.asarray(rho_n) - rho_minus) / q
    return {"f+": f_plus, "f-": f_minus, "f0": 0.5 * (f_plus + f_minus), "DD": f_plus - f_minus}


@dataclass
class GlobalIndices:
    """Finite-difference global CDFT indices, Hartree, Multiwfn conventions (manual 3.25.1)."""
    e_n: float
    e_plus: float
    e_minus: float
    homo_n: Optional[float] = None

    @property
    def ip(self) -> float:
        return self.e_minus - self.e_n

    @property
    def ea(self) -> float:
        return self.e_n - self.e_plus

    @property
    def chi(self) -> float:
        return 0.5 * (self.ip + self.ea)

    @property
    def mu(self) -> float:
        return -self.chi

    @property
    def eta(self) -> float:
        return self.ip - self.ea  # no 1/2 (fundamental gap), as in Multiwfn

    @property
    def softness(self) -> float:
        return 1.0 / self.eta if self.eta > 0 else float("nan")

    @property
    def omega(self) -> float:
        return self.mu ** 2 / (2.0 * self.eta) if self.eta > 0 else float("nan")

    @property
    def nucleophilicity(self) -> Optional[float]:
        return None if self.homo_n is None else self.homo_n - TCE_HOMO_HARTREE


def homo_energy(w: gw.FchkWavefunction) -> Optional[float]:
    e = w.alpha_orbital_energies
    if w.multiplicity != 1 or e is None or not 0 < w.n_alpha <= len(e):
        return None
    return float(e[w.n_alpha - 1])


# ----------------------------------------------------------------------
# Condensed indices: Multiwfn CDFT.txt
# ----------------------------------------------------------------------
@dataclass
class CondensedRow:
    index: int        # 1-based
    symbol: str
    q_n: float
    q_plus: float
    q_minus: float
    f_minus: float
    f_plus: float
    f_zero: float
    cdd: float


CDFT_HEADER = ["Atom", "q(N)", "q(N+1)", "q(N-1)", "f-", "f+", "f0", "CDD"]
_CDFT_ROW = re.compile(r"^\s*(\d+)\(\s*([A-Za-z]{1,3})\s*\)\s+(.*)$")
_CDFT_GLOBAL = re.compile(r"^\s*([^:]+?)\s*:\s+(-?\d+\.\d+)\s+Hartree")
CDFT_GLOBAL_KEYS = {
    "E(N)": "e_n", "E(N+1)": "e_plus", "E(N-1)": "e_minus", "E_HOMO(N)": "homo_n",
    "First vertical IP": "ip", "Vertical IP": "ip",
    "First vertical EA": "ea", "Vertical EA": "ea",
    "Mulliken electronegativity": "chi", "Chemical potential": "mu",
    "Hardness (=fundamental gap)": "eta", "Softness": "softness",
    "Electrophilicity index": "omega", "Nucleophilicity index": "nucleophilicity",
}


def parse_cdft_text(text: str) -> Tuple[List[CondensedRow], Dict[str, float]]:
    """Condensed rows and global values (Hartree) from Multiwfn's CDFT.txt."""
    lines = text.splitlines()
    start = next((i for i, l in enumerate(lines) if l.split()[:8] == CDFT_HEADER), None)
    if start is None:
        raise ValueError("No 'Atom q(N) q(N+1) q(N-1) f- f+ f0 CDD' table found: this is not a "
                         "CDFT.txt written by Multiwfn main function 22, option 2.")
    rows: List[CondensedRow] = []
    for line in lines[start + 1:]:
        m = _CDFT_ROW.match(line)
        if not m:
            break
        try:
            values = [float(x) for x in m.group(3).split()]
        except ValueError:
            raise ValueError(f"CDFT.txt: unreadable condensed row '{line.strip()}'.") from None
        if len(values) != 7:
            raise ValueError(f"CDFT.txt: row '{line.strip()}' has {len(values)} numbers, expected 7.")
        rows.append(CondensedRow(int(m.group(1)), m.group(2), *values))
    if not rows:
        raise ValueError("CDFT.txt: the condensed table is empty.")
    found: Dict[str, float] = {}
    for line in lines:
        m = _CDFT_GLOBAL.match(line)
        if m and m.group(1) in CDFT_GLOBAL_KEYS:
            found[CDFT_GLOBAL_KEYS[m.group(1)]] = float(m.group(2))
    return rows, found


def condensed_inconsistencies(rows: Sequence[CondensedRow], tol: float = CONDENSED_TOLERANCE) -> List[int]:
    """1-based atoms whose printed f-, f+, f0, CDD disagree with their own charges."""
    bad = []
    for r in rows:
        checks = (r.q_minus - r.q_n - r.f_minus, r.q_n - r.q_plus - r.f_plus,
                  0.5 * (r.q_minus - r.q_plus) - r.f_zero, r.f_plus - r.f_minus - r.cdd)
        if max(abs(c) for c in checks) > tol:
            bad.append(r.index)
    return bad


def check_rows_match(rows: Sequence[CondensedRow], w: gw.FchkWavefunction) -> None:
    symbols = w.atom_symbols()
    if len(rows) != len(symbols):
        raise ValueError(f"CDFT.txt lists {len(rows)} atoms but the N file has {len(symbols)}.")
    for r, s in zip(rows, symbols):
        if r.symbol.lower() != s.lower():
            raise ValueError(f"CDFT.txt atom {r.index} is {r.symbol} but atom {r.index} of the N file is {s}.")


# ----------------------------------------------------------------------
# Condensed indices: running Multiwfn
# ----------------------------------------------------------------------
def find_multiwfn(explicit: Optional[str] = None) -> Optional[str]:
    candidates = []
    if explicit:
        candidates.append(explicit)
    env = os.environ.get("Multiwfnpath")
    if env:
        candidates.append(os.path.join(env, "Multiwfn.exe"))
    for name in ("Multiwfn.exe", "Multiwfn"):
        found = shutil.which(name)
        if found:
            candidates.append(found)
    for c in candidates:
        if c and os.path.isfile(c):
            return c
    return None


def multiwfn_cdft_script(settings_found: bool) -> str:
    """Keystrokes for main function 22, option 2, with N.fch given on the command line.

    N.wfn is absent, so Multiwfn asks for the N, N+1 and N-1 files in turn
    (verified with Multiwfn 2026.7.15). Without a settings.ini Multiwfn first
    pauses for ENTER, so only then does the script start with a blank line.
    """
    lead = "" if settings_found else "\n"
    return lead + "22\n2\nN.fch\nN+1.fch\nN-1.fch\n0\nq\n"


def run_multiwfn_cdft(multiwfn: str, paths: Sequence[str], timeout: float = MULTIWFN_TIMEOUT_S) -> str:
    """CDFT.txt text from Multiwfn; RuntimeError when Multiwfn cannot produce it."""
    settings = Path(multiwfn).resolve().parent / "settings.ini"
    with tempfile.TemporaryDirectory() as tmp:
        for name, src in zip(("N.fch", "N+1.fch", "N-1.fch"), paths):
            shutil.copyfile(src, Path(tmp) / name)
        if settings.is_file():  # Multiwfn reads settings.ini from the current folder first
            shutil.copyfile(settings, Path(tmp) / "settings.ini")
        try:
            proc = subprocess.run([multiwfn, "N.fch"], input=multiwfn_cdft_script(settings.is_file()),
                                  text=True, cwd=tmp, capture_output=True, timeout=timeout,
                                  creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        except (OSError, subprocess.SubprocessError) as exc:
            raise RuntimeError(f"Multiwfn could not be run: {exc}") from exc
        out = Path(tmp) / "CDFT.txt"
        if not out.is_file():
            tail = " | ".join((proc.stdout or "").strip().splitlines()[-3:])
            raise RuntimeError(f"Multiwfn finished without writing CDFT.txt (exit {proc.returncode}): {tail}")
        return out.read_text(encoding="utf-8", errors="replace")
