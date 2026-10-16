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


# ----------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------
@dataclass
class LobeSurface:
    vertices: np.ndarray   # (V, 3) bohr
    triangles: np.ndarray  # (F, 3)
    normals: np.ndarray    # (V, 3)
    values: np.ndarray     # (V,) +isovalue or -isovalue


@dataclass
class FukuiResults:
    paths: Tuple[str, str, str]
    states: Tuple[gw.FchkWavefunction, gw.FchkWavefunction, gw.FchkWavefunction]
    p: int
    q: int
    fukui_isovalue: float
    surface_isovalue: float
    spacing: float
    origin: np.ndarray
    grids: Dict[str, np.ndarray]                 # a.u., on the cube grid
    surface: ms.SurfaceMesh                      # rho(N) = surface_isovalue, bohr
    surface_values: Dict[str, np.ndarray]        # a.u., one per vertex
    extrema: Dict[str, List[ms.SurfaceExtremum]]
    lobes: Dict[str, Optional[LobeSurface]]
    global_indices: Optional[GlobalIndices]
    condensed: List[CondensedRow]
    condensed_source: str
    notes: List[str] = field(default_factory=list)
    elapsed: float = 0.0

    @property
    def n_state(self) -> gw.FchkWavefunction:
        return self.states[0]

    def integral(self, key: str) -> float:
        return float(self.grids[key].sum() * self.spacing ** 3)


class FukuiAnalyser:
    """Headless pipeline: three fchk -> validation -> grids -> surface -> condensed indices."""

    def run(self, n_path: str, plus_path: str, minus_path: str,
            fukui_isovalue: float = DEFAULT_FUKUI_ISOVALUE,
            surface_isovalue: float = DEFAULT_SURFACE_ISOVALUE,
            spacing: float = DEFAULT_SPACING_BOHR,
            condensed_source: str = CONDENSED_MULTIWFN,
            multiwfn_path: Optional[str] = None, cdft_path: Optional[str] = None,
            progress: Optional[ProgressCallback] = None) -> FukuiResults:
        started = time.perf_counter()

        def report(fraction: float, message: str):
            if progress is not None:
                try:
                    progress(fraction, message)
                except Exception:
                    pass  # reporting must never abort an analysis

        if condensed_source not in CONDENSED_SOURCES:
            raise ValueError(f"Unknown condensed-index source '{condensed_source}'.")
        paths = (str(n_path), str(plus_path), str(minus_path))
        report(0.02, "Reading the three formatted checkpoints...")
        states = tuple(gw.read_fchk(p) for p in paths)
        n, plus, minus = states
        p, q, notes = validate_states(n, plus, minus)

        globals_ = None
        if p == 1 and q == 1:
            if None in (n.total_energy, plus.total_energy, minus.total_energy):
                notes.append("A 'Total Energy' entry is missing; global indices were not computed.")
            else:
                globals_ = GlobalIndices(n.total_energy, plus.total_energy, minus.total_energy, homo_energy(n))
                if not globals_.eta > 0:
                    notes.append("I - A is not positive; check the charges and multiplicities of the states.")
        else:
            notes.append(f"p = {p}, q = {q}: degenerate-frontier-orbital Fukui functions "
                         "(Multiwfn manual 3.25.4). Global indices need the N+1 and N-1 energies "
                         "and were not computed.")
        degeneracy = frontier_degeneracy(n)
        if degeneracy is not None and (degeneracy[0], degeneracy[1]) != (p, q) and max(degeneracy) > 1:
            notes.append(f"The N state has LUMO degeneracy {degeneracy[0]} and HOMO degeneracy "
                         f"{degeneracy[1]} (orbitals within {DEGENERACY_EV} eV) but the files give "
                         f"p = {p}, q = {q}. Standard f+/f- can then break the molecular symmetry; "
                         "use N+p / N-q states (manual 3.25.4) or orbital-weighted Fukui functions "
                         "(Multiwfn main function 22, options 4-7).")
        if has_diffuse_functions(n.basis_set):
            notes.append(f"The basis set {n.basis_set} has diffuse functions. They improve the electron "
                         "affinity, but with self-interaction-prone functionals the N+1 density can become "
                         "over-diffuse and f+ / the dual descriptor lose chemical meaning (manual 4.22.3.1).")

        def densities(points):
            return [gw.electron_density(w, points) for w in states]

        def fukui_at(points):
            return fukui_from_densities(*densities(points), p, q)

        def envelope(points):
            rn, rp, rm = densities(points)
            f = fukui_from_densities(rn, rp, rm, p, q)
            return np.maximum.reduce([rn / surface_isovalue,
                                      np.abs(f["f+"]) / fukui_isovalue,
                                      np.abs(f["f-"]) / fukui_isovalue])

        report(0.08, "Sizing the grid...")
        radii = np.array([AtomicData.vdw_radii.get(s, AtomicData.DEFAULT_VDW_RADIUS)
                          for s in n.atom_symbols()]) / BOHR_RADIUS_ANGSTROM
        origin, shape = ms.grid_box(n.coordinates, radii, envelope, 1.0, spacing)
        points = ms.grid_points(origin, spacing, shape)
        rho = []
        for k, (label, w) in enumerate(zip(("N", "N+1", "N-1"), states)):
            report(0.15 + 0.15 * k, f"Electron density of the {label} state on {len(points)} points...")
            rho.append(gw.electron_density(w, points).reshape(shape))
        grids = fukui_from_densities(*rho, p, q)

        report(0.62, f"rho(N) = {surface_isovalue:g} a.u. surface and the functions on it...")
        density_n = lambda pts: gw.electron_density(n, pts)  # noqa: E731
        surface = ms.build_isosurface(rho[0], origin, spacing, surface_isovalue,
                                      density_n, ms.density_outward(density_n))
        surface_values = fukui_at(surface.vertices)
        extrema = {k: ms.find_surface_extrema(v, surface.vertices, surface.neighbours)
                   for k, v in surface_values.items()}

        report(0.75, f"Isosurfaces at +/-{fukui_isovalue:g} a.u....")
        lobes = {k: self.lobe_surface(grids[k], origin, spacing, fukui_isovalue,
                                      lambda pts, key=k: fukui_at(pts)[key]) for k in FUNCTIONS}
        for k, lobe in lobes.items():
            if lobe is None:
                notes.append(f"{k} never reaches +/-{fukui_isovalue:g} a.u.; lower the Fukui isovalue to see it.")

        report(0.85, "Condensed indices...")
        condensed, source = self._condensed(condensed_source, paths, n, p, q, multiwfn_path, cdft_path, notes)
        report(1.0, "Done.")
        return FukuiResults(paths, states, p, q, fukui_isovalue, surface_isovalue, spacing, origin,
                            grids, surface, surface_values, extrema, lobes, globals_, condensed, source,
                            notes, time.perf_counter() - started)

    @staticmethod
    def lobe_surface(grid, origin, spacing, isovalue, exact) -> Optional[LobeSurface]:
        """The +isovalue and -isovalue surfaces of one function, as one mesh."""
        verts, tris, norms, vals = [], [], [], []
        offset = 0
        for sign in (1.0, -1.0):
            if (sign * grid).max() < isovalue:
                continue
            signed = lambda pts, s=sign: s * exact(pts)  # noqa: E731
            mesh = ms.build_isosurface(sign * grid, origin, spacing, isovalue, signed, ms.density_outward(signed))
            verts.append(mesh.vertices)
            tris.append(mesh.triangles + offset)
            norms.append(mesh.vertex_normals())
            vals.append(np.full(len(mesh.vertices), sign * isovalue))
            offset += len(mesh.vertices)
        if not verts:
            return None
        return LobeSurface(np.vstack(verts), np.vstack(tris), np.vstack(norms), np.concatenate(vals))

    @staticmethod
    def _condensed(source, paths, n, p, q, multiwfn_path, cdft_path, notes):
        if source == CONDENSED_NONE:
            return [], source
        if (p, q) != (1, 1):
            notes.append("Condensed indices for N+p / N-q states are not supported here; use Multiwfn "
                         "main function 22 option -3 then option 2 (manual 3.25.4).")
            return [], CONDENSED_NONE
        try:
            if source == CONDENSED_IMPORT:
                if not cdft_path:
                    raise RuntimeError("no CDFT.txt file was given")
                text = Path(cdft_path).read_text(encoding="utf-8", errors="replace")
            else:
                exe = find_multiwfn(multiwfn_path)
                if exe is None:
                    raise RuntimeError("Multiwfn.exe was not found (set its path, or Multiwfnpath, or PATH)")
                text = run_multiwfn_cdft(exe, paths)
        except (OSError, RuntimeError) as exc:
            notes.append(f"Condensed indices were not computed: {exc}.")
            return [], CONDENSED_NONE
        rows, _ = parse_cdft_text(text)
        check_rows_match(rows, n)
        bad = condensed_inconsistencies(rows)
        if bad:
            notes.append(f"CDFT.txt rows for atoms {bad} do not follow f- = q(N-1) - q(N), "
                         "f+ = q(N) - q(N+1); check the file.")
        negative = [r.index for r in rows if r.f_plus < 0 or r.f_minus < 0]
        if negative:
            notes.append(f"Negative condensed f+ or f- on atoms {negative}. Hirshfeld values are usually "
                         "non-negative; negative values can come from orbital relaxation or the partition "
                         "and should not be read as reactivity.")
        return rows, source

    # --- reporting ------------------------------------------------------
    @staticmethod
    def nearest_atom_label(w: gw.FchkWavefunction, position_bohr: np.ndarray) -> str:
        i = int(np.argmin(np.linalg.norm(w.coordinates - position_bohr, axis=1)))
        return f"{w.atom_symbols()[i]}{i + 1}"

    def extrema_rows(self, results: FukuiResults, key: str) -> List[Dict[str, object]]:
        rows = []
        for e in results.extrema[key]:
            rows.append({"row": len(rows) + 1, "kind": e.kind, "scope": e.scope, "value": e.value,
                         "position_angstrom": e.position * BOHR_RADIUS_ANGSTROM,
                         "nearest_atom": self.nearest_atom_label(results.n_state, e.position)})
        return rows

    @staticmethod
    def site_ranking(rows: Sequence[CondensedRow], top: int = 5) -> List[Tuple[str, List[CondensedRow]]]:
        def best(key, reverse):
            return sorted(rows, key=key, reverse=reverse)[:top]
        return [
            ("electrophilic attack: largest f-", best(lambda r: r.f_minus, True)),
            ("nucleophilic attack: largest f+", best(lambda r: r.f_plus, True)),
            ("radical attack: largest f0", best(lambda r: r.f_zero, True)),
            ("electrophilic attack: most negative CDD", best(lambda r: r.cdd, False)),
            ("nucleophilic attack: most positive CDD", best(lambda r: r.cdd, True)),
        ]

    def report_text(self, results: FukuiResults, written: Sequence[str] = ()) -> str:
        n, plus, minus = results.states
        p, q = results.p, results.q
        lines = [
            "# gQTEA Fukui Function Analysis",
            "# Indicators of likely reactive regions -- not predictions of reaction outcome.",
            f"# file_N                 {Path(results.paths[0]).name}",
            f"# file_N+{p}               {Path(results.paths[1]).name}",
            f"# file_N-{q}               {Path(results.paths[2]).name}",
            f"# method / basis         {n.method} / {n.basis_set}",
            f"# electrons              {n.n_electrons} / {plus.n_electrons} / {minus.n_electrons}",
            f"# charge multiplicity    N: {n.charge} {n.multiplicity}   N+{p}: {plus.charge} {plus.multiplicity}"
            f"   N-{q}: {minus.charge} {minus.multiplicity}",
            f"# density                {gw.DENSITY_SECTION}",
            f"# definitions            f+ = [rho(N+{p}) - rho(N)]/{p}; f- = [rho(N) - rho(N-{q})]/{q}; "
            "f0 = (f+ + f-)/2; DD = f+ - f-",
            f"# grid_spacing_bohr      {results.spacing:g}",
            f"# fukui_isovalue_au      {results.fukui_isovalue:g}",
            f"# surface                rho(N) = {results.surface_isovalue:g} a.u., "
            f"{len(results.surface.vertices)} vertices",
            "",
        ]
        g = results.global_indices
        if g is not None:
            def both(value):
                return f"{value:12.6f} Hartree  {value * HARTREE_TO_EV:10.4f} eV"
            lines += [
                "GLOBAL INDICES (hardness = I - A, without 1/2, as in Multiwfn)",
                f"  E(N)                 {g.e_n:.6f} Hartree",
                f"  E(N+1)               {g.e_plus:.6f} Hartree",
                f"  E(N-1)               {g.e_minus:.6f} Hartree",
                f"  Vertical IP          {both(g.ip)}",
                f"  Vertical EA          {both(g.ea)}",
                f"  Electronegativity    {both(g.chi)}",
                f"  Chemical potential   {both(g.mu)}",
                f"  Hardness             {both(g.eta)}",
                f"  Softness             {g.softness:12.6f} Hartree^-1",
                f"  Electrophilicity     {both(g.omega)}",
            ]
            if g.nucleophilicity is not None:
                lines.append(f"  Nucleophilicity      {both(g.nucleophilicity)}   "
                             "(E_HOMO(TCE) at B3LYP/6-31G*, as in Multiwfn)")
            lines.append("")
        lines.append("GRID INTEGRALS (ideally 1 for f+, f-, f0 and 0 for DD; a finite uniform grid undercounts)")
        lines.append("  " + "   ".join(f"{k} {results.integral(k):.4f}" for k in FUNCTIONS))
        lines.append("")
        if results.condensed:
            s = g.softness if g is not None else float("nan")
            lines += [f"CONDENSED INDICES  (source: {results.condensed_source}; Hirshfeld charges, e; "
                      "s = S x f in e/Hartree)",
                      "  atom      q(N)   q(N+1)   q(N-1)       f-       f+       f0      CDD       s-       s+       s0"]
            for r in results.condensed:
                lines.append(f"  {r.symbol:>2s}{r.index:<4d} {r.q_n:8.4f} {r.q_plus:8.4f} {r.q_minus:8.4f} "
                             f"{r.f_minus:8.4f} {r.f_plus:8.4f} {r.f_zero:8.4f} {r.cdd:8.4f} "
                             f"{s * r.f_minus:8.4f} {s * r.f_plus:8.4f} {s * r.f_zero:8.4f}")
            lines += ["", "SITE RANKING (consider only atoms that can actually react)"]
            for title, best in self.site_ranking(results.condensed):
                lines.append(f"  {title:42s} " + ", ".join(f"{r.symbol}{r.index}" for r in best))
            lines.append("")
        for key in FUNCTIONS:
            lines += [f"SURFACE EXTREMA OF {FUNCTION_TITLES[key]}  on rho(N) = {results.surface_isovalue:g} a.u.",
                      "   row  type     scope        value (a.u.)    x (A)     y (A)     z (A)   atom"]
            for r in self.extrema_rows(results, key):
                x, y, z = r["position_angstrom"]
                lines.append(f"  {r['row']:4d}  {r['kind']:7s}  {r['scope']:6s}  {r['value']:16.8f}  "
                             f"{x:8.4f}  {y:8.4f}  {z:8.4f}   {r['nearest_atom']}")
            lines.append("")
        for note in results.notes:
            lines.append(f"NOTE: {note}")
        for path in written:
            lines.append(f"Written: {path}")
        return "\n".join(lines) + "\n"

    def write_outputs(self, results: FukuiResults, output_dir: str, write_cubes: bool = True) -> List[str]:
        stem = Path(results.paths[0]).stem
        out = Path(output_dir)
        written: List[str] = []
        n = results.n_state
        if write_cubes:
            for key in FUNCTIONS:
                cube = cf.CubeData(f"{stem} {FUNCTION_TITLES[key]}", "gQTEA Fukui Function Analysis, a.u.",
                                   n.atomic_numbers, n.nuclear_charges, n.coordinates,
                                   results.origin, results.spacing, results.grids[key])
                written.append(cf.write_cube(out / (stem + CUBE_SUFFIXES[key]), cube))
        report = out / (stem + REPORT_SUFFIX)
        try:
            report.parent.mkdir(parents=True, exist_ok=True)
            report.write_text(self.report_text(results, written), encoding="utf-8")
        except OSError as exc:
            raise ValueError(f"Could not write {report}: {exc}") from exc
        return [str(report)] + written
