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

Everything is computed here, with no external program. Real-space functions
come from the exact densities (gaussianWavefunction, numba-accelerated when
installed). Atom-condensed indices come from Hirshfeld charges of the three
states (hirshfeldCharges: built-in free-atom densities, H-Ar), which reproduce
Multiwfn's main function 22 values to ~1e-4 e. Global indices follow
Multiwfn's conventions (hardness = I - A, without the 1/2).

These quantities indicate *likely* reactive regions; they are not
predictions of reaction outcome.
"""
from __future__ import annotations

import asyncio
import os
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
import hirshfeldCharges as hc
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

FUNCTIONS = ("f+", "f-", "f0", "DD")
FUNCTION_TITLES = {
    "f+": "f+  (nucleophilic attack)",
    "f-": "f-  (electrophilic attack)",
    "f0": "f0  (radical attack)",
    "DD": "dual descriptor  (+ nucleophilic / - electrophilic)",
}
CUBE_SUFFIXES = {"f+": "_fplus.cub", "f-": "_fminus.cub", "f0": "_fzero.cub", "DD": "_dual.cub"}
REPORT_SUFFIX = "_Fukui.txt"

CONDENSED_SOURCE = f"Hirshfeld (built in; sphericalized free atoms {hc.pad.LEVEL})"

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
# Condensed indices (Hirshfeld charges of the three states)
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


def condensed_rows(symbols: Sequence[str], q_n, q_plus, q_minus, p: int = 1, q: int = 1) -> List[CondensedRow]:
    """f+_A = [q_A(N) - q_A(N+p)]/p, f-_A = [q_A(N-q) - q_A(N)]/q, f0 = mean, CDD = f+ - f-."""
    rows = []
    for i, symbol in enumerate(symbols):
        f_plus = (q_n[i] - q_plus[i]) / p
        f_minus = (q_minus[i] - q_n[i]) / q
        rows.append(CondensedRow(i + 1, symbol, float(q_n[i]), float(q_plus[i]), float(q_minus[i]),
                                 float(f_minus), float(f_plus), float(0.5 * (f_plus + f_minus)),
                                 float(f_plus - f_minus)))
    return rows


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
            condensed: bool = True,
            progress: Optional[ProgressCallback] = None) -> FukuiResults:
        started = time.perf_counter()

        def report(fraction: float, message: str):
            if progress is not None:
                try:
                    progress(fraction, message)
                except Exception:
                    pass  # reporting must never abort an analysis

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

        rows, source = [], ""
        if condensed:
            report(0.85, "Condensed indices: Hirshfeld charges of the three states...")
            rows, source = self._condensed(states, p, q, notes)
        report(1.0, "Done.")
        return FukuiResults(paths, states, p, q, fukui_isovalue, surface_isovalue, spacing, origin,
                            grids, surface, surface_values, extrema, lobes, globals_, rows, source,
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
    def _condensed(states, p, q, notes) -> Tuple[List[CondensedRow], str]:
        """Hirshfeld charges of N, N+p, N-q on one shared grid; a note (never an abort) on failure."""
        n = states[0]
        try:
            grid = hc.molecular_grid(n.coordinates, n.atom_symbols())
            charges = [hc.hirshfeld_charges(w, grid) for w in states]
        except ValueError as exc:
            notes.append(f"Condensed indices were not computed: {exc}")
            return [], ""
        rows = condensed_rows(n.atom_symbols(), *charges, p=p, q=q)
        negative = [r.index for r in rows if r.f_plus < 0 or r.f_minus < 0]
        if negative:
            notes.append(f"Negative condensed f+ or f- on atoms {negative}. Hirshfeld values are usually "
                         "non-negative; negative values can come from orbital relaxation or the partition "
                         "and should not be read as reactivity.")
        return rows, CONDENSED_SOURCE

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


# ----------------------------------------------------------------------
# Toga UI
# ----------------------------------------------------------------------
DISPLAY_LOBES = "Isosurfaces (+/- isovalue)"
DISPLAY_MAPPED = "Mapped on rho(N) surface"
DISPLAY_MODES = (DISPLAY_LOBES, DISPLAY_MAPPED)
# Atom drawing styles of the 3D view: label -> molecularViewer.atom_display_style.
# "CPK" is ball-and-stick here, as in the 3D Molecular Viewer; space-filling is "vdW".
ATOM_STYLES = {"Line": "Line style", "Ball and stick (CPK)": "CPK style", "Space-filling (vdW)": "vdW style"}
DEFAULT_ATOM_STYLE = "Line"
BACKGROUND_COLOURS = {
    "Black": (0.0, 0.0, 0.0),
    "White": (1.0, 1.0, 1.0),
    "Dark grey": (0.25, 0.25, 0.25),
    "Light grey": (0.85, 0.85, 0.85),
}
IMAGE_FILE_TYPES = ["png", "tif", "tiff", "bmp"]


class FukuiFunctionAnalyzeUI:
    LABEL_WIDTH = 230
    HINT_COLOR = "#666666"

    def __init__(self, *args, **kwargs):
        # gqteaWinToga passes the Button that was pressed; ignore it.
        self.analyser = FukuiAnalyser()
        self.results: Optional[FukuiResults] = None
        self.written: List[str] = []
        self.viewer: Optional[ESPSurfaceViewer] = None
        self._syncing = False
        self.window = toga.Window(title="Fukui Function Analysis", size=(900, 700), on_close=self._on_close)
        self.layout_main_window()
        self.window.show()

    # --- layout -----------------------------------------------------------
    def _label(self, text):
        return toga.Label(text, style=Pack(width=self.LABEL_WIDTH, margin=(6, 6, 0, 0)))

    def _row(self, *children):
        return toga.Box(children=list(children), style=Pack(direction=ROW, margin_bottom=4))

    def _heading(self, text):
        return toga.Label(text, style=Pack(font_weight="bold", margin=(10, 0, 4, 0)))

    def _hint(self, text):
        return toga.Label(text, style=Pack(color=self.HINT_COLOR, font_size=9, margin=(0, 0, 4, 0)))

    def _file_row(self, label, attr, handler):
        field_ = toga.TextInput(placeholder=".fchk / .fch", style=Pack(flex=1))
        setattr(self, attr, field_)
        return self._row(self._label(label), field_,
                         toga.Button("Browse", on_press=handler, style=Pack(width=90, margin_left=5)))

    def layout_main_window(self):
        self.output_dir_input = toga.TextInput(placeholder="Blank = folder of the N file", style=Pack(flex=1))
        header = toga.Box(children=[
            self._file_row("N state (optimized geometry):", "n_input", self.browse_n),
            self._file_row("N+1 state (anion of a neutral):", "plus_input", self.browse_plus),
            self._file_row("N-1 state (cation of a neutral):", "minus_input", self.browse_minus),
            self._row(self._label("Output folder:"), self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
        ], style=Pack(direction=COLUMN))

        self.fukui_iso_input = toga.TextInput(placeholder=f"{DEFAULT_FUKUI_ISOVALUE:g}", style=Pack(width=160))
        self.surface_iso_input = toga.TextInput(placeholder=f"{DEFAULT_SURFACE_ISOVALUE:g}", style=Pack(width=160))
        self.spacing_input = toga.TextInput(placeholder=f"{DEFAULT_SPACING_BOHR:g}", style=Pack(width=160))
        self.condensed_switch = toga.Switch("Compute atom-condensed indices (Hirshfeld charges)", value=True)
        self.cubes_switch = toga.Switch("Write cube files (f+, f-, f0, dual descriptor)", value=True)
        setup = toga.Box(children=[
            self._heading("Real-space functions"),
            self._row(self._label("Fukui isovalue (a.u.):"), self.fukui_iso_input),
            self._hint("Blank = 0.005. Isosurfaces are drawn at +value and -value."),
            self._row(self._label("rho(N) surface for mapping (a.u.):"), self.surface_iso_input),
            self._hint("Blank = 0.01, as in Multiwfn manual 4.12.4 (f is too small on rho = 0.001)."),
            self._row(self._label("Grid spacing (bohr):"), self.spacing_input),
            self._hint("Blank = 0.25 bohr. Smaller = smoother and slower."),
            self._heading("Atom-condensed indices"),
            self.condensed_switch,
            self._hint(f"Built-in Hirshfeld partition, free atoms {hc.pad.LEVEL} (elements H-Ar). "
                       "Other elements: the run continues without them, with a note."),
            self.cubes_switch,
        ], style=Pack(direction=COLUMN, margin=10))

        self.function_selection = toga.Selection(items=list(FUNCTIONS), on_change=self.on_display_change,
                                                 style=Pack(width=140))
        self.display_selection = toga.Selection(items=list(DISPLAY_MODES), on_change=self.on_display_change,
                                                style=Pack(width=240))
        self.extrema_switch = toga.Switch("Show extrema markers", value=True, on_change=self.on_display_change)
        self.atom_style_selection = toga.Selection(items=list(ATOM_STYLES), on_change=self.on_atom_view_change,
                                                   style=Pack(width=200))
        self.atom_numbers_switch = toga.Switch("Show atom numbers", value=False,
                                               on_change=self.on_atom_view_change)
        self.background_selection = toga.Selection(items=list(BACKGROUND_COLOURS),
                                                   on_change=self.on_background_change, style=Pack(width=140))
        self.export_button = toga.Button("Export image...", on_press=self.export_image,
                                         style=Pack(width=130, margin_left=10))
        self.summary_output = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_fukui,
                                                      style=Pack(flex=1, font_family="monospace"))
        results = toga.Box(children=[
            self._row(self._label("Function:"), self.function_selection),
            self._row(self._label("Display:"), self.display_selection, self.extrema_switch),
            self._row(self._label("Atoms:"), self.atom_style_selection, self.atom_numbers_switch),
            self._row(self._label("Canvas background:"), self.background_selection, self.export_button),
            self.summary_output,
        ], style=Pack(direction=COLUMN, margin=10, flex=1))

        self.tabs = toga.OptionContainer(content=[
            ("Setup", toga.ScrollContainer(content=setup, style=Pack(flex=1))),
            ("Results", results),
        ], style=Pack(flex=1))
        self.run_button = toga.Button("Run", on_press=self.workflow, style=Pack(width=110, margin=5))
        self.show_button = toga.Button("Show 3D", on_press=self.show_surface, style=Pack(width=110, margin=5))
        self.clear_button = toga.Button("Clear", on_press=self.clear_output, style=Pack(width=110, margin=5))
        self.close_button = toga.Button("Close", on_press=self.closeTopLevel, style=Pack(width=110, margin=5))
        self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
        self.progress_bar = toga.ProgressBar(max=100, style=Pack(margin=(0, 5, 5, 5)))
        self.window.content = toga.Box(children=[
            header,
            self._row(self.run_button, self.show_button, self.clear_button, self.close_button),
            self.tabs, self.status_label, self.progress_bar,
        ], style=Pack(direction=COLUMN, margin=10, flex=1))

    # --- dialogs ------------------------------------------------------------
    async def _show_error(self, title: str, message: str):
        await self.window.dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        await self.window.dialog(toga.InfoDialog(title, message))

    async def _browse_into(self, target, title, file_types):
        try:
            path = await self.window.dialog(toga.OpenFileDialog(title=title, file_types=file_types))
            if path:
                target.value = str(path)
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))

    async def browse_n(self, widget):
        await self._browse_into(self.n_input, "Formatted checkpoint of the N state", ["fchk", "fch"])

    async def browse_plus(self, widget):
        await self._browse_into(self.plus_input, "Formatted checkpoint of the N+1 state", ["fchk", "fch"])

    async def browse_minus(self, widget):
        await self._browse_into(self.minus_input, "Formatted checkpoint of the N-1 state", ["fchk", "fch"])

    async def browse_output_dir(self, widget):
        try:
            folder = await self.window.dialog(toga.SelectFolderDialog(title="Choose the output folder"))
            if folder:
                self.output_dir_input.value = str(folder)
        except Exception as exc:
            await self._show_error("Folder Selection Error", str(exc))

    # --- run ----------------------------------------------------------------
    def _collect_inputs(self) -> Dict[str, object]:
        paths = {}
        for key, widget, label in (("n_path", self.n_input, "N"), ("plus_path", self.plus_input, "N+1"),
                                   ("minus_path", self.minus_input, "N-1")):
            path = (widget.value or "").strip()
            if not path:
                raise ValueError(f"Please select the {label} state .fchk file.")
            if not os.path.isfile(path):
                raise ValueError(f"File not found: {path}")
            paths[key] = path
        return {
            **paths,
            "fukui_isovalue": parse_positive_float(self.fukui_iso_input.value, DEFAULT_FUKUI_ISOVALUE,
                                                   "Fukui isovalue", 1.0),
            "surface_isovalue": parse_positive_float(self.surface_iso_input.value, DEFAULT_SURFACE_ISOVALUE,
                                                     "rho surface isovalue", 1.0),
            "spacing": parse_positive_float(self.spacing_input.value, DEFAULT_SPACING_BOHR,
                                            "Grid spacing", 0.6),
            "condensed": bool(self.condensed_switch.value),
        }

    def _apply_progress(self, fraction: float, message: str):
        self.progress_bar.value = max(0.0, min(100.0, 100.0 * fraction))
        self.status_label.text = message

    def _set_buttons(self, enabled: bool):
        for b in (self.run_button, self.show_button, self.clear_button):
            b.enabled = enabled

    async def workflow(self, widget):
        try:
            params = self._collect_inputs()
            output_dir = resolve_output_dir(params["n_path"], self.output_dir_input.value)
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        self._set_buttons(False)
        self.progress_bar.value = 0
        loop = asyncio.get_running_loop()

        def on_progress(fraction, message):
            loop.call_soon_threadsafe(self._apply_progress, fraction, message)

        try:
            results = await asyncio.to_thread(self.analyser.run, progress=on_progress, **params)
            self.written = self.analyser.write_outputs(results, output_dir, bool(self.cubes_switch.value))
            self.results = results
            self.summary_output.value = self.analyser.report_text(results, self.written)
            self.tabs.current_tab = 1
            # File name only: a long path would widen the window, and Toga never shrinks it.
            self.status_label.text = (f"Done in {results.elapsed:.1f} s. Report: {Path(self.written[0]).name} "
                                      "(full paths at the end of the report)")
            self.progress_bar.value = 100
            if self.viewer is not None and self.viewer.viewer_window_is_open():
                self._push_to_viewer()
        except Exception as exc:
            self.status_label.text = "Calculation failed."
            self.progress_bar.value = 0
            await self._show_error("Fukui Calculation Error", str(exc))
        finally:
            self._set_buttons(True)

    # --- display --------------------------------------------------------
    def selected_function(self) -> str:
        return self.function_selection.value if self.function_selection.value in FUNCTIONS else "f-"

    def surface_markers(self, key: str) -> List[SurfaceMarker]:
        return [SurfaceMarker(r["position_angstrom"], r["kind"], r["scope"], str(r["row"]))
                for r in self.analyser.extrema_rows(self.results, key)]

    def _push_to_viewer(self):
        if self.results is None or self.viewer is None:
            return
        key = self.selected_function()
        n = self.results.n_state
        self.viewer.set_molecule(n.atom_symbols(), n.coordinates * BOHR_RADIUS_ANGSTROM)
        self._apply_atom_view()
        self.viewer.set_display_options(unit=UNIT_AU, show_extrema=bool(self.extrema_switch.value))
        if self.display_selection.value == DISPLAY_MAPPED:
            s = self.results.surface
            self.viewer.set_surface(
                s.vertices * BOHR_RADIUS_ANGSTROM, s.triangles, s.vertex_normals(),
                self.results.surface_values[key], self.results.surface_isovalue, self.surface_markers(key),
                caption=f"{key} (a.u.)   mapped on rho(N) = {self.results.surface_isovalue:g} a.u.")
            return
        lobe = self.results.lobes[key]
        if lobe is None:
            self.viewer.clear_surface()
            self.status_label.text = f"{key} has no isosurface at +/-{self.results.fukui_isovalue:g} a.u."
            return
        iso = self.results.fukui_isovalue
        self.viewer.set_surface(lobe.vertices * BOHR_RADIUS_ANGSTROM, lobe.triangles, lobe.normals,
                                lobe.values, iso, (), half_range=iso,
                                caption=f"{key} isosurfaces at +/-{iso:g} a.u.: blue positive, red negative")

    def selected_atom_style(self) -> str:
        return ATOM_STYLES.get(self.atom_style_selection.value, ATOM_STYLES[DEFAULT_ATOM_STYLE])

    def _apply_atom_view(self):
        self.viewer.set_atom_view(self.selected_atom_style(), bool(self.atom_numbers_switch.value))

    def on_atom_view_change(self, widget):
        # A view setting, not a result: applies with or without a run.
        if self.viewer is not None:
            self._apply_atom_view()

    def on_display_change(self, widget):
        if self._syncing or self.results is None or self.viewer is None:
            return
        self._push_to_viewer()

    async def show_surface(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Run the calculation first.")
            return
        if self.viewer is None:
            self.viewer = ESPSurfaceViewer()
            self.viewer.window_title = "Fukui Function Viewer"
            self.viewer.loop = asyncio.get_running_loop()
            self.viewer.set_background_color(self.selected_background())
        self._push_to_viewer()
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)

    def selected_background(self):
        return BACKGROUND_COLOURS.get(self.background_selection.value, BACKGROUND_COLOURS["Black"])

    def on_background_change(self, widget):
        if self.viewer is not None:
            self.viewer.set_background_color(self.selected_background())

    def default_image_name(self) -> str:
        if self.results is None:
            return "Fukui.png"
        tag = {"f+": "fplus", "f-": "fminus", "f0": "fzero", "DD": "dual"}[self.selected_function()]
        return f"{Path(self.results.paths[0]).stem}_{tag}.png"

    async def _ask_image_path(self, suggested: str):
        return await self.window.dialog(toga.SaveFileDialog(
            title="Export the Fukui function image", suggested_filename=suggested, file_types=IMAGE_FILE_TYPES))

    async def export_image(self, widget):
        """Capture first, then ask for the file name (same contract as the ESP tool)."""
        if self.viewer is None or not self.viewer.viewer_window_is_open():
            await self._show_info("No Image", "Open the 3D view with 'Show 3D' first.")
            return
        try:
            image = await self.viewer.capture_frame()
        except RuntimeError as exc:
            await self._show_error("Export Image Error", str(exc))
            return
        try:
            path = await self._ask_image_path(self.default_image_name())
            if not path:
                return
            saved = save_canvas_image(image, str(path))
        except (ValueError, OSError) as exc:
            await self._show_error("Export Image Error", str(exc))
            return
        self.status_label.text = f"Image saved: {saved}"

    def clear_output(self, widget):
        self.results = None
        self.written = []
        self.summary_output.value = HelpGqteaWin.help_fukui
        if self.viewer is not None:
            self.viewer.clear_surface()
        self.progress_bar.value = 0
        self.status_label.text = "Ready."

    def _on_close(self, window, **kwargs) -> bool:
        if self.viewer is not None:
            self.viewer.request_close()
        return True

    def closeTopLevel(self, widget):
        self._on_close(self.window)
        self.window.close()
