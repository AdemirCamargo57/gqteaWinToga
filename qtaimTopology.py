"""QTAIM topology of the electron density: critical points and gradient paths.

Pure numpy/scipy, no Toga; every length is in bohr. The algorithms and their
defaults are Multiwfn's (topology analysis, main function 2, manual 3.14),
read from Multiwfn 2026.7.15 itself:

* critical points (CPs) by Newton iterations from starting points at the
  nuclei (mode 2), atom-pair midpoints (3), triangle centres (4) and
  tetrahedron centres (5); combinations whose atoms are farther apart than
  ``vdw_cutoff_scale`` x the sum of their vdW radii are skipped;
* a CP is new only if it is at least ``min_cp_distance`` from every CP found
  before (first found wins, modes in ascending order);
* its type (rank, signature) comes from the signs of the Hessian eigenvalues,
  and the Poincare-Hopf sum n(3,-3) - n(3,-1) + n(3,+1) - n(3,+3) must be 1
  for a molecule;
* bond paths climb the gradient from each (3,-1) displaced along the
  eigenvector of its positive eigenvalue; ring paths descend from each (3,+1)
  along the eigenvector of its negative eigenvalue to a (3,+3).

The density and its derivatives come from gaussianWavefunction.density_derivatives
(the compiled densityNumba kernel when numba is installed).
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from itertools import combinations
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

import gaussianWavefunction as gw
from help import AtomicData

BOHR_RADIUS_ANGSTROM = 0.529177210544  # CODATA 2022; equal to coordinateConverter's (tested)

NCP, BCP, RCP, CCP = "(3,-3)", "(3,-1)", "(3,+1)", "(3,+3)"
CP_TYPES = (NCP, BCP, RCP, CCP)
CP_NAMES = {NCP: "nuclear", BCP: "bond", RCP: "ring", CCP: "cage"}

MODE_NUCLEI, MODE_PAIRS, MODE_TRIANGLES, MODE_TETRAHEDRA = 2, 3, 4, 5
MODE_LABELS = {MODE_NUCLEI: "nuclei", MODE_PAIRS: "atom pairs",
               MODE_TRIANGLES: "triangles", MODE_TETRAHEDRA: "tetrahedra"}
MODE_SIZES = {MODE_NUCLEI: 1, MODE_PAIRS: 2, MODE_TRIANGLES: 3, MODE_TETRAHEDRA: 4}

MIN_CP_DENSITY = 1e-8      # a.u.; Newton points below it are far-field noise, not CPs
ZERO_EIGENVALUE = 1e-10    # |lambda| below it = degenerate CP (rank < 3)
NNA_DISTANCE = 0.3         # bohr; a (3,-3) farther from every nucleus is a non-nuclear attractor

ProgressCallback = Callable[[float], None]


def _report(progress: Optional[ProgressCallback], fraction: float):
    if progress is not None:
        try:
            progress(min(1.0, max(0.0, fraction)))
        except Exception:
            pass  # reporting must never abort an analysis


# ----------------------------------------------------------------------
# Parameters (Multiwfn topology menus -1 and -2)
# ----------------------------------------------------------------------
@dataclass
class CPSearchParameters:
    max_iterations: int = 120
    step_scale: float = 1.0
    gradient_tolerance: float = 1e-6
    displacement_tolerance: float = 1e-7
    min_cp_distance: float = 0.03          # bohr
    vdw_cutoff_scale: float = 1.5
    singular_determinant: float = 5e-22
    trust_radius: Optional[float] = None   # bohr; None = steps are not capped

    def validate(self) -> "CPSearchParameters":
        for f in fields(self):
            value = getattr(self, f.name)
            if value is None and f.name == "trust_radius":
                continue
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"CP search parameter '{f.name}' must be a positive number, not {value}.")
        if int(self.max_iterations) != self.max_iterations:
            raise ValueError(f"CP search parameter 'max_iterations' must be an integer, not {self.max_iterations}.")
        return self


@dataclass
class PathParameters:
    max_points: int = 451
    step: float = 0.03            # bohr
    stop_distance: float = 0.05   # bohr

    def validate(self) -> "PathParameters":
        for f in fields(self):
            value = getattr(self, f.name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"Path parameter '{f.name}' must be a positive number, not {value}.")
        if int(self.max_points) != self.max_points or self.max_points < 2:
            raise ValueError(f"Path parameter 'max_points' must be an integer of at least 2, not {self.max_points}.")
        return self


# ----------------------------------------------------------------------
# Critical points
# ----------------------------------------------------------------------
@dataclass
class CPProperties:
    """Values at a CP in a.u. (Multiwfn's definitions; V = -G - K, H = G + V = -K)."""
    rho: float
    laplacian: float
    gradient_norm: float
    g: float               # Lagrangian kinetic energy density G(r)
    k: float               # Hamiltonian kinetic energy density K(r) = G - lap/4
    v: float               # potential energy density V(r) = -G - K
    h: float               # total energy density H(r) = G + V = -K
    eigenvalues: np.ndarray    # lambda1 <= lambda2 <= lambda3
    eigenvectors: np.ndarray   # columns
    ellipticity: float     # lambda1/lambda2 - 1
    eta: float             # |lambda1|/lambda3

    @property
    def v_over_g(self) -> float:
        return abs(self.v) / self.g if self.g else float("nan")

    @property
    def h_over_rho(self) -> float:
        return self.h / self.rho if self.rho else float("nan")


@dataclass
class CriticalPoint:
    position: np.ndarray            # (3,) bohr
    kind: str                       # one of CP_TYPES
    properties: CPProperties
    assignment: Optional[Tuple] = None   # ("atom", i) | ("nna",) | ("pair", i, j); 0-based atoms, -1 = NNA
    index: int = 0                  # 1-based, after ordering


def classify(eigenvalues: Sequence[float], zero: float = ZERO_EIGENVALUE) -> Optional[str]:
    """(3, signature) from the eigenvalue signs; None for a degenerate point (a zero eigenvalue)."""
    lam = np.asarray(eigenvalues, dtype=float)
    if np.any(np.abs(lam) < zero):
        return None
    return CP_TYPES[int(np.sum(lam > 0))]


def _safe_ratio(a: float, b: float) -> float:
    return a / b if b != 0 else float("nan")


def cp_properties(dd: gw.DensityDerivatives, i: int) -> CPProperties:
    """Properties of point i of a ``density_derivatives(order=2)`` result."""
    hess = dd.hessian[i]
    lam, vec = np.linalg.eigh(hess)
    lap = float(np.trace(hess))
    g = float(dd.kinetic_g[i])
    k = g - 0.25 * lap
    v = -g - k
    return CPProperties(
        rho=float(dd.rho[i]), laplacian=lap, gradient_norm=float(np.linalg.norm(dd.gradient[i])),
        g=g, k=k, v=v, h=g + v, eigenvalues=lam, eigenvectors=vec,
        ellipticity=_safe_ratio(lam[0], lam[1]) - 1.0, eta=_safe_ratio(abs(lam[0]), lam[2]),
    )


def poincare_hopf(counts: Dict[str, int]) -> Tuple[int, bool]:
    """n(3,-3) - n(3,-1) + n(3,+1) - n(3,+3), and whether it is 1 (isolated molecule)."""
    total = counts.get(NCP, 0) - counts.get(BCP, 0) + counts.get(RCP, 0) - counts.get(CCP, 0)
    return total, total == 1


# ----------------------------------------------------------------------
# Starting points
# ----------------------------------------------------------------------
def vdw_radii_bohr(atomic_numbers: Sequence[int]) -> np.ndarray:
    symbols = [gw.SYMBOL_BY_ATOMIC_NUMBER.get(int(z), "X") for z in atomic_numbers]
    radii = [AtomicData.vdw_radii.get(s, AtomicData.DEFAULT_VDW_RADIUS) for s in symbols]
    return np.array(radii, dtype=float) / BOHR_RADIUS_ANGSTROM


def seeds_nuclei(coords: np.ndarray) -> np.ndarray:
    return np.asarray(coords, dtype=float).reshape(-1, 3).copy()


def atom_combinations(coords: np.ndarray, radii: np.ndarray, size: int, scale: float) -> List[Tuple[int, ...]]:
    """Sorted index tuples of ``size`` atoms, every pair within scale x (r_i + r_j)."""
    coords = np.asarray(coords, dtype=float).reshape(-1, 3)
    n = len(coords)
    if size == 1:
        return [(i,) for i in range(n)]
    if n < size:
        return []
    tree = cKDTree(coords)
    candidates = tree.query_pairs(2.0 * scale * float(np.max(radii)), output_type="ndarray")
    neighbours: List[set] = [set() for _ in range(n)]
    for i, j in candidates:
        if np.linalg.norm(coords[i] - coords[j]) <= scale * (radii[i] + radii[j]):
            neighbours[i].add(j)
            neighbours[j].add(i)
    result: List[Tuple[int, ...]] = []

    def extend(clique: Tuple[int, ...], allowed: set):
        if len(clique) == size:
            result.append(clique)
            return
        for k in sorted(x for x in allowed if x > clique[-1]):
            extend(clique + (k,), allowed & neighbours[k])

    for i in range(n):
        extend((i,), neighbours[i])
    return result


def seeds_combinations(coords: np.ndarray, radii: np.ndarray, size: int, scale: float) -> np.ndarray:
    """Centres (midpoint, triangle, tetrahedron centre) of the atom combinations."""
    coords = np.asarray(coords, dtype=float).reshape(-1, 3)
    combos = atom_combinations(coords, radii, size, scale)
    if not combos:
        return np.empty((0, 3))
    return coords[np.array(combos)].mean(axis=1)


def seeds_for_mode(mode: int, coords: np.ndarray, radii: np.ndarray, scale: float) -> np.ndarray:
    if mode not in MODE_SIZES:
        raise ValueError(f"Unknown CP search mode {mode}; use 2 (nuclei), 3 (pairs), 4 (triangles) "
                         "or 5 (tetrahedra).")
    if mode == MODE_NUCLEI:
        return seeds_nuclei(coords)
    return seeds_combinations(coords, radii, MODE_SIZES[mode], scale)


# ----------------------------------------------------------------------
# Newton search
# ----------------------------------------------------------------------
def newton_search(wfn: gw.FchkWavefunction, seeds: np.ndarray, params: CPSearchParameters,
                  nuclear: bool = False, use_numba: Optional[bool] = None,
                  progress: Optional[ProgressCallback] = None) -> List[np.ndarray]:
    """Converged Newton points from every seed, in seed order (batched over live seeds).

    Step = -step_scale H^-1 grad rho (capped at trust_radius if set). A seed is
    dropped when |det H| < singular_determinant, when rho < MIN_CP_DENSITY or
    after max_iterations. It has converged when the step is shorter than
    displacement_tolerance and |grad rho| (at the start of the step) is below
    gradient_tolerance; the gradient test is waived for nuclear seeds.
    """
    pos = np.array(seeds, dtype=float).reshape(-1, 3)
    n = len(pos)
    live = np.ones(n, dtype=bool)
    done = np.zeros(n, dtype=bool)
    for iteration in range(int(params.max_iterations)):
        idx = np.flatnonzero(live)
        if idx.size == 0:
            break
        dd = gw.density_derivatives(wfn, pos[idx], order=2, use_numba=use_numba)
        det = np.linalg.det(dd.hessian)
        ok = (np.abs(det) >= params.singular_determinant) & (dd.rho >= MIN_CP_DENSITY)
        live[idx[~ok]] = False
        idx = idx[ok]
        if idx.size == 0:
            break
        grad = dd.gradient[ok]
        step = -params.step_scale * np.linalg.solve(dd.hessian[ok], grad[..., None])[..., 0]
        length = np.linalg.norm(step, axis=1)
        if params.trust_radius is not None:
            over = length > params.trust_radius
            step[over] *= (params.trust_radius / length[over])[:, None]
            length[over] = params.trust_radius
        pos[idx] += step
        converged = length < params.displacement_tolerance
        if not nuclear:
            converged &= np.linalg.norm(grad, axis=1) < params.gradient_tolerance
        done[idx[converged]] = True
        live[idx[converged]] = False
        _report(progress, (iteration + 1) / params.max_iterations)
    _report(progress, 1.0)
    return [pos[i].copy() for i in np.flatnonzero(done)]


def merge_critical_points(candidates: Sequence[np.ndarray], existing: Sequence[np.ndarray],
                          min_distance: float) -> List[np.ndarray]:
    """Candidates at least ``min_distance`` from every existing and earlier-accepted point."""
    kept: List[np.ndarray] = []
    known = [np.asarray(p, dtype=float) for p in existing]
    for c in candidates:
        c = np.asarray(c, dtype=float)
        if all(np.linalg.norm(c - p) >= min_distance for p in known):
            kept.append(c)
            known.append(c)
    return kept


def assign_nuclear_cps(cps: Sequence[CriticalPoint], coords: np.ndarray) -> None:
    """Each (3,-3) gets its nearest nucleus, or ('nna',) beyond NNA_DISTANCE."""
    coords = np.asarray(coords, dtype=float).reshape(-1, 3)
    for cp in cps:
        if cp.kind != NCP:
            continue
        d = np.linalg.norm(coords - cp.position, axis=1)
        i = int(np.argmin(d))
        cp.assignment = ("atom", i) if d[i] <= NNA_DISTANCE else ("nna",)


def _ordered(cps: List[CriticalPoint]) -> List[CriticalPoint]:
    """(3,-3) by assigned atom (NNAs last), then (3,-1), (3,+1), (3,+3) in discovery order."""
    def ncp_key(item):
        position, cp = item
        a = cp.assignment
        return (0, a[1], position) if a and a[0] == "atom" else (1, 0, position)

    ncps = [cp for _, cp in sorted(((k, c) for k, c in enumerate(cps) if c.kind == NCP), key=ncp_key)]
    ordered = ncps + [c for t in CP_TYPES[1:] for c in cps if c.kind == t]
    for number, cp in enumerate(ordered, start=1):
        cp.index = number
    return ordered


def find_critical_points(wfn: gw.FchkWavefunction, modes: Sequence[int], params: CPSearchParameters,
                         use_numba: Optional[bool] = None,
                         progress: Optional[ProgressCallback] = None) -> Tuple[List[CriticalPoint], List[str]]:
    """All CPs from the given seed modes (run in ascending order), ordered and numbered."""
    params.validate()
    modes = sorted(set(modes))
    if not modes:
        raise ValueError("Select at least one CP search mode.")
    radii = vdw_radii_bohr(wfn.atomic_numbers)
    positions: List[np.ndarray] = []
    for k, mode in enumerate(modes):
        seeds = seeds_for_mode(mode, wfn.coordinates, radii, params.vdw_cutoff_scale)
        if len(seeds) == 0:
            continue
        found = newton_search(wfn, seeds, params, nuclear=(mode == MODE_NUCLEI), use_numba=use_numba,
                              progress=lambda f, k=k: _report(progress, (k + f) / len(modes)))
        positions += merge_critical_points(found, positions, params.min_cp_distance)
    notes: List[str] = []
    cps: List[CriticalPoint] = []
    if positions:
        dd = gw.density_derivatives(wfn, np.array(positions), order=2, use_numba=use_numba)
        degenerate = 0
        for i, p in enumerate(positions):
            props = cp_properties(dd, i)
            kind = classify(props.eigenvalues)
            if kind is None:
                degenerate += 1
                continue
            cps.append(CriticalPoint(position=p, kind=kind, properties=props))
        if degenerate:
            notes.append(f"{degenerate} converged point(s) had a zero Hessian eigenvalue (degenerate CP) "
                         "and are not reported.")
    assign_nuclear_cps(cps, wfn.coordinates)
    _report(progress, 1.0)
    return _ordered(cps), notes


def cp_counts(cps: Sequence[CriticalPoint]) -> Dict[str, int]:
    return {t: sum(1 for c in cps if c.kind == t) for t in CP_TYPES}
