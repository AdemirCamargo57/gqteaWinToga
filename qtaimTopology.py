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


# ----------------------------------------------------------------------
# Gradient paths
# ----------------------------------------------------------------------
MIN_PATH_GRADIENT = 1e-12   # |grad rho| below it: the path has stalled (no direction to follow)


@dataclass
class TopologyPath:
    kind: str              # "bond" (3,-1) -> (3,-3) or "ring" (3,+1) -> (3,+3)
    start_cp: int          # 1-based CP index
    end_cp: int            # 1-based CP index
    points: np.ndarray     # (m, 3) bohr, from the start CP's position to the end CP's position

    @property
    def length(self) -> float:
        return float(np.sum(np.linalg.norm(np.diff(self.points, axis=0), axis=1)))


@dataclass(frozen=True)
class Plane:
    """A plane through ``origin`` spanned by the orthonormal ``u`` and ``v`` (bohr)."""
    origin: np.ndarray
    u: np.ndarray
    v: np.ndarray

    @property
    def normal(self) -> np.ndarray:
        return np.cross(self.u, self.v)

    def to_plane(self, points: np.ndarray) -> np.ndarray:
        """(m, 2) in-plane coordinates of 3D points (bohr)."""
        d = np.asarray(points, dtype=float).reshape(-1, 3) - self.origin
        return np.stack([d @ self.u, d @ self.v], axis=1)

    def distance(self, points: np.ndarray) -> np.ndarray:
        """Signed distance of 3D points from the plane (bohr)."""
        return (np.asarray(points, dtype=float).reshape(-1, 3) - self.origin) @ self.normal

    def to_space(self, st: np.ndarray) -> np.ndarray:
        st = np.asarray(st, dtype=float).reshape(-1, 2)
        return self.origin + st[:, :1] * self.u + st[:, 1:] * self.v


def trace_paths(wfn: gw.FchkWavefunction, starts: np.ndarray, direction: float, targets: np.ndarray,
                params: PathParameters, plane: Optional[Plane] = None, max_points: Optional[int] = None,
                bounds: Optional[Tuple[float, float, float, float]] = None, keep_unfinished: bool = False,
                use_numba: Optional[bool] = None) -> List[Optional[Tuple[np.ndarray, int]]]:
    """Follow +grad rho (direction +1) or -grad rho (-1) from each start, all paths in one batch.

    Classical RK4 on the unit gradient with the fixed step ``params.step``, so
    points are evenly spaced, as in Multiwfn. A path ends when it comes within
    ``stop_distance`` of a target (the target's position is appended): the
    result is ``(points (m, 3), target index)``, points starting at the start.
    Reaching ``max_points`` (default params.max_points), stalling (|grad| <
    MIN_PATH_GRADIENT) or leaving ``bounds`` (plane coordinates, bohr) gives
    None, or ``(points, -1)`` with ``keep_unfinished``. With ``plane`` the
    gradient is projected onto it, so a start in the plane stays in it.
    """
    starts = np.asarray(starts, dtype=float).reshape(-1, 3)
    targets = np.asarray(targets, dtype=float).reshape(-1, 3)
    limit = int(max_points or params.max_points)
    tree = cKDTree(targets) if len(targets) else None
    h = float(params.step)
    paths: List[List[np.ndarray]] = [[s.copy()] for s in starts]
    results: List[Optional[Tuple[np.ndarray, int]]] = [None] * len(starts)
    pos = starts.copy()
    active = np.arange(len(starts))

    def unit(x):
        g = gw.density_derivatives(wfn, x, order=1, use_numba=use_numba).gradient
        if plane is not None:
            g = np.outer(g @ plane.u, plane.u) + np.outer(g @ plane.v, plane.v)
        norm = np.linalg.norm(g, axis=1)
        safe = np.where(norm > MIN_PATH_GRADIENT, norm, 1.0)
        return direction * g / safe[:, None], norm

    def finish(i, target):
        pts = np.array(paths[i])
        if target >= 0:
            results[i] = (np.vstack([pts, targets[target]]), int(target))
        elif keep_unfinished:
            results[i] = (pts, -1)

    def arrived(points):
        if tree is None:
            return np.full(len(points), -1)
        d, k = tree.query(points)
        return np.where(d < params.stop_distance, k, -1)

    hit = arrived(pos)
    for i in np.flatnonzero(hit >= 0):
        finish(i, hit[i])
    active = active[hit < 0]
    while active.size:
        p = pos[active]
        k1, n1 = unit(p)
        k2, _ = unit(p + 0.5 * h * k1)
        k3, _ = unit(p + 0.5 * h * k2)
        k4, _ = unit(p + h * k3)
        new = p + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        stalled = n1 <= MIN_PATH_GRADIENT
        pos[active] = new
        hit = arrived(new)
        outside = np.zeros(len(active), dtype=bool)
        if bounds is not None and plane is not None:
            st = plane.to_plane(new)
            outside = (st[:, 0] < bounds[0]) | (st[:, 0] > bounds[1]) | (st[:, 1] < bounds[2]) | (st[:, 1] > bounds[3])
        still = np.ones(len(active), dtype=bool)
        for n, i in enumerate(active):
            if stalled[n]:
                finish(i, -1)
                still[n] = False
                continue
            paths[i].append(new[n].copy())
            if hit[n] >= 0:
                finish(i, hit[n])
                still[n] = False
            elif outside[n] or len(paths[i]) >= limit:
                finish(i, -1)
                still[n] = False
        active = active[still]
    return results


def bond_paths(wfn: gw.FchkWavefunction, cps: Sequence[CriticalPoint], params: PathParameters,
               use_numba: Optional[bool] = None) -> Tuple[List[TopologyPath], List[str]]:
    """Both halves of every bond path; each complete BCP gets its ("pair", i, j) assignment.

    Each (3,-1) is displaced by +-step along the eigenvector of lambda3 and the
    gradient is climbed to a (3,-3). Atoms are 0-based, -1 for a non-nuclear
    attractor, sorted with NNAs last.
    """
    params.validate()
    bcps = [c for c in cps if c.kind == BCP]
    ncps = [c for c in cps if c.kind == NCP]
    if not bcps or not ncps:
        return [], []
    starts = []
    for c in bcps:
        e3 = c.properties.eigenvectors[:, 2]
        starts += [c.position + params.step * e3, c.position - params.step * e3]
    traced = trace_paths(wfn, np.array(starts), +1.0, np.array([c.position for c in ncps]), params,
                         use_numba=use_numba)
    paths: List[TopologyPath] = []
    incomplete = 0
    for n, bcp in enumerate(bcps):
        halves = traced[2 * n:2 * n + 2]
        if any(h is None for h in halves) or halves[0][1] == halves[1][1]:
            incomplete += 1
            continue
        ends = []
        for pts, target in halves:
            ncp = ncps[target]
            paths.append(TopologyPath("bond", bcp.index, ncp.index, np.vstack([bcp.position, pts])))
            a = ncp.assignment
            ends.append(a[1] if a and a[0] == "atom" else -1)
        bcp.assignment = ("pair", *sorted(ends, key=lambda x: (x < 0, x)))
    notes = []
    if incomplete:
        notes.append(f"{incomplete} bond critical point(s) have no complete bond path within "
                     f"{params.max_points} points of {params.step:g} bohr; they have no atom pair.")
    return paths, notes


def ring_paths(wfn: gw.FchkWavefunction, cps: Sequence[CriticalPoint], params: PathParameters,
               use_numba: Optional[bool] = None) -> List[TopologyPath]:
    """(3,+1) -> (3,+3) paths, descending along the eigenvector of lambda1.

    Only traced when the molecule has cage points; a half that does not reach
    one runs off to infinity, which is normal, and is dropped silently.
    """
    params.validate()
    rcps = [c for c in cps if c.kind == RCP]
    ccps = [c for c in cps if c.kind == CCP]
    if not rcps or not ccps:
        return []
    starts = []
    for c in rcps:
        e1 = c.properties.eigenvectors[:, 0]
        starts += [c.position + params.step * e1, c.position - params.step * e1]
    traced = trace_paths(wfn, np.array(starts), -1.0, np.array([c.position for c in ccps]), params,
                         use_numba=use_numba)
    paths = []
    for n, result in enumerate(traced):
        if result is not None:
            rcp = rcps[n // 2]
            paths.append(TopologyPath("ring", rcp.index, ccps[result[1]].index,
                                      np.vstack([rcp.position, result[0]])))
    return paths
