"""Periodic pre-optimization of a packed solvation box.

Random insertion leaves molecules just outside the builder's clash cutoff, and
when the target density is not reached, the molecules it did place are often
crowded. This module relaxes such a box before it is used as a starting
structure, while keeping every molecule whole and inside the periodic cell.

Pieces
------
* :class:`PeriodicCell` -- orthorhombic cell, minimum image, wrapping.
* :class:`PackedSystem` -- atoms plus **which molecule each atom belongs to**.
  The builders know that from construction, so molecular connectivity never
  has to be guessed from distances across a periodic boundary.
* :class:`OverlapPenalty` -- the intermolecular term, a Packmol-style penalty
  ``K * (d^2 - sigma^2)^2`` for ``d < sigma``, with ``sigma = vdw_scale *
  (r_i + r_j)`` (+ the solute clearance for solute-solvent pairs): exactly the
  builder's own acceptance rule, so an energy of zero means no clash is left.
  It has a finite range, which makes periodic boundaries exact: neighbour pairs
  come from ``cKDTree(..., boxsize=L)`` and every distance is recomputed with
  the minimum image. Pairs inside one molecule are excluded.
* Two models: :class:`RigidBodyModel` (6 degrees of freedom per molecule, the
  internal geometry is preserved to round-off) and :class:`FlexibleModel`
  (every atom moves; bonds and angles are held near the *input* geometry by
  :meth:`molecularPreOptimizer.ForceField.from_geometry`).
* :func:`pre_optimize` -- L-BFGS-B, then molecule-aware wrapping
  (:func:`wrap_molecules`): each molecule moves as a whole so that its centre
  of mass lies in ``[0, L)``; no molecule is ever split across a face.

Without periodicity (the builder's minimum-image switch off) there are no
images, and a harmonic wall keeps atoms inside ``[pad, L - pad]`` -- a purely
repulsive potential would otherwise push molecules out of the box.

Pure numpy/scipy, no Toga. Failures are ``ValueError``; non-convergence is
``success=False`` with the partially relaxed structure, never an exception.
"""
from dataclasses import dataclass, field as dataclass_field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from molecularPreOptimizer import ForceField, detect_bonds

MODE_RIGID = "rigid"
MODE_FLEXIBLE = "flexible"
MODES = (MODE_RIGID, MODE_FLEXIBLE)

SOLUTE_KIND = "solute"

DEFAULT_MAX_ITERATIONS = 500
PENALTY_CONSTANT = 1.0
# The minimizer works against sigma enlarged by this fraction, while every
# report (clash count, penalty) uses the true sigma. The penalty vanishes
# quadratically at d = sigma, so without a margin L-BFGS stops with a few
# contacts a hair inside the cutoff; with it they settle just outside.
OPTIMIZATION_MARGIN = 0.01
WALL_CONSTANT = 10.0

# Below this rotation angle the left-Jacobian coefficients use their Taylor
# series; the error is ~theta^4 / 720, far below double precision noise.
SMALL_ANGLE = 1e-3


# ----------------------------------------------------------------------
# Cell and wrapping
# ----------------------------------------------------------------------
class PeriodicCell:
    """Orthorhombic cell with side lengths ``a, b, c`` in Angstrom."""

    def __init__(self, lengths: Sequence[float], periodic: bool = True):
        self.lengths = np.asarray(lengths, dtype=float).reshape(3)
        if not np.all(np.isfinite(self.lengths)) or np.any(self.lengths <= 0.0):
            raise ValueError("Cell lengths must be positive numbers.")
        self.periodic = bool(periodic)

    def minimum_image(self, delta: np.ndarray) -> np.ndarray:
        if not self.periodic:
            return delta
        return delta - self.lengths * np.round(delta / self.lengths)

    def wrap_points(self, points: np.ndarray) -> np.ndarray:
        """Points folded into ``[0, L)`` -- per atom, for neighbour searches only."""
        wrapped = points - self.lengths * np.floor(points / self.lengths)
        # floor() can leave a tiny negative x at exactly L after rounding.
        wrapped = np.where(wrapped >= self.lengths, wrapped - self.lengths, wrapped)
        return np.where(wrapped < 0.0, 0.0, wrapped)


def molecule_centres(coordinates: np.ndarray, masses: np.ndarray,
                     slices: Sequence[Tuple[int, int]]) -> np.ndarray:
    """Mass-weighted centre of each molecule (geometric if its mass is zero)."""
    centres = np.empty((len(slices), 3))
    for index, (start, stop) in enumerate(slices):
        weights = masses[start:stop]
        total = float(np.sum(weights))
        if total > 0.0:
            centres[index] = weights @ coordinates[start:stop] / total
        else:
            centres[index] = coordinates[start:stop].mean(axis=0)
    return centres


def wrap_molecules(coordinates: np.ndarray, masses: np.ndarray,
                   slices: Sequence[Tuple[int, int]], cell: PeriodicCell) -> np.ndarray:
    """Translate each molecule as a whole so its centre of mass is in ``[0, L)``.

    The shift is a whole number of cell vectors per molecule, so intramolecular
    vectors -- and therefore the molecule itself -- are untouched. Atoms may
    stick out past a face; that is the convention that keeps molecules whole.
    """
    coordinates = np.array(coordinates, dtype=float, copy=True)
    if not cell.periodic:
        return coordinates
    centres = molecule_centres(coordinates, masses, slices)
    shifts = -cell.lengths * np.floor(centres / cell.lengths)
    for (start, stop), shift in zip(slices, shifts):
        coordinates[start:stop] += shift
    return coordinates


def make_whole(coordinates: np.ndarray, slices: Sequence[Tuple[int, int]],
               cell: PeriodicCell) -> np.ndarray:
    """Undo per-atom wrapping: each atom is placed at the minimum image of the
    previous atom of its molecule. Correct whenever consecutive atoms of a
    molecule are less than half a cell apart (any molecule smaller than half
    the box, which the builders already require)."""
    coordinates = np.array(coordinates, dtype=float, copy=True)
    if not cell.periodic:
        return coordinates
    for start, stop in slices:
        for atom in range(start + 1, stop):
            delta = cell.minimum_image(coordinates[atom] - coordinates[atom - 1])
            coordinates[atom] = coordinates[atom - 1] + delta
    return coordinates


# ----------------------------------------------------------------------
# System
# ----------------------------------------------------------------------
@dataclass
class PackedSystem:
    elements: List[str]
    coordinates: np.ndarray
    molecule_slices: List[Tuple[int, int]]
    molecule_kinds: List[str]
    vdw_radii: np.ndarray
    masses: np.ndarray

    def __post_init__(self):
        self.coordinates = np.asarray(self.coordinates, dtype=float).reshape(-1, 3)
        self.vdw_radii = np.asarray(self.vdw_radii, dtype=float).reshape(-1)
        self.masses = np.asarray(self.masses, dtype=float).reshape(-1)
        count = len(self.elements)
        if not (len(self.coordinates) == len(self.vdw_radii) == len(self.masses) == count):
            raise ValueError("Elements, coordinates, radii and masses must have one entry per atom.")
        if len(self.molecule_kinds) != len(self.molecule_slices):
            raise ValueError("Every molecule needs a kind.")
        expected = 0
        for start, stop in self.molecule_slices:
            if start != expected or stop <= start:
                raise ValueError("Molecule slices must be contiguous, non-empty and in order.")
            expected = stop
        if expected != count:
            raise ValueError("Molecule slices must cover every atom exactly once.")

    @classmethod
    def from_molecules(cls, molecules: Sequence[Tuple[str, Sequence[str], np.ndarray]],
                       radius_of: Callable[[str], float],
                       mass_of: Callable[[str], float]) -> "PackedSystem":
        """``molecules`` = ``[(kind, elements, (n, 3) coordinates), ...]``."""
        elements: List[str] = []
        blocks, slices, kinds = [], [], []
        for kind, mol_elements, mol_coords in molecules:
            start = len(elements)
            elements.extend(mol_elements)
            blocks.append(np.asarray(mol_coords, dtype=float).reshape(-1, 3))
            slices.append((start, len(elements)))
            kinds.append(kind)
        coordinates = np.concatenate(blocks) if blocks else np.zeros((0, 3))
        return cls(
            elements=elements,
            coordinates=coordinates,
            molecule_slices=slices,
            molecule_kinds=kinds,
            vdw_radii=np.array([radius_of(e) for e in elements], dtype=float),
            masses=np.array([mass_of(e) for e in elements], dtype=float),
        )

    @property
    def atom_count(self) -> int:
        return len(self.elements)

    def molecule_ids(self) -> np.ndarray:
        ids = np.empty(self.atom_count, dtype=np.int64)
        for index, (start, stop) in enumerate(self.molecule_slices):
            ids[start:stop] = index
        return ids


# ----------------------------------------------------------------------
# Intermolecular penalty
# ----------------------------------------------------------------------
class OverlapPenalty:
    def __init__(self, system: PackedSystem, cell: PeriodicCell, vdw_scale: float,
                 solute_clearance: float = 0.0, wall_padding: float = 0.0,
                 penalty_constant: float = PENALTY_CONSTANT, wall_constant: float = WALL_CONSTANT,
                 sigma_factor: float = 1.0, wall_bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None):
        if vdw_scale <= 0.0:
            raise ValueError("vdW scale must be positive.")
        self.cell = cell
        self.vdw_scale = float(vdw_scale)
        self.sigma_factor = float(sigma_factor)
        self.radii = system.vdw_radii
        self.molecule_id = system.molecule_ids()
        is_solute_molecule = np.array([kind == SOLUTE_KIND for kind in system.molecule_kinds], dtype=bool)
        self.is_solute = is_solute_molecule[self.molecule_id] if len(self.molecule_id) else np.zeros(0, bool)
        self.solute_clearance = max(0.0, float(solute_clearance))
        self.wall_padding = max(0.0, float(wall_padding))
        self.k = float(penalty_constant)
        self.k_wall = float(wall_constant)

        largest = 2.0 * float(np.max(self.radii)) * self.vdw_scale if len(self.radii) else 0.0
        has_mixed = bool(np.any(self.is_solute)) and bool(np.any(~self.is_solute))
        self.search_radius = self.sigma_factor * (largest + (self.solute_clearance if has_mixed else 0.0))
        if cell.periodic and self.search_radius >= 0.5 * float(np.min(cell.lengths)):
            raise ValueError(
                f"The contact range ({self.search_radius:.2f} A) must be shorter than half the "
                f"smallest cell length ({0.5 * float(np.min(cell.lengths)):.2f} A) for the "
                "minimum-image convention to hold. Use a larger box or a smaller vdW scale."
            )
        # Non-periodic walls: per-atom (low, high) bounds, (N, 3) each, +-inf
        # where an atom has no wall on that side. By default every atom is
        # kept in [pad, L - pad]; the shared-wall box passes its own bounds
        # (no wall at the interface between its two regions).
        self.wall_low = self.wall_high = None
        if not cell.periodic:
            if wall_bounds is None:
                inner = cell.lengths - 2.0 * self.wall_padding
                if np.any(inner <= 0.0):
                    raise ValueError("The wall padding leaves no room inside the box.")
                count = len(self.radii)
                self.wall_low = np.full((count, 3), self.wall_padding)
                self.wall_high = np.tile(cell.lengths - self.wall_padding, (count, 1))
            else:
                self.wall_low = np.asarray(wall_bounds[0], dtype=float)
                self.wall_high = np.asarray(wall_bounds[1], dtype=float)
                if self.wall_low.shape != (len(self.radii), 3) or self.wall_high.shape != self.wall_low.shape:
                    raise ValueError("Wall bounds must give (low, high) for every atom and axis.")
                if np.any(self.wall_high - self.wall_low <= 0.0):
                    raise ValueError("The wall padding leaves no room inside the box.")

    def contacts(self, coordinates: np.ndarray):
        """Intermolecular pairs closer than the search radius.

        Returns ``(i, j, delta, d2, sigma)`` with ``delta`` the minimum-image
        vector from j to i.
        """
        if len(coordinates) < 2 or self.search_radius <= 0.0:
            empty = np.zeros(0, dtype=np.int64)
            return empty, empty, np.zeros((0, 3)), np.zeros(0), np.zeros(0)
        if self.cell.periodic:
            tree = cKDTree(self.cell.wrap_points(coordinates), boxsize=self.cell.lengths)
        else:
            tree = cKDTree(coordinates)
        pairs = tree.query_pairs(self.search_radius, output_type="ndarray")
        if len(pairs) == 0:
            empty = np.zeros(0, dtype=np.int64)
            return empty, empty, np.zeros((0, 3)), np.zeros(0), np.zeros(0)
        i, j = pairs[:, 0], pairs[:, 1]
        inter = self.molecule_id[i] != self.molecule_id[j]
        i, j = i[inter], j[inter]
        delta = self.cell.minimum_image(coordinates[i] - coordinates[j])
        d2 = np.einsum("ij,ij->i", delta, delta)
        sigma = self.vdw_scale * (self.radii[i] + self.radii[j])
        if self.solute_clearance > 0.0:
            sigma = sigma + self.solute_clearance * (self.is_solute[i] != self.is_solute[j])
        if self.sigma_factor != 1.0:
            sigma = self.sigma_factor * sigma
        return i, j, delta, d2, sigma

    def energy_and_gradient(self, coordinates: np.ndarray) -> Tuple[float, np.ndarray]:
        gradient = np.zeros_like(coordinates)
        i, j, delta, d2, sigma = self.contacts(coordinates)
        energy = 0.0
        if len(i):
            violation = d2 - sigma * sigma
            active = violation < 0.0
            if np.any(active):
                v = violation[active]
                energy += self.k * float(np.sum(v * v))
                contribution = (4.0 * self.k * v)[:, None] * delta[active]
                np.add.at(gradient, i[active], contribution)
                np.add.at(gradient, j[active], -contribution)
        if not self.cell.periodic:
            below = np.minimum(coordinates - self.wall_low, 0.0)
            above = np.maximum(coordinates - self.wall_high, 0.0)
            energy += self.k_wall * float(np.sum(below * below) + np.sum(above * above))
            gradient += 2.0 * self.k_wall * (below + above)
        return energy, gradient

    def report(self, coordinates: np.ndarray) -> Dict[str, float]:
        """Clash count (pairs with d < sigma), closest contact and its d/sigma."""
        i, j, delta, d2, sigma = self.contacts(coordinates)
        if len(i) == 0:
            return {"clashes": 0, "min_distance": float("nan"), "min_ratio": float("nan")}
        distance = np.sqrt(d2)
        ratio = distance / sigma
        return {
            "clashes": int(np.count_nonzero(distance < sigma)),
            "min_distance": float(np.min(distance)),
            "min_ratio": float(np.min(ratio)),
        }


class SlabRestraint:
    """Flat-bottom harmonic restraint on each molecule's centre of mass along one axis.

    ``bounds[m] = (low, high)``; zero energy while the centre of mass of
    molecule ``m`` is inside, ``k (low - z)^2`` / ``k (z - high)^2`` outside.
    The shared-wall box uses it to keep each region's molecules in their own
    slab during pre-optimization, so the interface composition is preserved.
    The force on the centre of mass is shared by the atoms in proportion to
    their mass (geometric weights for a massless molecule).
    """

    def __init__(self, system: PackedSystem, bounds, axis: int = 2, force_constant: float = WALL_CONSTANT):
        self.bounds = np.asarray(bounds, dtype=float).reshape(-1, 2)
        if len(self.bounds) != len(system.molecule_slices):
            raise ValueError("Slab restraint bounds are needed for every molecule.")
        if np.any(self.bounds[:, 1] <= self.bounds[:, 0]):
            raise ValueError("Each slab must have high > low.")
        self.axis = int(axis)
        self.k = float(force_constant)
        self.molecule_id = system.molecule_ids()
        weights = system.masses.copy()
        totals = np.zeros(len(system.molecule_slices))
        np.add.at(totals, self.molecule_id, weights)
        for index, (start, stop) in enumerate(system.molecule_slices):
            if totals[index] <= 0.0:
                weights[start:stop] = 1.0
                totals[index] = float(stop - start)
        self.fraction = weights / totals[self.molecule_id]   # d(COM)/d(x_atom)
        self.count = len(system.molecule_slices)

    def centres(self, coordinates: np.ndarray) -> np.ndarray:
        com = np.zeros(self.count)
        np.add.at(com, self.molecule_id, self.fraction * coordinates[:, self.axis])
        return com

    def energy_and_gradient(self, coordinates: np.ndarray) -> Tuple[float, np.ndarray]:
        com = self.centres(coordinates)
        below = np.minimum(com - self.bounds[:, 0], 0.0)
        above = np.maximum(com - self.bounds[:, 1], 0.0)
        energy = self.k * float(np.sum(below * below) + np.sum(above * above))
        d_com = 2.0 * self.k * (below + above)
        gradient = np.zeros_like(coordinates)
        gradient[:, self.axis] = d_com[self.molecule_id] * self.fraction
        return energy, gradient

    def outside(self, coordinates: np.ndarray) -> int:
        com = self.centres(coordinates)
        return int(np.count_nonzero((com < self.bounds[:, 0] - 1e-6) | (com > self.bounds[:, 1] + 1e-6)))


class CompositeTerm:
    """Sum of energy terms that all take and return atom coordinates."""

    def __init__(self, terms):
        self.terms = list(terms)

    def energy_and_gradient(self, coordinates: np.ndarray) -> Tuple[float, np.ndarray]:
        energy = 0.0
        gradient = np.zeros_like(coordinates)
        for term in self.terms:
            e, g = term.energy_and_gradient(coordinates)
            energy += e
            gradient += g
        return energy, gradient


# ----------------------------------------------------------------------
# Models
# ----------------------------------------------------------------------
def _cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    # np.cross on (n, 3) arrays, written out: it is called in the minimizer's
    # inner loop and avoids np.cross's per-call overhead.
    return np.stack([
        a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1],
        a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2],
        a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0],
    ], axis=1)


class RigidBodyModel:
    """Each molecule = centroid translation ``t`` + rotation vector ``w``.

    ``x_i = R(w) b_i + t`` with ``b_i`` the atom's offset from its molecule's
    centroid in the starting structure. The gradient follows from the
    left-Jacobian of SO(3): a small change ``dw`` rotates ``R b`` by
    ``J_l(w) dw``, so ``dE/dw = J_l(w)^T tau`` with the torque
    ``tau = sum (x_i - t) x g_i``; ``J_l^T`` applied to a vector is
    ``tau - A (w x tau) + B (w x (w x tau))``, ``A = (1 - cos th)/th^2``,
    ``B = (th - sin th)/th^3``. Pinned by a finite-difference test.
    """

    mode = MODE_RIGID

    def __init__(self, system: PackedSystem, penalty: OverlapPenalty):
        self.system = system
        self.penalty = penalty
        self.molecule_id = system.molecule_ids()
        self.count = len(system.molecule_slices)
        centroids = np.array([system.coordinates[s:e].mean(axis=0) for s, e in system.molecule_slices])
        self.start_centroids = centroids
        self.body = system.coordinates - centroids[self.molecule_id]

    def initial_parameters(self) -> np.ndarray:
        return np.concatenate([self.start_centroids.ravel(), np.zeros(3 * self.count)])

    def _split(self, params):
        params = np.asarray(params, dtype=float)
        return params[: 3 * self.count].reshape(-1, 3), params[3 * self.count:].reshape(-1, 3)

    def coordinates(self, params) -> np.ndarray:
        t, w = self._split(params)
        rotations = Rotation.from_rotvec(w).as_matrix()
        rotated = np.einsum("nij,nj->ni", rotations[self.molecule_id], self.body)
        return rotated + t[self.molecule_id]

    def energy_and_gradient(self, params) -> Tuple[float, np.ndarray]:
        t, w = self._split(params)
        rotations = Rotation.from_rotvec(w).as_matrix()
        arms = np.einsum("nij,nj->ni", rotations[self.molecule_id], self.body)
        x = arms + t[self.molecule_id]
        energy, g = self.penalty.energy_and_gradient(x)

        grad_t = np.zeros((self.count, 3))
        np.add.at(grad_t, self.molecule_id, g)
        torque = np.zeros((self.count, 3))
        np.add.at(torque, self.molecule_id, _cross(arms, g))

        theta2 = np.einsum("ij,ij->i", w, w)
        theta = np.sqrt(theta2)
        small = theta < SMALL_ANGLE
        safe = np.where(small, 1.0, theta)
        a_coef = np.where(small, 0.5 - theta2 / 24.0, (1.0 - np.cos(safe)) / safe ** 2)
        b_coef = np.where(small, 1.0 / 6.0 - theta2 / 120.0, (safe - np.sin(safe)) / safe ** 3)
        w_x_tau = _cross(w, torque)
        grad_w = torque - a_coef[:, None] * w_x_tau + b_coef[:, None] * _cross(w, w_x_tau)
        return energy, np.concatenate([grad_t.ravel(), grad_w.ravel()])


class FlexibleModel:
    """Every atom moves; bonded terms keep each molecule near its input shape.

    One :meth:`ForceField.from_geometry` per molecule *kind*, built on the
    first molecule of that kind and tiled over the others (the copies are
    rigid rotations of the same template, so their internal references are
    identical). Intramolecular vectors are taken as they are -- molecules are
    kept whole, so no minimum image is needed or wanted inside one.
    """

    mode = MODE_FLEXIBLE

    def __init__(self, system: PackedSystem, penalty: OverlapPenalty):
        self.system = system
        self.penalty = penalty
        templates: Dict[str, ForceField] = {}
        parts = []
        for kind, (start, stop) in zip(system.molecule_kinds, system.molecule_slices):
            if kind not in templates:
                try:
                    templates[kind] = ForceField.from_geometry(
                        system.elements[start:stop], system.coordinates[start:stop]
                    )
                except ValueError as error:
                    raise ValueError(f"Flexible pre-optimization is not possible: {error} "
                                     "Use the rigid-body mode instead.")
            parts.append((templates[kind], start))
        self.field = ForceField.combine(parts)

    def initial_parameters(self) -> np.ndarray:
        return self.system.coordinates.ravel().copy()

    def coordinates(self, params) -> np.ndarray:
        return np.asarray(params, dtype=float).reshape(-1, 3)

    def energy_and_gradient(self, params) -> Tuple[float, np.ndarray]:
        x = self.coordinates(params)
        e_field, g_field = self.field.energy_and_gradient(params)
        e_pen, g_pen = self.penalty.energy_and_gradient(x)
        return e_field + e_pen, g_field + g_pen.ravel()


# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------
@dataclass
class PreOptimizationResult:
    success: bool
    message: str
    mode: str
    periodic: bool
    iterations: int
    evaluations: int
    coordinates: np.ndarray          # molecule-wrapped, what is written out
    energy_initial: float
    energy_final: float
    report_initial: Dict[str, float]
    report_final: Dict[str, float]
    max_internal_deviation: float    # largest change of an intramolecular distance (A)
    connectivity_preserved: Optional[bool]
    cancelled: bool = False
    notes: List[str] = dataclass_field(default_factory=list)

    def summary_lines(self) -> List[str]:
        def fmt(value, spec):
            return "n/a" if value is None or not np.isfinite(value) else format(value, spec)

        if self.connectivity_preserved is None:
            connectivity = "not checked (missing covalent radii)"
        else:
            connectivity = "preserved" if self.connectivity_preserved else "CHANGED"
        lines = [
            "--- Pre-optimization ---\n",
            f"Mode: {'rigid bodies' if self.mode == MODE_RIGID else 'flexible molecules'}"
            f" ({'periodic' if self.periodic else 'non-periodic, wall-confined'})\n",
            f"Status: {'cancelled' if self.cancelled else ('converged' if self.success else 'not converged')}"
            f" -- {self.message}\n",
            f"Iterations: {self.iterations} (energy evaluations: {self.evaluations})\n",
            f"Overlap penalty: {self.energy_initial:.6g} -> {self.energy_final:.6g}\n",
            f"Clashing atom pairs (d < sigma): {self.report_initial['clashes']} -> {self.report_final['clashes']}\n",
            f"Closest intermolecular contact (A): {fmt(self.report_initial['min_distance'], '.3f')} -> "
            f"{fmt(self.report_final['min_distance'], '.3f')}\n",
            f"Closest contact / sigma: {fmt(self.report_initial['min_ratio'], '.3f')} -> "
            f"{fmt(self.report_final['min_ratio'], '.3f')}\n",
            f"Largest change of an intramolecular distance (A): {self.max_internal_deviation:.2e}\n",
            f"Molecular connectivity: {connectivity}\n",
        ]
        lines.extend(f"NOTE: {note}\n" for note in self.notes)
        return lines


def _internal_deviation(before: np.ndarray, after: np.ndarray, slices) -> float:
    worst = 0.0
    for start, stop in slices:
        if stop - start < 2:
            continue
        a = before[start:stop]
        b = after[start:stop]
        da = np.sqrt(np.sum((a[:, None] - a[None]) ** 2, axis=2))
        db = np.sqrt(np.sum((b[:, None] - b[None]) ** 2, axis=2))
        worst = max(worst, float(np.max(np.abs(da - db))))
    return worst


def _connectivity_preserved(system: PackedSystem, after: np.ndarray) -> Optional[bool]:
    try:
        for start, stop in system.molecule_slices:
            elements = system.elements[start:stop]
            if detect_bonds(elements, system.coordinates[start:stop]) != detect_bonds(elements, after[start:stop]):
                return False
    except ValueError:
        return None
    return True


class PreOptimizationCancelled(Exception):
    pass


def pre_optimize(system: PackedSystem, cell: PeriodicCell, vdw_scale: float,
                 mode: str = MODE_RIGID, max_iterations: int = DEFAULT_MAX_ITERATIONS,
                 solute_clearance: float = 0.0, wall_padding: float = 0.0,
                 callback: Optional[Callable[[np.ndarray, int, float], None]] = None,
                 should_stop: Optional[Callable[[], bool]] = None,
                 wall_bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None,
                 slab_bounds=None, slab_axis: int = 2) -> PreOptimizationResult:
    """Relax intermolecular overlaps; return the molecule-wrapped structure.

    ``callback(wrapped_coordinates, iteration, energy)`` is called after every
    minimizer iteration (the caller throttles what it does with it; a callback
    that raises is ignored). ``should_stop()`` returning True ends the run
    early with the structure reached so far (``cancelled=True``).
    """
    if mode not in MODES:
        raise ValueError(f"Unknown pre-optimization mode '{mode}'.")
    if int(max_iterations) <= 0:
        raise ValueError("Maximum pre-optimization iterations must be a positive integer.")
    if system.atom_count == 0:
        raise ValueError("There are no atoms to pre-optimize.")

    # ``reporter`` measures against the true sigma; the model minimizes
    # ``penalty``, whose sigma carries OPTIMIZATION_MARGIN.
    reporter = OverlapPenalty(system, cell, vdw_scale, solute_clearance, wall_padding,
                              wall_bounds=wall_bounds)
    penalty = OverlapPenalty(system, cell, vdw_scale, solute_clearance, wall_padding,
                             sigma_factor=1.0 + OPTIMIZATION_MARGIN, wall_bounds=wall_bounds)
    # ``slab_bounds`` (one (low, high) per molecule along ``slab_axis``) keeps
    # each molecule's centre of mass in its own region -- the shared-wall box.
    slab = SlabRestraint(system, slab_bounds, slab_axis) if slab_bounds is not None else None
    term = CompositeTerm([penalty, slab]) if slab is not None else penalty
    model = RigidBodyModel(system, term) if mode == MODE_RIGID else FlexibleModel(system, term)

    start = system.coordinates
    energy_initial = reporter.energy_and_gradient(start)[0]
    report_initial = reporter.report(start)

    state = {"iterations": 0, "evaluations": 0, "cancelled": False}

    def objective(params):
        state["evaluations"] += 1
        return model.energy_and_gradient(params)

    def on_iteration(intermediate_result):
        # The parameter name matters: SciPy passes the full OptimizeResult
        # (with .fun) only to a callback whose argument is called this.
        state["iterations"] += 1
        params = getattr(intermediate_result, "x", intermediate_result)
        if callback is not None:
            try:
                x = model.coordinates(params)
                wrapped = wrap_molecules(x, system.masses, system.molecule_slices, cell)
                callback(wrapped, state["iterations"], float(getattr(intermediate_result, "fun", np.nan)))
            except Exception:
                pass  # reporting must never abort the optimization
        if should_stop is not None and should_stop():
            state["cancelled"] = True
            raise StopIteration

    notes: List[str] = []
    slab_clear = slab is None or slab.energy_and_gradient(start)[0] == 0.0
    if report_initial["clashes"] == 0 and energy_initial == 0.0 and slab_clear:
        params = model.initial_parameters()
        success, message = True, "no overlaps to relax"
    else:
        outcome = minimize(
            objective,
            model.initial_parameters(),
            jac=True,
            method="L-BFGS-B",
            callback=on_iteration,
            options={"maxiter": int(max_iterations)},
        )
        params = outcome.x
        success = bool(outcome.success)
        message = outcome.message
        if isinstance(message, bytes):  # older SciPy builds report bytes
            message = message.decode("utf-8", errors="replace")
        message = str(message)
        if state["cancelled"]:
            success, message = False, "stopped by the user"

    final = model.coordinates(params)
    if not np.all(np.isfinite(final)):
        final = start.copy()
        success = False
        notes.append("The minimizer produced non-finite coordinates; the starting structure was kept.")
    energy_final = reporter.energy_and_gradient(final)[0]
    report_final = reporter.report(final)
    if slab is not None:
        outside = slab.outside(final)
        if outside:
            notes.append(f"{outside} molecules ended with their centre of mass outside their region.")
    if report_final["clashes"]:
        notes.append(
            f"{report_final['clashes']} atom pairs are still closer than the clash cutoff; "
            "the box may be too dense for these molecules (lower the density or the vdW scale)."
        )

    wrapped = wrap_molecules(final, system.masses, system.molecule_slices, cell)
    return PreOptimizationResult(
        success=success,
        message=message,
        mode=mode,
        periodic=cell.periodic,
        iterations=state["iterations"],
        evaluations=state["evaluations"],
        coordinates=wrapped,
        energy_initial=float(energy_initial),
        energy_final=float(energy_final),
        report_initial=report_initial,
        report_final=report_final,
        max_internal_deviation=_internal_deviation(start, final, system.molecule_slices),
        connectivity_preserved=_connectivity_preserved(system, final),
        cancelled=state["cancelled"],
        notes=notes,
    )
