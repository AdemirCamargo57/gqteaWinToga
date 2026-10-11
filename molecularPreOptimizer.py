"""Geometry pre-optimization for structures drawn on the molecular design canvas.

The design tab hands this module a :class:`~molecularSketch.MoleculeSketch` --
a flat 2D drawing -- and gets back a three-dimensional structure in the same
``Frame`` shape molecularViewer already speaks, ready to be displayed, saved as
``.xyz`` and used as the starting point for a real CPMD / cp.x / ORCA run.

Three stages:

1. :func:`build_3d_seed` turns sketch coordinates into Angstrom and, when the
   drawing is flat, lifts it off its plane by a small, deterministic amount.
2. :class:`ForceField` assembles a compact UFF-flavoured potential from the
   sketch's own connectivity, with analytic gradients for every term.
3. :func:`optimize_sketch` minimizes it with ``scipy.optimize``, aligns the
   result back onto the drawing, and reports what happened in an
   :class:`OptimizationResult` -- including the geometry converted back to
   sketch coordinates, so the design canvas can adopt it.

**This is a pre-optimizer, not a calibrated force field.** Energies are in
loosely kcal/mol-like units and are meaningful only relative to each other
within one run. It is tuned for main-group, mostly organic molecules: angles
come from a VSEPR steric-number count, so a steric number of 6 is treated as
all-90-degree and an octahedron will come out distorted. Use a real
electronic-structure optimization for anything quantitative.

Nothing here imports Toga: failures are raised as ``ValueError`` with a
specific message and the UI layer turns them into dialogs.
"""
from collections import deque
from dataclasses import dataclass, field as dataclass_field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize

from help import AtomicData
from molecularSketch import RESONANCE_ORDER, MoleculeSketch

# The design canvas is drawn at 40 px to the Angstrom, so a typical 1.5 A bond
# is a comfortable 60 px on screen.
PIXELS_PER_ANGSTROM = 40.0

# Sketch coordinates run x right, y down, z away from the viewer; Angstrom
# coordinates run x right, y up, z towards the viewer. Going between them is a
# 180-degree turn about x -- (x, -y, -z) -- and deliberately not a flip of y
# alone, which is a reflection and would hand back the mirror image of any
# chiral 3D structure.
SKETCH_TO_ANGSTROM_AXES = np.array([1.0, -1.0, -1.0])

# A perfectly planar start leaves several angle and torsion gradients exactly
# zero, so the minimizer would never leave the plane. The displacement is drawn
# from a fixed-seed generator: the same sketch always gives the same structure,
# which is what makes the optimization reproducible and testable. It is applied
# only to a drawing that is flat to within SEED_PLANARITY_TOLERANCE Angstrom, so
# re-optimizing a structure that is already three-dimensional starts exactly
# where it is.
SEED_OUT_OF_PLANE_AMPLITUDE = 0.20
SEED_RANDOM_STATE = 20260920
SEED_PLANARITY_TOLERANCE = 0.05

# Bond length multipliers by bond order, applied to the sum of the covalent
# radii. Calibrated against experimental carbon-carbon bonds: single 1.52 A
# (exp. 1.54), resonance 1.40 A (exp. 1.397 in benzene), double 1.34 A
# (exp. 1.33 in ethene), triple 1.22 A (exp. 1.20 in acetylene). A delocalized
# bond sits between a single and a double, which is exactly where order 1.5
# puts it.
BOND_ORDER_LENGTH_SCALE = {
    1.0: 1.00,
    RESONANCE_ORDER: 0.92,
    2.0: 0.88,
    3.0: 0.80,
}

# Force constants. Angles use a potential in cos(theta) rather than theta,
# which removes the 1/sin(theta) singularity that a harmonic-in-theta bend hits
# whenever a fragment passes through linear -- exactly what a CO2 drawn by hand
# does on its way to 180 degrees.
BOND_FORCE_CONSTANT = 300.0
ANGLE_FORCE_CONSTANT = 150.0
LINEAR_ANGLE_FORCE_CONSTANT = 80.0
TORSION_BARRIER = 25.0
NONBONDED_EPSILON = 0.10
NONBONDED_SIGMA_SCALE = 1.80

# Ideal bond angle for each VSEPR steric number (sigma bonds + lone pairs).
IDEAL_ANGLE_BY_STERIC_NUMBER = {2: 180.0, 3: 120.0, 4: 109.47, 5: 109.47, 6: 90.0}
DEFAULT_IDEAL_ANGLE = 109.47

# Below these the cross products defining a torsion plane are meaningless, so
# the term contributes nothing rather than dividing by ~0.
TORSION_DEGENERACY_FLOOR = 1e-10


@dataclass
class OptimizationResult:
    """Outcome of one pre-optimization run.

    ``coordinates`` is always populated and always finite, even when
    ``success`` is False: a partially relaxed structure is more useful to the
    user than nothing, and the UI offers it with the failure message attached.
    """

    success: bool
    elements: List[str]
    coordinates: np.ndarray
    energy: float
    gradient_norm: float
    iterations: int
    message: str
    warnings: List[str] = dataclass_field(default_factory=list)
    # The same structure in sketch coordinates, centred where the drawing was,
    # in the sketch's atom order -- what the design canvas writes back.
    sketch_positions: Optional[np.ndarray] = None

    @property
    def atom_count(self) -> int:
        return len(self.elements)

    def to_frame(self) -> List[Tuple[str, Tuple[float, float, float]]]:
        """The structure in the ``Frame`` shape molecularViewer renders."""
        return [
            (element, (float(x), float(y), float(z)))
            for element, (x, y, z) in zip(self.elements, self.coordinates)
        ]

    def to_xyz_text(self, comment: str = "Generated by the gQTEA molecular design tool") -> str:
        lines = [str(self.atom_count), comment]
        for element, (x, y, z) in zip(self.elements, self.coordinates):
            lines.append(f"{element:<3s} {x: .8f} {y: .8f} {z: .8f}")
        return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------
# Stage 1: seeding
# ----------------------------------------------------------------------
def sketch_to_angstrom(points) -> np.ndarray:
    """Sketch coordinates (pixels, y down, z away) to Angstrom (y up, z towards)."""
    return np.asarray(points, dtype=float).reshape(-1, 3) * SKETCH_TO_ANGSTROM_AXES / PIXELS_PER_ANGSTROM


def angstrom_to_sketch(points) -> np.ndarray:
    """The inverse of :func:`sketch_to_angstrom`."""
    return np.asarray(points, dtype=float).reshape(-1, 3) * SKETCH_TO_ANGSTROM_AXES * PIXELS_PER_ANGSTROM


def reference_bond_length(element_i: str, element_j: str, order: float = 1.0) -> Optional[float]:
    """Force-field reference length in Angstrom, or ``None`` without radii.

    The sum of the covalent radii scaled by the bond order -- the same length
    the bond-stretch term relaxes towards, which is also how far from its
    anchor the design canvas places an atom grown out of another.
    """
    radius_i = AtomicData.covalent_radii.get(element_i)
    radius_j = AtomicData.covalent_radii.get(element_j)
    if radius_i is None or radius_j is None:
        return None
    return (radius_i + radius_j) * BOND_ORDER_LENGTH_SCALE.get(float(order), 1.0)


def _plane_normal(centred: np.ndarray) -> np.ndarray:
    """Unit normal of the best-fit plane, signed so the result is deterministic."""
    _, _, vectors = np.linalg.svd(centred, full_matrices=True)
    normal = vectors[-1]
    if normal[np.argmax(np.abs(normal))] < 0.0:
        normal = -normal
    return normal


def build_3d_seed(sketch: MoleculeSketch) -> Tuple[List[str], np.ndarray]:
    """Turn the drawing into a starting 3D geometry, centred on the origin.

    Sketch coordinates become Angstrom through a rotation (see
    :data:`SKETCH_TO_ANGSTROM_AXES`), so the structure is neither mirrored nor
    upside down against what the user sees. A drawing that is still flat --
    anything drawn without rotating, or drawn flat on a rotated view -- is
    lifted off its own plane by a small deterministic displacement; one that
    is already three-dimensional is used exactly as drawn.
    """
    elements = [atom.element for atom in sketch.atoms]
    if not elements:
        return [], np.zeros((0, 3), dtype=float)

    coordinates = sketch_to_angstrom(sketch.positions())
    coordinates -= coordinates.mean(axis=0)

    normal = _plane_normal(coordinates)
    if np.max(np.abs(coordinates @ normal)) < SEED_PLANARITY_TOLERANCE:
        generator = np.random.default_rng(SEED_RANDOM_STATE)
        out_of_plane = generator.normal(0.0, SEED_OUT_OF_PLANE_AMPLITUDE, size=len(elements))
        coordinates = coordinates + out_of_plane[:, None] * normal
        coordinates -= coordinates.mean(axis=0)
    return elements, coordinates


def align_onto(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rotate centred ``mobile`` to best fit centred ``target`` (Kabsch).

    Proper rotations only: the determinant correction stops the fit from
    "improving" the overlap by reflecting the structure into its mirror image.
    """
    if len(mobile) < 2:
        return mobile.copy()
    covariance = mobile.T @ target
    left, _, right = np.linalg.svd(covariance)
    handedness = np.sign(np.linalg.det(right.T @ left.T)) or 1.0
    correction = np.diag([1.0, 1.0, handedness])
    rotation = right.T @ correction @ left.T
    return mobile @ rotation.T


def ideal_angle_degrees(element: str, connections: int, bond_order_sum: int) -> float:
    """VSEPR ideal angle at an atom: sigma bonds plus lone pairs.

    Counting lone pairs is what distinguishes water (two bonds, two lone pairs,
    steric number 4, bent) from carbon dioxide (two double bonds, no lone pairs,
    steric number 2, linear). An element with no valence-electron entry falls
    back to counting bonds alone.

    ``bond_order_sum`` may be fractional when resonance bonds are involved, and
    the floor division handles that without a special case: an aromatic carbon
    (1.5 + 1.5 + 1 = 4, no lone pairs, three connections) is trigonal, and a
    pyridine nitrogen (1.5 + 1.5 = 3, one lone pair, two connections) is
    trigonal too rather than linear.
    """
    valence_electrons = AtomicData.valence_electrons.get(element)
    if valence_electrons is None:
        lone_pairs = 0
    else:
        lone_pairs = max(0, int((valence_electrons - bond_order_sum) // 2))
    steric_number = connections + lone_pairs
    return IDEAL_ANGLE_BY_STERIC_NUMBER.get(steric_number, DEFAULT_IDEAL_ANGLE)


# ----------------------------------------------------------------------
# Stage 2: the force field
# ----------------------------------------------------------------------
class ForceField:
    """A compact UFF-flavoured potential built from a sketch's connectivity.

    Four terms, each with an analytic gradient (checked against finite
    differences in tests/test_molecularPreOptimizer.py):

    * bond stretch, harmonic about a bond-order-dependent reference length;
    * angle bend, harmonic in ``cos(theta)`` about a VSEPR ideal angle;
    * torsion, a twofold barrier on double bonds only -- enough to keep alkenes
      and conjugated rings flat, without the parameter zoo a general torsion
      treatment would need;
    * non-bonded, purely repulsive, between atoms at least four bonds apart.
    """

    def __init__(
        self,
        elements: Sequence[str],
        bond_terms: Sequence[Tuple[int, int, float]],
        angle_terms: Sequence[Tuple[int, int, int, float, bool]],
        torsion_terms: Sequence[Tuple[int, int, int, int, float]],
        nonbonded_pairs: Sequence[Tuple[int, int, float]],
    ):
        self.elements = list(elements)
        self.bond_terms = list(bond_terms)
        self.angle_terms = list(angle_terms)
        self.torsion_terms = list(torsion_terms)
        self.nonbonded_pairs = list(nonbonded_pairs)

    # ------------------------------------------------------------------
    # Assembly
    # ------------------------------------------------------------------
    @classmethod
    def from_sketch(cls, sketch: MoleculeSketch, elements: Optional[Sequence[str]] = None) -> "ForceField":
        atom_ids = sketch.ordered_atom_ids()
        index_of = {atom_id: index for index, atom_id in enumerate(atom_ids)}
        if elements is None:
            elements = [atom.element for atom in sketch.atoms]

        bonds = [
            (index_of[bond.atom_i], index_of[bond.atom_j], float(bond.order))
            for bond in sketch.bonds
        ]

        bond_terms = []
        for i, j, order in bonds:
            length = reference_bond_length(elements[i], elements[j], order)
            if length is None:
                raise ValueError(
                    f"No covalent radius is known for {elements[i]} or {elements[j]}."
                )
            bond_terms.append((i, j, length))

        neighbours: Dict[int, List[int]] = {index: [] for index in range(len(elements))}
        order_sum: Dict[int, int] = {index: 0 for index in range(len(elements))}
        for i, j, order in bonds:
            neighbours[i].append(j)
            neighbours[j].append(i)
            order_sum[i] += order
            order_sum[j] += order

        angle_terms = []
        for centre, attached in neighbours.items():
            if len(attached) < 2:
                continue
            ideal = ideal_angle_degrees(elements[centre], len(attached), order_sum[centre])
            is_linear = ideal >= 179.999
            cos_ideal = float(np.cos(np.radians(ideal)))
            for position, first in enumerate(attached):
                for second in attached[position + 1:]:
                    angle_terms.append((first, centre, second, cos_ideal, is_linear))

        # Torsions exist only across multiple and delocalized bonds: that is
        # where planarity is a real constraint rather than a parameter guess.
        # Resonance is included, so an aromatic ring is flattened by every one
        # of its bonds instead of being free to pucker.
        torsion_terms = []
        for i, j, order in bonds:
            if order < RESONANCE_ORDER:
                continue
            for a in neighbours[i]:
                if a == j:
                    continue
                for b in neighbours[j]:
                    if b == i or b == a:
                        continue
                    torsion_terms.append((a, i, j, b, TORSION_BARRIER))

        nonbonded_pairs = []
        separation = cls._bond_separations(neighbours, len(elements))
        for i in range(len(elements)):
            for j in range(i + 1, len(elements)):
                if separation[i][j] <= 3:
                    continue
                radius_i = AtomicData.covalent_radii.get(elements[i], 1.0)
                radius_j = AtomicData.covalent_radii.get(elements[j], 1.0)
                nonbonded_pairs.append(
                    (i, j, NONBONDED_SIGMA_SCALE * (radius_i + radius_j))
                )

        return cls(elements, bond_terms, angle_terms, torsion_terms, nonbonded_pairs)

    @staticmethod
    def _bond_separations(neighbours: Dict[int, List[int]], count: int) -> List[List[int]]:
        """Shortest path length in bonds between every pair (BFS from each atom).

        Unreachable pairs get a large value, so separate fragments repel each
        other instead of being free to overlap.
        """
        unreachable = count + 10
        separations = [[unreachable] * count for _ in range(count)]
        for start in range(count):
            separations[start][start] = 0
            queue = deque([start])
            while queue:
                current = queue.popleft()
                for neighbour in neighbours[current]:
                    if separations[start][neighbour] == unreachable:
                        separations[start][neighbour] = separations[start][current] + 1
                        queue.append(neighbour)
        return separations

    # ------------------------------------------------------------------
    # Energy and gradient
    # ------------------------------------------------------------------
    def energy_and_gradient(self, flat_coordinates: np.ndarray) -> Tuple[float, np.ndarray]:
        coordinates = np.asarray(flat_coordinates, dtype=float).reshape(-1, 3)
        gradient = np.zeros_like(coordinates)
        energy = 0.0
        energy += self._bond_energy(coordinates, gradient)
        energy += self._angle_energy(coordinates, gradient)
        energy += self._torsion_energy(coordinates, gradient)
        energy += self._nonbonded_energy(coordinates, gradient)
        return float(energy), gradient.ravel()

    def energy(self, flat_coordinates: np.ndarray) -> float:
        return self.energy_and_gradient(flat_coordinates)[0]

    def _bond_energy(self, coordinates: np.ndarray, gradient: np.ndarray) -> float:
        if not self.bond_terms:
            return 0.0
        i_index = np.array([term[0] for term in self.bond_terms])
        j_index = np.array([term[1] for term in self.bond_terms])
        reference = np.array([term[2] for term in self.bond_terms])

        separation = coordinates[i_index] - coordinates[j_index]
        distance = np.linalg.norm(separation, axis=1)
        safe_distance = np.where(distance > 0.0, distance, 1.0)
        stretch = distance - reference

        energy = float(BOND_FORCE_CONSTANT * np.sum(stretch ** 2))
        contribution = (
            (2.0 * BOND_FORCE_CONSTANT * stretch / safe_distance)[:, None] * separation
        )
        np.add.at(gradient, i_index, contribution)
        np.add.at(gradient, j_index, -contribution)
        return energy

    def _angle_energy(self, coordinates: np.ndarray, gradient: np.ndarray) -> float:
        if not self.angle_terms:
            return 0.0
        first = np.array([term[0] for term in self.angle_terms])
        centre = np.array([term[1] for term in self.angle_terms])
        second = np.array([term[2] for term in self.angle_terms])
        cos_ideal = np.array([term[3] for term in self.angle_terms])
        is_linear = np.array([term[4] for term in self.angle_terms], dtype=bool)

        u = coordinates[first] - coordinates[centre]
        v = coordinates[second] - coordinates[centre]
        norm_u = np.linalg.norm(u, axis=1)
        norm_v = np.linalg.norm(v, axis=1)
        safe_u = np.where(norm_u > 0.0, norm_u, 1.0)
        safe_v = np.where(norm_v > 0.0, norm_v, 1.0)
        cosine = np.sum(u * v, axis=1) / (safe_u * safe_v)

        # A linear target uses (1 + cos) so that its minimum stays first-order
        # steep; every other target is harmonic in cos(theta).
        energy_terms = np.where(
            is_linear,
            LINEAR_ANGLE_FORCE_CONSTANT * (1.0 + cosine),
            ANGLE_FORCE_CONSTANT * (cosine - cos_ideal) ** 2,
        )
        d_energy_d_cos = np.where(
            is_linear,
            LINEAR_ANGLE_FORCE_CONSTANT,
            2.0 * ANGLE_FORCE_CONSTANT * (cosine - cos_ideal),
        )

        d_cos_d_u = v / (safe_u * safe_v)[:, None] - (cosine / safe_u ** 2)[:, None] * u
        d_cos_d_v = u / (safe_u * safe_v)[:, None] - (cosine / safe_v ** 2)[:, None] * v

        grad_first = d_energy_d_cos[:, None] * d_cos_d_u
        grad_second = d_energy_d_cos[:, None] * d_cos_d_v
        np.add.at(gradient, first, grad_first)
        np.add.at(gradient, second, grad_second)
        np.add.at(gradient, centre, -(grad_first + grad_second))
        return float(np.sum(energy_terms))

    def _torsion_energy(self, coordinates: np.ndarray, gradient: np.ndarray) -> float:
        if not self.torsion_terms:
            return 0.0
        energy = 0.0
        for a, i, j, b, barrier in self.torsion_terms:
            b1 = coordinates[i] - coordinates[a]
            b2 = coordinates[j] - coordinates[i]
            b3 = coordinates[b] - coordinates[j]

            n1 = np.cross(b1, b2)
            n2 = np.cross(b2, b3)
            n1_squared = float(np.dot(n1, n1))
            n2_squared = float(np.dot(n2, n2))
            b2_length = float(np.linalg.norm(b2))
            if n1_squared < TORSION_DEGENERACY_FLOOR or n2_squared < TORSION_DEGENERACY_FLOOR:
                continue  # three atoms in a line define no plane to twist about
            if b2_length < TORSION_DEGENERACY_FLOOR:
                continue

            phi = float(
                np.arctan2(float(np.dot(np.cross(n1, n2), b2 / b2_length)), float(np.dot(n1, n2)))
            )
            energy += 0.5 * barrier * (1.0 - np.cos(2.0 * phi))
            d_energy_d_phi = barrier * np.sin(2.0 * phi)

            # The two end atoms move phi directly; the two central atoms pick up
            # a projection of each. These coefficients are sign-convention
            # sensitive (they depend on which way b1 and b3 point), so they are
            # pinned by the finite-difference gradient test rather than trusted.
            d_phi_d_a = -(b2_length / n1_squared) * n1
            d_phi_d_b = (b2_length / n2_squared) * n2
            p = float(np.dot(b1, b2)) / (b2_length ** 2)
            q = float(np.dot(b3, b2)) / (b2_length ** 2)
            d_phi_d_i = -(1.0 + p) * d_phi_d_a + q * d_phi_d_b
            d_phi_d_j = p * d_phi_d_a - (1.0 + q) * d_phi_d_b

            gradient[a] += d_energy_d_phi * d_phi_d_a
            gradient[i] += d_energy_d_phi * d_phi_d_i
            gradient[j] += d_energy_d_phi * d_phi_d_j
            gradient[b] += d_energy_d_phi * d_phi_d_b
        return float(energy)

    def _nonbonded_energy(self, coordinates: np.ndarray, gradient: np.ndarray) -> float:
        if not self.nonbonded_pairs:
            return 0.0
        i_index = np.array([pair[0] for pair in self.nonbonded_pairs])
        j_index = np.array([pair[1] for pair in self.nonbonded_pairs])
        sigma = np.array([pair[2] for pair in self.nonbonded_pairs])

        separation = coordinates[i_index] - coordinates[j_index]
        distance = np.linalg.norm(separation, axis=1)
        safe_distance = np.where(distance > 1e-8, distance, 1e-8)
        ratio = (sigma / safe_distance) ** 12

        energy = float(NONBONDED_EPSILON * np.sum(ratio))
        # dE/dr = -12 * eps * (sigma/r)^12 / r, projected along the separation.
        contribution = (
            (-12.0 * NONBONDED_EPSILON * ratio / safe_distance ** 2)[:, None] * separation
        )
        np.add.at(gradient, i_index, contribution)
        np.add.at(gradient, j_index, -contribution)
        return energy


# ----------------------------------------------------------------------
# Stage 3: minimization
# ----------------------------------------------------------------------
def optimize_sketch(sketch: MoleculeSketch, max_iterations: int = 500) -> OptimizationResult:
    """Pre-optimize a drawn structure and report what happened.

    Raises ``ValueError`` when the sketch cannot produce a geometry at all
    (no atoms, a stray unbonded atom, an element with no covalent radius);
    validation *warnings* are carried into the result instead, since an unusual
    valence is still something the user may have meant.

    A run that hits the iteration limit is reported with ``success=False`` and
    the minimizer's own message -- never as an exception -- and its partially
    relaxed coordinates are returned all the same.
    """
    errors, warnings = sketch.validate()
    if errors:
        raise ValueError(" ".join(errors))

    elements, seed = build_3d_seed(sketch)
    field = ForceField.from_sketch(sketch, elements)

    outcome = minimize(
        field.energy_and_gradient,
        seed.ravel(),
        jac=True,
        method="L-BFGS-B",
        options={"maxiter": int(max_iterations)},
    )

    coordinates = np.asarray(outcome.x, dtype=float).reshape(-1, 3)
    if not np.all(np.isfinite(coordinates)):
        coordinates = seed
    coordinates = coordinates - coordinates.mean(axis=0)

    # Settle the result onto what was drawn (minus the seed's jitter), so the
    # canvas shows the molecule relaxing in place rather than jumping or
    # spinning, and put it back where the drawing was on the canvas.
    drawn_positions = sketch.positions()
    drawn = sketch_to_angstrom(drawn_positions)
    coordinates = align_onto(coordinates, drawn - drawn.mean(axis=0))
    sketch_positions = angstrom_to_sketch(coordinates) + drawn_positions.mean(axis=0)

    message = outcome.message
    if isinstance(message, bytes):  # older SciPy builds report bytes
        message = message.decode("utf-8", errors="replace")

    return OptimizationResult(
        success=bool(outcome.success),
        elements=elements,
        coordinates=coordinates,
        energy=float(outcome.fun),
        gradient_norm=float(np.linalg.norm(outcome.jac)),
        iterations=int(getattr(outcome, "nit", 0)),
        message=str(message),
        warnings=list(warnings),
        sketch_positions=sketch_positions,
    )
