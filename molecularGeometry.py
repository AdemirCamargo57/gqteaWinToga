"""Internal coordinates of a drawn molecule: measure, edit and constrain them.

The design canvas ([molecularDesign.py](molecularDesign.py)) measures a bond
length, a bond angle or a dihedral angle on 2, 3 or 4 picked atoms, sets one
to a typed value, and freezes some during pre-optimization. Everything about
those three quantities lives here, so the canvas, the edit dialog and the
pre-optimizer cannot disagree about a definition or a sign:

* :func:`measure` -- the value, in Angstrom or degrees;
* :func:`value_and_gradient` -- the value and its analytic derivative with
  respect to every coordinate, used to hold frozen parameters exactly
  (pinned against finite differences in tests/test_molecularGeometry.py);
* :func:`project_onto_constraints` -- move a geometry the least distance that
  puts every frozen parameter back on its value (Gauss-Newton, minimum norm);
* :func:`set_parameter` -- change one parameter by moving one side of the
  molecule rigidly, the way a model kit does it, refusing a ring bond (where
  no side can move alone) rather than distorting the molecule.

Coordinates are an ``(N, 3)`` array in Angstrom; atoms are 0-based indices.
The dihedral uses the IUPAC sign convention, the same as the pre-optimizer's
torsion term: looking along j->k, positive when i must turn clockwise to
eclipse l. No Toga, no canvas pixels.
"""
from collections import deque
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

DISTANCE, ANGLE, DIHEDRAL = "distance", "angle", "dihedral"
KIND_BY_ATOM_COUNT = {2: DISTANCE, 3: ANGLE, 4: DIHEDRAL}
ATOM_COUNT = {kind: count for count, kind in KIND_BY_ATOM_COUNT.items()}
KIND_LABELS = {DISTANCE: "Bond length", ANGLE: "Bond angle", DIHEDRAL: "Dihedral angle"}
UNITS = {DISTANCE: "Å", ANGLE: "°", DIHEDRAL: "°"}
DECIMALS = {DISTANCE: 4, ANGLE: 2, DIHEDRAL: 2}

# Values a typed edit may set. A bond shorter than 0.3 A or longer than 20 A
# is a typing error, not a structure.
MIN_DISTANCE = 0.3
MAX_DISTANCE = 20.0

# An angle this close to 0 or 180 degrees has no defined plane: it cannot be
# held by the constraint solver (its derivative is singular there), and it
# leaves the dihedral through it undefined.
LINEAR_TOLERANCE_DEGREES = 1.0
DIHEDRAL_LINEAR_LIMIT = 175.0

PROJECTION_TOLERANCE = 1e-10
PROJECTION_MAX_ITERATIONS = 60


class GeometryEditError(ValueError):
    """A parameter edit that cannot be made, with the reason as the message."""


@dataclass(frozen=True)
class Constraint:
    """A frozen parameter on 0-based atom indices; ``value`` in Angstrom or degrees."""

    kind: str
    atoms: Tuple[int, ...]
    value: float


# ----------------------------------------------------------------------
# Definitions
# ----------------------------------------------------------------------
def kind_for(atoms: Sequence) -> str:
    try:
        return KIND_BY_ATOM_COUNT[len(atoms)]
    except KeyError:
        raise ValueError("Pick 2 atoms for a bond length, 3 for an angle or 4 for a dihedral.") from None


def canonical_atoms(kind: str, atoms: Sequence) -> Tuple:
    """The same parameter read backwards is the same parameter: 1-2 is 2-1,
    1-2-3 is 3-2-1 and 1-2-3-4 is 4-3-2-1. Returns the smaller ordering."""
    atoms = tuple(atoms)
    return min(atoms, tuple(reversed(atoms)))


def format_value(kind: str, value: float) -> str:
    return f"{value:.{DECIMALS[kind]}f} {UNITS[kind]}" if kind == DISTANCE else \
        f"{value:.{DECIMALS[kind]}f}{UNITS[kind]}"


def wrap_degrees(value: float) -> float:
    """Into (-180, 180]."""
    wrapped = (float(value) + 180.0) % 360.0 - 180.0
    return 180.0 if wrapped == -180.0 else wrapped


def parse_parameter_value(kind: str, text: str) -> float:
    """A typed value, checked for the parameter's kind; ``ValueError`` says why not."""
    cleaned = (text or "").strip().replace(UNITS[kind], "").replace("deg", "").strip()
    if not cleaned:
        raise ValueError(f"Enter a value in {UNITS[kind]}.")
    try:
        value = float(cleaned.replace(",", "."))
    except ValueError:
        raise ValueError(f"'{text.strip()}' is not a number.") from None
    if not np.isfinite(value):
        raise ValueError("The value must be a finite number.")
    if kind == DISTANCE and not MIN_DISTANCE <= value <= MAX_DISTANCE:
        raise ValueError(f"A bond length must be between {MIN_DISTANCE} and {MAX_DISTANCE} Å.")
    if kind == ANGLE and not 0.0 < value <= 180.0:
        raise ValueError("A bond angle must be greater than 0° and at most 180°.")
    if kind == DIHEDRAL:
        value = wrap_degrees(value)
    return value


# ----------------------------------------------------------------------
# Values and analytic gradients
# ----------------------------------------------------------------------
def measure(kind: str, coordinates, atoms: Sequence[int]) -> float:
    """The parameter's value: Angstrom for a distance, degrees for the angles."""
    value, _ = value_and_gradient(kind, coordinates, atoms, gradient=False)
    return value


def value_and_gradient(kind: str, coordinates, atoms: Sequence[int], gradient: bool = True):
    """``(value, d value / d coordinates)``; the gradient is ``(N, 3)`` and is
    taken with respect to radians for the angles, the unit the constraint
    solver works in, while the value is reported in degrees."""
    coordinates = np.asarray(coordinates, dtype=float)
    positions = coordinates[list(atoms)]
    if kind == DISTANCE:
        value, parts = _distance(*positions)
    elif kind == ANGLE:
        value, parts = _angle(*positions)
    elif kind == DIHEDRAL:
        value, parts = _dihedral(*positions)
    else:
        raise ValueError(f"Unknown parameter kind: {kind}")
    if not gradient:
        return (value if kind == DISTANCE else float(np.degrees(value))), None
    full = np.zeros_like(coordinates)
    for atom, part in zip(atoms, parts):
        full[atom] += part
    return (value if kind == DISTANCE else float(np.degrees(value))), full


def _distance(p0, p1):
    separation = p0 - p1
    length = float(np.linalg.norm(separation))
    if length == 0.0:
        return 0.0, (np.zeros(3), np.zeros(3))
    unit = separation / length
    return length, (unit, -unit)


def _angle(p0, p1, p2):
    u, v = p0 - p1, p2 - p1
    norm_u, norm_v = float(np.linalg.norm(u)), float(np.linalg.norm(v))
    if norm_u == 0.0 or norm_v == 0.0:
        return 0.0, (np.zeros(3),) * 3
    cross = float(np.linalg.norm(np.cross(u, v)))
    dot = float(np.dot(u, v))
    theta = float(np.arctan2(cross, dot))
    sine = cross / (norm_u * norm_v)
    if sine < 1e-12:
        return theta, (np.zeros(3),) * 3  # linear: no defined direction
    cosine = dot / (norm_u * norm_v)
    d_u = -(v / (norm_u * norm_v) - cosine * u / norm_u ** 2) / sine
    d_v = -(u / (norm_u * norm_v) - cosine * v / norm_v ** 2) / sine
    return theta, (d_u, -(d_u + d_v), d_v)


def _dihedral(p0, p1, p2, p3):
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    n1_sq, n2_sq = float(np.dot(n1, n1)), float(np.dot(n2, n2))
    b2_length = float(np.linalg.norm(b2))
    if b2_length == 0.0 or n1_sq < 1e-20 or n2_sq < 1e-20:
        return 0.0, (np.zeros(3),) * 4
    phi = float(np.arctan2(np.dot(np.cross(n1, n2), b2 / b2_length), np.dot(n1, n2)))
    d0 = -(b2_length / n1_sq) * n1
    d3 = (b2_length / n2_sq) * n2
    p = float(np.dot(b1, b2)) / b2_length ** 2
    q = float(np.dot(b3, b2)) / b2_length ** 2
    d1 = -(1.0 + p) * d0 + q * d3
    d2 = p * d0 - (1.0 + q) * d3
    return phi, (d0, d1, d2, d3)


def is_degenerate(kind: str, coordinates, atoms: Sequence[int]) -> Optional[str]:
    """Why the parameter is undefined or cannot be held here, or ``None``."""
    coordinates = np.asarray(coordinates, dtype=float)
    if len(set(atoms)) != len(atoms):
        return "the same atom appears twice"
    if kind == ANGLE:
        theta = measure(ANGLE, coordinates, atoms)
        if theta < LINEAR_TOLERANCE_DEGREES or theta > 180.0 - LINEAR_TOLERANCE_DEGREES:
            return "the three atoms are (nearly) in a line, so the angle has no plane"
    if kind == DIHEDRAL:
        for triple in (atoms[:3], atoms[1:]):
            theta = measure(ANGLE, coordinates, triple)
            if theta > DIHEDRAL_LINEAR_LIMIT or theta < 180.0 - DIHEDRAL_LINEAR_LIMIT:
                return ("atoms " + "-".join(str(a + 1) for a in triple) +
                        " are (nearly) in a line, so the dihedral is undefined")
    return None


# ----------------------------------------------------------------------
# Holding frozen parameters
# ----------------------------------------------------------------------
def constraint_residuals(coordinates, constraints: Sequence[Constraint]):
    """``(residuals, jacobian)``: residual in Angstrom or radians, one row per constraint."""
    coordinates = np.asarray(coordinates, dtype=float).reshape(-1, 3)
    residuals = np.zeros(len(constraints))
    jacobian = np.zeros((len(constraints), coordinates.size))
    for row, constraint in enumerate(constraints):
        value, gradient = value_and_gradient(constraint.kind, coordinates, constraint.atoms)
        if constraint.kind == DISTANCE:
            residuals[row] = value - constraint.value
        else:
            difference = value - constraint.value
            if constraint.kind == DIHEDRAL:
                difference = wrap_degrees(difference)
            residuals[row] = np.radians(difference)
        jacobian[row] = gradient.ravel()
    return residuals, jacobian


def project_onto_constraints(coordinates, constraints: Sequence[Constraint],
                             tolerance: float = PROJECTION_TOLERANCE) -> np.ndarray:
    """The nearest geometry (Gauss-Newton, minimum-norm steps) holding every constraint.

    Raises ``ValueError`` when no such geometry is reached -- constraints
    that contradict each other, such as three bond lengths no triangle has.
    """
    x = np.asarray(coordinates, dtype=float).reshape(-1, 3).copy()
    if not constraints:
        return x
    residuals, jacobian = constraint_residuals(x, constraints)
    for _ in range(PROJECTION_MAX_ITERATIONS):
        error = float(np.max(np.abs(residuals)))
        if error < tolerance:
            return x
        step = -jacobian.T @ np.linalg.lstsq(jacobian @ jacobian.T, residuals, rcond=None)[0]
        scale = 1.0
        for _ in range(12):  # halve the step until the violation shrinks
            trial = x + scale * step.reshape(-1, 3)
            trial_residuals, trial_jacobian = constraint_residuals(trial, constraints)
            if np.max(np.abs(trial_residuals)) < error:
                break
            scale *= 0.5
        else:
            break
        x, residuals, jacobian = trial, trial_residuals, trial_jacobian
    raise ValueError(
        "The frozen parameters cannot all hold at once (largest violation "
        f"{float(np.max(np.abs(residuals))):.3g}). Remove or change one of them."
    )


def tangent_component(vector, coordinates, constraints: Sequence[Constraint]) -> np.ndarray:
    """``vector`` minus its part along the constraint normals: the directions
    in which a geometry can move without breaking a frozen parameter."""
    vector = np.asarray(vector, dtype=float).ravel()
    if not constraints:
        return vector
    _, jacobian = constraint_residuals(coordinates, constraints)
    coefficients = np.linalg.lstsq(jacobian @ jacobian.T, jacobian @ vector, rcond=None)[0]
    return vector - jacobian.T @ coefficients


# ----------------------------------------------------------------------
# Editing a parameter by moving one side
# ----------------------------------------------------------------------
def adjacency(atom_count: int, bonds: Iterable[Tuple[int, int]]) -> Dict[int, Set[int]]:
    neighbours: Dict[int, Set[int]] = {index: set() for index in range(atom_count)}
    for i, j in bonds:
        neighbours[i].add(j)
        neighbours[j].add(i)
    return neighbours


def side_of_bond(neighbours: Dict[int, Set[int]], fixed: int, moving: int) -> Optional[Set[int]]:
    """Atoms on ``moving``'s side once the fixed-moving bond is cut, or ``None``
    when the bond is in a ring (the two sides are still connected)."""
    side = {moving}
    queue = deque([moving])
    while queue:
        current = queue.popleft()
        for neighbour in neighbours[current]:
            if current == moving and neighbour == fixed:
                continue  # the cut bond
            if neighbour == fixed:
                return None
            if neighbour not in side:
                side.add(neighbour)
                queue.append(neighbour)
    return side


def fragment_of(neighbours: Dict[int, Set[int]], atom: int) -> Set[int]:
    seen = {atom}
    queue = deque([atom])
    while queue:
        for neighbour in neighbours[queue.popleft()]:
            if neighbour not in seen:
                seen.add(neighbour)
                queue.append(neighbour)
    return seen


def rotation_matrix(axis, angle_radians: float) -> np.ndarray:
    """Right-handed rotation about a unit ``axis`` (Rodrigues)."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    cross = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + np.sin(angle_radians) * cross + (1.0 - np.cos(angle_radians)) * cross @ cross


def _any_perpendicular(vector) -> np.ndarray:
    trial = np.array([0.0, 0.0, 1.0]) if abs(vector[2]) < 0.9 * np.linalg.norm(vector) else np.array([1.0, 0.0, 0.0])
    perpendicular = np.cross(vector, trial)
    return perpendicular / np.linalg.norm(perpendicular)


def moving_atoms(kind: str, atom_count: int, bonds, atoms: Sequence[int]) -> Set[int]:
    """The atoms an edit of this parameter moves; ``GeometryEditError`` if none can.

    The side of the *last* picked atom moves, cut at a bond: i-j for a length,
    j-k for an angle i-j-k, and the central bond j-k for a dihedral i-j-k-l.
    Every consecutive pair must be bonded, except a distance between two
    separate molecules, where the second molecule moves as a whole.
    """
    neighbours = adjacency(atom_count, bonds)
    names = "-".join(str(a + 1) for a in atoms)
    for first, second in zip(atoms, atoms[1:]):
        if second not in neighbours[first]:
            if kind == DISTANCE and second not in fragment_of(neighbours, first):
                return fragment_of(neighbours, second)
            raise GeometryEditError(
                f"Atoms {first + 1} and {second + 1} are not bonded. Only a {KIND_LABELS[kind].lower()} "
                f"along bonds ({names}) can be edited"
                + (", or the distance between two separate molecules." if kind == DISTANCE else ".")
            )
    pivot, moving = (atoms[0], atoms[1]) if kind == DISTANCE else (atoms[-3], atoms[-2]) \
        if kind == DIHEDRAL else (atoms[1], atoms[2])
    side = side_of_bond(neighbours, pivot, moving)
    if side is None:
        raise GeometryEditError(
            f"The bond {pivot + 1}-{moving + 1} is part of a ring, so no side of the molecule can "
            f"move on its own to change {names}. Ring geometry can only change by pre-optimizing."
        )
    if kind == DIHEDRAL:
        side.discard(moving)  # atom k is on the rotation axis: it does not move
    return side


def set_parameter(kind: str, coordinates, bonds, atoms: Sequence[int], value: float) -> np.ndarray:
    """New coordinates with the parameter at ``value``; one side moves rigidly.

    A length moves the last atom's side along the bond; an angle turns it
    about the vertex, in the plane of the three atoms; a dihedral turns it
    about the central bond. Nothing else moves.
    """
    coordinates = np.asarray(coordinates, dtype=float).reshape(-1, 3).copy()
    atoms = [int(a) for a in atoms]
    # A linear angle can still be bent (about any perpendicular axis); only a
    # dihedral through a linear angle, or a repeated atom, is undefined.
    if len(set(atoms)) != len(atoms):
        problem = "the same atom appears twice"
    else:
        problem = is_degenerate(kind, coordinates, atoms) if kind == DIHEDRAL else None
    if problem:
        raise GeometryEditError(f"The {KIND_LABELS[kind].lower()} cannot be set: {problem}.")
    side = sorted(moving_atoms(kind, len(coordinates), bonds, atoms))
    current = measure(kind, coordinates, atoms)

    if kind == DISTANCE:
        i, j = atoms
        direction = coordinates[j] - coordinates[i]
        length = np.linalg.norm(direction)
        if length == 0.0:
            raise GeometryEditError("The two atoms are on top of each other.")
        coordinates[side] += (value - current) * direction / length
        return coordinates

    if kind == ANGLE:
        i, j, k = atoms
        u, v = coordinates[i] - coordinates[j], coordinates[k] - coordinates[j]
        axis = np.cross(u, v)
        if np.linalg.norm(axis) < 1e-9 * np.linalg.norm(u) * np.linalg.norm(v):
            axis = _any_perpendicular(u)  # linear: any plane through the bond will do
        # Turning v about u x v by +a opens the angle by +a.
        rotation = rotation_matrix(axis, np.radians(value - current))
    else:
        j, k = atoms[1], atoms[2]
        axis = coordinates[k] - coordinates[j]
        # Turning the l side about j->k by +a changes the IUPAC dihedral by +a.
        rotation = rotation_matrix(axis, np.radians(wrap_degrees(value - current)))
    # Both rotations pass through atom j (the vertex, or the start of the central bond).
    centre = coordinates[atoms[1]].copy()
    coordinates[side] = (coordinates[side] - centre) @ rotation.T + centre
    return coordinates
