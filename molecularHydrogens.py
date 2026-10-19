"""Add the hydrogens a drawn molecule implies, or say why not.

:func:`plan_hydrogens` reads elements, Cartesian coordinates (Angstrom) and
bonds with orders, and returns the hydrogens to add -- each with the atom it
belongs to and a position -- plus the atoms it refused to complete and why.
The design canvas applies the plan as one undo step. Toga-free, tested
headless in tests/test_molecularHydrogens.py.

**How many.** An atom gets ``valence - bond order sum`` hydrogens, from
:data:`HYDROGEN_VALENCES`. Existing hydrogens are bonds like any other, so a
complete atom gets none and running the tool twice adds nothing. Nothing is
guessed where the drawing does not say:

* more bonds than the element's valence (ammonium N, a 5-bonded C) -- likely
  a charge the canvas does not record: flagged, untouched;
* an element with several valences and a bond sum between them (S with 3
  bonds, P with 4): flagged;
* a resonance-bonded atom whose sum is fractional and below its valence
  (a carboxylate O, 1.5): H or charge cannot be told apart: flagged;
* an aromatic N or P with exactly two resonance bonds: pyridine-like (no H)
  or pyrrole-like (N-H) cannot be told apart: flagged, no H added;
* an element without a valence rule (metals, noble gases): flagged.

A resonance-bonded atom whose sum is *above* its valence by less than one
(a ring-fusion carbon, 1.5 x 3 = 4.5) is simply complete: delocalization
over-counts, it is not a charge.

**Where.** VSEPR: the steric number is neighbours + new hydrogens + lone
pairs (valence electrons from help.py), and the new hydrogens take the free
directions of the linear, trigonal or tetrahedral arrangement around the
existing bonds, at the pre-optimizer's reference X-H length. A trigonal
centre is kept in the plane of its neighbour's substituents (so CH2=CH- comes
out planar) and the first tetrahedral position is anti to it (a zigzag). The
positions are a good start, not an optimized geometry: pre-optimize after.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from help import AtomicData
from molecularPreOptimizer import reference_bond_length

# Neutral valences for adding hydrogens; several values = hypervalent states.
HYDROGEN_VALENCES: Dict[str, Tuple[int, ...]] = {
    "B": (3,), "C": (4,), "N": (3,), "O": (2,), "F": (1,),
    "Al": (3,), "Si": (4,), "P": (3, 5), "S": (2, 4, 6), "Cl": (1,),
    "Ge": (4,), "As": (3, 5), "Se": (2, 4, 6), "Br": (1,),
    "Sn": (4,), "Sb": (3, 5), "Te": (2, 4, 6), "I": (1,),
}
DEFAULT_X_H_LENGTH = 1.09
TETRAHEDRAL = np.radians(109.4712206)
ORDER_TOLERANCE = 1e-6
RESONANCE = 1.5


@dataclass
class HydrogenPlan:
    """Hydrogens to add as ``(parent index, position)`` and the atoms left alone."""

    additions: List[Tuple[int, np.ndarray]] = field(default_factory=list)
    issues: List[str] = field(default_factory=list)

    @property
    def count(self) -> int:
        return len(self.additions)


def hydrogens_needed(element: str, orders: Sequence[float]) -> Tuple[int, Optional[str]]:
    """``(number of H, problem)``; a problem means none are added to this atom."""
    valences = HYDROGEN_VALENCES.get(element)
    if element == "H":
        return 0, None
    if valences is None:
        return 0, f"no valence rule for {element}"
    total = float(sum(orders))
    resonance = sum(1 for order in orders if abs(order - RESONANCE) < ORDER_TOLERANCE)
    fractional = abs(total - round(total)) > ORDER_TOLERANCE

    if resonance:
        missing = valences[0] - total
        if missing <= ORDER_TOLERANCE:
            if element in ("N", "P") and resonance == 2 and len(orders) == 2 and abs(missing) < ORDER_TOLERANCE:
                return 0, ("aromatic with two resonance bonds: pyridine-like (no H) or pyrrole-like "
                           "(N-H) cannot be told apart; add the H yourself if it is N-H")
            if missing > -1.0 + ORDER_TOLERANCE:
                return 0, None  # delocalization over-counts; the atom is complete
            return 0, f"{total:g} bonds exceed the valence of {valences[0]}"
        if fractional:
            return 0, (f"{total:g} bonds over resonance: whether it carries an H or a charge "
                       "cannot be told")
        return int(round(missing)), None

    for valence in valences:
        if total <= valence + ORDER_TOLERANCE:
            missing = valence - total
            if valence != valences[0] and missing > ORDER_TOLERANCE:
                return 0, (f"{total:g} bonds lie between the valences "
                           f"{', '.join(str(v) for v in valences)}: ambiguous")
            return int(round(missing)), None
    return 0, (f"{total:g} bonds exceed the usual valence of "
               f"{'/'.join(str(v) for v in valences)}: a charge is likely")


def plan_hydrogens(elements: Sequence[str], coordinates, bonds: Sequence[Tuple[int, int, float]],
                   scope: Optional[Set[int]] = None) -> HydrogenPlan:
    """Which hydrogens to add, and where, for the atoms in ``scope`` (default all)."""
    coordinates = np.asarray(coordinates, dtype=float).reshape(-1, 3)
    neighbours: List[List[Tuple[int, float]]] = [[] for _ in elements]
    for i, j, order in bonds:
        neighbours[i].append((j, float(order)))
        neighbours[j].append((i, float(order)))

    plan = HydrogenPlan()
    for atom, element in enumerate(elements):
        if scope is not None and atom not in scope:
            continue
        orders = [order for _, order in neighbours[atom]]
        count, problem = hydrogens_needed(element, orders)
        if problem:
            plan.issues.append(f"atom {atom + 1} ({element}): {problem}")
            continue
        if count <= 0:
            continue
        directions = _free_directions(atom, element, count, sum(orders), coordinates, neighbours)
        length = reference_bond_length(element, "H") or DEFAULT_X_H_LENGTH
        for direction in directions[:count]:
            plan.additions.append((atom, coordinates[atom] + length * direction))
    return plan


def _unit(vector) -> Optional[np.ndarray]:
    norm = float(np.linalg.norm(vector))
    return None if norm < 1e-8 else np.asarray(vector, dtype=float) / norm


def _perpendicular(vector, preferred=None) -> np.ndarray:
    """A unit vector perpendicular to ``vector``, as close to ``preferred`` as possible."""
    for candidate in ([preferred] if preferred is not None else []) + [
            np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])]:
        projected = _unit(candidate - np.dot(candidate, vector) * vector)
        if projected is not None:
            return projected
    raise ValueError("no perpendicular")  # unreachable for a unit vector


def _free_directions(atom, element, count, order_sum, coordinates, neighbours) -> List[np.ndarray]:
    """Unit vectors from ``atom`` for its new hydrogens, best first."""
    valence_electrons = AtomicData.valence_electrons.get(element, 0)
    lone_pairs = max(0, int((valence_electrons - (order_sum + count)) // 2))
    bonded = [j for j, _ in neighbours[atom]]
    steric = min(4, max(2, len(bonded) + count + lone_pairs))
    vectors = [v for v in (_unit(coordinates[j] - coordinates[atom]) for j in bonded) if v is not None]

    # Orientation reference: a substituent of the first neighbour, so the new
    # atoms come out planar (trigonal) or anti (tetrahedral) to it; failing
    # that, the screen plane (x-y in Angstrom).
    reference = None
    if bonded:
        others = [k for k, _ in neighbours[bonded[0]] if k != atom]
        if others:
            reference = _unit(coordinates[others[0]] - coordinates[bonded[0]])

    if not vectors:
        if steric == 2:
            return [np.array([1.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])]
        if steric == 3:
            return [np.array([np.cos(a), np.sin(a), 0.0]) for a in np.radians([90.0, 210.0, 330.0])]
        half = TETRAHEDRAL / 2.0
        return [np.array([np.sin(half), -np.cos(half), 0.0]), np.array([-np.sin(half), -np.cos(half), 0.0]),
                np.array([0.0, np.cos(half), np.sin(half)]), np.array([0.0, np.cos(half), -np.sin(half)])]

    b1 = vectors[0]
    if len(vectors) == 1:
        if steric == 2:
            return [-b1]
        p = _perpendicular(b1, reference if reference is not None else np.array([0.0, 1.0, 0.0]))
        if steric == 3:
            angle = np.radians(120.0)
            return [np.cos(angle) * b1 - np.sin(angle) * p, np.cos(angle) * b1 + np.sin(angle) * p]
        q = np.cross(b1, p)
        result = []
        for turn in np.radians([180.0, 60.0, -60.0]):  # anti to the reference first
            radial = np.cos(turn) * p + np.sin(turn) * q
            result.append(np.cos(TETRAHEDRAL) * b1 + np.sin(TETRAHEDRAL) * radial)
        return result

    total = _unit(-np.sum(vectors, axis=0))
    if len(vectors) == 2:
        normal = _unit(np.cross(vectors[0], vectors[1]))
        if total is None:  # the two bonds are opposite (linear)
            total = _perpendicular(vectors[0])
        if normal is None:
            normal = _unit(np.cross(vectors[0], total))
        if steric <= 3:
            return [total]
        half = TETRAHEDRAL / 2.0
        return [np.cos(half) * total + np.sin(half) * normal, np.cos(half) * total - np.sin(half) * normal]

    if total is None:  # three coplanar bonds: out of their plane
        total = _unit(np.cross(vectors[0], vectors[1]))
        if total is None:
            total = _perpendicular(vectors[0])
    return [total]
