"""Data model for the interactive molecular design canvas.

This is the single source of truth for what the user has drawn: a flat graph of
:class:`SketchAtom` nodes at canvas pixel positions and :class:`SketchBond`
edges carrying a bond order. It owns the rules the Design tab depends on --
bond-order cycling, hit-testing and structure validation -- and deliberately
knows nothing about Toga, OpenGL or geometry optimization, so the whole model
is testable headless (see tests/test_molecularSketch.py).

Two conventions matter throughout:

* **Atoms are addressed by a stable integer id, never by list position.**
  Deleting an atom would otherwise renumber every bond, which is exactly the
  class of silent off-by-one this project has been bitten by before. Ids are
  handed out in increasing order and never reused.
* **Positions are sketch coordinates**, not Angstrom: canvas pixels at the
  default view (x right, y down) plus a depth ``z`` that grows *away* from the
  viewer, in the same pixel units. Every atom is a 3D point; a structure drawn
  without rotating simply has ``z = 0`` everywhere. How those points reach the
  screen is molecularCamera's job, and the conversion to Angstrom belongs to
  molecularPreOptimizer, the only module that needs to know the scale.
"""
from dataclasses import dataclass, replace
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from help import AtomicData

# Where an atom appears on screen: ``(screen_x, screen_y, depth)``, with depth
# growing away from the viewer. Hit-testing takes a mapping of atom id to one of
# these, so it works the same on a rotated view as on the flat default one.
ScreenPositions = Mapping[int, Tuple[float, float, float]]

# Bond order used for a newly drawn bond, and the ceiling the cycle wraps at
# for a pair of elements that both support it.
SINGLE_BOND = 1
DEFAULT_MAX_BOND_ORDER = 3

# A delocalized (resonance) bond genuinely *is* order 1.5, so it is stored as a
# number rather than as a separate flag. That is what lets valence sums,
# reference bond lengths and the planarity torsion all keep working
# arithmetically, with no parallel code path for "is this aromatic".
RESONANCE_ORDER = 1.5

BOND_TYPE_NAMES = {
    1.0: "single",
    RESONANCE_ORDER: "resonance",
    2.0: "double",
    3.0: "triple",
}


def bond_type_name(order: float) -> str:
    """Human-readable name of a bond order, for status lines and summaries."""
    return BOND_TYPE_NAMES.get(float(order), f"order {format_order(order)}")


def format_order(value: float) -> str:
    """Render a bond order or valence sum without a pointless trailing zero."""
    number = float(value)
    if number.is_integer():
        return str(int(number))
    return f"{number:g}"


@dataclass
class SketchAtom:
    """One drawn atom in sketch coordinates; ``atom_id`` is stable."""

    atom_id: int
    element: str
    x: float
    y: float
    z: float = 0.0


@dataclass
class SketchBond:
    """One drawn bond between two atom ids.

    ``order`` is 1, 2 or 3 for a normal bond and :data:`RESONANCE_ORDER` (1.5)
    for a delocalized one.
    """

    atom_i: int
    atom_j: int
    order: float = SINGLE_BOND

    @property
    def is_resonance(self) -> bool:
        return self.order == RESONANCE_ORDER

    @property
    def type_name(self) -> str:
        return bond_type_name(self.order)

    def connects(self, atom_id: int) -> bool:
        return atom_id in (self.atom_i, self.atom_j)

    def other_end(self, atom_id: int) -> int:
        return self.atom_j if atom_id == self.atom_i else self.atom_i

    def same_pair(self, atom_i: int, atom_j: int) -> bool:
        return {self.atom_i, self.atom_j} == {atom_i, atom_j}


# Frozen-parameter kinds, keyed by how many atoms define them.
CONSTRAINT_KINDS = {2: "distance", 3: "angle", 4: "dihedral"}


@dataclass
class SketchConstraint:
    """A bond length, angle or dihedral frozen during pre-optimization.

    ``atom_ids`` are stable sketch ids in pick order; ``value`` is in Angstrom
    for a distance and degrees for the angles -- physical units, because the
    sketch's own pixel scale is not the molecule's.
    """

    kind: str
    atom_ids: Tuple[int, ...]
    value: float

    def same_parameter(self, atom_ids: Sequence[int]) -> bool:
        """True for the same atoms in the same or exactly reversed order."""
        ids = tuple(atom_ids)
        return ids in (self.atom_ids, tuple(reversed(self.atom_ids)))


def normalise_element(symbol: str) -> str:
    """Return ``symbol`` in standard capitalisation ("cL" -> "Cl").

    Raises ``ValueError`` when nothing is left after stripping, so a blank
    element can never reach the model.
    """
    cleaned = (symbol or "").strip()
    if not cleaned:
        raise ValueError("An element symbol is required to place an atom.")
    return cleaned[0].upper() + cleaned[1:].lower()


def point_to_segment_distance(
    px: float, py: float, x1: float, y1: float, x2: float, y2: float
) -> float:
    """Distance from a point to a finite segment (not to its infinite line).

    A point beyond either end measures to that end, which is what keeps a click
    far past a short bond from selecting it.
    """
    dx = x2 - x1
    dy = y2 - y1
    length_squared = dx * dx + dy * dy
    if length_squared == 0.0:
        return float(((px - x1) ** 2 + (py - y1) ** 2) ** 0.5)
    t = ((px - x1) * dx + (py - y1) * dy) / length_squared
    t = max(0.0, min(1.0, t))
    nearest_x = x1 + t * dx
    nearest_y = y1 + t * dy
    return float(((px - nearest_x) ** 2 + (py - nearest_y) ** 2) ** 0.5)


class MoleculeSketch:
    """The atoms and bonds the user has drawn on the design canvas."""

    # Click tolerances in screen pixels. An atom is drawn at ATOM_DISC_RADIUS,
    # so the atom target is slightly generous and the bond target deliberately
    # tighter, since the caller tests atoms first.
    ATOM_DISC_RADIUS = 11.0
    ATOM_HIT_RADIUS = 14.0
    BOND_HIT_TOLERANCE = 6.0

    def __init__(self):
        self._atoms: Dict[int, SketchAtom] = {}
        self._bonds: List[SketchBond] = []
        self._constraints: List[SketchConstraint] = []
        self._next_atom_id = 1

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------
    @property
    def atoms(self) -> List[SketchAtom]:
        """Atoms in insertion order; this order is the exported atom order."""
        return list(self._atoms.values())

    @property
    def bonds(self) -> List[SketchBond]:
        return list(self._bonds)

    @property
    def atom_count(self) -> int:
        return len(self._atoms)

    @property
    def bond_count(self) -> int:
        return len(self._bonds)

    @property
    def is_empty(self) -> bool:
        return not self._atoms

    def get_atom(self, atom_id: int) -> SketchAtom:
        try:
            return self._atoms[atom_id]
        except KeyError:
            raise ValueError(f"Unknown atom id {atom_id}.") from None

    def ordered_atom_ids(self) -> List[int]:
        """Atom ids in insertion order -- the single id-to-index mapping."""
        return list(self._atoms.keys())

    def bond_index_pairs(self) -> List[Tuple[int, int]]:
        """Bonds as 0-based ``(i, j)`` positions, low index first.

        This is the shape the viewer wants: it indexes atoms by position in the
        frame, so the stable ids used inside the sketch have to be mapped out
        exactly once, here.
        """
        index_of = {atom_id: index for index, atom_id in enumerate(self.ordered_atom_ids())}
        pairs = []
        for bond in self._bonds:
            i = index_of.get(bond.atom_i)
            j = index_of.get(bond.atom_j)
            if i is None or j is None:
                continue
            pairs.append((min(i, j), max(i, j)))
        return pairs

    def positions(self) -> np.ndarray:
        """``(N, 3)`` sketch coordinates in :meth:`ordered_atom_ids` order."""
        return np.array(
            [[atom.x, atom.y, atom.z] for atom in self._atoms.values()], dtype=float
        ).reshape(-1, 3)

    def update_positions(self, positions: Sequence[Sequence[float]]):
        """Move every atom at once, in :meth:`ordered_atom_ids` order.

        This is how an optimized geometry is written back onto the drawing.
        Ids, elements and bonds are untouched; a wrong shape or a non-finite
        value is refused before anything moves.
        """
        array = np.asarray(positions, dtype=float)
        if array.shape != (self.atom_count, 3):
            raise ValueError(
                f"Expected positions for {self.atom_count} atoms as (x, y, z) rows, "
                f"got an array of shape {array.shape}."
            )
        if not np.all(np.isfinite(array)):
            raise ValueError("Atom positions must be finite numbers.")
        for atom, (x, y, z) in zip(self._atoms.values(), array):
            atom.x, atom.y, atom.z = float(x), float(y), float(z)

    def atom_number(self, atom_id: int) -> int:
        """1-based position of an atom, matching how the rest of gQTEA labels atoms."""
        return self.ordered_atom_ids().index(atom_id) + 1

    def describe_atom(self, atom_id: int) -> str:
        atom = self.get_atom(atom_id)
        return f"atom {self.atom_number(atom_id)} ({atom.element})"

    # ------------------------------------------------------------------
    # Editing
    # ------------------------------------------------------------------
    def add_atom(self, element: str, x: float, y: float, z: float = 0.0) -> SketchAtom:
        atom = SketchAtom(
            atom_id=self._next_atom_id,
            element=normalise_element(element),
            x=float(x),
            y=float(y),
            z=float(z),
        )
        self._atoms[atom.atom_id] = atom
        self._next_atom_id += 1
        return atom

    def replace_element(self, atom_id: int, element: str) -> str:
        """Change an atom's element in place; returns the old element.

        The atom keeps its id, position, bonds (orders included) and frozen
        parameters -- only what it *is* changes, so picks, the selection and
        constraints that name it stay valid. Whether the result is chemically
        usual is :meth:`validate`'s business, as for any drawn atom.
        """
        atom = self.get_atom(atom_id)
        previous = atom.element
        atom.element = normalise_element(element)
        return previous

    def remove_atom(self, atom_id: int) -> bool:
        """Remove an atom, every bond touching it and every frozen parameter
        that names it. False if it is unknown."""
        if atom_id not in self._atoms:
            return False
        del self._atoms[atom_id]
        self._bonds = [bond for bond in self._bonds if not bond.connects(atom_id)]
        self._constraints = [c for c in self._constraints if atom_id not in c.atom_ids]
        return True

    def add_bond(self, atom_i: int, atom_j: int, order: int = SINGLE_BOND) -> Optional[SketchBond]:
        """Bond two atoms.

        Returns the new bond, or ``None`` when the pair is already bonded -- the
        caller reports that rather than silently changing an existing order,
        because raising an order is the job of :meth:`cycle_bond_order`.
        """
        if atom_i == atom_j:
            raise ValueError("An atom cannot be bonded to itself.")
        for atom_id in (atom_i, atom_j):
            if atom_id not in self._atoms:
                raise ValueError(f"Unknown atom id {atom_id}.")
        if self.find_bond(atom_i, atom_j) is not None:
            return None
        bond = SketchBond(atom_i=atom_i, atom_j=atom_j, order=float(order))
        self._bonds.append(bond)
        return bond

    def find_bond(self, atom_i: int, atom_j: int) -> Optional[SketchBond]:
        for bond in self._bonds:
            if bond.same_pair(atom_i, atom_j):
                return bond
        return None

    def remove_bond(self, atom_i: int, atom_j: int) -> bool:
        bond = self.find_bond(atom_i, atom_j)
        if bond is None:
            return False
        self._bonds.remove(bond)
        return True

    def remove_bond_object(self, bond: SketchBond) -> bool:
        if bond not in self._bonds:
            return False
        self._bonds.remove(bond)
        return True

    def clear(self):
        self._atoms.clear()
        self._bonds.clear()
        self._constraints.clear()

    def replace_contents(
        self,
        elements: Sequence[str],
        positions: Sequence[Sequence[float]],
        bonds: Sequence[Tuple[int, int, float]] = (),
    ):
        """Replace the whole drawing with a structure read from a file.

        ``positions`` are sketch coordinates, one ``(x, y, z)`` row per
        element; ``bonds`` are ``(i, j, order)`` with **0-based** positions
        into ``elements``. Everything is checked before the current drawing is
        touched, so a bad structure leaves the sketch exactly as it was. New
        atoms get fresh ids -- ids are never reused, even across a reload.
        """
        checked = self._checked_structure(elements, positions, bonds)
        self.clear()
        self._append_structure(*checked)

    def add_structure(
        self,
        elements: Sequence[str],
        positions: Sequence[Sequence[float]],
        bonds: Sequence[Tuple[int, int, float]] = (),
    ) -> List[int]:
        """Add a structure beside what is already drawn; returns the new atom ids.

        Same arguments and the same all-or-nothing checks as
        :meth:`replace_contents`; bonds index the *new* atoms only.
        """
        return self._append_structure(*self._checked_structure(elements, positions, bonds))

    def _append_structure(self, symbols, array, checked_bonds) -> List[int]:
        ids = [self.add_atom(symbol, x, y, z).atom_id for symbol, (x, y, z) in zip(symbols, array)]
        # Appended directly: the pairs are already known to be unique, and
        # add_bond's linear duplicate search would make a large import quadratic.
        self._bonds.extend(SketchBond(ids[i], ids[j], order) for i, j, order in checked_bonds)
        return ids

    def _checked_structure(self, elements, positions, bonds):
        """Validate a structure for :meth:`replace_contents` / :meth:`add_structure`."""
        symbols = [normalise_element(element) for element in elements]
        array = np.asarray(positions, dtype=float)
        if array.shape != (len(symbols), 3):
            raise ValueError(
                f"Expected positions for {len(symbols)} atoms as (x, y, z) rows, "
                f"got an array of shape {array.shape}."
            )
        if not np.all(np.isfinite(array)):
            raise ValueError("Atom positions must be finite numbers.")

        allowed_orders = set(BOND_TYPE_NAMES)
        seen = set()
        checked_bonds = []
        for i, j, order in bonds:
            i, j, order = int(i), int(j), float(order)
            if not (0 <= i < len(symbols) and 0 <= j < len(symbols)):
                raise ValueError(
                    f"A bond refers to atom {max(i, j) + 1}, but the structure has "
                    f"{len(symbols)} atoms."
                )
            if i == j:
                raise ValueError(f"Atom {i + 1} cannot be bonded to itself.")
            if order not in allowed_orders:
                raise ValueError(
                    f"Bond {i + 1}-{j + 1} has order {format_order(order)}; the canvas "
                    "supports single, double, triple and resonance (1.5) bonds."
                )
            pair = (min(i, j), max(i, j))
            if pair in seen:
                raise ValueError(f"Atoms {pair[0] + 1} and {pair[1] + 1} are bonded twice.")
            seen.add(pair)
            checked_bonds.append((i, j, order))
        return symbols, array, checked_bonds

    # ------------------------------------------------------------------
    # Frozen parameters
    # ------------------------------------------------------------------
    @property
    def constraints(self) -> List[SketchConstraint]:
        return list(self._constraints)

    def find_constraint(self, atom_ids: Sequence[int]) -> Optional[SketchConstraint]:
        for constraint in self._constraints:
            if constraint.same_parameter(atom_ids):
                return constraint
        return None

    def set_constraint(self, atom_ids: Sequence[int], value: float) -> SketchConstraint:
        """Freeze the parameter on these atoms at ``value`` (replacing an existing one)."""
        ids = tuple(int(atom_id) for atom_id in atom_ids)
        kind = CONSTRAINT_KINDS.get(len(ids))
        if kind is None:
            raise ValueError("A frozen parameter needs 2, 3 or 4 atoms.")
        if len(set(ids)) != len(ids):
            raise ValueError("A frozen parameter cannot use the same atom twice.")
        for atom_id in ids:
            self.get_atom(atom_id)
        if not np.isfinite(float(value)):
            raise ValueError("A frozen value must be a finite number.")
        existing = self.find_constraint(ids)
        if existing is not None:
            self._constraints.remove(existing)
        constraint = SketchConstraint(kind, ids, float(value))
        self._constraints.append(constraint)
        return constraint

    def remove_constraint(self, constraint: SketchConstraint) -> bool:
        if constraint not in self._constraints:
            return False
        self._constraints.remove(constraint)
        return True

    def clear_constraints(self):
        self._constraints.clear()

    def describe_parameter(self, atom_ids: Sequence[int]) -> str:
        """'1-2 (C-O)' style label of a parameter, with 1-based atom numbers."""
        numbers = "-".join(str(self.atom_number(atom_id)) for atom_id in atom_ids)
        elements = "-".join(self.get_atom(atom_id).element for atom_id in atom_ids)
        return f"{numbers} ({elements})"

    def fragment_containing(self, atom_id: int) -> List[int]:
        """The ids of the molecule (connected fragment) an atom belongs to."""
        self.get_atom(atom_id)
        for component in self.fragments():
            if atom_id in component:
                return component
        return [atom_id]

    # ------------------------------------------------------------------
    # Bond order
    # ------------------------------------------------------------------
    def max_bond_order_for(self, element_i: str, element_j: str) -> int:
        """Highest order a bond between these two elements may reach.

        The cap is the lower of the two elements' capabilities, so O-H stays
        single while C=O reaches double. An element with no entry is assumed to
        support a triple bond.
        """
        limits = [
            AtomicData.max_bond_order.get(normalise_element(element), DEFAULT_MAX_BOND_ORDER)
            for element in (element_i, element_j)
        ]
        return max(SINGLE_BOND, min(limits))

    def max_bond_order_of(self, bond: SketchBond) -> int:
        return self.max_bond_order_for(
            self.get_atom(bond.atom_i).element, self.get_atom(bond.atom_j).element
        )

    def bond_order_cycle(self, maximum: int) -> Tuple[float, ...]:
        """The orders a bond with this cap steps through, in click order.

        Resonance sits last, and only where a double bond is possible: a pair
        capped at single (O-H) never reaches it, so the cycle cannot produce a
        delocalized bond that makes no chemical sense.
        """
        orders: List[float] = [float(order) for order in range(1, int(maximum) + 1)]
        if maximum >= 2:
            orders.append(RESONANCE_ORDER)
        return tuple(orders)

    def cycle_bond_order(self, bond: SketchBond) -> float:
        """Advance a bond one step around its cycle.

        Single -> double -> triple -> resonance -> single for a pair that
        supports them all, with triple dropped when the elements cap at two.
        Clicking a bond that is already at the end of its cycle returns it to a
        single bond rather than doing nothing, so one gesture both raises and
        lowers the order and no state is unreachable.
        """
        cycle = self.bond_order_cycle(self.max_bond_order_of(bond))
        try:
            position = cycle.index(float(bond.order))
        except ValueError:
            # An order outside the cycle (a narrowed cap, a hand-set value):
            # step back to single rather than stranding the bond.
            bond.order = float(SINGLE_BOND)
            return bond.order
        bond.order = cycle[(position + 1) % len(cycle)]
        return bond.order

    def bond_order_sum(self, atom_id: int) -> float:
        """Total bond order on an atom -- a double counts twice, a resonance 1.5."""
        return float(sum(bond.order for bond in self._bonds if bond.connects(atom_id)))

    # ------------------------------------------------------------------
    # Hit-testing
    # ------------------------------------------------------------------
    def _screen_position(
        self, atom: SketchAtom, screen: Optional[ScreenPositions]
    ) -> Optional[Tuple[float, float, float]]:
        """Where ``atom`` is on screen; its own x/y/z on the unrotated view."""
        if screen is None:
            return atom.x, atom.y, atom.z
        return screen.get(atom.atom_id)

    def atom_at(
        self, x: float, y: float, screen: Optional[ScreenPositions] = None
    ) -> Optional[SketchAtom]:
        """The atom a click at screen point ``(x, y)`` selects, or ``None``.

        ``screen`` maps atom ids to projected ``(sx, sy, depth)``; without it
        the atoms' own x/y are the screen positions (the default view). A click
        inside an atom's drawn disc beats one merely within the hit tolerance,
        and among equals the atom in front wins, then the nearest -- so after a
        rotation the atom the user can see is the one that gets picked.
        """
        best: Optional[SketchAtom] = None
        best_key = None
        for atom in self._atoms.values():
            position = self._screen_position(atom, screen)
            if position is None:
                continue
            sx, sy, depth = position
            distance = ((sx - x) ** 2 + (sy - y) ** 2) ** 0.5
            if distance > self.ATOM_HIT_RADIUS:
                continue
            key = (distance > self.ATOM_DISC_RADIUS, depth, distance)
            if best_key is None or key < best_key:
                best = atom
                best_key = key
        return best

    def bond_at(
        self, x: float, y: float, screen: Optional[ScreenPositions] = None
    ) -> Optional[SketchBond]:
        """Nearest bond within ``BOND_HIT_TOLERANCE`` of a screen point."""
        best: Optional[SketchBond] = None
        best_distance = self.BOND_HIT_TOLERANCE
        for bond in self._bonds:
            start = self._atoms.get(bond.atom_i)
            end = self._atoms.get(bond.atom_j)
            if start is None or end is None:
                continue
            start_position = self._screen_position(start, screen)
            end_position = self._screen_position(end, screen)
            if start_position is None or end_position is None:
                continue
            distance = point_to_segment_distance(
                x, y, start_position[0], start_position[1], end_position[0], end_position[1]
            )
            if distance <= best_distance:
                best = bond
                best_distance = distance
        return best

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    def neighbours(self, atom_id: int) -> List[int]:
        return [bond.other_end(atom_id) for bond in self._bonds if bond.connects(atom_id)]

    def _has_resonance_neighbour(self, bond: SketchBond) -> bool:
        """Whether another resonance bond shares either end of this one.

        One shared end is enough: a carboxylate or an amide is delocalized
        without being a ring, and only a resonance bond with no resonance
        neighbour at all is meaningless.
        """
        for other in self._bonds:
            if other is bond or not other.is_resonance:
                continue
            if other.connects(bond.atom_i) or other.connects(bond.atom_j):
                return True
        return False

    def fragments(self) -> List[List[int]]:
        """Connected components, each a list of atom ids in insertion order."""
        unvisited = set(self._atoms)
        components: List[List[int]] = []
        for atom_id in self.ordered_atom_ids():
            if atom_id not in unvisited:
                continue
            component: List[int] = []
            stack = [atom_id]
            unvisited.discard(atom_id)
            while stack:
                current = stack.pop()
                component.append(current)
                for neighbour in self.neighbours(current):
                    if neighbour in unvisited:
                        unvisited.discard(neighbour)
                        stack.append(neighbour)
            components.append(component)
        return components

    def validate(self) -> Tuple[List[str], List[str]]:
        """Check the drawn structure.

        Returns ``(errors, warnings)``. Errors block pre-optimization because
        there is no geometry to produce; warnings are reported and the run goes
        ahead, because unusual valences and multi-fragment systems are
        chemically real and the tool does not police what the user draws.
        """
        errors: List[str] = []
        warnings: List[str] = []

        if not self._atoms:
            errors.append("The sketch has no atoms.")
            return errors, warnings

        for atom in self._atoms.values():
            if atom.element not in AtomicData.covalent_radii:
                errors.append(
                    f"No covalent radius is known for {atom.element} "
                    f"({self.describe_atom(atom.atom_id)}), so it cannot be optimized."
                )

        # A lone atom on its own is a perfectly good structure; a lone atom
        # beside a molecule is almost always a stray click.
        if self.atom_count > 1:
            for atom in self._atoms.values():
                if not self.neighbours(atom.atom_id):
                    errors.append(
                        f"{self.describe_atom(atom.atom_id).capitalize()} is not bonded "
                        "to anything. Bond it or delete it before optimizing."
                    )

        for atom in self._atoms.values():
            limit = AtomicData.max_valence.get(atom.element)
            if limit is None:
                continue
            total = self.bond_order_sum(atom.atom_id)
            if total > limit:
                warnings.append(
                    f"Unusual valence on {self.describe_atom(atom.atom_id)}: "
                    f"{format_order(total)} bonds where {limit} is typical."
                )

        for bond in self._bonds:
            if not bond.is_resonance or self._has_resonance_neighbour(bond):
                continue
            warnings.append(
                "Isolated resonance bond between "
                f"{self.describe_atom(bond.atom_i)} and {self.describe_atom(bond.atom_j)}: "
                "delocalization over a single bond has no meaning. Give it a resonance "
                "neighbour, or make it a single or double bond."
            )

        components = self.fragments()
        if len(components) > 1:
            warnings.append(
                f"The sketch contains {len(components)} disconnected fragments; "
                "they will be optimized together in one structure."
            )

        return errors, warnings

    # ------------------------------------------------------------------
    # Undo support
    # ------------------------------------------------------------------
    def snapshot(self) -> Tuple[List[SketchAtom], List[SketchBond], int, List[SketchConstraint]]:
        """A deep-enough copy of the current state for :meth:`restore`.

        Frozen parameters are part of it, so Undo brings them back with the
        atoms they refer to.
        """
        return (
            [replace(atom) for atom in self._atoms.values()],
            [replace(bond) for bond in self._bonds],
            self._next_atom_id,
            [replace(constraint) for constraint in self._constraints],
        )

    def restore(self, state):
        atoms, bonds, next_atom_id = state[:3]
        constraints = state[3] if len(state) > 3 else []
        self._atoms = {atom.atom_id: replace(atom) for atom in atoms}
        self._bonds = [replace(bond) for bond in bonds]
        self._next_atom_id = next_atom_id
        self._constraints = [replace(constraint) for constraint in constraints]

    def copy(self) -> "MoleculeSketch":
        """An independent duplicate, e.g. to optimize while the user keeps drawing."""
        duplicate = MoleculeSketch()
        duplicate.restore(self.snapshot())
        return duplicate
