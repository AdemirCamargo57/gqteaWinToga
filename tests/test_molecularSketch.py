"""Tests for the molecularSketch data model.

Pure model tests: no Toga, no OpenGL, no window. The sketch is the single
source of truth for what the user has drawn, so every rule the Design tab
relies on (stable ids, bond-order cycling, hit-testing, validation) is
pinned here.
"""
import pytest

from molecularSketch import (
    RESONANCE_ORDER,
    MoleculeSketch,
    SketchAtom,
    SketchBond,
    bond_type_name,
    format_order,
)


def _water() -> MoleculeSketch:
    """O with two bonded hydrogens, in canvas pixel coordinates."""
    sketch = MoleculeSketch()
    oxygen = sketch.add_atom("O", 100.0, 100.0)
    h1 = sketch.add_atom("H", 60.0, 140.0)
    h2 = sketch.add_atom("H", 140.0, 140.0)
    sketch.add_bond(oxygen.atom_id, h1.atom_id)
    sketch.add_bond(oxygen.atom_id, h2.atom_id)
    return sketch


# ------------------------------------------------------------------
# Atoms
# ------------------------------------------------------------------
def test_add_atom_returns_atom_with_element_and_position():
    sketch = MoleculeSketch()

    atom = sketch.add_atom("C", 10.0, 20.0)

    assert isinstance(atom, SketchAtom)
    assert atom.element == "C"
    assert (atom.x, atom.y) == (10.0, 20.0)
    assert sketch.atom_count == 1


def test_element_symbols_are_normalised_to_standard_capitalisation():
    sketch = MoleculeSketch()

    assert sketch.add_atom("c", 0.0, 0.0).element == "C"
    assert sketch.add_atom("cL", 0.0, 0.0).element == "Cl"
    assert sketch.add_atom("  na ", 0.0, 0.0).element == "Na"


def test_add_atom_rejects_a_blank_element():
    sketch = MoleculeSketch()

    with pytest.raises(ValueError, match="element"):
        sketch.add_atom("   ", 0.0, 0.0)


def test_atom_ids_are_unique_and_never_reused_after_deletion():
    sketch = MoleculeSketch()
    first = sketch.add_atom("C", 0.0, 0.0)
    second = sketch.add_atom("N", 10.0, 0.0)

    sketch.remove_atom(first.atom_id)
    third = sketch.add_atom("O", 20.0, 0.0)

    assert len({first.atom_id, second.atom_id, third.atom_id}) == 3
    assert third.atom_id != first.atom_id


def test_removing_an_atom_removes_every_bond_that_touched_it():
    sketch = _water()
    oxygen_id = sketch.atoms[0].atom_id

    removed = sketch.remove_atom(oxygen_id)

    assert removed is True
    assert sketch.atom_count == 2
    assert sketch.bond_count == 0


def test_removing_an_atom_leaves_the_other_atom_ids_unchanged():
    """Ids, not list indices: surviving bonds must still resolve."""
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 40.0, 0.0)
    c = sketch.add_atom("C", 80.0, 0.0)
    sketch.add_bond(b.atom_id, c.atom_id)

    sketch.remove_atom(a.atom_id)

    assert [atom.atom_id for atom in sketch.atoms] == [b.atom_id, c.atom_id]
    assert sketch.find_bond(b.atom_id, c.atom_id) is not None


def test_removing_an_unknown_atom_reports_false():
    sketch = MoleculeSketch()

    assert sketch.remove_atom(999) is False


# ------------------------------------------------------------------
# Bonds
# ------------------------------------------------------------------
def test_add_bond_creates_a_single_bond():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 40.0, 0.0)

    bond = sketch.add_bond(a.atom_id, b.atom_id)

    assert isinstance(bond, SketchBond)
    assert bond.order == 1
    assert sketch.bond_count == 1


def test_add_bond_returns_none_when_the_pair_is_already_bonded():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 40.0, 0.0)
    sketch.add_bond(a.atom_id, b.atom_id)

    assert sketch.add_bond(a.atom_id, b.atom_id) is None
    assert sketch.add_bond(b.atom_id, a.atom_id) is None
    assert sketch.bond_count == 1


def test_add_bond_rejects_an_atom_bonded_to_itself():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)

    with pytest.raises(ValueError, match="itself"):
        sketch.add_bond(a.atom_id, a.atom_id)


def test_add_bond_rejects_an_unknown_atom_id():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)

    with pytest.raises(ValueError, match="Unknown atom"):
        sketch.add_bond(a.atom_id, 4242)


def test_find_bond_is_order_independent():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("O", 40.0, 0.0)
    sketch.add_bond(a.atom_id, b.atom_id)

    assert sketch.find_bond(a.atom_id, b.atom_id) is sketch.find_bond(b.atom_id, a.atom_id)


def test_remove_bond_reports_whether_anything_was_removed():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 40.0, 0.0)
    sketch.add_bond(a.atom_id, b.atom_id)

    assert sketch.remove_bond(a.atom_id, b.atom_id) is True
    assert sketch.remove_bond(a.atom_id, b.atom_id) is False


# ------------------------------------------------------------------
# Bond-order cycling
# ------------------------------------------------------------------
def test_cycling_a_carbon_carbon_bond_passes_through_resonance_before_single():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 40.0, 0.0)
    bond = sketch.add_bond(a.atom_id, b.atom_id)

    assert sketch.cycle_bond_order(bond) == 2
    assert sketch.cycle_bond_order(bond) == 3
    assert sketch.cycle_bond_order(bond) == RESONANCE_ORDER
    assert sketch.cycle_bond_order(bond) == 1
    assert bond.order == 1


def test_an_oxygen_hydrogen_bond_never_leaves_single():
    sketch = MoleculeSketch()
    o = sketch.add_atom("O", 0.0, 0.0)
    h = sketch.add_atom("H", 40.0, 0.0)
    bond = sketch.add_bond(o.atom_id, h.atom_id)

    assert sketch.cycle_bond_order(bond) == 1


def test_a_carbon_oxygen_bond_skips_triple_but_still_reaches_resonance():
    sketch = MoleculeSketch()
    c = sketch.add_atom("C", 0.0, 0.0)
    o = sketch.add_atom("O", 40.0, 0.0)
    bond = sketch.add_bond(c.atom_id, o.atom_id)

    assert sketch.cycle_bond_order(bond) == 2
    assert sketch.cycle_bond_order(bond) == RESONANCE_ORDER
    assert sketch.cycle_bond_order(bond) == 1


def test_max_bond_order_is_the_lower_of_the_two_elements():
    sketch = MoleculeSketch()

    assert sketch.max_bond_order_for("C", "C") == 3
    assert sketch.max_bond_order_for("C", "O") == 2
    assert sketch.max_bond_order_for("C", "H") == 1
    assert sketch.max_bond_order_for("N", "N") == 3


def test_an_unlisted_element_falls_back_to_triple():
    sketch = MoleculeSketch()

    assert sketch.max_bond_order_for("C", "Fm") == 3


# ------------------------------------------------------------------
# Hit-testing
# ------------------------------------------------------------------
def test_atom_at_finds_an_atom_under_the_cursor():
    sketch = _water()

    hit = sketch.atom_at(102.0, 98.0)

    assert hit is not None
    assert hit.element == "O"


def test_atom_at_returns_none_on_empty_space():
    sketch = _water()

    assert sketch.atom_at(400.0, 400.0) is None


def test_atom_at_returns_the_nearest_of_two_overlapping_atoms():
    sketch = MoleculeSketch()
    near = sketch.add_atom("C", 100.0, 100.0)
    sketch.add_atom("N", 108.0, 100.0)

    assert sketch.atom_at(101.0, 100.0).atom_id == near.atom_id


def test_bond_at_finds_a_bond_by_its_midpoint():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 100.0, 0.0)
    bond = sketch.add_bond(a.atom_id, b.atom_id)

    assert sketch.bond_at(50.0, 2.0) is bond


def test_bond_at_ignores_a_point_beyond_the_ends_of_the_bond():
    """Point-segment distance, not point-line: past the end is a miss."""
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 100.0, 0.0)
    sketch.add_bond(a.atom_id, b.atom_id)

    assert sketch.bond_at(180.0, 0.0) is None


def test_bond_at_returns_none_far_from_the_bond_axis():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 100.0, 0.0)
    sketch.add_bond(a.atom_id, b.atom_id)

    assert sketch.bond_at(50.0, 60.0) is None


# ------------------------------------------------------------------
# Hit-testing through a projected (rotated) view
# ------------------------------------------------------------------
def test_atom_at_uses_projected_screen_positions_when_given():
    sketch = MoleculeSketch()
    atom = sketch.add_atom("C", 0.0, 0.0)
    screen = {atom.atom_id: (300.0, 200.0, 0.0)}

    assert sketch.atom_at(300.0, 200.0, screen).atom_id == atom.atom_id
    assert sketch.atom_at(0.0, 0.0, screen) is None


def test_atom_at_prefers_the_atom_in_front_when_two_overlap_on_screen():
    """After a rotation two atoms can sit on top of each other; the one drawn
    in front is the one the user sees, so it is the one a click selects."""
    sketch = MoleculeSketch()
    back = sketch.add_atom("C", 0.0, 0.0, 50.0)
    front = sketch.add_atom("O", 0.0, 0.0, -50.0)
    screen = {
        back.atom_id: (100.0, 100.0, 50.0),
        front.atom_id: (103.0, 100.0, -50.0),
    }

    assert sketch.atom_at(101.0, 100.0, screen).atom_id == front.atom_id


def test_atom_at_still_prefers_a_disc_hit_over_a_nearer_depth_near_miss():
    """An atom the click is inside beats one merely within the tolerance,
    whatever their depths."""
    sketch = MoleculeSketch()
    clicked = sketch.add_atom("C", 0.0, 0.0, 50.0)
    nearby = sketch.add_atom("O", 0.0, 0.0, -50.0)
    screen = {
        clicked.atom_id: (100.0, 100.0, 50.0),
        nearby.atom_id: (113.0, 100.0, -50.0),
    }

    assert sketch.atom_at(100.0, 100.0, screen).atom_id == clicked.atom_id


def test_bond_at_uses_projected_screen_positions_when_given():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 100.0, 0.0)
    bond = sketch.add_bond(a.atom_id, b.atom_id)
    screen = {a.atom_id: (0.0, 300.0, 0.0), b.atom_id: (0.0, 400.0, 0.0)}

    assert sketch.bond_at(1.0, 350.0, screen) is bond
    assert sketch.bond_at(50.0, 0.0, screen) is None


# ------------------------------------------------------------------
# 3D positions
# ------------------------------------------------------------------
def test_atoms_default_to_depth_zero():
    assert MoleculeSketch().add_atom("C", 1.0, 2.0).z == 0.0


def test_add_atom_accepts_a_depth():
    assert MoleculeSketch().add_atom("C", 1.0, 2.0, -7.5).z == -7.5


def test_positions_are_returned_in_atom_order():
    sketch = MoleculeSketch()
    sketch.add_atom("C", 1.0, 2.0, 3.0)
    sketch.add_atom("O", 4.0, 5.0, 6.0)

    assert sketch.positions().tolist() == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]


def test_positions_of_an_empty_sketch_have_three_columns():
    assert MoleculeSketch().positions().shape == (0, 3)


def test_update_positions_writes_back_in_atom_order_and_keeps_ids():
    sketch = _water()
    ids_before = sketch.ordered_atom_ids()
    new = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]

    sketch.update_positions(new)

    assert sketch.positions().tolist() == new
    assert sketch.ordered_atom_ids() == ids_before
    assert sketch.bond_count == 2


def test_update_positions_survives_a_deleted_atom():
    """Ids are not positions: writing back after a deletion must follow the
    current atom order, not the ids."""
    sketch = _water()
    sketch.remove_atom(sketch.ordered_atom_ids()[1])

    sketch.update_positions([[0.0, 0.0, 0.0], [9.0, 9.0, 9.0]])

    assert sketch.atoms[1].element == "H"
    assert (sketch.atoms[1].x, sketch.atoms[1].y, sketch.atoms[1].z) == (9.0, 9.0, 9.0)


def test_update_positions_rejects_the_wrong_atom_count():
    with pytest.raises(ValueError, match="3 atoms"):
        _water().update_positions([[0.0, 0.0, 0.0]])


def test_update_positions_rejects_non_finite_values():
    with pytest.raises(ValueError, match="finite"):
        _water().update_positions([[0.0, 0.0, 0.0], [1.0, float("nan"), 0.0], [2.0, 0.0, 0.0]])


def test_a_copy_is_independent_of_the_original():
    sketch = _water()
    duplicate = sketch.copy()

    duplicate.update_positions([[0.0, 0.0, 0.0]] * 3)
    duplicate.add_atom("C", 0.0, 0.0)

    assert sketch.atom_count == 3
    assert sketch.positions()[0].tolist() == [100.0, 100.0, 0.0]
    assert duplicate.ordered_atom_ids()[:3] == sketch.ordered_atom_ids()


def test_undo_snapshots_keep_the_depth():
    sketch = MoleculeSketch()
    sketch.add_atom("C", 1.0, 2.0, 3.0)
    state = sketch.snapshot()
    sketch.update_positions([[0.0, 0.0, 0.0]])

    sketch.restore(state)

    assert sketch.atoms[0].z == 3.0


# ------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------
def test_a_well_formed_molecule_has_no_errors_and_no_warnings():
    errors, warnings = _water().validate()

    assert errors == []
    assert warnings == []


def test_an_empty_sketch_is_an_error():
    errors, _ = MoleculeSketch().validate()

    assert any("no atoms" in message.lower() for message in errors)


def test_an_unbonded_atom_beside_a_molecule_is_an_error():
    sketch = _water()
    sketch.add_atom("C", 300.0, 300.0)

    errors, _ = sketch.validate()

    assert any("not bonded" in message.lower() for message in errors)


def test_a_single_lone_atom_is_a_valid_structure():
    sketch = MoleculeSketch()
    sketch.add_atom("Ne", 10.0, 10.0)

    errors, _ = sketch.validate()

    assert errors == []


def test_an_element_without_a_covalent_radius_is_an_error():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("Xx", 40.0, 0.0)
    sketch.add_bond(a.atom_id, b.atom_id)

    errors, _ = sketch.validate()

    assert any("Xx" in message for message in errors)


def test_an_over_valent_carbon_is_a_warning_not_an_error():
    sketch = MoleculeSketch()
    carbon = sketch.add_atom("C", 100.0, 100.0)
    for index in range(5):
        hydrogen = sketch.add_atom("H", 100.0 + 40.0 * index, 160.0)
        sketch.add_bond(carbon.atom_id, hydrogen.atom_id)

    errors, warnings = sketch.validate()

    assert errors == []
    assert any("valence" in message.lower() for message in warnings)


def test_bond_orders_count_towards_the_valence_warning():
    """Three double bonds on one carbon is 6 > 4, even though it is 3 bonds."""
    sketch = MoleculeSketch()
    carbon = sketch.add_atom("C", 100.0, 100.0)
    for index in range(3):
        oxygen = sketch.add_atom("O", 100.0 + 40.0 * index, 160.0)
        bond = sketch.add_bond(carbon.atom_id, oxygen.atom_id)
        sketch.cycle_bond_order(bond)

    _, warnings = sketch.validate()

    assert any("valence" in message.lower() for message in warnings)


def test_two_disconnected_molecules_are_a_warning():
    sketch = _water()
    a = sketch.add_atom("C", 400.0, 400.0)
    b = sketch.add_atom("O", 440.0, 400.0)
    sketch.add_bond(a.atom_id, b.atom_id)

    errors, warnings = sketch.validate()

    assert errors == []
    assert any("fragment" in message.lower() for message in warnings)


def test_validation_names_the_offending_atom_by_its_export_number():
    sketch = MoleculeSketch()
    first = sketch.add_atom("C", 0.0, 0.0)
    carbon = sketch.add_atom("C", 40.0, 0.0)
    sketch.add_bond(first.atom_id, carbon.atom_id)
    for index in range(4):
        hydrogen = sketch.add_atom("H", 40.0 + 30.0 * index, 60.0)
        sketch.add_bond(carbon.atom_id, hydrogen.atom_id)

    _, warnings = sketch.validate()

    assert any("atom 2 (C)" in message for message in warnings)


# ------------------------------------------------------------------
# Undo support and export ordering
# ------------------------------------------------------------------
def test_restoring_a_snapshot_undoes_later_edits():
    sketch = _water()
    saved = sketch.snapshot()
    sketch.add_atom("C", 500.0, 500.0)

    sketch.restore(saved)

    assert sketch.atom_count == 3
    assert sketch.bond_count == 2


def test_a_snapshot_is_not_affected_by_later_edits():
    sketch = _water()
    saved = sketch.snapshot()

    sketch.atoms[0].x = 999.0
    sketch.restore(saved)

    assert sketch.atoms[0].x == 100.0


def test_ordered_atom_ids_follow_insertion_order():
    sketch = _water()

    assert sketch.ordered_atom_ids() == [atom.atom_id for atom in sketch.atoms]


def test_clear_empties_the_sketch():
    sketch = _water()

    sketch.clear()

    assert sketch.atom_count == 0
    assert sketch.bond_count == 0


def test_the_model_never_imports_toga():
    """The sketch must stay GUI-free so it can be tested headless."""
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-c", "import molecularSketch, sys; print('toga' in sys.modules)"],
        capture_output=True,
        text=True,
    )

    assert result.stdout.strip() == "False", result.stderr


def test_bond_index_pairs_are_zero_based_positions_not_ids():
    """The viewer indexes atoms by position, so ids must be mapped out."""
    sketch = MoleculeSketch()
    first = sketch.add_atom("C", 0.0, 0.0)
    second = sketch.add_atom("C", 40.0, 0.0)
    third = sketch.add_atom("O", 80.0, 0.0)
    sketch.add_bond(second.atom_id, third.atom_id)
    sketch.remove_atom(first.atom_id)

    assert sketch.bond_index_pairs() == [(0, 1)]


def test_bond_index_pairs_are_ordered_low_to_high():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 40.0, 0.0)
    sketch.add_bond(b.atom_id, a.atom_id)

    assert sketch.bond_index_pairs() == [(0, 1)]


# ------------------------------------------------------------------
# Resonance bonds
# ------------------------------------------------------------------
def _benzene() -> MoleculeSketch:
    """Six carbons in a ring, every bond a resonance bond."""
    sketch = MoleculeSketch()
    carbons = [sketch.add_atom("C", 100.0 + 30.0 * i, 100.0 + 20.0 * (i % 2)) for i in range(6)]
    for index in range(6):
        bond = sketch.add_bond(carbons[index].atom_id, carbons[(index + 1) % 6].atom_id)
        bond.order = RESONANCE_ORDER
    return sketch


def test_a_resonance_bond_is_order_one_and_a_half():
    assert RESONANCE_ORDER == 1.5


def test_an_oxygen_hydrogen_bond_never_reaches_resonance():
    sketch = MoleculeSketch()
    o = sketch.add_atom("O", 0.0, 0.0)
    h = sketch.add_atom("H", 40.0, 0.0)
    bond = sketch.add_bond(o.atom_id, h.atom_id)

    orders = {sketch.cycle_bond_order(bond) for _ in range(6)}

    assert orders == {1}


def test_resonance_is_offered_wherever_a_double_bond_is():
    sketch = MoleculeSketch()

    assert RESONANCE_ORDER in sketch.bond_order_cycle(3)
    assert RESONANCE_ORDER in sketch.bond_order_cycle(2)
    assert RESONANCE_ORDER not in sketch.bond_order_cycle(1)


def test_resonance_sits_last_in_the_cycle():
    sketch = MoleculeSketch()

    assert sketch.bond_order_cycle(3) == (1.0, 2.0, 3.0, RESONANCE_ORDER)
    assert sketch.bond_order_cycle(2) == (1.0, 2.0, RESONANCE_ORDER)


def test_an_unexpected_order_cycles_back_to_single():
    """Defensive: an order outside the cycle must not strand the bond."""
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 40.0, 0.0)
    bond = sketch.add_bond(a.atom_id, b.atom_id)
    bond.order = 2.25

    assert sketch.cycle_bond_order(bond) == 1


def test_resonance_counts_as_one_and_a_half_towards_the_valence():
    sketch = _benzene()
    carbon_id = sketch.atoms[0].atom_id

    assert sketch.bond_order_sum(carbon_id) == 3.0


def test_an_aromatic_carbon_with_a_substituent_is_not_over_valent():
    """1.5 + 1.5 + 1 is exactly 4: benzene must not warn about its own carbons."""
    sketch = _benzene()
    hydrogen = sketch.add_atom("H", 300.0, 300.0)
    sketch.add_bond(sketch.atoms[0].atom_id, hydrogen.atom_id)

    _, warnings = sketch.validate()

    assert not any("valence" in message.lower() for message in warnings)


def test_bond_type_names_cover_every_order():
    assert bond_type_name(1) == "single"
    assert bond_type_name(2) == "double"
    assert bond_type_name(3) == "triple"
    assert bond_type_name(RESONANCE_ORDER) == "resonance"


def test_format_order_drops_a_trailing_zero_but_keeps_a_half():
    assert format_order(4.0) == "4"
    assert format_order(3.5) == "3.5"


def test_a_valence_warning_reads_as_a_whole_number():
    sketch = MoleculeSketch()
    carbon = sketch.add_atom("C", 100.0, 100.0)
    for index in range(5):
        hydrogen = sketch.add_atom("H", 100.0 + 30.0 * index, 160.0)
        sketch.add_bond(carbon.atom_id, hydrogen.atom_id)

    _, warnings = sketch.validate()

    assert any("5 bonds where 4 is typical" in message for message in warnings)


# ------------------------------------------------------------------
# Resonance validation
# ------------------------------------------------------------------
def test_an_isolated_resonance_bond_is_a_warning():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 40.0, 0.0)
    bond = sketch.add_bond(a.atom_id, b.atom_id)
    bond.order = RESONANCE_ORDER

    errors, warnings = sketch.validate()

    assert errors == []
    assert any("resonance" in message.lower() for message in warnings)


def test_a_resonance_ring_raises_no_resonance_warning():
    _, warnings = _benzene().validate()

    assert not any("resonance" in message.lower() for message in warnings)


def test_a_resonance_chain_raises_no_warning():
    """A carboxylate: O~C~O is delocalized without being a ring."""
    sketch = MoleculeSketch()
    carbon = sketch.add_atom("C", 100.0, 100.0)
    for x in (60.0, 140.0):
        oxygen = sketch.add_atom("O", x, 160.0)
        bond = sketch.add_bond(carbon.atom_id, oxygen.atom_id)
        bond.order = RESONANCE_ORDER

    _, warnings = sketch.validate()

    assert not any("resonance" in message.lower() for message in warnings)


def test_the_resonance_warning_names_both_atoms():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("N", 40.0, 0.0)
    bond = sketch.add_bond(a.atom_id, b.atom_id)
    bond.order = RESONANCE_ORDER

    _, warnings = sketch.validate()

    message = next(m for m in warnings if "resonance" in m.lower())
    assert "atom 1 (C)" in message
    assert "atom 2 (N)" in message


def test_a_resonance_bond_with_a_neighbour_at_one_end_is_enough():
    sketch = MoleculeSketch()
    atoms = [sketch.add_atom("C", 40.0 * i, 0.0) for i in range(3)]
    for index in range(2):
        bond = sketch.add_bond(atoms[index].atom_id, atoms[index + 1].atom_id)
        bond.order = RESONANCE_ORDER

    _, warnings = sketch.validate()

    assert not any("resonance" in message.lower() for message in warnings)


# ------------------------------------------------------------------
# Replacing the drawing with an imported structure
# ------------------------------------------------------------------
def test_replace_contents_loads_atoms_and_bonds_in_order():
    sketch = _water()
    sketch.replace_contents(
        ["c", "O"], [[0.0, 0.0, 0.0], [48.0, 0.0, 4.0]], [(0, 1, 2.0)]
    )

    assert [atom.element for atom in sketch.atoms] == ["C", "O"]
    assert sketch.positions()[1].tolist() == [48.0, 0.0, 4.0]
    assert sketch.bond_index_pairs() == [(0, 1)]
    assert sketch.bonds[0].order == 2.0


def test_replace_contents_never_reuses_atom_ids():
    sketch = _water()
    old_ids = set(sketch.ordered_atom_ids())

    sketch.replace_contents(["H", "H"], [[0.0, 0.0, 0.0], [30.0, 0.0, 0.0]])

    assert not old_ids & set(sketch.ordered_atom_ids())


@pytest.mark.parametrize("elements, positions, bonds, message", [
    (["C", "C"], [[0.0, 0.0, 0.0]], [], "Expected positions for 2 atoms"),
    (["C"], [[0.0, float("nan"), 0.0]], [], "finite"),
    (["C", "C"], [[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]], [(0, 2, 1.0)], "refers to atom 3"),
    (["C", "C"], [[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]], [(1, 1, 1.0)], "itself"),
    (["C", "C"], [[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]], [(0, 1, 4.0)], "order 4"),
    (["C", "C"], [[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]], [(0, 1, 1.0), (1, 0, 2.0)], "bonded twice"),
    ([""], [[0.0, 0.0, 0.0]], [], "element symbol is required"),
])
def test_a_refused_structure_leaves_the_drawing_untouched(elements, positions, bonds, message):
    sketch = _water()
    before = sketch.snapshot()

    with pytest.raises(ValueError, match=message):
        sketch.replace_contents(elements, positions, bonds)

    assert sketch.snapshot() == before


def test_add_structure_keeps_the_existing_drawing():
    sketch = _water()
    new_ids = sketch.add_structure(["C", "O"], [[300.0, 0.0, 0.0], [348.0, 0.0, 0.0]], [(0, 1, 2.0)])

    assert sketch.atom_count == 5 and sketch.bond_count == 3
    assert sketch.find_bond(*new_ids).order == 2.0


# ------------------------------------------------------------------
# Frozen parameters
# ------------------------------------------------------------------
def test_a_parameter_is_the_same_read_backwards():
    sketch = _water()
    o, h1, h2 = sketch.ordered_atom_ids()
    sketch.set_constraint([h1, o, h2], 104.5)

    assert sketch.find_constraint([h2, o, h1]).value == 104.5
    sketch.set_constraint([h2, o, h1], 100.0)  # replaces, never duplicates
    assert len(sketch.constraints) == 1 and sketch.constraints[0].value == 100.0
    assert sketch.constraints[0].kind == "angle"


@pytest.mark.parametrize("atoms, message", [
    ([1], "2, 3 or 4 atoms"), ([1, 1], "same atom twice"), ([1, 99], "Unknown atom id"),
])
def test_bad_constraints_are_refused(atoms, message):
    with pytest.raises(ValueError, match=message):
        _water().set_constraint(atoms, 1.0)


def test_deleting_an_atom_removes_its_constraints():
    sketch = _water()
    o, h1, h2 = sketch.ordered_atom_ids()
    sketch.set_constraint([o, h1], 0.96)
    sketch.set_constraint([o, h2], 0.96)

    sketch.remove_atom(h1)

    assert [c.atom_ids for c in sketch.constraints] == [(o, h2)]


def test_undo_snapshots_and_copies_carry_constraints():
    sketch = _water()
    o, h1, _ = sketch.ordered_atom_ids()
    state = sketch.snapshot()
    sketch.set_constraint([o, h1], 0.96)

    assert sketch.copy().constraints == sketch.constraints
    sketch.restore(state)
    assert sketch.constraints == []


def test_replacing_an_element_keeps_id_position_bonds_and_constraints():
    sketch = _water()
    o, h1, h2 = sketch.ordered_atom_ids()
    sketch.find_bond(o, h1).order = 2.0
    sketch.set_constraint([h1, o, h2], 104.5)
    before_positions = sketch.positions()

    previous = sketch.replace_element(o, "s")

    assert previous == "O"
    assert sketch.get_atom(o).element == "S"
    assert sketch.ordered_atom_ids() == [o, h1, h2]
    assert (sketch.positions() == before_positions).all()
    assert sketch.find_bond(o, h1).order == 2.0 and sketch.find_bond(o, h2) is not None
    assert sketch.constraints[0].atom_ids == (h1, o, h2)


def test_replacing_an_unknown_atom_is_refused():
    with pytest.raises(ValueError, match="Unknown atom id"):
        _water().replace_element(99, "C")


def test_the_selected_molecule_is_its_connected_fragment():
    sketch = _water()
    sketch.add_structure(["C"], [[400.0, 0.0, 0.0]])
    o = sketch.ordered_atom_ids()[0]

    assert sorted(sketch.fragment_containing(o)) == sorted(sketch.ordered_atom_ids()[:3])
    assert sketch.describe_parameter(sketch.ordered_atom_ids()[:2]) == "1-2 (O-H)"
