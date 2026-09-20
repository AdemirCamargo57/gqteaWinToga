"""Tests for the molecularSketch data model.

Pure model tests: no Toga, no OpenGL, no window. The sketch is the single
source of truth for what the user has drawn, so every rule the Design tab
relies on (stable ids, bond-order cycling, hit-testing, validation) is
pinned here.
"""
import pytest

from molecularSketch import MoleculeSketch, SketchAtom, SketchBond


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
def test_cycling_a_carbon_carbon_bond_goes_single_double_triple_then_back():
    sketch = MoleculeSketch()
    a = sketch.add_atom("C", 0.0, 0.0)
    b = sketch.add_atom("C", 40.0, 0.0)
    bond = sketch.add_bond(a.atom_id, b.atom_id)

    assert sketch.cycle_bond_order(bond) == 2
    assert sketch.cycle_bond_order(bond) == 3
    assert sketch.cycle_bond_order(bond) == 1
    assert bond.order == 1


def test_an_oxygen_hydrogen_bond_never_leaves_single():
    sketch = MoleculeSketch()
    o = sketch.add_atom("O", 0.0, 0.0)
    h = sketch.add_atom("H", 40.0, 0.0)
    bond = sketch.add_bond(o.atom_id, h.atom_id)

    assert sketch.cycle_bond_order(bond) == 1


def test_a_carbon_oxygen_bond_caps_at_double():
    sketch = MoleculeSketch()
    c = sketch.add_atom("C", 0.0, 0.0)
    o = sketch.add_atom("O", 40.0, 0.0)
    bond = sketch.add_bond(c.atom_id, o.atom_id)

    assert sketch.cycle_bond_order(bond) == 2
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
