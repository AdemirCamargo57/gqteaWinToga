"""Tests for molecularGeometry: measuring, editing and constraining parameters.

The analytic gradients are what hold frozen parameters exactly, so they are
pinned against finite differences; the edits are checked to land on the
requested value while leaving every other bond length unchanged.
"""
import subprocess
import sys

import numpy as np
import pytest

from molecularGeometry import (
    ANGLE,
    DIHEDRAL,
    DISTANCE,
    Constraint,
    GeometryEditError,
    canonical_atoms,
    constraint_residuals,
    format_value,
    is_degenerate,
    measure,
    moving_atoms,
    parse_parameter_value,
    project_onto_constraints,
    set_parameter,
    tangent_component,
    value_and_gradient,
    wrap_degrees,
)

# Butane-like zigzag chain C1-C2-C3-C4 with an H on C1 and one on C4.
CHAIN = np.array([
    [0.0, 0.0, 0.0], [1.52, 0.0, 0.0], [2.03, 1.43, 0.0], [3.55, 1.43, 0.0],
    [-0.36, -1.03, 0.0], [3.91, 2.46, 0.0],
])
CHAIN_BONDS = [(0, 1), (1, 2), (2, 3), (0, 4), (3, 5)]


def _ring(size=6, radius=1.4):
    angles = np.radians(np.arange(size) * 360.0 / size)
    return np.c_[radius * np.cos(angles), radius * np.sin(angles), np.zeros(size)], \
        [(k, (k + 1) % size) for k in range(size)]


def test_the_module_does_not_import_toga():
    result = subprocess.run(
        [sys.executable, "-c", "import molecularGeometry, sys; print('toga' in sys.modules)"],
        capture_output=True, text=True,
    )
    assert result.stdout.strip() == "False", result.stderr


# ------------------------------------------------------------------
# Values
# ------------------------------------------------------------------
def test_simple_values():
    square = np.array([[0.0, 0, 0], [1.0, 0, 0], [1.0, 1.0, 0], [2.0, 1.0, 0]])
    assert measure(DISTANCE, square, (0, 1)) == pytest.approx(1.0)
    assert measure(ANGLE, square, (0, 1, 2)) == pytest.approx(90.0)
    assert measure(DIHEDRAL, square, (0, 1, 2, 3)) == pytest.approx(180.0)


def test_dihedral_sign_follows_iupac():
    """Looking along 2->3, atom 4 clockwise from atom 1 is positive."""
    points = np.array([[1.0, 0, 0], [0.0, 0, 0], [0.0, 0, 1.0], [0.0, 1.0, 1.0]])
    assert measure(DIHEDRAL, points, (0, 1, 2, 3)) == pytest.approx(90.0)
    points[3] = [0.0, -1.0, 1.0]
    assert measure(DIHEDRAL, points, (0, 1, 2, 3)) == pytest.approx(-90.0)


@pytest.mark.parametrize("kind, atoms", [(DISTANCE, (0, 2)), (ANGLE, (0, 1, 3)), (DIHEDRAL, (0, 1, 2, 3))])
def test_gradients_match_finite_differences(kind, atoms):
    coordinates = np.random.default_rng(7).normal(size=(5, 3))
    _, gradient = value_and_gradient(kind, coordinates, atoms)
    numeric = np.zeros_like(coordinates)
    step = 1e-6
    scale = 1.0 if kind == DISTANCE else np.pi / 180.0
    for atom in range(5):
        for axis in range(3):
            plus, minus = coordinates.copy(), coordinates.copy()
            plus[atom, axis] += step
            minus[atom, axis] -= step
            numeric[atom, axis] = (measure(kind, plus, atoms) - measure(kind, minus, atoms)) / (2 * step) * scale
    assert np.allclose(gradient, numeric, atol=1e-7)


def test_formatting_and_parsing():
    assert format_value(DISTANCE, 1.52) == "1.5200 Å"
    assert format_value(ANGLE, 109.4712) == "109.47°"
    assert parse_parameter_value(DISTANCE, " 1,45 Å") == pytest.approx(1.45)
    assert parse_parameter_value(DIHEDRAL, "270") == pytest.approx(-90.0)
    assert wrap_degrees(-180.0) == 180.0
    assert canonical_atoms(DIHEDRAL, (4, 3, 2, 1)) == (1, 2, 3, 4)


@pytest.mark.parametrize("kind, text, message", [
    (DISTANCE, "", "Enter a value"),
    (DISTANCE, "abc", "not a number"),
    (DISTANCE, "0.1", "between 0.3 and 20"),
    (ANGLE, "0", "greater than 0"),
    (ANGLE, "181", "at most 180"),
    (DIHEDRAL, "inf", "finite"),
])
def test_bad_values_are_refused(kind, text, message):
    with pytest.raises(ValueError, match=message):
        parse_parameter_value(kind, text)


def test_degenerate_parameters_are_named():
    line = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0], [2.0, 1.0, 0]])
    assert "in a line" in is_degenerate(ANGLE, line, (0, 1, 2))
    assert "atoms 1-2-3" in is_degenerate(DIHEDRAL, line, (0, 1, 2, 3))
    assert is_degenerate(ANGLE, CHAIN, (0, 1, 2)) is None


# ------------------------------------------------------------------
# Editing
# ------------------------------------------------------------------
def _bond_lengths(coordinates, bonds):
    return [measure(DISTANCE, coordinates, bond) for bond in bonds]


@pytest.mark.parametrize("kind, atoms, value", [
    (DISTANCE, (1, 2), 1.80), (ANGLE, (0, 1, 2), 95.0), (DIHEDRAL, (0, 1, 2, 3), 65.0),
    (DIHEDRAL, (0, 1, 2, 3), -120.0), (ANGLE, (1, 2, 3), 180.0),
])
def test_an_edit_lands_on_the_value_and_keeps_other_bonds(kind, atoms, value):
    edited = set_parameter(kind, CHAIN, CHAIN_BONDS, atoms, value)
    assert measure(kind, edited, atoms) == pytest.approx(value, abs=1e-9)
    others = [bond for bond in CHAIN_BONDS if set(bond) != set(atoms)]
    assert np.allclose(_bond_lengths(edited, others), _bond_lengths(CHAIN, others))


def test_only_the_last_atoms_side_moves():
    edited = set_parameter(DIHEDRAL, CHAIN, CHAIN_BONDS, (0, 1, 2, 3), 60.0)
    assert np.array_equal(edited[[0, 1, 2, 4]], CHAIN[[0, 1, 2, 4]])
    assert moving_atoms(DIHEDRAL, 6, CHAIN_BONDS, (0, 1, 2, 3)) == {3, 5}


def test_a_linear_angle_can_be_bent():
    line = np.array([[0.0, 0, 0], [1.2, 0, 0], [2.4, 0, 0]])
    edited = set_parameter(ANGLE, line, [(0, 1), (1, 2)], (0, 1, 2), 120.0)
    assert measure(ANGLE, edited, (0, 1, 2)) == pytest.approx(120.0)


def test_a_ring_bond_is_refused():
    ring, bonds = _ring()
    with pytest.raises(GeometryEditError, match="part of a ring"):
        set_parameter(DISTANCE, ring, bonds, (0, 1), 1.6)
    with pytest.raises(GeometryEditError, match="part of a ring"):
        set_parameter(DIHEDRAL, ring, bonds, (5, 0, 1, 2), 30.0)


def test_unbonded_atoms_in_one_molecule_are_refused():
    with pytest.raises(GeometryEditError, match="not bonded"):
        set_parameter(DISTANCE, CHAIN, CHAIN_BONDS, (0, 2), 3.0)
    with pytest.raises(GeometryEditError, match="not bonded"):
        set_parameter(ANGLE, CHAIN, CHAIN_BONDS, (0, 1, 3), 90.0)


def test_a_distance_between_two_molecules_moves_the_second():
    pair = np.vstack([CHAIN[:2], CHAIN[:2] + [0.0, 5.0, 0.0]])
    bonds = [(0, 1), (2, 3)]
    edited = set_parameter(DISTANCE, pair, bonds, (1, 2), 4.0)
    assert measure(DISTANCE, edited, (1, 2)) == pytest.approx(4.0)
    assert np.array_equal(edited[:2], pair[:2])
    assert measure(DISTANCE, edited, (2, 3)) == pytest.approx(1.52)


def test_an_undefined_dihedral_cannot_be_set():
    line = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0], [2.0, 1.0, 0]])
    with pytest.raises(GeometryEditError, match="undefined"):
        set_parameter(DIHEDRAL, line, [(0, 1), (1, 2), (2, 3)], (0, 1, 2, 3), 60.0)


# ------------------------------------------------------------------
# Holding constraints
# ------------------------------------------------------------------
def test_projection_satisfies_every_constraint():
    constraints = [Constraint(DISTANCE, (0, 1), 1.7), Constraint(ANGLE, (0, 1, 2), 100.0),
                   Constraint(DIHEDRAL, (0, 1, 2, 3), 75.0)]
    projected = project_onto_constraints(CHAIN, constraints)
    residuals, _ = constraint_residuals(projected, constraints)
    assert np.max(np.abs(residuals)) < 1e-10


def test_contradictory_constraints_are_reported():
    constraints = [Constraint(DISTANCE, (0, 1), 1.0), Constraint(DISTANCE, (1, 2), 1.0),
                   Constraint(DISTANCE, (0, 2), 3.0)]  # no triangle has these sides
    with pytest.raises(ValueError, match="cannot all hold"):
        project_onto_constraints(CHAIN, constraints)


def test_the_tangent_component_does_not_change_constraints():
    constraints = [Constraint(DISTANCE, (0, 1), 1.52), Constraint(DIHEDRAL, (0, 1, 2, 3), 180.0)]
    vector = np.random.default_rng(3).normal(size=CHAIN.size)
    _, jacobian = constraint_residuals(CHAIN, constraints)
    assert np.allclose(jacobian @ tangent_component(vector, CHAIN, constraints), 0.0, atol=1e-10)
