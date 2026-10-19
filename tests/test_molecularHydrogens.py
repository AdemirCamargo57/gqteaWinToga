"""Tests for molecularHydrogens: how many hydrogens, where, and when to refuse."""
import subprocess
import sys

import numpy as np
import pytest

from molecularGeometry import ANGLE, DIHEDRAL, DISTANCE, measure
from molecularHydrogens import hydrogens_needed, plan_hydrogens
from molecularPreOptimizer import reference_bond_length


def _apply(elements, coordinates, bonds, scope=None):
    """The structure after the plan, as (elements, coordinates, bonds)."""
    plan = plan_hydrogens(elements, coordinates, bonds, scope)
    elements = list(elements) + ["H"] * plan.count
    coordinates = np.vstack([np.asarray(coordinates, dtype=float).reshape(-1, 3)]
                            + [position for _, position in plan.additions])
    start = len(elements) - plan.count
    bonds = list(bonds) + [(parent, start + k, 1.0) for k, (parent, _) in enumerate(plan.additions)]
    return plan, elements, coordinates, bonds


def _ring(size=6, radius=1.39):
    angles = np.radians(np.arange(size) * 360.0 / size)
    return np.c_[radius * np.cos(angles), radius * np.sin(angles), np.zeros(size)]


def test_the_module_does_not_import_toga():
    result = subprocess.run(
        [sys.executable, "-c", "import molecularHydrogens, sys; print('toga' in sys.modules)"],
        capture_output=True, text=True,
    )
    assert result.stdout.strip() == "False", result.stderr


@pytest.mark.parametrize("element, orders, expected", [
    ("C", [], 4), ("C", [1.0], 3), ("C", [2.0], 2), ("C", [3.0], 1), ("C", [1.5, 1.5], 1),
    ("C", [1.5, 1.5, 1.5], 0), ("N", [1.0], 2), ("N", [2.0], 1), ("O", [1.0], 1), ("O", [2.0], 0),
    ("S", [1.0], 1), ("S", [2.0, 2.0], 0), ("P", [], 3), ("F", [], 1), ("H", [1.0], 0),
])
def test_hydrogen_counts(element, orders, expected):
    assert hydrogens_needed(element, orders) == (expected, None)


@pytest.mark.parametrize("element, orders, phrase", [
    ("N", [1.0] * 4, "a charge is likely"),
    ("C", [1.0] * 5, "a charge is likely"),
    ("S", [1.0] * 3, "between the valences"),
    ("P", [1.0] * 4, "between the valences"),
    ("O", [1.5], "cannot be told"),
    ("N", [1.5, 1.5], "pyridine-like"),
    ("Fe", [], "no valence rule"),
])
def test_ambiguous_valences_are_flagged_not_guessed(element, orders, phrase):
    count, problem = hydrogens_needed(element, orders)
    assert count == 0
    assert phrase in problem


def test_ethanol_gets_tetrahedral_hydrogens_at_the_reference_length():
    plan, elements, coordinates, bonds = _apply(
        ["C", "C", "O"], [[0, 0, 0], [1.52, 0, 0], [2.0, 1.35, 0]], [(0, 1, 1.0), (1, 2, 1.0)]
    )
    assert plan.count == 6 and plan.issues == []
    for parent, start in zip([0, 0, 0, 1, 1, 2], range(3, 9)):
        assert measure(DISTANCE, coordinates, (parent, start)) == pytest.approx(
            reference_bond_length(elements[parent], "H"))
        other = 1 if parent != 1 else 0
        assert measure(ANGLE, coordinates, (other, parent, start)) == pytest.approx(109.47, abs=0.1)


def test_ethylene_comes_out_planar():
    _, _, coordinates, _ = _apply(["C", "C"], [[0, 0, 0], [1.34, 0, 0]], [(0, 1, 2.0)])
    assert np.allclose(coordinates[:, 2], 0.0, atol=1e-9)
    assert abs(measure(DIHEDRAL, coordinates, (2, 0, 1, 4))) in (pytest.approx(0.0, abs=1e-6),
                                                                 pytest.approx(180.0, abs=1e-6))


def test_benzene_hydrogens_point_out_of_the_ring():
    plan, _, coordinates, _ = _apply(["C"] * 6, _ring(), [(k, (k + 1) % 6, 1.5) for k in range(6)])
    assert plan.count == 6
    radial = np.linalg.norm(coordinates[6:], axis=1)
    assert np.allclose(radial, 1.39 + reference_bond_length("C", "H"))


def test_water_from_a_lone_oxygen():
    plan, _, coordinates, _ = _apply(["O"], [[0.0, 0.0, 0.0]], [])
    assert plan.count == 2
    assert measure(ANGLE, coordinates, (1, 0, 2)) == pytest.approx(109.47, abs=0.01)


def test_running_twice_adds_nothing_more():
    _, elements, coordinates, bonds = _apply(["C", "N"], [[0, 0, 0], [1.47, 0, 0]], [(0, 1, 1.0)])
    again = plan_hydrogens(elements, coordinates, bonds)
    assert again.count == 0 and again.issues == []


def test_flagged_atoms_are_left_alone_while_others_are_completed():
    plan = plan_hydrogens(["C", "C", "O", "O"], [[0, 0, 0], [1.5, 0, 0], [2.1, 1.1, 0], [2.1, -1.1, 0]],
                          [(0, 1, 1.0), (1, 2, 1.5), (1, 3, 1.5)])
    assert sorted({parent for parent, _ in plan.additions}) == [0]
    assert [issue.split(":")[0] for issue in plan.issues] == ["atom 3 (O)", "atom 4 (O)"]


def test_only_atoms_in_scope_get_hydrogens():
    plan = plan_hydrogens(["C", "C"], [[0, 0, 0], [5.0, 0, 0]], [], scope={1})
    assert {parent for parent, _ in plan.additions} == {1}
