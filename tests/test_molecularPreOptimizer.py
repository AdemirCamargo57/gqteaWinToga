"""Tests for the molecularPreOptimizer force field and geometry optimizer.

Pure numeric tests: no Toga, no OpenGL. The load-bearing test here is the
analytic-gradient check against finite differences -- every other assertion in
this file only means something if the gradients the minimizer is handed are
actually the gradients of the energy it is minimizing.

Geometry assertions use loose tolerances on purpose: this is a pre-optimizer
producing a reasonable starting structure, not a calibrated force field.
"""
import numpy as np
import pytest

from molecularPreOptimizer import (
    PIXELS_PER_ANGSTROM,
    ForceField,
    OptimizationResult,
    build_3d_seed,
    optimize_sketch,
)
from molecularSketch import MoleculeSketch


# ------------------------------------------------------------------
# Builders for the small molecules used below
# ------------------------------------------------------------------
def _diatomic(element_a: str, element_b: str, order: int = 1) -> MoleculeSketch:
    sketch = MoleculeSketch()
    a = sketch.add_atom(element_a, 100.0, 100.0)
    b = sketch.add_atom(element_b, 300.0, 100.0)
    bond = sketch.add_bond(a.atom_id, b.atom_id)
    bond.order = order
    return sketch


def _methane() -> MoleculeSketch:
    sketch = MoleculeSketch()
    carbon = sketch.add_atom("C", 200.0, 200.0)
    for x, y in ((200.0, 140.0), (260.0, 200.0), (200.0, 260.0), (140.0, 200.0)):
        hydrogen = sketch.add_atom("H", x, y)
        sketch.add_bond(carbon.atom_id, hydrogen.atom_id)
    return sketch


def _water() -> MoleculeSketch:
    sketch = MoleculeSketch()
    oxygen = sketch.add_atom("O", 200.0, 200.0)
    for x, y in ((140.0, 260.0), (260.0, 260.0)):
        hydrogen = sketch.add_atom("H", x, y)
        sketch.add_bond(oxygen.atom_id, hydrogen.atom_id)
    return sketch


def _carbon_dioxide() -> MoleculeSketch:
    sketch = MoleculeSketch()
    carbon = sketch.add_atom("C", 200.0, 200.0)
    for x in (140.0, 260.0):
        oxygen = sketch.add_atom("O", x, 210.0)
        bond = sketch.add_bond(carbon.atom_id, oxygen.atom_id)
        bond.order = 2
    return sketch


def _ethene() -> MoleculeSketch:
    sketch = MoleculeSketch()
    c1 = sketch.add_atom("C", 180.0, 200.0)
    c2 = sketch.add_atom("C", 240.0, 200.0)
    double = sketch.add_bond(c1.atom_id, c2.atom_id)
    double.order = 2
    for carbon, positions in (
        (c1, ((140.0, 170.0), (140.0, 230.0))),
        (c2, ((280.0, 170.0), (280.0, 230.0))),
    ):
        for x, y in positions:
            hydrogen = sketch.add_atom("H", x, y)
            sketch.add_bond(carbon.atom_id, hydrogen.atom_id)
    return sketch


def _butadiene() -> MoleculeSketch:
    """C1=C2-C3=C4 with hydrogens: exercises bonds, angles, torsions and
    non-bonded pairs (H on C1 is four bonds from C4) in one structure."""
    sketch = MoleculeSketch()
    carbons = [sketch.add_atom("C", 120.0 + 60.0 * i, 200.0 + 20.0 * (i % 2)) for i in range(4)]
    orders = (2, 1, 2)
    for index, order in enumerate(orders):
        bond = sketch.add_bond(carbons[index].atom_id, carbons[index + 1].atom_id)
        bond.order = order
    hydrogen_positions = {
        0: ((80.0, 170.0), (80.0, 240.0)),
        1: ((180.0, 260.0),),
        2: ((240.0, 150.0),),
        3: ((320.0, 170.0), (320.0, 240.0)),
    }
    for carbon_index, positions in hydrogen_positions.items():
        for x, y in positions:
            hydrogen = sketch.add_atom("H", x, y)
            sketch.add_bond(carbons[carbon_index].atom_id, hydrogen.atom_id)
    return sketch


def _bond_length(result: OptimizationResult, i: int, j: int) -> float:
    return float(np.linalg.norm(result.coordinates[i] - result.coordinates[j]))


def _angle_degrees(result: OptimizationResult, i: int, j: int, k: int) -> float:
    v1 = result.coordinates[i] - result.coordinates[j]
    v2 = result.coordinates[k] - result.coordinates[j]
    cosine = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def _max_out_of_plane(coordinates: np.ndarray) -> float:
    """Largest distance of any atom from the best-fit plane through them all."""
    centred = coordinates - coordinates.mean(axis=0)
    _, _, vectors = np.linalg.svd(centred)
    normal = vectors[-1]
    return float(np.max(np.abs(centred @ normal)))


# ------------------------------------------------------------------
# Seeding
# ------------------------------------------------------------------
def test_seed_converts_canvas_pixels_to_angstrom():
    sketch = _diatomic("C", "C")

    _, coordinates = build_3d_seed(sketch)

    separation = np.linalg.norm(coordinates[0][:2] - coordinates[1][:2])
    assert separation == pytest.approx(200.0 / PIXELS_PER_ANGSTROM)


def test_seed_flips_the_canvas_y_axis():
    """Canvas y grows downward; the 3D structure must not come out mirrored."""
    sketch = MoleculeSketch()
    top = sketch.add_atom("C", 100.0, 50.0)
    bottom = sketch.add_atom("O", 100.0, 250.0)
    sketch.add_bond(top.atom_id, bottom.atom_id)

    _, coordinates = build_3d_seed(sketch)

    assert coordinates[0][1] > coordinates[1][1]


def test_seed_is_not_perfectly_planar():
    """A flat start leaves several angle gradients exactly zero."""
    _, coordinates = build_3d_seed(_methane())

    assert np.ptp(coordinates[:, 2]) > 0.0


def test_seed_is_deterministic():
    first = build_3d_seed(_methane())[1]
    second = build_3d_seed(_methane())[1]

    assert np.array_equal(first, second)


def test_seed_reports_elements_in_sketch_order():
    elements, _ = build_3d_seed(_water())

    assert elements == ["O", "H", "H"]


# ------------------------------------------------------------------
# The force field and its gradients
# ------------------------------------------------------------------
def test_analytic_gradient_matches_finite_differences():
    """The one test that makes every other result trustworthy."""
    sketch = _butadiene()
    elements, coordinates = build_3d_seed(sketch)
    field = ForceField.from_sketch(sketch, elements)
    flat = coordinates.ravel().copy()

    _, analytic = field.energy_and_gradient(flat)

    step = 1e-6
    numeric = np.zeros_like(flat)
    for index in range(flat.size):
        forward = flat.copy()
        backward = flat.copy()
        forward[index] += step
        backward[index] -= step
        numeric[index] = (
            field.energy_and_gradient(forward)[0] - field.energy_and_gradient(backward)[0]
        ) / (2.0 * step)

    assert np.allclose(analytic, numeric, rtol=1e-4, atol=1e-6)


def test_the_force_field_includes_every_term_for_butadiene():
    sketch = _butadiene()
    elements, _ = build_3d_seed(sketch)

    field = ForceField.from_sketch(sketch, elements)

    assert field.bond_terms, "expected bond stretch terms"
    assert field.angle_terms, "expected angle bend terms"
    assert field.torsion_terms, "expected a torsion term on each double bond"
    assert field.nonbonded_pairs, "expected non-bonded pairs at least four bonds apart"


def test_atoms_three_bonds_apart_are_not_counted_as_non_bonded():
    sketch = _ethene()
    elements, _ = build_3d_seed(sketch)

    field = ForceField.from_sketch(sketch, elements)

    assert field.nonbonded_pairs == []


def test_the_energy_is_translation_invariant():
    sketch = _methane()
    elements, coordinates = build_3d_seed(sketch)
    field = ForceField.from_sketch(sketch, elements)

    here, _ = field.energy_and_gradient(coordinates.ravel())
    there, _ = field.energy_and_gradient((coordinates + np.array([3.0, -2.0, 1.0])).ravel())

    assert here == pytest.approx(there)


# ------------------------------------------------------------------
# Optimized geometries
# ------------------------------------------------------------------
def test_a_single_bond_relaxes_to_the_sum_of_the_covalent_radii():
    result = optimize_sketch(_diatomic("C", "C"))

    assert result.success
    assert _bond_length(result, 0, 1) == pytest.approx(1.52, abs=0.01)


def test_higher_bond_orders_give_shorter_bonds():
    single = _bond_length(optimize_sketch(_diatomic("C", "C", order=1)), 0, 1)
    double = _bond_length(optimize_sketch(_diatomic("C", "C", order=2)), 0, 1)
    triple = _bond_length(optimize_sketch(_diatomic("C", "C", order=3)), 0, 1)

    assert triple < double < single


def test_methane_relaxes_to_four_equal_bonds():
    result = optimize_sketch(_methane())

    lengths = [_bond_length(result, 0, index) for index in range(1, 5)]
    assert result.success
    assert max(lengths) - min(lengths) < 0.01
    assert lengths[0] == pytest.approx(1.07, abs=0.02)


def test_methane_relaxes_to_tetrahedral_angles():
    result = optimize_sketch(_methane())

    angles = [
        _angle_degrees(result, i, 0, j)
        for i in range(1, 5)
        for j in range(i + 1, 5)
    ]
    assert all(angle == pytest.approx(109.47, abs=2.0) for angle in angles)


def test_water_is_bent_not_linear():
    """Lone pairs count: two bonds on oxygen still give a tetrahedral angle."""
    result = optimize_sketch(_water())

    assert _angle_degrees(result, 1, 0, 2) == pytest.approx(109.47, abs=3.0)


def test_carbon_dioxide_is_linear():
    """Two double bonds and no lone pairs on carbon give a steric number of 2."""
    result = optimize_sketch(_carbon_dioxide())

    assert _angle_degrees(result, 1, 0, 2) == pytest.approx(180.0, abs=3.0)


def test_ethene_relaxes_flat():
    result = optimize_sketch(_ethene())

    assert result.success
    assert _max_out_of_plane(result.coordinates) < 0.05


def test_optimization_lowers_the_energy_from_the_seed():
    sketch = _butadiene()
    elements, coordinates = build_3d_seed(sketch)
    field = ForceField.from_sketch(sketch, elements)
    seed_energy, _ = field.energy_and_gradient(coordinates.ravel())

    result = optimize_sketch(sketch)

    assert result.energy < seed_energy


def test_optimization_is_deterministic():
    first = optimize_sketch(_methane())
    second = optimize_sketch(_methane())

    assert np.array_equal(first.coordinates, second.coordinates)


def test_a_lone_atom_optimizes_trivially():
    sketch = MoleculeSketch()
    sketch.add_atom("Ne", 100.0, 100.0)

    result = optimize_sketch(sketch)

    assert result.success
    assert result.coordinates.shape == (1, 3)


# ------------------------------------------------------------------
# Results, failures and reporting
# ------------------------------------------------------------------
def test_the_result_converts_to_the_frame_shape_the_viewer_uses():
    result = optimize_sketch(_water())

    frame = result.to_frame()

    assert len(frame) == 3
    element, position = frame[0]
    assert element == "O"
    assert len(position) == 3
    assert all(isinstance(value, float) for value in position)


def test_hitting_the_iteration_limit_reports_failure_without_raising():
    result = optimize_sketch(_butadiene(), max_iterations=1)

    assert result.success is False
    assert result.message
    assert result.coordinates.shape[0] == 10


def test_partial_coordinates_are_returned_even_when_optimization_fails():
    result = optimize_sketch(_butadiene(), max_iterations=1)

    assert np.all(np.isfinite(result.coordinates))


def test_an_empty_sketch_is_rejected_with_a_clear_message():
    with pytest.raises(ValueError, match="no atoms"):
        optimize_sketch(MoleculeSketch())


def test_a_structure_with_a_validation_error_is_rejected():
    sketch = _water()
    sketch.add_atom("C", 500.0, 500.0)

    with pytest.raises(ValueError, match="not bonded"):
        optimize_sketch(sketch)


def test_validation_warnings_are_carried_into_the_result():
    sketch = MoleculeSketch()
    carbon = sketch.add_atom("C", 200.0, 200.0)
    for index in range(5):
        hydrogen = sketch.add_atom("H", 140.0 + 30.0 * index, 260.0)
        sketch.add_bond(carbon.atom_id, hydrogen.atom_id)

    result = optimize_sketch(sketch)

    assert any("valence" in warning.lower() for warning in result.warnings)


def test_the_optimizer_reports_its_gradient_norm_and_iteration_count():
    result = optimize_sketch(_methane())

    assert result.gradient_norm >= 0.0
    assert result.iterations > 0


def test_the_optimizer_never_imports_toga():
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-c", "import molecularPreOptimizer, sys; print('toga' in sys.modules)"],
        capture_output=True,
        text=True,
    )

    assert result.stdout.strip() == "False", result.stderr
