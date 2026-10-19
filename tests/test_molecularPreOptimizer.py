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

import molecularPreOptimizer
from molecularGeometry import constraint_residuals, measure
from molecularPreOptimizer import (
    PIXELS_PER_ANGSTROM,
    ForceField,
    OptimizationResult,
    angstrom_to_sketch,
    build_3d_seed,
    ideal_angle_degrees,
    optimize_sketch,
    reference_bond_length,
    sketch_to_angstrom,
)
from molecularSketch import RESONANCE_ORDER, MoleculeSketch


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


# ------------------------------------------------------------------
# Resonance bonds
# ------------------------------------------------------------------
def _resonance_benzene() -> MoleculeSketch:
    """Six carbons in a ring, every bond delocalized."""
    import math

    sketch = MoleculeSketch()
    carbons = []
    for index in range(6):
        angle = math.radians(60 * index)
        carbons.append(
            sketch.add_atom("C", 200.0 + 60.0 * math.cos(angle), 200.0 + 60.0 * math.sin(angle))
        )
    for index in range(6):
        bond = sketch.add_bond(carbons[index].atom_id, carbons[(index + 1) % 6].atom_id)
        bond.order = RESONANCE_ORDER
    return sketch


def test_a_resonance_bond_lands_between_a_single_and_a_double():
    single = _bond_length(optimize_sketch(_diatomic("C", "C", order=1)), 0, 1)
    resonance = _bond_length(optimize_sketch(_diatomic("C", "C", order=RESONANCE_ORDER)), 0, 1)
    double = _bond_length(optimize_sketch(_diatomic("C", "C", order=2)), 0, 1)

    assert double < resonance < single


def test_an_aromatic_carbon_carbon_bond_is_about_one_point_four_angstrom():
    result = optimize_sketch(_diatomic("C", "C", order=RESONANCE_ORDER))

    assert _bond_length(result, 0, 1) == pytest.approx(1.40, abs=0.02)


def test_resonance_bonds_get_a_planarity_torsion():
    """Without this an aromatic ring is free to pucker."""
    sketch = _resonance_benzene()
    elements, _ = build_3d_seed(sketch)

    field = ForceField.from_sketch(sketch, elements)

    assert field.torsion_terms


def test_a_resonance_benzene_relaxes_flat_with_six_equal_bonds():
    result = optimize_sketch(_resonance_benzene())

    lengths = [_bond_length(result, index, (index + 1) % 6) for index in range(6)]
    assert result.success
    assert max(lengths) - min(lengths) < 0.01
    assert lengths[0] == pytest.approx(1.40, abs=0.03)
    assert _max_out_of_plane(result.coordinates) < 0.05


def test_an_aromatic_ring_keeps_its_one_hundred_and_twenty_degree_angles():
    result = optimize_sketch(_resonance_benzene())

    angles = [_angle_degrees(result, (i - 1) % 6, i, (i + 1) % 6) for i in range(6)]
    assert all(angle == pytest.approx(120.0, abs=2.0) for angle in angles)


def test_a_pyridine_nitrogen_keeps_a_trigonal_angle():
    """Two resonance bonds and one lone pair is steric number 3, not 2."""
    assert ideal_angle_degrees("N", connections=2, bond_order_sum=3.0) == pytest.approx(120.0)


def test_an_aromatic_carbon_is_trigonal():
    assert ideal_angle_degrees("C", connections=3, bond_order_sum=4.0) == pytest.approx(120.0)


def test_the_gradient_is_still_exact_with_resonance_bonds():
    sketch = _resonance_benzene()
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


# ------------------------------------------------------------------
# Recalibrated reference lengths
# ------------------------------------------------------------------
def test_reference_lengths_match_experimental_carbon_carbon_bonds():
    lengths = {
        order: _bond_length(optimize_sketch(_diatomic("C", "C", order=order)), 0, 1)
        for order in (1, RESONANCE_ORDER, 2, 3)
    }

    assert lengths[1] == pytest.approx(1.52, abs=0.02)
    assert lengths[RESONANCE_ORDER] == pytest.approx(1.40, abs=0.02)
    assert lengths[2] == pytest.approx(1.34, abs=0.02)
    assert lengths[3] == pytest.approx(1.22, abs=0.02)


# ------------------------------------------------------------------
# 3D sketches and writing the geometry back onto the canvas
# ------------------------------------------------------------------
def _signed_volume(points: np.ndarray, a: int, b: int, c: int, d: int) -> float:
    """Triple product of the three arms from ``a`` -- its sign is the handedness."""
    return float(np.dot(points[b] - points[a], np.cross(points[c] - points[a], points[d] - points[a])))


def _chiral_centre() -> MoleculeSketch:
    """CHFClBr drawn as a real 3D tetrahedron (sketch units, z away from the viewer)."""
    sketch = MoleculeSketch()
    carbon = sketch.add_atom("C", 200.0, 200.0, 0.0)
    for element, (x, y, z) in (
        ("H", (200.0, 156.0, -15.0)),
        ("F", (241.0, 220.0, -15.0)),
        ("Cl", (159.0, 220.0, -15.0)),
        ("Br", (200.0, 200.0, 50.0)),
    ):
        substituent = sketch.add_atom(element, x, y, z)
        sketch.add_bond(carbon.atom_id, substituent.atom_id)
    return sketch


def test_the_sketch_to_angstrom_conversion_round_trips():
    points = np.array([[10.0, -20.0, 30.0], [400.0, 5.0, -60.0]])

    assert np.allclose(angstrom_to_sketch(sketch_to_angstrom(points)), points)


def test_the_conversion_is_a_rotation_not_a_mirror():
    """Flipping y alone would turn a 3D structure into its mirror image."""
    points = _chiral_centre().positions()

    converted = sketch_to_angstrom(points)

    assert np.sign(_signed_volume(converted, 0, 1, 2, 3)) == np.sign(
        _signed_volume(points, 0, 1, 2, 3)
    )


def test_a_three_dimensional_sketch_is_seeded_without_any_jitter():
    sketch = _chiral_centre()

    _, coordinates = build_3d_seed(sketch)

    expected = sketch_to_angstrom(sketch.positions())
    assert np.allclose(coordinates, expected - expected.mean(axis=0))


def test_a_flat_sketch_tilted_out_of_the_screen_is_still_lifted_off_its_plane():
    """A planar drawing made after rotating the view has non-zero z, but it is
    still planar and still needs the out-of-plane nudge."""
    sketch = _methane()
    tilt = np.radians(35.0)
    rotation = np.array(
        [[1.0, 0.0, 0.0], [0.0, np.cos(tilt), -np.sin(tilt)], [0.0, np.sin(tilt), np.cos(tilt)]]
    )
    sketch.update_positions(sketch.positions() @ rotation.T)

    _, coordinates = build_3d_seed(sketch)

    assert _max_out_of_plane(coordinates) > 0.01


def test_reference_bond_length_is_the_scaled_sum_of_covalent_radii():
    from help import AtomicData

    single = AtomicData.covalent_radii["C"] + AtomicData.covalent_radii["O"]

    assert reference_bond_length("C", "O") == pytest.approx(single)
    assert reference_bond_length("C", "O", 2) == pytest.approx(single * 0.88)


def test_reference_bond_length_is_none_for_an_element_without_a_radius():
    assert reference_bond_length("C", "Xx") is None


def test_the_result_carries_sketch_positions_for_every_atom():
    sketch = _methane()

    result = optimize_sketch(sketch)

    assert result.sketch_positions.shape == (5, 3)
    assert np.all(np.isfinite(result.sketch_positions))


def test_written_back_positions_keep_the_drawing_where_it_was():
    sketch = _methane()

    result = optimize_sketch(sketch)

    assert np.allclose(result.sketch_positions.mean(axis=0), sketch.positions().mean(axis=0))


def test_sketch_positions_and_angstrom_coordinates_describe_the_same_structure():
    result = optimize_sketch(_ethene())

    converted = sketch_to_angstrom(result.sketch_positions)

    assert np.allclose(converted - converted.mean(axis=0), result.coordinates)


def test_the_result_is_aligned_onto_the_drawing():
    """The optimized structure settles in place instead of spinning: it sits
    closer to what was drawn than the same structure in any rotated pose."""
    sketch = _ethene()
    drawn = sketch_to_angstrom(sketch.positions())
    drawn -= drawn.mean(axis=0)

    result = optimize_sketch(sketch)

    aligned_rmsd = np.sqrt(np.mean(np.sum((result.coordinates - drawn) ** 2, axis=1)))
    turned = result.coordinates @ np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]).T
    turned_rmsd = np.sqrt(np.mean(np.sum((turned - drawn) ** 2, axis=1)))
    assert aligned_rmsd < turned_rmsd


def test_reoptimizing_a_written_back_structure_barely_moves_it():
    sketch = _methane()
    sketch.update_positions(optimize_sketch(sketch).sketch_positions)
    before = sketch.positions()

    again = optimize_sketch(sketch)

    assert np.max(np.abs(again.sketch_positions - before)) < 0.05 * PIXELS_PER_ANGSTROM


def test_optimizing_a_chiral_centre_keeps_its_handedness():
    sketch = _chiral_centre()
    drawn = sketch.positions()

    result = optimize_sketch(sketch)

    assert np.sign(_signed_volume(result.sketch_positions, 0, 1, 2, 3)) == np.sign(
        _signed_volume(drawn, 0, 1, 2, 3)
    )


# ----------------------------------------------------------------------
# Frozen parameters
# ----------------------------------------------------------------------
def _butane_skeleton() -> MoleculeSketch:
    sketch = MoleculeSketch()
    ids = [sketch.add_atom("C", 100.0 + 60.0 * k, 100.0 + (25.0 if k % 2 else 0.0)).atom_id for k in range(4)]
    for first, second in zip(ids, ids[1:]):
        sketch.add_bond(first, second)
    return sketch


def test_frozen_parameters_hold_at_every_evaluated_geometry(monkeypatch):
    sketch = _butane_skeleton()
    ids = sketch.ordered_atom_ids()
    sketch.set_constraint(ids, 60.0)
    sketch.set_constraint(ids[:2], 1.60)
    sketch.set_constraint(ids[:3], 100.0)
    evaluated = []
    original = molecularPreOptimizer.ForceField.energy_and_gradient

    def recording(self, flat):
        evaluated.append(np.array(flat).reshape(-1, 3))
        return original(self, flat)

    monkeypatch.setattr(molecularPreOptimizer.ForceField, "energy_and_gradient", recording)
    result = optimize_sketch(sketch)

    constraints = molecularPreOptimizer.sketch_constraints(sketch)
    assert result.success and result.constraint_count == 3
    assert result.constraint_deviation < 1e-8
    assert max(np.abs(constraint_residuals(x, constraints)[0]).max() for x in evaluated) < 1e-9
    assert measure("dihedral", result.coordinates, (0, 1, 2, 3)) == pytest.approx(60.0, abs=1e-6)
    assert measure("distance", result.coordinates, (0, 1)) == pytest.approx(1.60, abs=1e-8)
    assert measure("angle", result.coordinates, (0, 1, 2)) == pytest.approx(100.0, abs=1e-6)
    # Unconstrained parts still relax: the free C-C bond goes to its reference.
    assert measure("distance", result.coordinates, (2, 3)) == pytest.approx(
        reference_bond_length("C", "C"), abs=1e-3)


def test_a_constrained_run_is_also_a_minimum():
    sketch = _butane_skeleton()
    ids = sketch.ordered_atom_ids()
    sketch.set_constraint(ids, 180.0)
    result = optimize_sketch(sketch)
    assert result.success
    assert result.gradient_norm < 1e-3  # the tangent gradient, as for L-BFGS-B


def test_without_frozen_parameters_the_result_is_unchanged():
    first, second = _butane_skeleton(), _butane_skeleton()
    ids = second.ordered_atom_ids()
    second.set_constraint(ids[:2], 1.0)
    second.clear_constraints()
    assert np.array_equal(optimize_sketch(first).coordinates, optimize_sketch(second).coordinates)
    assert optimize_sketch(first).constraint_count == 0


def test_a_linear_frozen_angle_is_reported_as_unsupported():
    sketch = _butane_skeleton()
    ids = sketch.ordered_atom_ids()
    sketch.set_constraint(ids[:3], 180.0)
    with pytest.raises(ValueError, match="cannot be held by the pre-optimizer: angle 1-2-3"):
        optimize_sketch(sketch)


def test_contradictory_frozen_parameters_are_reported():
    sketch = _butane_skeleton()
    ids = sketch.ordered_atom_ids()
    sketch.add_bond(ids[0], ids[2])
    sketch.set_constraint(ids[:2], 1.0)
    sketch.set_constraint(ids[1:3], 1.0)
    sketch.set_constraint([ids[0], ids[2]], 3.0)
    with pytest.raises(ValueError, match="cannot all hold"):
        optimize_sketch(sketch)
