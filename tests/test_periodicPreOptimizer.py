import itertools
import pathlib
import subprocess
import sys

import numpy as np
import pytest

import periodicPreOptimizer as ppo
from help import AtomicData

SRC = pathlib.Path(__file__).resolve().parents[1]

WATER_ELEMENTS = ["O", "H", "H"]
WATER = np.array([[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]])


def radius(e):
    return AtomicData.vdw_radii.get(e, 1.7)


def mass(e):
    return AtomicData.atomic_masses.get(e, 0.0)


def water_system(centres, rotations=None, kinds=None):
    molecules = []
    for index, centre in enumerate(centres):
        coords = WATER.copy()
        if rotations is not None:
            coords = coords @ rotations[index].T
        kind = kinds[index] if kinds else "solvent"
        molecules.append((kind, WATER_ELEMENTS, coords + np.asarray(centre, dtype=float)))
    return ppo.PackedSystem.from_molecules(molecules, radius, mass)


def random_rotations(n, seed=1):
    from scipy.spatial.transform import Rotation
    return Rotation.random(n, random_state=seed).as_matrix()


def finite_difference(fn, x, h=1e-6):
    grad = np.zeros_like(x)
    for k in range(len(x)):
        step = np.zeros_like(x)
        step[k] = h
        grad[k] = (fn(x + step)[0] - fn(x - step)[0]) / (2 * h)
    return grad


def crowded_system(cell_length=9.0, n=12, seed=3):
    rng = np.random.default_rng(seed)
    centres = rng.uniform(0.0, cell_length, size=(n, 3))
    return water_system(centres, random_rotations(n, seed)), ppo.PeriodicCell([cell_length] * 3)


# ----------------------------------------------------------------------
# Cell and wrapping
# ----------------------------------------------------------------------
def test_minimum_image_matches_brute_force_image_scan():
    cell = ppo.PeriodicCell([7.0, 9.0, 11.0])
    rng = np.random.default_rng(0)
    for _ in range(200):
        delta = rng.uniform(-30, 30, size=3)
        images = [delta + cell.lengths * np.array(n) for n in itertools.product(range(-6, 7), repeat=3)]
        best = min(images, key=lambda v: float(v @ v))
        assert np.linalg.norm(cell.minimum_image(delta[None])[0]) == pytest.approx(np.linalg.norm(best), abs=1e-12)


def test_wrap_points_stay_in_cell():
    cell = ppo.PeriodicCell([5.0, 5.0, 5.0])
    pts = np.array([[-1e-17, 5.0, 10.0], [-5.0, 4.9999999999, 123.4]])
    wrapped = cell.wrap_points(pts)
    assert np.all(wrapped >= 0.0) and np.all(wrapped < cell.lengths)


def test_wrap_molecules_keeps_straddling_molecule_whole():
    cell = ppo.PeriodicCell([10.0, 10.0, 10.0])
    # A water whose O sits just inside x = 10 and whose H atoms poke past it.
    system = water_system([[9.9, 5.0, 5.0], [-3.0, 14.0, 25.5]],
                          rotations=[np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1.0]]), np.eye(3)])
    wrapped = ppo.wrap_molecules(system.coordinates, system.masses, system.molecule_slices, cell)
    centres = ppo.molecule_centres(wrapped, system.masses, system.molecule_slices)
    assert np.all(centres >= 0.0) and np.all(centres < 10.0)
    for start, stop in system.molecule_slices:
        before = system.coordinates[start:stop]
        after = wrapped[start:stop]
        shift = after - before
        assert np.allclose(shift, shift[0])                      # one rigid translation
        assert np.allclose(shift[0] / 10.0, np.round(shift[0] / 10.0))  # by whole cell vectors
    # The straddling molecule was not split: its H atoms are still beyond x = 10.
    assert wrapped[:3, 0].max() > 10.0


def test_make_whole_undoes_per_atom_wrapping():
    cell = ppo.PeriodicCell([8.0, 8.0, 8.0])
    system = water_system([[7.9, 0.05, 4.0]], rotations=[random_rotations(1)[0]])
    broken = cell.wrap_points(system.coordinates)
    whole = ppo.make_whole(broken, system.molecule_slices, cell)
    diffs = whole - system.coordinates
    assert np.allclose(diffs, diffs[0])


def test_packed_system_validates_slices():
    with pytest.raises(ValueError):
        ppo.PackedSystem(["O", "H"], np.zeros((2, 3)), [(0, 1)], ["a"], np.ones(2), np.ones(2))


# ----------------------------------------------------------------------
# Penalty
# ----------------------------------------------------------------------
def test_penalty_is_zero_without_clash_and_ignores_intramolecular_pairs():
    system = water_system([[2.0, 2.0, 2.0], [7.0, 7.0, 7.0]])
    penalty = ppo.OverlapPenalty(system, ppo.PeriodicCell([12.0] * 3), vdw_scale=0.8)
    energy, gradient = penalty.energy_and_gradient(system.coordinates)
    assert energy == 0.0 and not np.any(gradient)
    assert penalty.report(system.coordinates)["clashes"] == 0


def test_penalty_sees_clash_across_the_boundary():
    system = water_system([[0.3, 5.0, 5.0], [9.7, 5.0, 5.0]])
    periodic = ppo.OverlapPenalty(system, ppo.PeriodicCell([10.0] * 3), 0.8)
    assert periodic.energy_and_gradient(system.coordinates)[0] > 0.0
    assert periodic.report(system.coordinates)["clashes"] > 0
    flat = ppo.OverlapPenalty(system, ppo.PeriodicCell([10.0] * 3, periodic=False), 0.8)
    assert flat.report(system.coordinates)["clashes"] == 0


def test_penalty_is_invariant_under_a_whole_box_shift():
    system, cell = crowded_system()
    penalty = ppo.OverlapPenalty(system, cell, 0.8)
    e0 = penalty.energy_and_gradient(system.coordinates)[0]
    for shift in ([3.3, -1.2, 40.0], [-9.0, 9.0, 0.5]):
        e1 = penalty.energy_and_gradient(system.coordinates + np.array(shift))[0]
        assert e1 == pytest.approx(e0, rel=1e-10)
    assert e0 > 0.0


def test_solute_clearance_only_applies_to_solute_solvent_pairs():
    system = water_system([[2.0, 2.0, 2.0], [5.0, 2.0, 2.0]], kinds=["solute", "solvent"])
    cell = ppo.PeriodicCell([14.0] * 3)
    plain = ppo.OverlapPenalty(system, cell, 0.8)
    padded = ppo.OverlapPenalty(system, cell, 0.8, solute_clearance=1.5)
    assert plain.report(system.coordinates)["clashes"] == 0
    assert padded.report(system.coordinates)["clashes"] > 0


def test_contact_range_must_respect_half_box():
    system = water_system([[1.0, 1.0, 1.0]])
    with pytest.raises(ValueError, match="half the"):
        ppo.OverlapPenalty(system, ppo.PeriodicCell([4.0, 4.0, 4.0]), 0.9)


def test_penalty_gradient_matches_finite_differences():
    system, cell = crowded_system()
    penalty = ppo.OverlapPenalty(system, cell, 0.8, solute_clearance=0.0)
    x = system.coordinates.ravel()
    fn = lambda flat: (penalty.energy_and_gradient(flat.reshape(-1, 3))[0], None)
    analytic = penalty.energy_and_gradient(system.coordinates)[1].ravel()
    assert np.allclose(analytic, finite_difference(fn, x), rtol=1e-5, atol=1e-5)


def test_wall_gradient_matches_finite_differences():
    system = water_system([[0.5, 5.0, 9.6], [5.0, 5.0, 5.0]])
    penalty = ppo.OverlapPenalty(system, ppo.PeriodicCell([10.0] * 3, periodic=False), 0.8, wall_padding=0.5)
    fn = lambda flat: (penalty.energy_and_gradient(flat.reshape(-1, 3))[0], None)
    analytic = penalty.energy_and_gradient(system.coordinates)[1].ravel()
    assert penalty.energy_and_gradient(system.coordinates)[0] > 0.0
    assert np.allclose(analytic, finite_difference(fn, system.coordinates.ravel()), rtol=1e-5, atol=1e-5)


# ----------------------------------------------------------------------
# Models
# ----------------------------------------------------------------------
@pytest.mark.parametrize("w_scale", [0.0, 1e-5, 0.7, 2.5])
def test_rigid_gradient_matches_finite_differences(w_scale):
    system, cell = crowded_system()
    model = ppo.RigidBodyModel(system, ppo.OverlapPenalty(system, cell, 0.8))
    params = model.initial_parameters()
    rng = np.random.default_rng(5)
    params[3 * model.count:] = w_scale * rng.normal(size=3 * model.count)
    analytic = model.energy_and_gradient(params)[1]
    numeric = finite_difference(model.energy_and_gradient, params)
    assert np.allclose(analytic, numeric, rtol=1e-5, atol=1e-5)


def test_rigid_coordinates_preserve_internal_geometry():
    system, cell = crowded_system()
    model = ppo.RigidBodyModel(system, ppo.OverlapPenalty(system, cell, 0.8))
    params = model.initial_parameters() + np.random.default_rng(2).normal(size=6 * model.count)
    x = model.coordinates(params)
    assert ppo._internal_deviation(system.coordinates, x, system.molecule_slices) < 1e-12


def test_flexible_gradient_matches_finite_differences():
    system, cell = crowded_system(n=6)
    model = ppo.FlexibleModel(system, ppo.OverlapPenalty(system, cell, 0.8))
    params = model.initial_parameters() + np.random.default_rng(9).normal(scale=0.05, size=3 * system.atom_count)
    analytic = model.energy_and_gradient(params)[1]
    numeric = finite_difference(model.energy_and_gradient, params)
    assert np.allclose(analytic, numeric, rtol=1e-5, atol=1e-4)


def test_flexible_input_geometry_is_a_bonded_minimum():
    system = water_system([[3.0, 3.0, 3.0], [9.0, 9.0, 9.0]])
    model = ppo.FlexibleModel(system, ppo.OverlapPenalty(system, ppo.PeriodicCell([14.0] * 3), 0.8))
    _, gradient = model.energy_and_gradient(model.initial_parameters())
    assert np.max(np.abs(gradient)) < 1e-8


# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------
def _molecule_centres_inside(result, system, cell):
    centres = ppo.molecule_centres(result.coordinates, system.masses, system.molecule_slices)
    return np.all(centres >= 0.0) and np.all(centres < cell.lengths)


def test_rigid_pre_optimization_removes_clashes_and_keeps_molecules_whole():
    system, cell = crowded_system(cell_length=10.0, n=12)
    result = ppo.pre_optimize(system, cell, 0.8, mode=ppo.MODE_RIGID)
    assert result.report_initial["clashes"] > 0
    assert result.report_final["clashes"] == 0
    assert result.energy_final < 1e-6 * max(1.0, result.energy_initial)
    assert result.max_internal_deviation < 1e-10
    assert result.connectivity_preserved is True
    assert _molecule_centres_inside(result, system, cell)


def test_flexible_pre_optimization_keeps_connectivity():
    system, cell = crowded_system(cell_length=10.0, n=12)
    result = ppo.pre_optimize(system, cell, 0.8, mode=ppo.MODE_FLEXIBLE)
    assert result.report_final["clashes"] < result.report_initial["clashes"]
    assert result.connectivity_preserved is True
    assert result.max_internal_deviation < 0.1
    assert _molecule_centres_inside(result, system, cell)


def test_non_periodic_run_stays_inside_the_walls():
    system = water_system([[1.5, 5.0, 5.0], [2.3, 5.0, 5.0], [8.5, 5.0, 5.0]])
    cell = ppo.PeriodicCell([10.0] * 3, periodic=False)
    result = ppo.pre_optimize(system, cell, 0.8, wall_padding=0.5)
    assert result.report_final["clashes"] == 0
    assert np.all(result.coordinates > 0.5 - 0.05) and np.all(result.coordinates < 9.5 + 0.05)


def test_clash_free_box_is_returned_unchanged():
    system = water_system([[2.0, 2.0, 2.0], [7.0, 7.0, 7.0]])
    cell = ppo.PeriodicCell([12.0] * 3)
    result = ppo.pre_optimize(system, cell, 0.8)
    assert result.success and result.iterations == 0
    assert np.allclose(result.coordinates, system.coordinates)


def test_callback_and_cancel():
    system, cell = crowded_system(cell_length=10.0, n=12)
    seen = []

    def callback(coords, iteration, energy):
        seen.append(iteration)
        raise RuntimeError("reporting failures must be ignored")

    result = ppo.pre_optimize(system, cell, 0.8, callback=callback, should_stop=lambda: len(seen) >= 3)
    assert result.cancelled and not result.success
    assert seen == [1, 2, 3]
    assert "stopped" in result.message


def test_bad_mode_and_iterations():
    system, cell = crowded_system()
    with pytest.raises(ValueError):
        ppo.pre_optimize(system, cell, 0.8, mode="soft")
    with pytest.raises(ValueError):
        ppo.pre_optimize(system, cell, 0.8, max_iterations=0)


def test_summary_lines():
    system, cell = crowded_system(cell_length=10.0, n=12)
    text = "".join(ppo.pre_optimize(system, cell, 0.8).summary_lines())
    assert "rigid bodies" in text and "periodic" in text
    assert "Clashing atom pairs" in text and "Molecular connectivity: preserved" in text


# ----------------------------------------------------------------------
# Shared-wall support: slab restraint and per-atom walls
# ----------------------------------------------------------------------
def test_slab_restraint_gradient_matches_finite_differences():
    system, cell = crowded_system()
    bounds = np.tile([3.0, 5.0], (len(system.molecule_slices), 1))
    restraint = ppo.SlabRestraint(system, bounds, axis=2)
    energy, analytic = restraint.energy_and_gradient(system.coordinates)
    assert energy > 0.0
    fn = lambda flat: (restraint.energy_and_gradient(flat.reshape(-1, 3))[0], None)
    assert np.allclose(analytic.ravel(), finite_difference(fn, system.coordinates.ravel()), rtol=1e-6, atol=1e-6)


def test_composite_rigid_gradient_with_slab():
    system, cell = crowded_system()
    bounds = np.tile([2.0, 6.0], (len(system.molecule_slices), 1))
    term = ppo.CompositeTerm([ppo.OverlapPenalty(system, cell, 0.8), ppo.SlabRestraint(system, bounds)])
    model = ppo.RigidBodyModel(system, term)
    params = model.initial_parameters() + np.random.default_rng(4).normal(scale=0.3, size=6 * model.count)
    numeric = finite_difference(model.energy_and_gradient, params)
    assert np.allclose(model.energy_and_gradient(params)[1], numeric, rtol=1e-5, atol=1e-5)


def test_slab_restraint_keeps_two_regions_apart():
    # Two crowded layers of water across a z interface at 6 A in a 12 A periodic cell.
    rng = np.random.default_rng(8)
    lower = np.column_stack([rng.uniform(0, 12, 10), rng.uniform(0, 12, 10), rng.uniform(4.5, 5.9, 10)])
    upper = np.column_stack([rng.uniform(0, 12, 10), rng.uniform(0, 12, 10), rng.uniform(6.1, 7.5, 10)])
    system = water_system(np.vstack([lower, upper]), random_rotations(20, 8))
    cell = ppo.PeriodicCell([12.0, 12.0, 12.0])
    bounds = [(0.0, 6.0)] * 10 + [(6.0, 12.0)] * 10
    result = ppo.pre_optimize(system, cell, 0.9, slab_bounds=bounds)
    centres = ppo.molecule_centres(result.coordinates, system.masses, system.molecule_slices)
    assert result.report_final["clashes"] == 0
    assert np.all(centres[:10, 2] <= 6.0 + 1e-3) and np.all(centres[10:, 2] >= 6.0 - 1e-3)
    assert not any("outside their region" in note for note in result.notes)


def test_per_atom_wall_bounds():
    system = water_system([[1.0, 5.0, 5.0], [5.0, 5.0, 9.4]])
    cell = ppo.PeriodicCell([10.0] * 3, periodic=False)
    low = np.full((6, 3), -np.inf)
    high = np.full((6, 3), np.inf)
    low[:3, 0] = 2.0      # first molecule: wall at x = 2
    high[3:, 2] = 9.0     # second molecule: wall at z = 9
    penalty = ppo.OverlapPenalty(system, cell, 0.8, wall_bounds=(low, high))
    energy, gradient = penalty.energy_and_gradient(system.coordinates)
    assert energy > 0.0
    assert np.all(gradient[:3, 0] < 0.0) and np.all(gradient[:3, 1:] == 0.0)
    fn = lambda flat: (penalty.energy_and_gradient(flat.reshape(-1, 3))[0], None)
    assert np.allclose(gradient.ravel(), finite_difference(fn, system.coordinates.ravel()), rtol=1e-5, atol=1e-5)
    with pytest.raises(ValueError, match="every atom"):
        ppo.OverlapPenalty(system, cell, 0.8, wall_bounds=(low[:2], high[:2]))


def test_module_does_not_import_toga():
    code = "import sys, periodicPreOptimizer; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=str(SRC))
    assert out.stdout.strip() == "False", out.stderr
