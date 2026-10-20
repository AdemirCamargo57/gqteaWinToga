"""QTAIM topology: critical points and gradient paths (qtaimTopology)."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import gaussianWavefunction as gw
import qtaimTopology as qt
from qtaim_fixtures import fchk, make_h2

SRC = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def h2():
    wfn = make_h2()
    cps, notes = qt.find_critical_points(wfn, (2, 3), qt.CPSearchParameters())
    return wfn, cps, notes


@pytest.fixture(scope="module")
def water():
    wfn = gw.read_fchk(fchk("water_631gd_6d.fchk"))
    cps, notes = qt.find_critical_points(wfn, (2, 3, 4, 5), qt.CPSearchParameters())
    return wfn, cps, notes


class TestSyntheticH2:
    def test_helper_density_holds_two_electrons(self):
        import espIntegrals as ei
        wfn = make_h2()
        assert np.trace(wfn.cart_density @ ei.overlap_matrix(wfn)) == pytest.approx(2.0, rel=1e-7)  # 8-decimal STO-3G coefficients

    def test_two_ncps_and_one_bcp_at_the_midpoint(self, h2):
        _, cps, notes = h2
        assert [c.kind for c in cps] == [qt.NCP, qt.NCP, qt.BCP]
        assert not notes
        bcp = cps[2]
        assert np.allclose(bcp.position, [0.0, 0.0, 0.7], atol=1e-6)
        lam = bcp.properties.eigenvalues
        assert lam[0] < 0 and lam[1] < 0 and lam[2] > 0
        assert abs(bcp.properties.ellipticity) < 1e-6   # cylindrical symmetry

    def test_ncps_are_assigned_to_atoms_in_order(self, h2):
        _, cps, _ = h2
        assert [c.assignment for c in cps[:2]] == [("atom", 0), ("atom", 1)]
        assert [c.index for c in cps] == [1, 2, 3]

    def test_properties_follow_multiwfn_definitions(self, h2):
        p = h2[1][2].properties
        assert p.k == pytest.approx(p.g - p.laplacian / 4)
        assert p.v == pytest.approx(-p.g - p.k)
        assert p.h == pytest.approx(-p.k)
        assert p.h == pytest.approx(p.g + p.v)
        assert p.eta == pytest.approx(abs(p.eigenvalues[0]) / p.eigenvalues[2])
        assert p.ellipticity == pytest.approx(p.eigenvalues[0] / p.eigenvalues[1] - 1)
        assert p.gradient_norm < 1e-6
        # s functions only: grad(phi) = 0 at the midpoint by symmetry, so G = 0 exactly
        assert p.g == pytest.approx(0.0, abs=1e-14)
        assert np.isnan(p.v_over_g)


class TestWater:
    def test_counts_match_multiwfn(self, water):
        _, cps, notes = water
        counts = qt.cp_counts(cps)
        assert counts == {qt.NCP: 3, qt.BCP: 2, qt.RCP: 0, qt.CCP: 0}
        assert qt.poincare_hopf(counts) == (1, True)

    def test_bcp_position_matches_multiwfn(self, water):
        # Multiwfn 2026.7.15 CP 5: (0, 1.144645804772, -0.675954276022) bohr
        _, cps, _ = water
        bcps = [c for c in cps if c.kind == qt.BCP]
        target = np.array([0.0, 1.144645804772, -0.675954276022])
        assert min(np.linalg.norm(c.position - target) for c in bcps) < 1e-6

    def test_bcp_properties_match_multiwfn(self, water):
        _, cps, _ = water
        bcp = min((c for c in cps if c.kind == qt.BCP), key=lambda c: c.position[1])
        p = bcp.properties
        assert p.rho == pytest.approx(0.3567996402, rel=1e-8)
        assert p.g == pytest.approx(0.06335348440, rel=1e-8)
        assert p.k == pytest.approx(0.5563370013, rel=1e-8)
        assert p.v == pytest.approx(-0.6196904857, rel=1e-8)
        assert p.laplacian == pytest.approx(-1.971934067, rel=1e-8)
        assert p.eigenvalues == pytest.approx([-1.787374161, -1.746326200, 1.561766294], rel=1e-8)
        assert p.ellipticity == pytest.approx(0.023505, abs=1e-6)
        assert p.eta == pytest.approx(1.144457, abs=1e-6)
        assert p.h == pytest.approx(-0.5563370013, rel=1e-8)
        assert p.v_over_g == pytest.approx(0.6196904857 / 0.06335348440, rel=1e-8)
        assert p.h_over_rho == pytest.approx(-0.5563370013 / 0.3567996402, rel=1e-8)


class TestPieces:
    def test_classify(self):
        assert qt.classify([-1, -2, -3]) == qt.NCP
        assert qt.classify([-1, -2, 3]) == qt.BCP
        assert qt.classify([-1, 2, 3]) == qt.RCP
        assert qt.classify([1, 2, 3]) == qt.CCP
        assert qt.classify([0.0, 2, 3]) is None

    def test_combinations_respect_the_vdw_cutoff(self):
        coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], dtype=float)
        radii = np.full(4, 0.4)
        assert qt.atom_combinations(coords, radii, 2, 1.5) == [(0, 1), (1, 2), (2, 3)]
        assert qt.atom_combinations(coords, radii, 3, 1.5) == []
        assert len(qt.atom_combinations(coords, radii, 3, 3.0)) == 2   # cutoff 2.4: (0,1,2), (1,2,3)

    def test_combinations_are_cliques_sorted_and_deterministic(self):
        rng = np.random.default_rng(7)
        coords = rng.random((12, 3)) * 4
        radii = np.full(12, 1.0)
        a = qt.atom_combinations(coords, radii, 4, 1.0)
        assert a == sorted(a) and a == qt.atom_combinations(coords, radii, 4, 1.0)
        for combo in a:
            for i, j in zip(*np.triu_indices(4, 1)):
                assert np.linalg.norm(coords[combo[i]] - coords[combo[j]]) <= 2.0
        brute = [c for c in __import__("itertools").combinations(range(12), 4)
                 if all(np.linalg.norm(coords[x] - coords[y]) <= 2.0
                        for x, y in __import__("itertools").combinations(c, 2))]
        assert a == brute

    def test_seed_mode_validation(self):
        with pytest.raises(ValueError, match="Unknown CP search mode"):
            qt.seeds_for_mode(7, np.zeros((2, 3)), np.ones(2), 1.5)

    def test_merge_first_found_wins(self):
        existing = [np.zeros(3)]
        cand = [np.array([0.01, 0, 0]), np.array([1.0, 0, 0]), np.array([1.02, 0, 0])]
        kept = qt.merge_critical_points(cand, existing, 0.03)
        assert len(kept) == 1 and np.allclose(kept[0], [1.0, 0, 0])

    def test_poincare_hopf(self):
        assert qt.poincare_hopf({qt.NCP: 12, qt.BCP: 12, qt.RCP: 1}) == (1, True)
        assert qt.poincare_hopf({qt.NCP: 3, qt.BCP: 1}) == (2, False)

    def test_nna_assignment(self):
        props = qt.CPProperties(1, 0, 0, 0, 0, 0, 0, np.array([-1, -1, -1.0]), np.eye(3), 0, 0)
        cp = qt.CriticalPoint(np.array([0.5, 0.5, 0.5]), qt.NCP, props)
        qt.assign_nuclear_cps([cp], np.array([[0, 0, 0], [1, 1, 1.0]]))
        assert cp.assignment == ("nna",)

    @pytest.mark.parametrize("field, value", [("max_iterations", 0), ("step_scale", -1.0),
                                              ("min_cp_distance", float("nan")), ("trust_radius", 0.0),
                                              ("max_iterations", 2.5)])
    def test_search_parameters_validate(self, field, value):
        p = qt.CPSearchParameters(**{field: value})
        with pytest.raises(ValueError, match=field):
            p.validate()

    @pytest.mark.parametrize("field, value", [("max_points", 1), ("step", 0.0), ("stop_distance", -1.0)])
    def test_path_parameters_validate(self, field, value):
        with pytest.raises(ValueError, match=field):
            qt.PathParameters(**{field: value}).validate()

    def test_no_modes_is_an_error(self):
        with pytest.raises(ValueError, match="at least one"):
            qt.find_critical_points(make_h2(), (), qt.CPSearchParameters())

    def test_progress_that_raises_is_harmless(self):
        def bad(_):
            raise RuntimeError("boom")
        cps, _ = qt.find_critical_points(make_h2(), (2, 3), qt.CPSearchParameters(), progress=bad)
        assert len(cps) == 3

    def test_bohr_constant_matches_coordinate_converter(self):
        from coordinateConverter import BOHR_RADIUS_ANGSTROM
        assert qt.BOHR_RADIUS_ANGSTROM == BOHR_RADIUS_ANGSTROM

    def test_module_does_not_import_toga(self):
        code = "import sys, qtaimTopology; print('toga' in sys.modules)"
        out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
        assert out.stdout.strip() == "False"
