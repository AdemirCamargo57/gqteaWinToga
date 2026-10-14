"""Tests for the exact ESP integrals."""
import numpy as np
import pytest
from scipy.integrate import quad

import espIntegrals as ei
import gaussianWavefunction as gw
from esp_fixtures import REFERENCE_FIXTURES, fchk_path, load_reference


class TestBoys:
    @pytest.mark.parametrize("t", [0.0, 1e-12, 0.3, 5.0, 30.0, 80.0])
    def test_against_quadrature(self, t):
        f = ei.boys(8, np.array([t]))
        for n in range(9):
            ref = quad(lambda u: u ** (2 * n) * np.exp(-t * u * u), 0.0, 1.0, epsabs=1e-15)[0]
            assert f[n, 0] == pytest.approx(ref, rel=1e-10, abs=1e-15)


class TestOverlap:
    @pytest.mark.parametrize("name", REFERENCE_FIXTURES)
    def test_trace_ps_is_the_electron_count(self, name):
        w = gw.read_fchk(fchk_path(name))
        s = ei.overlap_matrix(w)
        assert np.trace(w.cart_density @ s) == pytest.approx(w.n_electrons, abs=1e-6)
        # every basis function as written in the fchk is normalised
        assert np.allclose(np.diag(w.transform @ s @ w.transform.T), 1.0, atol=1e-7)


class TestPotential:
    @pytest.mark.parametrize("name", REFERENCE_FIXTURES)
    def test_matches_cubegen(self, name):
        """Gaussian 09 cubegen at 300 points 0.5-8 bohr from the nuclei."""
        w = gw.read_fchk(fchk_path(name))
        points, _, esp_ref = load_reference(name)
        esp = ei.electrostatic_potential(w, points)
        assert np.max(np.abs(esp - esp_ref)) < 1e-6

    def test_threads_do_not_change_values(self):
        w = gw.read_fchk(fchk_path("water_ccpvtz"))
        points, _, _ = load_reference("water_ccpvtz")
        serial = ei.electrostatic_potential(w, points, workers=1, use_numba=False)
        threaded = ei.electrostatic_potential(w, points, workers=6, use_numba=False)
        assert np.allclose(serial, threaded, rtol=0, atol=1e-13)

    def test_far_field_is_the_dipole_potential(self):
        """Neutral water far away: V -> mu . r / r^3 with the fchk dipole."""
        w = gw.read_fchk(fchk_path("water_631gd_6d"))
        r = 150.0
        directions = np.eye(3)
        esp = ei.electrostatic_potential(w, directions * r)
        expected = directions @ w.dipole / r ** 2
        assert np.allclose(esp, expected, atol=2e-6)
        assert abs(esp[2]) > 1e-5  # water's dipole lies along z

    def test_progress_reaches_one(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))
        points, _, _ = load_reference("eteno_b3lyp_631g")
        seen = []
        ei.electrostatic_potential(w, points, workers=1, progress=seen.append)
        assert seen and seen[-1] == pytest.approx(1.0)

    def test_raising_progress_callback_is_ignored(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))
        points, _, esp_ref = load_reference("eteno_b3lyp_631g")

        def boom(_):
            raise RuntimeError("UI went away")

        esp = ei.electrostatic_potential(w, points, workers=2, progress=boom)
        assert np.max(np.abs(esp - esp_ref)) < 1e-6


@pytest.mark.skipif(not ei.numba_available(), reason="numba not installed (optional)")
class TestNumba:
    def test_boys_matches_the_scipy_version(self):
        import espNumba
        ts = np.r_[0.0, 1e-12, np.linspace(0.01, 150.0, 600)]
        ref = ei.boys(8, ts)
        f = np.zeros(9)
        for i, t in enumerate(ts):
            espNumba._boys_into(8, t, f)
            assert np.allclose(f, ref[:, i], rtol=1e-12, atol=0)

    @pytest.mark.parametrize("name", REFERENCE_FIXTURES)
    def test_kernel_matches_numpy(self, name):
        w = gw.read_fchk(fchk_path(name))
        points, _, esp_ref = load_reference(name)
        fast = ei.electrostatic_potential(w, points, use_numba=True)
        slow = ei.electrostatic_potential(w, points, use_numba=False)
        assert np.max(np.abs(fast - slow)) < 1e-10
        assert np.max(np.abs(fast - esp_ref)) < 1e-6


def test_forcing_numba_without_it_is_an_error(monkeypatch):
    monkeypatch.setattr(ei, "espNumba", None)
    w = gw.read_fchk(fchk_path("water_631gd_6d"))
    with pytest.raises(ValueError, match="numba"):
        ei.electrostatic_potential(w, np.zeros((1, 3)) + 5.0, use_numba=True)
