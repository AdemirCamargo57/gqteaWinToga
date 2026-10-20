"""Analytic AO and electron-density derivatives (gradient, Hessian, G) for QTAIM."""
from pathlib import Path

import numpy as np
import pytest

import gaussianWavefunction as gw

FCHK = Path(__file__).parent / "data" / "fchk"
# 6D, 5D, pure f (cc-pVTZ), pure g (cc-pVQZ), Cartesian 10F, SP shells, open shell
FIXTURES = ["water_631gd_6d.fchk", "water_631gd_5d.fchk", "water_ccpvtz.fchk", "water_ccpvqz.fchk",
            "h2co_6311g2df_10f.fchk", "eteno_b3lyp_631g.fchk", "ch3_ub3lyp_631gd.fchk"]
POINTS = np.array([[0.3, -0.2, 0.5], [1.1, 0.4, -0.7], [-0.6, 1.3, 0.2], [0.0, 0.0, 1.5]])


def random_points(wfn, n, seed):
    rng = np.random.default_rng(seed)
    lo = wfn.coordinates.min(axis=0) - 3.0
    hi = wfn.coordinates.max(axis=0) + 3.0
    return lo + (hi - lo) * rng.random((n, 3))


@pytest.fixture(scope="module", params=FIXTURES)
def wfn(request):
    return gw.read_fchk(FCHK / request.param)


class TestAODerivatives:
    def test_derivative_order_table(self):
        assert gw.DERIVATIVE_ORDER[:4] == ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        assert gw.SECOND_PAIRS == ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))

    def test_shapes(self, wfn):
        n_cart = sum(s.n_cart for s in wfn.shells)
        assert gw.evaluate_cartesian_ao_derivatives(wfn, POINTS, order=1, use_numba=False).shape == (4, 4, n_cart)
        assert gw.evaluate_cartesian_ao_derivatives(wfn, POINTS, order=2, use_numba=False).shape == (4, 10, n_cart)

    def test_value_matches_evaluate_cartesian_aos(self, wfn):
        d = gw.evaluate_cartesian_ao_derivatives(wfn, POINTS, order=1, use_numba=False)
        assert np.allclose(d[:, 0], gw.evaluate_cartesian_aos(wfn, POINTS, use_numba=False), rtol=1e-13, atol=0)

    def test_first_and_second_derivatives_match_finite_differences(self, wfn):
        d = gw.evaluate_cartesian_ao_derivatives(wfn, POINTS, order=2, use_numba=False)
        h = 1e-5
        for axis in range(3):
            e = np.zeros(3)
            e[axis] = h
            plus = gw.evaluate_cartesian_ao_derivatives(wfn, POINTS + e, order=1, use_numba=False)
            minus = gw.evaluate_cartesian_ao_derivatives(wfn, POINTS - e, order=1, use_numba=False)
            scale = np.abs(d[:, 1 + axis]).max()
            assert np.allclose((plus[:, 0] - minus[:, 0]) / (2 * h), d[:, 1 + axis], atol=1e-7 * scale, rtol=1e-6)
            for k, (a, b) in enumerate(gw.SECOND_PAIRS):
                if a != axis:
                    continue
                fd = (plus[:, 1 + b] - minus[:, 1 + b]) / (2 * h)
                scale = np.abs(d[:, 4 + k]).max()
                assert np.allclose(fd, d[:, 4 + k], atol=1e-6 * scale, rtol=1e-5)

    def test_order_must_be_one_or_two(self, wfn):
        with pytest.raises(ValueError, match="order"):
            gw.evaluate_cartesian_ao_derivatives(wfn, POINTS, order=3, use_numba=False)


class TestDensityDerivatives:
    def test_rho_matches_electron_density(self, wfn):
        dd = gw.density_derivatives(wfn, POINTS, use_numba=False)
        assert np.allclose(dd.rho, gw.electron_density(wfn, POINTS, use_numba=False), rtol=1e-12, atol=0)

    def test_gradient_and_hessian_match_finite_differences(self, wfn):
        dd = gw.density_derivatives(wfn, POINTS, use_numba=False)
        h = 1e-5
        for axis in range(3):
            e = np.zeros(3)
            e[axis] = h
            p = gw.density_derivatives(wfn, POINTS + e, order=1, use_numba=False)
            m = gw.density_derivatives(wfn, POINTS - e, order=1, use_numba=False)
            assert np.allclose((p.rho - m.rho) / (2 * h), dd.gradient[:, axis], rtol=1e-6, atol=1e-9)
            assert np.allclose((p.gradient - m.gradient) / (2 * h), dd.hessian[:, axis, :], rtol=1e-5, atol=1e-8)

    def test_hessian_is_symmetric_and_laplacian_is_its_trace(self, wfn):
        dd = gw.density_derivatives(wfn, POINTS, use_numba=False)
        assert np.array_equal(dd.hessian, np.transpose(dd.hessian, (0, 2, 1)))
        assert np.allclose(dd.laplacian, np.trace(dd.hessian, axis1=1, axis2=2))

    def test_order_one_has_no_hessian(self, wfn):
        dd = gw.density_derivatives(wfn, POINTS, order=1, use_numba=False)
        assert dd.hessian is None
        with pytest.raises(ValueError, match="order=2"):
            _ = dd.laplacian

    def test_kinetic_g_is_half_sum_of_gradient_products(self, wfn):
        d = gw.evaluate_cartesian_ao_derivatives(wfn, POINTS, order=1, use_numba=False)
        g = 0.5 * sum(np.einsum("pa,ab,pb->p", d[:, k], wfn.cart_density, d[:, k]) for k in (1, 2, 3))
        dd = gw.density_derivatives(wfn, POINTS, use_numba=False)
        assert np.allclose(dd.kinetic_g, g, rtol=1e-12, atol=0)
        assert np.all(dd.kinetic_g >= 0)

    def test_unknown_order_is_rejected(self, wfn):
        with pytest.raises(ValueError, match="order"):
            gw.density_derivatives(wfn, POINTS, order=0, use_numba=False)

    def test_chunking_does_not_change_results(self, wfn):
        pts = random_points(wfn, 37, seed=1)
        a = gw.density_derivatives(wfn, pts, chunk=5, use_numba=False)
        b = gw.density_derivatives(wfn, pts, chunk=1000, use_numba=False)
        assert np.allclose(a.rho, b.rho, rtol=1e-13, atol=0)
        assert np.allclose(a.hessian, b.hessian, rtol=1e-12, atol=1e-15)


class TestNumbaKernels:
    @pytest.fixture(autouse=True)
    def _need_numba(self):
        pytest.importorskip("numba")

    def test_ao_derivatives_agree_with_numpy(self, wfn):
        pts = random_points(wfn, 300, seed=3)
        for order in (1, 2):
            a = gw.evaluate_cartesian_ao_derivatives(wfn, pts, order=order, use_numba=False)
            b = gw.evaluate_cartesian_ao_derivatives(wfn, pts, order=order, use_numba=True)
            assert np.max(np.abs(a - b)) <= 1e-12 * np.max(np.abs(a))

    def test_fused_density_derivatives_agree_with_numpy(self, wfn):
        pts = np.vstack([random_points(wfn, 300, seed=4), wfn.coordinates + 1e-3])
        a = gw.density_derivatives(wfn, pts, use_numba=False)
        b = gw.density_derivatives(wfn, pts, use_numba=True)
        for x, y in ((a.rho, b.rho), (a.gradient, b.gradient), (a.hessian, b.hessian), (a.kinetic_g, b.kinetic_g)):
            assert np.max(np.abs(x - y)) <= 1e-12 * np.max(np.abs(x))
        c = gw.density_derivatives(wfn, pts, order=1, use_numba=True)
        assert c.hessian is None
        assert np.max(np.abs(c.gradient - a.gradient)) <= 1e-12 * np.max(np.abs(a.gradient))

    def test_numba_kernel_module_list_unchanged(self):
        from test_solventPacking import NUMBA_KERNEL_MODULES
        assert NUMBA_KERNEL_MODULES == {"espNumba.py", "densityNumba.py", "solventPackingNumba.py", "volumeNumba.py"}
