"""Dipole and second-moment integrals (espIntegrals.moment_matrices) and the ESP of an explicit density."""
from pathlib import Path

import numpy as np
import pytest

import espIntegrals as ei
import gaussianWavefunction as gw
import hirshfeldCharges as hc

FCHK = Path(__file__).resolve().parent / "data" / "fchk"
NAMES = ["water_631gd_5d.fchk", "water_631gd_6d.fchk", "water_ccpvtz.fchk", "water_ccpvqz.fchk",
         "h2co_6311g2df_10f.fchk", "ch3_ub3lyp_631gd.fchk"]


@pytest.fixture(scope="module", params=NAMES)
def case(request):
    w = gw.read_fchk(FCHK / request.param)
    return w, ei.moment_matrices(w)


def test_keys_and_shapes(case):
    w, m = case
    n = w.transform.shape[1]
    assert set(m) == {"s", "x", "y", "z", "xx", "yy", "zz"}
    assert all(v.shape == (n, n) and np.allclose(v, v.T, atol=1e-13) for v in m.values())


def test_overlap_key_equals_overlap_matrix(case):
    w, m = case
    assert np.allclose(m["s"], ei.overlap_matrix(w), atol=1e-13)


def test_moments_match_becke_quadrature(case):
    w, m = case
    g = hc.molecular_grid(w.coordinates, w.atom_symbols())
    ao = gw.evaluate_cartesian_aos(w, g.points)
    p = g.points
    for key, f in {"x": p[:, 0], "y": p[:, 1], "z": p[:, 2], "xx": p[:, 0] ** 2, "yy": p[:, 1] ** 2,
                   "zz": p[:, 2] ** 2}.items():
        ref = (ao * (g.weights * f)[:, None]).T @ ao
        assert np.allclose(m[key], ref, atol=5e-6 * np.abs(ref).max()), key


def test_density_dipole_matches_the_fchk(case):
    w, m = case
    electronic = np.array([np.sum(w.cart_density * m[k]) for k in "xyz"])
    dipole = -electronic + w.nuclear_charges @ w.coordinates
    assert np.allclose(dipole, w.dipole, atol=2e-4)


def test_esp_of_an_explicit_density_equals_the_default():
    w = gw.read_fchk(FCHK / "water_631gd_5d.fchk")
    pts = w.coordinates[0] + np.array([[1.5, 0.3, -0.2], [0.0, 2.0, 1.0], [-1.0, -1.0, 2.5]])
    a = ei.electronic_potential(ei.build_hermite_batches(w), pts)
    b = ei.electronic_potential(ei.build_hermite_batches(w, density=w.density), pts)
    assert np.array_equal(a, b)
    half = ei.electronic_potential(ei.build_hermite_batches(w, density=0.5 * w.density), pts)
    assert np.allclose(half, 0.5 * a, rtol=1e-12)
