"""MO coefficients and orbital values read from Gaussian formatted checkpoints."""
import numpy as np
import pytest

import espIntegrals
import gaussianWavefunction as gw
from esp_fixtures import fchk_path

FIXTURES = ("water_631gd_5d", "water_631gd_6d", "water_ccpvtz", "water_ccpvqz",
            "h2co_6311g2df_10f", "eteno_b3lyp_631g", "benzene_631gd", "ch3_ub3lyp_631gd")


@pytest.fixture(scope="module", params=FIXTURES)
def wfn(request):
    return gw.read_fchk(fchk_path(request.param))


def points_around(w, n=400, seed=7):
    rng = np.random.default_rng(seed)
    centres = w.coordinates[rng.integers(0, w.n_atoms, n)]
    return centres + rng.normal(scale=1.5, size=(n, 3))


def occupied(w):
    sets = [(w.alpha_mo[:, :w.n_alpha], 1.0 if w.is_unrestricted else 2.0)]
    if w.is_unrestricted:
        sets.append((w.beta_mo[:, :w.n_beta], 1.0))
    return sets


def test_sizes(wfn):
    assert wfn.n_basis == wfn.density.shape[0]
    assert wfn.alpha_mo.shape == (wfn.n_basis, wfn.n_mo)
    assert wfn.alpha_orbital_energies.shape == (wfn.n_mo,)
    if wfn.is_unrestricted:
        assert wfn.beta_mo.shape == wfn.alpha_mo.shape and wfn.beta_orbital_energies.shape == (wfn.n_mo,)
    else:
        assert wfn.beta_mo is None and wfn.beta_orbital_energies is None


def test_occupied_mos_rebuild_the_density_matrix(wfn):
    density = sum(n * c @ c.T for c, n in occupied(wfn))
    assert np.allclose(density, wfn.density, atol=2e-7)


def test_mos_are_orthonormal(wfn):
    s = wfn.transform @ espIntegrals.overlap_matrix(wfn) @ wfn.transform.T  # fchk-basis overlap
    for c in (wfn.alpha_mo, wfn.beta_mo):
        if c is not None:
            assert np.allclose(c.T @ s @ c, np.eye(wfn.n_mo), atol=1e-6)


def test_occupied_orbitals_sum_to_the_density(wfn):
    pts = points_around(wfn)
    rho = sum(n * (gw.orbital_values(wfn, c, pts) ** 2).sum(axis=1) for c, n in occupied(wfn))
    assert np.allclose(rho, gw.electron_density(wfn, pts), rtol=1e-6, atol=1e-10)


def test_a_vector_gives_one_column_and_matches_the_matrix_call():
    w = gw.read_fchk(fchk_path("water_631gd_5d"))
    pts = points_around(w, 50)
    many = gw.orbital_values(w, w.alpha_mo[:, 2:5], pts, chunk=17)
    one = gw.orbital_values(w, w.alpha_mo[:, 3], pts)
    assert one.shape == (50, 1) and np.allclose(one[:, 0], many[:, 1], rtol=1e-13, atol=1e-15)


@pytest.mark.skipif(not gw.numba_available(), reason="numba not installed")
def test_numba_and_numpy_orbital_values_agree():
    w = gw.read_fchk(fchk_path("water_ccpvqz"))
    pts = points_around(w, 300)
    a = gw.orbital_values(w, w.alpha_mo[:, :8], pts, use_numba=True)
    b = gw.orbital_values(w, w.alpha_mo[:, :8], pts, use_numba=False)
    assert np.max(np.abs(a - b)) <= 1e-12 * np.max(np.abs(b))


def test_mo_coefficients_by_spin():
    r = gw.read_fchk(fchk_path("water_631gd_5d"))
    u = gw.read_fchk(fchk_path("ch3_ub3lyp_631gd"))
    assert gw.mo_coefficients(r, "alpha") is r.alpha_mo and gw.mo_coefficients(u, "beta") is u.beta_mo
    with pytest.raises(ValueError, match="no beta MO coefficients"):
        gw.mo_coefficients(r, "beta")
    with pytest.raises(ValueError, match="spin must be"):
        gw.mo_coefficients(r, "gamma")
    r.alpha_mo = None
    with pytest.raises(ValueError, match="no alpha MO coefficients"):
        gw.mo_coefficients(r, "alpha")


def test_basis_of_a_td_job():
    """A TD fchk names the TD level between method and basis: 'SP  RPBE1PBE TD-FC  6-31G(d)'."""
    from nto_fixtures import fixture
    w = gw.read_fchk(fixture("uracil")[0])
    assert (w.method, w.basis_set) == ("RPBE1PBE", "6-31G(d)")
    assert gw.read_fchk(fchk_path("water_631gd_5d")).basis_set == "6-31G(d)"


def test_wrong_coefficient_count_is_reported(tmp_path):
    text = fchk_path("water_631gd_5d").read_text()
    lines = text.splitlines()
    k = next(i for i, l in enumerate(lines) if l.startswith("Alpha MO coefficients"))
    n = int(lines[k].split("N=")[1])
    lines[k] = lines[k].replace(f"N={n:>12d}", f"N={n - 5:>12d}")
    end = k + 1  # rewrite the section with n - 5 values, keeping the rest of the file
    while end < len(lines) and lines[end].startswith(" "):
        end += 1
    flat = " ".join(lines[k + 1:end]).split()[: n - 5]
    rows = ["".join(f"{float(x):16.8E}" for x in flat[i:i + 5]) for i in range(0, len(flat), 5)]
    bad = lines[:k + 1] + rows + lines[end:]
    path = tmp_path / "bad"
    path.write_text("\n".join(bad) + "\n")
    with pytest.raises(ValueError, match="Alpha MO coefficients"):
        gw.read_fchk(path)
