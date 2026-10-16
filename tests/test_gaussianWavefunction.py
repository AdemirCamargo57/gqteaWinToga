"""Tests for the fchk reader, basis expansion and electron density."""
import numpy as np
import pytest

import gaussianWavefunction as gw
from esp_fixtures import REFERENCE_FIXTURES, fchk_path, load_reference

MINI_FCHK = """Title line
SP        RHF                                                         STO-3G
Number of atoms                            I                2
Charge                                     I                0
Total Energy                               R     -1.117506712920000E+00
Atomic numbers                             I   N=           2
           1           1
Nuclear charges                            R   N=           2
  1.00000000E+00  1.00000000E+00
"""


class TestSectionReader:
    def test_scalars_and_arrays(self):
        s = gw.read_fchk_sections(MINI_FCHK)
        assert s["Number of atoms"] == 2
        assert s["Total Energy"] == pytest.approx(-1.11750671292)
        assert s["Atomic numbers"].tolist() == [1, 1]
        assert s["Nuclear charges"].tolist() == [1.0, 1.0]

    def test_truncated_array_is_rejected(self):
        text = MINI_FCHK.replace("           1           1\n", "           1\n")
        with pytest.raises(ValueError, match="Atomic numbers.*truncated"):
            gw.read_fchk_sections(text)

    def test_too_short_file(self):
        with pytest.raises(ValueError, match="too short"):
            gw.read_fchk_sections("x\n")

    def test_missing_section_is_named(self, tmp_path):
        bad = tmp_path / "bad.fchk"
        bad.write_text(MINI_FCHK)
        with pytest.raises(ValueError, match="Current cartesian coordinates"):
            gw.read_fchk(bad)

    def test_unreadable_file(self, tmp_path):
        with pytest.raises(ValueError, match="Cannot read"):
            gw.read_fchk(tmp_path / "absent.fchk")


class TestEtenoFile:
    def test_molecule_and_basis(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))
        assert w.atom_symbols() == ["C", "H", "H", "C", "H", "H"]
        assert w.n_electrons == 16 and w.charge == 0 and w.multiplicity == 1
        assert w.density.shape == (26, 26)
        # 6-31G: C = 1s + 2 SP shells, H = 2 s shells; each SP becomes s + p.
        assert len(w.shells) == 18
        assert np.allclose(w.density, w.density.T)
        assert np.abs(w.dipole).max() < 1e-10  # centrosymmetric: zero dipole

    def test_coordinates_are_bohr(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))
        cc = np.linalg.norm(w.coordinates[0] - w.coordinates[3])
        assert cc * gw_bohr() == pytest.approx(1.336, abs=0.01)  # C=C in Angstrom


def gw_bohr():
    return 0.529177210544


class TestPureToCartesian:
    @pytest.mark.parametrize("l", [1, 2, 3, 4])
    def test_rows_are_orthonormal(self, l):
        t = gw.pure_to_cartesian(l)
        comps = gw.cartesian_components(l)
        s = np.array([[gw._same_centre_overlap(a, b) for b in comps] for a in comps])
        assert t.shape == (2 * l + 1, len(comps))
        assert np.allclose(t @ s @ t.T, np.eye(2 * l + 1), atol=1e-12)

    def test_d_shell_gaussian_order(self):
        # Rows D0, D+1, D-1, D+2, D-2 over XX, YY, ZZ, XY, XZ, YZ.
        t = gw.pure_to_cartesian(2)
        nz = [tuple(np.nonzero(np.abs(row) > 1e-12)[0]) for row in t]
        assert nz == [(0, 1, 2), (4,), (5,), (0, 1), (3,)]
        assert t[3, 0] == pytest.approx(-t[3, 1])  # x^2 - y^2
        assert t[1, 4] > 0 and t[2, 5] > 0 and t[4, 3] > 0


class TestDensity:
    @pytest.mark.parametrize("name", REFERENCE_FIXTURES)
    def test_matches_cubegen_in_the_surface_band(self, name):
        """cubegen's density carries its own screening error (worst seen
        2.5e-4 relative); the ESP test is the stringent check of the basis."""
        w = gw.read_fchk(fchk_path(name))
        points, rho_ref, _ = load_reference(name)
        band = (rho_ref > 1e-5) & (rho_ref < 1e-2)
        assert band.sum() > 30
        rho = gw.electron_density(w, points[band])
        assert np.max(np.abs(rho - rho_ref[band]) / rho_ref[band]) < 1e-3

    def test_chunking_does_not_change_values(self):
        w = gw.read_fchk(fchk_path("water_ccpvtz"))
        points, _, _ = load_reference("water_ccpvtz")
        assert np.allclose(gw.electron_density(w, points, chunk=7), gw.electron_density(w, points))


class TestCdftFields:
    """Fields added for the Fukui tool: energy, method/basis, electron counts, orbital energies."""

    def test_energy_method_basis_and_orbitals(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))
        assert w.method == "RB3LYP" and w.basis_set == "6-31G"
        assert w.total_energy == pytest.approx(-78.572035869, abs=1e-8)
        assert (w.n_alpha, w.n_beta) == (8, 8) and w.other_densities == []
        assert w.alpha_orbital_energies[w.n_alpha - 1] == pytest.approx(-0.26819427, abs=1e-8)

    def test_open_shell_counts(self):
        w = gw.read_fchk(fchk_path("ch3_ub3lyp_631gd"))
        assert w.method == "UB3LYP" and (w.n_alpha, w.n_beta) == (5, 4)


class TestNumbaDensity:
    """densityNumba must reproduce the numpy density path on every basis-set kind."""

    @pytest.mark.skipif(not gw.numba_available(), reason="numba not installed (optional)")
    @pytest.mark.parametrize("name", REFERENCE_FIXTURES + ["benzene_631gd"])
    def test_kernel_matches_numpy(self, name):
        w = gw.read_fchk(fchk_path(name))
        points, _, _ = load_reference(name) if name in REFERENCE_FIXTURES else (
            np.random.default_rng(1).normal(scale=3.0, size=(500, 3)), None, None)
        points = np.vstack([points, w.coordinates, w.coordinates + 0.05])  # include nuclei
        fast = gw.electron_density(w, points, use_numba=True)
        slow = gw.electron_density(w, points, use_numba=False)
        assert np.max(np.abs(fast - slow)) <= 1e-12 * np.max(np.abs(slow))

    @pytest.mark.skipif(not gw.numba_available(), reason="numba not installed (optional)")
    def test_aos_match_numpy_including_far_points(self):
        w = gw.read_fchk(fchk_path("water_ccpvqz"))
        points = np.array([[0.0, 0.0, 0.0], [0.3, -0.2, 0.1], [4.0, 1.0, -2.0], [40.0, 0.0, 0.0]])
        fast = gw.evaluate_cartesian_aos(w, points, use_numba=True)
        slow = gw.evaluate_cartesian_aos(w, points, use_numba=False)
        assert np.allclose(fast, slow, rtol=1e-12, atol=1e-20)

    def test_forcing_numba_without_it_is_an_error(self, monkeypatch):
        monkeypatch.setattr(gw, "densityNumba", None)
        w = gw.read_fchk(fchk_path("water_631gd_6d"))
        with pytest.raises(ValueError, match="numba"):
            gw.electron_density(w, np.zeros((1, 3)), use_numba=True)
        assert gw.electron_density(w, np.zeros((1, 3))).shape == (1,)  # automatic -> numpy
