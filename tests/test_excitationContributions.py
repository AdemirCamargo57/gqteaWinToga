"""Atom and fragment contributions to hole and electron (Mulliken-like, Hirshfeld) and IFCT, against
Multiwfn 2026.7.15 (tests/data/excitation) and the Multiwfn manual (section 4.18.8.1)."""
from functools import lru_cache

import numpy as np
import pytest

import excitationAnalysis as xa
import excitationData as ed
import gaussianWavefunction as gw
import excitation_fixtures as xf


@lru_cache(maxsize=None)
def setup(name, state_number, threshold=0.0):
    fchk, log = xf.fixture_paths(name)
    wfn = gw.read_fchk(fchk)
    state = ed.read_excitation_file(log).state(state_number)
    factor = 1.0 if wfn.is_unrestricted else 2.0
    coeffs = xa.coefficient_matrices(state, wfn)
    dm = xa.density_matrices(coeffs, wfn, factor, cross_term_threshold=threshold)
    return wfn, coeffs, dm, factor, xa.basis_moments(wfn)["s"]


@lru_cache(maxsize=None)
def grid_of(name):
    return xa.molecular_grid(setup(name, 1)[0])


def mulliken(name, state, screen=0.001):
    wfn, coeffs, dm, factor, s = setup(name, state)
    return xa.mulliken_contributions(dm, s, wfn, pair_screen=screen, coeffs=coeffs, factor=factor)


def reference(name, state, section):
    return xf.sections(xf.multiwfn_text(name, state))[section]


class TestMulliken:
    @pytest.mark.parametrize("name,state", [("D-pi-A", 1), ("D-pi-A", 2), ("4-nitroaniline", 2), ("h2co_cis", 1),
                                            ("ch3", 1), ("h2co_td5050", 6)])
    def test_matches_multiwfn_with_its_pair_screen(self, name, state):
        ours = mulliken(name, state)
        ref = xf.parse_atom_contributions(reference(name, state, "mulliken"))
        assert set(ref) == set(range(1, len(ours.hole) + 1))
        for atom, (hole, ele, overlap, diff) in ref.items():
            k = atom - 1
            assert abs(100 * ours.hole[k] - hole) < 6e-3, atom
            assert abs(100 * ours.electron[k] - ele) < 6e-3, atom
            assert abs(100 * ours.overlap[k] - overlap) < 1.2e-2, atom
            assert abs(100 * ours.difference[k] - diff) < 1.2e-2, atom

    def test_exact_partition_sums_to_one_and_is_close_to_the_screened_one(self):
        exact, screened = mulliken("D-pi-A", 2, 0.0), mulliken("D-pi-A", 2)
        assert abs(exact.hole.sum() - 1) < 1e-12 and abs(exact.electron.sum() - 1) < 1e-12
        assert np.abs(exact.hole - screened.hole).max() < 2e-3

    def test_overlap_is_zero_for_negative_contributions(self):
        c = xa.AtomContributions("mulliken", np.array([0.6, 0.5, -0.1]), np.array([0.5, -0.2, 0.7]), True)
        assert np.allclose(c.overlap, [np.sqrt(0.3), 0.0, 0.0])
        assert np.allclose(c.difference, [-0.1, -0.7, 0.8])


class TestHirshfeld:
    @staticmethod
    def multiwfn_scheme(wfn, dm, radcut=10.0):
        """Multiwfn's atmcontri_holeele_Hirshfeld: w_A rho on atom A's own (single-centre) grid, no Becke
        weights, points beyond radcut = 10 bohr from A dropped, then renormalised."""
        raw = np.zeros((2, wfn.n_atoms))
        for a, symbol in enumerate(wfn.atom_symbols()):
            g = xa.hc.molecular_grid(wfn.coordinates[a:a + 1], [symbol])
            keep = np.linalg.norm(g.points - wfn.coordinates[a], axis=1) <= radcut
            pts, w = g.points[keep], g.weights[keep]
            share = xa.hc.hirshfeld_shares(pts, wfn.coordinates, wfn.atom_symbols())[:, a]
            hole, electron, _ = xa.point_densities(wfn, dm, pts)
            raw[:, a] = [(w * share) @ hole, (w * share) @ electron]
        return raw[0] / raw[0].sum(), raw[1] / raw[1].sum()

    # Agreement to the printed 2 decimals (%), except for the diffuse electron of N-phenylpyrrole S4
    # (6-31+G(d)): there the stockholder share far from the molecule is a ratio of free-atom density
    # tails, which differ between Multiwfn's built-in densities and proatomDensities. Investigated and
    # documented in tests/data/excitation/README.md; not an algorithmic difference.
    @pytest.mark.parametrize("name,state,tolerance", [("D-pi-A", 1, 0.01), ("4-nitroaniline", 2, 0.01),
                                                      ("h2co_cis", 1, 0.01), ("N-phenylpyrrole", 4, 0.06)])
    def test_matches_multiwfn(self, name, state, tolerance):
        wfn, coeffs, dm, factor, s = setup(name, state, 0.01)
        ref = xf.parse_atom_contributions(reference(name, state, "hirshfeld"))
        hole, electron = self.multiwfn_scheme(wfn, dm)
        for atom, (h, e, overlap, diff) in ref.items():
            assert abs(100 * hole[atom - 1] - h) < tolerance, atom
            assert abs(100 * electron[atom - 1] - e) < tolerance, atom

    @pytest.mark.parametrize("name,state,tolerance", [("D-pi-A", 1, 0.01), ("4-nitroaniline", 2, 0.01),
                                                      ("h2co_cis", 1, 0.01), ("N-phenylpyrrole", 4, 0.06)])
    def test_becke_grid_agrees_with_the_multiwfn_scheme(self, name, state, tolerance):
        # The Becke molecular grid keeps the far tail of a diffuse electron (N-phenylpyrrole S4) that the
        # 10-bohr single-centre grids drop; otherwise the two quadratures agree to the printed digits.
        wfn, coeffs, dm, factor, s = setup(name, state, 0.01)
        ours = xa.hirshfeld_contributions(wfn, dm, grid_of(name))
        hole, electron = self.multiwfn_scheme(wfn, dm)
        assert np.abs(100 * (ours.hole - hole)).max() < tolerance
        assert np.abs(100 * (ours.electron - electron)).max() < tolerance

    def test_ecp_is_refused(self):
        wfn = setup("h2co_cis", 1)[0]
        ecp = gw.read_fchk(wfn.path)
        ecp.nuclear_charges = ecp.nuclear_charges - np.array([2.0, 2.0, 0.0, 0.0])
        with pytest.raises(ValueError, match="effective core potentials"):
            xa.hirshfeld_contributions(ecp, setup("h2co_cis", 1)[2], grid_of("h2co_cis"))

    def test_element_beyond_argon_is_refused(self):
        heavy = gw.read_fchk(setup("h2co_cis", 1)[0].path)
        heavy.atomic_numbers = np.array([35, 8, 1, 1])
        with pytest.raises(ValueError, match="Br"):
            xa.hirshfeld_contributions(heavy, setup("h2co_cis", 1)[2], grid_of("h2co_cis"))


class TestFragments:
    def test_parse_with_rest(self):
        fragments, notes = xa.parse_fragments("1-3; 5", 6)
        assert [f.atoms for f in fragments] == [(0, 1, 2), (4,), (3, 5)]
        assert [f.name for f in fragments] == ["Fragment 1", "Fragment 2", "Rest"] and fragments[2].automatic
        assert len(notes) == 1

    def test_blank_means_no_fragments(self):
        assert xa.parse_fragments("  ", 6) == ([], [])

    def test_overlap_is_refused(self):
        with pytest.raises(ValueError, match="Atom 3 is in fragments 1 and 2"):
            xa.parse_fragments("1-3; 3-4", 6)

    def test_out_of_range_is_refused(self):
        with pytest.raises(ValueError, match="atom 9 does not exist"):
            xa.parse_fragments("1-3; 9", 6)

    def test_bad_syntax_is_refused(self):
        with pytest.raises(ValueError):
            xa.parse_fragments("1-x", 6)


class TestIFCT:
    @pytest.mark.parametrize("name,state", [("4-nitroaniline", 2), ("D-pi-A", 1), ("D-pi-A", 2)])
    def test_mulliken_ifct_matches_multiwfn(self, name, state):
        wfn = setup(name, state)[0]
        fragments, _ = xa.parse_fragments(xf.FRAGMENTS[name], wfn.n_atoms)
        result = xa.ifct(mulliken(name, state), fragments)
        ref = xf.parse_ifct(reference(name, state, "ifct_mulliken"))
        for k, (hole, ele) in ref["fragments"].items():
            assert abs(100 * result.hole[k - 1] - hole) < 6e-3 and abs(100 * result.electron[k - 1] - ele) < 6e-3
        for k, change in ref["population_change"].items():
            assert abs(result.population_change[k - 1] - change) < 6e-5
        for (i, j), q in ref["transfer"].items():
            assert abs(result.matrix[i - 1, j - 1] - q) < 6e-5
        assert abs(result.intrinsic_ct - ref["intrinsic_ct"]) < 6e-3
        assert abs(result.intrinsic_le - ref["intrinsic_le"]) < 6e-3

    def test_manual_values_for_4_nitroaniline(self):
        wfn = setup("4-nitroaniline", 2)[0]
        result = xa.ifct(mulliken("4-nitroaniline", 2), xa.parse_fragments(xf.FRAGMENTS["4-nitroaniline"],
                                                                          wfn.n_atoms)[0])
        assert np.allclose(100 * result.hole, [28.17, 65.89, 5.94], atol=0.01)
        assert np.allclose(100 * result.electron, [3.02, 29.52, 67.46], atol=0.01)
        assert abs(result.intrinsic_ct - 75.688) < 0.01

    def test_sums_and_two_fragment_apparent_ct(self):
        atoms = xa.AtomContributions("mulliken", np.array([0.7, 0.3]), np.array([0.2, 0.8]), False)
        result = xa.ifct(atoms, [xa.Fragment("A", (0,)), xa.Fragment("B", (1,))])
        assert abs(result.matrix.sum() - 1) < 1e-12
        assert abs(result.intrinsic_ct + result.intrinsic_le - 100) < 1e-10
        assert abs(result.apparent_ct - 50.0) < 1e-10           # |0.2 - 0.7|
        assert result.net_transfers()[0][:2] == (0, 1)

    def test_negative_fragment_is_clamped(self):
        atoms = xa.AtomContributions("mulliken", np.array([1.05, -0.05]), np.array([0.5, 0.5]), True)
        result = xa.ifct(atoms, [xa.Fragment("A", (0,)), xa.Fragment("B", (1,))])
        assert result.clamped and result.hole[1] == 0.0

    def test_atom_ct_matrix(self):
        atoms = xa.AtomContributions("mulliken", np.array([0.6, 0.4]), np.array([-0.1, 1.1]), True)
        assert np.allclose(xa.atom_ct_matrix(atoms), [[0.0, 0.66], [0.0, 0.44]])
