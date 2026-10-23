"""Electron excitation analysis (excitationAnalysis): Multiwfn's hole-electron framework, checked against
Multiwfn 2026.7.15 (tests/data/excitation/*_multiwfn.txt) and against basis-free identities."""
import subprocess
import sys
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import excitationAnalysis as xa
import excitationData as ed
import gaussianWavefunction as gw
import excitation_fixtures as xf

SRC = Path(__file__).resolve().parent.parent
A = xa.BOHR_RADIUS_ANGSTROM

# Multiwfn's uniform-grid values differ from its 1.5x coarser grid by at most these (README); the analytic
# values are compared with a margin of a few times that, plus Multiwfn's printed rounding.
TOL = {"centroid": 2e-4, "length": 2e-3, "dipole": 2e-4, "integral": 2e-4}


@lru_cache(maxsize=None)
def loaded(name):
    fchk, log = xf.fixture_paths(name)
    wfn = gw.read_fchk(fchk)
    data = ed.read_excitation_file(log)
    ed.check_against_wavefunction(data, wfn)
    return wfn, data


@lru_cache(maxsize=None)
def moments_of(name):
    return xa.basis_moments(loaded(name)[0])


@lru_cache(maxsize=None)
def core(name, state_number, threshold=0.0):
    wfn, data = loaded(name)
    state = data.state(state_number)
    factor = 1.0 if wfn.is_unrestricted else 2.0
    coeffs = xa.coefficient_matrices(state, wfn)
    dm = xa.density_matrices(coeffs, wfn, factor, cross_term_threshold=threshold)
    return SimpleNamespace(wfn=wfn, data=data, state=state, factor=factor, coeffs=coeffs, dm=dm,
                           moments=xa.moment_indices(dm, moments_of(name)))


def reference(name, state, section="indices"):
    return xf.parse_hole_electron(xf.sections(xf.multiwfn_text(name, state))[section])


class TestDensityMatrices:
    def test_cis_densities_integrate_to_the_sum_of_squares(self):
        r = core("h2co_cis", 1)
        x2, y2 = r.state.norm_terms()
        assert y2 == 0.0
        assert abs(r.moments.hole_integral - 2 * x2) < 1e-8
        assert abs(r.moments.electron_integral - 2 * x2) < 1e-8

    def test_td_densities_integrate_to_x2_minus_y2(self):
        r = core("h2co_td5050", 6)
        x2, y2 = r.state.norm_terms()
        assert y2 > 0
        assert abs(r.moments.hole_integral - 2 * (x2 - y2)) < 1e-8

    @pytest.mark.parametrize("name,state", [("h2co_cis", 1), ("h2co_td5050", 6), ("ch3", 1)])
    def test_transition_density_integrates_to_zero(self, name, state):
        assert abs(core(name, state).moments.transition_integral) < 1e-10

    def test_cis_matrices_are_symmetric_and_positive(self):
        dm = core("h2co_cis", 1).dm
        s = moments_of("h2co_cis")["s"]
        for p in (dm.hole, dm.electron):
            assert np.allclose(p, p.T, atol=1e-14)
            # P S has the non-negative spectrum of a density
            assert np.linalg.eigvals(p @ s).real.min() > -1e-10

    def test_unrestricted_uses_both_spins_and_no_factor(self):
        r = core("ch3", 1)
        assert set(r.coeffs) == {ed.SPIN_ALPHA, ed.SPIN_BETA} and r.factor == 1.0
        x2, y2 = r.state.norm_terms()
        assert abs(r.moments.hole_integral - (x2 - y2)) < 1e-8

    def test_out_of_range_pair_is_reported(self):
        wfn, _ = loaded("h2co_cis")
        bad = ed.parse_excitation_text(f"Excited State 1 1 5.0\n 2 -> {wfn.n_mo + 3} 0.7\n", "x.txt")
        with pytest.raises(ValueError, match="Line 2"):
            xa.coefficient_matrices(bad.state(1), wfn)


class TestMomentsAgainstMultiwfn:
    CASES = [("D-pi-A", 1), ("D-pi-A", 2), ("D-pi-A", 3), ("N-phenylpyrrole", 4), ("4-nitroaniline", 2),
             ("h2co_cis", 1), ("ch3", 1), ("h2co_ryd", 3)]

    @pytest.mark.parametrize("name,state", CASES)
    def test_centroids_rmsd_and_indices(self, name, state):
        m, ref = core(name, state, 0.01).moments, reference(name, state)
        assert np.allclose(m.hole_centroid * A, ref["hole_centroid"], atol=TOL["centroid"])
        assert np.allclose(m.electron_centroid * A, ref["electron_centroid"], atol=TOL["centroid"])
        assert abs(m.d_index * A - ref["d"]) < TOL["length"]
        assert np.allclose(np.append(m.sigma_hole, m.sigma_hole_norm) * A, ref["sigma_hole"], atol=TOL["length"])
        assert np.allclose(np.append(m.sigma_electron, m.sigma_electron_norm) * A, ref["sigma_electron"],
                           atol=TOL["length"])
        assert np.allclose(np.append(m.delta_sigma, m.delta_sigma_index) * A, ref["delta_sigma"], atol=TOL["length"])
        assert np.allclose(np.concatenate([m.h, [m.h_ct, m.h_index]]) * A, ref["h"], atol=TOL["length"])
        assert abs(m.t_index * A - ref["t"]) < TOL["length"]

    @pytest.mark.parametrize("name,state", CASES)
    def test_transition_dipole_and_dipole_change(self, name, state):
        m, ref = core(name, state, 0.01).moments, reference(name, state)
        assert np.allclose(m.transition_dipole, ref["transition_dipole"], atol=TOL["dipole"] * 5)
        assert np.allclose(m.dipole_change, ref["dipole_change"], atol=TOL["dipole"] * 5)

    @pytest.mark.parametrize("name,state", CASES)
    def test_integrals(self, name, state):
        m, ref = core(name, state, 0.01).moments, reference(name, state)
        assert abs(m.hole_integral - ref["hole_integral"]) < TOL["integral"]
        assert abs(m.electron_integral - ref["electron_integral"]) < TOL["integral"]

    def test_exact_and_screened_cross_terms_differ_little(self):
        exact, screened = core("D-pi-A", 2).moments, core("D-pi-A", 2, 0.01).moments
        assert abs(exact.d_index - screened.d_index) * A < 0.02


class TestTransitionDipoleAgainstLog:
    @pytest.mark.parametrize("name", ["D-pi-A", "N-phenylpyrrole", "4-nitroaniline", "uracil", "h2co_cis"])
    def test_magnitude_and_oscillator_strength(self, name):
        wfn, data = loaded(name)
        for s in data.states[:3]:
            if s.multiplicity == 3:
                continue
            m = core(name, s.number).moments
            log = np.array(data.transition_dipoles[s.number])
            # Gaussian prints coefficients with 5 decimals: compare the magnitude to 1 % (+ printed rounding)
            assert abs(np.linalg.norm(m.transition_dipole) - np.linalg.norm(log)) < 0.01 * np.linalg.norm(log) + 2e-3
            assert abs(m.oscillator(s.energy_ev) - s.oscillator) < 0.01 * s.oscillator + 2e-4


class TestOrbitalIndices:
    def test_mo_contributions_match_multiwfn(self):
        r = core("D-pi-A", 1)
        ours = {(c.mo, ""): (100 * c.hole, 100 * c.electron) for c in xa.mo_contributions(r.coeffs, r.factor)}
        ref = xf.parse_mo_contributions(xf.sections(xf.multiwfn_text("D-pi-A", 1))["mo"])
        for key, (hole, ele) in ref.items():
            got = ours.get(key, (0.0, 0.0))
            assert abs(got[0] - hole) < 6e-4 and abs(got[1] - ele) < 6e-4, key

    def test_unrestricted_mo_contributions_carry_the_spin(self):
        r = core("ch3", 1)
        contributions = xa.mo_contributions(r.coeffs, r.factor)
        ref = xf.parse_mo_contributions(xf.sections(xf.multiwfn_text("ch3", 1))["mo"])
        ours = {(c.mo, "A" if c.spin == ed.SPIN_ALPHA else "B"): (100 * c.hole, 100 * c.electron)
                for c in contributions}
        for key, (hole, ele) in ref.items():
            got = ours.get(key, (0.0, 0.0))
            assert abs(got[0] - hole) < 6e-4 and abs(got[1] - ele) < 6e-4, key

    @pytest.mark.parametrize("name", ["D-pi-A", "N-phenylpyrrole", "h2co_cis", "h2co_td5050", "ch3"])
    def test_major_transitions_match_multiwfn(self, name):
        wfn, data = loaded(name)
        ref = xf.parse_major_transitions(xf.sections(xf.multiwfn_text(name))["major"])
        for s in data.states:
            assert xa.major_transitions(s, wfn) == ref[s.number], s.number

    def test_ghost_hunter_matches_multiwfn(self):
        for name, state in (("D-pi-A", 1), ("4-nitroaniline", 2)):
            r = core(name, state, 0.01)
            ref = reference(name, state)["ghost"]
            g = xa.ghost_hunter(r.state, r.wfn, ref_d := r.moments.d_index)
            got = np.array([g.value, g.term1, g.term2]) * xa.HARTREE_TO_EV
            assert abs(got[1] - ref[1]) < 1.5e-3
            assert np.allclose(got, ref, atol=0.02), (got, ref, ref_d)

    @pytest.mark.parametrize("name", ["N-phenylpyrrole", "D-pi-A", "h2co_td5050", "ch3"])
    def test_delta_r_matches_multiwfn(self, name):
        wfn, data = loaded(name)
        ref = xf.parse_delta_r(xf.sections(xf.multiwfn_text(name))["delta_r"])
        for s in data.states:
            assert abs(xa.delta_r(s, wfn, moments_of(name)).value * A - ref[s.number]) < 2e-6, s.number

    def test_delta_r_manual_values(self):
        wfn, data = loaded("N-phenylpyrrole")
        values = [xa.delta_r(s, wfn, moments_of("N-phenylpyrrole")).value * A for s in data.states]
        assert np.allclose(values, [0.793368, 1.846333, 2.455982, 3.105966, 3.752463], atol=2e-6)

    def test_one_pair_gives_d_equal_to_delta_r(self):
        wfn, _ = loaded("h2co_cis")
        text = f"Excited State 1 1 5.0\n {wfn.n_alpha} -> {wfn.n_alpha + 2} 0.70710678\n"
        state = ed.parse_excitation_text(text, "one.txt").state(1)
        coeffs = xa.coefficient_matrices(state, wfn)
        m = xa.moment_indices(xa.density_matrices(coeffs, wfn, 2.0), moments_of("h2co_cis"))
        assert abs(m.d_index - xa.delta_r(state, wfn, moments_of("h2co_cis")).value) < 1e-8


@lru_cache(maxsize=None)
def grid_of(name):
    return xa.molecular_grid(loaded(name)[0])


class TestGridIndices:
    CASES = [("D-pi-A", 1), ("D-pi-A", 2), ("N-phenylpyrrole", 4), ("4-nitroaniline", 2), ("h2co_cis", 1),
             ("ch3", 1), ("h2co_ryd", 3)]

    @pytest.mark.parametrize("name,state", CASES)
    def test_sm_sr_hdi_edi_match_multiwfn(self, name, state):
        r = core(name, state, 0.01)
        g = xa.grid_indices(r.wfn, r.dm, grid_of(name))
        ref, coarse = reference(name, state), reference(name, state, "indices_coarse")
        for key, ours, scale in (("sm", g.sm, 1), ("sr", g.sr, 1), ("hdi", g.hdi, 100), ("edi", g.edi, 100)):
            grid_error = abs(ref[key] - coarse[key])
            assert abs(ours - ref[key]) < 3 * grid_error + 1e-3 * scale, (key, ours, ref[key], coarse[key])

    @pytest.mark.parametrize("name,state", [("h2co_cis", 1), ("D-pi-A", 2), ("h2co_ryd", 3)])
    def test_becke_integrals_equal_the_analytic_ones(self, name, state):
        r = core(name, state)
        g = xa.grid_indices(r.wfn, r.dm, grid_of(name))
        assert abs(g.hole_integral - r.moments.hole_integral) < 1e-4
        assert abs(g.electron_integral - r.moments.electron_integral) < 1e-4

    @pytest.mark.parametrize("name", ["N-phenylpyrrole", "D-pi-A", "h2co_cis", "h2co_ryd", "ch3"])
    def test_lambda_matches_multiwfn(self, name):
        wfn, data = loaded(name)
        refs = xf.sections(xf.multiwfn_text(name))
        same_grid, default = xf.parse_lambda(refs["lambda_75x434"]), xf.parse_lambda(refs["lambda"])
        for s in data.states:
            value = xa.lambda_index(s, wfn, grid_of(name)).value
            # Multiwfn on 75 x 434 points per atom, as here (its radial mapping differs slightly)
            assert abs(value - same_grid[s.number]) < 5e-4, s.number
            # Multiwfn's automatic 15 x 170 grid: |phi| has cusps at the nodes, so it is ~1e-2 off
            assert abs(value - default[s.number]) < 1e-2, s.number

    def test_one_pair_gives_lambda_equal_to_sr(self):
        wfn, _ = loaded("h2co_cis")
        text = f"Excited State 1 1 5.0\n {wfn.n_alpha} -> {wfn.n_alpha + 2} 0.70710678\n"
        state = ed.parse_excitation_text(text, "one.txt").state(1)
        dm = xa.density_matrices(xa.coefficient_matrices(state, wfn), wfn, 2.0)
        g = xa.grid_indices(wfn, dm, grid_of("h2co_cis"))
        assert abs(xa.lambda_index(state, wfn, grid_of("h2co_cis")).value - g.sr) < 1e-6


def test_numerical_modules_do_not_import_toga():
    code = "import sys, excitationAnalysis; print('toga' in sys.modules, 'molecularViewer' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False False"
