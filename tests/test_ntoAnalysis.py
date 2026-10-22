"""Natural transition orbitals: SVD of the transition matrix, validated against Multiwfn 2026.7.15."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import espIntegrals
import excitationData as ed
import gaussianWavefunction as gw
import ntoAnalysis as na
from nto_fixtures import REFERENCE_STATES, fixture, multiwfn_reference

SRC = Path(__file__).resolve().parent.parent
ANALYSER = na.NTOAnalyser()


def run(name, state):
    fchk, log = fixture(name)
    return ANALYSER.run(str(fchk), str(log), state)


@pytest.fixture(scope="module")
def uracil3():
    return run("uracil", 3)


def fchk_overlap(wfn):
    return wfn.transform @ espIntegrals.overlap_matrix(wfn) @ wfn.transform.T


def sample_points(wfn, n=1500, seed=3):
    rng = np.random.default_rng(seed)
    centres = wfn.coordinates[rng.integers(0, wfn.n_atoms, n)]
    return centres + rng.normal(scale=1.2, size=(n, 3))


# --------------------------------------------------------------------------- #
# Transition matrix and SVD                                                    #
# --------------------------------------------------------------------------- #
def make_state(*lines):
    configs = tuple(ed.Configuration(o, v, spin, back, coef, k + 1)
                    for k, (o, v, spin, back, coef) in enumerate(lines))
    return ed.ExcitedState(1, "Singlet-A", 1, 5.0, None, None, None, configs)


class TestTransitionMatrix:
    def test_excitations_only_and_one_spin(self):
        state = make_state((2, 3, "alpha", False, 0.6), (1, 4, "alpha", False, -0.3),
                           (2, 4, "alpha", True, 0.05), (1, 3, "beta", False, 0.2))
        t = na.transition_matrix(state, "alpha", n_occ=2, n_mo=4)
        assert t.tolist() == [[0.0, -0.3], [0.6, 0.0]]
        assert na.transition_matrix(state, "beta", n_occ=2, n_mo=4).tolist() == [[0.2, 0.0], [0.0, 0.0]]

    def test_out_of_range_pair_names_the_line(self):
        state = make_state((2, 3, "alpha", False, 0.6), (3, 4, "alpha", False, 0.1))
        with pytest.raises(ValueError, match="Line 2: 3 -> 4"):
            na.transition_matrix(state, "alpha", n_occ=2, n_mo=4)


class TestSVD:
    def test_identities(self, uracil3):
        s = uracil3.sets["alpha"]
        wfn = uracil3.wavefunction
        t = na.transition_matrix(uracil3.state, "alpha", s.n_occ, wfn.n_mo)
        sigma = np.sqrt(s.weights / uracil3.factor)
        k = len(sigma)
        assert np.allclose(s.hole_mo[:, :k] @ np.diag(sigma) @ s.particle_mo[:, :k].T, t, atol=1e-12)
        assert np.allclose(s.hole_mo.T @ s.hole_mo, np.eye(s.n_occ), atol=1e-12)
        assert np.allclose(s.particle_mo.T @ s.particle_mo, np.eye(s.particle_mo.shape[0]), atol=1e-12)
        overlap = fchk_overlap(wfn)
        assert np.allclose(s.holes.T @ overlap @ s.holes, np.eye(s.n_occ), atol=1e-6)
        assert np.allclose(s.particles.T @ overlap @ s.particles, np.eye(s.particles.shape[1]), atol=1e-6)
        assert np.all(np.diff(s.weights) <= 0)

    def test_sign_rule(self, uracil3):
        u = uracil3.sets["alpha"].hole_mo
        for k in range(len(uracil3.sets["alpha"].weights)):
            assert u[np.argmax(np.abs(u[:, k])), k] > 0

    def test_composition_sums_to_one(self, uracil3):
        s = uracil3.sets["alpha"]
        hole = s.composition("hole", 1, threshold=0.0)
        assert sum(f for _, f in hole) == pytest.approx(1.0)
        assert all(1 <= mo <= s.n_occ for mo, _ in hole)
        particle = s.composition("particle", 1, threshold=0.05)
        assert all(mo > s.n_occ and f >= 0.05 for mo, f in particle)
        assert [f for _, f in particle] == sorted((f for _, f in particle), reverse=True)


# --------------------------------------------------------------------------- #
# Published values and physics                                                 #
# --------------------------------------------------------------------------- #
class TestResults:
    def test_uracil_s3_matches_the_manual(self, uracil3):
        # Multiwfn manual section 4.18.6
        w = uracil3.sets["alpha"].weights
        assert np.allclose(w[:5], [0.865529, 0.134025, 0.000582, 0.000121, 0.000063], atol=5e-7)
        assert w[:10].sum() == pytest.approx(1.000387, abs=5e-7)  # Multiwfn sums the ten largest
        assert uracil3.deviation == pytest.approx(0.000007, abs=5e-7)
        assert (uracil3.sum_excitation, uracil3.sum_deexcitation) == (pytest.approx(0.500203, abs=5e-7),
                                                                       pytest.approx(0.000196, abs=5e-7))
        assert uracil3.restricted and uracil3.factor == 2.0 and uracil3.expected_norm == 0.5

    def test_cis_weights_are_at_most_one(self):
        r = run("h2co_cis", 3)
        w = r.sets["alpha"].weights
        assert w.max() <= 1.0 + 1e-9
        assert w.sum() == pytest.approx(1.0, abs=abs(r.deviation) * 2 + 1e-6)
        assert not any("de-excitation" in n for n in r.notes)

    def test_td_note_on_deexcitations(self):
        r = run("uracil", 2)
        assert any("de-excitation" in n for n in r.notes)
        assert r.weight_sum() > 1.0

    def test_open_shell_has_two_spin_sets(self):
        r = run("ch3", 1)
        assert not r.restricted and r.factor == 1.0 and r.expected_norm == 1.0
        assert set(r.sets) == {"alpha", "beta"}
        assert r.weight_sum() == pytest.approx(1.0, abs=0.01)
        pairs = r.listed_pairs()
        assert pairs[0].spin == "beta"  # Multiwfn: beta sum 0.98, alpha 0.02
        assert [p.weight for p in pairs] == sorted((p.weight for p in pairs), reverse=True)

    def test_listed_pairs_threshold(self, uracil3):
        pairs = uracil3.listed_pairs()
        assert all(p.weight >= na.LISTED_WEIGHT for p in pairs) and pairs[0].pair == 1
        assert len(uracil3.listed_pairs(minimum=10.0)) == 1  # always at least the largest

    def test_missing_coefficients_give_a_note(self, tmp_path):
        fchk, log = fixture("h2co_cis")
        text = Path(log).read_text().splitlines()
        # keep only the dominant coefficients, as Gaussian prints without IOp(9/40=4)
        kept = [l for l in text if not (("->" in l or "<-" in l) and abs(float(l.split()[-1])) < 0.1)]
        path = tmp_path / "default_print.log"
        path.write_text("\n".join(kept) + "\n")
        r = ANALYSER.run(str(fchk), str(path), 1)
        assert any("IOp(9/40=4)" in n for n in r.notes)

    def test_no_dominant_pair_note(self):
        r = run("benzene", 1)  # two degenerate pairs of ~0.5
        assert any("No NTO pair dominates" in n for n in r.notes)

    def test_unknown_state(self):
        fchk, log = fixture("uracil")
        with pytest.raises(ValueError, match="state 9"):
            ANALYSER.run(str(fchk), str(log), 9)

    def test_mismatched_files(self):
        with pytest.raises(ValueError, match="different atoms"):
            ANALYSER.run(str(fixture("uracil")[0]), str(fixture("h2co_cis")[1]), 1)

    def test_load_is_cached(self):
        a = ANALYSER.load(*map(str, fixture("uracil")))
        b = ANALYSER.load(*map(str, fixture("uracil")))
        assert a[0] is b[0] and a[1] is b[1]


# --------------------------------------------------------------------------- #
# Multiwfn                                                                     #
# --------------------------------------------------------------------------- #
DEGENERACY = 1e-6     # relative weight gap below which pairs form one subspace
COMPARED_WEIGHT = 1e-4


def groups(weights):
    """Indices of consecutive pairs with (relatively) equal weights."""
    out, current = [], [0]
    for k in range(1, len(weights)):
        if abs(weights[k] - weights[current[-1]]) <= DEGENERACY * max(weights[k], 1e-300):
            current.append(k)
        else:
            out.append(current)
            current = [k]
    out.append(current)
    return out


def assert_same_functions(ours, ref, label):
    """Columns of ``ours`` and ``ref`` (values on points) span the same functions, sign/rotation free."""
    scale = np.abs(ref).max()
    r, *_ = np.linalg.lstsq(ours, ref, rcond=None)
    assert np.allclose(r.T @ r, np.eye(r.shape[1]), atol=1e-6), f"{label}: map is not orthogonal"
    assert np.abs(ours @ r - ref).max() <= 1e-6 * scale, f"{label}: residual too large"


def reference_orbitals(ref, spin, n_occ, pairs):
    c = gw.mo_coefficients(ref, spin)
    energies = ref.alpha_orbital_energies if spin == "alpha" else ref.beta_orbital_energies
    holes = [n_occ - k for k in pairs]           # pair k (1-based) -> 0-based MO n_occ - k
    particles = [n_occ + k - 1 for k in pairs]
    return c[:, holes], c[:, particles], energies[holes], energies[particles]


CASES = [(name, state) for name, states in REFERENCE_STATES.items() for state in states]


@pytest.mark.parametrize("name, state", CASES, ids=[f"{n}-S{s}" for n, s in CASES])
def test_ntos_match_multiwfn(name, state, tmp_path):
    results = run(name, state)
    ref = gw.read_fchk(multiwfn_reference(name, state, tmp_path))
    pts = sample_points(results.wavefunction)
    for spin, nto in results.sets.items():
        n_pairs = len(nto.weights)
        energies = ref.alpha_orbital_energies if spin == "alpha" else ref.beta_orbital_energies
        ref_weights = energies[nto.n_occ:nto.n_occ + n_pairs]
        assert np.allclose(nto.weights, ref_weights, atol=2e-8), f"{spin} weights"
        assert np.allclose(energies[nto.n_occ - n_pairs:nto.n_occ][::-1], ref_weights, atol=2e-8)
        compared = [k for k in range(n_pairs) if nto.weights[k] >= COMPARED_WEIGHT]
        for group in groups(nto.weights[:len(compared)]):
            pairs = [k + 1 for k in group]
            h_ref, p_ref, _, _ = reference_orbitals(ref, spin, nto.n_occ, pairs)
            label = f"{name} state {state} {spin} pairs {pairs}"
            for ours, theirs, kind in ((nto.holes[:, group], h_ref, "hole"),
                                       (nto.particles[:, group], p_ref, "particle")):
                assert_same_functions(gw.orbital_values(results.wavefunction, ours, pts),
                                      gw.orbital_values(ref, theirs, pts), f"{label} {kind}")


def test_the_benzene_reference_has_a_degenerate_pair():
    """The subspace branch of the comparison is exercised (pairs 1-2 of benzene S1)."""
    weights = run("benzene", 1).sets["alpha"].weights
    assert groups(weights[:3])[0] == [0, 1]


def test_numerical_module_does_not_import_toga():
    code = "import sys, ntoAnalysis; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
