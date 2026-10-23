"""ExcitationAnalyser.run and the report (excitationReport)."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import excitationAnalysis as xa
import excitationClassification as xc
import excitationReport as xr
import excitation_fixtures as xf

SRC = Path(__file__).resolve().parent.parent
ANALYSER = xa.ExcitationAnalyser()
SECTIONS = ["NORMALIZATION (as Multiwfn)", "MAJOR MO TRANSITIONS", "MO CONTRIBUTIONS TO HOLE AND ELECTRON",
            "HOLE-ELECTRON INDICES", "ORBITAL-BASED INDICES", "TRANSITION ELECTRIC DIPOLE MOMENT",
            "ATOM CONTRIBUTIONS", "HOLE-ELECTRON COULOMB ATTRACTION", "CLASSIFICATION (descriptor-based)",
            "LIMITATIONS", "REFERENCES", "NOTES"]

_CACHE = {}


def run(name, state, **kwargs):
    key = (name, state, tuple(sorted(kwargs.items())))
    if key not in _CACHE:
        fchk, log = xf.fixture_paths(name)
        _CACHE[key] = ANALYSER.run(fchk, log, state, **kwargs)
    return _CACHE[key]


class TestRun:
    def test_progress_is_monotonic_and_ends_at_one(self):
        fchk, log = xf.fixture_paths("h2co_cis")
        seen = []
        xa.ExcitationAnalyser().run(fchk, log, 1, progress=lambda f, msg: seen.append((f, msg)))
        fractions = [f for f, _ in seen]
        assert fractions == sorted(fractions) and fractions[-1] == 1.0 and all(m for _, m in seen)

    def test_a_raising_progress_callback_is_ignored(self):
        fchk, log = xf.fixture_paths("h2co_cis")

        def boom(fraction, message):
            raise RuntimeError("ignored")

        assert xa.ExcitationAnalyser().run(fchk, log, 1, progress=boom).state.number == 1

    def test_bad_fragments_fail_before_heavy_work(self):
        fchk, log = xf.fixture_paths("h2co_cis")
        analyser = xa.ExcitationAnalyser()
        with pytest.raises(ValueError, match="does not exist"):
            analyser.run(fchk, log, 1, fragments_text="1-2; 40")
        assert analyser._geometry is None

    def test_unknown_partition(self):
        fchk, log = xf.fixture_paths("h2co_cis")
        with pytest.raises(ValueError, match="partition"):
            xa.ExcitationAnalyser().run(fchk, log, 1, partition="becke")

    def test_classification_of_the_manual_examples(self):
        assert run("D-pi-A", 1).classification.verdict == xc.VERDICT_LE        # n -> pi* on NO2 (manual 4.18.1)
        assert run("D-pi-A", 2).classification.verdict == xc.VERDICT_CT        # amino -> nitro CT
        assert run("4-nitroaniline", 2).classification.verdict == xc.VERDICT_CT
        assert run("N-phenylpyrrole", 5).classification.verdict == xc.VERDICT_CT
        assert run("h2co_ryd", 3).classification.verdict == xc.VERDICT_RYDBERG
        assert run("h2co_ryd", 4).classification.verdict == xc.VERDICT_RYDBERG

    def test_fragments_give_ifct_and_a_fragment_assignment(self):
        r = run("4-nitroaniline", 2, fragments_text=xf.FRAGMENTS["4-nitroaniline"])
        # exact pairs here: 75.54 %; Multiwfn's |w_i w_j| < 0.001 screening gives its 75.688 %
        # (reproduced in test_excitationContributions with pair_screen=0.001)
        assert r.ifct is not None and abs(r.ifct.intrinsic_ct - 75.54) < 0.01
        assert "Charge transfer between fragments" in r.classification.fragment_assignment

    def test_hirshfeld_falls_back_to_mulliken_with_a_note(self, tmp_path):
        fchk, log = xf.fixture_paths("h2co_cis")
        analyser = xa.ExcitationAnalyser()
        wfn, data, notes = analyser.load(fchk, log)
        wfn.nuclear_charges = wfn.nuclear_charges - np.array([2.0, 2.0, 0.0, 0.0])  # pretend ECPs on C and O
        r = analyser.run(fchk, log, 1, partition=xa.PARTITION_HIRSHFELD)
        assert r.partition == xa.PARTITION_MULLIKEN
        assert any("Hirshfeld partition not possible" in n for n in r.notes)

    def test_diffuse_basis_recommends_hirshfeld(self):
        assert any("Hirshfeld partition is recommended" in n for n in run("N-phenylpyrrole", 1).notes)

    def test_td_state_notes_the_de_excitations(self):
        assert any("de-excitations" in n for n in run("h2co_td5050", 6).notes)

    def test_exact_coulomb_only_on_request(self):
        off = run("h2co_cis", 1)
        assert off.coulomb.exact is None and off.coulomb.mode == xa.COULOMB_OFF
        fast = run("h2co_cis", 1, coulomb_mode=xa.COULOMB_FAST)
        precise = run("h2co_cis", 1, coulomb_mode=xa.COULOMB_PRECISE)
        assert (fast.coulomb.mode, precise.coulomb.mode) == (xa.COULOMB_FAST, xa.COULOMB_PRECISE)
        assert abs(fast.coulomb.exact - precise.coulomb.exact) < 2e-4 * precise.coulomb.exact

    def test_unknown_coulomb_mode(self):
        fchk, log = xf.fixture_paths("h2co_cis")
        with pytest.raises(ValueError, match="Coulomb mode"):
            xa.ExcitationAnalyser().run(fchk, log, 1, coulomb_mode="exact")

    def test_report_names_the_coulomb_mode(self):
        assert "not computed (tick" in xr.report_text(run("h2co_cis", 1))
        assert "fast mode" in xr.report_text(run("h2co_cis", 1, coulomb_mode=xa.COULOMB_FAST))
        assert "precise mode" in xr.report_text(run("h2co_cis", 1, coulomb_mode=xa.COULOMB_PRECISE))


class TestReport:
    def test_sections_in_order(self):
        text = xr.report_text(run("D-pi-A", 2))
        positions = [text.index(s) for s in SECTIONS]
        assert positions == sorted(positions)

    def test_units_and_verdict(self):
        text = xr.report_text(run("D-pi-A", 2))
        assert "D index" in text and " A " in text and "a.u." in text and "eV" in text
        assert f"Verdict: {xc.VERDICT_CT}" in text
        assert "Orbital character: The orbital character" in text

    def test_references_follow_the_analyses(self):
        plain = xr.report_text(run("D-pi-A", 1))
        with_fragments = xr.report_text(run("4-nitroaniline", 2, fragments_text=xf.FRAGMENTS["4-nitroaniline"]))
        assert "J. Comput. Chem. 33, 580" in plain and "Carbon 165, 461" in plain and "9, 3118" in plain
        assert "3.21.8" not in plain.split("REFERENCES")[1]
        assert "3.21.8" in with_fragments.split("REFERENCES")[1]
        assert "Q[R, S]" in with_fragments and "intrinsic CT" in with_fragments

    def test_unrestricted_lists_spin_labelled_mos(self):
        text = xr.report_text(run("ch3", 1))
        block = text.split("MO CONTRIBUTIONS")[1].split("HOLE-ELECTRON INDICES")[0]
        assert "A " in block or "B " in block

    def test_write_report_round_trips(self, tmp_path):
        r = run("h2co_cis", 1)
        path = xa.ExcitationAnalyser.write_report(r, str(tmp_path))
        assert Path(path).name == "h2co_cis_ccpvdz_state1_excitation.txt"
        assert Path(path).read_text(encoding="utf-8") == xr.report_text(r)

    def test_write_failure_is_a_value_error(self, tmp_path):
        blocker = tmp_path / "file"
        blocker.write_text("")
        with pytest.raises(ValueError, match="Could not write"):
            xa.ExcitationAnalyser.write_report(run("h2co_cis", 1), str(blocker / "sub"))


def test_modules_do_not_import_toga():
    code = "import sys, excitationReport; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
