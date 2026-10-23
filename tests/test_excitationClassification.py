"""Descriptor-based classification of excitations (excitationClassification)."""
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import excitationAnalysis as xa
import excitationClassification as xc
import gaussianWavefunction as gw
import excitation_fixtures as xf

SRC = Path(__file__).resolve().parent.parent

LE = xc.Descriptors(d=0.5, t=-1.0, sr=0.6, delta_sigma=0.1, delta_r=0.8, lambda_=0.6, normalization_deviation=0.0,
                    hole_integral=1.0, electron_integral=1.0, has_diffuse=False, diffuse_electron=0.0, energy_ev=4.0)


class TestRules:
    def test_local_excitation(self):
        c = xc.classify(LE)
        assert c.verdict == xc.VERDICT_LE and c.rule.startswith("Rule 5")
        assert c.confidence.startswith("consistent")

    def test_charge_transfer(self):
        c = xc.classify(replace(LE, d=3.5, t=0.6, delta_r=4.0, lambda_=0.25))
        assert c.verdict == xc.VERDICT_CT and c.rule.startswith("Rule 2")

    def test_ct_boundary_is_inclusive_for_d_and_strict_for_t(self):
        assert xc.classify(replace(LE, d=2.0, t=0.01, delta_r=2.5)).verdict == xc.VERDICT_CT
        assert xc.classify(replace(LE, d=2.5, t=0.0)).verdict == xc.VERDICT_INTERMEDIATE

    def test_rydberg_needs_a_diffuse_electron(self):
        ryd = replace(LE, sr=0.25, delta_sigma=2.0, has_diffuse=True, diffuse_electron=0.8)
        assert xc.classify(ryd).verdict == xc.VERDICT_RYDBERG
        no_basis = xc.classify(replace(ryd, has_diffuse=False, diffuse_electron=0.0))
        assert no_basis.verdict != xc.VERDICT_RYDBERG
        assert any("cannot be assigned in this basis" in w for w in no_basis.warnings)

    def test_multidirectional_ct(self):
        c = xc.classify(replace(LE, delta_sigma=1.5, sr=0.5))
        assert c.verdict == xc.VERDICT_MULTI_CT and "IFCT" in c.confidence

    def test_intermediate_names_the_conflict(self):
        c = xc.classify(replace(LE, d=1.5, t=0.4))
        assert c.verdict == xc.VERDICT_INTERMEDIATE and "partial charge transfer" in c.rule
        assert c.confidence.startswith("none")

    def test_incomplete_coefficients(self):
        c = xc.classify(replace(LE, normalization_deviation=0.03))
        assert c.verdict == xc.VERDICT_INCOMPLETE and "IOp(9/40=4)" in c.rule

    def test_disagreeing_orbital_indices_lower_the_confidence(self):
        c = xc.classify(replace(LE, delta_r=2.4))
        assert c.verdict == xc.VERDICT_LE and c.confidence.startswith("low") and "Delta_r" in c.confidence
        c = xc.classify(replace(LE, lambda_=0.2))
        assert c.confidence.startswith("low") and "Lambda" in c.confidence

    def test_evidence_quotes_values_and_sources(self):
        c = xc.classify(LE)
        text = " ".join(c.evidence)
        assert "D = 0.50 A" in text and "Guido" in text and "Peach" in text and "heuristic" in text

    def test_ghost_and_functional_warnings(self):
        c = xc.classify(replace(LE, d=3.5, t=0.6, delta_r=4.0, energy_ev=1.0, ghost_ev=1.5))
        assert any("necessary, not a sufficient" in w for w in c.warnings)
        assert any("range-separated" in w for w in c.warnings)

    @pytest.mark.parametrize("method,safe", [("RCAM-B3LYP", True), ("UWB97XD", True), ("RM062X", True),
                                             ("RHF", True), ("RB3LYP", False), ("RPBE1PBE", False), ("", False)])
    def test_long_range_safe_methods(self, method, safe):
        assert xc.long_range_safe(method) is safe

    def test_warnings_acknowledge_a_range_separated_functional(self):
        ct = replace(LE, d=3.5, t=0.6, delta_r=4.0, energy_ev=1.0, ghost_ev=1.5)
        safe = xc.classify(replace(ct, method="RCAM-B3LYP")).warnings
        unsafe = xc.classify(replace(ct, method="RB3LYP")).warnings
        assert any("is of that kind" in w for w in safe) and any("not expected" in w for w in safe)
        assert any("not recognised as one" in w for w in unsafe)

    def test_orbital_character_is_never_assigned(self):
        assert "not assigned" in xc.classify(LE).orbital_note


class TestFragments:
    def ifct(self, hole, electron):
        atoms = xa.AtomContributions("mulliken", np.array(hole), np.array(electron), False)
        return xa.ifct(atoms, [xa.Fragment(f"Fragment {k + 1}", (k,)) for k in range(len(hole))])

    def test_ct_assignment_names_the_largest_net_transfer(self):
        text = xc.classify(LE, ifct=self.ifct([0.3, 0.65, 0.05], [0.03, 0.3, 0.67])).fragment_assignment
        assert "Charge transfer" in text and "Fragment 2 -> Fragment 3" in text and "fragment definition" in text

    def test_le_assignment_names_the_main_fragment(self):
        text = xc.classify(LE, ifct=self.ifct([0.9, 0.1], [0.85, 0.15])).fragment_assignment
        assert "Local excitation" in text and "Fragment 1" in text

    def test_no_fragments_no_assignment(self):
        assert xc.classify(LE).fragment_assignment is None


class TestReferences:
    def test_volumes_and_pages(self):
        for key, needle in (("multiwfn", "33, 580"), ("multiwfn", "161, 082503"), ("hole_electron", "165, 461"),
                            ("le_bahers", "7, 2498"), ("delta_r", "9, 3118"), ("lambda", "128, 044118"),
                            ("ghost", "38, 2151"), ("ghost_fix", "154, 204102"), ("hirshfeld", "44, 129"),
                            ("becke", "88, 2547"), ("mulliken", "23, 1833"), ("exciton", "143, 244905")):
            assert needle in xc.REFERENCES[key]

    def test_references_for_keeps_order_and_drops_duplicates(self):
        refs = xc.references_for(["lambda", "multiwfn", "lambda"])
        assert refs == [xc.REFERENCES["multiwfn"], xc.REFERENCES["lambda"]]


class TestDiffuseFunctions:
    @pytest.mark.parametrize("name,expected", [("N-phenylpyrrole", True), ("h2co_ryd", True), ("D-pi-A", False),
                                               ("h2co_cis", False)])
    def test_detection(self, name, expected):
        assert xc.has_diffuse_functions(gw.read_fchk(xf.fixture_paths(name)[0])) is expected


def test_does_not_import_toga():
    code = "import sys, excitationClassification; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
