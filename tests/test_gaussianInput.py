"""Tests for the Gaussian 16 input builder core (gaussianInput, no Toga)."""
import re
import subprocess
import sys
from pathlib import Path

import pytest

import fukuiInput as fi
import gaussianInput as gi
import gaussianKeywords as gk
from gaussianInput import RouteKeyword as K

SRC = Path(__file__).resolve().parents[1]


def test_numerical_modules_do_not_import_toga():
    code = "import sys, gaussianInput, gaussianKeywords; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=str(SRC))
    assert out.stdout.strip() == "False", out.stderr


class TestKeywordAlgebra:
    @pytest.mark.parametrize("token,name,options", [
        ("NoSymm", "NoSymm", []),
        ("Opt=Tight", "Opt", ["Tight"]),
        ("Opt=(TS,CalcFC)", "Opt", ["TS", "CalcFC"]),
        ("SCF=(Tight, XQC)", "SCF", ["Tight", "XQC"]),
        ("TD(NStates=4)", "TD", ["NStates=4"]),
        ("TD(50-50,NStates=4)", "TD", ["50-50", "NStates=4"]),
        ("IOp(9/40=4)", "IOp", ["9/40=4"]),
        ("SCRF=(SMD,Solvent=Water)", "SCRF", ["SMD", "Solvent=Water"]),
        ("Integral(Grid=UltraFine)", "Integral", ["Grid=UltraFine"]),
    ])
    def test_parse(self, token, name, options):
        assert gi.parse_keyword(token) == K(name, options)

    @pytest.mark.parametrize("keyword,text", [
        (K("NoSymm"), "NoSymm"),
        (K("Opt", ["Tight"]), "Opt=Tight"),
        (K("Opt", ["TS", "CalcFC"]), "Opt=(TS,CalcFC)"),
        (K("SCRF", ["Solvent=Water"]), "SCRF=(Solvent=Water)"),
        (K("IOp", ["9/40=4"]), "IOp(9/40=4)"),
        (K("TD", ["NStates=10"]), "TD=(NStates=10)"),
    ])
    def test_format(self, keyword, text):
        assert gi.format_keyword(keyword) == text

    def test_key_and_option_key_ignore_case_and_value(self):
        assert K("opt").key == "OPT"
        assert gi.option_key("MaxCycles=50") == "MAXCYCLES" and gi.option_key("tight") == "TIGHT"

    def test_new_keywords_are_appended(self):
        merged, errors = gi.merge_keywords([K("Opt")], "Int=UltraFine 6D")
        assert [gi.format_keyword(k) for k in merged] == ["Opt", "Int=UltraFine", "6D"] and errors == []

    def test_options_of_a_panel_keyword_are_merged(self):
        merged, errors = gi.merge_keywords([K("Opt", ["TS", "CalcFC"])], "Opt=NoEigenTest")
        assert gi.format_keyword(merged[0]) == "Opt=(TS,CalcFC,NoEigenTest)" and errors == []

    def test_a_repeated_option_is_kept_once(self):
        merged, errors = gi.merge_keywords([K("Opt", ["Tight"])], "opt=tight")
        assert merged == [K("Opt", ["Tight"])] and errors == []

    def test_a_valued_option_with_another_value_is_an_error_naming_both_places(self):
        _, errors = gi.merge_keywords([K("Opt", ["MaxCycles=100"])], "Opt=MaxCycles=50")
        assert errors == ["Additional Keywords: Opt=MaxCycles=50 sets MaxCycles=50, but the Job Type panel "
                          "sets MaxCycles=100."]

    def test_unvalued_options_are_not_checked_against_each_other(self):
        merged, errors = gi.merge_keywords([], "Int=UltraFine Int=Fine")
        assert merged == [K("Int", ["UltraFine", "Fine"])] and errors == []

    def test_a_valued_conflict_between_two_additional_keywords(self):
        _, errors = gi.merge_keywords([], "SCF=MaxCycle=64 SCF=MaxCycle=128")
        assert errors == ["Additional Keywords: SCF=MaxCycle=128 sets MaxCycle=128, but an earlier "
                          "Additional Keyword sets MaxCycle=64."]

    def test_the_model_chemistry_is_refused_in_additional_keywords(self):
        _, errors = gi.merge_keywords([], "B3LYP/6-31G")
        assert errors and "Method panel" in errors[0]

    def test_iop_with_a_slash_is_not_a_model_chemistry(self):
        merged, errors = gi.merge_keywords([], "IOp(9/40=4)")
        assert merged == [K("IOp", ["9/40=4"])] and errors == []


from dataclasses import fields, replace


def step(**changes) -> gi.GaussianStepSpec:
    return replace(gi.GaussianStepSpec(), **changes)


def route(**changes) -> str:
    return " ".join(gi.route_lines(step(**changes)))


class TestStepSpec:
    def test_every_field_has_a_label(self):
        assert set(gi.FIELD_LABELS) == {f.name for f in fields(gi.GaussianStepSpec)}

    def test_every_choice_default_is_one_of_its_choices(self):
        default = gi.GaussianStepSpec()
        for name in gi.CHOICES:
            assert getattr(default, name) in gi.choices(name, default), name
        assert default.method_name in gi.choices("method_name", default)

    def test_method_choices_follow_the_family(self):
        assert gi.method_choices(gk.FAMILY_SEMI)[0] == "PM7"
        assert gi.method_choices(gk.FAMILY_HF) == (gk.NO_CHOICE,)


class TestModel:
    def test_default_route(self):
        assert route() == "#P B3LYP/6-31G(d)"

    @pytest.mark.parametrize("changes,model", [
        (dict(basis_diffuse="+", basis_polarization="(d,p)"), "B3LYP/6-31+G(d,p)"),
        (dict(basis="6-311G", basis_diffuse="++", basis_polarization="(2df,2pd)"), "B3LYP/6-311++G(2df,2pd)"),
        (dict(basis="6-31G", basis_polarization="None"), "B3LYP/6-31G"),
        (dict(basis="STO-3G", basis_diffuse="+"), "B3LYP/STO-3G"),          # lists apply to 6-31G/6-311G only
        (dict(basis="aug-cc-pVTZ"), "B3LYP/aug-cc-pVTZ"),
        (dict(basis=gk.OTHER, basis_other="def2-TZVPD"), "B3LYP/def2-TZVPD"),
        (dict(spin="Unrestricted (U)"), "UB3LYP/6-31G(d)"),
        (dict(spin="Restricted-open (RO)", family=gk.FAMILY_HF, method_name=gk.NO_CHOICE), "ROHF/6-31G(d)"),
        (dict(method_name="M06-2X"), "M062X/6-31G(d)"),
        (dict(method_name=gk.OTHER, method_other="revPBE"), "revPBE/6-31G(d)"),
        (dict(family=gk.FAMILY_SEMI, method_name="PM6"), "PM6"),
        (dict(family=gk.FAMILY_COMPOUND, method_name="CBS-QB3"), "CBS-QB3"),
        (dict(family=gk.FAMILY_CCSDT, method_name=gk.NO_CHOICE, basis="cc-pVDZ"), "CCSD(T)/cc-pVDZ"),
        (dict(state=gk.STATE_CIS, n_states=3, basis="cc-pVDZ"), "CIS(NStates=3)/cc-pVDZ"),
        (dict(state=gk.STATE_ZINDO, n_states=10, state_spin="Singlets"), "ZINDO(NStates=10,Singlets)"),
    ])
    def test_model_token(self, changes, model):
        assert gi.model_token(step(**changes)) == model


class TestJobKeywords:
    @pytest.mark.parametrize("changes,expected", [
        (dict(job_type=gk.JOB_OPT), "#P B3LYP/6-31G(d) Opt"),
        (dict(job_type=gk.JOB_OPT, opt_target=gk.OPT_TS, force_constants=gk.FC_CALCFC, opt_tight=True,
              opt_max_cycles=200), "#P B3LYP/6-31G(d) Opt=(TS,CalcFC,Tight,MaxCycles=200)"),
        (dict(job_type=gk.JOB_FREQ, freq_raman=True), "#P B3LYP/6-31G(d) Freq=Raman"),
        (dict(job_type=gk.JOB_OPT_FREQ, freq_vcd=True, freq_anharmonic=True),
         "#P B3LYP/6-31G(d) Opt Freq=(VCD,Anharmonic)"),
        (dict(job_type=gk.JOB_IRC, irc_direction="Forward only", irc_max_points=20, irc_step_size=5),
         "#P B3LYP/6-31G(d) IRC=(Forward,CalcFC,MaxPoints=20,StepSize=5)"),
        (dict(job_type=gk.JOB_SCAN), "#P B3LYP/6-31G(d) Opt=ModRedundant"),
        (dict(job_type=gk.JOB_OPT, modredundant="B 1 2 F"), "#P B3LYP/6-31G(d) Opt=ModRedundant"),
        (dict(job_type=gk.JOB_STABILITY, stable_opt=True), "#P B3LYP/6-31G(d) Stable=Opt"),
        (dict(job_type=gk.JOB_NMR, nmr_method="CSGT", nmr_spinspin=True), "#P B3LYP/6-31G(d) NMR=(CSGT,SpinSpin)"),
        (dict(job_type=gk.JOB_NMR), "#P B3LYP/6-31G(d) NMR"),
    ])
    def test_job_type_route(self, changes, expected):
        assert route(**changes) == expected

    def test_energy_writes_no_job_keyword(self):
        assert gi.job_keywords(step()) == []


class TestOtherPanels:
    @pytest.mark.parametrize("changes,expected", [
        (dict(state=gk.STATE_TD, n_states=10), "#P B3LYP/6-31G(d) TD=(NStates=10)"),
        (dict(state=gk.STATE_TDA, n_states=6, root=2, state_spin="Triplets"),
         "#P B3LYP/6-31G(d) TDA=(NStates=6,Root=2,Triplets)"),
        (dict(state=gk.STATE_TD, print_all_coefficients=True), "#P B3LYP/6-31G(d) TD=(NStates=3) IOp(9/40=4)"),
        (dict(dispersion="GD3BJ"), "#P B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ"),
        (dict(scf_tight=True, scf_qc=True, nosymm=True, polar=True), "#P B3LYP/6-31G(d) SCF=(Tight,QC) NoSymm Polar"),
        (dict(geometry_source=gk.GEOM_ALLCHECK), "#P B3LYP/6-31G(d) Geom=AllCheck"),
        (dict(connectivity=True), "#P B3LYP/6-31G(d) Geom=Connectivity"),
        (dict(guess="Read checkpoint", guess_mix=True), "#P B3LYP/6-31G(d) Guess=(Read,Mix)"),
        (dict(pop_level="Full", pop_mk=True, pop_nbo=True), "#P B3LYP/6-31G(d) Pop=(Full,MK,NBO)"),
        (dict(solvation="SMD", solvent="Acetonitrile"), "#P B3LYP/6-31G(d) SCRF=(SMD,Solvent=Acetonitrile)"),
        (dict(solvation="CPCM", solvent=gk.OTHER, solvent_other="Heptane"),
         "#P B3LYP/6-31G(d) SCRF=(CPCM,Solvent=Heptane)"),
        (dict(print_level="#N"), "#N B3LYP/6-31G(d)"),
    ])
    def test_route(self, changes, expected):
        assert route(**changes) == expected

    def test_additional_keywords_are_merged(self):
        assert route(job_type=gk.JOB_OPT, opt_target=gk.OPT_TS, additional_keywords="Opt=NoEigenTest 6D") == \
            "#P B3LYP/6-31G(d) Opt=(TS,NoEigenTest) 6D"

    def test_route_keywords_report_merge_errors(self):
        _, errors = gi.route_keywords(step(job_type=gk.JOB_OPT, opt_max_cycles=100,
                                           additional_keywords="Opt=MaxCycles=50"))
        assert errors and "Job Type panel" in errors[0]

    def test_a_long_route_is_wrapped_below_80_columns(self):
        lines = gi.route_lines(step(job_type=gk.JOB_OPT_FREQ, opt_tight=True, freq_raman=True,
                                    state=gk.STATE_TD, n_states=10, print_all_coefficients=True,
                                    dispersion="GD3BJ", scf_tight=True, nosymm=True, pop_level="Full",
                                    solvation="SMD", solvent="DiMethylSulfoxide"))
        assert len(lines) > 1 and all(len(line) <= fi.ROUTE_WIDTH for line in lines)


class TestModRedundant:
    def test_scan_and_freeze_lines(self):
        assert gi.parse_modredundant("B 1 2 S 10 0.05\n\nd 1 2 3 4 f\nA 1 2 3 F", n_atoms=4) == [
            ("B", (1, 2), "S"), ("D", (1, 2, 3, 4), "F"), ("A", (1, 2, 3), "F")]

    @pytest.mark.parametrize("text,fragment", [
        ("X 1 2 F", "start with B, A or D"),
        ("B 1 F", "needs 2 atom numbers"),
        ("A 1 1 2 F", "different atoms"),
        ("B 0 2 F", "different atoms"),
        ("B 1 9 F", "atom 9 does not exist (the molecule has 3 atoms)"),
        ("B 1 2", "end with F (freeze) or S <steps> <step size> (scan)"),
        ("B 1 2 S 0 0.1", "steps must be a positive integer"),
        ("B 1 2 S 5 0", "step size must be a non-zero number"),
        ("B 1 2 S 5 x", "step size must be a non-zero number"),
    ])
    def test_problems_name_the_line(self, text, fragment):
        with pytest.raises(ValueError, match=re.escape(fragment)):
            gi.parse_modredundant("B 1 2 F\n" + text, n_atoms=3)
        with pytest.raises(ValueError, match="line 2"):
            gi.parse_modredundant("B 1 2 F\n" + text, n_atoms=3)

    def test_atom_range_is_not_checked_without_a_molecule(self):
        assert gi.parse_modredundant("B 1 99 F") == [("B", (1, 99), "F")]

    def test_helper_row_builds_a_line(self):
        assert gi.modredundant_line("Dihedral", "1, 2 3 4", "Scan", "36", "10.0") == "D 1 2 3 4 S 36 10.0"
        assert gi.modredundant_line("Bond", "1 2", "Freeze", "", "") == "B 1 2 F"

    def test_helper_row_refuses_what_the_parser_refuses(self):
        with pytest.raises(ValueError, match="needs 3 atom numbers"):
            gi.modredundant_line("Angle", "1 2", "Freeze", "", "")
