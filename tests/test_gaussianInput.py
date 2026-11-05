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


class TestChkBasis:
    def test_route(self):
        assert route(basis=gk.BASIS_CHECKPOINT) == "#P B3LYP ChkBasis"
        assert route(basis=gk.BASIS_CHECKPOINT, state=gk.STATE_TD) == "#P B3LYP ChkBasis TD=(NStates=3)"

    def test_a_method_without_a_basis_writes_no_chkbasis(self):
        assert route(basis=gk.BASIS_CHECKPOINT, family=gk.FAMILY_SEMI, method_name="PM7") == "#P PM7"

    def test_it_reads_a_checkpoint(self, tmp_path):
        text = problems(tmp_path, step(basis=gk.BASIS_CHECKPOINT, multiplicity=1, write_chk=False))
        assert "ChkBasis reads a checkpoint" in text
        gi.GaussianInputBuilder().prepare(job_spec(tmp_path, step(basis=gk.BASIS_CHECKPOINT, multiplicity=1,
                                                                    old_chk="old.chk")))


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


WATER = ("O 0.000000 0.000000 0.117300\n"
         "H 0.000000 0.757200 -0.469200\n"
         "H 0.000000 -0.757200 -0.469200")
WATER_LINES = fi.format_geometry(fi.parse_geometry(WATER))


def job_spec(tmp_path, *steps, **changes) -> gi.GaussianJobSpec:
    steps = list(steps) or [step(multiplicity=1)]
    return replace(gi.GaussianJobSpec(stem="water", output_dir=str(tmp_path), geometry=WATER, steps=steps),
                   **changes)


def problems(tmp_path, *steps, **changes) -> str:
    with pytest.raises(ValueError) as info:
        gi.GaussianInputBuilder().prepare(job_spec(tmp_path, *steps, **changes))
    return str(info.value)


class TestInputText:
    def test_single_step(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path, step(job_type=gk.JOB_OPT_FREQ, multiplicity=1)))
        assert builder.input_text(job) == "\n".join(
            ["%Chk=water.chk", "%Mem=2GB", "%NProcShared=4", "#P B3LYP/6-31G(d) Opt Freq", "",
             "water", "", "0 1", *WATER_LINES, ""]) + "\n"
        assert job.input_path == tmp_path / "water.gjf"
        assert job.script_path == tmp_path / "water_formchk.sh" and job.chk_files == ["water.chk"]

    def test_linked_opt_then_td(self, tmp_path):
        first = step(job_type=gk.JOB_OPT, multiplicity=1)
        second = replace(gi.linked_step(first), state=gk.STATE_TD, n_states=10, print_all_coefficients=True)
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path, first, second))
        assert builder.input_text(job) == "\n".join(
            ["%Chk=water.chk", "%Mem=2GB", "%NProcShared=4", "#P B3LYP/6-31G(d) Opt", "",
             "water", "", "0 1", *WATER_LINES, "",
             "--Link1--",
             "%Chk=water.chk", "%Mem=2GB", "%NProcShared=4",
             "#P B3LYP/6-31G(d) TD=(NStates=10) IOp(9/40=4) Geom=AllCheck Guess=Read", ""]) + "\n"
        assert any("same checkpoint" in w for w in job.warnings)
        assert any("Density=Current" in w for w in job.warnings)

    def test_check_step_writes_title_and_charge_but_no_coordinates(self, tmp_path):
        second = step(geometry_source=gk.GEOM_CHECK, charge=1, multiplicity=2, guess="Read checkpoint")
        builder = gi.GaussianInputBuilder()
        text = builder.input_text(builder.prepare(job_spec(tmp_path, step(multiplicity=1), second)))
        block = text.split("--Link1--\n")[1]
        assert block.endswith("Geom=Check Guess=Read\n\nwater\n\n1 2\n\n")

    def test_sections_after_the_geometry(self, tmp_path):
        s = step(job_type=gk.JOB_SCAN, multiplicity=1, modredundant="B 1 2 S 5 0.05", connectivity=True,
                 basis="Gen", additional_input="O H 0\n6-31G(d)\n****")
        builder = gi.GaussianInputBuilder()
        text = builder.input_text(builder.prepare(job_spec(tmp_path, s)))
        assert text.endswith("\n".join([*WATER_LINES, "", " 1 2 1.0 3 1.0", " 2", " 3", "",
                                        "B 1 2 S 5 0.05", "", "O H 0", "6-31G(d)", "****", ""]) + "\n")
        assert "Geom=Connectivity" in text and "Opt=ModRedundant" in text

    def test_typed_title_and_link0(self, tmp_path):
        s = step(multiplicity=1, title="my water", chk="w1.chk", old_chk="start.chk", memory="8gb", nproc=16,
                 extra_link0="%KJob L502")
        builder = gi.GaussianInputBuilder()
        text = builder.input_text(builder.prepare(job_spec(tmp_path, s)))
        assert text.startswith("%OldChk=start.chk\n%Chk=w1.chk\n%Mem=8GB\n%NProcShared=16\n%KJob L502\n")
        assert "\nmy water\n" in text

    def test_no_checkpoint_means_no_script(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, step(multiplicity=1, write_chk=False)))
        assert job.script_path is None and any("no formchk script" in w for w in job.warnings)

    def test_files_are_written_with_lf_endings(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path))
        assert builder.write_files(job) == [job.input_path, job.script_path]
        assert b"\r" not in job.input_path.read_bytes()
        script = job.script_path.read_text()
        assert 'convert "water.chk" "water.fchk"' in script and "3D Molecular Viewer" in script

    def test_unencodable_text_leaves_an_existing_file_untouched(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path))
        job.steps[0].title = "água"
        job.input_path.write_text("keep")
        with pytest.raises(ValueError, match="Could not write"):
            builder.write_files(job)
        assert job.input_path.read_text() == "keep"

    def test_write_error_is_a_value_error(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path))
        job.input_path = tmp_path / "missing" / "x.gjf"
        with pytest.raises(ValueError, match="Could not write"):
            builder.write_files(job)


class TestValidation:
    @pytest.mark.parametrize("changes,fragment", [
        (dict(multiplicity=None), "Step 1: Multiplicity: enter it; it is not assumed (10 electrons allow 1, 3, 5, ...)"),
        (dict(multiplicity=2), "Multiplicity 2 is impossible"),
        (dict(multiplicity=3, spin="Restricted (R)"), "R (restricted) prefix needs multiplicity 1"),
        (dict(multiplicity=1, method_name=gk.OTHER, method_other="B3 LYP"), "Other method: enter one keyword"),
        (dict(multiplicity=1, basis=gk.OTHER, basis_other=""), "Other basis set: enter one keyword"),
        (dict(multiplicity=1, basis="GenECP"), "Basis set GenECP: put the basis-set section in the Add. Inp. panel"),
        (dict(multiplicity=1, state=gk.STATE_TD, n_states=0), "Number of states must be at least 1"),
        (dict(multiplicity=1, state=gk.STATE_TD, n_states=3, root=4), "Root) must be between 1 and 3"),
        (dict(multiplicity=1, state=gk.STATE_TD, family=gk.FAMILY_MP2, method_name=gk.NO_CHOICE),
         "TD-SCF needs a Hartree-Fock or DFT method"),
        (dict(multiplicity=1, state=gk.STATE_TDA, family=gk.FAMILY_HF, method_name=gk.NO_CHOICE), "TDA needs a DFT method"),
        (dict(multiplicity=1, state=gk.STATE_TD, state_spin="Singlets", spin="Unrestricted (U)"),
         "need a closed-shell reference"),
        (dict(multiplicity=1, print_all_coefficients=True), "needs an excited state"),
        (dict(multiplicity=1, title="a#b"), "Title: Gaussian does not allow #"),
        (dict(multiplicity=1, chk="C:\\jobs\\w.chk"), "%Chk 'C:\\jobs\\w.chk': give a file name relative"),
        (dict(multiplicity=1, chk="my w.chk"), "%Chk 'my w.chk': use a file name without blanks"),
        (dict(multiplicity=1, extra_link0="KJob L502"), "must start with %"),
        (dict(multiplicity=1, memory="2"), "%Mem = '2'"),
        (dict(multiplicity=1, guess="Read checkpoint"), "Guess=Read reads a checkpoint, but no earlier step writes"),
        (dict(multiplicity=1, geometry_source=gk.GEOM_ALLCHECK), "Geom=AllCheck reads a checkpoint"),
        (dict(multiplicity=1, job_type=gk.JOB_OPT, force_constants="Read from checkpoint (ReadFC)"),
         "Opt=ReadFC reads a checkpoint"),
        (dict(multiplicity=1, geometry_source=gk.GEOM_CHECK, connectivity=True, old_chk="a.chk"),
         "Write connectivity needs coordinates"),
        (dict(multiplicity=1, modredundant="B 1 2 F"), "Scan/freeze lines need the Optimization"),
        (dict(multiplicity=1, job_type=gk.JOB_OPT, modredundant="B 1 2 S 5 0.1"), "Scan lines (S) need the Scan job type"),
        (dict(multiplicity=1, job_type=gk.JOB_SCAN), "The Scan job type needs at least one scan line"),
        (dict(multiplicity=1, job_type=gk.JOB_OPT, modredundant="B 1 7 F"), "atom 7 does not exist"),
        (dict(multiplicity=1, solvation="SMD", solvent=gk.OTHER, solvent_other="n,n"), "Other solvent"),
        (dict(multiplicity=1, additional_keywords="B3LYP/STO-3G"), "Method panel"),
        (dict(multiplicity=1, job_type="Dance"), "Job type: 'Dance' is not one of the choices"),
    ])
    def test_step_problem(self, tmp_path, changes, fragment):
        assert fragment in problems(tmp_path, step(**changes))

    @pytest.mark.parametrize("changes,label", [
        (dict(title="água"), "Title"),
        (dict(extra_link0="%Chk=á.chk"), "Extra Link 0 lines"),
        (dict(additional_input="Ω"), "Additional input"),
        (dict(additional_keywords="Int=ÚltraFine"), "Additional Keywords"),
        (dict(method_name=gk.OTHER, method_other="B3ßLYP"), "Other method"),
        (dict(basis=gk.OTHER, basis_other="é"), "Other basis set"),
        (dict(solvation="SMD", solvent=gk.OTHER, solvent_other="água"), "Other solvent"),
    ])
    def test_non_ascii_text_is_refused(self, tmp_path, changes, label):
        text = problems(tmp_path, step(multiplicity=1, **changes))
        assert f"{label}: use ASCII characters only (Gaussian input files are ASCII)." in text

    def test_job_problems(self, tmp_path):
        text = problems(tmp_path, step(multiplicity=1), stem="bad name", output_dir=str(tmp_path / "nope"),
                        geometry="")
        assert "Job name 'bad name'" in text and "does not exist" in text and "geometry is empty" in text

    def test_every_problem_is_listed_on_its_own_line(self, tmp_path):
        lines = problems(tmp_path, step(multiplicity=None, memory="x"), step(multiplicity=None)).splitlines()
        assert len(lines) >= 3 and all(line.startswith("- ") for line in lines)
        assert any(line.startswith("- Step 2: ") for line in lines)

    def test_an_old_checkpoint_lets_step_one_read(self, tmp_path):
        s = step(geometry_source=gk.GEOM_ALLCHECK, old_chk="start.chk")
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, s, geometry=""))
        assert job.atoms is None

    def test_a_later_step_reads_the_checkpoint_an_earlier_step_writes(self, tmp_path):
        first = step(multiplicity=1, chk="a.chk")
        second = step(geometry_source=gk.GEOM_ALLCHECK, chk="a.chk", guess="Read checkpoint")
        gi.GaussianInputBuilder().prepare(job_spec(tmp_path, first, second))
        with pytest.raises(ValueError, match="no earlier step writes b.chk"):
            gi.GaussianInputBuilder().prepare(job_spec(tmp_path, first, replace(second, chk="b.chk")))

    def test_ts_without_force_constants_warns(self, tmp_path):
        ts = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, step(multiplicity=1, job_type=gk.JOB_OPT,
                                                                       opt_target=gk.OPT_TS)))
        assert any("TS" in w for w in ts.warnings)

    def test_irc_without_force_constants_is_refused(self, tmp_path):
        # Gaussian 16 C.01: "You must use ReadCartesianFC, FCCards, or CalcFC with IRC."
        with pytest.raises(ValueError, match="IRC needs force constants: give CalcFC, CalcAll or RCFC"):
            gi.GaussianInputBuilder().prepare(job_spec(tmp_path, step(multiplicity=1, job_type=gk.JOB_IRC,
                                                                      irc_force_constants=gk.FC_DEFAULT)))

    def test_irc_with_fccards_prepares(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, step(
            multiplicity=1, job_type=gk.JOB_IRC, irc_force_constants=gk.FC_DEFAULT, additional_keywords="IRC=FCCards")))
        assert "FCCards" in gi.GaussianInputBuilder().input_text(job)

    def test_linked_step_drops_job_keywords_from_additional_keywords(self):
        first = step(multiplicity=1, additional_keywords="Opt=NoEigenTest Int=UltraFine Freq=Raman 6D")
        assert gi.linked_step(first).additional_keywords == "Int=UltraFine 6D"

    def test_default_job_has_no_warnings(self, tmp_path):
        assert gi.GaussianInputBuilder().prepare(job_spec(tmp_path)).warnings == []


TS_FREQ = step(multiplicity=1, job_type=gk.JOB_OPT_FREQ, opt_target=gk.OPT_TS, force_constants=gk.FC_CALCFC)
TS_ONLY = replace(TS_FREQ, job_type=gk.JOB_OPT)
RCFC_WARNING = "Step 2: RCFC needs force constants from a frequency step earlier in the job"


def irc_after(first, **changes) -> gi.GaussianStepSpec:
    return replace(gi.linked_step(first), job_type=gk.JOB_IRC, **changes)


class TestIrcForceConstants:
    """Gaussian 16 C.01 refuses IRC=ReadFC; IRC=RCFC works after a TS step that ran Freq."""

    def test_rcfc_after_a_freq_step_prepares_without_a_warning(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path, TS_FREQ, irc_after(TS_FREQ, irc_force_constants=gk.FC_RCFC)))
        assert "IRC=RCFC Geom=AllCheck Guess=Read" in builder.input_text(job)
        assert not any("RCFC needs" in w for w in job.warnings)

    def test_rcfc_without_an_earlier_freq_step_warns(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, TS_ONLY,
                                                         irc_after(TS_ONLY, irc_force_constants=gk.FC_RCFC)))
        assert any(w.startswith(RCFC_WARNING) for w in job.warnings)

    def test_read_cartesian_fc_is_accepted_and_warns_the_same_way(self, tmp_path):
        typed = dict(irc_force_constants=gk.FC_DEFAULT, additional_keywords="IRC=ReadCartesianFC")
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path, TS_FREQ, irc_after(TS_FREQ, **typed)))
        assert "IRC=ReadCartesianFC" in builder.input_text(job) and not any("RCFC needs" in w for w in job.warnings)
        job = builder.prepare(job_spec(tmp_path, TS_ONLY, irc_after(TS_ONLY, **typed)))
        assert any(w.startswith(RCFC_WARNING) for w in job.warnings)

    def test_rcfc_reads_the_checkpoint(self, tmp_path):
        text = problems(tmp_path, step(multiplicity=1, job_type=gk.JOB_IRC, irc_force_constants=gk.FC_RCFC))
        assert "IRC=RCFC reads a checkpoint" in text

    def test_typed_readfc_is_refused_with_its_own_message(self, tmp_path):
        text = problems(tmp_path, TS_FREQ, irc_after(TS_FREQ, irc_force_constants=gk.FC_DEFAULT,
                                                     additional_keywords="IRC=ReadFC"))
        assert "Step 2: IRC=ReadFC: Gaussian 16 does not accept ReadFC for IRC" in text
        assert "IRC needs force constants" not in text

    def test_readfc_is_no_irc_choice_but_stays_an_opt_choice(self):
        assert not any("ReadFC" in c for c in gi.choices("irc_force_constants"))
        assert "Read from checkpoint (ReadFC)" in gi.choices("force_constants")

    def test_reopened_irc_rcfc_restores_the_control(self):
        spec, _ = gi.read_input_text("# hf/sto-3g irc=(rcfc,maxpoints=5)\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n")
        s = spec.steps[0]
        assert (s.job_type, s.irc_force_constants, s.irc_max_points) == (gk.JOB_IRC, gk.FC_RCFC, 5)
        assert s.additional_keywords == ""

    def test_reopened_irc_readfc_is_refused(self, tmp_path):
        text = (f"%Chk=water.chk\n# B3LYP/6-31G(d) Opt=(TS,CalcFC) Freq\n\nwater\n\n0 1\n{WATER}\n\n"
                "--Link1--\n%Chk=water.chk\n# B3LYP/6-31G(d) IRC=ReadFC Geom=AllCheck Guess=Read\n\n")
        spec, _ = gi.read_input_text(text, "water", str(tmp_path))
        assert spec.steps[1].irc_force_constants == gk.FC_DEFAULT and "ReadFC" in spec.steps[1].additional_keywords
        with pytest.raises(ValueError, match="Gaussian 16 does not accept ReadFC for IRC"):
            gi.GaussianInputBuilder().prepare(spec)


class TestFinalReviewFixes:
    def test_shared_checkpoint_note_only_with_a_formchk_script(self, tmp_path):
        first = step(job_type=gk.JOB_OPT, multiplicity=1)
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, first, gi.linked_step(first), write_formchk=False))
        assert not any("same checkpoint" in w for w in job.warnings)
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, first, gi.linked_step(first)))
        assert any("same checkpoint" in w and "formchk script" in w for w in job.warnings)

    def test_density_warning_only_when_the_step_writes_a_checkpoint(self, tmp_path):
        td = step(multiplicity=1, state=gk.STATE_TD)
        assert any("Density=Current" in w for w in gi.GaussianInputBuilder().prepare(job_spec(tmp_path, td)).warnings)
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, replace(td, write_chk=False)))
        assert not any("Density=Current" in w for w in job.warnings)

    @pytest.mark.parametrize("changes", [dict(method_other="B3ßLYP"), dict(basis_other="é"),
                                         dict(solvation="SMD", solvent_other="água")])
    def test_unused_other_fields_are_not_checked_for_ascii(self, tmp_path, changes):
        gi.GaussianInputBuilder().prepare(job_spec(tmp_path, step(multiplicity=1, **changes)))

    @pytest.mark.parametrize("changes,label", [
        (dict(job_type=gk.JOB_OPT, modredundant="B 1 2 F é"), "Scan/freeze lines"),
        (dict(memory="2µGB"), "%Mem"),
    ])
    def test_more_fields_are_checked_for_ascii(self, tmp_path, changes, label):
        text = problems(tmp_path, step(multiplicity=1, **changes))
        assert f"{label}: use ASCII characters only (Gaussian input files are ASCII)." in text

    @pytest.mark.parametrize("name", ["a*b.chk", "a;b.chk", "a$b.chk", "a'b.chk"])
    def test_checkpoint_names_keep_to_the_listed_characters(self, tmp_path, name):
        text = problems(tmp_path, step(multiplicity=1, chk=name))
        assert f"%Chk '{name}': use a file name without blanks (ASCII letters, digits, . _ -)." in text

    def test_checkpoint_chain_is_case_sensitive(self, tmp_path):
        first = step(multiplicity=1, chk="A.chk")
        second = step(geometry_source=gk.GEOM_ALLCHECK, chk="a.chk", guess="Read checkpoint")
        assert "no earlier step writes a.chk" in problems(tmp_path, first, second)
        gi.GaussianInputBuilder().prepare(job_spec(tmp_path, first, replace(second, chk="A.chk")))


class TestSteps:
    def test_linked_step_continues_from_the_checkpoint(self):
        first = step(job_type=gk.JOB_OPT_FREQ, multiplicity=1, title="t", modredundant="B 1 2 F",
                     additional_input="x", memory="8GB", chk="a.chk", old_chk="o.chk")
        new = gi.linked_step(first)
        assert (new.job_type, new.geometry_source, new.guess) == (gk.JOB_ENERGY, gk.GEOM_ALLCHECK, "Read checkpoint")
        assert (new.title, new.modredundant, new.additional_input, new.old_chk) == ("", "", "", "")
        assert (new.memory, new.chk, new.family, new.method_name) == ("8GB", "a.chk", gk.FAMILY_DFT, "B3LYP")
        assert (new.charge, new.multiplicity) == (0, None)

    def test_step_label(self):
        assert gi.step_label(1, step(job_type=gk.JOB_OPT)) == "1: Optimization B3LYP/6-31G(d)"
        assert gi.step_label(2, step(geometry_source=gk.GEOM_ALLCHECK, state=gk.STATE_TD)) == \
            "2: TD-SCF B3LYP/6-31G(d), Geom=AllCheck"


PHENOL_GJF = SRC / "tests" / "data" / "fukui" / "inputs" / "phenol_N.gjf"


def phenol_geometry() -> str:
    lines = PHENOL_GJF.read_text().splitlines()
    start = lines.index("0 1") + 1
    return "\n".join(line for line in lines[start:] if line.strip())


def fukui_job_spec(tmp_path, **changes) -> gi.GaussianJobSpec:
    s = gi.fukui_defaults(step(multiplicity=1, fukui_plus_multiplicity=2, fukui_minus_multiplicity=2))
    return gi.GaussianJobSpec(stem="phenol", output_dir=str(tmp_path), geometry=phenol_geometry(),
                              steps=[replace(s, **changes)])


class TestFukui:
    def test_output_is_the_fukui_builders_byte_for_byte(self, tmp_path):
        old = fi.FukuiInputBuilder()
        old_job = old.prepare(fi.FukuiJobSpec(stem="phenol", output_dir=str(tmp_path), geometry=phenol_geometry(),
                                              charge=0, multiplicity=1, plus_multiplicity=2, minus_multiplicity=2,
                                              memory="2GB", nproc=4))
        new = gi.GaussianInputBuilder()
        job = new.prepare(fukui_job_spec(tmp_path))
        assert new.input_text(job) == old.input_text(old_job)
        assert new.script_text(job) == old.script_text(old_job)
        assert (job.input_path, job.script_path) == (old_job.input_path, old_job.script_path)
        assert new.write_files(job) == [old_job.input_path, old_job.script_path]

    def test_panels_map_onto_the_fukui_spec(self, tmp_path):
        spec = fukui_job_spec(tmp_path, method_name="M06-2X", basis="Def2SVP", dispersion="GD3BJ", pop_mk=True,
                              solvation="SMD", force_constants=gk.FC_CALCFC, opt_tight=True, opt_max_cycles=50,
                              additional_keywords="Int=UltraFine", print_level="#N", memory="3GB", nproc=8,
                              charge=1, multiplicity=2, fukui_p=2, fukui_q=1, fukui_plus_multiplicity=3,
                              fukui_minus_multiplicity=1, fukui_guess_read=False)
        f = gi.fukui_spec(spec)
        assert (f.method, f.basis, f.print_level, f.memory, f.nproc) == ("M062X", "Def2SVP", "#N", "3GB", 8)
        assert f.keywords == "EmpiricalDispersion=GD3BJ SCF=Tight NoSymm Pop=MK SCRF=(SMD,Solvent=Water) Int=UltraFine"
        assert f.opt_keyword == "Opt=(CalcFC,Tight,MaxCycles=50)"
        assert (f.charge, f.multiplicity, f.plus_multiplicity, f.minus_multiplicity) == (1, 2, 3, 1)
        assert (f.p, f.q, f.guess_read, f.stem, f.geometry) == (2, 1, False, "phenol", phenol_geometry())

    def test_fukui_must_be_the_only_step(self, tmp_path):
        spec = fukui_job_spec(tmp_path)
        spec.steps.append(step(multiplicity=1))
        with pytest.raises(ValueError, match="must be the only step"):
            gi.GaussianInputBuilder().prepare(spec)

    def test_a_checkpoint_force_constant_option_names_the_allowed_ones(self, tmp_path):
        with pytest.raises(ValueError) as info:
            gi.GaussianInputBuilder().prepare(fukui_job_spec(tmp_path, force_constants="Read from checkpoint (ReadFC)"))
        text = str(info.value)
        assert "Force constants must be Default, Calculate at first point (CalcFC) or Calculate at all points " \
               "(CalcAll)" in text and "Read from checkpoint (ReadFC)" in text and "choose the job type" not in text

    @pytest.mark.parametrize("changes,label", [
        (dict(guess="Read checkpoint"), "Initial guess"),
        (dict(state=gk.STATE_TD), "Electronic state"),
        (dict(title="x"), "Title"),
    ])
    def test_controls_the_fukui_job_sets_itself_are_refused(self, tmp_path, changes, label):
        with pytest.raises(ValueError, match=f"sets {label} itself"):
            gi.GaussianInputBuilder().prepare(fukui_job_spec(tmp_path, **changes))

    def test_additional_keywords_merge_and_opt_is_refused_with_the_new_wording(self, tmp_path):
        assert gi.fukui_spec(fukui_job_spec(tmp_path, additional_keywords="SCF=XQC")).keywords == "SCF=(Tight,XQC) NoSymm"
        with pytest.raises(ValueError, match="Optimization options of the Job Type panel"):
            gi.GaussianInputBuilder().prepare(fukui_job_spec(tmp_path, additional_keywords="Opt=Tight"))

    def test_fukui_problems_come_from_fukui_input(self, tmp_path):
        with pytest.raises(ValueError, match=r"N\+1 state: enter its spin multiplicity"):
            gi.GaussianInputBuilder().prepare(fukui_job_spec(tmp_path, fukui_plus_multiplicity=None))

    def test_fukui_preview_shows_conflicting_additional_keywords(self, tmp_path):
        text = gi.route_preview(fukui_job_spec(tmp_path, additional_keywords="Opt=MaxCycles=5 Opt=MaxCycles=7"), 0)
        assert text.splitlines()[0].startswith("job 1:")
        assert "earlier Additional Keyword" in text

    def test_fukui_defaults(self):
        s = gi.fukui_defaults(step(guess="Read checkpoint", title="t", spin="Unrestricted (U)", chk="a.chk",
                                   force_constants="Read from checkpoint (ReadFC)"))
        assert s.job_type == gk.JOB_FUKUI and s.scf_tight and s.nosymm
        assert (s.guess, s.title, s.spin, s.chk, s.force_constants) == ("Default", "", "Default", "", gk.FC_DEFAULT)


class TestPreview:
    def test_fukui_route_preview_matches_the_former_window(self, tmp_path):
        assert gi.route_preview(fukui_job_spec(tmp_path), 0) == (
            "job 1:     #P B3LYP/6-31G(d) SCF=Tight NoSymm Opt\n"
            "jobs 2-3:  #P B3LYP/6-31G(d) SCF=Tight NoSymm Geom=Check Guess=Read")

    def test_route_preview_shows_merge_problems(self, tmp_path):
        spec = job_spec(tmp_path, step(job_type=gk.JOB_OPT, opt_max_cycles=9, additional_keywords="Opt=MaxCycles=5"))
        text = gi.route_preview(spec, 0)
        assert text.startswith("#P B3LYP/6-31G(d) Opt=(MaxCycles=9)") and "Job Type panel" in text

    def test_molecule_summary(self, tmp_path):
        assert gi.molecule_summary(job_spec(tmp_path, step()), 0) == \
            "3 atoms, charge 0, 10 electrons, multiplicity ? (enter it: 10 electrons allow 1, 3, 5, ...)"
        assert gi.molecule_summary(job_spec(tmp_path, step(charge=1, multiplicity=1)), 0).endswith(
            "multiplicity 1   <- impossible: " + gi.multiplicity_problem(9, 1))
        assert "checkpoint" in gi.molecule_summary(job_spec(tmp_path, step(geometry_source=gk.GEOM_ALLCHECK)), 0)
        assert "geometry" in gi.molecule_summary(job_spec(tmp_path, geometry=""), 0)

    def test_fukui_summary_lists_the_three_states(self, tmp_path):
        text = gi.molecule_summary(fukui_job_spec(tmp_path, fukui_plus_multiplicity=None), 0)
        assert "N+1" in text and "charge -1" in text and "51 electrons" in text and "phenol_Np1.chk" in text
        assert "multiplicity ? (enter it: 51 electrons allow 2, 4, 6, ...)" in text


import math

from fukui_fixtures import N

FUKUI_INPUTS = SRC / "tests" / "data" / "fukui" / "inputs"


@pytest.mark.parametrize("windows,linux", [(r"C:\jobs\my mol", "/mnt/c/jobs/my mol"),
                                           (r"D:\x", "/mnt/d/x"), ("/home/u/jobs", "/home/u/jobs")])
def test_wsl_path(windows, linux):
    assert gi.wsl_path(windows) == linux


class TestGeometryFiles:
    def test_gaussian_input_gives_charge_and_multiplicity(self):
        text, charge, mult = gi.read_geometry_file(FUKUI_INPUTS / "phenol_Np1.gjf")
        assert (charge, mult) == (-1, 2) and len(fi.parse_geometry(text)) == 13

    def test_fchk_gives_angstrom_coordinates(self):
        text, charge, mult = gi.read_geometry_file(N)
        atoms = fi.parse_geometry(text)
        reference = fi.parse_geometry(gi.read_geometry_file(FUKUI_INPUTS / "phenol_N.gjf")[0])
        assert (charge, mult) == (0, 1) and [a[0] for a in atoms] == [a[0] for a in reference]
        assert math.dist(atoms[0][1:], atoms[1][1:]) == pytest.approx(math.dist(reference[0][1:],
                                                                               reference[1][1:]), abs=1e-5)

    def test_xyz_has_no_charge(self, tmp_path):
        path = tmp_path / "w.xyz"
        path.write_text("3\nwater\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692\n")
        text, charge, mult = gi.read_geometry_file(path)
        assert charge is None and mult is None and len(fi.parse_geometry(text)) == 3

    OUTPUTS = SRC / "tests" / "data" / "gaussian_output"

    @pytest.mark.parametrize("name", ["h2o_opt.log", "h2o_opt_freq.out"])
    def test_gaussian_output_gives_the_optimized_geometry(self, name):
        g = gi.load_geometry_file(self.OUTPUTS / name)
        assert (g.charge, g.multiplicity) == (0, 1) and len(fi.parse_geometry(g.text)) == 3
        assert g.source.startswith("the optimized geometry of job 1 of") and g.notes == []
        assert gi.read_geometry_file(self.OUTPUTS / name) == (g.text, 0, 1)

    def test_gaussian_output_notes(self):
        g = gi.load_geometry_file(self.OUTPUTS / "h2o_opt_then_failed.log")
        assert len(g.notes) == 1 and g.source == "the optimized geometry of job 1 of 2"

    def test_unconverged_output_is_refused(self):
        with pytest.raises(ValueError, match="no successfully optimized geometry"):
            gi.load_geometry_file(self.OUTPUTS / "h2co_opt_maxcycles.log")

    def test_other_files_have_no_notes(self):
        g = gi.load_geometry_file(FUKUI_INPUTS / "phenol_Np1.gjf")
        assert g.notes == [] and g.source == "" and (g.charge, g.multiplicity) == (-1, 2)

    def test_output_suffixes_are_offered(self):
        assert {"log", "out"} <= set(gi.GEOMETRY_FILE_TYPES)

    def test_unsupported_file(self, tmp_path):
        path = tmp_path / "x.log"
        path.write_text("")
        with pytest.raises(ValueError, match="log"):
            gi.read_geometry_file(path)
