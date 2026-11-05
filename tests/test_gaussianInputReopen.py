"""Open input...: Gaussian inputs read back into the panels (gaussianInput, no Toga)."""
from dataclasses import replace
from pathlib import Path

import pytest

import fukuiInput as fi
import gaussianInput as gi
import gaussianKeywords as gk

DATA = Path(__file__).resolve().parent / "data"
FIXTURES = sorted(p for p in DATA.rglob("*") if p.suffix.lower() in gi.GAUSSIAN_INPUT_SUFFIXES)
# the layout the builder writes (and the reader returns), so a written file reads back to the same text
WATER = "\n".join(fi.format_geometry(fi.parse_geometry(
    "O 0.000000 0.000000 0.117300\nH 0.000000 0.757200 -0.469200\nH 0.000000 -0.757200 -0.469200")))


def step(**changes):
    return replace(gi.GaussianStepSpec(), **changes)


def spec_of(tmp_path, *steps, geometry=WATER):
    return gi.GaussianJobSpec(stem="water", output_dir=str(tmp_path), geometry=geometry, steps=list(steps))


ROUND_TRIPS = {
    "energy": [step(multiplicity=1)],
    "opt freq": [step(job_type=gk.JOB_OPT_FREQ, multiplicity=1, freq_raman=True, freq_anharmonic=True)],
    "ts": [step(job_type=gk.JOB_OPT, multiplicity=1, opt_target=gk.OPT_TS, force_constants=gk.FC_CALCFC,
                opt_tight=True, opt_max_cycles=200, additional_keywords="Opt=NoEigenTest 6D")],
    "irc": [step(job_type=gk.JOB_IRC, multiplicity=1, irc_direction="Reverse only", irc_max_points=30,
                 irc_step_size=8)],
    "ts freq then irc rcfc": [
        step(job_type=gk.JOB_OPT_FREQ, multiplicity=1, opt_target=gk.OPT_TS, force_constants=gk.FC_CALCFC),
        replace(gi.linked_step(step(job_type=gk.JOB_OPT_FREQ, multiplicity=1)), job_type=gk.JOB_IRC,
                irc_force_constants=gk.FC_RCFC, irc_max_points=5)],
    "scan": [step(job_type=gk.JOB_SCAN, multiplicity=1, modredundant="B 1 2 S 5 0.05\nA 2 1 3 F",
                  connectivity=True)],
    "stability": [step(job_type=gk.JOB_STABILITY, multiplicity=1, stable_opt=True)],
    "nmr": [step(job_type=gk.JOB_NMR, multiplicity=1, nmr_method="CSGT", nmr_spinspin=True)],
    "td": [step(multiplicity=1, state=gk.STATE_TD, n_states=10, root=2, state_spin="Singlets",
                print_all_coefficients=True, dispersion="GD3BJ")],
    "cis": [step(multiplicity=1, state=gk.STATE_CIS, n_states=4, family=gk.FAMILY_HF, method_name=gk.NO_CHOICE,
                 basis="cc-pVDZ")],
    "zindo": [step(multiplicity=1, state=gk.STATE_ZINDO, n_states=8, family=gk.FAMILY_HF,
                   method_name=gk.NO_CHOICE)],
    "semi": [step(multiplicity=1, family=gk.FAMILY_SEMI, method_name="PM6")],
    "compound": [step(multiplicity=1, family=gk.FAMILY_COMPOUND, method_name="G4")],
    "ccsd(t)": [step(multiplicity=1, family=gk.FAMILY_CCSDT, method_name=gk.NO_CHOICE, basis="aug-cc-pVTZ")],
    "open shell": [step(charge=1, multiplicity=2, spin="Unrestricted (U)", method_name="M06-2X",
                        basis="6-311G", basis_diffuse="++", basis_polarization="(2df,2pd)")],
    "other names": [step(multiplicity=1, method_name=gk.OTHER, method_other="revPBE", basis=gk.OTHER,
                         basis_other="def2TZVPD", solvation="CPCM", solvent=gk.OTHER, solvent_other="Heptane")],
    "general": [step(multiplicity=1, print_level="#N", scf_tight=True, scf_qc=True, nosymm=True, polar=True,
                     pop_level="Full", pop_mk=True, pop_nbo=True, solvation="SMD", solvent="Acetonitrile",
                     title="my water", memory="8GB", nproc=16, chk="w1.chk", old_chk="start.chk",
                     extra_link0="%KJob L502")],
    "gen basis": [step(multiplicity=1, basis="Gen", additional_input="O H 0\n6-31G(d)\n****")],
    "no checkpoint": [step(multiplicity=1, write_chk=False)],
    "linked": [step(job_type=gk.JOB_OPT, multiplicity=1),
               replace(gi.linked_step(step(job_type=gk.JOB_OPT, multiplicity=1)), state=gk.STATE_TD, n_states=6,
                       guess_mix=True),
               step(geometry_source=gk.GEOM_CHECK, charge=-1, multiplicity=2, guess="Read checkpoint"),
               step(multiplicity=3, title="triplet")],
}


@pytest.mark.parametrize("name", sorted(ROUND_TRIPS))
def test_a_file_the_builder_wrote_reads_back_to_its_spec(tmp_path, name):
    builder = gi.GaussianInputBuilder()
    spec = spec_of(tmp_path, *ROUND_TRIPS[name])
    text = builder.input_text(builder.prepare(spec))
    back, _ = gi.read_input_text(text, spec.stem, spec.output_dir)
    assert back == spec


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_every_project_input_reopens_prepares_and_keeps_its_keywords(tmp_path, path):
    spec, _ = gi.read_input_file(path)
    blocks = gi._split_blocks(path.read_text(errors="replace"))      # one per linked step
    assert len(blocks) == len(spec.steps)
    for block, step_spec in zip(blocks, spec.steps):
        route = " ".join(l.strip() for l in block if l.strip().startswith("#")).split(None, 1)[1]
        names = {gi.keyword_text_name(t).upper() for t in fi.route_tokens(route) if "/" not in t}
        emitted = {gi.keyword_text_name(t).upper() for t in fi.route_tokens(" ".join(gi.route_lines(step_spec)))}
        assert names <= emitted | {gi.keyword_text_name(gi.method_keyword(step_spec)).upper()}
    builder = gi.GaussianInputBuilder()
    job = builder.prepare(replace(spec, output_dir=str(tmp_path)))
    again, _ = gi.read_input_text(builder.input_text(job), spec.stem, str(tmp_path))
    assert again == replace(spec, output_dir=str(tmp_path))


class TestReading:
    def test_star_basis_sets_and_lower_case(self):
        spec, _ = gi.read_input_text("#p pbe1pbe/6-31g** opt\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n")
        s = spec.steps[0]
        assert (s.method_name, s.basis, s.basis_polarization, s.job_type, s.print_level) == \
            ("PBE0 (PBE1PBE)", "6-31G", "(d,p)", gk.JOB_OPT, "#P")

    def test_a_route_without_a_method_is_hf_sto3g_with_a_note(self):
        spec, notes = gi.read_input_text("# opt\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n")
        s = spec.steps[0]
        assert (s.family, s.basis) == (gk.FAMILY_HF, "STO-3G") and any("HF/STO-3G" in n for n in notes)

    def test_absolute_checkpoint_path_keeps_the_file_name(self):
        spec, notes = gi.read_input_text("%chk=C:\\jobs\\w.chk\n# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n",
                                         stem="other")
        assert spec.steps[0].chk == "w.chk" and any("C:\\jobs\\w.chk" in n for n in notes)

    def test_forbidden_title_characters_are_replaced(self):
        spec, notes = gi.read_input_text("# hf/sto-3g\n\nB3LYP/6-31G* opted\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n")
        assert spec.steps[0].title == "B3LYP/6 31G* opted" and any("title" in n for n in notes)

    def test_connectivity_section_is_consumed_with_a_note(self):
        text = "# hf/sto-3g geom=connectivity\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n 1 2 1.0\n 2\n\nextra\n\n"
        spec, notes = gi.read_input_text(text)
        assert spec.steps[0].connectivity and spec.steps[0].additional_input == "extra"
        assert any("connectivity" in n for n in notes)

    def test_fukui_file_reopens_as_three_steps_with_a_note(self, tmp_path):
        s = gi.fukui_defaults(step(multiplicity=1, fukui_plus_multiplicity=2, fukui_minus_multiplicity=2))
        builder = gi.GaussianInputBuilder()
        spec = spec_of(tmp_path, s)
        spec.stem = "water"
        text = builder.input_text(builder.prepare(spec))
        back, notes = gi.read_input_text(text, "water_fukui")
        assert len(back.steps) == 3 and any("Fukui" in n for n in notes)
        assert back.steps[1].geometry_source == gk.GEOM_CHECK and back.steps[1].old_chk == "water_N.chk"

    @pytest.mark.parametrize("text,fragment", [
        ("", "empty"),
        ("%mem=1GB\n\nt\n", "no route section"),
        ("# hf/sto-3g\n\nt\n\n0 1 0 1 0 1\nH 0 0 0\nH 0 0 0.74\n\n", "per-fragment charges"),
        ("# hf/sto-3g\n\nt\n\n0 1\nO\nH 1 0.96\n\n", "Only Cartesian coordinates"),
        ("# hf/sto-3g\n\nt\n\n", "ends before the molecule specification"),
        ("# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n--Link1--\n# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.80\n\n",
         "Step 2 has its own geometry"),
        ("# hf/sto-3g geom=check\n\nt\n\n", "Geom=Check needs a title"),
        ("mem=1GB\n# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\n\n", "is not a Link 0"),
    ])
    def test_refusals(self, text, fragment):
        with pytest.raises(ValueError, match=fragment.replace("(", "\\(").replace(")", "\\)")):
            gi.read_input_text(text)

    def test_read_input_file_names_the_file(self, tmp_path):
        path = tmp_path / "bad.gjf"
        path.write_text("nothing here\n")
        with pytest.raises(ValueError, match="bad.gjf"):
            gi.read_input_file(path)
        with pytest.raises(ValueError, match="not a Gaussian input"):
            gi.read_input_file(tmp_path / "x.xyz")

    def test_read_input_file_takes_stem_and_folder_from_the_path(self, tmp_path):
        path = tmp_path / "my mol.com"
        path.write_text("# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n")
        spec, _ = gi.read_input_file(path)
        assert (spec.stem, spec.output_dir) == ("my_mol", str(tmp_path.resolve()))


class TestReviewFixes:
    H2 = "0 1\nH 0 0 0\nH 0 0 0.74\n\n"

    def test_comment_line_between_link0_and_route(self):
        spec, _ = gi.read_input_text("%chk=a.chk\n! a note\n# hf/sto-3g\n\nt\n\n" + self.H2)
        assert spec.steps[0].chk == "a.chk" and spec.steps[0].method_name == gk.NO_CHOICE

    def test_comment_line_inside_the_atom_list(self):
        spec, _ = gi.read_input_text("# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\n! second atom\nH 0 0 0.74\n\n")
        assert len(spec.geometry.splitlines()) == 2 and spec.steps[0].additional_input == ""

    def test_trailing_comment_is_still_stripped(self):
        spec, _ = gi.read_input_text("# hf/sto-3g ! method\n\nt\n\n" + self.H2)
        assert spec.steps[0].additional_keywords == ""

    def test_a_second_job_keyword_is_kept(self):
        spec, notes = gi.read_input_text("# hf/sto-3g nmr freq\n\nt\n\n" + self.H2)
        s = spec.steps[0]
        assert s.job_type == gk.JOB_FREQ and "nmr" in s.additional_keywords.lower()
        assert {"NMR", "FREQ"} <= {gi.keyword_text_name(t).upper() for t in fi.route_tokens(" ".join(gi.route_lines(s)))}
        assert any("nmr" in n for n in notes)

    def test_opt_with_irc_keeps_the_irc_options(self):
        spec, _ = gi.read_input_text("# hf/sto-3g opt irc=calcfc\n\nt\n\n" + self.H2)
        s = spec.steps[0]
        assert s.job_type == gk.JOB_OPT and "irc=calcfc" in s.additional_keywords
        route = " ".join(gi.route_lines(s)).upper()
        assert "IRC" in route and "CALCFC" in route

    def test_unmodelled_cis_options_are_dropped_with_a_plain_note(self):
        spec, notes = gi.read_input_text("# cis(Direct,NStates=4)/sto-3g\n\nt\n\n" + self.H2)
        assert spec.steps[0].n_states == 4
        assert "Step 1: option(s) Direct of CIS have no panel and were dropped; add them back by hand if they matter." in notes

    def test_freeze_flag_coordinate_line_is_refused(self):
        with pytest.raises(ValueError, match="expected 'Element x y z'"):
            gi.read_input_text("# hf/sto-3g\n\nt\n\n0 1\nC -1 0.0 0.0 0.0\nH 0 0 0.9\n\n")

    def test_allcheck_step_has_no_title_or_charge(self):
        text = ("# hf/sto-3g\n\nt\n\n" + self.H2 + "--Link1--\n# hf/sto-3g geom=allcheck guess=read\n\n")
        spec, _ = gi.read_input_text(text)
        s = spec.steps[1]
        assert s.geometry_source == gk.GEOM_ALLCHECK and s.title == "" and s.additional_input == ""
        assert s.guess == "Read checkpoint" and s.multiplicity is None
