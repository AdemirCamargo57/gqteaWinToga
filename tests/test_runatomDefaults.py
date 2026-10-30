"""Suggested runatom inputs per element: library cases, ranking, starter and templates."""

import subprocess
import sys
from pathlib import Path

import pytest

import runatomDefaults as rd
import runatomInput as ri
from help import AtomicData
from runatom_fixtures import LIBRARY_DIR

ROOTS = [LIBRARY_DIR]
LIBRARY_SYMBOLS = {
    "H", "Li", "B", "C", "N", "O", "F", "Na", "Al", "Si", "P", "S", "Cl", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Fe", "Ni", "Cu", "Zn", "Ga", "As", "Se", "Sr", "Y", "Zr", "Nb", "Ru", "Rh", "Ag", "Cd", "Ba", "La",
    "Ta", "W", "Os", "Ir", "Pt", "Au", "Tl", "Pb", "Bi",
}


def test_module_does_not_import_toga():
    code = "import sys, runatomDefaults; assert 'toga' not in sys.modules, 'toga imported'"
    subprocess.run([sys.executable, "-c", code], check=True, cwd=Path(__file__).resolve().parents[1])


def test_xc_label_never_prints_minus_zero():
    assert ri.xc_label(-0.0) == ri.xc_label(0.0) == "Ceperley-Alder (LDA)"
    assert ri.xc_label(5.0) == "PBE (GGA)"
    assert ri.xc_label(0.7) == "exfact 0.7"


class TestLibraryIndex:
    def test_the_library_covers_45_elements(self):
        assert set(rd.library_elements(ROOTS)) == LIBRARY_SYMBOLS

    def test_the_example_and_history_folders_are_not_cases(self):
        folders = [f for found in rd.library_elements(ROOTS).values() for f in found]
        assert not any("000-Ex" in f.parts or "History" in f.parts for f in folders)

    def test_a_missing_root_is_skipped(self, tmp_path):
        assert rd.library_elements([tmp_path / "missing"]) == {}

    def test_titanium_has_two_described_variants(self):
        variants, notes = rd.element_variants("Ti", ROOTS)
        assert notes == []
        assert [v.folder.name for v in variants] == ["022-Ti-ca-sp-vgrp", "022-Ti-ca-sp-vgrp_alt"]
        first = variants[0]
        assert first.z == 22 and first.irel == 2 and first.valence[:2] == (300, 310)
        assert "Ceperley-Alder (LDA)" in first.label and "scalar-relativistic" in first.label
        assert "semicore 3s 3p" in first.label and "+1 ion" in first.label
        assert "+1 ion" not in variants[1].label            # the _alt case is the neutral atom

    def test_the_unreadable_vanadium_case_is_reported(self):
        variants, notes = rd.element_variants("V", ROOTS)
        assert variants == []
        assert len(notes) == 1 and "V_ps.adat" in notes[0] and "line 2" in notes[0]


class TestRanking:
    @pytest.mark.parametrize("symbol, best", [
        ("H", "001-H-gpbe--bm"),           # PBE first
        ("Cu", "029-Cu-gpbe-n-bm"),
        ("Ca", "020-Ca-gpw-nsp-campos"),   # no PBE: PW91 before LDA
        ("Ti", "022-Ti-ca-sp-vgrp"),       # LDA only: deterministic tie-break
    ])
    def test_pbe_first_then_pw91_then_lda(self, symbol, best):
        variants, _ = rd.element_variants(symbol, ROOTS)
        assert rd.rank_variants(variants)[0].folder.name == best

    def test_ranking_keeps_every_variant(self):
        variants, _ = rd.element_variants("Cu", ROOTS)
        assert sorted(v.folder for v in rd.rank_variants(variants)) == sorted(v.folder for v in variants)


ALL_Z = range(1, rd.MAX_DEFAULTS_Z + 1)


def _config_text(z):
    return " ".join(f"{n}{'spdf'[l]}{int(e)}" for n, l, e in rd.ground_state_configuration(z))


class TestConfiguration:
    @pytest.mark.parametrize("z", ALL_Z, ids=rd.element_symbol)
    def test_every_configuration_adds_up_to_z(self, z):
        assert sum(e for _, _, e in rd.ground_state_configuration(z)) == z

    @pytest.mark.parametrize("symbol, outer", [
        ("Fe", "4s2 3d6"), ("Cr", "4s1 3d5"), ("Cu", "4s1 3d10"), ("Nb", "5s1 4d4"), ("Mo", "5s1 4d5"),
        ("Ru", "5s1 4d7"), ("Rh", "5s1 4d8"), ("Ag", "5s1 4d10"), ("La", "6s2 5d1"), ("Ce", "6s2 4f1 5d1"),
        ("Gd", "6s2 4f7 5d1"), ("Pt", "6s1 4f14 5d9"), ("Au", "6s1 4f14 5d10"), ("Ac", "7s2 6d1"),
        ("Th", "7s2 6d2"), ("Pa", "7s2 5f2 6d1"), ("U", "7s2 5f3 6d1"),
    ])
    def test_nist_ground_configurations(self, symbol, outer):
        assert _config_text(AtomicData.atomic_numbers[symbol]).endswith(outer)

    def test_palladium_has_no_5s_electron(self):
        assert _config_text(46).endswith("4p6 4d10")

    def test_only_listed_shells_change(self):
        madelung_cr = "1s2 2s2 2p6 3s2 3p6 4s1 3d5"
        assert _config_text(24) == madelung_cr


class TestCoreValence:
    @pytest.mark.parametrize("symbol, valence", [
        ("H", [100]), ("Li", [200]), ("Be", [200]), ("C", [200, 210]),
        ("Na", [200, 210, 300]), ("Mg", [200, 210, 300]), ("K", [300, 310, 400]),
        ("Ti", [300, 310, 400, 320]), ("Cr", [300, 310, 400, 320]), ("Mn", [400, 320]),
        ("Cu", [400, 320]), ("Zn", [400, 320]), ("Ga", [400, 410]), ("Mo", [400, 410, 500, 420]),
        ("Pd", [420]), ("La", [500, 510, 600, 520]), ("Ce", [600, 430, 520]), ("Hf", [500, 510, 600, 520]),
        ("Pb", [600, 610]), ("Ac", [600, 610, 700, 620]), ("Th", [700, 620]), ("U", [700, 530, 620]),
    ])
    def test_valence_shells(self, symbol, valence):
        assert [s.nnlz for s in _valence(symbol)] == valence

    @pytest.mark.parametrize("z", ALL_Z, ids=rd.element_symbol)
    def test_every_starter_passes_runatoms_all_electron_checks(self, z):
        symbol = rd.element_symbol(z)
        ae = rd.starter_ae_input(symbol)
        fields = rd.starter_fields(symbol)
        ncores, nvales, nang = (int(v) for v in fields["ps_partition"].split())
        assert ri.check_ae(ae) == []
        assert ncores + nvales == len(ae.shells) and 1 <= nvales <= ri.MAX_VALENCE and nang <= ri.MAX_NANG
        assert ae.exfact == rd.PBE and ae.irel == 2

    @pytest.mark.parametrize("z", ALL_Z, ids=rd.element_symbol)
    def test_every_guess_is_runatoms_own_converged_eigenvalue(self, z):
        # Slater's rules miss open d and f shells by up to 100x; the table never needs them up to U
        ae = rd.starter_ae_input(rd.element_symbol(z))
        assert {s.nnlz: s.ee for s in ae.shells} == rd.STARTER_EIGENVALUES[z]

    def test_three_lanthanides_start_as_ions(self):
        # runatom never converges the 4f state of the neutral Dy, Ho and Er atoms
        assert {rd.element_symbol(z) for z in rd.STARTER_XION} == {"Dy", "Ho", "Er"}
        ae = rd.starter_ae_input("Dy")
        assert ae.xion == 1.0 and {s.nnlz: s.wwnl for s in ae.shells}[600] == 1.0
        assert rd.generation_tag("Dy") == "s1f10"


class TestStarterFields:
    def test_titanium(self):
        fields = rd.starter_fields("Ti")
        assert fields["symbol"] == "Ti" and fields["generation"] == "s2d2" and fields["title"] == "Ti"
        assert fields["ae_atom"] == "22  0  5" and fields["ps_partition"] == "3  4  3"
        assert fields["ae_mesh"] == f"{rd.STARTER_RMAX:g}  0  40  2"

    def test_the_pseudization_choices_are_left_blank(self):
        fields = rd.starter_fields("Mg")
        for key in ("ps_bessel", "ps_keyps", "ps_rinner", "ps_rcloc", "ps_rc", "ps_betas", "ps_npf",
                    "ps_local", "ps_qfunc", "ps_fixes", "tests"):
            assert fields[key] == "", key
        assert set(fields) == set(ri.FORM_FIELDS)

    def test_blank_fields_are_named_when_the_form_is_read(self):
        with pytest.raises(ValueError, match="keyps  ifpcor  rinner  rpcor: expected 4 values, got 0"):
            ri.case_from_fields(rd.starter_fields("Mg"))

    @pytest.mark.parametrize("symbol", ["Np", "Xx"])
    def test_beyond_uranium_or_unknown_is_refused(self, symbol):
        with pytest.raises(ValueError):
            rd.starter_fields(symbol)


class TestTemplates:
    @pytest.mark.parametrize("symbol, donor", [
        ("Mg", "Ca"), ("Be", "Ca"), ("Ge", "Si"), ("Sn", "Pb"), ("Mo", "W"), ("V", "Nb"),
        ("Ce", "La"), ("U", "La"), ("I", "Cl"), ("Rb", "K"), ("Co", "Rh"), ("Pd", "Pt"), ("Hf", "Zr"),
    ])
    def test_nearest_period_of_the_same_group_heavier_on_a_tie(self, symbol, donor):
        assert rd.template_donors(symbol, ROOTS)[0] == donor

    @pytest.mark.parametrize("symbol", ["Mn", "Tc", "Re", "He", "Ne", "Rn"])
    def test_groups_7_and_18_have_no_donor(self, symbol):
        assert rd.template_donors(symbol, ROOTS) == []

    def _template(self, symbol):
        options, _ = rd.suggestions(symbol, ROOTS)
        assert options[0].kind == "template"
        return options[0]

    def test_magnesium_takes_calciums_pseudopotential_tab(self):
        template = self._template("Mg")
        donor_folder = template.folder
        donor = ri.read_case_folder(donor_folder)
        assert donor_folder.parent.name == "020-Ca"
        donor_fields = ri.ps_fields(donor.ps)
        for key in ("ps_rc", "ps_rcloc", "ps_rinner", "ps_keyps", "ps_bessel", "ps_logd"):
            assert template.fields[key] == donor_fields[key], key
        case = ri.case_from_fields(template.fields)
        assert case.ae.z == 12 and case.ae.exfact == rd.PBE and case.symbol == "Mg"
        assert case.ps.nvales == donor.ps.nvales - 1 and case.ps.nang == donor.ps.nang     # Ca's empty 3d dropped

    @pytest.mark.parametrize("symbol, valence", [
        # Ca: 3s2 3p6 4s2 3d0, one period up; Mg has no d electron, so the empty d channel is dropped
        ("Mg", [(200, 2), (210, 6), (300, 2)]),
        # Pb (5d10 6s2 6p2) one period down: Sn's filled 4d joins the valence, as Pb's 5d does
        ("Sn", [(420, 10), (500, 2), (510, 2)]),
        # Nb: 4s2 4p6 4d4 5s1 5p0, one period down: V's own 3d3 4s2 occupations, the empty 4p dropped
        ("V", [(300, 2), (310, 6), (320, 3), (400, 2)]),
    ])
    def test_the_valence_list_follows_the_donor_with_the_targets_occupations(self, symbol, valence):
        case = ri.case_from_fields(self._template(symbol).fields)
        assert [(s.nnlz, s.wwnl) for s in case.ae.shells[case.ps.ncores:]] == valence
        assert sum(max(s.wwnl, 0) for s in case.ae.shells) == case.ae.z      # the target's own electrons
        assert ri.check_ae(case.ae) == []

    def test_projectors_follow_their_valence_shells(self):
        # keyee is a position in the valence list (its eigenvalue is the reference energy), whatever its l;
        # a projector on a dropped shell keeps that shell's energy from the donor input as eeread
        template = self._template("V")
        donor = ri.read_case_folder(template.folder)            # Nb: 4s 4p 4d 5s 5p(empty)
        dropped_ee = donor.ae.shells[donor.ps.ncores + 4].ee
        expected = [(b.l, 0, dropped_ee, b.iptype) if b.keyee == 5 else (b.l, b.keyee, b.eeread, b.iptype)
                    for b in donor.ps.betas]
        case = ri.case_from_fields(template.fields)
        assert [(b.l, b.keyee, b.eeread, b.iptype) for b in case.ps.betas] == expected
        assert any("5p" in w and "eeread" in w for w in template.warnings)

    def test_the_template_all_electron_input_has_the_starters_shells(self):
        for symbol in ("Be", "Mg", "V", "Sn", "Mo", "Pd", "Th"):
            case = ri.case_from_fields(self._template(symbol).fields)
            starter = rd.starter_ae_input(symbol)
            assert sorted((s.nnlz, s.wwnl, s.ee) for s in case.ae.shells) == \
                sorted((s.nnlz, s.wwnl, s.ee) for s in starter.shells), symbol

    def test_switching_to_ifqopt_3_raises_nqf_to_psqf2s_minimum(self):
        # Bi's LDA case uses ifqopt 2 with nqf 5; ifqopt 3 needs at least 6
        template = self._template("Sb")
        ifqopt, nqf, _ = template.fields["ps_qfunc"].split()
        assert ifqopt == "3" and int(nqf) >= 6
        assert any("nqf" in w for w in template.warnings)

    def test_the_warning_names_the_donor_and_says_the_radii_are_unchanged(self):
        warnings = " ".join(self._template("Mg").warnings)
        assert "Ca" in warnings and "020-Ca-" in warnings and "unchanged" in warnings
        assert "not a tested pseudopotential" in warnings

    def test_gga_templates_use_ifqopt_3(self):
        for symbol in ("Mg", "Ge", "Sn", "Mo", "Ce"):
            fields = self._template(symbol).fields
            if fields["ps_qfunc"]:
                assert fields["ps_qfunc"].split()[0] == "3", symbol

    def test_a_shell_the_donor_cannot_provide_gives_a_radii_only_template(self):
        # Ce's occupied 4f valence shell has no counterpart among La's valence shells (5s 5p 5d 6s 6p)
        template = self._template("Ce")
        assert template.fields["ps_rc"] == ri.ps_fields(ri.read_case_folder(template.folder).ps)["ps_rc"]
        assert template.fields["ps_betas"] == template.fields["ps_local"] == ""
        assert any("4f" in w and "projectors" in w for w in template.warnings)

    @pytest.mark.parametrize("symbol", ["Be", "Mg", "Rb", "Ge", "Sn", "Mo", "Hf", "V", "Co", "Pd",
                                        "In", "Sb", "Te", "I", "Cs", "Hg", "Po", "Ce", "Th", "U"])
    def test_check_problems_are_never_silent(self, symbol):
        option = rd.suggestions(symbol, ROOTS)[0][0]
        try:
            problems = ri.check_case(ri.case_from_fields(option.fields))
        except ValueError as exc:
            problems = [str(exc)]
        for problem in problems:
            assert any(problem in w for w in option.warnings), (symbol, problem)


class TestSuggestions:
    def test_a_library_element_gets_its_ranked_cases_unmodified(self):
        options, _ = rd.suggestions("Cu", ROOTS)
        assert [o.kind for o in options] == ["library"] * 4
        assert options[0].folder.name == "029-Cu-gpbe-n-bm" and options[0].warnings == ()
        assert options[0].fields == ri.case_to_fields(ri.read_case_folder(options[0].folder))

    def test_a_library_element_without_pbe_is_warned(self):
        options, _ = rd.suggestions("Ti", ROOTS)
        assert options[0].kind == "library" and "No PBE case for Ti" in options[0].warnings[0]

    def test_an_element_without_a_case_gets_a_template_then_the_starter(self):
        options, _ = rd.suggestions("Mg", ROOTS)
        assert [o.kind for o in options] == ["template", "starter"]
        assert options[0].label.startswith("Template from 020-Ca-")
        assert options[1].label == rd.STARTER_LABEL and options[1].fields == rd.starter_fields("Mg")

    def test_vanadium_reports_its_unreadable_case_and_borrows_from_niobium(self):
        options, notes = rd.suggestions("V", ROOTS)
        assert options[0].kind == "template" and options[0].folder.parent.name == "041-Nb"
        assert any("V_ps.adat" in w for w in options[0].warnings) and notes

    def test_group_7_gets_the_starter_with_a_warning(self):
        options, _ = rd.suggestions("Mn", ROOTS)
        assert [o.kind for o in options] == ["starter"]
        assert "group 7" in options[0].warnings[0]

    def test_the_ionic_starters_say_so(self):
        options, _ = rd.suggestions("Ho", ROOTS)
        assert any("+1 ion" in w for w in options[0].warnings)

    def test_without_a_library_every_element_gets_the_starter(self, tmp_path):
        options, _ = rd.suggestions("Ti", [tmp_path])
        assert [o.kind for o in options] == ["starter"]
        assert "library" in options[0].warnings[0].lower()

    def test_beyond_uranium_is_refused(self):
        with pytest.raises(ValueError, match="92"):
            rd.suggestions("Np", ROOTS)


def _valence(symbol):
    ae = rd.starter_ae_input(symbol)
    ncores = int(rd.starter_fields(symbol)["ps_partition"].split()[0])
    return ae.shells[ncores:]
