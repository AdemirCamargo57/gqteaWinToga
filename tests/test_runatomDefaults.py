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


def _valence(symbol):
    ae = rd.starter_ae_input(symbol)
    ncores = int(rd.starter_fields(symbol)["ps_partition"].split()[0])
    return ae.shells[ncores:]
