"""Suggested runatom inputs per element: library cases, ranking, starter and templates."""

import subprocess
import sys
from pathlib import Path

import pytest

import runatomDefaults as rd
import runatomInput as ri
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
