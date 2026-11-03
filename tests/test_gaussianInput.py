"""Tests for the Gaussian 16 input builder core (gaussianInput, no Toga)."""
import subprocess
import sys
from pathlib import Path

import pytest

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
