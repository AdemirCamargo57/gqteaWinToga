"""runatomInput: fixed-column .adat files exactly as runatom (rwsubs.f) reads them."""

import random
import subprocess
import sys
from pathlib import Path

import pytest

import runatomInput as ri

SRC = Path(__file__).resolve().parents[1]


def test_module_does_not_import_toga():
    code = "import sys, runatomInput; assert 'toga' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], cwd=SRC, check=True)


class TestReadField:
    def test_blanks_inside_an_integer_field_are_ignored(self):
        # gfortran BLANK='NULL': "1 0 0" in an i5 field is 100 (verified with gfortran 11.5)
        assert ri.read_record("1 0 0 0 3", "6i5", 5) == [100, 3, 0, 0, 0]

    def test_f_field_without_decimal_point_is_scaled(self):
        assert ri.read_record("   2.5        2         1.8", "3f10.5") == [2.5, 2e-5, 1.8]

    def test_explicit_decimal_point_wins(self):
        assert ri.read_field("f", 10, 5, "      2.  ") == 2.0
        assert ri.read_field("f", 5, 0, "  22.") == 22.0
        assert ri.read_field("f", 14, 6, "      -.55710 ") == -0.5571

    def test_fortran_exponents(self):
        assert ri.read_field("e", 10, 1, "   1.0d-10") == 1e-10
        assert ri.read_field("f", 10, 5, "1.23-100  ") == 1.23e-100

    def test_empty_and_short_records_read_as_zero(self):
        assert ri.read_record("    7", "2i5,f10.5") == [7, 0, 0.0]

    def test_sign_without_exponent_is_a_read_error(self):
        # Work/023-V-gpw-nsp-campos/V_ps.adat line 2: runatom stops here
        with pytest.raises(ValueError, match="floating point"):
            ri.read_record("  2.21   -6.0       2.0       80", "3f10.5,i5")

    def test_text_in_an_integer_field_is_a_read_error(self):
        with pytest.raises(ValueError):
            ri.read_field("i", 5, 0, "  1.0")


class TestWriteField:
    def test_every_float_is_written_with_a_decimal_point(self):
        for value in (2.0, 0.0, -6.0, 82.0, 1e-10, 150.0, -9000.0):
            assert "." in ri.write_field("f", 10, 5, value)

    def test_written_floats_read_back(self):
        rng = random.Random(7)
        for _ in range(2000):
            typed = round(rng.uniform(-1000, 1000), 4)                 # what people type: exact
            assert ri.read_field("f", 14, 6, ri.write_field("f", 14, 6, typed)) == typed
            tiny = rng.choice([1, -1]) * 10 ** rng.uniform(-9, -4)    # needs an exponent: 7 digits kept
            assert ri.read_field("f", 14, 6, ri.write_field("f", 14, 6, tiny)) == pytest.approx(tiny, rel=1e-6)

    def test_exact_value_may_fill_the_field(self):
        assert ri.write_field("f", 5, 0, 100.0) == "100.0"

    def test_integers_are_right_justified_and_checked(self):
        assert ri.write_field("i", 5, 0, 3) == "    3"
        with pytest.raises(ValueError):
            ri.write_field("i", 4, 0, 12345)

    def test_text_longer_than_the_field_is_refused(self):
        with pytest.raises(ValueError):
            ri.write_field("a", 20, 0, "a" * 21)


class TestFormatLine:
    def test_label_never_inside_the_columns_read(self):
        (line,) = ri.format_line("6i5", [1, 0, 0, 0, 3], "ifae,ifpsp,ifprt,ifplw,ilogd")
        assert line.index("ifae") >= ri.LABEL_COLUMN
        assert ri.read_record(line, "6i5", 5) == [1, 0, 0, 0, 3]

    def test_label_after_the_full_format_even_for_a_partial_record(self):
        (line,) = ri.format_line("5f10.5", [1.0, 2.0], "rinner")
        assert line.index("rinner") >= 52
        assert ri.read_record(line, "5f10.5") == [1.0, 2.0, 0.0, 0.0, 0.0]

    def test_format_reversion_continues_on_a_new_line(self):
        values = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]  # nang = 4: 7 rinner values, format 5f10.5
        lines = ri.format_line("5f10.5", values, "rinner")
        assert len(lines) == 2
        assert ri.read_record(lines[0], "5f10.5") + ri.read_record(lines[1], "5f10.5", 2) == values


from runatom_fixtures import LIBRARY_DIR

LIBRARY_INPUTS = sorted(p for p in LIBRARY_DIR.rglob("*.adat") if "History" not in p.parts)
UNREADABLE = {"V_ps.adat": 2}           # upstream defect, see PROVENANCE.md
CASE_DIRS = sorted(p.parent for p in LIBRARY_DIR.rglob("Makefile") if p.parent.parent != LIBRARY_DIR)


def _ti_case():
    return ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-test")


class TestAdatFiles:
    def test_library_is_large_enough_to_mean_something(self):
        assert len(LIBRARY_INPUTS) > 400

    @pytest.mark.parametrize("path", LIBRARY_INPUTS, ids=lambda p: f"{p.parent.name}/{p.name}")
    def test_every_library_input_round_trips(self, path):
        text = path.read_text(encoding="latin-1")
        if path.name in UNREADABLE:
            with pytest.raises(ri.AdatError) as info:
                ri.parse_adat(text)
            assert info.value.line_number == UNREADABLE[path.name]
            return
        parsed = ri.parse_adat(text)
        assert ri.parse_adat(ri.format_adat(parsed)) == parsed

    def test_written_file_uses_lf_and_ends_with_newline(self):
        text = ri.format_adat(_ti_case().ps)
        assert "\r" not in text and text.endswith("\n")

    def test_npf_line_only_when_a_projector_is_smooth(self):
        ps = _ti_case().ps
        assert "npf,ptryc" in ri.format_adat(ps)
        plain = ri.replace(ps, betas=[ri.replace(b, iptype=1) for b in ps.betas])
        text = ri.format_adat(plain)
        assert "npf,ptryc" not in text
        assert ri.parse_adat(text).lloc == ps.lloc   # the local line was not shifted

    def test_per_l_rinner_line_when_rinner_is_zero(self):
        case = ri.read_case_folder(LIBRARY_DIR / "026-Fe" / "026-Fe-gpw-n-campos")
        assert case.ps.rinner1 == 0.0 and len(case.ps.rinner) == 2 * case.ps.nang - 1


class TestCaseFolders:
    @pytest.mark.parametrize("folder", CASE_DIRS, ids=lambda p: p.name)
    def test_every_library_case_reads(self, folder):
        if (folder / "V_ps.adat").exists():
            with pytest.raises(ri.AdatError):
                ri.read_case_folder(folder)
            return
        case = ri.read_case_folder(folder)
        assert case.ae.controls.ifae == 1 and case.ps.controls.ifpsp == 2

    def test_tests_named_cf1_cf2_are_found(self):
        case = ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-tests")
        assert [t.name for t in case.tests] == ["s2d2", "s1d2"]

    def test_tests_become_occupation_changes(self):
        case = _ti_case()
        assert case.symbol == "ti" and case.generation_name == "s2d1"
        (test,) = case.tests
        assert test.name == "s2d2" and test.xion == 0.0 and test.occupations == {320: 2.0}
        derived = ri.test_ae_input(case, test)
        assert {s.nnlz: s.wwnl for s in derived.shells}[320] == 2.0

    def test_tests_with_their_own_shell_list_keep_their_input(self):
        case = ri.read_case_folder(LIBRARY_DIR / "011-Na" / "011-Na-ca-n-campos")
        own = [t for t in case.tests if t.ae_input is not None]
        assert own, "expected at least one test with its own shell list"
        assert all(t.occupations == {} for t in own)
        assert ri.test_ae_input(case, own[0]) is own[0].ae_input

    def test_written_case_reads_back(self, tmp_path):
        case = _ti_case()
        written = ri.write_case(case, tmp_path)
        assert sorted(p.name for p in written) == ["ti_ae_s2d1.adat", "ti_ae_s2d2.adat",
                                                   "ti_ps.adat", "ti_test_s2d2.adat"]
        (tmp_path / "Makefile").write_text(ri.makefile_text(case), newline="\n")
        again = ri.read_case_folder(tmp_path)
        assert again.ae == case.ae and again.ps == case.ps
        assert [(t.name, t.xion, t.occupations) for t in again.tests] == [("s2d2", 0.0, {320: 2.0})]

    def test_makefile_lists_every_stage(self):
        text = ri.makefile_text(_ti_case())
        assert "ATOM= ti" in text and "CFG= s2d1" in text and "CFT1= s2d2" in text
        assert "${PROG} ${TST1}.adat" in text
