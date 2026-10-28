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
