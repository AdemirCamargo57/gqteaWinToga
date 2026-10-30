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

    def test_a_tuned_occupation_test_keeps_its_all_electron_input(self):
        # Ca test 1: the same shells as the generation, but its own damp (and title)
        case = ri.read_case_folder(LIBRARY_DIR / "020-Ca" / "020-Ca-gpw-nsp-campos")
        test = next(t for t in case.tests if t.name == "1")
        assert test.ae_input is not None and test.test_input is not None
        assert test.ae_input.controls.damp != case.ae.controls.damp
        assert ri.test_ae_input(case, test) is test.ae_input
        assert "1" not in ri.unrepresentable_tests(case)
        rows = [line.split()[0] for line in ri.case_to_fields(case)["tests"].splitlines()]
        assert "1" in rows

    def test_an_untuned_test_stays_derived(self):
        (test,) = _ti_case().tests
        assert test.ae_input is None

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
        assert "# It assumes a USPP source tree" in text and "Bin/runatom.x" in text.split("ATOM=")[0]


class TestFields:
    def test_typed_values_are_parsed_not_copied(self):
        assert ri.parse_line_values("1 0 3", "iii", "ifprt ifplw ilogd") == [1, 0, 3]
        assert ri.parse_line_values("1.0d-10, 1e-9  0.5 0", "fffi", "x") == [1e-10, 1e-9, 0.5, 0]
        assert ri.parse_line_values("1.8 1.8 2", "f*", "rc") == [1.8, 1.8, 2.0]

    def test_wrong_count_and_type_name_the_field(self):
        with pytest.raises(ValueError, match="ifprt ifplw ilogd: expected 3 values, got 2"):
            ri.parse_line_values("1 0", "iii", "ifprt ifplw ilogd")
        with pytest.raises(ValueError, match="is not an integer"):
            ri.parse_line_values("1.5 0 3", "iii", "ifprt ifplw ilogd")

    @pytest.mark.parametrize("folder", [p for p in CASE_DIRS if not (p / "V_ps.adat").exists()],
                             ids=lambda p: p.name)
    def test_fields_round_trip_every_library_case(self, folder):
        case = ri.read_case_folder(folder)
        skipped = set(ri.unrepresentable_tests(case))
        kept = [t for t in case.tests if t.name not in skipped]
        expected = ri.replace(
            case, ps=ri.replace(case.ps, controls=ri.replace(case.ps.controls, title=case.ae.controls.title)),
            tests=[ri.replace(t, test_input=None, ae_input=None) for t in kept])
        from_fields = ri.case_from_fields(ri.case_to_fields(case))
        assert from_fields == expected                  # the fields hold name, xion and occupations only
        # the window adds the library's own inputs back to every unchanged row
        assert ri.attach_imported_tests(from_fields, case).tests == kept

    def test_attached_tests_keep_the_library_inputs_only_while_unchanged(self):
        case = ri.read_case_folder(LIBRARY_DIR / "020-Ca" / "020-Ca-gpw-nsp-campos")
        fields = ri.case_to_fields(case)
        rows = fields["tests"].splitlines()
        first = rows[0].split()
        rows[0] = "  ".join([first[0], "0.5"] + first[2:])              # xion changed in row 1 only
        changed = ri.attach_imported_tests(ri.case_from_fields(dict(fields, tests="\n".join(rows))), case)
        assert changed.tests[0].ae_input is None and changed.tests[0].test_input is None
        assert changed.tests[1:] == case.tests[1:]

    def test_no_library_inputs_without_an_import_or_after_a_shell_change(self):
        case = ri.read_case_folder(LIBRARY_DIR / "020-Ca" / "020-Ca-gpw-nsp-campos")
        derived = ri.case_from_fields(ri.case_to_fields(case))
        assert ri.attach_imported_tests(derived, None) == derived
        other_shells = ri.replace(derived, ae=ri.replace(derived.ae, shells=derived.ae.shells[:-1]))
        assert ri.attach_imported_tests(other_shells, case) == other_shells

    @pytest.mark.parametrize("key, value", [
        ("ae_atom", "20  0  5"),                  # exfact 4 -> 5
        ("ae_atom", "20.5  0.5  4"),              # z
        ("ae_mesh", "0  0  0  0"),                # rmax, aasf, bbsf and irel
        ("ae_mesh", "90  6  40  2"),              # rmax only
    ], ids=["exfact", "z", "mesh-irel", "rmax"])
    def test_a_changed_generation_gets_derived_tests(self, key, value):
        # runatom refuses a test whose all-electron input disagrees with the generation it tests
        case = ri.read_case_folder(LIBRARY_DIR / "020-Ca" / "020-Ca-gpw-nsp-campos")
        changed = ri.case_from_fields(dict(ri.case_to_fields(case), **{key: value}))
        attached = ri.attach_imported_tests(changed, case)
        assert all(t.ae_input is None and t.test_input is None for t in attached.tests)

    def test_a_changed_partition_gets_derived_tests(self):
        case = ri.read_case_folder(LIBRARY_DIR / "020-Ca" / "020-Ca-gpw-nsp-campos")
        derived = ri.case_from_fields(ri.case_to_fields(case))
        for ps in (ri.replace(derived.ps, ncores=derived.ps.ncores + 1, nvales=derived.ps.nvales - 1),
                   ri.replace(derived.ps, nang=derived.ps.nang - 1)):
            attached = ri.attach_imported_tests(ri.replace(derived, ps=ps), case)
            assert all(t.ae_input is None and t.test_input is None for t in attached.tests)

    def test_tests_with_another_shell_list_are_not_attached(self):
        case = ri.read_case_folder(LIBRARY_DIR / "011-Na" / "011-Na-ca-n-campos")
        own = next(t for t in case.tests if t.name in ri.unrepresentable_tests(case))
        typed = ri.replace(ri.case_from_fields(ri.case_to_fields(case)), tests=[ri.TestConfig(own.name, own.xion, {})])
        assert ri.attach_imported_tests(typed, case).tests[0].ae_input is None

    def test_tests_with_their_own_shell_list_are_named_not_lost_silently(self):
        case = ri.read_case_folder(LIBRARY_DIR / "011-Na" / "011-Na-ca-n-campos")
        skipped = ri.unrepresentable_tests(case)
        assert skipped and all(name not in ri.case_to_fields(case)["tests"].split() for name in skipped)


class TestChecks:
    def test_mesh_size_matches_runatom(self):
        # Ti example: rinit printed "value of mesh generated in rinit is   851"
        assert ri.radial_mesh_size(22.0, 200.0, 6.0, 59.0) == 851

    @pytest.mark.parametrize("folder", [p for p in CASE_DIRS if not (p / "V_ps.adat").exists()],
                             ids=lambda p: p.name)
    def test_no_library_case_is_refused(self, folder):
        # runatom accepted every one of these, so a refusal is a wrong rule
        assert ri.check_case(ri.read_case_folder(folder)) == []

    def _ti(self, **ps_changes):
        case = _ti_case()
        return ri.replace(case, ps=ri.replace(case.ps, **ps_changes))

    def test_gga_requires_ifqopt_3_even_for_minus_one(self):
        case = _ti_case()
        gga = ri.replace(case, ae=ri.replace(case.ae, exfact=5.0), ps=ri.replace(case.ps, ifqopt=-1))
        assert any("ifqopt = 3" in p for p in ri.check_case(gga))

    def test_occupations_must_add_up(self):
        case = _ti_case()
        bad = ri.replace(case, ae=ri.replace(case.ae, xion=0.0))
        assert any("add up to" in p for p in ri.check_case(bad))

    def test_projector_rules(self):
        betas = _ti_case().ps.betas
        assert any("non-decreasing" in p for p in ri.check_case(self._ti(betas=list(reversed(betas)))))
        assert any("keyee < 0" in p for p in ri.check_case(self._ti(betas=[ri.replace(betas[0], keyee=-1)] + betas[1:])))
        assert any("npf" in p for p in ri.check_case(self._ti(npf=4)))

    def test_mesh_limit(self):
        case = _ti_case()
        assert any("radial mesh" in p for p in ri.check_case(ri.replace(case, ae=ri.replace(case.ae, rmax=1e6))))

    def test_log_derivative_index_matches_runatom(self):
        # runatom.exe on the starter H input with rlogd 0: "klogd = 2 not allowed with mesh = 417"
        assert ri.radial_mesh_size(1.0, 0.0, 0.0, 0.0) == 417
        assert ri.log_derivative_index(1.0, 0.0, 0.0, 0.0) == 2

    def test_rlogd_must_give_a_usable_log_derivative_point(self):
        # rinit stops when klogd < 10 or klogd == mesh, whether or not log derivatives are asked for
        case = _ti_case()
        for rlogd in (0.0, 1e-6):
            ae = ri.replace(case.ae, controls=ri.replace(case.ae.controls, rlogd=rlogd))
            assert any("rlogd" in p for p in ri.check_ae(ae))
        ps = ri.replace(case.ps, controls=ri.replace(case.ps.controls, rlogd=0.0), rinner1=0.0,
                        rinner=[0.0] * 5)
        assert any("rlogd" in p for p in ri.check_ps(ps, case.ae))

    def test_long_file_names_are_refused(self):
        case = ri.replace(_ti_case(), generation_name="a" * 40)
        assert any("40 characters" in p for p in ri.check_case(case))

    def test_formatted_file_name_must_fit_reform(self):
        case = ri.replace(_ti_case(), symbol="x" * 19)        # x..x_ps_formatted.uspp: 37 characters
        problems = ri.check_case(case)
        assert len(problems) == 1 and "36" in problems[0] and "_formatted.uspp" in problems[0]
        assert ri.check_case(ri.replace(_ti_case(), symbol="x" * 18)) == []

    def test_a_tuned_test_is_checked_on_its_own_input(self):
        case = ri.read_case_folder(LIBRARY_DIR / "020-Ca" / "020-Ca-gpw-nsp-campos")
        test = next(t for t in case.tests if t.ae_input is not None)
        bad = ri.replace(test, ae_input=ri.replace(test.ae_input, xion=test.ae_input.xion + 1.0))
        problems = ri.check_case(ri.replace(case, tests=[bad]))
        assert any(p.startswith(f"Test {test.name}:") and "add up to" in p for p in problems)

    def test_unknown_shell_in_a_test(self):
        case = _ti_case()
        case.tests[0].occupations[430] = 1.0
        assert any("430" in p for p in ri.check_case(case))
