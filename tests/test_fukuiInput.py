"""Tests for the Gaussian Fukui input builder core (fukuiInput, no Toga)."""
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

import fukuiInput as fi
from fukui_fixtures import N

PHENOL_INPUT = Path(__file__).resolve().parent / "data" / "fukui" / "inputs" / "phenol_N.gjf"


def phenol_geometry() -> str:
    lines = PHENOL_INPUT.read_text().splitlines()
    start = lines.index("0 1") + 1
    return "\n".join(l for l in lines[start:] if l.strip())


def make_spec(tmp_path, **changes) -> fi.FukuiJobSpec:
    spec = fi.FukuiJobSpec(stem="phenol", output_dir=str(tmp_path), geometry=phenol_geometry(),
                           charge=0, multiplicity=1, plus_multiplicity=2, minus_multiplicity=2,
                           memory="2GB", nproc=4)
    return replace(spec, **changes)


def jobs_of(text: str):
    """The three linked jobs of an input as lists of lines."""
    return [block.strip("\n").split("\n") for block in text.split("--Link1--")]


def charge_index(job_lines):
    blanks = [i for i, l in enumerate(job_lines) if not l.strip()]
    return blanks[1] + 1   # route, blank, title, blank, charge line


def charge_line(job_lines):
    return job_lines[charge_index(job_lines)]


def molecule_lines(job_lines):
    """Lines after the charge/multiplicity line, up to the blank line or the end of the job."""
    out = []
    for line in job_lines[charge_index(job_lines) + 1:]:
        if not line.strip():
            break
        out.append(line)
    return out


def route_of(job_lines):
    return next(l for l in job_lines if l.startswith("#"))


# ----------------------------------------------------------------------
# Field parsers
# ----------------------------------------------------------------------
class TestParsers:
    @pytest.mark.parametrize("text,expected", [("2GB", "2GB"), ("500mb", "500MB"), (" 16 GB ", "16GB"),
                                               ("100MW", "100MW")])
    def test_memory_accepts_integer_with_unit(self, text, expected):
        assert fi.parse_memory(text) == expected

    def test_memory_bare_number_is_refused_because_gaussian_reads_words(self):
        with pytest.raises(ValueError, match="words"):
            fi.parse_memory("2000")

    @pytest.mark.parametrize("text", ["", "1.5GB", "0GB", "2 gigabytes", "-1GB"])
    def test_memory_rejects_other_forms(self, text):
        with pytest.raises(ValueError, match="%Mem"):
            fi.parse_memory(text)

    def test_nproc(self):
        assert fi.parse_nproc("8") == 8
        for bad in ("0", "-2", "two", "", "1.5"):
            with pytest.raises(ValueError, match="%NProcShared"):
                fi.parse_nproc(bad)

    @pytest.mark.parametrize("stem", ["phenol", "2-criegee", "mol_1.b"])
    def test_stem_accepted(self, stem):
        assert fi.validate_stem(stem) == stem

    @pytest.mark.parametrize("stem", ["", "my molecule", "a/b", "é", "-lead", "name%"])
    def test_stem_rejected(self, stem):
        with pytest.raises(ValueError, match="Job name"):
            fi.validate_stem(stem)

    @pytest.mark.parametrize("text,expected", [("Opt", "Opt"), (" opt=tight ", "opt=tight"),
                                               ("Opt=(Tight, CalcFC)", "Opt=(Tight, CalcFC)")])
    def test_optimization_keyword(self, text, expected):
        assert fi.parse_opt_keyword(text) == expected

    @pytest.mark.parametrize("text", ["", "Freq", "Opt Freq", "Opt=ReadFC", "Opt=ModRedundant", "Opt=TS"])
    def test_optimization_keyword_rejected(self, text):
        with pytest.raises(ValueError, match="Optimization"):
            fi.parse_opt_keyword(text)


class TestGeometry:
    def test_cartesian_lines(self):
        atoms = fi.parse_geometry("C 0 0 0\nh 1.09 0 0\n\n")
        assert [a[0] for a in atoms] == ["C", "H"] and atoms[1][1] == pytest.approx(1.09)

    def test_atomic_numbers_and_numbered_labels(self):
        atoms = fi.parse_geometry("8 0 0 0\nH1 0.96 0 0\nCL12 0 2 0")
        assert [a[0] for a in atoms] == ["O", "H", "Cl"]

    def test_full_xyz_file_header_is_skipped(self):
        atoms = fi.parse_geometry("2\nwater fragment\nO 0 0 0\nH 0.96 0 0\n")
        assert [a[0] for a in atoms] == ["O", "H"]

    def test_empty_geometry(self):
        with pytest.raises(ValueError, match="geometry"):
            fi.parse_geometry("  \n")

    @pytest.mark.parametrize("text,fragment", [
        ("C\nH 1 1.09", "line 1.*Z-matrix"),     # Z-matrix
        ("C 0 0 0\nH 0 1.0 0 0", "line 2"),      # five columns (freeze flag)
        ("C 0 0 x", "line 1"),
        ("C 0 0 nan", "line 1"),
    ])
    def test_malformed_lines_name_the_line(self, text, fragment):
        with pytest.raises(ValueError, match=fragment):
            fi.parse_geometry(text)

    def test_unknown_element(self):
        with pytest.raises(ValueError, match="Xx"):
            fi.parse_geometry("Xx 0 0 0")

    def test_overlapping_atoms(self):
        with pytest.raises(ValueError, match="atoms 1 and 2"):
            fi.parse_geometry("C 0 0 0\nH 0.05 0 0")

    def test_distance_check_can_be_skipped_for_live_counts(self):
        assert len(fi.parse_geometry("C 0 0 0\nH 0.05 0 0", check_distances=False)) == 2

    def test_format_is_fixed_width_angstrom(self):
        assert fi.format_geometry([("C", 0.020821, -1.857334, 0.0)]) == [
            "C       0.02082100     -1.85733400      0.00000000"]  # the fixture .gjf layout


class TestRoute:
    def test_default_route(self):
        assert fi.build_route("#P", "B3LYP", "6-31G(d)", "SCF=Tight NoSymm") == [
            "#P B3LYP/6-31G(d) SCF=Tight NoSymm"]

    def test_long_route_is_wrapped(self):
        lines = fi.build_route("#P", "B3LYP", "6-31G(d)", " ".join(f"IOp(3/{i}=1)" for i in range(20, 40)))
        assert len(lines) > 1 and all(len(l) <= fi.ROUTE_WIDTH for l in lines) and lines[0].startswith("#P ")

    def test_keyword_tokens_keep_parentheses(self):
        assert fi.route_tokens("SCF=(Tight, XQC) NoSymm") == ["SCF=(Tight, XQC)", "NoSymm"]

    @pytest.mark.parametrize("keywords,word", [
        ("Opt", "Optimization"), ("opt=tight", "Optimization"), ("IRC=CalcFC", "IRC"), ("Scan", "Scan"),
        ("Geom=Check", "Geom"), ("Guess=Mix", "Guess"), ("Units=AU", "Units"),
        ("SCRF=(PCM,Read)", "extra input"), ("Pseudo=Read", "Pseudo"), ("ExtraBasis", "ExtraBasis"),
    ])
    def test_forbidden_keywords(self, keywords, word):
        errors, _ = fi.check_route("B3LYP", "6-31G(d)", keywords)
        assert errors and word.lower() in " ".join(errors).lower()

    @pytest.mark.parametrize("basis", ["Gen", "GenECP", "gen"])
    def test_general_basis_is_refused(self, basis):
        errors, _ = fi.check_route("B3LYP", basis, "")
        assert errors and "Gen" in errors[0]

    def test_method_with_basis_attached(self):
        errors, _ = fi.check_route("B3LYP/6-31G(d)", "", "")
        assert errors and "Basis" in errors[0]

    def test_stable_opt_is_not_an_optimisation(self):
        errors, _ = fi.check_route("B3LYP", "6-31G(d)", "Stable=Opt")
        assert errors == []

    @pytest.mark.parametrize("method,keywords", [("MP2", ""), ("CCSD(T)", ""), ("B3LYP", "TD=(NStates=5)"),
                                                 ("UMP2", "")])
    def test_post_scf_is_a_warning(self, method, keywords):
        errors, warnings = fi.check_route(method, "6-31G(d)", keywords)
        assert errors == [] and any("Total SCF Density" in w for w in warnings)

    def test_freq_and_diffuse_are_warnings(self):
        errors, warnings = fi.check_route("B3LYP", "6-311+G(d,p)", "Freq")
        assert errors == [] and len(warnings) == 2


# ----------------------------------------------------------------------
# The linked job
# ----------------------------------------------------------------------
class TestPrepare:
    def test_states_charges_and_electrons(self, tmp_path):
        job = fi.FukuiInputBuilder().prepare(make_spec(tmp_path))
        assert [(s.label, s.tag, s.charge, s.multiplicity, s.n_electrons) for s in job.states] == [
            ("N", "N", 0, 1, 50), ("N+1", "Np1", -1, 2, 51), ("N-1", "Nm1", 1, 2, 49)]

    def test_p_and_q_states(self, tmp_path):
        job = fi.FukuiInputBuilder().prepare(make_spec(tmp_path, p=2, q=2, plus_multiplicity=3,
                                                       minus_multiplicity=3))
        assert [(s.label, s.tag, s.charge) for s in job.states] == [
            ("N", "N", 0), ("N+2", "Np2", -2), ("N-2", "Nm2", 2)]

    def test_paths_are_relative_for_linux(self, tmp_path):
        job = fi.FukuiInputBuilder().prepare(make_spec(tmp_path))
        assert job.input_path == tmp_path / "phenol_fukui.gjf"
        assert job.script_path == tmp_path / "phenol_formchk.sh"
        assert [s.chk for s in job.states] == ["phenol_N.chk", "phenol_Np1.chk", "phenol_Nm1.chk"]
        assert [s.fchk for s in job.states] == [tmp_path / f"phenol_{t}.fchk" for t in ("N", "Np1", "Nm1")]

    def test_folder_with_spaces_is_fine_with_relative_names(self, tmp_path):
        folder = tmp_path / "my jobs"
        folder.mkdir()
        assert fi.FukuiInputBuilder().prepare(make_spec(folder)).input_path.parent == folder

    def test_missing_output_folder(self, tmp_path):
        with pytest.raises(ValueError, match="Output folder"):
            fi.FukuiInputBuilder().prepare(make_spec(tmp_path / "missing"))
        with pytest.raises(ValueError, match="Output folder"):
            fi.FukuiInputBuilder().prepare(make_spec(tmp_path, output_dir=" "))

    @pytest.mark.parametrize("field,pattern", [("multiplicity", r"N state"),
                                               ("plus_multiplicity", r"N\+1 state"),
                                               ("minus_multiplicity", r"N-1 state")])
    def test_missing_multiplicity_is_asked_for_not_assumed(self, tmp_path, field, pattern):
        with pytest.raises(ValueError, match=pattern + r".*not assumed"):
            fi.FukuiInputBuilder().prepare(make_spec(tmp_path, **{field: None}))

    @pytest.mark.parametrize("changes,pattern", [
        ({"multiplicity": 2}, r"N state"),            # 50 electrons: odd multiplicity only
        ({"plus_multiplicity": 1}, r"N\+1 state"),    # 51 electrons: even only
        ({"minus_multiplicity": 3}, r"N-1 state"),    # 49 electrons: even only
        ({"multiplicity": 0}, r"N state"),
    ])
    def test_multiplicity_problems(self, tmp_path, changes, pattern):
        with pytest.raises(ValueError, match=pattern):
            fi.FukuiInputBuilder().prepare(make_spec(tmp_path, **changes))

    def test_n_minus_one_needs_an_electron(self, tmp_path):
        spec = make_spec(tmp_path, geometry="H 0 0 0", multiplicity=2, plus_multiplicity=1,
                         minus_multiplicity=1)
        with pytest.raises(ValueError, match="no electrons"):
            fi.FukuiInputBuilder().prepare(spec)

    def test_all_problems_are_reported_together(self, tmp_path):
        spec = make_spec(tmp_path, memory="2000", nproc=0, keywords="Opt", opt_keyword="Freq")
        with pytest.raises(ValueError) as info:
            fi.FukuiInputBuilder().prepare(spec)
        text = str(info.value)
        assert "%Mem" in text and "%NProcShared" in text and "Keyword Opt" in text and "Freq" in text

    def test_bad_print_level(self, tmp_path):
        with pytest.raises(ValueError, match="Print level"):
            fi.FukuiInputBuilder().prepare(make_spec(tmp_path, print_level="#X"))

    def test_restricted_prefix_with_open_shell_ion_warns(self, tmp_path):
        job = fi.FukuiInputBuilder().prepare(make_spec(tmp_path, method="RB3LYP"))
        assert any("RB3LYP" in w and "N+1" in w for w in job.warnings)

    def test_unusual_ion_multiplicity_warns(self, tmp_path):
        job = fi.FukuiInputBuilder().prepare(make_spec(tmp_path, plus_multiplicity=4))
        assert any("N+1" in w and "multiplicity 4" in w for w in job.warnings)

    def test_no_warnings_by_default(self, tmp_path):
        assert fi.FukuiInputBuilder().prepare(make_spec(tmp_path)).warnings == []


class TestInputText:
    def text(self, tmp_path, **changes):
        builder = fi.FukuiInputBuilder()
        return builder.input_text(builder.prepare(make_spec(tmp_path, **changes)))

    def test_three_linked_jobs(self, tmp_path):
        text = self.text(tmp_path)
        assert text.count("--Link1--") == 2
        assert all(line == "--Link1--" for line in text.split("\n") if "Link1" in line)
        assert len(jobs_of(text)) == 3

    def test_job_one_optimizes_the_n_state(self, tmp_path):
        job = jobs_of(self.text(tmp_path))[0]
        assert job[:4] == ["%Chk=phenol_N.chk", "%Mem=2GB", "%NProcShared=4",
                           "#P B3LYP/6-31G(d) SCF=Tight NoSymm Opt"]
        assert charge_line(job) == "0 1" and len(molecule_lines(job)) == 13

    def test_optimization_options_are_kept(self, tmp_path):
        job = jobs_of(self.text(tmp_path, opt_keyword="Opt=(Tight,CalcFC)"))[0]
        assert route_of(job).endswith(" Opt=(Tight,CalcFC)")

    def test_ions_read_the_optimized_geometry_into_their_own_checkpoints(self, tmp_path):
        for job, tag in zip(jobs_of(self.text(tmp_path))[1:], ("Np1", "Nm1")):
            assert job[:4] == ["%OldChk=phenol_N.chk", f"%Chk=phenol_{tag}.chk", "%Mem=2GB", "%NProcShared=4"]
            assert route_of(job) == "#P B3LYP/6-31G(d) SCF=Tight NoSymm Geom=Check Guess=Read"
            assert "Opt" not in route_of(job)
            assert molecule_lines(job) == []   # Geom=Check: no coordinates

    def test_guess_read_can_be_switched_off(self, tmp_path):
        for job in jobs_of(self.text(tmp_path, guess_read=False))[1:]:
            assert route_of(job).endswith(" Geom=Check") and "Guess" not in route_of(job)

    def test_charge_and_multiplicity_lines_match_the_fixture_jobs(self, tmp_path):
        # tests/data/fukui/inputs/phenol_{N,Np1,Nm1}.gjf: 0 1, -1 2, 1 2
        assert [charge_line(j) for j in jobs_of(self.text(tmp_path))] == ["0 1", "-1 2", "1 2"]

    def test_titles_avoid_characters_gaussian_forbids(self, tmp_path):
        text = self.text(tmp_path, stem="my_mol-1")
        titles = [j[j.index("") + 1] for j in jobs_of(text)]
        assert titles == ["my mol 1 Fukui N state optimization charge 0 multiplicity 1",
                          "my mol 1 Fukui N plus 1 state charge minus 1 multiplicity 2",
                          "my mol 1 Fukui N minus 1 state charge 1 multiplicity 2"]
        assert all(not set(t) & set("@#!-_\\") for t in titles)

    def test_file_ends_with_a_blank_line(self, tmp_path):
        assert self.text(tmp_path).endswith("\n\n")

    def test_job_one_reads_back_through_the_coordinate_converter(self, tmp_path):
        from coordinateConverter import CoordinateConverter
        reference = fi.parse_geometry(phenol_geometry())
        path = tmp_path / "job1.gjf"
        path.write_text("\n".join(jobs_of(self.text(tmp_path))[0]) + "\n\n")
        conv = CoordinateConverter()
        conv.set_input_file(str(path))
        atoms = conv.read_atoms("Gaussian")
        assert [a["symbol"] for a in atoms] == [r[0] for r in reference]
        for axis, column in (("x", 1), ("y", 2), ("z", 3)):
            assert [a[axis] for a in atoms] == pytest.approx([r[column] for r in reference], abs=1e-8)


class TestWriteFiles:
    def written(self, tmp_path, **changes):
        builder = fi.FukuiInputBuilder()
        job = builder.prepare(make_spec(tmp_path, **changes))
        return job, builder.write_files(job)

    def test_files_use_linux_line_endings(self, tmp_path):
        job, written = self.written(tmp_path)
        assert written == [job.input_path, job.script_path]
        for path in written:
            assert b"\r" not in path.read_bytes()
        assert b"\n--Link1--\n" in job.input_path.read_bytes()

    def test_script_converts_every_checkpoint(self, tmp_path):
        job, _ = self.written(tmp_path)
        script = job.script_path.read_text()
        assert script.startswith("#!/bin/bash\n")
        assert 'cd "$(dirname "$0")"' in script and "command -v formchk" in script
        for tag in ("N", "Np1", "Nm1"):
            assert f'convert "phenol_{tag}.chk" "phenol_{tag}.fchk"' in script
        assert script.rstrip().endswith("exit $status")

    def test_script_is_valid_bash(self, tmp_path):
        bash = shutil.which("bash")
        if bash is None:
            pytest.skip("no bash")
        job, _ = self.written(tmp_path)
        # Raw bytes: text=True would turn the piped LF into CRLF on Windows.
        done = subprocess.run([bash, "-n", "-"], input=job.script_path.read_bytes(), capture_output=True)
        assert done.returncode == 0, done.stderr.decode(errors="replace")

    def test_write_error_is_a_value_error(self, tmp_path):
        builder = fi.FukuiInputBuilder()
        job = builder.prepare(make_spec(tmp_path))
        job.input_path.mkdir()  # a folder where the file should go
        with pytest.raises(ValueError, match="Could not write"):
            builder.write_files(job)


class TestSiblings:
    def test_plus_and_minus_found_next_to_n(self, tmp_path):
        for s in ("N", "Np1", "Nm1"):
            shutil.copy(N, tmp_path / f"phenol_{s}.fchk")
        assert fi.sibling_state_files(tmp_path / "phenol_N.fchk") == (
            tmp_path / "phenol_Np1.fchk", tmp_path / "phenol_Nm1.fchk")

    def test_fch_extension_and_missing_siblings(self, tmp_path):
        (tmp_path / "m_N.fch").write_text("")
        (tmp_path / "m_Np1.fch").write_text("")
        assert fi.sibling_state_files(tmp_path / "m_N.fch") == (tmp_path / "m_Np1.fch", None)

    def test_not_an_n_file(self, tmp_path):
        assert fi.sibling_state_files(tmp_path / "phenol.fchk") == (None, None)


def test_module_does_not_import_toga():
    code = "import sys, fukuiInput; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(Path(__file__).resolve().parents[1]))
    assert out.stdout.strip() == "False", out.stderr
