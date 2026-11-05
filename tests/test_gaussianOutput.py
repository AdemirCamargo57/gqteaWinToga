"""Optimized geometries from Gaussian 16 outputs (gaussianOutput, no Toga)."""
import math
import subprocess
import sys
from pathlib import Path

import pytest

import gaussianOutput as go

SRC = Path(__file__).resolve().parents[1]
DATA = SRC / "tests" / "data" / "gaussian_output"


def test_module_does_not_import_toga():
    out = subprocess.run([sys.executable, "-c", "import sys, gaussianOutput; print('toga' in sys.modules)"],
                         capture_output=True, text=True, cwd=str(SRC))
    assert out.stdout.strip() == "False", out.stderr


def text_of(name: str) -> str:
    return (DATA / name).read_text()


def converged_block_start(lines):
    done = next(k for k, line in enumerate(lines) if "Optimization completed" in line)
    return max(k for k in range(done) if "Input orientation:" in lines[k])


class TestFixtures:
    def test_converged_optimization(self):
        g = go.read_optimized_geometry(DATA / "h2o_opt.log")
        assert [a[0] for a in g.atoms] == ["O", "H", "H"]
        assert (g.charge, g.multiplicity, g.job_number, g.job_count) == (0, 1, 1, 1)
        assert g.orientation == "Input orientation" and g.notes == []
        oh = math.dist(g.atoms[0][1:], g.atoms[1][1:])
        assert 0.95 < oh < 0.98          # HF/3-21G water: 0.967 A (the start was 1.00 A)

    def test_the_converged_geometry_is_the_last_block_before_completion(self):
        text = text_of("h2o_opt.log")
        lines = text.splitlines()
        start = converged_block_start(lines)
        g = go.parse_optimized_geometry(text, "h2o_opt.log")
        for row, atom in zip(lines[start + 5:start + 8], g.atoms):
            assert atom[1:] == tuple(float(v) for v in row.split()[3:6])

    def test_opt_freq_out_ignores_the_freq_job(self):
        g = go.read_optimized_geometry(DATA / "h2o_opt_freq.out")
        assert (g.job_number, g.job_count) == (1, 2) and g.notes == []

    def test_a_later_failed_optimization_is_reported(self):
        g = go.read_optimized_geometry(DATA / "h2o_opt_then_failed.log")
        assert (g.job_number, g.job_count, g.charge, g.multiplicity) == (1, 2, 0, 1)
        assert len(g.notes) == 1 and "job 2" in g.notes[0] and "Number of steps exceeded" in g.notes[0]

    def test_maxcycles_is_refused(self):
        with pytest.raises(ValueError, match="no successfully optimized geometry.*Number of steps exceeded"):
            go.read_optimized_geometry(DATA / "h2co_opt_maxcycles.log")

    def test_opt_then_td_gives_the_optimization(self):
        g = go.read_optimized_geometry(DATA / "h2co_opt_td.log")
        assert [a[0] for a in g.atoms] == ["C", "O", "H", "H"] and (g.job_number, g.job_count) == (1, 2)
        assert g.notes == []

    def test_unreadable_file(self, tmp_path):
        with pytest.raises(ValueError, match="could not be opened"):
            go.read_optimized_geometry(tmp_path / "missing.log")


class TestSynthetic:
    def test_cut_off_before_completion(self):
        text = text_of("h2o_opt.log")
        with pytest.raises(ValueError, match="ends before the optimization completed"):
            go.parse_optimized_geometry(text[:text.index("Optimization completed")], "cut.log")

    def test_completed_but_no_normal_termination_is_a_note(self):
        text = text_of("h2o_opt.log")
        g = go.parse_optimized_geometry(text[:text.index("Normal termination")], "x.log")
        assert len(g.notes) == 1 and "Normal termination" in g.notes[0]

    def test_not_a_gaussian_output(self):
        with pytest.raises(ValueError, match="not a Gaussian output"):
            go.parse_optimized_geometry("* O   R   A   *\nORCA TERMINATED NORMALLY\n", "orca.out")

    def test_no_optimization(self):
        text = text_of("h2co_opt_td.log")
        second_job = text[text.index("Normal termination"):].split("\n", 1)[1]
        with pytest.raises(ValueError, match="contains no geometry optimization"):
            go.parse_optimized_geometry(second_job, "td.log")

    def test_ghost_atoms_are_refused(self):
        lines = text_of("h2o_opt.log").splitlines()
        start = converged_block_start(lines)
        row = lines[start + 6].split()
        lines[start + 6] = f"      {row[0]}          0           0     {row[3]}    {row[4]}    {row[5]}"
        with pytest.raises(ValueError, match="ghost or dummy"):
            go.parse_optimized_geometry("\n".join(lines), "x.log")

    def test_relaxed_scan_note(self):
        text = text_of("h2o_opt.log").replace(" Optimization completed.",
                                              " Optimization completed.\n Optimization completed.", 1)
        g = go.parse_optimized_geometry(text, "scan.log")
        assert any("relaxed scan" in n for n in g.notes)

    def test_standard_orientation_when_no_input_orientation(self):
        text = text_of("h2o_opt.log").replace("Input orientation:", "Other heading:")
        g = go.parse_optimized_geometry(text, "x.log")
        assert g.orientation == "Standard orientation" and len(g.atoms) == 3

    @pytest.mark.parametrize("route,optimizes", [
        ("#P HF/3-21G Opt", True), ("#p opt=(maxcycles=2) hf/3-21g", True), ("# B3LYP/6-31G Opt(TS,CalcFC)", True),
        ("#P Geom=AllCheck Guess=TCheck SCRF=Check GenChk RHF/3-21G Freq", False),
        ("#P HF/3-21G Stable=Opt", False), ("#P HF/3-21G IOp(2/17=4)", False), ("#P HF/3-21G", False),
    ])
    def test_which_routes_optimize(self, route, optimizes):
        assert go._Job(1, [], route, [], True).optimizes is optimizes
