"""runatomOutput: the .out report and the .uspp file of runatom."""

import struct
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import runatomOutput as ro
from runatom_fixtures import REFERENCE_CASES, case_info, reference_bytes, reference_text

SRC = Path(__file__).resolve().parents[1]
TI = "000-Ex-ca-sp-test"


def test_module_does_not_import_toga():
    code = "import sys, runatomOutput; assert 'toga' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], cwd=SRC, check=True)


class TestParseOut:
    def test_ti_generation_report(self):
        report = ro.parse_out(reference_text(TI, "ti_ps.out"))
        assert [r.nlm for r in report.comparison] == [300, 310, 400, 320]
        assert report.comparison[0].all_electron == pytest.approx(-5.530317, abs=1e-6)
        assert [r.l for r in report.comparison] == [0, 1, 0, 2]
        assert sorted(report.bessel) == [0, 1, 2] and len(report.bessel[0]) == 4
        assert report.zv == 12.0
        assert report.errors == []
        assert any("max q eigval" in w for w in report.warnings)

    def test_ti_transferability_test(self):
        report = ro.parse_out(reference_text(TI, "ti_test_s2d2.out"))
        assert report.comparison[1].diff == pytest.approx(-0.000707, abs=1e-6)
        assert report.total_energy is not None

    def test_program_error_is_collected_with_its_detail(self):
        text = "something\n ***error in subroutine readin\n nang =  5 but idim5 =  4\n"
        assert ro.parse_out(text).errors == ["***error in subroutine readin nang =  5 but idim5 =  4"]


class TestGhosts:
    def test_no_ghost_in_ti(self):
        assert ro.ghost_candidates(ro.parse_out(reference_text(TI, "ti_ps.out"))) == []

    def test_extra_bound_state_is_a_ghost(self):
        report = ro.parse_out(reference_text(TI, "ti_ps.out"))
        report.bessel[2] = [-3.0] + report.bessel[2]     # one more d state than the AE atom has
        assert ro.ghost_candidates(report) == [2]


class TestReadUspp:
    @pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
    def test_every_reference_pseudopotential_reads(self, name, tmp_path):
        atom = case_info(name)["atom"]
        path = tmp_path / "p.uspp"
        path.write_bytes(reference_bytes(name, f"{atom}_ps.uspp"))
        data = ro.read_uspp(path)
        assert data.version == (7, 3, 6)
        assert data.beta.shape == (data.nbeta, data.kkbeta)
        assert data.qfunc.shape == (data.nbeta, data.nbeta, data.kkbeta)
        assert np.allclose(data.qqq, data.qqq.T)
        assert data.r.shape == (data.mesh,) and data.r[0] == 0.0

    def test_ti_header(self, tmp_path):
        path = tmp_path / "ti.uspp"
        path.write_bytes(reference_bytes(TI, "ti_ps.uspp"))
        data = ro.read_uspp(path)
        assert (data.title, data.zv, data.mesh, data.nbeta, data.kkbeta) == ("titanium", 12.0, 851, 6, 585)
        assert list(data.l) == [0, 0, 1, 1, 2, 2] and list(data.iptype) == [2] * 6

    def test_zero_filled_file_is_refused(self, tmp_path):
        path = tmp_path / "zeros.uspp"           # like the upstream Nb file
        path.write_bytes(bytes(10240))
        with pytest.raises(ValueError, match="record 1"):
            ro.read_uspp(path)

    def test_mismatched_record_markers_are_refused(self):
        with pytest.raises(ValueError, match="not a gfortran unformatted file"):
            ro.uspp_records(struct.pack("<i", 4) + b"abcd" + struct.pack("<i", 5))


def test_summary_flags_large_test_errors_and_reports_the_file(tmp_path):
    path = tmp_path / "ti.uspp"
    path.write_bytes(reference_bytes(TI, "ti_ps.uspp"))
    text = ro.job_summary([("Pseudopotential generation", 0.4, ro.parse_out(reference_text(TI, "ti_ps.out"))),
                           ("Transferability test (s2d2)", 0.2, ro.parse_out(reference_text(TI, "ti_test_s2d2.out")))],
                          ro.read_uspp(path))
    assert "no extra bound state" in text and "zv = 12" in text
    assert "more than 1 mRy" not in text      # every Ti test difference is below 1 mRy
