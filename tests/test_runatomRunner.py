"""runatomRunner: stages, child processes, progress, cancel, timeout and errors (fake program)."""

import json
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

import runatomInput as ri
import runatomRunner as rr
from runatom_fixtures import LIBRARY_DIR

SRC = Path(__file__).resolve().parents[1]
FAKE = [sys.executable, "-I", str(Path(__file__).resolve().parent / "fake_runatom.py")]


def _ti():
    return ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-test")


def _job(tmp_path, mode, monkeypatch, **kwargs):
    monkeypatch.setenv("FAKE_RUNATOM_MODE", mode)
    case = _ti()
    ri.write_case(case, tmp_path)
    return rr.run_job(FAKE, rr.plan_stages(case), tmp_path, **kwargs)


def test_module_does_not_import_toga():
    code = "import sys, runatomRunner; assert 'toga' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], cwd=SRC, check=True)


def test_stages_reproduce_the_library_makefile():
    stages = rr.plan_stages(_ti())
    assert [s.arguments() for s in stages] == [
        ["ti_ae_s2d1.adat", "ti_ae_s2d1.out", "ti_ae_s2d1.ae", "ti_ae_s2d1.atwf", "ti_ae_s2d1.logd", "dummy"],
        ["ti_ps.adat", "ti_ps.out", "ti_ae_s2d1.ae", "ti_ps.atwf", "ti_ps.logd", "ti_ps.uspp"],
        ["ti_ae_s2d2.adat", "ti_ae_s2d2.out", "ti_ae_s2d2.ae", "ti_ae_s2d2.atwf", "ti_ae_s2d2.logd", "dummy"],
        ["ti_test_s2d2.adat", "ti_test_s2d2.out", "ti_ae_s2d2.ae", "ti_test_s2d2.atwf", "ti_test_s2d2.logd", "ti_ps.uspp"],
    ]


def test_successful_job_runs_every_stage_in_the_job_folder(tmp_path, monkeypatch):
    results = _job(tmp_path, "ok", monkeypatch)
    assert [r.stage.kind for r in results] == ["ae", "pseudo", "ae", "test"]
    calls = [json.loads(l) for l in (tmp_path / "fake_calls.jsonl").read_text().splitlines()]
    assert all(Path(c["cwd"]) == tmp_path for c in calls)
    assert results[1].report.comparison[0].nlm == 200
    assert (tmp_path / "ti_ps.uspp").exists()


def test_progress_reports_stages_and_scf_cycles(tmp_path, monkeypatch):
    seen = []
    _job(tmp_path, "ok", monkeypatch, progress=seen.append)
    assert {i["stage_index"] for i in seen} == {0, 1, 2, 3}
    cycles = [i for i in seen if i["iteration"] is not None]
    assert cycles[0]["iteration"] == 1 and cycles[0]["delta"] == pytest.approx(0.1)
    assert rr.progress_text(cycles[0]).startswith("Stage 1/4: All-electron (s2d1) - SCF cycle 1")


def test_a_failing_progress_callback_never_stops_the_run(tmp_path, monkeypatch):
    def broken(info):
        raise RuntimeError("UI gone")
    assert len(_job(tmp_path, "ok", monkeypatch, progress=broken)) == 4


@pytest.mark.parametrize("mode, kind, text", [
    ("error", "program_error", "nang =  5 but idim5 =  4"),
    ("runtime", "input_format", "ti_ae_s2d1.adat"),
    ("negative", "negative_density", "rinner"),
    ("stop0", "missing_output", "stopped before the end"),
])
def test_failures_are_classified(tmp_path, monkeypatch, mode, kind, text):
    with pytest.raises(rr.RunatomError) as info:
        _job(tmp_path, mode, monkeypatch)
    assert info.value.kind == kind and text in str(info.value)
    assert info.value.stage_label == "All-electron (s2d1)"


def test_stale_outputs_are_removed_before_a_stage(tmp_path, monkeypatch):
    (tmp_path / "ti_ae_s2d1.ae").write_bytes(b"old")
    with pytest.raises(rr.RunatomError):
        _job(tmp_path, "stop0", monkeypatch)
    assert not (tmp_path / "ti_ae_s2d1.ae").exists()


def test_cancel_kills_the_program_and_removes_its_outputs(tmp_path, monkeypatch):
    cancel = threading.Event()
    threading.Timer(0.5, cancel.set).start()
    started = time.monotonic()
    with pytest.raises(rr.RunatomCancelled):
        _job(tmp_path, "hang", monkeypatch, cancel=cancel)
    assert time.monotonic() - started < 5
    assert not (tmp_path / "ti_ae_s2d1.out").exists() and not (tmp_path / "ti_ae_s2d1.ae").exists()


def test_timeout_stops_a_hung_program(tmp_path, monkeypatch):
    with pytest.raises(rr.RunatomError) as info:
        _job(tmp_path, "hang", monkeypatch, timeout=0.5)
    assert info.value.kind == "timeout"


def test_names_longer_than_fortran_accepts_are_refused(tmp_path):
    stage = rr.Stage("ae", "x", "a" * 40, "a.ae", "dummy")
    with pytest.raises(ValueError, match="40"):
        rr.run_stage(FAKE, stage, tmp_path)


class TestFindRunatom:
    def test_explicit_path_first(self, tmp_path, monkeypatch):
        exe = tmp_path / "my.exe"
        exe.write_bytes(b"")
        monkeypatch.setattr(rr, "VENDOR_EXE", tmp_path / "missing.exe")
        assert rr.find_runatom(str(exe)) == str(exe)

    def test_environment_variable(self, tmp_path, monkeypatch):
        exe = tmp_path / "env.exe"
        exe.write_bytes(b"")
        monkeypatch.setenv(rr.ENV_VAR, str(exe))
        assert rr.find_runatom() == str(exe)

    def test_frozen_build_looks_in_the_bundle(self, tmp_path, monkeypatch):
        exe = tmp_path / rr.FROZEN_SUBDIR / "runatom.exe"
        exe.parent.mkdir()
        exe.write_bytes(b"")
        monkeypatch.delenv(rr.ENV_VAR, raising=False)
        monkeypatch.setattr(sys, "frozen", True, raising=False)
        monkeypatch.setattr(sys, "_MEIPASS", str(tmp_path), raising=False)
        assert rr.find_runatom() == str(exe)
        assert rr.library_dir() == tmp_path / rr.FROZEN_SUBDIR / "Work"

    def test_nothing_found(self, tmp_path, monkeypatch):
        monkeypatch.delenv(rr.ENV_VAR, raising=False)
        monkeypatch.setattr(rr, "VENDOR_EXE", tmp_path / "missing.exe")
        monkeypatch.setattr(rr.shutil, "which", lambda name: None)
        assert rr.find_runatom() is None
