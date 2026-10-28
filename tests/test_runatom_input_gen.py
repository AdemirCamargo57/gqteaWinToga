"""The runatom window, headless: stub widgets, fake program, no Toga window."""

import asyncio
import sys
import threading
from pathlib import Path

import pytest

import runatomInput as ri
import runatomRunner as rr
import runatom_input_gen as gen
from runatom_fixtures import LIBRARY_DIR

FAKE = [sys.executable, "-I", str(Path(__file__).resolve().parent / "fake_runatom.py")]
FAKE_REFORM = [sys.executable, "-I", str(Path(__file__).resolve().parent / "fake_reform.py")]


class Value:
    def __init__(self, value=""):
        self.value = value
        self.enabled = True


class Label:
    text = ""


class Spinner:
    running = False

    def start(self):
        self.running = True

    def stop(self):
        self.running = False


class Tabs:
    current_tab = 0


class Window:
    def __init__(self):
        self.dialogs = []

    async def dialog(self, d):
        self.dialogs.append(d)


def _window(tmp_path):
    w = gen.RunatomGenerator.__new__(gen.RunatomGenerator)
    w._closed, w._busy, w._cancel, w._last_out = False, False, None, None
    case = ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-test")
    w.inputs = {k: Value(v) for k, v in ri.case_to_fields(case).items()}
    w.folder_input, w.exe_input = Value(str(tmp_path)), Value("")
    w.makefile_switch = Value(False)
    w.results_text = Value("")
    w.status_label, w.activity, w.tabs, w.main_window = Label(), Spinner(), Tabs(), Window()
    for name in ("check_button", "save_button", "run_button", "cancel_button"):
        setattr(w, name, Value())
    w.errors = []                     # a Toga dialog keeps no readable message: record _error instead

    async def record(title, message):
        w.errors.append((title, message))

    w._error = record
    return w, case


def test_every_form_field_has_a_widget():
    assert set(ri.FORM_FIELDS) == set(gen.NAME_FIELDS) | {k for _, keys in gen.TAB_FIELDS for k in keys}


def test_collect_case_reads_the_widgets(tmp_path):
    w, case = _window(tmp_path)
    assert w._collect_case() == ri.case_from_fields(ri.case_to_fields(case))


def test_run_writes_inputs_runs_every_stage_and_shows_results(tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "ok")
    monkeypatch.setenv("FAKE_REFORM_MODE", "ok")
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: "runatom.exe")
    monkeypatch.setattr(gen.rr, "runatom_command", lambda exe: FAKE)
    monkeypatch.setattr(gen.rr, "find_reform", lambda explicit=None: "reform.exe")
    monkeypatch.setattr(gen.rr, "reform_command", lambda exe: FAKE_REFORM)
    w, _ = _window(tmp_path)
    asyncio.run(w.on_run(None))
    assert (tmp_path / "ti_ps.adat").exists() and (tmp_path / "ti_ps.uspp").exists()
    assert (tmp_path / "ti_ps_formatted.uspp").exists()
    assert "*ti_ps_formatted.uspp FORMATTED" in w.results_text.value
    assert "Pseudopotential generation: finished" in w.results_text.value
    assert w.tabs.current_tab == len(gen.TAB_FIELDS)
    assert w.errors == [] and w.status_label.text.startswith("Finished: 4 stages")
    assert not w._busy and w.run_button.enabled and not w.activity.running


def test_missing_reform_leaves_the_binary_result_and_says_so(tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "ok")
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: "runatom.exe")
    monkeypatch.setattr(gen.rr, "runatom_command", lambda exe: FAKE)
    monkeypatch.setattr(gen.rr, "find_reform", lambda explicit=None: None)
    w, _ = _window(tmp_path)
    asyncio.run(w.on_run(None))
    assert "reform.exe was not found" in w.results_text.value and w.errors == []
    assert not (tmp_path / "ti_ps_formatted.uspp").exists()


def test_program_error_is_shown_in_a_dialog(tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "error")
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: "runatom.exe")
    monkeypatch.setattr(gen.rr, "runatom_command", lambda exe: FAKE)
    w, _ = _window(tmp_path)
    asyncio.run(w.on_run(None))
    ((title, message),) = w.errors
    assert title == "runatom stopped" and "idim5" in message and w.run_button.enabled


def test_invalid_inputs_never_reach_the_program(tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr(gen.rr, "run_job", lambda *a, **k: called.append(a))
    w, _ = _window(tmp_path)
    w.inputs["ae_atom"].value = "22 0.0 0"            # occupations no longer add up
    asyncio.run(w.on_run(None))
    assert called == [] and "add up to" in w.errors[0][1]


def test_missing_program_is_explained(tmp_path, monkeypatch):
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: None)
    w, _ = _window(tmp_path)
    asyncio.run(w.on_run(None))
    assert "runatom.exe was not found" in w.errors[0][1]


def test_progress_is_thinned_and_marshalled(tmp_path):
    w, _ = _window(tmp_path)
    posted = []

    class Loop:
        def call_soon_threadsafe(self, fn, *args):
            posted.append(args)

    report = w._progress_reporter(Loop())
    info = {"stage_index": 0, "stage_count": 4, "stage": "All-electron (s2d1)", "iteration": None, "delta": None}
    report(dict(info, message="starting"))
    for k in range(50):
        report(dict(info, message="SCF", iteration=k, delta=0.1))
    assert 1 <= len(posted) <= 3 and posted[0][0]["message"] == "starting"


def test_an_unexpected_error_never_leaves_the_window_busy(tmp_path, monkeypatch):
    def boom(*args, **kwargs):
        raise KeyError("unexpected")
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: "runatom.exe")
    monkeypatch.setattr(gen.rr, "run_job", boom)
    w, _ = _window(tmp_path)
    asyncio.run(w.on_run(None))
    assert not w._busy and w.run_button.enabled and not w.activity.running
    assert w.errors and "unexpected" in w.errors[0][1]


def test_a_bare_drive_as_runatom_path_does_not_break_the_cpmd_step(tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "ok")
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: "runatom.exe")
    monkeypatch.setattr(gen.rr, "runatom_command", lambda exe: FAKE)
    monkeypatch.setattr(gen.rr, "find_reform", lambda explicit=None: None)
    w, _ = _window(tmp_path)
    w.exe_input.value = "C:\\"
    asyncio.run(w.on_run(None))
    assert "Pseudopotential generation: finished" in w.results_text.value and not w._busy


def test_closing_the_window_cancels_a_run(tmp_path):
    w, _ = _window(tmp_path)
    w._cancel = threading.Event()
    assert w._on_close(None) is True
    assert w._closed and w._cancel.is_set()
