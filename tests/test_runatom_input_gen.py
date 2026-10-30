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
        self.answer = True            # what a QuestionDialog returns

    async def dialog(self, d):
        self.dialogs.append(d)
        return self.answer if isinstance(d, gen.toga.QuestionDialog) else None


class Selection:
    def __init__(self):
        self.items, self.value, self.enabled, self.writes = [], None, True, 0

    def __setattr__(self, name, value):
        if name in ("items", "value") and hasattr(self, "writes"):
            self.writes += 1
        object.__setattr__(self, name, value)


class Controller(gen.RunatomElementController):
    def __init__(self, on_select):
        super().__init__(on_select)
        self.shown, self.closed = [], 0

    def show_element(self, symbol):
        super().show_element(symbol)
        self.shown.append(symbol)

    def close_picker(self):
        self.closed += 1


def _window(tmp_path):
    w = gen.RunatomGenerator.__new__(gen.RunatomGenerator)
    w._closed, w._busy, w._cancel, w._last_out = False, False, None, None
    w._imported_case = None
    w._init_defaults_state()
    w.element_controller = Controller(w.element_selected)
    case = ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-test")
    w.inputs = {k: Value(v) for k, v in ri.case_to_fields(case).items()}
    w._baseline = w._current_fields()
    w.folder_input, w.exe_input = Value(str(tmp_path)), Value("")
    w.makefile_switch = Value(False)
    w.results_text = Value("")
    w.status_label, w.activity, w.tabs, w.main_window = Label(), Spinner(), Tabs(), Window()
    for name in ("check_button", "save_button", "run_button", "cancel_button",
                 "import_case_button", "import_adat_button", "browse_button", "element_button", "undo_button"):
        setattr(w, name, Value())
    w.undo_button.enabled = False
    w.element_label, w.element_note = Label(), Label()
    w.choice_selection = Selection()
    w.errors, w.warnings = [], []     # a Toga dialog keeps no readable message: record _error and _warn

    async def record(title, message):
        w.errors.append((title, message))

    async def warn(title, message):
        w.warnings.append((title, message))

    w._error, w._warn = record, warn
    w._start = asyncio.run            # the handlers' async part, run to completion in the test
    return w, case


def test_every_form_field_has_a_widget():
    assert set(ri.FORM_FIELDS) == set(gen.NAME_FIELDS) | {k for _, keys in gen.TAB_FIELDS for k in keys}


class FolderDialog:
    def __init__(self, title, initial_directory=None):
        self.title, self.initial_directory = title, initial_directory


def test_import_library_case_opens_in_the_library(tmp_path, monkeypatch):
    w, _ = _window(tmp_path)
    monkeypatch.setattr(gen.toga, "SelectFolderDialog", FolderDialog)
    asyncio.run(w.on_import_case(None))
    assert w.main_window.dialogs[-1].initial_directory == rr.library_dir() == LIBRARY_DIR


def test_import_library_case_without_a_library_opens_the_default_folder(tmp_path, monkeypatch):
    w, _ = _window(tmp_path)
    monkeypatch.setattr(gen.toga, "SelectFolderDialog", FolderDialog)
    monkeypatch.setattr(gen.rr, "library_dir", lambda: tmp_path / "missing")
    asyncio.run(w.on_import_case(None))
    assert w.main_window.dialogs[-1].initial_directory is None


def test_collect_case_reads_the_widgets(tmp_path):
    w, case = _window(tmp_path)
    assert w._collect_case() == ri.case_from_fields(ri.case_to_fields(case))


def test_an_imported_case_keeps_its_tests_own_inputs_while_their_rows_are_unchanged(tmp_path):
    w, _ = _window(tmp_path)
    folder = LIBRARY_DIR / "020-Ca" / "020-Ca-gpw-nsp-campos"
    w._load_case(folder)
    imported = ri.read_case_folder(folder)
    assert "keep the library's own settings" in w.status_label.text
    assert w._collect_case().tests == imported.tests
    rows = w.inputs["tests"].value.splitlines()
    w.inputs["tests"].value = "\n".join(rows[:-1] + [rows[-1] + "  999=0"])     # last row changed
    tests = w._collect_case().tests
    assert tests[:-1] == imported.tests[:-1] and tests[-1].ae_input is None and tests[-1].test_input is None


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


def _fakes(monkeypatch, reform=True):
    monkeypatch.setenv("FAKE_REFORM_MODE", "ok")
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: "runatom.exe")
    monkeypatch.setattr(gen.rr, "runatom_command", lambda exe: FAKE)
    monkeypatch.setattr(gen.rr, "find_reform", lambda explicit=None: "reform.exe" if reform else None)
    monkeypatch.setattr(gen.rr, "reform_command", lambda exe: FAKE_REFORM)


def test_a_failed_rerun_shows_only_its_own_results(tmp_path, monkeypatch):
    _fakes(monkeypatch)
    w, _ = _window(tmp_path)
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "ok")
    asyncio.run(w.on_run(None))
    assert "Transferability test (s2d2): finished" in w.results_text.value
    formatted = tmp_path / "ti_ps_formatted.uspp"
    formatted.write_text("written by run 1")
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "fail_test")
    asyncio.run(w.on_run(None))
    text = w.results_text.value
    assert "Pseudopotential generation: finished" in text and "Stopped at Transferability test (s2d2)" in text
    assert "Transferability test (s2d2): finished" not in text
    assert "*ti_ps_formatted.uspp FORMATTED" in text
    assert formatted.read_bytes() == b"    7    3    6\nformatted\n"      # rewritten by run 2
    assert w._last_out == tmp_path / "ti_test_s2d2.out"
    ((title, message),) = w.errors
    assert title == "runatom stopped" and "could not be converged" in message
    assert w.status_label.text.startswith("Stopped at") and not w._busy


def test_a_run_that_fails_at_once_clears_the_old_results(tmp_path, monkeypatch):
    _fakes(monkeypatch)
    w, _ = _window(tmp_path)
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "ok")
    asyncio.run(w.on_run(None))
    assert w._last_out == tmp_path / "ti_ps.out"
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "error")
    asyncio.run(w.on_run(None))
    assert "FORMATTED" not in w.results_text.value and "finished" not in w.results_text.value
    assert "Stopped at All-electron (s2d1)" in w.results_text.value
    assert not (tmp_path / "ti_ps_formatted.uspp").exists() and not (tmp_path / "ti_ps.uspp").exists()
    assert w._last_out == tmp_path / "ti_ae_s2d1.out"


def test_busy_window_disables_imports_and_browse(tmp_path):
    w, _ = _window(tmp_path)
    w._set_busy(True, "Running")
    assert not any(b.enabled for b in (w.import_case_button, w.import_adat_button, w.browse_button))
    w._set_busy(False, "done")
    assert all(b.enabled for b in (w.import_case_button, w.import_adat_button, w.browse_button))


def test_a_missing_library_at_startup_is_reported(tmp_path):
    w, _ = _window(tmp_path)
    w._load_case(tmp_path / "no-such-case", quiet=True)
    assert "Library not found" in w.status_label.text


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
    monkeypatch.setattr(gen.rr, "run_case_job", lambda *a, **k: called.append(a))
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
    monkeypatch.setattr(gen.rr, "run_case_job", boom)
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


# --- element defaults ------------------------------------------------------
def _pick(w, symbol):
    w.element_controller.select_element(symbol)        # what a periodic-table click does


def test_closing_the_window_closes_the_periodic_table(tmp_path):
    w, _ = _window(tmp_path)
    w._on_close(None)
    assert w.element_controller.closed == 1


def test_picking_a_library_element_loads_its_pbe_case(tmp_path):
    w, _ = _window(tmp_path)
    _pick(w, "H")
    folder = LIBRARY_DIR / "001-H" / "001-H-gpbe--bm"
    assert w._current_fields() == ri.case_to_fields(ri.read_case_folder(folder))
    assert w._imported_case == ri.read_case_folder(folder)
    assert w.choice_selection.value.startswith("001-H-gpbe--bm")
    assert len(w.choice_selection.items) == 3 and w.element_label.text == "H"
    assert w.element_controller.shown[-1] == "H" and w.warnings == []
    assert "001-H-gpbe--bm" in w.element_note.text and w.undo_button.enabled


def test_a_library_element_without_pbe_warns(tmp_path):
    w, _ = _window(tmp_path)
    _pick(w, "Ti")
    ((title, message),) = w.warnings
    assert title == "Suggested values for Ti" and "No PBE case for Ti" in message


def test_picking_an_element_without_a_case_fills_a_template(tmp_path):
    w, _ = _window(tmp_path)
    _pick(w, "Mg")
    fields = w._current_fields()
    assert fields["ae_atom"] == "12  0  5" and fields["symbol"] == "Mg"
    assert fields["ps_rc"] == ri.ps_fields(ri.read_case_folder(LIBRARY_DIR / "020-Ca" / "020-Ca-gpw-nsp-campos").ps)["ps_rc"]
    assert w._imported_case is None
    ((title, message),) = w.warnings
    assert title == "Values copied from Ca" and "unchanged" in message
    assert w.choice_selection.items == ["Template from 020-Ca-gpw-nsp-campos", gen.rd.STARTER_LABEL]


def test_choosing_the_starter_blanks_the_radii(tmp_path):
    w, _ = _window(tmp_path)
    _pick(w, "Mg")
    w.choice_selection.value = gen.rd.STARTER_LABEL
    w.on_choice_change(w.choice_selection)
    assert w._current_fields()["ps_rc"] == "" and w._current_fields()["ae_atom"] == "12  0  5"
    assert len(w._undo_stack) == 2


def test_group_7_gets_the_starter_and_a_warning(tmp_path):
    w, _ = _window(tmp_path)
    _pick(w, "Mn")
    assert w.choice_selection.items == [gen.rd.STARTER_LABEL]
    assert "group 7" in w.warnings[0][1]
    asyncio.run(w.on_check(None))
    assert w.errors and "keyps" in w.errors[0][1]            # the blank field is named, nothing runs


def test_writing_the_choices_does_not_fill_twice(tmp_path, monkeypatch):
    w, _ = _window(tmp_path)
    loads = []
    real = w._load_case
    monkeypatch.setattr(w, "_load_case", lambda *a, **k: (loads.append(a), real(*a, **k)))
    _pick(w, "Cu")
    w.on_choice_change(w.choice_selection)                    # the echo of the window's own write
    assert len(loads) == 1


def test_an_edited_form_asks_before_it_is_replaced(tmp_path):
    w, _ = _window(tmp_path)
    w.inputs["ae_atom"].value = "22  0  5"                     # a hand edit
    w.main_window.answer = False
    _pick(w, "H")
    assert w.inputs["ae_atom"].value == "22  0  5" and w._undo_stack == []
    assert w.element_controller.shown == [] and len(w.main_window.dialogs) == 1
    w.main_window.answer = True
    _pick(w, "H")
    assert w.inputs["ae_atom"].value != "22  0  5" and len(w._undo_stack) == 1


def test_an_unedited_form_is_replaced_without_asking(tmp_path):
    w, _ = _window(tmp_path)
    _pick(w, "H")
    _pick(w, "Mg")
    assert w.main_window.dialogs == []


def test_a_declined_choice_puts_the_selection_back(tmp_path):
    w, _ = _window(tmp_path)
    _pick(w, "Mg")
    w.inputs["ps_rc"].value = "1.9"
    w.main_window.answer = False
    w.choice_selection.value = gen.rd.STARTER_LABEL
    w.on_choice_change(w.choice_selection)
    assert w.inputs["ps_rc"].value == "1.9"
    assert w.choice_selection.value == "Template from 020-Ca-gpw-nsp-campos"


def test_undo_defaults_restores_everything(tmp_path):
    w, _ = _window(tmp_path)
    start = w._current_fields()
    w.inputs["title"].value = "edited"                          # an edit that a confirmed pick replaced
    _pick(w, "H")
    _pick(w, "Mg")
    w.on_undo(None)
    assert w._current_fields()["symbol"] == ri.read_case_folder(LIBRARY_DIR / "001-H" / "001-H-gpbe--bm").symbol
    assert w._imported_case is not None and w.element_controller.shown[-1] == "H"
    assert w.choice_selection.value.startswith("001-H-gpbe--bm") and w.undo_button.enabled
    w.on_undo(None)
    assert w._current_fields() == dict(start, title="edited") and w._imported_case is None
    assert not w.undo_button.enabled and w._form_edited()      # the hand edit is an edit again
    w.on_undo(None)                                             # nothing left: no change, no error
    assert w._current_fields()["title"] == "edited"


def test_undo_keeps_at_most_ten_steps(tmp_path):
    w, _ = _window(tmp_path)
    for symbol in ["H", "Mg", "C", "Ge", "O", "Si"] * 2:
        _pick(w, symbol)
    assert len(w._undo_stack) == gen.UNDO_DEPTH


def test_fields_stay_editable_after_a_fill(tmp_path):
    w, _ = _window(tmp_path)
    _pick(w, "H")
    w.inputs["ae_atom"].value = "1  0  4"
    assert w._collect_case().ae.exfact == 4.0 and w._form_edited()


def test_without_a_library_every_element_gets_the_starter(tmp_path, monkeypatch):
    w, _ = _window(tmp_path)
    monkeypatch.setattr(gen.rr, "library_dirs", lambda: [tmp_path / "missing"])
    _pick(w, "Ti")
    assert w.choice_selection.items == [gen.rd.STARTER_LABEL]
    assert "library" in w.warnings[0][1].lower()


def test_elements_beyond_uranium_are_refused(tmp_path):
    w, _ = _window(tmp_path)
    _pick(w, "Np")
    assert w.errors and "92" in w.errors[0][1]
    assert "Np" not in gen.SELECTABLE_ELEMENTS and "U" in gen.SELECTABLE_ELEMENTS


def test_a_run_disables_the_element_row(tmp_path):
    w, _ = _window(tmp_path)
    _pick(w, "H")
    w._set_busy(True, "Running")
    assert not (w.element_button.enabled or w.choice_selection.enabled or w.undo_button.enabled)
    _pick(w, "Mg")                                              # ignored while busy
    assert w._current_fields()["symbol"] != "Mg"
    w._set_busy(False, "done")
    assert w.element_button.enabled and w.choice_selection.enabled and w.undo_button.enabled


def test_an_imported_case_shows_its_element(tmp_path):
    w, _ = _window(tmp_path)
    w._load_case(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-test")
    assert w.element_label.text == "Ti" and w.element_controller.shown[-1] == "Ti"
    assert w.choice_selection.value == "Imported: 000-Ex-ca-sp-test"
    assert w.choice_selection.items[-1] == "Imported: 000-Ex-ca-sp-test"
    w._load_case(LIBRARY_DIR / "020-Ca" / "020-Ca-gpw-nsp-campos")
    assert w.choice_selection.value.startswith("020-Ca-gpw-nsp-campos")
    assert not w._form_edited()
