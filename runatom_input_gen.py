"""Inputs > Vanderbilt USPP Generator: prepare, check and run runatom (USPP 7.3.6).

The window only moves text between widgets and a dict of field values.
runatomInput turns that dict into checked, column-exact input files,
runatomRunner runs runatom on them in a worker thread, and runatomOutput reads
the results.
"""

import asyncio
import os
import threading
import time
from pathlib import Path

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import runatomInput as ri
import runatomOutput as ro
import runatomRunner as rr
from help import HelpGqteaWin

STATUS_INTERVAL = 0.1          # s between two status-line updates during a run
LABEL_WIDTH = 240
HINT_COLOR = "#666666"
MULTILINE_HEIGHT = 110
DEFAULT_CASE = ("000-Ex", "000-Ex-ca-sp-test")
MAKEFILE_HINT = ("The Makefile assumes a USPP source tree: the job folder three levels below its root "
                 "(ROOT=../../..) and the program built as Bin/runatom.x.")
LIBRARY_MISSING = "Library not found: start from an empty form or import a case."
NAME_FIELDS = ("symbol", "generation", "title")
TAB_FIELDS = (
    ("All-electron", ("ae_print", "ae_logd", "ae_conv", "ae_atom", "ae_mesh", "ae_shells")),
    ("Pseudopotential", ("ps_print", "ps_logd", "ps_conv", "ps_partition", "ps_bessel", "ps_keyps",
                         "ps_rinner", "ps_rcloc", "ps_rc", "ps_betas", "ps_npf", "ps_local",
                         "ps_qfunc", "ps_fixes")),
    ("Tests", ("tests",)),
)


class RunatomGenerator:
    def __init__(self, *args, **kwargs):
        # gqteaWinToga passes the Button that was pressed; ignore it.
        self._closed, self._busy, self._cancel, self._last_out = False, False, None, None
        self._imported_case = None      # the library case last imported: its tests keep their own inputs
        self.inputs = {}
        self.main_window = toga.Window(title="Vanderbilt USPP Generator (runatom)", size=(920, 700))
        self.main_window.on_close = self._on_close
        self.layout_main_window()
        self.main_window.show()
        self._load_case(rr.library_dir().joinpath(*DEFAULT_CASE), quiet=True)

    # --- layout -----------------------------------------------------------
    def _row(self, label, *widgets):
        return toga.Box(style=Pack(direction=ROW, margin_bottom=2), children=[
            toga.Label(label, style=Pack(width=LABEL_WIDTH, margin=(4, 6, 0, 0))), *widgets])

    def _field_row(self, key):
        label, _, multiline, hint = ri.FORM_FIELDS[key]
        style = Pack(flex=1, height=MULTILINE_HEIGHT) if multiline else Pack(flex=1)
        self.inputs[key] = toga.MultilineTextInput(style=style) if multiline else toga.TextInput(style=style)
        box = toga.Box(style=Pack(direction=COLUMN, margin_bottom=4), children=[self._row(label, self.inputs[key])])
        if hint:
            box.add(toga.Label(hint, style=Pack(color=HINT_COLOR, font_size=9, margin_left=LABEL_WIDTH)))
        return box

    def layout_main_window(self):
        self.status_label = toga.Label("", style=Pack(flex=1, margin=(6, 6)))
        self.activity = toga.ActivityIndicator(style=Pack(margin=(6, 6)))
        self.folder_input = toga.TextInput(placeholder="folder for the input and output files", style=Pack(flex=1))
        self.exe_input = toga.TextInput(placeholder="blank = the runatom.exe shipped with gqteaWinToga",
                                        style=Pack(flex=1), on_change=self._update_exe_status)
        self.exe_status = toga.Label("", style=Pack(color=HINT_COLOR, font_size=9, margin_left=LABEL_WIDTH))
        self.makefile_switch = toga.Switch("Also write a Makefile (to run the case by hand under Linux/WSL)", value=False)
        self.makefile_hint = toga.Label(MAKEFILE_HINT, style=Pack(color=HINT_COLOR, font_size=9, margin_top=2))
        self.browse_button = toga.Button("Browse...", on_press=self.on_browse_folder)
        self.import_case_button = toga.Button("Import library case...", on_press=self.on_import_case)
        self.import_adat_button = toga.Button("Import .adat file...", on_press=self.on_import_adat)
        header = toga.Box(style=Pack(direction=COLUMN, margin=8), children=[
            self._row("Job folder", self.folder_input, self.browse_button),
            self._row("runatom program", self.exe_input), self.exe_status,
            toga.Box(style=Pack(direction=ROW, margin_top=4), children=[
                self.import_case_button, self.import_adat_button, self.makefile_switch]),
            self.makefile_hint])
        names = toga.Box(style=Pack(direction=COLUMN, margin=(0, 8)),
                         children=[self._field_row(k) for k in NAME_FIELDS])
        pages = []
        for title, keys in TAB_FIELDS:
            page = toga.Box(style=Pack(direction=COLUMN, margin=8), children=[self._field_row(k) for k in keys])
            pages.append((title, toga.ScrollContainer(content=page, horizontal=False)))
        self.results_text = toga.MultilineTextInput(readonly=True, style=Pack(flex=1, font_family="monospace"))
        pages.append(("Results", toga.Box(style=Pack(direction=COLUMN, margin=8), children=[self.results_text])))
        help_text = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_runatom, style=Pack(flex=1))
        pages.append(("Help", toga.Box(style=Pack(direction=COLUMN, margin=8), children=[help_text])))
        self.tabs = toga.OptionContainer(content=pages, style=Pack(flex=1))
        self.check_button = toga.Button("Check", on_press=self.on_check)
        self.save_button = toga.Button("Save inputs", on_press=self.on_save)
        self.run_button = toga.Button("Run", on_press=self.on_run)
        self.cancel_button = toga.Button("Cancel", on_press=self.on_cancel, enabled=False)
        bar = toga.Box(style=Pack(direction=ROW, margin=8), children=[
            self.check_button, self.save_button, self.run_button, self.cancel_button,
            toga.Button("Open folder", on_press=self.on_open_folder),
            toga.Button("Open .out", on_press=self.on_open_out), self.activity, self.status_label])
        self.main_window.content = toga.Box(style=Pack(direction=COLUMN),
                                            children=[header, names, self.tabs, bar])
        self._update_exe_status(None)

    # --- widgets <-> fields ------------------------------------------------
    def _collect_case(self) -> ri.CaseInputs:
        fields = {key: widget.value or "" for key, widget in self.inputs.items()}
        return ri.attach_imported_tests(ri.case_from_fields(fields), self._imported_case)

    def _fill(self, fields):
        for key, text in fields.items():
            widget = self.inputs.get(key)
            if widget is not None and widget.value != text:
                widget.value = text

    def _load_case(self, folder: Path, quiet: bool = False):
        try:
            case = ri.read_case_folder(folder)
        except (OSError, ValueError, KeyError) as exc:
            if not quiet:
                raise ValueError(f"{folder.name}: {exc}") from None
            self._show_status(LIBRARY_MISSING)
            return
        self._imported_case = case
        self._fill(ri.case_to_fields(case))
        skipped = ri.unrepresentable_tests(case)
        note = (f" Tests with their own shell list were not imported: {', '.join(skipped)}."
                if skipped else "")
        self._show_status(f"Loaded {folder.name}. Its tests keep the library's own settings "
                          f"while their row on the Tests tab is unchanged.{note}")

    def _update_exe_status(self, widget):
        exe = rr.find_runatom(self.exe_input.value)
        text = f"using {exe}" if exe else "runatom.exe not found: enter its path"
        if self.exe_status.text != text:
            self.exe_status.text = text

    # --- status, busy, progress --------------------------------------------
    def _show_status(self, message: str):
        if not self._closed and self.status_label.text != message:
            self.status_label.text = message

    def _set_busy(self, busy: bool, message: str):
        self._busy = busy
        if self._closed:
            return
        (self.activity.start if busy else self.activity.stop)()
        for button in (self.check_button, self.save_button, self.run_button,
                       self.import_case_button, self.import_adat_button, self.browse_button):
            button.enabled = not busy
        self.cancel_button.enabled = busy
        self._show_status(message)

    def _apply_progress(self, info: dict):
        self._show_status(rr.progress_text(info))

    def _progress_reporter(self, loop):
        """progress(info) for the worker thread: thinned, and handed to the window's loop."""
        last = [0.0]

        def report(info: dict):
            now = time.monotonic()
            if info["message"] != "starting" and now - last[0] < STATUS_INTERVAL:
                return
            last[0] = now
            if not self._closed:
                try:
                    loop.call_soon_threadsafe(self._apply_progress, info)
                except RuntimeError:     # the event loop has closed (application exit)
                    pass

        return report

    async def _error(self, title: str, message: str):
        if not self._closed:
            await self.main_window.dialog(toga.ErrorDialog(title, message))

    # --- actions ----------------------------------------------------------
    def _job_folder(self) -> Path:
        text = (self.folder_input.value or "").strip()
        if not text:
            raise ValueError("Choose the job folder first (Browse...).")
        return Path(text)

    def _prepared_case(self) -> ri.CaseInputs:
        case = self._collect_case()
        problems = ri.check_case(case)
        if problems:
            raise ValueError("\n".join(problems))
        return case

    def _write_inputs(self, case, folder: Path):
        written = ri.write_case(case, folder)
        if self.makefile_switch.value:
            (folder / "Makefile").write_text(ri.makefile_text(case), encoding="ascii", newline="\n")
            written.append(folder / "Makefile")
        return written

    async def on_check(self, widget):
        try:
            self._prepared_case()
        except ValueError as exc:
            await self._error("Check the inputs", str(exc))
            return
        self._show_status("The inputs pass every check runatom makes before it starts.")

    async def on_save(self, widget):
        try:
            case, folder = self._prepared_case(), self._job_folder()
            written = self._write_inputs(case, folder)
        except (ValueError, OSError) as exc:
            await self._error("Inputs not saved", str(exc))
            return
        self._show_status(f"Saved {len(written)} files in {folder}")

    async def on_run(self, widget):
        if self._busy:
            return
        try:
            case, folder = self._prepared_case(), self._job_folder()
            exe = rr.find_runatom(self.exe_input.value)
            if exe is None:
                raise ValueError("runatom.exe was not found. Enter its path, or reinstall gqteaWinToga, which ships it.")
            reform = rr.find_reform(self._reform_beside_runatom())
            self._write_inputs(case, folder)
        except (ValueError, OSError) as exc:
            await self._error("Cannot run runatom", str(exc))
            return
        stages = rr.plan_stages(case)
        self._last_out = None
        self._show_results("")
        self._cancel = threading.Event()
        self._set_busy(True, f"Running {len(stages)} stages...")
        loop = asyncio.get_running_loop()
        message = "Run ended"
        try:
            outcome = await asyncio.to_thread(
                rr.run_case_job, rr.runatom_command(exe), rr.reform_command(reform) if reform else None,
                stages, folder, rr.formatted_name(case), self._progress_reporter(loop), self._cancel)
            self._show_results(self._outcome_text(outcome, folder))
            seconds = sum(r.seconds for r in outcome.results)
            if outcome.error is None:
                self._last_out = folder / f"{ri.file_stems(case)['ps']}.out"
                message = f"Finished: {len(outcome.results)} stages in {seconds:.1f} s"
            else:
                self._last_out = folder / f"{outcome.failed_stage.stem}.out"
                message = f"Stopped at: {outcome.error.stage_label}"
                await self._error("runatom stopped", f"{outcome.error}\n\nThe .out file of that stage is in {folder}.")
        except rr.RunatomCancelled:
            message = "Run cancelled; the stage that was running left no output."
        except OSError as exc:
            message = "runatom could not be started"
            await self._error("runatom could not be started", str(exc))
        except Exception as exc:        # nothing may escape the handler or leave the window busy
            message = "The run failed unexpectedly"
            await self._error("runatom run failed unexpectedly", f"{type(exc).__name__}: {exc}")
        finally:
            self._set_busy(False, message)

    def _reform_beside_runatom(self):
        """reform.exe next to a typed runatom path, or None."""
        typed = (self.exe_input.value or "").strip()
        return str(Path(typed).with_name("reform.exe")) if typed and Path(typed).name else None

    @staticmethod
    def _outcome_text(outcome: rr.JobOutcome, folder: Path) -> str:
        """The Results tab: the completed stages, the CPMD file or why there is none, then the stop."""
        uspp, parts = None, []
        generated = next((r.stage for r in outcome.results if r.stage.kind == "pseudo"), None)
        if generated is not None:
            try:
                uspp = ro.read_uspp(folder / generated.uspp_file)
            except (OSError, ValueError) as exc:
                parts.append(f"The pseudopotential file could not be read back: {exc}")
        summary = ro.job_summary([(r.stage.label, r.seconds, r.report) for r in outcome.results], uspp)
        if outcome.formatted is not None:
            name = outcome.formatted.name
            parts.append(f"CPMD: {name} is the formatted pseudopotential. Copy it to the CPMD run directory "
                         "or CPMD's pseudopotential folder (CPMD_PP_LIBRARY_PATH or PP_LIBRARY_PATH, depending "
                         "on the CPMD version) and use the &ATOMS species line\n"
                         f"  {rr.cpmd_atom_line(name)}\nfollowed, as for every species, by the LMAX line "
                         "and the number of atoms (see the Help tab).")
        elif outcome.formatted_note:
            parts.append(outcome.formatted_note)
        if outcome.error is not None:
            parts.append(f"Stopped at {outcome.failed_stage.label}: {outcome.error.message}")
        return "\n\n".join(p for p in [summary] + parts if p)

    def _show_results(self, text: str):
        """The Results tab's text (written only on change); a non-empty text brings the tab forward."""
        if self._closed:
            return
        if self.results_text.value != text:
            self.results_text.value = text
        current = self.tabs.current_tab                  # a Toga OptionItem, or an index in the tests
        if text and getattr(current, "index", current) != len(TAB_FIELDS):
            self.tabs.current_tab = len(TAB_FIELDS)

    def on_cancel(self, widget):
        if self._cancel is not None:
            self._cancel.set()
            self._show_status("Cancelling...")

    async def on_browse_folder(self, widget):
        path = await self.main_window.dialog(toga.SelectFolderDialog("Job folder"))
        if path:
            self.folder_input.value = str(path)

    async def on_import_case(self, widget):
        path = await self.main_window.dialog(toga.SelectFolderDialog(
            "Library case folder (a folder with a Makefile)", initial_directory=rr.library_dir()))
        if path:
            try:
                self._load_case(Path(path))
            except ValueError as exc:
                await self._error("Case not imported", str(exc))

    async def on_import_adat(self, widget):
        path = await self.main_window.dialog(toga.OpenFileDialog("Input file", file_types=["adat"]))
        if not path:
            return
        try:
            data = ri.parse_adat(Path(path).read_text(encoding="latin-1"))
        except (OSError, ValueError) as exc:
            await self._error("File not imported", str(exc))
            return
        if isinstance(data, ri.AEInput):
            self._fill(ri.ae_fields(data))
        elif isinstance(data, ri.PSInput):
            self._fill(ri.ps_fields(data))
        else:
            await self._error("File not imported", "Test inputs are written from the Pseudopotential tab; "
                                                    "import the test's all-electron file and list it on the Tests tab.")
            return
        self._show_status(f"Imported {Path(path).name}")

    def on_open_folder(self, widget):
        text = (self.folder_input.value or "").strip()
        if text and hasattr(os, "startfile") and Path(text).is_dir():
            os.startfile(text)

    def on_open_out(self, widget):
        if self._last_out is not None and self._last_out.is_file() and hasattr(os, "startfile"):
            os.startfile(str(self._last_out))

    def _on_close(self, window, **kwargs):
        self._closed = True        # a run still going finishes without touching the window
        if self._cancel is not None:
            self._cancel.set()
        return True

    def closeTopLevel(self, widget):
        self._on_close(self.main_window)
        self.main_window.close()
