"""Tests for toolLoader (the main-menu catalogue and the on-demand tool loader)
and for the main window's use of it in gqteaWinToga."""

import ast
import asyncio
import inspect
import subprocess
import sys
import threading
from pathlib import Path

import pytest

from toolLoader import (
    PREREQUISITE_MODULES,
    TOOL_TABS,
    ToolLoader,
    ToolSpec,
    all_tools,
    tool_modules,
)

SRC = Path(__file__).resolve().parents[1]


# -- the catalogue ----------------------------------------------------------

def test_tabs_are_the_five_menu_tabs_in_order():
    assert [title for title, _ in TOOL_TABS] == ["Geometry", "Inputs", "Structural", "Thermo", "Tools"]


def test_labels_are_unique_and_every_tab_has_tools():
    labels = [spec.label for spec in all_tools()]
    assert len(labels) == len(set(labels))
    assert all(specs for _, specs in TOOL_TABS)


def test_every_entry_opens_a_class_that_accepts_the_button():
    """Toga calls the handler with the Button; the tool class gets it as before."""
    for spec in all_tools():
        module = __import__(spec.module)
        tool = getattr(module, spec.attribute)
        assert callable(tool), spec
        inspect.signature(tool).bind(object())  # raises TypeError if it cannot


def test_frozen_build_imports_list_every_tool_and_prerequisite():
    """PyInstaller only bundles modules it sees imported literally."""
    tree = ast.parse((SRC / "toolLoader.py").read_text(encoding="utf-8"))
    function = next(node for node in ast.walk(tree)
                    if isinstance(node, ast.FunctionDef) and node.name == "_frozen_build_imports")
    imported = {alias.name for node in ast.walk(function) if isinstance(node, ast.Import)
                for alias in node.names}
    assert imported == set(tool_modules()) | set(PREREQUISITE_MODULES)


def _modules_after_import(statement):
    code = (f"import sys; sys.path.insert(0, {str(SRC)!r}); {statement}; "
            "print('\\n'.join(sorted(sys.modules)))")
    output = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                            check=True, cwd=str(SRC)).stdout
    return set(output.split())


def test_importing_the_entry_point_imports_no_tool():
    loaded = _modules_after_import("import gqteaWinToga")
    assert not loaded & set(tool_modules())
    assert not loaded & {"numba", "OpenGL", "scipy", "matplotlib", "displayPlots"}


def test_tool_loader_does_not_import_toga():
    assert "toga" not in _modules_after_import("import toolLoader")


# -- ToolLoader -------------------------------------------------------------

class FakeImporter:
    """Records the order of first imports; ``gate`` holds every import until released."""

    def __init__(self, failing=(), gate=None):
        self.order = []
        self.failing = set(failing)
        self.gate = gate
        self.started = threading.Event()

    def __call__(self, name):
        if name not in self.order:
            self.order.append(name)
        self.started.set()
        if self.gate is not None:
            self.gate.wait(5)
        if name in self.failing:
            raise ImportError(f"no module {name}")
        return f"module:{name}"


def _wait(loader):
    for _ in range(500):
        if not loader.is_busy():
            return
        threading.Event().wait(0.01)
    raise AssertionError("loader did not finish")


def test_preload_imports_prerequisites_first_then_tools_in_order():
    importer = FakeImporter()
    loader = ToolLoader(["a", "b", "c"], prerequisites=["pre"], import_module=importer)
    loader.start_preload()
    _wait(loader)
    assert importer.order == ["pre", "a", "b", "c"]
    assert loader.loaded_count() == loader.total_count() == 3
    assert loader.request("b").result() == "module:b"


def test_a_requested_tool_jumps_the_preload_queue():
    gate = threading.Event()
    importer = FakeImporter(gate=gate)
    loader = ToolLoader(["a", "b", "c", "d"], prerequisites=[], import_module=importer)
    loader.start_preload()
    assert importer.started.wait(5)  # "a" is being imported
    future = loader.request("d")
    gate.set()
    assert future.result(5) == "module:d"
    _wait(loader)
    assert importer.order == ["a", "d", "b", "c"]


def test_a_failed_import_is_reported_and_retried_on_the_next_request():
    importer = FakeImporter(failing={"b"})
    reports = []
    loader = ToolLoader(["a", "b"], prerequisites=[], import_module=importer,
                        on_loaded=lambda name, error: reports.append((name, type(error).__name__)))
    loader.start_preload()
    _wait(loader)
    assert reports == [("a", "NoneType"), ("b", "ImportError")]
    assert (loader.loaded_count(), loader.failed_count()) == (1, 1)
    with pytest.raises(ImportError):
        loader.request("b").result(5)  # retried, fails again
    importer.failing.clear()
    assert loader.request("b").result(5) == "module:b"
    assert (loader.loaded_count(), loader.failed_count()) == (2, 0)


def test_on_loaded_is_not_called_for_prerequisites_and_cannot_stop_the_loader():
    importer = FakeImporter()
    calls = []

    def explode(name, error):
        calls.append(name)
        raise RuntimeError("reporting failed")

    loader = ToolLoader(["a", "b"], prerequisites=["pre"], import_module=importer, on_loaded=explode)
    loader.start_preload()
    _wait(loader)
    assert calls == ["a", "b"] and loader.loaded_count() == 2


def test_stop_drops_the_queue_but_a_request_still_works():
    gate = threading.Event()
    importer = FakeImporter(gate=gate)
    loader = ToolLoader(["a", "b", "c"], prerequisites=[], import_module=importer)
    loader.start_preload()
    assert importer.started.wait(5)
    loader.stop()
    gate.set()
    _wait(loader)
    assert importer.order == ["a"]
    assert loader.request("c").result(5) == "module:c"


def test_the_real_prerequisite_pins_the_agg_backend():
    assert PREREQUISITE_MODULES[0] == "displayPlots"


# -- the main window --------------------------------------------------------

import gqteaWinToga as gq  # noqa: E402


class FakeLabel:
    def __init__(self):
        self._text = ""
        self.writes = 0

    @property
    def text(self):
        return self._text

    @text.setter
    def text(self, value):
        self._text = value
        self.writes += 1


class FakeIndicator:
    def __init__(self):
        self.is_running = False

    def start(self):
        self.is_running = True

    def stop(self):
        self.is_running = False


class FakeWindow:
    def __init__(self):
        self.dialogs = []

    async def dialog(self, dialog):
        self.dialogs.append(dialog)


class FakeDialog:
    def __init__(self, title, message):
        self.title, self.message = title, message


class FakeButton:
    def __init__(self):
        self.enabled = True
        self.enabled_history = []

    def __setattr__(self, name, value):
        if name == "enabled" and hasattr(self, "enabled_history"):
            self.enabled_history.append(value)
        object.__setattr__(self, name, value)


class OpenedTool:
    instances = []

    def __init__(self, *args, **kwargs):
        OpenedTool.instances.append((args, kwargs))


def _app(loader):
    app = gq.gqteaWin.__new__(gq.gqteaWin)
    app.tool_loader = loader
    app._opening = []
    app.status_label = FakeLabel()
    app.loading_indicator = FakeIndicator()
    app._main_window = FakeWindow()  # behind toga.App.main_window, which validates
    return app


def _module_with_tool():
    module = type(sys)("fake_tool_module")
    module.OpenedTool = OpenedTool
    return module


def test_open_tool_loads_the_module_then_calls_the_class_with_the_button():
    OpenedTool.instances.clear()
    loader = ToolLoader(["fake"], prerequisites=[], import_module=lambda name: _module_with_tool())
    app = _app(loader)
    button = FakeButton()
    asyncio.run(app.open_tool(ToolSpec("Fake tool", "fake", "OpenedTool"), button))
    assert OpenedTool.instances == [((button,), {})]
    assert button.enabled_history == [False, True] and app._opening == []
    assert app.status_label.text == "All 1 analysis tools are ready."
    assert not app.loading_indicator.is_running


def test_open_tool_on_a_loaded_module_does_not_touch_the_button():
    OpenedTool.instances.clear()
    loader = ToolLoader(["fake"], prerequisites=[], import_module=lambda name: _module_with_tool())
    loader.request("fake").result(5)
    app = _app(loader)
    button = FakeButton()
    asyncio.run(app.tool_handler(ToolSpec("Fake tool", "fake", "OpenedTool"))(button))
    assert OpenedTool.instances == [((button,), {})]
    assert button.enabled_history == []


def test_open_tool_reports_an_import_failure_in_a_dialog(monkeypatch, capsys):
    monkeypatch.setattr(gq.toga, "ErrorDialog", FakeDialog)
    OpenedTool.instances.clear()

    def fail(name):
        raise ImportError("No module named 'glfw'")

    app = _app(ToolLoader(["fake"], prerequisites=[], import_module=fail))
    button = FakeButton()
    asyncio.run(app.open_tool(ToolSpec("Fake tool", "fake", "OpenedTool"), button))
    assert OpenedTool.instances == []
    [dialog] = app.main_window.dialogs
    assert dialog.title == "Tool could not be loaded"
    assert "Fake tool" in dialog.message and "No module named 'glfw'" in dialog.message
    assert button.enabled and app._opening == []
    assert "could not be loaded" in app.status_label.text
    assert "ImportError" in capsys.readouterr().err


def test_load_status_texts():
    loader = ToolLoader(["a", "b"], prerequisites=[], import_module=lambda name: name)
    app = _app(loader)
    assert app.load_status() == ("Loading analysis tools in the background: 0 of 2 ready. "
                                 "You can already open any tool.", True)
    app._opening.append("Plots")
    assert app.load_status() == ("Opening Plots...", True)
    app._opening.clear()
    loader.start_preload()
    _wait(loader)
    assert app.load_status() == ("All 2 analysis tools are ready.", False)


def test_update_load_status_writes_the_label_only_when_it_changes():
    app = _app(ToolLoader(["a"], prerequisites=[], import_module=lambda name: name))
    app.update_load_status()
    app.update_load_status()
    assert app.status_label.writes == 1 and app.loading_indicator.is_running


def test_a_report_after_the_loop_closed_is_ignored(monkeypatch):
    app = _app(ToolLoader(["a"], prerequisites=[]))
    loop = asyncio.new_event_loop()
    loop.close()
    monkeypatch.setattr(gq.gqteaWin, "loop", property(lambda self: loop))
    app._tool_loaded_from_thread("a", None)  # must not raise
