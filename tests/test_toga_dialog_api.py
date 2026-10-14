"""Deprecated Toga 0.5 dialog calls stay out of the converted modules.

``Window.info_dialog`` / ``error_dialog`` / ``open_file_dialog`` (and friends)
print a DeprecationWarning; the codebase uses
``await window.dialog(toga.XxxDialog(...))`` instead. ``Window.dialog`` is a
coroutine, so a call site in sync code must be async or be scheduled -- a bare
rename would build a coroutine that never runs and the dialog would silently
never appear. These tests drive the async paths with a recording window.
"""

import asyncio
import importlib
import inspect
import re

import pytest
import toga

DEPRECATED = re.compile(
    r"main_window\.(info|error|question|confirm|stack_trace|"
    r"open_file|save_file|select_folder)_dialog\("
)


@pytest.mark.parametrize(
    "module_name", ["coordinateConverter", "cpxForces", "cpTrajec2xyz", "plotter"]
)
def test_no_deprecated_dialog_calls(module_name):
    source = inspect.getsource(importlib.import_module(module_name))
    assert DEPRECATED.findall(source) == []


class _RecordingWindow:
    def __init__(self, open_result=None):
        self.dialogs = []
        self.open_result = open_result

    async def dialog(self, dialog):
        self.dialogs.append(dialog)
        if isinstance(dialog, toga.OpenFileDialog):
            return self.open_result
        return None


class _StubInput:
    def __init__(self, value=""):
        self.value = value
        self.text = ""


@pytest.mark.parametrize(
    "module_name, class_name, path_attr",
    [
        ("cpxForces", "ForceXYZtoForce", "cp_force_textInput"),
        ("cpTrajec2xyz", "TrajectoryConverter", "cp_trajec_textInput"),
    ],
)
@pytest.mark.filterwarnings("error::DeprecationWarning")
def test_conversion_failure_shows_error_dialog(tmp_path, module_name, class_name, path_attr):
    cls = getattr(importlib.import_module(module_name), class_name)
    tool = cls.__new__(cls)
    cls.__init__(tool)
    tool.main_window = _RecordingWindow()
    tool.input_textInput = _StubInput(str(tmp_path / "missing.in"))
    setattr(tool, path_attr, _StubInput(str(tmp_path / "missing.traj")))
    tool.output_file_label = _StubInput()

    assert inspect.iscoroutinefunction(tool.perform_conversion)
    asyncio.run(tool.perform_conversion(None))

    assert len(tool.main_window.dialogs) == 1
    assert isinstance(tool.main_window.dialogs[0], toga.ErrorDialog)


@pytest.mark.parametrize(
    "module_name, class_name",
    [("cpxForces", "ForceXYZtoForce"), ("cpTrajec2xyz", "TrajectoryConverter")],
)
@pytest.mark.filterwarnings("error::DeprecationWarning")
def test_select_file_uses_open_file_dialog(tmp_path, module_name, class_name):
    cls = getattr(importlib.import_module(module_name), class_name)
    tool = cls.__new__(cls)
    cls.__init__(tool)
    chosen = tmp_path / "cp.in"
    tool.main_window = _RecordingWindow(open_result=chosen)
    target = _StubInput()

    result = asyncio.run(tool.select_file(None, "Pick", target))

    assert isinstance(tool.main_window.dialogs[0], toga.OpenFileDialog)
    assert result == str(chosen) and target.value == str(chosen)

    tool.main_window = _RecordingWindow(open_result=None)
    assert asyncio.run(tool.select_file(None, "Pick", _StubInput())) is None
    assert isinstance(tool.main_window.dialogs[-1], toga.InfoDialog)


@pytest.mark.filterwarnings("error::DeprecationWarning")
def test_plotter_sync_helper_schedules_info_dialog():
    from plotter import PlotterBase

    plotter = PlotterBase.__new__(PlotterBase)
    plotter.main_window = _RecordingWindow()
    plotter.data = []
    plotter.x_axis = [1.0]
    plotter.unit_selection = _StubInput("parsecs")

    async def run():
        plotter.compute_x_axis(1.0)  # sync code path, unsupported unit
        await asyncio.sleep(0)       # let the scheduled dialog run
        await asyncio.sleep(0)

    asyncio.run(run())
    assert plotter.x_axis == []
    assert len(plotter.main_window.dialogs) == 1
    assert isinstance(plotter.main_window.dialogs[0], toga.InfoDialog)
