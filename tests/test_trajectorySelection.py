"""Atom/molecule selection in the All bond / angle / dihedral analysis windows (headless).

The windows are built with ``__new__`` + stub widgets and a SelectionViewer
whose GL window is never opened; the Solute atom indices field is the
selection field (trajectorySelection.TrajectorySelectionUI).
"""
import asyncio
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import toga

import trajectorySelection as ts
import viewerSelectionControls as vsc
from allAnglesAnalysis import AllAnglesAnalysis, allAnglesAnalysisUI
from allBondAnalysis import AllBondAnalysis, allBondAnalysisUI
from allDihedralAnalysis import AllDihedralAnalysis, allDihedralAnalysisUI

SRC = Path(__file__).resolve().parents[1]

# Two waters: O1 H2 H3 and O4 H5 H6 (1-based), bonded within each molecule only.
DIMER = [
    ("O", (0.0, 0.0, 0.0)), ("H", (0.96, 0.0, 0.0)), ("H", (-0.24, 0.93, 0.0)),
    ("O", (2.9, 0.0, 0.0)), ("H", (3.4, 0.8, 0.0)), ("H", (3.4, -0.8, 0.0)),
]
# Hydrogen peroxide H-O-O-H (atoms 1-4, one dihedral) beside a water (atoms 5-7).
PEROXIDE_AND_WATER = [
    ("O", (0.0, 0.0, 0.0)), ("O", (1.45, 0.0, 0.0)), ("H", (-0.3, 0.9, 0.1)), ("H", (1.75, -0.2, 0.9)),
    ("O", (6.0, 0.0, 0.0)), ("H", (6.96, 0.0, 0.0)), ("H", (5.76, 0.93, 0.0)),
]
TOOLS = [
    (allBondAnalysisUI, AllBondAnalysis),
    (allAnglesAnalysisUI, AllAnglesAnalysis),
    (allDihedralAnalysisUI, AllDihedralAnalysis),
]


def write_dimer(path, frames=2, atoms=DIMER):
    lines = []
    for _ in range(frames):
        lines += [str(len(atoms)), "test structure"]
        lines += [f"{e} {x:.6f} {y:.6f} {z:.6f}" for e, (x, y, z) in atoms]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


class _StubLoop:
    def is_closed(self):
        return False

    def call_soon_threadsafe(self, callback, *args):
        callback(*args)


class _DialogWindow:
    def __init__(self):
        self.dialogs = []
        self.closed = False

    async def dialog(self, dialog):
        self.dialogs.append(dialog)

    def close(self):
        self.closed = True


def _message(dialog):
    return dialog._impl.message


def make_ui(ui_class, tmp_path, solute=""):
    """A tool window without Toga widgets; the solute field is the selection field."""
    ui = ui_class.__new__(ui_class)
    ui.trajec = write_dimer(tmp_path / "TRAJEC.xyz")
    ui.output_dir = str(tmp_path)
    ui.main_window = _DialogWindow()
    ui.textInput_file = SimpleNamespace(value=ui.trajec)
    ui.textInput_max_distance = SimpleNamespace(value="")
    ui.textInput_cell = SimpleNamespace(value="")
    ui.textInput_output = SimpleNamespace(value="combined.txt")
    ui.textInput_solute_indices = SimpleNamespace(value=solute)
    ui.selection_input = ui.textInput_solute_indices
    ui.selection_summary_label = SimpleNamespace(text="No atoms selected.")
    ui.selection_status_label = SimpleNamespace(text="")
    ui.show_selected_switch = SimpleNamespace(value=False)
    ui.multi_line_text = SimpleNamespace(value="")
    return ui


def load(ui, analysis_class):
    """What Browse does after reading the first frame."""
    elements, coordinates = analysis_class(ui.trajec).read_first_frame()
    ui.load_selection_structure(elements, coordinates)
    ui.structure_viewer.loop = _StubLoop()  # a Ctrl/Shift-click is marshalled onto the Toga loop
    return ui.structure_viewer


@pytest.fixture(params=TOOLS, ids=["bond", "angle", "dihedral"])
def tool(request):
    return request.param


# ----------------------------------------------------------------------
# The shared controls accept an existing field
# ----------------------------------------------------------------------
class _FakeWidget:
    def __init__(self, *args, **kwargs):
        self.args, self.kwargs = args, kwargs
        self.on_confirm = kwargs.get("on_confirm")


@pytest.fixture
def fake_toga(monkeypatch):
    monkeypatch.setattr(vsc, "toga", SimpleNamespace(TextInput=_FakeWidget, Label=_FakeWidget,
                                                     Button=_FakeWidget, Switch=_FakeWidget))
    monkeypatch.setattr(vsc, "Pack", lambda **kwargs: kwargs)


class _Host(vsc.SelectionControls):
    def selection_viewer(self):
        return None


def test_default_controls_create_their_own_field(fake_toga):
    host = _Host()
    host.build_selection_widgets()
    assert isinstance(host.selection_input, _FakeWidget)
    assert host.selection_input.on_confirm == host.apply_selection_field


def test_controls_can_reuse_an_existing_field(fake_toga):
    host = _Host()
    field = SimpleNamespace(value="1-3", on_confirm=None)
    host.build_selection_widgets(selection_input=field)
    assert host.selection_input is field and field.value == "1-3"
    assert field.on_confirm == host.apply_selection_field


def test_cube_tab_assign_is_the_shared_writer():
    import cubeTab
    import widgetWrites
    assert cubeTab._assign is widgetWrites.assign_if_changed
    assert vsc._assign is widgetWrites.assign_if_changed


# ----------------------------------------------------------------------
# Selection <-> Solute atom indices
# ----------------------------------------------------------------------
def test_clicks_fill_the_solute_field(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path)
    viewer = load(ui, tool[1])
    viewer.toggle_selected_atom(4)                      # Ctrl+click on H5
    assert ui.textInput_solute_indices.value == "5"
    assert ui.selection_summary_label.text == "1 atom: H"
    viewer.toggle_selected_molecule(0)                  # Shift+click on O1
    assert ui.textInput_solute_indices.value == "1-3, 5"
    assert ui.selection_status_label.text.startswith("Selection: 4 atoms")


def test_typed_field_and_enter_select_and_normalise(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path)
    viewer = load(ui, tool[1])
    ui.textInput_solute_indices.value = "6 4,5"
    asyncio.run(ui.apply_selection_field(None))
    assert viewer.get_selected_atoms() == [3, 4, 5]
    assert ui.textInput_solute_indices.value == "4-6"


def test_expand_and_clear(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path)
    viewer = load(ui, tool[1])
    viewer.set_selected_atoms([4])
    ui.expand_selection(None)
    assert ui.textInput_solute_indices.value == "4-6"
    ui.clear_atom_selection(None)
    assert ui.textInput_solute_indices.value == "" and viewer.get_selected_atoms() == []
    assert ui.selection_summary_label.text == "No atoms selected."


def test_browse_keeps_the_typed_indices_as_the_selection(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path, solute="1-3")
    viewer = load(ui, tool[1])
    assert viewer.get_selected_atoms() == [0, 1, 2]
    assert ui.textInput_solute_indices.value == "1-3"
    viewer = load(ui, tool[1])                          # a second Browse reuses the viewer
    assert viewer.get_selected_atoms() == [0, 1, 2] and ui.textInput_solute_indices.value == "1-3"


def test_browse_keeps_indices_the_new_file_does_not_have(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path, solute="1-3 9")
    viewer = load(ui, tool[1])
    assert viewer.get_selected_atoms() == []
    assert ui.textInput_solute_indices.value == "1-3 9"   # left as typed
    assert "Atom 9" in ui.selection_status_label.text
    assert ui.selection_summary_label.text == "No atoms selected."


def test_enter_before_a_trajectory_is_chosen(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path, solute="1-3")
    asyncio.run(ui.apply_selection_field(None))
    assert isinstance(ui.main_window.dialogs[-1], toga.ErrorDialog)
    assert _message(ui.main_window.dialogs[-1]) == vsc.NOTHING_LOADED_MESSAGE
    assert ui.textInput_solute_indices.value == "1-3"


def test_a_3d_failure_never_breaks_browse(tool, tmp_path, monkeypatch):
    ui = make_ui(tool[0], tmp_path, solute="2")

    def broken():
        raise RuntimeError("glfw3.dll not found")
    monkeypatch.setattr(ts, "_selection_viewer_class", broken)
    elements, coordinates = tool[1](ui.trajec).read_first_frame()
    ui.load_selection_structure(elements, coordinates)
    assert ui.structure_viewer is None and ui.textInput_solute_indices.value == "2"
    assert "glfw3.dll" in ui.selection_status_label.text


# ----------------------------------------------------------------------
# The analysis reads the field, as before
# ----------------------------------------------------------------------
def test_read_params_uses_the_clicked_selection(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path)
    viewer = load(ui, tool[1])
    viewer.toggle_selected_molecule(4)
    assert asyncio.run(ui.read_params()) is True
    assert ui.solute_atom_indices == [3, 4, 5]


def test_typed_text_without_enter_or_3d_still_works(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path, solute="4-6")
    assert asyncio.run(ui.read_params()) is True
    assert ui.solute_atom_indices == [3, 4, 5]
    ui.textInput_solute_indices.value = ""
    assert asyncio.run(ui.read_params()) is True and ui.solute_atom_indices is None


def test_workflow_writes_the_selected_scope(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path)
    ui.trajec = write_dimer(tmp_path / "peroxide.xyz", atoms=PEROXIDE_AND_WATER)
    viewer = load(ui, tool[1])
    viewer.toggle_selected_molecule(0)                  # Shift+click on the peroxide
    assert ui.textInput_solute_indices.value == "1-4"
    asyncio.run(ui.workflow(None))
    assert ui.main_window.dialogs == []
    text = (tmp_path / "combined.txt").read_text()
    assert "# atom_scope solute_atoms" in text and "# solute_atom_indices 1 2 3 4" in text


# ----------------------------------------------------------------------
# Show 3D and closing
# ----------------------------------------------------------------------
def test_show_3d_needs_a_trajectory(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path)
    asyncio.run(ui.show_structure(None))
    assert _message(ui.main_window.dialogs[-1]) == ts.NO_TRAJECTORY_MESSAGE


def test_show_3d_reports_a_window_that_cannot_open(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path)
    viewer = load(ui, tool[1])

    async def refuse(timeout=5.0):
        return False, "Only one 3D window at a time."
    viewer.ensure_viewer_window = refuse
    asyncio.run(ui.show_structure(None))
    assert "Only one 3D window" in _message(ui.main_window.dialogs[-1])


def test_closing_the_tool_closes_its_3d_window(tool, tmp_path):
    ui = make_ui(tool[0], tmp_path)
    ui.closeTopLevel(None)                              # no viewer yet: nothing to close
    assert ui.main_window.closed
    ui = make_ui(tool[0], tmp_path)
    viewer = load(ui, tool[1])
    requests = []
    viewer.request_close = lambda: requests.append(True)
    assert ui.on_window_close(ui.main_window) is True
    ui.closeTopLevel(None)
    assert requests == [True, True]


# ----------------------------------------------------------------------
# Viewer and imports
# ----------------------------------------------------------------------
def test_selection_changes_reach_the_listener_on_the_loop():
    from selectionViewer import SelectionViewer
    viewer = SelectionViewer()
    calls = []
    viewer.toggle_selected_atom(0)                      # no listener, no loop: nothing happens
    viewer.selection_listener = lambda: calls.append(viewer.get_selected_atoms())
    viewer.loop = _StubLoop()
    viewer.toggle_selected_atom(1)
    assert calls == [[0, 1]]


def test_the_fukui_viewer_is_a_selection_viewer():
    from espViewer import ESPSurfaceViewer
    from selectionViewer import SelectionViewer
    assert issubclass(ESPSurfaceViewer, SelectionViewer)
    assert ESPSurfaceViewer().selection_listener is None


def test_the_tools_do_not_import_the_3d_viewer():
    code = ("import sys, allBondAnalysis, allAnglesAnalysis, allDihedralAnalysis; "
            "print(sorted(m for m in ('molecularViewer', 'selectionViewer', 'OpenGL', 'numba', 'cubeTab') "
            "if m in sys.modules))")
    out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "[]"
