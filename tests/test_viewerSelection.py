"""Atom/molecule selection in the 3D Molecular Viewer (headless).

Pure helpers in viewerSelection/indexRanges first, then the viewer's
selection state and click routing, rendering filters, and the UI plumbing
(MolecularViewerUI built with __new__ + stub widgets; no window is created).
"""
import subprocess
import sys
from pathlib import Path

import pytest

import viewerSelection as vs
from indexRanges import format_index_ranges, parse_solute_index_ranges

SRC = Path(__file__).resolve().parents[1]

# Two waters: O0 H1 H2 and O3 H4 H5 (bonded within each molecule only).
DIMER_BONDS = [(0, 1), (0, 2), (3, 4), (3, 5)]
DIMER_ELEMENTS = ["O", "H", "H", "O", "H", "H"]


# ----------------------------------------------------------------------
# Pure helpers
# ----------------------------------------------------------------------
def test_format_index_ranges_compresses_runs():
    assert format_index_ranges([1, 2, 3, 7, 9, 10]) == "1-3, 7, 9-10"
    assert format_index_ranges([5, 3, 4, 4]) == "3-5"
    assert format_index_ranges([]) == ""
    assert format_index_ranges([295, 296, 297, 1], sep=",") == "1,295-297"


@pytest.mark.parametrize(
    "values", [[1], [1, 3, 5], list(range(1, 50)), [2, 3, 10, 11, 12, 40]]
)
def test_format_round_trips_through_the_parser(values):
    assert parse_solute_index_ranges(format_index_ranges(values)) == values


def test_connected_atoms_returns_the_seed_molecule():
    assert vs.connected_atoms(4, DIMER_BONDS, 6) == [3, 4, 5]
    assert vs.connected_atoms(0, DIMER_BONDS, 6) == [0, 1, 2]


def test_connected_atoms_of_an_unbonded_atom_is_itself():
    assert vs.connected_atoms(2, [(0, 1)], 4) == [2]
    assert vs.connected_atoms(0, [], 1) == [0]


def test_connected_atoms_rejects_an_index_outside_the_frame():
    with pytest.raises(IndexError):
        vs.connected_atoms(6, DIMER_BONDS, 6)


def test_expand_to_molecules_adds_every_partner():
    assert vs.expand_to_molecules([1, 6], DIMER_BONDS, 7) == [0, 1, 2, 6]
    assert vs.expand_to_molecules([], DIMER_BONDS, 6) == []


def test_visible_atoms_without_a_selection_filter_is_everything():
    bonds = [(0, 1), (1, 2)]
    for selection in (None, set()):
        atoms, shown = vs.visible_atoms_and_bonds(3, bonds, selection)
        assert list(atoms) == [0, 1, 2]
        assert shown is bonds  # no copy on the common path


def test_visible_atoms_keep_only_selected_atoms_and_inner_bonds():
    atoms, shown = vs.visible_atoms_and_bonds(4, [(0, 1), (1, 2), (2, 3)], {1, 2, 9})
    assert atoms == [1, 2]  # 9 does not exist in this frame
    assert shown == [(1, 2)]


def test_selection_summary():
    assert vs.selection_summary(DIMER_ELEMENTS, []) == "No atoms selected."
    assert vs.selection_summary(DIMER_ELEMENTS, [3]) == "1 atom: O"
    assert vs.selection_summary(DIMER_ELEMENTS, list(range(6))) == "6 atoms: H4 O2"
    assert vs.selection_summary(DIMER_ELEMENTS[:2], [0, 5]) == "2 atoms: O (1 not in this frame)"


def test_viewer_selection_imports_neither_toga_nor_opengl():
    code = "import sys, viewerSelection; print('toga' in sys.modules, 'OpenGL' in sys.modules)"
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True
    )
    assert result.stdout.split() == ["False", "False"]


# ----------------------------------------------------------------------
# XYZ export helpers
# ----------------------------------------------------------------------
FRAMES = [
    [("O", (0.0, 0.0, 0.0)), ("H", (0.96, 0.0, 0.0)), ("H", (-0.24, 0.93, 0.0)), ("O", (2.9, 0.0, 0.0))],
    [("O", (0.1, 0.0, 0.0)), ("H", (1.0, 0.1, 0.0)), ("H", (-0.2, 0.9, 0.1)), ("O", (2.8, 0.1, 0.0))],
]


def _old_viewer_frame_text(frame, comment):
    """What save_current_frame_xyz wrote before selection existed."""
    lines = [f"{len(frame)}\n", f"{comment}\n"]
    lines += [f"{element:<3s} {x: .8f} {y: .8f} {z: .8f}\n" for element, (x, y, z) in frame]
    return "".join(lines)


def test_xyz_frame_text_matches_the_viewer_format():
    comment = "Frame 0 exported from molecularViewer"
    assert vs.xyz_frame_text(FRAMES[0], None, comment) == _old_viewer_frame_text(FRAMES[0], comment)


def test_selected_frame_keeps_the_given_source_order():
    lines = vs.xyz_frame_text(FRAMES[1], [0, 3], "c").splitlines()
    assert lines[0] == "2" and lines[1] == "c"
    assert lines[2].split() == ["O", "0.10000000", "0.00000000", "0.00000000"]
    assert lines[3].split() == ["O", "2.80000000", "0.10000000", "0.00000000"]


def test_comment_is_forced_onto_one_line():
    assert vs.xyz_frame_text(FRAMES[0], [0], "a\nb\r\n").splitlines()[1] == "a b"


def test_selection_comment_names_frame_source_and_atoms():
    assert vs.selection_comment(7, "TRAJEC.xyz", [0, 1, 2, 294, 295, 296]) == (
        "Frame 7 of TRAJEC.xyz: 6 selected atoms (source atoms 1-3,295-297)"
    )
    assert vs.selection_comment(0, "x.xyz", [4]) == "Frame 0 of x.xyz: 1 selected atom (source atoms 5)"


def test_export_problem_for_an_empty_selection():
    assert vs.selection_export_problem(FRAMES, []) == vs.EMPTY_SELECTION_MESSAGE


def test_export_problem_names_the_first_short_frame():
    frames = FRAMES + [FRAMES[0][:2]]
    message = vs.selection_export_problem(frames, [0, 3], first_index=10)
    assert message.startswith("Frame 12 has 2 atoms")
    assert "atom 4" in message
    assert vs.selection_export_problem(FRAMES, [0, 3]) is None


def test_write_xyz_frames_writes_every_frame_in_order(tmp_path):
    from molecularViewer import read_xyz_frames

    path = tmp_path / "sel.xyz"
    reports = []
    count = vs.write_xyz_frames(
        str(path), FRAMES, [0, 3], lambda k: f"Frame {k}",
        progress=lambda done, total: reports.append((done, total)),
    )
    assert count == 2
    frames, _ = read_xyz_frames(str(path))  # the viewer's own reader reads it back
    assert [[atom[0] for atom in frame] for frame in frames] == [["O", "O"], ["O", "O"]]
    assert frames[1][1][1] == pytest.approx(FRAMES[1][3][1])
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[1] == "Frame 0" and lines[5] == "Frame 1"
    assert reports[-1] == (2, 2)


def test_write_xyz_frames_without_selection_writes_every_atom(tmp_path):
    from molecularViewer import read_xyz_frames

    path = tmp_path / "all.xyz"
    vs.write_xyz_frames(str(path), FRAMES, None, lambda k: "c")
    frames, _ = read_xyz_frames(str(path))
    assert [len(frame) for frame in frames] == [4, 4]


def test_first_index_numbers_the_comment(tmp_path):
    path = tmp_path / "one.xyz"
    vs.write_xyz_frames(str(path), [FRAMES[1]], None, lambda k: f"Frame {k}", first_index=5)
    assert path.read_text(encoding="utf-8").splitlines()[1] == "Frame 5"


# ----------------------------------------------------------------------
# Viewer selection state and click routing
# ----------------------------------------------------------------------
import asyncio
import gc

import glfw
import numpy as np

import molecularViewer as mv
from molecularViewer import MolecularViewer, MolecularViewerUI

DIMER_FRAME_0 = [
    ("O", (0.0, 0.0, 0.0)), ("H", (0.96, 0.0, 0.0)), ("H", (-0.24, 0.93, 0.0)),
    ("O", (2.9, 0.0, 0.0)), ("H", (3.4, 0.8, 0.0)), ("H", (3.4, -0.8, 0.0)),
]
# H1 has hopped to the second oxygen (0.9 A from O3, 2.0 A from O0).
DIMER_FRAME_1 = [DIMER_FRAME_0[0], ("H", (2.0, 0.0, 0.0))] + DIMER_FRAME_0[2:]


def _dimer_viewer():
    viewer = MolecularViewer()
    viewer.frames = [list(DIMER_FRAME_0), list(DIMER_FRAME_1)]
    viewer.molecule_data = viewer.frames[0]
    viewer.calculate_bonds(0)
    return viewer


def test_ctrl_toggles_one_atom_and_notifies():
    viewer = _dimer_viewer()
    calls = []
    viewer._notify_selection_changed = lambda: calls.append(viewer.get_selected_atoms())

    assert viewer.toggle_selected_atom(4) is True
    viewer.toggle_selected_atom(1)
    viewer.toggle_selected_atom(4)

    assert viewer.get_selected_atoms() == [1]
    assert calls == [[4], [1, 4], [1]]


def test_shift_selects_and_then_removes_the_whole_molecule():
    viewer = _dimer_viewer()
    viewer.toggle_selected_molecule(4)
    assert viewer.get_selected_atoms() == [3, 4, 5]
    viewer.toggle_selected_molecule(5)
    assert viewer.get_selected_atoms() == []


def test_shift_on_a_partly_selected_molecule_completes_it():
    viewer = _dimer_viewer()
    viewer.toggle_selected_atom(4)
    viewer.toggle_selected_molecule(3)
    assert viewer.get_selected_atoms() == [3, 4, 5]


def test_molecule_selection_is_frozen_when_connectivity_changes():
    viewer = _dimer_viewer()
    viewer.toggle_selected_molecule(0)
    viewer.set_recompute_bonds_per_frame(True)
    viewer.set_frame_index(1)  # H1 now bonded to O3

    assert viewer.get_selected_atoms() == [0, 1, 2]  # identities kept
    # A new Shift+click uses the bonds drawn for the displayed frame.
    assert viewer.molecule_atoms(3) == [1, 3, 4, 5]


def test_click_routing_by_modifier():
    viewer = _dimer_viewer()  # label-only mode is the default
    viewer._route_atom_click(4, glfw.MOD_CONTROL)
    assert viewer.get_selected_atoms() == [4]
    viewer._route_atom_click(0, glfw.MOD_SHIFT)
    assert viewer.get_selected_atoms() == [0, 1, 2, 4]
    viewer._route_atom_click(5, glfw.MOD_CONTROL | glfw.MOD_SHIFT)  # Ctrl wins
    assert viewer.get_selected_atoms() == [0, 1, 2, 4, 5]
    assert viewer.get_identified_atoms() == []

    viewer._route_atom_click(3, 0)  # plain click: unchanged behaviour
    assert viewer.get_identified_atoms() == [3]
    assert viewer.get_selected_atoms() == [0, 1, 2, 4, 5]


def _ortho(left, right, bottom, top, near, far):
    return np.array(
        [[2.0 / (right - left), 0, 0, 0], [0, 2.0 / (top - bottom), 0, 0],
         [0, 0, -2.0 / (far - near), 0],
         [-(right + left) / (right - left), -(top + bottom) / (top - bottom),
          -(far + near) / (far - near), 1.0]],
        dtype=float,
    )


def _picking_viewer():
    viewer = MolecularViewer()
    viewer.frames = []
    viewer.molecule_data = [("C", (5.0, 0.0, 0.0)), ("O", (-5.0, 0.0, 0.0)), ("N", (0.0, 0.0, 0.0))]
    viewer._pick_view = (np.eye(4), _ortho(-10, 10, -10, 10, -10, 10), (0, 0, 200, 200))
    return viewer


def test_mouse_modifiers_held_at_press_decide(monkeypatch):
    viewer = _picking_viewer()
    monkeypatch.setattr(mv.glfw, "get_cursor_pos", lambda window: (150.0, 100.0))  # atom C

    viewer.on_mouse_button(None, glfw.MOUSE_BUTTON_LEFT, glfw.PRESS, glfw.MOD_CONTROL)
    viewer.on_mouse_button(None, glfw.MOUSE_BUTTON_LEFT, glfw.RELEASE, 0)
    assert viewer.get_selected_atoms() == [0]

    viewer.on_mouse_button(None, glfw.MOUSE_BUTTON_LEFT, glfw.PRESS, 0)
    viewer.on_mouse_button(None, glfw.MOUSE_BUTTON_LEFT, glfw.RELEASE, glfw.MOD_CONTROL)
    assert viewer.get_selected_atoms() == [0]
    assert viewer.get_identified_atoms() == [0]  # plain press -> label


def test_set_clear_and_expand():
    viewer = _dimer_viewer()
    assert viewer.set_selected_atoms([4, 1]) is True
    assert viewer.get_selected_atoms() == [1, 4]
    assert viewer.set_selected_atoms([1, 4]) is False  # unchanged
    assert viewer.expand_selection_to_molecules() is True
    assert viewer.get_selected_atoms() == [0, 1, 2, 3, 4, 5]
    assert viewer.clear_selection() is True
    assert viewer.clear_selection() is False


def test_parse_selection_text():
    viewer = _dimer_viewer()
    assert viewer.parse_selection_text("1-3, 5") == [0, 1, 2, 4]
    assert viewer.parse_selection_text("  ") == []
    with pytest.raises(ValueError, match="Atom 7"):
        viewer.parse_selection_text("2 7")
    with pytest.raises(ValueError):
        viewer.parse_selection_text("3-1")


def test_loading_a_designed_molecule_clears_the_selection():
    viewer = _dimer_viewer()
    viewer.set_selected_atoms([0, 1])
    viewer.load_designed_molecule(list(DIMER_FRAME_0))
    assert viewer.get_selected_atoms() == []


def test_loading_a_trajectory_clears_the_selection(tmp_path):
    path = tmp_path / "t.xyz"
    path.write_text("2\nc\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    viewer = _dimer_viewer()
    viewer.set_selected_atoms([0, 1])
    try:
        asyncio.run(viewer.load_xyz_trajectory(str(path)))
    finally:
        gc.unfreeze()  # the loader freezes the heap (see test_viewerTrajectoryLoading)
    assert viewer.get_selected_atoms() == []


# ----------------------------------------------------------------------
# Show only selected: drawing filter, picking, camera
# ----------------------------------------------------------------------
def test_visible_atoms_follow_the_show_only_switch():
    viewer = _dimer_viewer()
    viewer.set_selected_atoms([3, 4, 5])
    atoms, bonds = viewer.visible_atoms_and_bonds(6, viewer.bonds)
    assert list(atoms) == list(range(6))

    viewer.set_show_selected_only(True)
    atoms, bonds = viewer.visible_atoms_and_bonds(6, viewer.bonds)
    assert atoms == [3, 4, 5]
    assert bonds == [(3, 4), (3, 5)]


def test_show_only_with_an_empty_selection_shows_everything():
    viewer = _dimer_viewer()
    viewer.set_show_selected_only(True)
    atoms, _ = viewer.visible_atoms_and_bonds(6, viewer.bonds)
    assert list(atoms) == list(range(6))
    assert viewer.camera_focus_frame() is None


def test_outlines_only_when_everything_is_shown():
    viewer = _dimer_viewer()
    viewer.set_selected_atoms([0, 9])
    assert viewer.selection_outline_atoms(6) == {0}
    viewer.set_show_selected_only(True)
    assert viewer.selection_outline_atoms(6) == set()


def test_hidden_atoms_cannot_be_picked():
    viewer = _picking_viewer()
    viewer.set_selected_atoms([0])
    viewer.set_show_selected_only(True)
    assert viewer._pick_atom_at(50, 100) is None  # O is hidden
    assert viewer._pick_atom_at(150, 100) == 0     # C is shown


def test_vectorised_pick_matches_the_per_atom_projection():
    rng = np.random.default_rng(7)
    viewer = MolecularViewer()
    viewer.frames = []
    viewer.molecule_data = [("C", tuple(p)) for p in rng.uniform(-9, 9, size=(400, 3))]
    modelview = np.eye(4)
    projection = _ortho(-10, 10, -10, 10, -10, 10)
    viewport = (0, 0, 200, 200)
    viewer._pick_view = (modelview, projection, viewport)

    def reference(px, py):  # the per-atom loop _pick_atom_at used to run
        hits = []
        for index, (_, position) in enumerate(viewer.molecule_data):
            sx, sy, sz = viewer._project_atom(position, modelview, projection, viewport)
            if 0.0 <= sz <= 1.0 and np.hypot(sx - px, sy - (200 - py)) <= viewer.pick_pixel_threshold:
                hits.append((sz, index))
        return min(hits)[1] if hits else None

    for px, py in rng.uniform(0, 200, size=(200, 2)):
        assert viewer._pick_atom_at(px, py) == reference(px, py)


def test_camera_frames_the_selected_atoms_and_refocuses():
    viewer = _dimer_viewer()
    viewer.set_selected_atoms([3, 4, 5])
    assert viewer.camera_focus_frame() is None  # switch off: frame everything

    viewer.get_active_scene_radius()  # fixes the radius
    viewer.consume_projection_update_request()
    viewer.set_show_selected_only(True)
    focus = viewer.camera_focus_frame()
    assert [atom[0] for atom in focus] == ["O", "H", "H"]
    assert viewer.get_system_center(focus)[0] == pytest.approx((2.9 + 3.4 + 3.4) / 3)
    assert viewer._fixed_scene_radius is None
    assert viewer.consume_projection_update_request() is True

    viewer.consume_projection_update_request()
    viewer.toggle_selected_atom(5)  # a change while showing only the selection refocuses
    assert viewer.consume_projection_update_request() is True


# ----------------------------------------------------------------------
# Display-tab plumbing and export (stub widgets; no window)
# ----------------------------------------------------------------------
from types import SimpleNamespace

import toga


class _StubLoop:
    def is_closed(self):
        return False

    def call_soon_threadsafe(self, callback, *args):
        callback(*args)

    def create_task(self, coro):
        coro.close()


class _DialogWindow:
    def __init__(self, answer=None):
        self.answer = answer
        self.dialogs = []

    async def dialog(self, dialog):
        self.dialogs.append(dialog)
        return self.answer if isinstance(dialog, toga.SaveFileDialog) else None


def _message(dialog):
    return dialog._impl.message


def _selection_ui(tmp_path, answer_name="out.xyz"):
    ui = MolecularViewerUI.__new__(MolecularViewerUI)
    MolecularViewer.__init__(ui)
    ui.loop = _StubLoop()
    ui.loading_label = None
    ui.textInput_file = SimpleNamespace(value=str(tmp_path / "TRAJEC.xyz"))
    ui.selection_input = SimpleNamespace(value="")
    ui.selection_summary_label = SimpleNamespace(text="")
    ui.show_selected_switch = SimpleNamespace(value=False)
    ui.export_selected_switch = SimpleNamespace(value=False)
    ui.save_trajectory_button = SimpleNamespace(enabled=True)
    ui.main_window = _DialogWindow(str(tmp_path / answer_name))
    ui.frames = [list(DIMER_FRAME_0), list(DIMER_FRAME_1)]
    ui.molecule_data = ui.frames[0]
    ui.calculate_bonds(0)
    return ui


def test_canvas_selection_fills_the_field_and_summary(tmp_path):
    ui = _selection_ui(tmp_path)
    ui._route_atom_click(4, glfw.MOD_CONTROL)
    assert ui.selection_input.value == "5"
    assert ui.selection_summary_label.text == "1 atom: H"
    ui._route_atom_click(0, glfw.MOD_SHIFT)
    assert ui.selection_input.value == "1-3, 5"


def test_typed_ranges_apply_and_normalise(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.selection_input.value = "6,4 5"
    asyncio.run(ui.apply_selection_field(None))
    assert ui.get_selected_atoms() == [3, 4, 5]
    assert ui.selection_input.value == "4-6"


def test_a_bad_field_reports_and_keeps_the_selection(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_selected_atoms([0])
    ui.selection_input.value = "9"
    asyncio.run(ui.apply_selection_field(None))
    assert isinstance(ui.main_window.dialogs[-1], toga.ErrorDialog)
    assert "Atom 9" in _message(ui.main_window.dialogs[-1])
    assert ui.get_selected_atoms() == [0]


def test_clear_and_expand_buttons(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_selected_atoms([4])
    ui.expand_selection(None)
    assert ui.get_selected_atoms() == [3, 4, 5]
    ui.clear_atom_selection(None)
    assert ui.get_selected_atoms() == []
    assert ui.selection_input.value == ""
    assert ui.selection_summary_label.text == "No atoms selected."


def test_switches_drive_the_viewer_flags(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.show_selected_switch.value = True
    ui.toggle_show_selected_only(ui.show_selected_switch)
    ui.export_selected_switch.value = True
    ui.toggle_export_selected_only(ui.export_selected_switch)
    assert ui.show_selected_only is True and ui.export_selected_only is True


def test_unticked_current_frame_export_is_byte_identical(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_frame_index(1)
    asyncio.run(ui.save_current_frame_xyz(None))

    expected = tmp_path / "expected.xyz"
    with open(expected, "w", encoding="utf-8") as f:  # the pre-selection writer, verbatim
        f.write(f"{len(DIMER_FRAME_1)}\n")
        f.write("Frame 1 exported from molecularViewer\n")
        for element, (x, y, z) in DIMER_FRAME_1:
            f.write(f"{element:<3s} {x: .8f} {y: .8f} {z: .8f}\n")
    assert (tmp_path / "out.xyz").read_bytes() == expected.read_bytes()
    assert isinstance(ui.main_window.dialogs[-1], toga.InfoDialog)


def test_selected_current_frame_export(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_selected_atoms([3, 4, 5])
    ui.export_selected_only = True
    asyncio.run(ui.save_current_frame_xyz(None))
    lines = (tmp_path / "out.xyz").read_text(encoding="utf-8").splitlines()
    assert lines[0] == "3"
    assert lines[1] == "Frame 0 of TRAJEC.xyz: 3 selected atoms (source atoms 4-6)"
    assert [line.split()[0] for line in lines[2:]] == ["O", "H", "H"]


def test_empty_selection_export_writes_nothing(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.export_selected_only = True
    asyncio.run(ui.save_current_frame_xyz(None))
    asyncio.run(ui.save_trajectory_xyz(None))
    assert not (tmp_path / "out.xyz").exists()
    assert len(ui.main_window.dialogs) == 2  # no file dialog was opened
    assert all(isinstance(d, toga.ErrorDialog) for d in ui.main_window.dialogs)
    assert all(_message(d) == vs.EMPTY_SELECTION_MESSAGE for d in ui.main_window.dialogs)


def test_selected_trajectory_export(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_selected_atoms([0, 1, 2])
    ui.export_selected_only = True
    asyncio.run(ui.save_trajectory_xyz(None))
    frames, _ = mv.read_xyz_frames(str(tmp_path / "out.xyz"))
    assert [len(frame) for frame in frames] == [3, 3]
    assert frames[1][1] == ("H", (2.0, 0.0, 0.0))  # the hopped H stays in, by identity
    lines = (tmp_path / "out.xyz").read_text(encoding="utf-8").splitlines()
    assert lines[6] == "Frame 1 of TRAJEC.xyz: 3 selected atoms (source atoms 1-3)"
    assert ui.save_trajectory_button.enabled is True
    assert isinstance(ui.main_window.dialogs[-1], toga.InfoDialog)


def test_full_trajectory_export_when_unticked(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_selected_atoms([0])  # ignored: the checkbox is off
    asyncio.run(ui.save_trajectory_xyz(None))
    frames, _ = mv.read_xyz_frames(str(tmp_path / "out.xyz"))
    assert [len(frame) for frame in frames] == [6, 6]


def test_trajectory_export_refuses_a_changing_atom_count(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.frames.append(list(DIMER_FRAME_0[:3]))
    ui.set_selected_atoms([4])
    ui.export_selected_only = True
    asyncio.run(ui.save_trajectory_xyz(None))
    assert not (tmp_path / "out.xyz").exists()
    assert isinstance(ui.main_window.dialogs[-1], toga.ErrorDialog)
    assert _message(ui.main_window.dialogs[-1]).startswith("Frame 2 has 3 atoms")


def test_trajectory_export_refuses_to_overwrite_the_loaded_file(tmp_path):
    ui = _selection_ui(tmp_path, answer_name="TRAJEC.xyz")
    asyncio.run(ui.save_trajectory_xyz(None))
    assert not (tmp_path / "TRAJEC.xyz").exists()
    assert isinstance(ui.main_window.dialogs[-1], toga.ErrorDialog)


def test_trajectory_export_needs_a_loaded_structure(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.frames = []  # only the built-in placeholder remains in molecule_data
    asyncio.run(ui.save_trajectory_xyz(None))
    assert len(ui.main_window.dialogs) == 1
    assert isinstance(ui.main_window.dialogs[0], toga.ErrorDialog)
