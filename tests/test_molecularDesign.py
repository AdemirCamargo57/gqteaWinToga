"""Tests for the molecular design controller.

The controller holds every rule and every piece of state and touches no widget,
so these tests construct it with ``__new__`` plus ``_init_state`` and attach
stub views. Nothing here creates a Toga window.

Three things are under test: the canvas state machine (click to place, drag to
bond or to grow a bonded atom, click to cycle a bond type, right-click to
delete), the 3D view (left-drag on empty space rotates, right-drag pans, the
wheel zooms, and none is an edit), and the plumbing around it -- broadcasting to several
views at once, surviving a closed window, writing an optimized geometry back
onto the drawing, and handing the finished structure to the 3D viewer.
"""
import asyncio
from pathlib import Path

import numpy as np
import pytest

import molecularDesign
import molecularFileFormats
from molecularDesign import MolecularDesignController, blend_color
from molecularFileFormats import parse_structure_text
from molecularGeometry import measure
from molecularPreOptimizer import PIXELS_PER_ANGSTROM, OptimizationResult, reference_bond_length
from molecularSketch import RESONANCE_ORDER


class _StubView:
    """Records every hook the controller broadcasts."""

    def __init__(self):
        self.element = None
        self.redraw_count = 0
        self.status = ""

    def on_element_changed(self, symbol):
        self.element = symbol

    def on_sketch_changed(self):
        self.redraw_count += 1

    def on_status_changed(self, message):
        self.status = message


class _StubViewer:
    def __init__(self, ready=True, error=None):
        self.loaded_frame = None
        self.loaded_bonds = None
        self.ensure_calls = 0
        self.main_window = None
        self._ready = ready
        self._error = error

    def load_designed_molecule(self, frame, bonds=None):
        self.loaded_frame = frame
        self.loaded_bonds = bonds

    async def ensure_viewer_window(self, timeout=5.0):
        self.ensure_calls += 1
        return self._ready, self._error


class _FailingViewer(_StubViewer):
    def load_designed_molecule(self, frame, bonds=None):
        raise ValueError("the structure could not be loaded")


def _controller(viewer=None) -> MolecularDesignController:
    controller = MolecularDesignController.__new__(MolecularDesignController)
    controller._init_state(viewer=viewer)
    # Never read the real keyboard: tests set controller.held to simulate keys.
    controller.held = frozenset()
    controller.modifier_state = lambda: controller.held
    controller.test_view = _StubView()
    controller.register_view(controller.test_view)
    return controller


def _click(controller, x, y):
    """A press and release at the same point: no drag."""
    controller.on_canvas_press(None, x, y)
    controller.on_canvas_release(None, x, y)


def _drag(controller, from_x, from_y, to_x, to_y):
    controller.on_canvas_press(None, from_x, from_y)
    controller.on_canvas_drag(None, (from_x + to_x) / 2.0, (from_y + to_y) / 2.0)
    controller.on_canvas_drag(None, to_x, to_y)
    controller.on_canvas_release(None, to_x, to_y)


def _right_click(controller, x, y):
    """A right press and release at the same point: delete, not rotate."""
    controller.on_canvas_alt_press(None, x, y)
    controller.on_canvas_alt_release(None, x, y)


def _right_drag(controller, from_x, from_y, to_x, to_y):
    controller.on_canvas_alt_press(None, from_x, from_y)
    controller.on_canvas_alt_drag(None, (from_x + to_x) / 2.0, (from_y + to_y) / 2.0)
    controller.on_canvas_alt_drag(None, to_x, to_y)
    controller.on_canvas_alt_release(None, to_x, to_y)


def _two_bonded_carbons(controller):
    _click(controller, 100.0, 100.0)
    _click(controller, 200.0, 100.0)
    _drag(controller, 100.0, 100.0, 200.0, 100.0)
    return controller.sketch.bonds[0]


def _successful_result() -> OptimizationResult:
    return OptimizationResult(
        success=True,
        elements=["O", "H", "H"],
        coordinates=np.zeros((3, 3)),
        energy=-1.25,
        gradient_norm=0.002,
        iterations=17,
        message="CONVERGENCE: NORM OF PROJECTED GRADIENT",
        warnings=[],
    )


# ------------------------------------------------------------------
# Views and broadcasting
# ------------------------------------------------------------------
def test_carbon_is_selected_by_default():
    assert _controller().selected_element == "C"


def test_selecting_an_element_changes_what_a_click_places():
    controller = _controller()

    controller.select_element("N")
    _click(controller, 50.0, 50.0)

    assert controller.sketch.atoms[0].element == "N"


def test_a_selection_reaches_every_open_view():
    """The periodic-table window and the canvas window are separate windows:
    choosing an element in one must immediately reach the other."""
    controller = _controller()
    table_view = _StubView()
    canvas_view = _StubView()
    controller.register_view(table_view)
    controller.register_view(canvas_view)

    controller.select_element("S")

    assert table_view.element == "S"
    assert canvas_view.element == "S"


def test_an_element_chosen_in_one_view_places_atoms_in_another():
    controller = _controller()
    canvas_view = _StubView()
    controller.register_view(canvas_view)

    controller.select_element("P")
    _click(controller, 80.0, 80.0)

    assert controller.sketch.atoms[0].element == "P"


def test_edits_redraw_every_open_view():
    controller = _controller()
    second_view = _StubView()
    controller.register_view(second_view)

    _click(controller, 120.0, 80.0)

    assert controller.test_view.redraw_count > 0
    assert second_view.redraw_count > 0


def test_status_messages_reach_every_open_view():
    controller = _controller()
    second_view = _StubView()
    controller.register_view(second_view)

    _click(controller, 120.0, 80.0)

    assert "Placed" in controller.test_view.status
    assert "Placed" in second_view.status


def test_a_closed_view_stops_receiving_updates():
    controller = _controller()
    closing_view = _StubView()
    controller.register_view(closing_view)
    controller.unregister_view(closing_view)
    redraws_at_close = closing_view.redraw_count

    _click(controller, 120.0, 80.0)

    assert closing_view.redraw_count == redraws_at_close
    assert controller.test_view.redraw_count > redraws_at_close


def test_the_drawing_survives_closing_every_view():
    """Windows are disposable views; a closed window must not lose the work."""
    controller = _controller()
    _two_bonded_carbons(controller)

    controller.unregister_view(controller.test_view)

    assert controller.sketch.atom_count == 2
    assert controller.sketch.bond_count == 1


def test_a_reopened_view_is_given_the_current_state():
    controller = _controller()
    controller.select_element("O")
    _two_bonded_carbons(controller)

    reopened = _StubView()
    controller.register_view(reopened)

    assert reopened.element == "O"
    assert reopened.redraw_count > 0


# ------------------------------------------------------------------
# Placing atoms
# ------------------------------------------------------------------
def test_clicking_empty_space_places_an_atom():
    controller = _controller()

    _click(controller, 120.0, 80.0)

    assert controller.sketch.atom_count == 1
    assert (controller.sketch.atoms[0].x, controller.sketch.atoms[0].y) == (120.0, 80.0)


def test_dragging_across_empty_space_places_nothing():
    controller = _controller()

    _drag(controller, 20.0, 20.0, 200.0, 200.0)

    assert controller.sketch.atom_count == 0


def test_clicking_an_existing_atom_places_nothing_on_top_of_it():
    controller = _controller()
    _click(controller, 100.0, 100.0)

    _click(controller, 102.0, 101.0)

    assert controller.sketch.atom_count == 1


# ------------------------------------------------------------------
# Drawing bonds
# ------------------------------------------------------------------
def test_dragging_from_one_atom_to_another_creates_a_single_bond():
    controller = _controller()
    _click(controller, 100.0, 100.0)
    _click(controller, 200.0, 100.0)

    _drag(controller, 100.0, 100.0, 200.0, 100.0)

    assert controller.sketch.bond_count == 1
    assert controller.sketch.bonds[0].order == 1


def test_dragging_from_an_atom_back_to_itself_creates_no_bond():
    controller = _controller()
    _click(controller, 100.0, 100.0)

    _drag(controller, 100.0, 100.0, 103.0, 103.0)

    assert controller.sketch.bond_count == 0


def test_dragging_from_an_atom_to_empty_space_grows_a_bonded_atom():
    controller = _controller()
    _click(controller, 100.0, 100.0)
    controller.select_element("O")

    _drag(controller, 100.0, 100.0, 300.0, 100.0)

    assert controller.sketch.atom_count == 2
    assert controller.sketch.bond_count == 1
    assert controller.sketch.atoms[1].element == "O"


def test_a_grown_atom_sits_at_the_reference_bond_length_in_the_drag_direction():
    controller = _controller()
    _click(controller, 100.0, 100.0)

    _drag(controller, 100.0, 100.0, 100.0, 400.0)

    grown = controller.sketch.atoms[1]
    expected = reference_bond_length("C", "C") * PIXELS_PER_ANGSTROM
    assert grown.x == pytest.approx(100.0)
    assert grown.y == pytest.approx(100.0 + expected)


def test_a_grown_atom_keeps_the_reference_length_when_zoomed():
    """Zooming changes how far apart atoms look, not how far apart they are."""
    controller = _controller()
    _click(controller, 100.0, 100.0)
    controller.zoom_in(None)
    sx, sy, _ = controller.screen_positions()[controller.sketch.atoms[0].atom_id]

    _drag(controller, sx, sy, sx + 300.0, sy)

    anchor, grown = controller.sketch.positions()
    expected = reference_bond_length("C", "C") * PIXELS_PER_ANGSTROM
    assert np.linalg.norm(grown - anchor) == pytest.approx(expected)


def test_a_grown_atom_is_placed_at_the_depth_of_its_anchor():
    controller = _controller()
    controller.sketch.add_atom("C", 100.0, 100.0, 37.0)

    _drag(controller, 100.0, 100.0, 300.0, 100.0)

    assert controller.sketch.atoms[1].z == pytest.approx(37.0)


def test_an_element_without_a_radius_grows_at_a_default_length():
    controller = _controller()
    _click(controller, 100.0, 100.0)
    controller.select_element("Xx")

    _drag(controller, 100.0, 100.0, 300.0, 100.0)

    expected = molecularDesign.DEFAULT_GROWN_BOND_LENGTH * PIXELS_PER_ANGSTROM
    assert controller.sketch.atoms[1].x == pytest.approx(100.0 + expected)


def test_growing_an_atom_is_undone_in_one_step():
    controller = _controller()
    _click(controller, 100.0, 100.0)
    _drag(controller, 100.0, 100.0, 300.0, 100.0)

    controller.undo(None)

    assert controller.sketch.atom_count == 1
    assert controller.sketch.bond_count == 0


def test_dragging_between_two_already_bonded_atoms_leaves_the_type_alone():
    controller = _controller()
    bond = _two_bonded_carbons(controller)
    bond.order = 2

    _drag(controller, 100.0, 100.0, 200.0, 100.0)

    assert controller.sketch.bond_count == 1
    assert bond.order == 2


def test_a_rubber_band_is_tracked_while_dragging_from_an_atom():
    controller = _controller()
    _click(controller, 100.0, 100.0)

    controller.on_canvas_press(None, 100.0, 100.0)
    controller.on_canvas_drag(None, 160.0, 140.0)

    assert controller.drag_preview == (100.0, 100.0, 160.0, 140.0)


def test_the_rubber_band_is_cleared_on_release():
    controller = _controller()
    _click(controller, 100.0, 100.0)

    _drag(controller, 100.0, 100.0, 300.0, 300.0)

    assert controller.drag_preview is None


# ------------------------------------------------------------------
# Bond types, including resonance
# ------------------------------------------------------------------
def test_clicking_a_bond_raises_its_order():
    controller = _controller()
    bond = _two_bonded_carbons(controller)

    _click(controller, 150.0, 100.0)

    assert bond.order == 2


def test_a_fourth_click_makes_a_carbon_carbon_bond_a_resonance_bond():
    controller = _controller()
    bond = _two_bonded_carbons(controller)

    for _ in range(3):
        _click(controller, 150.0, 100.0)

    assert bond.order == RESONANCE_ORDER


def test_clicking_past_resonance_returns_the_bond_to_single():
    controller = _controller()
    bond = _two_bonded_carbons(controller)

    for _ in range(4):
        _click(controller, 150.0, 100.0)

    assert bond.order == 1


def test_creating_a_resonance_bond_is_reported_by_name():
    controller = _controller()
    _two_bonded_carbons(controller)

    for _ in range(3):
        _click(controller, 150.0, 100.0)

    assert "resonance" in controller.status_message.lower()


def test_a_capped_bond_says_why_it_cannot_rise():
    controller = _controller()
    controller.select_element("O")
    _click(controller, 100.0, 100.0)
    controller.select_element("H")
    _click(controller, 200.0, 100.0)
    _drag(controller, 100.0, 100.0, 200.0, 100.0)

    _click(controller, 150.0, 100.0)

    assert controller.sketch.bonds[0].order == 1
    assert "single" in controller.status_message.lower()


def test_clicking_a_bond_never_places_an_atom():
    controller = _controller()
    _two_bonded_carbons(controller)

    _click(controller, 150.0, 100.0)

    assert controller.sketch.atom_count == 2


def test_a_resonance_bond_can_be_deleted_like_any_other():
    controller = _controller()
    _two_bonded_carbons(controller)
    for _ in range(3):
        _click(controller, 150.0, 100.0)

    _right_click(controller, 150.0, 100.0)

    assert controller.sketch.bond_count == 0


# ------------------------------------------------------------------
# Deleting
# ------------------------------------------------------------------
def test_right_clicking_an_atom_deletes_it_and_its_bonds():
    controller = _controller()
    _two_bonded_carbons(controller)

    _right_click(controller, 100.0, 100.0)

    assert controller.sketch.atom_count == 1
    assert controller.sketch.bond_count == 0


def test_right_clicking_a_bond_deletes_only_the_bond():
    controller = _controller()
    _two_bonded_carbons(controller)

    _right_click(controller, 150.0, 100.0)

    assert controller.sketch.atom_count == 2
    assert controller.sketch.bond_count == 0


def test_right_clicking_empty_space_does_nothing():
    controller = _controller()
    _two_bonded_carbons(controller)

    _right_click(controller, 500.0, 500.0)

    assert controller.sketch.atom_count == 2
    assert controller.sketch.bond_count == 1


# ------------------------------------------------------------------
# Undo and clear
# ------------------------------------------------------------------
def test_undo_removes_the_last_placed_atom():
    controller = _controller()
    _click(controller, 100.0, 100.0)
    _click(controller, 200.0, 100.0)

    controller.undo(None)

    assert controller.sketch.atom_count == 1


def test_undo_restores_a_deleted_atom():
    controller = _controller()
    _two_bonded_carbons(controller)

    _right_click(controller, 100.0, 100.0)
    controller.undo(None)

    assert controller.sketch.atom_count == 2
    assert controller.sketch.bond_count == 1


def test_undo_restores_a_bond_type():
    controller = _controller()
    _two_bonded_carbons(controller)
    for _ in range(3):
        _click(controller, 150.0, 100.0)

    controller.undo(None)

    assert controller.sketch.bonds[0].order == 3


def test_undo_on_an_empty_history_is_harmless():
    controller = _controller()

    controller.undo(None)

    assert controller.sketch.atom_count == 0


def test_clearing_the_canvas_empties_the_sketch():
    controller = _controller()
    _two_bonded_carbons(controller)

    controller.clear_canvas(None)

    assert controller.sketch.atom_count == 0


def test_clearing_the_canvas_can_be_undone():
    controller = _controller()
    _two_bonded_carbons(controller)

    controller.clear_canvas(None)
    controller.undo(None)

    assert controller.sketch.atom_count == 2


# ------------------------------------------------------------------
# Pre-optimization reporting
# ------------------------------------------------------------------
def test_a_successful_result_reports_energy_and_atom_count():
    summary = _controller().describe_result(_successful_result())

    assert "3 atoms" in summary
    assert "-1.25" in summary
    assert "17" in summary


def test_a_failed_result_reports_the_minimizer_message():
    result = _successful_result()
    result.success = False
    result.message = "STOP: TOTAL NO. of ITERATIONS REACHED LIMIT"

    summary = _controller().describe_result(result)

    assert "ITERATIONS REACHED LIMIT" in summary


def test_validation_warnings_appear_in_the_summary():
    result = _successful_result()
    result.warnings = ["Isolated resonance bond between atom 1 (C) and atom 2 (C)."]

    summary = _controller().describe_result(result)

    assert "Isolated resonance bond" in summary


def test_applying_a_result_stores_it_and_reports_success():
    controller = _controller()

    controller.apply_optimization_result(_successful_result())

    assert controller.last_result is not None
    assert "converged" in controller.status_message.lower()


def test_editing_the_sketch_invalidates_a_previous_result():
    """A structure that no longer matches the drawing must not be sent on."""
    controller = _controller()
    controller.apply_optimization_result(_successful_result())

    _click(controller, 300.0, 300.0)

    assert controller.last_result is None


# ------------------------------------------------------------------
# Transfer to the 3D viewer
# ------------------------------------------------------------------
def test_transfer_hands_over_the_frame_and_the_drawn_connectivity():
    viewer = _StubViewer()
    controller = _controller(viewer=viewer)
    _two_bonded_carbons(controller)
    controller.apply_optimization_result(_successful_result())

    transferred, _ = controller.transfer_to_viewer()

    assert transferred is True
    assert viewer.loaded_frame[0][0] == "O"
    assert viewer.loaded_bonds == [(0, 1)]


def test_transfer_before_optimizing_reports_why_it_cannot():
    viewer = _StubViewer()
    controller = _controller(viewer=viewer)

    transferred, message = controller.transfer_to_viewer()

    assert transferred is False
    assert viewer.loaded_frame is None
    assert "optimize" in message.lower()


def test_transfer_reports_a_viewer_that_refuses_the_structure():
    controller = _controller(viewer=_FailingViewer())
    controller.apply_optimization_result(_successful_result())

    transferred, message = controller.transfer_to_viewer()

    assert transferred is False
    assert "could not be loaded" in message


def test_sending_opens_the_3d_viewer_window():
    viewer = _StubViewer()
    controller = _controller(viewer=viewer)
    controller.apply_optimization_result(_successful_result())

    asyncio.run(controller.send_to_viewer(None))

    assert viewer.ensure_calls == 1
    assert "3D viewer" in controller.status_message


def test_a_viewer_that_cannot_open_is_reported_clearly():
    viewer = _StubViewer(ready=False, error="GLFW window creation failed")
    controller = _controller(viewer=viewer)
    controller.apply_optimization_result(_successful_result())

    asyncio.run(controller.send_to_viewer(None))

    assert "GLFW window creation failed" in controller.status_message


def test_a_failed_transfer_never_tries_to_open_the_viewer():
    viewer = _StubViewer()
    controller = _controller(viewer=viewer)

    asyncio.run(controller.send_to_viewer(None))

    assert viewer.ensure_calls == 0
    assert "optimize" in controller.status_message.lower()


def test_sending_without_a_viewer_attached_is_reported():
    controller = _controller(viewer=None)
    controller.apply_optimization_result(_successful_result())

    asyncio.run(controller.send_to_viewer(None))

    assert "viewer" in controller.status_message.lower()


# ------------------------------------------------------------------
# The 3D view: rotate (left-drag), pan (right-drag), zoom (wheel) -- the 3D
# viewer's mouse controls -- none of them an edit
# ------------------------------------------------------------------
def test_dragging_empty_space_rotates_the_view():
    controller = _controller()
    _two_bonded_carbons(controller)

    _drag(controller, 300.0, 300.0, 380.0, 300.0)

    assert not np.allclose(controller.camera.rotation, np.eye(3))
    assert controller.sketch.atom_count == 2          # no atom placed at the release point
    assert "rotated" in controller.status_message


def test_dragging_a_bond_rotates_the_view_and_keeps_its_type():
    controller = _controller()
    bond = _two_bonded_carbons(controller)

    _drag(controller, 150.0, 100.0, 230.0, 160.0)

    assert not np.allclose(controller.camera.rotation, np.eye(3))
    assert bond.order == 1


def test_dragging_from_an_atom_still_draws_a_bond():
    controller = _controller()
    _click(controller, 100.0, 100.0)
    _click(controller, 200.0, 100.0)

    _drag(controller, 100.0, 100.0, 200.0, 100.0)

    assert controller.sketch.bond_count == 1
    assert np.allclose(controller.camera.rotation, np.eye(3))


def test_right_dragging_pans_the_view():
    controller = _controller()
    _click(controller, 100.0, 100.0)

    _right_drag(controller, 300.0, 300.0, 340.0, 280.0)

    sx, sy, _ = controller.screen_positions()[controller.sketch.atoms[0].atom_id]
    assert (sx, sy) == pytest.approx((140.0, 80.0))
    assert np.allclose(controller.camera.rotation, np.eye(3))


def test_panning_moves_the_view_not_the_atoms():
    controller = _controller()
    _two_bonded_carbons(controller)
    before = controller.sketch.positions()
    undo_depth = len(controller.undo_stack)

    _right_drag(controller, 300.0, 300.0, 340.0, 280.0)

    assert np.array_equal(controller.sketch.positions(), before)
    assert len(controller.undo_stack) == undo_depth


def test_right_dragging_over_an_atom_pans_instead_of_deleting():
    controller = _controller()
    _two_bonded_carbons(controller)

    _right_drag(controller, 100.0, 100.0, 180.0, 160.0)

    assert controller.sketch.atom_count == 2
    assert controller.sketch.bond_count == 1
    sx, sy, _ = controller.screen_positions()[controller.sketch.atoms[0].atom_id]
    assert (sx, sy) == pytest.approx((180.0, 160.0))


def test_rotating_is_not_an_edit():
    """A finished optimization still matches the drawing after a rotation."""
    controller = _controller()
    _two_bonded_carbons(controller)
    controller.apply_optimization_result(_successful_result())
    before = controller.sketch.positions()
    undo_depth = len(controller.undo_stack)

    _drag(controller, 300.0, 300.0, 380.0, 340.0)

    assert np.array_equal(controller.sketch.positions(), before)
    assert len(controller.undo_stack) == undo_depth
    assert controller.last_result is not None


def test_rotation_turns_about_the_molecule_centre():
    controller = _controller()
    _two_bonded_carbons(controller)  # centred on (150, 100)

    _drag(controller, 300.0, 300.0, 390.0, 250.0)

    centre = controller.sketch.positions().mean(axis=0)
    projected = controller.camera.project(centre)[0]
    assert projected[:2] == pytest.approx((150.0, 100.0))


def test_a_bond_can_still_be_clicked_after_rotating():
    controller = _controller()
    bond = _two_bonded_carbons(controller)
    _drag(controller, 300.0, 300.0, 340.0, 330.0)
    a, b = (controller.screen_positions()[atom_id] for atom_id in (bond.atom_i, bond.atom_j))

    _click(controller, (a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)

    assert bond.order == 2


def test_a_click_after_rotating_places_the_atom_on_the_plane_through_the_centre():
    controller = _controller()
    _two_bonded_carbons(controller)
    _drag(controller, 300.0, 300.0, 380.0, 330.0)
    centre = controller.sketch.positions().mean(axis=0)

    _click(controller, 400.0, 350.0)

    placed = controller.sketch.positions()[-1]
    projected = controller.camera.project(placed)[0]
    assert projected[:2] == pytest.approx((400.0, 350.0))
    assert projected[2] == pytest.approx(controller.camera.depth_of(centre))


def test_a_rotated_view_lets_you_build_out_of_the_drawing_plane():
    """The point of the 3D canvas: atoms placed after a rotation are not all
    at the same depth as the original drawing."""
    controller = _controller()
    _click(controller, 280.0, 210.0)
    _drag(controller, 300.0, 300.0, 457.0, 300.0)  # about 90 degrees

    _drag(controller, 280.0, 210.0, 450.0, 210.0)

    assert abs(controller.sketch.atoms[1].z) > 10.0


def test_the_wheel_zooms_about_the_pointer():
    controller = _controller()
    _click(controller, 100.0, 100.0)
    _click(controller, 300.0, 200.0)
    held = controller.sketch.atoms[1].atom_id

    controller.on_canvas_wheel(1, 300.0, 200.0)       # one notch away from the user: in

    assert controller.camera.zoom == pytest.approx(molecularDesign.WHEEL_ZOOM_STEP)
    assert controller.screen_positions()[held][:2] == pytest.approx((300.0, 200.0))  # stays under the pointer
    assert "Zoom 110%" in controller.status_message
    controller.on_canvas_wheel(-2, 300.0, 200.0)      # two notches back: out
    assert controller.camera.zoom == pytest.approx(1.0 / molecularDesign.WHEEL_ZOOM_STEP)


def test_the_wheel_is_not_an_edit_and_stops_at_the_zoom_limits():
    controller = _controller()
    _two_bonded_carbons(controller)
    controller.apply_optimization_result(_successful_result())
    before, undo_depth = controller.sketch.positions(), len(controller.undo_stack)

    controller.on_canvas_wheel(100, 200.0, 200.0)
    assert controller.camera.zoom == controller.camera.MAX_ZOOM
    controller.on_canvas_wheel(-200, 200.0, 200.0)
    assert controller.camera.zoom == controller.camera.MIN_ZOOM

    assert np.array_equal(controller.sketch.positions(), before)
    assert len(controller.undo_stack) == undo_depth
    assert controller.last_result is not None


def test_a_smooth_scrolling_wheel_zooms_by_a_fraction_of_a_notch():
    controller = _controller()
    redraws = controller.test_view.redraw_count
    controller.on_canvas_wheel(0, 200.0, 200.0)       # nothing turned: nothing happens
    assert controller.camera.zoom == 1.0 and controller.test_view.redraw_count == redraws

    controller.on_canvas_wheel(0.5, 200.0, 200.0)
    assert controller.camera.zoom == pytest.approx(molecularDesign.WHEEL_ZOOM_STEP ** 0.5)


def test_zoom_buttons_change_the_zoom_and_reset_restores_it():
    controller = _controller()

    controller.zoom_in(None)
    zoomed = controller.camera.zoom
    controller.zoom_out(None)
    controller.zoom_out(None)
    controller.reset_view(None)

    assert zoomed > 1.0
    assert controller.camera.zoom == pytest.approx(1.0)
    assert np.allclose(controller.camera.rotation, np.eye(3))


def test_zooming_keeps_the_canvas_centre_fixed():
    controller = _controller()
    centre_x, centre_y = controller.view_centre
    controller.sketch.add_atom("C", centre_x, centre_y)

    controller.zoom_in(None)

    sx, sy, _ = controller.screen_positions()[controller.sketch.atoms[0].atom_id]
    assert (sx, sy) == pytest.approx((centre_x, centre_y))


# ------------------------------------------------------------------
# Writing the optimized geometry back onto the canvas
# ------------------------------------------------------------------
def _write_back_result(positions) -> OptimizationResult:
    result = _successful_result()
    result.elements = ["C", "C"]
    result.coordinates = np.zeros((2, 3))
    result.sketch_positions = np.asarray(positions, dtype=float)
    return result


def test_an_optimized_geometry_is_written_back_onto_the_drawing():
    controller = _controller()
    _two_bonded_carbons(controller)
    optimized = [[110.0, 100.0, 5.0], [170.0, 100.0, -5.0]]

    controller.apply_optimization_result(_write_back_result(optimized))

    assert controller.sketch.positions().tolist() == optimized
    assert controller.last_result is not None


def test_the_write_back_can_be_undone_to_recover_the_drawing():
    controller = _controller()
    _two_bonded_carbons(controller)
    drawn = controller.sketch.positions()

    controller.apply_optimization_result(_write_back_result([[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]]))
    controller.undo(None)

    assert np.array_equal(controller.sketch.positions(), drawn)
    assert controller.last_result is None


def test_the_write_back_redraws_the_canvas():
    controller = _controller()
    _two_bonded_carbons(controller)
    redraws = controller.test_view.redraw_count

    controller.apply_optimization_result(_write_back_result([[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]]))

    assert controller.test_view.redraw_count > redraws


def test_a_result_for_a_different_atom_count_is_not_written_back():
    controller = _controller()
    _click(controller, 100.0, 100.0)
    drawn = controller.sketch.positions()

    controller.apply_optimization_result(_write_back_result([[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]]))

    assert np.array_equal(controller.sketch.positions(), drawn)


def test_pre_optimizing_updates_the_canvas_geometry():
    controller = _controller()
    _two_bonded_carbons(controller)  # drawn 100 px = 2.5 A apart

    asyncio.run(controller.pre_optimize(None))

    first, second = controller.sketch.positions()
    expected = reference_bond_length("C", "C") * PIXELS_PER_ANGSTROM
    assert np.linalg.norm(second - first) == pytest.approx(expected, rel=0.02)
    assert controller.last_result is not None


def test_an_edit_made_while_optimizing_discards_the_stale_result(monkeypatch):
    controller = _controller()
    _two_bonded_carbons(controller)
    real_optimize = molecularDesign.optimize_sketch

    def optimize_then_edit(sketch):
        result = real_optimize(sketch)
        _click(controller, 400.0, 400.0)  # the user keeps drawing meanwhile
        return result

    monkeypatch.setattr(molecularDesign, "optimize_sketch", optimize_then_edit)

    asyncio.run(controller.pre_optimize(None))

    assert controller.sketch.atom_count == 3
    assert controller.last_result is None
    assert "changed" in controller.status_message.lower()


# ------------------------------------------------------------------
# Saving the displayed structure
# ------------------------------------------------------------------
class _StubDialogWindow:
    """Answers a SaveFileDialog with ``save_path``; records every dialog."""

    def __init__(self, save_path=None):
        self.save_path = save_path
        self.dialogs = []

    async def dialog(self, dialog):
        self.dialogs.append(dialog)
        if dialog.__class__.__name__ == "SaveFileDialog":
            return self.save_path
        return None


def _parse_xyz(text):
    lines = text.splitlines()
    count = int(lines[0])
    atoms = [line.split() for line in lines[2:2 + count]]
    elements = [atom[0] for atom in atoms]
    coordinates = np.array([[float(value) for value in atom[1:4]] for atom in atoms])
    return count, lines[1], elements, coordinates


def test_the_frame_is_exported_as_elements_and_angstrom_coordinates():
    controller = _controller()
    controller.select_element("O")
    _click(controller, 100.0, 100.0)
    controller.select_element("C")
    _click(controller, 200.0, 100.0)  # 100 px = 2.5 A to the right

    count, comment, elements, coordinates = _parse_xyz(controller.frame_xyz_text())

    assert count == 2
    assert "Angstrom" in comment
    assert elements == ["O", "C"]
    assert np.allclose(coordinates, [[-1.25, 0.0, 0.0], [1.25, 0.0, 0.0]])


def test_the_exported_frame_keeps_the_drawing_handedness():
    controller = _controller()
    # The second atom is one Angstrom down the screen and one away from the viewer.
    _click(controller, 100.0, 100.0)
    _click(controller, 100.0, 140.0)
    controller.sketch.update_positions([[100.0, 100.0, 0.0], [100.0, 140.0, 40.0]])

    _, _, _, coordinates = _parse_xyz(controller.frame_xyz_text())

    # Same axes as the pre-optimizer: y up and z towards the viewer.
    assert coordinates[1, 1] - coordinates[0, 1] == pytest.approx(-1.0)
    assert coordinates[1, 2] - coordinates[0, 2] == pytest.approx(-1.0)


def test_the_exported_frame_is_the_optimized_geometry_after_pre_optimizing():
    controller = _controller()
    _two_bonded_carbons(controller)  # drawn 2.5 A apart

    asyncio.run(controller.pre_optimize(None))
    _, comment, _, coordinates = _parse_xyz(controller.frame_xyz_text())

    assert "pre-optimized" in comment
    assert np.linalg.norm(coordinates[1] - coordinates[0]) == pytest.approx(
        reference_bond_length("C", "C"), rel=0.02
    )
    assert np.allclose(coordinates, controller.last_result.coordinates, atol=1e-6)


def test_the_exported_frame_ignores_the_view_rotation():
    controller = _controller()
    _two_bonded_carbons(controller)
    before = controller.frame_xyz_text()

    _drag(controller, 100.0, 300.0, 300.0, 200.0)   # empty space: rotates the view
    _right_drag(controller, 100.0, 300.0, 300.0, 200.0)   # pans it

    assert not np.allclose(controller.camera.rotation, np.eye(3))
    assert controller.frame_xyz_text() == before


def test_writing_the_frame_creates_an_xyz_file(tmp_path):
    controller = _controller()
    _two_bonded_carbons(controller)
    output = tmp_path / "frame.xyz"

    saved, message = controller.write_frame(str(output))

    assert saved
    assert str(output) in message
    assert output.read_text(encoding="utf-8") == controller.frame_xyz_text()


def test_writing_an_empty_canvas_is_refused(tmp_path):
    controller = _controller()
    output = tmp_path / "frame.xyz"

    saved, message = controller.write_frame(str(output))

    assert not saved
    assert "Draw a structure" in message
    assert not output.exists()


def test_a_write_failure_is_reported_not_raised(tmp_path):
    controller = _controller()
    _two_bonded_carbons(controller)

    saved, message = controller.write_frame(str(tmp_path / "missing" / "frame.xyz"))

    assert not saved
    assert "could not be saved" in message


def test_saving_is_not_an_edit(tmp_path):
    controller = _controller()
    _two_bonded_carbons(controller)
    controller.apply_optimization_result(_write_back_result([[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]]))
    undo_depth = len(controller.undo_stack)

    controller.write_frame(str(tmp_path / "frame.xyz"))

    assert controller.last_result is not None
    assert len(controller.undo_stack) == undo_depth


def test_the_save_frame_button_writes_the_chosen_file(tmp_path):
    viewer = _StubViewer()
    output = tmp_path / "chosen.xyz"
    viewer.main_window = _StubDialogWindow(save_path=output)
    controller = _controller(viewer=viewer)
    _two_bonded_carbons(controller)

    assert asyncio.run(controller.save_frame(None)) is True

    assert output.read_text(encoding="utf-8") == controller.frame_xyz_text()
    assert "saved" in controller.status_message


def test_cancelling_the_save_dialog_writes_nothing(tmp_path):
    viewer = _StubViewer()
    viewer.main_window = _StubDialogWindow(save_path=None)
    controller = _controller(viewer=viewer)
    _two_bonded_carbons(controller)

    assert asyncio.run(controller.save_frame(None)) is False
    assert list(tmp_path.iterdir()) == []


def test_saving_an_empty_canvas_never_opens_the_file_dialog():
    viewer = _StubViewer()
    viewer.main_window = _StubDialogWindow(save_path="unused.xyz")
    controller = _controller(viewer=viewer)

    assert asyncio.run(controller.save_frame(None)) is False
    assert [type(d).__name__ for d in viewer.main_window.dialogs] == ["ErrorDialog"]


# ------------------------------------------------------------------
# Drawing helpers
# ------------------------------------------------------------------
def test_drawing_order_paints_back_to_front():
    controller = _controller()
    controller.sketch.add_atom("C", 100.0, 100.0, -20.0)
    controller.sketch.add_atom("O", 200.0, 100.0, 20.0)

    order = controller.drawing_order()

    assert [item.element for kind, item, _ in order if kind == "atom"] == ["O", "C"]


def test_a_bond_is_painted_before_the_atom_in_front_of_it():
    controller = _controller()
    front = controller.sketch.add_atom("C", 100.0, 100.0, -20.0)
    back = controller.sketch.add_atom("O", 200.0, 100.0, 20.0)
    controller.sketch.add_bond(front.atom_id, back.atom_id)

    kinds = [kind for kind, _, _ in controller.drawing_order()]

    assert kinds == ["atom", "bond", "atom"]


def test_a_flat_drawing_paints_bonds_under_every_atom():
    controller = _controller()
    _two_bonded_carbons(controller)

    kinds = [kind for kind, _, _ in controller.drawing_order()]

    assert kinds == ["bond", "atom", "atom"]


def test_depth_fades_are_zero_for_a_flat_drawing():
    controller = _controller()
    _two_bonded_carbons(controller)

    assert all(fade == 0.0 for _, _, fade in controller.drawing_order())


def test_the_furthest_atom_is_faded_most():
    controller = _controller()
    controller.sketch.add_atom("C", 100.0, 100.0, -20.0)
    controller.sketch.add_atom("O", 200.0, 100.0, 20.0)

    fades = {item.element: fade for kind, item, fade in controller.drawing_order() if kind == "atom"}

    assert fades["C"] == pytest.approx(0.0)
    assert fades["O"] == pytest.approx(1.0)


def test_blend_color_mixes_towards_the_background():
    assert blend_color("#000000", "#ffffff", 0.0) == "#000000"
    assert blend_color("#000000", "#ffffff", 1.0) == "#ffffff"
    assert blend_color("#ff0000", "#0000ff", 0.5) == "#800080"


# ------------------------------------------------------------------
# Importing structure files and saving in other formats
# ------------------------------------------------------------------
ETHYLENE_GJF = Path(__file__).parent / "data" / "fchk" / "inputs" / "eteno_b3lyp_631g.gjf"


class _StubFileDialogWindow(_StubDialogWindow):
    """Also answers an OpenFileDialog with ``open_path``."""

    def __init__(self, save_path=None, open_path=None):
        super().__init__(save_path)
        self.open_path = open_path

    async def dialog(self, dialog):
        if dialog.__class__.__name__ == "OpenFileDialog":
            self.dialogs.append(dialog)
            return self.open_path
        return await super().dialog(dialog)


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_importing_a_file_replaces_the_drawing_with_its_structure():
    controller = _controller()
    _click(controller, 100.0, 100.0)

    imported, message, notes = controller.import_structure(ETHYLENE_GJF)

    assert imported
    assert [atom.element for atom in controller.sketch.atoms] == ["C", "H", "H", "C", "H", "H"]
    assert sorted(bond.order for bond in controller.sketch.bonds) == [1.0, 1.0, 1.0, 1.0, 2.0]
    assert "Imported 6 atoms and 5 bonds" in message and "Undo" in message
    assert controller.status_message == message
    assert notes == []


def test_an_imported_structure_keeps_its_geometry_on_export():
    controller = _controller()
    controller.import_structure(ETHYLENE_GJF)

    exported = parse_structure_text(controller.frame_xyz_text(), "XYZ").coordinates
    original = molecularFileFormats.read_structure_file(ETHYLENE_GJF).coordinates

    assert np.allclose(exported, original - original.mean(axis=0), atol=1e-8)


def test_an_import_is_one_undo_step_and_clears_the_optimization():
    controller = _controller()
    _two_bonded_carbons(controller)
    controller.apply_optimization_result(_write_back_result([[0.0, 0.0, 0.0], [60.0, 0.0, 0.0]]))
    drawn = controller.sketch.snapshot()

    controller.import_structure(ETHYLENE_GJF)
    assert controller.last_result is None

    controller.undo(None)
    assert controller.sketch.snapshot() == drawn


def test_an_imported_molecule_is_centred_and_the_view_reset():
    controller = _controller()
    _drag(controller, 100.0, 300.0, 300.0, 200.0)         # rotate
    _right_drag(controller, 100.0, 300.0, 300.0, 200.0)   # pan
    controller.on_canvas_wheel(3, 200.0, 200.0)           # zoom

    controller.import_structure(ETHYLENE_GJF)

    assert np.allclose(controller.sketch.positions().mean(axis=0)[:2], controller.view_centre)
    assert np.allclose(controller.camera.rotation, np.eye(3))
    assert controller.camera.zoom == 1.0


def test_a_large_molecule_is_zoomed_out_to_fit(tmp_path):
    chain = "\n".join(f"C {1.5 * k:.3f} 0 0" for k in range(30))
    path = _write(tmp_path, "chain.xyz", f"30\nchain\n{chain}\n")
    controller = _controller()

    controller.import_structure(path)

    screen = controller.camera.project(controller.sketch.positions())
    assert controller.camera.zoom < 1.0
    assert screen[:, 0].min() >= 0.0 and screen[:, 0].max() <= molecularDesign.CANVAS_WIDTH


def test_a_malformed_file_is_reported_and_changes_nothing(tmp_path):
    path = _write(tmp_path, "bad.xyz", "3\ncomment\nO 0 0 0\n")
    controller = _controller()
    _click(controller, 100.0, 100.0)
    before = controller.sketch.snapshot()

    imported, message, _ = controller.import_structure(path)

    assert not imported
    assert "'bad.xyz' could not be imported as XYZ" in message
    assert controller.sketch.snapshot() == before
    assert len(controller.undo_stack) == 1


def test_an_unsupported_file_type_is_reported(tmp_path):
    imported, message, _ = _controller().import_structure(_write(tmp_path, "run.log", "text"))

    assert not imported
    assert "Gaussian output files are not read" in message


def test_too_many_atoms_are_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(molecularDesign, "MAX_IMPORT_ATOMS", 2)
    path = _write(tmp_path, "water.xyz", "3\n\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n")
    controller = _controller()

    imported, message, _ = controller.import_structure(path)

    assert not imported
    assert "at most 2" in message and controller.sketch.is_empty


def test_the_import_button_loads_the_chosen_file():
    viewer = _StubViewer()
    viewer.main_window = _StubFileDialogWindow(open_path=ETHYLENE_GJF)
    controller = _controller(viewer=viewer)

    assert asyncio.run(controller.import_file(None)) is True

    assert controller.sketch.atom_count == 6
    kinds = [type(d).__name__ for d in viewer.main_window.dialogs]
    assert kinds == ["OpenFileDialog", "InfoDialog"]


def test_a_failed_import_shows_an_error_dialog(tmp_path):
    viewer = _StubViewer()
    viewer.main_window = _StubFileDialogWindow(open_path=_write(tmp_path, "bad.mol", "x\n"))
    controller = _controller(viewer=viewer)

    assert asyncio.run(controller.import_file(None)) is False
    assert [type(d).__name__ for d in viewer.main_window.dialogs] == ["OpenFileDialog", "ErrorDialog"]


def test_cancelling_the_import_dialog_changes_nothing():
    viewer = _StubViewer()
    viewer.main_window = _StubFileDialogWindow(open_path=None)
    controller = _controller(viewer=viewer)

    assert asyncio.run(controller.import_file(None)) is False
    assert controller.sketch.is_empty


def test_xyz_stays_the_default_export_format():
    controller = _controller()
    _two_bonded_carbons(controller)

    assert controller.export_format == "XYZ"
    assert controller.frame_text() == controller.frame_xyz_text()


@pytest.mark.parametrize("key, name", [
    ("Gaussian", "frame.gjf"), ("MOL", "frame.mol"), ("SDF", "frame.sdf"),
    ("MOL2", "frame.mol2"), ("PDB", "frame.pdb"),
])
def test_the_frame_is_written_in_the_chosen_format(tmp_path, key, name):
    controller = _controller()
    controller.import_structure(ETHYLENE_GJF)
    output = tmp_path / name

    saved, message = controller.write_frame(str(output), key)

    assert saved and str(output) in message
    back = parse_structure_text(output.read_text(encoding="utf-8"), key)
    assert back.elements == ["C", "H", "H", "C", "H", "H"]
    assert len(back.bonds) == 5


def test_the_exported_bond_orders_are_the_drawn_ones(tmp_path):
    controller = _controller()
    bond = _two_bonded_carbons(controller)
    controller._cycle_bond(bond)  # double

    saved, _ = controller.write_frame(str(tmp_path / "frame.mol"), "MOL")

    assert saved
    assert parse_structure_text((tmp_path / "frame.mol").read_text(), "MOL").bonds == [(0, 1, 2.0)]


def test_choosing_a_format_is_reported_and_unknown_ones_refused():
    controller = _controller()

    controller.set_export_format("PDB")

    assert controller.export_format == "PDB"
    assert "PDB" in controller.status_message
    with pytest.raises(ValueError):
        controller.set_export_format("FCHK")


def test_save_frame_offers_the_selected_format(tmp_path):
    viewer = _StubViewer()
    output = tmp_path / "chosen.mol2"
    viewer.main_window = _StubDialogWindow(save_path=output)
    controller = _controller(viewer=viewer)
    _two_bonded_carbons(controller)
    controller.set_export_format("MOL2")

    assert asyncio.run(controller.save_frame(None)) is True

    native = viewer.main_window.dialogs[0]._impl.native  # the WinForms dialog
    assert native.Title == "Save Frame as Tripos MOL2"
    assert native.FileName == "designed_molecule.mol2"
    assert native.Filter.startswith("mol2 files (*.mol2)|*.mol2")
    assert output.read_text(encoding="utf-8") == controller.frame_text("MOL2", name="chosen")


def test_a_typed_supported_extension_wins_over_the_selection(tmp_path):
    viewer = _StubViewer()
    output = tmp_path / "typed.pdb"
    viewer.main_window = _StubDialogWindow(save_path=output)
    controller = _controller(viewer=viewer)
    _two_bonded_carbons(controller)

    assert asyncio.run(controller.save_frame(None)) is True
    assert output.read_text(encoding="utf-8").startswith("COMPND")


def test_an_unknown_extension_gets_the_selected_format(tmp_path):
    viewer = _StubViewer()
    output = tmp_path / "frame.txt"
    viewer.main_window = _StubDialogWindow(save_path=output)
    controller = _controller(viewer=viewer)
    _two_bonded_carbons(controller)

    assert asyncio.run(controller.save_frame(None)) is True
    assert output.read_text(encoding="utf-8") == controller.frame_xyz_text()


def test_a_format_that_cannot_hold_the_structure_is_reported(tmp_path, monkeypatch):
    monkeypatch.setattr(molecularFileFormats, "MOL_MAX_COUNT", 1)
    controller = _controller()
    _two_bonded_carbons(controller)

    saved, message = controller.write_frame(str(tmp_path / "frame.mol"), "MOL")

    assert not saved
    assert "could not be saved as MOL" in message
    assert not (tmp_path / "frame.mol").exists()


def test_the_save_format_list_starts_with_xyz_and_covers_every_writer():
    labels = molecularDesign.DesignCanvasWindow.EXPORT_FORMAT_LABELS
    assert list(labels.values())[0] == "XYZ"
    assert set(labels.values()) == {fmt.key for fmt in molecularFileFormats.WRITABLE_FORMATS}


# ------------------------------------------------------------------
# Selecting and moving one molecule
# ------------------------------------------------------------------
def _two_molecules(controller):
    """Ethylene's heavy atoms at (100..160, 100) and an O-H pair lower down."""
    _two_bonded_carbons(controller)                 # atoms at x=100, 200, y=100
    controller.select_element("O")
    _click(controller, 100.0, 300.0)
    controller.select_element("H")
    _click(controller, 140.0, 300.0)
    _drag(controller, 100.0, 300.0, 140.0, 300.0)
    controller.select_element("C")
    return controller.sketch.ordered_atom_ids()


def _shift_click(controller, x, y):
    controller.held = frozenset({"shift"})
    _click(controller, x, y)
    controller.held = frozenset()


def _ctrl_drag(controller, *points):
    controller.held = frozenset({"ctrl"})
    _drag(controller, *points)
    controller.held = frozenset()


def _ctrl_right_drag(controller, *points):
    controller.held = frozenset({"ctrl"})
    _right_drag(controller, *points)
    controller.held = frozenset()


def _distances(positions):
    return np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=2)


def test_shift_click_on_an_atom_selects_its_whole_molecule():
    controller = _controller()
    ids = _two_molecules(controller)

    _shift_click(controller, 100.0, 300.0)

    assert controller.selected_atoms == set(ids[2:])
    assert "selected" in controller.status_message
    assert controller.sketch.atom_count == 4  # nothing placed, nothing bonded


def test_shift_click_on_a_bond_selects_its_molecule_and_on_empty_space_clears():
    controller = _controller()
    ids = _two_molecules(controller)

    _shift_click(controller, 150.0, 100.0)  # the C-C bond
    assert controller.selected_atoms == set(ids[:2])

    _shift_click(controller, 400.0, 400.0)
    assert controller.selected_atoms == set()
    assert controller.sketch.atom_count == 4


def test_ctrl_drag_moves_only_the_selected_molecule_rigidly():
    controller = _controller()
    _two_molecules(controller)
    _shift_click(controller, 100.0, 300.0)
    before = controller.sketch.positions()
    undo_depth = len(controller.undo_stack)

    _ctrl_drag(controller, 120.0, 300.0, 170.0, 260.0)

    after = controller.sketch.positions()
    assert np.array_equal(after[:2], before[:2])
    assert np.allclose(after[2:], before[2:] + [50.0, -40.0, 0.0])
    assert len(controller.undo_stack) == undo_depth + 1  # the whole drag is one undo step
    controller.undo()
    assert np.array_equal(controller.sketch.positions(), before)


def test_ctrl_drag_moves_by_screen_pixels_on_a_rotated_zoomed_view():
    controller = _controller()
    _two_molecules(controller)
    _shift_click(controller, 100.0, 300.0)
    _drag(controller, 300.0, 200.0, 360.0, 230.0)   # empty space: rotates the view
    assert not np.allclose(controller.camera.rotation, np.eye(3))
    controller.zoom_in()
    o_screen = controller.screen_positions()[controller.sketch.ordered_atom_ids()[2]]

    _ctrl_drag(controller, o_screen[0], o_screen[1], o_screen[0] + 30.0, o_screen[1] + 10.0)

    moved = controller.screen_positions()[controller.sketch.ordered_atom_ids()[2]]
    assert np.allclose(np.subtract(moved, o_screen), [30.0, 10.0, 0.0])


def test_ctrl_right_drag_turns_the_selected_molecule_about_its_centre():
    controller = _controller()
    _two_molecules(controller)
    _shift_click(controller, 200.0, 100.0)
    before = controller.sketch.positions()

    _ctrl_right_drag(controller, 300.0, 300.0, 380.0, 340.0)

    after = controller.sketch.positions()
    assert np.allclose(_distances(after[:2]), _distances(before[:2]))
    assert np.allclose(after[:2].mean(axis=0), before[:2].mean(axis=0))
    assert not np.allclose(after[:2], before[:2])
    assert np.array_equal(after[2:], before[2:])
    assert np.allclose(controller.camera.rotation, np.eye(3))  # the view did not turn


def test_a_ctrl_right_click_never_deletes():
    controller = _controller()
    _two_molecules(controller)
    controller.held = frozenset({"ctrl"})
    _right_click(controller, 100.0, 100.0)
    assert controller.sketch.atom_count == 4


def test_moving_without_a_selection_changes_nothing():
    controller = _controller()
    _two_molecules(controller)
    before = controller.sketch.positions()

    _ctrl_drag(controller, 120.0, 300.0, 170.0, 260.0)

    assert np.array_equal(controller.sketch.positions(), before)
    assert "Shift+click a molecule" in controller.status_message


def test_a_frozen_parameter_between_molecules_blocks_the_move():
    controller = _controller()
    ids = _two_molecules(controller)
    controller.sketch.set_constraint([ids[0], ids[2]], 4.0)
    _shift_click(controller, 100.0, 300.0)
    before = controller.sketch.positions()

    _ctrl_drag(controller, 120.0, 300.0, 170.0, 260.0)

    assert np.array_equal(controller.sketch.positions(), before)
    assert "links the selection to another molecule" in controller.status_message


def test_a_selection_bonded_to_outside_atoms_blocks_the_move():
    controller = _controller()
    ids = _two_molecules(controller)
    _shift_click(controller, 100.0, 300.0)
    controller.sketch.add_bond(ids[0], ids[2])

    _ctrl_drag(controller, 120.0, 300.0, 170.0, 260.0)

    assert "bonded to atoms outside it" in controller.status_message


# ------------------------------------------------------------------
# Measuring
# ------------------------------------------------------------------
def _chain(controller):
    """C1-C2-C3-C4 zigzag, one sketch pixel = 1/40 A."""
    controller.import_structure(_chain_file())
    return controller.sketch.ordered_atom_ids()


def _chain_file(directory=None):
    import tempfile
    directory = Path(directory or tempfile.mkdtemp())
    path = directory / "chain.xyz"
    path.write_text("4\nchain\nC 0 0 0\nC 1.52 0 0\nC 2.03 1.43 0\nC 3.55 1.43 0.6\n", encoding="utf-8")
    return path


def _click_atom(controller, atom_id):
    x, y, _ = controller.screen_positions()[atom_id]
    _click(controller, x, y)


def test_clicking_atoms_in_order_measures_length_angle_and_dihedral():
    controller = _controller()
    ids = _chain(controller)

    _click_atom(controller, ids[0])
    _click_atom(controller, ids[1])
    assert controller.measurement_text() == "Bond length 1-2 (C-C): 1.5200 Å"
    _click_atom(controller, ids[2])
    assert controller.measurement_text().startswith("Bond angle 1-2-3 (C-C-C): ")
    assert controller.measurement_text().endswith("°")
    _click_atom(controller, ids[3])
    kind, value = controller.current_measurement()
    assert kind == "dihedral" and value == pytest.approx(
        measure("dihedral", molecularFileFormats.read_structure_file(_chain_file()).coordinates, (0, 1, 2, 3)))
    assert controller.sketch.atom_count == 4  # picking never draws


def test_a_picked_atom_clicked_again_is_unpicked_and_a_fifth_pick_restarts():
    controller = _controller()
    ids = _chain(controller)
    for atom_id in ids:
        _click_atom(controller, atom_id)
    _click_atom(controller, ids[1])  # already picked: unpicked
    assert controller.picked_atoms == [ids[0], ids[2], ids[3]]
    _click_atom(controller, ids[1])  # four picks again
    _click(controller, 500.0, 40.0)  # an extra atom, away from the chain
    extra = controller.sketch.ordered_atom_ids()[-1]
    _click_atom(controller, extra)  # a fifth atom starts a new measurement
    assert controller.picked_atoms == [extra]


def test_picks_of_deleted_atoms_are_dropped():
    controller = _controller()
    ids = _chain(controller)
    _click_atom(controller, ids[0])
    _click_atom(controller, ids[1])
    x, y, _ = controller.screen_positions()[ids[1]]

    _right_click(controller, x, y)

    assert controller.valid_picks() == [ids[0]]
    assert controller.current_measurement() is None


# ------------------------------------------------------------------
# Editing a parameter
# ------------------------------------------------------------------
def _pick(controller, ids):
    controller.picked_atoms = list(ids)


def test_a_bond_length_is_set_to_the_typed_value(tmp_path):
    controller = _controller()
    ids = _chain(controller)
    _pick(controller, ids[1:3])

    applied, message = controller.apply_parameter_edit("1.6000")

    assert applied and "set to 1.6000 Å" in message
    assert controller.current_measurement()[1] == pytest.approx(1.6, abs=1e-9)
    controller.undo()
    assert controller.current_measurement()[1] == pytest.approx(
        np.linalg.norm(np.subtract([2.03, 1.43, 0], [1.52, 0, 0])), abs=1e-9)


@pytest.mark.parametrize("picks, value", [((0, 1, 2), 100.0), ((0, 1, 2, 3), -65.0)])
def test_angles_and_dihedrals_are_set_and_other_bonds_keep_their_length(picks, value):
    controller = _controller()
    ids = _chain(controller)
    _pick(controller, [ids[k] for k in picks])
    lengths = [measure("distance", controller._angstrom_positions(), pair) for pair in ((0, 1), (1, 2), (2, 3))]

    applied, _ = controller.apply_parameter_edit(str(value))

    assert applied
    assert controller.current_measurement()[1] == pytest.approx(value, abs=1e-6)
    after = [measure("distance", controller._angstrom_positions(), pair) for pair in ((0, 1), (1, 2), (2, 3))]
    assert np.allclose(after, lengths)


def test_the_edit_dialog_information_names_atoms_value_and_unit():
    controller = _controller()
    ids = _chain(controller)
    _pick(controller, ids[:3])

    info = controller.parameter_edit_info()

    assert info["label"] == "Bond angle" and info["atoms"] == "1-2-3 (C-C-C)"
    assert info["unit"] == "°" and info["frozen"] is False


@pytest.mark.parametrize("text, message", [("abc", "not a number"), ("0.05", "between 0.3")])
def test_a_bad_value_is_refused_and_nothing_moves(text, message):
    controller = _controller()
    ids = _chain(controller)
    _pick(controller, ids[:2])
    before = controller.sketch.positions()

    applied, reply = controller.apply_parameter_edit(text)

    assert not applied and message in reply
    assert np.array_equal(controller.sketch.positions(), before)


def test_a_ring_bond_cannot_be_edited():
    controller = _controller()
    controller.import_structure(Path(__file__).parent / "data" / "fchk" / "inputs" / "benzene_631gd.gjf")
    carbons = [a.atom_id for a in controller.sketch.atoms if a.element == "C"][:2]
    _pick(controller, carbons)

    with pytest.raises(ValueError, match="part of a ring"):
        controller.parameter_edit_info()
    assert controller.apply_parameter_edit("1.5")[0] is False


def test_an_edit_that_would_change_a_frozen_parameter_is_refused():
    controller = _controller()
    ids = _chain(controller)
    controller.sketch.set_constraint([ids[0], ids[3]], 3.9)  # end-to-end distance
    _pick(controller, ids)

    applied, message = controller.apply_parameter_edit("180")

    assert not applied and "would change frozen parameters" in message


def test_editing_with_freeze_freezes_the_new_value():
    controller = _controller()
    ids = _chain(controller)
    _pick(controller, ids[:2])

    applied, _ = controller.apply_parameter_edit("1.45", freeze=True)

    assert applied
    assert controller.sketch.find_constraint(ids[:2]).value == pytest.approx(1.45)
    assert "frozen at 1.4500" in controller.measurement_text()


def test_editing_needs_two_to_four_picks():
    controller = _controller()
    _chain(controller)
    applied, message = controller.apply_parameter_edit("1.5")
    assert not applied and "Click 2 atoms" in message


# ------------------------------------------------------------------
# Frozen parameters and constrained pre-optimization
# ------------------------------------------------------------------
def test_freezing_keeps_the_current_value_and_the_optimized_geometry():
    controller = _controller()
    ids = _chain(controller)
    controller.last_result = _successful_result()
    _pick(controller, ids)

    frozen, message = controller.freeze_picked()

    assert frozen and "Frozen: dihedral angle 1-2-3-4" in message
    kind, atoms, value = controller.constraint_rows()[0]
    assert (kind, atoms) == ("Dihedral", "1-2-3-4") and value.endswith("°")
    assert controller.last_result is not None  # nothing moved


def test_a_linear_angle_cannot_be_frozen():
    controller = _controller()
    controller.import_structure(_write_line())
    _pick(controller, controller.sketch.ordered_atom_ids())

    frozen, message = controller.freeze_picked()

    assert not frozen and "linear angle cannot be held" in message
    assert controller.sketch.constraints == []


def _write_line():
    import tempfile
    path = Path(tempfile.mkdtemp()) / "co2.xyz"
    path.write_text("3\n\nO -1.16 0 0\nC 0 0 0\nO 1.16 0 0\n", encoding="utf-8")
    return path


def test_constraints_can_be_listed_and_removed():
    controller = _controller()
    ids = _chain(controller)
    for picks in (ids[:2], ids[:3], ids):
        _pick(controller, picks)
        controller.freeze_picked()
    assert [row[0] for row in controller.constraint_rows()] == ["Length", "Angle", "Dihedral"]

    assert controller.remove_constraint_at(1)
    assert len(controller.sketch.constraints) == 2
    assert not controller.remove_constraint_at(5)
    controller.remove_all_constraints()
    assert controller.sketch.constraints == []
    controller.undo()
    assert len(controller.sketch.constraints) == 2


def test_pre_optimizing_holds_the_frozen_parameters():
    controller = _controller()
    ids = _chain(controller)
    _pick(controller, ids)
    controller.apply_parameter_edit("75", freeze=True)
    _pick(controller, ids[:2])
    controller.apply_parameter_edit("1.65", freeze=True)

    asyncio.run(controller.pre_optimize(None))

    assert controller.last_result.constraint_count == 2
    _pick(controller, ids)
    assert controller.current_measurement()[1] == pytest.approx(75.0, abs=1e-6)
    _pick(controller, ids[:2])
    assert controller.current_measurement()[1] == pytest.approx(1.65, abs=1e-8)
    assert "Frozen parameters held: 2" in controller.describe_result(controller.last_result)


def test_an_unsupported_constraint_is_reported_by_pre_optimize():
    viewer = _StubViewer()
    viewer.main_window = _StubDialogWindow()
    controller = _controller(viewer=viewer)
    ids = _chain(controller)
    controller.sketch.set_constraint(ids[:3], 180.0)

    asyncio.run(controller.pre_optimize(None))

    assert controller.last_result is None
    dialog = viewer.main_window.dialogs[-1]
    assert type(dialog).__name__ == "ErrorDialog"


# ------------------------------------------------------------------
# Hydrogens
# ------------------------------------------------------------------
def test_hydrogens_complete_the_valences_in_one_undo_step():
    controller = _controller()
    _two_bonded_carbons(controller)

    added, message, issues = controller.add_hydrogens()

    assert added and "Added 6 hydrogens" in message and issues == []
    assert sum(atom.element == "H" for atom in controller.sketch.atoms) == 6
    assert controller.add_hydrogens()[0] is False  # nothing more to add
    controller.undo()
    assert controller.sketch.atom_count == 2


def test_hydrogens_go_to_the_selected_molecule_only():
    controller = _controller()
    ids = _two_molecules(controller)
    _shift_click(controller, 100.0, 300.0)

    added, message, _ = controller.add_hydrogens()

    assert added and "to the selected molecule" in message
    new_h = controller.sketch.ordered_atom_ids()[4:]
    assert len(new_h) == 1
    assert controller.sketch.neighbours(new_h[0]) == [ids[2]]
    assert new_h[0] in controller.selected_atoms


def test_ambiguous_atoms_are_reported_in_a_dialog():
    viewer = _StubViewer()
    viewer.main_window = _StubDialogWindow()
    controller = _controller(viewer=viewer)
    controller.select_element("N")
    _click(controller, 100.0, 100.0)
    for x, y in ((100.0, 40.0), (160.0, 100.0), (100.0, 160.0), (40.0, 100.0)):
        controller.select_element("C")
        _drag(controller, 100.0, 100.0, x, y)

    asyncio.run(controller.add_hydrogens_button(None))

    dialog = viewer.main_window.dialogs[-1]
    assert type(dialog).__name__ == "InfoDialog"
    assert controller.sketch.atom_count == 5 + 12  # the carbons completed, N left alone


# ------------------------------------------------------------------
# Background colour
# ------------------------------------------------------------------
def test_the_background_accepts_presets_and_hex_and_refuses_the_rest():
    controller = _controller()
    controller.test_view.on_background_changed = lambda color: setattr(controller.test_view, "bg", color)

    assert controller.set_canvas_background("Black") == "#000000"
    assert controller.test_view.bg == "#000000"
    assert controller.set_canvas_background("A0B1C2") == "#a0b1c2"
    with pytest.raises(ValueError, match="not a colour"):
        controller.set_canvas_background("mauve")
    assert controller.canvas_background == "#a0b1c2"
    assert controller.undo_stack == []  # a view setting, not an edit


def test_dark_backgrounds_get_light_bonds():
    assert molecularDesign.canvas_palette("#000000")["bond"] != molecularDesign.BOND_COLOR
    assert molecularDesign.canvas_palette("#ffffff")["bond"] == molecularDesign.BOND_COLOR


def test_reading_the_keyboard_never_fails():
    assert molecularDesign.read_keyboard_modifiers() <= {"shift", "ctrl"}


# ------------------------------------------------------------------
# Bringing a molecule back from the 3D viewer
# ------------------------------------------------------------------
class _ViewerWithFrame(_StubViewer):
    def __init__(self, frame, bonds):
        super().__init__()
        self.frame, self.bonds = frame, bonds

    def current_frame_with_bonds(self):
        return list(self.frame), list(self.bonds)

    def load_designed_molecule(self, frame, bonds=None):
        super().load_designed_molecule(frame, bonds)
        self.frame, self.bonds = frame, bonds


WATER_FRAME = [("O", (0.0, 0.0, 0.117)), ("H", (0.0, 0.757, -0.469)), ("H", (0.0, -0.757, -0.469))]


def test_the_viewer_molecule_comes_back_with_its_bonds():
    controller = _controller(viewer=_ViewerWithFrame(WATER_FRAME, [(0, 1), (0, 2)]))

    loaded, message, notes = controller.bring_from_viewer()

    assert loaded and "Loaded 3 atoms and 2 bonds from the 3D viewer" in message
    assert [atom.element for atom in controller.sketch.atoms] == ["O", "H", "H"]
    assert controller.sketch.bond_index_pairs() == [(0, 1), (0, 2)]
    exported = parse_structure_text(controller.frame_xyz_text(), "XYZ").coordinates
    original = np.array([xyz for _, xyz in WATER_FRAME])
    assert np.allclose(exported, original - original.mean(axis=0))
    assert "no bond orders" in notes[0]


def test_bond_orders_survive_a_round_trip_through_the_viewer():
    viewer = _ViewerWithFrame([], [])
    controller = _controller(viewer=viewer)
    bond = _two_bonded_carbons(controller)
    controller._cycle_bond(bond)                       # double
    asyncio.run(controller.pre_optimize(None))
    assert controller.transfer_to_viewer()[0]
    controller.clear_canvas()

    loaded, _, notes = controller.bring_from_viewer()

    assert loaded and controller.sketch.bonds[0].order == 2.0
    assert "restored" in notes[0]


def test_adding_beside_keeps_the_drawing_and_does_not_overlap():
    controller = _controller(viewer=_ViewerWithFrame(WATER_FRAME, [(0, 1), (0, 2)]))
    _two_bonded_carbons(controller)
    before = controller.sketch.positions()

    loaded, message, _ = controller.bring_from_viewer("add")

    after = controller.sketch.positions()
    assert loaded and message.startswith("Added 3 atoms")
    assert np.array_equal(after[:2], before)
    assert after[2:, 0].min() > before[:, 0].max()
    screen = controller.camera.project(after)
    assert screen[:, 0].min() >= 0.0 and screen[:, 0].max() <= molecularDesign.CANVAS_WIDTH


class _AnsweringWindow(_StubDialogWindow):
    def __init__(self, answer):
        super().__init__()
        self.answer = answer

    async def dialog(self, dialog):
        self.dialogs.append(dialog)
        return self.answer if type(dialog).__name__ == "QuestionDialog" else None


@pytest.mark.parametrize("answer, atoms", [(False, 2), (True, 5)])
def test_a_drawing_is_never_overwritten_by_the_viewer_molecule(answer, atoms, monkeypatch):
    viewer = _ViewerWithFrame(WATER_FRAME, [(0, 1), (0, 2)])
    viewer.main_window = _AnsweringWindow(answer)
    controller = _controller(viewer=viewer)
    monkeypatch.setattr(controller, "open_design_canvas", lambda widget=None: None)
    _two_bonded_carbons(controller)

    asyncio.run(controller.import_from_viewer(None))

    assert type(viewer.main_window.dialogs[0]).__name__ == "QuestionDialog"
    assert controller.sketch.atom_count == atoms


def test_an_empty_canvas_takes_the_viewer_molecule_without_asking(monkeypatch):
    viewer = _ViewerWithFrame(WATER_FRAME, [(0, 1), (0, 2)])
    viewer.main_window = _AnsweringWindow(False)
    controller = _controller(viewer=viewer)
    monkeypatch.setattr(controller, "open_design_canvas", lambda widget=None: None)

    assert asyncio.run(controller.import_from_viewer(None)) is True
    assert [type(d).__name__ for d in viewer.main_window.dialogs] == ["InfoDialog"]


# ------------------------------------------------------------------
# Drawing the overlays (recording context, no window)
# ------------------------------------------------------------------
class _RecordingContext:
    def __init__(self):
        self.calls = []

    def clear(self):
        self.calls.append(("clear",))

    def _record(self, name, kwargs):
        context = self

        class _Path:
            def __enter__(self_inner):
                return self_inner

            def __exit__(self_inner, *exc):
                return False

            def __getattr__(self_inner, method):
                return lambda *args, **kw: context.calls.append((name, method, kwargs.get("color"), args))
        return _Path()

    def Fill(self, **kwargs):
        return self._record("fill", kwargs)

    def Stroke(self, **kwargs):
        return self._record("stroke", kwargs)


def test_the_canvas_draws_background_selection_picks_and_measurement():
    controller = _controller()
    ids = _two_molecules(controller)
    _shift_click(controller, 100.0, 300.0)
    _pick(controller, ids[:2])
    controller.set_canvas_background("Black")
    window = molecularDesign.DesignCanvasWindow.__new__(molecularDesign.DesignCanvasWindow)
    window.controller = controller
    window.canvas = type("Canvas", (), {"context": _RecordingContext()})()

    window._draw_sketch()

    calls = window.canvas.context.calls
    assert calls[1] == ("fill", "rect", "#000000", (0, 0, 560, 420))
    assert any(call[2] == molecularDesign.PICK_COLOR and call[1] == "write_text" and call[3][0] == "#2"
               for call in calls if len(call) > 2)
    texts = [call[3][0] for call in calls if len(call) > 2 and call[1] == "write_text"]
    assert any(text.startswith("Bond length 1-2 (C-C)") for text in texts)
    halo = molecularDesign.blend_color(molecularDesign.SELECTION_COLOR, "#000000", 0.45)
    assert sum(1 for call in calls if len(call) > 2 and call[2] == halo) == 2


# ------------------------------------------------------------------
# Responsiveness: what a drag event costs
# ------------------------------------------------------------------
def test_the_frame_primitives_render_to_a_bitmap():
    """The fast path paints exactly these primitives; they must all render."""
    from molecularCanvasRender import render_scene

    controller = _controller()
    ids = _two_molecules(controller)
    _shift_click(controller, 100.0, 300.0)
    _pick(controller, ids[:3])
    window = molecularDesign.DesignCanvasWindow.__new__(molecularDesign.DesignCanvasWindow)
    window.controller = controller

    primitives = window.scene_primitives()

    kinds = {primitive[0] for primitive in primitives}
    assert {"rect", "disc", "line", "ring", "polyline", "text"} <= kinds
    assert render_scene(primitives, (700, 525), 1.25).size == (700, 525)


def test_a_ctrl_drag_follows_the_pointer_from_the_first_pixel():
    controller = _controller()
    _two_molecules(controller)
    _shift_click(controller, 100.0, 300.0)
    before = controller.sketch.positions()
    controller.held = frozenset({"ctrl"})

    controller.on_canvas_press(None, 120.0, 300.0)
    controller.on_canvas_drag(None, 122.0, 301.0)  # inside the click threshold

    assert np.allclose(controller.sketch.positions()[2:], before[2:] + [2.0, 1.0, 0.0])
    controller.on_canvas_release(None, 122.0, 301.0)


def test_a_ctrl_right_drag_turns_from_the_first_pixel():
    controller = _controller()
    _two_molecules(controller)
    _shift_click(controller, 100.0, 300.0)
    before = controller.sketch.positions()
    controller.held = frozenset({"ctrl"})

    controller.on_canvas_alt_press(None, 300.0, 300.0)
    controller.on_canvas_alt_drag(None, 302.0, 300.0)

    assert not np.allclose(controller.sketch.positions()[2:], before[2:])
    controller.on_canvas_alt_release(None, 302.0, 300.0)


def test_plain_drags_keep_the_click_threshold():
    controller = _controller()
    _click(controller, 100.0, 100.0)
    anchor = controller.camera.anchor.copy()

    controller.on_canvas_press(None, 300.0, 300.0)
    controller.on_canvas_drag(None, 302.0, 301.0)
    assert np.array_equal(controller.camera.rotation, np.eye(3))  # no turn yet: could still be a click
    controller.on_canvas_release(None, 302.0, 301.0)
    assert controller.sketch.atom_count == 2                       # ...and it was one

    controller.on_canvas_alt_press(None, 300.0, 300.0)
    controller.on_canvas_alt_drag(None, 302.0, 301.0)
    assert np.array_equal(controller.camera.anchor, anchor)        # no pan yet: could still be a delete
    controller.on_canvas_alt_release(None, 302.0, 301.0)


class _StubCanvas:
    def __init__(self):
        self.context = _RecordingContext()
        self.redraws = 0

    def redraw(self):
        self.redraws += 1


def _headless_window(controller):
    window = molecularDesign.DesignCanvasWindow.__new__(molecularDesign.DesignCanvasWindow)
    window.controller = controller
    window.canvas = _StubCanvas()
    window._bitmap_painting = True  # no native control here: falls back to Toga drawing
    window._bitmap_active = False
    window.element_label = type("Label", (), {"text": ""})()
    window.status_label = type("Label", (), {"text": ""})()
    window.panel_refreshes = 0
    window._refresh_panel = lambda: setattr(window, "panel_refreshes", window.panel_refreshes + 1)
    return window


def test_the_side_panel_waits_for_the_end_of_a_drag():
    controller = _controller()
    _two_molecules(controller)
    _shift_click(controller, 100.0, 300.0)
    window = _headless_window(controller)
    controller.register_view(window)  # registering paints once: start counting after it
    window.panel_refreshes = 0
    window.canvas.redraws = 0
    controller.held = frozenset({"ctrl"})

    controller.on_canvas_press(None, 120.0, 300.0)
    for step in range(1, 11):
        controller.on_canvas_drag(None, 120.0 + 3.0 * step, 300.0)
    assert window.canvas.redraws == 10       # every move is painted...
    assert window.panel_refreshes == 0       # ...but the panel is left alone
    controller.on_canvas_release(None, 150.0, 300.0)
    assert window.panel_refreshes >= 1
    assert not controller.gesture_active


def test_the_native_wheel_zooms_at_the_pointer_in_canvas_units():
    """WinForms reports the wheel in physical pixels; Toga's own scaling turns
    them into the canvas units every other gesture uses."""
    from types import SimpleNamespace

    controller = _controller()
    controller.sketch.add_atom("C", 300.0, 200.0)
    window = _headless_window(controller)
    window.canvas._impl = SimpleNamespace(scale_out=lambda value: value / 1.25)   # 125 % display scaling
    event = SimpleNamespace(Delta=2 * molecularDesign.WHEEL_DELTA, X=375, Y=250, Handled=False)

    window._on_canvas_wheel(None, event)

    assert controller.camera.zoom == pytest.approx(molecularDesign.WHEEL_ZOOM_STEP ** 2)
    sx, sy, _ = controller.screen_positions()[controller.sketch.atoms[0].atom_id]
    assert (sx, sy) == pytest.approx((300.0, 200.0))
    assert event.Handled is True        # the window must not zoom a second time


def test_without_a_native_control_the_wheel_is_simply_not_wired():
    controller = _controller()
    window = _headless_window(controller)
    window.window = None
    window._attach_wheel_zoom()         # no native canvas or form: nothing to attach, no error


def test_without_a_native_control_the_toga_canvas_draws():
    controller = _controller()
    _two_bonded_carbons(controller)
    window = _headless_window(controller)

    window.on_sketch_changed()

    assert window._bitmap_painting is False  # tried once, then gave up for good
    assert window.canvas.context.calls[0] == ("clear",)
    assert window.canvas.redraws == 1


# ------------------------------------------------------------------
# Replacing an atom's element
# ------------------------------------------------------------------
def test_in_replace_mode_a_click_replaces_the_atom_in_place():
    controller = _controller()
    bond = _two_bonded_carbons(controller)
    first, second = controller.sketch.ordered_atom_ids()
    controller._cycle_bond(bond)                         # double
    controller.sketch.set_constraint([first, second], 1.34)
    controller.picked_atoms = [first, second]
    controller.selected_atoms = {first, second}
    before = controller.sketch.positions()
    undo_depth = len(controller.undo_stack)

    controller.set_replace_mode(True)
    controller.select_element("N")
    _click(controller, 100.0, 100.0)

    assert controller.sketch.get_atom(first).element == "N"
    assert controller.sketch.ordered_atom_ids() == [first, second]
    assert np.array_equal(controller.sketch.positions(), before)
    assert controller.sketch.bonds[0].order == 2.0
    assert controller.picked_atoms == [first, second] and controller.selected_atoms == {first, second}
    assert controller.sketch.constraints[0].atom_ids == (first, second)
    assert len(controller.undo_stack) == undo_depth + 1
    assert "C -> N" in controller.status_message
    controller.undo()
    assert controller.sketch.get_atom(first).element == "C"


def test_replace_mode_off_a_click_still_picks():
    controller = _controller()
    _two_bonded_carbons(controller)
    controller.select_element("N")

    _click(controller, 100.0, 100.0)

    assert controller.sketch.atoms[0].element == "C"
    assert controller.picked_atoms == [controller.sketch.ordered_atom_ids()[0]]


def test_in_replace_mode_drawing_still_works():
    controller = _controller()
    _two_bonded_carbons(controller)
    controller.set_replace_mode(True)
    controller.select_element("O")

    _drag(controller, 200.0, 100.0, 260.0, 160.0)   # grow a bonded atom
    _click(controller, 400.0, 300.0)               # place one

    assert [atom.element for atom in controller.sketch.atoms] == ["C", "C", "O", "O"]
    assert controller.sketch.bond_count == 2


def test_replacing_with_the_same_element_changes_nothing():
    controller = _controller()
    _two_bonded_carbons(controller)
    controller.set_replace_mode(True)
    undo_depth = len(controller.undo_stack)

    replaced, message = controller.replace_atom(controller.sketch.ordered_atom_ids()[0])

    assert not replaced and "already C" in message
    assert len(controller.undo_stack) == undo_depth


@pytest.mark.parametrize("element, order, phrase", [
    ("H", 2.0, "double bond to atom 2 is above the single bond H usually forms"),
    ("F", 1.0, None),
    ("Og", 1.0, "no covalent radius"),
    ("He", 1.0, "He now has 1 bond where 0 is typical"),
])
def test_unusual_replacements_are_reported_and_bonds_kept(element, order, phrase):
    controller = _controller()
    bond = _two_bonded_carbons(controller)
    bond.order = order
    controller.set_replace_mode(True)
    controller.select_element(element)

    replaced, message = controller.replace_atom(controller.sketch.ordered_atom_ids()[0])

    assert replaced and controller.sketch.bonds[0].order == order
    if phrase:
        assert phrase in message
    else:
        assert "Note:" not in message


def test_an_over_valent_replacement_is_reported():
    controller = _controller()
    controller.select_element("C")
    _click(controller, 200.0, 200.0)
    for x, y in ((200.0, 140.0), (260.0, 200.0), (200.0, 260.0)):
        _drag(controller, 200.0, 200.0, x, y)
    controller.set_replace_mode(True)
    controller.select_element("O")

    _, message = controller.replace_atom(controller.sketch.ordered_atom_ids()[0])

    assert "O now has 3 bonds where 2 is typical" in message


def test_selecting_an_element_in_replace_mode_says_what_a_click_does():
    controller = _controller()
    controller.set_replace_mode(True)
    controller.select_element("Cl")
    assert "replace it with Cl" in controller.status_message


# ------------------------------------------------------------------
# Atom numbers
# ------------------------------------------------------------------
def _numbered_window(controller):
    window = _headless_window(controller)
    controller.register_view(window)
    return window


def test_atom_numbers_are_drawn_only_when_switched_on():
    controller = _controller()
    _chain(controller)
    window = _numbered_window(controller)

    assert not [p for p in window.scene_primitives() if p[0] == "label"]
    undo_depth = len(controller.undo_stack)
    controller.set_show_atom_numbers(True)
    labels = [p for p in window.scene_primitives() if p[0] == "label"]
    assert [label[3] for label in labels] == ["1", "2", "3", "4"]
    assert len(controller.undo_stack) == undo_depth  # a view setting, not an edit
    assert "Atom numbers shown" in controller.status_message


def test_numbers_follow_the_sketch_numbering_after_a_deletion():
    controller = _controller()
    ids = _chain(controller)
    controller.set_show_atom_numbers(True)
    window = _numbered_window(controller)
    x, y, _ = controller.screen_positions()[ids[0]]
    _right_click(controller, x, y)

    labels = [p for p in window.scene_primitives() if p[0] == "label"]
    assert [label[3] for label in labels] == ["1", "2", "3"]


def test_numbers_are_laid_out_once_per_gesture(monkeypatch):
    controller = _controller()
    _two_molecules(controller)
    controller.set_show_atom_numbers(True)
    _shift_click(controller, 100.0, 300.0)
    window = _numbered_window(controller)
    calls = []
    original = molecularDesign.place_atom_labels
    monkeypatch.setattr(molecularDesign, "place_atom_labels",
                        lambda *args, **kwargs: calls.append(1) or original(*args, **kwargs))
    window.scene_primitives()  # a settled frame lays the numbers out
    calls.clear()

    controller.held = frozenset({"ctrl"})
    controller.on_canvas_press(None, 120.0, 300.0)
    for step in range(1, 6):
        controller.on_canvas_drag(None, 120.0 + 5.0 * step, 300.0)
    assert calls == []                      # offsets reused while dragging
    labels = [p for p in window.scene_primitives() if p[0] == "label"]
    atoms = controller.screen_positions()
    o_id = controller.sketch.ordered_atom_ids()[2]
    assert labels[2][1] - atoms[o_id][0] == pytest.approx(window._label_offsets[1][2][0])
    controller.on_canvas_release(None, 145.0, 300.0)
    assert calls                            # laid out again once it ends


def test_the_window_switches_follow_the_controller():
    controller = _controller()
    window = _headless_window(controller)
    window.replace_switch = type("Switch", (), {"text": "", "value": False})()
    controller.register_view(window)

    controller.select_element("Br")
    controller.set_replace_mode(True)

    assert window.replace_switch.text == "Click an atom to replace it with Br"
    assert window.replace_switch.value is True


def test_labels_are_only_assigned_a_new_text():
    class _Label:
        assignments = 0

        def __init__(self):
            self._text = "same"

        @property
        def text(self):
            return self._text

        @text.setter
        def text(self, value):
            self.assignments += 1
            self._text = value

    label = _Label()
    molecularDesign._set_label_text(label, "same")
    assert label.assignments == 0
    molecularDesign._set_label_text(label, "new")
    assert label.assignments == 1
