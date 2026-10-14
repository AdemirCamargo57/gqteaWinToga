"""Tests for the molecular design controller.

The controller holds every rule and every piece of state and touches no widget,
so these tests construct it with ``__new__`` plus ``_init_state`` and attach
stub views. Nothing here creates a Toga window.

Three things are under test: the canvas state machine (click to place, drag to
bond or to grow a bonded atom, click to cycle a bond type, right-click to
delete), the 3D view (left-drag on empty space pans, right-drag rotates, and
neither is an edit), and the plumbing around it -- broadcasting to several
views at once, surviving a closed window, writing an optimized geometry back
onto the drawing, and handing the finished structure to the 3D viewer.
"""
import asyncio

import numpy as np
import pytest

import molecularDesign
from molecularDesign import MolecularDesignController, blend_color
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
# The 3D view: pan, rotate, zoom -- none of them an edit
# ------------------------------------------------------------------
def test_dragging_empty_space_pans_the_view():
    controller = _controller()
    _click(controller, 100.0, 100.0)

    _drag(controller, 300.0, 300.0, 340.0, 280.0)

    sx, sy, _ = controller.screen_positions()[controller.sketch.atoms[0].atom_id]
    assert (sx, sy) == pytest.approx((140.0, 80.0))


def test_panning_moves_the_view_not_the_atoms():
    controller = _controller()
    _two_bonded_carbons(controller)
    before = controller.sketch.positions()
    undo_depth = len(controller.undo_stack)

    _drag(controller, 300.0, 300.0, 340.0, 280.0)

    assert np.array_equal(controller.sketch.positions(), before)
    assert len(controller.undo_stack) == undo_depth


def test_right_dragging_rotates_the_view():
    controller = _controller()
    _two_bonded_carbons(controller)

    _right_drag(controller, 300.0, 300.0, 380.0, 300.0)

    assert not np.allclose(controller.camera.rotation, np.eye(3))


def test_right_dragging_over_an_atom_rotates_instead_of_deleting():
    controller = _controller()
    _two_bonded_carbons(controller)

    _right_drag(controller, 100.0, 100.0, 180.0, 160.0)

    assert controller.sketch.atom_count == 2
    assert controller.sketch.bond_count == 1


def test_rotating_is_not_an_edit():
    """A finished optimization still matches the drawing after a rotation."""
    controller = _controller()
    _two_bonded_carbons(controller)
    controller.apply_optimization_result(_successful_result())
    before = controller.sketch.positions()
    undo_depth = len(controller.undo_stack)

    _right_drag(controller, 300.0, 300.0, 380.0, 340.0)

    assert np.array_equal(controller.sketch.positions(), before)
    assert len(controller.undo_stack) == undo_depth
    assert controller.last_result is not None


def test_rotation_turns_about_the_molecule_centre():
    controller = _controller()
    _two_bonded_carbons(controller)  # centred on (150, 100)

    _right_drag(controller, 300.0, 300.0, 390.0, 250.0)

    centre = controller.sketch.positions().mean(axis=0)
    projected = controller.camera.project(centre)[0]
    assert projected[:2] == pytest.approx((150.0, 100.0))


def test_a_bond_can_still_be_clicked_after_rotating():
    controller = _controller()
    bond = _two_bonded_carbons(controller)
    _right_drag(controller, 300.0, 300.0, 340.0, 330.0)
    a, b = (controller.screen_positions()[atom_id] for atom_id in (bond.atom_i, bond.atom_j))

    _click(controller, (a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)

    assert bond.order == 2


def test_a_click_after_rotating_places_the_atom_on_the_plane_through_the_centre():
    controller = _controller()
    _two_bonded_carbons(controller)
    _right_drag(controller, 300.0, 300.0, 380.0, 330.0)
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
    _right_drag(controller, 300.0, 300.0, 457.0, 300.0)  # about 90 degrees

    _drag(controller, 280.0, 210.0, 450.0, 210.0)

    assert abs(controller.sketch.atoms[1].z) > 10.0


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

    _right_drag(controller, 100.0, 300.0, 300.0, 200.0)

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
