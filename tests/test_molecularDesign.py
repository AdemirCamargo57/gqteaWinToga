"""Tests for the molecular design controller.

The controller holds every rule and every piece of state and touches no widget,
so these tests construct it with ``__new__`` plus ``_init_state`` and attach
stub views. Nothing here creates a Toga window.

Two things are under test: the canvas state machine (click to place, drag to
bond, click to cycle a bond type, right-click to delete) and the plumbing
around it -- broadcasting to several views at once, surviving a closed window,
and handing the finished structure to the 3D viewer.
"""
import asyncio

import numpy as np
import pytest

from molecularDesign import MolecularDesignController
from molecularPreOptimizer import OptimizationResult
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


def test_dragging_from_an_atom_to_empty_space_creates_no_bond():
    controller = _controller()
    _click(controller, 100.0, 100.0)
    _click(controller, 200.0, 100.0)

    _drag(controller, 100.0, 100.0, 400.0, 400.0)

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

    controller.on_canvas_alt_press(None, 150.0, 100.0)

    assert controller.sketch.bond_count == 0


# ------------------------------------------------------------------
# Deleting
# ------------------------------------------------------------------
def test_right_clicking_an_atom_deletes_it_and_its_bonds():
    controller = _controller()
    _two_bonded_carbons(controller)

    controller.on_canvas_alt_press(None, 100.0, 100.0)

    assert controller.sketch.atom_count == 1
    assert controller.sketch.bond_count == 0


def test_right_clicking_a_bond_deletes_only_the_bond():
    controller = _controller()
    _two_bonded_carbons(controller)

    controller.on_canvas_alt_press(None, 150.0, 100.0)

    assert controller.sketch.atom_count == 2
    assert controller.sketch.bond_count == 0


def test_right_clicking_empty_space_does_nothing():
    controller = _controller()
    _two_bonded_carbons(controller)

    controller.on_canvas_alt_press(None, 500.0, 500.0)

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

    controller.on_canvas_alt_press(None, 100.0, 100.0)
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
