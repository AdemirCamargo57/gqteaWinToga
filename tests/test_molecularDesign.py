"""Tests for the molecular design panel's interaction rules.

The panel is constructed with ``__new__`` plus ``_init_state`` and given stub
widgets, so no Toga window is ever created -- the same headless convention
tests/test_molecularViewer.py uses for MolecularViewerUI. What is under test is
the canvas state machine (click to place, drag to bond, click to raise a bond
order, right-click to delete) and the reporting around pre-optimization.
"""
import numpy as np
import pytest

from help import AtomicData
from molecularDesign import PERIODIC_TABLE_LAYOUT, MolecularDesignPanel
from molecularPreOptimizer import OptimizationResult


class _StubCanvas:
    def __init__(self):
        self.redraw_count = 0

    def redraw(self):
        self.redraw_count += 1


class _StubLabel:
    def __init__(self, text=""):
        self.text = text


class _StubViewer:
    def __init__(self):
        self.loaded_frame = None
        self.loaded_bonds = None

    def load_designed_molecule(self, frame, bonds=None):
        self.loaded_frame = frame
        self.loaded_bonds = bonds


def _panel(viewer=None) -> MolecularDesignPanel:
    panel = MolecularDesignPanel.__new__(MolecularDesignPanel)
    panel._init_state(viewer=viewer)
    panel.canvas = _StubCanvas()
    panel.status_label = _StubLabel()
    panel.selected_element_label = _StubLabel()
    return panel


def _click(panel, x, y):
    """A press and release at the same point: no drag."""
    panel.on_canvas_press(None, x, y)
    panel.on_canvas_release(None, x, y)


def _drag(panel, from_x, from_y, to_x, to_y):
    panel.on_canvas_press(None, from_x, from_y)
    panel.on_canvas_drag(None, (from_x + to_x) / 2.0, (from_y + to_y) / 2.0)
    panel.on_canvas_drag(None, to_x, to_y)
    panel.on_canvas_release(None, to_x, to_y)


def _two_bonded_carbons(panel):
    _click(panel, 100.0, 100.0)
    _click(panel, 200.0, 100.0)
    _drag(panel, 100.0, 100.0, 200.0, 100.0)
    return panel.sketch.bonds[0]


# ------------------------------------------------------------------
# Element selection
# ------------------------------------------------------------------
def test_carbon_is_selected_by_default():
    assert _panel().selected_element == "C"


def test_selecting_an_element_changes_what_a_click_places():
    panel = _panel()

    panel.select_element("N")
    _click(panel, 50.0, 50.0)

    assert panel.selected_element == "N"
    assert panel.sketch.atoms[0].element == "N"


def test_selecting_an_element_is_reported_in_the_label():
    panel = _panel()

    panel.select_element("Fe")

    assert "Fe" in panel.selected_element_label.text


# ------------------------------------------------------------------
# Placing atoms
# ------------------------------------------------------------------
def test_clicking_empty_space_places_an_atom():
    panel = _panel()

    _click(panel, 120.0, 80.0)

    assert panel.sketch.atom_count == 1
    assert (panel.sketch.atoms[0].x, panel.sketch.atoms[0].y) == (120.0, 80.0)


def test_dragging_across_empty_space_places_nothing():
    panel = _panel()

    _drag(panel, 20.0, 20.0, 200.0, 200.0)

    assert panel.sketch.atom_count == 0


def test_placing_an_atom_redraws_the_canvas():
    panel = _panel()

    _click(panel, 120.0, 80.0)

    assert panel.canvas.redraw_count > 0


def test_clicking_an_existing_atom_places_nothing_on_top_of_it():
    panel = _panel()
    _click(panel, 100.0, 100.0)

    _click(panel, 102.0, 101.0)

    assert panel.sketch.atom_count == 1


# ------------------------------------------------------------------
# Drawing bonds
# ------------------------------------------------------------------
def test_dragging_from_one_atom_to_another_creates_a_single_bond():
    panel = _panel()
    _click(panel, 100.0, 100.0)
    _click(panel, 200.0, 100.0)

    _drag(panel, 100.0, 100.0, 200.0, 100.0)

    assert panel.sketch.bond_count == 1
    assert panel.sketch.bonds[0].order == 1


def test_dragging_from_an_atom_back_to_itself_creates_no_bond():
    panel = _panel()
    _click(panel, 100.0, 100.0)

    _drag(panel, 100.0, 100.0, 103.0, 103.0)

    assert panel.sketch.bond_count == 0


def test_dragging_from_an_atom_to_empty_space_creates_no_bond():
    panel = _panel()
    _click(panel, 100.0, 100.0)
    _click(panel, 200.0, 100.0)

    _drag(panel, 100.0, 100.0, 400.0, 400.0)

    assert panel.sketch.bond_count == 0


def test_dragging_between_two_already_bonded_atoms_leaves_the_order_alone():
    panel = _panel()
    bond = _two_bonded_carbons(panel)
    bond.order = 2

    _drag(panel, 100.0, 100.0, 200.0, 100.0)

    assert panel.sketch.bond_count == 1
    assert bond.order == 2


def test_a_rubber_band_is_tracked_while_dragging_from_an_atom():
    panel = _panel()
    _click(panel, 100.0, 100.0)

    panel.on_canvas_press(None, 100.0, 100.0)
    panel.on_canvas_drag(None, 160.0, 140.0)

    assert panel.drag_preview == (100.0, 100.0, 160.0, 140.0)


def test_the_rubber_band_is_cleared_on_release():
    panel = _panel()
    _click(panel, 100.0, 100.0)

    _drag(panel, 100.0, 100.0, 300.0, 300.0)

    assert panel.drag_preview is None


# ------------------------------------------------------------------
# Bond order
# ------------------------------------------------------------------
def test_clicking_a_bond_raises_its_order():
    panel = _panel()
    bond = _two_bonded_carbons(panel)

    _click(panel, 150.0, 100.0)

    assert bond.order == 2


def test_clicking_a_bond_at_maximum_order_returns_it_to_single():
    panel = _panel()
    bond = _two_bonded_carbons(panel)

    _click(panel, 150.0, 100.0)
    _click(panel, 150.0, 100.0)
    _click(panel, 150.0, 100.0)

    assert bond.order == 1


def test_clicking_a_bond_reports_the_new_order():
    panel = _panel()
    _two_bonded_carbons(panel)

    _click(panel, 150.0, 100.0)

    assert "double" in panel.status_label.text.lower()


def test_clicking_a_capped_bond_says_why_it_cannot_rise():
    panel = _panel()
    panel.select_element("O")
    _click(panel, 100.0, 100.0)
    panel.select_element("H")
    _click(panel, 200.0, 100.0)
    _drag(panel, 100.0, 100.0, 200.0, 100.0)

    _click(panel, 150.0, 100.0)

    assert panel.sketch.bonds[0].order == 1
    assert "single" in panel.status_label.text.lower()


def test_clicking_a_bond_never_places_an_atom():
    panel = _panel()
    _two_bonded_carbons(panel)

    _click(panel, 150.0, 100.0)

    assert panel.sketch.atom_count == 2


# ------------------------------------------------------------------
# Deleting
# ------------------------------------------------------------------
def test_right_clicking_an_atom_deletes_it_and_its_bonds():
    panel = _panel()
    _two_bonded_carbons(panel)

    panel.on_canvas_alt_press(None, 100.0, 100.0)

    assert panel.sketch.atom_count == 1
    assert panel.sketch.bond_count == 0


def test_right_clicking_a_bond_deletes_only_the_bond():
    panel = _panel()
    _two_bonded_carbons(panel)

    panel.on_canvas_alt_press(None, 150.0, 100.0)

    assert panel.sketch.atom_count == 2
    assert panel.sketch.bond_count == 0


def test_right_clicking_empty_space_does_nothing():
    panel = _panel()
    _two_bonded_carbons(panel)

    panel.on_canvas_alt_press(None, 500.0, 500.0)

    assert panel.sketch.atom_count == 2
    assert panel.sketch.bond_count == 1


# ------------------------------------------------------------------
# Undo and clear
# ------------------------------------------------------------------
def test_undo_removes_the_last_placed_atom():
    panel = _panel()
    _click(panel, 100.0, 100.0)
    _click(panel, 200.0, 100.0)

    panel.undo(None)

    assert panel.sketch.atom_count == 1


def test_undo_restores_a_deleted_atom():
    panel = _panel()
    _two_bonded_carbons(panel)

    panel.on_canvas_alt_press(None, 100.0, 100.0)
    panel.undo(None)

    assert panel.sketch.atom_count == 2
    assert panel.sketch.bond_count == 1


def test_undo_restores_a_bond_order():
    panel = _panel()
    _two_bonded_carbons(panel)
    _click(panel, 150.0, 100.0)

    panel.undo(None)

    assert panel.sketch.bonds[0].order == 1


def test_undo_on_an_empty_history_is_harmless():
    panel = _panel()

    panel.undo(None)

    assert panel.sketch.atom_count == 0


def test_clearing_the_canvas_empties_the_sketch():
    panel = _panel()
    _two_bonded_carbons(panel)

    panel.clear_canvas(None)

    assert panel.sketch.atom_count == 0
    assert panel.sketch.bond_count == 0


def test_clearing_the_canvas_can_be_undone():
    panel = _panel()
    _two_bonded_carbons(panel)

    panel.clear_canvas(None)
    panel.undo(None)

    assert panel.sketch.atom_count == 2


# ------------------------------------------------------------------
# Pre-optimization reporting
# ------------------------------------------------------------------
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


def test_a_successful_result_reports_energy_and_atom_count():
    summary = _panel().describe_result(_successful_result())

    assert "3 atoms" in summary
    assert "-1.25" in summary
    assert "17" in summary


def test_a_failed_result_reports_the_minimizer_message():
    result = _successful_result()
    result.success = False
    result.message = "STOP: TOTAL NO. of ITERATIONS REACHED LIMIT"

    summary = _panel().describe_result(result)

    assert "ITERATIONS REACHED LIMIT" in summary


def test_validation_warnings_appear_in_the_summary():
    result = _successful_result()
    result.warnings = ["Unusual valence on atom 1 (C): 5 bonds where 4 is typical."]

    summary = _panel().describe_result(result)

    assert "Unusual valence" in summary


def test_applying_a_result_stores_it_and_reports_success():
    panel = _panel()

    panel.apply_optimization_result(_successful_result())

    assert panel.last_result is not None
    assert "converged" in panel.status_label.text.lower()


def test_sending_to_the_viewer_hands_over_the_optimized_frame():
    viewer = _StubViewer()
    panel = _panel(viewer=viewer)
    panel.apply_optimization_result(_successful_result())

    panel.send_to_viewer(None)

    assert viewer.loaded_frame is not None
    assert viewer.loaded_frame[0][0] == "O"


def test_sending_to_the_viewer_before_optimizing_reports_why_it_cannot():
    viewer = _StubViewer()
    panel = _panel(viewer=viewer)

    panel.send_to_viewer(None)

    assert viewer.loaded_frame is None
    assert "optimize" in panel.status_label.text.lower()


def test_editing_the_sketch_invalidates_a_previous_result():
    """A structure that no longer matches the drawing must not be sent on."""
    panel = _panel()
    panel.apply_optimization_result(_successful_result())

    _click(panel, 300.0, 300.0)

    assert panel.last_result is None


# ------------------------------------------------------------------
# The periodic table
# ------------------------------------------------------------------
def test_the_periodic_table_lists_all_118_elements():
    symbols = [symbol for row in PERIODIC_TABLE_LAYOUT for symbol in row if symbol]

    assert len(symbols) == 118


def test_every_periodic_table_symbol_is_a_real_element():
    symbols = [symbol for row in PERIODIC_TABLE_LAYOUT for symbol in row if symbol]

    unknown = [symbol for symbol in symbols if symbol not in AtomicData.atomic_numbers]
    assert unknown == []


def test_the_periodic_table_has_no_duplicate_symbols():
    symbols = [symbol for row in PERIODIC_TABLE_LAYOUT for symbol in row if symbol]

    assert len(set(symbols)) == len(symbols)


def test_the_periodic_table_rows_are_all_eighteen_columns_wide():
    assert all(len(row) == 18 for row in PERIODIC_TABLE_LAYOUT)


def test_hydrogen_and_helium_sit_at_the_ends_of_the_first_row():
    assert PERIODIC_TABLE_LAYOUT[0][0] == "H"
    assert PERIODIC_TABLE_LAYOUT[0][17] == "He"


def test_sending_to_the_viewer_also_hands_over_the_drawn_connectivity():
    """The viewer must not have to re-guess bonds it was already told about."""
    viewer = _StubViewer()
    panel = _panel(viewer=viewer)
    _two_bonded_carbons(panel)
    panel.apply_optimization_result(_successful_result())

    panel.send_to_viewer(None)

    assert viewer.loaded_bonds == [(0, 1)]
