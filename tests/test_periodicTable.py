"""Tests for the periodic-table element picker.

The layout data is pure and tested directly; the window half is tested through
a stub controller, so no Toga window is created.
"""
import pytest

from help import AtomicData
from periodicTable import (
    BLOCK_COLORS,
    PERIODIC_TABLE_LAYOUT,
    PeriodicTableWindow,
    element_block,
    table_symbols,
)


class _StubController:
    def __init__(self):
        self.selected_element = "C"
        self.selected_calls = []
        self.registered = []
        self.unregistered = []

    def select_element(self, symbol):
        self.selected_element = symbol
        self.selected_calls.append(symbol)

    def register_view(self, view):
        self.registered.append(view)

    def unregister_view(self, view):
        self.unregistered.append(view)


class _StubButton:
    def __init__(self, text):
        self.text = text
        self.style = type("Style", (), {"background_color": None})()


class _StubLabel:
    def __init__(self, text=""):
        self.text = text


def _headless_window(controller=None) -> PeriodicTableWindow:
    """A PeriodicTableWindow with its state but no Toga window."""
    controller = controller or _StubController()
    view = PeriodicTableWindow.__new__(PeriodicTableWindow)
    view.controller = controller
    view.window = None
    view.element_buttons = {symbol: _StubButton(symbol) for symbol in table_symbols()}
    view.selected_label = _StubLabel()
    return view


# ------------------------------------------------------------------
# Layout data
# ------------------------------------------------------------------
def test_the_periodic_table_lists_all_118_elements():
    assert len(table_symbols()) == 118


def test_every_periodic_table_symbol_is_a_real_element():
    unknown = [s for s in table_symbols() if s not in AtomicData.atomic_numbers]

    assert unknown == []


def test_the_periodic_table_has_no_duplicate_symbols():
    symbols = table_symbols()

    assert len(set(symbols)) == len(symbols)


def test_the_periodic_table_rows_are_all_eighteen_columns_wide():
    assert all(len(row) == 18 for row in PERIODIC_TABLE_LAYOUT)


def test_hydrogen_and_helium_sit_at_the_ends_of_the_first_row():
    assert PERIODIC_TABLE_LAYOUT[0][0] == "H"
    assert PERIODIC_TABLE_LAYOUT[0][17] == "He"


def test_elements_are_classified_into_blocks():
    assert element_block("Na") == "s"
    assert element_block("O") == "p"
    assert element_block("Fe") == "d"
    assert element_block("U") == "f"


def test_every_element_has_a_block_colour():
    assert all(element_block(symbol) in BLOCK_COLORS for symbol in table_symbols())


# ------------------------------------------------------------------
# Selecting an element
# ------------------------------------------------------------------
def test_pressing_an_element_button_tells_the_controller():
    controller = _StubController()
    view = _headless_window(controller)

    view.element_pressed("N")

    assert controller.selected_calls == ["N"]


def test_the_window_follows_a_selection_made_elsewhere():
    """The controller drives the label, so a canvas-side change is reflected."""
    view = _headless_window()

    view.on_element_changed("Fe")

    assert "Fe" in view.selected_label.text


def test_the_selected_button_is_highlighted_and_the_others_are_not():
    view = _headless_window()

    view.on_element_changed("O")

    assert view.element_buttons["O"].style.background_color != BLOCK_COLORS["p"]
    assert view.element_buttons["C"].style.background_color == BLOCK_COLORS["p"]


def test_changing_the_selection_releases_the_previous_button():
    view = _headless_window()

    view.on_element_changed("O")
    view.on_element_changed("N")

    assert view.element_buttons["O"].style.background_color == BLOCK_COLORS["p"]
    assert view.element_buttons["N"].style.background_color != BLOCK_COLORS["p"]
