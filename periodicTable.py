"""The periodic-table element picker, in its own window.

Split out of molecularDesign so the design canvas and the element picker can be
two independent windows, and because an element picker is useful on its own:
nothing here knows about sketches, bonds or geometry. It reports a chosen
symbol to a controller and, in the other direction, follows a selection made
anywhere else through :meth:`PeriodicTableWindow.on_element_changed`.

That second direction is the point of the split: the picker and the canvas
never reference each other, they only talk to the controller, so a click here
reaches the canvas window immediately without the two windows being coupled.
"""
from typing import AbstractSet, Dict, Optional

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, LEFT, ROW

# The layout data lives in the Toga-free periodicTableData; re-exported here
# for the window and for everything that imported it from this module.
from periodicTableData import (  # noqa: F401
    BLOCK_COLORS,
    PERIODIC_TABLE_LAYOUT,
    SELECTED_ELEMENT_COLOR,
    element_block,
    group_and_period,
    table_symbols,
)

DEFAULT_HINT = ("Click an element to select it. The design canvas starts using it "
                "immediately, so you can leave this window open while you draw.")
MARKED_LEGEND = "Bold: a tested case is in the USPP library."
BUTTON_SIZE = 32


def button_style_kwargs(symbol: str, marked: Optional[AbstractSet[str]]) -> dict:
    """The Pack style of one element button; marked elements are bold."""
    style = {"width": BUTTON_SIZE, "height": BUTTON_SIZE, "font_size": 7,
             "background_color": BLOCK_COLORS[element_block(symbol)]}
    if marked and symbol in marked:
        style["font_weight"] = "bold"
    return style


def button_enabled(symbol: str, selectable: Optional[AbstractSet[str]]) -> bool:
    """Every element is pressable unless a selectable set leaves it out."""
    return selectable is None or symbol in selectable


class PeriodicTableWindow:
    """A window of 118 element buttons that drives a design controller.

    ``hint`` replaces the instruction line, ``marked`` draws those elements in
    bold with a legend, ``selectable`` disables every element outside it. With
    none of them the window is the Design tab's picker, unchanged.
    """

    TITLE = "Periodic Table"
    WINDOW_SIZE = (660, 480)
    BUTTON_SIZE = BUTTON_SIZE
    HINT_COLOR = "#666666"

    def __init__(self, controller, *, hint: Optional[str] = None,
                 marked: Optional[AbstractSet[str]] = None,
                 selectable: Optional[AbstractSet[str]] = None):
        self.controller = controller
        self.hint = DEFAULT_HINT if hint is None else hint
        self.marked = marked
        self.selectable = selectable
        self.element_buttons: Dict[str, toga.Button] = {}
        self.selected_label = toga.Label(
            "", style=Pack(font_size=11, font_weight="bold", margin=(8, 0, 4, 0))
        )

        self.window = toga.Window(
            title=self.TITLE,
            size=self.WINDOW_SIZE,
            on_close=self._on_close,
        )
        self.window.content = self._build()

        # Registering hands us the current selection straight away, so the
        # window opens already showing whatever the canvas is set to.
        self.controller.register_view(self)
        self.window.show()

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------
    def _build(self) -> toga.Box:
        content = toga.Box(style=Pack(direction=COLUMN, margin=12))
        content.add(
            toga.Label(
                self.hint,
                style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 8, 0)),
            )
        )

        for row in PERIODIC_TABLE_LAYOUT:
            row_box = toga.Box(style=Pack(direction=ROW))
            for symbol in row:
                if symbol is None:
                    row_box.add(
                        toga.Box(style=Pack(width=self.BUTTON_SIZE, height=self.BUTTON_SIZE))
                    )
                    continue
                button = toga.Button(
                    symbol,
                    on_press=self._make_handler(symbol),
                    enabled=button_enabled(symbol, self.selectable),
                    style=Pack(**button_style_kwargs(symbol, self.marked)),
                )
                self.element_buttons[symbol] = button
                row_box.add(button)
            content.add(row_box)

        if self.marked:
            content.add(
                toga.Label(
                    MARKED_LEGEND,
                    style=Pack(font_size=9, color=self.HINT_COLOR, margin=(6, 0, 0, 0)),
                )
            )
        content.add(self.selected_label)
        return content

    def _make_handler(self, symbol: str):
        def handler(widget):
            self.element_pressed(symbol)

        return handler

    # ------------------------------------------------------------------
    # Controller conversation
    # ------------------------------------------------------------------
    def element_pressed(self, symbol: str):
        """A button was pressed: tell the controller and let it broadcast back."""
        self.controller.select_element(symbol)

    def on_element_changed(self, symbol: str):
        """Follow the selection, whoever made it."""
        self.selected_label.text = f"Selected element: {symbol}"
        for button_symbol, button in self.element_buttons.items():
            button.style.background_color = (
                SELECTED_ELEMENT_COLOR
                if button_symbol == symbol
                else BLOCK_COLORS[element_block(button_symbol)]
            )

    # ------------------------------------------------------------------
    # Lifetime
    # ------------------------------------------------------------------
    def _on_close(self, window, **kwargs) -> bool:
        """The user closed the window: stop receiving updates."""
        self.controller.unregister_view(self)
        self.controller.periodic_table_window = None
        return True

    def close(self):
        """Close programmatically. Toga does not fire on_close for this."""
        self.controller.unregister_view(self)
        if self.window is not None:
            self.window.close()
            self.window = None
