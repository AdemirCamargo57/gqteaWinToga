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
from typing import Dict, List, Optional

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, LEFT, ROW

# The periodic table as it is drawn: 18 columns per row, ``None`` for a gap.
# The two f-block rows sit below a blank spacer row, the way a printed table
# lays them out. A test asserts this contains all 118 elements exactly once.
PERIODIC_TABLE_LAYOUT: List[List[Optional[str]]] = [
    ["H", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, "He"],
    ["Li", "Be", None, None, None, None, None, None, None, None, None, None, "B", "C", "N", "O", "F", "Ne"],
    ["Na", "Mg", None, None, None, None, None, None, None, None, None, None, "Al", "Si", "P", "S", "Cl", "Ar"],
    ["K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"],
    ["Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"],
    ["Cs", "Ba", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"],
    ["Fr", "Ra", "Ac", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"],
    [None] * 18,
    [None, None, None, "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", None],
    [None, None, None, "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", None],
]

# Block colours, so the table reads as a periodic table rather than a wall of
# identical squares.
BLOCK_COLORS = {
    "s": "#ffd9d9",
    "p": "#d9e8ff",
    "d": "#e2f5df",
    "f": "#f3e2ff",
}
SELECTED_ELEMENT_COLOR = "#ffe066"


def table_symbols() -> List[str]:
    """Every element in the layout, in reading order."""
    return [symbol for row in PERIODIC_TABLE_LAYOUT for symbol in row if symbol]


def element_block(symbol: str) -> str:
    """Which block of the periodic table an element belongs to."""
    for row_index, row in enumerate(PERIODIC_TABLE_LAYOUT):
        if symbol not in row:
            continue
        column = row.index(symbol)
        if row_index >= 8:
            return "f"
        if column <= 1 or symbol == "He":
            return "s"
        if column >= 12:
            return "p"
        return "d"
    return "p"


class PeriodicTableWindow:
    """A window of 118 element buttons that drives a design controller."""

    TITLE = "Periodic Table"
    WINDOW_SIZE = (660, 480)
    BUTTON_SIZE = 32
    HINT_COLOR = "#666666"

    def __init__(self, controller):
        self.controller = controller
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
                "Click an element to select it. The design canvas starts using it "
                "immediately, so you can leave this window open while you draw.",
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
                    style=Pack(
                        width=self.BUTTON_SIZE,
                        height=self.BUTTON_SIZE,
                        font_size=7,
                        background_color=BLOCK_COLORS[element_block(symbol)],
                    ),
                )
                self.element_buttons[symbol] = button
                row_box.add(button)
            content.add(row_box)

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
