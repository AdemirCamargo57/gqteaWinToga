"""Interactive molecular design panel for the Molecular Viewer's Design tab.

The user picks an element from a periodic table and draws a molecule on a 2D
canvas; :mod:`molecularPreOptimizer` turns that drawing into a 3D structure,
which is handed to the viewer for display and saving.

Canvas gestures, all on the left button unless stated:

========================================  ==========================================
Gesture                                   Result
========================================  ==========================================
Click empty space                         Place an atom of the selected element
Press an atom, drag to another, release   Create a single bond between them
Click a bond                              Raise its order, wrapping back to single
Right-click an atom or bond               Delete it
========================================  ==========================================

The panel owns no chemistry of its own: what a click means is decided here, but
*whether it is allowed* comes from :class:`~molecularSketch.MoleculeSketch`.
The split is what lets every rule in this file be tested without a window --
the panel is built by ``_init_state`` plus ``build``, and the tests construct
the first half only (see tests/test_molecularDesign.py).
"""
import asyncio
from typing import Dict, List, Optional, Tuple

import toga
from toga.constants import Baseline
from toga.fonts import Font
from toga.style import Pack
from toga.style.pack import CENTER, COLUMN, LEFT, ROW

from help import AtomicData, HelpGqteaWin
from molecularPreOptimizer import OptimizationResult, optimize_sketch
from molecularSketch import MoleculeSketch, SketchBond

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

BOND_ORDER_NAMES = {1: "single", 2: "double", 3: "triple"}

# Canvas palette. Atom colours match molecularViewer's own, so a molecule looks
# the same in the sketch and in the 3D window.
CANVAS_BACKGROUND = "#fbfbfb"
BOND_COLOR = "#303030"
RUBBER_BAND_COLOR = "#0a84ff"
ATOM_OUTLINE_COLOR = "#303030"
ATOM_LABEL_COLOR = "#101010"
DEFAULT_ATOM_COLOR = "#b0b0b0"
ATOM_CANVAS_COLORS = {
    "H": "#ffffff",
    "C": "#33dddd",
    "N": "#3333ff",
    "O": "#ff3333",
    "F": "#99e699",
    "P": "#ff8000",
    "S": "#dddd00",
    "Cl": "#33dd33",
    "Br": "#996633",
    "I": "#8000cc",
}


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


class MolecularDesignPanel:
    """Builds and drives the Design tab of the Molecular Viewer window."""

    CANVAS_WIDTH = 470
    CANVAS_HEIGHT = 340
    # Pointer movement past this many pixels turns a click into a drag, the
    # same discrimination the 3D window makes between picking and rotating.
    DRAG_THRESHOLD = 4.0
    ATOM_RADIUS = 11.0
    DOUBLE_BOND_OFFSET = 3.0
    UNDO_DEPTH = 50

    ELEMENT_BUTTON_SIZE = 30
    HINT_COLOR = "#666666"

    def __init__(self, viewer=None):
        self._init_state(viewer=viewer)

    def _init_state(self, viewer=None):
        """Every piece of state the interaction rules need, and no widgets.

        Split out from ``__init__`` so the tests can exercise the whole state
        machine without constructing a Toga window.
        """
        self.viewer = viewer
        self.sketch = MoleculeSketch()
        self.selected_element = "C"
        self.undo_stack: List[tuple] = []
        self.last_result: Optional[OptimizationResult] = None

        # Rubber band shown while dragging a bond out of an atom, as
        # (from_x, from_y, to_x, to_y) in canvas pixels.
        self.drag_preview: Optional[Tuple[float, float, float, float]] = None
        self._press_point: Optional[Tuple[float, float]] = None
        self._press_target: Tuple[str, object] = ("empty", None)
        self._dragged = False
        self._optimization_running = False

        # Widgets, populated by build().
        self.canvas = None
        self.status_label = None
        self.selected_element_label = None
        self.optimize_button = None
        self.send_button = None
        self.element_buttons: Dict[str, toga.Button] = {}

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------
    def build(self) -> toga.Box:
        """Assemble the tab: periodic table, canvas, actions, help."""
        tab_box = toga.Box(style=Pack(direction=COLUMN, margin=12))

        self.status_label = toga.Label(
            "Pick an element, then click the canvas to place an atom.",
            style=Pack(font_size=9, color=self.HINT_COLOR, margin=(6, 0, 6, 0)),
        )
        self.selected_element_label = toga.Label(
            f"Selected element: {self.selected_element}",
            style=Pack(font_size=11, font_weight="bold", margin=(4, 0, 6, 0)),
        )

        tab_box.add(self._build_periodic_table())
        tab_box.add(self.selected_element_label)
        tab_box.add(toga.Divider(style=Pack(margin=(6, 0, 6, 0))))
        tab_box.add(self._build_canvas_row())
        tab_box.add(self._build_action_row())
        tab_box.add(self.status_label)
        tab_box.add(self._build_help())

        self._redraw()
        return tab_box

    def _build_periodic_table(self) -> toga.Box:
        table_box = toga.Box(style=Pack(direction=COLUMN))
        table_box.add(
            toga.Label(
                "Click an element to select it:",
                style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 4, 0)),
            )
        )
        for row in PERIODIC_TABLE_LAYOUT:
            row_box = toga.Box(style=Pack(direction=ROW))
            for symbol in row:
                if symbol is None:
                    row_box.add(
                        toga.Box(style=Pack(width=self.ELEMENT_BUTTON_SIZE, height=self.ELEMENT_BUTTON_SIZE))
                    )
                    continue
                button = toga.Button(
                    symbol,
                    on_press=self._make_element_handler(symbol),
                    style=Pack(
                        width=self.ELEMENT_BUTTON_SIZE,
                        height=self.ELEMENT_BUTTON_SIZE,
                        font_size=7,
                        background_color=BLOCK_COLORS[element_block(symbol)],
                    ),
                )
                self.element_buttons[symbol] = button
                row_box.add(button)
            table_box.add(row_box)
        self._highlight_selected_button()
        return table_box

    def _make_element_handler(self, symbol: str):
        def handler(widget):
            self.select_element(symbol)

        return handler

    def _build_canvas_row(self) -> toga.Box:
        canvas_box = toga.Box(style=Pack(direction=COLUMN))
        canvas_box.add(
            toga.Label(
                "Click empty space to place an atom - drag from one atom to another to bond "
                "them - click a bond to raise its order - right-click to delete.",
                style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 4, 0)),
            )
        )
        self.canvas = toga.Canvas(
            on_press=self.on_canvas_press,
            on_drag=self.on_canvas_drag,
            on_release=self.on_canvas_release,
            on_alt_press=self.on_canvas_alt_press,
            style=Pack(
                width=self.CANVAS_WIDTH,
                height=self.CANVAS_HEIGHT,
                background_color=CANVAS_BACKGROUND,
            ),
        )
        canvas_box.add(self.canvas)
        return canvas_box

    def _build_action_row(self) -> toga.Box:
        actions = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(8, 0, 4, 0)))
        self.optimize_button = toga.Button(
            "Pre-optimize geometry",
            on_press=self.pre_optimize,
            style=Pack(flex=1, margin=(0, 6, 0, 0)),
        )
        self.send_button = toga.Button(
            "Send to 3D viewer",
            on_press=self.send_to_viewer,
            style=Pack(flex=1, margin=(0, 6, 0, 0)),
        )
        undo_button = toga.Button("Undo", on_press=self.undo, style=Pack(width=70, margin=(0, 6, 0, 0)))
        clear_button = toga.Button("Clear", on_press=self.clear_canvas, style=Pack(width=70))
        for widget in (self.optimize_button, self.send_button, undo_button, clear_button):
            actions.add(widget)
        return actions

    def _build_help(self) -> toga.Box:
        help_box = toga.Box(style=Pack(direction=COLUMN, margin=(8, 0, 0, 0)))
        help_box.add(
            toga.Label("Help", style=Pack(font_size=11, font_weight="bold", margin=(8, 0, 6, 0)))
        )
        help_box.add(
            toga.MultilineTextInput(
                value=HelpGqteaWin.help_molecular_design,
                readonly=True,
                style=Pack(height=120, flex=1),
            )
        )
        return help_box

    # ------------------------------------------------------------------
    # Element selection
    # ------------------------------------------------------------------
    def select_element(self, symbol: str):
        self.selected_element = symbol
        if self.selected_element_label is not None:
            self.selected_element_label.text = f"Selected element: {symbol}"
        self._highlight_selected_button()
        self.set_status(f"{symbol} selected. Click the canvas to place an atom.")

    def _highlight_selected_button(self):
        for symbol, button in self.element_buttons.items():
            button.style.background_color = (
                SELECTED_ELEMENT_COLOR
                if symbol == self.selected_element
                else BLOCK_COLORS[element_block(symbol)]
            )

    # ------------------------------------------------------------------
    # Canvas interaction
    # ------------------------------------------------------------------
    def on_canvas_press(self, widget, x, y, **kwargs):
        """Classify what is under the pointer once, on press."""
        self._press_point = (x, y)
        self._dragged = False
        self.drag_preview = None

        atom = self.sketch.atom_at(x, y)
        if atom is not None:
            self._press_target = ("atom", atom.atom_id)
            return
        bond = self.sketch.bond_at(x, y)
        if bond is not None:
            self._press_target = ("bond", bond)
            return
        self._press_target = ("empty", None)

    def on_canvas_drag(self, widget, x, y, **kwargs):
        if self._press_point is None:
            return
        press_x, press_y = self._press_point
        if not self._dragged:
            moved = (x - press_x) ** 2 + (y - press_y) ** 2
            if moved <= self.DRAG_THRESHOLD ** 2:
                return
            self._dragged = True

        kind, target = self._press_target
        if kind != "atom":
            return
        anchor = self.sketch.get_atom(target)
        self.drag_preview = (anchor.x, anchor.y, x, y)
        self._redraw()

    def on_canvas_release(self, widget, x, y, **kwargs):
        kind, target = self._press_target
        dragged = self._dragged

        self._press_point = None
        self._press_target = ("empty", None)
        self._dragged = False
        self.drag_preview = None

        if kind == "atom":
            self._release_from_atom(target, x, y, dragged)
        elif kind == "bond":
            if not dragged:
                self._cycle_bond(target)
        elif not dragged:
            self._place_atom(x, y)

        self._redraw()

    def on_canvas_alt_press(self, widget, x, y, **kwargs):
        """Right-click deletes whatever is under the pointer."""
        atom = self.sketch.atom_at(x, y)
        if atom is not None:
            description = self.sketch.describe_atom(atom.atom_id)
            self._push_undo()
            self.sketch.remove_atom(atom.atom_id)
            self.set_status(f"Deleted {description}.")
            self._redraw()
            return

        bond = self.sketch.bond_at(x, y)
        if bond is not None:
            self._push_undo()
            self.sketch.remove_bond_object(bond)
            self.set_status("Bond deleted.")
            self._redraw()
            return

        self.set_status("Right-click an atom or a bond to delete it.")

    # ------------------------------------------------------------------
    # Edits
    # ------------------------------------------------------------------
    def _place_atom(self, x: float, y: float):
        self._push_undo()
        atom = self.sketch.add_atom(self.selected_element, x, y)
        self.set_status(
            f"Placed {self.sketch.describe_atom(atom.atom_id)}. "
            f"{self.sketch.atom_count} atoms, {self.sketch.bond_count} bonds."
        )

    def _release_from_atom(self, atom_id: int, x: float, y: float, dragged: bool):
        if not dragged:
            self.set_status("Drag from this atom to another atom to bond them.")
            return

        other = self.sketch.atom_at(x, y)
        if other is None or other.atom_id == atom_id:
            self.set_status("Release the pointer on another atom to create a bond.")
            return

        if self.sketch.find_bond(atom_id, other.atom_id) is not None:
            self.set_status(
                "Those atoms are already bonded. Click the bond itself to change its order."
            )
            return

        first = self.sketch.describe_atom(atom_id)
        second = self.sketch.describe_atom(other.atom_id)
        self._push_undo()
        self.sketch.add_bond(atom_id, other.atom_id)
        self.set_status(f"Single bond created between {first} and {second}.")

    def _cycle_bond(self, bond: SketchBond):
        element_i = self.sketch.get_atom(bond.atom_i).element
        element_j = self.sketch.get_atom(bond.atom_j).element
        maximum = self.sketch.max_bond_order_of(bond)

        self._push_undo()
        new_order = self.sketch.cycle_bond_order(bond)

        if maximum == 1:
            self.set_status(
                f"{element_i}-{element_j} bonds are limited to a single bond, so it stays single."
            )
            return
        name = BOND_ORDER_NAMES.get(new_order, str(new_order))
        wrapped = " (back to the start of the cycle)" if new_order == 1 else ""
        self.set_status(f"{element_i}-{element_j} bond is now {name}{wrapped}.")

    def _push_undo(self):
        """Snapshot before a change, and drop any geometry that no longer matches."""
        self.undo_stack.append(self.sketch.snapshot())
        if len(self.undo_stack) > self.UNDO_DEPTH:
            self.undo_stack.pop(0)
        self.last_result = None

    def undo(self, widget=None):
        if not self.undo_stack:
            self.set_status("Nothing to undo.")
            return
        self.sketch.restore(self.undo_stack.pop())
        self.last_result = None
        self.set_status(
            f"Undone. {self.sketch.atom_count} atoms, {self.sketch.bond_count} bonds."
        )
        self._redraw()

    def clear_canvas(self, widget=None):
        if self.sketch.is_empty:
            self.set_status("The canvas is already empty.")
            return
        self._push_undo()
        self.sketch.clear()
        self.set_status("Canvas cleared. Undo will bring the structure back.")
        self._redraw()

    # ------------------------------------------------------------------
    # Pre-optimization
    # ------------------------------------------------------------------
    async def pre_optimize(self, widget=None):
        """Validate, optimize off the UI thread, and report the outcome."""
        if self._optimization_running:
            return
        if self.sketch.is_empty:
            await self._show_error(
                "Nothing to optimize", "Draw a structure on the canvas first."
            )
            return

        self._optimization_running = True
        if self.optimize_button is not None:
            self.optimize_button.enabled = False
        self.set_status("Pre-optimizing geometry...")
        try:
            result = await asyncio.to_thread(optimize_sketch, self.sketch)
        except ValueError as exc:
            self.set_status("Optimization refused: the structure is not ready.")
            await self._show_error("Cannot optimize this structure", str(exc))
            return
        except Exception as exc:  # pragma: no cover - defensive, never expected
            self.set_status("Optimization failed.")
            await self._show_error("Optimization failed", str(exc))
            return
        finally:
            self._optimization_running = False
            if self.optimize_button is not None:
                self.optimize_button.enabled = True

        self.apply_optimization_result(result)
        title = "Pre-optimization complete" if result.success else "Pre-optimization did not converge"
        await self._show_info(title, self.describe_result(result))

    def apply_optimization_result(self, result: OptimizationResult):
        """Record a finished optimization and report it in the status line."""
        self.last_result = result
        if result.success:
            self.set_status(
                f"Optimization converged - {result.atom_count} atoms, "
                f"gradient norm {result.gradient_norm:.4f}. Send it to the 3D viewer to see it."
            )
        else:
            self.set_status(
                f"Optimization stopped early - {result.message} "
                "The partially relaxed structure is still available."
            )

    def describe_result(self, result: OptimizationResult) -> str:
        """The multi-line report shown in the dialog after an optimization."""
        lines: List[str] = []
        if result.success:
            lines.append(
                f"Converged for {result.atom_count} atoms in {result.iterations} iterations."
            )
        else:
            lines.append(
                f"Did not converge for {result.atom_count} atoms "
                f"after {result.iterations} iterations."
            )
            lines.append(f"Minimizer message: {result.message}")
            lines.append("The partially relaxed structure is still available.")

        lines.append(f"Final energy: {result.energy:.2f} (relative units)")
        lines.append(f"Gradient norm: {result.gradient_norm:.4f}")

        if result.warnings:
            lines.append("")
            lines.append("Warnings:")
            lines.extend(f"  - {warning}" for warning in result.warnings)

        lines.append("")
        lines.append(
            "This is a pre-optimized starting geometry, not a converged "
            "electronic-structure result."
        )
        return "\n".join(lines)

    def send_to_viewer(self, widget=None) -> bool:
        """Hand the optimized structure to the 3D viewer for display and saving."""
        if self.last_result is None:
            self.set_status(
                "Pre-optimize the structure before sending it to the 3D viewer."
            )
            return False
        if self.viewer is None:
            self.set_status("No viewer is attached to this design panel.")
            return False

        # The sketch still matches last_result: any edit clears it in _push_undo.
        self.viewer.load_designed_molecule(
            self.last_result.to_frame(), bonds=self.sketch.bond_index_pairs()
        )
        self.set_status(
            f"{self.last_result.atom_count} atoms sent to the 3D viewer. Use "
            "'Display Molecule/Trajectory' to see it, or 'Save Current Frame XYZ' to save it."
        )
        return True

    # ------------------------------------------------------------------
    # Reporting helpers
    # ------------------------------------------------------------------
    def set_status(self, message: str):
        if self.status_label is not None:
            self.status_label.text = message

    def _dialog_window(self):
        return getattr(self.viewer, "main_window", None)

    async def _show_error(self, title: str, message: str):
        window = self._dialog_window()
        if window is not None:
            await window.dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        window = self._dialog_window()
        if window is not None:
            await window.dialog(toga.InfoDialog(title, message))

    # ------------------------------------------------------------------
    # Drawing
    # ------------------------------------------------------------------
    def _redraw(self):
        if self.canvas is None:
            return
        self._draw_sketch()
        self.canvas.redraw()

    def _draw_sketch(self):
        context = getattr(self.canvas, "context", None)
        if context is None:  # the stub canvas used in tests draws nothing
            return
        context.clear()

        for bond in self.sketch.bonds:
            self._draw_bond(context, bond)

        if self.drag_preview is not None:
            from_x, from_y, to_x, to_y = self.drag_preview
            with context.Stroke(color=RUBBER_BAND_COLOR, line_width=1.5, line_dash=[4.0, 3.0]) as stroke:
                stroke.move_to(from_x, from_y)
                stroke.line_to(to_x, to_y)

        label_font = Font(family="sans-serif", size=9)
        for atom in self.sketch.atoms:
            color = ATOM_CANVAS_COLORS.get(atom.element, DEFAULT_ATOM_COLOR)
            with context.Fill(color=color) as fill:
                fill.ellipse(atom.x, atom.y, self.ATOM_RADIUS, self.ATOM_RADIUS)
            with context.Stroke(color=ATOM_OUTLINE_COLOR, line_width=1.0) as stroke:
                stroke.ellipse(atom.x, atom.y, self.ATOM_RADIUS, self.ATOM_RADIUS)
            with context.Fill(color=ATOM_LABEL_COLOR) as fill:
                fill.write_text(
                    atom.element,
                    atom.x - 2.0 * len(atom.element) - 1.0,
                    atom.y,
                    label_font,
                    Baseline.MIDDLE,
                )

    def _draw_bond(self, context, bond: SketchBond):
        """Draw a bond as one, two or three parallel lines."""
        start = self.sketch.get_atom(bond.atom_i)
        end = self.sketch.get_atom(bond.atom_j)
        dx = end.x - start.x
        dy = end.y - start.y
        length = (dx * dx + dy * dy) ** 0.5
        if length == 0.0:
            return
        # Unit normal, so parallel lines are offset sideways from the bond axis.
        normal_x = -dy / length
        normal_y = dx / length

        offsets = {
            1: (0.0,),
            2: (-self.DOUBLE_BOND_OFFSET, self.DOUBLE_BOND_OFFSET),
            3: (-2.0 * self.DOUBLE_BOND_OFFSET, 0.0, 2.0 * self.DOUBLE_BOND_OFFSET),
        }.get(bond.order, (0.0,))

        with context.Stroke(color=BOND_COLOR, line_width=1.8) as stroke:
            for offset in offsets:
                stroke.move_to(start.x + normal_x * offset, start.y + normal_y * offset)
                stroke.line_to(end.x + normal_x * offset, end.y + normal_y * offset)
