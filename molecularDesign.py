"""Interactive molecular design: the controller and the design canvas window.

The feature is three windows working on one piece of state:

* the viewer's **Design tab**, a launcher with two buttons and a status line;
* the **Periodic Table** window ([periodicTable.py](periodicTable.py));
* the **Design Canvas** window, defined here.

:class:`MolecularDesignController` owns everything that matters -- the sketch,
the selected element, the undo stack, the last optimization -- and touches no
widget at all. The windows are **disposable views** over it, registered with
``register_view`` and dropped on close.

That split is forced, not stylistic: Toga states that *a closed window cannot
be reused*, so reopening one means constructing a new one. If the sketch lived
in the canvas window, closing that window would destroy the user's drawing.
Keeping state in the controller means a window can be closed and reopened
freely, and it is also what lets a click in the periodic-table window reach the
canvas immediately -- the two windows never reference each other, they only
talk to the controller, which broadcasts to every live view.

Canvas gestures, all on the left button unless stated:

========================================  ==========================================
Gesture                                   Result
========================================  ==========================================
Click empty space                         Place an atom of the selected element
Press an atom, drag to another, release   Create a single bond between them
Click a bond                              Step its type: single, double, triple,
                                          resonance, back to single
Right-click an atom or bond               Delete it
========================================  ==========================================

The controller owns no chemistry: what a gesture *means* is decided here, but
whether it is allowed comes from :class:`~molecularSketch.MoleculeSketch`.
"""
import asyncio
from typing import List, Optional, Tuple

import toga
from toga.constants import Baseline
from toga.fonts import Font
from toga.style import Pack
from toga.style.pack import CENTER, COLUMN, ROW

from help import HelpGqteaWin
from molecularPreOptimizer import OptimizationResult, optimize_sketch
from molecularSketch import RESONANCE_ORDER, MoleculeSketch, SketchBond, bond_type_name
from periodicTable import PeriodicTableWindow

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


class MolecularDesignController:
    """State and rules for molecular design. Never touches a widget."""

    # Pointer movement past this many pixels turns a click into a drag, the
    # same discrimination the 3D window makes between picking and rotating.
    DRAG_THRESHOLD = 4.0
    ATOM_RADIUS = 11.0
    BOND_LINE_OFFSET = 3.0
    UNDO_DEPTH = 50

    def __init__(self, viewer=None):
        self._init_state(viewer=viewer)

    def _init_state(self, viewer=None):
        """All state, no widgets -- the half the tests construct on its own."""
        self.viewer = viewer
        self.sketch = MoleculeSketch()
        self.selected_element = "C"
        self.undo_stack: List[tuple] = []
        self.last_result: Optional[OptimizationResult] = None
        self.status_message = "Pick an element, then click the canvas to place an atom."

        self.views: List[object] = []
        self.periodic_table_window: Optional[PeriodicTableWindow] = None
        self.canvas_window: Optional["DesignCanvasWindow"] = None

        # Rubber band shown while dragging a bond out of an atom, as
        # (from_x, from_y, to_x, to_y) in canvas pixels.
        self.drag_preview: Optional[Tuple[float, float, float, float]] = None
        self._press_point: Optional[Tuple[float, float]] = None
        self._press_target: Tuple[str, object] = ("empty", None)
        self._dragged = False
        self._optimization_running = False

        # Launcher-tab widgets, the only ones the controller holds directly.
        self.launcher_status_label = None
        self.launcher_element_label = None

    # ------------------------------------------------------------------
    # Views
    # ------------------------------------------------------------------
    def register_view(self, view):
        """Attach a window and bring it up to date immediately.

        A view registered later (a window reopened after being closed) is given
        the current selection and told to draw, so it never shows stale state.
        """
        if view in self.views:
            return
        self.views.append(view)
        self._tell(view, "on_element_changed", self.selected_element)
        self._tell(view, "on_status_changed", self.status_message)
        self._tell(view, "on_sketch_changed")

    def unregister_view(self, view):
        if view in self.views:
            self.views.remove(view)

    def _broadcast(self, hook: str, *args):
        for view in list(self.views):
            self._tell(view, hook, *args)

    @staticmethod
    def _tell(view, hook: str, *args):
        handler = getattr(view, hook, None)
        if handler is not None:
            handler(*args)

    # ------------------------------------------------------------------
    # Windows
    # ------------------------------------------------------------------
    def open_periodic_table(self, widget=None):
        """Open the element picker, or re-focus it if it is already open."""
        if self.periodic_table_window is not None and self.periodic_table_window.window is not None:
            self.periodic_table_window.window.show()
            return self.periodic_table_window
        self.periodic_table_window = PeriodicTableWindow(self)
        self.set_status("Periodic table opened. Click an element to select it.")
        return self.periodic_table_window

    def open_design_canvas(self, widget=None):
        """Open the drawing canvas, or re-focus it if it is already open."""
        if self.canvas_window is not None and self.canvas_window.window is not None:
            self.canvas_window.window.show()
            return self.canvas_window
        self.canvas_window = DesignCanvasWindow(self)
        self.set_status(
            f"Design canvas opened. {self.sketch.atom_count} atoms, "
            f"{self.sketch.bond_count} bonds."
        )
        return self.canvas_window

    def close_windows(self):
        """Close both design windows, e.g. when the viewer window closes."""
        for window in (self.periodic_table_window, self.canvas_window):
            if window is not None:
                window.close()
        self.periodic_table_window = None
        self.canvas_window = None

    def build_launcher_tab(self) -> toga.Box:
        """The viewer's Design tab: two buttons, a status line and the help."""
        tab_box = toga.Box(style=Pack(direction=COLUMN, margin=12))
        tab_box.add(
            toga.Label(
                "Molecular design uses two windows: a periodic table to choose an "
                "element, and a canvas to draw on. Leave both open side by side.",
                style=Pack(font_size=9, color="#666666", margin=(0, 0, 8, 0)),
            )
        )

        buttons = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 8, 0)))
        buttons.add(
            toga.Button(
                "Open periodic table",
                on_press=self.open_periodic_table,
                style=Pack(flex=1, margin=(0, 6, 0, 0)),
            )
        )
        buttons.add(
            toga.Button(
                "Open design canvas",
                on_press=self.open_design_canvas,
                style=Pack(flex=1),
            )
        )
        tab_box.add(buttons)

        self.launcher_element_label = toga.Label(
            f"Selected element: {self.selected_element}",
            style=Pack(font_size=11, font_weight="bold", margin=(4, 0, 4, 0)),
        )
        self.launcher_status_label = toga.Label(
            self.status_message,
            style=Pack(font_size=9, color="#666666", margin=(0, 0, 8, 0)),
        )
        tab_box.add(self.launcher_element_label)
        tab_box.add(self.launcher_status_label)

        tab_box.add(
            toga.Label("Help", style=Pack(font_size=11, font_weight="bold", margin=(8, 0, 6, 0)))
        )
        tab_box.add(
            toga.MultilineTextInput(
                value=HelpGqteaWin.help_molecular_design,
                readonly=True,
                style=Pack(height=220, flex=1),
            )
        )

        # The launcher is itself a view, so it tracks the selection and status.
        self.register_view(_LauncherView(self))
        return tab_box

    # ------------------------------------------------------------------
    # Element selection
    # ------------------------------------------------------------------
    def select_element(self, symbol: str):
        """Choose the element the next canvas click will place.

        Broadcast rather than pushed: this is what carries a click in the
        periodic-table window over to the canvas window.
        """
        self.selected_element = symbol
        self._broadcast("on_element_changed", symbol)
        self.set_status(f"{symbol} selected. Click the canvas to place an atom.")

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
            name = bond_type_name(bond.order)
            self._push_undo()
            self.sketch.remove_bond_object(bond)
            self.set_status(f"Deleted a {name} bond.")
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
                "Those atoms are already bonded. Click the bond itself to change its type."
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

        name = bond_type_name(new_order)
        if new_order == RESONANCE_ORDER:
            suffix = " - click again to return it to single"
        elif new_order == 1:
            suffix = " (back to the start of the cycle)"
        else:
            suffix = ""
        self.set_status(f"{element_i}-{element_j} bond is now {name}{suffix}.")

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

    # ------------------------------------------------------------------
    # Handover to the 3D viewer
    # ------------------------------------------------------------------
    def transfer_to_viewer(self) -> Tuple[bool, str]:
        """Load the optimized structure into the viewer. No dialogs, no window.

        Split from :meth:`send_to_viewer` so the transfer itself stays testable
        without an event loop or a GL context.
        """
        if self.last_result is None:
            return False, "Pre-optimize the structure before sending it to the 3D viewer."
        if self.viewer is None:
            return False, "No 3D viewer is attached to this design session."

        try:
            # The sketch still matches last_result: any edit clears it in _push_undo.
            self.viewer.load_designed_molecule(
                self.last_result.to_frame(), bonds=self.sketch.bond_index_pairs()
            )
        except Exception as exc:
            return False, f"The 3D viewer could not accept the structure: {exc}"

        return True, f"{self.last_result.atom_count} atoms transferred to the 3D viewer."

    async def send_to_viewer(self, widget=None):
        """Transfer the structure and open (or re-focus) the 3D viewer window."""
        transferred, message = self.transfer_to_viewer()
        if not transferred:
            self.set_status(message)
            await self._show_error("Cannot send to the 3D viewer", message)
            return False

        self.set_status("Opening the 3D viewer...")
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            self.set_status(f"The 3D viewer could not be opened: {error}")
            await self._show_error(
                "The 3D viewer could not be opened",
                f"{error}\n\nThe structure has been transferred, so you can also use "
                "'Display Molecule/Trajectory' or save it with 'Save Current Frame XYZ'.",
            )
            return False

        atom_count = self.last_result.atom_count if self.last_result else 0
        self.set_status(
            f"{atom_count} atoms shown in the 3D viewer. "
            "Use 'Save Current Frame XYZ' to save the structure."
        )
        return True

    # ------------------------------------------------------------------
    # Reporting helpers
    # ------------------------------------------------------------------
    def set_status(self, message: str):
        self.status_message = message
        self._broadcast("on_status_changed", message)

    def _redraw(self):
        self._broadcast("on_sketch_changed")

    def _dialog_window(self):
        """Prefer the canvas window: that is where the user is working."""
        if self.canvas_window is not None and self.canvas_window.window is not None:
            return self.canvas_window.window
        return getattr(self.viewer, "main_window", None)

    async def _show_error(self, title: str, message: str):
        window = self._dialog_window()
        if window is not None:
            await window.dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        window = self._dialog_window()
        if window is not None:
            await window.dialog(toga.InfoDialog(title, message))


class _LauncherView:
    """Keeps the Design tab's two labels in step with the controller."""

    def __init__(self, controller: MolecularDesignController):
        self.controller = controller

    def on_element_changed(self, symbol: str):
        if self.controller.launcher_element_label is not None:
            self.controller.launcher_element_label.text = f"Selected element: {symbol}"

    def on_status_changed(self, message: str):
        if self.controller.launcher_status_label is not None:
            self.controller.launcher_status_label.text = message


class DesignCanvasWindow:
    """The drawing window: canvas, actions and a status line."""

    TITLE = "Molecular Design Canvas"
    CANVAS_WIDTH = 560
    CANVAS_HEIGHT = 420
    WINDOW_SIZE = (600, 620)
    HINT_COLOR = "#666666"

    def __init__(self, controller: MolecularDesignController):
        self.controller = controller
        self.canvas = None
        self.status_label = toga.Label(
            "", style=Pack(font_size=9, color=self.HINT_COLOR, margin=(6, 0, 0, 0))
        )
        self.element_label = toga.Label(
            "", style=Pack(font_size=11, font_weight="bold", margin=(0, 0, 6, 0))
        )

        self.window = toga.Window(title=self.TITLE, size=self.WINDOW_SIZE, on_close=self._on_close)
        self.window.content = self._build()

        self.controller.register_view(self)
        self.window.show()

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------
    def _build(self) -> toga.Box:
        content = toga.Box(style=Pack(direction=COLUMN, margin=12))
        content.add(self.element_label)
        content.add(
            toga.Label(
                "Click empty space to place an atom - drag from one atom to another to "
                "bond them - click a bond to step single, double, triple, resonance - "
                "right-click to delete.",
                style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 6, 0)),
            )
        )

        self.canvas = toga.Canvas(
            on_press=self.controller.on_canvas_press,
            on_drag=self.controller.on_canvas_drag,
            on_release=self.controller.on_canvas_release,
            on_alt_press=self.controller.on_canvas_alt_press,
            style=Pack(
                width=self.CANVAS_WIDTH,
                height=self.CANVAS_HEIGHT,
                background_color=CANVAS_BACKGROUND,
            ),
        )
        content.add(self.canvas)
        content.add(self._build_actions())
        content.add(self.status_label)
        return content

    def _build_actions(self) -> toga.Box:
        actions = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(8, 0, 0, 0)))
        actions.add(
            toga.Button(
                "Pre-optimize geometry",
                on_press=self.controller.pre_optimize,
                style=Pack(flex=1, margin=(0, 6, 0, 0)),
            )
        )
        actions.add(
            toga.Button(
                "Send to 3D viewer",
                on_press=self.controller.send_to_viewer,
                style=Pack(flex=1, margin=(0, 6, 0, 0)),
            )
        )
        actions.add(
            toga.Button(
                "Periodic table",
                on_press=self.controller.open_periodic_table,
                style=Pack(width=110, margin=(0, 6, 0, 0)),
            )
        )
        actions.add(
            toga.Button("Undo", on_press=self.controller.undo, style=Pack(width=64, margin=(0, 6, 0, 0)))
        )
        actions.add(toga.Button("Clear", on_press=self.controller.clear_canvas, style=Pack(width=64)))
        return actions

    # ------------------------------------------------------------------
    # Controller hooks
    # ------------------------------------------------------------------
    def on_element_changed(self, symbol: str):
        self.element_label.text = f"Selected element: {symbol}"

    def on_status_changed(self, message: str):
        self.status_label.text = message

    def on_sketch_changed(self):
        if self.canvas is None:
            return
        self._draw_sketch()
        self.canvas.redraw()

    # ------------------------------------------------------------------
    # Lifetime
    # ------------------------------------------------------------------
    def _on_close(self, window, **kwargs) -> bool:
        self.controller.unregister_view(self)
        self.controller.canvas_window = None
        return True

    def close(self):
        self.controller.unregister_view(self)
        if self.window is not None:
            self.window.close()
            self.window = None

    # ------------------------------------------------------------------
    # Drawing
    # ------------------------------------------------------------------
    def _draw_sketch(self):
        context = getattr(self.canvas, "context", None)
        if context is None:
            return
        context.clear()

        sketch = self.controller.sketch
        for bond in sketch.bonds:
            self._draw_bond(context, bond)

        if self.controller.drag_preview is not None:
            from_x, from_y, to_x, to_y = self.controller.drag_preview
            with context.Stroke(
                color=RUBBER_BAND_COLOR, line_width=1.5, line_dash=[4.0, 3.0]
            ) as stroke:
                stroke.move_to(from_x, from_y)
                stroke.line_to(to_x, to_y)

        radius = self.controller.ATOM_RADIUS
        label_font = Font(family="sans-serif", size=9)
        for atom in sketch.atoms:
            color = ATOM_CANVAS_COLORS.get(atom.element, DEFAULT_ATOM_COLOR)
            with context.Fill(color=color) as fill:
                fill.ellipse(atom.x, atom.y, radius, radius)
            with context.Stroke(color=ATOM_OUTLINE_COLOR, line_width=1.0) as stroke:
                stroke.ellipse(atom.x, atom.y, radius, radius)
            with context.Fill(color=ATOM_LABEL_COLOR) as fill:
                fill.write_text(
                    atom.element,
                    atom.x - 2.0 * len(atom.element) - 1.0,
                    atom.y,
                    label_font,
                    Baseline.MIDDLE,
                )

    def _draw_bond(self, context, bond: SketchBond):
        """One, two or three parallel lines -- or, for resonance, one solid
        line with a dashed line beside it."""
        sketch = self.controller.sketch
        start = sketch.get_atom(bond.atom_i)
        end = sketch.get_atom(bond.atom_j)
        dx = end.x - start.x
        dy = end.y - start.y
        length = (dx * dx + dy * dy) ** 0.5
        if length == 0.0:
            return
        # Unit normal, so parallel lines are offset sideways from the bond axis.
        normal_x = -dy / length
        normal_y = dx / length
        offset = self.controller.BOND_LINE_OFFSET

        if bond.is_resonance:
            self._stroke_line(context, start, end, normal_x, normal_y, -offset, dashed=False)
            self._stroke_line(context, start, end, normal_x, normal_y, offset, dashed=True)
            return

        offsets = {
            1: (0.0,),
            2: (-offset, offset),
            3: (-2.0 * offset, 0.0, 2.0 * offset),
        }.get(int(bond.order), (0.0,))
        for line_offset in offsets:
            self._stroke_line(context, start, end, normal_x, normal_y, line_offset, dashed=False)

    @staticmethod
    def _stroke_line(context, start, end, normal_x, normal_y, offset, dashed: bool):
        dash = [5.0, 4.0] if dashed else None
        with context.Stroke(color=BOND_COLOR, line_width=1.8, line_dash=dash) as stroke:
            stroke.move_to(start.x + normal_x * offset, start.y + normal_y * offset)
            stroke.line_to(end.x + normal_x * offset, end.y + normal_y * offset)
