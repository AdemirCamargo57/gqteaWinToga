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

The canvas is **three-dimensional**: every atom has a depth, a
:class:`~molecularCamera.ViewCamera` projects the drawing onto the screen, and
the view can be rotated, panned and zoomed while building. The view is never
an edit -- it pushes no undo entry and leaves a finished optimization valid.

Canvas gestures, all on the left button unless stated:

========================================  ==========================================
Gesture                                   Result
========================================  ==========================================
Click empty space                         Place an atom on the screen-facing plane
                                          through the molecule's centre
Press an atom, drag to another, release   Create a single bond between them
Press an atom, drag to empty space        Grow a new atom bonded to it, at the
                                          anchor's depth and the reference length
Click a bond                              Step its type: single, double, triple,
                                          resonance, back to single
Drag empty space                          Pan the view
Right-drag                                Rotate the view about the molecule
Right-click (no drag) an atom or bond     Delete it
========================================  ==========================================

Pre-optimizing writes the optimized geometry **back onto the drawing** (one
undo step), so the canvas always shows the structure that would be sent to
the 3D viewer, and the next edit builds on the relaxed geometry.

The controller owns no chemistry: what a gesture *means* is decided here, but
whether it is allowed comes from :class:`~molecularSketch.MoleculeSketch`.
"""
import asyncio
from typing import Dict, List, Optional, Tuple

import numpy as np
import toga
from toga.constants import Baseline
from toga.fonts import Font
from toga.style import Pack
from toga.style.pack import CENTER, COLUMN, ROW

from help import HelpGqteaWin
from molecularCamera import ViewCamera
from molecularPreOptimizer import (
    PIXELS_PER_ANGSTROM,
    OptimizationResult,
    optimize_sketch,
    reference_bond_length,
)
from molecularSketch import RESONANCE_ORDER, MoleculeSketch, SketchBond, bond_type_name
from periodicTable import PeriodicTableWindow

# Canvas size in pixels; the view zooms about its centre.
CANVAS_WIDTH = 560
CANVAS_HEIGHT = 420

# Length, in Angstrom, of a bond grown from an atom when either element has no
# covalent radius to derive a reference length from.
DEFAULT_GROWN_BOND_LENGTH = 1.5

ZOOM_STEP = 1.25

# Depth cue: the furthest atom is drawn this much smaller and this far faded
# towards the background; the nearest is drawn as is.
DEPTH_SHRINK = 0.25
DEPTH_FADE = 0.55

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


def blend_color(color: str, background: str, fraction: float) -> str:
    """Mix ``fraction`` of ``background`` into ``color`` (both ``#rrggbb``)."""
    fraction = min(1.0, max(0.0, float(fraction)))
    mixed = []
    for start in (1, 3, 5):
        front = int(color[start:start + 2], 16)
        back = int(background[start:start + 2], 16)
        mixed.append(int(round(front + (back - front) * fraction)))
    return "#" + "".join(f"{channel:02x}" for channel in mixed)


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
        self._last_drag_point: Optional[Tuple[float, float]] = None
        self._optimization_running = False

        # The 3D view, and the right-button gesture that rotates it.
        self.camera = ViewCamera()
        self.view_centre = (CANVAS_WIDTH / 2.0, CANVAS_HEIGHT / 2.0)
        self._alt_press_point: Optional[Tuple[float, float]] = None
        self._alt_last_point: Optional[Tuple[float, float]] = None
        self._alt_dragged = False

        # Bumped by every edit, so an optimization that finishes after the
        # user has changed the drawing can tell its result is stale.
        self._edit_generation = 0

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
    # The 3D view
    # ------------------------------------------------------------------
    def screen_positions(self) -> Dict[int, Tuple[float, float, float]]:
        """Every atom's projected ``(screen_x, screen_y, depth)``, by atom id."""
        projected = self.camera.project(self.sketch.positions())
        return {
            atom_id: (float(sx), float(sy), float(depth))
            for atom_id, (sx, sy, depth) in zip(self.sketch.ordered_atom_ids(), projected)
        }

    def drawing_order(self) -> List[Tuple[str, object, float]]:
        """What to paint, back to front, as ``(kind, item, fade)`` triples.

        ``kind`` is ``"bond"`` or ``"atom"``; ``fade`` runs from 0 for the
        nearest depth to 1 for the furthest, and is 0 throughout a flat
        drawing. A bond sits at the mean depth of its ends, so it is painted
        over the atom behind it and under the atom in front -- and at equal
        depth bonds go first, so a flat drawing looks exactly as it always did.
        """
        screen = self.screen_positions()
        if not screen:
            return []
        depths = [depth for _, _, depth in screen.values()]
        nearest = min(depths)
        span = max(depths) - nearest

        def fade(depth: float) -> float:
            return 0.0 if span < 1e-9 else (depth - nearest) / span

        items = []
        for bond in self.sketch.bonds:
            depth = (screen[bond.atom_i][2] + screen[bond.atom_j][2]) / 2.0
            items.append((depth, 0, "bond", bond))
        for atom in self.sketch.atoms:
            items.append((screen[atom.atom_id][2], 1, "atom", atom))
        items.sort(key=lambda item: (-item[0], item[1]))
        return [(kind, item, fade(depth)) for depth, _, kind, item in items]

    def zoom_in(self, widget=None):
        self._zoom(ZOOM_STEP)

    def zoom_out(self, widget=None):
        self._zoom(1.0 / ZOOM_STEP)

    def _zoom(self, factor: float):
        self.camera.zoom_by(factor, *self.view_centre)
        self.set_status(f"Zoom {self.camera.zoom * 100.0:.0f}%.")
        self._redraw()

    def reset_view(self, widget=None):
        self.camera.reset()
        self.set_status("View reset: unrotated, unpanned, 100% zoom.")
        self._redraw()

    def _placement_depth(self) -> float:
        """Depth of the screen-facing plane through the molecule's centre."""
        if self.sketch.is_empty:
            return 0.0
        return self.camera.depth_of(self.sketch.positions().mean(axis=0))

    # ------------------------------------------------------------------
    # Canvas interaction
    # ------------------------------------------------------------------
    def on_canvas_press(self, widget, x, y, **kwargs):
        """Classify what is under the pointer once, on press."""
        self._press_point = (x, y)
        self._last_drag_point = (x, y)
        self._dragged = False
        self.drag_preview = None

        screen = self.screen_positions()
        atom = self.sketch.atom_at(x, y, screen)
        if atom is not None:
            self._press_target = ("atom", atom.atom_id)
            return
        bond = self.sketch.bond_at(x, y, screen)
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
        if kind == "atom":
            anchor_x, anchor_y, _ = self.screen_positions()[target]
            self.drag_preview = (anchor_x, anchor_y, x, y)
            self._redraw()
        elif kind == "empty":
            # Panning starts from the press point, so the movement made while
            # still inside the drag threshold is not lost.
            last_x, last_y = self._last_drag_point
            self.camera.pan_by(x - last_x, y - last_y)
            self._last_drag_point = (x, y)
            self._redraw()

    def on_canvas_release(self, widget, x, y, **kwargs):
        kind, target = self._press_target
        dragged = self._dragged

        self._press_point = None
        self._last_drag_point = None
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
        """Start a right-button gesture: a drag rotates, a click deletes."""
        self._alt_press_point = (x, y)
        self._alt_last_point = (x, y)
        self._alt_dragged = False

    def on_canvas_alt_drag(self, widget, x, y, **kwargs):
        if self._alt_press_point is None:
            return
        if not self._alt_dragged:
            press_x, press_y = self._alt_press_point
            if (x - press_x) ** 2 + (y - press_y) ** 2 <= self.DRAG_THRESHOLD ** 2:
                return
            self._alt_dragged = True
            # Turn about the molecule's centre, wherever it is on screen.
            if not self.sketch.is_empty:
                self.camera.set_pivot(self.sketch.positions().mean(axis=0))

        last_x, last_y = self._alt_last_point
        self.camera.rotate_by_drag(x - last_x, y - last_y)
        self._alt_last_point = (x, y)
        self._redraw()

    def on_canvas_alt_release(self, widget, x, y, **kwargs):
        press = self._alt_press_point
        dragged = self._alt_dragged
        self._alt_press_point = None
        self._alt_last_point = None
        self._alt_dragged = False

        if press is None:
            return
        if dragged:
            self.set_status("View rotated. Reset view returns to the original orientation.")
            return
        self._delete_at(*press)

    def _delete_at(self, x: float, y: float):
        """Delete whatever is under a right-click."""
        screen = self.screen_positions()
        atom = self.sketch.atom_at(x, y, screen)
        if atom is not None:
            description = self.sketch.describe_atom(atom.atom_id)
            self._push_undo()
            self.sketch.remove_atom(atom.atom_id)
            self.set_status(f"Deleted {description}.")
            self._redraw()
            return

        bond = self.sketch.bond_at(x, y, screen)
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
        position = self.camera.unproject(x, y, self._placement_depth())
        self._push_undo()
        atom = self.sketch.add_atom(self.selected_element, *position)
        self.set_status(
            f"Placed {self.sketch.describe_atom(atom.atom_id)}. "
            f"{self.sketch.atom_count} atoms, {self.sketch.bond_count} bonds."
        )

    def _release_from_atom(self, atom_id: int, x: float, y: float, dragged: bool):
        if not dragged:
            self.set_status(
                "Drag from this atom to another atom to bond them, "
                "or to empty space to grow a new bonded atom."
            )
            return

        other = self.sketch.atom_at(x, y, self.screen_positions())
        if other is None:
            self._grow_atom(atom_id, x, y)
            return
        if other.atom_id == atom_id:
            self.set_status(
                "Release the pointer on another atom to bond them, "
                "or on empty space to grow a new bonded atom."
            )
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

    def _grow_atom(self, anchor_id: int, x: float, y: float):
        """Add an atom bonded to ``anchor_id``, towards the release point.

        It lands at the anchor's depth, at the single-bond reference length
        the pre-optimizer would relax to, so a structure grown this way is
        already close to its optimized geometry. Rotating the view first is
        how a branch is grown out of the drawing plane.
        """
        anchor = self.sketch.get_atom(anchor_id)
        anchor_x, anchor_y, depth = self.screen_positions()[anchor_id]
        dx, dy = x - anchor_x, y - anchor_y
        distance = (dx * dx + dy * dy) ** 0.5
        if distance == 0.0:
            return

        length = reference_bond_length(anchor.element, self.selected_element)
        if length is None:
            length = DEFAULT_GROWN_BOND_LENGTH
        on_screen = length * PIXELS_PER_ANGSTROM * self.camera.zoom
        position = self.camera.unproject(
            anchor_x + dx / distance * on_screen, anchor_y + dy / distance * on_screen, depth
        )

        self._push_undo()
        atom = self.sketch.add_atom(self.selected_element, *position)
        self.sketch.add_bond(anchor_id, atom.atom_id)
        self.set_status(
            f"Placed {self.sketch.describe_atom(atom.atom_id)} bonded to "
            f"{self.sketch.describe_atom(anchor_id)}."
        )

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
        self._edit_generation += 1

    def undo(self, widget=None):
        if not self.undo_stack:
            self.set_status("Nothing to undo.")
            return
        self.sketch.restore(self.undo_stack.pop())
        self.last_result = None
        self._edit_generation += 1
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
        # Optimize a copy: the canvas stays live while the minimizer runs on a
        # worker thread, and an edit made meanwhile must not race with it.
        generation = self._edit_generation
        try:
            result = await asyncio.to_thread(optimize_sketch, self.sketch.copy())
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

        if generation != self._edit_generation:
            self.set_status(
                "The drawing changed while it was being optimized, so that result was "
                "discarded. Pre-optimize again to relax the current structure."
            )
            return

        self.apply_optimization_result(result)
        title = "Pre-optimization complete" if result.success else "Pre-optimization did not converge"
        await self._show_info(title, self.describe_result(result))

    def apply_optimization_result(self, result: OptimizationResult):
        """Adopt a finished optimization: write it onto the drawing and report it.

        The optimized geometry replaces the drawn one as a single undo step, so
        Undo brings the hand-drawn layout back. ``last_result`` is set *after*
        that push, since every push clears it. A partially relaxed result is
        written back too -- it is still closer to a sensible structure than the
        drawing was.
        """
        positions = result.sketch_positions
        written_back = positions is not None and np.shape(positions) == (self.sketch.atom_count, 3)
        if written_back:
            self._push_undo()
            self.sketch.update_positions(positions)
        self.last_result = result

        on_canvas = (
            " The canvas shows the optimized geometry (Undo restores your drawing)."
            if written_back
            else ""
        )
        if result.success:
            self.set_status(
                f"Optimization converged - {result.atom_count} atoms, "
                f"gradient norm {result.gradient_norm:.4f}.{on_canvas} "
                "Send it to the 3D viewer to see it."
            )
        else:
            self.set_status(
                f"Optimization stopped early - {result.message} "
                f"The partially relaxed structure is still available.{on_canvas}"
            )
        self._redraw()

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
            "The canvas now shows this geometry; Undo restores your drawing. "
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
    """The drawing window: 3D canvas, edit and view actions, and a status line."""

    TITLE = "Molecular Design Canvas"
    CANVAS_WIDTH = CANVAS_WIDTH
    CANVAS_HEIGHT = CANVAS_HEIGHT
    WINDOW_SIZE = (600, 680)
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
                "Click empty space to place an atom - drag atom to atom to bond them, or "
                "atom to empty space to grow a bonded atom - click a bond to change its type.",
                style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 2, 0)),
            )
        )
        content.add(
            toga.Label(
                "Right-drag to rotate - drag empty space to pan - right-click to delete.",
                style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 6, 0)),
            )
        )

        self.canvas = toga.Canvas(
            on_press=self.controller.on_canvas_press,
            on_drag=self.controller.on_canvas_drag,
            on_release=self.controller.on_canvas_release,
            on_alt_press=self.controller.on_canvas_alt_press,
            on_alt_drag=self.controller.on_canvas_alt_drag,
            on_alt_release=self.controller.on_canvas_alt_release,
            style=Pack(
                width=self.CANVAS_WIDTH,
                height=self.CANVAS_HEIGHT,
                background_color=CANVAS_BACKGROUND,
            ),
        )
        content.add(self.canvas)
        content.add(self._build_actions())
        content.add(self._build_view_actions())
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

    def _build_view_actions(self) -> toga.Box:
        """Zoom and reset: the view controls the mouse cannot reach on its own."""
        view = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(6, 0, 0, 0)))
        view.add(toga.Label("View:", style=Pack(width=44, color=self.HINT_COLOR)))
        view.add(
            toga.Button(
                "Zoom +", on_press=self.controller.zoom_in, style=Pack(width=80, margin=(0, 6, 0, 0))
            )
        )
        view.add(
            toga.Button(
                "Zoom -", on_press=self.controller.zoom_out, style=Pack(width=80, margin=(0, 6, 0, 0))
            )
        )
        view.add(toga.Button("Reset view", on_press=self.controller.reset_view, style=Pack(width=100)))
        return view

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
        """Paint the projected drawing back to front, with a depth cue."""
        context = getattr(self.canvas, "context", None)
        if context is None:
            return
        context.clear()

        screen = self.controller.screen_positions()
        order = self.controller.drawing_order()
        radii = {
            item.atom_id: self._atom_radius(fade) for kind, item, fade in order if kind == "atom"
        }
        label_font = Font(family="sans-serif", size=9)
        preview_drawn = False
        for kind, item, fade in order:
            if kind == "bond":
                self._draw_bond(context, item, screen, radii, fade)
                continue
            # The rubber band sits over every bond behind the first atom and
            # under the atoms, as it did on the flat canvas.
            if not preview_drawn:
                self._draw_drag_preview(context)
                preview_drawn = True
            self._draw_atom(context, item, screen[item.atom_id], fade, label_font)
        if not preview_drawn:
            self._draw_drag_preview(context)

    def _draw_drag_preview(self, context):
        if self.controller.drag_preview is None:
            return
        from_x, from_y, to_x, to_y = self.controller.drag_preview
        with context.Stroke(color=RUBBER_BAND_COLOR, line_width=1.5, line_dash=[4.0, 3.0]) as stroke:
            stroke.move_to(from_x, from_y)
            stroke.line_to(to_x, to_y)

    def _atom_radius(self, fade: float) -> float:
        return self.controller.ATOM_RADIUS * (1.0 - DEPTH_SHRINK * fade)

    def _draw_atom(self, context, atom, position, fade: float, label_font):
        """A disc that shrinks and fades towards the background with depth."""
        sx, sy, _ = position
        radius = self._atom_radius(fade)
        tint = DEPTH_FADE * fade
        color = ATOM_CANVAS_COLORS.get(atom.element, DEFAULT_ATOM_COLOR)
        with context.Fill(color=blend_color(color, CANVAS_BACKGROUND, tint)) as fill:
            fill.ellipse(sx, sy, radius, radius)
        outline = blend_color(ATOM_OUTLINE_COLOR, CANVAS_BACKGROUND, tint)
        with context.Stroke(color=outline, line_width=1.0) as stroke:
            stroke.ellipse(sx, sy, radius, radius)
        with context.Fill(color=blend_color(ATOM_LABEL_COLOR, CANVAS_BACKGROUND, tint)) as fill:
            fill.write_text(
                atom.element,
                sx - 2.0 * len(atom.element) - 1.0,
                sy,
                label_font,
                Baseline.MIDDLE,
            )

    def _draw_bond(self, context, bond: SketchBond, screen, radii, fade: float):
        """One, two or three parallel lines -- or, for resonance, one solid
        line with a dashed line beside it -- offset sideways on screen.

        The lines run from disc edge to disc edge rather than centre to centre:
        a bond shares an end with each of its atoms, so no painting order can
        stop a centre-to-centre line crossing the disc of the atom it is
        painted after.
        """
        centre_i = screen[bond.atom_i][:2]
        centre_j = screen[bond.atom_j][:2]
        dx = centre_j[0] - centre_i[0]
        dy = centre_j[1] - centre_i[1]
        length = (dx * dx + dy * dy) ** 0.5
        radius_i = radii.get(bond.atom_i, self.controller.ATOM_RADIUS)
        radius_j = radii.get(bond.atom_j, self.controller.ATOM_RADIUS)
        if length <= radius_i + radius_j:
            return  # seen (nearly) end-on: the atoms hide it
        unit_x, unit_y = dx / length, dy / length
        start = (centre_i[0] + unit_x * radius_i, centre_i[1] + unit_y * radius_i)
        end = (centre_j[0] - unit_x * radius_j, centre_j[1] - unit_y * radius_j)
        # Unit normal, so parallel lines are offset sideways from the bond axis.
        normal_x = -unit_y
        normal_y = unit_x
        offset = self.controller.BOND_LINE_OFFSET
        color = blend_color(BOND_COLOR, CANVAS_BACKGROUND, DEPTH_FADE * fade)

        if bond.is_resonance:
            self._stroke_line(context, start, end, normal_x, normal_y, -offset, False, color)
            self._stroke_line(context, start, end, normal_x, normal_y, offset, True, color)
            return

        offsets = {
            1: (0.0,),
            2: (-offset, offset),
            3: (-2.0 * offset, 0.0, 2.0 * offset),
        }.get(int(bond.order), (0.0,))
        for line_offset in offsets:
            self._stroke_line(context, start, end, normal_x, normal_y, line_offset, False, color)

    @staticmethod
    def _stroke_line(context, start, end, normal_x, normal_y, offset, dashed: bool, color: str):
        dash = [5.0, 4.0] if dashed else None
        with context.Stroke(color=color, line_width=1.8, line_dash=dash) as stroke:
            stroke.move_to(start[0] + normal_x * offset, start[1] + normal_y * offset)
            stroke.line_to(end[0] + normal_x * offset, end[1] + normal_y * offset)
