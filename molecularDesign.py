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
Click an atom                              Pick it for a measurement (#1-#4, in order);
                                          click again to unpick. In replace mode:
                                          give it the selected element instead
Click a bond                              Step its type: single, double, triple,
                                          resonance, back to single
Drag empty space                          Pan the view
Right-drag                                Rotate the view about the molecule
Right-click (no drag) an atom or bond     Delete it
Shift+click an atom or bond               Select its whole molecule
Shift+click empty space                   Clear the selection
Ctrl+drag                                 Move the selected molecule rigidly
Ctrl+right-drag                           Turn the selected molecule about its centre
========================================  ==========================================

Toga reports no modifier keys, so Shift/Ctrl are read from WinForms at event
time (:func:`read_keyboard_modifiers`, injectable as ``modifier_state``).
Picks feed the measurement, the Edit value dialog and Freeze; lengths,
angles and dihedrals are defined once in
[molecularGeometry.py](molecularGeometry.py), hydrogens are planned by
[molecularHydrogens.py](molecularHydrogens.py).

Pre-optimizing writes the optimized geometry **back onto the drawing** (one
undo step), so the canvas always shows the structure that would be sent to
the 3D viewer, and the next edit builds on the relaxed geometry.

**Import structure...** replaces the drawing with a structure read from a file
(one undo step), and **Save Frame** writes the canvas in the format chosen
beside it -- XYZ by default. What each format can hold, and every reader and
writer, live in [molecularFileFormats.py](molecularFileFormats.py), which knows
nothing about the canvas; only the pixel/Angstrom conversion happens here.

The controller owns no chemistry: what a gesture *means* is decided here, but
whether it is allowed comes from :class:`~molecularSketch.MoleculeSketch`.
"""
import asyncio
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

import numpy as np
import toga
from toga.constants import Baseline
from toga.fonts import Font
from toga.style import Pack
from toga.style.pack import CENTER, COLUMN, ROW

from help import AtomicData, HelpGqteaWin
from molecularCamera import ViewCamera
from molecularFileFormats import (
    FORMATS_BY_KEY,
    READABLE_FORMATS,
    WRITABLE_FORMATS,
    StructureData,
    export_notes,
    format_for_path,
    read_structure_file,
    structure_to_text,
    write_xyz_text,
)
from molecularGeometry import (
    ANGLE,
    DIHEDRAL,
    DISTANCE,
    KIND_LABELS,
    LINEAR_TOLERANCE_DEGREES,
    UNITS,
    GeometryEditError,
    format_value,
    is_degenerate,
    kind_for,
    measure,
    moving_atoms,
    parse_parameter_value,
    set_parameter,
    wrap_degrees,
)
from molecularCanvasRender import pil_to_net_bitmap, place_atom_labels, render_scene, text_box_size
from molecularHydrogens import plan_hydrogens
from molecularPreOptimizer import (
    PIXELS_PER_ANGSTROM,
    OptimizationResult,
    angstrom_to_sketch,
    optimize_sketch,
    reference_bond_length,
    sketch_to_angstrom,
)
from molecularSketch import RESONANCE_ORDER, MoleculeSketch, SketchBond, bond_type_name, format_order
from periodicTable import PeriodicTableWindow

# Canvas size in pixels; the view zooms about its centre.
CANVAS_WIDTH = 560
CANVAS_HEIGHT = 420

# The canvas is for building molecules; a larger system belongs in the 3D
# viewer, which can open the same file directly.
MAX_IMPORT_ATOMS = 1000

# Format the Save Frame button writes until the user picks another.
DEFAULT_EXPORT_FORMAT = "XYZ"

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


# The selected molecule's halo, and the measurement picks with their order tags.
SELECTION_COLOR = "#ff9500"
PICK_COLOR = "#d6336c"

# Up to four picked atoms: a length, an angle or a dihedral.
MAX_PICKS = 4

# Background presets for the canvas, first = default.
CANVAS_BACKGROUNDS = {
    "Off-white (default)": CANVAS_BACKGROUND,
    "White": "#ffffff",
    "Light grey": "#d9d9d9",
    "Pale blue": "#e8f0fb",
    "Cream": "#fff8e7",
    "Dark grey": "#3a3a3a",
    "Black": "#000000",
}
_HEX_COLOR_RE = re.compile(r"^#?([0-9a-fA-F]{6})$")


def parse_hex_color(text: str) -> str:
    """'#a0b1c2' or 'a0b1c2' -> '#a0b1c2'; ``ValueError`` for anything else."""
    match = _HEX_COLOR_RE.match((text or "").strip())
    if not match:
        raise ValueError(f"'{(text or '').strip()}' is not a colour: use #rrggbb, e.g. #ffffff.")
    return "#" + match.group(1).lower()


def is_dark(color: str) -> bool:
    """Relative luminance below one half (sRGB weights, no gamma: good enough here)."""
    red, green, blue = (int(color[k:k + 2], 16) for k in (1, 3, 5))
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue < 128.0


def canvas_palette(background: str) -> Dict[str, str]:
    """Bond, outline and overlay-text colours that stay visible on ``background``."""
    if is_dark(background):
        return {"bond": "#d8d8d8", "outline": "#d8d8d8", "text": "#f0f0f0", "number": "#9cc3ff"}
    return {"bond": BOND_COLOR, "outline": ATOM_OUTLINE_COLOR, "text": ATOM_LABEL_COLOR, "number": "#1f5fbf"}


def read_keyboard_modifiers() -> FrozenSet[str]:
    """The modifier keys held right now: a subset of {'shift', 'ctrl'}.

    Toga's canvas events carry only x and y, so the modifiers are read from
    WinForms at the moment of the event. Anywhere that is not available
    (another backend, a headless test) no modifier is ever reported, and
    every gesture behaves as it did before modifiers existed.
    """
    try:
        from System.Windows.Forms import Control, Keys
    except Exception:
        return frozenset()
    state = int(Control.ModifierKeys)
    held = set()
    if state & int(Keys.Shift):
        held.add("shift")
    if state & int(Keys.Control):
        held.add("ctrl")
    return frozenset(held)


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

        # File format Save Frame writes (a molecularFileFormats key).
        self.export_format = DEFAULT_EXPORT_FORMAT

        # Which modifier keys are held; injectable so tests never read the
        # real keyboard.
        self.modifier_state = read_keyboard_modifiers
        # The molecule chosen with Shift+click (atom ids), and the atoms
        # picked in order for a measurement. Neither is an edit.
        self.selected_atoms: Set[int] = set()
        self.picked_atoms: List[int] = []
        # The Ctrl gesture under way: "translate" (left) or "rotate" (right).
        self._transform: Optional[str] = None
        self._transform_started = False
        self._alt_transform = False
        # View setting: the canvas background colour.
        self.canvas_background = CANVAS_BACKGROUND
        # While on, clicking an atom replaces its element instead of picking it.
        self.replace_mode = False
        # View setting: each atom's number drawn beside it.
        self.show_atom_numbers = False
        # The last structure handed to the 3D viewer, with its bond orders, so
        # bringing it back restores them (the viewer itself keeps none).
        self._sent_structure: Optional[StructureData] = None
        self.parameter_window: Optional["ParameterEditWindow"] = None

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

    def open_parameter_editor(self, widget=None):
        """Open the dialog that sets the picked parameter to a typed value."""
        try:
            info = self.parameter_edit_info()
        except ValueError as exc:
            self.set_status(str(exc))
            asyncio.ensure_future(self._show_error("Cannot edit this parameter", str(exc)))
            return None
        if self.parameter_window is not None and self.parameter_window.window is not None:
            self.parameter_window.close()
        self.parameter_window = ParameterEditWindow(self, info)
        return self.parameter_window

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
        """Close the design windows, e.g. when the viewer window closes."""
        for window in (self.periodic_table_window, self.canvas_window, self.parameter_window):
            if window is not None:
                window.close()
        self.periodic_table_window = None
        self.canvas_window = None
        self.parameter_window = None

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
        tab_box.add(
            toga.Button(
                "Bring the 3D molecule to the design canvas",
                on_press=self.import_from_viewer,
                style=Pack(margin=(0, 0, 8, 0)),
            )
        )

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
        if self.replace_mode:
            self.set_status(f"{symbol} selected. Click an atom to replace it with {symbol}, "
                            "or empty space to place one.")
        else:
            self.set_status(f"{symbol} selected. Click the canvas to place an atom.")

    # ------------------------------------------------------------------
    # Replacing an atom's element
    # ------------------------------------------------------------------
    def set_replace_mode(self, enabled: bool):
        """While on, clicking an atom gives it the selected element instead of
        picking it for a measurement. Not an edit; every other gesture is unchanged."""
        self.replace_mode = bool(enabled)
        self._broadcast("on_replace_mode_changed", self.replace_mode)
        if self.replace_mode:
            self.set_status(f"Replace mode: click an atom to make it {self.selected_element}. "
                            "Choose another element in the periodic table at any time.")
        else:
            self.set_status("Replace mode off: clicking an atom picks it for measuring again.")

    def replace_atom(self, atom_id: int) -> Tuple[bool, str]:
        """Give an atom the selected element, keeping its position and bonds.

        One undo step. The atom keeps its id, so picks, the selection and
        frozen parameters that name it stay attached. Anything unusual about
        the result -- a valence above the typical one, a bond order the new
        element does not usually form, no covalent radius -- is reported, and
        the bonds are kept as drawn rather than changed behind the user's back.
        """
        element = self.selected_element
        atom = self.sketch.get_atom(atom_id)
        if atom.element == element:
            message = f"{self.sketch.describe_atom(atom_id).capitalize()} is already {element}."
            self.set_status(message)
            return False, message
        self._push_undo()
        previous = self.sketch.replace_element(atom_id, element)
        bonds = len(self.sketch.neighbours(atom_id))
        message = (f"Atom {self.sketch.atom_number(atom_id)} replaced: {previous} -> {element}; "
                   f"its position and {bonds} bond{'s' if bonds != 1 else ''} are kept.")
        notes = self.replacement_notes(atom_id)
        if notes:
            message += " Note: " + " ".join(notes)
        self.set_status(message + " Undo restores it.")
        self._redraw()
        return True, message

    def replacement_notes(self, atom_id: int) -> List[str]:
        """What is unusual about an atom after a replacement (empty when nothing is)."""
        atom = self.sketch.get_atom(atom_id)
        element = atom.element
        notes = []
        if element not in AtomicData.covalent_radii:
            notes.append(f"{element} has no covalent radius, so the structure cannot be pre-optimized.")
        limit = AtomicData.max_valence.get(element)
        total = self.sketch.bond_order_sum(atom_id)
        if limit is not None and total > limit:
            noun = "bond" if total == 1 else "bonds"
            notes.append(f"{element} now has {format_order(total)} {noun} where {limit} is typical.")
        for bond in self.sketch.bonds:
            if not bond.connects(atom_id):
                continue
            cap = self.sketch.max_bond_order_of(bond)
            if bond.order > cap:
                other = self.sketch.atom_number(bond.other_end(atom_id))
                notes.append(f"The {bond_type_name(bond.order)} bond to atom {other} is above the "
                             f"{bond_type_name(cap)} bond {element} usually forms; click it to change it.")
        return notes

    # ------------------------------------------------------------------
    # Atom numbers
    # ------------------------------------------------------------------
    def set_show_atom_numbers(self, enabled: bool):
        """Show each atom's number beside it. A view setting: not an edit."""
        self.show_atom_numbers = bool(enabled)
        self.set_status("Atom numbers shown." if self.show_atom_numbers else "Atom numbers hidden.")
        self._redraw()

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

    def drawing_order(self, screen=None) -> List[Tuple[str, object, float]]:
        """What to paint, back to front, as ``(kind, item, fade)`` triples.

        ``kind`` is ``"bond"`` or ``"atom"``; ``fade`` runs from 0 for the
        nearest depth to 1 for the furthest, and is 0 throughout a flat
        drawing. A bond sits at the mean depth of its ends, so it is painted
        over the atom behind it and under the atom in front -- and at equal
        depth bonds go first, so a flat drawing looks exactly as it always did.
        ``screen`` may pass in :meth:`screen_positions` already computed.
        """
        screen = self.screen_positions() if screen is None else screen
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
    def _modifiers(self) -> FrozenSet[str]:
        try:
            return frozenset(self.modifier_state())
        except Exception:  # a keyboard query must never break a click
            return frozenset()

    def on_canvas_press(self, widget, x, y, **kwargs):
        """Classify what is under the pointer once, on press.

        Ctrl starts moving the selected molecule and Shift selects one; with
        no modifier the gesture is whatever is under the pointer.
        """
        self._press_point = (x, y)
        self._last_drag_point = (x, y)
        self._dragged = False
        self.drag_preview = None

        modifiers = self._modifiers()
        if "ctrl" in modifiers:
            self._press_target = ("translate", None)
            self._transform_started = False
            return
        screen = self.screen_positions()
        if "shift" in modifiers:
            atom = self.sketch.atom_at(x, y, screen)
            bond = None if atom is not None else self.sketch.bond_at(x, y, screen)
            target = atom.atom_id if atom is not None else (bond.atom_i if bond is not None else None)
            self._press_target = ("select", target)
            return
        atom = self.sketch.atom_at(x, y, screen)
        if atom is not None:
            self._press_target = ("atom", atom.atom_id)
            return
        bond = self.sketch.bond_at(x, y, screen)
        if bond is not None:
            self._press_target = ("bond", bond)
            return
        self._press_target = ("empty", None)

    @property
    def gesture_active(self) -> bool:
        """A mouse button is down on the canvas. Views skip work that cannot
        change mid-gesture (the side panel) until it ends."""
        return self._press_point is not None or self._alt_press_point is not None

    def on_canvas_drag(self, widget, x, y, **kwargs):
        if self._press_point is None:
            return
        press_x, press_y = self._press_point
        kind, target = self._press_target
        if not self._dragged:
            moved = (x - press_x) ** 2 + (y - press_y) ** 2
            # A Ctrl drag has no click meaning to protect, so it follows the
            # pointer from the first pixel instead of after the threshold.
            if kind != "translate" and moved <= self.DRAG_THRESHOLD ** 2:
                return
            self._dragged = True

        if kind == "translate":
            last_x, last_y = self._last_drag_point
            self._last_drag_point = (x, y)
            self._translate_selection(x - last_x, y - last_y)
        elif kind == "atom":
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

        if kind == "translate":
            self._finish_transform("Molecule moved")
        elif kind == "select":
            if not dragged:
                self.select_molecule(target)
        elif kind == "atom":
            self._release_from_atom(target, x, y, dragged)
        elif kind == "bond":
            if not dragged:
                self._cycle_bond(target)
        elif not dragged:
            self._place_atom(x, y)

        self._redraw()

    def on_canvas_alt_press(self, widget, x, y, **kwargs):
        """Start a right-button gesture: a drag rotates, a click deletes.

        With Ctrl held, a drag turns the selected molecule instead of the view.
        """
        self._alt_press_point = (x, y)
        self._alt_last_point = (x, y)
        self._alt_dragged = False
        self._alt_transform = "ctrl" in self._modifiers()
        self._transform_started = False

    def on_canvas_alt_drag(self, widget, x, y, **kwargs):
        if self._alt_press_point is None:
            return
        if not self._alt_dragged:
            press_x, press_y = self._alt_press_point
            within = (x - press_x) ** 2 + (y - press_y) ** 2 <= self.DRAG_THRESHOLD ** 2
            # The threshold tells a delete-click from a view rotation; a Ctrl
            # turn never deletes, so it starts at once.
            if within and not self._alt_transform:
                return
            self._alt_dragged = True
            # Turn about the molecule's centre, wherever it is on screen.
            if not self._alt_transform and not self.sketch.is_empty:
                self.camera.set_pivot(self.sketch.positions().mean(axis=0))

        last_x, last_y = self._alt_last_point
        self._alt_last_point = (x, y)
        if self._alt_transform:
            self._rotate_selection(x - last_x, y - last_y)
            return
        self.camera.rotate_by_drag(x - last_x, y - last_y)
        self._redraw()

    def on_canvas_alt_release(self, widget, x, y, **kwargs):
        press = self._alt_press_point
        dragged = self._alt_dragged
        transform = self._alt_transform
        self._alt_press_point = None
        self._alt_last_point = None
        self._alt_dragged = False
        self._alt_transform = False

        if press is None:
            return
        if transform:
            # Ctrl+right never deletes: a short Ctrl+right-click does nothing.
            if dragged:
                self._finish_transform("Molecule rotated")
            return
        if dragged:
            self.set_status("View rotated. Reset view returns to the original orientation.")
            self._redraw()  # the gesture is over: views catch up what they skipped (label layout)
            return
        self._delete_at(*press)

    # ------------------------------------------------------------------
    # Selecting and moving one molecule
    # ------------------------------------------------------------------
    def select_molecule(self, atom_id: Optional[int]):
        """Select the whole molecule containing ``atom_id``; ``None`` clears the selection."""
        if atom_id is None or atom_id not in self.sketch.ordered_atom_ids():
            had = bool(self.selected_atoms)
            self.selected_atoms = set()
            self.set_status("Selection cleared." if had else
                            "Shift+click an atom or a bond to select its molecule.")
            self._redraw()
            return
        self.selected_atoms = set(self.sketch.fragment_containing(atom_id))
        count = len(self.sketch.fragments())
        self.set_status(
            f"Molecule of {len(self.selected_atoms)} atoms selected"
            + (f" (one of {count})" if count > 1 else "")
            + ". Ctrl+drag moves it, Ctrl+right-drag turns it; Shift+click empty space to deselect."
        )
        self._redraw()

    def clear_selection(self, widget=None):
        self.select_molecule(None)

    def valid_selection(self) -> List[int]:
        """Selected ids that still exist, in sketch order."""
        return [atom_id for atom_id in self.sketch.ordered_atom_ids() if atom_id in self.selected_atoms]

    def _transform_refusal(self) -> Optional[str]:
        """Why the selection cannot be moved rigidly right now, or ``None``."""
        selected = set(self.valid_selection())
        if not selected:
            return "Shift+click a molecule to select it before moving it."
        for bond in self.sketch.bonds:
            if (bond.atom_i in selected) != (bond.atom_j in selected):
                return ("The selection is now bonded to atoms outside it, so moving it would "
                        "stretch that bond. Shift+click the molecule again to reselect it.")
        for constraint in self.sketch.constraints:
            inside = [atom_id in selected for atom_id in constraint.atom_ids]
            if any(inside) and not all(inside):
                return (f"The frozen {constraint.kind} "
                        f"{self.sketch.describe_parameter(constraint.atom_ids)} links the selection to "
                        "another molecule; remove it to move the molecule.")
        return None

    def _begin_transform(self) -> bool:
        """First movement of a Ctrl gesture: check it is allowed and record one undo step."""
        if self._transform_started:
            return True
        refusal = self._transform_refusal()
        if refusal:
            self.set_status(refusal)
            return False
        self._push_undo()
        self._transform_started = True
        return True

    def _translate_selection(self, dx: float, dy: float):
        if not self._begin_transform():
            return
        delta = self.camera.screen_to_sketch_delta(dx, dy)
        self._move_selection(lambda points: points + delta)

    def _rotate_selection(self, dx: float, dy: float):
        if not self._begin_transform():
            return
        rotation = self.camera.object_rotation(dx, dy)
        ids = self.valid_selection()
        centre = self.sketch.positions()[self._indices(ids)].mean(axis=0)
        self._move_selection(lambda points: (points - centre) @ rotation.T + centre)

    def _move_selection(self, transform):
        """Apply a rigid ``transform`` to the selected atoms only."""
        positions = self.sketch.positions()
        indices = self._indices(self.valid_selection())
        positions[indices] = transform(positions[indices])
        self.sketch.update_positions(positions)
        self._redraw()

    def _finish_transform(self, verb: str):
        if self._transform_started:
            self.set_status(f"{verb}; its shape is unchanged. Undo puts it back.")
            self._redraw()  # the gesture is over: views catch up what they skipped
        self._transform_started = False

    def _indices(self, atom_ids) -> List[int]:
        index_of = {atom_id: index for index, atom_id in enumerate(self.sketch.ordered_atom_ids())}
        return [index_of[atom_id] for atom_id in atom_ids]

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
            if self.replace_mode:
                self.replace_atom(atom_id)
            else:
                self.pick_atom(atom_id)
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

    def _push_undo(self, snapshot=None, keep_result: bool = False):
        """Snapshot before a change, and drop any geometry that no longer matches.

        ``snapshot`` is a state taken earlier, for a change that is only
        recorded once it is known to have succeeded (an import).
        ``keep_result`` is for a change that moves no atom (freezing or
        releasing a parameter): the optimized geometry still matches the
        drawing, but a run still in flight is stale and will be discarded.
        """
        self.undo_stack.append(self.sketch.snapshot() if snapshot is None else snapshot)
        if len(self.undo_stack) > self.UNDO_DEPTH:
            self.undo_stack.pop(0)
        if not keep_result:
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
        if result.constraint_count:
            lines.append(
                f"Frozen parameters held: {result.constraint_count} "
                f"(largest deviation {result.constraint_deviation:.1e})"
            )

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
            frame = self.last_result.to_frame()
            self.viewer.load_designed_molecule(frame, bonds=self.sketch.bond_index_pairs())
        except Exception as exc:
            return False, f"The 3D viewer could not accept the structure: {exc}"
        bonds = self.current_structure().bonds
        self._sent_structure = StructureData(
            [element for element, _ in frame], np.array([xyz for _, xyz in frame], dtype=float), bonds=bonds
        )

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
    # Measuring and editing geometric parameters
    # ------------------------------------------------------------------
    def valid_picks(self) -> List[int]:
        """Picked atom ids that still exist, in pick order."""
        existing = set(self.sketch.ordered_atom_ids())
        self.picked_atoms = [atom_id for atom_id in self.picked_atoms if atom_id in existing]
        return list(self.picked_atoms)

    def pick_atom(self, atom_id: int):
        """Add an atom to the measurement, or take it out if already picked.

        Picks are ordered -- the second of three is the angle's vertex -- and a
        fifth pick starts a new measurement with that atom.
        """
        picks = self.valid_picks()
        if atom_id in picks:
            picks.remove(atom_id)
        elif len(picks) >= MAX_PICKS:
            picks = [atom_id]
        else:
            picks.append(atom_id)
        self.picked_atoms = picks
        self.set_status(self.measurement_status())
        self._redraw()

    def clear_picks(self, widget=None):
        self.picked_atoms = []
        self.set_status("Measurement cleared.")
        self._redraw()

    def _angstrom_positions(self) -> np.ndarray:
        """Every atom in Angstrom, in sketch order -- a rotation of the drawing, so
        lengths, angles and dihedral signs are those of the molecule itself."""
        return sketch_to_angstrom(self.sketch.positions())

    def current_measurement(self) -> Optional[Tuple[str, float]]:
        """``(kind, value)`` of the picked atoms, or ``None`` for fewer than two."""
        picks = self.valid_picks()
        if len(picks) < 2:
            return None
        kind = kind_for(picks)
        return kind, measure(kind, self._angstrom_positions(), self._indices(picks))

    def measurement_text(self) -> str:
        """'Bond length 1-2 (C-C): 1.5230 Å' for the picks; '' for fewer than two."""
        measurement = self.current_measurement()
        if measurement is None:
            return ""
        kind, value = measurement
        text = f"{KIND_LABELS[kind]} {self.sketch.describe_parameter(self.picked_atoms)}: "
        if kind == DIHEDRAL and is_degenerate(DIHEDRAL, self._angstrom_positions(),
                                              self._indices(self.picked_atoms)):
            text += "undefined (three atoms in a line)"
        else:
            text += format_value(kind, value)
        frozen = self.sketch.find_constraint(self.picked_atoms)
        if frozen is not None:
            text += f", frozen at {format_value(kind, frozen.value)}"
        return text

    def measurement_status(self) -> str:
        picks = self.valid_picks()
        if not picks:
            return ("Click atoms in order to measure: 2 for a bond length, 3 for a bond angle, "
                    "4 for a dihedral angle.")
        if len(picks) == 1:
            return (f"Picked {self.sketch.describe_atom(picks[0])} as #1. Click 1 to 3 more atoms to "
                    "measure, or drag from an atom to bond.")
        return self.measurement_text() + "."

    def _picked_parameter(self) -> Tuple[str, List[int]]:
        picks = self.valid_picks()
        if len(picks) < 2:
            raise ValueError("Click 2 atoms (bond length), 3 (bond angle) or 4 (dihedral angle) "
                             "in order first.")
        return kind_for(picks), picks

    def parameter_edit_info(self) -> Dict[str, object]:
        """What the edit dialog shows. Raises ``ValueError`` (with the reason) when
        the picked parameter cannot be edited at all, before anything is typed."""
        kind, picks = self._picked_parameter()
        coordinates = self._angstrom_positions()
        indices = self._indices(picks)
        moving_atoms(kind, len(coordinates), self.sketch.bond_index_pairs(), indices)
        if kind == DIHEDRAL:
            reason = is_degenerate(DIHEDRAL, coordinates, indices)
            if reason:
                raise GeometryEditError(f"The dihedral cannot be set: {reason}.")
        value = measure(kind, coordinates, indices)
        return {
            "kind": kind,
            "label": KIND_LABELS[kind],
            "atoms": self.sketch.describe_parameter(picks),
            "value": value,
            "value_text": format_value(kind, value),
            "unit": UNITS[kind],
            "frozen": self.sketch.find_constraint(picks) is not None,
        }

    def apply_parameter_edit(self, text: str, freeze: bool = False) -> Tuple[bool, str]:
        """Set the picked parameter to the typed value. No dialogs.

        One side of the molecule moves rigidly (molecularGeometry.set_parameter).
        Refused, with the reason, for a bad number, a ring bond, unbonded atoms,
        an undefined dihedral, or an edit that would change another frozen
        parameter. A parameter that is frozen -- or ``freeze`` -- is frozen at
        the new value. One undo step.
        """
        try:
            info = self.parameter_edit_info()
            kind = info["kind"]
            value = parse_parameter_value(kind, text)
        except ValueError as exc:
            return False, str(exc)
        picks = self.valid_picks()
        indices = self._indices(picks)
        coordinates = self._angstrom_positions()
        if freeze or info["frozen"]:
            problem = self._constraint_problem(kind, indices, value, coordinates, check_shape=False)
            if problem:
                return False, f"The new value cannot be frozen: {problem}."
        try:
            edited = set_parameter(kind, coordinates, self.sketch.bond_index_pairs(), indices, value)
        except GeometryEditError as exc:
            return False, str(exc)
        broken = self._changed_constraints(coordinates, edited, exclude=picks)
        if broken:
            return False, ("This edit would change frozen parameters: " + "; ".join(broken)
                           + ". Remove them first, or edit another parameter.")

        moved = np.any(edited != coordinates, axis=1)
        positions = self.sketch.positions()
        positions[moved] = angstrom_to_sketch(edited[moved])
        self._push_undo()
        self.sketch.update_positions(positions)
        frozen_note = ""
        if freeze or info["frozen"]:
            self.sketch.set_constraint(picks, value)
            frozen_note = " and frozen"
        message = (f"{info['label']} {info['atoms']} set{frozen_note} to {format_value(kind, value)} "
                   f"(was {info['value_text']}); {int(moved.sum())} atoms moved. Undo restores it.")
        self.set_status(message)
        self._redraw()
        return True, message

    def _changed_constraints(self, before, after, exclude) -> List[str]:
        """Frozen parameters (other than ``exclude``) whose value differs between two geometries."""
        changed = []
        for constraint in self.sketch.constraints:
            if constraint.same_parameter(exclude):
                continue
            indices = self._indices(constraint.atom_ids)
            old = measure(constraint.kind, before, indices)
            new = measure(constraint.kind, after, indices)
            difference = abs(new - old) if constraint.kind != DIHEDRAL else abs(wrap_degrees(new - old))
            if difference > (1e-6 if constraint.kind == DISTANCE else 1e-4):
                changed.append(f"{constraint.kind} {self.sketch.describe_parameter(constraint.atom_ids)}")
        return changed

    # ------------------------------------------------------------------
    # Frozen parameters
    # ------------------------------------------------------------------
    def _constraint_problem(self, kind, indices, value, coordinates, check_shape=True) -> Optional[str]:
        """Why the pre-optimizer could not hold this parameter, or ``None``."""
        if kind == ANGLE and not LINEAR_TOLERANCE_DEGREES < value < 180.0 - LINEAR_TOLERANCE_DEGREES:
            return ("a (nearly) linear angle cannot be held by the pre-optimizer - the constraint has "
                    "no defined direction there; the angle term already keeps a linear centre straight")
        if kind == DIHEDRAL and check_shape:
            return is_degenerate(DIHEDRAL, coordinates, indices)
        return None

    def freeze_picked(self, widget=None) -> Tuple[bool, str]:
        """Freeze the picked parameter at its current value."""
        try:
            kind, picks = self._picked_parameter()
        except ValueError as exc:
            self.set_status(str(exc))
            return False, str(exc)
        coordinates = self._angstrom_positions()
        indices = self._indices(picks)
        value = measure(kind, coordinates, indices)
        label = f"{KIND_LABELS[kind].lower()} {self.sketch.describe_parameter(picks)}"
        problem = self._constraint_problem(kind, indices, value, coordinates)
        if problem:
            message = f"The {label} cannot be frozen: {problem}."
            self.set_status(message)
            return False, message
        self._push_undo(keep_result=True)
        self.sketch.set_constraint(picks, value)
        message = (f"Frozen: {label} at {format_value(kind, value)}. Pre-optimization will hold it; "
                   "remove it from the Frozen parameters list.")
        self.set_status(message)
        self._redraw()
        return True, message

    async def freeze_button(self, widget=None):
        frozen, message = self.freeze_picked()
        if not frozen:
            await self._show_error("Cannot freeze this parameter", message)

    # Short names for the Frozen parameters table, whose columns are narrow.
    CONSTRAINT_KIND_NAMES = {DISTANCE: "Length", ANGLE: "Angle", DIHEDRAL: "Dihedral"}

    def constraint_rows(self) -> List[Tuple[str, str, str]]:
        """``(kind, atom numbers, value)`` per frozen parameter, in the order frozen."""
        return [
            (self.CONSTRAINT_KIND_NAMES[c.kind],
             "-".join(str(self.sketch.atom_number(atom_id)) for atom_id in c.atom_ids),
             format_value(c.kind, c.value))
            for c in self.sketch.constraints
        ]

    def remove_constraint_at(self, index: int) -> bool:
        constraints = self.sketch.constraints
        if not 0 <= index < len(constraints):
            self.set_status("Select a frozen parameter in the list to remove it.")
            return False
        constraint = constraints[index]
        self._push_undo(keep_result=True)
        self.sketch.remove_constraint(constraint)
        self.set_status(f"No longer frozen: {constraint.kind} "
                        f"{self.sketch.describe_parameter(constraint.atom_ids)}.")
        self._redraw()
        return True

    def remove_all_constraints(self, widget=None):
        if not self.sketch.constraints:
            self.set_status("Nothing is frozen.")
            return
        count = len(self.sketch.constraints)
        self._push_undo(keep_result=True)
        self.sketch.clear_constraints()
        self.set_status(f"{count} frozen parameters removed. Undo brings them back.")
        self._redraw()

    # ------------------------------------------------------------------
    # Hydrogens
    # ------------------------------------------------------------------
    def add_hydrogens(self) -> Tuple[bool, str, List[str]]:
        """Complete valences with hydrogens -- the selected molecule only, if one is
        selected. Returns ``(added, message, issues)``; flagged atoms are left alone."""
        if self.sketch.is_empty:
            return False, "Draw or import a structure first.", []
        ids = self.sketch.ordered_atom_ids()
        structure = self.current_structure()
        selected = self.valid_selection()
        scope = set(self._indices(selected)) if selected else None
        plan = plan_hydrogens(structure.elements, self._angstrom_positions(), structure.bonds, scope)
        where = " to the selected molecule" if selected else ""
        if not plan.count:
            message = (f"No hydrogens were added{where}: every atom is complete"
                       + (" or flagged below." if plan.issues else "."))
            self.set_status(message)
            return False, message, plan.issues
        self._push_undo()
        for parent, position in plan.additions:
            hydrogen = self.sketch.add_atom("H", *angstrom_to_sketch(position)[0])
            self.sketch.add_bond(ids[parent], hydrogen.atom_id)
            if selected:
                self.selected_atoms.add(hydrogen.atom_id)
        message = f"Added {plan.count} hydrogens{where}. Undo removes them."
        if plan.issues:
            message += f" {len(plan.issues)} atoms were left unchanged (see the details)."
        self.set_status(message)
        self._redraw()
        return True, message, plan.issues

    async def add_hydrogens_button(self, widget=None):
        added, message, issues = self.add_hydrogens()
        if issues:
            details = "\n".join(f"  - {issue}" for issue in issues)
            await self._show_info(
                "Add hydrogens",
                f"{message}\n\nNot completed - these need your decision:\n{details}\n\n"
                "Add hydrogens to them by hand, or change their bonds, if needed.",
            )
        elif not added:
            await self._show_info("Add hydrogens", message)

    # ------------------------------------------------------------------
    # Canvas background
    # ------------------------------------------------------------------
    def set_canvas_background(self, color: str) -> str:
        """A preset name or ``#rrggbb``; ``ValueError`` otherwise. A view setting, not an edit."""
        value = CANVAS_BACKGROUNDS.get(color) or parse_hex_color(color)
        self.canvas_background = value
        self._broadcast("on_background_changed", value)
        self.set_status(f"Canvas background set to {value}.")
        self._redraw()
        return value

    # ------------------------------------------------------------------
    # Bringing a molecule back from the 3D viewer
    # ------------------------------------------------------------------
    def structure_from_viewer(self) -> StructureData:
        """The molecule the 3D viewer shows now, in Angstrom, with its bonds.

        The viewer keeps connectivity but no bond orders, so its bonds come in
        as single -- except for a structure this canvas sent there itself,
        whose drawn orders are restored. Raises ``ValueError`` with the reason.
        """
        if self.viewer is None:
            raise ValueError("No 3D viewer is attached to this design session.")
        frame, pairs = self.viewer.current_frame_with_bonds()
        if not frame:
            raise ValueError("The 3D viewer shows no molecule.")
        elements = [str(element) for element, _ in frame]
        coordinates = np.array([position for _, position in frame], dtype=float).reshape(-1, 3)
        bonds = sorted((min(i, j), max(i, j), 1.0) for i, j in pairs)
        sent = self._sent_structure
        if (sent is not None and sent.elements == elements
                and np.allclose(sent.coordinates, coordinates, atol=1e-6)):
            bonds = list(sent.bonds)
            notes = ["This is the structure sent from the canvas: its bond orders were restored."]
        elif len(elements) > 1:
            notes = [f"The 3D viewer stores no bond orders: its {len(bonds)} bonds came in as single "
                     "bonds. Click a bond to set its type before pre-optimizing."]
        else:
            notes = []
        return StructureData(elements, coordinates, bonds=bonds, notes=notes)

    def bring_from_viewer(self, mode: str = "replace") -> Tuple[bool, str, List[str]]:
        """Load the viewer's molecule: ``mode`` "replace" or "add" (beside the drawing). No dialogs."""
        try:
            structure = self.structure_from_viewer()
        except ValueError as exc:
            return False, str(exc), []
        return self._load_structure(structure, "the 3D viewer", mode)

    async def import_from_viewer(self, widget=None):
        """Bring the 3D viewer's molecule to the canvas without overwriting a drawing.

        On a non-empty canvas the user is asked first: Yes adds the molecule
        beside the drawing, No changes nothing (Clear first to replace it).
        """
        mode = "replace"
        if not self.sketch.is_empty:
            window = self._dialog_window()
            add = window is not None and await window.dialog(toga.QuestionDialog(
                "The canvas already has a drawing",
                f"The design canvas holds {self.sketch.atom_count} atoms. Add the 3D viewer's molecule "
                "beside them?\n\nYes: add it beside your drawing - nothing is lost.\n"
                "No: cancel. To replace the drawing instead, press Clear on the canvas first "
                "(Undo brings it back).",
            ))
            if not add:
                self.set_status("Nothing was brought from the 3D viewer; the drawing is unchanged.")
                return False
            mode = "add"
        loaded, message, notes = self.bring_from_viewer(mode)
        if not loaded:
            self.set_status(message)
            await self._show_error("Cannot bring the molecule from the 3D viewer", message)
            return False
        self.open_design_canvas()
        self.set_status(message)
        await self._show_info("Molecule brought from the 3D viewer", "\n\n".join([message] + notes))
        return True

    # ------------------------------------------------------------------
    # Importing a structure file
    # ------------------------------------------------------------------
    def import_structure(self, path) -> Tuple[bool, str, List[str]]:
        """Replace the drawing with the structure in ``path``. No dialogs.

        Returns ``(ok, message, notes)``; on failure the drawing is untouched
        and ``message`` says what is wrong with the file. The import is one
        undo step, so Undo brings the previous drawing back. The molecule is
        centred on the canvas, and the view is reset and zoomed out if that is
        needed to show all of it. Split from :meth:`import_file` so it is
        testable without an event loop, like :meth:`transfer_to_viewer`.
        """
        try:
            structure = read_structure_file(path)
        except ValueError as exc:
            return False, str(exc), []
        name = Path(path).name
        label = format_for_path(path).label
        return self._load_structure(structure, f"{name} ({label})", "replace", name=name, verb="Imported")

    def _load_structure(self, structure: StructureData, source: str, mode: str,
                        name: Optional[str] = None, verb: str = "Loaded") -> Tuple[bool, str, List[str]]:
        """Put a structure on the canvas: replacing the drawing, or beside it."""
        name = name or source
        adding = mode == "add" and not self.sketch.is_empty
        total = structure.atom_count + (self.sketch.atom_count if adding else 0)
        if total > MAX_IMPORT_ATOMS:
            return False, (
                f"'{name}' holds {structure.atom_count} atoms; the design canvas holds at most "
                f"{MAX_IMPORT_ATOMS}. Open large systems in the 3D viewer with "
                "'Display Molecule/Trajectory' instead."
            ), []

        coordinates = structure.coordinates - structure.coordinates.mean(axis=0)
        positions = angstrom_to_sketch(coordinates) + np.array([*self.view_centre, 0.0])
        if adding:
            # Beside the drawing: 2 A to the right of it, level with its middle.
            existing = self.sketch.positions()
            positions += np.array([
                existing[:, 0].max() - positions[:, 0].min() + 2.0 * PIXELS_PER_ANGSTROM,
                existing[:, 1].mean() - positions[:, 1].mean(),
                existing[:, 2].mean() - positions[:, 2].mean(),
            ])
        had_drawing = not self.sketch.is_empty
        before = self.sketch.snapshot()
        try:
            if adding:
                self.sketch.add_structure(structure.elements, positions, structure.bonds)
            else:
                self.sketch.replace_contents(structure.elements, positions, structure.bonds)
        except ValueError as exc:  # both change nothing when they refuse
            return False, f"'{name}' could not be placed on the canvas: {exc}", []
        self._push_undo(before)
        self._fit_view()

        message = (f"{'Added' if adding else verb} {structure.atom_count} atoms and {len(structure.bonds)} bonds "
                   f"from {source}.")
        if had_drawing:
            message += " Undo restores the previous drawing."
        self.set_status(message)
        self._redraw()
        return True, message, list(structure.notes)

    def _fit_view(self):
        """Reset the view onto the whole drawing, zooming out -- never in -- to fit it."""
        self.camera.reset()
        if self.sketch.is_empty:
            return
        centre = np.array(self.view_centre)
        points = self.sketch.positions()[:, :2]
        low, high = points.min(axis=0), points.max(axis=0)
        self.camera.pan_by(*(centre - (low + high) / 2.0))
        extent = (high - low) / 2.0
        room = centre - 2.0 * self.ATOM_RADIUS
        scales = [room[axis] / extent[axis] for axis in (0, 1) if extent[axis] > 0.0]
        if scales and min(scales) < 1.0:
            self.camera.zoom_by(min(scales), *centre)

    async def import_file(self, widget=None):
        """Ask for a structure file and load it onto the canvas."""
        window = self._dialog_window()
        if window is None:
            return False
        extensions = [ext.lstrip(".") for fmt in READABLE_FORMATS for ext in fmt.extensions]
        selected = await window.dialog(
            toga.OpenFileDialog(title="Import structure", file_types=extensions)
        )
        if isinstance(selected, (list, tuple)):
            selected = selected[0] if selected else None
        if not selected:
            return False

        imported, message, notes = self.import_structure(selected)
        if not imported:
            self.set_status("The file could not be imported; the drawing is unchanged.")
            await self._show_error("Cannot import this file", message)
            return False
        details = "\n\n".join([message] + notes)
        await self._show_info("Structure imported", details)
        return True

    # ------------------------------------------------------------------
    # Saving the displayed structure
    # ------------------------------------------------------------------
    def current_structure(self) -> StructureData:
        """The structure on the canvas in Angstrom, with its bonds and their orders.

        Built from the sketch itself, so it is exactly what the canvas shows --
        the optimized geometry after a pre-optimization, the drawing otherwise.
        The view (rotation, pan, zoom) is not an edit and is not applied; the
        coordinates are centred on the origin, like those sent to the 3D viewer.
        """
        elements = [atom.element for atom in self.sketch.atoms]
        coordinates = sketch_to_angstrom(self.sketch.positions())
        if elements:
            coordinates = coordinates - coordinates.mean(axis=0)
        index_of = {atom_id: index for index, atom_id in enumerate(self.sketch.ordered_atom_ids())}
        bonds = sorted(
            (min(index_of[b.atom_i], index_of[b.atom_j]), max(index_of[b.atom_i], index_of[b.atom_j]),
             float(b.order))
            for b in self.sketch.bonds
        )
        state = "pre-optimized" if self.last_result is not None else "as drawn, not optimized"
        title = f"Generated by the gQTEA molecular design tool ({state}); coordinates in Angstrom"
        return StructureData(elements, coordinates, bonds=bonds, title=title)

    def frame_xyz_text(self) -> str:
        """The structure on the canvas as XYZ text: element and x, y, z in Angstrom."""
        return write_xyz_text(self.current_structure())

    def frame_text(self, format_key: Optional[str] = None, name: str = "designed_molecule") -> str:
        """The structure on the canvas as a ``format_key`` file (default: the selected format)."""
        format_key = format_key or self.export_format
        if format_key == "XYZ":
            return self.frame_xyz_text()
        return structure_to_text(self.current_structure(), format_key, name)

    def set_export_format(self, format_key: str):
        """Choose the format Save Frame writes; XYZ until changed."""
        fmt = FORMATS_BY_KEY.get(format_key)
        if fmt is None or not fmt.writable:
            raise ValueError(f"Structures cannot be saved as {format_key}.")
        self.export_format = format_key
        self.set_status(f"Save Frame will write {fmt.label} files ({', '.join(fmt.extensions)}).")

    def _format_for_saving(self, output_path) -> str:
        """The selected format, unless the chosen name ends in another supported one.

        A typed ``frame.pdb`` is then written as PDB, so a file's content always
        matches its extension; a name with any other extension gets the
        selected format, as every name got XYZ before formats could be chosen.
        """
        try:
            return format_for_path(output_path, "write").key
        except ValueError:
            return self.export_format

    def write_frame(self, output_path: str, format_key: Optional[str] = None) -> Tuple[bool, str]:
        """Write the displayed structure to ``output_path``. No dialogs.

        Split from :meth:`save_frame` so the export stays testable without an
        event loop, the same seam as :meth:`transfer_to_viewer`.
        """
        if self.sketch.is_empty:
            return False, "Draw a structure on the canvas before saving it."
        format_key = format_key or self.export_format
        try:
            text = self.frame_text(format_key, name=Path(output_path).stem)
        except ValueError as exc:
            return False, f"The structure could not be saved as {format_key}: {exc}"
        try:
            with open(output_path, "w", encoding="utf-8") as handle:
                handle.write(text)
        except OSError as exc:
            return False, f"The structure could not be saved: {exc}"
        message = f"{self.sketch.atom_count} atoms saved to {output_path}"
        if format_key != "XYZ":
            message += f" ({FORMATS_BY_KEY[format_key].label})"
        return True, message

    async def save_frame(self, widget=None):
        """Ask for a file name and save the displayed structure in the selected format."""
        if self.sketch.is_empty:
            message = "Draw a structure on the canvas before saving it."
            self.set_status(message)
            await self._show_error("Nothing to save", message)
            return False

        window = self._dialog_window()
        if window is None:
            return False
        fmt = FORMATS_BY_KEY[self.export_format]
        selected = await window.dialog(
            toga.SaveFileDialog(
                title=f"Save Frame as {fmt.label}",
                suggested_filename=f"designed_molecule{fmt.extensions[0]}",
                file_types=[ext.lstrip(".") for ext in fmt.extensions],
            )
        )
        if not selected:
            return False

        format_key = self._format_for_saving(selected)
        saved, message = self.write_frame(str(selected), format_key)
        self.set_status(message)
        if not saved:
            await self._show_error("Save Error", message)
            return False
        notes = export_notes(self.current_structure(), format_key)
        await self._show_info("Frame Saved", "\n\n".join([f"Structure saved to:\n{selected}"] + notes))
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
            _set_label_text(self.controller.launcher_status_label, message)


class DesignCanvasWindow:
    """The drawing window: 3D canvas, edit and view actions, and a status line."""

    TITLE = "Molecular Design Canvas"
    CANVAS_WIDTH = CANVAS_WIDTH
    CANVAS_HEIGHT = CANVAS_HEIGHT
    # Canvas plus the side panel; Toga already widened the old 600 px window
    # to ~850 for the hint line, so the panel costs no extra width.
    WINDOW_SIZE = (860, 680)
    PANEL_WIDTH = 250
    HINT_COLOR = "#666666"

    def __init__(self, controller: MolecularDesignController):
        self.controller = controller
        self.canvas = None
        self._panel_signature = None
        # Paint frames as one native bitmap until that fails once; then Toga draws.
        self._bitmap_painting = True
        self._bitmap_active = False
        self._frame_bitmap = None
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
        top = toga.Box(style=Pack(direction=ROW, align_items=CENTER))
        top.add(self.element_label)
        top.add(toga.Box(style=Pack(flex=1)))
        self.replace_switch = toga.Switch(
            self.replace_switch_text(self.controller.selected_element),
            value=self.controller.replace_mode,
            style=Pack(margin=(0, 0, 6, 0)),
        )
        # Wired after the value is set, so building the window is not reported as a choice.
        self.replace_switch.on_change = lambda widget, **kwargs: self.controller.set_replace_mode(widget.value)
        top.add(self.replace_switch)
        content.add(top)
        content.add(
            toga.Label(
                "Click empty space to place an atom - drag atom to atom to bond them, or "
                "atom to empty space to grow a bonded atom - click a bond to change its type.",
                style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 2, 0)),
            )
        )
        content.add(
            toga.Label(
                "Right-drag to rotate - drag empty space to pan - right-click to delete - "
                "click atoms in order to measure.",
                style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 2, 0)),
            )
        )
        content.add(
            toga.Label(
                "Shift+click a molecule to select it - Ctrl+drag moves it - "
                "Ctrl+right-drag turns it - Shift+click empty space to deselect.",
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
                background_color=self.controller.canvas_background,
            ),
        )
        middle = toga.Box(style=Pack(direction=ROW))
        middle.add(self.canvas)
        middle.add(self._build_side_panel())
        content.add(middle)
        content.add(self._build_actions())
        content.add(self._build_view_actions())
        content.add(self._build_file_actions())
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
        """Zoom and reset -- the view controls the mouse cannot reach on its own."""
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
        view.add(toga.Button("Reset view", on_press=self.controller.reset_view,
                             style=Pack(width=100, margin=(0, 16, 0, 0))))
        self.numbers_switch = toga.Switch("Atom numbers", value=self.controller.show_atom_numbers)
        self.numbers_switch.on_change = lambda widget, **kwargs: self.controller.set_show_atom_numbers(widget.value)
        view.add(self.numbers_switch)
        return view

    # Save-format Selection labels, each mapped to its molecularFileFormats key.
    EXPORT_FORMAT_LABELS = {f"{fmt.label} ({fmt.extensions[0]})": fmt.key for fmt in WRITABLE_FORMATS}

    def _build_file_actions(self) -> toga.Box:
        """Import a structure file; Save Frame in the format chosen beside it."""
        row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(6, 0, 0, 0)))
        row.add(toga.Label("File:", style=Pack(width=44, color=self.HINT_COLOR)))
        row.add(
            toga.Button(
                "Import structure...",
                on_press=self.controller.import_file,
                style=Pack(width=140, margin=(0, 6, 0, 0)),
            )
        )
        row.add(
            toga.Button(
                "From 3D viewer",
                on_press=self.controller.import_from_viewer,
                style=Pack(width=120, margin=(0, 6, 0, 0)),
            )
        )
        row.add(toga.Box(style=Pack(flex=1)))
        labels = list(self.EXPORT_FORMAT_LABELS)
        self.export_selection = toga.Selection(items=labels, style=Pack(width=180, margin=(0, 6, 0, 0)))
        current = [label for label, key in self.EXPORT_FORMAT_LABELS.items()
                   if key == self.controller.export_format]
        self.export_selection.value = current[0] if current else labels[0]
        # Wired after the value is set, so building the window is not reported as a choice.
        self.export_selection.on_change = self._on_export_format_change
        row.add(self.export_selection)
        row.add(toga.Button("Save Frame", on_press=self.controller.save_frame, style=Pack(width=110)))
        return row

    def _on_export_format_change(self, widget, **kwargs):
        key = self.EXPORT_FORMAT_LABELS.get(str(widget.value))
        if key is not None:
            self.controller.set_export_format(key)

    # ------------------------------------------------------------------
    # Side panel: measurement, frozen parameters, molecule, canvas
    # ------------------------------------------------------------------
    def _heading(self, text: str) -> toga.Label:
        return toga.Label(text, style=Pack(font_size=10, font_weight="bold", margin=(6, 0, 3, 0)))

    def _small_button(self, text, handler, width=78, last=False) -> toga.Button:
        return toga.Button(text, on_press=handler, style=Pack(width=width, margin=(0, 0 if last else 4, 0, 0)))

    def _build_side_panel(self) -> toga.Box:
        panel = toga.Box(style=Pack(direction=COLUMN, width=self.PANEL_WIDTH, margin=(0, 0, 0, 10)))
        c = self.controller

        panel.add(self._heading("Measure"))
        self.measure_label = toga.Label("", style=Pack(font_size=9, margin=(0, 0, 4, 0)))
        panel.add(self.measure_label)
        row = toga.Box(style=Pack(direction=ROW))
        row.add(self._small_button("Edit value...", c.open_parameter_editor, width=86))
        row.add(self._small_button("Freeze", c.freeze_button, width=70))
        row.add(self._small_button("Clear", c.clear_picks, width=70, last=True))
        panel.add(row)

        panel.add(self._heading("Frozen parameters"))
        self.constraint_table = toga.Table(
            headings=["Kind", "Atoms", "Value"], data=[], style=Pack(height=100, width=self.PANEL_WIDTH)
        )
        panel.add(self.constraint_table)
        row = toga.Box(style=Pack(direction=ROW, margin=(4, 0, 0, 0)))
        row.add(self._small_button("Remove", self._remove_selected_constraint, width=86))
        row.add(self._small_button("Remove all", c.remove_all_constraints, width=90, last=True))
        panel.add(row)

        panel.add(self._heading("Molecule"))
        self.selection_label = toga.Label("", style=Pack(font_size=9, margin=(0, 0, 4, 0)))
        panel.add(self.selection_label)
        row = toga.Box(style=Pack(direction=ROW))
        row.add(self._small_button("Add hydrogens", c.add_hydrogens_button, width=110))
        row.add(self._small_button("Deselect", c.clear_selection, width=80, last=True))
        panel.add(row)

        panel.add(self._heading("Canvas"))
        row = toga.Box(style=Pack(direction=ROW, align_items=CENTER))
        row.add(toga.Label("Background", style=Pack(width=74, font_size=9)))
        labels = list(CANVAS_BACKGROUNDS)
        self.background_selection = toga.Selection(items=labels, style=Pack(width=170))
        current = [label for label, value in CANVAS_BACKGROUNDS.items() if value == c.canvas_background]
        self.background_selection.value = current[0] if current else labels[0]
        self.background_selection.on_change = self._on_background_preset
        row.add(self.background_selection)
        panel.add(row)
        row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(4, 0, 0, 0)))
        row.add(toga.Label("Custom", style=Pack(width=74, font_size=9)))
        self.background_input = toga.TextInput(
            placeholder="#rrggbb, Enter", on_confirm=self._on_background_custom, style=Pack(width=170)
        )
        row.add(self.background_input)
        panel.add(row)
        self._refresh_panel()
        return panel

    def _remove_selected_constraint(self, widget=None):
        row = self.constraint_table.selection
        if row is None:
            self.controller.set_status("Select a frozen parameter in the list to remove it.")
            return
        self.controller.remove_constraint_at(self.constraint_table.data.index(row))

    def _on_background_preset(self, widget, **kwargs):
        if widget.value:
            self.controller.set_canvas_background(str(widget.value))

    def _on_background_custom(self, widget, **kwargs):
        try:
            self.controller.set_canvas_background(widget.value)
        except ValueError as exc:
            self.controller.set_status(str(exc))

    def _refresh_panel(self):
        """Bring the side panel up to date, touching only what changed: every
        Toga label assignment re-lays-out the whole window."""
        c = self.controller
        picks = c.valid_picks()
        if len(picks) >= 2:
            # "Bond length 1-2 (C-C):" on one line, the value on the next.
            measure_text = c.measurement_text().replace(": ", ":\n", 1)
        elif picks:
            measure_text = f"#1 {c.sketch.describe_atom(picks[0])}\nClick 1-3 more atoms."
        else:
            measure_text = "Click 2, 3 or 4 atoms in order:\nlength, angle or dihedral."
        _set_label_text(self.measure_label, measure_text)
        selected = c.valid_selection()
        _set_label_text(self.selection_label, (
            f"Selected: {len(selected)} atoms (Ctrl+drag / Ctrl+right-drag).\nAdd hydrogens: this molecule."
            if selected else "None selected (Shift+click one).\nAdd hydrogens: every molecule."
        ))
        rows = c.constraint_rows()
        if rows != self._panel_signature:
            self._panel_signature = rows
            self.constraint_table.data = rows

    def on_background_changed(self, color: str):
        if self.canvas is not None:
            self.canvas.style.background_color = color

    # ------------------------------------------------------------------
    # Controller hooks
    # ------------------------------------------------------------------
    def on_element_changed(self, symbol: str):
        _set_label_text(self.element_label, f"Selected element: {symbol}")
        switch = getattr(self, "replace_switch", None)
        if switch is not None and switch.text != self.replace_switch_text(symbol):
            switch.text = self.replace_switch_text(symbol)

    @staticmethod
    def replace_switch_text(symbol: str) -> str:
        return f"Click an atom to replace it with {symbol}"

    def on_replace_mode_changed(self, enabled: bool):
        switch = getattr(self, "replace_switch", None)
        if switch is not None and bool(switch.value) != enabled:
            switch.value = enabled

    def on_status_changed(self, message: str):
        _set_label_text(self.status_label, message)

    def on_sketch_changed(self):
        """Repaint at once; the side panel only when no mouse gesture is under
        way (nothing in it changes mid-drag, and Toga re-lays-out the window on
        every label update, which made dragging lag)."""
        if self.canvas is None:
            return
        self.on_paint_requested()
        if not self.controller.gesture_active:
            self._refresh_panel()

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
    # Text sizes in points, the unit Toga fonts use.
    LABEL_FONT_SIZE = 9
    TAG_FONT_SIZE = 9
    MEASURE_FONT_SIZE = 10

    def on_paint_requested(self):
        """Paint the current drawing now (also used by :meth:`on_sketch_changed`)."""
        if not self._paint_bitmap():
            self._draw_sketch()
            self.canvas.redraw()
        self._paint_now()

    def _paint_now(self):
        """Have WinForms paint the invalidated canvas immediately.

        A plain Invalidate posts a low-priority paint message, which Windows
        delivers only once no input is waiting -- so while the pointer keeps
        moving, the picture can trail the gesture. ``Update`` paints at once.
        """
        native = getattr(getattr(self.canvas, "_impl", None), "native", None)
        try:
            native.Update()
        except Exception:
            pass

    def _paint_bitmap(self) -> bool:
        """WinForms fast path: render the frame with Pillow into a native bitmap
        that :meth:`_blit_frame` copies onto the control -- one blit per paint.
        Returns False (and stops trying) wherever that is not possible; Toga's
        canvas then draws instead."""
        if not getattr(self, "_bitmap_painting", False):
            return False
        native = getattr(getattr(self.canvas, "_impl", None), "native", None)
        if native is None or not hasattr(native, "ClientSize"):
            self._bitmap_painting = False
            return False
        try:
            width, height = int(native.ClientSize.Width), int(native.ClientSize.Height)
            if width <= 0 or height <= 0:
                return False  # not laid out yet: Toga draws this one frame
            image = render_scene(self.scene_primitives(), (width, height), width / self.CANVAS_WIDTH)
            bitmap = pil_to_net_bitmap(image)
            if not self._bitmap_active:
                native.Paint += self._blit_frame
                self.canvas.context.clear()  # Toga's own shapes would paint over the frame
                self._bitmap_active = True
            previous, self._frame_bitmap = getattr(self, "_frame_bitmap", None), bitmap
            native.Invalidate()
            if previous is not None:
                previous.Dispose()  # only ever a bitmap this window made
            return True
        except Exception:
            self._stop_bitmap_painting(native)
            return False

    def _blit_frame(self, sender, event):
        """Paint handler: copy the current frame 1:1 (no scaling, no blending)."""
        bitmap = getattr(self, "_frame_bitmap", None)
        if bitmap is None:
            return
        try:
            from System.Drawing import GraphicsUnit, Rectangle
            from System.Drawing.Drawing2D import CompositingMode, InterpolationMode

            graphics = event.Graphics
            graphics.CompositingMode = CompositingMode.SourceCopy
            graphics.InterpolationMode = InterpolationMode.NearestNeighbor
            area = Rectangle(0, 0, bitmap.Width, bitmap.Height)
            graphics.DrawImage(bitmap, area, area, GraphicsUnit.Pixel)
        except Exception:
            self._stop_bitmap_painting(sender)

    def _stop_bitmap_painting(self, native):
        """Fall back to Toga's canvas drawing for good."""
        self._bitmap_painting = False
        if getattr(self, "_bitmap_active", False):
            try:
                native.Paint -= self._blit_frame
            except Exception:
                pass
        self._bitmap_active = False
        self._frame_bitmap = None

    def _draw_sketch(self):
        """Toga-canvas path: the same primitives, painted through the canvas context."""
        context = getattr(self.canvas, "context", None)
        if context is None:
            return
        context.clear()
        self._paint_with_toga(context, self.scene_primitives())

    def scene_primitives(self) -> List[tuple]:
        """The frame as molecularCanvasRender primitives, back to front.

        The molecule with a depth cue, the selected molecule's halo, the drag
        rubber band, then the measurement picks with their order tags and the
        measured value on top.
        """
        controller = self.controller
        background = controller.canvas_background
        palette = canvas_palette(background)
        primitives: List[tuple] = [("rect", 0, 0, self.CANVAS_WIDTH, self.CANVAS_HEIGHT, background)]

        screen = controller.screen_positions()
        order = controller.drawing_order(screen)
        radii = {
            item.atom_id: self._atom_radius(fade) for kind, item, fade in order if kind == "atom"
        }
        selected = set(controller.valid_selection())
        preview_drawn = False
        for kind, item, fade in order:
            if kind == "bond":
                in_selection = item.atom_i in selected and item.atom_j in selected
                self._bond_primitives(primitives, item, screen, radii, fade, palette, background, in_selection)
                continue
            # The rubber band sits over every bond behind the first atom and
            # under the atoms, as it did on the flat canvas.
            if not preview_drawn:
                self._preview_primitives(primitives)
                preview_drawn = True
            self._atom_primitives(primitives, item, screen[item.atom_id], fade, palette,
                                  background, item.atom_id in selected)
        if not preview_drawn:
            self._preview_primitives(primitives)
        obstacles = self._pick_primitives(primitives, screen, radii, palette)
        self._number_primitives(primitives, screen, radii, palette, obstacles)
        return primitives

    def _preview_primitives(self, primitives):
        if self.controller.drag_preview is None:
            return
        from_x, from_y, to_x, to_y = self.controller.drag_preview
        primitives.append(("line", from_x, from_y, to_x, to_y, RUBBER_BAND_COLOR, 1.5, (4.0, 3.0)))

    def _atom_radius(self, fade: float) -> float:
        return self.controller.ATOM_RADIUS * (1.0 - DEPTH_SHRINK * fade)

    def _atom_primitives(self, primitives, atom, position, fade, palette, background, selected):
        """A disc that shrinks and fades towards the background with depth; an
        atom of the selected molecule gets an orange halo."""
        sx, sy, _ = position
        radius = self._atom_radius(fade)
        tint = DEPTH_FADE * fade
        if selected:
            halo = blend_color(SELECTION_COLOR, background, 0.45 + 0.3 * fade)
            primitives.append(("disc", sx, sy, radius + 5.0, halo, None, 0.0))
        color = ATOM_CANVAS_COLORS.get(atom.element, DEFAULT_ATOM_COLOR)
        primitives.append(("disc", sx, sy, radius, blend_color(color, background, tint),
                           blend_color(palette["outline"], background, tint), 1.0))
        primitives.append(("text", sx - 2.0 * len(atom.element) - 1.0, sy, atom.element,
                           blend_color(ATOM_LABEL_COLOR, background, tint), self.LABEL_FONT_SIZE,
                           False, "middle"))

    def _bond_primitives(self, primitives, bond: SketchBond, screen, radii, fade, palette,
                         background, selected):
        """One, two or three parallel lines -- or, for resonance, one solid
        line with a dashed line beside it -- offset sideways on screen.

        The lines run from disc edge to disc edge rather than centre to centre:
        a bond shares an end with each of its atoms, so no painting order can
        stop a centre-to-centre line crossing the disc of the atom it is
        painted after. Bonds of the selected molecule are drawn in orange.
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
        normal_x, normal_y = -unit_y, unit_x
        offset = self.controller.BOND_LINE_OFFSET
        base = SELECTION_COLOR if selected else palette["bond"]
        color = blend_color(base, background, DEPTH_FADE * fade)

        if bond.is_resonance:
            lines = [(-offset, None), (offset, (5.0, 4.0))]
        else:
            offsets = {
                1: (0.0,),
                2: (-offset, offset),
                3: (-2.0 * offset, 0.0, 2.0 * offset),
            }.get(int(bond.order), (0.0,))
            lines = [(line_offset, None) for line_offset in offsets]
        for line_offset, dash in lines:
            primitives.append((
                "line",
                start[0] + normal_x * line_offset, start[1] + normal_y * line_offset,
                end[0] + normal_x * line_offset, end[1] + normal_y * line_offset,
                color, 1.8, dash,
            ))

    def _pick_primitives(self, primitives, screen, radii, palette) -> List[Tuple[float, float, float, float]]:
        """Picked atoms: a ring, an order tag (#1-#4), a dotted path in pick
        order, and the measured value in the top-left corner. Returns the
        boxes that text occupies, for atom numbers to keep clear of."""
        picks = self.controller.valid_picks()
        if not picks:
            return []
        boxes = []
        points = [tuple(screen[atom_id][:2]) for atom_id in picks]
        if len(points) > 1:
            primitives.append(("polyline", points, PICK_COLOR, 1.2, (2.0, 3.0)))
        for number, atom_id in enumerate(picks, start=1):
            sx, sy = screen[atom_id][:2]
            radius = radii.get(atom_id, self.controller.ATOM_RADIUS) + 4.0
            primitives.append(("ring", sx, sy, radius, PICK_COLOR, 2.0))
            tag = f"#{number}"
            primitives.append(("text", sx + radius, sy - radius, tag, PICK_COLOR,
                               self.TAG_FONT_SIZE, True, "bottom"))
            width, height = text_box_size(tag, self.TAG_FONT_SIZE)
            boxes.append((sx + radius, sy - radius - height, sx + radius + width, sy - radius))
        text = self.controller.measurement_text()
        if text:
            primitives.append(("text", 8.0, 8.0, text, palette["text"], self.MEASURE_FONT_SIZE, False, "top"))
            width, height = text_box_size(text, self.MEASURE_FONT_SIZE)
            boxes.append((8.0, 8.0, 8.0 + width, 8.0 + height))
        return boxes

    NUMBER_FONT_SIZE = 8

    def _number_primitives(self, primitives, screen, radii, palette, obstacles):
        """Each atom's number beside it, where it covers no atom, bond or other text.

        Placement (molecularCanvasRender.place_atom_labels) is redone for every
        settled frame. During a mouse gesture each label keeps its offset from
        its atom instead, so dragging a large numbered molecule stays smooth;
        the layout is redone when the gesture ends.
        """
        controller = self.controller
        if not controller.show_atom_numbers or not screen:
            return
        ids = controller.sketch.ordered_atom_ids()
        texts = [str(number) for number in range(1, len(ids) + 1)]
        centres = np.array([screen[atom_id][:2] for atom_id in ids], dtype=float)
        key = (tuple(ids), controller.camera.zoom)
        cached = getattr(self, "_label_offsets", None)
        if controller.gesture_active and cached is not None and cached[0] == key:
            positions = centres + cached[1]
        else:
            index_of = {atom_id: index for index, atom_id in enumerate(ids)}
            bonds = [(index_of[b.atom_i], index_of[b.atom_j]) for b in controller.sketch.bonds]
            picked = set(controller.valid_picks())
            # A picked atom is drawn with a ring 4 px outside its disc (2 px wide):
            # its number has to clear the ring, not just the disc.
            positions = place_atom_labels(
                centres,
                [radii.get(atom_id, controller.ATOM_RADIUS) + (5.0 if atom_id in picked else 0.0)
                 for atom_id in ids],
                bonds,
                [text_box_size(text, self.NUMBER_FONT_SIZE) for text in texts],
                obstacles,
            )
            self._label_offsets = (key, positions - centres)
        for text, (x, y) in zip(texts, positions):
            primitives.append(("label", float(x), float(y), text, palette["number"], self.NUMBER_FONT_SIZE, False))

    TOGA_BASELINES = {"top": Baseline.TOP, "middle": Baseline.MIDDLE, "bottom": Baseline.BOTTOM}

    def _paint_with_toga(self, context, primitives):
        """Paint primitives through a Toga canvas context (the portable path)."""
        for primitive in primitives:
            kind = primitive[0]
            if kind == "rect":
                _, x0, y0, x1, y1, fill = primitive
                with context.Fill(color=fill) as path:
                    path.rect(x0, y0, x1 - x0, y1 - y0)
            elif kind == "disc":
                _, cx, cy, radius, fill, outline, width = primitive
                with context.Fill(color=fill) as path:
                    path.ellipse(cx, cy, radius, radius)
                if outline:
                    with context.Stroke(color=outline, line_width=width) as path:
                        path.ellipse(cx, cy, radius, radius)
            elif kind == "ring":
                _, cx, cy, radius, color, width = primitive
                with context.Stroke(color=color, line_width=width) as path:
                    path.ellipse(cx, cy, radius, radius)
            elif kind in ("line", "polyline"):
                if kind == "line":
                    _, x1, y1, x2, y2, color, width, dash = primitive
                    points = [(x1, y1), (x2, y2)]
                else:
                    _, points, color, width, dash = primitive
                with context.Stroke(color=color, line_width=width,
                                    line_dash=list(dash) if dash else None) as path:
                    path.move_to(*points[0])
                    for point in points[1:]:
                        path.line_to(*point)
            elif kind == "text":
                _, x, y, text, color, size, bold, baseline = primitive
                with context.Fill(color=color) as path:
                    path.write_text(text, x, y, _toga_font(size, bool(bold)), self.TOGA_BASELINES[baseline])
            elif kind == "label":
                _, x, y, text, color, size, bold = primitive
                width, _ = text_box_size(text, size)
                with context.Fill(color=color) as path:
                    path.write_text(text, x - (width - 2.0) / 2.0, y, _toga_font(size, bool(bold)),
                                    Baseline.MIDDLE)


def _set_label_text(label, text: str):
    """Assign only a different text: Toga relayouts the window on every assignment."""
    if label.text != text:
        label.text = text


@lru_cache(maxsize=16)
def _toga_font(size: int, bold: bool) -> Font:
    return Font(family="sans-serif", size=size, weight="bold" if bold else "normal")


class ParameterEditWindow:
    """The 'Edit value' dialog: the picked atoms, the current value, a new value.

    A thin view: validation and the edit itself are
    :meth:`MolecularDesignController.apply_parameter_edit`. A refused value is
    explained inside the dialog, which stays open so it can be corrected.
    """

    WINDOW_SIZE = (440, 240)
    ERROR_COLOR = "#c0392b"
    HINTS = {
        DISTANCE: "The second atom's side of the bond moves along it.",
        ANGLE: "The third atom's side turns about the vertex (the second atom).",
        DIHEDRAL: "The fourth atom's side turns about the central bond.",
    }

    def __init__(self, controller: MolecularDesignController, info: Dict[str, object]):
        self.controller = controller
        self.info = info
        kind = info["kind"]
        self.window = toga.Window(title=f"Set {str(info['label']).lower()}", size=self.WINDOW_SIZE,
                                  on_close=self._on_close)
        content = toga.Box(style=Pack(direction=COLUMN, margin=12))
        content.add(toga.Label(f"{info['label']}  -  atoms {info['atoms']}",
                               style=Pack(font_weight="bold", margin=(0, 0, 6, 0))))
        content.add(toga.Label(f"Current value: {info['value_text']}", style=Pack(margin=(0, 0, 6, 0))))
        row = toga.Box(style=Pack(direction=ROW, align_items=CENTER))
        row.add(toga.Label(f"New value ({info['unit']}):", style=Pack(width=120)))
        decimals = 4 if kind == DISTANCE else 2
        self.value_input = toga.TextInput(value=f"{float(info['value']):.{decimals}f}",
                                          on_confirm=self.apply, style=Pack(width=120))
        row.add(self.value_input)
        content.add(row)
        self.freeze_switch = toga.Switch("Freeze at the new value during pre-optimization",
                                         value=bool(info["frozen"]), style=Pack(margin=(6, 0, 0, 0)))
        content.add(self.freeze_switch)
        content.add(toga.Label(self.HINTS[kind], style=Pack(font_size=9, color="#666666", margin=(6, 0, 0, 0))))
        self.error_label = toga.Label("", style=Pack(font_size=9, color=self.ERROR_COLOR, margin=(4, 0, 0, 0)))
        content.add(self.error_label)
        buttons = toga.Box(style=Pack(direction=ROW, margin=(8, 0, 0, 0)))
        buttons.add(toga.Box(style=Pack(flex=1)))
        buttons.add(toga.Button("Apply", on_press=self.apply, style=Pack(width=80, margin=(0, 6, 0, 0))))
        buttons.add(toga.Button("Cancel", on_press=self._cancel, style=Pack(width=80)))
        content.add(buttons)
        self.window.content = content
        self.window.show()

    def apply(self, widget=None, **kwargs) -> bool:
        applied, message = self.controller.apply_parameter_edit(
            self.value_input.value, bool(self.freeze_switch.value)
        )
        if not applied:
            self.error_label.text = message
            return False
        self.close()
        return True

    def _cancel(self, widget=None, **kwargs):
        self.close()

    def _on_close(self, window, **kwargs) -> bool:
        if self.controller.parameter_window is self:
            self.controller.parameter_window = None
        return True

    def close(self):
        if self.controller.parameter_window is self:
            self.controller.parameter_window = None
        if self.window is not None:
            self.window.close()
            self.window = None
