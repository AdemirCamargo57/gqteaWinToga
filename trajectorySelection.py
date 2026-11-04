"""Atom/molecule selection for the All bond / angle / dihedral analysis windows.

``TrajectorySelectionUI`` puts the shared selection controls
(viewerSelectionControls: Clear, Expand to molecules, Show only selected,
summary) under a window's *Solute atom indices* field and makes that field
the selection field: a Ctrl/Shift-click in the 3D view, Clear or Expand
writes the 1-based ranges into it, and typing ranges + Enter selects them.
The analysis still reads only the field's text, so a typed value works
without Enter or the 3D view, and a blank field still means all atoms.

The structure is the trajectory's first frame, loaded after Browse into a
``selectionViewer.SelectionViewer`` whose GL window opens only on Show 3D.
That module (and the OpenGL/numba chain behind molecularViewer) is imported
lazily, so these windows open as fast as before.

A host supplies ``main_window`` and calls ``build_selection_rows`` while
laying out, ``load_selection_structure`` after reading the first frame and
``close_selection_viewer`` (or ``on_window_close``) when it closes.
"""
import asyncio
import textwrap

import toga
from toga.style import Pack

from viewerSelectionControls import SELECTION_HINT, SelectionControls
from widgetWrites import assign_if_changed

SOLUTE_SELECTION_HINT = (
    "Show 3D opens the first frame of the trajectory. " + SELECTION_HINT
    + " The selection is the Solute atom indices field; blank = all atoms."
)
HINT_WIDTH = 100  # characters per hint line: the hint must not widen the window
STATUS_WIDTH = 600  # px; fixed, so a long message is cut off instead of widening the window
HINT_COLOR = "#666666"
NO_TRAJECTORY_MESSAGE = "Choose a trajectory with Browse first: its first frame is what you select atoms in."


def _selection_viewer_class():
    """SelectionViewer, imported on first use (it brings in OpenGL, GLFW and numba)."""
    from selectionViewer import SelectionViewer
    return SelectionViewer


def _running_loop():
    try:
        return asyncio.get_running_loop()
    except RuntimeError:
        return None  # headless: no loop to marshal selection changes onto


class TrajectorySelectionUI(SelectionControls):
    """Selection controls bound to a tool window's Solute atom indices field."""

    structure_viewer = None          # the SelectionViewer, once a trajectory is chosen
    selection_viewer_title = "Atom Selection"

    # --- layout -------------------------------------------------------------
    def build_selection_rows(self, solute_input, label_width: int):
        """The selection rows to add under the solute row; ``label_width`` aligns them with its field."""
        self.build_selection_widgets(selection_input=solute_input)
        self.show_3d_button = toga.Button(
            "Show 3D", on_press=self.show_structure, style=Pack(width=90, margin=(0, 8, 0, 0)))
        self.selection_status_label = toga.Label(
            "", style=Pack(width=STATUS_WIDTH, color=HINT_COLOR, font_size=9, margin=(0, 5, 4, 5)))

        def row(*children):
            indent = toga.Box(style=Pack(width=label_width + 10))  # the label and its margins
            return toga.Box(children=[indent, *children], style=Pack(direction="row", margin=(0, 0, 4, 0)))

        hint = toga.Label("\n".join(textwrap.wrap(SOLUTE_SELECTION_HINT, HINT_WIDTH)),
                          style=Pack(color=HINT_COLOR, font_size=9, margin=(0, 5, 6, 5)))
        return [
            row(self.show_3d_button, self.selection_clear_button, self.selection_expand_button),
            row(self.show_selected_switch),
            row(self.selection_summary_label),
            self.selection_status_label,
            hint,
        ]

    # --- SelectionControls host ------------------------------------------
    def selection_viewer(self):
        return self.structure_viewer

    def set_selection_status(self, text: str):
        label = getattr(self, "selection_status_label", None)
        if label is not None:
            assign_if_changed(label, "text", text)

    async def _show_error(self, title: str, message: str):
        await self.main_window.dialog(toga.ErrorDialog(title, message))

    # --- the structure ------------------------------------------------------
    def load_selection_structure(self, elements, coordinates):
        """After Browse: the first frame becomes the structure; the typed field is kept.

        Loading clears the viewer's selection, so the field's ranges are
        re-applied; ranges this file does not have leave the field as typed
        (the analysis reports them). A 3D failure never breaks Browse.
        """
        frame = [(str(e), tuple(float(c) for c in xyz)) for e, xyz in zip(elements, coordinates)]
        try:
            viewer = self.structure_viewer
            if viewer is None:
                viewer = _selection_viewer_class()()
                viewer.window_title = self.selection_viewer_title
            viewer.loop = _running_loop()
            viewer.selection_listener = None  # the load's own clear is not a selection edit
            viewer.load_designed_molecule(frame)
        except Exception as exc:
            self.set_selection_status(f"3D selection unavailable: {exc}")
            return
        self.structure_viewer = viewer
        text = (self.selection_input.value or "").strip()
        problem = None
        if text:
            try:
                viewer.set_selected_atoms(viewer.parse_selection_text(text))
            except ValueError as exc:
                problem = str(exc)
        viewer.selection_listener = self._sync_selection_widgets
        if problem is None:
            self._sync_selection_widgets()
        else:
            assign_if_changed(self.selection_summary_label, "text", "No atoms selected.")
            self.set_selection_status(f"Solute atom indices kept as typed: {problem}")

    async def show_structure(self, widget):
        viewer = self.structure_viewer
        if viewer is None:
            await self._show_error("Show 3D", NO_TRAJECTORY_MESSAGE)
            return
        viewer.loop = asyncio.get_running_loop()
        ready, error = await viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D view could not be opened", error)

    def close_selection_viewer(self):
        if self.structure_viewer is not None:
            self.structure_viewer.request_close()

    def on_window_close(self, window, **kwargs) -> bool:
        """The tool window's on_close: its 3D view goes with it."""
        self.close_selection_viewer()
        return True
