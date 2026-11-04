"""The atom/molecule selection controls, shared by every window that shows a MolecularViewer.

``SelectionControls`` is a mixin for a Toga window class. It builds the
widgets (``build_selection_widgets``) and handles them; the host lays them out
with its own row helpers and supplies two methods:

- ``selection_viewer()``: the MolecularViewer the controls act on, or None
  (nothing loaded yet);
- ``set_selection_status(text)``: show a one-line message;

plus the project's usual ``async _show_error(title, message)``. The selection
itself lives in the viewer (``viewerSelection`` has the pure helpers); the
3D Molecular Viewer's Display tab, the Fukui Function Analysis window and the
All bond / angle / dihedral analysis windows (``trajectorySelection``, which
passes its Solute atom indices field as ``selection_input``) use these
controls, so all of them behave and read the same.
"""

import toga
from toga.style import Pack

from indexRanges import format_index_ranges
from viewerSelection import selection_summary
# Write-only-on-change (every widget write re-lays out the window); the same
# function as cubeTab._assign, without importing cubeTab's numba/OpenGL chain.
from widgetWrites import assign_if_changed as _assign

SELECTION_HINT = (
    "In the 3D window, Ctrl+click selects or deselects an atom and Shift+click "
    "its whole molecule (the bonds drawn now). Selected atoms keep their identity "
    "in every frame. Or type 1-based atoms and ranges, e.g. 1-3, 295-297, and "
    "press Enter."
)
NOTHING_LOADED_MESSAGE = "Nothing is loaded yet: there are no atoms to select."


class SelectionControls:
    """Selected-atoms field, summary, Clear / Expand to molecules and Show only selected."""

    def build_selection_widgets(self, selection_input=None):
        """Create the widgets (not laid out). The switch reads the viewer's flag, if any.

        ``selection_input``: an existing TextInput to use as the selection
        field (the All-* analysis windows' Solute atom indices); by default a
        new one is created.
        """
        viewer = self.selection_viewer()
        if selection_input is None:
            selection_input = toga.TextInput(
                placeholder="e.g. 1-3, 295-297",
                style=Pack(flex=1, margin=(0, 8, 0, 0)),
                on_confirm=self.apply_selection_field,
            )
        else:
            selection_input.on_confirm = self.apply_selection_field
        self.selection_input = selection_input
        self.selection_summary_label = toga.Label("No atoms selected.", style=Pack(margin=(0, 0, 6, 0)))
        self.selection_clear_button = toga.Button(
            "Clear", on_press=self.clear_atom_selection, style=Pack(width=90, margin=(0, 8, 0, 0)))
        self.selection_expand_button = toga.Button(
            "Expand to molecules", on_press=self.expand_selection, style=Pack(width=170))
        # value goes to the constructor so building the controls fires no on_change.
        self.show_selected_switch = toga.Switch(
            "Show only selected atoms/molecules",
            value=bool(viewer is not None and viewer.show_selected_only),
            on_change=self.toggle_show_selected_only,
        )

    def _sync_selection_widgets(self):
        """Selection field and summary from the selection (written only on change)."""
        viewer = self.selection_viewer()
        if getattr(self, "selection_input", None) is None or viewer is None:
            return
        atoms = viewer.get_selected_atoms()
        elements = [element for element, _ in viewer.get_current_frame_data()]
        summary = selection_summary(elements, atoms)
        _assign(self.selection_input, "value", format_index_ranges(i + 1 for i in atoms))
        _assign(self.selection_summary_label, "text", summary)
        self.set_selection_status(f"Selection: {summary}")

    async def apply_selection_field(self, widget):
        viewer = self.selection_viewer()
        if viewer is None:
            await self._show_error("Selection Error", NOTHING_LOADED_MESSAGE)
            return
        try:
            atoms = viewer.parse_selection_text(self.selection_input.value or "")
        except ValueError as exc:
            await self._show_error("Selection Error", str(exc))
            return
        viewer.set_selected_atoms(atoms)
        self._sync_selection_widgets()  # normalise the typed text even if unchanged

    def clear_atom_selection(self, widget):
        viewer = self.selection_viewer()
        if viewer is None:
            return
        viewer.clear_selection()
        self._sync_selection_widgets()

    def expand_selection(self, widget):
        viewer = self.selection_viewer()
        if viewer is None or not viewer.get_selected_atoms():
            self.set_selection_status("Select atoms first, then expand them to whole molecules.")
            return
        viewer.expand_selection_to_molecules()
        self._sync_selection_widgets()

    def toggle_show_selected_only(self, widget):
        viewer = self.selection_viewer()
        if viewer is None:
            return
        viewer.set_show_selected_only(bool(self.show_selected_switch.value))
        if viewer.show_selected_only and not viewer.get_selected_atoms():
            self.set_selection_status("Show only selected is on; select atoms to apply it.")
