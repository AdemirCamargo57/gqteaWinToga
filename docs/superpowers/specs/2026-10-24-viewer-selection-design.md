# 3D Molecular Viewer: atom/molecule selection and selective XYZ export

Date: 2026-10-24. Status: approved by the user.

## 0. Purpose and user decisions

Inspect intermolecular interactions, especially hydrogen bonds evolving
along an MD trajectory, by selecting atoms or whole molecules on the 3D
canvas, optionally showing only them, and exporting only them as XYZ. The
trajectory data are never modified.

Decisions (user, 2026-10-24):

1. **Molecule** = connected component of the bonds *drawn* for the displayed
   frame (the same Max bond length and frame-0 / per-frame switch the user
   sees).
2. **Connectivity changes**: a molecule selection is frozen into atom
   indices at click time; those atoms are tracked and exported in every
   frame (a hopping proton stays selected).
3. **Gesture**: Ctrl+click toggles one atom, Shift+click toggles its
   molecule; a plain click keeps today's label/measure behaviour.
4. **No PBC** in molecule detection (the viewer's bonds are non-periodic);
   a molecule split by the cell boundary needs a second Shift+click.

## 1. Selection state (MolecularViewer, headless)

- `selected_atoms: Set[int]`, 0-based, guarded by `_state_lock`;
  `_selection_version` (int, bumped on every change; render caches key on it).
- `toggle_selected_atom(i)`; `toggle_selected_molecule(i)`: component of `i`
  in `self.bonds` via `scipy.sparse.csgraph.connected_components`; if the
  whole component is already selected it is removed, otherwise added.
- `set_selected_atoms(indices)`, `clear_selection()`,
  `expand_selection_to_molecules()`, `get_selected_atoms()` (sorted list).
- `_notify_selection_changed()`: no-op hook on the base class;
  `MolecularViewerUI` marshals it onto the Toga loop (same pattern as
  `_notify_picked_atoms_changed`).
- Loading a trajectory or a designed molecule clears the selection.
- Indices outside a frame's atom range are ignored when rendering.

## 2. Canvas gestures

`on_mouse_button` records `mods` at the left press. On a release without a
drag: Ctrl → `toggle_selected_atom`, Shift → `toggle_selected_molecule`,
neither → the existing `_handle_atom_click`. Ctrl takes precedence over
Shift. A click on empty space does nothing. Dragging is unchanged.

## 3. Rendering

- Selected atoms keep their element colour and get an **orange outline**
  (a larger sphere drawn with front faces culled, the Design-style outline
  technique); in Line style a small orange sphere; in Design style an orange
  ring/outline. Distinct from measurement yellow and pick/label cyan.
- **Show only selected** (`show_selected_only`, effective only with a
  non-empty selection): `render_frame` draws only selected atoms and bonds
  with both ends selected; indices stay the original ones, so labels,
  measurements and pick tags keep working. Outlines are suppressed in this
  mode. Hidden atoms are not pickable.
- In that mode the scene centre and radius are computed from the visible
  atoms (the view follows the selection); the scene-radius cache key gains
  the selection version.
- The visible-bond list is filtered per frame with set lookups (well under
  1 ms for thousands of bonds); a cache is added only if measurement shows
  a need.
- The pick-order tags and click labels of hidden atoms are not drawn; the
  measurement overlay is unchanged.

## 4. Export

- A Toga-free module-level writer: for each frame in order, the selected
  atoms in ascending original index; line 1 = selected atom count; line 2 =
  `Frame k of <source>: n selected atoms (source atoms 1-3,295-297)`; atom
  lines in the existing `{element:<3s} {x: .8f} {y: .8f} {z: .8f}` format.
- The **Export only selected atoms/molecules** checkbox
  (`export_selected_only`) applies to the header's *Save Current Frame XYZ*
  and the new *Save Trajectory XYZ...* (all loaded frames; a single-frame
  structure is a one-frame trajectory; with no structure loaded -- only the
  viewer's built-in placeholder -- it refuses). Unticked, both export every
  atom (the current-frame output is byte-identical to before).
- The trajectory export refuses an output path equal to the loaded file.
- Pure helpers live in a new Toga- and OpenGL-free module `viewerSelection.py`
  (subprocess test); `indexRanges.format_index_ranges` is the inverse of the
  existing range parser.
- Empty selection with the checkbox on → dialog, nothing written.
- A frame whose atom count does not cover the largest selected index →
  refused before writing, naming the first such frame.
- The trajectory writer runs in `asyncio.to_thread` from a snapshot of the
  frames list, with throttled status-line progress; the button is disabled
  while it runs. A write error is a dialog.

## 5. Interface (Display tab, new "Selection" section)

Gesture hint; a selection field showing/accepting 1-based ranges
(`parse_solute_index_ranges`, Enter applies; display compressed as
`1-3, 295-297`); a summary label (`6 atoms: H4 O2`); **Clear** and
**Expand to molecules** buttons; **Show only selected atoms/molecules**
switch; **Export only selected atoms/molecules** checkbox;
**Save Trajectory XYZ...** button. The existing Display-tab button
"Clear selection" (which clears measurement picks and labels) is relabelled
**Clear labels**; its behaviour is unchanged.

Display rather than Dynamics: selection applies equally to single
structures and trajectories; Dynamics is playback only.

## 6. Performance

Per frame only Python filtering of visible atoms/bonds is added, cached
while the selection is unchanged. `_pick_atom_at` (a per-atom Python loop on
the render thread) is timed on a large box and vectorised only if a click
costs more than ~10 ms. Export is timed on a synthetic large trajectory;
Numba/parallelism only if measurement justifies it.

## 7. Testing

`tests/test_viewerSelection.py`, headless: atom/molecule toggling;
molecule frozen across frames whose connectivity changes; Ctrl/Shift routing
with plain clicks unchanged; hidden atoms unpickable; visible atom/bond
lists; export counts, comment line, atom and frame order, single frame,
empty selection, atom-count mismatch, full export unchanged; range field
parsing; selection cleared on load; UI plumbing via `__new__` + stub widgets.

## 8. Not implemented

Frame-range export; selection that follows connectivity frame by frame;
PBC-aware molecule detection.
