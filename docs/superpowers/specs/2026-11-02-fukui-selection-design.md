# Fukui window: show only selected atoms/molecules — design

Date: 2026-11-02. Status: approved by the user in chat (approach A + 1).

## Goal

In **Tools > Analysis of Fukui Functions**, let the user show only some atoms
or molecules of a multi-molecule system (e.g. the acetone of
`C:\AJCamargo\TESTE\AguaPlusAcenta_*.fchk`, 1 acetone + 42 waters) together
with *their* part of the Fukui-function display, and export an image of
exactly that. The selection is a **visualization filter**: results, report,
cubes and the extrema file are never changed.

## User decisions

1. **Surface rule — nearest atom.** A surface triangle is drawn when the
   nearest atom of each of its three vertices is selected (the rule the report
   already uses for an extremum's "nearest atom"). Applies to the ±isovalue
   lobes and to the ρ(N)-mapped surface. Extremum markers are kept when their
   nearest atom is selected; their labels stay the extrema-file row numbers.
2. **Colour range — whole system.** The mapped display's ±max|f| is the full
   surface's, so filtered and unfiltered images share one colour scale. (Lobe
   displays already use the fixed ±isovalue.)
3. **Printing — image only.** *Export image...* captures the filtered canvas;
   the suggested name gains `_selected`. No filtered report or cubes.

## Design

- **Filtering lives in `espViewer.ESPSurfaceViewer`** (approach A). `set_surface`
  takes optional `vertex_atoms` (owner atom per vertex) and `marker_atoms`
  (owner per marker); the viewer keeps the full surface and re-derives the
  drawn sub-mesh in `_refocus_camera` — the base viewer's hook for "the set of
  shown atoms changed", called by every selection change while *Show only* is on
  and by `set_show_selected_only`. A Ctrl/Shift-click in the GL window
  therefore refilters on the render thread with no Toga round trip.
  `half_range` is resolved from the full values before filtering (decision 2).
  The scene radius is measured from the camera's centre (the selection's when
  only it is shown). `ESPSurfaceLayer` is unchanged, so the 3D Molecular
  Viewer's ESP tab is unaffected.
- **`set_molecule` keeps the selection** when the structure is the one already
  shown (`volumeData.same_structure`): `_push_to_viewer` runs on every display
  change, and a reload clears the selection.
- **Pure helpers** in the Toga-free `viewerSelection.py`: `nearest_atoms(points,
  atom_positions)` (cKDTree) and `submesh(triangles, keep)` (triangles whose
  three vertices are kept, renumbered onto the vertices they use).
- **Shared controls** (approach 1): the Display tab's selection handlers move to
  a mixin `viewerSelectionControls.SelectionControls` (field, summary, Clear,
  Expand to molecules, Show only selected), used by `MolecularViewerUI` and
  `FukuiFunctionAnalyzeUI`. The host supplies `selection_viewer()` and
  `set_selection_status(text)` and lays the widgets out with its own helpers.
  The Display tab keeps its attribute and method names; its tests pass
  unchanged. No "Export only selected" switch in the Fukui window (it writes
  no XYZ).
- **Fukui window:** a *Selection* section on the Results tab. The viewer object
  is created as soon as results exist (the GL window still opens only on
  *Show 3D*), so a selection can be typed first. Owner atoms are computed per
  displayed mesh and cached until the next run/clear. The viewer's selection
  hook is forwarded to the window (`selection_listener`, marshalled with
  `loop.call_soon_threadsafe`) to update the field, summary and status, which
  says when no surface of the current function lies near the selection.

## Testing

Headless: helper tests; `ESPSurfaceViewer` filtering (every drawn vertex owned
by a selected atom, the dropped triangles owned by others, markers, whole-system
half-range, show off = the full surface, selection kept across a same-structure
`set_molecule`); Fukui UI plumbing with stub widgets (phenol fixture); the
existing Display-tab selection tests unchanged. Validation on the TESTE files:
run the analysis headless, check the acetone filter numerically, then live in
the window: select atoms 1-10, Show only selected, compare the exported PNG with
a frame capture of the canvas.
