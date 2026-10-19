# ESP tab in the 3D Molecular Viewer — design

Date: 2026-10-19. Status: approved by the user. Branch: `feature/cube-tab` (builds on the Cube tab).

## 0. User decisions

1. The ESP tool moves from **Tools > Electrostatic Surface Potential (ESP)** into a new **ESP** tab of the 3D Molecular Viewer. The Tools button is removed.
2. **Option A**: the ESP surface, its extremum markers and its colour bar are drawn in the viewer's **own** 3D window, not in a separate ESP window.
3. New **Save ESP cube...**: V(r) in a.u. on the run's density grid, written as a Gaussian cube. A save dialog chooses the folder and name; the suggestion is `<fchk name>_ESP.cube`, in the output folder.
4. A switch (on by default) also writes the **electron-density cube on the same grid**, `<name>_density.cube` beside the ESP cube.
5. Every existing ESP capability is preserved. Fukui Functions (which also uses `ESPSurfaceViewer`) is unchanged. Nothing else changes.

## 1. Drawing code: `espSurfaceLayer.py` (new, OpenGL only)

Moved verbatim from `espViewer.py`:
- the pure colour helpers: `hue_colours`, `symmetric_range`, `esp_colours`, `gaussview_label`, `unit_factor`, `colour_bar_text`, `contrasting_text_colour`;
- the unit and marker constants;
- `SurfaceMarker` and `SurfaceDrawData`;
- the GL code of `_draw_surface`, `_draw_markers` and `draw_screen_overlay`.

New: `ESPSurfaceLayer` holds the surface and its display options under its own lock:
- `set_surface(..., centre)`, `clear()`, `get()`, `set_display_options()`, `colour_bar_text()`, `radius()`;
- `draw_scene(viewer)` and `draw_overlay(viewer)`. These take the viewer only for `draw_sphere`, `_draw_text_3d`, `_draw_text_2d`, `display` and `background_color`.

The module imports neither Toga, GLFW nor `molecularViewer`, so `molecularViewer` can import it without a cycle.

`espViewer.py` re-exports every moved name, so `espViewer.hue_colours`, `from espViewer import SurfaceMarker` and similar keep working. `ESPSurfaceViewer` keeps its public API (`set_molecule`, `set_surface`, `clear_surface`, `set_atom_view`, `set_display_options`, `get_surface`, `colour_bar_text`, `get_scene_radius`, `request_close`, and the `colour_unit`/`surface_opacity`/`show_extrema` attributes as properties) and delegates to an `ESPSurfaceLayer`. The existing `test_espViewer.py` must pass unchanged; it pins the Fukui viewer.

## 2. Viewer integration (`molecularViewer.py`, `MolecularViewerUI` only)

- `self.esp_layer = ESPSurfaceLayer()`.
- `render_scene_extras` draws the cube layer, then `esp_layer.draw_scene(self)`.
- New `draw_screen_overlay` override: `esp_layer.draw_overlay(self)`.
- `get_scene_radius` also covers the ESP surface: `radius + 0.3`, as in `ESPSurfaceViewer`.
- The ESP tab comes after Cube: Display, Dynamics, Measure, Design, Cube, ESP, Box & Performance.
- The base-class hooks stay empty.

## 3. The tab: `ESPTabController` (in `electrostaticSurfacePotential.py`, replaces `ElectrostaticSurfacePotentialUI`)

- **Widgets and handlers keep their names and behaviour.** The old window's header, Setup and Results sections become sections of one scrolling page:
  - fchk and output folder;
  - surface definition and mapped property;
  - Run / Show surface / Clear, plus progress and status;
  - Results: surface shown, colour-bar unit, opacity, extrema markers, canvas background, Export image;
  - Export cubes;
  - the summary text.

  `Close` disappears (the viewer window owns the tab), and so does `tabs.current_tab = 1`.
- **Structure.** On **Show surface** (and after a run while the 3D window is open):
  - If the viewer already shows exactly this fchk molecule, or shows nothing, the molecule is loaded or kept silently.
  - Otherwise a `QuestionDialog` asks whether to replace the displayed structure. **No** keeps it, and the surface is drawn in the fchk's coordinates.
  - The answer is remembered for the current results; a new run asks again.
- **Imports.** The module no longer imports `espViewer` or, at module level, `molecularViewer` (`save_canvas_image` is imported inside `export_image`), so `molecularViewer` can import it.
- **The extrema file** is still written automatically after every run.

## 4. ESP cube (headless, `ESPSurfaceAnalyser`)

- `ESPResults` gains additive fields: `grid_origin`, `grid_shape`, `density_grid` (the grid `run` already computes) and `cubegen_path` (the resolved cubegen when that engine ran).
- `esp_on_grid(results, progress)` → `(values (shape), notes)`:
  - built-in engine with progress;
  - the cubegen engine when the run used it, falling back to built-in with a note on failure;
  - a grid point closer than `ON_NUCLEUS_DISTANCE = 1e-6` bohr to a nucleus is evaluated `NUCLEUS_SHIFT = 1e-3` bohr away along +x, and a note says so. The file never holds `inf`/`nan`.
- `esp_cube(results, values)` and `density_cube(results)` → `cubeFile.CubeData`: atomic numbers, nuclear charges, bohr coordinates, the run's origin and spacing.
- `write_cubes(results, esp_path, with_density, progress)` → `(paths, notes)` through `cubeFile.write_cube`.
  - Density path: `<esp stem without a trailing _ESP>_density.cube` in the same folder.
  - `default_esp_cube_name(results)` = `<fchk stem>_ESP.cube`.
- **UI:** an `Also save the electron density cube (same grid)` switch (on by default) and a `Save ESP cube...` button.
  - The `SaveFileDialog` is suggested `<output folder>/<fchk stem>_ESP.cube`.
  - The work runs in `asyncio.to_thread` with the progress bar; buttons are disabled meanwhile.
  - Errors appear as dialogs.

Measured (built-in engine with numba, 0.25 bohr): ethylene, 134 k points, 1.3 s; benzene, 256 k points, 13.5 s.

## 5. Main menu

`gqteaWinToga.py`: the ESP import and the Tools entry are removed.

## 6. Tests

- `test_espViewer.py` unchanged.
- New `test_espSurfaceLayer.py`:
  - the layer on its own: set, clear, options, radius, caption;
  - no Toga/GLFW import (subprocess).
- `test_electrostaticSurfacePotential.py`:
  - the UI tests are ported to `ESPTabController` (the same assertions, minus the inner tab switch);
  - a stub viewer with `esp_layer`, `frames` and `load_designed_molecule`;
  - the structure question (asked only for a different structure);
  - cube export: shape, origin and spacing read back with `cubeFile.read_cube`; density values equal `density_grid`; ESP values equal a direct evaluation at sampled grid points to the `%13.5E` precision; default names; the density switch; the nucleus guard (a grid built with a nucleus on a grid point stays finite and gives a note).
- `test_molecularViewer.py`: tab order includes ESP; the UI draws the ESP layer and its overlay; the scene radius includes the ESP surface.
- A test that `gqteaWinToga.py` no longer lists the ESP tool.
- Full suite against the 4 known failures; a live smoke run of the tab.

## 7. Docs

USER_MANUAL (the ESP section now describes the viewer tab and the cube export), `help_esp_surface`, and CLAUDE.md (ESP entry, module chain, layer, tab).
