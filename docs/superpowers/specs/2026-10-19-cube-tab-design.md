# Cube tab for the 3D Molecular Viewer — design

Date: 2026-10-19. Status: approved by the user (decisions §0). Branch: `feature/cube-tab`.

## 0. User decisions

1. **Structure ownership.** The first cube loaded while a structure is displayed asks (Yes/No dialog) whether to replace the displayed structure with the cube's own atoms; with nothing displayed it replaces it silently. Later cubes never change the structure; if their atoms differ from the displayed ones the status line says so.
2. **Slices.** A coloured plane inside the 3D scene is essential. A 2D map in the interactive plot viewer is optional and not part of this work (it would add a figure type to `displayPlots`/`plotViewer`/`plotter`).
3. **Phasing.** Phase 1 = cube I/O, multi-cube management, ± isosurfaces, visibility, grid operations, cube and image export. Phase 2 (same branch) = the slice plane.
4. **Numba.** The heavy grid kernel is compiled with numba. As everywhere in the project, numba stays **optional at run time**: a numpy path gives the same mesh, and every import is guarded.
5. **Scope rule.** All existing behaviour is preserved; changes to existing modules are limited to what the tab needs. No unrelated refactoring.

## 1. What exists and is reused

| Existing | Reused for |
|---|---|
| `cubeFile.py` (strict bohr, axis-aligned, one spacing; used by the Fukui tool) | Kept byte-for-byte; the general reader/writer is **added** beside it |
| `molecularSurface.marching_tetrahedra` (vectorised, shared vertices) | Isosurfaces, run in grid-index space (origin 0, spacing 1) and then mapped through the cube's axes — exact for skewed and unequal-spacing grids, since the linear interpolation is affine-invariant |
| `MolecularViewer` hooks `render_scene_extras` / `draw_screen_overlay`, `get_scene_radius` | Drawing the volume layer; overridden in **`MolecularViewerUI`**, never in the base class (which must keep them empty) |
| `MolecularViewer.capture_frame`, `save_canvas_image` | Export image |
| `MolecularViewer.load_designed_molecule(frame)` | Adopting a cube's atoms as a single-frame structure |
| `MolecularViewer.ensure_viewer_window()` | Opening/focusing the GL window after a load |
| `help.AtomicData.atomic_numbers` | Element symbols from the cube's atomic numbers |
| `_tab_page`, `_form_row`, `_form_label`, `_section_heading`, `_hint` | Tab layout |

## 2. Modules

One-way chain: `molecularViewer → cubeTab → {cubeRender, volumeData → (volumeNumba), cubeFile, molecularSurface}`.

### 2.1 `cubeFile.py` — additive

- `CubeVolume` dataclass: `title`, `comment`, `atomic_numbers (n,)`, `nuclear_charges (n,)`, `coordinates (n,3) bohr`, `origin (3,) bohr`, `axes (3,3) bohr` (row k = step vector of grid index k), `values (nx,ny,nz) float64`, `orbital: Optional[int]`.
- `read_cube_volumes(path) -> List[CubeVolume]`:
  - negative atom count → orbital cube: the line after the atoms holds `NMO` and the MO numbers (they may wrap over several lines); values are stored MO-fastest per point and are split into **one `CubeVolume` per orbital**;
  - a negative voxel count on lines 4-6 → Å: origin, axes and atom coordinates are converted to bohr with `BOHR_RADIUS_ANGSTROM = 0.529177210544` (repeated here as in `hirshfeldCharges`; a test pins it equal to `coordinateConverter`'s);
  - any axes (skewed, unequal spacing); a singular axis matrix is refused;
  - values parsed with `np.array(text.split(), dtype=float64)` (measured 1.2 s for 8 M values; `np.fromstring(sep=" ")` took 14 s). If that fails, Fortran exponents written without `E` (`1.234-100`) are repaired by a regex and parsing is retried; anything else is a `ValueError` naming the file;
  - a value count that does not match the grid → `ValueError` ("looks truncated").
- `write_cube_volume(path, volume)`: bohr, positive atom count (an orbital is written as a plain cube, its MO number on the comment line), the same `%13.5E`, 6-per-line, z-fastest layout as `format_cube`. `OSError` → `ValueError`.
- `CubeData`, `read_cube`, `write_cube`, `format_cube` are untouched.

### 2.2 `volumeData.py` — no Toga, no GL

- `voxel_volume(v) = |det(axes)|`; `integral(v) = Σ values · dV`; `value_statistics(v)` (min, max, integral).
- `grid_mismatch(a, b) -> Optional[str]`: `None` when shapes are equal and origin and axes agree within `GRID_TOLERANCE = 1e-4` bohr (headers carry 6 decimals); otherwise a sentence naming the cause. `atoms_mismatch(a, b) -> Optional[str]` for the warning.
- `combine(a, b, op)` with `op in ("add", "subtract")` and `scale(a, factor)`: new `CubeVolume` (atoms and grid of `a`, title describing the operation). Incompatible grids → `ValueError` with the `grid_mismatch` sentence. Non-finite factor refused.
- `default_isovalue(v)`: orbital cube or values of both signs → 0.02; all values ≥ 0 → 0.0004 (GaussView's density default); then capped below `max|v|` (half of it) and rounded to one significant figure.
- `isosurface(v, isovalue, sign=+1, use_numba=None) -> Optional[IsoMesh]`: marching tetrahedra on `sign · values` at `+isovalue` in index space, vertices mapped `origin + idx @ axes`, triangles flipped when `det(axes) < 0` so orientation is preserved, area-weighted vertex normals. Returns `None` when the isovalue cuts nothing. `IsoMesh(vertices bohr, triangles, normals)`.
- `slice_plane(v, axis, index)` → `(points (n1,n2,3) bohr, values (n1,n2))` for the lattice plane `index` along grid axis 0/1/2; `slice_mesh(points)` → triangles of the plane; `slice_table_text(v, axis, index)` → a commented text table (`x y z value`, Å) for export.
- Colour helpers (pure): `parse_colour(text)` (preset names, `#rrggbb`) → RGB floats; `diverging_colours(values, half_range)` (blue–white–red, used for the symmetric range) and `rainbow_colours(values, lo, hi)` (red→green→blue, GaussView hue sweep, for min–max). They are re-implemented here rather than imported because `espViewer` imports GLFW and `molecularDesign.parse_hex_color` lives in a Toga module.
- `numba_available()` and the `use_numba=None` convention of `gaussianWavefunction` (None = if installed; True without numba → `ValueError`).

### 2.3 `volumeNumba.py` — the only new numba importer

- `tetra_cut_edges(values, isovalue) -> (edge_ids (T,3) int64, outward (T,3) float64)`: a two-pass `prange` scan over cells (pass 1 counts the triangles per x-slab, prefix sum, pass 2 fills), reproducing exactly the polygons, edge ids (`lin*7 + k`) and outward directions of `molecularSurface.marching_tetrahedra`.
- Importing it raises `ImportError` without numba; `volumeData` falls back.

### 2.4 `molecularSurface.py` — minimal split

`marching_tetrahedra` is split into its two halves with **no change of behaviour**: `_tetra_cut_edges_numpy(values, isovalue)` (the cell scan) and `mesh_from_cut_edges(tri_edges, outward, values, origin, spacing, isovalue, density, n_bisect)` (decode, interpolate, orient). `marching_tetrahedra` calls both, so the ESP and Fukui tools get identical meshes (their tests, which compare against cubegen/Multiwfn, guard it). The cube tab calls `mesh_from_cut_edges` after either cell scan, so there is **one** finishing implementation. Because `np.unique` sorts the edge ids, the vertex arrays of the numpy and numba paths are identical; triangles may come in a different order.

### 2.5 `cubeRender.py` — OpenGL, no Toga, no GLFW

- `RenderItem` (frozen): `kind` (`"surface"` | `"slice"`), `vertices` (Å, float32), `normals` (surface only), `indices` (uint32), `rgba` (per vertex, float32), `opacity`.
- `CubeSceneLayer`: `set_items(items)` / `clear()` under its own lock (Toga thread), `draw()` (render thread) and `radius_about(centre)`. `draw_order(items)` (pure, tested): opaque before translucent. Surfaces: two-sided lighting, `GL_COLOR_MATERIAL`, specular, blending with depth writes off when translucent — the `espViewer._draw_surface` recipe; slices: lighting off, vertex colours. All GL state restored in `finally`.

### 2.6 `cubeTab.py` — Toga

- `CubeEntry`: `name`, `volume`, `path` (None for an operation result), `saved`, and display settings `isovalue`, `show` (`"+ and −"` / `"+ only"` / `"− only"`), `positive_colour`, `negative_colour`, `opacity`, `visible`, plus a mesh cache keyed by `(isovalue, sign)`.
- `CubeTabController(viewer)`: owns all state and touches no widget in its headless methods, which are the test seams:
  - `load_file(path) -> (entries, notes)` (sync, run in `asyncio.to_thread`): names are unique (`stem`, `stem#MO12`, `stem (2)`).
  - `adopt_structure(entry)`: cube atoms → `viewer.load_designed_molecule(frame)` (Å, symbols from atomic numbers; Z without a symbol → `X`). A cube without atoms is refused when no structure is displayed.
  - `remove(name)`, `apply_operation(a, op, b, factor, result_name)`, `save_entry(name, path)`, `slice_text()`.
  - `scene_items()` → the `RenderItem` list for all visible entries and the slice; computing meshes runs in a thread; a generation counter discards a result overtaken by a newer change (the design tab's `_edit_generation` pattern).
- `build_tab()` → the tab's `toga.Box`; async handlers wrap the seams and report errors through `toga.ErrorDialog` / status lines. Nothing raises out of a handler.

### 2.7 `molecularViewer.py` — integration only

1. `import cubeTab` and one `("Cube", …)` entry after Design, before Box & Performance (the other five tabs keep their order).
2. `_build_cube_tab()` creates `self.cube_controller` and returns its tab.
3. `MolecularViewerUI` overrides `render_scene_extras` (draws `self.cube_layer`), and `get_scene_radius` (max of the base radius and the layer's). The base-class hooks stay empty.
4. Base class: `show_structure = True`; `main_loop` calls `render_frame` only when it is True, and `_pick_atom_at` returns `None` while it is False. With the default nothing changes.

## 3. The tab

Sections, top to bottom (inside the usual scroll container):

1. **Cube files** — `Load cube(s)...` (multi-select `.cube/.cub`), `Remove`; a `toga.Table` (Name, Grid, Min, Max, Integral) whose selection drives the next section; `Show molecule` and `Show volumetric data` switches.
2. **Isosurface (selected cube)** — Isovalue (Enter applies), Show (`+ and −` / `+ only` / `− only`), positive and negative colour (preset name or `#rrggbb`, Enter applies; defaults Blue / Red), opacity slider (0.1-1), Visible switch.
3. **Slice** — Cube, Axis (a/b/c — x/y/z for an orthogonal grid), index slider with its value label, Colour range (`Symmetric ±max|V|` / `Min to max`), Show slice switch, `Export slice data...`.
4. **Grid operations** — A, operation (`A + B`, `A − B`, `c × A`), B, factor c, result name (blank → automatic), `Compute`. The result appears in the table as unsaved.
5. **Export** — `Save cube...` (selected entry), `Export image...` (captures before asking for the file name, as the ESP tool does).

Loading opens/focuses the 3D window (`ensure_viewer_window`).

## 4. Errors and limits

- Every file/grid problem is a `ValueError` with a message naming the file or the grids; the tab shows it in a dialog.
- "No surface at this isovalue" is a status message, not an error.
- Known limits: translucent surfaces are not depth-sorted against each other; volumes are static during trajectory playback; very large grids (> ~250³) are slow to mesh and are held fully in memory; no PBC images of the volume.

## 5. Tests

- `test_cubeFile.py` (new cases): orbital cube with 1 and 3 MOs (wrapped MO line), Å header, skewed axes, unequal spacing, Fortran exponent repair, truncated data, singular axes, write → read round trip, Bohr constant equal to `coordinateConverter`'s; the existing tests unchanged.
- `test_volumeData.py`: Gaussian blob sphere radius within one grid step; p-like function → mirror-image ± lobes; same function on a skewed grid → same surface radius; `det(axes) < 0` keeps orientation outward; `grid_mismatch` matrix; `A − A = 0`, linear scaling, integral of a normalised Gaussian ≈ 1; default isovalues; slice indexing and table; colour parsing; numba vs numpy identical vertices and identical sorted triangles (skipped without numba); `use_numba=True` without numba → `ValueError`; Toga-free subprocess test.
- `test_molecularSurface.py`/ESP/Fukui tests unchanged and passing (guards the split).
- `test_cubeRender.py`: `draw_order`, item construction, radius; Toga-free subprocess test.
- `test_cubeTab.py`: controller headless with a stub viewer (load, unique names, adopt structure, remove, operations incl. refusal, save, slice text, scene items honour visibility switches, stale-generation discard) and the widget plumbing via `__new__` + stub widgets.
- `test_molecularViewer.py`: base hooks still no-ops; `show_structure` default True; picking returns None while hidden; tab order.
- Fixtures: `tests/data/cube/` — small cubegen cubes of ethylene (HOMO orbital cube and total density) produced by `tests/data/cube/make_reference.py` (needs G09W; never run by the suite).
- Full suite compared with the 4 known baseline failures; a manual run of the tab in the live app.

## 6. Optional enhancements (not in this work)

2D contour map in the plot viewer; mapping one cube onto another's isosurface (ESP on density); product, |A|, A², linear combination, resampling onto a common grid; gradient-based normals; grid coarsening for very large cubes; background-colour selector; generating cubes from `.fchk` through cubegen.
