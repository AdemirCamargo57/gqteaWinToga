# QTAIM topology analysis (QTAIM tab of the 3D Molecular Viewer) — design

Date: 2026-10-20. Status: approved by the user. Branch: `feature/qtaim-topology`.

## 0. User decisions

1. **Input:** Gaussian `.fchk`/`.fch` only, read by the existing `gaussianWavefunction.read_fchk` (s–g, SP, pure and Cartesian shells). The analysed density is `Total SCF Density`, as in the ESP and Fukui tools. `.wfn`/`.wfx` and correlated densities are out of v1.
2. **Placement:** a new **QTAIM** tab of the 3D Molecular Viewer. CPs and paths are drawn in the viewer's own GL window, as the ESP and Cube tabs draw theirs (two GLFW windows do not work in this process).
3. **v1 ships both views:** the 3D view and the 2D plane map.
4. **Reference program:** Multiwfn 2026.7.15 at
   `C:\Users\ajc\Downloads\04_Software_and_Tools\Scientific_Software\Multiwfn\Multiwfn_2026.7.15_bin_Win64\Multiwfn_2026.7.15_bin_Win64\Multiwfn.exe`.
   It is used **only** by a dev-time reference script, never at run time. New small fixtures may be generated with G09W.
5. **ECP / heavy atoms:** the analysis runs, and the report carries a `NOTE:` that NCPs and nearby properties of ECP atoms are unreliable, because the `.fchk` has no core density (Multiwfn adds one from its EDF library; we do not). An ECP atom is one whose `nuclear_charges` differs from its `atomic_numbers`.
6. **Numba wherever it helps.** Every evaluation of ρ and its derivatives goes through a compiled kernel when numba is installed. Numba stays **optional**, with a numpy fallback that gives the same numbers (agreement ≤ 1e-12 relative, as for the existing density kernel). The kernels live in the existing `densityNumba.py`, so the project's list of numba kernel modules (`NUMBA_KERNEL_MODULES` in `tests/test_solventPacking.py`) does not change. After implementation, the remaining hot spots are profiled; any Python loop that still dominates gets a kernel in the same module.

## 1. Multiwfn conventions to reproduce

Read from Multiwfn 2026.7.15 itself: the topology menus −1/−2, and a water 6-31G(d) run of modes 2+3, paths, then option 7 → `CPprop.txt`.

**CP search defaults (menu −1):**
- maximum iterations 120
- step-size scale 1.0
- gradient-norm convergence 1e-6
- displacement convergence 1e-7
- minimum distance between CPs 0.03 bohr
- vdW-sum cutoff factor 1.5 for seed combinations
- singular-Hessian determinant 5e-22
- no trust radius
- Newton algorithm

In nuclear-seed mode (2), the gradient-norm criterion is waived (manual §3.14.2).

**Path defaults (menu −2):** at most 451 points, step 0.03 bohr, stop within 0.05 bohr of a CP.

**Seed modes:**
- 2 = nuclei
- 3 = midpoints of atom pairs
- 4 = triangle centres
- 5 = tetrahedron centres

In modes 3–5 a combination is skipped when any pair in it is farther apart than 1.5 × (r_vdW,i + r_vdW,j).

**Properties at a CP** (water BCP check: G = 0.06335348440, K = 0.5563370013, V = −0.6196904857, H = −0.5563370013, ∇²ρ = −1.971934067, λ = −1.787374161, −1.746326200, 1.561766294, ε = 0.023505, η = 1.144457):
- ρ = Σ_ab P_ab χ_a χ_b
- ∇ρ = 2 Σ P_ab χ_a ∇χ_b
- ∂_i∂_jρ = 2 Σ P_ab (∂_iχ_a ∂_jχ_b + χ_a ∂_i∂_jχ_b)
- **G = ½ Σ P_ab ∇χ_a·∇χ_b** (Lagrangian kinetic energy density; exact from the density matrix)
- K = G − ¼∇²ρ (Hamiltonian kinetic energy density)
- **V = −G − K** (= ¼∇²ρ − 2G, the local virial theorem)
- **H = G + V = −K** (total energy density)
- eigenvalues λ₁ ≤ λ₂ ≤ λ₃ with eigenvectors
- ellipticity ε = λ₁/λ₂ − 1
- η = |λ₁|/λ₃
- derived ratios |V|/G and H/ρ

**Type:** (rank, signature) from the eigenvalue signs, giving (3,−3), (3,−1), (3,+1) or (3,+3). A converged point with a zero eigenvalue (rank < 3) is degenerate; it is not reported as a CP and is counted in a NOTE.

**Poincaré–Hopf:** n(3,−3) − n(3,−1) + n(3,+1) − n(3,+3) = 1 for an isolated molecule. A failure becomes a `NOTE:` ("some CPs are missing"), never an error.

**Bond paths:** from each (3,−1), step forwards and backwards along the λ₃ eigenvector, then climb ∇ρ until within 0.05 bohr of a (3,−3). **Ring paths** (Multiwfn option 9): from each (3,+1), step along the λ₁ eigenvector, then descend to a (3,+3). A path that reaches 451 points without arriving is discarded and counted in a NOTE. Multiwfn reports each half-path's length separately (water O–H: 0.29034 + 1.45472 bohr).

**NCP ↔ atom:** each (3,−3) belongs to its nearest nucleus. A (3,−3) farther than `NNA_DISTANCE` (0.3 bohr) from every nucleus is a non-nuclear attractor (`NNA`). Hydrogen NCPs legitimately sit ~0.065 bohr off the nucleus with Gaussian basis sets (water, Multiwfn).

**Colours:**
- 3D GUI: magenta (3,−3), orange (3,−1), yellow (3,+1), green (3,+3), brown bond paths, cyan labels.
- 2D plane map (manual §3.5.5): brown (3,−3), blue (3,−1), orange (3,+1), green (3,+3), bold dark-brown bond paths, deep-blue interbasin paths. CPs and path points more than 0.5 bohr from the plane are not drawn.

## 2. Modules and dependency chain

```
molecularViewer ─► qtaimTab ─► { qtaimLayer, qtaimAnalysis ─► qtaimTopology ─► gaussianWavefunction ─► (densityNumba) }
displayPlots.draw_plane_map ◄─ plotViewer, plotter (validation), qtaimTab (manifest)
```

Lengths are **bohr** in the numerical modules. Å appears only in the report, the tab and the views. `BOHR_RADIUS_ANGSTROM` is repeated in the numerical modules and tested equal to `coordinateConverter`'s, as `molecularFileFormats` does, because `coordinateConverter` imports Toga.

Only `qtaimTab` imports Toga. `qtaimLayer` imports OpenGL, but neither Toga, GLFW nor `molecularViewer`. Subprocess tests keep `toga` out of `qtaimTopology`, `qtaimAnalysis` and `qtaimLayer`.

### 2.1 `gaussianWavefunction.py` (extended)

- `evaluate_cartesian_ao_derivatives(wfn, points, order, use_numba=None)` returns `(n_points, n_deriv, n_cart)`, where n_deriv = 4 for order 1 (value, x, y, z) and 10 for order 2 (+ xx, yy, zz, xy, xz, yz).
- `DensityDerivatives` dataclass: `rho (n)`, `gradient (n,3)`, `hessian (n,3,3)` (order 2 only), `kinetic_g (n)`.
- `density_derivatives(wfn, points, order=2, chunk=…, use_numba=None)`. The contraction with `cart_density` is a BLAS product, as in `matrix_densities`: χP and ∂_kχP are formed once per chunk.

### 2.2 `densityNumba.py` (extended)

`_ao_derivatives_kernel` is the existing `_aos_kernel` taken to first and second derivatives:
- radial part and its derivatives per primitive;
- the same e⁻⁵⁰ shell screening;
- `prange` over points.

The public wrapper `cartesian_ao_derivatives(points, packed, order)` reuses the `PackedShells` layout.

### 2.3 `qtaimTopology.py` (new, Toga-free)

- `CPSearchParameters`, `PathParameters`: dataclasses holding the Multiwfn defaults above.
- `CriticalPoint`: `position` (bohr), `kind` ("(3,-3)"…), `eigenvalues`, `eigenvectors`, `properties` (`CPProperties`: rho, laplacian, G, K, V, H, ellipticity, eta, |V|/G, H/ρ), `assignment` (atom index, `NNA`, atom pair, or None).
- `seeds_nuclei`, `seeds_pairs`, `seeds_triangles`, `seeds_tetrahedra`: combinations come from a `cKDTree` neighbour list at the largest cutoff and are then filtered pair by pair, never by an O(N⁴) loop.
- `newton_search(wfn, seeds, params, nuclear=False, progress)`: batched over every live seed, one `density_derivatives` call per iteration. Step = −step_scale · H⁻¹∇ρ, capped at `trust_radius` when it is set. The singular test and the iteration limit drop a seed. Points with ρ < `MIN_CP_DENSITY` (1e-8 a.u.; far-field noise) are dropped.
- `merge_critical_points(found, existing, min_distance)`: deterministic; the first-found point wins, in Multiwfn's mode order.
- `classify`, `cp_properties`, `poincare_hopf`, `assign_nuclear_cps`.
- `trace_paths(wfn, starts, direction (+1 ascend / −1 descend), targets, params, plane=None)`: batched RK4 on ∇ρ/|∇ρ| with the fixed step. A path stops within `stop_distance` of a target CP. `plane` restricts the gradient to the plane, for the 2D lines.
- `bond_paths`, `ring_paths`: these build the starts from the eigenvectors and record each half-path's points and length.

### 2.4 `qtaimAnalysis.py` (new, Toga-free)

- `QTAIMAnalyser.run(fchk_path, modes, search, paths, progress)` → `QTAIMResults`. The results hold the wfn, CPs (ordered (3,−3), (3,−1), (3,+1), (3,+3), each in discovery order), bond paths, ring paths, notes and parameters.
- Every problem is a `ValueError` with a specific message.
- `progress(fraction, message)` may raise; the exception is swallowed, as elsewhere.
- `report_text(results, written)` and `write_outputs(results, output_dir)` write `<stem>_QTAIM.txt`. The UI writes it automatically after every run.
- `plane_map(results, plane, options, progress)` → `PlaneMap`. It holds:
  - the grid of ρ;
  - Bader contour levels 1, 2, 4, 8 × 10ⁿ from 0.001 to 20 a.u.;
  - in-plane gradient lines from seed rings around each nucleus in the plane;
  - interbasin lines from each in-plane BCP along its in-plane negative-curvature direction;
  - CPs and path points within 0.5 bohr of the plane, projected into plane coordinates (Å).

  `plane_manifest(map)` returns the plot-viewer manifest entry.

### 2.5 `qtaimLayer.py` (new, OpenGL only)

`QTAIMSceneLayer` follows the pattern of `ESPSurfaceLayer`/`CubeSceneLayer`:
- own lock;
- immutable draw items in Å;
- `set_scene`, `clear`, `get`, `set_display_options` (per CP type, bond paths, ring paths, labels, highlighted CP), `radius()`, `draw_scene(viewer)`.

It uses only `draw_sphere`, `_draw_text_3d` and GL lines.

### 2.6 `qtaimTab.py` (new, Toga)

`QTAIMTabController(viewer).build_tab()` builds one scrolling page:
- **File:** `.fchk` + Browse; optional output folder (blank = next to the file).
- **Search:** switches for modes 2–5, all on by default.
- **Advanced:** every `CPSearchParameters`/`PathParameters` field; a blank field means the default.
- **Run** + progress bar + status. Real progress, marshalled with `loop.call_soon_threadsafe`.
- **Results:** the counts, the Poincaré–Hopf line and a CP table (row, type, assignment, ρ, ∇²ρ, ε). Selecting a row highlights that CP.
- **Display:** switches per CP type, bond paths, ring paths, labels; Display button.
- **2D map:**
  - Plane selection: "Three atoms" (default: the first three non-collinear atoms, editable), XY, XZ or YZ with an offset in Å.
  - Margin (2 Å) and grid points (200).
  - Switches for contours, gradient lines and interbasin paths.
  - **Plot 2D map**.

Rules carried over from the ESP and Cube tabs:
- `_settle_structure` behaves as in the ESP tab.
- Widget writes go through write-only-when-changed.
- Handlers are synchronous and start their async work with a guard.
- A test asserts every `run()` and `plane_map()` parameter is supplied by `_collect_inputs`.

### 2.7 Integration

- `molecularViewer.py`: a `("QTAIM", …)` tab after ESP; `_build_qtaim_tab`; `MolecularViewerUI.render_scene_extras`/`get_scene_radius` also draw and size `qtaim_layer` (guarded with `getattr`). The base-class hooks stay empty.
- `displayPlots.draw_plane_map(axes, entry, font_style=None)`: one renderer for the `{"type": "plane_map"}` manifest shape, used by `plotViewer.py` and the static fallback. Its import stays above the `TkAgg` pin.
- `plotter.py`: `validate_plot_manifest`, `describe_json_figures` and `json_plot_static` learn `plane_map`, so the replot tool reopens what the project writes.
- The manifest is written to the output folder as `<stem>_QTAIM_map.json` and opened through `launch_plot_viewer`. No PNG is written (`save_png=False` recipe).
- `help.py`: `HelpGqteaWin.help_qtaim`. `USER_MANUAL.md`: *QTAIM tab* (workflow + formulas + limits). `CLAUDE.md`: tool entry.

### 2.8 Manifest shape `plane_map`

```
{"type": "plane_map", "title", "xlabel", "ylabel",
 "extent": [xmin, xmax, ymin, ymax],                 # Å, plane coordinates
 "grid": [[...], ...] or null, "contours": [levels] or [],
 "lines": [{"kind": "gradient"|"interbasin"|"bond_path"|"ring_path", "points": [[x, y], ...]}],
 "points": [{"type": "(3,-3)"|..., "x", "y", "label"}]}
```

## 3. Report `<stem>_QTAIM.txt` (Fukui layout)

```
# gQTEA QTAIM Topology Analysis
# file                   <name>.fchk
# method / basis         RB3LYP / 6-31G(d)
# density                Total SCF Density
# electrons              10   charge 0   multiplicity 1
# search_modes           nuclei, atom pairs, triangles, tetrahedra
# cp_search              max_iterations 120; step_scale 1; gradient 1e-06; displacement 1e-07;
#                        min_cp_distance 0.03 bohr; vdw_cutoff_scale 1.5; singular 5e-22; trust_radius none
# path_generation        max_points 451; step 0.03 bohr; stop_distance 0.05 bohr
# units                  positions and lengths in Angstrom; properties in a.u.

CRITICAL POINT SUMMARY
  (3,-3) 3   (3,-1) 2   (3,+1) 0   (3,+3) 0
  Poincare-Hopf  3 - 2 + 0 - 0 = 1   satisfied

CRITICAL POINTS
   row  type      x (A)      y (A)      z (A)   assignment
DENSITY PROPERTIES AT CRITICAL POINTS (a.u.)
   row  type   rho  laplacian  lambda1  lambda2  lambda3  ellipticity  eta
ENERGY DENSITIES AT CRITICAL POINTS (a.u.; V = -G - K, H = G + V = -K)
   row  type   G  K  V  H  |V|/G  H/rho
BOND PATHS
   path  BCP  atoms   length_1 (A)  length_2 (A)  total (A)  distance (A)  ratio
RING PATHS
   path  RCP  CCP   length (A)

NOTE: ...
Written: ...
```

Properties use `%.10E`-style precision, as Multiwfn does.

## 4. Testing

- **Derivatives:** finite-difference checks of ∇ρ, ∇∇ρ and G against `electron_density` and the AO values, on water 6D/5D, pure f/g (water cc-pVTZ/QZ) and Cartesian 10F. The numba and numpy paths must agree.
- **Topology:**
  - H₂ BCP at the midpoint by symmetry;
  - Poincaré–Hopf on every fixture;
  - seed counts with the vdW cutoff;
  - deterministic merging;
  - RK4 path end-points;
  - plane projection;
  - ECP NOTE (stubbed nuclear charges).
- **Multiwfn reference:** `tests/data/qtaim/make_reference.py` drives Multiwfn (`2`,`2`,`3`,`4`,`5`,`8`,`9`,`7`→`0`, `-5`→`1`) and parses `CPprop.txt` and the path summary into `tests/data/qtaim/<name>_multiwfn.txt`. It is not run by the suite.

  Fixtures:
  - water (6D), ethylene, benzene (RCP) and CH₃ (open shell), which exist already;
  - new G09W fixtures: water dimer (H-bond BCP) and bicyclo[1.1.1]pentane (CCP), with inputs in `tests/data/qtaim/inputs/`.

  Agreement:
  - same CP counts and assignments;
  - positions ≤ 1e-5 bohr;
  - ρ, ∇²ρ, G, K, V, H and eigenvalues ≤ 1e-6 relative (or 1e-10 absolute);
  - path lengths ≤ 1e-3 bohr.
- **Report:** layout, every section, NOTE lines.
- **UI:** headless via `QTAIMTabController.__new__` with stub widgets. Every `run()`/`plane_map()` argument comes from a control; no-flicker rules.
- **Plot pipeline:** `draw_plane_map` renders under Agg; `plotter` validates and describes a `plane_map`; the subprocess tests for `toga`-free imports.

## 5. Not in v1

- 3D interbasin surfaces;
- AIM basin integration (charges, volumes);
- mode 1 (user points) and mode 6 (random spheres);
- Shannon aromaticity and ring-normal curvature;
- ESP and other real-space functions at CPs;
- correlated or spin densities;
- `.wfn`/`.wfx`;
- periodic systems;
- EDF core densities.
