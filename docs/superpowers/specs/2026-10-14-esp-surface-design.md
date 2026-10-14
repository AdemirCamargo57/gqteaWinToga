# ESP mapped on an electron-density isosurface — design & implementation plan (DRAFT)

Status: **draft for review, not committed.** Open questions are in §10;
anything marked *[assumption]* or *[proposed]* has not been verified.

## 0. Decisions (user, 2026-10-14)

| # | Decision |
|---|---|
| D1 | Reference input: `C:\AJCamargo\QCMM_2sem_2026\2026_10_08\Eteno.fch` (ethylene, RB3LYP/6-31G, opt; 6 atoms, 26 basis functions, s/p/SP shells only; has `Total SCF Density` and a (zero, centrosymmetric) `Dipole Moment`). Visual reference: the GaussView 6.1 screenshot `Eteno_ESP.png` in the same folder (**not** the TIF). |
| D2 | Visualisation: subclass of the existing GLFW/OpenGL `MolecularViewer` (§4.4). |
| D3 | Extrema: copy Multiwfn exactly — strict 1st+2nd-ring rule, no automatic merge, vertex positions (§2, §6). |
| D4 | Colour bar: GaussView defaults (§3). |
| D5 | Density: `Total SCF Density` only. No density selector in the UI. |
| D6 | numba accepted as an **optional** speed-up: used when importable, pure numpy otherwise; listed as optional in `installation_requirements.txt`. |
| D7 | Gaussian 09W `cubegen` available at `C:\G09W\cubegen.exe` — used as test reference and optional ESP engine (§4.2). |

## 1. Starting point (verified)

- [electrostaticSurfacePotential.py](../../../electrostaticSurfacePotential.py) is a stub:
  `class ElectrostaticSurfacePotentialUI: pass`. It is already imported by
  [gqteaWinToga.py](../../../gqteaWinToga.py) (line 56) and registered on the **Tools** tab
  as "Electrostatic Surface Potential (ESP)". Clicking it today constructs an empty object.
- No `.fchk` reader, Gaussian-basis code, integral code, or isosurface code exists anywhere
  in the tree. Nothing in `help.py`, `readme.md` or `USER_MANUAL.md` mentions ESP.
- Installed in the venv: numpy 2.4.2, scipy 1.17.0, matplotlib 3.10.8, glfw 2.10.0,
  PyOpenGL 3.1.10, toga 0.5.3. **Not** installed: scikit-image (marching cubes), numba,
  PySCF/libcint (no Windows wheels anyway).
- [molecularViewer.py](../../../molecularViewer.py) renders with legacy fixed-function GL
  (display lists, `glBegin`) in a GLFW daemon thread. `main_loop` has **no hook** for extra
  geometry: it calls `render_frame` inside the molecule push-matrix, then
  `draw_cartesian_axes_overlay`. Coordinates are Å.
- `help.py` `AtomicData` has masses, covalent radii, valence data — **no van der Waals radii**.
- Bohr constant precedent: `coordinateConverter.BOHR_RADIUS_ANGSTROM = 0.529177210544`
  (CODATA 2022) — reuse it rather than adding another bohr constant.

## 2. Multiwfn methodology (verified from source)

Source read: `src/surfana.f90`, `src/fileIO.f90`, `src/define.f90` from the
`ruixingw/Multiwfn-mirror` GitHub mirror (exact Multiwfn version of the mirror not confirmed).
Method paper: T. Lu, F. Chen, *J. Mol. Graph. Model.* **2012**, 38, 314 (doi:10.1016/j.jmgm.2012.07.004).

| Item | Multiwfn behaviour | Where |
|---|---|---|
| Default surface | ρ isosurface, `surfisoval = 0.001` a.u. | surfana.f90:25 |
| Grid spacing for ESP | `grdspc = 0.25` bohr (0.2 for ALIE/LEA etc.) | surfana.f90:29, 213 |
| Grid box | atoms ± `vdwmulti`·r_vdW, `vdwmulti = 1.7`; manual notes 1.7 suffices for 0.001 and **lower isovalues need more** | :32, :308, :384-389 |
| Polygonisation | Marching **tetrahedra** (main-axis decomposition of each cube) | :621, `marchtetra` |
| Vertex placement | `nbisec = 3` bisections of the cut edge using the **exact** density, then linear interpolation | :2997-3031 |
| Vertex merge | connected vertices closer than `critmerge = 0.5 × grdspc` are merged to their midpoint | :30-31, :659-725 |
| Vertex elimination | after merging, vertices with exactly 2 neighbours, then with exactly 3, are removed (3-neighbour hole re-capped with one triangle) | :727-780 |
| Topology check | warns if V + F − E − 2 ≠ 0 (closed single surface) | :637-642 |
| ESP evaluation | exact ESP at every surviving vertex (optionally via Gaussian `cubegen`) | :855-990 |
| Global extrema | vertex with the lowest / highest value | :1090-1105 |
| **Local minimum** | value **strictly lower** than every 1st-ring neighbour **and** every 2nd-ring neighbour (neighbours of neighbours, excluding itself) | :1107-1129 |
| Local maximum | mirror image (strictly higher) | :1146-1166 |
| Distinct-extrema merging | **none automatic.** Only manual menu options: discard by index (18), merge selected extrema to their mean position (19). Positions are **vertex positions, not refined.** | :1373-1374, :2142-2198 |
| Reporting | `a.u.`, `eV`, `kcal/mol`, `X/Y/Z (Angstrom)`; global extremum marked `*` inside the local list | :1110-1125, `surfanalysis.txt` |
| Conversion constants | `au2kcal = 627.51`, `au2eV = 27.2113838`, `b2a = 0.529177210903` (CODATA 2018) | define.f90:28-29 |
| ECP | reads `Nuclear charges` from the fch (≠ atomic number with ECP) and uses them as nuclear charges | fileIO.f90:506 |
| fch basis | each Cartesian primitive gets its **own** normalisation `normgau(type, exp)` × the fch contraction coefficient; `SP` shells use `P(S=P) Contraction coefficients` for the p part; pure shells (5D/7F/9G/11H) are converted to Cartesian | fileIO.f90:519-560, 686-745 |
| fch Cartesian order | s,p,d same as Gaussian `pop=full` (d: XX,YY,ZZ,XY,XZ,YZ); **f order in fch differs** from wfn (fch: XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ); g order (ZZZZ … XXXX) | define.f90:220-240 comments |

Consequence for "numerical criteria for distinct local extrema": in Multiwfn two extrema of
the same type are distinct **by construction** whenever both pass the strict 2-ring test,
so they are always ≥ 3 edges apart (≈ 3 × 0.25 bohr ≈ 0.4 Å on a 0.25-bohr grid after
merging). There is no value tolerance and no distance merge.

## 3. GaussView 6.1 conventions (verified from `Eteno_ESP.png`, a GaussView 6.1 screenshot)

Reference image: `C:\AJCamargo\QCMM_2sem_2026\2026_10_08\Eteno_ESP.png` (use this, not the TIF).

- **Window title:** `Electron density from Total SCF Density ( isoval = 0.0004 ; [mapped with ESP] )`
  — the surface quantity (ρ isovalue) and the mapped quantity (ESP) are named separately.
- **Colour bar:** horizontal strip across the **top** of the view, min label at the left end,
  max label at the right end, no unit printed (GaussView values are a.u.).
- **Range:** **symmetric**, `-2.072e-2 … 2.072e-2` for Eteno at ρ = 0.0004. Working rule:
  R = max(|V_min|, |V_max|) over the displayed surface; checked once our Eteno V_min/V_max
  exist (GaussView interpolates ESP from a cube grid, so agreement is approximate, ~1e-4).
- **Label format:** `%.3e`-style with the exponent's leading zero dropped (`-2.072e-2`).
- **Colour order:** red (−R) → orange → yellow → green (0) → cyan → blue (+R): a linear
  **hue sweep 0° → 240°**. On Eteno, red sits over the C=C π cloud, cyan/blue toward the H atoms.
- **Rendering:** opaque surface, Phong-style specular highlight, **black background**.
- Our bar additionally prints the unit (`ESP (a.u.)` or `ESP (kcal/mol)`) and the surface
  definition (`on ρ = 0.000400 a.u.`), which the requirements ask for and GaussView omits.
- The default isovalue 0.0004 a.u. is the user's requirement (Multiwfn's default is 0.001).

## 4. Proposed architecture

One-way dependency chain, Toga only at the top, everything below headless-testable
(same pattern as the molecular-design modules):

```
electrostaticSurfacePotential.py  (UI + analyser orchestration + ESPSurfaceViewer)
   ├── molecularViewer.MolecularViewer   (subclassed for the 3D window)
   ├── molecularSurface.py   isosurface (marching tetrahedra) + extrema + export rows
   ├── espIntegrals.py       Boys function + McMurchie–Davidson ESP integrals
   └── gaussianWavefunction.py   .fchk parser, basis expansion, ρ(r) on points
          └── help.py (AtomicData: symbols, + new vdw_radii)
```

### 4.1 `gaussianWavefunction.py` — FCHK parsing and density *[proposed]*
- `read_fchk(path) -> FchkData`: generic section reader (`I`/`R`/`C` scalars and `N=` arrays),
  then extracts: title/route, charge, multiplicity, electron counts, `Atomic numbers`,
  `Nuclear charges`, `Current cartesian coordinates` (bohr), `Shell types`,
  `Number of primitives per shell`, `Shell to atom map`, `Primitive exponents`,
  `Contraction coefficients`, `P(S=P) Contraction coefficients`, `Total SCF Density`
  (lower-triangular; D5), and `Dipole Moment` when present (validation only).
  Accepts `.fchk` and `.fch`. Every failure is a `ValueError` naming the missing/short section.
- `BasisSet`: shells → Cartesian primitives with per-component normalisation (Multiwfn
  convention above), SP split, pure→Cartesian transformation `T` per shell so that
  `P_cart = T P_pure Tᵀ` (one code path for 5D/6D/7F/10F/…). Angular momentum up to g
  fully; h if Q-time allows *[proposed]*.
- `density_at(points)`: chunked AO evaluation with per-shell distance screening, then
  `ρ = Σ (χ P) ∘ χ`. Pure numpy/BLAS.

### 4.2 `espIntegrals.py` — exact ESP *[proposed]*
`V(r) = Σ_A Z_A/|r−R_A| − Σ_μν P_μν ∫ χ_μ χ_ν /|r−r'| dr'` (Z_A = fchk `Nuclear charges`).
- McMurchie–Davidson: each primitive pair → Hermite expansion `E^{ij}_t`; **contract the
  density matrix into Hermite coefficients D_tuv once per primitive pair** (point-independent),
  so per point the work is `2π/p Σ D_tuv R_tuv(p, P−C)`.
- Boys function via `scipy.special.gammainc`/`gamma` with downward recursion, series for
  small T, asymptotic form for large T.
- Vectorised over (batch of primitive pairs × chunk of points); pairs screened by
  `|D| · exp(−μ R_AB²) < 1e-12`; point chunks run in a `ThreadPoolExecutor` (numpy releases
  the GIL) to use several cores.
- **numba (D6):** the two hot kernels (Hermite `R_tuv` evaluation and AO values) get an
  `@njit(parallel=True)` twin selected at import time (`try: import numba`); both paths are
  run by the tests when numba is present and must agree to 1e-12.
- **cubegen engine (D7) — verified working:** `cubegen 0 potential=SCF Eteno.fch out.cub -5 h
  < points` with `GAUSS_EXEDIR=C:\G09W` reads free-format x y z points **in Angstrom** from stdin (corrected by the spike: a point typed 0.1173 lands on the water O at 0.1173 A)
  and writes `x y z V` lines after the atom block; 20 000 points on Eteno took 0.5 s
  (e.g. V at (4,0,0) = −0.004575585872 a.u., above the π plane). This is the route Multiwfn
  uses (surfana.f90:880-905). UI: "ESP engine: Built-in exact / Gaussian cubegen" + path
  field, auto-detected from `GAUSS_EXEDIR`, `PATH` and `C:\G09W`; on any cubegen failure
  (e.g. a G16 fchk that G09's cubegen cannot read) the run falls back to the built-in engine
  and says so in the summary and the file header (`# esp_method`). The same command with
  `density=SCF` gives ρ reference values.

### 4.3 `molecularSurface.py` — isosurface and extrema *[proposed]*
- Grid: uniform cubic, spacing `h` (default 0.25 bohr, Multiwfn's ESP default). Box = atoms ±
  margin, starting from 1.7·r_vdW and **grown until max ρ on every box face < isovalue/10**
  (needed for 0.0004). One grid serves all requested isovalues (box sized for the smallest).
- Marching tetrahedra (6 tetrahedra per cube on the main diagonal, 16-case table — small
  enough to vectorise in numpy, no new dependency); vertices de-duplicated by global edge id
  so the mesh is connected; 3 bisections with exact ρ then linear interpolation; Multiwfn's
  0.5·h merge and 2-/3-neighbour elimination.
- Mesh outputs: vertices (bohr), triangles, vertex adjacency (CSR), vertex normals (from ∇ρ,
  analytic or from the grid), area, enclosed volume, Euler characteristic per connected
  component (warn when a component is not closed or χ ≠ 2).
- `find_surface_extrema(values, adjacency)`: Multiwfn's strict 2-ring rule, vectorised over
  the CSR adjacency; global min/max = argmin/argmax over all vertices, flagged even if a tie
  keeps it out of the strict local list.
- Optional, **off by default** [Q3]: merge same-type extrema closer than `d_merge` Å.

### 4.4 `electrostaticSurfacePotential.py` *[proposed]*
- `ESPSurfaceAnalyser` (headless): `run(fchk_path, isovalues, spacing, density_key,
  progress=callable) -> ESPResults` — parse → density grid → per isovalue: mesh → ESP at
  vertices → extrema. Same `progress` callback + `loop.call_soon_threadsafe` marshalling
  used by `meanResidenceTime_statistical.py`; runs in `asyncio.to_thread`.
- `ElectrostaticSurfacePotentialUI(*args, **kwargs)`: `layout_main_window` / `async workflow`
  pattern; newer `margin=` style; `toga.InfoDialog`/`window.dialog(...)` errors, never raises.
  Two visually separate groups so surface definition and mapped values cannot be confused:
  - **Surface definition — electron density ρ (a.u. = e/bohr³):** isovalue(s) field
    (`0.0004, 0.001, 0.002`; blank → 0.000400), grid spacing (blank → 0.25 bohr);
    density is always `Total SCF Density` (D5), shown read-only.
  - **Mapped property — ESP V(r) (a.u. = Hartree/e):** colour-range mode, colour-bar units
    (a.u. / kcal/mol), optional manual min/max.
  - Run / Show surface / Clear / Close; the extrema file is written automatically after every Run (no separate Export button); output folder (blank → next to the fchk);
    `MultilineTextInput` summary; help text `HelpGqteaWin.help_esp_surface`.
- `ESPSurfaceViewer(MolecularViewer)`: loads the molecule as a single-frame trajectory (same
  trick as `load_designed_molecule`), draws the active surface with vertex arrays
  (`glVertexPointer`/`glColorPointer`/`glNormalPointer` + `glDrawElements`), optional
  transparency, extrema as small spheres (max/min colours) with labels, and a 2D colour bar
  (GLUT bitmap text, already used by the viewer) titled e.g.
  `ESP (a.u.)  on  ρ = 0.000400 a.u.` with min / 0 / max tick labels. A Selection switches
  between isovalues.

### 4.5 Touches to existing files
- `molecularViewer.py`: two **no-op hooks** on the base class — `render_scene_extras()`
  (inside the molecule push-matrix, after `render_frame`) and `draw_screen_overlay()`
  (after the axes overlay) — plus letting `get_active_scene_radius` include extra extent.
  Same no-op-hook pattern as `_notify_picked_atoms_changed`; no behaviour change for the viewer.
- `help.py`: `AtomicData.vdw_radii` (Bondi 1964 + Mantina 2009 for missing main-group
  elements), `HelpGqteaWin.help_esp_surface`.
- `gqteaWinToga.py`: none (already registered).
- `USER_MANUAL.md`, `readme.md`, `CLAUDE.md`: new tool section / capabilities / module notes.
- No new runtime dependency → no packaging change.

## 5. Extrema output file *[proposed]*

`<fchk stem>_ESP_extrema.txt` next to the input (or chosen folder), one file for all isovalues:

```
# gQTEA ESP Surface Extrema
# source_file            water.fchk
# density                Total SCF Density
# surface_definition     electron density isosurface, rho in a.u. (e/bohr^3)
# mapped_property        electrostatic potential V(r), a.u. (Hartree/e)
# esp_method             exact (McMurchie-Davidson) | cubegen
# grid_spacing_bohr      0.25
# vertex_merge_bohr      0.125
# extremum_criterion     strict 1st+2nd-ring vertex neighbours (Multiwfn convention)
# hartree_to_kcal_mol    627.5094740631
# coordinates_unit       Angstrom
# surface rho=0.000400  vertices=… area_A2=… V_min_au=… V_max_au=…
# row rho_isovalue_au type scope esp_au esp_kcal_mol x_A y_A z_A nearest_atom
1 0.000400 minimum global -0.061234 -38.4253 0.000000 0.000000 -1.912345 O1
2 0.000400 maximum global  0.054321  34.0868 1.401234 0.000000  0.812345 H2
…
```

`type` ∈ {minimum, maximum}; `scope` ∈ {global, local}; every global extremum is listed
once with `scope=global`. Column names never use "isovalue" for an ESP quantity.

## 6. Numerical criteria to document (proposed defaults)

1. Grid spacing h = 0.25 bohr (Multiwfn ESP default); box grown until face ρ_max < iso/10.
2. Vertex position: 3 exact-ρ bisections + linear interpolation (residual |ρ−iso| reported).
3. Vertex merge: edge-connected vertices closer than 0.5 h merged; 2-/3-neighbour vertices removed.
4. Extremum: strict inequality against all 1st- and 2nd-ring neighbours (ties → no extremum).
5. Global: arg-min / arg-max over all vertices of that surface.
6. Units: 1 Hartree/e = 627.5094740631 kcal/(mol·e) (CODATA 2018 E_h·N_A; Multiwfn uses
   627.51, a 8×10⁻⁷ relative difference); coordinates in Å via CODATA 2022 a₀.

## 7. Validation plan

| Layer | Test |
|---|---|
| Parser | section reader on synthetic fchk text; malformed-file rejection matrix; 5D vs 6D of the same job give the same ρ |
| Basis / ρ | ∫ρ dV = N_electrons (< 1e-3 relative); dipole from ρ + nuclei vs fchk `Dipole Moment` (< 1e-3 a.u.) |
| Integrals | s–s vs closed form; d/f pairs vs brute-force numerical quadrature; far field → Q/r + μ·r/r³; translation/rotation invariance |
| ESP end-to-end | vs `cubegen` at the same points (target < 1e-6 a.u.). Reference values are generated **once** with cubegen and committed under `tests/data/esp_reference/` (points + ρ + V files, plus the generator script), so the suite never needs Gaussian at test time |
| Isosurface | analytic Gaussian "atom" → sphere radius error < h/10, area vs 4πR² < 1 %, χ = 2; two separated spheres → two components |
| Extrema | synthetic fields on a sphere mesh with known extrema; plateau → none; global flagged |
| Export | header, units, a.u.→kcal/mol, Å coordinates, one global per type per surface |
| UI | headless via `__new__` + stub widgets (project pattern); viewer colour-map / colour-bar tick maths as pure functions |
| Performance | timed benchmark on a medium molecule (see risk R1) |

## 8. Implementation order

0. **Spike (throwaway):** time the numpy ESP kernel on ~20 k points for a ~20-atom 6-31G(d)
   molecule to confirm R1 before committing to pure numpy.
1. fchk parser + basis + ρ (tests: ∫ρ, dipole).
2. ESP integrals (tests: analytic, quadrature, far field, reference values).
3. Isosurface + vertex elimination (tests: sphere, topology).
4. Extrema + export (tests: synthetic fields, file format).
5. Analyser orchestration + progress.
6. Viewer hooks + `ESPSurfaceViewer` (manual check against GaussView screenshot).
7. Toga UI, help text, manual/readme/CLAUDE.md.

## 9. Risks and assumptions

- **R1 performance [assumption]:** pure-numpy ESP is O(primitive pairs × vertices). Estimate
  for a 20-atom 6-31G(d) molecule at 0.0004 a.u. / 0.25 bohr: ~20–40 k vertices, tens of
  seconds to a few minutes. Unmeasured — hence the spike. Fallbacks: cubegen backend,
  coarser spacing, or numba (new dependency).
- **R2 box size at 0.0004:** the surface sits far out (≈ 2.5–3 Å beyond H); diffuse
  (`+`/`aug-`) bases push it further. Handled by the adaptive box; costs grid memory.
- **R3 ECP:** with ECPs the density lacks core electrons and Z is the effective charge
  (Multiwfn does the same); fine outside the core, which is where the surface is. Documented.
- **R4 f/g ordering:** taken from Multiwfn's documented fch convention; only proven by the
  5D-vs-6D and ∫ρ/dipole tests on real Gaussian files.
- **R5 multiple surfaces:** at some isovalues a weakly bound complex gives several closed
  components; extrema are still per surface, the Euler check is reported per component.
- **R6 GaussView range rule:** the range is verified symmetric (±2.072e-2, §3); that R equals
  max(|V_min|, |V_max|) is still an inference, to be checked against our Eteno result.
- **R7 coverage of Eteno.fch:** s/p/SP only — covered by the extra G09W fixtures (§10).

## 10. Resolved questions

1. GaussView range → §3 (from `Eteno_ESP.png`).
2. Extra fixtures → **approved**: small G09W jobs from `C:\G09W` — water HF/6-31G(d) with
   `6D` and with `5D`, water/cc-pVTZ (pure f), formaldehyde 6-311G(2df) with `10F`
   (Cartesian f ordering), methyl radical UB3LYP/6-31G(d) (open shell) — committed with
   their cubegen ρ/ESP reference values under `tests/data/`.

Spec approved by the user on 2026-10-14.

## 11. Spike results (2026-10-14, throwaway code in the session scratchpad)

Measured, not estimated — these supersede R1 and the assumptions above:

- **ESP vs Gaussian 09 cubegen**, 300 points 0.5-8 bohr from the nuclei: max |dV| <= 3.7e-7 a.u.
  on all seven fixtures (5D, 6D, pure f, pure g, Cartesian 10F, SP shells, open shell).
- **tr(P S) = N_electrons** to 1e-8 and every fchk basis function normalised to 1e-8.
- **Density vs cubegen**: <= 2.5e-4 relative for 1e-5 < rho < 1e-2. cubegen's density is
  screened near nuclei (12.8 vs 290 a.u. at the water O nucleus); irrelevant for surfaces.
- **Eteno, rho = 0.0004**: 6655 vertices, V+F-E = 2, |rho - iso| <= 1.2e-6; 2 pi minima
  -0.020967 a.u. (-13.16 kcal/mol), 4 H maxima +0.01166 a.u.; max|V| = 2.097e-2 vs GaussView's
  bar 2.072e-2 (GaussView interpolates from a coarse cube) — R = max|V| confirmed (R6). ~1 s total.
- **Timing, benzene B3LYP/6-31G(d)** (102 bf), 20 000 points: numpy 19.6 s single thread,
  4.0-4.7 s with 12 threads, **0.9 s with the numba kernel** (agrees with numpy to 1e-11).
- numba 0.68.0 supports the venv's numpy 2.4.2 (no downgrade).
- The GL viewer renders the GaussView look (screenshot checked against `Eteno_ESP.png`).
