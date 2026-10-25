# Molecular Orbitals tab (3D Molecular Viewer) — design

Date: 2026-10-25. Status: approved by the user.

## 0. User decisions

| Question | Decision |
|---|---|
| Calculation types | SCF only: restricted (RHF/RKS), unrestricted (UHF/UKS, alpha and beta lists), restricted-open (ROHF/ROKS, one set with occupancies 2/1/0). Natural orbitals are out of scope. |
| Orbital selection | A table of every MO (No., label, occupancy, energy Ha/eV), one MO drawn at a time; HOMO / LUMO / Previous / Next buttons; spin Selection for U files. No in-tab diagram. |
| Display options | Isovalue (0.02), +/- colours (GaussView red/green), opacity; plus grid presets Coarse/Medium/Fine/Custom, Solid/Transparent/Mesh style, orbital density \|phi\|^2, energy-level diagram in the plot viewer. |
| Exports | Cube of the displayed MO; cubes of an MO range; Send to Cube tab; 3D image; MO energy table (.txt). |
| Tab | Label "Orbitals", placed after Cube: Display, Box, Dynamics, Measure, Design, Cube, Orbitals, ESP, QTAIM, NTO, Excitation. |

## 1. Approach

New modules, existing tabs untouched:

- `molecularOrbitals.py` — numerical, no Toga (subprocess test): catalogue, orbital specs, grids, cubes, energy table, level-diagram manifest.
- `orbitalTab.py` — `OrbitalTabController`, the only module of the two that imports Toga; imports `molecularViewer` only inside `export_3d_image`.

Chain: `molecularViewer → orbitalTab → {cubeTab, cubeRender, qtaimTab, molecularOrbitals → {gaussianWavefunction, molecularSurface, cubeFile, indexRanges, ntoAnalysis (grid helpers), help}}`.

Rejected: extracting a base class from `ntoTab` (rewrites a validated tab), loading MOs in the Cube tab (no orbital list, energies or HOMO/LUMO).

## 2. Parsing and catalogue

- `gaussianWavefunction.read_fchk` unchanged. A file without `Alpha MO coefficients` or `Alpha Orbital Energies` is a `ValueError` naming the missing section.
- `orbital_catalogue(wfn) -> OrbitalCatalogue` with `kind` in `RESTRICTED`, `UNRESTRICTED`, `RESTRICTED_OPEN`:
  - restricted (no beta set, n_alpha == n_beta): occupancy 2 for MO <= n_alpha, else 0;
  - restricted-open (no beta set, n_alpha != n_beta): 2 for MO <= n_beta, 1 up to n_alpha, else 0;
  - unrestricted: alpha 1 for MO <= n_alpha, beta 1 for MO <= n_beta.
- Each `OrbitalInfo`: 1-based `number`, `spin` ("alpha"/"beta"), `energy` (Ha), `occupancy`, `label` relative to that spin's HOMO/LUMO (`HOMO`, `HOMO-1`, `LUMO`, `LUMO+2`). eV = Ha × 27.211386245988 (CODATA 2018, as the Fukui/Excitation tools).
- HOMO = highest MO with occupancy > 0 in that spin list; LUMO = HOMO + 1 (absent if HOMO = n_mo).
- `parse_orbital_spec(text, catalogue, spin)` → sorted unique 1-based MO numbers. Accepts numbers and number ranges (`20-25, 30`, through `indexRanges.parse_solute_index_ranges`) and HOMO/LUMO tokens with offsets (`HOMO-2`, `LUMO+1`), with `:` as the range operator between tokens (`HOMO-2:LUMO+2`, `10:LUMO`). Out-of-range numbers are a `ValueError`.

## 3. Orbital evaluation

- `orbital_volumes(wfn, spin, numbers, spacing, isovalue, progress=None, use_numba=None) -> List[CubeVolume]`: all requested MOs on one grid. Box = `molecularSurface.grid_box` (atoms ± 1.7 r_vdW, each face grown until max |phi| over the MOs < isovalue/10); `ntoAnalysis.check_grid_size` (25 M points) before evaluation; slabs of 2 M points through `gaussianWavefunction.orbital_values` (one AO pass for all MOs; numba AO kernel when installed). No new numba kernel.
- Grid presets (cubegen's): Coarse 3, Medium 6 (default), Fine 12 points per bohr; Custom = a typed spacing in bohr.
- Density mode: |phi|^2 = square of the cached amplitude grid (no new evaluation), one surface, its own isovalue (default 0.0004 = 0.02^2).
- Signs are the fchk coefficients as written, so phases match GaussView/cubegen for the same file.

## 4. Visualisation

- `viewer.orbital_layer` (a `cubeRender.CubeSceneLayer`), added to `MolecularViewerUI.ANALYSIS_LAYERS`, `render_scene_extras`, `get_scene_radius` (getattr-guarded). The tab claims the canvas with `show_analysis_layer("orbital_layer")` only for a non-empty result.
- Each displayed MO is a `cubeTab.CubeEntry` (meshes cached per isovalue and sign). Colour/opacity/style changes never re-mesh; MO, spacing or isovalue changes do. The last 3 amplitude grids are kept (LRU) so Previous/Next is instant.
- Styles: Solid (opacity 1), Transparent (the opacity slider, default 0.5), Mesh. Mesh = `volumeData.contour_lines` (the surface cut by planes 0.4 bohr apart in x, y and z, GaussView's look) drawn as an unlit `cubeRender.line_item` (new `"lines"` kind, `GL_LINES`); `surface_item` is unchanged. (Revised during implementation: the approved `wireframe` flag drew every marching-tetrahedra triangle with `GL_LINE`, which in the live window looked like a solid surface with speckles.)
- Colours default to red (+) / green (−); colour names or `#rrggbb` (`cubeTab.canonical_colour`).
- Structure: the NTO tab's `_settle_structure` rule (silent load when the viewer has no frames or shows the same molecule; asked once otherwise).

## 5. Tab layout and behaviour

Input (fchk + Browse, output folder + Browse, Load) → summary line (method/basis, electrons, R/U/RO, HOMO, LUMO, gap) → Orbitals (spin Selection, enabled for U; table No./Label/Occ./E (Ha)/E (eV), scrolled to the HOMO; HOMO, LUMO, Previous, Next, Energy levels...) → Display (surface: phases or |phi|^2; isovalue; grid quality + custom spacing; style; opacity; +/- colours; Display; status; progress) → Export (cube of displayed MO; MO range + Export cubes; Send to Cube tab; Export 3D image...; Export energy table) → help text (`HelpGqteaWin.help_orbitals`).

Load parses in `asyncio.to_thread`. After the first Display, a row click, Previous/Next, HOMO/LUMO or spin change redraws automatically. Conventions of the NTO tab: widget writes through `cubeTab._assign`, synchronous `_syncing`-guarded handlers, `to_thread` + `call_soon_threadsafe` progress, coalesced `refresh_scene`, every error a dialog.

## 6. Outputs

- `<stem>_MO<n>[_alpha|_beta][_density].cube` through `cubeFile.write_cube_volume` (ASCII titles).
- Range export: one cube per MO on the shared grid (Cube-tab arithmetic compatible).
- `<stem>_MO_energies.txt`: commented header (file, method, basis, kind, electrons, HOMO/LUMO/gap per spin), then `mo spin occupancy label energy_hartree energy_ev`.
- `<stem>_MO_levels.json`: a plotViewer manifest with the new `{"type": "energy_levels"}` shape: `columns` = [{`label`, `levels`: [{`number`, `energy_ev`, `occupancy`, `label`}], `homo`, `lumo`}], `ylim`, `title`. Rendered by `displayPlots.draw_energy_levels` (the single renderer for plotViewer and plotter's PNG fallback; imported in plotViewer above the TkAgg pin); occupied solid, virtual dashed, degenerate levels side by side, gap annotated, default window HOMO−10 … LUMO+10 but at most 20 eV beyond the HOMO/LUMO (added after a CH3 check: its C 1s level at −277 eV squashed the frontier region). plotter validates it (`_check_energy_levels`) and describes it.

## 7. Testing

- `tests/test_molecularOrbitals.py`: R/U/RO catalogues; spec parser matrix; Σ occ |phi|^2 = rho on a grid; agreement with cubegen MO cubes (dev-time `tests/data/orbitals/make_reference.py`, gzip-stored, near-nucleus points excluded); shared range grid; cube round trip; energy table format; level manifest; no Toga import.
- `tests/test_orbitalTab.py`: headless (`__new__` + stub widgets): every `orbital_volumes` parameter supplied by `_collect_display_inputs`; selection/redraw rules; no-flicker writes; layer claim; exports.
- Plot shape: plotter validation, Agg render; tab order and `ANALYSIS_LAYERS` pin tests updated.
- New fixture: `tests/data/orbitals/ch3_rob3lyp_631gd.fchk` (ROB3LYP/6-31G(d), G09W), with its input.

## 8. Findings during implementation

- G09W's `formchk` writes the **alpha density** as `Total SCF Density` for a ROHF job (the new CH3 fixture integrates to 5 of 9 electrons). The Orbitals tab does not use it; density-based tabs (ESP, QTAIM, Fukui) would see only the alpha electrons of such a file. Pinned by a test, documented in the help text and the manual; not fixed in those tabs (out of scope).
- Agreement with cubegen is limited by cubegen's six printed digits and its near-nucleus approximation: ≤ 3e-5 a.u. at points more than 0.3 bohr from a nucleus, identical signs.

## 9. Not implemented

Natural orbitals, `.wfn`/`.molden`/`.log` input, several MOs at once, an in-tab diagram, MO composition, spin density.
