# Solvation box: periodic pre-optimization, system properties, construction canvas, acceleration

Status: approved 2026-10-15. Branch `feature/solvation-box-preopt`.

## 0. User decisions

| Question | Decision |
|---|---|
| Pre-optimization model | **Both, user selectable**: rigid bodies (default) and flexible (input-geometry-referenced bonded terms). |
| Charge input | **Per species**: solute charge + charge per solvent molecule; total = q_solute + N_solvent * q_solvent. |
| Spin multiplicity | **Always required**; never inferred. Parity and range validated once N_e is known. |
| Canvas placement | **Tabs** Setup / Construction / Summary in the existing window. |
| Phasing | single_solute_solvent_box.py first, validated; then mixtureSolventBox.py, then sharedWallBox.py; consolidation assessment last and separate. |

## 1. Baseline (measured before any change)

20 A cubic water box, target 266 molecules, 100 000 attempts, seed 7: 234 accepted, 115 s under cProfile.
About 90 % of the time is interpreter overhead, not arithmetic: every attempt scans every accepted molecule;
`round()` on numpy scalars inside `minimum_image_displacement` (39 M calls, 42 s) and per-pair vdW radius look-ups
(15 M calls) dominate. The first win is algorithmic (cell list + vectorisation); Numba is a second, measured layer.

## 2. Modules

One-way dependencies; only the `*UI` module and the canvas wrapper import Toga.

- `solventPacking.py` -- `ClashChecker`: per-atom vdW radii as arrays, accepted molecules in a cell list keyed
  by centre, candidate-vs-neighbour atom clash test vectorised with the minimum image. **Decisions are
  bit-identical to the legacy loops**: same centre pre-filter `d_centre > r_cand + r_exist + 2.5`, same float64
  operation order for `dx -= L*round(dx/L)` (numpy and Python both round half to even), same `cutoff = s*(ri+rj)`
  and `dist2 < cutoff^2` / `dist2 < (cutoff+extra)^2`. So a fixed seed reproduces the legacy box exactly.
- `solventPackingNumba.py` -- optional compiled clash kernel; the only importer of `numba`; numpy fallback.
  Adopted as default only if it beats the numpy path in the benchmark.
- `periodicPreOptimizer.py` -- `PeriodicCell` (lengths, periodic flag, `minimum_image`, `wrap_molecules`),
  `PackedSystem` (elements, coordinates, molecule id per atom, molecule templates), the intermolecular penalty,
  `RigidBodyOptimizer` and `FlexibleOptimizer`, and `PreOptimizationResult`.
- `systemProperties.py` -- electrons, allowed multiplicities, validation messages.
- `boxConstructionView.py` -- pure projection/draw-list (`build_draw_list`) and `ConstructionCanvas`
  (Toga wrapper with throttled redraws).
- `molecularPreOptimizer.py` -- adds `ForceField.from_geometry(elements, coordinates, bonds)` and caches term
  arrays in `__init__` (no numerical change for sketches).

## 3. Pre-optimization physics

Intermolecular penalty (both modes), Packmol-style:
`E = K * sum_{pairs in different molecules, d < sigma} (d^2 - sigma^2)^2`, `sigma = vdw_scale * (r_vdW,i + r_vdW,j)`
-- the builder's own acceptance criterion, so E = 0 exactly means "no clash left". Finite range, so PBC is clean:
neighbour pairs from `cKDTree(wrapped, boxsize=L).query_pairs(max sigma)`; distances recomputed with the minimum
image on the unwrapped coordinates; same-molecule pairs excluded by molecule id. Requires `max sigma < min(L)/2`.
The solute-solvent extra clearance is added to sigma for solute-solvent pairs.

Non-periodic (minimum-image switch off): no images, plus a harmonic wall on atoms leaving
`[pad, L - pad]` (a repulsion-only potential would otherwise eject molecules).

**Rigid**: 6 DOF per molecule, `x_i = R(w) b_i + t` with body coordinates `b_i` relative to the molecule's
current centroid. Gradient: `dE/dt = sum g_i`; `dE/dw_k = sum g_i . (dR/dw_k b_i)` with the exact rotation-vector
derivative (Gallego & Yezzi 2015), `[e_k]x` at w = 0. Internal geometry preserved to round-off.

**Flexible**: per molecule template, bonds detected once on the input geometry (`r <= r_cov,i + r_cov,j + 0.45`,
the viewer's rule); `ForceField.from_geometry` takes reference lengths/angles from that geometry, adds the
existing >= 4-bonds-apart intramolecular repulsion, no torsions (no bond orders). Terms are tiled over molecule
copies with index offsets. Intramolecular vectors never use the minimum image (molecules are kept whole).

Minimiser: `scipy.optimize.minimize(L-BFGS-B)`, `max_iterations` from the UI (default 500), iteration callback
for progress/canvas snapshots. Non-convergence is `success=False`, never an exception.

Molecule-aware wrapping: each molecule shifted as a whole by `-L * floor(COM / L)` (mass-weighted COM). No
molecule is split; atoms may protrude past a face.

## 4. System properties

`N_atoms = len(data)`; `N_e = sum Z - Q`, `Q = q_solute + N_solvent * q_solvent` (integer charges).
Multiplicity M: required integer >= 1 before the build; after the build valid iff `(N_e + M) % 2 == 1` and
`M <= N_e + 1` and `N_e >= 0`. Inconsistent -> dialog with allowed values, summary marks it, and
"Update properties" re-validates without rebuilding. Unknown element (no atomic number) -> error naming it.

## 5. UI and outputs

Tabs Setup / Construction / Summary. Build runs in `asyncio.to_thread`; progress and snapshots marshalled by
`loop.call_soon_threadsafe`. Canvas: box edges, solute (orange, larger), solvent (blue; solvent B green later),
depth fade; stages empty box -> solute -> insertion -> pre-opt iterations -> final. Redraw at most ~4/s and only
for a new snapshot; above 4000 atoms solvent drawn as one dot per molecule.

Outputs unchanged in name and format (`single_solute_solvent_box.xyz`, `.txt`, `single_solute_solvent.txt`,
`single_solute_box_cmass.xyz`); new `single_solute_solvent_box_preopt.xyz` (wrapped, cell on the comment line);
summary gains a properties block and pre-opt statistics.

## 6. Phase 2

- Mixture: same as phase 1 with per-species charges q_A, q_B, q_solute.
- Shared wall: cell `a x b x (c1 + c2)`; flat-bottom restraint keeps each molecule's COM in its own slab during
  pre-optimization; charges per box and species; one total multiplicity.

## 7. Validation

Tests: legacy-vs-new packing identical for fixed seeds; finite-difference gradients (rigid, flexible); minimum
image vs brute-force image scan; box-shift invariance; straddling molecule wrapped whole; COMs in [0, L);
rigid internal distances preserved; electron/multiplicity matrix; draw list headless; `_collect_inputs` supplies
every engine parameter. Benchmark script: legacy vs numpy vs numba at 15/20/30 A.
