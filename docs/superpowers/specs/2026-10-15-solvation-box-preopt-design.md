# Solvation box: periodic pre-optimization, system properties, construction canvas, acceleration

Status: approved 2026-10-15. Branch `feature/solvation-box-preopt`.

## 0. User decisions

| Question | Decision |
| --- | --- |
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

## 6a. Changes made during phase-1 implementation (deviations from the approved design)

1. **Pre-optimization target scale.** As designed (sigma = the *packing* vdW scale) the
   pre-optimization was a guaranteed no-op: packing already enforces exactly that criterion
   (a live run reported "no overlaps to relax"). Added *Target contact scale*
   (`preopt_contact_scale`, default `max(0.90, vdw_scale)`); the minimizer aims at
   `1.01 * sigma` while reports use sigma, so no pair is left a hair inside the cutoff.
2. **Canvas -> bitmap view.** `toga.Canvas` on WinForms replays every shape through pythonnet on
   each repaint; with the packing thread holding the GIL the window froze for up to 9.5 s.
   The view now rasterises with Pillow and copies pixels into a native `System.Drawing.Bitmap`
   (`toga.Image(pil)` alone costs ~170 ms). Rotation moved from mouse drag to buttons
   (ImageView has no pointer events). Repaints are throttled adaptively (<= ~25 % of the UI thread).
3. **Cancel button** added (not in the design), so a long build can be stopped; no files are
   written for a cancelled run.

## 6b. Measured results (phase 1)

Water, density 0.997, 100 000 attempts, seed 7; every engine produced the identical box.

| Box | Placed | Legacy | NumPy | Numba | Pre-opt rigid | Pre-opt flexible |
| --- | --- | --- | --- | --- | --- | --- |
| 15 A | 100/112 | 54.2 s | 8.6 s (6x) | 1.8 s (31x) | 0.04 s, 170 -> 0 clashes | 0.10 s, 170 -> 0 |
| 20 A | 234/266 | 64.7 s | 9.4 s (7x) | 1.9 s (34x) | 0.08 s, 384 -> 0 | 0.25 s, 384 -> 0 |
| 30 A | 747/899 | 128.8 s | 12.1 s (11x) | 4.1 s (32x) | 0.28 s, 1055 -> 0 | 0.76 s, 1055 -> 0 |
| 40 A | 1670/2133 | -- | 12.9 s | 4.6 s | 0.54 s, 2088 -> 0 | 1.64 s, 2088 -> 0 |

(pre-opt: packing scale 0.8, target 0.9; rigid internal change <= 1e-14 A, flexible <= 1e-6 A;
connectivity preserved in every run). Live window, 30 A box: build 7.4 s (numba), UI event-loop
lag max 77 ms / mean 17 ms (was 9.5 s / 1.4 s with toga.Canvas).

Numba verdict: worthwhile -- 3-5x over the numpy cell list on top of the algorithmic gain, kept
optional (`Automatic` engine uses it when importable), bit-identical decisions.

## 6c. Phase 2 (mixture, shared wall) -- implemented and measured

- The phase-1 feature code moved into `solvationBoxCommon.py` (engine + UI mixins with geometry hooks);
  single_solute_solvent_box.py was switched to it with its tests unchanged.
- Mixture: per-species charges (solute, A, B); kinds `solvent_A`/`solvent_B` coloured on the canvas.
- Shared wall: per-box charges, one multiplicity; pre-optimization over the combined cell with a
  `SlabRestraint` (flat-bottom on each molecule's COM z) keeping molecules in their own region, and
  per-atom non-periodic walls (outer faces only, none at the interface). Packed file keeps its original
  z shift; `*_preopt.xyz` is in cell coordinates. Dead code (`fill_box_sync`, dict grid, publish/apply
  progress) removed.
- Equivalence: new legacy, numpy and numba engines of all three builders were checked bit-for-bit
  against the original `main` modules (periodic/non-periodic, with/without solutes, cell-list regime).

Water/methanol 60:40, 50 000 attempts per box, fixed seeds (identical output for every engine):

| Builder | Box | Legacy | NumPy | Numba |
| --- | --- | --- | --- | --- |
| mixture | 20 A | 14.4 s | 3.0 s (5x) | 0.62 s (23x) |
| mixture | 30 A | 27.8 s | 3.5 s (8x) | 1.01 s (28x) |
| shared wall | 20 x 20 x (20+20) A | 60.9 s | 6.8 s (9x) | 2.0 s (30x) |
| shared wall | 26 x 26 x (26+26) A | 114.0 s | 7.6 s (15x) | 2.4 s (48x) |

Live windows (flexible pre-optimization, 26-28 A cells): every build completed with 0 clashing pairs left,
connectivity preserved, UI event-loop lag <= 0.11 s.

## 7. Validation

Tests: legacy-vs-new packing identical for fixed seeds; finite-difference gradients (rigid, flexible); minimum
image vs brute-force image scan; box-shift invariance; straddling molecule wrapped whole; COMs in [0, L);
rigid internal distances preserved; electron/multiplicity matrix; draw list headless; `_collect_inputs` supplies
every engine parameter. Benchmark script: legacy vs numpy vs numba at 15/20/30 A.
