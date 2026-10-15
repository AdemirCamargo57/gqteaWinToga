# Solvation-box builders: consolidation assessment

Status: assessment and recommendation only (no code changed by this document). Date: 2026-10-15.
Scope: `single_solute_solvent_box.py`, `mixtureSolventBox.py`, `sharedWallBox.py`, after the
pre-optimization / properties / construction-view work on branch `feature/solvation-box-preopt`.

## 1. Question

Should the three builders be merged into one Python module, or refactored to share a common
computational core? What is duplicated, what is genuinely module-specific, what could break, and in
what order should any change be made?

## 2. Where things stand

| Module | Lines | Owns |
| --- | --- | --- |
| `single_solute_solvent_box.py` | 1349 (~120 of them commented-out "previous version" blocks) | one solvent, optional centred solute, `*_cmass.xyz` + parameters file |
| `mixtureSolventBox.py` | 981 | two solvents with a composition target, optional centred solute |
| `sharedWallBox.py` | 1190 | two z-stacked regions, each a mixture with its own solute, seeds and parameters |
| `solvationBoxCommon.py` | 568 | **already shared**: pre-optimization driver, system properties, progress, the Charge/Pre-optimization/Advanced sections, Construction/Summary tabs, run pipeline |
| `solventPacking.py`, `periodicPreOptimizer.py`, `systemProperties.py`, `boxConstructionView.py` | 1488 | **already shared**: cell list + clash kernels, PBC pre-optimizer, electrons/multiplicity, construction view |

So the newer features are already shared once. What is still duplicated is the **original packing
code**, which each builder kept exactly as it was so that a fixed seed still gives the same box.

## 3. Duplicated logic (still three copies)

| Concern | single | mixture | shared wall | Same behaviour? |
| --- | --- | --- | --- | --- |
| `read_input` validation closure inside `read_params` | yes | yes | yes (two variants) | nearly; messages differ |
| XYZ reading | two strict readers (count check, line numbers) | `_read_xyz_file` (no count check) | `_read_xyz_file` (no count check) | **no** - single validates more |
| `DEFAULT_VDW_RADII` table | copy | copy | copy | duplicates `help.AtomicData.vdw_radii` |
| `get_vdw_radius` | `atomic_vdw_radii` -> `vdw_radii` -> table -> 1.70 | `vdw_radii` -> table -> 1.70 | same as mixture | same numbers today |
| `center_molecule` | raises on zero mass | falls back to the centroid | falls back to the centroid | **no** |
| bounding radius used by the centre pre-filter | includes the atom's vdW radius | geometric only | geometric only | **no** |
| centre placement margin | `extent + padding` | `max(padding, radius)` | per-face, per-region | **no** |
| target molecule count | `floor` | `round` | `round` | **no** |
| seeding | always (`None` reseeds from the OS) | only if a seed is given | only if a seed is given, per region | **no** |
| `random_rotate_molecule`, `translate_molecule`, `minimum_image_displacement` | copy | copy | copy | yes |
| `pair_cutoff` / `molecules_clash` | cutoff, then cutoff + extra | same | one test against `s(ri+rj)+extra` | same decisions (verified) |
| solvent-type choice | n/a | composition-steered `random.random()` | same, but skips the draw when a composition is 0 | **no** (RNG consumption differs) |
| density | copy | copy | per region | yes |
| legacy insertion loop | `fill_solvent_box_legacy` | `fill_solvent_box_legacy` | `fill_box_legacy` | structurally similar, policies differ |
| fast insertion loop (new) | `fill_solvent_box_fast` | `fill_solvent_box_fast` | `fill_box_fast` | structurally ~80 % identical |
| XYZ / summary writing | own formats | own formats | own formats + z shift | file formats differ |

The rows marked **no** are the core finding. The three loops do not merely repeat each other; they
encode **different packing policies**. Today those policies are pinned bit-for-bit by
`tests/test_*` (and were checked against the original `main` code before the merge of phase 2).
Merging them naively would change every seeded box a user has produced.

## 4. Genuinely module-specific requirements

- **Single solute**: one species; `*_cmass.xyz` and the `single_solute_solvent.txt` parameters/element
  count file; `floor` target; vdW-inclusive bounding radius.
- **Mixture**: two species and the adaptive composition choice; `round` target.
- **Shared wall**: two regions with independent density, composition, padding, vdW scale, seed and
  solute; Box 2 checked against Box 1 under the combined-cell minimum image; padding only on outer
  faces; the written structure shifted to min z = 0; pre-optimization adds a slab restraint and
  per-atom walls.
- **UIs**: three different forms (one cell vs two regions; one vs two solvents). They share sections,
  tabs and the pipeline already, but their parameter pages are genuinely different.

## 5. Options

### A. Merge into one module

One class with a mode switch (single / mixture / shared wall), one window that shows or hides sections.

- Pros: one file to search; the menu could offer one entry.
- Cons: the policy differences of section 3 become `if mode == ...` branches through the loop, the
  readers and the output; the file would be ~3000 lines; the three windows already differ enough that
  a combined form would be harder to use (Toga cannot hide widgets without re-parenting, see the
  plotter notes in CLAUDE.md); every change risks all three tools at once.
- **Not recommended.** It merges the parts that differ and does nothing extra for the parts that are
  already shared.

### B. Shared computational core + thin builders (recommended)

Keep three modules (three menu entries, three windows, three output conventions), and move the
packing loop and its helpers into one core module, with the differences expressed as **data**:

```text
solvationPackingCore.py   (new, no Toga)
  MoleculeTemplate        elements, centred coords, vdW radii, mass, bounding radius (policy-dependent)
  Region                  lower/upper bounds per axis, padding per face (shared wall: no padding at the interface)
  PackingPolicy           target rounding (floor | round), bounding radius incl. vdW?, placement margin rule,
                          seeding rule, species-choice rule, zero-mass centring rule
  pack_region(templates, composition, region, cell, existing: MoleculeGrid, policy, engine, rng, progress)
  read_xyz(path) -> strict reader (single's validation for all three)
```

Each builder becomes: its UI and parameter reading + `SINGLE_POLICY` / `MIXTURE_POLICY` /
`SHARED_WALL_POLICY` + its own output writer. The existing `solventPacking` (grid, kernels),
`periodicPreOptimizer`, `systemProperties`, `boxConstructionView` and `solvationBoxCommon` stay
where they are.

- Pros: one fast loop instead of three; the policies become visible and testable instead of being
  scattered through three files; bugs (e.g. the weaker XYZ validation of two builders) get fixed once;
  a fourth builder (e.g. a slab with a solid surface) would mainly need a new policy and `Region`.
- Cons: a careful migration (the RNG call order must be reproduced per policy); a short period where
  both old and new loops exist.

### C. Status quo

Keep the current state: shared extras, three packing implementations.

- Acceptable short-term: everything is tested and identical to the original behaviour.
- Cost: three fast loops to maintain (the main new duplication introduced by this work), three readers,
  three radius tables.

## 6. Compatibility risks (for B)

| Risk | Mitigation |
| --- | --- |
| Seeded boxes change | The bit-identity tests (legacy vs fast, all three builders, periodic/non-periodic, solute, cell-list regime) are the gate: a policy migration is done only when they still pass with `np.array_equal`. |
| RNG consumption order (type choice -> rotation -> centre; zero-composition shortcut; per-region seeding) | Encode it in the policy, not in shared code paths; a test per policy compares against the legacy loop. |
| Output formats | Add golden-file tests for the written `.xyz`/`.txt` headers before moving any writer. |
| Imports by name (`gqteaWinToga.py` imports the three `*UI` classes; tests use attribute names) | Keep class names, module names and public attributes; the core is additive. |
| Frozen builds (auto-py-to-exe) | A new module is picked up by import analysis; `solventPackingNumba` stays optional. Re-check one packaged build after the migration. |
| Behaviour that *should* change (e.g. two builders silently accept an XYZ whose atom count is wrong) | Make it a separate, announced change with its own test, not part of the refactor. |

## 7. Recommended migration (each step independently shippable)

1. **Freeze outputs.** Add golden-file tests for the three builders' written files (header lines,
   formats) next to the existing bit-identity tests.
2. **Shared pure helpers.** Move `read_xyz` (single's strict version), radius/mass lookup (using
   `help.AtomicData`, dropping the three `DEFAULT_VDW_RADII` copies), `random_rotate_molecule`,
   `translate_molecule` and the minimum image into the core. Builders delegate; all tests unchanged.
   The stricter reader for mixture/shared wall is the one intentional behaviour change - announce it.
3. **`PackingPolicy` + `pack_region`.** Implement the fast loop once; switch `single` to it with
   `SINGLE_POLICY` and require the identity tests to pass; then `mixture`; then `shared wall`
   (two calls, the second with Box 1's grid as `existing`).
4. **Retire the extra loops.** Delete the three `fill_*_fast` methods. Move the three legacy loops into
   `tests/` as reference implementations (they are the specification of each policy) and decide whether
   the UI should keep a "Legacy" engine choice; the default engine is unaffected.
5. **Housekeeping.** Remove the commented-out "previous version" blocks in
   `single_solute_solvent_box.py` (the git history keeps them) once the core is in place.
6. **Optional, later, deliberate:** unify policies where the difference is accidental (e.g. `floor` vs
   `round` for the target count, always vs optional seeding). Each such change alters seeded outputs and
   should be released as a documented behaviour change, not folded into the refactor.

Estimated effort: steps 1-2 small; step 3 is the substantive one (the loop plus three policies and their
tests); steps 4-5 small. Risk is controlled throughout by the identity tests that already exist.

## 8. Recommendation

Do **not** merge the three modules into one. Adopt **option B**: keep three thin builders (UI, parameters,
outputs) over one packing core with explicit `PackingPolicy` objects, migrated in the order of section 7 and
gated by the existing bit-identity tests. Until then, the current state (option C) is safe to ship: all new
functionality is already shared, and every builder still reproduces its original output exactly.
