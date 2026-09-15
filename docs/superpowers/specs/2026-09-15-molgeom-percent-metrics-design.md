# Percent-difference metrics, comparison plot and global error metrics for MolGeomComparator

Date: 2026-09-15
Status: approved, ready for implementation planning

## Purpose

[MolGeomComparator.py](../../../MolGeomComparator.py) currently reports, for every
geometric parameter present in two simulations, the two averages and their signed
difference in the parameter's own unit. That answers "how much did it move" but not
"how much did it move *relatively*", which is the figure a method-vs-method or
calculation-vs-experiment table is built from. It also produces no figure and no
whole-dataset figure of merit.

This work adds four things:

1. two relative-difference columns, each selected by its own control;
2. two control-panel switches plus two supporting fields;
3. a grouped bar chart with error bars and an overlaid percentage curve;
4. global MAE and RMSD header lines for the compared geometry category.

## Scope

Touched:

- `MolGeomComparator.py` — percent maths, new columns, header lines, plot call, UI controls
- `displayPlots.py` — a new shared figure type (`save_bar_comparison_plot`)
- `plotViewer.py` — rendering for that figure type
- `help.py` — `HelpGqteaWin.help_mol_geom_comparator`
- `USER_MANUAL.md` — the *Comparison of Molecular Geometric Parameters* section
- `tests/test_MolGeomComparator.py` — new coverage

Not touched: the three all-\* analysis tools, their file formats, the matching and
remapping logic, the occurrence filter, and the two unmatched output sections.

## 1. Relative-difference formulas

Two formulas, taken from the user's reference figures and made **signed** (the
figures take an absolute value; the sign is kept here so the direction of the shift
survives, matching the convention of the existing `difference_<label2>_minus_<label1>`
column — a reader can always take the absolute value, but cannot recover a sign that
was discarded).

Symmetric percent difference, for method vs method:

```
%Dr = 100 * (avg_2 - avg_1) / mean(avg_1, avg_2)
```

Percent error, for calculation vs a reference (file 1 is the reference):

```
%Error = 100 * (avg_2 - avg_1) / avg_1
```

Both are implemented as **module-level pure functions** next to `angular_difference`,
so they are unit-testable without constructing a calculator:

```python
def percent_difference(average_1: float, average_2: float, periodic: bool, floor: float) -> float
def percent_error(average_1: float, average_2: float, periodic: bool, floor: float) -> float
```

### Numerator

For a periodic kind (dihedral only) the numerator is `angular_difference(avg_1, avg_2)`,
the existing (-180, 180] wrapped difference, so a +179 deg to -179 deg shift counts as
+2 deg and not -358 deg. For bonds and bond angles the numerator is `avg_2 - avg_1`.

### Denominator and the seam

For the symmetric formula the denominator is the mean of the two averages. For a
periodic kind the **arithmetic** mean is wrong at the seam: +179 deg and -179 deg
average to 0 deg, which would make two nearly identical dihedrals report a division by
(almost) zero. Periodic kinds therefore use the **circular mean**, expressed through
the already-wrapped difference so no new trigonometry is introduced:

```python
def circular_mean(average_1: float, difference: float) -> float:
    """Midpoint of two angles on the circle, wrapped into (-180, 180]."""
```

so `circular_mean(179, angular_difference(179, -179)) == 180`, the physically correct
midpoint. Bonds and bond angles use `(avg_1 + avg_2) / 2`.

For `percent_error` the denominator is `avg_1` unchanged; it has no seam problem.

### Near-zero guard

A dihedral averaging 0.5 deg gives a denominator near zero, so a 0.4 deg shift would
print as 80 %. Both functions return `float("nan")` when `abs(denominator) < floor`.
The floor is a new field on `ParameterKind`:

| kind | `percent_floor` | rationale |
|---|---|---|
| `BOND_KIND` | `1e-6` A | pure division guard; a bond length is never near zero |
| `ANGLE_KIND` | `1.0` deg | a bond angle near zero is a degenerate geometry |
| `DIHEDRAL_KIND` | `1.0` deg | dihedrals genuinely sit near zero |

A `nan` row keeps its atoms, averages and absolute difference; only its percent cells
read `nan`. The count of such rows is reported in the file header and in the summary
text. `nan` rows are excluded from the percent aggregates and included in the
native-unit MAE/RMSD.

### Computation is unconditional

Both percentages are computed for every matched parameter and stored on
`MatchedParameter` as `percent_difference` and `percent_error`. The control-panel
switches decide only which **columns are written** and which curve the plot draws.
Keeping the maths unconditional keeps it testable and keeps selection a presentation
concern.

`MatchedParameter` gains two float fields (both may be `nan`). No existing field
changes.

## 2. Global error metrics

Computed over the matched set for the file's single geometry category (bond
distances, bond angles or dihedral angles — a comparison is always one of the three,
since `_validate_same_parameter_type` already rejects a mixed pair).

- `mae_<unit>` — mean of `abs(difference)` over all matched parameters
- `rmsd_<unit>` — square root of the mean of `difference**2` over all matched parameters
- `mae_percent_difference` / `rmsd_percent_difference` — the same two statistics over
  `abs(percent_difference)`, skipping `nan`; written only when that switch is on
- `mae_percent_error` / `rmsd_percent_error` — likewise, written only when that switch is on

`<unit>` is `angstrom` or `degrees`, from `ParameterKind.unit`. Two further lines name
what the metrics cover: `global_metric_category` (`ParameterKind.label`, i.e. bond
distance / bond angle / dihedral angle) and `global_metric_n` (the number of matched
parameters). Metric, category, sample size and unit are therefore all identifiable
from the file alone.

These are exposed as a `global_metrics()` method on `MolGeomCalculator` returning an
ordered mapping, so the writer and the summary text share one implementation.

## 3. Output file

The three-section layout, the two unmatched sections and every existing column are
unchanged. New header lines, written after `# minimum_occurrence_fraction` and before
the counts:

```
# formula_percent_difference 100*(avg_2 - avg_1)/((avg_1 + avg_2)/2)
# formula_percent_error 100*(avg_2 - avg_1)/avg_1
# percent_denominator_floor 1 degrees
# percent_undefined_rows 3
# global_metric_category dihedral angle
# global_metric_n 42
# mae_degrees 0.83100000
# rmsd_degrees 1.20400000
# mae_percent_difference 0.56100000
# rmsd_percent_difference 0.80420000
```

A `formula_*` line, and the corresponding aggregate lines, appear only when that
switch is on. `percent_denominator_floor` and `percent_undefined_rows` appear only
when at least one switch is on. For a periodic kind the
`formula_percent_difference` line reads
`100*wrap(avg_2 - avg_1)/circular_mean(avg_1, avg_2)`, so the file never
misdescribes its own numbers.

New matched-section columns, **appended** after `difference_...` and before
`source_row_1`:

| switch on | column name |
|---|---|
| symmetric | `percent_difference` |
| percent error | `percent_error_vs_<label_1>` |

`<label_1>` is the **first element of the existing `column_labels` property**, not the
raw `label_1`, so that two identically named files stay distinguishable here exactly as
they already are in the `average_*` columns — e.g. `percent_error_vs_isolated`. Names are
ASCII single tokens: the output is written with the platform default encoding, so a
literal `%Dr` label would risk a `UnicodeEncodeError` on Windows, and a `%` in a
column name trips `numpy.genfromtxt(..., names=True)`. The `formula_*` header lines
carry the exact equations instead.

Appending rather than reordering means a reader that indexes the earlier columns
positionally still works. With **both switches off** the file is byte-for-byte what
the tool writes today.

Percent cells are formatted `{value:>16.8f}`, matching the surrounding columns;
`float("nan")` formats as `nan` under that spec, which is the intended marker.

## 4. Control panel

Five additions to `layout_main_window`, placed between the minimum-occurrence field
and the output-folder field:

| widget | default | meaning |
|---|---|---|
| `switch_percent_difference` | **on** | write `percent_difference` and its aggregates |
| `switch_percent_error` | **off** | write `percent_error_vs_<label_1>` and its aggregates |
| `switch_show_plot` | **on** | open the comparison figure after a successful run |
| `textInput_plot_count` | blank -> 25 | how many parameters the figure shows |

The two percent switches are **independent**: either, both, or neither. Neither is a
supported state that reproduces today's output exactly. When both are on the plot
curve draws `%Dr` and its legend entry says so.

`textInput_plot_count` is read by a new `parse_optional_positive_int(text, default=25)`
helper alongside `parse_optional_fraction`, rejecting zero, negatives and non-integers
with a message naming the offending text.

`MolGeomComparatorUI` gains `DisplayPlots` as a base class. It already sets
`self.output_dir`, which is what the mixin writes its manifest into.

`read_params` reads the four new controls into `self.percent_modes`, `self.show_plot`
and `self.plot_count`. `percent_modes` is a tuple drawn from two module-level
constants, `PERCENT_DIFFERENCE_MODE = "percent_difference"` and
`PERCENT_ERROR_MODE = "percent_error"`, in that fixed order — so column order does not
depend on which switch the user happened to tick first, and an empty tuple is the
supported "neither" state. `workflow` passes
`percent_modes` to `MolGeomCalculator` and, after writing the results, calls
`save_bar_comparison_plot` + `display_plots()` when `show_plot` is on.

## 5. Comparison plot

One figure per run. Content:

- grouped bars, one pair per parameter: file 1 and file 2 averages, left axis in
  A or deg, labelled with the two user labels
- error bars: `+/- 1 standard deviation`, the `std_dev` column already parsed into
  `MatchedParameter.std_1` / `std_2`, drawn as caps on each bar
- x tick labels: `MatchedParameter.atom_label` (`C1-C2`, `H5-C2-C3-H6`), rotated so
  they stay readable
- percentage curve: line with markers on a **twin right axis** labelled in %, drawn
  from the selected formula, with a dashed horizontal zero reference line
- a single merged legend carrying both bar groups and the curve
- title stating the scope, e.g. `Top 25 of 187 matched bond distances`

Parameters are ranked by `abs(difference)`, widest shift first, and the top
`plot_count` are shown; `largest_shifts(count)` already does this ranking. When the
matched set is smaller than `plot_count` every parameter is shown and the title says
so.

`save_png=False`, per the project's *Interactive figures without PNG clutter* recipe,
so the figure is shown only through the interactive viewer and no PNG is left next to
the user's data. When no percent switch is on, the figure is drawn with bars only and
no twin axis.

### New shared figure type

`DisplayPlots.save_bar_comparison_plot(k, categories, groups, line, xlabel, ylabel,
title, save_png=True)` records a manifest entry:

```json
{
  "type": "bars",
  "categories": ["C1-C2", "C2-O3"],
  "groups": [
    {"label": "isolated", "values": [1.52, 1.43], "errors": [0.01, 0.02]},
    {"label": "solvated", "values": [1.54, 1.42], "errors": [0.02, 0.01]}
  ],
  "line": {"label": "%Dr", "values": [1.25, -0.84], "ylabel": "percent difference (%)"},
  "xlabel": "parameter",
  "ylabel": "bond distance (A)",
  "title": "Top 25 of 187 matched bond distances"
}
```

`line` is optional. `errors` is optional per group. `plotViewer.build_figures`
dispatches on `fig.get("type") == "bars"` **before** its existing `series` and `x`/`y`
branches. Existing manifests carry no `type` key, so they render exactly as now — the
extension is purely additive, which is what keeps this safe for every other tool that
shares the mixin.

The static-PNG branch inside `save_bar_comparison_plot` mirrors the viewer rendering,
so the `_display_static` fallback still produces a meaningful image for any future
caller that passes `save_png=True`.

## 6. Summary text

`summary_text` gains a metrics block after the matched/unmatched counts:

```
Global metrics over 42 matched bond distances:
  MAE  0.00831 A     RMSD  0.01204 A
  MAE  0.5610 %Dr    RMSD  0.8042 %Dr
  3 parameters have an undefined percentage (denominator below 1 deg)
```

and the *Largest shifts* table gains a percent column for each selected formula.

## 7. Testing

Tests first, following the existing file's fixture style (`tmp_path`, written
parameter files, the logic class driven directly).

Percent maths:
- symmetric formula against hand-computed values for a bond pair
- swapping the two files negates the symmetric percentage and leaves its magnitude
  unchanged
- percent error uses file 1 as denominator and is *not* symmetric
- dihedral numerator wraps across the seam
- circular mean at the seam: +179 / -179 gives a 180 deg denominator, not 0
- the near-zero guard returns `nan` for each kind at its own floor, and does not
  trigger for a normal bond, angle or dihedral

Output:
- each switch combination writes exactly the expected columns; neither switch on
  reproduces the pre-change file byte-for-byte
- `formula_*` lines match the formula actually applied, including the periodic variant
- `percent_undefined_rows` counts the guarded rows
- MAE and RMSD match hand-computed values in native units, and the percent aggregates
  skip `nan` rows
- `global_metric_category` and the unit suffix follow the parameter kind

Plot:
- the manifest entry has `type == "bars"`, the two groups in file order with their
  errors, and the line when a formula is selected
- no `line` key when neither switch is on
- ranking is by `abs(difference)` and honours `plot_count`; a matched set smaller than
  `plot_count` shows everything
- `glob("*.png") == []` in the output dir after a run
- an existing single-curve and an existing `series` manifest still render through
  `build_figures` unchanged

UI plumbing, headless via `__new__` plus stub widgets, as the file already does:
- `read_params` maps the switches onto `percent_modes`
- `parse_optional_positive_int` accepts blank, rejects zero, negatives and text
- `show_plot` off means no manifest is written

Run with `python -m pytest tests/ -q` from `venv/src/`.

## 8. Documentation

- `help.py`: extend `help_mol_geom_comparator` with the two formulas, the checkbox
  rules, the near-zero guard, the error-bar meaning and the MAE/RMSD lines.
- `USER_MANUAL.md`: update the *Comparison of Molecular Geometric Parameters* section
  (line 512) with the new controls and the annotated output header.
- `CLAUDE.md`: extend the `MolGeomComparator.py` tool note, and the `displayPlots.py` /
  `plotViewer.py` notes with the `bars` figure type.

## Decisions recorded

| Question | Decision |
|---|---|
| Sign of the percentages | Signed, file 2 minus file 1 |
| Angles and dihedrals | Computed, with a per-kind near-zero floor giving `nan` |
| Checkbox behaviour | Independent; either, both or neither |
| Error-bar source | `+/- 1` standard deviation from the source files |
| Global metric | MAE and RMSD, native units always, percent when selected |
| Plot scope | Top N by `abs(difference)`, N editable, default 25 |
| Column naming | ASCII tokens plus `formula_*` header lines |
| Plot layout | Grouped bars, twin right axis, marker line |
| Periodic denominator | Circular mean, so the +/-180 seam does not divide by zero |
