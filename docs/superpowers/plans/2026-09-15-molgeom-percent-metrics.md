# MolGeomComparator Percent Metrics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add signed relative-difference columns, their control-panel switches, a grouped bar chart with error bars and an overlaid percentage curve, and global MAE/RMSD header lines to the Molecular Geometry Comparator.

**Architecture:** The percent maths lives in module-level pure functions next to the existing `angular_difference`, computed unconditionally onto `MatchedParameter`; the switches decide only what is *written* and what the plot *draws*. The figure is a new, purely additive `{"type": "bars"}` entry in the shared plot manifest, so `MolGeomComparatorUI` gains the `DisplayPlots` mixin and every other tool is untouched.

**Tech Stack:** Python 3.13, Toga (WinForms backend), numpy, matplotlib (Agg in-process, TkAgg in the viewer subprocess), pytest.

**Spec:** [docs/superpowers/specs/2026-09-15-molgeom-percent-metrics-design.md](../specs/2026-09-15-molgeom-percent-metrics-design.md)

## Global Constraints

- Run every command from `venv/src/` with the venv activated (`..\Scripts\Activate.ps1`). Test command: `python -m pytest tests/ -q`.
- There is no linter and no build step. Plain Python source tree.
- The output file is written with `open(path, "w")` and the **platform default encoding**. Every string written to it must be **ASCII**. Column names must be single whitespace-free tokens.
- Percentages are **signed** (file 2 minus file 1), never absolute.
- Percent values are `float("nan")` when the denominator is below the kind's floor: `1e-6` for `BOND_KIND`, `1.0` for `ANGLE_KIND` and `DIHEDRAL_KIND`.
- Both percentages are computed for every matched parameter regardless of which switches are on.
- With **both switches off** the output file must be byte-for-byte what the tool writes today. `MolGeomCalculator`'s `percent_modes` therefore **defaults to `()`** — the back-compatible state — and the UI always passes its selection explicitly. `tests/test_MolGeomComparator.py:863` asserts the exact matched-column list and must keep passing untouched.
- New matched-section columns are **appended** after `difference_...` and before `source_row_1`. No existing column moves.
- Numeric cells in the matched section use `f"{value:>16.8f}"`; header metric values use `f"{value:.8f}"`.
- Existing plot manifest entries carry no `"type"` key and must keep rendering exactly as they do now.
- Tests never create a Toga window: drive `MolGeomCalculator` directly, or build the UI with `MolGeomComparatorUI.__new__` plus stub widgets.
- Commit after every task.

---

### Task 1: Percent maths as pure functions

**Files:**
- Modify: `MolGeomComparator.py` (the `ParameterKind` dataclass at lines 39-76, the three kind instances at lines 82-119, and the pure-helpers section after `angular_difference` at line 454)
- Test: `tests/test_MolGeomComparator.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `ParameterKind.percent_floor: float` — a new trailing field on the frozen dataclass.
  - `PERCENT_DIFFERENCE_MODE = "percent_difference"`, `PERCENT_ERROR_MODE = "percent_error"`, `PERCENT_MODES = (PERCENT_DIFFERENCE_MODE, PERCENT_ERROR_MODE)` — module-level constants, fixed order.
  - `circular_mean(average_1: float, difference: float) -> float`
  - `percent_difference(kind: ParameterKind, average_1: float, average_2: float) -> float`
  - `percent_error(kind: ParameterKind, average_1: float, average_2: float) -> float`

Both percent functions take the `kind` rather than loose `periodic`/`floor` arguments, so a caller cannot pair a periodic numerator with a non-periodic floor.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_MolGeomComparator.py`. Extend the existing import block at line 23 with the six new names, then append this section after the `angular_difference` tests (around line 455):

```python
# --------------------------------------------------------------------------- #
# Relative-difference formulas                                                  #
# --------------------------------------------------------------------------- #
def test_percent_difference_is_signed_and_relative_to_the_mean():
    # mean 1.5, difference +0.1 -> 100 * 0.1 / 1.5
    assert percent_difference(BOND_KIND, 1.45, 1.55) == pytest.approx(100 * 0.1 / 1.5)
    assert percent_difference(BOND_KIND, 1.55, 1.45) == pytest.approx(-100 * 0.1 / 1.5)


def test_percent_difference_is_symmetric_in_magnitude():
    """Swapping the two files must only flip the sign -- that is the point of
    the symmetric formula: neither simulation is privileged."""
    forward = percent_difference(BOND_KIND, 1.4321, 1.5987)
    backward = percent_difference(BOND_KIND, 1.5987, 1.4321)
    assert forward == pytest.approx(-backward)


def test_percent_error_uses_file_1_as_the_reference():
    # 100 * (1.55 - 1.45) / 1.45
    assert percent_error(BOND_KIND, 1.45, 1.55) == pytest.approx(100 * 0.1 / 1.45)


def test_percent_error_is_not_symmetric():
    """Unlike the symmetric formula, swapping the files changes the magnitude."""
    forward = percent_error(BOND_KIND, 1.45, 1.55)
    backward = percent_error(BOND_KIND, 1.55, 1.45)
    assert abs(forward) != pytest.approx(abs(backward))


def test_dihedral_percentages_use_the_wrapped_numerator():
    """+179 -> -179 is a +2 deg shift, not -358 deg."""
    value = percent_difference(DIHEDRAL_KIND, 179.0, -179.0)
    # numerator +2, circular-mean denominator 180
    assert value == pytest.approx(100 * 2.0 / 180.0)


def test_circular_mean_takes_the_midpoint_on_the_circle():
    assert circular_mean(179.0, angular_difference(179.0, -179.0)) == pytest.approx(180.0)
    assert circular_mean(10.0, angular_difference(10.0, 20.0)) == pytest.approx(15.0)


def test_symmetric_dihedral_denominator_does_not_collapse_at_the_seam(tmp_path):
    """Regression: the arithmetic mean of +179 and -179 is 0, which would make
    two nearly identical dihedrals divide by zero."""
    assert percent_difference(DIHEDRAL_KIND, 179.0, -179.0) == pytest.approx(1.1111111, abs=1e-6)


def test_percent_is_undefined_when_the_denominator_falls_below_the_floor():
    import math
    # A dihedral averaging a few tenths of a degree: floor is 1.0 deg.
    assert math.isnan(percent_difference(DIHEDRAL_KIND, 0.4, 0.8))
    assert math.isnan(percent_error(DIHEDRAL_KIND, 0.4, 0.8))
    # A bond never trips its 1e-6 A floor.
    assert not math.isnan(percent_difference(BOND_KIND, 1.4, 1.5))
    # A normal angle does not trip the 1.0 deg floor either.
    assert not math.isnan(percent_difference(ANGLE_KIND, 109.5, 111.0))


def test_each_kind_carries_its_own_percent_floor():
    assert BOND_KIND.percent_floor == pytest.approx(1e-6)
    assert ANGLE_KIND.percent_floor == pytest.approx(1.0)
    assert DIHEDRAL_KIND.percent_floor == pytest.approx(1.0)


def test_percent_modes_have_a_fixed_order():
    """Column order must not depend on which switch the user ticked first."""
    assert PERCENT_MODES == (PERCENT_DIFFERENCE_MODE, PERCENT_ERROR_MODE)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_MolGeomComparator.py -q -k "percent or circular"`
Expected: FAIL at import — `ImportError: cannot import name 'percent_difference' from 'MolGeomComparator'`.

- [ ] **Step 3: Add `percent_floor` to the kind table**

In `MolGeomComparator.py`, add the field to the `ParameterKind` dataclass (after `periodic: bool`, line 56):

```python
    periodic: bool
    percent_floor: float
```

and give each of the three instances its value (each is already built with keyword arguments, so append one line before the closing parenthesis):

```python
BOND_KIND = ParameterKind(
    ...
    periodic=False,
    # A bond length is never near zero, so this is only a division guard.
    percent_floor=1e-6,
)

ANGLE_KIND = ParameterKind(
    ...
    periodic=False,
    # A bond angle near zero is a degenerate geometry, not a measurement.
    percent_floor=1.0,
)

DIHEDRAL_KIND = ParameterKind(
    ...
    periodic=True,
    # Dihedrals genuinely sit near zero, where a percentage is meaningless.
    percent_floor=1.0,
)
```

- [ ] **Step 4: Write the pure functions**

Insert immediately after `angular_difference` (line 454) in `MolGeomComparator.py`:

```python
# The two relative-difference formulas the comparator can report. The order of
# PERCENT_MODES is the order the columns are written in, so it does not depend
# on which switch the user happened to tick first.
PERCENT_DIFFERENCE_MODE = "percent_difference"
PERCENT_ERROR_MODE = "percent_error"
PERCENT_MODES = (PERCENT_DIFFERENCE_MODE, PERCENT_ERROR_MODE)


def circular_mean(average_1: float, difference: float) -> float:
    """The midpoint of two angles on the circle, wrapped into (-180, 180].

    Taken through the already-wrapped ``difference`` rather than the two angles,
    so no new trigonometry is introduced. The arithmetic mean is wrong at the
    seam: +179 deg and -179 deg average to 0 deg, which would make two nearly
    identical dihedrals divide by (almost) zero.
    """
    midpoint = (average_1 + difference / 2.0 + 180.0) % 360.0 - 180.0
    if midpoint <= -180.0:
        midpoint += 360.0
    return midpoint


def _signed_difference(kind: ParameterKind, average_1: float, average_2: float) -> float:
    """``average_2 - average_1``, wrapped for a coordinate that lives on a circle."""
    if kind.periodic:
        return angular_difference(average_1, average_2)
    return average_2 - average_1


def _as_percent(numerator: float, denominator: float, floor: float) -> float:
    """``100 * numerator / denominator``, or nan when the denominator is too small.

    A dihedral averaging 0.5 deg would otherwise turn a 0.4 deg shift into 80 %,
    which reads as a dramatic result and means nothing.
    """
    if abs(denominator) < floor:
        return float("nan")
    return 100.0 * numerator / denominator


def percent_difference(kind: ParameterKind, average_1: float, average_2: float) -> float:
    """The symmetric percent difference, signed as ``file 2 - file 1``.

    ``100 * (avg_2 - avg_1) / mean(avg_1, avg_2)``. Symmetric because neither
    simulation is the reference: swapping the two files only flips the sign.
    """
    difference = _signed_difference(kind, average_1, average_2)
    if kind.periodic:
        denominator = circular_mean(average_1, difference)
    else:
        denominator = (average_1 + average_2) / 2.0
    return _as_percent(difference, denominator, kind.percent_floor)


def percent_error(kind: ParameterKind, average_1: float, average_2: float) -> float:
    """The percent error of file 2 against file 1, signed.

    ``100 * (avg_2 - avg_1) / avg_1``. Asymmetric on purpose: file 1 is the
    reference, which is what a comparison against experiment needs.
    """
    difference = _signed_difference(kind, average_1, average_2)
    return _as_percent(difference, average_1, kind.percent_floor)
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_MolGeomComparator.py -q`
Expected: PASS, and every pre-existing test in the file still passes.

- [ ] **Step 6: Commit**

```bash
git add MolGeomComparator.py tests/test_MolGeomComparator.py
git commit -m "Add signed percent-difference and percent-error formulas

Both are computed from the parameter kind, so the periodic numerator and
the near-zero floor can never be mismatched. Dihedrals use the circular
mean as the symmetric denominator: the arithmetic mean of +179 and -179
is zero, which would make two nearly identical angles divide by zero.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: Matched parameters carry both percentages

**Files:**
- Modify: `MolGeomComparator.py` (`MatchedParameter` dataclass lines 467-488, `MolGeomCalculator._compare` lines 611-660)
- Test: `tests/test_MolGeomComparator.py`

**Interfaces:**
- Consumes: `percent_difference`, `percent_error` from Task 1.
- Produces: `MatchedParameter.percent_difference: float` and `MatchedParameter.percent_error: float`, populated for every matched parameter by `_compare`, either of which may be `nan`.

Both are computed unconditionally. Which ones are *written* is decided in Task 4; keeping the maths unconditional keeps it testable and keeps selection a presentation concern.

- [ ] **Step 1: Write the failing tests**

Append to the *Comparison* section of `tests/test_MolGeomComparator.py` (after `test_largest_shifts_ranks_by_absolute_difference`, around line 610):

```python
def test_matched_parameters_carry_both_percentages(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
    )

    calculator.run()

    match = calculator.matched[0]
    assert match.percent_difference == pytest.approx(100 * 0.1 / 1.5)
    assert match.percent_error == pytest.approx(100 * 0.1 / 1.45)


def test_percentages_are_computed_even_when_no_mode_is_selected(tmp_path):
    """The switches gate the output columns, not the arithmetic."""
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
        percent_modes=(),
    )

    calculator.run()

    assert calculator.matched[0].percent_difference == pytest.approx(100 * 0.1 / 1.5)


def test_a_near_zero_dihedral_gets_an_undefined_percentage(tmp_path):
    import math
    calculator = make_calculator(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.40)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.80)],
        writer=write_dihedral_file,
    )

    calculator.run()

    match = calculator.matched[0]
    assert math.isnan(match.percent_difference)
    assert match.difference == pytest.approx(0.40)  # the absolute shift survives


def test_the_calculator_counts_undefined_percentages(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.40),
         dihedral_row(2, 3, 4, 5, "C", "C", "H", "H", 120.0)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.80),
         dihedral_row(2, 3, 4, 5, "C", "C", "H", "H", 122.0)],
        writer=write_dihedral_file,
        percent_modes=PERCENT_MODES,
    )

    calculator.run()

    assert calculator.percent_undefined == 1
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_MolGeomComparator.py -q -k "percentages or undefined"`
Expected: FAIL — `TypeError: __init__() got an unexpected keyword argument 'percent_modes'` and `AttributeError: 'MatchedParameter' object has no attribute 'percent_difference'`.

- [ ] **Step 3: Add the fields and the constructor argument**

In `MatchedParameter` (line 467), add two fields after `difference: float`:

```python
    difference: float
    # Both relative measures are computed for every match; the selected modes
    # decide only which of them reach the output file and the plot.
    percent_difference: float
    percent_error: float
```

In `MolGeomCalculator.__init__` (line 501), add the parameter and store it, plus the counter:

```python
    def __init__(
        self,
        file_1: str,
        file_2: str,
        label_1: str = "isolated",
        label_2: str = "solvated",
        min_occurrence: float = 0.0,
        # Default to none: the back-compatible output. tests/test_MolGeomComparator.py:863
        # asserts the exact matched-column list, and the UI always passes its
        # switch selection explicitly, so nothing relies on a different default.
        percent_modes: Sequence[str] = (),
    ) -> None:
        ...
        self.min_occurrence = min_occurrence
        # Kept in the canonical order, so column order never depends on the
        # order the user ticked the switches in.
        self.percent_modes = tuple(
            mode for mode in PERCENT_MODES if mode in set(percent_modes)
        )
        ...
        self.filtered_out_2 = 0
        self.percent_undefined = 0
```

- [ ] **Step 4: Populate them in `_compare`**

In `_compare` (line 611), inside the matching loop, replace the `MatchedParameter(...)` construction so the two percentages are computed and appended:

```python
            self.matched.append(
                MatchedParameter(
                    local_atoms=identity,
                    elements=canonical_elements(kind, row_1.local_atoms, row_1.elements),
                    average_1=row_1.average,
                    std_1=row_1.std_dev,
                    occurrence_1=row_1.occurrence,
                    source_row_1=row_1.source_row,
                    average_2=row_2.average,
                    std_2=row_2.std_dev,
                    occurrence_2=row_2.occurrence,
                    source_row_2=row_2.source_row,
                    difference=difference,
                    percent_difference=percent_difference(kind, row_1.average, row_2.average),
                    percent_error=percent_error(kind, row_1.average, row_2.average),
                )
            )
```

Then, after the `only_in_file_2` assignment and before the `if not self.matched:` guard, count the undefined rows:

```python
        # A row is "undefined" when a percentage the user asked for could not be
        # formed; with no mode selected nothing was asked for, so nothing counts.
        self.percent_undefined = sum(
            1
            for match in self.matched
            if any(not math.isfinite(self.percent_of(match, mode)) for mode in self.percent_modes)
        )
```

Add `import math` to the imports at the top of the file (after `import os`), and add this small accessor to `MolGeomCalculator`, just above `_compare`:

```python
    @staticmethod
    def percent_of(match: MatchedParameter, mode: str) -> float:
        """The percentage a mode names, so writer, plot and counter agree."""
        if mode == PERCENT_DIFFERENCE_MODE:
            return match.percent_difference
        return match.percent_error
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_MolGeomComparator.py -q`
Expected: PASS, all tests in the file.

- [ ] **Step 6: Commit**

```bash
git add MolGeomComparator.py tests/test_MolGeomComparator.py
git commit -m "Compute both percentages for every matched parameter

The maths runs unconditionally and the selected modes gate only the
output, which keeps the arithmetic testable without a UI and keeps
selection a presentation concern.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: Global MAE and RMSD

**Files:**
- Modify: `MolGeomComparator.py` (add `global_metrics` to `MolGeomCalculator`, after `largest_shifts` at line 544)
- Test: `tests/test_MolGeomComparator.py`

**Interfaces:**
- Consumes: `MatchedParameter.percent_difference` / `.percent_error` and `MolGeomCalculator.percent_modes`, `percent_of` from Task 2.
- Produces: `MolGeomCalculator.global_metrics() -> Dict[str, float]`, an insertion-ordered mapping whose keys are exactly the header names the writer emits: `mae_<unit>`, `rmsd_<unit>`, then `mae_<mode>` / `rmsd_<mode>` for each selected mode in `PERCENT_MODES` order. `<unit>` is `angstrom` or `degrees`.

- [ ] **Step 1: Write the failing tests**

Append to the *Comparison* section of `tests/test_MolGeomComparator.py`:

```python
# --------------------------------------------------------------------------- #
# Global error metrics                                                          #
# --------------------------------------------------------------------------- #
def test_mae_and_rmsd_are_computed_in_the_parameters_own_unit(tmp_path):
    """Differences of +0.10 and -0.20 A: MAE 0.15, RMSD sqrt((0.01+0.04)/2)."""
    import math
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.40), bond_row(2, 3, "C", "O", 1.50)],
        [bond_row(1, 2, "C", "C", 1.50), bond_row(2, 3, "C", "O", 1.30)],
    )

    calculator.run()
    metrics = calculator.global_metrics()

    assert metrics["mae_angstrom"] == pytest.approx(0.15)
    assert metrics["rmsd_angstrom"] == pytest.approx(math.sqrt((0.01 + 0.04) / 2))


def test_angle_metrics_are_reported_in_degrees(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [angle_row(1, 2, 3, "H", "C", "H", 109.0)],
        [angle_row(1, 2, 3, "H", "C", "H", 111.0)],
        writer=write_angle_file,
    )

    calculator.run()
    metrics = calculator.global_metrics()

    assert metrics["mae_degrees"] == pytest.approx(2.0)
    assert "mae_angstrom" not in metrics


def test_percent_aggregates_appear_only_for_the_selected_modes(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
        percent_modes=(PERCENT_ERROR_MODE,),
    )

    calculator.run()
    metrics = calculator.global_metrics()

    assert "mae_percent_error" in metrics
    assert "mae_percent_difference" not in metrics


def test_percent_aggregates_use_absolute_values(tmp_path):
    """+2 % and -2 % must give an MAE of 2 %, not 0 %."""
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.00), bond_row(2, 3, "C", "O", 1.00)],
        [bond_row(1, 2, "C", "C", 1.02), bond_row(2, 3, "C", "O", 0.98)],
        percent_modes=(PERCENT_ERROR_MODE,),
    )

    calculator.run()

    assert calculator.global_metrics()["mae_percent_error"] == pytest.approx(2.0)


def test_percent_aggregates_skip_undefined_rows(tmp_path):
    """A nan must not poison the aggregate of the rows that are well defined."""
    import math
    calculator = make_calculator(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.40),
         dihedral_row(2, 3, 4, 5, "C", "C", "H", "H", 100.0)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.80),
         dihedral_row(2, 3, 4, 5, "C", "C", "H", "H", 102.0)],
        writer=write_dihedral_file,
        percent_modes=(PERCENT_ERROR_MODE,),
    )

    calculator.run()
    metrics = calculator.global_metrics()

    assert math.isfinite(metrics["mae_percent_error"])
    assert metrics["mae_percent_error"] == pytest.approx(2.0)
    # The undefined row still counts towards the native-unit metric.
    assert metrics["mae_degrees"] == pytest.approx((0.40 + 2.0) / 2)


def test_metrics_keys_follow_the_canonical_mode_order(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
        percent_modes=(PERCENT_ERROR_MODE, PERCENT_DIFFERENCE_MODE),
    )

    calculator.run()

    assert list(calculator.global_metrics()) == [
        "mae_angstrom",
        "rmsd_angstrom",
        "mae_percent_difference",
        "rmsd_percent_difference",
        "mae_percent_error",
        "rmsd_percent_error",
    ]
```

Extend the test module's import block with `PERCENT_DIFFERENCE_MODE` and `PERCENT_ERROR_MODE`.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_MolGeomComparator.py -q -k "mae or rmsd or metrics or aggregates"`
Expected: FAIL — `AttributeError: 'MolGeomCalculator' object has no attribute 'global_metrics'`.

- [ ] **Step 3: Implement `global_metrics`**

Insert into `MolGeomCalculator` after `largest_shifts` (line 544):

```python
    def global_metrics(self) -> Dict[str, float]:
        """MAE and RMSD over the whole matched set, keyed by output header name.

        The native-unit pair is always present and covers every matched
        parameter. A percent pair is added for each selected mode, over the
        absolute percentages, skipping the rows whose percentage is undefined --
        a single nan would otherwise poison the aggregate of every well-defined
        row.
        """
        metrics: Dict[str, float] = {}
        unit = self.parsed_1.kind.unit

        differences = [abs(match.difference) for match in self.matched]
        metrics[f"mae_{unit}"] = sum(differences) / len(differences)
        metrics[f"rmsd_{unit}"] = math.sqrt(
            sum(value * value for value in differences) / len(differences)
        )

        for mode in self.percent_modes:
            values = [
                abs(self.percent_of(match, mode))
                for match in self.matched
                if math.isfinite(self.percent_of(match, mode))
            ]
            if not values:
                continue
            metrics[f"mae_{mode}"] = sum(values) / len(values)
            metrics[f"rmsd_{mode}"] = math.sqrt(
                sum(value * value for value in values) / len(values)
            )

        return metrics
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_MolGeomComparator.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add MolGeomComparator.py tests/test_MolGeomComparator.py
git commit -m "Add global MAE and RMSD over the matched parameter set

Always in the parameter's own unit, plus an absolute-percentage pair for
each selected mode. Undefined percentages are skipped there but still
count towards the native-unit metric, where they are perfectly well
defined.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: Output file — header lines and percent columns

**Files:**
- Modify: `MolGeomComparator.py` (`matched_columns` lines 670-680, `write_results` lines 691-729, `_format_matched_row` lines 752-768)
- Test: `tests/test_MolGeomComparator.py`

**Interfaces:**
- Consumes: `percent_modes`, `percent_undefined`, `percent_of` (Task 2), `global_metrics()` (Task 3).
- Produces: `MolGeomCalculator.percent_column_name(mode) -> str` and `MolGeomCalculator.percent_formula(mode) -> str`, both used by the writer and reused by the help text.

`_format_matched_row` becomes an **instance method** (it is currently a `@staticmethod`) because it now needs `self.percent_modes`. Update its call site inside `write_results` from `self._format_matched_row(row_index, match)` — that call already reads correctly as a bound call, so only the decorator and signature change.

- [ ] **Step 1: Write the failing tests**

Append to the *Output file* section of `tests/test_MolGeomComparator.py`:

```python
# --------------------------------------------------------------------------- #
# Percent columns and metric header lines                                       #
# --------------------------------------------------------------------------- #
def run_with_percentages(tmp_path, rows_1, rows_2, **kwargs):
    """run_and_write with both relative measures selected.

    The calculator defaults to none of them, so that the file it writes without
    being asked is the one it has always written.
    """
    kwargs.setdefault("percent_modes", PERCENT_MODES)
    return run_and_write(tmp_path, rows_1, rows_2, **kwargs)


def test_both_modes_add_both_columns_after_the_difference(tmp_path):
    _, _, text = run_with_percentages(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
    )

    columns = section_columns(text, "matched")
    difference_at = columns.index("difference_solvated_minus_isolated")
    assert columns[difference_at + 1] == "percent_difference"
    assert columns[difference_at + 2] == "percent_error_vs_isolated"
    assert columns[difference_at + 3] == "source_row_1"


def test_one_mode_adds_only_its_own_column(tmp_path):
    _, _, text = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
        percent_modes=(PERCENT_DIFFERENCE_MODE,),
    )

    columns = section_columns(text, "matched")
    assert "percent_difference" in columns
    assert not any(c.startswith("percent_error") for c in columns)


def test_no_mode_reproduces_the_original_column_list(tmp_path):
    """Back-compatibility: with neither switch on, nothing about the file changes."""
    _, _, text = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
        percent_modes=(),
    )

    columns = section_columns(text, "matched")
    assert columns == [
        "row", "atom_i", "atom_j", "element_i", "element_j",
        "average_isolated", "std_dev_isolated", "occurrence_isolated",
        "average_solvated", "std_dev_solvated", "occurrence_solvated",
        "difference_solvated_minus_isolated", "source_row_1", "source_row_2",
    ]
    header = header_values(text)
    assert not any(key.startswith(("formula_", "mae_", "rmsd_", "global_metric",
                                   "percent_")) for key in header)


def test_percent_values_land_in_their_columns(tmp_path):
    _, _, text = run_with_percentages(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
    )

    columns = section_columns(text, "matched")
    row = section_rows(text, "matched")[0]
    assert float(row[columns.index("percent_difference")]) == pytest.approx(
        100 * 0.1 / 1.5, abs=1e-6
    )
    assert float(row[columns.index("percent_error_vs_isolated")]) == pytest.approx(
        100 * 0.1 / 1.45, abs=1e-6
    )


def test_an_undefined_percentage_is_written_as_nan(tmp_path):
    _, _, text = run_with_percentages(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.40)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.80)],
        writer=write_dihedral_file,
    )

    columns = section_columns(text, "matched")
    row = section_rows(text, "matched")[0]
    assert row[columns.index("percent_difference")] == "nan"


def test_the_column_name_follows_a_disambiguated_label(tmp_path):
    """Two files labelled the same must not produce a duplicate column name."""
    _, _, text = run_with_percentages(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
        label_1="run",
        label_2="run",
    )

    columns = section_columns(text, "matched")
    assert "percent_error_vs_run_1" in columns


def test_header_records_the_formula_that_was_applied(tmp_path):
    _, _, text = run_with_percentages(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
    )

    header = header_values(text)
    assert header["formula_percent_difference"] == "100*(avg_2-avg_1)/((avg_1+avg_2)/2)"
    assert header["formula_percent_error"] == "100*(avg_2-avg_1)/avg_1"


def test_a_dihedral_file_describes_its_wrapped_formula(tmp_path):
    """The file must never misdescribe its own numbers."""
    _, _, text = run_with_percentages(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 179.0)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", -179.0)],
        writer=write_dihedral_file,
    )

    header = header_values(text)
    assert header["formula_percent_difference"] == (
        "100*wrap(avg_2-avg_1)/circular_mean(avg_1,avg_2)"
    )


def test_header_records_the_floor_and_the_undefined_row_count(tmp_path):
    _, _, text = run_with_percentages(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.40),
         dihedral_row(2, 3, 4, 5, "C", "C", "H", "H", 120.0)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.80),
         dihedral_row(2, 3, 4, 5, "C", "C", "H", "H", 122.0)],
        writer=write_dihedral_file,
    )

    header = header_values(text)
    assert header["percent_denominator_floor"] == "1 degrees"
    assert header["percent_undefined_rows"] == "1"


def test_header_records_the_metric_category_size_and_units(tmp_path):
    _, _, text = run_with_percentages(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.40), bond_row(2, 3, "C", "O", 1.50)],
        [bond_row(1, 2, "C", "C", 1.50), bond_row(2, 3, "C", "O", 1.30)],
    )

    header = header_values(text)
    assert header["global_metric_category"] == "bond distance"
    assert header["global_metric_n"] == "2"
    assert float(header["mae_angstrom"]) == pytest.approx(0.15)
    assert float(header["rmsd_angstrom"]) == pytest.approx(0.1581138, abs=1e-6)


def test_the_output_file_is_pure_ascii(tmp_path):
    """The writer uses the platform default encoding; a stray unicode label
    would raise UnicodeEncodeError on Windows."""
    _, path, _ = run_with_percentages(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
    )

    open(path, "rb").read().decode("ascii")  # raises if anything is non-ASCII
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_MolGeomComparator.py -q -k "column or header or nan or ascii"`
Expected: FAIL — `ValueError: 'percent_difference' is not in list` and missing header keys.

- [ ] **Step 3: Add the naming helpers and extend the column list**

Insert into `MolGeomCalculator`, just after the `column_labels` property (line 668):

```python
    def percent_column_name(self, mode: str) -> str:
        """The output column a mode writes into.

        ASCII, single token: the file is written with the platform default
        encoding, and a '%' in a name breaks numpy's names=True reader. The
        error column carries the reference label, taken from column_labels so
        two identically named files stay distinguishable.
        """
        if mode == PERCENT_DIFFERENCE_MODE:
            return "percent_difference"
        return f"percent_error_vs_{self.column_labels[0]}"

    def percent_formula(self, mode: str) -> str:
        """The equation actually applied, for the file's own header."""
        periodic = self.parsed_1.kind.periodic
        if mode == PERCENT_DIFFERENCE_MODE:
            if periodic:
                return "100*wrap(avg_2-avg_1)/circular_mean(avg_1,avg_2)"
            return "100*(avg_2-avg_1)/((avg_1+avg_2)/2)"
        if periodic:
            return "100*wrap(avg_2-avg_1)/avg_1"
        return "100*(avg_2-avg_1)/avg_1"
```

Then extend `matched_columns` (line 670) so the percent columns sit between the difference and the traceback columns:

```python
    def matched_columns(self) -> List[str]:
        kind = self.parsed_1.kind
        first, second = self.column_labels
        return (
            ["row"]
            + list(kind.atom_columns)
            + list(kind.element_columns)
            + [f"average_{first}", f"std_dev_{first}", f"occurrence_{first}"]
            + [f"average_{second}", f"std_dev_{second}", f"occurrence_{second}"]
            + [f"difference_{second}_minus_{first}"]
            # Appended, never inserted earlier: a reader that indexes the
            # columns above positionally keeps working.
            + [self.percent_column_name(mode) for mode in self.percent_modes]
            + ["source_row_1", "source_row_2"]
        )
```

- [ ] **Step 4: Write the header lines and the row cells**

In `write_results` (line 691), replace the block between the `minimum_occurrence_fraction` line and the `matched_parameters` line with:

```python
            out.write(f"# minimum_occurrence_fraction {self.min_occurrence:g}\n")
            self._write_percent_header(out)
            out.write(f"# matched_parameters {len(self.matched)}\n")
```

and add the method next to `_write_source_header` (line 731):

```python
    def _write_percent_header(self, out) -> None:
        """The formula, guard and global-metric block.

        Nothing is written when no mode is selected, so the file stays exactly
        what it was before relative differences existed.
        """
        kind = self.parsed_1.kind
        if not self.percent_modes:
            return

        for mode in self.percent_modes:
            out.write(f"# formula_{mode} {self.percent_formula(mode)}\n")
        out.write(
            f"# percent_denominator_floor {kind.percent_floor:g} {kind.unit}\n"
        )
        out.write(f"# percent_undefined_rows {self.percent_undefined}\n")
        out.write(f"# global_metric_category {kind.label}\n")
        out.write(f"# global_metric_n {len(self.matched)}\n")
        for name, value in self.global_metrics().items():
            out.write(f"# {name} {value:.8f}\n")
```

Note the unit token: `kind.unit` is `angstrom` or `degrees`, so the floor line reads `# percent_denominator_floor 1 degrees`.

Finally, turn `_format_matched_row` (line 752) into an instance method that emits the selected cells:

```python
    def _format_matched_row(self, row_index: int, match: MatchedParameter) -> str:
        fields = [f"{row_index:>6d}"]
        fields += [f"{atom:>8d}" for atom in match.local_atoms]
        fields += [f"{element:>8s}" for element in match.elements]
        fields += [
            f"{match.average_1:>16.8f}",
            f"{match.std_1:>16.8f}",
            f"{match.occurrence_1:>16.8f}",
            f"{match.average_2:>16.8f}",
            f"{match.std_2:>16.8f}",
            f"{match.occurrence_2:>16.8f}",
            f"{match.difference:>16.8f}",
        ]
        # An undefined percentage formats as 'nan' under this spec, which is
        # exactly the marker the header's undefined-row count refers to.
        fields += [
            f"{self.percent_of(match, mode):>16.8f}" for mode in self.percent_modes
        ]
        fields += [
            f"{match.source_row_1:>12d}",
            f"{match.source_row_2:>12d}",
        ]
        return " ".join(fields) + "\n"
```

Remove its `@staticmethod` decorator. Leave `_format_unmatched_row` a static method and unchanged — the unmatched sections have no second value to compare against.

- [ ] **Step 5: Run the whole suite**

Run: `python -m pytest tests/ -q`
Expected: PASS, including every pre-existing `test_MolGeomComparator.py` test. If `test_output_file_can_be_reread_column_by_column` or `test_matched_rows_hold_the_values_and_the_traceback_columns` fail, the percent columns were inserted in the wrong place — they must come *after* `difference_...`.

- [ ] **Step 6: Commit**

```bash
git add MolGeomComparator.py tests/test_MolGeomComparator.py
git commit -m "Write percent columns and global metrics into the output file

Columns are appended after the difference column so positional readers of
the earlier columns keep working, and with no mode selected the file is
byte-for-byte what it was. The formula header line states the wrapped
variant for dihedrals, so the file never misdescribes its own numbers.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: Summary text

**Files:**
- Modify: `MolGeomComparator.py` (`summary_text` lines 786-832)
- Test: `tests/test_MolGeomComparator.py`

**Interfaces:**
- Consumes: `global_metrics()`, `percent_undefined`, `percent_modes`, `percent_of`.
- Produces: no new API; `summary_text(output_file, shifts=10)` keeps its signature.

- [ ] **Step 1: Write the failing tests**

Append near `test_summary_lines_report_the_counts_and_the_output_path` (line 996):

```python
def test_summary_reports_the_global_metrics(tmp_path):
    calculator, path, _ = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.40), bond_row(2, 3, "C", "O", 1.50)],
        [bond_row(1, 2, "C", "C", 1.50), bond_row(2, 3, "C", "O", 1.30)],
    )

    summary = calculator.summary_text(path)

    assert "Global metrics over 2 matched bond distance" in summary
    assert "MAE" in summary and "RMSD" in summary
    assert "0.15000" in summary


def test_summary_mentions_undefined_percentages_only_when_there_are_some(tmp_path):
    calculator, path, _ = run_with_percentages(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.40)],
        [bond_row(1, 2, "C", "C", 1.50)],
    )

    assert "undefined" not in calculator.summary_text(path)


def test_summary_counts_undefined_percentages_when_present(tmp_path):
    calculator, path, _ = run_with_percentages(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.40),
         dihedral_row(2, 3, 4, 5, "C", "C", "H", "H", 120.0)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.80),
         dihedral_row(2, 3, 4, 5, "C", "C", "H", "H", 122.0)],
        writer=write_dihedral_file,
    )

    assert "1 parameter" in calculator.summary_text(path)
    assert "undefined" in calculator.summary_text(path)


def test_summary_has_no_metric_block_when_no_mode_is_selected(tmp_path):
    """The native-unit metrics are still useful, so they stay; only the
    percentage lines disappear."""
    calculator, path, _ = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.40)],
        [bond_row(1, 2, "C", "C", 1.50)],
        percent_modes=(),
    )

    summary = calculator.summary_text(path)

    assert "MAE" in summary
    assert "%" not in summary
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_MolGeomComparator.py -q -k "summary"`
Expected: FAIL — `assert 'Global metrics over 2 matched bond distance' in summary`.

- [ ] **Step 3: Extend `summary_text`**

In `summary_text` (line 786), after the `f"Only in file 2: ..."` entry and before the blank line preceding `Output file`, insert the metric block. Replace:

```python
            f"Only in file 2: {len(self.only_in_file_2)}",
            "",
            f"Output file:\n{output_file}",
```

with:

```python
            f"Only in file 2: {len(self.only_in_file_2)}",
            "",
        ]

        metrics = self.global_metrics()
        lines.append(f"Global metrics over {len(self.matched)} matched {kind.label}s:")
        lines.append(
            f"  MAE  {metrics['mae_' + kind.unit]:.5f} {unit}"
            f"     RMSD  {metrics['rmsd_' + kind.unit]:.5f} {unit}"
        )
        for mode in self.percent_modes:
            lines.append(
                f"  MAE  {metrics['mae_' + mode]:.4f} %"
                f"     RMSD  {metrics['rmsd_' + mode]:.4f} %"
                f"   ({mode})"
            )
        if self.percent_undefined:
            plural = "" if self.percent_undefined == 1 else "s"
            lines.append(
                f"  {self.percent_undefined} parameter{plural} have an undefined "
                f"percentage (denominator below {kind.percent_floor:g} {unit})"
            )

        lines += [
            "",
            f"Output file:\n{output_file}",
```

A mode whose aggregate was skipped (every row undefined) has no `mae_<mode>` key, so guard the loop body with `if ('mae_' + mode) in metrics:` before appending.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_MolGeomComparator.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add MolGeomComparator.py tests/test_MolGeomComparator.py
git commit -m "Report the global metrics in the tool's summary text

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: A shared bar-chart figure type

**Files:**
- Modify: `displayPlots.py` (add a method after `save_multiseries_plot`, line 106)
- Modify: `plotViewer.py` (`build_figures`, lines 33-56, and the module docstring's manifest description at lines 8-13)
- Create: `tests/test_displayPlots.py`

**Interfaces:**
- Consumes: nothing from earlier tasks. This task is independent of Tasks 1-5 and could be done first.
- Produces:
  ```python
  DisplayPlots.save_bar_comparison_plot(
      self, k, categories, groups, xlabel, ylabel, title,
      line=None, save_png=True,
  )
  ```
  where `categories` is a list of `str`; `groups` is a list of `(label, values, errors)` tuples with `errors` either a sequence or `None`; and `line` is `None` or a `(label, values, y2label)` tuple. It appends one manifest entry:
  ```json
  {"type": "bars", "categories": [...],
   "groups": [{"label": "...", "values": [...], "errors": [...]}],
   "line": {"label": "...", "values": [...], "ylabel": "..."},
   "xlabel": "...", "ylabel": "...", "title": "..."}
  ```
  `"errors"` is omitted when `None`; `"line"` is omitted when `None`.

  `plotViewer.build_figures` renders that entry. Its dispatch order becomes `type == "bars"` → `"series"` → single curve.

**Why a shared type rather than a local plot:** the mixin already carries the frozen-build `--plot-viewer` dispatch and the `save_png=False` no-clutter convention. Re-implementing a figure inside one tool would duplicate both. Because existing entries carry no `"type"` key, the change is purely additive.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_displayPlots.py`:

```python
"""Tests for the shared plotting mixin and the standalone viewer.

Run from ``venv/src/``::

    python -m pytest tests/ -q

The bar-comparison figure is the third manifest shape (after the single curve
and the labelled-series overlay). These tests pin its manifest contract and,
just as importantly, that the two older shapes still render unchanged -- every
analysis tool shares this module.
"""
import glob
import os

import pytest

from displayPlots import DisplayPlots


def fresh(tmp_path):
    """A DisplayPlots with per-instance buffers.

    saved_plot_files/saved_plot_data are CLASS attributes on the mixin, so a
    test that did not shadow them would see another test's figures.
    """
    plotter = DisplayPlots()
    plotter.output_dir = str(tmp_path)
    plotter.saved_plot_files = []
    plotter.saved_plot_data = []
    return plotter


CATEGORIES = ["C1-C2", "C2-O3"]
GROUPS = [
    ("isolated", [1.52, 1.43], [0.01, 0.02]),
    ("solvated", [1.54, 1.42], [0.02, 0.01]),
]


def test_bar_plot_records_a_typed_manifest_entry(tmp_path):
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES, GROUPS, "parameter", "bond distance (A)", "t",
        save_png=False,
    )

    entry = plotter.saved_plot_data[0]
    assert entry["type"] == "bars"
    assert entry["categories"] == CATEGORIES
    assert [group["label"] for group in entry["groups"]] == ["isolated", "solvated"]
    assert entry["groups"][0]["values"] == [1.52, 1.43]
    assert entry["groups"][0]["errors"] == [0.01, 0.02]
    assert "line" not in entry


def test_bar_plot_records_the_overlaid_line_when_given(tmp_path):
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES, GROUPS, "parameter", "bond distance (A)", "t",
        line=("%dr", [1.25, -0.84], "percent difference (%)"),
        save_png=False,
    )

    line = plotter.saved_plot_data[0]["line"]
    assert line["label"] == "%dr"
    assert line["values"] == [1.25, -0.84]
    assert line["ylabel"] == "percent difference (%)"


def test_bar_plot_omits_errors_that_were_not_supplied(tmp_path):
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES, [("isolated", [1.52, 1.43], None)],
        "parameter", "d (A)", "t", save_png=False,
    )

    assert "errors" not in plotter.saved_plot_data[0]["groups"][0]


def test_bar_plot_writes_no_png_when_asked_not_to(tmp_path):
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES, GROUPS, "parameter", "d (A)", "t", save_png=False,
    )

    assert plotter.saved_plot_files == []
    assert glob.glob(os.path.join(str(tmp_path), "*.png")) == []


def test_bar_plot_can_still_write_a_static_png(tmp_path):
    """The _display_static fallback needs a real image when one is asked for."""
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES, GROUPS, "parameter", "d (A)", "t",
        line=("%dr", [1.25, -0.84], "percent (%)"), save_png=True,
    )

    assert len(plotter.saved_plot_files) == 1
    assert len(glob.glob(os.path.join(str(tmp_path), "*.png"))) == 1


def test_values_are_json_safe_floats(tmp_path):
    """numpy scalars must not reach json.dump."""
    import json
    np = pytest.importorskip("numpy")
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES,
        [("isolated", np.array([1.52, 1.43]), np.array([0.01, 0.02]))],
        "parameter", "d (A)", "t",
        line=("%dr", np.array([1.25, -0.84]), "percent (%)"),
        save_png=False,
    )

    json.dumps(plotter.saved_plot_data)  # raises TypeError on numpy scalars


# --------------------------------------------------------------------------- #
# The viewer renders all three manifest shapes                                  #
# --------------------------------------------------------------------------- #
def test_viewer_renders_a_bar_figure_with_a_twin_axis():
    import matplotlib.pyplot as plt
    import plotViewer

    plt.close("all")
    plotViewer.build_figures([{
        "type": "bars",
        "categories": CATEGORIES,
        "groups": [
            {"label": "isolated", "values": [1.52, 1.43], "errors": [0.01, 0.02]},
            {"label": "solvated", "values": [1.54, 1.42], "errors": [0.02, 0.01]},
        ],
        "line": {"label": "%dr", "values": [1.25, -0.84],
                 "ylabel": "percent difference (%)"},
        "xlabel": "parameter", "ylabel": "bond distance (A)", "title": "t",
    }])

    figure = plt.figure(1)
    # Two axes: the bars' own, plus the twin carrying the percentage line.
    assert len(figure.axes) == 2
    assert figure.axes[0].get_xticklabels()[0].get_text() == "C1-C2"
    plt.close("all")


def test_viewer_renders_a_bar_figure_without_a_line():
    import matplotlib.pyplot as plt
    import plotViewer

    plt.close("all")
    plotViewer.build_figures([{
        "type": "bars",
        "categories": CATEGORIES,
        "groups": [{"label": "isolated", "values": [1.52, 1.43]}],
        "xlabel": "parameter", "ylabel": "d (A)", "title": "t",
    }])

    assert len(plt.figure(1).axes) == 1
    plt.close("all")


def test_viewer_still_renders_the_two_older_manifest_shapes():
    """Regression guard: every other tool writes these, and they carry no
    'type' key at all."""
    import matplotlib.pyplot as plt
    import plotViewer

    plt.close("all")
    plotViewer.build_figures([
        {"x": [0, 1], "y": [0, 1], "xlabel": "x", "ylabel": "y", "title": "curve"},
        {"series": [{"x": [0, 1], "y": [1, 0], "label": "a"}],
         "xlabel": "x", "ylabel": "y", "title": "overlay"},
    ])

    assert plt.figure(1).axes[0].get_title() == "curve"
    assert plt.figure(2).axes[0].get_title() == "overlay"
    assert len(plt.figure(2).axes[0].get_legend().get_texts()) == 1
    plt.close("all")
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_displayPlots.py -q`
Expected: FAIL — `AttributeError: 'DisplayPlots' object has no attribute 'save_bar_comparison_plot'`.

- [ ] **Step 3: Implement `save_bar_comparison_plot`**

Insert into `displayPlots.py` after `save_multiseries_plot` (line 106):

```python
    def save_bar_comparison_plot(self, k, categories, groups, plot_xlabel,
                                 plot_ylabel, plot_title, line=None, save_png=True):
        """Record a grouped bar chart with error bars and an optional overlaid line.

        `categories` names the bar positions; `groups` is a list of
        (label, values, errors) tuples, one bar per category per group, with
        `errors` either a sequence of symmetric error-bar half-lengths or None.
        `line` is an optional (label, values, y2label) tuple drawn on a twin
        right-hand axis -- used to put a relative measure (a percentage) beside
        absolute values that share no scale with it.

        Pass save_png=False to skip writing a static PNG (see save_plots); the
        raw data recorded below is what the interactive viewer consumes.
        """
        entry = {
            "type": "bars",
            "categories": [str(c) for c in categories],
            "groups": [],
            "xlabel": plot_xlabel,
            "ylabel": plot_ylabel,
            "title": plot_title,
        }
        for label, values, errors in groups:
            group = {"label": label, "values": [float(v) for v in values]}
            if errors is not None:
                group["errors"] = [float(e) for e in errors]
            entry["groups"].append(group)
        if line is not None:
            line_label, line_values, line_ylabel = line
            entry["line"] = {
                "label": line_label,
                "values": [float(v) for v in line_values],
                "ylabel": line_ylabel,
            }

        if save_png:
            temp_filename = tempfile.NamedTemporaryFile(delete=False, suffix=".png", dir=self.output_dir).name
            figure = plt.figure(k)
            axes = figure.gca()
            draw_bar_comparison(axes, entry)
            plt.tight_layout()
            plt.savefig(temp_filename)
            plt.close(figure)
            self.saved_plot_files.append(temp_filename)

        self.saved_plot_data.append(entry)
```

The drawing itself is shared with the viewer, so put it at module level in `displayPlots.py`, above the class:

```python
def draw_bar_comparison(axes, figure_entry, font_style=None):
    """Draw one "bars" manifest entry onto `axes`; return the twin axis or None.

    Public (no leading underscore) because plotViewer.py imports it: it is part
    of this module's interface, not an internal detail. Kept out of the class so
    the standalone viewer can reuse the exact same
    rendering: a figure must not look different depending on which process drew
    it. The percentage line sits on a twin right-hand axis because it shares no
    scale with the absolute values in angstroms or degrees.
    """
    font_style = font_style or {'color': 'darkred', 'weight': 'normal', 'size': 14}
    categories = figure_entry.get("categories", [])
    groups = figure_entry.get("groups", [])
    positions = range(len(categories))
    count = max(len(groups), 1)
    width = 0.8 / count

    handles = []
    for index, group in enumerate(groups):
        offset = (index - (count - 1) / 2.0) * width
        bars = axes.bar(
            [p + offset for p in positions],
            group.get("values", []),
            width=width,
            yerr=group.get("errors"),
            capsize=4,
            label=group.get("label", ""),
        )
        handles.append(bars)

    axes.set_xticks(list(positions))
    axes.set_xticklabels(categories, rotation=60, ha="right", fontsize=8)
    axes.set_xlabel(figure_entry.get("xlabel", ""), fontdict=font_style)
    axes.set_ylabel(figure_entry.get("ylabel", ""), fontdict=font_style)
    axes.set_title(figure_entry.get("title", ""), fontdict=font_style)
    axes.grid(True, axis="y", alpha=0.3)

    twin = None
    line = figure_entry.get("line")
    if line is not None:
        twin = axes.twinx()
        drawn = twin.plot(
            list(positions), line.get("values", []),
            color="black", marker="o", markersize=4, linewidth=1.5,
            label=line.get("label", ""),
        )
        twin.axhline(0.0, color="grey", linestyle="--", linewidth=0.8)
        twin.set_ylabel(line.get("ylabel", ""), fontdict=font_style)
        handles.extend(drawn)

    # One legend for everything: two y axes would otherwise produce two.
    axes.legend(handles=handles, loc="best", fontsize=9)
    return twin
```

- [ ] **Step 4: Teach the viewer the new shape**

In `plotViewer.py`, import the shared renderer and dispatch on it. Add after the `FONT_STYLE` definition (line 30):

```python
from displayPlots import draw_bar_comparison
```

and make `build_figures` (line 33) dispatch before its existing branches:

```python
def build_figures(figures):
    """Create one matplotlib figure per manifest entry (no blocking show).

    Three entry shapes: a grouped bar chart with error bars and an optional
    percentage line on a twin axis ("type": "bars"), several labelled curves
    overlaid with a legend ("series"), or a single curve ("x"/"y"). Entries
    written before the bar chart existed carry no "type" key, so the plain
    curve stays the default.
    """
    for i, fig in enumerate(figures, 1):
        figure = plt.figure(i)
        if fig.get("type") == "bars":
            draw_bar_comparison(figure.gca(), fig, FONT_STYLE)
            plt.tight_layout()
            continue
        if "series" in fig:
            ...
```

Leave the rest of the loop body exactly as it is. Update the module docstring's manifest description (lines 8-13) to mention the third shape.

Note: `plotViewer.py` pins `TkAgg` at import time and `displayPlots.py` pins `Agg` at import time. `plotViewer` importing `displayPlots` therefore re-pins `Agg` *after* `TkAgg` unless the import order is right — put the `from displayPlots import ...` line **before** `matplotlib.use("TkAgg")` so the viewer's pin is the last one to run. Verify with the test in Step 5; the viewer is also exercised by `tests/test_radialDistribution.py::test_display_plots_spawns_viewer_source_and_frozen`.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/ -q`
Expected: PASS, the whole suite. Then confirm the viewer's backend pin survived the new import:

Run: `python -c "import plotViewer, matplotlib; print(matplotlib.get_backend())"`
Expected: `TkAgg` (case-insensitive). If it prints `Agg`, the `displayPlots` import is below the `matplotlib.use("TkAgg")` line — move it above.

- [ ] **Step 6: Commit**

```bash
git add displayPlots.py plotViewer.py tests/test_displayPlots.py
git commit -m "Add a grouped bar chart with error bars to the shared plot manifest

A third manifest shape, keyed by an explicit "type" that older entries do
not carry, so every existing figure renders unchanged. The renderer is
shared between the in-process PNG fallback and the standalone viewer, so
a figure cannot look different depending on which process drew it.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 7: Control panel and the comparison plot

**Files:**
- Modify: `MolGeomComparator.py` (add `plot_dataset` to `MolGeomCalculator`; `parse_optional_positive_int` near `parse_optional_fraction` line 838; `MolGeomComparatorUI` class declaration line 867, `__init__` line 872, `layout_main_window` line 884, `read_params` line 1075, `workflow` line 1132)
- Test: `tests/test_MolGeomComparator.py`

**Interfaces:**
- Consumes: `percent_modes`, `percent_of` (Task 2), `largest_shifts` (existing), `save_bar_comparison_plot` (Task 6).
- Produces:
  - `parse_optional_positive_int(text: str, default: int = 25) -> int`
  - `MolGeomCalculator.plot_dataset(count: int, mode: Optional[str]) -> Dict` returning keys `categories` (list of str), `groups` (list of `(label, values, errors)` ready for `save_bar_comparison_plot`), `percent` (list of float, or `None` when `mode` is None), `shown` (int), `total` (int).
  - UI attributes `switch_percent_difference`, `switch_percent_error`, `switch_show_plot`, `textInput_plot_count`, and the parsed `self.percent_modes`, `self.show_plot`, `self.plot_count`.

- [ ] **Step 1: Write the failing tests**

Append a plot section to `tests/test_MolGeomComparator.py`, and extend the `TestUIPlumbing.bare_ui` stub set:

```python
# --------------------------------------------------------------------------- #
# The comparison figure                                                         #
# --------------------------------------------------------------------------- #
def test_plot_dataset_ranks_by_the_size_of_the_shift(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.50), bond_row(2, 3, "C", "O", 1.40),
         bond_row(3, 4, "O", "H", 0.96)],
        [bond_row(1, 2, "C", "C", 1.51), bond_row(2, 3, "C", "O", 1.60),
         bond_row(3, 4, "O", "H", 0.99)],
    )
    calculator.run()

    dataset = calculator.plot_dataset(2, PERCENT_DIFFERENCE_MODE)

    assert dataset["categories"] == ["C2-O3", "O3-H4"]  # 0.20 then 0.03
    assert dataset["shown"] == 2
    assert dataset["total"] == 3


def test_plot_dataset_carries_both_averages_and_their_error_bars(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.50, std_dev=0.01)],
        [bond_row(1, 2, "C", "C", 1.55, std_dev=0.02)],
        label_1="iso", label_2="sol",
    )
    calculator.run()

    dataset = calculator.plot_dataset(25, None)

    (label_1, values_1, errors_1), (label_2, values_2, errors_2) = dataset["groups"]
    assert (label_1, label_2) == ("iso", "sol")
    assert values_1 == [pytest.approx(1.50)]
    assert errors_1 == [pytest.approx(0.01)]
    assert values_2 == [pytest.approx(1.55)]
    assert errors_2 == [pytest.approx(0.02)]


def test_plot_dataset_has_no_percentages_without_a_mode(tmp_path):
    calculator = make_calculator(
        tmp_path, [bond_row(1, 2, "C", "C", 1.5)], [bond_row(1, 2, "C", "C", 1.6)]
    )
    calculator.run()

    assert calculator.plot_dataset(25, None)["percent"] is None


def test_plot_dataset_shows_everything_when_the_set_is_small(tmp_path):
    calculator = make_calculator(
        tmp_path, [bond_row(1, 2, "C", "C", 1.5)], [bond_row(1, 2, "C", "C", 1.6)]
    )
    calculator.run()

    dataset = calculator.plot_dataset(25, None)
    assert dataset["shown"] == 1 and dataset["total"] == 1


# --------------------------------------------------------------------------- #
# The plot-count field                                                          #
# --------------------------------------------------------------------------- #
def test_blank_plot_count_uses_the_default():
    assert parse_optional_positive_int("") == 25
    assert parse_optional_positive_int("   ") == 25


def test_plot_count_parses_a_whole_number():
    assert parse_optional_positive_int("40") == 40


def test_plot_count_rejects_zero_negatives_and_text():
    for bad in ("0", "-3", "lots", "12.5"):
        with pytest.raises(ValueError):
            parse_optional_positive_int(bad)
```

Add these to `TestUIPlumbing`, and add the four new stub widgets to `bare_ui`:

```python
        ui.textInput_output = FakeWidget()
        ui.switch_percent_difference = FakeWidget(True)
        ui.switch_percent_error = FakeWidget(False)
        ui.switch_show_plot = FakeWidget(False)   # off: tests never spawn a viewer
        ui.textInput_plot_count = FakeWidget()
        ui.multi_line_text = FakeWidget()
```

```python
    def test_read_params_maps_the_switches_onto_the_modes(self, tmp_path):
        ui = self.bare_ui()
        ui.file_1 = write_bond_file(tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.5)])
        ui.file_2 = write_bond_file(tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.6)])
        ui.switch_percent_difference.value = True
        ui.switch_percent_error.value = True

        assert asyncio.run(ui.read_params()) is True
        assert ui.percent_modes == (PERCENT_DIFFERENCE_MODE, PERCENT_ERROR_MODE)

    def test_read_params_allows_neither_switch(self, tmp_path):
        """Neither is a supported state, not an error: it reproduces the
        original output file."""
        ui = self.bare_ui()
        ui.file_1 = write_bond_file(tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.5)])
        ui.file_2 = write_bond_file(tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.6)])
        ui.switch_percent_difference.value = False
        ui.switch_percent_error.value = False

        assert asyncio.run(ui.read_params()) is True
        assert ui.percent_modes == ()

    def test_read_params_refuses_a_bad_plot_count(self, tmp_path):
        ui = self.bare_ui()
        ui.file_1 = write_bond_file(tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.5)])
        ui.file_2 = write_bond_file(tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.6)])
        ui.textInput_plot_count.value = "-5"

        assert asyncio.run(ui.read_params()) is False
        assert ui.main_window.dialogs == ["InfoDialog"]

    def test_workflow_draws_no_figure_when_the_plot_switch_is_off(self, tmp_path):
        ui = self.bare_ui()
        ui.file_1 = write_bond_file(tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.50)])
        ui.file_2 = write_bond_file(tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.55)])
        ui.output_dir = str(tmp_path)
        ui.saved_plot_data = []
        ui.saved_plot_files = []
        ui.display_plots = lambda *a, **k: None

        asyncio.run(ui.workflow(None))

        assert ui.saved_plot_data == []

    def test_workflow_draws_the_figure_and_writes_no_png(self, tmp_path):
        import glob
        ui = self.bare_ui()
        ui.file_1 = write_bond_file(tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.50)])
        ui.file_2 = write_bond_file(tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.55)])
        ui.output_dir = str(tmp_path)
        ui.switch_show_plot.value = True
        ui.saved_plot_data = []
        ui.saved_plot_files = []
        launched = []
        ui.display_plots = lambda *a, **k: launched.append(True)

        asyncio.run(ui.workflow(None))

        assert launched == [True]
        assert ui.saved_plot_data[0]["type"] == "bars"
        assert glob.glob(os.path.join(str(tmp_path), "*.png")) == []

    def test_the_ui_is_a_display_plots_subclass(self):
        from displayPlots import DisplayPlots
        assert issubclass(MolGeomComparatorUI, DisplayPlots)
```

Extend the test module imports with `parse_optional_positive_int`.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_MolGeomComparator.py -q -k "plot or switch or neither"`
Expected: FAIL — `ImportError: cannot import name 'parse_optional_positive_int'`.

- [ ] **Step 3: Add `plot_dataset` and the field parser**

Insert into `MolGeomCalculator` after `global_metrics`:

```python
    PLOT_LABELS = {
        PERCENT_DIFFERENCE_MODE: "%Δr (symmetric)",
        PERCENT_ERROR_MODE: "%Error",
    }
    PLOT_AXIS_LABELS = {
        PERCENT_DIFFERENCE_MODE: "percent difference (%)",
        PERCENT_ERROR_MODE: "percent error (%)",
    }

    def plot_dataset(self, count: int, mode: Optional[str]) -> Dict[str, object]:
        """The figure's data: the `count` parameters that moved most.

        Ranked by the size of the shift so the figure shows what actually
        changed; a comparison can match hundreds of parameters, and a bar per
        parameter would be unreadable. Returned as plain lists so it can be
        checked without a display.
        """
        shown = self.largest_shifts(count)
        first, second = self.column_labels
        dataset = {
            "categories": [match.atom_label for match in shown],
            "groups": [
                (first,
                 [match.average_1 for match in shown],
                 [match.std_1 for match in shown]),
                (second,
                 [match.average_2 for match in shown],
                 [match.std_2 for match in shown]),
            ],
            "percent": None,
            "shown": len(shown),
            "total": len(self.matched),
        }
        if mode is not None:
            dataset["percent"] = [self.percent_of(match, mode) for match in shown]
        return dataset

    def plot_title(self, dataset: Dict[str, object]) -> str:
        kind = self.parsed_1.kind
        if dataset["shown"] >= dataset["total"]:
            return f"{dataset['total']} matched {kind.label}s"
        return f"Top {dataset['shown']} of {dataset['total']} matched {kind.label}s"
```

The plot labels may use `%Δr`: they go into a JSON manifest (escaped to ASCII by `json.dump`) and into matplotlib, never into the platform-encoded output file.

Add the field parser next to `parse_optional_fraction` (line 838):

```python
def parse_optional_positive_int(text: str, default: int = 25) -> int:
    """Read an optional whole-number field; a blank field means `default`."""
    stripped = (text or "").strip()
    if not stripped:
        return default

    try:
        value = int(stripped)
    except ValueError:
        raise ValueError(
            f"'{stripped}' is not a whole number. Enter how many parameters the "
            "plot should show, or leave the field blank for the default "
            f"({default})."
        )

    if value < 1:
        raise ValueError(
            f"The number of parameters to plot must be at least 1 (got {value})."
        )
    return value
```

- [ ] **Step 4: Wire the controls into the UI**

In `MolGeomComparator.py`:

Import the mixin next to the `help` import (line 33):

```python
from displayPlots import DisplayPlots
from help import HelpGqteaWin
```

Change the class declaration (line 867) to `class MolGeomComparatorUI(DisplayPlots):`.

In `__init__` (line 872), give the instance its own plot buffers before building the window — the mixin's are **class** attributes, so two open comparator windows would otherwise share one list:

```python
    def __init__(self, *args) -> None:
        self.file_1 = None
        self.file_2 = None
        self.parsed_1 = None
        self.parsed_2 = None
        self.output_dir = os.getcwd()
        # Per-instance: DisplayPlots keeps these on the class, so two open
        # windows would otherwise append into the same figure list.
        self.saved_plot_files = []
        self.saved_plot_data = []
        self.layout_main_window(*args)
```

In `layout_main_window`, insert four rows between the minimum-occurrence row (line 946) and the output-folder row (line 948):

```python
        self.switch_percent_difference = toga.Switch(
            "Symmetric percent difference: 100*(avg_2-avg_1)/mean", value=True
        )
        form_row("Relative difference columns:", self.switch_percent_difference)

        self.switch_percent_error = toga.Switch(
            "Percent error vs file 1: 100*(avg_2-avg_1)/avg_1", value=False
        )
        form_row("", self.switch_percent_error)

        self.switch_show_plot = toga.Switch(
            "Open the comparison figure after comparing", value=True
        )
        form_row("Plot:", self.switch_show_plot)

        self.textInput_plot_count = toga.TextInput(
            placeholder="Default: 25 (the parameters that shifted most)",
            style=input_style,
        )
        form_row("Parameters to plot:", self.textInput_plot_count)
```

The two percent switches are independent — either, both, or neither — so they are plain `Switch`es with no `on_change` coupling.

In `read_params` (line 1075), after the `min_occurrence` block:

```python
        try:
            self.plot_count = parse_optional_positive_int(self.textInput_plot_count.value)
        except ValueError as exc:
            await self.warning_function("Error", str(exc))
            return False

        selected = set()
        if self.switch_percent_difference.value:
            selected.add(PERCENT_DIFFERENCE_MODE)
        if self.switch_percent_error.value:
            selected.add(PERCENT_ERROR_MODE)
        # Canonical order, so the columns do not depend on the switching order.
        self.percent_modes = tuple(mode for mode in PERCENT_MODES if mode in selected)
        self.show_plot = bool(self.switch_show_plot.value)
```

In `workflow` (line 1132), pass the modes to the calculator:

```python
        calculator = MolGeomCalculator(
            self.file_1,
            self.file_2,
            label_1=self.label_1,
            label_2=self.label_2,
            min_occurrence=self.min_occurrence,
            percent_modes=self.percent_modes,
        )
```

and, after `self.multi_line_text.value = calculator.summary_text(output_file)` at the end, draw the figure:

```python
        if self.show_plot:
            self._show_comparison_plot(calculator)

    def _show_comparison_plot(self, calculator: MolGeomCalculator) -> None:
        """Bars for both averages with their spread, and the selected percentage.

        save_png=False per the project's interactive-figures recipe: the figure
        is shown only through the viewer, leaving no image files next to the
        user's parameter files.
        """
        # With both modes selected the curve draws the symmetric one, which is
        # the first in PERCENT_MODES; its legend entry names it.
        mode = calculator.percent_modes[0] if calculator.percent_modes else None
        dataset = calculator.plot_dataset(self.plot_count, mode)
        kind = calculator.kind

        line = None
        if mode is not None:
            line = (
                calculator.PLOT_LABELS[mode],
                dataset["percent"],
                calculator.PLOT_AXIS_LABELS[mode],
            )

        self.save_bar_comparison_plot(
            1,
            dataset["categories"],
            dataset["groups"],
            "parameter",
            f"{kind.label} ({kind.unit_symbol})",
            calculator.plot_title(dataset),
            line=line,
            save_png=False,
        )
        self.display_plots()
```

Bars are labelled `isolated` / `solvated` from `column_labels`, error bars are `+/- 1` standard deviation from the source files, and the legend carries all three entries — `draw_bar_comparison` merges them onto one legend.

- [ ] **Step 5: Run the whole suite**

Run: `python -m pytest tests/ -q`
Expected: PASS.

- [ ] **Step 6: Check the real window opens**

Run: `python gqteaWinToga.py`, open **Tools -> Comparison of Molecular Geometric Parameters**, and confirm the four new controls appear and are readable. If two combined parameter files are at hand, run a comparison with the plot switch on and confirm the figure opens with grouped bars, error-bar caps, rotated atom labels, a right-hand percent axis and one legend. Close without committing anything if the layout needs adjusting.

- [ ] **Step 7: Commit**

```bash
git add MolGeomComparator.py tests/test_MolGeomComparator.py
git commit -m "Add the percent switches, plot controls and comparison figure

The two percent switches are independent: either, both, or neither, and
neither reproduces the original output file. The figure ranks parameters
by the size of the shift because a comparison can match hundreds of them,
and passes save_png=False so no image files are left beside the data.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 8: Documentation

**Files:**
- Modify: `help.py` (`HelpGqteaWin.help_mol_geom_comparator`, from line 298)
- Modify: `USER_MANUAL.md` (*Comparison of Molecular Geometric Parameters*, from line 512)
- Modify: `CLAUDE.md` (the `MolGeomComparator.py` tool note, and the `displayPlots.py` / `plotViewer.py` mixin notes)

**Interfaces:**
- Consumes: everything from Tasks 1-7. No code interfaces produced.

- [ ] **Step 1: Read the current text**

Read `help.py` from line 298 to the end of `help_mol_geom_comparator`, and `USER_MANUAL.md` lines 512 onwards, so the additions match the existing voice and heading depth.

- [ ] **Step 2: Extend the help text**

Append a section to `help_mol_geom_comparator` covering, in the string's existing style:

- the two formulas, written out, and that both are **signed** as file 2 minus file 1;
- that the switches are independent — either, both, or neither, and that neither leaves the output exactly as it was;
- that with both on, the plot curve draws the symmetric difference;
- the near-zero guard: a percentage whose denominator falls below `1e-6` A (bonds) or `1.0` deg (angles and dihedrals) is written `nan` and counted in the `# percent_undefined_rows` header, because a dihedral averaging a few tenths of a degree would otherwise report an enormous and meaningless percentage;
- that dihedral differences are wrapped into (-180, 180] and their symmetric denominator is the circular mean, so the +/-180 seam does not collapse it;
- the `# mae_*` / `# rmsd_*` header lines, what they cover (`# global_metric_category`, `# global_metric_n`) and their units;
- the figure: grouped bars for both averages, `+/- 1` standard deviation error bars taken from the source files, the percentage on the right-hand axis, and the *Parameters to plot* field selecting how many of the largest shifts are shown.

Keep it ASCII — the string is displayed in a Toga `MultilineTextInput`, but staying ASCII matches the rest of the file.

- [ ] **Step 3: Update the user manual**

In the *Comparison of Molecular Geometric Parameters* section of `USER_MANUAL.md`, add the four new controls to the controls list and show an annotated example of the new header block and the two new columns, following the formatting the surrounding tool sections already use.

- [ ] **Step 4: Update CLAUDE.md**

Extend the `MolGeomComparator.py` bullet with: the two signed formulas and the independent switches; the per-kind `percent_floor` and the `nan` guard; the circular-mean denominator for periodic kinds and *why* (the arithmetic mean of +179 and -179 is 0); that percentages are computed unconditionally and the modes gate only output; that new columns are appended after `difference_...` so positional readers keep working, and that no selected mode reproduces the original file; the `mae_*`/`rmsd_*` header block; and that the tool now mixes in `DisplayPlots` and passes `save_png=False`.

Extend the `displayPlots.py` bullet with `save_bar_comparison_plot` and the `{"type": "bars"}` manifest shape, noting that the renderer `draw_bar_comparison` is shared with `plotViewer.py` so a figure looks the same from either process, that older entries carry no `"type"` key and are unaffected, and that `plotViewer.py` must import `displayPlots` **before** its `matplotlib.use("TkAgg")` line or the viewer loses its interactive backend.

- [ ] **Step 5: Verify nothing broke**

Run: `python -m pytest tests/ -q`
Expected: PASS.

Run: `python -c "from help import HelpGqteaWin; print(len(HelpGqteaWin.help_mol_geom_comparator))"`
Expected: a number larger than before, and no traceback.

- [ ] **Step 6: Commit**

```bash
git add help.py USER_MANUAL.md CLAUDE.md
git commit -m "Document the percent metrics, the comparison figure and the bar manifest

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Verification

After Task 8, confirm the whole feature end to end:

- [ ] `python -m pytest tests/ -q` — the full suite passes, not just the two touched files.
- [ ] `python -c "import plotViewer, matplotlib; print(matplotlib.get_backend())"` prints `TkAgg`.
- [ ] `git status` is clean.
- [ ] Launch `python gqteaWinToga.py`, run a real comparison of two combined parameter files with both switches on, and check: the output file carries both percent columns and the metric header block; the figure opens with bars, caps, the right-hand percent axis and one legend; no `.png` file appears next to the parameter files.
- [ ] Re-run with both switches off and confirm the output file has no percent columns and no metric header lines.
