# Periodic-Table Element Defaults for the USPP Generator — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Status:** APPROVED by the user on 2026-10-30 (rev. 2, with the user's decisions folded in). Implemented on branch `feature/runatom-element-defaults`.

**Goal:** In *Tools > Vanderbilt USPP Generator*, the user picks an element (Z ≤ 92) in a periodic table and every `.adat` field is filled with suggested values. Every field stays editable, and the fill can be undone.

**Architecture:** Three sources of defaults, chosen per element, best first:
1. **Library case**: a tested case from the bundled `Work/` library (45 elements), loaded through the existing `read_case_folder` → `case_to_fields` → `_fill` path.
2. **Same-group template**: for an element with no library case, a starter All-electron input plus the Pseudopotential-tab values of the nearest library element **in the same group**. Projector shell indices are remapped. A warning is shown.
3. **Computed starter**: for the two groups with no library member (7 and 18), the All-electron input only, with the radii left blank and a warning.

A new Toga-free module, `runatomDefaults.py`, holds the rules. The periodic-table layout moves to a Toga-free `periodicTableData.py`, which `periodicTable.py` re-exports, and the picker window gains three optional keyword arguments.

**Tech stack:** Python 3.13, Toga 0.5 (WinForms), pytest. No new dependencies.

**Spec:** this plan; parent design [docs/superpowers/specs/2026-10-28-runatom-integration-design.md](../specs/2026-10-28-runatom-integration-design.md) (D6: the `Work/` library ships "as templates that the UI can import").

## Decisions (user, 2026-10-30)

| # | Question | Decision |
|---|---|---|
| Q1 | Default functional | **PBE.** Library variants are ranked PBE → PW91 → other GGA → LDA; the starter and templates use exfact = 5 |
| Q2 | Elements without a library case | **Copy the radii (and the rest of the Pseudopotential tab) from a neighbouring element in the same group, with a warning** |
| Q3 | GBRV | **Follow the recommendation:** out of scope here. A follow-up plan adds an optional, user-chosen extra library folder (no download, nothing bundled) |
| Q4 | Z range | **Z ≤ 92 only.** Heavier elements are disabled in the picker and refused by `runatomDefaults` |
| Q5 | Overwrite protection | **Confirmation dialog and an *Undo defaults* button** |
| Q6 | `attach_imported_tests` open item | **Fix it first** (Task 0) |
| Q7 | Semicore shells in the starter | **Yes**, for alkali and alkaline-earth metals and early transition metals (rule in §3.3) |

## Global Constraints

- `runatomInput.py`, `runatomDefaults.py`, `periodicTableData.py` and `help.py` import no Toga. A subprocess test pins each one.
- Never edit `third_party/uspp-736/` (vendored upstream code, checked by `SOURCE_SHA256.txt`).
- Every value reaches a file through `runatomInput.case_from_fields` → `format_adat`, never as copied text (the fixed-column rule in the `runatomInput` docstring).
- `MAX_DEFAULTS_Z = 92`.
- Startup behaviour is unchanged: the window opens on `000-Ex/000-Ex-ca-sp-test`.
- Widget writes happen only on change (`_fill` already does this). Selection handlers are synchronous and `_syncing`-guarded, and async work starts as a task (the Cube-tab rule in CLAUDE.md). Toga labels do not wrap, so long warnings go to a dialog and the note line stays one line.
- `from periodicTable import PERIODIC_TABLE_LAYOUT, BLOCK_COLORS, element_block, table_symbols` keeps working (the Design tab and its tests import these names).
- Existing tests keep passing. The baseline is the 4 known failures listed in CLAUDE.md.
- Commits only when the user authorises them.

---

## 1. Relevant files, classes and functions

| File | Role today | Change |
|---|---|---|
| [runatom_input_gen.py](../../runatom_input_gen.py) | `RunatomGenerator`: `layout_main_window`, `_field_row`, `_fill`, `_load_case`, `_collect_case`, `_set_busy`, `_on_close`/`closeTopLevel`; `DEFAULT_CASE`, `NAME_FIELDS`, `TAB_FIELDS` | Element row (button, choice `Selection`, *Undo defaults*, note), `RunatomElementController`, confirmation, warning dialog, undo stack |
| [runatomInput.py](../../runatomInput.py) | `FORM_FIELDS`, `case_to_fields`/`case_from_fields`, `ae_fields`/`ps_fields`, `read_case_folder`, `attach_imported_tests`, `check_case`/`check_ae`/`check_ps`, `XC_FUNCTIONALS`, `MAX_VALENCE`, `MAX_NANG` | Task 0 fix to `attach_imported_tests`; new `xc_label(exfact)` |
| `runatomDefaults.py` (**new**, no Toga) | — | Library index, PBE-first ranking, configurations, starter, same-group template |
| `periodicTableData.py` (**new**, no Toga) | — | `PERIODIC_TABLE_LAYOUT`, `BLOCK_COLORS`, `SELECTED_ELEMENT_COLOR`, `table_symbols`, `element_block` (moved), new `group_and_period` |
| [periodicTable.py](../../periodicTable.py) | layout data + `PeriodicTableWindow(controller)` | Re-exports the data; `PeriodicTableWindow(controller, *, hint=None, marked=None, selectable=None)` |
| [runatomRunner.py](../../runatomRunner.py) | `library_dirs()` / `library_dir()` | None (library roots) |
| [molecularDesign.py](../../molecularDesign.py) | `MolecularDesignController` (the controller contract the picker expects) | None |
| [help.py](../../help.py) | `AtomicData.atomic_numbers`, `HelpGqteaWin.help_runatom` | `AtomicData.ground_state_exceptions`; help section |
| [USER_MANUAL.md](../../USER_MANUAL.md), [CLAUDE.md](../../CLAUDE.md) | runatom sections | Document the feature |
| Tests | [test_runatom_input_gen.py](../../tests/test_runatom_input_gen.py) (`_window()` stub helper), [test_runatomInput.py](../../tests/test_runatomInput.py), [test_periodicTable.py](../../tests/test_periodicTable.py), [test_runatomExe.py](../../tests/test_runatomExe.py) (skip without exe), [runatom_fixtures.py](../../tests/runatom_fixtures.py) | New `test_runatomDefaults.py`, `test_runatomStarterExe.py`, `test_periodicTableData.py`; additions to the others |

Facts from reading the code and the library (they shape the design):
- `_load_case(folder)` already does what a library pick needs: it reads the case, sets `_imported_case` (library tests keep their own inputs), fills the fields and writes a status line.
- 77 library cases (outside `000-Ex`) pass `read_case_folder` + `check_case` with no problems. `023-V-gpw-nsp-campos` is the known upstream defect, so **V has no usable library case and is a template target** (donor Nb or Ta).
- Library variants exist for 45 elements. PBE (`gpbe`) variants exist only for H, C, N, O, S, Fe, Cu, Zn. Most others are PW91 (`gpw`, exfact 4) or CA LDA.
- `periodicTable.py` imports Toga at module level, so `runatomDefaults` cannot import the layout from it. Hence `periodicTableData.py`, following the `indexRanges.py` precedent (moved out of `allBondAnalysis`, re-exported).
- runatom limits: 26 shells, **5 valence states** (`idim4 = 5`), nang ≤ 4. INPUT_AE: `ee` needs to be within a factor of 2–3; mesh `0 0 0` = program defaults; irel = 2 recommended for every element.
- In a `.adat` generation input, a projector's `keyee > 0` (and the local potential's) is a **1-based index into the valence shell list**. Copying a donor's projectors therefore needs remapping (§3.4).

## 2. The periodic-table implementation to reuse

`periodicTable.PeriodicTableWindow` is a standalone picker. Its controller contract has five parts: `register_view(view)` (must call `view.on_element_changed(current)`), `unregister_view(view)`, `select_element(symbol)`, a `periodic_table_window` attribute that `_on_close` sets to `None`, and `view.on_element_changed(symbol)`. The runatom side supplies them with a ~30-line `RunatomElementController`. `MolecularDesignController` is not refactored: the shared code would be ~15 lines, and the Design tab is large and well tested.

Changes:
- **`periodicTableData.py`**: the layout data and helpers move here unchanged, plus `group_and_period(symbol) -> Tuple[int, int]`. Main-table rows 0–6 give period = row + 1 and group = column + 1. The f rows (8, 9) give period 6 or 7 and **group 3**, the La/Ac column, so lanthanides and actinides find group-3 donors. `periodicTable.py` does `from periodicTableData import *` for the moved names.
- **`PeriodicTableWindow(controller, *, hint=None, marked=None, selectable=None)`**:
  - `hint`: replaces the Design-specific sentence; `None` gives today's text.
  - `marked`: elements drawn **bold**, plus a legend line ("Bold: a tested case is in the USPP library"); `None` gives today's look.
  - `selectable`: elements outside the set are created with `enabled=False` (here: Z > 92); `None` enables everything, as today.
  - The style computation moves into a pure `button_style_kwargs(symbol, marked) -> dict`, so it can be tested without Toga. `on_element_changed` rewrites only `background_color`, so the bold weight survives selection changes.

The picker stays a separate window. The table is about 600 × 330 px, and embedding it would push the generator's tabs off a laptop screen. The Design tab uses the same separate-window pattern.

## 3. Source and structure of the default data

### 3.1 Precedence per element (Z ≤ 92)

| Element has… | Choices offered (Selection items) | Default |
|---|---|---|
| ≥ 1 readable library case | every library variant, ranked | best-ranked variant (PBE first) |
| no library case; a library element in its group | `Template from <donor folder>`; `Computed starter (radii blank)` | template |
| no library case and none in its group (group 7: Mn, Tc, Re; group 18: He–Rn) | `Computed starter (radii blank)` | starter |

### 3.2 `runatomDefaults.py` (new, no Toga): public interface

```python
MAX_DEFAULTS_Z = 92
PBE = 5.0
FUNCTIONAL_PREFERENCE = (5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0)   # PBE, PW91, other GGA, LDA

@dataclass(frozen=True)
class LibraryVariant:
    folder: Path
    z: int
    file_symbol: str           # the Makefile's ATOM ("ti" or "Ti"), kept so file names match the library
    exfact: float
    irel: int
    xion: float
    nlcc: bool                 # ps.ifpcor == 1
    core: Tuple[int, ...]      # nnlz of the ncores core shells
    valence: Tuple[int, ...]   # nnlz of the nvales valence shells, in file order
    n_tests: int
    label: str                 # "PBE (GGA), scalar-relativistic, semicore 3s 3p, NLCC, 3 tests"

@dataclass(frozen=True)
class Suggestion:
    kind: str                  # "library" | "template" | "starter"
    label: str                 # the Selection item
    fields: Dict[str, str]     # complete form-field dict (library: case_to_fields of the case)
    folder: Optional[Path]     # library case or donor case
    note: str                  # one line for the element row
    warnings: List[str]        # shown in a dialog after the fill (template, starter, non-PBE library)

def library_elements(roots) -> Dict[str, List[Path]]: ...
def element_variants(symbol, roots) -> Tuple[List[LibraryVariant], List[str]]: ...
def rank_variants(variants) -> List[LibraryVariant]: ...
def suggestions(symbol, roots) -> Tuple[List[Suggestion], List[str]]: ...   # best first; notes = skipped cases
def ground_state_configuration(z) -> List[Tuple[int, int, float]]: ...
def split_core_valence(z, config) -> Tuple[list, list, List[str]]: ...
def eigenvalue_guess(z, config, n, l) -> float: ...
def starter_ae_input(symbol) -> ri.AEInput: ...
def starter_fields(symbol) -> Dict[str, str]: ...
def template_donors(symbol, roots) -> List[str]: ...
def template_fields(symbol, donor: LibraryVariant, donor_case: ri.CaseInputs) -> Tuple[Dict[str, str], List[str]]: ...
```

- `rank_variants` sorts by `(FUNCTIONAL_PREFERENCE.index(v.exfact + 0.0), v.irel != 2, -v.n_tests, v.folder.name)`. The form's current functional plays no part (Q1).
- `suggestions` raises `ValueError(f"{symbol}: Z = {z} > {MAX_DEFAULTS_Z}; defaults are available up to uranium")` above U. The window never sends one, because those buttons are disabled.
- Parsing reads only the chosen element's folders, plus the donor's for a template: a few ms each. Nothing is parsed when the window opens. The bold set is a folder-name scan done when the picker first opens.

### 3.3 Computed starter (All-electron tab; also the AE half of a template)

`AtomicData.ground_state_exceptions` in `help.py` holds NIST ASD ground configurations, only where they differ from Madelung order, up to U:
`Cr 3d5 4s1, Cu 3d10 4s1, Nb 4d4 5s1, Mo 4d5 5s1, Ru 4d7 5s1, Rh 4d8 5s1, Pd 4d10 5s0, Ag 4d10 5s1, La 5d1 6s2, Ce 4f1 5d1 6s2, Gd 4f7 5d1 6s2, Pt 5d9 6s1, Au 5d10 6s1, Ac 6d1 7s2, Th 6d2 7s2, Pa 5f2 6d1 7s2, U 5f3 6d1 7s2`.
They are stored as `{symbol: {(n, l): electrons}}` overrides.

**Core/valence split** (`split_core_valence`, deterministic):
1. Valence = the outermost ns and np, an open (n−1)d, an open (n−2)f. Group 12 (Zn, Cd, Hg) keeps its filled d. The p block puts its filled (n−1)d in the core (as `031-Ga-ca--vgrp`).
2. **Semicore (Q7):** add (n−1)s and (n−1)p to the valence for groups 1–2 from period 3 (Na…Ra) and groups 3–6 from period 4 (Sc–Cr, Y–Mo, La, Hf–W, Ac). Li and Be are excluded (their (n−1) shell is 1s, and the library's Li case keeps it in the core). Lanthanides and actinides other than La/Ac get none.
3. If the valence exceeds `MAX_VALENCE` (5), move the deepest valence shells to the core until it fits, and return a note.

**Eigenvalue guesses:** Slater's rules. Groups `[1s][2s2p][3s3p][3d][4s4p][4d][4f][5s5p][5d][5f][6s6p][6d][7s7p]`. Screening is 0.35 from the same group (0.30 for 1s); 0.85 / 1.00 for s,p electrons one / two or more shells below; 1.00 for every electron below a d or f group. n* = 1, 2, 3, 3.7, 4.0, 4.2 for n = 1–6 and 4.4 assumed for n = 7. `ee = -(Z_eff / n*)**2` Ry.

**AE fields:**

| Field | Value |
|---|---|
| `symbol` / `title` | periodic symbol (`Mg`) |
| `generation` | tag of the outer (non-semicore) valence shells, in the library's style: `s2` (Mg), `s2d2` (Ti), `s2d1f1` (Gd); at most 12 characters |
| `ae_print` / `ae_logd` / `ae_conv` | `0 0 0` / `0 0 0 0` / `0 0 0 0` (no plots or log derivatives; program defaults) |
| `ae_atom` | `z 0 5` (neutral, PBE) |
| `ae_mesh` | `0 0 0 2` (program-default mesh, scalar-relativistic) |
| `ae_shells` | core shells then valence shells, each in filling order, with `eigenvalue_guess` |
| `ps_partition` | `ncores nvales nang`, nang = max valence l + 1 (a template replaces nang, §3.4) |
| `tests` | blank |

Starter-only Pseudopotential tab: the controls are copied from the 000-Ex template (`1 1 3`, `2.001 -6.0 2.0 80`, `1e-10 1e-09 0.5 0`), plus partition. Every other PS field is blank, so **Check** names the first one (the existing `parse_line_values` message).

### 3.4 Same-group template (Q2)

**Donor choice** (`template_donors`): library elements with the same `group_and_period` group, sorted by `(|Δperiod|, -Z)`. The nearest period wins; on a tie the heavier wins, since 4d and 5d metals are closer to each other than to 3d. Examples: Mg → Ca; Ge → Si; Sn → Pb; Mo → W; V → Nb; Ce–Lu → La; Th–U → La; Be → Ca; Rb → K; I → Cl. Each donor's variants are tried in `rank_variants` order, and the first whose projectors remap is taken. If none remaps, the best-ranked one is used for a **radii-only** copy (below).

**What is copied** from the donor's `ps_fields`: the `ps_print`, `ps_logd`, `ps_conv`, `ps_bessel`, `ps_keyps` (keyps, ifpcor, rinner, rpcor), `ps_rinner`, `ps_rcloc` and `ps_rc` values, unchanged (the user asked for copies, not scaling), and `nang`. With a successful remap it also copies `ps_betas`, `ps_npf`, `ps_local`, `ps_qfunc` and `ps_fixes`.

**Remap** (`_remap_keyee`): rank each donor valence shell among the donor valence shells of the same l by n (0 = lowest n). A donor `keyee = k > 0` (projector or local) becomes the 1-based index of the target valence shell with the same l and the same rank. `keyee = 0` (energy given by `eeread`/`eloc`) is copied as is. If any `keyee > 0` has no counterpart, the remap fails: only the radii fields, controls, bessel and nang are copied, and projectors, local, qfunc and fixes stay blank.

**Functional:** the template is PBE (exfact 5). If the donor variant is not GGA and its `ifqopt ≠ 3`, `ps_qfunc` is set to `3 <nqf> <qtryc>` when the donor's `nqf > 3`, and left blank otherwise. Both cases get a warning.

**Check after building:** `ri.check_case(ri.case_from_fields(fields))` runs on the template when every field is filled. Each problem it lists is added to the warnings. A template is never silently "fixed" beyond the remap and the ifqopt rule.

**Warning** (dialog title `Values copied from <Donor>`; note line `Template from <folder> — review before use`):

> Mg has no case in the USPP library. The All-electron tab was computed from its ground-state configuration (PBE, scalar-relativistic, semicore 2s 2p). The Pseudopotential tab — rc, rcloc, rinner, rpcor, the projectors and the local potential — was copied **unchanged** from Ca (020-Ca-gpw-nsp-campos, PW91, same group, one period below). These radii were tuned for Ca, not Mg: they are a starting point, not a tested pseudopotential. Review them, add transferability tests, and compare against all-electron results before production use.

(Plus remap-failure, functional and `check_case` lines when they apply.)

## 4. UI and data-flow changes

### 4.1 Layout (`layout_main_window`, top of the `names` box)

```
Element  [Choose element...]  Mg  [ Template from 020-Ca-gpw-nsp-campos ▾ ]  [Undo defaults]
         Template from 020-Ca-gpw-nsp-campos — review before use. Every field stays editable.
```

- `element_button` → `on_choose_element`.
- `element_label`: the element.
- `choice_selection` (`toga.Selection`, `on_change=self.on_choice_change`): items = the labels of `suggestions(...)`.
- `undo_button` (`Undo defaults`, `enabled=False` until there is something to undo).
- `element_note`: one line, hint colour.

`_set_busy` adds `element_button`, `choice_selection` and `undo_button` (undo stays disabled after a run when its stack is empty). `_on_close`/`closeTopLevel` call `element_controller.close_picker()`.

### 4.2 Controller (in `runatom_input_gen.py`)

```python
class RunatomElementController:
    """The periodic table's controller (see periodicTable.PeriodicTableWindow)."""

    def __init__(self, on_select):
        self.selected_element = None
        self.views = []
        self.periodic_table_window = None
        self._on_select = on_select              # RunatomGenerator.element_selected

    def register_view(self, view):
        if view not in self.views:
            self.views.append(view)
            if self.selected_element:
                view.on_element_changed(self.selected_element)

    def unregister_view(self, view):
        if view in self.views:
            self.views.remove(view)

    def show_element(self, symbol):              # follow a change made elsewhere (startup, import, undo)
        self.selected_element = symbol
        for view in list(self.views):
            view.on_element_changed(symbol)

    def select_element(self, symbol):            # a click in the periodic table
        self._on_select(symbol)

    def open_picker(self, *, hint, marked, selectable):
        if self.periodic_table_window is not None and self.periodic_table_window.window is not None:
            self.periodic_table_window.window.show()
        else:
            self.periodic_table_window = PeriodicTableWindow(self, hint=hint, marked=marked,
                                                             selectable=selectable)

    def close_picker(self):
        if self.periodic_table_window is not None:
            self.periodic_table_window.close()
            self.periodic_table_window = None
```

`select_element` does not highlight the button itself. The generator calls `show_element` only after the fill was confirmed, so a declined replacement leaves the highlight on the element still in the form.

### 4.3 Data flow

```
click "Mg" ─► controller.select_element ─► element_selected("Mg")  (sync) ─► _start(_apply_element("Mg"))
   _apply_element:
     options, notes = rd.suggestions("Mg", rr.library_dirs())       (ValueError → error dialog)
     if self._form_edited() and not await self._confirm_replace("Mg"): return
     self._push_undo()
     self._show_choices(options, options[0])                        (under _syncing)
     await self._apply_suggestion(options[0])

on_choice_change (sync, _syncing-guarded) ─► _start(_apply_choice(label))
     same confirmation + _push_undo, then _apply_suggestion(that option)

_apply_suggestion(s):
     library  → self._load_case(s.folder)                           (existing path, sets _imported_case)
     template / starter → self._imported_case = None; self._fill(s.fields)
     self._baseline = self._current_fields()
     controller.show_element(symbol); element_note = s.note (+ skipped-case notes, shortened)
     if s.warnings: await self._warn(f"Values copied from …" / "Defaults for …", "\n\n".join(s.warnings))

Undo defaults ─► on_undo (sync) ─► restore the last snapshot
```

- **Snapshot** (`_DefaultsSnapshot`, a dataclass): `fields`, `imported_case`, `element`, `choices` (labels), `choice`, `options`, `note`, `baseline`. `_push_undo` captures the state **before** a defaults fill. `UNDO_DEPTH = 10`, oldest dropped. Undo restores the fields via `_fill`, `_imported_case`, the Selection items and value (under `_syncing`), the element highlight and the note, and pops the snapshot. The button is enabled only while the stack is non-empty and no run is going. Only defaults fills push snapshots. Imports and hand edits do not, because the button is named *Undo defaults*.
- **Dirty check**: `_baseline` is the field dict after the last programmatic fill (startup, defaults, either import, undo). `_form_edited()` is `_current_fields() != _baseline`. The `QuestionDialog` ("Replace the current inputs with the suggested values for Mg? You can undo this with Undo defaults.") appears only then.
- **Imports and startup**: after `_load_case` (startup 000-Ex, Import library case) the element row follows: `show_element(symbol_of_z)`, with items = that element's suggestions plus `Imported: <folder>` selected when the folder is not one of them. Import .adat updates `_baseline` only.
- **Editing stays free**: no widget is made read-only and nothing re-fills on edit. The Symbol field can still be typed.
- **Library tests**: `_load_case` sets `_imported_case`; templates and starters set `None`. Together with Task 0, a library test keeps its own inputs only while the generation it was made for is unchanged.

## 5. Validation, fallback and error handling

| Situation | Behaviour |
|---|---|
| Library element, PBE variant exists | PBE variant loaded; no dialog; note "From <folder>. Every field stays editable." |
| Library element, no PBE variant (e.g. Ti, Ca, Au) | best GGA, else LDA, loaded; **warning**: "No PBE case for Ti in the library; loaded <label>. Changing exfact to 5 alone is not enough: GGA needs ifqopt 3, and the radii were tuned for this functional." |
| A variant unreadable (023-V) or z ≠ folder number | skipped; the reader's message is listed in the warning dialog when it changes the outcome (V becomes a template) |
| No library case, group has a donor | template (§3.4) + warning dialog; the starter is also offered |
| Remap fails | radii-only template; warning says projectors/local/qfunc were left blank |
| Group 7 or 18 | starter only; warning: "No element of group 7 is in the USPP library, so no radii could be copied. Enter rc, rcloc, rinner, the projectors and the local potential." |
| Z > 92 | button disabled in the picker; `suggestions` raises `ValueError` (defensive) → error dialog |
| Starter valence trimmed to 5 shells | note in the warning |
| Library folder missing (frozen build without data) | every element gets the starter (no donors); warning says `LIBRARY_MISSING`; picker marks nothing |
| User edited the form | confirmation dialog; No leaves everything (including the highlight) unchanged |
| Run in progress | element button, choice selection and undo disabled |
| Check / Save / Run | unchanged: blank or malformed fields are named by `case_from_fields`, and runatom's rules by `check_case`. Defaults never bypass them |
| `OSError`/`ValueError`/`KeyError` while reading a case | caught in `element_variants` → note; nothing raises out of an async handler (the window's existing `except Exception` pattern in `on_run` is mirrored) |

`ri.xc_label(exfact)` returns `XC_FUNCTIONALS.get(exfact + 0.0, f"exfact {exfact:g}")`, so the library's `-0.` prints as CA (LDA), never "-0".

## 6. Testing requirements

Headless throughout. `_window()` in `test_runatom_input_gen.py` gains the stub attributes `element_button`, `element_label`, `choice_selection` (records `items`, `value`), `undo_button`, `element_note`, `element_controller` (records `show_element`, `close_picker`), `_baseline`, `_syncing`, `_undo_stack`, plus a scripted `QuestionDialog` answer and a recorded `_warn`.

**Task 0, `test_runatomInput.py`:**
- An imported Ca case whose form changes only exfact (or irel, rmax, aasf, bbsf, z, ncores, nvales, nang) gets derived tests (`ae_input is None`); unchanged, the tests are the library's (the existing test still passes).

**`test_periodicTableData.py`** (new) **and `test_periodicTable.py`:**
- `periodicTableData` imports no Toga (subprocess); the old names still import from `periodicTable`.
- `group_and_period`: H (1, 1), Mg (2, 3), Ti (4, 4), Sn (14, 5), U (3, 7), Ce (3, 6), Rn (18, 6).
- Default window: today's hint, no bold, every button enabled (regression pin via `button_style_kwargs`).
- `marked={"Ti"}` → Ti bold, C not; bold survives `on_element_changed`. `selectable` without Np → Np disabled.

**`test_runatomDefaults.py`** (new):
- Imports no Toga (subprocess).
- `library_elements` = the 45 symbols (H, Li, B, C, N, O, F, Na, Al, Si, P, S, Cl, K, Ca, Sc, Ti, V, Cr, Fe, Ni, Cu, Zn, Ga, As, Se, Sr, Y, Zr, Nb, Ru, Rh, Ag, Cd, Ba, La, Ta, W, Os, Ir, Pt, Au, Tl, Pb, Bi); `000-Ex` and `History/` excluded.
- Ranking: H → `001-H-gpbe--bm`; Cu → `029-Cu-gpbe-n-bm`; Ti (LDA only) → `022-Ti-ca-sp-vgrp`, with the no-PBE warning; Ca → `020-Ca-gpw-nsp-campos` (PW91 before LDA).
- A library suggestion's fields are exactly `ri.case_to_fields(ri.read_case_folder(folder))`.
- V: the library variant is skipped with a note naming `V_ps.adat` line 2; the default is a template whose donor is Nb.
- `ground_state_configuration(z)` sums to z for z = 1..92; the 17 exceptions are pinned; Fe `[Ar] 3d6 4s2`.
- Starter for z = 1..92: `ri.check_ae` → `[]`; `ncores + nvales == len(shells)`; `1 ≤ nvales ≤ 5`; `nang ≤ 4`; core before valence; Ti valence `{300, 310, 400, 320}`, K `{300, 310, 400}`, Ga has no 3d in the valence; exfact 5.
- `eigenvalue_guess` is within a factor of 3 of the library authors' `ee` for every valence shell of every library AE input.
- Donors: Mg → Ca, Ge → Si, Sn → Pb, Mo → W, Ce → La, U → La, I → Cl; Mn, Ne → none.
- Remap: Mg from Ca keeps every projector's l and points `keyee > 0` at Mg's matching valence shell; a constructed donor with an f projector bound to a 4f shell that the target lacks → radii-only, projectors blank, warning present.
- Every template for Z ≤ 92 whose remap succeeds either passes `check_case` or lists every `check_case` problem in its warnings (never silent).
- Z = 93 → `ValueError`.

**`test_runatomStarterExe.py`** (new; `skipif(rr.find_runatom() is None)`):
- The starter AE stage (`rd.starter_ae_input` → `format_adat` → `rr.run_stage`) completes for **every Z from 1 to 92**. Failures are fixed in `eigenvalue_guess`/`split_core_valence`, never by skipping an element.
- Dev-time report, not a pass/fail test: `benchmarks/report_runatom_templates.py` runs AE + generation for every template element. Its table (finished / stopped at stage + message) is quoted in USER_MANUAL, so users know which templates generate as copied.

**`test_runatom_input_gen.py`** (additions):
- Picking H fills the PBE case, sets `_imported_case`, Selection items in ranked order, no warning.
- Picking Ti records the no-PBE warning.
- Picking Mg fills template fields (`ae_atom == "12  0  5"`, `ps_rc` = Ca's), `_imported_case is None`, the warning title names Ca; choosing *Computed starter* blanks `ps_rc`.
- Picking Mn: starter only, group-7 warning.
- Setting items under `_syncing` causes no second fill (count `_fill`/`_load_case` calls).
- Edited form + pick, answer No → fields, element highlight, undo stack unchanged; Yes → replaced and one snapshot pushed. An unedited form → no dialog.
- Undo restores fields, `_imported_case`, Selection and highlight; two picks then two undos return to the startup state; the button is disabled when the stack is empty; depth is capped at 10.
- After a fill, editing a field changes `_collect_case()`.
- Library missing → Mg gets the starter with the library warning.
- `_set_busy(True)` disables the three new controls; `_on_close` closes the picker.
- Startup shows element Ti with `Imported: 000-Ex-ca-sp-test`.
- `test_every_form_field_has_a_widget` is unchanged.

**Regression:** `python -m pytest tests/ -q` against the 4-failure baseline.

**Live check (manual, recorded in Task 6):** bold marks and disabled Z > 92 buttons render on WinForms; Mg/Ti/V/Mn/H picks; warning dialog; Undo; no window growth or flicker; Run on H (library) and Mg (template).

## 7. Risks, assumptions and remaining open points

**Risks**
- *Copied radii mistaken for a tested pseudopotential.* Mitigated by the warning dialog on every template fill, the note line, a distinct Selection label, `check_case` results in the warning, and the template report in the manual. Residual: a template that runs to completion can still be poor. Only the user's tests reveal that.
- *Unscaled radii across periods* (e.g. Sn from Pb, Ge from Si). This is the user's choice, and the warning says so. Covalent-radius scaling could be added later as an explicit option.
- *Remap brittleness.* Different valence shapes between donor and target are handled by the radii-only fallback, never by guessing.
- *Overwriting work.* Confirmation dialog plus a 10-deep *Undo defaults*.
- *Selection re-entrancy.* `_syncing` + synchronous handlers.
- *Bold / disabled buttons on WinForms.* If bold does not render, fall back to a marker suffix (`Ti•`) with the same `marked` API.
- *Slater guesses for 4f/5f elements.* Exposed by the Z 1–92 exe test; the fallback is the nearest library element's guesses (INPUT_AE's other suggestion).
- *Moving the layout data.* Re-exports keep the Design tab's imports; its tests are in the regression run.

**Assumptions**
- The bundled `Work/` library is trusted and suggested unmodified (spec D6).
- NIST ASD ground configurations are the reference for the neutral-atom starter.
- A template is PBE even when its donor is not. The warning states the donor's functional.

**Remaining open points** (defaults chosen; say if you want otherwise)
- Groups 7 and 18 have no library donor: they get the starter only (blank radii). Alternative: borrow from a same-period neighbour group (Mn from Cr/Fe).
- The template's `tests` field stays blank. Alternative: add a +1 ion test automatically.

---

## Tasks

### Task 0: Library tests attach only to the generation they were made for

**Files:** Modify `runatomInput.py` (`attach_imported_tests`); extend `tests/test_runatomInput.py`.

- [ ] Write the failing test: import `020-Ca-gpw-nsp-campos`, change `ae_atom` exfact 4 → 5 in the fields, `attach_imported_tests(case_from_fields(fields), imported)` → every test has `ae_input is None and test_input is None`. Repeat for irel, rmax, ncores/nvales, nang. Unchanged fields → the library tests (existing test).
- [ ] Run `python -m pytest tests/test_runatomInput.py -q`: FAIL.
- [ ] Implement:

```python
def _same_generation(a: CaseInputs, b: CaseInputs) -> bool:
    """The all-electron and partition values a library test's own inputs were made for."""
    keys = lambda c: (c.ae.z, c.ae.exfact + 0.0, c.ae.irel, c.ae.rmax, c.ae.aasf, c.ae.bbsf,
                      c.ps.ncores, c.ps.nvales, c.ps.nang)
    return _same_shells(a.ae, b.ae) and keys(a) == keys(b)
```

  and use `if imported is None or not _same_generation(case, imported): return case`.
- [ ] Run `tests/test_runatomInput.py` and `tests/test_runatom_input_gen.py`: PASS.
- [ ] Commit `fix(uspp): library tests keep their own inputs only for an unchanged generation`, and update the memory file (open item 1 closed).

### Task 1: Toga-free periodic-table data and picker options

**Files:** Create `periodicTableData.py`, `tests/test_periodicTableData.py`; modify `periodicTable.py`, `tests/test_periodicTable.py`.

**Produces:** `periodicTableData.{PERIODIC_TABLE_LAYOUT, BLOCK_COLORS, SELECTED_ELEMENT_COLOR, table_symbols, element_block, group_and_period}`; `periodicTable.button_style_kwargs(symbol, marked)`; `PeriodicTableWindow(controller, *, hint=None, marked=None, selectable=None)`; `periodicTable.DEFAULT_HINT`, `MARKED_LEGEND`.

- [ ] Write the failing tests (§6).
- [ ] Move the data verbatim; add `group_and_period`:

```python
def group_and_period(symbol: str) -> Tuple[int, int]:
    """(IUPAC group 1-18, period); lanthanides and actinides count as group 3."""
    for row_index, row in enumerate(PERIODIC_TABLE_LAYOUT):
        if symbol in row:
            if row_index >= 8:
                return 3, row_index - 2          # f rows: periods 6 and 7
            return row.index(symbol) + 1, row_index + 1
    raise KeyError(symbol)
```

- [ ] Add the three keyword arguments and `button_style_kwargs`; disabled buttons get `enabled=False`.
- [ ] Run `tests/test_periodicTable*.py` and `tests/test_molecularDesign.py`: PASS.
- [ ] Commit `refactor(periodic-table): Toga-free layout data; optional hint, marks and selectable set`.

### Task 2: Library index and PBE-first ranking

**Files:** Create `runatomDefaults.py`; add `ri.xc_label`; create `tests/test_runatomDefaults.py`.

**Produces:** `MAX_DEFAULTS_Z`, `FUNCTIONAL_PREFERENCE`, `LibraryVariant`, `library_elements`, `element_variants`, `rank_variants`, `ri.xc_label`.

- [ ] Write the failing tests: no-toga subprocess, the 45 symbols, the ranking cases, fields unchanged, V skipped.
- [ ] Implement:

```python
ELEMENT_FOLDER = re.compile(r"^(\d{3})-([A-Z][a-z]?)$")

def library_elements(roots):
    found: Dict[str, List[Path]] = {}
    for root in roots:
        if not Path(root).is_dir():
            continue
        for element_dir in sorted(Path(root).iterdir()):
            m = ELEMENT_FOLDER.match(element_dir.name)
            if not m or int(m.group(1)) != AtomicData.atomic_numbers.get(m.group(2)):
                continue
            cases = sorted(p.parent for p in element_dir.glob("*/Makefile"))
            if cases:
                found.setdefault(m.group(2), []).extend(cases)
    return found

def element_variants(symbol, roots):
    variants, notes = [], []
    for folder in library_elements(roots).get(symbol, []):
        try:
            case = ri.read_case_folder(folder)
        except (OSError, ValueError, KeyError) as exc:
            notes.append(f"{folder.name} skipped: {exc}")
            continue
        if round(case.ae.z) != AtomicData.atomic_numbers[symbol]:
            notes.append(f"{folder.name} skipped: z = {case.ae.z:g} is not {symbol}")
            continue
        variants.append(_describe(folder, case))
    return variants, notes

def rank_variants(variants):
    def rank(v):
        xc = v.exfact + 0.0
        order = FUNCTIONAL_PREFERENCE.index(xc) if xc in FUNCTIONAL_PREFERENCE else len(FUNCTIONAL_PREFERENCE)
        return (order, v.irel != 2, -v.n_tests, v.folder.name)
    return sorted(variants, key=rank)
```

  `_describe` builds the label from `ri.xc_label`, "scalar-relativistic" (irel 2), "semicore …" (valence shells whose n is below the highest valence n, s/p only), "NLCC" (`ifpcor == 1`), `f"{xion:+g} ion"` when ≠ 0, and `f"{n} test(s)"`.
- [ ] Run: PASS. Commit `feat(uspp): index library cases by element, PBE first`.

### Task 3: Computed starter (Z ≤ 92)

**Files:** `help.py` (`AtomicData.ground_state_exceptions`); `runatomDefaults.py`; tests.

**Produces:** `ground_state_configuration`, `split_core_valence`, `eigenvalue_guess`, `starter_ae_input`, `starter_fields`, `STARTER_PS_CONTROLS`.

- [ ] Write the failing tests: sums for z 1..92, the 17 exceptions, `check_ae` clean, the semicore cases (Ti, K, Na, Sc, Mo, La) and their absence for Li, Be, Ga, Ce; the factor-3 `ee` check.
- [ ] Implement per §3.3 (Madelung order `sorted(((n, l) for n in range(1, 8) for l in range(min(n, 4))), key=lambda nl: (sum(nl), nl[0]))`, fill 2(2l+1), apply exceptions, refill to z; Slater's rules; split; fields).
- [ ] Run: PASS. Record any library-valence disagreements in the task notes.
- [ ] Commit `feat(uspp): computed starter inputs up to uranium`.

### Task 4: Same-group template and `suggestions`

**Files:** `runatomDefaults.py`; tests.

**Consumes:** `periodicTableData.group_and_period`, Tasks 2–3.
**Produces:** `template_donors`, `_remap_keyee`, `template_fields`, `Suggestion`, `suggestions`.

- [ ] Write the failing tests: donors, remap success (Mg ← Ca), remap failure, the ifqopt rule, `check_case` problems surfaced, `suggestions` order for H / Ti / Mg / V / Mn, Z = 93 refused.
- [ ] Implement:

```python
def _remap_keyee(donor_valence, target_valence, keyee):
    """Donor valence index (1-based) -> target index of the same l and n-rank; 0 stays 0."""
    if keyee <= 0:
        return keyee
    n, l = divmod(donor_valence[keyee - 1] // 10, 10)
    rank = sorted(s for s in donor_valence if (s // 10) % 10 == l).index(donor_valence[keyee - 1])
    same_l = sorted(s for s in target_valence if (s // 10) % 10 == l)
    if rank >= len(same_l):
        raise LookupError(f"no {('s', 'p', 'd', 'f')[l]} valence shell of rank {rank} in the target")
    return target_valence.index(same_l[rank]) + 1
```

  `template_fields` starts from `starter_fields(symbol)`, overlays the copied donor `ps_fields` (§3.4), rewrites the `ps_betas`/`ps_local` keyee column through `_remap_keyee` (radii-only on `LookupError`), sets nang in `ps_partition`, applies the ifqopt rule, then collects the warnings (the §3.4 text with donor, folder, period offset and functional, plus remap/ifqopt/`check_case` lines).
- [ ] Run: PASS. Commit `feat(uspp): same-group templates for elements without a library case`.

### Task 5: Validate the starter with the real program

**Files:** Create `tests/test_runatomStarterExe.py`, `benchmarks/report_runatom_templates.py`; adjust the Task 3 rules if needed.

- [ ] Write the exe test for Z 1..92 (AE stage in `tmp_path`, completion marker + exit 0).
- [ ] Run it, fix failures in the guess/split rules, and run it again until all 92 complete.
- [ ] Write and run the report script; save its table for Task 7.
- [ ] Commit `test(uspp): the starter all-electron input converges for Z 1-92`.

### Task 6: Element row, confirmation, warnings and Undo in the window

**Files:** `runatom_input_gen.py`; `tests/test_runatom_input_gen.py`.

**Produces:** `RunatomElementController`, `element_selected`, `on_choose_element`, `on_choice_change`, `on_undo`, `_apply_suggestion`, `_confirm_replace`, `_warn`, `_push_undo`, `_DefaultsSnapshot`, `_form_edited`, `_current_fields`, `UNDO_DEPTH = 10`, `PICKER_HINT`.

- [ ] Extend `_window()`; write the failing tests (§6).
- [ ] Implement §4. `_current_fields()` is the comprehension now inline in `_collect_case`; reuse it there. Element from z for imports: `next(s for s, n in AtomicData.atomic_numbers.items() if n == round(z))`. `_warn(title, message)` = `main_window.dialog(toga.InfoDialog(...))` when not closed, recorded in tests like `_error`.
- [ ] Run the file, then the whole suite against the baseline.
- [ ] Live check (§6); record the results.
- [ ] Commit `feat(uspp): choose an element to fill the inputs with suggested values`.

### Task 7: Documentation

**Files:** `help.py` (`help_runatom`), `USER_MANUAL.md` (runatom: *Choosing an element*, with precedence, PBE preference, templates and their warning, the template report table, Undo, Z ≤ 92), `CLAUDE.md` (runatom bullet: `runatomDefaults`, `periodicTableData`, picker options, Task 0 rule).

- [ ] Write the texts; run `python -m pytest tests/ -q`.
- [ ] Commit `docs(uspp): element defaults`.

## Implementation notes (2026-10-30): where the implementation departs from rev. 2

Each change was forced by running the shipped `runatom.exe`:

1. **Eigenvalue guesses.** Slater's rules missed converged eigenvalues by 2× to 100× (Nd 4f: −12 Ry guessed against −0.13 Ry), and 39 of 92 starters failed. So the planned check "within a factor of 3 of the library's `ee`" was dropped: the library's own guesses are often 50× off and still converge. The guesses now come from the generated `runatomStarterData.py` (built by `tools/make_runatom_starter_data.py`, bootstrapping through Z, as INPUT_AE suggests). Slater remains only as the fallback for shells not in the table.
2. **`rlogd`.** `rlogd = 0` makes `rinit` stop (klogd < 10) even with ilogd 0. The starter uses rlogd 2.0, and `check_ae`/`check_ps` now report it (`log_derivative_index`).
3. **Mesh and damping.** The starter uses rmax 150 and damp 0.2. Twelve atoms (Co, Ag, Ce–Tb, Pa, U) need rmax 400 and damp 0.1–0.2, recorded per atom in `STARTER_SETTINGS`.
4. **Dy, Ho, Er** do not converge as neutral atoms with any setting tried (the 4f state). The starter uses the +1 ion (`STARTER_XION`), with a warning. The exe test runs all 92 starters.
5. **Template remap.** The rule "same l and rank" was replaced by **valence position**: `keyee` takes a state's eigenvalue whatever its l (Nb's d projector uses 5p). Donor valence shells the target holds no electron in are **dropped**, because Ca's empty 3d made Mg's all-electron stage fail; their references become keyee 0 with the donor's `ee`. So a template's all-electron input always has exactly the starter's shells.
6. **`psqf2`.** ifqopt 3 needs nqf ≥ 6 (≥ 7 with a constraint), so the template raises nqf. `check_ps` now mirrors this rule; it was missing.
7. **Live check.** The window size is unchanged (1025 × 700 logical at 125 %). Bold and disabled picker buttons render, and the Mg template, starter choice, Undo and picker close all work.

## Self-review

- Q1 → `FUNCTIONAL_PREFERENCE`, ranking tests, starter exfact 5. Q2 → §3.4, Task 4, warning dialog. Q3 → out of scope, recorded. Q4 → `MAX_DEFAULTS_Z`, `selectable`, refusal test. Q5 → §4.3 confirmation + undo stack, Task 6 tests. Q6 → Task 0. Q7 → §3.3 semicore rule, Task 3 tests.
- Names are consistent across tasks: `suggestions`, `Suggestion`, `rank_variants`, `template_fields`, `_remap_keyee`, `starter_ae_input`, `group_and_period`, `button_style_kwargs`, `_DefaultsSnapshot`, `UNDO_DEPTH`.
