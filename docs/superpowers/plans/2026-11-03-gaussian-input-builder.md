# Gaussian 16 Input Builder Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A GaussView-6.1-style Gaussian 16 input builder on the Inputs tab (core panels, linked steps, reopening `.gjf` files), with the Fukui linked job folded in as a Job Type.

**Architecture:** Three new modules in a one-way chain: `gaussianKeywords` (data tables), `gaussianInput` (Toga-free job model, validation, route building, input text, reopening, files), and `gaussianInputBuilder` (the Toga window). Every control's route form is defined once in a table and read back through the same table. The Fukui Job Type builds a `fukuiInput.FukuiJobSpec` and delegates to the existing `fukuiInput.FukuiInputBuilder`, so its output is unchanged.

**Tech Stack:** Python 3.13, Toga 0.5 (WinForms), pytest. No new dependencies.

**Spec:** [docs/superpowers/specs/2026-11-03-gaussian-input-builder-design.md](../specs/2026-11-03-gaussian-input-builder-design.md)

## Global Constraints

- Run every command from `venv/src/` with the venv active; tests: `python -m pytest tests/ -q`.
- Baseline: four tests fail on `main` before any change (`test_meanResidenceTime.py::TestAdvancedOptionsExposed::test_every_run_parameter_is_supplied_by_the_ui` and three table-styling tests in `test_plotter.py`). Compare against that, do not "fix" them.
- `gaussianKeywords.py` and `gaussianInput.py` must not import Toga (subprocess test, Task 2). Only `gaussianInputBuilder.py` imports Toga.
- `gaussianInputBuilder.py` must not import `cubeTab` (it would load OpenGL for an input builder): it defines its own 6-line write-on-change `_assign` with the same rule.
- Input files: ASCII, LF line endings (`newline="\n"`), relative `%Chk`/`%OldChk` names, for Gaussian 16 (Rev. C.01) under WSL or a Linux cluster.
- Multiplicity is required and never assumed: a blank field is `None` and becomes an error naming the allowed values.
- `prepare` collects every problem into one `ValueError` whose lines start with `- `; problems of a step start with `- Step n: `.
- Window handlers are synchronous and `_syncing`-guarded; widget writes go through `_assign` (write only on change).
- Atom numbers shown to the user are 1-based.
- `fukuiInput.FukuiInputBuilder` keeps producing byte-identical files. The only edit to `fukuiInput.py` is Task 5's extraction of the formchk script body into `formchk_script()`, pinned byte-for-byte by a test written before the change. (This refactor is not in the spec; it keeps one copy of the script instead of two.)
- Two deliberate refinements of the spec, applied by this plan: (1) routes are written in `fukuiInput.build_route` order (print level, model, then keywords), not GaussView's job-first order; (2) `Guess=TCheck` does not require a checkpoint (it falls back to a new guess by definition); `Guess=Read`, `Geom=Check/AllCheck` and `ReadFC` do.
- Commit after each task with a message ending in:
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`

## File Structure

| File | Responsibility |
|---|---|
| `gaussianKeywords.py` (new) | Data tables: labels shown in the window → Gaussian keywords/options; the panel that owns each route keyword. |
| `gaussianInput.py` (new) | `RouteKeyword` algebra (parse/format/merge), `GaussianStepSpec`/`GaussianJobSpec`, route building, ModRedundant lines, `GaussianInputBuilder` (prepare/input_text/script_text/write_files), Fukui delegation, reopening (`read_input_text`/`read_input_file`), `read_geometry_file`, `wsl_path`, live preview texts. |
| `gaussianInputBuilder.py` (new) | `GaussianInputBuilderUI`: panels, `STEP_WIDGETS` field table, step bar, dialogs. |
| `fukuiInput.py` (modify) | Extract `formchk_script()`; `FukuiInputBuilder.script_text` calls it. |
| `fukuiInputBuilder.py` (delete) | Replaced by the new window. |
| `toolLoader.py` (modify) | Inputs-tab entry and `_frozen_build_imports`. |
| `help.py` (modify) | `help_gaussian_input` (new), `help_fukui_input` rewritten for the Job Type, `help_fukui` pointer. |
| `tests/test_gaussianKeywords.py`, `tests/test_gaussianInput.py`, `tests/test_gaussianInputReopen.py`, `tests/test_gaussianInputBuilder.py` (new) | Tests. |
| `tests/test_fukuiInputBuilder.py` (delete) | Migrated into `tests/test_gaussianInputBuilder.py`. |
| `tests/data/gaussian_input/phenol_formchk.sh` (new) | Byte-for-byte pin of the Fukui formchk script. |
| `USER_MANUAL.md`, `readme.md`, `docs/index.html`, `CLAUDE.md` (modify) | Docs. |

---

### Task 1: Keyword tables (`gaussianKeywords.py`)

**Files:**
- Create: `gaussianKeywords.py`
- Test: `tests/test_gaussianKeywords.py`

**Interfaces:**
- Consumes: `fukuiInput.PRINT_LEVELS`, `fukuiInput.DEFAULT_PRINT_LEVEL`.
- Produces (names used by every later task): `OTHER`, `NO_CHOICE`; job types `JOB_ENERGY`, `JOB_OPT`, `JOB_FREQ`, `JOB_OPT_FREQ`, `JOB_IRC`, `JOB_SCAN`, `JOB_STABILITY`, `JOB_NMR`, `JOB_FUKUI`, `JOB_TYPES`, `OPTIMIZING_JOBS`, `FREQUENCY_JOBS`; `OPT_MINIMUM`, `OPT_TS`, `OPT_TARGETS`; `FC_DEFAULT`, `FC_CALCFC`, `FORCE_CONSTANTS`, `FUKUI_FORCE_CONSTANTS`, `READ_FORCE_CONSTANTS`; `IRC_DIRECTIONS`; `NMR_METHODS`; `FREQ_FLAG_FIELDS`; states `STATE_GROUND`, `STATE_TD`, `STATE_TDA`, `STATE_CIS`, `STATE_ZINDO`, `STATES`, `STATE_KEYWORDS`, `MODEL_STATES`, `STATE_SPINS`, `DEFAULT_N_STATES`, `COEFFICIENTS_IOP`; `SPIN_PREFIXES`; families `FAMILY_HF`, `FAMILY_DFT`, `FAMILY_SEMI`, `FAMILY_MP2`, `FAMILY_MP4`, `FAMILY_CCSD`, `FAMILY_CCSDT`, `FAMILY_COMPOUND`, `FAMILIES`, `FAMILY_KEYWORDS`, `FAMILY_METHODS`, `NO_BASIS_FAMILIES`, `POST_SCF_FAMILIES`, `DEFAULT_METHOD`; bases `POPLE_BASES`, `EXTENSIBLE_POPLE`, `OTHER_BASES`, `BASES`, `GENERAL_BASES`, `DIFFUSE`, `POLARIZATION`; `DISPERSION`; geometry `GEOM_COORDINATES`, `GEOM_ALLCHECK`, `GEOM_CHECK`, `GEOMETRY_SOURCES`; `GUESS_METHODS`, `READ_GUESSES`, `GUESS_FLAG_FIELDS`; `POP_LEVELS`, `POP_FLAG_FIELDS`; `SOLVATION_NONE`, `SOLVATION_MODELS`, `SOLVATION_ALIASES`, `SOLVENTS`; `MODRED_TYPES`, `MODRED_ACTIONS`; `PANEL_OF_KEYWORD`; `PRINT_LEVELS`, `DEFAULT_PRINT_LEVEL`; `reverse(table) -> Dict[str, str]` (upper-case keyword → label, `None` values skipped).

- [ ] **Step 1: Write the failing tests**

`tests/test_gaussianKeywords.py`:

```python
"""The data tables of the Gaussian 16 input builder (gaussianKeywords, no Toga)."""
import pytest

import gaussianKeywords as gk

LABEL_TABLES = {
    "OPT_TARGETS": gk.OPT_TARGETS, "FORCE_CONSTANTS": gk.FORCE_CONSTANTS,
    "IRC_DIRECTIONS": gk.IRC_DIRECTIONS, "NMR_METHODS": gk.NMR_METHODS, "STATE_SPINS": gk.STATE_SPINS,
    "SPIN_PREFIXES": gk.SPIN_PREFIXES, "DIFFUSE": gk.DIFFUSE, "POLARIZATION": gk.POLARIZATION,
    "DISPERSION": gk.DISPERSION, "GEOMETRY_SOURCES": gk.GEOMETRY_SOURCES, "GUESS_METHODS": gk.GUESS_METHODS,
    "POP_LEVELS": gk.POP_LEVELS, "SOLVATION_MODELS": gk.SOLVATION_MODELS,
    **{f"FAMILY_METHODS[{family}]": table for family, table in gk.FAMILY_METHODS.items()},
}


@pytest.mark.parametrize("name", sorted(LABEL_TABLES))
def test_route_values_are_unique_so_a_route_reads_back_to_one_label(name):
    values = [v.upper() for v in LABEL_TABLES[name].values() if v]
    assert len(values) == len(set(values)), name


@pytest.mark.parametrize("name", sorted(LABEL_TABLES))
def test_reverse_maps_every_route_value_to_its_label(name):
    table = LABEL_TABLES[name]
    back = gk.reverse(table)
    for label, value in table.items():
        if value:
            assert back[value.upper()] == label


def test_choice_lists_have_unique_entries():
    for values in (gk.JOB_TYPES, gk.STATES, gk.FAMILIES, gk.BASES, gk.SOLVENTS, gk.PRINT_LEVELS):
        assert len(values) == len(set(values))


def test_every_family_has_a_route_form():
    for family in gk.FAMILIES:
        assert (family in gk.FAMILY_KEYWORDS) != (family in gk.FAMILY_METHODS), family


def test_default_method_of_each_choice_family_is_one_of_its_methods():
    for family, method in gk.DEFAULT_METHOD.items():
        assert method in gk.FAMILY_METHODS[family]


def test_keywords_gaussian_spells_differently_from_the_label():
    assert gk.FAMILY_METHODS[gk.FAMILY_DFT]["M06-2X"] == "M062X"
    assert gk.FAMILY_METHODS[gk.FAMILY_DFT]["wB97X-D"] == "wB97XD"
    assert gk.FAMILY_METHODS[gk.FAMILY_DFT]["PBE0 (PBE1PBE)"] == "PBE1PBE"


def test_solvent_names_hold_no_comma_or_blank_because_they_become_one_route_option():
    for solvent in gk.SOLVENTS:
        assert "," not in solvent and " " not in solvent


def test_fukui_force_constants_exclude_reading_the_checkpoint():
    assert set(gk.FUKUI_FORCE_CONSTANTS) < set(gk.FORCE_CONSTANTS)
    assert not any(gk.FORCE_CONSTANTS[label] in gk.READ_FORCE_CONSTANTS for label in gk.FUKUI_FORCE_CONSTANTS)


def test_flag_fields_name_route_options():
    assert gk.FREQ_FLAG_FIELDS == {"freq_raman": "Raman", "freq_vcd": "VCD", "freq_anharmonic": "Anharmonic"}
    assert set(gk.GUESS_FLAG_FIELDS.values()) == {"Mix", "Always", "Only", "Save"}
    assert set(gk.POP_FLAG_FIELDS.values()) == {"MK", "Hirshfeld", "NBO", "NPA", "SaveNTO"}


def test_every_panel_keyword_has_an_owner():
    for key in ("OPT", "FREQ", "IRC", "STABLE", "NMR", "TD", "TDA", "EMPIRICALDISPERSION", "IOP", "SCF",
                "NOSYMM", "POLAR", "GEOM", "GUESS", "POP", "SCRF"):
        assert key in gk.PANEL_OF_KEYWORD
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_gaussianKeywords.py -q`
Expected: collection error, `ModuleNotFoundError: No module named 'gaussianKeywords'`.

- [ ] **Step 3: Write `gaussianKeywords.py`**

```python
"""
Data tables of the Gaussian 16 input builder: what every control writes into the route.

Modelled on GaussView 6.1's Gaussian Calculation Setup (its help page gausscalcopts.htm lists each
control's keyword). Keys are the labels the window shows; values are Gaussian keywords or options
(None = write nothing). ``gaussianInput`` builds routes from these tables and reads routes back
through ``reverse``, so a label and its keyword are defined only here. No Toga.
"""

from typing import Dict, Mapping, Optional

from fukuiInput import DEFAULT_PRINT_LEVEL, PRINT_LEVELS

OTHER = "Other..."
NO_CHOICE = "(none)"          # the method list of a family that has no sub-choice (HF, MP2, ...)

# --- Job Type -----------------------------------------------------------------
JOB_ENERGY = "Energy"
JOB_OPT = "Optimization"
JOB_FREQ = "Frequency"
JOB_OPT_FREQ = "Opt+Freq"
JOB_IRC = "IRC"
JOB_SCAN = "Scan (relaxed)"
JOB_STABILITY = "Stability"
JOB_NMR = "NMR"
JOB_FUKUI = "Fukui functions (N, N+p, N-q)"
JOB_TYPES = (JOB_ENERGY, JOB_OPT, JOB_FREQ, JOB_OPT_FREQ, JOB_IRC, JOB_SCAN, JOB_STABILITY, JOB_NMR,
             JOB_FUKUI)
OPTIMIZING_JOBS = (JOB_OPT, JOB_OPT_FREQ, JOB_SCAN)
FREQUENCY_JOBS = (JOB_FREQ, JOB_OPT_FREQ)

OPT_MINIMUM = "Minimum"
OPT_TS = "TS (Berny)"
OPT_TARGETS: Dict[str, Optional[str]] = {OPT_MINIMUM: None, OPT_TS: "TS"}

FC_DEFAULT = "Default"
FC_CALCFC = "Calculate at first point (CalcFC)"
FORCE_CONSTANTS: Dict[str, Optional[str]] = {
    FC_DEFAULT: None,
    FC_CALCFC: "CalcFC",
    "Calculate at all points (CalcAll)": "CalcAll",
    "Read from checkpoint (ReadFC)": "ReadFC",
}
FUKUI_FORCE_CONSTANTS = (FC_DEFAULT, FC_CALCFC, "Calculate at all points (CalcAll)")
READ_FORCE_CONSTANTS = ("ReadFC", "RCFC")          # options that read the checkpoint

IRC_DIRECTIONS: Dict[str, Optional[str]] = {"Both": None, "Forward only": "Forward", "Reverse only": "Reverse"}
NMR_METHODS: Dict[str, Optional[str]] = {"GIAO": None, "CSGT": "CSGT", "IGAIM": "IGAIM", "All": "All"}
FREQ_FLAG_FIELDS = {"freq_raman": "Raman", "freq_vcd": "VCD", "freq_anharmonic": "Anharmonic"}

# --- Method -------------------------------------------------------------------
STATE_GROUND = "Ground state"
STATE_TD = "TD-SCF"
STATE_TDA = "TDA"
STATE_CIS = "CIS"
STATE_ZINDO = "ZINDO"
STATES = (STATE_GROUND, STATE_TD, STATE_TDA, STATE_CIS, STATE_ZINDO)
STATE_KEYWORDS = {STATE_TD: "TD", STATE_TDA: "TDA", STATE_CIS: "CIS", STATE_ZINDO: "ZINDO"}
MODEL_STATES = (STATE_CIS, STATE_ZINDO)            # written as the method of the model token: CIS(...)/basis
STATE_SPINS: Dict[str, Optional[str]] = {"Default": None, "Singlets": "Singlets", "Triplets": "Triplets",
                                         "50-50": "50-50"}
DEFAULT_N_STATES = 3                               # Gaussian's own default
COEFFICIENTS_IOP = "9/40=4"                        # print every CI coefficient (the NTO/Excitation tabs need it)

SPIN_PREFIXES: Dict[str, str] = {"Default": "", "Restricted (R)": "R", "Unrestricted (U)": "U",
                                 "Restricted-open (RO)": "RO"}

FAMILY_HF = "Hartree-Fock"
FAMILY_DFT = "DFT"
FAMILY_SEMI = "Semi-empirical"
FAMILY_MP2 = "MP2"
FAMILY_MP4 = "MP4(SDTQ)"
FAMILY_CCSD = "CCSD"
FAMILY_CCSDT = "CCSD(T)"
FAMILY_COMPOUND = "Compound"
FAMILIES = (FAMILY_HF, FAMILY_DFT, FAMILY_SEMI, FAMILY_MP2, FAMILY_MP4, FAMILY_CCSD, FAMILY_CCSDT,
            FAMILY_COMPOUND)
FAMILY_KEYWORDS = {FAMILY_HF: "HF", FAMILY_MP2: "MP2", FAMILY_MP4: "MP4(SDTQ)", FAMILY_CCSD: "CCSD",
                   FAMILY_CCSDT: "CCSD(T)"}       # families without a sub-choice
FAMILY_METHODS: Dict[str, Dict[str, Optional[str]]] = {
    FAMILY_DFT: {"B3LYP": "B3LYP", "CAM-B3LYP": "CAM-B3LYP", "wB97X-D": "wB97XD", "M06-2X": "M062X",
                 "M06": "M06", "PBE0 (PBE1PBE)": "PBE1PBE", "PBEPBE": "PBEPBE", "B3PW91": "B3PW91",
                 "BP86": "BP86", "TPSSTPSS": "TPSSTPSS", "APFD": "APFD", "HSEH1PBE": "HSEH1PBE",
                 "LC-wHPBE": "LC-wHPBE", OTHER: None},
    FAMILY_SEMI: {"PM7": "PM7", "PM6": "PM6", "AM1": "AM1", "PM3": "PM3", "PDDG": "PDDG"},
    FAMILY_COMPOUND: {"CBS-QB3": "CBS-QB3", "CBS-4M": "CBS-4M", "G4": "G4", "G4MP2": "G4MP2", "W1BD": "W1BD"},
}
NO_BASIS_FAMILIES = (FAMILY_SEMI, FAMILY_COMPOUND)
POST_SCF_FAMILIES = (FAMILY_MP2, FAMILY_MP4, FAMILY_CCSD, FAMILY_CCSDT)
DEFAULT_METHOD = {FAMILY_DFT: "B3LYP", FAMILY_SEMI: "PM7", FAMILY_COMPOUND: "CBS-QB3"}

POPLE_BASES = ("STO-3G", "3-21G", "6-31G", "6-311G")
EXTENSIBLE_POPLE = ("6-31G", "6-311G")             # the ones the diffuse/polarization lists apply to
OTHER_BASES = ("cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ",
               "Def2SVP", "Def2TZVP", "Def2TZVPP", "Def2QZVP", "LanL2DZ", "SDD", "DGDZVP", "Gen", "GenECP")
BASES = POPLE_BASES + OTHER_BASES + (OTHER,)
GENERAL_BASES = ("Gen", "GenECP")
DIFFUSE: Dict[str, str] = {"None": "", "+": "+", "++": "++"}
POLARIZATION: Dict[str, str] = {"None": "", "(d)": "(d)", "(d,p)": "(d,p)", "(2d,p)": "(2d,p)",
                                "(2df,2pd)": "(2df,2pd)", "(3df,3pd)": "(3df,3pd)"}
DISPERSION: Dict[str, Optional[str]] = {"None": None, "GD3": "GD3", "GD3BJ": "GD3BJ", "GD2": "GD2"}

# --- Molecule -----------------------------------------------------------------
GEOM_COORDINATES = "Coordinates (the step-1 geometry)"
GEOM_ALLCHECK = "Checkpoint: everything (Geom=AllCheck)"
GEOM_CHECK = "Checkpoint: new charge and multiplicity (Geom=Check)"
GEOMETRY_SOURCES: Dict[str, Optional[str]] = {GEOM_COORDINATES: None, GEOM_ALLCHECK: "AllCheck",
                                              GEOM_CHECK: "Check"}
MODRED_TYPES = {"Bond": ("B", 2), "Angle": ("A", 3), "Dihedral": ("D", 4)}
MODRED_ACTIONS = {"Scan": "S", "Freeze": "F"}

# --- Guess / Pop / Solvation --------------------------------------------------
GUESS_METHODS: Dict[str, Optional[str]] = {"Default": None, "Core Hamiltonian": "Core",
                                           "Extended Huckel": "Huckel", "Read checkpoint": "Read",
                                           "Read if possible (TCheck)": "TCheck"}
READ_GUESSES = ("Read",)                           # TCheck falls back to a new guess, so it needs no checkpoint
GUESS_FLAG_FIELDS = {"guess_mix": "Mix", "guess_always": "Always", "guess_only": "Only", "guess_save": "Save"}

POP_LEVELS: Dict[str, Optional[str]] = {"Default": None, "Full": "Full", "Regular": "Regular",
                                        "Minimal": "Minimal", "None": "None"}
POP_FLAG_FIELDS = {"pop_mk": "MK", "pop_hirshfeld": "Hirshfeld", "pop_nbo": "NBO", "pop_npa": "NPA",
                   "pop_savento": "SaveNTO"}

SOLVATION_NONE = "None (gas phase)"
SOLVATION_MODELS: Dict[str, Optional[str]] = {SOLVATION_NONE: None, "IEFPCM": "IEFPCM", "CPCM": "CPCM",
                                              "SMD": "SMD"}
SOLVATION_ALIASES = {"PCM": "IEFPCM"}              # read back: SCRF=PCM is IEFPCM in Gaussian 16
SOLVENTS = ("Water", "Acetonitrile", "Methanol", "Ethanol", "DiMethylSulfoxide", "Acetone", "Chloroform",
            "DiChloroMethane", "TetraHydroFuran", "Toluene", "Benzene", "n-Hexane", "CarbonTetraChloride",
            "DiethylEther", OTHER)

# Route keyword (upper case) -> the panel that writes it, for messages.
PANEL_OF_KEYWORD = {
    "OPT": "Job Type", "FREQ": "Job Type", "IRC": "Job Type", "STABLE": "Job Type", "NMR": "Job Type",
    "TD": "Method", "TDA": "Method", "EMPIRICALDISPERSION": "Method", "IOP": "Method",
    "SCF": "General", "NOSYMM": "General", "POLAR": "General", "GEOM": "Molecule/General",
    "GUESS": "Guess", "POP": "Pop", "SCRF": "Solvation",
}

def reverse(table: Mapping[str, Optional[str]]) -> Dict[str, str]:
    """Route value (upper case) -> label; entries that write nothing are skipped."""
    return {value.upper(): label for label, value in table.items() if value}
```

`PRINT_LEVELS` and `DEFAULT_PRINT_LEVEL` are imported from `fukuiInput` so both builders share one list; the window reads them as `gk.PRINT_LEVELS`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_gaussianKeywords.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add gaussianKeywords.py tests/test_gaussianKeywords.py
git commit -m "feat(gaussian): keyword tables of the Gaussian 16 input builder"
```

---

### Task 2: Route keyword algebra (`gaussianInput.py`, part 1)

**Files:**
- Create: `gaussianInput.py`
- Test: `tests/test_gaussianInput.py`

**Interfaces:**
- Consumes: `fukuiInput.route_tokens`, `gaussianKeywords.PANEL_OF_KEYWORD`.
- Produces: `RouteKeyword(name: str, options: List[str])` with property `key -> str` (upper-case name); `option_key(option: str) -> str`; `parse_keyword(token: str) -> RouteKeyword`; `format_keyword(keyword: RouteKeyword) -> str`; `merge_keywords(panel: Sequence[RouteKeyword], extra_text: str) -> Tuple[List[RouteKeyword], List[str]]` (merged keywords, error messages).

- [ ] **Step 1: Write the failing tests**

`tests/test_gaussianInput.py` (this file grows in Tasks 3–6):

```python
"""Tests for the Gaussian 16 input builder core (gaussianInput, no Toga)."""
import subprocess
import sys
from pathlib import Path

import pytest

import gaussianInput as gi
import gaussianKeywords as gk
from gaussianInput import RouteKeyword as K

SRC = Path(__file__).resolve().parents[1]


def test_numerical_modules_do_not_import_toga():
    code = "import sys, gaussianInput, gaussianKeywords; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=str(SRC))
    assert out.stdout.strip() == "False", out.stderr


class TestKeywordAlgebra:
    @pytest.mark.parametrize("token,name,options", [
        ("NoSymm", "NoSymm", []),
        ("Opt=Tight", "Opt", ["Tight"]),
        ("Opt=(TS,CalcFC)", "Opt", ["TS", "CalcFC"]),
        ("SCF=(Tight, XQC)", "SCF", ["Tight", "XQC"]),
        ("TD(NStates=4)", "TD", ["NStates=4"]),
        ("TD(50-50,NStates=4)", "TD", ["50-50", "NStates=4"]),
        ("IOp(9/40=4)", "IOp", ["9/40=4"]),
        ("SCRF=(SMD,Solvent=Water)", "SCRF", ["SMD", "Solvent=Water"]),
        ("Integral(Grid=UltraFine)", "Integral", ["Grid=UltraFine"]),
    ])
    def test_parse(self, token, name, options):
        assert gi.parse_keyword(token) == K(name, options)

    @pytest.mark.parametrize("keyword,text", [
        (K("NoSymm"), "NoSymm"),
        (K("Opt", ["Tight"]), "Opt=Tight"),
        (K("Opt", ["TS", "CalcFC"]), "Opt=(TS,CalcFC)"),
        (K("SCRF", ["Solvent=Water"]), "SCRF=(Solvent=Water)"),
        (K("IOp", ["9/40=4"]), "IOp(9/40=4)"),
        (K("TD", ["NStates=10"]), "TD=(NStates=10)"),
    ])
    def test_format(self, keyword, text):
        assert gi.format_keyword(keyword) == text

    def test_key_and_option_key_ignore_case_and_value(self):
        assert K("opt").key == "OPT"
        assert gi.option_key("MaxCycles=50") == "MAXCYCLES" and gi.option_key("tight") == "TIGHT"

    def test_new_keywords_are_appended(self):
        merged, errors = gi.merge_keywords([K("Opt")], "Int=UltraFine 6D")
        assert [gi.format_keyword(k) for k in merged] == ["Opt", "Int=UltraFine", "6D"] and errors == []

    def test_options_of_a_panel_keyword_are_merged(self):
        merged, errors = gi.merge_keywords([K("Opt", ["TS", "CalcFC"])], "Opt=NoEigenTest")
        assert gi.format_keyword(merged[0]) == "Opt=(TS,CalcFC,NoEigenTest)" and errors == []

    def test_a_repeated_option_is_kept_once(self):
        merged, errors = gi.merge_keywords([K("Opt", ["Tight"])], "opt=tight")
        assert merged == [K("Opt", ["Tight"])] and errors == []

    def test_a_valued_option_with_another_value_is_an_error_naming_both_places(self):
        _, errors = gi.merge_keywords([K("Opt", ["MaxCycles=100"])], "Opt=MaxCycles=50")
        assert errors == ["Additional Keywords: Opt=MaxCycles=50 sets MaxCycles=50, but the Job Type panel "
                          "sets MaxCycles=100."]

    def test_unvalued_options_are_not_checked_against_each_other(self):
        merged, errors = gi.merge_keywords([], "Int=UltraFine Int=Fine")
        assert merged == [K("Int", ["UltraFine", "Fine"])] and errors == []

    def test_a_valued_conflict_between_two_additional_keywords(self):
        _, errors = gi.merge_keywords([], "SCF=MaxCycle=64 SCF=MaxCycle=128")
        assert errors == ["Additional Keywords: SCF=MaxCycle=128 sets MaxCycle=128, but an earlier "
                          "Additional Keyword sets MaxCycle=64."]

    def test_the_model_chemistry_is_refused_in_additional_keywords(self):
        _, errors = gi.merge_keywords([], "B3LYP/6-31G")
        assert errors and "Method panel" in errors[0]

    def test_iop_with_a_slash_is_not_a_model_chemistry(self):
        merged, errors = gi.merge_keywords([], "IOp(9/40=4)")
        assert merged == [K("IOp", ["9/40=4"])] and errors == []
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_gaussianInput.py -q`
Expected: `ModuleNotFoundError: No module named 'gaussianInput'`.

- [ ] **Step 3: Write the first part of `gaussianInput.py`**

```python
"""
Gaussian 16 input builder: the job model, validation, route, input text and reopening (no Toga).

``GaussianJobSpec`` (job-wide fields + one ``GaussianStepSpec`` per linked step) ->
``GaussianInputBuilder.prepare`` (every problem in one ValueError, one '- Step n: ...' line each) ->
``GaussianJob`` -> ``input_text`` / ``script_text`` / ``write_files``. ``read_input_text`` is the
inverse used by Open input. The Fukui Job Type is delegated to ``fukuiInput`` unchanged. Files are
ASCII with LF endings and relative checkpoint names: the job runs under Linux from the input's folder.
"""

import re
from dataclasses import dataclass, field
from typing import List, Sequence, Tuple

import fukuiInput as fi
import gaussianKeywords as gk

PARENTHESIS_KEYWORDS = ("IOP",)        # written as IOp(9/40=4), the form every Gaussian manual uses


# ----------------------------------------------------------------------
# Route keywords
# ----------------------------------------------------------------------
@dataclass
class RouteKeyword:
    """One route keyword and its options, e.g. Opt=(TS,CalcFC) -> ('Opt', ['TS', 'CalcFC'])."""
    name: str
    options: List[str] = field(default_factory=list)

    @property
    def key(self) -> str:
        return self.name.upper()


def option_key(option: str) -> str:
    """'MaxCycles=50' -> 'MAXCYCLES': two options with the same key set the same thing."""
    return option.split("=", 1)[0].strip().upper()


def _outer_parentheses(text: str) -> str:
    """'(a,b)' -> 'a,b' when the first '(' closes at the very end; anything else unchanged."""
    if not (text.startswith("(") and text.endswith(")")):
        return text
    depth = 0
    for position, ch in enumerate(text):
        depth += (ch == "(") - (ch == ")")
        if depth == 0 and position < len(text) - 1:
            return text
    return text[1:-1]


def _split_options(text: str) -> List[str]:
    """Split on commas and blanks outside parentheses."""
    options, current, depth = [], "", 0
    for ch in text:
        if depth == 0 and (ch == "," or ch.isspace()):
            if current:
                options.append(current)
            current = ""
            continue
        depth += (ch == "(") - (ch == ")")
        current += ch
    if current:
        options.append(current)
    return options


def keyword_text_name(token: str) -> str:
    """The keyword part of a route token as written: 'Opt=(TS)' -> 'Opt', 'IOp(9/40=4)' -> 'IOp'."""
    return re.split(r"[=(]", token, maxsplit=1)[0]


def parse_keyword(token: str) -> RouteKeyword:
    name = keyword_text_name(token)
    rest = token[len(name):]
    if rest.startswith("="):
        rest = rest[1:]
    return RouteKeyword(name, _split_options(_outer_parentheses(rest)))


def format_keyword(keyword: RouteKeyword) -> str:
    options = keyword.options
    if not options:
        return keyword.name
    if keyword.key in PARENTHESIS_KEYWORDS:
        return f"{keyword.name}({','.join(options)})"
    if len(options) == 1 and "=" not in options[0]:
        return f"{keyword.name}={options[0]}"
    return f"{keyword.name}=({','.join(options)})"


def merge_keywords(panel: Sequence[RouteKeyword], extra_text: str) -> Tuple[List[RouteKeyword], List[str]]:
    """The panels' keywords plus Additional Keywords (spec 5.1).

    A keyword already written by a panel takes the extra options; an option repeated with the same
    key is kept once, and a valued option (X=...) given another value is an error naming both places.
    """
    merged = [RouteKeyword(k.name, list(k.options)) for k in panel]
    panel_keys = {k.key for k in panel}
    errors: List[str] = []
    for token in fi.route_tokens(extra_text):
        if "/" in keyword_text_name(token):
            errors.append(f"Additional Keywords: '{token}' is a method/basis set; choose them in the "
                          "Method panel.")
            continue
        extra = parse_keyword(token)
        target = next((k for k in merged if k.key == extra.key), None)
        if target is None:
            merged.append(extra)
            continue
        for option in extra.options:
            same = [o for o in target.options if option_key(o) == option_key(option)]
            if not same:
                target.options.append(option)
            elif "=" in option and same[0].upper() != option.upper():
                owner = (f"the {gk.PANEL_OF_KEYWORD.get(extra.key, 'Job Type')} panel"
                         if extra.key in panel_keys else "an earlier Additional Keyword")
                errors.append(f"Additional Keywords: {token} sets {option}, but {owner} sets {same[0]}.")
    return merged, errors
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_gaussianInput.py tests/test_gaussianKeywords.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add gaussianInput.py tests/test_gaussianInput.py
git commit -m "feat(gaussian): route keyword parsing, formatting and merging"
```

---

### Task 3: Step spec and route building (`gaussianInput.py`, part 2)

**Files:**
- Modify: `gaussianInput.py` (append)
- Test: `tests/test_gaussianInput.py` (append)

**Interfaces:**
- Consumes: Task 1 tables; Task 2 `RouteKeyword`, `format_keyword`, `merge_keywords`; `fukuiInput.build_route`, `DEFAULT_MEMORY`, `DEFAULT_NPROC`, `DEFAULT_PRINT_LEVEL`.
- Produces: `GaussianStepSpec` (dataclass, every field defaulted, listed below); `FIELD_LABELS: Dict[str, str]` (one per `GaussianStepSpec` field); `CHOICES: Dict[str, Sequence[str]]`; `method_choices(family) -> Tuple[str, ...]`; `choices(field_name, step=None) -> Tuple[str, ...]`; `no_basis(step) -> bool`; `state_options(step) -> List[str]`; `method_keyword(step, with_prefix=True) -> str`; `basis_keyword(step) -> str`; `model_token(step) -> str`; `opt_options(step) -> List[str]`; `job_keywords(step) -> List[RouteKeyword]`; `method_keywords(step) -> List[RouteKeyword]`; `general_keywords(step) -> List[RouteKeyword]`; `panel_keywords(step) -> List[RouteKeyword]`; `route_keywords(step) -> Tuple[List[RouteKeyword], List[str]]`; `route_lines(step) -> List[str]`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_gaussianInput.py`)

```python
from dataclasses import fields, replace


def step(**changes) -> gi.GaussianStepSpec:
    return replace(gi.GaussianStepSpec(), **changes)


def route(**changes) -> str:
    return " ".join(gi.route_lines(step(**changes)))


class TestStepSpec:
    def test_every_field_has_a_label(self):
        assert set(gi.FIELD_LABELS) == {f.name for f in fields(gi.GaussianStepSpec)}

    def test_every_choice_default_is_one_of_its_choices(self):
        default = gi.GaussianStepSpec()
        for name in gi.CHOICES:
            assert getattr(default, name) in gi.choices(name, default), name
        assert default.method_name in gi.choices("method_name", default)

    def test_method_choices_follow_the_family(self):
        assert gi.method_choices(gk.FAMILY_SEMI)[0] == "PM7"
        assert gi.method_choices(gk.FAMILY_HF) == (gk.NO_CHOICE,)


class TestModel:
    def test_default_route(self):
        assert route() == "#P B3LYP/6-31G(d)"

    @pytest.mark.parametrize("changes,model", [
        (dict(basis_diffuse="+", basis_polarization="(d,p)"), "B3LYP/6-31+G(d,p)"),
        (dict(basis="6-311G", basis_diffuse="++", basis_polarization="(2df,2pd)"), "B3LYP/6-311++G(2df,2pd)"),
        (dict(basis="6-31G", basis_polarization="None"), "B3LYP/6-31G"),
        (dict(basis="STO-3G", basis_diffuse="+"), "B3LYP/STO-3G"),          # lists apply to 6-31G/6-311G only
        (dict(basis="aug-cc-pVTZ"), "B3LYP/aug-cc-pVTZ"),
        (dict(basis=gk.OTHER, basis_other="def2-TZVPD"), "B3LYP/def2-TZVPD"),
        (dict(spin="Unrestricted (U)"), "UB3LYP/6-31G(d)"),
        (dict(spin="Restricted-open (RO)", family=gk.FAMILY_HF, method_name=gk.NO_CHOICE), "ROHF/6-31G(d)"),
        (dict(method_name="M06-2X"), "M062X/6-31G(d)"),
        (dict(method_name=gk.OTHER, method_other="revPBE"), "revPBE/6-31G(d)"),
        (dict(family=gk.FAMILY_SEMI, method_name="PM6"), "PM6"),
        (dict(family=gk.FAMILY_COMPOUND, method_name="CBS-QB3"), "CBS-QB3"),
        (dict(family=gk.FAMILY_CCSDT, method_name=gk.NO_CHOICE, basis="cc-pVDZ"), "CCSD(T)/cc-pVDZ"),
        (dict(state=gk.STATE_CIS, n_states=3, basis="cc-pVDZ"), "CIS(NStates=3)/cc-pVDZ"),
        (dict(state=gk.STATE_ZINDO, n_states=10, state_spin="Singlets"), "ZINDO(NStates=10,Singlets)"),
    ])
    def test_model_token(self, changes, model):
        assert gi.model_token(step(**changes)) == model


class TestJobKeywords:
    @pytest.mark.parametrize("changes,expected", [
        (dict(job_type=gk.JOB_OPT), "#P B3LYP/6-31G(d) Opt"),
        (dict(job_type=gk.JOB_OPT, opt_target=gk.OPT_TS, force_constants=gk.FC_CALCFC, opt_tight=True,
              opt_max_cycles=200), "#P B3LYP/6-31G(d) Opt=(TS,CalcFC,Tight,MaxCycles=200)"),
        (dict(job_type=gk.JOB_FREQ, freq_raman=True), "#P B3LYP/6-31G(d) Freq=Raman"),
        (dict(job_type=gk.JOB_OPT_FREQ, freq_vcd=True, freq_anharmonic=True),
         "#P B3LYP/6-31G(d) Opt Freq=(VCD,Anharmonic)"),
        (dict(job_type=gk.JOB_IRC, irc_direction="Forward only", irc_max_points=20, irc_step_size=5),
         "#P B3LYP/6-31G(d) IRC=(Forward,CalcFC,MaxPoints=20,StepSize=5)"),
        (dict(job_type=gk.JOB_SCAN), "#P B3LYP/6-31G(d) Opt=ModRedundant"),
        (dict(job_type=gk.JOB_OPT, modredundant="B 1 2 F"), "#P B3LYP/6-31G(d) Opt=ModRedundant"),
        (dict(job_type=gk.JOB_STABILITY, stable_opt=True), "#P B3LYP/6-31G(d) Stable=Opt"),
        (dict(job_type=gk.JOB_NMR, nmr_method="CSGT", nmr_spinspin=True), "#P B3LYP/6-31G(d) NMR=(CSGT,SpinSpin)"),
        (dict(job_type=gk.JOB_NMR), "#P B3LYP/6-31G(d) NMR"),
    ])
    def test_job_type_route(self, changes, expected):
        assert route(**changes) == expected

    def test_energy_writes_no_job_keyword(self):
        assert gi.job_keywords(step()) == []


class TestOtherPanels:
    @pytest.mark.parametrize("changes,expected", [
        (dict(state=gk.STATE_TD, n_states=10), "#P B3LYP/6-31G(d) TD=(NStates=10)"),
        (dict(state=gk.STATE_TDA, n_states=6, root=2, state_spin="Triplets"),
         "#P B3LYP/6-31G(d) TDA=(NStates=6,Root=2,Triplets)"),
        (dict(state=gk.STATE_TD, print_all_coefficients=True), "#P B3LYP/6-31G(d) TD=(NStates=3) IOp(9/40=4)"),
        (dict(dispersion="GD3BJ"), "#P B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ"),
        (dict(scf_tight=True, scf_qc=True, nosymm=True, polar=True), "#P B3LYP/6-31G(d) SCF=(Tight,QC) NoSymm Polar"),
        (dict(geometry_source=gk.GEOM_ALLCHECK), "#P B3LYP/6-31G(d) Geom=AllCheck"),
        (dict(connectivity=True), "#P B3LYP/6-31G(d) Geom=Connectivity"),
        (dict(guess="Read checkpoint", guess_mix=True), "#P B3LYP/6-31G(d) Guess=(Read,Mix)"),
        (dict(pop_level="Full", pop_mk=True, pop_nbo=True), "#P B3LYP/6-31G(d) Pop=(Full,MK,NBO)"),
        (dict(solvation="SMD", solvent="Acetonitrile"), "#P B3LYP/6-31G(d) SCRF=(SMD,Solvent=Acetonitrile)"),
        (dict(solvation="CPCM", solvent=gk.OTHER, solvent_other="Heptane"),
         "#P B3LYP/6-31G(d) SCRF=(CPCM,Solvent=Heptane)"),
        (dict(print_level="#N"), "#N B3LYP/6-31G(d)"),
    ])
    def test_route(self, changes, expected):
        assert route(**changes) == expected

    def test_additional_keywords_are_merged(self):
        assert route(job_type=gk.JOB_OPT, opt_target=gk.OPT_TS, additional_keywords="Opt=NoEigenTest 6D") == \
            "#P B3LYP/6-31G(d) Opt=(TS,NoEigenTest) 6D"

    def test_route_keywords_report_merge_errors(self):
        _, errors = gi.route_keywords(step(job_type=gk.JOB_OPT, opt_max_cycles=100,
                                           additional_keywords="Opt=MaxCycles=50"))
        assert errors and "Job Type panel" in errors[0]

    def test_a_long_route_is_wrapped_below_80_columns(self):
        lines = gi.route_lines(step(job_type=gk.JOB_OPT_FREQ, opt_tight=True, freq_raman=True,
                                    state=gk.STATE_TD, n_states=10, print_all_coefficients=True,
                                    dispersion="GD3BJ", scf_tight=True, nosymm=True, pop_level="Full",
                                    solvation="SMD", solvent="DiMethylSulfoxide"))
        assert len(lines) > 1 and all(len(line) <= fi.ROUTE_WIDTH for line in lines)
```

Add `import fukuiInput as fi` to the imports at the top of the test file.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_gaussianInput.py -q`
Expected: the new tests fail with `AttributeError: module 'gaussianInput' has no attribute 'GaussianStepSpec'`.

- [ ] **Step 3: Implement** (append to `gaussianInput.py`; add `Dict, Optional, Tuple` to the `typing` import)

```python
# ----------------------------------------------------------------------
# One step's panels
# ----------------------------------------------------------------------
@dataclass
class GaussianStepSpec:
    """Every panel of one linked step, as the window shows it (labels, not keywords).

    A multiplicity of None means "not given": ``prepare`` asks for it rather than assuming one.
    """
    # Job Type
    job_type: str = gk.JOB_ENERGY
    opt_target: str = gk.OPT_MINIMUM
    force_constants: str = gk.FC_DEFAULT
    opt_tight: bool = False
    opt_max_cycles: Optional[int] = None
    freq_raman: bool = False
    freq_vcd: bool = False
    freq_anharmonic: bool = False
    irc_direction: str = "Both"
    irc_force_constants: str = gk.FC_CALCFC
    irc_max_points: Optional[int] = None
    irc_step_size: Optional[int] = None
    stable_opt: bool = False
    nmr_method: str = "GIAO"
    nmr_spinspin: bool = False
    fukui_p: int = 1
    fukui_q: int = 1
    fukui_plus_multiplicity: Optional[int] = None
    fukui_minus_multiplicity: Optional[int] = None
    fukui_guess_read: bool = True
    # Method
    state: str = gk.STATE_GROUND
    n_states: int = gk.DEFAULT_N_STATES
    root: int = 1
    state_spin: str = "Default"
    print_all_coefficients: bool = False
    spin: str = "Default"
    family: str = gk.FAMILY_DFT
    method_name: str = "B3LYP"
    method_other: str = ""
    basis: str = "6-31G"
    basis_diffuse: str = "None"
    basis_polarization: str = "(d)"
    basis_other: str = ""
    dispersion: str = "None"
    # Molecule
    geometry_source: str = gk.GEOM_COORDINATES
    charge: int = 0
    multiplicity: Optional[int] = None
    modredundant: str = ""
    # Title
    title: str = ""
    # Link 0
    memory: str = fi.DEFAULT_MEMORY
    nproc: int = fi.DEFAULT_NPROC
    write_chk: bool = True
    chk: str = ""
    old_chk: str = ""
    extra_link0: str = ""
    # General
    print_level: str = fi.DEFAULT_PRINT_LEVEL
    scf_tight: bool = False
    scf_qc: bool = False
    nosymm: bool = False
    polar: bool = False
    connectivity: bool = False
    # Guess
    guess: str = "Default"
    guess_mix: bool = False
    guess_always: bool = False
    guess_only: bool = False
    guess_save: bool = False
    # Pop
    pop_level: str = "Default"
    pop_mk: bool = False
    pop_hirshfeld: bool = False
    pop_nbo: bool = False
    pop_npa: bool = False
    pop_savento: bool = False
    # Solvation
    solvation: str = gk.SOLVATION_NONE
    solvent: str = "Water"
    solvent_other: str = ""
    # Add. Inp. and the Additional Keywords field
    additional_input: str = ""
    additional_keywords: str = ""


FIELD_LABELS: Dict[str, str] = {
    "job_type": "Job type", "opt_target": "Optimize to a", "force_constants": "Force constants",
    "opt_tight": "Tight convergence", "opt_max_cycles": "Maximum optimization cycles",
    "freq_raman": "Raman", "freq_vcd": "VCD", "freq_anharmonic": "Anharmonic",
    "irc_direction": "IRC direction", "irc_force_constants": "IRC force constants",
    "irc_max_points": "IRC MaxPoints", "irc_step_size": "IRC StepSize", "stable_opt": "Reoptimize",
    "nmr_method": "NMR method", "nmr_spinspin": "Spin-spin couplings",
    "fukui_p": "p (electrons added)", "fukui_q": "q (electrons removed)",
    "fukui_plus_multiplicity": "Multiplicity of N+p", "fukui_minus_multiplicity": "Multiplicity of N-q",
    "fukui_guess_read": "Ions start from the N-state orbitals",
    "state": "Electronic state", "n_states": "Number of states", "root": "State of interest (Root)",
    "state_spin": "States", "print_all_coefficients": "Print all excitation coefficients",
    "spin": "Spin", "family": "Method family", "method_name": "Method", "method_other": "Other method",
    "basis": "Basis set", "basis_diffuse": "Diffuse functions", "basis_polarization": "Polarization",
    "basis_other": "Other basis set", "dispersion": "Empirical dispersion",
    "geometry_source": "Geometry source", "charge": "Charge", "multiplicity": "Multiplicity",
    "modredundant": "Scan/freeze lines", "title": "Title", "memory": "%Mem", "nproc": "%NProcShared",
    "write_chk": "Write a checkpoint (%Chk)", "chk": "%Chk", "old_chk": "%OldChk",
    "extra_link0": "Extra Link 0 lines", "print_level": "Print level", "scf_tight": "SCF=Tight",
    "scf_qc": "SCF=QC", "nosymm": "NoSymm", "polar": "Polar", "connectivity": "Write connectivity",
    "guess": "Initial guess", "guess_mix": "Mix HOMO and LUMO", "guess_always": "Always generate",
    "guess_only": "Only do guess", "guess_save": "Save orbitals", "pop_level": "Population output",
    "pop_mk": "MK charges", "pop_hirshfeld": "Hirshfeld (CM5) charges", "pop_nbo": "NBO", "pop_npa": "NPA",
    "pop_savento": "Save NTOs", "solvation": "Solvation model", "solvent": "Solvent",
    "solvent_other": "Other solvent", "additional_input": "Additional input",
    "additional_keywords": "Additional Keywords",
}

# Fields chosen from a list; method_name depends on the family (method_choices).
CHOICES: Dict[str, Sequence[str]] = {
    "job_type": gk.JOB_TYPES, "opt_target": tuple(gk.OPT_TARGETS), "force_constants": tuple(gk.FORCE_CONSTANTS),
    "irc_direction": tuple(gk.IRC_DIRECTIONS), "irc_force_constants": tuple(gk.FORCE_CONSTANTS),
    "nmr_method": tuple(gk.NMR_METHODS), "state": gk.STATES, "state_spin": tuple(gk.STATE_SPINS),
    "spin": tuple(gk.SPIN_PREFIXES), "family": gk.FAMILIES, "basis": gk.BASES,
    "basis_diffuse": tuple(gk.DIFFUSE), "basis_polarization": tuple(gk.POLARIZATION),
    "dispersion": tuple(gk.DISPERSION), "geometry_source": tuple(gk.GEOMETRY_SOURCES),
    "print_level": gk.PRINT_LEVELS, "guess": tuple(gk.GUESS_METHODS), "pop_level": tuple(gk.POP_LEVELS),
    "solvation": tuple(gk.SOLVATION_MODELS), "solvent": gk.SOLVENTS,
}


def method_choices(family: str) -> Tuple[str, ...]:
    return tuple(gk.FAMILY_METHODS[family]) if family in gk.FAMILY_METHODS else (gk.NO_CHOICE,)


def choices(field_name: str, step: Optional[GaussianStepSpec] = None) -> Tuple[str, ...]:
    if field_name == "method_name":
        return method_choices((step or GaussianStepSpec()).family)
    return tuple(CHOICES[field_name])


# ----------------------------------------------------------------------
# Route of one step
# ----------------------------------------------------------------------
def no_basis(step: GaussianStepSpec) -> bool:
    """ZINDO, semi-empirical and compound methods take no basis set."""
    if step.state == gk.STATE_ZINDO:
        return True
    return step.state != gk.STATE_CIS and step.family in gk.NO_BASIS_FAMILIES


def state_options(step: GaussianStepSpec) -> List[str]:
    if step.state == gk.STATE_GROUND:
        return []
    options = [f"NStates={step.n_states}"]
    if step.root != 1:
        options.append(f"Root={step.root}")
    spin = gk.STATE_SPINS.get(step.state_spin)
    if spin:
        options.append(spin)
    return options


def method_keyword(step: GaussianStepSpec, with_prefix: bool = True) -> str:
    """'UB3LYP', 'CIS(NStates=3)', 'PM6', ... (the part before the '/')."""
    prefix = gk.SPIN_PREFIXES.get(step.spin, "") if with_prefix else ""
    if step.state in gk.MODEL_STATES:
        name, options = gk.STATE_KEYWORDS[step.state], state_options(step)
        return prefix + (f"{name}({','.join(options)})" if options else name)
    if step.family in gk.FAMILY_KEYWORDS:
        return prefix + gk.FAMILY_KEYWORDS[step.family]
    keyword = gk.FAMILY_METHODS.get(step.family, {}).get(step.method_name)
    return prefix + (keyword if keyword else step.method_other.strip())


def basis_keyword(step: GaussianStepSpec) -> str:
    if no_basis(step):
        return ""
    if step.basis == gk.OTHER:
        return step.basis_other.strip()
    if step.basis in gk.EXTENSIBLE_POPLE:
        diffuse = gk.DIFFUSE.get(step.basis_diffuse, "")
        polarization = gk.POLARIZATION.get(step.basis_polarization, "")
        return f"{step.basis[:-1]}{diffuse}G{polarization}"
    return step.basis


def model_token(step: GaussianStepSpec) -> str:
    method, basis = method_keyword(step), basis_keyword(step)
    return f"{method}/{basis}" if basis else method


def opt_options(step: GaussianStepSpec) -> List[str]:
    options = [o for o in (gk.OPT_TARGETS.get(step.opt_target), gk.FORCE_CONSTANTS.get(step.force_constants)) if o]
    if step.opt_tight:
        options.append("Tight")
    if step.opt_max_cycles:
        options.append(f"MaxCycles={step.opt_max_cycles}")
    if step.job_type == gk.JOB_SCAN or step.modredundant.strip():
        options.append("ModRedundant")
    return options


def _flags(step: GaussianStepSpec, table: Dict[str, str]) -> List[str]:
    return [option for name, option in table.items() if getattr(step, name)]


def job_keywords(step: GaussianStepSpec) -> List[RouteKeyword]:
    job, out = step.job_type, []
    if job in gk.OPTIMIZING_JOBS:
        out.append(RouteKeyword("Opt", opt_options(step)))
    if job in gk.FREQUENCY_JOBS:
        out.append(RouteKeyword("Freq", _flags(step, gk.FREQ_FLAG_FIELDS)))
    if job == gk.JOB_IRC:
        options = [o for o in (gk.IRC_DIRECTIONS.get(step.irc_direction),
                               gk.FORCE_CONSTANTS.get(step.irc_force_constants)) if o]
        if step.irc_max_points:
            options.append(f"MaxPoints={step.irc_max_points}")
        if step.irc_step_size:
            options.append(f"StepSize={step.irc_step_size}")
        out.append(RouteKeyword("IRC", options))
    if job == gk.JOB_STABILITY:
        out.append(RouteKeyword("Stable", ["Opt"] if step.stable_opt else []))
    if job == gk.JOB_NMR:
        method = gk.NMR_METHODS.get(step.nmr_method)
        out.append(RouteKeyword("NMR", ([method] if method else []) + (["SpinSpin"] if step.nmr_spinspin else [])))
    return out


def method_keywords(step: GaussianStepSpec) -> List[RouteKeyword]:
    """TD/TDA, dispersion and the coefficient IOp (CIS and ZINDO live in the model token)."""
    out = []
    if step.state in (gk.STATE_TD, gk.STATE_TDA):
        out.append(RouteKeyword(gk.STATE_KEYWORDS[step.state], state_options(step)))
    dispersion = gk.DISPERSION.get(step.dispersion)
    if dispersion:
        out.append(RouteKeyword("EmpiricalDispersion", [dispersion]))
    if step.print_all_coefficients:
        out.append(RouteKeyword("IOp", [gk.COEFFICIENTS_IOP]))
    return out


def general_keywords(step: GaussianStepSpec) -> List[RouteKeyword]:
    """General, Molecule (Geom), Guess, Pop and Solvation panels, in that order."""
    out = []
    scf = (["Tight"] if step.scf_tight else []) + (["QC"] if step.scf_qc else [])
    if scf:
        out.append(RouteKeyword("SCF", scf))
    if step.nosymm:
        out.append(RouteKeyword("NoSymm"))
    if step.polar:
        out.append(RouteKeyword("Polar"))
    source = gk.GEOMETRY_SOURCES.get(step.geometry_source)
    geom = ([source] if source else []) + (["Connectivity"] if step.connectivity else [])
    if geom:
        out.append(RouteKeyword("Geom", geom))
    guess = gk.GUESS_METHODS.get(step.guess)
    guess_options = ([guess] if guess else []) + _flags(step, gk.GUESS_FLAG_FIELDS)
    if guess_options:
        out.append(RouteKeyword("Guess", guess_options))
    level = gk.POP_LEVELS.get(step.pop_level)
    pop = ([level] if level else []) + _flags(step, gk.POP_FLAG_FIELDS)
    if pop:
        out.append(RouteKeyword("Pop", pop))
    model = gk.SOLVATION_MODELS.get(step.solvation)
    if model:
        solvent = step.solvent_other.strip() if step.solvent == gk.OTHER else step.solvent
        out.append(RouteKeyword("SCRF", [model, f"Solvent={solvent}"]))
    return out


def panel_keywords(step: GaussianStepSpec) -> List[RouteKeyword]:
    return job_keywords(step) + method_keywords(step) + general_keywords(step)


def route_keywords(step: GaussianStepSpec) -> Tuple[List[RouteKeyword], List[str]]:
    """The panels' keywords merged with Additional Keywords; second item = merge errors."""
    return merge_keywords(panel_keywords(step), step.additional_keywords)


def route_lines(step: GaussianStepSpec) -> List[str]:
    """The route section, wrapped below 80 columns (merge errors are ignored here; prepare reports them)."""
    keywords, _ = route_keywords(step)
    return fi.build_route(step.print_level, method_keyword(step), basis_keyword(step),
                          " ".join(format_keyword(k) for k in keywords))
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_gaussianInput.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add gaussianInput.py tests/test_gaussianInput.py
git commit -m "feat(gaussian): step spec and route building for every panel"
```

---

### Task 4: Scan / freeze (ModRedundant) lines

**Files:**
- Modify: `gaussianInput.py` (append)
- Test: `tests/test_gaussianInput.py` (append)

**Interfaces:**
- Consumes: `gk.MODRED_TYPES`, `gk.MODRED_ACTIONS`.
- Produces: `parse_modredundant(text: str, n_atoms: Optional[int] = None) -> List[Tuple[str, Tuple[int, ...], str]]` — `(letter, atoms, action)` per non-blank line, action `"S"` or `"F"`, `ValueError` naming the line; `modredundant_line(kind: str, atoms_text: str, action: str, steps_text: str, size_text: str) -> str` — the helper row's line, validated by the same parser.

- [ ] **Step 1: Write the failing tests** (append)

```python
class TestModRedundant:
    def test_scan_and_freeze_lines(self):
        assert gi.parse_modredundant("B 1 2 S 10 0.05\n\nd 1 2 3 4 f\nA 1 2 3 F", n_atoms=4) == [
            ("B", (1, 2), "S"), ("D", (1, 2, 3, 4), "F"), ("A", (1, 2, 3), "F")]

    @pytest.mark.parametrize("text,fragment", [
        ("X 1 2 F", "start with B, A or D"),
        ("B 1 F", "needs 2 atom numbers"),
        ("A 1 1 2 F", "different atoms"),
        ("B 0 2 F", "different atoms"),
        ("B 1 9 F", "atom 9 does not exist (the molecule has 3 atoms)"),
        ("B 1 2", "end with F (freeze) or S <steps> <step size> (scan)"),
        ("B 1 2 S 0 0.1", "steps must be a positive integer"),
        ("B 1 2 S 5 0", "step size must be a non-zero number"),
        ("B 1 2 S 5 x", "step size must be a non-zero number"),
    ])
    def test_problems_name_the_line(self, text, fragment):
        with pytest.raises(ValueError, match=re.escape(fragment)):
            gi.parse_modredundant("B 1 2 F\n" + text, n_atoms=3)
        with pytest.raises(ValueError, match="line 2"):
            gi.parse_modredundant("B 1 2 F\n" + text, n_atoms=3)

    def test_atom_range_is_not_checked_without_a_molecule(self):
        assert gi.parse_modredundant("B 1 99 F") == [("B", (1, 99), "F")]

    def test_helper_row_builds_a_line(self):
        assert gi.modredundant_line("Dihedral", "1, 2 3 4", "Scan", "36", "10.0") == "D 1 2 3 4 S 36 10.0"
        assert gi.modredundant_line("Bond", "1 2", "Freeze", "", "") == "B 1 2 F"

    def test_helper_row_refuses_what_the_parser_refuses(self):
        with pytest.raises(ValueError, match="needs 3 atom numbers"):
            gi.modredundant_line("Angle", "1 2", "Freeze", "", "")
```

Add `import re` to the imports of the test file.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_gaussianInput.py -q -k ModRedundant`
Expected: FAIL, `AttributeError: ... no attribute 'parse_modredundant'`.

- [ ] **Step 3: Implement** (append; add `import math` at the top of `gaussianInput.py`)

```python
# ----------------------------------------------------------------------
# Scan / freeze lines (Opt=ModRedundant)
# ----------------------------------------------------------------------
_MODRED_COUNTS = {letter: count for letter, count in gk.MODRED_TYPES.values()}


def parse_modredundant(text: str, n_atoms: Optional[int] = None) -> List[Tuple[str, Tuple[int, ...], str]]:
    """'B 1 2 S 10 0.05' / 'D 1 2 3 4 F' lines -> (letter, atoms, 'S'|'F'); ValueError names the line."""
    parsed = []
    for number, raw in enumerate((text or "").splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        where = f"Scan/freeze line {number} ('{line}')"
        tokens = line.split()
        letter = tokens[0].upper()
        if letter not in _MODRED_COUNTS:
            raise ValueError(f"{where}: start with B, A or D (bond, angle, dihedral).")
        count = _MODRED_COUNTS[letter]
        atom_tokens = tokens[1:1 + count]
        if len(atom_tokens) != count or not all(t.isdigit() for t in atom_tokens):
            raise ValueError(f"{where}: {letter} needs {count} atom numbers.")
        atoms = tuple(int(t) for t in atom_tokens)
        if len(set(atoms)) != count or min(atoms) < 1:
            raise ValueError(f"{where}: give {count} different atoms, numbered from 1.")
        if n_atoms is not None and max(atoms) > n_atoms:
            raise ValueError(f"{where}: atom {max(atoms)} does not exist (the molecule has {n_atoms} atoms).")
        rest = [t.upper() for t in tokens[1 + count:]]
        if rest == ["F"]:
            action = "F"
        elif len(rest) == 3 and rest[0] == "S":
            if not rest[1].isdigit() or int(rest[1]) < 1:
                raise ValueError(f"{where}: the number of steps must be a positive integer.")
            try:
                size = float(rest[2])
            except ValueError:
                size = 0.0
            if size == 0.0 or not math.isfinite(size):
                raise ValueError(f"{where}: the step size must be a non-zero number.")
            action = "S"
        else:
            raise ValueError(f"{where}: end with F (freeze) or S <steps> <step size> (scan).")
        parsed.append((letter, atoms, action))
    return parsed


def modredundant_line(kind: str, atoms_text: str, action: str, steps_text: str, size_text: str) -> str:
    """The Molecule panel's helper row -> one line, checked by the same rules as typed lines."""
    letter, _ = gk.MODRED_TYPES[kind]
    atoms = [t for t in re.split(r"[\s,]+", (atoms_text or "").strip()) if t]
    code = gk.MODRED_ACTIONS[action]
    words = [letter, *atoms, code] + ([steps_text.strip(), size_text.strip()] if code == "S" else [])
    line = " ".join(w for w in words if w)
    parse_modredundant(line)
    return line
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_gaussianInput.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add gaussianInput.py tests/test_gaussianInput.py
git commit -m "feat(gaussian): validated scan and freeze lines"
```

---

### Task 5: Job model, validation, input text, script and files

**Files:**
- Modify: `fukuiInput.py:425-449` (extract `formchk_script`), `gaussianInput.py` (append)
- Create: `tests/data/gaussian_input/phenol_formchk.sh`
- Test: `tests/test_fukuiInput.py` (append one test), `tests/test_gaussianInput.py` (append)

**Interfaces:**
- Consumes: Tasks 2–4; `fukuiInput.parse_memory`, `validate_stem`, `parse_geometry`, `format_geometry`, `build_route`, `_title_text`, `TITLE_FORBIDDEN`, `SCRIPT_SUFFIX`, `Atom`; `systemProperties.allowed_multiplicities`, `multiplicity_problem`, `nuclear_charge`; `molecularFileFormats.detect_bonds`.
- Produces:
  - `fukuiInput.formchk_script(script_name: str, header: Sequence[str], pairs: Sequence[Tuple[str, str]]) -> str`, `fukuiInput.FUKUI_SCRIPT_HEADER`.
  - `gaussianInput.GaussianJobSpec(stem, output_dir, geometry, steps, write_formchk)`, `PreparedStep`, `GaussianJob(spec, atoms, steps, input_path, script_path, chk_files, warnings, fukui)`, `INPUT_SUFFIX = ".gjf"`, `SCRIPT_HEADER`, `connectivity_lines(atoms) -> List[str]`, `linked_step(step) -> GaussianStepSpec`, `step_label(index: int, step) -> str` (1-based index), `GaussianInputBuilder` with `prepare(spec) -> GaussianJob`, `input_text(job) -> str`, `script_text(job) -> Optional[str]`, `write_files(job) -> List[Path]`. `prepare` calls `self._prepare_fukui(spec)` for a Fukui job; Task 6 adds that method (until then it raises `NotImplementedError`).

- [ ] **Step 1: Pin the Fukui formchk script before touching it**

```bash
mkdir tests\data\gaussian_input
python -c "import sys, tempfile; sys.path.insert(0, 'tests'); from pathlib import Path; import fukuiInput as fi; from test_fukuiInput import make_spec; b = fi.FukuiInputBuilder(); j = b.prepare(make_spec(Path(tempfile.mkdtemp()))); open('tests/data/gaussian_input/phenol_formchk.sh', 'w', newline='\n', encoding='ascii').write(b.script_text(j))"
```

Append to `tests/test_fukuiInput.py`:

```python
def test_fukui_script_is_unchanged_byte_for_byte(tmp_path):
    builder = fi.FukuiInputBuilder()
    job = builder.prepare(make_spec(tmp_path))
    pinned = Path(__file__).resolve().parent / "data" / "gaussian_input" / "phenol_formchk.sh"
    assert builder.script_text(job) == pinned.read_bytes().decode("ascii")
```

Run: `python -m pytest tests/test_fukuiInput.py -q` — Expected: all pass (the pin matches the current code).

- [ ] **Step 2: Extract `formchk_script` in `fukuiInput.py`**

Add above `class FukuiInputBuilder`:

```python
FUKUI_SCRIPT_HEADER = ("# gQTEA Fukui input builder: convert the three checkpoints to formatted",
                       "# checkpoints for Tools > Analysis of Fukui Functions. Run it where Gaussian")


def formchk_script(script_name: str, header: Sequence[str], pairs: Sequence[Tuple[str, str]]) -> str:
    """A bash script that runs formchk on each (chk, fchk) pair, run where Gaussian ran."""
    lines = ["#!/bin/bash", *header,
             f"# ran (WSL or the cluster), after the jobs finish:  bash {script_name}",
             'cd "$(dirname "$0")" || exit 1',
             "if ! command -v formchk >/dev/null 2>&1; then",
             "  echo \"formchk not found: set up Gaussian 16 first (e.g. 'module load gaussian', or source",
             "  \\$g16root/g16/bsd/g16.profile).\" >&2",
             "  exit 1",
             "fi",
             "status=0",
             "convert() {",
             '  if [ ! -f "$1" ]; then',
             '    echo "$1 not found: did that job run? Gaussian stops the remaining linked jobs when one '
             'fails; see the .log file." >&2',
             "    status=1",
             "    return",
             "  fi",
             '  formchk "$1" "$2" || { echo "formchk failed for $1" >&2; status=1; }',
             "}"]
    lines += [f'convert "{chk}" "{fchk}"' for chk, fchk in pairs]
    lines.append("exit $status")
    return "\n".join(lines) + "\n"
```

Replace the body of `FukuiInputBuilder.script_text` with:

```python
    def script_text(self, job: FukuiJob) -> str:
        """A bash script that converts the three checkpoints, run where Gaussian ran."""
        return formchk_script(job.script_path.name, FUKUI_SCRIPT_HEADER,
                              [(s.chk, s.fchk.name) for s in job.states])
```

Run: `python -m pytest tests/test_fukuiInput.py -q` — Expected: all pass, including the byte-for-byte pin.

- [ ] **Step 3: Write the failing job tests** (append to `tests/test_gaussianInput.py`)

```python
WATER = ("O 0.000000 0.000000 0.117300\n"
         "H 0.000000 0.757200 -0.469200\n"
         "H 0.000000 -0.757200 -0.469200")
WATER_LINES = fi.format_geometry(fi.parse_geometry(WATER))


def job_spec(tmp_path, *steps, **changes) -> gi.GaussianJobSpec:
    steps = list(steps) or [step(multiplicity=1)]
    return replace(gi.GaussianJobSpec(stem="water", output_dir=str(tmp_path), geometry=WATER, steps=steps),
                   **changes)


def problems(tmp_path, *steps, **changes) -> str:
    with pytest.raises(ValueError) as info:
        gi.GaussianInputBuilder().prepare(job_spec(tmp_path, *steps, **changes))
    return str(info.value)


class TestInputText:
    def test_single_step(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path, step(job_type=gk.JOB_OPT_FREQ, multiplicity=1)))
        assert builder.input_text(job) == "\n".join(
            ["%Chk=water.chk", "%Mem=2GB", "%NProcShared=4", "#P B3LYP/6-31G(d) Opt Freq", "",
             "water", "", "0 1", *WATER_LINES, ""]) + "\n"
        assert job.input_path == tmp_path / "water.gjf"
        assert job.script_path == tmp_path / "water_formchk.sh" and job.chk_files == ["water.chk"]

    def test_linked_opt_then_td(self, tmp_path):
        first = step(job_type=gk.JOB_OPT, multiplicity=1)
        second = replace(gi.linked_step(first), state=gk.STATE_TD, n_states=10, print_all_coefficients=True)
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path, first, second))
        assert builder.input_text(job) == "\n".join(
            ["%Chk=water.chk", "%Mem=2GB", "%NProcShared=4", "#P B3LYP/6-31G(d) Opt", "",
             "water", "", "0 1", *WATER_LINES, "",
             "--Link1--",
             "%Chk=water.chk", "%Mem=2GB", "%NProcShared=4",
             "#P B3LYP/6-31G(d) TD=(NStates=10) IOp(9/40=4) Geom=AllCheck Guess=Read", ""]) + "\n"
        assert any("same checkpoint" in w for w in job.warnings)
        assert any("Density=Current" in w for w in job.warnings)

    def test_check_step_writes_title_and_charge_but_no_coordinates(self, tmp_path):
        second = step(geometry_source=gk.GEOM_CHECK, charge=1, multiplicity=2, guess="Read checkpoint")
        builder = gi.GaussianInputBuilder()
        text = builder.input_text(builder.prepare(job_spec(tmp_path, step(multiplicity=1), second)))
        block = text.split("--Link1--\n")[1]
        assert block.endswith("Geom=Check Guess=Read\n\nwater\n\n1 2\n\n")

    def test_sections_after_the_geometry(self, tmp_path):
        s = step(job_type=gk.JOB_SCAN, multiplicity=1, modredundant="B 1 2 S 5 0.05", connectivity=True,
                 basis="Gen", additional_input="O H 0\n6-31G(d)\n****")
        builder = gi.GaussianInputBuilder()
        text = builder.input_text(builder.prepare(job_spec(tmp_path, s)))
        assert text.endswith("\n".join([*WATER_LINES, "", " 1 2 1.0 3 1.0", " 2", " 3", "",
                                        "B 1 2 S 5 0.05", "", "O H 0", "6-31G(d)", "****", ""]) + "\n")
        assert "Geom=Connectivity" in text and "Opt=ModRedundant" in text

    def test_typed_title_and_link0(self, tmp_path):
        s = step(multiplicity=1, title="my water", chk="w1.chk", old_chk="start.chk", memory="8gb", nproc=16,
                 extra_link0="%KJob L502")
        builder = gi.GaussianInputBuilder()
        text = builder.input_text(builder.prepare(job_spec(tmp_path, s)))
        assert text.startswith("%OldChk=start.chk\n%Chk=w1.chk\n%Mem=8GB\n%NProcShared=16\n%KJob L502\n")
        assert "\nmy water\n" in text

    def test_no_checkpoint_means_no_script(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, step(multiplicity=1, write_chk=False)))
        assert job.script_path is None and any("no formchk script" in w for w in job.warnings)

    def test_files_are_written_with_lf_endings(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path))
        assert builder.write_files(job) == [job.input_path, job.script_path]
        assert b"\r" not in job.input_path.read_bytes()
        script = job.script_path.read_text()
        assert 'convert "water.chk" "water.fchk"' in script and "3D Molecular Viewer" in script

    def test_write_error_is_a_value_error(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path))
        job.input_path = tmp_path / "missing" / "x.gjf"
        with pytest.raises(ValueError, match="Could not write"):
            builder.write_files(job)


class TestValidation:
    @pytest.mark.parametrize("changes,fragment", [
        (dict(multiplicity=None), "Step 1: Multiplicity: enter it; it is not assumed (10 electrons allow 1, 3, 5, ...)"),
        (dict(multiplicity=2), "Multiplicity 2 is impossible"),
        (dict(multiplicity=3, spin="Restricted (R)"), "R (restricted) prefix needs multiplicity 1"),
        (dict(multiplicity=1, method_name=gk.OTHER, method_other="B3 LYP"), "Other method: enter one keyword"),
        (dict(multiplicity=1, basis=gk.OTHER, basis_other=""), "Other basis set: enter one keyword"),
        (dict(multiplicity=1, basis="GenECP"), "Basis set GenECP: put the basis-set section in the Add. Inp. panel"),
        (dict(multiplicity=1, state=gk.STATE_TD, n_states=0), "Number of states must be at least 1"),
        (dict(multiplicity=1, state=gk.STATE_TD, n_states=3, root=4), "Root) must be between 1 and 3"),
        (dict(multiplicity=1, state=gk.STATE_TD, family=gk.FAMILY_MP2, method_name=gk.NO_CHOICE),
         "TD-SCF needs a Hartree-Fock or DFT method"),
        (dict(multiplicity=1, state=gk.STATE_TDA, family=gk.FAMILY_HF, method_name=gk.NO_CHOICE), "TDA needs a DFT method"),
        (dict(multiplicity=1, state=gk.STATE_TD, state_spin="Singlets", spin="Unrestricted (U)"),
         "need a closed-shell reference"),
        (dict(multiplicity=1, print_all_coefficients=True), "needs an excited state"),
        (dict(multiplicity=1, title="a#b"), "Title: Gaussian does not allow #"),
        (dict(multiplicity=1, chk="C:\\jobs\\w.chk"), "%Chk 'C:\\jobs\\w.chk': give a file name relative"),
        (dict(multiplicity=1, chk="my w.chk"), "%Chk 'my w.chk': use a file name without blanks"),
        (dict(multiplicity=1, extra_link0="KJob L502"), "must start with %"),
        (dict(multiplicity=1, memory="2"), "%Mem = '2'"),
        (dict(multiplicity=1, guess="Read checkpoint"), "Guess=Read reads a checkpoint, but no earlier step writes"),
        (dict(multiplicity=1, geometry_source=gk.GEOM_ALLCHECK), "Geom=AllCheck reads a checkpoint"),
        (dict(multiplicity=1, job_type=gk.JOB_OPT, force_constants="Read from checkpoint (ReadFC)"),
         "Opt=ReadFC reads a checkpoint"),
        (dict(multiplicity=1, geometry_source=gk.GEOM_CHECK, connectivity=True, old_chk="a.chk"),
         "Write connectivity needs coordinates"),
        (dict(multiplicity=1, modredundant="B 1 2 F"), "Scan/freeze lines need the Optimization"),
        (dict(multiplicity=1, job_type=gk.JOB_OPT, modredundant="B 1 2 S 5 0.1"), "Scan lines (S) need the Scan job type"),
        (dict(multiplicity=1, job_type=gk.JOB_SCAN), "The Scan job type needs at least one scan line"),
        (dict(multiplicity=1, job_type=gk.JOB_OPT, modredundant="B 1 7 F"), "atom 7 does not exist"),
        (dict(multiplicity=1, solvation="SMD", solvent=gk.OTHER, solvent_other="n,n"), "Other solvent"),
        (dict(multiplicity=1, additional_keywords="B3LYP/STO-3G"), "Method panel"),
        (dict(multiplicity=1, job_type="Dance"), "Job type: 'Dance' is not one of the choices"),
    ])
    def test_step_problem(self, tmp_path, changes, fragment):
        assert fragment in problems(tmp_path, step(**changes))

    def test_job_problems(self, tmp_path):
        text = problems(tmp_path, step(multiplicity=1), stem="bad name", output_dir=str(tmp_path / "nope"),
                        geometry="")
        assert "Job name 'bad name'" in text and "does not exist" in text and "geometry is empty" in text

    def test_every_problem_is_listed_on_its_own_line(self, tmp_path):
        lines = problems(tmp_path, step(multiplicity=None, memory="x"), step(multiplicity=None)).splitlines()
        assert len(lines) >= 3 and all(line.startswith("- ") for line in lines)
        assert any(line.startswith("- Step 2: ") for line in lines)

    def test_an_old_checkpoint_lets_step_one_read(self, tmp_path):
        s = step(geometry_source=gk.GEOM_ALLCHECK, old_chk="start.chk")
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, s, geometry=""))
        assert job.atoms is None

    def test_a_later_step_reads_the_checkpoint_an_earlier_step_writes(self, tmp_path):
        first = step(multiplicity=1, chk="a.chk")
        second = step(geometry_source=gk.GEOM_ALLCHECK, chk="a.chk", guess="Read checkpoint")
        gi.GaussianInputBuilder().prepare(job_spec(tmp_path, first, second))
        with pytest.raises(ValueError, match="no earlier step writes b.chk"):
            gi.GaussianInputBuilder().prepare(job_spec(tmp_path, first, replace(second, chk="b.chk")))

    def test_ts_and_irc_without_force_constants_warn(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        ts = builder.prepare(job_spec(tmp_path, step(multiplicity=1, job_type=gk.JOB_OPT, opt_target=gk.OPT_TS)))
        irc = builder.prepare(job_spec(tmp_path, step(multiplicity=1, job_type=gk.JOB_IRC,
                                                      irc_force_constants=gk.FC_DEFAULT)))
        assert any("TS" in w for w in ts.warnings) and any("IRC" in w for w in irc.warnings)

    def test_default_job_has_no_warnings(self, tmp_path):
        assert gi.GaussianInputBuilder().prepare(job_spec(tmp_path)).warnings == []


class TestSteps:
    def test_linked_step_continues_from_the_checkpoint(self):
        first = step(job_type=gk.JOB_OPT_FREQ, multiplicity=1, title="t", modredundant="B 1 2 F",
                     additional_input="x", memory="8GB", chk="a.chk", old_chk="o.chk")
        new = gi.linked_step(first)
        assert (new.job_type, new.geometry_source, new.guess) == (gk.JOB_ENERGY, gk.GEOM_ALLCHECK, "Read checkpoint")
        assert (new.title, new.modredundant, new.additional_input, new.old_chk) == ("", "", "", "")
        assert (new.memory, new.chk, new.family, new.method_name) == ("8GB", "a.chk", gk.FAMILY_DFT, "B3LYP")
        assert (new.charge, new.multiplicity) == (0, None)

    def test_step_label(self):
        assert gi.step_label(1, step(job_type=gk.JOB_OPT)) == "1: Optimization B3LYP/6-31G(d)"
        assert gi.step_label(2, step(geometry_source=gk.GEOM_ALLCHECK, state=gk.STATE_TD)) == \
            "2: TD-SCF B3LYP/6-31G(d), Geom=AllCheck"
```

- [ ] **Step 4: Run the tests to verify they fail**

Run: `python -m pytest tests/test_gaussianInput.py -q`
Expected: the new tests fail (`AttributeError: ... 'GaussianJobSpec'`).

- [ ] **Step 5: Implement** (append to `gaussianInput.py`; add `import dataclasses`, `from pathlib import Path, PureWindowsPath` and `from systemProperties import allowed_multiplicities, multiplicity_problem, nuclear_charge` to the imports)

```python
# ----------------------------------------------------------------------
# The job
# ----------------------------------------------------------------------
INPUT_SUFFIX = ".gjf"
SCRIPT_HEADER = ("# gQTEA Gaussian 16 input builder: convert the checkpoints to formatted",
                 "# checkpoints for the 3D Molecular Viewer's analysis tabs. Run it where Gaussian")
FORCE_CONSTANT_OPTIONS = {"CALCFC", "CALCALL", "READFC", "RCFC"}


@dataclass
class GaussianJobSpec:
    """The whole job: file names, the step-1 geometry and one GaussianStepSpec per linked step."""
    stem: str = fi.DEFAULT_STEM
    output_dir: str = ""
    geometry: str = ""
    steps: List[GaussianStepSpec] = field(default_factory=lambda: [GaussianStepSpec()])
    write_formchk: bool = True


@dataclass
class PreparedStep:
    """One validated step, ready to write."""
    link0: List[str]
    route: List[str]
    title: Optional[str]            # None with Geom=AllCheck (no title section)
    charge_line: Optional[str]      # None with Geom=AllCheck
    coordinates: bool
    connectivity: List[str]
    modredundant: List[str]
    additional_input: List[str]
    chk: Optional[str]


@dataclass
class GaussianJob:
    spec: GaussianJobSpec
    atoms: Optional[List[fi.Atom]]
    steps: List[PreparedStep]
    input_path: Path
    script_path: Optional[Path]
    chk_files: List[str]
    warnings: List[str] = field(default_factory=list)
    fukui: Optional[fi.FukuiJob] = None


def _relative_name(value: str, label: str) -> str:
    name = (value or "").strip()
    if not name or any(ch.isspace() for ch in name) or not name.isascii():
        raise ValueError(f"{label} '{name}': use a file name without blanks (ASCII letters, digits, . _ -).")
    if PureWindowsPath(name).drive or "/" in name or "\\" in name:
        raise ValueError(f"{label} '{name}': give a file name relative to the input's folder "
                         "(the job runs under Linux from there).")
    return name


def _one_keyword(text: str) -> bool:
    text = (text or "").strip()
    return bool(text) and len(text.split()) == 1 and "/" not in text and "," not in text


def _checkpoint_readers(keywords: Sequence[RouteKeyword]) -> List[str]:
    """What in the route reads the checkpoint: 'Geom=AllCheck', 'Guess=Read', 'Opt=ReadFC', ..."""
    readers = []
    for keyword in keywords:
        options = {o.upper() for o in keyword.options}
        if keyword.key == "GEOM":
            readers += [f"Geom={o}" for o in ("AllCheck", "Check") if o.upper() in options]
        elif keyword.key == "GUESS":
            readers += [f"Guess={o}" for o in gk.READ_GUESSES if o.upper() in options]
        elif keyword.key in ("OPT", "IRC"):
            readers += [f"{keyword.name}={o}" for o in gk.READ_FORCE_CONSTANTS if o.upper() in options]
    return readers


def connectivity_lines(atoms: Sequence[fi.Atom]) -> List[str]:
    """The Geom=Connectivity section, single bonds from the distances (the 3D viewer's rule)."""
    from molecularFileFormats import detect_bonds    # numpy + scipy, needed only for this option
    bonds = detect_bonds([a[0] for a in atoms], [a[1:] for a in atoms])
    partners: Dict[int, List[int]] = {i: [] for i in range(len(atoms))}
    for i, j, _ in bonds:
        partners[i].append(j)
    return [" ".join([f" {i + 1}"] + [f"{j + 1} 1.0" for j in sorted(partners[i])]) for i in range(len(atoms))]


def _trim_blank_lines(text: str) -> List[str]:
    lines = [line.rstrip() for line in (text or "").splitlines()]
    while lines and not lines[0]:
        lines.pop(0)
    while lines and not lines[-1]:
        lines.pop()
    return lines


def _prepare_step(step: GaussianStepSpec, stem: str, atoms, written: set,
                  write_formchk: bool) -> Tuple[List[str], List[str], Optional[PreparedStep]]:
    """(errors, warnings, prepared step or None) of one step; messages carry no step prefix."""
    errors: List[str] = []
    warnings: List[str] = []

    def attempt(fn, *args):
        try:
            return fn(*args)
        except ValueError as exc:
            errors.append(str(exc))
            return None

    for name, values in CHOICES.items():
        if getattr(step, name) not in values:
            errors.append(f"{FIELD_LABELS[name]}: '{getattr(step, name)}' is not one of the choices.")
    if step.family in gk.FAMILIES and step.method_name not in method_choices(step.family):
        errors.append(f"Method: '{step.method_name}' is not a {step.family} method.")
    if errors:
        return errors, warnings, None

    # Link 0
    memory = attempt(fi.parse_memory, step.memory)
    if not isinstance(step.nproc, int) or step.nproc < 1:
        errors.append(f"%NProcShared must be a positive integer (got {step.nproc}).")
    chk = attempt(_relative_name, step.chk.strip() or f"{stem}.chk", "%Chk") if step.write_chk else None
    old_chk = attempt(_relative_name, step.old_chk, "%OldChk") if step.old_chk.strip() else None
    extra_link0 = [line.strip() for line in step.extra_link0.splitlines() if line.strip()]
    errors += [f"Extra Link 0 line '{line}' must start with %." for line in extra_link0 if not line.startswith("%")]

    # Method and basis
    excited = step.state != gk.STATE_GROUND
    if (step.state not in gk.MODEL_STATES and step.family == gk.FAMILY_DFT and step.method_name == gk.OTHER
            and not _one_keyword(step.method_other)):
        errors.append("Other method: enter one keyword, e.g. B97D3 (no blanks, no '/').")
    if not no_basis(step):
        if step.basis == gk.OTHER and not _one_keyword(step.basis_other):
            errors.append("Other basis set: enter one keyword, e.g. def2TZVPD (no blanks, no '/').")
        if step.basis in gk.GENERAL_BASES and not step.additional_input.strip():
            errors.append(f"Basis set {step.basis}: put the basis-set section in the Add. Inp. panel.")
    if step.solvation != gk.SOLVATION_NONE and step.solvent == gk.OTHER and not _one_keyword(step.solvent_other):
        errors.append("Other solvent: enter one Gaussian solvent name, e.g. Heptane.")
    if excited:
        if not isinstance(step.n_states, int) or step.n_states < 1:
            errors.append("Number of states must be at least 1.")
        elif not isinstance(step.root, int) or not 1 <= step.root <= step.n_states:
            errors.append(f"State of interest (Root) must be between 1 and {step.n_states}.")
        if step.state == gk.STATE_TD and step.family not in (gk.FAMILY_HF, gk.FAMILY_DFT):
            errors.append("TD-SCF needs a Hartree-Fock or DFT method.")
        if step.state == gk.STATE_TDA and step.family != gk.FAMILY_DFT:
            errors.append("TDA needs a DFT method.")
    elif step.print_all_coefficients:
        errors.append("Print all excitation coefficients needs an excited state (TD-SCF, TDA, CIS or ZINDO).")

    # Geometry, charge, multiplicity
    source = step.geometry_source
    prefix = gk.SPIN_PREFIXES[step.spin]
    mult = step.multiplicity
    if step.connectivity and source != gk.GEOM_COORDINATES:
        errors.append("Write connectivity needs coordinates; this step reads its geometry from the checkpoint.")
    if source != gk.GEOM_ALLCHECK:
        electrons = None if atoms is None else nuclear_charge(a[0] for a in atoms) - step.charge
        if mult is None:
            allowed = (f" ({electrons} electrons allow "
                       f"{', '.join(str(m) for m in allowed_multiplicities(electrons, 3))}, ...)"
                       if electrons is not None and electrons >= 1 else "")
            errors.append(f"Multiplicity: enter it; it is not assumed{allowed}.")
        elif electrons is not None:
            problem = multiplicity_problem(electrons, mult) if electrons >= 1 else "no electrons are left"
            if problem:
                errors.append(f"Multiplicity {mult} is impossible ({problem}).")
    if mult is not None and mult > 1 and prefix == "R":
        errors.append(f"The R (restricted) prefix needs multiplicity 1 (got {mult}); use RO or U.")
    if gk.STATE_SPINS.get(step.state_spin) and (prefix in ("U", "RO") or (mult is not None and mult != 1)):
        errors.append(f"States = {step.state_spin}: Singlets, Triplets and 50-50 need a closed-shell reference "
                      "(multiplicity 1, no U or RO prefix).")
    title = step.title.strip() or fi._title_text(stem)
    bad = sorted(set(title) & set(fi.TITLE_FORBIDDEN))
    if bad:
        errors.append(f"Title: Gaussian does not allow {' '.join(bad)} in the title.")

    # Scan / freeze lines
    parsed: List[Tuple[str, Tuple[int, ...], str]] = []
    modred_ok = True
    if step.modredundant.strip():
        if step.job_type not in gk.OPTIMIZING_JOBS:
            errors.append("Scan/freeze lines need the Optimization, Opt+Freq or Scan job type.")
            modred_ok = False
        else:
            parsed = attempt(parse_modredundant, step.modredundant, len(atoms) if atoms else None)
            modred_ok = parsed is not None
            parsed = parsed or []
            if any(p[2] == "S" for p in parsed) and step.job_type != gk.JOB_SCAN:
                errors.append("Scan lines (S) need the Scan job type; Optimization takes freeze lines (F) only.")
    if step.job_type == gk.JOB_SCAN and modred_ok and not any(p[2] == "S" for p in parsed):
        errors.append("The Scan job type needs at least one scan line (S); add one in the Molecule panel.")

    # Route and checkpoints
    keywords, merge_errors = route_keywords(step)
    errors += merge_errors
    readers = _checkpoint_readers(keywords)
    if readers and old_chk is None and (chk is None or chk.lower() not in written):
        errors.append(f"{', '.join(readers)} reads a checkpoint, but no earlier step writes "
                      f"{chk or 'one'}; set %OldChk to an existing checkpoint.")
    options = {k.key: {o.upper() for o in k.options} for k in keywords}
    if "TS" in options.get("OPT", set()) and not FORCE_CONSTANT_OPTIONS & options["OPT"]:
        warnings.append("Optimization to a TS without force constants (CalcFC, CalcAll or ReadFC) often fails; "
                        "CalcFC is recommended.")
    if "IRC" in options and not FORCE_CONSTANT_OPTIONS & options["IRC"]:
        warnings.append("IRC without force constants: give CalcFC, CalcAll or ReadFC.")
    if write_formchk and (excited or step.family in gk.POST_SCF_FAMILIES) and "DENSITY" not in options:
        warnings.append("The .fchk will hold the SCF density, not the excited-state or correlated one; add "
                        "Density=Current to the Additional Keywords if an analysis needs it.")
    if errors:
        return errors, warnings, None

    link0 = ([f"%OldChk={old_chk}"] if old_chk else []) + ([f"%Chk={chk}"] if chk else []) + \
        [f"%Mem={memory}", f"%NProcShared={step.nproc}"] + extra_link0
    route = fi.build_route(step.print_level, method_keyword(step), basis_keyword(step),
                           " ".join(format_keyword(k) for k in keywords))
    all_check = source == gk.GEOM_ALLCHECK
    return errors, warnings, PreparedStep(
        link0=link0, route=route, title=None if all_check else title,
        charge_line=None if all_check else f"{step.charge} {mult}",
        coordinates=source == gk.GEOM_COORDINATES,
        connectivity=connectivity_lines(atoms) if step.connectivity else [],
        modredundant=[line.strip() for line in step.modredundant.splitlines() if line.strip()],
        additional_input=_trim_blank_lines(step.additional_input), chk=chk)


def linked_step(step: GaussianStepSpec) -> GaussianStepSpec:
    """A new step continuing from ``step``: same method and Link 0, geometry and orbitals from its checkpoint.

    Charge and multiplicity are cleared: Geom=AllCheck reads them from the checkpoint, and a step
    switched to another source must be given them (they are never assumed).
    """
    return dataclasses.replace(step, job_type=gk.JOB_ENERGY, geometry_source=gk.GEOM_ALLCHECK,
                               guess="Read checkpoint", title="", modredundant="", additional_input="",
                               connectivity=False, old_chk="", charge=0, multiplicity=None)


def step_label(index: int, step: GaussianStepSpec) -> str:
    """'2: TD-SCF B3LYP/6-31G(d), Geom=AllCheck' for the step list (index is 1-based)."""
    what = step.state if step.state != gk.STATE_GROUND and step.job_type == gk.JOB_ENERGY else step.job_type
    source = gk.GEOMETRY_SOURCES.get(step.geometry_source)
    return f"{index}: {what} {model_token(step)}" + (f", Geom={source}" if source else "")


class GaussianInputBuilder:
    """Validation, the input text, the formchk script and the files."""

    def __init__(self):
        self._fukui = fi.FukuiInputBuilder()

    def prepare(self, spec: GaussianJobSpec) -> GaussianJob:
        """Validate everything; ValueError lists every problem, one '- ' line each."""
        if any(s.job_type == gk.JOB_FUKUI for s in spec.steps):
            return self._prepare_fukui(spec)
        errors: List[str] = []
        warnings: List[str] = []
        stem = None
        try:
            stem = fi.validate_stem(spec.stem)
        except ValueError as exc:
            errors.append(str(exc))
        folder_text = (spec.output_dir or "").strip()
        folder = Path(folder_text).resolve() if folder_text else None
        if folder is None or not folder.is_dir():
            errors.append(f"Output folder '{folder_text}' does not exist; choose an existing folder."
                          if folder_text else "Output folder: choose where the input files are written.")
        if not spec.steps:
            errors.append("The job has no steps.")
        atoms = None
        if any(s.geometry_source == gk.GEOM_COORDINATES for s in spec.steps):
            try:
                atoms = fi.parse_geometry(spec.geometry)
            except ValueError as exc:
                errors.append(str(exc))
        name = stem or fi.DEFAULT_STEM
        written: set = set()
        prepared: List[PreparedStep] = []
        for index, step in enumerate(spec.steps, start=1):
            step_errors, step_warnings, item = _prepare_step(step, name, atoms, written, spec.write_formchk)
            errors += [f"Step {index}: {e}" for e in step_errors]
            warnings += [f"Step {index}: {w}" for w in step_warnings]
            if item is not None:
                prepared.append(item)
            if step.write_chk:
                written.add((step.chk.strip() or f"{name}.chk").lower())
        if errors:
            raise ValueError("\n".join(f"- {e}" for e in errors))
        chk_files = list(dict.fromkeys(p.chk for p in prepared if p.chk))
        if len([p for p in prepared if p.chk]) > len(chk_files):
            warnings.append("Several steps write the same checkpoint: it keeps only the last step's results, "
                            "and that is what the formchk script converts.")
        script_path = None
        if spec.write_formchk and chk_files:
            script_path = folder / f"{stem}{fi.SCRIPT_SUFFIX}"
        elif spec.write_formchk:
            warnings.append("No step writes a checkpoint (%Chk), so no formchk script is written.")
        return GaussianJob(spec, atoms, prepared, folder / f"{stem}{INPUT_SUFFIX}", script_path, chk_files, warnings)

    def _prepare_fukui(self, spec: GaussianJobSpec) -> GaussianJob:
        raise NotImplementedError("Task 6")

    def input_text(self, job: GaussianJob) -> str:
        if job.fukui is not None:
            return self._fukui.input_text(job.fukui)
        blocks = []
        for step in job.steps:
            lines = list(step.link0) + list(step.route) + [""]
            if step.title is not None:
                lines += [step.title, "", step.charge_line]
                if step.coordinates:
                    lines += fi.format_geometry(job.atoms)
                lines.append("")
            for section in (step.connectivity, step.modredundant, step.additional_input):
                if section:
                    lines += list(section) + [""]
            blocks.append("\n".join(lines))
        return "\n--Link1--\n".join(blocks) + "\n"

    def script_text(self, job: GaussianJob) -> Optional[str]:
        if job.fukui is not None:
            return self._fukui.script_text(job.fukui)
        if job.script_path is None:
            return None
        return fi.formchk_script(job.script_path.name, SCRIPT_HEADER,
                                 [(chk, Path(chk).stem + fi.FCHK_EXTENSION) for chk in job.chk_files])

    def write_files(self, job: GaussianJob) -> List[Path]:
        """Write the .gjf (and the .sh) with LF line endings; ValueError on any OS error."""
        if job.fukui is not None:
            return self._fukui.write_files(job.fukui)
        written = []
        files = [(job.input_path, self.input_text(job))]
        if job.script_path is not None:
            files.append((job.script_path, self.script_text(job)))
        for path, text in files:
            try:
                with open(path, "w", encoding="ascii", errors="strict", newline="\n") as handle:
                    handle.write(text)
            except (OSError, UnicodeEncodeError) as exc:
                raise ValueError(f"Could not write {path}: {exc}") from exc
            written.append(path)
        return written
```

Check against the tests while implementing:
- `test_check_step_writes_title_and_charge_but_no_coordinates` expects the block to end with `"...Guess=Read\n\nwater\n\n1 2\n\n"`: lines = route, `""`, title, `""`, `"1 2"`, `""` → joined, then the final `"\n"`. Water with charge +1 has 9 electrons, so multiplicity 2 is valid.
- `test_sections_after_the_geometry`: the connectivity of water is O bonded to both H (` 1 2 1.0 3 1.0`, ` 2`, ` 3`).

- [ ] **Step 6: Run the tests to verify they pass**

Run: `python -m pytest tests/test_gaussianInput.py tests/test_fukuiInput.py -q`
Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add fukuiInput.py gaussianInput.py tests/test_fukuiInput.py tests/test_gaussianInput.py tests/data/gaussian_input/phenol_formchk.sh
git commit -m "feat(gaussian): job validation, linked steps, input text and formchk script"
```

---

### Task 6: Fukui Job Type and live preview texts

**Files:**
- Modify: `gaussianInput.py` (append; replace the `_prepare_fukui` stub)
- Test: `tests/test_gaussianInput.py` (append)

**Interfaces:**
- Consumes: Tasks 3 and 5; `fukuiInput.FukuiJobSpec`, `FukuiInputBuilder`, `state_tags`, `build_route`, `route_tokens`.
- Produces: `FUKUI_FIXED: Tuple[str, ...]`; `fukui_defaults(step) -> GaussianStepSpec`; `fukui_opt_options(step) -> List[str]`; `fukui_keywords(step) -> Tuple[List[RouteKeyword], List[str]]`; `fukui_spec(spec: GaussianJobSpec) -> fi.FukuiJobSpec`; `GaussianInputBuilder._prepare_fukui(spec) -> GaussianJob` (with `fukui` set); `route_preview(spec, index: int) -> str` (0-based index); `molecule_summary(spec, index: int) -> str`.

- [ ] **Step 1: Write the failing tests** (append)

```python
PHENOL_GJF = SRC / "tests" / "data" / "fukui" / "inputs" / "phenol_N.gjf"


def phenol_geometry() -> str:
    lines = PHENOL_GJF.read_text().splitlines()
    start = lines.index("0 1") + 1
    return "\n".join(line for line in lines[start:] if line.strip())


def fukui_job_spec(tmp_path, **changes) -> gi.GaussianJobSpec:
    s = gi.fukui_defaults(step(multiplicity=1, fukui_plus_multiplicity=2, fukui_minus_multiplicity=2))
    return gi.GaussianJobSpec(stem="phenol", output_dir=str(tmp_path), geometry=phenol_geometry(),
                              steps=[replace(s, **changes)])


class TestFukui:
    def test_output_is_the_fukui_builders_byte_for_byte(self, tmp_path):
        old = fi.FukuiInputBuilder()
        old_job = old.prepare(fi.FukuiJobSpec(stem="phenol", output_dir=str(tmp_path), geometry=phenol_geometry(),
                                              charge=0, multiplicity=1, plus_multiplicity=2, minus_multiplicity=2,
                                              memory="2GB", nproc=4))
        new = gi.GaussianInputBuilder()
        job = new.prepare(fukui_job_spec(tmp_path))
        assert new.input_text(job) == old.input_text(old_job)
        assert new.script_text(job) == old.script_text(old_job)
        assert (job.input_path, job.script_path) == (old_job.input_path, old_job.script_path)
        assert new.write_files(job) == [old_job.input_path, old_job.script_path]

    def test_panels_map_onto_the_fukui_spec(self, tmp_path):
        spec = fukui_job_spec(tmp_path, method_name="M06-2X", basis="Def2SVP", dispersion="GD3BJ", pop_mk=True,
                              solvation="SMD", force_constants=gk.FC_CALCFC, opt_tight=True, opt_max_cycles=50,
                              additional_keywords="Int=UltraFine", print_level="#N", memory="3GB", nproc=8,
                              charge=1, multiplicity=2, fukui_p=2, fukui_q=1, fukui_plus_multiplicity=3,
                              fukui_minus_multiplicity=1, fukui_guess_read=False)
        f = gi.fukui_spec(spec)
        assert (f.method, f.basis, f.print_level, f.memory, f.nproc) == ("M062X", "Def2SVP", "#N", "3GB", 8)
        assert f.keywords == "EmpiricalDispersion=GD3BJ SCF=Tight NoSymm Pop=MK SCRF=(SMD,Solvent=Water) Int=UltraFine"
        assert f.opt_keyword == "Opt=(CalcFC,Tight,MaxCycles=50)"
        assert (f.charge, f.multiplicity, f.plus_multiplicity, f.minus_multiplicity) == (1, 2, 3, 1)
        assert (f.p, f.q, f.guess_read, f.stem, f.geometry) == (2, 1, False, "phenol", phenol_geometry())

    def test_fukui_must_be_the_only_step(self, tmp_path):
        spec = fukui_job_spec(tmp_path)
        spec.steps.append(step(multiplicity=1))
        with pytest.raises(ValueError, match="must be the only step"):
            gi.GaussianInputBuilder().prepare(spec)

    @pytest.mark.parametrize("changes,label", [
        (dict(force_constants="Read from checkpoint (ReadFC)"), "Force constants"),
        (dict(guess="Read checkpoint"), "Initial guess"),
        (dict(state=gk.STATE_TD), "Electronic state"),
        (dict(title="x"), "Title"),
    ])
    def test_controls_the_fukui_job_sets_itself_are_refused(self, tmp_path, changes, label):
        with pytest.raises(ValueError, match=f"sets {label} itself"):
            gi.GaussianInputBuilder().prepare(fukui_job_spec(tmp_path, **changes))

    def test_additional_keywords_merge_and_opt_is_refused_with_the_new_wording(self, tmp_path):
        assert gi.fukui_spec(fukui_job_spec(tmp_path, additional_keywords="SCF=XQC")).keywords == "SCF=(Tight,XQC) NoSymm"
        with pytest.raises(ValueError, match="Optimization options of the Job Type panel"):
            gi.GaussianInputBuilder().prepare(fukui_job_spec(tmp_path, additional_keywords="Opt=Tight"))

    def test_fukui_problems_come_from_fukui_input(self, tmp_path):
        with pytest.raises(ValueError, match="N\\+1 state: enter its spin multiplicity"):
            gi.GaussianInputBuilder().prepare(fukui_job_spec(tmp_path, fukui_plus_multiplicity=None))

    def test_fukui_defaults(self):
        s = gi.fukui_defaults(step(guess="Read checkpoint", title="t", spin="Unrestricted (U)", chk="a.chk",
                                   force_constants="Read from checkpoint (ReadFC)"))
        assert s.job_type == gk.JOB_FUKUI and s.scf_tight and s.nosymm
        assert (s.guess, s.title, s.spin, s.chk, s.force_constants) == ("Default", "", "Default", "", gk.FC_DEFAULT)


class TestPreview:
    def test_fukui_route_preview_matches_the_former_window(self, tmp_path):
        assert gi.route_preview(fukui_job_spec(tmp_path), 0) == (
            "job 1:     #P B3LYP/6-31G(d) SCF=Tight NoSymm Opt\n"
            "jobs 2-3:  #P B3LYP/6-31G(d) SCF=Tight NoSymm Geom=Check Guess=Read")

    def test_route_preview_shows_merge_problems(self, tmp_path):
        spec = job_spec(tmp_path, step(job_type=gk.JOB_OPT, opt_max_cycles=9, additional_keywords="Opt=MaxCycles=5"))
        text = gi.route_preview(spec, 0)
        assert text.startswith("#P B3LYP/6-31G(d) Opt=MaxCycles=9") and "Job Type panel" in text

    def test_molecule_summary(self, tmp_path):
        assert gi.molecule_summary(job_spec(tmp_path, step()), 0) == \
            "3 atoms, charge 0, 10 electrons, multiplicity ? (enter it: 10 electrons allow 1, 3, 5, ...)"
        assert gi.molecule_summary(job_spec(tmp_path, step(charge=1, multiplicity=1)), 0).endswith(
            "multiplicity 1   <- impossible: " + gi.multiplicity_problem(9, 1))
        assert "checkpoint" in gi.molecule_summary(job_spec(tmp_path, step(geometry_source=gk.GEOM_ALLCHECK)), 0)
        assert "geometry" in gi.molecule_summary(job_spec(tmp_path, geometry=""), 0)

    def test_fukui_summary_lists_the_three_states(self, tmp_path):
        text = gi.molecule_summary(fukui_job_spec(tmp_path, fukui_plus_multiplicity=None), 0)
        assert "N+1" in text and "charge -1" in text and "51 electrons" in text and "phenol_Np1.chk" in text
        assert "multiplicity ? (enter it: 51 electrons allow 2, 4, 6, ...)" in text
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_gaussianInput.py -q -k "Fukui or Preview"`
Expected: FAIL (`AttributeError: ... 'fukui_defaults'`).

- [ ] **Step 3: Implement** (append to `gaussianInput.py`, and replace the `_prepare_fukui` stub in `GaussianInputBuilder` with the method below)

```python
# ----------------------------------------------------------------------
# The Fukui Job Type (delegated to fukuiInput, unchanged)
# ----------------------------------------------------------------------
# Fields the Fukui linked job writes itself; they must hold fukui_defaults' values.
FUKUI_FIXED = ("opt_target", "force_constants", "state", "print_all_coefficients", "spin", "geometry_source",
               "modredundant", "title", "write_chk", "chk", "old_chk", "extra_link0", "connectivity", "guess",
               "guess_mix", "guess_always", "guess_only", "guess_save", "additional_input")


def fukui_defaults(step: GaussianStepSpec) -> GaussianStepSpec:
    """The step as the Fukui Job Type uses it: the former builder's SCF=Tight NoSymm, unused controls reset."""
    force_constants = step.force_constants if step.force_constants in gk.FUKUI_FORCE_CONSTANTS else gk.FC_DEFAULT
    return dataclasses.replace(
        step, job_type=gk.JOB_FUKUI, opt_target=gk.OPT_MINIMUM, force_constants=force_constants,
        state=gk.STATE_GROUND, print_all_coefficients=False, spin="Default", geometry_source=gk.GEOM_COORDINATES,
        modredundant="", title="", write_chk=True, chk="", old_chk="", extra_link0="", connectivity=False,
        guess="Default", guess_mix=False, guess_always=False, guess_only=False, guess_save=False,
        additional_input="", scf_tight=True, nosymm=True)


def fukui_opt_options(step: GaussianStepSpec) -> List[str]:
    """Options of the N-state optimization (towards a minimum; no ReadFC, no ModRedundant)."""
    options = [o for o in (gk.FORCE_CONSTANTS.get(step.force_constants),) if o]
    if step.opt_tight:
        options.append("Tight")
    if step.opt_max_cycles:
        options.append(f"MaxCycles={step.opt_max_cycles}")
    return options


def fukui_keywords(step: GaussianStepSpec) -> Tuple[List[RouteKeyword], List[str]]:
    """Keywords shared by the three Fukui jobs (panels + Additional Keywords, merged as in spec 5.1)."""
    return merge_keywords(method_keywords(step) + general_keywords(step), step.additional_keywords)


def fukui_spec(spec: GaussianJobSpec) -> fi.FukuiJobSpec:
    """The panels of the single Fukui step as the FukuiJobSpec that fukuiInput validates and writes."""
    step = spec.steps[0]
    shared = [format_keyword(k) for k in fukui_keywords(step)[0]]
    return fi.FukuiJobSpec(
        stem=(spec.stem or "").strip(), output_dir=(spec.output_dir or "").strip(), geometry=spec.geometry,
        charge=step.charge, multiplicity=step.multiplicity,
        plus_multiplicity=step.fukui_plus_multiplicity, minus_multiplicity=step.fukui_minus_multiplicity,
        memory=step.memory, nproc=step.nproc, print_level=step.print_level,
        method=method_keyword(step, with_prefix=False), basis=basis_keyword(step),
        keywords=" ".join(shared),
        opt_keyword=format_keyword(RouteKeyword("Opt", fukui_opt_options(step))),
        p=step.fukui_p, q=step.fukui_q, guess_read=step.fukui_guess_read)


# ----------------------------------------------------------------------
# Live texts of the window (never raise)
# ----------------------------------------------------------------------
def route_preview(spec: GaussianJobSpec, index: int) -> str:
    """The route line(s) of step ``index`` (0-based), with merge problems on extra lines."""
    step = spec.steps[index]
    try:
        if step.job_type == gk.JOB_FUKUI:
            f = fukui_spec(spec)
            args = (f.print_level, f.method, f.basis, f.keywords)
            ion_extra = ("Geom=Check", "Guess=Read") if f.guess_read else ("Geom=Check",)
            return ("job 1:     " + " ".join(fi.build_route(*args, (f.opt_keyword,))) + "\n"
                    "jobs 2-3:  " + " ".join(fi.build_route(*args, ion_extra)))
        keywords, errors = route_keywords(step)
        route = fi.build_route(step.print_level, method_keyword(step), basis_keyword(step),
                               " ".join(format_keyword(k) for k in keywords))
        return " ".join(route) + "".join(f"\n{e}" for e in errors)
    except (KeyError, ValueError) as exc:
        return f"Cannot build the route yet: {exc}"


def _spin_text(electrons: int, multiplicity: Optional[int]) -> str:
    if multiplicity is None:
        allowed = ", ".join(str(m) for m in allowed_multiplicities(electrons, 3))
        return f"multiplicity ? (enter it: {electrons} electrons allow {allowed}, ...)"
    problem = multiplicity_problem(electrons, multiplicity) if electrons >= 1 else "no electrons left"
    return f"multiplicity {multiplicity}" + (f"   <- impossible: {problem}" if problem else "")


def molecule_summary(spec: GaussianJobSpec, index: int) -> str:
    """Charge, electron count and multiplicity of step ``index`` (the Fukui job: of its three states)."""
    step = spec.steps[index]
    if step.geometry_source == gk.GEOM_ALLCHECK:
        return "Geometry, charge and multiplicity come from the checkpoint."
    if not (spec.geometry or "").strip():
        return "Enter or load a geometry to see the electron count."
    try:
        atoms = fi.parse_geometry(spec.geometry, check_distances=False)
    except ValueError as exc:
        return f"Cannot count electrons yet: {exc}"
    z_total = nuclear_charge(a[0] for a in atoms)
    if step.job_type != gk.JOB_FUKUI:
        electrons = z_total - step.charge
        return (f"{len(atoms)} atoms, charge {step.charge:+d}, {electrons} electrons, "
                f"{_spin_text(electrons, step.multiplicity)}").replace("charge +0", "charge 0")
    stem = (spec.stem or "").strip() or fi.DEFAULT_STEM
    lines = []
    multiplicities = (step.multiplicity, step.fukui_plus_multiplicity, step.fukui_minus_multiplicity)
    for (label, tag, _, k, sign), mult in zip(fi.state_tags(step.fukui_p, step.fukui_q), multiplicities):
        charge = step.charge + sign * k
        electrons = z_total - charge
        lines.append(f"{label:<5} charge {charge:+d}, {electrons} electrons -> {stem}_{tag}.chk, "
                     f"{_spin_text(electrons, mult)}")
    return "\n".join(lines).replace("charge +0", "charge 0")
```

The method that replaces the stub inside `GaussianInputBuilder`:

```python
    def _prepare_fukui(self, spec: GaussianJobSpec) -> GaussianJob:
        if len(spec.steps) != 1:
            raise ValueError("- The Fukui functions job type is already three linked jobs (N, N+p, N-q): "
                             "it must be the only step.")
        step = spec.steps[0]
        reference = fukui_defaults(step)
        errors = [f"The Fukui functions job type sets {FIELD_LABELS[name]} itself; choose the job type again "
                  "to reset it." for name in FUKUI_FIXED if getattr(step, name) != getattr(reference, name)]
        keywords, merge_errors = fukui_keywords(step)
        errors += merge_errors
        if any(k.key == "OPT" for k in keywords):
            errors.append("The Fukui functions job type: set the N-state optimization with the Optimization "
                          "options of the Job Type panel, not in Additional Keywords.")
        if errors:
            raise ValueError("\n".join(f"- {e}" for e in errors))
        job = self._fukui.prepare(fukui_spec(spec))
        return GaussianJob(spec, job.atoms, [], job.input_path, job.script_path, [s.chk for s in job.states],
                           list(job.warnings), fukui=job)
```

`gi.multiplicity_problem` in the test is the name `gaussianInput` imported from `systemProperties`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_gaussianInput.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add gaussianInput.py tests/test_gaussianInput.py
git commit -m "feat(gaussian): Fukui functions job type delegated to fukuiInput"
```

---

### Task 7: Reopening a Gaussian input (`read_input_text` / `read_input_file`)

**Files:**
- Modify: `gaussianInput.py` (append)
- Test: `tests/test_gaussianInputReopen.py` (create)

**Interfaces:**
- Consumes: everything above.
- Produces: `read_input_text(text: str, stem: str = fi.DEFAULT_STEM, output_dir: str = "") -> Tuple[GaussianJobSpec, List[str]]` (spec, notes); `read_input_file(path) -> Tuple[GaussianJobSpec, List[str]]` (stem from the file name, output folder = the file's folder); `GAUSSIAN_INPUT_SUFFIXES = (".gjf", ".com", ".gau")`. Both raise `ValueError` naming the step (and the file).

Rules (spec section 6): split on `--Link1--`; per step `%` lines → Link 0, `#` lines → route, then title and charge/multiplicity unless `Geom=AllCheck`, Cartesian geometry unless `Geom=Check`/`AllCheck`, the connectivity section if `Geom=Connectivity` (dropped: rebuilt on save, with a note), the ModRedundant section if `Opt=ModRedundant`, the rest → Add. Inp. A `%Chk` equal to `<stem>.chk` and a title equal to the default title come back blank (the builder's defaults), so a file the builder wrote reads back to the spec it came from. Titles with characters Gaussian forbids are cleaned, with a note. Absolute checkpoint paths keep the file name, with a note.

- [ ] **Step 1: Write the failing tests**

`tests/test_gaussianInputReopen.py`:

```python
"""Open input...: Gaussian inputs read back into the panels (gaussianInput, no Toga)."""
from dataclasses import replace
from pathlib import Path

import pytest

import fukuiInput as fi
import gaussianInput as gi
import gaussianKeywords as gk

DATA = Path(__file__).resolve().parent / "data"
FIXTURES = sorted(p for p in DATA.rglob("*") if p.suffix.lower() in gi.GAUSSIAN_INPUT_SUFFIXES)
WATER = ("O 0.000000 0.000000 0.117300\n"
         "H 0.000000 0.757200 -0.469200\n"
         "H 0.000000 -0.757200 -0.469200")


def step(**changes):
    return replace(gi.GaussianStepSpec(), **changes)


def spec_of(tmp_path, *steps, geometry=WATER):
    return gi.GaussianJobSpec(stem="water", output_dir=str(tmp_path), geometry=geometry, steps=list(steps))


ROUND_TRIPS = {
    "energy": [step(multiplicity=1)],
    "opt freq": [step(job_type=gk.JOB_OPT_FREQ, multiplicity=1, freq_raman=True, freq_anharmonic=True)],
    "ts": [step(job_type=gk.JOB_OPT, multiplicity=1, opt_target=gk.OPT_TS, force_constants=gk.FC_CALCFC,
                opt_tight=True, opt_max_cycles=200, additional_keywords="Opt=NoEigenTest 6D")],
    "irc": [step(job_type=gk.JOB_IRC, multiplicity=1, irc_direction="Reverse only", irc_max_points=30,
                 irc_step_size=8)],
    "irc default fc": [step(job_type=gk.JOB_IRC, multiplicity=1, irc_force_constants=gk.FC_DEFAULT)],
    "scan": [step(job_type=gk.JOB_SCAN, multiplicity=1, modredundant="B 1 2 S 5 0.05\nA 2 1 3 F",
                  connectivity=True)],
    "stability": [step(job_type=gk.JOB_STABILITY, multiplicity=1, stable_opt=True)],
    "nmr": [step(job_type=gk.JOB_NMR, multiplicity=1, nmr_method="CSGT", nmr_spinspin=True)],
    "td": [step(multiplicity=1, state=gk.STATE_TD, n_states=10, root=2, state_spin="Singlets",
                print_all_coefficients=True, dispersion="GD3BJ")],
    "cis": [step(multiplicity=1, state=gk.STATE_CIS, n_states=4, family=gk.FAMILY_HF, method_name=gk.NO_CHOICE,
                 basis="cc-pVDZ")],
    "zindo": [step(multiplicity=1, state=gk.STATE_ZINDO, n_states=8, family=gk.FAMILY_HF,
                   method_name=gk.NO_CHOICE)],
    "semi": [step(multiplicity=1, family=gk.FAMILY_SEMI, method_name="PM6")],
    "compound": [step(multiplicity=1, family=gk.FAMILY_COMPOUND, method_name="G4")],
    "ccsd(t)": [step(multiplicity=1, family=gk.FAMILY_CCSDT, method_name=gk.NO_CHOICE, basis="aug-cc-pVTZ")],
    "open shell": [step(charge=1, multiplicity=2, spin="Unrestricted (U)", method_name="M06-2X",
                        basis="6-311G", basis_diffuse="++", basis_polarization="(2df,2pd)")],
    "other names": [step(multiplicity=1, method_name=gk.OTHER, method_other="revPBE", basis=gk.OTHER,
                         basis_other="def2TZVPD", solvation="CPCM", solvent=gk.OTHER, solvent_other="Heptane")],
    "general": [step(multiplicity=1, print_level="#N", scf_tight=True, scf_qc=True, nosymm=True, polar=True,
                     pop_level="Full", pop_mk=True, pop_nbo=True, solvation="SMD", solvent="Acetonitrile",
                     title="my water", memory="8GB", nproc=16, chk="w1.chk", old_chk="start.chk",
                     extra_link0="%KJob L502")],
    "gen basis": [step(multiplicity=1, basis="Gen", additional_input="O H 0\n6-31G(d)\n****")],
    "no checkpoint": [step(multiplicity=1, write_chk=False)],
    "linked": [step(job_type=gk.JOB_OPT, multiplicity=1),
               replace(gi.linked_step(step(job_type=gk.JOB_OPT, multiplicity=1)), state=gk.STATE_TD, n_states=6,
                       guess_mix=True),
               step(geometry_source=gk.GEOM_CHECK, charge=-1, multiplicity=2, guess="Read checkpoint"),
               step(multiplicity=3, title="triplet")],
}


@pytest.mark.parametrize("name", sorted(ROUND_TRIPS))
def test_a_file_the_builder_wrote_reads_back_to_its_spec(tmp_path, name):
    builder = gi.GaussianInputBuilder()
    spec = spec_of(tmp_path, *ROUND_TRIPS[name])
    text = builder.input_text(builder.prepare(spec))
    back, _ = gi.read_input_text(text, spec.stem, spec.output_dir)
    assert back == spec


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_every_project_input_reopens_prepares_and_keeps_its_keywords(tmp_path, path):
    spec, _ = gi.read_input_file(path)
    original = path.read_text(errors="replace").splitlines()
    route = " ".join(l.strip() for l in original if l.strip().startswith("#")).split(None, 1)[1]
    names = {gi.keyword_text_name(t).upper() for t in fi.route_tokens(route) if "/" not in t}
    emitted = {gi.keyword_text_name(t).upper() for t in fi.route_tokens(" ".join(gi.route_lines(spec.steps[0])))}
    assert names <= emitted | {gi.keyword_text_name(gi.method_keyword(spec.steps[0])).upper()}
    builder = gi.GaussianInputBuilder()
    job = builder.prepare(replace(spec, output_dir=str(tmp_path)))
    again, _ = gi.read_input_text(builder.input_text(job), spec.stem, str(tmp_path))
    assert again == replace(spec, output_dir=str(tmp_path))


class TestReading:
    def test_star_basis_sets_and_lower_case(self):
        spec, _ = gi.read_input_text("#p pbe1pbe/6-31g** opt\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n")
        s = spec.steps[0]
        assert (s.method_name, s.basis, s.basis_polarization, s.job_type, s.print_level) == \
            ("PBE0 (PBE1PBE)", "6-31G", "(d,p)", gk.JOB_OPT, "#P")

    def test_a_route_without_a_method_is_hf_sto3g_with_a_note(self):
        spec, notes = gi.read_input_text("# opt\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n")
        s = spec.steps[0]
        assert (s.family, s.basis) == (gk.FAMILY_HF, "STO-3G") and any("HF/STO-3G" in n for n in notes)

    def test_absolute_checkpoint_path_keeps_the_file_name(self):
        spec, notes = gi.read_input_text("%chk=C:\\jobs\\w.chk\n# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n",
                                         stem="other")
        assert spec.steps[0].chk == "w.chk" and any("C:\\jobs\\w.chk" in n for n in notes)

    def test_forbidden_title_characters_are_replaced(self):
        spec, notes = gi.read_input_text("# hf/sto-3g\n\nB3LYP/6-31G* opted\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n")
        assert spec.steps[0].title == "B3LYP/6 31G* opted" and any("title" in n for n in notes)

    def test_connectivity_section_is_consumed_with_a_note(self):
        text = "# hf/sto-3g geom=connectivity\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n 1 2 1.0\n 2\n\nextra\n\n"
        spec, notes = gi.read_input_text(text)
        assert spec.steps[0].connectivity and spec.steps[0].additional_input == "extra"
        assert any("connectivity" in n for n in notes)

    def test_fukui_file_reopens_as_three_steps_with_a_note(self, tmp_path):
        s = gi.fukui_defaults(step(multiplicity=1, fukui_plus_multiplicity=2, fukui_minus_multiplicity=2))
        builder = gi.GaussianInputBuilder()
        spec = spec_of(tmp_path, s)
        spec.stem = "water"
        text = builder.input_text(builder.prepare(spec))
        back, notes = gi.read_input_text(text, "water_fukui")
        assert len(back.steps) == 3 and any("Fukui" in n for n in notes)
        assert back.steps[1].geometry_source == gk.GEOM_CHECK and back.steps[1].old_chk == "water_N.chk"

    @pytest.mark.parametrize("text,fragment", [
        ("", "empty"),
        ("%mem=1GB\n\nt\n", "no route section"),
        ("# hf/sto-3g\n\nt\n\n0 1 0 1 0 1\nH 0 0 0\nH 0 0 0.74\n\n", "per-fragment charges"),
        ("# hf/sto-3g\n\nt\n\n0 1\nO\nH 1 0.96\n\n", "Only Cartesian coordinates"),
        ("# hf/sto-3g\n\nt\n\n", "ends before the molecule specification"),
        ("# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n--Link1--\n# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.80\n\n",
         "Step 2 has its own geometry"),
        ("# hf/sto-3g geom=check\n\nt\n\n", "Geom=Check needs a title"),
        ("mem=1GB\n# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\n\n", "is not a Link 0"),
    ])
    def test_refusals(self, text, fragment):
        with pytest.raises(ValueError, match=fragment.replace("(", "\\(").replace(")", "\\)")):
            gi.read_input_text(text)

    def test_read_input_file_names_the_file(self, tmp_path):
        path = tmp_path / "bad.gjf"
        path.write_text("nothing here\n")
        with pytest.raises(ValueError, match="bad.gjf"):
            gi.read_input_file(path)
        with pytest.raises(ValueError, match="not a Gaussian input"):
            gi.read_input_file(tmp_path / "x.xyz")

    def test_read_input_file_takes_stem_and_folder_from_the_path(self, tmp_path):
        path = tmp_path / "my mol.com"
        path.write_text("# hf/sto-3g\n\nt\n\n0 1\nH 0 0 0\nH 0 0 0.74\n\n")
        spec, _ = gi.read_input_file(path)
        assert (spec.stem, spec.output_dir) == ("my_mol", str(tmp_path.resolve()))
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_gaussianInputReopen.py -q`
Expected: collection error `AttributeError: module 'gaussianInput' has no attribute 'GAUSSIAN_INPUT_SUFFIXES'`.

- [ ] **Step 3: Implement** (append to `gaussianInput.py`)

```python
# ----------------------------------------------------------------------
# Open input...: a Gaussian input back into the panels
# ----------------------------------------------------------------------
GAUSSIAN_INPUT_SUFFIXES = (".gjf", ".com", ".gau")
_POPLE_RE = re.compile(r"^(6-31|6-311)(\+{0,2})G(\*{1,2}|\([^)]*\))?$", re.IGNORECASE)
_STAR_POLARIZATION = {"*": "(d)", "**": "(d,p)"}
_INTEGER_RE = re.compile(r"[+-]?\d+")


def _int_option(option: str) -> Optional[int]:
    value = option.split("=", 1)[1].strip() if "=" in option else ""
    return int(value) if value.isdigit() and int(value) > 0 else None


def _bare_method(upper: str) -> Optional[dict]:
    """Panel fields of a method keyword without prefix (upper case), or None if no panel knows it."""
    families = {v.upper(): k for k, v in gk.FAMILY_KEYWORDS.items()}
    families["MP4"] = gk.FAMILY_MP4
    if upper in families:
        return {"family": families[upper], "method_name": gk.NO_CHOICE}
    for family, table in gk.FAMILY_METHODS.items():
        back = gk.reverse(table)
        if upper in back:
            return {"family": family, "method_name": back[upper]}
    name = keyword_text_name(upper)
    states = {"CIS": gk.STATE_CIS, "ZINDO": gk.STATE_ZINDO}
    if name in states:
        return {"state": states[name], "family": gk.FAMILY_HF, "method_name": gk.NO_CHOICE,
                "_options": parse_keyword(upper).options}
    return None


def _method_fields(text: str) -> Optional[dict]:
    """'UB3LYP' -> {'spin': 'Unrestricted (U)', 'family': 'DFT', 'method_name': 'B3LYP'}; None if unknown."""
    upper = text.upper()
    fields = _bare_method(upper)
    if fields is not None:
        return {"spin": "Default", **fields}
    for prefix in ("RO", "R", "U"):
        if upper.startswith(prefix):
            fields = _bare_method(upper[len(prefix):])
            if fields is not None:
                return {"spin": gk.reverse(gk.SPIN_PREFIXES)[prefix], **fields}
    return None


class _RouteReader:
    """Puts route tokens back into a GaussianStepSpec; what no panel models is kept as extra keywords."""

    def __init__(self, step: GaussianStepSpec, notes: List[str], where: str):
        self.step, self.notes, self.where = step, notes, where
        self.jobs: List[str] = []
        self.modredundant = False
        self.extra: List[str] = []
        self.handlers = {"OPT": self._opt, "FREQ": self._freq, "IRC": self._irc, "STABLE": self._stable,
                         "NMR": self._nmr, "TD": self._td, "TDA": self._td,
                         "EMPIRICALDISPERSION": self._dispersion, "IOP": self._iop, "SCF": self._scf,
                         "NOSYMM": self._switch("nosymm"), "POLAR": self._switch("polar"), "GEOM": self._geom,
                         "GUESS": self._guess, "POP": self._pop, "SCRF": self._scrf}

    def read(self, tokens: Sequence[str]) -> None:
        model_seen = False
        for token in tokens:
            if not model_seen and self._model(token):
                model_seen = True
                continue
            keyword = parse_keyword(token)
            handler = self.handlers.get(keyword.key)
            leftover = handler(keyword) if handler else None
            if leftover is None:
                self.extra.append(token)
            elif leftover:
                self.extra.append(format_keyword(RouteKeyword(keyword.name, leftover)))
        if not model_seen:
            self.step.family, self.step.method_name, self.step.basis = gk.FAMILY_HF, gk.NO_CHOICE, "STO-3G"
            self.notes.append(f"{self.where}: the route names no method; Gaussian's default HF/STO-3G is used.")

    def job_type(self, scan_lines: bool) -> str:
        jobs = self.jobs
        kinds = [j for j in ("OPT", "IRC", "STABLE", "NMR") if j in jobs]
        if len(kinds) > 1:
            self.notes.append(f"{self.where}: the route combines {', '.join(kinds)}; only {kinds[0]} is kept as "
                              "the job type.")
        if "OPT" in jobs:
            return gk.JOB_OPT_FREQ if "FREQ" in jobs else (gk.JOB_SCAN if scan_lines else gk.JOB_OPT)
        if "FREQ" in jobs:
            return gk.JOB_FREQ
        return {"IRC": gk.JOB_IRC, "STABLE": gk.JOB_STABILITY, "NMR": gk.JOB_NMR}.get(kinds[0] if kinds else "",
                                                                                    gk.JOB_ENERGY)

    # --- model chemistry ---
    def _model(self, token: str) -> bool:
        method_text, slash, basis_text = token.partition("/")
        if slash and method_text.count("(") != method_text.count(")"):
            return False                                     # IOp(9/40=4)
        fields = _method_fields(method_text)
        if fields is None:
            if not slash:
                return False
            fields = {"spin": "Default", "family": gk.FAMILY_DFT, "method_name": gk.OTHER,
                      "method_other": method_text}
        options = fields.pop("_options", [])
        for name, value in fields.items():
            setattr(self.step, name, value)
        dropped = self._state_options(options) if self.step.state in gk.MODEL_STATES else []
        if dropped:
            self.notes.append(f"{self.where}: option(s) {', '.join(dropped)} of {self.step.state} are not kept.")
        if slash:
            self._basis(basis_text)
        elif not no_basis(self.step):
            self.step.basis = "STO-3G"
            self.notes.append(f"{self.where}: the route names no basis set; Gaussian's default STO-3G is used.")
        return True

    def _basis(self, text: str) -> None:
        step = self.step
        match = _POPLE_RE.match(text)
        if match:
            polarization = (match.group(3) or "")
            polarization = _STAR_POLARIZATION.get(polarization, polarization).lower()
            back = {v.lower(): k for k, v in gk.POLARIZATION.items()}
            if polarization in back:
                step.basis = match.group(1) + "G"
                step.basis_diffuse = {v: k for k, v in gk.DIFFUSE.items()}[match.group(2)]
                step.basis_polarization = back[polarization]
                return
        known = {b.upper(): b for b in gk.POPLE_BASES + gk.OTHER_BASES}
        if text.upper() in known:
            step.basis = known[text.upper()]
        else:
            step.basis, step.basis_other = gk.OTHER, text

    def _state_options(self, options: Sequence[str]) -> List[str]:
        step, rest = self.step, []
        spins = gk.reverse(gk.STATE_SPINS)
        step.n_states, step.root, step.state_spin = gk.DEFAULT_N_STATES, 1, "Default"
        for option in options:
            key = option_key(option)
            if key == "NSTATES" and _int_option(option):
                step.n_states = _int_option(option)
            elif key == "ROOT" and _int_option(option):
                step.root = _int_option(option)
            elif option.upper() in spins:
                step.state_spin = spins[option.upper()]
            else:
                rest.append(option)
        return rest

    # --- keywords (each returns the options it does not model; None = keep the token as it is) ---
    def _opt(self, keyword: RouteKeyword) -> List[str]:
        self.jobs.append("OPT")
        rest, fc = [], gk.reverse(gk.FORCE_CONSTANTS)
        for option in keyword.options:
            key = option_key(option)
            if key == "TS":
                self.step.opt_target = gk.OPT_TS
            elif key in fc:
                self.step.force_constants = fc[key]
            elif key == "TIGHT":
                self.step.opt_tight = True
            elif key == "MAXCYCLES" and _int_option(option):
                self.step.opt_max_cycles = _int_option(option)
            elif key == "MODREDUNDANT":
                self.modredundant = True
            else:
                rest.append(option)
        return rest

    def _flags(self, keyword: RouteKeyword, table: Dict[str, str]) -> List[str]:
        back = {option.upper(): name for name, option in table.items()}
        rest = []
        for option in keyword.options:
            if option.upper() in back:
                setattr(self.step, back[option.upper()], True)
            else:
                rest.append(option)
        return rest

    def _freq(self, keyword: RouteKeyword) -> List[str]:
        self.jobs.append("FREQ")
        return self._flags(keyword, gk.FREQ_FLAG_FIELDS)

    def _irc(self, keyword: RouteKeyword) -> List[str]:
        self.jobs.append("IRC")
        step, rest = self.step, []
        step.irc_force_constants = gk.FC_DEFAULT
        directions, fc = gk.reverse(gk.IRC_DIRECTIONS), gk.reverse(gk.FORCE_CONSTANTS)
        for option in keyword.options:
            key = option_key(option)
            if key in directions:
                step.irc_direction = directions[key]
            elif key in fc:
                step.irc_force_constants = fc[key]
            elif key == "MAXPOINTS" and _int_option(option):
                step.irc_max_points = _int_option(option)
            elif key == "STEPSIZE" and _int_option(option):
                step.irc_step_size = _int_option(option)
            else:
                rest.append(option)
        return rest

    def _stable(self, keyword: RouteKeyword) -> List[str]:
        self.jobs.append("STABLE")
        rest = [o for o in keyword.options if o.upper() != "OPT"]
        self.step.stable_opt = len(rest) < len(keyword.options)
        return rest

    def _nmr(self, keyword: RouteKeyword) -> List[str]:
        self.jobs.append("NMR")
        methods, rest = gk.reverse(gk.NMR_METHODS), []
        for option in keyword.options:
            key = option.upper()
            if key == "GIAO":
                self.step.nmr_method = "GIAO"
            elif key in methods:
                self.step.nmr_method = methods[key]
            elif key == "SPINSPIN":
                self.step.nmr_spinspin = True
            else:
                rest.append(option)
        return rest

    def _td(self, keyword: RouteKeyword) -> List[str]:
        self.step.state = gk.STATE_TD if keyword.key == "TD" else gk.STATE_TDA
        return self._state_options(keyword.options)

    def _dispersion(self, keyword: RouteKeyword) -> Optional[List[str]]:
        back = gk.reverse(gk.DISPERSION)
        if len(keyword.options) == 1 and keyword.options[0].upper() in back:
            self.step.dispersion = back[keyword.options[0].upper()]
            return []
        return None

    def _iop(self, keyword: RouteKeyword) -> List[str]:
        rest = [o for o in keyword.options if o.replace(" ", "") != gk.COEFFICIENTS_IOP]
        self.step.print_all_coefficients = len(rest) < len(keyword.options)
        return rest

    def _scf(self, keyword: RouteKeyword) -> List[str]:
        rest = []
        for option in keyword.options:
            if option.upper() == "TIGHT":
                self.step.scf_tight = True
            elif option.upper() == "QC":
                self.step.scf_qc = True
            else:
                rest.append(option)
        return rest

    def _switch(self, name: str):
        def handler(keyword: RouteKeyword) -> Optional[List[str]]:
            if keyword.options:
                return None
            setattr(self.step, name, True)
            return []
        return handler

    def _geom(self, keyword: RouteKeyword) -> List[str]:
        sources, rest = gk.reverse(gk.GEOMETRY_SOURCES), []
        for option in keyword.options:
            key = option.upper()
            if key in sources:
                self.step.geometry_source = sources[key]
            elif key in ("CONNECTIVITY", "CONNECT"):
                self.step.connectivity = True
            else:
                rest.append(option)
        return rest

    def _guess(self, keyword: RouteKeyword) -> List[str]:
        methods = gk.reverse(gk.GUESS_METHODS)
        flags = {option.upper(): name for name, option in gk.GUESS_FLAG_FIELDS.items()}
        rest = []
        for option in keyword.options:
            key = option.upper()
            if key in methods:
                self.step.guess = methods[key]
            elif key in flags:
                setattr(self.step, flags[key], True)
            else:
                rest.append(option)
        return rest

    def _pop(self, keyword: RouteKeyword) -> List[str]:
        levels = gk.reverse(gk.POP_LEVELS)
        flags = {option.upper(): name for name, option in gk.POP_FLAG_FIELDS.items()}
        rest = []
        for option in keyword.options:
            key = option.upper()
            if key in levels:
                self.step.pop_level = levels[key]
            elif key in flags:
                setattr(self.step, flags[key], True)
            else:
                rest.append(option)
        return rest

    def _scrf(self, keyword: RouteKeyword) -> List[str]:
        models = gk.reverse(gk.SOLVATION_MODELS)
        solvents = {s.upper(): s for s in gk.SOLVENTS if s != gk.OTHER}
        step, rest = self.step, []
        step.solvation, step.solvent = "IEFPCM", "Water"   # Gaussian's defaults for a bare SCRF
        for option in keyword.options:
            key = gk.SOLVATION_ALIASES.get(option.upper(), option.upper())
            if key in models:
                step.solvation = models[key]
            elif option_key(option) == "SOLVENT" and "=" in option:
                name = option.split("=", 1)[1].strip()
                if name.upper() in solvents:
                    step.solvent = solvents[name.upper()]
                else:
                    step.solvent, step.solvent_other = gk.OTHER, name
            else:
                rest.append(option)
        return rest


def _split_blocks(text: str) -> List[List[str]]:
    blocks, current = [], []
    for raw in (text or "").splitlines():
        line = raw.split("!", 1)[0].rstrip()          # '!' starts a comment in Gaussian input
        if line.strip().lower() == "--link1--":
            blocks.append(current)
            current = []
        else:
            current.append(line)
    blocks.append(current)
    return [b for b in blocks if any(line.strip() for line in b)]


def _sections(lines: Sequence[str]) -> List[List[str]]:
    sections, current = [], []
    for line in lines:
        if line.strip():
            current.append(line)
        elif current:
            sections.append(current)
            current = []
    if current:
        sections.append(current)
    return sections


def _reopened_name(value: str, label: str, where: str, notes: List[str]) -> str:
    name = PureWindowsPath(value.strip()).name
    if name != value.strip():
        notes.append(f"{where}: {label}={value.strip()} becomes {name} (the job runs from the input's folder).")
    return name


def _read_charge(step: GaussianStepSpec, line: str, where: str) -> None:
    tokens = line.replace(",", " ").split()
    if len(tokens) != 2 or not all(_INTEGER_RE.fullmatch(t) for t in tokens):
        raise ValueError(f"{where}: '{line.strip()}' is not one charge and one multiplicity "
                         "(per-fragment charges are not supported).")
    step.charge, step.multiplicity = int(tokens[0]), int(tokens[1])


def _read_block(lines: Sequence[str], index: int, stem: str, notes: List[str]
                ) -> Tuple[GaussianStepSpec, Optional[str]]:
    """One linked step -> (its panels, its geometry text or None)."""
    where = f"Step {index}"
    sections = _sections(lines)
    head = sections[0]
    route_start = next((k for k, line in enumerate(head) if line.lstrip().startswith("#")), None)
    if route_start is None:
        raise ValueError(f"{where}: no route section (a line starting with '#') was found.")
    step = GaussianStepSpec(write_chk=False)

    extra_link0, seen = [], set()
    for line in head[:route_start]:
        text = line.strip()
        if not text.startswith("%"):
            raise ValueError(f"{where}: '{text}' before the route is not a Link 0 (%) line.")
        key, _, value = text[1:].partition("=")
        key, value = key.strip().lower(), value.strip()
        seen.add(key)
        if key == "mem":
            step.memory = value
        elif key in ("nprocshared", "nproc") and value.isdigit():
            step.nproc = int(value)
        elif key == "chk":
            step.write_chk = True
            name = _reopened_name(value, "%Chk", where, notes)
            step.chk = "" if name == f"{stem}.chk" else name
        elif key == "oldchk":
            step.old_chk = _reopened_name(value, "%OldChk", where, notes)
        else:
            extra_link0.append(text)
    step.extra_link0 = "\n".join(extra_link0)
    if "mem" not in seen:
        notes.append(f"{where}: no %Mem in the file; {fi.DEFAULT_MEMORY} is used.")

    route_text = " ".join(line.strip() for line in head[route_start:])
    match = re.match(r"#([PNTpnt]?)(?=\s|$)", route_text)
    step.print_level = "#" + (match.group(1).upper() if match else "")
    reader = _RouteReader(step, notes, where)
    reader.read(fi.route_tokens(route_text[match.end():] if match else route_text[1:]))

    rest = sections[1:]
    geometry = None
    title_lines: List[str] = []
    if step.geometry_source == gk.GEOM_CHECK:
        if len(rest) < 2:
            raise ValueError(f"{where}: Geom=Check needs a title and a charge/multiplicity line.")
        title_lines, charge_lines, rest = rest[0], rest[1], rest[2:]
        _read_charge(step, charge_lines[0], where)
    elif step.geometry_source == gk.GEOM_COORDINATES:
        if len(rest) < 2:
            raise ValueError(f"{where}: the file ends before the molecule specification.")
        title_lines, molecule, rest = rest[0], rest[1], rest[2:]
        _read_charge(step, molecule[0], where)
        try:
            atoms = fi.parse_geometry("\n".join(molecule[1:]), check_distances=False)
        except ValueError as exc:
            raise ValueError(f"{where}: {exc}") from exc
        geometry = "\n".join(fi.format_geometry(atoms))
    if title_lines:
        title = " ".join(line.strip() for line in title_lines)
        clean = " ".join("".join(" " if ch in fi.TITLE_FORBIDDEN else ch for ch in title).split())
        if clean != title:
            notes.append(f"{where}: the title had characters Gaussian does not allow; they became blanks.")
        step.title = "" if clean == fi._title_text(stem) else clean

    if step.connectivity:
        if rest and rest[0][0].split()[0].isdigit():
            rest = rest[1:]
            notes.append(f"{where}: the connectivity section is rebuilt from the distances when saving.")
        else:
            notes.append(f"{where}: Geom=Connectivity but no connectivity section follows.")
    scan_lines = False
    if reader.modredundant:
        if rest:
            step.modredundant = "\n".join(line.strip() for line in rest[0])
            rest = rest[1:]
            try:
                scan_lines = any(p[2] == "S" for p in parse_modredundant(step.modredundant))
            except ValueError:
                scan_lines = False      # prepare reports the line
        else:
            notes.append(f"{where}: Opt=ModRedundant but no coordinate lines follow.")
    step.job_type = reader.job_type(scan_lines)
    step.additional_input = "\n\n".join("\n".join(section) for section in rest)
    step.additional_keywords = " ".join(reader.extra)
    return step, geometry


def read_input_text(text: str, stem: str = fi.DEFAULT_STEM, output_dir: str = ""
                    ) -> Tuple[GaussianJobSpec, List[str]]:
    """A Gaussian input -> (job spec, notes). ValueError names the step."""
    blocks = _split_blocks(text)
    if not blocks:
        raise ValueError("The input is empty.")
    notes: List[str] = []
    steps: List[GaussianStepSpec] = []
    geometry = ""
    for index, lines in enumerate(blocks, start=1):
        step, step_geometry = _read_block(lines, index, stem, notes)
        if step_geometry is not None:
            if not geometry:
                geometry = step_geometry
            elif step_geometry != geometry:
                raise ValueError(f"Step {index} has its own geometry, different from the first one; the builder "
                                 "writes one geometry per job.")
        steps.append(step)
    if len(steps) == 3 and " Fukui " in f" {steps[0].title} ":
        notes.append("This is a Fukui job written by the builder; it is shown as three linked steps. Choose the "
                     "Fukui functions job type to rebuild it from its panels.")
    return GaussianJobSpec(stem=stem, output_dir=output_dir, geometry=geometry, steps=steps), notes


def read_input_file(path) -> Tuple[GaussianJobSpec, List[str]]:
    path = Path(path)
    if path.suffix.lower() not in GAUSSIAN_INPUT_SUFFIXES:
        raise ValueError(f"'{path.name}' is not a Gaussian input ({', '.join(GAUSSIAN_INPUT_SUFFIXES)}).")
    try:
        text = path.read_bytes().decode("utf-8-sig", errors="replace")
    except OSError as exc:
        raise ValueError(f"'{path.name}' could not be opened: {exc.strerror or exc}") from exc
    stem = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in path.stem) or fi.DEFAULT_STEM
    try:
        return read_input_text(text, stem, str(path.resolve().parent))
    except ValueError as exc:
        raise ValueError(f"{path.name}: {exc}") from exc
```

While implementing, check these against the tests:
- In `_read_block`, the first sections entry always exists because `_split_blocks` drops empty blocks.
- The "%mem=1GB\n\nt\n" refusal case: the head is `%mem=1GB` with no `#` line, so the message is "no route section".
- The Fukui file test: the ion steps have `%OldChk=water_N.chk` and `Geom=Check`, and their titles contain " Fukui ".
- `test_a_file_the_builder_wrote_reads_back_to_its_spec["linked"]`: step 3 (`Geom=Check`) reads `water.chk`, written by steps 1–2; step 4 writes the coordinates again, and its geometry equals step 1's.

If a fixture in `test_every_project_input_reopens_prepares_and_keeps_its_keywords` fails, fix the reader, not the test: the fixtures are real Gaussian inputs.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_gaussianInputReopen.py tests/test_gaussianInput.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add gaussianInput.py tests/test_gaussianInputReopen.py
git commit -m "feat(gaussian): reopen Gaussian inputs into the panels"
```

---

### Task 8: Geometry files and WSL paths (moved out of the Fukui window)

**Files:**
- Modify: `gaussianInput.py` (append)
- Test: `tests/test_gaussianInput.py` (append; these tests replace `TestGeometryFiles` and `test_wsl_path` of `tests/test_fukuiInputBuilder.py`, which Task 11 deletes)

**Interfaces:**
- Consumes: `molecularFileFormats.read_structure_file`, `molecularFileFormats.BOHR_RADIUS_ANGSTROM`, `gaussianWavefunction.read_fchk_sections` (lazy).
- Produces: `GEOMETRY_FILE_TYPES: List[str]`; `read_geometry_file(path) -> Tuple[str, Optional[int], Optional[int]]` (Å `El x y z` text, charge or None, multiplicity or None); `wsl_path(path) -> str`.

The old `fukuiInputBuilder.read_geometry_file` used `coordinateConverter`, which imports Toga at module level, so it cannot move as it is into a Toga-free module. The new one reads through `molecularFileFormats` (Toga-free; it converts `Units=Bohr` to Å instead of refusing it). HyperChem `.hin` is no longer offered: `molecularFileFormats` deliberately does not read it (its CLAUDE.md entry explains why).

- [ ] **Step 1: Write the failing tests** (append to `tests/test_gaussianInput.py`)

```python
import math

from fukui_fixtures import N

FUKUI_INPUTS = SRC / "tests" / "data" / "fukui" / "inputs"


@pytest.mark.parametrize("windows,linux", [(r"C:\jobs\my mol", "/mnt/c/jobs/my mol"),
                                           (r"D:\x", "/mnt/d/x"), ("/home/u/jobs", "/home/u/jobs")])
def test_wsl_path(windows, linux):
    assert gi.wsl_path(windows) == linux


class TestGeometryFiles:
    def test_gaussian_input_gives_charge_and_multiplicity(self):
        text, charge, mult = gi.read_geometry_file(FUKUI_INPUTS / "phenol_Np1.gjf")
        assert (charge, mult) == (-1, 2) and len(fi.parse_geometry(text)) == 13

    def test_fchk_gives_angstrom_coordinates(self):
        text, charge, mult = gi.read_geometry_file(N)
        atoms = fi.parse_geometry(text)
        reference = fi.parse_geometry(gi.read_geometry_file(FUKUI_INPUTS / "phenol_N.gjf")[0])
        assert (charge, mult) == (0, 1) and [a[0] for a in atoms] == [a[0] for a in reference]
        assert math.dist(atoms[0][1:], atoms[1][1:]) == pytest.approx(math.dist(reference[0][1:],
                                                                               reference[1][1:]), abs=1e-5)

    def test_xyz_has_no_charge(self, tmp_path):
        path = tmp_path / "w.xyz"
        path.write_text("3\nwater\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692\n")
        text, charge, mult = gi.read_geometry_file(path)
        assert charge is None and mult is None and len(fi.parse_geometry(text)) == 3

    def test_unsupported_file(self, tmp_path):
        path = tmp_path / "x.log"
        path.write_text("")
        with pytest.raises(ValueError, match="log"):
            gi.read_geometry_file(path)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_gaussianInput.py -q -k "wsl or GeometryFiles"`
Expected: FAIL (`AttributeError: ... 'wsl_path'`).

- [ ] **Step 3: Implement** (append to `gaussianInput.py`)

```python
# ----------------------------------------------------------------------
# Files the window reads, and paths it prints
# ----------------------------------------------------------------------
GEOMETRY_FILE_TYPES = ["xyz", "gjf", "com", "gau", "fchk", "fch", "pdb", "ent", "mol", "sdf", "mol2"]
FCHK_SUFFIXES = (".fchk", ".fch")


def read_geometry_file(path) -> Tuple[str, Optional[int], Optional[int]]:
    """(Angstrom 'El x y z' text, charge or None, multiplicity or None) from a structure file.

    A Gaussian input or formatted checkpoint also gives the charge and multiplicity.
    """
    import molecularFileFormats as mff        # numpy + scipy: loaded when a file is opened
    path = Path(path)
    if path.suffix.lower() in FCHK_SUFFIXES:
        import gaussianWavefunction as gw
        try:
            sections = gw.read_fchk_sections(path.read_text(encoding="utf-8", errors="replace"))
            numbers = [int(z) for z in sections["Atomic numbers"]]
            xyz = list(sections["Current cartesian coordinates"])
        except (OSError, KeyError) as exc:
            raise ValueError(f"{path.name} is not a readable formatted checkpoint: {exc}") from exc
        atoms = [(fi.element_symbol(str(z)), *(mff.BOHR_RADIUS_ANGSTROM * float(v) for v in xyz[3 * i:3 * i + 3]))
                 for i, z in enumerate(numbers)]
        return "\n".join(fi.format_geometry(atoms)), int(sections.get("Charge", 0)), \
            int(sections.get("Multiplicity", 1))
    structure = mff.read_structure_file(path)
    atoms = [(element, *(float(v) for v in xyz)) for element, xyz in zip(structure.elements, structure.coordinates)]
    text = "\n".join(fi.format_geometry(atoms))
    charge = mult = None
    if path.suffix.lower() in GAUSSIAN_INPUT_SUFFIXES:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            tokens = line.split()
            if len(tokens) == 2 and all(t.lstrip("+-").isdigit() for t in tokens):
                charge, mult = int(tokens[0]), int(tokens[1])
                break
    return text, charge, mult


def wsl_path(path) -> str:
    """'C:\\jobs\\phenol' -> '/mnt/c/jobs/phenol' (the folder as WSL sees it)."""
    p = PureWindowsPath(str(path))
    if not p.drive or not p.drive.endswith(":"):
        return str(path)
    return "/mnt/" + p.drive[0].lower() + "/" + "/".join(p.parts[1:])
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_gaussianInput.py -q`
Expected: all pass (including the subprocess test: both imports are lazy, and neither module imports Toga).

- [ ] **Step 5: Commit**

```bash
git add gaussianInput.py tests/test_gaussianInput.py
git commit -m "feat(gaussian): geometry files and WSL paths for the input builder"
```

---

### Task 9: The window — panels, field table, live refresh (`gaussianInputBuilder.py`, part 1)

**Files:**
- Create: `gaussianInputBuilder.py`
- Test: `tests/test_gaussianInputBuilder.py` (create)

**Interfaces:**
- Consumes: `gaussianInput` (all of Tasks 2–8), `gaussianKeywords`, `systemProperties.parse_integer(text, field_name, default=None)`, `help.HelpGqteaWin.help_gaussian_input` (added in Task 11; until then the Help tab uses `getattr(HelpGqteaWin, "help_gaussian_input", "")`, and Task 11 replaces it with the plain attribute).
- Produces: module constants `CHOICE, TEXT, BOOL, INT, OPTIONAL_INT`, `STEP_WIDGETS: Dict[str, Tuple[str, str]]` (field → (widget attribute, kind)), `CONTROL_ATTRIBUTES: Tuple[str, ...]` (every other widget the methods touch), `TAB_TITLES`, `PREVIEW_TAB`, `_assign(widget, attr, value)`; class `GaussianInputBuilderUI` with `_init_state()`, `layout_main_window()`, `_collect_step() -> GaussianStepSpec`, `_load_step(step)`, `_collect_inputs() -> GaussianJobSpec`, `on_field_change(widget)`, `refresh()`, `refresh_controls(step)`, and (Task 10) the action handlers.

Layout: route line and step bar on top, a `toga.OptionContainer` with the 12 tabs of `TAB_TITLES` (each panel in a `ScrollContainer`), then Additional Keywords, job name and output folder, the buttons, a status line. Window 940 × 720.

- [ ] **Step 1: Write the failing tests**

`tests/test_gaussianInputBuilder.py`:

```python
"""The Gaussian 16 input builder window, headless (stub widgets, no window)."""
import asyncio
from dataclasses import fields, replace
from pathlib import Path
from types import SimpleNamespace

import pytest

import fukuiInput as fi
import gaussianInput as gi
import gaussianInputBuilder as gb
import gaussianKeywords as gk

SRC = Path(__file__).resolve().parents[1]
WATER = ("O 0.000000 0.000000 0.117300\n"
         "H 0.000000 0.757200 -0.469200\n"
         "H 0.000000 -0.757200 -0.469200")


class Recorder(SimpleNamespace):
    """A stub widget that counts writes, to check the write-on-change rule."""
    def __setattr__(self, name, value):
        if name != "writes":
            object.__setattr__(self, "writes", getattr(self, "writes", 0) + 1)
        object.__setattr__(self, name, value)


def make_ui(tmp_path):
    ui = gb.GaussianInputBuilderUI.__new__(gb.GaussianInputBuilderUI)
    ui._init_state()
    for attr, _ in gb.STEP_WIDGETS.values():
        setattr(ui, attr, SimpleNamespace(value=None, items=[], enabled=True))
    for attr in gb.CONTROL_ATTRIBUTES:
        setattr(ui, attr, SimpleNamespace(value="", text="", items=[], enabled=True, current_tab=0))
    ui.formchk_switch.value = True
    ui.stem_input.value, ui.output_dir_input.value, ui.geometry_input.value = "water", str(tmp_path), WATER
    ui._load_step(gi.GaussianStepSpec(multiplicity=1))
    ui.refresh()
    ui.dialogs, ui.answers = [], []

    async def record(title, message):
        ui.dialogs.append((title, message))

    async def confirm(title, message):
        ui.dialogs.append((title, message))
        return ui.answers.pop(0) if ui.answers else True
    ui._show_error = ui._show_info = record
    ui._confirm = confirm
    return ui


def change(ui, attr, value):
    """Set a widget as the user would, then fire its on_change."""
    widget = getattr(ui, attr)
    widget.value = value
    ui.on_field_change(widget)


def non_default_step() -> gi.GaussianStepSpec:
    default = gi.GaussianStepSpec()
    values = {}
    for f in fields(default):
        _, kind = gb.STEP_WIDGETS[f.name]
        value = getattr(default, f.name)
        if kind == gb.BOOL:
            values[f.name] = not value
        elif kind == gb.INT:
            values[f.name] = value + 2
        elif kind == gb.OPTIONAL_INT:
            values[f.name] = 7
        elif kind == gb.TEXT:
            values[f.name] = value + "x"
        else:
            values[f.name] = next(c for c in gi.choices(f.name, default) if c != value)
    values["family"], values["method_name"] = gk.FAMILY_SEMI, "PM6"
    return gi.GaussianStepSpec(**values)


class TestFieldTable:
    def test_every_step_field_has_a_widget(self):
        assert set(gb.STEP_WIDGETS) == {f.name for f in fields(gi.GaussianStepSpec)}

    def test_widget_attributes_are_unique(self):
        attrs = [attr for attr, _ in gb.STEP_WIDGETS.values()] + list(gb.CONTROL_ATTRIBUTES)
        assert len(attrs) == len(set(attrs))

    def test_every_field_round_trips_through_the_widgets(self, tmp_path):
        ui = make_ui(tmp_path)
        step = non_default_step()
        default = gi.GaussianStepSpec()
        assert all(getattr(step, f.name) != getattr(default, f.name) for f in fields(step))
        ui._load_step(step)
        assert ui._collect_step() == step

    def test_blank_numbers(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.charge_input.value = ui.multiplicity_input.value = ui.opt_max_cycles_input.value = ""
        step = ui._collect_step()
        assert (step.charge, step.multiplicity, step.opt_max_cycles) == (0, None, None)

    @pytest.mark.parametrize("attr,value,label", [("charge_input", "0.5", "Charge"),
                                                  ("nproc_input", "x", "%NProcShared"),
                                                  ("fukui_plus_mult_input", "2.5", "Multiplicity of N\\+p")])
    def test_bad_numbers_name_the_field(self, tmp_path, attr, value, label):
        ui = make_ui(tmp_path)
        getattr(ui, attr).value = value
        with pytest.raises(ValueError, match=label):
            ui._collect_step()

    def test_assign_writes_only_on_change(self):
        widget = Recorder(value="a")
        widget.writes = 0
        gb._assign(widget, "value", "a")
        gb._assign(widget, "value", "b")
        assert widget.writes == 1 and widget.value == "b"


class TestRefresh:
    def test_route_and_summary(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT_FREQ)
        assert ui.route_label.text == "#P B3LYP/6-31G(d) Opt Freq"
        assert ui.summary_label.text == "3 atoms, charge 0, 10 electrons, multiplicity 1"

    def test_a_bad_number_shows_in_the_route_line(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "charge_input", "x")
        assert ui.route_label.text.startswith("Cannot build the route yet:")

    def test_family_change_offers_its_methods(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "family_selection", gk.FAMILY_SEMI)
        assert ui.method_selection.items == list(gi.method_choices(gk.FAMILY_SEMI))
        assert ui.method_selection.value == "PM7" and ui.route_label.text == "#P PM7"

    def test_choosing_fukui_resets_what_the_fukui_job_sets(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.guess_selection.value = "Read checkpoint"
        ui.title_input.value = "t"
        change(ui, "job_type_selection", gk.JOB_FUKUI)
        assert (ui.guess_selection.value, ui.title_input.value) == ("Default", "")
        assert ui.scf_tight_switch.value and ui.nosymm_switch.value
        assert ui.route_label.text.startswith("job 1:     #P B3LYP/6-31G(d) SCF=Tight NoSymm Opt")

    def test_changes_while_syncing_are_ignored(self, tmp_path):
        ui = make_ui(tmp_path)
        ui._syncing = True
        change(ui, "job_type_selection", gk.JOB_OPT)
        assert ui.route_label.text == "#P B3LYP/6-31G(d)"

    def test_step_list_follows_the_steps(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT)
        assert ui.step_selection.items == ["1: Optimization B3LYP/6-31G(d)"]
        assert ui.step_selection.value == "1: Optimization B3LYP/6-31G(d)"


class TestControls:
    @pytest.mark.parametrize("changes,enabled,disabled", [
        (dict(), ["job_type_selection", "basis_selection", "diffuse_selection"],
         ["opt_target_selection", "irc_direction_selection", "fukui_p_input", "n_states_input",
          "method_other_input", "basis_other_input", "solvent_selection", "remove_step_button"]),
        (dict(job_type=gk.JOB_OPT), ["opt_target_selection", "modredundant_input", "modred_add_button"],
         ["freq_raman_switch"]),
        (dict(job_type=gk.JOB_FUKUI), ["fukui_p_input", "force_constants_selection"],
         ["guess_selection", "title_input", "chk_input", "state_selection", "spin_selection",
          "additional_input_input", "add_step_button", "opt_target_selection"]),
        (dict(state=gk.STATE_TD), ["n_states_input", "coefficients_switch"], []),
        (dict(basis="cc-pVDZ"), [], ["diffuse_selection", "polarization_selection"]),
        (dict(family=gk.FAMILY_SEMI, method_name="PM6"), [], ["basis_selection"]),
        (dict(method_name=gk.OTHER, solvation="SMD", solvent=gk.OTHER),
         ["method_other_input", "solvent_selection", "solvent_other_input"], []),
        (dict(geometry_source=gk.GEOM_ALLCHECK), [], ["connectivity_switch", "charge_input"]),
    ])
    def test_enabled_controls(self, tmp_path, changes, enabled, disabled):
        ui = make_ui(tmp_path)
        ui.refresh_controls(replace(gi.GaussianStepSpec(), **changes))
        assert [a for a in enabled if not getattr(ui, a).enabled] == []
        assert [a for a in disabled if getattr(ui, a).enabled] == []

    def test_fukui_note(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.refresh_controls(gi.GaussianStepSpec(job_type=gk.JOB_FUKUI))
        assert "three linked jobs" in ui.fukui_note_label.text
        ui.refresh_controls(gi.GaussianStepSpec())
        assert ui.fukui_note_label.text == ""


def test_window_module_imports_toga_but_the_core_does_not():
    import subprocess, sys
    code = "import sys, gaussianInputBuilder; print('toga' in sys.modules, 'cubeTab' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=str(SRC))
    assert out.stdout.strip() == "True False", out.stderr
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_gaussianInputBuilder.py -q`
Expected: `ModuleNotFoundError: No module named 'gaussianInputBuilder'`.

- [ ] **Step 3: Write `gaussianInputBuilder.py`** (Task 10 adds the action handlers; the buttons here already name them, so add the Task 10 methods as one-line stubs `async def preview(self, widget): pass` etc. in this task — `browse_output_dir`, `load_geometry`, `open_input`, `preview`, `save_input`, `defaults` async; `on_step_select`, `add_linked_step`, `remove_step`, `move_step_up`, `move_step_down`, `add_modredundant_line`, `closeTopLevel` sync — and replace them in Task 10)

```python
"""
Inputs > Gaussian 16 Input Builder: the window over ``gaussianInput``.

Modelled on GaussView 6.1's Gaussian Calculation Setup: a live route line, one tab per panel, an
Additional Keywords field, and a step bar for linked jobs (--Link1--). The window only moves values
between widgets and ``GaussianStepSpec``s (through the STEP_WIDGETS table) and turns every
ValueError into a dialog; validation, the input text, reopening and the Fukui job live in
``gaussianInput``. Handlers are synchronous and ``_syncing``-guarded, and widgets are written only
when a value changes (every write re-lays out the window on WinForms).
"""

import dataclasses
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import gaussianInput as gi
import gaussianKeywords as gk
from help import HelpGqteaWin
from systemProperties import parse_integer

CHOICE, TEXT, BOOL, INT, OPTIONAL_INT = "choice", "text", "bool", "int", "optional_int"

# GaussianStepSpec field -> (widget attribute, kind). _collect_step and _load_step go through this
# table, and a test checks it names every field, so a new spec field without a control fails the suite.
STEP_WIDGETS: Dict[str, Tuple[str, str]] = {
    "job_type": ("job_type_selection", CHOICE),
    "opt_target": ("opt_target_selection", CHOICE),
    "force_constants": ("force_constants_selection", CHOICE),
    "opt_tight": ("opt_tight_switch", BOOL),
    "opt_max_cycles": ("opt_max_cycles_input", OPTIONAL_INT),
    "freq_raman": ("freq_raman_switch", BOOL),
    "freq_vcd": ("freq_vcd_switch", BOOL),
    "freq_anharmonic": ("freq_anharmonic_switch", BOOL),
    "irc_direction": ("irc_direction_selection", CHOICE),
    "irc_force_constants": ("irc_force_constants_selection", CHOICE),
    "irc_max_points": ("irc_max_points_input", OPTIONAL_INT),
    "irc_step_size": ("irc_step_size_input", OPTIONAL_INT),
    "stable_opt": ("stable_opt_switch", BOOL),
    "nmr_method": ("nmr_method_selection", CHOICE),
    "nmr_spinspin": ("nmr_spinspin_switch", BOOL),
    "fukui_p": ("fukui_p_input", INT),
    "fukui_q": ("fukui_q_input", INT),
    "fukui_plus_multiplicity": ("fukui_plus_mult_input", OPTIONAL_INT),
    "fukui_minus_multiplicity": ("fukui_minus_mult_input", OPTIONAL_INT),
    "fukui_guess_read": ("fukui_guess_read_switch", BOOL),
    "state": ("state_selection", CHOICE),
    "n_states": ("n_states_input", INT),
    "root": ("root_input", INT),
    "state_spin": ("state_spin_selection", CHOICE),
    "print_all_coefficients": ("coefficients_switch", BOOL),
    "spin": ("spin_selection", CHOICE),
    "family": ("family_selection", CHOICE),
    "method_name": ("method_selection", CHOICE),
    "method_other": ("method_other_input", TEXT),
    "basis": ("basis_selection", CHOICE),
    "basis_diffuse": ("diffuse_selection", CHOICE),
    "basis_polarization": ("polarization_selection", CHOICE),
    "basis_other": ("basis_other_input", TEXT),
    "dispersion": ("dispersion_selection", CHOICE),
    "geometry_source": ("geometry_source_selection", CHOICE),
    "charge": ("charge_input", INT),
    "multiplicity": ("multiplicity_input", OPTIONAL_INT),
    "modredundant": ("modredundant_input", TEXT),
    "title": ("title_input", TEXT),
    "memory": ("memory_input", TEXT),
    "nproc": ("nproc_input", INT),
    "write_chk": ("write_chk_switch", BOOL),
    "chk": ("chk_input", TEXT),
    "old_chk": ("old_chk_input", TEXT),
    "extra_link0": ("extra_link0_input", TEXT),
    "print_level": ("print_level_selection", CHOICE),
    "scf_tight": ("scf_tight_switch", BOOL),
    "scf_qc": ("scf_qc_switch", BOOL),
    "nosymm": ("nosymm_switch", BOOL),
    "polar": ("polar_switch", BOOL),
    "connectivity": ("connectivity_switch", BOOL),
    "guess": ("guess_selection", CHOICE),
    "guess_mix": ("guess_mix_switch", BOOL),
    "guess_always": ("guess_always_switch", BOOL),
    "guess_only": ("guess_only_switch", BOOL),
    "guess_save": ("guess_save_switch", BOOL),
    "pop_level": ("pop_level_selection", CHOICE),
    "pop_mk": ("pop_mk_switch", BOOL),
    "pop_hirshfeld": ("pop_hirshfeld_switch", BOOL),
    "pop_nbo": ("pop_nbo_switch", BOOL),
    "pop_npa": ("pop_npa_switch", BOOL),
    "pop_savento": ("pop_savento_switch", BOOL),
    "solvation": ("solvation_selection", CHOICE),
    "solvent": ("solvent_selection", CHOICE),
    "solvent_other": ("solvent_other_input", TEXT),
    "additional_input": ("additional_input_input", TEXT),
    "additional_keywords": ("additional_keywords_input", TEXT),
}

# Every other widget the handlers touch (job-wide fields, labels, step bar, helper row, outputs).
CONTROL_ATTRIBUTES = (
    "stem_input", "output_dir_input", "geometry_input", "formchk_switch", "route_label", "summary_label",
    "status_label", "fukui_note_label", "step_selection", "add_step_button", "remove_step_button",
    "step_up_button", "step_down_button", "load_geometry_button", "modred_type_selection", "modred_atoms_input",
    "modred_action_selection", "modred_steps_input", "modred_size_input", "modred_add_button",
    "preview_output", "messages_output", "tabs",
)

TAB_TITLES = ("Job Type", "Method", "Molecule", "Title", "Link 0", "General", "Guess", "Pop", "Solvation",
              "Add. Inp.", "Preview", "Help")
PREVIEW_TAB = TAB_TITLES.index("Preview")
INPUT_FILE_TYPES = [suffix.lstrip(".") for suffix in gi.GAUSSIAN_INPUT_SUFFIXES]
FUKUI_NOTE = ("The Fukui job writes three linked jobs itself: the title, %Chk/%OldChk, the guess, Add. Inp., "
              "connectivity, the spin prefix and the electronic state are set by it.")


def _assign(widget, attr: str, value) -> None:
    """Write a widget property only when it changes (on WinForms every write re-lays out the window)."""
    if getattr(widget, attr, None) != value:
        setattr(widget, attr, value)


def _item_names(selection) -> List[str]:
    """The strings a Toga Selection holds (its rows carry them as ``.value``)."""
    return [getattr(item, "value", item) for item in selection.items]


class GaussianInputBuilderUI:
    LABEL_WIDTH = 190
    HINT_COLOR = "#666666"

    def __init__(self, *args, **kwargs):
        # gqteaWinToga passes the Button that was pressed; ignore it.
        self._init_state()
        self.window = toga.Window(title="Gaussian 16 Input Builder", size=(940, 720))
        self.layout_main_window()
        self.window.show()

    def _init_state(self):
        """Everything but widgets (shared with the headless tests)."""
        self.builder = gi.GaussianInputBuilder()
        self.steps: List[gi.GaussianStepSpec] = [gi.GaussianStepSpec()]
        self.current_step = 0
        self.job: Optional[gi.GaussianJob] = None
        self._syncing = False
        self._step_labels: List[str] = []

    # --- widget factories -------------------------------------------------
    def _label(self, text, width=None):
        return toga.Label(text, style=Pack(width=width or self.LABEL_WIDTH, margin=(6, 6, 0, 0)))

    def _row(self, *children):
        return toga.Box(children=list(children), style=Pack(direction=ROW, margin_bottom=4))

    def _heading(self, text):
        return toga.Label(text, style=Pack(font_weight="bold", margin=(10, 0, 4, 0)))

    def _hint(self, text):
        return toga.Label(text, style=Pack(color=self.HINT_COLOR, font_size=9, margin=(0, 0, 4, 0)))

    def _panel(self, *children):
        return toga.ScrollContainer(content=toga.Box(children=list(children), style=Pack(direction=COLUMN, margin=10)),
                                    style=Pack(flex=1))

    def _field(self, label, widget, hint=None):
        return self._row(self._label(label), widget, *([self._hint("   " + hint)] if hint else []))

    def _text(self, attr, width=None, placeholder=""):
        widget = toga.TextInput(placeholder=placeholder, style=Pack(width=width) if width else Pack(flex=1),
                                on_change=self.on_field_change)
        setattr(self, attr, widget)
        return widget

    def _multiline(self, attr, height, placeholder=""):
        widget = toga.MultilineTextInput(placeholder=placeholder, on_change=self.on_field_change,
                                         style=Pack(flex=1, height=height, font_family="monospace"))
        setattr(self, attr, widget)
        return widget

    def _choice(self, attr, items: Sequence[str], width=240):
        widget = toga.Selection(items=list(items), style=Pack(width=width), on_change=self.on_field_change)
        setattr(self, attr, widget)
        return widget

    def _switch(self, attr, text):
        widget = toga.Switch(text, on_change=self.on_field_change, style=Pack(margin_right=12))
        setattr(self, attr, widget)
        return widget

    def _checkbox_row(self, *switches):
        return self._row(self._label(""), *switches)

    # --- layout -------------------------------------------------------------
    def layout_main_window(self):
        self._syncing = True          # widgets fire on_change while they are built
        self.route_label = toga.Label("", style=Pack(font_family="monospace", margin=(4, 0, 4, 0), flex=1))
        self.summary_label = toga.Label("", style=Pack(font_family="monospace", margin=(4, 0, 4, 0)))
        self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
        self.fukui_note_label = toga.Label("", style=Pack(color=self.HINT_COLOR, margin=(6, 0, 4, 0)))

        self.step_selection = toga.Selection(items=["1"], style=Pack(width=520), on_change=self.on_step_select)
        self.add_step_button = toga.Button("Add linked step", on_press=self.add_linked_step, style=Pack(margin_left=6))
        self.remove_step_button = toga.Button("Remove", on_press=self.remove_step, style=Pack(margin_left=6))
        self.step_up_button = toga.Button("Up", on_press=self.move_step_up, style=Pack(margin_left=6, width=50))
        self.step_down_button = toga.Button("Down", on_press=self.move_step_down, style=Pack(margin_left=6, width=60))
        header = toga.Box(children=[
            self._row(toga.Label("Step:", style=Pack(margin=(6, 6, 0, 0))), self.step_selection,
                      self.add_step_button, self.remove_step_button, self.step_up_button, self.step_down_button),
            self._row(toga.Label("Route:", style=Pack(margin=(4, 6, 0, 0))), self.route_label),
        ], style=Pack(direction=COLUMN))

        self.preview_output = toga.MultilineTextInput(readonly=True, style=Pack(flex=2, font_family="monospace"))
        self.messages_output = toga.MultilineTextInput(readonly=True, style=Pack(flex=1, font_family="monospace"))
        preview = toga.Box(children=[self.preview_output, self.messages_output],
                           style=Pack(direction=COLUMN, margin=10, flex=1))
        help_box = toga.MultilineTextInput(readonly=True, style=Pack(flex=1),
                                           value=getattr(HelpGqteaWin, "help_gaussian_input", ""))
        panels = (self._build_job_panel(), self._build_method_panel(), self._build_molecule_panel(),
                  self._build_title_panel(), self._build_link0_panel(), self._build_general_panel(),
                  self._build_guess_panel(), self._build_pop_panel(), self._build_solvation_panel(),
                  self._build_additional_input_panel(), preview, help_box)
        self.tabs = toga.OptionContainer(content=list(zip(TAB_TITLES, panels)), style=Pack(flex=1))

        self.stem_input = toga.TextInput(value=gi.fi.DEFAULT_STEM, style=Pack(width=220), on_change=self.on_field_change)
        self.output_dir_input = toga.TextInput(placeholder="Required: folder for the .gjf and .sh files",
                                               style=Pack(flex=1), on_change=self.on_field_change)
        footer = toga.Box(children=[
            self._row(self._label("Additional Keywords:"), self._text("additional_keywords_input",
                                                                      placeholder="e.g. Int=UltraFine 6D")),
            self._hint("Any keyword or option the panels do not have. Options of a keyword a panel writes are "
                       "merged into it (Opt=NoEigenTest)."),
            self._row(self._label("Job name (file stem):"), self.stem_input, self._label("   Output folder:", 120),
                      self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
            self._row(*(toga.Button(text, on_press=handler, style=Pack(width=120, margin=5)) for text, handler in (
                ("Open input...", self.open_input), ("Preview", self.preview), ("Save input", self.save_input),
                ("Defaults", self.defaults), ("Close", self.closeTopLevel)))),
            self.status_label,
        ], style=Pack(direction=COLUMN))
        self.window.content = toga.Box(children=[header, self.tabs, footer],
                                       style=Pack(direction=COLUMN, margin=10, flex=1))
        self._syncing = False
        self._load_step(self.steps[0])
        self.refresh()

    def _build_job_panel(self):
        return self._panel(
            self._field("Job type:", self._choice("job_type_selection", gk.JOB_TYPES, 280)),
            self._heading("Optimization"),
            self._field("Optimize to a:", self._choice("opt_target_selection", gk.OPT_TARGETS)),
            self._field("Force constants:", self._choice("force_constants_selection", gk.FORCE_CONSTANTS, 280)),
            self._checkbox_row(self._switch("opt_tight_switch", "Tight convergence")),
            self._field("Maximum cycles:", self._text("opt_max_cycles_input", 100, "default")),
            self._hint("Scan and freeze lines are entered in the Molecule panel."),
            self._heading("Frequency"),
            self._checkbox_row(self._switch("freq_raman_switch", "Raman"), self._switch("freq_vcd_switch", "VCD"),
                               self._switch("freq_anharmonic_switch", "Anharmonic")),
            self._heading("IRC"),
            self._field("Direction:", self._choice("irc_direction_selection", gk.IRC_DIRECTIONS)),
            self._field("Force constants:", self._choice("irc_force_constants_selection", gk.FORCE_CONSTANTS, 280)),
            self._field("MaxPoints:", self._text("irc_max_points_input", 100, "default")),
            self._field("StepSize:", self._text("irc_step_size_input", 100, "default"), "units of 0.01 bohr"),
            self._heading("Stability"),
            self._checkbox_row(self._switch("stable_opt_switch", "Reoptimize the wavefunction if unstable (Stable=Opt)")),
            self._heading("NMR"),
            self._field("Method:", self._choice("nmr_method_selection", gk.NMR_METHODS)),
            self._checkbox_row(self._switch("nmr_spinspin_switch", "Spin-spin couplings")),
            self._heading("Fukui functions (N, N+p, N-q)"),
            self._field("p (electrons added):", self._text("fukui_p_input", 80)),
            self._field("q (electrons removed):", self._text("fukui_q_input", 80)),
            self._field("Multiplicity of N+p:", self._text("fukui_plus_mult_input", 80, "required")),
            self._field("Multiplicity of N-q:", self._text("fukui_minus_mult_input", 80, "required")),
            self._checkbox_row(self._switch("fukui_guess_read_switch",
                                            "Start the ions from the optimized N-state orbitals (Guess=Read)")),
            self._hint("Job 1 optimizes the N state (Optimization options above, minimum only); jobs 2-3 are single "
                       "points with Geom=Check. Writes <name>_fukui.gjf and <name>_formchk.sh."),
            self.fukui_note_label)

    def _build_method_panel(self):
        return self._panel(
            self._field("Electronic state:", self._choice("state_selection", gk.STATES)),
            self._field("Number of states:", self._text("n_states_input", 80)),
            self._field("State of interest (Root):", self._text("root_input", 80)),
            self._field("States:", self._choice("state_spin_selection", gk.STATE_SPINS)),
            self._checkbox_row(self._switch("coefficients_switch", "Print all excitation coefficients for the NTO "
                                                                   "and Excitation tabs (IOp(9/40=4))")),
            self._heading("Model chemistry"),
            self._field("Spin:", self._choice("spin_selection", gk.SPIN_PREFIXES)),
            self._field("Method family:", self._choice("family_selection", gk.FAMILIES)),
            self._field("Method:", self._choice("method_selection", gi.method_choices(gk.FAMILY_DFT))),
            self._field("Other method:", self._text("method_other_input", 200), "when Method is Other..."),
            self._field("Basis set:", self._choice("basis_selection", gk.BASES)),
            self._row(self._label("Diffuse, polarization:"), self._choice("diffuse_selection", gk.DIFFUSE, 90),
                      self._choice("polarization_selection", gk.POLARIZATION, 130),
                      self._hint("   6-31G and 6-311G only")),
            self._field("Other basis set:", self._text("basis_other_input", 200), "when Basis set is Other..."),
            self._field("Empirical dispersion:", self._choice("dispersion_selection", gk.DISPERSION)))

    def _build_molecule_panel(self):
        self.load_geometry_button = toga.Button("Load geometry...", on_press=self.load_geometry, style=Pack(width=150))
        self.geometry_input = toga.MultilineTextInput(
            placeholder="Cartesian geometry, Angstrom:\nC   0.000000   0.000000   0.000000\n...",
            style=Pack(flex=1, height=170, font_family="monospace"), on_change=self.on_field_change)
        self.modred_type_selection = toga.Selection(items=list(gk.MODRED_TYPES), style=Pack(width=100))
        self.modred_atoms_input = toga.TextInput(placeholder="atoms, e.g. 1 2 3 4", style=Pack(width=150))
        self.modred_action_selection = toga.Selection(items=list(gk.MODRED_ACTIONS), style=Pack(width=90))
        self.modred_steps_input = toga.TextInput(placeholder="steps", style=Pack(width=60))
        self.modred_size_input = toga.TextInput(placeholder="step size", style=Pack(width=80))
        self.modred_add_button = toga.Button("Add", on_press=self.add_modredundant_line, style=Pack(width=60))
        return self._panel(
            self._row(self.load_geometry_button,
                      self._hint("   .xyz .gjf .com .fchk .mol .sdf .mol2 .pdb; a .gjf/.fchk also sets the charge "
                                 "and multiplicity")),
            self.geometry_input,
            self._hint("The geometry belongs to the whole job: step 1 writes it, and a later step can write it again."),
            self._field("Geometry source:", self._choice("geometry_source_selection", gk.GEOMETRY_SOURCES, 380)),
            self._field("Charge:", self._text("charge_input", 80)),
            self._field("Multiplicity:", self._text("multiplicity_input", 80, "required"),
                        "never assumed; the line below lists the allowed values"),
            self.summary_label,
            self._heading("Scan and freeze coordinates (Opt=ModRedundant)"),
            self._row(self.modred_type_selection, self.modred_atoms_input, self.modred_action_selection,
                      self.modred_steps_input, self.modred_size_input, self.modred_add_button),
            self._hint("Atoms numbered from 1. Scan: number of steps and step size (Angstrom or degrees). "
                       "Typed lines: B/A/D, the atoms, then F, or S steps size."),
            self._multiline("modredundant_input", 90, "B 1 2 S 10 0.05\nD 1 2 3 4 F"))

    def _build_title_panel(self):
        return self._panel(
            self._field("Title:", self._text("title_input", placeholder="blank = the job name")),
            self._hint("One line. Gaussian does not allow @ # ! - _ \\ in the title."))

    def _build_link0_panel(self):
        return self._panel(
            self._field("%Mem:", self._text("memory_input", 120), "integer + unit, e.g. 2GB, 500MB"),
            self._field("%NProcShared:", self._text("nproc_input", 120), "cores of the machine that runs Gaussian"),
            self._checkbox_row(self._switch("write_chk_switch", "Write a checkpoint (%Chk)")),
            self._field("%Chk:", self._text("chk_input", 220, "<name>.chk"), "file name, no folders"),
            self._field("%OldChk:", self._text("old_chk_input", 220, "optional"),
                        "an existing checkpoint to read (Guess=Read, Geom=Check, ReadFC)"),
            self._label("Extra Link 0 lines:"),
            self._multiline("extra_link0_input", 70, "%KJob L502"),
            self._hint("Relative file names: the job runs under Linux from the input's folder."))

    def _build_general_panel(self):
        self.formchk_switch = toga.Switch("Write a formchk script (<name>_formchk.sh) for the 3D viewer's tabs",
                                          value=True, on_change=self.on_field_change)
        return self._panel(
            self._field("Print level:", self._choice("print_level_selection", gk.PRINT_LEVELS, 100)),
            self._checkbox_row(self._switch("scf_tight_switch", "SCF=Tight"), self._switch("scf_qc_switch", "SCF=QC"),
                               self._switch("nosymm_switch", "NoSymm"), self._switch("polar_switch", "Polar")),
            self._checkbox_row(self._switch("connectivity_switch", "Write connectivity (Geom=Connectivity)")),
            self._checkbox_row(self.formchk_switch))

    def _build_guess_panel(self):
        return self._panel(
            self._field("Initial guess:", self._choice("guess_selection", gk.GUESS_METHODS)),
            self._checkbox_row(self._switch("guess_mix_switch", "Mix HOMO and LUMO"),
                               self._switch("guess_always_switch", "Always generate in optimizations")),
            self._checkbox_row(self._switch("guess_only_switch", "Only do guess (no SCF)"),
                               self._switch("guess_save_switch", "Save orbitals")))

    def _build_pop_panel(self):
        return self._panel(
            self._field("Population output:", self._choice("pop_level_selection", gk.POP_LEVELS)),
            self._checkbox_row(self._switch("pop_mk_switch", "MK charges"),
                               self._switch("pop_hirshfeld_switch", "Hirshfeld (CM5) charges")),
            self._checkbox_row(self._switch("pop_nbo_switch", "NBO"), self._switch("pop_npa_switch", "NPA"),
                               self._switch("pop_savento_switch", "Save NTOs")))

    def _build_solvation_panel(self):
        return self._panel(
            self._field("Solvation model:", self._choice("solvation_selection", gk.SOLVATION_MODELS)),
            self._field("Solvent:", self._choice("solvent_selection", gk.SOLVENTS)),
            self._field("Other solvent:", self._text("solvent_other_input", 200), "Gaussian's name, one word"))

    def _build_additional_input_panel(self):
        return self._panel(
            self._hint("Input sections written after the geometry (and after any connectivity and scan/freeze "
                       "lines): a Gen basis set, ECPs, ... Separate sections with a blank line."),
            self._multiline("additional_input_input", 320))

    # --- panels <-> GaussianStepSpec -----------------------------------------
    def _set_items(self, selection, items: Sequence[str]) -> None:
        if _item_names(selection) != list(items):
            selection.items = list(items)

    def _collect_step(self) -> gi.GaussianStepSpec:
        """The panels as a GaussianStepSpec; ValueError (naming the field) for a non-integer number."""
        defaults = gi.GaussianStepSpec()
        values = {}
        for name, (attr, kind) in STEP_WIDGETS.items():
            widget = getattr(self, attr)
            if kind == BOOL:
                values[name] = bool(widget.value)
            elif kind == CHOICE:
                values[name] = widget.value if widget.value is not None else getattr(defaults, name)
            elif kind == TEXT:
                values[name] = widget.value or ""
            elif kind == INT:
                values[name] = parse_integer(widget.value, gi.FIELD_LABELS[name], getattr(defaults, name))
            else:
                values[name] = parse_integer(widget.value, gi.FIELD_LABELS[name])
        return gi.GaussianStepSpec(**values)

    def _load_step(self, step: gi.GaussianStepSpec) -> None:
        self._syncing = True
        try:
            self._set_items(self.method_selection, gi.method_choices(step.family))
            for name, (attr, kind) in STEP_WIDGETS.items():
                value = getattr(step, name)
                if kind in (INT, OPTIONAL_INT):
                    value = "" if value is None else str(value)
                _assign(getattr(self, attr), "value", value)
        finally:
            self._syncing = False

    def _collect_inputs(self) -> gi.GaussianJobSpec:
        self.steps[self.current_step] = self._collect_step()
        return gi.GaussianJobSpec(stem=(self.stem_input.value or "").strip(),
                                  output_dir=(self.output_dir_input.value or "").strip(),
                                  geometry=self.geometry_input.value or "",
                                  steps=[dataclasses.replace(s) for s in self.steps],
                                  write_formchk=bool(self.formchk_switch.value))

    # --- live updates ---------------------------------------------------------
    def on_field_change(self, widget):
        """Every control (synchronous): keep the panels consistent, then refresh the live texts."""
        if self._syncing:
            return
        if widget is self.family_selection:
            self._syncing = True
            try:
                choices = gi.method_choices(widget.value)
                self._set_items(self.method_selection, choices)
                _assign(self.method_selection, "value", gk.DEFAULT_METHOD.get(widget.value, choices[0]))
            finally:
                self._syncing = False
        elif widget is self.job_type_selection and widget.value == gk.JOB_FUKUI:
            try:
                step = self._collect_step()
            except ValueError:
                step = dataclasses.replace(self.steps[self.current_step], job_type=gk.JOB_FUKUI)
            self._load_step(gi.fukui_defaults(step))
        self.refresh()

    def refresh(self) -> None:
        try:
            spec = self._collect_inputs()
        except ValueError as exc:
            _assign(self.route_label, "text", f"Cannot build the route yet: {exc}")
            return
        _assign(self.route_label, "text", gi.route_preview(spec, self.current_step))
        _assign(self.summary_label, "text", gi.molecule_summary(spec, self.current_step))
        self._refresh_step_list(spec)
        self.refresh_controls(spec.steps[self.current_step])

    def _refresh_step_list(self, spec: gi.GaussianJobSpec) -> None:
        labels = [gi.step_label(i + 1, s) for i, s in enumerate(spec.steps)]
        if labels == self._step_labels:
            return
        self._syncing = True
        try:
            self._set_items(self.step_selection, labels)
            _assign(self.step_selection, "value", labels[self.current_step])
        finally:
            self._syncing = False
        self._step_labels = labels

    def refresh_controls(self, step: gi.GaussianStepSpec) -> None:
        """Enable only what applies to this step (GaussView greys out the same way)."""
        job, fukui = step.job_type, step.job_type == gk.JOB_FUKUI
        optimizing = job in gk.OPTIMIZING_JOBS
        excited = step.state != gk.STATE_GROUND
        model_state = step.state in gk.MODEL_STATES
        basis_used = not gi.no_basis(step)
        coordinates = step.geometry_source == gk.GEOM_COORDINATES
        solvated = step.solvation != gk.SOLVATION_NONE
        rules = {
            "opt_target_selection": optimizing, "force_constants_selection": optimizing or fukui,
            "opt_tight_switch": optimizing or fukui, "opt_max_cycles_input": optimizing or fukui,
            "freq_raman_switch": job in gk.FREQUENCY_JOBS, "freq_vcd_switch": job in gk.FREQUENCY_JOBS,
            "freq_anharmonic_switch": job in gk.FREQUENCY_JOBS,
            "irc_direction_selection": job == gk.JOB_IRC, "irc_force_constants_selection": job == gk.JOB_IRC,
            "irc_max_points_input": job == gk.JOB_IRC, "irc_step_size_input": job == gk.JOB_IRC,
            "stable_opt_switch": job == gk.JOB_STABILITY,
            "nmr_method_selection": job == gk.JOB_NMR, "nmr_spinspin_switch": job == gk.JOB_NMR,
            "fukui_p_input": fukui, "fukui_q_input": fukui, "fukui_plus_mult_input": fukui,
            "fukui_minus_mult_input": fukui, "fukui_guess_read_switch": fukui,
            "state_selection": not fukui, "n_states_input": excited and not fukui,
            "root_input": excited and not fukui, "state_spin_selection": excited and not fukui,
            "coefficients_switch": excited and not fukui, "spin_selection": not fukui,
            "family_selection": not model_state,
            "method_selection": not model_state and step.family in gk.FAMILY_METHODS,
            "method_other_input": not model_state and step.family == gk.FAMILY_DFT and step.method_name == gk.OTHER,
            "basis_selection": basis_used,
            "diffuse_selection": basis_used and step.basis in gk.EXTENSIBLE_POPLE,
            "polarization_selection": basis_used and step.basis in gk.EXTENSIBLE_POPLE,
            "basis_other_input": basis_used and step.basis == gk.OTHER,
            "geometry_source_selection": not fukui,
            "charge_input": step.geometry_source != gk.GEOM_ALLCHECK,
            "multiplicity_input": step.geometry_source != gk.GEOM_ALLCHECK,
            "geometry_input": self.current_step == 0 or coordinates,
            "load_geometry_button": self.current_step == 0 or coordinates,
            "modredundant_input": optimizing, "modred_add_button": optimizing,
            "modred_type_selection": optimizing, "modred_atoms_input": optimizing,
            "modred_action_selection": optimizing, "modred_steps_input": optimizing,
            "modred_size_input": optimizing,
            "title_input": not fukui, "write_chk_switch": not fukui, "chk_input": not fukui and step.write_chk,
            "old_chk_input": not fukui, "extra_link0_input": not fukui,
            "connectivity_switch": not fukui and coordinates,
            "guess_selection": not fukui, "guess_mix_switch": not fukui, "guess_always_switch": not fukui,
            "guess_only_switch": not fukui, "guess_save_switch": not fukui,
            "solvent_selection": solvated, "solvent_other_input": solvated and step.solvent == gk.OTHER,
            "additional_input_input": not fukui,
            "add_step_button": not fukui, "remove_step_button": len(self.steps) > 1,
            "step_up_button": self.current_step > 0, "step_down_button": self.current_step < len(self.steps) - 1,
        }
        for attr, enabled in rules.items():
            _assign(getattr(self, attr), "enabled", enabled)
        _assign(self.fukui_note_label, "text", FUKUI_NOTE if fukui else "")
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_gaussianInputBuilder.py -q`
Expected: all pass. If `test_window_module_imports_toga_but_the_core_does_not` reports `cubeTab` loaded, a module this file imports pulled it in; remove that import.

The window is opened by hand in Task 11, once its button exists.

- [ ] **Step 5: Commit**

```bash
git add gaussianInputBuilder.py tests/test_gaussianInputBuilder.py
git commit -m "feat(gaussian): input builder window with GaussView-style panels"
```

---

### Task 10: The window — step bar and actions (`gaussianInputBuilder.py`, part 2)

**Files:**
- Modify: `gaussianInputBuilder.py` (replace the Task 9 stubs)
- Test: `tests/test_gaussianInputBuilder.py` (append; includes the tests migrated from `tests/test_fukuiInputBuilder.py`)

**Interfaces:**
- Consumes: Task 9 window; `gi.linked_step`, `gi.modredundant_line`, `gi.read_geometry_file`, `gi.read_input_file`, `gi.wsl_path`, `gi.fukui_spec`, `gi.GEOMETRY_FILE_TYPES`.
- Produces: sync `on_step_select(widget)`, `add_linked_step(widget)`, `remove_step(widget)`, `move_step_up(widget)`, `move_step_down(widget)`, `add_modredundant_line(widget)`, `apply_geometry_file(path)`, `apply_input_file(path) -> List[str]` (notes), `closeTopLevel(widget)`; async `browse_output_dir`, `load_geometry`, `open_input`, `preview`, `save_input`, `defaults`; dialog helpers `_show_error`, `_show_info`, `_confirm` (overridden in tests).

- [ ] **Step 1: Write the failing tests** (append)

```python
FUKUI_INPUTS = SRC / "tests" / "data" / "fukui" / "inputs"


def load_phenol(ui, ions=True):
    ui.apply_geometry_file(FUKUI_INPUTS / "phenol_N.gjf")
    change(ui, "job_type_selection", gk.JOB_FUKUI)
    if ions:
        ui.fukui_plus_mult_input.value = ui.fukui_minus_mult_input.value = "2"
        ui.refresh()


class TestSteps:
    def test_add_linked_step_continues_from_the_checkpoint(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT)
        ui.add_linked_step(None)
        assert len(ui.steps) == 2 and ui.current_step == 1
        assert ui.geometry_source_selection.value == gk.GEOM_ALLCHECK and ui.guess_selection.value == "Read checkpoint"
        assert ui.step_selection.items == ["1: Optimization B3LYP/6-31G(d)",
                                           "2: Energy B3LYP/6-31G(d), Geom=AllCheck"]
        assert ui.route_label.text == "#P B3LYP/6-31G(d) Geom=AllCheck Guess=Read"

    def test_switching_steps_stores_and_loads_the_panels(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.add_linked_step(None)
        change(ui, "state_selection", gk.STATE_TD)
        ui.step_selection.value = ui._step_labels[0]
        ui.on_step_select(ui.step_selection)
        assert ui.current_step == 0 and ui.state_selection.value == gk.STATE_GROUND
        assert ui.steps[1].state == gk.STATE_TD

    def test_a_step_with_a_bad_number_is_not_left(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.add_linked_step(None)
        ui.n_states_input.value = "x"
        ui.step_selection.value = ui._step_labels[0]
        ui.on_step_select(ui.step_selection)
        assert ui.current_step == 1 and "Fix step 2" in ui.status_label.text
        assert ui.step_selection.value == ui._step_labels[1]

    def test_remove_and_move(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.remove_step(None)
        assert len(ui.steps) == 1 and "only step" in ui.status_label.text
        ui.add_linked_step(None)
        ui.move_step_up(None)
        assert ui.current_step == 0 and ui.steps[0].geometry_source == gk.GEOM_ALLCHECK
        ui.move_step_down(None)
        assert ui.current_step == 1 and ui.steps[1].geometry_source == gk.GEOM_ALLCHECK
        ui.remove_step(None)
        assert len(ui.steps) == 1 and ui.current_step == 0

    def test_a_fukui_step_takes_no_linked_step(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_FUKUI)
        ui.add_linked_step(None)
        assert len(ui.steps) == 1 and "three linked jobs" in ui.status_label.text

    def test_helper_row_appends_a_line(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.modredundant_input.value = "B 1 2 F"
        ui.modred_type_selection.value, ui.modred_atoms_input.value = "Angle", "2 1 3"
        ui.modred_action_selection.value, ui.modred_steps_input.value, ui.modred_size_input.value = "Scan", "5", "2.0"
        ui.add_modredundant_line(None)
        assert ui.modredundant_input.value == "B 1 2 F\nA 2 1 3 S 5 2.0"
        ui.modred_atoms_input.value = "1"
        ui.add_modredundant_line(None)
        assert "needs 3 atom numbers" in ui.status_label.text


class TestGeometryAndFiles:
    def test_apply_fills_the_form_but_not_the_ion_spins(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.stem_input.value = fi.DEFAULT_STEM
        ui.output_dir_input.value = ""
        ui.fukui_plus_mult_input.value = ui.fukui_minus_mult_input.value = "4"
        ui.apply_geometry_file(FUKUI_INPUTS / "phenol_Np1.gjf")
        assert ui.charge_input.value == "-1" and ui.multiplicity_input.value == "2"
        assert ui.fukui_plus_mult_input.value == "" and ui.fukui_minus_mult_input.value == ""
        assert ui.stem_input.value == "phenol_Np1" and ui.output_dir_input.value == str(FUKUI_INPUTS)

    def test_apply_from_an_n_fchk_drops_the_state_suffix(self, tmp_path):
        from fukui_fixtures import N
        ui = make_ui(tmp_path)
        ui.stem_input.value = ""
        ui.apply_geometry_file(N)
        assert ui.stem_input.value == "phenol" and ui.output_dir_input.value == str(tmp_path)

    def test_apply_keeps_a_typed_job_name(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.stem_input.value = "mine"
        ui.apply_geometry_file(FUKUI_INPUTS / "phenol_N.gjf")
        assert ui.stem_input.value == "mine"

    def test_open_input_restores_every_step(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT)
        ui.add_linked_step(None)
        change(ui, "state_selection", gk.STATE_TD)
        asyncio.run(ui.save_input(None))
        fresh = make_ui(tmp_path)
        notes = fresh.apply_input_file(tmp_path / "water.gjf")
        assert notes == [] and len(fresh.steps) == 2 and fresh.current_step == 0
        assert fresh.job_type_selection.value == gk.JOB_OPT and fresh.steps[1].state == gk.STATE_TD
        assert fresh.geometry_input.value == "\n".join(fi.format_geometry(fi.parse_geometry(WATER)))

    def test_open_input_errors_go_to_a_dialog(self, tmp_path):
        ui = make_ui(tmp_path)
        bad = tmp_path / "bad.gjf"
        bad.write_text("nothing\n")
        ui._ask_open_path = lambda: _resolved(bad)
        asyncio.run(ui.open_input(None))
        assert ui.dialogs and ui.dialogs[0][0] == "Open Input Error" and "bad.gjf" in ui.dialogs[0][1]


async def _resolved(value):
    return value


class TestActions:
    def test_preview_fills_the_input_text(self, tmp_path):
        ui = make_ui(tmp_path)
        asyncio.run(ui.preview(None))
        assert ui.preview_output.value.startswith("%Chk=water.chk") and not ui.dialogs
        assert ui.tabs.current_tab == gb.PREVIEW_TAB

    def test_preview_errors_go_to_a_dialog(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.multiplicity_input.value = ""
        asyncio.run(ui.preview(None))
        assert ui.dialogs[0][0] == "Input Error" and "it is not assumed" in ui.dialogs[0][1]

    def test_save_writes_the_files_and_says_how_to_run_them(self, tmp_path):
        ui = make_ui(tmp_path)
        asyncio.run(ui.save_input(None))
        assert (tmp_path / "water.gjf").is_file() and (tmp_path / "water_formchk.sh").is_file()
        assert "g16 < water.gjf > water.log" in ui.messages_output.value
        assert "bash water_formchk.sh" in ui.messages_output.value

    def test_existing_files_are_kept_unless_confirmed(self, tmp_path):
        ui = make_ui(tmp_path)
        (tmp_path / "water.gjf").write_text("mine")
        ui.answers = [False]
        asyncio.run(ui.save_input(None))
        assert (tmp_path / "water.gjf").read_text() == "mine"
        ui.answers = [True]
        asyncio.run(ui.save_input(None))
        assert "%Chk=water.chk" in (tmp_path / "water.gjf").read_text()

    def test_defaults_ask_first_and_keep_the_molecule(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.add_linked_step(None)
        ui.answers = [False]
        asyncio.run(ui.defaults(None))
        assert len(ui.steps) == 2
        ui.answers = [True]
        asyncio.run(ui.defaults(None))
        assert ui.steps == [gi.GaussianStepSpec(multiplicity=1)] and ui.geometry_input.value == WATER
        ui.dialogs.clear()
        asyncio.run(ui.defaults(None))
        assert ui.dialogs == []          # nothing changed: no question


class TestFukuiInTheWindow:
    """Migrated from tests/test_fukuiInputBuilder.py: the Fukui builder now lives here as a Job Type."""

    def test_every_fukui_spec_field_is_read_from_a_widget(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        ui.stem_input.value, ui.memory_input.value, ui.nproc_input.value = "ph", "3GB", "3"
        ui.method_selection.value, ui.basis_selection.value = "M06-2X", "Def2SVP"
        ui.additional_keywords_input.value = "SCF=XQC"
        ui.opt_tight_switch.value = True
        ui.print_level_selection.value = "#N"
        ui.charge_input.value, ui.multiplicity_input.value = "1", "2"
        ui.fukui_p_input.value, ui.fukui_q_input.value = "2", "3"
        ui.fukui_plus_mult_input.value, ui.fukui_minus_mult_input.value = "3", "5"
        ui.fukui_guess_read_switch.value = False
        spec = gi.fukui_spec(ui._collect_inputs())
        defaults = fi.FukuiJobSpec("", "", "", 0, None, None, None)
        for f in fields(fi.FukuiJobSpec):
            assert getattr(spec, f.name) != getattr(defaults, f.name), f"{f.name} is not read from the window"

    def test_blank_multiplicities_are_not_filled_in(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui, ions=False)
        spec = gi.fukui_spec(ui._collect_inputs())
        assert (spec.multiplicity, spec.plus_multiplicity, spec.minus_multiplicity) == (1, None, None)

    def test_summary_and_route_list_the_three_states(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        text = ui.summary_label.text
        assert "N+1" in text and "charge -1" in text and "51 electrons" in text and "phenol_Np1.chk" in text
        assert ui.route_label.text == ("job 1:     #P B3LYP/6-31G(d) SCF=Tight NoSymm Opt\n"
                                       "jobs 2-3:  #P B3LYP/6-31G(d) SCF=Tight NoSymm Geom=Check Guess=Read")

    def test_preview_asks_for_missing_spins(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui, ions=False)
        asyncio.run(ui.preview(None))
        assert ui.dialogs[0][0] == "Input Error"
        assert "N+1 state: enter its spin multiplicity; it is not assumed" in ui.dialogs[0][1]

    def test_preview_shows_three_linked_jobs(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        asyncio.run(ui.preview(None))
        assert ui.preview_output.value.count("--Link1--") == 2 and not ui.dialogs

    def test_save_writes_both_files_and_says_how_to_run_them(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.stem_input.value = fi.DEFAULT_STEM
        load_phenol(ui)
        ui.output_dir_input.value = str(tmp_path)
        asyncio.run(ui.save_input(None))
        assert (tmp_path / "phenol_fukui.gjf").is_file() and (tmp_path / "phenol_formchk.sh").is_file()
        assert "g16 < phenol_fukui.gjf > phenol_fukui.log" in ui.messages_output.value
        assert "bash phenol_formchk.sh" in ui.messages_output.value
        assert "Analysis of Fukui Functions" in ui.messages_output.value
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_gaussianInputBuilder.py -q`
Expected: the new tests fail (the Task 9 stubs do nothing).

- [ ] **Step 3: Implement** (replace the Task 9 stubs in `GaussianInputBuilderUI`)

```python
    # --- dialogs ------------------------------------------------------------
    async def _show_error(self, title: str, message: str):
        await self.window.dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        await self.window.dialog(toga.InfoDialog(title, message))

    async def _confirm(self, title: str, message: str) -> bool:
        return bool(await self.window.dialog(toga.QuestionDialog(title, message)))

    async def _ask_open_path(self):
        return await self.window.dialog(toga.OpenFileDialog(title="Open a Gaussian input",
                                                            file_types=INPUT_FILE_TYPES))

    def _status(self, text: str) -> None:
        _assign(self.status_label, "text", text)

    # --- steps --------------------------------------------------------------
    def _store_current(self) -> bool:
        """Keep the panels in self.steps; False (and a status message) when a number does not parse."""
        try:
            self.steps[self.current_step] = self._collect_step()
            return True
        except ValueError as exc:
            self._status(f"Fix step {self.current_step + 1} first: {exc}")
            return False

    def _show_step(self, index: int) -> None:
        self.current_step = index
        self._step_labels = []                       # force the list (and its selection) to be rewritten
        self._load_step(self.steps[index])
        self.refresh()

    def on_step_select(self, widget):
        if self._syncing or widget.value not in self._step_labels:
            return
        index = self._step_labels.index(widget.value)
        if index == self.current_step:
            return
        if not self._store_current():
            self._syncing = True
            try:
                _assign(self.step_selection, "value", self._step_labels[self.current_step])
            finally:
                self._syncing = False
            return
        self._show_step(index)

    def add_linked_step(self, widget):
        if not self._store_current():
            return
        step = self.steps[self.current_step]
        if step.job_type == gk.JOB_FUKUI:
            self._status("The Fukui functions job type is already three linked jobs; it cannot take more steps.")
            return
        self.steps.insert(self.current_step + 1, gi.linked_step(step))
        self._show_step(self.current_step + 1)
        self._status(f"Step {self.current_step + 1} added: it reads the geometry and orbitals of the checkpoint.")

    def remove_step(self, widget):
        if len(self.steps) == 1:
            self._status("This is the only step; there is nothing to remove.")
            return
        del self.steps[self.current_step]
        self._show_step(min(self.current_step, len(self.steps) - 1))

    def _move_step(self, delta: int) -> None:
        target = self.current_step + delta
        if not 0 <= target < len(self.steps) or not self._store_current():
            return
        self.steps[self.current_step], self.steps[target] = self.steps[target], self.steps[self.current_step]
        self._show_step(target)

    def move_step_up(self, widget):
        self._move_step(-1)

    def move_step_down(self, widget):
        self._move_step(+1)

    # --- molecule -----------------------------------------------------------
    def add_modredundant_line(self, widget):
        try:
            line = gi.modredundant_line(self.modred_type_selection.value, self.modred_atoms_input.value or "",
                                        self.modred_action_selection.value, self.modred_steps_input.value or "",
                                        self.modred_size_input.value or "")
        except ValueError as exc:
            self._status(str(exc))
            return
        text = (self.modredundant_input.value or "").rstrip()
        self.modredundant_input.value = f"{text}\n{line}" if text else line
        self._status(f"Added: {line}")
        self.refresh()

    async def load_geometry(self, widget):
        try:
            path = await self.window.dialog(toga.OpenFileDialog(title="Geometry of the molecule",
                                                                file_types=gi.GEOMETRY_FILE_TYPES))
            if path:
                self.apply_geometry_file(path)
        except ValueError as exc:
            await self._show_error("Geometry Error", str(exc))
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))

    def apply_geometry_file(self, path):
        """Fill the geometry (and the charge/multiplicity a .gjf/.fchk gives); keep a typed job name.

        The Fukui ion multiplicities are cleared, never guessed: they belong to the previous molecule.
        """
        path = Path(path)
        text, charge, mult = gi.read_geometry_file(path)
        self._syncing = True
        try:
            _assign(self.geometry_input, "value", text)
            if charge is not None:
                _assign(self.charge_input, "value", str(charge))
            if mult is not None:
                _assign(self.multiplicity_input, "value", str(mult))
            _assign(self.fukui_plus_mult_input, "value", "")
            _assign(self.fukui_minus_mult_input, "value", "")
            if (self.stem_input.value or "").strip() in ("", gi.fi.DEFAULT_STEM):
                stem = path.stem[:-2] if path.stem.endswith("_N") else path.stem
                _assign(self.stem_input, "value", "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in stem))
            if not (self.output_dir_input.value or "").strip():
                _assign(self.output_dir_input, "value", str(path.resolve().parent))
        finally:
            self._syncing = False
        self._status(f"Geometry loaded from {path.name}.")
        self.refresh()

    async def browse_output_dir(self, widget):
        try:
            folder = await self.window.dialog(toga.SelectFolderDialog(title="Folder for the Gaussian input"))
            if folder:
                self.output_dir_input.value = str(folder)
        except Exception as exc:
            await self._show_error("Folder Selection Error", str(exc))

    # --- whole job ------------------------------------------------------------
    def apply_input_file(self, path) -> List[str]:
        """Open input...: every step of a Gaussian input into the panels; returns the reader's notes."""
        spec, notes = gi.read_input_file(path)
        self.steps = spec.steps
        self._syncing = True
        try:
            _assign(self.stem_input, "value", spec.stem)
            _assign(self.output_dir_input, "value", spec.output_dir)
            _assign(self.geometry_input, "value", spec.geometry)
            _assign(self.formchk_switch, "value", spec.write_formchk)
        finally:
            self._syncing = False
        self._show_step(0)
        _assign(self.messages_output, "value", "\n".join(f"NOTE: {n}" for n in notes) or "Read without notes.")
        self._status(f"{Path(path).name}: {len(spec.steps)} step(s) read into the panels.")
        return notes

    async def open_input(self, widget):
        try:
            path = await self._ask_open_path()
            if path:
                self.apply_input_file(path)
        except ValueError as exc:
            await self._show_error("Open Input Error", str(exc))
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))

    def _prepare(self) -> gi.GaussianJob:
        return self.builder.prepare(self._collect_inputs())

    def _warnings_text(self, job: gi.GaussianJob) -> str:
        return "\n".join(f"WARNING: {w}" for w in job.warnings) if job.warnings else "No warnings."

    def _run_instructions(self, job: gi.GaussianJob) -> str:
        gjf, log = job.input_path.name, job.input_path.with_suffix(".log").name
        lines = ["Written:", f"  {job.input_path}"] + ([f"  {job.script_path}"] if job.script_path else []) + [
            "",
            "Next, in a Linux shell with Gaussian 16 set up (WSL, or the cluster after copying the folder):",
            f"  cd '{gi.wsl_path(job.input_path.parent)}'      (WSL; on a cluster: the copied folder)",
            f"  g16 < {gjf} > {log}      (or submit it through the cluster's queue)"]
        if job.script_path:
            lines.append(f"  bash {job.script_path.name}")
        if job.fukui is not None:
            lines.append(f"Then open {job.fukui.states[0].fchk.name} in Tools > Analysis of Fukui Functions; the N+1 "
                         "and N-1 files are filled in automatically (copy the three .fchk back first if they were "
                         "made on a cluster).")
        elif job.script_path:
            lines.append("Then open the .fchk in the 3D Molecular Viewer (Orbitals, ESP, QTAIM, NTO, Excitation tabs).")
        return "\n".join(lines) + "\n\n" + self._warnings_text(job)

    async def preview(self, widget):
        try:
            job = self._prepare()
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        self.job = job
        self.preview_output.value = self.builder.input_text(job)
        self.messages_output.value = self._warnings_text(job)
        self.tabs.current_tab = PREVIEW_TAB
        self._status(f"Preview of {job.input_path.name} (not saved yet).")

    async def save_input(self, widget):
        try:
            job = self._prepare()
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        existing = [p.name for p in (job.input_path, job.script_path) if p is not None and p.exists()]
        if existing and not await self._confirm("Replace Files?", f"{', '.join(existing)} already exist in "
                                                f"{job.input_path.parent}. Replace them?"):
            return
        try:
            self.builder.write_files(job)
        except ValueError as exc:
            await self._show_error("Write Error", str(exc))
            return
        self.job = job
        self.preview_output.value = self.builder.input_text(job)
        self.messages_output.value = self._run_instructions(job)
        self.tabs.current_tab = PREVIEW_TAB
        self._status(f"Saved {job.input_path.name}" + (f" and {job.script_path.name}." if job.script_path else "."))

    async def defaults(self, widget):
        """Every step back to one default step; the molecule (geometry, charge, multiplicity) is kept."""
        self._store_current()
        first = self.steps[0]
        reference = gi.GaussianStepSpec(charge=first.charge, multiplicity=first.multiplicity)
        if self.steps != [reference] and not await self._confirm(
                "Reset to Defaults?", "Every panel of every step goes back to its default (the geometry, "
                                      "charge and multiplicity are kept). Continue?"):
            return
        self.steps = [reference]
        self._show_step(0)
        self._status("Panels reset to their defaults.")

    def closeTopLevel(self, widget):
        self.window.close()
```

`test_open_input_restores_every_step`: the saved file's first step has `%Chk=water.chk`, and the second step is `Geom=AllCheck`, so the reader produces no notes (the `%Mem` line is present in both steps).

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_gaussianInputBuilder.py tests/test_gaussianInput.py tests/test_gaussianInputReopen.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add gaussianInputBuilder.py tests/test_gaussianInputBuilder.py
git commit -m "feat(gaussian): step bar, open/preview/save/defaults in the input builder"
```

---

### Task 11: Menu entry, removal of the old window, help texts

**Files:**
- Modify: `toolLoader.py:49` and `toolLoader.py:108`, `help.py:2402-2404` and `help.py:2441-2488`, `fukuiInput.py:19`, `gaussianInputBuilder.py` (Help tab line), `gaussianUtilities.py:4`
- Delete: `fukuiInputBuilder.py`, `tests/test_fukuiInputBuilder.py`
- Test: `tests/test_toolLoader.py` (existing tests must pass), `tests/test_gaussianInputBuilder.py` (append one test)

**Interfaces:**
- Produces: `ToolSpec("Gaussian 16 Input Builder", "gaussianInputBuilder", "GaussianInputBuilderUI")` on the Inputs tab; `HelpGqteaWin.help_gaussian_input` (str, ends with the rewritten `help_fukui_input`).

- [ ] **Step 1: Write the failing test** (append to `tests/test_gaussianInputBuilder.py`)

```python
def test_the_inputs_tab_offers_the_builder_instead_of_the_fukui_window():
    from toolLoader import TOOL_TABS
    inputs = dict(TOOL_TABS)["Inputs"]
    assert any(spec.module == "gaussianInputBuilder" and spec.attribute == "GaussianInputBuilderUI"
               and spec.label == "Gaussian 16 Input Builder" for spec in inputs)
    assert not any(spec.module == "fukuiInputBuilder" for tab in dict(TOOL_TABS).values() for spec in tab)
    assert not (SRC / "fukuiInputBuilder.py").exists()


def test_help_text_covers_the_builder_and_the_fukui_job_type():
    from help import HelpGqteaWin
    text = HelpGqteaWin.help_gaussian_input
    assert text.startswith("GAUSSIAN 16 INPUT BUILDER") and "FUKUI FUNCTIONS JOB TYPE" in text
    assert "Gaussian Fukui Input Builder" not in text and "Gaussian Fukui Input Builder" not in HelpGqteaWin.help_fukui
```

Check the shape of `TOOL_TABS` first (`toolLoader.py:30-80`): it is a sequence of `(tab title, [ToolSpec, ...])` pairs, which `dict(...)` accepts. If it is a different shape, adapt the two lookups, not the assertions.

Run: `python -m pytest tests/test_gaussianInputBuilder.py -q -k "inputs_tab or help_text"` — Expected: FAIL.

- [ ] **Step 2: Replace the menu entry**

`toolLoader.py:49`:
```python
        ToolSpec("Gaussian 16 Input Builder", "gaussianInputBuilder", "GaussianInputBuilderUI"),
```
`toolLoader.py:108` (inside `_frozen_build_imports`): `import gqteaMDinputBuilder, fukuiInputBuilder` → `import gqteaMDinputBuilder, gaussianInputBuilder`.

- [ ] **Step 3: Delete the old window and its tests**

```bash
git rm fukuiInputBuilder.py tests/test_fukuiInputBuilder.py
```
(Every test of `tests/test_fukuiInputBuilder.py` has a counterpart: `test_wsl_path` and `TestGeometryFiles` in `tests/test_gaussianInput.py` (Task 8); the apply/collect/summary/action tests in `tests/test_gaussianInputBuilder.py` (Task 10).)

- [ ] **Step 4: Docstrings that named the old window**

`fukuiInput.py:19`: `No Toga here (``fukuiInputBuilder`` is the window).` → `No Toga here; the window is the Fukui functions job type of ``gaussianInputBuilder``.`
`gaussianUtilities.py:4`: `the Fukui input builder` → `the Gaussian 16 input builder`.

- [ ] **Step 5: Help texts in `help.py`**

Replace lines 2402–2404 of `help_fukui` with:
```
  Inputs > Gaussian 16 Input Builder, job type "Fukui functions", writes
  steps 1-3 for Gaussian 16 (one linked job: Opt on N, then Geom=Check
  single points on N+1 and N-1, plus a formchk script); choosing its
  <name>_N.fchk here fills in the other two.
```

Replace the whole `help_fukui_input` string (lines 2441–2488) with:
```python
    help_fukui_input = """FUKUI FUNCTIONS JOB TYPE
One Gaussian 16 input with three linked jobs (separated by --Link1--) for a
finite-difference Fukui calculation, plus a script that converts the
checkpoints for Tools > Analysis of Fukui Functions.

  job 1  N    Opt           charge Q      mult. M    %Chk=<name>_N.chk
  job 2  N+p  single point  charge Q - p  mult. M+   %OldChk=<name>_N.chk
                                                     %Chk=<name>_Np1.chk
  job 3  N-q  single point  charge Q + q  mult. M-   %OldChk=<name>_N.chk
                                                     %Chk=<name>_Nm1.chk

Job 1 optimizes the N state from the geometry of the Molecule panel. Jobs 2
and 3 read the OPTIMIZED geometry from <name>_N.chk (%OldChk + Geom=Check,
and by default Guess=Read for the orbitals) and are not re-optimized: the
Fukui functions are derivatives at fixed nuclei. %OldChk is copied into each
job's own %Chk, so <name>_N.chk is never overwritten.

WHAT YOU MUST SUPPLY (cannot be safely defaulted)
  Geometry       Molecule panel: starting structure of the N state.
  Charge of N    Molecule panel; blank = 0 (a .gjf or .fchk sets it).
  Multiplicities of N (Molecule panel), N+p and N-q (Job Type panel): ALL
                 REQUIRED. Spin states are never assumed; the Molecule panel
                 lists the values each electron count allows. Use the
                 lowest-energy spin state of each ion (closed-shell
                 neutral: 1, 2, 2; neutral radical: 2, 1, 1).
  %Mem / %NProcShared  Link 0 panel, for the machine that runs Gaussian.
  Output folder  where the .gjf and .sh files go.

DEFAULTS (editable)
  Choosing this job type switches SCF=Tight and NoSymm on (General panel):
  #P B3LYP/6-31G(d) SCF=Tight NoSymm in all jobs, Opt in job 1, %Mem=2GB,
  %NProcShared=4, p = q = 1, Guess=Read on. NoSymm keeps one orientation.
  Job 1's optimization takes the Optimization options of the Job Type
  panel (minimum only; force constants, Tight, MaxCycles). The title,
  %Chk/%OldChk, Guess, Add. Inp., connectivity, spin prefix and excited
  states are set by the job and greyed out. In Additional Keywords Opt,
  Scan, IRC, Geom=, Guess=, Units=, Gen/GenECP, ECPs and options that read
  extra input are refused; Freq, post-SCF methods (the analysis reads only
  the Total SCF Density) and diffuse basis sets warn. The Fukui job must be
  the only step.

FILES (Save input; relative names, LF line endings)
  <name>_fukui.gjf    the linked job.
  <name>_formchk.sh   formchk <name>_N.chk <name>_N.fchk, and so on.

RUNNING (in the input's folder; under WSL C:\\jobs is /mnt/c/jobs)
  g16 < <name>_fukui.gjf > <name>_fukui.log   (or the cluster's queue)
  bash <name>_formchk.sh
Then open <name>_N.fchk in Tools > Analysis of Fukui Functions: the N+1
and N-1 files are filled in. If one job fails, Gaussian does not run the
jobs after it; the script then reports which .chk is missing.
"""

    help_gaussian_input = """GAUSSIAN 16 INPUT BUILDER
Writes a Gaussian 16 input (Rev. C.01, under WSL or on a Linux cluster), in
the way of GaussView 6.1's Gaussian Calculation Setup: the route line at the
top follows every panel as you change it.

PANELS (of the selected step)
  Job Type   Energy; Optimization (minimum or TS (Berny), force constants,
             Tight, MaxCycles); Frequency (Raman, VCD, Anharmonic);
             Opt+Freq; IRC (direction, force constants, MaxPoints,
             StepSize); Scan (relaxed, Opt=ModRedundant); Stability
             (Stable=Opt); NMR (GIAO, CSGT, IGAIM, All; SpinSpin); Fukui
             functions (see below).
  Method     ground or excited state (TD-SCF, TDA, CIS, ZINDO: number of
             states, Root, Singlets/Triplets/50-50; "Print all excitation
             coefficients" writes IOp(9/40=4), which the NTO and Excitation
             tabs of the 3D viewer need), spin prefix R/U/RO, method family
             and method, basis set (6-31G and 6-311G take the diffuse and
             polarization lists), empirical dispersion.
  Molecule   geometry (Load geometry... or Cartesian Angstrom lines), where
             this step's geometry comes from, charge, multiplicity, scan and
             freeze lines (helper row: type, atoms, Scan/Freeze, steps, step
             size, Add; or type B/A/D lines).
  Title  Link 0 (%Mem, %NProcShared, %Chk, %OldChk, extra % lines)  General
  (print level, SCF=Tight, SCF=QC, NoSymm, Polar, connectivity, formchk
  script)  Guess  Pop  Solvation (IEFPCM, CPCM, SMD and the solvent)
  Add. Inp. (sections written after the geometry: a Gen basis set, ECPs).
  Additional Keywords: anything the panels do not have. Options of a
  keyword a panel already writes are merged into it (Opt=NoEigenTest);
  the same option with another value (MaxCycles=50 and 100) is an error.

WHAT YOU MUST SUPPLY
  Geometry; multiplicity (never assumed: the Molecule panel lists the
  values the electron count allows); %Mem and %NProcShared for the machine
  that runs Gaussian; the output folder.

LINKED STEPS (--Link1--)
  Add linked step copies the selected step and presets it to continue from
  the checkpoint: Geom=AllCheck, Guess=Read, job type Energy (for example an
  optimization, then a TD-SCF single point on the optimized geometry).
  Choose a step in the Step list to edit it; Remove, Up and Down change the
  list. Every step writes %Chk=<name>.chk unless you change it. A step that
  reads a checkpoint (Geom=Check or AllCheck, Guess=Read, ReadFC) needs one
  written by an earlier step, or named in %OldChk.

OPEN INPUT...
  Reads a .gjf/.com/.gau back into the panels, every step. Keywords no
  panel has go to Additional Keywords and extra sections to Add. Inp., so
  saving again keeps them. Z-matrices, per-fragment charges and freeze flags
  inside the coordinates are refused; notes list what was adapted (an
  absolute %Chk path, a title with characters Gaussian forbids).

FILES (Save input; relative names, LF line endings)
  <name>.gjf          the input (Fukui functions: <name>_fukui.gjf).
  <name>_formchk.sh   formchk of every checkpoint the job writes.

RUNNING (in the input's folder; under WSL C:\\jobs is /mnt/c/jobs)
  g16 < <name>.gjf > <name>.log   (or the cluster's queue)
  bash <name>_formchk.sh
Then open the .fchk in the 3D Molecular Viewer (Orbitals, ESP, QTAIM, NTO
and Excitation tabs).

""" + help_fukui_input
```

- [ ] **Step 6: The Help tab reads the attribute directly**

In `gaussianInputBuilder.py`, `value=getattr(HelpGqteaWin, "help_gaussian_input", "")` → `value=HelpGqteaWin.help_gaussian_input`.

- [ ] **Step 7: Run the tests**

Run: `python -m pytest tests/test_gaussianInputBuilder.py tests/test_toolLoader.py tests/test_fukuiInput.py -q`
Expected: all pass. Then `python -m pytest tests/ -q` — Expected: only the four baseline failures.

- [ ] **Step 8: Open the window by hand**

Run `python gqteaWinToga.py`, press **Inputs > Gaussian 16 Input Builder**. Check:
- the window opens at 940 × 720 with the step bar, route line, 12 tabs and the footer, without growing past the screen;
- Job type Optimization enables the Optimization rows and the route reads `#P B3LYP/6-31G(d) Opt`;
- Method family Semi-empirical lists PM7 first and greys out the basis set;
- Load geometry with `tests/data/fukui/inputs/phenol_N.gjf` fills charge 0 and multiplicity 1, and the summary line counts 50 electrons;
- Add linked step shows step 2 with Geom=AllCheck Guess=Read, and choosing step 1 in the list restores its panels;
- Job type Fukui functions greys out Title, Guess, %Chk and Add linked step;
- Close ends the window without errors in the console.

Report anything that does not look right before committing; fix it in this task.

- [ ] **Step 9: Commit**

```bash
git add toolLoader.py help.py fukuiInput.py gaussianUtilities.py gaussianInputBuilder.py tests/test_gaussianInputBuilder.py
git commit -m "feat(gaussian): Inputs tab opens the Gaussian 16 input builder; Fukui builder moved into it"
```

---

### Task 12: User documentation and project notes

**Files:**
- Modify: `USER_MANUAL.md` (TOC lines 45 and 139; lines 154 and 488; section lines 1001–1069; line 3334; section `### Fukui input builder messages`; line 7194 and the version history), `readme.md:10`, `readme.md:36` and `readme.md:92`, `docs/index.html:92`, `CLAUDE.md`, `gqteaWinToga.py:47`

This feature bumps the app version from **1.6.0 to 1.7.0** (user decision, 2026-11-03; Step 7b). The historical version-history entries (e.g. "Version 0.8.5 ... added the Gaussian Fukui Input Builder") stay as they are.

- [ ] **Step 1: USER_MANUAL table of contents**

Line 45: `    - [Gaussian Fukui Input Builder](#gaussian-fukui-input-builder)` →
```markdown
    - [Gaussian 16 Input Builder](#gaussian-16-input-builder)
```
Line 139: `  - [Fukui input builder messages](#fukui-input-builder-messages)` →
```markdown
  - [Gaussian 16 input builder messages](#gaussian-16-input-builder-messages)
```
Line 488: `and the Gaussian Fukui input builder.` → `and the Gaussian 16 input builder (which also writes the Fukui-function jobs).`
Line 3334: `[Inputs > Gaussian Fukui Input Builder](#gaussian-fukui-input-builder)` → `[Inputs > Gaussian 16 Input Builder](#fukui-functions-job-type), job type *Fukui functions*,`

- [ ] **Step 2: Replace the section `#### Gaussian Fukui Input Builder` (lines 1001–1069, up to the line before `### Structural Tools`)** with:

````markdown
#### Gaussian 16 Input Builder

**Inputs > Gaussian 16 Input Builder** prepares Gaussian 16 inputs in the way of GaussView 6.1's *Gaussian Calculation Setup*. The input is written for **Gaussian 16 (Rev. C.01) under Linux**, in WSL or on a computing cluster: relative file names and Linux (LF) line endings. Nothing is run from the window; it writes the files and lists the commands.

**Window.**
- **Step** list with **Add linked step**, **Remove**, **Up**, **Down**: the steps of a linked job (`--Link1--`). The panels always show the selected step.
- **Route**: the route section of the selected step, updated as you change any control.
- Tabs **Job Type**, **Method**, **Molecule**, **Title**, **Link 0**, **General**, **Guess**, **Pop**, **Solvation**, **Add. Inp.**, **Preview**, **Help**.
- **Additional Keywords**: any keyword or option the panels do not offer.
- **Job name** (file stem) and **Output folder**; buttons **Open input...**, **Preview**, **Save input**, **Defaults**, **Close**.

**Panels.**

| Panel | Controls | Writes |
|---|---|---|
| Job Type | Energy; Optimization (minimum or TS (Berny), force constants, tight convergence, maximum cycles); Frequency (Raman, VCD, anharmonic); Opt+Freq; IRC (direction, force constants, MaxPoints, StepSize); Scan (relaxed); Stability (reoptimize); NMR (GIAO, CSGT, IGAIM, All; spin-spin); Fukui functions | `Opt=(TS,CalcFC,Tight,MaxCycles=n)`, `Freq=(Raman,VCD,Anharmonic)`, `IRC=(...)`, `Opt=ModRedundant`, `Stable=Opt`, `NMR=(...)` |
| Method | electronic state (ground, TD-SCF, TDA, CIS, ZINDO; number of states, Root, Singlets/Triplets/50-50); *Print all excitation coefficients*; spin (default, R, U, RO); method family (Hartree-Fock, DFT, semi-empirical, MP2, MP4(SDTQ), CCSD, CCSD(T), compound) and method; basis set, with diffuse (`+`, `++`) and polarization lists for 6-31G and 6-311G; empirical dispersion | `UB3LYP/6-31+G(d,p)`, `TD=(NStates=10)`, `CIS(NStates=3)/cc-pVDZ`, `IOp(9/40=4)`, `EmpiricalDispersion=GD3BJ` |
| Molecule | geometry (**Load geometry...** or Cartesian `Element x y z` lines in Å), geometry source of the step, charge, multiplicity, scan and freeze lines | the molecule specification, `Geom=Check` / `Geom=AllCheck`, the ModRedundant section |
| Title | one line; blank = the job name | the title section |
| Link 0 | `%Mem`, `%NProcShared`, `%Chk` (on by default, `<name>.chk`), `%OldChk`, extra `%` lines | the Link 0 section |
| General | print level (`#P`, `#`, `#N`, `#T`), `SCF=Tight`, `SCF=QC`, `NoSymm`, `Polar`, *Write connectivity*, *Write a formchk script* | `SCF=(Tight,QC)`, `Geom=Connectivity` and its section |
| Guess | default, core, Huckel, read the checkpoint, read if possible (TCheck); mix HOMO/LUMO, always, only, save | `Guess=(Read,Mix)` |
| Pop | output level; MK, Hirshfeld (CM5), NBO, NPA, save NTOs | `Pop=(Full,MK)` |
| Solvation | none, IEFPCM, CPCM, SMD; solvent (Gaussian's names, or Other) | `SCRF=(SMD,Solvent=Water)` |
| Add. Inp. | free text | sections after the geometry: a `Gen` basis set, ECPs, … |

*Print all excitation coefficients* writes `IOp(9/40=4)`: without it Gaussian prints only the largest CI coefficients, and the 3D viewer's NTO and Excitation tabs report a normalisation deviation.

**What you must supply.** The geometry; the **multiplicity**, which is never assumed (a blank field stops Preview and Save with the values the electron count allows; the line under it shows them as you type); `%Mem` and `%NProcShared` for the machine that runs Gaussian; the output folder.

**Additional Keywords.** A keyword that no panel writes is added to the route. If a panel already writes the keyword, the options are **merged** into it: with *Optimize to a* TS and `Opt=NoEigenTest` here the route reads `Opt=(TS,NoEigenTest)`. The same option with two values (the panel's `MaxCycles=100` and `Opt=MaxCycles=50`) is an error that names both places. The method and basis set belong in the Method panel.

**Linked steps.** **Add linked step** copies the selected step and presets it to continue from the checkpoint: geometry source *Checkpoint: everything* (`Geom=AllCheck`), `Guess=Read`, job type Energy. A typical job is an optimization followed by a TD-SCF single point. Every step writes `%Chk=<name>.chk` unless you change it, so each step reads what the previous one wrote. A step that reads a checkpoint (`Geom=Check`, `Geom=AllCheck`, `Guess=Read`, `ReadFC`) must find one written by an earlier step, or named in its `%OldChk`. Steps that share a checkpoint leave only the last step's results in it, and a warning says so. The geometry belongs to the whole job: a later step with the source *Coordinates* writes the same geometry again.

**Scan and freeze coordinates.** In the Molecule panel, choose Bond, Angle or Dihedral, type the atoms (numbered from 1), choose Scan (number of steps and step size, in Å or degrees) or Freeze, and press **Add**; or type lines such as `B 1 2 S 10 0.05` or `D 1 2 3 4 F`. Scan lines need the job type *Scan (relaxed)*; freeze lines also work with Optimization and Opt+Freq.

**Open input...** reads a `.gjf`, `.com` or `.gau` file back into the panels, every step. Keywords no panel has go to Additional Keywords and extra sections to Add. Inp., so saving again keeps them. An absolute `%Chk` path keeps only its file name, and title characters Gaussian does not allow (`@ # ! - _ \`) become blanks; the Preview tab lists these notes. Z-matrices, per-fragment charges and freeze flags inside the coordinates are refused.

**Buttons.**
- **Preview** shows the input on the Preview tab, with any warnings.
- **Save input** writes `<name>.gjf` and, when *Write a formchk script* is on, `<name>_formchk.sh`, asking before it replaces existing files. The Preview tab then lists the commands to run.
- **Defaults** resets every panel of every step (it asks first); the geometry, charge and multiplicity are kept.

**Running the job.**
1. In a Linux shell with Gaussian 16 set up, go to the folder that holds the input. Under WSL, `C:\jobs\water` is `/mnt/c/jobs/water`. On a cluster, copy the folder there first.
2. Run `g16 < <name>.gjf > <name>.log`, or submit the job through the cluster's queue.
3. Run `bash <name>_formchk.sh` to convert every checkpoint to a `.fchk`, which the 3D Molecular Viewer's Orbitals, ESP, QTAIM, NTO and Excitation tabs read.

##### Fukui functions job type

The job type **Fukui functions (N, N+p, N−q)** prepares the three Gaussian calculations that [Tools > Analysis of Fukui Functions](#analysis-of-fukui-functions) needs, as **one** input file with three linked jobs:

| Job | State | Calculation | Charge | Multiplicity | Checkpoints |
|---|---|---|---|---|---|
| 1 | N | `Opt` | Q | M | `%Chk=<name>_N.chk` |
| 2 | N+p | single point, `Geom=Check Guess=Read` | Q − p | M₊ | `%OldChk=<name>_N.chk`, `%Chk=<name>_Np1.chk` |
| 3 | N−q | single point, `Geom=Check Guess=Read` | Q + q | M₋ | `%OldChk=<name>_N.chk`, `%Chk=<name>_Nm1.chk` |

Job 1 optimizes the N state from the geometry of the Molecule panel. Jobs 2 and 3 read the **optimized** geometry (and, by default, the orbitals) from `<name>_N.chk` and are **not** re-optimized, because the Fukui functions are derivatives at fixed nuclei. Gaussian copies the `%OldChk` file into each job's own `%Chk`, so `<name>_N.chk` keeps the optimized N state. All three jobs share the method, basis set and keywords of the panels.

- **Job Type panel**: p and q (1 for ordinary f⁺/f⁻; larger for degenerate frontier orbitals), the multiplicities of N+p and N−q, and *Start the ions from the optimized N-state orbitals (Guess=Read)*. The Optimization options apply to job 1: a minimum only, with default, CalcFC or CalcAll force constants, tight convergence and a maximum number of cycles.
- **Molecule panel**: the geometry, the charge and the multiplicity of the N state. The summary line lists the charge, electron count, `.chk` name and multiplicity of the three states. **All three multiplicities are required**: for ordinary organic molecules the lowest states are usually 1, 2, 2 for a closed-shell neutral and 2, 1, 1 for a neutral radical, but you must enter them. For transition-metal systems, compute several spin states of each ion and use the lowest in energy.
- Choosing this job type switches `SCF=Tight` and `NoSymm` on (General panel). The title, `%Chk`/`%OldChk`, Guess, Add. Inp., connectivity, the spin prefix and excited states are set by the job and greyed out. It must be the only step.
- In Additional Keywords `Opt`, `IRC`, `Scan`, `Geom=`, `Guess=`, `Units=`, `Gen`/`GenECP`, `ExtraBasis`, `Pseudo`/ECPs and options that read an extra input section are refused; `Freq`, post-SCF methods, `TD` and diffuse basis sets only warn (the analysis reads only the `Total SCF Density`).

**Save input** writes `<name>_fukui.gjf` and `<name>_formchk.sh`. Run `g16 < <name>_fukui.gjf > <name>_fukui.log`, then `bash <name>_formchk.sh`, which runs `formchk <name>_N.chk <name>_N.fchk` and the same for the other two, and reports any missing checkpoint. Copy the three `.fchk` files back if they were made on a cluster, and choose `<name>_N.fchk` in *Tools > Analysis of Fukui Functions*: the other two files are filled in automatically, and the analysis checks that the three states share the same geometry, basis set and method. If one linked job fails, Gaussian does not run the jobs after it; the script then reports which `.chk` is missing, and the reason is in the `.log` file.

````

- [ ] **Step 3: Replace `### Fukui input builder messages`** heading with `### Gaussian 16 input builder messages`, change the second row's message/cause/fix to:

```markdown
| `The Fukui functions job type: set the N-state optimization with the Optimization options of the Job Type panel, not in Additional Keywords.` | `Opt` was typed in Additional Keywords for a Fukui job. | Use the Optimization options of the Job Type panel. |
```

and insert these rows after the first row:

```markdown
| `Step 1: Multiplicity: enter it; it is not assumed (10 electrons allow 1, 3, 5, ...)` | The multiplicity field of that step is blank. | Enter the multiplicity of the state you want. |
| `Step 2: Geom=AllCheck reads a checkpoint, but no earlier step writes x.chk; set %OldChk to an existing checkpoint.` | The step reads a checkpoint that no earlier step writes. | Give the steps the same `%Chk`, or name an existing checkpoint in `%OldChk`. |
| `Additional Keywords: Opt=MaxCycles=50 sets MaxCycles=50, but the Job Type panel sets MaxCycles=100.` | The same option was given two values. | Keep one: change the panel or remove the keyword. |
| `Basis set Gen: put the basis-set section in the Add. Inp. panel.` | `Gen`/`GenECP` without its basis-set section. | Type the basis set (and ECPs) in Add. Inp. |
| `Step 1: '0 1 0 1 0 1' is not one charge and one multiplicity (per-fragment charges are not supported).` | Open input... met a fragment-charge line. | Edit the file to one charge and multiplicity, or set up that job by hand. |
```

- [ ] **Step 4: `readme.md`**

Line 36: `- Gaussian Fukui input builder (one Gaussian 16 linked job: ...)` →
`- Gaussian 16 input builder modelled on GaussView 6.1's calculation setup (job types, methods, basis sets, solvation, linked steps, reopening of .gjf files), including the Fukui-function linked job (optimization of N, then `Geom=Check` single points on N+1 and N-1; plus a formchk script)`
Line 92: `fukuiInput.py` and `fukuiInputBuilder.py` for the Gaussian Fukui input builder` → `` `gaussianKeywords.py`, `gaussianInput.py` and `gaussianInputBuilder.py` for the Gaussian 16 input builder (`fukuiInput.py` writes its Fukui job) ``

- [ ] **Step 5: `docs/index.html:92`**

`<li>Gaussian 16 Fukui input builder: the N, N+1 and N&minus;1 calculations in one linked job, plus a <code>formchk</code> script <a href="#fukui">(new)</a></li>` →
`<li>Gaussian 16 input builder in the style of GaussView: job types, methods, solvation, linked steps and the Fukui N, N+1 and N&minus;1 calculations in one linked job, plus a <code>formchk</code> script</li>`

- [ ] **Step 6: `CLAUDE.md`**

Replace the bullet that starts `- **[fukuiInput.py](fukuiInput.py) / [fukuiInputBuilder.py](fukuiInputBuilder.py)** — **Inputs > Gaussian Fukui Input Builder**` by a bullet for the new tool, and keep the Fukui facts under it:

```markdown
- **[gaussianKeywords.py](gaussianKeywords.py) / [gaussianInput.py](gaussianInput.py) / [gaussianInputBuilder.py](gaussianInputBuilder.py) + [fukuiInput.py](fukuiInput.py)** — **Inputs > Gaussian 16 Input Builder**, GaussView 6.1's Gaussian Calculation Setup (its help pages in `C:\G09W\help`, `gausscalcopts.htm`, list each control's keyword). Spec: [docs/superpowers/specs/2026-11-03-gaussian-input-builder-design.md](docs/superpowers/specs/2026-11-03-gaussian-input-builder-design.md); plan: [docs/superpowers/plans/2026-11-03-gaussian-input-builder.md](docs/superpowers/plans/2026-11-03-gaussian-input-builder.md); help `help_gaussian_input` (ends with `help_fukui_input`).
  - **Chain:** `toolLoader → gaussianInputBuilder → gaussianInput → {gaussianKeywords, fukuiInput, systemProperties, (molecularFileFormats, gaussianWavefunction: lazy)}`. Only the window imports Toga, and it must not import `cubeTab` (it has its own `_assign`); subprocess tests pin both.
  - **Tables, not code paths:** every control's label → route form is a `gaussianKeywords` table, and `_RouteReader` (Open input) reads routes back through `gk.reverse` of the same tables, so a new option is one table entry plus a reader branch. `GaussianStepSpec` holds labels, not keywords. `STEP_WIDGETS` (field → widget, kind) drives `_collect_step`/`_load_step`; a test requires it to name every spec field, and `FIELD_LABELS` likewise.
  - **Route algebra:** `RouteKeyword(name, options)`; `merge_keywords` merges Additional Keywords into a panel keyword (same option kept once, a valued option with another value is an error naming the panel). Routes are written in `fukuiInput.build_route` order (print level, model, keywords; wrapped below 80 columns). `IOp` is written `IOp(9/40=4)`.
  - **Linked steps:** `GaussianJobSpec.steps`; the geometry is job-wide. `prepare` checks the checkpoint chain (`Geom=Check/AllCheck`, `Guess=Read`, `ReadFC` need a `%Chk` of an earlier step or `%OldChk`; `TCheck` needs none). `linked_step` clears charge and multiplicity (never assumed). Every message is `- Step n: ...`.
  - **Reopening:** `read_input_text`/`read_input_file`; a file the builder wrote reads back to the same spec (test over ~20 job shapes), and every `.gjf` in `tests/data` reopens, prepares and keeps its keywords. A `%Chk` equal to `<stem>.chk` and the default title come back blank (the builder's defaults).
  - **Fukui functions job type:** `fukui_spec` maps the panels onto `fi.FukuiJobSpec` and `GaussianInputBuilder` delegates to `fi.FukuiInputBuilder`, so its files are byte-identical to the former window's (test). `fukui_defaults` resets the controls the Fukui job sets itself (`FUKUI_FIXED`; non-default values are refused) and switches `SCF=Tight NoSymm` on. The rest of the Fukui rules below are `fukuiInput`'s and unchanged; `formchk_script` is the one formchk-script writer for both job kinds (pinned by `tests/data/gaussian_input/phenol_formchk.sh`).
```

Then edit the existing Fukui bullet's first line so it reads `- **[fukuiInput.py](fukuiInput.py)** — the Fukui functions job type of the Gaussian 16 Input Builder: one Gaussian 16 (Rev. C.01, WSL or a Linux cluster) linked job ...` (remove `fukuiInputBuilder.py`, the `**Inputs > Gaussian Fukui Input Builder**` name, and in its *Split* sub-bullet replace "`fukuiInputBuilder` is only the window (`_collect_inputs` → spec, errors → dialogs). A test sets every widget ..." by "the window is `gaussianInputBuilder` (Fukui job type); `TestFukuiInTheWindow` sets every widget to a non-default value and asserts **every `FukuiJobSpec` field** comes from the window."). In its *Tests* line, replace `[tests/test_fukuiInputBuilder.py](tests/test_fukuiInputBuilder.py)` with `[tests/test_gaussianInputBuilder.py](tests/test_gaussianInputBuilder.py)`.

Also update the third CLAUDE.md reference (`grep -n fukuiInputBuilder CLAUDE.md` after the edits must print nothing).

- [ ] **Step 7: Check that nothing still names the old tool**

First update the two remaining docstrings: `fukuiFunctionAnalyze.py:842` "with the Fukui input builder's names" → "with the names the Fukui functions job type gives them", and the first line of `tests/test_fukuiInput.py` → `"""Tests for the core of the Fukui functions job type (fukuiInput, no Toga)."""`.

Run: `git grep -n -i "fukui input builder\|fukuiInputBuilder" -- . ":(exclude)docs/superpowers"`
Expected: only two lines. (1) `fukuiInput.py`'s `FUKUI_SCRIPT_HEADER` ("# gQTEA Fukui input builder: convert the three checkpoints ..."): it is text inside the written script, pinned byte-for-byte, so it stays. (2) The historical USER_MANUAL version entry ("Version 0.8.5 ... added the Gaussian Fukui Input Builder"): keep its words, but change its link anchor to `#fukui-functions-job-type`.

- [ ] **Step 7b: Version 1.7.0**

- `gqteaWinToga.py:47`: `"gQTEA-1.6.0 Molecular Analysis Toolkit",` → `"gQTEA-1.7.0 Molecular Analysis Toolkit",`
- `readme.md:10`: `gQTEA-1.6.0` → `gQTEA-1.7.0`.
- `USER_MANUAL.md:154` and `USER_MANUAL.md:7194` (the two `gQTEA-1.6.0 Molecular Analysis Toolkit` text blocks): `1.6.0` → `1.7.0`.
- `USER_MANUAL.md`, *Version, Contributors, and License*: insert this paragraph before the one that starts `Version 1.6.0 (October 2026)`:

```markdown
Version 1.7.0 (November 2026) replaces the Gaussian Fukui Input Builder on the **Inputs** tab with the [Gaussian 16 Input Builder](#gaussian-16-input-builder), modelled on GaussView 6.1's Gaussian Calculation Setup. It writes Gaussian 16 inputs for energies, optimizations (minimum or transition state), frequencies, IRC, relaxed scans, stability and NMR jobs, with ground or excited states (TD-SCF, TDA, CIS, ZINDO), Hartree-Fock, DFT, semi-empirical, MP2, coupled-cluster and compound methods, basis sets, dispersion, solvation, population analysis and guess options, and a live route line. Several steps can be linked in one input (`--Link1--`), for example an optimization followed by a TD-SCF single point, and **Open input...** reads an existing `.gjf` back into the panels. The Fukui-function calculation is now its job type [Fukui functions](#fukui-functions-job-type) and writes exactly the same files as before.
```

- `CLAUDE.md:7`: `current app version 1.6.0` → `current app version 1.7.0`.
- Leave `installation_requirements.txt:21` ("Verified working (gqteaWinToga 1.6.0) ...") and CLAUDE.md's "before 1.6.0" (runatom entry) unchanged: they record what was verified, and history.

Check: `git grep -n "1\.6\.0" -- . ":(exclude)docs/superpowers"` shows only those two lines and the `Version 1.6.0 (October 2026)` history paragraph.

- [ ] **Step 8: Commit**

```bash
git add USER_MANUAL.md readme.md docs/index.html CLAUDE.md fukuiFunctionAnalyze.py tests/test_fukuiInput.py gqteaWinToga.py
git commit -m "docs(gaussian): Gaussian 16 input builder and its Fukui job type"
```

---

### Task 13: Verification with real Gaussian 16 (WSL) and the full suite

**Files:**
- Create: `tools/check_gaussian_inputs.py` (dev-time; never run by the suite)
- Modify (only if Step 4 requires it): `gaussianInput.py`, `tests/test_gaussianInput.py`, `help.py`, `USER_MANUAL.md`, the spec
- Modify: `docs/superpowers/specs/2026-11-03-gaussian-input-builder-design.md` (status line, measured results)

G16 C.01 is installed in the user's WSL (`/usr/local/g16/g16`, checked 2026-11-03), as it was for the Fukui builder's verification.

- [ ] **Step 1: Full suite**

Run: `python -m pytest tests/ -q`
Expected: only the four baseline failures listed in Global Constraints.

- [ ] **Step 2: Write `tools/check_gaussian_inputs.py`**

```python
"""Dev-time check of the Gaussian 16 input builder against the real program (never run by the suite).

Writes one small input per job type into an empty folder:
    python tools/check_gaussian_inputs.py <folder>
then, in WSL with Gaussian 16 set up, in that folder:
    for f in *.gjf; do g16 < "$f" > "${f%.gjf}.log"; done
    grep -c "Normal termination" *.log
Each log must count one Normal termination per step (the Fukui job: 3).
"""
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import gaussianInput as gi           # noqa: E402
import gaussianKeywords as gk        # noqa: E402

WATER = ("O 0.000000 0.000000 0.117300\n"
         "H 0.000000 0.757200 -0.469200\n"
         "H 0.000000 -0.757200 -0.469200")
HCN_TS_GUESS = ("C 0.000000 0.000000 0.000000\n"
                "N 0.000000 0.000000 1.180000\n"
                "H 1.100000 0.000000 0.600000")
BASE = gi.GaussianStepSpec(family=gk.FAMILY_HF, method_name=gk.NO_CHOICE, basis="STO-3G", multiplicity=1,
                           memory="1GB", nproc=2)


def step(**changes):
    return replace(BASE, **changes)


TS = step(job_type=gk.JOB_OPT, opt_target=gk.OPT_TS, force_constants=gk.FC_CALCFC,
          additional_keywords="Opt=NoEigenTest")

JOBS = {
    "energy": (WATER, [step()], True),
    "optfreq": (WATER, [step(job_type=gk.JOB_OPT_FREQ, freq_raman=True)], False),
    "scan": (WATER, [step(job_type=gk.JOB_SCAN, modredundant="B 1 2 S 3 0.05")], False),
    "freeze": (WATER, [step(job_type=gk.JOB_OPT, modredundant="A 2 1 3 F")], False),
    "stability": (WATER, [step(job_type=gk.JOB_STABILITY, stable_opt=True)], False),
    "nmr": (WATER, [step(job_type=gk.JOB_NMR, nmr_spinspin=True)], False),
    "cation_uhf": (WATER, [step(charge=1, multiplicity=2, spin="Unrestricted (U)")], False),
    "cis": (WATER, [step(state=gk.STATE_CIS, n_states=3)], False),
    "pm6": (WATER, [step(family=gk.FAMILY_SEMI, method_name="PM6")], False),
    "cbs4m": (WATER, [step(family=gk.FAMILY_COMPOUND, method_name="CBS-4M")], False),
    "dft_smd": (WATER, [step(family=gk.FAMILY_DFT, method_name="B3LYP", basis="6-31G", basis_polarization="(d)",
                             dispersion="GD3BJ", solvation="SMD", solvent="Water", pop_mk=True, scf_tight=True,
                             nosymm=True, print_level="#N")], False),
    "linked_td": (WATER, [step(job_type=gk.JOB_OPT),
                          replace(gi.linked_step(step(job_type=gk.JOB_OPT)), state=gk.STATE_TD, n_states=3,
                                  print_all_coefficients=True)], True),
    "ts_irc": (HCN_TS_GUESS, [TS, replace(gi.linked_step(TS), job_type=gk.JOB_IRC, irc_max_points=5)], False),
    "irc_without_fc": (HCN_TS_GUESS, [TS, replace(gi.linked_step(TS), job_type=gk.JOB_IRC,
                                                  irc_force_constants=gk.FC_DEFAULT, irc_max_points=5)], False),
    "fukui": (WATER, [gi.fukui_defaults(step(fukui_plus_multiplicity=2, fukui_minus_multiplicity=2))], True),
}


def main(folder: str) -> None:
    out = Path(folder)
    out.mkdir(parents=True, exist_ok=True)
    builder = gi.GaussianInputBuilder()
    for name, (geometry, steps, formchk) in JOBS.items():
        job = builder.prepare(gi.GaussianJobSpec(stem=name, output_dir=str(out), geometry=geometry, steps=steps,
                                                 write_formchk=formchk))
        builder.write_files(job)
        print(f"{job.input_path.name}: {len(job.steps) or 3} step(s)"
              + "".join(f"\n  WARNING {w}" for w in job.warnings))


if __name__ == "__main__":
    main(sys.argv[1])
```

- [ ] **Step 3: Write the inputs and run them in WSL**

```powershell
python tools/check_gaussian_inputs.py "$env:TEMP\gaussian_input_check"
```
Expected: 15 files listed; `irc_without_fc` and `linked_td`/`fukui` print their warnings (IRC without force constants; Density=Current; shared checkpoint).

Then run them in WSL, in the same folder as WSL sees it (`gaussianInput.wsl_path` converts it):

```powershell
$wslDir = python -c "import os, gaussianInput as gi; print(gi.wsl_path(os.path.join(os.environ['TEMP'], 'gaussian_input_check')))"
wsl -e bash -lc "cd '$wslDir' && for f in *.gjf; do g16 < `"`$f`" > `"`${f%.gjf}.log`"; done; grep -c 'Normal termination' *.log; bash energy_formchk.sh && bash linked_td_formchk.sh && bash fukui_formchk.sh && ls *.fchk"
```

(If the PowerShell quoting gets in the way, write the same three bash lines to `run.sh` in that folder with LF endings and run `wsl -e bash -lc "cd '$wslDir' && bash run.sh"`.)

Expected: `Normal termination` counts — 1 for every single-step job, 2 for `linked_td`, `ts_irc` and (if Gaussian accepts it) `irc_without_fc`, 3 for `fukui`; the three scripts make `energy.fchk`, `linked_td.fchk`, `fukui_N.fchk`, `fukui_Np1.fchk`, `fukui_Nm1.fchk`.

- [ ] **Step 4: Act on what Gaussian did**

- A job that fails: read the end of its `.log`. A refused keyword or input layout is a builder bug: fix `gaussianInput`, add a test that pins the corrected text, rerun Steps 1–3. Do not change the job list to hide it.
- `irc_without_fc`: if G16 refuses the IRC without force constants, make it an **error** in `_prepare_step` (replace the `warnings.append("IRC without force constants ...")` with an `errors.append(...)` naming the G16 message), update `test_ts_and_irc_without_force_constants_warn` (IRC part → `pytest.raises`), `tools/check_gaussian_inputs.py` (drop that job) and the help/manual sentences on IRC. If it runs, keep the warning. Either way, record the observed behaviour in the spec §5 "Warnings" line.
- Open `linked_td.gjf` in the window (Open input...) and press Preview: the text must equal the file (spec section 6 round trip, checked here on a file Gaussian has run).

- [ ] **Step 5: Record and commit**

In the spec, change the status line to `Status: implemented <date>; verified with Gaussian 16 C.01 in WSL (tools/check_gaussian_inputs.py: <n> jobs, all Normal termination)` and add the IRC finding under §5.

```bash
git add tools/check_gaussian_inputs.py docs/superpowers/specs/2026-11-03-gaussian-input-builder-design.md
git commit -m "test(gaussian): every job type checked with Gaussian 16 in WSL"
```

(Include any files changed in Step 4 in this commit.)

