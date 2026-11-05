# Optimized geometries from outputs + Opt/excited-states job type: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The Gaussian 16 Input Builder loads the final converged geometry from a Gaussian `.log`/`.out`, and gains a "Geometry optimization + excited states" job type whose two `--Link1--` steps feed the Excitation/NTO tabs.

**Architecture:** A new Toga-free `gaussianOutput.py` parses Gaussian outputs into jobs and returns the last completed optimization. The composite job type follows the Fukui pattern: one `GaussianStepSpec` with extra `excited_*` fields is expanded by a pure `opt_excited_steps` into two ordinary steps, which go through the existing `_prepare_step`. A new Method-panel basis choice "Read from checkpoint (ChkBasis)" lets step 2 reuse a Gen/GenECP basis and reopen correctly.

**Tech Stack:** Python 3.13, Toga 0.5 (window only), pytest; Gaussian 16 C.01 in WSL for fixtures (dev-time only).

**Spec:** [docs/superpowers/specs/2026-11-05-gaussian-opt-excited-design.md](../specs/2026-11-05-gaussian-opt-excited-design.md)

## Global Constraints

- Only `gaussianInputBuilder.py` imports Toga; `gaussianOutput`, `gaussianInput` and `gaussianKeywords` stay Toga-free (subprocess tests).
- The builder window must not import `cubeTab`.
- Files written are ASCII, LF, with relative checkpoint names; target G16 C.01 (`%OldChk` is required).
- Multiplicity is never assumed.
- Widget writes go through `_assign`; handlers are synchronous and `_syncing`-guarded.
- Every `GaussianStepSpec` field has a `FIELD_LABELS` entry and a `STEP_WIDGETS` entry (existing tests).
- Every `.gjf` under `tests/data` must reopen, prepare and read back to the same spec (`tests/test_gaussianInputReopen.py`).
- Baseline: 4 tests fail on `main` already (CLAUDE.md); compare against that.
- Run tests from `venv/src` with `..\Scripts\python.exe -m pytest`.

## Files

| File | Change |
|------|--------|
| `gaussianOutput.py` | **new**: Gaussian output → `OptimizedGeometry` |
| `tests/test_gaussianOutput.py` | **new** |
| `tests/data/gaussian_output/` | **new** fixtures, `inputs/*.gjf`, `run_gaussian.py`, `README.md` |
| `gaussianKeywords.py` | `BASIS_CHECKPOINT`, `JOB_OPT_EXCITED`, excited-state tables |
| `gaussianInput.py` | ChkBasis; `GeometryFile`/`load_geometry_file`; spec fields; `opt_excited_*`; `_prepare_steps` refactor; previews; reopen note |
| `gaussianInputBuilder.py` | new widgets/section, rules, job-type choice, geometry notes dialog, run instructions |
| `excitationData.py` | `sibling_excitation_file` tries the stem without `_excited` |
| `tests/test_gaussianInput.py`, `tests/test_gaussianInputBuilder.py`, `tests/test_gaussianInputReopen.py`, `tests/test_excitationData.py` | new cases |
| `tools/check_gaussian_inputs.py` | `opt_excited`, `opt_excited_gen` jobs |
| `help.py`, `USER_MANUAL.md`, `CLAUDE.md` | docs |

---

### Task 1: Fixtures and the Gaussian output reader

**Files:**
- Create: `gaussianOutput.py`, `tests/test_gaussianOutput.py`, `tests/data/gaussian_output/{README.md, run_gaussian.py}`
- Move: `tests/data/gaussian_output/raw/*` → `tests/data/gaussian_output/` (`h2o_opt_freq.log` → `h2o_opt_freq.out`; `h2co_opt_td_excited.fchk` gzip → `.fchk.gz`); delete `raw/` and the `h2co_opt_td_gen.*` outputs (only its input is kept)

**Interfaces:**
- Produces: `gaussianOutput.read_optimized_geometry(path) -> OptimizedGeometry`; `parse_optimized_geometry(text, name) -> OptimizedGeometry`; `OptimizedGeometry(atoms: List[Tuple[str, float, float, float]], charge: Optional[int], multiplicity: Optional[int], job_number: int, job_count: int, orientation: str, notes: List[str])`.

Markers verified in the G16 C.01 logs (2026-11-05): a job ends with `Normal termination of Gaussian`; the next starts with `Initial command:` (`--Link1--`) or `Link1:  Proceeding to internal job step number` (Opt Freq). **The Freq internal job also prints `Optimization completed.`**, so only jobs whose route has an `Opt` keyword count. Failed: `Optimization stopped.` + `-- Number of steps exceeded,  NStep=   1`, then `Error termination`. Orientation blocks: title, dashes, 2 heading lines, dashes, rows `center Z type x y z`, dashes.

- [ ] **Step 1: Organise the fixtures** (shell moves above), write `README.md` (table: file, route, what it covers) and `run_gaussian.py` (copies `inputs/*.gjf` to `/tmp/gqout` through WSL, runs `g16 < x.gjf > x.log`, formchks `h2co_opt_td_excited.chk`, copies back, renames `h2o_opt_freq.log` → `.out`, gzips the fchk; never run by the suite).

- [ ] **Step 2: Write the failing tests** `tests/test_gaussianOutput.py`:

```python
"""Optimized geometries from Gaussian 16 outputs (gaussianOutput, no Toga)."""
import math
import subprocess
import sys
from pathlib import Path

import pytest

import gaussianOutput as go

SRC = Path(__file__).resolve().parents[1]
DATA = SRC / "tests" / "data" / "gaussian_output"


def test_module_does_not_import_toga():
    out = subprocess.run([sys.executable, "-c", "import sys, gaussianOutput; print('toga' in sys.modules)"],
                         capture_output=True, text=True, cwd=str(SRC))
    assert out.stdout.strip() == "False", out.stderr


def angle(a, b, c):
    u = [a[i] - b[i] for i in range(3)]
    v = [c[i] - b[i] for i in range(3)]
    dot = sum(x * y for x, y in zip(u, v))
    return math.degrees(math.acos(dot / (math.hypot(*u) * math.hypot(*v))))


class TestFixtures:
    def test_converged_optimization(self):
        g = go.read_optimized_geometry(DATA / "h2o_opt.log")
        assert [a[0] for a in g.atoms] == ["O", "H", "H"]
        assert (g.charge, g.multiplicity, g.job_number, g.job_count) == (0, 1, 1, 1)
        assert g.orientation == "Input orientation" and g.notes == []
        oh = math.dist(g.atoms[0][1:], g.atoms[1][1:])
        assert 0.95 < oh < 0.98          # HF/3-21G water: 0.967 A (the start was 1.00 A)

    def test_the_converged_geometry_is_the_last_block_before_completion(self):
        text = (DATA / "h2o_opt.log").read_text()
        g = go.parse_optimized_geometry(text, "h2o_opt.log")
        lines = text.splitlines()
        done = next(k for k, l in enumerate(lines) if "Optimization completed" in l)
        start = max(k for k in range(done) if "Input orientation:" in lines[k])
        x, y, z = (float(v) for v in lines[start + 5].split()[3:6])
        assert g.atoms[0][1:] == (x, y, z)

    def test_opt_freq_out_ignores_the_freq_job(self):
        g = go.read_optimized_geometry(DATA / "h2o_opt_freq.out")
        assert (g.job_number, g.job_count) == (1, 2) and g.notes == []

    def test_a_later_failed_optimization_is_reported(self):
        g = go.read_optimized_geometry(DATA / "h2o_opt_then_failed.log")
        assert (g.job_number, g.job_count, g.charge, g.multiplicity) == (1, 2, 0, 1)
        assert len(g.notes) == 1 and "job 2" in g.notes[0] and "Number of steps exceeded" in g.notes[0]

    def test_maxcycles_is_refused(self):
        with pytest.raises(ValueError, match="no successfully optimized geometry.*Number of steps exceeded"):
            go.read_optimized_geometry(DATA / "h2co_opt_maxcycles.log")

    def test_opt_then_td_gives_the_optimization(self):
        g = go.read_optimized_geometry(DATA / "h2co_opt_td.log")
        assert [a[0] for a in g.atoms] == ["C", "O", "H", "H"] and (g.job_number, g.job_count) == (1, 2)


class TestSynthetic:
    def test_cut_off_before_completion(self):
        text = (DATA / "h2o_opt.log").read_text()
        cut = text[:text.index("Optimization completed")]
        with pytest.raises(ValueError, match="ends before the optimization completed"):
            go.parse_optimized_geometry(cut, "cut.log")

    def test_completed_but_no_normal_termination_is_a_note(self):
        text = (DATA / "h2o_opt.log").read_text()
        g = go.parse_optimized_geometry(text[:text.index("Normal termination")], "x.log")
        assert any("Normal termination" in n for n in g.notes)

    def test_not_a_gaussian_output(self):
        with pytest.raises(ValueError, match="not a Gaussian output"):
            go.parse_optimized_geometry("* O   R   A   *\nORCA TERMINATED NORMALLY\n", "orca.out")

    def test_no_optimization(self):
        text = (DATA / "h2co_opt_td.log").read_text()
        second = text[text.index("Normal termination") :]
        with pytest.raises(ValueError, match="contains no geometry optimization"):
            go.parse_optimized_geometry(" Entering Gaussian System\n" + second.split("\n", 1)[1], "td.log")

    def test_ghost_atoms_are_refused(self):
        text = (DATA / "h2o_opt.log").read_text()
        lines = text.splitlines()
        done = next(k for k, l in enumerate(lines) if "Optimization completed" in l)
        start = max(k for k in range(done) if "Input orientation:" in lines[k])
        row = lines[start + 6].split()
        lines[start + 6] = f"      {row[0]}          0           0     {row[3]}    {row[4]}    {row[5]}"
        with pytest.raises(ValueError, match="ghost or dummy"):
            go.parse_optimized_geometry("\n".join(lines), "x.log")

    def test_relaxed_scan_note(self):
        text = (DATA / "h2o_opt.log").read_text()
        doubled = text.replace(" Optimization completed.", " Optimization completed.\n Optimization completed.", 1)
        g = go.parse_optimized_geometry(doubled, "scan.log")
        assert any("relaxed scan" in n for n in g.notes)
```

- [ ] **Step 3: Run them — expect `ModuleNotFoundError: gaussianOutput`**

Run: `..\Scripts\python.exe -m pytest tests/test_gaussianOutput.py -q`

- [ ] **Step 4: Implement `gaussianOutput.py`**

```python
"""
The final optimized geometry of a Gaussian output (.log/.out), for the Gaussian 16 Input Builder (no Toga).

A Gaussian output is a sequence of jobs (each --Link1-- step, and the internal Freq job of Opt Freq),
each ending with "Normal termination". A job optimizes when its route has an Opt keyword (the Freq job
of Opt Freq also prints "Optimization completed", so the route decides). The geometry imported is the
last orientation block before the last "Optimization completed" of the last optimizing job that has
one; an optimization that never completed is refused, never imported.
"""
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

from fukuiInput import element_symbol

GAUSSIAN_BANNERS = ("Entering Gaussian System", "Gaussian, Inc.")
NORMAL_TERMINATION = "Normal termination of Gaussian"
ERROR_TERMINATION = "Error termination"
COMPLETED = "Optimization completed"
STOPPED = "Optimization stopped"
ORIENTATIONS = ("Input orientation", "Standard orientation", "Z-Matrix orientation")   # preference order
CHARGE_RE = re.compile(r"^\s*Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
OPT_RE = re.compile(r"(?i)(?:^|[\s/])opt(?=$|[\s=(])")

Atom = Tuple[str, float, float, float]


@dataclass
class OptimizedGeometry:
    atoms: List[Atom]
    charge: Optional[int]
    multiplicity: Optional[int]
    job_number: int
    job_count: int
    orientation: str
    notes: List[str] = field(default_factory=list)


@dataclass
class _Job:
    number: int
    lines: List[str]
    route: str
    completions: List[int]
    normal: bool

    @property
    def optimizes(self) -> bool:
        return bool(OPT_RE.search(self.route.lstrip("#").split(None, 1)[-1] if self.route else ""))

    def failure(self) -> str:
        for k, line in enumerate(self.lines):
            if STOPPED in line:
                detail = self.lines[k + 1].strip().lstrip("-").strip() if k + 1 < len(self.lines) else ""
                return f"the optimization stopped ({detail})" if detail else "the optimization stopped"
        if any(ERROR_TERMINATION in line for line in self.lines):
            return "Gaussian ended with an error termination during the optimization"
        return "the file ends before the optimization completed (the job is still running or was cut off)"


def _route(lines: List[str]) -> str:
    start = next((k for k, line in enumerate(lines) if line.strip().startswith("#")), None)
    if start is None:
        return ""
    words = []
    for line in lines[start:]:
        if line.strip().startswith("---"):
            break
        words.append(line.strip())
    return " ".join(words)


def _jobs(lines: List[str]) -> List[_Job]:
    jobs, start = [], 0
    bounds = []
    for k, line in enumerate(lines):
        if NORMAL_TERMINATION in line:
            bounds.append((start, k + 1, True))
            start = k + 1
    if any(line.strip() for line in lines[start:]):
        bounds.append((start, len(lines), False))
    for number, (a, b, normal) in enumerate(bounds, start=1):
        part = lines[a:b]
        jobs.append(_Job(number, part, _route(part), [k for k, l in enumerate(part) if COMPLETED in l], normal))
    return jobs


def _block(lines: List[str], start: int, name: str) -> List[Atom]:
    atoms, k = [], start + 5          # title, dashes, two heading lines, dashes
    while k < len(lines) and not lines[k].strip().startswith("---"):
        parts = lines[k].split()
        try:
            z, xyz = int(parts[1]), tuple(float(v) for v in parts[-3:])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{name}, line {k + 1}: cannot read the orientation row '{lines[k].strip()}'.") from exc
        if z <= 0:
            raise ValueError(f"{name}: atom {parts[0]} of the optimized geometry is a ghost or dummy atom "
                             f"(atomic number {z}); the builder takes real atoms only.")
        symbol = element_symbol(str(z))
        if symbol is None:
            raise ValueError(f"{name}: atomic number {z} is not an element.")
        atoms.append((symbol, *xyz))
        k += 1
    if not atoms:
        raise ValueError(f"{name}, line {start + 1}: the orientation block holds no atoms.")
    return atoms


def _converged_block(job: _Job, name: str) -> Tuple[List[Atom], str]:
    """The orientation block of the geometry step that the last completion belongs to (Input preferred)."""
    lines, done = job.lines, job.completions[-1]
    titles = {k: t for k, line in enumerate(lines[:done]) for t in ORIENTATIONS if f"{t}:" in line}
    if not titles:
        raise ValueError(f"{name}: job {job.number} completed its optimization but printed no geometry.")
    last = max(titles)
    earlier = [k for k in titles if k < last and titles[k] == titles[last]]
    step = {t: k for k, t in titles.items() if k > (max(earlier) if earlier else -1)}
    title = next(t for t in ORIENTATIONS if t in step)
    return _block(lines, step[title], name), title


def parse_optimized_geometry(text: str, name: str) -> OptimizedGeometry:
    if not any(banner in text for banner in GAUSSIAN_BANNERS):
        raise ValueError(f"{name} is not a Gaussian output (no Gaussian banner); load a Gaussian .log/.out, "
                         "an input file or an .xyz.")
    jobs = _jobs(text.splitlines())
    optimizing = [j for j in jobs if j.optimizes]
    if not optimizing:
        raise ValueError(f"{name} contains no geometry optimization (no route has Opt); load its input file "
                         "or an .xyz instead.")
    done = [j for j in optimizing if j.completions]
    if not done:
        last = optimizing[-1]
        raise ValueError(f"{name}: no successfully optimized geometry was found: in job {last.number}, "
                         f"{last.failure()}. The last geometry is not imported because it is not optimized.")
    job = done[-1]
    atoms, orientation = _converged_block(job, name)
    charges = [m for m in (CHARGE_RE.match(l) for l in job.lines[:job.completions[-1]]) if m]
    charge, mult = (int(charges[-1].group(1)), int(charges[-1].group(2))) if charges else (None, None)
    notes = [f"job {j.number}'s optimization did not complete ({j.failure()}); the geometry of job {job.number}'s "
             "completed optimization was imported." for j in optimizing if j.number > job.number]
    if len(job.completions) > 1:
        notes.append(f"job {job.number} completed {len(job.completions)} optimizations (a relaxed scan); the "
                     "last one's geometry was imported.")
    if not job.normal:
        notes.append(f"job {job.number} did not end with Normal termination; its optimization had completed, "
                     "so that geometry was imported.")
    return OptimizedGeometry(atoms, charge, mult, job.number, len(jobs), orientation, notes)


def read_optimized_geometry(path) -> OptimizedGeometry:
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ValueError(f"'{path.name}' could not be opened: {exc.strerror or exc}") from exc
    return parse_optimized_geometry(text, path.name)
```

- [ ] **Step 5: Run the tests — expect all PASS**

- [ ] **Step 6: Commit** `feat(gaussian): read the optimized geometry of a Gaussian output (+ G16 fixtures)`

---

### Task 2: "Load geometry..." reads `.log`/`.out`

**Files:** Modify `gaussianInput.py` (end section), `gaussianInputBuilder.py` (`_build_molecule_panel`, `load_geometry`, `apply_geometry_file`, `_show_info`); Test: `tests/test_gaussianInput.py::TestGeometryFiles`, `tests/test_gaussianInputBuilder.py::TestGeometryAndFiles`

**Interfaces:**
- Consumes: `gaussianOutput.read_optimized_geometry`.
- Produces: `gi.GeometryFile(text, charge, multiplicity, notes, source)`; `gi.load_geometry_file(path) -> GeometryFile`; `gi.OUTPUT_SUFFIXES = (".log", ".out")`; `GaussianInputBuilderUI.apply_geometry_file(path) -> List[str]` (the notes); `GaussianInputBuilderUI._show_info(title, message)` (async).

- [ ] **Step 1: Failing tests**

```python
# tests/test_gaussianInput.py, class TestGeometryFiles
    OUTPUTS = SRC / "tests" / "data" / "gaussian_output"

    @pytest.mark.parametrize("name", ["h2o_opt.log", "h2o_opt_freq.out"])
    def test_gaussian_output_gives_the_optimized_geometry(self, name):
        g = gi.load_geometry_file(self.OUTPUTS / name)
        assert (g.charge, g.multiplicity) == (0, 1) and len(fi.parse_geometry(g.text)) == 3
        assert g.source.startswith("the optimized geometry of job 1 of")
        assert gi.read_geometry_file(self.OUTPUTS / name) == (g.text, 0, 1)

    def test_gaussian_output_notes(self):
        assert len(gi.load_geometry_file(self.OUTPUTS / "h2o_opt_then_failed.log").notes) == 1

    def test_unconverged_output_is_refused(self):
        with pytest.raises(ValueError, match="no successfully optimized geometry"):
            gi.load_geometry_file(self.OUTPUTS / "h2co_opt_maxcycles.log")

    def test_other_files_have_no_notes(self):
        g = gi.load_geometry_file(FUKUI_INPUTS / "phenol_Np1.gjf")
        assert g.notes == [] and g.source == "" and (g.charge, g.multiplicity) == (-1, 2)

# tests/test_gaussianInputBuilder.py, class TestGeometryAndFiles
    def test_apply_an_output_fills_geometry_charge_and_multiplicity(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.charge_input.value, ui.multiplicity_input.value = "1", "2"
        notes = ui.apply_geometry_file(OUTPUTS / "h2o_opt.log")
        assert notes == [] and ui.charge_input.value == "0" and ui.multiplicity_input.value == "1"
        assert ui.status_label.text.startswith("The optimized geometry of job 1 of 1 loaded from h2o_opt.log")

    def test_load_geometry_shows_the_notes(self, tmp_path, monkeypatch):
        ui = make_ui(tmp_path)

        async def ask(*a, **k):
            return OUTPUTS / "h2o_opt_then_failed.log"
        ui.window = SimpleNamespace(dialog=ask)
        asyncio.run(ui.load_geometry(None))
        assert ui.dialogs[0][0] == "Geometry Notes" and "job 2" in ui.dialogs[0][1]

    def test_load_geometry_refuses_an_unconverged_output(self, tmp_path):
        ui = make_ui(tmp_path)
        before = ui.geometry_input.value

        async def ask(*a, **k):
            return OUTPUTS / "h2co_opt_maxcycles.log"
        ui.window = SimpleNamespace(dialog=ask)
        asyncio.run(ui.load_geometry(None))
        assert ui.dialogs[0][0] == "Geometry Error" and ui.geometry_input.value == before
```
(add `OUTPUTS = SRC / "tests" / "data" / "gaussian_output"` at module level of the window test, and in `make_ui` stub `ui._show_info = record`.)

- [ ] **Step 2: Run — expect AttributeError `load_geometry_file`**

- [ ] **Step 3: Implement.** In `gaussianInput.py`:

```python
GEOMETRY_FILE_TYPES = ["xyz", "gjf", "com", "gau", "fchk", "fch", "log", "out", "pdb", "ent", "mol", "sdf", "mol2"]
OUTPUT_SUFFIXES = (".log", ".out")


@dataclass
class GeometryFile:
    """A loaded geometry: 'El x y z' text, charge/multiplicity when the file gives them, notes, and what it is."""
    text: str
    charge: Optional[int]
    multiplicity: Optional[int]
    notes: List[str] = field(default_factory=list)
    source: str = ""


def load_geometry_file(path) -> GeometryFile:
    """Any supported structure file; a Gaussian .log/.out gives its last completed optimization."""
    path = Path(path)
    if path.suffix.lower() in OUTPUT_SUFFIXES:
        import gaussianOutput as go
        g = go.read_optimized_geometry(path)
        return GeometryFile("\n".join(fi.format_geometry(g.atoms)), g.charge, g.multiplicity, list(g.notes),
                            f"the optimized geometry of job {g.job_number} of {g.job_count}")
    ...  # the former read_geometry_file body, returning GeometryFile(text, charge, mult)


def read_geometry_file(path) -> Tuple[str, Optional[int], Optional[int]]:
    g = load_geometry_file(path)
    return g.text, g.charge, g.multiplicity
```

In the window: hint `".xyz .gjf .com .fchk .log .out .mol .sdf .mol2 .pdb; a .gjf/.fchk/.log/.out also sets the charge and multiplicity (a .log/.out: its last completed optimization)"`; `apply_geometry_file` calls `gi.load_geometry_file`, uses `geometry.text/charge/multiplicity`, status `f"{source[:1].upper()}{source[1:]} loaded from {path.name}."` when `source` (else the old text) + the existing step note, writes `NOTE:` lines into `messages_output` when there are notes, and returns `geometry.notes`. `load_geometry`: `notes = self.apply_geometry_file(path)`; `if notes: await self._show_info("Geometry Notes", "\n".join(notes))`. Add `async def _show_info(self, title, message): await self.window.dialog(toga.InfoDialog(title, message))`.

- [ ] **Step 4: Run** `tests/test_gaussianInput.py tests/test_gaussianInputBuilder.py` — PASS (the old `test_unsupported_file` still passes: an empty `.log` is "not a Gaussian output").

- [ ] **Step 5: Commit** `feat(gaussian): Load geometry reads the optimized geometry of a .log/.out`

---

### Task 3: Basis set "Read from checkpoint (ChkBasis)"

**Files:** Modify `gaussianKeywords.py`, `gaussianInput.py` (`basis_keyword`, `method_keywords`, `_checkpoint_readers`, `_RouteReader`); Test `tests/test_gaussianInput.py`, `tests/test_gaussianInputReopen.py`

**Interfaces:** Produces `gk.BASIS_CHECKPOINT = "Read from checkpoint (ChkBasis)"` (in `gk.BASES` before `OTHER`).

- [ ] **Step 1: Failing tests**

```python
# tests/test_gaussianInput.py
class TestChkBasis:
    def test_route(self):
        assert route(basis=gk.BASIS_CHECKPOINT) == "#P B3LYP ChkBasis"

    def test_it_reads_a_checkpoint(self, tmp_path):
        assert "ChkBasis reads a checkpoint" in problems(tmp_path, step(basis=gk.BASIS_CHECKPOINT, multiplicity=1))
        gi.GaussianInputBuilder().prepare(job_spec(tmp_path, step(basis=gk.BASIS_CHECKPOINT, multiplicity=1,
                                                                    old_chk="old.chk")))

# tests/test_gaussianInputReopen.py, class TestReading
    def test_chkbasis_reads_back(self):
        spec, notes = gi.read_input_text("%OldChk=a.chk\n#P B3LYP ChkBasis Geom=AllCheck Guess=Read TD\n\n")
        assert spec.steps[0].basis == gk.BASIS_CHECKPOINT and notes == [no for no in notes if "basis" not in no]
        assert "ChkBasis" not in spec.steps[0].additional_keywords
```

- [ ] **Step 2: Run — expect AttributeError `BASIS_CHECKPOINT`**

- [ ] **Step 3: Implement**
  - `gk`: `BASIS_CHECKPOINT = ...`; `BASES = POPLE_BASES + OTHER_BASES + (BASIS_CHECKPOINT, OTHER)`; `PANEL_OF_KEYWORD["CHKBASIS"] = "Method"`.
  - `basis_keyword`: `if step.basis == gk.BASIS_CHECKPOINT: return ""` (after the `no_basis` test).
  - `method_keywords`: first item `RouteKeyword("ChkBasis")` when `step.basis == gk.BASIS_CHECKPOINT and not no_basis(step)`.
  - `_checkpoint_readers`: `elif keyword.key == "CHKBASIS": readers.append("ChkBasis")`.
  - `_RouteReader.read`: before the loop, `self.chk_basis = any(keyword_text_name(t).upper() == "CHKBASIS" for t in tokens)`; in `_model`, the no-slash branch becomes `elif self.chk_basis and not no_basis(self.step): self.step.basis = gk.BASIS_CHECKPOINT` before the STO-3G default; `handlers["CHKBASIS"] = self._switch_basis` returning `[]` (consumed) for an option-less token, `None` otherwise.

- [ ] **Step 4: Run** `tests/test_gaussianInput.py tests/test_gaussianInputReopen.py tests/test_gaussianInputBuilder.py` — PASS; `h2co_opt_td_gen.gjf` reopens with `basis == BASIS_CHECKPOINT` in step 2.

- [ ] **Step 5: Commit** `feat(gaussian): basis set Read from checkpoint (ChkBasis)`

---

### Task 4: The composite job type in `gaussianInput`

**Files:** Modify `gaussianKeywords.py`, `gaussianInput.py`; Test `tests/test_gaussianInput.py` (new class `TestOptExcited`)

**Interfaces:**
- Consumes: `gk.BASIS_CHECKPOINT` (Task 3).
- Produces (gk): `JOB_OPT_EXCITED = "Geometry optimization + excited states"` (in `JOB_TYPES` before `JOB_FUKUI`), `SAME_AS_OPTIMIZATION = "Same as optimization"`, `EXCITED_STATES = (STATE_TD, STATE_TDA, STATE_CIS)`, `EXCITED_FAMILIES = (SAME_AS_OPTIMIZATION, FAMILY_HF, FAMILY_DFT)`, `EXCITED_DFT_METHODS = tuple(FAMILY_METHODS[FAMILY_DFT])`, `EXCITED_BASES = (SAME_AS_OPTIMIZATION,) + tuple(b for b in BASES if b not in GENERAL_BASES + (BASIS_CHECKPOINT,))`, `DEFAULT_EXCITED_N_STATES = 10`, `EXCITED_CHK_SUFFIX = "_excited"`.
- Produces (gi): spec fields `opt_confirm_freq, excited_state, excited_n_states, excited_root, excited_state_spin, excited_family, excited_method_name, excited_method_other, excited_basis, excited_basis_diffuse, excited_basis_polarization, excited_basis_other, excited_coefficients, excited_keywords, excited_chk` (defaults in the spec §2.1); `opt_excited_defaults(step) -> GaussianStepSpec`; `opt_excited_steps(step, stem) -> Tuple[GaussianStepSpec, GaussianStepSpec]`; `optimization_chk(step, stem) -> str`; `excited_chk(step, stem) -> str`; `GaussianJob.excited_fchk: Optional[str]`; `OPT_EXCITED_STEP_LABELS = ("Optimization step", "Excited-state step")`.

- [ ] **Step 1: Failing tests** (`ACETALDEHYDE` not needed: water)

```python
def opt_excited(**changes):
    return gi.opt_excited_defaults(step(multiplicity=1, **changes))


class TestOptExcited:
    def test_two_linked_steps(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(job_spec(tmp_path, opt_excited()))
        assert builder.input_text(job) == "\n".join(
            ["%Chk=water.chk", "%Mem=2GB", "%NProcShared=4", "#P B3LYP/6-31G(d) Opt", "",
             "water", "", "0 1", *WATER_LINES, "",
             "--Link1--",
             "%OldChk=water.chk", "%Chk=water_excited.chk", "%Mem=2GB", "%NProcShared=4",
             "#P B3LYP/6-31G(d) TD=(NStates=10) IOp(9/40=4) Geom=AllCheck Guess=Read", ""]) + "\n"
        assert job.chk_files == ["water.chk", "water_excited.chk"] and job.excited_fchk == "water_excited.fchk"
        assert job.warnings == []

    def test_opt_freq_switch(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, opt_excited(opt_confirm_freq=True, freq_raman=True)))
        assert job.steps[0].route == ["#P B3LYP/6-31G(d) Opt Freq=Raman"]

    def test_own_method_and_basis_give_a_new_guess(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, opt_excited(
            excited_family=gk.FAMILY_DFT, excited_method_name="CAM-B3LYP", excited_basis="6-311G",
            excited_basis_diffuse="+", excited_basis_polarization="(d,p)", excited_n_states=6,
            excited_state_spin="Singlets", dispersion="GD3BJ")))
        assert job.steps[0].route == ["#P B3LYP/6-31G(d) Opt EmpiricalDispersion=GD3BJ"]
        assert job.steps[1].route == ["#P CAM-B3LYP/6-311+G(d,p) TD=(NStates=6,Singlets) IOp(9/40=4) Geom=AllCheck"]
        assert any("new guess" in w for w in job.warnings)

    def test_same_basis_with_another_functional_reads_the_guess(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, opt_excited(
            excited_family=gk.FAMILY_DFT, excited_method_name="CAM-B3LYP", excited_state=gk.STATE_TDA)))
        assert job.steps[1].route == ["#P CAM-B3LYP/6-31G(d) TDA=(NStates=10) IOp(9/40=4) Geom=AllCheck Guess=Read"]

    def test_gen_basis_is_read_from_the_checkpoint(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, opt_excited(
            basis="Gen", additional_input="O H 0\n6-31G(d)\n****")))
        assert job.steps[1].route == ["#P B3LYP ChkBasis TD=(NStates=10) IOp(9/40=4) Geom=AllCheck Guess=Read"]
        assert job.steps[1].additional_input == []

    def test_cis(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, opt_excited(excited_state=gk.STATE_CIS,
                                                                                  excited_n_states=4)))
        assert job.steps[1].route == ["#P CIS(NStates=4)/6-31G(d) IOp(9/40=4) Geom=AllCheck Guess=Read"]

    def test_shared_and_step_two_keywords(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, opt_excited(
            scf_tight=True, nosymm=True, solvation="SMD", solvent="Acetone", pop_nbo=True,
            additional_keywords="Opt=NoEigenTest Int=UltraFine", excited_keywords="Density=Current")))
        assert job.steps[0].route == ["#P B3LYP/6-31G(d) Opt=NoEigenTest SCF=Tight NoSymm Pop=NBO",
                                      "SCRF=(SMD,Solvent=Acetone) Int=UltraFine"]
        assert " ".join(job.steps[1].route) == ("#P B3LYP/6-31G(d) TD=(NStates=10) IOp(9/40=4) SCF=Tight NoSymm "
                                                "Geom=AllCheck Guess=Read SCRF=(SMD,Solvent=Acetone) Int=UltraFine "
                                                "Density=Current")

    def test_open_shell_reference(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, opt_excited(charge=1, multiplicity=2)))
        assert job.steps[1].charge_line is None and "Geom=AllCheck" in job.steps[1].route[0]

    @pytest.mark.parametrize("changes,fragment", [
        (dict(family=gk.FAMILY_COMPOUND, method_name="CBS-QB3"), "Optimization step: Compound methods"),
        (dict(spin="Restricted-open (RO)", charge=1, multiplicity=2), "Restricted-open (RO)"),
        (dict(family=gk.FAMILY_SEMI, method_name="PM7"), "Excited-state step: the optimization (PM7) has no basis"),
        (dict(family=gk.FAMILY_MP2, method_name=gk.NO_CHOICE), "Excited-state step: TD-SCF needs"),
        (dict(excited_state=gk.STATE_TDA, excited_family=gk.FAMILY_HF), "Excited-state step: TDA needs a DFT"),
        (dict(charge=1, multiplicity=2, excited_state_spin="Triplets"), "Excited-state step: States = Triplets"),
        (dict(excited_chk="WATER.chk"), "Excited-state %Chk 'WATER.chk' is the optimization's %Chk"),
        (dict(excited_chk="a b.chk"), "Excited-state step: %Chk 'a b.chk'"),
        (dict(excited_keywords="Freq"), "Excited-state keywords: Freq"),
        (dict(excited_root=11), "Excited-state step: State of interest (Root) must be between 1 and 10"),
        (dict(modredundant="B 1 2 S 3 0.05"), "Optimization step: Scan lines (S)"),
        (dict(multiplicity=None), "Optimization step: Multiplicity: enter it"),
    ])
    def test_problems(self, tmp_path, changes, fragment):
        assert fragment in problems(tmp_path, opt_excited(**changes))

    def test_semi_empirical_optimization_with_an_excited_basis(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, opt_excited(
            family=gk.FAMILY_SEMI, method_name="PM7", excited_family=gk.FAMILY_DFT, excited_basis="6-31G")))
        assert job.steps[0].route == ["#P PM7 Opt"]
        assert job.steps[1].route == ["#P B3LYP/6-31G(d) TD=(NStates=10) IOp(9/40=4) Geom=AllCheck"]
        assert any("PM7" in w and "new guess" in w for w in job.warnings)

    def test_must_be_the_only_step(self, tmp_path):
        assert "must be the only step" in problems(tmp_path, opt_excited(), step(multiplicity=1))

    @pytest.mark.parametrize("changes,label", [(dict(opt_target=gk.OPT_TS), "Optimize to a"),
                                               (dict(state=gk.STATE_TD), "Electronic state"),
                                               (dict(geometry_source=gk.GEOM_CHECK), "Geometry source"),
                                               (dict(old_chk="x.chk"), "%OldChk")])
    def test_controls_the_job_sets_itself(self, tmp_path, changes, label):
        assert f"sets {label} itself" in problems(tmp_path, replace(opt_excited(), **changes))

    def test_readfc_is_refused(self, tmp_path):
        assert "Force constants must be" in problems(tmp_path, replace(opt_excited(),
                                                                       force_constants="Read from checkpoint (ReadFC)"))

    def test_warnings(self, tmp_path):
        job = gi.GaussianInputBuilder().prepare(job_spec(tmp_path, opt_excited(excited_coefficients=False),
                                                          write_formchk=False))
        text = "\n".join(job.warnings)
        assert "IOp(9/40=4)" in text and "water_excited.fchk" in text

    def test_defaults(self):
        s = gi.opt_excited_defaults(step(state=gk.STATE_TD, geometry_source=gk.GEOM_CHECK, old_chk="x", write_chk=False,
                                         opt_target=gk.OPT_TS, force_constants="Read from checkpoint (ReadFC)"))
        assert (s.job_type, s.state, s.geometry_source, s.old_chk, s.write_chk, s.opt_target, s.force_constants) == \
            (gk.JOB_OPT_EXCITED, gk.STATE_GROUND, gk.GEOM_COORDINATES, "", True, gk.OPT_MINIMUM, gk.FC_DEFAULT)

    def test_preview_and_summary(self, tmp_path):
        spec = job_spec(tmp_path, opt_excited())
        assert gi.route_preview(spec, 0) == ("step 1:  #P B3LYP/6-31G(d) Opt\n"
                                             "step 2:  #P B3LYP/6-31G(d) TD=(NStates=10) IOp(9/40=4) Geom=AllCheck "
                                             "Guess=Read")
        assert gi.molecule_summary(spec, 0).endswith("step 2 reads the geometry, charge and multiplicity from water.chk")

    def test_reopened_file_is_two_linked_steps_with_a_note(self, tmp_path):
        builder = gi.GaussianInputBuilder()
        text = builder.input_text(builder.prepare(job_spec(tmp_path, opt_excited())))
        spec, notes = gi.read_input_text(text, "water", str(tmp_path))
        assert len(spec.steps) == 2 and spec.steps[1].chk == "water_excited.chk" and spec.steps[1].old_chk == "water.chk"
        assert any("Geometry optimization + excited states" in n for n in notes)
        again = builder.prepare(spec)
        assert builder.input_text(again) == text
```

- [ ] **Step 2: Run — expect AttributeError `opt_excited_defaults`**

- [ ] **Step 3: Implement in `gaussianInput.py`**

(a) The 15 spec fields after the Fukui fields; their `FIELD_LABELS`; `CHOICES` entries `excited_state: gk.EXCITED_STATES, excited_state_spin: tuple(gk.STATE_SPINS), excited_family: gk.EXCITED_FAMILIES, excited_method_name: gk.EXCITED_DFT_METHODS, excited_basis: gk.EXCITED_BASES, excited_basis_diffuse: tuple(gk.DIFFUSE), excited_basis_polarization: tuple(gk.POLARIZATION)`; `"excited_keywords"` and `"excited_chk"` added to `ASCII_FIELDS`.

(b) Refactor `prepare` (non-Fukui path) into `_prepare_steps(self, spec, labels=None, excited_analysis=False) -> GaussianJob`, where step messages are prefixed with `labels[i]` when given (else `f"Step {i}"`). `_prepare_step` gains `excited_analysis: bool = False`, which skips the `Density=Current` warning (the Excitation/NTO tabs use the MO coefficients). `prepare` dispatches `JOB_FUKUI` → `_prepare_fukui`, `JOB_OPT_EXCITED` → `_prepare_opt_excited`, else `_prepare_steps(spec)`.

(c) The composite:

```python
OPT_EXCITED_STEP_LABELS = ("Optimization step", "Excited-state step")
OPT_EXCITED_FIXED = ("opt_target", "force_constants", "state", "print_all_coefficients", "geometry_source",
                     "write_chk", "old_chk")
STEP_TWO_DROPPED = dict(job_type=gk.JOB_ENERGY, geometry_source=gk.GEOM_ALLCHECK, title="", modredundant="",
                        additional_input="", connectivity=False, guess_mix=False, guess_always=False,
                        guess_only=False, guess_save=False, pop_level="Default", pop_mk=False, pop_hirshfeld=False,
                        pop_nbo=False, pop_npa=False, pop_savento=False, write_chk=True)


def opt_excited_defaults(step):
    """The step as this job type uses it: optimization towards a minimum from the typed coordinates."""
    force_constants = step.force_constants if step.force_constants in gk.FUKUI_FORCE_CONSTANTS else gk.FC_DEFAULT
    return dataclasses.replace(step, job_type=gk.JOB_OPT_EXCITED, opt_target=gk.OPT_MINIMUM,
                               force_constants=force_constants, state=gk.STATE_GROUND, print_all_coefficients=False,
                               geometry_source=gk.GEOM_COORDINATES, write_chk=True, old_chk="")


def optimization_chk(step, stem):
    return step.chk.strip() or f"{stem}.chk"


def excited_chk(step, stem):
    return step.excited_chk.strip() or f"{stem}{gk.EXCITED_CHK_SUFFIX}.chk"


def _excited_model(step) -> dict:
    """Step 2's Method-panel fields (spec 2.3/2.4)."""
    fields = dict(state=step.excited_state, n_states=step.excited_n_states, root=step.excited_root,
                  state_spin=step.excited_state_spin, print_all_coefficients=step.excited_coefficients)
    if step.excited_family != gk.SAME_AS_OPTIMIZATION or step.excited_state == gk.STATE_CIS:
        fields["dispersion"] = "None"
    if step.excited_family != gk.SAME_AS_OPTIMIZATION:
        fields.update(family=step.excited_family, method_other=step.excited_method_other,
                      method_name=step.excited_method_name if step.excited_family == gk.FAMILY_DFT else gk.NO_CHOICE)
    if step.excited_basis == gk.SAME_AS_OPTIMIZATION:
        if step.basis in gk.GENERAL_BASES:
            fields["basis"] = gk.BASIS_CHECKPOINT
    else:
        fields.update(basis=step.excited_basis, basis_diffuse=step.excited_basis_diffuse,
                      basis_polarization=step.excited_basis_polarization, basis_other=step.excited_basis_other)
    return fields


def _reuses_orbitals(first, second) -> bool:
    if no_basis(first):
        return False
    if second.basis == gk.BASIS_CHECKPOINT:
        return True
    return basis_keyword(first).upper() == basis_keyword(second).upper()


def opt_excited_steps(step, stem):
    """(optimization, excited states): two ordinary steps for _prepare_steps (spec 2.2-2.4)."""
    first = dataclasses.replace(step, job_type=gk.JOB_OPT_FREQ if step.opt_confirm_freq else gk.JOB_OPT,
                                state=gk.STATE_GROUND, print_all_coefficients=False, write_chk=True)
    second = dataclasses.replace(step, **STEP_TWO_DROPPED, **_excited_model(step),
                                 old_chk=optimization_chk(step, stem), chk=excited_chk(step, stem),
                                 additional_keywords=" ".join(t for t in (_without_job_keywords(step.additional_keywords),
                                                                          step.excited_keywords.strip()) if t))
    second = dataclasses.replace(second, guess="Read checkpoint" if _reuses_orbitals(first, second) else "Default")
    return first, second
```

`_prepare_opt_excited(self, spec)`:
1. `len(spec.steps) != 1` → `ValueError("- The Geometry optimization + excited states job type is already two linked steps (optimization, excited states): it must be the only step.")`.
2. Fixed controls, as in `_prepare_fukui`: force constants → `"... Force constants must be Default, Calculate at first point (CalcFC) or Calculate at all points (CalcAll) (got ...)"`; others → `f"The Geometry optimization + excited states job type sets {FIELD_LABELS[name]} itself; choose the job type again to reset it."`.
3. Composite problems: Compound family → `"Optimization step: Compound methods (CBS-QB3, G4, ...) cannot be optimized here; choose Hartree-Fock, DFT, semi-empirical or MP2."`; RO prefix → `"Restricted-open (RO) references are not supported: Gaussian's TD/CIS and the Excitation and NTO tabs need a restricted or unrestricted reference; choose Default or U."`; `no_basis(first)` and excited basis Same → `f"Excited-state step: the optimization ({method_keyword(first)}) has no basis set; choose the excited-state basis set."`; excited chk equal to the optimization chk ignoring case → `f"Excited-state %Chk '{name}' is the optimization's %Chk; give the excited states their own checkpoint (blank = {stem}_excited.chk)."`; each `excited_keywords` token whose name is in `LINKED_DROPPED_KEYWORDS` → `f"Excited-state keywords: {token} belongs to the optimization step; this step computes the excited states only."`.
4. `job = self._prepare_steps(replace(spec, steps=[first, second]), OPT_EXCITED_STEP_LABELS, excited_analysis=True)`, with the step-3 errors raised together with `_prepare_steps`' (collect both, one `ValueError`).
5. Warnings appended: coefficients off → `"Excited-state step: without IOp(9/40=4) Gaussian prints only the larger excitation coefficients; the NTO and Excitation tabs then report a normalisation deviation."`; `not spec.write_formchk` → `f"No formchk script: the Excitation and NTO tabs need {stem}_excited.fchk; run formchk {excited chk} after the job."`; no orbital reuse → `f"Excited-state step: its basis set ({basis_keyword(second)}) differs from the optimization's ({basis_keyword(first)}), so it starts from a new guess, not from the optimization's orbitals."` or, when `no_basis(first)`, `f"Excited-state step: the optimization ({method_keyword(first)}) has no basis set, so it starts from a new guess."`.
6. `job.excited_fchk = Path(excited chk).stem + fi.FCHK_EXTENSION`.

(d) `route_preview`: for `JOB_OPT_EXCITED`, `first, second = opt_excited_steps(step, stem)` and return `"step 1:  " + route(first) + "\nstep 2:  " + route(second)` + the merge errors of both, one per line (helper `_route_and_errors(step) -> Tuple[str, List[str]]` shared with the plain path).

(e) `molecule_summary`: for `JOB_OPT_EXCITED` append `f"\nstep 2 reads the geometry, charge and multiplicity from {optimization_chk(step, stem)}"`.

(f) `read_input_text`: after reading, if `len(steps) == 2 and steps[1].geometry_source == gk.GEOM_ALLCHECK and steps[1].state != gk.STATE_GROUND and steps[1].chk.endswith(gk.EXCITED_CHK_SUFFIX + ".chk")`, add the note `"This is a Geometry optimization + excited states job written by the builder; it is shown as two linked steps. Choose that job type to rebuild it from its panels."`.

- [ ] **Step 4: Run** `tests/test_gaussianInput.py tests/test_gaussianInputReopen.py` — PASS (fix route expectations only if the keyword order in the file differs from the existing panel order, which is fixed: job, method, general).

- [ ] **Step 5: Commit** `feat(gaussian): Geometry optimization + excited states job type (two linked steps)`

---

### Task 5: The window

**Files:** Modify `gaussianInputBuilder.py`; Test `tests/test_gaussianInputBuilder.py` (new class `TestOptExcitedInTheWindow`)

**Interfaces:** Consumes everything of Task 4. Produces widget attributes named in `STEP_WIDGETS` below, `OPT_EXCITED_STATUS`, `OPT_EXCITED_NOTE`.

- [ ] **Step 1: Failing tests**

```python
class TestOptExcitedInTheWindow:
    def test_choosing_it_resets_what_it_sets(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "state_selection", gk.STATE_TD)
        change(ui, "old_chk_input", "x.chk")
        change(ui, "job_type_selection", gk.JOB_OPT_EXCITED)
        assert ui.state_selection.value == gk.STATE_GROUND and ui.old_chk_input.value == ""
        assert ui.status_label.text == gb.OPT_EXCITED_STATUS
        assert not ui.add_step_button.enabled and not ui.state_selection.enabled
        assert ui.excited_state_selection.enabled and ui.fukui_note_label.text == gb.OPT_EXCITED_NOTE

    def test_controls_follow_the_choices(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT_EXCITED)
        assert not ui.excited_method_selection.enabled and not ui.freq_raman_switch.enabled
        change(ui, "excited_family_selection", gk.FAMILY_DFT)
        change(ui, "opt_confirm_freq_switch", True)
        assert ui.excited_method_selection.enabled and ui.freq_raman_switch.enabled
        change(ui, "excited_state_selection", gk.STATE_CIS)
        assert not ui.excited_family_selection.enabled
        change(ui, "excited_basis_selection", "6-311G")
        assert ui.excited_diffuse_selection.enabled

    def test_route_line_shows_both_steps(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT_EXCITED)
        assert ui.route_label.text.startswith("step 1:  #P B3LYP/6-31G(d) Opt\nstep 2:  ")

    def test_save_says_where_the_excitation_files_are(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT_EXCITED)
        asyncio.run(ui.save_input(None))
        text = (tmp_path / "water.gjf").read_text()
        assert text.count("--Link1--") == 1 and "%Chk=water_excited.chk" in text
        assert "water_excited.fchk" in ui.messages_output.value and "water.log" in ui.messages_output.value
        assert "formchk water_excited.chk" in (tmp_path / "water_formchk.sh").read_text()

    def test_add_linked_step_is_refused(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT_EXCITED)
        ui.add_linked_step(None)
        assert len(ui.steps) == 1 and "already two linked steps" in ui.status_label.text
```

- [ ] **Step 2: Run — expect AttributeError `OPT_EXCITED_STATUS`**

- [ ] **Step 3: Implement**
  - `STEP_WIDGETS` += `opt_confirm_freq: opt_confirm_freq_switch BOOL`, `excited_state: excited_state_selection CHOICE`, `excited_n_states: excited_n_states_input INT`, `excited_root: excited_root_input INT`, `excited_state_spin: excited_state_spin_selection CHOICE`, `excited_family: excited_family_selection CHOICE`, `excited_method_name: excited_method_selection CHOICE`, `excited_method_other: excited_method_other_input TEXT`, `excited_basis: excited_basis_selection CHOICE`, `excited_basis_diffuse: excited_diffuse_selection CHOICE`, `excited_basis_polarization: excited_polarization_selection CHOICE`, `excited_basis_other: excited_basis_other_input TEXT`, `excited_coefficients: excited_coefficients_switch BOOL`, `excited_keywords: excited_keywords_input TEXT`, `excited_chk: excited_chk_input TEXT`.
  - Job Type panel: the section of spec §2.7 between NMR and Fukui (heading "Geometry optimization + excited states", the switch, the fields with the window's `_field`/`_choice`/`_text`/`_switch` helpers, and the hint "Step 1 optimizes with the Method panel's model chemistry and writes %Chk; step 2 reads the geometry, charge and multiplicity from it (Geom=AllCheck, %OldChk) and writes <name>_excited.chk for the Excitation and NTO tabs.").
  - Constants: `OPT_EXCITED_STATUS = "Geometry optimization + excited states: two linked steps. Choosing it reset the electronic state, the geometry source and %OldChk; set the excited states in the Job Type panel."`; `OPT_EXCITED_NOTE = "Step 1 is the optimization (Method panel); step 2 computes the excited states at its geometry with the settings above. The Method panel's electronic state stays Ground state."`
  - `on_field_change`: `elif widget is self.job_type_selection and widget.value == gk.JOB_OPT_EXCITED:` load `gi.opt_excited_defaults(step)` and `_status(OPT_EXCITED_STATUS)` (same shape as the Fukui branch).
  - `refresh_controls`: `opt_excited = job == gk.JOB_OPT_EXCITED`; force constants/tight/max cycles also when `opt_excited`; freq flags when `job in FREQUENCY_JOBS or (opt_excited and step.opt_confirm_freq)`; `state_selection`, `geometry_source_selection`, `write_chk_switch`, `old_chk_input`, `add_step_button` also `and not opt_excited`; modredundant widgets when `optimizing or opt_excited`; the new widgets: all `opt_excited`, plus `excited_family_selection: opt_excited and step.excited_state != gk.STATE_CIS`, `excited_method_selection: that and step.excited_family == gk.FAMILY_DFT`, `excited_method_other_input: that and step.excited_method_name == gk.OTHER`, diffuse/polarization `opt_excited and step.excited_basis in gk.EXTENSIBLE_POPLE`, `excited_basis_other_input: opt_excited and step.excited_basis == gk.OTHER`; `fukui_note_label` text `FUKUI_NOTE if fukui else OPT_EXCITED_NOTE if opt_excited else ""`.
  - `add_linked_step`: refuse `JOB_OPT_EXCITED` with `"The Geometry optimization + excited states job type is already two linked steps; it cannot take more."`.
  - `_run_instructions`: when `job.excited_fchk`: `f"Then open {job.excited_fchk} in the 3D Molecular Viewer's Excitation (or NTO) tab with {log} (the output of both steps; the excited states are its last block)."` instead of the generic line.

- [ ] **Step 4: Run** `tests/test_gaussianInputBuilder.py` — PASS (the table tests now cover the 15 fields).

- [ ] **Step 5: Commit** `feat(gaussian): Opt + excited states controls in the builder window`

---

### Task 6: Excitation-tab compatibility, real G16 verification

**Files:** Modify `excitationData.py`, `tools/check_gaussian_inputs.py`, `tests/data/gaussian_output/inputs/h2co_opt_td*.gjf` (+ regenerated `h2co_opt_td.log`, `h2co_opt_td_excited.fchk.gz`); Test `tests/test_excitationData.py`, `tests/test_gaussianOutput.py` (end to end)

- [ ] **Step 1: Failing tests**

```python
# tests/test_excitationData.py
def test_sibling_of_an_excited_checkpoint_is_the_jobs_output(tmp_path):
    fchk = tmp_path / "mol_excited.fchk"
    fchk.write_text("")
    (tmp_path / "mol.log").write_text("")
    assert ed.sibling_excitation_file(fchk) == tmp_path / "mol.log"
    (tmp_path / "mol_excited.out").write_text("")
    assert ed.sibling_excitation_file(fchk) == tmp_path / "mol_excited.out"


def test_excited_suffix_matches_the_builder():
    import gaussianKeywords as gk
    assert ed.EXCITED_SUFFIX == gk.EXCITED_CHK_SUFFIX

# tests/test_gaussianOutput.py
class TestBuilderJobs:
    @pytest.mark.parametrize("name,changes", [
        ("h2co_opt_td", {}),
        ("h2co_opt_td_gen", {"basis": "Gen", "additional_input": "C O H 0\n6-31G(d)\n****", "excited_n_states": 3}),
    ])
    def test_fixture_input_is_the_builders(self, tmp_path, name, changes):
        import gaussianInput as gi
        import gaussianKeywords as gk
        geometry = "\n".join(l for l in (DATA / "inputs" / f"{name}.gjf").read_text().splitlines()[7:11])
        step = gi.opt_excited_defaults(gi.GaussianStepSpec(multiplicity=1, memory="1GB", nproc=2, **changes))
        builder = gi.GaussianInputBuilder()
        job = builder.prepare(gi.GaussianJobSpec(stem=name, output_dir=str(tmp_path), geometry=geometry, steps=[step]))
        assert builder.input_text(job) == (DATA / "inputs" / f"{name}.gjf").read_text()

    def test_excitation_analysis_accepts_the_jobs_files(self, tmp_path):
        import gzip, shutil
        import excitationAnalysis as ea
        fchk = tmp_path / "h2co_opt_td_excited.fchk"
        with gzip.open(DATA / "h2co_opt_td_excited.fchk.gz", "rb") as src, open(fchk, "wb") as dst:
            shutil.copyfileobj(src, dst)
        shutil.copy(DATA / "h2co_opt_td.log", tmp_path / "h2co_opt_td.log")
        import excitationData as ed
        assert ed.sibling_excitation_file(fchk) == tmp_path / "h2co_opt_td.log"
        results = ea.ExcitationAnalyser().run(fchk, tmp_path / "h2co_opt_td.log", 1)
        assert results is not None
```

- [ ] **Step 2: Run — expect FAIL (sibling, fixture text differs from the hand-written input)**

- [ ] **Step 3: Implement** `excitationData.EXCITED_SUFFIX = "_excited"`; `sibling_excitation_file` tries the stem, then the stem without a trailing `_excited`, each with `.log` then `.out`.

- [ ] **Step 4: Regenerate the two `h2co_opt_td*` inputs from the builder** (write `builder.input_text(job)` exactly as the test builds it, same geometry lines), rerun `python tests/data/gaussian_output/run_gaussian.py`, confirm two Normal terminations each, then replace `h2co_opt_td.log` and `h2co_opt_td_excited.fchk.gz`.

- [ ] **Step 5: `tools/check_gaussian_inputs.py`** += `"opt_excited": (WATER, [gi.opt_excited_defaults(step(excited_n_states=3))], True)` and `"opt_excited_gen": (WATER, [gi.opt_excited_defaults(step(basis="Gen", additional_input="O H 0\nSTO-3G\n****", excited_n_states=3))], True)`; run it into an empty folder in WSL `/tmp`, run g16 on the two files: 2 Normal terminations each.

- [ ] **Step 6: Run** `tests/test_gaussianOutput.py tests/test_excitationData.py tests/test_gaussianInputReopen.py` — PASS

- [ ] **Step 7: Commit** `feat(gaussian): excited-state checkpoints open in the Excitation/NTO tabs; verified on G16`

---

### Task 7: Documentation and the whole suite

**Files:** `help.py` (`help_gaussian_input`), `USER_MANUAL.md` (*Gaussian 16 Input Builder*: Load geometry from outputs, the job type, the Excitation files; *Analyse an Electronic Excitation*: the `_excited.fchk` + `.log` pair; the 1.7.0 history paragraph gets a sentence), `CLAUDE.md` (the Gaussian builder entry: `gaussianOutput`, ChkBasis, the composite).

- [ ] **Step 1:** Write the docs (what each file is, the conversion `formchk <stem>_excited.chk` and which files the Excitation tab needs).
- [ ] **Step 2:** `..\Scripts\python.exe -m pytest tests/ -q` — only the 4 baseline failures.
- [ ] **Step 3: Commit** `docs(gaussian): optimized geometries from outputs and the Opt + excited states job type`
