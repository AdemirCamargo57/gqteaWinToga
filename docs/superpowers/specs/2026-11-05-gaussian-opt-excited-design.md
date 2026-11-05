# Gaussian 16 Input Builder: optimized geometries from outputs, and Opt + excited states

Date: 2026-11-05. Extends the Gaussian 16 Input Builder (spec `2026-11-03-gaussian-input-builder-design.md`).
Target: Gaussian 16 Rev. C.01 under WSL or Linux (`%OldChk` is required; G09W A.02 lacks it).

## 0. User decisions (2026-11-05)

| # | Question | Decision |
|---|----------|----------|
| 1 | How the new job appears | **Composite job type** (the Fukui pattern): one step in the step bar, dedicated excited-state controls, expanded into two `--Link1--` steps at Preview/Save; it must be the only step. |
| 2 | When step 2 reuses step 1's orbitals | **Same basis only**: `Guess=Read` only when step 2's basis is identical to step 1's and step 1 has a basis; otherwise a fresh default guess and a warning saying why. |
| 3 | Frequencies in step 1 | An **optional switch, off by default**: "Confirm the minimum with frequencies (Opt Freq)". |
| 4 | What step 2 shares | **Resources + General + Solvation**: %Mem, %NProcShared, print level, SCF/NoSymm/Polar, the solvation model and solvent, and the Additional Keywords minus job keywords; plus its own "Excited-state keywords" field. |
| 5 | Defaults | NStates = 10; IOp(9/40=4) on; Gen/GenECP in step 2 only through "Same as optimization" (`ChkBasis`). |

## 1. Optimized geometry from a Gaussian output (`.log`, `.out`)

### 1.1 Module

New Toga-free module `gaussianOutput.py` (subprocess test: `toga` stays out of `sys.modules`).

```python
@dataclass
class OptimizedGeometry:
    atoms: List[Tuple[str, float, float, float]]   # element, x, y, z in Angstrom, file order
    charge: Optional[int]
    multiplicity: Optional[int]
    job_number: int          # 1-based job of the file that holds the completed optimization
    job_count: int
    orientation: str         # "Input orientation" | "Standard orientation" | "Z-Matrix orientation"
    notes: List[str]

def read_optimized_geometry(path) -> OptimizedGeometry      # ValueError names the file
def parse_optimized_geometry(text: str, name: str) -> OptimizedGeometry
```

### 1.2 Rules

1. **A Gaussian output** contains the Gaussian banner (`Entering Gaussian System` or `Gaussian, Inc.`).
   Anything else (an ORCA `.out`, a text file) is refused: "<name> is not a Gaussian output".
2. **Jobs.** The text is split into jobs: each `--Link1--` step of the input and each internal job
   Gaussian chains itself (the Freq part of `Opt Freq`). The exact boundary lines are taken from real
   G16 C.01 outputs (§4) and pinned by tests.
3. **Completed optimization** = a job containing `Optimization completed.` The converged geometry is
   the **last orientation block before that line** in the job. Preference inside that geometry step:
   *Input orientation* (the input's frame), else *Standard orientation*, else *Z-Matrix orientation*.
   Atom order is the file's.
4. **Which one**: the last completed optimization of the **last job that has one**. A relaxed scan
   (several completions in one job) gives its last point, with a note.
5. **Charge and multiplicity**: the last `Charge = q Multiplicity = m` line of that job before the
   completion (Geom=Check steps print it too). Missing → `None` (the window leaves the fields as they are).
6. **Notes** (shown in the status line and the messages box; never silent):
   - a later job runs an optimization that did not complete → "job k's optimization did not complete;
     the geometry of job j's completed optimization was imported";
   - the job did not end with `Normal termination` after the optimization completed;
   - a relaxed scan (rule 4).
7. **Errors** (ValueError → "Geometry Error" dialog; the molecule panel is not changed):
   - no optimization in the file: "<name> contains no geometry optimization; load its input or an .xyz";
   - an optimization that never completed, with the reason found: `Number of steps exceeded`,
     `Error termination`, or the file ends during the optimization ("still running or cut off");
     the message says the last geometry is **not** imported;
   - an orientation block with atomic number ≤ 0 (ghost/dummy atom) or that cannot be parsed.

### 1.3 Integration

- `gaussianInput.GEOMETRY_FILE_TYPES` gains `log`, `out`. A new `load_geometry_file(path) ->
  GeometryFile(text, charge, multiplicity, notes)` holds the dispatch (`.log`/`.out` → `gaussianOutput`,
  the rest as now with `notes=[]`); `read_geometry_file` stays as the 3-tuple wrapper over it, so its
  callers and tests are unchanged. The window calls `load_geometry_file`.
- The Molecule panel's hint names the new suffixes; `apply_geometry_file` reports the notes and
  `Optimized geometry of job j loaded from <name>`.
- Every other format keeps its current reader. `molecularFileFormats` (Design canvas import) is unchanged.

## 2. Job type "Geometry optimization + excited states"

### 2.1 Model

`gk.JOB_OPT_EXCITED = "Geometry optimization + excited states"`, in `JOB_TYPES` before the Fukui type.

New `GaussianStepSpec` fields (each with a widget in `STEP_WIDGETS`, a `FIELD_LABELS` entry and, for
lists, a `CHOICES` entry; the existing table tests then cover them):

| Field | Default | Meaning |
|-------|---------|---------|
| `opt_confirm_freq` | False | step 1 writes `Opt Freq` |
| `excited_state` | TD-SCF | TD-SCF / TDA / CIS |
| `excited_n_states` | 10 | NStates |
| `excited_root` | 1 | Root |
| `excited_state_spin` | Default | Default / Singlets / Triplets / 50-50 |
| `excited_family` | Same as optimization | Same as optimization / Hartree-Fock / DFT |
| `excited_method_name` | B3LYP | DFT list (when the family is DFT) |
| `excited_method_other` | "" | when the method is Other... |
| `excited_basis` | Same as optimization | Same as optimization + `BASES` without Gen/GenECP |
| `excited_basis_diffuse`, `excited_basis_polarization`, `excited_basis_other` | None, (d), "" | as in the Method panel |
| `excited_coefficients` | True | `IOp(9/40=4)` |
| `excited_keywords` | "" | route keywords for step 2 only, merged like Additional Keywords |
| `excited_chk` | "" | blank = `<stem>_excited.chk` |

`GaussianInputBuilder.prepare` routes a spec whose step is `JOB_OPT_EXCITED` to `_prepare_opt_excited`,
which refuses more than one step, then builds the two steps with the pure function
`opt_excited_steps(step, stem) -> (step1, step2)` and passes each through the existing `_prepare_step`
(so Link 0, method, multiplicity, title, file-name and keyword-merge validation is shared).

### 2.2 Step 1 (optimization)

`step1 = replace(step, job_type=Optimization or Opt+Freq, state=Ground, write_chk=True)`.
Everything else comes from the panels: Method, Molecule (coordinates, charge, multiplicity, freeze
lines only), Title, Link 0, General, Guess, Pop, Solvation, Add. Inp., Additional Keywords.
`%Chk` = the Link 0 field, blank = `<stem>.chk`. `Opt` (and `Freq`) are route keywords.

### 2.3 Step 2 (excited states)

- Link 0: `%OldChk=<step-1 chk>`, `%Chk=<excited chk>`, `%Mem`, `%NProcShared`, the extra Link 0 lines.
- Route: print level; the excited state (`TD(NStates=..,Root=..,Singlets)`, `TDA(...)` or
  `CIS(...)/basis`); method and basis; `Geom=AllCheck` (geometry, charge, multiplicity and title from the
  checkpoint, so no title or molecule section); `Guess=Read` per §2.4; `IOp(9/40=4)`; SCF/NoSymm/Polar
  and SCRF copied from step 1; dispersion copied only when the method is "Same as optimization";
  Additional Keywords minus job keywords (`LINKED_DROPPED_KEYWORDS`) merged with Excited-state keywords.
- Not copied: Opt/Freq, freeze lines, connectivity, Guess panel flags (Mix, Always, Only, Save), Pop
  panel, Add. Inp. (except through `ChkBasis`, §2.4).
- Spin prefix: step 1's. The step-1 multiplicity is kept in the step-2 spec for the spin checks (not written).

### 2.4 Basis and wavefunction reuse

| Step 1 basis | Step 2 basis | Step 2 writes | Guess |
|--------------|--------------|---------------|-------|
| any listed / Other | Same as optimization | the same basis keyword | `Guess=Read` |
| any listed / Other | a different one | that basis | default (warning: new guess) |
| Gen / GenECP | Same as optimization | `ChkBasis` (basis and ECP from the checkpoint) | `Guess=Read` |
| none (semi-empirical) | Same as optimization | error: give step 2 a basis set | — |
| none (semi-empirical) | a basis | that basis | default (warning) |

"Same" compares the written basis keyword (case-insensitive). `ChkBasis` is confirmed on G16 (§4).

### 2.5 Validation

Errors (one `- ` line each, labelled "Optimization step:" or "Excited-state step:"):
- the job type is not the only step;
- the Method panel's Electronic state is not Ground state (the panel greys it out; the excited state
  is set in the Job Type tab);
- a Compound method for the optimization;
- an RO spin prefix (Gaussian TD and the Excitation/NTO tabs do not support it);
- TD-SCF with a family other than HF/DFT, TDA without DFT (existing rules, on the resolved step 2);
- Singlets/Triplets/50-50 with an open-shell step 1 (existing rule);
- the excited `%Chk` equal to step 1's `%Chk` (case-insensitive) or `%OldChk`-style errors on the names;
- scan lines (S) (Optimization takes freeze lines only — existing rule);
- everything `_prepare_step` already checks.

Warnings: coefficients off ("the NTO and Excitation tabs need IOp(9/40=4) for complete coefficients");
the formchk switch off ("the Excitation tab needs <stem>_excited.fchk"); a new guess (§2.4).

### 2.6 Text, preview, files

- `input_text` writes the two prepared steps joined by `--Link1--` (existing code).
- `route_preview` shows `step 1: ...` and `step 2: ...` lines, as for the Fukui job.
- `molecule_summary` adds "step 2 reads the geometry, charge and multiplicity from <stem>.chk".
- `chk_files` = [step-1 chk, excited chk]; the formchk script converts both (`<stem>.fchk`,
  `<stem>_excited.fchk`).
- Run instructions after Save end with: open `<stem>_excited.fchk` in the 3D Molecular Viewer's
  Excitation (or NTO) tab with `<stem>.log`.
- Open input reads such a file as two linked steps with a note (the Fukui precedent).

### 2.7 Window

- New Job Type section "Excited states (after the optimization)" with the controls of §2.1; enabled only
  for this job type. The Optimization section stays enabled (target fixed to Minimum, force constants
  without ReadFC — the Fukui rules).
- Choosing the job type calls `opt_excited_defaults(step)`: electronic state Ground, geometry source
  Coordinates, `%OldChk` blank, write `%Chk` on, opt target Minimum, ReadFC → Default; status message.
- "Add linked step" is disabled for this job type.

## 3. Compatibility with the Excitation and NTO tabs

- The tabs read a **formatted checkpoint** (`<stem>_excited.fchk`, made with `formchk` from the binary
  `.chk`, which they cannot read) **and** the **Gaussian output of the same job** (`<stem>.log`, the
  output of the whole linked job; its last excitation block and last orientation are step 2's).
- `excitationData.sibling_excitation_file(fchk)` also tries the fchk stem without a trailing
  `_excited`, so the tabs fill `<stem>.log` for `<stem>_excited.fchk`.
- Documented in `help_gaussian_input`, the run instructions and USER_MANUAL.

## 4. Fixtures and tests

Real G16 C.01 runs in WSL (from `/tmp`), stored in `tests/data/gaussian_output/` with their `.gjf`
inputs and a dev-time `run_gaussian.py` (never run by the suite):

| Fixture | Content |
|---------|---------|
| `h2o_opt.log` | converged optimization |
| `h2o_opt_freq.out` | Opt Freq (two internal jobs), `.out` suffix |
| `h2co_opt_maxcycles.log` | optimization stopped by MaxCycles |
| `h2o_opt_then_failed.log` | a converged optimization, then a linked optimization that fails |
| `h2co_opt_td.log` + `h2co_opt_td_excited.fchk` | the builder's own job of this type |

Synthetic cases in the tests: a converged log truncated before completion, a non-Gaussian `.out`,
ghost atoms. Tests:
- parser: every fixture and synthetic case (geometry, orientation, charge, multiplicity, notes, errors);
- `read_geometry_file` / `apply_geometry_file` with `.log` and `.out`;
- golden input text of the job type (two blocks, distinct checkpoints, `Geom=AllCheck`, `Guess=Read`
  present/absent, `ChkBasis`, `IOp(9/40=4)`, shared keywords, dropped job keywords);
- validation matrix of §2.5;
- window: every new spec field from a widget (existing table tests), choosing the job type, add-step disabled;
- `sibling_excitation_file` with `_excited`;
- end to end: `ExcitationAnalyser` loads `h2co_opt_td.log` + `h2co_opt_td_excited.fchk` (geometry check passes);
- `tools/check_gaussian_inputs.py` gains `opt_excited` and `opt_excited_gen` (ChkBasis), run on G16:
  one Normal termination per step.

Docs: `help_gaussian_input`, USER_MANUAL (*Gaussian 16 Input Builder*, *Analyse an Electronic
Excitation*), CLAUDE.md.

## 5. Not included

Excited-state optimizations (`Opt` + `TD` in one step), state-specific solvation, NTO/Density=Current
switches (available through Excited-state keywords), rebuilding the composite job from a reopened file,
reading geometries from IRC/scan outputs other than §1.2 rule 4.
