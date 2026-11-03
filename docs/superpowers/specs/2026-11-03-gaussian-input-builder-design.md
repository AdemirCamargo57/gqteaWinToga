# Gaussian 16 Input Builder — design

Status: design approved 2026-11-03, not implemented.

A general Gaussian 16 input builder on the main window's **Inputs** tab, modelled on GaussView 6.1's
*Gaussian Calculation Setup* dialog. It replaces the *Gaussian Fukui Input Builder* button: the Fukui
linked job becomes one of its Job Types.

Reference used: GaussView 6.1.1 as installed in `C:\G09W` and its help pages `help\howto_gsetup.htm`
(dialog layout, panels, buttons), `help\gausscalcopts.htm` (every control → keyword) and
`help\gausscalc.htm` (calculation schemes).

## 0. User decisions

| Question | Decision |
|---|---|
| Scope of the first version | **Core panels**: Job Type, Method, Title, Link 0, General, Guess, Pop, Solvation, Add. Inp., Preview, Additional Keywords. No ONIOM, PBC, QST2/QST3, molecule groups or MM methods. |
| Where the Fukui builder goes | **A Job Type** ("Fukui functions (N, N+p, N−q)") sharing the Method / Link 0 / Molecule panels; `fukuiInput.py` is reused unchanged. The old button and window are removed. |
| After writing the input | **Write files only** (`.gjf` with LF endings and relative names for G16 under WSL or a cluster, optional formchk `.sh`), then show the commands to run. No Submit / Quick Launch. |
| Calculation schemes | **None**: a *Defaults* button only; reopening a saved `.gjf` restores its settings. |
| Molecule source | **File or typed text** (no link to the 3D viewer or the Design canvas). |
| Reopening a `.gjf` | **Yes, best effort**: recognised keywords go back to their panels, everything else to Additional Keywords / Add. Inp. |
| Scan / freeze coordinates | **Helper row + text box**: the row appends validated ModRedundant lines to an editable box; typed lines are validated too. |
| Linked steps (`--Link1--`) | **Full linked steps**: a step list; each step has its own panel state; Open input restores every step. |
| Multiplicity | **Required, never assumed** (the Fukui builder's rule, presented in the design and not objected to). |
| Additional Keywords naming a panel keyword | **Merged** into the panel keyword's options; an error only when the same valued option gets two values (§5.1). |

## 1. Modules

One-way chain, like the rest of the project:

```
toolLoader -> gaussianInputBuilder (window; the only Toga import)
                 -> gaussianInput (spec, validation, input text, reopen, file writing)
                      -> gaussianKeywords (data tables)
                      -> fukuiInput (unchanged)
                      -> molecularFileFormats, systemProperties, help (existing)
```

| File | Role |
|---|---|
| `gaussianKeywords.py` | No Toga. Declarative tables: job types and their options, method families, functionals, semi-empirical and compound methods, basis sets (Pople base + diffuse + polarization; Dunning; Def2; ECP bases), dispersion, General / Guess / Pop options, solvation models and solvents. Each entry gives its label, its route form, and what the reverse mapping recognises. |
| `gaussianInput.py` | No Toga. `GaussianStepSpec` (one step's panel state), `GaussianJobSpec` (job name, output folder, step-1 molecule, steps), `GaussianInputBuilder.prepare(spec) -> GaussianJob` (all problems in one `ValueError`, one line each, "Step n: …"), `route_keywords(step)`, `input_text(job)`, `script_text(job)`, `write_files(job)`, `read_input_file(path) -> (GaussianJobSpec, notes)`, the ModRedundant line parser/formatter, `read_geometry_file` and `wsl_path` (moved from `fukuiInputBuilder.py`). |
| `gaussianInputBuilder.py` | `GaussianInputBuilderUI`: one `_build_*` method per panel, `_collect_inputs` → spec, `_load_step` spec → widgets, dialogs. |
| `fukuiInput.py` | Unchanged. The generic helpers (`route_tokens`, `keyword_name`, `strip_prefix`, `parse_memory`, `parse_nproc`, `validate_stem`, `parse_geometry`, `format_geometry`, `TITLE_FORBIDDEN`, …) stay here and are imported by `gaussianInput`. |
| `fukuiInputBuilder.py` | **Removed.** Its tests move to `tests/test_gaussianInputBuilder.py`. |
| `toolLoader.py` | `ToolSpec("Gaussian 16 Input Builder", "gaussianInputBuilder", "GaussianInputBuilderUI")` replaces the Fukui entry on the Inputs tab; `_frozen_build_imports` updated. |
| `help.py` | `HelpGqteaWin.help_gaussian_input` (contains the former `help_fukui_input` text as its Fukui section); references to "Inputs > Gaussian Fukui Input Builder" updated, including `help_fukui`. |

Lengths are Å throughout (Gaussian input coordinates). Files: LF endings, relative `%Chk`/`%OldChk` names,
written with `newline="\n"` like `fukuiInput.write_files`.

## 2. Window

```
┌ Gaussian 16 Input Builder ──────────────────────────────────────────────────┐
│ Step: [2: TD-SCF B3LYP/6-31+G(d), Geom=AllCheck ▾] [Add linked step]        │
│       [Remove] [▲] [▼]                                                       │
│ Route: #P TD=(NStates=10) B3LYP/6-31+G(d) Geom=AllCheck Guess=Read           │
│ [Job Type|Method|Molecule|Title|Link 0|General|Guess|Pop|Solvation|Add. Inp.|│
│  Preview|Help]                                                               │
│ ┌ selected panel (ScrollContainer) ─────────────────────────────────────┐   │
│ └───────────────────────────────────────────────────────────────────────┘   │
│ Additional Keywords: [______________________________________]               │
│ Job name [molecule]   Output folder [.................] [Browse]            │
│ [Open input...] [Preview] [Save input] [Defaults] [Close]                   │
│ status line                                                                 │
└─────────────────────────────────────────────────────────────────────────────┘
```

- The **route line** shows the selected step's route, rebuilt on every change and wrapped as
  `fukuiInput.build_route` wraps it (`ROUTE_WIDTH` 72).
- **Steps.** Selecting a step stores the panels into the current step, then loads the selected one.
  Widget writes go through write-on-change (`cubeTab._assign`) and the handlers are synchronous and
  `_syncing`-guarded (the viewer tabs' no-flicker rule), so loading a step does not replay as edits.
  *Add linked step* copies the selected step and presets: geometry source "Checkpoint, all"
  (`Geom=AllCheck`), `Guess=Read`, Job Type Energy, same %Mem / %NProcShared / %Chk. *Remove* refuses the
  last remaining step. ▲/▼ reorder; the checkpoint rules (§5) are re-checked on Preview/Save.
- **Buttons.** *Open input...* (§6). *Preview* shows the whole input and the warnings on the Preview tab.
  *Save input* writes `<stem>.gjf` (Fukui: `<stem>_fukui.gjf`) and, when ticked, `<stem>_formchk.sh`; it
  asks before replacing existing files, then shows the WSL / cluster commands (as the Fukui builder did:
  `cd '<wsl path>'`, `g16 < x.gjf > x.log`, `bash x_formchk.sh`). *Defaults* resets every step to one
  default step, asking first when anything differs from the defaults. *Close*.
- Every `ValueError` becomes an `ErrorDialog`; nothing raises out of an async handler.

## 3. Panels (per step) and the keywords they write

Everything not covered below goes through Additional Keywords (§5.1).

### Job Type

| Job type | Options | Writes |
|---|---|---|
| Energy | — | — |
| Optimization | target Minimum / TS (Berny); force constants default / CalcFC / CalcAll / ReadFC; Tight; MaxCycles | `Opt`, `Opt=(TS,CalcFC,Tight,MaxCycles=n,…)` (+ `ModRedundant` when freeze lines exist) |
| Frequency | Raman, VCD, Anharmonic | `Freq`, `Freq=(Raman,VCD,Anharmonic)` |
| Opt+Freq | both rows above | `Opt… Freq…` |
| IRC | direction Both / Forward / Reverse; MaxPoints; StepSize; force constants CalcFC / CalcAll / ReadFC | `IRC=(CalcFC,Forward,MaxPoints=n,StepSize=n)` |
| Scan (relaxed) | the scan/freeze lines of the Molecule panel | `Opt=ModRedundant` (+ Optimization options) |
| Stability | Reoptimize | `Stable`, `Stable=Opt` |
| NMR | GIAO / CSGT / IGAIM / All; SpinSpin | `NMR`, `NMR=(CSGT,SpinSpin)` |
| Fukui functions (N, N+p, N−q) | p, q, multiplicities of N+p and N−q, Guess=Read for the ions, Optimization options of the N state (Minimum only; default / CalcFC / CalcAll; Tight; MaxCycles) | delegated to `fukuiInput` (§7) |

### Method

- **State**: Ground; TD-SCF; TDA; CIS; ZINDO. Excited states take NStates, Root and Singlets / Triplets /
  50-50 (`TD=(NStates=10,Root=1,Singlets)`, `CIS=(…)`, `ZINDO=(…)`). Option *Print all excitation
  coefficients for the NTO / Excitation tabs* writes `IOp(9/40=4)`.
- **Spin**: Default / R / U / RO (method prefix).
- **Family**: HF; DFT (B3LYP, CAM-B3LYP, ωB97X-D `wB97XD`, M06-2X, M06, PBE0 `PBE1PBE`, PBEPBE, B3PW91,
  BP86, TPSSTPSS, APFD, HSEH1PBE, LC-ωHPBE `LC-wHPBE`, Other = typed); Semi-empirical (PM7, PM6, AM1, PM3,
  PDDG; no basis); MP2; MP4(SDTQ); CCSD; CCSD(T); Compound (CBS-QB3, CBS-4M, G4, G4MP2, W1BD; no basis).
- **Basis**: STO-3G, 3-21G, 6-31G, 6-311G with diffuse (none, +, ++) and polarization ((d), (d,p),
  (2d,p), (2df,2pd), (3df,3pd)) dropdowns; cc-pVDZ/TZ/QZ, aug-cc-pVDZ/TZ/QZ; Def2SVP, Def2TZVP, Def2TZVPP,
  Def2QZVP; LanL2DZ; SDD; DGDZVP; Gen; GenECP; Other = typed.
- **Empirical dispersion**: None / GD3 / GD3BJ / GD2 (`EmpiricalDispersion=GD3BJ`).

### Molecule

- **Step 1**: geometry box (Cartesian `El x y z`, Å; pasted `.xyz` accepted) and *Load geometry...*
  (every format `molecularFileFormats` reads, plus `.fchk`; a `.gjf`/`.fchk` also fills charge and
  multiplicity — `read_geometry_file`, moved from the Fukui window). Charge and multiplicity fields; a
  live line shows the electron count and the multiplicities it allows (as the Fukui window's state line).
- **Steps 2+**: *Geometry source* replaces the box: **Checkpoint, all** (`Geom=AllCheck`; no title,
  charge or coordinates), **Checkpoint, new charge/multiplicity** (`Geom=Check`; title and
  charge/multiplicity, no coordinates), **Same geometry as step 1** (coordinates written again).
- **Scan / freeze helper**: type Bond / Angle / Dihedral, atom numbers, action Scan / Freeze, steps, step
  size, *Add* → appends `B 1 2 S 10 0.05`, `D 1 2 3 4 F`, … to an editable ModRedundant box. The box is
  parsed on every Preview/Save (§5).

### Title

Free text (one line); blank → the job name. `fukuiInput.TITLE_FORBIDDEN` characters are an error.

### Link 0

%Mem (`parse_memory`), %NProcShared (`parse_nproc`), %Chk (on by default, `<stem>.chk`), %OldChk
(optional), Extra Link 0 lines (free text, each must start with `%`).

### General

Print level (`#P` default, `#`, `#N`, `#T`); `SCF=QC`; `SCF=Tight`; `NoSymm`; `Polar`; *Write
connectivity* (`Geom=Connectivity`, bonds from `molecularFileFormats.detect_bonds`, written as single
bonds; off by default); *Write a formchk script* (on by default; one `.sh` converting every distinct
`.chk` the job writes, for the viewer's ESP / QTAIM / Orbitals / NTO / Excitation tabs).

### Guess

Default / Core / Huckel / Read / TCheck; Mix; Always; Only; Save (`Guess=(Read,Mix)`).

### Pop

Output level Default / Full / Regular / Minimal / None; MK; Hirshfeld (CM5); NBO; NPA; SaveNTO
(`Pop=(Full,MK,NBO)`).

### Solvation

None / IEFPCM / CPCM / SMD; solvent list (Gaussian's names) + Other (`SCRF=(SMD,Solvent=Water)`).

### Add. Inp.

Free text written after the geometry, the connectivity block and the ModRedundant lines, separated by the
blank lines Gaussian expects; the step ends with a blank line.

### Preview / Help

The whole input text plus the warnings; `help_gaussian_input`.

## 4. Linked steps

- A job is an ordered list of `GaussianStepSpec`; only the job name, output folder and step-1 molecule
  are job-wide. Steps are separated by `--Link1--`.
- **Checkpoints**: every step defaults to `%Chk=<stem>.chk`; a step may set its own %Chk and %OldChk.
- **The Fukui Job Type** is already a 3-step chain and must be the only step (§5).
- **formchk script**: converts every distinct `.chk` named by a `%Chk`; a warning notes that steps
  sharing one `.chk` leave only the last step's result in it.

## 5. Validation

`prepare` collects every problem; each message names its step.

**Errors**

- Job name (`validate_stem`), output folder required.
- Charge and multiplicity integers; multiplicity required and checked with
  `systemProperties.multiplicity_problem` (not needed for `Geom=AllCheck` steps); `R` prefix with
  multiplicity > 1.
- Basis required, except semi-empirical, compound and ZINDO, which refuse one; typed Other method /
  functional / basis / solvent must be one token without `/`; Gen / GenECP with an empty Add. Inp.
- Excited states: NStates ≥ 1, Root ≤ NStates; TD / TDA need HF or DFT (TDA DFT only); Singlets /
  Triplets / 50-50 need a closed-shell reference; `IOp(9/40=4)` only with an excited state.
- Scan needs ≥ 1 Scan line; Freeze lines only with Optimization, Opt+Freq, Scan; atom numbers 1…N,
  distinct, 2 / 3 / 4 for B / A / D; steps a positive integer; step size finite and non-zero.
- Geom=Check / AllCheck, Guess=Read / TCheck, ReadFC / RCFC need a checkpoint written by an earlier step
  (its `%Chk`) or named in `%OldChk`; step 1 can read one only through `%OldChk`.
- Write connectivity together with a checkpoint geometry.
- `%Chk` / `%OldChk` not a relative name without blanks; Extra Link 0 lines not starting with `%`.
- Fukui together with other steps; Title with forbidden characters.
- Geometry: `fukuiInput.parse_geometry` (atoms closer than 0.4 Å refused).

**Warnings**

- TS optimization without CalcFC / ReadFC / CalcAll.
- IRC without force constants (exact G16 behaviour checked against the Gaussian manual during
  implementation; an error instead if G16 refuses it).
- Post-SCF or excited-state method with the formchk script: the `.fchk` holds the SCF density unless
  `Density=Current` is given.
- Several steps sharing one `.chk` (§4).

### 5.1 Additional Keywords

Tokens split with `fukuiInput.route_tokens`. A token whose keyword a panel already writes is **merged**:
its options join the panel keyword's options, duplicates dropped (`Opt=(TS,CalcFC)` + `Opt=NoEigenTest`
→ `Opt=(TS,CalcFC,NoEigenTest)`). The same valued option with two values (`MaxCycles=50` vs the panel's
100, `Solvent=Ethanol` vs Water) is an error naming both places. Non-valued options that exclude each
other are not checked (Gaussian reports them). Merging is what makes reopening lossless (§6).

## 6. Reopening a `.gjf` / `.com` / `.gau`

`read_input_file(path) -> (GaussianJobSpec, notes)`:

1. Split on `--Link1--` lines (case-insensitive).
2. Per step: `%` lines → Link 0 (unknown ones → Extra Link 0 lines); `#` lines up to a blank line → the
   route; then title and charge/multiplicity unless `Geom=AllCheck`; Cartesian geometry unless
   `Geom=Check`/`AllCheck` (read with `molecularFileFormats`' Gaussian reader); the connectivity block if
   `Geom=Connectivity`; ModRedundant lines if `Opt=ModRedundant`; the rest → Add. Inp.
3. Route tokens are mapped through the `gaussianKeywords` tables: `B3LYP/6-31+G(d)` is split into
   prefix, family, functional, Pople base, diffuse and polarization; unmodelled options and keywords go to
   Additional Keywords (merged back on save, §5.1).
4. Refused with a message: Z-matrix geometries, freeze flags or fragment tags in the coordinates,
   per-fragment charges.
5. A `_fukui.gjf` reopens as 3 ordinary linked steps, with a note; it is not turned back into the Fukui
   Job Type.

Round trip: spec → text → reopen gives the same spec; text → spec → text is stable after the first
normalisation (keyword order may change).

## 7. Fukui Job Type

`gaussianInput` builds a `fi.FukuiJobSpec` from the panels and hands it to `fi.FukuiInputBuilder`
(`prepare`, `input_text`, `script_text`, `write_files`), so the output is that module's, unchanged.

| `FukuiJobSpec` field | From |
|---|---|
| `stem`, `output_dir` | job name, output folder |
| `geometry`, `charge`, `multiplicity` | Molecule panel |
| `plus_multiplicity`, `minus_multiplicity`, `p`, `q`, `guess_read` | Fukui options of the Job Type panel |
| `opt_keyword` | Optimization options of the N state (`Opt`, `Opt=(Tight,CalcFC,MaxCycles=n)`) |
| `method` | bare family / functional token (no R/U/RO) |
| `basis` | Basis |
| `keywords` | dispersion + General (`SCF=Tight`, `NoSymm`, `SCF=QC`, `Polar`) + Pop + Solvation + Additional Keywords |
| `print_level` | General |
| `memory`, `nproc` | Link 0 |

Choosing Fukui switches `SCF=Tight` and `NoSymm` on (the former defaults) and disables, with a line
saying why: Guess and Add. Inp. panels, Title, %Chk / %OldChk, Write connectivity, Spin prefix, excited
states, Add linked step. `fi.check_route` and `prepare` keep refusing what they refuse today.

## 8. Testing

- `tests/test_gaussianKeywords.py`: unique labels and keywords; every option has a route form and a
  reverse mapping.
- `tests/test_gaussianInput.py`: the exact route of every Job Type and option; whole-input references in
  `tests/data/gaussian_input/` for representative jobs (including linked Opt → TD); one test per error
  and warning; checkpoint-chain rules; scan / freeze lines; the merge rule; round trips for every
  reference; every `.gjf` already in `tests/data` reopens without error and re-emits the same keywords
  and options; the Fukui Job Type's output equals `fi.FukuiInputBuilder` called directly.
- `tests/test_gaussianInputBuilder.py` (headless, `__new__` + stub widgets): every `GaussianStepSpec`
  field and every `FukuiJobSpec` field comes from the window (the second migrated from
  `test_fukuiInputBuilder.py`); step switching stores and loads; write-on-change; Defaults confirmation;
  Add linked step presets; the remaining Fukui-window tests migrated.
- Subprocess test: `toga` stays out of `sys.modules` after importing `gaussianInput` and
  `gaussianKeywords`.
- `tests/test_toolLoader.py` follows the new `TOOL_TABS` entry.
- Live: open the real window; if G16 under WSL is available, run water HF/STO-3G jobs for each Job Type
  plus a linked Opt → TD job and require Normal termination.

## 9. Docs

`help_gaussian_input`; USER_MANUAL: a new *Gaussian 16 Input Builder* section, the *Gaussian Fukui Input
Builder* section rewritten as its Fukui Job Type, every "Inputs > Gaussian Fukui Input Builder"
reference updated; `readme.md`; `docs/index.html`; a CLAUDE.md entry.

## 10. Not implemented

ONIOM, PBC, QST2 / QST3, molecule groups, MM methods, rigid (Z-matrix) scans, calculation schemes,
Submit / Quick Launch, the 3D viewer / Design canvas as a molecule source, per-fragment charges, Z-matrix
input, restoring a Fukui file as the Fukui Job Type.
