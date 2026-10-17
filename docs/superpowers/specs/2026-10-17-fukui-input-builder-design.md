# Gaussian Fukui Input Builder — design

Status: implemented 2026-10-17, then **revised the same day for Gaussian 16 (§8)**. §8 supersedes §0–§5 wherever they differ: the G09W single-point design in §0–§7 is kept as the record of the first version.

## 0. User decisions

| Question | Decision |
|---|---|
| Checkpoint conversion | **Both**: a `.bat` script written next to the input, and a *Convert .chk -> .fchk* button that runs formchk locally. |
| Where the jobs run | **This Windows PC (G09W)**: `.bat` script, CRLF line endings, absolute `%Chk` paths by default. |
| Geometry workflow | **Single points only** on an already optimized N-state geometry; `Opt` and similar keywords are refused. |
| Placement | **Inputs tab** button, plus auto-fill in *Analysis of Fukui Functions* (`<stem>_N.fchk` fills `_Np1`/`_Nm1`). |
| `%Chk` paths / default route | Not answered; the proposals were used: absolute paths (switchable), `#P B3LYP/6-31G(d) SCF=Tight NoSymm`. |

## 1. Modules

| File | Role |
|---|---|
| `fukuiInput.py` | No Toga. `FukuiJobSpec` → `FukuiInputBuilder.prepare` → `FukuiJob`; `input_text`, `script_text`, `write_files`; `checkpoint_pairs`, `convert_checkpoints`, `check_ready_for_analysis`, `sibling_state_files`. |
| `fukuiInputBuilder.py` | `FukuiInputBuilderUI` (tabs Setup / Molecule & states / Preview & convert / Help) and `read_geometry_file`. |
| `gaussianUtilities.py` | `find_gaussian_utility(name, explicit)`; `electrostaticSurfacePotential.find_cubegen` delegates to it. |
| `fukuiFunctionAnalyze.py` | `fill_sibling_states()` after *Browse* on the N file. |
| `gqteaWinToga.py` | `("Gaussian Fukui Input Builder", FukuiInputBuilderUI)` on the Inputs tab. |
| `help.py` | `help_fukui_input`; `help_fukui` points to the builder. |

## 2. Fields, defaults, validation

| Field | Default | Rule |
|---|---|---|
| Job name | loaded file's stem (minus `_N`), else `molecule` | `[A-Za-z0-9][A-Za-z0-9._-]*`, ASCII |
| Output folder | loaded file's folder | must exist; required |
| `%Mem` | `2GB` | integer + `KB/MB/GB/TB/KW/MW/GW/TW`; a bare number is refused (words) |
| `%NProcShared` | `min(4, cpu_count)` | positive integer; warns above `cpu_count` |
| Print level | `#P` | `#P # #N #T` |
| Method / basis / keywords | `B3LYP` / `6-31G(d)` / `SCF=Tight NoSymm` | §3 |
| Geometry | — | Cartesian `El x y z` Å; symbol, numbered label or Z; pasted `.xyz` header accepted; atoms closer than 0.4 Å refused |
| Charge (N) | 0 | integer; every state must keep ≥ 1 electron |
| Multiplicity (N) | lowest for the electron count | `systemProperties.multiplicity_problem` |
| Multiplicity (N+p, N−q) | M−k if ≥ 1, else M+k | same check; a NOTE-style warning when \|ΔM\| ≠ k |
| p, q | 1, 1 | positive integers |
| Ions from N orbitals | off | §7 |
| Absolute `%Chk` | on | folder without spaces / non-ASCII when on |

System-specific values that must be supplied or confirmed: geometry, charge and multiplicity of N, the ions' multiplicities, `%Mem`, `%NProcShared`, output folder, formchk location if not found.

## 3. Route rules

Errors: `Opt`, `IRC`, `IRCMax`, `Scan`, `Geom=`, `Guess=`, `Units=`, `Pseudo`, `ExtraBasis`, `ChkBasis`, basis `Gen`/`GenECP`, any keyword token containing `READ`, a method containing `/` or blanks. Warnings: `Freq`; post-SCF methods (MP*, CC*, QCISD, CI*, BD, CASSCF, SAC-CI, EOM) and `TD`/`TDA`/`CIS` keywords (the analysis reads only `Total SCF Density`); diffuse basis; `R` prefix with an open-shell state; no basis.

## 4. Linked job

Order N, N+p, N−q. Each job: `[%OldChk]`, `%Chk`, `%Mem`, `%NProcShared`, route (wrapped at 72 columns), blank, title, blank, `charge multiplicity`, geometry, blank. Jobs joined by a `--Link1--` line; the file ends with a blank line. Title: `<stem with forbidden chars as spaces> Fukui <N | N plus p | N minus q> state charge <c | minus c> multiplicity <m>` — no `@ # ! - _ \`.

## 5. Conversion and compatibility

`<stem>_formchk.bat`: `cd /d "%~dp0"`, `set "FORMCHK=…"`, `set "GAUSS_EXEDIR=…"`, then `"%FORMCHK%" "<chk>" "<stem>_<tag>.fchk"` with an `errorlevel` message per file, `pause`; `%` escaped as `%%`. The explicit output name forces `.fchk` (G09W's one-argument default is `.fch`; the analysis tool accepts both).

The button reads the `%Chk` lines from the saved `.gjf`, asks before replacing existing `.fchk`, runs formchk in `asyncio.to_thread` (600 s per file, `GAUSS_EXEDIR` set), reports each file, and then calls `check_ready_for_analysis`, which applies `fukuiFunctionAnalyze.validate_states` itself (imported lazily so `fukuiInput` stays Toga-free).

## 6. Errors and tests

Every user error is a `ValueError` shown as a dialog; `prepare` reports all problems at once. Writes go through `write_files` (`OSError` → `ValueError`); existing files are replaced only after confirmation. Tests: `test_fukuiInput.py` (parsers, geometry, route, states, text layout, CRLF files, script, conversion with a fake `subprocess.run`, one real-formchk test, readiness on the phenol fixtures, siblings, no-Toga), `test_fukuiInputBuilder.py` (headless UI: file loading, every spec field from a widget, suggestions, preview/save/convert), `test_gaussianUtilities.py`, two auto-fill tests in `test_fukuiFunctionAnalyze.py`.

## 7. Found while building

- **G09 Rev. A.02 has no `%OldChk`.** With it, job 2 stops with `Error parsing percent card` (verified with both a long absolute and a short relative path). The Guess=Read option stays, off by default, with a permanent warning when switched on.
- **`charge -1` in a title** violates the no-`-` rule; negative charges are written as `minus 1`.
- **End-to-end check** on the installed G09W: water HF/STO-3G → three Normal terminations → formchk via the button path and via the `.bat` → ready → `FukuiAnalyser.run` (condensed f⁻: O 0.642, H 0.179, H 0.179).

## 8. Revision: Gaussian 16, optimization + Geom=Check

User requirements (2026-10-17): Gaussian 16 Rev. C.01 under WSL or on a cluster; job 1 optimizes the N state; the Link1 jobs are single points on N+1 and N−1 at that optimized geometry, never re-optimized; `%OldChk` from the N checkpoint plus `Geom=Check`, each ion with its own `%Chk`; charges and multiplicities per electron count, **asking rather than assuming** unspecified spin states.

User answers: formchk via a `.sh` script only (the `.bat`, the Convert button and the Windows formchk search go); the single-point mode is replaced, not kept; `Guess=Read` on the ions; no `Freq` in job 1; worked example = phenol 0 1 / −1 2 / +1 2.

| Item | Now |
|---|---|
| Job 1 | `%Chk=<stem>_N.chk`, route + `Opt` (field *Optimization (job 1)*, default `Opt`; one Opt token, TS/QST/ModRedundant/ReadFC refused), charge Q, M, explicit Cartesian coordinates |
| Jobs 2, 3 | `%OldChk=<stem>_N.chk`, `%Chk=<stem>_Np1.chk` / `_Nm1.chk`, route + `Geom=Check Guess=Read` (Guess=Read switchable), charge Q−p / Q+q, M₊ / M₋, no coordinates |
| Multiplicities | all three required (`Optional[int]`, None → "enter its spin multiplicity; it is not assumed (… allow …)"); the suggestion functions were removed; loading a file clears the ion fields |
| Files | relative names, LF line endings: `<stem>_fukui.gjf`, `<stem>_formchk.sh` (bash; `command -v formchk`, per-file missing-chk message, exit status) |
| Removed | `absolute_chk`, the core-count warning, `GUESS_READ_WARNING`, `checkpoint_pairs`, `convert_checkpoints`, `check_ready_for_analysis` (the analysis tool performs the same check on Run) |
| Kept | route rules (Opt now points to the job-1 field), title rules, `%Mem` rule, `sibling_state_files` auto-fill, `gaussianUtilities` (still used by the ESP tool's cubegen search) |

Verified on the user's WSL (Oracle Linux 9.5, `/usr/local/g16`, `ES64L-G16RevC.01`):
- phenol B3LYP/6-31G(d): 3 Normal terminations, Geom=Check read twice, no Opt steps in jobs 2–3, `bash phenol_formchk.sh` exit 0, ion coordinates identical to the N fchk (0.0 bohr), `validate_states` → p = q = 1 with no notes, `FukuiAnalyser.run` → IP 0.3014, EA −0.0812, η 0.3827 Ha (the Multiwfn reference values);
- distorted water HF/STO-3G: O–H 1.092 / 1.012 Å → 0.989 / 0.990 Å after 5 Opt steps in job 1; both ions at exactly that geometry.
