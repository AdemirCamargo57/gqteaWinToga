# Vanderbilt USPP generator (`runatom.x`) inside gqteaWinToga — design

Date: 2026-10-28. Status: **approved by the user (decisions in §0), implementation by subagent-driven development** following [the plan](../plans/2026-10-28-runatom-integration.md).

The program the request calls "run.x" is `runatom.x`, the only executable in USPP 7.3.6
(`Source/runatom.f`, `program runatom`). This document uses that name; the Windows build is `runatom.exe`.

## 0. User decisions (2026-10-28)

| # | Question | Decision |
|---|---|---|
| D1 | Integration approach | **Keep the Fortran, ship a statically linked `runatom.exe`, run it as a child process** (§3). Approved. |
| D2 | Windows Fortran toolchain (none was installed) | **MSYS2 UCRT64 `mingw-w64-ucrt-x86_64-gcc-fortran`**, installed once with `winget install MSYS2.MSYS2` (approved by the user). |
| D3 | Ship the built programs | **`runatom.exe` (and `reform.exe`, §6.1) are committed** with SHA-256 and compiler version, **and packaged with gqteaWinToga** (user requirement): the packaged build must run them with nothing else installed. |
| D4 | Licence (GPL v2 + `AA-README`'s non-commercial clause) | **No conflict**: gqteaWinToga is non-commercial research software (user). The unmodified source, `LICENSE`, `AA-README` and our patch are shipped beside the programs (GPL §3a). |
| D5 | Target code for the pseudopotential | **CPMD, formatted file** (user). Produced with USPP's own `Utils/reform.f` (§6.1) and used in CPMD with the atom-line keyword `FORMATTED` (verified in CPMD's source, §6.1). The binary `.uspp` is kept as well. |
| D6 | Ship the `Work/` input library? | **Yes, as templates** that the UI can import. |
| D7 | 7th command-line argument (`keyee < 0`) | **Not supported**; such an input is refused with an explanation. |
| D8 | How to implement | **Subagent-driven development** of the plan (user). |

Still unknown, not blocking: the download URL and checksum of `uspp-736.tar.gz` (to be recorded in `third_party/uspp-736/PROVENANCE.md` when known).

## 1. What inspection of the code established

Everything in this section was checked against the files in
`C:\Users\ajc\Downloads\uspp-736-0\uspp-736` or measured by building and running them
(gfortran 11.5 in the Oracle Linux WSL distro, scratch copies only).

### 1.1 Source

- 13 fixed-form Fortran 77 files, 16 089 lines, about 8 250 non-comment lines and about 7 300 executable statements excluding `write`/`format`. There are about 130 subroutines and functions.
- **Numerics in the code:**
  - all-electron radial solvers: `schsol`, `difsol` (Schrödinger), `dirsol` (Dirac) and `koesol` (Koelling–Harmon scalar-relativistic), with self-consistency mixing (`mixer`);
  - xc functionals: CA/PZ, Wigner, Hedin–Lundqvist, Gunnarsson–Lundqvist, B88, LYP, P86, PW91 (`aesubs.f`) and PBE (`pbe.f`);
  - pseudization: `pssubs.f` and Troullier–Martins-like norm-conserving (`trou.f`);
  - generalized eigenproblem and descreening: `scgsubs.f`;
  - linear algebra: bundled EISPACK (`eissubs.f`: `rsg`, `tred2`, `tql2`) and LINPACK (`trou.f`: `dgefa`, `dgesl`);
  - splines and integration: `numsubs.f`, `interpol.f`;
  - ghost-state check: Hamiltonian diagonalised in a spherical-Bessel basis (`bessel.f`).
- **No external libraries.** BLAS-like routines (`daxpy`, `ddot`, `idamax`) are bundled.
- **Non-standard extensions used:** `getarg`/`iargc`, `idate`, `call exit(n)`, Hollerith formats (`7h rlogd=`), `real*8`, `dfloat`. gfortran accepts all of them with `-std=legacy`.
- **Fixed capacities** (`runatom.f:20-26`, repeated as literals in `pssubs.f` (`qfplot(1000,…)`) and `trou.f` (`ndm=1000`)):

  | Limit | Value |
  |---|---|
  | Radial mesh points | ≤ 1000 (checked in `initsubs.f:38`) |
  | Atomic shells | ≤ 26 |
  | Projectors (β) | ≤ 10 |
  | Valence states | ≤ 5 |
  | lmax | ≤ 3 |
  | References per l | ≤ 4 |
  | Polynomial terms (`nqf`/`npf`) | ≤ 20 |

### 1.2 Build (verified)

- The shipped Makefile uses `FC = f77`. With gfortran the unmodified source **fails to link**: `runatom.f:86` declares `external idate`, giving `undefined reference to 'idate_'`.
- Commenting out that one line (gfortran then uses its intrinsic `IDATE(VALUES)`: day, month, year) builds `runatom.x`. This is the **only** source change needed. It affects only the date printed in the `.out` header and in the `.uspp` file's `idmy` record.
- With that patch, the reference example `Work/000-Ex/000-Ex-ca-sp-test` reproduces the shipped `Sample_output`. The eigenvalue tables are identical to the printed 6 decimals (e.g. 300: −4.576045 / −4.575601). The remaining differences are 10⁻¹²-level convergence noise, the platform sensitivity documented in `Source/BENIGN_BUG`.
- `Source/COMPILER_DETAILS` says to keep the compiler's default (static vs. automatic) array allocation.

### 1.3 Runtime interface (from `runatom.f`, `rwsubs.f`, `aesubs.f`)

**Command line:** exactly six positional arguments, each `character*40`:

```
runatom input.adat output.out file.ae file.atwf file.logd file.uspp|dummy
```

A seventh argument (reference-wavefunction file) is optional (D7). **Arguments longer than 40 characters are silently truncated** (`getarg` into `character*40`), so absolute Windows paths cannot be passed. The program must run with its working directory set to the job folder and get short relative names.

**Files:**

| File | Mode | Written / read by |
|---|---|---|
| `.adat` (input) | formatted, **fixed columns** | read with `5i5`, `3f10.5`, … (§1.4) |
| `.out` | formatted text, 1 700–1 800 lines | unit 6 is opened on the 2nd argument |
| `.ae` | unformatted | written by the AE run, read by the generation/test runs |
| `.uspp` | unformatted, layout = `Utils/read_uspp.f` | written by generation (`ifpsp=2`), read by tests (`ifpsp=1`) |
| `.logd` | unformatted | log-derivative tables |
| `.atwf` | formatted, `(r, f)` pairs separated by `&` | written only when `ifplw = 1` or `2` |

**Progress:** the program writes stage messages and one line per SCF cycle to **stderr** (unit 0, unbuffered):
- `beginning the all electron calculation`
- `completed self-consistent cycle  NN delta = D`
- `all electron calculation completed`
- `begin pseudopotential test`
- `fourier transforming the qfunctions`
- `descreening the potential`
- `writing the pseudopotential to file`
- `beginning search for possible ghost states`
- `eigensolution in bessel basis with cutoff XX.XX`

So progress can be streamed from a pipe.

**Exit status:**
- 141 sites `call exit(n)`. All but 7 are `exit(1)`, each preceded by `***error in subroutine <name>` plus one or two explanatory lines **in the `.out` file**, not on stderr.
- About 8 bare `stop` statements (exit status **0** in gfortran), with or without a message on stderr (e.g. `SPLINT: xfit not properly ordered`, `rtflsp max iter execeeded`, `stop 99`).
- Fortran runtime errors (exit status 2) report on stderr, e.g. `At line 165 of file rwsubs.f (unit = 10…) Fortran runtime error: Bad value during floating point read`.
- **A zero exit status therefore does not prove success.** Success is "status 0 and the stage's expected output exists and no `***error` line is in the `.out`".

**Run time (measured, one process per stage):**

| Case | AE run | Generation run |
|---|---|---|
| Ti (example) | 0.04 s | 0.41 s |
| O PBE | 0.04 s | 0.19 s |
| Pt PW91 | 0.15 s | 0.36 s |
| Ba PW91 | 0.20 s | 0.68 s |

A full job (AE + generation + 3 × (AE + test)) is about 1–3 s.

**Output landmarks in `.out` (for parsing):**
- `***error in subroutine`
- `** mild warning: max q eigval somewhat large **` (q-matrix eigenvalue > 5; "severe" > 20; program stops > 50)
- the boxed `pseudopotential report`
- `comparison of all-electron and pseudo eigenvalues (Ry)` with `nlm all-elec pseudo diff` rows
- `eigenvalues of US hamiltonian` with `l= 0 …` rows (last block = highest Bessel cutoff, the ghost check)
- `total energy =` (test runs)
- `chart:` blocks (ASCII plots)

### 1.4 Input format (from `rwsubs.f` `readin`, lines 130–673)

- **Every line is read with a fixed-column format, with blanks ignored** (gfortran default `BLANK='NULL'`). Verified with gfortran:
  - `1 0 0 0 3` read as `6i5` gives `ifae=100, ifpsp=3`.
  - `2` read as `f10.5` gives `2.0e-5` (an `F` field with no decimal point is scaled by 10⁻ᵈ).
- Only the columns covered by the format are read, so trailing labels are harmless **only if they start after the format width**. `6i5` reads 30 columns although only 5 values are written.
- Line order and conditional lines:
  - **All-electron file** (`ifae=1`): control (`6i5`), log-derivative mesh (`3f10.5,i5`), convergence (`2e10.1,f10.5,i5`), title (`a20`), then `z,xion,exfact` (`f5.0,2f10.5`), `rmax,aasf,bbsf` (`3f10.5`), `ncspvs,irel` (`2i5`) and `ncspvs` lines of `nnlz,wwnl,ee` (`i4,f7.3,f14.6`).
  - **Generation file** (`ifae=0, ifpsp=2`): the same 4 lines, then:
    - `ncores,nvales,nang` (`3i5`)
    - `besrmax,besemin,besemax,besde` (`4f10.5`)
    - `keyps,ifpcor,rinner1,rpcor` (`2i5,2f10.5`)
    - **if `rinner1 = 0`:** `2·nang−1` values of `rinner` (`5f10.5`)
    - `nbeta,rcloc` (`i5,f10.5`)
    - `nang` values of `rc` (`4f10.5`)
    - `nbeta` lines of `lll,keyee,eeread,iptype` (`2i5,f10.5,i5`)
    - **if any `iptype = 2`:** `npf,ptryc` (`i5,f10.5`)
    - `lloc,keyee,eloc,iploctype` (`2i5,f10.5,i5`)
    - `ifqopt,nqf,qtryc,nfix` (`2i5,f10.5,i5`)
    - `nfix` lines of `ibfix,rfix,qfix` (`i5,2f10.5`)
  - **Test file** (`ifae=0, ifpsp=1`): the 4 control lines, `ncores,nvales,nang`, and the Bessel line.
- **Checks the program makes** (each ends the run with `***error`):
  - |z − Σwwnl − xion| ≤ 1e-3; negative `wwnl` means "skip the orbital" and counts as 0.
  - `ncspvs` ≤ 26; `irel` ∈ {0, 1, 2}; no pseudopotential work with `irel = 1`.
  - `ncores + nvales = ncspvs`; `nvales` ≤ 5; `nang` ≤ 4.
  - `keyps` = 3 is the only implemented value.
  - every `rinner` ≤ `rlogd`.
  - `nbeta` ≤ 10; β lines in non-descending l; ≤ 4 β per l; `nang = lmax + 1`.
  - `keyee` ≤ `nvales`; `iptype` ∈ 0..4.
  - `npf` > 4 and ≤ 20 when any `iptype = 2`.
  - `iploctype = 2` needs an `iptype = 2` β when `keyee > 0`.
  - `ifqopt` ∈ {<0, 0, 1, 2, 3}; `ifqopt = 0` wants `nqf = 3` (warning only); `ifqopt ≥ 1` needs 3 < `nqf` ≤ 20.
  - `nfix = 0` unless `ifqopt ≥ 2`; 0 ≤ `nfix` ≤ `nbeta`.
  - **any GGA (`exfact ≥ 1`) requires `ifqopt = 3`.** That includes `ifqopt = −1`, so `Doc/TUTORIAL`'s advice to start with `ifqopt=-1` fails for PBE/PW91.
  - each `rfix` ≤ every `rinner`; `ibfix` strictly increasing, ≥ 0 and ≤ `nbeta`.
- The 419 library `.adat` files outside `History/` were parsed with an emulation of these rules: 418 read cleanly. The 2 files in `082-Pb/History/` are pre-7.2 inputs for an older program version. `023-V-gpw-nsp-campos/V_ps.adat` line 2 fails (`Bad value during floating point read`, reproduced with the real program).

### 1.5 The existing `runatom_input_gen.py` (538 lines, no tests)

It writes the text the user types **verbatim** into the fixed-column file, so it silently produces misread input (§1.4). It also:
- always writes the `npf,ptryc` line (misaligns every later line when no β has `iptype = 2`);
- supports one β line regardless of `nbeta`, and has no `rpcor`, `nfix` fix lines, per-L `rinner` line or test inputs;
- hard-codes the title as `hydrogen <cfg>`;
- labels the z line `f7.2` (the program reads `f5.0`);
- copies the outdated `ifprt` range from `INPUT_AE`;
- writes a Makefile that refers to `a7.3.3` and has no test targets.

The module name and `RunatomGenerator` class are referenced by `toolLoader.TOOL_TABS` and `_frozen_build_imports()`.

## 2. Approaches compared

| | **A. Fortran kept, `runatom.exe` as a child process** | A′. Fortran as an in-process library (f2py / ctypes DLL) | **B. Port to Python** (numpy, optional numba) |
|---|---|---|---|
| **Performance** | Native: 0.04–0.7 s per stage plus ~20–50 ms process start. Already far below anything a user notices. | Same numerics, no process start. | ~7 300 Fortran statements, dominated by per-mesh-point radial ODE recurrences inside eigenvalue bisection and SCF loops. These cannot be vectorised, so pure numpy/Python is expected to be 50–500× slower (seconds to minutes). numba could recover most of it, but numba is **optional** in this project (CLAUDE.md: every numba path needs a numpy fallback giving the same numbers), so the fallback would be the slow path. |
| **Numerical reliability** | Identical algorithm and code paths to the reference. Verified: reproduces `Sample_output` to printed precision. | Same as A. | A re-implementation of ~130 routines, including EISPACK `rsg`/`tql2` (which `scipy.linalg.eigh` would replace with LAPACK: different eigenvector signs/order) and spline/derivative helpers. `Source/BENIGN_BUG` documents that the β transformation matrix `ttt` is "supersensitive" to convergence level and platform, so a port **cannot** be compared array-for-array with the original. It can only be compared through invariants, and every subtle divergence becomes a research question. |
| **Effort** | Small: one-line patch, build script, Python wrapper (input writer/validator, runner, output parser, UI). About 2–3 weeks including tests. | Moderate, and blocked by the issues in the next row. | Very large: several person-months for a domain expert, plus a validation campaign over the 78-case library. |
| **Robustness inside a GUI** | Crashes, `call exit(1)` and runaway loops are isolated in the child; cancel = kill the process. | **Unworkable as is:** 141 `call exit` and ~8 `stop` statements would **terminate the whole GUI process**. The code keeps state in `SAVE`/static arrays (not re-entrant) and cannot be cancelled. Making it a library means rewriting all error exits. | In-process Python: cancellable by flag checks, errors are exceptions. |
| **Maintenance** | Upstream source untouched apart from one documented patch. The Fortran is frozen (last release 2006, "no active support"), so there is nothing to track. | Patched fork of every error path. | A second implementation to own forever, with no upstream to compare against except the Fortran itself. |
| **Windows deployment** | One static exe (no `libgfortran`/`libquadmath`/`libgcc` DLLs) added to the PyInstaller bundle as a binary, plus GPL source. Unsigned exe may trigger antivirus heuristics (see risks). | DLL + compiler runtime; f2py on Windows needs meson + a Fortran compiler at build time and an ABI match with the Python/numpy build. | Pure Python: simplest to ship, but numba adds ~60 MB to the bundle and JIT warm-up on first use. |

## 3. Recommendation

**Approach A: keep the Fortran unchanged except for the `idate` patch; build a statically linked `runatom.exe` with MSYS2 UCRT64 gfortran; run it as a child process from a worker thread.** Approach B is not justified.

1. **The bottleneck B would remove does not exist.** Measured runs take 0.04–0.7 s, and the user's wait is dominated by their own trial-and-error edits, not by computation.
2. **Exactness is the requirement that matters.** Pseudopotentials are judged against the original program's outputs (log derivatives, eigenvalue tables, ghost checks). Only the original code gives those outputs by construction. The documented `ttt` instability means a port would *look* different even when correct, making validation itself ambiguous.
3. **A child process is the only way to run this code safely inside a GUI:** 149 process-terminating exits rule out an in-process library without rewriting the error handling.
4. **The changes are small and reversible.** The integration is pure Python around an unmodified upstream program; if Approach B were ever wanted, the same input writer, output parser and UI would be reused unchanged, with only the runner swapped.

## 4. Architecture

New, all flat modules (the project convention). None of the three numerical modules imports Toga (subprocess tests):

- `runatomInput.py`
  - dataclasses `Controls`, `Shell`, `AEInput`, `PSInput`, `Beta`, `QFix`, `TestInput`, `TestConfig` and `CaseInputs`;
  - the fixed-column writer `format_adat`;
  - a reader `parse_adat` that emulates gfortran's fixed-format input (blank = null, implied decimal, short-record padding);
  - `parse_line_values` (space-separated UI text to typed values);
  - `check_ae`/`check_ps`/`check_case`, which mirror §1.4 one rule per message and run **before** anything is executed, on the values as runatom will read them (after a format/parse round trip).
- `runatomRunner.py`
  - `find_runatom(explicit)`;
  - `Stage`/`RunatomJob` (the Makefile's dependency chain as data: AE → generation → per test: AE → test);
  - `run_job(job, progress, cancel)`: each stage is a `subprocess.Popen` in the job folder with relative names ≤ 40 characters, stderr streamed and parsed into progress, and a timeout;
  - error classification into `RunatomError(stage, kind, message, details)` / `RunatomCancelled`.
- `runatomOutput.py`
  - `parse_out(text)` → `OutReport`: errors, warnings, report block, eigenvalue comparison, Bessel-basis eigenvalues per l, total energy;
  - `ghost_candidates`;
  - `read_uspp(path)` → `UsppData`: the record layout of `Utils/read_uspp.f`, with 4-byte little-endian record markers validated.
- `runatom_input_gen.py` (existing module and class names kept) becomes the window:
  - Job-folder header (with the detected `runatom.exe`);
  - tabs **All-electron / Pseudopotential / Tests / Results / Help**;
  - a run bar with **Validate**, **Save inputs**, **Run**, **Cancel**, **Open folder**, an `ActivityIndicator` and a status line.
- `help.py` gains `HelpGqteaWin.help_runatom` (replacing the module's two help constants).

Chain, one-way: `runatom_input_gen → {runatomRunner → runatomInput, runatomOutput, help}`.

Vendored upstream under `third_party/uspp-736/`:
- `Source/*.f` unmodified, plus `LICENSE`, `README`, `BENIGN_BUG`, `COMPILER_DETAILS`;
- `Doc/`, `Utils/read_uspp.f`, `Work/` and `AA-*`;
- `patches/0001-gfortran-intrinsic-idate.patch`, `build_runatom.ps1` and `PROVENANCE.md`;
- `bin/runatom.exe` with `bin/runatom.exe.sha256`.

The build script copies `Source/` to a temporary folder, applies the patch and compiles. The vendored source stays byte-identical to upstream.

## 5. Execution design

- **Job folder:**
  - The user chooses it; default: a new folder `<symbol>_uspp` beside the last saved input.
  - Inputs are written there as `<sym>_ae_<cfg>.adat`, `<sym>_ps.adat` and `<sym>_test_<cft>.adat` (the library's naming, so the library Makefiles still work), with `newline="\n"`.
  - The runner refuses names over 40 characters or containing anything outside `[A-Za-z0-9._-]`.
  - Existing outputs from a previous run in the folder are deleted for the stages being rerun, so a stale `.uspp` can never be mistaken for a new one.
- **Stages:** each stage is one process, `[exe, inp, out, ae, atwf, logd, uspp_or_dummy]` with `cwd = job folder`, mirroring `Work/000-Ex/000-Ex-ca-sp-test/Makefile` exactly:
  1. AE reference;
  2. generation;
  3. for each test: AE test, then test run.
- **Threading:**
  - `workflow` awaits `asyncio.to_thread(run_job, …)`. The worker thread owns the process.
  - A reader thread pushes stderr lines into a queue; the worker polls `proc.poll()` every 50 ms, drains the queue into `progress(dict)` and checks the cancel `threading.Event` and the stage deadline.
  - `progress` is marshalled with `loop.call_soon_threadsafe` and thinned to 0.1 s (the `meanResidenceTime_classical` pattern).
  - Toga's loop is a `ProactorEventLoop`, so `asyncio.create_subprocess_exec` would also work. The thread pattern is kept because it is the project's tested pattern and is testable without an event loop.
- **Cancellation:**
  - **Cancel** sets the event; the worker `kill()`s the process (TerminateProcess on Windows), waits for it and raises `RunatomCancelled`.
  - Closing the window cancels too (`_closed` flag, as in the classical MRT tool).
  - Outputs of the killed stage are deleted.
- **Timeout:** per stage, default 120 s (more than 100× the slowest measured stage), configurable for experiments with `maxit` raised.
- **Windows specifics:**
  - `creationflags=CREATE_NO_WINDOW` (no console flash in the frozen GUI);
  - stdout to `DEVNULL` (the program writes its report to the `.out` file);
  - stderr decoded as latin-1 (never fails).
- **Error classification** (first match):

  | Kind | Condition | Message source |
  |---|---|---|
  | `cancelled` | the user pressed Cancel | — |
  | `timeout` | the stage deadline passed | — |
  | `input_format` | stderr has `Fortran runtime error` and `unit = 10` | names the input file and the reported line |
  | `program_error` | the `.out` contains `***error in subroutine X` | that line and the next two |
  | `negative_density` | stderr has `stopping in vhxc: negative densities` | hint: reduce `rinner` or add `nfix` constraints (`Doc/TUTORIAL`) |
  | `not_converged` | `.out` says `potential could not be converged` | — |
  | `runtime_error` | any other non-zero status | stderr tail |
  | `missing_output` | status 0 but an expected file is absent | — |

## 6. Output processing

After each stage, the runner calls `parse_out`. The Results tab shows:
- per stage: status, time, warnings;
- the eigenvalue comparison table, with each |diff| > 1 mRy flagged (the `Doc/TUTORIAL` threshold);
- the q-matrix warning lines;
- the ghost check;
- the `.uspp` summary from `read_uspp`: zv, nbeta, kkbeta, mesh and rc.

**Ghost check** (`ghost_candidates`): for each l present in the eigenvalue comparison, take the k all-electron valence eigenvalues and the last-cutoff US-Hamiltonian eigenvalues u. The check flags l when more than k values of u lie below `max(AE_l) + 0.1 Ry`. This is the `Doc/TUTORIAL` criterion: an extra bound state below the expected ones. On Ti it gives no ghost for l = 0, 1, 2.

The full `.out` stays one click away (**Open .out**), because the ASCII charts are the program's main diagnostic and are not re-plotted in this plan.

### 6.1 Formatted pseudopotential for CPMD (D5)

**What CPMD reads** (checked in CPMD's source, `CPMD-code/CPMD` on GitHub):

- `src/ratom_utils.mod.F90`: an atom line containing `FORMAT…` sets `tvan` (a Vanderbilt potential, formatted). `BINARY` sets `tvan` and `tbin`; both together are an error. The species line is therefore `*<file> FORMATTED`.
- `src/readvan_utils.mod.F90`, `readvan`: the formatted branch reads the records of `Utils/read_uspp.f`'s layout, in the same order.
  - Every `READ` is list-directed (`READ(ifn,*)`), except the title record, `'(A20,3F15.9/2I5,1PE19.11)'`.
  - `iptype`, `npf` and `ptryc` are read in one statement.
  - `nchi` (`natwf`) is read for version ≥ 7.

That is exactly what USPP's own `Utils/reform.f` writes with `itype = 1`:
- formats `6i5`, `(a20,3f15.9/2i5,1pe19.11)`, `(i5,2f15.9)`, `(2i5,f15.9)`, `(2i5,f9.5,2i5,f9.5)`, `(i5)`, `(1p4e19.11)`, `(i5/(1p4e19.11))`;
- `nchi` list-directed.

So the formatted file is produced by **`reform.exe`, built from the unmodified `Utils/reform.f`** with the same toolchain and flags as `runatom.exe` (no patch needed: verified, it compiles and converts the Ti example). It is **not** produced by a Python re-implementation.

**Running reform:** the runner calls it after a successful generation stage, in the job folder, with stdin `1\n<symbol>_ps.uspp\n<symbol>_ps_formatted.uspp\n`.
- File names are read as `a36`, so they must be at most 36 characters.
- The output is opened with `status='new'`: reform stops with "File exists", so the runner deletes an old one first.
- Success means exit status 0, the output file present, and `pseudopotential has been written out` on stdout.

**Precision and exponents:** `1p4e19.11` keeps 12 significant digits (relative rounding ≤ 5·10⁻¹²). For a 3-digit exponent Fortran drops the `E` (`1.23456789012-100`); gfortran's list-directed read accepts that form (verified), which is what CPMD uses.

**In the window:**
- the Results tab names the file and prints the CPMD line `*<symbol>_ps_formatted.uspp FORMATTED` (the file goes to CPMD's `PP_LIBRARY_PATH` or the run directory);
- a missing `reform.exe` gives a note, not a failure — the binary `.uspp` is still produced.

**Verification:** `runatomOutput.read_formatted_uspp` reads the file in `readvan`'s statement order, with list-directed semantics: each READ starts on a new line and may span lines. It must equal `read_uspp` of the binary file to 1e-11 relative. That is the CPMD-compatibility check, without CPMD installed.

## 7. Validation

- **Reference implementation:** the unmodified source plus the `idate` patch, built with Linux gfortran in WSL. It is the build verified against `Sample_output` (§1.2).
- `tests/data/runatom/make_reference.sh` (dev-time only, never run by the suite) runs 8 library cases and stores their `.out` and `.uspp` files. The cases are chosen so that every reader branch and functional family is exercised:

  | Case | What it covers |
  |---|---|
  | `000-Ex-ca-sp-test` (Ti) | LDA, semicore, 2 β per channel, `iptype=2`, test |
  | `001-H-gpbe--bm` | PBE, 1 β, `irel=0` |
  | `008-O-ca--vgrp` | LDA, test |
  | `026-Fe-gpw-n-campos` | PW91, NLCC, per-L `rinner`, `iptype` 3/4, 3 tests |
  | `029-Cu-gpbe-n-bm` | PBE, NLCC, fractional occupations |
  | `079-Au-ca--ts` | `ifqopt=0`, `lloc=−1`, per-L `rinner` |
  | `082-Pb-ca-d-vgrp` | `nfix=5` constraints, test |
  | `056-Ba-gpw-nsp-campos` | 6th row, `lloc=3`, 3 tests |

- **Input round trip** (no exe needed): for all 418 readable library files, `parse_adat(format_adat(parse_adat(f))) == parse_adat(f)`; `023-V` must be reported as unreadable with the line number.
- **Exe equivalence** (skipped when no `runatom.exe` is found, like the cubegen tests):
  - re-formatting a library input with `format_adat` and running it gives a `.uspp` **byte-identical** to running the original input with the same exe (except the `idmy` date record, which is masked).
  - `runatom.exe` vs. the Linux reference:
    - every eigenvalue-comparison and total-energy number agrees within 2e-6 Ry (one unit in the printed last place);
    - Bessel-basis eigenvalues agree within 1e-5 Ry;
    - `.uspp` scalar fields agree exactly; `r`/`rab` within 1e-12 relative; `vloc`, `vloc0` and `rsatom` within 1e-7 relative or 1e-10 absolute.
    - `beta`, `ddd`, `qqq` and `qfunc` are compared up to the sign of each β within 1e-6 relative, **except** the s channel of Pb, which `BENIGN_BUG` documents as unstable. There only the invariants are compared: the eigenvalues of the D and Q blocks per l, and the Bessel eigenvalues.
  - Rule (as for the QTAIM/NTO references): a disagreement is investigated in the build, never fixed by loosening a tolerance.
- **Historical check:** `Sample_output` (2006) eigenvalue tables equal to printed precision.
- **Formatted files for CPMD:**
  - `reform.exe` applied to each *reference* binary `.uspp` must reproduce the Linux `reform.x` output of that same file **byte for byte**, except line 1 (date). Formatted output is correctly rounded, so it is deterministic.
  - `read_formatted_uspp(formatted)` must match `read_uspp(binary)` of our own runs to 1e-11 relative.
- **Benchmarks:** `benchmarks/bench_runatom.py` times every library case end to end (process start included) and reports per-stage times. Acceptance:
  - each case < 3 s on the dev machine;
  - Toga event-loop lag during a run ≤ 50 ms p95, measured with the `toolLoader` method.

## 8. Packaging

- auto-py-to-exe / PyInstaller:
  - add `third_party/uspp-736/bin/runatom.exe` and `bin/reform.exe` as **binaries** with destination `uspp/` (D3: they ship with gqteaWinToga);
  - add `LICENSE`, `AA-README`, `Source/`, `patches/` and `Work/` as **data** under `uspp/`.
- `find_runatom` looks, in order, at:
  1. the explicit path field;
  2. `%GQTEA_RUNATOM%`;
  3. `<sys._MEIPASS>/uspp/runatom.exe` when `sys.frozen`;
  4. `third_party/uspp-736/bin/runatom.exe` beside the module;
  5. `runatom.exe` / `runatom.x` on `PATH`.
- The exe must be **statically linked**. The build step fails if `objdump -p runatom.exe` lists any DLL other than `KERNEL32.dll` and the `api-ms-win-crt-*` UCRT set (part of Windows 10+).
- `readme.md` gets a packaging subsection and `USER_MANUAL.md` a rewritten tool section plus a "Generating an ultrasoft pseudopotential" workflow.

## 9. Risks and compatibility issues

| Risk | Mitigation |
|---|---|
| 40-character argument truncation | Run in the job folder with validated short relative names; a test pins the limit. Unicode or space characters in the folder path are then harmless (CreateProcess handles the cwd). |
| Exit status 0 on failure (`stop`) | Success requires the expected outputs plus no `***error` line; tests cover a fake program that exits 0 without output. |
| Newer gfortran (MSYS2 ships 14/15) changes optimisation or rounding | Build at `-O1` (the Makefile's advice); validate against the reference with the tolerances above. If they fail, try `-O0` and `-ffp-contract=off` before anything else, and record the outcome. |
| `-std=legacy` warnings or new errors in newer gfortran (e.g. argument mismatch in old F77 calls) | Add `-fallow-argument-mismatch` only if the compiler requires it; record it in the build script. |
| CRLF input | Always write LF; the reader accepts both. |
| Unsigned exe flagged by antivirus / SmartScreen | Static, no packer; document it; code signing is a release decision for the user. |
| Non-converging transferability tests (reported in `Doc/TUTORIAL` by B. Meyer) | Timeout + Cancel, plus the hint "lower `damp` to 0.1–0.2". |
| Library defects found earlier (`023-V` unreadable, all-zero Nb `.uspp`, `C`/`Cr` named `-n-` without NLCC, `Cd` misnamed) | Imported as is; the importer reports the V error; the defects are listed in `third_party/uspp-736/PROVENANCE.md`. |
| `Pot/`/`Bin/` Makefile workflow on Windows | Not used. The Makefile generator option stays only for users who run on Linux/WSL, and is fixed (test targets, version). |
| Licence (D4) | User decision before any release that includes the exe. |

## 10. Out of scope

- Re-plotting the `.out` ASCII charts, `.logd` and `.atwf` in the plot viewer (a natural follow-up: the readers are simple).
- UPF export (Quantum ESPRESSO). The CPMD formatted file is in scope (§6.1).
- `keyee < 0` reference-wavefunction runs (D7).
- Porting the Fortran.
- Dirac (`irel = 1`) generation, which the program itself refuses.
