# Vanderbilt USPP Generator (runatom) Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let a user prepare, check and run the Vanderbilt ultrasoft-pseudopotential generator `runatom` (USPP 7.3.6) from **Inputs > Vanderbilt USPP Generator** in gqteaWinToga on Windows. The tool covers the all-electron run, pseudopotential generation and transferability tests, plus a results summary and a **formatted pseudopotential file for CPMD** (written by USPP's own `reform`). `runatom.exe` and `reform.exe` are packaged with gqteaWinToga. Running an external program by hand is no longer needed.

**Architecture:** The Fortran stays as it is (one documented one-line patch) and is built once into a statically linked `runatom.exe` with MSYS2 UCRT64 gfortran. Three new Toga-free modules sit around it:
- `runatomInput.py` writes and reads the fixed-column `.adat` files exactly as `rwsubs.f` reads them, and checks them with the program's own rules before anything runs;
- `runatomRunner.py` runs each stage as a child process in the job folder from a worker thread, with streamed progress, cancel and timeout;
- `runatomOutput.py` parses the `.out` report and reads the binary `.uspp` file.

`runatom_input_gen.py` is rewritten as the window.

**Tech Stack:** Python 3.13, numpy, Toga 0.5 (WinForms), pytest. MSYS2 UCRT64 gfortran (Windows build) and gfortran 11.5 in the Oracle Linux WSL distro (reference build). No numba.

**Spec:** [docs/superpowers/specs/2026-10-28-runatom-integration-design.md](../specs/2026-10-28-runatom-integration-design.md). Read §1 (what the code does), §5 (execution) and §7 (validation) before starting.

## Global Constraints

- Run everything from `venv/src/`; tests: `python -m pytest tests/ -q`. Baseline failures to ignore: `test_meanResidenceTime.py::TestAdvancedOptionsExposed::test_every_run_parameter_is_supplied_by_the_ui` and three table-styling tests in `test_plotter.py`.
- The user's "run.x" is `runatom.x`; the Windows build is named `runatom.exe`.
- The vendored upstream files under `third_party/uspp-736/` are **never edited**. The only source change (comment out `      external idate` in `runatom.f`) is applied to a temporary copy at build time.
- `runatomInput.py`, `runatomOutput.py` and `runatomRunner.py` must not import Toga (a subprocess test in each test file).
- Values typed by a user are **never** copied into an `.adat` file as text. They are parsed and re-written column-exactly, always with a decimal point in F/E fields.
- `.adat` files are written as ASCII with LF line endings (`newline="\n"`).
- runatom's file-name arguments are `character*40`. The program is always started with `cwd` = the job folder and relative names of at most 40 characters.
- A stage succeeds only when all of these hold: exit status 0, its completion marker is on stderr, and its expected files exist. A `stop` exits with status 0 even on failure.
- Do not delete or modify user files outside the job folder. In the job folder, delete only the outputs of the stage about to run.
- Commits: only once the user has authorised commits; otherwise leave the work uncommitted and report it.
- The user's decisions D1–D8 are recorded in spec §0:
  - non-commercial research software, so the GPL/`AA-README` terms are fine;
  - MSYS2 approved (installed 2026-10-28 with `winget install MSYS2.MSYS2`);
  - **`runatom.exe` and `reform.exe` must be packaged with gqteaWinToga**;
  - **CPMD formatted `.uspp`** output;
  - subagent-driven execution.
- The CPMD formatted file is written by `reform.exe` (unmodified `Utils/reform.f`), never by Python. The Python reader `read_formatted_uspp` exists only to verify it, in CPMD `readvan`'s statement order.

## File Structure

| File | Responsibility |
|---|---|
| `third_party/uspp-736/**` | Unmodified USPP 7.3.6 tree (Source, Doc, Utils, Work, AA-*), `.gitattributes`, `PROVENANCE.md`, `SOURCE_SHA256.txt`, `patches/0001-gfortran-intrinsic-idate.patch`, `build_runatom.ps1`, `bin/runatom.exe` and `bin/reform.exe` (+ `.sha256`, `BUILDINFO.txt`) |
| `runatomInput.py` | Fixed-format field I/O, input data model, `.adat` writer/reader, library case folders, UI field mapping, pre-run checks |
| `runatomOutput.py` | `.out` parsing, ghost check, binary `.uspp` reader, formatted (CPMD) `.uspp` reader for verification, results summary text |
| `runatomRunner.py` | Locating the programs and the library, stage planning, child-process execution, progress, cancel, timeout, error classification, the `reform` conversion for CPMD |
| `runatom_input_gen.py` | The Toga window (rewritten; module and `RunatomGenerator` names kept for `toolLoader`) |
| `help.py` | `HelpGqteaWin.help_runatom` |
| `toolLoader.py` | Button label only |
| `tests/runatom_fixtures.py` | Paths and readers for the reference data and the library |
| `tests/fake_runatom.py`, `tests/fake_reform.py` | Stand-in programs for runner/UI tests |
| `tests/data/runatom/` | `make_reference.sh`, `MANIFEST.txt`, 8 case folders of gzipped Linux-reference `.out`/`.uspp`/`_formatted.uspp` |
| `tests/test_runatomVendor.py`, `test_runatomReferenceData.py`, `test_runatomInput.py`, `test_runatomOutput.py`, `test_runatomRunner.py`, `test_runatomFormatted.py`, `test_runatomExe.py`, `test_runatom_input_gen.py` | Tests |
| `benchmarks/bench_runatom.py` | Library timing and event-loop lag |
| `readme.md`, `USER_MANUAL.md`, `CLAUDE.md` | Packaging, user docs, project notes |

---

### Task 1: Vendor the unmodified USPP 7.3.6 tree

**Files:**
- Create: `third_party/uspp-736/` (copy of `C:\Users\ajc\Downloads\uspp-736-0\uspp-736`)
- Create: `third_party/uspp-736/.gitattributes`, `PROVENANCE.md`, `SOURCE_SHA256.txt`, `patches/0001-gfortran-intrinsic-idate.patch`
- Test: `tests/test_runatomVendor.py`

**Interfaces:**
- Produces: `third_party/uspp-736/Work/<element>/<case>/` (the library used by every later task) and `third_party/uspp-736/Source/*.f` (the build input of Task 2).

- [ ] **Step 1: Copy the tree byte for byte**

The tree must be copied byte for byte, before git can touch line endings:

```powershell
python -I -c "import shutil,sys; shutil.copytree(sys.argv[1], sys.argv[2])" "C:\Users\ajc\Downloads\uspp-736-0\uspp-736" third_party\uspp-736
```

Then create `third_party/uspp-736/.gitattributes`:

```
# Upstream USPP files: keep every byte (fixed-column Fortran input, LF endings, binary data).
* -text
```

- [ ] **Step 2: Write the patch record**

`third_party/uspp-736/patches/0001-gfortran-intrinsic-idate.patch` is the change that `build_runatom.ps1` and `make_reference.sh` apply to a temporary copy:

```diff
--- a/Source/runatom.f
+++ b/Source/runatom.f
@@ -83,7 +83,7 @@
       character*40 flname(6)
       character*40 filewfref
 c
-      external idate
+c     external idate
 c
 c=======================================================================
 c     c h a n g e s  m a d e  t o  c r e a t e  v e r s i o n  3.0.0
```

- [ ] **Step 3: Write the provenance note**

Create `third_party/uspp-736/PROVENANCE.md`:

```markdown
# USPP 7.3.6 in gqteaWinToga

Upstream: Vanderbilt ultrasoft pseudopotential package uspp-736 (D. Vanderbilt,
Rutgers University, December 2006). Licence: GNU GPL v2 (`Source/LICENSE`),
with the additional condition in `AA-README` (no profit-making use, no inclusion
in commercial code without arrangement). Download URL and archive checksum:
**not recorded yet — add them when known.**

The files outside `patches/`, `bin/`, `build_runatom.ps1`, `PROVENANCE.md`,
`SOURCE_SHA256.txt` and `.gitattributes` are unmodified upstream files
(`SOURCE_SHA256.txt`, checked by tests/test_runatomVendor.py).

## Change made to build runatom.exe

`patches/0001-gfortran-intrinsic-idate.patch` comments out `external idate` in
`Source/runatom.f`. With gfortran the declaration hides the IDATE intrinsic and
the link fails (`undefined reference to 'idate_'`). The change only affects the
date printed in the .out header and stored in the .uspp file. It is applied to a
temporary copy by `build_runatom.ps1`; the file in `Source/` stays unmodified.

## Known defects of the upstream input library (left as shipped)

- `Work/023-V/023-V-gpw-nsp-campos/V_ps.adat` line 2: the minus sign of emin
  falls in the rlogd field; runatom stops with "Bad value during floating point read".
- `Work/041-Nb/041-Nb-gpw-nsp-campos/041-Nb-gpw-nsp-campos_us_pw91_13elec` (and
  the copy in the .tar next to it): 10240 bytes of zeros, not a pseudopotential.
- `006-C-gpw-n-campos` and `024-Cr-gpw-nsp-campos` have ifpcor = 0 although their
  names and READMEs claim a nonlinear core correction; `048-Cd-gpw-campos` has
  ifpcor = 1 but no flags field in its name.
- `039-Y-ca-sp-vgrp` is a retracted entry (README only).
- `Doc/TUTORIAL` suggests starting with ifqopt = -1, which runatom refuses for
  every GGA functional (rwsubs.f: ifqopt < 3 with exfact >= 1 is an error).
```

- [ ] **Step 4: Write the source manifest**

Run from `venv/src`:

```powershell
python -I -c "import hashlib,pathlib; r=pathlib.Path('third_party/uspp-736'); roots=['AA-QUICK-START','AA-README','AA-SUBMISSIONS','Makefile','Doc','Pot','Source','Utils','Work']; files=sorted(p for n in roots for p in ([r/n] if (r/n).is_file() else (r/n).rglob('*')) if p.is_file()); (r/'SOURCE_SHA256.txt').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(r).as_posix()}\n' for p in files), newline='\n')"
```

- [ ] **Step 5: Write the failing test**

Create `tests/test_runatomVendor.py`:

```python
"""The vendored USPP 7.3.6 tree is unmodified and the build patch still applies."""

import hashlib
from pathlib import Path

import pytest

VENDOR = Path(__file__).resolve().parents[1] / "third_party" / "uspp-736"
IDATE_LINE = "      external idate"


def _manifest():
    rows = (VENDOR / "SOURCE_SHA256.txt").read_text().splitlines()
    return [row.split("  ", 1) for row in rows if row.strip()]


def test_every_upstream_file_is_unmodified():
    manifest = _manifest()
    assert len(manifest) > 600
    for digest, name in manifest:
        assert hashlib.sha256((VENDOR / name).read_bytes()).hexdigest() == digest, name


def test_line_endings_are_protected_from_git():
    assert "* -text" in (VENDOR / ".gitattributes").read_text()


def test_the_idate_patch_targets_exactly_one_line():
    lines = (VENDOR / "Source" / "runatom.f").read_text(encoding="latin-1").splitlines()
    assert lines.count(IDATE_LINE) == 1
    patch = (VENDOR / "patches" / "0001-gfortran-intrinsic-idate.patch").read_text()
    assert f"-{IDATE_LINE}" in patch and "+c     external idate" in patch


@pytest.mark.parametrize("name", ["runatom.exe", "reform.exe"])
def test_built_program_matches_its_recorded_checksum(name):
    exe = VENDOR / "bin" / name
    if not exe.exists():
        pytest.skip(f"{name} has not been built yet (Task 2)")
    recorded = (VENDOR / "bin" / f"{name}.sha256").read_text().split()[0]
    assert hashlib.sha256(exe.read_bytes()).hexdigest() == recorded
```

- [ ] **Step 6: Run it**

Run: `python -m pytest tests/test_runatomVendor.py -v`
Expected: 3 passed, 2 skipped. (Before Step 1 it fails with `FileNotFoundError`.)

- [ ] **Step 7: Commit (if authorised)**

```bash
git add third_party/uspp-736 tests/test_runatomVendor.py
git commit -m "build(uspp): vendor unmodified USPP 7.3.6 with provenance and checksum manifest"
```

---

### Task 2: Build the static Windows `runatom.exe` and `reform.exe`

**Prerequisite (done 2026-10-28, approved by the user):** MSYS2 was installed in `C:\msys64` with `winget install --id MSYS2.MSYS2 -e`. Then, from PowerShell with `MSYSTEM=UCRT64`:
- `C:\msys64\usr\bin\bash.exe -lc "pacman -Syuu --noconfirm"`, run twice;
- `pacman -S --needed --noconfirm mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-binutils`.

Installed and verified: GNU Fortran 16.2.0 (Rev6, MSYS2) and GNU objdump 2.47. This is newer than the WSL reference's gfortran 11.5, which is why Task 10 compares with tolerances. Once the two programs are committed, nothing else in the plan needs MSYS2.

**Files:**
- Create: `third_party/uspp-736/build_runatom.ps1`
- Create (by running it): `third_party/uspp-736/bin/runatom.exe`, `bin/runatom.exe.sha256`, `bin/BUILDINFO.txt`
- Test: `tests/test_runatomVendor.py` (already has the checksum test)

**Interfaces:**
- Produces: `third_party/uspp-736/bin/runatom.exe`, run as `runatom.exe <inp> <out> <ae> <atwf> <logd> <uspp|dummy>` from the job folder.

- [ ] **Step 1: Write the build script**

Create `third_party/uspp-736/build_runatom.ps1`:

```powershell
<#
Builds bin\runatom.exe from the unmodified USPP 7.3.6 source in Source\.

Needs MSYS2 UCRT64 with gcc-fortran and binutils (see PROVENANCE.md).
The idate patch is applied to a temporary copy; Source\ is never changed.
The result is statically linked: it may import only KERNEL32.dll and the
Windows 10+ UCRT (api-ms-win-crt-*.dll), or the build fails.

Usage (from venv\src):
    powershell -ExecutionPolicy Bypass -File third_party\uspp-736\build_runatom.ps1 [-Msys C:\msys64]
#>
param([string]$Msys = "C:\msys64")
$ErrorActionPreference = "Stop"

$root = $PSScriptRoot
$bin = Join-Path $Msys "ucrt64\bin"
$gfortran = Join-Path $bin "gfortran.exe"
$objdump = Join-Path $bin "objdump.exe"
if (-not (Test-Path $gfortran)) { throw "gfortran not found at $gfortran: install mingw-w64-ucrt-x86_64-gcc-fortran" }
if (-not (Test-Path $objdump)) { throw "objdump not found at $objdump: install mingw-w64-ucrt-x86_64-binutils" }
$env:PATH = "$bin;$env:PATH"
$flags = @("-O1", "-std=legacy", "-fallow-argument-mismatch", "-static")

$work = Join-Path ([IO.Path]::GetTempPath()) ("uspp-build-" + [guid]::NewGuid())
New-Item -ItemType Directory $work | Out-Null
try {
    Copy-Item (Join-Path $root "Source\*.f") $work
    $main = Join-Path $work "runatom.f"
    $text = [IO.File]::ReadAllText($main)
    $pattern = "(?m)^      external idate(\r?)$"
    if ([regex]::Matches($text, $pattern).Count -ne 1) { throw "the idate patch no longer matches runatom.f" }
    [IO.File]::WriteAllText($main, [regex]::Replace($text, $pattern, 'c     external idate$1'))

    Copy-Item (Join-Path $root "Utils\reform.f") $work      # unmodified: it needs no patch

    function Assert-SelfContained([string]$exe) {
        $dlls = & $objdump -p $exe | Select-String "DLL Name:" | ForEach-Object { ($_ -split ":", 2)[1].Trim() }
        $foreign = @($dlls | Where-Object { $_ -ne "KERNEL32.dll" -and $_ -notlike "api-ms-win-crt-*" })
        if ($foreign.Count -gt 0) { throw "$exe is not self-contained; it imports: $($foreign -join ', ')" }
        return ($dlls -join ', ')
    }

    Push-Location $work
    try {
        $sources = Get-ChildItem -Filter *.f | Where-Object Name -ne "reform.f" | ForEach-Object Name
        & $gfortran @flags -o runatom.exe @sources
        if ($LASTEXITCODE -ne 0) { throw "gfortran failed on runatom with exit code $LASTEXITCODE" }
        & $gfortran @flags -o reform.exe reform.f
        if ($LASTEXITCODE -ne 0) { throw "gfortran failed on reform with exit code $LASTEXITCODE" }
        $imports = Assert-SelfContained "runatom.exe"
        $reformImports = Assert-SelfContained "reform.exe"
    } finally { Pop-Location }

    # Smoke run: the package's own example must give a pseudopotential and its CPMD formatted copy.
    $smoke = Join-Path $work "smoke"
    Copy-Item (Join-Path $root "Work\000-Ex\000-Ex-ca-sp-test") $smoke -Recurse
    Push-Location $smoke
    try {
        & (Join-Path $work "runatom.exe") ti_ae_s2d1.adat ti_ae_s2d1.out ti_ae_s2d1.ae ti_ae_s2d1.atwf ti_ae_s2d1.logd dummy 2>$null
        & (Join-Path $work "runatom.exe") ti_ps.adat ti_ps.out ti_ae_s2d1.ae ti_ps.atwf ti_ps.logd ti_ps.uspp 2>$null
        if (-not (Test-Path ti_ps.uspp)) { throw "smoke run produced no ti_ps.uspp; see $smoke" }
        "1`nti_ps.uspp`nti_ps_formatted.uspp`n" | & (Join-Path $work "reform.exe") | Out-Null
        if (-not (Test-Path ti_ps_formatted.uspp)) { throw "reform.exe produced no ti_ps_formatted.uspp; see $smoke" }
    } finally { Pop-Location }

    $out = Join-Path $root "bin"
    New-Item -ItemType Directory $out -Force | Out-Null
    $lines = @()
    foreach ($name in "runatom.exe", "reform.exe") {
        Copy-Item (Join-Path $work $name) $out -Force
        $hash = (Get-FileHash (Join-Path $out $name) -Algorithm SHA256).Hash.ToLower()
        [IO.File]::WriteAllText((Join-Path $out "$name.sha256"), "$hash  $name`n")
        $lines += "$name sha256: $hash"
    }
    $version = (& $gfortran --version | Select-Object -First 1)
    [IO.File]::WriteAllText((Join-Path $out "BUILDINFO.txt"),
        "compiler: $version`nflags: $($flags -join ' ')`n" +
        "runatom.exe: Source/*.f + patches/0001-gfortran-intrinsic-idate.patch; imports: $imports`n" +
        "reform.exe: Utils/reform.f (unmodified); imports: $reformImports`n" + ($lines -join "`n") + "`n")
    Write-Host ($lines -join "`n")
} finally {
    Remove-Item $work -Recurse -Force -ErrorAction SilentlyContinue
}
```

- [ ] **Step 2: Run the build**

Run: `powershell -ExecutionPolicy Bypass -File third_party\uspp-736\build_runatom.ps1`
Expected: `runatom.exe sha256: <hash>` and `reform.exe sha256: <hash>`.

If gfortran rejects the code, record the exact message in `BUILDINFO.txt`'s commit message. Fix it only by adding a compiler flag to `$flags`, never by editing the source. Then add the same flag to `make_reference.sh` in Task 3.

- [ ] **Step 3: Run the checksum test**

Run: `python -m pytest tests/test_runatomVendor.py -v`
Expected: 5 passed.

- [ ] **Step 4: Commit (if authorised)**

```bash
git add third_party/uspp-736/build_runatom.ps1 third_party/uspp-736/bin
git commit -m "build(uspp): static MSYS2 UCRT64 build of runatom.exe with import check and smoke run"
```

---

### Task 3: Linux reference outputs and test fixtures

The reference is the original program built with Linux gfortran plus the same one-line patch. This build was verified against the upstream `Sample_output`.

**Files:**
- Create: `tests/data/runatom/make_reference.sh`, `tests/data/runatom/README.md`
- Create (by running it): `tests/data/runatom/MANIFEST.txt`, `tests/data/runatom/<case>/{CASE.txt,*.out.gz,*_ps.uspp.gz}`
- Create: `tests/runatom_fixtures.py`
- Test: `tests/test_runatomReferenceData.py`

**Interfaces:**
- Produces (in `tests/runatom_fixtures.py`):
  - `LIBRARY_DIR: Path`, `REFERENCE_DIR: Path`, `REFERENCE_CASES: Dict[str, str]` (case name → `Work/` relative path);
  - `library_case_dir(name) -> Path`, `case_info(name) -> Dict[str, str]` (`atom`, `cfg`, `tests`);
  - `reference_names(name) -> List[str]` (names of stored files without `.gz`);
  - `reference_text(name, filename) -> str`, `reference_bytes(name, filename) -> bytes`.

- [ ] **Step 1: Write the reference script**

Create `tests/data/runatom/make_reference.sh`:

```bash
#!/bin/bash
# Dev-time only, never run by the test suite. Builds the reference runatom.x
# with Linux gfortran from third_party/uspp-736/Source plus the idate patch,
# runs 8 library cases and stores their .out and .uspp files (gzipped) here.
# Run from venv/src inside WSL:  bash tests/data/runatom/make_reference.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$(cd "$HERE/../../../third_party/uspp-736" && pwd)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
FLAGS="-O1 -std=legacy -fallow-argument-mismatch"

mkdir "$TMP/build"; cp "$SRC"/Source/*.f "$TMP/build/"
sed -i 's/^      external idate\(\r\?\)$/c     external idate\1/' "$TMP/build/runatom.f"
[ "$(grep -c '^c     external idate' "$TMP/build/runatom.f")" = 1 ]
(cd "$TMP/build" && gfortran $FLAGS -o runatom.x *.f)
X="$TMP/build/runatom.x"
mkdir "$TMP/reform"; cp "$SRC/Utils/reform.f" "$TMP/reform/"          # unmodified
(cd "$TMP/reform" && gfortran $FLAGS -o reform.x reform.f)
R="$TMP/reform/reform.x"

CASES="000-Ex/000-Ex-ca-sp-test 001-H/001-H-gpbe--bm 008-O/008-O-ca--vgrp
       026-Fe/026-Fe-gpw-n-campos 029-Cu/029-Cu-gpbe-n-bm 079-Au/079-Au-ca--ts
       082-Pb/082-Pb-ca-d-vgrp 056-Ba/056-Ba-gpw-nsp-campos"
{ echo "compiler: $(gfortran --version | head -1)"; echo "flags: $FLAGS";
  echo "source: third_party/uspp-736/Source + patches/0001-gfortran-intrinsic-idate.patch"; } > "$HERE/MANIFEST.txt"

for c in $CASES; do
  name=$(basename "$c"); dst="$HERE/$name"; rm -rf "$dst"; mkdir -p "$dst"
  run="$TMP/$name"; cp -r "$SRC/Work/$c" "$run"; cd "$run"
  atom=$(sed -n 's/^ATOM= *//p' Makefile | tr -d '\r')
  cfg=$(sed -n 's/^CFG= *//p' Makefile | tr -d '\r')
  tests=$(sed -n 's/^CFT[0-9]*= *//p' Makefile | tr -d '\r')
  "$X" "${atom}_ae_${cfg}.adat" "${atom}_ae_${cfg}.out" "${atom}_ae_${cfg}.ae" "${atom}_ae_${cfg}.atwf" "${atom}_ae_${cfg}.logd" dummy 2>/dev/null
  "$X" "${atom}_ps.adat" "${atom}_ps.out" "${atom}_ae_${cfg}.ae" "${atom}_ps.atwf" "${atom}_ps.logd" "${atom}_ps.uspp" 2>/dev/null
  done_tests=""
  for t in $tests; do
    [ -f "${atom}_test_${t}.adat" ] || continue
    "$X" "${atom}_ae_${t}.adat" "${atom}_ae_${t}.out" "${atom}_ae_${t}.ae" "${atom}_ae_${t}.atwf" "${atom}_ae_${t}.logd" dummy 2>/dev/null
    "$X" "${atom}_test_${t}.adat" "${atom}_test_${t}.out" "${atom}_ae_${t}.ae" "${atom}_test_${t}.atwf" "${atom}_test_${t}.logd" "${atom}_ps.uspp" 2>/dev/null
    done_tests="$done_tests $t"
  done
  printf '1\n%s\n%s\n' "${atom}_ps.uspp" "${atom}_ps_formatted.uspp" | "$R" > /dev/null
  for f in *.out "${atom}_ps.uspp" "${atom}_ps_formatted.uspp"; do gzip -9 -n -c "$f" > "$dst/$f.gz"; done
  { echo "case: $c"; echo "atom: $atom"; echo "cfg: $cfg"; echo "tests:$done_tests"; } > "$dst/CASE.txt"
  cd "$HERE"
done
echo "reference written to $HERE"
```

Create `tests/data/runatom/README.md`:

```markdown
Reference outputs of the original USPP 7.3.6 runatom, and of reform (itype 1,
the formatted file CPMD reads), from a Linux gfortran build (see MANIFEST.txt),
for 8 library cases chosen to exercise every input branch:
LDA/PBE/PW91, semicore, nonlinear core correction, per-L rinner, iptype 2/3/4,
ifqopt 0/2/3, nfix constraints and transferability tests.
Regenerate only with make_reference.sh (WSL with gfortran); never by hand.
The .uspp files carry the run date in record 1; comparisons mask it.
```

- [ ] **Step 2: Run it once in WSL**

Run (PowerShell): `wsl.exe -d OracleLinux_9_5 -- bash tests/data/runatom/make_reference.sh`
Expected: `reference written to ...`, plus 8 folders, each with `CASE.txt`, `*_ps.uspp.gz`, `*_ps_formatted.uspp.gz` and at least two `*.out.gz`.

- [ ] **Step 3: Write the fixtures module**

Create `tests/runatom_fixtures.py`:

```python
"""Paths and readers for the USPP library and the runatom reference outputs."""

import gzip
from pathlib import Path
from typing import Dict, List

SRC_DIR = Path(__file__).resolve().parents[1]
LIBRARY_DIR = SRC_DIR / "third_party" / "uspp-736" / "Work"
REFERENCE_DIR = Path(__file__).resolve().parent / "data" / "runatom"
REFERENCE_CASES = {
    "000-Ex-ca-sp-test": "000-Ex/000-Ex-ca-sp-test",
    "001-H-gpbe--bm": "001-H/001-H-gpbe--bm",
    "008-O-ca--vgrp": "008-O/008-O-ca--vgrp",
    "026-Fe-gpw-n-campos": "026-Fe/026-Fe-gpw-n-campos",
    "029-Cu-gpbe-n-bm": "029-Cu/029-Cu-gpbe-n-bm",
    "079-Au-ca--ts": "079-Au/079-Au-ca--ts",
    "082-Pb-ca-d-vgrp": "082-Pb/082-Pb-ca-d-vgrp",
    "056-Ba-gpw-nsp-campos": "056-Ba/056-Ba-gpw-nsp-campos",
}


def library_case_dir(name: str) -> Path:
    return LIBRARY_DIR / REFERENCE_CASES[name]


def case_info(name: str) -> Dict[str, str]:
    rows = (REFERENCE_DIR / name / "CASE.txt").read_text().splitlines()
    return {k.strip(): v.strip() for k, v in (row.split(":", 1) for row in rows if ":" in row)}


def reference_names(name: str) -> List[str]:
    return sorted(p.name[:-3] for p in (REFERENCE_DIR / name).glob("*.gz"))


def reference_bytes(name: str, filename: str) -> bytes:
    return gzip.decompress((REFERENCE_DIR / name / f"{filename}.gz").read_bytes())


def reference_text(name: str, filename: str) -> str:
    return reference_bytes(name, filename).decode("latin-1")
```

- [ ] **Step 4: Write the test**

Create `tests/test_runatomReferenceData.py`:

```python
"""The stored reference outputs are complete and came from the documented build."""

import pytest

from runatom_fixtures import REFERENCE_CASES, REFERENCE_DIR, case_info, reference_names, reference_text


def test_manifest_records_compiler_flags_and_patch():
    manifest = (REFERENCE_DIR / "MANIFEST.txt").read_text()
    assert "gfortran" in manifest and "-O1 -std=legacy" in manifest and "idate" in manifest


@pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
def test_every_case_has_ae_pseudo_and_test_outputs(name):
    info = case_info(name)
    names = reference_names(name)
    atom, cfg = info["atom"], info["cfg"]
    assert f"{atom}_ae_{cfg}.out" in names and f"{atom}_ps.out" in names and f"{atom}_ps.uspp" in names
    assert f"{atom}_ps_formatted.uspp" in names       # reform.x itype=1: the CPMD file
    for test in info["tests"].split():
        assert f"{atom}_test_{test}.out" in names
    assert "writing the pseudopotential" not in reference_text(name, f"{atom}_ps.out")  # stderr text, not .out
    assert "comparison of all-electron and pseudo eigenvalues" in reference_text(name, f"{atom}_ps.out")
```

- [ ] **Step 5: Run it**

Run: `python -m pytest tests/test_runatomReferenceData.py -v`
Expected: 9 passed.

- [ ] **Step 6: Commit (if authorised)**

```bash
git add tests/data/runatom tests/runatom_fixtures.py tests/test_runatomReferenceData.py
git commit -m "test(uspp): Linux reference outputs of 8 library cases and fixture helpers"
```

---

### Task 4: Fixed-column Fortran fields (`runatomInput.py`, part 1)

**Files:**
- Create: `runatomInput.py`
- Test: `tests/test_runatomInput.py`

**Interfaces:**
- Produces:
  - `parse_format(fmt) -> List[Tuple[str, int, int]]`
  - `read_field(kind, width, decimals, raw)`
  - `read_record(line, fmt, count=None) -> list`
  - `write_field(kind, width, decimals, value) -> str`
  - `format_line(fmt, values, label) -> List[str]` (more than one line when the values exceed the format: Fortran format reversion)
  - constants `LABEL_COLUMN`, `MAX_FORTRAN_NAME`, `SAFE_NAME`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_runatomInput.py`:

```python
"""runatomInput: fixed-column .adat files exactly as runatom (rwsubs.f) reads them."""

import random
import subprocess
import sys
from pathlib import Path

import pytest

import runatomInput as ri

SRC = Path(__file__).resolve().parents[1]


def test_module_does_not_import_toga():
    code = "import sys, runatomInput; assert 'toga' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], cwd=SRC, check=True)


class TestReadField:
    def test_blanks_inside_an_integer_field_are_ignored(self):
        # gfortran BLANK='NULL': "1 0 0" in an i5 field is 100 (verified with gfortran 11.5)
        assert ri.read_record("1 0 0 0 3", "6i5", 5) == [100, 3, 0, 0, 0]

    def test_f_field_without_decimal_point_is_scaled(self):
        assert ri.read_record("   2.5        2         1.8", "3f10.5") == [2.5, 2e-5, 1.8]

    def test_explicit_decimal_point_wins(self):
        assert ri.read_field("f", 10, 5, "      2.  ") == 2.0
        assert ri.read_field("f", 5, 0, "  22.") == 22.0
        assert ri.read_field("f", 14, 6, "      -.55710 ") == -0.5571

    def test_fortran_exponents(self):
        assert ri.read_field("e", 10, 1, "   1.0d-10") == 1e-10
        assert ri.read_field("f", 10, 5, "1.23-100  ") == 1.23e-100

    def test_empty_and_short_records_read_as_zero(self):
        assert ri.read_record("    7", "2i5,f10.5") == [7, 0, 0.0]

    def test_sign_without_exponent_is_a_read_error(self):
        # Work/023-V-gpw-nsp-campos/V_ps.adat line 2: runatom stops here
        with pytest.raises(ValueError, match="floating point"):
            ri.read_record("  2.21   -6.0       2.0       80", "3f10.5,i5")

    def test_text_in_an_integer_field_is_a_read_error(self):
        with pytest.raises(ValueError):
            ri.read_field("i", 5, 0, "  1.0")


class TestWriteField:
    def test_every_float_is_written_with_a_decimal_point(self):
        for value in (2.0, 0.0, -6.0, 82.0, 1e-10, 150.0, -9000.0):
            assert "." in ri.write_field("f", 10, 5, value)

    def test_written_floats_read_back(self):
        rng = random.Random(7)
        for _ in range(2000):
            typed = round(rng.uniform(-1000, 1000), 4)                 # what people type: exact
            assert ri.read_field("f", 14, 6, ri.write_field("f", 14, 6, typed)) == typed
            tiny = rng.choice([1, -1]) * 10 ** rng.uniform(-9, -4)    # needs an exponent: 7 digits kept
            assert ri.read_field("f", 14, 6, ri.write_field("f", 14, 6, tiny)) == pytest.approx(tiny, rel=1e-6)

    def test_exact_value_may_fill_the_field(self):
        assert ri.write_field("f", 5, 0, 100.0) == "100.0"

    def test_integers_are_right_justified_and_checked(self):
        assert ri.write_field("i", 5, 0, 3) == "    3"
        with pytest.raises(ValueError):
            ri.write_field("i", 4, 0, 12345)

    def test_text_longer_than_the_field_is_refused(self):
        with pytest.raises(ValueError):
            ri.write_field("a", 20, 0, "a" * 21)


class TestFormatLine:
    def test_label_never_inside_the_columns_read(self):
        (line,) = ri.format_line("6i5", [1, 0, 0, 0, 3], "ifae,ifpsp,ifprt,ifplw,ilogd")
        assert line.index("ifae") >= ri.LABEL_COLUMN
        assert ri.read_record(line, "6i5", 5) == [1, 0, 0, 0, 3]

    def test_format_reversion_continues_on_a_new_line(self):
        values = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]  # nang = 4: 7 rinner values, format 5f10.5
        lines = ri.format_line("5f10.5", values, "rinner")
        assert len(lines) == 2
        assert ri.read_record(lines[0], "5f10.5") + ri.read_record(lines[1], "5f10.5", 2) == values
```

- [ ] **Step 2: Run them to see them fail**

Run: `python -m pytest tests/test_runatomInput.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'runatomInput'`.

- [ ] **Step 3: Implement**

Create `runatomInput.py`:

```python
"""Input files (.adat) of the Vanderbilt USPP program runatom (USPP 7.3.6).

No Toga here (a subprocess test pins it). runatom reads every input line with a
fixed-column Fortran format (rwsubs.f, subroutine readin). gfortran ignores
blanks inside a numeric field and scales an F/E field that has no decimal point
by 10**-d: "1 0 0" in an i5 field is 100 and "2" in an f10.5 field is 2e-5.
Text typed by a user is therefore never copied into a file: values are parsed,
checked and written column-exactly here.
"""

import math
import re
from typing import List, Optional, Sequence, Tuple

LABEL_COLUMN = 38           # labels start here, or 2 columns after the data
MAX_FORTRAN_NAME = 40       # runatom's file-name arguments are character*40
SAFE_NAME = re.compile(r"^[A-Za-z0-9._-]+$")

_FORMAT_ITEM = re.compile(r"^(\d*)([iafe])(\d+)(?:\.(\d+))?$")
_FORTRAN_REAL = re.compile(r"^([+-]?)(\d*)(\.?)(\d*)(?:[ED]([+-]?\d+)|([+-]\d+))?$")


def parse_format(fmt: str) -> List[Tuple[str, int, int]]:
    """'2i5,f10.5' -> [('i', 5, 0), ('i', 5, 0), ('f', 10, 5)]."""
    fields = []
    for item in fmt.split(","):
        match = _FORMAT_ITEM.match(item.strip())
        if match is None:
            raise ValueError(f"unsupported Fortran format item {item!r}")
        repeat = int(match.group(1) or 1)
        fields.extend([(match.group(2), int(match.group(3)), int(match.group(4) or 0))] * repeat)
    return fields


def read_field(kind: str, width: int, decimals: int, raw: str):
    """One input field as gfortran reads it (BLANK='NULL': blanks ignored)."""
    if kind == "a":
        return raw
    text = raw.replace(" ", "")
    if kind == "i":
        if not text:
            return 0
        if not re.fullmatch(r"[+-]?\d+", text):
            raise ValueError(f"Bad integer {raw.strip()!r}")
        return int(text)
    if not text:
        return 0.0
    match = _FORTRAN_REAL.match(text.upper())
    if match is None or not (match.group(2) or match.group(4)):
        raise ValueError(f"Bad value during floating point read: {raw.strip()!r}")
    sign, whole, point, frac, exp_e, exp_sign = match.groups()
    exponent = int(exp_e or exp_sign or 0)
    if point:
        return float(f"{sign}{whole or '0'}.{frac or '0'}e{exponent}")
    return float(f"{sign}{whole}e{exponent - decimals}")   # implied decimal point


def read_record(line: str, fmt: str, count: Optional[int] = None) -> list:
    """The first `count` fields of one record (a short record is padded with blanks)."""
    fields = parse_format(fmt)
    if count is not None:
        fields = fields[:count]     # a READ stops when its list is exhausted
    line = line.rstrip("\r\n")
    values, position = [], 0
    for kind, width, decimals in fields:
        values.append(read_field(kind, width, decimals, line[position:position + width].ljust(width)))
        position += width
    return values


def _float_texts(value: float):
    """Spellings of a float, most precise first; every one has a decimal point."""
    text = repr(value)
    if "e" in text:
        mantissa, exponent = text.split("e")
        yield f"{mantissa if '.' in mantissa else mantissa + '.0'}E{int(exponent)}"
    else:
        yield text
    if value == 0.0 or 1e-3 <= abs(value) < 1e15:
        for decimals in range(12, 0, -1):
            yield f"{value:.{decimals}f}"
    for digits in range(12, 0, -1):
        yield f"{value:.{digits}E}"


def write_field(kind: str, width: int, decimals: int, value) -> str:
    """One output field; floats keep one leading blank unless only the exact value fits."""
    if kind == "a":
        text = str(value)
        if len(text) > width:
            raise ValueError(f"{text!r} is longer than {width} characters")
        return text.ljust(width)
    if kind == "i":
        if isinstance(value, float) and not value.is_integer():
            raise ValueError(f"{value} is not an integer")
        text = str(int(value))
        if len(text) > width:
            raise ValueError(f"{value} does not fit in {width} columns")
        return text.rjust(width)
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{value!r} is not a finite number")
    for text in _float_texts(number):
        if len(text) < width or (len(text) == width and float(text) == number):
            return text.rjust(width)
    raise ValueError(f"{value} does not fit in {width} columns")


def format_line(fmt: str, values: Sequence, label: str) -> List[str]:
    """Data written column-exactly, then '<label> (<fmt>)' after the columns read.

    More values than the format has fields continue on a new line, as Fortran
    format reversion does when it reads them.
    """
    fields = parse_format(fmt)
    lines = []
    for start in range(0, max(len(values), 1), len(fields)):
        chunk = values[start:start + len(fields)]
        lines.append("".join(write_field(k, w, d, v) for (k, w, d), v in zip(fields, chunk)))
    lines[0] = lines[0].ljust(max(LABEL_COLUMN, len(lines[0]) + 2)) + f"{label} ({fmt})"
    return lines
```

- [ ] **Step 4: Run them to see them pass**

Run: `python -m pytest tests/test_runatomInput.py -v`
Expected: all passed.

- [ ] **Step 5: Commit (if authorised)**

```bash
git add runatomInput.py tests/test_runatomInput.py
git commit -m "feat(uspp): fixed-column Fortran field reader/writer for runatom inputs"
```

---

### Task 5: Input data model, `.adat` writer/reader and library case folders (`runatomInput.py`, part 2)

**Files:**
- Modify: `runatomInput.py` (append)
- Test: `tests/test_runatomInput.py` (append)

**Interfaces:**
- Consumes: Task 4 functions; `tests/runatom_fixtures.LIBRARY_DIR`.
- Produces:
  - dataclasses `Controls`, `Shell`, `AEInput`, `Beta`, `QFix`, `PSInput`, `TestInput`, `TestConfig`, `CaseInputs`;
  - `LINE_FORMATS`; `AdatError(ValueError)` with `.line_number`;
  - `format_adat(obj) -> str`, `parse_adat(text) -> AEInput | PSInput | TestInput`;
  - `file_stems(case) -> {"ae": str, "ps": str, "tests": List[Tuple[str, str]]}`;
  - `test_ae_input(case, test) -> AEInput`, `test_input(case, test) -> TestInput`;
  - `case_files(case) -> Dict[str, str]`, `write_case(case, folder) -> List[Path]`;
  - `read_case_folder(folder) -> CaseInputs`, `makefile_text(case) -> str`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_runatomInput.py`:

```python
from runatom_fixtures import LIBRARY_DIR

LIBRARY_INPUTS = sorted(p for p in LIBRARY_DIR.rglob("*.adat") if "History" not in p.parts)
UNREADABLE = {"V_ps.adat": 2}           # upstream defect, see PROVENANCE.md
CASE_DIRS = sorted(p.parent for p in LIBRARY_DIR.rglob("Makefile") if p.parent.parent != LIBRARY_DIR)


def _ti_case():
    return ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-test")


class TestAdatFiles:
    def test_library_is_large_enough_to_mean_something(self):
        assert len(LIBRARY_INPUTS) > 400

    @pytest.mark.parametrize("path", LIBRARY_INPUTS, ids=lambda p: f"{p.parent.name}/{p.name}")
    def test_every_library_input_round_trips(self, path):
        text = path.read_text(encoding="latin-1")
        if path.name in UNREADABLE:
            with pytest.raises(ri.AdatError) as info:
                ri.parse_adat(text)
            assert info.value.line_number == UNREADABLE[path.name]
            return
        parsed = ri.parse_adat(text)
        assert ri.parse_adat(ri.format_adat(parsed)) == parsed

    def test_written_file_uses_lf_and_ends_with_newline(self):
        text = ri.format_adat(_ti_case().ps)
        assert "\r" not in text and text.endswith("\n")

    def test_npf_line_only_when_a_projector_is_smooth(self):
        ps = _ti_case().ps
        assert "npf,ptryc" in ri.format_adat(ps)
        plain = ri.replace(ps, betas=[ri.replace(b, iptype=1) for b in ps.betas])
        text = ri.format_adat(plain)
        assert "npf,ptryc" not in text
        assert ri.parse_adat(text).lloc == ps.lloc   # the local line was not shifted

    def test_per_l_rinner_line_when_rinner_is_zero(self):
        case = ri.read_case_folder(LIBRARY_DIR / "026-Fe" / "026-Fe-gpw-n-campos")
        assert case.ps.rinner1 == 0.0 and len(case.ps.rinner) == 2 * case.ps.nang - 1


class TestCaseFolders:
    @pytest.mark.parametrize("folder", CASE_DIRS, ids=lambda p: p.name)
    def test_every_library_case_reads(self, folder):
        if (folder / "V_ps.adat").exists():
            with pytest.raises(ri.AdatError):
                ri.read_case_folder(folder)
            return
        case = ri.read_case_folder(folder)
        assert case.ae.controls.ifae == 1 and case.ps.controls.ifpsp == 2

    def test_tests_named_cf1_cf2_are_found(self):
        case = ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-tests")
        assert [t.name for t in case.tests] == ["s2d2", "s1d2"]

    def test_tests_become_occupation_changes(self):
        case = _ti_case()
        assert case.symbol == "ti" and case.generation_name == "s2d1"
        (test,) = case.tests
        assert test.name == "s2d2" and test.xion == 0.0 and test.occupations == {320: 2.0}
        derived = ri.test_ae_input(case, test)
        assert {s.nnlz: s.wwnl for s in derived.shells}[320] == 2.0

    def test_written_case_reads_back(self, tmp_path):
        case = _ti_case()
        written = ri.write_case(case, tmp_path)
        assert sorted(p.name for p in written) == ["ti_ae_s2d1.adat", "ti_ae_s2d2.adat",
                                                   "ti_ps.adat", "ti_test_s2d2.adat"]
        (tmp_path / "Makefile").write_text(ri.makefile_text(case), newline="\n")
        again = ri.read_case_folder(tmp_path)
        assert again.ae == case.ae and again.ps == case.ps
        assert [(t.name, t.xion, t.occupations) for t in again.tests] == [("s2d2", 0.0, {320: 2.0})]

    def test_makefile_lists_every_stage(self):
        text = ri.makefile_text(_ti_case())
        assert "ATOM= ti" in text and "CFG= s2d1" in text and "CFT1= s2d2" in text
        assert "${PROG} ${TST1}.adat" in text
```

- [ ] **Step 2: Run them to see them fail**

Run: `python -m pytest tests/test_runatomInput.py -v -k "Adat or CaseFolders"`
Expected: FAIL with `AttributeError: module 'runatomInput' has no attribute 'read_case_folder'`.

- [ ] **Step 3: Implement**

Append to `runatomInput.py`. Move the `from dataclasses …`, `from pathlib …` and `from typing …` imports to the top of the file:

```python
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

LINE_FORMATS = {
    "controls": ("6i5", "ifae,ifpsp,ifprt,ifplw,ilogd"),
    "logd": ("3f10.5,i5", "rlogd,emin,emax,nnt"),
    "convergence": ("2e10.1,f10.5,i5", "thresh,tol,damp,maxit"),
    "title": ("a20", "title"),
    "atom": ("f5.0,2f10.5", "z,xion,exfact"),
    "mesh": ("3f10.5", "rmax,aasf,bbsf"),
    "nstates": ("2i5", "ncspvs,irel"),
    "shell": ("i4,f7.3,f14.6", "nnlz,wwnl,ee"),
    "partition": ("3i5", "ncores,nvales,nang"),
    "bessel": ("4f10.5", "besrmax,besemin,besemax,besde"),
    "keyps": ("2i5,2f10.5", "keyps,ifpcor,rinner,rpcor"),
    "rinner": ("5f10.5", "rinner for L=1..2*nang-1"),
    "nbeta": ("i5,f10.5", "nbeta,rcloc"),
    "rc": ("4f10.5", "rc for l=0..nang-1"),
    "beta": ("2i5,f10.5,i5", "lll,keyee,eeread,iptype"),
    "npf": ("i5,f10.5", "npf,ptryc"),
    "local": ("2i5,f10.5,i5", "lloc,keyee,eloc,iploctype"),
    "qfunc": ("2i5,f10.5,i5", "ifqopt,nqf,qtryc,nfix"),
    "qfix": ("i5,2f10.5", "ibfix,rfix,qfix"),
}


@dataclass
class Controls:
    """The four lines every .adat file starts with."""
    ifae: int
    ifpsp: int
    ifprt: int
    ifplw: int
    ilogd: int
    rlogd: float
    emin: float
    emax: float
    nnt: int
    thresh: float
    tol: float
    damp: float
    maxit: int
    title: str


@dataclass
class Shell:
    nnlz: int       # n, l, 0 as three digits (320 = 3d)
    wwnl: float     # electrons (negative: orbital skipped)
    ee: float       # starting guess for the eigenvalue (Ry)


@dataclass
class AEInput:
    controls: Controls
    z: float
    xion: float
    exfact: float
    rmax: float
    aasf: float
    bbsf: float
    irel: int
    shells: List[Shell]


@dataclass
class Beta:
    l: int
    keyee: int
    eeread: float
    iptype: int


@dataclass
class QFix:
    ibeta: int
    rfix: float
    qfix: float


@dataclass
class PSInput:
    controls: Controls
    ncores: int
    nvales: int
    nang: int
    besrmax: float
    besemin: float
    besemax: float
    besde: float
    keyps: int
    ifpcor: int
    rinner1: float
    rpcor: float
    rinner: List[float]      # 2*nang-1 values, written only when rinner1 == 0
    rcloc: float
    rc: List[float]          # nang values
    betas: List[Beta]
    npf: int                 # npf, ptryc written only when some beta has iptype 2
    ptryc: float
    lloc: int
    keyee_loc: int
    eloc: float
    iploctype: int
    ifqopt: int
    nqf: int
    qtryc: float
    fixes: List[QFix]


@dataclass
class TestInput:
    controls: Controls
    ncores: int
    nvales: int
    nang: int
    besrmax: float
    besemin: float
    besemax: float
    besde: float


@dataclass
class TestConfig:
    """A transferability test: the generation configuration with other occupations."""
    name: str
    xion: float
    occupations: Dict[int, float]                # nnlz -> electrons; others unchanged
    test_input: Optional[TestInput] = None       # None: derived from the generation input


@dataclass
class CaseInputs:
    symbol: str
    generation_name: str
    ae: AEInput
    ps: PSInput
    tests: List[TestConfig] = field(default_factory=list)


class AdatError(ValueError):
    def __init__(self, line_number: int, message: str):
        super().__init__(f"line {line_number}: {message}")
        self.line_number = line_number


def _line(key: str, values: Sequence) -> List[str]:
    fmt, label = LINE_FORMATS[key]
    try:
        return format_line(fmt, list(values), label)
    except ValueError as exc:
        raise ValueError(f"{label}: {exc}") from None


def _control_lines(c: Controls) -> List[str]:
    return (_line("controls", [c.ifae, c.ifpsp, c.ifprt, c.ifplw, c.ilogd])
            + _line("logd", [c.rlogd, c.emin, c.emax, c.nnt])
            + _line("convergence", [c.thresh, c.tol, c.damp, c.maxit])
            + _line("title", [c.title]))


def format_adat(data: Union[AEInput, PSInput, TestInput]) -> str:
    lines = _control_lines(data.controls)
    if isinstance(data, AEInput):
        lines += _line("atom", [data.z, data.xion, data.exfact])
        lines += _line("mesh", [data.rmax, data.aasf, data.bbsf])
        lines += _line("nstates", [len(data.shells), data.irel])
        for s in data.shells:
            lines += _line("shell", [s.nnlz, s.wwnl, s.ee])
    else:
        lines += _line("partition", [data.ncores, data.nvales, data.nang])
        lines += _line("bessel", [data.besrmax, data.besemin, data.besemax, data.besde])
    if isinstance(data, PSInput):
        lines += _line("keyps", [data.keyps, data.ifpcor, data.rinner1, data.rpcor])
        if data.rinner1 == 0.0:
            lines += _line("rinner", data.rinner)
        lines += _line("nbeta", [len(data.betas), data.rcloc])
        lines += _line("rc", data.rc)
        for b in data.betas:
            lines += _line("beta", [b.l, b.keyee, b.eeread, b.iptype])
        if any(b.iptype == 2 for b in data.betas):
            lines += _line("npf", [data.npf, data.ptryc])
        lines += _line("local", [data.lloc, data.keyee_loc, data.eloc, data.iploctype])
        lines += _line("qfunc", [data.ifqopt, data.nqf, data.qtryc, len(data.fixes)])
        for f in data.fixes:
            lines += _line("qfix", [f.ibeta, f.rfix, f.qfix])
    return "\n".join(lines) + "\n"


class _Reader:
    def __init__(self, text: str):
        self.lines = text.replace("\r\n", "\n").split("\n")
        if self.lines and self.lines[-1] == "":
            self.lines.pop()            # the newline that ends the last record
        self.index = 0

    def read(self, key: str, count: Optional[int] = None) -> list:
        fmt, label = LINE_FORMATS[key]
        per_line = len(parse_format(fmt))
        wanted = per_line if count is None else count
        values: list = []
        while len(values) < wanted or not values and wanted == 0:
            if self.index >= len(self.lines):
                raise AdatError(self.index + 1, f"end of file while reading {label}")
            line_number, line = self.index + 1, self.lines[self.index]
            self.index += 1
            try:
                values += read_record(line, fmt, min(per_line, wanted - len(values)))
            except ValueError as exc:
                raise AdatError(line_number, f"{label}: {exc}") from None
            if wanted == 0:
                break
        return values


def parse_adat(text: str) -> Union[AEInput, PSInput, TestInput]:
    """An .adat file read exactly as rwsubs.f reads it (extra lines are ignored, as there)."""
    r = _Reader(text)
    ifae, ifpsp, ifprt, ifplw, ilogd = r.read("controls", 5)
    rlogd, emin, emax, nnt = r.read("logd")
    thresh, tol, damp, maxit = r.read("convergence")
    (title,) = r.read("title")
    controls = Controls(ifae, ifpsp, ifprt, ifplw, ilogd, rlogd, emin, emax, nnt,
                        thresh, tol, damp, maxit, title.rstrip())
    if ifae == 1:
        if ifpsp != 0:
            raise AdatError(1, "an all-electron run and a pseudopotential run in one file is not "
                               "supported here; use separate files as the Work/ library does")
        z, xion, exfact = r.read("atom")
        rmax, aasf, bbsf = r.read("mesh")
        ncspvs, irel = r.read("nstates")
        shells = [Shell(*r.read("shell")) for _ in range(max(ncspvs, 0))]
        return AEInput(controls, z, xion, exfact, rmax, aasf, bbsf, irel, shells)
    if ifae != 0 or ifpsp not in (1, 2):
        raise AdatError(1, f"ifae = {ifae}, ifpsp = {ifpsp}: not an all-electron, "
                           "generation or test input")
    ncores, nvales, nang = r.read("partition")
    if not 1 <= nang <= 4:
        raise AdatError(r.index, f"nang = {nang}; runatom accepts 1 to 4")
    bessel = r.read("bessel")
    if ifpsp == 1:
        return TestInput(controls, ncores, nvales, nang, *bessel)
    keyps, ifpcor, rinner1, rpcor = r.read("keyps")
    rinner = r.read("rinner", 2 * nang - 1) if rinner1 == 0.0 else []
    nbeta, rcloc = r.read("nbeta")
    rc = r.read("rc", nang)
    betas = [Beta(*r.read("beta")) for _ in range(max(nbeta, 0))]
    npf, ptryc = r.read("npf") if any(b.iptype == 2 for b in betas) else (0, 0.0)
    lloc, keyee_loc, eloc, iploctype = r.read("local")
    ifqopt, nqf, qtryc, nfix = r.read("qfunc")
    fixes = [QFix(*r.read("qfix")) for _ in range(max(nfix, 0))]
    return PSInput(controls, ncores, nvales, nang, *bessel, keyps, ifpcor, rinner1, rpcor,
                   rinner, rcloc, rc, betas, npf, ptryc, lloc, keyee_loc, eloc, iploctype,
                   ifqopt, nqf, qtryc, fixes)


def file_stems(case: CaseInputs) -> Dict[str, object]:
    """File stems in the library's naming (Work/000-Ex Makefiles)."""
    s = case.symbol
    return {"ae": f"{s}_ae_{case.generation_name}", "ps": f"{s}_ps",
            "tests": [(f"{s}_ae_{t.name}", f"{s}_test_{t.name}") for t in case.tests]}


def test_ae_input(case: CaseInputs, test: TestConfig) -> AEInput:
    shells = [replace(s, wwnl=test.occupations.get(s.nnlz, s.wwnl)) for s in case.ae.shells]
    return replace(case.ae, xion=test.xion, shells=shells)


def test_input(case: CaseInputs, test: TestConfig) -> TestInput:
    if test.test_input is not None:
        return test.test_input
    ps = case.ps
    controls = replace(ps.controls, ifae=0, ifpsp=1, ifplw=0)
    return TestInput(controls, ps.ncores, ps.nvales, ps.nang,
                     ps.besrmax, ps.besemin, ps.besemax, ps.besde)


def case_files(case: CaseInputs) -> Dict[str, str]:
    stems = file_stems(case)
    files = {f"{stems['ae']}.adat": format_adat(case.ae), f"{stems['ps']}.adat": format_adat(case.ps)}
    for test, (ae_stem, test_stem) in zip(case.tests, stems["tests"]):
        files[f"{ae_stem}.adat"] = format_adat(test_ae_input(case, test))
        files[f"{test_stem}.adat"] = format_adat(test_input(case, test))
    return files


def write_case(case: CaseInputs, folder) -> List[Path]:
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    for name, text in case_files(case).items():
        path = folder / name
        path.write_text(text, encoding="ascii", newline="\n")
        written.append(path)
    return written


def _makefile_vars(text: str) -> Dict[str, str]:
    return {m.group(1): m.group(2).strip()
            for m in re.finditer(r"^([A-Z][A-Z0-9_]*)\s*=\s*(.*?)\s*$", text, re.M)}


def read_case_folder(folder) -> CaseInputs:
    """A Work/ library case (Makefile + .adat files) as CaseInputs."""
    folder = Path(folder)
    variables = _makefile_vars((folder / "Makefile").read_text(encoding="latin-1"))
    symbol, generation = variables["ATOM"], variables["CFG"]
    # test configurations: CFT, CFT1.. (most cases) or CF1, CF2 (000-Ex-ca-sp-tests)
    test_keys = sorted((k for k in variables if re.fullmatch(r"CF(T\d*|\d+)", k)),
                       key=lambda k: int(re.sub(r"\D", "", k) or 0))

    def load(name: str):
        try:
            return parse_adat((folder / name).read_text(encoding="latin-1"))
        except AdatError as exc:
            raise AdatError(exc.line_number, f"{name}: {exc}") from None

    ae, ps = load(f"{symbol}_ae_{generation}.adat"), load(f"{symbol}_ps.adat")
    reference = {s.nnlz: s.wwnl for s in ae.shells}
    tests = []
    for key in test_keys:
        name = variables[key]
        if not (folder / f"{symbol}_test_{name}.adat").exists():
            continue
        test_ae = load(f"{symbol}_ae_{name}.adat")
        if [s.nnlz for s in test_ae.shells] != [s.nnlz for s in ae.shells]:
            raise ValueError(f"test {name}: its shells differ from the generation configuration's")
        changed = {s.nnlz: s.wwnl for s in test_ae.shells if s.wwnl != reference[s.nnlz]}
        tests.append(TestConfig(name, test_ae.xion, changed, load(f"{symbol}_test_{name}.adat")))
    return CaseInputs(symbol, generation, ae, ps, tests)


def makefile_text(case: CaseInputs) -> str:
    """A Makefile in the Work/000-Ex template's style, for running the case by hand."""
    tests = case.tests
    head = [f"ATOM= {case.symbol}", f"CFG= {case.generation_name}"]
    head += [f"CFT{k}= {t.name}" for k, t in enumerate(tests, 1)]
    names = ["", "ROOT= ../../..", "PROG= ${ROOT}/Bin/runatom.x", "POT=  ${ROOT}/Pot", "",
             "AEG= ${ATOM}_ae_${CFG}", "PSP= ${ATOM}_ps"]
    names += [f"AET{k}= ${{ATOM}}_ae_${{CFT{k}}}" for k in range(1, len(tests) + 1)]
    names += [f"TST{k}= ${{ATOM}}_test_${{CFT{k}}}" for k in range(1, len(tests) + 1)]
    rules = ["", "default: pseudo", "",
             "${AEG}.ae: ${AEG}.adat",
             "\t${PROG} ${AEG}.adat ${AEG}.out ${AEG}.ae ${AEG}.atwf ${AEG}.logd dummy", "",
             "pseudo: ${PSP}.uspp", "${PSP}.uspp: ${PSP}.adat ${AEG}.ae",
             "\t${PROG} ${PSP}.adat ${PSP}.out ${AEG}.ae ${PSP}.atwf ${PSP}.logd ${PSP}.uspp", ""]
    rules.append("test: " + " ".join(f"${{TST{k}}}.out" for k in range(1, len(tests) + 1)))
    for k in range(1, len(tests) + 1):
        rules += ["", f"${{AET{k}}}.ae: ${{AET{k}}}.adat",
                  f"\t${{PROG}} ${{AET{k}}}.adat ${{AET{k}}}.out ${{AET{k}}}.ae ${{AET{k}}}.atwf ${{AET{k}}}.logd dummy",
                  f"${{TST{k}}}.out: ${{TST{k}}}.adat ${{AET{k}}}.ae ${{PSP}}.uspp",
                  f"\t${{PROG}} ${{TST{k}}}.adat ${{TST{k}}}.out ${{AET{k}}}.ae ${{TST{k}}}.atwf ${{TST{k}}}.logd ${{PSP}}.uspp"]
    rules += ["", "clean:", "\t- rm *.out *.ae *.atwf *.logd *.uspp", ""]
    header = ["# Written by gqteaWinToga (Vanderbilt USPP Generator) for USPP 7.3.6.",
              "# make: generate the pseudopotential; make test: run the transferability tests."]
    return "\n".join(header + head + names + rules)
```

- [ ] **Step 4: Run them to see them pass**

Run: `python -m pytest tests/test_runatomInput.py -v`
Expected: all passed (~430 parametrised cases).

If a library file fails to round-trip, compare the reader with `rwsubs.f` lines 130–673 and fix the reader. The program read that file, so the file is not the problem. `V_ps.adat` is the only expected reader failure.

- [ ] **Step 5: Commit (if authorised)**

```bash
git add runatomInput.py tests/test_runatomInput.py
git commit -m "feat(uspp): runatom input model, .adat writer/reader and library case folders"
```

---

### Task 6: UI field mapping and pre-run checks (`runatomInput.py`, part 3)

**Files:**
- Modify: `runatomInput.py` (append)
- Test: `tests/test_runatomInput.py` (append)

**Interfaces:**
- Consumes: Task 5.
- Produces:
  - `FORM_FIELDS: Dict[str, Tuple[label, kinds, multiline, hint]]` (kinds: string of `i`/`f`, `f*` = one or more floats, `s` = text, `t` = test rows);
  - `parse_line_values(text, kinds, name) -> list`, `parse_rows(text, kinds, name) -> List[list]`;
  - `case_to_fields(case) -> Dict[str, str]`, `case_from_fields(fields) -> CaseInputs`, `ae_fields(ae)`, `ps_fields(ps)`;
  - `radial_mesh_size(z, rmax, aasf, bbsf) -> int`;
  - `check_ae(ae) -> List[str]`, `check_ps(ps, ae) -> List[str]`, `check_case(case) -> List[str]`;
  - capacity constants `MAX_MESH` … `MAX_POLY_TERMS`, `XC_FUNCTIONALS`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_runatomInput.py`:

```python
class TestFields:
    def test_typed_values_are_parsed_not_copied(self):
        assert ri.parse_line_values("1 0 3", "iii", "ifprt ifplw ilogd") == [1, 0, 3]
        assert ri.parse_line_values("1.0d-10, 1e-9  0.5 0", "fffi", "x") == [1e-10, 1e-9, 0.5, 0]
        assert ri.parse_line_values("1.8 1.8 2", "f*", "rc") == [1.8, 1.8, 2.0]

    def test_wrong_count_and_type_name_the_field(self):
        with pytest.raises(ValueError, match="ifprt ifplw ilogd: expected 3 values, got 2"):
            ri.parse_line_values("1 0", "iii", "ifprt ifplw ilogd")
        with pytest.raises(ValueError, match="is not an integer"):
            ri.parse_line_values("1.5 0 3", "iii", "ifprt ifplw ilogd")

    @pytest.mark.parametrize("folder", [p for p in CASE_DIRS if not (p / "V_ps.adat").exists()],
                             ids=lambda p: p.name)
    def test_fields_round_trip_every_library_case(self, folder):
        case = ri.read_case_folder(folder)
        expected = ri.replace(
            case, ps=ri.replace(case.ps, controls=ri.replace(case.ps.controls, title=case.ae.controls.title)),
            tests=[ri.replace(t, test_input=None) for t in case.tests])
        assert ri.case_from_fields(ri.case_to_fields(case)) == expected


class TestChecks:
    def test_mesh_size_matches_runatom(self):
        # Ti example: rinit printed "value of mesh generated in rinit is   851"
        assert ri.radial_mesh_size(22.0, 200.0, 6.0, 59.0) == 851

    @pytest.mark.parametrize("folder", [p for p in CASE_DIRS if not (p / "V_ps.adat").exists()],
                             ids=lambda p: p.name)
    def test_no_library_case_is_refused(self, folder):
        # runatom accepted every one of these, so a refusal is a wrong rule
        assert ri.check_case(ri.read_case_folder(folder)) == []

    def _ti(self, **ps_changes):
        case = _ti_case()
        return ri.replace(case, ps=ri.replace(case.ps, **ps_changes))

    def test_gga_requires_ifqopt_3_even_for_minus_one(self):
        case = _ti_case()
        gga = ri.replace(case, ae=ri.replace(case.ae, exfact=5.0), ps=ri.replace(case.ps, ifqopt=-1))
        assert any("ifqopt = 3" in p for p in ri.check_case(gga))

    def test_occupations_must_add_up(self):
        case = _ti_case()
        bad = ri.replace(case, ae=ri.replace(case.ae, xion=0.0))
        assert any("add up to" in p for p in ri.check_case(bad))

    def test_projector_rules(self):
        betas = _ti_case().ps.betas
        assert any("non-decreasing" in p for p in ri.check_case(self._ti(betas=list(reversed(betas)))))
        assert any("keyee < 0" in p for p in ri.check_case(self._ti(betas=[ri.replace(betas[0], keyee=-1)] + betas[1:])))
        assert any("npf" in p for p in ri.check_case(self._ti(npf=4)))

    def test_mesh_limit(self):
        case = _ti_case()
        assert any("radial mesh" in p for p in ri.check_case(ri.replace(case, ae=ri.replace(case.ae, rmax=1e6))))

    def test_long_file_names_are_refused(self):
        case = ri.replace(_ti_case(), generation_name="a" * 40)
        assert any("40 characters" in p for p in ri.check_case(case))

    def test_unknown_shell_in_a_test(self):
        case = _ti_case()
        case.tests[0].occupations[430] = 1.0
        assert any("430" in p for p in ri.check_case(case))
```

- [ ] **Step 2: Run them to see them fail**

Run: `python -m pytest tests/test_runatomInput.py -v -k "Fields or Checks"`
Expected: FAIL with `AttributeError: ... 'parse_line_values'`.

- [ ] **Step 3: Implement**

Append to `runatomInput.py`:

```python
# Capacities compiled into runatom (runatom.f parameters; initsubs.f for the mesh)
MAX_MESH = 1000
MAX_SHELLS = 26
MAX_BETAS = 10
MAX_VALENCE = 5
MAX_NANG = 4
MAX_BETAS_PER_L = 4
MAX_POLY_TERMS = 20
TITLE_WIDTH = 20
OCCUPATION_TOLERANCE = 1e-3
XC_FUNCTIONALS = {0.0: "Ceperley-Alder (LDA)", -1.0: "Wigner (LDA)", -2.0: "Hedin-Lundqvist (LDA)",
                  -3.0: "Gunnarsson-Lundqvist (LDA)", 1.0: "CA + Becke88 + LYP (GGA)",
                  2.0: "CA + Becke88 (GGA)", 3.0: "CA + Becke88 + Perdew86 (GGA)",
                  4.0: "Perdew-Wang 91 (GGA)", 5.0: "PBE (GGA)"}

FORM_FIELDS = {
    "symbol": ("Element symbol (file names)", "s", False, "e.g. ti gives ti_ae_<configuration>.adat and ti_ps.adat"),
    "generation": ("Generation configuration", "s", False, "a short name such as s2d1"),
    "title": ("Title (at most 20 characters)", "s", False, ""),
    "ae_print": ("ifprt  ifplw  ilogd", "iii", False, "printout level, wavefunction plot file (0/1), number of l for log derivatives"),
    "ae_logd": ("rlogd  emin  emax  nnt", "fffi", False, "log-derivative radius (bohr) and energy mesh (Ry)"),
    "ae_conv": ("thresh  tol  damp  maxit", "fffi", False, "0 = program default (1e-6, 1e-5, 0.5, 250)"),
    "ae_atom": ("z  xion  exfact", "fff", False, "exfact: 0 CA, -1 Wigner, -2 HL, -3 GL, 1 BLYP, 2 B88, 3 BP86, 4 PW91, 5 PBE"),
    "ae_mesh": ("rmax  aasf  bbsf  irel", "fffi", False, "radial mesh (0 = 80, 6, 40); irel 0 non-relativistic, 2 scalar-relativistic"),
    "ae_shells": ("Shells: nnlz  wwnl  ee", "iff", True, "one shell per line, core first, valence last; e.g. 300 2.0 -3.0"),
    "ps_print": ("ifprt  ifplw  ilogd", "iii", False, ""),
    "ps_logd": ("rlogd  emin  emax  nnt", "fffi", False, "rinner must not exceed rlogd"),
    "ps_conv": ("thresh  tol  damp  maxit", "fffi", False, ""),
    "ps_partition": ("ncores  nvales  nang", "iii", False, "ncores + nvales = number of shells; nang = lmax + 1 of the projectors"),
    "ps_bessel": ("besrmax  besemin  besemax  besde", "ffff", False, "ghost check in a Bessel basis: box radius (bohr), cutoffs (Ry)"),
    "ps_keyps": ("keyps  ifpcor  rinner  rpcor", "iiff", False, "keyps 3; ifpcor 1 = nonlinear core correction; rinner 0 = per-L values below"),
    "ps_rinner": ("rinner per L (if rinner = 0)", "f*", False, "2*nang - 1 values; ignored when rinner is not 0"),
    "ps_rcloc": ("rcloc", "f", False, "radius for the local potential"),
    "ps_rc": ("rc per l", "f*", False, "nang values (s, p, d, f)"),
    "ps_betas": ("Projectors: l  keyee  eeread  iptype", "iifi", True, "one per line, l non-decreasing; keyee > 0 takes that valence eigenvalue"),
    "ps_npf": ("npf  ptryc", "if", False, "written only when a projector has iptype 2"),
    "ps_local": ("lloc  keyee  eloc  iploctype", "iifi", False, "lloc -1 = polynomial inside rcloc"),
    "ps_qfunc": ("ifqopt  nqf  qtryc", "iif", False, "GGA functionals need ifqopt 3"),
    "ps_fixes": ("Q constraints: ibfix  rfix  qfix", "iff", True, "optional; ifqopt 2 or 3 only"),
    "tests": ("Tests: name  xion  nnlz=occupation ...", "t", True, "e.g. s2d2 0.0 320=2  (unlisted shells keep the generation occupation)"),
}


def parse_line_values(text: str, kinds: str, name: str) -> list:
    """Space- or comma-separated values typed by a user, typed by `kinds`."""
    tokens = [t for t in re.split(r"[\s,]+", (text or "").strip()) if t]
    if kinds.endswith("*"):
        if not tokens:
            raise ValueError(f"{name}: enter at least one value")
        kinds = kinds[0] * len(tokens)
    if len(tokens) != len(kinds):
        raise ValueError(f"{name}: expected {len(kinds)} values, got {len(tokens)}")
    values = []
    for token, kind in zip(tokens, kinds):
        try:
            values.append(int(token) if kind == "i" else float(token.lower().replace("d", "e")))
        except ValueError:
            what = "an integer" if kind == "i" else "a number"
            raise ValueError(f"{name}: {token!r} is not {what}") from None
    return values


def parse_rows(text: str, kinds: str, name: str) -> List[list]:
    rows = [line for line in (text or "").splitlines() if line.strip()]
    return [parse_line_values(line, kinds, f"{name}, row {k}") for k, line in enumerate(rows, 1)]


def _number(value: float) -> str:
    text = repr(float(value))
    return text[:-2] if text.endswith(".0") else text


def _text(values) -> str:
    return "  ".join(str(v) if isinstance(v, int) else _number(v) for v in values)


def _controls_fields(prefix: str, c: Controls) -> Dict[str, str]:
    return {f"{prefix}_print": _text([c.ifprt, c.ifplw, c.ilogd]),
            f"{prefix}_logd": _text([c.rlogd, c.emin, c.emax, c.nnt]),
            f"{prefix}_conv": _text([c.thresh, c.tol, c.damp, c.maxit])}


def ae_fields(ae: AEInput) -> Dict[str, str]:
    fields = _controls_fields("ae", ae.controls)
    fields.update(ae_atom=_text([ae.z, ae.xion, ae.exfact]),
                  ae_mesh=_text([ae.rmax, ae.aasf, ae.bbsf, ae.irel]),
                  ae_shells="\n".join(_text([s.nnlz, s.wwnl, s.ee]) for s in ae.shells))
    return fields


def ps_fields(ps: PSInput) -> Dict[str, str]:
    fields = _controls_fields("ps", ps.controls)
    fields.update(
        ps_partition=_text([ps.ncores, ps.nvales, ps.nang]),
        ps_bessel=_text([ps.besrmax, ps.besemin, ps.besemax, ps.besde]),
        ps_keyps=_text([ps.keyps, ps.ifpcor, ps.rinner1, ps.rpcor]),
        ps_rinner=_text(ps.rinner), ps_rcloc=_text([ps.rcloc]), ps_rc=_text(ps.rc),
        ps_betas="\n".join(_text([b.l, b.keyee, b.eeread, b.iptype]) for b in ps.betas),
        ps_npf=_text([ps.npf, ps.ptryc]) if any(b.iptype == 2 for b in ps.betas) else "",
        ps_local=_text([ps.lloc, ps.keyee_loc, ps.eloc, ps.iploctype]),
        ps_qfunc=_text([ps.ifqopt, ps.nqf, ps.qtryc]),
        ps_fixes="\n".join(_text([f.ibeta, f.rfix, f.qfix]) for f in ps.fixes))
    return fields


def case_to_fields(case: CaseInputs) -> Dict[str, str]:
    fields = {"symbol": case.symbol, "generation": case.generation_name, "title": case.ae.controls.title}
    fields.update(ae_fields(case.ae))
    fields.update(ps_fields(case.ps))
    fields["tests"] = "\n".join(
        "  ".join([t.name, _number(t.xion)] + [f"{n}={_number(w)}" for n, w in t.occupations.items()])
        for t in case.tests)
    return fields


def _values(fields: Dict[str, str], key: str):
    label, kinds, multiline, _ = FORM_FIELDS[key]
    text = fields.get(key, "")
    return parse_rows(text, kinds, label) if multiline else parse_line_values(text, kinds, label)


def _controls_from(fields: Dict[str, str], prefix: str, ifae: int, ifpsp: int, title: str) -> Controls:
    ifprt, ifplw, ilogd = _values(fields, f"{prefix}_print")
    rlogd, emin, emax, nnt = _values(fields, f"{prefix}_logd")
    thresh, tol, damp, maxit = _values(fields, f"{prefix}_conv")
    return Controls(ifae, ifpsp, ifprt, ifplw, ilogd, rlogd, emin, emax, nnt, thresh, tol, damp, maxit, title)


def _parse_tests(text: str) -> List[TestConfig]:
    tests = []
    for k, line in enumerate((l for l in (text or "").splitlines() if l.strip()), 1):
        tokens = line.split()
        if len(tokens) < 2:
            raise ValueError(f"Tests, row {k}: write a name and xion, then nnlz=occupation pairs")
        (xion,) = parse_line_values(tokens[1], "f", f"Tests, row {k}, xion")
        occupations = {}
        for pair in tokens[2:]:
            nnlz, _, value = pair.partition("=")
            if not value:
                raise ValueError(f"Tests, row {k}: {pair!r} is not nnlz=occupation")
            (n,) = parse_line_values(nnlz, "i", f"Tests, row {k}, {pair}")
            (w,) = parse_line_values(value, "f", f"Tests, row {k}, {pair}")
            occupations[n] = w
        tests.append(TestConfig(tokens[0], xion, occupations))
    return tests


def case_from_fields(fields: Dict[str, str]) -> CaseInputs:
    """The window's text fields as CaseInputs; ValueError names the field at fault."""
    title = (fields.get("title") or "").rstrip()      # as parse_adat keeps it (a20, trailing blanks)
    z, xion, exfact = _values(fields, "ae_atom")
    rmax, aasf, bbsf, irel = _values(fields, "ae_mesh")
    shells = [Shell(*row) for row in _values(fields, "ae_shells")]
    ae = AEInput(_controls_from(fields, "ae", 1, 0, title), z, xion, exfact, rmax, aasf, bbsf, irel, shells)
    ncores, nvales, nang = _values(fields, "ps_partition")
    keyps, ifpcor, rinner1, rpcor = _values(fields, "ps_keyps")
    betas = [Beta(*row) for row in _values(fields, "ps_betas")]
    npf, ptryc = _values(fields, "ps_npf") if any(b.iptype == 2 for b in betas) else (0, 0.0)
    lloc, keyee_loc, eloc, iploctype = _values(fields, "ps_local")
    ifqopt, nqf, qtryc = _values(fields, "ps_qfunc")
    ps = PSInput(_controls_from(fields, "ps", 0, 2, title), ncores, nvales, nang, *_values(fields, "ps_bessel"),
                 keyps, ifpcor, rinner1, rpcor,
                 _values(fields, "ps_rinner") if rinner1 == 0.0 else [],
                 _values(fields, "ps_rcloc")[0], _values(fields, "ps_rc"), betas, npf, ptryc,
                 lloc, keyee_loc, eloc, iploctype, ifqopt, nqf, qtryc,
                 [QFix(*row) for row in _values(fields, "ps_fixes")])
    return CaseInputs((fields.get("symbol") or "").strip(), (fields.get("generation") or "").strip(),
                      ae, ps, _parse_tests(fields.get("tests", "")))


def radial_mesh_size(z: float, rmax: float, aasf: float, bbsf: float) -> int:
    """Number of radial points runatom will use (initsubs.f, subroutine rinit)."""
    rmax, aasf, bbsf = rmax or 80.0, aasf or 6.0, bbsf or 40.0
    a = math.exp(-aasf) / z
    mesh = 2 + int(math.log(rmax / a + 1.0) / (1.0 / bbsf))
    return mesh + 1 if mesh % 2 == 0 else mesh


def check_ae(ae: AEInput) -> List[str]:
    """runatom's own all-electron checks (rwsubs.f readin, initsubs.f rinit) plus this tool's limits."""
    problems = []
    if len(ae.controls.title) > TITLE_WIDTH or not ae.controls.title.isascii():
        problems.append(f"The title must be at most {TITLE_WIDTH} ASCII characters.")
    if ae.z <= 0:
        problems.append("z must be positive.")
    if ae.exfact not in XC_FUNCTIONALS:
        problems.append(f"exfact = {ae.exfact:g} is not one of {sorted(XC_FUNCTIONALS)}.")
    if ae.irel not in (0, 2):
        problems.append("irel must be 0 or 2: runatom generates and tests pseudopotentials "
                        "only without the Dirac equation (irel = 1).")
    if not 1 <= len(ae.shells) <= MAX_SHELLS:
        problems.append(f"Give between 1 and {MAX_SHELLS} shells (runatom's limit).")
    for s in ae.shells:
        n, l = s.nnlz // 100, (s.nnlz // 10) % 10
        if not (100 <= s.nnlz <= 999 and s.nnlz % 10 == 0 and l < n and l <= 3):
            problems.append(f"Shell {s.nnlz}: write n, l and 0 as three digits (320 = 3d), with l < n and l <= 3.")
    electrons = sum(max(s.wwnl, 0.0) for s in ae.shells)
    if abs(ae.z - electrons - ae.xion) > OCCUPATION_TOLERANCE:
        problems.append(f"The occupations add up to {electrons:g} electrons, but z - xion = {ae.z - ae.xion:g}.")
    if ae.z > 0:
        mesh = radial_mesh_size(ae.z, ae.rmax, ae.aasf, ae.bbsf)
        if mesh > MAX_MESH:
            problems.append(f"The radial mesh would have {mesh} points; runatom allows {MAX_MESH}. Lower rmax or bbsf.")
    c = ae.controls
    if min(c.thresh, c.tol, c.damp) < 0 or c.maxit < 0:
        problems.append("thresh, tol, damp and maxit cannot be negative.")
    return problems


def check_ps(ps: PSInput, ae: AEInput) -> List[str]:
    """runatom's generation checks (rwsubs.f readin, lines 300-673) plus this tool's limits."""
    problems = []
    if ps.keyps != 3:
        problems.append("keyps must be 3 (the only scheme runatom implements).")
    if ps.ncores + ps.nvales != len(ae.shells):
        problems.append(f"ncores + nvales = {ps.ncores + ps.nvales}, but the all-electron input has {len(ae.shells)} shells.")
    if not 1 <= ps.nvales <= MAX_VALENCE:
        problems.append(f"nvales must be between 1 and {MAX_VALENCE}.")
    if not 1 <= ps.nang <= MAX_NANG:
        return problems + [f"nang must be between 1 and {MAX_NANG}."]
    n_l = 2 * ps.nang - 1
    rinner = ps.rinner if ps.rinner1 == 0.0 else [ps.rinner1] * n_l
    if ps.rinner1 == 0.0 and len(ps.rinner) != n_l:
        problems.append(f"rinner = 0 needs {n_l} per-L values.")
    if rinner and max(rinner) > ps.controls.rlogd:
        problems.append(f"rinner ({max(rinner):g}) must not exceed rlogd ({ps.controls.rlogd:g}).")
    if len(ps.rc) != ps.nang:
        problems.append(f"Give {ps.nang} rc values, one per l.")
    if not 1 <= len(ps.betas) <= MAX_BETAS:
        problems.append(f"Give between 1 and {MAX_BETAS} projectors.")
    ls = [b.l for b in ps.betas]
    if ls != sorted(ls):
        problems.append("List the projectors with l non-decreasing (s, then p, then d).")
    for l in sorted(set(ls)):
        if ls.count(l) > MAX_BETAS_PER_L:
            problems.append(f"At most {MAX_BETAS_PER_L} projectors per l (l = {l} has {ls.count(l)}).")
    if ls and max(ls) + 1 != ps.nang:
        problems.append(f"nang = {ps.nang}, but the highest projector l is {max(ls)}: nang must be lmax + 1.")
    for k, b in enumerate(ps.betas, 1):
        if b.keyee < 0:
            problems.append(f"Projector {k}: keyee < 0 (wavefunctions from another reference "
                            "configuration) is not supported by this tool.")
        elif b.keyee > ps.nvales:
            problems.append(f"Projector {k}: keyee = {b.keyee} exceeds nvales = {ps.nvales}.")
        if not 0 <= b.iptype <= 4:
            problems.append(f"Projector {k}: iptype must be 0 to 4.")
    smooth = any(b.iptype == 2 for b in ps.betas)
    if smooth and not 4 < ps.npf <= MAX_POLY_TERMS:
        problems.append(f"npf must be between 5 and {MAX_POLY_TERMS} when a projector has iptype 2.")
    if ps.keyee_loc > ps.nvales:
        problems.append(f"Local potential: keyee = {ps.keyee_loc} exceeds nvales = {ps.nvales}.")
    if ps.keyee_loc > 0 and ps.iploctype == 2 and not smooth:
        problems.append("Local potential: iploctype 2 needs npf and ptryc, i.e. a projector with iptype 2.")
    if ps.ifqopt > 3:
        problems.append("ifqopt must be -1 (no pseudization), 0, 1, 2 or 3.")
    if 1 <= ps.ifqopt <= 3 and not 3 < ps.nqf <= MAX_POLY_TERMS:
        problems.append(f"ifqopt {ps.ifqopt} needs 3 < nqf <= {MAX_POLY_TERMS}.")
    if ps.fixes and ps.ifqopt < 2:
        problems.append("Q-function constraints need ifqopt 2 or 3.")
    if len(ps.fixes) > len(ps.betas):
        problems.append("There are more Q-function constraints than projectors.")
    previous = -1
    for k, f in enumerate(ps.fixes, 1):
        if f.ibeta <= previous or f.ibeta > len(ps.betas):
            problems.append(f"Constraint {k}: ibfix must increase and be at most nbeta.")
        previous = f.ibeta
        if rinner and f.rfix > min(rinner):
            problems.append(f"Constraint {k}: rfix must not exceed rinner.")
    if ae.exfact >= 1 and ps.ifqopt != 3:
        problems.append("GGA functionals (exfact >= 1) need ifqopt = 3: runatom refuses -1, 0, 1 and 2.")
    return problems


def check_case(case: CaseInputs) -> List[str]:
    """Every problem runatom would stop on, found before it runs; [] when ready."""
    problems = []
    names = [case.symbol, case.generation_name] + [t.name for t in case.tests]
    for name in names:
        if not SAFE_NAME.match(name or ""):
            problems.append(f"Name {name!r}: use only letters, digits, '.', '_' and '-'.")
    test_names = [t.name for t in case.tests]
    if len(set(test_names)) != len(test_names) or case.generation_name in test_names:
        problems.append("Test names must be different from each other and from the generation configuration.")
    if problems:
        return problems
    stems = file_stems(case)
    longest = max([stems["ae"], stems["ps"]] + [s for pair in stems["tests"] for s in pair], key=len)
    if len(longest) + len(".atwf") > MAX_FORTRAN_NAME:
        return [f"The file name {longest}.atwf is longer than {MAX_FORTRAN_NAME} characters, which "
                "runatom cannot read; shorten the symbol or the configuration names."]
    try:      # check the values runatom will read, after the column-exact round trip
        ae = parse_adat(format_adat(case.ae))
        ps = parse_adat(format_adat(case.ps))
    except ValueError as exc:
        return [str(exc)]
    problems += [f"All-electron: {p}" for p in check_ae(ae)]
    problems += [f"Pseudopotential: {p}" for p in check_ps(ps, ae)]
    known = {s.nnlz for s in ae.shells}
    for test in case.tests:
        unknown = sorted(set(test.occupations) - known)
        if unknown:
            problems.append(f"Test {test.name}: shells {unknown} are not in the all-electron shell list.")
            continue
        try:
            test_ae = parse_adat(format_adat(test_ae_input(case, test)))
        except ValueError as exc:
            problems.append(f"Test {test.name}: {exc}")
            continue
        problems += [f"Test {test.name}: {p}" for p in check_ae(test_ae)]
    return problems
```

- [ ] **Step 4: Run them to see them pass**

Run: `python -m pytest tests/test_runatomInput.py -v`
Expected: all passed.

`test_no_library_case_is_refused` is the guard against over-strict rules. runatom accepted every one of those cases, so a failure means a rule here differs from `rwsubs.f`. Re-read the matching lines and fix the rule, not the test.

- [ ] **Step 5: Commit (if authorised)**

```bash
git add runatomInput.py tests/test_runatomInput.py
git commit -m "feat(uspp): runatom form fields and pre-run checks mirroring rwsubs.f"
```

---

### Task 7: Output parsing and the `.uspp` reader (`runatomOutput.py`)

**Files:**
- Create: `runatomOutput.py`
- Test: `tests/test_runatomOutput.py`

**Interfaces:**
- Consumes: `tests/runatom_fixtures` (reference `.out`/`.uspp`).
- Produces:
  - `EigenRow(nlm, all_electron, pseudo, diff)` with `.l`;
  - `OutReport(errors, warnings, comparison, bessel, total_energy, zv)`, `parse_out(text) -> OutReport`;
  - `ghost_candidates(report, tolerance=GHOST_TOLERANCE_RY) -> List[int]`;
  - `uspp_records(data: bytes) -> List[bytes]`, `UsppData`, `read_uspp(path) -> UsppData`;
  - `job_summary(stages: List[Tuple[str, float, OutReport]], uspp: Optional[UsppData]) -> str`;
  - constants `EIGENVALUE_FLAG_RY = 1e-3`, `GHOST_TOLERANCE_RY = 0.1`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_runatomOutput.py`:

```python
"""runatomOutput: the .out report and the .uspp file of runatom."""

import struct
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import runatomOutput as ro
from runatom_fixtures import REFERENCE_CASES, case_info, reference_bytes, reference_text

SRC = Path(__file__).resolve().parents[1]
TI = "000-Ex-ca-sp-test"


def test_module_does_not_import_toga():
    code = "import sys, runatomOutput; assert 'toga' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], cwd=SRC, check=True)


class TestParseOut:
    def test_ti_generation_report(self):
        report = ro.parse_out(reference_text(TI, "ti_ps.out"))
        assert [r.nlm for r in report.comparison] == [300, 310, 400, 320]
        assert report.comparison[0].all_electron == pytest.approx(-5.530317, abs=1e-6)
        assert [r.l for r in report.comparison] == [0, 1, 0, 2]
        assert sorted(report.bessel) == [0, 1, 2] and len(report.bessel[0]) == 4
        assert report.zv == 12.0
        assert report.errors == []
        assert any("max q eigval" in w for w in report.warnings)

    def test_ti_transferability_test(self):
        report = ro.parse_out(reference_text(TI, "ti_test_s2d2.out"))
        assert report.comparison[1].diff == pytest.approx(-0.000707, abs=1e-6)
        assert report.total_energy is not None

    def test_program_error_is_collected_with_its_detail(self):
        text = "something\n ***error in subroutine readin\n nang =  5 but idim5 =  4\n"
        assert ro.parse_out(text).errors == ["***error in subroutine readin nang =  5 but idim5 =  4"]


class TestGhosts:
    def test_no_ghost_in_ti(self):
        assert ro.ghost_candidates(ro.parse_out(reference_text(TI, "ti_ps.out"))) == []

    def test_extra_bound_state_is_a_ghost(self):
        report = ro.parse_out(reference_text(TI, "ti_ps.out"))
        report.bessel[2] = [-3.0] + report.bessel[2]     # one more d state than the AE atom has
        assert ro.ghost_candidates(report) == [2]


class TestReadUspp:
    @pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
    def test_every_reference_pseudopotential_reads(self, name, tmp_path):
        atom = case_info(name)["atom"]
        path = tmp_path / "p.uspp"
        path.write_bytes(reference_bytes(name, f"{atom}_ps.uspp"))
        data = ro.read_uspp(path)
        assert data.version == (7, 3, 6)
        assert data.beta.shape == (data.nbeta, data.kkbeta)
        assert data.qfunc.shape == (data.nbeta, data.nbeta, data.kkbeta)
        assert np.allclose(data.qqq, data.qqq.T)
        assert data.r.shape == (data.mesh,) and data.r[0] == 0.0

    def test_ti_header(self, tmp_path):
        path = tmp_path / "ti.uspp"
        path.write_bytes(reference_bytes(TI, "ti_ps.uspp"))
        data = ro.read_uspp(path)
        assert (data.title, data.zv, data.mesh, data.nbeta, data.kkbeta) == ("titanium", 12.0, 851, 6, 585)
        assert list(data.l) == [0, 0, 1, 1, 2, 2] and list(data.iptype) == [2] * 6

    def test_zero_filled_file_is_refused(self, tmp_path):
        path = tmp_path / "zeros.uspp"           # like the upstream Nb file
        path.write_bytes(bytes(10240))
        with pytest.raises(ValueError, match="record 1"):
            ro.read_uspp(path)

    def test_mismatched_record_markers_are_refused(self):
        with pytest.raises(ValueError, match="not a gfortran unformatted file"):
            ro.uspp_records(struct.pack("<i", 4) + b"abcd" + struct.pack("<i", 5))


def test_summary_flags_large_test_errors_and_reports_the_file(tmp_path):
    path = tmp_path / "ti.uspp"
    path.write_bytes(reference_bytes(TI, "ti_ps.uspp"))
    text = ro.job_summary([("Pseudopotential generation", 0.4, ro.parse_out(reference_text(TI, "ti_ps.out"))),
                           ("Transferability test (s2d2)", 0.2, ro.parse_out(reference_text(TI, "ti_test_s2d2.out")))],
                          ro.read_uspp(path))
    assert "no extra bound state" in text and "zv = 12" in text
    assert "more than 1 mRy" not in text      # every Ti test difference is below 1 mRy
```

- [ ] **Step 2: Run them to see them fail**

Run: `python -m pytest tests/test_runatomOutput.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'runatomOutput'`.

- [ ] **Step 3: Implement**

Create `runatomOutput.py`:

```python
"""Results of the Vanderbilt USPP program runatom: the .out report and the .uspp file.

No Toga here (a subprocess test pins it). read_uspp follows Utils/read_uspp.f
(files of version 7.3.2 or later): gfortran sequential unformatted records, each
framed by its 4-byte little-endian length before and after; integers are 4 bytes.
"""

import re
import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

EIGENVALUE_FLAG_RY = 1e-3      # Doc/TUTORIAL: a test is "probably OK" below 1 mRy
GHOST_TOLERANCE_RY = 0.1       # above the highest AE eigenvalue of a channel

_ROW = re.compile(r"^\s*(\d{3})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
_LROW = re.compile(r"^\s*l=\s*(\d+)\s+(.*)$")
_TOTAL = re.compile(r"total energy =\s*(-?\d+\.\d+)")
_ZV = re.compile(r"\|\s*z =\s*[-\d.]+\s+zv =\s*(-?\d+\.\d+)")   # the boxed pseudopotential report


@dataclass(frozen=True)
class EigenRow:
    nlm: int
    all_electron: float
    pseudo: float
    diff: float

    @property
    def l(self) -> int:
        return (self.nlm // 10) % 10


@dataclass
class OutReport:
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    comparison: List[EigenRow] = field(default_factory=list)   # last table in the file
    bessel: Dict[int, List[float]] = field(default_factory=dict)  # last US-Hamiltonian block
    total_energy: Optional[float] = None
    zv: Optional[float] = None


def _last_index(lines: Sequence[str], marker: str) -> Optional[int]:
    hits = [i for i, line in enumerate(lines) if marker in line]
    return hits[-1] if hits else None


def parse_out(text: str) -> OutReport:
    lines = text.splitlines()
    report = OutReport()
    for i, line in enumerate(lines):
        low = line.lower()
        if "***error" in low or "*** error" in low:
            report.errors.append(" ".join(l.strip() for l in lines[i:i + 3] if l.strip()))
        elif "warning" in low and line.strip() not in report.warnings:
            report.warnings.append(line.strip())
    start = _last_index(lines, "comparison of all-electron and pseudo eigenvalues")
    if start is not None:
        for line in lines[start + 1:]:
            match = _ROW.match(line)
            if match:
                report.comparison.append(EigenRow(int(match.group(1)), *map(float, match.group(2, 3, 4))))
            elif report.comparison and line.strip():
                break
    start = _last_index(lines, "eigenvalues of US hamiltonian")
    if start is not None:
        for line in lines[start + 1:]:
            match = _LROW.match(line)
            if not match:
                break
            report.bessel[int(match.group(1))] = [float(v) for v in match.group(2).split()]
    totals = _TOTAL.findall(text)
    report.total_energy = float(totals[-1]) if totals else None
    zv = _ZV.search(text)
    report.zv = float(zv.group(1)) if zv else None
    return report


def ghost_candidates(report: OutReport, tolerance: float = GHOST_TOLERANCE_RY) -> List[int]:
    """Channels with more bound Bessel-basis states than all-electron valence states (Doc/TUTORIAL)."""
    ghosts = []
    for l in sorted({row.l for row in report.comparison}):
        expected = [row.all_electron for row in report.comparison if row.l == l]
        found = report.bessel.get(l, [])
        if sum(1 for e in found if e < max(expected) + tolerance) > len(expected):
            ghosts.append(l)
    return ghosts


def uspp_records(data: bytes) -> List[bytes]:
    records, position = [], 0
    while position < len(data):
        index = len(records) + 1
        if position + 4 > len(data):
            raise ValueError(f"record {index}: truncated length marker")
        (length,) = struct.unpack_from("<i", data, position)
        end = position + 4 + length
        if length < 0 or end + 4 > len(data):
            raise ValueError(f"record {index}: length {length} runs past the end of the file")
        (trailer,) = struct.unpack_from("<i", data, end)
        if trailer != length:
            raise ValueError(f"record {index}: leading and trailing lengths differ ({length} vs {trailer}); "
                             "not a gfortran unformatted file")
        records.append(data[position + 4:end])
        position = end + 4
    return records


class _Record:
    def __init__(self, payload: bytes, index: int):
        self.payload, self.index, self.position = payload, index, 0

    def _take(self, dtype: str, count: int) -> np.ndarray:
        size = np.dtype(dtype).itemsize * count
        if self.position + size > len(self.payload):
            raise ValueError(f"record {self.index} is shorter than the USPP 7.3.6 layout expects")
        values = np.frombuffer(self.payload, dtype=dtype, count=count, offset=self.position).copy()
        self.position += size
        return values

    def ints(self, count: int) -> np.ndarray:
        return self._take("<i4", count)

    def reals(self, count: int) -> np.ndarray:
        return self._take("<f8", count)

    def int(self) -> int:
        return int(self.ints(1)[0])

    def real(self) -> float:
        return float(self.reals(1)[0])

    def text(self, count: int) -> str:
        if self.position + count > len(self.payload):
            raise ValueError(f"record {self.index} is shorter than the USPP 7.3.6 layout expects")
        value = self.payload[self.position:self.position + count].decode("latin-1")
        self.position += count
        return value

    def end(self):
        if self.position != len(self.payload):
            raise ValueError(f"record {self.index} is longer than the USPP 7.3.6 layout expects")


@dataclass
class UsppData:
    version: Tuple[int, int, int]
    date: Tuple[int, int, int]
    title: str
    z: float
    zv: float
    exfact: float
    valence: List[Tuple[int, float, float]]      # (nnlz, occupation, eigenvalue Ry)
    mesh: int
    etot: float
    keyps: int
    ifpcor: int
    nang: int
    lloc: int
    eloc: float
    ifqopt: int
    nqf: int
    qtryc: float
    rinner: np.ndarray
    irel: int
    rc: np.ndarray
    kkbeta: int
    l: np.ndarray
    eee: np.ndarray
    beta: np.ndarray            # (nbeta, kkbeta)
    ddd0: np.ndarray            # (nbeta, nbeta), symmetric
    ddd: np.ndarray
    qqq: np.ndarray
    qfunc: np.ndarray           # (nbeta, nbeta, kkbeta)
    qfcoef: np.ndarray          # (nbeta, nbeta, 2*nang-1, nqf)
    iptype: np.ndarray
    npf: int
    ptryc: float
    rcloc: float
    vloc0: np.ndarray
    rpcor: Optional[float]
    rspsco: Optional[np.ndarray]
    vloc: np.ndarray
    rsatom: np.ndarray
    r: np.ndarray
    rab: np.ndarray
    wf: np.ndarray              # (nchi, mesh)

    @property
    def nbeta(self) -> int:
        return len(self.l)


def read_uspp(path) -> UsppData:
    records = [_Record(p, i) for i, p in enumerate(uspp_records(Path(path).read_bytes()), 1)]
    stream = iter(records)

    def record() -> _Record:
        try:
            return next(stream)
        except StopIteration:
            raise ValueError("the file ends before the USPP 7.3.6 layout does") from None

    r = record(); version = tuple(int(v) for v in r.ints(3)); date = tuple(int(v) for v in r.ints(3)); r.end()
    if version < (7, 3, 2):
        raise ValueError(f"file version {version} predates the 7.3.2 layout this reader follows")
    r = record(); title = r.text(20).strip(); z, zv, exfact = r.real(), r.real(), r.real()
    nval, mesh = r.int(), r.int(); etot = r.real(); r.end()
    r = record(); valence = [(r.int(), r.real(), r.real()) for _ in range(nval)]; r.end()
    r = record(); keyps, ifpcor = r.int(), r.int(); r.real(); r.end()
    r = record(); nang, lloc, eloc, ifqopt, nqf, qtryc = r.int(), r.int(), r.real(), r.int(), r.int(), r.real(); r.end()
    n_l = 2 * nang - 1
    r = record(); rinner = r.reals(n_l); r.end()
    r = record(); irel = r.int(); r.end()
    r = record(); rc = r.reals(nang); r.end()
    r = record(); nbeta, kkbeta = r.int(), r.int(); r.end()
    l, eee = np.zeros(nbeta, int), np.zeros(nbeta)
    beta = np.zeros((nbeta, kkbeta))
    ddd0, ddd, qqq = (np.zeros((nbeta, nbeta)) for _ in range(3))
    qfunc = np.zeros((nbeta, nbeta, kkbeta))
    qfcoef = np.zeros((nbeta, nbeta, n_l, nqf))
    for j in range(nbeta):
        r = record(); l[j], eee[j] = r.int(), r.real(); beta[j] = r.reals(kkbeta); r.end()
        for k in range(j, nbeta):
            r = record()
            ddd0[j, k] = ddd0[k, j] = r.real()
            ddd[j, k] = ddd[k, j] = r.real()
            qqq[j, k] = qqq[k, j] = r.real()
            qfunc[j, k] = qfunc[k, j] = r.reals(kkbeta)
            qfcoef[j, k] = qfcoef[k, j] = r.reals(nqf * n_l).reshape(n_l, nqf)
            r.end()
    r = record(); iptype = r.ints(nbeta); npf, ptryc = r.int(), r.real(); r.end()
    r = record(); rcloc = r.real(); vloc0 = r.reals(mesh); r.end()
    rpcor = rspsco = None
    if ifpcor > 0:
        r = record(); rpcor = r.real(); r.end()
        r = record(); rspsco = r.reals(mesh); r.end()
    arrays = []
    for _ in range(4):                       # vloc, rsatom, r, rab
        r = record(); arrays.append(r.reals(mesh)); r.end()
    r = record(); nchi = r.int(); r.end()
    r = record(); wf = r.reals(mesh * nchi).reshape(nchi, mesh); r.end()
    if next(stream, None) is not None:
        raise ValueError("the file has more records than the USPP 7.3.6 layout")
    return UsppData(version, date, title, z, zv, exfact, valence, mesh, etot, keyps, ifpcor, nang,
                    lloc, eloc, ifqopt, nqf, qtryc, rinner, irel, rc, kkbeta, l, eee, beta, ddd0,
                    ddd, qqq, qfunc, qfcoef, iptype, npf, ptryc, rcloc, vloc0, rpcor, rspsco,
                    *arrays, wf)


def job_summary(stages: List[Tuple[str, float, OutReport]], uspp: Optional[UsppData] = None) -> str:
    """The Results tab: per stage time, warnings, eigenvalue table and ghost check."""
    lines = []
    for label, seconds, report in stages:
        lines.append(f"{label}: finished in {seconds:.2f} s")
        lines += [f"  WARNING: {w}" for w in report.warnings]
        if report.comparison:
            lines.append("    nlm     all-elec       pseudo         diff   (Ry)")
            for row in report.comparison:
                flag = "  <-- more than 1 mRy" if abs(row.diff) > EIGENVALUE_FLAG_RY else ""
                lines.append(f"    {row.nlm:3d} {row.all_electron:12.6f} {row.pseudo:12.6f} {row.diff:12.6f}{flag}")
        if report.comparison and report.bessel:
            ghosts = ghost_candidates(report)
            lines.append("  ghost check: " + (f"possible ghost state for l = {', '.join(map(str, ghosts))} "
                                              "(see 'eigenvalues of US hamiltonian' in the .out file)"
                                              if ghosts else "no extra bound state below the expected eigenvalues"))
        if report.total_energy is not None:
            lines.append(f"  total energy = {report.total_energy:.7f} Ry")
        lines.append("")
    if uspp is not None:
        lines.append(f"Pseudopotential file: {uspp.title}, zv = {uspp.zv:g}, nbeta = {uspp.nbeta}, "
                     f"mesh = {uspp.mesh}, kkbeta = {uspp.kkbeta}, "
                     f"rc = {', '.join(f'{v:g}' for v in uspp.rc)} bohr, "
                     f"core correction = {'yes' if uspp.ifpcor else 'no'}")
    return "\n".join(lines)
```

- [ ] **Step 4: Run them to see them pass**

Run: `python -m pytest tests/test_runatomOutput.py -v`
Expected: all passed.

- [ ] **Step 5: Commit (if authorised)**

```bash
git add runatomOutput.py tests/test_runatomOutput.py
git commit -m "feat(uspp): parse runatom .out reports and read .uspp files"
```

---

### Task 8: Running runatom (`runatomRunner.py`)

**Files:**
- Create: `runatomRunner.py`
- Create: `tests/fake_runatom.py`
- Test: `tests/test_runatomRunner.py`

**Interfaces:**
- Consumes: `runatomInput.CaseInputs`, `file_stems`, `MAX_FORTRAN_NAME`; `runatomOutput.parse_out`, `OutReport`.
- Produces:
  - `find_runatom(explicit=None) -> Optional[str]`, `library_dir() -> Path`, `runatom_command(exe) -> List[str]`;
  - `Stage(kind, label, stem, ae_file, uspp_file)` with `.arguments()`, `.products()`, `.required()`;
  - `plan_stages(case) -> List[Stage]`, `StageResult(stage, seconds, report, returncode)`;
  - `RunatomError(kind, stage_label, message)` (kinds: `input_format`, `program_error`, `not_converged`, `negative_density`, `runtime_error`, `missing_output`, `timeout`), `RunatomCancelled`;
  - `run_stage(...)`, `run_job(command, stages, workdir, progress=None, cancel=None, timeout=DEFAULT_STAGE_TIMEOUT) -> List[StageResult]`;
  - `progress_text(info) -> str`, where `info` = `{"stage_index", "stage_count", "stage", "message", "iteration", "delta"}`.

- [ ] **Step 1: Write the fake program**

Create `tests/fake_runatom.py`:

```python
"""Stand-in for runatom in tests (never the real program).

Called like runatom: input out ae atwf logd uspp. FAKE_RUNATOM_MODE selects:
ok, error, stop0, runtime, negative, hang. Every call is logged to fake_calls.jsonl.
"""

import json
import os
import sys
import time

SAMPLE_OUT = """ comparison of all-electron and pseudo eigenvalues (Ry)

   nlm    all-elec      pseudo        diff

   200   -1.000000   -1.000001   -0.000001
   210   -0.500000   -0.500000    0.000000

 -----------------------------------------------
   eigenvalues of US hamiltonian
   l= 0      -1.000000   0.300000
   l= 1      -0.500000   0.700000
"""

mode = os.environ.get("FAKE_RUNATOM_MODE", "ok")
inp, out, ae, atwf, logd, uspp = sys.argv[1:7]
with open("fake_calls.jsonl", "a") as log:
    log.write(json.dumps({"argv": sys.argv[1:], "cwd": os.getcwd()}) + "\n")
kind = "ae" if uspp == "dummy" else ("test" if "_test_" in inp else "pseudo")
err = sys.stderr
err.write("beginning execution pseudopotential program version  7.3.6\n")
err.flush()
if mode == "runtime":
    err.write(f"At line 165 of file rwsubs.f (unit = 10, file = '{inp}')\n"
              "Fortran runtime error: Bad value during floating point read\n")
    sys.exit(2)
if mode == "negative":
    err.write(" *** stopping in vhxc: negative densities ***\n")
    sys.exit(1)
if mode == "stop0":
    sys.exit(0)
for cycle in range(1, 4):
    err.write(f"completed self-consistent cycle{cycle:4d} delta ={10.0 ** -cycle:15.7E}\n".replace("E", "D"))
    err.flush()
if mode == "hang":
    for name in (out, ae if kind == "ae" else uspp):
        open(name, "w").close()
    while True:
        time.sleep(0.05)
if mode == "error":
    with open(out, "w") as f:
        f.write(" ***error in subroutine readin\n nang =  5 but idim5 =  4\n")
    sys.exit(1)
err.write({"ae": " all electron calculation completed\n",
           "pseudo": " writing the pseudopotential to file\n",
           "test": " pseudopotential test converged\n"}[kind])
with open(out, "w") as f:
    f.write(" value of mesh generated in rinit is   851\n" if kind == "ae" else SAMPLE_OUT)
if kind == "ae":
    open(ae, "wb").write(b"ae")
if kind == "pseudo":
    open(uspp, "wb").write(b"uspp")
```

- [ ] **Step 2: Write the failing tests**

Create `tests/test_runatomRunner.py`:

```python
"""runatomRunner: stages, child processes, progress, cancel, timeout and errors (fake program)."""

import json
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

import runatomInput as ri
import runatomRunner as rr
from runatom_fixtures import LIBRARY_DIR

SRC = Path(__file__).resolve().parents[1]
FAKE = [sys.executable, "-I", str(Path(__file__).resolve().parent / "fake_runatom.py")]


def _ti():
    return ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-test")


def _job(tmp_path, mode, monkeypatch, **kwargs):
    monkeypatch.setenv("FAKE_RUNATOM_MODE", mode)
    case = _ti()
    ri.write_case(case, tmp_path)
    return rr.run_job(FAKE, rr.plan_stages(case), tmp_path, **kwargs)


def test_module_does_not_import_toga():
    code = "import sys, runatomRunner; assert 'toga' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], cwd=SRC, check=True)


def test_stages_reproduce_the_library_makefile():
    stages = rr.plan_stages(_ti())
    assert [s.arguments() for s in stages] == [
        ["ti_ae_s2d1.adat", "ti_ae_s2d1.out", "ti_ae_s2d1.ae", "ti_ae_s2d1.atwf", "ti_ae_s2d1.logd", "dummy"],
        ["ti_ps.adat", "ti_ps.out", "ti_ae_s2d1.ae", "ti_ps.atwf", "ti_ps.logd", "ti_ps.uspp"],
        ["ti_ae_s2d2.adat", "ti_ae_s2d2.out", "ti_ae_s2d2.ae", "ti_ae_s2d2.atwf", "ti_ae_s2d2.logd", "dummy"],
        ["ti_test_s2d2.adat", "ti_test_s2d2.out", "ti_ae_s2d2.ae", "ti_test_s2d2.atwf", "ti_test_s2d2.logd", "ti_ps.uspp"],
    ]


def test_successful_job_runs_every_stage_in_the_job_folder(tmp_path, monkeypatch):
    results = _job(tmp_path, "ok", monkeypatch)
    assert [r.stage.kind for r in results] == ["ae", "pseudo", "ae", "test"]
    calls = [json.loads(l) for l in (tmp_path / "fake_calls.jsonl").read_text().splitlines()]
    assert all(Path(c["cwd"]) == tmp_path for c in calls)
    assert results[1].report.comparison[0].nlm == 200
    assert (tmp_path / "ti_ps.uspp").exists()


def test_progress_reports_stages_and_scf_cycles(tmp_path, monkeypatch):
    seen = []
    _job(tmp_path, "ok", monkeypatch, progress=seen.append)
    assert {i["stage_index"] for i in seen} == {0, 1, 2, 3}
    cycles = [i for i in seen if i["iteration"] is not None]
    assert cycles[0]["iteration"] == 1 and cycles[0]["delta"] == pytest.approx(0.1)
    assert rr.progress_text(cycles[0]).startswith("Stage 1/4: All-electron (s2d1) - SCF cycle 1")


def test_a_failing_progress_callback_never_stops_the_run(tmp_path, monkeypatch):
    def broken(info):
        raise RuntimeError("UI gone")
    assert len(_job(tmp_path, "ok", monkeypatch, progress=broken)) == 4


@pytest.mark.parametrize("mode, kind, text", [
    ("error", "program_error", "nang =  5 but idim5 =  4"),
    ("runtime", "input_format", "ti_ae_s2d1.adat"),
    ("negative", "negative_density", "rinner"),
    ("stop0", "missing_output", "stopped before the end"),
])
def test_failures_are_classified(tmp_path, monkeypatch, mode, kind, text):
    with pytest.raises(rr.RunatomError) as info:
        _job(tmp_path, mode, monkeypatch)
    assert info.value.kind == kind and text in str(info.value)
    assert info.value.stage_label == "All-electron (s2d1)"


def test_stale_outputs_are_removed_before_a_stage(tmp_path, monkeypatch):
    (tmp_path / "ti_ae_s2d1.ae").write_bytes(b"old")
    with pytest.raises(rr.RunatomError):
        _job(tmp_path, "stop0", monkeypatch)
    assert not (tmp_path / "ti_ae_s2d1.ae").exists()


def test_cancel_kills_the_program_and_removes_its_outputs(tmp_path, monkeypatch):
    cancel = threading.Event()
    threading.Timer(0.5, cancel.set).start()
    started = time.monotonic()
    with pytest.raises(rr.RunatomCancelled):
        _job(tmp_path, "hang", monkeypatch, cancel=cancel)
    assert time.monotonic() - started < 5
    assert not (tmp_path / "ti_ae_s2d1.out").exists() and not (tmp_path / "ti_ae_s2d1.ae").exists()


def test_timeout_stops_a_hung_program(tmp_path, monkeypatch):
    with pytest.raises(rr.RunatomError) as info:
        _job(tmp_path, "hang", monkeypatch, timeout=0.5)
    assert info.value.kind == "timeout"


def test_names_longer_than_fortran_accepts_are_refused(tmp_path):
    stage = rr.Stage("ae", "x", "a" * 40, "a.ae", "dummy")
    with pytest.raises(ValueError, match="40"):
        rr.run_stage(FAKE, stage, tmp_path)


class TestFindRunatom:
    def test_explicit_path_first(self, tmp_path, monkeypatch):
        exe = tmp_path / "my.exe"
        exe.write_bytes(b"")
        monkeypatch.setattr(rr, "VENDOR_EXE", tmp_path / "missing.exe")
        assert rr.find_runatom(str(exe)) == str(exe)

    def test_environment_variable(self, tmp_path, monkeypatch):
        exe = tmp_path / "env.exe"
        exe.write_bytes(b"")
        monkeypatch.setenv(rr.ENV_VAR, str(exe))
        assert rr.find_runatom() == str(exe)

    def test_frozen_build_looks_in_the_bundle(self, tmp_path, monkeypatch):
        exe = tmp_path / rr.FROZEN_SUBDIR / "runatom.exe"
        exe.parent.mkdir()
        exe.write_bytes(b"")
        monkeypatch.delenv(rr.ENV_VAR, raising=False)
        monkeypatch.setattr(sys, "frozen", True, raising=False)
        monkeypatch.setattr(sys, "_MEIPASS", str(tmp_path), raising=False)
        assert rr.find_runatom() == str(exe)
        assert rr.library_dir() == tmp_path / rr.FROZEN_SUBDIR / "Work"

    def test_nothing_found(self, tmp_path, monkeypatch):
        monkeypatch.delenv(rr.ENV_VAR, raising=False)
        monkeypatch.setattr(rr, "VENDOR_EXE", tmp_path / "missing.exe")
        monkeypatch.setattr(rr.shutil, "which", lambda name: None)
        assert rr.find_runatom() is None
```

- [ ] **Step 3: Run them to see them fail**

Run: `python -m pytest tests/test_runatomRunner.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'runatomRunner'`.

- [ ] **Step 4: Implement**

Create `runatomRunner.py`:

```python
"""Running the Vanderbilt USPP program runatom as a child process.

No Toga here (a subprocess test pins it). runatom stops on errors with
`call exit(1)` (134 places) or `stop` (exit status 0) and takes its six file
names as character*40 command-line arguments. It therefore runs as a separate
process, with its working directory set to the job folder and short relative
names. A stage succeeds when the exit status is 0, its completion marker is on
stderr and its output files exist.
"""

import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional, Sequence

from runatomInput import MAX_FORTRAN_NAME, CaseInputs, file_stems
from runatomOutput import OutReport, parse_out

DEFAULT_STAGE_TIMEOUT = 120.0     # s; the slowest measured stage takes 0.7 s
POLL_INTERVAL = 0.05              # s between cancel/timeout checks
ENV_VAR = "GQTEA_RUNATOM"
FROZEN_SUBDIR = "uspp"
VENDOR_DIR = Path(__file__).resolve().parent / "third_party" / "uspp-736"
VENDOR_EXE = VENDOR_DIR / "bin" / "runatom.exe"
COMPLETION_MARKERS = {"ae": "all electron calculation completed",
                      "pseudo": "writing the pseudopotential to file",
                      "test": "pseudopotential test converged"}
_SCF = re.compile(r"completed self-consistent cycle\s*(\d+)\s+delta\s*=\s*([-+0-9.EeDd]+)")
_NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)


class RunatomError(RuntimeError):
    def __init__(self, kind: str, stage_label: str, message: str):
        super().__init__(f"{stage_label}: {message}")
        self.kind, self.stage_label, self.message = kind, stage_label, message


class RunatomCancelled(RuntimeError):
    pass


def _bundle_dir() -> Optional[Path]:
    if getattr(sys, "frozen", False):
        return Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent)) / FROZEN_SUBDIR
    return None


def _find_program(explicit: Optional[str], env_var: str, exe_name: str, vendor_path: Path,
                  path_names: Sequence[str]) -> Optional[str]:
    """Explicit path, the environment variable, the frozen bundle, third_party, PATH; None if absent."""
    candidates = []
    if explicit and explicit.strip():
        candidates.append(explicit.strip())
    if os.environ.get(env_var):
        candidates.append(os.environ[env_var])
    bundle = _bundle_dir()
    if bundle is not None:
        candidates.append(str(bundle / exe_name))
    candidates.append(str(vendor_path))
    for name in path_names:
        found = shutil.which(name)
        if found:
            candidates.append(found)
    return next((c for c in candidates if os.path.isfile(c)), None)


def find_runatom(explicit: Optional[str] = None) -> Optional[str]:
    return _find_program(explicit, ENV_VAR, "runatom.exe", VENDOR_EXE, ("runatom.exe", "runatom.x"))


def library_dir() -> Path:
    bundle = _bundle_dir()
    return bundle / "Work" if bundle is not None else VENDOR_DIR / "Work"


def runatom_command(exe: str) -> List[str]:
    return [exe]


@dataclass(frozen=True)
class Stage:
    kind: str           # "ae", "pseudo" or "test"
    label: str
    stem: str           # <stem>.adat, .out, .atwf, .logd
    ae_file: str        # .ae written (ae) or read (pseudo, test)
    uspp_file: str      # "dummy" for an all-electron stage

    def arguments(self) -> List[str]:
        return [f"{self.stem}.adat", f"{self.stem}.out", self.ae_file,
                f"{self.stem}.atwf", f"{self.stem}.logd", self.uspp_file]

    def required(self) -> List[str]:
        files = [f"{self.stem}.out"]
        files += [self.ae_file] if self.kind == "ae" else []
        files += [self.uspp_file] if self.kind == "pseudo" else []
        return files

    def products(self) -> List[str]:
        return self.required() + [f"{self.stem}.atwf", f"{self.stem}.logd"]


@dataclass
class StageResult:
    stage: Stage
    seconds: float
    report: OutReport
    returncode: int


def plan_stages(case: CaseInputs) -> List[Stage]:
    """The Work/000-Ex Makefile's chain: AE, generation, then per test AE + test."""
    stems = file_stems(case)
    ae_ref, uspp = f"{stems['ae']}.ae", f"{stems['ps']}.uspp"
    stages = [Stage("ae", f"All-electron ({case.generation_name})", stems["ae"], ae_ref, "dummy"),
              Stage("pseudo", "Pseudopotential generation", stems["ps"], ae_ref, uspp)]
    for test, (ae_stem, test_stem) in zip(case.tests, stems["tests"]):
        stages.append(Stage("ae", f"All-electron ({test.name})", ae_stem, f"{ae_stem}.ae", "dummy"))
        stages.append(Stage("test", f"Transferability test ({test.name})", test_stem, f"{ae_stem}.ae", uspp))
    return stages


def progress_text(info: dict) -> str:
    text = f"Stage {info['stage_index'] + 1}/{info['stage_count']}: {info['stage']}"
    if info.get("iteration") is not None:
        return f"{text} - SCF cycle {info['iteration']}, delta = {info['delta']:.2e}"
    return f"{text} - {info['message']}"


def _tail(stderr_text: str) -> str:
    lines = [l.strip() for l in stderr_text.splitlines() if l.strip() and "self-consistent cycle" not in l]
    return " | ".join(lines[-3:])


def classify(stage: Stage, returncode: int, stderr_text: str, out_text: str, workdir: Path) -> Optional[RunatomError]:
    def error(kind, message):
        return RunatomError(kind, stage.label, message)

    if "Fortran runtime error" in stderr_text:
        detail = " ".join(l.strip() for l in stderr_text.splitlines()
                          if l.startswith("At line") or "Fortran runtime error" in l)
        if "unit = 10" in stderr_text:
            return error("input_format", f"runatom could not read {stage.stem}.adat ({detail}); "
                                         "the file does not follow runatom's fixed columns.")
        return error("runtime_error", detail)
    if "negative densities" in stderr_text:
        return error("negative_density", "the pseudo charge density became negative while descreening; "
                                         "reduce rinner or add Q-function constraints (Doc/TUTORIAL).")
    lines = out_text.splitlines()
    for i, line in enumerate(lines):
        if "***error" in line or "*** error" in line:
            detail = " ".join(l.strip() for l in lines[i:i + 3] if l.strip())
            return error("not_converged" if "could not be converged" in detail else "program_error", detail)
    if returncode != 0:
        return error("runtime_error", f"runatom stopped with exit status {returncode}: {_tail(stderr_text)}")
    if COMPLETION_MARKERS[stage.kind] not in stderr_text:
        tail = _tail(stderr_text)
        return error("missing_output", "runatom stopped before the end of this stage" + (f": {tail}" if tail else "."))
    missing = [name for name in stage.required() if not (workdir / name).is_file()]
    if missing:
        return error("missing_output", "expected output not written: " + ", ".join(missing))
    return None


def _pump(stream, lines: "queue.Queue"):
    for raw in iter(stream.readline, b""):
        lines.put(raw.decode("latin-1"))
    lines.put(None)


def _stop(process: subprocess.Popen):
    if process.poll() is None:
        process.kill()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass


def _remove(workdir: Path, names: Sequence[str]):
    for name in names:
        try:
            (workdir / name).unlink()
        except FileNotFoundError:
            pass


def run_stage(command: Sequence[str], stage: Stage, workdir, progress: Optional[Callable[[dict], None]] = None,
              cancel: Optional[threading.Event] = None, timeout: float = DEFAULT_STAGE_TIMEOUT,
              stage_index: int = 0, stage_count: int = 1) -> StageResult:
    workdir = Path(workdir)
    arguments = stage.arguments()
    too_long = [a for a in arguments if len(a) > MAX_FORTRAN_NAME]
    if too_long:
        raise ValueError(f"file names longer than {MAX_FORTRAN_NAME} characters: {', '.join(too_long)}")

    def report(message: str, iteration: Optional[int] = None, delta: Optional[float] = None):
        if progress is None:
            return
        try:
            progress({"stage_index": stage_index, "stage_count": stage_count, "stage": stage.label,
                      "message": message, "iteration": iteration, "delta": delta})
        except Exception:      # reporting must never abort a run
            pass

    _remove(workdir, stage.products())
    report("starting")
    started = time.monotonic()
    process = subprocess.Popen(list(command) + arguments, cwd=str(workdir), stdin=subprocess.DEVNULL,
                               stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, creationflags=_NO_WINDOW)
    lines: "queue.Queue" = queue.Queue()
    pump = threading.Thread(target=_pump, args=(process.stderr, lines), daemon=True)
    pump.start()
    collected, finished_reading = [], False
    try:
        while not (finished_reading and process.poll() is not None):
            if cancel is not None and cancel.is_set():
                raise RunatomCancelled(f"{stage.label}: cancelled")
            if time.monotonic() - started > timeout:
                raise RunatomError("timeout", stage.label, f"no result after {timeout:.0f} s; the run was stopped.")
            try:
                line = lines.get(timeout=POLL_INTERVAL)
            except queue.Empty:
                continue
            if line is None:
                finished_reading = True
                continue
            collected.append(line)
            match = _SCF.search(line)
            if match:
                report(f"SCF cycle {match.group(1)}", int(match.group(1)),
                       float(match.group(2).upper().replace("D", "E")))
            elif line.strip():
                report(line.strip())
    except (RunatomCancelled, RunatomError):
        _stop(process)
        _remove(workdir, stage.products())
        raise
    finally:
        _stop(process)
        pump.join(timeout=5)       # the killed process closes the pipe; never close it under the reader
        process.stderr.close()
    seconds = time.monotonic() - started
    out_path = workdir / f"{stage.stem}.out"
    out_text = out_path.read_text(encoding="latin-1") if out_path.is_file() else ""
    problem = classify(stage, process.returncode, "".join(collected), out_text, workdir)
    if problem is not None:
        _remove(workdir, [n for n in stage.products() if n != f"{stage.stem}.out"])  # keep the report
        raise problem
    return StageResult(stage, seconds, parse_out(out_text), process.returncode)


def run_job(command: Sequence[str], stages: Sequence[Stage], workdir,
            progress: Optional[Callable[[dict], None]] = None, cancel: Optional[threading.Event] = None,
            timeout: float = DEFAULT_STAGE_TIMEOUT) -> List[StageResult]:
    """All stages in order; the first failure raises RunatomError (or RunatomCancelled)."""
    return [run_stage(command, stage, workdir, progress, cancel, timeout, index, len(stages))
            for index, stage in enumerate(stages)]
```

- [ ] **Step 5: Run them to see them pass**

Run: `python -m pytest tests/test_runatomRunner.py -v`
Expected: all passed.

- [ ] **Step 6: Commit (if authorised)**

```bash
git add runatomRunner.py tests/fake_runatom.py tests/test_runatomRunner.py
git commit -m "feat(uspp): run runatom stages as child processes with progress, cancel and timeout"
```

---

### Task 9: Formatted pseudopotential for CPMD (`reform`)

CPMD reads a Vanderbilt potential given as `*<file> FORMATTED` with `readvan`, in the record order of `reform.f`'s formatted output (spec §6.1, checked in CPMD's `ratom_utils.mod.F90` and `readvan_utils.mod.F90`). The file is written by `reform.exe`, the unmodified upstream converter.
- Python only **verifies** it: `read_formatted_uspp` reads it in `readvan`'s statement order.
- The runner normalises line endings to LF. MinGW gfortran writes CRLF, CPMD runs on Linux, and the Linux reference uses LF.

**Files:**
- Modify: `runatomRunner.py` (append), `runatomOutput.py` (append)
- Create: `tests/fake_reform.py`
- Test: `tests/test_runatomFormatted.py`

**Interfaces:**
- Consumes: `runatomRunner._find_program`, `_bundle_dir`, `_remove`, `_tail`, `_NO_WINDOW`, `RunatomError`, `DEFAULT_STAGE_TIMEOUT`, `VENDOR_DIR`; `runatomInput.file_stems`; `runatomOutput.UsppData`; Task 3 reference files `<atom>_ps.uspp`, `<atom>_ps_formatted.uspp`.
- Produces:
  - `runatomRunner`: `REFORM_ENV_VAR`, `VENDOR_REFORM`, `MAX_REFORM_NAME = 36`, `find_reform(explicit=None) -> Optional[str]`, `reform_command(exe) -> List[str]`, `formatted_name(case) -> str` (`<symbol>_ps_formatted.uspp`), `cpmd_atom_line(name) -> str`;
  - `write_formatted(command, workdir, uspp_name, formatted, timeout=DEFAULT_STAGE_TIMEOUT) -> Path` (raises `RunatomError` with kind `formatted` or `timeout`, label `"Formatted file for CPMD"`);
  - `runatomOutput`: `read_formatted_uspp(path) -> UsppData`.

- [ ] **Step 1: Write the fake converter**

Create `tests/fake_reform.py`:

```python
"""Stand-in for USPP's reform (itype 1) in tests; never the real program.

Reads itype, input and output names from stdin like reform.f. FAKE_REFORM_MODE:
ok (writes the output with CRLF line endings, as MinGW gfortran does), fail
(missing input), silent (exit 0 without the completion message).
"""

import os
import sys

mode = os.environ.get("FAKE_REFORM_MODE", "ok")
itype, source, target = (sys.stdin.readline().strip() for _ in range(3))
print(" selected itype =", itype)
if os.path.exists(target):        # reform.f opens its output with status='new'
    sys.stderr.write(f"Fortran runtime error: Cannot open file '{target}': File exists\n")
    sys.exit(2)
if mode == "fail" or not os.path.exists(source):
    sys.stderr.write(f"Fortran runtime error: Cannot open file '{source}': No such file or directory\n")
    sys.exit(2)
with open(target, "w", newline="\r\n") as f:
    f.write("    7    3    6\nformatted\n")
if mode != "silent":
    print("1pseudopotential has been written out")
```

- [ ] **Step 2: Write the failing tests**

Create `tests/test_runatomFormatted.py`:

```python
"""The CPMD formatted pseudopotential: reform conversion and the readvan-order reader."""

import sys
from pathlib import Path

import numpy as np
import pytest

import runatomInput as ri
import runatomOutput as ro
import runatomRunner as rr
from runatom_fixtures import LIBRARY_DIR, REFERENCE_CASES, case_info, reference_bytes

FAKE = [sys.executable, "-I", str(Path(__file__).resolve().parent / "fake_reform.py")]
REL = 1e-11          # 1p4e19.11 keeps 12 significant digits
ABS = {"z": 5e-10, "zv": 5e-10, "exfact": 5e-10, "ptryc": 5e-10,    # written f15.9 by reform
       "eloc": 5e-6, "qtryc": 5e-6}                                  # written f9.5 by reform
SCALARS = ("version", "date", "title", "z", "zv", "exfact", "valence", "mesh", "etot", "keyps", "ifpcor",
           "nang", "lloc", "eloc", "ifqopt", "nqf", "qtryc", "irel", "kkbeta", "npf", "ptryc", "rcloc", "rpcor")
ARRAYS = ("rinner", "rc", "l", "eee", "beta", "ddd0", "ddd", "qqq", "qfunc", "qfcoef", "iptype",
          "vloc0", "rspsco", "vloc", "rsatom", "r", "rab", "wf")


def _files(name, tmp_path):
    atom = case_info(name)["atom"]
    binary, formatted = tmp_path / f"{atom}_ps.uspp", tmp_path / f"{atom}_ps_formatted.uspp"
    binary.write_bytes(reference_bytes(name, binary.name))
    return binary, formatted


def _assert_same(a: ro.UsppData, b: ro.UsppData):
    for attr in SCALARS:
        x, y = getattr(a, attr), getattr(b, attr)
        if isinstance(x, float):
            assert x == pytest.approx(y, rel=REL, abs=ABS.get(attr, 0.0)), attr
        elif attr == "valence":          # (i5,2f15.9)
            assert [(n, pytest.approx(w, abs=5e-10), pytest.approx(e, abs=5e-10)) for n, w, e in x] == y, attr
        else:
            assert x == y, attr
    for attr in ARRAYS:
        x, y = getattr(a, attr), getattr(b, attr)
        if x is None or y is None:
            assert x is None and y is None, attr
        else:
            assert np.allclose(x, y, rtol=REL, atol=0), attr


class TestReader:
    def test_fortran_reals_as_list_directed_input_reads_them(self):
        assert ro._fortran_float("1.23456789012-100") == pytest.approx(1.23456789012e-100)
        assert ro._fortran_float("-4.50000000000+200") == -4.5e200
        assert ro._fortran_float("1.0D-03") == 1e-3 and ro._fortran_float("-7") == -7.0

    @pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
    def test_reference_formatted_file_holds_the_binary_one(self, name, tmp_path):
        binary, formatted = _files(name, tmp_path)
        formatted.write_bytes(reference_bytes(name, formatted.name))
        _assert_same(ro.read_formatted_uspp(formatted), ro.read_uspp(binary))

    def test_truncated_file_is_refused(self, tmp_path):
        binary, formatted = _files("000-Ex-ca-sp-test", tmp_path)
        text = reference_bytes("000-Ex-ca-sp-test", formatted.name).decode("latin-1")
        formatted.write_text("\n".join(text.splitlines()[:40]) + "\n", encoding="latin-1")
        with pytest.raises(ValueError, match="ends early"):
            ro.read_formatted_uspp(formatted)


class TestConversion:
    def test_names_and_cpmd_line(self):
        case = ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-test")
        assert rr.formatted_name(case) == "ti_ps_formatted.uspp"
        assert rr.cpmd_atom_line("ti_ps_formatted.uspp") == "*ti_ps_formatted.uspp FORMATTED"

    def test_conversion_replaces_an_old_file_and_writes_lf(self, tmp_path, monkeypatch):
        monkeypatch.setenv("FAKE_REFORM_MODE", "ok")
        (tmp_path / "ti_ps.uspp").write_bytes(b"binary")
        (tmp_path / "ti_ps_formatted.uspp").write_text("old")
        path = rr.write_formatted(FAKE, tmp_path, "ti_ps.uspp", "ti_ps_formatted.uspp")
        assert path.read_bytes() == b"    7    3    6\nformatted\n"

    @pytest.mark.parametrize("mode, text", [("fail", "No such file"), ("silent", "could not convert")])
    def test_failures_are_reported(self, tmp_path, monkeypatch, mode, text):
        monkeypatch.setenv("FAKE_REFORM_MODE", mode)
        (tmp_path / "ti_ps.uspp").write_bytes(b"binary")
        with pytest.raises(rr.RunatomError) as info:
            rr.write_formatted(FAKE, tmp_path, "ti_ps.uspp", "ti_ps_formatted.uspp")
        assert info.value.kind == "formatted" and text in str(info.value)
        assert not (tmp_path / "ti_ps_formatted.uspp").exists()

    def test_names_longer_than_reform_reads_are_refused(self, tmp_path):
        with pytest.raises(ValueError, match="36"):
            rr.write_formatted(FAKE, tmp_path, "a" * 37, "b.uspp")

    def test_frozen_build_finds_the_bundled_reform(self, tmp_path, monkeypatch):
        exe = tmp_path / rr.FROZEN_SUBDIR / "reform.exe"
        exe.parent.mkdir()
        exe.write_bytes(b"")
        monkeypatch.delenv(rr.REFORM_ENV_VAR, raising=False)
        monkeypatch.setattr(sys, "frozen", True, raising=False)
        monkeypatch.setattr(sys, "_MEIPASS", str(tmp_path), raising=False)
        assert rr.find_reform() == str(exe)


REFORM = rr.find_reform()


@pytest.mark.skipif(REFORM is None, reason="reform.exe is not available")
@pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
def test_reform_exe_reproduces_the_linux_reference_byte_for_byte(name, tmp_path):
    binary, formatted = _files(name, tmp_path)
    rr.write_formatted(rr.reform_command(REFORM), tmp_path, binary.name, formatted.name)
    assert formatted.read_bytes() == reference_bytes(name, formatted.name)
```

- [ ] **Step 3: Run them to see them fail**

Run: `python -m pytest tests/test_runatomFormatted.py -v`
Expected: FAIL with `AttributeError: module 'runatomOutput' has no attribute '_fortran_float'`.

- [ ] **Step 4: Implement the converter**

Append to `runatomRunner.py`:

```python
REFORM_ENV_VAR = "GQTEA_REFORM"
VENDOR_REFORM = VENDOR_DIR / "bin" / "reform.exe"
MAX_REFORM_NAME = 36          # reform.f reads its file names with format a36
REFORM_DONE = "pseudopotential has been written out"
FORMATTED_LABEL = "Formatted file for CPMD"


def find_reform(explicit: Optional[str] = None) -> Optional[str]:
    return _find_program(explicit, REFORM_ENV_VAR, "reform.exe", VENDOR_REFORM, ("reform.exe", "reform.x"))


def reform_command(exe: str) -> List[str]:
    return [exe]


def formatted_name(case: CaseInputs) -> str:
    return f"{file_stems(case)['ps']}_formatted.uspp"


def cpmd_atom_line(name: str) -> str:
    """The &ATOMS species line that makes CPMD read the file as a formatted Vanderbilt potential."""
    return f"*{name} FORMATTED"


def write_formatted(command: Sequence[str], workdir, uspp_name: str, formatted: str,
                    timeout: float = DEFAULT_STAGE_TIMEOUT) -> Path:
    """The CPMD formatted copy of a binary .uspp, written by USPP's reform (itype 1), with LF endings."""
    workdir = Path(workdir)
    for name in (uspp_name, formatted):
        if len(name) > MAX_REFORM_NAME:
            raise ValueError(f"reform reads file names of at most {MAX_REFORM_NAME} characters: {name}")
    target = workdir / formatted
    _remove(workdir, [formatted])                  # reform opens its output with status='new'
    try:
        done = subprocess.run(list(command), cwd=str(workdir), input=f"1\n{uspp_name}\n{formatted}\n".encode("ascii"),
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout,
                              creationflags=_NO_WINDOW)
    except subprocess.TimeoutExpired:
        _remove(workdir, [formatted])
        raise RunatomError("timeout", FORMATTED_LABEL, f"no result after {timeout:.0f} s; the conversion was stopped.")
    stdout, stderr = done.stdout.decode("latin-1"), done.stderr.decode("latin-1")
    if done.returncode != 0 or REFORM_DONE not in stdout or not target.is_file():
        _remove(workdir, [formatted])
        detail = _tail(stderr)
        raise RunatomError("formatted", FORMATTED_LABEL,
                           f"reform could not convert {uspp_name}" + (f": {detail}" if detail else "."))
    data = target.read_bytes()
    if b"\r\n" in data:            # MinGW gfortran writes CRLF; CPMD runs on Linux
        target.write_bytes(data.replace(b"\r\n", b"\n"))
    return target
```

- [ ] **Step 5: Implement the reader**

Append to `runatomOutput.py`:

```python
_EXPONENT_WITHOUT_E = re.compile(r"(?<=\d)([+-]\d)")


def _fortran_float(token: str) -> float:
    """A real as Fortran list-directed input reads it ('1.0D-03', and '1.23456789012-100',
    which is how an E format writes a 3-digit exponent)."""
    text = token.upper().replace("D", "E")
    if "E" not in text:
        text = _EXPONENT_WITHOUT_E.sub(r"E\1", text, count=1)
    return float(text)


class _ListReader:
    """List-directed READs: each starts on a new record and may continue over several."""

    def __init__(self, text: str):
        self.lines, self.index = text.splitlines(), 0

    def line(self) -> str:
        if self.index >= len(self.lines):
            raise ValueError("the formatted file ends early")
        self.index += 1
        return self.lines[self.index - 1]

    def read(self, count: int) -> List[str]:
        tokens: List[str] = []
        while len(tokens) < count:
            tokens += self.line().replace(",", " ").split()
        return tokens[:count]

    def ints(self, count: int) -> List[int]:
        return [int(t) for t in self.read(count)]

    def reals(self, count: int) -> np.ndarray:
        return np.array([_fortran_float(t) for t in self.read(count)])


def read_formatted_uspp(path) -> UsppData:
    """A formatted .uspp (reform itype 1) read in CPMD readvan's statement order.

    This is the CPMD-compatibility check: CPMD reads every record list-directed
    except the title, '(A20,3F15.9/2I5,1PE19.11)'.
    """
    r = _ListReader(Path(path).read_text(encoding="latin-1"))
    head = r.ints(6)
    version, date = tuple(head[:3]), tuple(head[3:])
    if version < (7, 3, 2):
        raise ValueError(f"file version {version} predates the 7.3.2 layout this reader follows")
    first, second = r.line(), r.line()
    title = first[:20].strip()
    z, zv, exfact = (_fortran_float(first[20 + 15 * k:35 + 15 * k]) for k in range(3))
    nval, mesh, etot = int(second[0:5]), int(second[5:10]), _fortran_float(second[10:29])
    t = r.read(3 * nval)
    valence = [(int(t[3 * i]), _fortran_float(t[3 * i + 1]), _fortran_float(t[3 * i + 2])) for i in range(nval)]
    t = r.read(3)
    keyps, ifpcor = int(t[0]), int(t[1])
    t = r.read(6)
    nang, lloc, eloc, ifqopt, nqf, qtryc = int(t[0]), int(t[1]), _fortran_float(t[2]), int(t[3]), int(t[4]), _fortran_float(t[5])
    n_l = 2 * nang - 1
    rinner = r.reals(n_l)
    (irel,) = r.ints(1)
    rc = r.reals(nang)
    nbeta, kkbeta = r.ints(2)
    l, eee = np.zeros(nbeta, int), np.zeros(nbeta)
    beta = np.zeros((nbeta, kkbeta))
    ddd0, ddd, qqq = (np.zeros((nbeta, nbeta)) for _ in range(3))
    qfunc = np.zeros((nbeta, nbeta, kkbeta))
    qfcoef = np.zeros((nbeta, nbeta, n_l, nqf))
    for j in range(nbeta):
        t = r.read(2 + kkbeta)
        l[j], eee[j] = int(t[0]), _fortran_float(t[1])
        beta[j] = [_fortran_float(x) for x in t[2:]]
        for k in range(j, nbeta):
            v = r.reals(3 + kkbeta + nqf * n_l)
            ddd0[j, k] = ddd0[k, j] = v[0]
            ddd[j, k] = ddd[k, j] = v[1]
            qqq[j, k] = qqq[k, j] = v[2]
            qfunc[j, k] = qfunc[k, j] = v[3:3 + kkbeta]
            qfcoef[j, k] = qfcoef[k, j] = v[3 + kkbeta:].reshape(n_l, nqf)
    t = r.read(nbeta + 2)
    iptype = np.array([int(x) for x in t[:nbeta]])
    npf, ptryc = int(t[nbeta]), _fortran_float(t[nbeta + 1])
    v = r.reals(1 + mesh)
    rcloc, vloc0 = float(v[0]), v[1:]
    rpcor = rspsco = None
    if ifpcor > 0:
        rpcor = float(r.reals(1)[0])
        rspsco = r.reals(mesh)
    vloc, rsatom, radii, rab = (r.reals(mesh) for _ in range(4))
    (nchi,) = r.ints(1)
    wf = r.reals(mesh * nchi).reshape(nchi, mesh)
    return UsppData(version, date, title, z, zv, exfact, valence, mesh, etot, keyps, ifpcor, nang,
                    lloc, eloc, ifqopt, nqf, qtryc, rinner, irel, rc, kkbeta, l, eee, beta, ddd0,
                    ddd, qqq, qfunc, qfcoef, iptype, npf, ptryc, rcloc, vloc0, rpcor, rspsco,
                    vloc, rsatom, radii, rab, wf)
```

- [ ] **Step 6: Run them to see them pass**

Run: `python -m pytest tests/test_runatomFormatted.py -v`
Expected: all passed. The byte-for-byte reform test runs once `reform.exe` exists (Task 2); before that it is skipped.

If the byte comparison fails only in the last digit of some numbers, compare the two gfortran versions' rounding of `1PE19.11` (16.2 vs 11.5) before anything else. Both should round correctly. Record the finding in the spec and never weaken the test without the user's agreement.

- [ ] **Step 7: Commit (if authorised)**

```bash
git add runatomRunner.py runatomOutput.py tests/fake_reform.py tests/test_runatomFormatted.py
git commit -m "feat(uspp): CPMD formatted pseudopotential via USPP reform, verified in readvan order"
```

---

### Task 10: Validate `runatom.exe` against the reference, and benchmark

**Files:**
- Create: `tests/test_runatomExe.py`
- Create: `benchmarks/bench_runatom.py`

**Interfaces:**
- Consumes: everything above; `runatom.exe` from Task 2 (tests skip when it is absent).

- [ ] **Step 1: Write the tests**

Create `tests/test_runatomExe.py`:

```python
"""runatom.exe against the Linux reference build (spec §7). Skipped without runatom.exe.

Rule: a disagreement is investigated in the build (flags, compiler), never fixed
by loosening a tolerance.
"""

import shutil
import time

import numpy as np
import pytest

import runatomInput as ri
import runatomOutput as ro
import runatomRunner as rr
from runatom_fixtures import LIBRARY_DIR, REFERENCE_CASES, case_info, library_case_dir, reference_bytes, reference_names, reference_text

EXE = rr.find_runatom()
pytestmark = pytest.mark.skipif(EXE is None, reason="runatom.exe is not available")

EIGEN_TOL = 2e-6          # Ry: one unit in the last printed place
BESSEL_TOL = 1e-5         # Ry
MESH_RTOL = 1e-12
POT_RTOL, POT_ATOL = 1e-7, 1e-10
PROJ_RTOL = 1e-6
UNSTABLE = {"082-Pb-ca-d-vgrp": {0}}   # Source/BENIGN_BUG: the s channel of Pb


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    """Every reference case run once with runatom.exe on the original library inputs."""
    done = {}
    for name in REFERENCE_CASES:
        work = tmp_path_factory.mktemp(name)
        shutil.copytree(library_case_dir(name), work, dirs_exist_ok=True)
        case = ri.read_case_folder(work)
        started = time.monotonic()
        rr.run_job(rr.runatom_command(EXE), rr.plan_stages(case), work)
        done[name] = (work, case, time.monotonic() - started)
    return done


@pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
def test_eigenvalue_tables_match_the_reference(runs, name):
    work = runs[name][0]
    for filename in (n for n in reference_names(name) if n.endswith(".out")):
        ours = ro.parse_out((work / filename).read_text(encoding="latin-1"))
        ref = ro.parse_out(reference_text(name, filename))
        assert [r.nlm for r in ours.comparison] == [r.nlm for r in ref.comparison], filename
        for a, b in zip(ours.comparison, ref.comparison):
            assert abs(a.all_electron - b.all_electron) <= EIGEN_TOL and abs(a.pseudo - b.pseudo) <= EIGEN_TOL, filename
        assert sorted(ours.bessel) == sorted(ref.bessel), filename
        for l in ref.bessel:
            assert np.allclose(ours.bessel[l], ref.bessel[l], rtol=0, atol=BESSEL_TOL), (filename, l)
        if ref.total_energy is not None:
            assert abs(ours.total_energy - ref.total_energy) <= EIGEN_TOL, filename


def _signs(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.where(np.einsum("ij,ij->i", a, b) < 0, -1.0, 1.0)


def _block_eigenvalues(matrix: np.ndarray, l: np.ndarray, channel: int) -> np.ndarray:
    index = np.flatnonzero(l == channel)
    return np.linalg.eigvalsh(matrix[np.ix_(index, index)])


@pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
def test_pseudopotential_file_matches_the_reference(runs, name, tmp_path):
    atom = case_info(name)["atom"]
    ours = ro.read_uspp(runs[name][0] / f"{atom}_ps.uspp")
    (tmp_path / "ref.uspp").write_bytes(reference_bytes(name, f"{atom}_ps.uspp"))
    ref = ro.read_uspp(tmp_path / "ref.uspp")
    for attr in ("version", "title", "z", "zv", "exfact", "mesh", "keyps", "ifpcor", "nang", "lloc",
                 "ifqopt", "nqf", "irel", "kkbeta", "npf"):
        assert getattr(ours, attr) == getattr(ref, attr), attr
    assert np.array_equal(ours.l, ref.l) and np.array_equal(ours.iptype, ref.iptype)
    assert np.allclose(ours.r, ref.r, rtol=MESH_RTOL, atol=0) and np.allclose(ours.rab, ref.rab, rtol=MESH_RTOL, atol=0)
    for attr in ("vloc", "vloc0", "rsatom"):
        assert np.allclose(getattr(ours, attr), getattr(ref, attr), rtol=POT_RTOL, atol=POT_ATOL), attr
    unstable = UNSTABLE.get(name, set())
    stable = np.array([l not in unstable for l in ref.l])
    s = _signs(ours.beta, ref.beta)
    scale = np.abs(ref.beta).max()
    assert np.allclose((s[:, None] * ours.beta)[stable], ref.beta[stable], rtol=PROJ_RTOL, atol=PROJ_RTOL * scale)
    pair = np.outer(stable, stable)
    for attr in ("ddd", "qqq"):
        mine = np.outer(s, s) * getattr(ours, attr)
        assert np.allclose(mine[pair], getattr(ref, attr)[pair], rtol=PROJ_RTOL, atol=1e-9), attr
    for channel in unstable:      # representation-independent check (BENIGN_BUG)
        for attr in ("ddd", "qqq"):
            assert np.allclose(_block_eigenvalues(getattr(ours, attr), ours.l, channel),
                               _block_eigenvalues(getattr(ref, attr), ref.l, channel), rtol=1e-6, atol=1e-9)


@pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
def test_reformatted_inputs_give_an_identical_pseudopotential(runs, name, tmp_path):
    work, case, _ = runs[name]
    ri.write_case(ri.replace(case, tests=[]), tmp_path)        # our writer, same values
    rr.run_job(rr.runatom_command(EXE), rr.plan_stages(ri.replace(case, tests=[])), tmp_path)
    stem = ri.file_stems(case)["ps"]
    original = ro.uspp_records((work / f"{stem}.uspp").read_bytes())
    rewritten = ro.uspp_records((tmp_path / f"{stem}.uspp").read_bytes())
    assert original[1:] == rewritten[1:]          # record 1 holds the run date


def test_sample_output_of_2006_is_reproduced(runs):
    work = runs["000-Ex-ca-sp-test"][0]
    sample = LIBRARY_DIR / "000-Ex" / "Sample_output"
    for filename in ("ti_ps.out", "ti_test_s2d2.out"):
        ours = ro.parse_out((work / filename).read_text(encoding="latin-1")).comparison
        old = ro.parse_out((sample / filename).read_text(encoding="latin-1")).comparison
        assert [r.nlm for r in ours] == [r.nlm for r in old]
        for a, b in zip(ours, old):
            assert abs(a.all_electron - b.all_electron) <= EIGEN_TOL and abs(a.pseudo - b.pseudo) <= EIGEN_TOL


def test_unreadable_library_input_is_reported_as_an_input_problem(tmp_path):
    shutil.copytree(LIBRARY_DIR / "023-V" / "023-V-gpw-nsp-campos", tmp_path, dirs_exist_ok=True)
    stage = rr.Stage("pseudo", "V generation", "V_ps", "V_ae_ref.ae", "V_ps.uspp")
    rr.run_stage(rr.runatom_command(EXE), rr.Stage("ae", "V AE", "V_ae_ref", "V_ae_ref.ae", "dummy"), tmp_path)
    with pytest.raises(rr.RunatomError) as info:
        rr.run_stage(rr.runatom_command(EXE), stage, tmp_path)
    assert info.value.kind == "input_format"


def test_every_reference_case_is_fast(runs):
    assert max(seconds for _, _, seconds in runs.values()) < 10.0
```

- [ ] **Step 2: Run them**

Run: `python -m pytest tests/test_runatomExe.py -v`
Expected: all passed. If any fail:
1. First compare `BUILDINFO.txt` with `tests/data/runatom/MANIFEST.txt`.
2. Then rebuild once with `-O0 -ffp-contract=off` added to both builds.
3. Record what changed in `BUILDINFO.txt` and in the spec's §9 risk table.
4. Never widen a tolerance to make a test pass.

The one admissible exception is the `Source/BENIGN_BUG` case: projectors of one channel differ by more than a sign, while the eigenvalues of that channel's D and Q blocks agree within the tolerances. That channel may then be added to `UNSTABLE`, with the measured numbers written into the spec.

- [ ] **Step 3: Write the benchmark**

Create `benchmarks/bench_runatom.py`:

```python
"""Time runatom.exe on every library case and measure event-loop lag while it runs.

Usage (from venv/src):  python benchmarks/bench_runatom.py [--lag]
"""

import asyncio
import shutil
import statistics
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import runatomInput as ri        # noqa: E402
import runatomRunner as rr       # noqa: E402


def library_cases():
    for makefile in sorted(rr.library_dir().rglob("Makefile")):
        folder = makefile.parent
        if folder.parent == rr.library_dir() or (folder / "V_ps.adat").exists():
            continue
        yield folder


def time_library(exe: str):
    times = []
    for folder in library_cases():
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copytree(folder, tmp, dirs_exist_ok=True)
            case = ri.read_case_folder(tmp)
            started = time.perf_counter()
            results = rr.run_job(rr.runatom_command(exe), rr.plan_stages(case), tmp)
            total = time.perf_counter() - started
        times.append(total)
        stages = ", ".join(f"{r.stage.kind} {r.seconds:.2f}" for r in results)
        print(f"{folder.name:34s} {total:6.2f} s  ({stages})")
    print(f"\n{len(times)} cases: median {statistics.median(times):.2f} s, max {max(times):.2f} s")


async def measure_lag(exe: str, folder: Path):
    loop = asyncio.get_running_loop()
    lags, done = [], asyncio.Event()

    async def ticker():
        while not done.is_set():
            before = loop.time()
            await asyncio.sleep(0.01)
            lags.append(loop.time() - before - 0.01)

    with tempfile.TemporaryDirectory() as tmp:
        shutil.copytree(folder, tmp, dirs_exist_ok=True)
        case = ri.read_case_folder(tmp)
        tick = asyncio.create_task(ticker())
        progress = lambda info: loop.call_soon_threadsafe(lambda: None)
        await asyncio.to_thread(rr.run_job, rr.runatom_command(exe), rr.plan_stages(case), tmp, progress)
        done.set()
        await tick
    lags.sort()
    print(f"event-loop lag during a run: p95 {1000 * lags[int(0.95 * len(lags))]:.1f} ms, max {1000 * lags[-1]:.1f} ms")


if __name__ == "__main__":
    exe = rr.find_runatom()
    if exe is None:
        sys.exit("runatom.exe not found")
    time_library(exe)
    if "--lag" in sys.argv:
        asyncio.run(measure_lag(exe, rr.library_dir() / "056-Ba" / "056-Ba-gpw-nsp-campos"))
```

- [ ] **Step 4: Run the benchmark and record the numbers**

Run: `python benchmarks/bench_runatom.py --lag`
Acceptance:
- every case < 3 s;
- event-loop lag p95 ≤ 50 ms.

Paste the summary lines into the spec's §7 under "Measured".

- [ ] **Step 5: Commit (if authorised)**

```bash
git add tests/test_runatomExe.py benchmarks/bench_runatom.py docs/superpowers/specs/2026-10-28-runatom-integration-design.md
git commit -m "test(uspp): validate runatom.exe against the Linux reference; library benchmark"
```

---

### Task 11: The window (`runatom_input_gen.py`) and help text

**Files:**
- Rewrite: `runatom_input_gen.py` (module name and `RunatomGenerator` class kept)
- Modify: `help.py` (add `help_runatom` to `HelpGqteaWin`)
- Modify: `toolLoader.py:50` (label only)
- Test: `tests/test_runatom_input_gen.py`

**Interfaces:**
- Consumes: `runatomInput` (`FORM_FIELDS`, `case_from_fields`, `case_to_fields`, `ae_fields`, `ps_fields`, `check_case`, `write_case`, `makefile_text`, `parse_adat`, `read_case_folder`, `file_stems`, `AEInput`, `PSInput`), `runatomRunner` (`find_runatom`, `library_dir`, `runatom_command`, `plan_stages`, `run_job`, `progress_text`, `RunatomError`, `RunatomCancelled`, and from Task 9 `find_reform`, `reform_command`, `formatted_name`, `write_formatted`, `cpmd_atom_line`) and `runatomOutput` (`read_uspp`, `job_summary`).
- Produces: `RunatomGenerator(*args)`, with methods the tests use: `_collect_case()`, `_fill(fields)`, `_progress_reporter(loop)`, `on_run(widget)`, `on_cancel(widget)`, `_on_close(window)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_runatom_input_gen.py`:

```python
"""The runatom window, headless: stub widgets, fake program, no Toga window."""

import asyncio
import sys
import threading
from pathlib import Path

import pytest

import runatomInput as ri
import runatomRunner as rr
import runatom_input_gen as gen
from runatom_fixtures import LIBRARY_DIR

FAKE = [sys.executable, "-I", str(Path(__file__).resolve().parent / "fake_runatom.py")]
FAKE_REFORM = [sys.executable, "-I", str(Path(__file__).resolve().parent / "fake_reform.py")]


class Value:
    def __init__(self, value=""):
        self.value = value
        self.enabled = True


class Label:
    text = ""


class Spinner:
    running = False

    def start(self):
        self.running = True

    def stop(self):
        self.running = False


class Tabs:
    current_tab = 0


class Window:
    def __init__(self):
        self.dialogs = []

    async def dialog(self, d):
        self.dialogs.append(d)


def _window(tmp_path):
    w = gen.RunatomGenerator.__new__(gen.RunatomGenerator)
    w._closed, w._busy, w._cancel, w._last_out = False, False, None, None
    case = ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-test")
    w.inputs = {k: Value(v) for k, v in ri.case_to_fields(case).items()}
    w.folder_input, w.exe_input = Value(str(tmp_path)), Value("")
    w.makefile_switch = Value(False)
    w.results_text = Value("")
    w.status_label, w.activity, w.tabs, w.main_window = Label(), Spinner(), Tabs(), Window()
    for name in ("check_button", "save_button", "run_button", "cancel_button"):
        setattr(w, name, Value())
    w.errors = []                     # a Toga dialog keeps no readable message: record _error instead

    async def record(title, message):
        w.errors.append((title, message))

    w._error = record
    return w, case


def test_every_form_field_has_a_widget():
    assert set(ri.FORM_FIELDS) == set(gen.NAME_FIELDS) | {k for _, keys in gen.TAB_FIELDS for k in keys}


def test_collect_case_reads_the_widgets(tmp_path):
    w, case = _window(tmp_path)
    assert w._collect_case() == ri.case_from_fields(ri.case_to_fields(case))


def test_run_writes_inputs_runs_every_stage_and_shows_results(tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "ok")
    monkeypatch.setenv("FAKE_REFORM_MODE", "ok")
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: "runatom.exe")
    monkeypatch.setattr(gen.rr, "runatom_command", lambda exe: FAKE)
    monkeypatch.setattr(gen.rr, "find_reform", lambda explicit=None: "reform.exe")
    monkeypatch.setattr(gen.rr, "reform_command", lambda exe: FAKE_REFORM)
    w, _ = _window(tmp_path)
    asyncio.run(w.on_run(None))
    assert (tmp_path / "ti_ps.adat").exists() and (tmp_path / "ti_ps.uspp").exists()
    assert (tmp_path / "ti_ps_formatted.uspp").exists()
    assert "*ti_ps_formatted.uspp FORMATTED" in w.results_text.value
    assert "Pseudopotential generation: finished" in w.results_text.value
    assert w.tabs.current_tab == len(gen.TAB_FIELDS)
    assert w.errors == [] and w.status_label.text.startswith("Finished: 4 stages")
    assert not w._busy and w.run_button.enabled and not w.activity.running


def test_missing_reform_leaves_the_binary_result_and_says_so(tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "ok")
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: "runatom.exe")
    monkeypatch.setattr(gen.rr, "runatom_command", lambda exe: FAKE)
    monkeypatch.setattr(gen.rr, "find_reform", lambda explicit=None: None)
    w, _ = _window(tmp_path)
    asyncio.run(w.on_run(None))
    assert "reform.exe was not found" in w.results_text.value and w.errors == []
    assert not (tmp_path / "ti_ps_formatted.uspp").exists()


def test_program_error_is_shown_in_a_dialog(tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_RUNATOM_MODE", "error")
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: "runatom.exe")
    monkeypatch.setattr(gen.rr, "runatom_command", lambda exe: FAKE)
    w, _ = _window(tmp_path)
    asyncio.run(w.on_run(None))
    ((title, message),) = w.errors
    assert title == "runatom stopped" and "idim5" in message and w.run_button.enabled


def test_invalid_inputs_never_reach_the_program(tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr(gen.rr, "run_job", lambda *a, **k: called.append(a))
    w, _ = _window(tmp_path)
    w.inputs["ae_atom"].value = "22 0.0 0"            # occupations no longer add up
    asyncio.run(w.on_run(None))
    assert called == [] and "add up to" in w.errors[0][1]


def test_missing_program_is_explained(tmp_path, monkeypatch):
    monkeypatch.setattr(gen.rr, "find_runatom", lambda explicit=None: None)
    w, _ = _window(tmp_path)
    asyncio.run(w.on_run(None))
    assert "runatom.exe was not found" in w.errors[0][1]


def test_progress_is_thinned_and_marshalled(tmp_path):
    w, _ = _window(tmp_path)
    posted = []

    class Loop:
        def call_soon_threadsafe(self, fn, *args):
            posted.append(args)

    report = w._progress_reporter(Loop())
    info = {"stage_index": 0, "stage_count": 4, "stage": "All-electron (s2d1)", "iteration": None, "delta": None}
    report(dict(info, message="starting"))
    for k in range(50):
        report(dict(info, message="SCF", iteration=k, delta=0.1))
    assert 1 <= len(posted) <= 3 and posted[0][0]["message"] == "starting"


def test_closing_the_window_cancels_a_run(tmp_path):
    w, _ = _window(tmp_path)
    w._cancel = threading.Event()
    assert w._on_close(None) is True
    assert w._closed and w._cancel.is_set()
```

- [ ] **Step 2: Run them to see them fail**

Run: `python -m pytest tests/test_runatom_input_gen.py -v`
Expected: FAIL with `AttributeError: module 'runatom_input_gen' has no attribute 'NAME_FIELDS'`.

- [ ] **Step 3: Rewrite the window**

Replace the whole of `runatom_input_gen.py`:

```python
"""Inputs > Vanderbilt USPP Generator: prepare, check and run runatom (USPP 7.3.6).

The window only moves text between widgets and a dict of field values.
runatomInput turns that dict into checked, column-exact input files,
runatomRunner runs runatom on them in a worker thread, and runatomOutput reads
the results.
"""

import asyncio
import os
import threading
import time
from pathlib import Path

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import runatomInput as ri
import runatomOutput as ro
import runatomRunner as rr
from help import HelpGqteaWin

STATUS_INTERVAL = 0.1          # s between two status-line updates during a run
LABEL_WIDTH = 240
HINT_COLOR = "#666666"
MULTILINE_HEIGHT = 110
DEFAULT_CASE = ("000-Ex", "000-Ex-ca-sp-test")
NAME_FIELDS = ("symbol", "generation", "title")
TAB_FIELDS = (
    ("All-electron", ("ae_print", "ae_logd", "ae_conv", "ae_atom", "ae_mesh", "ae_shells")),
    ("Pseudopotential", ("ps_print", "ps_logd", "ps_conv", "ps_partition", "ps_bessel", "ps_keyps",
                         "ps_rinner", "ps_rcloc", "ps_rc", "ps_betas", "ps_npf", "ps_local",
                         "ps_qfunc", "ps_fixes")),
    ("Tests", ("tests",)),
)


class RunatomGenerator:
    def __init__(self, *args, **kwargs):
        # gqteaWinToga passes the Button that was pressed; ignore it.
        self._closed, self._busy, self._cancel, self._last_out = False, False, None, None
        self.inputs = {}
        self.main_window = toga.Window(title="Vanderbilt USPP Generator (runatom)", size=(920, 700))
        self.main_window.on_close = self._on_close
        self.layout_main_window()
        self.main_window.show()
        self._load_case(rr.library_dir().joinpath(*DEFAULT_CASE), quiet=True)

    # --- layout -----------------------------------------------------------
    def _row(self, label, *widgets):
        return toga.Box(style=Pack(direction=ROW, margin_bottom=2), children=[
            toga.Label(label, style=Pack(width=LABEL_WIDTH, margin=(4, 6, 0, 0))), *widgets])

    def _field_row(self, key):
        label, _, multiline, hint = ri.FORM_FIELDS[key]
        style = Pack(flex=1, height=MULTILINE_HEIGHT) if multiline else Pack(flex=1)
        self.inputs[key] = toga.MultilineTextInput(style=style) if multiline else toga.TextInput(style=style)
        box = toga.Box(style=Pack(direction=COLUMN, margin_bottom=4), children=[self._row(label, self.inputs[key])])
        if hint:
            box.add(toga.Label(hint, style=Pack(color=HINT_COLOR, font_size=9, margin_left=LABEL_WIDTH)))
        return box

    def layout_main_window(self):
        self.status_label = toga.Label("", style=Pack(flex=1, margin=(6, 6)))
        self.activity = toga.ActivityIndicator(style=Pack(margin=(6, 6)))
        self.folder_input = toga.TextInput(placeholder="folder for the input and output files", style=Pack(flex=1))
        self.exe_input = toga.TextInput(placeholder="blank = the runatom.exe shipped with gqteaWinToga",
                                        style=Pack(flex=1), on_change=self._update_exe_status)
        self.exe_status = toga.Label("", style=Pack(color=HINT_COLOR, font_size=9, margin_left=LABEL_WIDTH))
        self.makefile_switch = toga.Switch("Also write a Makefile (to run the case by hand under Linux/WSL)", value=False)
        header = toga.Box(style=Pack(direction=COLUMN, margin=8), children=[
            self._row("Job folder", self.folder_input, toga.Button("Browse...", on_press=self.on_browse_folder)),
            self._row("runatom program", self.exe_input), self.exe_status,
            toga.Box(style=Pack(direction=ROW, margin_top=4), children=[
                toga.Button("Import library case...", on_press=self.on_import_case),
                toga.Button("Import .adat file...", on_press=self.on_import_adat),
                self.makefile_switch])])
        names = toga.Box(style=Pack(direction=COLUMN, margin=(0, 8)),
                         children=[self._field_row(k) for k in NAME_FIELDS])
        pages = []
        for title, keys in TAB_FIELDS:
            page = toga.Box(style=Pack(direction=COLUMN, margin=8), children=[self._field_row(k) for k in keys])
            pages.append((title, toga.ScrollContainer(content=page, horizontal=False)))
        self.results_text = toga.MultilineTextInput(readonly=True, style=Pack(flex=1, font_family="monospace"))
        pages.append(("Results", toga.Box(style=Pack(direction=COLUMN, margin=8), children=[self.results_text])))
        help_text = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_runatom, style=Pack(flex=1))
        pages.append(("Help", toga.Box(style=Pack(direction=COLUMN, margin=8), children=[help_text])))
        self.tabs = toga.OptionContainer(content=pages, style=Pack(flex=1))
        self.check_button = toga.Button("Check", on_press=self.on_check)
        self.save_button = toga.Button("Save inputs", on_press=self.on_save)
        self.run_button = toga.Button("Run", on_press=self.on_run)
        self.cancel_button = toga.Button("Cancel", on_press=self.on_cancel, enabled=False)
        bar = toga.Box(style=Pack(direction=ROW, margin=8), children=[
            self.check_button, self.save_button, self.run_button, self.cancel_button,
            toga.Button("Open folder", on_press=self.on_open_folder),
            toga.Button("Open .out", on_press=self.on_open_out), self.activity, self.status_label])
        self.main_window.content = toga.Box(style=Pack(direction=COLUMN),
                                            children=[header, names, self.tabs, bar])
        self._update_exe_status(None)

    # --- widgets <-> fields ------------------------------------------------
    def _collect_case(self) -> ri.CaseInputs:
        return ri.case_from_fields({key: widget.value or "" for key, widget in self.inputs.items()})

    def _fill(self, fields):
        for key, text in fields.items():
            widget = self.inputs.get(key)
            if widget is not None and widget.value != text:
                widget.value = text

    def _load_case(self, folder: Path, quiet: bool = False):
        try:
            case = ri.read_case_folder(folder)
        except (OSError, ValueError, KeyError) as exc:
            if not quiet:
                raise ValueError(f"{folder.name}: {exc}") from None
            return
        self._fill(ri.case_to_fields(case))
        self._show_status(f"Loaded {folder.name}. Test settings now follow the Pseudopotential tab.")

    def _update_exe_status(self, widget):
        exe = rr.find_runatom(self.exe_input.value)
        text = f"using {exe}" if exe else "runatom.exe not found: enter its path"
        if self.exe_status.text != text:
            self.exe_status.text = text

    # --- status, busy, progress --------------------------------------------
    def _show_status(self, message: str):
        if not self._closed and self.status_label.text != message:
            self.status_label.text = message

    def _set_busy(self, busy: bool, message: str):
        self._busy = busy
        if self._closed:
            return
        (self.activity.start if busy else self.activity.stop)()
        for button in (self.check_button, self.save_button, self.run_button):
            button.enabled = not busy
        self.cancel_button.enabled = busy
        self._show_status(message)

    def _apply_progress(self, info: dict):
        self._show_status(rr.progress_text(info))

    def _progress_reporter(self, loop):
        """progress(info) for the worker thread: thinned, and handed to the window's loop."""
        last = [0.0]

        def report(info: dict):
            now = time.monotonic()
            if info["message"] != "starting" and now - last[0] < STATUS_INTERVAL:
                return
            last[0] = now
            if not self._closed:
                try:
                    loop.call_soon_threadsafe(self._apply_progress, info)
                except RuntimeError:     # the event loop has closed (application exit)
                    pass

        return report

    async def _error(self, title: str, message: str):
        if not self._closed:
            await self.main_window.dialog(toga.ErrorDialog(title, message))

    # --- actions ----------------------------------------------------------
    def _job_folder(self) -> Path:
        text = (self.folder_input.value or "").strip()
        if not text:
            raise ValueError("Choose the job folder first (Browse...).")
        return Path(text)

    def _prepared_case(self) -> ri.CaseInputs:
        case = self._collect_case()
        problems = ri.check_case(case)
        if problems:
            raise ValueError("\n".join(problems))
        return case

    def _write_inputs(self, case, folder: Path):
        written = ri.write_case(case, folder)
        if self.makefile_switch.value:
            (folder / "Makefile").write_text(ri.makefile_text(case), encoding="ascii", newline="\n")
            written.append(folder / "Makefile")
        return written

    async def on_check(self, widget):
        try:
            self._prepared_case()
        except ValueError as exc:
            await self._error("Check the inputs", str(exc))
            return
        self._show_status("The inputs pass every check runatom makes before it starts.")

    async def on_save(self, widget):
        try:
            case, folder = self._prepared_case(), self._job_folder()
            written = self._write_inputs(case, folder)
        except (ValueError, OSError) as exc:
            await self._error("Inputs not saved", str(exc))
            return
        self._show_status(f"Saved {len(written)} files in {folder}")

    async def on_run(self, widget):
        if self._busy:
            return
        try:
            case, folder = self._prepared_case(), self._job_folder()
            exe = rr.find_runatom(self.exe_input.value)
            if exe is None:
                raise ValueError("runatom.exe was not found. Enter its path, or reinstall gqteaWinToga, which ships it.")
            self._write_inputs(case, folder)
        except (ValueError, OSError) as exc:
            await self._error("Cannot run runatom", str(exc))
            return
        stages = rr.plan_stages(case)
        self._cancel = threading.Event()
        self._set_busy(True, f"Running {len(stages)} stages...")
        loop = asyncio.get_running_loop()
        try:
            results = await asyncio.to_thread(rr.run_job, rr.runatom_command(exe), stages, folder,
                                              self._progress_reporter(loop), self._cancel)
        except rr.RunatomCancelled:
            self._set_busy(False, "Run cancelled; the stage that was running left no output.")
            return
        except rr.RunatomError as exc:
            self._set_busy(False, f"Stopped at: {exc.stage_label}")
            await self._error("runatom stopped", f"{exc}\n\nThe .out file of that stage is in {folder}.")
            return
        except OSError as exc:
            self._set_busy(False, "runatom could not be started")
            await self._error("runatom could not be started", str(exc))
            return
        stem = ri.file_stems(case)["ps"]
        self._last_out = folder / f"{stem}.out"
        try:
            uspp, note = ro.read_uspp(folder / f"{stem}.uspp"), ""
        except (OSError, ValueError) as exc:
            uspp, note = None, f"\n\nThe pseudopotential file could not be read back: {exc}"
        note += await self._write_cpmd_file(case, folder, f"{stem}.uspp")
        text = ro.job_summary([(r.stage.label, r.seconds, r.report) for r in results], uspp) + note
        if not self._closed:
            self.results_text.value = text
            self.tabs.current_tab = len(TAB_FIELDS)
        self._set_busy(False, f"Finished: {len(results)} stages in {sum(r.seconds for r in results):.1f} s")

    async def _write_cpmd_file(self, case, folder: Path, uspp_name: str) -> str:
        """Write the formatted copy for CPMD with reform.exe; return the Results text about it."""
        typed = (self.exe_input.value or "").strip()
        reform = rr.find_reform(str(Path(typed).with_name("reform.exe")) if typed else None)
        if reform is None:
            return ("\n\nreform.exe was not found, so no formatted file for CPMD was written "
                    "(the binary .uspp is complete).")
        name = rr.formatted_name(case)
        try:
            await asyncio.to_thread(rr.write_formatted, rr.reform_command(reform), folder, uspp_name, name)
        except (rr.RunatomError, ValueError, OSError) as exc:
            return f"\n\nThe formatted file for CPMD was not written: {exc}"
        return (f"\n\nCPMD: {name} is the formatted pseudopotential. Copy it to CPMD's PP_LIBRARY_PATH "
                f"(or the run directory) and use the &ATOMS species line\n  {rr.cpmd_atom_line(name)}")

    def on_cancel(self, widget):
        if self._cancel is not None:
            self._cancel.set()
            self._show_status("Cancelling...")

    async def on_browse_folder(self, widget):
        path = await self.main_window.dialog(toga.SelectFolderDialog("Job folder"))
        if path:
            self.folder_input.value = str(path)

    async def on_import_case(self, widget):
        path = await self.main_window.dialog(toga.SelectFolderDialog(
            "Library case folder (a folder with a Makefile)", initial_directory=rr.library_dir()))
        if path:
            try:
                self._load_case(Path(path))
            except ValueError as exc:
                await self._error("Case not imported", str(exc))

    async def on_import_adat(self, widget):
        path = await self.main_window.dialog(toga.OpenFileDialog("Input file", file_types=["adat"]))
        if not path:
            return
        try:
            data = ri.parse_adat(Path(path).read_text(encoding="latin-1"))
        except (OSError, ValueError) as exc:
            await self._error("File not imported", str(exc))
            return
        if isinstance(data, ri.AEInput):
            self._fill(ri.ae_fields(data))
        elif isinstance(data, ri.PSInput):
            self._fill(ri.ps_fields(data))
        else:
            await self._error("File not imported", "Test inputs are written from the Pseudopotential tab; "
                                                    "import the test's all-electron file and list it on the Tests tab.")
            return
        self._show_status(f"Imported {Path(path).name}")

    def on_open_folder(self, widget):
        text = (self.folder_input.value or "").strip()
        if text and hasattr(os, "startfile") and Path(text).is_dir():
            os.startfile(text)

    def on_open_out(self, widget):
        if self._last_out is not None and self._last_out.is_file() and hasattr(os, "startfile"):
            os.startfile(str(self._last_out))

    def _on_close(self, window, **kwargs):
        self._closed = True        # a run still going finishes without touching the window
        if self._cancel is not None:
            self._cancel.set()
        return True

    def closeTopLevel(self, widget):
        self._on_close(self.main_window)
        self.main_window.close()
```

- [ ] **Step 4: Add the help text**

Add to `HelpGqteaWin` in `help.py`, beside `help_fukui_input`:

```python
    help_runatom = """VANDERBILT USPP GENERATOR (runatom, USPP 7.3.6)

Generates an ultrasoft pseudopotential with D. Vanderbilt's program runatom,
which is shipped with gqteaWinToga (Fortran, GPL v2; source in
third_party/uspp-736). No separate installation is needed.

WORKFLOW
1. Choose a job folder (Browse...). Every input and output file is written there.
2. Start from a library case (Import library case...: the Work/ folder of USPP
   holds tested inputs for 45 elements) or from the Ti example loaded at start.
3. Edit the tabs. Type values separated by spaces; the tool writes them in the
   fixed columns runatom reads, so spacing does not matter.
4. Check runs every test runatom makes before it starts; Run writes the inputs and
   runs: all-electron reference -> pseudopotential generation -> for each test,
   all-electron test configuration -> transferability test.
5. Results shows each stage's time, warnings, the all-electron vs pseudo
   eigenvalue table (differences above 1 mRy are flagged), and a ghost-state
   check. Open .out shows the full report with its plots of log derivatives,
   wavefunctions and Q functions.

FILES (in the job folder)
<symbol>_ae_<config>.adat/.out/.ae   all-electron run
<symbol>_ps.adat/.out/.uspp          generation; .uspp is the pseudopotential
                                     (binary, gfortran little-endian records)
<symbol>_ps_formatted.uspp           the same pseudopotential, formatted, for
                                     CPMD (written by USPP's reform). In the
                                     CPMD input, &ATOMS line:
                                       *<symbol>_ps_formatted.uspp FORMATTED
<symbol>_test_<name>.adat/.out       transferability test

FIELDS (see also Doc/INPUT_AE, Doc/INPUT_GEN and Doc/TUTORIAL of USPP)
All-electron:  ifprt ifplw ilogd; rlogd emin emax nnt; thresh tol damp maxit;
               z xion exfact (0 CA, -1 Wigner, -2 HL, -3 GL, 1 BLYP, 2 B88,
               3 BP86, 4 PW91, 5 PBE); rmax aasf bbsf irel; one shell per line:
               nnlz wwnl ee (e.g. 320 1.0 -0.2 = 3d, 1 electron, guess -0.2 Ry).
Pseudopotential: ncores nvales nang; Bessel-basis ghost check; keyps (3)
               ifpcor rinner rpcor; rcloc; rc per l; one projector per line:
               l keyee eeread iptype; npf ptryc (when some iptype is 2);
               lloc keyee eloc iploctype; ifqopt nqf qtryc; Q constraints.
Tests:         one per line: name xion nnlz=occupation ...  Shells not listed
               keep the generation occupation.

RULES WORTH KNOWING
- GGA functionals (exfact 1-5) need ifqopt = 3; runatom refuses -1, 0, 1 and 2,
  although Doc/TUTORIAL suggests starting with -1.
- The radial mesh may have at most 1000 points (lower rmax or bbsf).
- File names are limited to 40 characters, so keep the symbol and configuration
  names short.
- A large q-matrix eigenvalue gives a "mild" (> 5) or "severe" (> 20) warning;
  runatom stops above 50.
- If a run reports negative densities, reduce rinner or add Q constraints.
- If a test does not converge, lower damp to 0.1-0.2.

Not supported: reference wavefunctions from another configuration (keyee < 0),
Dirac generation (irel 1), UPF export.
"""
```

- [ ] **Step 5: Update the button label**

In `toolLoader.py:50` change the label only:

```python
        ToolSpec("Vanderbilt USPP Generator (runatom)", "runatom_input_gen", "RunatomGenerator"),
```

- [ ] **Step 6: Run the tests**

Run: `python -m pytest tests/test_runatom_input_gen.py tests/test_toolLoader.py -v`
Expected: all passed.

- [ ] **Step 7: Check it in the real window**

Run `python gqteaWinToga.py`, open **Inputs > Vanderbilt USPP Generator (runatom)**, choose an empty folder and press **Run**. Check that:
- the Ti example runs 4 stages;
- the status line advances through them;
- the Results tab shows 4 sections with no flagged difference, and ends with the CPMD line `*ti_ps_formatted.uspp FORMATTED`;
- `ti_ps_formatted.uspp` exists in the folder with LF line endings;
- **Cancel** during a run of the Ba case stops it and leaves no `.uspp`;
- the window does not freeze.

Record what you checked in the commit message.

- [ ] **Step 8: Commit (if authorised)**

```bash
git add runatom_input_gen.py help.py toolLoader.py tests/test_runatom_input_gen.py
git commit -m "feat(uspp): run runatom from the Vanderbilt USPP Generator window"
```

---

### Task 12: Packaging and documentation

**Files:**
- Modify: `readme.md` (packaging section), `USER_MANUAL.md` (replace the *Vanderbilt runatom.x Input Builder* section and its TOC entry), `CLAUDE.md` (tool note + numba/optional lines untouched), `installation_requirements.txt` (no change; record that no Python package is added)
- Test: `tests/test_runatomVendor.py` (append one test)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_runatomVendor.py`:

```python
def test_readme_tells_packagers_to_bundle_runatom_and_its_source():
    readme = (VENDOR.parents[1] / "readme.md").read_text(encoding="utf-8")
    assert "third_party/uspp-736/bin/runatom.exe" in readme and "uspp/" in readme
    assert "third_party/uspp-736/bin/reform.exe" in readme and "FORMATTED" in readme
```

- [ ] **Step 2: Run it to see it fail**

Run: `python -m pytest tests/test_runatomVendor.py::test_readme_tells_packagers_to_bundle_runatom_and_its_source -v`
Expected: FAIL (assertion).

- [ ] **Step 3: Add the packaging section to `readme.md`**

Put it after "Tool modules in the packaged build":

```markdown
### runatom (Vanderbilt USPP generator) in the packaged build

The Vanderbilt USPP Generator runs two programs that **must ship with
gqteaWinToga**:
- `runatom.exe`, a statically linked build of USPP 7.3.6;
- `reform.exe`, USPP's converter, which writes the formatted file CPMD reads
  with the `&ATOMS` line `*<file> FORMATTED`.

Both are committed in `third_party/uspp-736/bin/`. Rebuild them with
`third_party/uspp-736/build_runatom.ps1`, which needs MSYS2 UCRT64 gfortran.
In auto-py-to-exe, under *Additional Files*, add:

- `third_party/uspp-736/bin/runatom.exe` with destination `uspp/`
- `third_party/uspp-736/bin/reform.exe` with destination `uspp/`
- `third_party/uspp-736/Work` with destination `uspp/Work` (the input library)
- `third_party/uspp-736/Source`, `third_party/uspp-736/Utils`,
  `third_party/uspp-736/patches`, `third_party/uspp-736/AA-README` and
  `third_party/uspp-736/Source/LICENSE` with destination `uspp/` (GPL v2
  requires shipping the source with the programs)

At run time each program is looked for in this order:
1. the path typed in the window (`reform.exe` next to it);
2. `%GQTEA_RUNATOM%` / `%GQTEA_REFORM%`;
3. `<bundle>/uspp/`;
4. the repository copy;
5. `PATH`.

Windows SmartScreen or antivirus software may question the unsigned programs;
they are not packed or obfuscated.
```

- [ ] **Step 4: Rewrite the user documentation**

In `USER_MANUAL.md`:
- rename the TOC entry and heading to *Vanderbilt USPP Generator (runatom)*;
- replace the section body with the workflow, fields and rules from `help_runatom` (Task 11 Step 4) written as prose, including *Using the pseudopotential in CPMD* (the `_formatted.uspp` file, `PP_LIBRARY_PATH`, the `*<file> FORMATTED` species line), plus a short *Limitations* list (keyee < 0, Dirac, UPF export, the library defects listed in `PROVENANCE.md`);
- update the two inventory lines at `USER_MANUAL.md:290` and `:488` from "Vanderbilt `runatom.x` input builders" to "Vanderbilt USPP generator (runs `runatom`)".

In `readme.md:37` and `:214`, likewise change "Vanderbilt `runatom.x` input builder"/"inputs" to "Vanderbilt USPP generator (runs `runatom`)".

- [ ] **Step 5: Add the project note to `CLAUDE.md`**

Add under *Tool-specific notes*:

```markdown
- **[runatom_input_gen.py](runatom_input_gen.py) + [runatomInput.py](runatomInput.py) / [runatomRunner.py](runatomRunner.py) / [runatomOutput.py](runatomOutput.py)** — **Inputs > Vanderbilt USPP Generator (runatom)**: prepares, checks and runs D. Vanderbilt's USPP 7.3.6 `runatom` (spec [docs/superpowers/specs/2026-10-28-runatom-integration-design.md](docs/superpowers/specs/2026-10-28-runatom-integration-design.md)). The Fortran is vendored **unmodified** in `third_party/uspp-736/` (`.gitattributes` `* -text`, `SOURCE_SHA256.txt` checked by a test) and built once into a static `bin/runatom.exe` by `build_runatom.ps1` (MSYS2 UCRT64), which comments out `external idate` in a temporary copy (otherwise gfortran fails to link: `idate_`). It runs as a **child process**, never in-process: 134 `call exit(1)` and several `stop` would kill the GUI, and its six file-name arguments are `character*40`, so `cwd` = job folder with relative names ≤ 40 characters. **Every input line is read with fixed columns and blanks ignored** ("1 0 0" in `i5` is 100, "2" in `f10.5` is 2e-5), so typed text is parsed and re-written by `runatomInput.format_adat`, never copied; `parse_adat` emulates gfortran's reader (round-trip test over all 418 readable library inputs; `023-V/V_ps.adat` is the upstream defect it must report on line 2). `check_case` mirrors `rwsubs.f`'s checks — notably GGA needs `ifqopt = 3` (−1 is refused too) — and `test_no_library_case_is_refused` guards against over-strict rules. Success of a stage = exit 0 **and** its stderr marker (`all electron calculation completed` / `writing the pseudopotential to file` / `pseudopotential test converged`) **and** its files: `stop` exits 0 on failure. Progress = stderr lines streamed by a reader thread; cancel kills the process; the window marshals progress with `call_soon_threadsafe`. `read_uspp` follows `Utils/read_uspp.f` (4-byte little-endian record markers; record 1 holds the run date, masked in byte comparisons). **CPMD output:** `reform.exe` (unmodified `Utils/reform.f`, built by the same script) writes `<symbol>_ps_formatted.uspp`, which CPMD reads with the species line `*<file> FORMATTED` (`ratom_utils`/`readvan`, checked in CPMD's source). The runner normalises its CRLF to LF, deletes an old output first (reform opens with `status='new'`) and passes names of at most 36 characters (`a36`). `read_formatted_uspp` re-reads it in `readvan`'s order as a check, never to write. Both programs ship in the bundle under `uspp/`. Validation: `tests/test_runatomExe.py` (skips without the exe) against Linux-gfortran references in `tests/data/runatom/` (regenerate only with `make_reference.sh` in WSL); the s channel of Pb is compared through invariants only (`Source/BENIGN_BUG`). Measured: every library case runs in < 1 s per stage.
```

- [ ] **Step 6: Run the full suite**

Run: `python -m pytest tests/ -q`
Expected: only the four baseline failures named in Global Constraints.

- [ ] **Step 7: Check the packaged build**

Build with auto-py-to-exe as documented, start the `.exe` and run the Ti example from the window. Also confirm that `find_runatom` reports `...\_internal\uspp\runatom.exe` in the window's program line. Record the outcome in the commit message.

- [ ] **Step 8: Commit (if authorised)**

```bash
git add readme.md USER_MANUAL.md CLAUDE.md tests/test_runatomVendor.py
git commit -m "docs(uspp): packaging, user manual and project notes for the USPP generator"
```
