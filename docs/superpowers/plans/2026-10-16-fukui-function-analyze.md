# Fukui Function Analysis (`fukuiFunctionAnalyze.py`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add **Tools > Fukui Function Analysis**. It reads three Gaussian `.fchk` files (N, N+1, N−1) of one molecule at one geometry and computes f⁺, f⁻, f⁰ and the dual descriptor on a grid and on the ρ(N) = 0.01 a.u. surface, with Multiwfn-rule extrema. It takes Hirshfeld-condensed indices from Multiwfn (run as an external engine, or imported from `CDFT.txt`) and computes the global CDFT indices. It writes a report and four cube files, and shows the functions in 3D.

**Architecture:**
- One new Toga-free helper, `cubeFile.py`.
- One new tool module, `fukuiFunctionAnalyze.py`. It holds a headless `FukuiAnalyser` and a Toga `FukuiFunctionAnalyzeUI`, the same split as `electrostaticSurfacePotential.py`.
- Small backward-compatible extensions to `gaussianWavefunction.py` (energy, method/basis, orbital energies) and to `espViewer.py` (an optional `caption=`).
- The real-space functions are computed natively. Condensed indices are **never** reimplemented: they come from Multiwfn main function 22, option 2.

**Tech Stack:** Python 3.13, numpy, scipy, toga-winforms 0.5.3, pytest. Existing modules reused: `gaussianWavefunction`, `molecularSurface`, `espViewer`/`molecularViewer`. Optional external programs: Multiwfn (tested with 2026.7.15) and Gaussian 09W, needed only to regenerate fixtures.

**Spec:** [docs/superpowers/specs/2026-10-16-fukui-function-design.md](../specs/2026-10-16-fukui-function-design.md). Read §1 (definitions), §2.2 (charge/multiplicity), §3.6 (what was verified) and §7 (tool decisions) first.

## Global Constraints

- **Branch.** Work on `feature/fukui-function`, created from `main` (Task 1). Run every command from `C:\AJCamargo\gqteaWinToga\venv\src` with the venv interpreter `..\Scripts\python.exe`.
- **What was executed before this plan was written (2026-10-16, scratch tree `$SCRATCH`):**
  - All of Tasks 1–7. Every non-UI line, and the stub-driven UI tests: 59 tests in `test_fukuiFunctionAnalyze.py` + `test_cubeFile.py` + the caption test, all passing. That includes the live Multiwfn-engine test, which reproduces `phenol_CDFT.txt` exactly.
  - The modified `gaussianWavefunction.py`/`espViewer.py` against the existing `test_gaussianWavefunction.py`, `test_espViewer.py` and `test_electrostaticSurfacePotential.py`, plus the two new field tests: 71 passed.
  - **Not executed:** building the real Toga window (`layout_main_window`), the live GL view, and the documentation edits (Task 8). Task 8 has the manual smoke test.
  - **Copy code exactly.** If something must change, change the test first and say why in the commit message.
- **Baseline: 4 pre-existing, unrelated failures.** Do not fix them here and do not count them as regressions: `tests/test_meanResidenceTime.py::TestAdvancedOptionsExposed::test_every_run_parameter_is_supplied_by_the_ui` and three table-styling tests in `tests/test_plotter.py`.
- **No new required dependency.**
  - Multiwfn is an optional external program, found through the typed path, `Multiwfnpath\Multiwfn.exe` or `PATH`. Missing or failing ⇒ a `NOTE:`, **never** an abort (the ESP tool's cubegen rule).
  - `cubeFile.py` must not import Toga (a subprocess test enforces it).
  - `fukuiFunctionAnalyze.py` may import Toga, like `electrostaticSurfacePotential.py`.
- **Units and defaults.**
  - Lengths in **bohr** inside the computation; Å only at the viewer and report boundary, via `coordinateConverter.BOHR_RADIUS_ANGSTROM`.
  - Fukui functions in a.u. (e/bohr³).
  - `HARTREE_TO_EV = 27.211386245988`.
  - Defaults: Fukui isovalue **0.005 a.u.**, mapping surface **ρ(N) = 0.01 a.u.**, grid **0.25 bohr**.
  - Hardness **η = I − A** (no ½), as in Multiwfn.
- **Multiwfn scripting rule (verified).**
  - Pass the input file on the command line (`[exe, "N.fch"]`) and copy `settings.ini` from the exe's folder into the working folder.
  - Only when no `settings.ini` exists, start the stdin script with one blank line (Multiwfn pauses for ENTER).
  - Do **not** rely on the `Multiwfnpath` environment variable for `settings.ini`: it did not work in this Windows build (spec §3.6).
- **UI errors go through dialogs** (`window.dialog(toga.ErrorDialog(...))`). Nothing raises out of an async handler.
- Commit messages end with: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`

**Scratch source of the verified files** (this session): `C:\Users\ajc\AppData\Local\Temp\claude\c--AJCamargo-gqteaWinToga-venv-src\dc50d304-21d9-487f-ba44-c3aa1d2a2200\scratchpad\impl\`, referred to below as `$SCRATCH`. Every source file is reproduced in full in this plan, and the plan text is authoritative if that folder is gone. The binary fixtures can be regenerated (Task 1).

## File Structure

| File | Status | Responsibility |
|---|---|---|
| `tests/data/fukui/phenol_{N,Np1,Nm1}.fchk`, `inputs/*.gjf` | create | G09W B3LYP/6-31G(d) single points at Multiwfn's `examples/phenol.xyz` geometry (≈1.8 MB) |
| `tests/data/fukui/phenol_CDFT.txt`, `make_reference.py`, `multiwfn_reference/phenol_grid_sample.txt` | create | Multiwfn 2026.7.15 reference: option 2 output, plus 300 sampled grid nodes of its f⁺/f⁻/f⁰/DD cubes |
| `tests/fukui_fixtures.py` | create | fixture paths (not a test module) |
| `gaussianWavefunction.py` | modify | store `total_energy`, `method`, `basis_set`, `n_alpha`, `n_beta`, `alpha_orbital_energies`, `other_densities` |
| `cubeFile.py` | create | Gaussian cube write/read (Toga-free) |
| `espViewer.py` | modify | optional `caption=` on `set_surface` (ESP behaviour unchanged) |
| `fukuiFunctionAnalyze.py` | create | formulas, state validation, `CDFT.txt` parser, Multiwfn engine, `FukuiAnalyser`, report/cubes, `FukuiFunctionAnalyzeUI` |
| `help.py` | modify | `HelpGqteaWin.help_fukui` |
| `gqteaWinToga.py` | modify | import + register under **Tools** after the ESP tool |
| `tests/test_cubeFile.py`, `tests/test_fukuiFunctionAnalyze.py` | create | tests |
| `tests/test_gaussianWavefunction.py`, `tests/test_espViewer.py` | modify | new cases |
| `USER_MANUAL.md`, `readme.md`, `CLAUDE.md` | modify | documentation |

---

### Task 1: Branch, fixtures and Multiwfn reference data

**Files:**
- Create: `tests/data/fukui/phenol_N.fchk`, `phenol_Np1.fchk`, `phenol_Nm1.fchk`
- Create: `tests/data/fukui/inputs/phenol_N.gjf`, `phenol_Np1.gjf`, `phenol_Nm1.gjf`
- Create: `tests/data/fukui/phenol_CDFT.txt`, `tests/data/fukui/make_reference.py`, `tests/data/fukui/multiwfn_reference/phenol_grid_sample.txt`
- Create: `tests/fukui_fixtures.py`

**Interfaces:**
- Produces: `fukui_fixtures.N`, `PLUS`, `MINUS` (Paths of the three `.fchk`), `CDFT` (Path of the Multiwfn output), `REFERENCE_DIR`. The sample-file columns are `x_bohr y_bohr z_bohr f+ f- f0 DD`, 300 rows, with `#` comments.

- [ ] **Step 1: Create the branch**

```bash
git checkout main && git checkout -b feature/fukui-function
```

- [ ] **Step 2: Copy the fixtures**

```bash
SCRATCH="/c/Users/ajc/AppData/Local/Temp/claude/c--AJCamargo-gqteaWinToga-venv-src/dc50d304-21d9-487f-ba44-c3aa1d2a2200/scratchpad/impl"
mkdir -p tests/data/fukui/inputs tests/data/fukui/multiwfn_reference
cp "$SCRATCH"/tests/data/fukui/phenol_*.fchk "$SCRATCH"/tests/data/fukui/phenol_CDFT.txt "$SCRATCH"/tests/data/fukui/make_reference.py tests/data/fukui/
cp "$SCRATCH"/tests/data/fukui/inputs/*.gjf tests/data/fukui/inputs/
cp "$SCRATCH"/tests/data/fukui/multiwfn_reference/phenol_grid_sample.txt tests/data/fukui/multiwfn_reference/
cp "$SCRATCH"/tests/fukui_fixtures.py tests/
ls -R tests/data/fukui
```

Expected: three `.fchk`, `phenol_CDFT.txt`, `make_reference.py`, `inputs/` with three `.gjf`, and `multiwfn_reference/phenol_grid_sample.txt` (302 lines).

**If the scratch folder no longer exists**, rebuild everything from the text below.

1. Write the three inputs. The coordinates are Multiwfn's `examples\phenol.xyz` (B3LYP/6-31G\* optimized).

`tests/data/fukui/inputs/phenol_N.gjf`:

```text
%chk=phenol_N.chk
%mem=500MB
#p B3LYP/6-31G(d) scf=tight

phenol N state

0 1
C       0.02082100     -1.85733400      0.00000000
C       1.22048400     -1.14011500      0.00000000
C       1.21812200      0.25309000      0.00000000
C       0.00000000      0.94117500      0.00000000
C      -1.20625800      0.23245100      0.00000000
C      -1.18949000     -1.16312800      0.00000000
H       0.03003900     -2.94340800      0.00000000
H       2.16996900     -1.66955400      0.00000000
H       2.14242000      0.82231900      0.00000000
H      -2.15282600      0.77102800      0.00000000
H      -2.13099000     -1.70650600      0.00000000
O       0.05161500      2.30900800      0.00000000
H      -0.85360600      2.65722900      0.00000000
```

2. `phenol_Np1.gjf` and `phenol_Nm1.gjf` are identical except for two things:
   - the `%chk=` name (`phenol_Np1.chk`, `phenol_Nm1.chk`) and the title (`phenol Np1 state`, `phenol Nm1 state`);
   - the charge/multiplicity line: **`-1 2`** for Np1 and **`1 2`** for Nm1.
3. Run them and format the checkpoints:

```bash
cd tests/data/fukui/inputs
export GAUSS_EXEDIR='C:\G09W'; export GAUSS_SCRDIR="$(cygpath -w "$TMP")"
for n in N Np1 Nm1; do /c/G09W/g09.exe phenol_$n.gjf phenol_$n.out; /c/G09W/formchk.exe phenol_$n.chk ../phenol_$n.fchk; done
grep -c "Normal termination" phenol_*.out   # expect 1 for each
rm -f phenol_*.chk phenol_*.out; cd ../../../..
```

4. Write `tests/data/fukui/make_reference.py` (below) and run `..\Scripts\python.exe tests\data\fukui\make_reference.py "<path to>\Multiwfn.exe"`.

`tests/data/fukui/make_reference.py`:

```python
"""Regenerate the Multiwfn reference data used by the Fukui tests.

NOT run by the test suite (it needs Multiwfn). Run it by hand, from the source
directory, only when a phenol_*.fchk fixture changes:

    python tests/data/fukui/make_reference.py path\\to\\Multiwfn.exe

It copies phenol_N/Np1/Nm1.fchk to N.fch/N+1.fch/N-1.fch in a temporary folder
and drives Multiwfn main function 22 (verified with Multiwfn 2026.7.15):
  option 2            -> CDFT.txt                  (saved as phenol_CDFT.txt)
  option 3, grid 1    -> low-quality grid; options 5-8 export f+.cub, f-.cub,
                         f0.cub, DD.cub
From the cubes it keeps 300 reproducible grid nodes and writes
multiwfn_reference/phenol_grid_sample.txt with columns
    x_bohr y_bohr z_bohr f+ f- f0 DD
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "multiwfn_reference"
SEED = 20261016
N_POINTS = 300
FILES = {"N.fch": "phenol_N.fchk", "N+1.fch": "phenol_Np1.fchk", "N-1.fch": "phenol_Nm1.fchk"}
PATHS = "N.fch\nN+1.fch\nN-1.fch\n"
SCRIPT = "22\n2\n" + PATHS + "3\n" + PATHS + "1\n5\n6\n7\n8\n0\n0\nq\n"  # N.fch is argv[1]


def read_cube(path: Path):
    lines = path.read_text().splitlines()
    head = lines[2].split()
    natoms, origin = int(head[0]), np.array([float(x) for x in head[1:4]])
    shape = tuple(int(lines[3 + k].split()[0]) for k in range(3))
    spacing = float(lines[3].split()[1])
    data = np.array([float(v) for l in lines[6 + natoms:] for v in l.split()]).reshape(shape)
    return origin, spacing, data


def main(multiwfn: str) -> None:
    settings = Path(multiwfn).resolve().parent / "settings.ini"
    lead = "" if settings.is_file() else "\n"  # without settings.ini Multiwfn pauses for ENTER
    with tempfile.TemporaryDirectory() as tmp:
        for target, source in FILES.items():
            shutil.copyfile(HERE / source, Path(tmp) / target)
        if settings.is_file():
            shutil.copyfile(settings, Path(tmp) / "settings.ini")
        subprocess.run([multiwfn, "N.fch"], input=lead + SCRIPT, text=True, cwd=tmp, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.copyfile(Path(tmp) / "CDFT.txt", HERE / "phenol_CDFT.txt")
        cubes = {k: read_cube(Path(tmp) / f"{k}.cub") for k in ("f+", "f-", "f0", "DD")}
    origin, spacing, first = cubes["f+"]
    rng = np.random.default_rng(SEED)
    idx = np.column_stack([rng.integers(n, size=N_POINTS) for n in first.shape])
    xyz = origin + spacing * idx
    OUT.mkdir(exist_ok=True)
    with (OUT / "phenol_grid_sample.txt").open("w") as handle:
        handle.write("# Multiwfn 22/3 low grid, options 5-8, sampled at 300 grid nodes\n")
        handle.write("# x_bohr y_bohr z_bohr f+ f- f0 DD   (a.u.)\n")
        for (i, j, k), (x, y, z) in zip(idx, xyz):
            values = [cubes[key][2][i, j, k] for key in ("f+", "f-", "f0", "DD")]
            handle.write(f"{x:.6f} {y:.6f} {z:.6f} " + " ".join(f"{v:.5e}" for v in values) + "\n")
    print("wrote phenol_CDFT.txt and multiwfn_reference/phenol_grid_sample.txt")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else shutil.which("Multiwfn.exe") or "Multiwfn.exe")
```

`tests/fukui_fixtures.py`:

```python
"""Paths of the Fukui test fixtures (not a test module).

tests/data/fukui/phenol_{N,Np1,Nm1}.fchk: G09W B3LYP/6-31G(d) scf=tight single points
at Multiwfn's examples/phenol.xyz geometry, charge/multiplicity 0 1, -1 2, 1 2
(inputs/*.gjf). phenol_CDFT.txt: Multiwfn 2026.7.15 main function 22 option 2 on
those three files. multiwfn_reference/*.cub: Multiwfn option 3 (low grid), options
5-8, on the same files -- regenerate with tests/data/fukui/make_reference.py.
"""
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data" / "fukui"
N, PLUS, MINUS = (DATA / f"phenol_{s}.fchk" for s in ("N", "Np1", "Nm1"))
CDFT = DATA / "phenol_CDFT.txt"
REFERENCE_DIR = DATA / "multiwfn_reference"
```

- [ ] **Step 3: Sanity-check the reference**

Run:

```bash
head -18 tests/data/fukui/phenol_CDFT.txt | tail -14
```

Expected: the header `Atom q(N) q(N+1) q(N-1) f- f+ f0 CDD`, then 13 rows. The first is `1(C )  -0.0587  -0.1185   0.0852   0.1439   0.0598   0.1018  -0.0841`, identical to the Multiwfn manual §4.22.1 and §4.7.3.

- [ ] **Step 4: Commit**

```bash
git add tests/data/fukui tests/fukui_fixtures.py
git commit -m "Fukui fixtures: phenol N/N+1/N-1 fchk (G09W B3LYP/6-31G(d)) and Multiwfn 2026.7.15 reference

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: `gaussianWavefunction` keeps energies, method, basis and orbital energies

**Files:**
- Modify: `gaussianWavefunction.py` (the `FchkWavefunction` dataclass and the end of `read_fchk`)
- Test: `tests/test_gaussianWavefunction.py` (append a class)

**Interfaces:**
- Produces new `FchkWavefunction` fields, all with defaults so existing constructors keep working:
  - `total_energy: Optional[float]` (Hartree, from `Total Energy`)
  - `method: str` (fchk line 2, token 2, e.g. `"RB3LYP"`)
  - `basis_set: str` (token 3, e.g. `"6-31G(d)"`)
  - `n_alpha: int`, `n_beta: int`
  - `alpha_orbital_energies: Optional[np.ndarray]` (Hartree)
  - `other_densities: List[str]`: every `Total … Density` section except `Total SCF Density`

- [ ] **Step 1: Write the failing test.** Append to `tests/test_gaussianWavefunction.py`:

```python


class TestCdftFields:
    """Fields added for the Fukui tool: energy, method/basis, electron counts, orbital energies."""

    def test_energy_method_basis_and_orbitals(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))
        assert w.method == "RB3LYP" and w.basis_set == "6-31G"
        assert w.total_energy == pytest.approx(-78.572035869, abs=1e-8)
        assert (w.n_alpha, w.n_beta) == (8, 8) and w.other_densities == []
        assert w.alpha_orbital_energies[w.n_alpha - 1] == pytest.approx(-0.26819427, abs=1e-8)

    def test_open_shell_counts(self):
        w = gw.read_fchk(fchk_path("ch3_ub3lyp_631gd"))
        assert w.method == "UB3LYP" and (w.n_alpha, w.n_beta) == (5, 4)
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `..\Scripts\python.exe -m pytest tests/test_gaussianWavefunction.py -q -k CdftFields`
Expected: FAIL with `AttributeError: 'FchkWavefunction' object has no attribute 'method'`.

- [ ] **Step 3: Implement.** In `gaussianWavefunction.py`, replace

```python
    dipole: Optional[np.ndarray] = None
    transform: np.ndarray = field(default=None, repr=False)     # fchk -> Cartesian
    cart_density: np.ndarray = field(default=None, repr=False)  # T^T P T
```

with

```python
    dipole: Optional[np.ndarray] = None
    transform: np.ndarray = field(default=None, repr=False)     # fchk -> Cartesian
    cart_density: np.ndarray = field(default=None, repr=False)  # T^T P T
    total_energy: Optional[float] = None        # Hartree, "Total Energy"
    method: str = ""                            # line 2, e.g. "RB3LYP"
    basis_set: str = ""                         # line 2, e.g. "6-31G(d)"
    n_alpha: int = 0
    n_beta: int = 0
    alpha_orbital_energies: Optional[np.ndarray] = None  # Hartree
    other_densities: List[str] = field(default_factory=list)  # e.g. ["Total MP2 Density"]
```

then replace

```python
    lines = text.splitlines()
    wfn = FchkWavefunction(
```

with

```python
    lines = text.splitlines()
    line2 = lines[1].split() if len(lines) > 1 else []
    energies = sections.get("Alpha Orbital Energies")
    wfn = FchkWavefunction(
```

and replace

```python
        dipole=np.asarray(sections["Dipole Moment"]) if "Dipole Moment" in sections else None,
    )
```

with

```python
        dipole=np.asarray(sections["Dipole Moment"]) if "Dipole Moment" in sections else None,
        total_energy=float(sections["Total Energy"]) if "Total Energy" in sections else None,
        method=line2[1] if len(line2) > 1 else "",
        basis_set=line2[2] if len(line2) > 2 else "",
        n_alpha=int(sections.get("Number of alpha electrons", 0)),
        n_beta=int(sections.get("Number of beta electrons", 0)),
        alpha_orbital_energies=np.asarray(energies, dtype=float) if energies is not None else None,
        other_densities=sorted(k for k in sections
                               if k.startswith("Total ") and k.endswith(" Density") and k != DENSITY_SECTION),
    )
```

- [ ] **Step 4: Run the fchk, ESP and viewer tests**

Run: `..\Scripts\python.exe -m pytest tests/test_gaussianWavefunction.py tests/test_electrostaticSurfacePotential.py -q`
Expected: all pass (this was 71 tests together with `test_espViewer.py` in the scratch run).

- [ ] **Step 5: Commit**

```bash
git add gaussianWavefunction.py tests/test_gaussianWavefunction.py
git commit -m "gaussianWavefunction: keep total energy, method/basis, electron counts and alpha orbital energies

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: `cubeFile.py`, Gaussian cube writer and reader

**Files:**
- Create: `cubeFile.py`
- Test: `tests/test_cubeFile.py`

**Interfaces:**
- Produces:
  - `CubeData(title, comment, atomic_numbers, nuclear_charges, coordinates_bohr, origin_bohr, spacing_bohr, values[nx, ny, nz])`
  - `format_cube(cube) -> str`
  - `write_cube(path, cube) -> str`: raises `ValueError` on `OSError`
  - `read_cube(path) -> CubeData`: raises `ValueError` on malformed, skewed, orbital (negative atom count) or truncated files

- [ ] **Step 1: Write the failing test** `tests/test_cubeFile.py`:

```python
"""Tests for the Gaussian cube writer/reader."""
import subprocess
import sys

import numpy as np
import pytest

import cubeFile as cf


def sample_cube(shape=(3, 4, 7)):
    rng = np.random.default_rng(7)
    return cf.CubeData("title line", "comment line", np.array([8, 1]), np.array([8.0, 1.0]),
                       np.array([[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]]), np.array([-4.0, -3.5, -2.0]),
                       0.25, rng.normal(scale=0.01, size=shape))


class TestCube:
    def test_round_trip_keeps_five_significant_digits(self, tmp_path):
        cube = sample_cube()
        back = cf.read_cube(cf.write_cube(tmp_path / "a.cub", cube))
        assert back.values.shape == cube.values.shape and back.spacing == 0.25
        assert np.allclose(back.origin, cube.origin) and np.array_equal(back.atomic_numbers, [8, 1])
        assert np.allclose(back.values, cube.values, rtol=1e-5, atol=1e-12)

    def test_layout_is_six_per_line_and_restarts_each_column(self):
        text = cf.format_cube(sample_cube((1, 2, 7)))
        data = text.splitlines()[6 + 2:]
        assert [len(line.split()) for line in data] == [6, 1, 6, 1]
        assert text.splitlines()[3].split() == ["1", "0.250000", "0.000000", "0.000000"]

    def test_truncated_file_is_a_value_error(self, tmp_path):
        path = tmp_path / "t.cub"
        path.write_text("\n".join(cf.format_cube(sample_cube()).splitlines()[:-2]) + "\n")
        with pytest.raises(ValueError, match="truncated"):
            cf.read_cube(path)

    def test_skewed_grid_is_rejected(self, tmp_path):
        lines = cf.format_cube(sample_cube()).splitlines()
        lines[4] = "    4    0.100000    0.250000    0.000000"
        path = tmp_path / "s.cub"
        path.write_text("\n".join(lines) + "\n")
        with pytest.raises(ValueError, match="axis-aligned"):
            cf.read_cube(path)

    def test_unwritable_path_is_a_value_error(self, tmp_path):
        (tmp_path / "file").write_text("x")
        with pytest.raises(ValueError, match="Could not write"):
            cf.write_cube(tmp_path / "file" / "a.cub", sample_cube())

    def test_module_does_not_import_toga(self):
        code = "import sys, cubeFile; print('toga' in sys.modules)"
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                             cwd=str(__import__('pathlib').Path(cf.__file__).parent), check=True)
        assert out.stdout.strip() == "False"
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `..\Scripts\python.exe -m pytest tests/test_cubeFile.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'cubeFile'`.

- [ ] **Step 3: Implement** `cubeFile.py`:

```python
"""Gaussian cube files (.cub/.cube): write and read volumetric grid data.

Toga-free. Lengths are bohr, as in the cube format itself (a positive voxel
count on lines 4-6 means bohr). Only orthogonal grids whose three axes are
the Cartesian axes with one common spacing are written or read -- that is
what every grid in this project is, and what Multiwfn writes by default.
Values are written 6 per line in %13.5E, z fastest, a new line for every
(x, y) column: the layout Gaussian's cubegen and Multiwfn use.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

AXIS_TOLERANCE = 1e-6  # bohr


@dataclass
class CubeData:
    title: str
    comment: str
    atomic_numbers: np.ndarray   # (n,)
    nuclear_charges: np.ndarray  # (n,)
    coordinates: np.ndarray      # (n, 3) bohr
    origin: np.ndarray           # (3,) bohr
    spacing: float               # bohr
    values: np.ndarray           # (nx, ny, nz)


def format_cube(cube: CubeData) -> str:
    nx, ny, nz = cube.values.shape
    out = [cube.title.replace("\n", " "), cube.comment.replace("\n", " "),
           f"{len(cube.atomic_numbers):5d}" + "".join(f"{c:12.6f}" for c in cube.origin)]
    for n, axis in zip((nx, ny, nz), np.eye(3) * cube.spacing):
        out.append(f"{n:5d}" + "".join(f"{c:12.6f}" for c in axis))
    for z, q, xyz in zip(cube.atomic_numbers, cube.nuclear_charges, cube.coordinates):
        out.append(f"{int(z):5d}{float(q):12.6f}" + "".join(f"{c:12.6f}" for c in xyz))
    for column in np.asarray(cube.values, dtype=float).reshape(nx * ny, nz):
        for k in range(0, nz, 6):
            out.append("".join(f"{v:13.5E}" for v in column[k:k + 6]))
    return "\n".join(out) + "\n"


def write_cube(path, cube: CubeData) -> str:
    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(format_cube(cube), encoding="ascii")
    except OSError as exc:
        raise ValueError(f"Could not write {path}: {exc}") from exc
    return str(path)


def read_cube(path) -> CubeData:
    """Parse a cube written by this module, cubegen or Multiwfn. ValueError on any problem."""
    path = Path(path)
    try:
        lines = path.read_text(encoding="ascii", errors="replace").splitlines()
    except OSError as exc:
        raise ValueError(f"Cannot read {path}: {exc}") from exc
    try:
        head = lines[2].split()
        natoms = int(head[0])
        origin = np.array([float(x) for x in head[1:4]])
        shape, axes = [], []
        for k in range(3):
            t = lines[3 + k].split()
            shape.append(int(t[0]))
            axes.append([float(x) for x in t[1:4]])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"{path.name}: malformed cube header.") from exc
    if natoms < 0:
        raise ValueError(f"{path.name}: orbital cubes (negative atom count) are not supported.")
    axes = np.array(axes)
    if min(shape) <= 0:
        raise ValueError(f"{path.name}: grid sizes must be positive (bohr units).")
    spacing = float(axes[0, 0])
    if np.abs(axes - np.eye(3) * spacing).max() > AXIS_TOLERANCE:
        raise ValueError(f"{path.name}: only axis-aligned grids with one spacing are supported.")
    atoms = [lines[6 + i].split() for i in range(natoms)]
    try:
        z = np.array([int(a[0]) for a in atoms])
        q = np.array([float(a[1]) for a in atoms])
        xyz = np.array([[float(c) for c in a[2:5]] for a in atoms]).reshape(-1, 3)
        data = np.array([float(v) for line in lines[6 + natoms:] for v in line.split()])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"{path.name}: malformed atom or data lines.") from exc
    expected = shape[0] * shape[1] * shape[2]
    if data.size != expected:
        raise ValueError(f"{path.name}: {data.size} values, expected {expected}; the file looks truncated.")
    return CubeData(lines[0], lines[1], z, q, xyz, origin, spacing, data.reshape(shape))
```

- [ ] **Step 4: Run the test**

Run: `..\Scripts\python.exe -m pytest tests/test_cubeFile.py -q`
Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
git add cubeFile.py tests/test_cubeFile.py
git commit -m "cubeFile: Toga-free Gaussian cube writer/reader (Multiwfn/cubegen layout)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: `espViewer`: optional caption for properties other than the ESP

**Files:**
- Modify: `espViewer.py` (`SurfaceDrawData`, `ESPSurfaceViewer.set_surface`, `ESPSurfaceViewer.colour_bar_text`)
- Test: `tests/test_espViewer.py` (append one function)

**Interfaces:**
- Produces: `ESPSurfaceViewer.set_surface(..., half_range=None, caption=None)`. A non-empty `caption` replaces the ESP caption in the colour bar. Omitting it gives exactly the old behaviour.

- [ ] **Step 1: Write the failing test.** Append to `tests/test_espViewer.py`:

```python
def test_custom_caption_replaces_the_esp_caption():
    v = ev.ESPSurfaceViewer()
    tri = np.array([[0, 1, 2]])
    verts = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    v.set_surface(verts, tri, np.tile([0, 0, 1.0], (3, 1)), [0.01, -0.01, 0.0], 0.01, caption="f- (a.u.)")
    assert v.colour_bar_text()[2] == "f- (a.u.)"
    v.set_surface(verts, tri, np.tile([0, 0, 1.0], (3, 1)), [0.01, -0.01, 0.0], 0.01)
    assert v.colour_bar_text()[2].startswith("ESP (a.u.)")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `..\Scripts\python.exe -m pytest tests/test_espViewer.py -q -k custom_caption`
Expected: FAIL with `TypeError: ... unexpected keyword argument 'caption'`.

- [ ] **Step 3: Implement.** In `espViewer.py`:

(a) In `SurfaceDrawData`, after `radius: float          # max distance of a vertex from the atoms' centre`, add

```python
    caption: Optional[str] = None  # replaces the ESP caption (other mapped properties)
```

(b) Change the signature line of `set_surface` from

```python
                    markers: Sequence[SurfaceMarker] = (), half_range: Optional[float] = None):
```

to

```python
                    markers: Sequence[SurfaceMarker] = (), half_range: Optional[float] = None,
                    caption: Optional[str] = None):
```

(c) In `set_surface`, after `radius=float(np.linalg.norm(verts - centre, axis=1).max()) if len(verts) else 0.0,` add `caption=caption,` inside the `SurfaceDrawData(...)` call.

(d) Replace the body of `colour_bar_text` (the method):

```python
        s = self.get_surface()
        return None if s is None else colour_bar_text(s.half_range, self.colour_unit, s.isovalue)
```

with

```python
        s = self.get_surface()
        if s is None:
            return None
        left, right, caption = colour_bar_text(s.half_range, self.colour_unit, s.isovalue)
        return left, right, s.caption or caption
```

- [ ] **Step 4: Run the viewer and ESP tests**

Run: `..\Scripts\python.exe -m pytest tests/test_espViewer.py tests/test_electrostaticSurfacePotential.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add espViewer.py tests/test_espViewer.py
git commit -m "espViewer: optional colour-bar caption so other mapped properties can reuse the viewer

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Fukui core: formulas, three-state validation, `CDFT.txt` parser, Multiwfn engine

**Files:**
- Create: `fukuiFunctionAnalyze.py` (first part: module docstring through the Multiwfn runner)
- Test: `tests/test_fukuiFunctionAnalyze.py` (first part)

**Interfaces:**
- Consumes: `gaussianWavefunction.read_fchk` and the Task 2 fields.
- Produces (module level):
  - constants `HARTREE_TO_EV`, `TCE_HOMO_HARTREE`, `DEGENERACY_EV`, `DEFAULT_FUKUI_ISOVALUE`, `DEFAULT_SURFACE_ISOVALUE`, `DEFAULT_SPACING_BOHR`, `FUNCTIONS = ("f+", "f-", "f0", "DD")`, `FUNCTION_TITLES`, `CUBE_SUFFIXES`, `REPORT_SUFFIX`, and `CONDENSED_MULTIWFN` / `CONDENSED_IMPORT` / `CONDENSED_NONE` / `CONDENSED_SOURCES`;
  - `parse_positive_float(text, default, name, upper)`, `resolve_output_dir(fchk, chosen)`;
  - `strip_reference_prefix(method)`, `same_basis(a, b)`, `validate_states(n, plus, minus) -> (p, q, notes)`, `frontier_degeneracy(w) -> (p, q) | None`, `has_diffuse_functions(basis)`;
  - `fukui_from_densities(rho_n, rho_plus, rho_minus, p=1, q=1) -> {"f+","f-","f0","DD"}`;
  - `GlobalIndices(e_n, e_plus, e_minus, homo_n)`, with properties `ip, ea, chi, mu, eta, softness, omega, nucleophilicity`; `homo_energy(w)`;
  - `CondensedRow(index, symbol, q_n, q_plus, q_minus, f_minus, f_plus, f_zero, cdd)`, `parse_cdft_text(text) -> (rows, globals_dict)`, `condensed_inconsistencies(rows) -> [1-based atoms]`, `check_rows_match(rows, wfn)`;
  - `find_multiwfn(explicit=None)`, `multiwfn_cdft_script(settings_found)`, `run_multiwfn_cdft(exe, paths, timeout) -> CDFT.txt text`, which raises `RuntimeError`.

- [ ] **Step 1: Write the failing test** `tests/test_fukuiFunctionAnalyze.py` (first part; later tasks append to it):

```python
"""Tests for the Fukui function analyser, Multiwfn interop and the UI plumbing."""
import asyncio
import inspect
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

import fukuiFunctionAnalyze as fk
import gaussianWavefunction as gw
from fukui_fixtures import CDFT, MINUS, N, PLUS, REFERENCE_DIR

CDFT_TEXT = CDFT.read_text()


class TestFormulas:
    def test_fukui_definitions(self):
        f = fk.fukui_from_densities(np.array([1.0]), np.array([1.4]), np.array([0.7]))
        assert f["f+"][0] == pytest.approx(0.4) and f["f-"][0] == pytest.approx(0.3)
        assert f["f0"][0] == pytest.approx(0.35) and f["DD"][0] == pytest.approx(0.1)

    def test_degenerate_states_divide_by_p_and_q(self):
        f = fk.fukui_from_densities(np.array([1.0]), np.array([1.6]), np.array([0.4]), p=2, q=3)
        assert f["f+"][0] == pytest.approx(0.3) and f["f-"][0] == pytest.approx(0.2)

    def test_global_indices_reproduce_multiwfn(self):
        # E(N), E(N+1), E(N-1), E_HOMO(N) and the expected values: phenol_CDFT.txt
        g = fk.GlobalIndices(-307.464870, -307.383627, -307.163447, -0.218914)
        assert g.ip == pytest.approx(0.301423, abs=2e-6)
        assert g.ea == pytest.approx(-0.081243, abs=2e-6)
        assert g.chi == pytest.approx(0.110090, abs=2e-6) and g.mu == pytest.approx(-0.110090, abs=2e-6)
        assert g.eta == pytest.approx(0.382666, abs=2e-6)
        assert g.softness == pytest.approx(2.613246, abs=2e-5)
        assert g.omega == pytest.approx(0.015836, abs=2e-6)
        assert g.nucleophilicity == pytest.approx(0.116284, abs=2e-6)
        assert g.eta * fk.HARTREE_TO_EV == pytest.approx(10.4129, abs=1e-4)

    def test_non_positive_gap_gives_nan_not_an_exception(self):
        g = fk.GlobalIndices(-1.0, -1.2, -0.9)  # IP 0.1 < EA 0.2
        assert g.eta < 0 and np.isnan(g.softness) and np.isnan(g.omega)

    @pytest.mark.parametrize("method,expected", [("RB3LYP", "B3LYP"), ("UB3LYP", "B3LYP"),
                                                 ("ROB3LYP", "B3LYP"), ("RHF", "HF"), ("UHF", "HF")])
    def test_reference_prefix(self, method, expected):
        assert fk.strip_reference_prefix(method) == expected

    @pytest.mark.parametrize("basis,diffuse", [("6-31G(d)", False), ("6-311+G(d,p)", True),
                                               ("aug-cc-pVTZ", True), ("def2-TZVP", False)])
    def test_diffuse_basis_detection(self, basis, diffuse):
        assert fk.has_diffuse_functions(basis) is diffuse


class TestParsing:
    def test_blank_is_default(self):
        assert fk.parse_positive_float("", 0.005, "x", 1.0) == 0.005

    @pytest.mark.parametrize("text", ["abc", "0", "-1", "2"])
    def test_bad_values(self, text):
        with pytest.raises(ValueError):
            fk.parse_positive_float(text, 0.005, "Fukui isovalue", 1.0)


class TestStates:
    @pytest.fixture(scope="class")
    def states(self):
        return gw.read_fchk(N), gw.read_fchk(PLUS), gw.read_fchk(MINUS)

    def test_new_fchk_fields(self, states):
        n, plus, _ = states
        assert n.method == "RB3LYP" and plus.method == "UB3LYP" and n.basis_set == "6-31G(d)"
        assert n.total_energy == pytest.approx(-307.464870, abs=1e-6)
        assert (n.n_alpha, n.n_beta) == (25, 25) and n.other_densities == []
        assert n.alpha_orbital_energies[n.n_alpha - 1] == pytest.approx(-0.218914, abs=1e-6)

    def test_phenol_triplet_is_valid(self, states):
        p, q, notes = fk.validate_states(*states)
        assert (p, q, notes) == (1, 1, [])

    def test_swapped_files_are_explained(self, states):
        n, plus, minus = states
        with pytest.raises(ValueError, match="swapped"):
            fk.validate_states(n, minus, plus)

    def test_moved_geometry_is_rejected(self, states):
        n, plus, minus = states
        moved = replace(plus, coordinates=plus.coordinates + 0.01)
        with pytest.raises(ValueError, match="single points"):
            fk.validate_states(n, moved, minus)

    def test_other_method_is_rejected(self, states):
        n, plus, minus = states
        with pytest.raises(ValueError, match="same method"):
            fk.validate_states(n, replace(plus, method="UPBE0"), minus)

    def test_other_basis_is_rejected(self, states):
        n, plus, minus = states
        with pytest.raises(ValueError, match="basis set"):
            fk.validate_states(n, replace(plus, shells=plus.shells[:-1]), minus)

    def test_unusual_multiplicity_is_a_note(self, states):
        n, plus, minus = states
        _, _, notes = fk.validate_states(n, replace(plus, multiplicity=4), minus)
        assert "multiplicity 4" in notes[0]

    def test_phenol_frontier_orbitals_are_not_degenerate(self, states):
        assert fk.frontier_degeneracy(states[0]) == (1, 1)
        assert fk.frontier_degeneracy(states[1]) is None  # open shell


class TestCDFT:
    def test_parse_condensed_table_and_globals(self):
        rows, g = fk.parse_cdft_text(CDFT_TEXT)
        assert len(rows) == 13 and rows[0].symbol == "C" and rows[11].symbol == "O"
        c1 = rows[0]
        assert (c1.q_n, c1.q_plus, c1.q_minus, c1.f_minus, c1.f_plus, c1.f_zero, c1.cdd) == \
               (-0.0587, -0.1185, 0.0852, 0.1439, 0.0598, 0.1018, -0.0841)
        assert g["e_n"] == -307.464870 and g["eta"] == 0.382666 and g["ip"] == 0.301423

    def test_older_vertical_ip_label_is_accepted(self):
        _, g = fk.parse_cdft_text(CDFT_TEXT.replace("First vertical IP", "Vertical IP"))
        assert g["ip"] == 0.301423

    def test_phenol_rows_are_self_consistent(self):
        rows, _ = fk.parse_cdft_text(CDFT_TEXT)
        assert fk.condensed_inconsistencies(rows) == []

    def test_inconsistent_row_is_found(self):
        rows, _ = fk.parse_cdft_text(CDFT_TEXT)
        rows[2] = replace(rows[2], f_plus=rows[2].f_plus + 0.01)
        assert fk.condensed_inconsistencies(rows) == [3]

    def test_not_a_cdft_file(self):
        with pytest.raises(ValueError, match="main function 22"):
            fk.parse_cdft_text("hello\nworld\n")

    def test_rows_must_match_the_molecule(self):
        rows, _ = fk.parse_cdft_text(CDFT_TEXT)
        n = gw.read_fchk(N)
        fk.check_rows_match(rows, n)
        with pytest.raises(ValueError, match="lists 12 atoms"):
            fk.check_rows_match(rows[:-1], n)

    def test_script_adds_a_blank_line_only_without_settings(self):
        assert fk.multiwfn_cdft_script(True).startswith("22\n2\nN.fch\nN+1.fch\nN-1.fch\n")
        assert fk.multiwfn_cdft_script(False).startswith("\n22\n")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `..\Scripts\python.exe -m pytest tests/test_fukuiFunctionAnalyze.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'fukuiFunctionAnalyze'`.

- [ ] **Step 3: Implement.** Create `fukuiFunctionAnalyze.py` with exactly this content. It imports Toga and the viewer at the top, as the ESP tool does; those are only used by Task 7's UI.

```python
"""Tools > Fukui Function Analysis.

Fukui functions f+(r), f-(r), f0(r) and the dual descriptor Df(r) from three
Gaussian formatted checkpoints of the same molecule at the same geometry and
level of theory: the N-electron state and the N+p / N-q states (p = q = 1
normally; p, q > 1 for (quasi-)degenerate frontier orbitals, Multiwfn manual
section 3.25.4).

    f+  = [rho(N+p) - rho(N)] / p        nucleophilic attack
    f-  = [rho(N) - rho(N-q)] / q        electrophilic attack
    f0  = (f+ + f-) / 2                  radical attack
    Df  = f+ - f-                        dual descriptor (> 0 nucleophilic,
                                         < 0 electrophilic attack)

Real-space functions are computed here, from the exact densities
(gaussianWavefunction). Atom-condensed indices come from Multiwfn main
function 22 option 2 (Hirshfeld charges of the three states), either run here
as an external engine or imported from an existing CDFT.txt -- this module
does not compute Hirshfeld charges itself. Global indices follow Multiwfn's
conventions (hardness = I - A, without the 1/2).

These quantities indicate *likely* reactive regions; they are not
predictions of reaction outcome.
"""
from __future__ import annotations

import asyncio
import os
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import cubeFile as cf
import gaussianWavefunction as gw
import molecularSurface as ms
from coordinateConverter import BOHR_RADIUS_ANGSTROM
from espViewer import ESPSurfaceViewer, SurfaceMarker, UNIT_AU
from help import AtomicData, HelpGqteaWin
from molecularViewer import save_canvas_image

HARTREE_TO_EV = 27.211386245988   # CODATA 2018; reproduces Multiwfn's eV columns
TCE_HOMO_HARTREE = -0.335198      # E_HOMO(TCE), B3LYP/6-31G*, as used by Multiwfn (manual 3.25.2)
DEGENERACY_EV = 0.01              # Multiwfn's automatic FMO degeneracy threshold (manual 4.22.3.1)
DEFAULT_FUKUI_ISOVALUE = 0.005    # a.u.; the manual's examples use 0.005-0.01
DEFAULT_SURFACE_ISOVALUE = 0.01   # rho isosurface for mapping f (manual 4.12.4)
DEFAULT_SPACING_BOHR = 0.25
COORDINATE_TOLERANCE_BOHR = 1e-5
CONDENSED_TOLERANCE = 2e-4        # CDFT.txt prints 4 decimals
MULTIWFN_TIMEOUT_S = 1800

FUNCTIONS = ("f+", "f-", "f0", "DD")
FUNCTION_TITLES = {
    "f+": "f+  (nucleophilic attack)",
    "f-": "f-  (electrophilic attack)",
    "f0": "f0  (radical attack)",
    "DD": "dual descriptor  (+ nucleophilic / - electrophilic)",
}
CUBE_SUFFIXES = {"f+": "_fplus.cub", "f-": "_fminus.cub", "f0": "_fzero.cub", "DD": "_dual.cub"}
REPORT_SUFFIX = "_Fukui.txt"

CONDENSED_MULTIWFN = "Run Multiwfn (Hirshfeld)"
CONDENSED_IMPORT = "Import Multiwfn CDFT.txt"
CONDENSED_NONE = "Skip condensed indices"
CONDENSED_SOURCES = (CONDENSED_MULTIWFN, CONDENSED_IMPORT, CONDENSED_NONE)

ProgressCallback = Callable[[float, str], None]


# ----------------------------------------------------------------------
# Input parsing
# ----------------------------------------------------------------------
def parse_positive_float(text: Optional[str], default: float, name: str, upper: float) -> float:
    """Blank -> default; otherwise a number in (0, upper)."""
    text = (text or "").strip()
    if not text:
        return default
    try:
        value = float(text)
    except ValueError:
        raise ValueError(f"{name}: '{text}' is not a number.") from None
    if not 0.0 < value < upper:
        raise ValueError(f"{name} must be greater than 0 and smaller than {upper:g}.")
    return value


def resolve_output_dir(fchk_path: str, chosen: Optional[str]) -> str:
    chosen = (chosen or "").strip()
    return chosen if chosen else str(Path(fchk_path).resolve().parent)


# ----------------------------------------------------------------------
# Three states: validation
# ----------------------------------------------------------------------
def strip_reference_prefix(method: str) -> str:
    """'RB3LYP', 'UB3LYP', 'ROB3LYP' -> 'B3LYP' (the closed- and open-shell states differ only there)."""
    m = method.upper()
    for prefix in ("RO", "R", "U"):
        if m.startswith(prefix) and len(m) > len(prefix):
            return m[len(prefix):]
    return m


def same_basis(a: gw.FchkWavefunction, b: gw.FchkWavefunction) -> bool:
    if len(a.shells) != len(b.shells):
        return False
    return all(sa.l == sb.l and sa.pure == sb.pure and sa.atom == sb.atom
               and sa.exponents.shape == sb.exponents.shape
               and np.allclose(sa.exponents, sb.exponents, rtol=1e-10, atol=0.0)
               for sa, sb in zip(a.shells, b.shells))


def validate_states(n: gw.FchkWavefunction, plus: gw.FchkWavefunction,
                    minus: gw.FchkWavefunction) -> Tuple[int, int, List[str]]:
    """(p, q, notes). ValueError when the three files cannot form a finite difference."""
    for label, w in (("N+1", plus), ("N-1", minus)):
        if w.n_atoms != n.n_atoms or not np.array_equal(w.atomic_numbers, n.atomic_numbers):
            raise ValueError(f"The {label} file holds different atoms from the N file; "
                             "all three files must describe the same molecule.")
        shift = float(np.abs(w.coordinates - n.coordinates).max())
        if shift > COORDINATE_TOLERANCE_BOHR:
            raise ValueError(
                f"The {label} geometry differs from the N geometry by up to "
                f"{shift * BOHR_RADIUS_ANGSTROM:.2e} A. All three states must use the optimized "
                "N-state geometry (the derivative is taken at fixed nuclei): run single points "
                "for the ions, do not re-optimize them.")
        if not same_basis(w, n):
            raise ValueError(f"The {label} file uses a different basis set from the N file.")
        if strip_reference_prefix(w.method) != strip_reference_prefix(n.method):
            raise ValueError(f"The {label} file was computed with {w.method}, the N file with "
                             f"{n.method}; all three states need the same method.")
    p = plus.n_electrons - n.n_electrons
    q = n.n_electrons - minus.n_electrons
    if p <= 0 or q <= 0:
        raise ValueError(
            f"Electron counts N = {n.n_electrons}, 'N+1' = {plus.n_electrons}, "
            f"'N-1' = {minus.n_electrons}: the N+1 file must hold more electrons than N and the "
            "N-1 file fewer. Were two files swapped?")
    notes: List[str] = []
    for label, w, k in (("N+1" if p == 1 else f"N+{p}", plus, p), ("N-1" if q == 1 else f"N-{q}", minus, q)):
        if abs(w.multiplicity - n.multiplicity) != k:
            notes.append(f"The {label} state has multiplicity {w.multiplicity} (N state: {n.multiplicity}). "
                         "Usually adding/removing electrons into/from a non-degenerate orbital changes "
                         f"the multiplicity by {k}; confirm this is the lowest-energy spin state.")
    for label, w in (("N", n), ("N+1", plus), ("N-1", minus)):
        if w.other_densities:
            notes.append(f"The {label} file also contains {', '.join(w.other_densities)}; "
                         f"the {gw.DENSITY_SECTION} is used here.")
    return p, q, notes


def frontier_degeneracy(w: gw.FchkWavefunction) -> Optional[Tuple[int, int]]:
    """(LUMO degeneracy p, HOMO degeneracy q) of a closed-shell state, Multiwfn's 0.01 eV rule."""
    e = w.alpha_orbital_energies
    if w.multiplicity != 1 or e is None or not 0 < w.n_alpha < len(e):
        return None
    e = e * HARTREE_TO_EV
    homo, lumo = e[w.n_alpha - 1], e[w.n_alpha]
    q = int(np.sum(np.abs(e[:w.n_alpha] - homo) < DEGENERACY_EV))
    p = int(np.sum(np.abs(e[w.n_alpha:] - lumo) < DEGENERACY_EV))
    return p, q


def has_diffuse_functions(basis_set: str) -> bool:
    b = basis_set.upper()
    return "+" in b or b.startswith("AUG-") or "D-AUG" in b


# ----------------------------------------------------------------------
# Formulas
# ----------------------------------------------------------------------
def fukui_from_densities(rho_n, rho_plus, rho_minus, p: int = 1, q: int = 1) -> Dict[str, np.ndarray]:
    f_plus = (np.asarray(rho_plus) - rho_n) / p
    f_minus = (np.asarray(rho_n) - rho_minus) / q
    return {"f+": f_plus, "f-": f_minus, "f0": 0.5 * (f_plus + f_minus), "DD": f_plus - f_minus}


@dataclass
class GlobalIndices:
    """Finite-difference global CDFT indices, Hartree, Multiwfn conventions (manual 3.25.1)."""
    e_n: float
    e_plus: float
    e_minus: float
    homo_n: Optional[float] = None

    @property
    def ip(self) -> float:
        return self.e_minus - self.e_n

    @property
    def ea(self) -> float:
        return self.e_n - self.e_plus

    @property
    def chi(self) -> float:
        return 0.5 * (self.ip + self.ea)

    @property
    def mu(self) -> float:
        return -self.chi

    @property
    def eta(self) -> float:
        return self.ip - self.ea  # no 1/2 (fundamental gap), as in Multiwfn

    @property
    def softness(self) -> float:
        return 1.0 / self.eta if self.eta > 0 else float("nan")

    @property
    def omega(self) -> float:
        return self.mu ** 2 / (2.0 * self.eta) if self.eta > 0 else float("nan")

    @property
    def nucleophilicity(self) -> Optional[float]:
        return None if self.homo_n is None else self.homo_n - TCE_HOMO_HARTREE


def homo_energy(w: gw.FchkWavefunction) -> Optional[float]:
    e = w.alpha_orbital_energies
    if w.multiplicity != 1 or e is None or not 0 < w.n_alpha <= len(e):
        return None
    return float(e[w.n_alpha - 1])


# ----------------------------------------------------------------------
# Condensed indices: Multiwfn CDFT.txt
# ----------------------------------------------------------------------
@dataclass
class CondensedRow:
    index: int        # 1-based
    symbol: str
    q_n: float
    q_plus: float
    q_minus: float
    f_minus: float
    f_plus: float
    f_zero: float
    cdd: float


CDFT_HEADER = ["Atom", "q(N)", "q(N+1)", "q(N-1)", "f-", "f+", "f0", "CDD"]
_CDFT_ROW = re.compile(r"^\s*(\d+)\(\s*([A-Za-z]{1,3})\s*\)\s+(.*)$")
_CDFT_GLOBAL = re.compile(r"^\s*([^:]+?)\s*:\s+(-?\d+\.\d+)\s+Hartree")
CDFT_GLOBAL_KEYS = {
    "E(N)": "e_n", "E(N+1)": "e_plus", "E(N-1)": "e_minus", "E_HOMO(N)": "homo_n",
    "First vertical IP": "ip", "Vertical IP": "ip",
    "First vertical EA": "ea", "Vertical EA": "ea",
    "Mulliken electronegativity": "chi", "Chemical potential": "mu",
    "Hardness (=fundamental gap)": "eta", "Softness": "softness",
    "Electrophilicity index": "omega", "Nucleophilicity index": "nucleophilicity",
}


def parse_cdft_text(text: str) -> Tuple[List[CondensedRow], Dict[str, float]]:
    """Condensed rows and global values (Hartree) from Multiwfn's CDFT.txt."""
    lines = text.splitlines()
    start = next((i for i, l in enumerate(lines) if l.split()[:8] == CDFT_HEADER), None)
    if start is None:
        raise ValueError("No 'Atom q(N) q(N+1) q(N-1) f- f+ f0 CDD' table found: this is not a "
                         "CDFT.txt written by Multiwfn main function 22, option 2.")
    rows: List[CondensedRow] = []
    for line in lines[start + 1:]:
        m = _CDFT_ROW.match(line)
        if not m:
            break
        try:
            values = [float(x) for x in m.group(3).split()]
        except ValueError:
            raise ValueError(f"CDFT.txt: unreadable condensed row '{line.strip()}'.") from None
        if len(values) != 7:
            raise ValueError(f"CDFT.txt: row '{line.strip()}' has {len(values)} numbers, expected 7.")
        rows.append(CondensedRow(int(m.group(1)), m.group(2), *values))
    if not rows:
        raise ValueError("CDFT.txt: the condensed table is empty.")
    found: Dict[str, float] = {}
    for line in lines:
        m = _CDFT_GLOBAL.match(line)
        if m and m.group(1) in CDFT_GLOBAL_KEYS:
            found[CDFT_GLOBAL_KEYS[m.group(1)]] = float(m.group(2))
    return rows, found


def condensed_inconsistencies(rows: Sequence[CondensedRow], tol: float = CONDENSED_TOLERANCE) -> List[int]:
    """1-based atoms whose printed f-, f+, f0, CDD disagree with their own charges."""
    bad = []
    for r in rows:
        checks = (r.q_minus - r.q_n - r.f_minus, r.q_n - r.q_plus - r.f_plus,
                  0.5 * (r.q_minus - r.q_plus) - r.f_zero, r.f_plus - r.f_minus - r.cdd)
        if max(abs(c) for c in checks) > tol:
            bad.append(r.index)
    return bad


def check_rows_match(rows: Sequence[CondensedRow], w: gw.FchkWavefunction) -> None:
    symbols = w.atom_symbols()
    if len(rows) != len(symbols):
        raise ValueError(f"CDFT.txt lists {len(rows)} atoms but the N file has {len(symbols)}.")
    for r, s in zip(rows, symbols):
        if r.symbol.lower() != s.lower():
            raise ValueError(f"CDFT.txt atom {r.index} is {r.symbol} but atom {r.index} of the N file is {s}.")


# ----------------------------------------------------------------------
# Condensed indices: running Multiwfn
# ----------------------------------------------------------------------
def find_multiwfn(explicit: Optional[str] = None) -> Optional[str]:
    candidates = []
    if explicit:
        candidates.append(explicit)
    env = os.environ.get("Multiwfnpath")
    if env:
        candidates.append(os.path.join(env, "Multiwfn.exe"))
    for name in ("Multiwfn.exe", "Multiwfn"):
        found = shutil.which(name)
        if found:
            candidates.append(found)
    for c in candidates:
        if c and os.path.isfile(c):
            return c
    return None


def multiwfn_cdft_script(settings_found: bool) -> str:
    """Keystrokes for main function 22, option 2, with N.fch given on the command line.

    N.wfn is absent, so Multiwfn asks for the N, N+1 and N-1 files in turn
    (verified with Multiwfn 2026.7.15). Without a settings.ini Multiwfn first
    pauses for ENTER, so only then does the script start with a blank line.
    """
    lead = "" if settings_found else "\n"
    return lead + "22\n2\nN.fch\nN+1.fch\nN-1.fch\n0\nq\n"


def run_multiwfn_cdft(multiwfn: str, paths: Sequence[str], timeout: float = MULTIWFN_TIMEOUT_S) -> str:
    """CDFT.txt text from Multiwfn; RuntimeError when Multiwfn cannot produce it."""
    settings = Path(multiwfn).resolve().parent / "settings.ini"
    with tempfile.TemporaryDirectory() as tmp:
        for name, src in zip(("N.fch", "N+1.fch", "N-1.fch"), paths):
            shutil.copyfile(src, Path(tmp) / name)
        if settings.is_file():  # Multiwfn reads settings.ini from the current folder first
            shutil.copyfile(settings, Path(tmp) / "settings.ini")
        try:
            proc = subprocess.run([multiwfn, "N.fch"], input=multiwfn_cdft_script(settings.is_file()),
                                  text=True, cwd=tmp, capture_output=True, timeout=timeout,
                                  creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        except (OSError, subprocess.SubprocessError) as exc:
            raise RuntimeError(f"Multiwfn could not be run: {exc}") from exc
        out = Path(tmp) / "CDFT.txt"
        if not out.is_file():
            tail = " | ".join((proc.stdout or "").strip().splitlines()[-3:])
            raise RuntimeError(f"Multiwfn finished without writing CDFT.txt (exit {proc.returncode}): {tail}")
        return out.read_text(encoding="utf-8", errors="replace")
```

- [ ] **Step 4: Run the tests**

Run: `..\Scripts\python.exe -m pytest tests/test_fukuiFunctionAnalyze.py -q`
Expected: all pass, 33 tests. Two are guaranteed by verified data: `test_global_indices_reproduce_multiwfn` (every printed digit of `phenol_CDFT.txt`) and `test_phenol_rows_are_self_consistent`.

- [ ] **Step 5: Commit**

```bash
git add fukuiFunctionAnalyze.py tests/test_fukuiFunctionAnalyze.py
git commit -m "Fukui core: finite-difference formulas, N/N+1/N-1 validation, CDFT.txt parser, scripted Multiwfn engine

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: `FukuiAnalyser`: grids, surface, extrema, isosurface lobes, condensed sources, report and cubes

**Files:**
- Modify: `fukuiFunctionAnalyze.py` (append the Results/Analyser section)
- Test: `tests/test_fukuiFunctionAnalyze.py` (append)

**Interfaces:**
- Consumes: Task 5's functions. `molecularSurface.grid_box`, `grid_points`, `build_isosurface`, `density_outward` and `find_surface_extrema`. `cubeFile.CubeData`/`write_cube`.
- Produces:
  - `LobeSurface(vertices, triangles, normals, values)`;
  - `FukuiResults(paths, states, p, q, fukui_isovalue, surface_isovalue, spacing, origin, grids, surface, surface_values, extrema, lobes, global_indices, condensed, condensed_source, notes, elapsed)`, with `.n_state` and `.integral(key)`;
  - `FukuiAnalyser.run(n_path, plus_path, minus_path, fukui_isovalue, surface_isovalue, spacing, condensed_source, multiwfn_path, cdft_path, progress)`;
  - `FukuiAnalyser.lobe_surface(...)`, `.nearest_atom_label(w, pos)`, `.extrema_rows(results, key)`, `.site_ranking(rows, top)`, `.report_text(results, written)`, `.write_outputs(results, output_dir, write_cubes=True) -> [report, cubes...]`.

**Design notes** (why the code looks like this):
- **Grid size.** One grid serves both the lobes and the mapping surface. `grid_box` is driven by a single envelope, max(ρ_N/ρ_iso, |f⁺|/f_iso, |f⁻|/f_iso), with an isovalue of 1. So every face ends with each quantity below a tenth of its isovalue. |f⁰| and |Δf| are bounded by these two, so no lobe is truncated.
- **Lobes.** The +iso and −iso surfaces are built from the grid of `sign·f`. `build_isosurface` refines each cut point by bisecting the *exact* `sign·f`, which `fukui_at` evaluates from the three densities.
- **Condensed sources.** Missing or failing Multiwfn degrades to a note. A malformed `CDFT.txt` is an error, because the user supplied that file.

- [ ] **Step 1: Write the failing test.** Append to `tests/test_fukuiFunctionAnalyze.py`:

```python
@pytest.fixture(scope="module")
def phenol():
    seen = []
    results = fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), spacing=0.3,
                                     condensed_source=fk.CONDENSED_IMPORT, cdft_path=str(CDFT),
                                     progress=lambda f, m: seen.append(f))
    return results, seen


class TestPhenol:
    def test_progress_is_monotone_to_one(self, phenol):
        _, seen = phenol
        assert seen[-1] == 1.0 and all(b >= a - 1e-12 for a, b in zip(seen, seen[1:]))

    def test_grid_matches_multiwfn(self):
        data = np.loadtxt(REFERENCE_DIR / "phenol_grid_sample.txt", comments="#")
        states = [gw.read_fchk(p) for p in (N, PLUS, MINUS)]
        mine = fk.fukui_from_densities(*[gw.electron_density(w, data[:, :3]) for w in states])
        for k, key in enumerate(fk.FUNCTIONS):
            assert np.max(np.abs(mine[key] - data[:, 3 + k])) < 2e-5

    def test_integrals_are_close_to_one_and_zero(self, phenol):
        r, _ = phenol
        for key in ("f+", "f-", "f0"):
            assert 0.9 < r.integral(key) < 1.02
        assert abs(r.integral("DD")) < 0.05

    def test_globals_and_condensed_rows(self, phenol):
        r, _ = phenol
        assert r.global_indices.eta == pytest.approx(0.382666, abs=2e-6)
        assert len(r.condensed) == 13 and r.condensed_source == fk.CONDENSED_IMPORT
        assert r.notes == []

    def test_ranking_matches_the_manual_for_phenol(self, phenol):
        # Multiwfn manual 4.7.3: para C1 most favoured for electrophilic attack among ring carbons
        r, _ = phenol
        ranking = dict(fk.FukuiAnalyser.site_ranking(r.condensed, top=13))
        ring = [row.index for row in ranking["electrophilic attack: most negative CDD"] if row.symbol == "C"]
        assert ring[0] == 1

    def test_surface_maxima_of_f_minus_lie_on_ring_and_oxygen(self, phenol):
        r, _ = phenol
        maxima = [e for e in r.extrema["f-"] if e.kind == "maximum"]
        assert maxima and maxima[0].scope == "global"
        assert fk.FukuiAnalyser.nearest_atom_label(r.n_state, maxima[0].position)[0] in "CO"

    def test_lobes_carry_both_signs_for_the_dual_descriptor(self, phenol):
        r, _ = phenol
        lobe = r.lobes["DD"]
        assert lobe is not None and set(np.unique(lobe.values)) == {-r.fukui_isovalue, r.fukui_isovalue}
        assert lobe.triangles.max() < len(lobe.vertices)

    def test_outputs(self, phenol, tmp_path):
        r, _ = phenol
        written = fk.FukuiAnalyser().write_outputs(r, str(tmp_path))
        names = sorted(p.split("\\")[-1].split("/")[-1] for p in written)
        assert names == ["phenol_N_Fukui.txt", "phenol_N_dual.cub", "phenol_N_fminus.cub",
                         "phenol_N_fplus.cub", "phenol_N_fzero.cub"]
        text = (tmp_path / "phenol_N_Fukui.txt").read_text()
        assert "not predictions of reaction outcome" in text and "SITE RANKING" in text
        assert "Hardness" in text and "SURFACE EXTREMA OF f-" in text

    def test_cube_round_trip(self, phenol, tmp_path):
        import cubeFile as cf
        r, _ = phenol
        fk.FukuiAnalyser().write_outputs(r, str(tmp_path))
        cube = cf.read_cube(tmp_path / "phenol_N_fminus.cub")
        assert np.allclose(cube.values, r.grids["f-"], rtol=1e-5, atol=1e-10)


class TestCondensedSources:
    def test_missing_multiwfn_continues_with_a_note(self, monkeypatch):
        monkeypatch.setattr(fk, "find_multiwfn", lambda explicit=None: None)
        r = fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), spacing=0.4)
        assert r.condensed == [] and "not found" in r.notes[-1]

    def test_failing_multiwfn_continues_with_a_note(self, monkeypatch):
        monkeypatch.setattr(fk, "find_multiwfn", lambda explicit=None: "C:/fake/Multiwfn.exe")

        def broken(*args, **kwargs):
            raise RuntimeError("Multiwfn finished without writing CDFT.txt (exit 24): test")
        monkeypatch.setattr(fk, "run_multiwfn_cdft", broken)
        r = fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), spacing=0.4)
        assert r.condensed == [] and "without writing CDFT.txt" in r.notes[-1]

    @pytest.mark.skipif(fk.find_multiwfn() is None, reason="Multiwfn not installed")
    def test_multiwfn_engine_reproduces_the_reference(self):
        text = fk.run_multiwfn_cdft(fk.find_multiwfn(), [str(N), str(PLUS), str(MINUS)])
        assert fk.parse_cdft_text(text)[0] == fk.parse_cdft_text(CDFT_TEXT)[0]

    def test_unknown_source(self):
        with pytest.raises(ValueError):
            fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), condensed_source="magic")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `..\Scripts\python.exe -m pytest tests/test_fukuiFunctionAnalyze.py -q -k "Phenol or CondensedSources"`
Expected: FAIL with `AttributeError: module 'fukuiFunctionAnalyze' has no attribute 'FukuiAnalyser'`.

- [ ] **Step 3: Implement.** Append to `fukuiFunctionAnalyze.py`:

```python
# ----------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------
@dataclass
class LobeSurface:
    vertices: np.ndarray   # (V, 3) bohr
    triangles: np.ndarray  # (F, 3)
    normals: np.ndarray    # (V, 3)
    values: np.ndarray     # (V,) +isovalue or -isovalue


@dataclass
class FukuiResults:
    paths: Tuple[str, str, str]
    states: Tuple[gw.FchkWavefunction, gw.FchkWavefunction, gw.FchkWavefunction]
    p: int
    q: int
    fukui_isovalue: float
    surface_isovalue: float
    spacing: float
    origin: np.ndarray
    grids: Dict[str, np.ndarray]                 # a.u., on the cube grid
    surface: ms.SurfaceMesh                      # rho(N) = surface_isovalue, bohr
    surface_values: Dict[str, np.ndarray]        # a.u., one per vertex
    extrema: Dict[str, List[ms.SurfaceExtremum]]
    lobes: Dict[str, Optional[LobeSurface]]
    global_indices: Optional[GlobalIndices]
    condensed: List[CondensedRow]
    condensed_source: str
    notes: List[str] = field(default_factory=list)
    elapsed: float = 0.0

    @property
    def n_state(self) -> gw.FchkWavefunction:
        return self.states[0]

    def integral(self, key: str) -> float:
        return float(self.grids[key].sum() * self.spacing ** 3)


class FukuiAnalyser:
    """Headless pipeline: three fchk -> validation -> grids -> surface -> condensed indices."""

    def run(self, n_path: str, plus_path: str, minus_path: str,
            fukui_isovalue: float = DEFAULT_FUKUI_ISOVALUE,
            surface_isovalue: float = DEFAULT_SURFACE_ISOVALUE,
            spacing: float = DEFAULT_SPACING_BOHR,
            condensed_source: str = CONDENSED_MULTIWFN,
            multiwfn_path: Optional[str] = None, cdft_path: Optional[str] = None,
            progress: Optional[ProgressCallback] = None) -> FukuiResults:
        started = time.perf_counter()

        def report(fraction: float, message: str):
            if progress is not None:
                try:
                    progress(fraction, message)
                except Exception:
                    pass  # reporting must never abort an analysis

        if condensed_source not in CONDENSED_SOURCES:
            raise ValueError(f"Unknown condensed-index source '{condensed_source}'.")
        paths = (str(n_path), str(plus_path), str(minus_path))
        report(0.02, "Reading the three formatted checkpoints...")
        states = tuple(gw.read_fchk(p) for p in paths)
        n, plus, minus = states
        p, q, notes = validate_states(n, plus, minus)

        globals_ = None
        if p == 1 and q == 1:
            if None in (n.total_energy, plus.total_energy, minus.total_energy):
                notes.append("A 'Total Energy' entry is missing; global indices were not computed.")
            else:
                globals_ = GlobalIndices(n.total_energy, plus.total_energy, minus.total_energy, homo_energy(n))
                if not globals_.eta > 0:
                    notes.append("I - A is not positive; check the charges and multiplicities of the states.")
        else:
            notes.append(f"p = {p}, q = {q}: degenerate-frontier-orbital Fukui functions "
                         "(Multiwfn manual 3.25.4). Global indices need the N+1 and N-1 energies "
                         "and were not computed.")
        degeneracy = frontier_degeneracy(n)
        if degeneracy is not None and (degeneracy[0], degeneracy[1]) != (p, q) and max(degeneracy) > 1:
            notes.append(f"The N state has LUMO degeneracy {degeneracy[0]} and HOMO degeneracy "
                         f"{degeneracy[1]} (orbitals within {DEGENERACY_EV} eV) but the files give "
                         f"p = {p}, q = {q}. Standard f+/f- can then break the molecular symmetry; "
                         "use N+p / N-q states (manual 3.25.4) or orbital-weighted Fukui functions "
                         "(Multiwfn main function 22, options 4-7).")
        if has_diffuse_functions(n.basis_set):
            notes.append(f"The basis set {n.basis_set} has diffuse functions. They improve the electron "
                         "affinity, but with self-interaction-prone functionals the N+1 density can become "
                         "over-diffuse and f+ / the dual descriptor lose chemical meaning (manual 4.22.3.1).")

        def densities(points):
            return [gw.electron_density(w, points) for w in states]

        def fukui_at(points):
            return fukui_from_densities(*densities(points), p, q)

        def envelope(points):
            rn, rp, rm = densities(points)
            f = fukui_from_densities(rn, rp, rm, p, q)
            return np.maximum.reduce([rn / surface_isovalue,
                                      np.abs(f["f+"]) / fukui_isovalue,
                                      np.abs(f["f-"]) / fukui_isovalue])

        report(0.08, "Sizing the grid...")
        radii = np.array([AtomicData.vdw_radii.get(s, AtomicData.DEFAULT_VDW_RADIUS)
                          for s in n.atom_symbols()]) / BOHR_RADIUS_ANGSTROM
        origin, shape = ms.grid_box(n.coordinates, radii, envelope, 1.0, spacing)
        points = ms.grid_points(origin, spacing, shape)
        rho = []
        for k, (label, w) in enumerate(zip(("N", "N+1", "N-1"), states)):
            report(0.15 + 0.15 * k, f"Electron density of the {label} state on {len(points)} points...")
            rho.append(gw.electron_density(w, points).reshape(shape))
        grids = fukui_from_densities(*rho, p, q)

        report(0.62, f"rho(N) = {surface_isovalue:g} a.u. surface and the functions on it...")
        density_n = lambda pts: gw.electron_density(n, pts)  # noqa: E731
        surface = ms.build_isosurface(rho[0], origin, spacing, surface_isovalue,
                                      density_n, ms.density_outward(density_n))
        surface_values = fukui_at(surface.vertices)
        extrema = {k: ms.find_surface_extrema(v, surface.vertices, surface.neighbours)
                   for k, v in surface_values.items()}

        report(0.75, f"Isosurfaces at +/-{fukui_isovalue:g} a.u....")
        lobes = {k: self.lobe_surface(grids[k], origin, spacing, fukui_isovalue,
                                      lambda pts, key=k: fukui_at(pts)[key]) for k in FUNCTIONS}
        for k, lobe in lobes.items():
            if lobe is None:
                notes.append(f"{k} never reaches +/-{fukui_isovalue:g} a.u.; lower the Fukui isovalue to see it.")

        report(0.85, "Condensed indices...")
        condensed, source = self._condensed(condensed_source, paths, n, p, q, multiwfn_path, cdft_path, notes)
        report(1.0, "Done.")
        return FukuiResults(paths, states, p, q, fukui_isovalue, surface_isovalue, spacing, origin,
                            grids, surface, surface_values, extrema, lobes, globals_, condensed, source,
                            notes, time.perf_counter() - started)

    @staticmethod
    def lobe_surface(grid, origin, spacing, isovalue, exact) -> Optional[LobeSurface]:
        """The +isovalue and -isovalue surfaces of one function, as one mesh."""
        verts, tris, norms, vals = [], [], [], []
        offset = 0
        for sign in (1.0, -1.0):
            if (sign * grid).max() < isovalue:
                continue
            signed = lambda pts, s=sign: s * exact(pts)  # noqa: E731
            mesh = ms.build_isosurface(sign * grid, origin, spacing, isovalue, signed, ms.density_outward(signed))
            verts.append(mesh.vertices)
            tris.append(mesh.triangles + offset)
            norms.append(mesh.vertex_normals())
            vals.append(np.full(len(mesh.vertices), sign * isovalue))
            offset += len(mesh.vertices)
        if not verts:
            return None
        return LobeSurface(np.vstack(verts), np.vstack(tris), np.vstack(norms), np.concatenate(vals))

    @staticmethod
    def _condensed(source, paths, n, p, q, multiwfn_path, cdft_path, notes):
        if source == CONDENSED_NONE:
            return [], source
        if (p, q) != (1, 1):
            notes.append("Condensed indices for N+p / N-q states are not supported here; use Multiwfn "
                         "main function 22 option -3 then option 2 (manual 3.25.4).")
            return [], CONDENSED_NONE
        try:
            if source == CONDENSED_IMPORT:
                if not cdft_path:
                    raise RuntimeError("no CDFT.txt file was given")
                text = Path(cdft_path).read_text(encoding="utf-8", errors="replace")
            else:
                exe = find_multiwfn(multiwfn_path)
                if exe is None:
                    raise RuntimeError("Multiwfn.exe was not found (set its path, or Multiwfnpath, or PATH)")
                text = run_multiwfn_cdft(exe, paths)
        except (OSError, RuntimeError) as exc:
            notes.append(f"Condensed indices were not computed: {exc}.")
            return [], CONDENSED_NONE
        rows, _ = parse_cdft_text(text)
        check_rows_match(rows, n)
        bad = condensed_inconsistencies(rows)
        if bad:
            notes.append(f"CDFT.txt rows for atoms {bad} do not follow f- = q(N-1) - q(N), "
                         "f+ = q(N) - q(N+1); check the file.")
        negative = [r.index for r in rows if r.f_plus < 0 or r.f_minus < 0]
        if negative:
            notes.append(f"Negative condensed f+ or f- on atoms {negative}. Hirshfeld values are usually "
                         "non-negative; negative values can come from orbital relaxation or the partition "
                         "and should not be read as reactivity.")
        return rows, source

    # --- reporting ------------------------------------------------------
    @staticmethod
    def nearest_atom_label(w: gw.FchkWavefunction, position_bohr: np.ndarray) -> str:
        i = int(np.argmin(np.linalg.norm(w.coordinates - position_bohr, axis=1)))
        return f"{w.atom_symbols()[i]}{i + 1}"

    def extrema_rows(self, results: FukuiResults, key: str) -> List[Dict[str, object]]:
        rows = []
        for e in results.extrema[key]:
            rows.append({"row": len(rows) + 1, "kind": e.kind, "scope": e.scope, "value": e.value,
                         "position_angstrom": e.position * BOHR_RADIUS_ANGSTROM,
                         "nearest_atom": self.nearest_atom_label(results.n_state, e.position)})
        return rows

    @staticmethod
    def site_ranking(rows: Sequence[CondensedRow], top: int = 5) -> List[Tuple[str, List[CondensedRow]]]:
        def best(key, reverse):
            return sorted(rows, key=key, reverse=reverse)[:top]
        return [
            ("electrophilic attack: largest f-", best(lambda r: r.f_minus, True)),
            ("nucleophilic attack: largest f+", best(lambda r: r.f_plus, True)),
            ("radical attack: largest f0", best(lambda r: r.f_zero, True)),
            ("electrophilic attack: most negative CDD", best(lambda r: r.cdd, False)),
            ("nucleophilic attack: most positive CDD", best(lambda r: r.cdd, True)),
        ]

    def report_text(self, results: FukuiResults, written: Sequence[str] = ()) -> str:
        n, plus, minus = results.states
        p, q = results.p, results.q
        lines = [
            "# gQTEA Fukui Function Analysis",
            "# Indicators of likely reactive regions -- not predictions of reaction outcome.",
            f"# file_N                 {Path(results.paths[0]).name}",
            f"# file_N+{p}               {Path(results.paths[1]).name}",
            f"# file_N-{q}               {Path(results.paths[2]).name}",
            f"# method / basis         {n.method} / {n.basis_set}",
            f"# electrons              {n.n_electrons} / {plus.n_electrons} / {minus.n_electrons}",
            f"# charge multiplicity    N: {n.charge} {n.multiplicity}   N+{p}: {plus.charge} {plus.multiplicity}"
            f"   N-{q}: {minus.charge} {minus.multiplicity}",
            f"# density                {gw.DENSITY_SECTION}",
            f"# definitions            f+ = [rho(N+{p}) - rho(N)]/{p}; f- = [rho(N) - rho(N-{q})]/{q}; "
            "f0 = (f+ + f-)/2; DD = f+ - f-",
            f"# grid_spacing_bohr      {results.spacing:g}",
            f"# fukui_isovalue_au      {results.fukui_isovalue:g}",
            f"# surface                rho(N) = {results.surface_isovalue:g} a.u., "
            f"{len(results.surface.vertices)} vertices",
            "",
        ]
        g = results.global_indices
        if g is not None:
            def both(value):
                return f"{value:12.6f} Hartree  {value * HARTREE_TO_EV:10.4f} eV"
            lines += [
                "GLOBAL INDICES (hardness = I - A, without 1/2, as in Multiwfn)",
                f"  E(N)                 {g.e_n:.6f} Hartree",
                f"  E(N+1)               {g.e_plus:.6f} Hartree",
                f"  E(N-1)               {g.e_minus:.6f} Hartree",
                f"  Vertical IP          {both(g.ip)}",
                f"  Vertical EA          {both(g.ea)}",
                f"  Electronegativity    {both(g.chi)}",
                f"  Chemical potential   {both(g.mu)}",
                f"  Hardness             {both(g.eta)}",
                f"  Softness             {g.softness:12.6f} Hartree^-1",
                f"  Electrophilicity     {both(g.omega)}",
            ]
            if g.nucleophilicity is not None:
                lines.append(f"  Nucleophilicity      {both(g.nucleophilicity)}   "
                             "(E_HOMO(TCE) at B3LYP/6-31G*, as in Multiwfn)")
            lines.append("")
        lines.append("GRID INTEGRALS (ideally 1 for f+, f-, f0 and 0 for DD; a finite uniform grid undercounts)")
        lines.append("  " + "   ".join(f"{k} {results.integral(k):.4f}" for k in FUNCTIONS))
        lines.append("")
        if results.condensed:
            s = g.softness if g is not None else float("nan")
            lines += [f"CONDENSED INDICES  (source: {results.condensed_source}; Hirshfeld charges, e; "
                      "s = S x f in e/Hartree)",
                      "  atom      q(N)   q(N+1)   q(N-1)       f-       f+       f0      CDD       s-       s+       s0"]
            for r in results.condensed:
                lines.append(f"  {r.symbol:>2s}{r.index:<4d} {r.q_n:8.4f} {r.q_plus:8.4f} {r.q_minus:8.4f} "
                             f"{r.f_minus:8.4f} {r.f_plus:8.4f} {r.f_zero:8.4f} {r.cdd:8.4f} "
                             f"{s * r.f_minus:8.4f} {s * r.f_plus:8.4f} {s * r.f_zero:8.4f}")
            lines += ["", "SITE RANKING (consider only atoms that can actually react)"]
            for title, best in self.site_ranking(results.condensed):
                lines.append(f"  {title:42s} " + ", ".join(f"{r.symbol}{r.index}" for r in best))
            lines.append("")
        for key in FUNCTIONS:
            lines += [f"SURFACE EXTREMA OF {FUNCTION_TITLES[key]}  on rho(N) = {results.surface_isovalue:g} a.u.",
                      "   row  type     scope        value (a.u.)    x (A)     y (A)     z (A)   atom"]
            for r in self.extrema_rows(results, key):
                x, y, z = r["position_angstrom"]
                lines.append(f"  {r['row']:4d}  {r['kind']:7s}  {r['scope']:6s}  {r['value']:16.8f}  "
                             f"{x:8.4f}  {y:8.4f}  {z:8.4f}   {r['nearest_atom']}")
            lines.append("")
        for note in results.notes:
            lines.append(f"NOTE: {note}")
        for path in written:
            lines.append(f"Written: {path}")
        return "\n".join(lines) + "\n"

    def write_outputs(self, results: FukuiResults, output_dir: str, write_cubes: bool = True) -> List[str]:
        stem = Path(results.paths[0]).stem
        out = Path(output_dir)
        written: List[str] = []
        n = results.n_state
        if write_cubes:
            for key in FUNCTIONS:
                cube = cf.CubeData(f"{stem} {FUNCTION_TITLES[key]}", "gQTEA Fukui Function Analysis, a.u.",
                                   n.atomic_numbers, n.nuclear_charges, n.coordinates,
                                   results.origin, results.spacing, results.grids[key])
                written.append(cf.write_cube(out / (stem + CUBE_SUFFIXES[key]), cube))
        report = out / (stem + REPORT_SUFFIX)
        try:
            report.parent.mkdir(parents=True, exist_ok=True)
            report.write_text(self.report_text(results, written), encoding="utf-8")
        except OSError as exc:
            raise ValueError(f"Could not write {report}: {exc}") from exc
        return [str(report)] + written
```

- [ ] **Step 4: Run the tests**

Run: `..\Scripts\python.exe -m pytest tests/test_fukuiFunctionAnalyze.py tests/test_cubeFile.py -q`
Expected: all pass, in about 25 s.
- `test_grid_matches_multiwfn` compares against Multiwfn's own cubes at 300 nodes, within 2×10⁻⁵ a.u. The measured maximum was 8×10⁻⁶.
- `test_multiwfn_engine_reproduces_the_reference` runs only when `find_multiwfn()` finds an installation; otherwise it is skipped.

- [ ] **Step 5: Commit**

```bash
git add fukuiFunctionAnalyze.py tests/test_fukuiFunctionAnalyze.py
git commit -m "FukuiAnalyser: f+/f-/f0/DD grids, rho(N)=0.01 surface extrema, lobes, condensed sources, report and cubes

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: UI, help text and registration

**Files:**
- Modify: `fukuiFunctionAnalyze.py` (append the Toga UI section)
- Modify: `help.py` (add `help_fukui` to `HelpGqteaWin`, after `help_esp_surface`)
- Modify: `gqteaWinToga.py` (import + Tools entry)
- Test: `tests/test_fukuiFunctionAnalyze.py` (append the UI tests)

**Interfaces:**
- Consumes: `FukuiAnalyser`, `ESPSurfaceViewer` with `caption=` (Task 4), `molecularViewer.save_canvas_image`.
- Produces: `FukuiFunctionAnalyzeUI(*args, **kwargs)`, plus `DISPLAY_LOBES`, `DISPLAY_MAPPED`, `DISPLAY_MODES` and `BACKGROUND_COLOURS`. Behaviour follows the ESP tool:
  - `workflow` runs `run` in `asyncio.to_thread`, with progress through `loop.call_soon_threadsafe`.
  - It always writes the report; it writes the cubes only when the switch is on.
  - It switches to the Results tab when done.
  - The capture-then-ask image export contract is unchanged.
- **UI contract test.** Every `run()` parameter must be supplied by `_collect_inputs`, so adding an engine option without a control fails the suite.

- [ ] **Step 1: Write the failing test.** Append to `tests/test_fukuiFunctionAnalyze.py`:

```python
def make_ui():
    """The UI class without a window: stub widgets, as in test_electrostaticSurfacePotential."""
    ui = fk.FukuiFunctionAnalyzeUI.__new__(fk.FukuiFunctionAnalyzeUI)
    ui.analyser = fk.FukuiAnalyser()
    ui.results = None
    ui.written = []
    ui.viewer = None
    ui._syncing = False
    for name in ("n_input", "plus_input", "minus_input", "output_dir_input", "fukui_iso_input",
                 "surface_iso_input", "spacing_input", "multiwfn_input", "cdft_input"):
        setattr(ui, name, SimpleNamespace(value=""))
    ui.condensed_selection = SimpleNamespace(value=fk.CONDENSED_NONE)
    ui.cubes_switch = SimpleNamespace(value=False)
    ui.function_selection = SimpleNamespace(value="f-")
    ui.display_selection = SimpleNamespace(value=fk.DISPLAY_LOBES)
    ui.extrema_switch = SimpleNamespace(value=True)
    ui.background_selection = SimpleNamespace(value="Black")
    ui.summary_output = SimpleNamespace(value="")
    ui.status_label = SimpleNamespace(text="")
    ui.progress_bar = SimpleNamespace(value=0)
    ui.tabs = SimpleNamespace(current_tab=0)
    for name in ("run_button", "show_button", "clear_button"):
        setattr(ui, name, SimpleNamespace(enabled=True))
    ui.errors = []

    async def record(title, message):
        ui.errors.append((title, message))
    ui._show_error = record
    ui._show_info = record
    return ui


def fill(ui, tmp_path):
    ui.n_input.value, ui.plus_input.value, ui.minus_input.value = str(N), str(PLUS), str(MINUS)
    ui.output_dir_input.value = str(tmp_path)
    ui.spacing_input.value = "0.4"


class StubViewer:
    def __init__(self):
        self.surfaces = []

    def set_molecule(self, symbols, xyz):
        self.symbols = list(symbols)

    def set_display_options(self, **kwargs):
        self.options = kwargs

    def set_surface(self, *args, **kwargs):
        self.surfaces.append((args, kwargs))

    def clear_surface(self):
        self.surfaces.append(None)

    def viewer_window_is_open(self):
        return False


class TestUI:
    def test_missing_files_are_an_input_error(self):
        ui = make_ui()
        asyncio.run(ui.workflow(None))
        assert ui.errors and ui.errors[0][0] == "Input Error" and "N state" in ui.errors[0][1]

    def test_import_without_a_file_is_an_input_error(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        ui.condensed_selection.value = fk.CONDENSED_IMPORT
        with pytest.raises(ValueError, match="CDFT.txt"):
            ui._collect_inputs()

    def test_every_run_parameter_is_supplied_by_the_ui(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        wanted = set(inspect.signature(fk.FukuiAnalyser.run).parameters) - {"self", "progress"}
        assert wanted == set(ui._collect_inputs())

    def test_workflow_writes_the_report_and_switches_tab(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        assert not ui.errors and (tmp_path / "phenol_N_Fukui.txt").exists()
        assert not list(tmp_path.glob("*.cub"))  # cubes switch off
        assert ui.tabs.current_tab == 1 and ui.run_button.enabled

    def test_viewer_gets_lobes_or_the_mapped_surface(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        ui.viewer = StubViewer()
        ui._push_to_viewer()
        args, kwargs = ui.viewer.surfaces[-1]
        assert "blue positive, red negative" in kwargs["caption"] and kwargs["half_range"] == 0.005
        ui.display_selection.value = fk.DISPLAY_MAPPED
        ui._push_to_viewer()
        args, kwargs = ui.viewer.surfaces[-1]
        assert kwargs["caption"].startswith("f- (a.u.)") and args[5]  # extrema markers
        assert ui.default_image_name() == "phenol_N_fminus.png"

    def test_clear_resets_results(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        ui.clear_output(None)
        assert ui.results is None
        asyncio.run(ui.show_surface(None))
        assert ui.errors[-1][0] == "No Results"
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `..\Scripts\python.exe -m pytest tests/test_fukuiFunctionAnalyze.py -q -k TestUI`
Expected: FAIL with `AttributeError: module 'fukuiFunctionAnalyze' has no attribute 'FukuiFunctionAnalyzeUI'`.

- [ ] **Step 3: Add the help text.** In `help.py`, directly after the closing `"""` of `help_esp_surface`, add:

```python
    help_fukui = """FUKUI FUNCTION ANALYSIS
Fukui functions show where a molecule's electron density changes most when
it gains or loses an electron. They INDICATE likely reactive regions; they do
not predict reaction outcomes (sterics, solvent, kinetics are not included).

  f+(r) = rho(N+1) - rho(N)        large: susceptible to NUCLEOPHILIC attack
  f-(r) = rho(N) - rho(N-1)        large: susceptible to ELECTROPHILIC attack
  f0(r) = [rho(N+1) - rho(N-1)]/2  large: susceptible to RADICAL attack
  dual descriptor DD = f+ - f-     > 0 nucleophilic, < 0 electrophilic attack
  (compare sites by how positive/negative DD is, not only by its sign)

WHAT TO PREPARE (Gaussian)
  1. Optimize the N state (the molecule as it reacts) and run a frequency job.
  2. At that SAME geometry run three single points with the SAME method and
     basis set: N, N+1 and N-1. Do not re-optimize the ions.
     Neutral closed shell: N 0 1, N+1 -1 2, N-1 1 2.
     Neutral radical:      N 0 2, N+1 -1 1, N-1 1 1 (check other spins for
     transition-metal systems).
  3. formchk each .chk to .fchk and select the three files here.

CONDENSED (per-atom) INDICES
  f-_A = q_A(N-1) - q_A(N), f+_A = q_A(N) - q_A(N+1), f0_A = average,
  CDD_A = f+_A - f-_A, from Hirshfeld charges. They come from Multiwfn
  (main function 22, option 2): run here when Multiwfn.exe is found, or
  imported from an existing CDFT.txt.

OUTPUTS (next to the N file unless another folder is chosen)
  <name>_Fukui.txt  global indices, condensed table, site ranking, surface
                    extrema of f+, f-, f0, DD on rho(N) = 0.01 a.u., notes
  <name>_fplus.cub, _fminus.cub, _fzero.cub, _dual.cub  (a.u., bohr grid)
  Show 3D: isosurfaces (blue +, red -) or the function mapped on the
  rho(N) surface with its extrema (cyan minimum, orange maximum).

CAVEATS
  Diffuse basis functions with B3LYP-like functionals can make f+ and DD
  over-diffuse; degenerate HOMO/LUMO need N+p/N-q states or orbital-weighted
  Fukui functions. Read the NOTES at the end of the report.
"""
```

- [ ] **Step 4: Implement the UI.** Append to `fukuiFunctionAnalyze.py`:

```python
# ----------------------------------------------------------------------
# Toga UI
# ----------------------------------------------------------------------
DISPLAY_LOBES = "Isosurfaces (+/- isovalue)"
DISPLAY_MAPPED = "Mapped on rho(N) surface"
DISPLAY_MODES = (DISPLAY_LOBES, DISPLAY_MAPPED)
BACKGROUND_COLOURS = {
    "Black": (0.0, 0.0, 0.0),
    "White": (1.0, 1.0, 1.0),
    "Dark grey": (0.25, 0.25, 0.25),
    "Light grey": (0.85, 0.85, 0.85),
}
IMAGE_FILE_TYPES = ["png", "tif", "tiff", "bmp"]


class FukuiFunctionAnalyzeUI:
    LABEL_WIDTH = 230
    HINT_COLOR = "#666666"

    def __init__(self, *args, **kwargs):
        # gqteaWinToga passes the Button that was pressed; ignore it.
        self.analyser = FukuiAnalyser()
        self.results: Optional[FukuiResults] = None
        self.written: List[str] = []
        self.viewer: Optional[ESPSurfaceViewer] = None
        self._syncing = False
        self.window = toga.Window(title="Fukui Function Analysis", size=(900, 700), on_close=self._on_close)
        self.layout_main_window()
        self.window.show()

    # --- layout -----------------------------------------------------------
    def _label(self, text):
        return toga.Label(text, style=Pack(width=self.LABEL_WIDTH, margin=(6, 6, 0, 0)))

    def _row(self, *children):
        return toga.Box(children=list(children), style=Pack(direction=ROW, margin_bottom=4))

    def _heading(self, text):
        return toga.Label(text, style=Pack(font_weight="bold", margin=(10, 0, 4, 0)))

    def _hint(self, text):
        return toga.Label(text, style=Pack(color=self.HINT_COLOR, font_size=9, margin=(0, 0, 4, 0)))

    def _file_row(self, label, attr, handler):
        field_ = toga.TextInput(placeholder=".fchk / .fch", style=Pack(flex=1))
        setattr(self, attr, field_)
        return self._row(self._label(label), field_,
                         toga.Button("Browse", on_press=handler, style=Pack(width=90, margin_left=5)))

    def layout_main_window(self):
        self.output_dir_input = toga.TextInput(placeholder="Blank = folder of the N file", style=Pack(flex=1))
        header = toga.Box(children=[
            self._file_row("N state (optimized geometry):", "n_input", self.browse_n),
            self._file_row("N+1 state (anion of a neutral):", "plus_input", self.browse_plus),
            self._file_row("N-1 state (cation of a neutral):", "minus_input", self.browse_minus),
            self._row(self._label("Output folder:"), self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
        ], style=Pack(direction=COLUMN))

        self.fukui_iso_input = toga.TextInput(placeholder=f"{DEFAULT_FUKUI_ISOVALUE:g}", style=Pack(width=160))
        self.surface_iso_input = toga.TextInput(placeholder=f"{DEFAULT_SURFACE_ISOVALUE:g}", style=Pack(width=160))
        self.spacing_input = toga.TextInput(placeholder=f"{DEFAULT_SPACING_BOHR:g}", style=Pack(width=160))
        self.condensed_selection = toga.Selection(items=list(CONDENSED_SOURCES), style=Pack(width=260))
        self.multiwfn_input = toga.TextInput(value=find_multiwfn() or "", placeholder="Multiwfn.exe",
                                             style=Pack(flex=1))
        self.cdft_input = toga.TextInput(placeholder="CDFT.txt from Multiwfn main function 22, option 2",
                                         style=Pack(flex=1))
        self.cubes_switch = toga.Switch("Write cube files (f+, f-, f0, dual descriptor)", value=True)
        setup = toga.Box(children=[
            self._heading("Real-space functions"),
            self._row(self._label("Fukui isovalue (a.u.):"), self.fukui_iso_input),
            self._hint("Blank = 0.005. Isosurfaces are drawn at +value and -value."),
            self._row(self._label("rho(N) surface for mapping (a.u.):"), self.surface_iso_input),
            self._hint("Blank = 0.01, as in Multiwfn manual 4.12.4 (f is too small on rho = 0.001)."),
            self._row(self._label("Grid spacing (bohr):"), self.spacing_input),
            self._hint("Blank = 0.25 bohr. Smaller = smoother and slower."),
            self._heading("Atom-condensed indices (Hirshfeld, via Multiwfn)"),
            self._row(self._label("Source:"), self.condensed_selection),
            self._row(self._label("Multiwfn executable:"), self.multiwfn_input),
            self._row(self._label("Existing CDFT.txt:"), self.cdft_input,
                      toga.Button("Browse", on_press=self.browse_cdft, style=Pack(width=90, margin_left=5))),
            self._hint("If Multiwfn is missing or fails, the run continues without condensed indices."),
            self.cubes_switch,
        ], style=Pack(direction=COLUMN, margin=10))

        self.function_selection = toga.Selection(items=list(FUNCTIONS), on_change=self.on_display_change,
                                                 style=Pack(width=140))
        self.display_selection = toga.Selection(items=list(DISPLAY_MODES), on_change=self.on_display_change,
                                                style=Pack(width=240))
        self.extrema_switch = toga.Switch("Show extrema markers", value=True, on_change=self.on_display_change)
        self.background_selection = toga.Selection(items=list(BACKGROUND_COLOURS),
                                                   on_change=self.on_background_change, style=Pack(width=140))
        self.export_button = toga.Button("Export image...", on_press=self.export_image,
                                         style=Pack(width=130, margin_left=10))
        self.summary_output = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_fukui,
                                                      style=Pack(flex=1, font_family="monospace"))
        results = toga.Box(children=[
            self._row(self._label("Function:"), self.function_selection),
            self._row(self._label("Display:"), self.display_selection, self.extrema_switch),
            self._row(self._label("Canvas background:"), self.background_selection, self.export_button),
            self.summary_output,
        ], style=Pack(direction=COLUMN, margin=10, flex=1))

        self.tabs = toga.OptionContainer(content=[
            ("Setup", toga.ScrollContainer(content=setup, style=Pack(flex=1))),
            ("Results", results),
        ], style=Pack(flex=1))
        self.run_button = toga.Button("Run", on_press=self.workflow, style=Pack(width=110, margin=5))
        self.show_button = toga.Button("Show 3D", on_press=self.show_surface, style=Pack(width=110, margin=5))
        self.clear_button = toga.Button("Clear", on_press=self.clear_output, style=Pack(width=110, margin=5))
        self.close_button = toga.Button("Close", on_press=self.closeTopLevel, style=Pack(width=110, margin=5))
        self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
        self.progress_bar = toga.ProgressBar(max=100, style=Pack(margin=(0, 5, 5, 5)))
        self.window.content = toga.Box(children=[
            header,
            self._row(self.run_button, self.show_button, self.clear_button, self.close_button),
            self.tabs, self.status_label, self.progress_bar,
        ], style=Pack(direction=COLUMN, margin=10, flex=1))

    # --- dialogs ------------------------------------------------------------
    async def _show_error(self, title: str, message: str):
        await self.window.dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        await self.window.dialog(toga.InfoDialog(title, message))

    async def _browse_into(self, target, title, file_types):
        try:
            path = await self.window.dialog(toga.OpenFileDialog(title=title, file_types=file_types))
            if path:
                target.value = str(path)
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))

    async def browse_n(self, widget):
        await self._browse_into(self.n_input, "Formatted checkpoint of the N state", ["fchk", "fch"])

    async def browse_plus(self, widget):
        await self._browse_into(self.plus_input, "Formatted checkpoint of the N+1 state", ["fchk", "fch"])

    async def browse_minus(self, widget):
        await self._browse_into(self.minus_input, "Formatted checkpoint of the N-1 state", ["fchk", "fch"])

    async def browse_cdft(self, widget):
        await self._browse_into(self.cdft_input, "Multiwfn CDFT.txt", ["txt"])

    async def browse_output_dir(self, widget):
        try:
            folder = await self.window.dialog(toga.SelectFolderDialog(title="Choose the output folder"))
            if folder:
                self.output_dir_input.value = str(folder)
        except Exception as exc:
            await self._show_error("Folder Selection Error", str(exc))

    # --- run ----------------------------------------------------------------
    def _collect_inputs(self) -> Dict[str, object]:
        paths = {}
        for key, widget, label in (("n_path", self.n_input, "N"), ("plus_path", self.plus_input, "N+1"),
                                   ("minus_path", self.minus_input, "N-1")):
            path = (widget.value or "").strip()
            if not path:
                raise ValueError(f"Please select the {label} state .fchk file.")
            if not os.path.isfile(path):
                raise ValueError(f"File not found: {path}")
            paths[key] = path
        source = self.condensed_selection.value or CONDENSED_MULTIWFN
        cdft = (self.cdft_input.value or "").strip() or None
        if source == CONDENSED_IMPORT and (cdft is None or not os.path.isfile(cdft)):
            raise ValueError("Please select an existing CDFT.txt, or choose another condensed-index source.")
        return {
            **paths,
            "fukui_isovalue": parse_positive_float(self.fukui_iso_input.value, DEFAULT_FUKUI_ISOVALUE,
                                                   "Fukui isovalue", 1.0),
            "surface_isovalue": parse_positive_float(self.surface_iso_input.value, DEFAULT_SURFACE_ISOVALUE,
                                                     "rho surface isovalue", 1.0),
            "spacing": parse_positive_float(self.spacing_input.value, DEFAULT_SPACING_BOHR,
                                            "Grid spacing", 0.6),
            "condensed_source": source,
            "multiwfn_path": (self.multiwfn_input.value or "").strip() or None,
            "cdft_path": cdft,
        }

    def _apply_progress(self, fraction: float, message: str):
        self.progress_bar.value = max(0.0, min(100.0, 100.0 * fraction))
        self.status_label.text = message

    def _set_buttons(self, enabled: bool):
        for b in (self.run_button, self.show_button, self.clear_button):
            b.enabled = enabled

    async def workflow(self, widget):
        try:
            params = self._collect_inputs()
            output_dir = resolve_output_dir(params["n_path"], self.output_dir_input.value)
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        self._set_buttons(False)
        self.progress_bar.value = 0
        loop = asyncio.get_running_loop()

        def on_progress(fraction, message):
            loop.call_soon_threadsafe(self._apply_progress, fraction, message)

        try:
            results = await asyncio.to_thread(self.analyser.run, progress=on_progress, **params)
            self.written = self.analyser.write_outputs(results, output_dir, bool(self.cubes_switch.value))
            self.results = results
            self.summary_output.value = self.analyser.report_text(results, self.written)
            self.tabs.current_tab = 1
            self.status_label.text = f"Done in {results.elapsed:.1f} s. Report: {self.written[0]}"
            self.progress_bar.value = 100
            if self.viewer is not None and self.viewer.viewer_window_is_open():
                self._push_to_viewer()
        except Exception as exc:
            self.status_label.text = "Calculation failed."
            self.progress_bar.value = 0
            await self._show_error("Fukui Calculation Error", str(exc))
        finally:
            self._set_buttons(True)

    # --- display --------------------------------------------------------
    def selected_function(self) -> str:
        return self.function_selection.value if self.function_selection.value in FUNCTIONS else "f-"

    def surface_markers(self, key: str) -> List[SurfaceMarker]:
        return [SurfaceMarker(r["position_angstrom"], r["kind"], r["scope"], str(r["row"]))
                for r in self.analyser.extrema_rows(self.results, key)]

    def _push_to_viewer(self):
        if self.results is None or self.viewer is None:
            return
        key = self.selected_function()
        n = self.results.n_state
        self.viewer.set_molecule(n.atom_symbols(), n.coordinates * BOHR_RADIUS_ANGSTROM)
        self.viewer.set_display_options(unit=UNIT_AU, show_extrema=bool(self.extrema_switch.value))
        if self.display_selection.value == DISPLAY_MAPPED:
            s = self.results.surface
            self.viewer.set_surface(
                s.vertices * BOHR_RADIUS_ANGSTROM, s.triangles, s.vertex_normals(),
                self.results.surface_values[key], self.results.surface_isovalue, self.surface_markers(key),
                caption=f"{key} (a.u.)   mapped on rho(N) = {self.results.surface_isovalue:g} a.u.")
            return
        lobe = self.results.lobes[key]
        if lobe is None:
            self.viewer.clear_surface()
            self.status_label.text = f"{key} has no isosurface at +/-{self.results.fukui_isovalue:g} a.u."
            return
        iso = self.results.fukui_isovalue
        self.viewer.set_surface(lobe.vertices * BOHR_RADIUS_ANGSTROM, lobe.triangles, lobe.normals,
                                lobe.values, iso, (), half_range=iso,
                                caption=f"{key} isosurfaces at +/-{iso:g} a.u.: blue positive, red negative")

    def on_display_change(self, widget):
        if self._syncing or self.results is None or self.viewer is None:
            return
        self._push_to_viewer()

    async def show_surface(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Run the calculation first.")
            return
        if self.viewer is None:
            self.viewer = ESPSurfaceViewer()
            self.viewer.window_title = "Fukui Function Viewer"
            self.viewer.loop = asyncio.get_running_loop()
            self.viewer.set_background_color(self.selected_background())
        self._push_to_viewer()
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)

    def selected_background(self):
        return BACKGROUND_COLOURS.get(self.background_selection.value, BACKGROUND_COLOURS["Black"])

    def on_background_change(self, widget):
        if self.viewer is not None:
            self.viewer.set_background_color(self.selected_background())

    def default_image_name(self) -> str:
        if self.results is None:
            return "Fukui.png"
        tag = {"f+": "fplus", "f-": "fminus", "f0": "fzero", "DD": "dual"}[self.selected_function()]
        return f"{Path(self.results.paths[0]).stem}_{tag}.png"

    async def _ask_image_path(self, suggested: str):
        return await self.window.dialog(toga.SaveFileDialog(
            title="Export the Fukui function image", suggested_filename=suggested, file_types=IMAGE_FILE_TYPES))

    async def export_image(self, widget):
        """Capture first, then ask for the file name (same contract as the ESP tool)."""
        if self.viewer is None or not self.viewer.viewer_window_is_open():
            await self._show_info("No Image", "Open the 3D view with 'Show 3D' first.")
            return
        try:
            image = await self.viewer.capture_frame()
        except RuntimeError as exc:
            await self._show_error("Export Image Error", str(exc))
            return
        try:
            path = await self._ask_image_path(self.default_image_name())
            if not path:
                return
            saved = save_canvas_image(image, str(path))
        except (ValueError, OSError) as exc:
            await self._show_error("Export Image Error", str(exc))
            return
        self.status_label.text = f"Image saved: {saved}"

    def clear_output(self, widget):
        self.results = None
        self.written = []
        self.summary_output.value = HelpGqteaWin.help_fukui
        if self.viewer is not None:
            self.viewer.clear_surface()
        self.progress_bar.value = 0
        self.status_label.text = "Ready."

    def _on_close(self, window, **kwargs) -> bool:
        if self.viewer is not None:
            self.viewer.request_close()
        return True

    def closeTopLevel(self, widget):
        self._on_close(self.window)
        self.window.close()
```

- [ ] **Step 5: Register the tool.** In `gqteaWinToga.py`, after `from electrostaticSurfacePotential import ElectrostaticSurfacePotentialUI`, add

```python
from fukuiFunctionAnalyze import FukuiFunctionAnalyzeUI
```

and in `general_tools`, after `("Electrostatic Surface Potential (ESP)", ElectrostaticSurfacePotentialUI),`, add

```python
            ("Fukui Function Analysis", FukuiFunctionAnalyzeUI),
```

- [ ] **Step 6: Run the tests**

Run: `..\Scripts\python.exe -m pytest tests/test_fukuiFunctionAnalyze.py tests/test_cubeFile.py tests/test_espViewer.py tests/test_gaussianWavefunction.py tests/test_electrostaticSurfacePotential.py -q`
Expected: all pass. The Fukui and cube files account for 58 tests (59 with the caption test).

When appending in Tasks 6 and 7, separate each appended part from the existing text by **two blank lines**. The result must then equal the verified files. Check with the following, if `$SCRATCH` still exists; no output means identical:

```bash
SCRATCH="/c/Users/ajc/AppData/Local/Temp/claude/c--AJCamargo-gqteaWinToga-venv-src/dc50d304-21d9-487f-ba44-c3aa1d2a2200/scratchpad/impl"
for f in fukuiFunctionAnalyze.py cubeFile.py tests/test_fukuiFunctionAnalyze.py tests/test_cubeFile.py; do diff -q "$SCRATCH/$f" "$f"; done
```

- [ ] **Step 7: Commit**

```bash
git add fukuiFunctionAnalyze.py help.py gqteaWinToga.py tests/test_fukuiFunctionAnalyze.py
git commit -m "Fukui Function Analysis tool: Toga UI (Setup/Results), 3D lobes or rho-surface map, help, Tools entry

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Documentation, full suite, and a live smoke test

**Files:**
- Modify: `USER_MANUAL.md`: add a tool section before `#### Energy Plots`, and a formulas section before `## Formulas and Methods (Solvation box builders)`.
- Modify: `readme.md`: the capability list, the module list, and the accepted-files list.
- Modify: `CLAUDE.md`: a tool-specific note after the ESP bullet.

- [ ] **Step 1: USER_MANUAL tool section.** Insert before the line `#### Energy Plots`:

```markdown
#### Fukui Function Analysis

**Tools > Fukui Function Analysis** shows where a molecule is most susceptible to nucleophilic (f⁺), electrophilic (f⁻) and radical (f⁰) attack, and the dual descriptor Δf = f⁺ − f⁻ (positive: nucleophilic attack; negative: electrophilic attack). The results indicate *likely* reactive regions. They are not predictions of reaction outcome: sterics, solvent and barriers are not included.

**Before you start (Gaussian).**
1. Optimize the molecule in the state in which it reacts (N electrons) and check that it is a minimum.
2. At that same geometry, run three single points with the same method and basis set: N, N+1 and N−1. Do not re-optimize the ions.

   | Molecule | N | N+1 | N−1 |
   |---|---|---|---|
   | neutral, closed shell | 0 1 | −1 2 | 1 2 |
   | neutral radical | 0 2 | −1 1 | 1 1 |

3. Run `formchk` on each `.chk`.

**Setup tab.**
- The three `.fchk` files and an optional output folder (blank = next to the N file).
- Fukui isovalue: blank = 0.005 a.u.
- ρ(N) mapping surface: blank = 0.01 a.u.
- Grid spacing: blank = 0.25 bohr.
- Condensed-index source:
  - *Run Multiwfn (Hirshfeld)* uses Multiwfn.exe, found from the typed path, the `Multiwfnpath` folder or `PATH`;
  - *Import Multiwfn CDFT.txt* reads a file from your own Multiwfn run (main function 22, option 2);
  - *Skip condensed indices*.

  If Multiwfn is missing or fails, the run still finishes and a NOTE explains why.
- Write cube files: on by default.

**Run** checks the three files (same atoms, geometry, basis set and method; N+1 has one electron more, N−1 one fewer), then computes everything. The Results tab shows the report, which is also saved as `<name>_Fukui.txt`. It contains:
- the global indices (IP, EA, electronegativity, hardness = I − A, softness, electrophilicity, nucleophilicity);
- grid integrals: about 1 for f⁺/f⁻/f⁰ and 0 for Δf, a quick check that the files are not swapped;
- the condensed table q(N), q(N+1), q(N−1), f⁻, f⁺, f⁰, CDD, s⁻, s⁺, s⁰;
- a site ranking for each kind of attack;
- the surface extrema of each function on ρ(N) = 0.01 a.u., with the nearest atom;
- NOTES, e.g. degenerate frontier orbitals, diffuse basis functions, negative condensed values, unusual multiplicities.

Cubes: `<name>_fplus.cub`, `_fminus.cub`, `_fzero.cub`, `_dual.cub` (a.u. on a bohr grid; open them in VMD or Multiwfn).

**Show 3D.** Pick a function and a display:
- *Isosurfaces* at ±isovalue: blue positive, red negative.
- *Mapped on ρ(N) surface*: colour = value, with extrema markers (cyan minimum, orange maximum).

Canvas background and **Export image...** work as in the ESP tool.

**Reading the results.**
- Compare only atoms that can actually react.
- For electrophilic attack, prefer large f⁻ and the most negative Δf. For nucleophilic attack, prefer the most positive Δf: f⁺ alone is the least reliable of the descriptors.
- If the descriptors disagree, compute the transition states.
- For degenerate HOMO/LUMO (benzene, C₆₀), use N+p/N−q states or Multiwfn's orbital-weighted Fukui functions (see the NOTE).

Theory and limitations: *Formulas and Methods (Fukui Function tool)* below, and `docs/superpowers/specs/2026-10-16-fukui-function-design.md`.
```

- [ ] **Step 2: USER_MANUAL formulas section.** Insert before `## Formulas and Methods (Solvation box builders)`:

```markdown
## Formulas and Methods (Fukui Function tool)

**Real-space functions** (a.u. = e/bohr³), with p = q = 1 for ordinary N±1 states:
f⁺ = [ρ_{N+p} − ρ_N]/p, f⁻ = [ρ_N − ρ_{N−q}]/q, f⁰ = (f⁺ + f⁻)/2, Δf = f⁺ − f⁻.

The densities are the `Total SCF Density` of each `.fchk`, evaluated exactly at every grid node. On phenol they agree with Multiwfn 2026.7.15's exported cubes to 8×10⁻⁶ a.u.

The grid grows until every face has ρ_N < ρ_surface/10 and |f⁺|, |f⁻| < isovalue/10. Isosurfaces use marching tetrahedra with bisection of the exact function. Surface extrema use Multiwfn's rule: strict comparison against the first and second vertex rings (as in the ESP tool).

**Condensed indices** (Hirshfeld charges q_A from Multiwfn main function 22, option 2):
f⁻_A = q_A(N−1) − q_A(N), f⁺_A = q_A(N) − q_A(N+1), f⁰_A = [q_A(N−1) − q_A(N+1)]/2, CDD_A = f⁺_A − f⁻_A.
Each imported row is checked against these relations to ±2×10⁻⁴ (`CDFT.txt` prints 4 decimals). Condensed local softness is s_A = S·f_A.

**Global indices** (Hartree, Multiwfn conventions):

| Index | Formula |
|---|---|
| vertical IP | I = E(N−1) − E(N) |
| vertical EA | A = E(N) − E(N+1) |
| electronegativity | χ = (I+A)/2 |
| chemical potential | μ = −χ |
| hardness | **η = I − A** (no ½) |
| softness | S = 1/η |
| electrophilicity | ω = μ²/2η |
| nucleophilicity | N_Nu = E_HOMO(N) − E_HOMO(TCE), with E_HOMO(TCE) = −0.335198 Ha |

E_HOMO(TCE) is the B3LYP/6-31G\* value. At another level, recompute TCE.

**Not provided** (use Multiwfn):
- Fukui potential and dual-descriptor potential;
- orbital-weighted Fukui functions;
- bond dual descriptor;
- ω_cubic and ε;
- condensed indices for N+p/N−q;
- averages on local atomic surfaces;
- orbital contributions;
- spin-polarized Fukui functions.
```

- [ ] **Step 3: readme.md.**
  - After the ESP capability line (`- Electrostatic surface potential (ESP) mapped on ...`), add `- Fukui functions f+, f-, f0 and the dual descriptor from Gaussian N/N+1/N-1 .fchk files, with surface extrema, cube export and Multiwfn (Hirshfeld) condensed indices`.
  - After the ESP module line (`- electrostaticSurfacePotential.py, ...`), add `` - `fukuiFunctionAnalyze.py` and `cubeFile.py` for the Fukui function tool (optional external program: Multiwfn, for condensed indices) ``.
  - In the accepted-files list, after the `.fchk` line, add `` - Multiwfn `CDFT.txt` (main function 22, option 2) and Gaussian cube files written by the Fukui tool ``.

- [ ] **Step 4: CLAUDE.md.** After the ESP bullet (the one that ends with `…GaussView's bar reads ±2.072e-2 because it interpolates from a coarse cube.`), add one bullet:

```markdown
- **[fukuiFunctionAnalyze.py](fukuiFunctionAnalyze.py) / [cubeFile.py](cubeFile.py)** — **Tools > Fukui Function Analysis**. Spec: [docs/superpowers/specs/2026-10-16-fukui-function-design.md](docs/superpowers/specs/2026-10-16-fukui-function-design.md).
  - **Inputs and checks.** Three Gaussian `.fchk` files: N, N+p, N−q. p and q come from the electron counts and are 1 for ordinary states. `validate_states` rejects different atoms, a geometry shift > 1e-5 bohr ("run single points, do not re-optimize the ions"), a different basis, a different method after `strip_reference_prefix` (R/U/RO), and swapped files. An unusual multiplicity is a NOTE, not an error.
  - **Real-space functions** f⁺, f⁻, f⁰ and DD are computed **natively** from `gaussianWavefunction.electron_density`. They match Multiwfn's cubes to ≤ 8e-6 a.u. (tested against `tests/data/fukui/multiwfn_reference/`).
  - **Condensed (Hirshfeld) indices are never reimplemented.** They come from Multiwfn main function 22 option 2: either `run_multiwfn_cdft`, the scripted external engine, or an imported `CDFT.txt`, read by `parse_cdft_text` (which also accepts the older `Vertical IP` label). Each row is checked against q(N±1) by `condensed_inconsistencies`.
  - **Multiwfn scripting rule** (verified with 2026.7.15): pass `N.fch` as argv, copy `settings.ini` into the working folder, and add a leading blank line only when no `settings.ini` exists. `Multiwfnpath` did **not** make the Windows build find `settings.ini`. A missing or failing Multiwfn is a NOTE, never an abort (the cubegen rule).
  - **Global indices** use Multiwfn's conventions: **η = I − A** without ½, and E_HOMO(TCE) = −0.335198 Ha.
  - **Grid and display.** The grid is sized by one envelope, max(ρ/ρ_iso, |f±|/f_iso). The surface is ρ(N) = 0.01 a.u. (manual §4.12.4), with Multiwfn-rule extrema. Lobes are ±isovalue meshes, shown through `ESPSurfaceViewer`'s optional `caption=` (Task 4 of the plan), which leaves the ESP tool unchanged.
  - **Not implemented** (spec §7): DDP/Fukui potential, orbital-weighted forms, BDD, condensed N+p/N−q, local-surface averages, spin-polarized CDFT.
```

- [ ] **Step 5: Full test suite**

Run: `..\Scripts\python.exe -m pytest tests/ -q`
Expected: only the 4 baseline failures; everything else passes.

- [ ] **Step 6: Live smoke test** (needs a desktop session; this step was **not** executed when the plan was written).
1. `..\Scripts\python.exe gqteaWinToga.py` → Tools → **Fukui Function Analysis**.
2. Select `tests\data\fukui\phenol_N.fchk`, `phenol_Np1.fchk` and `phenol_Nm1.fchk`, and an output folder in `%TEMP%`.
3. Press **Run**. Expected: about 10–15 s. Results tab, with Hardness 10.4129 eV. The condensed table matches `phenol_CDFT.txt` when Multiwfn is installed; otherwise there is a NOTE. The ranking for "electrophilic attack: most negative CDD" starts with `O12, C1`.
4. **Show 3D** with f⁻ and *Isosurfaces*: blue lobes on O12 and the ortho/para carbons, above and below the ring plane. Then *Mapped on ρ(N) surface*: an orange global-maximum marker near the ring or O. Switch the background to White and **Export image...**.
5. Confirm the window fits a laptop screen (900×700) and that Close also closes the 3D window.

- [ ] **Step 7: Commit**

```bash
git add USER_MANUAL.md readme.md CLAUDE.md
git commit -m "Docs: Fukui Function Analysis tool (manual, formulas, readme, CLAUDE.md)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Self-review (done when the plan was written)

- **Spec coverage.**
  - §1 definitions: Tasks 5–6 formulas and the report wording.
  - §2 states: `validate_states` (Task 5) and the help text (Task 7).
  - §3 Multiwfn: the engine and parser (Task 5) and the reference (Task 1).
  - §4 visualization: lobes, surface and extrema (Task 6), viewer (Tasks 4 and 7).
  - §5 limitations: NOTES (Tasks 5–6) and docs (Task 8).
  - §6 checklist: USER_MANUAL (Task 8).
  - §7 decisions: all tasks.
- **Placeholders.** None. All code is reproduced from executed files.
- **Names.** The same in every task: `FUNCTIONS`, `CONDENSED_*`, `DISPLAY_*`, `FukuiResults` fields, `run()` parameters (pinned by `test_every_run_parameter_is_supplied_by_the_ui`), and `caption=`.
