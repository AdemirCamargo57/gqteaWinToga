# ESP Mapped on Electron-Density Isosurfaces — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the empty `electrostaticSurfacePotential.py` stub into a working tool that reads a Gaussian `.fchk`/`.fch`, builds electron-density isosurfaces, maps the exact ESP onto them, finds Multiwfn-identical surface extrema, writes them to a file and shows a GaussView-style colour-mapped surface in 3D.

**Architecture:** Three Toga-free numerical modules (`gaussianWavefunction` → `espIntegrals` → `molecularSurface`), an OpenGL viewer subclass (`espViewer`, built on two new no-op hooks in `molecularViewer`), and the tool module (`electrostaticSurfacePotential`: headless `ESPSurfaceAnalyser` + Toga `ElectrostaticSurfacePotentialUI`). An optional `espNumba` kernel is used automatically when numba is installed.

**Tech Stack:** Python 3.13, numpy 2.4.2, scipy 1.17.0, toga-winforms 0.5.3, glfw + PyOpenGL (legacy fixed pipeline, as in `molecularViewer.py`), pytest. Optional: numba 0.68.0. Test reference data from Gaussian 09W `cubegen`.

**Spec:** [docs/superpowers/specs/2026-10-14-esp-surface-design.md](../specs/2026-10-14-esp-surface-design.md) — read §0 (decisions), §2 (Multiwfn rules), §3 (GaussView conventions) and §11 (measured spike results) first.

## Global Constraints

- Work on branch `feature/esp-surface` (already created; the spec is committed there). Run every command from `C:\AJCamargo\gqteaWinToga\venv\src` with the venv's interpreter `..\Scripts\python.exe`.
- **All code in this plan was executed before the plan was written** (scratch tree, 2026-10-14): 97 tests passed + 8 numba tests skipped without numba, 104 passed with numba; the 149 existing `test_molecularViewer.py`/`test_molecularDesign.py` tests pass against the modified viewer. Copy code **exactly**; if something must change, change the test first and say why in the commit message.
- **Baseline has 4 pre-existing, unrelated failures** — do not fix them in this plan, do not count them as regressions: `tests/test_meanResidenceTime.py::TestAdvancedOptionsExposed::test_every_run_parameter_is_supplied_by_the_ui` (imports a renamed module) and three in `tests/test_plotter.py` (`TestTableFigureCustomisation::test_each_column_gets_its_own_style`, `::test_styled_figure_renders_with_matplotlib`, `TestTableColumnSelection::test_styles_stay_attached_to_their_column`). Expected suite result after the plan: those 4 failed, everything else passed.
- No new **required** dependency. numba is optional (Task 8) and must never be imported outside `espNumba.py`.
- `gaussianWavefunction.py`, `espIntegrals.py`, `molecularSurface.py`, `espNumba.py` must not import Toga or OpenGL (a subprocess test enforces it for the first three).
- Units: lengths in **bohr** inside the numerical modules; Å only at the UI/viewer/file boundary via `coordinateConverter.BOHR_RADIUS_ANGSTROM = 0.529177210544`. ESP in a.u.; `HARTREE_TO_KCAL_MOL = 627.5094740631`.
- Default isovalue **0.000400 a.u.**; default grid spacing **0.25 bohr**; density section **`Total SCF Density`** only.
- Extrema rule is **Multiwfn's, unchanged**: strict comparison against 1st- and 2nd-ring vertex neighbours, no tolerance, no automatic merge, vertex positions.
- UI errors go through dialogs (`window.dialog(toga.ErrorDialog(...))`); nothing raises out of an async handler.
- Commit messages end with: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`

**Scratch source of the verified files** (this session): `C:\Users\ajc\AppData\Local\Temp\claude\c--AJCamargo-gqteaWinToga-venv-src\d12a792d-4775-412e-89f0-20172e75b918\scratchpad\final\` — referred to below as `$SCRATCH`. Every file is also reproduced in full in this plan; the plan text is authoritative if the scratch folder is gone.

## File Structure

| File | Status | Responsibility |
|---|---|---|
| `tests/data/fchk/*.fchk`, `tests/data/fchk/inputs/*.gjf` | create | 8 fixtures (user's Eteno + 7 G09W jobs) and the inputs that made them |
| `tests/data/esp_reference/*.txt`, `make_reference.py` | create | cubegen ρ/ESP at 300 points per fixture; generator script (needs Gaussian; never run by tests) |
| `tests/esp_fixtures.py` | create | fixture paths + reference loader (not a test module) |
| `gaussianWavefunction.py` | create | fchk parser, basis expansion, pure→Cartesian, ρ(r) |
| `espIntegrals.py` | create | Boys function, McMurchie–Davidson ESP, threads, overlap check, kcal constant |
| `molecularSurface.py` | create | grid box, marching tetrahedra, Multiwfn vertex elimination, mesh topology, extrema |
| `help.py` | modify | `AtomicData.vdw_radii`, `DEFAULT_VDW_RADIUS`, `HelpGqteaWin.help_esp_surface` |
| `molecularViewer.py` | modify | `window_title` attribute + two no-op hooks called by `main_loop` |
| `espViewer.py` | create | GaussView colour maths (pure) + `ESPSurfaceViewer(MolecularViewer)` |
| `electrostaticSurfacePotential.py` | replace stub | `ESPSurfaceAnalyser`, cubegen engine, extrema file, `ElectrostaticSurfacePotentialUI` |
| `espNumba.py` | create (optional task) | compiled ESP kernel |
| `USER_MANUAL.md`, `readme.md`, `CLAUDE.md`, `installation_requirements.txt` | modify | documentation |
| `gqteaWinToga.py` | commit user's edit | already imports and registers `ElectrostaticSurfacePotentialUI` (uncommitted user change; committed in Task 6) |

---

### Task 1: Test fixtures and cubegen reference data

**Files:**
- Create: `tests/data/fchk/{eteno_b3lyp_631g,water_631gd_6d,water_631gd_5d,water_ccpvtz,water_ccpvqz,h2co_6311g2df_10f,ch3_ub3lyp_631gd,benzene_631gd}.fchk`
- Create: `tests/data/fchk/inputs/*.gjf` (same 8 names)
- Create: `tests/data/esp_reference/make_reference.py`, `tests/data/esp_reference/*.txt` (8 files)

**Interfaces:**
- Produces: fixture files consumed by `tests/esp_fixtures.py` (Task 2). Reference file columns: `x_A y_A z_A rho_au esp_au`, 300 rows, `#` comments.

What each fixture exercises: `eteno_b3lyp_631g` = the user's `C:\AJCamargo\QCMM_2sem_2026\2026_10_08\Eteno.fch` (RB3LYP/6-31G, SP shells); `water_631gd_6d` Cartesian d; `water_631gd_5d` pure d; `water_ccpvtz` pure f; `water_ccpvqz` pure g; `h2co_6311g2df_10f` Cartesian f (fch f order); `ch3_ub3lyp_631gd` open shell; `benzene_631gd` only for a manual timing check. Total ≈ 1 MB.

- [x] **Step 1: Copy the fixtures**

```bash
SCRATCH="/c/Users/ajc/AppData/Local/Temp/claude/c--AJCamargo-gqteaWinToga-venv-src/d12a792d-4775-412e-89f0-20172e75b918/scratchpad/final"
mkdir -p tests/data/fchk/inputs tests/data/esp_reference
cp "$SCRATCH"/tests/data/fchk/*.fchk tests/data/fchk/
cp "$SCRATCH"/tests/data/fchk/inputs/*.gjf tests/data/fchk/inputs/
cp "$SCRATCH"/tests/data/esp_reference/*.txt tests/data/esp_reference/
cp "$SCRATCH"/tests/data/esp_reference/make_reference.py tests/data/esp_reference/
ls tests/data/fchk tests/data/esp_reference
```

Expected: 8 `.fchk`, 8 `.txt`, `make_reference.py`, and `inputs/` with 8 `.gjf`.

*If the scratch folder no longer exists:* run each `inputs/*.gjf` with `C:\G09W\g09.exe name.gjf name.out` (set `GAUSS_EXEDIR=C:\G09W`, `GAUSS_SCRDIR` to a temp folder), `C:\G09W\formchk.exe name.chk name.fchk`, copy `Eteno.fch` as `eteno_b3lyp_631g.fchk`, then Step 2's script. The `.gjf` files are:

`inputs/benzene_631gd.gjf`:

```text
%chk=benzene_631gd.chk
%mem=500MB
#p B3LYP/6-31G(d) scf=tight

benzene

0 1
C 1.397000 0.000000 0.0
C 0.698500 1.209837 0.0
C -0.698500 1.209837 0.0
C -1.397000 0.000000 0.0
C -0.698500 -1.209837 0.0
C 0.698500 -1.209837 0.0
H 2.481000 0.000000 0.0
H 1.240500 2.148609 0.0
H -1.240500 2.148609 0.0
H -2.481000 0.000000 0.0
H -1.240500 -2.148609 0.0
H 1.240500 -2.148609 0.0
```

`inputs/ch3_ub3lyp_631gd.gjf`:

```text
%chk=ch3_ub3lyp_631gd.chk
%mem=500MB
#p UB3LYP/6-31G(d) scf=tight

ch3_ub3lyp_631gd

0 2
C  0.000000  0.000000  0.000000
H  1.079000  0.000000  0.000000
H -0.539500  0.934400  0.000000
H -0.539500 -0.934400  0.000000
```

`inputs/eteno_b3lyp_631g.gjf`:

```text
%chk=C:\AJCamargo\QCMM_2sem_2026\2026_10_08\Eteno.chk
# opt b3lyp/6-31g geom=connectivity

Title Card Required

0 1
 C                  0.00000000   -0.67760000    0.00000000
 H                  0.92664718   -1.21260000    0.00000000
 H                 -0.92664718   -1.21260000    0.00000000
 C                 -0.00000000    0.67760000    0.00000000
 H                 -0.92664718    1.21260000    0.00000000
 H                  0.92664718    1.21260000   -0.00000000

 1 2 1.0 3 1.0 4 2.0
 2
 3
 4 5 1.0 6 1.0
 5
 6
```

`inputs/h2co_6311g2df_10f.gjf`:

```text
%chk=h2co_6311g2df_10f.chk
%mem=500MB
#p HF/6-311G(2df) 6D 10F scf=tight

h2co_6311g2df_10f

0 1
C  0.000000  0.000000 -0.529000
O  0.000000  0.000000  0.676000
H  0.000000  0.937000 -1.116000
H  0.000000 -0.937000 -1.116000
```

`inputs/water_631gd_5d.gjf`:

```text
%chk=water_631gd_5d.chk
%mem=500MB
#p HF/6-31G(d) 5D scf=tight

water_631gd_5d

0 1
O  0.000000  0.000000  0.117300
H  0.000000  0.757200 -0.469200
H  0.000000 -0.757200 -0.469200
```

`inputs/water_631gd_6d.gjf`:

```text
%chk=water_631gd_6d.chk
%mem=500MB
#p HF/6-31G(d) 6D scf=tight

water_631gd_6d

0 1
O  0.000000  0.000000  0.117300
H  0.000000  0.757200 -0.469200
H  0.000000 -0.757200 -0.469200
```

`inputs/water_ccpvqz.gjf`:

```text
%chk=water_ccpvqz.chk
%mem=500MB
#p HF/cc-pVQZ scf=tight

water_ccpvqz

0 1
O  0.000000  0.000000  0.117300
H  0.000000  0.757200 -0.469200
H  0.000000 -0.757200 -0.469200
```

`inputs/water_ccpvtz.gjf`:

```text
%chk=water_ccpvtz.chk
%mem=500MB
#p HF/cc-pVTZ scf=tight

water_ccpvtz

0 1
O  0.000000  0.000000  0.117300
H  0.000000  0.757200 -0.469200
H  0.000000 -0.757200 -0.469200
```


- [x] **Step 2: Create the reference generator** — `tests/data/esp_reference/make_reference.py` (already copied; this is its content):

```python
"""Regenerate the cubegen reference values used by the ESP tests.

NOT run by the test suite (it needs Gaussian). Run it by hand, from the
source directory, only when a fixture in tests/data/fchk changes:

    python tests/data/esp_reference/make_reference.py C:\\G09W\\cubegen.exe

For every tests/data/fchk/<name>.fchk it draws 300 reproducible points,
0.5-8 bohr from a random atom, asks cubegen for the SCF density and the ESP
there, and writes tests/data/esp_reference/<name>.txt with columns
    x_A  y_A  z_A  rho_au  esp_au
cubegen reads free-format point coordinates from stdin **in Angstrom**
(verified 2026-10-14 with G09W: a point typed as 0.1173 lands on the water
oxygen at 0.1173 A, not 0.1173 bohr), and echoes them back unchanged.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

BOHR_RADIUS_ANGSTROM = 0.529177210544
HERE = Path(__file__).resolve().parent
FCHK_DIR = HERE.parent / "fchk"
N_POINTS = 300
SEED = 20261014


def fchk_coordinates_bohr(path: Path) -> np.ndarray:
    lines = path.read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("Current cartesian coordinates"))
    count = int(lines[start].split()[-1])
    values = []
    k = start + 1
    while len(values) < count:
        values += [float(x) for x in lines[k].split()]
        k += 1
    return np.array(values).reshape(-1, 3)


def sample_points_angstrom(coords_bohr: np.ndarray) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    atoms = coords_bohr[rng.integers(len(coords_bohr), size=N_POINTS)]
    direction = rng.normal(size=(N_POINTS, 3))
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    distance = rng.uniform(0.5, 8.0, size=N_POINTS)
    return (atoms + direction * distance[:, None]) * BOHR_RADIUS_ANGSTROM


def run_cubegen(cubegen: str, fchk: Path, kind: str, points: np.ndarray) -> np.ndarray:
    env = dict(os.environ, GAUSS_EXEDIR=str(Path(cubegen).parent))
    with tempfile.TemporaryDirectory() as tmp:
        cube = Path(tmp) / "out.cub"
        stdin = "".join(f"{x:.10f} {y:.10f} {z:.10f}\n" for x, y, z in points)
        subprocess.run([cubegen, "0", f"{kind}=SCF", str(fchk), str(cube), "-5", "h"],
                       input=stdin, text=True, env=env, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        lines = cube.read_text().splitlines()
    natoms = abs(int(lines[2].split()[0]))
    data = np.array([[float(v) for v in l.split()] for l in lines[6 + natoms:] if l.strip()])
    if data.shape[0] != len(points) or not np.allclose(data[:, :3], points, atol=1e-8):
        raise RuntimeError(f"cubegen did not echo the requested points for {fchk.name}")
    return data[:, 3]


def main(cubegen: str) -> None:
    for fchk in sorted(FCHK_DIR.glob("*.fchk")):
        points = sample_points_angstrom(fchk_coordinates_bohr(fchk))
        rho = run_cubegen(cubegen, fchk, "density", points)
        esp = run_cubegen(cubegen, fchk, "potential", points)
        out = HERE / f"{fchk.stem}.txt"
        with out.open("w") as handle:
            handle.write(f"# cubegen reference for {fchk.name}: density=SCF and potential=SCF\n")
            handle.write("# x_A y_A z_A rho_au esp_au   (points sent to cubegen in Angstrom)\n")
            for (x, y, z), r, v in zip(points, rho, esp):
                handle.write(f"{x:.10f} {y:.10f} {z:.10f} {r:.12f} {v:.12f}\n")
        print(f"wrote {out.name}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else r"C:\G09W\cubegen.exe")
```

- [x] **Step 3: Check the reference data regenerates identically** (only possible with Gaussian installed; skip otherwise)

Run: `..\Scripts\python.exe tests/data/esp_reference/make_reference.py C:\G09W\cubegen.exe` then `git status --short tests/data/esp_reference`
Expected: 8 `wrote …txt` lines; files byte-identical to the copied ones (git shows them only as untracked, no diff after a second run).

- [x] **Step 4: Commit**

```bash
git add tests/data/fchk tests/data/esp_reference
git commit -m "Add Gaussian fchk fixtures and cubegen ESP reference data

Eight formatted checkpoints covering SP shells, 5D/6D, pure f and g,
Cartesian 10F and open shell, with Gaussian 09 cubegen density and ESP
at 300 points each (cubegen reads -5 points in Angstrom).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: `gaussianWavefunction.py` — fchk reader, basis and density

**Files:**
- Create: `tests/esp_fixtures.py`
- Create: `tests/test_gaussianWavefunction.py`
- Create: `gaussianWavefunction.py`

**Interfaces:**
- Consumes: `help.AtomicData.atomic_numbers` (exists).
- Produces (used by Tasks 3, 4, 6):
  - `read_fchk(path) -> FchkWavefunction` (raises `ValueError` naming the problem)
  - `FchkWavefunction` fields: `path, title, atomic_numbers, nuclear_charges, coordinates (n,3) bohr, shells: List[Shell], density, n_electrons, charge, multiplicity, dipole, transform, cart_density`; property `n_atoms`; method `atom_symbols() -> List[str]`
  - `Shell(l, pure, center, atom, exponents, coefficients)` with `n_cart`, `n_functions`
  - `electron_density(wfn, points (N,3) bohr, chunk=20000) -> (N,)`
  - `cartesian_components(l)`, `component_angular_norm(powers)`, `primitive_norm(alpha, l)`, `pure_to_cartesian(l)`, `read_fchk_sections(text)`, `DENSITY_SECTION = "Total SCF Density"`

- [x] **Step 1: Write the fixture helper** — `tests/esp_fixtures.py`:

```python
"""Shared paths and loaders for the ESP test fixtures (not a test module)."""
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent / "data"
FCHK_DIR = DATA / "fchk"
REFERENCE_DIR = DATA / "esp_reference"
BOHR_RADIUS_ANGSTROM = 0.529177210544

# Every fixture except benzene (kept for the integration/timing test only).
REFERENCE_FIXTURES = [
    "eteno_b3lyp_631g",     # user's file: RB3LYP/6-31G, SP shells
    "water_631gd_6d",       # Cartesian d
    "water_631gd_5d",       # pure d
    "water_ccpvtz",         # pure f
    "water_ccpvqz",         # pure g
    "h2co_6311g2df_10f",    # Cartesian f (fchk f order)
    "ch3_ub3lyp_631gd",     # open shell, UB3LYP
]


def fchk_path(name: str) -> Path:
    return FCHK_DIR / f"{name}.fchk"


def load_reference(name: str):
    """(points_bohr, rho_au, esp_au) from tests/data/esp_reference/<name>.txt."""
    data = np.loadtxt(REFERENCE_DIR / f"{name}.txt", comments="#")
    return data[:, :3] / BOHR_RADIUS_ANGSTROM, data[:, 3], data[:, 4]
```

- [x] **Step 2: Write the failing tests** — `tests/test_gaussianWavefunction.py`:

```python
"""Tests for the fchk reader, basis expansion and electron density."""
import numpy as np
import pytest

import gaussianWavefunction as gw
from esp_fixtures import REFERENCE_FIXTURES, fchk_path, load_reference

MINI_FCHK = """Title line
SP        RHF                                                         STO-3G
Number of atoms                            I                2
Charge                                     I                0
Total Energy                               R     -1.117506712920000E+00
Atomic numbers                             I   N=           2
           1           1
Nuclear charges                            R   N=           2
  1.00000000E+00  1.00000000E+00
"""


class TestSectionReader:
    def test_scalars_and_arrays(self):
        s = gw.read_fchk_sections(MINI_FCHK)
        assert s["Number of atoms"] == 2
        assert s["Total Energy"] == pytest.approx(-1.11750671292)
        assert s["Atomic numbers"].tolist() == [1, 1]
        assert s["Nuclear charges"].tolist() == [1.0, 1.0]

    def test_truncated_array_is_rejected(self):
        text = MINI_FCHK.replace("           1           1\n", "           1\n")
        with pytest.raises(ValueError, match="Atomic numbers.*truncated"):
            gw.read_fchk_sections(text)

    def test_too_short_file(self):
        with pytest.raises(ValueError, match="too short"):
            gw.read_fchk_sections("x\n")

    def test_missing_section_is_named(self, tmp_path):
        bad = tmp_path / "bad.fchk"
        bad.write_text(MINI_FCHK)
        with pytest.raises(ValueError, match="Current cartesian coordinates"):
            gw.read_fchk(bad)

    def test_unreadable_file(self, tmp_path):
        with pytest.raises(ValueError, match="Cannot read"):
            gw.read_fchk(tmp_path / "absent.fchk")


class TestEtenoFile:
    def test_molecule_and_basis(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))
        assert w.atom_symbols() == ["C", "H", "H", "C", "H", "H"]
        assert w.n_electrons == 16 and w.charge == 0 and w.multiplicity == 1
        assert w.density.shape == (26, 26)
        # 6-31G: C = 1s + 2 SP shells, H = 2 s shells; each SP becomes s + p.
        assert len(w.shells) == 18
        assert np.allclose(w.density, w.density.T)
        assert np.abs(w.dipole).max() < 1e-10  # centrosymmetric: zero dipole

    def test_coordinates_are_bohr(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))
        cc = np.linalg.norm(w.coordinates[0] - w.coordinates[3])
        assert cc * gw_bohr() == pytest.approx(1.336, abs=0.01)  # C=C in Angstrom


def gw_bohr():
    return 0.529177210544


class TestPureToCartesian:
    @pytest.mark.parametrize("l", [1, 2, 3, 4])
    def test_rows_are_orthonormal(self, l):
        t = gw.pure_to_cartesian(l)
        comps = gw.cartesian_components(l)
        s = np.array([[gw._same_centre_overlap(a, b) for b in comps] for a in comps])
        assert t.shape == (2 * l + 1, len(comps))
        assert np.allclose(t @ s @ t.T, np.eye(2 * l + 1), atol=1e-12)

    def test_d_shell_gaussian_order(self):
        # Rows D0, D+1, D-1, D+2, D-2 over XX, YY, ZZ, XY, XZ, YZ.
        t = gw.pure_to_cartesian(2)
        nz = [tuple(np.nonzero(np.abs(row) > 1e-12)[0]) for row in t]
        assert nz == [(0, 1, 2), (4,), (5,), (0, 1), (3,)]
        assert t[3, 0] == pytest.approx(-t[3, 1])  # x^2 - y^2
        assert t[1, 4] > 0 and t[2, 5] > 0 and t[4, 3] > 0


class TestDensity:
    @pytest.mark.parametrize("name", REFERENCE_FIXTURES)
    def test_matches_cubegen_in_the_surface_band(self, name):
        """cubegen's density carries its own screening error (worst seen
        2.5e-4 relative); the ESP test is the stringent check of the basis."""
        w = gw.read_fchk(fchk_path(name))
        points, rho_ref, _ = load_reference(name)
        band = (rho_ref > 1e-5) & (rho_ref < 1e-2)
        assert band.sum() > 30
        rho = gw.electron_density(w, points[band])
        assert np.max(np.abs(rho - rho_ref[band]) / rho_ref[band]) < 1e-3

    def test_chunking_does_not_change_values(self):
        w = gw.read_fchk(fchk_path("water_ccpvtz"))
        points, _, _ = load_reference("water_ccpvtz")
        assert np.allclose(gw.electron_density(w, points, chunk=7), gw.electron_density(w, points))
```

- [x] **Step 3: Run them to confirm they fail**

Run: `..\Scripts\python.exe -m pytest tests/test_gaussianWavefunction.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'gaussianWavefunction'`.

- [x] **Step 4: Implement** — `gaussianWavefunction.py`:

```python
"""Gaussian formatted checkpoint (.fchk / .fch) reader, basis set and density.

Pure numpy, no Toga: the electrostatic surface potential tool builds on it
and the tests drive it headless. Every length is in bohr.

Basis conventions follow Multiwfn's fch reader (fileIO.f90 readfch), and
were verified against Gaussian 09 cubegen on water 5D/6D, cc-pVTZ (pure f),
cc-pVQZ (pure g), formaldehyde 10F (Cartesian f), ethylene (SP shells) and
the methyl radical (open shell):

* each Cartesian primitive is normalised on its own, x^i y^j z^k
  exp(-a r^2) with N = (2a/pi)^(3/4) (4a)^(L/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!),
  and multiplied by the fchk contraction coefficient as is (the contracted
  functions then come out normalised: tr(P S) = N_electrons to 1e-8);
* an SP shell (type -1) is an s shell plus a p shell sharing exponents, the
  p part taking ``P(S=P) Contraction coefficients``;
* a pure shell (type <= -2) is mapped onto its Cartesian components by
  ``pure_to_cartesian``, so all integrals run in one Cartesian code path;
* ``Total SCF Density`` (lower triangle, row-wise) is the density used.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import comb, pi, sqrt
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from help import AtomicData

# Cartesian component order written to the fchk for each angular momentum
# (Multiwfn define.f90: s, p, d as Gaussian pop=full; f and g in the fch
# order, which is not the .wfn order). Entries are powers (i, j, k) of x, y, z.
FCHK_CARTESIAN_ORDER: Dict[int, Tuple[Tuple[int, int, int], ...]] = {
    0: ((0, 0, 0),),
    1: ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    2: ((2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)),
    3: ((3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0),
        (2, 0, 1), (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)),
    4: ((0, 0, 4), (0, 1, 3), (0, 2, 2), (0, 3, 1), (0, 4, 0),
        (1, 0, 3), (1, 1, 2), (1, 2, 1), (1, 3, 0), (2, 0, 2),
        (2, 1, 1), (2, 2, 0), (3, 0, 1), (3, 1, 0), (4, 0, 0)),
}
MAX_ANGULAR_MOMENTUM = max(FCHK_CARTESIAN_ORDER)
DENSITY_SECTION = "Total SCF Density"
REQUIRED_SECTIONS = (
    "Atomic numbers", "Nuclear charges", "Current cartesian coordinates",
    "Shell types", "Number of primitives per shell", "Shell to atom map",
    "Primitive exponents", "Contraction coefficients", DENSITY_SECTION,
)
SYMBOL_BY_ATOMIC_NUMBER = {z: s for s, z in AtomicData.atomic_numbers.items()}


def cartesian_components(l: int) -> Tuple[Tuple[int, int, int], ...]:
    return FCHK_CARTESIAN_ORDER[l]


def double_factorial(n: int) -> int:
    """n!! with (-1)!! = 0!! = 1."""
    result = 1
    while n > 1:
        result *= n
        n -= 2
    return result


def component_angular_norm(powers: Tuple[int, int, int]) -> float:
    """1/sqrt((2i-1)!!(2j-1)!!(2k-1)!!), the component part of a Cartesian norm."""
    i, j, k = powers
    return 1.0 / sqrt(double_factorial(2 * i - 1) * double_factorial(2 * j - 1)
                      * double_factorial(2 * k - 1))


def primitive_norm(alpha: np.ndarray, l: int) -> np.ndarray:
    """(2a/pi)^(3/4) (4a)^(l/2), the exponent part of a Cartesian norm."""
    return (2.0 * alpha / pi) ** 0.75 * (4.0 * alpha) ** (l / 2.0)


def _same_centre_overlap(p1, p2) -> float:
    """Overlap of two normalised Cartesian components with one centre and exponent."""
    value = 1.0
    for a, b in zip(p1, p2):
        if (a + b) % 2:
            return 0.0
        value *= double_factorial(a + b - 1)
    return value * component_angular_norm(p1) * component_angular_norm(p2)


def _solid_harmonic_monomials(l: int, m: int) -> Dict[Tuple[int, int, int], float]:
    """Real solid harmonic S_lm as {(i, j, k): coefficient} over x^i y^j z^k.

    Helgaker, Jorgensen & Olsen, Molecular Electronic-Structure Theory,
    eq. 6.4.48; the overall scale is irrelevant (rows are renormalised).
    """
    am = abs(m)
    vm = 0.0 if m >= 0 else 0.5
    terms: Dict[Tuple[int, int, int], float] = {}
    v_values = [vm + n for n in range(int(np.floor(am / 2.0 - vm)) + 1)]
    for t in range((l - am) // 2 + 1):
        for u in range(t + 1):
            for v in v_values:
                two_v = int(round(2 * v))
                coef = ((-1) ** int(round(t + v - vm)) * (0.25 ** t) * comb(l, t)
                        * comb(l - t, am + t) * comb(t, u) * comb(am, two_v))
                key = (2 * t + am - 2 * u - two_v, 2 * u + two_v, l - 2 * t - am)
                terms[key] = terms.get(key, 0.0) + coef
    return terms


def pure_to_cartesian(l: int) -> np.ndarray:
    """T with pure_p = sum_c T[p, c] cart_c (both normalised).

    Rows follow Gaussian's pure order 0, +1, -1, +2, -2, ...; columns the
    fchk Cartesian order of ``cartesian_components(l)``.
    """
    comps = cartesian_components(l)
    index = {c: n for n, c in enumerate(comps)}
    order = [0]
    for m in range(1, l + 1):
        order += [m, -m]
    overlap = np.array([[_same_centre_overlap(a, b) for b in comps] for a in comps])
    rows = []
    for m in order:
        row = np.zeros(len(comps))
        for powers, coef in _solid_harmonic_monomials(l, m).items():
            row[index[powers]] += coef / component_angular_norm(powers)
        row /= sqrt(row @ overlap @ row)
        rows.append(row)
    return np.array(rows)


def read_fchk_sections(text: str) -> Dict[str, object]:
    """Every integer/real scalar and array section of an fchk, keyed by name.

    A header line has the name in columns 1-43, the type letter in column 44
    and either a scalar or ``N=  count``; array values follow on lines that
    start with a blank. Character and logical sections are skipped.
    """
    lines = text.splitlines()
    if len(lines) < 3:
        raise ValueError("The file is too short to be a Gaussian formatted checkpoint.")
    sections: Dict[str, object] = {}
    n = 2  # line 1 = title, line 2 = job type / method / basis
    while n < len(lines):
        line = lines[n]
        n += 1
        if not line.strip() or line[0] == " " or len(line) < 45:
            continue
        name, kind, rest = line[:43].strip(), line[43], line[44:].split()
        if not rest:
            continue
        if rest[0].startswith("N="):
            count = int(rest[0][2:]) if len(rest[0]) > 2 else int(rest[1])
            tokens: List[str] = []
            while n < len(lines) and lines[n][:1] == " ":
                tokens.extend(lines[n].split())
                n += 1
            try:
                if kind == "I":
                    values = np.array([int(t) for t in tokens], dtype=np.int64)
                elif kind == "R":
                    values = np.array([float(t.replace("D", "E")) for t in tokens])
                else:
                    continue
            except ValueError as exc:
                raise ValueError(f"Section '{name}' contains a value that is not a number.") from exc
            if values.size != count:
                raise ValueError(
                    f"Section '{name}' declares {count} values but holds {values.size}; "
                    "the file looks truncated."
                )
            sections[name] = values
        elif kind == "I":
            sections[name] = int(rest[0])
        elif kind == "R":
            sections[name] = float(rest[0].replace("D", "E"))
    return sections


@dataclass
class Shell:
    l: int
    pure: bool
    center: np.ndarray       # bohr
    atom: int                # 0-based atom index
    exponents: np.ndarray
    coefficients: np.ndarray  # fchk contraction coefficients, as written

    @property
    def n_cart(self) -> int:
        return (self.l + 1) * (self.l + 2) // 2

    @property
    def n_functions(self) -> int:
        return 2 * self.l + 1 if self.pure else self.n_cart


@dataclass
class FchkWavefunction:
    path: str
    title: str
    atomic_numbers: np.ndarray
    nuclear_charges: np.ndarray   # effective charges (differ from Z with an ECP)
    coordinates: np.ndarray       # (n_atoms, 3) bohr
    shells: List[Shell]
    density: np.ndarray           # symmetric, in the fchk basis
    n_electrons: int
    charge: int
    multiplicity: int
    dipole: Optional[np.ndarray] = None
    transform: np.ndarray = field(default=None, repr=False)     # fchk -> Cartesian
    cart_density: np.ndarray = field(default=None, repr=False)  # T^T P T

    @property
    def n_atoms(self) -> int:
        return len(self.atomic_numbers)

    def atom_symbols(self) -> List[str]:
        return [SYMBOL_BY_ATOMIC_NUMBER.get(int(z), "X") for z in self.atomic_numbers]


def read_fchk(path) -> FchkWavefunction:
    """Parse a .fchk/.fch file. Every problem is a ValueError naming it."""
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ValueError(f"Cannot read {path}: {exc}") from exc
    sections = read_fchk_sections(text)
    missing = [name for name in REQUIRED_SECTIONS if name not in sections]
    if missing:
        raise ValueError(f"{path.name} is missing fchk section(s): {', '.join(missing)}.")

    coords = np.asarray(sections["Current cartesian coordinates"], dtype=float).reshape(-1, 3)
    z = np.asarray(sections["Atomic numbers"])
    if len(z) != len(coords):
        raise ValueError(f"{path.name}: {len(z)} atomic numbers but {len(coords)} coordinate triples.")
    shell_types = sections["Shell types"]
    nprim = sections["Number of primitives per shell"]
    shell_atom = sections["Shell to atom map"]
    exps = sections["Primitive exponents"]
    coefs = sections["Contraction coefficients"]
    sp_coefs = sections.get("P(S=P) Contraction coefficients")

    shells: List[Shell] = []
    start = 0
    for stype, count, atom in zip(shell_types.tolist(), nprim.tolist(), shell_atom.tolist()):
        sl = slice(start, start + count)
        start += count
        centre = coords[atom - 1]
        if stype == -1:
            if sp_coefs is None:
                raise ValueError(f"{path.name}: SP shell without 'P(S=P) Contraction coefficients'.")
            shells.append(Shell(0, False, centre, atom - 1, exps[sl].copy(), coefs[sl].copy()))
            shells.append(Shell(1, False, centre, atom - 1, exps[sl].copy(), sp_coefs[sl].copy()))
            continue
        l = abs(stype)
        if l > MAX_ANGULAR_MOMENTUM:
            raise ValueError(
                f"{path.name}: basis functions with l = {l} are not supported "
                f"(highest supported is g, l = {MAX_ANGULAR_MOMENTUM})."
            )
        shells.append(Shell(l, stype <= -2, centre, atom - 1, exps[sl].copy(), coefs[sl].copy()))
    if start != exps.size:
        raise ValueError(f"{path.name}: shell primitive counts do not match the exponent list.")

    nbf = sum(s.n_functions for s in shells)
    tri = sections[DENSITY_SECTION]
    if tri.size != nbf * (nbf + 1) // 2:
        raise ValueError(
            f"{path.name}: '{DENSITY_SECTION}' has {tri.size} elements, "
            f"expected {nbf * (nbf + 1) // 2} for {nbf} basis functions."
        )
    density = np.zeros((nbf, nbf))
    density[np.tril_indices(nbf)] = tri
    density = density + density.T - np.diag(np.diag(density))

    lines = text.splitlines()
    wfn = FchkWavefunction(
        path=str(path),
        title=lines[0].strip() if lines else "",
        atomic_numbers=z,
        nuclear_charges=np.asarray(sections["Nuclear charges"], dtype=float),
        coordinates=coords,
        shells=shells,
        density=density,
        n_electrons=int(sections.get("Number of electrons", 0)),
        charge=int(sections.get("Charge", 0)),
        multiplicity=int(sections.get("Multiplicity", 1)),
        dipole=np.asarray(sections["Dipole Moment"]) if "Dipole Moment" in sections else None,
    )
    wfn.transform = basis_transform(shells)
    wfn.cart_density = wfn.transform.T @ density @ wfn.transform
    return wfn


def basis_transform(shells: List[Shell]) -> np.ndarray:
    """Block-diagonal (n_functions, n_cartesian) map onto normalised Cartesians."""
    nbf = sum(s.n_functions for s in shells)
    ncart = sum(s.n_cart for s in shells)
    t = np.zeros((nbf, ncart))
    row = col = 0
    for s in shells:
        block = pure_to_cartesian(s.l) if s.pure else np.eye(s.n_cart)
        t[row:row + s.n_functions, col:col + s.n_cart] = block
        row += s.n_functions
        col += s.n_cart
    return t


def evaluate_cartesian_aos(wfn: FchkWavefunction, points: np.ndarray) -> np.ndarray:
    """(n_points, n_cartesian) values of every normalised Cartesian component."""
    points = np.asarray(points, dtype=float)
    out = np.empty((points.shape[0], sum(s.n_cart for s in wfn.shells)))
    col = 0
    for s in wfn.shells:
        d = points - s.center
        r2 = np.einsum("ij,ij->i", d, d)
        radial = np.exp(-np.outer(r2, s.exponents)) @ (s.coefficients * primitive_norm(s.exponents, s.l))
        for c in cartesian_components(s.l):
            out[:, col] = (radial * component_angular_norm(c)
                           * d[:, 0] ** c[0] * d[:, 1] ** c[1] * d[:, 2] ** c[2])
            col += 1
    return out


def electron_density(wfn: FchkWavefunction, points: np.ndarray, chunk: int = 20000) -> np.ndarray:
    """rho(r) = sum_ab P_ab chi_a(r) chi_b(r) at every point (a.u.)."""
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    rho = np.empty(points.shape[0])
    for a in range(0, points.shape[0], chunk):
        ao = evaluate_cartesian_aos(wfn, points[a:a + chunk])
        rho[a:a + chunk] = np.einsum("ij,ij->i", ao @ wfn.cart_density, ao)
    return rho
```

- [x] **Step 5: Run the tests**

Run: `..\Scripts\python.exe -m pytest tests/test_gaussianWavefunction.py -q`
Expected: `20 passed`.

- [x] **Step 6: Commit**

```bash
git add gaussianWavefunction.py tests/esp_fixtures.py tests/test_gaussianWavefunction.py
git commit -m "Add Gaussian fchk reader, basis expansion and electron density

Follows Multiwfn's fch conventions (per-primitive Cartesian norms, SP
split, pure->Cartesian, fch f/g order); density matches cubegen in the
surface band.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: `espIntegrals.py` — exact ESP

**Files:**
- Create: `tests/test_espIntegrals.py`
- Create: `espIntegrals.py`

**Interfaces:**
- Consumes: `gaussianWavefunction.{FchkWavefunction, cartesian_components, component_angular_norm, primitive_norm}`.
- Produces (used by Tasks 5, 6, 8):
  - `electrostatic_potential(wfn, points_bohr, batches=None, workers=None, progress=None, use_numba=None) -> (N,)` a.u.; `progress(fraction)`
  - `build_hermite_batches(wfn) -> List[HermiteBatch]` (`HermiteBatch(l, p, centers, coeffs)`), `electronic_potential(batches, points)`, `nuclear_potential(wfn, points)`
  - `boys(n_max, t)`, `hermite_coefficients(...)`, `hermite_indices(l)`, `overlap_matrix(wfn)`, `numba_available()`, `default_workers()`
  - `HARTREE_TO_KCAL_MOL = 627.5094740631`
  - module attribute `espNumba` (None until Task 8 adds `espNumba.py`)

- [x] **Step 1: Write the failing tests** — `tests/test_espIntegrals.py` (the `TestNumba` class is skipped until Task 8):

```python
"""Tests for the exact ESP integrals."""
import numpy as np
import pytest
from scipy.integrate import quad

import espIntegrals as ei
import gaussianWavefunction as gw
from esp_fixtures import REFERENCE_FIXTURES, fchk_path, load_reference


class TestBoys:
    @pytest.mark.parametrize("t", [0.0, 1e-12, 0.3, 5.0, 30.0, 80.0])
    def test_against_quadrature(self, t):
        f = ei.boys(8, np.array([t]))
        for n in range(9):
            ref = quad(lambda u: u ** (2 * n) * np.exp(-t * u * u), 0.0, 1.0, epsabs=1e-15)[0]
            assert f[n, 0] == pytest.approx(ref, rel=1e-10, abs=1e-15)


class TestOverlap:
    @pytest.mark.parametrize("name", REFERENCE_FIXTURES)
    def test_trace_ps_is_the_electron_count(self, name):
        w = gw.read_fchk(fchk_path(name))
        s = ei.overlap_matrix(w)
        assert np.trace(w.cart_density @ s) == pytest.approx(w.n_electrons, abs=1e-6)
        # every basis function as written in the fchk is normalised
        assert np.allclose(np.diag(w.transform @ s @ w.transform.T), 1.0, atol=1e-7)


class TestPotential:
    @pytest.mark.parametrize("name", REFERENCE_FIXTURES)
    def test_matches_cubegen(self, name):
        """Gaussian 09 cubegen at 300 points 0.5-8 bohr from the nuclei."""
        w = gw.read_fchk(fchk_path(name))
        points, _, esp_ref = load_reference(name)
        esp = ei.electrostatic_potential(w, points)
        assert np.max(np.abs(esp - esp_ref)) < 1e-6

    def test_threads_do_not_change_values(self):
        w = gw.read_fchk(fchk_path("water_ccpvtz"))
        points, _, _ = load_reference("water_ccpvtz")
        serial = ei.electrostatic_potential(w, points, workers=1, use_numba=False)
        threaded = ei.electrostatic_potential(w, points, workers=6, use_numba=False)
        assert np.allclose(serial, threaded, rtol=0, atol=1e-13)

    def test_far_field_is_the_dipole_potential(self):
        """Neutral water far away: V -> mu . r / r^3 with the fchk dipole."""
        w = gw.read_fchk(fchk_path("water_631gd_6d"))
        r = 150.0
        directions = np.eye(3)
        esp = ei.electrostatic_potential(w, directions * r)
        expected = directions @ w.dipole / r ** 2
        assert np.allclose(esp, expected, atol=2e-6)
        assert abs(esp[2]) > 1e-5  # water's dipole lies along z

    def test_progress_reaches_one(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))
        points, _, _ = load_reference("eteno_b3lyp_631g")
        seen = []
        ei.electrostatic_potential(w, points, workers=1, progress=seen.append)
        assert seen and seen[-1] == pytest.approx(1.0)

    def test_raising_progress_callback_is_ignored(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))
        points, _, esp_ref = load_reference("eteno_b3lyp_631g")

        def boom(_):
            raise RuntimeError("UI went away")

        esp = ei.electrostatic_potential(w, points, workers=2, progress=boom)
        assert np.max(np.abs(esp - esp_ref)) < 1e-6


@pytest.mark.skipif(not ei.numba_available(), reason="numba not installed (optional)")
class TestNumba:
    def test_boys_matches_the_scipy_version(self):
        import espNumba
        ts = np.r_[0.0, 1e-12, np.linspace(0.01, 150.0, 600)]
        ref = ei.boys(8, ts)
        f = np.zeros(9)
        for i, t in enumerate(ts):
            espNumba._boys_into(8, t, f)
            assert np.allclose(f, ref[:, i], rtol=1e-12, atol=0)

    @pytest.mark.parametrize("name", REFERENCE_FIXTURES)
    def test_kernel_matches_numpy(self, name):
        w = gw.read_fchk(fchk_path(name))
        points, _, esp_ref = load_reference(name)
        fast = ei.electrostatic_potential(w, points, use_numba=True)
        slow = ei.electrostatic_potential(w, points, use_numba=False)
        assert np.max(np.abs(fast - slow)) < 1e-10
        assert np.max(np.abs(fast - esp_ref)) < 1e-6


def test_forcing_numba_without_it_is_an_error(monkeypatch):
    monkeypatch.setattr(ei, "espNumba", None)
    w = gw.read_fchk(fchk_path("water_631gd_6d"))
    with pytest.raises(ValueError, match="numba"):
        ei.electrostatic_potential(w, np.zeros((1, 3)) + 5.0, use_numba=True)
```

- [x] **Step 2: Run them to confirm they fail**

Run: `..\Scripts\python.exe -m pytest tests/test_espIntegrals.py -q`
Expected: `ModuleNotFoundError: No module named 'espIntegrals'`.

- [x] **Step 3: Implement** — `espIntegrals.py`:

```python
"""Exact molecular electrostatic potential of a Gaussian-basis wavefunction.

    V(C) = sum_A Z_A / |C - R_A|  -  sum_ab P_ab  int chi_a chi_b / |r - C| dr

Z_A are the fchk ``Nuclear charges`` (effective charges with an ECP, as in
Multiwfn). The electronic term uses McMurchie-Davidson: every primitive
pair's share of the density is reduced **once** to Hermite coefficients
D_tuv (independent of C), so the potential at C is

    sum_pairs  sum_tuv  D_tuv R_tuv(p, P - C),

with R_tuv from the Boys function. Pure numpy/scipy; point chunks run in a
thread pool, since numpy and scipy ufuncs release the GIL. Verified against
Gaussian 09 cubegen to <= 4e-7 a.u. on every fixture in tests/data/fchk.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from math import pi
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.special import gamma, gammainc

from gaussianWavefunction import (FchkWavefunction, cartesian_components,
                                  component_angular_norm, primitive_norm)

try:  # optional compiled kernel (pip install numba)
    import espNumba
except ImportError:
    espNumba = None

# A primitive pair is dropped when its largest possible contribution,
# |prefactor| x max|P block|, is below this (a.u.).
PAIR_SCREEN = 1e-14
# Elements of one (pairs x points) work array; bounds memory per thread.
BLOCK_ELEMENTS = 2_000_000


def boys(n_max: int, t: np.ndarray) -> np.ndarray:
    """Boys functions F_0..F_n_max at t, stacked on a new leading axis.

    F_n(t) = int_0^1 u^(2n) exp(-t u^2) du. F_n_max comes from the
    regularised incomplete gamma function, the rest by the stable downward
    recursion F_(m-1) = (2 t F_m + exp(-t)) / (2m - 1).
    """
    t = np.asarray(t, dtype=float)
    out = np.empty((n_max + 1,) + t.shape)
    small = t < 1e-10
    a = n_max + 0.5
    ts = np.where(small, 1.0, t)
    top = gamma(a) * gammainc(a, ts) / (2.0 * ts ** a)
    out[n_max] = np.where(small, 1.0 / (2 * n_max + 1) - t / (2 * n_max + 3), top)
    emt = np.exp(-t)
    for m in range(n_max, 0, -1):
        out[m - 1] = (2.0 * t * out[m] + emt) / (2 * m - 1)
    return out


def hermite_coefficients(la: int, lb: int, xpa: np.ndarray, xpb: np.ndarray,
                         inv2p: np.ndarray) -> np.ndarray:
    """E[pair, i, j, t] (one Cartesian direction, without the K prefactor).

    E^{i+1,j}_t = 1/(2p) E^{ij}_{t-1} + X_PA E^{ij}_t + (t+1) E^{ij}_{t+1}
    E^{i,j+1}_t = 1/(2p) E^{ij}_{t-1} + X_PB E^{ij}_t + (t+1) E^{ij}_{t+1}
    """
    npair = xpa.shape[0]
    e = np.zeros((npair, la + 1, lb + 1, la + lb + 2))
    e[:, 0, 0, 0] = 1.0
    for i in range(la):
        for t in range(i + 2):
            prev = e[:, i, 0, t - 1] if t > 0 else 0.0
            e[:, i + 1, 0, t] = inv2p * prev + xpa * e[:, i, 0, t] + (t + 1) * e[:, i, 0, t + 1]
    for j in range(lb):
        for i in range(la + 1):
            for t in range(i + j + 2):
                prev = e[:, i, j, t - 1] if t > 0 else 0.0
                e[:, i, j + 1, t] = inv2p * prev + xpb * e[:, i, j, t] + (t + 1) * e[:, i, j, t + 1]
    return e[..., : la + lb + 1]


def hermite_indices(l: int) -> List[Tuple[int, int, int]]:
    """Every (t, u, v) with t + u + v <= l."""
    return [(t, u, v) for t in range(l + 1) for u in range(l + 1 - t) for v in range(l + 1 - t - u)]


@dataclass
class HermiteBatch:
    """All primitive pairs whose angular momenta sum to ``l``."""
    l: int
    p: np.ndarray        # (M,) total exponent
    centers: np.ndarray  # (M, 3) Gaussian product centre P
    coeffs: np.ndarray   # (M, n_tuv): 2 pi/p x K x contraction x density, summed over components


def _primitive_pairs(sa, sb):
    al = np.repeat(sa.exponents, sb.exponents.size)
    be = np.tile(sb.exponents, sa.exponents.size)
    ca = np.repeat(sa.coefficients * primitive_norm(sa.exponents, sa.l), sb.exponents.size)
    cb = np.tile(sb.coefficients * primitive_norm(sb.exponents, sb.l), sa.exponents.size)
    p = al + be
    ab = sa.center - sb.center
    k = np.exp(-al * be / p * (ab @ ab))
    centers = (al[:, None] * sa.center + be[:, None] * sb.center) / p[:, None]
    return al, be, p, ca * cb * k, centers


def build_hermite_batches(wfn: FchkWavefunction) -> List[HermiteBatch]:
    """Contract the density into point-independent Hermite coefficients."""
    shells = wfn.shells
    offsets = np.cumsum([0] + [s.n_cart for s in shells])
    acc: Dict[int, List[Tuple[np.ndarray, np.ndarray, np.ndarray]]] = {}
    for ia, sa in enumerate(shells):
        for ib in range(ia + 1):
            sb = shells[ib]
            block = wfn.cart_density[offsets[ia]:offsets[ia + 1], offsets[ib]:offsets[ib + 1]]
            block_max = np.abs(block).max()
            if block_max == 0.0:
                continue
            _, _, p, ckk, centers = _primitive_pairs(sa, sb)
            factor = 2.0 if ia != ib else 1.0  # the (b, a) block is the transpose
            pre = factor * 2.0 * pi / p * ckk
            keep = np.abs(pre) * block_max > PAIR_SCREEN
            if not keep.any():
                continue
            p, pre, centers = p[keep], pre[keep], centers[keep]
            inv2p = 0.5 / p
            ex, ey, ez = (hermite_coefficients(sa.l, sb.l, centers[:, d] - sa.center[d],
                                               centers[:, d] - sb.center[d], inv2p) for d in range(3))
            l = sa.l + sb.l
            tuv = np.array(hermite_indices(l))
            coeffs = np.zeros((p.size, len(tuv)))
            for a, pa in enumerate(cartesian_components(sa.l)):
                for b, pb in enumerate(cartesian_components(sb.l)):
                    dab = block[a, b] * component_angular_norm(pa) * component_angular_norm(pb)
                    if dab == 0.0:
                        continue
                    coeffs += dab * (ex[:, pa[0], pb[0], tuv[:, 0]] * ey[:, pa[1], pb[1], tuv[:, 1]]
                                     * ez[:, pa[2], pb[2], tuv[:, 2]])
            acc.setdefault(l, []).append((p, centers, coeffs * pre[:, None]))
    return [HermiteBatch(l, np.concatenate([x[0] for x in v]), np.concatenate([x[1] for x in v]),
                         np.concatenate([x[2] for x in v])) for l, v in sorted(acc.items())]


def hermite_integrals(l: int, p: np.ndarray, pc: np.ndarray) -> Dict[Tuple[int, int, int], np.ndarray]:
    """R^0_tuv(p, P - C) for t + u + v <= l; ``pc`` is (3, M, N).

    R^n_000 = (-2p)^n F_n(p |PC|^2), then e.g.
    R^n_{t+1,u,v} = t R^{n+1}_{t-1,u,v} + X_PC R^{n+1}_{t,u,v},
    keeping only two levels of n alive at a time.
    """
    x, y, z = pc
    f = boys(l, p[:, None] * (x * x + y * y + z * z))
    m2p = -2.0 * p[:, None]
    level = {(0, 0, 0): (m2p ** l) * f[l]}
    for n in range(l - 1, -1, -1):
        new = {}
        for (t, u, v) in hermite_indices(l - n):
            if t > 0:
                val = x * level[(t - 1, u, v)]
                if t > 1:
                    val = val + (t - 1) * level[(t - 2, u, v)]
            elif u > 0:
                val = y * level[(t, u - 1, v)]
                if u > 1:
                    val = val + (u - 1) * level[(t, u - 2, v)]
            elif v > 0:
                val = z * level[(t, u, v - 1)]
                if v > 1:
                    val = val + (v - 1) * level[(t, u, v - 2)]
            else:
                val = (m2p ** n) * f[n]
            new[(t, u, v)] = val
        level = new
    return level


def electronic_potential(batches: List[HermiteBatch], points: np.ndarray) -> np.ndarray:
    """int rho(r) / |r - C| dr at every point C (positive; a.u.)."""
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    out = np.zeros(points.shape[0])
    for batch in batches:
        tuv = hermite_indices(batch.l)
        npt = max(1, BLOCK_ELEMENTS // max(batch.p.size, 1))
        for a in range(0, points.shape[0], npt):
            c = points[a:a + npt]
            pc = batch.centers.T[:, :, None] - c.T[:, None, :]
            r = hermite_integrals(batch.l, batch.p, pc)
            stacked = np.stack([r[key] for key in tuv], axis=1)  # (M, n_tuv, N)
            out[a:a + npt] += np.einsum("mk,mkn->n", batch.coeffs, stacked)
    return out


def nuclear_potential(wfn: FchkWavefunction, points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    d = np.linalg.norm(points[:, None, :] - wfn.coordinates[None, :, :], axis=2)
    return (wfn.nuclear_charges[None, :] / d).sum(axis=1)


def default_workers() -> int:
    return max(1, min(12, (os.cpu_count() or 1)))


def numba_available() -> bool:
    return espNumba is not None


def electrostatic_potential(wfn: FchkWavefunction, points: np.ndarray,
                            batches: Optional[List[HermiteBatch]] = None,
                            workers: Optional[int] = None,
                            progress: Optional[Callable[[float], None]] = None,
                            use_numba: Optional[bool] = None) -> np.ndarray:
    """Total ESP (a.u., Hartree/e) at ``points`` (bohr).

    ``use_numba`` None means "if installed" (the compiled kernel is ~5x
    faster on benzene and agrees with numpy to 1e-10). ``progress(fraction)``
    is called as point chunks finish; a callback that raises is ignored.
    """
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    if batches is None:
        batches = build_hermite_batches(wfn)
    if use_numba is None:
        use_numba = numba_available()
    if use_numba and not numba_available():
        raise ValueError("numba is not installed.")
    electronic = espNumba.electronic_potential if use_numba else electronic_potential
    # The numba kernel is parallel itself, so its chunks run one at a time.
    workers = 1 if use_numba else (workers or default_workers())
    n_chunks = max(1, min(len(points), 20 if use_numba else workers * 4))
    chunks = np.array_split(np.arange(len(points)), n_chunks)
    out = np.empty(len(points))
    done = [0]

    def work(idx):
        out[idx] = nuclear_potential(wfn, points[idx]) - electronic(batches, points[idx])
        done[0] += len(idx)
        if progress is not None:
            try:
                progress(done[0] / len(points))
            except Exception:
                pass

    if workers == 1:
        for idx in chunks:
            work(idx)
    else:
        with ThreadPoolExecutor(workers) as pool:
            list(pool.map(work, chunks))
    return out


def overlap_matrix(wfn: FchkWavefunction) -> np.ndarray:
    """Overlap of the normalised Cartesian components (validation: tr(P S) = N)."""
    shells = wfn.shells
    offsets = np.cumsum([0] + [s.n_cart for s in shells])
    s_mat = np.zeros((offsets[-1], offsets[-1]))
    for ia, sa in enumerate(shells):
        for ib, sb in enumerate(shells):
            _, _, p, ckk, centers = _primitive_pairs(sa, sb)
            pre = ckk * (pi / p) ** 1.5
            e = [hermite_coefficients(sa.l, sb.l, centers[:, d] - sa.center[d],
                                      centers[:, d] - sb.center[d], 0.5 / p) for d in range(3)]
            for a, pa in enumerate(cartesian_components(sa.l)):
                for b, pb in enumerate(cartesian_components(sb.l)):
                    s_mat[offsets[ia] + a, offsets[ib] + b] = (
                        component_angular_norm(pa) * component_angular_norm(pb)
                        * np.sum(pre * e[0][:, pa[0], pb[0], 0] * e[1][:, pa[1], pb[1], 0]
                                 * e[2][:, pa[2], pb[2], 0]))
    return s_mat


# 1 Hartree/e in kcal/(mol e): CODATA 2018 E_h = 4.3597447222071e-18 J times
# N_A = 6.02214076e23 / mol, over 4184 J/kcal. Multiwfn uses 627.51.
HARTREE_TO_KCAL_MOL = 627.5094740631
```

- [x] **Step 4: Run the tests**

Run: `..\Scripts\python.exe -m pytest tests/test_espIntegrals.py -q`
Expected: `25 passed, 8 skipped` (the skips are `TestNumba`).

- [x] **Step 5: Commit**

```bash
git add espIntegrals.py tests/test_espIntegrals.py
git commit -m "Add exact ESP integrals (McMurchie-Davidson)

Density contracted once into point-independent Hermite coefficients,
threaded over point chunks; agrees with Gaussian 09 cubegen to <1e-6 a.u.
on every fixture, tr(PS)=N to 1e-6.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: `molecularSurface.py` — isosurface and Multiwfn extrema

**Files:**
- Create: `tests/test_molecularSurface.py`
- Create: `molecularSurface.py`

**Interfaces:**
- Consumes: `gaussianWavefunction.{read_fchk, electron_density}` (tests only).
- Produces (used by Task 6):
  - `DEFAULT_SPACING_BOHR = 0.25`, `MERGE_RATIO = 0.5`, `N_BISECT = 3`
  - `grid_points(origin, spacing, shape) -> (N,3)`; `grid_box(coords, radii_bohr, density, isovalue, spacing) -> (origin, shape)`
  - `density_outward(density) -> callable`; `marching_tetrahedra(values, origin, spacing, iso, density=None, n_bisect=3) -> (verts, tris)`
  - `eliminate_redundant_vertices(verts, tris, merge_distance) -> (verts, tris, csr)`; `build_isosurface(values, origin, spacing, iso, density, outward) -> SurfaceMesh`
  - `SurfaceMesh(vertices, triangles, neighbours, isovalue, spacing)` with `n_edges`, `euler_characteristic()`, `component_labels()`, `component_euler_characteristics()`, `area()`, `vertex_normals()`
  - `SurfaceExtremum(vertex, kind, scope, value, position)`; `find_surface_extrema(values, vertices, neighbours) -> List[SurfaceExtremum]` ordered: minima (global first, then by value ascending), then maxima (global first, then descending); `second_ring(csr)`

- [x] **Step 1: Write the failing tests** — `tests/test_molecularSurface.py`:

```python
"""Tests for isosurface extraction and Multiwfn-style surface extrema."""
import numpy as np
import pytest
from scipy import sparse

import gaussianWavefunction as gw
import molecularSurface as ms
from esp_fixtures import fchk_path

SPACING = 0.25


def gaussian_blobs(centres):
    centres = np.atleast_2d(centres)

    def density(p):
        d2 = ((p[:, None, :] - centres[None, :, :]) ** 2).sum(axis=2)
        return np.exp(-d2).sum(axis=1)
    return density


def mesh_for(density, isovalue, lo, hi):
    origin = np.array(lo, dtype=float)
    shape = tuple(int(round((h - l) / SPACING)) + 1 for l, h in zip(lo, hi))
    values = density(ms.grid_points(origin, SPACING, shape)).reshape(shape)
    return ms.build_isosurface(values, origin, SPACING, isovalue, density, ms.density_outward(density))


@pytest.fixture(scope="module")
def sphere():
    radius = 3.0
    return radius, mesh_for(gaussian_blobs([0, 0, 0]), np.exp(-radius ** 2), (-5, -5, -5), (5, 5, 5))


class TestSphere:
    def test_closed_single_surface(self, sphere):
        _, mesh = sphere
        assert mesh.euler_characteristic() == 2
        assert mesh.component_euler_characteristics() == [2]

    def test_vertices_lie_on_the_sphere(self, sphere):
        radius, mesh = sphere
        assert np.abs(np.linalg.norm(mesh.vertices, axis=1) - radius).max() < SPACING / 10

    def test_area(self, sphere):
        radius, mesh = sphere
        assert mesh.area() == pytest.approx(4 * np.pi * radius ** 2, rel=0.01)

    def test_normals_point_outward(self, sphere):
        _, mesh = sphere
        n = mesh.vertex_normals()
        assert np.all(np.einsum("ij,ij->i", n, mesh.vertices) > 0)

    def test_redundant_vertices_are_removed(self, sphere):
        radius, mesh = sphere
        density = gaussian_blobs([0, 0, 0])
        shape = (41, 41, 41)
        origin = np.array([-5.0, -5.0, -5.0])
        values = density(ms.grid_points(origin, SPACING, shape)).reshape(shape)
        raw, _ = ms.marching_tetrahedra(values, origin, SPACING, np.exp(-radius ** 2), density)
        # Lu & Chen report > 60 % of the marching-tetrahedra vertices removed.
        # (No minimum spacing is guaranteed afterwards: Multiwfn moves the
        # survivor to the midpoint of each merge, which can create new close
        # pairs, and it does not iterate.)
        assert len(mesh.vertices) < 0.5 * len(raw)


class TestTopology:
    def test_two_separate_spheres(self):
        density = gaussian_blobs([[-4, 0, 0], [4, 0, 0]])
        mesh = mesh_for(density, np.exp(-4.0), (-7, -3, -3), (7, 3, 3))
        assert mesh.component_euler_characteristics() == [2, 2]

    def test_no_surface_inside_the_grid(self):
        density = gaussian_blobs([0, 0, 0])
        shape = (9, 9, 9)
        values = density(ms.grid_points(np.zeros(3) - 1, SPACING, shape)).reshape(shape)
        with pytest.raises(ValueError, match="No isosurface"):
            ms.marching_tetrahedra(values, np.zeros(3) - 1, SPACING, 5.0)


class TestMolecularSurface:
    @pytest.fixture(scope="class")
    def eteno(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))

        def density(p):
            return gw.electron_density(w, p)

        radii = np.array([1.7, 1.2, 1.2, 1.7, 1.2, 1.2]) / 0.529177210544
        origin, shape = ms.grid_box(w.coordinates, radii, density, 0.0004, SPACING)
        values = density(ms.grid_points(origin, SPACING, shape)).reshape(shape)
        mesh = ms.build_isosurface(values, origin, SPACING, 0.0004, density, ms.density_outward(density))
        return w, density, origin, shape, mesh

    def test_grid_box_faces_are_below_a_tenth_of_the_isovalue(self, eteno):
        _, density, origin, shape, _ = eteno
        values = density(ms.grid_points(origin, SPACING, shape)).reshape(shape)
        faces = [values[0], values[-1], values[:, 0], values[:, -1], values[:, :, 0], values[:, :, -1]]
        assert max(f.max() for f in faces) < 0.0004 / 10

    def test_closed_and_on_the_isovalue(self, eteno):
        _, density, _, _, mesh = eteno
        assert mesh.component_euler_characteristics() == [2]
        assert np.abs(density(mesh.vertices) - 0.0004).max() < 1e-5


def chain(values):
    """A path graph 0-1-2-...: the smallest mesh that tells 1-ring from 2-ring."""
    n = len(values)
    rows = list(range(n - 1)) + list(range(1, n))
    cols = list(range(1, n)) + list(range(n - 1))
    nbr = sparse.csr_matrix((np.ones(len(rows), bool), (rows, cols)), shape=(n, n))
    return np.asarray(values, float), np.zeros((n, 3)), nbr


class TestExtrema:
    def test_second_ring_is_required(self):
        # vertex 3 (value 2) beats its direct neighbours (3, 6) but not vertex 1
        # (value 1) two steps away, so Multiwfn does not call it a minimum.
        values, verts, nbr = chain([5, 1, 3, 2, 6])
        minima = [e for e in ms.find_surface_extrema(values, verts, nbr) if e.kind == "minimum"]
        assert [(e.vertex, e.scope) for e in minima] == [(1, "global")]

    def test_ties_give_no_local_extremum(self):
        values, verts, nbr = chain([2, 1, 1, 2, 3])
        minima = [e for e in ms.find_surface_extrema(values, verts, nbr) if e.kind == "minimum"]
        # neither 1 is strictly below the other; the global minimum is still reported
        assert len(minima) == 1 and minima[0].scope == "global" and minima[0].value == 1.0

    def test_sphere_poles(self, sphere):
        _, mesh = sphere
        ex = ms.find_surface_extrema(mesh.vertices[:, 2], mesh.vertices, mesh.neighbours)
        assert [(e.kind, e.scope) for e in ex] == [("minimum", "global"), ("maximum", "global")]
        assert ex[0].position[2] < -2.9 and ex[1].position[2] > 2.9

    def test_two_maxima_two_minima(self, sphere):
        _, mesh = sphere
        x, y = mesh.vertices[:, 0], mesh.vertices[:, 1]
        ex = ms.find_surface_extrema(x * x - y * y + 0.01 * x, mesh.vertices, mesh.neighbours)
        kinds = [(e.kind, e.scope) for e in ex]
        assert kinds.count(("maximum", "global")) == 1 and kinds.count(("maximum", "local")) == 1
        assert kinds.count(("minimum", "global")) == 1 and kinds.count(("minimum", "local")) == 1
        # global first, then locals in order of extremity
        maxima = [e for e in ex if e.kind == "maximum"]
        assert maxima[0].value >= maxima[1].value


def test_numerical_modules_do_not_import_toga():
    """The fchk/integral/surface modules stay GUI-free (project convention)."""
    import os
    import subprocess
    import sys
    code = ("import sys, gaussianWavefunction, espIntegrals, molecularSurface; "
            "sys.exit(1 if 'toga' in sys.modules else 0)")
    src = os.path.dirname(os.path.abspath(ms.__file__))
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([src] + sys.path))
    assert subprocess.run([sys.executable, "-c", code], env=env, cwd=src).returncode == 0
```

- [x] **Step 2: Run them to confirm they fail**

Run: `..\Scripts\python.exe -m pytest tests/test_molecularSurface.py -q`
Expected: `ModuleNotFoundError: No module named 'molecularSurface'`.

- [x] **Step 3: Implement** — `molecularSurface.py`:

```python
"""Isosurface extraction (marching tetrahedra) and surface extrema.

Follows Multiwfn's quantitative molecular surface analysis (surfana.f90):
exact-function bisection of cut edges, merging of vertices closer than
0.5 x grid spacing, elimination of 2- and 3-connected vertices, and the
strict first+second-neighbour extremum test. Pure numpy/scipy, no Toga.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Set, Tuple

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

DEFAULT_SPACING_BOHR = 0.25
MERGE_RATIO = 0.5
N_BISECT = 3

# Kuhn (Freudenthal) decomposition of a cube into 6 tetrahedra sharing the
# main diagonal 0-7. Corner c has offset (c & 1, (c >> 1) & 1, (c >> 2) & 1).
_CORNER_OFFSETS = np.array([[c & 1, (c >> 1) & 1, (c >> 2) & 1] for c in range(8)])
_TETRAHEDRA = np.array([[0, 7, 1, 3], [0, 7, 3, 2], [0, 7, 2, 6],
                        [0, 7, 6, 4], [0, 7, 4, 5], [0, 7, 5, 1]])
# Every tetrahedron edge joins corners whose offsets differ by one of these.
_EDGE_TYPES = {(1, 0, 0): 0, (0, 1, 0): 1, (0, 0, 1): 2, (1, 1, 0): 3,
               (1, 0, 1): 4, (0, 1, 1): 5, (1, 1, 1): 6}
_EDGE_LUT = np.full(8, -1)
for _off, _k in _EDGE_TYPES.items():
    _EDGE_LUT[_off[0] + 2 * _off[1] + 4 * _off[2]] = _k


def density_outward(density: Callable[[np.ndarray], np.ndarray], eps: float = 1e-3):
    """-grad(rho) by central differences: the outward direction of a rho isosurface."""
    def outward(points: np.ndarray) -> np.ndarray:
        g = np.empty_like(points)
        for d in range(3):
            step = np.zeros(3); step[d] = eps
            g[:, d] = density(points - step) - density(points + step)
        return g
    return outward


@dataclass
class SurfaceMesh:
    vertices: np.ndarray      # (V, 3) bohr
    triangles: np.ndarray     # (F, 3) indices into vertices, outward oriented
    neighbours: sparse.csr_matrix  # (V, V) boolean vertex adjacency
    isovalue: float
    spacing: float

    @property
    def n_edges(self) -> int:
        return int(self.neighbours.nnz // 2)

    def euler_characteristic(self) -> int:
        return len(self.vertices) + len(self.triangles) - self.n_edges

    def component_labels(self) -> np.ndarray:
        """Connected-component label of every vertex (0-based)."""
        return connected_components(self.neighbours, directed=False)[1]

    def component_euler_characteristics(self) -> List[int]:
        """V + F - E of each connected component; 2 for every closed sphere-like surface."""
        labels = self.component_labels()
        n = int(labels.max()) + 1 if labels.size else 0
        v = np.bincount(labels, minlength=n)
        f = np.bincount(labels[self.triangles[:, 0]], minlength=n) if len(self.triangles) else np.zeros(n, int)
        rows, cols = self.neighbours.nonzero()
        e = np.bincount(labels[rows[rows < cols]], minlength=n)
        return [int(x) for x in v + f - e]

    def area(self) -> float:
        a, b, c = (self.vertices[self.triangles[:, k]] for k in range(3))
        return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())

    def vertex_normals(self) -> np.ndarray:
        a, b, c = (self.vertices[self.triangles[:, k]] for k in range(3))
        fn = np.cross(b - a, c - a)
        n = np.zeros_like(self.vertices)
        for k in range(3):
            np.add.at(n, self.triangles[:, k], fn)
        length = np.linalg.norm(n, axis=1)
        return n / np.where(length > 0, length, 1.0)[:, None]


def grid_points(origin: np.ndarray, spacing: float, shape: Tuple[int, int, int]) -> np.ndarray:
    axes = [origin[d] + spacing * np.arange(shape[d]) for d in range(3)]
    gx, gy, gz = np.meshgrid(*axes, indexing="ij")
    return np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])


def grid_box(coordinates: np.ndarray, radii: np.ndarray, density: Callable[[np.ndarray], np.ndarray],
             isovalue: float, spacing: float, max_rounds: int = 30):
    """Origin and shape of a cubic grid whose every face has rho < isovalue / 10.

    Starts from atoms +/- 1.7 x van der Waals radius (Multiwfn's vdwmulti) and
    pushes any face whose maximum density is still too high out by 1 bohr.
    """
    lo = (coordinates - 1.7 * radii[:, None]).min(axis=0)
    hi = (coordinates + 1.7 * radii[:, None]).max(axis=0)
    limit = isovalue / 10.0
    for _ in range(max_rounds):
        shape = tuple(int(np.ceil((hi[d] - lo[d]) / spacing)) + 1 for d in range(3))
        grown = False
        for d in range(3):
            for side in (0, 1):
                face_shape = list(shape)
                face_shape[d] = 1
                face_origin = lo.copy()
                if side:
                    face_origin[d] = lo[d] + (shape[d] - 1) * spacing
                pts = grid_points(face_origin, spacing, tuple(face_shape))
                if density(pts).max() >= limit:
                    grown = True
                    if side:
                        hi[d] += 1.0
                    else:
                        lo[d] -= 1.0
        if not grown:
            return lo, shape
    raise ValueError("Could not enclose the isosurface: density stays above "
                     f"{limit:.3g} a.u. on the grid boundary.")


def marching_tetrahedra(values: np.ndarray, origin: np.ndarray, spacing: float, isovalue: float,
                        density: Callable[[np.ndarray], np.ndarray] | None = None,
                        n_bisect: int = N_BISECT):
    """Raw triangulation of {rho = isovalue}: (vertices, triangles).

    Vertices are shared between triangles (one per cut grid edge), so the
    result is a connected mesh. With ``density`` the cut point on each edge
    is refined by ``n_bisect`` bisections of the exact function before the
    final linear interpolation (Multiwfn nbisec=3).
    """
    nx, ny, nz = values.shape
    inside = values >= isovalue
    cube = np.stack([inside[ox:nx - 1 + ox, oy:ny - 1 + oy, oz:nz - 1 + oz]
                     for ox, oy, oz in _CORNER_OFFSETS], axis=-1)
    boundary = np.argwhere(cube.any(axis=-1) & ~cube.all(axis=-1))
    if boundary.size == 0:
        raise ValueError(f"No isosurface at rho = {isovalue:g} a.u. inside the grid.")
    corner_in = cube[boundary[:, 0], boundary[:, 1], boundary[:, 2]]  # (C, 8)
    corner_idx = boundary[:, None, :] + _CORNER_OFFSETS[None, :, :]   # (C, 8, 3)

    def edge_id(ga, gb):
        lo_ = np.minimum(ga, gb)
        k = _EDGE_LUT[np.abs(gb - ga) @ np.array([1, 2, 4])]
        lin = (lo_[..., 0] * ny + lo_[..., 1]) * nz + lo_[..., 2]
        return lin * 7 + k

    tri_edges, tri_dirs = [], []
    for tet in _TETRAHEDRA:
        ins = corner_in[:, tet]                       # (C, 4)
        code = ins.astype(int) @ np.array([1, 2, 4, 8])
        g = corner_idx[:, tet, :]                     # (C, 4, 3)
        n_in = ins.sum(axis=1)
        for count in (1, 2, 3):
            sel = np.nonzero(n_in == count)[0]
            if sel.size == 0:
                continue
            for c in np.unique(code[sel]):
                rows = sel[code[sel] == c]
                bits = [(c >> k) & 1 for k in range(4)]
                ins_k = [k for k in range(4) if bits[k]]
                out_k = [k for k in range(4) if not bits[k]]
                pairs = [(i, o) for i in ins_k for o in out_k]
                if len(pairs) == 3:
                    polys = [pairs]
                else:  # 2 in, 2 out: quad (i0o0, i0o1, i1o1, i1o0) -> 2 triangles
                    (i0, i1), (o0, o1) = ins_k, out_k
                    quad = [(i0, o0), (i0, o1), (i1, o1), (i1, o0)]
                    polys = [[quad[0], quad[1], quad[2]], [quad[0], quad[2], quad[3]]]
                for poly in polys:
                    ids = np.stack([edge_id(g[rows, i], g[rows, o]) for i, o in poly], axis=1)
                    dirs = np.stack([g[rows, o] - g[rows, i] for i, o in poly], axis=1)
                    tri_edges.append(ids)
                    tri_dirs.append(dirs.sum(axis=1))
    tri_edges = np.concatenate(tri_edges)
    outward = np.concatenate(tri_dirs).astype(float)
    unique_ids, triangles = np.unique(tri_edges, return_inverse=True)
    triangles = triangles.reshape(-1, 3)

    # Decode each unique edge back to its two grid points (a inside, b outside).
    kinds = np.array(list(_EDGE_TYPES.keys()))
    lin, k = np.divmod(unique_ids, 7)
    ga = np.column_stack([lin // (ny * nz), (lin // nz) % ny, lin % nz])
    gb = ga + kinds[k]
    va = values[ga[:, 0], ga[:, 1], ga[:, 2]]
    vb = values[gb[:, 0], gb[:, 1], gb[:, 2]]
    pa = origin + spacing * ga
    pb = origin + spacing * gb
    if density is not None:
        for _ in range(n_bisect):
            mid = 0.5 * (pa + pb)
            vm = density(mid)
            same_as_a = (vm - isovalue) * (va - isovalue) > 0
            pa = np.where(same_as_a[:, None], mid, pa)
            va = np.where(same_as_a, vm, va)
            pb = np.where(same_as_a[:, None], pb, mid)
            vb = np.where(same_as_a, vb, vm)
    scl = (va - isovalue) / (va - vb)
    vertices = pa + scl[:, None] * (pb - pa)

    a, b, c = (vertices[triangles[:, k]] for k in range(3))
    flip = np.einsum("ij,ij->i", np.cross(b - a, c - a), outward) < 0
    triangles[flip] = triangles[flip][:, [0, 2, 1]]
    return vertices, triangles


def _adjacency_sets(n: int, triangles: np.ndarray) -> List[Set[int]]:
    adj: List[Set[int]] = [set() for _ in range(n)]
    for i, j, k in triangles.tolist():
        adj[i].update((j, k)); adj[j].update((i, k)); adj[k].update((i, j))
    return adj


def eliminate_redundant_vertices(vertices: np.ndarray, triangles: np.ndarray, merge_distance: float):
    """Multiwfn's vertex merge + 2-/3-neighbour elimination (surfana.f90:654-780).

    Returns the surviving vertices, re-indexed triangles and adjacency sets.
    """
    v = vertices.copy()
    n = len(v)
    adj = _adjacency_sets(n, triangles)
    dead = np.zeros(n, dtype=bool)
    parent = np.arange(n)
    crit2 = merge_distance ** 2
    for i in range(n):
        if dead[i]:
            continue
        for j in list(adj[i]):
            if dead[j] or j == i:
                continue
            d = v[i] - v[j]
            if d @ d < crit2:
                for k in adj[j]:
                    if dead[k] or k == i:
                        continue
                    if i not in adj[k]:
                        adj[k].add(i)
                        adj[i].add(k)
                parent[j] = i
                dead[j] = True
                v[i] = 0.5 * (v[i] + v[j])

    def root(x):
        while dead[x] and parent[x] != x:
            x = parent[x]
        return x

    tris = [[root(x) for x in t] for t in triangles.tolist()]
    tris = [t for t in tris if len(set(t)) == 3]

    def live(i):
        return [j for j in adj[i] if not dead[j]]

    for i in range(n):  # vertices with exactly two live neighbours
        if not dead[i] and len(live(i)) == 2:
            dead[i] = True
    for i in range(n):  # vertices with exactly three: drop and re-cap
        if not dead[i]:
            nb = live(i)
            if len(nb) == 3:
                dead[i] = True
                tris.append(nb)
    tris = np.array([t for t in tris if not any(dead[x] for x in t)], dtype=np.int64).reshape(-1, 3)
    keep = np.nonzero(~dead)[0]
    new_index = -np.ones(n, dtype=np.int64)
    new_index[keep] = np.arange(keep.size)
    rows, cols = [], []
    for i in keep:
        for j in adj[i]:
            if not dead[j]:
                rows.append(new_index[i]); cols.append(new_index[j])
    m = keep.size
    nbr = sparse.csr_matrix((np.ones(len(rows), dtype=bool), (rows, cols)), shape=(m, m))
    return v[keep], new_index[tris], nbr


def orient_triangles(vertices, triangles, density_gradient_sign):
    """Make every triangle's normal point along ``density_gradient_sign`` (outward)."""
    a, b, c = (vertices[triangles[:, k]] for k in range(3))
    centre = (a + b + c) / 3.0
    normal = np.cross(b - a, c - a)
    flip = np.einsum("ij,ij->i", normal, density_gradient_sign(centre)) < 0
    out = triangles.copy()
    out[flip] = out[flip][:, [0, 2, 1]]
    return out


def build_isosurface(values, origin, spacing, isovalue, density, outward_direction,
                     merge_ratio: float = MERGE_RATIO, n_bisect: int = N_BISECT) -> SurfaceMesh:
    verts, tris = marching_tetrahedra(values, origin, spacing, isovalue, density, n_bisect)
    verts, tris, nbr = eliminate_redundant_vertices(verts, tris, merge_ratio * spacing)
    tris = orient_triangles(verts, tris, outward_direction)
    return SurfaceMesh(verts, tris, nbr, isovalue, spacing)


# ----------------------------------------------------------------------
# Extrema
# ----------------------------------------------------------------------
@dataclass
class SurfaceExtremum:
    vertex: int
    kind: str      # "minimum" | "maximum"
    scope: str     # "global" | "local"
    value: float
    position: np.ndarray  # bohr


def second_ring(neighbours: sparse.csr_matrix) -> sparse.csr_matrix:
    """First + second neighbours of every vertex, the vertex itself excluded."""
    a = neighbours.astype(np.int32)
    ring = ((a + a @ a) > 0).tolil()
    ring.setdiag(False)
    return ring.tocsr()


def find_surface_extrema(values: np.ndarray, vertices: np.ndarray,
                         neighbours: sparse.csr_matrix) -> List[SurfaceExtremum]:
    """Multiwfn's rule: a minimum is strictly below every 1st- and 2nd-ring
    neighbour (maximum: strictly above). Ties give no extremum. The global
    minimum/maximum are argmin/argmax over all vertices and are reported as
    "global" even when a tie keeps them out of the strict local test."""
    ring = second_ring(neighbours)
    has = np.diff(ring.indptr) > 0
    vals = values[ring.indices]
    starts = ring.indptr[:-1][has]
    lo = np.full(values.size, np.inf)
    hi = np.full(values.size, -np.inf)
    lo[has] = np.minimum.reduceat(vals, starts)
    hi[has] = np.maximum.reduceat(vals, starts)
    gmin, gmax = int(np.argmin(values)), int(np.argmax(values))
    out: List[SurfaceExtremum] = []
    for kind, mask, g in (("minimum", has & (values < lo), gmin), ("maximum", has & (values > hi), gmax)):
        idx = set(np.nonzero(mask)[0].tolist()) | {g}
        order = sorted(idx, key=lambda i: values[i], reverse=(kind == "maximum"))
        for i in order:
            out.append(SurfaceExtremum(i, kind, "global" if i == g else "local",
                                       float(values[i]), vertices[i].copy()))
    return out
```

- [x] **Step 4: Run the tests**

Run: `..\Scripts\python.exe -m pytest tests/test_molecularSurface.py -q`
Expected: `14 passed`.

- [x] **Step 5: Commit**

```bash
git add molecularSurface.py tests/test_molecularSurface.py
git commit -m "Add marching-tetrahedra isosurface and Multiwfn surface extrema

Exact-density bisection, 0.5h vertex merge and 2/3-neighbour elimination
as in Multiwfn surfana.f90; strict 1st+2nd-ring extremum rule; adaptive
grid box for low isovalues.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Viewer hooks, vdW radii/help text, and `espViewer.py`

**Files:**
- Modify: `molecularViewer.py` (4 small edits, below)
- Modify: `help.py` (vdW table after `max_bond_order`; `help_esp_surface` at the end of `HelpGqteaWin`)
- Create: `tests/test_espViewer.py`
- Create: `espViewer.py`

**Interfaces:**
- Consumes: `espIntegrals.HARTREE_TO_KCAL_MOL`, `molecularViewer.MolecularViewer`.
- Produces (used by Task 6):
  - `MolecularViewer.window_title`, `MolecularViewer.render_scene_extras()`, `MolecularViewer.draw_screen_overlay()` (no-ops)
  - `AtomicData.vdw_radii` (Å), `AtomicData.DEFAULT_VDW_RADIUS = 2.00`, `HelpGqteaWin.help_esp_surface`
  - `espViewer`: `UNIT_AU`, `UNIT_KCAL`, `COLOUR_UNITS`, `hue_colours`, `symmetric_range`, `esp_colours`, `gaussview_label`, `unit_factor`, `colour_bar_text`, `SurfaceMarker(position, kind, scope, label)`, `ESPSurfaceViewer` with `set_molecule(symbols, coords_A)`, `set_surface(vertices_A, triangles, normals, values_au, isovalue, markers=(), half_range=None)`, `clear_surface()`, `set_display_options(unit=None, opacity=None, show_extrema=None)`, `get_surface()`, `colour_bar_text()`, `request_close()`

- [x] **Step 1: Write the failing tests** — `tests/test_espViewer.py`:

```python
"""Headless tests for the ESP viewer: colour maths and surface state (no GL window)."""
import numpy as np
import pytest

import espViewer as ev
from molecularViewer import MolecularViewer


class TestColours:
    def test_hue_sweep_end_points(self):
        rgb = ev.hue_colours(np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
        assert np.allclose(rgb, [[1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1]])

    def test_esp_zero_is_green_and_clipped(self):
        rgb = ev.esp_colours(np.array([-0.05, 0.0, 0.05, 0.2]), 0.05)
        assert np.allclose(rgb[0], [1, 0, 0]) and np.allclose(rgb[1], [0, 1, 0])
        assert np.allclose(rgb[2], [0, 0, 1]) and np.allclose(rgb[3], [0, 0, 1])

    def test_symmetric_range_is_max_abs(self):
        assert ev.symmetric_range(np.array([-0.02097, 0.01166])) == pytest.approx(0.02097)
        assert ev.symmetric_range(np.zeros(3)) == 1.0  # never a zero-width bar

    @pytest.mark.parametrize("value, text", [(-0.02072, "-2.072e-2"), (0.02072, "2.072e-2"),
                                             (-13.157, "-1.316e1"), (0.0, "0.000e0")])
    def test_gaussview_label(self, value, text):
        assert ev.gaussview_label(value) == text

    def test_colour_bar_text_names_both_quantities(self):
        left, right, caption = ev.colour_bar_text(0.02097, ev.UNIT_KCAL, 0.0004)
        assert left == "-1.316e1" and right == "1.316e1"
        assert "ESP (kcal/mol)" in caption and "rho = 0.000400 a.u." in caption

    def test_unknown_unit(self):
        with pytest.raises(ValueError):
            ev.unit_factor("eV")


class TestViewerState:
    @pytest.fixture
    def viewer(self):
        v = ev.ESPSurfaceViewer()
        v.set_molecule(["C", "C"], np.array([[0, 0, -0.67], [0, 0, 0.67]]))
        return v

    def test_base_viewer_hooks_are_no_ops(self):
        base = MolecularViewer()
        assert base.render_scene_extras() is None and base.draw_screen_overlay() is None
        assert base.window_title == "3D Molecular Viewer"

    def test_set_surface_builds_draw_data(self, viewer):
        verts = np.array([[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]])
        viewer.set_surface(verts, np.array([[0, 1, 2]]), np.eye(3), np.array([-0.02, 0.0, 0.01]), 0.0004)
        s = viewer.get_surface()
        assert s.vertices.dtype == np.float32 and s.indices.dtype == np.uint32
        assert s.half_range == pytest.approx(0.02)
        assert np.allclose(s.rgba[:, 3], 1.0)
        assert viewer.get_scene_radius() >= 3.0  # the camera frames the surface, not just the atoms
        assert viewer.colour_bar_text()[0] == "-2.000e-2"

    def test_opacity_updates_alpha(self, viewer):
        viewer.set_surface(np.zeros((3, 3)) + [[1, 0, 0]], np.array([[0, 1, 2]]), np.eye(3), np.zeros(3), 0.001)
        viewer.set_display_options(opacity=0.5)
        assert np.allclose(viewer.get_surface().rgba[:, 3], 0.5)

    def test_clear_surface(self, viewer):
        viewer.set_surface(np.eye(3), np.array([[0, 1, 2]]), np.eye(3), np.zeros(3), 0.001)
        viewer.clear_surface()
        assert viewer.get_surface() is None and viewer.colour_bar_text() is None
```

- [x] **Step 2: Run them to confirm they fail**

Run: `..\Scripts\python.exe -m pytest tests/test_espViewer.py -q`
Expected: `ModuleNotFoundError: No module named 'espViewer'`.

- [x] **Step 3: Add the viewer hooks** — four edits in `molecularViewer.py`:

(a) at the end of `MolecularViewer.__init__`, after `self.loop = None`:

```python
        # UI-loop pointer (set by subclass)
        self.loop = None

        # GLFW window title; subclasses that reuse the viewer set their own.
        self.window_title = "3D Molecular Viewer"
```

(b) in `init_glfw`, replace the literal title:

```python
        self.window = glfw.create_window(
            self.display[0], self.display[1], self.window_title, None, None
        )
```

(c) immediately **before** `def main_loop(self):` add:

```python
    def render_scene_extras(self):
        """Hook: extra geometry in the molecule's own frame (Angstrom).

        Called by main_loop on the render thread, inside the same matrix the
        atoms are drawn in, after the atoms, bonds and box. No-op here;
        espViewer.ESPSurfaceViewer draws its mapped surface from it.
        """

    def draw_screen_overlay(self):
        """Hook: 2D overlay drawn after the scene (e.g. a colour bar).

        Called on the render thread once per frame, after the axis triad.
        No-op here.
        """

```

(d) in `main_loop`, call the hooks. The two lines marked `# + new` are the only additions; do **not** copy the `# + new` markers themselves:

```python
                    self.render_frame(frame_data, bonds)
                    if self.show_box and not self.fast_playback_mode:
                        self.draw_box(frame_data)
                    self.render_scene_extras()          # + new
                    glPopMatrix()

                if not self.fast_playback_mode:
                    self.draw_cartesian_axes_overlay()
                self.draw_screen_overlay()              # + new
```

- [x] **Step 4: Extend `help.py`**

(a) In `class AtomicData`, directly after the closing `}` of `max_bond_order`:

```python
    # Van der Waals radii in Angstrom: Bondi, J. Phys. Chem. 1964, 68, 441,
    # completed for main-group elements by Mantina et al., J. Phys. Chem. A
    # 2009, 113, 5806. The ESP surface tool uses them only to start its grid
    # box (atoms +/- 1.7 r_vdW, as Multiwfn does); the box then grows until
    # the density on every face is below a tenth of the isovalue, so a rough
    # or missing radius costs time, never accuracy. Absent elements use
    # DEFAULT_VDW_RADIUS.
    vdw_radii = {
        "H": 1.20, "He": 1.40,
        "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52,
        "F": 1.47, "Ne": 1.54,
        "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80,
        "Cl": 1.75, "Ar": 1.88,
        "K": 2.75, "Ca": 2.31, "Ni": 1.63, "Cu": 1.40, "Zn": 1.39,
        "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.85, "Kr": 2.02,
        "Rb": 3.03, "Sr": 2.49, "Pd": 1.63, "Ag": 1.72, "Cd": 1.58,
        "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16,
        "Cs": 3.43, "Ba": 2.68, "Pt": 1.72, "Au": 1.66, "Hg": 1.55,
        "Tl": 1.96, "Pb": 2.02, "Bi": 2.07, "Po": 1.97, "At": 2.02, "Rn": 2.20,
        "Fr": 3.48, "Ra": 2.83, "U": 1.86,
    }
    DEFAULT_VDW_RADIUS = 2.00
```

(b) At the end of `class HelpGqteaWin` (after `help_molecular_design`), add:

```python
    help_esp_surface = """ELECTROSTATIC SURFACE POTENTIAL (ESP)

Maps the molecular electrostatic potential V(r) onto an electron-density
isosurface, shows it colour-coded in a 3D window and lists the surface
maxima and minima - the quantitative molecular surface analysis of
Multiwfn (Lu & Chen, J. Mol. Graph. Model. 2012, 38, 314), drawn the way
GaussView 6 draws a mapped ESP.

INPUT
  A Gaussian formatted checkpoint, .fchk or .fch (formchk name.chk name.fchk).
  The 'Total SCF Density' stored in it is the density used, for both the
  surface and the potential. Open-shell files work: the total density is
  what the ESP needs.

TWO DIFFERENT QUANTITIES - do not confuse them
  Surface definition : electron density rho(r), in a.u. (e/bohr^3). The
                       isovalue picks WHERE the surface is. Default 0.000400
                       (GaussView's default); 0.001 is the classic
                       Bader van der Waals surface.
  Mapped property    : electrostatic potential V(r), in a.u. (Hartree/e)
                       or kcal/mol. It is evaluated ON that surface and
                       is what the colours and the extrema report.

SETUP
  Isovalue(s)     One or several, separated by commas or spaces
                  (e.g. 0.0004, 0.001, 0.002). Blank = 0.000400.
  Grid spacing    Spacing of the density grid the surface is built from,
                  in bohr. Blank = 0.25 (Multiwfn's default for the ESP).
                  Smaller is smoother and slower.
  ESP engine      Built-in exact: analytic integrals in this program.
                  Gaussian cubegen: Gaussian's own program, if installed
                  (give its path). If cubegen fails the built-in engine is
                  used and the summary says so.

HOW THE SURFACE AND THE EXTREMA ARE FOUND (Multiwfn's procedure)
  1. rho is computed on a cubic grid that encloses the whole surface.
  2. The isosurface is triangulated by marching tetrahedra; each surface
     point is refined by 3 bisections of the exact density.
  3. Surface points closer than half the grid spacing are merged, and
     points joined to only 2 or 3 neighbours are removed.
  4. V(r) is computed exactly at every remaining surface point.
  5. A point is a LOCAL MINIMUM when its V is strictly lower than at every
     neighbouring point AND every neighbour of those neighbours (two rings);
     a LOCAL MAXIMUM likewise, strictly higher. Equal values give no
     extremum, and nearby extrema are never merged automatically.
     The GLOBAL minimum/maximum is the lowest/highest point of the surface.
  Extremum positions are surface-point positions (not interpolated).

OUTPUT  (next to the .fchk unless another folder is chosen)
  <name>_ESP_extrema.txt  every extremum of every surface: rho isovalue,
  minimum/maximum, global/local, V in a.u. and kcal/mol, x y z in Angstrom
  and the nearest atom. 1 a.u. = 627.5094740631 kcal/mol (CODATA 2018);
  Multiwfn prints with 627.51, a difference below 1e-6.

3D WINDOW
  Colour bar across the top, GaussView style: red = most negative,
  green = zero, blue = most positive, over a symmetric range +/- max|V|
  of the shown surface. Rotate with the left mouse button, zoom with the
  wheel. Markers: cyan = minimum, orange = maximum (larger = global); the
  number is the row in the extrema file.
"""
```

- [x] **Step 5: Implement** — `espViewer.py`:

```python
"""3D window for an ESP-mapped molecular surface, GaussView 6 style.

``ESPSurfaceViewer`` subclasses the project's GLFW/OpenGL MolecularViewer,
so camera, rotation, zoom, atoms, bonds and click-to-label all come from
there; it adds the colour-mapped surface (``render_scene_extras``) and the
colour bar (``draw_screen_overlay``).

GaussView conventions reproduced (from a GaussView 6.1 screenshot of the
ethylene ESP at rho = 0.0004): horizontal colour bar across the top, a
symmetric range +/-R with the min label at the left end and the max label at
the right end written like ``-2.072e-2``, colours a linear hue sweep red ->
orange -> yellow -> green (V = 0) -> cyan -> blue, opaque surface with a
specular highlight on a black background. Unlike GaussView, the bar also
prints the unit and the density isovalue that defines the surface.

The colour maths are module-level pure functions, tested headless.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import glfw
import numpy as np
from OpenGL.GL import *  # noqa: F401,F403 - same style as molecularViewer

from espIntegrals import HARTREE_TO_KCAL_MOL
from molecularViewer import MolecularViewer

UNIT_AU = "a.u."
UNIT_KCAL = "kcal/mol"
COLOUR_UNITS = (UNIT_AU, UNIT_KCAL)
HUE_SPAN_DEGREES = 240.0

MINIMUM_MARKER_COLOUR = (0.0, 1.0, 1.0)
MAXIMUM_MARKER_COLOUR = (1.0, 0.55, 0.0)


def hue_colours(fraction: np.ndarray) -> np.ndarray:
    """RGB (N, 3) for fractions in [0, 1]: 0 = red, 0.5 = green, 1 = blue.

    HSV with full saturation and value, hue = 240 degrees x fraction.
    """
    h = np.clip(np.asarray(fraction, dtype=float), 0.0, 1.0) * HUE_SPAN_DEGREES / 60.0
    rgb = np.zeros(h.shape + (3,))
    seg = np.minimum(h.astype(int), 3)
    x = h - seg
    rgb[seg == 0] = np.column_stack([np.ones_like(x), x, np.zeros_like(x)])[seg == 0]
    rgb[seg == 1] = np.column_stack([1 - x, np.ones_like(x), np.zeros_like(x)])[seg == 1]
    rgb[seg == 2] = np.column_stack([np.zeros_like(x), np.ones_like(x), x])[seg == 2]
    rgb[seg == 3] = np.column_stack([np.zeros_like(x), 1 - x, np.ones_like(x)])[seg == 3]
    return rgb


def symmetric_range(values: np.ndarray) -> float:
    """GaussView's default half-range: max |V| over the surface."""
    r = float(np.max(np.abs(values))) if np.size(values) else 0.0
    return r if r > 0 else 1.0


def esp_colours(values: np.ndarray, half_range: float) -> np.ndarray:
    """RGB for ESP values over the symmetric range [-R, +R] (clipped)."""
    return hue_colours((np.asarray(values, dtype=float) + half_range) / (2.0 * half_range))


def gaussview_label(value: float) -> str:
    """'-2.072e-2' style: 4 significant digits, exponent without padding."""
    mantissa, exponent = f"{value:.3e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def unit_factor(unit: str) -> float:
    if unit == UNIT_AU:
        return 1.0
    if unit == UNIT_KCAL:
        return HARTREE_TO_KCAL_MOL
    raise ValueError(f"Unknown colour-bar unit '{unit}'.")


def colour_bar_text(half_range_au: float, unit: str, isovalue: float) -> Tuple[str, str, str]:
    """(left label, right label, caption) for the colour bar."""
    factor = unit_factor(unit)
    caption = f"ESP ({unit})   mapped on electron density isosurface  rho = {isovalue:.6f} a.u."
    return gaussview_label(-half_range_au * factor), gaussview_label(half_range_au * factor), caption


@dataclass
class SurfaceMarker:
    position: np.ndarray  # Angstrom
    kind: str             # "minimum" | "maximum"
    scope: str            # "global" | "local"
    label: str


@dataclass
class SurfaceDrawData:
    vertices: np.ndarray   # (V, 3) float32, Angstrom
    normals: np.ndarray    # (V, 3) float32
    indices: np.ndarray    # (3F,) uint32
    values: np.ndarray     # (V,) ESP in a.u.
    rgba: np.ndarray       # (V, 4) float32
    half_range: float      # a.u.
    isovalue: float
    markers: List[SurfaceMarker]
    radius: float          # max distance of a vertex from the atoms' centre


class ESPSurfaceViewer(MolecularViewer):
    """MolecularViewer that also draws one ESP-mapped surface and its colour bar."""

    def __init__(self):
        super().__init__()
        self.window_title = "ESP Surface Viewer"
        self.colour_unit = UNIT_AU
        self.surface_opacity = 1.0
        self.show_extrema = True
        self._surface: Optional[SurfaceDrawData] = None

    # --- state (Toga thread) -------------------------------------------
    def set_molecule(self, symbols: Sequence[str], coordinates_angstrom: np.ndarray):
        frame = [(s, tuple(float(c) for c in xyz)) for s, xyz in zip(symbols, coordinates_angstrom)]
        self.load_designed_molecule(frame)

    def set_surface(self, vertices_angstrom, triangles, normals, values_au, isovalue,
                    markers: Sequence[SurfaceMarker] = (), half_range: Optional[float] = None):
        values = np.asarray(values_au, dtype=float)
        r = half_range if half_range else symmetric_range(values)
        verts = np.ascontiguousarray(vertices_angstrom, dtype=np.float32)
        centre = np.array(self.get_system_center(), dtype=float)
        data = SurfaceDrawData(
            vertices=verts,
            normals=np.ascontiguousarray(normals, dtype=np.float32),
            indices=np.ascontiguousarray(np.asarray(triangles).ravel(), dtype=np.uint32),
            values=values,
            rgba=self._rgba(values, r),
            half_range=r,
            isovalue=float(isovalue),
            markers=list(markers),
            radius=float(np.linalg.norm(verts - centre, axis=1).max()) if len(verts) else 0.0,
        )
        with self._state_lock:
            self._surface = data
        self.invalidate_scene_cache()
        self.request_projection_update()

    def clear_surface(self):
        with self._state_lock:
            self._surface = None
        self.invalidate_scene_cache()

    def set_display_options(self, unit: str = None, opacity: float = None, show_extrema: bool = None):
        if unit is not None:
            unit_factor(unit)  # validates
            self.colour_unit = unit
        if show_extrema is not None:
            self.show_extrema = bool(show_extrema)
        if opacity is not None:
            self.surface_opacity = float(min(1.0, max(0.05, opacity)))
            with self._state_lock:
                if self._surface is not None:
                    self._surface.rgba = self._rgba(self._surface.values, self._surface.half_range)

    def get_surface(self) -> Optional[SurfaceDrawData]:
        with self._state_lock:
            return self._surface

    def _rgba(self, values, half_range) -> np.ndarray:
        rgb = esp_colours(values, half_range)
        alpha = np.full((len(rgb), 1), self.surface_opacity)
        return np.ascontiguousarray(np.hstack([rgb, alpha]), dtype=np.float32)

    def colour_bar_text(self) -> Optional[Tuple[str, str, str]]:
        s = self.get_surface()
        return None if s is None else colour_bar_text(s.half_range, self.colour_unit, s.isovalue)

    def get_scene_radius(self, frame_data=None) -> float:
        base = super().get_scene_radius(frame_data)
        s = self.get_surface()
        return max(base, s.radius + 0.3) if s is not None else base

    def request_close(self):
        """Ask GLFW to close the window (glfwSetWindowShouldClose is thread-safe)."""
        if self.window is not None:
            try:
                glfw.set_window_should_close(self.window, True)
            except Exception:
                pass

    # --- drawing (render thread) ----------------------------------------
    def render_scene_extras(self):
        s = self.get_surface()
        if s is None:
            return
        self._draw_surface(s)
        if self.show_extrema:
            self._draw_markers(s)

    def _draw_surface(self, s: SurfaceDrawData):
        translucent = self.surface_opacity < 0.999
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (0.6, 0.6, 0.6, 1.0))
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 40.0)
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
        if translucent:
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glDepthMask(GL_FALSE)
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_NORMAL_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        try:
            glVertexPointer(3, GL_FLOAT, 0, s.vertices)
            glNormalPointer(GL_FLOAT, 0, s.normals)
            glColorPointer(4, GL_FLOAT, 0, s.rgba)
            glDrawElements(GL_TRIANGLES, int(s.indices.size), GL_UNSIGNED_INT, s.indices)
        finally:
            glDisableClientState(GL_COLOR_ARRAY)
            glDisableClientState(GL_NORMAL_ARRAY)
            glDisableClientState(GL_VERTEX_ARRAY)
            if translucent:
                glDepthMask(GL_TRUE)
                glDisable(GL_BLEND)
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (0.0, 0.0, 0.0, 1.0))
            glDisable(GL_COLOR_MATERIAL)

    def _draw_markers(self, s: SurfaceDrawData):
        for m in s.markers:
            colour = MINIMUM_MARKER_COLOUR if m.kind == "minimum" else MAXIMUM_MARKER_COLOUR
            radius = 0.18 if m.scope == "global" else 0.12
            glPushMatrix()
            glTranslatef(*[float(c) for c in m.position])
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, (*colour, 1.0))
            self.draw_sphere(radius)
            glPopMatrix()
            label_pos = tuple(float(c) + 0.2 for c in m.position)
            self._draw_text_3d(label_pos, m.label, color=colour)

    def draw_screen_overlay(self):
        text = self.colour_bar_text()
        if text is None:
            return
        left, right, caption = text
        w, h = self.display
        x0, x1 = 100.0, max(110.0, w - 100.0)
        y0, y1 = h - 34.0, h - 16.0
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, w, 0, h, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        glDisable(GL_LIGHTING)
        glDisable(GL_DEPTH_TEST)
        try:
            steps = 96
            colours = hue_colours(np.linspace(0.0, 1.0, steps + 1))
            glBegin(GL_QUAD_STRIP)
            for k in range(steps + 1):
                x = x0 + (x1 - x0) * k / steps
                glColor3f(*colours[k])
                glVertex2f(x, y0)
                glVertex2f(x, y1)
            glEnd()
            self._draw_text_2d(6.0, y0 + 2.0, left)
            self._draw_text_2d(x1 + 6.0, y0 + 2.0, right)
            self._draw_text_2d(x0, y0 - 22.0, caption)
        finally:
            glEnable(GL_DEPTH_TEST)
            glEnable(GL_LIGHTING)
            glPopMatrix()
            glMatrixMode(GL_PROJECTION)
            glPopMatrix()
            glMatrixMode(GL_MODELVIEW)
```

- [x] **Step 6: Run the new and the affected existing tests**

Run: `..\Scripts\python.exe -m pytest tests/test_espViewer.py tests/test_molecularViewer.py tests/test_molecularDesign.py -q`
Expected: all pass (`13` new + `149` existing = `162 passed`).

- [x] **Step 7: Commit**

```bash
git add molecularViewer.py help.py espViewer.py tests/test_espViewer.py
git commit -m "Add ESP surface viewer on two new MolecularViewer hooks

render_scene_extras/draw_screen_overlay are no-ops on the base viewer.
ESPSurfaceViewer draws the colour-mapped surface, extrema markers and a
GaussView-style colour bar. Adds Bondi/Mantina vdW radii and ESP help.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: `electrostaticSurfacePotential.py` — analyser, extrema file, Toga UI

**Files:**
- Create: `tests/test_electrostaticSurfacePotential.py`
- Replace: `electrostaticSurfacePotential.py` (currently a 5-line stub, **untracked** in git)
- Commit as-is: `gqteaWinToga.py` — the user's own uncommitted edit (2 lines: the import and the `("Electrostatic Surface Potential (ESP)", ElectrostaticSurfacePotentialUI)` entry in `general_tools`). Do not change it; just include it in this task's commit so the tool is reachable from the committed branch.

**Interfaces:**
- Consumes: everything produced by Tasks 2–5; `coordinateConverter.BOHR_RADIUS_ANGSTROM`.
- Produces:
  - `parse_isovalues(text) -> List[float]`, `parse_spacing(text) -> float`, `resolve_output_dir(fchk, chosen) -> str`
  - `find_cubegen(explicit=None) -> Optional[str]`, `cubegen_potential(fchk, points_bohr, cubegen, timeout=3600) -> ndarray` (raises `RuntimeError`)
  - `ESPSurfaceAnalyser.run(fchk_path, isovalues=(0.0004,), spacing=0.25, engine=ENGINE_BUILTIN, cubegen_path=None, progress=None) -> ESPResults`; `.extrema_rows(results)`, `.format_extrema_file(results)`, `.write_extrema_file(results, output_dir) -> str` (raises `ValueError`), `.summary_text(results, extrema_path=None)`, `.nearest_atom_label(wfn, position)`
  - `SurfaceResult(isovalue, mesh, esp, extrema)` with `half_range`; `ESPResults(fchk_path, wavefunction, spacing, engine, surfaces, notes, elapsed)`
  - `ElectrostaticSurfacePotentialUI(*args, **kwargs)` with `layout_main_window`, async `workflow`, async `show_surface`, `clear_output`, `closeTopLevel`

- [x] **Step 1: Write the failing tests** — `tests/test_electrostaticSurfacePotential.py`:

```python
"""Tests for the ESP surface analyser, the extrema file and the UI plumbing."""
import asyncio
from types import SimpleNamespace

import numpy as np
import pytest

import electrostaticSurfacePotential as esp
from esp_fixtures import fchk_path

A0 = 0.529177210544


class TestParsing:
    def test_blank_isovalue_is_the_default(self):
        assert esp.parse_isovalues("") == [0.0004] and esp.parse_isovalues(None) == [0.0004]

    def test_several_isovalues_keep_order_without_duplicates(self):
        assert esp.parse_isovalues("0.001, 0.0004;0.002 0.001") == [0.001, 0.0004, 0.002]

    @pytest.mark.parametrize("text", ["abc", "0", "-0.001", "2.5"])
    def test_bad_isovalue(self, text):
        with pytest.raises(ValueError):
            esp.parse_isovalues(text)

    def test_esp_like_value_is_explained(self):
        with pytest.raises(ValueError, match="electron-density"):
            esp.parse_isovalues("1.5")

    def test_spacing(self):
        assert esp.parse_spacing("") == 0.25 and esp.parse_spacing("0.2") == 0.2
        with pytest.raises(ValueError):
            esp.parse_spacing("1.0")

    def test_output_dir_defaults_to_the_fchk_folder(self, tmp_path):
        f = tmp_path / "x.fchk"
        assert esp.resolve_output_dir(str(f), "") == str(tmp_path.resolve())
        assert esp.resolve_output_dir(str(f), "D:/out") == "D:/out"


@pytest.fixture(scope="module")
def eteno_results():
    analyser = esp.ESPSurfaceAnalyser()
    seen = []
    results = analyser.run(str(fchk_path("eteno_b3lyp_631g")), [0.0004, 0.001],
                           progress=lambda f, m: seen.append(f))
    return analyser, results, seen


class TestEteno:
    def test_progress_is_monotone_to_one(self, eteno_results):
        _, _, seen = eteno_results
        assert seen[-1] == 1.0 and all(b >= a - 1e-12 for a, b in zip(seen, seen[1:]))

    def test_two_pi_minima_and_four_h_maxima_per_surface(self, eteno_results):
        analyser, results, _ = eteno_results
        for surface in results.surfaces:
            kinds = [(e.kind, e.scope) for e in surface.extrema]
            assert kinds.count(("minimum", "global")) == 1 and kinds.count(("minimum", "local")) == 1
            assert kinds.count(("maximum", "global")) == 1 and kinds.count(("maximum", "local")) == 3
            assert surface.mesh.component_euler_characteristics() == [2]
        rows = analyser.extrema_rows(results)
        assert {r["nearest_atom"][0] for r in rows if r["kind"] == "minimum"} == {"C"}
        assert {r["nearest_atom"][0] for r in rows if r["kind"] == "maximum"} == {"H"}

    def test_values_at_rho_0_0004(self, eteno_results):
        """Spike reference (Multiwfn-style); GaussView 6.1 shows a +/-2.072e-2 bar for this file."""
        _, results, _ = eteno_results
        s = results.surfaces[0]
        assert s.isovalue == 0.0004
        assert s.esp.min() == pytest.approx(-0.02097, abs=2e-4)
        assert s.esp.max() == pytest.approx(0.01166, abs=2e-4)
        assert s.half_range == pytest.approx(0.0207, abs=5e-4)
        # the minima sit over the pi cloud: off the molecular (yz) plane along x
        minima = [e for e in s.extrema if e.kind == "minimum"]
        assert all(abs(e.position[0] * A0) > 2.0 for e in minima)

    def test_tighter_isovalue_gives_a_smaller_surface_and_deeper_minimum(self, eteno_results):
        _, results, _ = eteno_results
        outer, inner = results.surfaces
        assert inner.mesh.area() < outer.mesh.area()
        assert inner.esp.min() < outer.esp.min()


class TestExtremaFile:
    def test_format(self, eteno_results, tmp_path):
        analyser, results, _ = eteno_results
        path = analyser.write_extrema_file(results, str(tmp_path))
        assert path.endswith("eteno_b3lyp_631g_ESP_extrema.txt")
        text = open(path, encoding="utf-8").read()
        assert "# surface_definition     electron density isosurface" in text
        assert "# mapped_property        electrostatic potential" in text
        assert "# hartree_to_kcal_mol    627.5094740631" in text
        assert text.count("# surface rho_isovalue_au=") == 2
        data = [l.split() for l in text.splitlines() if not l.startswith("#")]
        assert len(data) == 12  # 6 extrema on each of the two surfaces
        for row in data:
            assert row[2] in ("minimum", "maximum") and row[3] in ("global", "local")
            assert float(row[5]) == pytest.approx(float(row[4]) * 627.5094740631, abs=1e-3)
        globals_ = [(r[1], r[2]) for r in data if r[3] == "global"]
        assert sorted(globals_) == sorted({(r[1], r[2]) for r in data})  # one global per type per surface

    def test_unwritable_folder_is_a_value_error(self, eteno_results, tmp_path):
        analyser, results, _ = eteno_results
        blocker = tmp_path / "file"
        blocker.write_text("x")
        with pytest.raises(ValueError, match="Could not write"):
            analyser.write_extrema_file(results, str(blocker / "sub"))

    def test_summary_names_the_surface_and_the_file(self, eteno_results):
        analyser, results, _ = eteno_results
        text = analyser.summary_text(results, "C:/x/eteno_ESP_extrema.txt")
        assert "rho = 0.000400 a.u. (electron density isovalue)" in text
        assert "Extrema written to: C:/x/eteno_ESP_extrema.txt" in text


class TestEngines:
    def test_missing_cubegen_falls_back_with_a_note(self, monkeypatch):
        monkeypatch.setattr(esp, "find_cubegen", lambda explicit=None: None)
        r = esp.ESPSurfaceAnalyser().run(str(fchk_path("water_631gd_6d")), [0.001],
                                         engine=esp.ENGINE_CUBEGEN)
        assert r.engine == esp.ENGINE_BUILTIN and "not found" in r.notes[0]

    def test_failing_cubegen_falls_back_with_a_note(self, monkeypatch):
        monkeypatch.setattr(esp, "find_cubegen", lambda explicit=None: "C:/fake/cubegen.exe")

        def broken(*args, **kwargs):
            raise RuntimeError("cubegen failed: test")
        monkeypatch.setattr(esp, "cubegen_potential", broken)
        r = esp.ESPSurfaceAnalyser().run(str(fchk_path("water_631gd_6d")), [0.001],
                                         engine=esp.ENGINE_CUBEGEN)
        assert r.engine == esp.ENGINE_BUILTIN and "cubegen failed" in r.notes[0]

    @pytest.mark.skipif(esp.find_cubegen() is None, reason="Gaussian cubegen not installed")
    def test_cubegen_engine_agrees_with_the_built_in_one(self):
        f = str(fchk_path("water_631gd_6d"))
        a = esp.ESPSurfaceAnalyser().run(f, [0.001], engine=esp.ENGINE_BUILTIN).surfaces[0]
        b = esp.ESPSurfaceAnalyser().run(f, [0.001], engine=esp.ENGINE_CUBEGEN).surfaces[0]
        assert np.max(np.abs(a.esp - b.esp)) < 1e-6

    def test_unknown_engine(self):
        with pytest.raises(ValueError):
            esp.ESPSurfaceAnalyser().run(str(fchk_path("water_631gd_6d")), engine="magic")


def make_ui():
    """The UI class without a window: stub widgets, as in test_meanResidenceTime."""
    ui = esp.ElectrostaticSurfacePotentialUI.__new__(esp.ElectrostaticSurfacePotentialUI)
    ui.analyser = esp.ESPSurfaceAnalyser()
    ui.results = None
    ui.extrema_path = None
    ui.viewer = None
    ui._syncing = False
    for name in ("fchk_input", "output_dir_input", "isovalue_input", "spacing_input", "cubegen_input"):
        setattr(ui, name, SimpleNamespace(value=""))
    ui.engine_selection = SimpleNamespace(value=esp.ENGINE_BUILTIN)
    ui.surface_selection = SimpleNamespace(value=None, items=[])
    ui.unit_selection = SimpleNamespace(value="a.u.")
    ui.opacity_selection = SimpleNamespace(value="Opaque")
    ui.extrema_switch = SimpleNamespace(value=True)
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


class TestUI:
    def test_missing_file_is_an_input_error(self):
        ui = make_ui()
        asyncio.run(ui.workflow(None))
        assert ui.errors and ui.errors[0][0] == "Input Error"

    def test_collect_inputs(self):
        ui = make_ui()
        ui.fchk_input.value = str(fchk_path("eteno_b3lyp_631g"))
        ui.isovalue_input.value = "0.0004 0.001"
        p = ui._collect_inputs()
        assert p["isovalues"] == [0.0004, 0.001] and p["spacing"] == 0.25 and p["cubegen_path"] is None

    def test_every_run_parameter_is_supplied_by_the_ui(self):
        import inspect
        ui = make_ui()
        ui.fchk_input.value = str(fchk_path("eteno_b3lyp_631g"))
        wanted = set(inspect.signature(esp.ESPSurfaceAnalyser.run).parameters) - {"self", "progress"}
        assert wanted == set(ui._collect_inputs())

    def test_workflow_writes_the_file_and_fills_the_selection(self, tmp_path):
        ui = make_ui()
        ui.fchk_input.value = str(fchk_path("water_631gd_6d"))
        ui.output_dir_input.value = str(tmp_path)
        ui.isovalue_input.value = "0.001"
        asyncio.run(ui.workflow(None))
        assert not ui.errors
        assert (tmp_path / "water_631gd_6d_ESP_extrema.txt").exists()
        assert ui.surface_selection.items == ["rho = 0.001000 a.u."]
        assert ui.selected_surface().isovalue == 0.001
        assert ui.tabs.current_tab == 1 and ui.run_button.enabled
        markers = ui.surface_markers(ui.selected_surface())
        assert markers and all(m.label.isdigit() for m in markers)

    def test_clear_resets_results(self, tmp_path):
        ui = make_ui()
        ui.fchk_input.value = str(fchk_path("water_631gd_6d"))
        ui.output_dir_input.value = str(tmp_path)
        asyncio.run(ui.workflow(None))
        ui.clear_output(None)
        assert ui.results is None and ui.surface_selection.items == []
        asyncio.run(ui.show_surface(None))
        assert ui.errors[-1][0] == "No Results"
```

- [x] **Step 2: Run them to confirm they fail**

Run: `..\Scripts\python.exe -m pytest tests/test_electrostaticSurfacePotential.py -q`
Expected: `AttributeError`/`ImportError` failures — the stub has no `parse_isovalues`, `ESPSurfaceAnalyser`, …

- [x] **Step 3: Implement** — replace the whole of `electrostaticSurfacePotential.py` with:

```python
"""Electrostatic surface potential (ESP) mapped on electron-density isosurfaces.

``ESPSurfaceAnalyser`` (headless) reads a Gaussian .fchk/.fch, builds one
isosurface of the Total SCF Density per requested isovalue, evaluates the
exact ESP at every surface vertex and finds the surface extrema with
Multiwfn's rules (see molecularSurface). ``ElectrostaticSurfacePotentialUI``
is the Toga window; ``espViewer.ESPSurfaceViewer`` draws the result.

Two quantities, never to be confused: the *surface* is defined by the
electron density rho (isovalue, a.u. = e/bohr^3); the *mapped property* is
the electrostatic potential V (a.u. = Hartree/e, or kcal/mol).
"""

from __future__ import annotations

import asyncio
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import espIntegrals as ei
import gaussianWavefunction as gw
import molecularSurface as ms
from coordinateConverter import BOHR_RADIUS_ANGSTROM
from espIntegrals import HARTREE_TO_KCAL_MOL
from espViewer import COLOUR_UNITS, ESPSurfaceViewer, SurfaceMarker, UNIT_AU, symmetric_range
from help import AtomicData, HelpGqteaWin

DEFAULT_ISOVALUE = 0.0004
MIN_SPACING, MAX_SPACING = 0.05, 0.6
ENGINE_BUILTIN = "Built-in exact"
ENGINE_CUBEGEN = "Gaussian cubegen"
ENGINES = (ENGINE_BUILTIN, ENGINE_CUBEGEN)
DEFAULT_CUBEGEN_LOCATIONS = (r"C:\G09W\cubegen.exe", r"C:\G16W\cubegen.exe")
EXTREMA_FILE_SUFFIX = "_ESP_extrema.txt"

ProgressCallback = Callable[[float, str], None]


# ----------------------------------------------------------------------
# Input parsing
# ----------------------------------------------------------------------
def parse_isovalues(text: Optional[str]) -> List[float]:
    """'0.0004, 0.001' -> [0.0004, 0.001]; blank -> [DEFAULT_ISOVALUE].

    Order is kept, duplicates dropped. Every value must be a positive
    number below 1 a.u. (a density isovalue, not an ESP value).
    """
    tokens = (text or "").replace(",", " ").replace(";", " ").split()
    if not tokens:
        return [DEFAULT_ISOVALUE]
    values: List[float] = []
    for tok in tokens:
        try:
            v = float(tok)
        except ValueError as exc:
            raise ValueError(f"Isovalue '{tok}' is not a number.") from exc
        if not (0.0 < v < 1.0):
            raise ValueError(
                f"Isovalue {tok} is outside (0, 1) a.u. It is an electron-density "
                "value in e/bohr^3 (e.g. 0.0004 or 0.001), not an ESP value."
            )
        if v not in values:
            values.append(v)
    return values


def parse_spacing(text: Optional[str]) -> float:
    """Grid spacing in bohr; blank -> molecularSurface.DEFAULT_SPACING_BOHR."""
    if not (text or "").strip():
        return ms.DEFAULT_SPACING_BOHR
    try:
        h = float(text)
    except ValueError as exc:
        raise ValueError(f"Grid spacing '{text}' is not a number.") from exc
    if not (MIN_SPACING <= h <= MAX_SPACING):
        raise ValueError(f"Grid spacing must be between {MIN_SPACING} and {MAX_SPACING} bohr.")
    return h


def resolve_output_dir(fchk_path: str, chosen: Optional[str]) -> str:
    """Chosen folder, or the folder of the .fchk when left blank."""
    chosen = (chosen or "").strip()
    return chosen if chosen else str(Path(fchk_path).resolve().parent)


# ----------------------------------------------------------------------
# cubegen engine
# ----------------------------------------------------------------------
def find_cubegen(explicit: Optional[str] = None) -> Optional[str]:
    """Path of Gaussian's cubegen: explicit, $GAUSS_EXEDIR, PATH, usual folders."""
    candidates = []
    if explicit and explicit.strip():
        candidates.append(explicit.strip())
    exedir = os.environ.get("GAUSS_EXEDIR")
    if exedir:
        candidates += [os.path.join(d, "cubegen.exe") for d in exedir.split(os.pathsep)]
    which = shutil.which("cubegen")
    if which:
        candidates.append(which)
    candidates += list(DEFAULT_CUBEGEN_LOCATIONS)
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None


def cubegen_potential(fchk_path: str, points_bohr: np.ndarray, cubegen: str,
                      timeout: float = 3600.0) -> np.ndarray:
    """ESP from Gaussian cubegen at arbitrary points.

    ``cubegen 0 potential=SCF file.fchk out.cub -5 h < points``: with -5 the
    points are read from stdin as free-format x y z **in Angstrom** and the
    result lists x y z V after the atom block. Raises RuntimeError on any
    failure, so the caller can fall back to the built-in engine.
    """
    points_angstrom = np.asarray(points_bohr, dtype=float) * BOHR_RADIUS_ANGSTROM
    env = dict(os.environ, GAUSS_EXEDIR=str(Path(cubegen).parent))
    with tempfile.TemporaryDirectory() as tmp:
        cube = Path(tmp) / "esp.cub"
        stdin = "".join(f"{x:.10f} {y:.10f} {z:.10f}\n" for x, y, z in points_angstrom)
        try:
            subprocess.run([cubegen, "0", "potential=SCF", str(fchk_path), str(cube), "-5", "h"],
                           input=stdin, text=True, env=env, timeout=timeout, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            lines = cube.read_text().splitlines()
        except (OSError, subprocess.SubprocessError) as exc:
            raise RuntimeError(f"cubegen failed: {exc}") from exc
    try:
        natoms = abs(int(lines[2].split()[0]))
        values = np.array([float(l.split()[3]) for l in lines[6 + natoms:] if l.strip()])
    except (IndexError, ValueError) as exc:
        raise RuntimeError("cubegen wrote an output file that could not be read.") from exc
    if values.size != len(points_angstrom):
        raise RuntimeError(f"cubegen returned {values.size} values for {len(points_angstrom)} points.")
    return values


# ----------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------
@dataclass
class SurfaceResult:
    isovalue: float               # electron density, a.u.
    mesh: ms.SurfaceMesh          # bohr
    esp: np.ndarray               # a.u., one per vertex
    extrema: List[ms.SurfaceExtremum]

    @property
    def half_range(self) -> float:
        return symmetric_range(self.esp)


@dataclass
class ESPResults:
    fchk_path: str
    wavefunction: gw.FchkWavefunction
    spacing: float
    engine: str
    surfaces: List[SurfaceResult]
    notes: List[str] = field(default_factory=list)
    elapsed: float = 0.0


class ESPSurfaceAnalyser:
    """Headless pipeline: fchk -> density grid -> isosurfaces -> ESP -> extrema."""

    def run(self, fchk_path: str, isovalues: Sequence[float] = (DEFAULT_ISOVALUE,),
            spacing: float = ms.DEFAULT_SPACING_BOHR, engine: str = ENGINE_BUILTIN,
            cubegen_path: Optional[str] = None,
            progress: Optional[ProgressCallback] = None) -> ESPResults:
        started = time.perf_counter()

        def report(fraction: float, message: str):
            if progress is not None:
                try:
                    progress(fraction, message)
                except Exception:
                    pass  # reporting must never abort an analysis

        if engine not in ENGINES:
            raise ValueError(f"Unknown ESP engine '{engine}'.")
        isovalues = list(isovalues) or [DEFAULT_ISOVALUE]
        notes: List[str] = []

        report(0.02, "Reading the formatted checkpoint...")
        wfn = gw.read_fchk(fchk_path)

        def density(points):
            return gw.electron_density(wfn, points)

        report(0.05, "Building the density grid...")
        radii = np.array([AtomicData.vdw_radii.get(s, AtomicData.DEFAULT_VDW_RADIUS)
                          for s in wfn.atom_symbols()]) / BOHR_RADIUS_ANGSTROM
        origin, shape = ms.grid_box(wfn.coordinates, radii, density, min(isovalues), spacing)
        values = density(ms.grid_points(origin, spacing, shape)).reshape(shape)

        cubegen = None
        if engine == ENGINE_CUBEGEN:
            cubegen = find_cubegen(cubegen_path)
            if cubegen is None:
                notes.append("cubegen was not found; the built-in exact engine was used instead.")
                engine = ENGINE_BUILTIN
        batches = ei.build_hermite_batches(wfn) if engine == ENGINE_BUILTIN else None

        surfaces: List[SurfaceResult] = []
        span = 0.65 / len(isovalues)
        for k, iso in enumerate(isovalues):
            base = 0.30 + k * span
            report(base, f"rho = {iso:g} a.u.: building the isosurface...")
            mesh = ms.build_isosurface(values, origin, spacing, iso, density, ms.density_outward(density))
            closure = mesh.component_euler_characteristics()
            if any(c != 2 for c in closure):
                notes.append(f"rho = {iso:g}: surface components with V+F-E = {closure} "
                             "(2 expected for each closed surface); try a smaller grid spacing.")
            report(base + 0.2 * span, f"rho = {iso:g} a.u.: ESP at {len(mesh.vertices)} surface points...")
            esp = None
            if engine == ENGINE_CUBEGEN:
                try:
                    esp = cubegen_potential(fchk_path, mesh.vertices, cubegen)
                except RuntimeError as exc:
                    notes.append(f"{exc} The built-in exact engine was used instead.")
                    engine = ENGINE_BUILTIN
                    batches = ei.build_hermite_batches(wfn)
            if esp is None:
                esp = ei.electrostatic_potential(
                    wfn, mesh.vertices, batches,
                    progress=lambda f, b=base: report(b + (0.2 + 0.75 * f) * span,
                                                      f"rho = {iso:g} a.u.: ESP {f:.0%}"))
            extrema = ms.find_surface_extrema(esp, mesh.vertices, mesh.neighbours)
            surfaces.append(SurfaceResult(iso, mesh, esp, extrema))

        report(1.0, "Done.")
        return ESPResults(str(fchk_path), wfn, spacing, engine, surfaces, notes,
                          time.perf_counter() - started)

    # --- reporting ------------------------------------------------------
    @staticmethod
    def nearest_atom_label(wfn: gw.FchkWavefunction, position_bohr: np.ndarray) -> str:
        i = int(np.argmin(np.linalg.norm(wfn.coordinates - position_bohr, axis=1)))
        return f"{wfn.atom_symbols()[i]}{i + 1}"

    def extrema_rows(self, results: ESPResults) -> List[Dict[str, object]]:
        """One dict per extremum of every surface, numbered from 1 in file order."""
        rows = []
        for surface in results.surfaces:
            for e in surface.extrema:
                rows.append({
                    "row": len(rows) + 1,
                    "isovalue": surface.isovalue,
                    "kind": e.kind,
                    "scope": e.scope,
                    "esp_au": e.value,
                    "esp_kcal_mol": e.value * HARTREE_TO_KCAL_MOL,
                    "position_angstrom": e.position * BOHR_RADIUS_ANGSTROM,
                    "nearest_atom": self.nearest_atom_label(results.wavefunction, e.position),
                })
        return rows

    def format_extrema_file(self, results: ESPResults) -> str:
        wfn = results.wavefunction
        out = [
            "# gQTEA ESP Surface Extrema",
            f"# source_file            {Path(results.fchk_path).name}",
            f"# density                {gw.DENSITY_SECTION}",
            "# surface_definition     electron density isosurface, rho in a.u. (e/bohr^3)",
            "# mapped_property        electrostatic potential V(r), a.u. (Hartree/e)",
            f"# esp_method             {results.engine}",
            f"# grid_spacing_bohr      {results.spacing:g}",
            f"# vertex_merge_bohr      {ms.MERGE_RATIO * results.spacing:g}",
            "# extremum_criterion     strict 1st+2nd-ring vertex neighbours (Multiwfn convention)",
            f"# hartree_to_kcal_mol    {HARTREE_TO_KCAL_MOL}",
            "# coordinates_unit       Angstrom",
            f"# atoms                  {wfn.n_atoms}  charge {wfn.charge}  multiplicity {wfn.multiplicity}",
        ]
        for s in results.surfaces:
            out.append(
                f"# surface rho_isovalue_au={s.isovalue:.6f} vertices={len(s.mesh.vertices)} "
                f"triangles={len(s.mesh.triangles)} "
                f"area_A2={s.mesh.area() * BOHR_RADIUS_ANGSTROM ** 2:.4f} "
                f"esp_min_au={s.esp.min():.8f} esp_max_au={s.esp.max():.8f}"
            )
        out.append("# row rho_isovalue_au type scope esp_au esp_kcal_mol x_A y_A z_A nearest_atom")
        for r in self.extrema_rows(results):
            x, y, z = r["position_angstrom"]
            out.append(
                f"{r['row']:d} {r['isovalue']:.6f} {r['kind']} {r['scope']} "
                f"{r['esp_au']:.8f} {r['esp_kcal_mol']:.4f} {x:.6f} {y:.6f} {z:.6f} {r['nearest_atom']}"
            )
        return "\n".join(out) + "\n"

    def write_extrema_file(self, results: ESPResults, output_dir: str) -> str:
        path = Path(output_dir) / (Path(results.fchk_path).stem + EXTREMA_FILE_SUFFIX)
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(self.format_extrema_file(results), encoding="utf-8")
        except OSError as exc:
            raise ValueError(f"Could not write {path}: {exc}") from exc
        return str(path)

    def summary_text(self, results: ESPResults, extrema_path: Optional[str] = None) -> str:
        wfn = results.wavefunction
        lines = [
            "ELECTROSTATIC SURFACE POTENTIAL",
            f"File        : {results.fchk_path}",
            f"Molecule    : {wfn.n_atoms} atoms, {wfn.n_electrons} electrons, "
            f"charge {wfn.charge}, multiplicity {wfn.multiplicity}",
            f"Density     : {gw.DENSITY_SECTION}",
            f"ESP engine  : {results.engine}",
            f"Grid        : {results.spacing:g} bohr",
            f"Time        : {results.elapsed:.1f} s",
            "",
        ]
        rows = self.extrema_rows(results)
        for s in results.surfaces:
            lines += [
                f"SURFACE  rho = {s.isovalue:.6f} a.u. (electron density isovalue)",
                f"  vertices {len(s.mesh.vertices)}, area {s.mesh.area() * BOHR_RADIUS_ANGSTROM ** 2:.2f} A^2",
                f"  ESP on this surface: min {s.esp.min():.6f} a.u. ({s.esp.min() * HARTREE_TO_KCAL_MOL:.2f} kcal/mol), "
                f"max {s.esp.max():.6f} a.u. ({s.esp.max() * HARTREE_TO_KCAL_MOL:.2f} kcal/mol)",
                "   row  type     scope      ESP (a.u.)  ESP (kcal/mol)   x (A)     y (A)     z (A)   atom",
            ]
            for r in rows:
                if r["isovalue"] != s.isovalue:
                    continue
                x, y, z = r["position_angstrom"]
                lines.append(f"  {r['row']:4d}  {r['kind']:7s}  {r['scope']:6s}  {r['esp_au']:12.6f}  "
                             f"{r['esp_kcal_mol']:12.3f}  {x:8.4f}  {y:8.4f}  {z:8.4f}   {r['nearest_atom']}")
            lines.append("")
        for note in results.notes:
            lines.append(f"NOTE: {note}")
        if extrema_path:
            lines.append(f"Extrema written to: {extrema_path}")
        return "\n".join(lines)


# ----------------------------------------------------------------------
# Toga UI
# ----------------------------------------------------------------------
OPACITY_LABELS = {"Opaque": 1.0, "Transparent (50 %)": 0.5}


class ElectrostaticSurfacePotentialUI:
    LABEL_WIDTH = 230
    HINT_COLOR = "#666666"

    def __init__(self, *args, **kwargs):
        # gqteaWinToga passes the Button that was pressed; ignore it.
        self.analyser = ESPSurfaceAnalyser()
        self.results: Optional[ESPResults] = None
        self.extrema_path: Optional[str] = None
        self.viewer: Optional[ESPSurfaceViewer] = None
        self._syncing = False
        self.window = toga.Window(title="Electrostatic Surface Potential (ESP)", size=(900, 660),
                                  on_close=self._on_close)
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

    def layout_main_window(self):
        self.fchk_input = toga.TextInput(placeholder="Gaussian formatted checkpoint (.fchk / .fch)",
                                         style=Pack(flex=1))
        self.output_dir_input = toga.TextInput(placeholder="Blank = folder of the .fchk file",
                                               style=Pack(flex=1))
        header = toga.Box(children=[
            self._row(self._label("Formatted checkpoint:"), self.fchk_input,
                      toga.Button("Browse", on_press=self.browse_fchk, style=Pack(width=90, margin_left=5))),
            self._row(self._label("Output folder:"), self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
        ], style=Pack(direction=COLUMN))

        self.isovalue_input = toga.TextInput(placeholder=f"{DEFAULT_ISOVALUE:.6f}", style=Pack(width=220))
        self.spacing_input = toga.TextInput(placeholder=f"{ms.DEFAULT_SPACING_BOHR}", style=Pack(width=120))
        self.engine_selection = toga.Selection(items=list(ENGINES), style=Pack(width=220))
        self.cubegen_input = toga.TextInput(value=find_cubegen() or "", placeholder=DEFAULT_CUBEGEN_LOCATIONS[0],
                                            style=Pack(flex=1))
        setup = toga.Box(children=[
            self._heading("Surface definition - electron density isosurface rho(r)"),
            self._row(self._label("Isovalue(s) rho  (a.u. = e/bohr^3):"), self.isovalue_input),
            self._hint("Blank = 0.000400 (GaussView default). Several: 0.0004, 0.001, 0.002"),
            self._row(self._label("Grid spacing (bohr):"), self.spacing_input),
            self._hint("Blank = 0.25 bohr (Multiwfn default for the ESP). Smaller = smoother and slower."),
            self._row(self._label("Density:"), toga.Label(gw.DENSITY_SECTION, style=Pack(margin_top=6))),
            self._heading("Mapped property - electrostatic potential V(r)"),
            self._row(self._label("ESP engine:"), self.engine_selection),
            self._row(self._label("cubegen executable:"), self.cubegen_input),
            self._hint("Only used with the Gaussian cubegen engine; on failure the built-in engine runs."),
        ], style=Pack(direction=COLUMN, margin=10))

        self.surface_selection = toga.Selection(items=[], on_change=self.on_display_change, style=Pack(width=260))
        self.unit_selection = toga.Selection(items=list(COLOUR_UNITS), on_change=self.on_display_change,
                                             style=Pack(width=140))
        self.opacity_selection = toga.Selection(items=list(OPACITY_LABELS), on_change=self.on_display_change,
                                                style=Pack(width=180))
        self.extrema_switch = toga.Switch("Show extrema markers", value=True, on_change=self.on_display_change)
        self.summary_output = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_esp_surface,
                                                      style=Pack(flex=1, font_family="monospace"))
        results = toga.Box(children=[
            self._row(self._label("Surface shown:"), self.surface_selection),
            self._row(self._label("Colour-bar unit:"), self.unit_selection),
            self._row(self._label("Surface:"), self.opacity_selection, self.extrema_switch),
            self.summary_output,
        ], style=Pack(direction=COLUMN, margin=10, flex=1))

        self.tabs = toga.OptionContainer(content=[
            ("Setup", toga.ScrollContainer(content=setup, style=Pack(flex=1))),
            ("Results", results),
        ], style=Pack(flex=1))

        self.run_button = toga.Button("Run", on_press=self.workflow, style=Pack(width=110, margin=5))
        self.show_button = toga.Button("Show surface", on_press=self.show_surface, style=Pack(width=130, margin=5))
        self.clear_button = toga.Button("Clear", on_press=self.clear_output, style=Pack(width=110, margin=5))
        self.close_button = toga.Button("Close", on_press=self.closeTopLevel, style=Pack(width=110, margin=5))
        self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
        self.progress_bar = toga.ProgressBar(max=100, style=Pack(margin=(0, 5, 5, 5)))

        self.window.content = toga.Box(children=[
            header,
            self._row(self.run_button, self.show_button, self.clear_button, self.close_button),
            self.tabs,
            self.status_label,
            self.progress_bar,
        ], style=Pack(direction=COLUMN, margin=10, flex=1))

    # --- dialogs ------------------------------------------------------------
    async def _show_error(self, title: str, message: str):
        await self.window.dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        await self.window.dialog(toga.InfoDialog(title, message))

    async def browse_fchk(self, widget):
        try:
            path = await self.window.dialog(toga.OpenFileDialog(
                title="Select a Gaussian formatted checkpoint", file_types=["fchk", "fch"]))
            if path:
                self.fchk_input.value = str(path)
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))

    async def browse_output_dir(self, widget):
        try:
            folder = await self.window.dialog(toga.SelectFolderDialog(title="Choose the output folder"))
            if folder:
                self.output_dir_input.value = str(folder)
        except Exception as exc:
            await self._show_error("Folder Selection Error", str(exc))

    # --- run ----------------------------------------------------------------
    def _collect_inputs(self) -> Dict[str, object]:
        fchk = (self.fchk_input.value or "").strip()
        if not fchk:
            raise ValueError("Please select a .fchk or .fch file.")
        if not os.path.isfile(fchk):
            raise ValueError(f"File not found: {fchk}")
        return {
            "fchk_path": fchk,
            "isovalues": parse_isovalues(self.isovalue_input.value),
            "spacing": parse_spacing(self.spacing_input.value),
            "engine": self.engine_selection.value or ENGINE_BUILTIN,
            "cubegen_path": (self.cubegen_input.value or "").strip() or None,
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
            output_dir = resolve_output_dir(params["fchk_path"], self.output_dir_input.value)
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
            self.extrema_path = self.analyser.write_extrema_file(results, output_dir)
            self.results = results
            self.summary_output.value = self.analyser.summary_text(results, self.extrema_path)
            self._populate_surface_selection()
            self.tabs.current_tab = 1
            self.status_label.text = f"Done in {results.elapsed:.1f} s. Extrema: {self.extrema_path}"
            self.progress_bar.value = 100
            if self.viewer is not None and self.viewer.viewer_window_is_open():
                self._push_to_viewer()
        except Exception as exc:
            self.status_label.text = "Calculation failed."
            self.progress_bar.value = 0
            await self._show_error("ESP Calculation Error", str(exc))
        finally:
            self._set_buttons(True)

    # --- display --------------------------------------------------------
    @staticmethod
    def surface_label(isovalue: float) -> str:
        return f"rho = {isovalue:.6f} a.u."

    def _populate_surface_selection(self):
        self._syncing = True
        try:
            labels = [self.surface_label(s.isovalue) for s in self.results.surfaces]
            self.surface_selection.items = labels
            self.surface_selection.value = labels[0]
        finally:
            self._syncing = False

    def selected_surface(self) -> Optional[SurfaceResult]:
        if self.results is None:
            return None
        for s in self.results.surfaces:
            if self.surface_label(s.isovalue) == self.surface_selection.value:
                return s
        return self.results.surfaces[0]

    def surface_markers(self, surface: SurfaceResult) -> List[SurfaceMarker]:
        rows = [r for r in self.analyser.extrema_rows(self.results) if r["isovalue"] == surface.isovalue]
        return [SurfaceMarker(r["position_angstrom"], r["kind"], r["scope"], str(r["row"])) for r in rows]

    def _push_to_viewer(self):
        surface = self.selected_surface()
        if surface is None or self.viewer is None:
            return
        wfn = self.results.wavefunction
        self.viewer.set_molecule(wfn.atom_symbols(), wfn.coordinates * BOHR_RADIUS_ANGSTROM)
        self.viewer.set_display_options(unit=self.unit_selection.value or UNIT_AU,
                                        opacity=OPACITY_LABELS.get(self.opacity_selection.value, 1.0),
                                        show_extrema=bool(self.extrema_switch.value))
        self.viewer.set_surface(surface.mesh.vertices * BOHR_RADIUS_ANGSTROM, surface.mesh.triangles,
                                surface.mesh.vertex_normals(), surface.esp, surface.isovalue,
                                self.surface_markers(surface))

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
            self.viewer.loop = asyncio.get_running_loop()
        self._push_to_viewer()
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)

    def clear_output(self, widget):
        self.results = None
        self.extrema_path = None
        self.summary_output.value = HelpGqteaWin.help_esp_surface
        self._syncing = True
        try:
            self.surface_selection.items = []
        finally:
            self._syncing = False
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

- [x] **Step 4: Run the tests**

Run: `..\Scripts\python.exe -m pytest tests/test_electrostaticSurfacePotential.py -q`
Expected: `25 passed` (with `C:\G09W\cubegen.exe` present; on a machine without Gaussian `24 passed, 1 skipped`).

- [x] **Step 5: Run the whole suite**

Run: `..\Scripts\python.exe -m pytest tests/ -q`
Expected: only the 4 pre-existing failures listed in Global Constraints; every new test passes.

- [x] **Step 6: Commit**

```bash
git add electrostaticSurfacePotential.py gqteaWinToga.py tests/test_electrostaticSurfacePotential.py
git commit -m "Implement the ESP surface tool: analyser, extrema file and UI

Replaces the stub. Runs fchk -> density grid -> isosurfaces -> exact or
cubegen ESP -> Multiwfn extrema in a worker thread with staged progress,
writes <name>_ESP_extrema.txt and opens the GaussView-style viewer.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Documentation

**Files:**
- Modify: `USER_MANUAL.md` (TOC line after `- [3D Molecular Viewer](#3d-molecular-viewer)`; new section before `#### Energy Plots`)
- Modify: `readme.md` (General Tools list, Repository Layout, Requirements, Supported files)
- Modify: `installation_requirements.txt` (optional numba block)
- Modify: `CLAUDE.md` (Tool-specific notes: new bullet after the molecular-design bullet, before `## Adding a new tool`)

- [x] **Step 1: USER_MANUAL.md — table of contents.** After the line `    - [3D Molecular Viewer](#3d-molecular-viewer)` insert:

```markdown
    - [Electrostatic Surface Potential (ESP)](#electrostatic-surface-potential-esp)
```

- [x] **Step 2: USER_MANUAL.md — section.** Insert immediately before the line `#### Energy Plots`:

````markdown
#### Electrostatic Surface Potential (ESP)

Use **Tools > Electrostatic Surface Potential (ESP)** to map the molecular
electrostatic potential onto an electron-density isosurface, view it
colour-coded in 3D (GaussView 6 style) and list the surface maxima and
minima, following Multiwfn's quantitative molecular surface analysis.

**Input.** A Gaussian formatted checkpoint (`.fchk` or `.fch`), made from the
`.chk` with `formchk name.chk name.fchk`. The `Total SCF Density` stored in the
file is used for both the surface and the potential; open-shell files work.

**Two different quantities.**

| Quantity | What it is | Unit |
|---|---|---|
| Electron-density isovalue | Defines *where* the surface is: the surface is the set of points where ρ(r) equals this value. | a.u. = e/bohr³ |
| ESP V(r) | The property *mapped onto* that surface: the colours and the extrema. | a.u. = Hartree/e, or kcal/mol |

**Steps.**

1. **Browse** to the `.fchk`/`.fch` file. Optionally choose an output folder
   (blank = the folder of the `.fchk`).
2. *Setup* tab, **Surface definition**: type one or several isovalues
   (`0.0004, 0.001, 0.002`); blank uses **0.000400 a.u.** Leave the grid
   spacing blank for 0.25 bohr (smaller is smoother and slower).
3. **Mapped property**: keep **Built-in exact**, or choose **Gaussian cubegen**
   and give its path (found automatically in `C:\G09W` or `%GAUSS_EXEDIR%`).
   If cubegen fails, the built-in engine runs and the summary says so.
4. Press **Run**. The *Results* tab lists, for every surface, its area, the
   ESP range and every extremum; the file `<name>_ESP_extrema.txt` is written.
5. Press **Show surface** for the 3D window. Choose the **Surface shown**
   when several isovalues were computed, the **Colour-bar unit** (a.u. or
   kcal/mol), opaque or transparent, and whether to show the extrema markers.
   Rotate with the left mouse button, zoom with the wheel.

**The 3D window.** The colour bar runs across the top: red is the most
negative ESP, green zero, blue the most positive, over the symmetric range
±max|V| of the shown surface, labelled like GaussView (`-2.097e-2 … 2.097e-2`).
The caption repeats the unit and the density isovalue of the surface.
Cyan spheres are minima, orange spheres maxima (larger = global); the number
beside each is its row in the extrema file.

**How extrema are found.** A surface point is a *local minimum* when its ESP is
strictly lower than at every neighbouring surface point **and** every
neighbour of those neighbours; a *local maximum* likewise, strictly higher.
Equal values give no extremum, and nearby extrema are not merged. The *global*
minimum/maximum is the lowest/highest point of that surface. Before this test
the surface points closer than half the grid spacing are merged and points
with only two or three neighbours are removed — Multiwfn's procedure.

**Output file** `<name>_ESP_extrema.txt`: a commented header (source file,
density, engine, grid spacing, criteria, conversion factor, one `# surface`
line per isovalue with its area and ESP range), then one row per extremum:

```text
# row rho_isovalue_au type scope esp_au esp_kcal_mol x_A y_A z_A nearest_atom
1 0.000400 minimum global -0.02096700 -13.1570 -2.245500 -0.062900 0.009700 C4
```

1 a.u. = 627.5094740631 kcal/mol (CODATA 2018); coordinates in Å.

**Limits.** Basis functions up to g (l = 4). With an effective core potential
the stored density lacks the core electrons and the nuclear charges are the
effective ones (as in Multiwfn) — fine on a molecular surface, which lies far
outside the core. Installing the optional `numba` package makes the built-in
ESP about five times faster.
````

- [x] **Step 3: readme.md.**
  - In `### General Tools`, after `- 3D molecular viewer` add: `- Electrostatic surface potential (ESP) mapped on electron-density isosurfaces from Gaussian `.fchk`/`.fch`, with Multiwfn-style surface extrema`
  - In `## Repository Layout`, after the `molecularViewer.py` bullet add: `- `electrostaticSurfacePotential.py`, `gaussianWavefunction.py`, `espIntegrals.py`, `molecularSurface.py`, `espViewer.py` (and the optional `espNumba.py`) for the ESP surface tool`
  - In `## Requirements`, after the `PyMuPDF` bullet add: `- optional: `numba` (about 5x faster ESP in the ESP surface tool)`
  - In `## Supported Input and Output Files`, after the `.xyz` bullet add: `- Gaussian formatted checkpoint files (`.fchk`, `.fch`) and the ESP extrema file `<name>_ESP_extrema.txt``

- [x] **Step 4: installation_requirements.txt.** After the `PyMuPDF` block add:

```text
# Optional. Compiles the electrostatic-potential kernel of the ESP surface
# tool (espNumba.py), about 5x faster on benzene. Without it the tool uses
# its numpy code path and gives the same numbers. Leave it out of packaged
# builds unless llvmlite is bundled too.
# numba
```

- [x] **Step 5: CLAUDE.md.** Delete the temporary `## In-progress work: ESP surface tool` section near the top (it only pointed at this plan), then insert this bullet as the last item of `## Tool-specific notes` (just before `## Adding a new tool`):

````markdown
- **[electrostaticSurfacePotential.py](electrostaticSurfacePotential.py) / [gaussianWavefunction.py](gaussianWavefunction.py) / [espIntegrals.py](espIntegrals.py) / [molecularSurface.py](molecularSurface.py) / [espViewer.py](espViewer.py) / [espNumba.py](espNumba.py)** — ESP mapped on electron-density isosurfaces from a Gaussian `.fchk`/`.fch` (Tools tab). One-way chain `electrostaticSurfacePotential → {espViewer → molecularViewer, molecularSurface, espIntegrals → gaussianWavefunction → help}`; **only `electrostaticSurfacePotential.py` and `espViewer.py` (via molecularViewer) import Toga/GL**, and a subprocess test keeps `toga` out of the other three. Spec: [docs/superpowers/specs/2026-10-14-esp-surface-design.md](docs/superpowers/specs/2026-10-14-esp-surface-design.md). **Two quantities, never conflated:** the surface is a `Total SCF Density` isovalue ρ (a.u. = e/bohr³, default `DEFAULT_ISOVALUE = 0.0004`); the mapped property is V (a.u. = Hartree/e, or kcal/mol via `HARTREE_TO_KCAL_MOL = 627.5094740631`). `gaussianWavefunction` reproduces Multiwfn's fch conventions — every Cartesian primitive normalised on its own, SP shells split with `P(S=P)` coefficients, pure shells mapped by `pure_to_cartesian` (Gaussian order 0,+1,−1,…), `FCHK_CARTESIAN_ORDER` with the fch f order (XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ) — and is pinned by `tr(P S) = N` and by **Gaussian 09 cubegen reference values** in `tests/data/esp_reference/` (ESP agrees to < 1e-6 a.u. on 5D/6D, pure f/g, Cartesian 10F, SP and open-shell fixtures in `tests/data/fchk/`; regenerate with `make_reference.py`, never at test time). **cubegen reads `-5` stdin points in Å, not bohr** (verified) — `cubegen_potential` converts. cubegen's own *density* is screened near nuclei (12.8 vs 290 a.u. at a water O nucleus), so the density test only compares the 1e-5 < ρ < 1e-2 band. `espIntegrals` is McMurchie–Davidson with the density contracted once into point-independent Hermite coefficients (`build_hermite_batches`), threads over point chunks, and uses `espNumba` automatically when `numba` imports (~5× faster; `BOYS_SERIES_LIMIT = 60`, because at 30 the asymptotic F₈ is still 1e-6 off). `molecularSurface` copies Multiwfn `surfana.f90`: grid spacing 0.25 bohr; box from atoms ± 1.7 r_vdW (`AtomicData.vdw_radii`) **grown until every face has ρ < iso/10** (needed at 0.0004); Kuhn 6-tetrahedra marching tetrahedra with vertices shared by global edge id; 3 exact-ρ bisections; merge of edge-connected vertices closer than 0.5 h to their midpoint, then removal of 2- and 3-neighbour vertices (sequential, order-dependent, does **not** guarantee a minimum spacing afterwards — don't test for one). **Extrema = Multiwfn's strict 1st+2nd-ring rule** (`second_ring`), no tolerance, no automatic merging, vertex positions; the global min/max is argmin/argmax and is reported as `global` even when a tie keeps it out of the strict test. Output `<stem>_ESP_extrema.txt` (header + `row rho_isovalue_au type scope esp_au esp_kcal_mol x_A y_A z_A nearest_atom`). `espViewer.ESPSurfaceViewer` subclasses `MolecularViewer` through two **no-op hooks** added to the base main loop (`render_scene_extras` inside the molecule matrix, `draw_screen_overlay` after the axis triad) plus a `window_title` attribute; GaussView look verified against a GaussView 6.1 screenshot: top colour bar, symmetric ±max|V| (`symmetric_range`), `-2.072e-2` labels (`gaussview_label`), hue sweep red→green(0)→blue (`hue_colours`), opaque + specular, black background. The colour maths are pure functions tested headless; the GL calls are not. Eteno (ethylene B3LYP/6-31G) at 0.0004 gives 2 π minima −0.02097 a.u. and 4 H maxima +0.01166 a.u.; GaussView's bar reads ±2.072e-2 because it interpolates from a coarse cube.
````

- [x] **Step 6: Check the help text and manual agree** — open both; the defaults (0.000400, 0.25 bohr), the colours (cyan = minimum, orange = maximum) and the file name (`<name>_ESP_extrema.txt`) must match what `electrostaticSurfacePotential.py` and `espViewer.py` do.

- [x] **Step 7: Commit**

```bash
git add USER_MANUAL.md readme.md installation_requirements.txt CLAUDE.md
git commit -m "Document the ESP surface tool

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8 (optional speed-up): numba kernel

numba is accepted as an **optional** dependency (spec D6). Measured on benzene B3LYP/6-31G(d), 20 000 points: numpy 4.7 s (12 threads) vs numba 0.9 s, agreement 4e-12 a.u.

**Files:**
- Create: `espNumba.py`
- Test: `tests/test_espIntegrals.py` (`TestNumba`, already written in Task 3)

**Interfaces:**
- Consumes: `espIntegrals.hermite_indices`, `HermiteBatch` objects.
- Produces: `espNumba.electronic_potential(batches, points)`, `espNumba._boys_into(n_max, t, out)`, `espNumba.BOYS_SERIES_LIMIT = 60.0`. `espIntegrals` picks it up automatically through its guarded import.

- [x] **Step 1: Install numba into the venv** (numba 0.68.0 accepts numpy 2.4.2 — verified with `pip install --dry-run`, no downgrade)

Run: `..\Scripts\python.exe -m pip install numba==0.68.0`
Expected: installs `numba-0.68.0` and `llvmlite-0.50.0`; numpy stays 2.4.2 (`..\Scripts\python.exe -c "import numpy; print(numpy.__version__)"` → `2.4.2`).

- [x] **Step 2: Run the numba tests to confirm they fail for the right reason**

Run: `..\Scripts\python.exe -m pytest tests/test_espIntegrals.py -q`
Expected: `TestNumba` still **skipped** (8) because `espNumba` does not exist yet — `ei.numba_available()` is False.

- [x] **Step 3: Implement** — `espNumba.py`:

```python
"""Optional numba kernel for the electronic ESP (imported by espIntegrals).

Same McMurchie-Davidson sum as espIntegrals.electronic_potential, but one
compiled loop over (point, primitive pair) instead of numpy arrays of shape
(pairs, points, n_tuv). Importing this module raises ImportError when numba
is not installed; espIntegrals then keeps the numpy path.
"""

from __future__ import annotations

from math import exp, pi, sqrt
from typing import Dict, Tuple

import numba
import numpy as np

# Boys: power series below this t, asymptotic form above. At 30 the
# asymptotic F_8 is still 1e-6 off; at 60 it is exact to double precision.
BOYS_SERIES_LIMIT = 60.0


@numba.njit(cache=True)
def _boys_into(n_max, t, out):
    if t < BOYS_SERIES_LIMIT:
        # F_n(t) = exp(-t) sum_k (2t)^k / ((2n+1)(2n+3)...(2n+2k+1))
        term = 1.0 / (2 * n_max + 1)
        total = term
        k = 1
        while term > 1e-17 * total:
            term *= 2.0 * t / (2 * n_max + 2 * k + 1)
            total += term
            k += 1
        out[n_max] = exp(-t) * total
    else:
        df = 1.0  # (2n-1)!!
        for m in range(1, n_max + 1):
            df *= 2 * m - 1
        out[n_max] = df / 2.0 ** (n_max + 1) * sqrt(pi / t ** (2 * n_max + 1))
    emt = exp(-t)
    for m in range(n_max, 0, -1):
        out[m - 1] = (2.0 * t * out[m] + emt) / (2 * m - 1)


@numba.njit(parallel=True, cache=True)
def _potential_kernel(l, p, centers, coeffs, points, kind, src1, src2, coef, out):
    """out[i] += sum_m sum_k coeffs[m, k] R^0_k(p_m, P_m - C_i).

    The recursion is a precomputed instruction list over the full tuv
    index of order l, applied level by level n = l-1 .. 0:
    kind 0 = base (-2p)^n F_n, 1/2/3 = x/y/z step with sources src1 (-1)
    and src2 (-2, or -1 for none) scaled by coef.
    """
    n_tuv = coeffs.shape[1]
    for i in numba.prange(points.shape[0]):
        prev = np.zeros(n_tuv)
        cur = np.zeros(n_tuv)
        f = np.zeros(l + 1)
        acc = 0.0
        for m in range(p.shape[0]):
            x = centers[m, 0] - points[i, 0]
            y = centers[m, 1] - points[i, 1]
            z = centers[m, 2] - points[i, 2]
            _boys_into(l, p[m] * (x * x + y * y + z * z), f)
            m2p = -2.0 * p[m]
            prev[0] = m2p ** l * f[l]
            for n in range(l - 1, -1, -1):
                count = (l - n + 1) * (l - n + 2) * (l - n + 3) // 6
                for j in range(count):
                    s = kind[j]
                    if s == 0:
                        val = m2p ** n * f[n]
                    else:
                        d = x if s == 1 else (y if s == 2 else z)
                        val = d * prev[src1[j]]
                        if src2[j] >= 0:
                            val += coef[j] * prev[src2[j]]
                    cur[j] = val
                prev, cur = cur, prev
            for k in range(n_tuv):
                acc += coeffs[m, k] * prev[k]
        out[i] += acc


_PLANS: Dict[int, Tuple[np.ndarray, ...]] = {}


def recursion_plan(l: int):
    """Instruction arrays for _potential_kernel, cached per l.

    Entries are listed in hermite_indices order; that order is a prefix-
    closed list (every order-(l-n) index set is a prefix of the order-l
    one), so position j in a level is the same tuv in every level.
    """
    if l in _PLANS:
        return _PLANS[l]
    from espIntegrals import hermite_indices
    full = _graded_indices(l)
    pos = {key: j for j, key in enumerate(full)}
    kind = np.zeros(len(full), np.int64)
    src1 = np.zeros(len(full), np.int64)
    src2 = -np.ones(len(full), np.int64)
    coef = np.zeros(len(full))
    for j, (t, u, v) in enumerate(full):
        if t > 0:
            kind[j], src1[j] = 1, pos[(t - 1, u, v)]
            if t > 1:
                src2[j], coef[j] = pos[(t - 2, u, v)], t - 1
        elif u > 0:
            kind[j], src1[j] = 2, pos[(t, u - 1, v)]
            if u > 1:
                src2[j], coef[j] = pos[(t, u - 2, v)], u - 1
        elif v > 0:
            kind[j], src1[j] = 3, pos[(t, u, v - 1)]
            if v > 1:
                src2[j], coef[j] = pos[(t, u, v - 2)], v - 1
    order = np.array([pos[key] for key in hermite_indices(l)])
    _PLANS[l] = (kind, src1, src2, coef, order)
    return _PLANS[l]


def _graded_indices(l: int):
    """All (t, u, v) with t+u+v <= l, sorted by total order: prefix-closed."""
    out = []
    for total in range(l + 1):
        for t in range(total, -1, -1):
            for u in range(total - t, -1, -1):
                out.append((t, u, total - t - u))
    return out


def electronic_potential(batches, points: np.ndarray) -> np.ndarray:
    points = np.ascontiguousarray(points, dtype=float).reshape(-1, 3)
    out = np.zeros(points.shape[0])
    for b in batches:
        kind, src1, src2, coef, order = recursion_plan(b.l)
        # coeffs are stored in hermite_indices order; the kernel uses graded order
        coeffs = np.zeros((b.coeffs.shape[0], len(kind)))
        coeffs[:, order] = b.coeffs
        _potential_kernel(b.l, b.p, np.ascontiguousarray(b.centers), coeffs, points,
                          kind, src1, src2, coef, out)
    return out
```

- [x] **Step 4: Run the tests (first run compiles, ~10 s)**

Run: `..\Scripts\python.exe -m pytest tests/test_espIntegrals.py tests/test_electrostaticSurfacePotential.py -q`
Expected: `33 passed` for test_espIntegrals.py (no skips) and `25 passed` for test_electrostaticSurfacePotential.py.

- [x] **Step 5: Commit**

```bash
git add espNumba.py
git commit -m "Add optional numba kernel for the ESP (about 5x faster)

Used automatically when numba imports; agrees with the numpy path to
1e-10 a.u. Boys series up to t=60 (the asymptotic F_8 is 1e-6 off at 30).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: End-to-end check in the real application

**Files:** none (verification only).

- [x] **Step 1: Full suite**

Run: `..\Scripts\python.exe -m pytest tests/ -q`
Expected: only the 4 pre-existing failures from Global Constraints.

- [x] **Step 2: Launch the app and run the user's file**

Run: `..\Scripts\python.exe gqteaWinToga.py` → **Tools > Electrostatic Surface Potential (ESP)** → Browse to `C:\AJCamargo\QCMM_2sem_2026\2026_10_08\Eteno.fch`, isovalue blank, **Run**.
Expected: *Results* tab shows `SURFACE rho = 0.000400 a.u.` with 2 minima ≈ −0.0210 a.u. (≈ −13.16 kcal/mol, nearest atom C) and 4 maxima ≈ +0.0117 a.u. (nearest H); `Eteno_ESP_extrema.txt` appears next to the `.fch`.

- [x] **Step 3: Compare with GaussView**

Press **Show surface**, rotate to look onto the π face, and compare with `C:\AJCamargo\QCMM_2sem_2026\2026_10_08\Eteno_ESP.png`.
Expected: same shape and colour sweep (red over the C=C π cloud, green body, cyan/blue toward the H atoms), colour bar across the top labelled `-2.097e-2 … 2.097e-2` (GaussView: ±2.072e-2 — it interpolates from a coarse cube), caption `ESP (a.u.) mapped on electron density isosurface rho = 0.000400 a.u.`, cyan marker on the π minimum, orange markers by the H atoms. Switch the unit to kcal/mol: labels become `-1.316e1 … 1.316e1`. Close the tool window: the GL window closes too.

- [x] **Step 4: Several isovalues and the cubegen engine**

Isovalues `0.0004, 0.001`, engine **Gaussian cubegen**, Run. Expected: two `SURFACE` blocks, `ESP engine : Gaussian cubegen`, values equal to the built-in run to ~1e-6; the *Surface shown* selector switches the 3D surface.

- [x] **Step 5: Report** — no commit; report the observed numbers and any visual difference from GaussView to the user.
