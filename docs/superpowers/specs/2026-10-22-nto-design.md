# Natural Transition Orbitals (NTO tab) - design

Date: 2026-10-22. Branch: `feature/nto-analysis`. Status: approved by the user
on 2026-10-22 (placement, inputs, TDDFT treatment, validation source, design).

## 0. User decisions

| Question | Decision |
|---|---|
| Where | A new **NTO** tab of the 3D Molecular Viewer, after QTAIM and before Box & Performance, drawing in the viewer's own window (a second GLFW window does not work in this process). |
| Inputs | Gaussian `.fchk`/`.fch` (MOs and basis) **plus** a Gaussian `.log`/`.out` of a CIS / TDA / TD job, or Multiwfn's plain-text "Excited State" format. No ORCA/CP2K/GAMESS. |
| TD-DFT/TDHF | De-excitation (`<-`, Y) coefficients are **ignored** in T, as Multiwfn and Gaussian `Pop=NTO` do. No X+Y option. |
| Validation | Multiwfn 2026.7.15 (`C:\Users\ajc\Downloads\04_Software_and_Tools\Scientific_Software\Multiwfn\Multiwfn_2026.7.15_bin_Win64\Multiwfn_2026.7.15_bin_Win64\Multiwfn.exe`) through a dev-time `tests/data/nto/make_reference.py`, never run by the suite. |
| Design | Approach 1: own tab + own `CubeSceneLayer` instance, reusing the Cube tab's meshing and rendering; NTO `.fchk` export; cube names below; colours below. |

## 1. Sources checked

* Multiwfn manual 2026.7.10: section 3.21.A (inputs: `.fch` + Gaussian output with
  `IOp(9/40=4)`; closed-shell coefficients normalised to 0.5; plain-text format),
  3.21.6 (NTO theory), 4.18.6 (uracil example: S3 weights 0.865529, 0.134025,
  0.000582, ...; sum 1.000387).
* Multiwfn has **no official GitHub repository**; the official source is
  `http://sobereva.com/multiwfn/misc/Multiwfn_2026.10.18_src_Linux.zip`
  (GitHub repos are third-party mirrors/build recipes). Read: `excittrans.f90`
  `subroutine NTO` (3713-3894), `loadallexcinfo` (142-), `loadexccoeff` (1138-).
* Facts taken from that code (not assumed):
  * T(i, a-nocc) = w(i->a) for `->` lines only (`excdir==1`); `<-` ignored
    ("this treatment makes the result identical to that outputted by Gaussian").
  * Eigen-decomposition of T Tᵀ and Tᵀ T **separately** (DSYEV), eigenvalues x 2
    when there is no beta MO set (restricted), x 1 per spin otherwise.
  * Exported `.fch`: occupied NTO eigenvalues ascending (largest at the HOMO
    position), virtual descending (largest at the LUMO position); MO coefficients
    replaced, density untouched; orbital energies = eigenvalues.
  * Gaussian log: the **last** "Excitation energies and oscillator strengths:"
    block is used; a state's configuration lines are those containing `>` or `<`
    right after its header; multiplicity from "Singlet"/"Triplet" else 0;
    unrestricted indices carry `A`/`B`; an `A`->`B` (spin-flip) pair is refused.
  * Multiwfn prints "Sum of square of excitation coefficients", "Negative of the
    sum of square of de-excitation coefficients", their sum and the "Deviation to
    expected normalization value (0.5)".
* Running Multiwfn 2026.7.15 on its own `examples\excit\NTO\uracil.fch/.out`
  reproduces the manual (S3: 0.865529 0.134025 ...; sum 1.000387; deviation
  0.000007). Its exported `.fch` splits SP shells into s + p (56 vs 40 shell
  types; same functions) and writes values as `E16.8`.
* Real Gaussian lines can lack a space: `10 ->103  0.00017` - parse by regex.

## 2. Computation

**Inputs.** `FchkWavefunction` gains the MO coefficients (`Alpha MO coefficients`,
`Beta MO coefficients`, `n_basis x n_mo`, columns = MOs, fchk basis),
`Number of basis functions`, `Number of independent functions` (`n_mo`; can be
below `n_basis` with linear dependencies) and `Beta Orbital Energies`.

**Transition matrix** for spin s (alpha for restricted): `n_occ = n_alpha` (or
`n_beta`), `n_vir = n_mo - n_occ`, `T[i-1, a-n_occ-1] = w(i->a)` over that spin's
`->` lines. A repeated (i, a) line is an error naming the line.

**NTOs = SVD**, `T = U Σ Vᵀ` (`numpy.linalg.svd(full_matrices=True)`):
* weight `λ_k = f σ_k²`, f = 2 restricted, 1 unrestricted; descending;
  `k = 1..min(n_occ, n_vir)`.
* hole_k = `C_occ u_k`, particle_k = `C_vir v_k`, in the fchk basis.
* The SVD's eigenvalues equal Multiwfn's (σ² are the eigenvalues of both T Tᵀ
  and Tᵀ T) and it ties each particle to its hole (`v_k = Tᵀ u_k / σ_k`), which
  separate diagonalisations do not guarantee when weights are close.
* Sign: the one joint freedom (u_k, v_k) -> (-u_k, -v_k) is fixed by making the
  largest-magnitude element of u_k positive (first one on a tie). Null-space
  vectors (σ = 0) are kept only to complete the unitary sets for the fchk export.
* Composition: hole k over occupied MOs = `U[:, k]²`, particle k over virtual
  MOs = `V[:, k]²` (each sums to 1).
* Normalisation report, as Multiwfn: Σw² (`->`), −Σw'² (`<-`), their sum, the
  expected value (0.5 restricted, 1 unrestricted) and the deviation. A deviation
  above `NORMALIZATION_TOLERANCE = 0.01` adds a NOTE ("add IOp(9/40=4)").
  Σλ is reported; it can exceed 1 for full TD (Multiwfn's explanation is quoted).

**Interpretation** (help + manual): λ_k x 100 = % of the excitation carried by
pair k; one pair near 1 means the transition is that hole -> particle; if no pair
dominates, NTOs do not simplify the transition (hole-electron analysis would,
but is out of scope). Open shell: alpha and beta pairs are listed together,
sorted by weight; the largest says which spin carries the excitation.

**Grid.** Per displayed pair, one grid for hole and particle together:
`molecularSurface.grid_box` with `max(|hole|, |particle|)` as the function and
the isovalue (faces end up below isovalue/10), start box atoms ± 1.7 r_vdW,
spacing `DEFAULT_SPACING_BOHR = 0.2` (field; blank = default). Values are exact
at every point (`gaussianWavefunction.orbital_values`, no interpolation), in
slabs. `MAX_GRID_POINTS = 25_000_000` checked before evaluating (message gives
the smallest workable spacing, as the Fukui tool does).

**Performance.** The cost is the AO evaluation, already the compiled, screened
`densityNumba.cartesian_aos` when numba is installed; orbital values are then one
BLAS product per chunk, `(points x n_cart) @ (n_cart x k)`, with hole and particle
in one AO pass. The SVD is a small LAPACK call. **No new numba kernel** unless
profiling on benzene shows a Python loop above 20 % of the run (the QTAIM rule).
Meshes come from `volumeData.isosurface` (numba scan when installed).

## 3. Files read and checks

`excitationData.read_excitation_file(path)` detects the format: a Gaussian
output contains "Excitation energies and oscillator strengths:"; otherwise a
file with "Excited State" lines is the plain-text format; else a `ValueError`.

* Gaussian header: `Excited State   3:      Singlet-A"     6.0180 eV  206.02 nm  f=0.0000  <S**2>=0.000`.
* Plain text header: `Excited State   1   1   5.7945` (number, multiplicity, eV).
* Configuration line: `(\d+)([AB]?)\s*(->|<-)\s*(\d+)([AB]?)\s+(coefficient)`.
* Gaussian extras read for checks: last `NBasis=` / `NBsUse=`, `alpha electrons`
  / `beta electrons`, last "Standard orientation:" (else "Input orientation:")
  block.

`check_against_wavefunction(data, wfn)` raises `ValueError` (message names the
mismatch) when: the fchk has no MO coefficients; it is a restricted open-shell
reference (method starts with `RO`); restricted/unrestricted disagree between
the files; basis-function or electron counts differ; an orbital index is out of
range for its spin; the geometry differs by more than 1e-4 Å (message: "Use the
.fchk of the same job"). Several excitation blocks -> NOTE (last one used).

## 4. GUI: the NTO tab (`ntoTab.NTOTabController`)

Same conventions as the ESP/QTAIM tabs: `_settle_structure`, work in
`asyncio.to_thread` with `loop.call_soon_threadsafe` progress, widgets written
through `cubeTab._assign`, synchronous `_syncing`-guarded change handlers that
`_start` async work, dialogs for every failure, nothing raised out of a handler.

Sections:
1. **Input** - Formatted checkpoint (Browse; fills an empty excitation field from
   `<stem>.log` then `<stem>.out`), Excitation file (Browse), Output folder
   (blank = folder of the .fchk), **Load**.
2. **Excited state** - state Selection (`3: Singlet-A" 6.018 eV 206.0 nm f=0.0000`),
   **Run** (computes the NTOs and writes the report automatically), summary line
   (normalisation), progress bar, status.
3. **NTO pairs** - `toga.Table` (Pair, Spin, Weight, %, Cumulative %), the pairs
   with weight ≥ `LISTED_WEIGHT = 1e-4` (at least the top one), columns fitted
   with `qtaimTab.fit_table_columns`; selecting a row selects the pair.
4. **Display** - Pair Selection, Show (`Hole`, `Particle`, `Hole and particle`),
   Isovalue (default 0.02, the Cube tab's orbital default), Grid spacing,
   four colour fields (`Hole +` Blue, `Hole −` Red, `Particle +` Green,
   `Particle −` Yellow; `volumeData.parse_colour` names or `#rrggbb`), Opacity,
   **Display**. Changing show/colours/opacity re-colours cached meshes; a new
   isovalue re-meshes; a new pair or spacing recomputes the grid.
5. **Export** - **Export cubes** (selected pair; switch "all listed pairs"),
   **Export NTO .fchk**, **Send to Cube tab** (adds the pair's two volumes as
   unsaved Cube-tab entries), **Export 3D image** (capture first, then ask).
6. Report box (help text until a run).

Rendering: `MolecularViewerUI.nto_layer = CubeSceneLayer()` drawn in
`render_scene_extras` and included in `get_scene_radius` (getattr-guarded).

## 5. Outputs

* Report `<stem>_state<n>_NTO.txt`: files, method/basis, state line, normalisation
  block, weights of the listed pairs with their hole/particle MO numbers in the
  exported fchk and their dominant MO contributions (≥ 5 %), Σλ, notes.
* Cubes `<stem>_state<n>_NTO<k>_hole.cube` / `_particle.cube`; unrestricted:
  `<stem>_state<n>_NTO<k>_alpha_hole.cube` etc. Plain cubes (bohr), values =
  orbital amplitude (a.u.), title naming state, pair, kind, spin and weight.
* `<stem>_state<n>_NTO.fchk`: the input fchk text with only the MO-coefficient
  and orbital-energy sections replaced (Multiwfn's ordering); the SCF density is
  unchanged by occupied-occupied / virtual-virtual rotations, so the file stays
  valid for GaussView, Multiwfn, cubegen and this program.

## 6. Validation

Fixtures in `tests/data/nto/`: Multiwfn's uracil example (PBE0/6-31G*, TD, 3
states; attributed in `tests/data/nto/README.md`) plus G09W jobs with
`IOp(9/40=4)` (inputs kept): formaldehyde TD-B3LYP/6-31G(d) with
`TD(50-50)` (singlets + triplets, de-excitations), formaldehyde CIS/cc-pVDZ
(pure d shells, weights ≤ 1, Σλ ≈ 1), CH₃ UB3LYP/6-31G(d) TD (open shell, A/B
labels), benzene TD-B3LYP/6-31G(d) (degenerate e shells).

Basis-free checks: P_SCF = Σ n C Cᵀ; Cᵀ S C = I; Σ_occ |φ|² from `orbital_values`
equals `electron_density` at random points; U Σ Vᵀ = T; NTOs orthonormal;
uracil S3 weights equal the manual's printed digits.

Multiwfn references (`make_reference.py`): for each fixture and state, Multiwfn
18 -> 6 -> file -> state -> 2 (export `.fch`). The exported file is stored
gzip-compressed as `<fixture>_state<n>_multiwfn.fch.gz` (one per checked state)
and read back with our own fchk reader, so its basis (SP split into s + p) needs
no special handling. Its orbital energies are the NTO eigenvalues; hole k is MO
`n_occ + 1 - k` and particle k is MO `n_occ + k`. Comparison is basis-agnostic,
on a shared set of points around the molecule:
* weights: |Δλ| ≤ 2e-8 (Multiwfn writes 8 digits);
* non-degenerate pair (relative gap to its neighbours > 1e-6): max over points of
  |φ_ours − s φ_ref| with the best sign s = ±1 ≤ 1e-6 x max|φ|;
* degenerate group: least-squares map R from ours to Multiwfn's group; R must be
  orthogonal to 1e-6 and the residual ≤ 1e-6 x max|φ| (rotations within the
  subspace are arbitrary);
* only pairs with λ ≥ 1e-4 are compared (the rest are null-space vectors).

## 7. Not implemented / unsupported (help, manual, errors)

ORCA, CP2K, BDF, GAMESS and Firefly files; spin-flip TD; restricted open-shell
references; EOM-CC, CASSCF and other multi-reference methods; an X+Y transition
density; NTOs between two excited states; hole-electron analysis; NTO energies;
h functions (fchk reader limit); PBC.

## 8. Missing information / assumptions

* Gaussian prints coefficients to 5 decimals; agreement with Gaussian's internal
  NTOs is therefore ~1e-4, while agreement with Multiwfn (same printed input)
  is to round-off. Both facts go in the manual.
* The log and fchk must come from the same job (checked by geometry, basis and
  electron counts; the method is not compared because the fchk names e.g.
  `RPBE1PBE` and the log route `PBE1PBE`).
