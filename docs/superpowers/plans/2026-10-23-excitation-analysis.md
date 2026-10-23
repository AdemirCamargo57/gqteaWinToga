# Electron Excitation Analysis (Excitation tab) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A new **Excitation** tab in the 3D Molecular Viewer, right after NTO, that characterises a CIS/TDA/TD excited state from a Gaussian `.fchk` + `.log` the way Multiwfn's main function 18 does (hole-electron indices, MO/atom/fragment contributions, IFCT, Δr, Λ, ghost-hunter, Coulomb attraction), classifies it with documented criteria, draws hole/electron/Sr/CDD/T/C functions, exports cubes, heat maps and a referenced report — validated against Multiwfn 2026.7.15 but never calling it.

**Architecture:** Toga-free numerics in `excitationAnalysis.py` (matrices, indices, contributions, IFCT, Coulomb, display volumes, heat-map entries, the `ExcitationAnalyser`), `excitationClassification.py` (thresholds, rules, references) and `excitationReport.py` (report text), built on `excitationData`, `gaussianWavefunction`, `espIntegrals` (new moment integrals, density argument), `hirshfeldCharges` (grid + pro-atoms). One Toga module, `excitationTab.py`, follows `ntoTab.py` and draws through a third `cubeRender.CubeSceneLayer` (`viewer.excitation_layer`). Chain: `molecularViewer → excitationTab → {cubeTab, cubeRender, qtaimTab, excitationAnalysis → {excitationClassification, excitationReport, excitationData, gaussianWavefunction, espIntegrals, hirshfeldCharges, molecularSurface, cubeFile, allBondAnalysis}}`.

**Tech Stack:** Python 3.13, numpy/scipy, optional numba (existing kernels only), Toga 0.5 WinForms, PyOpenGL, matplotlib (plot viewer), pytest.

**Spec:** [docs/superpowers/specs/2026-10-23-excitation-analysis-design.md](../specs/2026-10-23-excitation-analysis-design.md)

## Global Constraints

- Run from `venv/src/` with `..\Scripts\python.exe`; suite `python -m pytest tests/ -q`. Baseline failures (not ours): `test_meanResidenceTime.py::TestAdvancedOptionsExposed::test_every_run_parameter_is_supplied_by_the_ui` and three table-styling tests in `test_plotter.py`. The user has an **uncommitted edit of `molecularViewer.py`** (tab order: Box second, NTO last) — never revert or commit it unless asked; add the Excitation tab after NTO in whatever order is on disk, and stage `molecularViewer.py` hunks selectively (`git add -p`) or ask.
- Lengths **bohr** in numerical modules; Å only at report/viewer boundaries (`cubeFile.BOHR_RADIUS_ANGSTROM`); eV = Hartree × 27.211386245988 (`HARTREE_TO_EV`, CODATA 2018, as `uvVisSpectrum` if it defines one — reuse it, else define once in `excitationAnalysis`).
- `excitationAnalysis`, `excitationClassification`, `excitationReport` must not import Toga or `molecularViewer` (subprocess test). `excitationTab` imports `molecularViewer` only inside methods. No new numba importer (`NUMBA_KERNEL_MODULES` unchanged).
- Closed-shell factor **2** on hole, electron, transition density and MO contributions (Multiwfn `wfntype==0`); Σw² expected 0.5 restricted, 1 unrestricted.
- Defaults: display spacing 0.2 bohr, isovalue 0.002 (transition density 0.001), `MAX_GRID_POINTS = 25_000_000`, `NORMALIZATION_TOLERANCE = 0.01`, Becke grid 75 × 434, partition Mulliken-like, hydrogens hidden, exact Coulomb off.
- Classification constants exactly: `CT_DISTANCE = 2.0` Å, `LOW_OVERLAP = 0.3`, t sign 0, `DIFFUSE_SPREAD = 1.0` Å (heuristic, labelled), `IFCT_CT_MAJORITY = 50.0` %, `DIFFUSE_EXPONENT = 0.05` (non-H) / `0.03` (H) bohr⁻².
- File names: `<stem>_state<n>_excitation.txt`; `<stem>_state<n>_<field>.cube` with field in `hole, electron, Sr, CDD, transdens, Chole, Cele`; `<stem>_state<n>_<kind>_heatmap.json` with kind in `atoms, fragments, ifct, atom_ct`.
- Colours: hole Blue, electron Green, Sr Orange, CDD + Green / − Blue, T + Cyan / − Magenta; centroid spheres hole Blue / electron Green, arrow Orange.
- Widget writes via `cubeTab._assign`; handlers synchronous and `_syncing`-guarded, async work via `_start`; dialogs via `window.dialog(toga.ErrorDialog/InfoDialog/QuestionDialog)`; nothing raises out of a handler.
- Validation-only knobs (`cross_term_threshold`, `pair_screen`, `use_numba`) are **constructor attributes** of `ExcitationAnalyser`, not `run()` parameters, so the "every run() parameter comes from the UI" test stays exact.
- Never loosen a tolerance to make a Multiwfn comparison pass: find the algorithmic difference (the source is in the scratchpad `mwfn_src/x/Multiwfn_2026.10.18_src_Linux/excittrans.f90`) and either fix it or document it in the spec and README.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

---

### Task 1: Fixtures and Multiwfn references

**Files:**
- Create: `tests/data/excitation/README.md`, `tests/data/excitation/{D-pi-A,N-phenylpyrrole,4-nitroaniline}.{fchk.gz,out}` (from `C:\Users\ajc\Downloads\Multiwfn_2026.7.15_bin_Win64\examples\excit`), `tests/data/excitation/inputs/h2co_td_augccpvdz.gjf` + `.log` + `.fchk.gz`, `tests/data/excitation/make_reference.py`, `tests/data/excitation/<fixture>_state<n>_multiwfn.txt`, `tests/excitation_fixtures.py`

**Interfaces:**
- Produces: `excitation_fixtures.FIXTURES: Dict[str, Tuple[str, str]]` (name → (fchk file, excitation file)); `fixture_paths(name) -> Tuple[Path, Path]` (fchk decompressed once into a session temp dir; NTO fixture names resolve into `tests/data/nto/`); `multiwfn_text(name, state) -> str` (`pytest.skip` when absent); `parse_multiwfn(text) -> Dict[str, object]` (regex parser of the printed lines below); `FRAGMENTS: Dict[str, str]` (`"4-nitroaniline": "11-13; 1-10; 14-16"`, `"D-pi-A": "21-23; 11-20; 1-10; 24-26"`).

- [ ] **Step 1:** Copy the three Multiwfn example pairs; gzip the fchk files (`gzip -9`). Write `inputs/h2co_td_augccpvdz.gjf`: `#p B3LYP/aug-cc-pVDZ TD(NStates=6) IOp(9/40=4)`, C2v formaldehyde (geometry of `tests/data/nto/inputs/h2co_cis_ccpvdz.gjf`), `%chk=h2co_td_augccpvdz.chk`. Run with `GAUSS_EXEDIR=C:\G09W`: `g09.exe` then `formchk.exe`; require "Normal termination"; gzip the fchk; delete the chk.
- [ ] **Step 2:** Write `make_reference.py <Multiwfn.exe>` on the QTAIM/NTO pattern (temp dir, copy `settings.ini` with `nthreads` = CPU count, `subprocess.run([exe, fchk], input=script, capture_output=True, text=True, timeout=...)`, store stdout). Per fixture and state (D-pi-A 1, 2, 3; N-phenylpyrrole 1-5; 4-nitroaniline 2; h2co_td5050 1 + a triplet; h2co_cis 1; ch3_utd 1; h2co_td_augccpvdz all Rydberg candidates), one run each for:
  - `18 / 1 / <log> / <state> / 1 / 4 / <spacing 0.08 bohr> / 0 / 0 / 0 / q` (check the real `setgrid` prompts interactively once and encode them; record a second run at 0.12 bohr for the grid-error estimate);
  - Coulomb: same with low-quality grid (`1`) then post-processing `18`;
  - `18 / 1 / <log> / <state> / 2 / 0 / 3 / 1 / 7 / 0 / 0` (MO contributions threshold 0; Mulliken atoms with hydrogens) and the same with `3 / 2` (Hirshfeld);
  - `18 / 8 / 1 / <log> / <state> / <n> / <fragments...> / 0` and with partition 2, only for fixtures in `FRAGMENTS`;
  - `18 / 4 / <log> / 1-<nstates>`; `18 / 14 / <log> / 1-<nstates> / n`; `18 / 15` with the log as the boot file.
  Write `<fixture>_state<n>_multiwfn.txt` = the concatenated stdout sections, each preceded by `### <section>`.
- [ ] **Step 3:** Run it; check D-pi-A S1 prints `Sr index ... 0.51` and `D index:   0.522` (manual 4.18.1), N-phenylpyrrole S1 `Delta_r = ... 0.793` Å and `lambda = 0.68` (4.18.4), 4-nitroaniline S2 `CT(%): 75.6` (4.18.8.1). Record in the README the measured 0.08 vs 0.12 bohr differences for every grid quantity (these are the grid-error tolerances of Tasks 4-5).
- [ ] **Step 4:** Write `tests/excitation_fixtures.py` (`fixture_paths` with a module-level temp dir via `tempfile.mkdtemp`, `gzip.open` copy; `parse_multiwfn` returns floats keyed `sr`, `sm`, `d`, `t`, `hdi`, `edi`, `hole_centroid`, `electron_centroid`, `sigma_hole`, `sigma_electron`, `h`, `h_ct`, `transition_dipole`, `ghost`, `coulomb`, `mo_contributions` {mo: (hole, ele)}, `mulliken`/`hirshfeld` {atom: (hole, ele, overlap, diff)}, `ifct_*`, `delta_r` {state: Å}, `lambda` {state: value}, `major` {state: text}).
- [ ] **Step 5:** README: provenance (Multiwfn examples redistributed under its licence, manual sections and printed values; G09W A.02 job), regeneration commands, tolerance table.
- [ ] **Step 6: Commit** `test: excitation-analysis fixtures and Multiwfn 2026.7.15 references`.

### Task 2: `excitationData` — transition dipoles from the log, neutral RO message

**Files:**
- Modify: `excitationData.py` (new field and parser; message wording in `check_against_wavefunction`)
- Test: `tests/test_excitationData.py` (new cases; the existing `match="restricted open-shell"` still passes)

**Interfaces:**
- Produces: `ExcitationFile.transition_dipoles: Dict[int, Tuple[float, float, float]]` (state number → (X, Y, Z) a.u. from the **last** "Ground to excited state transition electric dipole moments (Au):" table; empty for plain text or when absent).

- [ ] **Step 1: Failing tests**

```python
def test_transition_dipoles_are_read_from_the_last_table():
    data = ed.read_excitation_file(nto_dir / "uracil.out")
    assert set(data.transition_dipoles) == {1, 2, 3}
    x, y, z = data.transition_dipoles[1]
    assert all(isinstance(v, float) for v in (x, y, z))

def test_dipole_table_matches_printed_oscillator():   # f = 2/3 E |mu|^2 to the printed digits
    data = ed.read_excitation_file(excitation_dir / "N-phenylpyrrole.out")
    s = data.state(1); mu = np.array(data.transition_dipoles[1])
    assert abs(2 / 3 * s.energy_ev / HARTREE_TO_EV * mu @ mu - s.oscillator) < 2e-4

def test_plain_text_has_no_dipoles(): parse_excitation_text("Excited State 1 1 5.0\n 5 -> 6 0.70711\n", "x").transition_dipoles == {}

def test_ro_message_is_not_nto_specific():
    with pytest.raises(ValueError) as err: ...   # existing RO fixture path
    assert "NTO" not in str(err.value)
```

- [ ] **Step 2:** Run `python -m pytest tests/test_excitationData.py -q` → new tests FAIL.
- [ ] **Step 3:** Implement `_last_dipole_table(lines)`: find last index of the header; skip the column-title line; read rows matching `^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)`; stop at the first non-matching line. Set it in `parse_excitation_text` for Gaussian files. Change the RO message to "… is a restricted open-shell (…) reference; excited-state analyses support restricted closed-shell and unrestricted references only."
- [ ] **Step 4:** Run `tests/test_excitationData.py tests/test_ntoAnalysis.py tests/test_ntoTab.py` → PASS.
- [ ] **Step 5: Commit** `feat(excitationData): transition dipoles from the Gaussian log`.

### Task 3: `espIntegrals` — moment integrals and an arbitrary density for the ESP

**Files:**
- Modify: `espIntegrals.py` (new `moment_matrices`; `build_hermite_batches(wfn, density=None)`)
- Test: `tests/test_espMoments.py`

**Interfaces:**
- Produces: `moment_matrices(wfn) -> Dict[str, np.ndarray]` with keys `"s", "x", "y", "z", "xx", "yy", "zz"`, each (n_cart, n_cart) over normalised Cartesian components (same layout as `overlap_matrix`), origin (0, 0, 0) bohr. `build_hermite_batches(wfn, density: Optional[np.ndarray] = None)` — `density` is a symmetric matrix in the **fchk basis**; default `wfn.density`; internally `wfn.transform.T @ density @ wfn.transform`.

- [ ] **Step 1: Failing tests**

```python
def test_overlap_key_equals_overlap_matrix(): assert np.allclose(m["s"], ei.overlap_matrix(w), atol=1e-14)

@pytest.mark.parametrize("name", [the six fchk fixtures: 5D, 6D, pure f, pure g, 10F, open shell])
def test_moments_match_becke_quadrature(name):
    w = gw.read_fchk(...); m = ei.moment_matrices(w)
    g = hirshfeldCharges.molecular_grid(w.coordinates, w.atom_symbols())
    ao = gw.evaluate_cartesian_aos(w, g.points)            # (P, n_cart)
    for key, f in {"x": g.points[:, 0], "yy": g.points[:, 1] ** 2, "s": 1.0}.items():
        ref = (ao * (g.weights * f)[:, None]).T @ ao
        assert np.allclose(m[key], ref, atol=2e-6 * np.abs(ref).max())

def test_density_dipole_matches_fchk():   # electronic dipole: -tr(P x) + nuclear = "Dipole Moment"
    p_cart = w.cart_density
    mu = -np.array([np.sum(p_cart * m[k]) for k in "xyz"]) + w.nuclear_charges @ w.coordinates
    assert np.allclose(mu, w.dipole, atol=1e-4)

def test_esp_of_an_explicit_density_equals_default():
    a = ei.electronic_potential(ei.build_hermite_batches(w), pts)
    b = ei.electronic_potential(ei.build_hermite_batches(w, density=w.density), pts)
    assert np.array_equal(a, b)
```

- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3: Implement** with McMurchie-Davidson, raising the right-hand angular momentum (exact, no extra recursion): for a Cartesian pair, `<a|x|b> = <a|b+1_x> + B_x <a|b>` and `<a|x²|b> = <a|b+2_x> + 2 B_x <a|b+1_x> + B_x² <a|b>`, where `<a|b+k_x>` uses `hermite_coefficients(sa.l, sb.l + 2, ...)[:, i, j+k, 0]` in x and the plain `[..., 0]` in y, z, times `pre = ckk (π/p)^{3/2}` and the two `component_angular_norm` factors (as `overlap_matrix`). One loop over shell pairs fills all seven matrices. `build_hermite_batches`: `cart = wfn.cart_density if density is None else wfn.transform.T @ density @ wfn.transform`.
- [ ] **Step 4:** Run the new file plus `tests/test_espIntegrals*.py tests/test_electrostaticSurfacePotential.py` → PASS.
- [ ] **Step 5: Commit** `feat(espIntegrals): dipole and second-moment integrals; ESP of any density matrix`.

### Task 4: `excitationAnalysis` core — matrices, analytic indices, orbital indices

**Files:**
- Create: `excitationAnalysis.py`
- Test: `tests/test_excitationAnalysis.py`

**Interfaces:**
- Consumes: `excitationData` (`ExcitedState`, `Configuration`, `SPIN_ALPHA/BETA`, `check_against_wavefunction`), `gaussianWavefunction.mo_coefficients`, `espIntegrals.moment_matrices`.
- Produces:

```python
HARTREE_TO_EV = 27.211386245988
BOHR_RADIUS_ANGSTROM = cubeFile.BOHR_RADIUS_ANGSTROM
@dataclass
class SpinCoefficients:
    spin: str; n_occ: int
    x: np.ndarray   # (n_occ, n_vir) excitation coefficients
    y: np.ndarray   # (n_occ, n_vir) de-excitation coefficients
@dataclass
class DensityMatrices:          # fchk basis (n_basis, n_basis)
    hole: np.ndarray; electron: np.ndarray
    transition: np.ndarray      # unsymmetrised P_T; densities use (P + P^T)/2
@dataclass
class MomentIndices:
    hole_integral: float; electron_integral: float; transition_integral: float
    hole_centroid: np.ndarray; electron_centroid: np.ndarray      # (3,) bohr
    sigma_hole: np.ndarray; sigma_electron: np.ndarray            # (3,) bohr
    transition_dipole: np.ndarray                                 # (3,) a.u., -tr(P_T r)
    # derived properties: d_components (|R_e - R_h|), d_index, dipole_change,
    # sigma_hole_norm, sigma_electron_norm, delta_sigma (3,), delta_sigma_index,
    # h (3,), h_ct, h_index, t_index, oscillator(energy_ev) -> float
@dataclass
class MOContribution: spin: str; mo: int; hole: float; electron: float   # fractions
@dataclass
class GhostHunter: value: float; term1: float; term2: float              # Hartree
@dataclass
class OrbitalIndex: value: float; pairs: List[Tuple[ed.Configuration, float]]  # combined K configs, contribution
def coefficient_matrices(state, wfn) -> Dict[str, SpinCoefficients]
def density_matrices(coeffs, wfn, factor, cross_term_threshold=0.0) -> DensityMatrices
def basis_moments(wfn) -> Dict[str, np.ndarray]       # fchk basis: T M T^T
def moment_indices(dm, moments) -> MomentIndices
def mo_contributions(coeffs, factor) -> List[MOContribution]           # nonzero only, spin then MO order
def major_transitions(state, wfn, threshold=0.05) -> str               # Multiwfn's 'H-4 -> L 81.9%, ...'
def ghost_hunter(state, wfn, d_bohr) -> GhostHunter
def delta_r(state, wfn, moments) -> OrbitalIndex                       # bohr
```

- [ ] **Step 1: Failing tests** (fixtures via `excitation_fixtures.fixture_paths`; references via `parse_multiwfn`; tolerances from the README of Task 1)

```python
def test_cis_densities_integrate_to_sum_of_squares():   # h2co_cis: no de-excitations
    r = core("h2co_cis", 1)                               # integral = f * sum w^2 exactly (MOs orthonormal)
    x2, _ = r.state.norm_terms()
    assert abs(r.moments.hole_integral - 2 * x2) < 1e-10 and abs(r.moments.electron_integral - 2 * x2) < 1e-10
def test_transition_density_integrates_to_zero(): abs(transition_integral) < 1e-10
def test_hole_and_electron_matrices_are_symmetric_and_psd_for_cis(): eigvalsh >= -1e-12
@pytest.mark.parametrize("name,state", [("D-pi-A", 1), ("D-pi-A", 2), ("N-phenylpyrrole", 4), ("ch3_utd", 1)])
def test_moments_match_multiwfn(name, state):       # cross_term_threshold=0.01 as Multiwfn
    ours, ref = ..., parse_multiwfn(multiwfn_text(name, state))
    assert np.allclose(ours.hole_centroid * A, ref["hole_centroid"], atol=TOL["centroid"])
    ... d, sigma, h, h_ct, t, transition_dipole likewise
def test_transition_dipole_matches_log(): |mu| vs data.transition_dipoles within 2e-3 relative (5-decimal coefficients)
def test_mo_contributions_match_multiwfn(): D-pi-A S1 MO52 hole 0.96207, MO57 electron 0.85411 to 5e-6; sums 1 within the normalisation deviation
def test_major_transitions_match_multiwfn(): text equality with the 18/15 line of each state
def test_ghost_hunter_matches_multiwfn(): D-pi-A S1 three numbers to the printed 3 decimals (eV), using Multiwfn's D
def test_delta_r_matches_manual(): N-phenylpyrrole S1..S5 = 0.793368, 1.846333, 2.455982, 3.105966, 3.752463 A to 1e-6
def test_unrestricted_has_two_spins_and_no_factor(): ch3_utd -> keys {"alpha","beta"}, factor 1, integrals ~ 1
def test_single_pair_limits():                        # manual: D == delta_r when one pair dominates
    synthetic state with one configuration -> moment D == delta_r value to 1e-10
```

- [ ] **Step 2:** Run → FAIL (module missing).
- [ ] **Step 3: Implement.**
  - `coefficient_matrices`: per spin (alpha only when restricted), `n_occ = wfn.n_alpha/n_beta`, `n_vir = wfn.n_mo - n_occ`; fill `x[i-1, a-n_occ-1]` or `y[...]`; out-of-range → `ValueError` naming the line (as `ntoAnalysis.transition_matrix`).
  - `density_matrices`:
    ```python
    for s in coeffs.values():
        c = gw.mo_coefficients(wfn, s.spin); co, cv = c[:, :s.n_occ], c[:, s.n_occ:]
        xh, xe = _screened(s.x, thr), _screened(s.y, thr)   # cross terms only, see below
        hole += factor * co @ (_cross(s.x, thr, "hole") - _cross(s.y, thr, "hole")) @ co.T
        ele  += factor * cv @ (_cross(s.x, thr, "ele")  - _cross(s.y, thr, "ele"))  @ cv.T
        trans += factor * co @ (s.x + s.y) @ cv.T
    ```
    `_cross(m, thr, kind)`: exact `m @ m.T` (hole) / `m.T @ m` (ele) when `thr == 0`; otherwise diagonal (local) terms from all of `m` plus off-diagonal terms from `m * (abs(m) >= thr)` — i.e. `local = diag(sum(m², axis))`, `cross = mk @ mk.T - diag(diag(mk @ mk.T))` with `mk` the screened matrix (Multiwfn `skippair`, item 2 of spec section 1).
  - `basis_moments`: `{k: wfn.transform @ v @ wfn.transform.T for k, v in ei.moment_matrices(wfn).items()}`, cached on the analyser per wavefunction.
  - `moment_indices`: `n_h = tr(P_h S)`, `R_h = [tr(P_h x)/n_h, ...]`, `σ_h,λ = sqrt(tr(P_h λλ)/n_h − R_h,λ²)`; same for electron; `transition_integral = tr(P_T S)`; `transition_dipole = −[tr(P_T x), ...]` (use `np.sum(P * M)` for non-symmetric P with symmetric M). Properties: `d_components = |R_e − R_h|`, `d_index = norm`, `dipole_change = −(R_e − R_h)(n_h + n_e)/2`, `delta_sigma = σ_e − σ_h`, `delta_sigma_index = |σ_e| − |σ_h|`, `h = (σ_e + σ_h)/2`, `h_ct = norm(h * (R_e − R_h)) / d_index` (0 when d_index = 0), `h_index = (|σ_e| + |σ_h|)/2`, `t_index = d_index − h_ct`, `oscillator(E) = 2/3 (E/HARTREE_TO_EV) |μ|²`.
  - `mo_contributions`: hole_i = f Σ_a (x² − y²), ele_a = f Σ_i (x² − y²); MO numbers 1-based (virtual = n_occ + index + 1).
  - `major_transitions`: sort configurations by |w| descending (stable), percent = w²·(2 if restricted), keep > 5 %, label `H` / `H-n`, `L` / `L+m`; unrestricted `Ha`/`Hb`, `La`/`Lb` with the same offsets; `<-` written `i <- a -x.x%`; join with `", "`; percent `f"{p:.1f}%"` (Multiwfn `f5.1` trimmed).
  - `ghost_hunter`: excitation lines only, energies from `alpha_orbital_energies` / `beta_orbital_energies` by spin; `term2 = 1/d_bohr`.
  - `delta_r`: per configuration MO centroid `⟨i|r|i⟩ = c_iᵀ r c_i` (fchk-basis moments); K = w + w' of the same (i, a, spin) `<-` line; N = Σ over `->` lines of K² + Σ over `<-` lines of w'² (spec item 12); value = Σ_{->} K²/N |⟨a|r|a⟩ − ⟨i|r|i⟩|.
- [ ] **Step 4:** Run `tests/test_excitationAnalysis.py` → PASS.
- [ ] **Step 5: Commit** `feat(excitationAnalysis): hole-electron density matrices, analytic indices, Δr, ghost-hunter`.

### Task 5: Becke-grid quantities — Sm, Sr, HDI, EDI, Λ

**Files:**
- Modify: `excitationAnalysis.py`; `hirshfeldCharges.py` (extract `hirshfeld_shares`)
- Test: `tests/test_excitationAnalysis.py`, `tests/test_hirshfeldCharges.py` (unchanged results)

**Interfaces:**
- Produces: `hirshfeldCharges.hirshfeld_shares(points, coordinates, symbols) -> np.ndarray` (P, n_atoms) stockholder weights (the block now inside `hirshfeld_populations`, which calls it); `excitationAnalysis.GridIndices(sm, sr, hdi, edi, hole_integral, electron_integral)`; `grid_indices(wfn, dm, grid, use_numba=None) -> GridIndices`; `lambda_index(state, wfn, grid, use_numba=None) -> OrbitalIndex`; `point_densities(wfn, dm, points, use_numba=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]` (ρh, ρe, T via `gw.matrix_densities` with the symmetrised P_T).

- [ ] **Step 1: Failing tests**: `hirshfeld_populations` unchanged on its fixtures (existing tests); `grid_indices` on D-pi-A S1/S2, N-phenylpyrrole S4 vs Multiwfn Sm, Sr, HDI, EDI within the README grid tolerance; grid integrals within 1e-4 of the analytic ones for CIS; Λ on N-phenylpyrrole S1..S5 vs `0.684853 0.563804 0.530928 0.198710 0.235255` to 2e-3; single-pair synthetic state: Λ == Sr to grid precision (1e-4).
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3: Implement.** Chunks of `hirshfeldCharges.POINT_CHUNK`; densities clamped `np.maximum(ρ, 0)` (Multiwfn item 3) before `min`, `sqrt(ρhρe)`, `ρ²`; HDI/EDI = `100*sqrt(Σ w ρ²)`. Λ: MOs of the `->` lines only, `orbital_values` for the involved columns per chunk, `ovlp[i, a] += (w * |φ_i|) @ |φ_a|`, K and N as Δr.
- [ ] **Step 4:** Run → PASS (plus `tests/test_hirshfeldCharges.py`).
- [ ] **Step 5: Commit** `feat(excitationAnalysis): Sm, Sr, HDI, EDI and Λ on the Becke grid`.

### Task 6: Atom and fragment contributions, IFCT

**Files:**
- Modify: `excitationAnalysis.py`
- Test: `tests/test_excitationContributions.py`

**Interfaces:**
- Produces:

```python
PARTITION_MULLIKEN, PARTITION_HIRSHFELD = "mulliken", "hirshfeld"
@dataclass
class AtomContributions:
    method: str; hole: np.ndarray; electron: np.ndarray        # (n_atoms,) fractions, sum 1
    overlap: np.ndarray; difference: np.ndarray; negative: bool  # any raw value < 0
@dataclass(frozen=True)
class Fragment: name: str; atoms: Tuple[int, ...]; automatic: bool   # 0-based atoms
@dataclass
class IFCTResult:
    fragments: List[Fragment]; hole: np.ndarray; electron: np.ndarray   # clamped >= 0
    matrix: np.ndarray                  # Q[R, S] = hole_R * electron_S
    population_change: np.ndarray       # sum_S Q[S, R] - sum_S Q[R, S]
    intrinsic_ct: float; intrinsic_le: float; apparent_ct: Optional[float]   # percent
    clamped: bool
def basis_atoms(wfn) -> np.ndarray                  # fchk basis function -> 0-based atom
def mulliken_contributions(dm, s_matrix, wfn, pair_screen=0.0, coeffs=None, factor=2.0) -> AtomContributions
def hirshfeld_contributions(wfn, dm, grid, use_numba=None) -> AtomContributions
def parse_fragments(text: str, n_atoms: int) -> Tuple[List[Fragment], List[str]]   # notes
def ifct(atoms: AtomContributions, fragments: List[Fragment]) -> IFCTResult
def atom_ct_matrix(atoms: AtomContributions) -> np.ndarray   # max(hole_A * electron_B, 0)
```

- [ ] **Step 1: Failing tests**
  - Mulliken with `pair_screen=0.001` vs Multiwfn's printed atom table (D-pi-A S1 with hydrogens; 4-nitroaniline S2) to the printed 2 decimals in %; exact (`pair_screen=0`) sums to 1 and differs from screened by < 1e-3.
  - Hirshfeld vs Multiwfn within the README tolerance; ECP or element > Ar → `ValueError` naming the element (use a synthetic wavefunction copy with `atomic_numbers` edited).
  - `parse_fragments("1-3; 5", 6)` → fragments `(0,1,2)`, `(4,)`, automatic rest `(3,5)` and one note; overlap `"1-3; 3-4"` → `ValueError("Atom 3 is in fragments 1 and 2.")`; out of range → `ValueError`; blank → `([], [])`.
  - `ifct` on 4-nitroaniline S2 Mulliken vs manual: hole 28.17/65.89/5.94 %, electron 3.02/29.52/67.46 %, CT% 75.688 (to the regenerated reference's digits); Q sums to 1; two fragments → apparent CT = |Δp_1|·100; negative fragment value clamped and `clamped` set.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3: Implement.** `basis_atoms`: `np.repeat([s.atom for s in shells], [s.n_functions for s in shells])`. Mulliken exact: `raw = np.bincount(atoms, weights=np.diag(P @ S))` (P symmetric, S = fchk-basis overlap from `basis_moments`). With `pair_screen > 0`: rebuild P_hole/P_ele from `coeffs` keeping only configuration pairs with `|w_i w_j| >= pair_screen` (local terms always pass when `w² >= pair_screen`, as Multiwfn compares the product), then the same bincount. Renormalise to sum 1 (Multiwfn); `negative = (raw < 0).any()`; overlap = `sqrt(h e)` where both > 0 else 0. Hirshfeld: `hirshfeld_shares` × clamped-free ρ (Multiwfn does not clamp here) × grid weights; renormalise. Fragments: split on `;`, each through `allBondAnalysis.parse_solute_index_ranges(part)` (1-based → 0-based), names `"Fragment 1"`…, rest named `"Rest"`. IFCT exactly spec item 14.
- [ ] **Step 4:** Run → PASS.
- [ ] **Step 5: Commit** `feat(excitationAnalysis): atom/fragment contributions (Mulliken-like, Hirshfeld) and IFCT`.

### Task 7: Hole-electron Coulomb attraction

**Files:**
- Modify: `excitationAnalysis.py`
- Test: `tests/test_excitationAnalysis.py`

**Interfaces:**
- Produces: `@dataclass class CoulombEnergies: gaussian: float; exact: Optional[float]` (Hartree, positive = attraction magnitude as Multiwfn prints); `coulomb_gaussian(moments: MomentIndices) -> float`; `coulomb_exact(wfn, dm, grid, progress=None, use_numba=None) -> float`.

- [ ] **Step 1: Failing tests**: two isotropic unit Gaussians at distance d → `erf(d / (2σ))/d` with σ² summed (analytic: `erf(d/sqrt(2Σ))/d` for Σ = σ_e² + σ_h²) to 1e-10; d = 0 → `sqrt(2/(πΣ))`... (limit of erf/d); D-pi-A S1 `coulomb_gaussian` vs Multiwfn option 18 (manual 0.287031 a.u. at medium grid; use the regenerated value) to 1 %; `coulomb_exact` for a CIS fixture positive and within 30 % of the Gaussian value (sanity, not equality); progress callback that raises is ignored.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3: Implement.**
  ```python
  def coulomb_gaussian(m):
      var = m.sigma_electron ** 2 + m.sigma_hole ** 2; d = m.electron_centroid - m.hole_centroid
      def f(s):
          den = 1 + 2 * s * s * var
          return np.prod(den ** -0.5) * np.exp(-np.sum(s * s * d * d / den))
      return 2 / np.sqrt(np.pi) * quad(f, 0, np.inf, limit=200)[0]
  ```
  `coulomb_exact`: `batches = ei.build_hermite_batches(wfn, density=sym(dm.electron))`, `V = ei.electronic_potential(batches, grid.points)` (or `espNumba` through `ei.electrostatic_potential`-style dispatch: write a small `_electronic(batches, pts, use_numba)` that picks `espNumba.electronic_potential` when available), `E = Σ w ρ_hole V` in chunks with progress.
- [ ] **Step 4:** Run → PASS.
- [ ] **Step 5: Commit** `feat(excitationAnalysis): hole-electron Coulomb attraction (C functions and exact)`.

### Task 8: `excitationClassification.py`

**Files:**
- Create: `excitationClassification.py`
- Test: `tests/test_excitationClassification.py`

**Interfaces:**
- Produces:

```python
CT_DISTANCE = 2.0; LOW_OVERLAP = 0.3; DIFFUSE_SPREAD = 1.0; NORMALIZATION_TOLERANCE = 0.01
IFCT_CT_MAJORITY = 50.0; DIFFUSE_EXPONENT = 0.05; DIFFUSE_EXPONENT_H = 0.03
VERDICT_INCOMPLETE = "Undetermined (incomplete data)"
VERDICT_CT = "Single-direction charge transfer (CT)"
VERDICT_RYDBERG = "Rydberg-like excitation"
VERDICT_MULTI_CT = "Possible centrosymmetric / multidirectional CT"
VERDICT_LE = "Local excitation (LE)"
VERDICT_INTERMEDIATE = "Intermediate / undetermined"
@dataclass(frozen=True)
class Descriptors:
    d: float; t: float; sr: float; delta_sigma: float            # Angstrom except sr
    delta_r: Optional[float]; lambda_: Optional[float]           # Angstrom, dimensionless
    normalization_deviation: float; hole_integral: float; electron_integral: float
    has_diffuse: bool; energy_ev: float; ghost_ev: Optional[float]
@dataclass
class Classification:
    verdict: str; rule: str; evidence: List[str]; checks: List[str]; confidence: str
    fragment_assignment: Optional[str]; orbital_note: str; warnings: List[str]
def has_diffuse_functions(wfn) -> bool
def classify(d: Descriptors, ifct=None) -> Classification     # ifct: excitationAnalysis.IFCTResult or None
REFERENCES: Dict[str, str]     # keys: multiwfn, hole_electron, le_bahers, delta_r, lambda, ghost, ghost_fix, ifct, hirshfeld, becke, mulliken, exciton
def references_for(keys: Iterable[str]) -> List[str]
```

- [ ] **Step 1: Failing tests**: one synthetic `Descriptors` per rule (first-match order of spec section 3), boundary values (`d == 2.0` with `t > 0` → CT; `t == 0` → not separated), Rydberg without diffuse functions → intermediate + note "cannot be assigned in this basis", disagreeing Δr/Λ → `confidence` starts with `"low"`, IFCT with CT% 75.7 → `fragment_assignment` names the largest net transfer `"Fragment 2 -> Fragment 3"`, ghost (E < M_AC) → warning text contains "necessary, not a sufficient", Δr ≥ 2 → functional warning; `orbital_note` always says the orbital character is not assigned; every `REFERENCES` value contains the volume/page from spec section 5; `has_diffuse_functions` true for N-phenylpyrrole and h2co_td_augccpvdz, false for D-pi-A.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3: Implement** the rules verbatim from spec section 3 (evidence lines quote the value and the threshold, e.g. `"D = 3.48 A >= 2.0 A (CT distance; Guido et al. 2013 Δr criterion applied to D)"`; Δσ evidence says "heuristic threshold of this program"). `has_diffuse_functions`: any shell with `min(exponents)` below the element's threshold. Keep the module free of numpy-heavy work and of Toga.
- [ ] **Step 4:** Run → PASS.
- [ ] **Step 5: Commit** `feat: descriptor-based classification of electronic excitations`.

### Task 9: `ExcitationAnalyser.run`, results, report

**Files:**
- Modify: `excitationAnalysis.py` (results dataclass, analyser, `write_report`)
- Create: `excitationReport.py`
- Test: `tests/test_excitationAnalysis.py` (run/report), `tests/test_excitationReport.py`

**Interfaces:**
- Produces:

```python
@dataclass
class ExcitationResults:
    fchk_path: str; excitation_path: str
    wavefunction; excitations; state                  # repr=False
    restricted: bool; factor: float
    coefficients: Dict[str, SpinCoefficients]; densities: DensityMatrices
    moments: MomentIndices; grid: GridIndices
    mo_contributions: List[MOContribution]; major: str
    ghost: GhostHunter; delta_r: OrbitalIndex; lambda_: OrbitalIndex
    partition: str; atoms: AtomContributions; partition_note: Optional[str]
    fragments: List[Fragment]; ifct: Optional[IFCTResult]
    coulomb: CoulombEnergies; classification: Classification
    sum_excitation: float; sum_deexcitation: float; expected_norm: float
    log_dipole: Optional[Tuple[float, float, float]]
    notes: List[str]; elapsed: float = 0.0
    # properties: deviation, stem, has_diffuse
class ExcitationAnalyser:
    def __init__(self, cross_term_threshold=0.0, pair_screen=0.0, use_numba=None)
    def load(self, fchk_path, excitation_path)       # (wfn, data, notes), cached by path/mtime/size
    def run(self, fchk_path, excitation_path, state_number, partition=PARTITION_MULLIKEN,
            fragments_text="", exact_coulomb=False, progress=None) -> ExcitationResults
    @staticmethod
    def write_report(results, output_dir) -> str      # <stem>_state<n>_excitation.txt; OSError -> ValueError
def report_filename(stem, state_number) -> str
def resolve_output_dir(fchk_path, chosen) -> str      # reuse ntoAnalysis.resolve_output_dir (import it)
# excitationReport.py
def report_text(results, written: Sequence[str] = ()) -> str
```

- [ ] **Step 1: Failing tests**: `run` on D-pi-A S1 with progress recording messages in increasing fraction order; a raising progress callback is ignored; Hirshfeld failure (ECP synthetic) → Mulliken kept? No — `ValueError` from the partition becomes `partition_note` and the run falls back to Mulliken with a NOTE; bad fragments text → `ValueError` before any heavy work; report contains every section heading of spec section 5 in order, each descriptor line with its unit, the classification verdict, "LIMITATIONS" and only the references of the analyses run (no `ifct` reference without fragments); unrestricted report lists alpha/beta MO contributions; written file round-trips.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3: Implement.** `run` order and progress: read (0.05) → coefficients/densities (0.15) → moments + Δr + ghost (0.25) → grid (0.35) → Sm/Sr/HDI/EDI (0.5) → Λ (0.65) → partition (0.75) → IFCT (0.8) → Coulomb (0.85-0.98) → classification (1.0). Grid built once per wavefunction and cached with the moments. Notes: normalisation deviation > 0.01 (IOp(9/40=4)), TD de-excitations present (moments exact vs Multiwfn's clamped grid), Mulliken with diffuse functions or negative contributions → "use Hirshfeld", automatic rest fragment, several blocks. `report_text` sections: header (files, method/basis, reference, electrons, basis size, state title, multiplicity, method used: "analytic moments + Becke grid 75 x 434"), NORMALIZATION (as Multiwfn), MAJOR MO TRANSITIONS, MO CONTRIBUTIONS (≥ 1 %), HOLE-ELECTRON INDICES (table: name, value, unit, definition), ORBITAL-BASED INDICES, TRANSITION DIPOLE (calc vs log, f calc vs printed), ATOM CONTRIBUTIONS (method; non-H unless all), FRAGMENTS / IFCT (when defined: contributions, Q matrix, net transfers, Δp, CT%/LE%, apparent for 2), COULOMB (both, a.u. and eV, labelled "C functions (Multiwfn's option 18)" / "ρ_hole ρ_ele (manual definition)"), CLASSIFICATION, LIMITATIONS (unrelaxed densities; Gaussian 5-decimal coefficients; Mulliken with diffuse functions; thresholds partly heuristic; isosurface inspection needed), REFERENCES, NOTES, FILES.
- [ ] **Step 4:** Run → PASS.
- [ ] **Step 5: Commit** `feat(excitationAnalysis): analyser, results and the referenced report`.

### Task 10: Display volumes, cubes and marker meshes

**Files:**
- Modify: `excitationAnalysis.py` (volumes, cube writing); `cubeRender.py` (`sphere_item`, `arrow_item`)
- Test: `tests/test_excitationAnalysis.py`, `tests/test_cubeRender.py`

**Interfaces:**
- Produces: `FIELDS = ("hole", "electron", "Sr", "CDD", "transdens", "Chole", "Cele")`; `SIGNED_FIELDS = {"CDD", "transdens"}`; `DEFAULT_ISOVALUES = {"transdens": 0.001}` (others 0.002); `field_volumes(results, spacing=0.2, isovalue=0.002, progress=None, use_numba=None) -> Dict[str, cubeFile.CubeVolume]`; `cube_filename(stem, state, field) -> str`; `ExcitationAnalyser.write_cubes(results, output_dir, fields, spacing, isovalue, progress=None) -> List[str]`; `cubeRender.sphere_item(centre_bohr, radius_bohr, rgb, opacity=1.0, subdivisions=2) -> RenderItem`, `cubeRender.arrow_item(start_bohr, end_bohr, radius_bohr, rgb, opacity=1.0) -> Optional[RenderItem]` (None for zero length).

- [ ] **Step 1: Failing tests**: on D-pi-A S1 at 0.25 bohr, `Σ hole · h³` within 1e-2 of 1, `Σ transdens · h³` ≈ 0, `CDD == electron − hole`, `Sr == sqrt(max(h,0) max(e,0))`, `Chole` integrates to the hole integral (Multiwfn normalises C to ∫hole), C centred at the centroid; `MAX_GRID_POINTS` exceeded → `ValueError` with the smallest workable spacing (reuse `ntoAnalysis.check_grid_size`); cube names; `sphere_item` vertices at the radius (1e-12), outward normals; `arrow_item` spans start→end; zero arrow → None.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3: Implement.** Box: `molecularSurface.grid_box(coords, vdW radii, envelope, isovalue, spacing)` with `envelope = max(|ρh|, |ρe|, |T|)` from `point_densities`; slabs of `ntoAnalysis.SLAB_POINTS`; C functions evaluated analytically on the same points and scaled so `Σ C h³ = ∫hole` (Multiwfn); titles name state, field and units. Spheres: icosahedron subdivided twice, normalised; arrow: 16-sided cylinder (radius r, 80 % of the length) + cone (radius 2r) as one triangle mesh with flat normals; both through `cubeRender._item("surface", ...)`.
- [ ] **Step 4:** Run → PASS.
- [ ] **Step 5: Commit** `feat: excitation display volumes, cubes and centroid/arrow meshes`.

### Task 11: Heat-map plot shape

**Files:**
- Modify: `displayPlots.py` (`draw_heat_map`, `render_figures_static` dispatch), `plotViewer.py` (import on the `draw_bar_comparison` line, dispatch `"heat_map"`), `plotter.py` (`_check_heat_map`, `describe_json_figures`), `excitationAnalysis.py` (`heat_map_entry`, `write_heat_map`)
- Test: `tests/test_displayPlots.py`, `tests/test_plotter.py`, `tests/test_excitationAnalysis.py`

**Interfaces:**
- Produces: manifest entry `{"type": "heat_map", "matrix": [[...]], "row_labels": [...], "col_labels": [...], "title": str, "xlabel": str, "ylabel": str, "colorbar": str, "vmin": float|None, "vmax": float|None}`; `displayPlots.draw_heat_map(axes, entry, font_style=None)`; `HEATMAP_KINDS = {"atoms": "Atom composition", "fragments": "Fragment composition", "ifct": "IFCT matrix", "atom_ct": "Atom-atom CT matrix"}`; `heat_map_entry(results, kind, include_hydrogens=False) -> dict` (`ValueError` for `fragments`/`ifct` without fragments); `heatmap_filename(stem, state, kind)`; `ExcitationAnalyser.write_heat_map(results, kind, output_dir, include_hydrogens=False) -> str` (a one-figure manifest, as plotter's table type writes).

- [ ] **Step 1: Failing tests**: `draw_heat_map` on an Agg axes produces one `QuadMesh`/`AxesImage`, tick labels equal the labels, a colour bar labelled `colorbar`; `validate_plot_manifest` accepts a good entry and rejects ragged rows / label-count mismatch / non-finite values with `"Figure 1: ..."` messages; `describe_json_figures` says `"3 x 16 heat map"`; atoms kind rows `["Hole", "Electron", "Overlap"]`, columns `"22(O)"`-style labels without hydrogens by default; ifct kind is Q with fragment names; the subprocess tests that keep `toga` out of `displayPlots`/`plotViewer` still pass.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3: Implement** `draw_heat_map` with `axes.imshow(matrix, aspect="auto", cmap="viridis", vmin, vmax)`, ticks at every label (rotate column labels 90° when more than 12), `axes.figure.colorbar(image, ax=axes, label=colorbar)`, title/labels with `font_style` as `draw_bar_comparison` does. Keep `from displayPlots import draw_bar_comparison, ..., draw_heat_map` **above** `matplotlib.use("TkAgg")` in `plotViewer.py`.
- [ ] **Step 4:** Run the three test files → PASS (baseline plotter failures unchanged).
- [ ] **Step 5: Commit** `feat(plots): heat-map figure shape for the Excitation tab`.

### Task 12: `excitationTab.py` and viewer integration

**Files:**
- Create: `excitationTab.py`
- Modify: `molecularViewer.py` (import, `excitation_layer`/`excitation_controller` in `__init__`, `_build_excitation_tab`, tab entry after NTO, `render_scene_extras`, `get_scene_radius`, tab-order test expectation), `help.py` (`HelpGqteaWin.help_excitation`)
- Test: `tests/test_excitationTab.py`, `tests/test_molecularViewer.py` (tab order)

**Interfaces:**
- Consumes: everything above; `cubeTab.CubeEntry/_assign/_item_names/canonical_colour`, `qtaimTab.fit_table_columns/parse_optional_float`, `ntoTab` patterns.
- Produces: `ExcitationTabController(viewer)` with `build_tab() -> toga.Box`, `_collect_inputs() -> Dict` (exactly the `run()` keyword parameters: `fchk_path, excitation_path, state_number, partition, fragments_text, exact_coulomb`), `_collect_display_inputs() -> Dict` (`fields: Tuple[str, ...]`, `spacing`, `isovalue`), `SHOW_MODES` (label → fields: `Hole`→(hole,), `Electron`→(electron,), `Hole and electron`→(hole, electron), `Sr overlap`→(Sr,), `Charge density difference`→(CDD,), `Transition density`→(transdens,), `C_hole and C_ele`→(Chole, Cele)), `PARTITION_LABELS`, `descriptor_rows() -> List[Tuple[str, str, str]]`, `fragment_rows()`, `compute_scene_items() -> Tuple[List[RenderItem], List[str]]`, handlers `load_files, workflow, display_excitation, export_cubes, send_to_cube_tab, export_3d_image, plot_heat_map, clear_output, on_show_change, on_opacity_change, on_colour_confirm, on_grid_confirm, on_centroids_change`.

- [ ] **Step 1: Failing tests** (headless: `ExcitationTabController.__new__` + stub widgets/viewer as `tests/test_ntoTab.py` does): every `run()` parameter comes from `_collect_inputs` (introspect `inspect.signature(ExcitationAnalyser.run)` minus `self`, `progress`); every `field_volumes` display parameter from `_collect_display_inputs`; Load fills the state Selection with `ExcitedState.title()`; Run stores results even when the report write fails (dialog "Excitation Output Error"); `descriptor_rows` contains D, Sr, t, Δσ, Δr, Λ, M_AC, E_C with units; classification label text = verdict; Show change re-colours cached meshes without recomputing volumes (count `field_volumes` calls); isovalue change re-meshes, spacing change recomputes; centroids switch adds/removes the two spheres and the arrow; heat-map without fragments for `fragments` kind → error dialog; tab order in the viewer is the on-disk order with `"Excitation"` immediately after `"NTO"`; `excitationTab` imports neither `molecularViewer` at module level (subprocess) nor heavy modules unnecessarily.
- [ ] **Step 2:** Run → FAIL.
- [ ] **Step 3: Implement** by following `ntoTab.py` section by section (spec section 4): same `_label/_row/_heading/_hint`, `LABEL_WIDTH = 210`, Browse/sibling fill, `_load` in `asyncio.to_thread`, `_progress_callback`, `_set_buttons`, `_settle_structure`, coalesced `refresh_scene` drawing into `viewer.excitation_layer`, `_start`, colour fields per field/sign (`COLOUR_KEYS = (("hole",1), ("electron",1), ("Sr",1), ("CDD",1), ("CDD",-1), ("transdens",1), ("transdens",-1), ("Chole",1), ("Cele",1))`, Chole/Cele reuse hole/electron colours), opacity slider, `Export 3D image` (capture first, then ask; default name `<stem>_state<n>_excitation.png`), `plot_heat_map` (`write_heat_map` → `displayPlots.launch_plot_viewer`, info dialog with the path when it cannot start), Send to Cube tab (as NTO, Show unticked). Results section: `classification_label` (bold), `summary_label` (`"State n: D = … Å, Sr = …, t = … Å; Δr = … Å, Λ = …"`), descriptor Table (Quantity, Value, Unit; fractions 0.5/0.3/0.2), fragment Table (Fragment, Hole %, Electron %, Overlap %, Δ %), both with `fit_table_columns` on Resize. Help text `help_excitation` (inputs, IOp(9/40=4), what each index means, classification rules and their sources, partitions, fragments syntax, display, exports, limitations).
- [ ] **Step 4:** Run `tests/test_excitationTab.py tests/test_molecularViewer.py tests/test_ntoTab.py` → PASS.
- [ ] **Step 5: Commit** `feat: Excitation tab in the 3D Molecular Viewer` (stage only this task's `molecularViewer.py` hunks).

### Task 13: Documentation, full suite, live check

**Files:**
- Modify: `USER_MANUAL.md` (Basic Workflows: *Analyse an Electronic Excitation*; tool section *Electron Excitation Analysis (Excitation tab)* with step-by-step procedure, inputs, every descriptor with definition/unit/interpretation, classification table with sources, partitions, IFCT, Coulomb (both definitions), display/exports, validation table, limitations; *Example 11: D-pi-A and N-phenylpyrrole* with the fixtures' real numbers; *Formulas and Methods (Excitation tab)*; references), `CLAUDE.md` (a tool-specific entry in the style of the NTO entry: chain, Multiwfn facts and the code-vs-manual differences, conventions, validation, not implemented; mention `heat_map` in the displayPlots/plotter bullets; `excitation_layer` in the molecularViewer integration notes), `tests/data/excitation/README.md` (final tolerances)
- Test: full suite

- [ ] **Step 1:** Write the docs; quote numbers only from the fixtures' actual runs (regenerate them with a short script, never by hand).
- [ ] **Step 2:** Run `python -m pytest tests/ -q`; compare against the baseline failures (plus any test the user's tab-order edit already broke on disk — report, don't fix silently).
- [ ] **Step 3:** Live check in the real viewer (the `run` skill / launch `python gqteaWinToga.py`): load D-pi-A, run S2, confirm CT verdict, display Hole and electron with centroids, export a cube and the IFCT heat map, export a 3D image; record what was checked in the commit message.
- [ ] **Step 4: Commit** `docs: Excitation tab in USER_MANUAL, help and CLAUDE.md`.
