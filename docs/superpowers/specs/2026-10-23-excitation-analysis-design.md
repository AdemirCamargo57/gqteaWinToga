# Electron excitation analysis (Excitation tab) - design

Date: 2026-10-23. Branch: `feature/excitation-analysis`. Status: design approved by
the user on 2026-10-23 (scope, partition, integration, tab name, sections 1-5 of the
chat design); this file records it with the facts read from the Multiwfn source.

## 0. User decisions

| Question | Decision |
|---|---|
| Where | A new **Excitation** tab of the 3D Molecular Viewer, immediately after the NTO tab, drawing in the viewer's own window. |
| Inputs | Same as the NTO tab: Gaussian `.fchk`/`.fch` + the `.log`/`.out` of the same CIS/TDA/TD job (G09 and G16), or Multiwfn's plain-text "Excited State" format. No ORCA/CP2K/GAMESS. |
| Scope | Core (hole-electron indices, MO contributions, major MO transitions, transition dipole/f check, Δr, Λ, ghost-hunter index, classification, report) **plus** atom/fragment contributions and IFCT, 3D display + cubes, hole-electron Coulomb attraction, heat maps. |
| Partition | Mulliken-like and Hirshfeld, chosen by the user; Mulliken-like is the default; a note recommends Hirshfeld when the basis has diffuse functions or a contribution is negative. |
| Integration | Analytic (dipole and second-moment integrals) for every quantity linear in the densities; a Becke molecular grid for the non-linear ones. Multiwfn is used only for validation (dev-time script), never at run time. |
| Tab name | `Excitation`; its heading reads "Electron excitation analysis". |
| Approach | New Toga-free module `excitationAnalysis.py` + tab `excitationTab.py` following the NTO tab; the NTO tab is not refactored. |

## 1. Sources checked

* Multiwfn manual 2026.7.10, sections 3.21.A (inputs), 3.21.1 (hole-electron theory,
  Theories 1-6, classification table p. 245), 3.21.4 (Δr), 3.21.7 (ghost-hunter),
  3.21.8 (IFCT), 3.21.14 (Λ), 3.21.15 (major MO transitions); 4.18.1 (D-pi-A
  reference output), 4.18.4 (N-phenylpyrrole Δr and Λ), 4.18.8.1 (4-nitroaniline
  IFCT).
* Multiwfn source `Multiwfn_2026.10.18_src_Linux.zip` (sobereva.com; no official
  GitHub), `excittrans.f90`: `hole_electron` (1700-2414), `hole_ele_MOcontri`,
  `hole_ele_atmcontri_heatmap`, `atmcontri_holeele_Mulliken`,
  `atmcontri_holeele_Hirshfeld`, `delta_r` (3356), `lambda_excit` (3514), `IFCT`
  (4907), `majorMOtrans` (6118), `calc_holeele` (6656), `calc_transdens`,
  `loadexccoeff` (1138).

Facts taken from that code (where it differs from the manual text, the code is
what Multiwfn prints, and both are documented):

1. **Closed-shell factor.** Gaussian normalises closed-shell coefficients to
   Σw² = 0.5; Multiwfn multiplies hole, electron, their cross terms and the
   transition density by **2** (`wfntype==0`). MO contributions likewise x2.
2. **Hole/electron.** Local terms over all configurations (`->` added, `<-`
   subtracted); cross terms only between configurations of the same direction
   sharing the virtual (hole: i->l, j->l) or the occupied (electron: i->l, i->m)
   orbital, `<-` pairs subtracted. Equivalent to the matrices
   P_hole = C_occ (X Xᵀ − Y Yᵀ) C_occᵀ, P_ele = C_vir (Xᵀ X − Yᵀ Y) C_virᵀ, except
   that **Multiwfn drops configurations with |w| < 0.01 from the cross terms**
   (`skippair`, `calc_holeele(skipthres=0.01)`).
3. **Negative grid values are set to 0** before every integral (`where
   (holegrid<0) holegrid=0`), including the normalisation and the centroids.
4. **Transition density** T(r) = Σ (w + w') φ_i φ_a over all lines, x2 closed
   shell; transition dipole = −∫ r T(r) dr.
5. **Centroids, σ** are normalised by the integral of hole/electron; Δμ =
   −(R_ele − R_hole) x (∫hole + ∫ele)/2.
6. **H_CT = |(H_x D_x, H_y D_y, H_z D_z)| / D** (not H·u_CT as printed in the
   manual); t = D − H_CT; H = (|σ_ele| + |σ_hole|)/2; Δσ = |σ_ele| − |σ_hole|.
7. **HDI/EDI** = 100 sqrt(∫ρ²).
8. **Ghost-hunter** M_AC = Σ_{->} w² (ε_a − ε_i)/Σ_{->} w² − 1/D (D in bohr);
   de-excitations ignored.
9. **Coulomb attraction**: the manual defines E_C = ∫∫ ρ_hole ρ_ele/|r1−r2|; the
   code (option 18) integrates **C_ele(r1) C_hole(r2)/|r1−r2|** on the uniform grid
   (C = the normalised anisotropic Gaussians). We report both, labelled.
10. **Mulliken-like contributions**: terms with |w_i w_j| < 0.001 dropped, then the
    atomic values are **renormalised** to sum 1. Overlap = sqrt(hole x ele) only
    when both are positive, else 0. Hydrogens hidden by default.
11. **Hirshfeld contributions**: Multiwfn's built-in free-atom densities, its Becke
    grid, cross terms screened at 0.01, renormalised to sum 1.
12. **Δr** = Σ_{->} K² |⟨a|r|a⟩ − ⟨i|r|i⟩| / N with K = w + w' (the `<-` of the same
    pair) and N = Σ over **all lines** of their combined coefficient squared (so a
    `<-` line also adds w'²). Reproduced exactly.
13. **Λ** = Σ_{->} K² ∫|φ_i||φ_a| / N (same K, N) on a 15 x 170 Becke grid
    (radcut 18 bohr).
14. **IFCT**: fragment hole/electron = sums of atomic values; negatives set to 0;
    Q[R,S] = Θ_R,hole Θ_S,ele; Δp_R = Σ_S Q[S,R] − Σ_S Q[R,S]; intrinsic CT% =
    Σ_{R≠S} Q, LE% = Σ_R Q[R,R]; apparent CT% (two fragments) = |Δp_1|.
15. **Major transitions**: |w|² (x2 closed shell) > 5 %, sorted by |w|, labels
    `H-n -> L+m` (`Ha`/`Hb` unrestricted), `<-` shown with a minus sign.

## 2. Computation (`excitationAnalysis.py`, bohr, no Toga)

**Coefficients.** Per spin s: X_s[i, a] (`->`) and Y_s[i, a] (`<-`), shape
(n_occ_s, n_vir_s). Restricted: factor f = 2 on every density (so ∫hole = 1).

**Density matrices (AO basis, fchk Cartesian handling as `gaussianWavefunction`).**
P_hole = f Σ_s C_occ,s (X Xᵀ − Y Yᵀ) C_occ,sᵀ; P_ele = f Σ_s C_vir,s (Xᵀ X − Yᵀ Y)
C_vir,sᵀ; P_T = f Σ_s C_occ,s (X + Y) C_vir,sᵀ (unsymmetrised; densities use
(P + Pᵀ)/2). Parameter `cross_term_threshold` (default **0**, exact): when > 0,
configurations below it are removed from the cross terms only (Multiwfn = 0.01),
used by the validation tests.

**Analytic quantities** from new `espIntegrals.moment_matrices(wfn)`: S, ⟨x⟩,
⟨y⟩, ⟨z⟩, ⟨x²⟩, ⟨y²⟩, ⟨z²⟩ (McMurchie-Davidson, origin at 0, Cartesian
normalised components like `overlap_matrix`). Integrals ∫hole = tr(P_hole S) etc.;
centroids R = tr(P r)/∫; σ_λ² = tr(P λ²)/∫ − R_λ²; D, Δμ, Δσ, H, H_CT, t as in
section 1 item 5-6; transition dipole μ = −tr(P_T r), f_calc = 2/3 E |μ|².
Negative regions are **not** clamped for these (exact moments of the true
densities); the difference from Multiwfn's clamped grid values is reported only
for TD states (de-excitations), where ρ_hole can dip below 0.

**Becke-grid quantities** on `hirshfeldCharges.molecular_grid` (75 x 434):
Sm = ∫min(ρh, ρe), Sr = ∫sqrt(ρh ρe), HDI, EDI with ρ clamped at 0 (Multiwfn's
rule), plus grid ∫hole/∫ele as a quality check. Λ uses ∫|φ_i||φ_a| for the MOs of
the listed configurations, on the same grid.

**MO contributions** Θ_i = f Σ_a (X² − Y²), Θ_a = f Σ_i (X² − Y²); **major
transitions** as item 15; **ghost-hunter** item 8 with fchk orbital energies;
**Δr** item 12 with ⟨i|r|i⟩ = C_iᵀ r C_i.

**Atom contributions.** Mulliken-like: Θ_A = Σ_{μ∈A} (P S)_μμ for P_hole and P_ele
(symmetric P; equals Multiwfn's ½(Σ_{μ∈A,ν} + Σ_{μ,ν∈A}) form), then renormalised
to 1 (Multiwfn), with `pair_screen` (default 0; tests use 0.001). Hirshfeld:
Θ_A = ∫ w_A ρ on the Becke grid with `hirshfeldCharges` pro-atoms (H-Ar, no ECP;
other cases `ValueError` -> NOTE, Mulliken still available), renormalised.
Overlap = sqrt(h e) if both > 0 else 0; Δ = e − h.

**Fragments** typed `1-10; 11-13` (`;` separates fragments, each parsed by
`allBondAnalysis.parse_solute_index_ranges`); overlap between fragments is an
error; atoms left out form an automatic last fragment "rest" (NOTE). **IFCT** as
item 14, plus the atom-atom CT matrix Θ_A,hole Θ_B,ele (negatives -> 0).

**Coulomb attraction.** (a) E_C(C): Coulomb energy of the two normalised
anisotropic Gaussians, E = (2/√π) ∫₀^∞ Π_λ (1 + 2 s² Σ_λ)^(-1/2) exp(−s² d_λ²/(1 +
2 s² Σ_λ)) ds with Σ_λ = σ_ele,λ² + σ_hole,λ², d = R_ele − R_hole (1-D quadrature,
`scipy.integrate.quad`); always computed (cheap); compared with Multiwfn's
option 18. (b) E_C(ρ) exact = ∫ ρ_hole(r) V_ele(r) dr with V_ele from
`espIntegrals` given P_ele (new `density=` argument of `build_hermite_batches`),
integrated on the Becke grid; optional (switch), with progress.

**Open shell.** α and β handled per spin and summed into the densities (as
Multiwfn's combined index); MO contributions and major transitions list spin.
Restricted open-shell references refused (as NTO).

## 3. Classification (descriptor-based, explicit, may be "undetermined")

Thresholds (constants in the module, quoted in the report):

| Constant | Value | Justification |
|---|---|---|
| `CT_DISTANCE` | 2.0 Å | Guido et al. JCTC 9, 3118 (2013) Δr criterion; the manual shows D = Δr for one orbital pair, so the same length is applied to D. |
| `LOW_OVERLAP` | 0.3 | Peach et al. JCP 128, 044118 (2008) Λ value; the manual shows Sr = Λ for one orbital pair. |
| `t` sign | 0 | Manual: t < 0 = hole and electron not separated; clearly positive = separated. |
| `DIFFUSE_SPREAD` | 1.0 Å (Δσ) | No literature value; **our heuristic**, labelled as such, calibrated on the H₂CO/aug-cc-pVDZ Rydberg fixture. |
| `DIFFUSE_ELECTRON` | 50 % | Fraction of the electron on diffuse basis functions; **our heuristic** (rule 3). |
| `NORMALIZATION_TOLERANCE` | 0.01 | As the NTO tab; also applied to abs(∫hole − 1) and abs(∫ele − 1). |
| `IFCT_CT_MAJORITY` | 50 % | Intrinsic CT% above half of the excited electron. |

Rules, first match wins:
1. **Undetermined - incomplete data** if the normalisation or an integral fails the
   tolerance.
2. **Single-direction charge transfer (CT)**: D ≥ 2.0 Å and t > 0.
3. **Rydberg-like**: Δσ ≥ 1.0 Å, D < 2.0 Å and at least `DIFFUSE_ELECTRON` = 50 %
   of the electron (Mulliken population) on diffuse basis functions (shell exponent
   < `DIFFUSE_EXPONENT` = 0.05 bohr⁻² on a non-H atom or < 0.03 on H) - the manual's
   definition ("the electron mainly consists of very diffuse MOs"); 50 % is our
   heuristic. Without diffuse functions: NOTE "cannot be assigned in this basis".
   *Revised during implementation (2026-10-23)*: the first draft required Sr < 0.3,
   but H₂CO/aug-cc-pVDZ S4 (B₂, n→3p Rydberg) has Sr = 0.33 and would have been
   called centrosymmetric CT; the diffuse fraction separates the fixtures cleanly
   (Rydberg states 71-143 %, valence states 7-19 %). Sr stays in the evidence.
4. **Possible centrosymmetric / multidirectional CT**: D < 2.0, t < 0,
   Δσ ≥ 1.0, Sr ≥ 0.3 - needs fragment confirmation (IFCT CT% ≥ 50 %).
5. **Local excitation (LE)**: D < 2.0, t < 0, |Δσ| < 1.0, Sr ≥ 0.3.
6. **Intermediate / undetermined**: anything else, listing the indicators that
   disagree.

Independent checks reported next to the verdict: Δr ≥ 2.0 Å (CT-like), Λ < 0.3
(low overlap). Agreement → confidence "consistent"; disagreement → "low, the
orbital-based indices disagree" with the values.

Separate statements (never merged into the verdict):
* **Fragment-based assignment** (only with fragments): intrinsic CT% ≥ 50 % →
  CT, with the largest net transfers R → S; else LE located mainly on the fragment
  with the largest Q[R,R]. Labelled "depends on the fragment definition".
* **Orbital character** (π→π*, n→π*, d-d, MLCT/LMCT): **not assigned**. The report
  says it needs the hole/electron isosurfaces or NTOs, and lists the evidence it
  has (dominant MO pairs, atoms with the largest hole/electron contributions).
  MLCT/LMCT can be quantified by defining the metal as a fragment.
* **Ghost warning**: E < M_AC → "possible ghost CT state; a necessary, not a
  sufficient condition (Multiwfn manual 3.21.7)", only meaningful for one-dimensional CT.
* **Functional warning**: Δr ≥ 2.0 Å or Λ < 0.3 → pure/low-HF-exchange
  functionals are unreliable; use a range-separated hybrid (manual 3.21.4, Peach).

## 4. GUI: the Excitation tab (`excitationTab.ExcitationTabController`)

NTO-tab conventions: `_settle_structure`, `asyncio.to_thread` +
`call_soon_threadsafe` progress, `cubeTab._assign` writes, synchronous
`_syncing`-guarded handlers that `_start` async work, dialogs for every failure,
results kept when writing the report fails. Sections:

1. **Input** - Formatted checkpoint (Browse; fills an empty excitation field with
   `<stem>.log`/`.out`), Excitation file, Output folder, **Load**.
2. **Excited state and options** - state Selection, Partition (`Mulliken-like`,
   `Hirshfeld`), Fragments (TextInput, hint), Include hydrogens (switch, atom
   table/heat map), exact Coulomb as two mutually exclusive checkboxes *fast*
   and *precise* (added 2026-10-23 at the user's request; ticking one disables
   the other, neither = off; fast = 40 × 110 grid + skipping points with
   < 1e-6 of the hole, ~10× faster, ≤ 1.2e-4 relative), **Run**, **Clear**, status,
   progress bar.
3. **Results** - classification label (bold), summary line, descriptor
   `toga.Table` (Quantity, Value, Unit), fragment `toga.Table` (Fragment, Hole %,
   Electron %, Overlap %, Δ %), columns fitted with `qtaimTab.fit_table_columns`.
4. **Display** - Show (`Hole`, `Electron`, `Hole and electron`, `Sr overlap`,
   `Charge density difference`, `Transition density`, `C_hole and C_ele`),
   Isovalue (blank = 0.002, transition density 0.001 - Multiwfn's defaults),
   Grid spacing (0.2 bohr), colours (hole Blue, electron Green, Sr Orange,
   CDD + Green / − Blue, T + Cyan / − Magenta), Opacity, Centroids switch
   (spheres at R_hole/R_ele and an arrow hole → electron), **Display**.
5. **Export** - **Export cubes** (the shown quantities), **Send to Cube tab**,
   **Export 3D image** (capture first, then ask), Heat map Selection (`Atom
   composition`, `Fragment composition`, `IFCT matrix`, `Atom-atom CT matrix`) +
   **Plot heat map**.
6. Report box (help text until a run).

Rendering: `MolecularViewerUI.excitation_layer = CubeSceneLayer()`, drawn in
`render_scene_extras` and included in `get_scene_radius` (getattr-guarded).
Centroid spheres and the arrow are triangle meshes passed as ordinary
`RenderItem`s (no renderer change).

**Display grid**: one uniform grid (`molecularSurface.grid_box`, envelope
max(ρh, ρe, |T|) at the isovalue, 0.2 bohr default, `MAX_GRID_POINTS` 25 M),
all of ρh, ρe, T computed in one AO pass (`gaussianWavefunction.matrix_densities`);
Sr, CDD from them; C functions analytic. Cached per (results, spacing, isovalue).

## 5. Outputs

* Report `<stem>_state<n>_excitation.txt` (written after every run): source
  files, method/basis, reference, state line; normalisation block (as Multiwfn);
  major MO transitions; MO contributions (≥ 1 %); HOLE-ELECTRON INDICES (each
  with definition and unit); ORBITAL-BASED INDICES (Δr, Λ, M_AC and its two terms);
  TRANSITION DIPOLE (calculated vs the log, f calculated vs printed); ATOM and
  FRAGMENT contributions; IFCT; COULOMB; CLASSIFICATION (verdict, rule matched,
  evidence lines, independent checks, fragment-based assignment, orbital
  character statement, warnings); LIMITATIONS; REFERENCES; NOTES; FILES.
* Cubes `<stem>_state<n>_<quantity>.cube` (`hole`, `electron`, `Sr`, `CDD`,
  `transdens`, `Chole`, `Cele`), bohr, plain cubes.
* Heat maps `<stem>_state<n>_<kind>_heatmap.json`: a new plot-manifest shape
  `{"type": "heat_map", "matrix", "row_labels", "col_labels", "title", "xlabel",
  "ylabel", "colorbar", "vmin", "vmax"}` drawn by `displayPlots.draw_heat_map`
  (imported in `plotViewer.py` on the `draw_bar_comparison` line, above the
  `TkAgg` pin), validated by plotter's `_check_heat_map`, described as `R x C
  heat map`; opened with `launch_plot_viewer`.

**References** printed in the report (verified against the manual):
Multiwfn - T. Lu, F. Chen, J. Comput. Chem. 33, 580 (2012); T. Lu, J. Chem. Phys.
161, 082503 (2024). Hole-electron analysis - Z. Liu, T. Lu, Q. Chen, Carbon 165,
461 (2020). RMSD/H/t/C indices - T. Le Bahers, C. Adamo, I. Ciofini, J. Chem.
Theory Comput. 7, 2498 (2011). Δr - C. A. Guido, P. Cortona, B. Mennucci, C. Adamo,
J. Chem. Theory Comput. 9, 3118 (2013). Λ - M. J. G. Peach, P. Benfield, T.
Helgaker, D. J. Tozer, J. Chem. Phys. 128, 044118 (2008). Ghost-hunter - M.
Campetella et al., J. Comput. Chem. 38, 2151 (2017); corrected form J. Chem.
Phys. 154, 204102 (2021). IFCT - Multiwfn manual section 3.21.8 (no separate
publication; cite Multiwfn). Hirshfeld - F. L. Hirshfeld, Theor. Chim. Acta 44,
129 (1977). Becke grid - A. D. Becke, J. Chem. Phys. 88, 2547 (1988). Mulliken -
R. S. Mulliken, J. Chem. Phys. 23, 1833 (1955). Exciton binding discussion - J.
Chem. Phys. 143, 244905 (2015). Only the references for the analyses actually run
are printed.

## 6. Validation

Fixtures in `tests/data/excitation/` (gzip-compressed fchk where large;
attributed in its README):
* Multiwfn `examples\excit`: `D-pi-A` (CAM-B3LYP/6-31G(d), 5 states; manual
  4.18.1 prints S1-S3 indices), `N-phenylpyrrole` (CAM-B3LYP/6-31+G(d), diffuse;
  manual 4.18.4 prints Δr and Λ of S1-S5), `4-nitroaniline` (PBE0/6-311G*; manual
  4.18.8.1 prints IFCT of S2).
* The NTO fixtures (`tests/data/nto/`): CIS (no Y), TD(50-50) (triplets, Y), CH₃
  UB3LYP (open shell), benzene (degenerate), uracil.
* New G09W job `h2co_td_augccpvdz` (`TD(NStates=6) IOp(9/40=4)`, B3LYP/aug-cc-pVDZ):
  Rydberg states for the classification rule.

References: `tests/data/excitation/make_reference.py <Multiwfn.exe>` (dev-time,
never run by the suite) drives Multiwfn 2026.7.15 through stdin for each fixture
and state: 18 → 1 → state → 1 with a fine uniform grid (spacing 0.08 bohr,
extension as Multiwfn's default) → printed indices; option 18 Coulomb; 1 → 2
(MO contributions, threshold 0); 1 → 3 Mulliken and Hirshfeld (with hydrogens);
18 → 8 IFCT with the fixture's fragments; 18 → 4 Δr; 18 → 14 Λ; 18 → 15 major
transitions. Raw outputs stored as `<fixture>_state<n>_multiwfn.txt`.

Tolerances (Multiwfn prints limited digits; grid quantities carry Multiwfn's grid
error, measured by comparing its 0.08 and 0.12 bohr runs and recorded in the
README):
* orbital/analytic, same algorithm: Δr, MO contributions, major transitions,
  Mulliken contributions (with `pair_screen=0.001`), IFCT on those - to the
  printed digits;
* centroids, D, σ, H, t, Δμ, transition dipole (Multiwfn grid vs our analytic,
  with `cross_term_threshold=0.01`): within the measured grid error;
* Sm, Sr, HDI, EDI, Hirshfeld (both numerical): within the measured grid error;
* Λ (Multiwfn 15 x 170 vs our 75 x 434 grid): ≤ 2e-3;
* E_C(C) vs option 18: ≤ 1 % (Multiwfn skips C_ele < 1e-6 and the self point).
* Basis-free checks: tr(P_hole S) = tr(P_ele S) = 1 (CIS exactly; TD = Σ(X²−Y²)x f),
  tr(P_T S) = 0, |μ| and f against the Gaussian log's transition-dipole table and
  f to their printed precision, IFCT sums to 1, E_C(ρ) > 0, moment integrals
  against Becke-grid quadrature of the same AO products.
* Classification: each fixture state's verdict pinned (D-pi-A S1 LE, S2 CT;
  N-phenylpyrrole S1 LE, S4/S5 CT; a Rydberg state of H₂CO/aug-cc-pVDZ), plus
  synthetic descriptor sets for every rule, including "undetermined".
* UI headless (`__new__` + stub widgets): every `run()` and display parameter comes
  from `_collect_inputs` / `_collect_display_inputs`; no-flicker rules; Toga-free
  numerical modules (subprocess test).

## 7. Not implemented / unsupported (help, manual, errors)

ORCA, CP2K, BDF, GAMESS, Firefly; spin-flip TD; ROHF/ROKS references;
EOM-CC/CASSCF/multi-reference; relaxed densities (Gaussian `density=` /
Z-vector; the hole-electron densities are unrelaxed); transition magnetic dipole
and its density; transition dipole moment density cubes; TDM heat maps
(3.21.2), density-difference CT (3.21.3), transitions between excited states,
NTOs (own tab), atomic transition charges, excited-state natural orbitals,
CT spectrum, polarisation analysis, ECD/CPL g factors; Becke/Hirshfeld-I
partitions; Hirshfeld beyond Ar or with ECPs; h functions (fchk reader limit);
PBC.

## 8. Missing information / assumptions

* Gaussian prints coefficients to 5 decimals; IOp(9/40=4) (|w| > 1e-4) is required
  for reliable indices, and its absence is detected by the normalisation check.
* The classification thresholds other than D, Sr/Λ and t are our own and are
  labelled as heuristics; the classification is a summary of the descriptors, not
  a substitute for inspecting hole/electron isosurfaces (manual 4.18.4).
* Exact moments versus Multiwfn's clamped uniform-grid moments differ for TD
  states with de-excitations; the report states the method used.
