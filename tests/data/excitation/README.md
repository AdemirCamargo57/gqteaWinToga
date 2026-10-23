# Excitation-analysis test fixtures

Used by `tests/test_excitationAnalysis.py`, `tests/test_excitationContributions.py`,
`tests/test_excitationClassification.py`, `tests/test_excitationReport.py` and
`tests/test_excitationTab.py` through `tests/excitation_fixtures.py` (which also
reaches the NTO fixtures in `tests/data/nto/`: `h2co_td5050`, `h2co_cis`, `ch3`, `uracil`).

| Files | Origin | Covers |
|---|---|---|
| `D-pi-A.fchk.gz`, `D-pi-A.out`, `inputs/D-pi-A.gjf` | Multiwfn 2026.7.15 `examples\excit` (Gaussian, CAM-B3LYP/6-31G(d) `TD(nstates=5) IOp(9/40=4)`). Manual 4.18.1 prints S1: Sr 0.51896, D 0.522 A, t -0.998 A, HDI 22.68, EDI 17.11 | local (S1, n->pi*) and charge-transfer (S2) states; IFCT with four fragments |
| `N-phenylpyrrole.fch.gz`, `N-phenylpyrrole.out` | Multiwfn `examples\excit` (CAM-B3LYP/6-31+G(d), 5 states). Manual 4.18.4 prints Delta_r 0.793368 ... 3.752463 A and Lambda 0.684853 ... 0.235255 | diffuse functions (Mulliken unreliable), Delta_r / Lambda |
| `4-nitroaniline.fch.gz`, `4-nitroaniline.out`, `inputs/4-nitroaniline.gjf` | Multiwfn `examples\excit` (PBE0/6-311G*). Manual 4.18.8.1 prints the S2 IFCT (CT 75.688 %) | IFCT with three fragments |
| `h2co_td_augccpvdz.fchk.gz`, `h2co_td_augccpvdz.log`, `inputs/h2co_td_augccpvdz.gjf` | G09W A.02, B3LYP/aug-cc-pVDZ `TD(NStates=6) IOp(9/40=4)`, formaldehyde at the geometry of the NTO CIS fixture (2026-10-23) | Rydberg-like states (very diffuse electron) |
| `<fixture>_state<n>_multiwfn.txt.gz` | Multiwfn 2026.7.15 output written by `make_reference.py` | hole-electron indices (two grids), Coulomb attraction, MO / Mulliken / Hirshfeld contributions, IFCT |
| `<fixture>_states_multiwfn.txt.gz` | idem | Delta_r, Lambda (Multiwfn's automatic grid and 75 x 434), major MO transitions |

The Multiwfn examples are redistributed under Multiwfn's licence for testing only.

Regenerate (neither is run by the suite):

* the Gaussian job: `g09.exe inputs/h2co_td_augccpvdz.gjf` with `GAUSS_EXEDIR=C:\G09W`, then
  `formchk`, `gzip -9` the `.fchk`;
* the references: `python tests/data/excitation/make_reference.py <Multiwfn.exe> [fixture ...] [--states-only]`.

## What the references showed

* **Multiwfn's default grid box (atoms +/- 6 bohr) truncates diffuse electrons**: with it,
  N-phenylpyrrole S4's electron integrates to 0.955 and the H2CO Rydberg S3's to 0.802.
  The references of those two fixtures use 10 and 14 bohr. This program's analytic moments
  have no box at all.
* **Lambda converges slowly with the grid** (|phi| has cusps at the orbital nodes): Multiwfn's
  automatic 15 x 170 points per atom differ from the 75 x 434 result by up to ~6e-3; on the
  same 75 x 434 size the two programs agree to <= 3.5e-4.

## Tolerances

Grid quantities carry Multiwfn's grid error. The largest difference between each
reference's fine grid and its 1.5x coarser grid, over all fixtures and states:

| Quantity | max fine - coarse |
|---|---|
| integral of hole / electron | 1.1e-3 / 2.3e-4 |
| Sm / Sr | 1.0e-3 / 4.6e-3 |
| HDI / EDI | 0.10 / 0.11 |
| D index | 1.0e-3 A (printed to 3 decimals) |
| centroids | 2.5e-4 A |
| transition dipole | 5.5e-4 a.u. |

The tests compare the exact (analytic) or Becke-grid values with the fine-grid reference
within a few times these, Multiwfn's |w| < 0.01 cross-term screening reproduced
(`cross_term_threshold=0.01`). Orbital quantities computed the same way in both programs
(Delta_r, MO contributions, major transitions) are compared to the printed digits.
