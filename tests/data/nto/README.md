# NTO test fixtures

Used by `tests/test_excitationData.py`, `tests/test_ntoAnalysis.py` and
`tests/test_ntoTab.py` through `tests/nto_fixtures.py`.

| Files | Origin | Covers |
|---|---|---|
| `uracil.fch`, `uracil.out`, `inputs/uracil.gjf` | Multiwfn 2026.7.15 `examples\excit\NTO` (Gaussian 09 D.01, PBE0/6-31G*, TD, 3 states, `IOp(9/40=4)`). Redistributed under Multiwfn's licence; the manual (section 4.18.6) prints the S3 weights 0.865529, 0.134025, 0.000582, 0.000121, 0.000063, ... | the published reference case |
| `h2co_td5050_631gd.*` | G09W A.02, B3LYP/6-31G(d) `TD(50-50,NStates=4)` | singlets and triplets in one file, de-excitations (`<-`) |
| `h2co_cis_ccpvdz.*` | G09W A.02, `CIS(NStates=3)/cc-pVDZ` | CIS (no de-excitations), pure d shells |
| `ch3_utd_631gd.*` | G09W A.02, UB3LYP/6-31G(d) `TD(NStates=3)` | open shell: `A`/`B` labels, separate alpha and beta NTOs, `2.004-E"` state labels |
| `benzene_td_631gd.*` | G09W A.02, B3LYP/6-31G(d) `TD(NStates=4)` | degenerate NTO pairs (e1g -> e2u), `Singlet-?Sym` labels |
| `*_state<n>_multiwfn.fch.gz` | Multiwfn 2026.7.15, main function 18, subfunction 6, option 2 (`make_reference.py`) | reference NTOs; orbital energies = NTO eigenvalues |

Every Gaussian input carries `IOp(9/40=4)` (all configuration coefficients
above 1e-4 printed), as the Multiwfn manual requires for accurate NTOs.

Regenerate:

* the Gaussian files: `python tests/data/nto/run_gaussian.py C:\G09W\g09.exe`
  (only when an input changes);
* the Multiwfn references: `python tests/data/nto/make_reference.py <Multiwfn.exe>`.

Neither script is run by the test suite. Note that when there are more than 10
virtual orbitals Multiwfn's printed "Sum of all eigenvalues" is the sum of the
**10 largest** eigenvalues only (`excittrans.f90`, `subroutine NTO`).
