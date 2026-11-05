# Gaussian 16 output fixtures

Used by `tests/test_gaussianOutput.py`, `tests/test_gaussianInput.py` and `tests/test_gaussianInputBuilder.py`
(Load geometry from a `.log`/`.out`), and by the end-to-end test that the Excitation tab accepts the files of the
builder's *Geometry optimization + excited states* job. All were run with Gaussian 16 Rev. C.01 in WSL
(2026-11-05); regenerate them only with `run_gaussian.py` (never run by the suite).

| File | Input (`inputs/`) | Covers |
|---|---|---|
| `h2o_opt.log` | HF/3-21G `Opt`, water from 1.00 A O-H | one converged optimization |
| `h2o_opt_freq.out` | HF/3-21G `Opt Freq` | `.out` suffix; the internal Freq job also prints "Optimization completed" and must be ignored |
| `h2co_opt_maxcycles.log` | HF/3-21G `Opt=(MaxCycles=2)` | "Number of steps exceeded", error termination: refused |
| `h2o_opt_then_failed.log` | `Opt`, then `--Link1--` UHF cation `Opt=(MaxCycles=1) Geom=Check` | a later failed optimization: the earlier geometry is imported with a note |
| `h2co_opt_td.log`, `h2co_opt_td_excited.fchk.gz` | the builder's own job: B3LYP/6-31G(d) `Opt`, then `TD(NStates=10)` with `Geom=AllCheck Guess=Read IOp(9/40=4)`, `%OldChk`/`%Chk=h2co_opt_td_excited.chk` | two-step output; the Excitation tab reads the `_excited.fchk` with this log |
| (input only) `h2co_opt_td_gen.gjf` | `B3LYP/Gen Opt`, then `B3LYP ChkBasis TD` | `ChkBasis` runs in G16 (2 Normal terminations); the output is not kept |
