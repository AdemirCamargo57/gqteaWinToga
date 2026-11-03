"""Dev-time check of the Gaussian 16 input builder against the real program (never run by the suite).

Writes one small input per job type into an empty folder:
    python tools/check_gaussian_inputs.py <folder>
then, in WSL with Gaussian 16 set up, in that folder:
    for f in *.gjf; do g16 < "$f" > "${f%.gjf}.log"; done
    grep -c "Normal termination" *.log
Each log must count one Normal termination per step (the Fukui job: 3).
"""
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import gaussianInput as gi           # noqa: E402
import gaussianKeywords as gk        # noqa: E402

WATER = ("O 0.000000 0.000000 0.117300\n"
         "H 0.000000 0.757200 -0.469200\n"
         "H 0.000000 -0.757200 -0.469200")
HCN_TS_GUESS = ("C 0.000000 0.000000 0.000000\n"
                "N 0.000000 0.000000 1.180000\n"
                "H 1.100000 0.000000 0.600000")
BASE = gi.GaussianStepSpec(family=gk.FAMILY_HF, method_name=gk.NO_CHOICE, basis="STO-3G", multiplicity=1,
                           memory="1GB", nproc=2)


def step(**changes):
    return replace(BASE, **changes)


TS = step(job_type=gk.JOB_OPT, opt_target=gk.OPT_TS, force_constants=gk.FC_CALCFC,
          additional_keywords="Opt=NoEigenTest")
# IRC=RCFC reads the Cartesian force constants a Freq run stores; after an Opt-only TS step it stops in l123.
TS_FREQ = replace(TS, job_type=gk.JOB_OPT_FREQ)

JOBS = {
    "energy": (WATER, [step()], True),
    "optfreq": (WATER, [step(job_type=gk.JOB_OPT_FREQ, freq_raman=True)], False),
    # NoSymm: with symmetry on, the scanned O-H bond is tied to the other and the scan never converges.
    "scan": (WATER, [step(job_type=gk.JOB_SCAN, modredundant="B 1 2 S 3 0.05", nosymm=True)], False),
    "freeze": (WATER, [step(job_type=gk.JOB_OPT, modredundant="A 2 1 3 F")], False),
    "stability": (WATER, [step(job_type=gk.JOB_STABILITY, stable_opt=True)], False),
    "nmr": (WATER, [step(job_type=gk.JOB_NMR, nmr_spinspin=True)], False),
    "cation_uhf": (WATER, [step(charge=1, multiplicity=2, spin="Unrestricted (U)")], False),
    "cis": (WATER, [step(state=gk.STATE_CIS, n_states=3)], False),
    "pm6": (WATER, [step(family=gk.FAMILY_SEMI, method_name="PM6")], False),
    "cbs4m": (WATER, [step(family=gk.FAMILY_COMPOUND, method_name="CBS-4M")], False),
    "dft_smd": (WATER, [step(family=gk.FAMILY_DFT, method_name="B3LYP", basis="6-31G", basis_polarization="(d)",
                             dispersion="GD3BJ", solvation="SMD", solvent="Water", pop_mk=True, scf_tight=True,
                             nosymm=True, print_level="#N")], False),
    "linked_td": (WATER, [step(job_type=gk.JOB_OPT),
                          replace(gi.linked_step(step(job_type=gk.JOB_OPT)), state=gk.STATE_TD, n_states=3,
                                  print_all_coefficients=True)], True),
    "ts_irc": (HCN_TS_GUESS, [TS, replace(gi.linked_step(TS), job_type=gk.JOB_IRC, irc_max_points=5)], False),
    "irc_rcfc": (HCN_TS_GUESS, [TS_FREQ, replace(gi.linked_step(TS_FREQ), job_type=gk.JOB_IRC,
                                                 irc_force_constants=gk.FC_RCFC, irc_max_points=5)], False),
    "fukui": (WATER, [gi.fukui_defaults(step(fukui_plus_multiplicity=2, fukui_minus_multiplicity=2))], True),
}


def main(folder: str) -> None:
    out = Path(folder)
    out.mkdir(parents=True, exist_ok=True)
    builder = gi.GaussianInputBuilder()
    for name, (geometry, steps, formchk) in JOBS.items():
        job = builder.prepare(gi.GaussianJobSpec(stem=name, output_dir=str(out), geometry=geometry, steps=steps,
                                                 write_formchk=formchk))
        builder.write_files(job)
        print(f"{job.input_path.name}: {len(job.steps) or 3} step(s)"
              + "".join(f"\n  WARNING {w}" for w in job.warnings))


if __name__ == "__main__":
    main(sys.argv[1])
