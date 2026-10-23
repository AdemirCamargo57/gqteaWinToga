"""Regenerate the Multiwfn references for the excitation-analysis tests.

NOT run by the test suite (it needs Multiwfn). Run it by hand from the source
directory with the path of Multiwfn.exe (verified with Multiwfn 2026.7.15):

    python tests/data/excitation/make_reference.py path\\to\\Multiwfn.exe [fixture ...]

For every (fixture, state) of excitation_fixtures.REFERENCE_STATES it drives
Multiwfn's main function 18 through stdin and stores the printed output as
<fixture>_state<n>_multiwfn.txt.gz, one '### <section>' per run:

* indices / indices_coarse - subfunction 1, option 1, uniform grid of the
  fixture's spacing / 1.5 x that spacing (their difference estimates the grid error)
  and the fixture's box extension (Multiwfn's 6 bohr truncates diffuse electrons);
* coulomb - option 1 on the low-quality grid, then post-processing option 18;
* mo - option 2 (MO contributions, threshold 1e-4 %);
* mulliken / hirshfeld - option 3 with partition 1 / 2, hydrogens included (7);
* ifct_mulliken / ifct_hirshfeld - subfunction 8 for fixtures in FRAGMENTS.

Per fixture, <fixture>_states_multiwfn.txt.gz holds subfunction 4 (delta_r),
14 (lambda) and 15 (major MO transitions, booted from the Gaussian output) for
all states.
"""

import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))  # tests/, for excitation_fixtures

from excitation_fixtures import FIXTURES, FRAGMENTS, REFERENCE_STATES, fixture_paths  # noqa: E402

TIMEOUT = 7200


def _run(exe, boot, script, work):
    proc = subprocess.run([exe, boot], input=script, text=True, cwd=work, capture_output=True, timeout=TIMEOUT)
    return proc.stdout


def _write(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=9) as handle:
        handle.write(text)


def _n_states(log: Path) -> int:
    text = log.read_text(encoding="utf-8", errors="replace")
    block = text.rsplit("Excitation energies and oscillator strengths:", 1)[-1]
    return len(re.findall(r"^\s*Excited State\s+\d+:", block, re.M))


def _fragment_lines(text: str):
    parts = [p.strip().replace(" ", "") for p in text.split(";")]
    return [str(len(parts))] + parts


def main(exe: str, only=(), states_only=False) -> None:
    settings = Path(exe).resolve().parent / "settings.ini"
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        text = settings.read_text(encoding="utf-8", errors="replace")
        text = re.sub(r"nthreads=\s*\d+", f"nthreads= {os.cpu_count() or 4}", text)
        (work / "settings.ini").write_text(text, encoding="utf-8")
        for name in FIXTURES:
            if only and name not in only:
                continue
            fchk_src, log_src = fixture_paths(name)
            fchk, log = fchk_src.name, log_src.name
            shutil.copy(fchk_src, work / fchk)
            shutil.copy(log_src, work / log)
            if name in REFERENCE_STATES and not states_only:
                states, spacing, extension = REFERENCE_STATES[name]
                for state in states:
                    head = f"18\n1\n{log}\n{state}\n"
                    box = f"1\n-10\n{extension:g}\n"  # option 1, then the grid menu's extension distance
                    runs = {
                        "indices": head + box + f"4\n{spacing:g}\n0\n0\n0\nq\n",
                        "indices_coarse": head + box + f"4\n{1.5 * spacing:g}\n0\n0\n0\nq\n",
                        "coulomb": head + box + "1\n18\n0\n0\n0\nq\n",
                        "mo": head + "2\n0.0001\n0\n0\nq\n",
                        "mulliken": head + "3\n1\n7\n0\n0\n0\nq\n",
                        "hirshfeld": head + "3\n2\n7\n0\n0\n0\nq\n",
                    }
                    if name in FRAGMENTS:
                        frags = "\n".join(_fragment_lines(FRAGMENTS[name]))
                        for method, key in ((1, "ifct_mulliken"), (2, "ifct_hirshfeld")):
                            runs[key] = f"18\n8\n{method}\n{log}\n{state}\n{frags}\n0\n0\nq\n"
                    out = []
                    for key, script in runs.items():
                        print(f"{name} state {state}: {key}...", flush=True)
                        out.append(f"### {key}\n" + _run(exe, fchk, script, work))
                    _write(HERE / f"{name}_state{state}_multiwfn.txt.gz", "\n".join(out))
            n = _n_states(work / log)
            lam = f"18\n14\n{log}\n1-{n}\nn\n0\nq\n"
            out = [f"### delta_r\n" + _run(exe, fchk, f"18\n4\n{log}\n1-{n}\n0\nq\n", work),
                   f"### lambda\n" + _run(exe, fchk, lam, work),
                   f"### major\n" + _run(exe, log, "18\n15\nn\n0\nq\n", work)]
            # Lambda again on this program's Becke grid size (75 x 434) instead of Multiwfn's automatic 15 x 170
            (work / "settings.ini").write_text(re.sub(r"iautointgrid=\s*\d", "iautointgrid= 0", text), encoding="utf-8")
            out.append(f"### lambda_75x434\n" + _run(exe, fchk, lam, work))
            (work / "settings.ini").write_text(text, encoding="utf-8")
            _write(HERE / f"{name}_states_multiwfn.txt.gz", "\n".join(out))
            print(f"{name}: all-state references written", flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    args = [a for a in sys.argv[2:] if a != "--states-only"]
    main(sys.argv[1], tuple(args), states_only="--states-only" in sys.argv)
