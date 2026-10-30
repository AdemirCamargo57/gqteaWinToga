"""Regenerate runatomStarterData.py: runatom's converged eigenvalues of the starter atoms.

Usage (from venv/src):  python tools/make_runatom_starter_data.py [path\\to\\runatom.exe]

Doc/INPUT_AE asks for eigenvalue guesses "within about a factor of 2-3" and
suggests "results from a nearby atomic number as a guide". Slater's rules miss
open d and f shells by far more (Nd 4f: -12 Ry against -0.13 Ry), and runatom
then fails to converge. So the guesses are bootstrapped through the table:
atom Z starts from the converged eigenvalues of atom Z - 1 for the shells they
share and from Slater's rules for the others, and runs with the shipped
runatom.exe. When the starter's own settings (runatomDefaults.STARTER_RMAX,
STARTER_DAMP, STARTER_MAXIT) do not converge, gentler settings are tried only to
obtain the eigenvalues; tests/test_runatomStarterExe.py then checks that the
starter, with its own settings and these guesses, converges for every atom.
"""

import re
import sys
import tempfile
from pathlib import Path

SRC = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SRC))

import runatomDefaults as rd   # noqa: E402
import runatomInput as ri      # noqa: E402
import runatomRunner as rr     # noqa: E402

OUTPUT = SRC / "runatomStarterData.py"
# (rmax, bbsf, damp, maxit); a setting whose mesh exceeds runatom's 1000 points is skipped
FALLBACK_SETTINGS = [(400.0, 40.0, 0.1, 2000), (400.0, 40.0, 0.2, 2000), (400.0, 60.0, 0.1, 2000),
                     (400.0, 60.0, 0.2, 2000), (400.0, 40.0, 0.3, 2000), (400.0, 40.0, 0.5, 2000),
                     (400.0, 60.0, 0.5, 2000), (400.0, 40.0, 0.05, 4000)]
EIGENVALUE_ROW = re.compile(r"^\s+(\d[0-3]0)\s+([-\d.]+)\s+([-\d.]+)\s*$", re.M)
HEADER = '''"""Converged runatom eigenvalues of the computed starter atoms (GENERATED).

Do not edit: regenerate with tools/make_runatom_starter_data.py (needs runatom.exe).
Atoms Z = 1..92 in the starter's configuration (runatomDefaults; neutral except
STARTER_XION), PBE, scalar-relativistic, runatom 7.3.6. {z: {nnlz: eigenvalue in
Ry}}; used as the eigenvalue guesses (ee) of the computed starter.
"""

'''


def converged_eigenvalues(command, ae: ri.AEInput) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "starter_ae.adat").write_text(ri.format_adat(ae), encoding="ascii", newline="\n")
        rr.run_stage(command, rr.Stage("ae", "All-electron", "starter_ae", "starter_ae.ae", "dummy"), tmp)
        text = (Path(tmp) / "starter_ae.out").read_text(encoding="latin-1")
    table = text.split("the eigenvalues are")[-1].split("etot")[0]
    return {int(nnlz): float(energy) for nnlz, _, energy in EIGENVALUE_ROW.findall(table)}


def main():
    exe = rr.find_runatom(sys.argv[1] if len(sys.argv) > 1 else None)
    if exe is None:
        sys.exit("runatom.exe not found")
    command = rr.runatom_command(exe)
    table, used, ions, previous, unconverged = {}, {}, {}, {}, []
    for z in range(1, rd.MAX_DEFAULTS_Z + 1):
        symbol = rd.element_symbol(z)
        settings = [(rd.STARTER_RMAX, rd.STARTER_BBSF, rd.STARTER_DAMP, rd.STARTER_MAXIT)] + FALLBACK_SETTINGS
        failed, values = [], None
        # the neutral atom first; the +1 ion only when no setting converges it
        # (Dy, Ho, Er: the 4f state never converges in the neutral atom)
        for xion in (0.0, 1.0):
            for rmax, bbsf, damp, maxit in settings:
                ae = rd.starter_ae_input(symbol, guesses=previous, rmax=rmax, bbsf=bbsf, damp=damp,
                                         maxit=maxit, xion=xion)
                if ri.check_ae(ae):
                    continue
                try:
                    values = converged_eigenvalues(command, ae)
                    break
                except rr.RunatomError as exc:
                    failed.append(f"xion {xion:g}, rmax {rmax:g}, bbsf {bbsf:g}, damp {damp:g}: {exc}")
            if values is not None:
                break
        if values is None:
            print(f"{symbol}: NOT CONVERGED\n  " + "\n  ".join(failed))
            unconverged.append(symbol)
            continue
        missing = {s.nnlz for s in ae.shells} - set(values)
        if missing:
            sys.exit(f"{symbol}: no eigenvalue read for shells {sorted(missing)}")
        table[z] = {s.nnlz: round(values[s.nnlz], 4) for s in ae.shells}
        used[z] = (rmax, bbsf, damp, maxit)
        if xion:
            ions[z] = xion
        previous = table[z]
        print(f"{symbol}: converged (xion {xion:g}, rmax {rmax:g}, bbsf {bbsf:g}, damp {damp:g})")
    lines = [HEADER, "STARTER_EIGENVALUES = {"]
    for z, values in table.items():
        lines.append(f"    {z}: {{" + ", ".join(f"{k}: {v!r}" for k, v in values.items()) + "},"
                     f"   # {rd.element_symbol(z)}")
    lines += ["}", "", "# (rmax, bbsf, damp, maxit) where the starter's defaults do not converge",
              "STARTER_SETTINGS = {"]
    default = (rd.STARTER_RMAX, rd.STARTER_BBSF, rd.STARTER_DAMP, rd.STARTER_MAXIT)
    lines += [f"    {z}: {s!r},   # {rd.element_symbol(z)}" for z, s in used.items() if s != default]
    lines += ["}", "", "# atoms whose neutral configuration runatom does not converge: the starter",
              "# uses this ion charge (electrons taken from the outermost s shell)",
              "STARTER_XION = {"]
    lines += [f"    {z}: {x!r},   # {rd.element_symbol(z)}" for z, x in ions.items()]
    lines += ["}", "", "# atoms runatom did not converge with any setting tried",
              f"UNCONVERGED = {tuple(unconverged)!r}"]
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
