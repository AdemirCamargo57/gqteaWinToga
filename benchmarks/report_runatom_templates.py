"""Report how the same-group templates of the USPP generator fare (dev-time, not a test).

Usage (from venv/src):  python benchmarks/report_runatom_templates.py [runatom.exe]

For every element up to U without a library case: the donor case, whether the
projectors could be carried over, the problems runatom's own checks find, and,
when there are none, the outcome of the all-electron and generation stages with
the shipped runatom.exe. The table is quoted in USER_MANUAL.md, so users know
which templates generate as copied.
"""

import sys
import tempfile
from pathlib import Path

SRC = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SRC))

import runatomDefaults as rd   # noqa: E402
import runatomInput as ri      # noqa: E402
import runatomRunner as rr     # noqa: E402


def main():
    exe = rr.find_runatom(sys.argv[1] if len(sys.argv) > 1 else None)
    roots = rr.library_dirs()
    rows = []
    for z in range(1, rd.MAX_DEFAULTS_Z + 1):
        symbol = rd.element_symbol(z)
        option = rd.suggestions(symbol, roots)[0][0]
        if option.kind == "library":
            continue
        donor = option.folder.name if option.folder else "-"
        projectors = "yes" if option.fields["ps_betas"] else "no"
        try:
            case = ri.case_from_fields(option.fields)
            problems = ri.check_case(case)
        except ValueError as exc:
            case, problems = None, [str(exc)]
        if problems:
            outcome = "not run: " + problems[0][:70]
        elif exe is None:
            outcome = "not run: runatom.exe not found"
        else:
            with tempfile.TemporaryDirectory() as tmp:
                ri.write_case(case, tmp)
                try:
                    rr.run_job(rr.runatom_command(exe), rr.plan_stages(case)[:2], tmp, timeout=60)
                    outcome = "generated"
                except rr.RunatomError as exc:
                    outcome = f"stopped: {exc.stage_label}: {exc.message[:60]}"
        rows.append((symbol, option.kind, donor, projectors, outcome))
    width = max(len(r[2]) for r in rows)
    print(f"{'El':3} {'kind':9} {'donor':{width}} {'proj':4} outcome")
    for symbol, kind, donor, projectors, outcome in rows:
        print(f"{symbol:3} {kind:9} {donor:{width}} {projectors:4} {outcome}")
    generated = sum(r[4] == "generated" for r in rows)
    print(f"\n{len(rows)} elements without a library case; {generated} templates generate as copied.")


if __name__ == "__main__":
    main()
