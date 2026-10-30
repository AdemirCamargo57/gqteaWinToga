"""The computed starter's all-electron input converges with the shipped runatom.exe (Z = 1..92).

Skipped without the exe. runatomStarterData holds the eigenvalue guesses and
the per-atom settings; this is what keeps them honest.
"""

import pytest

import runatomDefaults as rd
import runatomInput as ri
import runatomRunner as rr

EXE = rr.find_runatom()
pytestmark = pytest.mark.skipif(EXE is None, reason="runatom.exe not built")


def test_no_atom_is_left_unconverged():
    assert rd.UNCONVERGED == ()


@pytest.mark.parametrize("z", range(1, rd.MAX_DEFAULTS_Z + 1), ids=rd.element_symbol)
def test_the_starter_all_electron_stage_converges(z, tmp_path):
    symbol = rd.element_symbol(z)
    case = ri.case_from_fields(dict(rd.starter_fields(symbol), **_placeholder_ps()))
    stem = ri.file_stems(case)["ae"]
    (tmp_path / f"{stem}.adat").write_text(ri.format_adat(case.ae), encoding="ascii", newline="\n")
    result = rr.run_stage(rr.runatom_command(EXE), rr.Stage("ae", "All-electron", stem, f"{stem}.ae", "dummy"),
                          tmp_path, timeout=60)
    assert result.returncode == 0 and (tmp_path / f"{stem}.ae").is_file()


def _placeholder_ps() -> dict:
    """Values for the starter's blank pseudopotential fields, so the form parses; only the AE input runs."""
    return {"ps_bessel": "10 30 40 5", "ps_keyps": "3 0 1.0 0", "ps_rcloc": "1.0", "ps_rc": "1.0",
            "ps_betas": "0 1 0 0", "ps_local": "0 0 0 1", "ps_qfunc": "0 0 0"}
