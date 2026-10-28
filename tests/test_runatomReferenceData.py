"""The stored reference outputs are complete and came from the documented build."""

import pytest

from runatom_fixtures import REFERENCE_CASES, REFERENCE_DIR, case_info, reference_names, reference_text


def test_manifest_records_compiler_flags_and_patch():
    manifest = (REFERENCE_DIR / "MANIFEST.txt").read_text()
    assert "gfortran" in manifest and "-O1 -std=legacy" in manifest and "idate" in manifest


@pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
def test_every_case_has_ae_pseudo_and_test_outputs(name):
    info = case_info(name)
    names = reference_names(name)
    atom, cfg = info["atom"], info["cfg"]
    assert f"{atom}_ae_{cfg}.out" in names and f"{atom}_ps.out" in names and f"{atom}_ps.uspp" in names
    assert f"{atom}_ps_formatted.uspp" in names       # reform.x itype=1: the CPMD file
    for test in info["tests"].split():
        assert f"{atom}_test_{test}.out" in names
    assert "writing the pseudopotential" not in reference_text(name, f"{atom}_ps.out")  # stderr text, not .out
    assert "comparison of all-electron and pseudo eigenvalues" in reference_text(name, f"{atom}_ps.out")
