"""The vendored USPP 7.3.6 tree is unmodified and the build patch still applies."""

import hashlib
from pathlib import Path

import pytest

VENDOR = Path(__file__).resolve().parents[1] / "third_party" / "uspp-736"
IDATE_LINE = "      external idate"


def _manifest():
    rows = (VENDOR / "SOURCE_SHA256.txt").read_text().splitlines()
    return [row.split("  ", 1) for row in rows if row.strip()]


def test_every_upstream_file_is_unmodified():
    manifest = _manifest()
    assert len(manifest) > 600
    for digest, name in manifest:
        assert hashlib.sha256((VENDOR / name).read_bytes()).hexdigest() == digest, name


def test_line_endings_are_protected_from_git():
    assert "* -text" in (VENDOR / ".gitattributes").read_text()


def test_the_idate_patch_targets_exactly_one_line():
    lines = (VENDOR / "Source" / "runatom.f").read_text(encoding="latin-1").splitlines()
    assert lines.count(IDATE_LINE) == 1
    patch = (VENDOR / "patches" / "0001-gfortran-intrinsic-idate.patch").read_text()
    assert f"-{IDATE_LINE}" in patch and "+c     external idate" in patch


@pytest.mark.parametrize("name", ["runatom.exe", "reform.exe"])
def test_built_program_matches_its_recorded_checksum(name):
    exe = VENDOR / "bin" / name
    if not exe.exists():
        pytest.skip(f"{name} has not been built yet (Task 2)")
    recorded = (VENDOR / "bin" / f"{name}.sha256").read_text().split()[0]
    assert hashlib.sha256(exe.read_bytes()).hexdigest() == recorded
