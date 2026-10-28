"""The CPMD formatted pseudopotential: reform conversion and the readvan-order reader."""

import sys
from pathlib import Path

import numpy as np
import pytest

import runatomInput as ri
import runatomOutput as ro
import runatomRunner as rr
from runatom_fixtures import LIBRARY_DIR, REFERENCE_CASES, case_info, reference_bytes

FAKE = [sys.executable, "-I", str(Path(__file__).resolve().parent / "fake_reform.py")]
REL = 1e-11          # 1p4e19.11 keeps 12 significant digits
ABS = {"z": 5e-10, "zv": 5e-10, "exfact": 5e-10, "ptryc": 5e-10,    # written f15.9 by reform
       "eloc": 5e-6, "qtryc": 5e-6}                                  # written f9.5 by reform
SCALARS = ("version", "date", "title", "z", "zv", "exfact", "valence", "mesh", "etot", "keyps", "ifpcor",
           "nang", "lloc", "eloc", "ifqopt", "nqf", "qtryc", "irel", "kkbeta", "npf", "ptryc", "rcloc", "rpcor")
ARRAYS = ("rinner", "rc", "l", "eee", "beta", "ddd0", "ddd", "qqq", "qfunc", "qfcoef", "iptype",
          "vloc0", "rspsco", "vloc", "rsatom", "r", "rab", "wf")


def _files(name, tmp_path):
    atom = case_info(name)["atom"]
    binary, formatted = tmp_path / f"{atom}_ps.uspp", tmp_path / f"{atom}_ps_formatted.uspp"
    binary.write_bytes(reference_bytes(name, binary.name))
    return binary, formatted


def _assert_same(a: ro.UsppData, b: ro.UsppData):
    for attr in SCALARS:
        x, y = getattr(a, attr), getattr(b, attr)
        if isinstance(x, float):
            assert x == pytest.approx(y, rel=REL, abs=ABS.get(attr, 0.0)), attr
        elif attr == "valence":          # (i5,2f15.9)
            assert [(n, pytest.approx(w, abs=5e-10), pytest.approx(e, abs=5e-10)) for n, w, e in x] == y, attr
        else:
            assert x == y, attr
    for attr in ARRAYS:
        x, y = getattr(a, attr), getattr(b, attr)
        if x is None or y is None:
            assert x is None and y is None, attr
        else:
            assert np.allclose(x, y, rtol=REL, atol=0), attr


class TestReader:
    def test_fortran_reals_as_list_directed_input_reads_them(self):
        assert ro._fortran_float("1.23456789012-100") == pytest.approx(1.23456789012e-100)
        assert ro._fortran_float("-4.50000000000+200") == -4.5e200
        assert ro._fortran_float("1.0D-03") == 1e-3 and ro._fortran_float("-7") == -7.0

    @pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
    def test_reference_formatted_file_holds_the_binary_one(self, name, tmp_path):
        binary, formatted = _files(name, tmp_path)
        formatted.write_bytes(reference_bytes(name, formatted.name))
        _assert_same(ro.read_formatted_uspp(formatted), ro.read_uspp(binary))

    def test_truncated_file_is_refused(self, tmp_path):
        binary, formatted = _files("000-Ex-ca-sp-test", tmp_path)
        text = reference_bytes("000-Ex-ca-sp-test", formatted.name).decode("latin-1")
        formatted.write_text("\n".join(text.splitlines()[:40]) + "\n", encoding="latin-1")
        with pytest.raises(ValueError, match="ends early"):
            ro.read_formatted_uspp(formatted)


class TestConversion:
    def test_names_and_cpmd_line(self):
        case = ri.read_case_folder(LIBRARY_DIR / "000-Ex" / "000-Ex-ca-sp-test")
        assert rr.formatted_name(case) == "ti_ps_formatted.uspp"
        assert rr.cpmd_atom_line("ti_ps_formatted.uspp") == "*ti_ps_formatted.uspp FORMATTED"

    def test_conversion_replaces_an_old_file_and_writes_lf(self, tmp_path, monkeypatch):
        monkeypatch.setenv("FAKE_REFORM_MODE", "ok")
        (tmp_path / "ti_ps.uspp").write_bytes(b"binary")
        (tmp_path / "ti_ps_formatted.uspp").write_text("old")
        path = rr.write_formatted(FAKE, tmp_path, "ti_ps.uspp", "ti_ps_formatted.uspp")
        assert path.read_bytes() == b"    7    3    6\nformatted\n"

    @pytest.mark.parametrize("mode, text", [("fail", "No such file"), ("silent", "could not convert")])
    def test_failures_are_reported(self, tmp_path, monkeypatch, mode, text):
        monkeypatch.setenv("FAKE_REFORM_MODE", mode)
        (tmp_path / "ti_ps.uspp").write_bytes(b"binary")
        with pytest.raises(rr.RunatomError) as info:
            rr.write_formatted(FAKE, tmp_path, "ti_ps.uspp", "ti_ps_formatted.uspp")
        assert info.value.kind == "formatted" and text in str(info.value)
        assert not (tmp_path / "ti_ps_formatted.uspp").exists()

    def test_names_longer_than_reform_reads_are_refused(self, tmp_path):
        with pytest.raises(ValueError, match="36"):
            rr.write_formatted(FAKE, tmp_path, "a" * 37, "b.uspp")

    def test_frozen_build_finds_the_bundled_reform(self, tmp_path, monkeypatch):
        exe = tmp_path / rr.FROZEN_SUBDIR / "reform.exe"
        exe.parent.mkdir()
        exe.write_bytes(b"")
        monkeypatch.delenv(rr.REFORM_ENV_VAR, raising=False)
        monkeypatch.setattr(sys, "frozen", True, raising=False)
        monkeypatch.setattr(sys, "_MEIPASS", str(tmp_path), raising=False)
        assert rr.find_reform() == str(exe)


REFORM = rr.find_reform()


@pytest.mark.skipif(REFORM is None, reason="reform.exe is not available")
@pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
def test_reform_exe_reproduces_the_linux_reference_byte_for_byte(name, tmp_path):
    binary, formatted = _files(name, tmp_path)
    rr.write_formatted(rr.reform_command(REFORM), tmp_path, binary.name, formatted.name)
    assert formatted.read_bytes() == reference_bytes(name, formatted.name)
