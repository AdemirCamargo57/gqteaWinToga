"""Tests for the shared Gaussian utility finder (cubegen, formchk, ...)."""
import subprocess
import sys
from pathlib import Path

import pytest

import gaussianUtilities as gu


@pytest.fixture
def no_system_gaussian(monkeypatch):
    """Hide any real Gaussian install: no GAUSS_EXEDIR, nothing on PATH, no default folders."""
    monkeypatch.delenv("GAUSS_EXEDIR", raising=False)
    monkeypatch.setattr(gu.shutil, "which", lambda name: None)
    monkeypatch.setattr(gu, "DEFAULT_GAUSSIAN_FOLDERS", ())


def make_exe(folder: Path, name: str) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / name
    path.write_text("")
    return path


class TestFindGaussianUtility:
    def test_explicit_path_wins(self, tmp_path, no_system_gaussian):
        exe = make_exe(tmp_path, "formchk.exe")
        assert gu.find_gaussian_utility("formchk", str(exe)) == str(exe)

    def test_nothing_found_gives_none(self, tmp_path, no_system_gaussian):
        assert gu.find_gaussian_utility("formchk", str(tmp_path / "missing.exe")) is None

    def test_gauss_exedir_with_several_folders(self, tmp_path, monkeypatch, no_system_gaussian):
        exe = make_exe(tmp_path / "b", "formchk.exe")
        monkeypatch.setenv("GAUSS_EXEDIR", f"{tmp_path / 'a'}{gu.os.pathsep}{tmp_path / 'b'}")
        assert gu.find_gaussian_utility("formchk") == str(exe)

    def test_path_lookup(self, tmp_path, monkeypatch, no_system_gaussian):
        exe = make_exe(tmp_path, "formchk.exe")
        monkeypatch.setattr(gu.shutil, "which", lambda name: str(exe) if name == "formchk" else None)
        assert gu.find_gaussian_utility("formchk") == str(exe)

    def test_default_folders(self, tmp_path, monkeypatch, no_system_gaussian):
        exe = make_exe(tmp_path / "G16W", "cubegen.exe")
        monkeypatch.setattr(gu, "DEFAULT_GAUSSIAN_FOLDERS", (str(tmp_path / "G09W"), str(tmp_path / "G16W")))
        assert gu.find_gaussian_utility("cubegen") == str(exe)

    def test_cubegen_finder_delegates(self, tmp_path, no_system_gaussian):
        import electrostaticSurfacePotential as esp
        exe = make_exe(tmp_path, "cubegen.exe")
        assert esp.find_cubegen(str(exe)) == str(exe)
        assert esp.find_cubegen(str(tmp_path / "nope.exe")) is None


def test_module_does_not_import_toga():
    code = "import sys, gaussianUtilities; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(Path(__file__).resolve().parents[1]))
    assert out.stdout.strip() == "False", out.stderr
