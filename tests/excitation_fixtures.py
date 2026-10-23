"""Fixtures of the excitation-analysis tests: .fchk + Gaussian output pairs and Multiwfn references."""

import gzip
import re
import shutil
import tempfile
from pathlib import Path
from typing import Dict, Tuple

import pytest

DATA = Path(__file__).resolve().parent / "data"
EXCITATION_DIR = DATA / "excitation"
NTO_DIR = DATA / "nto"

# name -> (folder, formatted checkpoint (maybe .gz), Gaussian output)
FIXTURES = {
    "D-pi-A": (EXCITATION_DIR, "D-pi-A.fchk.gz", "D-pi-A.out"),                    # CAM-B3LYP/6-31G(d), LE + CT
    "N-phenylpyrrole": (EXCITATION_DIR, "N-phenylpyrrole.fch.gz", "N-phenylpyrrole.out"),  # 6-31+G(d), diffuse
    "4-nitroaniline": (EXCITATION_DIR, "4-nitroaniline.fch.gz", "4-nitroaniline.out"),    # PBE0/6-311G*, IFCT
    "h2co_ryd": (EXCITATION_DIR, "h2co_td_augccpvdz.fchk.gz", "h2co_td_augccpvdz.log"),  # aug-cc-pVDZ, Rydberg
    "h2co_td5050": (NTO_DIR, "h2co_td5050_631gd.fchk", "h2co_td5050_631gd.log"),     # triplets, de-excitations
    "h2co_cis": (NTO_DIR, "h2co_cis_ccpvdz.fchk", "h2co_cis_ccpvdz.log"),            # CIS, no de-excitations
    "ch3": (NTO_DIR, "ch3_utd_631gd.fchk", "ch3_utd_631gd.log"),                      # open shell
    "uracil": (NTO_DIR, "uracil.fch", "uracil.out"),                                  # Multiwfn NTO example
}

# Hole-electron references written by make_reference.py:
# name -> (states, grid spacing in bohr, grid extension beyond the atoms in bohr).
# Multiwfn's default extension (6 bohr) cuts off diffuse electrons (N-phenylpyrrole S4 integrates to 0.955,
# the H2CO Rydberg S3 to 0.802), so the diffuse fixtures use a wider box.
REFERENCE_STATES = {
    "D-pi-A": ((1, 2, 3), 0.12, 6.0),
    "N-phenylpyrrole": ((1, 4), 0.15, 10.0),
    "4-nitroaniline": ((2,), 0.12, 6.0),
    "h2co_ryd": ((1, 2, 3, 4, 5, 6), 0.15, 14.0),
    "h2co_td5050": ((1, 2, 6), 0.08, 6.0),
    "h2co_cis": ((1, 3), 0.08, 6.0),
    "ch3": ((1, 3), 0.08, 6.0),
}

# IFCT fragments (1-based atom ranges, ';' between fragments) as in the Multiwfn manual
FRAGMENTS = {
    "4-nitroaniline": "11-13; 1-10; 14-16",   # amino, benzene, nitro (manual 4.18.8.1)
    "D-pi-A": "21-23; 11-20; 1-10; 24-26",    # nitro, benzene, benzene, amino (manual 4.18.1)
}

_TMP = None


def _tmp_dir() -> Path:
    global _TMP
    if _TMP is None:
        _TMP = Path(tempfile.mkdtemp(prefix="excitation_fixtures_"))
    return _TMP


def fixture_paths(name: str) -> Tuple[Path, Path]:
    """(fchk path, Gaussian output path); a .gz fchk is decompressed once per session."""
    folder, fchk, log = FIXTURES[name]
    source = folder / fchk
    if source.suffix == ".gz":
        target = _tmp_dir() / source.stem
        if not target.is_file():
            with gzip.open(source, "rb") as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst)
        source = target
    return source, folder / log


def reference_path(name: str, state=None) -> Path:
    suffix = f"_state{state}" if state is not None else "_states"
    return EXCITATION_DIR / f"{name}{suffix}_multiwfn.txt.gz"


def multiwfn_text(name: str, state=None) -> str:
    """Multiwfn's stored output for (fixture, state); state None = the all-state run (Δr, Λ, transitions)."""
    path = reference_path(name, state)
    if not path.is_file():
        pytest.skip(f"no Multiwfn reference {path.name} (run tests/data/excitation/make_reference.py)")
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def sections(text: str) -> Dict[str, str]:
    """The '### <name>' sections of a reference file."""
    out, name = {}, None
    for line in text.splitlines():
        if line.startswith("### "):
            name = line[4:].strip()
            out[name] = []
        elif name is not None:
            out[name].append(line)
    return {k: "\n".join(v) for k, v in out.items()}


_NUM = r"(-?\d+\.\d+)"


def _floats(pattern: str, text: str):
    m = re.search(pattern, text)
    return None if m is None else [float(v) for v in m.groups()]


def parse_hole_electron(text: str) -> Dict[str, object]:
    """Printed indices of option 1 (hole_electron), Angstrom where Multiwfn prints Angstrom."""
    one = lambda pat: (_floats(pat, text) or [None])[0]
    return {
        "hole_integral": one(r"Integral of hole:\s+" + _NUM),
        "electron_integral": one(r"Integral of electron:\s*" + _NUM),
        "transition_integral": one(r"Integral of transition density:\s+" + _NUM),
        "transition_dipole": _floats(r"Transition dipole moment in X/Y/Z:\s*" + r"\s*".join([_NUM] * 3), text),
        "sm": one(r"Sm index \(integral of Sm function\):\s+" + _NUM),
        "sr": one(r"Sr index \(integral of Sr function\):\s+" + _NUM),
        "hole_centroid": _floats(r"Centroid of hole in X/Y/Z:\s+" + r"\s+".join([_NUM] * 3), text),
        "electron_centroid": _floats(r"Centroid of electron in X/Y/Z:\s+" + r"\s+".join([_NUM] * 3), text),
        "d": one(r"D index:\s+" + _NUM),
        "dipole_change": _floats(r"X:\s+" + _NUM + r"\s+Y:\s+" + _NUM + r"\s+Z:\s+" + _NUM + r"\s+Norm", text),
        "sigma_hole": _floats(r"RMSD of hole in X/Y/Z:\s+" + r"\s+".join([_NUM] * 3) + r"\s+Norm:\s*" + _NUM, text),
        "sigma_electron": _floats(r"RMSD of electron in X/Y/Z:\s*" + r"\s+".join([_NUM] * 3) + r"\s+Norm:\s*" + _NUM,
                                  text),
        "delta_sigma": _floats(r"X:\s*" + _NUM + r"\s+Y:\s*" + _NUM + r"\s+Z:\s*" + _NUM + r"\s+Overall:\s*" + _NUM, text),
        "h": _floats(r"H_x:\s*" + _NUM + r"\s+H_y:\s*" + _NUM + r"\s+H_z:\s*" + _NUM + r"\s+H_CT:\s*" + _NUM
                     + r"\s+H index:\s*" + _NUM, text),
        "t": one(r"t index:\s*" + _NUM),
        "hdi": one(r"\(HDI\):\s+" + _NUM),
        "edi": one(r"\(EDI\):\s*" + _NUM),
        "ghost": _floats(r"Ghost-hunter index:\s*" + _NUM + r" eV, 1st term:\s*" + _NUM + r" eV, 2nd term:\s*"
                         + _NUM, text),
        "coulomb": one(r"Coulomb attractive energy:\s*" + _NUM),
    }


def parse_mo_contributions(text: str) -> Dict[Tuple[int, str], Tuple[float, float]]:
    """{(MO, spin label '' / 'A' / 'B'): (hole %, electron %)} of option 2."""
    out = {}
    for m in re.finditer(r"MO\s+(\d+)([AB]?), Occ:\s*\S+\s+Hole:\s*(-?\d+\.\d+) %\s+Electron:\s*(-?\d+\.\d+) %", text):
        out[(int(m.group(1)), m.group(2))] = (float(m.group(3)), float(m.group(4)))
    return out


def parse_atom_contributions(text: str) -> Dict[int, Tuple[float, float, float, float]]:
    """{atom (1-based): (hole %, electron %, overlap %, diff %)} - the last table printed (all atoms)."""
    out = {}
    pat = re.compile(r"^\s*(\d+)\((\S+)\s*\)\s+Hole:\s*(-?\d+\.\d+) %\s+Electron:\s*(-?\d+\.\d+) %\s+Overlap:\s*"
                     r"(-?\d+\.\d+) %\s+Diff\.:\s*(-?\d+\.\d+) %")
    for line in text.splitlines():
        m = pat.match(line)
        if m:
            out[int(m.group(1))] = tuple(float(m.group(k)) for k in range(3, 7))
    return out


def parse_ifct(text: str) -> Dict[str, object]:
    frag = {int(m.group(1)): (float(m.group(2)), float(m.group(3))) for m in
            re.finditer(r"^\s*(\d+)\s+Hole:\s*(-?\d+\.\d+) %\s+Electron:\s*(-?\d+\.\d+) %\s*$", text, re.M)}
    change = {int(m.group(1)): float(m.group(2)) for m in
              re.finditer(r"Variation of population number of fragment\s*(\d+):\s*(-?\d+\.\d+)", text)}
    redis = {int(m.group(1)): float(m.group(2)) for m in
             re.finditer(r"Intrafragment electron redistribution of fragment\s*(\d+):\s*(-?\d+\.\d+)", text)}
    transfer = {}
    for m in re.finditer(r"(\d+)\s+->\s*(\d+):\s*(-?\d+\.\d+)\s+\d+\s+<-\s*\d+:\s*(-?\d+\.\d+)", text):
        i, j = int(m.group(1)), int(m.group(2))
        transfer[(i, j)], transfer[(j, i)] = float(m.group(3)), float(m.group(4))
    one = lambda pat: (_floats(pat, text) or [None])[0]
    return {"fragments": frag, "population_change": change, "redistribution": redis, "transfer": transfer,
            "intrinsic_ct": one(r"Intrinsic charge transfer percentage, CT\(%\):\s*" + _NUM),
            "intrinsic_le": one(r"Intrinsic local excitation percentage, LE\(%\):\s*" + _NUM),
            "apparent_ct": one(r"Apparent charge transfer percentage, CT\(%\):\s*" + _NUM)}


def parse_delta_r(text: str) -> Dict[int, float]:
    """{state: Δr in Angstrom}."""
    return {int(m.group(1)): float(m.group(2)) for m in
            re.finditer(r"Excited state\s+(\d+):\s+Delta_r =\s*\S+ Bohr,\s*(-?\d+\.\d+) Angstrom", text)}


def parse_lambda(text: str) -> Dict[int, float]:
    return {int(m.group(1)): float(m.group(2)) for m in
            re.finditer(r"Excited state\s+(\d+):\s+lambda =\s*(-?\d+\.\d+)", text)}


def parse_major_transitions(text: str) -> Dict[int, str]:
    """{state: the line of MO transitions under its '#  n ...' header} (option 15)."""
    out, lines = {}, text.splitlines()
    for k, line in enumerate(lines):
        m = re.match(r"^\s*#\s*(\d+)\s+\S+ eV", line)
        if m and k + 1 < len(lines):
            out[int(m.group(1))] = lines[k + 1].strip()
    return out
