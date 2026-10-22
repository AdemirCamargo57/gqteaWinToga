"""Fixtures of the NTO tests: Gaussian .fchk + output pairs and Multiwfn references."""

import gzip
import shutil
from pathlib import Path

import pytest

NTO_DIR = Path(__file__).resolve().parent / "data" / "nto"

# name -> (formatted checkpoint, Gaussian output)
FIXTURES = {
    "uracil": ("uracil.fch", "uracil.out"),                       # Multiwfn example, PBE0/6-31G* TD
    "h2co_td5050": ("h2co_td5050_631gd.fchk", "h2co_td5050_631gd.log"),  # singlets + triplets, de-excitations
    "h2co_cis": ("h2co_cis_ccpvdz.fchk", "h2co_cis_ccpvdz.log"),         # CIS, pure d
    "ch3": ("ch3_utd_631gd.fchk", "ch3_utd_631gd.log"),                   # open shell (A/B)
    "benzene": ("benzene_td_631gd.fchk", "benzene_td_631gd.log"),         # degenerate NTO pairs
}

# states exported by Multiwfn (make_reference.py)
REFERENCE_STATES = {
    "uracil": (1, 2, 3),
    "h2co_td5050": (1, 2, 6),
    "h2co_cis": (1, 3),
    "ch3": (1, 3),
    "benzene": (1, 3),
}


def fixture(name: str):
    """(fchk path, output path) of a fixture."""
    fchk, log = FIXTURES[name]
    return NTO_DIR / fchk, NTO_DIR / log


def multiwfn_reference(name: str, state: int, tmp_path) -> Path:
    """Multiwfn's exported NTO .fch, decompressed into ``tmp_path``."""
    source = NTO_DIR / f"{name}_state{state}_multiwfn.fch.gz"
    if not source.is_file():
        pytest.skip(f"no Multiwfn reference {source.name} (run tests/data/nto/make_reference.py)")
    target = Path(tmp_path) / f"{name}_state{state}_multiwfn.fch"
    with gzip.open(source, "rb") as src, open(target, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return target
