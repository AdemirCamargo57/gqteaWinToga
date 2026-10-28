"""Paths and readers for the USPP library and the runatom reference outputs."""

import gzip
from pathlib import Path
from typing import Dict, List

SRC_DIR = Path(__file__).resolve().parents[1]
LIBRARY_DIR = SRC_DIR / "third_party" / "uspp-736" / "Work"
REFERENCE_DIR = Path(__file__).resolve().parent / "data" / "runatom"
REFERENCE_CASES = {
    "000-Ex-ca-sp-test": "000-Ex/000-Ex-ca-sp-test",
    "001-H-gpbe--bm": "001-H/001-H-gpbe--bm",
    "008-O-ca--vgrp": "008-O/008-O-ca--vgrp",
    "026-Fe-gpw-n-campos": "026-Fe/026-Fe-gpw-n-campos",
    "029-Cu-gpbe-n-bm": "029-Cu/029-Cu-gpbe-n-bm",
    "079-Au-ca--ts": "079-Au/079-Au-ca--ts",
    "082-Pb-ca-d-vgrp": "082-Pb/082-Pb-ca-d-vgrp",
    "056-Ba-gpw-nsp-campos": "056-Ba/056-Ba-gpw-nsp-campos",
}


def library_case_dir(name: str) -> Path:
    return LIBRARY_DIR / REFERENCE_CASES[name]


def case_info(name: str) -> Dict[str, str]:
    rows = (REFERENCE_DIR / name / "CASE.txt").read_text().splitlines()
    return {k.strip(): v.strip() for k, v in (row.split(":", 1) for row in rows if ":" in row)}


def reference_names(name: str) -> List[str]:
    return sorted(p.name[:-3] for p in (REFERENCE_DIR / name).glob("*.gz"))


def reference_bytes(name: str, filename: str) -> bytes:
    return gzip.decompress((REFERENCE_DIR / name / f"{filename}.gz").read_bytes())


def reference_text(name: str, filename: str) -> str:
    return reference_bytes(name, filename).decode("latin-1")
