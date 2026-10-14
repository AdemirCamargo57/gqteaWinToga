"""Shared paths and loaders for the ESP test fixtures (not a test module)."""
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent / "data"
FCHK_DIR = DATA / "fchk"
REFERENCE_DIR = DATA / "esp_reference"
BOHR_RADIUS_ANGSTROM = 0.529177210544

# Every fixture except benzene (kept for the integration/timing test only).
REFERENCE_FIXTURES = [
    "eteno_b3lyp_631g",     # user's file: RB3LYP/6-31G, SP shells
    "water_631gd_6d",       # Cartesian d
    "water_631gd_5d",       # pure d
    "water_ccpvtz",         # pure f
    "water_ccpvqz",         # pure g
    "h2co_6311g2df_10f",    # Cartesian f (fchk f order)
    "ch3_ub3lyp_631gd",     # open shell, UB3LYP
]


def fchk_path(name: str) -> Path:
    return FCHK_DIR / f"{name}.fchk"


def load_reference(name: str):
    """(points_bohr, rho_au, esp_au) from tests/data/esp_reference/<name>.txt."""
    data = np.loadtxt(REFERENCE_DIR / f"{name}.txt", comments="#")
    return data[:, :3] / BOHR_RADIUS_ANGSTROM, data[:, 3], data[:, 4]
