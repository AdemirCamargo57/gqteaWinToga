"""Shared QTAIM test data: a synthetic H2 wavefunction and fixture paths."""
from math import pi
from pathlib import Path

import numpy as np

import gaussianWavefunction as gw

FCHK = Path(__file__).parent / "data" / "fchk"
QTAIM = Path(__file__).parent / "data" / "qtaim"

# STO-3G hydrogen (exponents and coefficients of normalised primitives, as in an fchk)
STO3G_H_EXPONENTS = np.array([3.42525091, 0.62391373, 0.16885540])
STO3G_H_COEFFICIENTS = np.array([0.15432897, 0.53532814, 0.44463454])


def _contracted_overlap(distance: float) -> float:
    a, c = STO3G_H_EXPONENTS, STO3G_H_COEFFICIENTS
    n = (2 * a / pi) ** 0.75
    p = a[:, None] + a[None, :]
    s = (pi / p) ** 1.5 * np.exp(-(a[:, None] * a[None, :] / p) * distance ** 2)
    return float(np.sum(c[:, None] * c[None, :] * n[:, None] * n[None, :] * s))


def make_h2(distance: float = 1.4) -> gw.FchkWavefunction:
    """H2 along z (atoms at z = 0 and z = distance bohr), doubly occupied bonding MO."""
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, distance]])
    shells = [gw.Shell(0, False, coords[i], i, STO3G_H_EXPONENTS.copy(), STO3G_H_COEFFICIENTS.copy())
              for i in range(2)]
    s = _contracted_overlap(distance)
    c = np.full(2, 1.0 / np.sqrt(2.0 * (1.0 + s)))
    density = 2.0 * np.outer(c, c)
    wfn = gw.FchkWavefunction(
        path="H2.fchk", title="H2", atomic_numbers=np.array([1, 1]), nuclear_charges=np.array([1.0, 1.0]),
        coordinates=coords, shells=shells, density=density, n_electrons=2, charge=0, multiplicity=1,
        method="RHF", basis_set="STO-3G", n_alpha=1, n_beta=1)
    wfn.transform = gw.basis_transform(shells)
    wfn.cart_density = wfn.transform.T @ density @ wfn.transform
    return wfn


def fchk(name: str) -> Path:
    """A fixture by file name: tests/data/qtaim first, then tests/data/fchk."""
    path = QTAIM / name
    return path if path.is_file() else FCHK / name


REFERENCE_FIXTURES = ["water_631gd_6d.fchk", "eteno_b3lyp_631g.fchk", "benzene_631gd.fchk",
                      "ch3_ub3lyp_631gd.fchk", "water_dimer_631gd.fchk", "bicyclopentane_631gd.fchk"]


class ReferenceCP:
    def __init__(self, index, kind, position):
        self.index, self.kind, self.position = index, kind, np.array(position)
        self.props = None   # rho lap G K V H l1 l2 l3 ellipticity eta
        self.atoms = None   # frozenset of 0-based atoms (BCPs)


def load_reference(name: str):
    """(cps, paths) from tests/data/qtaim/<stem>_multiwfn.txt (see make_reference.py)."""
    cps, paths = [], []
    for line in (QTAIM / (Path(name).stem + "_multiwfn.txt")).read_text().splitlines():
        parts = line.split()
        if not parts or parts[0].startswith("#"):
            continue
        if parts[0] == "cp":
            cps.append(ReferenceCP(int(parts[1]), parts[2], [float(x) for x in parts[3:6]]))
        elif parts[0] == "prop":
            cps[-1].props = np.array([float(x) for x in parts[1:]])
        elif parts[0] == "atoms":
            cps[-1].atoms = frozenset(int(x) - 1 for x in parts[1:3])
        elif parts[0] == "path":
            paths.append((int(parts[1]), int(parts[2]), float(parts[3])))
    return cps, paths
