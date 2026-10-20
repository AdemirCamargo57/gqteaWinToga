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
