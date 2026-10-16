"""Hirshfeld (stockholder) atomic populations and charges from a Gaussian fchk density.

Toga-free. Lengths in bohr. Same scheme as Multiwfn's population analysis
(main function 7, option 1, built-in densities), which it reproduces to
~1e-4 e on the test set (H, Li, B-F, Na-Cl):

* pro-atoms: spherically averaged neutral free-atom densities of
  ``proatomDensities`` (UB3LYP/cc-pVQZ, H-Ar), interpolated linearly in
  log(rho); w_A(r) = rho_A^0(|r - R_A|) / sum_B rho_B^0(|r - R_B|);
* integration: Becke multicentre grid, per atom 75 Gauss-Chebyshev
  (2nd kind) radial points mapped as r = r_m (1+x)/(1-x) with r_m the
  covalent radius, times a 434-point Lebedev sphere (Multiwfn's defaults);
* populations are scaled to sum to the wavefunction's electron count before
  q_A = Z_A - p_A, as Multiwfn does for its "Final atomic charges".

With numba installed the Becke weights and the densities come from
densityNumba; without it the numpy paths give the same numbers.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np
from scipy.integrate import lebedev_rule

import gaussianWavefunction as gw
import proatomDensities as pad
from help import AtomicData

try:  # optional compiled kernel (pip install numba)
    import densityNumba
except ImportError:
    densityNumba = None

BOHR_RADIUS_ANGSTROM = 0.529177210544  # coordinateConverter imports Toga, so the constant is repeated
DEFAULT_RADIAL_POINTS = 75
DEFAULT_LEBEDEV_ORDER = 35  # 434 points
POINT_CHUNK = 20000

_R = np.asarray(pad.RADII_BOHR)
_LOG_RHO = {s: np.log(np.maximum(np.asarray(v), 1e-300)) for s, v in pad.DENSITIES.items()}


def supported_elements() -> List[str]:
    return list(pad.DENSITIES)


def unsupported_elements(symbols: Sequence[str]) -> List[str]:
    return sorted({s for s in symbols if s not in pad.DENSITIES})


def proatom_density(symbol: str, r: np.ndarray) -> np.ndarray:
    """Spherical free-atom density at distances r (bohr); 0 beyond the table."""
    log_rho = _LOG_RHO[symbol]
    return np.exp(np.interp(r, _R, log_rho, left=log_rho[0], right=-np.inf))


@dataclass
class MolecularGrid:
    points: np.ndarray   # (P, 3) bohr
    weights: np.ndarray  # (P,) quadrature x Becke cell weight


def _becke_weights_numpy(points, coords, centre):
    d = np.linalg.norm(points[:, None, :] - coords[None, :, :], axis=2)
    r = np.linalg.norm(coords[:, None] - coords[None, :], axis=2)
    cell = np.ones_like(d)
    for b in range(len(coords)):
        for c in range(len(coords)):
            if b != c:
                mu = (d[:, b] - d[:, c]) / r[b, c]
                for _ in range(3):
                    mu = 1.5 * mu - 0.5 * mu ** 3
                cell[:, b] *= 0.5 * (1.0 - mu)
    total = cell.sum(axis=1)
    return np.where(total > 0, cell[:, centre] / np.where(total > 0, total, 1.0), 0.0)


def molecular_grid(coordinates: np.ndarray, symbols: Sequence[str],
                   n_radial: int = DEFAULT_RADIAL_POINTS, lebedev_order: int = DEFAULT_LEBEDEV_ORDER,
                   use_numba: Optional[bool] = None) -> MolecularGrid:
    """Becke multicentre grid: integrates any smooth molecular function exactly enough for charges."""
    use_numba = gw._resolve_numba(use_numba) and densityNumba is not None
    coords = np.asarray(coordinates, dtype=float).reshape(-1, 3)
    sphere, w_sphere = lebedev_rule(lebedev_order)
    i = np.arange(1, n_radial + 1)
    theta = i * np.pi / (n_radial + 1)
    x = np.cos(theta)
    points, weights = [], []
    for a, symbol in enumerate(symbols):
        rm = AtomicData.covalent_radii.get(symbol, 1.5) / BOHR_RADIUS_ANGSTROM
        r = rm * (1 + x) / (1 - x)
        w_r = np.pi / (n_radial + 1) * np.sin(theta) * 2.0 * rm / (1 - x) ** 2 * r ** 2
        p = coords[a] + (r[:, None, None] * sphere.T[None, :, :]).reshape(-1, 3)
        w = (w_r[:, None] * w_sphere[None, :]).ravel()
        if len(coords) > 1:
            if use_numba:
                w = w * densityNumba.becke_cell_weights(p, coords, a)
            else:
                w = w * np.concatenate([_becke_weights_numpy(p[k:k + POINT_CHUNK], coords, a)
                                        for k in range(0, len(p), POINT_CHUNK)])
        points.append(p)
        weights.append(w)
    return MolecularGrid(np.vstack(points), np.concatenate(weights))


def hirshfeld_populations(wfn: gw.FchkWavefunction, grid: Optional[MolecularGrid] = None,
                          normalize: bool = True, use_numba: Optional[bool] = None) -> np.ndarray:
    """Electrons assigned to each atom. ValueError for elements without a free-atom density."""
    symbols = wfn.atom_symbols()
    missing = unsupported_elements(symbols)
    if missing:
        raise ValueError(f"No free-atom density for {', '.join(missing)}: Hirshfeld charges are "
                         f"available for {supported_elements()[0]}-{supported_elements()[-1]} only.")
    if not np.allclose(wfn.nuclear_charges, wfn.atomic_numbers):
        raise ValueError("The file uses effective core potentials: the density lacks the core "
                         "electrons, so Hirshfeld charges against all-electron free atoms would be wrong.")
    if grid is None:
        grid = molecular_grid(wfn.coordinates, symbols, use_numba=use_numba)
    populations = np.zeros(len(symbols))
    for k in range(0, len(grid.points), POINT_CHUNK):
        pts = grid.points[k:k + POINT_CHUNK]
        d = np.linalg.norm(pts[:, None, :] - wfn.coordinates[None, :, :], axis=2)
        pro = np.column_stack([proatom_density(s, d[:, a]) for a, s in enumerate(symbols)])
        total = pro.sum(axis=1)
        share = np.where(total[:, None] > 0, pro / np.where(total > 0, total, 1.0)[:, None], 0.0)
        rho = gw.electron_density(wfn, pts, use_numba=use_numba)
        populations += (grid.weights[k:k + POINT_CHUNK] * rho) @ share
    if normalize and populations.sum() > 0:
        populations *= wfn.n_electrons / populations.sum()
    return populations


def hirshfeld_charges(wfn: gw.FchkWavefunction, grid: Optional[MolecularGrid] = None,
                      normalize: bool = True, use_numba: Optional[bool] = None) -> np.ndarray:
    """q_A = Z_A - p_A (e). Atomic numbers, not ECP-reduced nuclear charges."""
    return wfn.atomic_numbers - hirshfeld_populations(wfn, grid, normalize, use_numba)
