"""Optional numba kernels for grid work on Gaussian wavefunctions.

``cartesian_aos`` evaluates every normalised Cartesian basis component at
every point (the expensive half of rho(r) = sum_ab P_ab chi_a chi_b; the
contraction with P stays a BLAS product in gaussianWavefunction).
``becke_cell_weights`` gives Becke's fuzzy-cell weight of one centre at many
points, for multicentre integration (hirshfeldCharges).

Importing this module raises ImportError when numba is not installed; the
callers then keep their numpy paths, which give the same numbers.
"""

from __future__ import annotations

from math import exp, sqrt

import numba
import numpy as np

# A shell is skipped at a point when even its most diffuse primitive is below
# exp(-SHELL_SCREEN_EXPONENT) there: < 1e-22 of its coefficient, far below
# double-precision relevance for densities of order 1e-6 and larger.
SHELL_SCREEN_EXPONENT = 50.0


@numba.njit(parallel=True, cache=True)
def _aos_kernel(points, centers, prim_start, prim_count, min_exponent, exps, coefs,
                comp_start, comp_count, powers, comp_norm, out):
    n_shell = centers.shape[0]
    for i in numba.prange(points.shape[0]):
        px, py, pz = points[i, 0], points[i, 1], points[i, 2]
        for s in range(n_shell):
            dx = px - centers[s, 0]
            dy = py - centers[s, 1]
            dz = pz - centers[s, 2]
            r2 = dx * dx + dy * dy + dz * dz
            c0 = comp_start[s]
            nc = comp_count[s]
            if min_exponent[s] * r2 > SHELL_SCREEN_EXPONENT:
                for c in range(nc):
                    out[i, c0 + c] = 0.0
                continue
            radial = 0.0
            for k in range(prim_start[s], prim_start[s] + prim_count[s]):
                radial += coefs[k] * exp(-exps[k] * r2)
            for c in range(nc):
                j = c0 + c
                value = radial * comp_norm[j]
                for _ in range(powers[j, 0]):
                    value *= dx
                for _ in range(powers[j, 1]):
                    value *= dy
                for _ in range(powers[j, 2]):
                    value *= dz
                out[i, j] = value


def cartesian_aos(points: np.ndarray, packed) -> np.ndarray:
    """(n_points, n_cartesian) AO values; ``packed`` from gaussianWavefunction.pack_shells."""
    points = np.ascontiguousarray(points, dtype=np.float64)
    out = np.empty((points.shape[0], packed.n_cart))
    _aos_kernel(points, packed.centers, packed.prim_start, packed.prim_count, packed.min_exponent,
                packed.exponents, packed.coefficients, packed.comp_start, packed.comp_count,
                packed.powers, packed.comp_norm, out)
    return out


@numba.njit(parallel=True, cache=True)
def _becke_kernel(points, coords, inv_r, centre, out):
    n_atom = coords.shape[0]
    for i in numba.prange(points.shape[0]):
        d = np.empty(n_atom)
        for a in range(n_atom):
            dx = points[i, 0] - coords[a, 0]
            dy = points[i, 1] - coords[a, 1]
            dz = points[i, 2] - coords[a, 2]
            d[a] = sqrt(dx * dx + dy * dy + dz * dz)
        total = 0.0
        mine = 0.0
        for b in range(n_atom):
            cell = 1.0
            for c in range(n_atom):
                if b == c:
                    continue
                mu = (d[b] - d[c]) * inv_r[b, c]
                for _ in range(3):
                    mu = 1.5 * mu - 0.5 * mu * mu * mu
                cell *= 0.5 * (1.0 - mu)
                if cell == 0.0:
                    break
            total += cell
            if b == centre:
                mine = cell
        out[i] = mine / total if total > 0.0 else 0.0


def becke_cell_weights(points: np.ndarray, coords: np.ndarray, centre: int) -> np.ndarray:
    """Becke (1988) fuzzy-cell weight of atom ``centre`` at each point (3 smoothing steps)."""
    coords = np.ascontiguousarray(coords, dtype=np.float64)
    r = np.linalg.norm(coords[:, None] - coords[None, :], axis=2)
    inv_r = np.where(r > 0, 1.0 / np.where(r > 0, r, 1.0), 0.0)
    out = np.empty(len(points))
    _becke_kernel(np.ascontiguousarray(points, dtype=np.float64), coords, inv_r, int(centre), out)
    return out
