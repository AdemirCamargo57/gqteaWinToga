"""Optional numba kernels for grid work on Gaussian wavefunctions.

``cartesian_aos`` evaluates every normalised Cartesian basis component at
every point (the expensive half of rho(r) = sum_ab P_ab chi_a chi_b; the
contraction with P stays a BLAS product in gaussianWavefunction).
``becke_cell_weights`` gives Becke's fuzzy-cell weight of one centre at many
points, for multicentre integration (hirshfeldCharges).
``cartesian_ao_derivatives`` adds first and second AO derivatives, and
``density_derivatives`` fuses them with the contraction (rho, grad rho, its
Hessian and G per point, P summed only over unscreened AOs) for QTAIM: it
measured 8-65x faster than numpy and 1.4-3.6x faster than the AO kernel plus a
BLAS contraction, for 50 to 40 000 points (benzene 6-31G(d), water cc-pVQZ).

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


@numba.njit(inline="always", cache=True)
def _ipow(t, n):
    """t**n for small integers n >= 0; 0 for n < 0 (vanishing derivative terms)."""
    if n < 0:
        return 0.0
    value = 1.0
    for _ in range(n):
        value *= t
    return value


@numba.njit(inline="always", cache=True)
def _axis_factors(t, n, alpha):
    """d^m/dt^m [t^n exp(-a t^2)] / exp(-a t^2) for m = 0, 1, 2 (see gaussianWavefunction)."""
    f0 = _ipow(t, n)
    f1 = n * _ipow(t, n - 1) - 2.0 * alpha * _ipow(t, n + 1)
    f2 = n * (n - 1) * _ipow(t, n - 2) - 2.0 * alpha * (2 * n + 1) * f0 + 4.0 * alpha * alpha * _ipow(t, n + 2)
    return f0, f1, f2


@numba.njit(inline="always", cache=True)
def _point_ao_derivatives(px, py, pz, centers, prim_start, prim_count, min_exponent, exps, coefs,
                          comp_start, comp_count, powers, comp_norm, n_deriv, buf, active):
    """AO derivatives of one point into buf (n_deriv, n_cart); returns the number of active components.

    ``active`` receives the indices of components in unscreened shells; screened
    shells are left as zeros in ``buf``. Derivative order as DERIVATIVE_ORDER.
    """
    n_active = 0
    for s in range(centers.shape[0]):
        dx = px - centers[s, 0]
        dy = py - centers[s, 1]
        dz = pz - centers[s, 2]
        r2 = dx * dx + dy * dy + dz * dz
        c0 = comp_start[s]
        nc = comp_count[s]
        if min_exponent[s] * r2 > SHELL_SCREEN_EXPONENT:
            for c in range(nc):
                for k in range(n_deriv):
                    buf[k, c0 + c] = 0.0
            continue
        for c in range(nc):
            j = c0 + c
            for k in range(n_deriv):
                buf[k, j] = 0.0
            nx, ny, nz = powers[j, 0], powers[j, 1], powers[j, 2]
            for p in range(prim_start[s], prim_start[s] + prim_count[s]):
                a = exps[p]
                w = coefs[p] * np.exp(-a * r2)
                x0, x1, x2 = _axis_factors(dx, nx, a)
                y0, y1, y2 = _axis_factors(dy, ny, a)
                z0, z1, z2 = _axis_factors(dz, nz, a)
                buf[0, j] += w * x0 * y0 * z0
                buf[1, j] += w * x1 * y0 * z0
                buf[2, j] += w * x0 * y1 * z0
                buf[3, j] += w * x0 * y0 * z1
                if n_deriv > 4:
                    buf[4, j] += w * x2 * y0 * z0
                    buf[5, j] += w * x0 * y2 * z0
                    buf[6, j] += w * x0 * y0 * z2
                    buf[7, j] += w * x1 * y1 * z0
                    buf[8, j] += w * x1 * y0 * z1
                    buf[9, j] += w * x0 * y1 * z1
            for k in range(n_deriv):
                buf[k, j] *= comp_norm[j]
            active[n_active] = j
            n_active += 1
    return n_active


@numba.njit(parallel=True, cache=True)
def _ao_derivatives_kernel(points, centers, prim_start, prim_count, min_exponent, exps, coefs,
                           comp_start, comp_count, powers, comp_norm, out):
    n_deriv = out.shape[1]
    n_cart = out.shape[2]
    for i in numba.prange(points.shape[0]):
        active = np.empty(n_cart, dtype=np.int64)
        _point_ao_derivatives(points[i, 0], points[i, 1], points[i, 2], centers, prim_start, prim_count,
                              min_exponent, exps, coefs, comp_start, comp_count, powers, comp_norm,
                              n_deriv, out[i], active)


def cartesian_ao_derivatives(points: np.ndarray, packed, order: int) -> np.ndarray:
    """(n_points, 4 or 10, n_cartesian) AO values and derivatives (gaussianWavefunction.DERIVATIVE_ORDER)."""
    points = np.ascontiguousarray(points, dtype=np.float64)
    out = np.empty((points.shape[0], 4 if order == 1 else 10, packed.n_cart))
    _ao_derivatives_kernel(points, packed.centers, packed.prim_start, packed.prim_count, packed.min_exponent,
                           packed.exponents, packed.coefficients, packed.comp_start, packed.comp_count,
                           packed.powers, packed.comp_norm, out)
    return out


@numba.njit(parallel=True, cache=True)
def _density_derivatives_kernel(points, centers, prim_start, prim_count, min_exponent, exps, coefs,
                                comp_start, comp_count, powers, comp_norm, density, n_deriv,
                                rho, grad, hess, kin):
    """rho, grad rho, Hessian and G per point, contracting P only over unscreened AOs."""
    n_cart = density.shape[0]
    for i in numba.prange(points.shape[0]):
        buf = np.empty((n_deriv, n_cart))
        active = np.empty(n_cart, dtype=np.int64)
        m = _point_ao_derivatives(points[i, 0], points[i, 1], points[i, 2], centers, prim_start, prim_count,
                                  min_exponent, exps, coefs, comp_start, comp_count, powers, comp_norm,
                                  n_deriv, buf, active)
        r = gx = gy = gz = g = 0.0
        hxx = hyy = hzz = hxy = hxz = hyz = 0.0
        for bb in range(m):
            b = active[bb]
            vp = dpx = dpy = dpz = 0.0
            for aa in range(m):
                a = active[aa]
                p = density[a, b]
                vp += p * buf[0, a]
                dpx += p * buf[1, a]
                dpy += p * buf[2, a]
                dpz += p * buf[3, a]
            r += vp * buf[0, b]
            gx += vp * buf[1, b]
            gy += vp * buf[2, b]
            gz += vp * buf[3, b]
            g += dpx * buf[1, b] + dpy * buf[2, b] + dpz * buf[3, b]
            if n_deriv > 4:
                hxx += dpx * buf[1, b] + vp * buf[4, b]
                hyy += dpy * buf[2, b] + vp * buf[5, b]
                hzz += dpz * buf[3, b] + vp * buf[6, b]
                hxy += dpx * buf[2, b] + vp * buf[7, b]
                hxz += dpx * buf[3, b] + vp * buf[8, b]
                hyz += dpy * buf[3, b] + vp * buf[9, b]
        rho[i] = r
        grad[i, 0] = 2.0 * gx
        grad[i, 1] = 2.0 * gy
        grad[i, 2] = 2.0 * gz
        kin[i] = 0.5 * g
        if n_deriv > 4:
            hess[i, 0, 0] = 2.0 * hxx
            hess[i, 1, 1] = 2.0 * hyy
            hess[i, 2, 2] = 2.0 * hzz
            hess[i, 0, 1] = hess[i, 1, 0] = 2.0 * hxy
            hess[i, 0, 2] = hess[i, 2, 0] = 2.0 * hxz
            hess[i, 1, 2] = hess[i, 2, 1] = 2.0 * hyz


def density_derivatives(points: np.ndarray, packed, density: np.ndarray, order: int):
    """(rho, grad, Hessian or None, G) at every point; ``density`` is the Cartesian density matrix."""
    points = np.ascontiguousarray(points, dtype=np.float64)
    n = points.shape[0]
    n_deriv = 4 if order == 1 else 10
    rho, grad, kin = np.empty(n), np.empty((n, 3)), np.empty(n)
    hess = np.empty((n, 3, 3)) if order == 2 else np.empty((1, 3, 3))
    _density_derivatives_kernel(points, packed.centers, packed.prim_start, packed.prim_count, packed.min_exponent,
                                packed.exponents, packed.coefficients, packed.comp_start, packed.comp_count,
                                packed.powers, packed.comp_norm, np.ascontiguousarray(density, dtype=np.float64),
                                n_deriv, rho, grad, hess, kin)
    return rho, grad, (hess if order == 2 else None), kin


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
