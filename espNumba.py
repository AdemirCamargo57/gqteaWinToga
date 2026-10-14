"""Optional numba kernel for the electronic ESP (imported by espIntegrals).

Same McMurchie-Davidson sum as espIntegrals.electronic_potential, but one
compiled loop over (point, primitive pair) instead of numpy arrays of shape
(pairs, points, n_tuv). Importing this module raises ImportError when numba
is not installed; espIntegrals then keeps the numpy path.
"""

from __future__ import annotations

from math import exp, pi, sqrt
from typing import Dict, Tuple

import numba
import numpy as np

# Boys: power series below this t, asymptotic form above. At 30 the
# asymptotic F_8 is still 1e-6 off; at 60 it is exact to double precision.
BOYS_SERIES_LIMIT = 60.0


@numba.njit(cache=True)
def _boys_into(n_max, t, out):
    if t < BOYS_SERIES_LIMIT:
        # F_n(t) = exp(-t) sum_k (2t)^k / ((2n+1)(2n+3)...(2n+2k+1))
        term = 1.0 / (2 * n_max + 1)
        total = term
        k = 1
        while term > 1e-17 * total:
            term *= 2.0 * t / (2 * n_max + 2 * k + 1)
            total += term
            k += 1
        out[n_max] = exp(-t) * total
    else:
        df = 1.0  # (2n-1)!!
        for m in range(1, n_max + 1):
            df *= 2 * m - 1
        out[n_max] = df / 2.0 ** (n_max + 1) * sqrt(pi / t ** (2 * n_max + 1))
    emt = exp(-t)
    for m in range(n_max, 0, -1):
        out[m - 1] = (2.0 * t * out[m] + emt) / (2 * m - 1)


@numba.njit(parallel=True, cache=True)
def _potential_kernel(l, p, centers, coeffs, points, kind, src1, src2, coef, out):
    """out[i] += sum_m sum_k coeffs[m, k] R^0_k(p_m, P_m - C_i).

    The recursion is a precomputed instruction list over the full tuv
    index of order l, applied level by level n = l-1 .. 0:
    kind 0 = base (-2p)^n F_n, 1/2/3 = x/y/z step with sources src1 (-1)
    and src2 (-2, or -1 for none) scaled by coef.
    """
    n_tuv = coeffs.shape[1]
    for i in numba.prange(points.shape[0]):
        prev = np.zeros(n_tuv)
        cur = np.zeros(n_tuv)
        f = np.zeros(l + 1)
        acc = 0.0
        for m in range(p.shape[0]):
            x = centers[m, 0] - points[i, 0]
            y = centers[m, 1] - points[i, 1]
            z = centers[m, 2] - points[i, 2]
            _boys_into(l, p[m] * (x * x + y * y + z * z), f)
            m2p = -2.0 * p[m]
            prev[0] = m2p ** l * f[l]
            for n in range(l - 1, -1, -1):
                count = (l - n + 1) * (l - n + 2) * (l - n + 3) // 6
                for j in range(count):
                    s = kind[j]
                    if s == 0:
                        val = m2p ** n * f[n]
                    else:
                        d = x if s == 1 else (y if s == 2 else z)
                        val = d * prev[src1[j]]
                        if src2[j] >= 0:
                            val += coef[j] * prev[src2[j]]
                    cur[j] = val
                prev, cur = cur, prev
            for k in range(n_tuv):
                acc += coeffs[m, k] * prev[k]
        out[i] += acc


_PLANS: Dict[int, Tuple[np.ndarray, ...]] = {}


def recursion_plan(l: int):
    """Instruction arrays for _potential_kernel, cached per l.

    Entries are listed in hermite_indices order; that order is a prefix-
    closed list (every order-(l-n) index set is a prefix of the order-l
    one), so position j in a level is the same tuv in every level.
    """
    if l in _PLANS:
        return _PLANS[l]
    from espIntegrals import hermite_indices
    full = _graded_indices(l)
    pos = {key: j for j, key in enumerate(full)}
    kind = np.zeros(len(full), np.int64)
    src1 = np.zeros(len(full), np.int64)
    src2 = -np.ones(len(full), np.int64)
    coef = np.zeros(len(full))
    for j, (t, u, v) in enumerate(full):
        if t > 0:
            kind[j], src1[j] = 1, pos[(t - 1, u, v)]
            if t > 1:
                src2[j], coef[j] = pos[(t - 2, u, v)], t - 1
        elif u > 0:
            kind[j], src1[j] = 2, pos[(t, u - 1, v)]
            if u > 1:
                src2[j], coef[j] = pos[(t, u - 2, v)], u - 1
        elif v > 0:
            kind[j], src1[j] = 3, pos[(t, u, v - 1)]
            if v > 1:
                src2[j], coef[j] = pos[(t, u, v - 2)], v - 1
    order = np.array([pos[key] for key in hermite_indices(l)])
    _PLANS[l] = (kind, src1, src2, coef, order)
    return _PLANS[l]


def _graded_indices(l: int):
    """All (t, u, v) with t+u+v <= l, sorted by total order: prefix-closed."""
    out = []
    for total in range(l + 1):
        for t in range(total, -1, -1):
            for u in range(total - t, -1, -1):
                out.append((t, u, total - t - u))
    return out


def electronic_potential(batches, points: np.ndarray) -> np.ndarray:
    points = np.ascontiguousarray(points, dtype=float).reshape(-1, 3)
    out = np.zeros(points.shape[0])
    for b in batches:
        kind, src1, src2, coef, order = recursion_plan(b.l)
        # coeffs are stored in hermite_indices order; the kernel uses graded order
        coeffs = np.zeros((b.coeffs.shape[0], len(kind)))
        coeffs[:, order] = b.coeffs
        _potential_kernel(b.l, b.p, np.ascontiguousarray(b.centers), coeffs, points,
                          kind, src1, src2, coef, out)
    return out
