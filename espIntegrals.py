"""Exact molecular electrostatic potential of a Gaussian-basis wavefunction.

    V(C) = sum_A Z_A / |C - R_A|  -  sum_ab P_ab  int chi_a chi_b / |r - C| dr

Z_A are the fchk ``Nuclear charges`` (effective charges with an ECP, as in
Multiwfn). The electronic term uses McMurchie-Davidson: every primitive
pair's share of the density is reduced **once** to Hermite coefficients
D_tuv (independent of C), so the potential at C is

    sum_pairs  sum_tuv  D_tuv R_tuv(p, P - C),

with R_tuv from the Boys function. Pure numpy/scipy; point chunks run in a
thread pool, since numpy and scipy ufuncs release the GIL. Verified against
Gaussian 09 cubegen to <= 4e-7 a.u. on every fixture in tests/data/fchk.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from math import pi
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.special import gamma, gammainc

from gaussianWavefunction import (FchkWavefunction, cartesian_components,
                                  component_angular_norm, primitive_norm)

try:  # optional compiled kernel (pip install numba)
    import espNumba
except ImportError:
    espNumba = None

# A primitive pair is dropped when its largest possible contribution,
# |prefactor| x max|P block|, is below this (a.u.).
PAIR_SCREEN = 1e-14
# Elements of one (pairs x points) work array; bounds memory per thread.
BLOCK_ELEMENTS = 2_000_000


def boys(n_max: int, t: np.ndarray) -> np.ndarray:
    """Boys functions F_0..F_n_max at t, stacked on a new leading axis.

    F_n(t) = int_0^1 u^(2n) exp(-t u^2) du. F_n_max comes from the
    regularised incomplete gamma function, the rest by the stable downward
    recursion F_(m-1) = (2 t F_m + exp(-t)) / (2m - 1).
    """
    t = np.asarray(t, dtype=float)
    out = np.empty((n_max + 1,) + t.shape)
    small = t < 1e-10
    a = n_max + 0.5
    ts = np.where(small, 1.0, t)
    top = gamma(a) * gammainc(a, ts) / (2.0 * ts ** a)
    out[n_max] = np.where(small, 1.0 / (2 * n_max + 1) - t / (2 * n_max + 3), top)
    emt = np.exp(-t)
    for m in range(n_max, 0, -1):
        out[m - 1] = (2.0 * t * out[m] + emt) / (2 * m - 1)
    return out


def hermite_coefficients(la: int, lb: int, xpa: np.ndarray, xpb: np.ndarray,
                         inv2p: np.ndarray) -> np.ndarray:
    """E[pair, i, j, t] (one Cartesian direction, without the K prefactor).

    E^{i+1,j}_t = 1/(2p) E^{ij}_{t-1} + X_PA E^{ij}_t + (t+1) E^{ij}_{t+1}
    E^{i,j+1}_t = 1/(2p) E^{ij}_{t-1} + X_PB E^{ij}_t + (t+1) E^{ij}_{t+1}
    """
    npair = xpa.shape[0]
    e = np.zeros((npair, la + 1, lb + 1, la + lb + 2))
    e[:, 0, 0, 0] = 1.0
    for i in range(la):
        for t in range(i + 2):
            prev = e[:, i, 0, t - 1] if t > 0 else 0.0
            e[:, i + 1, 0, t] = inv2p * prev + xpa * e[:, i, 0, t] + (t + 1) * e[:, i, 0, t + 1]
    for j in range(lb):
        for i in range(la + 1):
            for t in range(i + j + 2):
                prev = e[:, i, j, t - 1] if t > 0 else 0.0
                e[:, i, j + 1, t] = inv2p * prev + xpb * e[:, i, j, t] + (t + 1) * e[:, i, j, t + 1]
    return e[..., : la + lb + 1]


def hermite_indices(l: int) -> List[Tuple[int, int, int]]:
    """Every (t, u, v) with t + u + v <= l."""
    return [(t, u, v) for t in range(l + 1) for u in range(l + 1 - t) for v in range(l + 1 - t - u)]


@dataclass
class HermiteBatch:
    """All primitive pairs whose angular momenta sum to ``l``."""
    l: int
    p: np.ndarray        # (M,) total exponent
    centers: np.ndarray  # (M, 3) Gaussian product centre P
    coeffs: np.ndarray   # (M, n_tuv): 2 pi/p x K x contraction x density, summed over components


def _primitive_pairs(sa, sb):
    al = np.repeat(sa.exponents, sb.exponents.size)
    be = np.tile(sb.exponents, sa.exponents.size)
    ca = np.repeat(sa.coefficients * primitive_norm(sa.exponents, sa.l), sb.exponents.size)
    cb = np.tile(sb.coefficients * primitive_norm(sb.exponents, sb.l), sa.exponents.size)
    p = al + be
    ab = sa.center - sb.center
    k = np.exp(-al * be / p * (ab @ ab))
    centers = (al[:, None] * sa.center + be[:, None] * sb.center) / p[:, None]
    return al, be, p, ca * cb * k, centers


def build_hermite_batches(wfn: FchkWavefunction, density: Optional[np.ndarray] = None) -> List[HermiteBatch]:
    """Contract a density into point-independent Hermite coefficients.

    ``density`` is a symmetric matrix in the fchk basis (default: the SCF
    density of ``wfn``), e.g. the electron density of an excitation.
    """
    cart_density = wfn.cart_density if density is None else wfn.transform.T @ density @ wfn.transform
    shells = wfn.shells
    offsets = np.cumsum([0] + [s.n_cart for s in shells])
    acc: Dict[int, List[Tuple[np.ndarray, np.ndarray, np.ndarray]]] = {}
    for ia, sa in enumerate(shells):
        for ib in range(ia + 1):
            sb = shells[ib]
            block = cart_density[offsets[ia]:offsets[ia + 1], offsets[ib]:offsets[ib + 1]]
            block_max = np.abs(block).max()
            if block_max == 0.0:
                continue
            _, _, p, ckk, centers = _primitive_pairs(sa, sb)
            factor = 2.0 if ia != ib else 1.0  # the (b, a) block is the transpose
            pre = factor * 2.0 * pi / p * ckk
            keep = np.abs(pre) * block_max > PAIR_SCREEN
            if not keep.any():
                continue
            p, pre, centers = p[keep], pre[keep], centers[keep]
            inv2p = 0.5 / p
            ex, ey, ez = (hermite_coefficients(sa.l, sb.l, centers[:, d] - sa.center[d],
                                               centers[:, d] - sb.center[d], inv2p) for d in range(3))
            l = sa.l + sb.l
            tuv = np.array(hermite_indices(l))
            coeffs = np.zeros((p.size, len(tuv)))
            for a, pa in enumerate(cartesian_components(sa.l)):
                for b, pb in enumerate(cartesian_components(sb.l)):
                    dab = block[a, b] * component_angular_norm(pa) * component_angular_norm(pb)
                    if dab == 0.0:
                        continue
                    coeffs += dab * (ex[:, pa[0], pb[0], tuv[:, 0]] * ey[:, pa[1], pb[1], tuv[:, 1]]
                                     * ez[:, pa[2], pb[2], tuv[:, 2]])
            acc.setdefault(l, []).append((p, centers, coeffs * pre[:, None]))
    return [HermiteBatch(l, np.concatenate([x[0] for x in v]), np.concatenate([x[1] for x in v]),
                         np.concatenate([x[2] for x in v])) for l, v in sorted(acc.items())]


def hermite_integrals(l: int, p: np.ndarray, pc: np.ndarray) -> Dict[Tuple[int, int, int], np.ndarray]:
    """R^0_tuv(p, P - C) for t + u + v <= l; ``pc`` is (3, M, N).

    R^n_000 = (-2p)^n F_n(p |PC|^2), then e.g.
    R^n_{t+1,u,v} = t R^{n+1}_{t-1,u,v} + X_PC R^{n+1}_{t,u,v},
    keeping only two levels of n alive at a time.
    """
    x, y, z = pc
    f = boys(l, p[:, None] * (x * x + y * y + z * z))
    m2p = -2.0 * p[:, None]
    level = {(0, 0, 0): (m2p ** l) * f[l]}
    for n in range(l - 1, -1, -1):
        new = {}
        for (t, u, v) in hermite_indices(l - n):
            if t > 0:
                val = x * level[(t - 1, u, v)]
                if t > 1:
                    val = val + (t - 1) * level[(t - 2, u, v)]
            elif u > 0:
                val = y * level[(t, u - 1, v)]
                if u > 1:
                    val = val + (u - 1) * level[(t, u - 2, v)]
            elif v > 0:
                val = z * level[(t, u, v - 1)]
                if v > 1:
                    val = val + (v - 1) * level[(t, u, v - 2)]
            else:
                val = (m2p ** n) * f[n]
            new[(t, u, v)] = val
        level = new
    return level


def electronic_potential(batches: List[HermiteBatch], points: np.ndarray) -> np.ndarray:
    """int rho(r) / |r - C| dr at every point C (positive; a.u.)."""
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    out = np.zeros(points.shape[0])
    for batch in batches:
        tuv = hermite_indices(batch.l)
        npt = max(1, BLOCK_ELEMENTS // max(batch.p.size, 1))
        for a in range(0, points.shape[0], npt):
            c = points[a:a + npt]
            pc = batch.centers.T[:, :, None] - c.T[:, None, :]
            r = hermite_integrals(batch.l, batch.p, pc)
            stacked = np.stack([r[key] for key in tuv], axis=1)  # (M, n_tuv, N)
            out[a:a + npt] += np.einsum("mk,mkn->n", batch.coeffs, stacked)
    return out


def nuclear_potential(wfn: FchkWavefunction, points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    d = np.linalg.norm(points[:, None, :] - wfn.coordinates[None, :, :], axis=2)
    return (wfn.nuclear_charges[None, :] / d).sum(axis=1)


def default_workers() -> int:
    return max(1, min(12, (os.cpu_count() or 1)))


def numba_available() -> bool:
    return espNumba is not None


def electrostatic_potential(wfn: FchkWavefunction, points: np.ndarray,
                            batches: Optional[List[HermiteBatch]] = None,
                            workers: Optional[int] = None,
                            progress: Optional[Callable[[float], None]] = None,
                            use_numba: Optional[bool] = None) -> np.ndarray:
    """Total ESP (a.u., Hartree/e) at ``points`` (bohr).

    ``use_numba`` None means "if installed" (the compiled kernel is ~5x
    faster on benzene and agrees with numpy to 1e-10). ``progress(fraction)``
    is called as point chunks finish; a callback that raises is ignored.
    """
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    if batches is None:
        batches = build_hermite_batches(wfn)
    if use_numba is None:
        use_numba = numba_available()
    if use_numba and not numba_available():
        raise ValueError("numba is not installed.")
    electronic = espNumba.electronic_potential if use_numba else electronic_potential
    # The numba kernel is parallel itself, so its chunks run one at a time.
    workers = 1 if use_numba else (workers or default_workers())
    n_chunks = max(1, min(len(points), 20 if use_numba else workers * 4))
    chunks = np.array_split(np.arange(len(points)), n_chunks)
    out = np.empty(len(points))
    done = [0]

    def work(idx):
        out[idx] = nuclear_potential(wfn, points[idx]) - electronic(batches, points[idx])
        done[0] += len(idx)
        if progress is not None:
            try:
                progress(done[0] / len(points))
            except Exception:
                pass

    if workers == 1:
        for idx in chunks:
            work(idx)
    else:
        with ThreadPoolExecutor(workers) as pool:
            list(pool.map(work, chunks))
    return out


def overlap_matrix(wfn: FchkWavefunction) -> np.ndarray:
    """Overlap of the normalised Cartesian components (validation: tr(P S) = N)."""
    shells = wfn.shells
    offsets = np.cumsum([0] + [s.n_cart for s in shells])
    s_mat = np.zeros((offsets[-1], offsets[-1]))
    for ia, sa in enumerate(shells):
        for ib, sb in enumerate(shells):
            _, _, p, ckk, centers = _primitive_pairs(sa, sb)
            pre = ckk * (pi / p) ** 1.5
            e = [hermite_coefficients(sa.l, sb.l, centers[:, d] - sa.center[d],
                                      centers[:, d] - sb.center[d], 0.5 / p) for d in range(3)]
            for a, pa in enumerate(cartesian_components(sa.l)):
                for b, pb in enumerate(cartesian_components(sb.l)):
                    s_mat[offsets[ia] + a, offsets[ib] + b] = (
                        component_angular_norm(pa) * component_angular_norm(pb)
                        * np.sum(pre * e[0][:, pa[0], pb[0], 0] * e[1][:, pa[1], pb[1], 0]
                                 * e[2][:, pa[2], pb[2], 0]))
    return s_mat


MOMENT_KEYS = ("s", "x", "y", "z", "xx", "yy", "zz")


def moment_matrices(wfn: FchkWavefunction) -> Dict[str, np.ndarray]:
    """Overlap, dipole (x, y, z) and diagonal second-moment (xx, yy, zz) integrals.

    Same layout as ``overlap_matrix`` (normalised Cartesian components), origin
    at (0, 0, 0) bohr. Along one axis, with x = (x - B_x) + B_x on the right-hand
    Gaussian, <a|x|b> = S(i, j+1) + B_x S(i, j) and
    <a|x^2|b> = S(i, j+2) + 2 B_x S(i, j+1) + B_x^2 S(i, j), where S(i, j) is the
    one-dimensional overlap E^{ij}_0 sqrt(pi/p).
    """
    shells = wfn.shells
    offsets = np.cumsum([0] + [s.n_cart for s in shells])
    out = {key: np.zeros((offsets[-1], offsets[-1])) for key in MOMENT_KEYS}
    for ia, sa in enumerate(shells):
        for ib, sb in enumerate(shells):
            _, _, p, ckk, centers = _primitive_pairs(sa, sb)
            root = np.sqrt(pi / p)
            e = [hermite_coefficients(sa.l, sb.l + 2, centers[:, d] - sa.center[d],
                                      centers[:, d] - sb.center[d], 0.5 / p)[..., 0] * root[:, None, None]
                 for d in range(3)]  # (pairs, la+1, lb+3): 1-D overlaps
            for a, pa in enumerate(cartesian_components(sa.l)):
                for b, pb in enumerate(cartesian_components(sb.l)):
                    norm = component_angular_norm(pa) * component_angular_norm(pb)
                    s1 = [e[d][:, pa[d], pb[d]] for d in range(3)]
                    first = [e[d][:, pa[d], pb[d] + 1] + sb.center[d] * s1[d] for d in range(3)]
                    second = [e[d][:, pa[d], pb[d] + 2] + 2.0 * sb.center[d] * e[d][:, pa[d], pb[d] + 1]
                              + sb.center[d] ** 2 * s1[d] for d in range(3)]
                    row, col = offsets[ia] + a, offsets[ib] + b
                    out["s"][row, col] = norm * np.sum(ckk * s1[0] * s1[1] * s1[2])
                    for d, (k1, k2) in enumerate((("x", "xx"), ("y", "yy"), ("z", "zz"))):
                        others = [s1[o] for o in range(3) if o != d]
                        out[k1][row, col] = norm * np.sum(ckk * first[d] * others[0] * others[1])
                        out[k2][row, col] = norm * np.sum(ckk * second[d] * others[0] * others[1])
    return out


# 1 Hartree/e in kcal/(mol e): CODATA 2018 E_h = 4.3597447222071e-18 J times
# N_A = 6.02214076e23 / mol, over 4184 J/kcal. Multiwfn uses 627.51.
HARTREE_TO_KCAL_MOL = 627.5094740631
