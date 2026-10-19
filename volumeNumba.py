"""Optional numba kernel for isosurfaces of volumetric grids (the Cube tab).

``tetra_cut_edges`` is the compiled equivalent of
``molecularSurface.tetra_cut_edges``: the cell scan of marching tetrahedra.
It produces the same triangles -- the same cut edges, in the same vertex
order, with the same outward directions -- but scans the cells in parallel
without the (cells x 8) boolean stack the numpy path builds, so it is both
faster and much lighter on memory for large cubes. Triangles come out in
cell order rather than tetrahedron order; ``mesh_from_cut_edges`` sorts the
edge ids, so the vertex arrays of both paths are identical.

Importing this module raises ImportError when numba is not installed; the
caller (volumeData) then keeps the numpy path.
"""

from __future__ import annotations

import numba
import numpy as np

from molecularSurface import _CORNER_OFFSETS, _EDGE_LUT, _TETRAHEDRA

_OFFSETS = np.ascontiguousarray(_CORNER_OFFSETS, dtype=np.int64)
_TETS = np.ascontiguousarray(_TETRAHEDRA, dtype=np.int64)
_LUT = np.ascontiguousarray(_EDGE_LUT, dtype=np.int64)


@numba.njit(cache=True)
def _cell_corners(values, x, y, z, isovalue, offsets, inside):
    """Fill ``inside`` for the 8 corners; True when the cell is cut."""
    n_in = 0
    for c in range(8):
        inside[c] = values[x + offsets[c, 0], y + offsets[c, 1], z + offsets[c, 2]] >= isovalue
        if inside[c]:
            n_in += 1
    return 0 < n_in < 8


@numba.njit(cache=True)
def _emit(x, y, z, ny, nz, ca, cb, offsets, lut, edges, dirs, row, col):
    """Edge id and inside->outside step of the edge from corner ca (in) to cb (out)."""
    ax, ay, az = x + offsets[ca, 0], y + offsets[ca, 1], z + offsets[ca, 2]
    bx, by, bz = x + offsets[cb, 0], y + offsets[cb, 1], z + offsets[cb, 2]
    lx, ly, lz = min(ax, bx), min(ay, by), min(az, bz)
    k = lut[abs(bx - ax) + 2 * abs(by - ay) + 4 * abs(bz - az)]
    edges[row, col] = ((lx * ny + ly) * nz + lz) * 7 + k
    dirs[row, 0] += bx - ax
    dirs[row, 1] += by - ay
    dirs[row, 2] += bz - az


@numba.njit(parallel=True, cache=True)
def _count(values, isovalue, offsets, tets):
    nx, ny, nz = values.shape
    counts = np.zeros(max(nx - 1, 0), dtype=np.int64)
    for x in numba.prange(nx - 1):
        inside = np.empty(8, dtype=np.bool_)
        n = 0
        for y in range(ny - 1):
            for z in range(nz - 1):
                if not _cell_corners(values, x, y, z, isovalue, offsets, inside):
                    continue
                for t in range(6):
                    n_in = 0
                    for k in range(4):
                        if inside[tets[t, k]]:
                            n_in += 1
                    if n_in == 1 or n_in == 3:
                        n += 1
                    elif n_in == 2:
                        n += 2
        counts[x] = n
    return counts


@numba.njit(parallel=True, cache=True)
def _fill(values, isovalue, offsets, tets, lut, starts, edges, dirs):
    nx, ny, nz = values.shape
    for x in numba.prange(nx - 1):
        inside = np.empty(8, dtype=np.bool_)
        ins = np.empty(4, dtype=np.int64)
        outs = np.empty(4, dtype=np.int64)
        row = starts[x]
        for y in range(ny - 1):
            for z in range(nz - 1):
                if not _cell_corners(values, x, y, z, isovalue, offsets, inside):
                    continue
                for t in range(6):
                    ni = 0
                    no = 0
                    for k in range(4):
                        c = tets[t, k]
                        if inside[c]:
                            ins[ni] = c
                            ni += 1
                        else:
                            outs[no] = c
                            no += 1
                    if ni == 1:  # (i, o0), (i, o1), (i, o2)
                        for j in range(3):
                            _emit(x, y, z, ny, nz, ins[0], outs[j], offsets, lut, edges, dirs, row, j)
                        row += 1
                    elif ni == 3:  # (i0, o), (i1, o), (i2, o)
                        for j in range(3):
                            _emit(x, y, z, ny, nz, ins[j], outs[0], offsets, lut, edges, dirs, row, j)
                        row += 1
                    elif ni == 2:  # quad (i0o0, i0o1, i1o1, i1o0) -> (q0, q1, q2), (q0, q2, q3)
                        _emit(x, y, z, ny, nz, ins[0], outs[0], offsets, lut, edges, dirs, row, 0)
                        _emit(x, y, z, ny, nz, ins[0], outs[1], offsets, lut, edges, dirs, row, 1)
                        _emit(x, y, z, ny, nz, ins[1], outs[1], offsets, lut, edges, dirs, row, 2)
                        row += 1
                        _emit(x, y, z, ny, nz, ins[0], outs[0], offsets, lut, edges, dirs, row, 0)
                        _emit(x, y, z, ny, nz, ins[1], outs[1], offsets, lut, edges, dirs, row, 1)
                        _emit(x, y, z, ny, nz, ins[1], outs[0], offsets, lut, edges, dirs, row, 2)
                        row += 1


def tetra_cut_edges(values: np.ndarray, isovalue: float):
    """(tri_edges (T, 3) int64, outward (T, 3) float64); T = 0 when nothing is cut."""
    values = np.ascontiguousarray(values, dtype=np.float64)
    counts = _count(values, float(isovalue), _OFFSETS, _TETS)
    starts = np.zeros(counts.size, dtype=np.int64)
    if counts.size > 1:
        starts[1:] = np.cumsum(counts)[:-1]
    total = int(counts.sum())
    edges = np.empty((total, 3), dtype=np.int64)
    dirs = np.zeros((total, 3), dtype=np.int64)
    if total:
        _fill(values, float(isovalue), _OFFSETS, _TETS, _LUT, starts, edges, dirs)
    return edges, dirs.astype(float)
