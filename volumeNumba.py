"""Optional numba kernel for isosurfaces of volumetric grids (the Cube tab).

``tetra_cut_edges`` is the compiled equivalent of
``molecularSurface.tetra_cut_edges``: the cell scan of marching tetrahedra.
It produces the same triangles -- the same cut edges, in the same vertex
order, with the same outward directions -- but scans the cells in parallel
without the (cells x 8) boolean stack the numpy path builds, so it is both
faster and much lighter on memory for large cubes. Triangles come out in
cell order rather than tetrahedron order; ``mesh_from_cut_edges`` sorts the
edge ids, so the vertex arrays of both paths are identical.

``format_e13`` writes cube-file values as text, the ``%13.5E`` layout of
``cubeFile`` (6 per line, a new line for every (x, y) column), in parallel.
It produces the bytes Python's own formatting produces: the few values it
cannot settle exactly (a sixth significant digit within 1e-7 of a rounding
tie, non-finite values, exponents beyond +/-300) are flagged and formatted
by Python instead.

Importing this module raises ImportError when numba is not installed; the
callers (volumeData, cubeFile) then keep their Python/numpy paths.
"""

from __future__ import annotations

import math

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


# ----------------------------------------------------------------------
# Cube values as text: %13.5E, 6 per line, a new line for every column
# ----------------------------------------------------------------------
E13_WIDTH = 13
VALUES_PER_LINE = 6
# A sixth significant digit this close to a rounding tie goes to Python: the
# scaled value carries ~1e-9 absolute error, the tie test needs it exact.
TIE_TOLERANCE = 1e-7
MAX_DECIMAL_EXPONENT = 300  # beyond it 10**(5 - e) leaves the normal range


@numba.njit(cache=True)
def _scaled_mantissa(a, e):
    """a / 10**(e - 5): six significant digits before the decimal point."""
    k = 5 - e
    if k >= 0:
        return a * 10.0 ** k
    return a / 10.0 ** (-k)


@numba.njit(cache=True)
def _write_e13(out, pos, v):
    """Write v as Python's f"{v:13.5E}" into out[pos:pos + 13]; False when Python must do it."""
    if not math.isfinite(v):
        return False
    negative = math.copysign(1.0, v) < 0.0
    a = abs(v)
    digits = 0
    e = 0
    if a != 0.0:
        e = int(math.floor(math.log10(a)))
        if e > MAX_DECIMAL_EXPONENT or e < -MAX_DECIMAL_EXPONENT:
            return False
        m = _scaled_mantissa(a, e)
        if m < 100000.0:
            e -= 1
            m = _scaled_mantissa(a, e)
        elif m >= 1000000.0:
            e += 1
            m = _scaled_mantissa(a, e)
        if m < 100000.0 or m >= 1000000.0:
            return False
        whole = math.floor(m)
        fraction = m - whole
        if abs(fraction - 0.5) < TIE_TOLERANCE:
            return False
        digits = int(whole) + (1 if fraction > 0.5 else 0)
        if digits == 1000000:  # 9.999996 -> 1.00000E+01
            digits = 100000
            e += 1
    j = pos + E13_WIDTH
    exponent = abs(e)
    n_exp = 3 if exponent >= 100 else 2
    for _ in range(n_exp):
        j -= 1
        out[j] = 48 + exponent % 10
        exponent //= 10
    j -= 1
    out[j] = 45 if e < 0 else 43   # - or +
    j -= 1
    out[j] = 69                    # E
    for _ in range(5):
        j -= 1
        out[j] = 48 + digits % 10
        digits //= 10
    j -= 1
    out[j] = 46                    # .
    j -= 1
    out[j] = 48 + digits
    if negative:
        j -= 1
        out[j] = 45
    for k in range(pos, j):
        out[k] = 32
    return True


def column_bytes(nz: int) -> int:
    """Bytes of one (x, y) column: nz values plus one newline per line of up to six."""
    return nz * E13_WIDTH + (nz + VALUES_PER_LINE - 1) // VALUES_PER_LINE


@numba.njit(parallel=True, cache=True)
def _format_columns(columns, per_column, out, python_needed):
    n_columns, nz = columns.shape
    for c in numba.prange(n_columns):
        pos = c * per_column
        for k in range(nz):
            if not _write_e13(out, pos, columns[c, k]):
                python_needed[c, k] = True
            pos += E13_WIDTH
            if k % VALUES_PER_LINE == VALUES_PER_LINE - 1 or k == nz - 1:
                out[pos] = 10
                pos += 1


def format_e13(columns: np.ndarray) -> str:
    """(n_columns, nz) values as cube text, byte for byte what the Python loop writes."""
    columns = np.ascontiguousarray(columns, dtype=np.float64)
    n_columns, nz = columns.shape
    per_column = column_bytes(nz)
    out = np.empty(n_columns * per_column, dtype=np.uint8)
    python_needed = np.zeros(columns.shape, dtype=np.bool_)
    if columns.size:
        _format_columns(columns, per_column, out, python_needed)
    for c, k in np.argwhere(python_needed):
        start = c * per_column + k * E13_WIDTH + k // VALUES_PER_LINE
        out[start:start + E13_WIDTH] = np.frombuffer(f"{columns[c, k]:13.5E}".encode("ascii"), np.uint8)
    return out.tobytes().decode("ascii")
