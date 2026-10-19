"""Numerics of the 3D viewer's Cube tab: isosurfaces, slices and grid arithmetic.

Toga-free and GL-free, so all of it is tested headless. Volumes are
``cubeFile.CubeVolume`` objects with lengths in bohr; the grid point
(i, j, k) sits at ``origin + (i, j, k) @ axes``, so skewed cells and unequal
spacings need no special case.

Isosurfaces reuse molecularSurface's marching tetrahedra. The cell scan runs
in grid-index space (origin 0, spacing 1) -- through the compiled
``volumeNumba`` kernel when numba is installed, otherwise numpy, with
identical vertices -- and the vertices are then mapped through the cube's
axes. Linear interpolation along an edge commutes with that affine map, so
the result is exact for any lattice.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

import molecularSurface as ms
from cubeFile import BOHR_RADIUS_ANGSTROM, CubeVolume
from help import AtomicData

try:  # optional compiled cell scan; same vertices as the numpy path
    import volumeNumba
except ImportError:  # pragma: no cover - depends on the installation
    volumeNumba = None


def numba_available() -> bool:
    return volumeNumba is not None


def _resolve_numba(use_numba: Optional[bool]) -> bool:
    """None = "if installed" (the gaussianWavefunction convention)."""
    if use_numba is None:
        return numba_available()
    if use_numba and not numba_available():
        raise ValueError("use_numba=True but numba is not installed (pip install numba).")
    return bool(use_numba)


# ----------------------------------------------------------------------
# Isosurfaces
# ----------------------------------------------------------------------
@dataclass
class IsoMesh:
    vertices: np.ndarray   # (V, 3) bohr
    triangles: np.ndarray  # (F, 3) int64, oriented towards decreasing sign*value
    normals: np.ndarray    # (V, 3) unit, area-weighted


def vertex_normals(vertices: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    a, b, c = (vertices[triangles[:, k]] for k in range(3))
    face = np.cross(b - a, c - a)  # length = 2 x area, so the sum is area-weighted
    normals = np.zeros_like(vertices)
    for k in range(3):
        for d in range(3):
            normals[:, d] += np.bincount(triangles[:, k], weights=face[:, d], minlength=len(vertices))
    length = np.linalg.norm(normals, axis=1)
    return normals / np.where(length > 0, length, 1.0)[:, None]


def isosurface(volume: CubeVolume, isovalue: float, sign: int = 1,
               use_numba: Optional[bool] = None) -> Optional[IsoMesh]:
    """The surface ``sign * value = isovalue``; None when the isovalue cuts nothing.

    ``sign = -1`` gives the negative lobe of a signed function (an orbital,
    a density difference) at -isovalue.
    """
    if not np.isfinite(isovalue):
        raise ValueError("The isovalue must be a finite number.")
    if sign not in (1, -1):
        raise ValueError("sign must be +1 or -1.")
    values = np.ascontiguousarray(volume.values if sign == 1 else -volume.values, dtype=np.float64)
    if _resolve_numba(use_numba):
        edges, outward = volumeNumba.tetra_cut_edges(values, isovalue)
    else:
        try:
            edges, outward = ms.tetra_cut_edges(values, isovalue)
        except ValueError:  # no cell is cut
            return None
    if len(edges) == 0:
        return None
    index_vertices, triangles = ms.mesh_from_cut_edges(edges, outward, values, np.zeros(3), 1.0, isovalue)
    axes = np.asarray(volume.axes, dtype=float)
    vertices = np.asarray(volume.origin, dtype=float) + index_vertices @ axes
    if np.linalg.det(axes) < 0:  # a left-handed lattice mirrors the winding
        triangles = triangles[:, [0, 2, 1]]
    triangles = np.ascontiguousarray(triangles, dtype=np.int64)
    return IsoMesh(vertices, triangles, vertex_normals(vertices, triangles))


# ----------------------------------------------------------------------
# Grid compatibility and arithmetic
# ----------------------------------------------------------------------
GRID_TOLERANCE = 1e-4  # bohr; cube headers carry 6 decimals
OPERATIONS = ("add", "subtract")
RESULT_COMMENT = "Made by the gQTEA 3D Molecular Viewer (Cube tab)"


def _shape_text(shape) -> str:
    return " x ".join(str(n) for n in shape)


def grid_mismatch(a: CubeVolume, b: CubeVolume) -> Optional[str]:
    """None when the two grids are point-for-point the same, otherwise why not."""
    if a.values.shape != b.values.shape:
        return (f"The grids have different sizes: {_shape_text(a.values.shape)} and "
                f"{_shape_text(b.values.shape)} points.")
    shift = float(np.abs(np.asarray(a.origin) - np.asarray(b.origin)).max())
    if shift > GRID_TOLERANCE:
        return f"The grid origins differ by up to {shift:.6f} bohr."
    axes = float(np.abs(np.asarray(a.axes) - np.asarray(b.axes)).max())
    if axes > GRID_TOLERANCE:
        return f"The grid axes (step vectors) differ by up to {axes:.6f} bohr."
    return None


def atoms_mismatch(a: CubeVolume, b: CubeVolume) -> Optional[str]:
    """None when both cubes describe the same atoms at the same positions."""
    if not np.array_equal(np.asarray(a.atomic_numbers), np.asarray(b.atomic_numbers)):
        return "The two cubes list different atoms."
    if len(a.coordinates) and float(np.abs(np.asarray(a.coordinates) - np.asarray(b.coordinates)).max()) > GRID_TOLERANCE:
        return "The two cubes place their atoms at different positions."
    return None


def _result(a: CubeVolume, values: np.ndarray, title: str) -> CubeVolume:
    return CubeVolume(title, RESULT_COMMENT, np.array(a.atomic_numbers).copy(),
                      np.array(a.nuclear_charges).copy(), np.array(a.coordinates).copy(),
                      np.array(a.origin, dtype=float).copy(), np.array(a.axes, dtype=float).copy(),
                      values, None)


def combine(a: CubeVolume, b: CubeVolume, op: str, title: str = "") -> CubeVolume:
    """A + B or A - B on identical grids; the result keeps A's atoms and grid."""
    if op not in OPERATIONS:
        raise ValueError(f"Unknown grid operation '{op}'; expected one of {OPERATIONS}.")
    reason = grid_mismatch(a, b)
    if reason:
        raise ValueError(f"The cubes are not on the same grid. {reason}")
    values = a.values + b.values if op == "add" else a.values - b.values
    return _result(a, values, title or ("Sum of two cubes" if op == "add" else "Difference of two cubes"))


def scale(a: CubeVolume, factor: float, title: str = "") -> CubeVolume:
    factor = float(factor)
    if not np.isfinite(factor):
        raise ValueError("The scale factor must be a finite number.")
    return _result(a, a.values * factor, title or f"Cube scaled by {factor:g}")


# ----------------------------------------------------------------------
# Statistics and defaults
# ----------------------------------------------------------------------
def voxel_volume(volume: CubeVolume) -> float:
    return float(abs(np.linalg.det(np.asarray(volume.axes, dtype=float))))


def integral(volume: CubeVolume) -> float:
    """Riemann sum of the values over the grid (electrons, for a density in e/bohr^3)."""
    return float(volume.values.sum()) * voxel_volume(volume)


def value_statistics(volume: CubeVolume) -> Tuple[float, float, float]:
    return float(volume.values.min()), float(volume.values.max()), integral(volume)


ORBITAL_ISOVALUE = 0.02    # GaussView's default for orbitals
DENSITY_ISOVALUE = 0.0004  # GaussView's default for total densities
SIGNED_FRACTION = 0.01     # negative values beyond 1 % of max|v| make a cube "signed"


def default_isovalue(volume: CubeVolume) -> float:
    """A starting isovalue: orbital/signed data 0.02, densities 0.0004, capped below max|v|."""
    max_abs = float(np.abs(volume.values).max()) if volume.values.size else 0.0
    if max_abs == 0.0:
        return ORBITAL_ISOVALUE
    signed = volume.orbital is not None or float(volume.values.min()) < -SIGNED_FRACTION * max_abs
    value = ORBITAL_ISOVALUE if signed else DENSITY_ISOVALUE
    if value >= max_abs:
        value = 0.5 * max_abs
    return float(f"{value:.0e}")


# ----------------------------------------------------------------------
# Slices
# ----------------------------------------------------------------------
AXIS_NAMES = ("a", "b", "c")


def slice_plane(volume: CubeVolume, axis: int, index: int) -> Tuple[np.ndarray, np.ndarray]:
    """Points (n1, n2, 3) bohr and values (n1, n2) of lattice plane ``index`` along ``axis``."""
    if axis not in (0, 1, 2):
        raise ValueError(f"The slice axis must be 0, 1 or 2 (a, b or c), not {axis}.")
    n = volume.values.shape[axis]
    if not 0 <= index < n:
        raise ValueError(f"The slice plane must be between 1 and {n} along axis "
                         f"{AXIS_NAMES[axis]}, not {index + 1}.")
    u, w = [d for d in range(3) if d != axis]
    axes = np.asarray(volume.axes, dtype=float)
    shape = volume.values.shape
    base = np.asarray(volume.origin, dtype=float) + index * axes[axis]
    i = np.arange(shape[u])[:, None, None]
    j = np.arange(shape[w])[None, :, None]
    points = base + i * axes[u] + j * axes[w]
    return points, np.take(volume.values, index, axis=axis)


def slice_triangles(n1: int, n2: int) -> np.ndarray:
    """Two triangles per cell of an (n1, n2) point grid, indices in row-major order."""
    if n1 < 2 or n2 < 2:
        return np.zeros((0, 3), dtype=np.int64)
    i, j = np.meshgrid(np.arange(n1 - 1), np.arange(n2 - 1), indexing="ij")
    v00 = (i * n2 + j).ravel()
    v10, v01, v11 = v00 + n2, v00 + 1, v00 + n2 + 1
    return np.ascontiguousarray(np.vstack([np.column_stack([v00, v10, v11]),
                                           np.column_stack([v00, v11, v01])]), dtype=np.int64)


def slice_table_text(volume: CubeVolume, axis: int, index: int, name: str = "") -> str:
    points, values = slice_plane(volume, axis, index)
    shape = values.shape
    lines = [f"# Cube slice {name}".rstrip(),
             f"# plane {index + 1} of {volume.values.shape[axis]} along grid axis {AXIS_NAMES[axis]}",
             f"# points {shape[0]} x {shape[1]}",
             "# x_angstrom y_angstrom z_angstrom value"]
    xyz = points.reshape(-1, 3) * BOHR_RADIUS_ANGSTROM
    for (x, y, z), v in zip(xyz, values.ravel()):
        lines.append(f"{x:14.6f} {y:14.6f} {z:14.6f} {v:16.8E}")
    return "\n".join(lines) + "\n"


COLOUR_RANGE_MODES = ("symmetric", "minmax")


def colour_range(values: np.ndarray, mode: str) -> Tuple[float, float]:
    """(lo, hi) for slice colouring: symmetric +/-max|v| or the data's min..max."""
    values = np.asarray(values, dtype=float)
    if mode == "symmetric":
        r = float(np.abs(values).max()) if values.size else 0.0
        r = r if r > 0 else 1.0
        return -r, r
    if mode == "minmax":
        return float(values.min()), float(values.max())
    raise ValueError(f"Unknown colour range '{mode}'; expected one of {COLOUR_RANGE_MODES}.")


# ----------------------------------------------------------------------
# Colours (pure; espViewer imports GLFW and molecularDesign imports Toga,
# so their small colour helpers are not imported here)
# ----------------------------------------------------------------------
COLOUR_PRESETS = {
    "Blue": (0.0, 0.35, 1.0),
    "Red": (0.9, 0.1, 0.1),
    "Green": (0.1, 0.75, 0.2),
    "Yellow": (1.0, 0.85, 0.0),
    "Orange": (1.0, 0.55, 0.0),
    "Cyan": (0.0, 0.8, 0.9),
    "Magenta": (0.85, 0.1, 0.85),
    "Purple": (0.5, 0.2, 0.8),
    "White": (1.0, 1.0, 1.0),
    "Grey": (0.6, 0.6, 0.6),
}


def parse_colour(text: str) -> Tuple[float, float, float]:
    """A preset name (any case) or #rrggbb -> RGB in [0, 1]."""
    value = (text or "").strip()
    for name, rgb in COLOUR_PRESETS.items():
        if value.lower() == name.lower():
            return rgb
    if len(value) == 7 and value.startswith("#"):
        try:
            return tuple(int(value[k:k + 2], 16) / 255.0 for k in (1, 3, 5))
        except ValueError:
            pass
    raise ValueError(f"Unknown colour '{value}': use one of {', '.join(COLOUR_PRESETS)} or #rrggbb.")


def diverging_colours(values: np.ndarray, half_range: float) -> np.ndarray:
    """Blue (-R) -> white (0) -> red (+R), clipped; RGB (N, 3)."""
    t = np.clip(np.asarray(values, dtype=float) / (half_range if half_range > 0 else 1.0), -1.0, 1.0)
    pos, neg = np.clip(t, 0.0, 1.0), np.clip(-t, 0.0, 1.0)
    return np.column_stack([1.0 - neg, 1.0 - pos - neg, 1.0 - pos])


def rainbow_colours(values: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """GaussView hue sweep red (lo) -> green -> blue (hi); RGB (N, 3)."""
    values = np.asarray(values, dtype=float).ravel()
    fraction = np.full(values.shape, 0.5) if hi <= lo else (values - lo) / (hi - lo)
    h = np.clip(fraction, 0.0, 1.0) * 4.0
    seg = np.minimum(h.astype(int), 3)
    x = h - seg
    one, zero = np.ones_like(x), np.zeros_like(x)
    choices = [np.column_stack(c) for c in ((one, x, zero), (1 - x, one, zero),
                                            (zero, one, x), (zero, 1 - x, one))]
    return np.select([(seg == k)[:, None] for k in range(4)], choices)


# ----------------------------------------------------------------------
# Structure
# ----------------------------------------------------------------------
UNKNOWN_ELEMENT = "X"


def structure_frame(volume: CubeVolume):
    """The cube's atoms as a viewer frame: [(symbol, (x, y, z) Angstrom), ...]."""
    symbols = {z: s for s, z in AtomicData.atomic_numbers.items()}
    xyz = np.asarray(volume.coordinates, dtype=float).reshape(-1, 3) * BOHR_RADIUS_ANGSTROM
    return [(symbols.get(int(z), UNKNOWN_ELEMENT), tuple(float(c) for c in p))
            for z, p in zip(volume.atomic_numbers, xyz)]
