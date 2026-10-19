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
