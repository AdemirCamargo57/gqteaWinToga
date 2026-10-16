"""Gaussian cube files (.cub/.cube): write and read volumetric grid data.

Toga-free. Lengths are bohr, as in the cube format itself (a positive voxel
count on lines 4-6 means bohr). Only orthogonal grids whose three axes are
the Cartesian axes with one common spacing are written or read -- that is
what every grid in this project is, and what Multiwfn writes by default.
Values are written 6 per line in %13.5E, z fastest, a new line for every
(x, y) column: the layout Gaussian's cubegen and Multiwfn use.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

AXIS_TOLERANCE = 1e-6  # bohr


@dataclass
class CubeData:
    title: str
    comment: str
    atomic_numbers: np.ndarray   # (n,)
    nuclear_charges: np.ndarray  # (n,)
    coordinates: np.ndarray      # (n, 3) bohr
    origin: np.ndarray           # (3,) bohr
    spacing: float               # bohr
    values: np.ndarray           # (nx, ny, nz)


def format_cube(cube: CubeData) -> str:
    nx, ny, nz = cube.values.shape
    out = [cube.title.replace("\n", " "), cube.comment.replace("\n", " "),
           f"{len(cube.atomic_numbers):5d}" + "".join(f"{c:12.6f}" for c in cube.origin)]
    for n, axis in zip((nx, ny, nz), np.eye(3) * cube.spacing):
        out.append(f"{n:5d}" + "".join(f"{c:12.6f}" for c in axis))
    for z, q, xyz in zip(cube.atomic_numbers, cube.nuclear_charges, cube.coordinates):
        out.append(f"{int(z):5d}{float(q):12.6f}" + "".join(f"{c:12.6f}" for c in xyz))
    for column in np.asarray(cube.values, dtype=float).reshape(nx * ny, nz):
        for k in range(0, nz, 6):
            out.append("".join(f"{v:13.5E}" for v in column[k:k + 6]))
    return "\n".join(out) + "\n"


def write_cube(path, cube: CubeData) -> str:
    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(format_cube(cube), encoding="ascii")
    except OSError as exc:
        raise ValueError(f"Could not write {path}: {exc}") from exc
    return str(path)


def read_cube(path) -> CubeData:
    """Parse a cube written by this module, cubegen or Multiwfn. ValueError on any problem."""
    path = Path(path)
    try:
        lines = path.read_text(encoding="ascii", errors="replace").splitlines()
    except OSError as exc:
        raise ValueError(f"Cannot read {path}: {exc}") from exc
    try:
        head = lines[2].split()
        natoms = int(head[0])
        origin = np.array([float(x) for x in head[1:4]])
        shape, axes = [], []
        for k in range(3):
            t = lines[3 + k].split()
            shape.append(int(t[0]))
            axes.append([float(x) for x in t[1:4]])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"{path.name}: malformed cube header.") from exc
    if natoms < 0:
        raise ValueError(f"{path.name}: orbital cubes (negative atom count) are not supported.")
    axes = np.array(axes)
    if min(shape) <= 0:
        raise ValueError(f"{path.name}: grid sizes must be positive (bohr units).")
    spacing = float(axes[0, 0])
    if np.abs(axes - np.eye(3) * spacing).max() > AXIS_TOLERANCE:
        raise ValueError(f"{path.name}: only axis-aligned grids with one spacing are supported.")
    atoms = [lines[6 + i].split() for i in range(natoms)]
    try:
        z = np.array([int(a[0]) for a in atoms])
        q = np.array([float(a[1]) for a in atoms])
        xyz = np.array([[float(c) for c in a[2:5]] for a in atoms]).reshape(-1, 3)
        data = np.array([float(v) for line in lines[6 + natoms:] for v in line.split()])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"{path.name}: malformed atom or data lines.") from exc
    expected = shape[0] * shape[1] * shape[2]
    if data.size != expected:
        raise ValueError(f"{path.name}: {data.size} values, expected {expected}; the file looks truncated.")
    return CubeData(lines[0], lines[1], z, q, xyz, origin, spacing, data.reshape(shape))
