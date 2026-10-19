"""Gaussian cube files (.cub/.cube): write and read volumetric grid data.

Toga-free. Lengths are bohr, as in the cube format itself (a positive voxel
count on lines 4-6 means bohr). ``read_cube``/``write_cube`` handle only
orthogonal grids whose three axes are the Cartesian axes with one common
spacing -- every grid the Fukui tool makes, and what Multiwfn writes by
default. ``read_cube_volumes``/``write_cube_volume`` (the 3D viewer's Cube
tab) accept any cube: orbital cubes, Angstrom headers, skewed cells and
unequal spacings. Values are written 6 per line in %13.5E, z fastest, a new
line for every (x, y) column: the layout Gaussian's cubegen and Multiwfn use.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

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


# ----------------------------------------------------------------------
# General cubes (the 3D viewer's Cube tab). read_cube above stays strict
# for the Fukui tool; this reader accepts what other programs write:
# orbital cubes (negative atom count, one volume per orbital), Angstrom
# headers (negative voxel counts), skewed cells and unequal spacings.
# ----------------------------------------------------------------------
BOHR_RADIUS_ANGSTROM = 0.529177210544  # coordinateConverter imports Toga, so the constant is repeated
SINGULAR_AXES_VOLUME = 1e-12  # bohr^3
_FORTRAN_EXPONENT = re.compile(r"(?<=[0-9.])([+-][0-9]{2,3})(?=\s|$)")


@dataclass
class CubeVolume:
    title: str
    comment: str
    atomic_numbers: np.ndarray   # (n,)
    nuclear_charges: np.ndarray  # (n,)
    coordinates: np.ndarray      # (n, 3) bohr
    origin: np.ndarray           # (3,) bohr
    axes: np.ndarray             # (3, 3) bohr; row k = step of grid index k
    values: np.ndarray           # (nx, ny, nz)
    orbital: Optional[int] = None  # MO number when read from an orbital cube


def _parse_values(text: str, name: str) -> np.ndarray:
    try:
        return np.array(text.split(), dtype=np.float64)
    except ValueError:
        pass
    # Some Fortran writers drop the E of a three-digit exponent: 1.234-100.
    repaired = _FORTRAN_EXPONENT.sub(r"E\1", text)
    try:
        return np.array(repaired.split(), dtype=np.float64)
    except ValueError:
        for token in repaired.split():
            try:
                float(token)
            except ValueError:
                raise ValueError(f"{name}: non-numeric value '{token}' in the volumetric data.") from None
        raise


def read_cube_volumes(path) -> List[CubeVolume]:
    """Every volume in a cube file (one per orbital for an orbital cube), lengths in bohr."""
    path = Path(path)
    try:
        lines = path.read_text(encoding="ascii", errors="replace").splitlines()
    except OSError as exc:
        raise ValueError(f"Cannot read {path}: {exc}") from exc
    name = path.name
    try:
        head = lines[2].split()
        natoms = int(head[0])
        origin = np.array([float(x) for x in head[1:4]])
        counts, axes = [], []
        for k in range(3):
            t = lines[3 + k].split()
            counts.append(int(t[0]))
            axes.append([float(x) for x in t[1:4]])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"{name}: malformed cube header (lines 3-6).") from exc
    if 0 in counts:
        raise ValueError(f"{name}: a grid dimension is zero.")
    shape = tuple(abs(n) for n in counts)
    axes = np.array(axes)
    n_atoms = abs(natoms)
    try:
        atoms = [lines[6 + i].split() for i in range(n_atoms)]
        z = np.array([int(float(a[0])) for a in atoms], dtype=int)
        q = np.array([float(a[1]) for a in atoms])
        xyz = np.array([[float(c) for c in a[2:5]] for a in atoms]).reshape(-1, 3)
    except (IndexError, ValueError) as exc:
        raise ValueError(f"{name}: malformed atom lines.") from exc
    if any(n < 0 for n in counts):  # Gaussian convention: negative voxel counts mean Angstrom
        factor = 1.0 / BOHR_RADIUS_ANGSTROM
        origin, axes, xyz = origin * factor, axes * factor, xyz * factor
    if abs(np.linalg.det(axes)) < SINGULAR_AXES_VOLUME:
        raise ValueError(f"{name}: the three grid axes are singular (they do not span a volume).")

    start = 6 + n_atoms
    orbitals: List[Optional[int]] = [None]
    if natoms < 0:  # orbital cube: NMO and the MO numbers, possibly over several lines
        ids: List[int] = []
        try:
            while not ids or len(ids) < ids[0] + 1:
                ids.extend(int(t) for t in lines[start].split())
                start += 1
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{name}: malformed orbital line after the atoms.") from exc
        if ids[0] <= 0 or len(ids) != ids[0] + 1:
            raise ValueError(f"{name}: malformed orbital line after the atoms.")
        orbitals = list(ids[1:])

    data = _parse_values("\n".join(lines[start:]), name)
    expected = shape[0] * shape[1] * shape[2] * len(orbitals)
    if data.size != expected:
        raise ValueError(f"{name}: {data.size} values, expected {expected}; the file looks truncated.")
    data = data.reshape(shape + (len(orbitals),))
    return [CubeVolume(lines[0], lines[1], z, q, xyz.copy(), origin.copy(), axes.copy(),
                       np.ascontiguousarray(data[..., k]), mo)
            for k, mo in enumerate(orbitals)]


def format_cube_volume(volume: CubeVolume) -> str:
    """The format_cube layout with general axes; an orbital is written as a plain cube."""
    nx, ny, nz = volume.values.shape
    comment = volume.comment.replace("\n", " ")
    if volume.orbital is not None and f"MO {volume.orbital}" not in comment:
        comment = f"{comment} MO {volume.orbital}".strip()
    out = [volume.title.replace("\n", " "), comment,
           f"{len(volume.atomic_numbers):5d}" + "".join(f"{c:12.6f}" for c in volume.origin)]
    for n, axis in zip((nx, ny, nz), np.asarray(volume.axes, dtype=float)):
        out.append(f"{n:5d}" + "".join(f"{c:12.6f}" for c in axis))
    for z, q, xyz in zip(volume.atomic_numbers, volume.nuclear_charges, volume.coordinates):
        out.append(f"{int(z):5d}{float(q):12.6f}" + "".join(f"{c:12.6f}" for c in xyz))
    for column in np.asarray(volume.values, dtype=float).reshape(nx * ny, nz):
        for k in range(0, nz, 6):
            out.append("".join(f"{v:13.5E}" for v in column[k:k + 6]))
    return "\n".join(out) + "\n"


def write_cube_volume(path, volume: CubeVolume) -> str:
    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(format_cube_volume(volume), encoding="ascii")
    except OSError as exc:
        raise ValueError(f"Could not write {path}: {exc}") from exc
    return str(path)
