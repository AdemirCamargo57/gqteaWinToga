"""Fast clash checking for the random-insertion solvent-box builders.

The builders accept a candidate solvent molecule only if none of its atoms
comes closer than ``vdw_scale * (r_i + r_j)`` to an atom of the solute or of a
molecule accepted earlier, optionally under the orthorhombic minimum image.
The original code did that with pure-Python atom loops over **every** accepted
molecule on every attempt; profiling a 20 A water box showed ~90 % of the run
was interpreter overhead (``round()`` on numpy scalars alone: 39 M calls).

This module keeps the *decision* and changes only how it is reached:

* :class:`MoleculeGrid` -- accepted molecules in a cell list keyed by their
  centre, so a candidate is compared only with molecules in its own and the
  neighbouring cells instead of all of them;
* :func:`atoms_clash` -- the atom-pair test vectorised with numpy, or with the
  optional compiled kernel in :mod:`solventPackingNumba`.

**Every accept/reject decision is bit-identical to the legacy loops**, so a
fixed random seed still produces exactly the same box:

* the legacy centre pre-filter ``d_centre > r_cand + r_exist + 2.5`` is applied
  to the cell-list candidates unchanged (the cell list only returns a
  superset of the molecules it would keep);
* the minimum image is ``dx - L * round(dx / L)``; Python's ``round`` and
  ``np.round`` both round half to even and both are exact, so the products and
  differences are the same float64 operations;
* ``dist2 = dx*dx + dy*dy + dz*dz`` is written out in that order, and the
  cutoff is ``scale * (r_i + r_j)`` compared as ``dist2 < cutoff*cutoff`` (and
  ``dist2 < (cutoff + extra)**2`` when an extra clearance is requested).

Tests in tests/test_solventPacking.py compare whole boxes against the legacy
implementation. No Toga here.
"""
import math
from typing import List, Optional, Sequence

import numpy as np

try:  # optional compiled kernel (pip install numba)
    import solventPackingNumba
except Exception:  # pragma: no cover - depends on the environment
    solventPackingNumba = None

# The legacy builders add this margin (in Angstrom) to the two bounding radii
# when deciding whether two molecules are close enough to need an atom check.
CENTRE_FILTER_MARGIN = 2.5

ENGINE_AUTO = "auto"
ENGINE_NUMPY = "numpy"
ENGINE_NUMBA = "numba"
ENGINE_LEGACY = "legacy"
ENGINES = (ENGINE_AUTO, ENGINE_NUMPY, ENGINE_NUMBA, ENGINE_LEGACY)


def numba_available() -> bool:
    return solventPackingNumba is not None


def resolve_engine(engine: str) -> str:
    """``auto`` becomes ``numba`` when it imports, otherwise ``numpy``."""
    if engine not in ENGINES:
        raise ValueError(f"Unknown packing engine '{engine}'.")
    if engine == ENGINE_AUTO:
        return ENGINE_NUMBA if numba_available() else ENGINE_NUMPY
    if engine == ENGINE_NUMBA and not numba_available():
        raise ValueError("The numba engine was requested but numba is not installed.")
    return engine


def atoms_clash_numpy(cand_xyz: np.ndarray, cand_r: np.ndarray,
                      other_xyz: np.ndarray, other_r: np.ndarray,
                      vdw_scale: float, lengths: np.ndarray, periodic: bool,
                      extra: float = 0.0) -> bool:
    """True if any candidate atom overlaps any other atom (legacy rule)."""
    if len(other_xyz) == 0 or len(cand_xyz) == 0:
        return False
    dx = cand_xyz[:, None, 0] - other_xyz[None, :, 0]
    dy = cand_xyz[:, None, 1] - other_xyz[None, :, 1]
    dz = cand_xyz[:, None, 2] - other_xyz[None, :, 2]
    if periodic:
        a, b, c = float(lengths[0]), float(lengths[1]), float(lengths[2])
        if a > 0.0:
            dx = dx - a * np.round(dx / a)
        if b > 0.0:
            dy = dy - b * np.round(dy / b)
        if c > 0.0:
            dz = dz - c * np.round(dz / c)
    dist2 = dx * dx + dy * dy + dz * dz
    cutoff = vdw_scale * (cand_r[:, None] + other_r[None, :])
    if np.any(dist2 < cutoff * cutoff):
        return True
    if extra > 0.0:
        expanded = cutoff + extra
        return bool(np.any(dist2 < expanded * expanded))
    return False


def atoms_clash(cand_xyz, cand_r, other_xyz, other_r, vdw_scale, lengths, periodic,
                extra=0.0, engine=ENGINE_NUMPY) -> bool:
    if engine == ENGINE_NUMBA:
        if len(other_xyz) == 0 or len(cand_xyz) == 0:
            return False
        return bool(solventPackingNumba.atoms_clash(
            cand_xyz, cand_r, other_xyz, other_r, float(vdw_scale),
            np.asarray(lengths, dtype=float), bool(periodic), float(extra)))
    return atoms_clash_numpy(cand_xyz, cand_r, other_xyz, other_r, vdw_scale, lengths, periodic, extra)


class MoleculeGrid:
    """Accepted molecules in a uniform cell list over the box.

    Cells are at least ``min_cell`` wide, where ``min_cell`` is the largest
    centre-filter distance any query will use, so the 27 cells around a
    candidate always contain every molecule the legacy filter would keep.
    With fewer than three cells along an axis the neighbour set simply covers
    the whole axis (duplicates removed); when that is true of all three axes
    -- every box under about three filter distances -- the cell look-up is
    skipped and all molecules are candidates.

    Centres, radii and atoms live in contiguous, growing numpy arrays, so a
    query gathers its neighbours by index instead of rebuilding arrays from
    Python lists (which profiling showed was most of the remaining cost).
    """

    def __init__(self, lengths: Sequence[float], periodic: bool, min_cell: float):
        self.lengths = np.asarray(lengths, dtype=float)
        self.periodic = bool(periodic)
        min_cell = max(float(min_cell), 1e-6)
        self.shape = tuple(max(1, int(math.floor(length / min_cell))) for length in self.lengths)
        self.widths = self.lengths / np.array(self.shape, dtype=float)
        # -1/0/+1 neighbours reach every cell when an axis has at most 3 cells
        # (periodic, indices wrap) or at most 2 (non-periodic).
        reach = 3 if self.periodic else 2
        self.covers_all = all(n <= reach for n in self.shape)
        self.cells = {}
        self._neighbour_cache = {}
        self.count = 0
        self.atom_count = 0
        self.centres = np.zeros((64, 3))
        self.radii = np.zeros(64)
        self.mol_start = np.zeros(64, dtype=np.int64)
        self.mol_count = np.zeros(64, dtype=np.int64)
        self.atoms = np.zeros((256, 3))
        self.atom_r = np.zeros(256)

    def __len__(self):
        return self.count

    def _cell_of(self, centre) -> tuple:
        key = []
        for axis in range(3):
            value = float(centre[axis])
            if self.periodic:
                value = value % self.lengths[axis]
            index = int(math.floor(value / self.widths[axis]))
            key.append(min(self.shape[axis] - 1, max(0, index)))
        return tuple(key)

    def _neighbour_cells(self, key) -> List[tuple]:
        cached = self._neighbour_cache.get(key)
        if cached is not None:
            return cached
        axes = []
        for axis in range(3):
            n = self.shape[axis]
            values = set()
            for offset in (-1, 0, 1):
                index = key[axis] + offset
                if self.periodic:
                    values.add(index % n)
                elif 0 <= index < n:
                    values.add(index)
            axes.append(sorted(values))
        cells = [(i, j, k) for i in axes[0] for j in axes[1] for k in axes[2]]
        self._neighbour_cache[key] = cells
        return cells

    @staticmethod
    def _grow(array, needed):
        if needed <= len(array):
            return array
        size = len(array)
        while size < needed:
            size *= 2
        grown = np.zeros((size,) + array.shape[1:], dtype=array.dtype)
        grown[: len(array)] = array
        return grown

    def add(self, centre, radius: float, coords: np.ndarray, atom_radii: np.ndarray) -> int:
        index = self.count
        n = len(coords)
        self.centres = self._grow(self.centres, index + 1)
        self.radii = self._grow(self.radii, index + 1)
        self.mol_start = self._grow(self.mol_start, index + 1)
        self.mol_count = self._grow(self.mol_count, index + 1)
        self.atoms = self._grow(self.atoms, self.atom_count + n)
        self.atom_r = self._grow(self.atom_r, self.atom_count + n)
        self.centres[index] = centre
        self.radii[index] = radius
        self.mol_start[index] = self.atom_count
        self.mol_count[index] = n
        self.atoms[self.atom_count: self.atom_count + n] = coords
        self.atom_r[self.atom_count: self.atom_count + n] = atom_radii
        self.atom_count += n
        self.count += 1
        if not self.covers_all:
            self.cells.setdefault(self._cell_of(centre), []).append(index)
        return index

    def candidates(self, centre) -> np.ndarray:
        """Indices of molecules in the neighbouring cells (a superset).

        Unordered: the clash decision is "any overlap", so order is irrelevant.
        """
        if self.covers_all:
            return np.arange(self.count, dtype=np.int64)
        found = []
        cells = self.cells
        for cell in self._neighbour_cells(self._cell_of(centre)):
            members = cells.get(cell)
            if members:
                found.extend(members)
        return np.array(found, dtype=np.int64)

    def _centre_filter(self, centre, candidate_radius, found: np.ndarray) -> np.ndarray:
        centres = self.centres[found]
        dx = centre[0] - centres[:, 0]
        dy = centre[1] - centres[:, 1]
        dz = centre[2] - centres[:, 2]
        if self.periodic:
            a, b, c = self.lengths
            if a > 0.0:
                dx = dx - a * np.round(dx / a)
            if b > 0.0:
                dy = dy - b * np.round(dy / b)
            if c > 0.0:
                dz = dz - c * np.round(dz / c)
        distance = np.sqrt(dx * dx + dy * dy + dz * dz)
        keep = ~(distance > (candidate_radius + self.radii[found] + CENTRE_FILTER_MARGIN))
        return found[keep]

    def nearby(self, centre, candidate_radius: float) -> List[int]:
        """Molecules passing the legacy centre filter, in acceptance order."""
        found = np.sort(self.candidates(centre))
        if len(found) == 0:
            return []
        return [int(i) for i in self._centre_filter(centre, candidate_radius, found)]

    def _atom_indices(self, molecules: np.ndarray) -> np.ndarray:
        counts = self.mol_count[molecules]
        starts = self.mol_start[molecules]
        total = int(counts.sum())
        offsets = np.repeat(starts - (np.cumsum(counts) - counts), counts)
        return offsets + np.arange(total)

    def clashes(self, cand_xyz, cand_r, centre, candidate_radius, vdw_scale,
                extra: float = 0.0, engine: str = ENGINE_NUMPY) -> bool:
        found = self.candidates(centre)
        if len(found) == 0:
            return False
        if engine == ENGINE_NUMBA:
            return bool(solventPackingNumba.grid_clash(
                found, cand_xyz, cand_r, centre, float(candidate_radius), self.centres, self.radii,
                self.mol_start, self.mol_count, self.atoms, self.atom_r, float(vdw_scale),
                self.lengths, self.periodic, float(extra), CENTRE_FILTER_MARGIN,
            ))
        near = self._centre_filter(centre, candidate_radius, found)
        if len(near) == 0:
            return False
        atoms = self._atom_indices(near)
        return atoms_clash_numpy(cand_xyz, cand_r, self.atoms[atoms], self.atom_r[atoms], vdw_scale,
                                 self.lengths, self.periodic, extra)


class RotationStream:
    """Random rotation matrices drawn in blocks from numpy's global generator.

    ``Rotation.random(n)`` consumes the global ``np.random`` state exactly as
    ``n`` calls of ``Rotation.random()`` do, and ``x @ R.T`` reproduces
    ``Rotation.apply`` bit for bit, so drawing in blocks changes no result --
    it only removes ~20 us of scipy call overhead per insertion attempt.
    """

    def __init__(self, block: int = 512):
        self.block = int(block)
        self._matrices: Optional[np.ndarray] = None
        self._next = 0

    def next(self) -> np.ndarray:
        if self._matrices is None or self._next >= len(self._matrices):
            from scipy.spatial.transform import Rotation
            self._matrices = Rotation.random(self.block).as_matrix()
            self._next = 0
        matrix = self._matrices[self._next]
        self._next += 1
        return matrix
