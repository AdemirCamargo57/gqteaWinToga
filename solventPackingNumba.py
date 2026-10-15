"""Optional numba kernel for :mod:`solventPacking` -- the only module that imports numba.

Same arithmetic, in the same order, as :func:`solventPacking.atoms_clash_numpy`
(``np.rint`` rounds half to even like ``round``; numba does not contract
multiply-adds without ``fastmath``), so decisions are identical; it returns at
the first overlapping pair instead of building the full distance matrix.
"""
import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def _clash_kernel(cand_xyz, cand_r, other_xyz, other_r, vdw_scale, lengths, periodic, extra):
    a = lengths[0]
    b = lengths[1]
    c = lengths[2]
    for i in range(cand_xyz.shape[0]):
        ax = cand_xyz[i, 0]
        ay = cand_xyz[i, 1]
        az = cand_xyz[i, 2]
        ri = cand_r[i]
        for j in range(other_xyz.shape[0]):
            dx = ax - other_xyz[j, 0]
            dy = ay - other_xyz[j, 1]
            dz = az - other_xyz[j, 2]
            if periodic:
                if a > 0.0:
                    dx = dx - a * np.rint(dx / a)
                if b > 0.0:
                    dy = dy - b * np.rint(dy / b)
                if c > 0.0:
                    dz = dz - c * np.rint(dz / c)
            dist2 = dx * dx + dy * dy + dz * dz
            cutoff = vdw_scale * (ri + other_r[j])
            if dist2 < cutoff * cutoff:
                return True
            if extra > 0.0:
                expanded = cutoff + extra
                if dist2 < expanded * expanded:
                    return True
    return False


@njit(cache=True, nogil=True)
def _grid_clash_kernel(found, cand_xyz, cand_r, centre, cand_radius, centres, radii,
                       mol_start, mol_count, atoms, atom_r, vdw_scale, lengths, periodic, extra, margin):
    """Legacy centre pre-filter + atom test over the candidate molecules, in one call."""
    a = lengths[0]
    b = lengths[1]
    c = lengths[2]
    for f in range(found.shape[0]):
        k = found[f]
        dx = centre[0] - centres[k, 0]
        dy = centre[1] - centres[k, 1]
        dz = centre[2] - centres[k, 2]
        if periodic:
            if a > 0.0:
                dx = dx - a * np.rint(dx / a)
            if b > 0.0:
                dy = dy - b * np.rint(dy / b)
            if c > 0.0:
                dz = dz - c * np.rint(dz / c)
        distance = np.sqrt(dx * dx + dy * dy + dz * dz)
        if distance > (cand_radius + radii[k] + margin):
            continue
        start = mol_start[k]
        stop = start + mol_count[k]
        if _clash_kernel(cand_xyz, cand_r, atoms[start:stop], atom_r[start:stop],
                         vdw_scale, lengths, periodic, extra):
            return True
    return False


def grid_clash(found, cand_xyz, cand_r, centre, cand_radius, centres, radii, mol_start, mol_count,
               atoms, atom_r, vdw_scale, lengths, periodic, extra, margin):
    return _grid_clash_kernel(
        found, np.ascontiguousarray(cand_xyz, dtype=np.float64), cand_r,
        np.ascontiguousarray(centre, dtype=np.float64), cand_radius, centres, radii,
        mol_start, mol_count, atoms, atom_r, vdw_scale, lengths, periodic, extra, margin,
    )


def atoms_clash(cand_xyz, cand_r, other_xyz, other_r, vdw_scale, lengths, periodic, extra):
    return _clash_kernel(
        np.ascontiguousarray(cand_xyz, dtype=np.float64),
        np.ascontiguousarray(cand_r, dtype=np.float64),
        np.ascontiguousarray(other_xyz, dtype=np.float64),
        np.ascontiguousarray(other_r, dtype=np.float64),
        float(vdw_scale), np.ascontiguousarray(lengths, dtype=np.float64), bool(periodic), float(extra),
    )
