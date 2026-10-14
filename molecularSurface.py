"""Isosurface extraction (marching tetrahedra) and surface extrema.

Follows Multiwfn's quantitative molecular surface analysis (surfana.f90):
exact-function bisection of cut edges, merging of vertices closer than
0.5 x grid spacing, elimination of 2- and 3-connected vertices, and the
strict first+second-neighbour extremum test. Pure numpy/scipy, no Toga.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Set, Tuple

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

DEFAULT_SPACING_BOHR = 0.25
MERGE_RATIO = 0.5
N_BISECT = 3

# Kuhn (Freudenthal) decomposition of a cube into 6 tetrahedra sharing the
# main diagonal 0-7. Corner c has offset (c & 1, (c >> 1) & 1, (c >> 2) & 1).
_CORNER_OFFSETS = np.array([[c & 1, (c >> 1) & 1, (c >> 2) & 1] for c in range(8)])
_TETRAHEDRA = np.array([[0, 7, 1, 3], [0, 7, 3, 2], [0, 7, 2, 6],
                        [0, 7, 6, 4], [0, 7, 4, 5], [0, 7, 5, 1]])
# Every tetrahedron edge joins corners whose offsets differ by one of these.
_EDGE_TYPES = {(1, 0, 0): 0, (0, 1, 0): 1, (0, 0, 1): 2, (1, 1, 0): 3,
               (1, 0, 1): 4, (0, 1, 1): 5, (1, 1, 1): 6}
_EDGE_LUT = np.full(8, -1)
for _off, _k in _EDGE_TYPES.items():
    _EDGE_LUT[_off[0] + 2 * _off[1] + 4 * _off[2]] = _k


def density_outward(density: Callable[[np.ndarray], np.ndarray], eps: float = 1e-3):
    """-grad(rho) by central differences: the outward direction of a rho isosurface."""
    def outward(points: np.ndarray) -> np.ndarray:
        g = np.empty_like(points)
        for d in range(3):
            step = np.zeros(3); step[d] = eps
            g[:, d] = density(points - step) - density(points + step)
        return g
    return outward


@dataclass
class SurfaceMesh:
    vertices: np.ndarray      # (V, 3) bohr
    triangles: np.ndarray     # (F, 3) indices into vertices, outward oriented
    neighbours: sparse.csr_matrix  # (V, V) boolean vertex adjacency
    isovalue: float
    spacing: float

    @property
    def n_edges(self) -> int:
        return int(self.neighbours.nnz // 2)

    def euler_characteristic(self) -> int:
        return len(self.vertices) + len(self.triangles) - self.n_edges

    def component_labels(self) -> np.ndarray:
        """Connected-component label of every vertex (0-based)."""
        return connected_components(self.neighbours, directed=False)[1]

    def component_euler_characteristics(self) -> List[int]:
        """V + F - E of each connected component; 2 for every closed sphere-like surface."""
        labels = self.component_labels()
        n = int(labels.max()) + 1 if labels.size else 0
        v = np.bincount(labels, minlength=n)
        f = np.bincount(labels[self.triangles[:, 0]], minlength=n) if len(self.triangles) else np.zeros(n, int)
        rows, cols = self.neighbours.nonzero()
        e = np.bincount(labels[rows[rows < cols]], minlength=n)
        return [int(x) for x in v + f - e]

    def area(self) -> float:
        a, b, c = (self.vertices[self.triangles[:, k]] for k in range(3))
        return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())

    def vertex_normals(self) -> np.ndarray:
        a, b, c = (self.vertices[self.triangles[:, k]] for k in range(3))
        fn = np.cross(b - a, c - a)
        n = np.zeros_like(self.vertices)
        for k in range(3):
            np.add.at(n, self.triangles[:, k], fn)
        length = np.linalg.norm(n, axis=1)
        return n / np.where(length > 0, length, 1.0)[:, None]


def grid_points(origin: np.ndarray, spacing: float, shape: Tuple[int, int, int]) -> np.ndarray:
    axes = [origin[d] + spacing * np.arange(shape[d]) for d in range(3)]
    gx, gy, gz = np.meshgrid(*axes, indexing="ij")
    return np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])


def grid_box(coordinates: np.ndarray, radii: np.ndarray, density: Callable[[np.ndarray], np.ndarray],
             isovalue: float, spacing: float, max_rounds: int = 30):
    """Origin and shape of a cubic grid whose every face has rho < isovalue / 10.

    Starts from atoms +/- 1.7 x van der Waals radius (Multiwfn's vdwmulti) and
    pushes any face whose maximum density is still too high out by 1 bohr.
    """
    lo = (coordinates - 1.7 * radii[:, None]).min(axis=0)
    hi = (coordinates + 1.7 * radii[:, None]).max(axis=0)
    limit = isovalue / 10.0
    for _ in range(max_rounds):
        shape = tuple(int(np.ceil((hi[d] - lo[d]) / spacing)) + 1 for d in range(3))
        grown = False
        for d in range(3):
            for side in (0, 1):
                face_shape = list(shape)
                face_shape[d] = 1
                face_origin = lo.copy()
                if side:
                    face_origin[d] = lo[d] + (shape[d] - 1) * spacing
                pts = grid_points(face_origin, spacing, tuple(face_shape))
                if density(pts).max() >= limit:
                    grown = True
                    if side:
                        hi[d] += 1.0
                    else:
                        lo[d] -= 1.0
        if not grown:
            return lo, shape
    raise ValueError("Could not enclose the isosurface: density stays above "
                     f"{limit:.3g} a.u. on the grid boundary.")


def marching_tetrahedra(values: np.ndarray, origin: np.ndarray, spacing: float, isovalue: float,
                        density: Callable[[np.ndarray], np.ndarray] | None = None,
                        n_bisect: int = N_BISECT):
    """Raw triangulation of {rho = isovalue}: (vertices, triangles).

    Vertices are shared between triangles (one per cut grid edge), so the
    result is a connected mesh. With ``density`` the cut point on each edge
    is refined by ``n_bisect`` bisections of the exact function before the
    final linear interpolation (Multiwfn nbisec=3).
    """
    nx, ny, nz = values.shape
    inside = values >= isovalue
    cube = np.stack([inside[ox:nx - 1 + ox, oy:ny - 1 + oy, oz:nz - 1 + oz]
                     for ox, oy, oz in _CORNER_OFFSETS], axis=-1)
    boundary = np.argwhere(cube.any(axis=-1) & ~cube.all(axis=-1))
    if boundary.size == 0:
        raise ValueError(f"No isosurface at rho = {isovalue:g} a.u. inside the grid.")
    corner_in = cube[boundary[:, 0], boundary[:, 1], boundary[:, 2]]  # (C, 8)
    corner_idx = boundary[:, None, :] + _CORNER_OFFSETS[None, :, :]   # (C, 8, 3)

    def edge_id(ga, gb):
        lo_ = np.minimum(ga, gb)
        k = _EDGE_LUT[np.abs(gb - ga) @ np.array([1, 2, 4])]
        lin = (lo_[..., 0] * ny + lo_[..., 1]) * nz + lo_[..., 2]
        return lin * 7 + k

    tri_edges, tri_dirs = [], []
    for tet in _TETRAHEDRA:
        ins = corner_in[:, tet]                       # (C, 4)
        code = ins.astype(int) @ np.array([1, 2, 4, 8])
        g = corner_idx[:, tet, :]                     # (C, 4, 3)
        n_in = ins.sum(axis=1)
        for count in (1, 2, 3):
            sel = np.nonzero(n_in == count)[0]
            if sel.size == 0:
                continue
            for c in np.unique(code[sel]):
                rows = sel[code[sel] == c]
                bits = [(c >> k) & 1 for k in range(4)]
                ins_k = [k for k in range(4) if bits[k]]
                out_k = [k for k in range(4) if not bits[k]]
                pairs = [(i, o) for i in ins_k for o in out_k]
                if len(pairs) == 3:
                    polys = [pairs]
                else:  # 2 in, 2 out: quad (i0o0, i0o1, i1o1, i1o0) -> 2 triangles
                    (i0, i1), (o0, o1) = ins_k, out_k
                    quad = [(i0, o0), (i0, o1), (i1, o1), (i1, o0)]
                    polys = [[quad[0], quad[1], quad[2]], [quad[0], quad[2], quad[3]]]
                for poly in polys:
                    ids = np.stack([edge_id(g[rows, i], g[rows, o]) for i, o in poly], axis=1)
                    dirs = np.stack([g[rows, o] - g[rows, i] for i, o in poly], axis=1)
                    tri_edges.append(ids)
                    tri_dirs.append(dirs.sum(axis=1))
    tri_edges = np.concatenate(tri_edges)
    outward = np.concatenate(tri_dirs).astype(float)
    unique_ids, triangles = np.unique(tri_edges, return_inverse=True)
    triangles = triangles.reshape(-1, 3)

    # Decode each unique edge back to its two grid points (a inside, b outside).
    kinds = np.array(list(_EDGE_TYPES.keys()))
    lin, k = np.divmod(unique_ids, 7)
    ga = np.column_stack([lin // (ny * nz), (lin // nz) % ny, lin % nz])
    gb = ga + kinds[k]
    va = values[ga[:, 0], ga[:, 1], ga[:, 2]]
    vb = values[gb[:, 0], gb[:, 1], gb[:, 2]]
    pa = origin + spacing * ga
    pb = origin + spacing * gb
    if density is not None:
        for _ in range(n_bisect):
            mid = 0.5 * (pa + pb)
            vm = density(mid)
            same_as_a = (vm - isovalue) * (va - isovalue) > 0
            pa = np.where(same_as_a[:, None], mid, pa)
            va = np.where(same_as_a, vm, va)
            pb = np.where(same_as_a[:, None], pb, mid)
            vb = np.where(same_as_a, vb, vm)
    scl = (va - isovalue) / (va - vb)
    vertices = pa + scl[:, None] * (pb - pa)

    a, b, c = (vertices[triangles[:, k]] for k in range(3))
    flip = np.einsum("ij,ij->i", np.cross(b - a, c - a), outward) < 0
    triangles[flip] = triangles[flip][:, [0, 2, 1]]
    return vertices, triangles


def _adjacency_sets(n: int, triangles: np.ndarray) -> List[Set[int]]:
    adj: List[Set[int]] = [set() for _ in range(n)]
    for i, j, k in triangles.tolist():
        adj[i].update((j, k)); adj[j].update((i, k)); adj[k].update((i, j))
    return adj


def eliminate_redundant_vertices(vertices: np.ndarray, triangles: np.ndarray, merge_distance: float):
    """Multiwfn's vertex merge + 2-/3-neighbour elimination (surfana.f90:654-780).

    Returns the surviving vertices, re-indexed triangles and adjacency sets.
    """
    v = vertices.copy()
    n = len(v)
    adj = _adjacency_sets(n, triangles)
    dead = np.zeros(n, dtype=bool)
    parent = np.arange(n)
    crit2 = merge_distance ** 2
    for i in range(n):
        if dead[i]:
            continue
        for j in list(adj[i]):
            if dead[j] or j == i:
                continue
            d = v[i] - v[j]
            if d @ d < crit2:
                for k in adj[j]:
                    if dead[k] or k == i:
                        continue
                    if i not in adj[k]:
                        adj[k].add(i)
                        adj[i].add(k)
                parent[j] = i
                dead[j] = True
                v[i] = 0.5 * (v[i] + v[j])

    def root(x):
        while dead[x] and parent[x] != x:
            x = parent[x]
        return x

    tris = [[root(x) for x in t] for t in triangles.tolist()]
    tris = [t for t in tris if len(set(t)) == 3]

    def live(i):
        return [j for j in adj[i] if not dead[j]]

    for i in range(n):  # vertices with exactly two live neighbours
        if not dead[i] and len(live(i)) == 2:
            dead[i] = True
    for i in range(n):  # vertices with exactly three: drop and re-cap
        if not dead[i]:
            nb = live(i)
            if len(nb) == 3:
                dead[i] = True
                tris.append(nb)
    tris = np.array([t for t in tris if not any(dead[x] for x in t)], dtype=np.int64).reshape(-1, 3)
    keep = np.nonzero(~dead)[0]
    new_index = -np.ones(n, dtype=np.int64)
    new_index[keep] = np.arange(keep.size)
    rows, cols = [], []
    for i in keep:
        for j in adj[i]:
            if not dead[j]:
                rows.append(new_index[i]); cols.append(new_index[j])
    m = keep.size
    nbr = sparse.csr_matrix((np.ones(len(rows), dtype=bool), (rows, cols)), shape=(m, m))
    return v[keep], new_index[tris], nbr


def orient_triangles(vertices, triangles, density_gradient_sign):
    """Make every triangle's normal point along ``density_gradient_sign`` (outward)."""
    a, b, c = (vertices[triangles[:, k]] for k in range(3))
    centre = (a + b + c) / 3.0
    normal = np.cross(b - a, c - a)
    flip = np.einsum("ij,ij->i", normal, density_gradient_sign(centre)) < 0
    out = triangles.copy()
    out[flip] = out[flip][:, [0, 2, 1]]
    return out


def build_isosurface(values, origin, spacing, isovalue, density, outward_direction,
                     merge_ratio: float = MERGE_RATIO, n_bisect: int = N_BISECT) -> SurfaceMesh:
    verts, tris = marching_tetrahedra(values, origin, spacing, isovalue, density, n_bisect)
    verts, tris, nbr = eliminate_redundant_vertices(verts, tris, merge_ratio * spacing)
    tris = orient_triangles(verts, tris, outward_direction)
    return SurfaceMesh(verts, tris, nbr, isovalue, spacing)


# ----------------------------------------------------------------------
# Extrema
# ----------------------------------------------------------------------
@dataclass
class SurfaceExtremum:
    vertex: int
    kind: str      # "minimum" | "maximum"
    scope: str     # "global" | "local"
    value: float
    position: np.ndarray  # bohr


def second_ring(neighbours: sparse.csr_matrix) -> sparse.csr_matrix:
    """First + second neighbours of every vertex, the vertex itself excluded."""
    a = neighbours.astype(np.int32)
    ring = ((a + a @ a) > 0).tolil()
    ring.setdiag(False)
    return ring.tocsr()


def find_surface_extrema(values: np.ndarray, vertices: np.ndarray,
                         neighbours: sparse.csr_matrix) -> List[SurfaceExtremum]:
    """Multiwfn's rule: a minimum is strictly below every 1st- and 2nd-ring
    neighbour (maximum: strictly above). Ties give no extremum. The global
    minimum/maximum are argmin/argmax over all vertices and are reported as
    "global" even when a tie keeps them out of the strict local test."""
    ring = second_ring(neighbours)
    has = np.diff(ring.indptr) > 0
    vals = values[ring.indices]
    starts = ring.indptr[:-1][has]
    lo = np.full(values.size, np.inf)
    hi = np.full(values.size, -np.inf)
    lo[has] = np.minimum.reduceat(vals, starts)
    hi[has] = np.maximum.reduceat(vals, starts)
    gmin, gmax = int(np.argmin(values)), int(np.argmax(values))
    out: List[SurfaceExtremum] = []
    for kind, mask, g in (("minimum", has & (values < lo), gmin), ("maximum", has & (values > hi), gmax)):
        idx = set(np.nonzero(mask)[0].tolist()) | {g}
        order = sorted(idx, key=lambda i: values[i], reverse=(kind == "maximum"))
        for i in order:
            out.append(SurfaceExtremum(i, kind, "global" if i == g else "local",
                                       float(values[i]), vertices[i].copy()))
    return out
