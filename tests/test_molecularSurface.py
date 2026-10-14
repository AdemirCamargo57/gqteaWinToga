"""Tests for isosurface extraction and Multiwfn-style surface extrema."""
import numpy as np
import pytest
from scipy import sparse

import gaussianWavefunction as gw
import molecularSurface as ms
from esp_fixtures import fchk_path

SPACING = 0.25


def gaussian_blobs(centres):
    centres = np.atleast_2d(centres)

    def density(p):
        d2 = ((p[:, None, :] - centres[None, :, :]) ** 2).sum(axis=2)
        return np.exp(-d2).sum(axis=1)
    return density


def mesh_for(density, isovalue, lo, hi):
    origin = np.array(lo, dtype=float)
    shape = tuple(int(round((h - l) / SPACING)) + 1 for l, h in zip(lo, hi))
    values = density(ms.grid_points(origin, SPACING, shape)).reshape(shape)
    return ms.build_isosurface(values, origin, SPACING, isovalue, density, ms.density_outward(density))


@pytest.fixture(scope="module")
def sphere():
    radius = 3.0
    return radius, mesh_for(gaussian_blobs([0, 0, 0]), np.exp(-radius ** 2), (-5, -5, -5), (5, 5, 5))


class TestSphere:
    def test_closed_single_surface(self, sphere):
        _, mesh = sphere
        assert mesh.euler_characteristic() == 2
        assert mesh.component_euler_characteristics() == [2]

    def test_vertices_lie_on_the_sphere(self, sphere):
        radius, mesh = sphere
        assert np.abs(np.linalg.norm(mesh.vertices, axis=1) - radius).max() < SPACING / 10

    def test_area(self, sphere):
        radius, mesh = sphere
        assert mesh.area() == pytest.approx(4 * np.pi * radius ** 2, rel=0.01)

    def test_normals_point_outward(self, sphere):
        _, mesh = sphere
        n = mesh.vertex_normals()
        assert np.all(np.einsum("ij,ij->i", n, mesh.vertices) > 0)

    def test_redundant_vertices_are_removed(self, sphere):
        radius, mesh = sphere
        density = gaussian_blobs([0, 0, 0])
        shape = (41, 41, 41)
        origin = np.array([-5.0, -5.0, -5.0])
        values = density(ms.grid_points(origin, SPACING, shape)).reshape(shape)
        raw, _ = ms.marching_tetrahedra(values, origin, SPACING, np.exp(-radius ** 2), density)
        # Lu & Chen report > 60 % of the marching-tetrahedra vertices removed.
        # (No minimum spacing is guaranteed afterwards: Multiwfn moves the
        # survivor to the midpoint of each merge, which can create new close
        # pairs, and it does not iterate.)
        assert len(mesh.vertices) < 0.5 * len(raw)


class TestTopology:
    def test_two_separate_spheres(self):
        density = gaussian_blobs([[-4, 0, 0], [4, 0, 0]])
        mesh = mesh_for(density, np.exp(-4.0), (-7, -3, -3), (7, 3, 3))
        assert mesh.component_euler_characteristics() == [2, 2]

    def test_no_surface_inside_the_grid(self):
        density = gaussian_blobs([0, 0, 0])
        shape = (9, 9, 9)
        values = density(ms.grid_points(np.zeros(3) - 1, SPACING, shape)).reshape(shape)
        with pytest.raises(ValueError, match="No isosurface"):
            ms.marching_tetrahedra(values, np.zeros(3) - 1, SPACING, 5.0)


class TestMolecularSurface:
    @pytest.fixture(scope="class")
    def eteno(self):
        w = gw.read_fchk(fchk_path("eteno_b3lyp_631g"))

        def density(p):
            return gw.electron_density(w, p)

        radii = np.array([1.7, 1.2, 1.2, 1.7, 1.2, 1.2]) / 0.529177210544
        origin, shape = ms.grid_box(w.coordinates, radii, density, 0.0004, SPACING)
        values = density(ms.grid_points(origin, SPACING, shape)).reshape(shape)
        mesh = ms.build_isosurface(values, origin, SPACING, 0.0004, density, ms.density_outward(density))
        return w, density, origin, shape, mesh

    def test_grid_box_faces_are_below_a_tenth_of_the_isovalue(self, eteno):
        _, density, origin, shape, _ = eteno
        values = density(ms.grid_points(origin, SPACING, shape)).reshape(shape)
        faces = [values[0], values[-1], values[:, 0], values[:, -1], values[:, :, 0], values[:, :, -1]]
        assert max(f.max() for f in faces) < 0.0004 / 10

    def test_closed_and_on_the_isovalue(self, eteno):
        _, density, _, _, mesh = eteno
        assert mesh.component_euler_characteristics() == [2]
        assert np.abs(density(mesh.vertices) - 0.0004).max() < 1e-5


def chain(values):
    """A path graph 0-1-2-...: the smallest mesh that tells 1-ring from 2-ring."""
    n = len(values)
    rows = list(range(n - 1)) + list(range(1, n))
    cols = list(range(1, n)) + list(range(n - 1))
    nbr = sparse.csr_matrix((np.ones(len(rows), bool), (rows, cols)), shape=(n, n))
    return np.asarray(values, float), np.zeros((n, 3)), nbr


class TestExtrema:
    def test_second_ring_is_required(self):
        # vertex 3 (value 2) beats its direct neighbours (3, 6) but not vertex 1
        # (value 1) two steps away, so Multiwfn does not call it a minimum.
        values, verts, nbr = chain([5, 1, 3, 2, 6])
        minima = [e for e in ms.find_surface_extrema(values, verts, nbr) if e.kind == "minimum"]
        assert [(e.vertex, e.scope) for e in minima] == [(1, "global")]

    def test_ties_give_no_local_extremum(self):
        values, verts, nbr = chain([2, 1, 1, 2, 3])
        minima = [e for e in ms.find_surface_extrema(values, verts, nbr) if e.kind == "minimum"]
        # neither 1 is strictly below the other; the global minimum is still reported
        assert len(minima) == 1 and minima[0].scope == "global" and minima[0].value == 1.0

    def test_sphere_poles(self, sphere):
        _, mesh = sphere
        ex = ms.find_surface_extrema(mesh.vertices[:, 2], mesh.vertices, mesh.neighbours)
        assert [(e.kind, e.scope) for e in ex] == [("minimum", "global"), ("maximum", "global")]
        assert ex[0].position[2] < -2.9 and ex[1].position[2] > 2.9

    def test_two_maxima_two_minima(self, sphere):
        _, mesh = sphere
        x, y = mesh.vertices[:, 0], mesh.vertices[:, 1]
        ex = ms.find_surface_extrema(x * x - y * y + 0.01 * x, mesh.vertices, mesh.neighbours)
        kinds = [(e.kind, e.scope) for e in ex]
        assert kinds.count(("maximum", "global")) == 1 and kinds.count(("maximum", "local")) == 1
        assert kinds.count(("minimum", "global")) == 1 and kinds.count(("minimum", "local")) == 1
        # global first, then locals in order of extremity
        maxima = [e for e in ex if e.kind == "maximum"]
        assert maxima[0].value >= maxima[1].value


def test_numerical_modules_do_not_import_toga():
    """The fchk/integral/surface modules stay GUI-free (project convention)."""
    import os
    import subprocess
    import sys
    code = ("import sys, gaussianWavefunction, espIntegrals, molecularSurface; "
            "sys.exit(1 if 'toga' in sys.modules else 0)")
    src = os.path.dirname(os.path.abspath(ms.__file__))
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([src] + sys.path))
    assert subprocess.run([sys.executable, "-c", code], env=env, cwd=src).returncode == 0
