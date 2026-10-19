"""Tests for volumeData (the Cube tab's numerics) and its numba kernel."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import cubeFile as cf
import volumeData as vd

SRC = str(Path(vd.__file__).parent)


def volume_from(function, shape, origin=(0.0, 0.0, 0.0), axes=np.eye(3) * 0.2, orbital=None,
                atoms=((8, (0.0, 0.0, 0.0)),)):
    """A CubeVolume sampling ``function(points (n, 3) bohr)`` on the given lattice."""
    origin = np.asarray(origin, dtype=float)
    axes = np.asarray(axes, dtype=float)
    idx = np.stack(np.meshgrid(*[np.arange(n) for n in shape], indexing="ij"), axis=-1).reshape(-1, 3)
    points = origin + idx @ axes
    z = np.array([a[0] for a in atoms], dtype=int)
    xyz = np.array([a[1] for a in atoms], dtype=float).reshape(-1, 3)
    return cf.CubeVolume("t", "c", z, z.astype(float), xyz, origin, axes,
                         function(points).reshape(shape), orbital)


def blob(points):
    return np.exp(-(points ** 2).sum(axis=1) / 2.0)


def p_orbital(points):
    return points[:, 0] * np.exp(-(points ** 2).sum(axis=1))


BLOB_RADIUS = np.sqrt(2.0 * np.log(2.0))  # where blob = 0.5


@pytest.fixture(scope="module")
def blob_volume():
    return volume_from(blob, (41, 41, 41), origin=(-4.0, -4.0, -4.0))


class TestIsosurface:
    def test_blob_surface_is_a_sphere(self, blob_volume):
        mesh = vd.isosurface(blob_volume, 0.5, use_numba=False)
        radii = np.linalg.norm(mesh.vertices, axis=1)
        assert np.abs(radii - BLOB_RADIUS).max() < 0.2
        assert mesh.triangles.shape[1] == 3 and mesh.normals.shape == mesh.vertices.shape

    def test_normals_point_away_from_higher_values(self, blob_volume):
        mesh = vd.isosurface(blob_volume, 0.5, use_numba=False)
        assert (np.einsum("ij,ij->i", mesh.normals, mesh.vertices) > 0).mean() > 0.99

    def test_positive_and_negative_lobes_mirror_each_other(self):
        volume = volume_from(p_orbital, (41, 41, 41), origin=(-4.0, -4.0, -4.0))
        plus = vd.isosurface(volume, 0.1, sign=1, use_numba=False)
        minus = vd.isosurface(volume, 0.1, sign=-1, use_numba=False)
        assert plus.vertices[:, 0].min() > 0 and minus.vertices[:, 0].max() < 0
        assert np.sort(plus.vertices[:, 0]) == pytest.approx(np.sort(-minus.vertices[:, 0]), abs=1e-9)

    def test_skewed_grid_gives_the_same_sphere(self):
        h = 0.2
        axes = np.array([[h, 0.0, 0.0], [h / 2, h, 0.0], [0.0, h / 3, h]])
        volume = volume_from(blob, (60, 45, 45), origin=(-6.5, -4.5, -4.5), axes=axes)
        mesh = vd.isosurface(volume, 0.5, use_numba=False)
        radii = np.linalg.norm(mesh.vertices, axis=1)
        assert np.abs(radii - BLOB_RADIUS).max() < 0.2

    def test_left_handed_axes_keep_normals_outward(self):
        axes = np.diag([0.2, 0.2, -0.2])
        volume = volume_from(blob, (41, 41, 41), origin=(-4.0, -4.0, 4.0), axes=axes)
        mesh = vd.isosurface(volume, 0.5, use_numba=False)
        a, b, c = (mesh.vertices[mesh.triangles[:, k]] for k in range(3))
        face = np.cross(b - a, c - a)
        centre = (a + b + c) / 3
        assert (np.einsum("ij,ij->i", face, centre) > 0).mean() > 0.99

    def test_isovalue_outside_the_data_gives_none(self, blob_volume):
        assert vd.isosurface(blob_volume, 2.0, use_numba=False) is None
        assert vd.isosurface(blob_volume, 0.5, sign=-1, use_numba=False) is None

    def test_bad_arguments_are_refused(self, blob_volume):
        with pytest.raises(ValueError, match="isovalue"):
            vd.isosurface(blob_volume, float("nan"))
        with pytest.raises(ValueError, match="sign"):
            vd.isosurface(blob_volume, 0.5, sign=2)

    @pytest.mark.skipif(not vd.numba_available(), reason="numba not installed")
    def test_numba_and_numpy_give_the_same_mesh(self):
        volume = volume_from(p_orbital, (33, 30, 28), origin=(-3.2, -2.9, -2.7))
        for sign in (1, -1):
            a = vd.isosurface(volume, 0.08, sign=sign, use_numba=False)
            b = vd.isosurface(volume, 0.08, sign=sign, use_numba=True)
            assert np.array_equal(a.vertices, b.vertices)
            key = lambda t: sorted(map(tuple, t.tolist()))  # noqa: E731
            assert key(a.triangles) == key(b.triangles)

    @pytest.mark.skipif(not vd.numba_available(), reason="numba not installed")
    def test_numba_path_returns_none_when_nothing_is_cut(self, blob_volume):
        assert vd.isosurface(blob_volume, 2.0, use_numba=True) is None

    def test_forcing_numba_without_it_is_a_value_error(self, blob_volume, monkeypatch):
        monkeypatch.setattr(vd, "volumeNumba", None)
        assert not vd.numba_available()
        with pytest.raises(ValueError, match="numba"):
            vd.isosurface(blob_volume, 0.5, use_numba=True)
        assert vd.isosurface(blob_volume, 0.5) is not None  # None = "if installed" -> numpy


class TestModuleImports:
    @pytest.mark.parametrize("module", ["volumeData", "volumeNumba"])
    def test_does_not_import_toga(self, module):
        if module == "volumeNumba" and not vd.numba_available():
            pytest.skip("numba not installed")
        code = f"import sys, {module}; print('toga' in sys.modules)"
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=SRC, check=True)
        assert out.stdout.strip() == "False"
