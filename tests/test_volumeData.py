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


def small(values=None, origin=(0.0, 0.0, 0.0), axes=np.eye(3) * 0.5, shape=(3, 4, 5)):
    values = np.arange(np.prod(shape), dtype=float).reshape(shape) if values is None else values
    return cf.CubeVolume("t", "c", np.array([1]), np.array([1.0]), np.zeros((1, 3)),
                         np.asarray(origin, float), np.asarray(axes, float), values)


class TestGridChecks:
    def test_identical_grids_match(self):
        assert vd.grid_mismatch(small(), small()) is None

    def test_differences_within_tolerance_match(self):
        assert vd.grid_mismatch(small(), small(origin=(5e-5, 0, 0))) is None

    def test_shape_mismatch_names_both_sizes(self):
        message = vd.grid_mismatch(small(), small(shape=(3, 4, 6)))
        assert "3 x 4 x 5" in message and "3 x 4 x 6" in message

    def test_origin_mismatch(self):
        assert "origin" in vd.grid_mismatch(small(), small(origin=(0.01, 0, 0)))

    def test_axes_mismatch(self):
        assert "axes" in vd.grid_mismatch(small(), small(axes=np.eye(3) * 0.51))

    def test_atoms_mismatch(self):
        a = small()
        b = small()
        assert vd.atoms_mismatch(a, b) is None
        b.coordinates = np.array([[0.1, 0.0, 0.0]])
        assert "atoms" in vd.atoms_mismatch(a, b)
        b.atomic_numbers = np.array([1, 1])
        assert "atoms" in vd.atoms_mismatch(a, b)


class TestOperations:
    def test_sum_and_difference(self):
        a, b = small(), small(values=np.ones((3, 4, 5)))
        assert np.array_equal(vd.combine(a, b, "add").values, a.values + 1)
        assert np.array_equal(vd.combine(a, b, "subtract").values, a.values - 1)

    def test_a_minus_a_is_zero_and_keeps_the_grid_of_a(self):
        a = small(origin=(1.0, 2.0, 3.0))
        result = vd.combine(a, a, "subtract", title="A - A")
        assert not result.values.any() and result.title == "A - A"
        assert np.array_equal(result.origin, a.origin) and np.array_equal(result.axes, a.axes)
        assert result.orbital is None and result.values is not a.values

    def test_scale_is_linear(self):
        a = small()
        assert np.array_equal(vd.scale(a, 2.0).values, 2.0 * a.values)
        assert np.array_equal(vd.scale(a, -1.0).values, -a.values)

    def test_incompatible_grids_are_refused_with_the_reason(self):
        with pytest.raises(ValueError, match="3 x 4 x 6"):
            vd.combine(small(), small(shape=(3, 4, 6)), "add")

    def test_unknown_operation_and_bad_factor_are_refused(self):
        with pytest.raises(ValueError, match="operation"):
            vd.combine(small(), small(), "multiply")
        for factor in (float("nan"), float("inf")):
            with pytest.raises(ValueError, match="factor"):
                vd.scale(small(), factor)

    def test_operations_do_not_modify_their_inputs(self):
        a = small()
        before = a.values.copy()
        vd.combine(a, a, "add")
        vd.scale(a, 3.0)
        assert np.array_equal(a.values, before)


class TestStatistics:
    def test_voxel_volume_of_a_skewed_cell(self):
        axes = np.array([[0.2, 0.0, 0.0], [0.1, 0.3, 0.0], [0.0, 0.05, 0.4]])
        assert vd.voxel_volume(small(axes=axes)) == pytest.approx(0.2 * 0.3 * 0.4)

    def test_integral_of_a_normalised_gaussian_on_a_skewed_grid(self):
        def normalised(points):
            return np.exp(-(points ** 2).sum(axis=1)) / np.pi ** 1.5
        h = 0.15
        axes = np.array([[h, 0.0, 0.0], [h / 2, h, 0.0], [0.0, 0.0, h]])
        volume = volume_from(normalised, (90, 60, 60), origin=(-8.0, -4.5, -4.5), axes=axes)
        assert vd.integral(volume) == pytest.approx(1.0, abs=1e-3)

    def test_value_statistics(self):
        lo, hi, total = vd.value_statistics(small())
        assert (lo, hi) == (0.0, 59.0) and total == pytest.approx(59 * 60 / 2 * 0.125)

    def test_default_isovalues(self):
        density = small(values=np.linspace(0.0, 2.0, 60).reshape(3, 4, 5))
        orbital = small(values=np.linspace(-0.3, 0.3, 60).reshape(3, 4, 5))
        orbital.orbital = 7
        signed = small(values=np.linspace(-0.3, 0.5, 60).reshape(3, 4, 5))
        tiny = small(values=np.linspace(-0.004, 0.004, 60).reshape(3, 4, 5))
        assert vd.default_isovalue(density) == 0.0004
        assert vd.default_isovalue(orbital) == 0.02
        assert vd.default_isovalue(signed) == 0.02
        assert vd.default_isovalue(tiny) == 0.002
        assert vd.default_isovalue(small(values=np.zeros((3, 4, 5)))) == 0.02


class TestSlices:
    @pytest.mark.parametrize("axis", [0, 1, 2])
    def test_plane_values_and_points(self, axis):
        axes = np.array([[0.2, 0.0, 0.0], [0.1, 0.3, 0.0], [0.0, 0.05, 0.4]])
        volume = small(axes=axes, origin=(1.0, -1.0, 0.5))
        points, values = vd.slice_plane(volume, axis, 2)
        assert np.array_equal(values, np.take(volume.values, 2, axis=axis))
        assert points.shape == values.shape + (3,)
        index = [0, 0, 0]
        index[axis] = 2
        assert np.allclose(points[0, 0], volume.origin + np.array(index) @ axes)
        others = [d for d in range(3) if d != axis]
        assert np.allclose(points[1, 0] - points[0, 0], axes[others[0]])
        assert np.allclose(points[0, 1] - points[0, 0], axes[others[1]])

    def test_out_of_range_plane_and_axis_are_refused(self):
        with pytest.raises(ValueError, match="plane"):
            vd.slice_plane(small(), 0, 3)
        with pytest.raises(ValueError, match="axis"):
            vd.slice_plane(small(), 3, 0)

    def test_slice_triangles_cover_the_plane(self):
        triangles = vd.slice_triangles(3, 4)
        assert triangles.shape == (2 * 2 * 3, 3)
        assert set(triangles.ravel()) == set(range(12))

    def test_slice_table_is_in_angstrom_with_one_row_per_point(self):
        volume = small(origin=(1.0, 0.0, 0.0))
        text = vd.slice_table_text(volume, 2, 0, name="rho")
        rows = [line for line in text.splitlines() if not line.startswith("#")]
        assert len(rows) == 3 * 4
        first = [float(x) for x in rows[0].split()]
        assert first[:3] == pytest.approx([vd.BOHR_RADIUS_ANGSTROM, 0.0, 0.0], abs=1e-6)
        assert first[3] == 0.0
        assert "rho" in text and "angstrom" in text

    def test_colour_range_modes(self):
        values = np.array([-1.0, 0.5, 3.0])
        assert vd.colour_range(values, "symmetric") == (-3.0, 3.0)
        assert vd.colour_range(values, "minmax") == (-1.0, 3.0)
        assert vd.colour_range(np.zeros(3), "symmetric") == (-1.0, 1.0)
        with pytest.raises(ValueError):
            vd.colour_range(values, "other")


class TestColours:
    def test_presets_and_hex(self):
        assert vd.parse_colour("blue") == vd.COLOUR_PRESETS["Blue"]
        assert vd.parse_colour(" #FF8000 ") == pytest.approx((1.0, 128 / 255, 0.0))

    @pytest.mark.parametrize("text", ["", "teal-ish", "#12345", "#GG0000"])
    def test_bad_colours_are_refused(self, text):
        with pytest.raises(ValueError, match="colour"):
            vd.parse_colour(text)

    def test_diverging_colours(self):
        rgb = vd.diverging_colours(np.array([0.0, 2.0, -2.0, 9.0]), 2.0)
        assert rgb[0] == pytest.approx([1, 1, 1]) and rgb[1] == pytest.approx([1, 0, 0])
        assert rgb[2] == pytest.approx([0, 0, 1]) and rgb[3] == pytest.approx([1, 0, 0])

    def test_rainbow_colours_run_red_green_blue(self):
        rgb = vd.rainbow_colours(np.array([0.0, 0.5, 1.0]), 0.0, 1.0)
        assert rgb[0] == pytest.approx([1, 0, 0]) and rgb[1] == pytest.approx([0, 1, 0])
        assert rgb[2] == pytest.approx([0, 0, 1])
        assert vd.rainbow_colours(np.array([4.0]), 4.0, 4.0)[0] == pytest.approx([0, 1, 0])


class TestStructure:
    def test_frame_in_angstrom_with_symbols(self):
        volume = small()
        volume.atomic_numbers = np.array([8, 1, 0])
        volume.coordinates = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        frame = vd.structure_frame(volume)
        assert [s for s, _ in frame] == ["O", "H", "X"]
        assert frame[1][1] == pytest.approx((vd.BOHR_RADIUS_ANGSTROM, 0.0, 0.0))
