"""Tests for the Gaussian cube writer/reader."""
import subprocess
import sys

import numpy as np
import pytest

import cubeFile as cf


def sample_cube(shape=(3, 4, 7)):
    rng = np.random.default_rng(7)
    return cf.CubeData("title line", "comment line", np.array([8, 1]), np.array([8.0, 1.0]),
                       np.array([[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]]), np.array([-4.0, -3.5, -2.0]),
                       0.25, rng.normal(scale=0.01, size=shape))


class TestCube:
    def test_round_trip_keeps_five_significant_digits(self, tmp_path):
        cube = sample_cube()
        back = cf.read_cube(cf.write_cube(tmp_path / "a.cub", cube))
        assert back.values.shape == cube.values.shape and back.spacing == 0.25
        assert np.allclose(back.origin, cube.origin) and np.array_equal(back.atomic_numbers, [8, 1])
        assert np.allclose(back.values, cube.values, rtol=1e-5, atol=1e-12)

    def test_layout_is_six_per_line_and_restarts_each_column(self):
        text = cf.format_cube(sample_cube((1, 2, 7)))
        data = text.splitlines()[6 + 2:]
        assert [len(line.split()) for line in data] == [6, 1, 6, 1]
        assert text.splitlines()[3].split() == ["1", "0.250000", "0.000000", "0.000000"]

    def test_truncated_file_is_a_value_error(self, tmp_path):
        path = tmp_path / "t.cub"
        path.write_text("\n".join(cf.format_cube(sample_cube()).splitlines()[:-2]) + "\n")
        with pytest.raises(ValueError, match="truncated"):
            cf.read_cube(path)

    def test_skewed_grid_is_rejected(self, tmp_path):
        lines = cf.format_cube(sample_cube()).splitlines()
        lines[4] = "    4    0.100000    0.250000    0.000000"
        path = tmp_path / "s.cub"
        path.write_text("\n".join(lines) + "\n")
        with pytest.raises(ValueError, match="axis-aligned"):
            cf.read_cube(path)

    def test_unwritable_path_is_a_value_error(self, tmp_path):
        (tmp_path / "file").write_text("x")
        with pytest.raises(ValueError, match="Could not write"):
            cf.write_cube(tmp_path / "file" / "a.cub", sample_cube())

    def test_module_does_not_import_toga(self):
        code = "import sys, cubeFile; print('toga' in sys.modules)"
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                             cwd=str(__import__('pathlib').Path(cf.__file__).parent), check=True)
        assert out.stdout.strip() == "False"


# ----------------------------------------------------------------------
# General reader/writer (CubeVolume): orbital cubes, Angstrom headers,
# skewed and unequal-spacing grids. read_cube above stays strict.
# ----------------------------------------------------------------------
def cube_text(values, origin=(0.0, 0.0, 0.0), axes=np.eye(3) * 0.5, counts_sign=1,
              atoms=((8, 8.0, (0.0, 0.0, 0.0)),), mos=None, per_line=6):
    """Hand-written cube text; values (nx, ny, nz) or (nx, ny, nz, nmo) for an orbital cube."""
    values = np.asarray(values, dtype=float)
    nx, ny, nz = values.shape[:3]
    natoms = -len(atoms) if mos is not None else len(atoms)
    lines = ["title", "comment", f"{natoms:5d}" + "".join(f"{c:12.6f}" for c in origin)]
    for n, axis in zip((nx, ny, nz), np.asarray(axes, dtype=float)):
        lines.append(f"{counts_sign * n:5d}" + "".join(f"{c:12.6f}" for c in axis))
    for z, q, xyz in atoms:
        lines.append(f"{z:5d}{q:12.6f}" + "".join(f"{c:12.6f}" for c in xyz))
    if mos is not None:
        ids = [len(mos)] + list(mos)
        for k in range(0, len(ids), 10):
            lines.append("".join(f"{i:5d}" for i in ids[k:k + 10]))
    columns = values.reshape(nx * ny, -1)
    for column in columns:
        for k in range(0, len(column), per_line):
            lines.append("".join(f"{v:13.5E}" for v in column[k:k + per_line]))
    return "\n".join(lines) + "\n"


def write_text(tmp_path, text, name="v.cube"):
    path = tmp_path / name
    path.write_text(text)
    return path


class TestCubeVolume:
    def test_plain_cube_matches_strict_reader(self, tmp_path):
        path = cf.write_cube(tmp_path / "a.cub", sample_cube())
        (volume,) = cf.read_cube_volumes(path)
        strict = cf.read_cube(path)
        assert np.array_equal(volume.values, strict.values)
        assert np.allclose(volume.axes, np.eye(3) * 0.25) and np.allclose(volume.origin, strict.origin)
        assert np.array_equal(volume.atomic_numbers, [8, 1]) and volume.orbital is None
        assert np.allclose(volume.coordinates, strict.coordinates)

    def test_orbital_cube_with_one_mo(self, tmp_path):
        values = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4, 1) * 1e-3
        (volume,) = cf.read_cube_volumes(write_text(tmp_path, cube_text(values, mos=[9])))
        assert volume.orbital == 9
        assert np.allclose(volume.values, values[..., 0])

    def test_orbital_cube_with_three_mos_is_split(self, tmp_path):
        rng = np.random.default_rng(3)
        values = rng.normal(size=(2, 3, 5, 3))
        volumes = cf.read_cube_volumes(write_text(tmp_path, cube_text(values, mos=[5, 6, 7])))
        assert [v.orbital for v in volumes] == [5, 6, 7]
        for k, v in enumerate(volumes):
            assert np.allclose(v.values, values[..., k], rtol=1e-5)

    def test_mo_list_may_wrap_over_two_lines(self, tmp_path):
        values = np.ones((1, 1, 2, 12)) * np.arange(12)
        volumes = cf.read_cube_volumes(write_text(tmp_path, cube_text(values, mos=list(range(1, 13)))))
        assert [v.orbital for v in volumes] == list(range(1, 13))
        assert np.allclose(volumes[11].values, 11.0)

    def test_negative_voxel_counts_mean_angstrom(self, tmp_path):
        a = cf.BOHR_RADIUS_ANGSTROM
        text = cube_text(np.zeros((2, 2, 2)), origin=(a, 0.0, 0.0), axes=np.eye(3) * a,
                         counts_sign=-1, atoms=((1, 1.0, (0.0, 2 * a, 0.0)),))
        (volume,) = cf.read_cube_volumes(write_text(tmp_path, text))
        assert volume.values.shape == (2, 2, 2)
        assert np.allclose(volume.axes, np.eye(3), atol=1e-5)
        assert np.allclose(volume.origin, [1.0, 0.0, 0.0], atol=1e-5)
        assert np.allclose(volume.coordinates, [[0.0, 2.0, 0.0]], atol=1e-5)

    def test_skewed_and_unequal_axes_are_kept(self, tmp_path):
        axes = np.array([[0.2, 0.0, 0.0], [0.1, 0.3, 0.0], [0.0, 0.05, 0.4]])
        (volume,) = cf.read_cube_volumes(write_text(tmp_path, cube_text(np.zeros((2, 3, 4)), axes=axes)))
        assert np.allclose(volume.axes, axes)

    def test_singular_axes_are_rejected(self, tmp_path):
        axes = np.array([[0.2, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.0, 0.4]])
        with pytest.raises(ValueError, match="axes"):
            cf.read_cube_volumes(write_text(tmp_path, cube_text(np.zeros((2, 2, 2)), axes=axes)))

    def test_fortran_exponent_without_e_is_repaired(self, tmp_path):
        text = cube_text(np.full((1, 1, 3), 2.0)).replace("2.00000E+00", "1.23400-100", 1)
        (volume,) = cf.read_cube_volumes(write_text(tmp_path, text))
        assert volume.values[0, 0, 0] == pytest.approx(1.234e-100)
        assert volume.values[0, 0, 2] == 2.0

    def test_text_in_the_data_is_a_value_error(self, tmp_path):
        text = cube_text(np.ones((1, 1, 3))).replace("1.00000E+00", "abc", 1)
        with pytest.raises(ValueError, match="v.cube"):
            cf.read_cube_volumes(write_text(tmp_path, text))

    def test_truncated_volume_is_a_value_error(self, tmp_path):
        text = cube_text(np.ones((2, 2, 7)))
        lines = text.splitlines()[:-1]
        with pytest.raises(ValueError, match="truncated"):
            cf.read_cube_volumes(write_text(tmp_path, "\n".join(lines) + "\n"))

    def test_malformed_header_is_a_value_error(self, tmp_path):
        with pytest.raises(ValueError, match="header"):
            cf.read_cube_volumes(write_text(tmp_path, "t\nc\nx y z\n"))

    def test_missing_file_is_a_value_error(self, tmp_path):
        with pytest.raises(ValueError, match="Cannot read"):
            cf.read_cube_volumes(tmp_path / "none.cube")

    def test_round_trip_with_skewed_axes(self, tmp_path):
        rng = np.random.default_rng(11)
        volume = cf.CubeVolume("t", "c", np.array([6, 1]), np.array([6.0, 1.0]),
                               np.array([[0.0, 0.0, 0.0], [2.0, 0.1, 0.0]]), np.array([-3.0, -2.0, -1.0]),
                               np.array([[0.2, 0.0, 0.0], [0.1, 0.3, 0.0], [0.0, 0.0, 0.25]]),
                               rng.normal(size=(3, 4, 8)))
        (back,) = cf.read_cube_volumes(cf.write_cube_volume(tmp_path / "r.cube", volume))
        assert np.allclose(back.axes, volume.axes) and np.allclose(back.origin, volume.origin)
        assert np.allclose(back.coordinates, volume.coordinates)
        assert np.allclose(back.values, volume.values, rtol=1e-5, atol=1e-12)

    def test_orbital_volume_is_written_as_a_plain_cube(self, tmp_path):
        values = np.ones((1, 2, 2, 1))
        (volume,) = cf.read_cube_volumes(write_text(tmp_path, cube_text(values, mos=[4])))
        path = cf.write_cube_volume(tmp_path / "out.cube", volume)
        assert int(open(path).read().splitlines()[2].split()[0]) == 1
        (back,) = cf.read_cube_volumes(path)
        assert back.orbital is None and "MO 4" in back.comment

    def test_header_text_is_written_as_ascii(self, tmp_path):
        (volume,) = cf.read_cube_volumes(cf.write_cube(tmp_path / "a.cub", sample_cube()))
        volume.title = "rho(N) − rho(N−1) × 2 é"
        volume.comment = "line\nbreak"
        path = cf.write_cube_volume(tmp_path / "u.cube", volume)
        title, comment = open(path, encoding="ascii").read().splitlines()[:2]
        assert title == "rho(N) - rho(N-1) x 2 ?" and comment == "line break"

    def test_unwritable_volume_path_is_a_value_error(self, tmp_path):
        (tmp_path / "file").write_text("x")
        (volume,) = cf.read_cube_volumes(cf.write_cube(tmp_path / "a.cub", sample_cube()))
        with pytest.raises(ValueError, match="Could not write"):
            cf.write_cube_volume(tmp_path / "file" / "b.cube", volume)

    def test_bohr_constant_matches_coordinate_converter(self):
        source = (__import__('pathlib').Path(cf.__file__).parent / "coordinateConverter.py").read_text()
        assert f"BOHR_RADIUS_ANGSTROM = {cf.BOHR_RADIUS_ANGSTROM!r}" in source
