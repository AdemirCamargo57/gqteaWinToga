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
