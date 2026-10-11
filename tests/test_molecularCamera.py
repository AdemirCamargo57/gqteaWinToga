"""Tests for the design canvas view camera.

Pure numpy: the camera maps sketch coordinates to screen pixels and back, and
these tests pin the conventions the canvas relies on -- the default view is the
identity, depth grows away from the viewer, a drag rotates the front of the
molecule the way the pointer moves, and changing the pivot never makes the
picture jump.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from molecularCamera import ViewCamera


def _points():
    return np.array(
        [[100.0, 120.0, 0.0], [160.0, 80.0, 30.0], [40.0, 200.0, -45.0]], dtype=float
    )


# ------------------------------------------------------------------
# Default view
# ------------------------------------------------------------------
def test_the_default_view_is_the_identity():
    """A fresh canvas must behave exactly as the flat canvas did."""
    camera = ViewCamera()

    projected = camera.project(_points())

    assert np.allclose(projected, _points())


def test_unproject_inverts_project():
    camera = ViewCamera()
    camera.rotate_by_drag(37.0, -21.0)
    camera.pan_by(12.0, -8.0)
    camera.zoom_by(1.7, 280.0, 210.0)

    for point in _points():
        sx, sy, depth = camera.project(point[None, :])[0]
        assert np.allclose(camera.unproject(sx, sy, depth), point)


def test_project_accepts_an_empty_array():
    assert ViewCamera().project(np.zeros((0, 3))).shape == (0, 3)


# ------------------------------------------------------------------
# Rotation
# ------------------------------------------------------------------
def test_dragging_right_moves_the_front_of_the_molecule_right():
    camera = ViewCamera()
    front = np.array([[0.0, 0.0, -50.0]])  # negative depth = towards the viewer

    camera.rotate_by_drag(30.0, 0.0)

    assert camera.project(front)[0][0] > 0.0


def test_dragging_down_moves_the_front_of_the_molecule_down():
    camera = ViewCamera()
    front = np.array([[0.0, 0.0, -50.0]])

    camera.rotate_by_drag(0.0, 30.0)

    assert camera.project(front)[0][1] > 0.0


def test_rotation_preserves_distances():
    camera = ViewCamera()
    before = _points()
    for _ in range(200):
        camera.rotate_by_drag(13.0, -7.0)

    projected = camera.project(before)

    for i in range(3):
        for j in range(i + 1, 3):
            assert np.linalg.norm(projected[i] - projected[j]) == pytest.approx(
                np.linalg.norm(before[i] - before[j])
            )


def test_the_rotation_stays_orthonormal_after_many_drags():
    camera = ViewCamera()
    for _ in range(1000):
        camera.rotate_by_drag(3.1, 2.7)

    assert np.allclose(camera.rotation @ camera.rotation.T, np.eye(3), atol=1e-12)
    assert np.linalg.det(camera.rotation) == pytest.approx(1.0)


def test_rotation_turns_about_the_pivot():
    camera = ViewCamera()
    pivot = np.array([200.0, 150.0, 0.0])
    camera.set_pivot(pivot)

    camera.rotate_by_drag(45.0, 20.0)

    assert np.allclose(camera.project(pivot[None, :])[0][:2], pivot[:2])


def test_moving_the_pivot_does_not_move_the_picture():
    camera = ViewCamera()
    camera.rotate_by_drag(25.0, -15.0)
    before = camera.project(_points())

    camera.set_pivot(np.array([90.0, 40.0, 10.0]))

    assert np.allclose(camera.project(_points())[:, :2], before[:, :2])


# ------------------------------------------------------------------
# Pan and zoom
# ------------------------------------------------------------------
def test_panning_shifts_every_point_by_the_drag():
    camera = ViewCamera()

    camera.pan_by(15.0, -5.0)

    assert np.allclose(camera.project(_points())[:, :2], _points()[:, :2] + [15.0, -5.0])


def test_zoom_keeps_the_point_under_its_centre_fixed():
    camera = ViewCamera()
    centre = (280.0, 210.0)

    camera.zoom_by(2.0, *centre)

    assert np.allclose(camera.project(np.array([[280.0, 210.0, 0.0]]))[0][:2], centre)
    assert camera.zoom == pytest.approx(2.0)


def test_zoom_scales_on_screen_distances():
    camera = ViewCamera()

    camera.zoom_by(1.5, 0.0, 0.0)
    projected = camera.project(np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]]))

    assert projected[1][0] - projected[0][0] == pytest.approx(150.0)


def test_zoom_is_clamped():
    camera = ViewCamera()

    for _ in range(50):
        camera.zoom_by(2.0, 0.0, 0.0)

    assert camera.zoom == pytest.approx(ViewCamera.MAX_ZOOM)


def test_a_clamped_zoom_does_not_drift_the_picture():
    camera = ViewCamera()
    for _ in range(50):
        camera.zoom_by(2.0, 280.0, 210.0)
    before = camera.project(_points())

    camera.zoom_by(2.0, 280.0, 210.0)

    assert np.allclose(camera.project(_points()), before)


def test_reset_restores_the_identity():
    camera = ViewCamera()
    camera.rotate_by_drag(40.0, 40.0)
    camera.pan_by(10.0, 10.0)
    camera.zoom_by(3.0, 5.0, 5.0)
    camera.set_pivot(np.array([1.0, 2.0, 3.0]))

    camera.reset()

    assert np.allclose(camera.project(_points()), _points())


def test_the_camera_never_imports_toga():
    source_dir = Path(__file__).resolve().parents[1]
    code = (
        "import sys; sys.path.insert(0, %r); import molecularCamera; "
        "print('toga' in sys.modules)" % str(source_dir)
    )
    output = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert output.stdout.strip() == "False"
