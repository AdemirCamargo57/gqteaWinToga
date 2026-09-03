"""Tests for the molecularViewer atom-picking geometry.

These exercise the pure projection / nearest-atom math used by the
click-to-identify feature. They construct known OpenGL-style matrices and
feed them directly, so no GL context (and no GUI) is required — matching the
headless test convention in tests/conftest.py.
"""
import numpy as np
import pytest

from molecularViewer import MolecularViewer


def _ortho_matrix(left, right, bottom, top, near, far):
    """glOrtho matrix in the column-major layout glGetDoublev returns.

    Stored so that ``obj_vec @ M`` reproduces the OpenGL transform.
    """
    return np.array(
        [
            [2.0 / (right - left), 0.0, 0.0, 0.0],
            [0.0, 2.0 / (top - bottom), 0.0, 0.0],
            [0.0, 0.0, -2.0 / (far - near), 0.0],
            [
                -(right + left) / (right - left),
                -(top + bottom) / (top - bottom),
                -(far + near) / (far - near),
                1.0,
            ],
        ],
        dtype=float,
    )


def _viewer_with_view():
    viewer = MolecularViewer()
    viewer.frames = []
    viewer.molecule_data = [
        ("C", (5.0, 0.0, 0.0)),
        ("O", (-5.0, 0.0, 0.0)),
        ("N", (0.0, 0.0, 0.0)),
    ]
    identity = np.eye(4)
    proj = _ortho_matrix(-10, 10, -10, 10, -10, 10)
    viewer._pick_view = (identity, proj, (0, 0, 200, 200))
    return viewer


def test_project_atom_maps_to_expected_pixels():
    viewer = _viewer_with_view()
    identity = np.eye(4)
    proj = _ortho_matrix(-10, 10, -10, 10, -10, 10)

    sx, sy, sz = viewer._project_atom((5.0, 0.0, 0.0), identity, proj, (0, 0, 200, 200))

    assert sx == pytest.approx(150.0)
    assert sy == pytest.approx(100.0)
    assert 0.0 <= sz <= 1.0


def test_pick_selects_atom_under_cursor():
    viewer = _viewer_with_view()

    # Cursor uses GLFW top-left origin; window height is 200.
    # C -> window (150, 100) -> cursor (150, 100)
    assert viewer._pick_atom_at(150, 100) == 0
    # O -> window (50, 100) -> cursor (50, 100)
    assert viewer._pick_atom_at(50, 100) == 1
    # N -> window (100, 100) -> cursor (100, 100)
    assert viewer._pick_atom_at(100, 100) == 2


def test_pick_empty_space_returns_none():
    viewer = _viewer_with_view()

    assert viewer._pick_atom_at(10, 10) is None


def test_pick_returns_none_without_captured_view():
    viewer = MolecularViewer()
    viewer._pick_view = None

    assert viewer._pick_atom_at(100, 100) is None


def test_pick_prefers_nearer_atom_on_overlap():
    viewer = MolecularViewer()
    viewer.frames = []
    # Both atoms project to the same pixel but differ in depth.
    viewer.molecule_data = [
        ("C", (0.0, 0.0, -5.0)),
        ("O", (0.0, 0.0, 5.0)),
    ]
    identity = np.eye(4)
    proj = _ortho_matrix(-10, 10, -10, 10, -10, 10)
    viewer._pick_view = (identity, proj, (0, 0, 200, 200))

    # O (z=5) maps to the smaller window-z (nearer), so it wins.
    assert viewer._pick_atom_at(100, 100) == 1
