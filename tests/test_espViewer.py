"""Headless tests for the ESP viewer: colour maths and surface state (no GL window)."""
import numpy as np
import pytest

import espViewer as ev
from molecularViewer import MolecularViewer


class TestColours:
    def test_hue_sweep_end_points(self):
        rgb = ev.hue_colours(np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
        assert np.allclose(rgb, [[1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1]])

    def test_esp_zero_is_green_and_clipped(self):
        rgb = ev.esp_colours(np.array([-0.05, 0.0, 0.05, 0.2]), 0.05)
        assert np.allclose(rgb[0], [1, 0, 0]) and np.allclose(rgb[1], [0, 1, 0])
        assert np.allclose(rgb[2], [0, 0, 1]) and np.allclose(rgb[3], [0, 0, 1])

    def test_symmetric_range_is_max_abs(self):
        assert ev.symmetric_range(np.array([-0.02097, 0.01166])) == pytest.approx(0.02097)
        assert ev.symmetric_range(np.zeros(3)) == 1.0  # never a zero-width bar

    @pytest.mark.parametrize("value, text", [(-0.02072, "-2.072e-2"), (0.02072, "2.072e-2"),
                                             (-13.157, "-1.316e1"), (0.0, "0.000e0")])
    def test_gaussview_label(self, value, text):
        assert ev.gaussview_label(value) == text

    def test_colour_bar_text_names_both_quantities(self):
        left, right, caption = ev.colour_bar_text(0.02097, ev.UNIT_KCAL, 0.0004)
        assert left == "-1.316e1" and right == "1.316e1"
        assert "ESP (kcal/mol)" in caption and "rho = 0.000400 a.u." in caption

    def test_unknown_unit(self):
        with pytest.raises(ValueError):
            ev.unit_factor("eV")


class TestViewerState:
    @pytest.fixture
    def viewer(self):
        v = ev.ESPSurfaceViewer()
        v.set_molecule(["C", "C"], np.array([[0, 0, -0.67], [0, 0, 0.67]]))
        return v

    def test_base_viewer_hooks_are_no_ops(self):
        base = MolecularViewer()
        assert base.render_scene_extras() is None and base.draw_screen_overlay() is None
        assert base.window_title == "3D Molecular Viewer"

    def test_set_surface_builds_draw_data(self, viewer):
        verts = np.array([[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]])
        viewer.set_surface(verts, np.array([[0, 1, 2]]), np.eye(3), np.array([-0.02, 0.0, 0.01]), 0.0004)
        s = viewer.get_surface()
        assert s.vertices.dtype == np.float32 and s.indices.dtype == np.uint32
        assert s.half_range == pytest.approx(0.02)
        assert np.allclose(s.rgba[:, 3], 1.0)
        assert viewer.get_scene_radius() >= 3.0  # the camera frames the surface, not just the atoms
        assert viewer.colour_bar_text()[0] == "-2.000e-2"

    def test_opacity_updates_alpha(self, viewer):
        viewer.set_surface(np.zeros((3, 3)) + [[1, 0, 0]], np.array([[0, 1, 2]]), np.eye(3), np.zeros(3), 0.001)
        viewer.set_display_options(opacity=0.5)
        assert np.allclose(viewer.get_surface().rgba[:, 3], 0.5)

    def test_clear_surface(self, viewer):
        viewer.set_surface(np.eye(3), np.array([[0, 1, 2]]), np.eye(3), np.zeros(3), 0.001)
        viewer.clear_surface()
        assert viewer.get_surface() is None and viewer.colour_bar_text() is None
