"""Tests for espSurfaceLayer: the ESP surface state shared by the ESP tab and the Fukui viewer."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import espSurfaceLayer as el
import espViewer as ev

TRIANGLE = np.array([[0, 1, 2]])


def surface(layer, values=(-0.02, 0.0, 0.01), centre=(0.0, 0.0, 0.0), **kwargs):
    verts = np.array([[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]])
    layer.set_surface(verts, TRIANGLE, np.eye(3), np.array(values), 0.0004, centre=centre, **kwargs)


class TestLayer:
    def test_empty_layer(self):
        layer = el.ESPSurfaceLayer()
        assert layer.get() is None and layer.colour_bar_text() is None and layer.radius() == 0.0

    def test_set_surface_builds_draw_data(self):
        layer = el.ESPSurfaceLayer()
        surface(layer, markers=[el.SurfaceMarker(np.zeros(3), "minimum", "global", "1")])
        s = layer.get()
        assert s.vertices.dtype == np.float32 and s.indices.dtype == np.uint32
        assert s.half_range == pytest.approx(0.02) and np.allclose(s.rgba[:, 3], 1.0)
        assert s.radius == pytest.approx(3.0) and layer.radius() == pytest.approx(3.0)
        assert len(s.markers) == 1
        assert layer.colour_bar_text()[0] == "-2.000e-2"

    def test_radius_is_measured_from_the_given_centre(self):
        layer = el.ESPSurfaceLayer()
        surface(layer, centre=(3.0, 0.0, 0.0))
        assert layer.radius() == pytest.approx(np.sqrt(18.0))

    def test_display_options(self):
        layer = el.ESPSurfaceLayer()
        surface(layer)
        layer.set_display_options(unit=el.UNIT_KCAL, opacity=0.5, show_extrema=False)
        assert np.allclose(layer.get().rgba[:, 3], 0.5)
        assert layer.colour_unit == el.UNIT_KCAL and layer.show_extrema is False
        assert "kcal/mol" in layer.colour_bar_text()[2]
        layer.set_display_options(opacity=0.0)
        assert layer.surface_opacity == pytest.approx(0.05)
        with pytest.raises(ValueError):
            layer.set_display_options(unit="eV")

    def test_explicit_range_and_caption(self):
        layer = el.ESPSurfaceLayer()
        surface(layer, half_range=0.5, caption="f+ (a.u.)")
        assert layer.get().half_range == 0.5 and layer.colour_bar_text()[2] == "f+ (a.u.)"

    def test_clear(self):
        layer = el.ESPSurfaceLayer()
        surface(layer)
        layer.clear()
        assert layer.get() is None and layer.radius() == 0.0

    def test_drawing_an_empty_layer_makes_no_gl_call(self):
        layer = el.ESPSurfaceLayer()
        layer.draw_scene(None)  # nothing to draw: the viewer is never touched
        layer.draw_overlay(None)


class TestEspViewerReExports:
    def test_moved_names_are_still_importable_from_esp_viewer(self):
        for name in ("hue_colours", "symmetric_range", "esp_colours", "gaussview_label", "unit_factor",
                     "colour_bar_text", "contrasting_text_colour", "SurfaceMarker", "SurfaceDrawData",
                     "UNIT_AU", "UNIT_KCAL", "COLOUR_UNITS"):
            assert getattr(ev, name) is getattr(el, name)

    def test_esp_viewer_delegates_to_its_layer(self):
        v = ev.ESPSurfaceViewer()
        v.set_molecule(["C"], np.zeros((1, 3)))
        v.set_surface(np.array([[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]]), TRIANGLE, np.eye(3), np.zeros(3), 0.001)
        assert v.get_surface() is v.esp_layer.get()
        v.surface_opacity = 0.3
        assert v.esp_layer.surface_opacity == 0.3


def test_module_imports_neither_toga_nor_glfw_nor_the_viewer():
    code = ("import sys, espSurfaceLayer; "
            "print('toga' in sys.modules, 'glfw' in sys.modules, 'molecularViewer' in sys.modules)")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(Path(el.__file__).parent), check=True)
    assert out.stdout.strip() == "False False False"


def test_the_esp_colour_bar_is_drawn_at_the_top(monkeypatch):
    calls = []
    monkeypatch.setattr(el, "draw_colour_bar", lambda viewer, *args, **kwargs: calls.append((args, kwargs)))
    layer = el.ESPSurfaceLayer()
    surface(layer)
    layer.draw_overlay("viewer")
    assert calls == [(layer.colour_bar_text(), {"at_top": True})]
