"""Tests for espSurfaceLayer: the ESP surface state shared by the ESP tab and the Fukui viewer."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import espSurfaceLayer as el
import espViewer as ev
import volumeData as vd

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


GL_NAMES = ("glEnable", "glDisable", "glColorMaterial", "glMaterialfv", "glMaterialf", "glLightModeli",
            "glEnableClientState", "glDisableClientState", "glVertexPointer", "glNormalPointer",
            "glColorPointer", "glDrawElements", "glDrawArrays", "glBlendFunc", "glDepthMask")


@pytest.fixture
def gl_calls(monkeypatch):
    calls = []
    for name in GL_NAMES:
        monkeypatch.setattr(el, name, lambda *args, _n=name: calls.append((_n, args)))
    return calls


def sphere_surface(layer, radius=3.0):
    """An icosphere surface (Angstrom) with a value gradient, large enough for contour lines."""
    from cubeRender import _icosphere
    unit, triangles = _icosphere(3)
    layer.set_surface(radius * unit, triangles, unit, unit[:, 0] * 0.02, 0.0004)


class TestCanvasSurfaceStyle:
    """The 3D viewer's Display-tab surface style, applied when the layer draws."""

    def test_no_style_draws_as_before(self, gl_calls):
        layer = el.ESPSurfaceLayer()
        sphere_surface(layer)
        layer.set_display_options(opacity=0.5)
        assert layer.surface_style() is None and layer.styled_rgba(layer.get()) is layer.get().rgba
        layer.draw_scene(None)
        assert ("glDepthMask", (el.GL_FALSE,)) in gl_calls      # translucent, as surface_opacity says

    def test_solid_is_opaque_and_transparent_scales_the_tab_opacity(self):
        layer = el.ESPSurfaceLayer()
        sphere_surface(layer)
        layer.set_display_options(opacity=0.5)
        layer.set_surface_style(vd.SURFACE_SOLID)
        solid = layer.styled_rgba(layer.get())
        assert np.allclose(solid[:, 3], 1.0) and np.allclose(solid[:, :3], layer.get().rgba[:, :3])
        assert layer.styled_rgba(layer.get()) is solid                # computed once, not every frame
        layer.set_surface_style(vd.SURFACE_TRANSPARENT, 0.5)
        assert np.allclose(layer.styled_rgba(layer.get())[:, 3], 0.25)   # tab opacity 0.5 x (1 - 0.5)
        layer.set_display_options(opacity=1.0)
        assert np.allclose(layer.styled_rgba(layer.get())[:, 3], 0.5)
        layer.set_surface_style(vd.SURFACE_TRANSPARENT, 0.0)
        assert layer.styled_rgba(layer.get()) is layer.get().rgba      # nothing to change: no copy

    def test_mesh_draws_coloured_contour_lines(self, gl_calls):
        layer = el.ESPSurfaceLayer()
        sphere_surface(layer)
        layer.set_surface_style(vd.SURFACE_MESH)
        assert layer.styled_rgba(layer.get()) is None
        vertices, rgba = layer.mesh_lines(layer.get())
        assert len(vertices) > 50 and len(vertices) % 2 == 0 and rgba.shape == (len(vertices), 4)
        assert np.allclose(np.linalg.norm(vertices, axis=1), 3.0, atol=0.1)
        assert np.allclose(rgba[:, 3], 1.0) and rgba[:, 0].max() > 0.5 and rgba[:, 2].max() > 0.5  # red..blue
        assert layer.mesh_lines(layer.get())[0] is vertices           # cached per surface
        layer.draw_scene(None)
        draws = [args for name, args in gl_calls if name == "glDrawArrays"]
        assert draws == [(el.GL_LINES, 0, len(vertices))]
        assert not [name for name, _ in gl_calls if name == "glDrawElements"]

    def test_unknown_style(self):
        with pytest.raises(ValueError, match="surface style"):
            el.ESPSurfaceLayer().set_surface_style("Wire")


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
