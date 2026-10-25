"""Tests for cubeRender: render items, draw order and the scene layer (no GL context)."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import cubeRender as cr
import volumeData as vd


def triangle_mesh():
    vertices = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    triangles = np.array([[0, 1, 2]])
    return vd.IsoMesh(vertices, triangles, vd.vertex_normals(vertices, triangles))


class TestItems:
    def test_surface_item_is_in_angstrom_with_gl_types(self):
        item = cr.surface_item(triangle_mesh(), (1.0, 0.0, 0.0), 0.5)
        assert item.kind == "surface"
        assert item.vertices.dtype == np.float32 and item.vertices.flags["C_CONTIGUOUS"]
        assert item.indices.dtype == np.uint32 and item.indices.tolist() == [0, 1, 2]
        assert item.vertices[1] == pytest.approx([vd.BOHR_RADIUS_ANGSTROM, 0, 0])
        assert item.rgba.shape == (3, 4) and item.rgba[:, 3] == pytest.approx([0.5] * 3)
        assert item.normals.dtype == np.float32 and not item.opaque

    def test_opacity_is_clamped(self):
        assert cr.surface_item(triangle_mesh(), (1, 1, 1), 5.0).opaque
        assert cr.surface_item(triangle_mesh(), (1, 1, 1), 0.0).rgba[0, 3] == pytest.approx(cr.MIN_OPACITY)

    def test_slice_item(self):
        points = np.zeros((2, 3, 3))
        points[..., 0] = np.arange(2)[:, None]
        points[..., 1] = np.arange(3)[None, :]
        colours = np.tile([0.0, 1.0, 0.0], (6, 1))
        item = cr.slice_item(points, colours, 1.0)
        assert item.kind == "slice" and item.normals is None and item.opaque
        assert item.vertices.shape == (6, 3) and item.indices.size == 2 * 2 * 3
        assert item.vertices[3] == pytest.approx([vd.BOHR_RADIUS_ANGSTROM, 0, 0])

    def test_line_item(self):
        segments = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0], [0.0, 1.0, 2.0]]])
        item = cr.line_item(segments, (0.0, 1.0, 0.0))
        assert item.kind == "lines" and item.normals is None and item.opaque
        assert item.vertices.shape == (4, 3) and item.indices.tolist() == [0, 1, 2, 3]
        assert item.vertices[1] == pytest.approx([vd.BOHR_RADIUS_ANGSTROM, 0, 0])
        assert item.bounds[1] == pytest.approx([vd.BOHR_RADIUS_ANGSTROM, vd.BOHR_RADIUS_ANGSTROM,
                                                2 * vd.BOHR_RADIUS_ANGSTROM])
        empty = cr.line_item(np.zeros((0, 2, 3)), (1, 1, 1))
        assert empty.vertices.shape == (0, 3) and empty.indices.size == 0

    def test_bounds(self):
        item = cr.surface_item(triangle_mesh(), (1, 1, 1), 1.0)
        lo, hi = item.bounds
        assert lo == pytest.approx([0, 0, 0]) and hi == pytest.approx([vd.BOHR_RADIUS_ANGSTROM, 2 * vd.BOHR_RADIUS_ANGSTROM, 0])


class TestDrawOrder:
    def test_opaque_items_come_first_in_their_original_order(self):
        a = cr.surface_item(triangle_mesh(), (1, 0, 0), 0.4)
        b = cr.surface_item(triangle_mesh(), (0, 1, 0), 1.0)
        c = cr.surface_item(triangle_mesh(), (0, 0, 1), 0.7)
        d = cr.surface_item(triangle_mesh(), (1, 1, 0), 1.0)
        assert cr.draw_order([a, b, c, d]) == [b, d, a, c]


class TestLayer:
    def test_empty_layer(self):
        layer = cr.CubeSceneLayer()
        assert layer.items() == () and layer.radius_about((0, 0, 0)) == 0.0
        layer.draw()  # nothing to draw: no GL call is made

    def test_set_and_clear(self):
        layer = cr.CubeSceneLayer()
        item = cr.surface_item(triangle_mesh(), (1, 1, 1), 1.0)
        layer.set_items([item])
        assert layer.items() == (item,)
        layer.clear()
        assert layer.items() == ()

    def test_lines_are_drawn_unlit_as_gl_lines(self, monkeypatch):
        calls = []
        for name in ("glEnable", "glDisable", "glColorMaterial", "glMaterialfv", "glMaterialf", "glLightModeli",
                     "glEnableClientState", "glDisableClientState", "glVertexPointer", "glNormalPointer",
                     "glColorPointer", "glDrawElements", "glBlendFunc", "glDepthMask"):
            monkeypatch.setattr(cr, name, lambda *args, _n=name: calls.append((_n, args)))
        layer = cr.CubeSceneLayer()
        layer.set_items([cr.line_item(np.zeros((3, 2, 3)), (1, 1, 1))])
        layer.draw()
        draws = [args for name, args in calls if name == "glDrawElements"]
        assert len(draws) == 1 and draws[0][:2] == (cr.GL_LINES, 6)
        assert calls[0] == ("glDisable", (cr.GL_LIGHTING,)) and calls[-1] == ("glEnable", (cr.GL_LIGHTING,))
        assert "glNormalPointer" not in [name for name, _ in calls]

    def test_radius_covers_every_item(self):
        layer = cr.CubeSceneLayer()
        layer.set_items([cr.surface_item(triangle_mesh(), (1, 1, 1), 1.0)])
        a = vd.BOHR_RADIUS_ANGSTROM
        assert layer.radius_about((0.0, 0.0, 0.0)) == pytest.approx(np.hypot(a, 2 * a))


def test_module_imports_neither_toga_nor_glfw():
    code = "import sys, cubeRender; print('toga' in sys.modules, 'glfw' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(Path(cr.__file__).parent), check=True)
    assert out.stdout.strip() == "False False"


class TestColourBar:
    def test_bar_travels_with_the_items_and_is_drawn_at_the_bottom(self, monkeypatch):
        calls = []
        monkeypatch.setattr(cr, "draw_colour_bar", lambda viewer, *args, **kwargs: calls.append((args, kwargs)))
        layer = cr.CubeSceneLayer()
        layer.draw_overlay("viewer")
        assert calls == []  # no bar: nothing drawn
        bar = ("-2.000e-2", "2.000e-2", "ESP mapped on rho")
        layer.set_items([cr.surface_item(triangle_mesh(), (1, 0, 0), 1.0)], colour_bar=bar)
        assert layer.colour_bar() == bar
        layer.draw_overlay("viewer")
        assert calls == [(bar, {"at_top": False})]
        layer.set_items([])
        assert layer.colour_bar() is None
        layer.set_items([], colour_bar=bar)
        layer.clear()
        assert layer.colour_bar() is None and layer.items() == ()


class TestMarkers:
    def test_sphere_vertices_lie_on_the_radius_with_outward_normals(self):
        centre = np.array([1.0, -2.0, 0.5])
        item = cr.sphere_item(centre, 0.4, (1.0, 0.0, 0.0))
        bohr = item.vertices.astype(float) / vd.BOHR_RADIUS_ANGSTROM
        assert np.allclose(np.linalg.norm(bohr - centre, axis=1), 0.4, atol=1e-6)
        assert np.allclose(item.normals, (bohr - centre) / 0.4, atol=1e-5)
        assert item.kind == "surface" and item.indices.size % 3 == 0 and item.opaque

    def test_arrow_spans_start_to_end(self):
        start, end = np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 4.0])
        item = cr.arrow_item(start, end, 0.1, (1.0, 0.5, 0.0))
        bohr = item.vertices.astype(float) / vd.BOHR_RADIUS_ANGSTROM
        assert abs(bohr[:, 2].min()) < 1e-6 and abs(bohr[:, 2].max() - 4.0) < 1e-6
        assert np.abs(bohr[:, :2]).max() <= 0.2 + 1e-6          # cone radius 2 x 0.1
        assert np.allclose(np.linalg.norm(item.normals, axis=1), 1.0, atol=1e-5)

    def test_zero_length_arrow_is_none(self):
        assert cr.arrow_item(np.ones(3), np.ones(3), 0.1, (1.0, 1.0, 1.0)) is None
