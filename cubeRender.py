"""OpenGL drawing of the Cube tab's isosurfaces and slices.

GL only: no Toga, no GLFW, so the item construction and draw order are
tested headless. The Toga thread builds immutable ``RenderItem``s and hands
them to a ``CubeSceneLayer`` under its lock; the render thread calls
``CubeSceneLayer.draw()`` from ``MolecularViewerUI.render_scene_extras``,
inside the molecule's push-matrix, so the volume is drawn in the atoms' Å
frame and turns with the molecule.

Surfaces use the ``espViewer._draw_surface`` recipe (vertex arrays with
per-vertex colour, two-sided lighting, a specular highlight, blending with
depth writes off when translucent); slices are drawn unlit. Opaque items
are drawn before translucent ones; translucent items are not depth-sorted
against each other.
"""

from __future__ import annotations

import threading
from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from OpenGL.GL import *  # noqa: F401,F403 - same style as molecularViewer/espViewer

from espSurfaceLayer import draw_colour_bar
from volumeData import (BOHR_RADIUS_ANGSTROM, DEFAULT_TRANSPARENCY, MESH_LINE_SPACING_BOHR, SURFACE_SOLID,  # noqa: F401
                        SURFACE_STYLES, SURFACE_TRANSPARENT, IsoMesh, checked_surface_style, contour_lines,
                        slice_triangles)

MIN_OPACITY = 0.05
OPAQUE_ALPHA = 0.999


@dataclass(frozen=True, eq=False)
class RenderItem:
    kind: str                # "surface" | "slice" | "lines"
    vertices: np.ndarray     # (V, 3) float32, Angstrom
    indices: np.ndarray      # (3F,) uint32; (2S,) for lines
    rgba: np.ndarray         # (V, 4) float32
    normals: Optional[np.ndarray] = None  # (V, 3) float32, surfaces only
    bounds: Tuple[np.ndarray, np.ndarray] = field(default=None)
    isosurface: bool = False  # an isosurface (surface_item): follows the canvas surface style; markers do not

    @property
    def opaque(self) -> bool:
        return not len(self.rgba) or float(self.rgba[:, 3].min()) >= OPAQUE_ALPHA


def _clamp_opacity(opacity: float) -> float:
    return float(min(1.0, max(MIN_OPACITY, opacity)))


def _rgba(rgb: np.ndarray, n: int, opacity: float) -> np.ndarray:
    rgb = np.broadcast_to(np.asarray(rgb, dtype=float), (n, 3))
    return np.ascontiguousarray(np.hstack([rgb, np.full((n, 1), _clamp_opacity(opacity))]), dtype=np.float32)


def _item(kind, vertices_bohr, triangles, rgba, normals=None) -> RenderItem:
    vertices = np.ascontiguousarray(np.asarray(vertices_bohr, dtype=float) * BOHR_RADIUS_ANGSTROM,
                                    dtype=np.float32)
    bounds = ((vertices.min(axis=0), vertices.max(axis=0)) if len(vertices)
              else (np.zeros(3, np.float32), np.zeros(3, np.float32)))
    return RenderItem(
        kind=kind,
        vertices=vertices,
        indices=np.ascontiguousarray(np.asarray(triangles).ravel(), dtype=np.uint32),
        rgba=rgba,
        normals=None if normals is None else np.ascontiguousarray(normals, dtype=np.float32),
        bounds=bounds,
    )


def surface_item(mesh: IsoMesh, rgb, opacity: float) -> RenderItem:
    """One isosurface in a single colour."""
    item = _item("surface", mesh.vertices, mesh.triangles, _rgba(rgb, len(mesh.vertices), opacity),
                 mesh.normals)
    return replace(item, isosurface=True)


def slice_item(points_bohr: np.ndarray, colours: np.ndarray, opacity: float) -> RenderItem:
    """A slice plane: points (n1, n2, 3) bohr and one RGB per point (n1*n2, 3)."""
    n1, n2 = points_bohr.shape[:2]
    return _item("slice", points_bohr.reshape(-1, 3), slice_triangles(n1, n2),
                 _rgba(colours, n1 * n2, opacity))


def _icosphere(subdivisions: int) -> Tuple[np.ndarray, np.ndarray]:
    """Unit-sphere vertices and triangles: an icosahedron split ``subdivisions`` times."""
    t = (1.0 + 5 ** 0.5) / 2.0
    v = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
         (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    f = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2),
         (10, 7, 6), (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11),
         (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    vertices = [np.array(p, dtype=float) / np.linalg.norm(p) for p in v]
    for _ in range(subdivisions):
        cache, faces = {}, []

        def middle(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                p = vertices[a] + vertices[b]
                vertices.append(p / np.linalg.norm(p))
                cache[key] = len(vertices) - 1
            return cache[key]

        for a, b, c in f:
            ab, bc, ca = middle(a, b), middle(b, c), middle(c, a)
            faces += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        f = faces
    return np.array(vertices), np.array(f, dtype=np.uint32)


def sphere_item(centre_bohr, radius_bohr: float, rgb, opacity: float = 1.0, subdivisions: int = 2) -> RenderItem:
    """A lit sphere (e.g. a hole or electron centroid marker)."""
    unit, triangles = _icosphere(subdivisions)
    vertices = np.asarray(centre_bohr, dtype=float) + radius_bohr * unit
    return _item("surface", vertices, triangles, _rgba(rgb, len(vertices), opacity), unit)


def arrow_item(start_bohr, end_bohr, radius_bohr: float, rgb, opacity: float = 1.0,
               sides: int = 16) -> Optional[RenderItem]:
    """A shaft (80 % of the length) and a cone (radius 2x) from ``start`` to ``end``; None for zero length."""
    start, end = np.asarray(start_bohr, dtype=float), np.asarray(end_bohr, dtype=float)
    axis = end - start
    length = float(np.linalg.norm(axis))
    if length < 1e-9:
        return None
    w = axis / length
    helper = np.array([1.0, 0.0, 0.0]) if abs(w[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(w, helper)
    u /= np.linalg.norm(u)
    v = np.cross(w, u)
    angles = 2 * np.pi * np.arange(sides) / sides
    ring = np.outer(np.cos(angles), u) + np.outer(np.sin(angles), v)       # (sides, 3) unit radial
    neck = start + 0.8 * axis
    vertices = np.vstack([start + radius_bohr * ring, neck + radius_bohr * ring,  # shaft
                          neck + 2 * radius_bohr * ring, np.repeat(end[None], sides, axis=0)])  # cone
    slope = 2 * radius_bohr / (0.2 * length)
    cone_normal = ring + slope * w
    cone_normal /= np.linalg.norm(cone_normal, axis=1)[:, None]
    normals = np.vstack([ring, ring, cone_normal, cone_normal])
    triangles = []
    for k in range(sides):
        n = (k + 1) % sides
        triangles += [(k, n, sides + n), (k, sides + n, sides + k)]                        # shaft
        triangles += [(2 * sides + k, 2 * sides + n, 3 * sides + k)]                       # cone
    return _item("surface", vertices, np.array(triangles, dtype=np.uint32),
                 _rgba(rgb, len(vertices), opacity), normals)


def draw_order(items: Iterable[RenderItem]) -> List[RenderItem]:
    """Opaque items first, then translucent ones, each group in the given order."""
    items = list(items)
    return [i for i in items if i.opaque] + [i for i in items if not i.opaque]


def mesh_lines_item(item: RenderItem) -> RenderItem:
    """The Mesh style of one isosurface: its contour lines (planes MESH_LINE_SPACING_BOHR apart), in its colours."""
    mesh = IsoMesh(np.asarray(item.vertices, dtype=float), item.indices.reshape(-1, 3), None)
    segments, colours = contour_lines(mesh, np.zeros(3), MESH_LINE_SPACING_BOHR * BOHR_RADIUS_ANGSTROM,
                                      values=item.rgba[:, :3])
    vertices = np.ascontiguousarray(segments.reshape(-1, 3), dtype=np.float32)        # already Angstrom
    rgba = np.ascontiguousarray(np.hstack([colours.reshape(-1, 3), np.ones((len(vertices), 1))]), dtype=np.float32)
    bounds = ((vertices.min(axis=0), vertices.max(axis=0)) if len(vertices)
              else (np.zeros(3, np.float32), np.zeros(3, np.float32)))
    return RenderItem("lines", vertices, np.arange(len(vertices), dtype=np.uint32), rgba, None, bounds)


def styled_items(items: Sequence[RenderItem], style: str, transparency: float,
                 mesh_cache: Optional[Dict[int, Tuple[RenderItem, RenderItem]]] = None) -> List[RenderItem]:
    """The items as the canvas surface style draws them; only isosurfaces change, markers never.

    Solid: fully opaque. Transparent: alpha x (1 - transparency), so the
    tab's own opacity still scales its surfaces. Mesh: contour lines
    (``mesh_lines_item``), cached in ``mesh_cache`` by item identity.
    """
    out = []
    for item in items:
        if not item.isosurface:
            out.append(item)
        elif style == SURFACE_SOLID:
            out.append(item if item.opaque else replace(item, rgba=_with_alpha(item.rgba, 1.0)))
        elif style == SURFACE_TRANSPARENT:
            out.append(item if transparency <= 0.0 else
                       replace(item, rgba=_with_alpha(item.rgba, item.rgba[:, 3] * (1.0 - transparency))))
        else:
            cached = mesh_cache.get(id(item)) if mesh_cache is not None else None
            if cached is None or cached[0] is not item:
                cached = (item, mesh_lines_item(item))
                if mesh_cache is not None:
                    mesh_cache[id(item)] = cached
            out.append(cached[1])
    return out


def _with_alpha(rgba: np.ndarray, alpha) -> np.ndarray:
    out = np.array(rgba, dtype=np.float32, copy=True)
    out[:, 3] = np.maximum(MIN_OPACITY, alpha)
    return out


class CubeSceneLayer:
    """What the Cube tab wants drawn, shared between the Toga and render threads."""

    def __init__(self):
        self._lock = threading.Lock()
        self._items: Tuple[RenderItem, ...] = ()
        self._colour_bar: Optional[Tuple[str, str, str]] = None
        self._style: Optional[Tuple[str, float]] = None   # None: draw items as built (the layer's old behaviour)
        self._styled = None                               # (items, style, styled items) cache
        self._mesh_cache: Dict[int, Tuple[RenderItem, RenderItem]] = {}

    # --- Toga thread ---------------------------------------------------
    def set_items(self, items: Sequence[RenderItem], colour_bar: Optional[Tuple[str, str, str]] = None):
        """The scene, and the (left label, right label, caption) of a mapped surface's colour bar."""
        with self._lock:
            self._items = tuple(items)
            self._colour_bar = tuple(colour_bar) if colour_bar else None
            self._styled = None
            self._mesh_cache = {}

    def set_surface_style(self, style: Optional[str], transparency: float = DEFAULT_TRANSPARENCY):
        """The canvas surface style for this layer's isosurfaces (None: as built)."""
        checked = None if style is None else checked_surface_style(style, transparency)
        with self._lock:
            self._style = checked
            self._styled = None

    def surface_style(self) -> Optional[Tuple[str, float]]:
        with self._lock:
            return self._style

    def clear(self):
        self.set_items(())

    def items(self) -> Tuple[RenderItem, ...]:
        with self._lock:
            return self._items

    def colour_bar(self) -> Optional[Tuple[str, str, str]]:
        with self._lock:
            return self._colour_bar

    def radius_about(self, centre) -> float:
        """Distance from ``centre`` (Angstrom) to the farthest bounding-box corner of any item."""
        centre = np.asarray(centre, dtype=float)
        radius = 0.0
        for item in self.items():
            lo, hi = (np.asarray(b, dtype=float) for b in item.bounds)
            far = np.maximum(np.abs(lo - centre), np.abs(hi - centre))
            radius = max(radius, float(np.linalg.norm(far)))
        return radius

    # --- render thread -------------------------------------------------
    def draw_overlay(self, viewer):
        """The mapped surface's colour bar, at the bottom (the ESP tab's bar is at the top)."""
        bar = self.colour_bar()
        if bar is not None:
            draw_colour_bar(viewer, *bar, at_top=False)

    def drawn_items(self) -> Tuple[RenderItem, ...]:
        """The items with the surface style applied (computed once per items/style pair)."""
        with self._lock:
            items, style, styled, mesh_cache = self._items, self._style, self._styled, self._mesh_cache
        if style is None or not items:
            return items
        if styled is not None and styled[0] is items and styled[1] == style:
            return styled[2]
        result = tuple(styled_items(items, style[0], style[1], mesh_cache))
        with self._lock:
            if self._items is items and self._style == style:
                self._styled = (items, style, result)
        return result

    def draw(self):
        items = self.drawn_items()
        if not items:
            return
        for item in draw_order(items):
            if item.kind == "surface":
                self._draw_surface(item)
            elif item.kind == "lines":
                self._draw_lines(item)
            else:
                self._draw_slice(item)

    @staticmethod
    def _draw_surface(item: RenderItem):
        translucent = not item.opaque
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (0.6, 0.6, 0.6, 1.0))
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 40.0)
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
        if translucent:
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glDepthMask(GL_FALSE)
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_NORMAL_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        try:
            glVertexPointer(3, GL_FLOAT, 0, item.vertices)
            glNormalPointer(GL_FLOAT, 0, item.normals)
            glColorPointer(4, GL_FLOAT, 0, item.rgba)
            glDrawElements(GL_TRIANGLES, int(item.indices.size), GL_UNSIGNED_INT, item.indices)
        finally:
            glDisableClientState(GL_COLOR_ARRAY)
            glDisableClientState(GL_NORMAL_ARRAY)
            glDisableClientState(GL_VERTEX_ARRAY)
            if translucent:
                glDepthMask(GL_TRUE)
                glDisable(GL_BLEND)
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (0.0, 0.0, 0.0, 1.0))
            glDisable(GL_COLOR_MATERIAL)

    @staticmethod
    def _draw_lines(item: RenderItem):
        translucent = not item.opaque
        glDisable(GL_LIGHTING)
        if translucent:
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glDepthMask(GL_FALSE)
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        try:
            glVertexPointer(3, GL_FLOAT, 0, item.vertices)
            glColorPointer(4, GL_FLOAT, 0, item.rgba)
            glDrawElements(GL_LINES, int(item.indices.size), GL_UNSIGNED_INT, item.indices)
        finally:
            glDisableClientState(GL_COLOR_ARRAY)
            glDisableClientState(GL_VERTEX_ARRAY)
            if translucent:
                glDepthMask(GL_TRUE)
                glDisable(GL_BLEND)
            glEnable(GL_LIGHTING)

    @staticmethod
    def _draw_slice(item: RenderItem):
        translucent = not item.opaque
        glDisable(GL_LIGHTING)
        if translucent:
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glDepthMask(GL_FALSE)
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        try:
            glVertexPointer(3, GL_FLOAT, 0, item.vertices)
            glColorPointer(4, GL_FLOAT, 0, item.rgba)
            glDrawElements(GL_TRIANGLES, int(item.indices.size), GL_UNSIGNED_INT, item.indices)
        finally:
            glDisableClientState(GL_COLOR_ARRAY)
            glDisableClientState(GL_VERTEX_ARRAY)
            if translucent:
                glDepthMask(GL_TRUE)
                glDisable(GL_BLEND)
            glEnable(GL_LIGHTING)
