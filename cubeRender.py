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
from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
from OpenGL.GL import *  # noqa: F401,F403 - same style as molecularViewer/espViewer

from volumeData import BOHR_RADIUS_ANGSTROM, IsoMesh, slice_triangles

MIN_OPACITY = 0.05
OPAQUE_ALPHA = 0.999


@dataclass(frozen=True, eq=False)
class RenderItem:
    kind: str                # "surface" | "slice"
    vertices: np.ndarray     # (V, 3) float32, Angstrom
    indices: np.ndarray      # (3F,) uint32
    rgba: np.ndarray         # (V, 4) float32
    normals: Optional[np.ndarray] = None  # (V, 3) float32, surfaces only
    bounds: Tuple[np.ndarray, np.ndarray] = field(default=None)

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
    return _item("surface", mesh.vertices, mesh.triangles, _rgba(rgb, len(mesh.vertices), opacity),
                 mesh.normals)


def slice_item(points_bohr: np.ndarray, colours: np.ndarray, opacity: float) -> RenderItem:
    """A slice plane: points (n1, n2, 3) bohr and one RGB per point (n1*n2, 3)."""
    n1, n2 = points_bohr.shape[:2]
    return _item("slice", points_bohr.reshape(-1, 3), slice_triangles(n1, n2),
                 _rgba(colours, n1 * n2, opacity))


def draw_order(items: Iterable[RenderItem]) -> List[RenderItem]:
    """Opaque items first, then translucent ones, each group in the given order."""
    items = list(items)
    return [i for i in items if i.opaque] + [i for i in items if not i.opaque]


class CubeSceneLayer:
    """What the Cube tab wants drawn, shared between the Toga and render threads."""

    def __init__(self):
        self._lock = threading.Lock()
        self._items: Tuple[RenderItem, ...] = ()

    # --- Toga thread ---------------------------------------------------
    def set_items(self, items: Sequence[RenderItem]):
        with self._lock:
            self._items = tuple(items)

    def clear(self):
        self.set_items(())

    def items(self) -> Tuple[RenderItem, ...]:
        with self._lock:
            return self._items

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
    def draw(self):
        items = self.items()
        if not items:
            return
        for item in draw_order(items):
            if item.kind == "surface":
                self._draw_surface(item)
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
