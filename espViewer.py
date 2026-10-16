"""3D window for an ESP-mapped molecular surface, GaussView 6 style.

``ESPSurfaceViewer`` subclasses the project's GLFW/OpenGL MolecularViewer,
so camera, rotation, zoom, atoms, bonds and click-to-label all come from
there; it adds the colour-mapped surface (``render_scene_extras``) and the
colour bar (``draw_screen_overlay``).

GaussView conventions reproduced (from a GaussView 6.1 screenshot of the
ethylene ESP at rho = 0.0004): horizontal colour bar across the top, a
symmetric range +/-R with the min label at the left end and the max label at
the right end written like ``-2.072e-2``, colours a linear hue sweep red ->
orange -> yellow -> green (V = 0) -> cyan -> blue, opaque surface with a
specular highlight on a black background (the default; the tool can switch
the canvas to white or grey, and the bar's text then turns black so it stays
readable). Unlike GaussView, the bar also prints the unit and the density
isovalue that defines the surface.

The colour maths are module-level pure functions, tested headless.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import glfw
import numpy as np
from OpenGL.GL import *  # noqa: F401,F403 - same style as molecularViewer

from espIntegrals import HARTREE_TO_KCAL_MOL
from molecularViewer import MolecularViewer

UNIT_AU = "a.u."
UNIT_KCAL = "kcal/mol"
COLOUR_UNITS = (UNIT_AU, UNIT_KCAL)
HUE_SPAN_DEGREES = 240.0
# MolecularViewer.atom_display_style values: lines, ball-and-stick ("CPK"), space-filling ("vdW").
ATOM_DISPLAY_STYLES = ("Line style", "CPK style", "vdW style")

MINIMUM_MARKER_COLOUR = (0.0, 1.0, 1.0)
MAXIMUM_MARKER_COLOUR = (1.0, 0.55, 0.0)


def hue_colours(fraction: np.ndarray) -> np.ndarray:
    """RGB (N, 3) for fractions in [0, 1]: 0 = red, 0.5 = green, 1 = blue.

    HSV with full saturation and value, hue = 240 degrees x fraction.
    """
    h = np.clip(np.asarray(fraction, dtype=float), 0.0, 1.0) * HUE_SPAN_DEGREES / 60.0
    rgb = np.zeros(h.shape + (3,))
    seg = np.minimum(h.astype(int), 3)
    x = h - seg
    rgb[seg == 0] = np.column_stack([np.ones_like(x), x, np.zeros_like(x)])[seg == 0]
    rgb[seg == 1] = np.column_stack([1 - x, np.ones_like(x), np.zeros_like(x)])[seg == 1]
    rgb[seg == 2] = np.column_stack([np.zeros_like(x), np.ones_like(x), x])[seg == 2]
    rgb[seg == 3] = np.column_stack([np.zeros_like(x), 1 - x, np.ones_like(x)])[seg == 3]
    return rgb


def symmetric_range(values: np.ndarray) -> float:
    """GaussView's default half-range: max |V| over the surface."""
    r = float(np.max(np.abs(values))) if np.size(values) else 0.0
    return r if r > 0 else 1.0


def esp_colours(values: np.ndarray, half_range: float) -> np.ndarray:
    """RGB for ESP values over the symmetric range [-R, +R] (clipped)."""
    return hue_colours((np.asarray(values, dtype=float) + half_range) / (2.0 * half_range))


def gaussview_label(value: float) -> str:
    """'-2.072e-2' style: 4 significant digits, exponent without padding."""
    mantissa, exponent = f"{value:.3e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def unit_factor(unit: str) -> float:
    if unit == UNIT_AU:
        return 1.0
    if unit == UNIT_KCAL:
        return HARTREE_TO_KCAL_MOL
    raise ValueError(f"Unknown colour-bar unit '{unit}'.")


def colour_bar_text(half_range_au: float, unit: str, isovalue: float) -> Tuple[str, str, str]:
    """(left label, right label, caption) for the colour bar."""
    factor = unit_factor(unit)
    caption = f"ESP ({unit})   mapped on electron density isosurface  rho = {isovalue:.6f} a.u."
    return gaussview_label(-half_range_au * factor), gaussview_label(half_range_au * factor), caption


def contrasting_text_colour(background) -> Tuple[float, float, float]:
    """Black text on a light background, white on a dark one (Rec. 709 luma)."""
    r, g, b = (float(c) for c in background)
    return (0.0, 0.0, 0.0) if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.5 else (1.0, 1.0, 1.0)


@dataclass
class SurfaceMarker:
    position: np.ndarray  # Angstrom
    kind: str             # "minimum" | "maximum"
    scope: str            # "global" | "local"
    label: str


@dataclass
class SurfaceDrawData:
    vertices: np.ndarray   # (V, 3) float32, Angstrom
    normals: np.ndarray    # (V, 3) float32
    indices: np.ndarray    # (3F,) uint32
    values: np.ndarray     # (V,) ESP in a.u.
    rgba: np.ndarray       # (V, 4) float32
    half_range: float      # a.u.
    isovalue: float
    markers: List[SurfaceMarker]
    radius: float          # max distance of a vertex from the atoms' centre
    caption: Optional[str] = None  # replaces the ESP caption (other mapped properties)


class ESPSurfaceViewer(MolecularViewer):
    """MolecularViewer that also draws one ESP-mapped surface and its colour bar."""

    def __init__(self):
        super().__init__()
        self.window_title = "ESP Surface Viewer"
        self.colour_unit = UNIT_AU
        self.surface_opacity = 1.0
        self.show_extrema = True
        self._surface: Optional[SurfaceDrawData] = None

    # --- state (Toga thread) -------------------------------------------
    def set_molecule(self, symbols: Sequence[str], coordinates_angstrom: np.ndarray):
        frame = [(s, tuple(float(c) for c in xyz)) for s, xyz in zip(symbols, coordinates_angstrom)]
        self.load_designed_molecule(frame)

    def set_surface(self, vertices_angstrom, triangles, normals, values_au, isovalue,
                    markers: Sequence[SurfaceMarker] = (), half_range: Optional[float] = None,
                    caption: Optional[str] = None):
        values = np.asarray(values_au, dtype=float)
        r = half_range if half_range else symmetric_range(values)
        verts = np.ascontiguousarray(vertices_angstrom, dtype=np.float32)
        centre = np.array(self.get_system_center(), dtype=float)
        data = SurfaceDrawData(
            vertices=verts,
            normals=np.ascontiguousarray(normals, dtype=np.float32),
            indices=np.ascontiguousarray(np.asarray(triangles).ravel(), dtype=np.uint32),
            values=values,
            rgba=self._rgba(values, r),
            half_range=r,
            isovalue=float(isovalue),
            markers=list(markers),
            radius=float(np.linalg.norm(verts - centre, axis=1).max()) if len(verts) else 0.0,
            caption=caption,
        )
        with self._state_lock:
            self._surface = data
        self.invalidate_scene_cache()
        self.request_projection_update()

    def clear_surface(self):
        with self._state_lock:
            self._surface = None
        self.invalidate_scene_cache()

    def set_atom_view(self, style: str, show_numbers: bool):
        """How the atoms are drawn (one of ATOM_DISPLAY_STYLES) and whether 1-based numbers are shown."""
        if style not in ATOM_DISPLAY_STYLES:
            raise ValueError(f"Unknown atom style '{style}'; expected one of {ATOM_DISPLAY_STYLES}.")
        self.atom_display_style = style
        self.show_atom_numbers = bool(show_numbers)
        self.invalidate_scene_cache()
        self.request_projection_update()

    def set_display_options(self, unit: str = None, opacity: float = None, show_extrema: bool = None):
        if unit is not None:
            unit_factor(unit)  # validates
            self.colour_unit = unit
        if show_extrema is not None:
            self.show_extrema = bool(show_extrema)
        if opacity is not None:
            self.surface_opacity = float(min(1.0, max(0.05, opacity)))
            with self._state_lock:
                if self._surface is not None:
                    self._surface.rgba = self._rgba(self._surface.values, self._surface.half_range)

    def get_surface(self) -> Optional[SurfaceDrawData]:
        with self._state_lock:
            return self._surface

    def _rgba(self, values, half_range) -> np.ndarray:
        rgb = esp_colours(values, half_range)
        alpha = np.full((len(rgb), 1), self.surface_opacity)
        return np.ascontiguousarray(np.hstack([rgb, alpha]), dtype=np.float32)

    def colour_bar_text(self) -> Optional[Tuple[str, str, str]]:
        s = self.get_surface()
        if s is None:
            return None
        left, right, caption = colour_bar_text(s.half_range, self.colour_unit, s.isovalue)
        return left, right, s.caption or caption

    def get_scene_radius(self, frame_data=None) -> float:
        base = super().get_scene_radius(frame_data)
        s = self.get_surface()
        return max(base, s.radius + 0.3) if s is not None else base

    def request_close(self):
        """Ask GLFW to close the window (glfwSetWindowShouldClose is thread-safe)."""
        if self.window is not None:
            try:
                glfw.set_window_should_close(self.window, True)
            except Exception:
                pass

    # --- drawing (render thread) ----------------------------------------
    def render_scene_extras(self):
        s = self.get_surface()
        if s is None:
            return
        self._draw_surface(s)
        if self.show_extrema:
            self._draw_markers(s)

    def _draw_surface(self, s: SurfaceDrawData):
        translucent = self.surface_opacity < 0.999
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
            glVertexPointer(3, GL_FLOAT, 0, s.vertices)
            glNormalPointer(GL_FLOAT, 0, s.normals)
            glColorPointer(4, GL_FLOAT, 0, s.rgba)
            glDrawElements(GL_TRIANGLES, int(s.indices.size), GL_UNSIGNED_INT, s.indices)
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

    def _draw_markers(self, s: SurfaceDrawData):
        for m in s.markers:
            colour = MINIMUM_MARKER_COLOUR if m.kind == "minimum" else MAXIMUM_MARKER_COLOUR
            radius = 0.18 if m.scope == "global" else 0.12
            glPushMatrix()
            glTranslatef(*[float(c) for c in m.position])
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, (*colour, 1.0))
            self.draw_sphere(radius)
            glPopMatrix()
            label_pos = tuple(float(c) + 0.2 for c in m.position)
            self._draw_text_3d(label_pos, m.label, color=colour)

    def draw_screen_overlay(self):
        text = self.colour_bar_text()
        if text is None:
            return
        left, right, caption = text
        w, h = self.display
        x0, x1 = 100.0, max(110.0, w - 100.0)
        y0, y1 = h - 34.0, h - 16.0
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, w, 0, h, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        glDisable(GL_LIGHTING)
        glDisable(GL_DEPTH_TEST)
        try:
            steps = 96
            colours = hue_colours(np.linspace(0.0, 1.0, steps + 1))
            glBegin(GL_QUAD_STRIP)
            for k in range(steps + 1):
                x = x0 + (x1 - x0) * k / steps
                glColor3f(*colours[k])
                glVertex2f(x, y0)
                glVertex2f(x, y1)
            glEnd()
            ink = contrasting_text_colour(self.background_color)
            self._draw_text_2d(6.0, y0 + 2.0, left, color=ink)
            self._draw_text_2d(x1 + 6.0, y0 + 2.0, right, color=ink)
            self._draw_text_2d(x0, y0 - 22.0, caption, color=ink)
        finally:
            glEnable(GL_DEPTH_TEST)
            glEnable(GL_LIGHTING)
            glPopMatrix()
            glMatrixMode(GL_PROJECTION)
            glPopMatrix()
            glMatrixMode(GL_MODELVIEW)
