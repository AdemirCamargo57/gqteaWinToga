"""3D window for an ESP-mapped molecular surface, GaussView 6 style.

``ESPSurfaceViewer`` subclasses the project's GLFW/OpenGL MolecularViewer,
so camera, rotation, zoom, atoms, bonds and click-to-label all come from
there; it adds the colour-mapped surface (``render_scene_extras``) and the
colour bar (``draw_screen_overlay``). Both are drawn by an
``espSurfaceLayer.ESPSurfaceLayer``, which the 3D Molecular Viewer's ESP tab
shares; this window is what the Fukui tool opens.

The colour maths, ``SurfaceMarker`` and ``SurfaceDrawData`` live in
espSurfaceLayer and are re-exported here, so ``espViewer.hue_colours`` and
``from espViewer import SurfaceMarker`` keep working.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import glfw
import numpy as np

from espSurfaceLayer import (  # noqa: F401 - re-exported for existing importers
    COLOUR_UNITS,
    HUE_SPAN_DEGREES,
    MAXIMUM_MARKER_COLOUR,
    MINIMUM_MARKER_COLOUR,
    UNIT_AU,
    UNIT_KCAL,
    ESPSurfaceLayer,
    SurfaceDrawData,
    SurfaceMarker,
    colour_bar_text,
    contrasting_text_colour,
    esp_colours,
    gaussview_label,
    hue_colours,
    symmetric_range,
    unit_factor,
)
from molecularViewer import MolecularViewer

# MolecularViewer.atom_display_style values: lines, ball-and-stick ("CPK"), space-filling ("vdW").
ATOM_DISPLAY_STYLES = ("Line style", "CPK style", "vdW style")


class ESPSurfaceViewer(MolecularViewer):
    """MolecularViewer that also draws one ESP-mapped surface and its colour bar."""

    def __init__(self):
        super().__init__()
        self.window_title = "ESP Surface Viewer"
        self.esp_layer = ESPSurfaceLayer()

    # The display options live on the layer; these keep the old attributes.
    @property
    def colour_unit(self) -> str:
        return self.esp_layer.colour_unit

    @colour_unit.setter
    def colour_unit(self, value: str):
        self.esp_layer.colour_unit = value

    @property
    def surface_opacity(self) -> float:
        return self.esp_layer.surface_opacity

    @surface_opacity.setter
    def surface_opacity(self, value: float):
        self.esp_layer.surface_opacity = value

    @property
    def show_extrema(self) -> bool:
        return self.esp_layer.show_extrema

    @show_extrema.setter
    def show_extrema(self, value: bool):
        self.esp_layer.show_extrema = value

    # --- state (Toga thread) -------------------------------------------
    def set_molecule(self, symbols: Sequence[str], coordinates_angstrom: np.ndarray):
        frame = [(s, tuple(float(c) for c in xyz)) for s, xyz in zip(symbols, coordinates_angstrom)]
        self.load_designed_molecule(frame)

    def set_surface(self, vertices_angstrom, triangles, normals, values_au, isovalue,
                    markers: Sequence[SurfaceMarker] = (), half_range: Optional[float] = None,
                    caption: Optional[str] = None):
        self.esp_layer.set_surface(vertices_angstrom, triangles, normals, values_au, isovalue, markers,
                                   half_range, caption, centre=self.get_system_center())
        self.invalidate_scene_cache()
        self.request_projection_update()

    def clear_surface(self):
        self.esp_layer.clear()
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
        self.esp_layer.set_display_options(unit, opacity, show_extrema)

    def get_surface(self) -> Optional[SurfaceDrawData]:
        return self.esp_layer.get()

    def colour_bar_text(self) -> Optional[Tuple[str, str, str]]:
        return self.esp_layer.colour_bar_text()

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
        self.esp_layer.draw_scene(self)

    def draw_screen_overlay(self):
        self.esp_layer.draw_overlay(self)
