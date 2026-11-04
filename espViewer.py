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

import threading
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

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
from selectionViewer import SelectionViewer
from viewerSelection import submesh
from volumeData import same_structure

# MolecularViewer.atom_display_style values: lines, ball-and-stick ("CPK"), space-filling ("vdW").
ATOM_DISPLAY_STYLES = ("Line style", "CPK style", "vdW style")


@dataclass(frozen=True)
class _SurfaceSource:
    """The whole surface as given to set_surface; the drawn part is derived from it."""
    vertices: np.ndarray
    triangles: np.ndarray
    normals: np.ndarray
    values: np.ndarray
    isovalue: float
    markers: List[SurfaceMarker]
    half_range: float
    caption: Optional[str]
    vertex_atoms: Optional[np.ndarray]   # None: never filtered
    marker_atoms: Optional[np.ndarray]


class ESPSurfaceViewer(SelectionViewer):
    """MolecularViewer that also draws one ESP-mapped surface and its colour bar.

    ``selection_listener`` (from SelectionViewer) is set by the Fukui window.
    """

    def __init__(self):
        # Before super(): the base constructor may already touch the selection hooks.
        self._source_lock = threading.Lock()
        self._surface_source: Optional[_SurfaceSource] = None
        self._hidden_by_selection = False
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
        """Show this molecule; the structure already shown is left alone, selection included."""
        frame = [(s, tuple(float(c) for c in xyz)) for s, xyz in zip(symbols, coordinates_angstrom)]
        with self._state_lock:
            shown = list(self.frames[self.current_frame]) if self.frames else None
        if shown is not None and same_structure(shown, frame):
            return
        self.load_designed_molecule(frame)

    def set_surface(self, vertices_angstrom, triangles, normals, values_au, isovalue,
                    markers: Sequence[SurfaceMarker] = (), half_range: Optional[float] = None,
                    caption: Optional[str] = None, vertex_atoms=None, marker_atoms=None):
        """The surface to draw; with ``vertex_atoms`` (0-based owner atom of each vertex, and
        ``marker_atoms`` of each marker) only the selected atoms' part is drawn while the
        viewer shows only its selection."""
        values = np.asarray(values_au, dtype=float)
        markers = list(markers)
        if vertex_atoms is not None:
            vertex_atoms = np.asarray(vertex_atoms, dtype=np.intp)
            if vertex_atoms.shape != (len(values),):
                raise ValueError(f"vertex_atoms has {vertex_atoms.size} entries for {len(values)} vertices.")
            marker_atoms = np.asarray([] if marker_atoms is None else marker_atoms, dtype=np.intp)
            if marker_atoms.shape != (len(markers),):
                raise ValueError(f"marker_atoms has {marker_atoms.size} entries for {len(markers)} markers.")
        source = _SurfaceSource(
            np.asarray(vertices_angstrom, dtype=float), np.asarray(triangles).reshape(-1, 3),
            np.asarray(normals, dtype=float), values, float(isovalue), markers,
            # Resolved on the whole surface, so a filtered view keeps the full colour scale.
            float(half_range) if half_range else symmetric_range(values),
            caption, vertex_atoms, marker_atoms)
        with self._source_lock:
            self._surface_source = source
        self._show_surface_part()

    def clear_surface(self):
        with self._source_lock:
            self._surface_source = None
            self._hidden_by_selection = False
        self.esp_layer.clear()
        self.invalidate_scene_cache()

    # --- show only selected ----------------------------------------------
    def filtering_by_selection(self) -> bool:
        """True while only the selected atoms (and their part of the surface) are drawn."""
        return self._shown_selection() is not None

    def surface_hidden_by_selection(self) -> bool:
        """True when a surface exists but none of it belongs to the shown atoms."""
        with self._source_lock:
            return self._hidden_by_selection

    def _refocus_camera(self):
        # The base viewer's hook for "the set of shown atoms changed" (any
        # selection change while showing only it, and the switch itself),
        # possibly on the render thread after a Ctrl/Shift-click.
        self._show_surface_part()
        super()._refocus_camera()

    def _show_surface_part(self):
        """Give the layer the whole surface, or only the shown atoms' part of it."""
        with self._source_lock:
            source = self._surface_source
            if source is None:
                return
            shown = self._shown_selection() if source.vertex_atoms is not None else None
            vertices, triangles, normals, values, markers = (
                source.vertices, source.triangles, source.normals, source.values, source.markers)
            if shown is not None:
                selected = np.fromiter(shown, dtype=np.intp)
                used, triangles = submesh(source.triangles, np.isin(source.vertex_atoms, selected))
                vertices, normals, values = source.vertices[used], source.normals[used], source.values[used]
                markers = [m for m, atom in zip(source.markers, source.marker_atoms) if atom in shown]
            self._hidden_by_selection = shown is not None and not len(triangles)
            if self._hidden_by_selection:
                self.esp_layer.clear()
            else:
                focus = self.camera_focus_frame()
                self.esp_layer.set_surface(vertices, triangles, normals, values, source.isovalue, markers,
                                           source.half_range, source.caption,
                                           centre=self.get_system_center(focus))
        self.invalidate_scene_cache()
        self.request_projection_update()

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

    # --- drawing (render thread) ----------------------------------------
    def render_scene_extras(self):
        self.esp_layer.draw_scene(self)

    def draw_screen_overlay(self):
        self.esp_layer.draw_overlay(self)
