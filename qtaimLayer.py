"""QTAIM critical points and gradient paths drawn in the 3D Molecular Viewer's window.

``QTAIMSceneLayer`` holds an immutable ``QTAIMScene`` (Angstrom) under its own
lock and draws it from the render thread through ``MolecularViewerUI.
render_scene_extras``, inside the molecule's push-matrix. It imports OpenGL
but neither Toga, GLFW nor molecularViewer, so molecularViewer can import it.

Colours follow Multiwfn's 3D topology GUI (manual 4.2.1): magenta (3,-3),
orange (3,-1), yellow (3,+1), green (3,+3) spheres, brown bond paths and
cyan CP-index labels (dark teal on a light background).
"""

from __future__ import annotations

import threading
from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
from OpenGL.GL import *  # noqa: F401,F403 - same style as molecularViewer

from qtaimTopology import BOHR_RADIUS_ANGSTROM, CP_TYPES

CP_COLOURS = {"(3,-3)": (1.0, 0.0, 1.0), "(3,-1)": (1.0, 0.55, 0.0),
              "(3,+1)": (1.0, 1.0, 0.0), "(3,+3)": (0.0, 0.85, 0.0)}
BOND_PATH_COLOUR = (0.55, 0.33, 0.12)
RING_PATH_COLOUR = (0.2, 0.7, 0.4)
LABEL_COLOUR = (0.0, 1.0, 1.0)              # cyan, on a dark background (Multiwfn's)
LABEL_COLOUR_ON_LIGHT = (0.0, 0.42, 0.48)    # dark teal: cyan is unreadable on white
HIGHLIGHT_COLOUR = (1.0, 1.0, 1.0)
CP_RADIUS = 0.10          # Angstrom
HIGHLIGHT_RADIUS = 0.18   # Angstrom
LABEL_OFFSET = 0.15       # Angstrom, along each axis
LABEL_SCREEN_OFFSET = 0.2  # Angstrom, to the screen lower right (atom numbers sit upper right)
ON_TOP_DEPTH_RANGE = 0.01  # depth range of CPs painted over flat atoms: in front, still self-sorted
PATH_LINE_WIDTH = 2.5


def cp_label_colour(background) -> Tuple[float, float, float]:
    """CP-number colour that stays readable on ``background`` (r, g, b in [0, 1]).

    Same luminance rule as the Molecular Design Canvas's palette (is_dark).
    """
    r, g, b = (float(c) for c in background)
    return LABEL_COLOUR if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.5 else LABEL_COLOUR_ON_LIGHT


@dataclass(frozen=True)
class SceneCP:
    index: int
    kind: str
    position: Tuple[float, float, float]   # Angstrom


@dataclass(frozen=True)
class QTAIMScene:
    cps: Tuple[SceneCP, ...]
    bond_paths: Tuple[np.ndarray, ...]     # (m, 3) Angstrom, float32
    ring_paths: Tuple[np.ndarray, ...]
    centre: Tuple[float, float, float] = (0.0, 0.0, 0.0)


def scene_from_results(results, centre=(0.0, 0.0, 0.0)) -> QTAIMScene:
    """The drawable scene of a qtaimAnalysis.QTAIMResults (bohr -> Angstrom)."""
    a = BOHR_RADIUS_ANGSTROM
    cps = tuple(SceneCP(c.index, c.kind, tuple(float(x) for x in c.position * a)) for c in results.critical_points)
    bonds = tuple(np.ascontiguousarray(p.points * a, dtype=np.float32) for p in results.bond_paths)
    rings = tuple(np.ascontiguousarray(p.points * a, dtype=np.float32) for p in results.ring_paths)
    return QTAIMScene(cps, bonds, rings, tuple(float(x) for x in centre))


class QTAIMSceneLayer:
    """CP spheres, bond/ring paths and CP labels, drawable in any MolecularViewer."""

    def __init__(self):
        self._lock = threading.Lock()
        self._scene: Optional[QTAIMScene] = None
        self.show_types: Dict[str, bool] = {t: True for t in CP_TYPES}
        self.show_bond_paths = True
        self.show_ring_paths = True
        self.show_labels = True
        self.highlight: Optional[int] = None

    # --- state (Toga thread) -------------------------------------------
    def set_scene(self, scene: Optional[QTAIMScene]):
        with self._lock:
            self._scene = scene

    def clear(self):
        self.set_scene(None)
        self.highlight = None

    def get(self) -> Optional[QTAIMScene]:
        with self._lock:
            return self._scene

    def set_display_options(self, show_types: Optional[Dict[str, bool]] = None, bond_paths: Optional[bool] = None,
                            ring_paths: Optional[bool] = None, labels: Optional[bool] = None, highlight=...):
        """None leaves an option unchanged; ``highlight`` is a CP index, or None for no highlight."""
        if show_types is not None:
            unknown = set(show_types) - set(CP_TYPES)
            if unknown:
                raise ValueError(f"Unknown critical point type(s): {', '.join(sorted(unknown))}.")
            self.show_types = {**self.show_types, **{k: bool(v) for k, v in show_types.items()}}
        if bond_paths is not None:
            self.show_bond_paths = bool(bond_paths)
        if ring_paths is not None:
            self.show_ring_paths = bool(ring_paths)
        if labels is not None:
            self.show_labels = bool(labels)
        if highlight is not ...:
            self.highlight = highlight

    def visible_cps(self) -> Tuple[SceneCP, ...]:
        s = self.get()
        if s is None:
            return ()
        return tuple(c for c in s.cps if self.show_types.get(c.kind, True))

    def radius(self) -> float:
        """Largest distance of a CP or path point from the scene centre (Angstrom)."""
        s = self.get()
        if s is None:
            return 0.0
        points = [np.array([c.position for c in s.cps], dtype=float).reshape(-1, 3)]
        points += [np.asarray(p, dtype=float) for p in s.bond_paths + s.ring_paths]
        allp = np.vstack(points)
        return float(np.linalg.norm(allp - np.array(s.centre), axis=1).max()) if len(allp) else 0.0

    # --- drawing (render thread) ----------------------------------------
    def draw_scene(self, viewer):
        """Draw in the molecule's frame (call inside its push-matrix)."""
        s = self.get()
        if s is None:
            return
        paths = []
        if self.show_bond_paths:
            paths += [(p, BOND_PATH_COLOUR) for p in s.bond_paths]
        if self.show_ring_paths:
            paths += [(p, RING_PATH_COLOUR) for p in s.ring_paths]
        if paths:
            self._draw_paths(paths)
        # With flat atom discs (the viewer's Design style) a bond point of an X-H
        # bond lies inside the H disc, so (3,-1), (3,+1) and (3,+3) points and the
        # CP numbers are painted over the discs; a (3,-3) stays under its own atom.
        # The spheres keep the depth test (switching it off would paint their back
        # halves over their fronts) and are squeezed into the front of the range.
        flat = bool(getattr(viewer, "atoms_are_flat", lambda: False)())
        cps = self.visible_cps()
        for c in cps:
            highlighted = c.index == self.highlight
            colour = HIGHLIGHT_COLOUR if highlighted else CP_COLOURS.get(c.kind, (0.7, 0.7, 0.7))
            on_top = flat and c.kind != CP_TYPES[0]
            if on_top:
                glDepthRange(0.0, ON_TOP_DEPTH_RANGE)
            glPushMatrix()
            glTranslatef(*c.position)
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, (*colour, 1.0))
            viewer.draw_sphere(HIGHLIGHT_RADIUS if highlighted else CP_RADIUS)
            glPopMatrix()
            if on_top:
                glDepthRange(0.0, 1.0)
        # On flat discs a (3,-3)'s number would repeat its atom's number (both follow
        # atom order), right next to it, so only the other points are numbered there.
        labelled = [c for c in cps if not (flat and c.kind == CP_TYPES[0])]
        if self.show_labels and labelled:
            colour = cp_label_colour(getattr(viewer, "background_color", (0.0, 0.0, 0.0)))
            offset = np.full(3, LABEL_OFFSET)
            axes = getattr(viewer, "current_screen_axes", None)
            if axes is not None:
                right, up, _ = axes()
                offset = (np.asarray(right) - np.asarray(up)) * LABEL_SCREEN_OFFSET
            if flat:
                glDisable(GL_DEPTH_TEST)
            try:
                for c in labelled:
                    viewer._draw_text_3d(tuple(float(x) for x in np.asarray(c.position) + offset), str(c.index),
                                         color=colour)
            finally:
                if flat:
                    glEnable(GL_DEPTH_TEST)

    @staticmethod
    def _draw_paths(paths: Sequence[Tuple[np.ndarray, Tuple[float, float, float]]]):
        glDisable(GL_LIGHTING)
        glLineWidth(PATH_LINE_WIDTH)
        try:
            for points, colour in paths:
                glColor3f(*colour)
                glBegin(GL_LINE_STRIP)
                for x, y, z in points:
                    glVertex3f(float(x), float(y), float(z))
                glEnd()
        finally:
            glLineWidth(1.0)
            glEnable(GL_LIGHTING)
