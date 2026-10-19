"""View camera for the 3D molecular design canvas.

The sketch stores every atom in **sketch coordinates**: canvas pixels at the
default view (x to the right, y downward, as the canvas reports clicks) plus a
depth ``z`` that grows *away* from the viewer. With y pointing down that keeps
the frame right-handed, so rotating it never mirrors a chiral structure.

:class:`ViewCamera` maps those coordinates to the screen and back with an
orthographic projection::

    v      = R @ (p - pivot)
    screen = zoom * v[:2] + anchor        depth = v[2]

The default view is the identity (``R = I``, ``zoom = 1``, ``pivot = anchor =
0``), so a fresh canvas behaves exactly like the flat 2D canvas it replaced.
The view never changes stored coordinates: rotating, panning and zooming are
not edits, push no undo entry and leave a finished optimization valid.

Pure numpy, no Toga -- tested headless in tests/test_molecularCamera.py.
"""
import numpy as np


class ViewCamera:
    """Orthographic rotate / pan / zoom view over sketch coordinates."""

    # A full-width drag across the canvas (~560 px) turns the molecule about
    # 320 degrees: fine control without needing several drags for a half turn.
    ROTATION_PER_PIXEL = 0.01
    MIN_ZOOM = 0.2
    MAX_ZOOM = 8.0

    def __init__(self):
        self.reset()

    def reset(self):
        """Back to the identity view."""
        self.rotation = np.eye(3)
        self.pivot = np.zeros(3)
        self.anchor = np.zeros(2)
        self.zoom = 1.0

    # ------------------------------------------------------------------
    # Mapping
    # ------------------------------------------------------------------
    def project(self, points) -> np.ndarray:
        """``(N, 3)`` sketch points to ``(N, 3)`` rows of ``(sx, sy, depth)``."""
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        view = (points - self.pivot) @ self.rotation.T
        screen = np.empty_like(view)
        screen[:, :2] = self.zoom * view[:, :2] + self.anchor
        screen[:, 2] = view[:, 2]
        return screen

    def unproject(self, sx: float, sy: float, depth: float) -> np.ndarray:
        """The sketch point that projects to ``(sx, sy)`` at view ``depth``."""
        view = np.array(
            [(sx - self.anchor[0]) / self.zoom, (sy - self.anchor[1]) / self.zoom, depth],
            dtype=float,
        )
        return self.rotation.T @ view + self.pivot

    def depth_of(self, point) -> float:
        return float(self.project(point)[0][2])

    # ------------------------------------------------------------------
    # Changing the view
    # ------------------------------------------------------------------
    def set_pivot(self, point):
        """Rotate about ``point`` from now on, without moving the picture.

        The anchor is moved to wherever the new pivot currently projects, so
        every point keeps its screen position; only the depth origin shifts.
        """
        point = np.asarray(point, dtype=float).reshape(3)
        self.anchor = self.project(point)[0][:2].copy()
        self.pivot = point.copy()

    def rotate_by_drag(self, dx: float, dy: float):
        """Trackball rotation: the front of the molecule follows the pointer.

        A horizontal drag turns about the screen's vertical axis and a vertical
        drag about its horizontal axis. The increment is applied in screen
        space, so the response is the same however the molecule is oriented.
        """
        rotation = self.drag_increment(dx, dy) @ self.rotation
        # Re-orthonormalize so thousands of small drags cannot accumulate
        # into a shear or a scale.
        left, _, right = np.linalg.svd(rotation)
        self.rotation = left @ right

    def drag_increment(self, dx: float, dy: float) -> np.ndarray:
        """The screen-space rotation a drag of ``(dx, dy)`` pixels stands for."""
        about_vertical = -dx * self.ROTATION_PER_PIXEL
        about_horizontal = dy * self.ROTATION_PER_PIXEL

        cos_y, sin_y = np.cos(about_vertical), np.sin(about_vertical)
        turn_y = np.array([[cos_y, 0.0, sin_y], [0.0, 1.0, 0.0], [-sin_y, 0.0, cos_y]])
        cos_x, sin_x = np.cos(about_horizontal), np.sin(about_horizontal)
        turn_x = np.array([[1.0, 0.0, 0.0], [0.0, cos_x, -sin_x], [0.0, sin_x, cos_x]])
        return turn_x @ turn_y

    def object_rotation(self, dx: float, dy: float) -> np.ndarray:
        """The sketch-space rotation that turns an object the way a drag turns the view.

        Applied about the object's own centre, it moves the object on screen
        exactly as :meth:`rotate_by_drag` would move the whole picture.
        """
        return self.rotation.T @ self.drag_increment(dx, dy) @ self.rotation

    def screen_to_sketch_delta(self, dx: float, dy: float) -> np.ndarray:
        """The sketch displacement that moves a point ``(dx, dy)`` on screen, at constant depth."""
        return self.rotation.T @ np.array([dx / self.zoom, dy / self.zoom, 0.0])

    def pan_by(self, dx: float, dy: float):
        self.anchor = self.anchor + np.array([dx, dy], dtype=float)

    def zoom_by(self, factor: float, centre_x: float, centre_y: float):
        """Scale the view about a screen point, which stays where it is."""
        new_zoom = float(np.clip(self.zoom * factor, self.MIN_ZOOM, self.MAX_ZOOM))
        effective = new_zoom / self.zoom
        centre = np.array([centre_x, centre_y], dtype=float)
        self.anchor = centre + effective * (self.anchor - centre)
        self.zoom = new_zoom
