"""Live view of a solvation box while it is built and pre-optimized.

Three layers:

* :func:`build_draw_list` -- pure numpy: projects a :class:`ConstructionSnapshot`
  (cell + molecules tagged by kind) through an orthographic
  :class:`BoxCamera` and returns paint primitives back to front.
* :func:`render_image` -- rasterises a draw list into one Pillow image
  (2x supersampled, then downscaled), headless.
* :class:`ConstructionView` -- a ``toga.ImageView`` showing that image, with
  the rotate / zoom / reset actions. Toga is imported lazily inside it, so
  importing this module (as the tests do) never loads the GUI stack.

**Why a bitmap and not ``toga.Canvas``.** On WinForms a Toga canvas keeps its
drawing as Python objects and replays every one of them through pythonnet on
each repaint: ~120 ms for a 2300-atom box when idle, and *seconds* while a
packing worker thread competes for the GIL (each replayed call waits up to the
5 ms switch interval) -- the window froze for up to 9 s in a live test. Pillow
draws the same frame in C in a few milliseconds, and the ImageView repaints a
native bitmap with no Python at all.

Keeping the window responsive is also bounded by :meth:`ConstructionView.submit`:
it only stores the newest snapshot and renders it at most once per
``min_interval`` seconds (a timer flushes the last one); above ``max_atoms``
each solvent molecule is drawn as a single dot at its centroid.
"""
import time
from dataclasses import dataclass, field as dataclass_field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from molecularCamera import ViewCamera

CANVAS_BACKGROUND = "#fbfbfb"
BOX_EDGE_COLOR = "#555555"
CAPTION_COLOR = "#222222"
KIND_COLORS = {
    "solute": "#e6550d",
    "solvent": "#3182bd",
    "solvent_A": "#3182bd",
    "solvent_B": "#31a354",
}
DEFAULT_KIND_COLOR = "#7f7f7f"
LEGEND_LABELS = {
    "solute": "solute",
    "solvent": "solvent",
    "solvent_A": "solvent A",
    "solvent_B": "solvent B",
}

# Disc radii in Angstrom (scaled by the camera zoom): small enough that a
# packed box still shows its structure. The solute is drawn larger so it
# stands out among thousands of solvent atoms.
ATOM_DRAW_RADIUS = 0.30
HYDROGEN_DRAW_RADIUS = 0.20
MOLECULE_DOT_RADIUS = 0.55
SOLUTE_SCALE = 1.5
MIN_DISC_PIXELS = 1.0
DEPTH_FADE = 0.55
FADE_LEVELS = 6

DEFAULT_MAX_ATOMS = 4000
DEFAULT_MIN_INTERVAL = 0.25
# The next repaint waits at least this many times as long as the last one
# took, so painting can never use more than ~1/(1 + factor) of the UI thread.
PAINT_BUDGET_FACTOR = 3.0

# The default view: Angstrom y points up but screen y points down, so the
# frame is turned 180 degrees about x (not mirrored), then tilted so three
# faces of the box are visible.
_FLIP = np.diag([1.0, -1.0, -1.0])


def _rotation_x(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def _rotation_y(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


DEFAULT_VIEW = _rotation_x(np.radians(-20.0)) @ _rotation_y(np.radians(-30.0)) @ _FLIP


class BoxCamera(ViewCamera):
    """ViewCamera with a zoom range wide enough for 10-100 A boxes (px per A)."""

    MIN_ZOOM = 0.05
    MAX_ZOOM = 200.0


def blend_color(color: str, background: str, fraction: float) -> str:
    fraction = min(1.0, max(0.0, float(fraction)))
    mixed = []
    for start in (1, 3, 5):
        front = int(color[start:start + 2], 16)
        back = int(background[start:start + 2], 16)
        mixed.append(int(round(front + (back - front) * fraction)))
    return "#" + "".join(f"{channel:02x}" for channel in mixed)


def box_corners(lengths) -> np.ndarray:
    a, b, c = (float(v) for v in lengths)
    return np.array([[x, y, z] for x in (0.0, a) for y in (0.0, b) for z in (0.0, c)])


# Index pairs into box_corners() that differ in exactly one coordinate.
BOX_EDGES = [(i, j) for i in range(8) for j in range(i + 1, 8) if bin(i ^ j).count("1") == 1]


@dataclass
class ConstructionSnapshot:
    stage: str
    lengths: Tuple[float, float, float]
    # (kind, elements, (n, 3) coordinates in Angstrom)
    molecules: List[Tuple[str, Sequence[str], np.ndarray]] = dataclass_field(default_factory=list)
    caption: str = ""


def fit_camera(camera: BoxCamera, lengths, width: float, height: float, margin: float = 0.85):
    """Default oblique view with the whole box centred and filling the canvas."""
    camera.reset()
    camera.rotation = DEFAULT_VIEW.copy()
    lengths = np.asarray(lengths, dtype=float)
    camera.pivot = 0.5 * lengths
    camera.anchor = np.array([0.5 * width, 0.5 * height + 8.0])
    projected = (box_corners(lengths) - camera.pivot) @ camera.rotation.T
    span_x = float(np.ptp(projected[:, 0])) or 1.0
    span_y = float(np.ptp(projected[:, 1])) or 1.0
    zoom = margin * min(width / span_x, (height - 30.0) / span_y)
    camera.zoom = float(np.clip(zoom, camera.MIN_ZOOM, camera.MAX_ZOOM))


def _disc_radius(element: str, kind: str, zoom: float) -> float:
    base = HYDROGEN_DRAW_RADIUS if element == "H" else ATOM_DRAW_RADIUS
    if kind == "solute":
        base *= SOLUTE_SCALE
    return max(MIN_DISC_PIXELS, base * zoom)


def build_draw_list(snapshot: ConstructionSnapshot, camera: BoxCamera,
                    max_atoms: int = DEFAULT_MAX_ATOMS) -> List[tuple]:
    """Paint primitives, back to front.

    ``("lines", color, width, [(x1, y1, x2, y2), ...])`` -- the 12 box edges;
    ``("discs", color, [(x, y, r), ...])`` -- runs of same-coloured discs,
    in depth order (farthest first);
    ``("text", color, x, y, text)`` -- caption and legend.
    """
    primitives: List[tuple] = []
    corners = camera.project(box_corners(snapshot.lengths))
    segments = [(corners[i, 0], corners[i, 1], corners[j, 0], corners[j, 1]) for i, j in BOX_EDGES]
    primitives.append(("lines", BOX_EDGE_COLOR, 1.2, segments))

    total_atoms = sum(len(m[1]) for m in snapshot.molecules)
    as_dots = total_atoms > max_atoms
    points, radii, kinds = [], [], []
    for kind, elements, coords in snapshot.molecules:
        coords = np.asarray(coords, dtype=float).reshape(-1, 3)
        if as_dots and kind != "solute":
            points.append(coords.mean(axis=0)[None, :])
            radii.append([max(MIN_DISC_PIXELS, MOLECULE_DOT_RADIUS * camera.zoom)])
            kinds.append(kind)
            continue
        points.append(coords)
        radii.append([_disc_radius(e, kind, camera.zoom) for e in elements])
        kinds.extend([kind] * len(elements))

    if points:
        projected = camera.project(np.concatenate(points))
        radii = np.concatenate([np.asarray(r, dtype=float) for r in radii])
        depth = projected[:, 2]
        span = float(np.ptp(depth))
        fade = (depth - depth.min()) / span if span > 0.0 else np.zeros_like(depth)
        level = np.minimum(FADE_LEVELS - 1, (fade * FADE_LEVELS).astype(int))
        order = np.argsort(-depth, kind="stable")  # farthest first
        run_color, run = None, []
        for index in order:
            base = KIND_COLORS.get(kinds[index], DEFAULT_KIND_COLOR)
            tint = DEPTH_FADE * level[index] / max(1, FADE_LEVELS - 1)
            color = blend_color(base, CANVAS_BACKGROUND, tint)
            if color != run_color and run:
                primitives.append(("discs", run_color, run))
                run = []
            run_color = color
            run.append((float(projected[index, 0]), float(projected[index, 1]), float(radii[index])))
        if run:
            primitives.append(("discs", run_color, run))

    caption = snapshot.stage + (f" -- {snapshot.caption}" if snapshot.caption else "")
    primitives.append(("text", CAPTION_COLOR, 8.0, 6.0, caption))
    present = []
    for kind, _, _ in snapshot.molecules:
        if kind not in present:
            present.append(kind)
    for row, kind in enumerate(present):
        label = LEGEND_LABELS.get(kind, kind) + (" (centroids)" if as_dots and kind != "solute" else "")
        primitives.append(("text", KIND_COLORS.get(kind, DEFAULT_KIND_COLOR), 8.0, 24.0 + 16.0 * row, label))
    return primitives


def _rgb(color: str) -> Tuple[int, int, int]:
    return int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)


def render_image(primitives: List[tuple], width: int, height: int, supersample: int = 2):
    """One RGB Pillow image of a draw list (``supersample`` x for smooth edges)."""
    from PIL import Image, ImageDraw, ImageFont

    k = max(1, int(supersample))
    image = Image.new("RGB", (int(width) * k, int(height) * k), _rgb(CANVAS_BACKGROUND))
    draw = ImageDraw.Draw(image)
    try:
        font = ImageFont.load_default(size=11 * k)
    except TypeError:  # Pillow < 10.1 has no sized default font
        font = ImageFont.load_default()
    for primitive in primitives:
        kind = primitive[0]
        if kind == "lines":
            _, color, line_width, segments = primitive
            for x1, y1, x2, y2 in segments:
                draw.line((x1 * k, y1 * k, x2 * k, y2 * k), fill=_rgb(color), width=max(1, int(round(line_width * k))))
        elif kind == "discs":
            _, color, discs = primitive
            fill = _rgb(color)
            for x, y, r in discs:
                draw.ellipse((k * (x - r), k * (y - r), k * (x + r), k * (y + r)), fill=fill)
        elif kind == "text":
            _, color, x, y, text = primitive
            draw.text((x * k, y * k), text, fill=_rgb(color), font=font)
    if k > 1:
        image = image.resize((int(width), int(height)), Image.LANCZOS)
    return image


ROTATION_STEP = 15.0  # degrees per rotate-button press


class ConstructionView:
    """``toga.ImageView`` showing the newest snapshot, re-rendered at a bounded rate.

    :meth:`rotate`, :meth:`zoom` and :meth:`reset_view` back the window's view
    buttons. Must be used from the Toga thread (the builders hand their
    snapshots over with ``loop.call_soon_threadsafe``).
    """

    def __init__(self, width: int = 760, height: int = 430,
                 min_interval: float = DEFAULT_MIN_INTERVAL, max_atoms: int = DEFAULT_MAX_ATOMS):
        import toga
        from toga.style import Pack

        self.width = int(width)
        self.height = int(height)
        self.min_interval = float(min_interval)
        self.max_atoms = int(max_atoms)
        self.camera = BoxCamera()
        self.snapshot: Optional[ConstructionSnapshot] = None
        self._fitted_lengths = None
        self._last_paint = 0.0
        self._timer = None
        self.paint_count = 0
        self.widget = toga.ImageView(
            image=render_image([], self.width, self.height, supersample=1),
            style=Pack(flex=1, height=self.height, background_color=CANVAS_BACKGROUND),
        )

    # -- snapshots -----------------------------------------------------
    def submit(self, snapshot: ConstructionSnapshot, force: bool = False):
        self.snapshot = snapshot
        lengths = tuple(float(v) for v in snapshot.lengths)
        if lengths != self._fitted_lengths:
            fit_camera(self.camera, lengths, self.width, self.height)
            self._fitted_lengths = lengths
        now = time.monotonic()
        interval = self._effective_interval()
        if force or now - self._last_paint >= interval:
            self._cancel_timer()
            self.paint()
        elif self._timer is None:
            import asyncio
            delay = interval - (now - self._last_paint)
            self._timer = asyncio.get_event_loop().call_later(delay, self._flush)

    def _flush(self):
        self._timer = None
        self.paint()

    def _cancel_timer(self):
        if self._timer is not None:
            self._timer.cancel()
            self._timer = None

    # -- view ----------------------------------------------------------
    def rotate(self, horizontal_degrees: float = 0.0, vertical_degrees: float = 0.0):
        """Turn the box about the screen's vertical / horizontal axis."""
        if self.snapshot is not None:
            self.camera.set_pivot(0.5 * np.asarray(self.snapshot.lengths, dtype=float))
        per_pixel = self.camera.ROTATION_PER_PIXEL
        self.camera.rotate_by_drag(-np.radians(horizontal_degrees) / per_pixel,
                                   np.radians(vertical_degrees) / per_pixel)
        self.paint()

    def zoom(self, factor: float):
        self.camera.zoom_by(factor, 0.5 * self.width, 0.5 * self.height)
        self.paint()

    def reset_view(self):
        if self.snapshot is not None:
            fit_camera(self.camera, self.snapshot.lengths, self.width, self.height)
        self.paint()

    # -- painting ------------------------------------------------------
    def render(self):
        """The current frame as a Pillow image (headless)."""
        primitives = [] if self.snapshot is None else build_draw_list(self.snapshot, self.camera, self.max_atoms)
        return render_image(primitives, self.width, self.height)

    def paint(self):
        started = time.monotonic()
        image = self.render()
        if not self._set_native_bitmap(image):
            self.widget.image = image
        self.paint_count += 1
        finished = time.monotonic()
        self._last_paint = finished
        # Never let the view take more than about a quarter of the UI thread.
        self._paint_cost = finished - started

    _paint_cost = 0.0
    _native_bitmap = None

    def _effective_interval(self) -> float:
        return max(self.min_interval, PAINT_BUDGET_FACTOR * self._paint_cost)

    def _set_native_bitmap(self, image) -> bool:
        """WinForms fast path: copy the pixels straight into a .NET Bitmap.

        ``toga.Image(pil)`` encodes a PNG and pythonnet then converts it to a
        .NET byte array element by element -- ~170 ms per frame for a
        770x430 view, measured. Locking a ``System.Drawing.Bitmap`` and
        ``memmove``-ing the BGR rows into it takes about a millisecond. Any
        failure (another backend, no pythonnet, a stub widget in the tests)
        returns False and the caller uses Toga's own ``image`` setter.
        """
        native = getattr(getattr(self.widget, "_impl", None), "native", None)
        if native is None or not hasattr(native, "Image"):
            return False
        try:
            import ctypes
            from System.Drawing import Bitmap, Rectangle
            from System.Drawing.Imaging import ImageLockMode, PixelFormat

            rgb = image.convert("RGB")
            width, height = rgb.size
            bitmap = Bitmap(width, height, PixelFormat.Format24bppRgb)
            locked = bitmap.LockBits(Rectangle(0, 0, width, height), ImageLockMode.WriteOnly,
                                     PixelFormat.Format24bppRgb)
            try:
                stride = int(locked.Stride)
                # Rows padded to the bitmap's stride, in the BGR order it stores.
                raw = rgb.tobytes("raw", "BGR", stride, 1)
                ctypes.memmove(int(locked.Scan0.ToInt64()), raw, stride * height)
            finally:
                bitmap.UnlockBits(locked)
        except Exception:
            return False
        previous = self._native_bitmap
        native.Image = bitmap
        self._native_bitmap = bitmap
        if previous is not None:
            previous.Dispose()  # only bitmaps this view created; Toga's own image is left alone
        return True
