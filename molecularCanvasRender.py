"""Fast painting for the molecular design canvas: one bitmap per frame.

A Toga canvas on WinForms replays every drawn shape through pythonnet as an
anti-aliased GDI+ path on each repaint: ~15 ms for benzene and growing with
every atom, which made dragging, rotating and panning lag behind the pointer.
Here a frame is instead described as a list of plain **primitives**, drawn by
Pillow in C (2x supersampled for smooth edges) at the control's physical
pixel size, and handed to WinForms as one ``System.Drawing.Bitmap`` that the
canvas's own Paint handler copies 1:1 -- so a repaint is a single blit,
whatever the molecule. The same technique, for the same reason, as
[boxConstructionView.py](boxConstructionView.py). (Showing the frame as the
control's ``BackgroundImage`` instead costs ~13 ms per paint, measured,
against ~1 ms for the Paint-handler blit: do not switch to it.)

Primitives, all in logical canvas units (the coordinates Toga reports clicks
in), colours as ``#rrggbb``:

    ("rect", x0, y0, x1, y1, fill)
    ("disc", cx, cy, radius, fill, outline or None, outline_width)
    ("ring", cx, cy, radius, color, width)
    ("line", x1, y1, x2, y2, color, width, dash or None)
    ("polyline", [(x, y), ...], color, width, dash or None)
    ("text", x, y, text, color, size_points, bold, baseline)   baseline: top / middle / bottom
    ("label", cx, cy, text, color, size_points, bold)          centred on (cx, cy)

molecularDesign builds the list once per frame and paints it either here or,
where no native bitmap can be made, through Toga's own canvas context, so both
paths show the same picture. No Toga import; Pillow and pythonnet are imported
lazily.
"""
from functools import lru_cache
from typing import List, Optional, Sequence, Tuple

import numpy as np

SUPERSAMPLE = 2
POINTS_TO_PIXELS = 96.0 / 72.0
FONT_FILES = ("segoeui.ttf", "arial.ttf", "DejaVuSans.ttf")
BOLD_FONT_FILES = ("segoeuib.ttf", "arialbd.ttf", "DejaVuSans-Bold.ttf")
# Pillow text anchors for the Toga baselines the canvas uses.
ANCHORS = {"top": "lt", "middle": "lm", "bottom": "lb"}


@lru_cache(maxsize=32)
def _font(pixel_size: int, bold: bool):
    from PIL import ImageFont

    for name in BOLD_FONT_FILES if bold else FONT_FILES:
        try:
            return ImageFont.truetype(name, pixel_size)
        except OSError:
            continue
    try:
        return ImageFont.load_default(size=pixel_size)
    except TypeError:  # Pillow < 10.1 has a single bitmap font
        return ImageFont.load_default()


def dash_segments(x1: float, y1: float, x2: float, y2: float,
                  pattern: Sequence[float]) -> List[Tuple[Tuple[float, float], Tuple[float, float]]]:
    """The 'on' pieces of a dashed line, as (start, end) points."""
    dx, dy = x2 - x1, y2 - y1
    length = (dx * dx + dy * dy) ** 0.5
    if length == 0.0 or not pattern or sum(pattern) <= 0.0:
        return [((x1, y1), (x2, y2))]
    ux, uy = dx / length, dy / length
    segments, position, index = [], 0.0, 0
    while position < length:
        step = pattern[index % len(pattern)]
        end = min(length, position + step)
        if index % 2 == 0:
            segments.append(((x1 + ux * position, y1 + uy * position), (x1 + ux * end, y1 + uy * end)))
        position, index = end, index + 1
    return segments


def text_box_size(text: str, size_points: float) -> Tuple[float, float]:
    """Approximate (width, height) of a one-line label, in logical units.

    Deliberately a little generous: label placement only needs to keep the
    box clear, and both renderers (Pillow, Toga) draw within it.
    """
    pixels = size_points * POINTS_TO_PIXELS
    return 0.62 * pixels * len(text) + 2.0, pixels + 1.0


# Label directions in order of preference (screen y grows downwards): upper
# right first, as chemists write atom numbers, then the other diagonals and
# the axes, then the in-between directions.
_LABEL_ANGLES = np.radians([45, 0, 90, 135, 315, 180, 270, 225,
                            22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5])
LABEL_DIRECTIONS = np.c_[np.cos(_LABEL_ANGLES), -np.sin(_LABEL_ANGLES)]
_BOND_SAMPLES = np.linspace(0.08, 0.92, 12)


def place_atom_labels(centres, radii, bonds: Sequence[Tuple[int, int]], sizes,
                      obstacles: Optional[Sequence[Tuple[float, float, float, float]]] = None,
                      gap: float = 2.0) -> np.ndarray:
    """Centres for one label beside each atom, clear of atoms, bonds and each other.

    ``centres`` are the atoms' screen positions ``(N, 2)``, ``radii`` their
    drawn disc radii, ``bonds`` 0-based index pairs, ``sizes`` each label's
    ``(width, height)``; ``obstacles`` are boxes ``(x0, y0, x1, y1)`` to keep
    clear too (other text already on the canvas). For every atom, 16
    directions at two distances are scored -- a collision with an atom disc
    or a bond weighs most, one with a label already placed next, the farther
    ring a little, and the direction order breaks ties -- and the cheapest is
    taken. A label touches its own atom's disc at most at ``gap``. Greedy and
    deterministic; vectorized over candidates and nearby obstacles only.
    """
    from scipy.spatial import cKDTree

    centres = np.asarray(centres, dtype=float).reshape(-1, 2)
    count = len(centres)
    result = np.zeros((count, 2))
    if count == 0:
        return result
    radii = np.broadcast_to(np.asarray(radii, dtype=float), (count,))
    half = np.asarray(sizes, dtype=float).reshape(-1, 2) / 2.0
    pairs = np.asarray(bonds, dtype=int).reshape(-1, 2)
    bonds_of: List[List[int]] = [[] for _ in range(count)]
    for index, (i, j) in enumerate(pairs):
        bonds_of[i].append(index)
        bonds_of[j].append(index)
    # Points along every bond: a label box containing one sits on the bond.
    samples = (centres[pairs[:, 0], None, :] * (1.0 - _BOND_SAMPLES[None, :, None])
               + centres[pairs[:, 1], None, :] * _BOND_SAMPLES[None, :, None]) if len(pairs) else None

    placed = np.asarray(obstacles if obstacles else np.zeros((0, 4)), dtype=float).reshape(-1, 4)
    tree = cKDTree(centres)
    longest_bond = float(np.max(np.linalg.norm(centres[pairs[:, 0]] - centres[pairs[:, 1]], axis=1))) \
        if len(pairs) else 0.0
    reach = 2.0 * float(radii.max()) + gap + 3.0 * float(np.hypot(*half.max(axis=0))) + 4.0
    preference = np.arange(2 * len(LABEL_DIRECTIONS)) * 0.01

    for atom in range(count):
        half_w, half_h = half[atom]
        support = np.abs(LABEL_DIRECTIONS[:, 0]) * half_w + np.abs(LABEL_DIRECTIONS[:, 1]) * half_h
        near = radii[atom] + gap + support
        far = near + 2.0 * half_h + 2.0
        candidates = np.vstack([centres[atom] + near[:, None] * LABEL_DIRECTIONS,
                                centres[atom] + far[:, None] * LABEL_DIRECTIONS])
        cost = preference.copy()
        cost[len(LABEL_DIRECTIONS):] += 0.5
        local = [j for j in tree.query_ball_point(centres[atom], reach + longest_bond) if j != atom]
        # Mildly prefer pointing away from the surrounding atoms, so a ring
        # atom is numbered outside the ring and a chain end beyond the chain
        # (well below any collision cost, so it only chooses among free spots).
        if local:
            away = centres[atom] - centres[local].mean(axis=0)
            if np.linalg.norm(away) > 1e-6:
                away /= np.linalg.norm(away)
                alignment = np.tile(LABEL_DIRECTIONS @ away, 2)
                cost += 0.6 * (1.0 - alignment) / 2.0
        if local:
            others = centres[local]
            dx = np.maximum(np.abs(candidates[:, None, 0] - others[None, :, 0]) - half_w, 0.0)
            dy = np.maximum(np.abs(candidates[:, None, 1] - others[None, :, 1]) - half_h, 0.0)
            hits = dx * dx + dy * dy < (radii[local][None, :] + 1.0) ** 2
            cost += 10.0 * hits.sum(axis=1)

        if samples is not None:
            nearby_bonds = sorted({b for j in local + [atom] for b in bonds_of[j]})
            if nearby_bonds:
                points = samples[nearby_bonds].reshape(-1, 2)
                inside = ((np.abs(candidates[:, None, 0] - points[None, :, 0]) < half_w + 1.0)
                          & (np.abs(candidates[:, None, 1] - points[None, :, 1]) < half_h + 1.0))
                cost += 10.0 * inside.reshape(len(candidates), len(nearby_bonds), -1).any(axis=2).sum(axis=1)

        if len(placed):
            overlap = ((candidates[:, None, 0] - half_w < placed[None, :, 2])
                       & (candidates[:, None, 0] + half_w > placed[None, :, 0])
                       & (candidates[:, None, 1] - half_h < placed[None, :, 3])
                       & (candidates[:, None, 1] + half_h > placed[None, :, 1]))
            cost += 6.0 * overlap.sum(axis=1)

        best = candidates[int(np.argmin(cost))]
        result[atom] = best
        placed = np.vstack([placed, [best[0] - half_w, best[1] - half_h, best[0] + half_w, best[1] + half_h]])
    return result


def render_scene(primitives: Sequence[tuple], size: Tuple[int, int], scale: float,
                 supersample: int = SUPERSAMPLE):
    """A Pillow RGB image of ``size`` physical pixels; logical units times ``scale``."""
    from PIL import Image, ImageDraw

    k = scale * supersample
    width, height = size
    image = Image.new("RGB", (width * supersample, height * supersample), "#ffffff")
    draw = ImageDraw.Draw(image)

    def stroke_width(width_logical):
        return max(1, int(round(width_logical * k)))

    def line(points, color, width_logical, dash):
        width_px = stroke_width(width_logical)
        for (xa, ya), (xb, yb) in zip(points, points[1:]):
            if dash:
                pieces = dash_segments(xa * k, ya * k, xb * k, yb * k, [d * k for d in dash])
            else:
                pieces = [((xa * k, ya * k), (xb * k, yb * k))]
            for start, end in pieces:
                draw.line([start, end], fill=color, width=width_px)

    for primitive in primitives:
        kind = primitive[0]
        if kind == "rect":
            _, x0, y0, x1, y1, fill = primitive
            draw.rectangle([x0 * k, y0 * k, x1 * k, y1 * k], fill=fill)
        elif kind == "disc":
            _, cx, cy, radius, fill, outline, outline_width = primitive
            box = [(cx - radius) * k, (cy - radius) * k, (cx + radius) * k, (cy + radius) * k]
            draw.ellipse(box, fill=fill, outline=outline,
                         width=stroke_width(outline_width) if outline else 0)
        elif kind == "ring":
            _, cx, cy, radius, color, width_logical = primitive
            box = [(cx - radius) * k, (cy - radius) * k, (cx + radius) * k, (cy + radius) * k]
            draw.ellipse(box, outline=color, width=stroke_width(width_logical))
        elif kind == "line":
            _, x1, y1, x2, y2, color, width_logical, dash = primitive
            line([(x1, y1), (x2, y2)], color, width_logical, dash)
        elif kind == "polyline":
            _, points, color, width_logical, dash = primitive
            line(list(points), color, width_logical, dash)
        elif kind == "text":
            _, x, y, text, color, size_points, bold, baseline = primitive
            font = _font(max(1, int(round(size_points * POINTS_TO_PIXELS * k))), bool(bold))
            draw.text((x * k, y * k), text, fill=color, font=font, anchor=ANCHORS.get(baseline, "ls"))
        elif kind == "label":
            _, x, y, text, color, size_points, bold = primitive
            font = _font(max(1, int(round(size_points * POINTS_TO_PIXELS * k))), bool(bold))
            draw.text((x * k, y * k), text, fill=color, font=font, anchor="mm")
        else:
            raise ValueError(f"Unknown canvas primitive: {kind}")

    if supersample > 1:
        # reduce() is an integer box filter: the same average as
        # resize(BOX) at this factor, about 4x faster (1.2 vs 4.8 ms measured).
        image = image.reduce(supersample)
        if image.size != (width, height):
            image = image.resize((width, height), Image.BOX)
    return image


def pil_to_net_bitmap(image):
    """A 32-bit ``System.Drawing.Bitmap`` holding ``image``'s pixels (WinForms only).

    ``toga.Image(pil)`` would encode a PNG and convert it byte by byte through
    pythonnet (~170 ms per frame, measured for boxConstructionView); locking
    the bitmap and ``memmove``-ing the rows in takes about a millisecond.
    32 bits per pixel is GDI+'s native layout, so the blit needs no format
    conversion (0.9 ms against 1.4 ms for 24 bits, measured). Raises on any
    platform without pythonnet/System.Drawing.
    """
    import ctypes

    from System.Drawing import Bitmap, Rectangle
    from System.Drawing.Imaging import ImageLockMode, PixelFormat

    rgb = image.convert("RGB")
    width, height = rgb.size
    bitmap = Bitmap(width, height, PixelFormat.Format32bppRgb)
    locked = bitmap.LockBits(Rectangle(0, 0, width, height), ImageLockMode.WriteOnly,
                             PixelFormat.Format32bppRgb)
    try:
        stride = int(locked.Stride)
        raw = rgb.tobytes("raw", "BGRX", stride, 1)  # GDI+ byte order; X is ignored
        ctypes.memmove(int(locked.Scan0.ToInt64()), raw, stride * height)
    finally:
        bitmap.UnlockBits(locked)
    return bitmap
