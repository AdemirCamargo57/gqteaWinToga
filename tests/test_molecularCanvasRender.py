"""Tests for molecularCanvasRender: the Pillow frame behind the design canvas.

The native-bitmap blit needs WinForms and is checked live; everything that
decides what the frame looks like is tested here, headless.
"""
import subprocess
import sys

import numpy as np
import pytest

from molecularCanvasRender import dash_segments, place_atom_labels, render_scene, text_box_size


def test_the_module_does_not_import_toga():
    result = subprocess.run(
        [sys.executable, "-c", "import molecularCanvasRender, sys; print('toga' in sys.modules)"],
        capture_output=True, text=True,
    )
    assert result.stdout.strip() == "False", result.stderr


def test_dashes_alternate_on_and_off_along_the_line():
    pieces = dash_segments(0.0, 0.0, 10.0, 0.0, [3.0, 2.0])
    assert pieces == [((0.0, 0.0), (3.0, 0.0)), ((5.0, 0.0), (8.0, 0.0))]


def test_a_degenerate_dash_is_a_plain_line():
    assert dash_segments(1.0, 1.0, 1.0, 1.0, [3.0, 2.0]) == [((1.0, 1.0), (1.0, 1.0))]
    assert dash_segments(0.0, 0.0, 4.0, 0.0, []) == [((0.0, 0.0), (4.0, 0.0))]


@pytest.mark.parametrize("scale, size", [(1.0, (100, 80)), (1.25, (125, 100))])
def test_the_frame_has_the_physical_size_and_draws_in_logical_units(scale, size):
    primitives = [
        ("rect", 0, 0, 100, 80, "#102030"),
        ("disc", 50.0, 40.0, 10.0, "#ff0000", "#000000", 1.0),
    ]
    image = render_scene(primitives, size, scale)
    assert image.size == size
    assert image.getpixel((2, 2)) == (0x10, 0x20, 0x30)
    centre = (round(50 * scale), round(40 * scale))
    assert image.getpixel(centre) == (255, 0, 0)


def test_every_primitive_kind_renders():
    primitives = [
        ("rect", 0, 0, 120, 60, "#ffffff"),
        ("disc", 20.0, 20.0, 8.0, "#00ff00", None, 0.0),
        ("ring", 50.0, 20.0, 8.0, "#0000ff", 2.0),
        ("line", 0.0, 50.0, 120.0, 50.0, "#000000", 1.8, None),
        ("line", 0.0, 55.0, 120.0, 55.0, "#000000", 1.8, (5.0, 4.0)),
        ("polyline", [(70.0, 10.0), (90.0, 30.0), (110.0, 10.0)], "#d6336c", 1.2, (2.0, 3.0)),
        ("text", 5.0, 5.0, "C", "#101010", 9, False, "middle"),
        ("text", 60.0, 40.0, "#2", "#d6336c", 9, True, "bottom"),
    ]
    image = render_scene(primitives, (120, 60), 1.0)
    colours = {colour for _, colour in image.getcolors(maxcolors=100000)}
    assert (0, 255, 0) in colours and len(colours) > 10  # shapes, anti-aliased edges


def test_an_unknown_primitive_is_an_error():
    with pytest.raises(ValueError, match="Unknown canvas primitive"):
        render_scene([("hexagon", 0, 0)], (10, 10), 1.0)


def test_a_centred_label_renders_around_its_point():
    image = render_scene([("rect", 0, 0, 60, 30, "#ffffff"),
                          ("label", 30.0, 15.0, "88", "#000000", 12, True)], (60, 30), 1.0)
    dark = [(x, y) for x in range(60) for y in range(30) if sum(image.getpixel((x, y))) < 200]
    xs = [x for x, _ in dark]
    assert dark and abs((min(xs) + max(xs)) / 2.0 - 30.0) <= 2.0


# ------------------------------------------------------------------
# Atom-number placement
# ------------------------------------------------------------------
RADIUS = 11.0


def _benzene_screen():
    angles = np.radians(np.arange(6) * 60.0)
    ring = np.c_[200 + 56 * np.cos(angles), 200 + 56 * np.sin(angles)]
    hydrogens = np.c_[200 + 99 * np.cos(angles), 200 + 99 * np.sin(angles)]
    bonds = [(k, (k + 1) % 6) for k in range(6)] + [(k, k + 6) for k in range(6)]
    return np.vstack([ring, hydrogens]), bonds


def _boxes(positions, sizes):
    half = np.asarray(sizes) / 2.0
    return np.c_[positions - half, positions + half]


def _box_hits_disc(box, centre, radius):
    dx = max(box[0] - centre[0], 0.0, centre[0] - box[2])
    dy = max(box[1] - centre[1], 0.0, centre[1] - box[3])
    return dx * dx + dy * dy < radius * radius


def _box_hits_segment(box, start, end, samples=200):
    t = np.linspace(0.0, 1.0, samples)[:, None]
    points = np.asarray(start) * (1 - t) + np.asarray(end) * t
    inside = ((points[:, 0] > box[0]) & (points[:, 0] < box[2])
              & (points[:, 1] > box[1]) & (points[:, 1] < box[3]))
    return bool(inside.any())


def test_numbers_cover_no_atom_no_bond_and_no_other_number():
    centres, bonds = _benzene_screen()
    sizes = [text_box_size(str(k + 1), 8) for k in range(len(centres))]

    positions = place_atom_labels(centres, RADIUS, bonds, sizes)

    boxes = _boxes(positions, sizes)
    for box in boxes:
        assert not any(_box_hits_disc(box, centre, RADIUS) for centre in centres)
        assert not any(_box_hits_segment(box, centres[i], centres[j]) for i, j in bonds)
    for first in range(len(boxes)):
        for second in range(first + 1, len(boxes)):
            a, b = boxes[first], boxes[second]
            assert a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1]


def test_numbers_sit_next_to_their_atom_and_outside_the_ring():
    centres, bonds = _benzene_screen()
    sizes = [text_box_size(str(k + 1), 8) for k in range(len(centres))]

    positions = place_atom_labels(centres, RADIUS, bonds, sizes)

    distances = np.linalg.norm(positions - centres, axis=1)
    assert (distances < RADIUS + 25.0).all()  # beside the atom, not drifting away
    ring_centre = np.array([200.0, 200.0])
    assert (np.linalg.norm(positions - ring_centre, axis=1)
            > np.linalg.norm(centres - ring_centre, axis=1)).all()


def test_numbers_keep_clear_of_other_text():
    centres = np.array([[100.0, 100.0]])
    size = text_box_size("1", 8)
    free = place_atom_labels(centres, RADIUS, [], [size])
    blocking = [(free[0][0] - 20, free[0][1] - 20, free[0][0] + 20, free[0][1] + 20)]

    moved = place_atom_labels(centres, RADIUS, [], [size], obstacles=blocking)

    box = _boxes(moved, [size])[0]
    assert box[2] <= blocking[0][0] or box[0] >= blocking[0][2] \
        or box[3] <= blocking[0][1] or box[1] >= blocking[0][3]


def test_no_atoms_no_labels():
    assert place_atom_labels(np.zeros((0, 2)), RADIUS, [], np.zeros((0, 2))).shape == (0, 2)
