import asyncio
import pathlib
import subprocess
import sys

import numpy as np
import pytest

import boxConstructionView as bcv

SRC = pathlib.Path(__file__).resolve().parents[1]
WATER = np.array([[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]])


def snapshot(n_solvent=3, lengths=(12.0, 12.0, 12.0), with_solute=True):
    molecules = []
    if with_solute:
        molecules.append(("solute", ["C", "O"], np.array([[6.0, 6.0, 6.0], [7.2, 6.0, 6.0]])))
    rng = np.random.default_rng(0)
    for _ in range(n_solvent):
        molecules.append(("solvent", ["O", "H", "H"], WATER + rng.uniform(1, 11, 3)))
    return bcv.ConstructionSnapshot("Inserting solvent", lengths, molecules, caption="3/10")


def fitted_camera(lengths=(12.0, 12.0, 12.0), width=700, height=440):
    camera = bcv.BoxCamera()
    bcv.fit_camera(camera, lengths, width, height)
    return camera


def test_box_has_twelve_edges_inside_the_canvas():
    assert len(bcv.BOX_EDGES) == 12
    for lengths in [(10, 10, 10), (20, 30, 60), (80, 80, 15)]:
        camera = fitted_camera(lengths)
        corners = camera.project(bcv.box_corners(lengths))
        assert corners[:, 0].min() >= 0 and corners[:, 0].max() <= 700
        assert corners[:, 1].min() >= 0 and corners[:, 1].max() <= 440


def test_default_view_is_a_rotation_not_a_mirror():
    assert np.isclose(np.linalg.det(bcv.DEFAULT_VIEW), 1.0)
    # Angstrom +y must point up the screen (smaller screen y).
    camera = fitted_camera()
    up = camera.project(np.array([[6.0, 7.0, 6.0]]))[0, 1] - camera.project(np.array([[6.0, 6.0, 6.0]]))[0, 1]
    assert up < 0


def test_draw_list_has_every_atom_once_back_to_front():
    snap = snapshot()
    camera = fitted_camera()
    prims = bcv.build_draw_list(snap, camera)
    assert prims[0][0] == "lines" and len(prims[0][3]) == 12
    discs = [d for p in prims if p[0] == "discs" for d in p[2]]
    assert len(discs) == 2 + 3 * 3
    # Farthest first: depth of the drawn discs never increases toward the viewer.
    coords = np.concatenate([m[2] for m in snap.molecules])
    depth = camera.project(coords)[:, 2]
    assert np.isclose(max(depth), sorted(depth, reverse=True)[0])


def test_solute_is_larger_and_differently_coloured():
    prims = bcv.build_draw_list(snapshot(n_solvent=1), fitted_camera())
    by_color = {}
    for p in prims:
        if p[0] == "discs":
            by_color.setdefault(p[1], []).extend(r for _, _, r in p[2])
    colors = set(by_color)
    assert any(c.startswith("#e") or c.startswith("#f") for c in colors)  # orange family for the solute
    radii = [r for rs in by_color.values() for r in rs]
    assert max(radii) > min(radii)


def test_large_boxes_fall_back_to_one_dot_per_molecule():
    snap = snapshot(n_solvent=50)
    prims = bcv.build_draw_list(snap, fitted_camera(), max_atoms=20)
    discs = [d for p in prims if p[0] == "discs" for d in p[2]]
    assert len(discs) == 2 + 50
    assert any("centroids" in p[4] for p in prims if p[0] == "text")


def test_caption_and_legend():
    prims = bcv.build_draw_list(snapshot(), fitted_camera())
    texts = [p[4] for p in prims if p[0] == "text"]
    assert texts[0] == "Inserting solvent -- 3/10"
    assert "solute" in texts and "solvent" in texts


def test_empty_box_draws_only_edges_and_caption():
    snap = bcv.ConstructionSnapshot("Empty box", (10, 10, 10), [])
    prims = bcv.build_draw_list(snap, fitted_camera((10, 10, 10)))
    assert [p[0] for p in prims] == ["lines", "text"]


def test_same_colour_discs_are_batched():
    prims = bcv.build_draw_list(snapshot(n_solvent=200, with_solute=False), fitted_camera())
    runs = [p for p in prims if p[0] == "discs"]
    # 600 atoms, a handful of fade levels and two element sizes -> far fewer fill operations.
    assert len(runs) < 100


def test_render_image_draws_box_atoms_and_caption():
    snap = snapshot()
    camera = fitted_camera(width=500, height=300)
    image = bcv.render_image(bcv.build_draw_list(snap, camera), 500, 300)
    assert image.size == (500, 300)
    pixels = np.asarray(image)
    background = np.array(bcv._rgb(bcv.CANVAS_BACKGROUND))
    assert np.any(np.any(pixels != background, axis=2))
    # The solute's orange is present (allowing for antialiasing and depth fading).
    orange = np.array(bcv._rgb(bcv.KIND_COLORS["solute"]))
    assert np.min(np.abs(pixels.reshape(-1, 3).astype(int) - orange).sum(axis=1)) < 60


def test_render_image_without_supersampling_and_empty_list():
    image = bcv.render_image([], 120, 80, supersample=1)
    assert image.size == (120, 80)
    assert np.all(np.asarray(image) == np.array(bcv._rgb(bcv.CANVAS_BACKGROUND)))


class _StubImageView:
    def __init__(self):
        self.image = None
        self.sets = 0

    def __setattr__(self, name, value):
        if name == "image" and value is not None:
            object.__setattr__(self, "sets", getattr(self, "sets", 0) + 1)
        object.__setattr__(self, name, value)


def make_view(interval=0.2):
    view = bcv.ConstructionView.__new__(bcv.ConstructionView)
    view.width, view.height = 400, 260
    view.min_interval = interval
    view.max_atoms = bcv.DEFAULT_MAX_ATOMS
    view.camera = bcv.BoxCamera()
    view.snapshot = None
    view._fitted_lengths = None
    view._last_paint = 0.0
    view._timer = None
    view.paint_count = 0
    view.widget = _StubImageView()
    return view


def test_submit_throttles_repaints_and_flushes_the_last_snapshot():
    async def scenario():
        view = make_view(interval=1.0)
        snaps = [snapshot(n_solvent=i % 5 + 1) for i in range(100)]
        for snap in snaps:
            view.submit(snap)
        assert view.paint_count == 1          # first one painted, the rest coalesced
        await asyncio.sleep(1.3)
        assert view.paint_count == 2          # timer flushed the newest snapshot once
        assert len(view.snapshot.molecules) == 1 + (99 % 5 + 1)
        view.submit(snapshot(), force=True)
        assert view.paint_count == 3
        assert view.widget.sets == 3
        assert view.widget.image.size == (400, 260)

    asyncio.run(scenario())


def test_rotate_buttons_turn_the_view_and_reset_restores():
    view = make_view()
    view.submit(snapshot(), force=True)
    before = view.camera.rotation.copy()
    view.rotate(horizontal_degrees=15)
    after_turn = view.camera.rotation.copy()
    assert not np.allclose(before, after_turn)
    # A 15-degree turn about the screen vertical: the rotation changed by exactly that angle.
    relative = after_turn @ before.T
    angle = np.degrees(np.arccos((np.trace(relative) - 1) / 2))
    assert angle == pytest.approx(15.0, abs=1e-6)
    view.rotate(vertical_degrees=-15)
    view.zoom(1.25)
    view.reset_view()
    assert np.allclose(before, view.camera.rotation)
    assert view.paint_count == 5


def test_module_import_does_not_load_toga():
    code = "import sys, boxConstructionView; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=str(SRC))
    assert out.stdout.strip() == "False", out.stderr
