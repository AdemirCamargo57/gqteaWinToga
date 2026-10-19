"""Tests for the 3D viewer's Cube tab controller, headless (stub viewer and widgets)."""
import asyncio
from types import SimpleNamespace

import numpy as np
import pytest

import cubeFile as cf
import cubeTab as ct
import volumeData as vd
from cubeRender import CubeSceneLayer


def p_orbital_volume(shape=(21, 21, 21), h=0.4, orbital=None, atoms=((8, (0.0, 0.0, 0.0)), (1, (1.8, 0.0, 0.0)))):
    origin = -h * (np.array(shape) - 1) / 2
    idx = np.stack(np.meshgrid(*[np.arange(n) for n in shape], indexing="ij"), axis=-1).reshape(-1, 3)
    points = origin + idx * h
    values = points[:, 0] * np.exp(-(points ** 2).sum(axis=1) / 2)
    z = np.array([a[0] for a in atoms], dtype=int)
    xyz = np.array([a[1] for a in atoms], dtype=float).reshape(-1, 3)
    return cf.CubeVolume("p", "test", z, z.astype(float), xyz, origin, np.eye(3) * h,
                         values.reshape(shape), orbital)


def write(tmp_path, name, volume):
    return cf.write_cube_volume(tmp_path / name, volume)


class StubViewer:
    def __init__(self, frames=None):
        self.frames = list(frames or [])
        self.loaded = []
        self.cube_layer = CubeSceneLayer()
        self.show_structure = True
        self.invalidated = 0
        self.projection_requests = 0
        self.status = []

    def load_designed_molecule(self, frame, bonds=None):
        self.loaded.append(frame)
        self.frames = [list(frame)]

    def get_current_frame_data(self):
        return list(self.frames[0]) if self.frames else []

    def invalidate_scene_cache(self):
        self.invalidated += 1

    def request_projection_update(self):
        self.projection_requests += 1

    def set_status_message(self, message):
        self.status.append(message)


def controller(frames=None, use_numba=False):
    c = ct.CubeTabController(StubViewer(frames))
    c.use_numba = use_numba
    return c


def loaded(tmp_path, *volumes, frames=None):
    c = controller(frames)
    for k, volume in enumerate(volumes):
        c.add_entries(c.read_entries(write(tmp_path, f"cube{k}.cube", volume)))
    return c


class TestLoading:
    def test_plain_cube_is_named_after_the_file(self, tmp_path):
        entries = ct.CubeTabController.read_entries(write(tmp_path, "water.cube", p_orbital_volume()))
        assert [e.name for e in entries] == ["water"]
        entry = entries[0]
        assert entry.saved and entry.path.endswith("water.cube")
        assert entry.isovalue == vd.default_isovalue(entry.volume)
        assert entry.statistics == pytest.approx(vd.value_statistics(entry.volume))

    def test_orbital_cube_gives_one_entry_per_orbital(self, tmp_path):
        volume = p_orbital_volume()
        text = cf.format_cube_volume(volume).splitlines()
        # Turn it into an orbital cube with two MOs (values duplicated, the second negated).
        n_atoms = len(volume.atomic_numbers)
        text[2] = f"{-n_atoms:5d}" + text[2][5:]
        data = np.stack([volume.values, -volume.values], axis=-1).reshape(21 * 21, -1)
        body = ["    2    7    8"]
        for column in data:
            for k in range(0, len(column), 6):
                body.append("".join(f"{v:13.5E}" for v in column[k:k + 6]))
        path = tmp_path / "mo.cube"
        path.write_text("\n".join(text[:6 + n_atoms] + body) + "\n")
        entries = ct.CubeTabController.read_entries(path)
        assert [e.name for e in entries] == ["mo MO 7", "mo MO 8"]
        assert entries[0].isovalue == 0.02

    def test_density_cubes_show_only_the_positive_surface(self, tmp_path):
        density = p_orbital_volume()
        density.values = density.values ** 2
        (entry,) = ct.CubeTabController.read_entries(write(tmp_path, "rho.cube", density))
        assert entry.show == ct.SHOW_POSITIVE
        (signed,) = ct.CubeTabController.read_entries(write(tmp_path, "p.cube", p_orbital_volume()))
        assert signed.show == ct.SHOW_BOTH

    def test_names_are_unique(self, tmp_path):
        c = controller()
        path = write(tmp_path, "a.cube", p_orbital_volume())
        c.add_entries(c.read_entries(path))
        c.add_entries(c.read_entries(path))
        c.add_entries(c.read_entries(path))
        assert c.names() == ["a", "a (2)", "a (3)"]

    def test_first_cube_sets_the_slice_cube(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        assert c.slice_cube == "cube0" and c.slice_index == 10 and c.slice_axis == 2


class TestStructure:
    def test_question_only_when_a_structure_is_displayed(self):
        assert not controller().structure_question_needed()
        assert controller(frames=[[("C", (0.0, 0.0, 0.0))]]).structure_question_needed()

    def test_no_question_when_the_cube_holds_the_displayed_molecule(self, tmp_path):
        """Regression: the molecule on screen (e.g. from the ESP tab) is the cube's own one."""
        (entry,) = ct.CubeTabController.read_entries(write(tmp_path, "same.cube", p_orbital_volume()))
        c = controller(frames=[vd.structure_frame(entry.volume)])
        assert not c.structure_question_needed(entry)
        asked = []

        async def ask(title, message):
            asked.append(title)
            return True
        c._ask = ask
        assert asyncio.run(c._settle_structure(entry)) is True
        assert asked == [] and c.viewer.loaded == []  # kept as is, not reloaded
        assert c.structure_volume is entry.volume  # later cubes are compared with it

    def test_a_different_displayed_molecule_is_still_asked_about(self, tmp_path):
        (entry,) = ct.CubeTabController.read_entries(write(tmp_path, "other.cube", p_orbital_volume()))
        c = controller(frames=[[("C", (0.0, 0.0, 0.0))]])
        assert c.structure_question_needed(entry)
        asked = []

        async def ask(title, message):
            asked.append(message)
            return False
        c._ask = ask
        assert asyncio.run(c._settle_structure(entry)) is True
        assert len(asked) == 1 and "(1 atoms)" in asked[0] and c.viewer.loaded == []

    def test_esp_tab_molecule_then_its_cube_asks_nothing(self, tmp_path):
        """The reported sequence, with a real viewer: ESP run shown, then its ESP cube loaded."""
        import electrostaticSurfacePotential as esp
        from espSurfaceLayer import ESPSurfaceLayer
        from esp_fixtures import fchk_path
        from molecularViewer import MolecularViewer
        viewer = MolecularViewer()
        viewer.esp_layer = ESPSurfaceLayer()
        viewer.cube_layer = CubeSceneLayer()
        analyser = esp.ESPSurfaceAnalyser()
        results = analyser.run(str(fchk_path("water_631gd_6d")), [0.001], spacing=0.4)
        tab = esp.ESPTabController.__new__(esp.ESPTabController)
        tab.viewer, tab.results, tab._structure_settled_for = viewer, results, None
        asyncio.run(tab._settle_structure())
        assert len(viewer.frames) == 1
        paths, _ = analyser.write_cubes(results, tmp_path / "water_ESP.cube", with_density=False)
        c = ct.CubeTabController(viewer)
        (entry,) = c.read_entries(paths[0])
        assert not c.structure_question_needed(entry)

    def test_adopting_loads_the_cube_atoms_in_angstrom(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        c.adopt_structure(c.entries[0])
        (frame,) = c.viewer.loaded
        assert [s for s, _ in frame] == ["O", "H"]
        assert frame[1][1] == pytest.approx((1.8 * cf.BOHR_RADIUS_ANGSTROM, 0.0, 0.0))
        assert c.structure_volume is c.entries[0].volume
        assert not c.structure_question_needed()

    def test_a_cube_without_atoms_cannot_supply_a_structure(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume(atoms=()))
        with pytest.raises(ValueError, match="no atoms"):
            c.adopt_structure(c.entries[0])

    def test_later_cubes_with_other_atoms_give_a_note(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        c.adopt_structure(c.entries[0])
        notes = c.add_entries(c.read_entries(write(tmp_path, "b.cube", p_orbital_volume(atoms=((6, (0, 0, 0)),)))))
        assert len(notes) == 1 and "different atoms" in notes[0] and "unchanged" in notes[0]
        assert len(c.viewer.loaded) == 1


class TestScene:
    def test_both_lobes_in_their_colours(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        c.entries[0].isovalue = 0.1
        items, notes = c.compute_scene_items()
        assert [i.kind for i in items] == ["surface", "surface"] and notes == []
        assert items[0].rgba[0, :3] == pytest.approx(vd.COLOUR_PRESETS["Blue"])
        assert items[1].rgba[0, :3] == pytest.approx(vd.COLOUR_PRESETS["Red"])

    def test_show_modes_visibility_and_master_switch(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        entry = c.entries[0]
        entry.isovalue = 0.1
        entry.show = ct.SHOW_POSITIVE
        assert len(c.compute_scene_items()[0]) == 1
        entry.show = ct.SHOW_BOTH
        entry.visible = False
        assert c.compute_scene_items()[0] == []
        entry.visible = True
        c.show_volumes = False
        assert c.compute_scene_items()[0] == []

    def test_no_surface_is_a_note_not_an_error(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        c.entries[0].isovalue = 50.0
        items, notes = c.compute_scene_items()
        assert items == [] and len(notes) == 1 and "no surface" in notes[0]

    def test_opacity_reaches_the_items(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        c.entries[0].isovalue = 0.1
        c.entries[0].opacity = 0.5
        assert all(not i.opaque for i in c.compute_scene_items()[0])

    def test_meshes_are_cached_per_isovalue(self, tmp_path, monkeypatch):
        c = loaded(tmp_path, p_orbital_volume())
        c.entries[0].isovalue = 0.1
        calls = []
        real = vd.isosurface
        monkeypatch.setattr(ct.vd, "isosurface", lambda *a, **k: calls.append(a) or real(*a, **k))
        c.compute_scene_items()
        c.compute_scene_items()
        assert len(calls) == 2
        c.entries[0].isovalue = 0.12
        c.compute_scene_items()
        assert len(calls) == 4

    def test_slice_item(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        c.entries[0].visible = False
        c.set_slice(visible=True)
        (item,) = c.compute_scene_items()[0]
        assert item.kind == "slice" and item.vertices.shape == (21 * 21, 3)

    def test_refresh_pushes_items_to_the_viewer(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        c.entries[0].isovalue = 0.1
        asyncio.run(c.refresh_scene())
        assert len(c.viewer.cube_layer.items()) == 2
        assert c.viewer.invalidated and c.viewer.projection_requests

    def test_hiding_the_molecule(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        c.set_show_structure(False)
        assert c.viewer.show_structure is False


class TestSlices:
    def test_axis_change_moves_to_the_middle_plane(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume(shape=(21, 15, 11)))
        c.set_slice(axis=1)
        assert c.slice_index == 7 and c.slice_plane_count() == 15

    def test_index_is_validated(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        with pytest.raises(ValueError, match="plane"):
            c.set_slice(index=21)
        with pytest.raises(ValueError, match="range"):
            c.set_slice(range_mode="other")

    def test_slice_text(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume())
        text = c.slice_text()
        assert len([r for r in text.splitlines() if not r.startswith("#")]) == 21 * 21
        assert "cube0" in text

    def test_removing_the_slice_cube_moves_the_slice(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume(), p_orbital_volume())
        c.remove("cube0")
        assert c.names() == ["cube1"] and c.slice_cube == "cube1"
        c.remove("cube1")
        assert c.slice_cube is None and c.compute_scene_items()[0] == []
        with pytest.raises(ValueError, match="No cube"):
            c.remove("cube1")


class TestOperations:
    def test_difference_is_added_unsaved(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume(), p_orbital_volume())
        volume, name, notes = c.compute_operation("cube0", ct.OP_SUBTRACT, "cube1", "", "")
        entry = c.add_result(volume, name)
        assert entry.name == "cube0 − cube1" and not entry.saved and entry.path is None
        assert not entry.volume.values.any() and notes == []
        assert c.names()[-1] == entry.name

    def test_sum_and_scale(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume(), p_orbital_volume())
        a = c.entries[0].volume.values
        volume, name, _ = c.compute_operation("cube0", ct.OP_ADD, "cube1", "", "rho sum")
        assert name == "rho sum" and np.allclose(volume.values, 2 * a)
        volume, name, _ = c.compute_operation("cube0", ct.OP_SCALE, None, "-0.5", "")
        assert name == "-0.5 × cube0" and np.allclose(volume.values, -0.5 * a)

    def test_refusals(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume(), p_orbital_volume(shape=(21, 21, 19)))
        with pytest.raises(ValueError, match="same grid"):
            c.compute_operation("cube0", ct.OP_ADD, "cube1", "", "")
        with pytest.raises(ValueError, match="factor"):
            c.compute_operation("cube0", ct.OP_SCALE, None, "two", "")
        with pytest.raises(ValueError, match="cube B"):
            c.compute_operation("cube0", ct.OP_ADD, None, "", "")
        with pytest.raises(ValueError, match="operation"):
            c.compute_operation("cube0", "A / B", "cube1", "", "")

    def test_different_atoms_give_a_note(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume(), p_orbital_volume(atoms=((6, (0, 0, 0)),)))
        _, _, notes = c.compute_operation("cube0", ct.OP_ADD, "cube1", "", "")
        assert len(notes) == 1 and "atoms" in notes[0]


class TestTableColumns:
    def test_columns_are_refitted_whenever_the_native_table_is_resized(self):
        """Toga's WinForms Table splits its width between the columns once, at its first
        layout; on a hidden tab that is a few pixels, so every column stayed ~19 px."""
        class Event:
            def __init__(self):
                self.handlers = []

            def __iadd__(self, handler):
                self.handlers.append(handler)
                return self

        refits = []
        impl = SimpleNamespace(native=SimpleNamespace(Resize=Event()), _resize_columns=lambda: refits.append(1))
        table = SimpleNamespace(_impl=impl)
        ct.fit_table_columns_on_resize(table)
        for handler in impl.native.Resize.handlers:
            handler(None, None)
        assert refits == [1]

    def test_other_backends_are_left_alone(self):
        ct.fit_table_columns_on_resize(SimpleNamespace())  # no _impl.native: nothing to do, no error


class TestSaving:
    def test_save_marks_the_entry_saved(self, tmp_path):
        c = loaded(tmp_path, p_orbital_volume(), p_orbital_volume())
        volume, name, _ = c.compute_operation("cube0", ct.OP_SUBTRACT, "cube1", "", "")
        entry = c.add_result(volume, name)
        path = c.save_entry(entry.name, tmp_path / "out" / ct.default_cube_filename(entry.name))
        assert entry.saved and entry.path == path
        (back,) = cf.read_cube_volumes(path)
        assert back.values.shape == (21, 21, 21)

    def test_default_file_names(self):
        assert ct.default_cube_filename("cube0 − cube1") == "cube0_cube1.cube"
        assert ct.default_cube_filename("mo MO 7") == "mo_MO_7.cube"
        assert ct.default_cube_filename("***") == "cube.cube"


# ----------------------------------------------------------------------
# Widget plumbing (no Toga window: stub widgets on a real controller)
# ----------------------------------------------------------------------
def with_widgets(c):
    c.iso_cube_selection = SimpleNamespace(value=c.names()[0] if c.names() else None, items=c.names())
    c.isovalue_input = SimpleNamespace(value="0.15")
    c.show_selection = SimpleNamespace(value=ct.SHOW_NEGATIVE)
    c.positive_colour_input = SimpleNamespace(value="green")
    c.negative_colour_input = SimpleNamespace(value="#00FF00")
    c.opacity_slider = SimpleNamespace(value=0.4)
    c.visible_switch = SimpleNamespace(value=False)
    c.slice_cube_selection = SimpleNamespace(value=c.names()[0] if c.names() else None)
    c.slice_axis_selection = SimpleNamespace(value=ct.AXIS_LABELS[0])
    c.slice_slider = SimpleNamespace(value=3.0, min=0, max=20)
    c.slice_range_selection = SimpleNamespace(value=list(ct.RANGE_LABELS)[1])
    c.slice_switch = SimpleNamespace(value=True)
    c.status_label = SimpleNamespace(text="")
    return c


class TestWidgetPlumbing:
    def test_every_isosurface_setting_comes_from_the_widgets(self, tmp_path):
        c = with_widgets(loaded(tmp_path, p_orbital_volume()))
        assert c.apply_isosurface_fields() is None
        e = c.entries[0]
        assert (e.isovalue, e.show, e.opacity, e.visible) == (0.15, ct.SHOW_NEGATIVE, 0.4, False)
        assert (e.positive_colour, e.negative_colour) == ("Green", "#00FF00")

    @pytest.mark.parametrize("field, text, message", [
        ("isovalue_input", "-1", "positive"),
        ("isovalue_input", "abc", "positive"),
        ("positive_colour_input", "teal-ish", "colour"),
    ])
    def test_bad_isosurface_fields_change_nothing(self, tmp_path, field, text, message):
        c = with_widgets(loaded(tmp_path, p_orbital_volume()))
        before = (c.entries[0].isovalue, c.entries[0].positive_colour, c.entries[0].show)
        getattr(c, field).value = text
        error = c.apply_isosurface_fields()
        assert message in error
        assert (c.entries[0].isovalue, c.entries[0].positive_colour, c.entries[0].show) == before

    def test_every_slice_setting_comes_from_the_widgets(self, tmp_path):
        c = with_widgets(loaded(tmp_path, p_orbital_volume()))
        c.slice_axis_selection.value = ct.AXIS_LABELS[2]
        assert c.apply_slice_fields() is None
        assert (c.slice_cube, c.slice_axis, c.slice_index) == ("cube0", 2, 3)
        assert c.slice_range == "minmax" and c.slice_visible is True

    def test_a_new_axis_starts_at_its_middle_plane(self, tmp_path):
        c = with_widgets(loaded(tmp_path, p_orbital_volume(shape=(21, 15, 11))))
        c.slice_axis_selection.value = ct.AXIS_LABELS[1]
        assert c.apply_slice_fields() is None
        assert (c.slice_axis, c.slice_index) == (1, 7)

    def test_compute_handler_adds_the_result(self, tmp_path):
        c = with_widgets(loaded(tmp_path, p_orbital_volume(), p_orbital_volume()))
        c.op_a_selection = SimpleNamespace(value="cube0")
        c.op_selection = SimpleNamespace(value=ct.OP_SCALE)
        c.op_b_selection = SimpleNamespace(value=None)
        c.factor_input = SimpleNamespace(value="2")
        c.result_name_input = SimpleNamespace(value="")
        errors = []
        c._show_error = lambda title, message: errors.append(message) or asyncio.sleep(0)
        c.sync_widgets = lambda **kwargs: None
        asyncio.run(c.on_compute(None))
        assert errors == [] and c.names()[-1] == "2 × cube0"
        assert "2 × cube0" in c.status_label.text



# ----------------------------------------------------------------------
# Flicker: on WinForms every Toga widget write re-lays out the whole viewer
# window (~70 ms), so the tab must write only what changed, and handlers
# must not re-run for the tab's own programmatic updates.
# ----------------------------------------------------------------------
class Recorder:
    """A stub widget that counts every attribute write after construction."""

    def __init__(self, **values):
        object.__setattr__(self, "writes", [])
        for name, value in values.items():
            object.__setattr__(self, name, value)

    def __setattr__(self, name, value):
        self.writes.append(name)
        object.__setattr__(self, name, value)


def recorded_widgets(c):
    c.cube_table = Recorder(data=[])
    for name in ("iso_cube_selection", "slice_cube_selection", "op_a_selection", "op_b_selection"):
        setattr(c, name, Recorder(items=[], value=None))
    c.isovalue_input = Recorder(value="")
    c.show_selection = Recorder(value=ct.SHOW_BOTH)
    c.positive_colour_input = Recorder(value="Blue")
    c.negative_colour_input = Recorder(value="Red")
    c.opacity_slider = Recorder(value=1.0)
    c.visible_switch = Recorder(value=True)
    c.slice_axis_selection = Recorder(value=ct.AXIS_LABELS[2])
    c.slice_range_selection = Recorder(value=list(ct.RANGE_LABELS)[0])
    c.slice_switch = Recorder(value=False)
    c.slice_slider = Recorder(value=0.0, max=1, tick_count=None, enabled=False)
    c.slice_index_label = Recorder(text="No cube")
    c.status_label = Recorder(text=" ")
    return c


def all_writes(c):
    widgets = [w for w in vars(c).values() if isinstance(w, Recorder)]
    return sum(len(w.writes) for w in widgets)


def clear_writes(c):
    for w in vars(c).values():
        if isinstance(w, Recorder):
            object.__setattr__(w, "writes", [])


class TestNoFlicker:
    def test_a_sync_with_nothing_changed_writes_no_widget(self, tmp_path):
        c = recorded_widgets(loaded(tmp_path, p_orbital_volume(), p_orbital_volume()))
        c.sync_widgets(selected="cube1")
        assert all_writes(c) > 0  # the first sync fills the widgets
        clear_writes(c)
        c.sync_widgets(selected="cube1")
        assert all_writes(c) == 0

    def test_only_the_changed_field_is_written(self, tmp_path):
        c = recorded_widgets(loaded(tmp_path, p_orbital_volume()))
        c.sync_widgets(selected="cube0")
        clear_writes(c)
        c.entries[0].isovalue = 0.123
        c.sync_widgets()
        assert c.isovalue_input.writes == ["value"] and all_writes(c) == 1

    def test_repeating_a_status_message_writes_nothing(self, tmp_path):
        c = recorded_widgets(loaded(tmp_path, p_orbital_volume()))
        c.set_status("Reading x.cube...")
        c.set_status("Reading x.cube...")
        assert c.status_label.writes == ["text"]

    @pytest.mark.parametrize("handler", ["on_isosurface_change", "on_slice_change", "on_show_volumes_change"])
    def test_handlers_ignore_the_tabs_own_updates(self, tmp_path, handler):
        """Toga runs an async handler as a task after the write that fired it has
        returned - after _syncing is cleared - so the guard must be checked at once."""
        c = with_widgets(loaded(tmp_path, p_orbital_volume()))
        started = []
        c._start = lambda coroutine: started.append(coroutine.__name__) or coroutine.close()
        c._syncing = True
        getattr(c, handler)(SimpleNamespace(value=True))
        assert started == []
        c._syncing = False
        getattr(c, handler)(SimpleNamespace(value=True))
        assert len(started) == 1
