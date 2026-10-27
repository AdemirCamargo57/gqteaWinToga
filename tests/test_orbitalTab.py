"""Orbitals tab of the 3D Molecular Viewer: the controller (headless, stub widgets) and its rules."""
import asyncio
import inspect
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import molecularOrbitals as mo
import orbitalTab
import volumeData as vd
from cubeRender import CubeSceneLayer
from help import HelpGqteaWin
from orbitalTab import OrbitalTabController

SRC = Path(__file__).resolve().parent.parent
DATA = Path(__file__).resolve().parent / "data"
WATER = DATA / "fchk" / "water_631gd_5d.fchk"
CH3_U = DATA / "fchk" / "ch3_ub3lyp_631gd.fchk"
CH3_RO = DATA / "orbitals" / "ch3_rob3lyp_631gd.fchk"


def W(**kwargs):
    return SimpleNamespace(**kwargs)


@pytest.fixture(autouse=True)
def plain_dialogs(monkeypatch):
    """Dialogs as plain records, so tests can read what the user would be told."""
    for name in ("ErrorDialog", "InfoDialog", "QuestionDialog"):
        monkeypatch.setattr(orbitalTab.toga, name, lambda title, message, _k=name: W(kind=_k, title=title,
                                                                                        message=message))


class StubWindow:
    def __init__(self, answer=True):
        self.dialogs, self.answer = [], answer

    async def dialog(self, d):
        self.dialogs.append(d)
        return self.answer


class StubCubeController:
    def __init__(self):
        self.added, self.synced = [], []

    def add_result(self, volume, name):
        entry = SimpleNamespace(name=name, volume=volume, visible=True)
        self.added.append(entry)
        return entry

    def sync_widgets(self, selected=None):
        self.synced.append(selected)


class StubViewer:
    def __init__(self):
        self.orbital_layer = CubeSceneLayer()
        self.cube_controller = StubCubeController()
        self.frames, self.loaded, self.refreshed, self.opened = [], [], 0, 0
        self.main_window = StubWindow()
        self.window_open, self.frame = True, np.zeros((4, 6, 3), dtype=np.uint8)
        self.claims = []
        # The canvas surface style, as MolecularViewerUI keeps it (Solid on every layer at start).
        self.canvas_surface_style, self.canvas_transparency = vd.SURFACE_SOLID, vd.DEFAULT_TRANSPARENCY
        self.orbital_layer.set_surface_style(vd.SURFACE_SOLID, vd.DEFAULT_TRANSPARENCY)
        self.orbital_controller = None
        self.style_calls = []

    def apply_surface_style(self, style=None, transparency=None):
        self.style_calls.append({"style": style, "transparency": transparency})
        self.canvas_surface_style = style or self.canvas_surface_style
        self.canvas_transparency = self.canvas_transparency if transparency is None else transparency
        self.orbital_layer.set_surface_style(self.canvas_surface_style, self.canvas_transparency)
        if self.orbital_controller is not None:
            self.orbital_controller.sync_surface_widgets(self.canvas_surface_style, self.canvas_transparency)

    def show_analysis_layer(self, name):
        self.claims.append(name)

    def viewer_window_is_open(self):
        return self.window_open

    async def capture_frame(self):
        return self.frame

    def get_current_frame_data(self):
        return self.loaded[-1] if self.loaded else []

    def load_designed_molecule(self, frame):
        self.loaded.append(frame)
        self.frames = [frame]

    def invalidate_scene_cache(self):
        self.refreshed += 1

    def request_projection_update(self):
        pass

    async def ensure_viewer_window(self):
        self.opened += 1
        return True, ""


BUTTONS = ("load_button", "clear_button", "display_button", "homo_button", "lumo_button", "previous_button",
           "next_button", "levels_button", "export_cube_button", "send_button", "export_range_button",
           "export_table_button", "export_3d_button")


def make_controller(path=WATER, output_dir=""):
    c = OrbitalTabController(viewer=StubViewer())
    c.fchk_input = W(value=str(path))
    c.output_dir_input = W(value=str(output_dir))
    for button in BUTTONS:
        setattr(c, button, W(enabled=True))
    c.status_label, c.summary_label, c.selected_label = W(text="Ready."), W(text=""), W(text="")
    c.progress_bar = W(value=0)
    c.spin_selection = W(items=[orbitalTab.ONE_SET_ITEM], value=orbitalTab.ONE_SET_ITEM, enabled=False)
    c.orbital_table = W(data=[], selection=None)
    c.surface_selection = W(value=orbitalTab.SURFACE_PHASES)
    c.isovalue_label = W(text=orbitalTab.ISOVALUE_LABELS[orbitalTab.SURFACE_PHASES])
    c.isovalue_input = W(value="", placeholder="0.02")
    c.grid_selection = W(value="Coarse")   # fast tests
    c.spacing_input = W(value="", enabled=False)
    c.style_selection = W(value=orbitalTab.STYLE_SOLID)
    c.transparency_slider = W(value=vd.DEFAULT_TRANSPARENCY, enabled=False)
    c.viewer.orbital_controller = c
    c.colour_inputs = {sign: W(value=colour) for sign, colour in orbitalTab.DEFAULT_COLOURS.items()}
    c.range_input = W(value="")
    c.report_output = W(value="")
    return c


async def settle(c):
    while c._tasks:
        await asyncio.gather(*list(c._tasks))


@pytest.fixture
def loaded(tmp_path):
    c = make_controller(output_dir=tmp_path)
    asyncio.run(c.load_file(None))
    return c


@pytest.fixture
def shown(loaded):
    asyncio.run(loaded.display_orbital(None))
    return loaded


def act(c, handler, *args):
    """Call a synchronous widget handler on a running loop (as Toga does) and wait for what it starts."""
    async def go():
        handler(*args)
        await settle(c)
    asyncio.run(go())


class TestInputs:
    def test_every_orbital_volumes_parameter_is_supplied_by_the_ui(self, loaded):
        wanted = set(inspect.signature(mo.orbital_volumes).parameters) - {"catalogue", "progress", "use_numba"}
        assert set(loaded._collect_display_inputs()) == wanted

    def test_display_inputs(self, loaded):
        assert loaded._collect_display_inputs() == {"spin": "alpha", "numbers": [5], "spacing": pytest.approx(1 / 3),
                                                    "isovalue": 0.02}
        loaded.grid_selection.value = "Fine"
        assert loaded._collect_display_inputs()["spacing"] == pytest.approx(1 / 12)
        loaded.grid_selection.value = orbitalTab.GRID_CUSTOM
        assert loaded._collect_display_inputs()["spacing"] == pytest.approx(1 / 6)  # blank: Medium's
        loaded.spacing_input.value = "0.3"
        assert loaded._collect_display_inputs()["spacing"] == 0.3
        loaded.spacing_input.value = "-1"
        with pytest.raises(ValueError, match="grid spacing"):
            loaded._collect_display_inputs()

    def test_density_isovalue_becomes_the_amplitude_grid_isovalue(self, loaded):
        loaded.on_surface_change(W(value=orbitalTab.SURFACE_DENSITY))
        assert loaded._collect_display_inputs()["isovalue"] == pytest.approx(0.02)  # sqrt(0.0004)
        loaded.isovalue_input.value = "0.0009"
        assert loaded._collect_display_inputs()["isovalue"] == pytest.approx(0.03)

    def test_inputs_need_a_file(self):
        c = make_controller()
        with pytest.raises(ValueError, match="Load a .fchk"):
            c._collect_display_inputs()
        c.fchk_input.value = ""
        with pytest.raises(ValueError, match="select a .fchk"):
            c._fchk_path()
        c.fchk_input.value = str(DATA / "none.fchk")
        with pytest.raises(ValueError, match="File not found"):
            c._fchk_path()


class TestLoad:
    def test_load_lists_every_orbital_and_selects_the_homo(self, loaded):
        assert len(loaded.orbital_table.data) == 18
        assert loaded.orbital_table.data[4] == ("5", "HOMO", "2", "-0.497563", "-13.5394")
        assert loaded.orbital_table.data[5][:3] == ("6", "LUMO", "0")
        assert loaded.number == 5 and loaded.selected_label.text.startswith("Selected: MO 5 (HOMO), occupancy 2")
        assert loaded.summary_label.text.startswith("RHF/6-31G(d), 10 electrons, 18 MOs (restricted).")
        assert loaded.spin_selection.items == [orbitalTab.ONE_SET_ITEM] and not loaded.spin_selection.enabled
        assert "18 MOs read (restricted)" in loaded.status_label.text

    def test_unrestricted_files_offer_both_spins(self, tmp_path):
        c = make_controller(CH3_U, tmp_path)
        asyncio.run(c.load_file(None))
        assert c.spin_selection.items == ["Alpha", "Beta"] and c.spin_selection.enabled
        assert c.number == 5 and c.orbital_table.data[4][1:3] == ("HOMO", "1")
        c.on_spin_change(W(value="Beta"))
        assert c.spin == "beta" and c.number == 5           # the MO number is kept
        assert c.orbital_table.data[4][1:3] == ("LUMO", "0")
        assert "beta" in c.selected_label.text
        assert c._collect_display_inputs()["spin"] == "beta"

    def test_restricted_open_occupancies(self, tmp_path):
        c = make_controller(CH3_RO, tmp_path)
        asyncio.run(c.load_file(None))
        assert [row[2] for row in c.orbital_table.data[:6]] == ["2", "2", "2", "2", "1", "0"]
        assert "(restricted-open)" in c.summary_label.text

    def test_a_bad_file_shows_a_dialog(self, tmp_path):
        bad = tmp_path / "bad.fchk"
        bad.write_text("not a checkpoint\n")
        c = make_controller(bad)
        asyncio.run(c.load_file(None))
        assert c.catalogue is None and c.viewer.main_window.dialogs[0].title == "Orbitals Input Error"
        assert c.status_label.text == "The file was not loaded." and all(getattr(c, b).enabled for b in BUTTONS)


class TestSelection:
    def test_buttons_move_the_selection_and_clamp(self, loaded):
        loaded.go_lumo(None)
        assert loaded.number == 6
        loaded.go_next(None)
        loaded.go_next(None)
        assert loaded.number == 8
        loaded.go_previous(None)
        assert loaded.number == 7
        loaded.go_homo(None)
        assert loaded.number == 5
        loaded.select_orbital(99)
        assert loaded.number == 18
        loaded.go_next(None)
        assert loaded.number == 18
        loaded.select_orbital(-3)
        assert loaded.number == 1 and "MO 1 (HOMO-4)" in loaded.selected_label.text
        assert loaded._tasks == set()      # nothing displayed: no redraw

    def test_clicking_a_row_selects_it(self, loaded):
        loaded.on_orbital_select(W(selection=W(mo="3")))
        assert loaded.number == 3 and "MO 3 (HOMO-2)" in loaded.selected_label.text
        loaded._syncing = True
        loaded.on_orbital_select(W(selection=W(mo="7")))   # the tab's own write: ignored
        assert loaded.number == 3
        loaded._syncing = False
        loaded.on_orbital_select(W(selection=None))
        assert loaded.number == 3

    def test_buttons_before_loading_say_so(self):
        c = make_controller()
        act(c, c.go_next, None)
        assert c.viewer.main_window.dialogs[0].message == "Load a .fchk file first."


class TestDisplay:
    def test_display_loads_the_molecule_and_two_phases(self, shown):
        v = shown.viewer
        assert len(v.loaded) == 1 and len(v.loaded[0]) == 3 and v.opened == 1
        assert v.claims == ["orbital_layer"]
        items = v.orbital_layer.items()
        assert len(items) == 2 and all(i.kind == "surface" and i.opaque for i in items)
        assert np.allclose([i.rgba[0, :3] for i in items], [vd.parse_colour("Red"), vd.parse_colour("Green")])
        assert "2 isosurfaces shown. MO 5 (HOMO)" in shown.status_label.text and shown.progress_bar.value == 100
        asyncio.run(shown.display_orbital(None))   # same molecule: no second load, no question
        assert len(v.loaded) == 1 and v.main_window.dialogs == []

    def test_a_different_structure_is_asked_about(self, loaded):
        loaded.viewer.load_designed_molecule([("C", (0.0, 0.0, 0.0))])
        loaded.viewer.main_window.answer = False
        asyncio.run(loaded.display_orbital(None))
        assert loaded.viewer.main_window.dialogs[0].kind == "QuestionDialog"
        assert len(loaded.viewer.loaded) == 1   # kept the user's structure

    def test_once_displayed_every_selection_redraws(self, shown):
        act(shown, shown.go_next, None)
        assert shown._current[2] == 6 and "MO 6 (LUMO)" in shown.status_label.text
        act(shown, shown.on_orbital_select, W(selection=W(mo="2")))
        assert shown._current[2] == 2

    def test_fast_presses_are_coalesced_into_the_last_selection(self, shown, monkeypatch):
        computed = []
        original = shown.orbital_entries

        def spy(**kwargs):
            computed.append(kwargs["numbers"][0])
            return original(**kwargs)
        monkeypatch.setattr(shown, "orbital_entries", spy)

        async def go():
            for _ in range(4):
                shown.go_next(None)
            await settle(shown)
        asyncio.run(go())
        assert shown.number == 9 and shown._current[2] == 9
        assert computed[-1] == 9 and len(computed) <= 2   # the first request, then the latest one

    def test_style_goes_to_the_canvas_and_never_recomputes(self, shown):
        entry = shown._grids[shown._current]["phases"]
        built = shown.viewer.orbital_layer.items()
        act(shown, shown.on_style_change, W(value=orbitalTab.STYLE_MESH))
        assert shown.viewer.style_calls == [{"style": orbitalTab.STYLE_MESH, "transparency": None}]
        assert shown.viewer.orbital_layer.items() is built                # the tab built nothing new
        lines = shown.viewer.orbital_layer.drawn_items()                  # the layer draws contour lines
        assert [i.kind for i in lines] == ["lines", "lines"] and all(i.indices.size > 0 for i in lines)
        assert np.allclose([i.rgba[0, :3] for i in lines], [vd.parse_colour("Red"), vd.parse_colour("Green")])
        assert shown.transparency_slider.enabled is False and shown.style_selection.value == orbitalTab.STYLE_MESH
        act(shown, shown.on_style_change, W(value=orbitalTab.STYLE_TRANSPARENT))
        assert shown.transparency_slider.enabled is True
        assert all(i.rgba[0, 3] == pytest.approx(0.5) for i in shown.viewer.orbital_layer.drawn_items())
        act(shown, shown.on_transparency_change, W(value=0.7))
        assert shown.transparency_slider.value == pytest.approx(0.7)
        assert all(i.rgba[0, 3] == pytest.approx(0.3) for i in shown.viewer.orbital_layer.drawn_items())
        assert all(i.rgba[0, 3] == pytest.approx(1.0) for i in shown.viewer.orbital_layer.items())
        shown.colour_inputs[1].value = "blue"
        act(shown, shown.on_colour_confirm, shown.colour_inputs[1])
        assert np.allclose(shown.viewer.orbital_layer.items()[0].rgba[0, :3], vd.parse_colour("Blue"))
        assert shown._grids[shown._current]["phases"] is entry and len(shown._grids) == 1

    def test_the_display_tab_s_style_shows_here(self, loaded):
        loaded.sync_surface_widgets(orbitalTab.STYLE_TRANSPARENT, 0.25)
        assert loaded.style_selection.value == orbitalTab.STYLE_TRANSPARENT
        assert loaded.transparency_slider.value == 0.25 and loaded.transparency_slider.enabled is True
        assert loaded.viewer.style_calls == []        # showing it is not a change: nothing is sent back
        loaded._syncing = True
        loaded.on_style_change(W(value=orbitalTab.STYLE_MESH))   # the tab's own write: ignored
        loaded._syncing = False
        assert loaded.viewer.style_calls == []

    def test_without_a_canvas_style_the_tab_keeps_its_own(self):
        c = make_controller()
        c.viewer.apply_surface_style = None        # e.g. a viewer without the Display-tab setting
        c.on_style_change(W(value=orbitalTab.STYLE_TRANSPARENT))
        assert c.style == orbitalTab.STYLE_TRANSPARENT and c.transparency_slider.enabled is True

    def test_a_bad_colour_is_refused_and_put_back(self, shown):
        shown.colour_inputs[-1].value = "nonsense"
        act(shown, shown.on_colour_confirm, shown.colour_inputs[-1])
        assert shown.colour_inputs[-1].value == "Green" and shown.viewer.main_window.dialogs[0].title == "Colour"

    def test_density_draws_one_surface_at_its_own_isovalue(self, shown):
        act(shown, shown.on_surface_change, W(value=orbitalTab.SURFACE_DENSITY))
        items = shown.viewer.orbital_layer.items()
        assert len(items) == 1 and np.allclose(items[0].rgba[0, :3], vd.parse_colour("Red"))
        density = shown._grids[shown._current]["density"]
        assert density.isovalue == pytest.approx(0.0004)
        assert np.array_equal(density.volume.values, shown._grids[shown._current]["phases"].volume.values ** 2)
        assert shown.isovalue_label.text == orbitalTab.ISOVALUE_LABELS[orbitalTab.SURFACE_DENSITY]
        assert shown.isovalue_input.placeholder == "0.0004"
        assert len(shown._grids) == 1          # same grid: sqrt(0.0004) = 0.02
        shown.isovalue_input.value = "0.0001"
        act(shown, shown.on_surface_change, W(value=orbitalTab.SURFACE_PHASES))
        assert shown.isovalue_input.value == ""   # each surface keeps its own typed isovalue
        act(shown, shown.on_surface_change, W(value=orbitalTab.SURFACE_DENSITY))
        assert shown.isovalue_input.value == "0.0001"

    def test_the_last_three_grids_are_kept(self, shown):
        for number in (6, 7, 8):
            act(shown, shown.select_orbital, number)
        assert [key[2] for key in shown._grids] == [6, 7, 8]
        act(shown, shown.select_orbital, 7)
        assert [key[2] for key in shown._grids] == [6, 8, 7]   # reused, not recomputed

    def test_a_grid_that_is_too_large_is_refused(self, loaded):
        loaded.grid_selection.value = orbitalTab.GRID_CUSTOM
        loaded.spacing_input.value = "0.01"
        asyncio.run(loaded.display_orbital(None))
        dialog = loaded.viewer.main_window.dialogs[0]
        assert dialog.title == "Orbital Grid Error" and "million points" in dialog.message
        assert loaded.viewer.orbital_layer.items() == () and loaded.viewer.opened == 0

    def test_another_tab_taking_the_canvas_stops_the_redraws(self, shown):
        shown.viewer.orbital_layer.clear()      # what show_analysis_layer does for another tab
        shown.go_next(None)
        assert shown._tasks == set()

    def test_loading_another_file_while_displayed_shows_its_homo(self, shown):
        shown.fchk_input.value = str(CH3_RO)
        asyncio.run(shown.load_file(None))
        assert shown.catalogue.kind == mo.RESTRICTED_OPEN and shown._current[2] == 5
        assert shown.viewer.main_window.dialogs[0].kind == "QuestionDialog"   # another molecule


class TestExports:
    def test_export_cube(self, loaded, tmp_path):
        asyncio.run(loaded.export_cube(None))
        path = tmp_path / "water_631gd_5d_MO5.cube"
        assert path.is_file() and str(path) in loaded.status_label.text
        loaded.on_surface_change(W(value=orbitalTab.SURFACE_DENSITY))
        asyncio.run(loaded.export_cube(None))
        assert (tmp_path / "water_631gd_5d_MO5_density.cube").is_file()

    def test_unrestricted_cube_names_carry_the_spin(self, tmp_path):
        c = make_controller(CH3_U, tmp_path)
        asyncio.run(c.load_file(None))
        c.on_spin_change(W(value="Beta"))
        asyncio.run(c.export_cube(None))
        assert (tmp_path / "ch3_ub3lyp_631gd_MO5_beta.cube").is_file()

    def test_export_range(self, loaded, tmp_path):
        loaded.range_input.value = "HOMO-1:LUMO"
        asyncio.run(loaded.export_range(None))
        names = sorted(p.name for p in tmp_path.glob("*.cube"))
        assert names == ["water_631gd_5d_MO4.cube", "water_631gd_5d_MO5.cube", "water_631gd_5d_MO6.cube"]
        assert "3 cube files written" in loaded.status_label.text and loaded.progress_bar.value == 100

    def test_a_bad_range_shows_a_dialog(self, loaded, tmp_path):
        loaded.range_input.value = "SOMO"
        asyncio.run(loaded.export_range(None))
        assert loaded.viewer.main_window.dialogs[0].title == "Export Input Error"
        assert list(tmp_path.glob("*.cube")) == []

    def test_send_to_cube_tab(self, loaded):
        asyncio.run(loaded.send_to_cube_tab(None))
        assert loaded.viewer.main_window.dialogs[0].title == "Nothing to Send"
        asyncio.run(loaded.display_orbital(None))
        asyncio.run(loaded.send_to_cube_tab(None))
        added = loaded.viewer.cube_controller.added[0]
        assert added.name == "water_631gd_5d_MO5" and added.visible is False
        assert loaded.viewer.cube_controller.synced == ["water_631gd_5d_MO5"]

    def test_send_builds_an_unvisited_cube_tab(self, loaded):
        cube, built = loaded.viewer.cube_controller, []
        loaded.viewer.cube_controller = None

        def ensure_tab_built(title):
            built.append(title)
            loaded.viewer.cube_controller = cube
            return cube
        loaded.viewer.ensure_tab_built = ensure_tab_built
        asyncio.run(loaded.display_orbital(None))
        asyncio.run(loaded.send_to_cube_tab(None))
        assert built == ["Cube"] and cube.added[0].name == "water_631gd_5d_MO5"

    def test_export_3d_image(self, shown, tmp_path):
        async def ask(title, suggested, file_types):
            assert suggested == "water_631gd_5d_MO5.png"
            return tmp_path / "view.png"
        shown._ask_save_path = ask
        asyncio.run(shown.export_3d_image(None))
        assert (tmp_path / "view.png").is_file() and "6 x 4 px" in shown.status_label.text

    def test_export_3d_image_needs_the_window(self, loaded):
        loaded.viewer.window_open = False
        asyncio.run(loaded.export_3d_image(None))
        assert loaded.viewer.main_window.dialogs[0].title == "No Image"

    def test_export_energy_table(self, loaded, tmp_path):
        asyncio.run(loaded.export_energy_table(None))
        path = tmp_path / "water_631gd_5d_MO_energies.txt"
        assert path.is_file() and loaded.report_output.value == path.read_text(encoding="utf-8")

    def test_energy_levels_open_the_plot_viewer(self, loaded, tmp_path, monkeypatch):
        import displayPlots
        launched = []
        monkeypatch.setattr(displayPlots, "launch_plot_viewer", lambda path: launched.append(path) or True)
        asyncio.run(loaded.show_energy_levels(None))
        assert launched == [str(tmp_path / "water_631gd_5d_MO_levels.json")]
        monkeypatch.setattr(displayPlots, "launch_plot_viewer", lambda path: False)
        asyncio.run(loaded.show_energy_levels(None))
        assert loaded.viewer.main_window.dialogs[0].title == "Energy Levels Saved"

    def test_exports_need_a_file(self):
        c = make_controller()
        for handler in (c.export_cube, c.export_range, c.export_energy_table, c.show_energy_levels):
            asyncio.run(handler(None))
        assert [d.message for d in c.viewer.main_window.dialogs] == ["Load a .fchk file first."] * 4


class TestClear:
    def test_clear_resets_the_tab(self, shown):
        shown.clear_output(None)
        assert shown.catalogue is None and shown.number is None and shown._grids == {}
        assert shown.orbital_table.data == [] and shown.viewer.orbital_layer.items() == ()
        assert shown.summary_label.text == "No file loaded." and shown.report_output.value == HelpGqteaWin.help_orbitals


def test_the_tab_does_not_import_the_viewer_at_module_level():
    code = "import sys, orbitalTab; assert 'molecularViewer' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], cwd=SRC, check=True)


def test_select_table_row_never_raises():
    orbitalTab.select_table_row(W(), 3)                                  # no native control, no scroll
    rows = []
    orbitalTab.select_table_row(W(scroll_to_row=rows.append), 3)
    assert rows == [3]
    broken = W(_impl=W(native=W()), scroll_to_row=lambda i: 1 / 0)
    orbitalTab.select_table_row(broken, 3)
