"""QTAIM tab of the 3D Molecular Viewer: scene layer, controller (headless) and integration."""
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import qtaimAnalysis as qa
import qtaimLayer as ql
import qtaimTopology as qt
from qtaim_fixtures import fchk

SRC = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def water_results():
    return qa.QTAIMAnalyser().run(str(fchk("water_631gd_6d.fchk")))


@pytest.fixture(scope="module")
def cage_results():
    return qa.QTAIMAnalyser().run(str(fchk("bicyclopentane_631gd.fchk")))


class RecordingViewer:
    def __init__(self):
        self.spheres, self.texts = [], []

    def draw_sphere(self, radius):
        self.spheres.append(radius)

    def _draw_text_3d(self, position, text, color=(1, 1, 0)):
        self.texts.append((text, color))


@pytest.fixture
def no_gl(monkeypatch):
    calls = []
    for name in ("glPushMatrix", "glPopMatrix", "glTranslatef", "glMaterialfv", "glDisable", "glEnable",
                 "glLineWidth", "glColor3f", "glBegin", "glEnd", "glVertex3f", "glDepthRange"):
        monkeypatch.setattr(ql, name, lambda *a, _n=name: calls.append(_n))
    return calls


class TestSceneLayer:
    def test_scene_from_results_is_in_angstrom(self, water_results):
        scene = ql.scene_from_results(water_results)
        assert len(scene.cps) == 5 and len(scene.bond_paths) == 4 and scene.ring_paths == ()
        cp = water_results.critical_points[3]
        assert np.allclose(scene.cps[3].position, cp.position * qt.BOHR_RADIUS_ANGSTROM)
        assert scene.bond_paths[0].dtype == np.float32

    def test_visible_cps_follow_the_type_switches(self, cage_results):
        layer = ql.QTAIMSceneLayer()
        layer.set_scene(ql.scene_from_results(cage_results))
        assert len(layer.visible_cps()) == 31
        layer.set_display_options(show_types={qt.NCP: False})
        assert {c.kind for c in layer.visible_cps()} == {qt.BCP, qt.RCP, qt.CCP}
        layer.set_display_options(show_types={qt.BCP: False, qt.RCP: False, qt.CCP: False})
        assert layer.visible_cps() == ()
        with pytest.raises(ValueError, match="Unknown critical point type"):
            layer.set_display_options(show_types={"(3,0)": True})

    def test_radius(self):
        layer = ql.QTAIMSceneLayer()
        assert layer.radius() == 0.0
        layer.set_scene(ql.QTAIMScene((ql.SceneCP(1, qt.NCP, (1.0, 0.0, 0.0)),),
                                      (np.array([[0, 0, 0], [0, 3.0, 0]], dtype=np.float32),), (), (0.0, 0.0, 0.0)))
        assert layer.radius() == pytest.approx(3.0)

    def test_draw_scene_spheres_highlight_labels_and_paths(self, water_results, no_gl):
        layer = ql.QTAIMSceneLayer()
        layer.set_scene(ql.scene_from_results(water_results))
        layer.set_display_options(highlight=4)
        viewer = RecordingViewer()
        layer.draw_scene(viewer)
        assert viewer.spheres.count(ql.HIGHLIGHT_RADIUS) == 1 and viewer.spheres.count(ql.CP_RADIUS) == 4
        assert [t for t, _ in viewer.texts] == ["1", "2", "3", "4", "5"]
        assert no_gl.count("glBegin") == 4 and no_gl.count("glEnd") == 4
        last_end = len(no_gl) - 1 - no_gl[::-1].index("glEnd")
        assert "glEnable" in no_gl[last_end:]   # lighting restored after the paths

    def test_draw_scene_honours_the_switches(self, cage_results, no_gl):
        layer = ql.QTAIMSceneLayer()
        layer.set_scene(ql.scene_from_results(cage_results))
        layer.set_display_options(bond_paths=False, labels=False, show_types={qt.NCP: False})
        viewer = RecordingViewer()
        layer.draw_scene(viewer)
        assert no_gl.count("glBegin") == 3      # the three ring paths only
        assert viewer.texts == [] and len(viewer.spheres) == 18

    def test_empty_layer_draws_nothing(self, no_gl):
        viewer = RecordingViewer()
        ql.QTAIMSceneLayer().draw_scene(viewer)
        assert viewer.spheres == [] and no_gl == []

    def test_layer_does_not_import_toga_or_the_viewer(self):
        code = "import sys, qtaimLayer; print('toga' in sys.modules, 'molecularViewer' in sys.modules)"
        out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
        assert out.stdout.strip() == "False False"


# --------------------------------------------------------------------------- #
# The controller, headless: stub widgets and a stub viewer                      #
# --------------------------------------------------------------------------- #
import asyncio
import inspect
import json

import qtaimTab
from qtaimTab import QTAIMTabController


def W(**kwargs):
    return SimpleNamespace(**kwargs)


class StubWindow:
    def __init__(self, answer=True):
        self.dialogs, self.answer = [], answer

    async def dialog(self, d):
        self.dialogs.append(d)
        return self.answer


class StubViewer:
    def __init__(self):
        self.qtaim_layer = ql.QTAIMSceneLayer()
        self.frames, self.loaded, self.refreshed, self.opened = [], [], 0, 0
        self.main_window = StubWindow()
        self.background_color = (0.0, 0.0, 0.0)
        self.atom_display_style, self.show_atom_symbols, self.show_atom_numbers = "Line style", False, False
        self.window_open, self.frame = True, np.zeros((4, 6, 3), dtype=np.uint8)

    def set_background_color(self, rgb):
        self.background_color = tuple(rgb)

    def apply_view_settings(self, style=None, symbols=None, numbers=None):
        if style is not None:
            self.atom_display_style = style
        if symbols is not None:
            self.show_atom_symbols = symbols
        if numbers is not None:
            self.show_atom_numbers = numbers

    def viewer_window_is_open(self):
        return self.window_open

    async def capture_frame(self):
        if not self.window_open:
            raise RuntimeError("closed")
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

    def get_system_center(self):
        return (0.0, 0.0, 0.0)

    async def ensure_viewer_window(self):
        self.opened += 1
        return True, ""


def make_controller(fchk_path="", output_dir=""):
    c = QTAIMTabController(viewer=StubViewer())
    c.fchk_input, c.output_dir_input = W(value=str(fchk_path)), W(value=str(output_dir))
    c.mode_switches = {m: W(value=True) for m in qt.MODE_LABELS}
    c.search_inputs = {name: W(value="") for name, _ in qtaimTab.SEARCH_FIELDS}
    c.path_inputs = {name: W(value="") for name, _ in qtaimTab.PATH_FIELDS}
    for name in ("run_button", "display_button", "clear_button", "map_button"):
        setattr(c, name, W(enabled=True))
    c.status_label, c.summary_label = W(text="Ready."), W(text="No results.")
    c.progress_bar, c.cp_table = W(value=0), W(data=[], selection=None)
    c.report_output = W(value="")
    c.type_switches = {t: W(value=True) for t in qt.CP_TYPES}
    c.bond_path_switch, c.ring_path_switch, c.label_switch = W(value=True), W(value=True), W(value=True)
    c.plane_selection = W(value=qtaimTab.PLANE_THREE_ATOMS)
    c.plane_atoms_input, c.plane_offset_input = W(value=""), W(value="")
    c.margin_input, c.grid_points_input = W(value=""), W(value="30")
    c.contour_switch, c.gradient_switch, c.interbasin_switch = W(value=True), W(value=True), W(value=True)
    c.background_selection, c.background_input = W(value="Black"), W(value="")
    c.style_selection = W(value="Line style")
    c.symbols_switch, c.numbers_switch = W(value=False), W(value=False)
    c.export_3d_button, c.export_2d_button = W(enabled=True), W(enabled=True)
    return c


@pytest.fixture
def ran(tmp_path):
    c = make_controller(fchk("water_631gd_6d.fchk"), tmp_path)
    asyncio.run(c.workflow(None))
    return c


class TestInputs:
    def test_every_run_parameter_is_supplied_by_the_ui(self):
        c = make_controller(fchk("water_631gd_6d.fchk"))
        wanted = set(inspect.signature(qa.QTAIMAnalyser.run).parameters) - {"self", "use_numba", "progress"}
        assert set(c._collect_inputs()) == wanted

    def test_every_plane_map_parameter_is_supplied_by_the_ui(self, water_results):
        c = make_controller()
        c.results = water_results
        wanted = set(inspect.signature(qa.QTAIMAnalyser.plane_map).parameters) - {"self", "results", "use_numba",
                                                                                   "progress"}
        assert set(c._collect_map_inputs()) == wanted

    def test_blank_advanced_fields_are_multiwfn_defaults(self):
        params = make_controller(fchk("water_631gd_6d.fchk"))._collect_inputs()
        assert params["search"] == qt.CPSearchParameters() and params["paths"] == qt.PathParameters()
        assert params["modes"] == (2, 3, 4, 5)

    def test_typed_values_are_passed_through(self):
        c = make_controller(fchk("water_631gd_6d.fchk"))
        c.search_inputs["trust_radius"].value = "0.1"
        c.search_inputs["max_iterations"].value = "200"
        c.path_inputs["step"].value = "0,02"
        c.mode_switches[5].value = False
        params = c._collect_inputs()
        assert params["search"].trust_radius == 0.1 and params["search"].max_iterations == 200
        assert params["paths"].step == 0.02 and params["modes"] == (2, 3, 4)

    @pytest.mark.parametrize("field, text, message", [("max_iterations", "1.5", "Maximum iterations"),
                                                      ("step_scale", "-1", "Step-size scale"),
                                                      ("min_cp_distance", "abc", "Minimum distance")])
    def test_bad_values_name_the_field(self, field, text, message):
        c = make_controller(fchk("water_631gd_6d.fchk"))
        c.search_inputs[field].value = text
        with pytest.raises(ValueError, match=message):
            c._collect_inputs()

    def test_no_mode_and_no_file(self, tmp_path):
        c = make_controller(fchk("water_631gd_6d.fchk"))
        for s in c.mode_switches.values():
            s.value = False
        with pytest.raises(ValueError, match="at least one CP search mode"):
            c._collect_inputs()
        with pytest.raises(ValueError, match="Please select"):
            make_controller()._collect_inputs()
        with pytest.raises(ValueError, match="File not found"):
            make_controller(tmp_path / "x.fchk")._collect_inputs()

    def test_parse_atom_triple(self):
        assert qtaimTab.parse_atom_triple("1, 2, 3", 3) == (0, 1, 2)
        assert qtaimTab.parse_atom_triple("3 1 2", 3) == (2, 0, 1)
        for bad, message in (("1 1 2", "different"), ("1 2", "three"), ("1 2 9", "outside 1-3"),
                             ("a b c", "atom numbers")):
            with pytest.raises(ValueError, match=message):
                qtaimTab.parse_atom_triple(bad, 3)

    def test_map_inputs(self, water_results):
        c = make_controller()
        c.results = water_results
        plane = c._collect_map_inputs()["plane"]
        assert plane.label == "plane of O1, H2, H3"
        c.plane_selection.value, c.plane_offset_input.value = "XZ", "1.0"
        plane = c._collect_map_inputs()["plane"]
        assert plane.axis_names == ("x", "z") and plane.origin[1] == pytest.approx(1.0 / qt.BOHR_RADIUS_ANGSTROM)
        c.grid_points_input.value = "5"
        with pytest.raises(ValueError, match="20 to 1000"):
            c._collect_map_inputs()


class TestWorkflow:
    def test_run_writes_the_report_and_fills_the_table(self, ran, tmp_path):
        assert ran.report_path == str(tmp_path / "water_631gd_6d_QTAIM.txt")
        assert ran.report_output.value.startswith("# gQTEA QTAIM Topology Analysis")
        assert [r[2] for r in ran.cp_table.data] == ["O1", "H2", "H3", "O1-H2", "O1-H3"]
        assert "(3,-3) 3   (3,-1) 2" in ran.summary_label.text and "(satisfied)" in ran.summary_label.text
        assert ran.progress_bar.value == 100 and ran.run_button.enabled

    def test_input_error_shows_a_dialog(self):
        c = make_controller()
        asyncio.run(c.workflow(None))
        assert c.results is None and len(c.viewer.main_window.dialogs) == 1

    def test_display_loads_the_molecule_and_the_scene(self, ran):
        asyncio.run(ran.display_results(None))
        v = ran.viewer
        assert len(v.loaded) == 1 and len(v.loaded[0]) == 3 and v.opened == 1
        assert len(v.qtaim_layer.get().cps) == 5
        asyncio.run(ran.display_results(None))   # the same molecule: no second load, no question
        assert len(v.loaded) == 1 and v.main_window.dialogs == []

    def test_a_different_structure_is_asked_about(self, ran):
        ran.viewer.frames = [[("C", (0.0, 0.0, 0.0))]]
        ran.viewer.loaded = [[("C", (0.0, 0.0, 0.0))]]
        ran.viewer.main_window.answer = False
        asyncio.run(ran.display_results(None))
        assert len(ran.viewer.main_window.dialogs) == 1 and len(ran.viewer.loaded) == 1

    def test_display_switches_follow_only_a_shown_scene(self, ran):
        ran.type_switches[qt.NCP].value = False
        ran.on_display_change(None)
        assert ran.viewer.qtaim_layer.show_types[qt.NCP]          # nothing shown yet: no change
        asyncio.run(ran.display_results(None))
        assert not ran.viewer.qtaim_layer.show_types[qt.NCP]
        ran.type_switches[qt.NCP].value = True
        ran._syncing = True
        ran.on_display_change(None)
        assert not ran.viewer.qtaim_layer.show_types[qt.NCP]      # muted while syncing
        ran._syncing = False
        ran.on_display_change(None)
        assert ran.viewer.qtaim_layer.show_types[qt.NCP]

    def test_selecting_a_row_highlights_the_cp(self, ran):
        ran.on_cp_select(W(selection=W(number="4")))
        assert ran.viewer.qtaim_layer.highlight == 4
        assert ran.status_label.text.startswith("CP 4 (3,-1)")
        ran.on_cp_select(W(selection=None))
        assert ran.viewer.qtaim_layer.highlight is None

    def test_plot_map_writes_the_manifest_and_opens_the_viewer(self, ran, tmp_path, monkeypatch):
        import displayPlots
        launched = []
        monkeypatch.setattr(displayPlots, "launch_plot_viewer", lambda p: launched.append(p) or True)
        asyncio.run(ran.plot_map(None))
        path = tmp_path / "water_631gd_6d_QTAIM_map.json"
        assert launched == [str(path)]
        entry = json.loads(path.read_text(encoding="utf-8"))[0]
        assert entry["type"] == "plane_map" and len(entry["grid"]) == 30
        assert list(tmp_path.glob("*.png")) == []

    def test_plot_map_reports_when_the_viewer_cannot_start(self, ran, monkeypatch):
        import displayPlots
        monkeypatch.setattr(displayPlots, "launch_plot_viewer", lambda p: False)
        monkeypatch.setattr(qtaimTab.toga, "InfoDialog", lambda title, message: W(title=title, message=message))
        asyncio.run(ran.plot_map(None))
        assert "Plotter (JSON plot type)" in ran.viewer.main_window.dialogs[-1].message

    def test_clear_resets_everything(self, ran):
        asyncio.run(ran.display_results(None))
        ran.clear_output(None)
        assert ran.results is None and ran.cp_table.data == [] and ran.viewer.qtaim_layer.get() is None
        assert ran.summary_label.text == "No results." and ran.report_output.value.startswith("QTAIM TOPOLOGY")


class TestViewerIntegration:
    def test_ui_draws_and_frames_the_qtaim_layer(self, water_results, no_gl):
        from molecularViewer import MolecularViewer, MolecularViewerUI
        ui = MolecularViewerUI.__new__(MolecularViewerUI)
        MolecularViewer.__init__(ui)
        ui.frames = []
        ui.molecule_data = [("H", (0.0, 0.0, 0.0))]
        plain = ui.get_scene_radius()
        ui.qtaim_layer = ql.QTAIMSceneLayer()
        drawn = []
        ui.qtaim_layer.draw_scene = lambda viewer: drawn.append(viewer)
        ui.render_scene_extras()
        assert drawn == [ui]
        ui.qtaim_layer = ql.QTAIMSceneLayer()
        ui.qtaim_layer.set_scene(ql.QTAIMScene((ql.SceneCP(1, qt.NCP, (20.0, 0.0, 0.0)),), (), ()))
        assert ui.get_scene_radius() > max(plain, 20.0)


# --------------------------------------------------------------------------- #
# Canvas: background, molecule style, labels (shared with the Display tab)      #
# --------------------------------------------------------------------------- #
class Recorder:
    """Dialog stand-in that records title and message (Toga dialogs expose neither)."""
    def __init__(self, title, message=None, **kwargs):
        self.title, self.message = title, message


class TestCanvas:
    def test_presets_are_the_design_canvas_backgrounds(self):
        from molecularDesign import CANVAS_BACKGROUNDS
        assert qtaimTab.background_presets() == CANVAS_BACKGROUNDS

    def test_preset_background_is_applied(self):
        c = make_controller()
        c.background_selection.value = "White"
        c.on_background_change(None)
        assert c.viewer.background_color == pytest.approx((1.0, 1.0, 1.0))

    def test_custom_background(self):
        c = make_controller()
        c.background_input.value = "#e8f0fb"
        c.on_custom_background(None)
        assert c.viewer.background_color == pytest.approx((0xe8 / 255, 0xf0 / 255, 0xfb / 255))

    def test_bad_custom_background_shows_an_error_and_keeps_the_colour(self, monkeypatch):
        monkeypatch.setattr(qtaimTab.toga, "ErrorDialog", Recorder)
        c = make_controller()

        async def go():
            c.background_input.value = "teal"
            c.on_custom_background(None)
            await asyncio.gather(*c._tasks)
        c.viewer.loop = None
        asyncio.run(go())
        assert c.viewer.background_color == (0.0, 0.0, 0.0)
        assert "#rrggbb" in c.viewer.main_window.dialogs[-1].message

    def test_style_and_label_switches_drive_the_viewer(self):
        c = make_controller()
        c.style_selection.value = "Design style"
        c.on_style_change(None)
        c.symbols_switch.value, c.numbers_switch.value = True, False
        c.on_label_change(None)
        v = c.viewer
        assert v.atom_display_style == "Design style" and v.show_atom_symbols and not v.show_atom_numbers

    def test_handlers_are_muted_while_syncing(self):
        c = make_controller()
        c._syncing = True
        c.style_selection.value = "vdW style"
        c.on_style_change(None)
        c.background_selection.value = "White"
        c.on_background_change(None)
        assert c.viewer.atom_display_style == "Line style" and c.viewer.background_color == (0.0, 0.0, 0.0)

    def test_sync_canvas_widgets_reflects_the_viewer(self):
        c = make_controller()
        c.viewer.atom_display_style, c.viewer.show_atom_numbers = "CPK style", True
        c.sync_canvas_widgets()
        assert c.style_selection.value == "CPK style" and c.numbers_switch.value is True
        assert c._syncing is False

    def test_sync_before_the_tab_is_built_is_harmless(self):
        QTAIMTabController(viewer=StubViewer()).sync_canvas_widgets()

    def test_initial_background_label(self):
        assert qtaimTab.background_label((0.0, 0.0, 0.0)) == "Black"
        assert qtaimTab.background_label((0.1, 0.2, 0.3)) is None

    def test_cp_labels_stay_readable_on_any_background(self, monkeypatch):
        assert ql.cp_label_colour((0.0, 0.0, 0.0)) == ql.LABEL_COLOUR
        assert ql.cp_label_colour((1.0, 1.0, 1.0)) == ql.LABEL_COLOUR_ON_LIGHT
        layer = ql.QTAIMSceneLayer()
        layer.set_scene(ql.QTAIMScene((ql.SceneCP(1, qt.NCP, (0.0, 0.0, 0.0)),), (), ()))
        viewer = RecordingViewer()
        viewer.background_color = (1.0, 1.0, 1.0)
        import qtaimLayer
        for name in ("glPushMatrix", "glPopMatrix", "glTranslatef", "glMaterialfv"):
            monkeypatch.setattr(qtaimLayer, name, lambda *a: None)
        layer.draw_scene(viewer)
        assert viewer.texts == [("1", ql.LABEL_COLOUR_ON_LIGHT)]


# --------------------------------------------------------------------------- #
# Image export: the 3D view and the 2D map                                       #
# --------------------------------------------------------------------------- #
class TestExport:
    def test_3d_export_captures_then_saves(self, ran, tmp_path, monkeypatch):
        saved = []
        import molecularViewer
        monkeypatch.setattr(molecularViewer, "save_canvas_image", lambda image, path: saved.append(path) or path)
        asked = []

        async def ask(title, suggested, types):
            asked.append((suggested, tuple(types)))
            return str(tmp_path / "view.png")
        ran._ask_save_path = ask
        asyncio.run(ran.export_3d_image(None))
        assert asked == [("water_631gd_6d_QTAIM_3D.png", ("png", "tif", "tiff", "bmp"))]
        assert saved == [str(tmp_path / "view.png")] and "6 x 4 px" in ran.status_label.text

    def test_3d_export_needs_an_open_window(self, ran, monkeypatch):
        monkeypatch.setattr(qtaimTab.toga, "InfoDialog", Recorder)
        ran.viewer.window_open = False
        asyncio.run(ran.export_3d_image(None))
        assert "Display" in ran.viewer.main_window.dialogs[-1].message

    def test_2d_export_needs_a_map(self, ran, monkeypatch):
        monkeypatch.setattr(qtaimTab.toga, "InfoDialog", Recorder)
        asyncio.run(ran.export_2d_image(None))
        assert "Plot 2D map" in ran.viewer.main_window.dialogs[-1].message

    @pytest.mark.parametrize("suffix", [".png", ".svg", ".pdf", ".tif"])
    def test_2d_export_writes_the_map(self, ran, tmp_path, monkeypatch, suffix):
        import displayPlots
        monkeypatch.setattr(displayPlots, "launch_plot_viewer", lambda p: True)
        asyncio.run(ran.plot_map(None))

        async def ask(title, suggested, types):
            assert suggested == "water_631gd_6d_QTAIM_map.png" and "svg" in types
            return str(tmp_path / ("map" + suffix))
        ran._ask_save_path = ask
        asyncio.run(ran.export_2d_image(None))
        out = tmp_path / ("map" + suffix)
        assert out.is_file() and out.stat().st_size > 1000

    def test_save_plane_map_image_rules(self, tmp_path):
        import displayPlots
        from test_displayPlots import plane_map_entry
        assert displayPlots.save_plane_map_image(plane_map_entry(), str(tmp_path / "noext")).endswith("noext.png")
        with pytest.raises(ValueError, match="png, svg, pdf, tif"):
            displayPlots.save_plane_map_image(plane_map_entry(), str(tmp_path / "map.jpg"))

    def test_clear_forgets_the_map(self, ran, monkeypatch):
        import displayPlots
        monkeypatch.setattr(displayPlots, "launch_plot_viewer", lambda p: True)
        asyncio.run(ran.plot_map(None))
        assert ran.last_map_entry is not None
        ran.clear_output(None)
        assert ran.last_map_entry is None


# --------------------------------------------------------------------------- #
# The results table keeps readable columns (WinForms sizes them only once)      #
# --------------------------------------------------------------------------- #
class TestTableColumns:
    def _table(self, width):
        columns = [W(Width=5) for _ in qtaimTab.TABLE_HEADINGS]
        return W(_impl=W(native=W(ClientSize=W(Width=width), Columns=columns))), columns

    def test_columns_follow_the_fractions(self):
        table, columns = self._table(800)
        qtaimTab.fit_table_columns(table)
        assert [c.Width for c in columns] == [int(800 * f) for f in qtaimTab.TABLE_COLUMN_FRACTIONS]
        assert sum(c.Width for c in columns) <= 800

    def test_a_tiny_table_keeps_a_minimum_width(self):
        table, columns = self._table(60)
        qtaimTab.fit_table_columns(table)
        assert all(c.Width >= qtaimTab.MIN_COLUMN_WIDTH for c in columns)

    def test_no_native_widget_is_harmless(self):
        qtaimTab.fit_table_columns(W(data=[]))
        table, columns = self._table(0)
        qtaimTab.fit_table_columns(table)
        assert [c.Width for c in columns] == [5] * len(columns)

    def test_fill_results_fits_the_columns(self, ran):
        table, columns = self._table(600)
        table.data = []
        ran.cp_table = table
        ran._fill_results()
        assert columns[0].Width == int(600 * qtaimTab.TABLE_COLUMN_FRACTIONS[0]) and len(table.data) == 5


# --------------------------------------------------------------------------- #
# Design style (flat atom discs): CPs stay visible, labels avoid atom numbers   #
# --------------------------------------------------------------------------- #
class FlatViewer(RecordingViewer):
    background_color = (1.0, 1.0, 1.0)

    def atoms_are_flat(self):
        return True

    def current_screen_axes(self):
        return np.eye(3)


class TestFlatAtoms:
    def _scene(self):
        return ql.QTAIMScene((ql.SceneCP(1, qt.NCP, (0.0, 0.0, 0.0)), ql.SceneCP(2, qt.BCP, (1.0, 0.0, 0.0))),
                             (), ())

    def test_bond_ring_cage_points_are_drawn_over_the_discs(self, no_gl):
        layer = ql.QTAIMSceneLayer()
        layer.set_scene(self._scene())
        layer.set_display_options(labels=False)
        layer.draw_scene(FlatViewer())
        assert no_gl.count("glDepthRange") == 2           # pulled forward and back: the (3,-1) only
        assert "glDisable" not in no_gl                    # the spheres keep the depth test

    def test_labels_sit_at_the_screen_lower_right(self, no_gl):
        layer = ql.QTAIMSceneLayer()
        layer.set_scene(self._scene())
        positions = []
        viewer = FlatViewer()
        viewer._draw_text_3d = lambda pos, text, color=(1, 1, 0): positions.append(np.array(pos))
        layer.draw_scene(viewer)
        assert len(positions) == 1                        # the (3,-3) repeats its atom number: not labelled
        offset = positions[0] - np.array([1.0, 0.0, 0.0])
        assert offset[0] > 0 and offset[1] < 0 and offset[2] == pytest.approx(0.0)

    def test_the_viewer_reports_flat_atoms_and_screen_axes(self, monkeypatch):
        import molecularViewer as mv
        v = mv.MolecularViewer()
        assert not v.atoms_are_flat()
        v.atom_display_style = mv.DESIGN_STYLE
        assert v.atoms_are_flat()
        monkeypatch.setattr(mv, "glGetDoublev", lambda what: np.eye(4))
        right, up, toward = v.current_screen_axes()
        assert np.allclose(right, [1, 0, 0]) and np.allclose(up, [0, 1, 0])
