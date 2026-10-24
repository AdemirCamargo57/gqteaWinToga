"""NTO tab of the 3D Molecular Viewer: the controller (headless) and its integration."""
import asyncio
import inspect
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import ntoAnalysis as na
import ntoTab
from cubeRender import CubeSceneLayer
from help import HelpGqteaWin
from nto_fixtures import fixture
from ntoTab import NTOTabController

SRC = Path(__file__).resolve().parent.parent


def W(**kwargs):
    return SimpleNamespace(**kwargs)


@pytest.fixture(autouse=True)
def plain_dialogs(monkeypatch):
    """Dialogs as plain records, so tests can read what the user would be told."""
    for name in ("ErrorDialog", "InfoDialog", "QuestionDialog"):
        monkeypatch.setattr(ntoTab.toga, name, lambda title, message, _k=name: W(kind=_k, title=title,
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
        self.nto_layer = CubeSceneLayer()
        self.cube_controller = StubCubeController()
        self.frames, self.loaded, self.refreshed, self.opened = [], [], 0, 0
        self.main_window = StubWindow()
        self.window_open, self.frame = True, np.zeros((4, 6, 3), dtype=np.uint8)
        self.claims = []

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


def make_controller(name="uracil", output_dir="", excitation=True):
    fchk, log = fixture(name)
    c = NTOTabController(viewer=StubViewer())
    c.fchk_input = W(value=str(fchk))
    c.excitation_input = W(value=str(log) if excitation else "")
    c.output_dir_input = W(value=str(output_dir))
    for button in ("load_button", "run_button", "plot_spectrum_button", "clear_button", "display_button",
                   "export_cubes_button", "export_fchk_button", "send_button", "export_3d_button"):
        setattr(c, button, W(enabled=True))
    c.state_selection = W(items=[], value=None)
    c.hwhm_input = W(value="")
    c.status_label, c.summary_label = W(text="Ready."), W(text="No results.")
    c.progress_bar, c.pair_table = W(value=0), W(data=[], selection=None)
    c.pair_selection = W(items=[], value=None)
    c.show_selection = W(value=ntoTab.SHOW_BOTH)
    c.isovalue_input, c.spacing_input = W(value=""), W(value="0.4")
    c.colour_inputs = {key: W(value=colour) for key, colour in ntoTab.DEFAULT_COLOURS.items()}
    c.opacity_slider, c.all_pairs_switch = W(value=1.0), W(value=False)
    c.report_output = W(value="")
    return c


async def settle(c):
    while c._tasks:
        await asyncio.gather(*list(c._tasks))


def choose_state(c, number):
    c.state_selection.value = next(t for t, n in c._state_items.items() if n == number)


@pytest.fixture
def ran(tmp_path):
    c = make_controller(output_dir=tmp_path)
    asyncio.run(c.load_files(None))
    choose_state(c, 3)
    asyncio.run(c.workflow(None))
    return c


@pytest.fixture
def shown(ran):
    asyncio.run(ran.display_nto(None))
    return ran


class TestInputs:
    def test_every_run_parameter_is_supplied_by_the_ui(self):
        c = make_controller()
        asyncio.run(c.load_files(None))
        wanted = set(inspect.signature(na.NTOAnalyser.run).parameters) - {"self", "progress"}
        assert set(c._collect_inputs()) == wanted

    def test_every_pair_volumes_parameter_is_supplied_by_the_ui(self, ran):
        wanted = set(inspect.signature(na.pair_volumes).parameters) - {"results", "progress", "use_numba"}
        assert set(ran._collect_display_inputs()) == wanted

    def test_display_inputs(self, ran):
        params = ran._collect_display_inputs()
        assert params == {"spin": "alpha", "pair": 1, "spacing": 0.4, "isovalue": na.ORBITAL_ISOVALUE}
        ran.spacing_input.value = "-1"
        with pytest.raises(ValueError, match="grid spacing"):
            ran._collect_display_inputs()

    def test_the_fchk_fills_in_its_gaussian_output(self):
        c = make_controller(excitation=False)
        c.set_fchk_path(str(fixture("h2co_cis")[0]))
        assert c.excitation_input.value == str(fixture("h2co_cis")[1])
        c.set_fchk_path(str(fixture("uracil")[0]))  # a typed excitation file is kept
        assert c.excitation_input.value == str(fixture("h2co_cis")[1])

    def test_inputs_need_load_and_files(self, tmp_path):
        c = make_controller()
        with pytest.raises(ValueError, match="Press Load"):
            c._collect_inputs()
        c.excitation_input.value = ""
        with pytest.raises(ValueError, match="CIS/TD output"):
            c._collect_inputs()
        c.excitation_input.value = str(tmp_path / "none.log")
        with pytest.raises(ValueError, match="File not found"):
            c._collect_inputs()


class TestLoadAndRun:
    def test_load_lists_the_states(self):
        c = make_controller()
        asyncio.run(c.load_files(None))
        assert c.state_selection.items == ['1: Singlet-A" 4.7968 eV 258.47 nm f=0.0001',
                                           "2: Singlet-A' 5.4219 eV 228.67 nm f=0.1317",
                                           '3: Singlet-A" 6.0180 eV 206.02 nm f=0.0000']
        assert c.state_selection.value == c.state_selection.items[0]
        assert "3 excited states read (restricted, 128 MOs)" in c.status_label.text

    def test_mismatched_files_are_refused(self):
        c = make_controller()
        c.excitation_input.value = str(fixture("h2co_cis")[1])
        asyncio.run(c.load_files(None))
        assert c.loaded is None and c.state_selection.items == []
        assert "different atoms" in c.viewer.main_window.dialogs[0].message

    def test_run_writes_the_report_and_fills_the_table(self, ran, tmp_path):
        assert ran.report_path == str(tmp_path / "uracil_state3_NTO.txt")
        assert ran.report_output.value.startswith("# gQTEA Natural Transition Orbitals")
        assert ran.pair_table.data[0] == ("1", "alpha", "0.865529", "86.55", "86.55")
        assert ran.pair_selection.items[:2] == ["1 (86.55 %)", "2 (13.40 %)"]
        assert ran.pair_selection.value == "1 (86.55 %)"
        assert "86.55 %" in ran.summary_label.text and ran.progress_bar.value == 100

    def test_run_loads_when_needed(self, tmp_path):
        c = make_controller(output_dir=tmp_path)
        asyncio.run(c.workflow(None))  # no Load pressed: state 1 is run
        assert c.results.state.number == 1

    def test_open_shell_pairs_name_their_spin(self, tmp_path):
        c = make_controller("ch3", output_dir=tmp_path)
        asyncio.run(c.workflow(None))
        assert c.pair_selection.items[0] == "1 beta (98.08 %)"
        assert c._collect_display_inputs()["spin"] == "beta"

    def test_input_error_shows_a_dialog(self):
        c = make_controller()
        c.fchk_input.value = ""
        asyncio.run(c.workflow(None))
        assert c.results is None and len(c.viewer.main_window.dialogs) == 1


class TestDisplay:
    def test_display_loads_the_molecule_and_four_surfaces(self, shown):
        v = shown.viewer
        assert len(v.loaded) == 1 and len(v.loaded[0]) == 12 and v.opened == 1
        assert len(v.nto_layer.items()) == 4
        import volumeData as vd
        colours = [i.rgba[0, :3] for i in v.nto_layer.items()]
        assert np.allclose(colours, [vd.parse_colour(n) for n in ("Blue", "Red", "Green", "Yellow")])
        asyncio.run(shown.display_nto(None))   # same molecule: no second load, no question
        assert len(v.loaded) == 1 and v.main_window.dialogs == []

    def test_show_hole_only(self, shown):
        async def go():
            shown.show_selection.value = ntoTab.SHOW_HOLE
            shown.on_show_change(shown.show_selection)
            await settle(shown)
        asyncio.run(go())
        assert len(shown.viewer.nto_layer.items()) == 2

    def test_colour_and_opacity_changes_reuse_the_meshes(self, shown):
        entries = shown._orbitals[1]
        before = {k: dict(e.meshes) for k, e in entries.items()}

        async def go():
            shown.colour_inputs[("hole", 1)].value = "#00ff00"
            shown.on_colour_confirm(shown.colour_inputs[("hole", 1)])
            shown.opacity_slider.value = 0.5
            shown.on_opacity_change(shown.opacity_slider)
            await settle(shown)
        asyncio.run(go())
        items = shown.viewer.nto_layer.items()
        assert np.allclose(items[0].rgba[0], [0.0, 1.0, 0.0, 0.5])
        for kind, entry in entries.items():
            assert all(entry.meshes[key] is mesh for key, mesh in before[kind].items())

    def test_a_bad_colour_is_refused_and_restored(self, shown):
        async def go():
            shown.colour_inputs[("particle", -1)].value = "mauve"
            shown.on_colour_confirm(shown.colour_inputs[("particle", -1)])
            await settle(shown)
        asyncio.run(go())
        assert shown.colour_inputs[("particle", -1)].value == "Yellow"
        assert "Unknown colour" in shown.viewer.main_window.dialogs[-1].message

    def test_a_new_isovalue_recomputes_the_orbitals(self, shown):
        old = shown._orbitals

        async def go():
            shown.isovalue_input.value = "0.05"
            shown.on_grid_confirm(shown.isovalue_input)
            await settle(shown)
        asyncio.run(go())
        assert shown._orbitals is not old and shown._orbitals[1]["hole"].isovalue == 0.05

    def test_another_pair(self, shown):
        async def go():
            shown.pair_table.selection = W(pair="2", spin="alpha")
            shown.on_pair_select(shown.pair_table)
            await settle(shown)
        asyncio.run(go())
        assert shown.pair_selection.value == "2 (13.40 %)" and shown._orbitals[0][2] == 2

    def test_handlers_are_muted_while_syncing(self, shown):
        shown._syncing = True
        shown.show_selection.value = ntoTab.SHOW_PARTICLE
        shown.on_show_change(shown.show_selection)
        shown.on_pair_change(shown.pair_selection)
        assert shown.show == ntoTab.SHOW_BOTH and not shown._tasks

    def test_display_before_run(self):
        c = make_controller()
        asyncio.run(c.display_nto(None))
        assert "Run the analysis first" in c.viewer.main_window.dialogs[0].message


class TestExport:
    def test_export_the_chosen_pair(self, ran, tmp_path):
        asyncio.run(ran.export_cubes(None))
        assert sorted(p.name for p in tmp_path.glob("*.cube")) == ["uracil_state3_NTO1_hole.cube",
                                                                   "uracil_state3_NTO1_particle.cube"]

    def test_export_all_listed_pairs(self, ran, tmp_path):
        ran.all_pairs_switch.value = True
        ran.spacing_input.value = "0.6"
        asyncio.run(ran.export_cubes(None))
        assert len(list(tmp_path.glob("*.cube"))) == 2 * len(ran.results.listed_pairs())

    def test_export_fchk(self, ran, tmp_path):
        asyncio.run(ran.export_fchk(None))
        assert (tmp_path / "uracil_state3_NTO.fchk").is_file()

    def test_send_to_cube_tab(self, shown):
        asyncio.run(shown.send_to_cube_tab(None))
        cube = shown.viewer.cube_controller
        assert [e.name for e in cube.added] == ["uracil_state3_NTO1_hole", "uracil_state3_NTO1_particle"]
        assert not any(e.visible for e in cube.added) and cube.synced == ["uracil_state3_NTO1_hole"]

    def test_send_needs_a_displayed_pair(self, ran):
        asyncio.run(ran.send_to_cube_tab(None))
        assert ran.viewer.cube_controller.added == [] and ran.viewer.main_window.dialogs

    def test_export_3d_image_captures_first(self, shown, tmp_path, monkeypatch):
        asked = []

        async def ask(title, suggested, types):
            asked.append(suggested)
            return str(tmp_path / suggested)
        monkeypatch.setattr(shown, "_ask_save_path", ask)
        asyncio.run(shown.export_3d_image(None))
        assert asked == ["uracil_state3_NTO1.png"] and (tmp_path / "uracil_state3_NTO1.png").is_file()


class TestPlotSpectrum:
    @pytest.fixture
    def launched(self, monkeypatch):
        import displayPlots
        calls = []
        monkeypatch.setattr(displayPlots, "launch_plot_viewer", lambda path: calls.append(path) or True)
        return calls

    def read(self, path):
        import json
        return json.loads(Path(path).read_text(encoding="utf-8"))[0]

    def test_every_state_is_plotted_and_saved(self, launched, tmp_path):
        c = make_controller(output_dir=tmp_path)
        asyncio.run(c.load_files(None))
        choose_state(c, 3)                      # the selection does not limit the plot
        asyncio.run(c.plot_spectrum(None))
        path = tmp_path / "uracil_UVVis.json"
        assert launched == [str(path)]
        entry = self.read(path)
        assert entry["type"] == "uv_vis" and [s["number"] for s in entry["states"]] == [1, 2, 3]
        assert [s["oscillator"] for s in entry["states"]] == [0.0001, 0.1317, 0.0]
        assert entry["broadening"]["hwhm_ev"] == 0.333
        assert entry["source"] == {"fchk": "uracil.fch", "excitation_file": "uracil.out"}
        assert str(path) in c.status_label.text

    def test_loads_the_files_when_needed(self, launched, tmp_path):
        c = make_controller("h2co_td5050", output_dir=tmp_path)
        asyncio.run(c.plot_spectrum(None))     # no Load pressed
        assert len(self.read(tmp_path / "h2co_td5050_631gd_UVVis.json")["states"]) == 8
        assert c.state_selection.items          # and the states are listed as Load would

    def test_half_width_field(self, launched, tmp_path):
        c = make_controller(output_dir=tmp_path)
        c.hwhm_input.value = "0,2"
        asyncio.run(c.plot_spectrum(None))
        assert self.read(tmp_path / "uracil_UVVis.json")["broadening"]["hwhm_ev"] == 0.2
        c.hwhm_input.value = "-1"
        asyncio.run(c.plot_spectrum(None))
        assert "half-width" in c.viewer.main_window.dialogs[-1].message and len(launched) == 1

    def test_viewer_that_cannot_start(self, monkeypatch, tmp_path):
        import displayPlots
        monkeypatch.setattr(displayPlots, "launch_plot_viewer", lambda path: False)
        c = make_controller(output_dir=tmp_path)
        asyncio.run(c.plot_spectrum(None))
        message = c.viewer.main_window.dialogs[-1].message
        assert "uracil_UVVis.json" in message and "Plotter" in message

    def test_plain_text_file_is_refused(self, launched, tmp_path):
        text = tmp_path / "exc.txt"
        text.write_text("Excited State   1   1   5.0\n   29 -> 30   0.70\n")
        c = make_controller(output_dir=tmp_path)
        c.excitation_input.value = str(text)
        asyncio.run(c.plot_spectrum(None))
        assert "oscillator strengths" in c.viewer.main_window.dialogs[-1].message
        assert launched == [] and not list(tmp_path.glob("*.json"))

    def test_the_nto_results_and_view_are_untouched(self, launched, shown):
        results, orbitals, items = shown.results, shown._orbitals, shown.viewer.nto_layer.items()
        asyncio.run(shown.plot_spectrum(None))
        assert shown.results is results and shown._orbitals is orbitals
        assert shown.viewer.nto_layer.items() == items and len(launched) == 1


def test_clear(shown):
    shown.clear_output(None)
    assert shown.results is None and shown.pair_table.data == [] and shown.pair_selection.items == []
    assert shown.viewer.nto_layer.items() == () and shown.report_output.value == HelpGqteaWin.help_nto


def test_tab_does_not_import_the_viewer_at_module_level():
    code = "import sys, ntoTab; print('molecularViewer' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"


class TestCanvasOwnership:
    """Displaying replaces whatever another tab left on the canvas (ESP surface, hole/electron, ...)."""

    def test_displaying_claims_the_canvas(self, shown):
        assert shown.viewer.claims and set(shown.viewer.claims) == {"nto_layer"}
        assert shown.viewer.nto_layer.items()
