"""Excitation tab of the 3D Molecular Viewer: the controller (headless) and its integration."""
import asyncio
import inspect
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import excitationAnalysis as xa
import excitationClassification as xc
import excitationTab
import excitation_fixtures as xf
from cubeRender import CubeSceneLayer
from excitationTab import ExcitationTabController
from help import HelpGqteaWin

SRC = Path(__file__).resolve().parent.parent


def W(**kwargs):
    return SimpleNamespace(**kwargs)


@pytest.fixture(autouse=True)
def plain_dialogs(monkeypatch):
    """Dialogs as plain records, so tests can read what the user would be told."""
    for name in ("ErrorDialog", "InfoDialog", "QuestionDialog"):
        monkeypatch.setattr(excitationTab.toga, name, lambda title, message, _k=name: W(kind=_k, title=title,
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
        self.excitation_layer = CubeSceneLayer()
        self.cube_controller = StubCubeController()
        self.frames, self.loaded, self.refreshed, self.opened = [], [], 0, 0
        self.main_window = StubWindow()
        self.claims = []

    def show_analysis_layer(self, name):
        self.claims.append(name)

    def viewer_window_is_open(self):
        return False

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


def make_controller(name="h2co_cis", output_dir=""):
    fchk, log = xf.fixture_paths(name)
    c = ExcitationTabController(viewer=StubViewer())
    c.fchk_input, c.excitation_input = W(value=str(fchk)), W(value=str(log))
    c.output_dir_input = W(value=str(output_dir))
    for button in ("load_button", "run_button", "clear_button", "display_button", "export_cubes_button",
                   "send_button", "heatmap_button", "export_3d_button"):
        setattr(c, button, W(enabled=True))
    c.state_selection = W(items=[], value=None)
    c.partition_selection = W(value="Mulliken-like")
    c.fragments_input, c.hydrogens_switch = W(value=""), W(value=False)
    c.fast_coulomb_switch, c.precise_coulomb_switch = W(value=False, enabled=True), W(value=False, enabled=True)
    c.status_label, c.progress_bar = W(text="Ready."), W(value=0)
    c.classification_label, c.summary_label = W(text="No results."), W(text="")
    c.descriptor_table, c.part_table = W(data=[]), W(data=[])
    c.show_selection = W(value=excitationTab.DEFAULT_SHOW)
    c.isovalue_input, c.spacing_input = W(value=""), W(value="0.35")
    c.colour_inputs = {key: W(value=colour) for key, colour in excitationTab.DEFAULT_COLOURS.items()}
    c.opacity_slider, c.centroids_switch = W(value=1.0), W(value=True)
    c.heatmap_selection = W(value=xa.HEATMAP_KINDS["atoms"])
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
    choose_state(c, 1)
    asyncio.run(c.workflow(None))
    return c


@pytest.fixture
def shown(ran):
    asyncio.run(ran.display_excitation(None))
    return ran


class TestInputs:
    def test_every_run_parameter_is_supplied_by_the_ui(self):
        c = make_controller()
        asyncio.run(c.load_files(None))
        wanted = set(inspect.signature(xa.ExcitationAnalyser.run).parameters) - {"self", "progress"}
        assert set(c._collect_inputs()) == wanted

    def test_every_field_volumes_grid_parameter_is_supplied_by_the_ui(self, ran):
        wanted = set(inspect.signature(xa.field_volumes).parameters) - {"results", "progress", "use_numba"}
        assert set(ran._collect_display_inputs()) == wanted

    def test_options_reach_the_run(self):
        c = make_controller()
        asyncio.run(c.load_files(None))
        c.partition_selection.value, c.fragments_input.value, c.precise_coulomb_switch.value = "Hirshfeld", "1-2", True
        params = c._collect_inputs()
        assert (params["partition"], params["fragments_text"], params["coulomb_mode"]) == ("hirshfeld", "1-2",
                                                                                           xa.COULOMB_PRECISE)


class TestCoulombCheckboxes:
    def tick(self, c, widget, value):
        widget.value = value
        c.on_coulomb_mode_change(widget)

    def test_neither_ticked_means_no_exact_coulomb(self):
        c = make_controller()
        assert c.coulomb_mode() == xa.COULOMB_OFF
        assert c.fast_coulomb_switch.enabled and c.precise_coulomb_switch.enabled

    def test_ticking_fast_disables_precise_and_unticking_enables_it(self):
        c = make_controller()
        self.tick(c, c.fast_coulomb_switch, True)
        assert c.coulomb_mode() == xa.COULOMB_FAST and c.precise_coulomb_switch.enabled is False
        self.tick(c, c.fast_coulomb_switch, False)
        assert c.coulomb_mode() == xa.COULOMB_OFF and c.precise_coulomb_switch.enabled is True

    def test_ticking_precise_disables_fast_and_unticking_enables_it(self):
        c = make_controller()
        self.tick(c, c.precise_coulomb_switch, True)
        assert c.coulomb_mode() == xa.COULOMB_PRECISE and c.fast_coulomb_switch.enabled is False
        self.tick(c, c.precise_coulomb_switch, False)
        assert c.coulomb_mode() == xa.COULOMB_OFF and c.fast_coulomb_switch.enabled is True

    def test_the_chosen_mode_reaches_the_run(self, tmp_path):
        c = make_controller(output_dir=tmp_path)
        asyncio.run(c.load_files(None))
        self.tick(c, c.fast_coulomb_switch, True)
        asyncio.run(c.workflow(None))
        assert c.results.coulomb.mode == xa.COULOMB_FAST and c.results.coulomb.exact > 0
        assert any("exact, fast" in row[0] for row in c.descriptor_table.data)

    def test_a_bad_fragment_definition_is_reported_before_running(self):
        c = make_controller()
        asyncio.run(c.load_files(None))
        c.fragments_input.value = "1-2; 2-3"
        asyncio.run(c.workflow(None))
        assert c.results is None and "fragments 1 and 2" in c.viewer.main_window.dialogs[-1].message

    def test_blank_isovalue_means_per_field_defaults(self, ran):
        assert ran._collect_display_inputs() == {"spacing": 0.35, "isovalue": None}
        ran.isovalue_input.value = "0.004"
        assert ran._collect_display_inputs()["isovalue"] == 0.004

    def test_load_lists_the_states(self):
        c = make_controller()
        asyncio.run(c.load_files(None))
        assert len(c._state_items) == 3 and c.state_selection.value.startswith("1: ")


class TestRun:
    def test_results_tables_and_report(self, ran):
        assert ran.results is not None and ran.report_path and Path(ran.report_path).is_file()
        assert ran.classification_label.text == f"Classification: {ran.results.classification.verdict}"
        quantities = [row[0] for row in ran.descriptor_table.data]
        for needle in ("D index", "t index", "Sr index", "Delta sigma", "Delta_r", "Lambda", "Ghost-hunter",
                       "Coulomb attraction"):
            assert any(needle in q for q in quantities), needle
        assert "CLASSIFICATION" in ran.report_output.value

    def test_atom_rows_hide_hydrogens_unless_asked(self, ran):
        heavy = len(ran.part_table.data)
        ran.hydrogens_switch.value = True
        ran.on_hydrogens_change(ran.hydrogens_switch)
        assert len(ran.part_table.data) == heavy + 2   # formaldehyde: two H

    def test_fragments_fill_the_part_table(self, tmp_path):
        c = make_controller("4-nitroaniline", tmp_path)
        asyncio.run(c.load_files(None))
        choose_state(c, 2)
        c.fragments_input.value = xf.FRAGMENTS["4-nitroaniline"]
        asyncio.run(c.workflow(None))
        assert [row[0] for row in c.part_table.data][:3] == ["Fragment 1", "Fragment 2", "Fragment 3"]
        assert c.part_table.data[-1][0] == "Intrinsic CT / LE"
        assert c.results.classification.verdict == xc.VERDICT_CT

    def test_results_kept_when_the_report_cannot_be_written(self, tmp_path):
        blocker = tmp_path / "file"
        blocker.write_text("")
        c = make_controller(output_dir=blocker / "sub")
        asyncio.run(c.load_files(None))
        asyncio.run(c.workflow(None))
        assert c.results is not None
        assert c.viewer.main_window.dialogs[-1].title == "Excitation Output Error"
        assert "CLASSIFICATION" in c.report_output.value


class TestDisplay:
    def test_display_draws_surfaces_and_markers(self, shown):
        items = shown.viewer.excitation_layer.items()
        assert len(items) >= 4      # hole, electron, two centroid spheres (+ arrow)
        assert shown.viewer.opened == 1 and shown.viewer.loaded

    def test_show_change_recolours_without_recomputing(self, shown, monkeypatch):
        calls = []
        real = xa.field_volumes
        monkeypatch.setattr(xa, "field_volumes", lambda *a, **k: calls.append(1) or real(*a, **k))
        async def go():
            shown.show_selection.value = "Charge density difference"
            shown.on_show_change(shown.show_selection)
            await settle(shown)
        asyncio.run(go())
        assert calls == [] and shown.show == "Charge density difference"

    def test_spacing_change_recomputes(self, shown, monkeypatch):
        calls = []
        real = xa.field_volumes
        monkeypatch.setattr(xa, "field_volumes", lambda *a, **k: calls.append(1) or real(*a, **k))
        shown.spacing_input.value = "0.4"
        asyncio.run(shown._display())
        assert calls == [1]

    def test_centroids_switch(self, shown):
        before = len(shown.viewer.excitation_layer.items())
        async def go():
            shown.centroids_switch.value = False
            shown.on_centroids_change(shown.centroids_switch)
            await settle(shown)
        asyncio.run(go())
        assert len(shown.viewer.excitation_layer.items()) < before

    def test_signed_fields_draw_both_signs(self, shown):
        shown.show = "Transition density"
        items, notes = shown.compute_scene_items()
        surfaces = len(items) - (3 if shown.results.moments.d_index >= excitationTab.MIN_ARROW_LENGTH else 2)
        assert surfaces + sum("No " in n for n in notes) == 2

    def test_bad_colour_is_put_back(self, shown):
        widget = shown.colour_inputs[("hole", 1)]
        async def go():
            widget.value = "not a colour"
            shown.on_colour_confirm(widget)
            await settle(shown)
        asyncio.run(go())
        assert widget.value == "Blue" and shown.viewer.main_window.dialogs[-1].title == "Colour"

    def test_send_to_cube_tab(self, shown):
        asyncio.run(shown.send_to_cube_tab(None))
        added = shown.viewer.cube_controller.added
        assert [e.name for e in added] == ["h2co_cis_ccpvdz_state1_hole", "h2co_cis_ccpvdz_state1_electron"]
        assert not any(e.visible for e in added)

    def test_handlers_are_muted_while_syncing(self, shown):
        shown._syncing = True
        shown.show_selection.value = "Hole"
        shown.on_show_change(shown.show_selection)
        shown.on_grid_confirm(shown.spacing_input)
        assert shown.show == excitationTab.DEFAULT_SHOW and not shown._tasks

    def test_clear(self, shown):
        shown.clear_output(None)
        assert shown.results is None and not shown.viewer.excitation_layer.items()
        assert shown.report_output.value == HelpGqteaWin.help_excitation


class TestExports:
    def test_export_cubes_writes_the_shown_quantities(self, ran, tmp_path):
        asyncio.run(ran.export_cubes(None))
        names = sorted(p.name for p in Path(ran.output_dir()).glob("*.cube"))
        assert names == ["h2co_cis_ccpvdz_state1_electron.cube", "h2co_cis_ccpvdz_state1_hole.cube"]

    def test_heat_map_is_written_and_viewer_launched(self, ran, monkeypatch):
        import displayPlots
        launched = []
        monkeypatch.setattr(displayPlots, "launch_plot_viewer", lambda path: launched.append(path) or True)
        asyncio.run(ran.plot_heat_map(None))
        payload = json.loads(Path(launched[0]).read_text(encoding="utf-8"))
        assert payload[0]["type"] == "heat_map"

    def test_fragment_heat_map_without_fragments_is_an_error(self, ran):
        ran.heatmap_selection.value = xa.HEATMAP_KINDS["ifct"]
        asyncio.run(ran.plot_heat_map(None))
        assert ran.viewer.main_window.dialogs[-1].title == "Heat Map Error"

    def test_3d_image_needs_an_open_window(self, ran):
        asyncio.run(ran.export_3d_image(None))
        assert ran.viewer.main_window.dialogs[-1].title == "No Image"


def test_help_text_mentions_the_rules_and_iop():
    text = HelpGqteaWin.help_excitation
    assert "IOp(9/40=4)" in text and "D >= 2.0 A and t > 0" in text and "heuristics" in text


def test_tab_does_not_import_the_viewer_at_module_level():
    code = "import sys, excitationTab; print('molecularViewer' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"


class TestCanvasOwnership:
    """Displaying replaces whatever another tab left on the canvas (ESP surface, NTOs, ...)."""

    def test_displaying_claims_the_canvas(self, shown):
        assert shown.viewer.claims and set(shown.viewer.claims) == {"excitation_layer"}
        assert shown.viewer.excitation_layer.items()

    def test_an_empty_scene_claims_nothing(self, ran):
        ran._volumes = None  # no field grids: compute_scene_items has nothing to draw
        asyncio.run(ran.refresh_scene())
        assert not ran.viewer.excitation_layer.items()
        assert ran.viewer.claims == []
