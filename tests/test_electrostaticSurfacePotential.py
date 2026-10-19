"""Tests for the ESP surface analyser, the extrema file and the UI plumbing."""
import asyncio
from types import SimpleNamespace

import numpy as np
import pytest

import electrostaticSurfacePotential as esp
from esp_fixtures import fchk_path

A0 = 0.529177210544


class TestParsing:
    def test_blank_isovalue_is_the_default(self):
        assert esp.parse_isovalues("") == [0.0004] and esp.parse_isovalues(None) == [0.0004]

    def test_several_isovalues_keep_order_without_duplicates(self):
        assert esp.parse_isovalues("0.001, 0.0004;0.002 0.001") == [0.001, 0.0004, 0.002]

    @pytest.mark.parametrize("text", ["abc", "0", "-0.001", "2.5"])
    def test_bad_isovalue(self, text):
        with pytest.raises(ValueError):
            esp.parse_isovalues(text)

    def test_esp_like_value_is_explained(self):
        with pytest.raises(ValueError, match="electron-density"):
            esp.parse_isovalues("1.5")

    def test_spacing(self):
        assert esp.parse_spacing("") == 0.25 and esp.parse_spacing("0.2") == 0.2
        with pytest.raises(ValueError):
            esp.parse_spacing("1.0")

    def test_output_dir_defaults_to_the_fchk_folder(self, tmp_path):
        f = tmp_path / "x.fchk"
        assert esp.resolve_output_dir(str(f), "") == str(tmp_path.resolve())
        assert esp.resolve_output_dir(str(f), "D:/out") == "D:/out"


@pytest.fixture(scope="module")
def eteno_results():
    analyser = esp.ESPSurfaceAnalyser()
    seen = []
    results = analyser.run(str(fchk_path("eteno_b3lyp_631g")), [0.0004, 0.001],
                           progress=lambda f, m: seen.append(f))
    return analyser, results, seen


class TestEteno:
    def test_progress_is_monotone_to_one(self, eteno_results):
        _, _, seen = eteno_results
        assert seen[-1] == 1.0 and all(b >= a - 1e-12 for a, b in zip(seen, seen[1:]))

    def test_two_pi_minima_and_four_h_maxima_per_surface(self, eteno_results):
        analyser, results, _ = eteno_results
        for surface in results.surfaces:
            kinds = [(e.kind, e.scope) for e in surface.extrema]
            assert kinds.count(("minimum", "global")) == 1 and kinds.count(("minimum", "local")) == 1
            assert kinds.count(("maximum", "global")) == 1 and kinds.count(("maximum", "local")) == 3
            assert surface.mesh.component_euler_characteristics() == [2]
        rows = analyser.extrema_rows(results)
        assert {r["nearest_atom"][0] for r in rows if r["kind"] == "minimum"} == {"C"}
        assert {r["nearest_atom"][0] for r in rows if r["kind"] == "maximum"} == {"H"}

    def test_values_at_rho_0_0004(self, eteno_results):
        """Spike reference (Multiwfn-style); GaussView 6.1 shows a +/-2.072e-2 bar for this file."""
        _, results, _ = eteno_results
        s = results.surfaces[0]
        assert s.isovalue == 0.0004
        assert s.esp.min() == pytest.approx(-0.02097, abs=2e-4)
        assert s.esp.max() == pytest.approx(0.01166, abs=2e-4)
        assert s.half_range == pytest.approx(0.0207, abs=5e-4)
        # the minima sit over the pi cloud: off the molecular (yz) plane along x
        minima = [e for e in s.extrema if e.kind == "minimum"]
        assert all(abs(e.position[0] * A0) > 2.0 for e in minima)

    def test_tighter_isovalue_gives_a_smaller_surface_and_deeper_minimum(self, eteno_results):
        _, results, _ = eteno_results
        outer, inner = results.surfaces
        assert inner.mesh.area() < outer.mesh.area()
        assert inner.esp.min() < outer.esp.min()


class TestExtremaFile:
    def test_format(self, eteno_results, tmp_path):
        analyser, results, _ = eteno_results
        path = analyser.write_extrema_file(results, str(tmp_path))
        assert path.endswith("eteno_b3lyp_631g_ESP_extrema.txt")
        text = open(path, encoding="utf-8").read()
        assert "# surface_definition     electron density isosurface" in text
        assert "# mapped_property        electrostatic potential" in text
        assert "# hartree_to_kcal_mol    627.5094740631" in text
        assert text.count("# surface rho_isovalue_au=") == 2
        data = [l.split() for l in text.splitlines() if not l.startswith("#")]
        assert len(data) == 12  # 6 extrema on each of the two surfaces
        for row in data:
            assert row[2] in ("minimum", "maximum") and row[3] in ("global", "local")
            assert float(row[5]) == pytest.approx(float(row[4]) * 627.5094740631, abs=1e-3)
        globals_ = [(r[1], r[2]) for r in data if r[3] == "global"]
        assert sorted(globals_) == sorted({(r[1], r[2]) for r in data})  # one global per type per surface

    def test_unwritable_folder_is_a_value_error(self, eteno_results, tmp_path):
        analyser, results, _ = eteno_results
        blocker = tmp_path / "file"
        blocker.write_text("x")
        with pytest.raises(ValueError, match="Could not write"):
            analyser.write_extrema_file(results, str(blocker / "sub"))

    def test_summary_names_the_surface_and_the_file(self, eteno_results):
        analyser, results, _ = eteno_results
        text = analyser.summary_text(results, "C:/x/eteno_ESP_extrema.txt")
        assert "rho = 0.000400 a.u. (electron density isovalue)" in text
        assert "Extrema written to: C:/x/eteno_ESP_extrema.txt" in text


class TestEngines:
    def test_missing_cubegen_falls_back_with_a_note(self, monkeypatch):
        monkeypatch.setattr(esp, "find_cubegen", lambda explicit=None: None)
        r = esp.ESPSurfaceAnalyser().run(str(fchk_path("water_631gd_6d")), [0.001],
                                         engine=esp.ENGINE_CUBEGEN)
        assert r.engine == esp.ENGINE_BUILTIN and "not found" in r.notes[0]

    def test_failing_cubegen_falls_back_with_a_note(self, monkeypatch):
        monkeypatch.setattr(esp, "find_cubegen", lambda explicit=None: "C:/fake/cubegen.exe")

        def broken(*args, **kwargs):
            raise RuntimeError("cubegen failed: test")
        monkeypatch.setattr(esp, "cubegen_potential", broken)
        r = esp.ESPSurfaceAnalyser().run(str(fchk_path("water_631gd_6d")), [0.001],
                                         engine=esp.ENGINE_CUBEGEN)
        assert r.engine == esp.ENGINE_BUILTIN and "cubegen failed" in r.notes[0]

    @pytest.mark.skipif(esp.find_cubegen() is None, reason="Gaussian cubegen not installed")
    def test_cubegen_engine_agrees_with_the_built_in_one(self):
        f = str(fchk_path("water_631gd_6d"))
        a = esp.ESPSurfaceAnalyser().run(f, [0.001], engine=esp.ENGINE_BUILTIN).surfaces[0]
        b = esp.ESPSurfaceAnalyser().run(f, [0.001], engine=esp.ENGINE_CUBEGEN).surfaces[0]
        assert np.max(np.abs(a.esp - b.esp)) < 1e-6

    def test_unknown_engine(self):
        with pytest.raises(ValueError):
            esp.ESPSurfaceAnalyser().run(str(fchk_path("water_631gd_6d")), engine="magic")


def make_ui():
    """The ESP tab without a window: stub widgets, as in test_meanResidenceTime."""
    ui = esp.ESPTabController.__new__(esp.ESPTabController)
    ui.analyser = esp.ESPSurfaceAnalyser()
    ui.results = None
    ui.extrema_path = None
    ui.viewer = None
    ui._syncing = False
    ui._structure_settled_for = None
    for name in ("fchk_input", "output_dir_input", "isovalue_input", "spacing_input", "cubegen_input"):
        setattr(ui, name, SimpleNamespace(value=""))
    ui.engine_selection = SimpleNamespace(value=esp.ENGINE_BUILTIN)
    ui.surface_selection = SimpleNamespace(value=None, items=[])
    ui.unit_selection = SimpleNamespace(value="a.u.")
    ui.opacity_selection = SimpleNamespace(value="Opaque")
    ui.extrema_switch = SimpleNamespace(value=True)
    ui.background_selection = SimpleNamespace(value="Black")
    ui.summary_output = SimpleNamespace(value="")
    ui.status_label = SimpleNamespace(text="")
    ui.progress_bar = SimpleNamespace(value=0)
    ui.density_cube_switch = SimpleNamespace(value=True)
    for name in ("run_button", "show_button", "clear_button", "cube_button"):
        setattr(ui, name, SimpleNamespace(enabled=True))
    ui.errors = []
    ui.questions = []
    ui.answer = True

    async def record(title, message):
        ui.errors.append((title, message))

    async def ask(title, message):
        ui.questions.append(title)
        return ui.answer
    ui._show_error = record
    ui._show_info = record
    ui._ask = ask
    return ui


class TestUI:
    def test_missing_file_is_an_input_error(self):
        ui = make_ui()
        asyncio.run(ui.workflow(None))
        assert ui.errors and ui.errors[0][0] == "Input Error"

    def test_collect_inputs(self):
        ui = make_ui()
        ui.fchk_input.value = str(fchk_path("eteno_b3lyp_631g"))
        ui.isovalue_input.value = "0.0004 0.001"
        p = ui._collect_inputs()
        assert p["isovalues"] == [0.0004, 0.001] and p["spacing"] == 0.25 and p["cubegen_path"] is None

    def test_every_run_parameter_is_supplied_by_the_ui(self):
        import inspect
        ui = make_ui()
        ui.fchk_input.value = str(fchk_path("eteno_b3lyp_631g"))
        wanted = set(inspect.signature(esp.ESPSurfaceAnalyser.run).parameters) - {"self", "progress"}
        assert wanted == set(ui._collect_inputs())

    def test_workflow_writes_the_file_and_fills_the_selection(self, tmp_path):
        ui = make_ui()
        ui.fchk_input.value = str(fchk_path("water_631gd_6d"))
        ui.output_dir_input.value = str(tmp_path)
        ui.isovalue_input.value = "0.001"
        asyncio.run(ui.workflow(None))
        assert not ui.errors
        assert (tmp_path / "water_631gd_6d_ESP_extrema.txt").exists()
        assert ui.surface_selection.items == ["rho = 0.001000 a.u."]
        assert ui.selected_surface().isovalue == 0.001
        assert ui.run_button.enabled
        markers = ui.surface_markers(ui.selected_surface())
        assert markers and all(m.label.isdigit() for m in markers)

    def test_clear_resets_results(self, tmp_path):
        ui = make_ui()
        ui.fchk_input.value = str(fchk_path("water_631gd_6d"))
        ui.output_dir_input.value = str(tmp_path)
        asyncio.run(ui.workflow(None))
        ui.clear_output(None)
        assert ui.results is None and ui.surface_selection.items == []
        asyncio.run(ui.show_surface(None))
        assert ui.errors[-1][0] == "No Results"


class StubViewer:
    """Stands in for the 3D Molecular Viewer: an ESP layer, a structure, one canned frame."""

    def __init__(self, image=None, open_=True, error=None, frame=None):
        from espSurfaceLayer import ESPSurfaceLayer
        self.image = image
        self.open = open_
        self.error = error
        self.background = None
        self.captures = 0
        self.esp_layer = ESPSurfaceLayer()
        self.frame = list(frame or [])
        self.loaded = []

    @property
    def frames(self):
        """Loaded structures, as on MolecularViewer (its built-in placeholder is not one)."""
        return [list(self.frame)] if self.frame else []

    def get_current_frame_data(self):
        return list(self.frame)

    def load_designed_molecule(self, frame, bonds=None):
        self.loaded.append(frame)
        self.frame = list(frame)

    def get_system_center(self, frame_data=None):
        return tuple(np.mean([p for _, p in self.frame], axis=0)) if self.frame else (0.0, 0.0, 0.0)

    def invalidate_scene_cache(self):
        pass

    def request_projection_update(self):
        pass

    async def ensure_viewer_window(self, timeout=5.0):
        return True, None

    def viewer_window_is_open(self):
        return self.open

    def set_background_color(self, rgb):
        self.background = rgb

    async def capture_frame(self, timeout=5.0):
        self.captures += 1
        if self.error:
            raise RuntimeError(self.error)
        return self.image


class TestCanvasBackground:
    def test_black_and_white_are_offered_and_black_is_the_default(self):
        assert list(esp.BACKGROUND_COLOURS)[0] == "Black"
        assert esp.BACKGROUND_COLOURS["Black"] == (0.0, 0.0, 0.0)
        assert esp.BACKGROUND_COLOURS["White"] == (1.0, 1.0, 1.0)

    def test_change_is_applied_to_an_open_viewer_immediately(self):
        ui = make_ui()
        ui.viewer = StubViewer()
        ui.background_selection.value = "White"
        ui.on_background_change(ui.background_selection)
        assert ui.viewer.background == (1.0, 1.0, 1.0)

    def test_change_without_results_still_applies(self):
        ui = make_ui()  # no run yet: the background is a view setting, not a result
        ui.viewer = StubViewer()
        ui.background_selection.value = "Light grey"
        ui.on_background_change(None)
        assert ui.viewer.background == esp.BACKGROUND_COLOURS["Light grey"]

    def test_change_without_a_viewer_is_harmless(self):
        ui = make_ui()
        ui.background_selection.value = "White"
        ui.on_background_change(None)
        assert ui.selected_background() == (1.0, 1.0, 1.0)


class TestExportImage:
    def _ui_with_results(self, tmp_path):
        ui = make_ui()
        ui.fchk_input.value = str(fchk_path("water_631gd_6d"))
        ui.output_dir_input.value = str(tmp_path)
        ui.isovalue_input.value = "0.001"
        asyncio.run(ui.workflow(None))
        assert not ui.errors
        return ui

    def test_without_an_open_viewer_explains_what_to_do(self):
        ui = make_ui()
        asyncio.run(ui.export_image(None))
        assert ui.errors[-1][0] == "No Image"
        ui.viewer = StubViewer(open_=False)
        asyncio.run(ui.export_image(None))
        assert ui.errors[-1][0] == "No Image" and ui.viewer.captures == 0

    def test_default_name_names_the_file_and_the_isovalue(self, tmp_path):
        ui = self._ui_with_results(tmp_path)
        assert ui.default_image_name() == "water_631gd_6d_ESP_rho0.001000.png"
        assert make_ui().default_image_name() == "ESP_surface.png"

    def test_saves_exactly_the_captured_pixels(self, tmp_path):
        from PIL import Image

        ui = self._ui_with_results(tmp_path)
        image = np.random.default_rng(1).integers(0, 256, size=(6, 9, 3), dtype=np.uint8)
        ui.viewer = StubViewer(image)
        asked = []

        async def ask(suggested):
            asked.append(suggested)
            return tmp_path / "picture.png"
        ui._ask_image_path = ask
        asyncio.run(ui.export_image(None))
        assert not ui.errors and asked == ["water_631gd_6d_ESP_rho0.001000.png"]
        with Image.open(tmp_path / "picture.png") as saved:
            assert np.array_equal(np.asarray(saved.convert("RGB")), image)
        assert "picture.png" in ui.status_label.text and "9 x 6" in ui.status_label.text

    def test_the_frame_is_captured_before_the_dialog_opens(self, tmp_path):
        ui = make_ui()
        ui.viewer = StubViewer(np.zeros((2, 2, 3), np.uint8))
        order = []

        async def ask(suggested):
            order.append(("dialog", ui.viewer.captures))
            return None
        ui._ask_image_path = ask
        asyncio.run(ui.export_image(None))
        assert order == [("dialog", 1)] and not ui.errors  # cancelled: nothing written, no error

    def test_capture_failure_is_a_dialog(self):
        ui = make_ui()
        ui.viewer = StubViewer(error="The viewer window is minimised.")
        asyncio.run(ui.export_image(None))
        assert ui.errors[-1] == ("Export Image Error", "The viewer window is minimised.")

    def test_write_failure_is_a_dialog(self, tmp_path):
        ui = make_ui()
        ui.viewer = StubViewer(np.zeros((2, 2, 3), np.uint8))

        async def ask(suggested):
            return tmp_path / "picture.jpg"
        ui._ask_image_path = ask
        asyncio.run(ui.export_image(None))
        assert ui.errors[-1][0] == "Export Image Error"


# ----------------------------------------------------------------------
# ESP cube export (headless)
# ----------------------------------------------------------------------
import dataclasses  # noqa: E402

import cubeFile as cf  # noqa: E402
import espIntegrals as ei  # noqa: E402


@pytest.fixture(scope="module")
def water_results():
    return esp.ESPSurfaceAnalyser().run(str(fchk_path("water_631gd_6d")), [0.001], spacing=0.4)


class TestCubeExport:
    def test_results_keep_the_density_grid(self, water_results):
        r = water_results
        assert r.density_grid.shape == tuple(r.grid_shape) and len(r.grid_origin) == 3
        assert r.cubegen_path is None  # the built-in engine ran

    def test_default_names(self, tmp_path):
        assert esp.default_esp_cube_name("C:/runs/water.fchk") == "water_ESP.cube"
        assert esp.density_cube_path(tmp_path / "water_ESP.cube") == tmp_path / "water_density.cube"
        assert esp.density_cube_path(tmp_path / "mine.cube") == tmp_path / "mine_density.cube"

    def test_both_cubes_on_the_run_grid(self, water_results, tmp_path):
        analyser = esp.ESPSurfaceAnalyser()
        seen = []
        paths, notes = analyser.write_cubes(water_results, tmp_path / "water_ESP.cube", with_density=True,
                                            progress=seen.append)
        assert [p.name for p in map(__import__("pathlib").Path, paths)] == ["water_ESP.cube", "water_density.cube"]
        assert notes == [] and seen and seen[-1] == pytest.approx(1.0)
        v, rho = (cf.read_cube(p) for p in paths)
        r = water_results
        for cube in (v, rho):
            assert cube.values.shape == tuple(r.grid_shape) and cube.spacing == pytest.approx(r.spacing)
            assert np.allclose(cube.origin, r.grid_origin, atol=1e-6)
            assert np.array_equal(cube.atomic_numbers, r.wavefunction.atomic_numbers)
        assert np.allclose(rho.values, r.density_grid, rtol=1e-5, atol=1e-12)
        # ESP values: a direct evaluation at sampled grid points, to the file's %13.5E precision.
        rng = np.random.default_rng(5)
        idx = np.column_stack([rng.integers(n, size=25) for n in r.grid_shape])
        points = r.grid_origin + idx * r.spacing
        direct = ei.electrostatic_potential(r.wavefunction, points, ei.build_hermite_batches(r.wavefunction))
        assert np.allclose(v.values[idx[:, 0], idx[:, 1], idx[:, 2]], direct, rtol=2e-5, atol=1e-10)

    def test_density_cube_is_optional_and_the_suffix_is_added(self, water_results, tmp_path):
        paths, _ = esp.ESPSurfaceAnalyser().write_cubes(water_results, tmp_path / "only_esp", with_density=False)
        assert [str(p) for p in paths] == [str(tmp_path / "only_esp.cube")]
        assert not (tmp_path / "only_esp_density.cube").exists()

    def test_a_grid_point_on_a_nucleus_stays_finite(self, water_results):
        r = water_results
        nucleus = r.wavefunction.coordinates[0]
        origin = nucleus - r.spacing * np.floor((nucleus - r.grid_origin) / r.spacing)
        shifted = dataclasses.replace(r, grid_origin=origin)
        values, notes = esp.ESPSurfaceAnalyser().esp_on_grid(shifted)
        assert np.all(np.isfinite(values))
        assert len(notes) == 1 and "nucleus" in notes[0]

    def test_failing_cubegen_falls_back_for_the_cube(self, water_results, monkeypatch):
        def broken(*args, **kwargs):
            raise RuntimeError("cubegen failed: test")
        monkeypatch.setattr(esp, "cubegen_potential", broken)
        r = dataclasses.replace(water_results, engine=esp.ENGINE_CUBEGEN, cubegen_path="C:/fake/cubegen.exe")
        values, notes = esp.ESPSurfaceAnalyser().esp_on_grid(r)
        assert np.all(np.isfinite(values)) and "cubegen failed" in notes[0] and "built-in" in notes[0]

    def test_results_without_a_grid_are_refused(self, water_results):
        r = dataclasses.replace(water_results, grid_origin=None, grid_shape=None, density_grid=None)
        with pytest.raises(ValueError, match="grid"):
            esp.ESPSurfaceAnalyser().esp_on_grid(r)

    def test_unwritable_destination_is_a_value_error(self, water_results, tmp_path):
        (tmp_path / "file").write_text("x")
        with pytest.raises(ValueError, match="Could not write"):
            esp.ESPSurfaceAnalyser().write_cubes(water_results, tmp_path / "file" / "a.cube", with_density=False)


# ----------------------------------------------------------------------
# The tab in the 3D viewer: structure, display and cube export
# ----------------------------------------------------------------------
@pytest.fixture
def water_ui(tmp_path):
    ui = make_ui()
    ui.fchk_input.value = str(fchk_path("water_631gd_6d"))
    ui.output_dir_input.value = str(tmp_path)
    ui.isovalue_input.value = "0.001"
    ui.spacing_input.value = "0.4"
    asyncio.run(ui.workflow(None))
    assert not ui.errors
    return ui


class TestTabDisplay:
    def test_an_empty_viewer_gets_the_molecule_without_a_question(self, water_ui):
        water_ui.viewer = StubViewer()
        asyncio.run(water_ui.show_surface(None))
        assert water_ui.questions == [] and len(water_ui.viewer.loaded) == 1
        assert [s for s, _ in water_ui.viewer.loaded[0]] == ["O", "H", "H"]
        assert water_ui.viewer.esp_layer.get() is not None and water_ui.viewer.background == (0.0, 0.0, 0.0)

    def test_a_fresh_viewers_placeholder_molecule_is_not_asked_about(self, water_ui):
        """Regression: a new MolecularViewer holds a built-in 5-atom placeholder in
        molecule_data (no frames loaded); it must not trigger the replace question."""
        from espSurfaceLayer import ESPSurfaceLayer
        from molecularViewer import MolecularViewer
        viewer = MolecularViewer()
        viewer.esp_layer = ESPSurfaceLayer()
        assert viewer.frames == [] and len(viewer.get_current_frame_data()) == 5  # the placeholder
        water_ui.viewer = viewer
        asyncio.run(water_ui._settle_structure())
        assert water_ui.questions == []
        assert [s for s, _ in viewer.get_current_frame_data()] == ["O", "H", "H"]

    def test_the_same_molecule_is_kept_without_a_question(self, water_ui):
        water_ui.viewer = StubViewer(frame=water_ui.fchk_frame())
        asyncio.run(water_ui.show_surface(None))
        assert water_ui.questions == [] and water_ui.viewer.loaded == []

    def test_a_different_structure_is_asked_about_once_per_run(self, water_ui):
        other = [("C", (0.0, 0.0, 0.0))]
        water_ui.viewer = StubViewer(frame=other)
        water_ui.answer = False
        asyncio.run(water_ui.show_surface(None))
        assert water_ui.questions == ["Show the ESP molecule?"] and water_ui.viewer.frame == other
        assert water_ui.viewer.esp_layer.get() is not None  # drawn in the fchk coordinates anyway
        asyncio.run(water_ui.show_surface(None))
        assert len(water_ui.questions) == 1
        water_ui.answer = True
        asyncio.run(water_ui.workflow(None))  # a new run, with the surface on screen: asked again
        assert len(water_ui.questions) == 2 and len(water_ui.viewer.loaded) == 1

    def test_results_controls_follow_only_a_shown_surface(self, water_ui):
        water_ui.viewer = StubViewer()
        water_ui.opacity_selection.value = "Transparent (50 %)"
        water_ui.on_display_change(None)
        assert water_ui.viewer.esp_layer.get() is None  # not shown yet: nothing drawn
        asyncio.run(water_ui.show_surface(None))
        assert np.allclose(water_ui.viewer.esp_layer.get().rgba[:, 3], 0.5)
        water_ui.unit_selection.value = "kcal/mol"
        water_ui.on_display_change(None)
        assert water_ui.viewer.esp_layer.colour_unit == "kcal/mol"

    def test_clear_removes_the_surface(self, water_ui):
        water_ui.viewer = StubViewer()
        asyncio.run(water_ui.show_surface(None))
        water_ui.clear_output(None)
        assert water_ui.viewer.esp_layer.get() is None and water_ui.results is None


class TestCubeButton:
    def test_without_results(self):
        ui = make_ui()
        asyncio.run(ui.export_cubes(None))
        assert ui.errors[-1][0] == "No Results"

    def test_suggests_the_output_folder_and_writes_both_cubes(self, water_ui, tmp_path):
        asked = []

        async def ask(suggested):
            asked.append(suggested)
            return tmp_path / "chosen" / "my_ESP.cube"
        water_ui._ask_cube_path = ask
        asyncio.run(water_ui.export_cubes(None))
        assert not water_ui.errors
        assert asked == [str(tmp_path / "water_631gd_6d_ESP.cube")]
        assert (tmp_path / "chosen" / "my_ESP.cube").exists()
        assert (tmp_path / "chosen" / "my_density.cube").exists()
        assert "my_ESP.cube" in water_ui.status_label.text and water_ui.cube_button.enabled
        assert water_ui.progress_bar.value == 100

    def test_density_switch_off_and_cancel(self, water_ui, tmp_path):
        water_ui.density_cube_switch.value = False

        async def ask(suggested):
            return tmp_path / "v.cube"
        water_ui._ask_cube_path = ask
        asyncio.run(water_ui.export_cubes(None))
        assert (tmp_path / "v.cube").exists() and not (tmp_path / "v_density.cube").exists()

        async def cancel(suggested):
            return None
        water_ui._ask_cube_path = cancel
        water_ui.status_label.text = "unchanged"
        asyncio.run(water_ui.export_cubes(None))
        assert water_ui.status_label.text == "unchanged" and not water_ui.errors

    def test_write_failure_is_a_dialog(self, water_ui, tmp_path):
        (tmp_path / "file").write_text("x")

        async def ask(suggested):
            return tmp_path / "file" / "v.cube"
        water_ui._ask_cube_path = ask
        asyncio.run(water_ui.export_cubes(None))
        assert water_ui.errors[-1][0] == "ESP Cube Error" and water_ui.cube_button.enabled


def test_the_module_does_not_import_the_viewer():
    import subprocess
    import sys
    from pathlib import Path
    code = ("import sys, electrostaticSurfacePotential; "
            "print('molecularViewer' in sys.modules, 'espViewer' in sys.modules)")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(Path(esp.__file__).parent), check=True)
    assert out.stdout.strip() == "False False"
