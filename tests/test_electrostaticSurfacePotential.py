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
    """The UI class without a window: stub widgets, as in test_meanResidenceTime."""
    ui = esp.ElectrostaticSurfacePotentialUI.__new__(esp.ElectrostaticSurfacePotentialUI)
    ui.analyser = esp.ESPSurfaceAnalyser()
    ui.results = None
    ui.extrema_path = None
    ui.viewer = None
    ui._syncing = False
    for name in ("fchk_input", "output_dir_input", "isovalue_input", "spacing_input", "cubegen_input"):
        setattr(ui, name, SimpleNamespace(value=""))
    ui.engine_selection = SimpleNamespace(value=esp.ENGINE_BUILTIN)
    ui.surface_selection = SimpleNamespace(value=None, items=[])
    ui.unit_selection = SimpleNamespace(value="a.u.")
    ui.opacity_selection = SimpleNamespace(value="Opaque")
    ui.extrema_switch = SimpleNamespace(value=True)
    ui.summary_output = SimpleNamespace(value="")
    ui.status_label = SimpleNamespace(text="")
    ui.progress_bar = SimpleNamespace(value=0)
    ui.tabs = SimpleNamespace(current_tab=0)
    for name in ("run_button", "show_button", "clear_button"):
        setattr(ui, name, SimpleNamespace(enabled=True))
    ui.errors = []

    async def record(title, message):
        ui.errors.append((title, message))
    ui._show_error = record
    ui._show_info = record
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
        assert ui.tabs.current_tab == 1 and ui.run_button.enabled
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
