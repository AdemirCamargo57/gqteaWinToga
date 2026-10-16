"""Tests for the Fukui function analyser, Multiwfn interop and the UI plumbing."""
import asyncio
import inspect
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

import fukuiFunctionAnalyze as fk
import gaussianWavefunction as gw
from fukui_fixtures import CDFT, MINUS, N, PLUS, REFERENCE_DIR

CDFT_TEXT = CDFT.read_text()


class TestFormulas:
    def test_fukui_definitions(self):
        f = fk.fukui_from_densities(np.array([1.0]), np.array([1.4]), np.array([0.7]))
        assert f["f+"][0] == pytest.approx(0.4) and f["f-"][0] == pytest.approx(0.3)
        assert f["f0"][0] == pytest.approx(0.35) and f["DD"][0] == pytest.approx(0.1)

    def test_degenerate_states_divide_by_p_and_q(self):
        f = fk.fukui_from_densities(np.array([1.0]), np.array([1.6]), np.array([0.4]), p=2, q=3)
        assert f["f+"][0] == pytest.approx(0.3) and f["f-"][0] == pytest.approx(0.2)

    def test_global_indices_reproduce_multiwfn(self):
        # E(N), E(N+1), E(N-1), E_HOMO(N) and the expected values: phenol_CDFT.txt
        g = fk.GlobalIndices(-307.464870, -307.383627, -307.163447, -0.218914)
        assert g.ip == pytest.approx(0.301423, abs=2e-6)
        assert g.ea == pytest.approx(-0.081243, abs=2e-6)
        assert g.chi == pytest.approx(0.110090, abs=2e-6) and g.mu == pytest.approx(-0.110090, abs=2e-6)
        assert g.eta == pytest.approx(0.382666, abs=2e-6)
        assert g.softness == pytest.approx(2.613246, abs=2e-5)
        assert g.omega == pytest.approx(0.015836, abs=2e-6)
        assert g.nucleophilicity == pytest.approx(0.116284, abs=2e-6)
        assert g.eta * fk.HARTREE_TO_EV == pytest.approx(10.4129, abs=1e-4)

    def test_non_positive_gap_gives_nan_not_an_exception(self):
        g = fk.GlobalIndices(-1.0, -1.2, -0.9)  # IP 0.1 < EA 0.2
        assert g.eta < 0 and np.isnan(g.softness) and np.isnan(g.omega)

    @pytest.mark.parametrize("method,expected", [("RB3LYP", "B3LYP"), ("UB3LYP", "B3LYP"),
                                                 ("ROB3LYP", "B3LYP"), ("RHF", "HF"), ("UHF", "HF")])
    def test_reference_prefix(self, method, expected):
        assert fk.strip_reference_prefix(method) == expected

    @pytest.mark.parametrize("basis,diffuse", [("6-31G(d)", False), ("6-311+G(d,p)", True),
                                               ("aug-cc-pVTZ", True), ("def2-TZVP", False)])
    def test_diffuse_basis_detection(self, basis, diffuse):
        assert fk.has_diffuse_functions(basis) is diffuse


class TestParsing:
    def test_blank_is_default(self):
        assert fk.parse_positive_float("", 0.005, "x", 1.0) == 0.005

    @pytest.mark.parametrize("text", ["abc", "0", "-1", "2"])
    def test_bad_values(self, text):
        with pytest.raises(ValueError):
            fk.parse_positive_float(text, 0.005, "Fukui isovalue", 1.0)


class TestStates:
    @pytest.fixture(scope="class")
    def states(self):
        return gw.read_fchk(N), gw.read_fchk(PLUS), gw.read_fchk(MINUS)

    def test_new_fchk_fields(self, states):
        n, plus, _ = states
        assert n.method == "RB3LYP" and plus.method == "UB3LYP" and n.basis_set == "6-31G(d)"
        assert n.total_energy == pytest.approx(-307.464870, abs=1e-6)
        assert (n.n_alpha, n.n_beta) == (25, 25) and n.other_densities == []
        assert n.alpha_orbital_energies[n.n_alpha - 1] == pytest.approx(-0.218914, abs=1e-6)

    def test_phenol_triplet_is_valid(self, states):
        p, q, notes = fk.validate_states(*states)
        assert (p, q, notes) == (1, 1, [])

    def test_swapped_files_are_explained(self, states):
        n, plus, minus = states
        with pytest.raises(ValueError, match="swapped"):
            fk.validate_states(n, minus, plus)

    def test_moved_geometry_is_rejected(self, states):
        n, plus, minus = states
        moved = replace(plus, coordinates=plus.coordinates + 0.01)
        with pytest.raises(ValueError, match="single points"):
            fk.validate_states(n, moved, minus)

    def test_other_method_is_rejected(self, states):
        n, plus, minus = states
        with pytest.raises(ValueError, match="same method"):
            fk.validate_states(n, replace(plus, method="UPBE0"), minus)

    def test_other_basis_is_rejected(self, states):
        n, plus, minus = states
        with pytest.raises(ValueError, match="basis set"):
            fk.validate_states(n, replace(plus, shells=plus.shells[:-1]), minus)

    def test_unusual_multiplicity_is_a_note(self, states):
        n, plus, minus = states
        _, _, notes = fk.validate_states(n, replace(plus, multiplicity=4), minus)
        assert "multiplicity 4" in notes[0]

    def test_phenol_frontier_orbitals_are_not_degenerate(self, states):
        assert fk.frontier_degeneracy(states[0]) == (1, 1)
        assert fk.frontier_degeneracy(states[1]) is None  # open shell


class TestCDFT:
    def test_parse_condensed_table_and_globals(self):
        rows, g = fk.parse_cdft_text(CDFT_TEXT)
        assert len(rows) == 13 and rows[0].symbol == "C" and rows[11].symbol == "O"
        c1 = rows[0]
        assert (c1.q_n, c1.q_plus, c1.q_minus, c1.f_minus, c1.f_plus, c1.f_zero, c1.cdd) == \
               (-0.0587, -0.1185, 0.0852, 0.1439, 0.0598, 0.1018, -0.0841)
        assert g["e_n"] == -307.464870 and g["eta"] == 0.382666 and g["ip"] == 0.301423

    def test_older_vertical_ip_label_is_accepted(self):
        _, g = fk.parse_cdft_text(CDFT_TEXT.replace("First vertical IP", "Vertical IP"))
        assert g["ip"] == 0.301423

    def test_phenol_rows_are_self_consistent(self):
        rows, _ = fk.parse_cdft_text(CDFT_TEXT)
        assert fk.condensed_inconsistencies(rows) == []

    def test_inconsistent_row_is_found(self):
        rows, _ = fk.parse_cdft_text(CDFT_TEXT)
        rows[2] = replace(rows[2], f_plus=rows[2].f_plus + 0.01)
        assert fk.condensed_inconsistencies(rows) == [3]

    def test_not_a_cdft_file(self):
        with pytest.raises(ValueError, match="main function 22"):
            fk.parse_cdft_text("hello\nworld\n")

    def test_rows_must_match_the_molecule(self):
        rows, _ = fk.parse_cdft_text(CDFT_TEXT)
        n = gw.read_fchk(N)
        fk.check_rows_match(rows, n)
        with pytest.raises(ValueError, match="lists 12 atoms"):
            fk.check_rows_match(rows[:-1], n)

    def test_script_adds_a_blank_line_only_without_settings(self):
        assert fk.multiwfn_cdft_script(True).startswith("22\n2\nN.fch\nN+1.fch\nN-1.fch\n")
        assert fk.multiwfn_cdft_script(False).startswith("\n22\n")


@pytest.fixture(scope="module")
def phenol():
    seen = []
    results = fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), spacing=0.3,
                                     condensed_source=fk.CONDENSED_IMPORT, cdft_path=str(CDFT),
                                     progress=lambda f, m: seen.append(f))
    return results, seen


class TestPhenol:
    def test_progress_is_monotone_to_one(self, phenol):
        _, seen = phenol
        assert seen[-1] == 1.0 and all(b >= a - 1e-12 for a, b in zip(seen, seen[1:]))

    def test_grid_matches_multiwfn(self):
        data = np.loadtxt(REFERENCE_DIR / "phenol_grid_sample.txt", comments="#")
        states = [gw.read_fchk(p) for p in (N, PLUS, MINUS)]
        mine = fk.fukui_from_densities(*[gw.electron_density(w, data[:, :3]) for w in states])
        for k, key in enumerate(fk.FUNCTIONS):
            assert np.max(np.abs(mine[key] - data[:, 3 + k])) < 2e-5

    def test_integrals_are_close_to_one_and_zero(self, phenol):
        r, _ = phenol
        for key in ("f+", "f-", "f0"):
            assert 0.9 < r.integral(key) < 1.02
        assert abs(r.integral("DD")) < 0.05

    def test_globals_and_condensed_rows(self, phenol):
        r, _ = phenol
        assert r.global_indices.eta == pytest.approx(0.382666, abs=2e-6)
        assert len(r.condensed) == 13 and r.condensed_source == fk.CONDENSED_IMPORT
        assert r.notes == []

    def test_ranking_matches_the_manual_for_phenol(self, phenol):
        # Multiwfn manual 4.7.3: para C1 most favoured for electrophilic attack among ring carbons
        r, _ = phenol
        ranking = dict(fk.FukuiAnalyser.site_ranking(r.condensed, top=13))
        ring = [row.index for row in ranking["electrophilic attack: most negative CDD"] if row.symbol == "C"]
        assert ring[0] == 1

    def test_surface_maxima_of_f_minus_lie_on_ring_and_oxygen(self, phenol):
        r, _ = phenol
        maxima = [e for e in r.extrema["f-"] if e.kind == "maximum"]
        assert maxima and maxima[0].scope == "global"
        assert fk.FukuiAnalyser.nearest_atom_label(r.n_state, maxima[0].position)[0] in "CO"

    def test_lobes_carry_both_signs_for_the_dual_descriptor(self, phenol):
        r, _ = phenol
        lobe = r.lobes["DD"]
        assert lobe is not None and set(np.unique(lobe.values)) == {-r.fukui_isovalue, r.fukui_isovalue}
        assert lobe.triangles.max() < len(lobe.vertices)

    def test_outputs(self, phenol, tmp_path):
        r, _ = phenol
        written = fk.FukuiAnalyser().write_outputs(r, str(tmp_path))
        names = sorted(p.split("\\")[-1].split("/")[-1] for p in written)
        assert names == ["phenol_N_Fukui.txt", "phenol_N_dual.cub", "phenol_N_fminus.cub",
                         "phenol_N_fplus.cub", "phenol_N_fzero.cub"]
        text = (tmp_path / "phenol_N_Fukui.txt").read_text()
        assert "not predictions of reaction outcome" in text and "SITE RANKING" in text
        assert "Hardness" in text and "SURFACE EXTREMA OF f-" in text

    def test_cube_round_trip(self, phenol, tmp_path):
        import cubeFile as cf
        r, _ = phenol
        fk.FukuiAnalyser().write_outputs(r, str(tmp_path))
        cube = cf.read_cube(tmp_path / "phenol_N_fminus.cub")
        assert np.allclose(cube.values, r.grids["f-"], rtol=1e-5, atol=1e-10)


class TestCondensedSources:
    def test_missing_multiwfn_continues_with_a_note(self, monkeypatch):
        monkeypatch.setattr(fk, "find_multiwfn", lambda explicit=None: None)
        r = fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), spacing=0.4)
        assert r.condensed == [] and "not found" in r.notes[-1]

    def test_failing_multiwfn_continues_with_a_note(self, monkeypatch):
        monkeypatch.setattr(fk, "find_multiwfn", lambda explicit=None: "C:/fake/Multiwfn.exe")

        def broken(*args, **kwargs):
            raise RuntimeError("Multiwfn finished without writing CDFT.txt (exit 24): test")
        monkeypatch.setattr(fk, "run_multiwfn_cdft", broken)
        r = fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), spacing=0.4)
        assert r.condensed == [] and "without writing CDFT.txt" in r.notes[-1]

    @pytest.mark.skipif(fk.find_multiwfn() is None, reason="Multiwfn not installed")
    def test_multiwfn_engine_reproduces_the_reference(self):
        text = fk.run_multiwfn_cdft(fk.find_multiwfn(), [str(N), str(PLUS), str(MINUS)])
        assert fk.parse_cdft_text(text)[0] == fk.parse_cdft_text(CDFT_TEXT)[0]

    def test_unknown_source(self):
        with pytest.raises(ValueError):
            fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), condensed_source="magic")


def make_ui():
    """The UI class without a window: stub widgets, as in test_electrostaticSurfacePotential."""
    ui = fk.FukuiFunctionAnalyzeUI.__new__(fk.FukuiFunctionAnalyzeUI)
    ui.analyser = fk.FukuiAnalyser()
    ui.results = None
    ui.written = []
    ui.viewer = None
    ui._syncing = False
    for name in ("n_input", "plus_input", "minus_input", "output_dir_input", "fukui_iso_input",
                 "surface_iso_input", "spacing_input", "multiwfn_input", "cdft_input"):
        setattr(ui, name, SimpleNamespace(value=""))
    ui.condensed_selection = SimpleNamespace(value=fk.CONDENSED_NONE)
    ui.cubes_switch = SimpleNamespace(value=False)
    ui.function_selection = SimpleNamespace(value="f-")
    ui.display_selection = SimpleNamespace(value=fk.DISPLAY_LOBES)
    ui.extrema_switch = SimpleNamespace(value=True)
    ui.background_selection = SimpleNamespace(value="Black")
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


def fill(ui, tmp_path):
    ui.n_input.value, ui.plus_input.value, ui.minus_input.value = str(N), str(PLUS), str(MINUS)
    ui.output_dir_input.value = str(tmp_path)
    ui.spacing_input.value = "0.4"


class StubViewer:
    def __init__(self):
        self.surfaces = []

    def set_molecule(self, symbols, xyz):
        self.symbols = list(symbols)

    def set_display_options(self, **kwargs):
        self.options = kwargs

    def set_surface(self, *args, **kwargs):
        self.surfaces.append((args, kwargs))

    def clear_surface(self):
        self.surfaces.append(None)

    def viewer_window_is_open(self):
        return False


class TestUI:
    def test_missing_files_are_an_input_error(self):
        ui = make_ui()
        asyncio.run(ui.workflow(None))
        assert ui.errors and ui.errors[0][0] == "Input Error" and "N state" in ui.errors[0][1]

    def test_import_without_a_file_is_an_input_error(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        ui.condensed_selection.value = fk.CONDENSED_IMPORT
        with pytest.raises(ValueError, match="CDFT.txt"):
            ui._collect_inputs()

    def test_every_run_parameter_is_supplied_by_the_ui(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        wanted = set(inspect.signature(fk.FukuiAnalyser.run).parameters) - {"self", "progress"}
        assert wanted == set(ui._collect_inputs())

    def test_workflow_writes_the_report_and_switches_tab(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        assert not ui.errors and (tmp_path / "phenol_N_Fukui.txt").exists()
        assert not list(tmp_path.glob("*.cub"))  # cubes switch off
        assert ui.tabs.current_tab == 1 and ui.run_button.enabled

    def test_status_line_names_the_report_without_its_folder(self, tmp_path):
        # A full path in the status line widens the window, and Toga never shrinks it back.
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        assert "phenol_N_Fukui.txt" in ui.status_label.text
        assert str(tmp_path) not in ui.status_label.text

    def test_viewer_gets_lobes_or_the_mapped_surface(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        ui.viewer = StubViewer()
        ui._push_to_viewer()
        args, kwargs = ui.viewer.surfaces[-1]
        assert "blue positive, red negative" in kwargs["caption"] and kwargs["half_range"] == 0.005
        ui.display_selection.value = fk.DISPLAY_MAPPED
        ui._push_to_viewer()
        args, kwargs = ui.viewer.surfaces[-1]
        assert kwargs["caption"].startswith("f- (a.u.)") and args[5]  # extrema markers
        assert ui.default_image_name() == "phenol_N_fminus.png"

    def test_clear_resets_results(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        ui.clear_output(None)
        assert ui.results is None
        asyncio.run(ui.show_surface(None))
        assert ui.errors[-1][0] == "No Results"
