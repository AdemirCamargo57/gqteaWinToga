"""Tests for the Fukui function analyser (native Hirshfeld condensed indices) and the UI plumbing."""
import asyncio
import inspect
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

import fukuiFunctionAnalyze as fk
import gaussianWavefunction as gw
from fukui_fixtures import CDFT, MINUS, N, PLUS, REFERENCE_DIR



def multiwfn_condensed():
    """Reference only: Multiwfn 2026.7.15 CDFT.txt rows q(N) q(N+1) q(N-1) f- f+ f0 CDD (4 decimals)."""
    lines = CDFT.read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.split()[:2] == ["Atom", "q(N)"])
    rows = []
    for line in lines[start + 1:]:
        parts = line.split()
        if len(parts) < 8:
            break
        rows.append([float(v) for v in parts[-7:]])
    return np.array(rows)


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


class TestCondensed:
    def test_definitions_from_charges(self):
        r = fk.condensed_rows(["C"], np.array([-0.0587]), np.array([-0.1185]), np.array([0.0852]))[0]
        assert (r.index, r.symbol, r.q_n, r.q_plus, r.q_minus) == (1, "C", -0.0587, -0.1185, 0.0852)
        assert r.f_minus == pytest.approx(0.1439) and r.f_plus == pytest.approx(0.0598)
        assert r.f_zero == pytest.approx(0.10185) and r.cdd == pytest.approx(-0.0841)

    def test_degenerate_states_divide_by_p_and_q(self):
        r = fk.condensed_rows(["C"], np.array([0.0]), np.array([-0.4]), np.array([0.3]), p=2, q=3)[0]
        assert r.f_plus == pytest.approx(0.2) and r.f_minus == pytest.approx(0.1)

    def test_multiwfn_is_no_longer_driven(self):
        for name in ("run_multiwfn_cdft", "find_multiwfn", "parse_cdft_text", "subprocess", "CONDENSED_SOURCES"):
            assert not hasattr(fk, name), name


@pytest.fixture(scope="module")
def phenol():
    seen = []
    results = fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), spacing=0.3,
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
        assert len(r.condensed) == 13 and r.condensed_source.startswith("Hirshfeld")
        assert r.notes == []

    def test_condensed_indices_reproduce_multiwfn(self, phenol):
        r, _ = phenol
        mine = np.array([[c.q_n, c.q_plus, c.q_minus, c.f_minus, c.f_plus, c.f_zero, c.cdd] for c in r.condensed])
        assert np.max(np.abs(mine - multiwfn_condensed())) < 3e-4  # Multiwfn prints 4 decimals

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


class TestCubeGrid:
    """Exported cubes on their own, finer grid, computed from density-matrix differences."""

    def test_same_spacing_reuses_the_analysis_grids(self, phenol):
        r, _ = phenol
        analyser = fk.FukuiAnalyser()
        assert analyser.cube_grids(r, None) is r.grids
        assert analyser.cube_grids(r, r.spacing) is r.grids

    def test_density_matrix_differences_reproduce_the_density_differences(self, phenol):
        r, _ = phenol
        f_plus, f_minus = fk.FukuiAnalyser().cube_pair(r, r.spacing)
        assert f_plus.shape == r.grids["f+"].shape
        assert np.max(np.abs(f_plus - r.grids["f+"])) < 1e-10
        assert np.max(np.abs(f_minus - r.grids["f-"])) < 1e-10

    def test_finer_cube_is_exact_at_its_points_and_covers_the_box(self, phenol, tmp_path):
        import cubeFile as cf
        r, _ = phenol
        seen = []
        fk.FukuiAnalyser().write_outputs(r, str(tmp_path), cube_spacing=0.25,
                                         progress=lambda f, m: seen.append(f))
        assert seen[-1] == 1.0 and all(b >= a - 1e-12 for a, b in zip(seen, seen[1:]))
        cube = cf.read_cube(tmp_path / "phenol_N_dual.cub")
        assert cube.spacing == 0.25 and np.allclose(cube.origin, r.origin, atol=1e-6)
        assert cube.values.shape == r.cube_shape == fk.FukuiAnalyser.cube_shape(r, 0.25)
        extent = (np.array(r.grids["f+"].shape) - 1) * r.spacing
        assert np.all((np.array(cube.values.shape) - 1) * 0.25 >= extent - 1e-9)
        # Spot-check against the densities of the three states at those grid points.
        rng = np.random.default_rng(5)
        index = np.column_stack([rng.integers(0, n, 40) for n in cube.values.shape])
        points = r.origin + 0.25 * index
        direct = fk.fukui_from_densities(*[gw.electron_density(w, points) for w in r.states])["DD"]
        assert np.allclose(cube.values[tuple(index.T)], direct, rtol=2e-5, atol=1e-11)
        text = (tmp_path / "phenol_N_Fukui.txt").read_text()
        nx, ny, nz = r.cube_shape
        assert f"# cube_spacing_bohr      0.25  ({nx} x {ny} x {nz} points" in text

    def test_no_cubes_means_no_cube_line(self, phenol, tmp_path):
        r, _ = phenol
        fk.FukuiAnalyser().write_outputs(r, str(tmp_path), write_cubes=False, cube_spacing=0.25)
        assert r.cube_spacing is None and not list(tmp_path.glob("*.cub"))
        assert "cube_spacing_bohr" not in (tmp_path / "phenol_N_Fukui.txt").read_text()

    def test_too_many_points_asks_for_a_larger_spacing(self, phenol, tmp_path, monkeypatch):
        r, _ = phenol
        monkeypatch.setattr(fk, "MAX_CUBE_POINTS", 1000)
        with pytest.raises(ValueError, match="cube grid spacing of at least"):
            fk.FukuiAnalyser().write_outputs(r, str(tmp_path), cube_spacing=0.25)
        assert not list(tmp_path.glob("*"))  # refused before anything was written


class TestLobeParts:
    def test_positive_and_negative_parts_split_the_mesh(self, phenol):
        r, _ = phenol
        lobe = r.lobes["DD"]
        plus, minus = fk.lobe_part(lobe, 1.0), fk.lobe_part(lobe, -1.0)
        assert set(plus.values) == {r.fukui_isovalue} and set(minus.values) == {-r.fukui_isovalue}
        assert len(plus.vertices) + len(minus.vertices) == len(lobe.vertices)
        assert len(plus.triangles) + len(minus.triangles) == len(lobe.triangles)
        for part in (plus, minus):
            assert 0 <= part.triangles.min() and part.triangles.max() < len(part.vertices)
            assert len(part.normals) == len(part.vertices)
        # The kept triangles are the original ones (+ first, then -), renumbered.
        n_plus = len(plus.triangles)
        assert np.array_equal(plus.vertices[plus.triangles], lobe.vertices[lobe.triangles[:n_plus]])
        assert np.array_equal(minus.vertices[minus.triangles], lobe.vertices[lobe.triangles[n_plus:]])

    def test_both_or_missing(self):
        lobe = fk.LobeSurface(np.zeros((3, 3)), np.array([[0, 1, 2]]), np.zeros((3, 3)), np.full(3, 0.005))
        assert fk.lobe_part(lobe, None) is lobe
        assert fk.lobe_part(lobe, -1.0) is None and fk.lobe_part(None, 1.0) is None


class TestCondensedOptions:
    def test_condensed_indices_can_be_skipped(self):
        r = fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), spacing=0.4, condensed=False)
        assert r.condensed == [] and r.condensed_source == "" and r.notes == []

    def test_unsupported_element_continues_with_a_note(self, monkeypatch):
        def refuse(*args, **kwargs):
            raise ValueError("No free-atom density for Fe: Hirshfeld charges are available for H-Ar only.")
        monkeypatch.setattr(fk.hc, "hirshfeld_charges", refuse)
        r = fk.FukuiAnalyser().run(str(N), str(PLUS), str(MINUS), spacing=0.4)
        assert r.condensed == [] and "No free-atom density for Fe" in r.notes[-1]
        assert r.grids["f-"].size > 0  # real-space results are still there


def make_ui():
    """The UI class without a window: stub widgets, as in test_electrostaticSurfacePotential."""
    ui = fk.FukuiFunctionAnalyzeUI.__new__(fk.FukuiFunctionAnalyzeUI)
    ui.analyser = fk.FukuiAnalyser()
    ui.results = None
    ui.written = []
    ui.viewer = None
    ui._syncing = False
    for name in ("n_input", "plus_input", "minus_input", "output_dir_input", "fukui_iso_input",
                 "surface_iso_input", "spacing_input", "cube_spacing_input"):
        setattr(ui, name, SimpleNamespace(value=""))
    ui.condensed_switch = SimpleNamespace(value=False)
    ui.atom_style_selection = SimpleNamespace(value=fk.DEFAULT_ATOM_STYLE)
    ui.atom_numbers_switch = SimpleNamespace(value=False)
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

    def set_atom_view(self, style, show_numbers):
        self.atom_view = (style, show_numbers)

    def viewer_window_is_open(self):
        return False


class TestUI:
    def test_choosing_an_n_file_fills_in_the_ion_files(self, tmp_path):
        import shutil
        for tag, source in (("N", N), ("Np1", PLUS), ("Nm1", MINUS)):
            shutil.copy(source, tmp_path / f"phenol_{tag}.fchk")
        ui = make_ui()
        ui.n_input.value = str(tmp_path / "phenol_N.fchk")
        ui.fill_sibling_states()
        assert ui.plus_input.value == str(tmp_path / "phenol_Np1.fchk")
        assert ui.minus_input.value == str(tmp_path / "phenol_Nm1.fchk")

    def test_auto_fill_keeps_typed_ion_files(self, tmp_path):
        import shutil
        for tag, source in (("N", N), ("Np1", PLUS), ("Nm1", MINUS)):
            shutil.copy(source, tmp_path / f"phenol_{tag}.fchk")
        ui = make_ui()
        ui.n_input.value = str(tmp_path / "phenol_N.fchk")
        ui.plus_input.value = "mine.fchk"
        ui.fill_sibling_states()
        assert ui.plus_input.value == "mine.fchk" and ui.minus_input.value.endswith("phenol_Nm1.fchk")

    def test_missing_files_are_an_input_error(self):
        ui = make_ui()
        asyncio.run(ui.workflow(None))
        assert ui.errors and ui.errors[0][0] == "Input Error" and "N state" in ui.errors[0][1]

    def test_every_run_parameter_is_supplied_by_the_ui(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        wanted = set(inspect.signature(fk.FukuiAnalyser.run).parameters) - {"self", "progress"}
        assert wanted == set(ui._collect_inputs())

    def test_every_output_option_is_supplied_by_the_ui(self):
        ui = make_ui()
        wanted = set(inspect.signature(fk.FukuiAnalyser.write_outputs).parameters) - {
            "self", "results", "output_dir", "progress"}
        assert wanted == set(ui._collect_output_options())

    def test_cube_spacing_blank_is_the_fine_default_and_is_validated(self):
        ui = make_ui()
        assert ui._collect_output_options()["cube_spacing"] == fk.DEFAULT_CUBE_SPACING_BOHR == 0.1
        ui.cube_spacing_input.value = "0.9"
        with pytest.raises(ValueError, match="Cube grid spacing"):
            ui._collect_output_options()

    def test_workflow_writes_cubes_on_the_chosen_grid(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        ui.cubes_switch.value = True
        ui.cube_spacing_input.value = "0.35"
        asyncio.run(ui.workflow(None))
        assert not ui.errors and len(list(tmp_path.glob("*.cub"))) == 4
        assert ui.results.cube_spacing == 0.35 and "cube_spacing_bohr      0.35" in ui.summary_output.value
        assert ui.progress_bar.value == 100

    def test_a_write_failure_keeps_the_analysis(self, tmp_path, monkeypatch):
        ui = make_ui()
        fill(ui, tmp_path)

        def refuse(*args, **kwargs):
            raise ValueError("Could not write the report")
        monkeypatch.setattr(ui.analyser, "write_outputs", refuse)
        asyncio.run(ui.workflow(None))
        assert ui.errors == [("Fukui Output Error", "Could not write the report")]
        assert ui.results is not None and "GLOBAL INDICES" in ui.summary_output.value
        assert ui.run_button.enabled

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
        assert ui.viewer.atom_view == ("Line style", False)  # the default view is pushed too

    def test_display_offers_each_sign_on_its_own(self):
        assert fk.DISPLAY_MODES == (fk.DISPLAY_LOBES, fk.DISPLAY_POSITIVE, fk.DISPLAY_NEGATIVE, fk.DISPLAY_MAPPED)

    def test_positive_or_negative_isosurface_only(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        ui.viewer = StubViewer()
        ui.function_selection.value = "DD"
        ui.display_selection.value = fk.DISPLAY_POSITIVE
        ui._push_to_viewer()
        args, kwargs = ui.viewer.surfaces[-1]
        assert set(args[3]) == {0.005} and kwargs["caption"] == "DD positive isosurface at +0.005 a.u. (blue)"
        assert ui.default_image_name() == "phenol_N_dual_positive.png"
        ui.display_selection.value = fk.DISPLAY_NEGATIVE
        ui._push_to_viewer()
        args, kwargs = ui.viewer.surfaces[-1]
        assert set(args[3]) == {-0.005} and kwargs["caption"] == "DD negative isosurface at -0.005 a.u. (red)"
        assert ui.default_image_name() == "phenol_N_dual_negative.png"

    def test_a_missing_sign_clears_the_surface_and_says_so(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        ui.results.lobes["f-"] = fk.lobe_part(ui.results.lobes["f-"], 1.0)  # no negative part
        ui.viewer = StubViewer()
        ui.display_selection.value = fk.DISPLAY_NEGATIVE
        ui._push_to_viewer()
        assert ui.viewer.surfaces[-1] is None
        assert ui.status_label.text == "f- has no negative isosurface at -0.005 a.u."

    def test_atom_style_and_numbers_reach_the_viewer(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        ui.viewer = StubViewer()
        ui.atom_style_selection.value = "Space-filling (vdW)"
        ui.atom_numbers_switch.value = True
        ui.on_atom_view_change(None)
        assert ui.viewer.atom_view == ("vdW style", True)
        ui.atom_style_selection.value = "Ball and stick (CPK)"
        ui.atom_numbers_switch.value = False
        ui.on_atom_view_change(None)
        assert ui.viewer.atom_view == ("CPK style", False)

    def test_atom_view_applies_without_results(self):
        ui = make_ui()
        ui.viewer = StubViewer()
        ui.atom_style_selection.value = "Ball and stick (CPK)"
        ui.on_atom_view_change(None)
        assert ui.viewer.atom_view == ("CPK style", False)

    def test_three_atom_styles_line_is_the_default(self):
        assert list(fk.ATOM_STYLES) == ["Line", "Ball and stick (CPK)", "Space-filling (vdW)"]
        assert fk.DEFAULT_ATOM_STYLE == "Line"

    def test_clear_resets_results(self, tmp_path):
        ui = make_ui()
        fill(ui, tmp_path)
        asyncio.run(ui.workflow(None))
        ui.clear_output(None)
        assert ui.results is None
        asyncio.run(ui.show_surface(None))
        assert ui.errors[-1][0] == "No Results"
