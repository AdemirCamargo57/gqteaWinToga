"""UI plumbing of the Gaussian Fukui input builder, headless (stub widgets, no window)."""
import asyncio
import dataclasses
import math
from pathlib import Path
from types import SimpleNamespace

import pytest

import fukuiInput as fi
import fukuiInputBuilder as fb
from fukui_fixtures import N

INPUTS = Path(__file__).resolve().parent / "data" / "fukui" / "inputs"

TEXT_WIDGETS = ("stem_input", "output_dir_input", "memory_input", "nproc_input", "method_input", "basis_input",
                "keywords_input", "opt_input", "geometry_input", "charge_input", "multiplicity_input",
                "plus_mult_input", "minus_mult_input", "p_input", "q_input", "preview_output",
                "messages_output")


def make_ui(tmp_path):
    ui = fb.FukuiInputBuilderUI.__new__(fb.FukuiInputBuilderUI)
    ui.builder = fi.FukuiInputBuilder()
    ui.job = None
    for name in TEXT_WIDGETS:
        setattr(ui, name, SimpleNamespace(value=""))
    ui.stem_input.value = "phenol"
    ui.output_dir_input.value = str(tmp_path)
    ui.memory_input.value = fi.DEFAULT_MEMORY
    ui.nproc_input.value = "4"
    ui.method_input.value = fi.DEFAULT_METHOD
    ui.basis_input.value = fi.DEFAULT_BASIS
    ui.keywords_input.value = fi.DEFAULT_KEYWORDS
    ui.opt_input.value = fi.DEFAULT_OPT_KEYWORD
    ui.charge_input.value = "0"
    ui.p_input.value = ui.q_input.value = "1"
    ui.print_level_selection = SimpleNamespace(value=fi.DEFAULT_PRINT_LEVEL)
    ui.guess_read_switch = SimpleNamespace(value=True)
    ui.states_label = SimpleNamespace(text="")
    ui.route_label = SimpleNamespace(text="")
    ui.status_label = SimpleNamespace(text="")
    ui.tabs = SimpleNamespace(current_tab=0)
    ui.dialogs, ui.answers = [], []

    async def record(title, message):
        ui.dialogs.append((title, message))

    async def confirm(title, message):
        ui.dialogs.append((title, message))
        return ui.answers.pop(0) if ui.answers else True
    ui._show_error = ui._show_info = record
    ui._confirm = confirm
    return ui


def load_phenol(ui, ions=True):
    ui.apply_geometry_file(INPUTS / "phenol_N.gjf")
    if ions:
        ui.plus_mult_input.value = ui.minus_mult_input.value = "2"


@pytest.mark.parametrize("windows,linux", [(r"C:\jobs\my mol", "/mnt/c/jobs/my mol"),
                                           (r"D:\x", "/mnt/d/x"), ("/home/u/jobs", "/home/u/jobs")])
def test_wsl_path(windows, linux):
    assert fb.wsl_path(windows) == linux


class TestGeometryFiles:
    def test_gaussian_input_gives_charge_and_multiplicity(self):
        text, charge, mult = fb.read_geometry_file(INPUTS / "phenol_Np1.gjf")
        assert (charge, mult) == (-1, 2) and len(fi.parse_geometry(text)) == 13

    def test_fchk_gives_angstrom_coordinates(self):
        text, charge, mult = fb.read_geometry_file(N)
        atoms = fi.parse_geometry(text)
        reference = fi.parse_geometry(fb.read_geometry_file(INPUTS / "phenol_N.gjf")[0])
        assert (charge, mult) == (0, 1) and [a[0] for a in atoms] == [a[0] for a in reference]
        # Gaussian's standard orientation may rotate the frame: compare a distance, not coordinates.
        assert math.dist(atoms[0][1:], atoms[1][1:]) == pytest.approx(math.dist(reference[0][1:],
                                                                               reference[1][1:]), abs=1e-5)

    def test_xyz_has_no_charge(self, tmp_path):
        path = tmp_path / "w.xyz"
        path.write_text("3\nwater\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692\n")
        text, charge, mult = fb.read_geometry_file(path)
        assert charge is None and mult is None and len(fi.parse_geometry(text)) == 3

    def test_unsupported_file(self, tmp_path):
        path = tmp_path / "x.log"
        path.write_text("")
        with pytest.raises(ValueError, match="log"):
            fb.read_geometry_file(path)

    def test_apply_fills_the_form_but_not_the_ion_spins(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.stem_input.value = fi.DEFAULT_STEM
        ui.output_dir_input.value = ""
        ui.plus_mult_input.value = ui.minus_mult_input.value = "4"   # from a previous molecule
        ui.apply_geometry_file(INPUTS / "phenol_Np1.gjf")
        assert ui.charge_input.value == "-1" and ui.multiplicity_input.value == "2"
        assert ui.plus_mult_input.value == "" and ui.minus_mult_input.value == ""
        assert ui.stem_input.value == "phenol_Np1" and ui.output_dir_input.value == str(INPUTS)

    def test_apply_from_an_n_fchk_drops_the_state_suffix(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.stem_input.value = ""
        ui.apply_geometry_file(N)
        assert ui.stem_input.value == "phenol" and ui.output_dir_input.value == str(tmp_path)

    def test_apply_keeps_a_typed_job_name(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.stem_input.value = "mine"
        load_phenol(ui)
        assert ui.stem_input.value == "mine"


class TestCollect:
    def test_every_spec_field_is_read_from_a_widget(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        ui.stem_input.value, ui.memory_input.value, ui.nproc_input.value = "ph", "3GB", "3"
        ui.method_input.value, ui.basis_input.value, ui.keywords_input.value = "PBE0", "def2SVP", "SCF=XQC"
        ui.opt_input.value = "Opt=Tight"
        ui.print_level_selection.value = "#N"
        ui.charge_input.value, ui.multiplicity_input.value = "1", "2"
        ui.p_input.value, ui.q_input.value = "2", "3"
        ui.plus_mult_input.value, ui.minus_mult_input.value = "3", "5"
        ui.guess_read_switch.value = False
        spec = ui._collect_inputs()
        defaults = fi.FukuiJobSpec("", "", "", 0, None, None, None)
        for f in dataclasses.fields(fi.FukuiJobSpec):
            assert getattr(spec, f.name) != getattr(defaults, f.name), f"{f.name} is not read from the window"

    def test_blank_multiplicities_are_not_filled_in(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui, ions=False)
        spec = ui._collect_inputs()
        assert (spec.multiplicity, spec.plus_multiplicity, spec.minus_multiplicity) == (1, None, None)

    def test_blank_charge_is_neutral(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        ui.charge_input.value = ""
        assert ui._collect_inputs().charge == 0

    @pytest.mark.parametrize("widget,value,fragment", [("charge_input", "0.5", "Charge"),
                                                       ("nproc_input", "x", "%NProcShared"),
                                                       ("p_input", "two", "p "),
                                                       ("plus_mult_input", "2.5", "N\\+p")])
    def test_bad_numbers(self, tmp_path, widget, value, fragment):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        getattr(ui, widget).value = value
        with pytest.raises(ValueError, match=fragment):
            ui._collect_inputs()


class TestStatesSummary:
    def test_summary_lists_the_three_states(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        ui.refresh_states(None)
        text = ui.states_label.text
        assert "N+1" in text and "charge -1" in text and "51 electrons" in text and "phenol_Np1.chk" in text
        assert "multiplicity 2" in text and "suggest" not in text
        assert ui.route_label.text == ("job 1:     #P B3LYP/6-31G(d) SCF=Tight NoSymm Opt\n"
                                       "jobs 2-3:  #P B3LYP/6-31G(d) SCF=Tight NoSymm Geom=Check Guess=Read")

    def test_summary_asks_for_missing_spins(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui, ions=False)
        ui.refresh_states(None)
        assert "multiplicity ? (enter it: 51 electrons allow 2, 4, 6, ...)" in ui.states_label.text

    def test_summary_without_geometry(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.refresh_states(None)
        assert "geometry" in ui.states_label.text.lower()

    def test_summary_reports_an_impossible_multiplicity(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        ui.plus_mult_input.value = "1"
        ui.refresh_states(None)
        assert "impossible" in ui.states_label.text


class TestActions:
    def test_preview_fills_the_input_text(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        asyncio.run(ui.preview(None))
        assert ui.preview_output.value.count("--Link1--") == 2 and not ui.dialogs

    def test_preview_asks_for_missing_spins(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui, ions=False)
        asyncio.run(ui.preview(None))
        assert ui.dialogs[0][0] == "Input Error"
        assert "N+1 state: enter its spin multiplicity; it is not assumed" in ui.dialogs[0][1]
        assert "N-1 state" in ui.dialogs[0][1]

    def test_preview_errors_go_to_a_dialog(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        ui.keywords_input.value = "Opt"
        asyncio.run(ui.preview(None))
        assert ui.dialogs and ui.dialogs[0][0] == "Input Error" and "Opt" in ui.dialogs[0][1]

    def test_preview_shows_warnings(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        ui.keywords_input.value = "Freq"
        asyncio.run(ui.preview(None))
        assert "Freq" in ui.messages_output.value

    def test_save_writes_both_files_and_says_how_to_run_them(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        asyncio.run(ui.save_input(None))
        assert (tmp_path / "phenol_fukui.gjf").is_file() and (tmp_path / "phenol_formchk.sh").is_file()
        assert "g16 < phenol_fukui.gjf > phenol_fukui.log" in ui.messages_output.value
        assert "bash phenol_formchk.sh" in ui.messages_output.value

    def test_existing_files_are_kept_unless_confirmed(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        (tmp_path / "phenol_fukui.gjf").write_text("mine")
        ui.answers = [False]
        asyncio.run(ui.save_input(None))
        assert (tmp_path / "phenol_fukui.gjf").read_text() == "mine"
        ui.answers = [True]
        asyncio.run(ui.save_input(None))
        assert "--Link1--" in (tmp_path / "phenol_fukui.gjf").read_text()
