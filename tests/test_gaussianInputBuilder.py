"""The Gaussian 16 input builder window, headless (stub widgets, no window)."""
import asyncio
from dataclasses import fields, replace
from pathlib import Path
from types import SimpleNamespace

import pytest

import fukuiInput as fi
import gaussianInput as gi
import gaussianInputBuilder as gb
import gaussianKeywords as gk

SRC = Path(__file__).resolve().parents[1]
WATER = ("O 0.000000 0.000000 0.117300\n"
         "H 0.000000 0.757200 -0.469200\n"
         "H 0.000000 -0.757200 -0.469200")


class Recorder(SimpleNamespace):
    """A stub widget that counts writes, to check the write-on-change rule."""
    def __setattr__(self, name, value):
        if name != "writes":
            object.__setattr__(self, "writes", getattr(self, "writes", 0) + 1)
        object.__setattr__(self, name, value)


def make_ui(tmp_path):
    ui = gb.GaussianInputBuilderUI.__new__(gb.GaussianInputBuilderUI)
    ui._init_state()
    for attr, _ in gb.STEP_WIDGETS.values():
        setattr(ui, attr, SimpleNamespace(value=None, items=[], enabled=True))
    for attr in gb.CONTROL_ATTRIBUTES:
        setattr(ui, attr, SimpleNamespace(value="", text="", items=[], enabled=True, current_tab=0))
    ui.formchk_switch.value = True
    ui.stem_input.value, ui.output_dir_input.value, ui.geometry_input.value = "water", str(tmp_path), WATER
    ui._load_step(gi.GaussianStepSpec(multiplicity=1))
    ui.refresh()
    ui.dialogs, ui.answers = [], []

    async def record(title, message):
        ui.dialogs.append((title, message))

    async def confirm(title, message):
        ui.dialogs.append((title, message))
        return ui.answers.pop(0) if ui.answers else True
    ui._show_error = ui._show_info = record
    ui._confirm = confirm
    return ui


def change(ui, attr, value):
    """Set a widget as the user would, then fire its on_change."""
    widget = getattr(ui, attr)
    widget.value = value
    ui.on_field_change(widget)


def non_default_step() -> gi.GaussianStepSpec:
    default = gi.GaussianStepSpec()
    values = {}
    for f in fields(default):
        _, kind = gb.STEP_WIDGETS[f.name]
        value = getattr(default, f.name)
        if kind == gb.BOOL:
            values[f.name] = not value
        elif kind == gb.INT:
            values[f.name] = value + 2
        elif kind == gb.OPTIONAL_INT:
            values[f.name] = 7
        elif kind == gb.TEXT:
            values[f.name] = value + "x"
        else:
            values[f.name] = next(c for c in gi.choices(f.name, default) if c != value)
    values["family"], values["method_name"] = gk.FAMILY_SEMI, "PM6"
    return gi.GaussianStepSpec(**values)


class TestFieldTable:
    def test_every_step_field_has_a_widget(self):
        assert set(gb.STEP_WIDGETS) == {f.name for f in fields(gi.GaussianStepSpec)}

    def test_widget_attributes_are_unique(self):
        attrs = [attr for attr, _ in gb.STEP_WIDGETS.values()] + list(gb.CONTROL_ATTRIBUTES)
        assert len(attrs) == len(set(attrs))

    def test_every_field_round_trips_through_the_widgets(self, tmp_path):
        ui = make_ui(tmp_path)
        step = non_default_step()
        default = gi.GaussianStepSpec()
        assert all(getattr(step, f.name) != getattr(default, f.name) for f in fields(step))
        ui._load_step(step)
        assert ui._collect_step() == step

    def test_blank_numbers(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.charge_input.value = ui.multiplicity_input.value = ui.opt_max_cycles_input.value = ""
        step = ui._collect_step()
        assert (step.charge, step.multiplicity, step.opt_max_cycles) == (0, None, None)

    @pytest.mark.parametrize("attr,value,label", [("charge_input", "0.5", "Charge"),
                                                  ("nproc_input", "x", "%NProcShared"),
                                                  ("fukui_plus_mult_input", "2.5", "Multiplicity of N\\+p")])
    def test_bad_numbers_name_the_field(self, tmp_path, attr, value, label):
        ui = make_ui(tmp_path)
        getattr(ui, attr).value = value
        with pytest.raises(ValueError, match=label):
            ui._collect_step()

    def test_assign_writes_only_on_change(self):
        widget = Recorder(value="a")
        widget.writes = 0
        gb._assign(widget, "value", "a")
        gb._assign(widget, "value", "b")
        assert widget.writes == 1 and widget.value == "b"


class TestRefresh:
    def test_route_and_summary(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT_FREQ)
        assert ui.route_label.text == "#P B3LYP/6-31G(d) Opt Freq"
        assert ui.summary_label.text == "3 atoms, charge 0, 10 electrons, multiplicity 1"

    def test_a_bad_number_shows_in_the_route_line(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "charge_input", "x")
        assert ui.route_label.text.startswith("Cannot build the route yet:")

    def test_family_change_offers_its_methods(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "family_selection", gk.FAMILY_SEMI)
        assert ui.method_selection.items == list(gi.method_choices(gk.FAMILY_SEMI))
        assert ui.method_selection.value == "PM7" and ui.route_label.text == "#P PM7"

    def test_choosing_fukui_resets_what_the_fukui_job_sets(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.guess_selection.value = "Read checkpoint"
        ui.title_input.value = "t"
        change(ui, "job_type_selection", gk.JOB_FUKUI)
        assert (ui.guess_selection.value, ui.title_input.value) == ("Default", "")
        assert ui.scf_tight_switch.value and ui.nosymm_switch.value
        assert ui.route_label.text.startswith("job 1:     #P B3LYP/6-31G(d) SCF=Tight NoSymm Opt")

    def test_changes_while_syncing_are_ignored(self, tmp_path):
        ui = make_ui(tmp_path)
        ui._syncing = True
        change(ui, "job_type_selection", gk.JOB_OPT)
        assert ui.route_label.text == "#P B3LYP/6-31G(d)"

    def test_step_list_follows_the_steps(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT)
        assert ui.step_selection.items == ["1: Optimization B3LYP/6-31G(d)"]
        assert ui.step_selection.value == "1: Optimization B3LYP/6-31G(d)"


class TestControls:
    @pytest.mark.parametrize("changes,enabled,disabled", [
        (dict(), ["job_type_selection", "basis_selection", "diffuse_selection"],
         ["opt_target_selection", "irc_direction_selection", "fukui_p_input", "n_states_input",
          "method_other_input", "basis_other_input", "solvent_selection", "remove_step_button"]),
        (dict(job_type=gk.JOB_OPT), ["opt_target_selection", "modredundant_input", "modred_add_button"],
         ["freq_raman_switch"]),
        (dict(job_type=gk.JOB_FUKUI), ["fukui_p_input", "force_constants_selection"],
         ["guess_selection", "title_input", "chk_input", "state_selection", "spin_selection",
          "additional_input_input", "add_step_button", "opt_target_selection"]),
        (dict(state=gk.STATE_TD), ["n_states_input", "coefficients_switch"], []),
        (dict(basis="cc-pVDZ"), [], ["diffuse_selection", "polarization_selection"]),
        (dict(family=gk.FAMILY_SEMI, method_name="PM6"), [], ["basis_selection"]),
        (dict(method_name=gk.OTHER, solvation="SMD", solvent=gk.OTHER),
         ["method_other_input", "solvent_selection", "solvent_other_input"], []),
        (dict(geometry_source=gk.GEOM_ALLCHECK), [], ["connectivity_switch", "charge_input"]),
    ])
    def test_enabled_controls(self, tmp_path, changes, enabled, disabled):
        ui = make_ui(tmp_path)
        ui.refresh_controls(replace(gi.GaussianStepSpec(), **changes))
        assert [a for a in enabled if not getattr(ui, a).enabled] == []
        assert [a for a in disabled if getattr(ui, a).enabled] == []

    def test_fukui_note(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.refresh_controls(gi.GaussianStepSpec(job_type=gk.JOB_FUKUI))
        assert "three linked jobs" in ui.fukui_note_label.text
        ui.refresh_controls(gi.GaussianStepSpec())
        assert ui.fukui_note_label.text == ""


def test_window_module_imports_toga_but_the_core_does_not():
    import subprocess, sys
    code = "import sys, gaussianInputBuilder; print('toga' in sys.modules, 'cubeTab' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=str(SRC))
    assert out.stdout.strip() == "True False", out.stderr


FUKUI_INPUTS = SRC / "tests" / "data" / "fukui" / "inputs"


def load_phenol(ui, ions=True):
    ui.apply_geometry_file(FUKUI_INPUTS / "phenol_N.gjf")
    change(ui, "job_type_selection", gk.JOB_FUKUI)
    if ions:
        ui.fukui_plus_mult_input.value = ui.fukui_minus_mult_input.value = "2"
        ui.refresh()


class TestSteps:
    def test_add_linked_step_continues_from_the_checkpoint(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT)
        ui.add_linked_step(None)
        assert len(ui.steps) == 2 and ui.current_step == 1
        assert ui.geometry_source_selection.value == gk.GEOM_ALLCHECK and ui.guess_selection.value == "Read checkpoint"
        assert ui.step_selection.items == ["1: Optimization B3LYP/6-31G(d)",
                                           "2: Energy B3LYP/6-31G(d), Geom=AllCheck"]
        assert ui.route_label.text == "#P B3LYP/6-31G(d) Geom=AllCheck Guess=Read"

    def test_switching_steps_stores_and_loads_the_panels(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.add_linked_step(None)
        change(ui, "state_selection", gk.STATE_TD)
        ui.step_selection.value = ui._step_labels[0]
        ui.on_step_select(ui.step_selection)
        assert ui.current_step == 0 and ui.state_selection.value == gk.STATE_GROUND
        assert ui.steps[1].state == gk.STATE_TD

    def test_a_step_with_a_bad_number_is_not_left(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.add_linked_step(None)
        ui.n_states_input.value = "x"
        ui.step_selection.value = ui._step_labels[0]
        ui.on_step_select(ui.step_selection)
        assert ui.current_step == 1 and "Fix step 2" in ui.status_label.text
        assert ui.step_selection.value == ui._step_labels[1]

    def test_remove_and_move(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.remove_step(None)
        assert len(ui.steps) == 1 and "only step" in ui.status_label.text
        ui.add_linked_step(None)
        ui.move_step_up(None)
        assert ui.current_step == 0 and ui.steps[0].geometry_source == gk.GEOM_ALLCHECK
        ui.move_step_down(None)
        assert ui.current_step == 1 and ui.steps[1].geometry_source == gk.GEOM_ALLCHECK
        ui.remove_step(None)
        assert len(ui.steps) == 1 and ui.current_step == 0

    def test_a_fukui_step_takes_no_linked_step(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_FUKUI)
        ui.add_linked_step(None)
        assert len(ui.steps) == 1 and "three linked jobs" in ui.status_label.text

    def test_helper_row_appends_a_line(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.modredundant_input.value = "B 1 2 F"
        ui.modred_type_selection.value, ui.modred_atoms_input.value = "Angle", "2 1 3"
        ui.modred_action_selection.value, ui.modred_steps_input.value, ui.modred_size_input.value = "Scan", "5", "2.0"
        ui.add_modredundant_line(None)
        assert ui.modredundant_input.value == "B 1 2 F\nA 2 1 3 S 5 2.0"
        ui.modred_atoms_input.value = "1"
        ui.add_modredundant_line(None)
        assert "needs 3 atom numbers" in ui.status_label.text


class TestGeometryAndFiles:
    def test_apply_fills_the_form_but_not_the_ion_spins(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.stem_input.value = fi.DEFAULT_STEM
        ui.output_dir_input.value = ""
        ui.fukui_plus_mult_input.value = ui.fukui_minus_mult_input.value = "4"
        ui.apply_geometry_file(FUKUI_INPUTS / "phenol_Np1.gjf")
        assert ui.charge_input.value == "-1" and ui.multiplicity_input.value == "2"
        assert ui.fukui_plus_mult_input.value == "" and ui.fukui_minus_mult_input.value == ""
        assert ui.stem_input.value == "phenol_Np1" and ui.output_dir_input.value == str(FUKUI_INPUTS)

    def test_apply_from_an_n_fchk_drops_the_state_suffix(self, tmp_path):
        from fukui_fixtures import N
        ui = make_ui(tmp_path)
        ui.stem_input.value = ""
        ui.apply_geometry_file(N)
        assert ui.stem_input.value == "phenol" and ui.output_dir_input.value == str(tmp_path)

    def test_apply_keeps_a_typed_job_name(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.stem_input.value = "mine"
        ui.apply_geometry_file(FUKUI_INPUTS / "phenol_N.gjf")
        assert ui.stem_input.value == "mine"

    def test_open_input_restores_every_step(self, tmp_path):
        ui = make_ui(tmp_path)
        change(ui, "job_type_selection", gk.JOB_OPT)
        ui.add_linked_step(None)
        change(ui, "state_selection", gk.STATE_TD)
        asyncio.run(ui.save_input(None))
        fresh = make_ui(tmp_path)
        notes = fresh.apply_input_file(tmp_path / "water.gjf")
        assert notes == [] and len(fresh.steps) == 2 and fresh.current_step == 0
        assert fresh.job_type_selection.value == gk.JOB_OPT and fresh.steps[1].state == gk.STATE_TD
        assert fresh.geometry_input.value == "\n".join(fi.format_geometry(fi.parse_geometry(WATER)))

    def test_open_input_errors_go_to_a_dialog(self, tmp_path):
        ui = make_ui(tmp_path)
        bad = tmp_path / "bad.gjf"
        bad.write_text("nothing\n")
        ui._ask_open_path = lambda: _resolved(bad)
        asyncio.run(ui.open_input(None))
        assert ui.dialogs and ui.dialogs[0][0] == "Open Input Error" and "bad.gjf" in ui.dialogs[0][1]


async def _resolved(value):
    return value


class TestActions:
    def test_preview_fills_the_input_text(self, tmp_path):
        ui = make_ui(tmp_path)
        asyncio.run(ui.preview(None))
        assert ui.preview_output.value.startswith("%Chk=water.chk") and not ui.dialogs
        assert ui.tabs.current_tab == gb.PREVIEW_TAB

    def test_preview_errors_go_to_a_dialog(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.multiplicity_input.value = ""
        asyncio.run(ui.preview(None))
        assert ui.dialogs[0][0] == "Input Error" and "it is not assumed" in ui.dialogs[0][1]

    def test_save_writes_the_files_and_says_how_to_run_them(self, tmp_path):
        ui = make_ui(tmp_path)
        asyncio.run(ui.save_input(None))
        assert (tmp_path / "water.gjf").is_file() and (tmp_path / "water_formchk.sh").is_file()
        assert "g16 < water.gjf > water.log" in ui.messages_output.value
        assert "bash water_formchk.sh" in ui.messages_output.value

    def test_existing_files_are_kept_unless_confirmed(self, tmp_path):
        ui = make_ui(tmp_path)
        (tmp_path / "water.gjf").write_text("mine")
        ui.answers = [False]
        asyncio.run(ui.save_input(None))
        assert (tmp_path / "water.gjf").read_text() == "mine"
        ui.answers = [True]
        asyncio.run(ui.save_input(None))
        assert "%Chk=water.chk" in (tmp_path / "water.gjf").read_text()

    def test_defaults_ask_first_and_keep_the_molecule(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.add_linked_step(None)
        ui.answers = [False]
        asyncio.run(ui.defaults(None))
        assert len(ui.steps) == 2
        ui.answers = [True]
        asyncio.run(ui.defaults(None))
        assert ui.steps == [gi.GaussianStepSpec(multiplicity=1)] and ui.geometry_input.value == WATER
        ui.dialogs.clear()
        asyncio.run(ui.defaults(None))
        assert ui.dialogs == []          # nothing changed: no question


class TestFukuiInTheWindow:
    """Migrated from tests/test_fukuiInputBuilder.py: the Fukui builder now lives here as a Job Type."""

    def test_every_fukui_spec_field_is_read_from_a_widget(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        ui.stem_input.value, ui.memory_input.value, ui.nproc_input.value = "ph", "3GB", "3"
        ui.method_selection.value, ui.basis_selection.value = "M06-2X", "Def2SVP"
        ui.additional_keywords_input.value = "SCF=XQC"
        ui.opt_tight_switch.value = True
        ui.print_level_selection.value = "#N"
        ui.charge_input.value, ui.multiplicity_input.value = "1", "2"
        ui.fukui_p_input.value, ui.fukui_q_input.value = "2", "3"
        ui.fukui_plus_mult_input.value, ui.fukui_minus_mult_input.value = "3", "5"
        ui.fukui_guess_read_switch.value = False
        spec = gi.fukui_spec(ui._collect_inputs())
        defaults = fi.FukuiJobSpec("", "", "", 0, None, None, None)
        for f in fields(fi.FukuiJobSpec):
            assert getattr(spec, f.name) != getattr(defaults, f.name), f"{f.name} is not read from the window"

    def test_blank_multiplicities_are_not_filled_in(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui, ions=False)
        spec = gi.fukui_spec(ui._collect_inputs())
        assert (spec.multiplicity, spec.plus_multiplicity, spec.minus_multiplicity) == (1, None, None)

    def test_summary_and_route_list_the_three_states(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.stem_input.value = fi.DEFAULT_STEM
        load_phenol(ui)
        text = ui.summary_label.text
        assert "N+1" in text and "charge -1" in text and "51 electrons" in text and "phenol_Np1.chk" in text
        assert ui.route_label.text == ("job 1:     #P B3LYP/6-31G(d) SCF=Tight NoSymm Opt\n"
                                       "jobs 2-3:  #P B3LYP/6-31G(d) SCF=Tight NoSymm Geom=Check Guess=Read")

    def test_preview_asks_for_missing_spins(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui, ions=False)
        asyncio.run(ui.preview(None))
        assert ui.dialogs[0][0] == "Input Error"
        assert "N+1 state: enter its spin multiplicity; it is not assumed" in ui.dialogs[0][1]

    def test_preview_shows_three_linked_jobs(self, tmp_path):
        ui = make_ui(tmp_path)
        load_phenol(ui)
        asyncio.run(ui.preview(None))
        assert ui.preview_output.value.count("--Link1--") == 2 and not ui.dialogs

    def test_save_writes_both_files_and_says_how_to_run_them(self, tmp_path):
        ui = make_ui(tmp_path)
        ui.stem_input.value = fi.DEFAULT_STEM
        load_phenol(ui)
        ui.output_dir_input.value = str(tmp_path)
        asyncio.run(ui.save_input(None))
        assert (tmp_path / "phenol_fukui.gjf").is_file() and (tmp_path / "phenol_formchk.sh").is_file()
        assert "g16 < phenol_fukui.gjf > phenol_fukui.log" in ui.messages_output.value
        assert "bash phenol_formchk.sh" in ui.messages_output.value
        assert "Analysis of Fukui Functions" in ui.messages_output.value
