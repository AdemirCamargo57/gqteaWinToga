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
