"""
Inputs > Gaussian 16 Input Builder: the window over ``gaussianInput``.

Modelled on GaussView 6.1's Gaussian Calculation Setup: a live route line, one tab per panel, an
Additional Keywords field, and a step bar for linked jobs (--Link1--). The window only moves values
between widgets and ``GaussianStepSpec``s (through the STEP_WIDGETS table) and turns every
ValueError into a dialog; validation, the input text, reopening and the Fukui job live in
``gaussianInput``. Handlers are synchronous and ``_syncing``-guarded, and widgets are written only
when a value changes (every write re-lays out the window on WinForms).
"""

import dataclasses
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import gaussianInput as gi
import gaussianKeywords as gk
from help import HelpGqteaWin
from systemProperties import parse_integer

CHOICE, TEXT, BOOL, INT, OPTIONAL_INT = "choice", "text", "bool", "int", "optional_int"

# GaussianStepSpec field -> (widget attribute, kind). _collect_step and _load_step go through this
# table, and a test checks it names every field, so a new spec field without a control fails the suite.
STEP_WIDGETS: Dict[str, Tuple[str, str]] = {
    "job_type": ("job_type_selection", CHOICE),
    "opt_target": ("opt_target_selection", CHOICE),
    "force_constants": ("force_constants_selection", CHOICE),
    "opt_tight": ("opt_tight_switch", BOOL),
    "opt_max_cycles": ("opt_max_cycles_input", OPTIONAL_INT),
    "freq_raman": ("freq_raman_switch", BOOL),
    "freq_vcd": ("freq_vcd_switch", BOOL),
    "freq_anharmonic": ("freq_anharmonic_switch", BOOL),
    "irc_direction": ("irc_direction_selection", CHOICE),
    "irc_force_constants": ("irc_force_constants_selection", CHOICE),
    "irc_max_points": ("irc_max_points_input", OPTIONAL_INT),
    "irc_step_size": ("irc_step_size_input", OPTIONAL_INT),
    "stable_opt": ("stable_opt_switch", BOOL),
    "nmr_method": ("nmr_method_selection", CHOICE),
    "nmr_spinspin": ("nmr_spinspin_switch", BOOL),
    "fukui_p": ("fukui_p_input", INT),
    "fukui_q": ("fukui_q_input", INT),
    "fukui_plus_multiplicity": ("fukui_plus_mult_input", OPTIONAL_INT),
    "fukui_minus_multiplicity": ("fukui_minus_mult_input", OPTIONAL_INT),
    "fukui_guess_read": ("fukui_guess_read_switch", BOOL),
    "state": ("state_selection", CHOICE),
    "n_states": ("n_states_input", INT),
    "root": ("root_input", INT),
    "state_spin": ("state_spin_selection", CHOICE),
    "print_all_coefficients": ("coefficients_switch", BOOL),
    "spin": ("spin_selection", CHOICE),
    "family": ("family_selection", CHOICE),
    "method_name": ("method_selection", CHOICE),
    "method_other": ("method_other_input", TEXT),
    "basis": ("basis_selection", CHOICE),
    "basis_diffuse": ("diffuse_selection", CHOICE),
    "basis_polarization": ("polarization_selection", CHOICE),
    "basis_other": ("basis_other_input", TEXT),
    "dispersion": ("dispersion_selection", CHOICE),
    "geometry_source": ("geometry_source_selection", CHOICE),
    "charge": ("charge_input", INT),
    "multiplicity": ("multiplicity_input", OPTIONAL_INT),
    "modredundant": ("modredundant_input", TEXT),
    "title": ("title_input", TEXT),
    "memory": ("memory_input", TEXT),
    "nproc": ("nproc_input", INT),
    "write_chk": ("write_chk_switch", BOOL),
    "chk": ("chk_input", TEXT),
    "old_chk": ("old_chk_input", TEXT),
    "extra_link0": ("extra_link0_input", TEXT),
    "print_level": ("print_level_selection", CHOICE),
    "scf_tight": ("scf_tight_switch", BOOL),
    "scf_qc": ("scf_qc_switch", BOOL),
    "nosymm": ("nosymm_switch", BOOL),
    "polar": ("polar_switch", BOOL),
    "connectivity": ("connectivity_switch", BOOL),
    "guess": ("guess_selection", CHOICE),
    "guess_mix": ("guess_mix_switch", BOOL),
    "guess_always": ("guess_always_switch", BOOL),
    "guess_only": ("guess_only_switch", BOOL),
    "guess_save": ("guess_save_switch", BOOL),
    "pop_level": ("pop_level_selection", CHOICE),
    "pop_mk": ("pop_mk_switch", BOOL),
    "pop_hirshfeld": ("pop_hirshfeld_switch", BOOL),
    "pop_nbo": ("pop_nbo_switch", BOOL),
    "pop_npa": ("pop_npa_switch", BOOL),
    "pop_savento": ("pop_savento_switch", BOOL),
    "solvation": ("solvation_selection", CHOICE),
    "solvent": ("solvent_selection", CHOICE),
    "solvent_other": ("solvent_other_input", TEXT),
    "additional_input": ("additional_input_input", TEXT),
    "additional_keywords": ("additional_keywords_input", TEXT),
}

# Every other widget the handlers touch (job-wide fields, labels, step bar, helper row, outputs).
CONTROL_ATTRIBUTES = (
    "stem_input", "output_dir_input", "geometry_input", "formchk_switch", "route_label", "summary_label",
    "status_label", "fukui_note_label", "step_selection", "add_step_button", "remove_step_button",
    "step_up_button", "step_down_button", "load_geometry_button", "modred_type_selection", "modred_atoms_input",
    "modred_action_selection", "modred_steps_input", "modred_size_input", "modred_add_button",
    "preview_output", "messages_output", "tabs",
)

TAB_TITLES = ("Job Type", "Method", "Molecule", "Title", "Link 0", "General", "Guess", "Pop", "Solvation",
              "Add. Inp.", "Preview", "Help")
PREVIEW_TAB = TAB_TITLES.index("Preview")
INPUT_FILE_TYPES = [suffix.lstrip(".") for suffix in gi.GAUSSIAN_INPUT_SUFFIXES]
FUKUI_NOTE = ("The Fukui job writes three linked jobs itself: the title, %Chk/%OldChk, the guess, Add. Inp., "
              "connectivity, the spin prefix and the electronic state are set by it.")


def _assign(widget, attr: str, value) -> None:
    """Write a widget property only when it changes (on WinForms every write re-lays out the window)."""
    if getattr(widget, attr, None) != value:
        setattr(widget, attr, value)


def _item_names(selection) -> List[str]:
    """The strings a Toga Selection holds (its rows carry them as ``.value``)."""
    return [getattr(item, "value", item) for item in selection.items]


class GaussianInputBuilderUI:
    LABEL_WIDTH = 190
    HINT_COLOR = "#666666"

    def __init__(self, *args, **kwargs):
        # gqteaWinToga passes the Button that was pressed; ignore it.
        self._init_state()
        self.window = toga.Window(title="Gaussian 16 Input Builder", size=(940, 720))
        self.layout_main_window()
        self.window.show()

    def _init_state(self):
        """Everything but widgets (shared with the headless tests)."""
        self.builder = gi.GaussianInputBuilder()
        self.steps: List[gi.GaussianStepSpec] = [gi.GaussianStepSpec()]
        self.current_step = 0
        self.job: Optional[gi.GaussianJob] = None
        self._syncing = False
        self._step_labels: List[str] = []

    # --- widget factories -------------------------------------------------
    def _label(self, text, width=None):
        return toga.Label(text, style=Pack(width=width or self.LABEL_WIDTH, margin=(6, 6, 0, 0)))

    def _row(self, *children):
        return toga.Box(children=list(children), style=Pack(direction=ROW, margin_bottom=4))

    def _heading(self, text):
        return toga.Label(text, style=Pack(font_weight="bold", margin=(10, 0, 4, 0)))

    def _hint(self, text):
        return toga.Label(text, style=Pack(color=self.HINT_COLOR, font_size=9, margin=(0, 0, 4, 0)))

    def _panel(self, *children):
        return toga.ScrollContainer(content=toga.Box(children=list(children), style=Pack(direction=COLUMN, margin=10)),
                                    style=Pack(flex=1))

    def _field(self, label, widget, hint=None):
        return self._row(self._label(label), widget, *([self._hint("   " + hint)] if hint else []))

    def _text(self, attr, width=None, placeholder=""):
        widget = toga.TextInput(placeholder=placeholder, style=Pack(width=width) if width else Pack(flex=1),
                                on_change=self.on_field_change)
        setattr(self, attr, widget)
        return widget

    def _multiline(self, attr, height, placeholder=""):
        widget = toga.MultilineTextInput(placeholder=placeholder, on_change=self.on_field_change,
                                         style=Pack(flex=1, height=height, font_family="monospace"))
        setattr(self, attr, widget)
        return widget

    def _choice(self, attr, items: Sequence[str], width=240):
        widget = toga.Selection(items=list(items), style=Pack(width=width), on_change=self.on_field_change)
        setattr(self, attr, widget)
        return widget

    def _switch(self, attr, text):
        widget = toga.Switch(text, on_change=self.on_field_change, style=Pack(margin_right=12))
        setattr(self, attr, widget)
        return widget

    def _checkbox_row(self, *switches):
        return self._row(self._label(""), *switches)

    # --- layout -------------------------------------------------------------
    def layout_main_window(self):
        self._syncing = True          # widgets fire on_change while they are built
        self.route_label = toga.Label("", style=Pack(font_family="monospace", margin=(4, 0, 4, 0), flex=1))
        self.summary_label = toga.Label("", style=Pack(font_family="monospace", margin=(4, 0, 4, 0)))
        self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
        self.fukui_note_label = toga.Label("", style=Pack(color=self.HINT_COLOR, margin=(6, 0, 4, 0)))

        self.step_selection = toga.Selection(items=["1"], style=Pack(width=520), on_change=self.on_step_select)
        self.add_step_button = toga.Button("Add linked step", on_press=self.add_linked_step, style=Pack(margin_left=6))
        self.remove_step_button = toga.Button("Remove", on_press=self.remove_step, style=Pack(margin_left=6))
        self.step_up_button = toga.Button("Up", on_press=self.move_step_up, style=Pack(margin_left=6, width=50))
        self.step_down_button = toga.Button("Down", on_press=self.move_step_down, style=Pack(margin_left=6, width=60))
        header = toga.Box(children=[
            self._row(toga.Label("Step:", style=Pack(margin=(6, 6, 0, 0))), self.step_selection,
                      self.add_step_button, self.remove_step_button, self.step_up_button, self.step_down_button),
            self._row(toga.Label("Route:", style=Pack(margin=(4, 6, 0, 0))), self.route_label),
        ], style=Pack(direction=COLUMN))

        self.preview_output = toga.MultilineTextInput(readonly=True, style=Pack(flex=2, font_family="monospace"))
        self.messages_output = toga.MultilineTextInput(readonly=True, style=Pack(flex=1, font_family="monospace"))
        preview = toga.Box(children=[self.preview_output, self.messages_output],
                           style=Pack(direction=COLUMN, margin=10, flex=1))
        help_box = toga.MultilineTextInput(readonly=True, style=Pack(flex=1),
                                           value=getattr(HelpGqteaWin, "help_gaussian_input", ""))
        panels = (self._build_job_panel(), self._build_method_panel(), self._build_molecule_panel(),
                  self._build_title_panel(), self._build_link0_panel(), self._build_general_panel(),
                  self._build_guess_panel(), self._build_pop_panel(), self._build_solvation_panel(),
                  self._build_additional_input_panel(), preview, help_box)
        self.tabs = toga.OptionContainer(content=list(zip(TAB_TITLES, panels)), style=Pack(flex=1))

        self.stem_input = toga.TextInput(value=gi.fi.DEFAULT_STEM, style=Pack(width=220), on_change=self.on_field_change)
        self.output_dir_input = toga.TextInput(placeholder="Required: folder for the .gjf and .sh files",
                                               style=Pack(flex=1), on_change=self.on_field_change)
        footer = toga.Box(children=[
            self._row(self._label("Additional Keywords:"), self._text("additional_keywords_input",
                                                                      placeholder="e.g. Int=UltraFine 6D")),
            self._hint("Any keyword or option the panels do not have. Options of a keyword a panel writes are "
                       "merged into it (Opt=NoEigenTest)."),
            self._row(self._label("Job name (file stem):"), self.stem_input, self._label("   Output folder:", 120),
                      self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
            self._row(*(toga.Button(text, on_press=handler, style=Pack(width=120, margin=5)) for text, handler in (
                ("Open input...", self.open_input), ("Preview", self.preview), ("Save input", self.save_input),
                ("Defaults", self.defaults), ("Close", self.closeTopLevel)))),
            self.status_label,
        ], style=Pack(direction=COLUMN))
        self.window.content = toga.Box(children=[header, self.tabs, footer],
                                       style=Pack(direction=COLUMN, margin=10, flex=1))
        self._syncing = False
        self._load_step(self.steps[0])
        self.refresh()

    def _build_job_panel(self):
        return self._panel(
            self._field("Job type:", self._choice("job_type_selection", gk.JOB_TYPES, 280)),
            self._heading("Optimization"),
            self._field("Optimize to a:", self._choice("opt_target_selection", gk.OPT_TARGETS)),
            self._field("Force constants:", self._choice("force_constants_selection", gk.FORCE_CONSTANTS, 280)),
            self._checkbox_row(self._switch("opt_tight_switch", "Tight convergence")),
            self._field("Maximum cycles:", self._text("opt_max_cycles_input", 100, "default")),
            self._hint("Scan and freeze lines are entered in the Molecule panel."),
            self._heading("Frequency"),
            self._checkbox_row(self._switch("freq_raman_switch", "Raman"), self._switch("freq_vcd_switch", "VCD"),
                               self._switch("freq_anharmonic_switch", "Anharmonic")),
            self._heading("IRC"),
            self._field("Direction:", self._choice("irc_direction_selection", gk.IRC_DIRECTIONS)),
            self._field("Force constants:", self._choice("irc_force_constants_selection", gk.FORCE_CONSTANTS, 280)),
            self._field("MaxPoints:", self._text("irc_max_points_input", 100, "default")),
            self._field("StepSize:", self._text("irc_step_size_input", 100, "default"), "units of 0.01 bohr"),
            self._heading("Stability"),
            self._checkbox_row(self._switch("stable_opt_switch", "Reoptimize the wavefunction if unstable (Stable=Opt)")),
            self._heading("NMR"),
            self._field("Method:", self._choice("nmr_method_selection", gk.NMR_METHODS)),
            self._checkbox_row(self._switch("nmr_spinspin_switch", "Spin-spin couplings")),
            self._heading("Fukui functions (N, N+p, N-q)"),
            self._field("p (electrons added):", self._text("fukui_p_input", 80)),
            self._field("q (electrons removed):", self._text("fukui_q_input", 80)),
            self._field("Multiplicity of N+p:", self._text("fukui_plus_mult_input", 80, "required")),
            self._field("Multiplicity of N-q:", self._text("fukui_minus_mult_input", 80, "required")),
            self._checkbox_row(self._switch("fukui_guess_read_switch",
                                            "Start the ions from the optimized N-state orbitals (Guess=Read)")),
            self._hint("Job 1 optimizes the N state (Optimization options above, minimum only); jobs 2-3 are single "
                       "points with Geom=Check. Writes <name>_fukui.gjf and <name>_formchk.sh."),
            self.fukui_note_label)

    def _build_method_panel(self):
        return self._panel(
            self._field("Electronic state:", self._choice("state_selection", gk.STATES)),
            self._field("Number of states:", self._text("n_states_input", 80)),
            self._field("State of interest (Root):", self._text("root_input", 80)),
            self._field("States:", self._choice("state_spin_selection", gk.STATE_SPINS)),
            self._checkbox_row(self._switch("coefficients_switch", "Print all excitation coefficients for the NTO "
                                                                   "and Excitation tabs (IOp(9/40=4))")),
            self._heading("Model chemistry"),
            self._field("Spin:", self._choice("spin_selection", gk.SPIN_PREFIXES)),
            self._field("Method family:", self._choice("family_selection", gk.FAMILIES)),
            self._field("Method:", self._choice("method_selection", gi.method_choices(gk.FAMILY_DFT))),
            self._field("Other method:", self._text("method_other_input", 200), "when Method is Other..."),
            self._field("Basis set:", self._choice("basis_selection", gk.BASES)),
            self._row(self._label("Diffuse, polarization:"), self._choice("diffuse_selection", gk.DIFFUSE, 90),
                      self._choice("polarization_selection", gk.POLARIZATION, 130),
                      self._hint("   6-31G and 6-311G only")),
            self._field("Other basis set:", self._text("basis_other_input", 200), "when Basis set is Other..."),
            self._field("Empirical dispersion:", self._choice("dispersion_selection", gk.DISPERSION)))

    def _build_molecule_panel(self):
        self.load_geometry_button = toga.Button("Load geometry...", on_press=self.load_geometry, style=Pack(width=150))
        self.geometry_input = toga.MultilineTextInput(
            placeholder="Cartesian geometry, Angstrom:\nC   0.000000   0.000000   0.000000\n...",
            style=Pack(flex=1, height=170, font_family="monospace"), on_change=self.on_field_change)
        self.modred_type_selection = toga.Selection(items=list(gk.MODRED_TYPES), style=Pack(width=100))
        self.modred_atoms_input = toga.TextInput(placeholder="atoms, e.g. 1 2 3 4", style=Pack(width=150))
        self.modred_action_selection = toga.Selection(items=list(gk.MODRED_ACTIONS), style=Pack(width=90))
        self.modred_steps_input = toga.TextInput(placeholder="steps", style=Pack(width=60))
        self.modred_size_input = toga.TextInput(placeholder="step size", style=Pack(width=80))
        self.modred_add_button = toga.Button("Add", on_press=self.add_modredundant_line, style=Pack(width=60))
        return self._panel(
            self._row(self.load_geometry_button,
                      self._hint("   .xyz .gjf .com .fchk .mol .sdf .mol2 .pdb; a .gjf/.fchk also sets the charge "
                                 "and multiplicity")),
            self.geometry_input,
            self._hint("The geometry belongs to the whole job: step 1 writes it, and a later step can write it again."),
            self._field("Geometry source:", self._choice("geometry_source_selection", gk.GEOMETRY_SOURCES, 380)),
            self._field("Charge:", self._text("charge_input", 80)),
            self._field("Multiplicity:", self._text("multiplicity_input", 80, "required"),
                        "never assumed; the line below lists the allowed values"),
            self.summary_label,
            self._heading("Scan and freeze coordinates (Opt=ModRedundant)"),
            self._row(self.modred_type_selection, self.modred_atoms_input, self.modred_action_selection,
                      self.modred_steps_input, self.modred_size_input, self.modred_add_button),
            self._hint("Atoms numbered from 1. Scan: number of steps and step size (Angstrom or degrees). "
                       "Typed lines: B/A/D, the atoms, then F, or S steps size."),
            self._multiline("modredundant_input", 90, "B 1 2 S 10 0.05\nD 1 2 3 4 F"))

    def _build_title_panel(self):
        return self._panel(
            self._field("Title:", self._text("title_input", placeholder="blank = the job name")),
            self._hint("One line. Gaussian does not allow @ # ! - _ \\ in the title."))

    def _build_link0_panel(self):
        return self._panel(
            self._field("%Mem:", self._text("memory_input", 120), "integer + unit, e.g. 2GB, 500MB"),
            self._field("%NProcShared:", self._text("nproc_input", 120), "cores of the machine that runs Gaussian"),
            self._checkbox_row(self._switch("write_chk_switch", "Write a checkpoint (%Chk)")),
            self._field("%Chk:", self._text("chk_input", 220, "<name>.chk"), "file name, no folders"),
            self._field("%OldChk:", self._text("old_chk_input", 220, "optional"),
                        "an existing checkpoint to read (Guess=Read, Geom=Check, ReadFC)"),
            self._label("Extra Link 0 lines:"),
            self._multiline("extra_link0_input", 70, "%KJob L502"),
            self._hint("Relative file names: the job runs under Linux from the input's folder."))

    def _build_general_panel(self):
        self.formchk_switch = toga.Switch("Write a formchk script (<name>_formchk.sh) for the 3D viewer's tabs",
                                          value=True, on_change=self.on_field_change)
        return self._panel(
            self._field("Print level:", self._choice("print_level_selection", gk.PRINT_LEVELS, 100)),
            self._checkbox_row(self._switch("scf_tight_switch", "SCF=Tight"), self._switch("scf_qc_switch", "SCF=QC"),
                               self._switch("nosymm_switch", "NoSymm"), self._switch("polar_switch", "Polar")),
            self._checkbox_row(self._switch("connectivity_switch", "Write connectivity (Geom=Connectivity)")),
            self._checkbox_row(self.formchk_switch))

    def _build_guess_panel(self):
        return self._panel(
            self._field("Initial guess:", self._choice("guess_selection", gk.GUESS_METHODS)),
            self._checkbox_row(self._switch("guess_mix_switch", "Mix HOMO and LUMO"),
                               self._switch("guess_always_switch", "Always generate in optimizations")),
            self._checkbox_row(self._switch("guess_only_switch", "Only do guess (no SCF)"),
                               self._switch("guess_save_switch", "Save orbitals")))

    def _build_pop_panel(self):
        return self._panel(
            self._field("Population output:", self._choice("pop_level_selection", gk.POP_LEVELS)),
            self._checkbox_row(self._switch("pop_mk_switch", "MK charges"),
                               self._switch("pop_hirshfeld_switch", "Hirshfeld (CM5) charges")),
            self._checkbox_row(self._switch("pop_nbo_switch", "NBO"), self._switch("pop_npa_switch", "NPA"),
                               self._switch("pop_savento_switch", "Save NTOs")))

    def _build_solvation_panel(self):
        return self._panel(
            self._field("Solvation model:", self._choice("solvation_selection", gk.SOLVATION_MODELS)),
            self._field("Solvent:", self._choice("solvent_selection", gk.SOLVENTS)),
            self._field("Other solvent:", self._text("solvent_other_input", 200), "Gaussian's name, one word"))

    def _build_additional_input_panel(self):
        return self._panel(
            self._hint("Input sections written after the geometry (and after any connectivity and scan/freeze "
                       "lines): a Gen basis set, ECPs, ... Separate sections with a blank line."),
            self._multiline("additional_input_input", 320))

    # --- panels <-> GaussianStepSpec -----------------------------------------
    def _set_items(self, selection, items: Sequence[str]) -> None:
        if _item_names(selection) != list(items):
            selection.items = list(items)

    def _collect_step(self) -> gi.GaussianStepSpec:
        """The panels as a GaussianStepSpec; ValueError (naming the field) for a non-integer number."""
        defaults = gi.GaussianStepSpec()
        values = {}
        for name, (attr, kind) in STEP_WIDGETS.items():
            widget = getattr(self, attr)
            if kind == BOOL:
                values[name] = bool(widget.value)
            elif kind == CHOICE:
                values[name] = widget.value if widget.value is not None else getattr(defaults, name)
            elif kind == TEXT:
                values[name] = widget.value or ""
            elif kind == INT:
                values[name] = parse_integer(widget.value, gi.FIELD_LABELS[name], getattr(defaults, name))
            else:
                values[name] = parse_integer(widget.value, gi.FIELD_LABELS[name])
        return gi.GaussianStepSpec(**values)

    def _load_step(self, step: gi.GaussianStepSpec) -> None:
        self._syncing = True
        try:
            self._set_items(self.method_selection, gi.method_choices(step.family))
            for name, (attr, kind) in STEP_WIDGETS.items():
                value = getattr(step, name)
                if kind in (INT, OPTIONAL_INT):
                    value = "" if value is None else str(value)
                _assign(getattr(self, attr), "value", value)
        finally:
            self._syncing = False

    def _collect_inputs(self) -> gi.GaussianJobSpec:
        self.steps[self.current_step] = self._collect_step()
        return gi.GaussianJobSpec(stem=(self.stem_input.value or "").strip(),
                                  output_dir=(self.output_dir_input.value or "").strip(),
                                  geometry=self.geometry_input.value or "",
                                  steps=[dataclasses.replace(s) for s in self.steps],
                                  write_formchk=bool(self.formchk_switch.value))

    # --- live updates ---------------------------------------------------------
    def on_field_change(self, widget):
        """Every control (synchronous): keep the panels consistent, then refresh the live texts."""
        if self._syncing:
            return
        if widget is self.family_selection:
            self._syncing = True
            try:
                choices = gi.method_choices(widget.value)
                self._set_items(self.method_selection, choices)
                _assign(self.method_selection, "value", gk.DEFAULT_METHOD.get(widget.value, choices[0]))
            finally:
                self._syncing = False
        elif widget is self.job_type_selection and widget.value == gk.JOB_FUKUI:
            try:
                step = self._collect_step()
            except ValueError:
                step = dataclasses.replace(self.steps[self.current_step], job_type=gk.JOB_FUKUI)
            self._load_step(gi.fukui_defaults(step))
        self.refresh()

    def refresh(self) -> None:
        try:
            spec = self._collect_inputs()
        except ValueError as exc:
            _assign(self.route_label, "text", f"Cannot build the route yet: {exc}")
            return
        _assign(self.route_label, "text", gi.route_preview(spec, self.current_step))
        _assign(self.summary_label, "text", gi.molecule_summary(spec, self.current_step))
        self._refresh_step_list(spec)
        self.refresh_controls(spec.steps[self.current_step])

    def _refresh_step_list(self, spec: gi.GaussianJobSpec) -> None:
        labels = [gi.step_label(i + 1, s) for i, s in enumerate(spec.steps)]
        if labels == self._step_labels:
            return
        self._syncing = True
        try:
            self._set_items(self.step_selection, labels)
            _assign(self.step_selection, "value", labels[self.current_step])
        finally:
            self._syncing = False
        self._step_labels = labels

    def refresh_controls(self, step: gi.GaussianStepSpec) -> None:
        """Enable only what applies to this step (GaussView greys out the same way)."""
        job, fukui = step.job_type, step.job_type == gk.JOB_FUKUI
        optimizing = job in gk.OPTIMIZING_JOBS
        excited = step.state != gk.STATE_GROUND
        model_state = step.state in gk.MODEL_STATES
        basis_used = not gi.no_basis(step)
        coordinates = step.geometry_source == gk.GEOM_COORDINATES
        solvated = step.solvation != gk.SOLVATION_NONE
        rules = {
            "opt_target_selection": optimizing, "force_constants_selection": optimizing or fukui,
            "opt_tight_switch": optimizing or fukui, "opt_max_cycles_input": optimizing or fukui,
            "freq_raman_switch": job in gk.FREQUENCY_JOBS, "freq_vcd_switch": job in gk.FREQUENCY_JOBS,
            "freq_anharmonic_switch": job in gk.FREQUENCY_JOBS,
            "irc_direction_selection": job == gk.JOB_IRC, "irc_force_constants_selection": job == gk.JOB_IRC,
            "irc_max_points_input": job == gk.JOB_IRC, "irc_step_size_input": job == gk.JOB_IRC,
            "stable_opt_switch": job == gk.JOB_STABILITY,
            "nmr_method_selection": job == gk.JOB_NMR, "nmr_spinspin_switch": job == gk.JOB_NMR,
            "fukui_p_input": fukui, "fukui_q_input": fukui, "fukui_plus_mult_input": fukui,
            "fukui_minus_mult_input": fukui, "fukui_guess_read_switch": fukui,
            "state_selection": not fukui, "n_states_input": excited and not fukui,
            "root_input": excited and not fukui, "state_spin_selection": excited and not fukui,
            "coefficients_switch": excited and not fukui, "spin_selection": not fukui,
            "family_selection": not model_state,
            "method_selection": not model_state and step.family in gk.FAMILY_METHODS,
            "method_other_input": not model_state and step.family == gk.FAMILY_DFT and step.method_name == gk.OTHER,
            "basis_selection": basis_used,
            "diffuse_selection": basis_used and step.basis in gk.EXTENSIBLE_POPLE,
            "polarization_selection": basis_used and step.basis in gk.EXTENSIBLE_POPLE,
            "basis_other_input": basis_used and step.basis == gk.OTHER,
            "geometry_source_selection": not fukui,
            "charge_input": step.geometry_source != gk.GEOM_ALLCHECK,
            "multiplicity_input": step.geometry_source != gk.GEOM_ALLCHECK,
            "geometry_input": self.current_step == 0 or coordinates,
            "load_geometry_button": self.current_step == 0 or coordinates,
            "modredundant_input": optimizing, "modred_add_button": optimizing,
            "modred_type_selection": optimizing, "modred_atoms_input": optimizing,
            "modred_action_selection": optimizing, "modred_steps_input": optimizing,
            "modred_size_input": optimizing,
            "title_input": not fukui, "write_chk_switch": not fukui, "chk_input": not fukui and step.write_chk,
            "old_chk_input": not fukui, "extra_link0_input": not fukui,
            "connectivity_switch": not fukui and coordinates,
            "guess_selection": not fukui, "guess_mix_switch": not fukui, "guess_always_switch": not fukui,
            "guess_only_switch": not fukui, "guess_save_switch": not fukui,
            "solvent_selection": solvated, "solvent_other_input": solvated and step.solvent == gk.OTHER,
            "additional_input_input": not fukui,
            "add_step_button": not fukui, "remove_step_button": len(self.steps) > 1,
            "step_up_button": self.current_step > 0, "step_down_button": self.current_step < len(self.steps) - 1,
        }
        for attr, enabled in rules.items():
            _assign(getattr(self, attr), "enabled", enabled)
        _assign(self.fukui_note_label, "text", FUKUI_NOTE if fukui else "")

    # --- dialogs ------------------------------------------------------------
    async def _show_error(self, title: str, message: str):
        await self.window.dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        await self.window.dialog(toga.InfoDialog(title, message))

    async def _confirm(self, title: str, message: str) -> bool:
        return bool(await self.window.dialog(toga.QuestionDialog(title, message)))

    async def _ask_open_path(self):
        return await self.window.dialog(toga.OpenFileDialog(title="Open a Gaussian input",
                                                            file_types=INPUT_FILE_TYPES))

    def _status(self, text: str) -> None:
        _assign(self.status_label, "text", text)

    # --- steps --------------------------------------------------------------
    def _store_current(self) -> bool:
        """Keep the panels in self.steps; False (and a status message) when a number does not parse."""
        try:
            self.steps[self.current_step] = self._collect_step()
            return True
        except ValueError as exc:
            self._status(f"Fix step {self.current_step + 1} first: {exc}")
            return False

    def _show_step(self, index: int) -> None:
        self.current_step = index
        self._step_labels = []                       # force the list (and its selection) to be rewritten
        self._load_step(self.steps[index])
        self.refresh()

    def on_step_select(self, widget):
        if self._syncing or widget.value not in self._step_labels:
            return
        index = self._step_labels.index(widget.value)
        if index == self.current_step:
            return
        if not self._store_current():
            self._syncing = True
            try:
                _assign(self.step_selection, "value", self._step_labels[self.current_step])
            finally:
                self._syncing = False
            return
        self._show_step(index)

    def add_linked_step(self, widget):
        if not self._store_current():
            return
        step = self.steps[self.current_step]
        if step.job_type == gk.JOB_FUKUI:
            self._status("The Fukui functions job type is already three linked jobs; it cannot take more steps.")
            return
        self.steps.insert(self.current_step + 1, gi.linked_step(step))
        self._show_step(self.current_step + 1)
        self._status(f"Step {self.current_step + 1} added: it reads the geometry and orbitals of the checkpoint.")

    def remove_step(self, widget):
        if len(self.steps) == 1:
            self._status("This is the only step; there is nothing to remove.")
            return
        del self.steps[self.current_step]
        self._show_step(min(self.current_step, len(self.steps) - 1))

    def _move_step(self, delta: int) -> None:
        target = self.current_step + delta
        if not 0 <= target < len(self.steps) or not self._store_current():
            return
        self.steps[self.current_step], self.steps[target] = self.steps[target], self.steps[self.current_step]
        self._show_step(target)

    def move_step_up(self, widget):
        self._move_step(-1)

    def move_step_down(self, widget):
        self._move_step(+1)

    # --- molecule -----------------------------------------------------------
    def add_modredundant_line(self, widget):
        try:
            line = gi.modredundant_line(self.modred_type_selection.value, self.modred_atoms_input.value or "",
                                        self.modred_action_selection.value, self.modred_steps_input.value or "",
                                        self.modred_size_input.value or "")
        except ValueError as exc:
            self._status(str(exc))
            return
        text = (self.modredundant_input.value or "").rstrip()
        self.modredundant_input.value = f"{text}\n{line}" if text else line
        self._status(f"Added: {line}")
        self.refresh()

    async def load_geometry(self, widget):
        try:
            path = await self.window.dialog(toga.OpenFileDialog(title="Geometry of the molecule",
                                                                file_types=gi.GEOMETRY_FILE_TYPES))
            if path:
                self.apply_geometry_file(path)
        except ValueError as exc:
            await self._show_error("Geometry Error", str(exc))
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))

    def apply_geometry_file(self, path):
        """Fill the geometry (and the charge/multiplicity a .gjf/.fchk gives); keep a typed job name.

        The Fukui ion multiplicities are cleared, never guessed: they belong to the previous molecule.
        """
        path = Path(path)
        text, charge, mult = gi.read_geometry_file(path)
        self._syncing = True
        try:
            _assign(self.geometry_input, "value", text)
            if charge is not None:
                _assign(self.charge_input, "value", str(charge))
            if mult is not None:
                _assign(self.multiplicity_input, "value", str(mult))
            _assign(self.fukui_plus_mult_input, "value", "")
            _assign(self.fukui_minus_mult_input, "value", "")
            if (self.stem_input.value or "").strip() in ("", gi.fi.DEFAULT_STEM):
                stem = path.stem[:-2] if path.stem.endswith("_N") else path.stem
                _assign(self.stem_input, "value", "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in stem))
            if not (self.output_dir_input.value or "").strip():
                _assign(self.output_dir_input, "value", str(path.resolve().parent))
        finally:
            self._syncing = False
        self._status(f"Geometry loaded from {path.name}.")
        self.refresh()

    async def browse_output_dir(self, widget):
        try:
            folder = await self.window.dialog(toga.SelectFolderDialog(title="Folder for the Gaussian input"))
            if folder:
                self.output_dir_input.value = str(folder)
        except Exception as exc:
            await self._show_error("Folder Selection Error", str(exc))

    # --- whole job ------------------------------------------------------------
    def apply_input_file(self, path) -> List[str]:
        """Open input...: every step of a Gaussian input into the panels; returns the reader's notes."""
        spec, notes = gi.read_input_file(path)
        self.steps = spec.steps
        self._syncing = True
        try:
            _assign(self.stem_input, "value", spec.stem)
            _assign(self.output_dir_input, "value", spec.output_dir)
            _assign(self.geometry_input, "value", spec.geometry)
            _assign(self.formchk_switch, "value", spec.write_formchk)
        finally:
            self._syncing = False
        self._show_step(0)
        _assign(self.messages_output, "value", "\n".join(f"NOTE: {n}" for n in notes) or "Read without notes.")
        self._status(f"{Path(path).name}: {len(spec.steps)} step(s) read into the panels.")
        return notes

    async def open_input(self, widget):
        try:
            path = await self._ask_open_path()
            if path:
                self.apply_input_file(path)
        except ValueError as exc:
            await self._show_error("Open Input Error", str(exc))
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))

    def _prepare(self) -> gi.GaussianJob:
        return self.builder.prepare(self._collect_inputs())

    def _warnings_text(self, job: gi.GaussianJob) -> str:
        return "\n".join(f"WARNING: {w}" for w in job.warnings) if job.warnings else "No warnings."

    def _run_instructions(self, job: gi.GaussianJob) -> str:
        gjf, log = job.input_path.name, job.input_path.with_suffix(".log").name
        lines = ["Written:", f"  {job.input_path}"] + ([f"  {job.script_path}"] if job.script_path else []) + [
            "",
            "Next, in a Linux shell with Gaussian 16 set up (WSL, or the cluster after copying the folder):",
            f"  cd '{gi.wsl_path(job.input_path.parent)}'      (WSL; on a cluster: the copied folder)",
            f"  g16 < {gjf} > {log}      (or submit it through the cluster's queue)"]
        if job.script_path:
            lines.append(f"  bash {job.script_path.name}")
        if job.fukui is not None:
            lines.append(f"Then open {job.fukui.states[0].fchk.name} in Tools > Analysis of Fukui Functions; the N+1 "
                         "and N-1 files are filled in automatically (copy the three .fchk back first if they were "
                         "made on a cluster).")
        elif job.script_path:
            lines.append("Then open the .fchk in the 3D Molecular Viewer (Orbitals, ESP, QTAIM, NTO, Excitation tabs).")
        return "\n".join(lines) + "\n\n" + self._warnings_text(job)

    async def preview(self, widget):
        try:
            job = self._prepare()
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        self.job = job
        self.preview_output.value = self.builder.input_text(job)
        self.messages_output.value = self._warnings_text(job)
        self.tabs.current_tab = PREVIEW_TAB
        self._status(f"Preview of {job.input_path.name} (not saved yet).")

    async def save_input(self, widget):
        try:
            job = self._prepare()
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        existing = [p.name for p in (job.input_path, job.script_path) if p is not None and p.exists()]
        if existing and not await self._confirm("Replace Files?", f"{', '.join(existing)} already exist in "
                                                f"{job.input_path.parent}. Replace them?"):
            return
        try:
            self.builder.write_files(job)
        except ValueError as exc:
            await self._show_error("Write Error", str(exc))
            return
        self.job = job
        self.preview_output.value = self.builder.input_text(job)
        self.messages_output.value = self._run_instructions(job)
        self.tabs.current_tab = PREVIEW_TAB
        self._status(f"Saved {job.input_path.name}" + (f" and {job.script_path.name}." if job.script_path else "."))

    async def defaults(self, widget):
        """Every step back to one default step; the molecule (geometry, charge, multiplicity) is kept."""
        self._store_current()
        first = self.steps[0]
        reference = gi.GaussianStepSpec(charge=first.charge, multiplicity=first.multiplicity)
        if self.steps != [reference] and not await self._confirm(
                "Reset to Defaults?", "Every panel of every step goes back to its default (the geometry, "
                                      "charge and multiplicity are kept). Continue?"):
            return
        self.steps = [reference]
        self._show_step(0)
        self._status("Panels reset to their defaults.")

    def closeTopLevel(self, widget):
        self.window.close()
