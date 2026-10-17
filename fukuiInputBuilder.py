"""
Inputs > Gaussian Fukui Input Builder: the window over ``fukuiInput``.

The window only reads widgets into a ``FukuiJobSpec`` and turns every
``ValueError`` into a dialog; validation, the linked-job text and the formchk
script all live in ``fukuiInput`` (no Toga there). Spin multiplicities are
never filled in: a blank field reaches ``prepare`` as None and is asked for.
"""

from pathlib import Path, PureWindowsPath
from typing import Optional, Tuple

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import fukuiInput as fi
from coordinateConverter import BOHR_RADIUS_ANGSTROM, CoordinateConverter
from help import HelpGqteaWin
from systemProperties import allowed_multiplicities, multiplicity_problem, nuclear_charge, parse_integer

GEOMETRY_FILE_TYPES = ["xyz", "gjf", "com", "gau", "fchk", "fch", "pdb", "ent", "mol", "mol2", "hin"]
FCHK_SUFFIXES = (".fchk", ".fch")


def read_geometry_file(path) -> Tuple[str, Optional[int], Optional[int]]:
    """(Angstrom 'El x y z' text, charge or None, multiplicity or None) from a structure file.

    A Gaussian input or formatted checkpoint also gives the charge and multiplicity.
    """
    import gaussianWavefunction as gw     # numpy-only, but not needed until a .fchk is loaded
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in FCHK_SUFFIXES:
        try:
            sections = gw.read_fchk_sections(path.read_text(encoding="utf-8", errors="replace"))
            numbers = [int(z) for z in sections["Atomic numbers"]]
            xyz = list(sections["Current cartesian coordinates"])
        except (OSError, KeyError) as exc:
            raise ValueError(f"{path.name} is not a readable formatted checkpoint: {exc}") from exc
        atoms = [(fi.element_symbol(str(z)), *(BOHR_RADIUS_ANGSTROM * float(v) for v in xyz[3 * i:3 * i + 3]))
                 for i, z in enumerate(numbers)]
        return "\n".join(fi.format_geometry(atoms)), int(sections.get("Charge", 0)), \
            int(sections.get("Multiplicity", 1))
    converter = CoordinateConverter()
    fmt = converter.guess_input_format_from_filename(str(path))
    if fmt is None:
        raise ValueError(f"Cannot read '{path.suffix}' files; use one of: {', '.join(GEOMETRY_FILE_TYPES)}.")
    converter.set_input_file(str(path))
    if converter.input_unit != "Angstrom":
        raise ValueError(f"{path.name} declares bohr coordinates; convert it to Angstrom first "
                         "(Tools > Coordinate converter).")
    atoms = [(a["symbol"], a["x"], a["y"], a["z"]) for a in converter.read_atoms(fmt)]
    text = "\n".join(fi.format_geometry(atoms))
    fi.parse_geometry(text, check_distances=False)   # unknown elements fail here, naming the line
    charge = mult = None
    if fmt == "Gaussian":
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            tokens = line.split()
            if len(tokens) == 2 and all(t.lstrip("+-").isdigit() for t in tokens):
                charge, mult = int(tokens[0]), int(tokens[1])
                break
    return text, charge, mult


def wsl_path(path) -> str:
    """'C:\\jobs\\phenol' -> '/mnt/c/jobs/phenol' (the folder as WSL sees it)."""
    p = PureWindowsPath(str(path))
    if not p.drive or not p.drive.endswith(":"):
        return str(path)
    return "/mnt/" + p.drive[0].lower() + "/" + "/".join(p.parts[1:])


class FukuiInputBuilderUI:
    LABEL_WIDTH = 200
    HINT_COLOR = "#666666"

    def __init__(self, *args, **kwargs):
        # gqteaWinToga passes the Button that was pressed; ignore it.
        self.builder = fi.FukuiInputBuilder()
        self.job: Optional[fi.FukuiJob] = None
        self.window = toga.Window(title="Gaussian Fukui Input Builder", size=(860, 660))
        self.layout_main_window()
        self.window.show()

    # --- layout -----------------------------------------------------------
    def _label(self, text):
        return toga.Label(text, style=Pack(width=self.LABEL_WIDTH, margin=(6, 6, 0, 0)))

    def _row(self, *children):
        return toga.Box(children=list(children), style=Pack(direction=ROW, margin_bottom=4))

    def _heading(self, text):
        return toga.Label(text, style=Pack(font_weight="bold", margin=(10, 0, 4, 0)))

    def _hint(self, text):
        return toga.Label(text, style=Pack(color=self.HINT_COLOR, font_size=9, margin=(0, 0, 4, 0)))

    def _text(self, value="", placeholder="", width=None, on_change=None):
        style = Pack(width=width) if width else Pack(flex=1)
        return toga.TextInput(value=value, placeholder=placeholder, style=style, on_change=on_change)

    def layout_main_window(self):
        # Labels first: widgets below fire on_change -> refresh_states while being built.
        self.states_label = toga.Label("", style=Pack(font_family="monospace", margin=(4, 0, 4, 0)))
        self.route_label = toga.Label("", style=Pack(font_family="monospace", margin=(6, 0, 0, 0), flex=1))
        self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))

        refresh = self.refresh_states
        self.stem_input = self._text(fi.DEFAULT_STEM, width=260, on_change=refresh)
        self.output_dir_input = self._text(placeholder="Required: folder for the .gjf and .sh files")
        self.memory_input = self._text(fi.DEFAULT_MEMORY, width=120)
        self.nproc_input = self._text(str(fi.DEFAULT_NPROC), width=120)
        self.print_level_selection = toga.Selection(items=list(fi.PRINT_LEVELS), style=Pack(width=120),
                                                    on_change=refresh)
        self.method_input = self._text(fi.DEFAULT_METHOD, width=200, on_change=refresh)
        self.basis_input = self._text(fi.DEFAULT_BASIS, width=200, on_change=refresh)
        self.keywords_input = self._text(fi.DEFAULT_KEYWORDS, on_change=refresh)
        self.opt_input = self._text(fi.DEFAULT_OPT_KEYWORD, width=200, on_change=refresh)
        setup = toga.Box(children=[
            self._hint("For Gaussian 16 (Rev. C.01) under WSL or on a Linux cluster: relative file names, "
                       "LF line endings; run the job from the folder that holds the input."),
            self._heading("Job"),
            self._row(self._label("Job name (file stem):"), self.stem_input),
            self._hint("Files: <name>_fukui.gjf, <name>_N.chk / _Np1.chk / _Nm1.chk, <name>_formchk.sh."),
            self._row(self._label("Output folder:"), self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
            self._heading("Link 0 (set these for the machine that runs Gaussian)"),
            self._row(self._label("%Mem:"), self.memory_input,
                      self._hint("   integer + unit, e.g. 2GB, 500MB (a bare number means words)")),
            self._row(self._label("%NProcShared:"), self.nproc_input,
                      self._hint("   cores of the WSL machine or the cluster node")),
            self._heading("Route section"),
            self._row(self._label("Print level:"), self.print_level_selection),
            self._row(self._label("Method:"), self.method_input,
                      self._hint("   no R/U prefix needed: Gaussian picks R or U per state")),
            self._row(self._label("Basis set:"), self.basis_input),
            self._row(self._label("Extra keywords (all jobs):"), self.keywords_input),
            self._row(self._label("Optimization (job 1):"), self.opt_input,
                      self._hint("   e.g. Opt, Opt=Tight, Opt=(Tight,CalcFC)")),
            self._hint("Job 1 optimizes the N state; jobs 2-3 are single points with Geom=Check. "
                       "Scan, IRC, Geom=, Guess=, Gen and ECPs are refused."),
            self._row(self._label("Route:"), self.route_label),
        ], style=Pack(direction=COLUMN, margin=10))

        self.geometry_input = toga.MultilineTextInput(
            placeholder="Starting geometry of the N state, Cartesian, Angstrom:\n"
                        "C   0.000000   0.000000   0.000000\n...",
            style=Pack(flex=1, height=180, font_family="monospace"), on_change=refresh)
        self.charge_input = self._text("0", width=80, on_change=refresh)
        self.multiplicity_input = self._text(placeholder="required", width=80, on_change=refresh)
        self.p_input = self._text("1", width=80, on_change=refresh)
        self.q_input = self._text("1", width=80, on_change=refresh)
        self.plus_mult_input = self._text(placeholder="required", width=80, on_change=refresh)
        self.minus_mult_input = self._text(placeholder="required", width=80, on_change=refresh)
        self.guess_read_switch = toga.Switch("Start the ions from the optimized N-state orbitals (Guess=Read)",
                                             value=True, on_change=refresh)
        molecule = toga.Box(children=[
            self._row(toga.Button("Load geometry...", on_press=self.load_geometry, style=Pack(width=150)),
                      self._hint("   .xyz .gjf .fchk .pdb .mol .mol2 .hin; a .gjf/.fchk also sets the N "
                                 "charge and multiplicity")),
            self.geometry_input,
            self._hint("Job 1 optimizes this geometry; jobs 2-3 read the optimized one (Geom=Check)."),
            self._heading("Electron-count states"),
            self._row(self._label("Charge of the N state:"), self.charge_input),
            self._row(self._label("Multiplicity of the N state:"), self.multiplicity_input),
            self._row(self._label("Multiplicity of N+p:"), self.plus_mult_input,
                      self._label("   Multiplicity of N-q:"), self.minus_mult_input),
            self._hint("All three multiplicities are required: spin states are not assumed. Use the "
                       "lowest-energy spin state of each ion."),
            self._row(self._label("p (electrons added):"), self.p_input,
                      self._label("   q (electrons removed):"), self.q_input),
            self._hint("1 and 1 for ordinary f+ / f-; larger only for degenerate frontier orbitals."),
            self.states_label,
            self.guess_read_switch,
            self._hint("Both ions read %OldChk=<name>_N.chk; each writes its own .chk, so the N "
                       "checkpoint is kept."),
        ], style=Pack(direction=COLUMN, margin=10))

        self.preview_output = toga.MultilineTextInput(readonly=True, style=Pack(flex=2, font_family="monospace"))
        self.messages_output = toga.MultilineTextInput(readonly=True, style=Pack(flex=1, font_family="monospace"))
        output = toga.Box(children=[self.preview_output, self.messages_output],
                          style=Pack(direction=COLUMN, margin=10, flex=1))

        help_box = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_fukui_input,
                                           style=Pack(flex=1))
        self.tabs = toga.OptionContainer(content=[
            ("Setup", toga.ScrollContainer(content=setup, style=Pack(flex=1))),
            ("Molecule & states", toga.ScrollContainer(content=molecule, style=Pack(flex=1))),
            ("Preview", output),
            ("Help", help_box),
        ], style=Pack(flex=1))
        buttons = self._row(
            toga.Button("Preview", on_press=self.preview, style=Pack(width=110, margin=5)),
            toga.Button("Save input", on_press=self.save_input, style=Pack(width=110, margin=5)),
            toga.Button("Close", on_press=self.closeTopLevel, style=Pack(width=110, margin=5)))
        self.print_level_selection.value = fi.DEFAULT_PRINT_LEVEL
        self.refresh_states(None)
        self.window.content = toga.Box(children=[self.tabs, buttons, self.status_label],
                                       style=Pack(direction=COLUMN, margin=10, flex=1))

    # --- dialogs ------------------------------------------------------------
    async def _show_error(self, title: str, message: str):
        await self.window.dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        await self.window.dialog(toga.InfoDialog(title, message))

    async def _confirm(self, title: str, message: str) -> bool:
        return bool(await self.window.dialog(toga.QuestionDialog(title, message)))

    async def browse_output_dir(self, widget):
        try:
            folder = await self.window.dialog(toga.SelectFolderDialog(title="Folder for the Fukui inputs"))
            if folder:
                self.output_dir_input.value = str(folder)
        except Exception as exc:
            await self._show_error("Folder Selection Error", str(exc))

    async def load_geometry(self, widget):
        try:
            path = await self.window.dialog(toga.OpenFileDialog(title="Starting geometry of the N state",
                                                                file_types=GEOMETRY_FILE_TYPES))
            if path:
                self.apply_geometry_file(path)
        except ValueError as exc:
            await self._show_error("Geometry Error", str(exc))
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))

    # --- form <-> spec ------------------------------------------------------
    def apply_geometry_file(self, path):
        """Fill the geometry (and the N charge/multiplicity a .gjf/.fchk gives); keep a typed job name.

        The ion multiplicities are cleared, never guessed: they belong to the previous molecule.
        """
        path = Path(path)
        text, charge, mult = read_geometry_file(path)
        self.geometry_input.value = text
        if charge is not None:
            self.charge_input.value = str(charge)
        if mult is not None:
            self.multiplicity_input.value = str(mult)
        self.plus_mult_input.value = self.minus_mult_input.value = ""
        if (self.stem_input.value or "").strip() in ("", fi.DEFAULT_STEM):
            stem = path.stem[:-2] if path.stem.endswith("_N") else path.stem
            self.stem_input.value = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in stem)
        if not (self.output_dir_input.value or "").strip():
            self.output_dir_input.value = str(path.resolve().parent)
        self.status_label.text = f"Geometry loaded from {path.name}."
        self.refresh_states(None)

    def _multiplicities(self):
        """(M, M+, M-): integers, or None where the field is blank (never filled in)."""
        values = []
        for widget, label in ((self.multiplicity_input, "the N state"), (self.plus_mult_input, "N+p"),
                              (self.minus_mult_input, "N-q")):
            values.append(parse_integer(widget.value, f"Multiplicity of {label}"))
        return tuple(values)

    def _collect_inputs(self) -> fi.FukuiJobSpec:
        charge = parse_integer(self.charge_input.value, "Charge of the N state", 0)
        p = parse_integer(self.p_input.value, "p (electrons added)", 1)
        q = parse_integer(self.q_input.value, "q (electrons removed)", 1)
        m, m_plus, m_minus = self._multiplicities()
        return fi.FukuiJobSpec(
            stem=(self.stem_input.value or "").strip(),
            output_dir=(self.output_dir_input.value or "").strip(),
            geometry=self.geometry_input.value or "",
            charge=charge, multiplicity=m, plus_multiplicity=m_plus, minus_multiplicity=m_minus,
            memory=self.memory_input.value or "",
            nproc=fi.parse_nproc(self.nproc_input.value),
            print_level=self.print_level_selection.value or fi.DEFAULT_PRINT_LEVEL,
            method=self.method_input.value or "", basis=self.basis_input.value or "",
            keywords=self.keywords_input.value or "",
            opt_keyword=self.opt_input.value or "",
            p=p, q=q,
            guess_read=bool(self.guess_read_switch.value))

    def refresh_states(self, widget):
        """Live route lines and the charge / multiplicity / electrons / .chk of each state."""
        if not hasattr(self, "states_label"):
            return
        try:
            args = (self.print_level_selection.value or fi.DEFAULT_PRINT_LEVEL, self.method_input.value or "",
                    self.basis_input.value or "", self.keywords_input.value or "")
            ion_extra = ("Geom=Check", "Guess=Read") if self.guess_read_switch.value else ("Geom=Check",)
            self.route_label.text = (
                "job 1:     " + " ".join(fi.build_route(*args, ((self.opt_input.value or "").strip(),))) + "\n"
                "jobs 2-3:  " + " ".join(fi.build_route(*args, ion_extra)))
        except AttributeError:
            return   # still building the window
        try:
            charge = parse_integer(self.charge_input.value, "Charge of the N state", 0)
            p = parse_integer(self.p_input.value, "p", 1)
            q = parse_integer(self.q_input.value, "q", 1)
            if (self.geometry_input.value or "").strip() == "":
                self.states_label.text = "Enter or load a geometry to see the electron counts of the states."
                return
            mults = self._multiplicities()
            atoms = fi.parse_geometry(self.geometry_input.value, check_distances=False)
            z_total = nuclear_charge(a[0] for a in atoms)
        except (ValueError, AttributeError) as exc:
            self.states_label.text = f"Cannot count electrons yet: {exc}"
            return
        stem = (self.stem_input.value or "").strip() or fi.DEFAULT_STEM
        lines = []
        for (label, tag, _, k, sign), mult in zip(fi.state_tags(p, q), mults):
            state_charge = charge + sign * k
            n_e = z_total - state_charge
            if mult is None:
                allowed = ", ".join(str(m) for m in allowed_multiplicities(n_e, 3))
                spin = f"multiplicity ? (enter it: {n_e} electrons allow {allowed}, ...)"
            else:
                spin = f"multiplicity {mult}"
                problem = multiplicity_problem(n_e, mult) if n_e >= 1 else "no electrons left"
                if problem:
                    spin += f"   <- impossible: {problem}"
            lines.append(f"{label:<5} charge {state_charge:+d}, {n_e} electrons -> {stem}_{tag}.chk, {spin}")
        self.states_label.text = "\n".join(lines).replace("charge +0", "charge 0")

    # --- actions ------------------------------------------------------------
    def _prepare(self) -> fi.FukuiJob:
        return self.builder.prepare(self._collect_inputs())

    def _warnings_text(self, job: fi.FukuiJob) -> str:
        return "\n".join(f"WARNING: {w}" for w in job.warnings) if job.warnings else "No warnings."

    async def preview(self, widget):
        try:
            job = self._prepare()
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        self.job = job
        self.preview_output.value = self.builder.input_text(job)
        self.messages_output.value = self._warnings_text(job)
        self.tabs.current_tab = 2
        self.status_label.text = f"Preview of {job.input_path.name} (not saved yet)."

    async def save_input(self, widget):
        try:
            job = self._prepare()
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        existing = [p.name for p in (job.input_path, job.script_path) if p.exists()]
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
        gjf, sh = job.input_path.name, job.script_path.name
        log = job.input_path.with_suffix(".log").name
        self.messages_output.value = (
            f"Written:\n  {job.input_path}\n  {job.script_path}\n\n"
            "Next, in a Linux shell with Gaussian 16 set up (WSL, or the cluster after copying the folder):\n"
            f"  cd '{wsl_path(job.input_path.parent)}'      (WSL; on a cluster: the copied folder)\n"
            f"  g16 < {gjf} > {log}      (or submit it through the cluster's queue)\n"
            f"  bash {sh}\n"
            f"Then open {job.states[0].fchk.name} in Tools > Analysis of Fukui Functions; the N+1 and N-1 "
            "files are filled in automatically (copy the three .fchk back first if they were made on a "
            "cluster).\n\n" + self._warnings_text(job))
        self.tabs.current_tab = 2
        self.status_label.text = f"Saved {gjf} and {sh}."

    def closeTopLevel(self, widget):
        self.window.close()
