"""The 3D Molecular Viewer's Excitation tab (``ExcitationTabController``).

Electron excitation analysis of one CIS/TDA/TD state from a Gaussian .fchk and
the output of the same job, in the style of Multiwfn's main function 18:
hole-electron indices, MO/atom/fragment contributions, IFCT, Delta_r, Lambda,
ghost-hunter, Coulomb attraction and a descriptor-based classification. The
numerical work is ``excitationAnalysis.ExcitationAnalyser`` (headless). Hole,
electron, Sr, charge density difference, transition density and the C functions
are drawn in the viewer's own 3D window through ``viewer.excitation_layer`` (a
third ``cubeRender.CubeSceneLayer``), with spheres at the hole and electron
centroids and an arrow between them. Conventions as in the NTO tab: widget
writes through ``cubeTab._assign``, synchronous ``_syncing``-guarded change
handlers that start their async work with ``_start``.
"""

from __future__ import annotations

import asyncio
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import excitationAnalysis as xa
import excitationData as ed
import volumeData as vd
from cubeRender import arrow_item, sphere_item, surface_item
from cubeTab import CubeEntry, _assign, _item_names, canonical_colour
from help import HelpGqteaWin
from qtaimTab import fit_table_columns, parse_optional_float

A = xa.BOHR_RADIUS_ANGSTROM
SHOW_MODES = {
    "Hole and electron": ("hole", "electron"),
    "Hole": ("hole",),
    "Electron": ("electron",),
    "Sr overlap": ("Sr",),
    "Charge density difference": ("CDD",),
    "Transition density": ("transdens",),
    "C_hole and C_ele": ("Chole", "Cele"),
}
DEFAULT_SHOW = "Hole and electron"
PARTITION_LABELS = {"Mulliken-like": xa.PARTITION_MULLIKEN, "Hirshfeld": xa.PARTITION_HIRSHFELD}
COLOUR_KEYS = (("hole", 1), ("electron", 1), ("Sr", 1), ("CDD", 1), ("CDD", -1), ("transdens", 1),
               ("transdens", -1))
DEFAULT_COLOURS = {("hole", 1): "Blue", ("electron", 1): "Green", ("Sr", 1): "Orange", ("CDD", 1): "Green",
                   ("CDD", -1): "Blue", ("transdens", 1): "Cyan", ("transdens", -1): "Magenta"}
COLOUR_LABELS = {("hole", 1): "Hole", ("electron", 1): "Electron", ("Sr", 1): "Sr", ("CDD", 1): "CDD +",
                 ("CDD", -1): "CDD −", ("transdens", 1): "Transition +", ("transdens", -1): "Transition −"}
COLOUR_ALIASES = {"Chole": ("hole", 1), "Cele": ("electron", 1)}   # C functions take the hole/electron colours
CENTROID_RADIUS = 0.25 / A     # bohr (Multiwfn's VMD example draws 0.25 A spheres)
ARROW_RADIUS = 0.06 / A
MIN_ARROW_LENGTH = 0.1 / A
DESCRIPTOR_HEADINGS = ["Quantity", "Value", "Unit"]
DESCRIPTOR_ACCESSORS = ["quantity", "value", "unit"]
DESCRIPTOR_FRACTIONS = (0.55, 0.28, 0.17)
PART_HEADINGS = ["Atom / fragment", "Hole %", "Electron %", "Overlap %", "Δ %"]
PART_ACCESSORS = ["part", "hole", "electron", "overlap", "difference"]
PART_FRACTIONS = (0.32, 0.17, 0.17, 0.17, 0.17)
IMAGE_FILE_TYPES = ["png", "tif", "tiff", "bmp"]
MIN_OPACITY = 0.1


def colour_key(name: str, sign: int) -> Tuple[str, int]:
    return COLOUR_ALIASES.get(name, (name, sign))


class ExcitationTabController:
    """The 3D Molecular Viewer's Excitation tab."""

    LABEL_WIDTH = 210
    HINT_COLOR = "#666666"

    def __init__(self, viewer):
        self.viewer = viewer
        self.analyser = xa.ExcitationAnalyser()
        self.loaded = None                       # (wavefunction, ExcitationFile, notes) of loaded_paths
        self.loaded_paths: Optional[Tuple[str, str]] = None
        self.results: Optional[xa.ExcitationResults] = None
        self.report_path: Optional[str] = None
        self.show = DEFAULT_SHOW
        self.colours: Dict[Tuple[str, int], str] = dict(DEFAULT_COLOURS)
        self.opacity = 1.0
        self.show_centroids = True
        self.use_numba: Optional[bool] = None    # None = if installed
        self._state_items: Dict[str, int] = {}
        self._volumes: Optional[Tuple[tuple, Dict[str, CubeEntry]]] = None  # (key, {field: entry})
        self._structure_settled_for = None
        self._syncing = False
        self._refreshing = False
        self._refresh_pending = False
        self._tasks = set()

    # --- layout -----------------------------------------------------------
    def _label(self, text):
        return toga.Label(text, style=Pack(width=self.LABEL_WIDTH, margin=(6, 6, 0, 0)))

    def _row(self, *children):
        return toga.Box(children=list(children), style=Pack(direction=ROW, margin_bottom=4))

    def _heading(self, text):
        return toga.Label(text, style=Pack(font_weight="bold", margin=(10, 0, 4, 0)))

    def _hint(self, text):
        return toga.Label(text, style=Pack(color=self.HINT_COLOR, font_size=9, margin=(0, 0, 4, 0)))

    def build_tab(self) -> toga.Box:
        self._syncing = True
        try:
            self.fchk_input = toga.TextInput(placeholder="Gaussian formatted checkpoint (.fchk / .fch)",
                                             style=Pack(flex=1))
            self.excitation_input = toga.TextInput(placeholder="Gaussian CIS/TD output (.log / .out) or plain text",
                                                   style=Pack(flex=1))
            self.output_dir_input = toga.TextInput(placeholder="Blank = folder of the .fchk file", style=Pack(flex=1))
            self.load_button = toga.Button("Load", on_press=self.load_files, style=Pack(width=110, margin=5))
            self.state_selection = toga.Selection(items=[], style=Pack(flex=1))
            self.partition_selection = toga.Selection(items=list(PARTITION_LABELS), style=Pack(width=200))
            self.partition_selection.value = "Mulliken-like"
            self.fragments_input = toga.TextInput(placeholder="e.g. 11-13; 1-10; 14-16 (blank = no fragments)",
                                                  style=Pack(flex=1))
            self.hydrogens_switch = toga.Switch("Include hydrogens in atom tables and heat maps", value=False,
                                                on_change=self.on_hydrogens_change, style=Pack(margin=(4, 0, 4, 0)))
            self.fast_coulomb_switch = toga.Switch("Exact Coulomb attraction - fast (about 10x faster, ~1e-4 "
                                                   "relative error)", value=False,
                                                   on_change=self.on_coulomb_mode_change,
                                                   style=Pack(margin=(4, 0, 2, 0)))
            self.precise_coulomb_switch = toga.Switch("Exact Coulomb attraction - precise (slower, full grid)",
                                                      value=False, on_change=self.on_coulomb_mode_change,
                                                      style=Pack(margin=(2, 0, 4, 0)))
            self.run_button = toga.Button("Run", on_press=self.workflow, style=Pack(width=110, margin=5))
            self.clear_button = toga.Button("Clear", on_press=self.clear_output, style=Pack(width=110, margin=5))
            self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
            self.progress_bar = toga.ProgressBar(max=100, style=Pack(margin=(0, 5, 5, 5)))
            self.classification_label = toga.Label("No results.", style=Pack(font_weight="bold", margin=(4, 0, 2, 0)))
            self.summary_label = toga.Label("", style=Pack(margin=(0, 0, 4, 0)))
            self.descriptor_table = toga.Table(headings=DESCRIPTOR_HEADINGS, accessors=DESCRIPTOR_ACCESSORS, data=[],
                                               style=Pack(height=230))
            self.part_table = toga.Table(headings=PART_HEADINGS, accessors=PART_ACCESSORS, data=[],
                                         style=Pack(height=170))

            self.show_selection = toga.Selection(items=list(SHOW_MODES), on_change=self.on_show_change,
                                                 style=Pack(width=230))
            self.show_selection.value = self.show
            self.isovalue_input = toga.TextInput(placeholder="0.002 (transition density 0.001)",
                                                 on_confirm=self.on_grid_confirm, style=Pack(width=230))
            self.spacing_input = toga.TextInput(placeholder=f"{xa.DEFAULT_SPACING_BOHR:g}",
                                                on_confirm=self.on_grid_confirm, style=Pack(width=120))
            self.colour_inputs = {key: toga.TextInput(value=DEFAULT_COLOURS[key], on_confirm=self.on_colour_confirm,
                                                      style=Pack(width=110)) for key in COLOUR_KEYS}
            self.opacity_slider = toga.Slider(min=MIN_OPACITY, max=1.0, value=1.0, on_change=self.on_opacity_change,
                                              style=Pack(width=200))
            self.centroids_switch = toga.Switch("Centroids and CT arrow", value=True,
                                                on_change=self.on_centroids_change, style=Pack(margin=(4, 0, 4, 0)))
            self.display_button = toga.Button("Display", on_press=self.display_excitation,
                                              style=Pack(width=110, margin=5))

            self.export_cubes_button = toga.Button("Export cubes", on_press=self.export_cubes,
                                                   style=Pack(width=140, margin=5))
            self.send_button = toga.Button("Send to Cube tab", on_press=self.send_to_cube_tab,
                                           style=Pack(width=140, margin=5))
            self.export_3d_button = toga.Button("Export 3D image...", on_press=self.export_3d_image,
                                                style=Pack(width=150, margin=5))
            self.heatmap_selection = toga.Selection(items=list(xa.HEATMAP_KINDS.values()), style=Pack(width=200))
            self.heatmap_button = toga.Button("Plot heat map", on_press=self.plot_heat_map,
                                              style=Pack(width=140, margin=5))
            self.report_output = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_excitation,
                                                         style=Pack(height=320, font_family="monospace"))
        finally:
            self._syncing = False

        colour_rows = [self._row(self._label(" / ".join(COLOUR_LABELS[k] for k in keys) + ":"),
                                 *[w for k in keys for w in (self.colour_inputs[k],
                                                             toga.Label(" ", style=Pack(width=8)))])
                       for keys in ((("hole", 1), ("electron", 1), ("Sr", 1)), (("CDD", 1), ("CDD", -1)),
                                    (("transdens", 1), ("transdens", -1)))]
        children = [
            self._hint("Electron excitation analysis (hole-electron, IFCT, Delta_r, Lambda) of a CIS/TDA/TD state, "
                       "as Multiwfn's main function 18; drawn in this viewer's 3D window."),
            self._heading("Input"),
            self._row(self._label("Formatted checkpoint:"), self.fchk_input,
                      toga.Button("Browse", on_press=self.browse_fchk, style=Pack(width=90, margin_left=5))),
            self._row(self._label("Excitation file:"), self.excitation_input,
                      toga.Button("Browse", on_press=self.browse_excitation, style=Pack(width=90, margin_left=5))),
            self._hint("The Gaussian output of the same job, run with IOp(9/40=4); .log/.out beside the .fchk is "
                       "filled in."),
            self._row(self._label("Output folder:"), self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
            self._row(self.load_button),
            self._heading("Excited state and options"),
            self._row(self._label("State:"), self.state_selection),
            self._row(self._label("Atom partition:"), self.partition_selection),
            self._hint("Mulliken-like is fast; use Hirshfeld (H-Ar, all-electron) with diffuse basis functions."),
            self._row(self._label("Fragments (IFCT):"), self.fragments_input),
            self._hint("1-based atoms, ';' between fragments; atoms left out form the fragment 'Rest'."),
            self.hydrogens_switch,
            self.fast_coulomb_switch,
            self.precise_coulomb_switch,
            self._hint("Tick one to also compute the Coulomb attraction of the hole and electron densities "
                       "(the C-function value is always given)."),
            self._row(self.run_button, self.clear_button),
            self.status_label,
            self.progress_bar,
            self._heading("Results"),
            self.classification_label,
            self.summary_label,
            self.descriptor_table,
            self.part_table,
            self._heading("Display"),
            self._row(self._label("Show:"), self.show_selection),
            self._row(self._label("Isovalue (a.u.):"), self.isovalue_input),
            self._row(self._label("Grid spacing (bohr):"), self.spacing_input),
            *colour_rows,
            self._hint("C_hole / C_ele use the hole / electron colours. Colours: Blue, Red, Green, Yellow, Orange, "
                       "Cyan, Magenta, Purple, White, Grey or #rrggbb (Enter applies)."),
            self._row(self._label("Opacity:"), self.opacity_slider),
            self.centroids_switch,
            self._row(self.display_button),
            self._heading("Export"),
            self._row(self.export_cubes_button, self.send_button, self.export_3d_button),
            self._row(self._label("Heat map:"), self.heatmap_selection, self.heatmap_button),
            self._hint("The report is written after every run; cubes and heat maps go to the output folder."),
            self.report_output,
        ]
        self._attach_column_fitting()
        return toga.Box(children=children, style=Pack(direction=COLUMN, margin=10))

    def _attach_column_fitting(self):
        for table, fractions in ((self.descriptor_table, DESCRIPTOR_FRACTIONS), (self.part_table, PART_FRACTIONS)):
            native = getattr(getattr(table, "_impl", None), "native", None)
            if native is None:
                continue
            try:
                native.Resize += lambda sender, event, t=table, f=fractions: fit_table_columns(t, f)
            except Exception:
                pass

    # --- dialogs ------------------------------------------------------------
    def _window(self):
        return self.viewer.main_window

    async def _show_error(self, title: str, message: str):
        await self._window().dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        await self._window().dialog(toga.InfoDialog(title, message))

    async def _ask(self, title: str, message: str) -> bool:
        return bool(await self._window().dialog(toga.QuestionDialog(title, message)))

    async def _ask_save_path(self, title: str, suggested: str, file_types):
        return await self._window().dialog(toga.SaveFileDialog(title=title, suggested_filename=suggested,
                                                               file_types=list(file_types)))

    async def browse_fchk(self, widget):
        try:
            path = await self._window().dialog(toga.OpenFileDialog(
                title="Select a Gaussian formatted checkpoint", file_types=["fchk", "fch"]))
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))
            return
        if path:
            self.set_fchk_path(str(path))

    def set_fchk_path(self, path: str):
        """Set the .fchk and, when the excitation field is empty, its <stem>.log / <stem>.out."""
        self.fchk_input.value = path
        if not (self.excitation_input.value or "").strip():
            sibling = ed.sibling_excitation_file(path)
            if sibling is not None:
                self.excitation_input.value = str(sibling)

    async def browse_excitation(self, widget):
        try:
            path = await self._window().dialog(toga.OpenFileDialog(
                title="Select the Gaussian CIS/TD output", file_types=["log", "out", "txt"]))
            if path:
                self.excitation_input.value = str(path)
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))

    async def browse_output_dir(self, widget):
        try:
            folder = await self._window().dialog(toga.SelectFolderDialog(title="Choose the output folder"))
            if folder:
                self.output_dir_input.value = str(folder)
        except Exception as exc:
            await self._show_error("Folder Selection Error", str(exc))

    def _start(self, coroutine):
        """Run a handler's async work on the UI loop, keeping the task referenced until done."""
        loop = getattr(self.viewer, "loop", None) or asyncio.get_event_loop()
        task = loop.create_task(coroutine)
        self._tasks.add(task)
        task.add_done_callback(self._tasks.discard)
        return task

    # --- inputs ---------------------------------------------------------------
    def _paths(self) -> Tuple[str, str]:
        fchk = (self.fchk_input.value or "").strip()
        excitation = (self.excitation_input.value or "").strip()
        if not fchk:
            raise ValueError("Please select a .fchk or .fch file.")
        if not excitation:
            raise ValueError("Please select the Gaussian CIS/TD output (.log/.out) of the same job.")
        for path in (fchk, excitation):
            if not os.path.isfile(path):
                raise ValueError(f"File not found: {path}")
        return fchk, excitation

    def _collect_inputs(self) -> Dict[str, object]:
        """Exactly the keyword parameters of ExcitationAnalyser.run (a test pins this)."""
        fchk, excitation = self._paths()
        if self.loaded_paths != (fchk, excitation):
            raise ValueError("Press Load to read these files and list their excited states.")
        number = self._state_items.get(self.state_selection.value)
        if number is None:
            raise ValueError("Choose an excited state.")
        partition = PARTITION_LABELS.get(self.partition_selection.value, xa.PARTITION_MULLIKEN)
        fragments = (self.fragments_input.value or "").strip()
        xa.parse_fragments(fragments, self.loaded[0].n_atoms)  # report a bad definition before running
        return {"fchk_path": fchk, "excitation_path": excitation, "state_number": number, "partition": partition,
                "fragments_text": fragments, "coulomb_mode": self.coulomb_mode()}

    def coulomb_mode(self) -> str:
        """The exact-Coulomb mode chosen with the two mutually exclusive checkboxes."""
        if self.fast_coulomb_switch.value:
            return xa.COULOMB_FAST
        if self.precise_coulomb_switch.value:
            return xa.COULOMB_PRECISE
        return xa.COULOMB_OFF

    def on_coulomb_mode_change(self, widget, **kwargs):
        """Ticking one Coulomb checkbox disables the other; unticking it enables the other again."""
        other = self.precise_coulomb_switch if widget is self.fast_coulomb_switch else self.fast_coulomb_switch
        other.enabled = not bool(widget.value)

    def _collect_display_inputs(self) -> Dict[str, object]:
        """The grid parameters of excitationAnalysis.field_volumes; isovalue None = each field's default."""
        if self.results is None:
            raise ValueError("Run the analysis first.")
        spacing = parse_optional_float(self.spacing_input.value, xa.DEFAULT_SPACING_BOHR, "The grid spacing")
        text = (self.isovalue_input.value or "").strip()
        isovalue = parse_optional_float(text, xa.DEFAULT_ISOVALUE, "The isovalue") if text else None
        return {"spacing": spacing, "isovalue": isovalue}

    def _collect_colours(self) -> Dict[Tuple[str, int], str]:
        colours = {}
        for key, widget in self.colour_inputs.items():
            try:
                colours[key] = canonical_colour(widget.value)
            except ValueError as exc:
                raise ValueError(f"{COLOUR_LABELS[key]}: {exc}") from None
        return colours

    def output_dir(self) -> str:
        return xa.resolve_output_dir(self.results.fchk_path, self.output_dir_input.value)

    # --- progress -------------------------------------------------------------
    def _apply_progress(self, fraction: float, message: str):
        _assign(self.progress_bar, "value", max(0.0, min(100.0, 100.0 * fraction)))
        _assign(self.status_label, "text", message)

    def _progress_callback(self, message: Optional[str] = None):
        loop = asyncio.get_running_loop()

        def on_progress(fraction, text=None):
            loop.call_soon_threadsafe(self._apply_progress, fraction, text or message or "")
        return on_progress

    def _set_buttons(self, enabled: bool):
        for b in (self.load_button, self.run_button, self.clear_button, self.display_button,
                  self.export_cubes_button, self.send_button, self.heatmap_button):
            b.enabled = enabled

    # --- load and run -----------------------------------------------------------
    async def load_files(self, widget):
        try:
            paths = self._paths()
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        if not await self._load(paths):
            return
        wfn, data, notes = self.loaded
        kind = "unrestricted" if data.unrestricted else "restricted"
        message = f"{len(data.states)} excited states read ({kind}, {wfn.n_atoms} atoms). Choose one and press Run."
        _assign(self.status_label, "text", " ".join([message] + notes))

    async def _load(self, paths: Tuple[str, str]) -> bool:
        self._set_buttons(False)
        _assign(self.status_label, "text", "Reading the files...")
        try:
            self.loaded = await asyncio.to_thread(self.analyser.load, *paths)
        except (ValueError, MemoryError) as exc:
            self.loaded, self.loaded_paths = None, None
            self._fill_states([])
            _assign(self.status_label, "text", "The files were not loaded.")
            await self._show_error("Excitation Input Error", str(exc))
            return False
        finally:
            self._set_buttons(True)
        self.loaded_paths = paths
        self._fill_states(self.loaded[1].states)
        return True

    def _fill_states(self, states: List[ed.ExcitedState]):
        previous = self._state_items.get(self.state_selection.value)
        self._state_items = {s.title(): s.number for s in states}
        self._syncing = True
        try:
            if list(self._state_items) != _item_names(self.state_selection):
                self.state_selection.items = list(self._state_items)
            if self._state_items:
                keep = next((t for t, n in self._state_items.items() if n == previous), None)
                _assign(self.state_selection, "value", keep or next(iter(self._state_items)))
        finally:
            self._syncing = False

    async def workflow(self, widget):
        try:
            paths = self._paths()
            if self.loaded_paths != paths and not await self._load(paths):
                return
            params = self._collect_inputs()
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        self._set_buttons(False)
        self.progress_bar.value = 0
        try:
            results = await asyncio.to_thread(self.analyser.run, progress=self._progress_callback(), **params)
        except Exception as exc:
            _assign(self.status_label, "text", "Calculation failed.")
            self.progress_bar.value = 0
            self._set_buttons(True)
            await self._show_error("Excitation Analysis Error", str(exc))
            return
        was_displayed = self._scene_is_displayed()
        self.results = results  # kept even if the report cannot be written
        self._volumes = None
        self._fill_results()
        try:
            self.report_path = self.analyser.write_report(results, self.output_dir())
            self.report_output.value = Path(self.report_path).read_text(encoding="utf-8")
            _assign(self.status_label, "text", f"Done in {results.elapsed:.2g} s. Report: {self.report_path}")
            self.progress_bar.value = 100
        except (ValueError, OSError) as exc:
            import excitationReport
            self.report_output.value = excitationReport.report_text(results)
            await self._show_error("Excitation Output Error", str(exc))
        finally:
            self._set_buttons(True)
        if was_displayed:
            await self._display()  # a scene is on screen: show the new state

    # --- results ------------------------------------------------------------
    def summary_text(self) -> str:
        r, m = self.results, self.results.moments
        return (f"State {r.state.number}: D = {m.d_index * A:.2f} A, t = {m.t_index * A:.2f} A, Sr = {r.grid.sr:.2f}, "
                f"Delta sigma = {m.delta_sigma_index * A:.2f} A; Delta_r = {r.delta_r.value * A:.2f} A, "
                f"Lambda = {r.lambda_.value:.2f}. Confidence: {r.classification.confidence.split(':')[0]}.")

    def descriptor_rows(self) -> List[Tuple[str, str, str]]:
        r, m, g = self.results, self.results.moments, self.results.grid
        ev = xa.HARTREE_TO_EV
        rows = [
            ("Excitation energy", f"{r.state.energy_ev:.4f}", "eV"),
            ("D index (hole-electron centroid distance)", f"{m.d_index * A:.3f}", "A"),
            ("t index (separation)", f"{m.t_index * A:.3f}", "A"),
            ("Sr index (overlap)", f"{g.sr:.4f}", ""),
            ("Sm index (overlap)", f"{g.sm:.4f}", ""),
            ("Delta sigma (electron - hole RMSD)", f"{m.delta_sigma_index * A:.3f}", "A"),
            ("RMSD of hole / electron", f"{m.sigma_hole_norm * A:.3f} / {m.sigma_electron_norm * A:.3f}", "A"),
            ("H index", f"{m.h_index * A:.3f}", "A"),
            ("HDI / EDI (delocalization)", f"{g.hdi:.2f} / {g.edi:.2f}", ""),
            ("Delta_r", f"{r.delta_r.value * A:.3f}", "A"),
            ("Lambda", f"{r.lambda_.value:.4f}", ""),
            ("Ghost-hunter index", f"{r.ghost.value * ev:.3f}", "eV"),
            ("Transition dipole |mu|", f"{float((m.transition_dipole ** 2).sum()) ** 0.5:.4f}", "a.u."),
            ("Oscillator strength (calc / printed)",
             f"{m.oscillator(r.state.energy_ev):.4f} / " + (f"{r.state.oscillator:.4f}" if r.state.oscillator
                                                            is not None else "-"), ""),
            ("Coulomb attraction (C functions)", f"{r.coulomb.gaussian * ev:.3f}", "eV"),
        ]
        if r.coulomb.exact is not None:
            rows.append((f"Coulomb attraction (exact, {r.coulomb.mode})", f"{r.coulomb.exact * ev:.3f}", "eV"))
        if r.has_diffuse:
            rows.append(("Electron on diffuse functions", f"{100 * r.diffuse_electron:.1f}", "%"))
        rows.append(("Integral of hole / electron", f"{m.hole_integral:.5f} / {m.electron_integral:.5f}", ""))
        return rows

    def part_rows(self) -> List[Tuple[str, str, str, str, str]]:
        r = self.results
        if r.ifct is not None:
            f = r.ifct
            rows = []
            for x, h, e in zip(f.fragments, f.hole, f.electron):
                o = (h * e) ** 0.5
                rows.append((x.name, f"{100 * h:.2f}", f"{100 * e:.2f}", f"{100 * o:.2f}", f"{100 * (e - h):.2f}"))
            rows.append(("Intrinsic CT / LE", f"{f.intrinsic_ct:.2f}", f"{f.intrinsic_le:.2f}", "", ""))
            return rows
        hydrogens = bool(self.hydrogens_switch.value)
        wfn, atoms = r.wavefunction, r.atoms
        symbols = wfn.atom_symbols()
        return [(f"{k + 1}({symbols[k]})", f"{100 * atoms.hole[k]:.2f}", f"{100 * atoms.electron[k]:.2f}",
                 f"{100 * atoms.overlap[k]:.2f}", f"{100 * atoms.difference[k]:.2f}")
                for k in range(wfn.n_atoms) if hydrogens or int(wfn.atomic_numbers[k]) != 1]

    def _fill_results(self):
        self._syncing = True
        try:
            _assign(self.classification_label, "text", f"Classification: {self.results.classification.verdict}")
            _assign(self.summary_label, "text", self.summary_text())
            self.descriptor_table.data = self.descriptor_rows()
            self.part_table.data = self.part_rows()
            fit_table_columns(self.descriptor_table, DESCRIPTOR_FRACTIONS)
            fit_table_columns(self.part_table, PART_FRACTIONS)
        finally:
            self._syncing = False

    def on_hydrogens_change(self, widget, **kwargs):
        if self._syncing or self.results is None:
            return
        self.part_table.data = self.part_rows()

    # --- 3D display -------------------------------------------------------------
    def fchk_frame(self):
        wfn = self.results.wavefunction
        return [(s, tuple(float(c) for c in xyz)) for s, xyz in zip(wfn.atom_symbols(), wfn.coordinates * A)]

    async def _settle_structure(self):
        """Show the fchk molecule, asking first if the viewer shows a different structure."""
        if self.results is None or self._structure_settled_for is self.results.wavefunction:
            return
        frame = self.fchk_frame()
        current = self.viewer.get_current_frame_data() if self.viewer.frames else []
        if not current:
            self.viewer.load_designed_molecule(frame)
        elif not vd.same_structure(current, frame):
            if await self._ask("Show the molecule of this analysis?",
                               f"The 3D viewer shows another structure ({len(current)} atoms).\n\n"
                               f"Yes: replace it with the {len(frame)} atoms of {Path(self.results.fchk_path).name}.\n"
                               "No: keep it; the surfaces are drawn in the .fchk coordinates."):
                self.viewer.load_designed_molecule(frame)
        self._structure_settled_for = self.results.wavefunction

    def _scene_is_displayed(self) -> bool:
        layer = getattr(self.viewer, "excitation_layer", None)
        return layer is not None and bool(layer.items())

    def field_entries(self, spacing: float, isovalue: Optional[float]) -> Dict[str, CubeEntry]:
        """Every field as a Cube-tab entry (mesh cache), computed once per grid. Thread-safe."""
        box_isovalue = isovalue if isovalue is not None else min(xa.default_isovalue(f) for f in xa.FIELDS)
        key = (id(self.results), spacing, box_isovalue)
        if self._volumes is None or self._volumes[0] != key:
            volumes = xa.field_volumes(self.results, spacing, box_isovalue, use_numba=self.use_numba)
            self._volumes = (key, {name: CubeEntry.from_volume(volume.title, volume, saved=False)
                                   for name, volume in volumes.items()})
        entries = self._volumes[1]
        for name, entry in entries.items():
            entry.isovalue = isovalue if isovalue is not None else xa.default_isovalue(name)
        return entries

    def compute_scene_items(self):
        """(render items, notes) of the current Show, colours, opacity and centroids. Runs in a worker thread."""
        if self._volumes is None:
            return [], []
        entries = self._volumes[1]
        items, notes = [], []
        for name in SHOW_MODES[self.show]:
            entry = entries[name]
            for sign in ((1, -1) if name in xa.SIGNED_FIELDS else (1,)):
                mesh = entry.mesh(sign, self.use_numba)
                if mesh is None:
                    notes.append(f"No {'+' if sign > 0 else '−'}{entry.isovalue:g} surface for {name}.")
                    continue
                if mesh.clipped:
                    notes.append(f"The {name} surface reaches the edge of the grid; raise the isovalue.")
                items.append(surface_item(mesh, vd.parse_colour(self.colours[colour_key(name, sign)]), self.opacity))
        if self.show_centroids:
            m = self.results.moments
            items.append(sphere_item(m.hole_centroid, CENTROID_RADIUS, vd.parse_colour(self.colours[("hole", 1)])))
            items.append(sphere_item(m.electron_centroid, CENTROID_RADIUS,
                                     vd.parse_colour(self.colours[("electron", 1)])))
            if m.d_index >= MIN_ARROW_LENGTH:
                items.append(arrow_item(m.hole_centroid, m.electron_centroid, ARROW_RADIUS,
                                        vd.parse_colour("Orange")))
        return items, notes

    async def refresh_scene(self):
        """Rebuild the drawn surfaces; changes arriving meanwhile are coalesced into one rerun."""
        if self._refreshing:
            self._refresh_pending = True
            return
        self._refreshing = True
        try:
            while True:
                self._refresh_pending = False
                items, notes = await asyncio.to_thread(self.compute_scene_items)
                if items:
                    # Replace another tab's result (an opaque ESP surface would hide these).
                    claim = getattr(self.viewer, "show_analysis_layer", None)
                    if claim is not None:
                        claim("excitation_layer")
                self.viewer.excitation_layer.set_items(items)
                self.viewer.invalidate_scene_cache()
                self.viewer.request_projection_update()
                if not self._refresh_pending:
                    break
        finally:
            self._refreshing = False
        _assign(self.status_label, "text", " ".join([f"{self.show} displayed."] + notes))

    async def _display(self) -> bool:
        try:
            params = self._collect_display_inputs()
            colours = self._collect_colours()
        except ValueError as exc:
            await self._show_error("Display Input Error", str(exc))
            return False
        self.colours = colours
        self.show = self.show_selection.value or DEFAULT_SHOW
        _assign(self.status_label, "text", "Computing hole, electron and transition density on the grid...")
        try:
            await asyncio.to_thread(self.field_entries, **params)
        except (ValueError, MemoryError) as exc:
            await self._show_error("Excitation Grid Error", str(exc))
            return False
        await self._settle_structure()
        await self.refresh_scene()
        return True

    async def display_excitation(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Run the analysis first.")
            return
        if not await self._display():
            return
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)

    # Synchronous, guarded handlers (see the Cube tab: an async handler would run
    # after _syncing is cleared and replay the tab's own writes as edits).
    def on_grid_confirm(self, widget, **kwargs):
        if self._syncing or self.results is None or not self._scene_is_displayed():
            return
        self._start(self._display())

    def on_show_change(self, widget, **kwargs):
        if self._syncing:
            return
        self.show = widget.value or DEFAULT_SHOW
        if self._scene_is_displayed():
            self._start(self.refresh_scene())

    def on_opacity_change(self, widget, **kwargs):
        if self._syncing:
            return
        self.opacity = max(MIN_OPACITY, min(1.0, float(widget.value)))
        if self._scene_is_displayed():
            self._start(self.refresh_scene())

    def on_centroids_change(self, widget, **kwargs):
        if self._syncing:
            return
        self.show_centroids = bool(widget.value)
        if self._scene_is_displayed():
            self._start(self.refresh_scene())

    def on_colour_confirm(self, widget, **kwargs):
        if self._syncing:
            return
        key = next(k for k, w in self.colour_inputs.items() if w is widget)
        try:
            self.colours[key] = canonical_colour(widget.value)
        except ValueError as exc:
            self._syncing = True
            try:
                _assign(widget, "value", self.colours[key])  # put the colour in use back
            finally:
                self._syncing = False
            self._start(self._show_error("Colour", f"{COLOUR_LABELS[key]}: {exc}"))
            return
        if self._scene_is_displayed():
            self._start(self.refresh_scene())

    # --- export -------------------------------------------------------------------
    async def export_cubes(self, widget):
        """The quantities of the current Show as cube files in the output folder."""
        if self.results is None:
            await self._show_info("No Results", "Run the analysis first.")
            return
        try:
            params = self._collect_display_inputs()
        except ValueError as exc:
            await self._show_error("Export Input Error", str(exc))
            return
        fields = SHOW_MODES[self.show_selection.value or DEFAULT_SHOW]
        isovalue = params["isovalue"] if params["isovalue"] is not None else min(
            xa.default_isovalue(f) for f in xa.FIELDS)
        self._set_buttons(False)
        self.progress_bar.value = 0
        try:
            written = await asyncio.to_thread(
                self.analyser.write_cubes, self.results, self.output_dir(), fields, params["spacing"], isovalue,
                self._progress_callback("Writing the cube files..."), self.use_numba)
        except (ValueError, MemoryError) as exc:
            await self._show_error("Export Cubes Error", str(exc))
            return
        finally:
            self._set_buttons(True)
        self.progress_bar.value = 100
        _assign(self.status_label, "text", f"{len(written)} cube files written to {Path(written[0]).parent}.")

    async def send_to_cube_tab(self, widget):
        """Add the displayed quantities to the Cube tab (slices, arithmetic, saving)."""
        cube = getattr(self.viewer, "cube_controller", None)
        if self.results is None or self._volumes is None or cube is None:
            await self._show_info("Nothing to Send", "Display the analysis first; the shown quantities are sent.")
            return
        names = []
        for field in SHOW_MODES[self.show]:
            entry = self._volumes[1][field]
            name = xa.cube_filename(self.results.stem, self.results.state.number, field)[:-len(".cube")]
            added = cube.add_result(entry.volume, name)
            added.visible = False  # already drawn here; tick Show in the Cube tab to draw it there too
            names.append(added.name)
        cube.sync_widgets(selected=names[0])
        _assign(self.status_label, "text", f"Sent to the Cube tab (Show unticked): {', '.join(names)}.")

    async def plot_heat_map(self, widget):
        """The chosen composition / charge-transfer matrix as a heat map in the plot viewer."""
        from displayPlots import launch_plot_viewer  # pins Agg: import only when a map is plotted

        if self.results is None:
            await self._show_info("No Results", "Run the analysis first.")
            return
        labels = {label: kind for kind, label in xa.HEATMAP_KINDS.items()}
        kind = labels.get(self.heatmap_selection.value, "atoms")
        try:
            path = xa.ExcitationAnalyser.write_heat_map(self.results, kind, self.output_dir(),
                                                        bool(self.hydrogens_switch.value))
        except ValueError as exc:
            await self._show_error("Heat Map Error", str(exc))
            return
        _assign(self.status_label, "text", f"Heat map: {path}")
        if not launch_plot_viewer(path):
            await self._show_info("Heat Map Saved", f"The plot viewer could not be started. The heat map was saved "
                                                    f"as\n{path}\nand can be opened with Plotter (JSON plot type).")

    async def export_3d_image(self, widget):
        """Save the 3D canvas exactly as shown: captured first, then the file is asked for."""
        from molecularViewer import save_canvas_image  # molecularViewer imports this module

        if not self.viewer.viewer_window_is_open():
            await self._show_info("No Image", "Open the 3D view with Display first; the image is taken from "
                                              "that window as shown.")
            return
        try:
            image = await self.viewer.capture_frame()
        except RuntimeError as exc:
            await self._show_error("Export Image Error", str(exc))
            return
        stem = "excitation"
        if self.results is not None:
            stem = f"{self.results.stem}_state{self.results.state.number}_excitation"
        try:
            path = await self._ask_save_path("Export the 3D view", f"{stem}.png", IMAGE_FILE_TYPES)
            if not path:
                return
            saved = save_canvas_image(image, str(path))
        except (ValueError, OSError) as exc:
            await self._show_error("Export Image Error", str(exc))
            return
        height, width = image.shape[:2]
        _assign(self.status_label, "text", f"Image saved ({width} x {height} px): {saved}")

    def clear_output(self, widget):
        self.results = None
        self.report_path = None
        self._volumes = None
        self._structure_settled_for = None
        self._syncing = True
        try:
            self.descriptor_table.data = []
            self.part_table.data = []
            _assign(self.classification_label, "text", "No results.")
            _assign(self.summary_label, "text", "")
        finally:
            self._syncing = False
        self.report_output.value = HelpGqteaWin.help_excitation
        self.viewer.excitation_layer.clear()
        self.viewer.invalidate_scene_cache()
        self.progress_bar.value = 0
        _assign(self.status_label, "text", "Ready.")
