"""The 3D Molecular Viewer's NTO tab (``NTOTabController``).

Natural transition orbitals of one excited state from a Gaussian .fchk and
the CIS/TDA/TD output of the same job. The numerical work is
``ntoAnalysis.NTOAnalyser`` (headless). The hole and particle isosurfaces are
drawn in the viewer's own 3D window through ``viewer.nto_layer`` (a second
``cubeRender.CubeSceneLayer``); meshes are cached per orbital, isovalue and
sign by the Cube tab's ``CubeEntry``, so colours, opacity and Show never
re-mesh. Widget writes go through ``cubeTab._assign`` (every write re-lays out
the viewer window) and the change handlers are synchronous and
``_syncing``-guarded, as in the Cube and QTAIM tabs.
"""

from __future__ import annotations

import asyncio
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import excitationData as ed
import ntoAnalysis as na
import uvVisSpectrum as uv
import volumeData as vd
from cubeRender import surface_item
from cubeTab import CubeEntry, _assign, _item_names, canonical_colour
from help import HelpGqteaWin
from qtaimTab import fit_table_columns, parse_optional_float

SHOW_HOLE, SHOW_PARTICLE, SHOW_BOTH = "Hole", "Particle", "Hole and particle"
SHOW_MODES = {SHOW_HOLE: ("hole",), SHOW_PARTICLE: ("particle",), SHOW_BOTH: ("hole", "particle")}
COLOUR_KEYS = (("hole", 1), ("hole", -1), ("particle", 1), ("particle", -1))
DEFAULT_COLOURS = {("hole", 1): "Blue", ("hole", -1): "Red", ("particle", 1): "Green", ("particle", -1): "Yellow"}
COLOUR_LABELS = {("hole", 1): "Hole +", ("hole", -1): "Hole −", ("particle", 1): "Particle +",
                 ("particle", -1): "Particle −"}
TABLE_HEADINGS = ["Pair", "Spin", "Weight", "%", "Cumulative %"]
TABLE_ACCESSORS = ["pair", "spin", "weight", "percent", "cumulative"]
TABLE_COLUMN_FRACTIONS = (0.12, 0.16, 0.24, 0.20, 0.26)
IMAGE_FILE_TYPES = ["png", "tif", "tiff", "bmp"]
MIN_OPACITY = 0.1


def pair_label(results: na.NTOResults, pair: na.NTOPair) -> str:
    spin = "" if results.restricted else f" {pair.spin}"
    return f"{pair.pair}{spin} ({100 * pair.weight:.2f} %)"


class NTOTabController:
    """The 3D Molecular Viewer's NTO tab."""

    LABEL_WIDTH = 210
    HINT_COLOR = "#666666"

    def __init__(self, viewer):
        self.viewer = viewer
        self.analyser = na.NTOAnalyser()
        self.loaded = None                       # (wavefunction, ExcitationFile, notes) of loaded_paths
        self.loaded_paths: Optional[Tuple[str, str]] = None
        self.results: Optional[na.NTOResults] = None
        self.report_path: Optional[str] = None
        self.show = SHOW_BOTH
        self.colours: Dict[Tuple[str, int], str] = dict(DEFAULT_COLOURS)
        self.opacity = 1.0
        self.use_numba: Optional[bool] = None    # None = if installed
        self._state_items: Dict[str, int] = {}   # state Selection label -> state number
        self._pair_items: Dict[str, na.NTOPair] = {}
        self._orbitals: Optional[Tuple[tuple, Dict[str, CubeEntry]]] = None  # (key, {"hole", "particle"})
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
            self.run_button = toga.Button("Run", on_press=self.workflow, style=Pack(width=110, margin=5))
            self.plot_spectrum_button = toga.Button("Plot Spectrum", on_press=self.plot_spectrum,
                                                    style=Pack(width=110, margin=5))
            self.hwhm_input = toga.TextInput(placeholder=f"{uv.DEFAULT_HWHM_EV:g}", style=Pack(width=120))
            self.clear_button = toga.Button("Clear", on_press=self.clear_output, style=Pack(width=110, margin=5))
            self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
            self.progress_bar = toga.ProgressBar(max=100, style=Pack(margin=(0, 5, 5, 5)))
            self.summary_label = toga.Label("No results.", style=Pack(margin=(4, 0, 4, 0)))
            self.pair_table = toga.Table(headings=TABLE_HEADINGS, accessors=TABLE_ACCESSORS, data=[],
                                         on_select=self.on_pair_select, style=Pack(height=170))

            self.pair_selection = toga.Selection(items=[], on_change=self.on_pair_change, style=Pack(width=200))
            self.show_selection = toga.Selection(items=list(SHOW_MODES), on_change=self.on_show_change,
                                                 style=Pack(width=200))
            self.show_selection.value = self.show
            self.isovalue_input = toga.TextInput(placeholder=f"{na.ORBITAL_ISOVALUE:g}",
                                                 on_confirm=self.on_grid_confirm, style=Pack(width=120))
            self.spacing_input = toga.TextInput(placeholder=f"{na.DEFAULT_SPACING_BOHR:g}",
                                                on_confirm=self.on_grid_confirm, style=Pack(width=120))
            self.colour_inputs = {key: toga.TextInput(value=DEFAULT_COLOURS[key], on_confirm=self.on_colour_confirm,
                                                      style=Pack(width=120)) for key in COLOUR_KEYS}
            self.opacity_slider = toga.Slider(min=MIN_OPACITY, max=1.0, value=1.0, on_change=self.on_opacity_change,
                                              style=Pack(width=200))
            self.display_button = toga.Button("Display", on_press=self.display_nto, style=Pack(width=110, margin=5))

            self.all_pairs_switch = toga.Switch("All listed pairs", value=False, style=Pack(margin=(8, 10, 0, 5)))
            self.export_cubes_button = toga.Button("Export cubes", on_press=self.export_cubes,
                                                   style=Pack(width=140, margin=5))
            self.export_fchk_button = toga.Button("Export NTO .fchk", on_press=self.export_fchk,
                                                  style=Pack(width=140, margin=5))
            self.send_button = toga.Button("Send to Cube tab", on_press=self.send_to_cube_tab,
                                           style=Pack(width=140, margin=5))
            self.export_3d_button = toga.Button("Export 3D image...", on_press=self.export_3d_image,
                                                style=Pack(width=150, margin=5))
            self.report_output = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_nto,
                                                         style=Pack(height=300, font_family="monospace"))
        finally:
            self._syncing = False

        colour_rows = [self._row(self._label(f"{COLOUR_LABELS[a]} / {COLOUR_LABELS[b]}:"),
                                 self.colour_inputs[a], toga.Label("  ", style=Pack(width=10)), self.colour_inputs[b])
                       for a, b in ((COLOUR_KEYS[0], COLOUR_KEYS[1]), (COLOUR_KEYS[2], COLOUR_KEYS[3]))]
        children = [
            self._hint("Natural transition orbitals of a CIS/TDA/TD state; hole and particle are drawn in this "
                       "viewer's 3D window."),
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
            self._heading("Excited state"),
            self._row(self._label("State:"), self.state_selection),
            self._row(self._label("Peak half-width (eV):"), self.hwhm_input),
            self._hint("Plot Spectrum draws every listed state (lambda, f) and the broadened absorption, "
                       "as GaussView does."),
            self._row(self.run_button, self.plot_spectrum_button, self.clear_button),
            self.status_label,
            self.progress_bar,
            self._heading("NTO pairs"),
            self.summary_label,
            self.pair_table,
            self._heading("Display"),
            self._row(self._label("Pair:"), self.pair_selection),
            self._row(self._label("Show:"), self.show_selection),
            self._row(self._label("Isovalue (a.u.):"), self.isovalue_input),
            self._row(self._label("Grid spacing (bohr):"), self.spacing_input),
            *colour_rows,
            self._hint("Colours: Blue, Red, Green, Yellow, Orange, Cyan, Magenta, Purple, White, Grey or #rrggbb "
                       "(Enter applies)."),
            self._row(self._label("Opacity:"), self.opacity_slider),
            self._row(self.display_button),
            self._heading("Export"),
            self._row(self.export_cubes_button, self.all_pairs_switch),
            self._row(self.export_fchk_button, self.send_button, self.export_3d_button),
            self._hint("Cubes and the NTO .fchk are written to the output folder; their names give the state, "
                       "pair, spin and hole/particle."),
            self.report_output,
        ]
        self._attach_column_fitting()
        return toga.Box(children=children, style=Pack(direction=COLUMN, margin=10))

    def _attach_column_fitting(self):
        native = getattr(getattr(self.pair_table, "_impl", None), "native", None)
        if native is None:
            return
        try:
            native.Resize += lambda sender, event: fit_table_columns(self.pair_table, TABLE_COLUMN_FRACTIONS)
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
        fchk, excitation = self._paths()
        if self.loaded_paths != (fchk, excitation):
            raise ValueError("Press Load to read these files and list their excited states.")
        number = self._state_items.get(self.state_selection.value)
        if number is None:
            raise ValueError("Choose an excited state.")
        return {"fchk_path": fchk, "excitation_path": excitation, "state_number": number}

    def _collect_display_inputs(self) -> Dict[str, object]:
        if self.results is None:
            raise ValueError("Run the analysis first.")
        pair = self._pair_items.get(self.pair_selection.value)
        if pair is None:
            raise ValueError("Choose an NTO pair.")
        return {
            "spin": pair.spin,
            "pair": pair.pair,
            "spacing": parse_optional_float(self.spacing_input.value, na.DEFAULT_SPACING_BOHR, "The grid spacing"),
            "isovalue": parse_optional_float(self.isovalue_input.value, na.ORBITAL_ISOVALUE, "The isovalue"),
        }

    def _collect_colours(self) -> Dict[Tuple[str, int], str]:
        colours = {}
        for key, widget in self.colour_inputs.items():
            try:
                colours[key] = canonical_colour(widget.value)
            except ValueError as exc:
                raise ValueError(f"{COLOUR_LABELS[key]}: {exc}") from None
        return colours

    def output_dir(self) -> str:
        return na.resolve_output_dir(self.results.fchk_path, self.output_dir_input.value)

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
        for b in (self.load_button, self.run_button, self.plot_spectrum_button, self.clear_button, self.display_button,
                  self.export_cubes_button, self.export_fchk_button, self.send_button):
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
        message = f"{len(data.states)} excited states read ({kind}, {wfn.n_mo} MOs). Choose one and press Run."
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
            await self._show_error("NTO Input Error", str(exc))
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
            await self._show_error("NTO Calculation Error", str(exc))
            return
        was_displayed = self._scene_is_displayed()
        self.results = results  # kept even if the report cannot be written
        self._orbitals = None
        self._fill_results()
        try:
            self.report_path = self.analyser.write_report(results, self.output_dir())
            self.report_output.value = Path(self.report_path).read_text(encoding="utf-8")
            _assign(self.status_label, "text", f"Done in {results.elapsed:.2g} s. Report: {self.report_path}")
            self.progress_bar.value = 100
        except (ValueError, OSError) as exc:
            self.report_output.value = na.report_text(results)
            await self._show_error("NTO Output Error", str(exc))
        finally:
            self._set_buttons(True)
        if was_displayed:
            await self._display()  # a scene is on screen: show the new state's first pair

    async def plot_spectrum(self, widget):
        """Every listed excited state as a UV-Vis spectrum, saved as <stem>_UVVis.json and opened in
        the plot viewer. Independent of the NTO run: it reads only the loaded excitation file."""
        from displayPlots import launch_plot_viewer  # pins Agg: import only when a spectrum is plotted

        try:
            paths = self._paths()
            hwhm = parse_optional_float(self.hwhm_input.value, uv.DEFAULT_HWHM_EV, "The peak half-width")
        except ValueError as exc:
            await self._show_error("Spectrum Input Error", str(exc))
            return
        if self.loaded_paths != paths and not await self._load(paths):
            return
        fchk, excitation = paths
        folder = na.resolve_output_dir(fchk, self.output_dir_input.value)
        source = {"fchk": Path(fchk).name, "excitation_file": Path(excitation).name}
        try:
            entry = uv.spectrum_entry(self.loaded[1].states, hwhm, source=source)
            path = uv.write_spectrum(entry, Path(folder) / uv.spectrum_filename(Path(fchk).stem))
        except ValueError as exc:
            await self._show_error("Spectrum Error", str(exc))
            return
        _assign(self.status_label, "text", f"Spectrum of {len(entry['states'])} states: {path}")
        if not launch_plot_viewer(path):
            await self._show_info("Spectrum Saved", f"The plot viewer could not be started. The spectrum was "
                                                    f"saved as\n{path}\nand can be opened with Plotter "
                                                    "(JSON plot type).")

    def summary_text(self) -> str:
        r = self.results
        top = r.listed_pairs()[0]
        spin = "" if r.restricted else f" {top.spin}"
        return (f"State {r.state.number}: largest pair {top.pair}{spin} carries {100 * top.weight:.2f} %; "
                f"weights sum to {r.weight_sum():.4f}; normalisation deviation {r.deviation:+.6f}.")

    def table_rows(self):
        rows, cumulative = [], 0.0
        for p in self.results.listed_pairs():
            cumulative += p.weight
            rows.append((str(p.pair), p.spin, f"{p.weight:.6f}", f"{100 * p.weight:.2f}", f"{100 * cumulative:.2f}"))
        return rows

    def _fill_results(self):
        pairs = self.results.listed_pairs()
        self._pair_items = {pair_label(self.results, p): p for p in pairs}
        self._syncing = True
        try:
            _assign(self.summary_label, "text", self.summary_text())
            self.pair_table.data = self.table_rows()
            fit_table_columns(self.pair_table, TABLE_COLUMN_FRACTIONS)
            if list(self._pair_items) != _item_names(self.pair_selection):
                self.pair_selection.items = list(self._pair_items)
            _assign(self.pair_selection, "value", next(iter(self._pair_items)))
        finally:
            self._syncing = False

    # --- 3D display -------------------------------------------------------------
    def fchk_frame(self):
        """The molecule as a viewer frame: [(symbol, (x, y, z) Angstrom), ...]."""
        wfn = self.results.wavefunction
        return [(s, tuple(float(c) for c in xyz))
                for s, xyz in zip(wfn.atom_symbols(), wfn.coordinates * na.BOHR_RADIUS_ANGSTROM)]

    async def _settle_structure(self):
        """Show the fchk molecule, asking first if the viewer shows a different structure."""
        if self.results is None or self._structure_settled_for is self.results.wavefunction:
            return
        frame = self.fchk_frame()
        # A fresh viewer keeps a 5-atom placeholder in molecule_data with no frames: never worth a question.
        current = self.viewer.get_current_frame_data() if self.viewer.frames else []
        if not current:
            self.viewer.load_designed_molecule(frame)
        elif not vd.same_structure(current, frame):
            if await self._ask("Show the NTO molecule?",
                               f"The 3D viewer shows another structure ({len(current)} atoms).\n\n"
                               f"Yes: replace it with the {len(frame)} atoms of {Path(self.results.fchk_path).name}.\n"
                               "No: keep it; the orbitals are drawn in the .fchk coordinates."):
                self.viewer.load_designed_molecule(frame)
        self._structure_settled_for = self.results.wavefunction

    def _scene_is_displayed(self) -> bool:
        layer = getattr(self.viewer, "nto_layer", None)
        return layer is not None and bool(layer.items())

    def orbital_entries(self, spin: str, pair: int, spacing: float, isovalue: float) -> Dict[str, CubeEntry]:
        """The pair's hole and particle as Cube-tab entries (mesh cache), computed once per grid. Thread-safe."""
        key = (id(self.results), spin, pair, spacing, isovalue)
        if self._orbitals is not None and self._orbitals[0] == key:
            return self._orbitals[1]
        volumes = na.pair_volumes(self.results, spin, pair, spacing, isovalue, use_numba=self.use_numba)
        entries = {}
        for kind, volume in zip(na.KINDS, volumes):
            entry = CubeEntry.from_volume(volume.title, volume, saved=False)
            entry.isovalue = isovalue
            entries[kind] = entry
        self._orbitals = (key, entries)
        return entries

    def compute_scene_items(self):
        """(render items, notes) of the current pair, Show, colours and opacity. Runs in a worker thread."""
        if self._orbitals is None:
            return [], []
        entries = self._orbitals[1]
        items, notes = [], []
        for kind in SHOW_MODES[self.show]:
            entry = entries[kind]
            for sign in (1, -1):
                mesh = entry.mesh(sign, self.use_numba)
                if mesh is None:
                    notes.append(f"No {'+' if sign > 0 else '−'}{entry.isovalue:g} surface for the {kind}.")
                    continue
                if mesh.clipped:
                    notes.append(f"The {kind} surface reaches the edge of the grid; raise the isovalue.")
                items.append(surface_item(mesh, vd.parse_colour(self.colours[(kind, sign)]), self.opacity))
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
                self.viewer.nto_layer.set_items(items)
                self.viewer.invalidate_scene_cache()
                self.viewer.request_projection_update()
                if not self._refresh_pending:
                    break
        finally:
            self._refreshing = False
        shown = sum(1 for i in items if i.kind == "surface")
        _assign(self.status_label, "text", " ".join([f"{shown} isosurface{'s' if shown != 1 else ''} shown."]
                                                    + notes))

    async def _display(self) -> bool:
        try:
            params = self._collect_display_inputs()
            colours = self._collect_colours()
        except ValueError as exc:
            await self._show_error("Display Input Error", str(exc))
            return False
        self.colours = colours
        self.show = self.show_selection.value or SHOW_BOTH
        _assign(self.status_label, "text", "Computing the orbitals on the grid...")
        try:
            await asyncio.to_thread(self.orbital_entries, **params)
        except (ValueError, MemoryError) as exc:
            await self._show_error("NTO Grid Error", str(exc))
            return False
        await self._settle_structure()
        await self.refresh_scene()
        return True

    async def display_nto(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Run the analysis first.")
            return
        if not await self._display():
            return
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)

    # Synchronous, guarded handlers: see the Cube tab (an async handler would
    # run after _syncing is cleared and replay the tab's own writes as edits).
    def on_pair_select(self, widget, **kwargs):
        if self._syncing or self.results is None:
            return
        row = getattr(widget, "selection", None)
        if row is None:
            return
        wanted = (getattr(row, "spin", None), int(getattr(row, "pair", 0)))
        label = next((l for l, p in self._pair_items.items() if (p.spin, p.pair) == wanted), None)
        if label is None or label == self.pair_selection.value:
            return
        self._syncing = True
        try:
            _assign(self.pair_selection, "value", label)
        finally:
            self._syncing = False
        if self._scene_is_displayed():
            self._start(self._display())

    def on_pair_change(self, widget, **kwargs):
        if self._syncing or self.results is None or not self._scene_is_displayed():
            return
        self._start(self._display())

    def on_grid_confirm(self, widget, **kwargs):
        if self._syncing or self.results is None or not self._scene_is_displayed():
            return
        self._start(self._display())

    def on_show_change(self, widget, **kwargs):
        if self._syncing:
            return
        self.show = widget.value or SHOW_BOTH
        if self._scene_is_displayed():
            self._start(self.refresh_scene())

    def on_opacity_change(self, widget, **kwargs):
        if self._syncing:
            return
        self.opacity = max(MIN_OPACITY, min(1.0, float(widget.value)))
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
        if self.results is None:
            await self._show_info("No Results", "Run the analysis first.")
            return
        try:
            params = self._collect_display_inputs()
        except ValueError as exc:
            await self._show_error("Export Input Error", str(exc))
            return
        if self.all_pairs_switch.value:
            pairs = self.results.listed_pairs()
        else:
            pairs = [self._pair_items[self.pair_selection.value]]
        self._set_buttons(False)
        self.progress_bar.value = 0
        try:
            written = await asyncio.to_thread(
                self.analyser.write_cubes, self.results, self.output_dir(), pairs, params["spacing"],
                params["isovalue"], self._progress_callback("Writing the cube files..."), self.use_numba)
        except (ValueError, MemoryError) as exc:
            await self._show_error("Export Cubes Error", str(exc))
            return
        finally:
            self._set_buttons(True)
        self.progress_bar.value = 100
        _assign(self.status_label, "text", f"{len(written)} cube files written to {Path(written[0]).parent}.")

    async def export_fchk(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Run the analysis first.")
            return
        path = Path(self.output_dir()) / na.nto_fchk_filename(self.results.stem, self.results.state.number)
        try:
            written = await asyncio.to_thread(self.analyser.write_nto_fchk, self.results, str(path))
        except ValueError as exc:
            await self._show_error("Export NTO .fchk Error", str(exc))
            return
        _assign(self.status_label, "text", f"NTO .fchk written: {written} (orbital energies = NTO weights).")

    async def send_to_cube_tab(self, widget):
        """Add the displayed pair's hole and particle to the Cube tab (slices, arithmetic, saving)."""
        cube = getattr(self.viewer, "cube_controller", None)
        if self.results is None or self._orbitals is None or cube is None:
            await self._show_info("Nothing to Send", "Display a pair first; its two orbitals are sent.")
            return
        _, spin, pair, _, _ = self._orbitals[0]
        names = []
        for kind, entry in self._orbitals[1].items():
            name = na.cube_filename(self.results.stem, self.results.state.number, pair, kind,
                                    None if self.results.restricted else spin)[:-len(".cube")]
            added = cube.add_result(entry.volume, name)
            added.visible = False  # already drawn here; tick Show in the Cube tab to draw it there too
            names.append(added.name)
        cube.sync_widgets(selected=names[0])
        _assign(self.status_label, "text", f"Sent to the Cube tab (Show unticked): {', '.join(names)}.")

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
        stem = "NTO"
        if self.results is not None and self._orbitals is not None:
            _, spin, pair, _, _ = self._orbitals[0]
            spin_part = "" if self.results.restricted else f"_{spin}"
            stem = f"{self.results.stem}_state{self.results.state.number}_NTO{pair}{spin_part}"
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
        self._orbitals = None
        self._pair_items = {}
        self._structure_settled_for = None
        self._syncing = True
        try:
            self.pair_table.data = []
            self.pair_selection.items = []
            _assign(self.summary_label, "text", "No results.")
        finally:
            self._syncing = False
        self.report_output.value = HelpGqteaWin.help_nto
        self.viewer.nto_layer.clear()
        self.viewer.invalidate_scene_cache()
        self.progress_bar.value = 0
        _assign(self.status_label, "text", "Ready.")
