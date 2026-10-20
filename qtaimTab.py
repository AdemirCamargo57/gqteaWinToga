"""The 3D Molecular Viewer's QTAIM tab (``QTAIMTabController``).

It owns the run, the results and the widgets. The critical points and paths
are drawn in the viewer's own 3D window through ``viewer.qtaim_layer`` (a
``qtaimLayer.QTAIMSceneLayer``); the 2D plane map opens in the interactive
plot viewer from a ``{"type": "plane_map"}`` manifest written beside the
report. The numerical work is ``qtaimAnalysis.QTAIMAnalyser`` (headless).
"""

from __future__ import annotations

import asyncio
import json
import math
import os
from pathlib import Path
from typing import Dict, Optional, Tuple

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import qtaimAnalysis as qa
import qtaimTopology as qt
from cubeTab import _assign  # write-only-on-change: every widget write re-lays out the viewer
from help import HelpGqteaWin
from qtaimLayer import scene_from_results
from qtaimTopology import BOHR_RADIUS_ANGSTROM, CP_NAMES, CP_TYPES
from volumeData import same_structure  # shared with the Cube and ESP tabs

PLANE_THREE_ATOMS = "Three atoms"
PLANE_MODES = (PLANE_THREE_ATOMS, "XY", "XZ", "YZ")

# (field of qtaimTopology.CPSearchParameters / PathParameters, label); blank field = Multiwfn default
SEARCH_FIELDS = (
    ("max_iterations", "Maximum iterations:"),
    ("step_scale", "Step-size scale factor:"),
    ("gradient_tolerance", "Gradient-norm convergence:"),
    ("displacement_tolerance", "Displacement convergence:"),
    ("min_cp_distance", "Minimum distance between CPs (bohr):"),
    ("vdw_cutoff_scale", "Skip atoms beyond vdW sum x:"),
    ("singular_determinant", "Singular Hessian if |det| below:"),
    ("trust_radius", "Trust radius (bohr):"),
)
PATH_FIELDS = (
    ("max_points", "Maximum points per path:"),
    ("step", "Path step (bohr):"),
    ("stop_distance", "Stop within this of a CP (bohr):"),
)
INTEGER_FIELDS = {"max_iterations", "max_points"}
TABLE_HEADINGS = ["Row", "Type", "Assignment", "rho (a.u.)", "Laplacian", "Ellipticity"]
TABLE_ACCESSORS = ["number", "kind", "assignment", "rho", "laplacian", "ellipticity"]


def parse_optional_float(text: Optional[str], default: Optional[float], name: str,
                         positive: bool = True) -> Optional[float]:
    """Blank -> ``default``; otherwise a finite number (and > 0 when ``positive``)."""
    text = (text or "").strip()
    if not text:
        return default
    try:
        value = float(text.replace(",", "."))
    except ValueError:
        raise ValueError(f"{name} must be a number, not '{text}'.") from None
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} must be a {'positive ' if positive else ''}finite number, not {text}.")
    return value


def parse_optional_int(text: Optional[str], default: int, name: str, minimum: int = 1) -> int:
    text = (text or "").strip()
    if not text:
        return default
    try:
        value = int(text)
    except ValueError:
        raise ValueError(f"{name} must be a whole number, not '{text}'.") from None
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}, not {value}.")
    return value


def parse_atom_triple(text: str, n_atoms: int) -> Tuple[int, int, int]:
    """'1, 2, 3' or '1 2 3' (1-based) -> three different 0-based atom indices."""
    parts = (text or "").replace(",", " ").split()
    try:
        atoms = [int(p) for p in parts]
    except ValueError:
        raise ValueError(f"Plane atoms must be three atom numbers such as 1, 2, 3, not '{text}'.") from None
    if len(atoms) != 3 or len(set(atoms)) != 3:
        raise ValueError(f"Plane atoms must be three different atom numbers, not '{text}'.")
    bad = [a for a in atoms if not 1 <= a <= n_atoms]
    if bad:
        raise ValueError(f"Atom number(s) {', '.join(map(str, bad))} outside 1-{n_atoms}.")
    return tuple(a - 1 for a in atoms)


class QTAIMTabController:
    """The 3D Molecular Viewer's QTAIM tab."""

    LABEL_WIDTH = 250
    HINT_COLOR = "#666666"

    def __init__(self, viewer):
        self.viewer = viewer
        self.analyser = qa.QTAIMAnalyser()
        self.results: Optional[qa.QTAIMResults] = None
        self.report_path: Optional[str] = None
        self._syncing = False
        self._structure_settled_for: Optional[qa.QTAIMResults] = None
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

    @staticmethod
    def _default_text(params, name) -> str:
        value = getattr(params, name)
        return "none" if value is None else f"{value:g}"

    def build_tab(self) -> toga.Box:
        self.fchk_input = toga.TextInput(placeholder="Gaussian formatted checkpoint (.fchk / .fch)", style=Pack(flex=1))
        self.output_dir_input = toga.TextInput(placeholder="Blank = folder of the .fchk file", style=Pack(flex=1))
        self.mode_switches = {m: toga.Switch(f"{qt.MODE_LABELS[m]} (mode {m})", value=True,
                                             style=Pack(margin_right=12)) for m in qt.MODE_LABELS}
        search, paths = qt.CPSearchParameters(), qt.PathParameters()
        self.search_inputs = {name: toga.TextInput(placeholder=self._default_text(search, name), style=Pack(width=140))
                              for name, _ in SEARCH_FIELDS}
        self.path_inputs = {name: toga.TextInput(placeholder=self._default_text(paths, name), style=Pack(width=140))
                            for name, _ in PATH_FIELDS}

        self.run_button = toga.Button("Run", on_press=self.workflow, style=Pack(width=110, margin=5))
        self.display_button = toga.Button("Display", on_press=self.display_results, style=Pack(width=110, margin=5))
        self.clear_button = toga.Button("Clear", on_press=self.clear_output, style=Pack(width=110, margin=5))
        self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
        self.progress_bar = toga.ProgressBar(max=100, style=Pack(margin=(0, 5, 5, 5)))
        self.summary_label = toga.Label("No results.", style=Pack(margin=(4, 0, 4, 0)))
        self.cp_table = toga.Table(headings=TABLE_HEADINGS, accessors=TABLE_ACCESSORS, data=[],
                                   on_select=self.on_cp_select, style=Pack(height=220))

        self.type_switches = {t: toga.Switch(f"{t} {CP_NAMES[t]}", value=True, on_change=self.on_display_change,
                                             style=Pack(margin_right=10)) for t in CP_TYPES}
        self.bond_path_switch = toga.Switch("Bond paths", value=True, on_change=self.on_display_change,
                                            style=Pack(margin_right=10))
        self.ring_path_switch = toga.Switch("Ring paths", value=True, on_change=self.on_display_change,
                                            style=Pack(margin_right=10))
        self.label_switch = toga.Switch("CP numbers", value=True, on_change=self.on_display_change)

        self.plane_selection = toga.Selection(items=list(PLANE_MODES), style=Pack(width=160))
        self.plane_atoms_input = toga.TextInput(placeholder="e.g. 1, 2, 3 (blank = first three non-collinear)",
                                                style=Pack(flex=1))
        self.plane_offset_input = toga.TextInput(placeholder="0", style=Pack(width=120))
        self.margin_input = toga.TextInput(placeholder=f"{qa.PlaneMapOptions().margin_angstrom:g}",
                                           style=Pack(width=120))
        self.grid_points_input = toga.TextInput(placeholder=f"{qa.PlaneMapOptions().points}", style=Pack(width=120))
        self.contour_switch = toga.Switch("Density contours", value=True, style=Pack(margin_right=10))
        self.gradient_switch = toga.Switch("Gradient lines", value=True, style=Pack(margin_right=10))
        self.interbasin_switch = toga.Switch("Interbasin paths", value=True)
        self.map_button = toga.Button("Plot 2D map", on_press=self.plot_map, style=Pack(width=130, margin=5))

        self.report_output = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_qtaim,
                                                     style=Pack(height=300, font_family="monospace"))

        children = [
            self._hint("Critical points and paths are drawn in this viewer's 3D window; the 2D map opens in the "
                       "plot viewer."),
            self._heading("Input"),
            self._row(self._label("Formatted checkpoint:"), self.fchk_input,
                      toga.Button("Browse", on_press=self.browse_fchk, style=Pack(width=90, margin_left=5))),
            self._row(self._label("Output folder:"), self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
            self._row(self._label("Density:"), toga.Label(qa.gw.DENSITY_SECTION, style=Pack(margin_top=6))),
            self._heading("Critical point search - starting points"),
            self._row(*self.mode_switches.values()),
            self._heading("Advanced (blank = Multiwfn default)"),
        ]
        children += [self._row(self._label(label), self.search_inputs[name]) for name, label in SEARCH_FIELDS]
        children += [self._row(self._label(label), self.path_inputs[name]) for name, label in PATH_FIELDS]
        children += [
            self._row(self.run_button, self.display_button, self.clear_button),
            self.status_label,
            self.progress_bar,
            self._heading("Results"),
            self.summary_label,
            self.cp_table,
            self._row(*self.type_switches.values()),
            self._row(self.bond_path_switch, self.ring_path_switch, self.label_switch),
            self._heading("2D plane map"),
            self._row(self._label("Plane:"), self.plane_selection),
            self._row(self._label("Plane atoms (three atoms):"), self.plane_atoms_input),
            self._row(self._label("Offset along the normal (A, XY/XZ/YZ):"), self.plane_offset_input),
            self._row(self._label("Margin around the atoms (A):"), self.margin_input),
            self._row(self._label("Grid points per axis:"), self.grid_points_input),
            self._row(self.contour_switch, self.gradient_switch, self.interbasin_switch),
            self._row(self.map_button),
            self.report_output,
        ]
        return toga.Box(children=children, style=Pack(direction=COLUMN, margin=10))

    # --- dialogs ------------------------------------------------------------
    def _window(self):
        return self.viewer.main_window

    async def _show_error(self, title: str, message: str):
        await self._window().dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        await self._window().dialog(toga.InfoDialog(title, message))

    async def _ask(self, title: str, message: str) -> bool:
        return bool(await self._window().dialog(toga.QuestionDialog(title, message)))

    async def browse_fchk(self, widget):
        try:
            path = await self._window().dialog(toga.OpenFileDialog(
                title="Select a Gaussian formatted checkpoint", file_types=["fchk", "fch"]))
            if path:
                self.fchk_input.value = str(path)
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
    def _parameters(self, cls, fields, inputs, title):
        defaults = cls()
        values = {}
        for name, label in fields:
            text = inputs[name].value
            name_text = f"{title}: '{label.rstrip(':')}'"
            if name in INTEGER_FIELDS:
                values[name] = parse_optional_int(text, getattr(defaults, name), name_text)
            else:
                values[name] = parse_optional_float(text, getattr(defaults, name), name_text)
        return cls(**values).validate()

    def _collect_inputs(self) -> Dict[str, object]:
        fchk = (self.fchk_input.value or "").strip()
        if not fchk:
            raise ValueError("Please select a .fchk or .fch file.")
        if not os.path.isfile(fchk):
            raise ValueError(f"File not found: {fchk}")
        modes = tuple(m for m, switch in self.mode_switches.items() if switch.value)
        if not modes:
            raise ValueError("Select at least one CP search mode (starting points).")
        return {
            "fchk_path": fchk,
            "modes": modes,
            "search": self._parameters(qt.CPSearchParameters, SEARCH_FIELDS, self.search_inputs, "CP search"),
            "paths": self._parameters(qt.PathParameters, PATH_FIELDS, self.path_inputs, "Paths"),
        }

    def _collect_map_inputs(self) -> Dict[str, object]:
        wfn = self.results.wavefunction
        mode = self.plane_selection.value or PLANE_THREE_ATOMS
        if mode == PLANE_THREE_ATOMS:
            text = (self.plane_atoms_input.value or "").strip()
            atoms = parse_atom_triple(text, wfn.n_atoms) if text else qt.default_plane_atoms(wfn.coordinates)
            plane = qt.plane_from_atoms(wfn.coordinates, *atoms, symbols=wfn.atom_symbols())
        else:
            offset = parse_optional_float(self.plane_offset_input.value, 0.0, "The plane offset", positive=False)
            plane = qt.axis_plane(mode, offset / BOHR_RADIUS_ANGSTROM)
        defaults = qa.PlaneMapOptions()
        options = qa.PlaneMapOptions(
            margin_angstrom=parse_optional_float(self.margin_input.value, defaults.margin_angstrom,
                                                 "The map margin", positive=False),
            points=parse_optional_int(self.grid_points_input.value, defaults.points, "Grid points per axis"),
            contours=bool(self.contour_switch.value),
            gradient_lines=bool(self.gradient_switch.value),
            interbasin=bool(self.interbasin_switch.value),
        ).validate()
        return {"plane": plane, "options": options}

    # --- run ----------------------------------------------------------------
    def _apply_progress(self, fraction: float, message: str):
        _assign(self.progress_bar, "value", max(0.0, min(100.0, 100.0 * fraction)))
        _assign(self.status_label, "text", message)

    def _set_buttons(self, enabled: bool):
        for b in (self.run_button, self.display_button, self.clear_button, self.map_button):
            b.enabled = enabled

    def _progress_callback(self):
        loop = asyncio.get_running_loop()

        def on_progress(fraction, message):
            loop.call_soon_threadsafe(self._apply_progress, fraction, message)
        return on_progress

    def _scene_is_displayed(self) -> bool:
        layer = getattr(self.viewer, "qtaim_layer", None)
        return layer is not None and layer.get() is not None

    async def workflow(self, widget):
        try:
            params = self._collect_inputs()
            output_dir = qa.resolve_output_dir(params["fchk_path"], self.output_dir_input.value)
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        self._set_buttons(False)
        self.progress_bar.value = 0
        try:
            results = await asyncio.to_thread(self.analyser.run, progress=self._progress_callback(), **params)
        except Exception as exc:
            self.status_label.text = "Calculation failed."
            self.progress_bar.value = 0
            self._set_buttons(True)
            await self._show_error("QTAIM Calculation Error", str(exc))
            return
        self.results = results  # kept even if the report cannot be written
        self._fill_results()
        try:
            self.report_path = self.analyser.write_outputs(results, output_dir)[0]
            self.report_output.value = Path(self.report_path).read_text(encoding="utf-8")
            self.status_label.text = f"Done in {results.elapsed:.1f} s. Report: {self.report_path}"
            self.progress_bar.value = 100
        except (ValueError, OSError) as exc:
            self.report_output.value = self.analyser.report_text(results)
            await self._show_error("QTAIM Output Error", str(exc))
        finally:
            self._set_buttons(True)
        if self._scene_is_displayed():
            await self._display()  # a scene is on screen: replace it with the new run

    def summary_text(self) -> str:
        counts = self.results.counts()
        total, satisfied = self.results.poincare_hopf()
        return ("   ".join(f"{t} {counts[t]}" for t in CP_TYPES)
                + f"     Poincare-Hopf = {total} " + ("(satisfied)" if satisfied else "(NOT satisfied - CPs missing)"))

    def table_rows(self):
        symbols = self.results.wavefunction.atom_symbols()
        return [(str(c.index), c.kind, qa.assignment_label(c, symbols), f"{c.properties.rho:.6f}",
                 f"{c.properties.laplacian:.6f}", f"{c.properties.ellipticity:.4f}")
                for c in self.results.critical_points]

    def _fill_results(self):
        self._syncing = True
        try:
            _assign(self.summary_label, "text", self.summary_text())
            self.cp_table.data = self.table_rows()
        finally:
            self._syncing = False

    # --- 3D display -------------------------------------------------------------
    def fchk_frame(self):
        """The run's molecule as a viewer frame: [(symbol, (x, y, z) Angstrom), ...]."""
        wfn = self.results.wavefunction
        return [(s, tuple(float(c) for c in xyz))
                for s, xyz in zip(wfn.atom_symbols(), wfn.coordinates * BOHR_RADIUS_ANGSTROM)]

    async def _settle_structure(self):
        """Show the fchk molecule, asking first if the viewer shows a different structure."""
        if self.results is None or self._structure_settled_for is self.results:
            return
        frame = self.fchk_frame()
        # A fresh viewer keeps a 5-atom placeholder in molecule_data with no frames: never worth a question.
        current = self.viewer.get_current_frame_data() if self.viewer.frames else []
        if not current:
            self.viewer.load_designed_molecule(frame)
        elif not same_structure(current, frame):
            if await self._ask("Show the QTAIM molecule?",
                               f"The 3D viewer shows another structure ({len(current)} atoms).\n\n"
                               f"Yes: replace it with the {len(frame)} atoms of {Path(self.results.fchk_path).name}.\n"
                               "No: keep it; the critical points are drawn in the .fchk coordinates."):
                self.viewer.load_designed_molecule(frame)
        self._structure_settled_for = self.results

    def display_options(self) -> Dict[str, object]:
        return {"show_types": {t: bool(s.value) for t, s in self.type_switches.items()},
                "bond_paths": bool(self.bond_path_switch.value),
                "ring_paths": bool(self.ring_path_switch.value),
                "labels": bool(self.label_switch.value)}

    def _refresh_view(self):
        self.viewer.invalidate_scene_cache()
        self.viewer.request_projection_update()

    def _push_to_viewer(self):
        layer = self.viewer.qtaim_layer
        layer.set_display_options(**self.display_options())
        layer.set_scene(scene_from_results(self.results, centre=self.viewer.get_system_center()))
        self._refresh_view()

    async def _display(self):
        await self._settle_structure()
        self._push_to_viewer()

    async def display_results(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Run the analysis first.")
            return
        await self._display()
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)

    def on_display_change(self, widget, **kwargs):
        # Only a scene already on screen follows the switches; nothing is recomputed.
        if self._syncing or self.results is None or not self._scene_is_displayed():
            return
        self.viewer.qtaim_layer.set_display_options(**self.display_options())
        self._refresh_view()

    def on_cp_select(self, widget, **kwargs):
        if self._syncing or self.results is None:
            return
        row = getattr(widget, "selection", None)
        index = int(row.number) if row is not None else None
        self.viewer.qtaim_layer.set_display_options(highlight=index)
        if self._scene_is_displayed():
            self._refresh_view()
        if index is not None:
            cp = self.results.cp(index)
            p = cp.properties
            x, y, z = cp.position * BOHR_RADIUS_ANGSTROM
            _assign(self.status_label, "text",
                    f"CP {index} {cp.kind} at ({x:.4f}, {y:.4f}, {z:.4f}) A: rho {p.rho:.6f}, lap {p.laplacian:.6f}, "
                    f"G {p.g:.6f}, V {p.v:.6f}, H {p.h:.6f} a.u.")

    # --- 2D map -----------------------------------------------------------------
    def map_path(self) -> Path:
        folder = qa.resolve_output_dir(self.results.fchk_path, self.output_dir_input.value)
        return Path(folder) / (Path(self.results.fchk_path).stem + qa.MAP_SUFFIX)

    async def plot_map(self, widget):
        from displayPlots import launch_plot_viewer  # pins Agg: import only when a map is plotted

        if self.results is None:
            await self._show_info("No Results", "Run the analysis first; the map shows its critical points.")
            return
        try:
            params = self._collect_map_inputs()
        except ValueError as exc:
            await self._show_error("2D Map Input Error", str(exc))
            return
        self._set_buttons(False)
        self.progress_bar.value = 0
        try:
            plane_map = await asyncio.to_thread(self.analyser.plane_map, self.results,
                                                progress=self._progress_callback(), **params)
            path = self.map_path()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps([qa.plane_manifest(plane_map)]), encoding="utf-8")
        except (ValueError, OSError, MemoryError) as exc:
            self.status_label.text = "The 2D map was not made."
            self.progress_bar.value = 0
            await self._show_error("2D Map Error", str(exc))
            return
        finally:
            self._set_buttons(True)
        self.progress_bar.value = 100
        self.status_label.text = f"2D map: {path}"
        if not launch_plot_viewer(str(path)):
            await self._show_info("2D Map Saved", f"The plot viewer could not be started. The map was saved as\n"
                                                  f"{path}\nand can be opened with Plotter (JSON plot type).")

    def clear_output(self, widget):
        self.results = None
        self.report_path = None
        self._structure_settled_for = None
        self._syncing = True
        try:
            self.cp_table.data = []
            _assign(self.summary_label, "text", "No results.")
        finally:
            self._syncing = False
        self.report_output.value = HelpGqteaWin.help_qtaim
        if self.viewer is not None:
            self.viewer.qtaim_layer.clear()
            self.viewer.invalidate_scene_cache()
        self.progress_bar.value = 0
        self.status_label.text = "Ready."
