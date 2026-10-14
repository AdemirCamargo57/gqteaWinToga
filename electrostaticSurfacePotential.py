"""Electrostatic surface potential (ESP) mapped on electron-density isosurfaces.

``ESPSurfaceAnalyser`` (headless) reads a Gaussian .fchk/.fch, builds one
isosurface of the Total SCF Density per requested isovalue, evaluates the
exact ESP at every surface vertex and finds the surface extrema with
Multiwfn's rules (see molecularSurface). ``ElectrostaticSurfacePotentialUI``
is the Toga window; ``espViewer.ESPSurfaceViewer`` draws the result.

Two quantities, never to be confused: the *surface* is defined by the
electron density rho (isovalue, a.u. = e/bohr^3); the *mapped property* is
the electrostatic potential V (a.u. = Hartree/e, or kcal/mol).
"""

from __future__ import annotations

import asyncio
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import espIntegrals as ei
import gaussianWavefunction as gw
import molecularSurface as ms
from coordinateConverter import BOHR_RADIUS_ANGSTROM
from espIntegrals import HARTREE_TO_KCAL_MOL
from espViewer import COLOUR_UNITS, ESPSurfaceViewer, SurfaceMarker, UNIT_AU, symmetric_range
from help import AtomicData, HelpGqteaWin

DEFAULT_ISOVALUE = 0.0004
MIN_SPACING, MAX_SPACING = 0.05, 0.6
ENGINE_BUILTIN = "Built-in exact"
ENGINE_CUBEGEN = "Gaussian cubegen"
ENGINES = (ENGINE_BUILTIN, ENGINE_CUBEGEN)
DEFAULT_CUBEGEN_LOCATIONS = (r"C:\G09W\cubegen.exe", r"C:\G16W\cubegen.exe")
EXTREMA_FILE_SUFFIX = "_ESP_extrema.txt"

ProgressCallback = Callable[[float, str], None]


# ----------------------------------------------------------------------
# Input parsing
# ----------------------------------------------------------------------
def parse_isovalues(text: Optional[str]) -> List[float]:
    """'0.0004, 0.001' -> [0.0004, 0.001]; blank -> [DEFAULT_ISOVALUE].

    Order is kept, duplicates dropped. Every value must be a positive
    number below 1 a.u. (a density isovalue, not an ESP value).
    """
    tokens = (text or "").replace(",", " ").replace(";", " ").split()
    if not tokens:
        return [DEFAULT_ISOVALUE]
    values: List[float] = []
    for tok in tokens:
        try:
            v = float(tok)
        except ValueError as exc:
            raise ValueError(f"Isovalue '{tok}' is not a number.") from exc
        if not (0.0 < v < 1.0):
            raise ValueError(
                f"Isovalue {tok} is outside (0, 1) a.u. It is an electron-density "
                "value in e/bohr^3 (e.g. 0.0004 or 0.001), not an ESP value."
            )
        if v not in values:
            values.append(v)
    return values


def parse_spacing(text: Optional[str]) -> float:
    """Grid spacing in bohr; blank -> molecularSurface.DEFAULT_SPACING_BOHR."""
    if not (text or "").strip():
        return ms.DEFAULT_SPACING_BOHR
    try:
        h = float(text)
    except ValueError as exc:
        raise ValueError(f"Grid spacing '{text}' is not a number.") from exc
    if not (MIN_SPACING <= h <= MAX_SPACING):
        raise ValueError(f"Grid spacing must be between {MIN_SPACING} and {MAX_SPACING} bohr.")
    return h


def resolve_output_dir(fchk_path: str, chosen: Optional[str]) -> str:
    """Chosen folder, or the folder of the .fchk when left blank."""
    chosen = (chosen or "").strip()
    return chosen if chosen else str(Path(fchk_path).resolve().parent)


# ----------------------------------------------------------------------
# cubegen engine
# ----------------------------------------------------------------------
def find_cubegen(explicit: Optional[str] = None) -> Optional[str]:
    """Path of Gaussian's cubegen: explicit, $GAUSS_EXEDIR, PATH, usual folders."""
    candidates = []
    if explicit and explicit.strip():
        candidates.append(explicit.strip())
    exedir = os.environ.get("GAUSS_EXEDIR")
    if exedir:
        candidates += [os.path.join(d, "cubegen.exe") for d in exedir.split(os.pathsep)]
    which = shutil.which("cubegen")
    if which:
        candidates.append(which)
    candidates += list(DEFAULT_CUBEGEN_LOCATIONS)
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None


def cubegen_potential(fchk_path: str, points_bohr: np.ndarray, cubegen: str,
                      timeout: float = 3600.0) -> np.ndarray:
    """ESP from Gaussian cubegen at arbitrary points.

    ``cubegen 0 potential=SCF file.fchk out.cub -5 h < points``: with -5 the
    points are read from stdin as free-format x y z **in Angstrom** and the
    result lists x y z V after the atom block. Raises RuntimeError on any
    failure, so the caller can fall back to the built-in engine.
    """
    points_angstrom = np.asarray(points_bohr, dtype=float) * BOHR_RADIUS_ANGSTROM
    env = dict(os.environ, GAUSS_EXEDIR=str(Path(cubegen).parent))
    with tempfile.TemporaryDirectory() as tmp:
        cube = Path(tmp) / "esp.cub"
        stdin = "".join(f"{x:.10f} {y:.10f} {z:.10f}\n" for x, y, z in points_angstrom)
        try:
            subprocess.run([cubegen, "0", "potential=SCF", str(fchk_path), str(cube), "-5", "h"],
                           input=stdin, text=True, env=env, timeout=timeout, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            lines = cube.read_text().splitlines()
        except (OSError, subprocess.SubprocessError) as exc:
            raise RuntimeError(f"cubegen failed: {exc}") from exc
    try:
        natoms = abs(int(lines[2].split()[0]))
        values = np.array([float(l.split()[3]) for l in lines[6 + natoms:] if l.strip()])
    except (IndexError, ValueError) as exc:
        raise RuntimeError("cubegen wrote an output file that could not be read.") from exc
    if values.size != len(points_angstrom):
        raise RuntimeError(f"cubegen returned {values.size} values for {len(points_angstrom)} points.")
    return values


# ----------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------
@dataclass
class SurfaceResult:
    isovalue: float               # electron density, a.u.
    mesh: ms.SurfaceMesh          # bohr
    esp: np.ndarray               # a.u., one per vertex
    extrema: List[ms.SurfaceExtremum]

    @property
    def half_range(self) -> float:
        return symmetric_range(self.esp)


@dataclass
class ESPResults:
    fchk_path: str
    wavefunction: gw.FchkWavefunction
    spacing: float
    engine: str
    surfaces: List[SurfaceResult]
    notes: List[str] = field(default_factory=list)
    elapsed: float = 0.0


class ESPSurfaceAnalyser:
    """Headless pipeline: fchk -> density grid -> isosurfaces -> ESP -> extrema."""

    def run(self, fchk_path: str, isovalues: Sequence[float] = (DEFAULT_ISOVALUE,),
            spacing: float = ms.DEFAULT_SPACING_BOHR, engine: str = ENGINE_BUILTIN,
            cubegen_path: Optional[str] = None,
            progress: Optional[ProgressCallback] = None) -> ESPResults:
        started = time.perf_counter()

        def report(fraction: float, message: str):
            if progress is not None:
                try:
                    progress(fraction, message)
                except Exception:
                    pass  # reporting must never abort an analysis

        if engine not in ENGINES:
            raise ValueError(f"Unknown ESP engine '{engine}'.")
        isovalues = list(isovalues) or [DEFAULT_ISOVALUE]
        notes: List[str] = []

        report(0.02, "Reading the formatted checkpoint...")
        wfn = gw.read_fchk(fchk_path)

        def density(points):
            return gw.electron_density(wfn, points)

        report(0.05, "Building the density grid...")
        radii = np.array([AtomicData.vdw_radii.get(s, AtomicData.DEFAULT_VDW_RADIUS)
                          for s in wfn.atom_symbols()]) / BOHR_RADIUS_ANGSTROM
        origin, shape = ms.grid_box(wfn.coordinates, radii, density, min(isovalues), spacing)
        values = density(ms.grid_points(origin, spacing, shape)).reshape(shape)

        cubegen = None
        if engine == ENGINE_CUBEGEN:
            cubegen = find_cubegen(cubegen_path)
            if cubegen is None:
                notes.append("cubegen was not found; the built-in exact engine was used instead.")
                engine = ENGINE_BUILTIN
        batches = ei.build_hermite_batches(wfn) if engine == ENGINE_BUILTIN else None

        surfaces: List[SurfaceResult] = []
        span = 0.65 / len(isovalues)
        for k, iso in enumerate(isovalues):
            base = 0.30 + k * span
            report(base, f"rho = {iso:g} a.u.: building the isosurface...")
            mesh = ms.build_isosurface(values, origin, spacing, iso, density, ms.density_outward(density))
            closure = mesh.component_euler_characteristics()
            if any(c != 2 for c in closure):
                notes.append(f"rho = {iso:g}: surface components with V+F-E = {closure} "
                             "(2 expected for each closed surface); try a smaller grid spacing.")
            report(base + 0.2 * span, f"rho = {iso:g} a.u.: ESP at {len(mesh.vertices)} surface points...")
            esp = None
            if engine == ENGINE_CUBEGEN:
                try:
                    esp = cubegen_potential(fchk_path, mesh.vertices, cubegen)
                except RuntimeError as exc:
                    notes.append(f"{exc} The built-in exact engine was used instead.")
                    engine = ENGINE_BUILTIN
                    batches = ei.build_hermite_batches(wfn)
            if esp is None:
                esp = ei.electrostatic_potential(
                    wfn, mesh.vertices, batches,
                    progress=lambda f, b=base: report(b + (0.2 + 0.75 * f) * span,
                                                      f"rho = {iso:g} a.u.: ESP {f:.0%}"))
            extrema = ms.find_surface_extrema(esp, mesh.vertices, mesh.neighbours)
            surfaces.append(SurfaceResult(iso, mesh, esp, extrema))

        report(1.0, "Done.")
        return ESPResults(str(fchk_path), wfn, spacing, engine, surfaces, notes,
                          time.perf_counter() - started)

    # --- reporting ------------------------------------------------------
    @staticmethod
    def nearest_atom_label(wfn: gw.FchkWavefunction, position_bohr: np.ndarray) -> str:
        i = int(np.argmin(np.linalg.norm(wfn.coordinates - position_bohr, axis=1)))
        return f"{wfn.atom_symbols()[i]}{i + 1}"

    def extrema_rows(self, results: ESPResults) -> List[Dict[str, object]]:
        """One dict per extremum of every surface, numbered from 1 in file order."""
        rows = []
        for surface in results.surfaces:
            for e in surface.extrema:
                rows.append({
                    "row": len(rows) + 1,
                    "isovalue": surface.isovalue,
                    "kind": e.kind,
                    "scope": e.scope,
                    "esp_au": e.value,
                    "esp_kcal_mol": e.value * HARTREE_TO_KCAL_MOL,
                    "position_angstrom": e.position * BOHR_RADIUS_ANGSTROM,
                    "nearest_atom": self.nearest_atom_label(results.wavefunction, e.position),
                })
        return rows

    def format_extrema_file(self, results: ESPResults) -> str:
        wfn = results.wavefunction
        out = [
            "# gQTEA ESP Surface Extrema",
            f"# source_file            {Path(results.fchk_path).name}",
            f"# density                {gw.DENSITY_SECTION}",
            "# surface_definition     electron density isosurface, rho in a.u. (e/bohr^3)",
            "# mapped_property        electrostatic potential V(r), a.u. (Hartree/e)",
            f"# esp_method             {results.engine}",
            f"# grid_spacing_bohr      {results.spacing:g}",
            f"# vertex_merge_bohr      {ms.MERGE_RATIO * results.spacing:g}",
            "# extremum_criterion     strict 1st+2nd-ring vertex neighbours (Multiwfn convention)",
            f"# hartree_to_kcal_mol    {HARTREE_TO_KCAL_MOL}",
            "# coordinates_unit       Angstrom",
            f"# atoms                  {wfn.n_atoms}  charge {wfn.charge}  multiplicity {wfn.multiplicity}",
        ]
        for s in results.surfaces:
            out.append(
                f"# surface rho_isovalue_au={s.isovalue:.6f} vertices={len(s.mesh.vertices)} "
                f"triangles={len(s.mesh.triangles)} "
                f"area_A2={s.mesh.area() * BOHR_RADIUS_ANGSTROM ** 2:.4f} "
                f"esp_min_au={s.esp.min():.8f} esp_max_au={s.esp.max():.8f}"
            )
        out.append("# row rho_isovalue_au type scope esp_au esp_kcal_mol x_A y_A z_A nearest_atom")
        for r in self.extrema_rows(results):
            x, y, z = r["position_angstrom"]
            out.append(
                f"{r['row']:d} {r['isovalue']:.6f} {r['kind']} {r['scope']} "
                f"{r['esp_au']:.8f} {r['esp_kcal_mol']:.4f} {x:.6f} {y:.6f} {z:.6f} {r['nearest_atom']}"
            )
        return "\n".join(out) + "\n"

    def write_extrema_file(self, results: ESPResults, output_dir: str) -> str:
        path = Path(output_dir) / (Path(results.fchk_path).stem + EXTREMA_FILE_SUFFIX)
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(self.format_extrema_file(results), encoding="utf-8")
        except OSError as exc:
            raise ValueError(f"Could not write {path}: {exc}") from exc
        return str(path)

    def summary_text(self, results: ESPResults, extrema_path: Optional[str] = None) -> str:
        wfn = results.wavefunction
        lines = [
            "ELECTROSTATIC SURFACE POTENTIAL",
            f"File        : {results.fchk_path}",
            f"Molecule    : {wfn.n_atoms} atoms, {wfn.n_electrons} electrons, "
            f"charge {wfn.charge}, multiplicity {wfn.multiplicity}",
            f"Density     : {gw.DENSITY_SECTION}",
            f"ESP engine  : {results.engine}",
            f"Grid        : {results.spacing:g} bohr",
            f"Time        : {results.elapsed:.1f} s",
            "",
        ]
        rows = self.extrema_rows(results)
        for s in results.surfaces:
            lines += [
                f"SURFACE  rho = {s.isovalue:.6f} a.u. (electron density isovalue)",
                f"  vertices {len(s.mesh.vertices)}, area {s.mesh.area() * BOHR_RADIUS_ANGSTROM ** 2:.2f} A^2",
                f"  ESP on this surface: min {s.esp.min():.6f} a.u. ({s.esp.min() * HARTREE_TO_KCAL_MOL:.2f} kcal/mol), "
                f"max {s.esp.max():.6f} a.u. ({s.esp.max() * HARTREE_TO_KCAL_MOL:.2f} kcal/mol)",
                "   row  type     scope      ESP (a.u.)  ESP (kcal/mol)   x (A)     y (A)     z (A)   atom",
            ]
            for r in rows:
                if r["isovalue"] != s.isovalue:
                    continue
                x, y, z = r["position_angstrom"]
                lines.append(f"  {r['row']:4d}  {r['kind']:7s}  {r['scope']:6s}  {r['esp_au']:12.6f}  "
                             f"{r['esp_kcal_mol']:12.3f}  {x:8.4f}  {y:8.4f}  {z:8.4f}   {r['nearest_atom']}")
            lines.append("")
        for note in results.notes:
            lines.append(f"NOTE: {note}")
        if extrema_path:
            lines.append(f"Extrema written to: {extrema_path}")
        return "\n".join(lines)


# ----------------------------------------------------------------------
# Toga UI
# ----------------------------------------------------------------------
OPACITY_LABELS = {"Opaque": 1.0, "Transparent (50 %)": 0.5}


class ElectrostaticSurfacePotentialUI:
    LABEL_WIDTH = 230
    HINT_COLOR = "#666666"

    def __init__(self, *args, **kwargs):
        # gqteaWinToga passes the Button that was pressed; ignore it.
        self.analyser = ESPSurfaceAnalyser()
        self.results: Optional[ESPResults] = None
        self.extrema_path: Optional[str] = None
        self.viewer: Optional[ESPSurfaceViewer] = None
        self._syncing = False
        self.window = toga.Window(title="Electrostatic Surface Potential (ESP)", size=(900, 660),
                                  on_close=self._on_close)
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

    def layout_main_window(self):
        self.fchk_input = toga.TextInput(placeholder="Gaussian formatted checkpoint (.fchk / .fch)",
                                         style=Pack(flex=1))
        self.output_dir_input = toga.TextInput(placeholder="Blank = folder of the .fchk file",
                                               style=Pack(flex=1))
        header = toga.Box(children=[
            self._row(self._label("Formatted checkpoint:"), self.fchk_input,
                      toga.Button("Browse", on_press=self.browse_fchk, style=Pack(width=90, margin_left=5))),
            self._row(self._label("Output folder:"), self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
        ], style=Pack(direction=COLUMN))

        self.isovalue_input = toga.TextInput(placeholder=f"{DEFAULT_ISOVALUE:.6f}", style=Pack(width=220))
        self.spacing_input = toga.TextInput(placeholder=f"{ms.DEFAULT_SPACING_BOHR}", style=Pack(width=120))
        self.engine_selection = toga.Selection(items=list(ENGINES), style=Pack(width=220))
        self.cubegen_input = toga.TextInput(value=find_cubegen() or "", placeholder=DEFAULT_CUBEGEN_LOCATIONS[0],
                                            style=Pack(flex=1))
        setup = toga.Box(children=[
            self._heading("Surface definition - electron density isosurface rho(r)"),
            self._row(self._label("Isovalue(s) rho  (a.u. = e/bohr^3):"), self.isovalue_input),
            self._hint("Blank = 0.000400 (GaussView default). Several: 0.0004, 0.001, 0.002"),
            self._row(self._label("Grid spacing (bohr):"), self.spacing_input),
            self._hint("Blank = 0.25 bohr (Multiwfn default for the ESP). Smaller = smoother and slower."),
            self._row(self._label("Density:"), toga.Label(gw.DENSITY_SECTION, style=Pack(margin_top=6))),
            self._heading("Mapped property - electrostatic potential V(r)"),
            self._row(self._label("ESP engine:"), self.engine_selection),
            self._row(self._label("cubegen executable:"), self.cubegen_input),
            self._hint("Only used with the Gaussian cubegen engine; on failure the built-in engine runs."),
        ], style=Pack(direction=COLUMN, margin=10))

        self.surface_selection = toga.Selection(items=[], on_change=self.on_display_change, style=Pack(width=260))
        self.unit_selection = toga.Selection(items=list(COLOUR_UNITS), on_change=self.on_display_change,
                                             style=Pack(width=140))
        self.opacity_selection = toga.Selection(items=list(OPACITY_LABELS), on_change=self.on_display_change,
                                                style=Pack(width=180))
        self.extrema_switch = toga.Switch("Show extrema markers", value=True, on_change=self.on_display_change)
        self.summary_output = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_esp_surface,
                                                      style=Pack(flex=1, font_family="monospace"))
        results = toga.Box(children=[
            self._row(self._label("Surface shown:"), self.surface_selection),
            self._row(self._label("Colour-bar unit:"), self.unit_selection),
            self._row(self._label("Surface:"), self.opacity_selection, self.extrema_switch),
            self.summary_output,
        ], style=Pack(direction=COLUMN, margin=10, flex=1))

        self.tabs = toga.OptionContainer(content=[
            ("Setup", toga.ScrollContainer(content=setup, style=Pack(flex=1))),
            ("Results", results),
        ], style=Pack(flex=1))

        self.run_button = toga.Button("Run", on_press=self.workflow, style=Pack(width=110, margin=5))
        self.show_button = toga.Button("Show surface", on_press=self.show_surface, style=Pack(width=130, margin=5))
        self.clear_button = toga.Button("Clear", on_press=self.clear_output, style=Pack(width=110, margin=5))
        self.close_button = toga.Button("Close", on_press=self.closeTopLevel, style=Pack(width=110, margin=5))
        self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
        self.progress_bar = toga.ProgressBar(max=100, style=Pack(margin=(0, 5, 5, 5)))

        self.window.content = toga.Box(children=[
            header,
            self._row(self.run_button, self.show_button, self.clear_button, self.close_button),
            self.tabs,
            self.status_label,
            self.progress_bar,
        ], style=Pack(direction=COLUMN, margin=10, flex=1))

    # --- dialogs ------------------------------------------------------------
    async def _show_error(self, title: str, message: str):
        await self.window.dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        await self.window.dialog(toga.InfoDialog(title, message))

    async def browse_fchk(self, widget):
        try:
            path = await self.window.dialog(toga.OpenFileDialog(
                title="Select a Gaussian formatted checkpoint", file_types=["fchk", "fch"]))
            if path:
                self.fchk_input.value = str(path)
        except Exception as exc:
            await self._show_error("File Selection Error", str(exc))

    async def browse_output_dir(self, widget):
        try:
            folder = await self.window.dialog(toga.SelectFolderDialog(title="Choose the output folder"))
            if folder:
                self.output_dir_input.value = str(folder)
        except Exception as exc:
            await self._show_error("Folder Selection Error", str(exc))

    # --- run ----------------------------------------------------------------
    def _collect_inputs(self) -> Dict[str, object]:
        fchk = (self.fchk_input.value or "").strip()
        if not fchk:
            raise ValueError("Please select a .fchk or .fch file.")
        if not os.path.isfile(fchk):
            raise ValueError(f"File not found: {fchk}")
        return {
            "fchk_path": fchk,
            "isovalues": parse_isovalues(self.isovalue_input.value),
            "spacing": parse_spacing(self.spacing_input.value),
            "engine": self.engine_selection.value or ENGINE_BUILTIN,
            "cubegen_path": (self.cubegen_input.value or "").strip() or None,
        }

    def _apply_progress(self, fraction: float, message: str):
        self.progress_bar.value = max(0.0, min(100.0, 100.0 * fraction))
        self.status_label.text = message

    def _set_buttons(self, enabled: bool):
        for b in (self.run_button, self.show_button, self.clear_button):
            b.enabled = enabled

    async def workflow(self, widget):
        try:
            params = self._collect_inputs()
            output_dir = resolve_output_dir(params["fchk_path"], self.output_dir_input.value)
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        self._set_buttons(False)
        self.progress_bar.value = 0
        loop = asyncio.get_running_loop()

        def on_progress(fraction, message):
            loop.call_soon_threadsafe(self._apply_progress, fraction, message)

        try:
            results = await asyncio.to_thread(self.analyser.run, progress=on_progress, **params)
            self.extrema_path = self.analyser.write_extrema_file(results, output_dir)
            self.results = results
            self.summary_output.value = self.analyser.summary_text(results, self.extrema_path)
            self._populate_surface_selection()
            self.tabs.current_tab = 1
            self.status_label.text = f"Done in {results.elapsed:.1f} s. Extrema: {self.extrema_path}"
            self.progress_bar.value = 100
            if self.viewer is not None and self.viewer.viewer_window_is_open():
                self._push_to_viewer()
        except Exception as exc:
            self.status_label.text = "Calculation failed."
            self.progress_bar.value = 0
            await self._show_error("ESP Calculation Error", str(exc))
        finally:
            self._set_buttons(True)

    # --- display --------------------------------------------------------
    @staticmethod
    def surface_label(isovalue: float) -> str:
        return f"rho = {isovalue:.6f} a.u."

    def _populate_surface_selection(self):
        self._syncing = True
        try:
            labels = [self.surface_label(s.isovalue) for s in self.results.surfaces]
            self.surface_selection.items = labels
            self.surface_selection.value = labels[0]
        finally:
            self._syncing = False

    def selected_surface(self) -> Optional[SurfaceResult]:
        if self.results is None:
            return None
        for s in self.results.surfaces:
            if self.surface_label(s.isovalue) == self.surface_selection.value:
                return s
        return self.results.surfaces[0]

    def surface_markers(self, surface: SurfaceResult) -> List[SurfaceMarker]:
        rows = [r for r in self.analyser.extrema_rows(self.results) if r["isovalue"] == surface.isovalue]
        return [SurfaceMarker(r["position_angstrom"], r["kind"], r["scope"], str(r["row"])) for r in rows]

    def _push_to_viewer(self):
        surface = self.selected_surface()
        if surface is None or self.viewer is None:
            return
        wfn = self.results.wavefunction
        self.viewer.set_molecule(wfn.atom_symbols(), wfn.coordinates * BOHR_RADIUS_ANGSTROM)
        self.viewer.set_display_options(unit=self.unit_selection.value or UNIT_AU,
                                        opacity=OPACITY_LABELS.get(self.opacity_selection.value, 1.0),
                                        show_extrema=bool(self.extrema_switch.value))
        self.viewer.set_surface(surface.mesh.vertices * BOHR_RADIUS_ANGSTROM, surface.mesh.triangles,
                                surface.mesh.vertex_normals(), surface.esp, surface.isovalue,
                                self.surface_markers(surface))

    def on_display_change(self, widget):
        if self._syncing or self.results is None or self.viewer is None:
            return
        self._push_to_viewer()

    async def show_surface(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Run the calculation first.")
            return
        if self.viewer is None:
            self.viewer = ESPSurfaceViewer()
            self.viewer.loop = asyncio.get_running_loop()
        self._push_to_viewer()
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)

    def clear_output(self, widget):
        self.results = None
        self.extrema_path = None
        self.summary_output.value = HelpGqteaWin.help_esp_surface
        self._syncing = True
        try:
            self.surface_selection.items = []
        finally:
            self._syncing = False
        if self.viewer is not None:
            self.viewer.clear_surface()
        self.progress_bar.value = 0
        self.status_label.text = "Ready."

    def _on_close(self, window, **kwargs) -> bool:
        if self.viewer is not None:
            self.viewer.request_close()
        return True

    def closeTopLevel(self, widget):
        self._on_close(self.window)
        self.window.close()
