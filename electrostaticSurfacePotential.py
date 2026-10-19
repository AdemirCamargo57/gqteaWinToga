"""Electrostatic surface potential (ESP) mapped on electron-density isosurfaces.

``ESPSurfaceAnalyser`` (headless) reads a Gaussian .fchk/.fch, builds one
isosurface of the Total SCF Density per requested isovalue, evaluates the
exact ESP at every surface vertex and finds the surface extrema with
Multiwfn's rules (see molecularSurface), and writes the ESP (and density)
on its grid as Gaussian cubes. ``ESPTabController`` is the 3D Molecular
Viewer's ESP tab; the viewer draws the result through its
``espSurfaceLayer.ESPSurfaceLayer``. This module must not import
molecularViewer or espViewer at module level: molecularViewer imports it.

Two quantities, never to be confused: the *surface* is defined by the
electron density rho (isovalue, a.u. = e/bohr^3); the *mapped property* is
the electrostatic potential V (a.u. = Hartree/e, or kcal/mol).
"""

from __future__ import annotations

import asyncio
import os
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import toga
from scipy.spatial import cKDTree
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import cubeFile as cf
import espIntegrals as ei
import gaussianWavefunction as gw
import molecularSurface as ms
from coordinateConverter import BOHR_RADIUS_ANGSTROM
from espIntegrals import HARTREE_TO_KCAL_MOL
from espSurfaceLayer import COLOUR_UNITS, UNIT_AU, SurfaceMarker, symmetric_range
from gaussianUtilities import find_gaussian_utility
from help import AtomicData, HelpGqteaWin
from volumeData import same_structure  # shared with the Cube tab

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
    return find_gaussian_utility("cubegen", explicit)


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
    # The density grid the surfaces were built from (bohr); the ESP cube uses it too.
    grid_origin: Optional[np.ndarray] = None
    grid_shape: Optional[Tuple[int, int, int]] = None
    density_grid: Optional[np.ndarray] = None
    cubegen_path: Optional[str] = None  # the cubegen that ran, when that engine was used


ESP_CUBE_SUFFIX = "_ESP.cube"
DENSITY_CUBE_SUFFIX = "_density.cube"
CUBE_FILE_TYPES = ["cube"]
ON_NUCLEUS_DISTANCE = 1e-6  # bohr: a grid point this close to a nucleus has V = infinity
NUCLEUS_SHIFT = 1e-3        # bohr: such a point is evaluated this far away along +x


def default_esp_cube_name(fchk_path: str) -> str:
    return Path(fchk_path).stem + ESP_CUBE_SUFFIX


def density_cube_path(esp_cube_path) -> Path:
    """'<name>_ESP.cube' -> '<name>_density.cube' beside it; any other name gets '_density' appended."""
    path = Path(esp_cube_path)
    stem = path.stem[:-len("_ESP")] if path.stem.endswith("_ESP") else path.stem
    return path.with_name(stem + DENSITY_CUBE_SUFFIX)


def _ascii(text: str) -> str:
    """Cube header lines are ASCII (cubeFile.write_cube writes ascii)."""
    return text.encode("ascii", "replace").decode("ascii")


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
                          time.perf_counter() - started, grid_origin=np.asarray(origin, dtype=float),
                          grid_shape=tuple(int(n) for n in shape), density_grid=values,
                          cubegen_path=cubegen if engine == ENGINE_CUBEGEN else None)

    # --- ESP and density cubes --------------------------------------------
    def esp_on_grid(self, results: ESPResults, progress: Optional[Callable[[float], None]] = None):
        """(V on every point of the run's density grid, shape grid_shape; notes).

        Uses the run's engine (cubegen falls back to the built-in engine with
        a note). V is infinite on a nucleus, so a grid point within
        ON_NUCLEUS_DISTANCE of one is evaluated NUCLEUS_SHIFT away along +x
        and a note says so: the cube never holds inf or nan.
        """
        if results.grid_origin is None or results.grid_shape is None:
            raise ValueError("These results hold no density grid; run the calculation again.")
        wfn = results.wavefunction
        points = ms.grid_points(results.grid_origin, results.spacing, results.grid_shape)
        distance, _ = cKDTree(wfn.coordinates).query(points)
        on_nucleus = distance < ON_NUCLEUS_DISTANCE
        points[on_nucleus, 0] += NUCLEUS_SHIFT
        notes: List[str] = []
        values = None
        if results.engine == ENGINE_CUBEGEN and results.cubegen_path:
            try:
                values = cubegen_potential(results.fchk_path, points, results.cubegen_path)
            except RuntimeError as exc:
                notes.append(f"{exc} The built-in exact engine was used for the cube.")
        if values is None:
            values = ei.electrostatic_potential(wfn, points, ei.build_hermite_batches(wfn), progress=progress)
        if on_nucleus.any():
            notes.append(f"{int(on_nucleus.sum())} grid point(s) lie on a nucleus, where V is infinite; "
                         f"they hold V evaluated {NUCLEUS_SHIFT:g} bohr away along +x.")
        if not np.all(np.isfinite(values)):
            raise ValueError("The ESP is not finite at some grid points.")
        if progress is not None:
            try:
                progress(1.0)
            except Exception:
                pass
        return np.asarray(values, dtype=float).reshape(results.grid_shape), notes

    def _cube(self, results: ESPResults, title: str, comment: str, values: np.ndarray) -> cf.CubeData:
        wfn = results.wavefunction
        return cf.CubeData(_ascii(title), _ascii(comment), wfn.atomic_numbers, wfn.nuclear_charges,
                           wfn.coordinates, results.grid_origin, results.spacing, values)

    def esp_cube(self, results: ESPResults, values: np.ndarray) -> cf.CubeData:
        stem = Path(results.fchk_path).stem
        return self._cube(results, f"{stem} electrostatic potential V(r)",
                          f"gQTEA ESP, a.u. (Hartree/e), {results.engine} engine", values)

    def density_cube(self, results: ESPResults) -> cf.CubeData:
        stem = Path(results.fchk_path).stem
        return self._cube(results, f"{stem} electron density rho(r)",
                          f"gQTEA ESP, {gw.DENSITY_SECTION}, a.u. (e/bohr^3)", results.density_grid)

    def write_cubes(self, results: ESPResults, esp_path, with_density: bool = True,
                    progress: Optional[Callable[[float], None]] = None):
        """Write the ESP cube (and the density cube on the same grid); (paths, notes)."""
        path = Path(esp_path)
        if not path.suffix:
            path = path.with_suffix(".cube")
        values, notes = self.esp_on_grid(results, progress)
        paths = [cf.write_cube(path, self.esp_cube(results, values))]
        if with_density:
            paths.append(cf.write_cube(density_cube_path(path), self.density_cube(results)))
        return paths, notes

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
# Canvas background; the first entry is the default (GaussView's black).
BACKGROUND_COLOURS = {
    "Black": (0.0, 0.0, 0.0),
    "White": (1.0, 1.0, 1.0),
    "Dark grey": (0.25, 0.25, 0.25),
    "Light grey": (0.85, 0.85, 0.85),
}
IMAGE_FILE_TYPES = ["png", "tif", "tiff", "bmp"]


class ESPTabController:
    """The 3D Molecular Viewer's ESP tab (formerly the Tools > ESP window).

    It owns the run, the results and the widgets; the surface, its extremum
    markers and the colour bar are drawn in the viewer's own 3D window through
    ``viewer.esp_layer`` (an ``espSurfaceLayer.ESPSurfaceLayer``).
    """

    LABEL_WIDTH = 230
    HINT_COLOR = "#666666"

    def __init__(self, viewer):
        self.viewer = viewer
        self.analyser = ESPSurfaceAnalyser()
        self.results: Optional[ESPResults] = None
        self.extrema_path: Optional[str] = None
        self._syncing = False
        # The results whose structure question has been answered (asked once per run).
        self._structure_settled_for: Optional[ESPResults] = None

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
        self.fchk_input = toga.TextInput(placeholder="Gaussian formatted checkpoint (.fchk / .fch)",
                                         style=Pack(flex=1))
        self.output_dir_input = toga.TextInput(placeholder="Blank = folder of the .fchk file",
                                               style=Pack(flex=1))
        self.isovalue_input = toga.TextInput(placeholder=f"{DEFAULT_ISOVALUE:.6f}", style=Pack(width=220))
        self.spacing_input = toga.TextInput(placeholder=f"{ms.DEFAULT_SPACING_BOHR}", style=Pack(width=120))
        self.engine_selection = toga.Selection(items=list(ENGINES), style=Pack(width=220))
        self.cubegen_input = toga.TextInput(value=find_cubegen() or "", placeholder=DEFAULT_CUBEGEN_LOCATIONS[0],
                                            style=Pack(flex=1))

        self.run_button = toga.Button("Run", on_press=self.workflow, style=Pack(width=110, margin=5))
        self.show_button = toga.Button("Show surface", on_press=self.show_surface, style=Pack(width=130, margin=5))
        self.clear_button = toga.Button("Clear", on_press=self.clear_output, style=Pack(width=110, margin=5))
        self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
        self.progress_bar = toga.ProgressBar(max=100, style=Pack(margin=(0, 5, 5, 5)))

        self.surface_selection = toga.Selection(items=[], on_change=self.on_display_change, style=Pack(width=260))
        self.unit_selection = toga.Selection(items=list(COLOUR_UNITS), on_change=self.on_display_change,
                                             style=Pack(width=140))
        self.opacity_selection = toga.Selection(items=list(OPACITY_LABELS), on_change=self.on_display_change,
                                                style=Pack(width=180))
        self.extrema_switch = toga.Switch("Show extrema markers", value=True, on_change=self.on_display_change)
        self.background_selection = toga.Selection(items=list(BACKGROUND_COLOURS),
                                                   on_change=self.on_background_change, style=Pack(width=140))
        self.export_button = toga.Button("Export image...", on_press=self.export_image,
                                         style=Pack(width=130, margin_left=10))

        self.density_cube_switch = toga.Switch("Also save the electron density cube (same grid)", value=True)
        self.cube_button = toga.Button("Save ESP cube...", on_press=self.export_cubes,
                                       style=Pack(width=160, margin=5))
        self.summary_output = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_esp_surface,
                                                      style=Pack(height=300, font_family="monospace"))

        return toga.Box(children=[
            self._hint("The surface, its extrema and the colour bar are drawn in this viewer's 3D window."),
            self._heading("Input"),
            self._row(self._label("Formatted checkpoint:"), self.fchk_input,
                      toga.Button("Browse", on_press=self.browse_fchk, style=Pack(width=90, margin_left=5))),
            self._row(self._label("Output folder:"), self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
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
            self._row(self.run_button, self.show_button, self.clear_button),
            self.status_label,
            self.progress_bar,
            self._heading("Results"),
            self._row(self._label("Surface shown:"), self.surface_selection),
            self._row(self._label("Colour-bar unit:"), self.unit_selection),
            self._row(self._label("Surface:"), self.opacity_selection, self.extrema_switch),
            self._row(self._label("Canvas background:"), self.background_selection, self.export_button),
            self._heading("ESP cube"),
            self._row(self.density_cube_switch),
            self._row(self.cube_button),
            self._hint("V(r) in a.u. on the run's density grid (its box and spacing). The suggested name is "
                       "<fchk name>_ESP.cube; the density cube is written beside it as <name>_density.cube."),
            self.summary_output,
        ], style=Pack(direction=COLUMN, margin=10))

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
        for b in (self.run_button, self.show_button, self.clear_button, self.cube_button):
            b.enabled = enabled

    def _surface_is_displayed(self) -> bool:
        layer = getattr(self.viewer, "esp_layer", None)
        return layer is not None and layer.get() is not None

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
            self.status_label.text = f"Done in {results.elapsed:.1f} s. Extrema: {self.extrema_path}"
            self.progress_bar.value = 100
            if self._surface_is_displayed():
                await self._display_results()  # a surface is on screen: replace it with the new run
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
        # Only a loaded structure counts (XYZ, Design, a cube, a previous ESP
        # run): a fresh MolecularViewer keeps a 5-atom placeholder in
        # molecule_data with no frames, and that is never worth a question.
        current = self.viewer.get_current_frame_data() if self.viewer.frames else []
        if not current:
            self.viewer.load_designed_molecule(frame)
        elif not same_structure(current, frame):
            if await self._ask("Show the ESP molecule?",
                               f"The 3D viewer shows another structure ({len(current)} atoms).\n\n"
                               f"Yes: replace it with the {len(frame)} atoms of "
                               f"{Path(self.results.fchk_path).name}.\n"
                               "No: keep it; the surface is drawn in the .fchk coordinates."):
                self.viewer.load_designed_molecule(frame)
        self._structure_settled_for = self.results

    def _push_to_viewer(self):
        surface = self.selected_surface()
        if surface is None or self.viewer is None:
            return
        layer = self.viewer.esp_layer
        layer.set_display_options(unit=self.unit_selection.value or UNIT_AU,
                                  opacity=OPACITY_LABELS.get(self.opacity_selection.value, 1.0),
                                  show_extrema=bool(self.extrema_switch.value))
        layer.set_surface(surface.mesh.vertices * BOHR_RADIUS_ANGSTROM, surface.mesh.triangles,
                          surface.mesh.vertex_normals(), surface.esp, surface.isovalue,
                          self.surface_markers(surface), centre=self.viewer.get_system_center())
        self.viewer.invalidate_scene_cache()
        self.viewer.request_projection_update()

    async def _display_results(self):
        await self._settle_structure()
        self.viewer.set_background_color(self.selected_background())
        self._push_to_viewer()

    def on_display_change(self, widget):
        # Only a surface already on screen follows the Results controls.
        if self._syncing or self.results is None or self.viewer is None or not self._surface_is_displayed():
            return
        self._push_to_viewer()

    async def show_surface(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Run the calculation first.")
            return
        await self._display_results()
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)

    # --- canvas background and image export ------------------------------
    def selected_background(self):
        return BACKGROUND_COLOURS.get(self.background_selection.value, BACKGROUND_COLOURS["Black"])

    def on_background_change(self, widget):
        # A view setting, not a result: applies with or without a run.
        if self.viewer is not None:
            self.viewer.set_background_color(self.selected_background())

    def default_image_name(self) -> str:
        surface = self.selected_surface()
        if surface is None:
            return "ESP_surface.png"
        return f"{Path(self.results.fchk_path).stem}_ESP_rho{surface.isovalue:.6f}.png"

    async def _ask_image_path(self, suggested: str):
        return await self._window().dialog(toga.SaveFileDialog(
            title="Export the ESP surface image", suggested_filename=suggested,
            file_types=IMAGE_FILE_TYPES))

    async def export_image(self, widget):
        """Save the 3D canvas exactly as shown: same view, colours and background.

        The frame is captured first, the instant the button is pressed, and
        only then is the file name asked for -- so the dialog can neither
        cover the canvas nor let the view drift before the capture.
        """
        from molecularViewer import save_canvas_image  # molecularViewer imports this module

        if self.viewer is None or not self.viewer.viewer_window_is_open():
            await self._show_info("No Image", "Open the 3D view with 'Show surface' first; "
                                              "the image is taken from that window as shown.")
            return
        try:
            image = await self.viewer.capture_frame()
        except RuntimeError as exc:
            await self._show_error("Export Image Error", str(exc))
            return
        try:
            path = await self._ask_image_path(self.default_image_name())
            if not path:
                return
            saved = save_canvas_image(image, str(path))
        except (ValueError, OSError) as exc:
            await self._show_error("Export Image Error", str(exc))
            return
        height, width = image.shape[:2]
        self.status_label.text = f"Image saved ({width} x {height} px): {saved}"

    # --- ESP cube export --------------------------------------------------
    def suggested_cube_path(self) -> str:
        """<output folder>/<fchk name>_ESP.cube (the save dialog opens in that folder)."""
        folder = resolve_output_dir(self.results.fchk_path, self.output_dir_input.value)
        return str(Path(folder) / default_esp_cube_name(self.results.fchk_path))

    async def _ask_cube_path(self, suggested: str):
        return await self._window().dialog(toga.SaveFileDialog(
            title="Save the ESP cube", suggested_filename=suggested, file_types=CUBE_FILE_TYPES))

    async def export_cubes(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Run the calculation first; the cube uses its density grid.")
            return
        path = await self._ask_cube_path(self.suggested_cube_path())
        if not path:
            return
        with_density = bool(self.density_cube_switch.value)
        self._set_buttons(False)
        self.progress_bar.value = 0
        loop = asyncio.get_running_loop()

        def on_progress(fraction):
            loop.call_soon_threadsafe(self._apply_progress, fraction, f"ESP on the grid: {fraction:.0%}")

        try:
            paths, notes = await asyncio.to_thread(self.analyser.write_cubes, self.results, str(path),
                                                   with_density, on_progress)
        except (ValueError, MemoryError) as exc:
            self.status_label.text = "The ESP cube was not saved."
            self.progress_bar.value = 0
            await self._show_error("ESP Cube Error", str(exc))
            return
        finally:
            self._set_buttons(True)
        self.progress_bar.value = 100
        self.status_label.text = " ".join(["Saved: " + ", ".join(str(p) for p in paths) + "."] + notes)

    def clear_output(self, widget):
        self.results = None
        self.extrema_path = None
        self._structure_settled_for = None
        self.summary_output.value = HelpGqteaWin.help_esp_surface
        self._syncing = True
        try:
            self.surface_selection.items = []
        finally:
            self._syncing = False
        if self.viewer is not None:
            self.viewer.esp_layer.clear()
            self.viewer.invalidate_scene_cache()
        self.progress_bar.value = 0
        self.status_label.text = "Ready."
