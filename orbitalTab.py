"""The 3D Molecular Viewer's Orbitals tab (``OrbitalTabController``).

Molecular orbitals of a Gaussian .fchk (RHF/RKS, UHF/UKS, ROHF/ROKS), one at
a time, as GaussView's MO editor: a table of every MO with its label,
occupancy and energy, HOMO / LUMO / Previous / Next buttons, and the + and -
phases (or |phi|^2) as isosurfaces in the viewer's own 3D window through
``viewer.orbital_layer`` (a ``cubeRender.CubeSceneLayer``). The numerical
work is ``molecularOrbitals`` (headless). Each computed orbital is a Cube-tab
``CubeEntry`` (meshes cached per isovalue and sign), the last few grids are
kept so stepping back and forth is instant, and colours never re-mesh. The
Style and Transparency controls mirror the viewer's canvas-wide surface style
(its Display tab), which the layer applies when it draws. Widget writes go
through ``cubeTab._assign`` (every write
re-lays out the viewer window) and the change handlers are synchronous and
``_syncing``-guarded, as in the Cube, QTAIM and NTO tabs.
"""

from __future__ import annotations

import asyncio
import math
import os
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Optional

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

import molecularOrbitals as mo
import volumeData as vd
from cubeRender import surface_item
from cubeTab import CubeEntry, _assign, _item_names, canonical_colour, cube_controller_of
from help import HelpGqteaWin
from qtaimTab import fit_table_columns, parse_optional_float

SURFACE_PHASES = "Orbital phases (+/−)"
SURFACE_DENSITY = "Orbital density |φ|²"
SURFACES = (SURFACE_PHASES, SURFACE_DENSITY)
ISOVALUE_LABELS = {SURFACE_PHASES: "Isovalue (a.u.):", SURFACE_DENSITY: "|φ|² isovalue (a.u.):"}
DEFAULT_ISOVALUES = {SURFACE_PHASES: mo.ORBITAL_ISOVALUE, SURFACE_DENSITY: mo.DENSITY_ISOVALUE}
STYLE_SOLID, STYLE_TRANSPARENT, STYLE_MESH = vd.SURFACE_SOLID, vd.SURFACE_TRANSPARENT, vd.SURFACE_MESH
STYLES = vd.SURFACE_STYLES
GRID_CUSTOM = "Custom"
GRID_ITEMS = [*mo.GRID_PRESETS, GRID_CUSTOM]
SIGNS = (1, -1)
DEFAULT_COLOURS = {1: "Red", -1: "Green"}           # GaussView's MO colours
COLOUR_LABELS = {1: "Positive (+)", -1: "Negative (−)"}
SPIN_ITEMS = {mo.ALPHA: "Alpha", mo.BETA: "Beta"}
ONE_SET_ITEM = "One MO set (alpha = beta)"
TABLE_HEADINGS = ["MO", "Label", "Occ.", "Energy (Ha)", "Energy (eV)"]
TABLE_ACCESSORS = ["mo", "label", "occupancy", "energy_ha", "energy_ev"]
TABLE_COLUMN_FRACTIONS = (0.13, 0.24, 0.13, 0.24, 0.24)   # 0.98: room for the vertical scroll bar
IMAGE_FILE_TYPES = ["png", "tif", "tiff", "bmp"]
MAX_TRANSPARENCY = vd.MAX_TRANSPARENCY
GRID_CACHE_SIZE = 3      # amplitude grids kept (about 45 MB each for benzene at Medium)


def select_table_row(table, index: int) -> None:
    """Select and show row ``index``. Toga 0.5's Table has no selection setter, so the
    native ListView is used where there is one (WinForms); best effort, never raises."""
    native = getattr(getattr(table, "_impl", None), "native", None)
    if native is not None:
        try:
            native.SelectedIndices.Clear()
            native.Items[index].Selected = True
            native.Items[index].Focused = True
        except Exception:
            pass
    scroll = getattr(table, "scroll_to_row", None)
    if scroll is not None:
        try:
            scroll(index)
        except Exception:
            pass


class OrbitalTabController:
    """The 3D Molecular Viewer's Orbitals tab."""

    LABEL_WIDTH = 210
    HINT_COLOR = "#666666"

    def __init__(self, viewer):
        self.viewer = viewer
        self.catalogue: Optional[mo.OrbitalCatalogue] = None
        self.spin = mo.ALPHA
        self.number: Optional[int] = None          # the selected MO (1-based)
        self.surface = SURFACE_PHASES
        self.surface_isovalue = mo.ORBITAL_ISOVALUE  # as typed: phi for phases, |phi|^2 for density
        # Mirrors of the viewer's canvas surface style (set in the Display tab or here).
        self.style = getattr(viewer, "canvas_surface_style", STYLE_SOLID)
        self.transparency = float(getattr(viewer, "canvas_transparency", vd.DEFAULT_TRANSPARENCY))
        self.colours: Dict[int, str] = dict(DEFAULT_COLOURS)
        self.use_numba: Optional[bool] = None       # None = if installed
        self._isovalue_text = {surface: "" for surface in SURFACES}
        self._grids: "OrderedDict[tuple, Dict[str, CubeEntry]]" = OrderedDict()  # LRU of computed orbitals
        self._current: Optional[tuple] = None       # key of the orbital on screen (or last computed)
        self._structure_settled_for = None
        self._syncing = False
        self._refreshing = False
        self._refresh_pending = False
        self._displaying = False
        self._display_pending = False
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

    def _button(self, text, handler, width=110):
        return toga.Button(text, on_press=handler, style=Pack(width=width, margin=5))

    def build_tab(self) -> toga.Box:
        self._syncing = True
        try:
            self.fchk_input = toga.TextInput(placeholder="Gaussian formatted checkpoint (.fchk / .fch)",
                                             style=Pack(flex=1))
            self.output_dir_input = toga.TextInput(placeholder="Blank = folder of the .fchk file", style=Pack(flex=1))
            self.load_button = self._button("Load", self.load_file)
            self.clear_button = self._button("Clear", self.clear_output)
            self.status_label = toga.Label("Ready.", style=Pack(margin=6, flex=1))
            self.progress_bar = toga.ProgressBar(max=100, style=Pack(margin=(0, 5, 5, 5)))

            self.summary_label = toga.Label("No file loaded.", style=Pack(margin=(4, 0, 4, 0)))
            self.spin_selection = toga.Selection(items=[ONE_SET_ITEM], on_change=self.on_spin_change,
                                                 style=Pack(width=230))
            self.spin_selection.enabled = False
            self.orbital_table = toga.Table(headings=TABLE_HEADINGS, accessors=TABLE_ACCESSORS, data=[],
                                            on_select=self.on_orbital_select, style=Pack(height=220))
            self.homo_button = self._button("HOMO", self.go_homo, 90)
            self.lumo_button = self._button("LUMO", self.go_lumo, 90)
            self.previous_button = self._button("◀ Previous", self.go_previous, 100)
            self.next_button = self._button("Next ▶", self.go_next, 100)
            self.levels_button = self._button("Energy levels...", self.show_energy_levels, 140)
            self.selected_label = toga.Label("No orbital selected.", style=Pack(margin=(4, 0, 4, 0)))

            self.surface_selection = toga.Selection(items=list(SURFACES), on_change=self.on_surface_change,
                                                    style=Pack(width=230))
            self.surface_selection.value = self.surface
            self.isovalue_label = self._label(ISOVALUE_LABELS[self.surface])
            self.isovalue_input = toga.TextInput(placeholder=f"{DEFAULT_ISOVALUES[self.surface]:g}",
                                                 on_confirm=self.on_grid_confirm, style=Pack(width=120))
            self.grid_selection = toga.Selection(items=GRID_ITEMS, on_change=self.on_grid_change,
                                                 style=Pack(width=120))
            self.grid_selection.value = mo.DEFAULT_GRID
            self.spacing_input = toga.TextInput(placeholder=f"{mo.preset_spacing(mo.DEFAULT_GRID):.4g} bohr",
                                                on_confirm=self.on_grid_confirm, style=Pack(width=110, margin_left=10))
            self.spacing_input.enabled = False
            # Mirrors of the canvas surface style (the Display tab's): changing either changes both.
            self.style_selection = toga.Selection(items=list(STYLES), on_change=self.on_style_change,
                                                  style=Pack(width=120))
            self.style_selection.value = self.style
            self.transparency_slider = toga.Slider(min=0.0, max=MAX_TRANSPARENCY, value=self.transparency,
                                                   on_change=self.on_transparency_change, style=Pack(width=200))
            self.transparency_slider.enabled = self.style == STYLE_TRANSPARENT
            self.colour_inputs = {sign: toga.TextInput(value=DEFAULT_COLOURS[sign], on_confirm=self.on_colour_confirm,
                                                       style=Pack(width=120)) for sign in SIGNS}
            self.display_button = self._button("Display", self.display_orbital)

            self.export_cube_button = self._button("Export cube", self.export_cube, 140)
            self.send_button = self._button("Send to Cube tab", self.send_to_cube_tab, 140)
            self.export_3d_button = self._button("Export 3D image...", self.export_3d_image, 150)
            self.range_input = toga.TextInput(placeholder="e.g. HOMO-2:LUMO+2 or 20-25, 30", style=Pack(flex=1))
            self.export_range_button = self._button("Export cubes", self.export_range, 140)
            self.export_table_button = self._button("Export energy table", self.export_energy_table, 160)
            self.report_output = toga.MultilineTextInput(readonly=True, value=HelpGqteaWin.help_orbitals,
                                                         style=Pack(height=300, font_family="monospace"))
        finally:
            self._syncing = False

        children = [
            toga.Label("Molecular Orbitals", style=Pack(font_size=11, font_weight="bold", margin=(0, 0, 4, 0))),
            self._hint("Orbitals of a Gaussian SCF calculation (RHF/RKS, UHF/UKS, ROHF/ROKS); the + and − phases "
                       "are drawn in this viewer's 3D window, as in GaussView."),
            self._heading("Input"),
            self._row(self._label("Formatted checkpoint:"), self.fchk_input,
                      toga.Button("Browse", on_press=self.browse_fchk, style=Pack(width=90, margin_left=5))),
            self._row(self._label("Output folder:"), self.output_dir_input,
                      toga.Button("Browse", on_press=self.browse_output_dir, style=Pack(width=90, margin_left=5))),
            self._row(self.load_button, self.clear_button),
            self.status_label,
            self.progress_bar,
            self._heading("Orbitals"),
            self.summary_label,
            self._row(self._label("Spin:"), self.spin_selection),
            self.orbital_table,
            self._row(self.homo_button, self.lumo_button, self.previous_button, self.next_button,
                      self.levels_button),
            self.selected_label,
            self._hint("Click a row or use the buttons. Once an orbital is displayed, every choice redraws. "
                       "Energy levels opens the diagram in the plot viewer."),
            self._heading("Display"),
            self._row(self._label("Surface:"), self.surface_selection),
            self._row(self.isovalue_label, self.isovalue_input),
            self._row(self._label("Grid quality:"), self.grid_selection, self.spacing_input),
            self._hint("Coarse, Medium, Fine = 3, 6, 12 points per bohr (cubegen's grids); Custom uses the typed "
                       "spacing in bohr. Enter applies a typed value."),
            self._row(self._label("Style:"), self.style_selection),
            self._row(self._label("Transparency:"), self.transparency_slider),
            self._hint("Style and Transparency are the canvas-wide settings of the Display tab: they apply to "
                       "every surface drawn."),
            self._row(self._label(f"{COLOUR_LABELS[1]} / {COLOUR_LABELS[-1]}:"), self.colour_inputs[1],
                      toga.Label("  ", style=Pack(width=10)), self.colour_inputs[-1]),
            self._hint("Colours: Red, Green, Blue, Yellow, Orange, Cyan, Magenta, Purple, White, Grey or #rrggbb "
                       "(Enter applies). |φ|² uses the positive colour."),
            self._row(self.display_button),
            self._heading("Export"),
            self._row(self.export_cube_button, self.send_button, self.export_3d_button),
            self._row(self._label("Orbital range:"), self.range_input, self.export_range_button),
            self._hint("One cube per orbital, all on one grid (so the Cube tab can add or subtract them). "
                       "Cubes follow the Surface choice."),
            self._row(self.export_table_button),
            self._hint("Files are written to the output folder: <name>_MO<n>[_alpha|_beta][_density].cube, "
                       "<name>_MO_energies.txt, <name>_MO_levels.json."),
            self.report_output,
        ]
        self._attach_column_fitting()
        return toga.Box(children=children, style=Pack(direction=COLUMN, margin=10))

    def _attach_column_fitting(self):
        native = getattr(getattr(self.orbital_table, "_impl", None), "native", None)
        if native is None:
            return
        try:
            native.Resize += lambda sender, event: fit_table_columns(self.orbital_table, TABLE_COLUMN_FRACTIONS)
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
            self.fchk_input.value = str(path)

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
    def _fchk_path(self) -> str:
        path = (self.fchk_input.value or "").strip()
        if not path:
            raise ValueError("Please select a .fchk or .fch file.")
        if not os.path.isfile(path):
            raise ValueError(f"File not found: {path}")
        return path

    def _is_loaded(self, path: str) -> bool:
        return self.catalogue is not None and Path(self.catalogue.path) == Path(path)

    def _spacing(self) -> float:
        grid = self.grid_selection.value or mo.DEFAULT_GRID
        if grid == GRID_CUSTOM:
            return parse_optional_float(self.spacing_input.value, mo.preset_spacing(mo.DEFAULT_GRID),
                                        "The grid spacing")
        return mo.preset_spacing(grid)

    def _typed_isovalue(self) -> float:
        return parse_optional_float(self.isovalue_input.value, DEFAULT_ISOVALUES[self.surface], "The isovalue")

    def _collect_display_inputs(self) -> Dict[str, object]:
        """The ``molecularOrbitals.orbital_volumes`` arguments of the selected orbital.

        The grid is always the amplitude's: for |phi|^2 the typed isovalue is
        that of the density, so the grid uses its square root (the same surface).
        """
        if self.catalogue is None:
            raise ValueError("Load a .fchk file first.")
        if self.number is None:
            raise ValueError("Choose an orbital.")
        isovalue = self._typed_isovalue()
        return {
            "spin": self.spin,
            "numbers": [self.number],
            "spacing": self._spacing(),
            "isovalue": math.sqrt(isovalue) if self.surface == SURFACE_DENSITY else isovalue,
        }

    def _collect_colours(self) -> Dict[int, str]:
        colours = {}
        for sign, widget in self.colour_inputs.items():
            try:
                colours[sign] = canonical_colour(widget.value)
            except ValueError as exc:
                raise ValueError(f"{COLOUR_LABELS[sign]}: {exc}") from None
        return colours

    def output_dir(self) -> str:
        return mo.resolve_output_dir(self.catalogue.path, self.output_dir_input.value)

    @property
    def density_mode(self) -> bool:
        return self.surface == SURFACE_DENSITY

    # --- progress -------------------------------------------------------------
    def _apply_progress(self, fraction: float, message: str):
        _assign(self.progress_bar, "value", max(0.0, min(100.0, 100.0 * fraction)))
        _assign(self.status_label, "text", message)

    def _progress_callback(self, message: str):
        loop = asyncio.get_running_loop()

        def on_progress(fraction):
            loop.call_soon_threadsafe(self._apply_progress, fraction, message)
        return on_progress

    def _set_buttons(self, enabled: bool):
        for b in (self.load_button, self.clear_button, self.display_button, self.homo_button, self.lumo_button,
                  self.previous_button, self.next_button, self.levels_button, self.export_cube_button,
                  self.send_button, self.export_range_button, self.export_table_button):
            b.enabled = enabled

    # --- load ---------------------------------------------------------------------
    async def load_file(self, widget):
        try:
            path = self._fchk_path()
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        await self._load(path)

    async def _load(self, path: str) -> bool:
        self._set_buttons(False)
        _assign(self.status_label, "text", "Reading the file...")
        try:
            catalogue = await asyncio.to_thread(mo.load_orbitals, path)
        except (ValueError, MemoryError) as exc:
            _assign(self.status_label, "text", "The file was not loaded.")
            await self._show_error("Orbitals Input Error", str(exc))
            return False
        finally:
            self._set_buttons(True)
        was_displayed = self._scene_is_displayed()
        self.catalogue = catalogue
        self._grids.clear()
        self._current = None
        self._structure_settled_for = None
        self.spin = mo.ALPHA
        self.number = catalogue.homo[mo.ALPHA] or 1
        self._fill_catalogue()
        _assign(self.status_label, "text", f"{catalogue.n_mo} MOs read ({catalogue.kind}). The HOMO is selected; "
                                           "press Display to draw it.")
        if was_displayed:
            await self._display()  # an orbital is on screen: show the new file's HOMO
        return True

    def table_rows(self):
        return [(str(o.number), o.label, f"{o.occupancy:g}", f"{o.energy:.6f}", f"{o.energy_ev:.4f}")
                for o in self.catalogue.orbitals[self.spin]]

    def selected_text(self) -> str:
        if self.catalogue is None or self.number is None:
            return "No orbital selected."
        o = self.catalogue.orbital(self.spin, self.number)
        return (f"Selected: MO {o.number} ({o.label}){self.catalogue.spin_text(self.spin)}, occupancy "
                f"{o.occupancy:g}, {o.energy:.6f} Ha ({o.energy_ev:.4f} eV).")

    def _fill_catalogue(self):
        c = self.catalogue
        items = [SPIN_ITEMS[s] for s in c.spins] if c.unrestricted else [ONE_SET_ITEM]
        self._syncing = True
        try:
            _assign(self.summary_label, "text", c.summary())
            if items != _item_names(self.spin_selection):
                self.spin_selection.items = items
            _assign(self.spin_selection, "value", items[0])
            _assign(self.spin_selection, "enabled", c.unrestricted)
        finally:
            self._syncing = False
        self._fill_table()

    def _fill_table(self):
        self._syncing = True
        try:
            self.orbital_table.data = self.table_rows()
            fit_table_columns(self.orbital_table, TABLE_COLUMN_FRACTIONS)
        finally:
            self._syncing = False
        self._show_selected()

    def _show_selected(self):
        self._syncing = True
        try:
            _assign(self.selected_label, "text", self.selected_text())
            if self.number is not None:
                select_table_row(self.orbital_table, self.number - 1)
        finally:
            self._syncing = False

    # --- orbital selection ------------------------------------------------------
    def select_orbital(self, number: int):
        """Make MO ``number`` the selected one (clamped to the file); redraw if an orbital is shown."""
        if self.catalogue is None:
            self._start(self._show_info("No File", "Load a .fchk file first."))
            return
        number = max(1, min(self.catalogue.n_mo, int(number)))
        if number == self.number:
            return
        self.number = number
        self._show_selected()
        if self._scene_is_displayed():
            self._start(self._display())

    def go_homo(self, widget):
        if self.catalogue is not None:
            self.select_orbital(self.catalogue.homo[self.spin] or 1)
        else:
            self.select_orbital(1)

    def go_lumo(self, widget):
        if self.catalogue is None:
            self.select_orbital(1)
            return
        lumo = self.catalogue.lumo(self.spin)
        if lumo is None:
            _assign(self.status_label, "text", "This MO set has no virtual orbital.")
            return
        self.select_orbital(lumo)

    def go_previous(self, widget):
        self.select_orbital((self.number or 2) - 1)

    def go_next(self, widget):
        self.select_orbital((self.number or 0) + 1)

    # Synchronous, guarded handlers: see the Cube tab (an async handler would
    # run after _syncing is cleared and replay the tab's own writes as edits).
    def on_orbital_select(self, widget, **kwargs):
        if self._syncing or self.catalogue is None:
            return
        row = getattr(widget, "selection", None)
        try:
            number = int(getattr(row, "mo", ""))
        except (TypeError, ValueError):
            return
        if number == self.number:
            return
        self.number = number
        _assign(self.selected_label, "text", self.selected_text())
        if self._scene_is_displayed():
            self._start(self._display())

    def on_spin_change(self, widget, **kwargs):
        if self._syncing or self.catalogue is None or not self.catalogue.unrestricted:
            return
        spin = next((s for s, label in SPIN_ITEMS.items() if label == widget.value), mo.ALPHA)
        if spin == self.spin:
            return
        self.spin = spin  # the MO number is kept, so alpha n and beta n can be compared
        self._fill_table()
        if self._scene_is_displayed():
            self._start(self._display())

    def on_surface_change(self, widget, **kwargs):
        if self._syncing:
            return
        surface = widget.value if widget.value in SURFACES else SURFACE_PHASES
        if surface == self.surface:
            return
        self._isovalue_text[self.surface] = self.isovalue_input.value or ""
        self.surface = surface
        self._syncing = True
        try:
            _assign(self.isovalue_label, "text", ISOVALUE_LABELS[surface])
            _assign(self.isovalue_input, "value", self._isovalue_text[surface])
            _assign(self.isovalue_input, "placeholder", f"{DEFAULT_ISOVALUES[surface]:g}")
        finally:
            self._syncing = False
        if self._scene_is_displayed():
            self._start(self._display())

    def on_grid_change(self, widget, **kwargs):
        if self._syncing:
            return
        _assign(self.spacing_input, "enabled", widget.value == GRID_CUSTOM)
        if self._scene_is_displayed():
            self._start(self._display())

    def on_grid_confirm(self, widget, **kwargs):
        if self._syncing or not self._scene_is_displayed():
            return
        self._start(self._display())

    # Style and Transparency set the viewer's canvas-wide surface style, which
    # the layers apply when they draw: nothing here is recomputed.
    def on_style_change(self, widget, **kwargs):
        if self._syncing or widget.value not in STYLES:
            return
        self._apply_canvas_style(style=widget.value)

    def on_transparency_change(self, widget, **kwargs):
        if self._syncing:
            return
        self._apply_canvas_style(transparency=float(widget.value))

    def _apply_canvas_style(self, style=None, transparency=None):
        apply = getattr(self.viewer, "apply_surface_style", None)
        if apply is not None:
            apply(style=style, transparency=transparency)  # calls sync_surface_widgets back
        else:
            self.sync_surface_widgets(self.style if style is None else style,
                                      self.transparency if transparency is None else transparency)

    def sync_surface_widgets(self, style: str, transparency: float):
        """Show the canvas surface style (set here or in the Display tab) in this tab's widgets."""
        self.style, self.transparency = style, float(transparency)
        self._syncing = True
        try:
            _assign(self.style_selection, "value", style)
            _assign(self.transparency_slider, "value", float(transparency))
            _assign(self.transparency_slider, "enabled", style == STYLE_TRANSPARENT)
        finally:
            self._syncing = False

    def on_colour_confirm(self, widget, **kwargs):
        if self._syncing:
            return
        sign = next(s for s, w in self.colour_inputs.items() if w is widget)
        try:
            self.colours[sign] = canonical_colour(widget.value)
        except ValueError as exc:
            self._syncing = True
            try:
                _assign(widget, "value", self.colours[sign])  # put the colour in use back
            finally:
                self._syncing = False
            self._start(self._show_error("Colour", f"{COLOUR_LABELS[sign]}: {exc}"))
            return
        if self._scene_is_displayed():
            self._start(self.refresh_scene())

    # --- 3D display -------------------------------------------------------------
    def fchk_frame(self):
        """The molecule as a viewer frame: [(symbol, (x, y, z) Angstrom), ...]."""
        wfn = self.catalogue.wavefunction
        return [(s, tuple(float(c) for c in xyz))
                for s, xyz in zip(wfn.atom_symbols(), wfn.coordinates * mo.BOHR_RADIUS_ANGSTROM)]

    async def _settle_structure(self):
        """Show the fchk molecule, asking first if the viewer shows a different structure."""
        wfn = self.catalogue.wavefunction
        if self._structure_settled_for is wfn:
            return
        frame = self.fchk_frame()
        # A fresh viewer keeps a 5-atom placeholder in molecule_data with no frames: never worth a question.
        current = self.viewer.get_current_frame_data() if self.viewer.frames else []
        if not current:
            self.viewer.load_designed_molecule(frame)
        elif not vd.same_structure(current, frame):
            if await self._ask("Show the orbital's molecule?",
                               f"The 3D viewer shows another structure ({len(current)} atoms).\n\n"
                               f"Yes: replace it with the {len(frame)} atoms of {Path(self.catalogue.path).name}.\n"
                               "No: keep it; the orbitals are drawn in the .fchk coordinates."):
                self.viewer.load_designed_molecule(frame)
        self._structure_settled_for = wfn

    def _scene_is_displayed(self) -> bool:
        layer = getattr(self.viewer, "orbital_layer", None)
        return layer is not None and bool(layer.items())

    def orbital_entries(self, spin: str, numbers, spacing: float, isovalue: float,
                        progress=None) -> Dict[str, CubeEntry]:
        """The orbital as a Cube-tab entry (mesh cache), computed once per grid and kept in a small LRU.
        Runs in a worker thread."""
        number = numbers[0]
        key = (id(self.catalogue), spin, number, spacing, isovalue)
        if key in self._grids:
            self._grids.move_to_end(key)
        else:
            volume = mo.orbital_volumes(self.catalogue, spin, [number], spacing, isovalue, progress,
                                        self.use_numba)[0]
            entry = CubeEntry.from_volume(volume.title, volume, saved=False)
            entry.isovalue = isovalue
            self._grids[key] = {"phases": entry}
            while len(self._grids) > GRID_CACHE_SIZE:
                self._grids.popitem(last=False)
        self._current = key
        return self._grids[key]

    def surface_entry(self) -> Optional[CubeEntry]:
        """The current orbital's entry for the Surface choice (|phi|^2 is derived once, on demand)."""
        if self._current is None or self._current not in self._grids:
            return None
        entries = self._grids[self._current]
        if not self.density_mode:
            return entries["phases"]
        if "density" not in entries:
            _, spin, number, _, _ = self._current
            volume = mo.density_volume(self.catalogue, entries["phases"].volume, spin, number)
            entries["density"] = CubeEntry.from_volume(volume.title, volume, saved=False)
        entry = entries["density"]
        entry.isovalue = self.surface_isovalue
        return entry

    def compute_scene_items(self):
        """(render items, notes) of the current orbital, Surface and colours. Runs in a worker thread.

        The items are plain opaque isosurfaces: the canvas surface style
        (Solid/Transparent/Mesh) is applied by the layer when it draws.
        """
        entry = self.surface_entry()
        if entry is None:
            return [], []
        items, notes = [], []
        what = "|φ|²" if self.density_mode else "orbital"
        for sign in ((1,) if self.density_mode else SIGNS):
            mesh = entry.mesh(sign, self.use_numba)
            if mesh is None:
                notes.append(f"No {'+' if sign > 0 else '−'}{entry.isovalue:g} surface for this {what}.")
                continue
            if mesh.clipped:
                notes.append("The surface reaches the edge of the grid; raise the isovalue.")
            items.append(surface_item(mesh, vd.parse_colour(self.colours[sign]), 1.0))
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
                        claim("orbital_layer")
                self.viewer.orbital_layer.set_items(items)
                self.viewer.invalidate_scene_cache()
                self.viewer.request_projection_update()
                if not self._refresh_pending:
                    break
        finally:
            self._refreshing = False
        shown = len(items)
        head = self.selected_text().replace("Selected: ", "") if shown else ""
        _assign(self.status_label, "text", " ".join(
            [f"{shown} isosurface{'s' if shown != 1 else ''} shown." + (f" {head}" if head else "")] + notes))

    async def _display(self) -> bool:
        """Compute (or reuse) the selected orbital's grid and draw it. Requests arriving while one
        runs (fast Next presses) are coalesced: the latest selection is drawn once."""
        if self._displaying:
            self._display_pending = True
            return True
        self._displaying = True
        try:
            while True:
                self._display_pending = False
                ok = await self._display_once()
                if not ok or not self._display_pending:
                    return ok
        finally:
            self._displaying = False

    async def _display_once(self) -> bool:
        try:
            params = self._collect_display_inputs()
            colours = self._collect_colours()
            surface_isovalue = self._typed_isovalue()
        except ValueError as exc:
            await self._show_error("Display Input Error", str(exc))
            return False
        self.colours = colours
        self.surface_isovalue = surface_isovalue
        message = f"Computing MO {self.number} on the grid..."
        _assign(self.status_label, "text", message)
        try:
            await asyncio.to_thread(self.orbital_entries, progress=self._progress_callback(message), **params)
        except (ValueError, MemoryError) as exc:
            _assign(self.status_label, "text", "The orbital was not computed.")
            await self._show_error("Orbital Grid Error", str(exc))
            return False
        _assign(self.progress_bar, "value", 100.0)
        await self._settle_structure()
        await self.refresh_scene()
        return True

    async def display_orbital(self, widget):
        try:
            path = self._fchk_path()
        except ValueError as exc:
            await self._show_error("Input Error", str(exc))
            return
        if not self._is_loaded(path) and not await self._load(path):
            return
        if not await self._display():
            return
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)

    # --- export -------------------------------------------------------------------
    def _current_names(self) -> str:
        """File stem of the orbital computed last: <name>_MO<n>[_alpha|_beta][_density]."""
        _, spin, number, _, _ = self._current
        spin_part = spin if self.catalogue.unrestricted else None  # only unrestricted files have two sets
        return mo.cube_filename(self.catalogue.stem, number, spin_part, self.density_mode)[:-len(".cube")]

    async def _compute_selected(self) -> Optional[CubeEntry]:
        """The selected orbital's entry for the Surface choice, computing its grid if needed."""
        try:
            params = self._collect_display_inputs()
            self.surface_isovalue = self._typed_isovalue()
        except ValueError as exc:
            await self._show_error("Export Input Error", str(exc))
            return None
        try:
            await asyncio.to_thread(self.orbital_entries, progress=self._progress_callback(
                f"Computing MO {self.number} on the grid..."), **params)
            return await asyncio.to_thread(self.surface_entry)
        except (ValueError, MemoryError) as exc:
            await self._show_error("Orbital Grid Error", str(exc))
            return None

    async def export_cube(self, widget):
        """The selected orbital (or its |phi|^2) on the chosen grid, as one cube file."""
        if self.catalogue is None:
            await self._show_info("No File", "Load a .fchk file first.")
            return
        self._set_buttons(False)
        try:
            entry = await self._compute_selected()
            if entry is None:
                return
            path = Path(self.output_dir()) / f"{self._current_names()}.cube"
            written = await asyncio.to_thread(mo.write_volume, entry.volume, path, self.use_numba)
        except ValueError as exc:
            await self._show_error("Export Cube Error", str(exc))
            return
        finally:
            self._set_buttons(True)
        _assign(self.status_label, "text", f"Cube written: {written}")

    async def export_range(self, widget):
        """One cube per orbital of the typed range, all on one grid."""
        if self.catalogue is None:
            await self._show_info("No File", "Load a .fchk file first.")
            return
        try:
            numbers = mo.parse_orbital_spec(self.range_input.value, self.catalogue, self.spin)
            params = self._collect_display_inputs()
        except ValueError as exc:
            await self._show_error("Export Input Error", str(exc))
            return
        self._set_buttons(False)
        _assign(self.progress_bar, "value", 0.0)
        try:
            written = await asyncio.to_thread(
                mo.write_cubes, self.catalogue, self.output_dir(), self.spin, numbers, params["spacing"],
                params["isovalue"], self.density_mode,
                self._progress_callback(f"Writing {len(numbers)} cube files..."), self.use_numba)
        except (ValueError, MemoryError) as exc:
            await self._show_error("Export Cubes Error", str(exc))
            return
        finally:
            self._set_buttons(True)
        _assign(self.progress_bar, "value", 100.0)
        _assign(self.status_label, "text", f"{len(written)} cube files written to {Path(written[0]).parent}.")

    async def send_to_cube_tab(self, widget):
        """Add the displayed orbital to the Cube tab (slices, arithmetic, saving)."""
        entry = self.surface_entry() if self.catalogue is not None else None
        # Asked for only when there is something to send: it builds an unvisited Cube tab.
        cube = cube_controller_of(self.viewer) if entry is not None else None
        if cube is None:
            await self._show_info("Nothing to Send", "Display an orbital first; it is sent as shown.")
            return
        added = cube.add_result(entry.volume, self._current_names())
        added.visible = False  # already drawn here; tick it in the Cube tab to draw it there too
        cube.sync_widgets(selected=added.name)
        _assign(self.status_label, "text", f"Sent to the Cube tab (unticked): {added.name}.")

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
        stem = self._current_names() if self.catalogue is not None and self._current is not None else "orbital"
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

    async def export_energy_table(self, widget):
        if self.catalogue is None:
            await self._show_info("No File", "Load a .fchk file first.")
            return
        try:
            path = mo.write_energy_table(self.catalogue, self.output_dir())
        except ValueError as exc:
            await self._show_error("Export Energy Table Error", str(exc))
            return
        self.report_output.value = mo.energy_table_text(self.catalogue)
        _assign(self.status_label, "text", f"Energy table written: {path}")

    async def show_energy_levels(self, widget):
        """The MO energy-level diagram, saved as <stem>_MO_levels.json and opened in the plot viewer."""
        from displayPlots import launch_plot_viewer  # pins Agg: import only when a diagram is shown

        if self.catalogue is None:
            await self._show_info("No File", "Load a .fchk file first.")
            return
        try:
            path = mo.write_levels(mo.energy_levels_entry(self.catalogue),
                                   Path(self.output_dir()) / mo.levels_filename(self.catalogue.stem))
        except ValueError as exc:
            await self._show_error("Energy Levels Error", str(exc))
            return
        _assign(self.status_label, "text", f"Energy levels: {path}")
        if not launch_plot_viewer(path):
            await self._show_info("Energy Levels Saved", f"The plot viewer could not be started. The diagram was "
                                                         f"saved as\n{path}\nand can be opened with Plotter "
                                                         "(JSON plot type).")

    def clear_output(self, widget):
        self.catalogue = None
        self.number = None
        self.spin = mo.ALPHA
        self._grids.clear()
        self._current = None
        self._structure_settled_for = None
        self._syncing = True
        try:
            self.orbital_table.data = []
            if _item_names(self.spin_selection) != [ONE_SET_ITEM]:
                self.spin_selection.items = [ONE_SET_ITEM]
            _assign(self.spin_selection, "enabled", False)
            _assign(self.summary_label, "text", "No file loaded.")
            _assign(self.selected_label, "text", "No orbital selected.")
        finally:
            self._syncing = False
        self.report_output.value = HelpGqteaWin.help_orbitals
        self.viewer.orbital_layer.clear()
        self.viewer.invalidate_scene_cache()
        _assign(self.progress_bar, "value", 0.0)
        _assign(self.status_label, "text", "Ready.")
