"""The 3D Molecular Viewer's Cube tab: Gaussian cube files as isosurfaces and slices.

``CubeTabController`` owns the tab's state -- the loaded cubes
(``CubeEntry``), their display settings, the slice and the structure the
cubes brought -- and its headless methods (``read_entries``, ``add_entries``,
``adopt_structure``, ``compute_scene_items``, ``compute_operation``,
``save_entry``, ...) are the test seams; ``build_tab`` creates the widgets and
the ``on_*`` handlers wrap the seams in dialogs, the same split as the
Design tab's ``MolecularDesignController``.

What is drawn goes to ``viewer.cube_layer`` (``cubeRender.CubeSceneLayer``),
which ``MolecularViewerUI.render_scene_extras`` draws. Reading files,
meshing and grid arithmetic run in ``asyncio.to_thread`` so the window stays
responsive; meshes are cached per (isovalue, sign) so a colour, opacity or
slice change does not re-mesh.
"""

from __future__ import annotations

import asyncio
import functools
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import toga
from toga.style import Pack
from toga.style.pack import CENTER, COLUMN, LEFT, ROW

import cubeFile as cf
import volumeData as vd
from cubeRender import slice_item, surface_item
from help import HelpGqteaWin

SHOW_BOTH = "+ and −"
SHOW_POSITIVE = "+ only"
SHOW_NEGATIVE = "− only"
SHOW_MODES = {SHOW_BOTH: (1, -1), SHOW_POSITIVE: (1,), SHOW_NEGATIVE: (-1,)}

OP_ADD = "A + B"
OP_SUBTRACT = "A − B"
OP_SCALE = "c × A"
OPERATION_LABELS = {OP_ADD: "add", OP_SUBTRACT: "subtract", OP_SCALE: "scale"}

AXIS_LABELS = ("a (x)", "b (y)", "c (z)")
RANGE_LABELS = {"Symmetric ±max|V|": "symmetric", "Min to max": "minmax"}

DEFAULT_POSITIVE_COLOUR = "Blue"
DEFAULT_NEGATIVE_COLOUR = "Red"
MIN_OPACITY = 0.1
CUBE_FILE_TYPES = ["cube", "cub"]
IMAGE_FILE_TYPES = ["png", "tif", "tiff", "bmp"]
DEFAULT_SLICE_AXIS = 2


def canonical_colour(text: str) -> str:
    """The stored form of a colour: the preset's own name, or #RRGGBB in capitals."""
    vd.parse_colour(text)  # validates
    value = text.strip()
    for name in vd.COLOUR_PRESETS:
        if value.lower() == name.lower():
            return name
    return value.upper()


def _assign(widget, attr: str, value) -> None:
    """Write a widget property only when it actually changes.

    On WinForms every Toga widget write re-lays out the whole viewer window
    (~70 ms measured, seven tabs of controls). Loading one cube used to write
    ~25 mostly unchanged values in a row, which showed as the window
    flickering for seconds. Numbers are compared with a tolerance because a
    slider reads back a rounded value.
    """
    current = getattr(widget, attr, None)
    if isinstance(value, float) and isinstance(current, (int, float)) and not isinstance(current, bool):
        if math.isclose(current, value, rel_tol=1e-9, abs_tol=1e-9):
            return
    elif current == value:
        return
    setattr(widget, attr, value)


def _item_names(selection) -> List[str]:
    """The strings a Toga Selection holds (its rows carry them as ``.value``)."""
    return [getattr(item, "value", item) for item in selection.items]


def rgb_to_hex(rgb) -> str:
    """(r, g, b) in [0, 1] -> '#RRGGBB'."""
    return "#" + "".join(f"{int(round(255 * min(1.0, max(0.0, float(c))))):02X}" for c in rgb)


def pick_colour(initial, owner=None) -> Optional[Tuple[float, float, float]]:
    """The standard Windows colour dialog, opened on ``initial`` (RGB in [0, 1]).

    Returns the chosen RGB, or None when cancelled. Raises RuntimeError where
    WinForms is not available; the typed-colour field still works there.
    The dialog is modal (``owner`` is the Toga window's native Form).
    """
    try:
        import clr  # pythonnet; load the assemblies here rather than rely on toga_winforms having done it
        clr.AddReference("System.Drawing")
        clr.AddReference("System.Windows.Forms")
        from System.Drawing import Color
        from System.Windows.Forms import ColorDialog, DialogResult
    except Exception as exc:
        raise RuntimeError("The colour picker needs Windows; type a colour name or #rrggbb "
                           "in the field and press Enter instead.") from exc
    dialog = ColorDialog()
    dialog.FullOpen = True
    dialog.AnyColor = True
    r, g, b = (int(round(255 * min(1.0, max(0.0, float(c))))) for c in initial)
    dialog.Color = Color.FromArgb(r, g, b)
    result = dialog.ShowDialog(owner) if owner is not None else dialog.ShowDialog()
    if result != DialogResult.OK:
        return None
    chosen = dialog.Color
    return chosen.R / 255.0, chosen.G / 255.0, chosen.B / 255.0


# (heading, width in px; None = takes the remaining width)
CUBE_LIST_COLUMNS = (("Name", None), ("Grid", 110), ("Min", 80), ("Max", 80), ("Integral", 80), ("Show", 46))
CUBE_LIST_HEIGHT = 130


class CubeList:
    """The Cube files list: a header, then one row per cube - name, grid and
    statistics - ending in a **Show** checkbox (a toga.Switch, which is a
    native CheckBox on WinForms).

    It replaces a toga.Table, which cannot hold checkboxes (on WinForms it is
    a virtual ListView). Rows are rebuilt only when the cubes or their
    statistics change, as one new content box swapped into the scroll
    container (one layout pass); a visibility change only writes that box.
    ``on_toggle(name, shown)`` is called for every box the user ticks.
    """

    def __init__(self, on_toggle):
        self.on_toggle = on_toggle
        self.names: List[str] = []
        self.rows: Optional[List[Tuple[str, ...]]] = None
        self.switches: Dict[str, toga.Switch] = {}
        self.writes = 0  # widget writes made by show(): rebuilds and box changes
        self.container = toga.ScrollContainer(horizontal=False,
                                              style=Pack(height=CUBE_LIST_HEIGHT, margin=(0, 0, 6, 0)))
        self.container.content = self._build([], [], [])

    @staticmethod
    def _cell(text: str, width: Optional[int], bold: bool = False) -> toga.Label:
        weight = "bold" if bold else "normal"
        if width is None:
            return toga.Label(text, style=Pack(flex=1, margin=(2, 6, 2, 0), font_weight=weight))
        return toga.Label(text, style=Pack(width=width, margin=(2, 6, 2, 0), font_weight=weight))

    def _build(self, names, rows, visible) -> toga.Box:
        box = toga.Box(style=Pack(direction=COLUMN, margin=(2, 4, 2, 4)))
        header = toga.Box(style=Pack(direction=ROW))
        for title, width in CUBE_LIST_COLUMNS:
            header.add(self._cell(title, width, bold=True))
        box.add(header)
        switches = {}
        for name, row, shown in zip(names, rows, visible):
            line = toga.Box(style=Pack(direction=ROW, align_items=CENTER))
            for text, (_, width) in zip(row, CUBE_LIST_COLUMNS[:-1]):
                line.add(self._cell(text, width))
            switch = toga.Switch("", value=bool(shown), on_change=functools.partial(self._toggled, name),
                                 style=Pack(width=CUBE_LIST_COLUMNS[-1][1]))
            switches[name] = switch
            line.add(switch)
            box.add(line)
        self.switches = switches
        return box

    def _toggled(self, name: str, widget, **kwargs):
        self.on_toggle(name, bool(widget.value))

    def show(self, names: List[str], rows, visible: List[bool]):
        """Bring the list up to date, writing only what changed."""
        rows = [tuple(r) for r in rows]
        if list(names) != self.names or rows != self.rows:
            self.container.content = self._build(names, rows, visible)  # one layout pass
            self.names, self.rows = list(names), rows
            self.writes += 1
            return
        for name, shown in zip(names, visible):
            switch = self.switches[name]
            if switch.value != bool(shown):
                switch.value = bool(shown)
                self.writes += 1


def default_cube_filename(name: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    return f"{stem or 'cube'}.cube"


@dataclass(eq=False)
class CubeEntry:
    name: str
    volume: cf.CubeVolume
    path: Optional[str] = None
    saved: bool = True
    isovalue: float = vd.ORBITAL_ISOVALUE
    show: str = SHOW_BOTH
    positive_colour: str = DEFAULT_POSITIVE_COLOUR
    negative_colour: str = DEFAULT_NEGATIVE_COLOUR
    opacity: float = 1.0
    visible: bool = True
    statistics: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    meshes: Dict[Tuple[float, int], Optional[vd.IsoMesh]] = field(default_factory=dict)

    @classmethod
    def from_volume(cls, name: str, volume: cf.CubeVolume, path: Optional[str] = None,
                    saved: bool = True) -> "CubeEntry":
        # A density has no negative lobe to draw; asking for one would only
        # produce a "no surface" note on every refresh.
        show = SHOW_BOTH if vd.is_signed(volume) else SHOW_POSITIVE
        return cls(name, volume, path, saved, vd.default_isovalue(volume), show,
                   statistics=vd.value_statistics(volume))

    def mesh(self, sign: int, use_numba: Optional[bool] = None) -> Optional[vd.IsoMesh]:
        """The isosurface at sign * isovalue, cached; only the current isovalue is kept."""
        key = (self.isovalue, sign)
        meshes = self.meshes
        if key not in meshes:
            meshes = {k: m for k, m in meshes.items() if k[0] == self.isovalue}
            meshes[key] = vd.isosurface(self.volume, self.isovalue, sign, use_numba)
            self.meshes = meshes
        return meshes[key]

    def table_row(self) -> Tuple[str, str, str, str, str]:
        lo, hi, total = self.statistics
        name = self.name if self.saved else f"{self.name} (unsaved)"
        return (name, " x ".join(str(n) for n in self.volume.values.shape),
                f"{lo:.4g}", f"{hi:.4g}", f"{total:.4g}")


class CubeTabController:
    def __init__(self, viewer):
        self.viewer = viewer
        self.entries: List[CubeEntry] = []
        self.show_volumes = True
        self.slice_cube: Optional[str] = None
        self.slice_axis = DEFAULT_SLICE_AXIS
        self.slice_index = 0
        self.slice_range = "symmetric"
        self.slice_visible = False
        self.structure_volume: Optional[cf.CubeVolume] = None  # whose atoms the viewer shows
        self.use_numba: Optional[bool] = None  # None = if installed
        self._syncing = False
        self._refreshing = False
        self._refresh_pending = False
        self._last_status = None   # last message sent to the viewer's status line
        self._tasks = set()        # handler work in flight (see _start)
        self._swatch_hex = None    # colour last written to the background swatch
        # Widgets (build_tab)
        self.cube_list: Optional[CubeList] = None
        self.background_input = None
        self.background_swatch = None
        self.show_molecule_switch = None
        self.show_volumes_switch = None
        self.iso_cube_selection = None
        self.isovalue_input = None
        self.show_selection = None
        self.positive_colour_input = None
        self.negative_colour_input = None
        self.opacity_slider = None
        self.visible_switch = None
        self.slice_cube_selection = None
        self.slice_axis_selection = None
        self.slice_slider = None
        self.slice_index_label = None
        self.slice_range_selection = None
        self.slice_switch = None
        self.op_a_selection = None
        self.op_selection = None
        self.op_b_selection = None
        self.factor_input = None
        self.result_name_input = None
        self.status_label = None

    # ------------------------------------------------------------------
    # Entries
    # ------------------------------------------------------------------
    def names(self) -> List[str]:
        return [e.name for e in self.entries]

    def entry(self, name: Optional[str]) -> CubeEntry:
        for e in self.entries:
            if e.name == name:
                return e
        raise ValueError(f"No cube named '{name}' is loaded.")

    def unique_name(self, base: str) -> str:
        names = set(self.names())
        if base not in names:
            return base
        k = 2
        while f"{base} ({k})" in names:
            k += 1
        return f"{base} ({k})"

    @staticmethod
    def read_entries(path) -> List[CubeEntry]:
        """Read one cube file (thread-safe; touches no state). ValueError on any problem."""
        volumes = cf.read_cube_volumes(path)
        stem = Path(path).stem
        return [CubeEntry.from_volume(stem if v.orbital is None else f"{stem} MO {v.orbital}", v, str(path))
                for v in volumes]

    def add_entries(self, entries: List[CubeEntry]) -> List[str]:
        """Add read entries under unique names; notes for cubes whose atoms differ from the shown ones."""
        notes = []
        for e in entries:
            e.name = self.unique_name(e.name)
            self.entries.append(e)
            if self.structure_volume is not None:
                mismatch = vd.atoms_mismatch(self.structure_volume, e.volume)
                if mismatch:
                    notes.append(f"{e.name}: {mismatch.replace('The two cubes', 'It and the displayed cube')} "
                                 "The displayed structure is unchanged.")
        if self.slice_cube is None and self.entries:
            self.set_slice(cube=self.entries[0].name, axis=DEFAULT_SLICE_AXIS)
        return notes

    def add_result(self, volume: cf.CubeVolume, name: str) -> CubeEntry:
        entry = CubeEntry.from_volume(self.unique_name(name), volume, None, saved=False)
        self.entries.append(entry)
        return entry

    def remove(self, name: str):
        entry = self.entry(name)
        self.entries.remove(entry)
        if self.slice_cube == name:
            self.slice_cube = None
            if self.entries:
                self.set_slice(cube=self.entries[0].name)

    # ------------------------------------------------------------------
    # Structure
    # ------------------------------------------------------------------
    def displays_cube_molecule(self, entry: CubeEntry) -> bool:
        """True when the viewer already shows this cube's own atoms (e.g. put there by the ESP tab)."""
        frame = vd.structure_frame(entry.volume)
        return bool(frame) and bool(self.viewer.frames) and vd.same_structure(
            self.viewer.get_current_frame_data(), frame)

    def structure_question_needed(self, entry: Optional[CubeEntry] = None) -> bool:
        """Ask before replacing a structure the user already has on screen - unless
        that structure is the cube's own molecule, which needs no replacing."""
        if self.structure_volume is not None or not self.viewer.frames:
            return False
        return entry is None or not self.displays_cube_molecule(entry)

    def adopt_structure(self, entry: CubeEntry):
        frame = vd.structure_frame(entry.volume)
        if not frame:
            raise ValueError(f"{entry.name} lists no atoms, so it has no structure to show. "
                             "Load the molecule's XYZ file first.")
        self.viewer.load_designed_molecule(frame)
        self.structure_volume = entry.volume

    def set_show_structure(self, show: bool):
        self.viewer.show_structure = bool(show)

    # ------------------------------------------------------------------
    # Slice
    # ------------------------------------------------------------------
    def slice_plane_count(self) -> int:
        if self.slice_cube is None:
            return 0
        return int(self.entry(self.slice_cube).volume.values.shape[self.slice_axis])

    def set_slice(self, cube: Optional[str] = None, axis: Optional[int] = None, index: Optional[int] = None,
                  range_mode: Optional[str] = None, visible: Optional[bool] = None):
        """Validate everything first, then apply; a new cube or axis starts at its middle plane."""
        new_cube = self.slice_cube if cube is None else cube
        new_axis = self.slice_axis if axis is None else int(axis)
        if new_axis not in (0, 1, 2):
            raise ValueError(f"The slice axis must be a, b or c, not {axis}.")
        new_index = 0
        if new_cube is not None:
            n = int(self.entry(new_cube).volume.values.shape[new_axis])
            moved = new_cube != self.slice_cube or new_axis != self.slice_axis
            if index is None:
                new_index = n // 2 if moved else min(self.slice_index, n - 1)
            else:
                new_index = int(index)
                if not 0 <= new_index < n:
                    raise ValueError(f"The slice plane must be between 1 and {n}, not {new_index + 1}.")
        new_range = self.slice_range if range_mode is None else range_mode
        if new_range not in vd.COLOUR_RANGE_MODES:
            raise ValueError(f"Unknown colour range '{range_mode}'.")
        self.slice_cube, self.slice_axis, self.slice_index = new_cube, new_axis, new_index
        self.slice_range = new_range
        if visible is not None:
            self.slice_visible = bool(visible)

    def slice_text(self) -> str:
        if self.slice_cube is None:
            raise ValueError("Load a cube first; there is no slice to export.")
        volume = self.entry(self.slice_cube).volume
        return vd.slice_table_text(volume, self.slice_axis, self.slice_index, name=self.slice_cube)

    def default_slice_filename(self) -> str:
        stem = default_cube_filename(self.slice_cube or "cube")[:-len(".cube")]
        return f"{stem}_slice_{vd.AXIS_NAMES[self.slice_axis]}{self.slice_index + 1}.dat"

    # ------------------------------------------------------------------
    # Scene
    # ------------------------------------------------------------------
    def compute_scene_items(self):
        """(render items, notes) for everything visible. Runs in a worker thread."""
        items, notes = [], []
        if not self.show_volumes:
            return items, notes
        for e in list(self.entries):
            if not e.visible:
                continue
            missing = []
            for sign in SHOW_MODES[e.show]:
                mesh = e.mesh(sign, self.use_numba)
                if mesh is None:
                    missing.append(sign)
                    continue
                colour = vd.parse_colour(e.positive_colour if sign > 0 else e.negative_colour)
                items.append(surface_item(mesh, colour, e.opacity))
            if missing:
                where = "±" if len(missing) == 2 else ("+" if missing[0] > 0 else "−")
                lo, hi, _ = e.statistics
                notes.append(f"{e.name}: no surface at {where}{e.isovalue:g} "
                             f"(values run from {lo:.4g} to {hi:.4g}).")
        if self.slice_visible and self.slice_cube in self.names():
            volume = self.entry(self.slice_cube).volume
            points, values = vd.slice_plane(volume, self.slice_axis, self.slice_index)
            lo, hi = vd.colour_range(values, self.slice_range)
            flat = values.ravel()
            colours = (vd.diverging_colours(flat, hi) if self.slice_range == "symmetric"
                       else vd.rainbow_colours(flat, lo, hi))
            items.append(slice_item(points, colours, 1.0))
        return items, notes

    def apply_scene(self, items):
        self.viewer.cube_layer.set_items(items)
        self.viewer.invalidate_scene_cache()
        self.viewer.request_projection_update()

    @staticmethod
    def scene_summary(items, notes) -> str:
        surfaces = sum(1 for i in items if i.kind == "surface")
        text = f"{surfaces} isosurface{'s' if surfaces != 1 else ''} shown"
        if any(i.kind == "slice" for i in items):
            text += " and a slice"
        return " ".join([text + "."] + list(notes))

    async def refresh_scene(self):
        """Recompute and push the scene; changes made meanwhile are coalesced into one rerun."""
        if self._refreshing:
            self._refresh_pending = True
            return
        self._refreshing = True
        try:
            while True:
                self._refresh_pending = False
                try:
                    items, notes = await asyncio.to_thread(self.compute_scene_items)
                except (ValueError, MemoryError) as exc:
                    self.set_status(f"Could not draw the cubes: {exc}")
                    return
                if not self._refresh_pending:
                    self.apply_scene(items)
                    self.set_status(self.scene_summary(items, notes))
                    return
        finally:
            self._refreshing = False

    # ------------------------------------------------------------------
    # Operations and export
    # ------------------------------------------------------------------
    def compute_operation(self, a_name, op_label, b_name, factor_text, result_name):
        """(volume, name, notes) for A + B, A − B or c × A. Thread-safe; adds nothing."""
        if op_label not in OPERATION_LABELS:
            raise ValueError(f"Unknown grid operation '{op_label}'.")
        op = OPERATION_LABELS[op_label]
        if not a_name:
            raise ValueError("Choose cube A for this operation.")
        a = self.entry(a_name)
        notes = []
        if op == "scale":
            try:
                factor = float(str(factor_text).strip())
            except ValueError:
                raise ValueError("The factor c must be a number, e.g. 2 or -0.5.") from None
            default = f"{factor:g} × {a.name}"
            name = (result_name or "").strip() or default
            volume = vd.scale(a.volume, factor, title=name)
        else:
            if not b_name:
                raise ValueError("Choose cube B for this operation.")
            b = self.entry(b_name)
            default = f"{a.name} {'+' if op == 'add' else '−'} {b.name}"
            name = (result_name or "").strip() or default
            volume = vd.combine(a.volume, b.volume, op, title=name)
            mismatch = vd.atoms_mismatch(a.volume, b.volume)
            if mismatch:
                notes.append(f"{mismatch} The result keeps the atoms of {a.name}.")
        return volume, name, notes

    def save_entry(self, name: str, path) -> str:
        entry = self.entry(name)
        written = cf.write_cube_volume(path, entry.volume)
        entry.path, entry.saved = written, True
        return written

    # ------------------------------------------------------------------
    # Widgets -> state
    # ------------------------------------------------------------------
    def apply_isosurface_fields(self) -> Optional[str]:
        """Copy every isosurface widget into the selected entry; an error message or None."""
        try:
            entry = self.entry(self.iso_cube_selection.value)
            try:
                isovalue = float(str(self.isovalue_input.value).strip())
            except ValueError:
                isovalue = float("nan")
            if not (isovalue > 0 and isovalue < float("inf")):
                raise ValueError("The isovalue must be a positive number; the negative surface is "
                                 "drawn at minus this value.")
            show = self.show_selection.value
            if show not in SHOW_MODES:
                raise ValueError(f"Unknown surface choice '{show}'.")
            positive = canonical_colour(self.positive_colour_input.value)
            negative = canonical_colour(self.negative_colour_input.value)
            opacity = min(1.0, max(MIN_OPACITY, float(self.opacity_slider.value)))
        except ValueError as exc:
            return str(exc)
        entry.isovalue, entry.show = isovalue, show
        entry.positive_colour, entry.negative_colour = positive, negative
        entry.opacity, entry.visible = opacity, bool(self.visible_switch.value)
        return None

    def apply_slice_fields(self) -> Optional[str]:
        """Copy every slice widget into the slice state; an error message or None."""
        try:
            cube = self.slice_cube_selection.value
            if not cube:
                raise ValueError("Load a cube first; there is nothing to slice.")
            label = self.slice_axis_selection.value
            if label not in AXIS_LABELS:
                raise ValueError(f"Unknown slice axis '{label}'.")
            axis = AXIS_LABELS.index(label)
            moved = cube != self.slice_cube or axis != self.slice_axis
            index = None if moved else int(round(float(self.slice_slider.value)))
            range_label = self.slice_range_selection.value
            if range_label not in RANGE_LABELS:
                raise ValueError(f"Unknown colour range '{range_label}'.")
            self.set_slice(cube=cube, axis=axis, index=index, range_mode=RANGE_LABELS[range_label],
                           visible=bool(self.slice_switch.value))
        except ValueError as exc:
            return str(exc)
        return None

    # ------------------------------------------------------------------
    # State -> widgets
    # ------------------------------------------------------------------
    def set_status(self, message: str):
        if self.status_label is not None:
            _assign(self.status_label, "text", message)
        setter = getattr(self.viewer, "set_status_message", None)
        if setter is not None and message != self._last_status:
            setter(message)
        self._last_status = message

    def sync_widgets(self, selected: Optional[str] = None):
        """Bring the table and selections up to date with the state, writing only what
        changed (see _assign); handlers are muted meanwhile."""
        if self.cube_list is None:
            return
        self._syncing = True
        try:
            names = self.names()
            self.cube_list.show(names, [e.table_row() for e in self.entries], [e.visible for e in self.entries])
            for selection in (self.iso_cube_selection, self.slice_cube_selection,
                              self.op_a_selection, self.op_b_selection):
                if _item_names(selection) != names:
                    current = selection.value
                    selection.items = names
                    if current in names:
                        _assign(selection, "value", current)
            if selected in names:
                _assign(self.iso_cube_selection, "value", selected)
            if self.slice_cube in names:
                _assign(self.slice_cube_selection, "value", self.slice_cube)
            self._load_isosurface_fields()
            self._load_slice_fields()
        finally:
            self._syncing = False

    def _load_isosurface_fields(self):
        try:
            entry = self.entry(self.iso_cube_selection.value)
        except ValueError:
            return
        was = self._syncing
        self._syncing = True
        try:
            _assign(self.isovalue_input, "value", f"{entry.isovalue:g}")
            _assign(self.show_selection, "value", entry.show)
            _assign(self.positive_colour_input, "value", entry.positive_colour)
            _assign(self.negative_colour_input, "value", entry.negative_colour)
            _assign(self.opacity_slider, "value", float(entry.opacity))
            _assign(self.visible_switch, "value", entry.visible)
        finally:
            self._syncing = was

    def _load_slice_fields(self):
        was = self._syncing
        self._syncing = True
        try:
            n = self.slice_plane_count()
            _assign(self.slice_axis_selection, "value", AXIS_LABELS[self.slice_axis])
            _assign(self.slice_range_selection, "value",
                    next(k for k, v in RANGE_LABELS.items() if v == self.slice_range))
            _assign(self.slice_switch, "value", self.slice_visible)
            _assign(self.slice_slider, "max", float(max(n - 1, 1)))
            _assign(self.slice_slider, "tick_count", max(n, 2))
            _assign(self.slice_slider, "value", float(self.slice_index))
            _assign(self.slice_slider, "enabled", n > 1)
            _assign(self.slice_index_label, "text", f"Plane {self.slice_index + 1} of {n}" if n else "No cube")
        finally:
            self._syncing = was

    # ------------------------------------------------------------------
    # Dialog helpers
    # ------------------------------------------------------------------
    def _window(self):
        return self.viewer.main_window

    async def _show_error(self, title: str, message: str):
        await self._window().dialog(toga.ErrorDialog(title, message))

    async def _ask(self, title: str, message: str) -> bool:
        return bool(await self._window().dialog(toga.QuestionDialog(title, message)))

    async def _ask_save_path(self, title: str, suggested: str, file_types):
        return await self._window().dialog(toga.SaveFileDialog(
            title=title, suggested_filename=suggested, file_types=file_types))

    # ------------------------------------------------------------------
    # Handlers
    # ------------------------------------------------------------------
    async def on_load(self, widget):
        selected = await self._window().dialog(toga.OpenFileDialog(
            "Open Gaussian cube file(s)", multiple_select=True, file_types=CUBE_FILE_TYPES))
        if not selected:
            return
        paths = list(selected) if isinstance(selected, (list, tuple)) else [selected]
        last, notes = None, []
        for path in paths:
            self.set_status(f"Reading {Path(str(path)).name}...")
            try:
                entries = await asyncio.to_thread(self.read_entries, str(path))
            except (ValueError, MemoryError) as exc:
                await self._show_error("Cube File Error", str(exc))
                continue
            if not await self._settle_structure(entries[0]):
                continue
            notes += self.add_entries(entries)
            last = entries[-1].name
        if last is None:
            self.set_status("No cube was loaded.")
            return
        self.sync_widgets(selected=last)
        await self.refresh_scene()
        if notes:
            self.set_status(" ".join(notes))
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)

    async def _settle_structure(self, entry: CubeEntry) -> bool:
        """The first cube decides the displayed structure (asking before replacing one)."""
        if self.entries or self.structure_volume is not None:
            return True
        has_atoms = bool(len(entry.volume.atomic_numbers))
        if not self.viewer.frames:
            if not has_atoms:
                await self._show_error("No Structure", f"{entry.name} lists no atoms and no molecule is "
                                                       "displayed. Load the molecule's XYZ file first.")
                return False
            self.adopt_structure(entry)
            return True
        if self.displays_cube_molecule(entry):
            self.structure_volume = entry.volume  # already on screen: nothing to ask or reload
            return True
        shown = len(self.viewer.get_current_frame_data())
        if has_atoms and await self._ask(
                "Show the cube's molecule?",
                f"The 3D viewer shows another structure ({shown} atoms); "
                f"{entry.name} stores {len(entry.volume.atomic_numbers)} atoms.\n\n"
                "Yes: replace the displayed structure with the cube's atoms.\n"
                "No: keep the displayed structure; the volume is drawn in the cube's own coordinates."):
            self.adopt_structure(entry)
        return True

    async def on_remove(self, widget):
        name = self.iso_cube_selection.value
        if not name:
            return
        entry = self.entry(name)
        if not entry.saved and not await self._ask(
                "Remove an unsaved cube?", f"{name} has not been saved. Remove it anyway?"):
            return
        self.remove(name)
        self.sync_widgets()
        await self.refresh_scene()

    def on_cube_toggled(self, name: str, shown: bool):
        """A Show box in the cube list: draw or hide that cube (others are untouched)."""
        if self._syncing:
            return
        try:
            entry = self.entry(name)
        except ValueError:
            return
        if entry.visible == shown:
            return
        entry.visible = shown
        # The Isosurface section's Visible switch is the same setting for the cube it edits.
        if self.visible_switch is not None and self.iso_cube_selection is not None \
                and self.iso_cube_selection.value == name:
            was = self._syncing
            self._syncing = True
            try:
                _assign(self.visible_switch, "value", shown)
            finally:
                self._syncing = was
        self._start(self.refresh_scene())

    async def on_display(self, widget):
        """Reopen the 3D window (or bring it forward) and redraw every ticked cube."""
        self.set_status("Opening the 3D viewer...")
        ready, error = await self.viewer.ensure_viewer_window()
        if not ready:
            await self._show_error("The 3D viewer could not be opened", error)
            return
        await self.refresh_scene()

    # --- canvas background (a view setting: no cube is re-meshed or hidden) ---
    def apply_background(self, rgb):
        hex_colour = rgb_to_hex(rgb)
        self.viewer.set_background_color(tuple(float(c) for c in rgb))
        if self.background_input is not None:
            _assign(self.background_input, "value", hex_colour)
        if self.background_swatch is not None and hex_colour != self._swatch_hex:
            self.background_swatch.style.background_color = hex_colour
            self._swatch_hex = hex_colour

    def on_background_confirm(self, widget, **kwargs):
        """A colour typed in the field (name or #rrggbb) and Enter."""
        try:
            rgb = vd.parse_colour(widget.value)
        except ValueError as exc:
            self.apply_background(self.viewer.background_color)  # show the colour in use again
            self._start(self._show_error("Background colour", str(exc)))
            return
        self.apply_background(rgb)

    def on_background_pick(self, widget, **kwargs):
        """The Windows colour picker, opened on the colour in use."""
        try:
            owner = self._window()._impl.native
        except Exception:
            owner = None
        try:
            rgb = pick_colour(tuple(float(c) for c in self.viewer.background_color), owner)
        except RuntimeError as exc:
            self._start(self._show_error("Background colour", str(exc)))
            return
        if rgb is not None:
            self.apply_background(rgb)

    def on_iso_cube_change(self, widget, **kwargs):
        if not self._syncing:
            self._load_isosurface_fields()

    # The widget handlers below are deliberately *synchronous*. Toga runs an
    # async handler as a task that starts only after the widget write that
    # fired it has returned - by then sync_widgets has cleared _syncing - so
    # the tab's own updates used to come back as "user changes" and redo the
    # whole sync and scene (measured: two syncs and two scene builds per cube
    # load). Checking the guard here, at the moment Toga fires the handler,
    # and only then starting the async work (_start) keeps them out.
    def _start(self, coroutine):
        """Run a handler's async work on the UI loop, keeping the task referenced until done."""
        loop = getattr(self.viewer, "loop", None) or asyncio.get_event_loop()
        task = loop.create_task(coroutine)
        self._tasks.add(task)
        task.add_done_callback(self._tasks.discard)
        return task

    def on_isosurface_change(self, widget, **kwargs):
        if self._syncing or not self.entries:
            return
        self._start(self._isosurface_changed())

    async def _isosurface_changed(self):
        error = self.apply_isosurface_fields()
        if error:
            self._load_isosurface_fields()  # put the last good values back
            await self._show_error("Isosurface", error)
            return
        self.sync_widgets()
        await self.refresh_scene()

    def on_show_molecule_change(self, widget, **kwargs):
        self.set_show_structure(widget.value)

    def on_show_volumes_change(self, widget, **kwargs):
        if self._syncing:
            return
        self.show_volumes = bool(widget.value)
        self._start(self.refresh_scene())

    def on_slice_change(self, widget, **kwargs):
        if self._syncing or not self.entries:
            return
        self._start(self._slice_changed())

    async def _slice_changed(self):
        error = self.apply_slice_fields()
        self._load_slice_fields()
        if error:
            await self._show_error("Slice", error)
            return
        await self.refresh_scene()

    async def on_export_slice(self, widget):
        try:
            text = self.slice_text()
        except ValueError as exc:
            await self._show_error("Export Slice", str(exc))
            return
        path = await self._ask_save_path("Export slice data", self.default_slice_filename(), ["dat", "txt"])
        if not path:
            return
        try:
            Path(str(path)).write_text(text, encoding="utf-8")
        except OSError as exc:
            await self._show_error("Export Slice", f"Could not write {path}: {exc}")
            return
        self.set_status(f"Slice data saved: {path}")

    async def on_compute(self, widget):
        args = (self.op_a_selection.value, self.op_selection.value, self.op_b_selection.value,
                self.factor_input.value, self.result_name_input.value)
        try:
            volume, name, notes = await asyncio.to_thread(self.compute_operation, *args)
        except (ValueError, MemoryError) as exc:
            await self._show_error("Grid Operation", str(exc))
            return
        entry = self.add_result(volume, name)
        self.sync_widgets(selected=entry.name)
        await self.refresh_scene()
        self.set_status(" ".join([f"{entry.name} added (not saved yet)."] + notes))

    async def on_save_cube(self, widget):
        name = self.iso_cube_selection.value
        if not name:
            await self._show_error("Save Cube", "Load or compute a cube first.")
            return
        path = await self._ask_save_path("Save cube file", default_cube_filename(name), CUBE_FILE_TYPES)
        if not path:
            return
        try:
            written = await asyncio.to_thread(self.save_entry, name, str(path))
        except ValueError as exc:
            await self._show_error("Save Cube", str(exc))
            return
        self.sync_widgets()
        self.set_status(f"{name} saved: {written}")

    async def on_export_image(self, widget):
        """Capture the 3D window as shown, then ask for the file name (the ESP tool's order)."""
        from molecularViewer import save_canvas_image  # molecularViewer imports this module
        if not self.viewer.viewer_window_is_open():
            await self._show_error("No Image", "Open the 3D view first; the image is taken from "
                                               "that window exactly as shown.")
            return
        try:
            image = await self.viewer.capture_frame()
        except RuntimeError as exc:
            await self._show_error("Export Image", str(exc))
            return
        stem = default_cube_filename(self.iso_cube_selection.value or "cube")[:-len(".cube")]
        path = await self._ask_save_path("Export the 3D view as an image", f"{stem}.png", IMAGE_FILE_TYPES)
        if not path:
            return
        try:
            saved = save_canvas_image(image, str(path))
        except (ValueError, OSError) as exc:
            await self._show_error("Export Image", str(exc))
            return
        height, width = image.shape[:2]
        self.set_status(f"Image saved ({width} x {height} px): {saved}")

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------
    def build_tab(self) -> toga.Box:
        v = self.viewer
        field_width = v.FIELD_WIDTH
        tab = toga.Box(style=Pack(direction=COLUMN, margin=12))
        tab.add(v._hint("Gaussian cube files (.cube/.cub) as isosurfaces and slices. "
                        "Type a value and press Enter to apply it."))

        # Cube files
        tab.add(v._section_heading("Cube files"))
        tab.add(v._form_row(
            toga.Button("Load cube(s)...", on_press=self.on_load, style=Pack(flex=1, margin=(0, 8, 0, 0))),
            toga.Button("Remove selected", on_press=self.on_remove, style=Pack(flex=1, margin=(0, 8, 0, 0))),
            toga.Button("Display", on_press=self.on_display, style=Pack(flex=1)),
        ))
        self.cube_list = CubeList(self.on_cube_toggled)
        tab.add(self.cube_list.container)
        # Toga labels do not wrap: keep each hint line short enough for the window.
        tab.add(v._hint("Tick Show to draw a cube; ticked cubes are drawn together. "
                        "Display reopens the 3D window."))
        tab.add(v._hint("Remove selected, Isosurface and Save act on the cube chosen under Isosurface > Cube."))
        self.show_molecule_switch = toga.Switch("Show molecule", value=True, on_change=self.on_show_molecule_change,
                                                style=Pack(margin=(0, 16, 0, 0)))
        self.show_volumes_switch = toga.Switch("Show volumetric data", value=True,
                                               on_change=self.on_show_volumes_change)
        tab.add(v._form_row(self.show_molecule_switch, self.show_volumes_switch))
        current = rgb_to_hex(self.viewer.background_color)
        self.background_swatch = toga.Box(style=Pack(width=24, height=18, margin=(0, 8, 0, 0),
                                                     background_color=current))
        self._swatch_hex = current
        self.background_input = toga.TextInput(value=current, on_confirm=self.on_background_confirm,
                                               style=Pack(width=field_width, margin=(0, 8, 0, 0)))
        tab.add(v._form_row(v._form_label("Background:"), self.background_swatch, self.background_input,
                            toga.Button("Background colour...", on_press=self.on_background_pick,
                                        style=Pack(width=150))))

        # Isosurface
        tab.add(v._section_heading("Isosurface"))
        self.iso_cube_selection = toga.Selection(items=[], on_change=self.on_iso_cube_change, style=Pack(flex=1))
        tab.add(v._form_row(v._form_label("Cube:"), self.iso_cube_selection))
        self.isovalue_input = toga.TextInput(on_confirm=self.on_isosurface_change,
                                             style=Pack(width=field_width, margin=(0, 16, 0, 0)))
        self.show_selection = toga.Selection(items=list(SHOW_MODES), on_change=self.on_isosurface_change,
                                             style=Pack(width=110))
        tab.add(v._form_row(v._form_label("Isovalue:"), self.isovalue_input,
                            v._form_label("Show:", width=50), self.show_selection))
        self.positive_colour_input = toga.TextInput(value=DEFAULT_POSITIVE_COLOUR,
                                                    on_confirm=self.on_isosurface_change,
                                                    style=Pack(width=field_width, margin=(0, 16, 0, 0)))
        self.negative_colour_input = toga.TextInput(value=DEFAULT_NEGATIVE_COLOUR,
                                                    on_confirm=self.on_isosurface_change,
                                                    style=Pack(width=field_width))
        tab.add(v._form_row(v._form_label("Positive colour:"), self.positive_colour_input,
                            v._form_label("Negative colour:", width=104), self.negative_colour_input))
        tab.add(v._hint("Colours: " + ", ".join(vd.COLOUR_PRESETS) + ", or #rrggbb. The negative "
                        "surface is drawn at minus the isovalue."))
        self.opacity_slider = toga.Slider(min=MIN_OPACITY, max=1.0, value=1.0,
                                          on_release=self.on_isosurface_change,
                                          style=Pack(flex=1, margin=(0, 16, 0, 0)))
        self.visible_switch = toga.Switch("Visible", value=True, on_change=self.on_isosurface_change)
        tab.add(v._form_row(v._form_label("Opacity:"), self.opacity_slider, self.visible_switch))

        # Slice
        tab.add(v._section_heading("Slice"))
        self.slice_cube_selection = toga.Selection(items=[], on_change=self.on_slice_change, style=Pack(flex=1))
        tab.add(v._form_row(v._form_label("Cube:"), self.slice_cube_selection))
        self.slice_axis_selection = toga.Selection(items=list(AXIS_LABELS), value=AXIS_LABELS[DEFAULT_SLICE_AXIS],
                                                   on_change=self.on_slice_change,
                                                   style=Pack(width=field_width, margin=(0, 16, 0, 0)))
        self.slice_range_selection = toga.Selection(items=list(RANGE_LABELS), on_change=self.on_slice_change,
                                                    style=Pack(width=150))
        tab.add(v._form_row(v._form_label("Grid axis:"), self.slice_axis_selection,
                            v._form_label("Colours:", width=60), self.slice_range_selection))
        self.slice_slider = toga.Slider(min=0, max=1, value=0, on_change=self.on_slice_change,
                                        enabled=False, style=Pack(flex=1, margin=(0, 8, 0, 0)))
        self.slice_index_label = toga.Label("No cube", style=Pack(width=110, text_align=LEFT))
        tab.add(v._form_row(v._form_label("Plane:"), self.slice_slider, self.slice_index_label))
        self.slice_switch = toga.Switch("Show slice", value=False, on_change=self.on_slice_change,
                                        style=Pack(margin=(0, 16, 0, 0)))
        tab.add(v._form_row(self.slice_switch,
                            toga.Button("Export slice data...", on_press=self.on_export_slice, style=Pack(flex=1))))
        tab.add(v._hint("A slice is a plane of the cube's own grid; a, b and c are its axes (x, y and z "
                        "for an orthogonal grid). Symmetric colours run blue (−), white (0), red (+); "
                        "min to max runs red, green, blue."))

        # Grid operations
        tab.add(v._section_heading("Grid operations"))
        self.op_a_selection = toga.Selection(items=[], style=Pack(flex=1, margin=(0, 8, 0, 0)))
        self.op_selection = toga.Selection(items=list(OPERATION_LABELS), style=Pack(width=80, margin=(0, 8, 0, 0)))
        self.op_b_selection = toga.Selection(items=[], style=Pack(flex=1))
        tab.add(v._form_row(v._form_label("A:", width=20), self.op_a_selection, self.op_selection,
                            v._form_label("B:", width=20), self.op_b_selection))
        self.factor_input = toga.TextInput(value="1", style=Pack(width=60, margin=(0, 16, 0, 0)))
        self.result_name_input = toga.TextInput(placeholder="automatic", style=Pack(flex=1, margin=(0, 8, 0, 0)))
        tab.add(v._form_row(v._form_label("Factor c:", width=60), self.factor_input,
                            v._form_label("Result name:", width=86), self.result_name_input,
                            toga.Button("Compute", on_press=self.on_compute, style=Pack(width=90))))
        tab.add(v._hint("A + B and A − B need identical grids (same points, origin and axes). "
                        "The result stays in memory, marked unsaved, until you save it."))

        # Export
        tab.add(v._section_heading("Export"))
        tab.add(v._form_row(
            toga.Button("Save selected cube...", on_press=self.on_save_cube, style=Pack(flex=1, margin=(0, 8, 0, 0))),
            toga.Button("Export image...", on_press=self.on_export_image, style=Pack(flex=1)),
        ))
        self.status_label = toga.Label(" ", style=Pack(font_size=9, margin=(4, 0, 0, 0)))
        tab.add(self.status_label)

        tab.add(v._section_heading("Help"))
        tab.add(toga.MultilineTextInput(value=HelpGqteaWin.help_cube_tab, readonly=True,
                                        style=Pack(height=220, flex=1)))
        return tab
