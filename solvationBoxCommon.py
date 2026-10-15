"""Shared pre-optimization, system-property and construction-view plumbing
for the three solvation-box builders.

The builders differ in geometry -- one cell with one solute
(single_solute_solvent_box.py), one cell with a two-solvent mixture
(mixtureSolventBox.py), two stacked regions sharing a face (sharedWallBox.py)
-- but the features added on top of their packing are the same. They live
here once, as two mixins:

* :class:`SolvationExtrasEngine` -- headless (never touches a widget, runs in a
  worker thread): progress throttling, the pre-optimization driver, system
  properties and their summary lines. Geometry enters only through hooks the
  builders override: :meth:`~SolvationExtrasEngine.packing_vdw_scale`,
  :meth:`~SolvationExtrasEngine.preopt_solute_clearance`,
  :meth:`~SolvationExtrasEngine.preopt_wall_padding`,
  :meth:`~SolvationExtrasEngine.preopt_wall_bounds`,
  :meth:`~SolvationExtrasEngine.preopt_slab_bounds`,
  :meth:`~SolvationExtrasEngine.box_atoms` and
  :meth:`~SolvationExtrasEngine.net_charge`.
* :class:`SolvationExtrasUI` -- the Toga half: the Charge and spin,
  Pre-optimization and Advanced sections, the Construction and Summary tabs,
  the run bar, progress handling on the Toga thread, the build pipeline
  (pack -> pre-optimize -> properties) and *Update properties*. Builders supply
  :meth:`~SolvationExtrasUI.read_charges` and the charge-field list.

A builder's own packing code (``build_box``) reports progress dicts through
:meth:`SolvationExtrasEngine._report`::

    {"stage": "insertion", "attempts", "max_attempts", "accepted", "target",
     "molecules": [(kind, elements, coords), ...],
     optional "label", "base", "span" (percent of the packing bar)}
    {"stage": "preopt", "iteration", "max_iterations", "energy", "molecules"}

The existing packing algorithms are deliberately *not* moved here; whether
the three builders should share a packing core is assessed separately
(docs/superpowers/specs/2026-10-15-solvation-box-consolidation-assessment.md).
"""
import asyncio
import textwrap
import time

import toga
from toga.style import Pack
from toga.style.pack import CENTER, COLUMN, LEFT, ROW

from boxConstructionView import ROTATION_STEP, ConstructionSnapshot, ConstructionView
from periodicPreOptimizer import (
    DEFAULT_MAX_ITERATIONS, MODE_FLEXIBLE, MODE_RIGID, PackedSystem, PeriodicCell, pre_optimize,
)
from solventPacking import ENGINE_AUTO, ENGINE_LEGACY, ENGINE_NUMBA, ENGINE_NUMPY, numba_available
from systemProperties import compute_properties, parse_integer, parse_multiplicity

# Seconds between progress reports from the worker thread. The construction
# view repaints at its own, lower rate; this only bounds cross-thread traffic.
PROGRESS_REPORT_INTERVAL = 0.1
# Share of the progress bar given to packing when a pre-optimization follows.
PACKING_PROGRESS_SHARE = 70.0
# Default contact scale for the pre-optimization. Packing already enforces
# d >= vdw_scale * (r_i + r_j), so relaxing to the *same* scale would find
# nothing to do; the pre-optimization pushes contacts out to this larger
# scale (or the packing scale, if that is larger).
DEFAULT_PREOPT_CONTACT_SCALE = 0.90


class SolvationExtrasEngine:
    """Headless half. Call :meth:`init_solvation_extras` from ``__init__``."""

    def init_solvation_extras(self):
        self.packing_engine = ENGINE_AUTO
        self.engine_used = ""
        self.build_seconds = 0.0

        self.preoptimize = False
        self.preopt_mode = MODE_RIGID
        self.preopt_max_iterations = DEFAULT_MAX_ITERATIONS
        self.preopt_contact_scale = None   # None -> effective_preopt_contact_scale() default
        self.preopt_result = None
        self.preopt_error = ""
        self.preopt_data = []
        self.preopt_seconds = 0.0

        self.multiplicity = None           # required; never inferred
        self.system_properties = None

        # (kind, elements, (n, 3) coordinates) per placed molecule, in the
        # order of the box's atoms -- the pre-optimizer's connectivity (never
        # guessed across a periodic boundary) and the construction view's input.
        self.molecules = []
        self.cancel_requested = False
        self._last_report = 0.0

    # -- progress ------------------------------------------------------
    def _report(self, progress, payload, force=False):
        if progress is None:
            return
        now = time.monotonic()
        if not force and now - self._last_report < PROGRESS_REPORT_INTERVAL:
            return
        self._last_report = now
        try:
            progress(payload)
        except Exception:
            pass  # reporting must never abort a build

    # -- hooks (geometry-specific) -------------------------------------
    def packing_vdw_scale(self) -> float:
        return float(self.vdw_scale)

    def preopt_solute_clearance(self) -> float:
        return 0.0

    def preopt_wall_padding(self) -> float:
        return max(0.0, float(getattr(self, "spacing", 0.0) or 0.0))

    def preopt_wall_bounds(self, system):
        return None

    def preopt_slab_bounds(self, system):
        return None

    def box_atoms(self):
        return self.data

    def net_charge(self) -> int:
        raise NotImplementedError

    # -- pre-optimization ------------------------------------------------
    def effective_preopt_contact_scale(self) -> float:
        if self.preopt_contact_scale is not None:
            return float(self.preopt_contact_scale)
        return max(DEFAULT_PREOPT_CONTACT_SCALE, self.packing_vdw_scale())

    def packed_system(self):
        return PackedSystem.from_molecules(self.molecules, self.get_vdw_radius, self.get_atomic_mass)

    def preopt_molecules(self, system, coordinates):
        return [
            (kind, system.elements[start:stop], coordinates[start:stop])
            for kind, (start, stop) in zip(system.molecule_kinds, system.molecule_slices)
        ]

    def run_preoptimization(self, progress=None):
        """Relax the packed box (worker thread). A setup problem -- e.g. a box
        too small for the minimum image -- is recorded in ``preopt_error``
        instead of raised, so the packed box is still written."""
        self.preopt_result = None
        self.preopt_error = ""
        self.preopt_data = []
        if not self.molecules:
            self.preopt_error = "there are no molecules to pre-optimize"
            return
        system = self.packed_system()
        cell = PeriodicCell(self.lattice_vectors, periodic=self.use_periodic_minimum_image)

        def callback(wrapped, iteration, energy):
            self._report(progress, {
                "stage": "preopt",
                "iteration": iteration,
                "max_iterations": self.preopt_max_iterations,
                "energy": energy,
                "molecules": self.preopt_molecules(system, wrapped),
            })

        started = time.perf_counter()
        try:
            result = pre_optimize(
                system, cell, self.effective_preopt_contact_scale(),
                mode=self.preopt_mode,
                max_iterations=self.preopt_max_iterations,
                solute_clearance=self.preopt_solute_clearance(),
                wall_padding=self.preopt_wall_padding(),
                wall_bounds=None if cell.periodic else self.preopt_wall_bounds(system),
                slab_bounds=self.preopt_slab_bounds(system),
                callback=callback,
                should_stop=lambda: self.cancel_requested,
            )
        except ValueError as error:
            self.preopt_error = str(error)
            return
        finally:
            self.preopt_seconds = time.perf_counter() - started
        self.preopt_result = result
        self.preopt_data = [
            [element, float(x), float(y), float(z)]
            for element, (x, y, z) in zip(system.elements, result.coordinates)
        ]

    # -- properties ------------------------------------------------------
    def compute_system_properties(self):
        self.system_properties = compute_properties(
            [atom[0] for atom in self.box_atoms()], self.net_charge(), self.multiplicity
        )
        return self.system_properties

    def extras_summary_lines(self, preopt_path):
        lines = []
        if self.system_properties is not None:
            lines.append("\n")
            lines.extend(self.system_properties.summary_lines())
        if self.preoptimize:
            lines.append("\n")
            if self.preopt_result is not None:
                lines.append(
                    f"Target contact scale: {self.effective_preopt_contact_scale():.3f} "
                    f"(packing vdW scale {self.packing_vdw_scale():.3f})\n"
                )
                lines.extend(self.preopt_result.summary_lines())
                lines.append(f"Pre-optimization time: {self.preopt_seconds:.2f} s\n")
                lines.append(f"Saved pre-optimized XYZ: {preopt_path}\n")
            else:
                lines.append(f"Pre-optimization not performed: {self.preopt_error}\n")
        return lines

    @staticmethod
    def write_xyz(path, atoms, comment):
        with open(path, "w") as f:
            f.write(f"{len(atoms)}\n")
            f.write(f"{comment}\n")
            for atom in atoms:
                f.write(f"{atom[0]:<3s}{atom[1]:>14.7f}{atom[2]:>14.7f}{atom[3]:>14.7f}\n")

    def write_preopt_xyz(self, path):
        a, b, c = self.lattice_vectors
        self.write_xyz(
            path, self.preopt_data,
            f"Pre-optimized solvation box ({self.preopt_result.mode}, molecules wrapped whole "
            f"into the cell); cell a b c = {a:.6f} {b:.6f} {c:.6f} Angstrom",
        )


class SolvationExtrasUI:
    """Toga half; mix in *before* the builder's logic class."""

    ENGINE_LABELS = {
        "Automatic (numba if installed)": ENGINE_AUTO,
        "NumPy cell list": ENGINE_NUMPY,
        "Numba kernel": ENGINE_NUMBA,
        "Legacy (reference, slow)": ENGINE_LEGACY,
    }
    PREOPT_MODE_LABELS = {
        "Rigid molecules": MODE_RIGID,
        "Flexible molecules": MODE_FLEXIBLE,
    }
    SETUP_TAB, CONSTRUCTION_TAB, SUMMARY_TAB = 0, 1, 2
    HINT_COLOR = "#666666"
    # (label, placeholder, attribute) for each charge field -- set per builder.
    CHARGE_FIELDS = ()
    # Where the multiplicity field is, for the message shown when it is missing
    # or impossible -- a long Setup tab can hide it below the fold.
    MULTIPLICITY_LOCATION = "Setup tab, Charge and spin section"

    # -- layout helpers ----------------------------------------------------
    def _extras_styles(self):
        return {
            "section": Pack(font_size=12, font_weight="bold", text_align=LEFT, margin=(8, 0, 4, 0)),
            "label": Pack(margin=(0, 6, 0, 0), text_align=LEFT, width=170),
            "input": Pack(flex=1, margin=0),
            "row": Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0)),
        }

    def _extras_field(self, parent, label_text, placeholder, attr_name):
        styles = self._extras_styles()
        row = toga.Box(style=styles["row"])
        text_input = toga.TextInput(placeholder=placeholder, style=styles["input"])
        setattr(self, attr_name, text_input)
        row.add(toga.Label(label_text, style=styles["label"]))
        row.add(text_input)
        parent.add(row)

    def _hint(self, text):
        # Toga labels do not wrap on their own; break the text to the window width.
        return toga.Label(textwrap.fill(text, 125),
                          style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 4, 0)))

    def _two_columns(self, parent):
        row = toga.Box(style=Pack(direction=ROW, margin=(0, 0, 4, 0)))
        left = toga.Box(style=Pack(direction=COLUMN, flex=1, margin=(0, 8, 0, 0)))
        right = toga.Box(style=Pack(direction=COLUMN, flex=1))
        row.add(left)
        row.add(right)
        parent.add(row)
        return left, right

    # -- sections ----------------------------------------------------------
    def build_charge_spin_section(self, parent, hint):
        parent.add(toga.Label("Charge and spin", style=self._extras_styles()["section"]))
        left, right = self._two_columns(parent)
        for index, (label, placeholder, attr) in enumerate(self.CHARGE_FIELDS):
            self._extras_field(left if index % 2 == 0 else right, label, placeholder, attr)
        self.build_multiplicity_field(right)
        parent.add(self._hint(hint))

    def build_multiplicity_field(self, parent):
        self._extras_field(parent, "Total spin multiplicity:", "Required: 2S+1, e.g. 1", "textInput_multiplicity")

    def build_preopt_section(self, parent):
        parent.add(toga.Label("Pre-optimization", style=self._extras_styles()["section"]))
        row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0)))
        self.switch_preoptimize = toga.Switch(
            "Pre-optimize after building", style=Pack(margin=(0, 12, 0, 0)), on_change=self.update_preopt_controls
        )
        self.selection_preopt_mode = toga.Selection(
            items=list(self.PREOPT_MODE_LABELS), style=Pack(width=180, margin=(0, 12, 0, 0))
        )
        row.add(self.switch_preoptimize)
        row.add(self.selection_preopt_mode)
        parent.add(row)
        left, right = self._two_columns(parent)
        self._extras_field(left, "Target contact scale:",
                           f"Default {DEFAULT_PREOPT_CONTACT_SCALE:.2f} or the vdW scale",
                           "textInput_preopt_contact_scale")
        self._extras_field(right, "Max iterations:", f"Default {DEFAULT_MAX_ITERATIONS}",
                           "textInput_preopt_iterations")
        parent.add(self._hint(self.preopt_hint()))

    def preopt_hint(self):
        return (
            "Pushes intermolecular contacts out to d >= scale x (r_i + r_j) under periodic boundaries "
            "(minimum image on, otherwise inside soft walls). Packing already guarantees the vdW scale, so "
            "the target must be larger to have an effect. Rigid keeps every molecule's geometry exactly; "
            "flexible lets bonds and angles relax around the input geometry. Molecules are wrapped whole "
            "into the cell."
        )

    def build_engine_section(self, parent):
        styles = self._extras_styles()
        parent.add(toga.Label("Advanced", style=styles["section"]))
        row = toga.Box(style=styles["row"])
        row.add(toga.Label("Packing engine:", style=styles["label"]))
        self.selection_engine = toga.Selection(items=list(self.ENGINE_LABELS), style=Pack(width=260))
        if not numba_available():
            self.selection_engine.items = [k for k, v in self.ENGINE_LABELS.items() if v != ENGINE_NUMBA]
        row.add(self.selection_engine)
        parent.add(row)
        parent.add(self._hint(
            "All engines produce the same box for the same seed; Legacy is the original, slower loop."
        ))

    def build_construction_tab(self):
        page = toga.Box(style=Pack(direction=COLUMN, margin=(6, 6, 6, 6), flex=1))
        self.construction_view = ConstructionView(width=770, height=430)
        page.add(self.construction_view.widget)
        view_row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(6, 0, 0, 0)))
        view = self.construction_view
        step = ROTATION_STEP
        for text, action in [
            ("Turn left", lambda w: view.rotate(horizontal_degrees=-step)),
            ("Turn right", lambda w: view.rotate(horizontal_degrees=step)),
            ("Tilt up", lambda w: view.rotate(vertical_degrees=-step)),
            ("Tilt down", lambda w: view.rotate(vertical_degrees=step)),
            ("Zoom +", lambda w: view.zoom(1.25)),
            ("Zoom -", lambda w: view.zoom(0.8)),
            ("Reset view", lambda w: view.reset_view()),
        ]:
            view_row.add(toga.Button(text, on_press=action, style=Pack(width=92, margin=(0, 6, 0, 0))))
        page.add(view_row)
        return page

    def build_summary_tab(self, intro_text):
        page = toga.Box(style=Pack(direction=COLUMN, margin=(6, 6, 6, 6), flex=1))
        self.multi_line_text.value = intro_text
        page.add(self.multi_line_text)
        row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(6, 0, 0, 0)))
        self.btn_update_properties = toga.Button(
            "Update properties", on_press=self.update_properties, style=Pack(width=150, margin=(0, 8, 0, 0))
        )
        row.add(self.btn_update_properties)
        row.add(toga.Label("Re-checks charge and multiplicity on the last box without rebuilding it.",
                           style=Pack(color=self.HINT_COLOR)))
        page.add(row)
        return page

    def build_run_bar(self, start_label):
        run_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(8, 0, 6, 0)))
        self.start_button = toga.Button(start_label, style=Pack(margin=(0, 8, 0, 0), width=180),
                                        on_press=self.workflow)
        self.btn_cancel = toga.Button("Cancel", style=Pack(margin=(0, 8, 0, 0), width=120), on_press=self.cancel_run)
        self.btn_cancel.enabled = False
        self.btn_close = toga.Button("Close", style=Pack(width=120), on_press=self.closeTopLevel)
        run_box.add(self.start_button)
        run_box.add(self.btn_cancel)
        run_box.add(self.btn_close)
        return run_box

    def build_tabs(self, setup_box, intro_text):
        self.tabs = toga.OptionContainer(
            content=[
                ("Setup", toga.ScrollContainer(content=setup_box, horizontal=False)),
                ("Construction", self.build_construction_tab()),
                ("Summary", self.build_summary_tab(intro_text)),
            ],
            style=Pack(flex=1),
        )
        return self.tabs

    # -- control state -----------------------------------------------------
    def update_preopt_controls(self, widget=None):
        enabled = bool(self.switch_preoptimize.value) and not getattr(self, "is_running", False)
        self.selection_preopt_mode.enabled = enabled
        self.textInput_preopt_iterations.enabled = enabled
        self.textInput_preopt_contact_scale.enabled = enabled

    def set_extras_running(self, is_running):
        """Enable/disable the controls this mixin owns."""
        self.is_running = is_running
        self.btn_cancel.enabled = is_running
        self.btn_update_properties.enabled = not is_running
        self.switch_preoptimize.enabled = not is_running
        self.selection_engine.enabled = not is_running
        self.textInput_multiplicity.enabled = not is_running
        if is_running:
            self.cancel_requested = False
        self.update_preopt_controls()

    def cancel_run(self, widget):
        self.cancel_requested = True
        self.btn_cancel.enabled = False
        self.status_label.text = "Cancel requested. Stopping after the current step..."

    # -- inputs ------------------------------------------------------------
    def read_charges(self):
        """Parse the builder's charge fields into attributes; raise ValueError."""
        raise NotImplementedError

    async def read_extras_params(self) -> bool:
        self.preoptimize = bool(self.switch_preoptimize.value)
        self.preopt_mode = self.PREOPT_MODE_LABELS.get(self.selection_preopt_mode.value, MODE_RIGID)
        self.packing_engine = self.ENGINE_LABELS.get(self.selection_engine.value, ENGINE_AUTO)
        try:
            iterations = parse_integer(
                self.textInput_preopt_iterations.value, "Maximum pre-optimization iterations",
                default=DEFAULT_MAX_ITERATIONS,
            )
            if iterations <= 0:
                raise ValueError("Maximum pre-optimization iterations must be a positive integer.")
            self.preopt_max_iterations = iterations
            contact_text = self.textInput_preopt_contact_scale.value.strip()
            if contact_text:
                try:
                    contact = float(contact_text)
                except ValueError:
                    raise ValueError(f"Target contact scale must be a number (got '{contact_text}').")
                if contact <= 0.0:
                    raise ValueError("Target contact scale must be positive.")
                self.preopt_contact_scale = contact
            else:
                self.preopt_contact_scale = None
            self.read_charges()
        except ValueError as error:
            await self.main_window.dialog(toga.ErrorDialog("Error", str(error)))
            return False
        return await self.read_multiplicity()

    async def read_multiplicity(self) -> bool:
        try:
            self.multiplicity = parse_multiplicity(self.textInput_multiplicity.value)
        except ValueError as error:
            await self.point_to_multiplicity(str(error))
            return False
        return True

    async def point_to_multiplicity(self, message):
        """Show the Setup tab, say where the field is, then put the cursor in it
        (focusing also scrolls the Setup page to it)."""
        try:
            self.tabs.current_tab = self.SETUP_TAB
        except Exception:
            pass
        await self.main_window.dialog(toga.ErrorDialog(
            "Spin multiplicity",
            f"{message}\n\nEnter it in 'Total spin multiplicity' ({self.MULTIPLICITY_LOCATION}).",
        ))
        try:
            self.textInput_multiplicity.focus()
        except Exception:
            pass

    def extras_parameter_text(self, width=32):
        text = f"{'Total spin multiplicity:':<{width}} {self.multiplicity}\n"
        if self.preoptimize:
            text += (
                f"{'Pre-optimization:':<{width}} {self.preopt_mode}, contact scale "
                f"{self.effective_preopt_contact_scale():.3f}, at most {self.preopt_max_iterations} iterations\n"
            )
        return text

    # -- progress (always on the Toga thread) ------------------------------
    def _packing_share(self):
        return PACKING_PROGRESS_SHARE if self.preoptimize else 100.0

    def _snapshot(self, stage, molecules, caption=""):
        return ConstructionSnapshot(stage, tuple(self.lattice_vectors), molecules, caption)

    def _apply_progress(self, payload):
        if payload.get("stage") == "insertion":
            attempts, accepted, target = payload["attempts"], payload["accepted"], payload["target"]
            max_attempts = max(1, payload.get("max_attempts", getattr(self, "max_attempts", 1)))
            local = max(attempts / max_attempts, 1.0 if target <= 0 else accepted / max(1, target))
            local = max(0.0, min(1.0, local))
            base, span = payload.get("base", 0.0), payload.get("span", 100.0)
            self.progress_bar.max = 100
            self.progress_bar.value = self._packing_share() * (base + span * local) / 100.0
            label = payload.get("label", "Insertion")
            self.status_label.text = (
                f"{label}: attempts {attempts}/{max_attempts} | "
                f"accepted {accepted}/{target} solvent molecules"
            )
            self.construction_view.submit(self._snapshot(
                "Inserting solvent", payload["molecules"], f"{label.lower()}: {accepted}/{target} molecules"
                if "label" in payload else f"{accepted}/{target} molecules"
            ))
        elif payload.get("stage") == "preopt":
            iteration, limit = payload["iteration"], max(1, payload["max_iterations"])
            share = self._packing_share()
            self.progress_bar.value = share + (100.0 - share) * min(1.0, iteration / limit)
            self.status_label.text = f"Pre-optimizing: iteration {iteration} | penalty {payload['energy']:.4g}"
            self.construction_view.submit(self._snapshot(
                "Pre-optimizing", payload["molecules"], f"iteration {iteration}"
            ))

    # -- pipeline ----------------------------------------------------------
    async def run_construction(self, initial_molecules, packed_caption):
        """Pack (worker thread), optionally pre-optimize, compute properties.

        Returns False if the run was cancelled (nothing must be written)."""
        loop = asyncio.get_running_loop()

        def progress(payload):
            loop.call_soon_threadsafe(self._apply_progress, payload)

        self.preopt_result = None
        self.preopt_error = ""
        self.system_properties = None
        self.tabs.current_tab = self.CONSTRUCTION_TAB
        self.construction_view.submit(
            self._snapshot("Solute placed" if initial_molecules else "Empty box", initial_molecules), force=True
        )

        await asyncio.to_thread(self.build_box, progress)
        if self.cancel_requested:
            self.report_cancelled()
            return False
        self.construction_view.submit(self._snapshot("Packed box", list(self.molecules), packed_caption()),
                                      force=True)

        if self.preoptimize:
            self.status_label.text = "Pre-optimizing..."
            await asyncio.to_thread(self.run_preoptimization, progress)
            if self.cancel_requested:
                self.report_cancelled()
                return False
            if self.preopt_result is not None:
                system = self.packed_system()
                self.construction_view.submit(self._snapshot(
                    "Pre-optimized box (molecules wrapped whole)",
                    self.preopt_molecules(system, self.preopt_result.coordinates),
                    f"{self.preopt_result.report_final['clashes']} clashing pairs left",
                ), force=True)

        self.compute_system_properties()
        return True

    def report_cancelled(self):
        self.status_label.text = "Cancelled"
        self.multi_line_text.value += "\nRun cancelled: no output files were written.\n"

    async def report_problems(self):
        if self.preoptimize and self.preopt_result is None and self.preopt_error:
            await self.main_window.dialog(toga.InfoDialog(
                "Pre-optimization not performed",
                f"{self.preopt_error}\n\nThe packed box was written without pre-optimization.",
            ))
        props = self.system_properties
        if props is not None and not props.multiplicity_valid:
            await self.main_window.dialog(toga.InfoDialog(
                "Check the spin multiplicity",
                f"Multiplicity {props.multiplicity} is not possible here: {props.message}.\n\n"
                f"Enter a valid value in 'Total spin multiplicity' ({self.MULTIPLICITY_LOCATION}) and "
                "press 'Update properties' on the Summary tab; the box does not need to be rebuilt.",
            ))

    async def update_properties(self, widget):
        if not self.box_atoms():
            await self.main_window.dialog(toga.InfoDialog("No box yet", "Build a solvent box first."))
            return
        try:
            self.read_charges()
        except ValueError as error:
            await self.main_window.dialog(toga.ErrorDialog("Error", str(error)))
            return
        if not await self.read_multiplicity():
            return
        try:
            self.compute_system_properties()
            self.multi_line_text.value = self.write_summary()
        except ValueError as error:
            await self.main_window.dialog(toga.ErrorDialog("Error", str(error)))
            return
        except OSError as error:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Error saving summary TXT file: {error}"))
            return
        await self.report_problems()
