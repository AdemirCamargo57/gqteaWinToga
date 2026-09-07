# Refactored bond angle analysis module for GQTEA at 2025 November
# Combines bond angle calculation, distribution function, and free energy computation
import math
import os

import numpy as np
import toga
from toga.style import Pack
from toga.constants import LEFT  # for label text alignment

from help import HelpGqteaWin
from framesCounter import FramesCounter
from displayPlots import DisplayPlots
from analysisCommon import GeometryAnalysisBase


class BondAngleAnalyser(GeometryAnalysisBase, FramesCounter, DisplayPlots):
    """Analyze the i-j-k bond angle (vertex at j) from an XYZ/CPMD trajectory.

    The shared machinery (trajectory reader, text status reporting, guarded
    writes, minimum image, histogram, statistics, summary layout) lives in
    :class:`analysisCommon.GeometryAnalysisBase`, shared with bond.py and
    dihedralAngle.py. Only what is specific to an angle is here.
    """

    tool_title = "BOND ANGLE ANALYSIS"

    def _reset_analysis_state(self) -> None:
        """Clear analysis results so failed runs cannot reuse stale data."""
        self.angles = []
        self.bin_centers = []
        self.histogram = []
        self.histogram_sum = 0
        self.out_of_range_count = 0
        self.frames_used = 0
        self.skipped_degenerate = 0
        self.stats = {}
        self.free_energy_pairs = []
        self.min_y = None
        self.min_y_idx = None
        self.min_angle_value = None
        self.elmt1 = None
        self.elmt2 = None
        self.elmt3 = None
        self.cell_lengths = None
        self.apply_jacobian = False
        self.smooth_fe = False
        self.angle_file = None
        self.fd = None
        self.fe = None
        self.fs = None

    async def read_params(self, widget) -> bool:
        """Read and validate user inputs from the UI.

        Returns True when every parameter is present and valid; otherwise shows
        an error dialog and returns False. Callers must check the return value
        before running the calculation.
        """
        if not self.require_trajectory():
            return await self.fail(
                "Please select an input TRAJEC.xyz file with the Browse button "
                "before running the analysis."
            )

        ok, self.cell_lengths = await self.read_cell_lengths(self.textInput_cell_lengths)
        if not ok:
            return False

        self.max_angle = await self.read_value(
            self.textInput_max_angle, "max_angle (deg)", float
        )
        if self.max_angle is None:
            return False

        self.dt = await self.read_value(self.textInput_time_step, "time step (a.u.)", float)
        if self.dt is None:
            return False

        self.sampling = await self.read_value(
            self.textInput_sampling_interval, "sampling (frames)", int
        )
        if self.sampling is None:
            return False

        self.sim_temp = await self.read_value(
            self.textInput_temperature, "simulation temperature (K)", float
        )
        if self.sim_temp is None:
            return False

        self.atom_labels = await self.read_atom_labels(self.textInput_atom_labels, 3)
        if self.atom_labels is None:
            return False

        self.bin_width = await self.read_value(
            self.textInput_bin_width, "bin width (deg)", float
        )
        if self.bin_width is None:
            return False

        if self.max_angle <= 0:
            return await self.fail("Maximum angle must be greater than zero.")
        if self.max_angle > 180.0:
            return await self.fail("Maximum angle cannot exceed 180 degrees.")
        if self.dt <= 0:
            return await self.fail("Time step must be greater than zero.")
        if self.sampling <= 0:
            return await self.fail("Sampling interval must be a positive integer.")
        if self.sim_temp <= 0:
            return await self.fail("Simulation temperature must be greater than zero.")
        if self.bin_width <= 0:
            return await self.fail("Histogram bin width must be greater than zero.")
        if self.bin_width >= self.max_angle:
            return await self.fail("Histogram bin width must be smaller than maximum angle.")

        self.apply_jacobian = bool(
            getattr(self, "switch_use_jacobian", None) and self.switch_use_jacobian.value
        )
        self.smooth_fe = bool(
            getattr(self, "switch_smooth_fe", None) and self.switch_smooth_fe.value
        )

        pbc_text = "none" if self.cell_lengths is None else str(self.cell_lengths)
        self.multi_line_text.value = (
            f"Maximum angle (deg)  --> {self.max_angle}\n"
            f"Time step (a.u.)     --> {self.dt}\n"
            f"Sampling (frames)    --> {self.sampling}\n"
            f"Temperature (K)      --> {self.sim_temp}\n"
            f"Atom labels (i j k)  --> {self.atom_labels}\n"
            f"Bin width (deg)      --> {self.bin_width}\n"
            f"Cell lengths (Å)     --> {pbc_text}\n"
            f"sin(theta) Jacobian  --> {self.apply_jacobian}\n"
            f"Smooth free energy   --> {self.smooth_fe}\n"
            f"Show plots at the end --> {self.switch_show_plots.value}\n"
            f"Save CSV outputs      --> {self.switch_save_csv.value}\n"
        )
        return True

    async def bond_angle(self) -> bool:
        """Compute the i-j-k angle (vertex at j) per frame and write the series."""
        self._ensure_output_dir()
        self.angles = []

        idx1, idx2, idx3 = self.atom_labels
        rows = [idx1 - 1, idx2 - 1, idx3 - 1]

        def status(done, total):
            return (
                f"Calculating bond angles ...\n"
                f"  processed {done} / {total} frames"
            )

        self._set_status(status(0, getattr(self, "total_frame_number", 0)))
        frames = await self.read_selected_atoms(rows, status=status)
        if frames is None:
            return False

        self.elmt1, self.elmt2, self.elmt3 = (frames[0][n][0] for n in range(3))

        # Vectorised: (N, 3, 3) coordinates -> N angles. The two arms run from
        # the vertex j out to i and k, each through the minimum image.
        coords = np.array([[a[1:4], b[1:4], c[1:4]] for a, b, c in frames], dtype=float)
        vec1 = self.minimum_image(coords[:, 0, :] - coords[:, 1, :])
        vec2 = self.minimum_image(coords[:, 2, :] - coords[:, 1, :])
        norm1 = np.linalg.norm(vec1, axis=1)
        norm2 = np.linalg.norm(vec2, axis=1)

        # An arm of zero length has no defined angle; drop those frames.
        valid = (norm1 > 0.0) & (norm2 > 0.0)
        self.skipped_degenerate = int((~valid).sum())
        if not valid.any():
            return await self.fail(
                "Every frame has a degenerate geometry (an atom sits on the vertex)."
            )

        cosines = np.clip(
            (vec1[valid] * vec2[valid]).sum(axis=1) / (norm1[valid] * norm2[valid]),
            -1.0, 1.0,
        )
        angles = np.degrees(np.arccos(cosines))

        times = (np.arange(len(frames)) * self._time_increment_ps())[valid]
        self.angles = [[float(t), float(v)] for t, v in zip(times, angles)]
        self.frames_used = len(self.angles)
        self._set_status(status(self.frames_used, getattr(self, "total_frame_number", 0)))

        # Minimum image is only exact for arms shorter than half the smallest
        # box length; past that the wrapped vector may not be the real one.
        if self.cell_lengths is not None:
            longest = float(max(norm1[valid].max(), norm2[valid].max()))
            if longest > self.min_image_limit:
                await self.warning_function(
                    f"An arm of the angle is {longest:.3f} Å long, more than half the "
                    f"smallest cell length ({self.min_image_limit:.3f} Å). The "
                    f"minimum-image convention may not give the true geometry there."
                )

        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}"
        self.angle_file = os.path.join(self.output_dir, f"angles_{tag}.dat")
        self.fd = os.path.join(self.output_dir, f"angles_distribution_{tag}.dat")
        self.fe = os.path.join(self.output_dir, f"angle_free_energy_{tag}.dat")
        self.fs = os.path.join(self.output_dir, f"summary_{tag}.txt")

        body = "".join(f"{t:>12.6f}{v:>12.6f}\n" for t, v in self.angles)
        if not await self._write_text(self.angle_file, body):
            return False

        # Figures go to the shared interactive viewer only; no static PNGs.
        self.save_plots(
            1, [t for t, _ in self.angles], [v for _, v in self.angles],
            "Simulation time (ps)", "Bond angle (°)", "Bond Angles over Time",
            save_png=False
        )
        return True

    async def distribution_function(self) -> bool:
        """Histogram (%) of bond angles and basic statistics."""
        if not self.angles:
            return await self.fail("No angles computed.")

        self._set_status("Building the bond angle distribution ...")
        values = np.asarray([v for _, v in self.angles], dtype=float)

        self.bin_centers, self.histogram, total, self.out_of_range_count = \
            self.histogram_percent(values, self.bin_width, 0.0, self.max_angle)
        self.num_bins = len(self.bin_centers)

        if total == 0:
            return await self.fail(
                "The histogram is empty. Increase 'Maximum angle for Distribution'.\n"
                f"Suggested minimum: {round(min(180.0, float(values.max()) * 1.1), 1)} degrees."
            )
        self.histogram_sum = total

        if self.out_of_range_count:
            # Excluded rather than clamped into the end bin: clamping silently
            # inflates the last bar and biases the free energy derived from it.
            pct = 100.0 * self.out_of_range_count / values.size
            await self.warning_function(
                f"{self.out_of_range_count} of {values.size} angles ({pct:.1f}%) are "
                f"out of range (>= max angle = {self.max_angle}°) and were excluded "
                f"from the distribution and the free energy. The largest is "
                f"{float(values.max()):.3f}°."
            )

        stats = self.basic_stats(values)
        self.stats = {
            "average_angle": stats["average"],
            "variance": stats["variance"],
            "std_dev": stats["std_dev"],
            "largest_angle": stats["largest"],
            "smallest_angle": stats["smallest"],
            "bin_centers": self.bin_centers,
            "histogram_percentage": self.histogram,
        }

        body = "".join(
            f"{center:>12.6f}{pct:>12.6f}\n"
            for center, pct in zip(self.bin_centers, self.histogram)
        )
        if not await self._write_text(self.fd, body):
            return False

        self.save_plots(
            2, self.bin_centers, self.histogram,
            "Bond angle (°)", "Bond angle distribution (%)",
            "Bond Angle Distribution", save_png=False
        )
        return True

    async def free_energy(self) -> bool:
        """Compute free energy from the histogram (kcal/mol).

        By default this is ``-RT ln P(θ)``. With the Jacobian switch on it
        becomes ``-RT ln[P(θ)/sin θ]``, which removes the sin θ volume-element
        bias of sampling an angle in three dimensions. Both are defined only up
        to an additive constant.
        """
        if not getattr(self, "histogram_sum", 0):
            self.stats = {}
            return False

        self._set_status("Computing the free energy ...")
        jacobian = (lambda deg: math.sin(math.radians(deg))) if self.apply_jacobian else None
        self.free_energy_pairs = self.free_energy_from_histogram(
            self.bin_centers, self.histogram, jacobian=jacobian
        )

        if not self.free_energy_pairs:
            self.stats = {}
            return await self.fail(
                "Free energy could not be computed from an empty distribution."
            )

        xs = [x for x, _ in self.free_energy_pairs]
        ys = [y for _, y in self.free_energy_pairs]
        self.min_y = float(min(ys))
        self.min_angle_value = float(xs[ys.index(self.min_y)])
        self.min_y_idx = self.min_angle_value  # kept for backwards compatibility

        body = "".join(f"{xi:>12.6f}{fe:>12.6f}\n" for xi, fe in self.free_energy_pairs)
        if not await self._write_text(self.fe, body):
            return False

        label = "Bond angle PMF (kcal/mol)" if self.apply_jacobian \
            else "Bond angle free energy (kcal/mol)"
        self.save_plots(
            3, xs, ys, "Bond angle (°)", label, "Bond Angle Free Energy",
            save_png=False
        )
        return True

    async def export_csv(self) -> bool:
        """Optional CSV export with headers."""
        idx1, idx2, idx3 = self.atom_labels
        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}"

        body = "time_ps,bond_angle_deg\n" + "".join(
            f"{t:.6f},{v:.6f}\n" for t, v in self.angles
        )
        if not await self._write_text(
            os.path.join(self.output_dir, f"angles_{tag}.csv"), body
        ):
            return False

        body = "bond_angle_deg,distribution_pct\n" + "".join(
            f"{x:.6f},{pct:.6f}\n" for x, pct in zip(self.bin_centers, self.histogram)
        )
        if not await self._write_text(
            os.path.join(self.output_dir, f"angles_distribution_{tag}.csv"), body
        ):
            return False

        if self.free_energy_pairs:
            column = "pmf_kcal_mol" if self.apply_jacobian else "free_energy_kcal_mol"
            body = f"bond_angle_deg,{column}\n" + "".join(
                f"{x:.6f},{fe:.6f}\n" for x, fe in self.free_energy_pairs
            )
            if not await self._write_text(
                os.path.join(self.output_dir, f"angle_free_energy_{tag}.csv"), body
            ):
                return False
        return True

    async def save_summary(self) -> bool:
        """Write the closing summary: the parameters used and the results."""
        if not getattr(self, "stats", None):
            return False

        idx1, idx2, idx3 = self.atom_labels
        pbc_text = "none" if self.cell_lengths is None else \
            " ".join(f"{v:g}" for v in self.cell_lengths)
        fe_label = "Lowest PMF" if self.apply_jacobian else "Lowest free energy"

        results = [
            ("Largest angle", f"{self.stats['largest_angle']:.4f} degrees"),
            ("Smallest angle", f"{self.stats['smallest_angle']:.4f} degrees"),
            ("Average angle", f"{self.stats['average_angle']:.4f} degrees"),
            ("Variance", f"{self.stats['variance']:.6f}"),
            ("Std deviation", f"{self.stats['std_dev']:.6f}"),
            (fe_label, f"{self.min_y:.4f} kcal/mol at "
                       f"{self.min_angle_value:.2f} degrees"),
            ("Excluded angles", f"{self.out_of_range_count} frame(s) (>= max angle)"),
        ]
        if self.skipped_degenerate:
            results.append(
                ("Degenerate frames", f"{self.skipped_degenerate} skipped")
            )

        summary = self.render_summary(
            f"{self.elmt1}{idx1}-{self.elmt2}{idx2}-{self.elmt3}{idx3}",
            parameters=[
                ("Maximum angle", f"{self.max_angle} degrees"),
                ("Bin width", f"{self.bin_width} degrees"),
                ("Time step", f"{self.dt} a.u."),
                ("Sampling interval", f"{self.sampling} frames"),
                ("Temperature", f"{self.sim_temp} K"),
                ("Cell (a b c)", pbc_text),
                ("Jacobian sin(theta)", "yes" if self.apply_jacobian else "no"),
                ("Smoothed FE", "yes" if self.smooth_fe else "no"),
            ],
            results=results,
        )
        if not await self._write_text(self.fs, summary):
            return False
        self.multi_line_text.value = summary
        return True


class BondAngleUI(BondAngleAnalyser):
    """UI for the Bond Angle module (mirrors the Bond Length UI conventions)."""

    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget) -> None:
        self.main_window = toga.Window(
            title="Bond Angle Analysis from TRAJEC.xyz",
            size=(760, 700),
        )

        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style   = Pack(margin=(5, 5), text_align=LEFT, width=260)
        input_style   = Pack(flex=1, margin=(5, 5))
        button_style  = Pack(margin=5, width=110)
        row_style     = Pack(direction="row", margin=(0, 0, 5, 0))

        main_box = toga.Box(style=Pack(direction="column", margin=20))

        # Header row. There is no frame-count label and no progress bar: loading
        # and calculation progress are reported as text in multi_line_text, which
        # ends the run showing only the summary.
        box_1 = toga.Box(style=Pack(direction="row", margin=(0, 0, 10, 0)))
        main_box.add(box_1)
        box_1.add(toga.Label("Bond Angle Analysis", style=heading_style))

        input_fields = [
            ("Maximum angle for Distribution Function (deg):", "Enter the maximum angle in (°) (suggestion: 180).", "textInput_max_angle"),
            ("Simulation Time Step (a.u.):", "Enter the time step used in the simulation in a.u. (e.g., 5).", "textInput_time_step"),
            ("Sampling Interval (frames):", "Enter sampling interval in a.u. (e.g., 5)", "textInput_sampling_interval"),
            ("Simulation Temperature (K):", "Enter the temperature used in the simulation in K (e.g., 300)", "textInput_temperature"),
            ("Atom Labels (i j k):", "Enter 3 atom labels, separated by a space (e.g., 1 2 3)", "textInput_atom_labels"),
            ("Histogram Bin Width (deg):", "Enter the bin width for the histogram (suggestion 1)", "textInput_bin_width"),
            ("Cell Lengths a b c (Å, optional):",
             "Leave blank for no periodic boundaries (e.g., 12.4 12.4 12.4)",
             "textInput_cell_lengths"),
        ]

        for label_text, placeholder, attr_name in input_fields:
            row = toga.Box(style=row_style)
            row.add(toga.Label(label_text, style=label_style))
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            row.add(text_input)
            main_box.add(row)

        file_row = toga.Box(style=row_style)
        file_row.add(toga.Label("Select Trajectory File:", style=label_style))
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select TRAJEC.xyz file", style=input_style
        )
        browse_button = toga.Button("Browse", on_press=self.frames_counter, style=button_style)
        file_row.add(self.textInput_file)
        file_row.add(browse_button)
        main_box.add(file_row)

        switches_row = toga.Box(style=row_style)
        self.switch_show_plots = toga.Switch("Show plots at the end", value=True, style=Pack(margin=(4, 10, 4, 0)))
        self.switch_save_csv = toga.Switch("Save CSV outputs", value=True, style=Pack(margin=(4, 10, 4, 10)))
        self.switch_use_jacobian = toga.Switch("Use sin(θ) Jacobian", value=False, style=Pack(margin=(4, 10, 4, 10)))
        self.switch_smooth_fe = toga.Switch("Smooth free energy", value=True, style=Pack(margin=(4, 10, 4, 10)))
        switches_row.add(self.switch_show_plots)
        switches_row.add(self.switch_save_csv)
        switches_row.add(self.switch_use_jacobian)
        switches_row.add(self.switch_smooth_fe)
        main_box.add(switches_row)

        self.multi_line_text = toga.MultilineTextInput(style=Pack(flex=1, margin=(10, 0), font_size=12))
        self.multi_line_text.value = HelpGqteaWin.help_bond_angle
        main_box.add(self.multi_line_text)

        button_row = toga.Box(style=Pack(direction="row", margin=(10, 0, 0, 0)))
        self.btn_execute = toga.Button("Analyze", style=button_style, on_press=self.workflow)
        self.btn_close = toga.Button("Close", style=button_style, on_press=self.close_bond_angle_window)
        button_row.add(self.btn_execute)
        button_row.add(self.btn_close)
        main_box.add(button_row)

        self.main_window.content = main_box
        self.main_window.show()

    async def workflow(self, widget) -> None:
        self._reset_analysis_state()
        if not await self.read_params(widget):
            return
        if not await self.bond_angle():
            return
        if not await self.distribution_function():
            return
        if not await self.free_energy():
            return
        await self.save_summary()

        if self.switch_save_csv.value:
            await self.export_csv()
        if self.switch_show_plots.value:
            self.display_plots()
        # save_summary has already replaced the status text with the closing
        # summary of the parameters used and the results obtained.

    def close_bond_angle_window(self, widget) -> None:
        self.main_window.close()
