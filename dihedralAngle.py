# Refactored dihedral angle analysis module for GQTEA at 2025 November
# Combines dihedral calculation, distribution function, and free energy computation
import os

import numpy as np
import toga
from toga.style import Pack
from toga.constants import LEFT  # for label text alignment

from help import HelpGqteaWin
from framesCounter import FramesCounter
from displayPlots import DisplayPlots
from analysisCommon import GeometryAnalysisBase


class DihedralAngleAnalyser(GeometryAnalysisBase, FramesCounter, DisplayPlots):
    """Analyze the i-j-k-l dihedral angle from an XYZ/CPMD trajectory.

    The shared machinery (trajectory reader, text status reporting, guarded
    writes, minimum image, histogram, statistics, summary layout) lives in
    :class:`analysisCommon.GeometryAnalysisBase`, shared with bond.py and
    bondAngle.py. Only what is specific to a dihedral is here.

    Unlike a distance (4*pi*r^2 dr) or an angle (sin(theta) d(theta)), a
    dihedral has a **uniform** volume element: d(phi) carries no geometric
    weighting. There is therefore no Jacobian correction to apply and no switch
    for one -- ``-RT ln P(phi)`` is already the potential of mean force.
    """

    tool_title = "DIHEDRAL ANGLE ANALYSIS"

    def _reset_analysis_state(self) -> None:
        """Clear analysis results so failed runs cannot reuse stale data."""
        self.dihedral = []
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
        self.min_dihedral_value = None
        self.elmt1 = None
        self.elmt2 = None
        self.elmt3 = None
        self.elmt4 = None
        self.cell_lengths = None
        self.wrap_180 = False
        self.smooth_fe = False
        self.hist_lower = 0.0
        self.hist_upper = 360.0
        self.fts = None
        self.fdist = None
        self.ffe = None
        self.fsummary = None

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

        self.atom_labels = await self.read_atom_labels(self.textInput_atom_labels, 4)
        if self.atom_labels is None:
            return False

        self.bin_width = await self.read_value(
            self.textInput_bin_width, "bin width (deg)", float
        )
        if self.bin_width is None:
            return False

        if self.max_angle <= 0:
            return await self.fail("Maximum angle must be greater than zero.")
        if self.max_angle > 360.0:
            return await self.fail("Maximum angle cannot exceed 360 degrees.")
        if self.dt <= 0:
            return await self.fail("Time step must be greater than zero.")
        if self.sampling <= 0:
            return await self.fail("Sampling interval must be a positive integer.")
        if self.sim_temp <= 0:
            return await self.fail("Simulation temperature must be greater than zero.")
        if self.bin_width <= 0:
            return await self.fail("Histogram bin width must be greater than zero.")

        self.wrap_180 = bool(
            getattr(self, "switch_wrap_180", None) and self.switch_wrap_180.value
        )
        self.smooth_fe = bool(
            getattr(self, "switch_smooth_fe", None) and self.switch_smooth_fe.value
        )

        # The histogram must span the range the angles actually live in. When
        # wrapping is on the values are in (-180, 180], so a [0, max_angle]
        # histogram would give every negative angle a negative bin index and
        # silently drop half the data.
        if self.wrap_180:
            self.hist_lower, self.hist_upper = -180.0, 180.0
        else:
            self.hist_lower, self.hist_upper = 0.0, self.max_angle

        if self.bin_width >= (self.hist_upper - self.hist_lower):
            return await self.fail(
                "Histogram bin width must be smaller than the histogram range "
                f"({self.hist_upper - self.hist_lower:g} degrees)."
            )

        pbc_text = "none" if self.cell_lengths is None else str(self.cell_lengths)
        self.multi_line_text.value = (
            f"Maximum dihedral (deg) --> {self.max_angle}\n"
            f"Time step (a.u.)       --> {self.dt}\n"
            f"Sampling (frames)      --> {self.sampling}\n"
            f"Temperature (K)        --> {self.sim_temp}\n"
            f"Atom labels (i j k l)  --> {self.atom_labels}\n"
            f"Bin width (deg)        --> {self.bin_width}\n"
            f"Cell lengths (Å)       --> {pbc_text}\n"
            f"Wrap to [-180, 180]    --> {self.wrap_180}\n"
            f"Histogram range (deg)  --> [{self.hist_lower:g}, {self.hist_upper:g}]\n"
            f"Smooth free energy     --> {self.smooth_fe}\n"
            f"Show plots at the end  --> {self.switch_show_plots.value}\n"
            f"Save CSV outputs       --> {self.switch_save_csv.value}\n"
        )
        return True

    async def dihedral_angle(self, widget) -> bool:
        """Compute the i-j-k-l dihedral per frame and write the time series."""
        self._ensure_output_dir()
        self.dihedral = []

        idx1, idx2, idx3, idx4 = self.atom_labels
        rows = [idx1 - 1, idx2 - 1, idx3 - 1, idx4 - 1]

        def status(done, total):
            return (
                f"Calculating dihedral angles ...\n"
                f"  processed {done} / {total} frames"
            )

        self._set_status(status(0, getattr(self, "total_frame_number", 0)))
        frames = await self.read_selected_atoms(rows, status=status)
        if frames is None:
            return False

        self.elmt1, self.elmt2, self.elmt3, self.elmt4 = (
            frames[0][n][0] for n in range(4)
        )

        # Vectorised stable atan2 formulation over all frames at once. Each of
        # the three connecting vectors goes through the minimum image, so a
        # quadruplet split across a periodic boundary still gives the real
        # geometry.
        coords = np.array(
            [[a[1:4], b[1:4], c[1:4], d[1:4]] for a, b, c, d in frames], dtype=float
        )
        ab = self.minimum_image(coords[:, 1, :] - coords[:, 0, :])
        bc = self.minimum_image(coords[:, 2, :] - coords[:, 1, :])
        cd = self.minimum_image(coords[:, 3, :] - coords[:, 2, :])

        n1 = np.cross(ab, bc)
        n2 = np.cross(bc, cd)
        nbc = np.linalg.norm(bc, axis=1)
        n1_norm = np.linalg.norm(n1, axis=1)
        n2_norm = np.linalg.norm(n2, axis=1)

        # Collinear atoms define no plane, so the dihedral is undefined there.
        valid = (nbc > 0.0) & (n1_norm > 0.0) & (n2_norm > 0.0)
        self.skipped_degenerate = int((~valid).sum())
        if not valid.any():
            return await self.fail(
                "Every frame has a degenerate geometry (collinear atoms), so no "
                "dihedral plane is defined."
            )

        n1u = n1[valid] / n1_norm[valid, None]
        n2u = n2[valid] / n2_norm[valid, None]
        bcu = bc[valid] / nbc[valid, None]

        x = (n1u * n2u).sum(axis=1)
        y = (np.cross(n1u, bcu) * n2u).sum(axis=1)
        angles = np.degrees(np.arctan2(y, x))  # (-180, 180]

        if not self.wrap_180:
            angles = np.where(angles < 0.0, angles + 360.0, angles)  # [0, 360)

        times = (np.arange(len(frames)) * self._time_increment_ps())[valid]
        self.dihedral = [[float(t), float(v)] for t, v in zip(times, angles)]
        self.frames_used = len(self.dihedral)
        self._set_status(status(self.frames_used, getattr(self, "total_frame_number", 0)))

        if self.cell_lengths is not None:
            longest = float(max(
                np.linalg.norm(ab[valid], axis=1).max(),
                nbc[valid].max(),
                np.linalg.norm(cd[valid], axis=1).max(),
            ))
            if longest > self.min_image_limit:
                await self.warning_function(
                    f"A bond of the dihedral is {longest:.3f} Å long, more than half "
                    f"the smallest cell length ({self.min_image_limit:.3f} Å). The "
                    f"minimum-image convention may not give the true geometry there."
                )

        tag = (f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_"
               f"{self.elmt3}{idx3}_{self.elmt4}{idx4}")
        self.fts = os.path.join(self.output_dir, f"dihedral_{tag}.dat")
        self.fdist = os.path.join(self.output_dir, f"dihedral_distribution_{tag}.dat")
        self.ffe = os.path.join(self.output_dir, f"dihedral_free_energy_{tag}.dat")
        self.fsummary = os.path.join(self.output_dir, f"summary_{tag}.txt")

        body = "".join(f"{t:>12.6f}{v:>12.6f}\n" for t, v in self.dihedral)
        if not await self._write_text(self.fts, body):
            return False

        # Figures go to the shared interactive viewer only; no static PNGs.
        self.save_plots(
            1, [t for t, _ in self.dihedral], [v for _, v in self.dihedral],
            "Simulation time (ps)", "Dihedral angle (°)", "Dihedral Angle Over Time",
            save_png=False
        )
        return True

    async def distribution_function(self) -> bool:
        """Histogram (%) of dihedral angles and basic statistics."""
        if not self.dihedral:
            return await self.fail("No dihedral angles computed.")

        self._set_status("Building the dihedral distribution ...")
        values = np.asarray([v for _, v in self.dihedral], dtype=float)

        self.bin_centers, self.histogram, total, self.out_of_range_count = \
            self.histogram_percent(
                values, self.bin_width, self.hist_lower, self.hist_upper
            )
        self.num_bins = len(self.bin_centers)

        if total == 0:
            return await self.fail(
                "The histogram is empty. Check the maximum angle and the "
                "'Wrap dihedral to [-180, 180]' setting."
            )
        self.histogram_sum = total

        if self.out_of_range_count:
            pct = 100.0 * self.out_of_range_count / values.size
            await self.warning_function(
                f"{self.out_of_range_count} of {values.size} dihedral angles "
                f"({pct:.1f}%) are out of range (outside "
                f"[{self.hist_lower:g}, {self.hist_upper:g}]°) and were excluded from "
                f"the distribution and the free energy."
            )

        stats = self.basic_stats(values)
        self.stats = {
            "average_dihedral": stats["average"],
            "variance": stats["variance"],
            "std_dev": stats["std_dev"],
            "largest_dihedral": stats["largest"],
            "smallest_dihedral": stats["smallest"],
            "bin_centers": self.bin_centers,
            "histogram_percentage": self.histogram,
        }

        body = "".join(
            f"{center:>12.6f}{pct:>12.6f}\n"
            for center, pct in zip(self.bin_centers, self.histogram)
        )
        if not await self._write_text(self.fdist, body):
            return False

        self.save_plots(
            2, self.bin_centers, self.histogram,
            "Dihedral angle (°)", "Dihedral distribution (%)",
            "Dihedral Angle Distribution", save_png=False
        )
        return True

    async def dihedral_free_energy(self) -> bool:
        """Compute free energy from the histogram (kcal/mol).

        ``-RT ln P(phi)``, with no Jacobian: a dihedral has a uniform measure,
        so this is already the potential of mean force. Defined up to an
        additive constant.
        """
        if not getattr(self, "histogram_sum", 0):
            self.stats = {}
            return False

        self._set_status("Computing the free energy ...")
        self.free_energy_pairs = self.free_energy_from_histogram(
            self.bin_centers, self.histogram, jacobian=None
        )

        if not self.free_energy_pairs:
            self.stats = {}
            return await self.fail(
                "Free energy could not be computed from an empty distribution."
            )

        xs = [x for x, _ in self.free_energy_pairs]
        ys = [y for _, y in self.free_energy_pairs]
        self.min_y = float(min(ys))
        self.min_dihedral_value = float(xs[ys.index(self.min_y)])
        self.min_y_idx = self.min_dihedral_value  # kept for backwards compatibility

        body = "".join(f"{xi:>12.6f}{fe:>12.6f}\n" for xi, fe in self.free_energy_pairs)
        if not await self._write_text(self.ffe, body):
            return False

        self.save_plots(
            3, xs, ys, "Dihedral angle (°)", "Free energy (kcal/mol)",
            "Dihedral Free Energy", save_png=False
        )
        return True

    async def export_csv(self) -> bool:
        """Optional CSV export with headers."""
        idx1, idx2, idx3, idx4 = self.atom_labels
        tag = (f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_"
               f"{self.elmt3}{idx3}_{self.elmt4}{idx4}")

        body = "time_ps,dihedral_deg\n" + "".join(
            f"{t:.6f},{v:.6f}\n" for t, v in self.dihedral
        )
        if not await self._write_text(
            os.path.join(self.output_dir, f"dihedral_{tag}.csv"), body
        ):
            return False

        body = "dihedral_deg,distribution_pct\n" + "".join(
            f"{x:.6f},{pct:.6f}\n" for x, pct in zip(self.bin_centers, self.histogram)
        )
        if not await self._write_text(
            os.path.join(self.output_dir, f"dihedral_distribution_{tag}.csv"), body
        ):
            return False

        if self.free_energy_pairs:
            body = "dihedral_deg,free_energy_kcal_mol\n" + "".join(
                f"{x:.6f},{fe:.6f}\n" for x, fe in self.free_energy_pairs
            )
            if not await self._write_text(
                os.path.join(self.output_dir, f"dihedral_free_energy_{tag}.csv"), body
            ):
                return False
        return True

    async def save_summary(self) -> bool:
        """Write the closing summary: the parameters used and the results."""
        if not getattr(self, "stats", None):
            return False

        idx1, idx2, idx3, idx4 = self.atom_labels
        pbc_text = "none" if self.cell_lengths is None else \
            " ".join(f"{v:g}" for v in self.cell_lengths)

        results = [
            ("Largest dihedral", f"{self.stats['largest_dihedral']:.4f} degrees"),
            ("Smallest dihedral", f"{self.stats['smallest_dihedral']:.4f} degrees"),
            ("Average dihedral", f"{self.stats['average_dihedral']:.4f} degrees"),
            ("Variance", f"{self.stats['variance']:.6f}"),
            ("Std deviation", f"{self.stats['std_dev']:.6f}"),
            ("Lowest free energy", f"{self.min_y:.4f} kcal/mol at "
                                   f"{self.min_dihedral_value:.2f} degrees"),
            ("Excluded angles", f"{self.out_of_range_count} frame(s) (out of range)"),
        ]
        if self.skipped_degenerate:
            results.append(("Degenerate frames", f"{self.skipped_degenerate} skipped"))

        summary = self.render_summary(
            f"{self.elmt1}{idx1}-{self.elmt2}{idx2}-"
            f"{self.elmt3}{idx3}-{self.elmt4}{idx4}",
            parameters=[
                ("Maximum angle", f"{self.max_angle} degrees"),
                ("Bin width", f"{self.bin_width} degrees"),
                ("Time step", f"{self.dt} a.u."),
                ("Sampling interval", f"{self.sampling} frames"),
                ("Temperature", f"{self.sim_temp} K"),
                ("Cell (a b c)", pbc_text),
                ("Wrap to [-180,180]", "yes" if self.wrap_180 else "no"),
                ("Histogram range", f"[{self.hist_lower:g}, {self.hist_upper:g}] degrees"),
                ("Jacobian", "not applicable (uniform measure)"),
                ("Smoothed FE", "yes" if self.smooth_fe else "no"),
            ],
            results=results,
        )
        if not await self._write_text(self.fsummary, summary):
            return False
        self.multi_line_text.value = summary
        return True


class DihedralUI(DihedralAngleAnalyser):
    """UI for the Dihedral Angle module (harmonized with Bond/BondAngle)."""

    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget) -> None:
        self.main_window = toga.Window(
            title="Dihedral Angle Analysis from TRAJEC.xyz",
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
        box_1.add(toga.Label("Dihedral Angle Analysis", style=heading_style))

        input_fields = [
            ("Maximum angle for Distribution (deg):", "Enter e.g. 360", "textInput_max_angle"),
            ("Simulation Time Step (a.u.):", "Enter time step", "textInput_time_step"),
            ("Sampling Interval (frames):", "Enter sampling interval", "textInput_sampling_interval"),
            ("Simulation Temperature (K):", "Enter temperature", "textInput_temperature"),
            ("Atom Labels (i j k l):", "Enter four labels, e.g., 1 2 3 4", "textInput_atom_labels"),
            ("Histogram Bin Width (deg):", "Enter e.g. 1.0", "textInput_bin_width"),
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
        self.switch_wrap_180 = toga.Switch("Wrap dihedral to [-180, 180]", value=False, style=Pack(margin=(4, 10, 4, 10)))
        self.switch_smooth_fe = toga.Switch("Smooth free energy", value=True, style=Pack(margin=(4, 10, 4, 10)))
        switches_row.add(self.switch_show_plots)
        switches_row.add(self.switch_save_csv)
        switches_row.add(self.switch_wrap_180)
        switches_row.add(self.switch_smooth_fe)
        main_box.add(switches_row)

        self.multi_line_text = toga.MultilineTextInput(style=Pack(flex=1, margin=(10, 0), font_size=12))
        self.multi_line_text.value = HelpGqteaWin.help_dihedral_angle
        main_box.add(self.multi_line_text)

        button_row = toga.Box(style=Pack(direction="row", margin=(10, 0, 0, 0)))
        self.btn_execute = toga.Button("Analyze", style=button_style, on_press=self.workflow)
        self.btn_close = toga.Button("Close", style=button_style, on_press=self.closeTopLevel)
        button_row.add(self.btn_execute)
        button_row.add(self.btn_close)
        main_box.add(button_row)

        self.main_window.content = main_box
        self.main_window.show()

    async def workflow(self, widget) -> None:
        self._reset_analysis_state()
        if not await self.read_params(widget):
            return
        if not await self.dihedral_angle(widget):
            return
        if not await self.distribution_function():
            return
        if not await self.dihedral_free_energy():
            return
        await self.save_summary()

        if self.switch_save_csv.value:
            await self.export_csv()
        if self.switch_show_plots.value:
            self.display_plots()
        # save_summary has already replaced the status text with the closing
        # summary of the parameters used and the results obtained.

    def closeTopLevel(self, widget) -> None:
        self.main_window.close()
