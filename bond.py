# Refactored bond analysis module for GQTEA at 2025 November
import os

import numpy as np
import toga
from toga.style import Pack
from toga.constants import LEFT  # for label text alignment

from help import HelpGqteaWin
from framesCounter import FramesCounter
from displayPlots import DisplayPlots
from analysisCommon import GeometryAnalysisBase


class BondAnalyser(GeometryAnalysisBase, FramesCounter, DisplayPlots):
    """Numerical analysis for bond length, distribution and free energy.

    The shared machinery (trajectory reader, text status reporting, guarded
    writes, minimum image, histogram, statistics, summary layout) lives in
    :class:`analysisCommon.GeometryAnalysisBase`, which bondAngle.py and
    dihedralAngle.py use as well. Only what is specific to a distance is here.
    """

    tool_title = "BOND LENGTH ANALYSIS"

    def _reset_analysis_state(self) -> None:
        """Clear analysis results so failed runs cannot reuse stale data."""
        self.bond_lengths = []
        self.bin_centers = []
        self.histogram = []
        self.histogram_sum = 0
        self.out_of_range_count = 0
        self.frames_used = 0
        self.stats = {}
        self.free_energy_pairs = []
        self.average_bond = None
        self.min_y = None
        self.min_y_idx = None
        self.min_bond_length = None
        self.elmt1 = None
        self.elmt2 = None
        self.cell_lengths = None
        self.apply_jacobian = False
        self.smooth_fe = False
        self.fb = None
        self.fe = None
        self.fd = None
        self.fs = None

    async def read_params(self, widget) -> bool:
        """Read and validate user inputs from the UI.

        Returns True when every parameter is present and valid; otherwise shows
        an error dialog and returns False. Callers must check the return value
        before running the calculation.
        """
        # A trajectory must be loaded first: validating the atom labels needs
        # the atom count, and the status text needs the frame count.
        if not self.require_trajectory():
            return await self.fail(
                "Please select an input TRAJEC.xyz file with the Browse button "
                "before running the analysis."
            )

        # Cell lengths are optional and read first: when they are given, the
        # minimum-image scheme caps the usable maximum r at half the smallest
        # box length, which constrains the max_r check below.
        ok, self.cell_lengths = await self.read_cell_lengths(self.textInput_cell_lengths)
        if not ok:
            return False

        self.max_r = await self.read_value(self.textInput_max_r, "max_r", float)
        if self.max_r is None:
            return False

        self.dt = await self.read_value(self.textInput_time_step, "time step (a.u.)", float)
        if self.dt is None:
            return False

        # sampling is number of frames → integer
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

        self.atom_labels = await self.read_atom_labels(self.textInput_atom_labels, 2)
        if self.atom_labels is None:
            return False

        self.bin_width = await self.read_value(self.textInput_bin_width, "bin width (Å)", float)
        if self.bin_width is None:
            return False

        if self.max_r <= 0:
            return await self.fail("Maximum r must be greater than zero.")
        if self.dt <= 0:
            return await self.fail("Time step must be greater than zero.")
        if self.sampling <= 0:
            return await self.fail("Sampling interval must be a positive integer.")
        if self.sim_temp <= 0:
            return await self.fail("Simulation temperature must be greater than zero.")
        if self.bin_width <= 0:
            return await self.fail("Histogram bin width must be greater than zero.")
        if self.bin_width >= self.max_r:
            return await self.fail("Histogram bin width must be smaller than maximum r.")

        # Minimum image cannot resolve separations beyond half the smallest box
        # length, so max_r must stay inside that bound.
        if self.cell_lengths is not None and self.max_r > self.min_image_limit:
            return await self.fail(
                f"Maximum r ({self.max_r}) must not exceed half the smallest cell "
                f"length ({self.min_image_limit:.4f} Å) for the periodic "
                f"minimum-image method to be valid."
            )

        self.apply_jacobian = bool(
            getattr(self, "switch_jacobian", None) and self.switch_jacobian.value
        )
        self.smooth_fe = bool(
            getattr(self, "switch_smooth_fe", None) and self.switch_smooth_fe.value
        )

        pbc_text = "none" if self.cell_lengths is None else str(self.cell_lengths)
        self.multi_line_text.value = (
            f"Maximum r --> {self.max_r}\n"
            f"Time step (a.u.) --> {self.dt}\n"
            f"Interval sampling (frames) --> {self.sampling}\n"
            f"Simulation temperature (K) --> {self.sim_temp}\n"
            f"Selected atom labels --> {self.atom_labels}\n"
            f"Histogram bin width (Å) --> {self.bin_width}\n"
            f"Periodic cell lengths (Å) --> {pbc_text}\n"
            f"Jacobian r² correction --> {self.apply_jacobian}\n"
            f"Smooth free energy --> {self.smooth_fe}\n"
            f"Show plots at the end --> {self.switch_show_plots.value}\n"
            f"Save CSV outputs --> {self.switch_save_csv.value}\n"
        )
        return True

    async def bond_length(self, widget) -> bool:
        """Compute bond lengths per frame and write the time-series file."""
        self._ensure_output_dir()
        self.bond_lengths = []

        idx1, idx2 = self.atom_labels
        rows = [idx1 - 1, idx2 - 1]

        def status(done, total):
            return (
                f"Calculating bond lengths ...\n"
                f"  processed {done} / {total} frames"
            )

        self._set_status(status(0, getattr(self, "total_frame_number", 0)))
        frames = await self.read_selected_atoms(rows, status=status)
        if frames is None:
            return False

        self.elmt1, self.elmt2 = frames[0][0][0], frames[0][1][0]

        # Vectorised: (N, 2, 3) coordinates -> N separations.
        coords = np.array([[a[1:4], b[1:4]] for a, b in frames], dtype=float)
        delta = self.minimum_image(coords[:, 1, :] - coords[:, 0, :])
        distances = np.sqrt((delta * delta).sum(axis=1))

        times = np.arange(len(frames)) * self._time_increment_ps()
        self.bond_lengths = [[float(t), float(r)] for t, r in zip(times, distances)]
        self.frames_used = len(frames)
        self._set_status(status(self.frames_used, getattr(self, "total_frame_number", 0)))

        pair_tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}"
        self.fb = os.path.join(self.output_dir, f"bond_{pair_tag}.dat")
        self.fe = os.path.join(self.output_dir, f"free_energy_{pair_tag}.dat")
        self.fd = os.path.join(self.output_dir, f"bond_distribution_{pair_tag}.dat")
        self.fs = os.path.join(self.output_dir, f"summary_{pair_tag}.txt")

        body = "".join(f"{t:>12.6f}{r:>12.6f}\n" for t, r in self.bond_lengths)
        if not await self._write_text(self.fb, body):
            return False

        # Figures go to the shared interactive viewer only. save_png is False so
        # no static PNG image files are left next to the trajectory.
        self.save_plots(
            1, [t for t, _ in self.bond_lengths], [r for _, r in self.bond_lengths],
            "Simulation time (ps)", "Bond length (Å)", "Bond lengths", save_png=False
        )
        return True

    async def distribution_function(self) -> bool:
        """Histogram (%) of bond lengths and basic statistics."""
        if not self.bond_lengths:
            return await self.fail("No bond lengths computed.")

        self._set_status("Building the bond length distribution ...")
        bonds = np.asarray([r for _, r in self.bond_lengths], dtype=float)

        self.bin_centers, self.histogram, total, self.out_of_range_count = \
            self.histogram_percent(bonds, self.bin_width, 0.0, self.max_r)
        self.num_bins = len(self.bin_centers)

        if total == 0:
            return await self.fail(
                "The histogram is empty. Increase 'Maximum r for Distribution'.\n"
                f"Suggested minimum: {round(float(bonds.max()) * 1.1, 2)} Å."
            )
        self.histogram_sum = total

        if self.out_of_range_count:
            # The histogram is normalised over the surviving samples only, so
            # the distribution and the free energy derived from it would be
            # silently biased if this went unreported.
            pct = 100.0 * self.out_of_range_count / bonds.size
            await self.warning_function(
                f"{self.out_of_range_count} of {bonds.size} bond lengths "
                f"({pct:.1f}%) are out of range (>= max_r = {self.max_r} Å) and were "
                f"excluded from the distribution and the free energy. The longest is "
                f"{float(bonds.max()):.3f} Å."
            )

        stats = self.basic_stats(bonds)
        self.average_bond = stats["average"]
        self.stats = {
            "average_bond": stats["average"],
            "variance": stats["variance"],
            "std_dev": stats["std_dev"],
            "largest_bond": stats["largest"],
            "smallest_bond": stats["smallest"],
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
            "Bond length (Å)", "Bond length distribution (%)",
            "Bond length distribution function", save_png=False
        )
        return True

    async def free_energy(self) -> bool:
        """Compute free energy from the histogram (kcal/mol).

        By default this is ``-RT ln P(r)``. With the Jacobian switch on it
        becomes the potential of mean force ``W(r) = -RT ln[P(r)/r²]``, which
        removes the 4πr² volume-element bias of sampling a separation in three
        dimensions. Both are defined only up to an additive constant.
        """
        if not getattr(self, "histogram_sum", 0):
            self.stats = {}
            return False

        self._set_status("Computing the free energy ...")
        jacobian = (lambda r: r * r) if self.apply_jacobian else None
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
        self.min_bond_length = float(xs[ys.index(self.min_y)])
        self.min_y_idx = self.min_bond_length  # kept for backwards compatibility

        body = "".join(f"{xi:>12.6f}{fe:>12.6f}\n" for xi, fe in self.free_energy_pairs)
        if not await self._write_text(self.fe, body):
            return False

        label = "Bond potential of mean force (kcal/mol)" if self.apply_jacobian \
            else "Bond free energy (kcal/mol)"
        self.save_plots(
            3, xs, ys, "Bond length (Å)", label,
            "Bond length Free Energy Distribution", save_png=False
        )
        return True

    async def export_csv(self) -> bool:
        """Optional CSV export with headers."""
        idx1, idx2 = self.atom_labels
        pair_tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}"

        body = "time_ps,bond_length_A\n" + "".join(
            f"{t:.6f},{r:.6f}\n" for t, r in self.bond_lengths
        )
        if not await self._write_text(
            os.path.join(self.output_dir, f"bond_{pair_tag}.csv"), body
        ):
            return False

        body = "bond_length_A,distribution_pct\n" + "".join(
            f"{x:.6f},{pct:.6f}\n" for x, pct in zip(self.bin_centers, self.histogram)
        )
        if not await self._write_text(
            os.path.join(self.output_dir, f"bond_distribution_{pair_tag}.csv"), body
        ):
            return False

        if self.free_energy_pairs:
            column = "potential_of_mean_force_kcal_mol" if self.apply_jacobian \
                else "free_energy_kcal_mol"
            body = f"bond_length_A,{column}\n" + "".join(
                f"{x:.6f},{fe:.6f}\n" for x, fe in self.free_energy_pairs
            )
            if not await self._write_text(
                os.path.join(self.output_dir, f"free_energy_{pair_tag}.csv"), body
            ):
                return False
        return True

    async def save_summary(self) -> bool:
        """Write the closing summary: the parameters used and the results."""
        if not getattr(self, "stats", None):
            return False

        idx1, idx2 = self.atom_labels
        pbc_text = "none" if self.cell_lengths is None else \
            " ".join(f"{v:g}" for v in self.cell_lengths)
        fe_label = "Lowest PMF" if self.apply_jacobian else "Lowest free energy"

        summary = self.render_summary(
            f"{self.elmt1}{idx1}-{self.elmt2}{idx2}",
            parameters=[
                ("Maximum r", f"{self.max_r} Angstrom"),
                ("Bin width", f"{self.bin_width} Angstrom"),
                ("Time step", f"{self.dt} a.u."),
                ("Sampling interval", f"{self.sampling} frames"),
                ("Temperature", f"{self.sim_temp} K"),
                ("Cell (a b c)", pbc_text),
                ("Jacobian r2 (PMF)", "yes" if self.apply_jacobian else "no"),
                ("Smoothed FE", "yes" if self.smooth_fe else "no"),
            ],
            results=[
                ("Largest bond", f"{self.stats['largest_bond']:.4f} Angstrom"),
                ("Smallest bond", f"{self.stats['smallest_bond']:.4f} Angstrom"),
                ("Average bond", f"{self.stats['average_bond']:.4f} Angstrom"),
                ("Variance", f"{self.stats['variance']:.6f}"),
                ("Std deviation", f"{self.stats['std_dev']:.6f}"),
                (fe_label, f"{self.min_y:.4f} kcal/mol at "
                           f"{self.min_bond_length:.2f} Angstrom"),
                ("Excluded bonds", f"{self.out_of_range_count} frame(s) (r >= max r)"),
            ],
        )
        if not await self._write_text(self.fs, summary):
            return False
        self.multi_line_text.value = summary
        return True


class BondUI(BondAnalyser):
    """User Interface for the Bond Analysis module."""

    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget) -> None:
        self.main_window = toga.Window(
            title="Bond Length Analysis from Molecular Dynamics Simulations",
            size=(720, 700),
        )

        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style   = Pack(margin=(5, 5), text_align=LEFT, width=240)
        input_style   = Pack(flex=1, margin=(5, 5))
        button_style  = Pack(margin=5, width=110)
        row_style     = Pack(direction="row", margin=(0, 0, 5, 0))

        main_box = toga.Box(style=Pack(direction="column", margin=20))

        # Header row. There is no frame-count label and no progress bar: loading
        # and calculation progress are reported as text in multi_line_text, which
        # ends the run showing only the summary.
        box_1 = toga.Box(style=Pack(direction="row", margin=(0, 0, 10, 0)))
        main_box.add(box_1)
        box_1.add(toga.Label("Bond Length Analysis", style=heading_style))

        input_fields = [
            ("Maximum r for Distribution Function:", "Enter maximum r value", "textInput_max_r"),
            ("Simulation Time Step (a.u.):", "Enter time step", "textInput_time_step"),
            ("Sampling Interval (frames):", "Enter sampling interval", "textInput_sampling_interval"),
            ("Simulation Temperature (K):", "Enter temperature", "textInput_temperature"),
            ("Atom Labels (e.g., 2 3):",
             "Enter two atom labels, separated by a space (e.g., 1 2)", "textInput_atom_labels"),
            ("Histogram Bin Width (Å):", "Enter bin width", "textInput_bin_width"),
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
        self.switch_jacobian = toga.Switch(
            "Jacobian r² correction (PMF)", value=False, style=Pack(margin=(4, 10, 4, 10))
        )
        self.switch_smooth_fe = toga.Switch(
            "Smooth free energy", value=False, style=Pack(margin=(4, 10, 4, 10))
        )
        switches_row.add(self.switch_show_plots)
        switches_row.add(self.switch_save_csv)
        switches_row.add(self.switch_jacobian)
        switches_row.add(self.switch_smooth_fe)
        main_box.add(switches_row)

        self.multi_line_text = toga.MultilineTextInput(style=Pack(flex=1, margin=(10, 0), font_size=12))
        self.multi_line_text.value = HelpGqteaWin.help_bond_analysis
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
        if not await self.bond_length(widget):
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

    def closeTopLevel(self, widget) -> None:
        self.main_window.close()
