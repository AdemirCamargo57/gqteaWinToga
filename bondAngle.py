# Refactored bond angle analysis module for GQTEA at 2025 november
# Combines bond angle calculation, distribution function, and free energy computation
import os
import numpy as np
import toga
from statistics import variance, stdev
from toga.style import Pack
from toga.constants import LEFT  # for label text alignment

from help import HelpGqteaWin
from framesCounter import FramesCounter
from displayPlots import DisplayPlots


class BondAngleAnalyser(FramesCounter, DisplayPlots):
    """Analyze bond angles from a trajectory in XYZ/CPMD-like format."""

    atufs: float = 0.02418884326505  # Atomic time unit (fs)

    def _smoothing_window(self, n_points: int) -> int:
        """Choose a small odd-sized moving-average window for FE smoothing."""
        if n_points < 5:
            return 1

        window = max(5, min(11, n_points // 12))
        if window % 2 == 0:
            window += 1
        return min(window, n_points if n_points % 2 == 1 else n_points - 1)

    def _smooth_series(self, values: list[float]) -> list[float]:
        """Apply a gentle moving-average smoother with edge padding."""
        n_points = len(values)
        window = self._smoothing_window(n_points)
        if window <= 1:
            return values[:]

        pad = window // 2
        padded = np.pad(np.asarray(values, dtype=float), pad_width=pad, mode="edge")
        kernel = np.ones(window, dtype=float) / window
        return np.convolve(padded, kernel, mode="valid").tolist()

    async def warning_function(self, msg: str) -> None:
        await self.main_window.dialog(toga.InfoDialog("Error", f"{msg}"))

    def _ensure_output_dir(self) -> None:
        """Guarantee output_dir exists; default to TRAJEC folder or CWD."""
        if not getattr(self, "output_dir", None):
            if hasattr(self, "trajec") and self.trajec:
                self.output_dir = os.path.dirname(os.path.abspath(self.trajec)) or os.getcwd()
            else:
                self.output_dir = os.getcwd()
        os.makedirs(self.output_dir, exist_ok=True)

    def _time_increment_ps(self) -> float:
        """Time between stored frames in picoseconds."""
        # dt [a.u.] * sampling [frames] * atufs [fs/a.u.] -> fs; /1000 -> ps
        return self.dt * self.sampling * self.atufs / 1000.0

    def _update_progress(self, current: int, total: int) -> None:
        if total <= 0:
            return
        if current % 10 == 0 or current == total:
            self.progress_bar.value = (current / total) * 100

    async def read_params(self, widget) -> bool:
        """Read/validate user inputs and echo them to the multiline box."""

        async def read_input(text_input: toga.TextInput, field_name: str, expected_type):
            value = text_input.value.strip()
            if not value:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Please input a valid value for {field_name}.")
                )
                return None
            try:
                if expected_type == list:
                    labels = [int(label) for label in value.split()]
                    if len(labels) != 3:
                        raise ValueError("Please input exactly three atom labels.")
                    return labels
                else:
                    return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        self.max_angle = await read_input(self.textInput_max_angle, "max_angle (deg)", float)
        if self.max_angle is None:
            return False

        self.dt = await read_input(self.textInput_time_step, "time step (a.u.)", float)
        if self.dt is None:
            return False

        self.sampling = await read_input(self.textInput_sampling_interval, "sampling (frames)", int)
        if self.sampling is None:
            return False

        self.sim_temp = await read_input(self.textInput_temperature, "simulation temperature (K)", float)
        if self.sim_temp is None:
            return False

        self.atom_labels = await read_input(self.textInput_atom_labels, "atom labels (i j k)", list)
        if self.atom_labels is None:
            return False

        self.bin_width = await read_input(self.textInput_bin_width, "bin width (deg)", float)
        if self.bin_width is None:
            return False

        if self.max_angle <= 0:
            await self.warning_function("Maximum angle must be greater than zero.")
            return False
        if self.dt <= 0:
            await self.warning_function("Time step must be greater than zero.")
            return False
        if self.sampling <= 0:
            await self.warning_function("Sampling interval must be a positive integer.")
            return False
        if self.sim_temp <= 0:
            await self.warning_function("Simulation temperature must be greater than zero.")
            return False
        if self.bin_width <= 0:
            await self.warning_function("Histogram bin width must be greater than zero.")
            return False
        if self.bin_width >= self.max_angle:
            await self.warning_function("Histogram bin width must be smaller than maximum angle.")
            return False

        if any(label <= 0 for label in self.atom_labels):
            await self.warning_function("Atom labels must be positive integers starting at 1.")
            return False

        if getattr(self, "num_atoms", 0):
            if any(label > self.num_atoms for label in self.atom_labels):
                await self.warning_function(
                    f"Atom labels must be between 1 and {self.num_atoms}."
                )
                return False

        # UI feedback with switches
        update_text = (
            f"Maximum angle (deg)  --> {self.max_angle}\n"
            f"Time step (a.u.)     --> {self.dt}\n"
            f"Sampling (frames)    --> {self.sampling}\n"
            f"Temperature (K)      --> {self.sim_temp}\n"
            f"Atom labels (i j k)  --> {self.atom_labels}\n"
            f"Bin width (deg)      --> {self.bin_width}\n"
            f"Show plots at the end --> {self.switch_show_plots.value}\n"
            f"Save CSV outputs      --> {self.switch_save_csv.value}\n"
            f"Use sin(theta) Jacobian --> {self.switch_use_jacobian.value}\n"
        )
        self.multi_line_text.value = update_text
        return True

    async def bond_angle(self) -> None:
        """Compute bond angles (i-j-k; vertex at j) for each frame and save time series."""
        self._ensure_output_dir()
        self.angles = []

        idx1, idx2, idx3 = map(int, self.atom_labels)
        time_increment_ps = self._time_increment_ps()

        frame_number = 0
        line_number = 0

        try:
            with open(self.trajec, "r") as f:
                self.progress_bar.start()

                while True:
                    title_line = f.readline()
                    if not title_line:
                        break
                    line_number += 1

                    comment_line = f.readline()
                    if not comment_line:
                        break
                    line_number += 1

                    atom_data = []
                    for _ in range(self.num_atoms):
                        line = f.readline()
                        if not line:
                            await self.warning_function("Unexpected end of file while reading atoms.")
                            self.progress_bar.stop()
                            return
                        tokens = line.strip().split()
                        if len(tokens) < 4:
                            await self.warning_function(f"Line {line_number} has invalid atom line format.")
                            self.progress_bar.stop()
                            return
                        element = tokens[0]
                        x, y, z = map(float, tokens[1:4])
                        atom_data.append((element, x, y, z))
                        line_number += 1

                    try:
                        a = np.array(atom_data[idx1 - 1][1:4], dtype=float)
                        b = np.array(atom_data[idx2 - 1][1:4], dtype=float)  # vertex
                        c = np.array(atom_data[idx3 - 1][1:4], dtype=float)

                        elmt1 = atom_data[idx1 - 1][0]
                        elmt2 = atom_data[idx2 - 1][0]
                        elmt3 = atom_data[idx3 - 1][0]
                    except IndexError:
                        await self.warning_function(f"Atom index out of range at line {line_number}.")
                        self.progress_bar.stop()
                        return

                    vec1 = a - b
                    vec2 = c - b
                    n1 = np.linalg.norm(vec1)
                    n2 = np.linalg.norm(vec2)
                    if n1 == 0.0 or n2 == 0.0:
                        # Degenerate geometry; skip this frame
                        frame_number += 1
                        self._update_progress(frame_number, self.total_frame_number)
                        continue

                    cos_angle = float(np.dot(vec1, vec2) / (n1 * n2))
                    cos_angle = float(np.clip(cos_angle, -1.0, 1.0))
                    angle_deg = float(np.degrees(np.arccos(cos_angle)))
                    # Clip to [0, max_angle]
                    angle_deg = max(0.0, min(self.max_angle, angle_deg))

                    sim_time = frame_number * time_increment_ps
                    self.angles.append([sim_time, angle_deg])

                    frame_number += 1
                    self._update_progress(frame_number, self.total_frame_number)

                self.progress_bar.value = 100
                self.progress_bar.stop()

        except Exception as e:
            await self.warning_function(f"Error calculating bond angles: {e}")
            return

        # Save element labels for filenames
        self.elmt1, self.elmt2, self.elmt3 = elmt1, elmt2, elmt3

        # Save time series
        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}"
        self.angle_file = os.path.join(self.output_dir, f"angles_{tag}.dat")
        with open(self.angle_file, "w") as f:
            for t, ang in self.angles:
                f.write(f"{t:>12.6f}{ang:>12.6f}\n")

        # Plot
        xs = [row[0] for row in self.angles]
        ys = [row[1] for row in self.angles]
        self.save_plots(1, xs, ys, "Simulation time (ps)", "Bond angle (°)", "Bond Angles over Time")

    async def distribution_function(self) -> None:
        """Histogram (%) of bond angles and basic statistics."""
        if not self.angles:
            await self.warning_function("No angles computed.")
            return

        angles_only = [ang for _, ang in self.angles]

        self.num_bins = max(1, int(self.max_angle / self.bin_width))
        self.bin_centers = [self.bin_width * (i + 0.5) for i in range(self.num_bins)]
        counts = [0] * self.num_bins

        for ang in angles_only:
            # bin index
            idx = int(ang / self.bin_width)
            if idx == self.num_bins:  # edge case when ang == max_angle
                idx -= 1
            if 0 <= idx < self.num_bins:
                counts[idx] += 1

        total = sum(counts)
        if total == 0:
            msg = (
                "The histogram is empty. Increase 'Maximum angle for Distribution'.\n"
                f"Suggested minimum: {round(min(180.0, 2 * max(angles_only)), 1)} degrees."
            )
            await self.warning_function(msg)
            return

        histogram_pct = [(c / total) * 100.0 for c in counts]
        self.histogram_sum = total
        self.histogram = histogram_pct

        # Stats (safe for N=1)
        avg = sum(angles_only) / len(angles_only)
        largest = max(angles_only)
        smallest = min(angles_only)
        var_val = variance(angles_only) if len(angles_only) > 1 else 0.0
        std_val = stdev(angles_only) if len(angles_only) > 1 else 0.0

        self.stats = {
            "average_angle": avg,
            "variance": var_val,
            "std_dev": std_val,
            "largest_angle": largest,
            "smallest_angle": smallest,
            "bin_centers": self.bin_centers,
            "histogram_percentage": histogram_pct,
        }

        # Save distribution (per-triplet)
        idx1, idx2, idx3 = map(int, self.atom_labels)
        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}"
        self.fd = os.path.join(self.output_dir, f"angles_distribution_{tag}.dat")
        with open(self.fd, "w") as f:
            for center, pct in zip(self.bin_centers, histogram_pct):
                f.write(f"{center:>12.6f}{pct:>12.6f}\n")

        # Plot histogram
        self.save_plots(
            2, self.bin_centers, self.histogram,
            "Bond angle (°)", "Bond angle distribution (%)",
            "Bond Angle Distribution"
        )

    def free_energy(self) -> None:
        """Compute free energy (kcal/mol) from P(θ) or P(θ)/sinθ."""
        if not getattr(self, "histogram_sum", 0):
            self.stats = {}
            return

        R = 0.001987204  # kcal/(mol*K)
        T = self.sim_temp
        use_jacobian = bool(self.switch_use_jacobian.value)

        raw_probabilities = []
        for theta_deg, pct in zip(self.bin_centers, self.histogram):
            if pct <= 0.0:
                continue
            p = pct / 100.0
            if use_jacobian:
                # Divide by sin(theta) to remove geometric bias (ensure nonzero)
                p = p / max(np.sin(np.radians(theta_deg)), 1e-12)
            raw_probabilities.append([theta_deg, max(p, 1e-12)])

        if not raw_probabilities:
            self.stats = {}
            return

        xs = [theta for theta, _ in raw_probabilities]
        smoothed_probabilities = self._smooth_series([p for _, p in raw_probabilities])

        self.free_energy_pairs = []
        for theta_deg, p_smooth in zip(xs, smoothed_probabilities):
            FE = -R * T * np.log(max(p_smooth, 1e-12))
            self.free_energy_pairs.append([theta_deg, FE])

        ys = [y for _, y in self.free_energy_pairs]

        self.min_y = float(min(ys))
        self.min_y_idx = float(xs[ys.index(self.min_y)])

        # Save FE (per-triplet)
        idx1, idx2, idx3 = map(int, self.atom_labels)
        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}"
        self.fe = os.path.join(self.output_dir, f"angle_free_energy_{tag}.dat")
        with open(self.fe, "w") as f:
            for xi, FE in self.free_energy_pairs:
                f.write(f"{xi:>12.6f}{FE:>12.6f}\n")

        # Plot FE
        self.save_plots(
            3, xs, ys,
            "Bond angle (°)", "Bond angle free energy (kcal/mol)",
            "Bond Angle Free Energy"
        )

    def export_csv(self) -> None:
        """CSV export with headers for angles, distribution, and FE."""
        idx1, idx2, idx3 = map(int, self.atom_labels)
        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}"

        # Angles time series
        csv_angles = os.path.join(self.output_dir, f"angles_{tag}.csv")
        with open(csv_angles, "w") as f:
            f.write("time_ps,bond_angle_deg\n")
            for t, ang in self.angles:
                f.write(f"{t:.6f},{ang:.6f}\n")

        # Distribution
        csv_dist = os.path.join(self.output_dir, f"angles_distribution_{tag}.csv")
        with open(csv_dist, "w") as f:
            f.write("bond_angle_deg,distribution_pct\n")
            for x, pct in zip(self.bin_centers, self.histogram):
                f.write(f"{x:.6f},{pct:.6f}\n")

        # Free energy
        if hasattr(self, "free_energy_pairs") and self.free_energy_pairs:
            csv_fe = os.path.join(self.output_dir, f"angle_free_energy_{tag}.csv")
            with open(csv_fe, "w") as f:
                f.write("bond_angle_deg,free_energy_kcal_mol\n")
                for x, fe in self.free_energy_pairs:
                    f.write(f"{x:.6f},{fe:.6f}\n")

    def save_summary(self) -> None:
        """Write summary to per-triplet file and mirror it into the multiline box."""
        if not getattr(self, "stats", None):
            return

        idx1, idx2, idx3 = map(int, self.atom_labels)
        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}"
        self.fs = os.path.join(self.output_dir, f"summary_{tag}.txt")

        with open(self.fs, "w") as f:
            f.write(f"Number of frames:........... {self.total_frame_number}\n")
            f.write(f"Number of atoms:............ {self.num_atoms}\n")
            f.write(f"Selected atoms:............. {self.elmt1}{idx1}-{self.elmt2}{idx2}-{self.elmt3}{idx3}\n")
            f.write(f"The largest bond angle:..... {self.stats['largest_angle']:.4f} degrees\n")
            f.write(f"The smallest bond angle:.... {self.stats['smallest_angle']:.4f} degrees\n")
            f.write(f"Bond angle average:......... {self.stats['average_angle']:.4f} degrees\n")
            f.write(f"Bond angle variance:........ {self.stats['variance']:.4f}\n")
            f.write(f"Bond angle std deviation:... {self.stats['std_dev']:.4f}\n")
            if hasattr(self, "min_y") and hasattr(self, "min_y_idx"):
                f.write(
                    f"Lowest free energy:......... {self.min_y:.4f} kcal/mol "
                    f"at {self.min_y_idx:.2f} degrees\n"
                )

        with open(self.fs, "r") as f:
            self.multi_line_text.value = f.read()


class BondAngleUI(BondAngleAnalyser):
    """UI for the Bond Angle module (mirrors Bond UI conventions)."""

    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget) -> None:
        self.main_window = toga.Window(
            title="Bond Angle Analysis from TRAJEC.xyz",
            size=(740, 680),
        )

        # Styles (use margin, direction strings)
        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style   = Pack(margin=(5, 5), text_align=LEFT, width=260)
        input_style   = Pack(flex=1, margin=(5, 5))
        button_style  = Pack(margin=5, width=110)
        row_style     = Pack(direction="row", margin=(0, 0, 5, 0))
        col_style     = Pack(direction="column", margin=(0, 0, 5, 0))

        # Main container
        main_box = toga.Box(style=Pack(direction="column", margin=20))

        # Header row (no colored background)
        box_1 = toga.Box(style=Pack(direction="row", margin=(0, 0, 10, 0)))
        main_box.add(box_1)

        box_1a = toga.Box(style=Pack(width=460))
        box_1b = toga.Box(style=Pack(width=180))
        box_1c = toga.Box(style=Pack(width=100))
        box_1.add(box_1a)
        box_1.add(box_1b)
        box_1.add(box_1c)

        title_label = toga.Label("Bond Angle Analysis", style=heading_style)
        self.progress_label = toga.Label(" ", style=heading_style)
        self.status_label = toga.Label(" ", style=Pack(margin=(0, 0, 0, 10)))  # will print output dir

        box_1a.add(title_label)
        box_1b.add(self.progress_label)
        box_1c.add(self.status_label)

        # Inputs
        input_fields = [
            ("Maximum angle for Distribution Function (deg):", "Enter e.g. 180", "textInput_max_angle"),
            ("Simulation Time Step (a.u.):", "Enter time step", "textInput_time_step"),
            ("Sampling Interval (frames):", "Enter sampling interval", "textInput_sampling_interval"),
            ("Simulation Temperature (K):", "Enter temperature", "textInput_temperature"),
            ("Atom Labels (i j k):", "Enter three labels, e.g., 1 2 3", "textInput_atom_labels"),
            ("Histogram Bin Width (deg):", "Enter e.g. 1.0", "textInput_bin_width"),
        ]

        for label_text, placeholder, attr_name in input_fields:
            row = toga.Box(style=row_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            row.add(label)
            row.add(text_input)
            main_box.add(row)

        # File selection
        file_row = toga.Box(style=row_style)
        file_label = toga.Label("Select Trajectory File:", style=label_style)
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select TRAJEC.xyz file", style=input_style
        )
        browse_button = toga.Button("Browse", on_press=self.frames_counter, style=button_style)
        file_row.add(file_label)
        file_row.add(self.textInput_file)
        file_row.add(browse_button)
        main_box.add(file_row)

        # Switches: show plots, save csv, use Jacobian
        switches_row = toga.Box(style=row_style)
        self.switch_show_plots = toga.Switch("Show plots at the end", value=True, style=Pack(margin=(4, 10, 4, 0)))
        self.switch_save_csv = toga.Switch("Save CSV outputs", value=True, style=Pack(margin=(4, 10, 4, 10)))
        self.switch_use_jacobian = toga.Switch("Use sin(θ) Jacobian", value=False, style=Pack(margin=(4, 10, 4, 10)))
        switches_row.add(self.switch_show_plots)
        switches_row.add(self.switch_save_csv)
        switches_row.add(self.switch_use_jacobian)

        # Progress bar
        progress_col = toga.Box(style=col_style)
        self.progress_bar = toga.ProgressBar(max=100)
        progress_col.add(self.progress_bar)

        main_box.add(switches_row)
        main_box.add(progress_col)

        # Help / Output
        self.multi_line_text = toga.MultilineTextInput(style=Pack(flex=1, margin=(10, 0), font_size=12))
        self.multi_line_text.value = HelpGqteaWin.help_bond_angle
        main_box.add(self.multi_line_text)

        # Buttons
        button_row = toga.Box(style=Pack(direction="row", margin=(10, 0, 0, 0)))
        self.btn_execute = toga.Button("Analyze", style=button_style, on_press=self.workflow)
        self.btn_close = toga.Button("Close", style=button_style, on_press=self.close_bond_angle_window)
        button_row.add(self.btn_execute)
        button_row.add(self.btn_close)
        main_box.add(button_row)

        self.main_window.content = main_box
        self.main_window.show()

    async def workflow(self, widget) -> None:
        if not await self.read_params(widget):
            return
        if not hasattr(self, "trajec") or not self.trajec:
            await self.warning_function("No trajectory file selected.")
            return

        await self.bond_angle()
        await self.distribution_function()
        self.free_energy()
        self.save_summary()

        if self.switch_save_csv.value:
            self.export_csv()
        if self.switch_show_plots.value:
            self.display_plots()

        # Print directory status
        out_dir = getattr(self, "output_dir", os.getcwd())
        msg = f"Outputs saved to: {out_dir}"
        self.status_label.text = msg
        self.multi_line_text.value = (self.multi_line_text.value.rstrip() + "\n\n" + msg + "\n")

    def close_bond_angle_window(self, widget) -> None:
        self.main_window.close()


# Previous version of bondAngle.py before refactoring
# import os, toga
# import numpy as np
# from statistics import stdev, variance
# from toga.style import Pack
# from toga.style.pack import COLUMN, ROW, LEFT,CENTER
# from gqteaHelp import HelpGqteaWin
# from framesCounter import FramesCounter
# from displayPlots import DisplayPlots


# class BondAngleAnalyser(FramesCounter,DisplayPlots):
#     """Class for analyzing bond angles from a CPMD trajectory file."""

#     async def warning_function(self, msg):
#         await self.main_window.dialog(toga.InfoDialog("Error", f"{msg}"))

#     async def read_params(self, widget):
#         #self.atufs = 0.02418884326505  # Atomic time unit in femtoseconds.

#         async def read_input(text_input, field_name, expected_type):
#             value = text_input.value.strip()
#             if not value:
#                 await self.main_window.dialog(
#                     toga.InfoDialog("Error", f"Please input a valid value for {field_name}.")
#                 )
#                 return None
#             try:
#                 if expected_type == list:
#                     # Assuming atom labels are space-separated integers
#                     labels = [int(label) for label in value.split()]
#                     if len(labels) != 3:
#                         raise ValueError("Please input exactly two atom labels.")
#                     return labels
#                 else:
#                     return expected_type(value)
#             except ValueError as e:
#                 await self.main_window.dialog(
#                     toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
#                 )
#                 return None

#         # Read and validate inputs
#         self.max_angle = await read_input(self.textInput_max_angle, "max_angle", float)
#         if self.max_angle is None:
#             return

#         self.dt = await read_input(self.textInput_time_step, "time step", float)
#         if self.dt is None:
#             return

#         self.sampling = await read_input(self.textInput_sampling_interval, "sampling", float)
#         if self.sampling is None:
#             return

#         self.sim_temp = await read_input(self.textInput_temperature, "simulation temperature", float)
#         if self.sim_temp is None:
#             return

#         self.atom_labels = await read_input(self.textInput_atom_labels, "atom labels", list)
#         if self.atom_labels is None:
#             return

#         self.bin_width = await read_input(self.textInput_bin_width, "bin width", float)
#         if self.bin_width is None:
#             return

#         update_text = (
#             f"Maximum angle          --> {self.max_angle}\n"
#             f"Time step              --> {self.dt}\n"
#             f"Interval sampling      --> {self.sampling}\n"
#             f"Simulation temperature --> {self.sim_temp}\n"
#             f"Selected atom labels   --> {self.atom_labels}\n"
#             f"Histogram bin width    --> {self.bin_width}\n"
#         )
#         self.multi_line_text.value = update_text

#     async def bond_angle(self):
#         """Calculates the bond angle for each frame."""
#         idx1,idx2,idx3 = int(self.atom_labels[0]), int(self.atom_labels[1]), int(self.atom_labels[2])

#         time_step = self.dt
#         sampling_interval = self.sampling
#         atufs = 0.02418884326505
#         time_between_frames = time_step * sampling_interval * atufs / 1000.0

#         self.angles = []
#         line_number = 0
#         frame_number = 1.0

#         try:
#             with open(self.trajec, "r") as f:
#                 self.progress_bar.start()

#                 while True:
#                     title_line = f.readline()
#                     if not title_line:
#                         break
#                     line_number += 1
#                     comment_line = f.readline()
#                     if not comment_line:
#                         break
#                     line_number += 1
#                     atom_data = []

#                     for _ in range(self.num_atoms):
#                         line = f.readline()
#                         if not line:
#                             await self.warning_function(f"Unexpected end of file!")
#                             break
#                         tokens = line.strip().split()
#                         if len(tokens) <4:
#                             await self.warning_function(f"Line {line_number} has invalid atom line format!")
#                         element = tokens[0]
#                         x, y, z = map(float, tokens[1:4])
#                         atom_data.append([element, x, y, z])
                        
#                     #Get positons of the selected atoms

#                     try:
#                         a = np.array(atom_data[idx1 - 1][1:4])
#                         b = np.array(atom_data[idx2 - 1][1:4])
#                         c = np.array(atom_data[idx3 - 1][1:4])

#                         elmt1, elmt2, elmt3 = atom_data[idx1 - 1][0], atom_data[idx2 - 1][0], atom_data[idx3 - 1][0]
                        
                        
#                         vec1 = a - b
#                         vec2 = c - b
#                         vec1_normalized = vec1 / np.linalg.norm(vec1)
#                         vec2_normalized = vec2 / np.linalg.norm(vec2)
                            
#                     except:
#                         self.warning_function(
#                             f"Atom index out of range: line {line_number}"
#                         )

#                     cos_angle = np.dot(vec1_normalized, vec2_normalized)
#                     angle_rad = np.arccos(cos_angle)
#                     angle_deg = np.degrees(angle_rad)

#                     sim_time = frame_number * time_between_frames
#                     self.angles.append((sim_time, angle_deg))

#                     frame_number += 1.0

#                     if (frame_number % 400) == 0:
#                         progress_bar_increment = (frame_number / self.total_frame_number)*100
#                         self.progress_bar.value = progress_bar_increment

#                     if frame_number == self.total_frame_number:
#                         self.progress_bar.value = 100

#                 self.progress_bar.stop()

#             self.elmt1, self.elmt2,self.elmt3 = elmt1, elmt2, elmt3 #To be used later

#             self.save_angles(self.elmt1, idx1, self.elmt2, idx2, elmt3, idx3)

#             # Plotting bond angles
#             y = [sublist[1] for sublist in self.angles]
#             x = [sublist[0] for sublist in self.angles]

#             #Save the temporary plot file
#             plot_xlabel = "Simulation time (ps)"
#             plot_ylabel = "Bond angle (°)"
#             plot_title = "Bond Angles over Time"
#             plot_number = 1
#             self.save_plots(plot_number, x ,y, plot_xlabel, plot_ylabel, plot_title)

#         except:
#             await self.warning_function(self, f"Error calculating bond angles")

#     def save_angles(self, elm1, idx1, elm2, idx2, elm3, idx3):
#         """Saves the calculated angles to a file."""
#         angle_file = os.path.join(
#             self.output_dir, f"angles_{elm1}{idx1}_{elm2}{idx2}_{elm3}{idx3}.dat"
#         )
#         self.angle_file = angle_file

#         with open(angle_file, "w") as f:
#             for sim_time, angle in self.angles:
#                 f.write(f"{sim_time:>10.5f}{angle:>10.5f}\n")

#     async def distribution_function(self):
#         """Calculates the distribution function of the bond angles."""
#         # max_angle = self.max_angle
#         # bin_width = self.bin_width

#         self.num_bins = int(self.max_angle / self.bin_width)
#         self.bin_centers = [self.bin_width * (i + 0.5) for i in range(self.num_bins)]
#         histogram = [0 for _ in range(self.num_bins)]

#         angles_only = [sublist[1] for sublist in self.angles]

#         for angle in angles_only:
#             bin_index = int(angle / self.bin_width)
#             if bin_index < self.num_bins:
#                 histogram[bin_index] += 1

#         self.histogram_sum = sum(histogram)

#         if self.histogram_sum > 0:
#             histogram_percentage = [(value / self.histogram_sum) * 100 for value in histogram]

#         else:
#             msg = (f"The histrogram is empty. Try to increse the maximum bond angle for distribution\n"
#                   f"Try something greater than {max(angles_only)}")
#             await self.warning_function(msg)
#             return  

#         self.histogram = histogram_percentage

#         # Compute statistics
#         average_angle = sum(angles_only) / len(angles_only)
#         var = variance(angles_only)
#         std_dev = stdev(angles_only)
#         largest_angle = max(angles_only)
#         smallest_angle = min(angles_only)

#         self.stats = {
#             "average_angle": average_angle,
#             "variance": var,
#             "std_dev": std_dev,
#             "largest_angle": largest_angle,
#             "smallest_angle": smallest_angle,
#             "bin_centers": self.bin_centers,
#             "histogram_percentage": histogram_percentage,
#         }

#         # Save distribution to file
#         angle_dist_file = os.path.join(self.output_dir, "angles_distribution.dat")
#         with open(angle_dist_file, "w") as f:
#             for center, percent in zip(self.bin_centers, histogram_percentage):
#                 f.write(f"{center:>10.5f}{percent:>10.5f}\n")

#         # Save the bond length plot file and plot
#         x = self.bin_centers
#         y = self.histogram
#         plot_xlabel = "Bond Angles (°)"
#         plot_ylabel = "Bond Angle Distribution (%)"
#         plot_title = "Bond Angle Distribution Plot"
#         plot_number = 2
#         self.save_plots(plot_number, x ,y, plot_xlabel, plot_ylabel, plot_title)

#     def free_energy(self):
#         """Calculates the free energy based on the bond angle distribution."""

#         if not hasattr(self, 'histogram_sum') or self.histogram_sum == 0:
#             self.stats = {}
#             return
                
#         R = 0.001987204  # kcal/(mol*K)

#         free_energy = []
#         for xi, yi in zip(self.bin_centers,self.histogram):
#             if yi > 0.0:
#                 FE = -R * self.sim_temp * np.log(yi / 100.0)
#                 free_energy.append([xi, FE])

#         y = [sublist[1] for sublist in free_energy]
#         x = [sublist[0] for sublist in free_energy]
#         self.min_y = min(y)
#         self.min_y_idx = x[y.index(self.min_y)]

#         # Save free energy to file
#         fe_file = os.path.join(self.output_dir, "angle_free_energy.dat")
#         with open(fe_file, "w") as f:
#             for xi, FE in free_energy:
#                 f.write(f"{xi:>10.5f}{FE:>10.5f}\n")

#         #Save the temporary plot file
#         plot_xlabel = "Bond angles (°)"
#         plot_ylabel = "Bond Angle Free Energy (kcal/mol)"
#         plot_title = "Bond Angle Free Energy plot"
#         plot_number = 3
#         self.save_plots(plot_number, x ,y, plot_xlabel, plot_ylabel, plot_title)

#     def save_summary(self):
#         """Saves the summary of the analysis."""

#         if not self.stats:
#             return

#         idx1, idx2, idx3 = self.atom_labels
#         summary_file = os.path.join(
#             self.output_dir, f"summary_{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}.txt"
#         )
#         self.fs = summary_file

#         with open(summary_file, "w") as f:
#             f.write(f"Number of frames:........... {self.total_frame_number}\n")
#             f.write(f"Number of atoms:............ {self.num_atoms}\n")
#             f.write(f"Selected atoms:............. {self.elmt1}{idx1}-{self.elmt2}{idx2}-{self.elmt3}{idx3}\n")
#             f.write(f"The largest bond angle:..... {self.stats['largest_angle']:.4f} degrees\n")
#             f.write(f"The smallest bond angle:.... {self.stats['smallest_angle']:.4f} degrees\n")
#             f.write(f"Bond angle average:......... {self.stats['average_angle']:.4f} degrees\n")
#             f.write(f"Bond angle variance:........ {self.stats['variance']:.4f}\n")
#             f.write(f"Bond angle std deviation:... {self.stats['std_dev']:.4f}\n")
#             f.write(f"Lowest free energy:......... {self.min_y:.4f} kcal/mol "
#                 f"at {self.min_y_idx:.2f} degrees\n"
#             )

#         with open(summary_file, "r") as f:
#             summary_content = f.read()
#             self.multi_line_text.value = summary_content

# class BondAngleUI(BondAngleAnalyser):

#     def __init__(self, *args):
#         self.layout_main_window(*args)

#     def layout_main_window(self, widget):
#         # Create the main window
#         self.main_window = toga.Window(
#             title="Bond Angle Analysis from TRAJEC.xyz",
#             size=(700, 600),
#         )

#         # Define common styles
#         heading_style = Pack(font_size=16, font_weight="bold", padding=(0, 0, 10, 0))
#         label_style = Pack(padding=(5, 5), text_align=LEFT, width=200)
#         input_style = Pack(flex=1, padding=(5, 5))
#         button_style = Pack(padding=5, width=100)
#         box_style = Pack(direction=ROW, alignment=CENTER, padding=(0, 0, 5, 0))

#         # Main container
#         main_box = toga.Box(style=Pack(direction=COLUMN, padding=20))

#         box_1 = toga.Box(style=Pack(direction=ROW, 
#                         alignment=CENTER, padding=(0,0,10,0))
#         )

#         main_box.add(box_1)
#         box_1a = toga.Box(style=Pack(width=400))
#         box_1b = toga.Box(style=Pack(width=150))
#         box_1c = toga.Box(style=Pack(width=100))
#         box_1.add(box_1a)
#         box_1.add(box_1b)
#         box_1.add(box_1c)

#         # Title
#         title_label = toga.Label("Bond Angle Analysis", 
#                                  style=heading_style
#         )
#         empty_label = toga.Label(" ",style=heading_style
#         )
#         self.progress_label = toga.Label(" ",style=heading_style
#         )

#         box_1a.add(title_label)
#         box_1b.add(empty_label)
#         box_1c.add(self.progress_label)

#         input_fields = [
#             (
#                 "Maximum angle for Distribution Function:",
#                 "Enter 200 for default value",
#                 "textInput_max_angle",
#             ),
#             (
#                 "Simulation Time Step (a.u.):", 
#                 "Enter time step value", 
#                 "textInput_time_step"
#             ),

#             (
#                 "Sampling Interval (frames):",
#                 "Enter sampling interval",
#                 "textInput_sampling_interval",
#             ),
#             (
#                 "Simulation Temperature (K):",
#                 "Enter temperature",
#                 "textInput_temperature",
#             ),
#             (
#                 "Atom Labels (e.g., 2 3 4):",
#                 "Enter three atom labels, separated by a space (e.g., 1 2 3)",
#                 "textInput_atom_labels",
#             ),
#             (
#                 "Histogram Bin Width (°):",
#                 "Enter 1 for default value",
#                 "textInput_bin_width",
#             ),
#         ]

#         for label_text, placeholder, attr_name in input_fields:
#             box = toga.Box(style=box_style)
#             label = toga.Label(label_text, style=label_style)
#             text_input = toga.TextInput(placeholder=placeholder, style=input_style)
#             setattr(self, attr_name, text_input)
#             box.add(label)
#             box.add(text_input)
#             main_box.add(box)

#         # File selection button
#         file_box = toga.Box(style=box_style)
#         file_label = toga.Label("Select Trajectory File:", style=label_style)
#         self.textInput_file = toga.TextInput(
#             placeholder="Click Browse to select TRAJEC.xyz file", style=input_style
#         )
#         browse_button = toga.Button(
#             "Browse", on_press = self.frames_counter, style=button_style
#         )

#         progress_box = toga.Box(style=Pack(direction=COLUMN))
#         self.progress_bar = toga.ProgressBar(max=100)
#         progress_box.add(self.progress_bar)

#         file_box.add(file_label)
#         file_box.add(self.textInput_file)
#         file_box.add(browse_button)
#         main_box.add(file_box)
#         main_box.add(progress_box)

#         # Multi-line text input for help or output
#         self.multi_line_text = toga.MultilineTextInput(
#             style=Pack(flex=1, padding=(10, 0), font_size=12)
#         )
#         self.multi_line_text.value = HelpGqteaWin.help_bond_angle
#         main_box.add(self.multi_line_text)

#         # Buttons at the bottom
#         button_box = toga.Box(
#             style=Pack(direction=ROW, alignment=CENTER, padding_top=10)
#         )
#         self.btn_execute = toga.Button(
#             "Analyze", style=button_style, on_press=self.workflow
#         )
#         self.btn_close = toga.Button(
#             "Close", style=button_style, on_press=self.close_bond_angle_window
#         )
#         button_box.add(self.btn_execute)
#         button_box.add(self.btn_close)
#         main_box.add(button_box)

#         # Set the content of the main window
#         self.main_window.content = main_box
#         self.main_window.show()

#     async def workflow(self, widget):
#         await self.read_params(widget)
#         if not hasattr(self, "trajec") or not self.trajec:
#             await self.warning_function("No trajectory file selected.")
#             return
#         await self.bond_angle()
#         await self.distribution_function()
#         self.free_energy()
#         self.save_summary()
#         self.display_plots()


#     def close_bond_angle_window(self, widget):
#         self.main_window.close()

