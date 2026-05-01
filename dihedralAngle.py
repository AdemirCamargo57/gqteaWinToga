# Refactoring dihedralAngle.py to improve structure, error handling, and code clarity at 2025-11-08
import os
import numpy as np
import toga
from statistics import variance, stdev
from toga.style import Pack
from toga.constants import LEFT  # for label text alignment

from help import HelpGqteaWin
from framesCounter import FramesCounter
from displayPlots import DisplayPlots


class DihedralAngleAnalyser(FramesCounter, DisplayPlots):
    """Analyze dihedral (torsion) angles i-j-k-l from XYZ/CPMD-like trajectories."""

    atufs: float = 0.02418884326505  # Atomic time unit in fs

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
        """Ensure output_dir exists; default to TRAJEC folder or CWD."""
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

    async def read_params(self, widget) -> None:
        """Read and validate inputs, echo to multiline box."""

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
                    if len(labels) != 4:
                        raise ValueError("Please input exactly four atom labels.")
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
            return

        self.dt = await read_input(self.textInput_time_step, "time step (a.u.)", float)
        if self.dt is None:
            return

        self.sampling = await read_input(self.textInput_sampling_interval, "sampling (frames)", int)
        if self.sampling is None:
            return

        self.sim_temp = await read_input(self.textInput_temperature, "simulation temperature (K)", float)
        if self.sim_temp is None:
            return

        self.atom_labels = await read_input(self.textInput_atom_labels, "atom labels (i j k l)", list)
        if self.atom_labels is None:
            return

        self.bin_width = await read_input(self.textInput_bin_width, "bin width (deg)", float)
        if self.bin_width is None:
            return

        # UI feedback
        update_text = (
            f"Maximum dihedral (deg) --> {self.max_angle}\n"
            f"Time step (a.u.)       --> {self.dt}\n"
            f"Sampling (frames)      --> {self.sampling}\n"
            f"Temperature (K)        --> {self.sim_temp}\n"
            f"Atom labels (i j k l)  --> {self.atom_labels}\n"
            f"Bin width (deg)        --> {self.bin_width}\n"
            f"Show plots at the end  --> {self.switch_show_plots.value}\n"
            f"Save CSV outputs       --> {self.switch_save_csv.value}\n"
            f"Wrap to [-180, 180]    --> {self.switch_wrap_180.value}\n"
        )
        self.multi_line_text.value = update_text

    @staticmethod
    def _dihedral_deg(a, b, c, d) -> float:
        """Compute dihedral angle (degrees) using stable atan2 formulation.
        Returns angle in (-180, 180]."""
        ab = b - a
        bc = c - b
        cd = d - c

        # Norms and guards
        nbc = np.linalg.norm(bc)
        if nbc == 0.0:
            return None

        n1 = np.cross(ab, bc)
        n2 = np.cross(bc, cd)

        n1_norm = np.linalg.norm(n1)
        n2_norm = np.linalg.norm(n2)
        if n1_norm == 0.0 or n2_norm == 0.0:
            return None

        n1u = n1 / n1_norm
        n2u = n2 / n2_norm
        bcu = bc / nbc

        x = np.dot(n1u, n2u)
        y = np.dot(np.cross(n1u, bcu), n2u)
        angle_rad = np.arctan2(y, x)
        return float(np.degrees(angle_rad))  # (-180, 180]

    async def dihedral_angle(self, widget) -> None:
        """Compute dihedral angles for each frame and save the time series."""
        self._ensure_output_dir()
        self.dihedral = []

        idx1, idx2, idx3, idx4 = map(int, self.atom_labels)
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
                        tok = line.strip().split()
                        if len(tok) < 4:
                            await self.warning_function(f"Line {line_number} has invalid atom line format.")
                            self.progress_bar.stop()
                            return
                        element = tok[0]
                        x, y, z = map(float, tok[1:4])
                        atom_data.append((element, x, y, z))
                        line_number += 1

                    try:
                        a = np.array(atom_data[idx1 - 1][1:4], dtype=float)
                        b = np.array(atom_data[idx2 - 1][1:4], dtype=float)
                        c = np.array(atom_data[idx3 - 1][1:4], dtype=float)
                        d = np.array(atom_data[idx4 - 1][1:4], dtype=float)

                        elmt1 = atom_data[idx1 - 1][0]
                        elmt2 = atom_data[idx2 - 1][0]
                        elmt3 = atom_data[idx3 - 1][0]
                        elmt4 = atom_data[idx4 - 1][0]
                    except IndexError:
                        await self.warning_function(f"Atom index out of range at line {line_number}.")
                        self.progress_bar.stop()
                        return

                    angle = self._dihedral_deg(a, b, c, d)
                    if angle is None:
                        # Degenerate geometry; skip this frame
                        frame_number += 1
                        self._update_progress(frame_number, self.total_frame_number)
                        continue

                    # Optionally wrap into [-180, 180]; otherwise map to [0, 360)
                    if self.switch_wrap_180.value:
                        angle_deg = angle  # already (-180, 180]
                    else:
                        angle_deg = angle if angle >= 0.0 else angle + 360.0  # [0, 360)

                    sim_time = frame_number * time_increment_ps
                    self.dihedral.append([sim_time, float(angle_deg)])

                    frame_number += 1
                    self._update_progress(frame_number, self.total_frame_number)

                self.progress_bar.value = 100
                self.progress_bar.stop()

        except Exception as e:
            await self.warning_function(f"Error calculating dihedral angles: {e}")
            return

        # Save tags for filenames
        self.elmt1, self.elmt2, self.elmt3, self.elmt4 = elmt1, elmt2, elmt3, elmt4

        # Time series file
        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}_{self.elmt4}{idx4}"
        self.ftimes = os.path.join(self.output_dir, f"dihedral_{tag}.dat")
        with open(self.ftimes, "w") as f:
            for t, ang in self.dihedral:
                f.write(f"{t:>12.6f}{ang:>12.6f}\n")

        # Plot time series
        xs = [row[0] for row in self.dihedral]
        ys = [row[1] for row in self.dihedral]
        self.save_plots(1, xs, ys, "Simulation time (ps)", "Dihedral angle (°)", "Dihedral Angle Over Time")

    def distribution_function(self) -> None:
        """Histogram (%) and stats."""
        if not self.dihedral:
            return

        angles_only = [ang for _, ang in self.dihedral]

        self.num_bins = max(1, int(self.max_angle / self.bin_width))
        self.bin_centers = [self.bin_width * (i + 0.5) for i in range(self.num_bins)]
        counts = [0] * self.num_bins

        for ang in angles_only:
            idx = int(ang / self.bin_width)
            if idx == self.num_bins:  # when ang == max_angle
                idx -= 1
            if 0 <= idx < self.num_bins:
                counts[idx] += 1

        total = sum(counts)
        if total == 0:
            # Leave quietly; free_energy() will no-op.
            self.histogram = [0.0] * self.num_bins
            return

        histogram_pct = [(c / total) * 100.0 for c in counts]
        self.histogram = histogram_pct

        # Stats (safe for N=1)
        avg = sum(angles_only) / len(angles_only)
        largest = max(angles_only)
        smallest = min(angles_only)
        var_val = variance(angles_only) if len(angles_only) > 1 else 0.0
        std_val = stdev(angles_only) if len(angles_only) > 1 else 0.0

        # Write distribution to file
        idx1, idx2, idx3, idx4 = map(int, self.atom_labels)
        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}_{self.elmt4}{idx4}"
        self.fdist = os.path.join(self.output_dir, f"dihedral_distribution_{tag}.dat")
        with open(self.fdist, "w") as f:
            for center, pct in zip(self.bin_centers, histogram_pct):
                f.write(f"{center:>12.6f}{pct:>12.6f}\n")

        # Persist stats for summary
        self.stats = {
            "average_dihedral": avg,
            "variance": var_val,
            "std_dev": std_val,
            "largest_dihedral": largest,
            "smallest_dihedral": smallest,
        }

        # Plot histogram
        self.save_plots(
            2, self.bin_centers, self.histogram,
            "Dihedral angle (°)", "Dihedral distribution (%)",
            "Dihedral Angle Distribution"
        )

    def dihedral_free_energy(self) -> None:
        """Compute free energy (kcal/mol) from P(φ)."""
        if not hasattr(self, "histogram") or not self.histogram:
            self.stats = {}
            return

        R = 0.001987204  # kcal/(mol*K)
        T = self.sim_temp

        raw_probabilities = []
        for phi_deg, pct in zip(self.bin_centers, self.histogram):
            if pct <= 0.0:
                continue
            p = pct / 100.0
            raw_probabilities.append([phi_deg, max(p, 1e-12)])

        if not raw_probabilities:
            return

        xs = [phi for phi, _ in raw_probabilities]
        smoothed_probabilities = self._smooth_series([p for _, p in raw_probabilities])

        self.free_energy_pairs = []
        for phi_deg, p_smooth in zip(xs, smoothed_probabilities):
            FE = -R * T * np.log(max(p_smooth, 1e-12))
            self.free_energy_pairs.append([phi_deg, FE])

        ys = [y for _, y in self.free_energy_pairs]

        self.min_y = float(min(ys))
        self.min_y_idx = float(xs[ys.index(self.min_y)])

        # Save FE
        idx1, idx2, idx3, idx4 = map(int, self.atom_labels)
        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}_{self.elmt4}{idx4}"
        self.ffe = os.path.join(self.output_dir, f"dihedral_free_energy_{tag}.dat")
        with open(self.ffe, "w") as f:
            for x, fe in self.free_energy_pairs:
                f.write(f"{x:>12.6f}{fe:>12.6f}\n")

        # Append summary + show in textbox
        self.fsummary = os.path.join(self.output_dir, f"summary_{tag}.txt")
        with open(self.fsummary, "w") as fs:
            fs.write(f"Number of frames:........... {self.total_frame_number}\n")
            fs.write(f"Number of atoms:............ {self.num_atoms}\n")
            fs.write(
                "Selected atoms:............ "
                f"{self.elmt1}{idx1}-{self.elmt2}{idx2}-{self.elmt3}{idx3}-{self.elmt4}{idx4}\n"
            )
            if hasattr(self, "stats") and self.stats:
                fs.write(f"The largest dihedral angle:..... {self.stats['largest_dihedral']:.4f} deg\n")
                fs.write(f"The smallest dihedral angle:.... {self.stats['smallest_dihedral']:.4f} deg\n")
                fs.write(f"Dihedral angle average:......... {self.stats['average_dihedral']:.4f} deg\n")
                fs.write(f"Dihedral angle variance:........ {self.stats['variance']:.4f}\n")
                fs.write(f"Dihedral angle std deviation:... {self.stats['std_dev']:.4f}\n")
            fs.write(
                f"Lowest free energy:......... {self.min_y:.4f} kcal/mol "
                f"at {self.min_y_idx:.2f} deg\n"
            )

        with open(self.fsummary, "r") as fs:
            self.multi_line_text.value = fs.read()

        # Plot FE
        self.save_plots(
            3, xs, ys,
            "Dihedral angle (°)", "Free energy (kcal/mol)",
            "Dihedral Free Energy"
        )

    def export_csv(self) -> None:
        """CSV export with headers for timeseries, distribution, FE."""
        idx1, idx2, idx3, idx4 = map(int, self.atom_labels)
        tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}_{self.elmt3}{idx3}_{self.elmt4}{idx4}"

        # Time series
        csv_ts = os.path.join(self.output_dir, f"dihedral_{tag}.csv")
        with open(csv_ts, "w") as f:
            f.write("time_ps,dihedral_deg\n")
            for t, ang in self.dihedral:
                f.write(f"{t:.6f},{ang:.6f}\n")

        # Distribution
        if hasattr(self, "bin_centers") and hasattr(self, "histogram"):
            csv_dist = os.path.join(self.output_dir, f"dihedral_distribution_{tag}.csv")
            with open(csv_dist, "w") as f:
                f.write("dihedral_deg,distribution_pct\n")
                for x, pct in zip(self.bin_centers, self.histogram):
                    f.write(f"{x:.6f},{pct:.6f}\n")

        # Free energy
        if hasattr(self, "free_energy_pairs") and self.free_energy_pairs:
            csv_fe = os.path.join(self.output_dir, f"dihedral_free_energy_{tag}.csv")
            with open(csv_fe, "w") as f:
                f.write("dihedral_deg,free_energy_kcal_mol\n")
                for x, fe in self.free_energy_pairs:
                    f.write(f"{x:.6f},{fe:.6f}\n")


class DihedralUI(DihedralAngleAnalyser):
    """UI for the Dihedral Angle module (harmonized with Bond/BondAngle)."""

    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget) -> None:
        # Window
        self.main_window = toga.Window(
            title="Dihedral Angle Analysis from TRAJEC.xyz",
            size=(760, 700),
        )

        # Styles (use margin + direction strings)
        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style   = Pack(margin=(5, 5), text_align=LEFT, width=270)
        input_style   = Pack(flex=1, margin=(5, 5))
        button_style  = Pack(margin=5, width=110)
        row_style     = Pack(direction="row", margin=(0, 0, 5, 0))
        col_style     = Pack(direction="column", margin=(0, 0, 5, 0))

        # Main container
        main_box = toga.Box(style=Pack(direction="column", margin=20))

        # Header row
        box_1 = toga.Box(style=Pack(direction="row", margin=(0, 0, 10, 0)))
        main_box.add(box_1)

        box_1a = toga.Box(style=Pack(width=480))
        box_1b = toga.Box(style=Pack(width=180))
        box_1c = toga.Box(style=Pack(width=100))
        box_1.add(box_1a)
        box_1.add(box_1b)
        box_1.add(box_1c)

        title_label = toga.Label("Dihedral Angle Analysis", style=heading_style)
        self.progress_label = toga.Label(" ", style=heading_style)
        self.status_label = toga.Label(" ", style=Pack(margin=(0, 0, 0, 10)))  # will print output dir

        box_1a.add(title_label)
        box_1b.add(self.progress_label)
        box_1c.add(self.status_label)

        # Inputs
        input_fields = [
            ("Maximum angle for Distribution (deg):", "Enter e.g. 360", "textInput_max_angle"),
            ("Simulation Time Step (a.u.):", "Enter time step", "textInput_time_step"),
            ("Sampling Interval (frames):", "Enter sampling interval", "textInput_sampling_interval"),
            ("Simulation Temperature (K):", "Enter temperature", "textInput_temperature"),
            ("Atom Labels (i j k l):", "Enter four labels, e.g., 1 2 3 4", "textInput_atom_labels"),
            ("Histogram Bin Width (deg):", "Enter e.g. 5.0", "textInput_bin_width"),
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

        # Switches
        switches_row = toga.Box(style=row_style)
        self.switch_show_plots = toga.Switch("Show plots at the end", value=True,  style=Pack(margin=(4, 10, 4, 0)))
        self.switch_save_csv   = toga.Switch("Save CSV outputs",     value=True,  style=Pack(margin=(4, 10, 4, 10)))
        self.switch_wrap_180   = toga.Switch("Wrap dihedral to [-180, 180]", value=False, style=Pack(margin=(4, 10, 4, 10)))
        switches_row.add(self.switch_show_plots)
        switches_row.add(self.switch_save_csv)
        switches_row.add(self.switch_wrap_180)

        # Progress bar
        progress_col = toga.Box(style=col_style)
        self.progress_bar = toga.ProgressBar(max=100)
        progress_col.add(self.progress_bar)

        main_box.add(switches_row)
        main_box.add(progress_col)

        # Help / Output
        self.multi_line_text = toga.MultilineTextInput(style=Pack(flex=1, margin=(10, 0), font_size=12))
        self.multi_line_text.value = HelpGqteaWin.help_dihedral_angle
        main_box.add(self.multi_line_text)

        # Buttons
        button_row = toga.Box(style=Pack(direction="row", margin=(10, 0, 0, 0)))
        self.btn_execute = toga.Button("Analyze", style=button_style, on_press=self.workflow)
        self.btn_close   = toga.Button("Close",   style=button_style, on_press=self.closeTopLevel)
        button_row.add(self.btn_execute)
        button_row.add(self.btn_close)
        main_box.add(button_row)

        self.main_window.content = main_box
        self.main_window.show()

    async def workflow(self, widget) -> None:
        await self.read_params(widget)
        if not hasattr(self, "trajec") or not self.trajec:
            await self.warning_function("No trajectory file selected.")
            return

        await self.dihedral_angle(widget)
        self.distribution_function()
        self.dihedral_free_energy()

        if self.switch_save_csv.value:
            self.export_csv()
        if self.switch_show_plots.value:
            self.display_plots()

        # Print output directory
        out_dir = getattr(self, "output_dir", os.getcwd())
        msg = f"Outputs saved to: {out_dir}"
        self.status_label.text = msg
        self.multi_line_text.value = (self.multi_line_text.value.rstrip() + "\n\n" + msg + "\n")

    def closeTopLevel(self, widget) -> None:
        self.main_window.close()



# old code for reference
# import os, toga 
# import numpy as np
# from statistics import stdev, variance
# from toga.style import Pack
# from toga.style.pack import COLUMN, ROW, LEFT,CENTER
# from gqteaHelp import HelpGqteaWin
# import matplotlib.pyplot as plt
# from framesCounter import FramesCounter
# from displayPlots import DisplayPlots


# class DihedralAngleAnalyser(FramesCounter, DisplayPlots):
    
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
#                     if len(labels) != 4:
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
#             f"Maximum Dihedral Angle --> {self.max_angle}\n"
#             f"Time Step              --> {self.dt}\n"
#             f"Interval Sampling      --> {self.sampling}\n"
#             f"Simulation Temperature --> {self.sim_temp}\n"
#             f"Selected Atom Labels   --> {self.atom_labels}\n"
#             f"Histogram Bin Width    --> {self.bin_width}\n"
#         )
#         self.multi_line_text.value = update_text
   
#     async def dihedral_angle(self,widget):

#         atufs=0.02418884326505
#         self.dihedral = []

#         # Time between collected frames     
#         time_btf = self.dt * self.sampling * atufs / 1000.0
#         num_line = 0
#         t = 1.0
#         frame_number = 0
#         idx1 = self.atom_labels[0]
#         idx2 = self.atom_labels[1]
#         idx3 = self.atom_labels[2]
#         idx4 = self.atom_labels[3]

#         with open(self.trajec) as f:
#             self.progress_bar.start()
#             while True:
#                 title_line = f.readline()
#                 if not title_line:
#                     break
#                 comment_line = f.readline()
#                 if not comment_line:
#                     break                

#                 num_line += 2
#                 try:
#                    for n in range(self.num_atoms):
#                        line = f.readline().split()
#                        num_line += 1
#                        line[1] = float(line[1])       # Positions  x in angstroms
#                        line[2] = float(line[2])       # Positions  y in angstroms
#                        line[3] = float(line[3])       # Positions  z in angstroms
#                        if idx1 == n + 1:
#                            elm1 = line[0]
#                            coords1 = [line[1],line[2],line[3]]
#                        elif idx2 == n + 1:
#                            elm2 = line[0]
#                            coords2 = [line[1],line[2],line[3]]
#                        elif idx3 == n + 1:
#                            elm3 = line[0]
#                            coords3 = [line[1],line[2],line[3]]
#                        elif idx4 == n + 1:
#                            elm4 = line[0]
#                            coords4 = [line[1],line[2],line[3]]
                        
#                 except:
#                     msg = f"Error while reading line {num_line} in TRAJEC.xyz"
#                     await self.warning_function(msg)
#                     break

#                 a = np.array(coords1)
#                 b = np.array(coords2)
#                 c = np.array(coords3)
#                 d = np.array(coords4)

#                 ab = b - a
#                 bc = c - b
#                 cd = d - c

#                 normal1 = np.cross(ab,bc)
#                 normal2 = np.cross(bc,cd)

#                 x = np.dot(normal1, normal2)
#                 y = np.dot(np.cross(normal1, bc/np.linalg.norm(bc)), normal2)
#                 angle_rad = np.arctan2(y, x)
#                 angle_deg = np.degrees(angle_rad)
        
#                 sim_time = t*time_btf

#                 if angle_deg >= 0:
#                     self.dihedral.append([sim_time, angle_deg])
#                 else:
#                     self.dihedral.append([sim_time, 360.0 + angle_deg])

#                 t += 1.0
#                 frame_number += 1

#                 if (frame_number % 400) == 0:
#                     progress_bar_increment = (frame_number / self.total_frame_number)*100
#                     self.progress_bar.value = progress_bar_increment

#                 if frame_number == self.total_frame_number:
#                     self.progress_bar.value = 100

#         self.progress_bar.stop()
            

#         f0 = f"{self.output_dir}/dihedral_{elm1}{idx1}_{elm2}{idx2}_{elm3}{idx3}_{elm4}{idx4}.dat"  
#         self.f1 = f"{self.output_dir}/dihedral_free_energy_{elm1}{idx1}_{elm2}{idx2}_{elm3}{idx3}_{elm4}{idx4}.dat"
#         self.f2 = f"{self.output_dir}/dihedral_dist_{elm1}{idx1}_{elm2}{idx2}_{elm3}{idx3}_{elm4}{idx4}.dat"
#         self.fs = f"{self.output_dir}/summary_{elm1}{idx1}-{elm2}{idx2}-{elm3}{idx3}-{elm4}{idx4}.txt"

       
#         with open(self.fs, "w") as fs:
#             fs.write(f"Number of frames:...................{self.total_frame_number}\n")
#             fs.write(f"File number of lines:...............{num_line}\n")
#             fs.write(f"Number of atoms:....................{self.num_atoms}\n")            

#         with open(f0,"w") as f:
#             for k in range(len(self.dihedral)):
#                 p0 = self.dihedral[k][0]
#                 p1 = self.dihedral[k][1]
#                 f.write(f"{p0:>15.5f}{p1:>15.5f}\n")

#         with open(self.fs,"a") as f:
#             f.write(f"Selected atoms:.....................{elm1}{idx1}-{elm2}{idx2}-{elm3}{idx3}-{elm4}{idx4}\n") 
        

#         #Plotting bond angles
#         x = [sublist[0] for sublist in self.dihedral]
#         y = [sublist[1] for sublist in self.dihedral]

#         # Save the bond length plot file and plot
#         plot_xlabel = "Simulation Time (ps)"
#         plot_ylabel = "Dihedral Angle (°)"
#         plot_title = "Dihedral Angle Over Time"
#         plot_number = 1
#         self.save_plots(plot_number, x ,y, plot_xlabel, plot_ylabel, plot_title)        
    
#     def distribution_function(self):
#         # bin_width = self.bin_width

#         num_bins = int(self.max_angle / self.bin_width)
#         self.bin_centers = [self.bin_width * (i + 0.5) for i in range(num_bins)]
#         histogram = [0 for _ in range(num_bins)]

#         angles_only = [sublist[1] for sublist in self.dihedral]

#         for angle in angles_only:
#             bin_index = int(angle / self.bin_width)
#             if bin_index < num_bins:
#                 histogram[bin_index] += 1

#         total_counts = sum(histogram)
#         histogram_percentage = [(count / total_counts) * 100 for count in histogram]

#         self.histogram = histogram_percentage
#         self.num_bins = num_bins
    
#         # Compute statistics
#         average_dihedral = sum(angles_only) / len(angles_only)
#         largest_dihedral = max(angles_only)
#         smallest_dihedral = min(angles_only)
#         var = variance(angles_only)
#         std = stdev(angles_only)
       
#         with open(self.fs, "a") as f:
#             f.write(f"The largest dihedral angle:.........{largest_dihedral:.4f} degrees\n")
#             f.write(f"The smallest dihedral angle:........{smallest_dihedral:.4f} degrees\n")
#             f.write(f"Dihedral angle average:.............{average_dihedral:.4f} degrees\n")
#             f.write(f"Dihedral angle variance:............{var:.4f}\n")
#             f.write(f"Dihedral angle stdev:...............{std:.4f}\n")

#         x = self.bin_centers
#         y = self.histogram

#         with open(self.f2, "w") as f:
#             for n in range(len(y)):
#                 f.write(f"{x[n]:>14.5f}{y[n]:>14.5f}\n")

#         # Save the bond length plot file and plot
#         plot_xlabel = "Dihedral angle (°)"
#         plot_ylabel = "Dihedral probability distribution"
#         plot_title = "Dihhedral Angle Distribution Function"
#         plot_number = 2
#         self.save_plots(plot_number, x ,y, plot_xlabel, plot_ylabel, plot_title)        
    
#     def dihedral_free_energy(self):
#         """Calculate the free energy based on the dihedral angle distribuition."""

#         free_energy = []
#         R = 0.001987204
#         T = self.sim_temp

#         for xi, yi in zip(self.bin_centers, self.histogram):
#             if yi > 0.0:
#                 FE = -R* self.sim_temp * np.log(yi / 100.0)
#                 free_energy.append([xi, FE])

#         with open(self.f1,"w") as f1:
#             for n in range(len(free_energy)):
#                 p0,p1 = free_energy[n][0],free_energy[n][1]
#                 f1.write(f"{p0:>14.7f}{p1:>14.7f}\n")
        
#         y = [sublist[1] for sublist in free_energy]
#         x = [sublist[0] for sublist in free_energy]

#         min_y = min(y)
#         min_y_idx = y.index(min_y)

#         with open(self.fs, "a") as fs:
#             tmp = f"Lowest free energy..................{min_y:.4f} kcal/mol at {x[min_y_idx]:.2f} degrees\n"
#             fs.write(tmp)

#         # Save the Dihedral angle plot file
#         plot_xlabel = "Dihedral Angle (°)"
#         plot_ylabel = "Free Energy (kcal/mol)"
#         plot_title = "Free Energy Distribution Function"
#         plot_number = 3
#         self.save_plots(plot_number, x , y, plot_xlabel, plot_ylabel, plot_title)
        
#         # Iputing the the calculation sammay into multi_line_text
#         with open(self.fs,"r") as fs:
#             summary_content = fs.read()
#             self.multi_line_text.value = summary_content


# class DihedralUI(DihedralAngleAnalyser):
#     def __init__(self, *args):
#         self.layout_main_window(*args)

#     def layout_main_window(self, widget):
#         # Create the main window
#         self.main_window = toga.Window(
#             title="DIHEDRAL ANGLE ANALYSIS",
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
#         title_label = toga.Label("Dihedral Angle Analysis", 
#                                  style=heading_style
#         )
#         empty_label = toga.Label(" ",style=heading_style
#         )
#         self.progress_label = toga.Label(" ",style=heading_style
#         )

#         box_1a.add(title_label)
#         box_1b.add(empty_label)
#         box_1c.add(self.progress_label)

#         # Input fields with labels
#         input_fields = [
#             (
#                 "Maximum Angle for Distribution Function:",
#                 "Enter 400 for default value",
#                 "textInput_max_angle",
#             ),
#             ("Simulation Time Step (a.u.):", "Enter simulation time step", "textInput_time_step"),
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
#                 "Atom Labels (e.g., 2 3 4 5):",
#                 "Enter four atom labels, separated by a space (e.g., 1 2 3 4)",
#                 "textInput_atom_labels",
#             ),
#             (
#                 "Histogram Bin Width (º):",
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

#         progress_box = toga.Box(style=Pack(direction=COLUMN, padding= (0,0,0,0)))
#         self.progress_bar = toga.ProgressBar(max=100)
#         progress_box.add(self.progress_bar)

#         file_box.add(file_label)
#         file_box.add(self.textInput_file)
#         file_box.add(browse_button)
#         main_box.add(file_box)
#         main_box.add(progress_box)

#         # Multi-line text input for help or output
#         self.multi_line_text = toga.MultilineTextInput(
#             style=Pack(flex=1, padding=(5, 0, 0, 0), font_size=12)
#         )
#         self.multi_line_text.value = HelpGqteaWin.help_dihedral_angle
#         main_box.add(self.multi_line_text)

#         # Buttons at the bottom
#         button_box = toga.Box(
#             style=Pack(direction=ROW, alignment=CENTER, padding_top=5)
#         )
#         self.btn_execute = toga.Button(
#             "Analyze", style=button_style, on_press=self.workflow
#         )
#         self.btn_close = toga.Button(
#             "Close", style=button_style, on_press=self.closeTopLevel
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
#         await self.dihedral_angle(widget)
#         self.distribution_function()
#         self.dihedral_free_energy()

#         # Display the plots within the application
#         self.display_plots()

#     def closeTopLevel(self, widget):
#         self.main_window.close()
