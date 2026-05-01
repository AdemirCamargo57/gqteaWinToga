# Refactored bond analysis module for GQTEA at 2025 November
import os
import numpy as np
import toga
from statistics import variance, stdev
from toga.style import Pack
from toga.constants import LEFT  # for label text alignment

from help import HelpGqteaWin
from framesCounter import FramesCounter
from displayPlots import DisplayPlots


class BondAnalyser(FramesCounter, DisplayPlots):
    """Numerical analysis for bond length, distribution and free energy."""

    atufs: float = 0.02418884326505  # Atomic time unit in femtoseconds.

    async def warning_function(self, msg: str) -> None:
        await self.main_window.dialog(toga.InfoDialog("Error", f"{msg}"))

    def _ensure_output_dir(self) -> None:
        """Ensure self.output_dir exists; default to the TRAJEC directory if missing."""
        if not getattr(self, "output_dir", None):
            if hasattr(self, "trajec") and self.trajec:
                self.output_dir = os.path.dirname(os.path.abspath(self.trajec)) or os.getcwd()
            else:
                self.output_dir = os.getcwd()
        os.makedirs(self.output_dir, exist_ok=True)

    def _reset_analysis_state(self) -> None:
        """Clear analysis results so failed runs cannot reuse stale data."""
        self.bond_lengths = []
        self.bin_centers = []
        self.histogram = []
        self.histogram_sum = 0
        self.stats = {}
        self.free_energy_pairs = []
        self.average_bond = None
        self.min_y = None
        self.min_y_idx = None
        self.elmt1 = None
        self.elmt2 = None
        self.fb = None
        self.fe = None
        self.fd = None
        self.fs = None

    async def read_params(self, widget) -> bool:
        """Read and validate user inputs from the UI."""

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
                    if len(labels) != 2:
                        raise ValueError("Please input exactly two atom labels.")
                    return labels
                else:
                    return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        self.max_r = await read_input(self.textInput_max_r, "max_r", float)
        if self.max_r is None:
            return False

        self.dt = await read_input(self.textInput_time_step, "time step (a.u.)", float)
        if self.dt is None:
            return False

        # sampling is number of frames → integer
        self.sampling = await read_input(self.textInput_sampling_interval, "sampling (frames)", int)
        if self.sampling is None:
            return False

        self.sim_temp = await read_input(self.textInput_temperature, "simulation temperature (K)", float)
        if self.sim_temp is None:
            return False

        self.atom_labels = await read_input(self.textInput_atom_labels, "atom labels", list)
        if self.atom_labels is None:
            return False

        self.bin_width = await read_input(self.textInput_bin_width, "bin width (Å)", float)
        if self.bin_width is None:
            return False

        if self.max_r <= 0:
            await self.warning_function("Maximum r must be greater than zero.")
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
        if self.bin_width >= self.max_r:
            await self.warning_function("Histogram bin width must be smaller than maximum r.")
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

        # UI feedback (read switches via .value)
        update_text = (
            f"Maximum r --> {self.max_r}\n"
            f"Time step (a.u.) --> {self.dt}\n"
            f"Interval sampling (frames) --> {self.sampling}\n"
            f"Simulation temperature (K) --> {self.sim_temp}\n"
            f"Selected atom labels --> {self.atom_labels}\n"
            f"Histogram bin width (Å) --> {self.bin_width}\n"
            f"Show plots at the end --> {self.switch_show_plots.value}\n"
            f"Save CSV outputs --> {self.switch_save_csv.value}\n"
        )
        self.multi_line_text.value = update_text
        return True

    def _time_increment_ps(self) -> float:
        """Increment of simulation time between stored frames in ps."""
        return self.dt * self.sampling * self.atufs / 1000.0

    def _update_progress(self, current: int, total: int) -> None:
        if total <= 0:
            return
        if current % 10 == 0 or current == total:
            self.progress_bar.value = (current / total) * 100

    async def bond_length(self, widget) -> bool:
        """Compute bond lengths per frame and write time series file."""

        self._ensure_output_dir()
        self.bond_lengths = []
        time_increment_ps = self._time_increment_ps()

        idx1 = int(self.atom_labels[0])
        idx2 = int(self.atom_labels[1])

        try:
            with open(self.trajec, "r") as f:
                self.progress_bar.start()
                frame_number = 0
                while True:
                    title_line = f.readline()
                    if not title_line:
                        break
                    comment_line = f.readline()
                    if not comment_line:
                        break

                    atom_data = []
                    for _ in range(self.num_atoms):
                        line = f.readline()
                        if not line:
                            raise ValueError("Unexpected end of file.")
                        tokens = line.strip().split()
                        if len(tokens) < 4:
                            raise ValueError("Invalid atom line format.")
                        element = tokens[0]
                        x, y, z = map(float, tokens[1:4])
                        atom_data.append((element, x, y, z))

                    try:
                        elmt1, x1, y1, z1 = atom_data[idx1 - 1]
                        elmt2, x2, y2, z2 = atom_data[idx2 - 1]
                    except IndexError:
                        await self.warning_function(
                            f"Atom index out of range in frame {frame_number}."
                        )
                        self.progress_bar.stop()
                        return False

                    dx, dy, dz = (x2 - x1), (y2 - y1), (z2 - z1)
                    bond_len = float(np.sqrt(dx * dx + dy * dy + dz * dz))
                    sim_time = frame_number * time_increment_ps
                    self.bond_lengths.append([sim_time, bond_len])

                    frame_number += 1
                    self._update_progress(frame_number, self.total_frame_number)

                self.progress_bar.value = 100
                self.progress_bar.stop()

        except Exception as e:
            await self.warning_function(f"Error reading TRAJEC.xyz file: {e}")
            return False

        if not self.bond_lengths:
            await self.warning_function("No frames were read from the trajectory file.")
            return False

        # Save element names
        self.elmt1, self.elmt2 = elmt1, elmt2

        # File basenames
        self.fb = os.path.join(self.output_dir, f"bond_{elmt1}{idx1}_{elmt2}{idx2}.dat")
        self.fe = os.path.join(self.output_dir, f"free_energy_{elmt1}{idx1}_{elmt2}{idx2}.dat")
        self.fd = os.path.join(self.output_dir, f"bond_distribution_{elmt1}{idx1}_{elmt2}{idx2}.dat")
        self.fs = os.path.join(self.output_dir, f"summary_{elmt1}{idx1}_{elmt2}{idx2}.txt")

        # Save bond series (.dat)
        with open(self.fb, "w") as out_bond:
            for t, r in self.bond_lengths:
                out_bond.write(f"{t:>12.6f}{r:>12.6f}\n")

        # Plot series
        xs = [row[0] for row in self.bond_lengths]
        ys = [row[1] for row in self.bond_lengths]
        self.save_plots(1, xs, ys, "Simulation time (ps)", "Bond length (Å)", "Bond lengths")
        return True

    async def distribution_function(self) -> bool:
        """Histogram (%) of bond lengths and basic statistics."""
        if not self.bond_lengths:
            await self.warning_function("No bond lengths computed.")
            return False

        bonds_only = [r for _, r in self.bond_lengths]

        self.num_bins = max(1, int(self.max_r / self.bin_width))
        self.bin_centers = [self.bin_width * (i + 0.5) for i in range(self.num_bins)]
        counts = [0] * self.num_bins

        for r in bonds_only:
            idx = int(r / self.bin_width)
            if 0 <= idx < self.num_bins:
                counts[idx] += 1

        total = sum(counts)
        if total == 0:
            msg = (
                "The histogram is empty. Increase 'Maximum r for Distribution'.\n"
                f"Suggested minimum: {round(2 * max(bonds_only), 2)} Å."
            )
            await self.warning_function(msg)
            return False

        histogram_pct = [(c / total) * 100.0 for c in counts]
        self.histogram_sum = total
        self.histogram = histogram_pct

        # Stats
        avg = sum(bonds_only) / len(bonds_only)
        largest = max(bonds_only)
        smallest = min(bonds_only)
        var_val = variance(bonds_only) if len(bonds_only) > 1 else 0.0
        std_val = stdev(bonds_only) if len(bonds_only) > 1 else 0.0
        self.average_bond = avg

        self.stats = {
            "average_bond": avg,
            "variance": var_val,
            "std_dev": std_val,
            "largest_bond": largest,
            "smallest_bond": smallest,
            "bin_centers": self.bin_centers,
            "histogram_percentage": histogram_pct,
        }

        # Save distribution (.dat)
        with open(self.fd, "w") as f:
            for center, pct in zip(self.bin_centers, histogram_pct):
                f.write(f"{center:>12.6f}{pct:>12.6f}\n")

        # Plot histogram
        self.save_plots(
            2, self.bin_centers, self.histogram,
            "Bond length (Å)", "Bond length distribution (%)",
            "Bond length distribution function"
        )
        return True

    async def free_energy(self) -> bool:
        """Compute free energy from the histogram (kcal/mol)."""
        if not getattr(self, "histogram_sum", 0):
            self.stats = {}
            return False

        R = 0.001987204  # kcal/(mol*K)
        T = self.sim_temp
        self.free_energy_pairs = []  # store for CSV export

        for x_i, p_pct in zip(self.bin_centers, self.histogram):
            if p_pct <= 0.0:
                continue
            p = p_pct / 100.0
            FE = -R * T * np.log(p)
            self.free_energy_pairs.append([x_i, FE])

        if not self.free_energy_pairs:
            await self.warning_function("Free energy could not be computed from an empty distribution.")
            self.stats = {}
            return False

        xs = [x for x, _ in self.free_energy_pairs]
        ys = [y for _, y in self.free_energy_pairs]

        self.min_y = float(min(ys))
        self.min_y_idx = float(xs[ys.index(self.min_y)])

        # Save FE (.dat)
        with open(self.fe, "w") as f:
            for xi, FE in self.free_energy_pairs:
                f.write(f"{xi:>12.6f}{FE:>12.6f}\n")

        # Plot FE
        self.save_plots(
            3, xs, ys,
            "Bond length (Å)", "Bond free energy (kcal/mol)",
            "Bond length Free Energy Distribution"
        )
        return True

    def export_csv(self) -> None:
        """Optional CSV export with headers."""
        idx1, idx2 = self.atom_labels
        pair_tag = f"{self.elmt1}{idx1}_{self.elmt2}{idx2}"

        # 1) Time series
        csv_bond = os.path.join(self.output_dir, f"bond_{pair_tag}.csv")
        with open(csv_bond, "w") as f:
            f.write("time_ps,bond_length_A\n")
            for t, r in self.bond_lengths:
                f.write(f"{t:.6f},{r:.6f}\n")

        # 2) Distribution
        csv_dist = os.path.join(self.output_dir, f"bond_distribution_{pair_tag}.csv")
        with open(csv_dist, "w") as f:
            f.write("bond_length_A,distribution_pct\n")
            for x, pct in zip(self.bin_centers, self.histogram):
                f.write(f"{x:.6f},{pct:.6f}\n")

        # 3) Free energy
        if hasattr(self, "free_energy_pairs"):
            csv_fe = os.path.join(self.output_dir, f"free_energy_{pair_tag}.csv")
            with open(csv_fe, "w") as f:
                f.write("bond_length_A,free_energy_kcal_mol\n")
                for x, fe in self.free_energy_pairs:
                    f.write(f"{x:.6f},{fe:.6f}\n")

    def save_summary(self) -> None:
        """Write a consolidated summary once, at the end."""
        if not getattr(self, "stats", None):
            return

        idx1, idx2 = self.atom_labels
        with open(self.fs, "w") as f:
            f.write(f"Number of frames:............ {self.total_frame_number}\n")
            f.write(f"Number of atoms:............. {self.num_atoms}\n")
            f.write(f"Selected atoms:.............. {self.elmt1}{idx1}-{self.elmt2}{idx2}\n")
            f.write(f"The largest bond length:..... {self.stats['largest_bond']:.4f} Å\n")
            f.write(f"The smallest bond length:.... {self.stats['smallest_bond']:.4f} Å\n")
            f.write(f"Bond length average:......... {self.stats['average_bond']:.4f} Å\n")
            f.write(f"Bond length variance:........ {self.stats['variance']:.4f}\n")
            f.write(f"Bond length std deviation:... {self.stats['std_dev']:.4f}\n")
            f.write(
                f"Lowest free energy:.......... {self.min_y:.4f} kcal/mol "
                f"at {self.min_y_idx:.2f} Å\n"
            )

        with open(self.fs, "r") as f:
            self.multi_line_text.value = f.read()


class BondUI(BondAnalyser):
    """User Interface for the Bond Analysis module."""

    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget) -> None:
        # Window
        self.main_window = toga.Window(
            title="Bond Length Analysis from Molecular Dynamics Simulations",
            size=(720, 660),
        )

        # Styles
        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style   = Pack(margin=(5, 5), text_align=LEFT, width=240)
        input_style   = Pack(flex=1, margin=(5, 5))
        button_style  = Pack(margin=5, width=110)
        row_style     = Pack(direction="row", margin=(0, 0, 5, 0))
        col_style     = Pack(direction="column", margin=(0, 0, 5, 0))

        # Main container
        main_box = toga.Box(style=Pack(direction="column", margin=20))

        # Header row — no green background
        box_1 = toga.Box(style=Pack(direction="row", margin=(0, 0, 10, 0)))
        main_box.add(box_1)

        box_1a = toga.Box(style=Pack(width=420))
        box_1b = toga.Box(style=Pack(width=200))
        box_1c = toga.Box(style=Pack(width=100))
        box_1.add(box_1a)
        box_1.add(box_1b)
        box_1.add(box_1c)

        title_label = toga.Label("Bond Length Analysis", style=heading_style)
        self.progress_label = toga.Label(" ", style=heading_style)
        self.status_label = toga.Label(" ", style=Pack(margin=(0, 0, 0, 10)))  # prints output dir

        box_1a.add(title_label)
        box_1b.add(self.progress_label)
        box_1c.add(self.status_label)

        # Inputs
        input_fields = [
            ("Maximum r for Distribution Function:", "Enter maximum r value", "textInput_max_r"),
            ("Simulation Time Step (a.u.):", "Enter time step", "textInput_time_step"),
            ("Sampling Interval (frames):", "Enter sampling interval", "textInput_sampling_interval"),
            ("Simulation Temperature (K):", "Enter temperature", "textInput_temperature"),
            ("Atom Labels (e.g., 2 3):",
             "Enter two atom labels, separated by a space (e.g., 1 2)", "textInput_atom_labels"),
            ("Histogram Bin Width (Å):", "Enter bin width", "textInput_bin_width"),
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

        # Switches (use value=, not is_on=)
        switches_row = toga.Box(style=row_style)
        self.switch_show_plots = toga.Switch("Show plots at the end", value=True, style=Pack(margin=(4, 10, 4, 0)))
        self.switch_save_csv = toga.Switch("Save CSV outputs", value=True, style=Pack(margin=(4, 10, 4, 10)))
        switches_row.add(self.switch_show_plots)
        switches_row.add(self.switch_save_csv)

        # Progress bar
        progress_col = toga.Box(style=col_style)
        self.progress_bar = toga.ProgressBar(max=100)
        progress_col.add(self.progress_bar)

        main_box.add(switches_row)
        main_box.add(progress_col)

        # Help / Output
        self.multi_line_text = toga.MultilineTextInput(style=Pack(flex=1, margin=(10, 0), font_size=12))
        self.multi_line_text.value = HelpGqteaWin.help_bond_analysis
        main_box.add(self.multi_line_text)

        # Buttons
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
        if not hasattr(self, "trajec") or not self.trajec:
            await self.warning_function("No trajectory file selected.")
            return

        if not await self.bond_length(widget):
            return
        if not await self.distribution_function():
            return
        if not await self.free_energy():
            return
        self.save_summary()

        if self.switch_save_csv.value:
            self.export_csv()
        if self.switch_show_plots.value:
            self.display_plots()

        # PRINT the output directory (status label + multiline box)
        out_dir = getattr(self, "output_dir", os.getcwd())
        msg = f"Outputs saved to: {out_dir}"
        self.status_label.text = msg
        self.multi_line_text.value = (self.multi_line_text.value.rstrip() + "\n\n" + msg + "\n")

    def closeTopLevel(self, widget) -> None:
        self.main_window.close()


# Previous version of bond.py for reference
# import toga,os
# import numpy as np
# from toga.style import Pack
# from toga.style.pack import COLUMN, ROW, LEFT, CENTER
# from gqteaHelp import HelpGqteaWin
# from statistics import variance, stdev
# from framesCounter import FramesCounter
# from displayPlots import DisplayPlots


# class BondAnalyser(FramesCounter,DisplayPlots):

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
#                     if len(labels) != 2:
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
#         self.max_r = await read_input(self.textInput_max_r, "max_r", float)
#         if self.max_r is None:
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
#             f"Maximum r --> {self.max_r}\n"
#             f"Time step --> {self.dt}\n"
#             f"Interval sampling --> {self.sampling}\n"
#             f"Simulation temperature --> {self.sim_temp}\n"
#             f"Selected atom labels --> {self.atom_labels}\n"
#             f"Histogram bin width --> {self.bin_width}\n"
#         )
#         self.multi_line_text.value = update_text

#     async def bond_length(self, widget):
#         """Calculates the bond length between two given atom indices for each frame and
#         saves them in the output file like bond_C2_C1.dat. The first column is the simulation
#         time and the second column is the bond length.
#         """
           
#         self.atufs = 0.02418884326505  # Atomic time unit in femtoseconds.
#         self.bond_lengths = []

#         time_of_simulation = self.dt * self.sampling * self.atufs / 1000.0

#         idx1 = int(self.atom_labels[0])
#         idx2 = int(self.atom_labels[1])

#         try:
#             with open(self.trajec, "r") as f:
#                 self.progress_bar.start()
#                 frame_number = 0
#                 while True:
#                     title_line = f.readline()
#                     if not title_line:
#                         break
#                     comment_line = f.readline()
#                     if not comment_line:
#                         break

#                     atom_data = []
#                     for _ in range(self.num_atoms):
#                         line = f.readline()
#                         if not line:
#                             raise ValueError("Unexpected end of file")
#                         tokens = line.strip().split()
#                         if len(tokens) < 4:
#                             raise ValueError("Invalid atom line format")
#                         element = tokens[0]
#                         x, y, z = map(float, tokens[1:4])
#                         atom_data.append((element, x, y, z))

#                     # Get positions of the selected atoms
#                     try:
#                         elmt1, pos_x1, pos_y1, pos_z1 = atom_data[idx1 - 1]
#                         elmt2, pos_x2, pos_y2, pos_z2 = atom_data[idx2 - 1]
#                     except IndexError:
#                         await self.warning_function(
#                             f"Atom index out of range in frame {frame_number}"
#                         )
#                         return

#                     delta_x = pos_x2 - pos_x1
#                     delta_y = pos_y2 - pos_y1
#                     delta_z = pos_z2 - pos_z1

#                     bond_length = np.sqrt(delta_x ** 2 + delta_y ** 2 + delta_z ** 2)
#                     sim_time = frame_number * time_of_simulation
#                     self.bond_lengths.append([sim_time, bond_length])

#                     frame_number += 1

#                     if (frame_number % 400) == 0:
#                         progress_bar_increment = (frame_number / self.total_frame_number)*100
#                         self.progress_bar.value = progress_bar_increment

#                     if frame_number == self.total_frame_number:
#                         self.progress_bar.value = 100

#                 self.progress_bar.stop()

#         except Exception as e:
#             await self.warning_function(f"Error reading TRAJEC.xyz file: {e}")
#             return

#         # Prepare output files
#         self.elmt1 = elmt1 = atom_data[idx1 - 1][0]
#         self.elmt2 = elmt2 = atom_data[idx2 - 1][0]

#         self.fb = f"{self.output_dir}/bond_{elmt1}{idx1}_{elmt2}{idx2}.dat"
#         self.fe = f"{self.output_dir}/free_energy_{elmt1}{idx1}_{elmt2}{idx2}.dat"
#         self.fd = f"{self.output_dir}/bond_distribution_{elmt1}{idx1}_{elmt2}{idx2}.dat"
#         self.fs = f"{self.output_dir}/summary_{elmt1}{idx1}_{elmt2}{idx2}.txt"

#         # Write bond lengths to file
#         with open(self.fb, "w") as out_bond:
#             for sim_time, bond_length in self.bond_lengths:
#                 out_bond.write(f"{sim_time:>10.5f}{bond_length:>10.5f}\n")

#         # Write summary
#         with open(self.fs, "w") as f:
#             f.write(f"Number of frames:...................{self.total_frame_number}\n")
#             f.write(f"Number of atoms:....................{self.num_atoms}\n")
#             f.write(f"Selected atoms:.....................{elmt1}{idx1}--{elmt2}{idx2}\n")

#         # Plotting bond lengths
#         y = [sublist[1] for sublist in self.bond_lengths]
#         x = [sublist[0] for sublist in self.bond_lengths]

#         # Save the bond length plot file and plot
#         plot_xlabel = "Simulation time (ps)"
#         plot_ylabel = "Bond length (Å)"
#         plot_title = "Bond lengths"
#         plot_number = 1
#         self.save_plots(plot_number, x ,y, plot_xlabel, plot_ylabel, plot_title)
  
#     async def distribution_function(self):
#         """Calculate the distribution function of the bond length."""

#         self.num_bins = int(self.max_r / self.bin_width)
#         self.bin_centers = [self.bin_width * (i + 0.5) for i in range(self.num_bins)]
#         histogram = [0 for _ in range(self.num_bins)]

#         bonds_only = [sublist[1] for sublist in self.bond_lengths]

#         for bond in bonds_only:
#             bin_index = int(bond / self.bin_width)
#             if bin_index < self.num_bins:
#                 histogram[bin_index] += 1

#         self.histogram_sum = sum(histogram)

#         if self.histogram_sum > 0:
#             histogram_percentage = [(value / self.histogram_sum) * 100 for value in histogram]

#         else:
#             msg = (f"The histrogram is empty. Try to increse the maximum r for distribution\n"
#                   f" At least {round(2 * max(bonds_only), 0)} Å")
#             await self.warning_function(msg)
#             return                  

#         self.histogram = histogram_percentage

#         # Compute statistics
#         average_bond = sum(bonds_only) / len(bonds_only)
#         largest_bond = max(bonds_only)
#         smallest_bond = min(bonds_only)
#         var = variance(bonds_only)
#         std = stdev(bonds_only)
#         self.average_bond = average_bond

#         self.stats = {
#             "average_bond": average_bond,
#             "variance": var,
#             "std_dev": std,
#             "largest_bond": largest_bond,
#             "smallest_bond": smallest_bond,
#             "bin_centers": self.bin_centers,
#             "histogram_percentage": histogram_percentage,

#         }

#         # Save distribution to file
#         bond_dist_file = os.path.join(self.output_dir, "bond_distribution.dat")
#         with open(bond_dist_file, "w") as f:
#             for center, percent in zip(self.bin_centers, histogram_percentage):
#                 f.write(f"{center:>10.5f}{percent:>10.5f}\n")

#         # Save the bond length plot file
#         x = self.bin_centers
#         y = self.histogram
#         plot_xlabel = "Bond length (Å)"
#         plot_ylabel = "Bond length distribution (%)"
#         plot_title = "Bond length distribution function"
#         plot_number = 2
#         self.save_plots(plot_number,x ,y, plot_xlabel, plot_ylabel, plot_title)

#     async def free_energy(self):
#         """ Calculate the free energies based on the bond length distribution."""
#         if self.histogram_sum == 0:
#             self.stats = {}
#             return
        
#         temperature = self.sim_temp
#         R = 0.001987204  # Gas constant in kcal/(mol*K)
#         free_energy = []

#         for xi, yi in zip(self.bin_centers, self.histogram):
#             if yi > 0.0:
#                 FE = -R * temperature * np.log(yi / 100.0)
#                 free_energy.append([xi, FE])

#         y = [sublist[1] for sublist in free_energy]
#         x = [sublist[0] for sublist in free_energy]

#         self.min_y = min(y)
#         self.min_y_idx = x[y.index(self.min_y)]

#         # Save free energy to file
#         fe_file = os.path.join(self.output_dir, "bond_free_energy.dat")
#         with open(fe_file, "w") as f:
#             for xi, FE in free_energy:
#                 f.write(f"{xi:>10.5f}{FE:>10.5f}\n")

#         #Save the temporary plot file
#         plot_xlabel = "Bond Length (Å)"
#         plot_ylabel = "Bond Free Energy (kcal/mol)"
#         plot_title = "Bond length Free Energy Distribution"
#         plot_number = 3
#         self.save_plots(plot_number, x ,y, plot_xlabel, plot_ylabel, plot_title)

#     def save_summary(self):
#         """Saves the summary of the analysis."""
#         if not self.stats:
#             return
        
#         idx1, idx2 = self.atom_labels
#         summary_file = os.path.join(
#             self.output_dir, f"summary_{self.elmt1}{idx1}_{self.elmt2}{idx2}.txt"
#         )

#         with open(summary_file, "w") as f:
#             f.write(f"Number of frames:............ {self.total_frame_number}\n")
#             f.write(f"Number of atoms:............. {self.num_atoms}\n")
#             f.write(f"Selected atoms:.............. {self.elmt1}{idx1}-{self.elmt2}{idx2}\n")
#             f.write(f"The largest bond length:..... {self.stats['largest_bond']:.4f} Å\n")
#             f.write(f"The smallest bond length:.... {self.stats['smallest_bond']:.4f} Å\n")
#             f.write(f"Bond length average:......... {self.stats['average_bond']:.4f} Å\n")
#             f.write(f"Bond length variance:........ {self.stats['variance']:.4f}\n")
#             f.write(f"Bond length std deviation:... {self.stats['std_dev']:.4f}\n")
#             f.write(f"Lowest free energy:.......... {self.min_y:.4f} kcal/mol "
#                 f"at {self.min_y_idx:.2f} Å\n"
#             )

#         with open(summary_file, "r") as f:
#             summary_content = f.read()
#             self.multi_line_text.value = summary_content


# class BondUI(BondAnalyser):
#     def __init__(self, *args):
#         self.layout_main_window(*args)

#     def layout_main_window(self, widget):
#         # Create the main window
#         self.main_window = toga.Window(
#             title="Bond Length Analysis from Molecular Dynamics Simulations",
#             size=(700, 600),
#         )

#         # Define common styles
#         heading_style = Pack(font_size=18, font_weight="bold", padding=(0, 0, 10, 0))
#         label_style = Pack(padding=(5, 5), text_align=LEFT, width=200)
#         input_style = Pack(flex=1, padding=(5, 5))
#         button_style = Pack(padding=5, width=100)
#         box_style = Pack(direction=ROW, alignment=CENTER, padding=(0, 0, 5, 0))

#         # Main container
#         main_box = toga.Box(style=Pack(direction=COLUMN, padding=20))

#         box_1 = toga.Box(style=Pack(direction=ROW, 
#                         alignment=CENTER, padding=(0,0,10,0),
#                         background_color = "green")
#         )

#         main_box.add(box_1)
#         box_1a = toga.Box(style=Pack(width=400))
#         box_1b = toga.Box(style=Pack(width=150))
#         box_1c = toga.Box(style=Pack(width=100))
#         box_1.add(box_1a)
#         box_1.add(box_1b)
#         box_1.add(box_1c)

#         # Title
#         title_label = toga.Label("Bond Length Analysis", 
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
#                 "Maximum r for Distribution Function:",
#                 "Enter maximum r value",
#                 "textInput_max_r",
#             ),
            
#             (
#                 "Simulation Time Step (a.u.):", 
#                 "Enter time step", 
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
#                 "Atom Labels (e.g., 2 3):",
#                 "Enter two atom labels, separated by a space (e.g., 1 2)",
#                 "textInput_atom_labels",
#             ),
            
#             (
#                 "Histogram Bin Width (Å):",
#                 "Enter bin width",
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
#             "Browse", on_press=self.frames_counter, style=button_style
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
#         self.multi_line_text.value = HelpGqteaWin.help_bond_analysis
#         main_box.add(self.multi_line_text)

#         # Buttons at the bottom
#         button_box = toga.Box(
#             style=Pack(direction=ROW, alignment=CENTER, padding_top=10)
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
        
#         await self.bond_length(widget)
#         await self.distribution_function()
#         await self.free_energy()
#         self.save_summary()
#         self.display_plots()

#     def closeTopLevel(self, widget):
#         self.main_window.close()

