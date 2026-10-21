import math

import toga
import numpy as np
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from help import HelpGqteaWin
from framesCounter import FramesCounter
from displayPlots import DisplayPlots

ATU_FS = 0.02418884326505  # fs per atomic unit of time

TIME_STEP_UNIT_AU = "a.u."
TIME_STEP_UNIT_FS = "fs"
TIME_STEP_UNITS = (TIME_STEP_UNIT_AU, TIME_STEP_UNIT_FS)   # a.u. first: the CPMD convention

TOLERANCE_UNIT_FRAMES = "frames"
TOLERANCE_UNIT_FS = "fs"
TOLERANCE_UNITS = (TOLERANCE_UNIT_FRAMES, TOLERANCE_UNIT_FS)
DEFAULT_TOLERANCE = 0

# t / dt is computed in floating point, so a tolerance that is an exact
# multiple of the frame interval can come out a hair below the whole number
# (3 * dt / dt = 2.9999999999999996); this much is added before rounding down.
FRAME_ROUNDING_SLACK = 1e-9


def frame_interval_fs(time_step: float, time_step_unit: str, sampling: int) -> float:
    """Time between two sampled frames, in fs: simulation time step (fs) x sampling interval."""
    if time_step_unit == TIME_STEP_UNIT_AU:
        step_fs = time_step * ATU_FS
    elif time_step_unit == TIME_STEP_UNIT_FS:
        step_fs = time_step
    else:
        raise ValueError(f"Unknown time step unit '{time_step_unit}'; use {' or '.join(TIME_STEP_UNITS)}.")
    return step_fs * sampling


def tolerance_in_frames(value: float, unit: str, interval_fs: float) -> int:
    """The tolerance as a whole number of sampled frames.

    The tolerance is the largest number of consecutive sampled frames a
    molecule may spend outside the shell without that excursion counting as
    an exit. In frames the value must already be a whole number. In fs it is
    converted with ``frames = time / interval_fs`` and **rounded down**: an
    excursion of n frames lasts n x interval_fs, and it is tolerated when that
    fits within the time given, i.e. n <= time / interval_fs. So n frames and
    n x interval_fs fs always give the same result, and a time between two
    multiples of the interval acts as the lower one (e.g. 2.5 frames -> 2).
    """
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"The tolerance must be zero or a positive number, not {value:g}.")
    if unit == TOLERANCE_UNIT_FRAMES:
        if not value.is_integer():
            raise ValueError(
                f"A tolerance in frames must be a whole number, not {value:g}. "
                "Choose fs as the unit to give it as a time."
            )
        return int(value)
    if unit == TOLERANCE_UNIT_FS:
        if not (math.isfinite(interval_fs) and interval_fs > 0):
            raise ValueError("The time between sampled frames must be greater than zero.")
        return int(math.floor(value / interval_fs + FRAME_ROUNDING_SLACK))
    raise ValueError(f"Unknown tolerance unit '{unit}'; use {' or '.join(TOLERANCE_UNITS)}.")


class MRTAnalyzer(FramesCounter,DisplayPlots):

    @staticmethod
    def apply_tolerance_to_occupancy(occupancy: np.ndarray, tolerance_frames: int) -> np.ndarray:
        if tolerance_frames <= 0:
            return occupancy.copy()

        corrected = occupancy.copy()
        n_objects, n_frames = corrected.shape

        for obj in range(n_objects):
            arr = corrected[obj]
            i = 0
            while i < n_frames:
                if arr[i] == 1:
                    i += 1
                    continue

                start_zero = i
                while i < n_frames and arr[i] == 0:
                    i += 1
                end_zero = i - 1
                zero_len = end_zero - start_zero + 1

                left_is_one = start_zero - 1 >= 0 and arr[start_zero - 1] == 1
                right_is_one = i < n_frames and arr[i] == 1

                if left_is_one and right_is_one and zero_len <= tolerance_frames:
                    arr[start_zero:end_zero + 1] = 1

        return corrected

    @staticmethod
    def count_exits(occupancy_row, tolerance_frames: int) -> int:
        """Exits from the shell in one atom's tolerance-corrected occupancy (1 = inside).

        An exit is an inside -> outside step whose excursion lasts more than
        ``tolerance_frames`` sampled frames. Excursions between two in-shell
        frames that were no longer have already been filled in by
        :meth:`apply_tolerance_to_occupancy`, so every such step left is an
        exit. An excursion still under way when the trajectory ends counts
        only if it already lasts more than the tolerance. With a tolerance of
        0 every inside -> outside step is an exit.
        """
        row = [int(value) for value in occupancy_row]
        n_frames = len(row)
        exits = 0
        i = 0
        while i < n_frames - 1:
            if row[i] == 1 and row[i + 1] == 0:
                j = i + 1
                while j < n_frames and row[j] == 0:
                    j += 1
                if j < n_frames or (j - i - 1) > tolerance_frames:
                    exits += 1
                i = j
            else:
                i += 1
        return exits

    async def warning_function(self, msg):
        await self.main_window.dialog(toga.InfoDialog("Error", f"{msg}"))

    async def read_params(self, widget):

        async def read_input(text_input, field_name, expected_type):
            value = text_input.value.strip()
            if not value:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Please input a valid value for {field_name}.")
                )
                return None
            try:
                if expected_type == list:
                    # Assuming atom labels are space-separated integers
                    labels = [int(label) for label in value.split()]
                    return labels
                else:
                    return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        if not hasattr(self, "trajec") or not self.trajec:
            await self.warning_function(
                "Please select a TRAJEC.xyz file before starting the MRT calculation."
            )
            return False

        if not hasattr(self, "total_frame_number") or self.total_frame_number <= 0:
            await self.warning_function(
                "Please load the trajectory with Browse before starting the MRT calculation."
            )
            return False

        # Read and validate inputs
        self.inner_radius = 0.0
        self.outer_radius = await read_input(self.textInput_outer_radius, "radius for the outer shell", float)
        if self.outer_radius is None:
            return False

        self.center_label = await read_input(self.textInput_center_label, "atom label for the center of the shell", int)
        if self.center_label is None:
            return False

        self.atom_symbol = await read_input(self.textInput_atom_symbol, "atom symbol for mean residence time", str)
        if self.atom_symbol is None:
            return False
        self.atom_symbol = self.atom_symbol.strip()

        self.time_step = await read_input(self.textInput_time_step, "simulation time step", float)
        if self.time_step is None:
            return False
        self.time_step_unit = self.selection_time_step_unit.value or TIME_STEP_UNIT_AU

        self.sampling = await read_input(self.textInput_sampling, "sampling interval between two colected frames", int)
        if self.sampling is None:
            return False

        # Blank means the default tolerance, 0: every exit is counted.
        self.tolerance_unit = self.selection_tolerance_unit.value or TOLERANCE_UNIT_FRAMES
        tolerance_text = self.textInput_tolerance.value.strip()
        try:
            self.tolerance_value = float(tolerance_text) if tolerance_text else float(DEFAULT_TOLERANCE)
        except ValueError:
            await self.warning_function(
                f"Invalid format for the tolerance: '{tolerance_text}' is not a number."
            )
            return False

        self.exclude_atom_list = await read_input(self.textInput_atom_labels, "atom labels list to exclude", list)
        if self.exclude_atom_list is None:
            return False

        # Adjust for zero-based indexing (assuming user input is one-based)
        if self.exclude_atom_list != [0]:
            self.exclude_atom_list = [label for label in self.exclude_atom_list]
        else:
            self.exclude_atom_list = []

        if self.outer_radius < 0:
            await self.warning_function("The cutoff radius must be non-negative.")
            return False

        if self.center_label < 1 or self.center_label > self.num_atoms:
            await self.warning_function(
                f"The shell center atom label must be between 1 and {self.num_atoms}."
            )
            return False

        if self.time_step <= 0:
            await self.warning_function("The simulation time step must be greater than zero.")
            return False

        if self.sampling <= 0:
            await self.warning_function("The sampling interval must be greater than zero.")
            return False

        self.frame_interval_fs = frame_interval_fs(self.time_step, self.time_step_unit, self.sampling)
        try:
            self.tolerance_frames = tolerance_in_frames(
                self.tolerance_value, self.tolerance_unit, self.frame_interval_fs
            )
        except ValueError as e:
            await self.warning_function(str(e))
            return False

        invalid_excluded_labels = [
            label for label in self.exclude_atom_list if label < 1 or label > self.num_atoms
        ]
        if invalid_excluded_labels:
            await self.warning_function(
                f"Excluded atom labels must be between 1 and {self.num_atoms}."
            )
            return False

        self.multi_line_text.value = self.parameters_text() + "\nRunning the MRT calculation...\n"
        return True

    def tolerance_text(self) -> str:
        """The tolerance as entered and as applied, e.g. '25 fs -> 2 frames (24.19 fs)'."""
        applied = (f"{self.tolerance_frames} frame{'s' if self.tolerance_frames != 1 else ''} "
                   f"({self.tolerance_frames * self.frame_interval_fs:.6g} fs)")
        if self.tolerance_unit == TOLERANCE_UNIT_FS:
            return f"{self.tolerance_value:g} fs -> {applied}"
        return applied

    def parameters_text(self) -> str:
        """The parameters of the run, as shown in the text box."""
        exits = ("every exit from the shell is counted" if self.tolerance_frames == 0 else
                 f"excursions of up to {self.tolerance_frames} sampled frames are not counted as exits")
        return (
            f"Cutoff radius --> {self.outer_radius} Å\n"
            f"Atom label for the center of the shell --> {self.center_label}\n"
            f"Atom symbol for mean residence time --> {self.atom_symbol}\n"
            f"Simulation time step --> {self.time_step:g} {self.time_step_unit}\n"
            f"Sampling interval between two collected frames --> {self.sampling}\n"
            f"Time between sampled frames --> {self.frame_interval_fs:.6g} fs\n"
            f"Tolerance --> {self.tolerance_text()}: {exits}\n"
            f"Atom labels list to exclude --> {self.exclude_atom_list}\n"
        )

    #Get atom symbol at the center of the shell
    def shell_center_symbol(self):
        with open(self.trajec, "r") as f:
            first_line = f.readline().split()
            if not first_line:
                raise ValueError("The trajectory file is empty.")
            num_atoms = int(first_line[0])
            f.readline()
            shell_center_symb = None
            for atom in range(num_atoms):
                line = f.readline().split()
                if len(line) < 4:
                    raise ValueError("Invalid coordinate line while reading the shell center atom.")
                if (atom + 1) == self.center_label:
                    shell_center_symb = line[0]
                    break
        if shell_center_symb is None:
            raise ValueError("The selected shell center atom label was not found in the trajectory.")
        return shell_center_symb

    def _validate_analysis_state(self):
        required_attributes = [
            "trajec",
            "output_dir",
            "inner_radius",
            "outer_radius",
            "center_label",
            "atom_symbol",
            "time_step",
            "time_step_unit",
            "sampling",
            "tolerance_frames",
            "exclude_atom_list",
            "total_frame_number",
        ]
        missing_attributes = [
            attr for attr in required_attributes if not hasattr(self, attr)
        ]
        if missing_attributes:
            missing_attrs = ", ".join(missing_attributes)
            raise ValueError(
                f"Missing required analysis data: {missing_attrs}. Please reload the inputs and trajectory."
            )

    def calc_mrt(self):
        """
        Calculate the mean residence time (MRT) for atoms within a specified shell.

        This function reads trajectory data to identify and store the labels of atoms
        that reside within the defined shell at least once during the simulation.

        Returns:
            None: The results are stored as class attributes, including a list of atom labels
            that enter the shell and the number of frames processed.
        """
        self._validate_analysis_state()
        inner_r = self.inner_radius
        outer_r = self.outer_radius
        acs = self.center_label  # Atom index at the center of the shell
        symbol = self.atom_symbol  # Element symbol for MRT calculation
        # Time between two sampled frames, in ps (time step in a.u. or fs x sampling interval)
        sim_time = frame_interval_fs(self.time_step, self.time_step_unit, self.sampling) / 1000.0
        self.sim_time = sim_time  # To be used in other functions
        shell_atom_list = []  # Labels of atoms that enter the shell at least once

        num_frames = 0  # Counter for the number of frames processed
        with open(self.trajec, "r") as f:
            self.progress_bar.start()

            while True:
                tmp_frame = []  # Temporary list to store atoms of the current frame

                num_atoms = f.readline().split()
                if not num_atoms:
                    break
                num_atoms = int(num_atoms[0])  # Read the number of atoms per frame

                line = f.readline()  # Read the comment line
                if not line:
                    break

                for _ in range(num_atoms):
                    line = f.readline().split()
                    if len(line) != 4:
                        raise ValueError(
                            f"Invalid atom record found in frame {num_frames + 1}."
                        )
                    tmp_frame.append(
                        [line[0], float(line[1]), float(line[2]), float(line[3])]
                    )

                coords1 = np.array(
                    [tmp_frame[acs - 1][1], tmp_frame[acs - 1][2], tmp_frame[acs - 1][3]]
                )

                for n in range(num_atoms):
                    if (
                        n != (acs - 1)
                        and tmp_frame[n][0] == symbol
                        and (n + 1) not in self.exclude_atom_list
                    ):
                        coords2 = np.array([tmp_frame[n][1], tmp_frame[n][2], tmp_frame[n][3]])
                        vec = coords2 - coords1
                        vec_norm = np.linalg.norm(vec)
                        if inner_r <= vec_norm <= outer_r:
                            if (n + 1) not in shell_atom_list:
                                shell_atom_list.append(n + 1)

                num_frames += 1

                if num_frames % 300 == 0:
                    self.progress_bar.value = min(
                        100, (num_frames / self.total_frame_number) * 100
                    )

        self.progress_bar.value = 100
        self.progress_bar.stop()

        self.shell_atom_list = shell_atom_list


        with open(self.trajec, "r") as f:
            with open(f"{self.output_dir}/mrt.dat", "w") as mrt:
                header = [f"{self.atom_symbol}{s}" for s in shell_atom_list]
                mrt.write(f"Frame  Time(ps)  {'    '.join(header)}\n")

                self.progress_bar.start()

                num_frames = 0
                occupancy_rows = []
                while True:
                    num_atoms_line = f.readline().split()
                    if not num_atoms_line:
                        break
                    num_atoms = int(num_atoms_line[0])

                    tmp_frame = []
                    tmp_list = [0 for i in shell_atom_list]
                    line = f.readline()  # Read the comment line
                    if not line:
                        break

                    for n in range(num_atoms):
                        line = f.readline().split()
                        if len(line) != 4:
                            raise ValueError(
                                f"Invalid atom record found in frame {num_frames + 1}."
                            )
                        tmp_frame.append(
                            [line[0], float(line[1]), float(line[2]), float(line[3])]
                        )

                    coords1 = np.array(
                        [tmp_frame[acs - 1][1], tmp_frame[acs - 1][2], tmp_frame[acs - 1][3]]
                    )

                    for index, n in enumerate(shell_atom_list):
                        coords2 = np.array(
                            [tmp_frame[n - 1][1], tmp_frame[n - 1][2], tmp_frame[n - 1][3]]
                        )
                        vec = coords2 - coords1
                        vec_norm = np.linalg.norm(vec)
                        if vec_norm >= inner_r and vec_norm <= outer_r:
                            tmp_list[index] = 1

                    num_frames += 1
                    occupancy_rows.append(tmp_list)

                    if num_frames % 300 == 0:
                        self.progress_bar.value = min(
                            100, (num_frames / self.total_frame_number) * 100
                        )

                if shell_atom_list and occupancy_rows:
                    occupancy_matrix = np.array(occupancy_rows, dtype=int).T
                    corrected_occupancy = self.apply_tolerance_to_occupancy(
                        occupancy_matrix,
                        self.tolerance_frames,
                    )
                    corrected_rows = corrected_occupancy.T.tolist()
                else:
                    corrected_rows = occupancy_rows

                for frame_index, tmp_list in enumerate(corrected_rows, start=1):
                    Time = sim_time * frame_index
                    tmp_str_list = [str(s) for s in tmp_list]
                    mrt.write(f"{frame_index:<6d} {Time:>8.4f}   {'      '.join(tmp_str_list)}\n")

                self.progress_bar.value = 100
                self.progress_bar.stop()

    async def write_mrt_total(self):
        with open(f"{self.output_dir}/mrt_total.dat", "w") as f:
            with open(f"{self.output_dir}/mrt.dat", "r") as mrt:
                # Read and discard the first line (header)
                mrt.readline()
                
                # Now iterate through the remaining lines and write them
                f.write("Time(ps)  Total coordination number\n")
                for line in mrt: # Iterating directly over the file object is memory efficient
                    line = line.split()
                    try:
                        # Check if the line has enough columns to avoid IndexError
                        if len(line) > 2:  # Ensure there are enough columns
                            time = float(line[1])
                            total_coordination_number = sum(int(s) for s in line[2:])
                            f.write(f"{time:>8.4f}  {total_coordination_number:>10d}\n")
                    except ValueError as e:
                        # Handle any conversion errors (e.g., if the line is not formatted as expected)
                        # If the line doesn't have enough columns, just write it as is
                        await self.main_window.dialog(
                            toga.ErrorDialog(
                                "Error",
                                f"Error processing line: {line} - {e}",
                            )
                            )

    def mrt_summary(self):
        symbol = self.atom_symbol # Atom symbol for mrt
        mrt_dat = []
        with open(f"{self.output_dir}/mrt.dat", "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.split()
                line = line[2:]
                line = [int(s) for s in line]
                mrt_dat.append(line)

        with open(f"{self.output_dir}/mrt_summary.dat", "w") as f:
            rt = [0 for i in range(len(self.shell_atom_list))]  # Residence time for each atom
            # mrt.dat already holds the tolerance-corrected occupancy
            exits = [
                self.count_exits([row[k] for row in mrt_dat], self.tolerance_frames)
                for k in range(len(self.shell_atom_list))
            ]

            for n in range(len(mrt_dat)):
                for k in range(len(self.shell_atom_list)):
                    if mrt_dat[n][k] == 1:
                        rt[k] += 1

            total_time_sim = self.sim_time * self.total_frame_number  # Total time simulation
            f.write(f"Total simulation time.........{total_time_sim:>10.4f}\n")
            f.write(f"Time between frames (fs)......{self.sim_time * 1000.0:>10.4f}\n")
            f.write(f"Tolerance (frames)............{self.tolerance_frames:>10d}\n")
            f.write(f"Tolerance (fs)................{self.tolerance_frames * self.sim_time * 1000.0:>10.4f}\n\n")
            self.mrt_results = []
            for n in range(len(rt)):
                tot_rt = rt[n] * self.sim_time
                if exits[n] != 0:
                    tmp_rt = rt[n] * self.sim_time / exits[n]
                else:
                    tmp_rt = rt[n] * self.sim_time  # never seen to leave: a lower bound
                self.mrt_results.append((self.shell_atom_list[n], tot_rt, tmp_rt, exits[n]))

                f.write(
                    f"{symbol}{self.shell_atom_list[n]} total rt.................{tot_rt:>10.4f}\n"
                )
                f.write(
                    f"{symbol}{self.shell_atom_list[n]} mrt......................{tmp_rt:>10.4f}\n"
                )
                f.write(
                    f"{symbol}{self.shell_atom_list[n]} exits....................{exits[n]:>10d}\n\n"
                )

        center_atom_symbol = self.shell_center_symbol()
        self.multi_line_text.value = self.results_text(center_atom_symbol)

    def results_text(self, center_atom_symbol: str) -> str:
        """The closing message: parameters, per-atom results (times in ps) and the output files."""
        lines = [
            f"MEAN RESIDENCE TIME ANALYSIS FOR {center_atom_symbol}{self.center_label} SITE FINISHED SUCCESSFULLY!",
            "",
            self.parameters_text(),
            "Results (times in ps; mrt = total rt / exits, or total rt when no exit was seen):",
        ]
        if not self.mrt_results:
            lines.append(f"  No {self.atom_symbol} atom entered the shell.")
        for label, tot_rt, mrt, exits in self.mrt_results:
            lines.append(f"  {self.atom_symbol}{label:<6} total rt {tot_rt:>10.4f}   mrt {mrt:>10.4f}   exits {exits:>6d}")
        lines += ["", f"Files written to {self.output_dir}: mrt.dat, mrt_total.dat, mrt_summary.dat"]
        return "\n".join(lines) + "\n"


    async def mrt_plot(self,widget):
        x_axis = []
        y_axis = []
        with open(f"{self.output_dir}/mrt.dat", "r") as f:
            lines = f.readlines()[1:]  # Skip header line

            for line in lines:
                line = line.split()  # Ensure line is split into columns

                try:
                    x_axis.append(float(line[1]))  # Convert time to float
                    data = [int(s) for s in line[2:]]  # Convert data points to int
                    y_axis.append(sum(data))  # Sum data for y-axis value
                except (ValueError, IndexError) as e:
                    await self.main_window.dialog(
                        toga.ErrorDialog("Error", f"Error processing line: {line} - {e}")
                    )

        # Plot labels and title
        plot_xlabel = "Time (ps)"
        plot_ylabel = "Coordination number"
        plot_title = "Mean Residence Time"

        # Save the plot with provided data
        self.save_plots(1, x_axis, y_axis, plot_xlabel, plot_ylabel, plot_title)


class ClassicalMeanResidenceTimeUI(MRTAnalyzer):
    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self,widget):
        # Create the main window
        self.main_window = toga.Window(
            title="Mean Residence Time (MRT)",
            size=(700, 600),
        )

        # Define common styles
        heading_style = Pack(
            font_size=18, font_weight="bold", text_align=LEFT, margin=(0, 0, 10, 0)
        )
        label_style = Pack(
            margin=(0, 0, 5, 5), text_align=LEFT, width=200
        )
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=100)
        box_style = Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0))

        # Main container
        main_box = toga.Box(style=Pack(direction=COLUMN, margin=20))

        box_title = toga.Box(
            style=Pack(
                direction=ROW,
                align_items=CENTER,
                margin=(0, 0, 10, 0),
            )
        )

        sub_box_title = toga.Box(style=Pack(width=410))
        sub_box_empty = toga.Box(style=Pack(width=140))
        sub_box_progress = toga.Box(style=Pack(width=100))

        box_title.add(sub_box_title)
        box_title.add(sub_box_empty)
        box_title.add(sub_box_progress)

        # Title
        title_label = toga.Label("Mean Residence Time (MRT)", style=heading_style)
        empty_label = toga.Label(" ", style=heading_style)
        self.progress_label = toga.Label(" ", style=heading_style)

        sub_box_title.add(title_label)
        sub_box_empty.add(empty_label)
        sub_box_progress.add(self.progress_label)

        main_box.add(box_title)

        # Input fields with labels
        input_fields = [
            (
                "Cutoff radius (\u00c5):",
                "Enter the cutoff radius",
                "textInput_outer_radius",
            ),
            (
                "Atom label for the shell center:",
                "Enter the atom label for the center of the shell",
                "textInput_center_label",
            ),
            (
                "Atom symbol for MRT:",
                "Enter the atom symbol for mean residence time",
                "textInput_atom_symbol",
            ),
            (
                "Atom label list to exclude:",
                "Enter the list of atom labels to exclude from the MRT analysis",
                "textInput_atom_labels",
            ),
        ]

        for label_text, placeholder, attr_name in input_fields:
            box = toga.Box(style=box_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            box.add(label)
            box.add(text_input)
            main_box.add(box)

        time_box = toga.Box(style=box_style)
        time_label = toga.Label("Simulation time step:", style=label_style)
        self.textInput_time_step = toga.TextInput(
            placeholder="Time step",
            style=Pack(width=90, margin=(5, 5)),
        )
        self.selection_time_step_unit = toga.Selection(
            items=list(TIME_STEP_UNITS), style=Pack(width=60, margin=(5, 5))
        )
        self.selection_time_step_unit.value = TIME_STEP_UNIT_AU
        sampling_label = toga.Label(
            "Sampling interval (steps):",
            style=Pack(margin=(0, 0, 5, 10), text_align=LEFT, width=150),
        )
        self.textInput_sampling = toga.TextInput(
            placeholder="MD steps per frame",
            style=Pack(width=110, margin=(5, 5)),
        )
        time_box.add(time_label)
        time_box.add(self.textInput_time_step)
        time_box.add(self.selection_time_step_unit)
        time_box.add(sampling_label)
        time_box.add(self.textInput_sampling)
        main_box.add(time_box)

        tolerance_box = toga.Box(style=box_style)
        self.tolerance_label = toga.Label("", style=label_style)
        self.textInput_tolerance = toga.TextInput(value=str(DEFAULT_TOLERANCE), style=input_style)
        self.selection_tolerance_unit = toga.Selection(
            items=list(TOLERANCE_UNITS), style=Pack(width=80, margin=(5, 5))
        )
        self.selection_tolerance_unit.value = TOLERANCE_UNIT_FRAMES
        # Wired after the value is set, then run once to write the label.
        self.selection_tolerance_unit.on_change = self.on_tolerance_unit_change
        self.on_tolerance_unit_change(self.selection_tolerance_unit)
        tolerance_box.add(self.tolerance_label)
        tolerance_box.add(self.textInput_tolerance)
        tolerance_box.add(self.selection_tolerance_unit)
        main_box.add(tolerance_box)

        # File selection button
        file_box = toga.Box(style=box_style)
        file_label = toga.Label("Select input file:", style=label_style)
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select TRAJEC.xyz file",
            style=input_style,
        )
        browse_button = toga.Button(
            "Browse", on_press=self.frames_counter, style=button_style
        )
        file_box.add(file_label)
        file_box.add(self.textInput_file)
        file_box.add(browse_button)

        # Progress bar
        self.progress_bar = toga.ProgressBar(max=100)
        main_box.add(file_box)
        main_box.add(self.progress_bar)

        # Multi-line text input for output
        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(5, 0), font_size=11)
        )
        self.multi_line_text.value = HelpGqteaWin.help_mrt
        main_box.add(self.multi_line_text)

        # Buttons at the bottom
        button_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin_top=5)
        )
        self.btn_execute = toga.Button(
            "MRT calculation", style=button_style, on_press=self.workflow
        )
        self.btn_close = toga.Button(
            "Close", style=button_style, on_press=self.closeTopLevel
        )
        button_box.add(self.btn_execute)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        # Set the content of the main window
        self.main_window.content = main_box
        self.main_window.show()

    TOLERANCE_LABELS = {
        TOLERANCE_UNIT_FRAMES: ("Tolerance (frames):",
                                "Frames outside the shell not counted as an exit (0 = count every exit)"),
        TOLERANCE_UNIT_FS: ("Tolerance (fs):",
                            "Time outside the shell not counted as an exit (0 = count every exit)"),
    }

    def on_tolerance_unit_change(self, widget, **kwargs):
        """Name the tolerance field after the chosen unit."""
        label, placeholder = self.TOLERANCE_LABELS.get(
            self.selection_tolerance_unit.value, self.TOLERANCE_LABELS[TOLERANCE_UNIT_FRAMES]
        )
        self.tolerance_label.text = label
        self.textInput_tolerance.placeholder = placeholder

    async def workflow(self, widget):
        params_are_valid = await self.read_params(widget)
        if not params_are_valid:
            return

        try:
            self.calc_mrt()
            self.mrt_summary()
            await self.mrt_plot(widget)
            self.display_plots()
            await self.write_mrt_total()
            # mrt_summary has left the results in the text box
        except Exception as e:
            self.progress_bar.stop()
            await self.warning_function(str(e))

    def closeTopLevel(self, widget):
        self.main_window.close()
