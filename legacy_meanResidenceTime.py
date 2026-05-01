import toga
import numpy as np
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from help import HelpGqteaWin
from framesCounter import FramesCounter
from displayPlots import DisplayPlots


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

        self.sampling = await read_input(self.textInput_sampling, "sampling interval between two colected frames", int)
        if self.sampling is None:
            return False

        self.tolerance_frames = await read_input(self.textInput_tolerance_frames, "tolerance frames", int)
        if self.tolerance_frames is None:
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

        if self.tolerance_frames < 0:
            await self.warning_function("Tolerance frames cannot be negative.")
            return False

        invalid_excluded_labels = [
            label for label in self.exclude_atom_list if label < 1 or label > self.num_atoms
        ]
        if invalid_excluded_labels:
            await self.warning_function(
                f"Excluded atom labels must be between 1 and {self.num_atoms}."
            )
            return False

        update_text = (
            f"Cutoff radius --> {self.outer_radius}\n"
            f"Atom label for the center of the shell --> {self.center_label}\n"
            f"Atom symbol for mean residence time --> {self.atom_symbol}\n"
            f"Time step --> {self.time_step}\n"
            f"Sampling interval between two colected frames --> {self.sampling}\n"
            f"Tolerance frames --> {self.tolerance_frames}\n"
            f"Atom labels list to exclude --> {self.exclude_atom_list}\n"
        )
        self.multi_line_text.value = update_text
        return True

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
        atufs = 0.02418884326505
        inner_r = self.inner_radius
        outer_r = self.outer_radius
        step = self.time_step  # Simulation time step
        sampling = self.sampling  # Number of frames skipped between collected frames
        acs = self.center_label  # Atom index at the center of the shell
        symbol = self.atom_symbol  # Element symbol for MRT calculation
        sim_time = (atufs * step * sampling) / 1000.0  # Simulation time in ps
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
            num_exchange = [0 for i in range(len(self.shell_atom_list))]

            for n in range(len(mrt_dat)):
                for k in range(len(self.shell_atom_list)):
                    if mrt_dat[n][k] == 1:
                        rt[k] += 1
                    if n < len(mrt_dat) - 1:
                        if mrt_dat[n][k] == 0 and mrt_dat[n + 1][k] == 1:
                            num_exchange[k] += 1

            total_time_sim = self.sim_time * self.total_frame_number  # Total time simulation
            f.write(f"Total simulation time.........{total_time_sim:>10.4f}\n")
            for n in range(len(rt)):
                tot_rt = rt[n] * self.sim_time
                if num_exchange[n] != 0:
                    tmp_rt = rt[n] * self.sim_time / num_exchange[n]
                else:
                    tmp_rt = rt[n] * self.sim_time

                f.write(
                    f"{symbol}{self.shell_atom_list[n]} total rt.................{tot_rt:>10.4f}\n"
                )
                f.write(
                    f"{symbol}{self.shell_atom_list[n]} mrt......................{tmp_rt:>10.4f}\n"
                )
                f.write(
                    f"{symbol}{self.shell_atom_list[n]} exchange.................{num_exchange[n]:>10d}\n\n"
                )

        center_atom_symbol = self.shell_center_symbol()
        self.multi_line_text.value = f"\n\n MEAN RESIDENCE TIME ANALYSIS FOR {center_atom_symbol}{self.center_label} SITE FINISHED SUCCSSEFULLY!\n"


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


class LegacyMeanResidenceTimeUI(MRTAnalyzer):
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
            placeholder="Enter the simulation time step",
            style=Pack(width=110, margin=(5, 5)),
        )
        sampling_label = toga.Label(
            "Sampling interval:",
            style=Pack(margin=(0, 0, 5, 10), text_align=LEFT, width=120),
        )
        self.textInput_sampling = toga.TextInput(
            placeholder="Enter the sampling interval for MRT",
            style=Pack(width=110, margin=(5, 5)),
        )
        time_box.add(time_label)
        time_box.add(self.textInput_time_step)
        time_box.add(sampling_label)
        time_box.add(self.textInput_sampling)
        main_box.add(time_box)

        tolerance_box = toga.Box(style=box_style)
        tolerance_label = toga.Label("Tolerance frames:", style=label_style)
        self.textInput_tolerance_frames = toga.TextInput(
            placeholder="Enter the tolerance frames",
            style=input_style,
        )
        tolerance_box.add(tolerance_label)
        tolerance_box.add(self.textInput_tolerance_frames)
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
        self.multi_line_text.value = (
            "\nMean Residence Time calculation using the TRAJEC.xyz file"
        )
        main_box.add(self.multi_line_text)

        # Buttons at the bottom
        button_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin_top=5)
        )
        self.btn_execute = toga.Button(
            "MRT calculation", style=button_style, on_press=self.workflow
        )
        self.btn_help = toga.Button(
            "Help", style=button_style, on_press=self.open_window_help 
        )
        self.btn_close = toga.Button(
            "Close", style=button_style, on_press=self.closeTopLevel
        )
        button_box.add(self.btn_execute)
        button_box.add(self.btn_help)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        # Set the content of the main window
        self.main_window.content = main_box
        self.main_window.show()

    def open_window_help(self, widget):
        window = toga.Window(
            title="Instructions for Conducting Mean Residence Time (MRT) Analysis"
        )
        main_box = toga.Box(style=Pack(direction=COLUMN, flex=1))
        multi_line_text = toga.MultilineTextInput(
            style=Pack(font_size=11, margin=(5, 5), flex=1)
        )
        multi_line_text.value = HelpGqteaWin.help_mrt

        main_box.add(multi_line_text)
        window.content = main_box
        window.show()

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
            self.multi_line_text.value = (
                f"\n\n MEAN RESIDENCE TIME ANALYSIS FINISHED SUCCSSEFULLY!\n"
            )
        except Exception as e:
            self.progress_bar.stop()
            await self.warning_function(str(e))

    def closeTopLevel(self, widget):
        self.main_window.close()
