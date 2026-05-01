import toga, os, asyncio
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, CENTER, LEFT, RIGHT
from help import HelpGqteaWin, AtomicData


class InputSHBuilder:

    async def read_params(self, widget):
        # self.atufs = 0.02418884326505  # Atomic time unit in femtoseconds.

        async def read_input(text_input, field_name, expected_type):
            value = text_input.value.strip()
            if not value:
                await self.main_window.dialog(
                    toga.InfoDialog(
                        "Error", f"Please input a valid value for {field_name}."
                    )
                )
                return None
            try:
                if expected_type == list:
                    # Assuming atom labels are space-separated integers
                    labels = [int(label) for label in value.split()]
                    if len(labels) != 6:
                        raise ValueError("Please input exactly six cell parameters.")
                    return labels
                else:
                    return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        # Read and validate inputs
        self.start_frame = await read_input(
            self.textInput_start_frame, "start frame for extraction", int
        )
        if self.start_frame is None:
            return

        self.skip_frame = await read_input(
            self.textInput_skipped_frames,
            "number of frames to be skipped between two collected frames",
            int,
        )
        if self.skip_frame is None:
            return

        self.stop_frame = await read_input(
            self.textInput_stop_frame, "final frame to be collected", int
        )
        if self.stop_frame is None:
            return
        self.num_states = await read_input(
            self.textInput_num_states, "number of states used in the spectrum", int
        )
        if self.num_states is None:
            return

        self.start_state_SHTDFT = await read_input(
            self.textInput_start_state_SHTDFT, "start state for SHTDFT", int
        )
        if self.start_state_SHTDFT is None:
            return

        self.prefix = await read_input(self.textInput_prefix, "prefix file names", str)
        if self.prefix is None:
            return

        self.dt = await read_input(
            self.textInput_dt, "Time step used in the MD run", float
        )
        if self.dt is None:
            return

        self.charge = await read_input(
            self.textInput_charge, "total charge on the system", float
        )
        if self.charge is None:
            return

        self.cell_parm = await read_input(
            self.textInput_cell_parm, "a b c cosα cosβ cosγ spaced by white space", list
        )
        if self.cell_parm is None:
            return

        self.cutoff = await read_input(
            self.textInput_cutoff, "energy cutoff in Ry", float
        )
        if self.cutoff is None:
            return

        self.dual = await read_input(self.textInput_dual, "Dual for ρ expansion", float)
        if self.dual is None:
            return

        self.cell_symmetry = await read_input(
            self.textInput_cell_symmetry, "Cell symmetry: enter 1 for cubic symm", int
        )
        if self.cell_symmetry is None:
            return
        self.max_num_steps = await read_input(
            self.textInput_max_num_steps, "max number of steps for SHTDDFT", int
        )
        if self.max_num_steps is None:
            return
        await self.main_window.dialog(
            toga.InfoDialog("Info", "All SH parameters has been read successfully!")
        )

    async def open_geometry(self, widget):

        try:
            self.geom = await self.main_window.dialog(toga.OpenFileDialog("Open file"))

            if self.geom is not None:
                self.textInput_geom.value = f"{self.geom}"
            else:
                await self.main_window.dialog(
                    toga.InfoDialog("Warning", "No file was selected!")
                )
                return

        except ValueError:
            await self.main_window.dialog(
                toga.InfoDialog("Error", "Open file was canceled!")
            )
            return

        self.output_dir = os.path.dirname(self.geom)

        # Calling function to read GEOMETRY.xyz

        self.read_geometry()

    async def open_trajectory(self, widget):
        try:
            self.trajec = await self.main_window.dialog(
                toga.OpenFileDialog("Open TRAJECTORY file")
            )

            if self.trajec is not None:
                self.textInput_trajectory.value = f"{self.trajec}"
            else:
                await self.main_window.dialog(
                    toga.InfoDialog("Warning", "No file was selected!")
                )
                return

        except ValueError:
            await self.main_window.dialog(
                toga.InfoDialog("Error", "Open file was canceled!")
            )
            return

        # Read each line to validate the trajectory file
        await self.validate_trajectory()

        await self.frames_counter()

    def read_geometry(self):
        atoms = []
        atom_counts = {}
        geom = []

        with open(self.geom, "r") as f:
            lines = f.readlines()
            self.num_atoms = int(lines[0])

        for line in lines[2:]:  # Skip the first two lines
            geom_data = line.split()
            geom.append(geom_data)

        for atom_data in geom:
            atom = atom_data[0]
            atoms.append(atom)
            if atom in atom_counts:
                atom_counts[atom] += 1
            else:
                atom_counts[atom] = 1

        self.atoms = atoms
        self.atom_counts = atom_counts

        # Open new window to atomic system info
        window = toga.Window(title=" ")
        main_box = toga.Box(style=Pack(direction=COLUMN, flex=1))
        multi_line_text = toga.MultilineTextInput(
            style=Pack(font_size=11, padding=(5, 5), flex=1)
        )
        multi_line_text.value = f"Total number of atoms: \t\t {self.num_atoms:>4d}\n"

        for key, value in atom_counts.items():
            multi_line_text.value += f"Number of {key} atoms: \t\t {value:>4d} \n"

        multi_line_text.value += f" \n"

        if len(self.info_text) != 0:
            multi_line_text.value += self.info_text

        for line in lines:
            multi_line_text.value += f"{line}"

        main_box.add(multi_line_text)
        window.content = main_box
        window.show()

    async def validate_trajectory(self):
        """Counts the number of frames in TRAJECTORY file"""

        try:
            with open(self.trajec, "r") as f:
                num_lines = 0
                while True:
                    line = f.readline().split()
                    if line == "\n" or not line:
                        break
                    # Check if the first element of the line is numeric
                    if not line[0].isdigit():
                        await self.main_window.dialog(
                            toga.InfoDialog(
                                "Error",
                                f"Error while reading line {num_lines}: {line[0]}",
                            )
                        )
                        break

                    if num_lines % 100000 == 0:
                        self.progress_label.text = f"{num_lines}"
                        await asyncio.sleep(0)  # Yield control back to the event loop

                    num_lines += 1
                self.progress_label.text = f"Number of lines: {num_lines}"
                self.total_num_frames = num_lines / self.num_atoms

        except Exception as e:
            await self.main_window.dialog(
                toga.InfoDialog("Error", f"An unexpected error occurred: {str(e)}")
            )

    async def frames_counter(self):
        try:
            with open(self.trajec, "r") as f:
                self.progress_bar.start()
                num_lines = 0
                num_frames = 0

                while True:
                    atom = []
                    line = f.readline().split()
                    num_lines += 1
                    if line == "\n" or not line:
                        break
                    atom.append(line)
                    for n in range(self.num_atoms - 1):
                        line = f.readline().split()
                        num_lines += 1
                        if line == "\n" or not line:
                            break
                        atom.append(line)
                    num_frames += 1
                    # Check if the first element of all atoms is the same
                    if not all([x[0] == atom[0][0] for x in atom]):
                        await self.main_window.dialog(
                            toga.InfoDialog(
                                "Error",
                                f"Reading problem on frame {num_frames} at line {num_lines}",
                            )
                        )
                        break

                    if num_frames % 1000 == 0:
                        progress_bar_update = (num_frames / self.total_num_frames) * 100
                        self.progress_bar.value = progress_bar_update

                    if num_frames >= (self.total_num_frames - 1000):
                        self.progress_bar.value = 100

                self.progress_label.text = f"Number of frames: {num_frames}"
                # await asyncio.sleep(0)

                await self.main_window.dialog(
                    toga.InfoDialog(
                        "Info",
                        f"Total number of lines:\t {num_lines}\n"
                        f"Total number of frames:\t {num_frames}\n",
                    )
                )
            self.progress_bar.stop()

        except Exception as e:
            await self.main_window.dialog(
                toga.InfoDialog("Error", f"An unexpected error occurred: {str(e)}")
            )

    async def extract_frames(self):

        self.progress_bar.start()

        start_frame = self.start_frame
        skip_frames = self.skip_frame
        stop_frame = self.stop_frame

        save_dir = open(f"{self.output_dir}/selected_frames.dat", "w")

        traj = open(self.trajec, "r")

        counter = 1
        step = start_frame
        self.num_select_frames = 0
        while counter < start_frame:
            for n in range(self.num_atoms):
                line = traj.readline()
            counter += 1
            if counter == start_frame:
                break

        while counter <= stop_frame:
            while counter <= step:
                if counter == step:
                    for n in range(self.num_atoms):
                        line = traj.readline()
                        save_dir.write(line)
                    self.num_select_frames += 1
                    self.progress_bar.value = round(
                        (
                            self.num_select_frames
                            / ((self.stop_frame - self.start_frame) / self.skip_frame)
                        )
                        * 100
                    )
                else:
                    for n in range(self.num_atoms):
                        line = traj.readline()
                counter += 1
            step += skip_frames

        self.progress_bar.stop()

        save_dir.close()
        traj.close()
        await self.main_window.dialog(
            toga.InfoDialog(
                "Info:",
                f"It was selected {self.num_select_frames:>4d} frames\n"
                f"The selected frames were stored into selected_frames.dat file\n",
            )
        )

    async def spectraSHInputs(self):
        """
        Generate all inputs necessary for SH dynamics: GEOMETRY,
        SPECTRA, and SHTDDFT. The GEOMETRY.xyz file is generated for
         visualization purpose
        """
        tmpList = []
        Bhor_unit = 0.52917720859  # Bohr unit conversion

        # ---New Modification---05-11-2023------------
        # Creating the directories to save the input files
        for n in range(self.num_select_frames):
            directory_path = f"{self.output_dir}/{self.prefix}-{n+1:0>4d}"
            try:
                # Create the directory
                os.makedirs(
                    directory_path, exist_ok=True
                )  # exist_ok=True will not raise an error if the directory already exists

            except OSError as error:
                await self.main_window.dialog(
                    toga.InfoDialog(
                        f"Creation of the directory {directory_path} failed due to: {error}"
                    )
                )

        await self.main_window.dialog(
            toga.InfoDialog(
                "Info",
                f"{self.num_select_frames} new folders were created successfully ",
            )
        )

        # ----End modification--------------------------

        tmp = open(f"{self.output_dir}/selected_frames.dat", "r")
        for n in range(self.num_select_frames):

            geom = open(
                f"{self.output_dir}/{self.prefix}-{n+1:0>4d}/{self.prefix}-GEOMETRY-{n+1:0>4d}",
                "w",
            )
            spectra = open(
                f"{self.output_dir}/{self.prefix}-{n+1:0>4d}/{self.prefix}-SPECTRA-{n+1:0>4d}.inp",
                "w",
            )
            tddft = open(
                f"{self.output_dir}/{self.prefix}-{n+1:0>4d}/{self.prefix}-SHTDDFT-{n+1:0>4d}.inp",
                "w",
            )
            geom_xyz = open(
                f"{self.output_dir}/{self.prefix}-{n+1:0>4d}/{self.prefix}-GEOMETRY-{n+1:0>4d}.xyz",
                "w",
            )
            geom_xyz.write(f"{self.num_atoms:<4d}\n")
            geom_xyz.write(f"File generated by gqteaWin program\n")
            for k in range(self.num_atoms):
                line = tmp.readline().split()
                line = [float(i) for i in line]
                sybl = self.atoms[k]
                geom_xyz.write(
                    f"{sybl:<4s}{Bhor_unit*line[1]:10.4f}{Bhor_unit*line[2]:10.4f}{Bhor_unit*line[3]:10.4f}\n"
                )
                geom.write(
                    f"{line[1]:20.14f}{line[2]:20.14f}{line[3]:20.14f}{line[4]:20.14f}{line[5]:20.14f}{line[6]:20.14f}\n"
                )
                tmpList.append(
                    [
                        sybl,
                        Bhor_unit * line[1],
                        Bhor_unit * line[2],
                        Bhor_unit * line[3],
                    ]
                )

            spectra.write("&INFO\n")
            tddft.write("&INFO\n")
            spectra.write(
                " CPMD INPUT FILE FOR SPECTRA CALCULATION GENERATED BY gqtea PROGRAM\n"
            )
            tddft.write(
                " CPMD INPUT FILE FOR SHTDDFT DYNAMICS GENERATED BY gqtea PROGRAM\n"
            )
            spectra.write("&END\n")
            tddft.write("&END\n")
            spectra.write("\n")
            spectra.write("&CPMD\n")
            tddft.write("&CPMD\n")
            spectra.write("\n")
            tddft.write("\n")
            spectra.write(" ELECTRONIC SPECTRA\n")
            tddft.write(" MOLECULAR DYNAMICS BO\n")
            spectra.write(" DIAGONALIZATION LANCZOS\n")
            tddft.write(" TDDFT\n")
            spectra.write(" COMPRESS WRITE32\n")
            tddft.write(" RESTART  COORDINATES VELOCITIES GEOFILE  LINRES LATEST\n")
            spectra.write("\n")
            tddft.write("\n")
            spectra.write(" MEMORY BIG\n")
            tddft.write(" MEMORY BIG\n")
            spectra.write("\n")
            tddft.write(" STORE\n")
            tddft.write("  200\n")
            tddft.write(" TIMESTEP\n")
            tddft.write(f" {self.dt}\n")
            tddft.write(" MAXSTEP\n")
            tddft.write(f" {self.max_num_steps}\n")
            tddft.write("\n")
            tddft.write(" TRAJECTORY XYZ\n")
            tddft.write(" 5\n")
            tddft.write(" NOSE IONS\n")
            tddft.write("  500.0 2000.0\n")
            tddft.write("\n")
            spectra.write("&END\n")
            tddft.write("&END\n")
            spectra.write("\n")
            tddft.write("\n")
            spectra.write("&TDDFT\n")
            tddft.write("&TDDFT\n")
            spectra.write("\n")
            spectra.write(" STATES SINGLET\n")
            tddft.write(" STATES SINGLET\n")
            spectra.write(f" {self.num_states}\n")
            tddft.write(f" {self.num_states}\n")
            spectra.write("\n")
            tddft.write(" T-SHTDDFT\n")
            tddft.write("\n")
            tddft.write(" FORCE STATE\n")
            tddft.write(f" {self.start_state_SHTDFT}\n")
            spectra.write(" TAMM-DANCOFF\n")
            tddft.write(" TAMM-DANCOFF\n")
            spectra.write(" DAVIDSON PARAMETER\n")
            tddft.write(" DAVIDSON PARAMETER\n")
            spectra.write("  150 1.D-7 50\n")
            tddft.write("  150 1.D-7 50\n")
            spectra.write("\n")
            tddft.write("\n")
            spectra.write("&END\n")
            tddft.write("&END\n")
            spectra.write("\n")
            tddft.write("\n")
            spectra.write("&DFT\n")
            tddft.write("&DFT\n")
            spectra.write(" XC_DRIVER\n")
            tddft.write(" XC_DRIVER\n")
            spectra.write(" FUNCTIONAL GGA_XC_PBE\n")
            tddft.write(" FUNCTIONAL GGA_XC_PBE\n")
            spectra.write("&END\n")
            tddft.write("&END\n")
            spectra.write("\n")
            tddft.write("\n")
            spectra.write("&SYSTEM\n")
            tddft.write("&SYSTEM\n")
            spectra.write(" CHARGE\n")
            tddft.write(" CHARGE\n")
            spectra.write(f" {self.charge}\n")
            tddft.write(f" {self.charge}\n")
            spectra.write(" SYMMETRY\n")
            tddft.write(" SYMMETRY\n")
            spectra.write(f"  {self.cell_symmetry}\n")
            tddft.write(f"  {self.cell_symmetry}\n")
            spectra.write(" ANGSTROM\n")
            tddft.write(" ANGSTROM\n")
            spectra.write(" CELL\n")
            tddft.write(" CELL\n")
            spectra.write(
                f" {self.cell_parm[0]:>7.4f}{self.cell_parm[1]/self.cell_parm[0]:>7.4f}{self.cell_parm[2]/self.cell_parm[0]:>7.4f}{self.cell_parm[3]:>7.4f}{self.cell_parm[4]:>7.4f}{self.cell_parm[5]:>7.4f}\n"
            )
            tddft.write(
                f" {self.cell_parm[0]:>7.4f}{self.cell_parm[1]/self.cell_parm[0]:>7.4f}{self.cell_parm[2]/self.cell_parm[0]:>7.4f}{self.cell_parm[3]:>7.4f}{self.cell_parm[4]:>7.4f}{self.cell_parm[5]:>7.4f}\n"
            )
            spectra.write("\n")
            tddft.write("\n")
            spectra.write(" CUTOFF\n")
            tddft.write(" CUTOFF\n")
            spectra.write(f" {self.cutoff}\n")
            tddft.write(f" {self.cutoff}\n")
            spectra.write(" DUAL\n")
            tddft.write(" DUAL\n")
            spectra.write(f" {self.dual}\n")
            tddft.write(f" {self.dual}\n")
            spectra.write("\n")
            tddft.write("\n")
            spectra.write("&END\n")
            tddft.write("&END\n")
            spectra.write("\n")
            tddft.write("\n")
            spectra.write("&ATOMS\n")
            tddft.write("&ATOMS\n")

            for key, value in self.atom_counts.items():
                spectra.write(f"*{key}_MT_PBE.psp\n")
                tddft.write(f"*{key}_MT_PBE.psp\n")

                spectra.write(f"LMAX={AtomicData.lang[key]}\n")
                tddft.write(f"LMAX={AtomicData.lang[key]}\n")

                spectra.write(f" {value}\n")
                tddft.write(f" {value}\n")

                for n in range(self.num_atoms):
                    if tmpList[n][0] == key:
                        p1 = tmpList[n][1]
                        p2 = tmpList[n][2]
                        p3 = tmpList[n][3]
                        spectra.write(f"{p1:>14.7f}{p2:>14.7f}{p3:>14.7f}\n")
                        tddft.write(f"{p1:>14.7f}{p2:>14.7f}{p3:>14.7f}\n")

            spectra.write("\n")
            tddft.write("\n")
            spectra.write("&END\n")
            tddft.write("&END\n")

            geom.close()
            geom_xyz.close()
            spectra.close()
            tddft.close()
        tmp.close()


class SurfaceHoppingUI(InputSHBuilder):
    def __init__(self, *args):
        self.info_text = ""
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        # Create the main window
        self.main_window = toga.Window(
            title="Surface Hopping Molecular Dynamics Input Files",
            size=(700, 600),
        )

        # Define common styles
        heading_style = Pack(font_size=18, font_weight="bold", padding=(0, 0, 10, 0))
        label_style = Pack(padding=(5, 5), text_align=LEFT, width=200)
        input_style = Pack(flex=1)  
        button_style = Pack(padding=5, width=150)
        box_style = Pack(direction=ROW, alignment=CENTER, padding=(0, 0, 5, 0))

        # Main container
        main_box = toga.Box(style=Pack(direction=COLUMN, padding=20))

        box_1 = toga.Box(
            style=Pack(
                direction=ROW,
                alignment=CENTER,
                padding=(0, 0, 10, 0),
            )
        )

        main_box.add(box_1)
        box_title = toga.Box(style=Pack(width=400))
        box_progress_label = toga.Box(style=Pack(flex=1))
        box_1.add(box_title)
        box_1.add(box_progress_label)

        # Title
        title_label = toga.Label("SH Molecular Dynamics Inputs", style=heading_style)
        self.progress_label = toga.Label(
            " ",
            style=Pack(
                font_size=18,
                font_weight="bold",
                text_align=RIGHT,
                padding=(0, 0, 10, 0),
                flex=1,
            ),
        )

        box_title.add(title_label)
        box_progress_label.add(self.progress_label)

        # Input fields with labels
        input_fields = [
            (
                "Starting Frame:",
                "Enter the initial frame for extraction",
                "textInput_start_frame",
            ),
            (
                "Number of skipped frames",
                "Enter the number of frames to be skippd between two extracted frames",
                "textInput_skipped_frames",
            ),
            (
                "Final frame to be extracted:",
                "Enter the final frame to be extracted",
                "textInput_stop_frame",
            ),
            (
                "Number of states:",
                "Enter the number of states in the electronic spectrum",
                "textInput_num_states",
            ),
            (
                "Initial state to start SHTDDFT:",
                "Enter the initial state for the SHTDDFT dynamics",
                "textInput_start_state_SHTDFT",
            ),
            (
                "Prefix file name:",
                "Enter prefix file name, e.g., VitC",
                "textInput_prefix",
            ),
            (
                "Molecular Dynamics Time Step:",
                "Enter the time step:, e.g., 5",
                "textInput_dt",
            ),
            (
                "Total charge on the system:",
                "Enter the charge on the system:, e.g., 0",
                "textInput_charge",
            ),
            (
                "Periodic box parameters:",
                "a b c cosα cosβ cosγ separated by white space",
                "textInput_cell_parm",
            ),
            (
                "Energy Cutoff (Ry):",
                "Enter energy cutoff (Ry), e.g., 25",
                "textInput_cutoff",
            ),
            (
                "Dual for ρ expansion:",
                "Enter the Dual to expand ρ in planewaves, e.g., 8",
                "textInput_dual",
            ),
            (
                "Cell symmetry:",
                "Enter the cell symmetry, e.g., 1",
                "textInput_cell_symmetry",
            ),
            (
                "Maximum number of steps:",
                "Enter the maximum number of steps for SHTDDFT",
                "textInput_max_num_steps",
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

        # File selection buttons
        box_geom = toga.Box(style=box_style)
        label_geom = toga.Label("Select GEOMETRY.xyz input file:", style=label_style)
        self.textInput_geom = toga.TextInput(
            placeholder="Click Browse to select GEOMETRY.xyz input file",
            style=Pack(flex=1, font_size=10, color="blue"),
        )

        box_trajectory = toga.Box(style=box_style)
        label_trajectory = toga.Label(
            "Select TRAJECTORY input file:", style=label_style
        )
        self.textInput_trajectory = toga.TextInput(
            placeholder="Click Browse to select TRAJECTORY input file",
            style=Pack(flex=1, font_size=10, color="blue"),
        )

        browse_btn_geom = toga.Button(
            "Browse", on_press=self.open_geometry, style=button_style
        )

        browse_btn_trajectory = toga.Button(
            "Browse", on_press=self.open_trajectory, style=button_style
        )

        box_geom.add(label_geom)
        box_trajectory.add(label_trajectory)
        box_geom.add(self.textInput_geom)
        box_trajectory.add(self.textInput_trajectory)
        box_geom.add(browse_btn_geom)
        box_trajectory.add(browse_btn_trajectory)
        main_box.add(box_geom)
        main_box.add(box_trajectory)

        # Progress bar creation
        box_progress_bar = toga.Box(
            style=Pack(direction=COLUMN, height=5, padding=(0, 5, 10, 5))
        )
        self.progress_bar = toga.ProgressBar(max=100)
        box_progress_bar.add(self.progress_bar)
        main_box.add(box_progress_bar)

        # Buttons at the bottom
        button_box = toga.Box(
            style=Pack(direction=ROW, alignment=CENTER, padding_top=10)
        )

        self.btn_read_params= toga.Button(
            "Read SH parms", style=button_style, on_press=self.read_params
        )

        self.btn_execute = toga.Button(
            "SH Input Builder", style=button_style, on_press=self.workflow
        )

        self.btn_help = toga.Button(
            "Help", style=button_style, on_press=self.open_window_help
        )

        self.btn_close = toga.Button(
            "Close", style=button_style, on_press=self.closeTopLevel
        )
        button_box.add(self.btn_read_params)
        button_box.add(self.btn_execute)
        button_box.add(self.btn_help)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        # Set the content of the main window
        self.main_window.content = main_box
        self.main_window.show()

    def open_window_help(self, widget):

        window = toga.Window(title=f"Instructions to use surface hopping builder ")
        main_box = toga.Box(style=Pack(direction=COLUMN, flex=1))
        multi_line_text = toga.MultilineTextInput(
            style=Pack(font_size=11, padding=(5, 5), flex=1)
        )
        multi_line_text.value = HelpGqteaWin.help_spectraSH

        main_box.add(multi_line_text)

        window.content = main_box

        window.show()

    async def workflow(self, widget):
        await self.extract_frames()
        await self.spectraSHInputs()

    def closeTopLevel(self, widget):
        self.main_window.close()
