import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from help import HelpGqteaWin, AtomicData
from framesCounter import FramesCounter


class CpmdInputBuilder(FramesCounter):

    async def warning_function(self, msg):
        await self.main_window.dialog(toga.InfoDialog("Error", f"{msg}"))

    async def read_params(self, widget):

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
                        raise ValueError("Please input exactly two atom labels.")
                    return labels
                else:
                    return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        # Read and validate inputs
        self.prefix_name = await read_input(
            self.textInput_prefix, "prefix file name", str
        )
        if self.prefix_name is None:
            return

        self.charge = await read_input(
            self.textInput_charge, "charge on the system", int
        )
        if self.charge is None:
            return

        self.cell_parm = await read_input(self.textInput_cell_parm, "sampling", list)
        if self.cell_parm is None:
            return

        self.cutoff = await read_input(
            self.textInput_cutoff, "simulation temperature", float
        )
        if self.cutoff is None:
            return

        self.dual = await read_input(self.textInput_dual, "atom labels", int)
        if self.dual is None:
            return

        symmetry_value = str(self.selection_symmetry.value).strip()
        if not symmetry_value:
            await self.main_window.dialog(
                toga.InfoDialog("Error", "Please select a valid value for cell symmetry.")
            )
            return
        self.symmetry = int(symmetry_value.split()[0])

        update_text = (
            f"Prefix name --> {self.prefix_name}\n"
            f"Charge on the systema --> {self.charge}\n"
            f"Cell parameters --> {self.cell_parm}\n"
            f"Energy cutoff --> {self.cutoff}\n"
            f"Dual --> {self.dual}\n"
            f"Cell symmetry --> {self.symmetry}\n"
        )
        self.multi_line_text.value = update_text

    def cpmd_input_builder(self):
        """Generate CPMD input files based on user input and display the contents."""

        # Parse user input for cell parameters

        lat_a, lat_b, lat_c = map(float, self.cell_parm[:3])
        cos_alfa, cos_beta, cos_gamma = map(float, self.cell_parm[3:])

        cutoff = self.cutoff
        dual = self.dual
        charge = self.charge
        symm = self.symmetry

        # Define the paths for the CPMD input files
        file_paths = {
            "wfn": f"{self.output_dir}/{self.prefix_name}_wfnopt.inp",
            "equi": f"{self.output_dir}/{self.prefix_name}_eq.inp",
            "sim": f"{self.output_dir}/{self.prefix_name}_sim.inp",
            "bomd": f"{self.output_dir}/{self.prefix_name}_bomd.inp",
        }

        # Create and write to the CPMD input files
        files = {key: open(path, "w") for key, path in file_paths.items()}
        for f in files.values():
            f.write(f"&INFO\nInput file created by gqteaWin\n \n&END\n\n")

        # Write specific sections to each file
        files["equi"].write(
            f"&CPMD\n MOLECULAR DYNAMICS CP\n RESTART WAVEFUNCTION COORDINATES LATEST\n PRINT FORCES ON\n MEMORY BIG\n TRAJECTORY SAMPLE XYZ\n 5\n RESTFILE\n 1\n STORE\n 200\n MAXSTEP\n 1000\n TIMESTEP\n 3.0\n NOSE IONS MASSIVE\n 300.0 2000.0\n NOSE ELECTRONS\n 0.007 15000.0\n NOSE PARAMETERS\n 3 3 3 6.0D0 15 4\n SPLINE POINTS\n 5000\n&END\n\n"
        )
        files["sim"].write(
            f"&CPMD\n MOLECULAR DYNAMICS CP\n RESTART WAVEFUNCTION COORDINATES VELOCITIES NOSEP NOSEE LATEST\n MEMORY BIG\n TRAJECTORY SAMPLE XYZ\n 5\n RESTFILE\n 1\n STORE\n 200\n MAXSTEP\n 300000\n TIMESTEP\n 5.0\n NOSE IONS MASSIVE\n 300.0 2000.0\n NOSE ELECTRONS\n 0.007 15000.0\n NOSE PARAMETERS\n 3 3 3 6.0D0 15 4\n SPLINE POINTS\n 5000\n&END\n\n"
        )
        files["wfn"].write(
            f"&CPMD\n OPTIMIZE WAVEFUNCTION\n CONVERGENCE ORBITALS\n 1.0d-6\n CENTER MOLECULE ON\n PRINT FORCES ON\n&END\n\n"
        )
        files["bomd"].write(
            f"&CPMD\n MOLECULAR DYNAMICS BO\n CONVERGENCE ORBITALS\n 1.0d-5\n CENTER MOLECULE ON\n TRAJECTORY SAMPLE XYZ\n 3\n RESTFILE\n 1\n STORE\n 200\n MAXSTEP\n 300000\n TIMESTEP\n 3.0\n&END\n\n"
        )
        # Write the DFT section
        for f in files.values():
            f.write(f"&DFT\n XC_DRIVER\n FUNCTIONAL GGA_XC_PBE\n&END\n\n")
        # Write the system information
        for f in files.values():
            f.write(
                f"&SYSTEM\n CHARGE\n {charge}\n SYMMETRY\n {symm}\n ANGSTROM\n CELL\n {lat_a:>6.2f} {lat_b/lat_a:>6.2f} {lat_c/lat_a:>6.2f} {cos_alfa:>6.2f} {cos_beta:>6.2f} {cos_gamma:>6.2f}\n CUTOFF\n {cutoff}\n DUAL\n {dual}\n&END\n\n"
            )
        # Write the atoms information
        for f in files.values():
            f.write(f"&ATOMS\n")
        for key, value in self.element_count.items():
            for f in files.values():
                f.write(
                    f"*{key}_VDB_PBE.psp FORMATTED\n LMAX={AtomicData.lang[key]}\n {value}\n"
                )
                for n in range(self.num_atoms):
                    if self.coords[n][0] == key:
                        p1, p2, p3 = self.coords[n][1:]
                        f.write(f"{p1:>14.7f}{p2:>14.7f}{p3:>14.7f}\n")
                f.write(f"\n")
        for f in files.values():
            f.write(f"&END")
        
        # Close all file streams
        for f in files.values():
            f.close()
        self.print_cpmd_inputs()

    def print_cpmd_inputs(self):
        """Print the generated CPMD input files."""
        for file_type in ["wfnopt", "eq", "sim", "bomd"]:
            file_path = f"{self.output_dir}/{self.prefix_name}_{file_type}.inp"
            self.multi_line_text.value += (
                f"\n>>> CPMD input file for {file_type} <<<\n\n"
            )

            with open(file_path, "r") as f:
                lines = f.readlines()
                for line in lines:
                    self.multi_line_text.value += line


class CpmdInputUI(CpmdInputBuilder):
    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        # Create the main window
        self.main_window = toga.Window(
            title="CPMD Input Files",
            size=(700, 600),
        )

        # Define common styles
        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style = Pack(margin=(5, 5), text_align=LEFT, width=200)
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=100)
        box_style = Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0))

        # Main container
        main_box = toga.Box(style=Pack(direction=COLUMN, margin=20))

        box_1 = toga.Box(
            style=Pack(
                direction=ROW,
                align_items=CENTER,
                margin=(0, 0, 10, 0))
        )

        main_box.add(box_1)
        box_1a = toga.Box(style=Pack(width=400))
        box_1b = toga.Box(style=Pack(width=150))
        box_1c = toga.Box(style=Pack(width=100))
        box_1.add(box_1a)
        box_1.add(box_1b)
        box_1.add(box_1c)

        # Title
        title_label = toga.Label(
            "CPMD Input Files Builder", style=heading_style
        )
        empty_label = toga.Label(" ", style=heading_style)
        self.progress_label = toga.Label(" ", style=heading_style)

        box_1a.add(title_label)
        box_1b.add(empty_label)
        box_1c.add(self.progress_label)

        # Input fields with labels
        input_fields = [
            (
                "Prefix file name:",
                "Input prefix file name, e.g., VitC",
                "textInput_prefix",
            ),
            (
                "Total charge on the system:",
                "Input the charge on the system:, e.g., 0",
                "textInput_charge",
            ),
            (
                "Periodic box parameters:",
                "a b c cosα cosβ cosγ separated by white space",
                "textInput_cell_parm",
            ),
            (
                "Energy Cutoff (Ry):",
                "Input energy cutoff (Ry), e.g., 25",
                "textInput_cutoff",
            ),
            (
                "Dual for ρ expansion:",
                "Input the Dual to expand ρ in planewaves, e.g., 8",
                "textInput_dual",
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

        symmetry_box = toga.Box(style=box_style)
        symmetry_label = toga.Label("Cell symmetry:", style=label_style)
        self.selection_symmetry = toga.Selection(
            items=[
                "1 - CUBIC a=b=c α=β=γ=90°",
                "6 - TETRAGONAL a=b≠c α=β=γ=90°",
                "8 - ORTHORHOMBIC a≠b≠c α=β=γ=90°",
            ],
            style=input_style,
        )
        self.selection_symmetry.value = "1 - CUBIC a=b=c α=β=γ=90°"
        symmetry_box.add(symmetry_label)
        symmetry_box.add(self.selection_symmetry)
        main_box.add(symmetry_box)

        # File selection button
        file_box = toga.Box(style=box_style)
        file_label = toga.Label("Select input file:", style=label_style)
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select input file in xyz format",
            style=input_style,
        )
        browse_button = toga.Button(
            "Browse", on_press=self.open_geometry_xyz, style=button_style
        )
        progress_box = toga.Box(style=Pack(direction=COLUMN))
        self.progress_bar = toga.ProgressBar(max=100)
        progress_box.add(self.progress_bar)

        file_box.add(file_label)
        file_box.add(self.textInput_file)
        file_box.add(browse_button)
        main_box.add(file_box)
        main_box.add(progress_box)

        # Multi-line text input for help or output
        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(5, 0), font_size=11)
        )
        self.multi_line_text.value = HelpGqteaWin.help_cpmd_input
        main_box.add(self.multi_line_text)

        # Buttons at the bottom
        button_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin_top=5)
        )
        self.btn_execute = toga.Button(
            "Input Builder", style=button_style, on_press=self.workflow
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

    async def workflow(self, widget):
        await self.read_params(widget)
        if not hasattr(self, "trajec") or not self.trajec:
            await self.warning_function("No trajectory file selected.")
            return

        self.cpmd_input_builder()

    def closeTopLevel(self, widget):
        self.main_window.close()
