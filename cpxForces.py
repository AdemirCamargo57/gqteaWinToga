import os
import numpy as np
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, CENTER

class ForceXYZtoForce:
    """
    Converts cp.x trajectory files in XYZ format to a force file (in atomic units).
    Processes file selection, extracts metadata (atom count and labels) from the cp.x input file,
    and computes the force magnitude for each atom in the trajectory.
    """
    def __init__(self):
        # File paths and conversion parameters
        self.input_file = None            # cp.x input file path
        self.cp_force_file = None         # cp.x trajectory (*.for) file path in XYZ format
        self.atom_labels = []             # List of atomic labels extracted from the input file
        self.nat = 0                      # Number of atoms
        self.bohr_to_angstrom = 0.5291772109  # (Unused conversion factor)

    async def select_file(self, widget, title, target_input):
        """
        Opens a file dialog with a specified title and updates the target TextInput.

        Parameters:
            widget: The widget that triggers the dialog.
            title (str): Title of the file dialog.
            target_input (toga.TextInput): Widget to display the selected file path.

        Returns:
            str or None: The selected file path, or None if no file was chosen.
        """
        try:
            file_path = await self.main_window.open_file_dialog(title=title)
            if file_path:
                target_input.value = file_path
                return file_path
            else:
                await self.main_window.info_dialog("Warning", "No file was selected!")
                return None
        except Exception as e:
            await self.main_window.error_dialog("Error", f"Failed to open file: {e}")
            return None

    async def select_cp_input_file(self, widget):
        """Opens a dialog to select the cp.x input file."""
        self.input_file = await self.select_file(widget, "Select cp.x Input File", self.input_textInput)

    async def select_cp_force_xyz_file(self, widget):
        """Opens a dialog to select the cp.x trajectory (*.for) file in XYZ format."""
        self.cp_force_file = await self.select_file(widget, "Select cp.x *.for File (XYZ)", self.cp_force_textInput)

    def read_cp_input(self):
        """
        Reads the cp.x input file to determine the number of atoms and extract atomic labels.
        Expects a line containing 'nat' and an 'ATOMIC_POSITIONS' section.

        Raises:
            FileNotFoundError: If the input file does not exist.
            ValueError: If parsing the atom count fails.
        """
        cp_input_path = self.input_textInput.value
        if not cp_input_path or not os.path.exists(cp_input_path):
            raise FileNotFoundError(f"cp.x input file not found: {cp_input_path}")
        
        self.atom_labels.clear()
        with open(cp_input_path, 'r') as f:
            for line in f:
                tokens = line.strip().split()
                if not tokens:
                    continue

                if 'nat' in tokens:
                    try:
                        # Assumes the number of atoms is the third token (stripping any trailing comma)
                        self.nat = int(tokens[2].strip(','))
                    except (IndexError, ValueError) as err:
                        raise ValueError(f"Error parsing number of atoms from: {line}") from err

                if 'ATOMIC_POSITIONS' in tokens:
                    # Read exactly self.nat lines following the header
                    for _ in range(self.nat):
                        pos_line = f.readline().strip().split()
                        if pos_line:
                            self.atom_labels.append(pos_line[0])
                    break

    def write_force(self):
        """
        Processes the cp.x trajectory file and writes an output file with computed force magnitudes.
        Each frame is written in XYZ format, where the force magnitude is computed as sqrt(x^2 + y^2 + z^2).

        Output:
            The converted file is saved as 'cp_force.for' in the same directory as the trajectory file.

        Raises:
            FileNotFoundError: If the trajectory file does not exist.
            ValueError: If any frame contains incomplete coordinate data.
        """
        cp_force_path = self.cp_force_textInput.value
        if not cp_force_path or not os.path.exists(cp_force_path):
            raise FileNotFoundError(f"cp.x trajectory file not found: {cp_force_path}")

        output_dir = os.path.dirname(cp_force_path)
        output_file = os.path.join(output_dir, 'cp_force.for')

        with open(cp_force_path, 'r') as f_in, open(output_file, 'w') as f_out:
            while True:
                title_line = f_in.readline().strip()
                if not title_line:
                    break  # End of file

                # Write frame header: number of atoms and a comment
                f_out.write(f"{self.nat}\n")
                f_out.write("Trajectory generated by cp2xyz converter\n")

                # Process each atom's coordinate line in the current frame
                for i in range(self.nat):
                    line = f_in.readline().strip().split()
                    if len(line) < 3:
                        raise ValueError("Incomplete atomic coordinate data in trajectory.")
                    try:
                        x, y, z = map(float, line[:3])
                    except ValueError as err:
                        raise ValueError("Non-numeric coordinate encountered.") from err

                    # Compute the force magnitude (Euclidean norm)
                    force = np.hypot(np.hypot(x, y), z)
                    f_out.write(f"{self.atom_labels[i]}  {force:14.7f}\n")

    def perform_conversion(self, widget):
        """
        Executes the conversion process:
         1. Reads the cp.x input file to extract the number of atoms and labels.
         2. Processes the trajectory file to compute force magnitudes.
        On success, updates the output file label; on failure, displays an error dialog.
        """
        try:
            self.read_cp_input()
            self.write_force()
            output_dir = os.path.dirname(self.cp_force_textInput.value)
            output_file = os.path.join(output_dir, 'cp_force.for')
            self.output_file_label.text = f"Output File: {output_file}"
        except Exception as e:
            self.main_window.error_dialog("Conversion Error", str(e))


class CPforcesUI(ForceXYZtoForce):
    """
    GUI application for converting cp.x trajectory files (XYZ format) to atomic force units.
    Inherits conversion functionality from ForceXYZtoForce.
    """
    def __init__(self, app):
        super().__init__()
        self.app = app
        self.main_window = toga.Window(
            title="Convert cp.x force file (*.for) in trajectory force file in a.u.", size=(700, 300)
        )
        self.layout_main_window()
        self.main_window.show()

    def layout_main_window(self):
        """
        Constructs the main window layout with file selection inputs, an output label, and a conversion button.
        """
        main_box_style = Pack(direction=COLUMN, margin=10)
        box_style = Pack(direction=ROW, margin=5)
        label_style = Pack(width=120, margin=(0, 5))
        input_style = Pack(flex=1, margin=(0, 5))
        button_style = Pack(margin=5)

        # cp.x Input File Selection
        input_file_box = toga.Box(style=box_style)
        input_file_label = toga.Label("cp.x Input File:", style=label_style)
        self.input_textInput = toga.TextInput(
            readonly=False,
            placeholder="Select cp.x input file",
            style=input_style
        )
        browse_button = toga.Button("Browse", on_press=self.select_cp_input_file, style=button_style)
        input_file_box.add(input_file_label)
        input_file_box.add(self.input_textInput)
        input_file_box.add(browse_button)

        # cp.x Trajectory (*.for) File Selection
        cp_for_box = toga.Box(style=box_style)
        cp_for_label = toga.Label("cp *.for file:", style=label_style)
        self.cp_force_textInput = toga.TextInput(
            readonly=False,
            placeholder="Select cp.x force trajectory (*.for) file in xyz format",
            style=input_style
        )
        cp_for_browse_button = toga.Button("Browse", on_press=self.select_cp_force_xyz_file, style=button_style)
        cp_for_box.add(cp_for_label)
        cp_for_box.add(self.cp_force_textInput)
        cp_for_box.add(cp_for_browse_button)

        # Output File Display
        self.output_file_label = toga.Label("Output File: ", style=Pack(margin=(5, 5)))

        # Convert Button
        convert_button = toga.Button(
            "Convert",
            on_press=self.perform_conversion,
            style=Pack(margin_left=300, align_items=CENTER, width=100)
        )

        # Assemble the main layout
        main_box = toga.Box(style=main_box_style)
        main_box.add(input_file_box)
        main_box.add(cp_for_box)
        main_box.add(self.output_file_label)
        main_box.add(convert_button)
        self.main_window.content = main_box
