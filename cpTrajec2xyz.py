import os
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, CENTER

class TrajectoryConverter:
    """
    A class to convert cp.x trajectory files into XYZ format.
    This class handles file selection, reading metadata from the cp.x input file,
    and writing the trajectory data into a new XYZ file.
    """
    def __init__(self):
        # Initialize file paths and data storage
        self.input_file = None           # Path for the cp.x input file
        self.cp_trajec_file = None       # Path for the cp.x *.pos trajectory file
        self.atom_labels = []            # List to store atomic labels read from input file
        self.nat = 0                     # Number of atoms
        self.bohr_to_angstrom = 0.5291772109  # Conversion factor from Bohr to Angstrom

    async def select_file(self, widget, title, target_input):
        """
        Open a file dialog with the specified title and update the provided TextInput widget.

        Parameters:
            widget: The widget that triggers the file dialog.
            title (str): The title for the file dialog.
            target_input (toga.TextInput): The TextInput widget where the file path will be set.
        Returns:
            The selected file path or None if no file was selected.
        """
        try:
            file_path = await self.main_window.open_file_dialog(title=title)
            if file_path:
                target_input.value = file_path  # Update the UI with the selected file path
                return file_path
            else:
                await self.main_window.info_dialog("Warning", "No file was selected!")
                return None
        except Exception as e:
            await self.main_window.error_dialog("Error", f"Failed to open file: {e}")
            return None

    async def select_cp_input_file(self, widget):
        """
        Open a file dialog to select the cp.x input file.
        """
        self.input_file = await self.select_file(
            widget, "Select cp.x Input File", self.input_textInput
        )

    async def select_cp_trajec_file(self, widget):
        """
        Open a file dialog to select the cp.x *.pos trajectory file.
        """
        self.cp_trajec_file = await self.select_file(
            widget, "Select cp.x *.pos Trajectory", self.cp_trajec_textInput
        )

    def read_cp_input(self):
        """
        Read the cp.x input file to extract the number of atoms (nat)
        and the atomic labels from the ATOMIC_POSITIONS section.
        """
        cp_input_path = self.input_textInput.value
        if not cp_input_path or not os.path.exists(cp_input_path):
            raise FileNotFoundError(f"cp.x input file not found: {cp_input_path}")
        
        # Clear any previously stored atomic labels
        self.atom_labels.clear()

        with open(cp_input_path, 'r') as f:
            for line in f:
                tokens = line.strip().split()
                if not tokens:
                    continue  # Skip empty lines

                # Extract number of atoms from a line containing 'nat'
                if 'nat' in tokens:
                    try:
                        self.nat = int(tokens[2].strip(','))  # Extract the number of atoms
                    except (IndexError, ValueError):
                        raise ValueError("Failed to parse number of atoms from line: " + line)

                # When the ATOMIC_POSITIONS section is found, read the next nat lines
                if 'ATOMIC_POSITIONS' in tokens:
                    for _ in range(self.nat):
                        pos_line = f.readline().strip().split()
                        if pos_line:
                            # Store the atomic label (first token) from each line
                            self.atom_labels.append(pos_line[0])
                    # Once the atomic positions are processed, exit the loop
                    break

    def write_xyz_trajec(self):
        """
        Convert the cp.x *.pos trajectory file to an XYZ formatted file.
        The output file is saved as 'cp_trajec.xyz' in the same directory as the input trajectory.
        """
        cp_trajec_path = self.cp_trajec_textInput.value
        if not cp_trajec_path or not os.path.exists(cp_trajec_path):
            raise FileNotFoundError(f"cp.x trajectory file not found: {cp_trajec_path}")

        # Define the output file path
        output_dir = os.path.dirname(cp_trajec_path)
        output_file = os.path.join(output_dir, 'cp_trajec.xyz')

        with open(cp_trajec_path, 'r') as f_in, open(output_file, 'w') as f_out:
            # Process each frame in the trajectory file
            while True:
                # Read the title line (which signals the start of a new frame)
                title_line = f_in.readline().strip()
                if not title_line:
                    break  # End of file reached

                # Write the number of atoms and a comment line in XYZ format
                f_out.write(f"{self.nat}\n")
                f_out.write("Trajectory generated by cp2xyz converter\n")

                # Write each atomic position line in the frame
                for i in range(self.nat):
                    atom_line = f_in.readline().strip().split()
                    if len(atom_line) < 3:
                        raise ValueError("Incomplete atomic coordinate data in trajectory.")
                    x, y, z = atom_line[:3]
                    # Convert atomic positions from Bohr to Angstrom
                    x = float(x) * self.bohr_to_angstrom
                    y = float(y) * self.bohr_to_angstrom
                    z = float(z) * self.bohr_to_angstrom

                    # Write the atomic label and its coordinates
                    f_out.write(f"{self.atom_labels[i]}  {x:>14.7f}  {y:>14.7f}  {z:>14.7}\n")

    def perform_conversion(self, widget):
        """
        Convert the cp.x *.pos trajectory file to XYZ format by:
            1. Reading the cp.x input file for metadata (number of atoms and atomic labels).
            2. Converting the trajectory data to XYZ format.
        Updates the output file label upon successful conversion.
        """
        try:
            self.read_cp_input()
            self.write_xyz_trajec()
            # Determine the output file location and update the label for user feedback
            cp_trajec_path = self.cp_trajec_textInput.value
            output_dir = os.path.dirname(cp_trajec_path)
            output_file = os.path.join(output_dir, 'cp_trajec.xyz')
            self.output_file_label.text = f"Output File: {output_file}"
        except Exception as e:
            # Display an error dialog if conversion fails
            self.main_window.error_dialog("Conversion Error", str(e))

class CPtraj2xyzUI(TrajectoryConverter):
    """
    The GUI application that provides an interface for converting cp.x trajectory files to XYZ format.
    Inherits conversion functionality from TrajectoryConverter.
    """
    def __init__(self, app):
        super().__init__()
        self.app = app
        # Create the main window with a title and a specified size
        self.main_window = toga.Window(title="Converting cp.x trajectory files to XYZ format", size=(700, 300))
        self.layout_main_window()
        self.main_window.show()

    def layout_main_window(self):
        """
        Layout the main window with file input fields, labels, and a conversion button.
        """
        # Define style properties for the layout using Toga's Pack layout
        main_box_style = Pack(direction=COLUMN, margin=10)
        box_style = Pack(direction=ROW, margin=5)
        label_style = Pack(width=120, margin=(0, 5))
        input_style = Pack(flex=1, margin=(0, 5))
        button_style = Pack(margin=5)

        # --- cp.x Input File Selection ---
        input_file_box = toga.Box(style=box_style)
        input_file_label = toga.Label("cp.x Input File:", style=label_style)
        self.input_textInput = toga.TextInput(
            readonly=False, 
            placeholder="Select cp.x input file", 
            style=input_style
        )
        browse_button = toga.Button(
            "Browse",
            on_press=self.select_cp_input_file,
            style=button_style
        )
        input_file_box.add(input_file_label)
        input_file_box.add(self.input_textInput)
        input_file_box.add(browse_button)

        # --- cp.x *.pos Trajectory File Selection ---
        cp_pos_trajec_box = toga.Box(style=box_style)
        cp_pos_trajec_label = toga.Label("Select *.pos file", style=label_style)
        self.cp_trajec_textInput = toga.TextInput(
            readonly=False, 
            placeholder="Select *.pos trajectory file from cp.x run", 
            style=input_style
        )
        cp_pos_trajec_browse_button = toga.Button(
            "Browse",
            on_press=self.select_cp_trajec_file,
            style=button_style
        )
        cp_pos_trajec_box.add(cp_pos_trajec_label)
        cp_pos_trajec_box.add(self.cp_trajec_textInput)
        cp_pos_trajec_box.add(cp_pos_trajec_browse_button)

        # --- Output File Display ---
        self.output_file_label = toga.Label(
            "Output File: ",
            style=Pack(margin=(5, 5))
        )

        # --- Convert Button ---
        convert_button = toga.Button(
            "Convert",
            on_press=self.perform_conversion,
            style=Pack(margin_left=300, align_items=CENTER, width=100)
        )

        # --- Main Layout Box ---
        main_box = toga.Box(style=main_box_style)
        main_box.add(input_file_box)
        main_box.add(cp_pos_trajec_box)
        main_box.add(self.output_file_label)
        main_box.add(convert_button)

        # Set the main window content to the assembled layout
        self.main_window.content = main_box

