import os
import asyncio
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, CENTER

class XYZ2BohrConverter:
    """
    A converter for XYZ files with coordinates in Ångströms to Bohr units.
    Conversion is performed by multiplying coordinates by a conversion factor.
    """
    ANGSTROM_TO_BOHR = 1.889725989  # 1 Ångström = 1.889725989 Bohr

    def __init__(self):
        self.input_file = None   # Path to the input XYZ file
        self.output_file = None  # Path to the output XYZ file

    async def convert_xyz_angstrom_to_bohr(self, input_filename, output_filename):
        """
        Convert coordinates in the XYZ file from Ångströms to Bohr units.

        Parameters:
            input_filename (str): Path of the input file.
            output_filename (str): Path to save the converted file.

        Returns:
            bool: True if conversion succeeded, False otherwise.
        """
        try:
            with open(input_filename, 'r') as infile:
                lines = infile.readlines()
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Error reading input file: {e}"))
            return False

        if len(lines) < 2:
            await self.main_window.dialog(toga.ErrorDialog("Error", "The input file does not have the expected XYZ format."))
            return False

        # The first two lines are the atom count and a comment.
        atom_count_line = lines[0]
        comment_line = lines[1]
        converted_lines = [atom_count_line, comment_line]

        # Process each atomic coordinate line.
        for line in lines[2:]:
            parts = line.strip().split()
            if len(parts) == 4:
                converted_lines.append(line)
                continue
            # Check if the line has the expected format.
            elif len(parts) != 4:
                # If the line does not have 4 parts, it may be a comment or an empty line.
                await self.main_window.dialog(toga.ErrorDialog("Error", f"Unexpected line format: {line.strip()}"))
                return False

            atom_label = parts[0]
            try:
                # Convert coordinate strings to floats.
                x_ang, y_ang, z_ang = map(float, parts[1:4])
                # Convert the coordinates from Ångströms to Bohr.
                x_bohr = x_ang * self.ANGSTROM_TO_BOHR
                y_bohr = y_ang * self.ANGSTROM_TO_BOHR
                z_bohr = z_ang * self.ANGSTROM_TO_BOHR
                converted_lines.append(f"{atom_label} {x_bohr:.6f} {y_bohr:.6f} {z_bohr:.6f}\n")
            except ValueError:
                # If conversion fails, output the line unchanged.
                converted_lines.append(line)

        try:
            with open(output_filename, 'w') as outfile:
                outfile.writelines(converted_lines)
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Error writing output file: {e}"))
            return False

        return True

class XYZ2BohrConverterUI(XYZ2BohrConverter):
    """
    Toga-based UI for converting XYZ files from Ångströms to Bohr units.
    Provides a graphical interface for file selection and conversion.
    """
    def __init__(self, *args):
        super().__init__()
        self.main_window = toga.Window(
            title="XYZ to Bohr Converter", size=(700, 300)
        )
        self.build_ui()
        self.main_window.show()

    async def select_input_file(self, widget, target_input, dialog_title="Select Input XYZ File"):
        """
        Opens a file dialog to select a file and sets the chosen path in target_input.
        """
        try:
            file_path = await self.main_window.dialog(toga.OpenFileDialog(title=dialog_title))
            if file_path:
                target_input.value = file_path
                self.file_path = file_path  # Save file path for determining output directory.
                return file_path
            else:
                await self.main_window.dialog(toga.InfoDialog("Warning", "No file was selected!"))
                return None
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to open file: {e}"))
            return None

    async def on_convert_press(self, widget):
        """
        Async handler for the 'Convert' button click.
        Retrieves file paths from text inputs, validates them, and performs the conversion.
        """
        input_filename = self.input_text_input.value.strip()
        output_filename = self.output_text_input.value.strip()

        if not input_filename:
            await self.main_window.dialog(toga.ErrorDialog("Input Error", "Please select an input file."))
            return

        if not output_filename:
            await self.main_window.dialog(toga.ErrorDialog("Input Error", "Please specify an output file name."))
            return

        # Determine the output directory (falling back to current directory if no file_path exists).
        output_dir = os.path.dirname(self.file_path) if hasattr(self, 'file_path') and self.file_path else os.getcwd()
        output_filepath = os.path.join(output_dir, output_filename)

        success = await self.convert_xyz_angstrom_to_bohr(input_filename, output_filepath)
        if success:
            self.output_file_label.text = f"Conversion successful. Output saved as: {output_filepath}"
        else:
            self.output_file_label.text = "Conversion failed. Check the console for details."

    def build_ui(self):
        """
        Build and set up the main window layout with file selection fields, labels, and a conversion button.
        """
        main_box_style = Pack(direction=COLUMN, padding=10)
        row_style = Pack(direction=ROW, padding=5)
        label_style = Pack(width=150, padding=(0, 5))
        input_style = Pack(flex=1, padding=(0, 5))
        button_style = Pack(padding=5)

        # --- Input File Selection ---
        input_box = toga.Box(style=row_style)
        input_label = toga.Label("Input XYZ file:", style=label_style)
        self.input_text_input = toga.TextInput(placeholder="Select input file", style=input_style)
        input_browse_button = toga.Button(
            "Browse",
            on_press=lambda w: asyncio.create_task(self.select_input_file(w, self.input_text_input)),
            style=button_style
        )
        input_box.add(input_label)
        input_box.add(self.input_text_input)
        input_box.add(input_browse_button)

        # --- Output File Name ---
        output_box = toga.Box(style=row_style)
        output_label = toga.Label("Output file name:", style=label_style)
        self.output_text_input = toga.TextInput(placeholder="Enter output file name", style=input_style)
        output_box.add(output_label)
        output_box.add(self.output_text_input)

        # --- Status Label for Output ---
        self.output_file_label = toga.Label("Output File: ", style=Pack(padding=5))

        # --- Convert Button ---
        convert_button = toga.Button(
            "Convert",
            on_press=lambda w: asyncio.create_task(self.on_convert_press(w)),
            style=Pack(padding_left=300, alignment=CENTER, width=100)
        )

        # --- Assemble the UI layout ---
        main_box = toga.Box(style=main_box_style)
        main_box.add(input_box)
        main_box.add(output_box)
        main_box.add(self.output_file_label)
        main_box.add(convert_button)

        self.main_window.content = main_box

