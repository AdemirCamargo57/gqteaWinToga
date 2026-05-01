import os
import re
import asyncio
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, CENTER


class CpmdInputToXYZConverter:
    """
    Converter for CPMD input files to XYZ format.

    Features
    --------
    1. Reads CELL information from the &SYSTEM block.
    2. Reads atomic coordinates from the &ATOMS block.
    3. Detects whether coordinates are in ANGSTROM or BOHR.
    4. Writes a standard XYZ file.
    5. Writes the box lattice information in the XYZ comment line.
    """

    BOHR_TO_ANGSTROM = 0.529177210903

    def __init__(self):
        self.input_file = None
        self.output_file = None
        self.lines = []
        self.cell_values = []
        self.cell_comment = ""
        self.use_angstrom = False
        self.atoms = []

    async def convert_cpmd_input_to_xyz(self, input_filename, output_filename):
        """
        Convert a CPMD input file into an XYZ file.

        Parameters
        ----------
        input_filename : str
            Path to the CPMD input file.
        output_filename : str
            Path to save the XYZ file.

        Returns
        -------
        bool
            True if conversion succeeded, False otherwise.
        """
        try:
            with open(input_filename, "r", encoding="utf-8") as infile:
                self.lines = infile.readlines()
        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Error reading input file: {e}")
            )
            return False

        if not self.lines:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", "The input file is empty.")
            )
            return False

        try:
            self.use_angstrom = self._detect_angstrom_mode()
            self.cell_values = self._extract_cell_values()
            self.cell_comment = self._build_cell_comment(self.cell_values)
            self.atoms = self._extract_atoms()
        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Error processing CPMD input: {e}")
            )
            return False

        if not self.atoms:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", "No atoms were found in the &ATOMS block.")
            )
            return False

        try:
            with open(output_filename, "w", encoding="utf-8") as outfile:
                outfile.write(f"{len(self.atoms)}\n")
                outfile.write(f"{self.cell_comment}\n")
                for symbol, x, y, z in self.atoms:
                    outfile.write(f"{symbol:<3s} {x:16.8f} {y:16.8f} {z:16.8f}\n")
        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Error writing output file: {e}")
            )
            return False

        return True

    def _detect_angstrom_mode(self):
        """
        Detect if ANGSTROM is present in the &SYSTEM block.
        """
        inside_system = False

        for raw_line in self.lines:
            line = raw_line.strip().upper()

            if line == "&SYSTEM":
                inside_system = True
                continue

            if inside_system and line == "&END":
                break

            if inside_system and line == "ANGSTROM":
                return True

        return False

    def _extract_cell_values(self):
        """
        Extract CELL values from the &SYSTEM block.
        """
        inside_system = False

        for i, raw_line in enumerate(self.lines):
            line = raw_line.strip().upper()

            if line == "&SYSTEM":
                inside_system = True
                continue

            if inside_system and line == "&END":
                break

            if inside_system and line == "CELL":
                if i + 1 >= len(self.lines):
                    raise ValueError("CELL keyword found, but no values were provided.")

                next_line = self.lines[i + 1].strip()
                values = self._extract_floats_from_line(next_line)

                if not values:
                    raise ValueError("CELL keyword found, but the following line has no numeric values.")

                return values

        raise ValueError("No CELL section was found in the &SYSTEM block.")

    def _extract_atoms(self):
        """
        Extract atom coordinates from the &ATOMS block.
        """
        inside_atoms = False
        parsed_atoms = []

        i = 0
        nlines = len(self.lines)

        while i < nlines:
            stripped = self.lines[i].strip()

            if stripped.upper() == "&ATOMS":
                inside_atoms = True
                i += 1
                continue

            if inside_atoms and stripped.upper() == "&END":
                break

            if inside_atoms and stripped.startswith("*"):
                symbol = self._extract_element_symbol_from_pseudopotential_line(stripped)

                if i + 2 >= nlines:
                    raise ValueError("Unexpected end of file while reading atomic block.")

                line_index = self._find_next_nonempty_line_index(i + 1)
                if line_index is None:
                    raise ValueError(f"Could not find LMAX for element {symbol}.")

                atom_count_index = self._find_next_nonempty_line_index(line_index + 1)
                if atom_count_index is None:
                    raise ValueError(f"Could not find atom count for element {symbol}.")

                atom_count_str = self.lines[atom_count_index].strip()
                try:
                    atom_count = int(atom_count_str)
                except ValueError:
                    raise ValueError(f"Invalid atom count '{atom_count_str}' for element {symbol}.")

                coord_line_index = atom_count_index + 1
                read_count = 0

                while coord_line_index < nlines and read_count < atom_count:
                    coord_line = self.lines[coord_line_index].strip()

                    if coord_line == "":
                        coord_line_index += 1
                        continue

                    values = self._extract_floats_from_line(coord_line)

                    if len(values) < 3:
                        raise ValueError(
                            f"Invalid coordinate line for element {symbol}: '{coord_line}'"
                        )

                    x, y, z = values[:3]

                    if not self.use_angstrom:
                        x *= self.BOHR_TO_ANGSTROM
                        y *= self.BOHR_TO_ANGSTROM
                        z *= self.BOHR_TO_ANGSTROM

                    parsed_atoms.append((symbol, x, y, z))
                    read_count += 1
                    coord_line_index += 1

                if read_count != atom_count:
                    raise ValueError(
                        f"Expected {atom_count} atoms for element {symbol}, but found {read_count}."
                    )

                i = coord_line_index
                continue

            i += 1

        return parsed_atoms

    def _extract_element_symbol_from_pseudopotential_line(self, line):
        """
        Extract atomic symbol from lines like:
            *C_VDB_PBE.psp FORMATTED
            *H_VDB_PBE.psp FORMATTED
        """
        cleaned = line.strip().lstrip("*")
        match = re.match(r"([A-Za-z]{1,2})", cleaned)

        if not match:
            raise ValueError(f"Could not determine atomic symbol from line: {line}")

        symbol = match.group(1)
        return symbol[0].upper() + symbol[1:].lower()

    def _extract_floats_from_line(self, line):
        """
        Extract floating-point values from a line.
        """
        float_strings = re.findall(r"[-+]?\d*\.?\d+(?:[DdEe][-+]?\d+)?", line)
        values = []

        for item in float_strings:
            values.append(float(item.replace("D", "E").replace("d", "e")))

        return values

    def _find_next_nonempty_line_index(self, start_index):
        """
        Find the next non-empty line index.
        """
        for index in range(start_index, len(self.lines)):
            if self.lines[index].strip() != "":
                return index
        return None

    def _build_cell_comment(self, cell_values):
        """
        Build the XYZ comment line with CELL information.
        """
        unit_label = "ANGSTROM" if self.use_angstrom else "BOHR->ANGSTROM"
        values_as_text = " ".join(f"{value:.6f}" for value in cell_values)
        return f"CELL: {values_as_text} | UNITS: {unit_label}"


class CpmdInputToXYZUI(CpmdInputToXYZConverter):
    """
    Toga-based UI for converting CPMD input files to XYZ files.
    """

    def __init__(self, *args):
        super().__init__()
        self.main_window = toga.Window(
            title="CPMD Input to XYZ Converter", size=(760, 200)
        )
        self.build_ui()
        self.main_window.show()

    async def select_input_file(self, widget, target_input, dialog_title="Select CPMD Input File"):
        """
        Open a file dialog to select the CPMD input file.
        """
        try:
            file_path = await self.main_window.dialog(
                toga.OpenFileDialog(title=dialog_title)
            )
            if file_path:
                target_input.value = str(file_path)
                self.file_path = str(file_path)
                return str(file_path)
            else:
                await self.main_window.dialog(
                    toga.InfoDialog("Warning", "No file was selected!")
                )
                return None
        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Failed to open file: {e}")
            )
            return None

    async def on_convert_press(self, widget):
        """
        Async handler for the Convert button.
        """
        input_filename = self.input_text_input.value.strip()
        output_filename = self.output_text_input.value.strip()

        if not input_filename:
            await self.main_window.dialog(
                toga.ErrorDialog("Input Error", "Please select an input file.")
            )
            return

        if not output_filename:
            await self.main_window.dialog(
                toga.ErrorDialog("Input Error", "Please specify an output file name.")
            )
            return

        output_dir = os.path.dirname(self.file_path) if hasattr(self, "file_path") and self.file_path else os.getcwd()
        output_filepath = os.path.join(output_dir, output_filename)

        self.output_file_label.text = "Processing CPMD input file..."
        await asyncio.sleep(0.01)

        success = await self.convert_cpmd_input_to_xyz(input_filename, output_filepath)

        if success:
            self.output_file_label.text = (
                f"Conversion successful. Output saved as: {output_filepath}"
            )
        else:
            self.output_file_label.text = "Conversion failed. Check the messages for details."

    def build_ui(self):
        """
        Build the main window layout.
        """
        main_box_style = Pack(direction=COLUMN, margin=10)
        row_style = Pack(direction=ROW, margin=5)
        label_style = Pack(width=160, margin=(0, 5))
        input_style = Pack(flex=1, margin=(0, 5))
        button_style = Pack(margin=5)

        # --- Input file selection ---
        input_box = toga.Box(style=row_style)
        input_label = toga.Label("Input CPMD file:", style=label_style)
        self.input_text_input = toga.TextInput(
            placeholder="Select CPMD input file", style=input_style
        )
        input_browse_button = toga.Button(
            "Browse",
            on_press=lambda w: asyncio.create_task(
                self.select_input_file(w, self.input_text_input)
            ),
            style=button_style
        )
        input_box.add(input_label)
        input_box.add(self.input_text_input)
        input_box.add(input_browse_button)

        # --- Output file name ---
        output_box = toga.Box(style=row_style)
        output_label = toga.Label("Output XYZ file name:", style=label_style)
        self.output_text_input = toga.TextInput(
            placeholder="Enter output file name", style=input_style
        )
        output_box.add(output_label)
        output_box.add(self.output_text_input)

        # --- Status label ---
        self.output_file_label = toga.Label("Output File: ", style=Pack(margin=5))

        # --- Convert button ---
        convert_button = toga.Button(
            "Convert",
            on_press=lambda w: asyncio.create_task(self.on_convert_press(w)),
            style=Pack(margin_left=300, width=100)
        )

        # --- Main layout ---
        main_box = toga.Box(style=main_box_style)
        main_box.add(input_box)
        main_box.add(output_box)
        main_box.add(self.output_file_label)
        main_box.add(convert_button)

        self.main_window.content = main_box
