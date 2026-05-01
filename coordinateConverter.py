import os
from typing import List, Dict, Optional

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW


Atom = Dict[str, float]


class CoordinateConverter:
    """Coordinate file converter based on a normalized atom list.

    Supported formats:
        XYZ, PDB, MOL, MOL2, Gaussian, ENT, HIN

    Notes:
        * Only atom labels and Cartesian coordinates are preserved.
        * Bond orders, charges, residue metadata, and force-field specific
          information are not reconstructed when absent in the source file.
        * ENT is treated as a PDB-style format.
    """

    EXTENSION_TO_FORMAT = {
        ".xyz": "XYZ",
        ".pdb": "PDB",
        ".ent": "ENT",
        ".mol": "MOL",
        ".mol2": "MOL2",
        ".gjf": "Gaussian",
        ".com": "Gaussian",
        ".gau": "Gaussian",
        ".hin": "HIN",
    }

    FORMAT_TO_EXTENSION = {
        "XYZ": "xyz",
        "PDB": "pdb",
        "MOL": "mol",
        "MOL2": "mol2",
        "Gaussian": "gjf",
        "ENT": "ent",
        "HIN": "hin",
    }

    def __init__(self):
        self.input_format = "XYZ"
        self.output_format = "PDB"
        self.input_file = ""
        self.output_file = ""
        self.supported_formats = ["XYZ", "PDB", "MOL", "MOL2", "Gaussian", "ENT", "HIN"]

    # ------------------------------------------------------------------
    # GUI helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _normalize_dialog_result(dialog_result) -> str:
        """Convert Toga dialog result into a filesystem path string."""
        if not dialog_result:
            return ""

        # Some backends may return a list/tuple.
        if isinstance(dialog_result, (list, tuple)):
            if not dialog_result:
                return ""
            dialog_result = dialog_result[0]

        # Toga may return a Path-like object.
        path_value = getattr(dialog_result, "path", dialog_result)
        return str(path_value)

    def guess_input_format_from_filename(self, filename: str) -> Optional[str]:
        ext = os.path.splitext(filename)[1].lower()
        return self.EXTENSION_TO_FORMAT.get(ext)

    def set_input_file(self, filepath: str):
        self.input_file = filepath
        guessed = self.guess_input_format_from_filename(filepath)
        if guessed and guessed in self.supported_formats:
            self.input_format = guessed

    def ensure_output_file(self):
        if not self.output_file:
            self.update_output_filename()

    async def select_input_file(self, widget):
        """Open a file dialog to select the input file."""
        try:
            file_dialog = toga.OpenFileDialog(title="Select Input File")
            selected = await self.main_window.dialog(file_dialog)
            input_file = self._normalize_dialog_result(selected)

            if not input_file:
                await self.main_window.info_dialog("Warning", "No file was selected.")
                return

            self.set_input_file(input_file)
            self.input_file_input.value = self.input_file
            self.sync_input_selection()
            self.update_output_options()
            self.update_output_filename()
        except Exception as exc:
            await self.main_window.error_dialog("Error", f"Failed to open file: {exc}")

    async def select_output_file(self, widget):
        """Allow the user to override the automatically generated output filename."""
        if not self.input_file:
            await self.main_window.error_dialog("Error", "Please select an input file first.")
            return

        try:
            suggested = os.path.basename(self.output_file) if self.output_file else f"converted.{self.FORMAT_TO_EXTENSION[self.output_format]}"
            save_dialog = toga.SaveFileDialog(title="Save Converted File", suggested_filename=suggested)
            selected = await self.main_window.dialog(save_dialog)
            output_file = self._normalize_dialog_result(selected)
            if output_file:
                self.output_file = output_file
                self.output_file_label.text = f"Output File: {self.output_file}"
        except Exception as exc:
            await self.main_window.error_dialog("Error", f"Failed to select output file: {exc}")

    def update_output_filename(self):
        """Update the output filename from input file and selected output format."""
        if not self.input_file:
            self.output_file = ""
            self.output_file_label.text = "Output File: "
            return

        base, _ = os.path.splitext(self.input_file)
        extension = self.FORMAT_TO_EXTENSION.get(self.output_format, self.output_format.lower())
        self.output_file = f"{base}_converted.{extension}"
        self.output_file_label.text = f"Output File: {self.output_file}"

    async def perform_conversion(self, widget):
        """Perform the file format conversion."""
        if not self.input_file:
            await self.main_window.error_dialog("Error", "Please select an input file.")
            return

        if self.input_format == self.output_format:
            await self.main_window.error_dialog(
                "Error",
                "Input and output formats are the same. Please select different formats.",
            )
            return

        try:
            self.ensure_output_file()
            atoms = self.read_atoms(self.input_format)
            if not atoms:
                raise ValueError("No atoms could be read from the input file.")
            self.write_atoms(atoms, self.output_format)
            await self.main_window.info_dialog(
                "Success",
                f"File converted successfully and saved as:\n{self.output_file}",
            )
        except Exception as exc:
            await self.main_window.error_dialog("Error", f"Conversion failed: {exc}")

    # ------------------------------------------------------------------
    # General conversion engine
    # ------------------------------------------------------------------
    def read_atoms(self, fmt: str) -> List[Atom]:
        readers = {
            "XYZ": self.read_xyz,
            "PDB": self.read_pdb,
            "ENT": self.read_ent,
            "MOL": self.read_mol,
            "MOL2": self.read_mol2,
            "Gaussian": self.read_gaussian,
            "HIN": self.read_hin,
        }
        reader = readers.get(fmt)
        if reader is None:
            raise ValueError(f"Unsupported input format: {fmt}")
        return reader()

    def write_atoms(self, atoms: List[Atom], fmt: str):
        writers = {
            "XYZ": self.write_xyz,
            "PDB": self.write_pdb,
            "ENT": self.write_ent,
            "MOL": self.write_mol,
            "MOL2": self.write_mol2,
            "Gaussian": self.write_gaussian,
            "HIN": self.write_hin,
        }
        writer = writers.get(fmt)
        if writer is None:
            raise ValueError(f"Unsupported output format: {fmt}")
        writer(atoms)

    @staticmethod
    def _make_atom(symbol: str, x: float, y: float, z: float) -> Atom:
        return {
            "symbol": CoordinateConverter.normalize_symbol(symbol),
            "x": float(x),
            "y": float(y),
            "z": float(z),
        }

    @staticmethod
    def normalize_symbol(symbol: str) -> str:
        symbol = (symbol or "X").strip()
        if not symbol:
            return "X"
        if len(symbol) == 1:
            return symbol.upper()
        return symbol[0].upper() + symbol[1:].lower()

    # ------------------------------------------------------------------
    # Readers
    # ------------------------------------------------------------------
    def read_xyz(self) -> List[Atom]:
        atoms: List[Atom] = []
        try:
            with open(self.input_file, "r", encoding="utf-8") as file_handle:
                lines = [line.rstrip("\n") for line in file_handle]

            if len(lines) < 2:
                raise ValueError("XYZ file is too short.")

            expected_atoms = int(lines[0].strip())
            for index, line in enumerate(lines[2:], start=3):
                if not line.strip():
                    continue
                tokens = line.split()
                if len(tokens) < 4:
                    raise ValueError(f"Invalid XYZ line {index}: {line}")
                atoms.append(self._make_atom(tokens[0], tokens[1], tokens[2], tokens[3]))

            if expected_atoms != len(atoms):
                raise ValueError(
                    f"XYZ atom count mismatch: header says {expected_atoms}, read {len(atoms)} atoms."
                )
        except Exception as exc:
            raise ValueError(f"Failed to read XYZ file: {exc}") from exc
        return atoms

    def read_pdb(self) -> List[Atom]:
        atoms: List[Atom] = []
        try:
            with open(self.input_file, "r", encoding="utf-8") as file_handle:
                for line in file_handle:
                    if line.startswith(("ATOM", "HETATM")):
                        symbol = line[76:78].strip() or line[12:16].strip()
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        atoms.append(self._make_atom(symbol, x, y, z))
        except Exception as exc:
            raise ValueError(f"Failed to read PDB file: {exc}") from exc
        return atoms

    def read_ent(self) -> List[Atom]:
        return self.read_pdb()

    def read_hin(self) -> List[Atom]:
        atoms: List[Atom] = []
        try:
            with open(self.input_file, "r", encoding="utf-8") as file_handle:
                for line in file_handle:
                    if not line.startswith("atom"):
                        continue
                    tokens = line.strip().split()
                    if len(tokens) < 7:
                        raise ValueError(f"Invalid HIN line: {line.strip()}")
                    atoms.append(self._make_atom(tokens[2], tokens[4], tokens[5], tokens[6]))
        except Exception as exc:
            raise ValueError(f"Failed to read HIN file: {exc}") from exc
        return atoms

    def read_mol(self) -> List[Atom]:
        atoms: List[Atom] = []
        try:
            with open(self.input_file, "r", encoding="utf-8") as file_handle:
                lines = [line.rstrip("\n") for line in file_handle]

            if len(lines) < 4:
                raise ValueError("MOL file is too short.")

            counts_line = lines[3]
            atom_count = int(counts_line[0:3].strip())
            atom_block = lines[4:4 + atom_count]
            if len(atom_block) != atom_count:
                raise ValueError("MOL atom block is incomplete.")

            for line in atom_block:
                x = float(line[0:10].strip())
                y = float(line[10:20].strip())
                z = float(line[20:30].strip())
                symbol = line[31:34].strip()
                atoms.append(self._make_atom(symbol, x, y, z))
        except Exception as exc:
            raise ValueError(f"Failed to read MOL file: {exc}") from exc
        return atoms

    def read_mol2(self) -> List[Atom]:
        atoms: List[Atom] = []
        try:
            with open(self.input_file, "r", encoding="utf-8") as file_handle:
                in_atom_section = False
                for line in file_handle:
                    striped = line.strip()
                    if striped.startswith("@<TRIPOS>ATOM"):
                        in_atom_section = True
                        continue
                    if striped.startswith("@<TRIPOS>") and in_atom_section:
                        break
                    if in_atom_section and striped:
                        tokens = striped.split()
                        if len(tokens) < 6:
                            raise ValueError(f"Invalid MOL2 atom line: {striped}")
                        symbol = tokens[5].split(".")[0]
                        atoms.append(self._make_atom(symbol, tokens[2], tokens[3], tokens[4]))
        except Exception as exc:
            raise ValueError(f"Failed to read MOL2 file: {exc}") from exc
        return atoms

    def read_gaussian(self) -> List[Atom]:
        atoms: List[Atom] = []
        try:
            with open(self.input_file, "r", encoding="utf-8") as file_handle:
                lines = [line.rstrip("\n") for line in file_handle]

            charge_mult_index = None
            for index, line in enumerate(lines):
                tokens = line.split()
                if len(tokens) == 2:
                    try:
                        int(tokens[0])
                        int(tokens[1])
                        charge_mult_index = index
                        break
                    except ValueError:
                        continue

            if charge_mult_index is None:
                raise ValueError("Could not locate the charge/multiplicity line in Gaussian input.")

            for line in lines[charge_mult_index + 1:]:
                striped = line.strip()
                if not striped:
                    break
                tokens = striped.split()
                if len(tokens) < 4:
                    raise ValueError(f"Invalid Gaussian coordinate line: {line}")
                atoms.append(self._make_atom(tokens[0], tokens[1], tokens[2], tokens[3]))
        except Exception as exc:
            raise ValueError(f"Failed to read Gaussian file: {exc}") from exc
        return atoms

    # ------------------------------------------------------------------
    # Writers
    # ------------------------------------------------------------------
    def write_xyz(self, atoms: List[Atom]):
        try:
            with open(self.output_file, "w", encoding="utf-8") as file_handle:
                file_handle.write(f"{len(atoms)}\n")
                file_handle.write("Converted file\n")
                for atom in atoms:
                    file_handle.write(
                        f"{atom['symbol']:<2} {atom['x']:16.8f} {atom['y']:16.8f} {atom['z']:16.8f}\n"
                    )
        except Exception as exc:
            raise ValueError(f"Failed to write XYZ file: {exc}") from exc

    def write_pdb(self, atoms: List[Atom]):
        try:
            with open(self.output_file, "w", encoding="utf-8") as file_handle:
                file_handle.write("REMARK   Converted file\n")
                for index, atom in enumerate(atoms, start=1):
                    file_handle.write(
                        f"ATOM  {index:5d} {atom['symbol']:<4} MOL A{1:4d}    "
                        f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
                        f"  1.00  0.00          {atom['symbol']:>2}\n"
                    )
                file_handle.write("END\n")
        except Exception as exc:
            raise ValueError(f"Failed to write PDB file: {exc}") from exc

    def write_ent(self, atoms: List[Atom]):
        self.write_pdb(atoms)

    def write_hin(self, atoms: List[Atom]):
        try:
            with open(self.output_file, "w", encoding="utf-8") as file_handle:
                file_handle.write("mol 1 Converted file\n")
                for index, atom in enumerate(atoms, start=1):
                    file_handle.write(
                        f"atom {index:5d} {atom['symbol']:<2} 0 {atom['x']:.6f} {atom['y']:.6f} {atom['z']:.6f}\n"
                    )
                file_handle.write("endmol 1\n")
        except Exception as exc:
            raise ValueError(f"Failed to write HIN file: {exc}") from exc

    def write_mol(self, atoms: List[Atom]):
        try:
            with open(self.output_file, "w", encoding="utf-8") as file_handle:
                file_handle.write("Converted file\n")
                file_handle.write("CoordinateConverter\n")
                file_handle.write("Generated from Cartesian coordinates\n")
                file_handle.write(f"{len(atoms):3d}{0:3d}  0  0  0  0            999 V2000\n")
                for atom in atoms:
                    file_handle.write(
                        f"{atom['x']:10.4f}{atom['y']:10.4f}{atom['z']:10.4f} {atom['symbol']:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n"
                    )
                file_handle.write("M  END\n")
        except Exception as exc:
            raise ValueError(f"Failed to write MOL file: {exc}") from exc

    def write_mol2(self, atoms: List[Atom]):
        try:
            with open(self.output_file, "w", encoding="utf-8") as file_handle:
                file_handle.write("@<TRIPOS>MOLECULE\n")
                file_handle.write("Converted file\n")
                file_handle.write(f"{len(atoms)} 0 0 0 0\n")
                file_handle.write("SMALL\n")
                file_handle.write("NO_CHARGES\n\n")
                file_handle.write("@<TRIPOS>ATOM\n")
                for index, atom in enumerate(atoms, start=1):
                    atom_type = atom['symbol']
                    file_handle.write(
                        f"{index:7d} {atom['symbol']}{index:<3d} "
                        f"{atom['x']:10.4f} {atom['y']:10.4f} {atom['z']:10.4f} "
                        f"{atom_type:<6} 1 MOL 0.0000\n"
                    )
                file_handle.write("@<TRIPOS>BOND\n")
        except Exception as exc:
            raise ValueError(f"Failed to write MOL2 file: {exc}") from exc

    def write_gaussian(self, atoms: List[Atom]):
        try:
            with open(self.output_file, "w", encoding="utf-8") as file_handle:
                file_handle.write("%chk=converted.chk\n")
                file_handle.write("#p hf/3-21g\n\n")
                file_handle.write("Converted structure\n\n")
                file_handle.write("0 1\n")
                for atom in atoms:
                    file_handle.write(
                        f"{atom['symbol']:<2} {atom['x']:16.8f} {atom['y']:16.8f} {atom['z']:16.8f}\n"
                    )
                file_handle.write("\n")
        except Exception as exc:
            raise ValueError(f"Failed to write Gaussian file: {exc}") from exc


class CoordinateConverterUI(CoordinateConverter):
    def __init__(self, app):
        super().__init__()
        self.app = app
        self.main_window = toga.Window(title="Coordinate Converter", size=(760, 320))
        self.layout_main_window()
        self.initialize_defaults()
        self.main_window.show()

    def initialize_defaults(self):
        self.input_format_selection.value = self.input_format
        self.update_output_options()
        self.output_format_selection.value = self.output_format
        self.update_output_filename()

    def layout_main_window(self):
        main_box_style = Pack(direction=COLUMN, margin=10)
        box_style = Pack(direction=ROW, margin=5)
        label_style = Pack(width=130, margin=(0, 5))
        input_style = Pack(flex=1, margin=(0, 5))
        button_style = Pack(margin=5)

        input_format_box = toga.Box(style=box_style)
        input_format_label = toga.Label("Input Format:", style=label_style)
        self.input_format_selection = toga.Selection(
            items=self.supported_formats,
            style=input_style,
            on_change=self.on_input_format_change,
        )
        input_format_box.add(input_format_label)
        input_format_box.add(self.input_format_selection)

        output_format_box = toga.Box(style=box_style)
        output_format_label = toga.Label("Output Format:", style=label_style)
        self.output_format_selection = toga.Selection(
            items=self.get_output_formats(),
            style=input_style,
            on_change=self.on_output_format_change,
        )
        output_format_box.add(output_format_label)
        output_format_box.add(self.output_format_selection)

        input_file_box = toga.Box(style=box_style)
        input_file_label = toga.Label("Input File:", style=label_style)
        self.input_file_input = toga.TextInput(
            readonly=True,
            placeholder="Select input file",
            style=input_style,
        )
        browse_button = toga.Button("Browse", on_press=self.select_input_file, style=button_style)
        input_file_box.add(input_file_label)
        input_file_box.add(self.input_file_input)
        input_file_box.add(browse_button)

        output_file_box = toga.Box(style=box_style)
        output_file_title = toga.Label("Output File:", style=label_style)
        self.output_file_label = toga.Label("", style=Pack(flex=1, margin=(5, 5)))
        choose_output_button = toga.Button(
            "Choose...",
            on_press=self.select_output_file,
            style=button_style,
        )
        output_file_box.add(output_file_title)
        output_file_box.add(self.output_file_label)
        output_file_box.add(choose_output_button)

        button_row = toga.Box(style=Pack(direction=ROW, margin_top=10))
        convert_button = toga.Button("Convert", on_press=self.perform_conversion, style=button_style)
        button_row.add(convert_button)

        main_box = toga.Box(style=main_box_style)
        main_box.add(input_format_box)
        main_box.add(output_format_box)
        main_box.add(input_file_box)
        main_box.add(output_file_box)
        main_box.add(button_row)

        self.main_window.content = main_box

    def _selection_value(self, widget) -> str:
        value = widget.value
        return str(value) if value is not None else ""

    def sync_input_selection(self):
        try:
            self.input_format_selection.value = self.input_format
        except Exception:
            pass

    def on_input_format_change(self, widget):
        selected = self._selection_value(widget)
        if selected:
            self.input_format = selected
        self.update_output_options()
        self.update_output_filename()

    def on_output_format_change(self, widget):
        selected = self._selection_value(widget)
        if selected:
            self.output_format = selected
        self.update_output_filename()

    def update_output_options(self):
        formats = [fmt for fmt in self.supported_formats if fmt != self.input_format]
        self.output_format_selection.items = formats
        if self.output_format not in formats and formats:
            self.output_format = formats[0]
        if formats:
            try:
                self.output_format_selection.value = self.output_format
            except Exception:
                pass

    def get_output_formats(self):
        return [fmt for fmt in self.supported_formats if fmt != self.input_format]
