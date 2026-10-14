import os
import re
from typing import List, Dict, Optional

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW


Atom = Dict[str, float]

# CODATA 2022 recommended value of the Bohr radius a0 (NIST, May 2024):
# a0 = 0.529177210544(82) Angstrom, relative standard uncertainty 1.6e-10.
BOHR_RADIUS_ANGSTROM = 0.529177210544

UNIT_ANGSTROM = "Angstrom"
UNIT_BOHR = "Bohr"
COORDINATE_UNITS = (UNIT_ANGSTROM, UNIT_BOHR)

# Decimal places written by the free-format writers (XYZ, Gaussian, HIN).
# Ten places keep an 8-decimal input bit-recoverable after an Angstrom -> bohr
# -> Angstrom (or bohr -> Angstrom -> bohr) round trip.
FREE_FORMAT_DECIMALS = 10

# Matches a Gaussian route "Units=AU", "Units=Bohr", "Units(AU,Deg)", ...
_GAUSSIAN_UNITS_RE = re.compile(r"\bunits\s*(?:=\s*\(?|\()\s*([^)\s]*)", re.IGNORECASE)


def convert_length(value: float, from_unit: str, to_unit: str) -> float:
    """Convert one Cartesian coordinate between Angstrom and bohr.

    A single correctly-rounded division (Angstrom -> bohr) or multiplication
    (bohr -> Angstrom) by ``BOHR_RADIUS_ANGSTROM`` is used, never a stored
    reciprocal, so the only rounding is the one the operation itself makes.
    Equal units return the value untouched.
    """
    for unit in (from_unit, to_unit):
        if unit not in COORDINATE_UNITS:
            raise ValueError(f"Unsupported coordinate unit: {unit}")
    value = float(value)
    if from_unit == to_unit:
        return value
    if from_unit == UNIT_ANGSTROM:
        return value / BOHR_RADIUS_ANGSTROM
    return value * BOHR_RADIUS_ANGSTROM


def convert_atom_units(atoms: List[Atom], from_unit: str, to_unit: str) -> List[Atom]:
    """Return a new atom list with x/y/z converted; the input list is not modified."""
    converted = []
    for atom in atoms:
        new_atom = dict(atom)
        for axis in ("x", "y", "z"):
            new_atom[axis] = convert_length(atom[axis], from_unit, to_unit)
        converted.append(new_atom)
    return converted


class CoordinateConverter:
    """Coordinate file converter based on a normalized atom list.

    Supported formats:
        XYZ, PDB, MOL, MOL2, Gaussian, ENT, HIN

    Supported coordinate units:
        Angstrom and bohr (atomic units), in both directions. PDB, ENT, MOL,
        MOL2 and HIN fix Angstrom in their specification, so bohr is only
        accepted on the XYZ and Gaussian side of a conversion.

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

    # Formats whose specification defines coordinates in Angstrom.
    ANGSTROM_ONLY_FORMATS = ("PDB", "ENT", "MOL", "MOL2", "HIN")

    def __init__(self):
        self.input_format = "XYZ"
        self.output_format = "PDB"
        self.input_unit = UNIT_ANGSTROM
        self.output_unit = UNIT_ANGSTROM
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

    async def _show_info(self, title: str, message: str):
        await self.main_window.dialog(toga.InfoDialog(title, message))

    async def _show_error(self, title: str, message: str):
        await self.main_window.dialog(toga.ErrorDialog(title, message))

    def guess_input_format_from_filename(self, filename: str) -> Optional[str]:
        ext = os.path.splitext(filename)[1].lower()
        return self.EXTENSION_TO_FORMAT.get(ext)

    def set_input_file(self, filepath: str):
        self.input_file = filepath
        guessed = self.guess_input_format_from_filename(filepath)
        if guessed and guessed in self.supported_formats:
            self.input_format = guessed
        if self.input_format == "Gaussian":
            # Gaussian declares its unit in the route section; XYZ declares none,
            # so an XYZ keeps whatever unit the user selected.
            declared = self.guess_gaussian_unit(filepath)
            if declared:
                self.input_unit = declared
        if not self.unit_selectable(self.input_format):
            self.input_unit = UNIT_ANGSTROM

    def unit_selectable(self, fmt: str) -> bool:
        """True when ``fmt`` may hold bohr coordinates (XYZ, Gaussian)."""
        return fmt not in self.ANGSTROM_ONLY_FORMATS

    @staticmethod
    def guess_gaussian_unit(filepath: str) -> Optional[str]:
        """Return the unit declared by a Gaussian route section, or None if unreadable.

        ``Units=AU`` / ``Units=Bohr`` (also in the ``Units(AU,Deg)`` form) mean
        bohr; a route without the keyword means Gaussian's default, Angstrom.
        """
        try:
            with open(filepath, "r", encoding="utf-8") as file_handle:
                lines = [line.strip() for line in file_handle]
        except (OSError, UnicodeDecodeError):
            return None

        route_parts = []
        in_route = False
        for line in lines:
            if line.startswith("#"):
                in_route = True
            if in_route:
                if not line:
                    break
                route_parts.append(line)

        for match in _GAUSSIAN_UNITS_RE.finditer(" ".join(route_parts)):
            options = {option.strip().lower() for option in match.group(1).split(",")}
            if options & {"au", "bohr"}:
                return UNIT_BOHR
        return UNIT_ANGSTROM

    def default_output_path(self) -> str:
        """``<input>_converted[_<unit>].<ext>``; the unit tag appears only when units change."""
        if not self.input_file:
            return ""
        base, _ = os.path.splitext(self.input_file)
        extension = self.FORMAT_TO_EXTENSION.get(self.output_format, self.output_format.lower())
        suffix = "_converted"
        if self.input_unit != self.output_unit:
            suffix += f"_{self.output_unit.lower()}"
        return f"{base}{suffix}.{extension}"

    def ensure_output_file(self):
        if not self.output_file:
            self.output_file = self.default_output_path()

    def validate_conversion(self):
        """Raise ``ValueError`` describing the first problem with the current settings."""
        if not self.input_file:
            raise ValueError("Please select an input file.")
        for label, unit in (("input", self.input_unit), ("output", self.output_unit)):
            if unit not in COORDINATE_UNITS:
                raise ValueError(f"Unsupported {label} coordinate unit: {unit}")
        for label, fmt, unit in (
            ("Input", self.input_format, self.input_unit),
            ("Output", self.output_format, self.output_unit),
        ):
            if unit != UNIT_ANGSTROM and not self.unit_selectable(fmt):
                raise ValueError(
                    f"{label} format {fmt} stores coordinates in Angstrom by specification; "
                    f"use XYZ or Gaussian for {unit} coordinates."
                )
        if self.input_format == self.output_format and self.input_unit == self.output_unit:
            raise ValueError(
                "Input and output formats and units are the same. "
                "Please select a different format or a different unit."
            )
        self.ensure_output_file()
        same_file = os.path.normcase(os.path.abspath(self.output_file)) == os.path.normcase(
            os.path.abspath(self.input_file)
        )
        if same_file:
            raise ValueError("The output file would overwrite the input file. Please choose another name.")

    def convert_file(self) -> int:
        """Read, convert units, write. Returns the number of atoms written.

        Headless: touches no widget, raises ``ValueError`` on any failure.
        """
        self.validate_conversion()
        atoms = self.read_atoms(self.input_format)
        if not atoms:
            raise ValueError("No atoms could be read from the input file.")
        atoms = convert_atom_units(atoms, self.input_unit, self.output_unit)
        self.write_atoms(atoms, self.output_format)
        return len(atoms)

    async def select_input_file(self, widget):
        """Open a file dialog to select the input file."""
        try:
            file_dialog = toga.OpenFileDialog(title="Select Input File")
            selected = await self.main_window.dialog(file_dialog)
            input_file = self._normalize_dialog_result(selected)

            if not input_file:
                await self._show_info("Warning", "No file was selected.")
                return

            self.set_input_file(input_file)
            self.input_file_input.value = self.input_file
            self.sync_input_selection()
            self.update_output_options()
            self.update_output_filename()
        except Exception as exc:
            await self._show_error("Error", f"Failed to open file: {exc}")

    async def select_output_file(self, widget):
        """Allow the user to override the automatically generated output filename."""
        if not self.input_file:
            await self._show_error("Error", "Please select an input file first.")
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
            await self._show_error("Error", f"Failed to select output file: {exc}")

    def update_output_filename(self):
        """Update the output filename from input file and selected output format."""
        self.output_file = self.default_output_path()
        self.output_file_label.text = f"Output File: {self.output_file}"

    async def perform_conversion(self, widget):
        """Perform the file format and coordinate unit conversion."""
        try:
            atom_count = self.convert_file()
        except Exception as exc:
            await self._show_error("Error", f"Conversion failed: {exc}")
            return

        await self._show_info(
            "Success",
            f"Converted {atom_count} atoms from {self.input_format} ({self.input_unit}) "
            f"to {self.output_format} ({self.output_unit}).\n\n"
            f"Saved as:\n{self.output_file}",
        )

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
                file_handle.write(f"Converted file (coordinates in {self.output_unit.lower()})\n")
                for atom in atoms:
                    file_handle.write(f"{atom['symbol']:<2} {self._free_format_xyz(atom)}\n")
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
                        f"atom {index:5d} {atom['symbol']:<2} 0 "
                        f"{atom['x']:.{FREE_FORMAT_DECIMALS}f} {atom['y']:.{FREE_FORMAT_DECIMALS}f} "
                        f"{atom['z']:.{FREE_FORMAT_DECIMALS}f}\n"
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
                # Gaussian reads Angstrom unless the route says otherwise.
                units_keyword = " Units=Bohr" if self.output_unit == UNIT_BOHR else ""
                file_handle.write(f"#p hf/3-21g{units_keyword}\n\n")
                file_handle.write("Converted structure\n\n")
                file_handle.write("0 1\n")
                for atom in atoms:
                    file_handle.write(f"{atom['symbol']:<2} {self._free_format_xyz(atom)}\n")
                file_handle.write("\n")
        except Exception as exc:
            raise ValueError(f"Failed to write Gaussian file: {exc}") from exc

    @staticmethod
    def _free_format_xyz(atom: Atom) -> str:
        width = FREE_FORMAT_DECIMALS + 8
        return " ".join(f"{atom[axis]:{width}.{FREE_FORMAT_DECIMALS}f}" for axis in ("x", "y", "z"))


class CoordinateConverterUI(CoordinateConverter):
    def __init__(self, app):
        super().__init__()
        self.app = app
        self.main_window = toga.Window(title="Coordinate Converter", size=(760, 360))
        self.layout_main_window()
        self.initialize_defaults()
        self.main_window.show()

    def initialize_defaults(self):
        self.input_format_selection.value = self.input_format
        self.update_output_options()
        self.output_format_selection.value = self.output_format
        self.sync_unit_controls()
        self.update_output_filename()

    def layout_main_window(self):
        main_box_style = Pack(direction=COLUMN, margin=10)
        box_style = Pack(direction=ROW, margin=5)
        label_style = Pack(width=130, margin=(0, 5))
        unit_label_style = Pack(margin=(5, 5))
        input_style = Pack(flex=1, margin=(0, 5))
        unit_style = Pack(width=120, margin=(0, 5))
        button_style = Pack(margin=5)

        input_format_box = toga.Box(style=box_style)
        input_format_label = toga.Label("Input Format:", style=label_style)
        self.input_format_selection = toga.Selection(
            items=self.supported_formats,
            style=input_style,
            on_change=self.on_input_format_change,
        )
        self.input_unit_selection = toga.Selection(
            items=list(COORDINATE_UNITS),
            style=unit_style,
            on_change=self.on_input_unit_change,
        )
        input_format_box.add(input_format_label)
        input_format_box.add(self.input_format_selection)
        input_format_box.add(toga.Label("Unit:", style=unit_label_style))
        input_format_box.add(self.input_unit_selection)

        output_format_box = toga.Box(style=box_style)
        output_format_label = toga.Label("Output Format:", style=label_style)
        self.output_format_selection = toga.Selection(
            items=self.get_output_formats(),
            style=input_style,
            on_change=self.on_output_format_change,
        )
        self.output_unit_selection = toga.Selection(
            items=list(COORDINATE_UNITS),
            style=unit_style,
            on_change=self.on_output_unit_change,
        )
        output_format_box.add(output_format_label)
        output_format_box.add(self.output_format_selection)
        output_format_box.add(toga.Label("Unit:", style=unit_label_style))
        output_format_box.add(self.output_unit_selection)

        unit_hint = toga.Label(
            f"1 bohr = {BOHR_RADIUS_ANGSTROM} Å (CODATA 2022). "
            "PDB, ENT, MOL, MOL2 and HIN are always in Å; bohr is available for XYZ and Gaussian.",
            style=Pack(margin=(0, 10), color="#555555"),
        )

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
        main_box.add(unit_hint)
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
        self.sync_unit_controls()

    def sync_unit_controls(self):
        """Show the current units and lock Angstrom-only formats to Angstrom.

        Setting ``.value`` fires ``on_*_unit_change``, which only stores the
        same value back, so this cannot recurse.
        """
        for fmt, unit_attr, selection in (
            (self.input_format, "input_unit", self.input_unit_selection),
            (self.output_format, "output_unit", self.output_unit_selection),
        ):
            selectable = self.unit_selectable(fmt)
            if not selectable:
                setattr(self, unit_attr, UNIT_ANGSTROM)
            selection.value = getattr(self, unit_attr)
            selection.enabled = selectable

    def on_input_format_change(self, widget):
        selected = self._selection_value(widget)
        if selected:
            self.input_format = selected
        self.update_output_options()
        self.sync_unit_controls()
        self.update_output_filename()

    def on_output_format_change(self, widget):
        selected = self._selection_value(widget)
        if selected:
            self.output_format = selected
        self.sync_unit_controls()
        self.update_output_filename()

    def on_input_unit_change(self, widget):
        selected = self._selection_value(widget)
        if selected:
            self.input_unit = selected
        self.update_output_filename()

    def on_output_unit_change(self, widget):
        selected = self._selection_value(widget)
        if selected:
            self.output_unit = selected
        self.update_output_filename()

    def update_output_options(self):
        # Every format is a valid target: XYZ -> XYZ is a pure unit conversion.
        formats = self.get_output_formats()
        if self.output_format not in formats and formats:
            self.output_format = formats[0]
        if formats:
            try:
                self.output_format_selection.value = self.output_format
            except Exception:
                pass

    def get_output_formats(self):
        return list(self.supported_formats)
