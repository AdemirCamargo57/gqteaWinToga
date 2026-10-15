import os
import math
import random
import time
import numpy as np
from scipy.spatial.transform import Rotation
import toga
from toga.style import Pack
from toga.style.pack import LEFT, CENTER, ROW, COLUMN
import asyncio
from help import AtomicData, HelpGqteaWin
from solventPacking import (
    ENGINE_LEGACY, ENGINE_NUMPY, MoleculeGrid, RotationStream, atoms_clash, resolve_engine,
)
from systemProperties import parse_integer
from solvationBoxCommon import SolvationExtrasEngine, SolvationExtrasUI


class SingleSoluteSolventBox(SolvationExtrasEngine):
    """Build an orthorhombic solvent box containing one optional solute.

    Revised algorithm:
    - center solvent at origin
    - optionally center solute at box center
    - derive the target number of solvent molecules from the target density
    - insert solvent molecules at random positions with optional random rotations
    - reject insertions that clash with the solute or with previously accepted solvent molecules
    - save the accepted structure and report insertion statistics
    """

    AMU_TO_GRAMS = 1.66053906660e-24

    # Fallback radii in Å used if AtomicData does not provide them.
    DEFAULT_VDW_RADII = {
        "H": 1.20,
        "He": 1.40,
        "Li": 1.82,
        "Be": 1.53,
        "B": 1.92,
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "F": 1.47,
        "Ne": 1.54,
        "Na": 2.27,
        "Mg": 1.73,
        "Al": 1.84,
        "Si": 2.10,
        "P": 1.80,
        "S": 1.80,
        "Cl": 1.75,
        "Ar": 1.88,
        "K": 2.75,
        "Ca": 2.31,
        "Br": 1.85,
        "I": 1.98,
    }

    def __init__(self):
        self.solvent_coords = []
        self.solute_coords = []
        self.centered_solvent_coords = []
        self.centered_solute_coords = []
        self.positioned_solute_coords = []
        self.lattice_vectors = []
        self.spacing = 0.0              # retained for backward UI compatibility; used as optional safety padding
        self.data = []
        self.output_dir = ""
        self.min_distance = None
        self.target_density = 1.0
        self.max_attempts = 50000
        self.vdw_scale = 0.80
        self.random_seed = None
        self.solvent_insertion_count = 0
        self.solvent_attempts_used = 0
        self.failed_insertions = 0
        self.final_density = 0.0
        self.num_target_solvent_molecules = 0
        self.solvent_extent = 0.0
        self.solvent_mass_amu = 0.0
        self.solute_mass_amu = 0.0
        self.solvent_bounding_radius = 0.0
        self.solute_bounding_radius = 0.0
        self.use_periodic_minimum_image = True

        # Packing options (set from the widgets by read_params; plain
        # attributes so the engine runs headless and in a worker thread).
        self.rotate_solvent = True
        self.include_solute = False
        # Charges per species; the multiplicity (in the extras) is never inferred.
        self.solute_charge = 0
        self.solvent_charge = 0
        # Pre-optimization, properties, progress, molecules (solvationBoxCommon).
        self.init_solvation_extras()

    async def read_params(self, widget):
        """Read and validate every input; True when the run can start."""
        async def read_input(text_input, field_name, expected_type, required=True, default=None):
            value = text_input.value.strip()
            if not value:
                if required:
                    await self.main_window.dialog(toga.ErrorDialog(
                        "Error", f"Please input a valid value for {field_name}."
                    ))
                    return None
                return default
            try:
                if expected_type == list:
                    labels = [float(label) for label in value.split()]
                    if len(labels) != 3:
                        await self.main_window.dialog(toga.InfoDialog(
                            "Info", "Please input exactly three lattice vectors."
                        ))
                        return None
                    return labels
                return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(toga.ErrorDialog(
                    "Error", f"Invalid format for {field_name}: {e}"
                ))
                return None

        self.lattice_vectors = await read_input(
            self.textInput_lattice_vectors, "box lattice vectors", list
        )
        if self.lattice_vectors is None:
            return False
        if any(length <= 0.0 for length in self.lattice_vectors):
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "All three box lattice vectors must be positive."
            ))
            self.lattice_vectors = []
            return False

        self.spacing = await read_input(
            self.textInput_spacing, "legacy spacing / extra padding", float, required=False, default=0.0
        )
        if self.spacing is None:
            return False

        self.target_density = await read_input(
            self.textInput_target_density, "target density", float
        )
        if self.target_density is None or self.target_density <= 0.0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "Target density must be a positive number."
            ))
            return False

        self.max_attempts = await read_input(
            self.textInput_max_attempts, "maximum insertion attempts", int, required=False, default=50000
        )
        if self.max_attempts is None or self.max_attempts <= 0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "Maximum insertion attempts must be a positive integer."
            ))
            return False

        self.vdw_scale = await read_input(
            self.textInput_vdw_scale, "vdW scaling factor", float, required=False, default=0.80
        )
        if self.vdw_scale is None or self.vdw_scale <= 0.0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "vdW scaling factor must be positive."
            ))
            return False

        seed_value = await read_input(
            self.textInput_seed, "random seed", int, required=False, default=None
        )
        if seed_value is None and self.textInput_seed.value.strip():
            return False  # read_input already explained the bad value
        self.random_seed = seed_value
        self.use_periodic_minimum_image = bool(self.switch_minimum_image.value)
        self.rotate_solvent = bool(self.switch.value)
        self.include_solute = bool(self.switch_include_solute.value)

        if self.include_solute:
            self.min_distance = await read_input(
                self.textInput_min_distance, "minimum distance", float, required=False, default=0.0
            )
            if self.min_distance is None or self.min_distance < 0.0:
                await self.main_window.dialog(toga.ErrorDialog(
                    "Error", "Minimum distance must be zero or positive."
                ))
                return False
        else:
            self.min_distance = 0.0

        # Pre-optimization, charge/spin and engine options.
        if not await self.read_extras_params():
            return False

        update_text = (
            f"{'Box lattice vectors:':<32} {self.lattice_vectors}\n"
            f"{'Extra wall padding (Å):':<32} {self.spacing:.3f}\n"
            f"{'Target density (g/cm³):':<32} {self.target_density:.4f}\n"
            f"{'Max insertion attempts:':<32} {self.max_attempts}\n"
            f"{'vdW scale:':<32} {self.vdw_scale:.3f}\n"
        )
        if self.random_seed is not None:
            update_text += f"{'Random seed:':<32} {self.random_seed}\n"
        update_text += f"{'Periodic minimum-image clash:':<32} {'ON' if self.use_periodic_minimum_image else 'OFF'}\n"
        if self.include_solute:
            update_text += f"{'Extra solute clearance (Å):':<32} {self.min_distance:.3f}\n"
            update_text += f"{'Solute charge:':<32} {self.solute_charge:+d}\n"
        update_text += f"{'Solvent molecule charge:':<32} {self.solvent_charge:+d}\n"
        update_text += self.extras_parameter_text()
        self.multi_line_text.value += update_text
        return True

    async def read_solvent_xyz(self, widget):
        self.solvent_coords = []

        try:
            solvent_file = await self.main_window.dialog(
                toga.OpenFileDialog(title="Open Solvent XYZ file")
            )

            if not solvent_file:
                await self.main_window.dialog(
                    toga.InfoDialog("Warning", "No file was selected!")
                )
                return

            # Normalize path
            self.solvent_file = str(solvent_file)
            self.textInput_solvent.value = self.solvent_file
            self.output_dir = os.path.dirname(self.solvent_file)

        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Failed to open file: {e}")
            )
            return

        try:
            with open(self.solvent_file, "r", encoding="utf-8") as f:
                lines = f.readlines()

            if not lines:
                await self.main_window.dialog(
                    toga.ErrorDialog("Error", "The file is empty!")
                )
                return

            if len(lines) < 3:
                await self.main_window.dialog(
                    toga.ErrorDialog(
                        "Error",
                        "Invalid XYZ file: it must contain at least 3 lines "
                        "(number of atoms, comment line, and coordinates)."
                    )
                )
                return

            try:
                self.num_solvent_atoms = int(lines[0].strip())
            except ValueError:
                await self.main_window.dialog(
                    toga.ErrorDialog(
                        "Error",
                        "Invalid XYZ file: the first line must be an integer "
                        "with the number of atoms."
                    )
                )
                return

            coord_lines = [line for line in lines[2:] if line.strip()]

            if len(coord_lines) != self.num_solvent_atoms:
                await self.main_window.dialog(
                    toga.ErrorDialog(
                        "Error",
                        f"XYZ atom count mismatch: first line says "
                        f"{self.num_solvent_atoms}, but {len(coord_lines)} "
                        f"coordinate lines were found."
                    )
                )
                return

            for i, line in enumerate(coord_lines, start=3):
                tokens = line.split()

                if len(tokens) < 4:
                    await self.main_window.dialog(
                        toga.ErrorDialog(
                            "Error",
                            f"Invalid line {i} in solvent XYZ file:\n{line}"
                        )
                    )
                    return

                try:
                    atom = [
                        tokens[0],
                        float(tokens[1]),
                        float(tokens[2]),
                        float(tokens[3]),
                    ]
                except ValueError:
                    await self.main_window.dialog(
                        toga.ErrorDialog(
                            "Error",
                            f"Invalid numeric coordinates at line {i}:\n{line}"
                        )
                    )
                    return

                self.solvent_coords.append(atom)

            self.multi_line_text.value += (
                "\n--- Solvent molecule loaded successfully ---\n"
                f"File: {self.solvent_file}\n"
                f"Number of atoms: {self.num_solvent_atoms}\n\n"
                "Solvent Molecule:\n" + "".join(lines)
            )

        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Failed to read file: {e}")
            )
            return

    # async def read_solvent_xyz(self, widget):
    #     self.solvent_coords = []
    #     try:
    #         self.solvent_file = await self.main_window.dialog(toga.OpenFileDialog(
    #             title="Open Solvent XYZ file"
    #         ))
    #         if self.solvent_file:
    #             self.textInput_solvent.value = f"{self.solvent_file}"
    #         else:
    #             await self.main_window.dialog(toga.InfoDialog("Warning", "No file was selected!"))
    #             return
    #     except Exception as e:
    #         await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to open file: {e}"))
    #         return

    #     self.output_dir = os.path.dirname(self.solvent_file)

    #     try:
    #         with open(self.solvent_file, "r") as f:
    #             lines = f.readlines()
    #             if not lines:
    #                 await self.main_window.dialog(toga.ErrorDialog("Error", "The file is empty!"))
    #                 return
    #             self.num_solvent_atoms = int(lines[0].strip())
    #             for line in lines[2:]:
    #                 tokens = line.strip().split()
    #                 if len(tokens) >= 4:
    #                     atom = [tokens[0], float(tokens[1]), float(tokens[2]), float(tokens[3])]
    #                     self.solvent_coords.append(atom)
    #                 else:
    #                     await self.main_window.dialog(toga.ErrorDialog(
    #                         "Error", f"Invalid line in solvent XYZ file: {line}"
    #                     ))
    #                     return
    #             self.multi_line_text.value = "\nSolvent Molecule:\n" + "".join(lines)
    #     except Exception as e:
    #         await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to read file: {e}"))
    #         return

    async def read_solute_xyz(self, widget):
        self.solute_coords = []

        try:
            solute_file = await self.main_window.dialog(
                toga.OpenFileDialog(title="Open Solute XYZ file")
            )

            if not solute_file:
                await self.main_window.dialog(
                    toga.InfoDialog("Warning", "No file was selected!")
                )
                return

            self.solute_file = str(solute_file)
            self.textInput_solute.value = self.solute_file

            # Optional, but useful if the user loads the solute first
            if not getattr(self, "output_dir", ""):
                self.output_dir = os.path.dirname(self.solute_file)

        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Failed to open file: {e}")
            )
            return

        try:
            with open(self.solute_file, "r", encoding="utf-8") as f:
                lines = f.readlines()

            if not lines:
                await self.main_window.dialog(
                    toga.ErrorDialog("Error", "The file is empty!")
                )
                return

            if len(lines) < 3:
                await self.main_window.dialog(
                    toga.ErrorDialog(
                        "Error",
                        "Invalid XYZ file: it must contain at least 3 lines "
                        "(number of atoms, comment line, and coordinates)."
                    )
                )
                return

            try:
                self.num_solute_atoms = int(lines[0].strip())
            except ValueError:
                await self.main_window.dialog(
                    toga.ErrorDialog(
                        "Error",
                        "First line of XYZ must be an integer with the number of atoms."
                    )
                )
                return

            if self.num_solute_atoms <= 0:
                await self.main_window.dialog(
                    toga.ErrorDialog(
                        "Error",
                        "The number of atoms in the XYZ file must be greater than zero."
                    )
                )
                return

            coord_lines = [line for line in lines[2:] if line.strip()]

            if len(coord_lines) != self.num_solute_atoms:
                await self.main_window.dialog(
                    toga.ErrorDialog(
                        "Error",
                        f"XYZ atom count mismatch: first line says "
                        f"{self.num_solute_atoms}, but {len(coord_lines)} "
                        f"coordinate lines were found."
                    )
                )
                return

            for i, line in enumerate(coord_lines, start=3):
                tokens = line.split()

                if len(tokens) < 4:
                    await self.main_window.dialog(
                        toga.ErrorDialog(
                            "Error",
                            f"Invalid line {i} in solute XYZ file:\n{line}"
                        )
                    )
                    return

                symbol = tokens[0]

                try:
                    x = float(tokens[1])
                    y = float(tokens[2])
                    z = float(tokens[3])
                except ValueError:
                    await self.main_window.dialog(
                        toga.ErrorDialog(
                            "Error",
                            f"Invalid numeric coordinates at line {i}:\n{line}"
                        )
                    )
                    return

                self.solute_coords.append([symbol, x, y, z])

            if not self.solute_coords:
                await self.main_window.dialog(
                    toga.ErrorDialog(
                        "Error",
                        "No atom coordinates were found in the XYZ file."
                    )
                )
                return

            xs = [atom[1] for atom in self.solute_coords]
            ys = [atom[2] for atom in self.solute_coords]
            zs = [atom[3] for atom in self.solute_coords]

            xmin, xmax = min(xs), max(xs)
            ymin, ymax = min(ys), max(ys)
            zmin, zmax = min(zs), max(zs)

            dx = xmax - xmin
            dy = ymax - ymin
            dz = zmax - zmin

            self.solute_min = (xmin, ymin, zmin)
            self.solute_max = (xmax, ymax, zmax)
            self.solute_dims = (dx, dy, dz)
            self.solute_center = (
                0.5 * (xmin + xmax),
                0.5 * (ymin + ymax),
                0.5 * (zmin + zmax),
            )

            # Keep consistent with the rest of the class.
            # Your class uses self.spacing as the "extra wall padding" input. :contentReference[oaicite:1]{index=1}
            padding = max(0.0, getattr(self, "spacing", 0.0))
            self.orthorhombic_box = (
                dx + 2.0 * padding,
                dy + 2.0 * padding,
                dz + 2.0 * padding,
            )

            dims_str = (
                "\n--- Solute bounding box (Å) ---\n"
                f"File = {self.solute_file}\n"
                f"Number of atoms = {self.num_solute_atoms}\n"
                f"min = ({xmin:.4f}, {ymin:.4f}, {zmin:.4f})\n"
                f"max = ({xmax:.4f}, {ymax:.4f}, {zmax:.4f})\n"
                f"size (dx, dy, dz) = ({dx:.4f}, {dy:.4f}, {dz:.4f})\n"
                f"center = ({self.solute_center[0]:.4f}, "
                f"{self.solute_center[1]:.4f}, {self.solute_center[2]:.4f})\n"
                f"orthorhombic box with padding {padding:.2f} Å -> "
                f"a,b,c = ({self.orthorhombic_box[0]:.4f}, "
                f"{self.orthorhombic_box[1]:.4f}, "
                f"{self.orthorhombic_box[2]:.4f})\n"
            )

            self.multi_line_text.value += dims_str + "\nSolute Molecule:\n" + "".join(lines)

        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Failed to read file: {e}")
            )
            return


    # async def read_solute_xyz(self, widget):
    #     self.solute_coords = []
    #     try:
    #         self.solute_file = await self.main_window.dialog(toga.OpenFileDialog(
    #             title="Open Solute XYZ file"
    #         ))
    #         if self.solute_file:
    #             self.textInput_solute.value = f"{self.solute_file}"
    #         else:
    #             await self.main_window.dialog(toga.InfoDialog("Warning", "No file was selected!"))
    #             return
    #     except Exception as e:
    #         await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to open file: {e}"))
    #         return

    #     try:
    #         with open(self.solute_file, "r") as f:
    #             lines = f.readlines()
    #             if not lines:
    #                 await self.main_window.dialog(toga.ErrorDialog("Error", "The file is empty!"))
    #                 return
    #             try:
    #                 self.num_solute_atoms = int(lines[0].strip())
    #             except Exception:
    #                 await self.main_window.dialog(toga.ErrorDialog(
    #                     "Error", "First line of XYZ must be the number of atoms."
    #                 ))
    #                 return

    #             for line in lines[2:]:
    #                 tokens = line.strip().split()
    #                 if len(tokens) >= 4:
    #                     atom = [tokens[0], float(tokens[1]), float(tokens[2]), float(tokens[3])]
    #                     self.solute_coords.append(atom)
    #                 else:
    #                     await self.main_window.dialog(toga.ErrorDialog(
    #                         "Error", f"Invalid line in solute XYZ file: {line}"
    #                     ))
    #                     return

    #             if not self.solute_coords:
    #                 await self.main_window.dialog(toga.ErrorDialog(
    #                     "Error", "No atom coordinates were found in the XYZ file."
    #                 ))
    #                 return

    #             xs = [a[1] for a in self.solute_coords]
    #             ys = [a[2] for a in self.solute_coords]
    #             zs = [a[3] for a in self.solute_coords]
    #             xmin, xmax = min(xs), max(xs)
    #             ymin, ymax = min(ys), max(ys)
    #             zmin, zmax = min(zs), max(zs)
    #             dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
    #             self.solute_min = (xmin, ymin, zmin)
    #             self.solute_max = (xmax, ymax, zmax)
    #             self.solute_dims = (dx, dy, dz)
    #             self.solute_center = (0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax))

    #             padding = getattr(self, "box_padding", 1.0)
    #             self.orthorhombic_box = (dx + 2.0 * padding, dy + 2.0 * padding, dz + 2.0 * padding)

    #             dims_str = (
    #                 "\n--- Solute bounding box (Å) ---\n"
    #                 f"min = ({xmin:.4f}, {ymin:.4f}, {zmin:.4f})\n"
    #                 f"max = ({xmax:.4f}, {ymax:.4f}, {zmax:.4f})\n"
    #                 f"size (dx, dy, dz) = ({dx:.4f}, {dy:.4f}, {dz:.4f})\n"
    #                 f"center = ({self.solute_center[0]:.4f}, {self.solute_center[1]:.4f}, {self.solute_center[2]:.4f})\n"
    #                 f"orthorhombic box with padding {padding:.2f} Å → "
    #                 f"a,b,c = ({self.orthorhombic_box[0]:.4f}, {self.orthorhombic_box[1]:.4f}, {self.orthorhombic_box[2]:.4f})\n"
    #             )
    #             self.multi_line_text.value = dims_str + "\nSolute Molecule:\n" + "".join(lines)
    #     except Exception as e:
    #         await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to read file: {e}"))
    #         return

    def get_atomic_mass(self, symbol):
        return AtomicData.atomic_masses.get(symbol, 0.0)

    def get_vdw_radius(self, symbol):
        if hasattr(AtomicData, "atomic_vdw_radii"):
            radius = AtomicData.atomic_vdw_radii.get(symbol)
            if radius is not None:
                return float(radius)
        if hasattr(AtomicData, "vdw_radii"):
            radius = AtomicData.vdw_radii.get(symbol)
            if radius is not None:
                return float(radius)
        return self.DEFAULT_VDW_RADII.get(symbol, 1.70)

    def molecule_mass(self, coords):
        return sum(self.get_atomic_mass(atom[0]) for atom in coords)

    def center_molecule(self, coords):
        total_mass = 0.0
        center_x = center_y = center_z = 0.0
        for atom in coords:
            mass = self.get_atomic_mass(atom[0])
            total_mass += mass
            center_x += mass * atom[1]
            center_y += mass * atom[2]
            center_z += mass * atom[3]
        if total_mass <= 0.0:
            raise ValueError("Total molecular mass is zero. Check atomic symbols.")
        center_x /= total_mass
        center_y /= total_mass
        center_z /= total_mass
        centered = []
        for atom in coords:
            centered.append([atom[0], atom[1] - center_x, atom[2] - center_y, atom[3] - center_z])
        return centered

    def get_solvent_extent(self):
        max_extent = 0.0
        for atom in self.centered_solvent_coords:
            distance = math.sqrt(atom[1] ** 2 + atom[2] ** 2 + atom[3] ** 2)
            max_extent = max(max_extent, distance)
        return max_extent

    def get_bounding_radius(self, coords):
        radius = 0.0
        for atom in coords:
            atom_radius = math.sqrt(atom[1] ** 2 + atom[2] ** 2 + atom[3] ** 2) + self.get_vdw_radius(atom[0])
            radius = max(radius, atom_radius)
        return radius

    def center_of_mass(self):
        self.centered_solvent_coords = self.center_molecule(self.solvent_coords)
        self.solvent_extent = self.get_solvent_extent()
        self.solvent_mass_amu = self.molecule_mass(self.centered_solvent_coords)
        self.solvent_bounding_radius = self.get_bounding_radius(self.centered_solvent_coords)

        self.centered_solute_coords = []
        self.positioned_solute_coords = []
        self.solute_mass_amu = 0.0
        self.solute_bounding_radius = 0.0

        if self.include_solute and self.solute_coords:
            self.centered_solute_coords = self.center_molecule(self.solute_coords)
            self.solute_mass_amu = self.molecule_mass(self.centered_solute_coords)
            self.solute_bounding_radius = self.get_bounding_radius(self.centered_solute_coords)
            a, b, c = self.lattice_vectors
            self.positioned_solute_coords = self.translate_molecule(
                self.centered_solute_coords, np.array([a / 2.0, b / 2.0, c / 2.0], dtype=float)
            )

    def random_rotate_molecule(self, coordinates):
        atoms = [atom[0] for atom in coordinates]
        coords = np.array([atom[1:] for atom in coordinates], dtype=float)
        rotated = Rotation.random().apply(coords)
        return [[atoms[i], rotated[i][0], rotated[i][1], rotated[i][2]] for i in range(len(atoms))]

    def translate_molecule(self, coordinates, shift):
        sx, sy, sz = float(shift[0]), float(shift[1]), float(shift[2])
        return [[atom[0], atom[1] + sx, atom[2] + sy, atom[3] + sz] for atom in coordinates]

    def molecule_fits_in_box(self, coords):
        a, b, c = self.lattice_vectors
        wall_padding = max(0.0, self.spacing)
        for atom in coords:
            if not (wall_padding <= atom[1] <= a - wall_padding):
                return False
            if not (wall_padding <= atom[2] <= b - wall_padding):
                return False
            if not (wall_padding <= atom[3] <= c - wall_padding):
                return False
        return True

    def minimum_image_displacement(self, dx, dy, dz):
        """Apply the minimum-image convention for an orthorhombic periodic box."""
        if not getattr(self, "use_periodic_minimum_image", False):
            return dx, dy, dz
        if len(self.lattice_vectors) != 3:
            return dx, dy, dz

        a, b, c = self.lattice_vectors
        if a > 0.0:
            dx -= a * round(dx / a)
        if b > 0.0:
            dy -= b * round(dy / b)
        if c > 0.0:
            dz -= c * round(dz / c)
        return dx, dy, dz

    def pair_cutoff(self, symbol_i, symbol_j, extra=0.0):
        return self.vdw_scale * (self.get_vdw_radius(symbol_i) + self.get_vdw_radius(symbol_j)) + extra

    def molecules_clash(self, coords_a, coords_b, extra=0.0):
        for atom_a in coords_a:
            ax, ay, az = atom_a[1], atom_a[2], atom_a[3]
            for atom_b in coords_b:
                cutoff = self.pair_cutoff(atom_a[0], atom_b[0], extra=0.0)
                dx = ax - atom_b[1]
                dy = ay - atom_b[2]
                dz = az - atom_b[3]
                dx, dy, dz = self.minimum_image_displacement(dx, dy, dz)
                dist2 = dx * dx + dy * dy + dz * dz
                if dist2 < cutoff * cutoff:
                    return True
                if extra > 0.0:
                    expanded_cutoff = cutoff + extra
                    if dist2 < expanded_cutoff * expanded_cutoff:
                        return True
        return False

    def solvent_too_close_to_solute(self, solvent_coords):
        if not self.positioned_solute_coords:
            return False
        extra = max(0.0, self.min_distance)
        return self.molecules_clash(solvent_coords, self.positioned_solute_coords, extra=extra)

    def solvent_too_close_to_accepted(self, solvent_coords, accepted_molecules, candidate_center):
        candidate_radius = self.solvent_bounding_radius
        for existing_center, existing_coords in accepted_molecules:
            dx = candidate_center[0] - existing_center[0]
            dy = candidate_center[1] - existing_center[1]
            dz = candidate_center[2] - existing_center[2]
            dx, dy, dz = self.minimum_image_displacement(dx, dy, dz)
            center_dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if center_dist > (candidate_radius + self.solvent_bounding_radius + 2.5):
                continue
            if self.molecules_clash(solvent_coords, existing_coords, extra=0.0):
                return True
        return False

    def estimate_target_num_solvent_molecules(self):
        a, b, c = self.lattice_vectors
        volume_cm3 = (a * b * c) * 1e-24
        target_total_mass_grams = self.target_density * volume_cm3
        target_total_mass_amu = target_total_mass_grams / self.AMU_TO_GRAMS
        available_for_solvent_amu = target_total_mass_amu - self.solute_mass_amu
        if self.solvent_mass_amu <= 0.0:
            raise ValueError("Solvent molecular mass is zero. Check atomic symbols.")
        n_solvent = max(0, int(math.floor(available_for_solvent_amu / self.solvent_mass_amu)))
        return n_solvent

    def random_center_in_box(self):
        a, b, c = self.lattice_vectors
        margin = self.solvent_extent + max(0.0, self.spacing)
        if 2.0 * margin >= min(a, b, c):
            raise ValueError(
                "The solvent molecule is too large for the chosen box dimensions and wall padding."
            )
        x = random.uniform(margin, a - margin)
        y = random.uniform(margin, b - margin)
        z = random.uniform(margin, c - margin)
        return np.array([x, y, z], dtype=float)

    def reset_run_state(self):
        self.data = []
        self.solvent_insertion_count = 0
        self.solvent_attempts_used = 0
        self.failed_insertions = 0
        self.final_density = 0.0
        self.num_target_solvent_molecules = 0

    def insertion_progress_percent(self, attempts, accepted, target):
        """Return determinate progress for the random insertion loop.

        The loop can finish either by accepting the target number of solvent
        molecules or by exhausting max_attempts. The progress bar therefore
        tracks the closer completion condition so a successful early finish
        still reaches 100%.
        """
        attempt_fraction = attempts / max(1, self.max_attempts)
        accepted_fraction = 1.0 if target <= 0 else accepted / max(1, target)
        progress_fraction = max(attempt_fraction, accepted_fraction)
        return max(0.0, min(100.0, 100.0 * progress_fraction))

    def update_progress_display(self, attempts, accepted, target):
        progress_percent = self.insertion_progress_percent(attempts, accepted, target)
        if hasattr(self, "progress_bar"):
            # Keep the Toga widget in a stable 0-100 range; some backends do
            # not repaint reliably when max is changed during a long callback.
            self.progress_bar.max = 100
            self.progress_bar.value = progress_percent
        if hasattr(self, "status_label"):
            attempt_percent = 100.0 * attempts / max(1, self.max_attempts)
            accepted_percent = 100.0 if target <= 0 else 100.0 * accepted / max(1, target)
            self.status_label.text = (
                f"Insertion progress: {progress_percent:.1f}% | "
                f"attempts {attempts}/{self.max_attempts} ({attempt_percent:.1f}%) | "
                f"accepted {accepted}/{target} solvent molecules ({accepted_percent:.1f}%)"
            )

    # ------------------------------------------------------------------
    # Packing. Both loops below make the same accept/reject decision for
    # every attempt and consume the random generators identically, so a
    # fixed seed gives exactly the same box from either (pinned by
    # tests/test_single_solute_solvent_box.py). They run in a worker thread
    # and never touch a widget: progress goes through ``progress(dict)``.
    # ------------------------------------------------------------------
    def _insertion_payload(self, attempts, accepted, molecules):
        return {
            "stage": "insertion",
            "attempts": attempts,
            "max_attempts": self.max_attempts,
            "accepted": accepted,
            "target": self.num_target_solvent_molecules,
            "molecules": list(molecules),
        }

    def _solute_molecule(self):
        if not (self.include_solute and self.positioned_solute_coords):
            return None
        return (
            "solute",
            [atom[0] for atom in self.positioned_solute_coords],
            np.array([atom[1:] for atom in self.positioned_solute_coords], dtype=float),
        )

    def _finish_packing(self, box, molecules, accepted, attempts):
        self.solvent_insertion_count = accepted
        self.solvent_attempts_used = attempts
        self.failed_insertions = max(0, attempts - accepted)
        self.data = box
        self.molecules = molecules
        self.final_density = self.calculate_density()

    def fill_solvent_box_legacy(self, progress=None):
        """The original pure-Python insertion loop, kept as the reference
        implementation (the 'Legacy' engine, the tests and the benchmark)."""
        random.seed(self.random_seed)
        np.random.seed(self.random_seed if self.random_seed is not None else None)
        self.reset_run_state()

        box = []
        accepted_molecules = []
        molecules = []
        solute = self._solute_molecule()
        if self.positioned_solute_coords:
            box.extend([[atom[0], atom[1], atom[2], atom[3]] for atom in self.positioned_solute_coords])
        if solute is not None:
            molecules.append(solute)
        solvent_elements = [atom[0] for atom in self.centered_solvent_coords]

        self.num_target_solvent_molecules = self.estimate_target_num_solvent_molecules()
        accepted = 0
        attempts = 0
        self._report(progress, self._insertion_payload(attempts, accepted, molecules), force=True)

        while accepted < self.num_target_solvent_molecules and attempts < self.max_attempts:
            if self.cancel_requested:
                break
            attempts += 1
            if self.rotate_solvent:
                rotated_coords = self.random_rotate_molecule(self.centered_solvent_coords)
            else:
                rotated_coords = [atom[:] for atom in self.centered_solvent_coords]

            candidate_center = self.random_center_in_box()
            translated_coords = self.translate_molecule(rotated_coords, candidate_center)

            accepted_candidate = True
            if not self.molecule_fits_in_box(translated_coords):
                accepted_candidate = False
            elif self.include_solute and self.solvent_too_close_to_solute(translated_coords):
                accepted_candidate = False
            elif self.solvent_too_close_to_accepted(translated_coords, accepted_molecules, candidate_center):
                accepted_candidate = False

            if accepted_candidate:
                accepted_molecules.append((candidate_center, translated_coords))
                box.extend(translated_coords)
                molecules.append((
                    "solvent", solvent_elements,
                    np.array([atom[1:] for atom in translated_coords], dtype=float),
                ))
                accepted += 1

            self._report(progress, self._insertion_payload(attempts, accepted, molecules))

        self._finish_packing(box, molecules, accepted, attempts)
        self._report(progress, self._insertion_payload(attempts, accepted, molecules), force=True)

    def fill_solvent_box_fast(self, engine=ENGINE_NUMPY, progress=None):
        """Same decisions as :meth:`fill_solvent_box_legacy`, reached through
        a cell list and vectorised (or numba-compiled) clash tests."""
        random.seed(self.random_seed)
        np.random.seed(self.random_seed if self.random_seed is not None else None)
        self.reset_run_state()

        lengths = np.array(self.lattice_vectors, dtype=float)
        periodic = bool(self.use_periodic_minimum_image)
        solvent_elements = [atom[0] for atom in self.centered_solvent_coords]
        template = np.array([atom[1:] for atom in self.centered_solvent_coords], dtype=float)
        solvent_radii = np.array([self.get_vdw_radius(e) for e in solvent_elements], dtype=float)
        bounding = self.solvent_bounding_radius

        box = []
        molecules = []
        solute = self._solute_molecule()
        if self.positioned_solute_coords:
            box.extend([[atom[0], atom[1], atom[2], atom[3]] for atom in self.positioned_solute_coords])
        solute_xyz = solute_r = None
        if solute is not None:
            molecules.append(solute)
            solute_xyz = solute[2]
            solute_r = np.array([self.get_vdw_radius(e) for e in solute[1]], dtype=float)
        extra = max(0.0, self.min_distance or 0.0)

        wall = max(0.0, self.spacing)
        low = wall
        high = lengths - wall
        grid = MoleculeGrid(lengths, periodic, min_cell=2.0 * bounding + 2.5)
        rotations = RotationStream()

        self.num_target_solvent_molecules = self.estimate_target_num_solvent_molecules()
        accepted = 0
        attempts = 0
        self._report(progress, self._insertion_payload(attempts, accepted, molecules), force=True)

        while accepted < self.num_target_solvent_molecules and attempts < self.max_attempts:
            if self.cancel_requested:
                break
            attempts += 1
            if self.rotate_solvent:
                coords = template @ rotations.next().T
            else:
                coords = template.copy()
            center = self.random_center_in_box()
            xyz = coords + center

            ok = bool(np.all(xyz >= low) and np.all(xyz <= high))
            if ok and solute_xyz is not None:
                ok = not atoms_clash(xyz, solvent_radii, solute_xyz, solute_r, self.vdw_scale,
                                     lengths, periodic, extra, engine)
            if ok:
                ok = not grid.clashes(xyz, solvent_radii, center, bounding, self.vdw_scale, 0.0, engine)
            if ok:
                grid.add(center, bounding, xyz, solvent_radii)
                molecules.append(("solvent", solvent_elements, xyz))
                box.extend([[e, float(x), float(y), float(z)] for e, (x, y, z) in zip(solvent_elements, xyz)])
                accepted += 1

            self._report(progress, self._insertion_payload(attempts, accepted, molecules))

        self._finish_packing(box, molecules, accepted, attempts)
        self._report(progress, self._insertion_payload(attempts, accepted, molecules), force=True)

    def build_box(self, progress=None):
        """Pack the box with the selected engine (runs in a worker thread)."""
        engine = resolve_engine(self.packing_engine)
        started = time.perf_counter()
        if engine == ENGINE_LEGACY:
            self.fill_solvent_box_legacy(progress)
        else:
            self.fill_solvent_box_fast(engine, progress)
        self.build_seconds = time.perf_counter() - started
        self.engine_used = engine

    # ------------------------------------------------------------------
    # Hooks for the shared pre-optimization and properties (solvationBoxCommon)
    # ------------------------------------------------------------------
    def preopt_solute_clearance(self):
        return self.min_distance if self._solute_molecule() is not None else 0.0

    def net_charge(self):
        solute = self.solute_charge if self._solute_molecule() is not None else 0
        return int(solute) + int(self.solvent_insertion_count) * int(self.solvent_charge)

    def calculate_density(self):
        total_mass_amu = 0.0
        for atom in self.data:
            total_mass_amu += self.get_atomic_mass(atom[0])
        total_mass_grams = total_mass_amu * self.AMU_TO_GRAMS
        a, b, c = self.lattice_vectors
        volume_cm3 = (a * b * c) * 1e-24
        if volume_cm3 <= 0.0:
            return 0.0
        return total_mass_grams / volume_cm3

    def element_count(self):
        atom_list = [atom[0] for atom in self.data]
        num_atoms = len(atom_list)
        element_count = {}
        for element in atom_list:
            element_count[element] = element_count.get(element, 0) + 1
        with open(f"{self.output_dir}/single_solute_solvent.txt", "a") as f:
            f.write(f"Total number of atoms: {num_atoms}\n")
            for element, count in sorted(element_count.items()):
                f.write(f"{element}: {count}\n")
        return None

    def center_box(self):
        if not self.data:
            return None
        coords = [[atom[0], atom[1], atom[2], atom[3]] for atom in self.data]
        total_mass = 0.0
        center_x = center_y = center_z = 0.0
        for atom in coords:
            mass = self.get_atomic_mass(atom[0])
            total_mass += mass
            center_x += mass * atom[1]
            center_y += mass * atom[2]
            center_z += mass * atom[3]
        if total_mass <= 0.0:
            return None
        center_x /= total_mass
        center_y /= total_mass
        center_z /= total_mass
        for atom in coords:
            atom[1] -= center_x
            atom[2] -= center_y
            atom[3] -= center_z
        with open(f"{self.output_dir}/single_solute_box_cmass.xyz", "w") as f:
            f.write(f"{len(coords)}\n")
            f.write("Solvation box with 1 solute molecule plus solvent translated to the center of mass\n")
            for atom in coords:
                f.write(f"{atom[0]:<3s}{atom[1]:>14.7f}{atom[2]:>14.7f}{atom[3]:>14.7f}\n")
        return None

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    OUTPUT_LABEL = "single_solute_solvent_box"

    def output_paths(self):
        label = self.OUTPUT_LABEL
        return {
            "xyz": f"{self.output_dir}/{label}.xyz",
            "summary": f"{self.output_dir}/{label}.txt",
            "preopt": f"{self.output_dir}/{label}_preopt.xyz",
            "parameters": f"{self.output_dir}/single_solute_solvent.txt",
        }

    def build_summary_text(self):
        paths = self.output_paths()
        summary_lines = [
            "Random-insertion solvent box build completed.\n",
            f"Target solvent molecules: {self.num_target_solvent_molecules}\n",
            f"Inserted solvent molecules: {self.solvent_insertion_count}\n",
            f"Insertion attempts used: {self.solvent_attempts_used}\n",
            f"Final density: {self.final_density:.4f} g/cm³\n",
            f"Saved XYZ: {paths['xyz']}\n",
            f"Saved TXT summary: {paths['summary']}\n",
        ]
        if self.solvent_insertion_count < self.num_target_solvent_molecules:
            summary_lines.append(
                "Warning: the target number of solvent molecules was not fully reached. "
                "Try a larger box, a lower target density, a smaller vdW scale, or more attempts.\n"
            )
        if self.use_periodic_minimum_image:
            summary_lines.append("Periodic minimum-image clash detection: ON\n")
        else:
            summary_lines.append("Periodic minimum-image clash detection: OFF\n")
        if self.engine_used:
            summary_lines.append(f"Packing engine: {self.engine_used} ({self.build_seconds:.2f} s)\n")

        summary_lines.extend(self.extras_summary_lines(paths["preopt"]))
        return "".join(summary_lines)

    def write_summary(self):
        summary_text = self.build_summary_text()
        with open(self.output_paths()["summary"], "w", encoding="utf-8") as f:
            f.write(summary_text)
        return summary_text

    def write_parameters_file(self):
        with open(self.output_paths()["parameters"], "w") as f:
            f.write(f"Target density: {self.target_density:.4f} g/cm³\n")
            f.write(f"Final density: {self.final_density:.4f} g/cm³\n")
            f.write(f"Box lattice vectors:\n")
            f.write(f"  a: {self.lattice_vectors[0]:.3f} Å\n")
            f.write(f"  b: {self.lattice_vectors[1]:.3f} Å\n")
            f.write(f"  c: {self.lattice_vectors[2]:.3f} Å\n")
            f.write(f"Extra wall padding: {self.spacing:.3f} Å\n")
            f.write(f"vdW scale: {self.vdw_scale:.3f}\n")
            f.write(f"Target solvent molecules: {self.num_target_solvent_molecules}\n")
            f.write(f"Inserted solvent molecules: {self.solvent_insertion_count}\n")
            f.write(f"Insertion attempts used: {self.solvent_attempts_used}\n")
            f.write(f"Failed insertion attempts: {self.failed_insertions}\n")
            f.write(f"Periodic minimum-image clash detection: {'ON' if self.use_periodic_minimum_image else 'OFF'}\n")
            if self.include_solute:
                f.write(f"Extra solute clearance: {self.min_distance:.3f} Å\n")

    def write_outputs(self):
        """Every output file, in the original order; returns the summary text.

        File names and formats of the original outputs are unchanged; the
        pre-optimized structure goes to its own ``*_preopt.xyz``. Raises
        ``OSError``; the UI turns it into a dialog.
        """
        paths = self.output_paths()
        self.write_xyz(paths["xyz"], self.data,
                        "Solvent box generated by random-insertion solvent box builder")
        if self.preopt_result is not None:
            self.write_preopt_xyz(paths["preopt"])
        summary_text = self.write_summary()
        self.write_parameters_file()
        self.element_count()
        self.center_box()
        return summary_text

    async def save_and_display_results(self):
        try:
            summary_text = self.write_outputs()
        except OSError as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Error saving output files: {e}\n"))
            return False
        self.multi_line_text.value = summary_text
        return True


class SingleSoluteSolventBoxUI(SolvationExtrasUI, SingleSoluteSolventBox):
    """Setup / Construction / Summary tabs over a fixed run bar."""

    CHARGE_FIELDS = (
        ("Solute charge:", "Integer, default 0", "textInput_solute_charge"),
        ("Solvent molecule charge:", "Integer per molecule, default 0", "textInput_solvent_charge"),
    )

    def __init__(self, *args):
        super().__init__()
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        self.main_window = toga.Window(
            title="Single Solute Solvent Box Builder",
            size=(820, 720),
        )

        heading_style = Pack(font_size=18, font_weight="bold", text_align=LEFT, margin=(0, 0, 8, 0))
        self._section_style = Pack(font_size=12, font_weight="bold", text_align=LEFT, margin=(8, 0, 4, 0))
        self._label_style = Pack(margin=(0, 6, 0, 0), text_align=LEFT, width=170)
        self._input_style = Pack(flex=1, margin=0)
        self._button_style = Pack(margin=(0, 8, 0, 0), width=120)
        self._row_style = Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0))

        main_box = toga.Box(style=Pack(direction=COLUMN, margin=12, flex=1))
        main_box.add(toga.Label("Single Solute Solvent Box Builder", style=heading_style))

        # The status line and the summary box are created first: building the
        # tabs sets Selection/Switch values, whose handlers may report status.
        self.status_label = toga.Label("Idle", style=Pack(margin=(0, 0, 2, 0), text_align=LEFT))
        self.progress_bar = toga.ProgressBar(max=100, value=0, style=Pack(margin=(0, 0, 6, 0)))
        self.multi_line_text = toga.MultilineTextInput(
            readonly=True,
            style=Pack(flex=1, margin=0, font_size=11),
        )

        help_text = getattr(HelpGqteaWin, "single_solute_solvent_box", "Random insertion solvent box builder.")
        intro = (
            f"{help_text}\n\n"
            "Algorithm: target-density random insertion with random rotations and "
            "explicit solvent-solute / solvent-solvent vdW clash rejection, optionally followed by a "
            "periodic pre-optimization.\n"
        )
        main_box.add(self.build_tabs(self._build_setup_tab(), intro))
        main_box.add(self.build_run_bar("Build Solvent Box"))
        main_box.add(self.status_label)
        main_box.add(self.progress_bar)

        self.main_window.content = main_box
        self.update_solute_controls()
        self.update_preopt_controls()
        self.main_window.show()

    # -- layout helpers ----------------------------------------------------
    def _add_field(self, parent, label_text, placeholder, attr_name, button=None):
        row = toga.Box(style=self._row_style)
        label = toga.Label(label_text, style=self._label_style)
        text_input = toga.TextInput(placeholder=placeholder, style=self._input_style)
        setattr(self, attr_name, text_input)
        row.add(label)
        row.add(text_input)
        if button is not None:
            row.add(button)
        parent.add(row)

    def _build_setup_tab(self):
        setup = toga.Box(style=Pack(direction=COLUMN, margin=(8, 8, 8, 8)))
        section = self._section_style

        setup.add(toga.Label("System and packing parameters", style=section))
        params_row = toga.Box(style=Pack(direction=ROW, margin=(0, 0, 4, 0)))
        left_col = toga.Box(style=Pack(direction=COLUMN, flex=1, margin=(0, 8, 0, 0)))
        right_col = toga.Box(style=Pack(direction=COLUMN, flex=1))
        self._add_field(left_col, "Box lattice vectors:", "Enter a b c separated by spaces", "textInput_lattice_vectors")
        self._add_field(left_col, "Extra wall padding (Å):", "Optional padding from walls; 0.0 is fine", "textInput_spacing")
        self._add_field(left_col, "Target density (g/cm³):", "Example: 0.997", "textInput_target_density")
        self._add_field(right_col, "Max attempts:", "Example: 50000 or 100000", "textInput_max_attempts")
        self._add_field(right_col, "vdW scale:", "Example: 0.80", "textInput_vdw_scale")
        self._add_field(right_col, "Random seed:", "Optional integer", "textInput_seed")
        params_row.add(left_col)
        params_row.add(right_col)
        setup.add(params_row)

        setup.add(toga.Label("Packing options", style=section))
        options_row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0)))
        self.switch = toga.Switch("Randomly rotate solvent", style=Pack(margin=(0, 12, 0, 0)))
        self.switch.value = True
        self.switch_density = toga.Switch("Calculate density", style=Pack(margin=(0, 12, 0, 0)))
        self.switch_density.value = True
        self.switch_minimum_image = toga.Switch("Periodic minimum-image clash detection", style=Pack(margin=(0, 12, 0, 0)))
        self.switch_minimum_image.value = True
        self.switch_include_solute = toga.Switch(
            "Include solute", style=Pack(margin=(0, 0, 0, 0)), on_change=self.update_solute_controls
        )
        options_row.add(self.switch)
        options_row.add(self.switch_density)
        options_row.add(self.switch_minimum_image)
        options_row.add(self.switch_include_solute)
        setup.add(options_row)
        self._add_field(
            setup,
            "Extra solute clearance (Å):",
            "Additional solute-solvent clearance beyond vdW overlap",
            "textInput_min_distance",
        )

        setup.add(toga.Label("Input structures", style=section))
        self.btn_browse_solvent = toga.Button("Browse", on_press=self.read_solvent_xyz, style=self._button_style)
        self._add_field(setup, "Select Solvent:", "Click Browse to select solvent in XYZ format",
                        "textInput_solvent", button=self.btn_browse_solvent)
        self.btn_browse_solute = toga.Button("Browse", on_press=self.read_solute_xyz, style=self._button_style)
        self._add_field(setup, "Select Solute:", "Click Browse to select solute in XYZ format",
                        "textInput_solute", button=self.btn_browse_solute)

        self.build_charge_spin_section(setup, (
            "Net charge = solute charge + (inserted solvent molecules x solvent charge). The multiplicity "
            "is never guessed: it is checked against the electron count once the box is built."
        ))
        self.build_preopt_section(setup)
        self.build_engine_section(setup)
        return setup

    # -- inputs and control state -------------------------------------------
    def read_charges(self):
        self.solute_charge = parse_integer(self.textInput_solute_charge.value, "Solute charge", default=0)
        self.solvent_charge = parse_integer(
            self.textInput_solvent_charge.value, "Solvent molecule charge", default=0
        )

    def update_solute_controls(self, widget=None):
        enabled = bool(self.switch_include_solute.value) and not getattr(self, "is_running", False)
        self.textInput_min_distance.enabled = enabled
        self.textInput_solute.enabled = enabled
        self.btn_browse_solute.enabled = enabled
        self.textInput_solute_charge.enabled = enabled

    def set_running_state(self, is_running):
        self.set_extras_running(is_running)
        self.start_button.enabled = not is_running
        self.btn_close.enabled = not is_running
        self.btn_browse_solvent.enabled = not is_running
        for switch in [self.switch, self.switch_density, self.switch_minimum_image, self.switch_include_solute]:
            switch.enabled = not is_running
        for field in [
            self.textInput_lattice_vectors,
            self.textInput_spacing,
            self.textInput_target_density,
            self.textInput_max_attempts,
            self.textInput_vdw_scale,
            self.textInput_seed,
            self.textInput_solvent,
            self.textInput_solvent_charge,
        ]:
            field.enabled = not is_running
        self.update_solute_controls()
        if is_running:
            self.status_label.text = "Running insertion..."
            self.progress_bar.value = 0
            self.multi_line_text.value = "Preparing random insertion...\n"

    # -- run ---------------------------------------------------------------
    async def workflow(self, widget):
        if not await self.read_params(widget):
            return
        if not self.lattice_vectors:
            await self.main_window.dialog(toga.ErrorDialog("Error", "Please input valid box lattice vectors."))
            return

        if not self.solvent_coords:
            await self.main_window.dialog(toga.ErrorDialog("Error", "Please select a solvent molecule file."))
            return

        if self.include_solute and not self.solute_coords:
            await self.main_window.dialog(toga.ErrorDialog("Error", "Please select a solute molecule file."))
            return

        try:
            self.set_running_state(True)
            self.center_of_mass()
            initial = [m for m in [self._solute_molecule()] if m is not None]
            completed = await self.run_construction(
                initial,
                lambda: f"{self.solvent_insertion_count}/{self.num_target_solvent_molecules} molecules",
            )
            if not completed:
                return
            if not await self.save_and_display_results():
                self.status_label.text = "Failed"
                return

            if self.switch_density.value:
                self.multi_line_text.value += (
                    f"\nTarget density: {self.target_density:.4f} g/cm³\n"
                    f"Final density: {self.final_density:.4f} g/cm³\n"
                )

            self.status_label.text = "Completed"
            self.progress_bar.max = 100
            self.progress_bar.value = 100
            self.tabs.current_tab = self.SUMMARY_TAB
            await self.report_problems()
        except Exception as e:
            self.status_label.text = "Failed"
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to build solvent box: {e}"))
            return
        finally:
            self.set_running_state(False)

    def closeTopLevel(self, widget):
        self.main_window.close()
