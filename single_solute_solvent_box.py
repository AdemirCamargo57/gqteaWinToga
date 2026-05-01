import os
import math
import random
import numpy as np
from scipy.spatial.transform import Rotation
import toga
from toga.style import Pack
from toga.style.pack import LEFT, CENTER, ROW, COLUMN
import asyncio
from help import AtomicData, HelpGqteaWin


class SingleSoluteSolventBox:
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

    async def read_params(self, widget):
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
            return

        self.spacing = await read_input(
            self.textInput_spacing, "legacy spacing / extra padding", float, required=False, default=0.0
        )
        if self.spacing is None:
            return

        self.target_density = await read_input(
            self.textInput_target_density, "target density", float
        )
        if self.target_density is None or self.target_density <= 0.0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "Target density must be a positive number."
            ))
            return

        self.max_attempts = await read_input(
            self.textInput_max_attempts, "maximum insertion attempts", int, required=False, default=50000
        )
        if self.max_attempts is None or self.max_attempts <= 0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "Maximum insertion attempts must be a positive integer."
            ))
            return

        self.vdw_scale = await read_input(
            self.textInput_vdw_scale, "vdW scaling factor", float, required=False, default=0.80
        )
        if self.vdw_scale is None or self.vdw_scale <= 0.0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "vdW scaling factor must be positive."
            ))
            return

        seed_value = await read_input(
            self.textInput_seed, "random seed", int, required=False, default=None
        )
        self.random_seed = seed_value
        self.use_periodic_minimum_image = bool(self.switch_minimum_image.value)

        if self.switch_include_solute.value:
            self.min_distance = await read_input(
                self.textInput_min_distance, "minimum distance", float, required=False, default=0.0
            )
            if self.min_distance is None or self.min_distance < 0.0:
                await self.main_window.dialog(toga.ErrorDialog(
                    "Error", "Minimum distance must be zero or positive."
                ))
                return
        else:
            self.min_distance = 0.0

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
        if self.switch_include_solute.value:
            update_text += f"{'Extra solute clearance (Å):':<32} {self.min_distance:.3f}\n"
        self.multi_line_text.value += update_text

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

        if self.switch_include_solute.value and self.solute_coords:
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

    def update_progress_display(self, attempts, accepted, target):
        if hasattr(self, "progress_bar"):
            self.progress_bar.max = max(1, self.max_attempts)
            self.progress_bar.value = min(attempts, self.max_attempts)
        if hasattr(self, "status_label"):
            percent = 100.0 * attempts / max(1, self.max_attempts)
            self.status_label.text = (
                f"Insertion progress: {attempts}/{self.max_attempts} attempts "
                f"({percent:.1f}%) | accepted {accepted}/{target} solvent molecules"
            )

    async def fill_solvent_box_async(self):
        random.seed(self.random_seed)
        np.random.seed(self.random_seed if self.random_seed is not None else None)
        self.reset_run_state()

        box = []
        accepted_molecules = []

        if self.positioned_solute_coords:
            box.extend([[atom[0], atom[1], atom[2], atom[3]] for atom in self.positioned_solute_coords])

        self.num_target_solvent_molecules = self.estimate_target_num_solvent_molecules()
        accepted = 0
        attempts = 0
        update_every = max(100, min(2000, self.max_attempts // 100 if self.max_attempts > 0 else 100))

        self.update_progress_display(attempts, accepted, self.num_target_solvent_molecules)
        await asyncio.sleep(0)

        while accepted < self.num_target_solvent_molecules and attempts < self.max_attempts:
            attempts += 1
            if self.switch.value:
                rotated_coords = self.random_rotate_molecule(self.centered_solvent_coords)
            else:
                rotated_coords = [atom[:] for atom in self.centered_solvent_coords]

            candidate_center = self.random_center_in_box()
            translated_coords = self.translate_molecule(rotated_coords, candidate_center)

            accepted_candidate = True
            if not self.molecule_fits_in_box(translated_coords):
                accepted_candidate = False
            elif self.switch_include_solute.value and self.solvent_too_close_to_solute(translated_coords):
                accepted_candidate = False
            elif self.solvent_too_close_to_accepted(translated_coords, accepted_molecules, candidate_center):
                accepted_candidate = False

            if accepted_candidate:
                accepted_molecules.append((candidate_center, translated_coords))
                box.extend(translated_coords)
                accepted += 1

            if attempts % update_every == 0 or accepted == self.num_target_solvent_molecules or attempts == self.max_attempts:
                self.update_progress_display(attempts, accepted, self.num_target_solvent_molecules)
                if hasattr(self, "multi_line_text"):
                    self.multi_line_text.value = (
                        "Running random insertion...\n"
                        f"Target solvent molecules: {self.num_target_solvent_molecules}\n"
                        f"Accepted so far: {accepted}\n"
                        f"Attempts used: {attempts}\n"
                    )
                await asyncio.sleep(0)

        self.solvent_insertion_count = accepted
        self.solvent_attempts_used = attempts
        self.failed_insertions = max(0, attempts - accepted)
        self.data = box
        self.final_density = self.calculate_density()
        self.update_progress_display(attempts, accepted, self.num_target_solvent_molecules)

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

    async def save_and_display_results(self):
        label = "single_solute_solvent_box"
        xyz_path = f"{self.output_dir}/{label}.xyz"
        try:
            with open(xyz_path, "w") as f:
                f.write(f"{len(self.data)}\n")
                f.write("Solvent box generated by random-insertion solvent box builder\n")
                for atom in self.data:
                    f.write(f"{atom[0]:<3s}{atom[1]:>14.7f}{atom[2]:>14.7f}{atom[3]:>14.7f}\n")
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Error saving {label}: {e}\n"))
            return None

        summary_lines = [
            "Random-insertion solvent box build completed.\n",
            f"Target solvent molecules: {self.num_target_solvent_molecules}\n",
            f"Inserted solvent molecules: {self.solvent_insertion_count}\n",
            f"Insertion attempts used: {self.solvent_attempts_used}\n",
            f"Final density: {self.final_density:.4f} g/cm³\n",
            f"Saved XYZ: {xyz_path}\n",
            f"Saved TXT summary: {self.output_dir}/{label}.txt\n",
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

        summary_text = "".join(summary_lines)
        self.multi_line_text.value = summary_text

        try:
            summary_txt_path = f"{self.output_dir}/{label}.txt"
            with open(summary_txt_path, "w", encoding="utf-8") as f:
                f.write(summary_text)
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Error saving summary TXT file: {e}\n"))
            return None
        if hasattr(self, "status_label"):
            self.status_label.text = "Completed"
        if hasattr(self, "progress_bar"):
            self.progress_bar.value = self.max_attempts if self.solvent_insertion_count < self.num_target_solvent_molecules else min(self.solvent_attempts_used, self.max_attempts)
        return None


class SingleSoluteSolventBoxUI(SingleSoluteSolventBox):
    def __init__(self, *args):
        super().__init__()
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        self.main_window = toga.Window(
            title="Single Solute Solvent Box Builder",
            size=(760, 620),
        )

        heading_style = Pack(font_size=18, font_weight="bold", text_align=LEFT, margin=(0, 0, 8, 0))
        section_style = Pack(font_size=12, font_weight="bold", text_align=LEFT, margin=(8, 0, 4, 0))
        label_style = Pack(margin=(0, 6, 0, 0), text_align=LEFT, width=150)
        input_style = Pack(flex=1, margin=0)
        button_style = Pack(margin=(0, 8, 0, 0), width=120)
        row_style = Pack(direction=ROW, alignment=CENTER, margin=(0, 0, 4, 0))
        main_box = toga.Box(style=Pack(direction=COLUMN, margin=12))

        def add_field(parent, label_text, placeholder, attr_name, button=None):
            row = toga.Box(style=row_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            row.add(label)
            row.add(text_input)
            if button is not None:
                row.add(button)
            parent.add(row)

        title_label = toga.Label("Single Solute Solvent Box Builder", style=heading_style)
        main_box.add(title_label)

        main_box.add(toga.Label("System and packing parameters", style=section_style))
        params_row = toga.Box(style=Pack(direction=ROW, margin=(0, 0, 4, 0)))

        left_col = toga.Box(style=Pack(direction=COLUMN, flex=1, margin=(0, 8, 0, 0)))
        right_col = toga.Box(style=Pack(direction=COLUMN, flex=1))

        add_field(left_col, "Box lattice vectors:", "Enter a b c separated by spaces", "textInput_lattice_vectors")
        add_field(left_col, "Extra wall padding (Å):", "Optional padding from walls; 0.0 is fine", "textInput_spacing")
        add_field(left_col, "Target density (g/cm³):", "Example: 0.997", "textInput_target_density")

        add_field(right_col, "Max attempts:", "Example: 50000 or 100000", "textInput_max_attempts")
        add_field(right_col, "vdW scale:", "Example: 0.80", "textInput_vdw_scale")
        add_field(right_col, "Random seed:", "Optional integer", "textInput_seed")

        params_row.add(left_col)
        params_row.add(right_col)
        main_box.add(params_row)

        main_box.add(toga.Label("Packing options", style=section_style))
        options_row = toga.Box(style=Pack(direction=ROW, alignment=CENTER, margin=(0, 0, 4, 0)))
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
        main_box.add(options_row)

        add_field(
            main_box,
            "Extra solute clearance (Å):",
            "Additional solute-solvent clearance beyond vdW overlap",
            "textInput_min_distance",
        )

        main_box.add(toga.Label("Input structures", style=section_style))
        self.btn_browse_solvent = toga.Button("Browse", on_press=self.read_solvent_xyz, style=button_style)
        add_field(
            main_box,
            "Select Solvent:",
            "Click Browse to select solvent in XYZ format",
            "textInput_solvent",
            button=self.btn_browse_solvent,
        )

        self.btn_browse_solute = toga.Button("Browse", on_press=self.read_solute_xyz, style=button_style)
        add_field(
            main_box,
            "Select Solute:",
            "Click Browse to select solute in XYZ format",
            "textInput_solute",
            button=self.btn_browse_solute,
        )

        main_box.add(toga.Label("Run controls", style=section_style))
        run_box = toga.Box(style=Pack(direction=ROW, alignment=CENTER, margin=(0, 0, 6, 0)))
        self.start_button = toga.Button("Build Solvent Box", style=button_style, on_press=self.workflow)
        self.btn_close = toga.Button("Close", style=Pack(width=120), on_press=self.closeTopLevel)
        run_box.add(self.start_button)
        run_box.add(self.btn_close)
        main_box.add(run_box)

        main_box.add(toga.Label("Run status and summary", style=section_style))
        self.status_label = toga.Label("Idle", style=Pack(margin=(0, 0, 2, 0), text_align=LEFT))
        main_box.add(self.status_label)
        self.progress_bar = toga.ProgressBar(max=100, value=0, style=Pack(margin=(0, 0, 6, 0)))
        main_box.add(self.progress_bar)
        self.multi_line_text = toga.MultilineTextInput(
            readonly=True,
            style=Pack(flex=1, margin=0, font_size=11),
        )
        help_text = getattr(HelpGqteaWin, "single_solute_solvent_box", "Random insertion solvent box builder.")
        self.multi_line_text.value = (
            f"{help_text}\n\n"
            "Algorithm: target-density random insertion with random rotations and "
            "explicit solvent-solute / solvent-solvent vdW clash rejection.\n"
        )
        main_box.add(self.multi_line_text)

        self.main_window.content = main_box
        self.update_solute_controls()
        self.main_window.show()



    def update_solute_controls(self, widget=None):
        enabled = bool(self.switch_include_solute.value)
        self.textInput_min_distance.enabled = enabled
        self.textInput_solute.enabled = enabled
        self.btn_browse_solute.enabled = enabled

    def set_running_state(self, is_running):
        self.start_button.enabled = not is_running
        self.btn_close.enabled = not is_running
        self.btn_browse_solvent.enabled = not is_running
        self.btn_browse_solute.enabled = (not is_running) and bool(self.switch_include_solute.value)
        self.switch.enabled = not is_running
        self.switch_density.enabled = not is_running
        self.switch_minimum_image.enabled = not is_running
        self.switch_include_solute.enabled = not is_running

        for field in [
            self.textInput_lattice_vectors,
            self.textInput_spacing,
            self.textInput_target_density,
            self.textInput_max_attempts,
            self.textInput_vdw_scale,
            self.textInput_seed,
            self.textInput_min_distance,
            self.textInput_solvent,
            self.textInput_solute,
        ]:
            if field is self.textInput_min_distance:
                field.enabled = (not is_running) and bool(self.switch_include_solute.value)
            elif field is self.textInput_solute:
                field.enabled = (not is_running) and bool(self.switch_include_solute.value)
            else:
                field.enabled = not is_running

        if is_running:
            self.status_label.text = "Running insertion..."
            self.progress_bar.value = 0
            self.multi_line_text.value = "Preparing random insertion...\n"
        else:
            self.update_solute_controls()


    async def workflow(self, widget):
        await self.read_params(widget)
        if not self.lattice_vectors:
            await self.main_window.dialog(toga.ErrorDialog("Error", "Please input valid box lattice vectors."))
            return

        if not self.solvent_coords:
            await self.main_window.dialog(toga.ErrorDialog("Error", "Please select a solvent molecule file."))
            return

        if self.switch_include_solute.value and not self.solute_coords:
            await self.main_window.dialog(toga.ErrorDialog("Error", "Please select a solute molecule file."))
            return

        try:
            self.set_running_state(True)
            self.center_of_mass()
            await self.fill_solvent_box_async()
            await self.save_and_display_results()
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to build solvent box: {e}"))
            return
        finally:
            self.set_running_state(False)

        output_file = f"{self.output_dir}/single_solute_solvent.txt"
        with open(output_file, "w") as f:
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
            if self.switch_include_solute.value:
                f.write(f"Extra solute clearance: {self.min_distance:.3f} Å\n")

        if self.switch_density.value:
            self.multi_line_text.value += (
                f"\nTarget density: {self.target_density:.4f} g/cm³\n"
                f"Final density: {self.final_density:.4f} g/cm³\n"
            )

        self.element_count()
        self.center_box()

    def closeTopLevel(self, widget):
        self.main_window.close()
