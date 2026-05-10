
import os
import math
import random
import asyncio
import numpy as np
from scipy.spatial.transform import Rotation
import toga
from toga.style import Pack
from toga.style.pack import LEFT, CENTER, ROW, COLUMN
from help import AtomicData, HelpGqteaWin


class MixtureSolvationBox:
    """Build an orthorhombic solvent mixture box with an optional centered solute.

    Revised design:
    - density-driven random insertion instead of grid filling
    - optional solute centered in the box
    - explicit solvent-solvent and solvent-solute clash rejection
    - optional periodic minimum-image clash detection
    - compact Toga UI similar to the latest single-solute builder
    """

    AMU_TO_GRAMS = 1.66053906660e-24

    DEFAULT_VDW_RADII = {
        "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70,
        "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
        "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
        "K": 2.75, "Ca": 2.31, "Br": 1.85, "I": 1.98,
    }

    def __init__(self):
        self.coords_A = []
        self.coords_B = []
        self.solute_coords = []

        self.centered_A = []
        self.centered_B = []
        self.centered_solute = []
        self.positioned_solute_coords = []

        self.lattice_vectors = []
        self.spacing = 0.0  # optional wall padding
        self.output_dir = ""

        self.target_density = 1.0
        self.max_attempts = 50000
        self.vdw_scale = 0.80
        self.random_seed = None

        self.min_distance = 0.0  # extra solute clearance beyond vdW cutoff
        self.composition_A = 0.5
        self.composition_B = 0.5

        self.data = []  # flat list of atoms in final box
        self.accepted_molecules = []  # list of dicts with center/type/coords

        self.count_solvent_A = 0
        self.count_solvent_B = 0
        self.num_target_solvent_molecules = 0
        self.solvent_attempts_used = 0
        self.failed_insertions = 0
        self.final_density = 0.0

        self.mass_A_amu = 0.0
        self.mass_B_amu = 0.0
        self.solute_mass_amu = 0.0
        self.radius_A = 0.0
        self.radius_B = 0.0
        self.radius_solute = 0.0

        self.use_periodic_minimum_image = True

    # ---------- IO helpers ----------

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
            return False

        self.spacing = await read_input(
            self.textInput_spacing, "extra wall padding", float, required=False, default=0.0
        )
        if self.spacing is None or self.spacing < 0.0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "Extra wall padding must be zero or positive."
            ))
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
        self.random_seed = seed_value

        comp_A = await read_input(
            self.textInput_composition_A, "composition of solvent A (%)", float, required=False, default=50.0
        )
        comp_B = await read_input(
            self.textInput_composition_B, "composition of solvent B (%)", float, required=False, default=50.0
        )
        if comp_A is None or comp_B is None:
            return False
        if comp_A < 0.0 or comp_B < 0.0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "Solvent compositions must be zero or positive."
            ))
            return False
        if abs((comp_A + comp_B) - 100.0) > 1.0e-8:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "The compositions of solvent A and B must sum to 100%."
            ))
            return False

        self.composition_A = comp_A / 100.0
        self.composition_B = comp_B / 100.0
        self.use_periodic_minimum_image = bool(self.switch_minimum_image.value)

        if self.switch_include_solute.value:
            self.min_distance = await read_input(
                self.textInput_min_distance, "extra solute clearance", float, required=False, default=0.0
            )
            if self.min_distance is None or self.min_distance < 0.0:
                await self.main_window.dialog(toga.ErrorDialog(
                    "Error", "Extra solute clearance must be zero or positive."
                ))
                return False
        else:
            self.min_distance = 0.0

        update_text = (
            f"{'Box lattice vectors:':<36} {self.lattice_vectors}\n"
            f"{'Extra wall padding (Å):':<36} {self.spacing:.3f}\n"
            f"{'Target density (g/cm³):':<36} {self.target_density:.4f}\n"
            f"{'Max insertion attempts:':<36} {self.max_attempts}\n"
            f"{'vdW scale:':<36} {self.vdw_scale:.3f}\n"
            f"{'Composition solvent A (%):':<36} {100.0 * self.composition_A:.2f}\n"
            f"{'Composition solvent B (%):':<36} {100.0 * self.composition_B:.2f}\n"
            f"{'Periodic minimum-image clash:':<36} {'ON' if self.use_periodic_minimum_image else 'OFF'}\n"
        )
        if self.random_seed is not None:
            update_text += f"{'Random seed:':<36} {self.random_seed}\n"
        if self.switch_include_solute.value:
            update_text += f"{'Extra solute clearance (Å):':<36} {self.min_distance:.3f}\n"

        self.multi_line_text.value = update_text
        return True

    async def _read_xyz_file(self, file_title, target_attr, input_widget):
        coords = []
        try:
            file_path = await self.main_window.dialog(toga.OpenFileDialog(title=file_title))
            if file_path:
                input_widget.value = f"{file_path}"
            else:
                await self.main_window.dialog(toga.InfoDialog("Warning", "No file was selected!"))
                return None
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to open file: {e}"))
            return None

        if not self.output_dir:
            self.output_dir = os.path.dirname(file_path)

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                lines = f.readlines()
            if not lines:
                await self.main_window.dialog(toga.ErrorDialog("Error", "The file is empty!"))
                return None

            num_atoms = int(lines[0].strip())
            for line in lines[2:]:
                tokens = line.strip().split()
                if len(tokens) >= 4:
                    coords.append([tokens[0], float(tokens[1]), float(tokens[2]), float(tokens[3])])
                elif line.strip():
                    await self.main_window.dialog(toga.ErrorDialog(
                        "Error", f"Invalid line in XYZ file: {line}"
                    ))
                    return None

            setattr(self, target_attr, coords)
            self.multi_line_text.value = f"\n{file_title}:\n" + "".join(lines)
            return num_atoms
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to read file: {e}"))
            return None

    async def read_solvent_A_xyz(self, widget):
        await self._read_xyz_file("Open Solvent A XYZ file", "coords_A", self.textInput_solvent_A)

    async def read_solvent_B_xyz(self, widget):
        await self._read_xyz_file("Open Solvent B XYZ file", "coords_B", self.textInput_solvent_B)

    async def read_solute_xyz(self, widget):
        await self._read_xyz_file("Open Solute XYZ file", "solute_coords", self.textInput_solute)

    # ---------- geometry helpers ----------

    def get_atomic_mass(self, symbol):
        return AtomicData.atomic_masses.get(symbol, 0.0)

    def get_vdw_radius(self, symbol):
        if hasattr(AtomicData, "vdw_radii") and isinstance(AtomicData.vdw_radii, dict):
            return AtomicData.vdw_radii.get(symbol, self.DEFAULT_VDW_RADII.get(symbol, 1.70))
        return self.DEFAULT_VDW_RADII.get(symbol, 1.70)

    def molecule_mass(self, coords):
        return sum(self.get_atomic_mass(atom[0]) for atom in coords)

    def center_molecule(self, coords):
        if not coords:
            return []
        centered = [[atom[0], atom[1], atom[2], atom[3]] for atom in coords]
        total_mass = 0.0
        center_x = center_y = center_z = 0.0
        for atom in centered:
            mass = self.get_atomic_mass(atom[0])
            total_mass += mass
            center_x += mass * atom[1]
            center_y += mass * atom[2]
            center_z += mass * atom[3]
        if total_mass <= 0.0:
            arr = np.array([[a[1], a[2], a[3]] for a in centered], dtype=float)
            com = arr.mean(axis=0)
        else:
            com = np.array([center_x / total_mass, center_y / total_mass, center_z / total_mass], dtype=float)
        for atom in centered:
            atom[1] -= com[0]
            atom[2] -= com[1]
            atom[3] -= com[2]
        return centered

    def molecule_bounding_radius(self, coords):
        if not coords:
            return 0.0
        max_r = 0.0
        for atom in coords:
            r = math.sqrt(atom[1] ** 2 + atom[2] ** 2 + atom[3] ** 2)
            if r > max_r:
                max_r = r
        return max_r

    def random_rotate_molecule(self, coordinates):
        atoms = [atom[0] for atom in coordinates]
        coords = [atom[1:] for atom in coordinates]
        rotated = Rotation.random().apply(coords)
        return [[atoms[i], rotated[i][0], rotated[i][1], rotated[i][2]] for i in range(len(atoms))]

    def minimum_image_displacement(self, dx, dy, dz):
        if not self.use_periodic_minimum_image:
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

    def molecule_inside_box(self, coords):
        a, b, c = self.lattice_vectors
        wall_padding = self.spacing
        for atom in coords:
            if not (wall_padding <= atom[1] <= a - wall_padding):
                return False
            if not (wall_padding <= atom[2] <= b - wall_padding):
                return False
            if not (wall_padding <= atom[3] <= c - wall_padding):
                return False
        return True

    def prepare_molecules(self):
        self.centered_A = self.center_molecule(self.coords_A)
        self.centered_B = self.center_molecule(self.coords_B)
        self.mass_A_amu = self.molecule_mass(self.centered_A)
        self.mass_B_amu = self.molecule_mass(self.centered_B)
        self.radius_A = self.molecule_bounding_radius(self.centered_A)
        self.radius_B = self.molecule_bounding_radius(self.centered_B)

        if self.switch_include_solute.value and self.solute_coords:
            self.centered_solute = self.center_molecule(self.solute_coords)
            self.solute_mass_amu = self.molecule_mass(self.centered_solute)
            self.radius_solute = self.molecule_bounding_radius(self.centered_solute)
            a, b, c = self.lattice_vectors
            self.positioned_solute_coords = [
                [atom[0], atom[1] + a / 2.0, atom[2] + b / 2.0, atom[3] + c / 2.0]
                for atom in self.centered_solute
            ]
        else:
            self.centered_solute = []
            self.positioned_solute_coords = []
            self.solute_mass_amu = 0.0
            self.radius_solute = 0.0

    def estimate_target_molecule_count(self):
        a, b, c = self.lattice_vectors
        volume_cm3 = (a * b * c) * 1.0e-24
        total_mass_grams = self.target_density * volume_cm3
        total_mass_amu = total_mass_grams / self.AMU_TO_GRAMS

        avg_mass = self.composition_A * self.mass_A_amu + self.composition_B * self.mass_B_amu
        if avg_mass <= 0.0:
            self.num_target_solvent_molecules = 0
            return

        correction = 0.0
        if self.positioned_solute_coords:
            correction = self.solute_mass_amu

        estimate = int(max(0, round((total_mass_amu - correction) / avg_mass)))
        self.num_target_solvent_molecules = estimate

    def choose_solvent_type(self):
        total_placed = self.count_solvent_A + self.count_solvent_B
        if total_placed <= 0:
            return "A" if random.random() < self.composition_A else "B"

        current_fraction_A = self.count_solvent_A / total_placed if total_placed else 0.0
        current_fraction_B = self.count_solvent_B / total_placed if total_placed else 0.0

        need_A = max(0.0, self.composition_A - current_fraction_A)
        need_B = max(0.0, self.composition_B - current_fraction_B)
        norm = need_A + need_B
        if norm <= 1.0e-12:
            return "A" if random.random() < self.composition_A else "B"
        prob_A = need_A / norm
        return "A" if random.random() < prob_A else "B"

    def random_center_for_radius(self, radius):
        a, b, c = self.lattice_vectors
        wall_padding = max(self.spacing, radius)
        if wall_padding * 2 >= min(a, b, c):
            return None
        x = random.uniform(wall_padding, a - wall_padding)
        y = random.uniform(wall_padding, b - wall_padding)
        z = random.uniform(wall_padding, c - wall_padding)
        return np.array([x, y, z], dtype=float)

    def translate_molecule(self, coords, center):
        return [[atom[0], atom[1] + center[0], atom[2] + center[1], atom[3] + center[2]] for atom in coords]

    def solvent_too_close_to_accepted(self, solvent_coords, accepted_molecules, candidate_center, candidate_radius):
        for entry in accepted_molecules:
            existing_center = entry["center"]
            existing_radius = entry["radius"]
            dx = candidate_center[0] - existing_center[0]
            dy = candidate_center[1] - existing_center[1]
            dz = candidate_center[2] - existing_center[2]
            dx, dy, dz = self.minimum_image_displacement(dx, dy, dz)
            center_dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if center_dist > (candidate_radius + existing_radius + 2.5):
                continue
            if self.molecules_clash(solvent_coords, entry["coords"], extra=0.0):
                return True
        return False

    async def fill_solvent_box_async(self):
        self.accepted_molecules = []
        self.data = []
        self.count_solvent_A = 0
        self.count_solvent_B = 0

        if self.random_seed is not None:
            random.seed(self.random_seed)
            np.random.seed(self.random_seed)

        if self.positioned_solute_coords:
            self.data.extend([[atom[0], atom[1], atom[2], atom[3]] for atom in self.positioned_solute_coords])

        attempts = 0
        accepted = 0
        self.update_progress_display(attempts, accepted, self.num_target_solvent_molecules)
        await asyncio.sleep(0.001)

        while attempts < self.max_attempts and accepted < self.num_target_solvent_molecules:
            attempts += 1
            if attempts % 25 == 0:
                self.update_progress_display(attempts, accepted, self.num_target_solvent_molecules)
                await asyncio.sleep(0.001)

            solvent_type = self.choose_solvent_type()
            template = self.centered_A if solvent_type == "A" else self.centered_B
            radius = self.radius_A if solvent_type == "A" else self.radius_B

            if not template:
                continue

            rotated = self.random_rotate_molecule(template) if self.switch_rotate.value else [
                [atom[0], atom[1], atom[2], atom[3]] for atom in template
            ]

            center = self.random_center_for_radius(radius)
            if center is None:
                break

            translated = self.translate_molecule(rotated, center)

            if not self.molecule_inside_box(translated):
                continue

            if self.positioned_solute_coords:
                if self.molecules_clash(translated, self.positioned_solute_coords, extra=self.min_distance):
                    continue

            if self.solvent_too_close_to_accepted(translated, self.accepted_molecules, center, radius):
                continue

            self.accepted_molecules.append({
                "type": solvent_type,
                "center": center,
                "radius": radius,
                "coords": translated,
            })
            self.data.extend(translated)

            accepted += 1
            if solvent_type == "A":
                self.count_solvent_A += 1
            else:
                self.count_solvent_B += 1

            if attempts % 25 == 0 or accepted == self.num_target_solvent_molecules:
                self.update_progress_display(attempts, accepted, self.num_target_solvent_molecules)
                self.multi_line_text.value = (
                    "Running random insertion...\n"
                    f"Target solvent molecules: {self.num_target_solvent_molecules}\n"
                    f"Accepted so far: {accepted}\n"
                    f"Attempts used: {attempts}\n"
                    f"Solvent A molecules: {self.count_solvent_A}\n"
                    f"Solvent B molecules: {self.count_solvent_B}\n"
                )
                await asyncio.sleep(0.001)

        self.solvent_attempts_used = attempts
        self.failed_insertions = max(0, attempts - accepted)
        self.final_density = self.calculate_density()

    def calculate_density(self):
        total_mass_amu = 0.0
        for atom in self.data:
            total_mass_amu += self.get_atomic_mass(atom[0])
        total_mass_grams = total_mass_amu * self.AMU_TO_GRAMS
        a, b, c = self.lattice_vectors
        volume_cm3 = (a * b * c) * 1.0e-24
        if volume_cm3 <= 0.0:
            return 0.0
        return total_mass_grams / volume_cm3

    def update_progress_display(self, attempts, accepted, target):
        if hasattr(self, "status_label"):
            attempt_fraction = 0.0 if self.max_attempts <= 0 else min(1.0, attempts / self.max_attempts)
            self.status_label.text = (
                f"Insertion progress: {attempts}/{self.max_attempts} attempts ({100.0 * attempt_fraction:.1f}%) | "
                f"accepted {accepted}/{target} solvent molecules"
            )
        if hasattr(self, "progress_bar"):
            self.progress_bar.max = max(1, self.max_attempts)
            self.progress_bar.value = min(attempts, self.max_attempts)

    async def save_and_display_results(self):
        label = "mixture_solvent_box"
        xyz_path = f"{self.output_dir}/{label}.xyz"
        txt_path = f"{self.output_dir}/{label}.txt"

        try:
            with open(xyz_path, "w", encoding="utf-8") as f:
                f.write(f"{len(self.data)}\n")
                f.write("Mixture solvent box generated by density-driven random insertion\n")
                for atom in self.data:
                    f.write(f"{atom[0]:<3s}{atom[1]:>14.7f}{atom[2]:>14.7f}{atom[3]:>14.7f}\n")
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Error saving {label}: {e}\n"))
            return None

        summary_lines = [
            "Random-insertion solvent mixture box build completed.\n",
            f"Target solvent molecules: {self.num_target_solvent_molecules}\n",
            f"Inserted solvent A molecules: {self.count_solvent_A}\n",
            f"Inserted solvent B molecules: {self.count_solvent_B}\n",
            f"Insertion attempts used: {self.solvent_attempts_used}\n",
            f"Failed insertion attempts: {self.failed_insertions}\n",
            f"Final density: {self.final_density:.4f} g/cm³\n",
            f"Composition solvent A (%): {100.0 * self.composition_A:.2f}\n",
            f"Composition solvent B (%): {100.0 * self.composition_B:.2f}\n",
            f"Periodic minimum-image clash detection: {'ON' if self.use_periodic_minimum_image else 'OFF'}\n",
            f"Centered solute inserted: {'ON' if bool(self.positioned_solute_coords) else 'OFF'}\n",
            f"Saved XYZ: {xyz_path}\n",
            f"Saved TXT summary: {txt_path}\n",
            f"a = {self.lattice_vectors[0]:.3f} Å\n",
            f"b = {self.lattice_vectors[1]:.3f} Å\n",
            f"c = {self.lattice_vectors[2]:.3f} Å\n",
        ]
        if self.positioned_solute_coords:
            summary_lines.append(f"Extra solute clearance (Å): {self.min_distance:.3f}\n")
        if (self.count_solvent_A + self.count_solvent_B) < self.num_target_solvent_molecules:
            summary_lines.append(
                "Warning: the target number of solvent molecules was not fully reached. "
                "Try a larger box, a lower target density, a smaller vdW scale, or more attempts.\n"
            )

        summary_text = "".join(summary_lines)
        self.multi_line_text.value = summary_text

        try:
            with open(txt_path, "w", encoding="utf-8") as f:
                f.write(summary_text)
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Error saving summary TXT file: {e}\n"))
            return None

        if hasattr(self, "status_label"):
            self.status_label.text = "Completed"
        if hasattr(self, "progress_bar"):
            self.progress_bar.value = self.max_attempts if (self.count_solvent_A + self.count_solvent_B) < self.num_target_solvent_molecules else min(self.solvent_attempts_used, self.max_attempts)
        return None


class MixtureSolventBoxUI(MixtureSolvationBox):
    def __init__(self, *args):
        super().__init__()
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        self.main_window = toga.Window(
            title="Solvent Mixture Box Builder",
            size=(820, 650),
        )

        heading_style = Pack(font_size=18, font_weight="bold", text_align=LEFT, margin=(0, 0, 8, 0))
        section_style = Pack(font_size=12, font_weight="bold", text_align=LEFT, margin=(8, 0, 4, 0))
        label_style = Pack(margin=(0, 6, 0, 0), text_align=LEFT, width=170)
        input_style = Pack(flex=1, margin=0)
        button_style = Pack(margin=(0, 8, 0, 0), width=120)
        row_style = Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0))
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

        title_label = toga.Label("Solvent Mixture Box Builder", style=heading_style)
        main_box.add(title_label)

        main_box.add(toga.Label("System and packing parameters", style=section_style))
        params_row = toga.Box(style=Pack(direction=ROW, margin=(0, 0, 4, 0)))
        left_col = toga.Box(style=Pack(direction=COLUMN, flex=1, margin=(0, 8, 0, 0)))
        right_col = toga.Box(style=Pack(direction=COLUMN, flex=1))

        add_field(left_col, "Box lattice vectors:", "Enter a b c separated by spaces", "textInput_lattice_vectors")
        add_field(left_col, "Extra wall padding (Å):", "Optional padding from walls; 0.0 is fine", "textInput_spacing")
        add_field(left_col, "Target density (g/cm³):", "Example: 0.95", "textInput_target_density")

        add_field(right_col, "Max attempts:", "Example: 50000 or 100000", "textInput_max_attempts")
        add_field(right_col, "vdW scale:", "Example: 0.80", "textInput_vdw_scale")
        add_field(right_col, "Random seed:", "Optional integer", "textInput_seed")

        params_row.add(left_col)
        params_row.add(right_col)
        main_box.add(params_row)

        main_box.add(toga.Label("Mixture composition", style=section_style))
        comp_row = toga.Box(style=Pack(direction=ROW, margin=(0, 0, 4, 0)))
        comp_left = toga.Box(style=Pack(direction=COLUMN, flex=1, margin=(0, 8, 0, 0)))
        comp_right = toga.Box(style=Pack(direction=COLUMN, flex=1))
        add_field(comp_left, "Solvent A (%):", "Example: 70", "textInput_composition_A")
        add_field(comp_right, "Solvent B (%):", "Example: 30", "textInput_composition_B")
        comp_row.add(comp_left)
        comp_row.add(comp_right)
        main_box.add(comp_row)

        main_box.add(toga.Label("Packing options", style=section_style))
        options_row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0)))
        self.switch_rotate = toga.Switch("Randomly rotate solvent", style=Pack(margin=(0, 12, 0, 0)))
        self.switch_rotate.value = True
        self.switch_density = toga.Switch("Calculate density", style=Pack(margin=(0, 12, 0, 0)))
        self.switch_density.value = True
        self.switch_minimum_image = toga.Switch("Periodic minimum-image clash detection", style=Pack(margin=(0, 12, 0, 0)))
        self.switch_minimum_image.value = True
        self.switch_include_solute = toga.Switch("Insert centered solute", style=Pack(margin=0), on_change=self.update_solute_controls)
        self.switch_include_solute.value = False
        options_row.add(self.switch_rotate)
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
        self.btn_browse_solvent_A = toga.Button("Browse", on_press=self.read_solvent_A_xyz, style=button_style)
        add_field(
            main_box,
            "Select Solvent A:",
            "Click Browse to select solvent A in XYZ format",
            "textInput_solvent_A",
            button=self.btn_browse_solvent_A,
        )

        self.btn_browse_solvent_B = toga.Button("Browse", on_press=self.read_solvent_B_xyz, style=button_style)
        add_field(
            main_box,
            "Select Solvent B:",
            "Click Browse to select solvent B in XYZ format",
            "textInput_solvent_B",
            button=self.btn_browse_solvent_B,
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
        run_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 6, 0)))
        self.start_button = toga.Button("Build Mixture Box", style=button_style, on_press=self.workflow)
        self.btn_close = toga.Button("Close", style=Pack(width=120), on_press=self.closeTopLevel)
        run_box.add(self.start_button)
        run_box.add(self.btn_close)
        main_box.add(run_box)

        main_box.add(toga.Label("Run status and summary", style=section_style))
        self.status_label = toga.Label("Idle", style=Pack(margin=(0, 0, 2, 0), text_align=LEFT))
        main_box.add(self.status_label)
        self.progress_bar = toga.ProgressBar(max=100, value=0, style=Pack(margin=(0, 0, 6, 0)))
        main_box.add(self.progress_bar)
        self.multi_line_text = toga.MultilineTextInput(readonly=True, style=Pack(flex=1, margin=0, font_size=11))
        help_text = getattr(HelpGqteaWin, "mixture_solvent_box", "Random insertion solvent mixture box builder.")
        self.multi_line_text.value = (
            f"{help_text}\n\n"
            "Algorithm: target-density random insertion for solvent A/B, "
            "optional centered solute insertion, optional periodic minimum-image clash detection, "
            "and summary output in XYZ/TXT files.\n"
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
        self.btn_browse_solvent_A.enabled = not is_running
        self.btn_browse_solvent_B.enabled = not is_running
        self.btn_browse_solute.enabled = (not is_running) and bool(self.switch_include_solute.value)
        self.switch_rotate.enabled = not is_running
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
            self.textInput_composition_A,
            self.textInput_composition_B,
            self.textInput_min_distance,
            self.textInput_solvent_A,
            self.textInput_solvent_B,
            self.textInput_solute,
        ]:
            if field is self.textInput_min_distance:
                field.enabled = (not is_running) and bool(self.switch_include_solute.value)
            elif field is self.textInput_solute:
                field.enabled = (not is_running) and bool(self.switch_include_solute.value)
            else:
                field.enabled = not is_running

        if is_running:
            self.status_label.text = "Running insertion... please wait"
            self.progress_bar.value = 0
            self.multi_line_text.value = "Preparing random insertion...\n"
        else:
            self.update_solute_controls()

    async def workflow(self, widget):
        params_ok = await self.read_params(widget)
        if not params_ok:
            return

        if not self.coords_A:
            await self.main_window.dialog(toga.ErrorDialog("Error", "Please select solvent A molecule file."))
            return
        if not self.coords_B:
            await self.main_window.dialog(toga.ErrorDialog("Error", "Please select solvent B molecule file."))
            return
        if self.switch_include_solute.value and not self.solute_coords:
            await self.main_window.dialog(toga.ErrorDialog("Error", "Please select a solute molecule file."))
            return

        try:
            self.set_running_state(True)
            self.prepare_molecules()
            self.estimate_target_molecule_count()
            await self.fill_solvent_box_async()
            await self.save_and_display_results()
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to build solvent mixture box: {e}"))
            return
        finally:
            self.set_running_state(False)

    def closeTopLevel(self, widget):
        self.main_window.close()
