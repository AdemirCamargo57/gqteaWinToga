import os
import math
import random
import asyncio
import time
import numpy as np
from scipy.spatial.transform import Rotation
import toga
from toga.style import Pack
from toga.style.pack import LEFT, CENTER, ROW, COLUMN
from help import AtomicData, HelpGqteaWin


class SharedWallBox:
    """Build two adjacent solvent boxes sharing one z-normal face.

    Box 1 occupies z = 0..c1 and Box 2 occupies z = c1..c1+c2. The output XYZ
    contains only atoms; no wall or boundary marker is written at the interface.
    Periodic minimum-image clash checks use the full combined a x b x (c1+c2)
    lattice so the final structure is suitable for interface simulations.
    """

    AMU_TO_GRAMS = 1.66053906660e-24

    DEFAULT_VDW_RADII = {
        "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70,
        "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
        "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
        "K": 2.75, "Ca": 2.31, "Br": 1.85, "I": 1.98,
    }

    def __init__(self):
        self.lattice_vectors = []
        self.shared_z = 0.0
        self.output_dir = ""
        self.use_periodic_minimum_image = True
        self.combined_data = []
        self.boxes = [self.new_box_state(1), self.new_box_state(2)]
        self.cancel_requested = False
        self.is_running = False
        self.latest_progress = None
        self.last_applied_progress_key = None
        self.last_progress_text_update = 0.0
        self.last_progress_text_box = None
        self.last_text_update = 0.0

    def new_box_state(self, index):
        return {
            "index": index,
            "coords_A": [],
            "coords_B": [],
            "solute_coords": [],
            "centered_A": [],
            "centered_B": [],
            "centered_solute": [],
            "positioned_solute_coords": [],
            "spacing": 0.0,
            "target_density": 1.0,
            "max_attempts": 50000,
            "vdw_scale": 0.80,
            "random_seed": None,
            "min_distance": 0.0,
            "composition_A": 1.0,
            "composition_B": 0.0,
            "use_solute": False,
            "rotate": True,
            "data": [],
            "accepted_molecules": [],
            "count_solvent_A": 0,
            "count_solvent_B": 0,
            "num_target_solvent_molecules": 0,
            "solvent_attempts_used": 0,
            "failed_insertions": 0,
            "final_density": 0.0,
            "mass_A_amu": 0.0,
            "mass_B_amu": 0.0,
            "solute_mass_amu": 0.0,
            "radius_A": 0.0,
            "radius_B": 0.0,
            "radius_solute": 0.0,
        }

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
                    if len(labels) != 2:
                        await self.main_window.dialog(toga.InfoDialog(
                            "Info", "Please input exactly two shared lattice values: a b."
                        ))
                        return None
                    return labels
                return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(toga.ErrorDialog(
                    "Error", f"Invalid format for {field_name}: {e}"
                ))
                return None

        shared_ab = await read_input(self.textInput_shared_ab, "shared a b lattice values", list)
        if shared_ab is None:
            return False
        a, b = shared_ab
        c1 = await read_input(self.textInput_c1, "Box 1 height c1", float)
        c2 = await read_input(self.textInput_c2, "Box 2 height c2", float)
        if c1 is None or c2 is None:
            return False
        if a <= 0.0 or b <= 0.0 or c1 <= 0.0 or c2 <= 0.0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", "All lattice dimensions must be positive."
            ))
            return False

        self.shared_z = c1
        self.lattice_vectors = [a, b, c1 + c2]
        self.use_periodic_minimum_image = bool(self.switch_minimum_image.value)

        for box in self.boxes:
            ok = await self.read_box_params(box)
            if not ok:
                return False

        self.multi_line_text.value = self.build_parameter_summary()
        return True

    async def read_box_params(self, box):
        index = box["index"]

        async def read_input(attr_name, field_name, expected_type, required=True, default=None):
            text_input = getattr(self, attr_name)
            value = text_input.value.strip()
            if not value:
                if required:
                    await self.main_window.dialog(toga.ErrorDialog(
                        "Error", f"Please input a valid value for {field_name}."
                    ))
                    return None
                return default
            try:
                return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(toga.ErrorDialog(
                    "Error", f"Invalid format for {field_name}: {e}"
                ))
                return None

        box["spacing"] = await read_input(f"textInput_box{index}_spacing", f"Box {index} extra wall padding", float, required=False, default=0.0)
        box["target_density"] = await read_input(f"textInput_box{index}_target_density", f"Box {index} target density", float)
        box["max_attempts"] = await read_input(f"textInput_box{index}_max_attempts", f"Box {index} maximum insertion attempts", int, required=False, default=50000)
        box["vdw_scale"] = await read_input(f"textInput_box{index}_vdw_scale", f"Box {index} vdW scaling factor", float, required=False, default=0.80)
        box["random_seed"] = await read_input(f"textInput_box{index}_seed", f"Box {index} random seed", int, required=False, default=None)
        comp_A = await read_input(f"textInput_box{index}_composition_A", f"Box {index} solvent A composition", float, required=False, default=100.0)
        comp_B = await read_input(f"textInput_box{index}_composition_B", f"Box {index} solvent B composition", float, required=False, default=0.0)
        box["min_distance"] = await read_input(f"textInput_box{index}_min_distance", f"Box {index} extra solute clearance", float, required=False, default=0.0)

        if None in [box["spacing"], box["target_density"], box["max_attempts"], box["vdw_scale"], comp_A, comp_B, box["min_distance"]]:
            return False
        if box["spacing"] < 0.0 or box["target_density"] < 0.0 or box["max_attempts"] <= 0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", f"Box {index} padding, density, and max attempts are outside valid ranges."
            ))
            return False
        if box["vdw_scale"] <= 0.0 or box["min_distance"] < 0.0:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", f"Box {index} vdW scale must be positive and solute clearance cannot be negative."
            ))
            return False
        if comp_A < 0.0 or comp_B < 0.0 or abs((comp_A + comp_B) - 100.0) > 1.0e-8:
            await self.main_window.dialog(toga.ErrorDialog(
                "Error", f"Box {index} solvent A and B compositions must be non-negative and sum to 100%."
            ))
            return False

        box["composition_A"] = comp_A / 100.0
        box["composition_B"] = comp_B / 100.0
        box["rotate"] = bool(getattr(self, f"switch_box{index}_rotate").value)
        box["use_solute"] = bool(getattr(self, f"switch_box{index}_include_solute").value)
        if not box["use_solute"]:
            box["min_distance"] = 0.0
        return True

    async def _read_xyz_file(self, file_title, box_index, target_key, input_widget):
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

            self.boxes[box_index - 1][target_key] = coords
            self.multi_line_text.value = f"\n{file_title}:\n" + "".join(lines)
            return num_atoms
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to read file: {e}"))
            return None

    async def read_box1_solvent_A_xyz(self, widget):
        await self._read_xyz_file("Open Box 1 Solvent A XYZ file", 1, "coords_A", self.textInput_box1_solvent_A)

    async def read_box1_solvent_B_xyz(self, widget):
        await self._read_xyz_file("Open Box 1 Solvent B XYZ file", 1, "coords_B", self.textInput_box1_solvent_B)

    async def read_box1_solute_xyz(self, widget):
        await self._read_xyz_file("Open Box 1 Solute XYZ file", 1, "solute_coords", self.textInput_box1_solute)

    async def read_box2_solvent_A_xyz(self, widget):
        await self._read_xyz_file("Open Box 2 Solvent A XYZ file", 2, "coords_A", self.textInput_box2_solvent_A)

    async def read_box2_solvent_B_xyz(self, widget):
        await self._read_xyz_file("Open Box 2 Solvent B XYZ file", 2, "coords_B", self.textInput_box2_solvent_B)

    async def read_box2_solute_xyz(self, widget):
        await self._read_xyz_file("Open Box 2 Solute XYZ file", 2, "solute_coords", self.textInput_box2_solute)

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

    def pair_cutoff(self, symbol_i, symbol_j, vdw_scale, extra=0.0):
        return vdw_scale * (self.get_vdw_radius(symbol_i) + self.get_vdw_radius(symbol_j)) + extra

    def molecules_clash(self, coords_a, coords_b, vdw_scale, extra=0.0):
        for atom_a in coords_a:
            ax, ay, az = atom_a[1], atom_a[2], atom_a[3]
            for atom_b in coords_b:
                cutoff = self.pair_cutoff(atom_a[0], atom_b[0], vdw_scale, extra=extra)
                dx = ax - atom_b[1]
                dy = ay - atom_b[2]
                dz = az - atom_b[3]
                dx, dy, dz = self.minimum_image_displacement(dx, dy, dz)
                dist2 = dx * dx + dy * dy + dz * dz
                if dist2 < cutoff * cutoff:
                    return True
        return False

    def translate_molecule(self, coords, center):
        return [[atom[0], atom[1] + center[0], atom[2] + center[1], atom[3] + center[2]] for atom in coords]

    def slab_bounds(self, box):
        if box["index"] == 1:
            return 0.0, self.shared_z
        return self.shared_z, self.lattice_vectors[2]

    def random_center_for_radius(self, box, radius):
        a, b, total_c = self.lattice_vectors
        z_min, z_max = self.slab_bounds(box)
        xy_padding = max(box["spacing"], radius)
        z_low_padding = box["spacing"] if z_min == 0.0 else 0.0
        z_high_padding = box["spacing"] if z_max == total_c else 0.0

        if xy_padding * 2.0 >= min(a, b):
            return None
        if z_low_padding + z_high_padding >= (z_max - z_min):
            return None

        x = random.uniform(xy_padding, a - xy_padding)
        y = random.uniform(xy_padding, b - xy_padding)
        z = random.uniform(z_min + z_low_padding, z_max - z_high_padding)
        return np.array([x, y, z], dtype=float)

    def molecule_inside_slab(self, coords, box):
        a, b, total_c = self.lattice_vectors
        z_min, z_max = self.slab_bounds(box)
        x_min = box["spacing"]
        x_max = a - box["spacing"]
        y_min = box["spacing"]
        y_max = b - box["spacing"]
        z_lower = z_min + (box["spacing"] if z_min == 0.0 else 0.0)
        z_upper = z_max - (box["spacing"] if z_max == total_c else 0.0)
        for atom in coords:
            if not (x_min <= atom[1] <= x_max):
                return False
            if not (y_min <= atom[2] <= y_max):
                return False
            if not (z_lower <= atom[3] <= z_upper):
                return False
        return True

    def solvent_too_close_to_entries(self, solvent_coords, entries, candidate_center, candidate_radius, vdw_scale):
        for entry in entries:
            existing_center = entry["center"]
            existing_radius = entry["radius"]
            dx = candidate_center[0] - existing_center[0]
            dy = candidate_center[1] - existing_center[1]
            dz = candidate_center[2] - existing_center[2]
            dx, dy, dz = self.minimum_image_displacement(dx, dy, dz)
            center_dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if center_dist > (candidate_radius + existing_radius + 2.5):
                continue
            if self.molecules_clash(solvent_coords, entry["coords"], vdw_scale, extra=0.0):
                return True
        return False

    def grid_cell_size(self):
        radii = []
        for box in self.boxes:
            for key in ["radius_A", "radius_B", "radius_solute"]:
                if box.get(key, 0.0) > 0.0:
                    radii.append(box[key])
        max_radius = max(radii) if radii else 2.0
        return max(4.0, min(12.0, max_radius + 2.5))

    def make_entry_grid(self, entries, cell_size):
        a, b, c = self.lattice_vectors
        shape = (
            max(1, int(math.ceil(a / cell_size))),
            max(1, int(math.ceil(b / cell_size))),
            max(1, int(math.ceil(c / cell_size))),
        )
        grid = {
            "cells": {},
            "shape": shape,
            "widths": (a / shape[0], b / shape[1], c / shape[2]),
            "max_radius": 0.0,
        }
        for entry in entries:
            self.add_entry_to_grid(grid, entry)
        return grid

    def entry_cell_key(self, center, grid):
        shape = grid["shape"]
        key = []
        for axis, length in enumerate(self.lattice_vectors):
            value = float(center[axis])
            if self.use_periodic_minimum_image and length > 0.0:
                value = value % length
            index = int(math.floor(value / grid["widths"][axis]))
            index = max(0, min(shape[axis] - 1, index))
            key.append(index)
        return tuple(key)

    def add_entry_to_grid(self, grid, entry):
        key = self.entry_cell_key(entry["center"], grid)
        grid["cells"].setdefault(key, []).append(entry)
        grid["max_radius"] = max(grid["max_radius"], entry["radius"])

    def iter_nearby_entries(self, grid, center, search_radius):
        if not grid["cells"]:
            return []
        base_key = self.entry_cell_key(center, grid)
        ranges = [
            int(math.ceil(search_radius / width)) if width > 0.0 else 0
            for width in grid["widths"]
        ]
        nearby = []
        seen = set()
        shape = grid["shape"]
        for dx in range(-ranges[0], ranges[0] + 1):
            ix = base_key[0] + dx
            if self.use_periodic_minimum_image:
                ix %= shape[0]
            elif ix < 0 or ix >= shape[0]:
                continue
            for dy in range(-ranges[1], ranges[1] + 1):
                iy = base_key[1] + dy
                if self.use_periodic_minimum_image:
                    iy %= shape[1]
                elif iy < 0 or iy >= shape[1]:
                    continue
                for dz in range(-ranges[2], ranges[2] + 1):
                    iz = base_key[2] + dz
                    if self.use_periodic_minimum_image:
                        iz %= shape[2]
                    elif iz < 0 or iz >= shape[2]:
                        continue
                    for entry in grid["cells"].get((ix, iy, iz), []):
                        entry_id = id(entry)
                        if entry_id in seen:
                            continue
                        seen.add(entry_id)
                        nearby.append(entry)
        return nearby

    def solvent_too_close_to_grid(self, solvent_coords, grid, candidate_center, candidate_radius, vdw_scale):
        if not grid["cells"]:
            return False
        search_radius = candidate_radius + grid["max_radius"] + 2.5
        nearby_entries = self.iter_nearby_entries(grid, candidate_center, search_radius)
        return self.solvent_too_close_to_entries(
            solvent_coords, nearby_entries, candidate_center, candidate_radius, vdw_scale
        )

    def make_entry(self, kind, center, radius, coords, box_index):
        return {
            "kind": kind,
            "box": box_index,
            "center": np.array(center, dtype=float),
            "radius": radius,
            "coords": coords,
        }

    # ---------- packing helpers ----------

    def prepare_box_molecules(self, box):
        box["centered_A"] = self.center_molecule(box["coords_A"])
        box["centered_B"] = self.center_molecule(box["coords_B"])
        box["mass_A_amu"] = self.molecule_mass(box["centered_A"])
        box["mass_B_amu"] = self.molecule_mass(box["centered_B"])
        box["radius_A"] = self.molecule_bounding_radius(box["centered_A"])
        box["radius_B"] = self.molecule_bounding_radius(box["centered_B"])

        if box["use_solute"] and box["solute_coords"]:
            box["centered_solute"] = self.center_molecule(box["solute_coords"])
            box["solute_mass_amu"] = self.molecule_mass(box["centered_solute"])
            box["radius_solute"] = self.molecule_bounding_radius(box["centered_solute"])
            a, b, _ = self.lattice_vectors
            z_min, z_max = self.slab_bounds(box)
            center = np.array([a / 2.0, b / 2.0, (z_min + z_max) / 2.0], dtype=float)
            box["positioned_solute_coords"] = self.translate_molecule(box["centered_solute"], center)
        else:
            box["centered_solute"] = []
            box["positioned_solute_coords"] = []
            box["solute_mass_amu"] = 0.0
            box["radius_solute"] = 0.0

    def estimate_target_molecule_count(self, box):
        a, b, _ = self.lattice_vectors
        z_min, z_max = self.slab_bounds(box)
        volume_cm3 = (a * b * (z_max - z_min)) * 1.0e-24
        total_mass_grams = box["target_density"] * volume_cm3
        total_mass_amu = total_mass_grams / self.AMU_TO_GRAMS

        avg_mass = box["composition_A"] * box["mass_A_amu"] + box["composition_B"] * box["mass_B_amu"]
        if avg_mass <= 0.0:
            box["num_target_solvent_molecules"] = 0
            return

        correction = box["solute_mass_amu"] if box["positioned_solute_coords"] else 0.0
        box["num_target_solvent_molecules"] = int(max(0, round((total_mass_amu - correction) / avg_mass)))

    def choose_solvent_type(self, box):
        if box["composition_A"] <= 0.0:
            return "B"
        if box["composition_B"] <= 0.0:
            return "A"

        total_placed = box["count_solvent_A"] + box["count_solvent_B"]
        if total_placed <= 0:
            return "A" if random.random() < box["composition_A"] else "B"

        current_fraction_A = box["count_solvent_A"] / total_placed
        current_fraction_B = box["count_solvent_B"] / total_placed
        need_A = max(0.0, box["composition_A"] - current_fraction_A)
        need_B = max(0.0, box["composition_B"] - current_fraction_B)
        norm = need_A + need_B
        if norm <= 1.0e-12:
            return "A" if random.random() < box["composition_A"] else "B"
        return "A" if random.random() < (need_A / norm) else "B"

    def calculate_box_density(self, box):
        total_mass_amu = sum(self.get_atomic_mass(atom[0]) for atom in box["data"])
        total_mass_grams = total_mass_amu * self.AMU_TO_GRAMS
        a, b, _ = self.lattice_vectors
        z_min, z_max = self.slab_bounds(box)
        volume_cm3 = (a * b * (z_max - z_min)) * 1.0e-24
        if volume_cm3 <= 0.0:
            return 0.0
        return total_mass_grams / volume_cm3

    def local_progress_percent(self, attempts, accepted, target, max_attempts):
        attempt_fraction = attempts / max(1, max_attempts)
        accepted_fraction = 1.0 if target <= 0 else accepted / max(1, target)
        return 100.0 * max(attempt_fraction, accepted_fraction)

    def progress_update_interval(self, box):
        """Throttle UI work during packing.

        The worker publishes progress by attempt count, including long
        rejection-heavy stretches where no molecules are accepted. The Toga UI
        still applies only changed snapshots, so this can stay responsive
        without repainting every candidate.
        """
        return max(100, min(500, box["max_attempts"] // 200 if box["max_attempts"] > 0 else 100))

    def progress_report_interval(self, box):
        """Report progress every 10% of max insertion attempts."""
        return max(1, math.ceil(box["max_attempts"] / 10))

    def event_yield_interval(self, box):
        """Yield often enough for Cancel/Close clicks without repainting every time."""
        return max(50, min(250, box["max_attempts"] // 500 if box["max_attempts"] > 0 else 50))

    def update_progress_display(self, box, attempts, accepted, base_percent, span_percent):
        local_percent = max(0.0, min(100.0, self.local_progress_percent(
            attempts, accepted, box["num_target_solvent_molecules"], box["max_attempts"]
        )))
        total_percent = base_percent + span_percent * local_percent / 100.0
        if hasattr(self, "progress_bar"):
            self.progress_bar.max = 100
            self.progress_bar.value = max(0.0, min(100.0, total_percent))
        if hasattr(self, "status_label"):
            self.status_label.text = (
                f"Packing Box {box['index']}: {local_percent:.1f}% | "
                f"attempts {attempts}/{box['max_attempts']} | "
                f"accepted {accepted}/{box['num_target_solvent_molecules']}"
            )

    def publish_progress(self, box, attempts, accepted, base_percent, span_percent):
        local_percent = max(0.0, min(100.0, self.local_progress_percent(
            attempts, accepted, box["num_target_solvent_molecules"], box["max_attempts"]
        )))
        total_percent = base_percent + span_percent * local_percent / 100.0
        self.latest_progress = {
            "box_index": box["index"],
            "attempts": attempts,
            "accepted": accepted,
            "target": box["num_target_solvent_molecules"],
            "max_attempts": box["max_attempts"],
            "local_percent": local_percent,
            "total_percent": max(0.0, min(100.0, total_percent)),
            "count_A": box["count_solvent_A"],
            "count_B": box["count_solvent_B"],
        }

    def apply_latest_progress(self, force_text=False):
        progress = self.latest_progress
        if not progress:
            return
        progress_key = (progress["box_index"], progress["attempts"], progress["accepted"])
        if progress_key == self.last_applied_progress_key:
            return
        self.last_applied_progress_key = progress_key
        self.progress_bar.max = 100
        self.progress_bar.value = progress["total_percent"]
        self.status_label.text = (
            f"Packing Box {progress['box_index']}: {progress['local_percent']:.1f}% | "
            f"attempts {progress['attempts']}/{progress['max_attempts']} | "
            f"accepted {progress['accepted']}/{progress['target']}"
        )
        now = time.monotonic()
        box_changed = progress["box_index"] != self.last_progress_text_box
        should_update_text = force_text or box_changed or (now - self.last_progress_text_update >= 0.75)
        if not should_update_text:
            return
        self.last_progress_text_update = now
        self.last_progress_text_box = progress["box_index"]
        self.multi_line_text.value = (
            "Running shared-wall random insertion...\n"
            "The packing backend is running outside the UI thread.\n"
            f"Current box: {progress['box_index']}\n"
            f"Target solvent molecules: {progress['target']}\n"
            f"Accepted so far: {progress['accepted']}\n"
            f"Attempts used: {progress['attempts']}\n"
            f"Solvent A molecules: {progress['count_A']}\n"
            f"Solvent B molecules: {progress['count_B']}\n"
        )

    def fill_box_sync(self, box, existing_entries, base_percent, span_percent):
        box["data"] = []
        box["accepted_molecules"] = []
        box["count_solvent_A"] = 0
        box["count_solvent_B"] = 0

        if box["random_seed"] is not None:
            random.seed(box["random_seed"])
            np.random.seed(box["random_seed"])

        solute_entry = None
        if box["positioned_solute_coords"]:
            box["data"].extend([[atom[0], atom[1], atom[2], atom[3]] for atom in box["positioned_solute_coords"]])
            z_min, z_max = self.slab_bounds(box)
            center = [self.lattice_vectors[0] / 2.0, self.lattice_vectors[1] / 2.0, (z_min + z_max) / 2.0]
            solute_entry = self.make_entry("solute", center, box["radius_solute"], box["positioned_solute_coords"], box["index"])

        attempts = 0
        accepted = 0
        update_every = self.progress_update_interval(box)
        cell_size = self.grid_cell_size()
        existing_grid = self.make_entry_grid(existing_entries, cell_size)
        accepted_grid = self.make_entry_grid([], cell_size)
        self.publish_progress(box, attempts, accepted, base_percent, span_percent)

        while attempts < box["max_attempts"] and accepted < box["num_target_solvent_molecules"]:
            if self.cancel_requested:
                break
            attempts += 1
            if attempts % update_every == 0:
                self.publish_progress(box, attempts, accepted, base_percent, span_percent)

            solvent_type = self.choose_solvent_type(box)
            template = box["centered_A"] if solvent_type == "A" else box["centered_B"]
            radius = box["radius_A"] if solvent_type == "A" else box["radius_B"]
            if not template:
                continue

            rotated = self.random_rotate_molecule(template) if box["rotate"] else [
                [atom[0], atom[1], atom[2], atom[3]] for atom in template
            ]
            center = self.random_center_for_radius(box, radius)
            if center is None:
                break

            translated = self.translate_molecule(rotated, center)
            if not self.molecule_inside_slab(translated, box):
                continue
            if box["positioned_solute_coords"]:
                if self.molecules_clash(translated, box["positioned_solute_coords"], box["vdw_scale"], extra=box["min_distance"]):
                    continue
            if self.solvent_too_close_to_grid(translated, existing_grid, center, radius, box["vdw_scale"]):
                continue
            if self.solvent_too_close_to_grid(translated, accepted_grid, center, radius, box["vdw_scale"]):
                continue

            entry = self.make_entry(f"solvent_{solvent_type}", center, radius, translated, box["index"])
            box["accepted_molecules"].append(entry)
            self.add_entry_to_grid(accepted_grid, entry)
            box["data"].extend(translated)
            accepted += 1
            if solvent_type == "A":
                box["count_solvent_A"] += 1
            else:
                box["count_solvent_B"] += 1

            self.publish_progress(box, attempts, accepted, base_percent, span_percent)

        box["solvent_attempts_used"] = attempts
        box["failed_insertions"] = max(0, attempts - accepted)
        box["final_density"] = self.calculate_box_density(box)
        self.publish_progress(box, attempts, accepted, base_percent, span_percent)

        entries = []
        if solute_entry is not None:
            entries.append(solute_entry)
        entries.extend(box["accepted_molecules"])
        return entries

    def build_shared_wall_box_sync(self):
        self.combined_data = []
        global_entries = []
        for box in self.boxes:
            self.prepare_box_molecules(box)
            self.estimate_target_molecule_count(box)

        box1_entries = self.fill_box_sync(self.boxes[0], global_entries, 0.0, 50.0)
        if self.cancel_requested:
            return False
        global_entries.extend(box1_entries)
        box2_entries = self.fill_box_sync(self.boxes[1], global_entries, 50.0, 50.0)
        if self.cancel_requested:
            return False
        global_entries.extend(box2_entries)

        self.combined_data = []
        for box in self.boxes:
            self.combined_data.extend(box["data"])

        if self.combined_data:
            min_z = min(atom[3] for atom in self.combined_data)
            for atom in self.combined_data:
                atom[3] -= min_z

        return True

    async def run_packing_worker(self):
        completed = await self.build_shared_wall_box_async()
        return completed

    async def fill_box_async(self, box, existing_entries, base_percent, span_percent):
        box["data"] = []
        box["accepted_molecules"] = []
        box["count_solvent_A"] = 0
        box["count_solvent_B"] = 0

        if box["random_seed"] is not None:
            random.seed(box["random_seed"])
            np.random.seed(box["random_seed"])

        solute_entry = None
        if box["positioned_solute_coords"]:
            box["data"].extend([[atom[0], atom[1], atom[2], atom[3]] for atom in box["positioned_solute_coords"]])
            z_min, z_max = self.slab_bounds(box)
            center = [self.lattice_vectors[0] / 2.0, self.lattice_vectors[1] / 2.0, (z_min + z_max) / 2.0]
            solute_entry = self.make_entry("solute", center, box["radius_solute"], box["positioned_solute_coords"], box["index"])

        attempts = 0
        accepted = 0
        report_every = self.progress_report_interval(box)
        yield_every = self.event_yield_interval(box)
        self.update_progress_display(box, attempts, accepted, base_percent, span_percent)
        await asyncio.sleep(0)

        while attempts < box["max_attempts"] and accepted < box["num_target_solvent_molecules"]:
            if self.cancel_requested:
                break
            attempts += 1
            if attempts % yield_every == 0:
                await asyncio.sleep(0)
                if self.cancel_requested:
                    break
            if attempts % report_every == 0:
                self.update_progress_display(box, attempts, accepted, base_percent, span_percent)
                await asyncio.sleep(0)
                if self.cancel_requested:
                    break

            solvent_type = self.choose_solvent_type(box)
            template = box["centered_A"] if solvent_type == "A" else box["centered_B"]
            radius = box["radius_A"] if solvent_type == "A" else box["radius_B"]
            if not template:
                continue

            rotated = self.random_rotate_molecule(template) if box["rotate"] else [
                [atom[0], atom[1], atom[2], atom[3]] for atom in template
            ]
            center = self.random_center_for_radius(box, radius)
            if center is None:
                break

            translated = self.translate_molecule(rotated, center)
            if not self.molecule_inside_slab(translated, box):
                continue
            if box["positioned_solute_coords"]:
                if self.molecules_clash(translated, box["positioned_solute_coords"], box["vdw_scale"], extra=box["min_distance"]):
                    continue

            if self.solvent_too_close_to_entries(translated, existing_entries, center, radius, box["vdw_scale"]):
                continue
            if self.solvent_too_close_to_entries(translated, box["accepted_molecules"], center, radius, box["vdw_scale"]):
                continue

            entry = self.make_entry(f"solvent_{solvent_type}", center, radius, translated, box["index"])
            box["accepted_molecules"].append(entry)
            box["data"].extend(translated)
            accepted += 1
            if solvent_type == "A":
                box["count_solvent_A"] += 1
            else:
                box["count_solvent_B"] += 1

            if attempts % report_every == 0 or accepted == box["num_target_solvent_molecules"]:
                self.update_progress_display(box, attempts, accepted, base_percent, span_percent)
                now = time.monotonic()
                if now - self.last_text_update >= 0.75:
                    self.last_text_update = now
                    self.multi_line_text.value = (
                        "Running shared-wall random insertion...\n"
                        "Box 1 and Box 2 are adjacent along z; no wall atoms are written.\n"
                        "Periodic minimum-image clash checks use the full combined box.\n"
                        f"Current box: {box['index']}\n"
                        f"Target solvent molecules: {box['num_target_solvent_molecules']}\n"
                        f"Accepted so far: {accepted}\n"
                        f"Attempts used: {attempts}\n"
                        f"Solvent A molecules: {box['count_solvent_A']}\n"
                        f"Solvent B molecules: {box['count_solvent_B']}\n"
                    )
                await asyncio.sleep(0)

        box["solvent_attempts_used"] = attempts
        box["failed_insertions"] = max(0, attempts - accepted)
        box["final_density"] = self.calculate_box_density(box)
        self.update_progress_display(box, attempts, accepted, base_percent, span_percent)
        if self.cancel_requested and hasattr(self, "status_label"):
            self.status_label.text = "Cancelling after current packing step..."

        entries = []
        if solute_entry is not None:
            entries.append(solute_entry)
        entries.extend(box["accepted_molecules"])
        return entries

    async def build_shared_wall_box_async(self):
        self.combined_data = []
        global_entries = []
        for box in self.boxes:
            self.prepare_box_molecules(box)
            self.estimate_target_molecule_count(box)

        box1_entries = await self.fill_box_async(self.boxes[0], global_entries, 0.0, 50.0)
        if self.cancel_requested:
            return False
        global_entries.extend(box1_entries)
        box2_entries = await self.fill_box_async(self.boxes[1], global_entries, 50.0, 50.0)
        if self.cancel_requested:
            return False
        global_entries.extend(box2_entries)

        self.combined_data = []
        for box in self.boxes:
            self.combined_data.extend(box["data"])

        if self.combined_data:
            min_z = min(atom[3] for atom in self.combined_data)
            for atom in self.combined_data:
                atom[3] -= min_z

        return True

    # ---------- validation and output ----------

    async def validate_selected_files(self):
        for box in self.boxes:
            index = box["index"]
            if box["target_density"] > 0.0:
                if box["composition_A"] > 0.0 and not box["coords_A"]:
                    await self.main_window.dialog(toga.ErrorDialog("Error", f"Please select Box {index} solvent A molecule file."))
                    return False
                if box["composition_B"] > 0.0 and not box["coords_B"]:
                    await self.main_window.dialog(toga.ErrorDialog("Error", f"Please select Box {index} solvent B molecule file."))
                    return False
            if box["use_solute"] and not box["solute_coords"]:
                await self.main_window.dialog(toga.ErrorDialog("Error", f"Please select Box {index} solute molecule file."))
                return False
        return True

    def build_parameter_summary(self):
        lines = [
            "Shared-wall solvent box parameters\n",
            f"Combined lattice vectors: a={self.lattice_vectors[0]:.3f}, b={self.lattice_vectors[1]:.3f}, c={self.lattice_vectors[2]:.3f} Angstrom\n",
            f"Shared interface at z={self.shared_z:.3f} Angstrom\n",
            f"Periodic minimum-image clash detection: {'ON' if self.use_periodic_minimum_image else 'OFF'}\n",
        ]
        for box in self.boxes:
            z_min, z_max = self.slab_bounds(box)
            lines.extend([
                f"\nBox {box['index']} z range: {z_min:.3f} to {z_max:.3f} Angstrom\n",
                f"Target density (g/cm^3): {box['target_density']:.4f}\n",
                f"Max insertion attempts: {box['max_attempts']}\n",
                f"vdW scale: {box['vdw_scale']:.3f}\n",
                f"Composition solvent A (%): {100.0 * box['composition_A']:.2f}\n",
                f"Composition solvent B (%): {100.0 * box['composition_B']:.2f}\n",
                f"Centered solute: {'ON' if box['use_solute'] else 'OFF'}\n",
            ])
        return "".join(lines)

    def build_result_summary(self, xyz_path, txt_path):
        lines = [
            "Shared-wall solvent box build completed.\n",
            "The XYZ file contains two adjacent z-stacked regions with no wall atoms at the interface.\n",
            f"Saved XYZ: {xyz_path}\n",
            f"Saved TXT summary: {txt_path}\n",
            f"Combined lattice: a={self.lattice_vectors[0]:.3f}, b={self.lattice_vectors[1]:.3f}, c={self.lattice_vectors[2]:.3f} Angstrom\n",
            f"Shared interface z: {self.shared_z:.3f} Angstrom\n",
            f"Periodic minimum-image clash detection: {'ON' if self.use_periodic_minimum_image else 'OFF'}\n",
        ]
        for box in self.boxes:
            z_min, z_max = self.slab_bounds(box)
            total_inserted = box["count_solvent_A"] + box["count_solvent_B"]
            lines.extend([
                f"\nBox {box['index']} summary\n",
                f"z range: {z_min:.3f} to {z_max:.3f} Angstrom\n",
                f"Target solvent molecules: {box['num_target_solvent_molecules']}\n",
                f"Inserted solvent A molecules: {box['count_solvent_A']}\n",
                f"Inserted solvent B molecules: {box['count_solvent_B']}\n",
                f"Insertion attempts used: {box['solvent_attempts_used']}\n",
                f"Failed insertion attempts: {box['failed_insertions']}\n",
                f"Final density: {box['final_density']:.4f} g/cm^3\n",
                f"Composition solvent A (%): {100.0 * box['composition_A']:.2f}\n",
                f"Composition solvent B (%): {100.0 * box['composition_B']:.2f}\n",
                f"Centered solute inserted: {'ON' if bool(box['positioned_solute_coords']) else 'OFF'}\n",
            ])
            if box["positioned_solute_coords"]:
                lines.append(f"Extra solute clearance (Angstrom): {box['min_distance']:.3f}\n")
            if total_inserted < box["num_target_solvent_molecules"]:
                lines.append(
                    "Warning: the target number of solvent molecules was not fully reached. "
                    "Try a larger slab, lower density, smaller vdW scale, or more attempts.\n"
                )
        return "".join(lines)

    async def save_and_display_results(self):
        label = "shared_wall_box"
        xyz_path = f"{self.output_dir}/{label}.xyz"
        txt_path = f"{self.output_dir}/{label}.txt"

        try:
            with open(xyz_path, "w", encoding="utf-8") as f:
                f.write(f"{len(self.combined_data)}\n")
                f.write("Shared-wall double solvent box generated by density-driven random insertion\n")
                for atom in self.combined_data:
                    f.write(f"{atom[0]:<3s}{atom[1]:>14.7f}{atom[2]:>14.7f}{atom[3]:>14.7f}\n")
        except Exception as e:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Error saving {label}: {e}\n"))
            return None

        summary_text = self.build_result_summary(xyz_path, txt_path)
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
            self.progress_bar.max = 100
            self.progress_bar.value = 100
        return None


class SharedWallBoxUI(SharedWallBox):
    def __init__(self, *args):
        super().__init__()
        self._controls = []
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        self.main_window = toga.Window(
            title="Shared-Wall Double Solvent Box Builder",
            size=(760, 720),
        )

        heading_style = Pack(font_size=18, font_weight="bold", text_align=LEFT, margin=(0, 0, 8, 0))
        section_style = Pack(font_size=12, font_weight="bold", text_align=LEFT, margin=(8, 0, 4, 0))
        label_style = Pack(margin=(0, 6, 0, 0), text_align=LEFT, width=150)
        input_style = Pack(flex=1, margin=0)
        button_style = Pack(margin=(0, 8, 0, 0), width=120)
        row_style = Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0))
        content_box = toga.Box(style=Pack(direction=COLUMN, margin=12))

        def add_field(parent, label_text, placeholder, attr_name, button=None):
            row = toga.Box(style=row_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            self._controls.append(text_input)
            row.add(label)
            row.add(text_input)
            if button is not None:
                row.add(button)
            parent.add(row)

        title_label = toga.Label("Shared-Wall Double Box Builder", style=heading_style)
        content_box.add(title_label)

        content_box.add(toga.Label("Shared geometry", style=section_style))
        geometry_row = toga.Box(style=Pack(direction=ROW, margin=(0, 0, 4, 0)))
        geom_left = toga.Box(style=Pack(direction=COLUMN, flex=1, margin=(0, 8, 0, 0)))
        geom_right = toga.Box(style=Pack(direction=COLUMN, flex=1))
        add_field(geom_left, "Shared a b:", "a b", "textInput_shared_ab")
        add_field(geom_left, "Box 1 c:", "c1", "textInput_c1")
        add_field(geom_right, "Box 2 c:", "c2", "textInput_c2")
        self.switch_minimum_image = toga.Switch("Periodic minimum-image clash detection", style=Pack(margin=(0, 8, 0, 0)))
        self.switch_minimum_image.value = True
        self._controls.append(self.switch_minimum_image)
        geom_right.add(self.switch_minimum_image)
        geometry_row.add(geom_left)
        geometry_row.add(geom_right)
        content_box.add(geometry_row)

        self.add_box_controls(content_box, 1, "Box 1 lower z region", button_style, section_style, row_style, label_style, input_style)
        self.add_box_controls(content_box, 2, "Box 2 upper z region", button_style, section_style, row_style, label_style, input_style)

        content_box.add(toga.Label("Run controls", style=section_style))
        run_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 6, 0)))
        self.start_button = toga.Button("Build Shared-Wall Box", style=Pack(margin=(0, 8, 0, 0), width=180), on_press=self.workflow)
        self.btn_cancel = toga.Button("Cancel", style=button_style, on_press=self.cancel_run)
        self.btn_cancel.enabled = False
        self.btn_close = toga.Button("Close", style=Pack(width=120), on_press=self.closeTopLevel)
        self._controls.append(self.start_button)
        run_box.add(self.start_button)
        run_box.add(self.btn_cancel)
        run_box.add(self.btn_close)
        content_box.add(run_box)

        content_box.add(toga.Label("Run status and summary", style=section_style))
        self.progress_bar = toga.ProgressBar(max=100, value=0, style=Pack(margin=(0, 0, 6, 0)))
        content_box.add(self.progress_bar)
        self.multi_line_text = toga.MultilineTextInput(readonly=True, style=Pack(height=180, margin=0, font_size=11))
        help_text = getattr(HelpGqteaWin, "shared_wall_box", "Build two z-stacked solvent regions that share one face.")
        self.multi_line_text.value = (
            f"{help_text}\n\n"
            "Algorithm: density-driven random insertion in Box 1 and Box 2, "
            "separate solvent compositions and solutes per box, shared z interface, "
            "and full combined-box periodic minimum-image clash detection.\n"
        )
        content_box.add(self.multi_line_text)

        self.main_window.content = toga.ScrollContainer(content=content_box, style=Pack(flex=1))
        self.update_solute_controls()
        self.update_solvent_controls()
        self.main_window.show()

    def add_box_controls(self, parent, index, title, button_style, section_style, row_style, label_style, input_style):
        def add_field(local_parent, label_text, placeholder, attr_name, button=None, on_change=None):
            row = toga.Box(style=row_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(placeholder=placeholder, style=input_style, on_change=on_change)
            setattr(self, attr_name, text_input)
            self._controls.append(text_input)
            row.add(label)
            row.add(text_input)
            if button is not None:
                row.add(button)
            local_parent.add(row)

        parent.add(toga.Label(title, style=section_style))
        params_row = toga.Box(style=Pack(direction=ROW, margin=(0, 0, 4, 0)))
        left_col = toga.Box(style=Pack(direction=COLUMN, flex=1, margin=(0, 8, 0, 0)))
        right_col = toga.Box(style=Pack(direction=COLUMN, flex=1))

        add_field(left_col, "Wall padding:", "0.0", f"textInput_box{index}_spacing")
        add_field(left_col, "Density:", "0.95", f"textInput_box{index}_target_density", on_change=self.update_solvent_controls)
        add_field(left_col, "Solvent A (%):", "Example: 100", f"textInput_box{index}_composition_A")
        add_field(left_col, "Solvent B (%):", "Example: 0", f"textInput_box{index}_composition_B")
        add_field(right_col, "Max attempts:", "50000", f"textInput_box{index}_max_attempts")
        add_field(right_col, "vdW scale:", "Example: 0.80", f"textInput_box{index}_vdw_scale")
        add_field(right_col, "Random seed:", "Optional integer", f"textInput_box{index}_seed")
        add_field(right_col, "Solute gap:", "0.0", f"textInput_box{index}_min_distance")
        params_row.add(left_col)
        params_row.add(right_col)
        parent.add(params_row)

        options_row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0)))
        rotate_switch = toga.Switch(f"Rotate Box {index}", style=Pack(margin=(0, 12, 0, 0)))
        rotate_switch.value = True
        solute_switch = toga.Switch(f"Box {index} solute", style=Pack(margin=0), on_change=self.update_solute_controls)
        solute_switch.value = False
        setattr(self, f"switch_box{index}_rotate", rotate_switch)
        setattr(self, f"switch_box{index}_include_solute", solute_switch)
        self._controls.extend([rotate_switch, solute_switch])
        options_row.add(rotate_switch)
        options_row.add(solute_switch)
        parent.add(options_row)

        browse_A = toga.Button("Browse", on_press=getattr(self, f"read_box{index}_solvent_A_xyz"), style=button_style)
        browse_B = toga.Button("Browse", on_press=getattr(self, f"read_box{index}_solvent_B_xyz"), style=button_style)
        browse_solute = toga.Button("Browse", on_press=getattr(self, f"read_box{index}_solute_xyz"), style=button_style)
        setattr(self, f"btn_box{index}_solvent_A", browse_A)
        setattr(self, f"btn_box{index}_solvent_B", browse_B)
        setattr(self, f"btn_box{index}_solute", browse_solute)
        self._controls.extend([browse_A, browse_B, browse_solute])

        add_field(parent, "Solvent A:", "Required if A > 0", f"textInput_box{index}_solvent_A", button=browse_A)
        add_field(parent, "Solvent B:", "Required if B > 0", f"textInput_box{index}_solvent_B", button=browse_B)
        add_field(parent, "Solute:", "Required if solute is enabled", f"textInput_box{index}_solute", button=browse_solute)

    def update_solvent_controls(self, widget=None):
        for index in [1, 2]:
            try:
                density_str = getattr(self, f"textInput_box{index}_target_density").value.strip()
                density = float(density_str) if density_str else 1.0
            except ValueError:
                density = 1.0

            zero_density = (density == 0.0)
            editable = not zero_density and not self.is_running

            getattr(self, f"textInput_box{index}_composition_A").enabled = editable
            getattr(self, f"textInput_box{index}_composition_B").enabled = editable
            if hasattr(self, f"btn_box{index}_solvent_A"):
                getattr(self, f"btn_box{index}_solvent_A").enabled = editable
            if hasattr(self, f"btn_box{index}_solvent_B"):
                getattr(self, f"btn_box{index}_solvent_B").enabled = editable
            if hasattr(self, f"textInput_box{index}_solvent_A"):
                getattr(self, f"textInput_box{index}_solvent_A").enabled = editable
            if hasattr(self, f"textInput_box{index}_solvent_B"):
                getattr(self, f"textInput_box{index}_solvent_B").enabled = editable

    def update_solute_controls(self, widget=None):
        for index in [1, 2]:
            enabled = bool(getattr(self, f"switch_box{index}_include_solute").value)
            editable = enabled and not self.is_running
            getattr(self, f"textInput_box{index}_min_distance").enabled = editable
            if hasattr(self, f"textInput_box{index}_solute"):
                getattr(self, f"textInput_box{index}_solute").enabled = editable
            if hasattr(self, f"btn_box{index}_solute"):
                getattr(self, f"btn_box{index}_solute").enabled = editable

    def set_file_buttons_enabled(self, enabled):
        for index in [1, 2]:
            getattr(self, f"btn_box{index}_solvent_A").enabled = enabled
            getattr(self, f"btn_box{index}_solvent_B").enabled = enabled
            getattr(self, f"btn_box{index}_solute").enabled = (
                enabled and bool(getattr(self, f"switch_box{index}_include_solute").value)
            )

    def set_running_state(self, is_running):
        self.is_running = is_running
        for control in self._controls:
            control.enabled = not is_running
        self.start_button.enabled = not is_running
        self.btn_cancel.enabled = is_running
        self.btn_close.enabled = True
        self.switch_minimum_image.enabled = not is_running
        for index in [1, 2]:
            getattr(self, f"switch_box{index}_rotate").enabled = not is_running
            getattr(self, f"switch_box{index}_include_solute").enabled = not is_running
        self.set_file_buttons_enabled(not is_running)
        if is_running:
            self.cancel_requested = False
            self.latest_progress = None
            self.last_applied_progress_key = None
            self.last_progress_text_update = 0.0
            self.last_progress_text_box = None
            self.last_text_update = 0.0
            if hasattr(self, "status_label"):
                self.status_label.text = "Running shared-wall insertion... please wait"
            self.progress_bar.max = 100
            self.progress_bar.value = 0
            self.multi_line_text.value = "Preparing shared-wall random insertion...\n"
        else:
            self.update_solute_controls()
            self.update_solvent_controls()

    def cancel_run(self, widget):
        if not self.is_running:
            return
        self.cancel_requested = True
        self.btn_cancel.enabled = False
        if hasattr(self, "status_label"):
            self.status_label.text = "Cancel requested. Stopping after the current packing step..."
        self.multi_line_text.value += "\nCancel requested. No output will be written for this run.\n"

    async def workflow(self, widget):
        params_ok = await self.read_params(widget)
        if not params_ok:
            return
        files_ok = await self.validate_selected_files()
        if not files_ok:
            return

        try:
            self.set_running_state(True)
            completed = await self.run_packing_worker()
            if completed:
                await self.save_and_display_results()
            else:
                if hasattr(self, "status_label"):
                    self.status_label.text = "Cancelled"
                self.multi_line_text.value += "\nRun cancelled before output files were written.\n"
        except Exception as e:
            if hasattr(self, "status_label"):
                self.status_label.text = "Failed"
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to build shared-wall solvent box: {e}"))
            return
        finally:
            self.set_running_state(False)

    def closeTopLevel(self, widget):
        if self.is_running:
            self.cancel_run(widget)
            if hasattr(self, "status_label"):
                self.status_label.text = "Cancel requested. Close will be available after the run stops."
            return
        self.main_window.close()
