import asyncio
import os
from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple

import numpy as np
import toga
from toga.constants import LEFT
from toga.style import Pack


@dataclass
class DihedralQuadruplet:
    """A dihedral quadruplet detected from the first-frame connectivity graph."""

    atom_i: int
    atom_j: int
    atom_k: int
    atom_l: int
    element_i: str
    element_j: str
    element_k: str
    element_l: str
    first_frame_dihedral: float


class AllDihedralAnalysis:
    """Calculate statistics for all connected interatomic dihedrals in an XYZ trajectory."""

    DEFAULT_CONNECTION_DISTANCE = 1.5

    def __init__(
        self,
        trajectory_file: Optional[str] = None,
        max_connection_distance: float = DEFAULT_CONNECTION_DISTANCE,
        solute_atom_indices: Optional[List[int]] = None,
    ) -> None:
        self.trajectory_file = trajectory_file
        self.max_connection_distance = float(max_connection_distance)
        self.solute_atom_indices = solute_atom_indices
        self.num_atoms = 0
        self.num_frames = 0
        self.elements: List[str] = []
        self.connected_atom_pairs: List[Tuple[int, int]] = []
        self.dihedral_quadruplets: List[DihedralQuadruplet] = []
        self.statistics: List[Tuple[float, float, float]] = []

    def _read_frames(self) -> Iterable[Tuple[List[str], np.ndarray]]:
        """Yield frames from a standard multi-frame XYZ trajectory."""
        if not self.trajectory_file:
            raise ValueError("No trajectory file was provided.")

        with open(self.trajectory_file, "r") as xyz_file:
            frame_index = 0
            expected_atoms = None

            while True:
                atom_count_line = xyz_file.readline()
                if not atom_count_line:
                    break
                if not atom_count_line.strip():
                    continue

                try:
                    atom_count = int(atom_count_line.strip().split()[0])
                except (ValueError, IndexError) as exc:
                    raise ValueError(
                        f"Invalid atom-count line at frame {frame_index + 1}."
                    ) from exc

                if expected_atoms is None:
                    expected_atoms = atom_count
                elif atom_count != expected_atoms:
                    raise ValueError(
                        "Inconsistent number of atoms: "
                        f"frame {frame_index + 1} has {atom_count}, expected {expected_atoms}."
                    )

                comment_line = xyz_file.readline()
                if not comment_line:
                    raise ValueError(f"Missing comment line at frame {frame_index + 1}.")

                elements = []
                coordinates = np.zeros((atom_count, 3), dtype=float)

                for atom_index in range(atom_count):
                    atom_line = xyz_file.readline()
                    if not atom_line:
                        raise ValueError(
                            "Unexpected end of file while reading "
                            f"frame {frame_index + 1}, atom {atom_index + 1}."
                        )

                    parts = atom_line.split()
                    if len(parts) < 4:
                        raise ValueError(
                            "Invalid atom line while reading "
                            f"frame {frame_index + 1}, atom {atom_index + 1}."
                        )

                    elements.append(parts[0])
                    try:
                        coordinates[atom_index] = [float(value) for value in parts[1:4]]
                    except ValueError as exc:
                        raise ValueError(
                            "Invalid coordinates while reading "
                            f"frame {frame_index + 1}, atom {atom_index + 1}."
                        ) from exc

                frame_index += 1
                yield elements, coordinates

    def read_first_frame(self) -> Tuple[List[str], np.ndarray]:
        """Read and store metadata from the first trajectory frame."""
        try:
            elements, coordinates = next(iter(self._read_frames()))
        except StopIteration as exc:
            raise ValueError("The trajectory file does not contain any frames.") from exc

        self.num_atoms = len(elements)
        self.elements = elements
        self._validate_solute_atom_indices()
        return elements, coordinates

    def _validate_solute_atom_indices(self) -> None:
        if self.solute_atom_indices is None:
            return

        if not self.solute_atom_indices:
            raise ValueError("Solute atom index list cannot be empty.")

        normalized_indices = []
        seen_indices = set()
        for atom_index in self.solute_atom_indices:
            if atom_index < 0 or atom_index >= self.num_atoms:
                raise ValueError(
                    "Solute atom indices must be between "
                    f"1 and {self.num_atoms}."
                )
            if atom_index in seen_indices:
                raise ValueError("Solute atom indices cannot contain duplicates.")
            seen_indices.add(atom_index)
            normalized_indices.append(atom_index)

        self.solute_atom_indices = normalized_indices

    def identify_connected_atom_pairs(self, coordinates: np.ndarray) -> List[Tuple[int, int]]:
        """Infer connected solute atom pairs from first-frame distances."""
        if self.max_connection_distance <= 0:
            raise ValueError("Maximum connection distance must be greater than zero.")

        cutoff_sq = self.max_connection_distance * self.max_connection_distance
        candidate_indices = (
            self.solute_atom_indices
            if self.solute_atom_indices is not None
            else list(range(self.num_atoms))
        )
        candidate_indices = sorted(candidate_indices)
        pairs: List[Tuple[int, int]] = []

        for index_i, atom_i in enumerate(candidate_indices):
            for atom_j in candidate_indices[index_i + 1:]:
                diff = coordinates[atom_i] - coordinates[atom_j]
                distance_sq = float(np.dot(diff, diff))
                if distance_sq <= cutoff_sq:
                    pairs.append((atom_i, atom_j))

        self.connected_atom_pairs = pairs
        return pairs

    def _dihedral_degrees(
        self,
        coordinates: np.ndarray,
        atom_i: int,
        atom_j: int,
        atom_k: int,
        atom_l: int,
    ) -> float:
        atom_a = coordinates[atom_i]
        atom_b = coordinates[atom_j]
        atom_c = coordinates[atom_k]
        atom_d = coordinates[atom_l]

        ab = atom_b - atom_a
        bc = atom_c - atom_b
        cd = atom_d - atom_c

        norm_bc = float(np.linalg.norm(bc))
        if norm_bc == 0.0:
            raise ValueError("Cannot calculate dihedral because central atoms overlap.")

        normal_1 = np.cross(ab, bc)
        normal_2 = np.cross(bc, cd)
        normal_1_norm = float(np.linalg.norm(normal_1))
        normal_2_norm = float(np.linalg.norm(normal_2))
        if normal_1_norm == 0.0 or normal_2_norm == 0.0:
            raise ValueError("Cannot calculate dihedral because one plane is degenerate.")

        normal_1_unit = normal_1 / normal_1_norm
        normal_2_unit = normal_2 / normal_2_norm
        bc_unit = bc / norm_bc

        x_value = float(np.dot(normal_1_unit, normal_2_unit))
        y_value = float(np.dot(np.cross(normal_1_unit, bc_unit), normal_2_unit))
        return float(np.degrees(np.arctan2(y_value, x_value)))

    def identify_dihedral_quadruplets(self, coordinates: np.ndarray) -> List[DihedralQuadruplet]:
        """Build all i-j-k-l dihedrals from connected paths around each central bond j-k."""
        if not self.connected_atom_pairs:
            raise ValueError("No connected atom pairs were identified.")

        adjacency = {atom_index: [] for atom_index in range(self.num_atoms)}
        for atom_i, atom_j in self.connected_atom_pairs:
            adjacency[atom_i].append(atom_j)
            adjacency[atom_j].append(atom_i)

        quadruplets: List[DihedralQuadruplet] = []
        for atom_j, atom_k in self.connected_atom_pairs:
            neighbors_j = sorted(atom for atom in adjacency[atom_j] if atom != atom_k)
            neighbors_k = sorted(atom for atom in adjacency[atom_k] if atom != atom_j)

            for atom_i in neighbors_j:
                for atom_l in neighbors_k:
                    if atom_i == atom_l:
                        continue
                    try:
                        first_frame_dihedral = self._dihedral_degrees(
                            coordinates, atom_i, atom_j, atom_k, atom_l
                        )
                    except ValueError:
                        continue

                    quadruplets.append(
                        DihedralQuadruplet(
                            atom_i=atom_i,
                            atom_j=atom_j,
                            atom_k=atom_k,
                            atom_l=atom_l,
                            element_i=self.elements[atom_i],
                            element_j=self.elements[atom_j],
                            element_k=self.elements[atom_k],
                            element_l=self.elements[atom_l],
                            first_frame_dihedral=first_frame_dihedral,
                        )
                    )

        if not quadruplets:
            raise ValueError("No connected dihedral quadruplets were identified.")

        self.dihedral_quadruplets = quadruplets
        return quadruplets

    def _calculate_dihedrals_for_frame(
        self,
        coordinates: np.ndarray,
        quadruplet_indices: np.ndarray,
    ) -> np.ndarray:
        atoms_a = coordinates[quadruplet_indices[:, 0]]
        atoms_b = coordinates[quadruplet_indices[:, 1]]
        atoms_c = coordinates[quadruplet_indices[:, 2]]
        atoms_d = coordinates[quadruplet_indices[:, 3]]

        ab = atoms_b - atoms_a
        bc = atoms_c - atoms_b
        cd = atoms_d - atoms_c

        norm_bc = np.linalg.norm(bc, axis=1)
        if np.any(norm_bc == 0.0):
            raise ValueError("Cannot calculate dihedral because central atoms overlap.")

        normal_1 = np.cross(ab, bc)
        normal_2 = np.cross(bc, cd)
        normal_1_norm = np.linalg.norm(normal_1, axis=1)
        normal_2_norm = np.linalg.norm(normal_2, axis=1)
        if np.any(normal_1_norm == 0.0) or np.any(normal_2_norm == 0.0):
            raise ValueError("Cannot calculate dihedral because one plane is degenerate.")

        normal_1_unit = normal_1 / normal_1_norm[:, np.newaxis]
        normal_2_unit = normal_2 / normal_2_norm[:, np.newaxis]
        bc_unit = bc / norm_bc[:, np.newaxis]

        x_values = np.sum(normal_1_unit * normal_2_unit, axis=1)
        y_values = np.sum(np.cross(normal_1_unit, bc_unit) * normal_2_unit, axis=1)
        return np.degrees(np.arctan2(y_values, x_values))

    async def compute_dihedral_statistics_async(
        self,
        ui_update_interval: int = 500,
    ) -> List[Tuple[float, float, float]]:
        """Compute average, population variance, and standard deviation for each dihedral."""
        if not self.dihedral_quadruplets:
            raise ValueError("No dihedral quadruplets were identified.")

        quadruplet_indices = np.array(
            [
                (
                    quadruplet.atom_i,
                    quadruplet.atom_j,
                    quadruplet.atom_k,
                    quadruplet.atom_l,
                )
                for quadruplet in self.dihedral_quadruplets
            ],
            dtype=int,
        )
        num_quadruplets = len(quadruplet_indices)
        means = np.zeros(num_quadruplets, dtype=float)
        m2 = np.zeros(num_quadruplets, dtype=float)
        frame_count = 0

        for elements, coordinates in self._read_frames():
            if len(elements) != self.num_atoms:
                raise ValueError("A frame has a different number of atoms than the first frame.")

            frame_count += 1
            dihedrals = self._calculate_dihedrals_for_frame(coordinates, quadruplet_indices)

            delta = dihedrals - means
            means += delta / frame_count
            delta_after_update = dihedrals - means
            m2 += delta * delta_after_update

            if frame_count % ui_update_interval == 0:
                await asyncio.sleep(0)

        if frame_count == 0:
            raise ValueError("No frames were read from the trajectory file.")

        await asyncio.sleep(0)

        variances = m2 / frame_count
        std_devs = np.sqrt(variances)
        self.num_frames = frame_count
        self.statistics = [
            (float(mean), float(variance), float(std_dev))
            for mean, variance, std_dev in zip(means, variances, std_devs)
        ]
        return self.statistics

    async def analyze_async(self, ui_update_interval: int = 500) -> List[Tuple[float, float, float]]:
        """Run the full all-dihedral analysis."""
        _, first_coordinates = self.read_first_frame()
        self.identify_connected_atom_pairs(first_coordinates)
        self.identify_dihedral_quadruplets(first_coordinates)
        return await self.compute_dihedral_statistics_async(ui_update_interval=ui_update_interval)

    def write_results(self, output_file: str, mapping_file: Optional[str] = None) -> Tuple[str, str]:
        """Write dihedral statistics and row-to-quadruplet mapping files."""
        if not self.statistics:
            raise ValueError("No statistics are available to write.")

        output_file = os.path.abspath(output_file)
        output_dir = os.path.dirname(output_file) or os.getcwd()
        os.makedirs(output_dir, exist_ok=True)

        if mapping_file is None:
            base, _ = os.path.splitext(output_file)
            mapping_file = f"{base}_quadruplets.txt"
        mapping_file = os.path.abspath(mapping_file)

        with open(output_file, "w") as stats_file:
            stats_file.write(
                "# average_dihedral_degrees variance_degrees2 standard_deviation_degrees\n"
            )
            for average, variance, std_dev in self.statistics:
                stats_file.write(f"{average:>16.8f} {variance:>16.8f} {std_dev:>16.8f}\n")

        with open(mapping_file, "w") as quadruplet_file:
            quadruplet_file.write(f"# frames_used {self.num_frames}\n")
            quadruplet_file.write("# dihedral_convention signed_degrees_minus180_to_180\n")
            if self.solute_atom_indices is None:
                quadruplet_file.write("# atom_scope all_atoms\n")
            else:
                solute_labels = " ".join(str(index + 1) for index in self.solute_atom_indices)
                quadruplet_file.write("# atom_scope solute_atoms\n")
                quadruplet_file.write(f"# solute_atom_indices {solute_labels}\n")
            quadruplet_file.write(
                "# row atom_i atom_j atom_k atom_l element_i element_j element_k element_l first_frame_dihedral_degrees\n"
            )
            for row_index, quadruplet in enumerate(self.dihedral_quadruplets, start=1):
                quadruplet_file.write(
                    f"{row_index:>6d} "
                    f"{quadruplet.atom_i + 1:>8d} {quadruplet.atom_j + 1:>8d} "
                    f"{quadruplet.atom_k + 1:>8d} {quadruplet.atom_l + 1:>8d} "
                    f"{quadruplet.element_i:>8s} {quadruplet.element_j:>8s} "
                    f"{quadruplet.element_k:>8s} {quadruplet.element_l:>8s} "
                    f"{quadruplet.first_frame_dihedral:>16.8f}\n"
                )

        return output_file, mapping_file


class allDihedralAnalysisUI:
    """Toga frontend for all connected dihedral statistics."""

    def __init__(self, *args) -> None:
        self.trajec = None
        self.output_dir = os.getcwd()
        self.layout_main_window(*args)

    async def warning_function(self, title: str, message: str) -> None:
        await self.main_window.dialog(toga.InfoDialog(title, message))

    def layout_main_window(self, widget) -> None:
        self.main_window = toga.Window(
            title="All Dihedral Angle Analysis from Molecular Dynamics Simulations",
            size=(760, 600),
        )

        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style = Pack(margin=(5, 5), text_align=LEFT, width=260)
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=120)
        row_style = Pack(direction="row", margin=(0, 0, 5, 0))

        main_box = toga.Box(style=Pack(direction="column", margin=20))

        title_row = toga.Box(style=Pack(direction="row", margin=(0, 0, 10, 0)))
        title_box = toga.Box(style=Pack(width=700))
        title_label = toga.Label("All Dihedral Angle Analysis", style=heading_style)
        title_box.add(title_label)
        title_row.add(title_box)
        main_box.add(title_row)

        file_row = toga.Box(style=row_style)
        file_label = toga.Label("Select trajectory.xyz file:", style=label_style)
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select trajectory.xyz file",
            style=input_style,
        )
        browse_button = toga.Button("Browse", on_press=self.open_file_dialog, style=button_style)
        file_row.add(file_label)
        file_row.add(self.textInput_file)
        file_row.add(browse_button)
        main_box.add(file_row)

        distance_row = toga.Box(style=row_style)
        distance_label = toga.Label("Maximum connection distance (A):", style=label_style)
        self.textInput_max_distance = toga.TextInput(
            value=str(AllDihedralAnalysis.DEFAULT_CONNECTION_DISTANCE),
            placeholder="Default: 1.5",
            style=input_style,
        )
        distance_row.add(distance_label)
        distance_row.add(self.textInput_max_distance)
        main_box.add(distance_row)

        solute_row = toga.Box(style=row_style)
        solute_label = toga.Label("Solute atom indices:", style=label_style)
        self.textInput_solute_indices = toga.TextInput(
            placeholder="Example: 1 2 3 4; leave blank to use all atoms",
            style=input_style,
        )
        solute_row.add(solute_label)
        solute_row.add(self.textInput_solute_indices)
        main_box.add(solute_row)

        output_row = toga.Box(style=row_style)
        output_label = toga.Label("Output txt filename:", style=label_style)
        self.textInput_output = toga.TextInput(
            value="all_dihedral_analysis.txt",
            placeholder="all_dihedral_analysis.txt",
            style=input_style,
        )
        output_row.add(output_label)
        output_row.add(self.textInput_output)
        main_box.add(output_row)

        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(10, 0), font_size=12)
        )
        self.multi_line_text.value = (
            "This module infers connected atom pairs from the first frame of an XYZ "
            "trajectory, builds every connected i-j-k-l dihedral path, and computes "
            "the average signed dihedral, population variance, and standard deviation "
            "in degrees. Provide solute atom indices to exclude solvent atoms."
        )
        main_box.add(self.multi_line_text)

        button_row = toga.Box(style=Pack(direction="row", margin=(10, 0, 0, 0)))
        self.btn_execute = toga.Button("Analyze", style=button_style, on_press=self.workflow)
        self.btn_close = toga.Button("Close", style=button_style, on_press=self.closeTopLevel)
        button_row.add(self.btn_execute)
        button_row.add(self.btn_close)
        main_box.add(button_row)

        self.main_window.content = main_box
        self.main_window.show()

    async def open_file_dialog(self, widget) -> None:
        try:
            selected_file = await self.main_window.dialog(
                toga.OpenFileDialog("Open trajectory.xyz file")
            )
        except ValueError:
            await self.warning_function("Error", "Open file was canceled.")
            return

        if selected_file is None:
            await self.warning_function("Warning", "No file was selected.")
            return

        self.trajec = str(selected_file)
        self.output_dir = os.path.dirname(os.path.abspath(self.trajec)) or os.getcwd()
        self.textInput_file.value = self.trajec

        try:
            analyzer = AllDihedralAnalysis(self.trajec)
            analyzer.read_first_frame()
        except Exception as exc:
            await self.warning_function("Error", f"Failed to read trajectory file: {exc}")
            return

        self.multi_line_text.value = (
            f"Selected file: {self.trajec}\n"
            f"Number of atoms in first frame: {analyzer.num_atoms}\n"
            f"Output directory: {self.output_dir}\n"
        )

    async def read_params(self) -> bool:
        if not self.trajec:
            await self.warning_function("Error", "No trajectory file selected.")
            return False

        try:
            self.max_connection_distance = float(self.textInput_max_distance.value.strip())
        except ValueError:
            await self.warning_function(
                "Error", "Maximum connection distance must be a valid number."
            )
            return False

        if self.max_connection_distance <= 0:
            await self.warning_function(
                "Error", "Maximum connection distance must be greater than zero."
            )
            return False

        solute_index_text = self.textInput_solute_indices.value.strip()
        self.solute_atom_indices = None
        if solute_index_text:
            try:
                solute_indices = [int(value) for value in solute_index_text.split()]
            except ValueError:
                await self.warning_function(
                    "Error", "Solute atom indices must be positive integers separated by spaces."
                )
                return False

            if any(index <= 0 for index in solute_indices):
                await self.warning_function(
                    "Error", "Solute atom indices must be positive integers starting at 1."
                )
                return False

            self.solute_atom_indices = [index - 1 for index in solute_indices]

        output_name = self.textInput_output.value.strip()
        if not output_name:
            await self.warning_function("Error", "Please provide an output filename.")
            return False

        if os.path.isabs(output_name):
            self.output_file = output_name
        else:
            self.output_file = os.path.join(self.output_dir, output_name)

        return True

    async def workflow(self, widget) -> None:
        if not await self.read_params():
            return

        try:
            analyzer = AllDihedralAnalysis(
                trajectory_file=self.trajec,
                max_connection_distance=self.max_connection_distance,
                solute_atom_indices=self.solute_atom_indices,
            )

            self.multi_line_text.value = (
                "Preparing all-dihedral analysis...\n"
                "Reading the first trajectory frame and checking the selected solute atoms."
            )
            await asyncio.sleep(0)

            _, first_coordinates = analyzer.read_first_frame()

            self.multi_line_text.value = (
                "First frame loaded.\n"
                "Identifying connected atom pairs using the maximum connection distance.\n"
                "Only solute-solute connectivity is considered when solute atom indices are provided."
            )
            await asyncio.sleep(0)

            analyzer.identify_connected_atom_pairs(first_coordinates)

            self.multi_line_text.value = (
                f"Connected pairs identified: {len(analyzer.connected_atom_pairs)}\n"
                "Building all connected i-j-k-l dihedral paths around each central bond."
            )
            await asyncio.sleep(0)

            analyzer.identify_dihedral_quadruplets(first_coordinates)

            self.multi_line_text.value = (
                f"Dihedral quadruplets identified: {len(analyzer.dihedral_quadruplets)}\n"
                "Computing signed dihedral values across all trajectory frames.\n"
                "Calculating the average signed dihedral, population variance, and "
                "standard deviation for each selected quadruplet. Please wait until "
                "the final summary appears."
            )
            await asyncio.sleep(0)

            await analyzer.compute_dihedral_statistics_async(ui_update_interval=500)

            self.multi_line_text.value = (
                "Dihedral statistics completed.\n"
                "Writing the statistics file and the atom-quadruplet mapping file."
            )
            await asyncio.sleep(0)

            output_file, mapping_file = analyzer.write_results(self.output_file)
        except Exception as exc:
            await self.warning_function("Error", f"All dihedral analysis failed: {exc}")
            return

        preview_quadruplets = []
        for row_index, quadruplet in enumerate(analyzer.dihedral_quadruplets[:10], start=1):
            preview_quadruplets.append(
                f"{row_index:>3d}: {quadruplet.element_i}{quadruplet.atom_i + 1}-"
                f"{quadruplet.element_j}{quadruplet.atom_j + 1}-"
                f"{quadruplet.element_k}{quadruplet.atom_k + 1}-"
                f"{quadruplet.element_l}{quadruplet.atom_l + 1} "
                f"first-frame dihedral = {quadruplet.first_frame_dihedral:.6f} deg"
            )
        quadruplet_preview = "\n".join(preview_quadruplets)
        if len(analyzer.dihedral_quadruplets) > 10:
            quadruplet_preview += "\n..."

        if analyzer.solute_atom_indices is None:
            atom_scope = "All atoms"
        else:
            solute_labels = " ".join(str(index + 1) for index in analyzer.solute_atom_indices)
            atom_scope = f"Solute atoms only: {solute_labels}"

        self.multi_line_text.value = (
            f"Analysis completed.\n"
            f"Trajectory: {self.trajec}\n"
            f"Frames processed: {analyzer.num_frames}\n"
            f"Atoms per frame: {analyzer.num_atoms}\n"
            f"Atom scope: {atom_scope}\n"
            f"Connected pairs: {len(analyzer.connected_atom_pairs)}\n"
            f"Dihedral quadruplets: {len(analyzer.dihedral_quadruplets)}\n"
            f"Maximum connection distance: {self.max_connection_distance:.6f} A\n\n"
            f"Statistics file:\n{output_file}\n\n"
            f"Quadruplet mapping file:\n{mapping_file}\n\n"
            f"First mapped quadruplets:\n{quadruplet_preview}"
        )

    def closeTopLevel(self, widget) -> None:
        self.main_window.close()
