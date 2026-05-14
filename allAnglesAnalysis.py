import asyncio
import os
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, List, Optional, Tuple

import numpy as np
import toga
from toga.constants import LEFT
from toga.style import Pack


@dataclass
class AngleTriplet:
    """A bond angle triplet detected from the first-frame connectivity graph."""

    atom_i: int
    atom_j: int
    atom_k: int
    element_i: str
    element_j: str
    element_k: str
    first_frame_angle: float


class AllAnglesAnalysis:
    """Calculate statistics for all connected interatomic angles in an XYZ trajectory."""

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
        self.angle_triplets: List[AngleTriplet] = []
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
        pairs: List[Tuple[int, int]] = []

        for index_i, atom_i in enumerate(candidate_indices):
            for atom_j in candidate_indices[index_i + 1:]:
                diff = coordinates[atom_i] - coordinates[atom_j]
                distance_sq = float(np.dot(diff, diff))
                if distance_sq <= cutoff_sq:
                    pairs.append((atom_i, atom_j))

        self.connected_atom_pairs = pairs
        return pairs

    def _angle_degrees(self, coordinates: np.ndarray, atom_i: int, atom_j: int, atom_k: int) -> float:
        vector_ji = coordinates[atom_i] - coordinates[atom_j]
        vector_jk = coordinates[atom_k] - coordinates[atom_j]
        norm_ji = float(np.linalg.norm(vector_ji))
        norm_jk = float(np.linalg.norm(vector_jk))

        if norm_ji == 0.0 or norm_jk == 0.0:
            raise ValueError(
                "Cannot calculate angle because two atoms occupy the same coordinates."
            )

        cosine = float(np.dot(vector_ji, vector_jk) / (norm_ji * norm_jk))
        cosine = max(-1.0, min(1.0, cosine))
        return float(np.degrees(np.arccos(cosine)))

    def identify_angle_triplets(self, coordinates: np.ndarray) -> List[AngleTriplet]:
        """Build all i-j-k angles where i and k are connected to central atom j."""
        if not self.connected_atom_pairs:
            raise ValueError("No connected atom pairs were identified.")

        adjacency = {atom_index: [] for atom_index in range(self.num_atoms)}
        for atom_i, atom_j in self.connected_atom_pairs:
            adjacency[atom_i].append(atom_j)
            adjacency[atom_j].append(atom_i)

        triplets: List[AngleTriplet] = []
        for atom_j in sorted(adjacency):
            neighbors = sorted(adjacency[atom_j])
            for atom_i, atom_k in combinations(neighbors, 2):
                triplets.append(
                    AngleTriplet(
                        atom_i=atom_i,
                        atom_j=atom_j,
                        atom_k=atom_k,
                        element_i=self.elements[atom_i],
                        element_j=self.elements[atom_j],
                        element_k=self.elements[atom_k],
                        first_frame_angle=self._angle_degrees(
                            coordinates, atom_i, atom_j, atom_k
                        ),
                    )
                )

        if not triplets:
            raise ValueError("No connected angle triplets were identified.")

        self.angle_triplets = triplets
        return triplets

    def _calculate_angles_for_frame(self, coordinates: np.ndarray, triplet_indices: np.ndarray) -> np.ndarray:
        coords_i = coordinates[triplet_indices[:, 0]]
        coords_j = coordinates[triplet_indices[:, 1]]
        coords_k = coordinates[triplet_indices[:, 2]]

        vectors_ji = coords_i - coords_j
        vectors_jk = coords_k - coords_j
        norms_ji = np.linalg.norm(vectors_ji, axis=1)
        norms_jk = np.linalg.norm(vectors_jk, axis=1)

        if np.any(norms_ji == 0.0) or np.any(norms_jk == 0.0):
            raise ValueError(
                "Cannot calculate angle because two atoms occupy the same coordinates."
            )

        cosines = np.sum(vectors_ji * vectors_jk, axis=1) / (norms_ji * norms_jk)
        cosines = np.clip(cosines, -1.0, 1.0)
        return np.degrees(np.arccos(cosines))

    async def compute_angle_statistics_async(
        self,
        ui_update_interval: int = 500,
    ) -> List[Tuple[float, float, float]]:
        """Compute average, population variance, and standard deviation for each angle."""
        if not self.angle_triplets:
            raise ValueError("No angle triplets were identified.")

        triplet_indices = np.array(
            [
                (triplet.atom_i, triplet.atom_j, triplet.atom_k)
                for triplet in self.angle_triplets
            ],
            dtype=int,
        )
        num_triplets = len(triplet_indices)
        means = np.zeros(num_triplets, dtype=float)
        m2 = np.zeros(num_triplets, dtype=float)
        frame_count = 0

        for elements, coordinates in self._read_frames():
            if len(elements) != self.num_atoms:
                raise ValueError("A frame has a different number of atoms than the first frame.")

            frame_count += 1
            angles = self._calculate_angles_for_frame(coordinates, triplet_indices)

            delta = angles - means
            means += delta / frame_count
            delta_after_update = angles - means
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
        """Run the full all-angle analysis."""
        _, first_coordinates = self.read_first_frame()
        self.identify_connected_atom_pairs(first_coordinates)
        self.identify_angle_triplets(first_coordinates)
        return await self.compute_angle_statistics_async(ui_update_interval=ui_update_interval)

    def write_results(self, output_file: str, mapping_file: Optional[str] = None) -> Tuple[str, str]:
        """Write angle statistics and row-to-triplet mapping files."""
        if not self.statistics:
            raise ValueError("No statistics are available to write.")

        output_file = os.path.abspath(output_file)
        output_dir = os.path.dirname(output_file) or os.getcwd()
        os.makedirs(output_dir, exist_ok=True)

        if mapping_file is None:
            base, _ = os.path.splitext(output_file)
            mapping_file = f"{base}_triplets.txt"
        mapping_file = os.path.abspath(mapping_file)

        with open(output_file, "w") as stats_file:
            stats_file.write(
                "# average_angle_degrees variance_degrees2 standard_deviation_degrees\n"
            )
            for average, variance, std_dev in self.statistics:
                stats_file.write(f"{average:>16.8f} {variance:>16.8f} {std_dev:>16.8f}\n")

        with open(mapping_file, "w") as triplet_file:
            triplet_file.write(f"# frames_used {self.num_frames}\n")
            if self.solute_atom_indices is None:
                triplet_file.write("# atom_scope all_atoms\n")
            else:
                solute_labels = " ".join(str(index + 1) for index in self.solute_atom_indices)
                triplet_file.write("# atom_scope solute_atoms\n")
                triplet_file.write(f"# solute_atom_indices {solute_labels}\n")
            triplet_file.write(
                "# row atom_i atom_j atom_k element_i element_j element_k first_frame_angle_degrees\n"
            )
            for row_index, triplet in enumerate(self.angle_triplets, start=1):
                triplet_file.write(
                    f"{row_index:>6d} "
                    f"{triplet.atom_i + 1:>8d} {triplet.atom_j + 1:>8d} {triplet.atom_k + 1:>8d} "
                    f"{triplet.element_i:>8s} {triplet.element_j:>8s} {triplet.element_k:>8s} "
                    f"{triplet.first_frame_angle:>16.8f}\n"
                )

        return output_file, mapping_file


class allAnglesAnalysisUI:
    """Toga frontend for all connected bond-angle statistics."""

    def __init__(self, *args) -> None:
        self.trajec = None
        self.output_dir = os.getcwd()
        self.layout_main_window(*args)

    async def warning_function(self, title: str, message: str) -> None:
        await self.main_window.dialog(toga.InfoDialog(title, message))

    def layout_main_window(self, widget) -> None:
        self.main_window = toga.Window(
            title="All Bond Angle Analysis from Molecular Dynamics Simulations",
            size=(740, 580),
        )

        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style = Pack(margin=(5, 5), text_align=LEFT, width=250)
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=120)
        row_style = Pack(direction="row", margin=(0, 0, 5, 0))

        main_box = toga.Box(style=Pack(direction="column", margin=20))

        title_row = toga.Box(style=Pack(direction="row", margin=(0, 0, 10, 0)))
        title_box = toga.Box(style=Pack(width=680))
        title_label = toga.Label("All Bond Angle Analysis", style=heading_style)
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
            value=str(AllAnglesAnalysis.DEFAULT_CONNECTION_DISTANCE),
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
            value="all_angles_analysis.txt",
            placeholder="all_angles_analysis.txt",
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
            "trajectory, builds every connected i-j-k bond angle, and computes the "
            "average angle, population variance, and standard deviation in degrees. "
            "Provide solute atom indices to exclude solvent atoms from the calculation."
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
            analyzer = AllAnglesAnalysis(self.trajec)
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
            analyzer = AllAnglesAnalysis(
                trajectory_file=self.trajec,
                max_connection_distance=self.max_connection_distance,
                solute_atom_indices=self.solute_atom_indices,
            )

            self.multi_line_text.value = (
                "Preparing all-angle analysis...\n"
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
                "Building all connected i-j-k angle triplets around each central atom."
            )
            await asyncio.sleep(0)

            analyzer.identify_angle_triplets(first_coordinates)

            self.multi_line_text.value = (
                f"Angle triplets identified: {len(analyzer.angle_triplets)}\n"
                "Computing angle values across all trajectory frames.\n"
                "Calculating the average angle, population variance, and standard deviation "
                "for each selected triplet. Please wait until the final summary appears."
            )
            await asyncio.sleep(0)

            await analyzer.compute_angle_statistics_async(ui_update_interval=500)

            self.multi_line_text.value = (
                "Angle statistics completed.\n"
                "Writing the statistics file and the atom-triplet mapping file."
            )
            await asyncio.sleep(0)

            output_file, mapping_file = analyzer.write_results(self.output_file)
        except Exception as exc:
            await self.warning_function("Error", f"All angle analysis failed: {exc}")
            return

        preview_triplets = []
        for row_index, triplet in enumerate(analyzer.angle_triplets[:10], start=1):
            preview_triplets.append(
                f"{row_index:>3d}: {triplet.element_i}{triplet.atom_i + 1}-"
                f"{triplet.element_j}{triplet.atom_j + 1}-"
                f"{triplet.element_k}{triplet.atom_k + 1} "
                f"first-frame angle = {triplet.first_frame_angle:.6f} deg"
            )
        triplet_preview = "\n".join(preview_triplets)
        if len(analyzer.angle_triplets) > 10:
            triplet_preview += "\n..."

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
            f"Angle triplets: {len(analyzer.angle_triplets)}\n"
            f"Maximum connection distance: {self.max_connection_distance:.6f} A\n\n"
            f"Statistics file:\n{output_file}\n\n"
            f"Triplet mapping file:\n{mapping_file}\n\n"
            f"First mapped triplets:\n{triplet_preview}"
        )

    def closeTopLevel(self, widget) -> None:
        self.main_window.close()
