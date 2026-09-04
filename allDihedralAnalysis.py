import asyncio
import os
from dataclasses import dataclass
from typing import Callable, Iterable, List, Optional, Tuple

import numpy as np
import toga
from toga.constants import LEFT
from toga.style import Pack

# Shared, unit-tested range parser (e.g. "1-5 14-16 18 20"); reused across the
# all-bond / all-angle / all-dihedral tools so the syntax stays identical.
from allBondAnalysis import parse_solute_index_ranges


@dataclass
class DihedralQuadruplet:
    """A dihedral i-j-k-l that is connected (bonds i-j, j-k, k-l all within cutoff)
    in at least one frame.

    ``first_present_dihedral`` is the signed dihedral in the first frame the
    quadruplet appears; ``frames_present`` is how many frames it stayed connected.
    """

    atom_i: int
    atom_j: int
    atom_k: int
    atom_l: int
    element_i: str
    element_j: str
    element_k: str
    element_l: str
    first_present_dihedral: float
    frames_present: int


class AllDihedralAnalysis:
    """Calculate statistics for all connected interatomic dihedrals in an XYZ trajectory."""

    DEFAULT_CONNECTION_DISTANCE = 1.7

    def __init__(
        self,
        trajectory_file: Optional[str] = None,
        max_connection_distance: float = DEFAULT_CONNECTION_DISTANCE,
        solute_atom_indices: Optional[List[int]] = None,
        cell_lengths: Optional[Iterable[float]] = None,
    ) -> None:
        self.trajectory_file = trajectory_file
        self.max_connection_distance = float(max_connection_distance)
        self.solute_atom_indices = solute_atom_indices
        self.cell_lengths = None if cell_lengths is None else np.asarray(cell_lengths, dtype=float)
        self.num_atoms = 0
        self.num_frames = 0
        self.elements: List[str] = []
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

    def _candidate_pairs(self) -> np.ndarray:
        """All in-scope ``i<j`` atom-index pairs used for per-frame connectivity."""
        if self.solute_atom_indices is not None:
            scope = np.array(sorted(self.solute_atom_indices), dtype=int)
        else:
            scope = np.arange(self.num_atoms, dtype=int)
        upper_i, upper_j = np.triu_indices(len(scope), k=1)
        return np.column_stack((scope[upper_i], scope[upper_j]))

    def _prepare_run(self) -> np.ndarray:
        """Validate the cutoff / cell and return the candidate pair index array."""
        if self.max_connection_distance <= 0:
            raise ValueError("Maximum connection distance must be greater than zero.")
        if self.cell_lengths is not None:
            if self.cell_lengths.shape != (3,):
                raise ValueError("Cell lengths must be exactly three values: a b c.")
            if np.any(self.cell_lengths <= 0):
                raise ValueError("Cell lengths must be positive.")
            half_box = float(np.min(self.cell_lengths)) / 2.0
            if self.max_connection_distance > half_box:
                raise ValueError(
                    f"Maximum connection distance ({self.max_connection_distance}) must not "
                    f"exceed half the smallest cell length ({half_box:.4f} Angstrom) for the "
                    f"minimum-image convention to be valid."
                )
        pair_indices = self._candidate_pairs()
        if len(pair_indices) == 0:
            raise ValueError("At least two in-scope atoms are required to form a dihedral.")
        return pair_indices

    def _minimum_image(self, deltas: np.ndarray) -> np.ndarray:
        """Apply the orthorhombic minimum-image convention to bond vectors (PBC)."""
        if self.cell_lengths is not None:
            deltas = deltas - self.cell_lengths * np.round(deltas / self.cell_lengths)
        return deltas

    def _bonded_adjacency(self, coordinates: np.ndarray, pair_indices: np.ndarray):
        """Return (bonded_pairs, adjacency) for the pairs within the cutoff this frame."""
        deltas = self._minimum_image(coordinates[pair_indices[:, 0]] - coordinates[pair_indices[:, 1]])
        distances = np.sqrt(np.sum(deltas * deltas, axis=1))
        bonded = pair_indices[distances <= self.max_connection_distance]
        adjacency: dict = {}
        for atom_i, atom_j in bonded:
            adjacency.setdefault(int(atom_i), []).append(int(atom_j))
            adjacency.setdefault(int(atom_j), []).append(int(atom_i))
        return bonded, adjacency

    def _calculate_dihedrals(self, coordinates: np.ndarray, quad_indices: np.ndarray) -> np.ndarray:
        """Vectorized signed i-j-k-l dihedrals in degrees (PBC-aware); NaN for
        degenerate geometries (overlapping central atoms or a collinear plane)."""
        ab = self._minimum_image(coordinates[quad_indices[:, 1]] - coordinates[quad_indices[:, 0]])
        bc = self._minimum_image(coordinates[quad_indices[:, 2]] - coordinates[quad_indices[:, 1]])
        cd = self._minimum_image(coordinates[quad_indices[:, 3]] - coordinates[quad_indices[:, 2]])

        normal_1 = np.cross(ab, bc)
        normal_2 = np.cross(bc, cd)
        norm_bc = np.linalg.norm(bc, axis=1)
        norm_1 = np.linalg.norm(normal_1, axis=1)
        norm_2 = np.linalg.norm(normal_2, axis=1)
        degenerate = (norm_bc == 0.0) | (norm_1 == 0.0) | (norm_2 == 0.0)

        with np.errstate(invalid="ignore", divide="ignore"):
            normal_1_unit = normal_1 / norm_1[:, np.newaxis]
            normal_2_unit = normal_2 / norm_2[:, np.newaxis]
            bc_unit = bc / norm_bc[:, np.newaxis]
            x_values = np.sum(normal_1_unit * normal_2_unit, axis=1)
            y_values = np.sum(np.cross(normal_1_unit, bc_unit) * normal_2_unit, axis=1)
            dihedrals = np.degrees(np.arctan2(y_values, x_values))
        dihedrals[degenerate] = np.nan
        return dihedrals

    def _accumulate_frame(self, coordinates: np.ndarray, pair_indices: np.ndarray, acc: dict) -> None:
        """Detect this frame's connectivity, build its i-j-k-l dihedrals around each
        central bond, and fold each value into a per-quadruplet Welford accumulator
        keyed by ``(i, j, k, l)`` (central bond j<k, matching the original builder)."""
        bonded, adjacency = self._bonded_adjacency(coordinates, pair_indices)
        if len(bonded) == 0:
            return
        q_i: List[int] = []
        q_j: List[int] = []
        q_k: List[int] = []
        q_l: List[int] = []
        for atom_j, atom_k in bonded:  # candidate pairs are stored with atom_j < atom_k
            atom_j, atom_k = int(atom_j), int(atom_k)
            neighbours_j = sorted(a for a in adjacency[atom_j] if a != atom_k)
            neighbours_k = sorted(a for a in adjacency[atom_k] if a != atom_j)
            for atom_i in neighbours_j:
                for atom_l in neighbours_k:
                    if atom_i == atom_l:
                        continue
                    q_i.append(atom_i)
                    q_j.append(atom_j)
                    q_k.append(atom_k)
                    q_l.append(atom_l)
        if not q_j:
            return
        dihedrals = self._calculate_dihedrals(coordinates, np.array([q_i, q_j, q_k, q_l]).T)
        for atom_i, atom_j, atom_k, atom_l, value in zip(q_i, q_j, q_k, q_l, dihedrals):
            if np.isnan(value):
                continue
            key = (atom_i, atom_j, atom_k, atom_l)
            entry = acc.get(key)
            if entry is None:
                acc[key] = [1, float(value), 0.0, float(value)]  # count, mean, M2, first
            else:
                entry[0] += 1
                delta = value - entry[1]
                entry[1] += delta / entry[0]
                entry[2] += delta * (value - entry[1])

    def _finalize(self, acc: dict, frame_count: int) -> List[Tuple[float, float, float]]:
        """Turn the per-quadruplet accumulators into ``dihedral_quadruplets`` + ``statistics``."""
        if frame_count == 0:
            raise ValueError("No frames were read from the trajectory file.")
        if not acc:
            raise ValueError(
                "No connected dihedral quadruplets were identified in any frame "
                "(no i-j-k-l with all three bonds within the maximum connection distance)."
            )
        self.num_frames = frame_count
        self.dihedral_quadruplets = []
        self.statistics = []
        for (atom_i, atom_j, atom_k, atom_l) in sorted(acc):
            count, mean, m2, first_value = acc[(atom_i, atom_j, atom_k, atom_l)]
            variance = m2 / count
            self.dihedral_quadruplets.append(
                DihedralQuadruplet(
                    atom_i=atom_i, atom_j=atom_j, atom_k=atom_k, atom_l=atom_l,
                    element_i=self.elements[atom_i],
                    element_j=self.elements[atom_j],
                    element_k=self.elements[atom_k],
                    element_l=self.elements[atom_l],
                    first_present_dihedral=float(first_value),
                    frames_present=int(count),
                )
            )
            self.statistics.append((float(mean), float(variance), float(np.sqrt(variance))))
        return self.statistics

    def compute_dihedral_statistics(
        self, progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> List[Tuple[float, float, float]]:
        """Per-frame connectivity + per-quadruplet dihedral mean/variance/std over present frames."""
        pair_indices = self._prepare_run()
        acc: dict = {}
        frame_count = 0
        for elements, coordinates in self._read_frames():
            if len(elements) != self.num_atoms:
                raise ValueError("A frame has a different number of atoms than the first frame.")
            frame_count += 1
            self._accumulate_frame(coordinates, pair_indices, acc)
            if progress_callback:
                progress_callback(frame_count, 0)
        return self._finalize(acc, frame_count)

    async def compute_dihedral_statistics_async(
        self,
        progress_callback: Optional[Callable[[int, int], None]] = None,
        ui_update_interval: int = 500,
    ) -> List[Tuple[float, float, float]]:
        """Async-friendly variant that periodically lets the Toga UI repaint."""
        pair_indices = self._prepare_run()
        acc: dict = {}
        frame_count = 0
        for elements, coordinates in self._read_frames():
            if len(elements) != self.num_atoms:
                raise ValueError("A frame has a different number of atoms than the first frame.")
            frame_count += 1
            self._accumulate_frame(coordinates, pair_indices, acc)
            if progress_callback and frame_count % ui_update_interval == 0:
                progress_callback(frame_count, 0)
                await asyncio.sleep(0)
        stats = self._finalize(acc, frame_count)
        if progress_callback:
            progress_callback(frame_count, frame_count)
            await asyncio.sleep(0)
        return stats

    def analyze(
        self, progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> List[Tuple[float, float, float]]:
        """Run the full all-dihedral analysis."""
        self.read_first_frame()
        return self.compute_dihedral_statistics(progress_callback=progress_callback)

    async def analyze_async(
        self,
        progress_callback: Optional[Callable[[int, int], None]] = None,
        ui_update_interval: int = 500,
    ) -> List[Tuple[float, float, float]]:
        """Run the full all-dihedral analysis while keeping the Toga UI responsive."""
        self.read_first_frame()
        return await self.compute_dihedral_statistics_async(
            progress_callback=progress_callback, ui_update_interval=ui_update_interval
        )

    def write_results(self, output_file: str) -> str:
        """Write a single self-describing dihedral-analysis file (metadata header +
        one row per quadruplet: identity, first-present dihedral, stats, occurrence)."""
        if not self.statistics:
            raise ValueError("No statistics are available to write.")

        output_file = os.path.abspath(output_file)
        os.makedirs(os.path.dirname(output_file) or os.getcwd(), exist_ok=True)

        with open(output_file, "w") as out:
            out.write("# gQTEA All Dihedral Angle Analysis\n")
            out.write(f"# frames_used {self.num_frames}\n")
            out.write(f"# max_connection_distance {self.max_connection_distance:g}\n")
            out.write("# dihedral_convention signed_degrees_minus180_to_180\n")
            if self.solute_atom_indices is None:
                out.write("# atom_scope all_atoms\n")
            else:
                solute_labels = " ".join(str(index + 1) for index in sorted(self.solute_atom_indices))
                out.write("# atom_scope solute_atoms\n")
                out.write(f"# solute_atom_indices {solute_labels}\n")
            if self.cell_lengths is None:
                out.write("# periodic_boundary none\n")
            else:
                a, b, c = self.cell_lengths
                out.write(f"# cell_lengths {a:g} {b:g} {c:g}\n")
            out.write(
                "# row atom_i atom_j atom_k atom_l element_i element_j element_k element_l "
                "first_present_dihedral_degrees average_dihedral_degrees variance_degrees2 "
                "standard_deviation_degrees occurrence_fraction frames_present\n"
            )
            for row_index, (quad, (average, variance, std_dev)) in enumerate(
                zip(self.dihedral_quadruplets, self.statistics), start=1
            ):
                occurrence = quad.frames_present / self.num_frames
                out.write(
                    f"{row_index:>6d} "
                    f"{quad.atom_i + 1:>8d} {quad.atom_j + 1:>8d} "
                    f"{quad.atom_k + 1:>8d} {quad.atom_l + 1:>8d} "
                    f"{quad.element_i:>8s} {quad.element_j:>8s} "
                    f"{quad.element_k:>8s} {quad.element_l:>8s} "
                    f"{quad.first_present_dihedral:>16.8f} "
                    f"{average:>16.8f} {variance:>16.8f} {std_dev:>16.8f} "
                    f"{occurrence:>16.8f} {quad.frames_present:>10d}\n"
                )

        return output_file


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
            placeholder="Default: 1.7 (leave blank to use it)",
            style=input_style,
        )
        distance_row.add(distance_label)
        distance_row.add(self.textInput_max_distance)
        main_box.add(distance_row)

        cell_row = toga.Box(style=row_style)
        cell_label = toga.Label("Cell lattices a b c (A):", style=label_style)
        self.textInput_cell = toga.TextInput(
            placeholder="Example: 12.5 12.5 12.5; leave blank for no PBC",
            style=input_style,
        )
        cell_row.add(cell_label)
        cell_row.add(self.textInput_cell)
        main_box.add(cell_row)

        solute_row = toga.Box(style=row_style)
        solute_label = toga.Label("Solute atom indices:", style=label_style)
        self.textInput_solute_indices = toga.TextInput(
            placeholder="Example: 1-5 14-16 18 20; leave blank to use all atoms",
            style=input_style,
        )
        solute_row.add(solute_label)
        solute_row.add(self.textInput_solute_indices)
        main_box.add(solute_row)

        output_row = toga.Box(style=row_style)
        output_label = toga.Label("Output txt filename:", style=label_style)
        self.textInput_output = toga.TextInput(
            value="all_dihedral_analysis_combined.txt",
            placeholder="all_dihedral_analysis_combined.txt",
            style=input_style,
        )
        output_row.add(output_label)
        output_row.add(self.textInput_output)
        main_box.add(output_row)

        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(10, 0), font_size=12)
        )
        self.multi_line_text.value = (
            "This module re-evaluates connectivity every frame of an XYZ trajectory, builds "
            "every connected i-j-k-l dihedral path, and computes the average signed dihedral, "
            "population variance, standard deviation (degrees), and occurrence (fraction of "
            "frames the quadruplet is connected) for each. Enter cell lattices a b c to apply "
            "the minimum-image convention (PBC); leave blank for an isolated system. Solute "
            "atom indices accept ranges (e.g. 1-5 14-16 18 20) to exclude solvent atoms."
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

        # Maximum connection distance: a blank field falls back to the 1.7 A default.
        distance_text = self.textInput_max_distance.value.strip()
        if not distance_text:
            self.max_connection_distance = AllDihedralAnalysis.DEFAULT_CONNECTION_DISTANCE
        else:
            try:
                self.max_connection_distance = float(distance_text)
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

        # Cell lattices for PBC: optional. Blank = no periodic boundaries.
        cell_text = self.textInput_cell.value.strip()
        self.cell_lengths = None
        if cell_text:
            try:
                lengths = [float(value) for value in cell_text.split()]
            except ValueError:
                await self.warning_function(
                    "Error", "Cell lattices must be three numbers separated by spaces: a b c."
                )
                return False
            if len(lengths) != 3 or any(length <= 0 for length in lengths):
                await self.warning_function(
                    "Error", "Enter exactly three positive cell lattices: a b c."
                )
                return False
            self.cell_lengths = lengths

        # Solute atom indices: optional, range syntax allowed (e.g. 1-5 14-16 18 20).
        solute_index_text = self.textInput_solute_indices.value.strip()
        self.solute_atom_indices = None
        if solute_index_text:
            try:
                solute_indices = parse_solute_index_ranges(solute_index_text)
            except ValueError as exc:
                await self.warning_function("Error", str(exc))
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
                cell_lengths=self.cell_lengths,
            )

            self.multi_line_text.value = (
                "Preparing all-dihedral analysis...\n"
                "Reading the first trajectory frame and checking the selected solute atoms."
            )
            await asyncio.sleep(0)

            analyzer.read_first_frame()

            self.multi_line_text.value = (
                f"First frame loaded ({analyzer.num_atoms} atoms).\n"
                "Re-evaluating connectivity every frame, building all connected i-j-k-l "
                "dihedral paths, and computing the average signed dihedral, population "
                "variance, standard deviation, and occurrence for each quadruplet.\n"
                "Please wait until the final summary appears."
            )
            await asyncio.sleep(0)

            await analyzer.compute_dihedral_statistics_async(ui_update_interval=500)

            self.multi_line_text.value = (
                "Dihedral statistics completed.\n"
                "Writing the combined dihedral-analysis file."
            )
            await asyncio.sleep(0)

            output_file = analyzer.write_results(self.output_file)
        except Exception as exc:
            await self.warning_function("Error", f"All dihedral analysis failed: {exc}")
            return

        preview_quadruplets = []
        for row_index, quadruplet in enumerate(analyzer.dihedral_quadruplets[:10], start=1):
            occurrence = quadruplet.frames_present / analyzer.num_frames
            preview_quadruplets.append(
                f"{row_index:>3d}: {quadruplet.element_i}{quadruplet.atom_i + 1}-"
                f"{quadruplet.element_j}{quadruplet.atom_j + 1}-"
                f"{quadruplet.element_k}{quadruplet.atom_k + 1}-"
                f"{quadruplet.element_l}{quadruplet.atom_l + 1} "
                f"first-present dihedral = {quadruplet.first_present_dihedral:.6f} deg, "
                f"occurrence = {occurrence:.1%}"
            )
        quadruplet_preview = "\n".join(preview_quadruplets)
        if len(analyzer.dihedral_quadruplets) > 10:
            quadruplet_preview += "\n..."

        if analyzer.solute_atom_indices is None:
            atom_scope = "All atoms"
        else:
            solute_labels = " ".join(str(index + 1) for index in sorted(analyzer.solute_atom_indices))
            atom_scope = f"Solute atoms only: {solute_labels}"

        if analyzer.cell_lengths is None:
            pbc_state = "off (isolated system)"
        else:
            a, b, c = analyzer.cell_lengths
            pbc_state = f"on, minimum-image cell = {a:g} {b:g} {c:g} A"

        self.multi_line_text.value = (
            f"Analysis completed.\n"
            f"Trajectory: {self.trajec}\n"
            f"Frames processed: {analyzer.num_frames}\n"
            f"Atoms per frame: {analyzer.num_atoms}\n"
            f"Atom scope: {atom_scope}\n"
            f"Periodic boundaries: {pbc_state}\n"
            f"Dihedral quadruplets (present in >=1 frame): {len(analyzer.dihedral_quadruplets)}\n"
            f"Maximum connection distance: {self.max_connection_distance:.6f} A\n\n"
            f"Output file:\n{output_file}\n\n"
            f"First mapped quadruplets:\n{quadruplet_preview}"
        )

    def closeTopLevel(self, widget) -> None:
        self.main_window.close()
