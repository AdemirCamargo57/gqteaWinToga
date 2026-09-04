import asyncio
import os
from dataclasses import dataclass
from typing import Callable, Iterable, List, Optional, Tuple

import numpy as np
import toga
from toga.constants import LEFT
from toga.style import Pack


def parse_solute_index_ranges(text: str) -> List[int]:
    """Expand compact solute-index syntax into an explicit 1-based list.

    Accepts whitespace- and/or comma-separated tokens, where each token is
    either a single positive integer (``18``) or an inclusive range
    (``14-16``). Returns a sorted, de-duplicated list of 1-based indices; an
    empty/blank string returns ``[]``. Raises ``ValueError`` on non-positive
    values, descending ranges, or malformed tokens.
    """
    indices = set()
    for token in text.replace(",", " ").split():
        if "-" in token:
            parts = token.split("-")
            if len(parts) != 2 or not parts[0] or not parts[1]:
                raise ValueError(f"Invalid index range: '{token}'.")
            try:
                low, high = int(parts[0]), int(parts[1])
            except ValueError as exc:
                raise ValueError(f"Invalid index range: '{token}'.") from exc
            if low <= 0 or high <= 0:
                raise ValueError("Solute atom indices must be positive integers starting at 1.")
            if high < low:
                raise ValueError(f"Range bounds must be ascending: '{token}'.")
            indices.update(range(low, high + 1))
        else:
            try:
                value = int(token)
            except ValueError as exc:
                raise ValueError(f"Invalid solute atom index: '{token}'.") from exc
            if value <= 0:
                raise ValueError("Solute atom indices must be positive integers starting at 1.")
            indices.add(value)
    return sorted(indices)


@dataclass
class BondPair:
    """An atom pair that is within the connection cutoff in at least one frame.

    ``first_bonded_distance`` is the separation in the first frame where the
    pair became bonded; ``frames_bonded`` is how many frames it stayed within
    the cutoff (its occurrence count).
    """

    atom_i: int
    atom_j: int
    element_i: str
    element_j: str
    first_bonded_distance: float
    frames_bonded: int


class AllBondAnalysis:
    """Calculate statistics for all connected interatomic distances in an XYZ trajectory."""

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
        self.connected_atom_pairs: List[BondPair] = []
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
        """All in-scope ``i<j`` atom-index pairs (solute subset, or every atom).

        Returned as an ``(P, 2)`` int array with ``atom_i < atom_j`` in every
        row; connectivity itself is decided per frame in the streaming pass.
        """
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
            raise ValueError("At least two in-scope atoms are required to form a pair.")
        return pair_indices

    def _pair_distances(self, coordinates: np.ndarray, pair_indices: np.ndarray) -> np.ndarray:
        """Vectorized distances for every candidate pair, applying PBC if a cell is set."""
        deltas = coordinates[pair_indices[:, 0]] - coordinates[pair_indices[:, 1]]
        if self.cell_lengths is not None:
            deltas -= self.cell_lengths * np.round(deltas / self.cell_lengths)
        return np.sqrt(np.sum(deltas * deltas, axis=1))

    def _finalize(self, pair_indices, counts, means, m2, first_dist, frame_count):
        """Turn the per-pair accumulators into ``connected_atom_pairs`` + ``statistics``."""
        if frame_count == 0:
            raise ValueError("No frames were read from the trajectory file.")

        bonded_ever = counts > 0
        if not np.any(bonded_ever):
            raise ValueError(
                "No connected atom pairs were identified in any frame "
                "(no pair came within the maximum connection distance)."
            )

        self.num_frames = frame_count
        self.connected_atom_pairs = []
        self.statistics = []
        for row in np.flatnonzero(bonded_ever):
            atom_i, atom_j = int(pair_indices[row, 0]), int(pair_indices[row, 1])
            count = int(counts[row])
            variance = float(m2[row] / count)
            self.connected_atom_pairs.append(
                BondPair(
                    atom_i=atom_i,
                    atom_j=atom_j,
                    element_i=self.elements[atom_i],
                    element_j=self.elements[atom_j],
                    first_bonded_distance=float(first_dist[row]),
                    frames_bonded=count,
                )
            )
            self.statistics.append((float(means[row]), variance, float(np.sqrt(variance))))
        return self.statistics

    @staticmethod
    def _accumulate(distances, bonded, counts, means, m2, first_dist):
        """One masked Welford update for the pairs bonded in the current frame."""
        counts[bonded] += 1
        newly = bonded & (counts == 1)
        first_dist[newly] = distances[newly]
        delta = distances[bonded] - means[bonded]
        means[bonded] += delta / counts[bonded]
        m2[bonded] += delta * (distances[bonded] - means[bonded])

    def compute_bond_statistics(
        self,
        progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> List[Tuple[float, float, float]]:
        """Per-frame connectivity + per-pair mean/variance/std over its bonded frames."""
        pair_indices = self._prepare_run()
        num_pairs = len(pair_indices)
        counts = np.zeros(num_pairs, dtype=int)
        means = np.zeros(num_pairs, dtype=float)
        m2 = np.zeros(num_pairs, dtype=float)
        first_dist = np.full(num_pairs, np.nan, dtype=float)
        frame_count = 0

        for elements, coordinates in self._read_frames():
            if len(elements) != self.num_atoms:
                raise ValueError("A frame has a different number of atoms than the first frame.")
            frame_count += 1
            distances = self._pair_distances(coordinates, pair_indices)
            bonded = distances <= self.max_connection_distance
            self._accumulate(distances, bonded, counts, means, m2, first_dist)
            if progress_callback:
                progress_callback(frame_count, 0)

        return self._finalize(pair_indices, counts, means, m2, first_dist, frame_count)

    async def compute_bond_statistics_async(
        self,
        progress_callback: Optional[Callable[[int, int], None]] = None,
        ui_update_interval: int = 500,
    ) -> List[Tuple[float, float, float]]:
        """Async-friendly variant that periodically lets the Toga UI repaint."""
        pair_indices = self._prepare_run()
        num_pairs = len(pair_indices)
        counts = np.zeros(num_pairs, dtype=int)
        means = np.zeros(num_pairs, dtype=float)
        m2 = np.zeros(num_pairs, dtype=float)
        first_dist = np.full(num_pairs, np.nan, dtype=float)
        frame_count = 0

        for elements, coordinates in self._read_frames():
            if len(elements) != self.num_atoms:
                raise ValueError("A frame has a different number of atoms than the first frame.")
            frame_count += 1
            distances = self._pair_distances(coordinates, pair_indices)
            bonded = distances <= self.max_connection_distance
            self._accumulate(distances, bonded, counts, means, m2, first_dist)
            if progress_callback and frame_count % ui_update_interval == 0:
                progress_callback(frame_count, 0)
                await asyncio.sleep(0)

        stats = self._finalize(pair_indices, counts, means, m2, first_dist, frame_count)
        if progress_callback:
            progress_callback(frame_count, frame_count)
            await asyncio.sleep(0)
        return stats

    def analyze(
        self,
        progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> List[Tuple[float, float, float]]:
        """Run the full all-bond distance analysis."""
        self.read_first_frame()
        return self.compute_bond_statistics(progress_callback=progress_callback)

    async def analyze_async(
        self,
        progress_callback: Optional[Callable[[int, int], None]] = None,
        ui_update_interval: int = 500,
    ) -> List[Tuple[float, float, float]]:
        """Run the full analysis while keeping the Toga UI responsive."""
        self.read_first_frame()
        return await self.compute_bond_statistics_async(
            progress_callback=progress_callback,
            ui_update_interval=ui_update_interval,
        )

    def write_results(self, output_file: str) -> str:
        """Write a single self-describing bond-analysis file.

        Combines what used to be two files: a commented metadata header
        (frames used, cutoff, atom scope, periodic-boundary state) followed by
        one row per connected pair carrying its identity, first-bonded
        distance, distance statistics, and occurrence. Returns the file path.
        """
        if not self.statistics:
            raise ValueError("No statistics are available to write.")

        output_file = os.path.abspath(output_file)
        os.makedirs(os.path.dirname(output_file) or os.getcwd(), exist_ok=True)

        with open(output_file, "w") as out:
            out.write("# gQTEA All Bond Distance Analysis\n")
            out.write(f"# frames_used {self.num_frames}\n")
            out.write(f"# max_connection_distance {self.max_connection_distance:g}\n")
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
                "# row atom_i atom_j element_i element_j first_bonded_distance "
                "average variance standard_deviation occurrence_fraction frames_bonded\n"
            )
            for row_index, (pair, (average, variance, std_dev)) in enumerate(
                zip(self.connected_atom_pairs, self.statistics), start=1
            ):
                occurrence = pair.frames_bonded / self.num_frames
                out.write(
                    f"{row_index:>6d} "
                    f"{pair.atom_i + 1:>8d} {pair.atom_j + 1:>8d} "
                    f"{pair.element_i:>8s} {pair.element_j:>8s} "
                    f"{pair.first_bonded_distance:>16.8f} "
                    f"{average:>16.8f} {variance:>16.8f} {std_dev:>16.8f} "
                    f"{occurrence:>16.8f} {pair.frames_bonded:>10d}\n"
                )

        return output_file


class allBondAnalysisUI:
    """Toga frontend for all connected bond distance statistics."""

    def __init__(self, *args) -> None:
        self.trajec = None
        self.output_dir = os.getcwd()
        self.layout_main_window(*args)

    async def warning_function(self, title: str, message: str) -> None:
        await self.main_window.dialog(toga.InfoDialog(title, message))

    def layout_main_window(self, widget) -> None:
        self.main_window = toga.Window(
            title="All Bond Distance Analysis from Molecular Dynamics Simulations",
            size=(720, 560),
        )

        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style = Pack(margin=(5, 5), text_align=LEFT, width=240)
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=120)
        row_style = Pack(direction="row", margin=(0, 0, 5, 0))

        main_box = toga.Box(style=Pack(direction="column", margin=20))

        title_row = toga.Box(style=Pack(direction="row", margin=(0, 0, 10, 0)))
        title_box = toga.Box(style=Pack(width=660))
        title_label = toga.Label("All Bond Distance Analysis", style=heading_style)
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
            value="all_bond_analysis_combined.txt",
            placeholder="all_bond_analysis_combined.txt",
            style=input_style,
        )
        output_row.add(output_label)
        output_row.add(self.textInput_output)
        main_box.add(output_row)

        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(10, 0), font_size=12)
        )
        self.multi_line_text.value = (
            "This module re-evaluates connectivity every frame of an XYZ trajectory using "
            "the maximum connection distance, then computes the average distance, population "
            "variance, standard deviation, and occurrence (fraction of frames bonded) for each "
            "pair. Enter cell lattices a b c to apply the minimum-image convention (PBC); leave "
            "them blank for an isolated system. Provide solute atom indices (ranges allowed, "
            "e.g. 1-5 14-16 18 20) to exclude solvent atoms from the calculation."
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
            analyzer = AllBondAnalysis(self.trajec)
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
            self.max_connection_distance = AllBondAnalysis.DEFAULT_CONNECTION_DISTANCE
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
            analyzer = AllBondAnalysis(
                trajectory_file=self.trajec,
                max_connection_distance=self.max_connection_distance,
                solute_atom_indices=self.solute_atom_indices,
                cell_lengths=self.cell_lengths,
            )

            self.multi_line_text.value = (
                "Preparing all-bond analysis...\n"
                "Reading the first trajectory frame and checking the selected solute atoms."
            )
            await asyncio.sleep(0)

            analyzer.read_first_frame()

            self.multi_line_text.value = (
                f"First frame loaded ({analyzer.num_atoms} atoms).\n"
                "Re-evaluating connectivity every frame and computing the average distance, "
                "population variance, standard deviation, and occurrence for each pair.\n"
                "Please wait until the final summary appears."
            )
            await asyncio.sleep(0)

            await analyzer.compute_bond_statistics_async(
                progress_callback=None,
                ui_update_interval=500,
            )

            self.multi_line_text.value = (
                "Bond statistics completed.\n"
                "Writing the combined bond-analysis file."
            )
            await asyncio.sleep(0)

            output_file = analyzer.write_results(self.output_file)
        except Exception as exc:
            await self.warning_function("Error", f"All bond analysis failed: {exc}")
            return

        preview_pairs = []
        for row_index, pair in enumerate(analyzer.connected_atom_pairs[:10], start=1):
            occurrence = pair.frames_bonded / analyzer.num_frames
            preview_pairs.append(
                f"{row_index:>3d}: {pair.element_i}{pair.atom_i + 1}-"
                f"{pair.element_j}{pair.atom_j + 1} "
                f"first-bonded distance = {pair.first_bonded_distance:.6f} A, "
                f"occurrence = {occurrence:.1%}"
            )
        pair_preview = "\n".join(preview_pairs)
        if len(analyzer.connected_atom_pairs) > 10:
            pair_preview += "\n..."

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
            f"Connected pairs (bonded in >=1 frame): {len(analyzer.connected_atom_pairs)}\n"
            f"Maximum connection distance: {self.max_connection_distance:.6f} A\n\n"
            f"Output file:\n{output_file}\n\n"
            f"First mapped pairs:\n{pair_preview}"
        )

    def closeTopLevel(self, widget) -> None:
        self.main_window.close()
