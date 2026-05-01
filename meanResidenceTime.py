import asyncio
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Dict, Any

import numpy as np
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW


@dataclass
class MRTResults:
    xyz_file: str
    n_frames: int
    n_atoms: int
    dt: float
    dt_unit: str
    cutoff: float
    tolerance_frames: int
    reference_mode: str
    reference_definition: str
    observed_mode: str
    observed_definition: str
    n_observed_objects: int
    occupancy_mean: float
    n_continuous_events: int
    direct_continuous_mrt: float
    continuous_mrt: float
    intermittent_mrt: float
    time_axis: np.ndarray
    continuous_survival: np.ndarray
    intermittent_correlation: np.ndarray
    intermittent_relaxation: np.ndarray
    event_durations: np.ndarray
    object_labels: List[str]


class MeanResidenceTimeCalculator:
    """
    Computational engine for mean residence time (MRT) analysis from XYZ trajectories.

    Supported reference modes:
        - single_atom
        - geometric_center
        - center_of_mass

    Supported observed modes:
        - single_atom
        - atom_list_individual
        - group_list_geometric
        - group_list_com
    """

    ATOMIC_MASSES = {
        "H": 1.00784, "He": 4.002602,
        "Li": 6.94, "Be": 9.0121831, "B": 10.81, "C": 12.011, "N": 14.007,
        "O": 15.999, "F": 18.998403163, "Ne": 20.1797,
        "Na": 22.98976928, "Mg": 24.305, "Al": 26.9815385, "Si": 28.085,
        "P": 30.973761998, "S": 32.06, "Cl": 35.45, "Ar": 39.948,
        "K": 39.0983, "Ca": 40.078, "Sc": 44.955908, "Ti": 47.867,
        "V": 50.9415, "Cr": 51.9961, "Mn": 54.938044, "Fe": 55.845,
        "Co": 58.933194, "Ni": 58.6934, "Cu": 63.546, "Zn": 65.38,
        "Ga": 69.723, "Ge": 72.630, "As": 74.921595, "Se": 78.971,
        "Br": 79.904, "Kr": 83.798,
        "Rb": 85.4678, "Sr": 87.62, "Y": 88.90584, "Zr": 91.224,
        "Nb": 92.90637, "Mo": 95.95, "Tc": 98.0, "Ru": 101.07,
        "Rh": 102.90550, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.414,
        "In": 114.818, "Sn": 118.710, "Sb": 121.760, "Te": 127.60,
        "I": 126.90447, "Xe": 131.293,
        "Cs": 132.90545196, "Ba": 137.327, "La": 138.90547, "Ce": 140.116,
        "Pr": 140.90766, "Nd": 144.242, "Pm": 145.0, "Sm": 150.36,
        "Eu": 151.964, "Gd": 157.25, "Tb": 158.92535, "Dy": 162.500,
        "Ho": 164.93033, "Er": 167.259, "Tm": 168.93422, "Yb": 173.045,
        "Lu": 174.9668, "Hf": 178.49, "Ta": 180.94788, "W": 183.84,
        "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.084,
        "Au": 196.966569, "Hg": 200.592, "Tl": 204.38, "Pb": 207.2,
        "Bi": 208.98040
    }

    def __init__(self):
        pass

    @staticmethod
    def _clean_symbol(symbol: str) -> str:
        if not symbol:
            raise ValueError("Empty atomic symbol found in XYZ file.")
        symbol = symbol.strip()
        return symbol[0].upper() + symbol[1:].lower()

    @staticmethod
    def parse_atom_list(text: str) -> List[int]:
        """
        Parse 1-based atom indices from:
            '1,2,3'
            '1 2 3'
            '1, 2 3'
        Returns 0-based indices.
        """
        if not text or not text.strip():
            raise ValueError("Atom list is empty.")

        raw = text.replace(",", " ").split()
        atoms = []
        for item in raw:
            try:
                idx = int(item)
            except ValueError as exc:
                raise ValueError(f"Invalid atom index: '{item}'.") from exc
            if idx < 1:
                raise ValueError("Atom indices must start at 1.")
            atoms.append(idx - 1)

        if len(atoms) == 0:
            raise ValueError("No valid atom indices found.")
        return atoms

    @classmethod
    def parse_group_list(cls, text: str) -> List[List[int]]:
        """
        Parse group definitions using ';' between groups and ',' or whitespace within each group.

        Example:
            '1,2,3; 4,5,6; 10 11 12'
        """
        if not text or not text.strip():
            raise ValueError("Group list is empty.")

        groups = []
        for chunk in text.split(";"):
            chunk = chunk.strip()
            if not chunk:
                continue
            atoms = cls.parse_atom_list(chunk)
            groups.append(atoms)

        if len(groups) == 0:
            raise ValueError("No valid groups found.")
        return groups

    @staticmethod
    def _validate_indices(indices: List[int], n_atoms: int, label: str = "indices") -> None:
        for idx in indices:
            if idx < 0 or idx >= n_atoms:
                raise ValueError(
                    f"Atom index {idx + 1} in {label} is out of range. "
                    f"Trajectory has {n_atoms} atoms."
                )

    @classmethod
    def _validate_group_indices(cls, groups: List[List[int]], n_atoms: int, label: str = "groups") -> None:
        for i, group in enumerate(groups, start=1):
            if len(group) == 0:
                raise ValueError(f"Group {i} in {label} is empty.")
            cls._validate_indices(group, n_atoms, label=f"{label} group {i}")

    def read_xyz(self, filepath: str) -> Tuple[List[str], np.ndarray]:
        """
        Read a multi-frame XYZ trajectory.

        Returns
        -------
        symbols : list[str]
            Atomic symbols from the first frame.
        coords : np.ndarray
            Shape = (n_frames, n_atoms, 3)
        """
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"XYZ file not found: {filepath}")

        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        if not lines:
            raise ValueError("The XYZ file is empty.")

        frames = []
        symbols_ref = None
        line_idx = 0
        n_lines = len(lines)

        while line_idx < n_lines:
            line = lines[line_idx].strip()
            if not line:
                line_idx += 1
                continue

            try:
                n_atoms = int(line)
            except ValueError as exc:
                raise ValueError(
                    f"Expected number of atoms at line {line_idx + 1}, found: '{lines[line_idx]}'"
                ) from exc

            if line_idx + 1 >= n_lines:
                raise ValueError("Incomplete XYZ frame: missing comment line.")

            comment_line = lines[line_idx + 1]  # noqa: F841
            start = line_idx + 2
            end = start + n_atoms

            if end > n_lines:
                raise ValueError(
                    f"Incomplete XYZ frame starting at line {line_idx + 1}: "
                    f"expected {n_atoms} atomic lines, but file ended early."
                )

            frame_symbols = []
            frame_coords = np.zeros((n_atoms, 3), dtype=float)

            for i, row in enumerate(lines[start:end]):
                parts = row.split()
                if len(parts) < 4:
                    raise ValueError(
                        f"Invalid XYZ coordinate line at line {start + i + 1}: '{row}'"
                    )

                symbol = self._clean_symbol(parts[0])
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid coordinates at line {start + i + 1}: '{row}'"
                    ) from exc

                frame_symbols.append(symbol)
                frame_coords[i] = [x, y, z]

            if symbols_ref is None:
                symbols_ref = frame_symbols
            else:
                if len(frame_symbols) != len(symbols_ref):
                    raise ValueError("The number of atoms changes between frames.")
                if frame_symbols != symbols_ref:
                    raise ValueError(
                        "Atomic ordering or symbols change between frames. "
                        "This script requires a consistent atom order through the trajectory."
                    )

            frames.append(frame_coords)
            line_idx = end

        if not frames:
            raise ValueError("No valid XYZ frames were read.")

        coords = np.array(frames, dtype=float)
        return symbols_ref, coords

    def get_masses(self, symbols: List[str]) -> np.ndarray:
        masses = []
        for sym in symbols:
            if sym not in self.ATOMIC_MASSES:
                raise ValueError(
                    f"Atomic mass for element '{sym}' is not available in the internal mass table."
                )
            masses.append(self.ATOMIC_MASSES[sym])
        return np.array(masses, dtype=float)

    @staticmethod
    def geometric_center(frame_coords: np.ndarray, indices: List[int]) -> np.ndarray:
        return np.mean(frame_coords[indices], axis=0)

    @staticmethod
    def center_of_mass(frame_coords: np.ndarray, indices: List[int], masses: np.ndarray) -> np.ndarray:
        local_masses = masses[indices]
        weighted = frame_coords[indices] * local_masses[:, None]
        return np.sum(weighted, axis=0) / np.sum(local_masses)

    def build_reference_positions(
        self,
        coords: np.ndarray,
        masses: np.ndarray,
        mode: str,
        definition: str
    ) -> np.ndarray:
        n_frames, n_atoms, _ = coords.shape

        if mode == "single_atom":
            atom = self.parse_atom_list(definition)
            if len(atom) != 1:
                raise ValueError("Reference mode 'single_atom' requires exactly one atom index.")
            self._validate_indices(atom, n_atoms, label="reference atom")
            idx = atom[0]
            return coords[:, idx, :]

        atoms = self.parse_atom_list(definition)
        self._validate_indices(atoms, n_atoms, label="reference group")

        ref_positions = np.zeros((n_frames, 3), dtype=float)

        if mode == "geometric_center":
            for i in range(n_frames):
                ref_positions[i] = self.geometric_center(coords[i], atoms)
            return ref_positions

        if mode == "center_of_mass":
            for i in range(n_frames):
                ref_positions[i] = self.center_of_mass(coords[i], atoms, masses)
            return ref_positions

        raise ValueError(f"Unknown reference mode: {mode}")

    def build_observed_positions(
        self,
        coords: np.ndarray,
        masses: np.ndarray,
        mode: str,
        definition: str
    ) -> Tuple[np.ndarray, List[str]]:
        """
        Returns
        -------
        obs_positions : np.ndarray
            Shape = (n_objects, n_frames, 3)
        labels : list[str]
        """
        n_frames, n_atoms, _ = coords.shape

        if mode == "single_atom":
            atoms = self.parse_atom_list(definition)
            if len(atoms) != 1:
                raise ValueError("Observed mode 'single_atom' requires exactly one atom index.")
            self._validate_indices(atoms, n_atoms, label="observed atom")
            idx = atoms[0]
            return coords[:, idx, :][None, :, :], [f"Atom {idx + 1}"]

        if mode == "atom_list_individual":
            atoms = self.parse_atom_list(definition)
            self._validate_indices(atoms, n_atoms, label="observed atom list")

            positions = np.zeros((len(atoms), n_frames, 3), dtype=float)
            labels = []
            for k, idx in enumerate(atoms):
                positions[k] = coords[:, idx, :]
                labels.append(f"Atom {idx + 1}")
            return positions, labels

        if mode in ("group_list_geometric", "group_list_com"):
            groups = self.parse_group_list(definition)
            self._validate_group_indices(groups, n_atoms, label="observed groups")

            positions = np.zeros((len(groups), n_frames, 3), dtype=float)
            labels = []

            for g, group in enumerate(groups):
                labels.append("Group " + str(g + 1) + " (" + ",".join(str(i + 1) for i in group) + ")")
                for i in range(n_frames):
                    if mode == "group_list_geometric":
                        positions[g, i] = self.geometric_center(coords[i], group)
                    else:
                        positions[g, i] = self.center_of_mass(coords[i], group, masses)

            return positions, labels

        raise ValueError(f"Unknown observed mode: {mode}")

    @staticmethod
    def build_occupancy(
        reference_positions: np.ndarray,
        observed_positions: np.ndarray,
        cutoff: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        reference_positions: (n_frames, 3)
        observed_positions:  (n_objects, n_frames, 3)

        Returns
        -------
        distances : np.ndarray, shape (n_objects, n_frames)
        occupancy : np.ndarray, shape (n_objects, n_frames), dtype=int
        """
        if cutoff <= 0:
            raise ValueError("Cutoff radius must be positive.")

        # Broadcasting:
        # observed_positions -> (n_objects, n_frames, 3)
        # reference_positions -> (1, n_frames, 3)
        diffs = observed_positions - reference_positions[None, :, :]
        distances = np.linalg.norm(diffs, axis=2)
        occupancy = (distances <= cutoff).astype(int)
        return distances, occupancy

    @staticmethod
    def apply_tolerance_to_occupancy(occupancy: np.ndarray, tolerance_frames: int) -> np.ndarray:
        """
        Fill short zero-gaps of length <= tolerance_frames between runs of ones.
        """
        if tolerance_frames <= 0:
            return occupancy.copy()

        corrected = occupancy.copy()
        n_objects, n_frames = corrected.shape

        for obj in range(n_objects):
            arr = corrected[obj]
            i = 0
            while i < n_frames:
                if arr[i] == 1:
                    i += 1
                    continue

                start_zero = i
                while i < n_frames and arr[i] == 0:
                    i += 1
                end_zero = i - 1
                zero_len = end_zero - start_zero + 1

                left_is_one = (start_zero - 1 >= 0 and arr[start_zero - 1] == 1)
                right_is_one = (i < n_frames and arr[i] == 1)

                if left_is_one and right_is_one and zero_len <= tolerance_frames:
                    arr[start_zero:end_zero + 1] = 1

        return corrected

    @staticmethod
    def extract_event_durations(occupancy: np.ndarray, dt: float) -> np.ndarray:
        durations = []
        n_objects, n_frames = occupancy.shape

        for obj in range(n_objects):
            arr = occupancy[obj]
            i = 0
            while i < n_frames:
                if arr[i] == 1:
                    start = i
                    while i < n_frames and arr[i] == 1:
                        i += 1
                    end = i - 1
                    n_block = end - start + 1
                    durations.append(n_block * dt)
                else:
                    i += 1

        if len(durations) == 0:
            return np.array([], dtype=float)
        return np.array(durations, dtype=float)

    @staticmethod
    def compute_continuous_survival_from_events(
        event_durations: np.ndarray,
        dt: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        S(t_l) = fraction of event durations >= t_l
        """
        if dt <= 0:
            raise ValueError("Time between saved frames must be positive.")

        if len(event_durations) == 0:
            time_axis = np.array([0.0], dtype=float)
            survival = np.array([0.0], dtype=float)
            return time_axis, survival

        max_duration = float(np.max(event_durations))
        n_lags = int(np.floor(max_duration / dt)) + 1
        time_axis = np.arange(n_lags, dtype=float) * dt
        survival = np.zeros_like(time_axis)

        n_events = len(event_durations)
        for lag, t in enumerate(time_axis):
            survival[lag] = np.sum(event_durations >= t) / n_events

        return time_axis, survival

    @staticmethod
    def compute_intermittent_functions(
        occupancy: np.ndarray,
        dt: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        """
        Computes:
            C(t) = <h(0)h(t)> / <h>
            R(t) = (C(t) - <h>) / (1 - <h>)

        Returns
        -------
        time_axis, C, R, h_mean
        """
        if dt <= 0:
            raise ValueError("Time between saved frames must be positive.")

        n_objects, n_frames = occupancy.shape
        h_mean = float(np.mean(occupancy))

        time_axis = np.arange(n_frames, dtype=float) * dt
        C = np.zeros(n_frames, dtype=float)
        R = np.zeros(n_frames, dtype=float)

        for lag in range(n_frames):
            valid_origins = occupancy[:, :n_frames - lag]
            shifted = occupancy[:, lag:]

            denominator = np.sum(valid_origins)
            if denominator == 0:
                C[lag] = 0.0
            else:
                numerator = np.sum(valid_origins * shifted)
                C[lag] = numerator / denominator

        # If h_mean is 1.0, occupancy never changes. Then R(t) is undefined by the usual formula.
        # In that special case, set R(t)=1 for all t.
        if np.isclose(h_mean, 1.0):
            R[:] = 1.0
        else:
            R = (C - h_mean) / (1.0 - h_mean)

        # Numerical cleanup
        C = np.clip(C, 0.0, 1.0)
        R = np.clip(R, -1.0, 1.0)

        return time_axis, C, R, h_mean

    @staticmethod
    def discrete_integral(y: np.ndarray, dt: float) -> float:
        if len(y) == 0:
            return 0.0
        return float(np.sum(y) * dt)

    def run(
        self,
        xyz_file: str,
        dt: float,
        dt_unit: str,
        cutoff: float,
        tolerance_frames: int,
        reference_mode: str,
        reference_definition: str,
        observed_mode: str,
        observed_definition: str
    ) -> MRTResults:
        symbols, coords = self.read_xyz(xyz_file)
        masses = self.get_masses(symbols)

        n_frames, n_atoms, _ = coords.shape

        reference_positions = self.build_reference_positions(
            coords=coords,
            masses=masses,
            mode=reference_mode,
            definition=reference_definition
        )

        observed_positions, object_labels = self.build_observed_positions(
            coords=coords,
            masses=masses,
            mode=observed_mode,
            definition=observed_definition
        )

        _, occupancy = self.build_occupancy(
            reference_positions=reference_positions,
            observed_positions=observed_positions,
            cutoff=cutoff
        )

        corrected_occupancy = self.apply_tolerance_to_occupancy(
            occupancy=occupancy,
            tolerance_frames=tolerance_frames
        )

        event_durations = self.extract_event_durations(corrected_occupancy, dt=dt)
        t_surv, S = self.compute_continuous_survival_from_events(event_durations, dt=dt)
        t_corr, C, R, h_mean = self.compute_intermittent_functions(corrected_occupancy, dt=dt)

        # Use the full time axis for intermittent integration.
        # For continuous, integrate the survival built from events.
        direct_cont_mrt = float(np.mean(event_durations)) if len(event_durations) > 0 else 0.0
        continuous_mrt = self.discrete_integral(S, dt)
        intermittent_mrt = self.discrete_integral(R, dt)

        return MRTResults(
            xyz_file=xyz_file,
            n_frames=n_frames,
            n_atoms=n_atoms,
            dt=dt,
            dt_unit=dt_unit,
            cutoff=cutoff,
            tolerance_frames=tolerance_frames,
            reference_mode=reference_mode,
            reference_definition=reference_definition,
            observed_mode=observed_mode,
            observed_definition=observed_definition,
            n_observed_objects=corrected_occupancy.shape[0],
            occupancy_mean=h_mean,
            n_continuous_events=len(event_durations),
            direct_continuous_mrt=direct_cont_mrt,
            continuous_mrt=continuous_mrt,
            intermittent_mrt=intermittent_mrt,
            time_axis=t_corr,
            continuous_survival=S,
            intermittent_correlation=C,
            intermittent_relaxation=R,
            event_durations=event_durations,
            object_labels=object_labels
        )

    @staticmethod
    def export_curves(results: MRTResults, output_dir: str | None = None) -> Dict[str, str]:
        """
        Export data curves to .dat files in the same folder as the XYZ file unless output_dir is given.
        """
        xyz_path = Path(results.xyz_file)
        base_dir = Path(output_dir) if output_dir else xyz_path.parent
        base_name = xyz_path.stem

        files = {}

        continuous_file = base_dir / f"{base_name}_continuous_survival.dat"
        with continuous_file.open("w", encoding="utf-8") as f:
            f.write("# time  S(t)\n")
            dt = results.dt
            for i, value in enumerate(results.continuous_survival):
                f.write(f"{i * dt:20.10f} {value:20.10f}\n")
        files["continuous_survival"] = str(continuous_file)

        intermittent_file = base_dir / f"{base_name}_intermittent_correlation.dat"
        with intermittent_file.open("w", encoding="utf-8") as f:
            f.write("# time  C(t)  R(t)\n")
            for i, (cval, rval) in enumerate(
                zip(results.intermittent_correlation, results.intermittent_relaxation)
            ):
                f.write(f"{i * results.dt:20.10f} {cval:20.10f} {rval:20.10f}\n")
        files["intermittent_correlation"] = str(intermittent_file)

        durations_file = base_dir / f"{base_name}_mrt_event_durations.dat"
        with durations_file.open("w", encoding="utf-8") as f:
            f.write("# event_index  duration\n")
            for i, duration in enumerate(results.event_durations, start=1):
                f.write(f"{i:10d} {duration:20.10f}\n")
        files["event_durations"] = str(durations_file)

        summary_file = base_dir / f"{base_name}_mrt_summary.txt"
        with summary_file.open("w", encoding="utf-8") as f:
            f.write(MeanResidenceTimeCalculator.format_results_text(results))
        files["summary"] = str(summary_file)

        return files

    @staticmethod
    def format_results_text(results: MRTResults) -> str:
        lines = [
            "Mean Residence Time (MRT) Analysis",
            "=" * 60,
            f"XYZ file: {results.xyz_file}",
            f"Number of frames: {results.n_frames}",
            f"Number of atoms: {results.n_atoms}",
            f"Time between saved frames: {results.dt} {results.dt_unit}",
            f"Cutoff radius: {results.cutoff}",
            f"Tolerance frames: {results.tolerance_frames}",
            "",
            f"Reference mode: {results.reference_mode}",
            f"Reference definition: {results.reference_definition}",
            f"Observed mode: {results.observed_mode}",
            f"Observed definition: {results.observed_definition}",
            f"Number of observed objects: {results.n_observed_objects}",
            "",
            f"Average occupancy <h>: {results.occupancy_mean:.8f}",
            f"Number of continuous residence events: {results.n_continuous_events}",
            f"Direct continuous MRT (mean event duration): {results.direct_continuous_mrt:.8f} {results.dt_unit}",
            f"Continuous MRT from S(t): {results.continuous_mrt:.8f} {results.dt_unit}",
            f"Intermittent MRT from R(t): {results.intermittent_mrt:.8f} {results.dt_unit}",
            "",
            "Observed objects:",
        ]

        for label in results.object_labels:
            lines.append(f"  - {label}")

        lines.extend([
            "",
            "Notes:",
            "1. Distances are computed directly from XYZ coordinates.",
            "2. Standard XYZ files do not contain periodic box vectors.",
            "3. Therefore, this script does not apply minimum-image periodic boundary corrections.",
            "4. If your trajectory is periodic, a future extension should read the box and apply PBC-aware distances."
        ])
        return "\n".join(lines)


class MeanResidenceTimeUI:
    """
    Toga user interface for mean residence time calculations.

    Typical integration:
        MeanResidenceTimeUI(self.main_window)
    """

    REFERENCE_MODES = [
        "single_atom",
        "geometric_center",
        "center_of_mass",
    ]

    OBSERVED_MODES = [
        "single_atom",
        "atom_list_individual",
        "group_list_geometric",
        "group_list_com",
    ]

    def __init__(self, main_window):
        self.main_window = main_window
        self.app = main_window.app if hasattr(main_window, "app") else None
        self.calculator = MeanResidenceTimeCalculator()
        self.results = None

        self.window = toga.Window(
            title="Mean Residence Time (MRT)",
            size=(1100, 400)
        )

        self._build_ui()
        self.window.show()    
    
    def _build_ui(self):
        # File selection
        self.xyz_path_input = toga.TextInput(
            placeholder="Select the XYZ trajectory file...",
            readonly=False,
            style=Pack(flex=1, margin=5)
        )

        self.browse_button = toga.Button(
            "Browse XYZ",
            on_press=self.browse_xyz_file,
            style=Pack(width=140, margin=5)
        )

        file_box = toga.Box(
            children=[
                toga.Label("XYZ trajectory file:", style=Pack(width=150, margin=8)),
                self.xyz_path_input,
                self.browse_button
            ],
            style=Pack(direction=ROW, margin_bottom=5)
        )

        # Numerical controls
        self.dt_input = toga.TextInput(
            value="1.0",
            placeholder="e.g. 1.0",
            style=Pack(width=120, margin=5)
        )

        self.dt_unit_selection = toga.Selection(
            items=["fs", "ps", "ns", "a.u."],
            style=Pack(width=100, margin=5)
        )
        self.dt_unit_selection.value = "fs"

        self.cutoff_input = toga.TextInput(
            value="3.5",
            placeholder="Cutoff radius",
            style=Pack(width=120, margin=5)
        )

        self.tolerance_input = toga.TextInput(
            value="0",
            placeholder="Tolerance frames",
            style=Pack(width=120, margin=5)
        )

        numeric_box = toga.Box(
            children=[
                toga.Label("Δt between saved frames:", style=Pack(width=170, margin=8)),
                self.dt_input,
                self.dt_unit_selection,
                toga.Label("Cutoff radius (Å):", style=Pack(width=130, margin=8)),
                self.cutoff_input,
                toga.Label("Tolerance frames:", style=Pack(width=120, margin=8)),
                self.tolerance_input,
            ],
            style=Pack(direction=ROW, margin_bottom=5)
        )

        # Reference controls
        self.reference_mode_selection = toga.Selection(
            items=self.REFERENCE_MODES,
            style=Pack(width=220, margin=5)
        )
        self.reference_mode_selection.value = self.REFERENCE_MODES[0]

        self.reference_definition_input = toga.TextInput(
            placeholder="Reference definition",
            style=Pack(flex=1, margin=5)
        )

        reference_box = toga.Box(
            children=[
                toga.Label("Reference mode:", style=Pack(width=120, margin=8)),
                self.reference_mode_selection,
            ],
            style=Pack(direction=ROW, margin_bottom=3)
        )

        reference_def_box = toga.Box(
            children=[
                toga.Label("Reference definition:", style=Pack(width=140, margin=8)),
                self.reference_definition_input,
            ],
            style=Pack(direction=ROW, margin_bottom=8)
        )

        # Observed controls
        self.observed_mode_selection = toga.Selection(
            items=self.OBSERVED_MODES,
            style=Pack(width=220, margin=5)
        )
        self.observed_mode_selection.value = self.OBSERVED_MODES[1]

        self.observed_definition_input = toga.TextInput(
            placeholder="Observed definition",
            style=Pack(flex=1, margin=5)
        )

        observed_box = toga.Box(
            children=[
                toga.Label("Observed mode:", style=Pack(width=120, margin=8)),
                self.observed_mode_selection,
            ],
            style=Pack(direction=ROW, margin_bottom=3)
        )

        observed_def_box = toga.Box(
            children=[
                toga.Label("Observed definition:", style=Pack(width=140, margin=8)),
                self.observed_definition_input,
            ],
            style=Pack(direction=ROW, margin_bottom=8)
        )

        # Buttons
        self.run_button = toga.Button(
            "Run MRT",
            on_press=self.run_mrt_calculation,
            style=Pack(width=120, margin=5)
        )

        self.export_button = toga.Button(
            "Export Results",
            on_press=self.export_results,
            style=Pack(width=140, margin=5)
        )

        self.clear_button = toga.Button(
            "Clear Output",
            on_press=self.clear_output,
            style=Pack(width=120, margin=5)
        )

        self.status_label = toga.Label(
            "Ready.",
            style=Pack(margin=8)
        )

        self.progress_bar = toga.ProgressBar(
            max=100,
            style=Pack(margin_left=5, margin_right=5, margin_bottom=5)
        )
        self.progress_bar.value = 0

        button_box = toga.Box(
            children=[
                self.run_button,
                self.export_button,
                self.clear_button,
                self.status_label
            ],
            style=Pack(direction=ROW, margin_bottom=5)
        )

        main_box = toga.Box(
            children=[
                file_box,
                numeric_box,
                toga.Divider(style=Pack(margin_top=5, margin_bottom=5)),
                reference_box,
                reference_def_box,
                observed_box,
                observed_def_box,
                button_box,
                self.progress_bar,
            ],
            style=Pack(direction=COLUMN, margin=10, flex=1)
        )

        self.window.content = main_box

    async def show_window(self):
        self.window.show()

    def update_help_texts(self, widget):
        return

    async def browse_xyz_file(self, widget):
        try:
            file_path = await self.window.dialog(
                toga.OpenFileDialog(
                    title="Select XYZ trajectory file",
                    file_types=["xyz"]
                )
            )
            if file_path:
                self.xyz_path_input.value = str(file_path)
                self.status_label.text = "XYZ file selected."
        except Exception as exc:
            await self.window.error_dialog("File Selection Error", str(exc))

    def _collect_inputs(self) -> Dict[str, Any]:
        xyz_file = self.xyz_path_input.value.strip() if self.xyz_path_input.value else ""
        if not xyz_file:
            raise ValueError("Please select an XYZ trajectory file.")

        try:
            dt = float(self.dt_input.value.strip())
        except Exception as exc:
            raise ValueError("The time between saved frames must be a valid number.") from exc

        if dt <= 0:
            raise ValueError("The time between saved frames must be positive.")

        dt_unit = str(self.dt_unit_selection.value).strip()

        try:
            cutoff = float(self.cutoff_input.value.strip())
        except Exception as exc:
            raise ValueError("The cutoff radius must be a valid number.") from exc

        if cutoff <= 0:
            raise ValueError("The cutoff radius must be positive.")

        try:
            tolerance_frames = int(self.tolerance_input.value.strip())
        except Exception as exc:
            raise ValueError("Tolerance frames must be an integer.") from exc

        if tolerance_frames < 0:
            raise ValueError("Tolerance frames cannot be negative.")

        reference_mode = str(self.reference_mode_selection.value).strip()
        observed_mode = str(self.observed_mode_selection.value).strip()

        reference_definition = (
            self.reference_definition_input.value.strip()
            if self.reference_definition_input.value else ""
        )
        observed_definition = (
            self.observed_definition_input.value.strip()
            if self.observed_definition_input.value else ""
        )

        if not reference_definition:
            raise ValueError("Please provide the reference definition.")
        if not observed_definition:
            raise ValueError("Please provide the observed definition.")

        return {
            "xyz_file": xyz_file,
            "dt": dt,
            "dt_unit": dt_unit,
            "cutoff": cutoff,
            "tolerance_frames": tolerance_frames,
            "reference_mode": reference_mode,
            "reference_definition": reference_definition,
            "observed_mode": observed_mode,
            "observed_definition": observed_definition
        }

    async def run_mrt_calculation(self, widget):
        try:
            params = self._collect_inputs()
        except Exception as exc:
            await self.window.error_dialog("Input Error", str(exc))
            return

        self.status_label.text = "Reading trajectory and computing MRT..."
        self.progress_bar.value = 0
        self.progress_bar.start()
        self.run_button.enabled = False
        self.export_button.enabled = False
        self.clear_button.enabled = False
        await asyncio.sleep(0)

        try:
            worker = asyncio.create_task(asyncio.to_thread(self.calculator.run, **params))
            progress_value = 3.0
            progress_direction = 1.0

            while not worker.done():
                self.progress_bar.value = progress_value
                progress_value += progress_direction * 3.0

                if progress_value >= 92.0:
                    progress_value = 92.0
                    progress_direction = -1.0
                elif progress_value <= 8.0:
                    progress_value = 8.0
                    progress_direction = 1.0

                await asyncio.sleep(0.08)

            results = await worker
            self.results = results
            self.status_label.text = (
                "Done. "
                f"Objects: {results.n_observed_objects} | "
                f"Direct MRT: {results.direct_continuous_mrt:.4f} {results.dt_unit} | "
                f"Continuous MRT: {results.continuous_mrt:.4f} {results.dt_unit}"
            )
            self.progress_bar.value = 100
        except Exception as exc:
            self.status_label.text = "Calculation failed."
            self.progress_bar.value = 0
            await self.window.error_dialog("MRT Calculation Error", str(exc))
        finally:
            self.progress_bar.stop()
            self.run_button.enabled = True
            self.export_button.enabled = True
            self.clear_button.enabled = True

    async def export_results(self, widget):
        if self.results is None:
            await self.window.info_dialog("No Results", "Please run the MRT calculation first.")
            return

        try:
            files = self.calculator.export_curves(self.results)
            msg = (
                "Results exported successfully:\n\n"
                f"Summary: {files['summary']}\n"
                f"Continuous survival: {files['continuous_survival']}\n"
                f"Intermittent correlation: {files['intermittent_correlation']}\n"
                f"Event durations: {files['event_durations']}"
            )
            self.status_label.text = "Results exported."
            await self.window.info_dialog("Export Successful", msg)
        except Exception as exc:
            await self.window.error_dialog("Export Error", str(exc))

    def clear_output(self, widget):
        self.status_label.text = "Status cleared."
        self.progress_bar.value = 0


# Example integration helper
# You can adapt this to the style used in your main gqteaWinToga launcher.
async def open_mean_residence_time_window(app, main_window):
    mrt_ui = MeanResidenceTimeUI(app=app, main_window=main_window)
    await mrt_ui.show_window()
