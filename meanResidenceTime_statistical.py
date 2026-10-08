import asyncio
import os
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Dict, Any

import numpy as np
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

from allBondAnalysis import parse_solute_index_ranges
from displayPlots import DisplayPlots
from help import AtomicData, HelpGqteaWin


@dataclass
class MRTResults:
    xyz_file: str
    n_frames: int
    n_atoms: int
    dt: float
    dt_unit: str
    cutoff: float
    tolerance_frames: int
    cell_lengths: Tuple[float, float, float] | None
    reference_mode: str
    reference_definition: str
    observed_mode: str
    observed_definition: str
    n_observed_objects: int
    occupancy_mean: float

    # --- continuous statistics -------------------------------------------
    # These are two DIFFERENT quantities and coincide only for exponential
    # kinetics; see the module docstring of format_results_text().
    n_continuous_events: int
    n_censored_events: int
    mean_residence_time: float           # <T>, mean duration of complete events
    mean_residence_uncertainty: float
    imm_survival_time: float             # integral of the origin-averaged S(t)
    event_durations: np.ndarray
    survival_time_axis: np.ndarray
    origin_survival: np.ndarray

    # --- intermittent statistics -----------------------------------------
    intermittent_mrt: float
    intermittent_lags_used: int
    integration_rule: str
    max_lag: int
    time_axis: np.ndarray
    intermittent_correlation: np.ndarray
    intermittent_relaxation: np.ndarray

    object_labels: List[str]
    tolerance_scan: List[Dict[str, Any]]


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

    # help.py's AtomicData is the project's reference-data hub. The module used
    # to carry its own 83-element copy of the periodic table; this covers all
    # 118 (so actinides no longer raise) and there is one table to maintain.
    ATOMIC_MASSES = AtomicData.atomic_masses

    def __init__(self):
        pass

    @staticmethod
    def resolve_output_dir(xyz_file: str, chosen: str | None) -> str:
        """Where the .dat/.txt files go: the chosen folder, else the input's."""
        if chosen and chosen.strip():
            candidate = chosen.strip()
            if not os.path.isdir(candidate):
                raise ValueError(f"Output directory does not exist: {candidate}")
            return candidate
        return str(Path(xyz_file).parent)

    @staticmethod
    def _clean_symbol(symbol: str) -> str:
        if not symbol:
            raise ValueError("Empty atomic symbol found in XYZ file.")
        symbol = symbol.strip()
        return symbol[0].upper() + symbol[1:].lower()

    @staticmethod
    def parse_atom_list(text: str) -> List[int]:
        """
        Parse 1-based atom indices, with or without ranges:
            '1,2,3'
            '1 2 3'
            '1, 2 3'
            '3-7,10,15-17'      ->  3,4,5,6,7,10,15,16,17
            '10,3-5,1'          ->  10,3,4,5,1

        Ranges are inclusive. Input order is preserved (the tool labels the
        observed objects in the order they were typed), so the expansion is
        not sorted or de-duplicated.

        Returns 0-based indices.
        """
        if not text or not text.strip():
            raise ValueError("Atom list is empty.")

        # Shared with the all-* connectivity tools: one tested implementation
        # of the range syntax for the whole project.
        atoms = parse_solute_index_ranges(text, sort_unique=False)

        if len(atoms) == 0:
            raise ValueError("No valid atom indices found.")
        return [idx - 1 for idx in atoms]

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
    def parse_cell_lengths(text: str | None) -> Tuple[float, float, float] | None:
        """Parse 'a b c' (or 'a,b,c') orthorhombic cell edges in Angstrom.

        An empty or whitespace-only string means "no periodic boundary", the
        same convention bond.py and radialDistribution.py use.
        """
        if text is None or not text.strip():
            return None

        parts = text.replace(",", " ").split()
        if len(parts) != 3:
            raise ValueError(
                "Cell lengths must be three numbers 'a b c', or blank for no "
                "periodic boundary."
            )

        try:
            a, b, c = (float(p) for p in parts)
        except ValueError as exc:
            raise ValueError(f"Invalid cell lengths: '{text}'.") from exc

        if a <= 0 or b <= 0 or c <= 0:
            raise ValueError("Cell lengths must all be positive.")
        return (a, b, c)

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

        # Vectorised over frames: the per-frame Python loop this replaces was
        # ~86x slower on a 3000-frame trajectory, with identical results.
        if mode == "geometric_center":
            return self._geometric_center_series(coords, atoms)

        if mode == "center_of_mass":
            return self._center_of_mass_series(coords, atoms, masses)

        raise ValueError(f"Unknown reference mode: {mode}")

    @staticmethod
    def _geometric_center_series(coords: np.ndarray, indices: List[int]) -> np.ndarray:
        """Geometric centre of ``indices`` for every frame -> (n_frames, 3)."""
        return coords[:, indices, :].mean(axis=1)

    @staticmethod
    def _center_of_mass_series(
        coords: np.ndarray, indices: List[int], masses: np.ndarray
    ) -> np.ndarray:
        """Mass-weighted centre of ``indices`` for every frame -> (n_frames, 3)."""
        local = masses[indices]
        total = float(local.sum())
        if total <= 0:
            raise ValueError(
                "Total mass of the selected atoms is zero; check the selection."
            )
        return np.tensordot(coords[:, indices, :], local, axes=([1], [0])) / total

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

            geometric = (mode == "group_list_geometric")
            positions = np.stack([
                self._geometric_center_series(coords, group) if geometric
                else self._center_of_mass_series(coords, group, masses)
                for group in groups
            ])
            labels = [
                "Group " + str(g + 1) + " (" + ",".join(str(i + 1) for i in group) + ")"
                for g, group in enumerate(groups)
            ]

            return positions, labels

        raise ValueError(f"Unknown observed mode: {mode}")

    @staticmethod
    def build_occupancy(
        reference_positions: np.ndarray,
        observed_positions: np.ndarray,
        cutoff: float,
        cell_lengths: Tuple[float, float, float] | None = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        reference_positions: (n_frames, 3)
        observed_positions:  (n_objects, n_frames, 3)
        cell_lengths:        orthorhombic (a, b, c), or None for no PBC.

        With a cell, distances use the orthorhombic minimum-image convention
        ``delta -= box * round(delta / box)``, which is exact only while the
        cutoff stays within half the shortest edge; that bound is enforced
        here.  Coordinates need not be wrapped into the box -- minimum image
        tolerates unwrapped input.

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

        if cell_lengths is not None:
            box = np.asarray(cell_lengths, dtype=float)
            if box.shape != (3,):
                raise ValueError("Cell lengths must be three numbers 'a b c'.")
            if np.any(box <= 0):
                raise ValueError("Cell lengths must all be positive.")
            half_shortest = float(np.min(box)) / 2.0
            if cutoff > half_shortest:
                raise ValueError(
                    f"Cutoff radius {cutoff:.4f} exceeds half the shortest cell "
                    f"edge ({half_shortest:.4f}). The minimum-image convention is "
                    f"only valid up to that bound; reduce the cutoff or use a "
                    f"larger cell."
                )
            diffs = diffs - box * np.round(diffs / box)

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
    def _runs(occupancy: np.ndarray) -> List[Tuple[int, int, int]]:
        """Every maximal run of ones, as (object, start, length)."""
        runs: List[Tuple[int, int, int]] = []
        n_objects, n_frames = occupancy.shape

        for obj in range(n_objects):
            arr = occupancy[obj]
            i = 0
            while i < n_frames:
                if arr[i] == 1:
                    start = i
                    while i < n_frames and arr[i] == 1:
                        i += 1
                    runs.append((obj, start, i - start))
                else:
                    i += 1
        return runs

    @classmethod
    def extract_event_durations(
        cls,
        occupancy: np.ndarray,
        dt: float,
        censor_boundary: bool = False
    ) -> np.ndarray:
        """Durations of residence events.

        A run touching frame 0 or the final frame was already in progress when
        the trajectory started, or had not finished when it ended.  Recording
        it at its truncated length biases the mean residence time downwards,
        and biases it worst for the longest-lived species.  With
        ``censor_boundary=True`` such runs are dropped (right-censoring);
        ``count_boundary_events`` reports how many.
        """
        n_frames = occupancy.shape[1]
        durations = [
            length * dt
            for _, start, length in cls._runs(occupancy)
            if not (censor_boundary and (start == 0 or start + length == n_frames))
        ]
        return np.array(durations, dtype=float)

    @classmethod
    def count_boundary_events(cls, occupancy: np.ndarray) -> int:
        """Number of runs clipped by the start or the end of the trajectory."""
        n_frames = occupancy.shape[1]
        return sum(
            1
            for _, start, length in cls._runs(occupancy)
            if start == 0 or start + length == n_frames
        )

    @classmethod
    def compute_origin_averaged_survival(
        cls,
        occupancy: np.ndarray,
        dt: float,
        max_lag: int | None = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Continuous survival function in the Impey-Madden-McDonald sense.

        S(t) is averaged over every *time origin* at which an object is inside
        the shell -- not over events -- so a long event contributes as many
        origins as it has frames:

            S(l) = sum_i max(0, n_i - l) / sum_i n_i

        for runs of n_i frames.  That closed form is exactly equivalent to
        scanning all origins explicitly (and is what the unit tests check
        against), but costs O(events) per lag instead of O(frames^2).

        Its integral is the mean *residual* time, <T^2>/(2<T>), which equals
        the mean residence time <T> only for exponential kinetics.
        """
        if dt <= 0:
            raise ValueError("Time between saved frames must be positive.")

        lengths = np.array([length for _, _, length in cls._runs(occupancy)], dtype=float)

        if lengths.size == 0:
            n_lags = max(1, int(max_lag) if max_lag else 1)
            return np.arange(n_lags, dtype=float) * dt, np.zeros(n_lags, dtype=float)

        if max_lag is None:
            n_lags = int(lengths.max())
        else:
            n_lags = max(1, min(int(max_lag), occupancy.shape[1]))

        lags = np.arange(n_lags, dtype=float)
        # survivors[l] = sum_i max(0, n_i - l)
        survivors = np.maximum(0.0, lengths[:, None] - lags[None, :]).sum(axis=0)
        survival = survivors / lengths.sum()

        return lags * dt, survival

    @staticmethod
    def default_max_lag(n_frames: int) -> int:
        """Lag ceiling for correlation functions: a tenth of the trajectory.

        The longest lags average over a handful of time origins and are pure
        noise; integrating them into a residence time is meaningless.
        """
        return max(1, min(n_frames, n_frames // 10))

    @classmethod
    def compute_intermittent_functions(
        cls,
        occupancy: np.ndarray,
        dt: float,
        max_lag: int | None = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        """
        Computes, up to ``max_lag`` lags:
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

        if max_lag is None:
            n_lags = cls.default_max_lag(n_frames)
        else:
            n_lags = max(1, min(int(max_lag), n_frames))

        time_axis = np.arange(n_lags, dtype=float) * dt
        C = np.zeros(n_lags, dtype=float)

        for lag in range(n_lags):
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
            R = np.ones(n_lags, dtype=float)
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

    @classmethod
    def integrate_correlation(
        cls,
        y: np.ndarray,
        dt: float,
        rule: str = "zero_crossing"
    ) -> Tuple[float, int, str]:
        """Turn a decaying correlation function into a time constant.

        rule
        ----
        ``"zero_crossing"``  integrate up to the first negative value (the
                             point beyond which the curve is noise).
        ``"exponential"``    fit ln y = a - t/tau over the positive part and
                             return tau.  Falls back to ``zero_crossing`` when
                             the curve does not decay.
        ``"full"``           integrate the whole axis (the old behaviour).

        Returns ``(value, n_lags_used, rule_actually_applied)``.
        """
        y = np.asarray(y, dtype=float)

        if rule not in ("zero_crossing", "exponential", "full"):
            raise ValueError(
                f"Unknown integration rule '{rule}'. "
                f"Use 'zero_crossing', 'exponential' or 'full'."
            )

        if y.size == 0:
            return 0.0, 0, rule

        if rule == "full":
            return cls.discrete_integral(y, dt), int(y.size), "full"

        if rule == "exponential":
            positive = y > 0
            if positive.sum() >= 3:
                t = np.arange(y.size, dtype=float) * dt
                slope, _ = np.polyfit(t[positive], np.log(y[positive]), 1)
                if slope < 0:
                    return float(-1.0 / slope), int(positive.sum()), "exponential"
            # did not decay -- fall through to the robust rule
            rule = "zero_crossing"

        negative = np.flatnonzero(y < 0)
        n_used = int(negative[0]) if negative.size else int(y.size)
        return cls.discrete_integral(y[:n_used], dt), n_used, "zero_crossing"

    @staticmethod
    def block_average_uncertainty(values: np.ndarray, n_blocks: int = 5) -> float:
        """Standard error of the mean from contiguous block averages.

        Block averaging is used rather than a naive SEM because successive
        residence events are not independent samples.  Returns NaN when there
        is not enough data to form at least two blocks.
        """
        values = np.asarray(values, dtype=float)
        if values.size < 2:
            return float("nan")

        n_blocks = max(2, min(int(n_blocks), values.size))
        block_means = np.array(
            [block.mean() for block in np.array_split(values, n_blocks) if block.size]
        )
        if block_means.size < 2:
            return float("nan")

        return float(np.std(block_means, ddof=1) / np.sqrt(block_means.size))

    @classmethod
    def scan_tolerance(
        cls,
        occupancy: np.ndarray,
        dt: float,
        tolerances: List[int],
        max_lag: int | None = None
    ) -> List[Dict[str, Any]]:
        """Recompute the residence statistics for a range of t* values.

        The IMM result is known to be extremely sensitive to the tolerance time
        (Laage & Hynes, J. Phys. Chem. B 2008), so the sensitivity is reported
        rather than hidden behind a single chosen value.
        """
        rows: List[Dict[str, Any]] = []

        for tol in tolerances:
            corrected = cls.apply_tolerance_to_occupancy(occupancy, int(tol))
            durations = cls.extract_event_durations(corrected, dt, censor_boundary=True)
            _, survival = cls.compute_origin_averaged_survival(corrected, dt, max_lag)

            rows.append({
                "tolerance_frames": int(tol),
                "n_events": int(durations.size),
                "n_censored": cls.count_boundary_events(corrected),
                "mean_residence_time": float(durations.mean()) if durations.size else 0.0,
                "imm_survival_time": cls.discrete_integral(survival, dt),
            })

        return rows

    # Modes that need atomic masses; every other mode is purely geometric and
    # must not fail just because an element is missing from the mass table.
    MASS_DEPENDENT_MODES = frozenset({"center_of_mass", "group_list_com"})

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
        observed_definition: str,
        cell_lengths: Tuple[float, float, float] | None = None,
        max_lag: int | None = None,
        censor_boundary: bool = True,
        integration_rule: str = "zero_crossing",
        tolerance_scan: List[int] | None = None,
        progress=None
    ) -> MRTResults:
        """Run the full analysis.

        ``progress`` is an optional callable ``progress(fraction, message)``
        invoked as each stage starts, so a caller can show real progress rather
        than an animated placeholder.  It runs on whatever thread calls run(),
        and any exception it raises is swallowed: reporting must never abort an
        analysis.
        """
        def report(fraction: float, message: str):
            if progress is None:
                return
            try:
                progress(fraction, message)
            except Exception:
                pass

        report(0.0, "Reading trajectory...")
        symbols, coords = self.read_xyz(xyz_file)

        needs_masses = (
            reference_mode in self.MASS_DEPENDENT_MODES
            or observed_mode in self.MASS_DEPENDENT_MODES
        )
        masses = self.get_masses(symbols) if needs_masses else np.zeros(len(symbols))

        n_frames, n_atoms, _ = coords.shape

        report(0.45, "Building reference positions...")
        reference_positions = self.build_reference_positions(
            coords=coords,
            masses=masses,
            mode=reference_mode,
            definition=reference_definition
        )

        report(0.55, "Building observed positions...")
        observed_positions, object_labels = self.build_observed_positions(
            coords=coords,
            masses=masses,
            mode=observed_mode,
            definition=observed_definition
        )

        report(0.65, "Computing shell occupancy...")
        _, occupancy = self.build_occupancy(
            reference_positions=reference_positions,
            observed_positions=observed_positions,
            cutoff=cutoff,
            cell_lengths=cell_lengths
        )

        corrected_occupancy = self.apply_tolerance_to_occupancy(
            occupancy=occupancy,
            tolerance_frames=tolerance_frames
        )

        effective_max_lag = (
            self.default_max_lag(n_frames) if max_lag is None
            else max(1, min(int(max_lag), n_frames))
        )

        # --- continuous: two distinct quantities, not one reported twice ----
        report(0.75, "Extracting residence events...")
        event_durations = self.extract_event_durations(
            corrected_occupancy, dt=dt, censor_boundary=censor_boundary
        )
        n_censored = self.count_boundary_events(corrected_occupancy)

        mean_residence = float(np.mean(event_durations)) if event_durations.size else 0.0
        residence_uncertainty = self.block_average_uncertainty(event_durations)

        report(0.82, "Computing survival function...")
        t_surv, S_origin = self.compute_origin_averaged_survival(
            corrected_occupancy, dt=dt, max_lag=effective_max_lag
        )
        imm_survival_time = self.discrete_integral(S_origin, dt)

        # --- intermittent: bounded axis, bounded integral -------------------
        report(0.88, "Computing intermittent correlation...")
        t_corr, C, R, h_mean = self.compute_intermittent_functions(
            corrected_occupancy, dt=dt, max_lag=effective_max_lag
        )
        intermittent_mrt, lags_used, rule_used = self.integrate_correlation(
            R, dt=dt, rule=integration_rule
        )

        if tolerance_scan:
            report(0.95, f"Scanning {len(tolerance_scan)} tolerance values...")
        scan_rows = (
            self.scan_tolerance(occupancy, dt, tolerance_scan, effective_max_lag)
            if tolerance_scan else []
        )

        report(1.0, "Done.")
        return MRTResults(
            xyz_file=xyz_file,
            n_frames=n_frames,
            n_atoms=n_atoms,
            dt=dt,
            dt_unit=dt_unit,
            cutoff=cutoff,
            tolerance_frames=tolerance_frames,
            cell_lengths=cell_lengths,
            reference_mode=reference_mode,
            reference_definition=reference_definition,
            observed_mode=observed_mode,
            observed_definition=observed_definition,
            n_observed_objects=corrected_occupancy.shape[0],
            occupancy_mean=h_mean,
            n_continuous_events=int(event_durations.size),
            n_censored_events=n_censored,
            mean_residence_time=mean_residence,
            mean_residence_uncertainty=residence_uncertainty,
            imm_survival_time=imm_survival_time,
            event_durations=event_durations,
            survival_time_axis=t_surv,
            origin_survival=S_origin,
            intermittent_mrt=intermittent_mrt,
            intermittent_lags_used=lags_used,
            integration_rule=rule_used,
            max_lag=effective_max_lag,
            time_axis=t_corr,
            intermittent_correlation=C,
            intermittent_relaxation=R,
            object_labels=object_labels,
            tolerance_scan=scan_rows
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
            f.write("# Origin-averaged continuous survival (Impey-Madden-McDonald)\n")
            f.write(f"# integral = {results.imm_survival_time:.10f} {results.dt_unit}\n")
            f.write("# time  S(t)\n")
            for t, value in zip(results.survival_time_axis, results.origin_survival):
                f.write(f"{t:20.10f} {value:20.10f}\n")
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
    def _abbreviate(text: str, limit: int = 120) -> str:
        """Shorten a long definition string for display.

        Selecting every solvent hydrogen produces a definition thousands of
        characters long; echoing it in full buries the rest of the summary.
        """
        text = text.strip()
        if len(text) <= limit:
            return text
        n_items = len([p for p in text.replace(";", ",").split(",") if p.strip()])
        return f"{text[:limit].rstrip().rstrip(',')} ... ({n_items} entries in total)"

    @staticmethod
    def _abbreviate_labels(labels: List[str], head: int = 5, tail: int = 2) -> List[str]:
        """List a few object labels rather than all of them."""
        if len(labels) <= head + tail + 1:
            return list(labels)
        hidden = len(labels) - head - tail
        return (
            list(labels[:head])
            + [f"... {hidden} more ..."]
            + list(labels[-tail:])
        )

    @classmethod
    def format_results_text(cls, results: MRTResults) -> str:
        unit = results.dt_unit

        if results.cell_lengths is None:
            pbc_line = "Periodic boundary: none (distances taken directly from the coordinates)"
        else:
            a, b, c = results.cell_lengths
            pbc_line = (
                f"Periodic boundary: orthorhombic minimum image, "
                f"a={a:.4f} b={b:.4f} c={c:.4f} A"
            )

        uncertainty = results.mean_residence_uncertainty
        uncertainty_text = (
            "n/a" if uncertainty != uncertainty else f"+/- {uncertainty:.8f}"
        )

        lines = [
            "Mean Residence Time (MRT) Analysis",
            "=" * 72,
            "",
            "PARAMETERS",
            "-" * 72,
            f"XYZ file: {results.xyz_file}",
            f"Number of frames: {results.n_frames}",
            f"Number of atoms: {results.n_atoms}",
            f"Time between saved frames: {results.dt} {unit}",
            f"Cutoff radius: {results.cutoff} A",
            f"Tolerance frames (t*): {results.tolerance_frames}",
            pbc_line,
            f"Maximum lag used for correlations: {results.max_lag} frames",
            "",
            f"Reference mode: {results.reference_mode}",
            f"Reference definition: {cls._abbreviate(results.reference_definition)}",
            f"Observed mode: {results.observed_mode}",
            f"Observed definition: {cls._abbreviate(results.observed_definition)}",
            f"Number of observed objects: {results.n_observed_objects}",
            "",
            "RESULTS",
            "-" * 72,
            f"Average occupancy <h>: {results.occupancy_mean:.8f}",
            f"Complete residence events: {results.n_continuous_events}",
            f"Censored events (clipped by the start or end of the run): {results.n_censored_events}",
            "",
            f"Mean residence time <T>: {results.mean_residence_time:.8f} {unit}  {uncertainty_text}",
            f"IMM survival integral:   {results.imm_survival_time:.8f} {unit}",
            f"Intermittent time from R(t): {results.intermittent_mrt:.8f} {unit}",
            f"  (integration rule: {results.integration_rule}, "
            f"{results.intermittent_lags_used} of {results.max_lag} lags used)",
            "",
        ]

        if results.tolerance_scan:
            lines.extend([
                "TOLERANCE (t*) SENSITIVITY",
                "-" * 72,
                f"{'t* (frames)':>12} {'events':>8} {'<T>':>16} {'IMM integral':>16}",
            ])
            for row in results.tolerance_scan:
                lines.append(
                    f"{row['tolerance_frames']:>12} {row['n_events']:>8} "
                    f"{row['mean_residence_time']:>16.8f} {row['imm_survival_time']:>16.8f}"
                )
            lines.append("")

        lines.append(f"Observed objects ({len(results.object_labels)}):")
        for label in cls._abbreviate_labels(results.object_labels):
            lines.append(f"  - {label}")

        lines.extend([
            "",
            "NOTES ON THE DEFINITIONS",
            "-" * 72,
            "1. Mean residence time <T> is the average duration of a complete",
            "   residence event. Events still running when the trajectory started",
            "   or ended are censored (excluded), because counting them at their",
            "   truncated length biases <T> downwards, worst for the longest-lived",
            "   species. The censored count is reported above.",
            "2. The IMM survival integral is the integral of the survival function",
            "   averaged over time origins, in the sense of Impey, Madden and",
            "   McDonald, J. Phys. Chem. 87 (1983) 5071. Mathematically this is the",
            "   mean residual time <T^2>/(2<T>); it equals <T> only for exponential",
            "   kinetics and is larger for broad duration distributions. The two",
            "   numbers are reported separately because they answer different",
            "   questions.",
            "3. R(t) is the rescaled intermittent correlation (C(t)-<h>)/(1-<h>).",
            "   It is integrated only over the lag range stated above; the longest",
            "   lags average over very few time origins and are noise.",
            "4. The tolerance time t* has a strong effect on the result (Laage and",
            "   Hynes, J. Phys. Chem. B 112 (2008) 7697). Run a t* scan before",
            "   quoting a residence time from a single value.",
            "5. Without cell lengths no periodic boundary correction is applied, so",
            "   a molecule crossing the box edge registers as leaving and returning.",
        ])
        return "\n".join(lines)


class MeanResidenceTimeUI(DisplayPlots):
    """
    Toga user interface for mean residence time calculations.

    Registered in gqteaWinToga as a Button's on_press, so __init__ receives the
    Button and ignores it.
    """

    # Human-readable label -> the identifier the calculator expects. The raw
    # snake_case names used to be shown to the user verbatim.
    REFERENCE_MODE_LABELS = {
        "Single atom": "single_atom",
        "Geometric centre of a group": "geometric_center",
        "Centre of mass of a group": "center_of_mass",
    }

    OBSERVED_MODE_LABELS = {
        "Single atom": "single_atom",
        "Each atom in a list, tracked separately": "atom_list_individual",
        "Groups, geometric centre": "group_list_geometric",
        "Groups, centre of mass": "group_list_com",
    }

    # mode -> (what to type, a worked example). Drives both the placeholder and
    # the inline hint, which is what makes the free-text boxes usable at all.
    MODE_HELP = {
        "single_atom": (
            "One 1-based atom index.",
            "295",
        ),
        "geometric_center": (
            "Two or more 1-based atom indices, separated by commas or spaces; "
            "ranges like 1-3 are allowed. Their geometric centre is the "
            "reference point in every frame.",
            "1-3",
        ),
        "center_of_mass": (
            "Two or more 1-based atom indices, separated by commas or spaces; "
            "ranges like 1-3 are allowed. Their mass-weighted centre is the "
            "reference point in every frame.",
            "1-3",
        ),
        "atom_list_individual": (
            "A list of 1-based atom indices, with ranges allowed: 3-7,10,15-17 "
            "means 3,4,5,6,7,10,15,16,17. Each atom is followed separately, so "
            "the statistics pool over all of them.",
            "3-7,10,15-17",
        ),
        "group_list_geometric": (
            "Groups separated by ';', atoms within a group by commas, spaces or "
            "ranges. Each group is followed by its geometric centre.",
            "1-3; 4-6",
        ),
        "group_list_com": (
            "Groups separated by ';', atoms within a group by commas, spaces or "
            "ranges. Each group is followed by its centre of mass.",
            "1-3; 4-6",
        ),
    }

    # The rules integrate_correlation() understands, and how they are labelled.
    INTEGRATION_RULES = ("zero_crossing", "exponential", "full")

    INTEGRATION_RULE_LABELS = {
        "Truncate at first zero crossing": "zero_crossing",
        "Fit a single exponential": "exponential",
        "Integrate the whole lag range": "full",
    }

    @staticmethod
    def parse_optional_int(text: str | None) -> int | None:
        """Blank means 'decide automatically'; anything else must be a positive
        integer. Used for the maximum-lag box."""
        if text is None or not str(text).strip():
            return None
        raw = str(text).strip()
        try:
            value = int(raw)
        except ValueError as exc:
            raise ValueError(
                f"Maximum lag must be a whole number of frames, or blank: got '{raw}'."
            ) from exc
        if value < 1:
            raise ValueError("Maximum lag must be at least 1 frame.")
        return value

    @classmethod
    def mode_hint(cls, mode: str) -> str:
        """One-line description of what the definition box expects."""
        return cls.MODE_HELP.get(
            mode, ("Enter 1-based atom indices.", "1")
        )[0]

    @classmethod
    def mode_example(cls, mode: str) -> str:
        return cls.MODE_HELP.get(
            mode, ("Enter 1-based atom indices.", "1")
        )[1]

    def __init__(self, *args, **kwargs):
        # gqteaWinToga registers this class itself as a Button's on_press, so
        # __init__ is handed the Button. Accept and ignore whatever is passed,
        # exactly like every other tool in the project.
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

        # Optional orthorhombic cell: blank means no periodic boundary, the
        # same convention as bond.py and radialDistribution.py.
        self.cell_lengths_input = toga.TextInput(
            placeholder="blank = no PBC, or: 12.4 12.4 12.4",
            style=Pack(width=220, margin=5)
        )

        self.tolerance_scan_input = toga.TextInput(
            placeholder="optional, e.g. 0 1 2 5 10",
            style=Pack(width=200, margin=5)
        )

        periodic_box = toga.Box(
            children=[
                toga.Label("Cell lengths a b c (Å):", style=Pack(width=170, margin=8)),
                self.cell_lengths_input,
                toga.Label("Scan t* values:", style=Pack(width=120, margin=8)),
                self.tolerance_scan_input,
            ],
            style=Pack(direction=ROW, margin_bottom=5)
        )

        # Reference controls. The Selection is created empty and given its
        # items/value further down, once the hint labels its on_change writes
        # to exist.
        self.reference_mode_selection = toga.Selection(style=Pack(width=280, margin=5))

        self.reference_definition_input = toga.TextInput(
            placeholder="Reference definition",
            style=Pack(flex=1, margin=5)
        )

        self.reference_hint_label = toga.Label(
            "", style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 6, 8))
        )

        # Observed controls
        self.observed_definition_input = toga.TextInput(
            placeholder="Observed definition",
            style=Pack(flex=1, margin=5)
        )
        self.observed_hint_label = toga.Label(
            "", style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 6, 8))
        )

        # The hint labels must exist before the Selections get a value, because
        # assigning .value fires on_change -> update_help_texts.
        self.reference_mode_selection.items = list(self.REFERENCE_MODE_LABELS)
        self.reference_mode_selection.on_change = self.update_help_texts
        self.reference_mode_selection.value = "Single atom"

        self.observed_mode_selection = toga.Selection(
            items=list(self.OBSERVED_MODE_LABELS),
            style=Pack(width=280, margin=5),
            on_change=self.update_help_texts,
        )
        self.observed_mode_selection.value = "Each atom in a list, tracked separately"

        definitions_box = toga.Box(
            children=[
                toga.Box(
                    children=[
                        toga.Label("Reference mode:", style=Pack(width=150, margin=8)),
                        self.reference_mode_selection,
                    ],
                    style=Pack(direction=ROW),
                ),
                toga.Box(
                    children=[
                        toga.Label("Reference definition:", style=Pack(width=150, margin=8)),
                        self.reference_definition_input,
                    ],
                    style=Pack(direction=ROW),
                ),
                self.reference_hint_label,
                toga.Divider(style=Pack(margin_top=6, margin_bottom=6)),
                toga.Box(
                    children=[
                        toga.Label("Observed mode:", style=Pack(width=150, margin=8)),
                        self.observed_mode_selection,
                    ],
                    style=Pack(direction=ROW),
                ),
                toga.Box(
                    children=[
                        toga.Label("Observed definition:", style=Pack(width=150, margin=8)),
                        self.observed_definition_input,
                    ],
                    style=Pack(direction=ROW),
                ),
                self.observed_hint_label,
            ],
            style=Pack(direction=COLUMN, margin=10),
        )

        # --- Setup tab ---------------------------------------------------
        self.output_dir_input = toga.TextInput(
            placeholder="blank = next to the trajectory file",
            style=Pack(flex=1, margin=5),
        )
        output_dir_button = toga.Button(
            "Choose...",
            on_press=self.browse_output_dir,
            style=Pack(width=110, margin=5),
        )
        # --- Advanced: the knobs that decide how the correlation is turned
        # into a time constant. Defaults reproduce the standard behaviour.
        self.integration_rule_selection = toga.Selection(
            items=list(self.INTEGRATION_RULE_LABELS),
            style=Pack(width=280, margin=5),
        )
        self.integration_rule_selection.value = "Truncate at first zero crossing"

        self.max_lag_input = toga.TextInput(
            placeholder="blank = a tenth of the trajectory",
            style=Pack(width=220, margin=5),
        )

        self.censor_switch = toga.Switch(
            "Censor events clipped by the start/end of the run",
            style=Pack(margin=5),
        )
        self.censor_switch.value = True

        advanced_box = toga.Box(
            children=[
                toga.Divider(style=Pack(margin_top=8, margin_bottom=6)),
                toga.Label(
                    "Advanced",
                    style=Pack(font_size=11, font_weight="bold", margin=(0, 0, 4, 8)),
                ),
                toga.Box(
                    children=[
                        toga.Label("Integration rule:", style=Pack(width=170, margin=8)),
                        self.integration_rule_selection,
                    ],
                    style=Pack(direction=ROW, margin_bottom=3),
                ),
                toga.Label(
                    "How R(t) becomes a time constant. Use the exponential fit when "
                    "R(t) has not decayed to zero within the lag range.",
                    style=Pack(font_size=9, color=self.HINT_COLOR, margin=(0, 0, 6, 8)),
                ),
                toga.Box(
                    children=[
                        toga.Label("Maximum lag (frames):", style=Pack(width=170, margin=8)),
                        self.max_lag_input,
                    ],
                    style=Pack(direction=ROW, margin_bottom=3),
                ),
                self.censor_switch,
            ],
            style=Pack(direction=COLUMN),
        )

        setup_box = toga.Box(
            children=[
                numeric_box,
                periodic_box,
                toga.Box(
                    children=[
                        toga.Label("Output folder:", style=Pack(width=170, margin=8)),
                        self.output_dir_input,
                        output_dir_button,
                    ],
                    style=Pack(direction=ROW, margin_bottom=5),
                ),
                advanced_box,
            ],
            style=Pack(direction=COLUMN, margin=10),
        )

        # --- Results tab -------------------------------------------------
        self.summary_output = toga.MultilineTextInput(
            readonly=True,
            value="Run an analysis to see the summary here.",
            style=Pack(flex=1, margin=5, font_family="monospace", font_size=9),
        )
        self.plots_button = toga.Button(
            "Show plots",
            on_press=self.show_plots,
            style=Pack(width=140, margin=5),
        )
        results_box = toga.Box(
            children=[
                toga.Box(
                    children=[self.plots_button, self.export_button_placeholder()],
                    style=Pack(direction=ROW),
                ),
                self.summary_output,
            ],
            style=Pack(direction=COLUMN, margin=10, flex=1),
        )

        # --- Help tab ----------------------------------------------------
        help_box = toga.Box(
            children=[
                toga.MultilineTextInput(
                    readonly=True,
                    value=HelpGqteaWin.help_mrt_advanced,
                    style=Pack(flex=1, margin=5, font_size=10),
                )
            ],
            style=Pack(direction=COLUMN, margin=10, flex=1),
        )

        # --- Actions, status, assembly -----------------------------------
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
        self.status_label = toga.Label("Ready.", style=Pack(margin=8, flex=1))
        self.progress_bar = toga.ProgressBar(
            max=100,
            style=Pack(margin_left=5, margin_right=5, margin_bottom=5)
        )
        self.progress_bar.value = 0

        button_box = toga.Box(
            children=[self.run_button, self.export_button, self.clear_button],
            style=Pack(direction=ROW, margin_bottom=5)
        )

        tabs = toga.OptionContainer(
            content=[
                ("Setup", setup_box),
                ("Definitions", definitions_box),
                ("Results", results_box),
                ("Help", help_box),
            ],
            style=Pack(flex=1),
        )

        main_box = toga.Box(
            children=[
                file_box,
                button_box,
                toga.Divider(style=Pack(margin_top=5, margin_bottom=5)),
                tabs,
                toga.Divider(style=Pack(margin_top=5, margin_bottom=2)),
                self.status_label,
                self.progress_bar,
            ],
            style=Pack(direction=COLUMN, margin=10, flex=1)
        )

        self.window.content = main_box
        self.update_help_texts(None)

    @staticmethod
    def export_button_placeholder():
        """Spacer keeping the Results toolbar aligned with the header row."""
        return toga.Label("", style=Pack(flex=1))

    async def _show_error(self, title: str, message: str):
        await self.window.dialog(toga.ErrorDialog(title, message))

    async def _show_info(self, title: str, message: str):
        await self.window.dialog(toga.InfoDialog(title, message))

    HINT_COLOR = "#666666"

    def selected_reference_mode(self) -> str:
        return self.REFERENCE_MODE_LABELS.get(
            self.reference_mode_selection.value, "single_atom"
        )

    def selected_observed_mode(self) -> str:
        return self.OBSERVED_MODE_LABELS.get(
            self.observed_mode_selection.value, "atom_list_individual"
        )

    def update_help_texts(self, widget):
        """Rewrite the placeholder and inline hint whenever a mode changes.

        The definition boxes are free text whose syntax depends on the mode, so
        without this the user has to guess. (This method used to be a stub that
        returned immediately and was never wired to anything.)
        """
        if getattr(self, "reference_hint_label", None) is None:
            return

        ref_mode = self.selected_reference_mode()
        obs_mode = self.selected_observed_mode()

        self.reference_definition_input.placeholder = (
            f"e.g. {self.mode_example(ref_mode)}"
        )
        self.reference_hint_label.text = (
            f"{self.mode_hint(ref_mode)}  Example: {self.mode_example(ref_mode)}"
        )

        self.observed_definition_input.placeholder = (
            f"e.g. {self.mode_example(obs_mode)}"
        )
        self.observed_hint_label.text = (
            f"{self.mode_hint(obs_mode)}  Example: {self.mode_example(obs_mode)}"
        )

    async def browse_output_dir(self, widget):
        try:
            folder = await self.window.dialog(
                toga.SelectFolderDialog(title="Choose the output folder")
            )
            if folder:
                self.output_dir_input.value = str(folder)
        except Exception as exc:
            await self._show_error("Folder Selection Error", str(exc))

    def save_result_plots(self, results: MRTResults):
        """Record the three MRT figures for the interactive viewer.

        Follows the project's 'interactive figures without PNG clutter' recipe:
        save_png=False everywhere, so nothing is left next to the user's data.
        """
        unit = results.dt_unit

        self.save_plots(
            1,
            results.survival_time_axis,
            results.origin_survival,
            f"time ({unit})",
            "S(t)",
            "Continuous survival S(t), origin-averaged",
            save_png=False,
        )

        self.save_multiseries_plot(
            2,
            [
                (results.time_axis, results.intermittent_correlation, "C(t)"),
                (results.time_axis, results.intermittent_relaxation, "R(t)"),
            ],
            f"time ({unit})",
            "correlation",
            "Intermittent correlation C(t) and rescaled R(t)",
            save_png=False,
        )

        if results.event_durations.size:
            counts, edges = np.histogram(results.event_durations, bins="auto")
            centres = 0.5 * (edges[:-1] + edges[1:])
            self.save_plots(
                3,
                centres,
                counts,
                f"event duration ({unit})",
                "count",
                "Residence event duration distribution",
                save_png=False,
            )

    async def show_plots(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Run the MRT calculation first.")
            return
        try:
            self.saved_plot_data = []
            self.saved_plot_files = []
            self.output_dir = self.calculator.resolve_output_dir(
                self.results.xyz_file, self.output_dir_input.value
            )
            self.save_result_plots(self.results)
            self.display_plots()
        except Exception as exc:
            await self._show_error("Plot Error", str(exc))

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
            await self._show_error("File Selection Error", str(exc))

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

        # The Selections show human-readable labels; the calculator wants the
        # identifiers behind them.
        reference_mode = self.selected_reference_mode()
        observed_mode = self.selected_observed_mode()

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

        # Blank cell lengths mean "no periodic boundary"; parse_cell_lengths
        # raises a clear message for anything else that is not three numbers.
        cell_lengths = self.calculator.parse_cell_lengths(
            self.cell_lengths_input.value
        )

        if cell_lengths is not None and cutoff > min(cell_lengths) / 2.0:
            raise ValueError(
                f"The cutoff radius ({cutoff:g} Å) must not exceed half the "
                f"shortest cell edge ({min(cell_lengths) / 2.0:g} Å), otherwise "
                f"the minimum-image convention is invalid."
            )

        scan_text = (self.tolerance_scan_input.value or "").replace(",", " ").strip()
        tolerance_scan = None
        if scan_text:
            try:
                tolerance_scan = [int(v) for v in scan_text.split()]
            except ValueError as exc:
                raise ValueError(
                    "The t* scan must be a list of integers, e.g. '0 1 2 5 10'."
                ) from exc
            if any(v < 0 for v in tolerance_scan):
                raise ValueError("t* scan values cannot be negative.")

        return {
            "xyz_file": xyz_file,
            "dt": dt,
            "dt_unit": dt_unit,
            "cutoff": cutoff,
            "tolerance_frames": tolerance_frames,
            "reference_mode": reference_mode,
            "reference_definition": reference_definition,
            "observed_mode": observed_mode,
            "observed_definition": observed_definition,
            "cell_lengths": cell_lengths,
            "tolerance_scan": tolerance_scan,
            "max_lag": self.parse_optional_int(self.max_lag_input.value),
            "censor_boundary": bool(self.censor_switch.value),
            "integration_rule": self.INTEGRATION_RULE_LABELS.get(
                self.integration_rule_selection.value, "zero_crossing"
            ),
            "output_dir": self.calculator.resolve_output_dir(
                xyz_file, self.output_dir_input.value
            ),
        }

    def _apply_progress(self, fraction: float, message: str):
        """Move the bar and the status line. Runs on the Toga event loop."""
        self.progress_bar.value = max(0.0, min(100.0, fraction * 100.0))
        self.status_label.text = message

    async def run_mrt_calculation(self, widget):
        try:
            params = self._collect_inputs()
            # output_dir is the UI's business, not the calculator's.
            self.output_dir = params.pop("output_dir")
        except Exception as exc:
            await self._show_error("Input Error", str(exc))
            return

        self.status_label.text = "Reading trajectory and computing MRT..."
        self.progress_bar.value = 0
        self.run_button.enabled = False
        self.export_button.enabled = False
        self.clear_button.enabled = False
        await asyncio.sleep(0)

        try:
            # run() reports each stage as it starts. It executes on a worker
            # thread, so hop back to the event loop before touching widgets.
            loop = asyncio.get_running_loop()

            def on_progress(fraction: float, message: str):
                loop.call_soon_threadsafe(self._apply_progress, fraction, message)

            worker = asyncio.create_task(
                asyncio.to_thread(self.calculator.run, progress=on_progress, **params)
            )
            while not worker.done():
                await asyncio.sleep(0.05)

            results = await worker
            self.results = results
            self.summary_output.value = self.calculator.format_results_text(results)
            self.status_label.text = (
                "Done. "
                f"Objects: {results.n_observed_objects} | "
                f"<T>: {results.mean_residence_time:.4f} {results.dt_unit} | "
                f"IMM: {results.imm_survival_time:.4f} {results.dt_unit} | "
                f"events: {results.n_continuous_events} "
                f"({results.n_censored_events} censored)"
            )
            self.progress_bar.value = 100
        except Exception as exc:
            self.status_label.text = "Calculation failed."
            self.progress_bar.value = 0
            await self._show_error("MRT Calculation Error", str(exc))
        finally:
            # The bar is determinate now (start()/stop() drove the old
            # indeterminate animation), so there is nothing to stop.
            self.run_button.enabled = True
            self.export_button.enabled = True
            self.clear_button.enabled = True

    async def export_results(self, widget):
        if self.results is None:
            await self._show_info("No Results", "Please run the MRT calculation first.")
            return

        try:
            files = self.calculator.export_curves(
                self.results,
                output_dir=self.calculator.resolve_output_dir(
                    self.results.xyz_file, self.output_dir_input.value
                ),
            )
            msg = (
                "Results exported successfully:\n\n"
                f"Summary: {files['summary']}\n"
                f"Continuous survival: {files['continuous_survival']}\n"
                f"Intermittent correlation: {files['intermittent_correlation']}\n"
                f"Event durations: {files['event_durations']}"
            )
            self.status_label.text = "Results exported."
            await self._show_info("Export Successful", msg)
        except Exception as exc:
            await self._show_error("Export Error", str(exc))

    def clear_output(self, widget):
        # Discard the results too: without this, Export after Clear silently
        # re-wrote the previous run's curves under the current file's name.
        self.results = None
        self.summary_output.value = "Run an analysis to see the summary here."
        self.status_label.text = "Output cleared."
        self.progress_bar.value = 0

