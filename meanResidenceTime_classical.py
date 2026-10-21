import asyncio
import itertools
import math
import os
import time
from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple

import toga
import numpy as np
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from help import HelpGqteaWin
from framesCounter import FramesCounter
from displayPlots import DisplayPlots

ATU_FS = 0.02418884326505  # fs per atomic unit of time

TIME_STEP_UNIT_AU = "a.u."
TIME_STEP_UNIT_FS = "fs"
TIME_STEP_UNITS = (TIME_STEP_UNIT_AU, TIME_STEP_UNIT_FS)   # a.u. first: the CPMD convention

TOLERANCE_UNIT_FRAMES = "frames"
TOLERANCE_UNIT_FS = "fs"
TOLERANCE_UNITS = (TOLERANCE_UNIT_FRAMES, TOLERANCE_UNIT_FS)
DEFAULT_TOLERANCE = 0

# t / dt is computed in floating point, so a tolerance that is an exact
# multiple of the frame interval can come out a hair below the whole number
# (3 * dt / dt = 2.9999999999999996); this much is added before rounding down.
FRAME_ROUNDING_SLACK = 1e-9


def frame_interval_fs(time_step: float, time_step_unit: str, sampling: int) -> float:
    """Time between two sampled frames, in fs: simulation time step (fs) x sampling interval."""
    if time_step_unit == TIME_STEP_UNIT_AU:
        step_fs = time_step * ATU_FS
    elif time_step_unit == TIME_STEP_UNIT_FS:
        step_fs = time_step
    else:
        raise ValueError(f"Unknown time step unit '{time_step_unit}'; use {' or '.join(TIME_STEP_UNITS)}.")
    return step_fs * sampling


def tolerance_in_frames(value: float, unit: str, interval_fs: float) -> int:
    """The tolerance as a whole number of sampled frames.

    The tolerance is the largest number of consecutive sampled frames a
    molecule may spend outside the shell without that excursion counting as
    an exit. In frames the value must already be a whole number. In fs it is
    converted with ``frames = time / interval_fs`` and **rounded down**: an
    excursion of n frames lasts n x interval_fs, and it is tolerated when that
    fits within the time given, i.e. n <= time / interval_fs. So n frames and
    n x interval_fs fs always give the same result, and a time between two
    multiples of the interval acts as the lower one (e.g. 2.5 frames -> 2).
    """
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"The tolerance must be zero or a positive number, not {value:g}.")
    if unit == TOLERANCE_UNIT_FRAMES:
        if not value.is_integer():
            raise ValueError(
                f"A tolerance in frames must be a whole number, not {value:g}. "
                "Choose fs as the unit to give it as a time."
            )
        return int(value)
    if unit == TOLERANCE_UNIT_FS:
        if not (math.isfinite(interval_fs) and interval_fs > 0):
            raise ValueError("The time between sampled frames must be greater than zero.")
        return int(math.floor(value / interval_fs + FRAME_ROUNDING_SLACK))
    raise ValueError(f"Unknown tolerance unit '{unit}'; use {' or '.join(TOLERANCE_UNITS)}.")


# The trajectory is read in blocks of this many frames: it bounds the memory
# held at once and how long one block's bulk string operations keep the GIL
# from the window's thread while the analysis runs in a worker.
READ_CHUNK_FRAMES = 100
COUNT_BLOCK_BYTES = 1 << 24    # bytes read at a time when counting frames
STATUS_INTERVAL = 0.25         # s between two status-line updates during a job

# progress(message, force): the status line; force bypasses the throttling.
ProgressCallback = Callable[[str, bool], None]


def _notify(progress: Optional[ProgressCallback], message: str, force: bool = False) -> None:
    """Report to the status line; a failing callback never stops the analysis."""
    if progress is None:
        return
    try:
        progress(message, force)
    except Exception:
        pass


def count_frames(path: str, progress: Optional[Callable[[float], None]] = None) -> Tuple[int, int]:
    """(atoms per frame, number of frames) of an .xyz trajectory, from its line count.

    Newline bytes are counted instead of parsing the lines (0.4 s instead of
    2.5 s for a 351 MB file); blank lines at the end of the file are ignored.
    The lines themselves are checked when the trajectory is read for the
    analysis (:func:`scan_trajectory`). ``progress(fraction)`` follows the
    bytes read.
    """
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        fields = f.readline().split()
        if not fields:
            raise ValueError("The trajectory file is empty.")
        try:
            num_atoms = int(fields[0])
        except ValueError:
            raise ValueError("The first line of the trajectory must hold the number of atoms.") from None
        if num_atoms <= 0:
            raise ValueError("The number of atoms on the first line must be greater than zero.")
        f.seek(0)
        newlines = trailing = done = 0
        while block := f.read(COUNT_BLOCK_BYTES):
            newlines += block.count(b"\n")
            content_end = len(block.rstrip())
            # Newlines after the last non-blank character of the file so far
            trailing = block.count(b"\n", content_end) if content_end else trailing + block.count(b"\n")
            done += len(block)
            if progress is not None:
                progress(done / size)
    lines = newlines - trailing + 1
    stride = num_atoms + 2
    if lines % stride:
        raise ValueError(
            f"The trajectory has {lines} lines, which is not a whole number of {stride}-line frames "
            f"({num_atoms} atoms plus 2 header lines): the last frame is incomplete or the number "
            "of atoms changes."
        )
    return num_atoms, lines // stride


@dataclass
class ShellScan:
    """What :func:`scan_trajectory` found."""
    shell_atom_list: List[int]   # 1-based labels of the atoms that entered the shell, in order of first entry
    occupancy: np.ndarray        # int8 (atoms, frames): 1 = inside the shell, before any tolerance
    centre_symbol: str


def scan_trajectory(path: str, centre_label: int, symbol: str, exclude_labels: Sequence[int],
                    outer_radius: float, inner_radius: float = 0.0, n_frames: Optional[int] = None,
                    progress: Optional[Callable[[int], None]] = None,
                    chunk_frames: int = READ_CHUNK_FRAMES) -> ShellScan:
    """Which atoms of ``symbol`` are inside the shell round atom ``centre_label``, frame by frame.

    One pass over the trajectory. Frame 1 is checked line by line and fixes
    the candidates (every atom of ``symbol`` except the centre and the
    excluded labels). In later frames every atom-count line is checked, but
    only the centre's and the candidates' lines are split and converted, in
    bulk per block of frames; the other lines are never parsed, so a
    malformed coordinate on an atom the analysis does not use goes
    unnoticed (the bond tool makes the same trade-off). A parsed line must
    have four fields and the element it has in frame 1.

    The distance is ``np.linalg.norm`` of the centre -> atom vector, the same
    float64 value the per-pair loop this replaced computed, and an atom is
    inside when ``inner_radius <= d <= outer_radius``. The shell list is in
    order of first entry (atom order within a frame), which is the column
    order of mrt.dat.

    ``n_frames`` is the count from :func:`count_frames` (counted here when
    None): a file that ends early or holds more than that raises
    ``ValueError``. ``progress(frames_read)`` is called after each block.
    """
    if n_frames is None:
        n_frames = count_frames(path)[1]
    exclude = set(exclude_labels)
    with open(path, "r") as f:
        fields = f.readline().split()
        if not fields:
            raise ValueError("The trajectory file is empty.")
        num_atoms = int(fields[0])
        if not 1 <= centre_label <= num_atoms:
            raise ValueError(f"The shell center atom label must be between 1 and {num_atoms}.")
        f.seek(0)
        stride = num_atoms + 2
        centre = centre_label - 1

        symbols = None
        frames_read = 0
        inside_blocks = []
        while frames_read < n_frames:
            n_block = min(chunk_frames, n_frames - frames_read)
            lines = list(itertools.islice(f, n_block * stride))
            if len(lines) < n_block * stride:
                raise ValueError(f"Invalid atom record found in frame {frames_read + len(lines) // stride + 1}.")
            for j in range(n_block):
                header = lines[j * stride].split()
                frame = frames_read + j + 1
                if not header or header[0] != fields[0]:
                    found = f"'{header[0]}'" if header else "a blank line"
                    raise ValueError(f"Frame {frame} should start with the atom count {num_atoms}, "
                                     f"but {found} was found.")

            if symbols is None:
                symbols = []
                for line in lines[2:stride]:
                    if len(line.split()) != 4:
                        raise ValueError("Invalid atom record found in frame 1.")
                    symbols.append(line.split()[0])
                candidates = [i for i, s in enumerate(symbols)
                              if i != centre and s == symbol and (i + 1) not in exclude]
                needed = [centre] + candidates
                expected = np.array([symbols[i] for i in needed], dtype=object)

            picked = []
            for i in needed:
                picked.extend(lines[2 + i::stride])
            tokens = " ".join(picked).split()
            if len(tokens) != 4 * len(picked):
                for k, line in enumerate(picked):
                    if len(line.split()) != 4:
                        raise ValueError(f"Invalid atom record found in frame {frames_read + k % n_block + 1}.")
            table = np.array(tokens, dtype=object).reshape(len(needed), n_block, 4)
            moved = np.argwhere(table[:, :, 0] != expected[:, None])
            if moved.size:
                a, j = moved[0]
                raise ValueError(
                    f"In frame {frames_read + j + 1} atom {needed[a] + 1} is {table[a, j, 0]}, but it is "
                    f"{expected[a]} in frame 1: the atom order must be the same in every frame."
                )
            xyz = table[:, :, 1:].astype(np.float64)    # float() of each field, as the old loop did
            distance = np.linalg.norm(xyz[1:] - xyz[0], axis=-1)
            inside_blocks.append((distance >= inner_radius) & (distance <= outer_radius))
            frames_read += n_block
            if progress is not None:
                progress(frames_read)

        if any(line.strip() for line in f):
            raise ValueError("The trajectory holds more frames than were counted when it was loaded; "
                             "load it again with Browse.")

    if symbols is None:   # no frame at all
        return ShellScan([], np.zeros((0, 0), np.int8), "")
    inside = np.concatenate(inside_blocks, axis=1)
    entered = np.flatnonzero(inside.any(axis=1))
    first_entry = inside.argmax(axis=1)
    entered = entered[np.lexsort((entered, first_entry[entered]))]
    return ShellScan(
        shell_atom_list=[candidates[k] + 1 for k in entered],
        occupancy=inside[entered].astype(np.int8),
        centre_symbol=symbols[centre],
    )


def _excursions(occupancy: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Every stay outside the shell that begins with an inside -> outside step.

    Returns (row, start, end) arrays: frames start .. end-1 of that row are
    outside, frame start-1 is inside, and end is the frame the atom is back
    inside, or n_frames if it never comes back. Leading outside frames are
    not an excursion.
    """
    n_frames = occupancy.shape[1]
    step = np.diff((occupancy == 1).astype(np.int8), axis=1)
    out_rows, out_at = np.nonzero(step == -1)     # inside at t, outside at t+1
    in_rows, in_at = np.nonzero(step == 1)        # outside at t, inside at t+1
    starts = out_at + 1
    # np.nonzero is row-major, so these keys are sorted: the first return at
    # or after each start in the same row is that excursion's end.
    return_keys = in_rows * (n_frames + 1) + in_at + 1
    k = np.searchsorted(return_keys, out_rows * (n_frames + 1) + starts)
    found = k < return_keys.size
    found[found] = in_rows[k[found]] == out_rows[found]
    ends = np.full(starts.shape, n_frames)
    ends[found] = in_at[k[found]] + 1
    return out_rows, starts, ends


def exits_per_atom(occupancy: np.ndarray, tolerance_frames: int) -> np.ndarray:
    """Number of exits of every row of a tolerance-corrected occupancy (see MRTAnalyzer.count_exits)."""
    occupancy = np.asarray(occupancy)
    if occupancy.shape[1] < 2:
        return np.zeros(occupancy.shape[0], dtype=np.int64)
    rows, starts, ends = _excursions(occupancy)
    counted = (ends < occupancy.shape[1]) | (ends - starts > tolerance_frames)
    return np.bincount(rows[counted], minlength=occupancy.shape[0])


class MRTAnalyzer(FramesCounter,DisplayPlots):

    @staticmethod
    def apply_tolerance_to_occupancy(occupancy: np.ndarray, tolerance_frames: int) -> np.ndarray:
        """Fill in every gap of at most ``tolerance_frames`` frames between two in-shell frames.

        Gaps at the start or the end of a row are left as they are.
        """
        corrected = occupancy.copy()
        if tolerance_frames <= 0 or corrected.shape[1] < 2:
            return corrected
        n_frames = corrected.shape[1]
        rows, starts, ends = _excursions(corrected)
        bridged = (ends < n_frames) & (ends - starts <= tolerance_frames)
        # +1 where a bridged gap starts, -1 where it ends: the running sum marks the gap
        marks = np.zeros((corrected.shape[0], n_frames + 1), dtype=np.int32)
        marks[rows[bridged], starts[bridged]] = 1
        marks[rows[bridged], ends[bridged]] = -1
        corrected[np.cumsum(marks, axis=1)[:, :n_frames] > 0] = 1
        return corrected

    @staticmethod
    def count_exits(occupancy_row, tolerance_frames: int) -> int:
        """Exits from the shell in one atom's tolerance-corrected occupancy (1 = inside).

        An exit is an inside -> outside step whose excursion lasts more than
        ``tolerance_frames`` sampled frames. Excursions between two in-shell
        frames that were no longer have already been filled in by
        :meth:`apply_tolerance_to_occupancy`, so every such step left is an
        exit. An excursion still under way when the trajectory ends counts
        only if it already lasts more than the tolerance. With a tolerance of
        0 every inside -> outside step is an exit.
        """
        row = np.asarray([int(value) for value in occupancy_row], dtype=np.int64).reshape(1, -1)
        return int(exits_per_atom(row, tolerance_frames)[0])

    async def warning_function(self, msg):
        await self.main_window.dialog(toga.InfoDialog("Error", f"{msg}"))

    async def read_params(self, widget):

        async def read_input(text_input, field_name, expected_type):
            value = text_input.value.strip()
            if not value:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Please input a valid value for {field_name}.")
                )
                return None
            try:
                if expected_type == list:
                    # Assuming atom labels are space-separated integers
                    labels = [int(label) for label in value.split()]
                    return labels
                else:
                    return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        if not hasattr(self, "trajec") or not self.trajec:
            await self.warning_function(
                "Please select a TRAJEC.xyz file before starting the MRT calculation."
            )
            return False

        if not hasattr(self, "total_frame_number") or self.total_frame_number <= 0:
            await self.warning_function(
                "Please load the trajectory with Browse before starting the MRT calculation."
            )
            return False

        # Read and validate inputs
        self.inner_radius = 0.0
        self.outer_radius = await read_input(self.textInput_outer_radius, "radius for the outer shell", float)
        if self.outer_radius is None:
            return False

        self.center_label = await read_input(self.textInput_center_label, "atom label for the center of the shell", int)
        if self.center_label is None:
            return False

        self.atom_symbol = await read_input(self.textInput_atom_symbol, "atom symbol for mean residence time", str)
        if self.atom_symbol is None:
            return False
        self.atom_symbol = self.atom_symbol.strip()

        self.time_step = await read_input(self.textInput_time_step, "simulation time step", float)
        if self.time_step is None:
            return False
        self.time_step_unit = self.selection_time_step_unit.value or TIME_STEP_UNIT_AU

        self.sampling = await read_input(self.textInput_sampling, "sampling interval between two colected frames", int)
        if self.sampling is None:
            return False

        # Blank means the default tolerance, 0: every exit is counted.
        self.tolerance_unit = self.selection_tolerance_unit.value or TOLERANCE_UNIT_FRAMES
        tolerance_text = self.textInput_tolerance.value.strip()
        try:
            self.tolerance_value = float(tolerance_text) if tolerance_text else float(DEFAULT_TOLERANCE)
        except ValueError:
            await self.warning_function(
                f"Invalid format for the tolerance: '{tolerance_text}' is not a number."
            )
            return False

        self.exclude_atom_list = await read_input(self.textInput_atom_labels, "atom labels list to exclude", list)
        if self.exclude_atom_list is None:
            return False

        # Adjust for zero-based indexing (assuming user input is one-based)
        if self.exclude_atom_list != [0]:
            self.exclude_atom_list = [label for label in self.exclude_atom_list]
        else:
            self.exclude_atom_list = []

        if self.outer_radius < 0:
            await self.warning_function("The cutoff radius must be non-negative.")
            return False

        if self.center_label < 1 or self.center_label > self.num_atoms:
            await self.warning_function(
                f"The shell center atom label must be between 1 and {self.num_atoms}."
            )
            return False

        if self.time_step <= 0:
            await self.warning_function("The simulation time step must be greater than zero.")
            return False

        if self.sampling <= 0:
            await self.warning_function("The sampling interval must be greater than zero.")
            return False

        self.frame_interval_fs = frame_interval_fs(self.time_step, self.time_step_unit, self.sampling)
        try:
            self.tolerance_frames = tolerance_in_frames(
                self.tolerance_value, self.tolerance_unit, self.frame_interval_fs
            )
        except ValueError as e:
            await self.warning_function(str(e))
            return False

        invalid_excluded_labels = [
            label for label in self.exclude_atom_list if label < 1 or label > self.num_atoms
        ]
        if invalid_excluded_labels:
            await self.warning_function(
                f"Excluded atom labels must be between 1 and {self.num_atoms}."
            )
            return False

        self.multi_line_text.value = self.parameters_text() + "\nRunning the MRT calculation...\n"
        return True

    def tolerance_text(self) -> str:
        """The tolerance as entered and as applied, e.g. '25 fs -> 2 frames (24.19 fs)'."""
        applied = (f"{self.tolerance_frames} frame{'s' if self.tolerance_frames != 1 else ''} "
                   f"({self.tolerance_frames * self.frame_interval_fs:.6g} fs)")
        if self.tolerance_unit == TOLERANCE_UNIT_FS:
            return f"{self.tolerance_value:g} fs -> {applied}"
        return applied

    def parameters_text(self) -> str:
        """The parameters of the run, as shown in the text box."""
        exits = ("every exit from the shell is counted" if self.tolerance_frames == 0 else
                 f"excursions of up to {self.tolerance_frames} sampled frames are not counted as exits")
        return (
            f"Cutoff radius --> {self.outer_radius} Å\n"
            f"Atom label for the center of the shell --> {self.center_label}\n"
            f"Atom symbol for mean residence time --> {self.atom_symbol}\n"
            f"Simulation time step --> {self.time_step:g} {self.time_step_unit}\n"
            f"Sampling interval between two collected frames --> {self.sampling}\n"
            f"Time between sampled frames --> {self.frame_interval_fs:.6g} fs\n"
            f"Tolerance --> {self.tolerance_text()}: {exits}\n"
            f"Atom labels list to exclude --> {self.exclude_atom_list}\n"
        )

    def _validate_analysis_state(self):
        required_attributes = [
            "trajec",
            "output_dir",
            "inner_radius",
            "outer_radius",
            "center_label",
            "atom_symbol",
            "time_step",
            "time_step_unit",
            "sampling",
            "tolerance_frames",
            "exclude_atom_list",
            "total_frame_number",
        ]
        missing_attributes = [
            attr for attr in required_attributes if not hasattr(self, attr)
        ]
        if missing_attributes:
            missing_attrs = ", ".join(missing_attributes)
            raise ValueError(
                f"Missing required analysis data: {missing_attrs}. Please reload the inputs and trajectory."
            )

    def run_analysis(self, progress: Optional[ProgressCallback] = None):
        """The whole calculation after read_params: mrt.dat, mrt_summary.dat, mrt_total.dat and the figure.

        Reads the trajectory once and touches no widget, so the window runs
        it in a worker thread; ``progress(message, force)`` feeds the status line.
        """
        self.calc_mrt(progress)
        _notify(progress, "Writing mrt_summary.dat and mrt_total.dat...", True)
        self.write_summary()
        self.write_mrt_total()
        _notify(progress, "Preparing the figure...", True)
        self.record_plot()

    def calc_mrt(self, progress: Optional[ProgressCallback] = None):
        """Read the trajectory, apply the tolerance and write mrt.dat.

        Sets ``shell_atom_list`` (1-based labels of the atoms of the chosen
        symbol that entered the shell at least once, in order of first
        entry), ``occupancy`` (int8, one row per such atom, one column per
        frame, 1 = inside, tolerance-corrected), ``center_atom_symbol`` and
        ``sim_time`` (ps between sampled frames).
        """
        self._validate_analysis_state()
        # Time between two sampled frames, in ps (time step in a.u. or fs x sampling interval)
        self.sim_time = frame_interval_fs(self.time_step, self.time_step_unit, self.sampling) / 1000.0
        total = self.total_frame_number
        _notify(progress, "Reading the trajectory...", True)
        scan = scan_trajectory(
            self.trajec, self.center_label, self.atom_symbol, self.exclude_atom_list,
            self.outer_radius, self.inner_radius, n_frames=total,
            progress=lambda frames: _notify(progress, f"Reading the trajectory... frame {frames} of {total}"),
        )
        self.shell_atom_list = scan.shell_atom_list
        self.center_atom_symbol = scan.centre_symbol
        _notify(progress, "Applying the tolerance...", True)
        self.occupancy = self.apply_tolerance_to_occupancy(scan.occupancy, self.tolerance_frames)
        _notify(progress, "Writing mrt.dat...", True)
        self.write_mrt_dat()

    def write_mrt_dat(self):
        header = [f"{self.atom_symbol}{s}" for s in self.shell_atom_list]
        lines = [f"Frame  Time(ps)  {'    '.join(header)}\n"]
        for frame_index, row in enumerate(self.occupancy.T.tolist(), start=1):
            Time = self.sim_time * frame_index
            lines.append(f"{frame_index:<6d} {Time:>8.4f}   {'      '.join(str(s) for s in row)}\n")
        with open(f"{self.output_dir}/mrt.dat", "w") as mrt:
            mrt.writelines(lines)

    def frame_times(self) -> List[float]:
        """The time of every frame, in ps, as written in mrt.dat (4 decimals)."""
        return [float(f"{self.sim_time * i:.4f}") for i in range(1, self.occupancy.shape[1] + 1)]

    def write_mrt_total(self):
        """mrt_total.dat: time and the number of atoms inside the shell in every frame."""
        with open(f"{self.output_dir}/mrt_total.dat", "w") as f:
            f.write("Time(ps)  Total coordination number\n")
            if self.shell_atom_list:
                totals = self.occupancy.sum(axis=0).tolist()
                for time_ps, total in zip(self.frame_times(), totals):
                    f.write(f"{time_ps:>8.4f}  {total:>10d}\n")

    def write_summary(self):
        """mrt_summary.dat, and ``mrt_results``: (label, total rt, mrt, exits) per atom, times in ps."""
        symbol = self.atom_symbol # Atom symbol for mrt
        rt = self.occupancy.sum(axis=1).tolist()   # frames inside the shell, per atom
        # self.occupancy is already tolerance-corrected
        exits = exits_per_atom(self.occupancy, self.tolerance_frames).tolist()

        with open(f"{self.output_dir}/mrt_summary.dat", "w") as f:
            total_time_sim = self.sim_time * self.total_frame_number  # Total time simulation
            f.write(f"Total simulation time.........{total_time_sim:>10.4f}\n")
            f.write(f"Time between frames (fs)......{self.sim_time * 1000.0:>10.4f}\n")
            f.write(f"Tolerance (frames)............{self.tolerance_frames:>10d}\n")
            f.write(f"Tolerance (fs)................{self.tolerance_frames * self.sim_time * 1000.0:>10.4f}\n\n")
            self.mrt_results = []
            for n in range(len(rt)):
                tot_rt = rt[n] * self.sim_time
                if exits[n] != 0:
                    tmp_rt = rt[n] * self.sim_time / exits[n]
                else:
                    tmp_rt = rt[n] * self.sim_time  # never seen to leave: a lower bound
                self.mrt_results.append((self.shell_atom_list[n], tot_rt, tmp_rt, exits[n]))

                f.write(
                    f"{symbol}{self.shell_atom_list[n]} total rt.................{tot_rt:>10.4f}\n"
                )
                f.write(
                    f"{symbol}{self.shell_atom_list[n]} mrt......................{tmp_rt:>10.4f}\n"
                )
                f.write(
                    f"{symbol}{self.shell_atom_list[n]} exits....................{exits[n]:>10d}\n\n"
                )

    def mrt_summary(self):
        """write_summary, then the results in the text box (the synchronous path)."""
        self.write_summary()
        self.multi_line_text.value = self.results_text(self.center_atom_symbol)

    def results_text(self, center_atom_symbol: str) -> str:
        """The closing message: parameters, per-atom results (times in ps) and the output files."""
        lines = [
            f"MEAN RESIDENCE TIME ANALYSIS FOR {center_atom_symbol}{self.center_label} SITE FINISHED SUCCESSFULLY!",
            "",
            self.parameters_text(),
            "Results (times in ps; mrt = total rt / exits, or total rt when no exit was seen):",
        ]
        if not self.mrt_results:
            lines.append(f"  No {self.atom_symbol} atom entered the shell.")
        for label, tot_rt, mrt, exits in self.mrt_results:
            lines.append(f"  {self.atom_symbol}{label:<6} total rt {tot_rt:>10.4f}   mrt {mrt:>10.4f}   exits {exits:>6d}")
        lines += ["", f"Files written to {self.output_dir}: mrt.dat, mrt_total.dat, mrt_summary.dat"]
        return "\n".join(lines) + "\n"


    def record_plot(self):
        """Record the coordination-number figure (time, atoms inside) for display_plots."""
        x_axis = self.frame_times()
        y_axis = self.occupancy.sum(axis=0).tolist()

        # Plot labels and title
        plot_xlabel = "Time (ps)"
        plot_ylabel = "Coordination number"
        plot_title = "Mean Residence Time"

        # Save the plot with provided data
        self.save_plots(1, x_axis, y_axis, plot_xlabel, plot_ylabel, plot_title)


class ClassicalMeanResidenceTimeUI(MRTAnalyzer):
    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self,widget):
        # Create the main window
        self.main_window = toga.Window(
            title="Mean Residence Time (MRT)",
            size=(700, 600),
        )

        # Define common styles
        heading_style = Pack(
            font_size=18, font_weight="bold", text_align=LEFT, margin=(0, 0, 10, 0)
        )
        label_style = Pack(
            margin=(0, 0, 5, 5), text_align=LEFT, width=200
        )
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=100)
        box_style = Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0))

        # Main container
        main_box = toga.Box(style=Pack(direction=COLUMN, margin=20))

        # Title
        main_box.add(toga.Label("Mean Residence Time (MRT)", style=heading_style))

        # Input fields with labels
        input_fields = [
            (
                "Cutoff radius (\u00c5):",
                "Enter the cutoff radius",
                "textInput_outer_radius",
            ),
            (
                "Atom label for the shell center:",
                "Enter the atom label for the center of the shell",
                "textInput_center_label",
            ),
            (
                "Atom symbol for MRT:",
                "Enter the atom symbol for mean residence time",
                "textInput_atom_symbol",
            ),
            (
                "Atom label list to exclude:",
                "Enter the list of atom labels to exclude from the MRT analysis",
                "textInput_atom_labels",
            ),
        ]

        for label_text, placeholder, attr_name in input_fields:
            box = toga.Box(style=box_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            box.add(label)
            box.add(text_input)
            main_box.add(box)

        time_box = toga.Box(style=box_style)
        time_label = toga.Label("Simulation time step:", style=label_style)
        self.textInput_time_step = toga.TextInput(
            placeholder="Time step",
            style=Pack(width=90, margin=(5, 5)),
        )
        self.selection_time_step_unit = toga.Selection(
            items=list(TIME_STEP_UNITS), style=Pack(width=60, margin=(5, 5))
        )
        self.selection_time_step_unit.value = TIME_STEP_UNIT_AU
        sampling_label = toga.Label(
            "Sampling interval (steps):",
            style=Pack(margin=(0, 0, 5, 10), text_align=LEFT, width=150),
        )
        self.textInput_sampling = toga.TextInput(
            placeholder="MD steps per frame",
            style=Pack(width=110, margin=(5, 5)),
        )
        time_box.add(time_label)
        time_box.add(self.textInput_time_step)
        time_box.add(self.selection_time_step_unit)
        time_box.add(sampling_label)
        time_box.add(self.textInput_sampling)
        main_box.add(time_box)

        tolerance_box = toga.Box(style=box_style)
        self.tolerance_label = toga.Label("", style=label_style)
        self.textInput_tolerance = toga.TextInput(value=str(DEFAULT_TOLERANCE), style=input_style)
        self.selection_tolerance_unit = toga.Selection(
            items=list(TOLERANCE_UNITS), style=Pack(width=80, margin=(5, 5))
        )
        self.selection_tolerance_unit.value = TOLERANCE_UNIT_FRAMES
        # Wired after the value is set, then run once to write the label.
        self.selection_tolerance_unit.on_change = self.on_tolerance_unit_change
        self.on_tolerance_unit_change(self.selection_tolerance_unit)
        tolerance_box.add(self.tolerance_label)
        tolerance_box.add(self.textInput_tolerance)
        tolerance_box.add(self.selection_tolerance_unit)
        main_box.add(tolerance_box)

        # File selection button
        file_box = toga.Box(style=box_style)
        file_label = toga.Label("Select input file:", style=label_style)
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select TRAJEC.xyz file",
            style=input_style,
        )
        self.browse_button = toga.Button(
            "Browse", on_press=self.frames_counter, style=button_style
        )
        file_box.add(file_label)
        file_box.add(self.textInput_file)
        file_box.add(self.browse_button)
        main_box.add(file_box)

        # Activity: a spinner while a job runs, and what it is doing. Starting
        # or stopping the spinner only shows/hides it, so the window is not
        # re-laid out.
        activity_box = toga.Box(style=box_style)
        self.activity_indicator = toga.ActivityIndicator(style=Pack(margin=(0, 8, 0, 5)))
        self.status_label = toga.Label("Choose a TRAJEC.xyz file with Browse.", style=Pack(flex=1))
        activity_box.add(self.activity_indicator)
        activity_box.add(self.status_label)
        main_box.add(activity_box)

        # Multi-line text input for output
        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(5, 0), font_size=11)
        )
        self.multi_line_text.value = HelpGqteaWin.help_mrt
        main_box.add(self.multi_line_text)

        # Buttons at the bottom
        button_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin_top=5)
        )
        self.btn_execute = toga.Button(
            "MRT calculation", style=button_style, on_press=self.workflow
        )
        self.btn_close = toga.Button(
            "Close", style=button_style, on_press=self.closeTopLevel
        )
        button_box.add(self.btn_execute)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        # Set the content of the main window
        self.main_window.content = main_box
        self.main_window.on_close = self._on_close
        self.main_window.show()

    TOLERANCE_LABELS = {
        TOLERANCE_UNIT_FRAMES: ("Tolerance (frames):",
                                "Frames outside the shell not counted as an exit (0 = count every exit)"),
        TOLERANCE_UNIT_FS: ("Tolerance (fs):",
                            "Time outside the shell not counted as an exit (0 = count every exit)"),
    }

    def on_tolerance_unit_change(self, widget, **kwargs):
        """Name the tolerance field after the chosen unit."""
        label, placeholder = self.TOLERANCE_LABELS.get(
            self.selection_tolerance_unit.value, self.TOLERANCE_LABELS[TOLERANCE_UNIT_FRAMES]
        )
        self.tolerance_label.text = label
        self.textInput_tolerance.placeholder = placeholder

    # ------------------------------------------------------------------
    # Activity: spinner, status line, buttons
    # ------------------------------------------------------------------
    def _show_status(self, message: str):
        """Write the status line (only when it changes: every Label write re-lays out the window)."""
        if getattr(self, "_closed", False):
            return
        if self.status_label.text != message:
            self.status_label.text = message

    def _set_busy(self, busy: bool, message: str):
        """Start or stop the spinner, and keep Browse and Run from starting a second job."""
        self._busy = busy
        if getattr(self, "_closed", False):
            return
        if busy:
            self.activity_indicator.start()
        else:
            self.activity_indicator.stop()
        self.browse_button.enabled = not busy
        self.btn_execute.enabled = not busy
        self._show_status(message)

    def _progress_reporter(self, loop) -> ProgressCallback:
        """A progress(message, force) callback for the worker thread.

        Messages are thinned to one per STATUS_INTERVAL (unless forced) and
        handed to the window's event loop, which alone touches widgets.
        """
        last = [0.0]

        def report(message: str, force: bool = False):
            now = time.monotonic()
            if not force and now - last[0] < STATUS_INTERVAL:
                return
            last[0] = now
            try:
                loop.call_soon_threadsafe(self._show_status, message)
            except RuntimeError:    # the event loop has closed (application exit)
                pass

        return report

    async def frames_counter(self, widget):
        """Browse: choose the trajectory and count its frames in a worker thread.

        Replaces FramesCounter.frames_counter for this tool only: the frames
        are counted from the newlines (count_frames), with the spinner running.
        """
        await self.open_file_dialog(widget)
        if not getattr(self, "trajec", None):
            return
        self.total_frame_number = 0     # a failed count must not leave the previous file's
        loop = asyncio.get_running_loop()
        report = self._progress_reporter(loop)
        self._set_busy(True, "Counting frames...")
        try:
            self.num_atoms, frames = await asyncio.to_thread(
                count_frames, self.trajec,
                lambda fraction: report(f"Counting frames... {fraction:.0%}"),
            )
        except (OSError, ValueError) as e:
            self._set_busy(False, "The trajectory could not be read.")
            await self.warning_function(f"Error reading TRAJEC.xyz file: {e}")
            return
        self._set_busy(False, f"{frames} frames of {self.num_atoms} atoms. Fill in the fields and "
                              "press 'MRT calculation'.")
        self.total_frame_number = frames

        with open(self.trajec, "r") as f:
            frame_lines = [f.readline() for _ in range(self.num_atoms + 2)]
        self.multi_line_text.value = (
            f"Number of atoms  -->  {self.num_atoms}\n"
            f"Number of lines  -->  {frames * (self.num_atoms + 2)}\n"
            f"Number of frames -->  {frames}\n"
            f"  TRAJEC.xyz FIRST FRAME\n"
            + "".join(frame_lines)
        )

    async def workflow(self, widget):
        params_are_valid = await self.read_params(widget)
        if not params_are_valid:
            return

        loop = asyncio.get_running_loop()
        started = time.monotonic()
        self._set_busy(True, "Reading the trajectory...")
        try:
            await asyncio.to_thread(self.run_analysis, self._progress_reporter(loop))
        except Exception as e:
            self._set_busy(False, "The MRT calculation stopped with an error.")
            if not getattr(self, "_closed", False):
                await self.warning_function(str(e))
            return
        self._set_busy(False, f"Finished in {time.monotonic() - started:.1f} s: "
                              f"{len(self.shell_atom_list)} {self.atom_symbol} atoms entered the shell.")
        if getattr(self, "_closed", False):
            return
        self.multi_line_text.value = self.results_text(self.center_atom_symbol)
        try:
            self.display_plots()
        except Exception as e:
            await self.warning_function(f"The figure could not be shown: {e}")

    def _on_close(self, window, **kwargs):
        self._closed = True   # a job still running finishes without touching the window
        return True

    def closeTopLevel(self, widget):
        self._closed = True
        self.main_window.close()
