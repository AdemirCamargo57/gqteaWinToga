"""Shared infrastructure for the single-item geometry analysis tools.

[bond.py](bond.py) (one atom pair), [bondAngle.py](bondAngle.py) (one i-j-k
triplet) and [dihedralAngle.py](dihedralAngle.py) (one i-j-k-l quadruplet) all
do the same thing to a different number of atoms: stream ``TRAJEC.xyz``, pull
out a handful of selected atoms per frame, turn them into one scalar per frame,
histogram it, derive a free energy from the histogram, and write ``.dat`` /
``.csv`` / summary files.

Everything identical between those three lives here so the three tools cannot
drift apart in behaviour, user feedback, output formatting or error handling.
Anything specific to a coordinate (how the scalar is computed, its units, its
histogram range, whether it has a Jacobian factor) stays in the tool.

These tools report progress **as text** in their ``MultilineTextInput``: none
of them has a progress bar or a frame-count label. ``frames_counter`` below
overrides the shared ``FramesCounter.frames_counter`` for that reason; the
mixin itself is left untouched, so every other tool still uses ``progress_label``.
"""
import asyncio
import math
import os

import numpy as np
import toga


class GeometryAnalysisBase:
    """Mixin holding the machinery common to bond / angle / dihedral tools.

    Combine it with ``FramesCounter`` and ``DisplayPlots``, e.g.::

        class BondAnalyser(GeometryAnalysisBase, FramesCounter, DisplayPlots):
            tool_title = "BOND LENGTH ANALYSIS"

    It must come before ``FramesCounter`` in the MRO so its text-reporting
    ``frames_counter`` wins.
    """

    atufs: float = 0.02418884326505   # Atomic time unit in femtoseconds.
    gas_constant: float = 0.001987204  # kcal/(mol*K)

    # Frames between status-text refreshes. Each refresh also yields to the Toga
    # event loop so the window keeps repainting during a long read.
    progress_stride: int = 200

    # Overridden per tool.
    tool_title: str = "GEOMETRY ANALYSIS"

    # ------------------------------------------------------------------ #
    # Dialogs and output directory                                        #
    # ------------------------------------------------------------------ #
    async def warning_function(self, msg: str) -> None:
        await self.main_window.dialog(toga.InfoDialog("Error", f"{msg}"))

    async def fail(self, message: str) -> bool:
        """Report a validation error and return False, for ``return await self.fail(...)``."""
        await self.warning_function(message)
        return False

    def _ensure_output_dir(self) -> None:
        """Ensure self.output_dir exists; default to the trajectory's folder."""
        if not getattr(self, "output_dir", None):
            if getattr(self, "trajec", None):
                self.output_dir = os.path.dirname(os.path.abspath(self.trajec)) or os.getcwd()
            else:
                self.output_dir = os.getcwd()
        os.makedirs(self.output_dir, exist_ok=True)

    async def _write_text(self, path: str, text: str) -> bool:
        """Write ``text`` to ``path``, reporting OS errors instead of raising.

        Without this an unwritable output directory raises straight out of an
        async Toga handler, which fails silently for the user.
        """
        try:
            with open(path, "w") as f:
                f.write(text)
        except OSError as e:
            await self.warning_function(
                f"Could not write {os.path.basename(path)}: {e}"
            )
            return False
        return True

    # ------------------------------------------------------------------ #
    # Text status reporting (no progress bar, no frame-count label)        #
    # ------------------------------------------------------------------ #
    def _status_header(self) -> str:
        labels = " ".join(str(v) for v in getattr(self, "atom_labels", []))
        return (
            f"{self.tool_title}\n"
            f"Trajectory: {os.path.basename(self.trajec)}\n"
            f"Selected atoms: {labels}\n\n"
        )

    def _set_status(self, message: str) -> None:
        """Report progress as text in the status box."""
        self.multi_line_text.value = self._status_header() + message

    async def frames_counter(self, widget):
        """Load the trajectory and count its frames, reporting progress as text.

        Overrides the shared ``FramesCounter.frames_counter`` so these tools can
        report in the ``MultilineTextInput`` instead of a progress label.
        """
        await self.open_file_dialog(widget)
        if not getattr(self, "trajec", None):
            return

        traj_name = os.path.basename(self.trajec)
        self.multi_line_text.value = f"Reading trajectory: {traj_name} ...\n"
        frame_count = 0
        try:
            with open(self.trajec, "r") as f:
                while True:
                    title_line = f.readline()
                    if not title_line.strip():
                        break
                    if not f.readline():  # comment line
                        break
                    for _ in range(self.num_atoms):
                        atom_line = f.readline().split()
                        if not atom_line or len(atom_line) != 4:
                            raise ValueError(
                                f"malformed atom line in frame {frame_count}"
                            )
                    frame_count += 1
                    if frame_count % 400 == 0:
                        self.multi_line_text.value = (
                            f"Reading trajectory: {traj_name} ...\n"
                            f"  frames: {frame_count}"
                        )
                        await asyncio.sleep(0)
            self.total_frame_number = frame_count
            self.multi_line_text.value = (
                f"Trajectory loaded: {traj_name}\n"
                f"  atoms: {self.num_atoms}   frames: {frame_count}\n"
            )
        except Exception as e:
            await self.main_window.dialog(
                toga.InfoDialog("Error", f"Error reading TRAJEC.xyz file: {e}")
            )

    # ------------------------------------------------------------------ #
    # Parameter reading helpers                                           #
    # ------------------------------------------------------------------ #
    async def read_value(self, text_input, field_name: str, expected_type):
        """Read one required field, reporting blanks and bad formats."""
        value = text_input.value.strip()
        if not value:
            await self.warning_function(
                f"Please input a valid value for {field_name}."
            )
            return None
        try:
            return expected_type(value)
        except ValueError as e:
            await self.warning_function(f"Invalid format for {field_name}: {e}")
            return None

    async def read_atom_labels(self, text_input, count: int):
        """Read exactly ``count`` distinct 1-based atom labels."""
        value = text_input.value.strip()
        if not value:
            await self.warning_function(
                f"Please input {count} atom labels."
            )
            return None
        try:
            labels = [int(v) for v in value.replace(",", " ").split()]
        except ValueError as e:
            await self.warning_function(f"Invalid format for atom labels: {e}")
            return None
        if len(labels) != count:
            await self.warning_function(
                f"Please input exactly {count} atom labels."
            )
            return None
        if any(label <= 0 for label in labels):
            await self.warning_function(
                "Atom labels must be positive integers starting at 1."
            )
            return None
        num_atoms = getattr(self, "num_atoms", 0)
        if num_atoms and any(label > num_atoms for label in labels):
            await self.warning_function(
                f"Atom labels must be between 1 and {num_atoms}."
            )
            return None
        if len(set(labels)) != count:
            await self.warning_function("The atom labels must all be different.")
            return None
        return labels

    async def read_cell_lengths(self, text_input):
        """Read the optional orthorhombic cell.

        Returns ``(ok, lengths)``: a blank field is valid and yields
        ``(True, None)`` meaning no periodic boundaries.
        """
        text = text_input.value.strip()
        if not text:
            return True, None
        try:
            lengths = [float(v) for v in text.replace(",", " ").split()]
        except ValueError as e:
            await self.warning_function(f"Invalid format for cell lengths a b c: {e}")
            return False, None
        if len(lengths) != 3:
            await self.warning_function(
                "Enter exactly three cell lengths: a b c (or leave blank)."
            )
            return False, None
        if any(length <= 0 for length in lengths):
            await self.warning_function("Cell lengths must be positive.")
            return False, None
        return True, lengths

    def require_trajectory(self) -> bool:
        """True when a trajectory has been loaded with Browse."""
        return bool(
            getattr(self, "trajec", None)
            and hasattr(self, "num_atoms")
            and hasattr(self, "total_frame_number")
        )

    # ------------------------------------------------------------------ #
    # Geometry helpers                                                    #
    # ------------------------------------------------------------------ #
    def _time_increment_ps(self) -> float:
        """Increment of simulation time between stored frames, in ps."""
        return self.dt * self.sampling * self.atufs / 1000.0

    def minimum_image(self, delta):
        """Apply the orthorhombic minimum-image convention to a displacement.

        ``delta`` is a length-3 sequence. With ``self.cell_lengths`` unset the
        displacement is returned unchanged. The convention is exact only for
        true separations up to ``min(a, b, c) / 2``; callers that cannot bound
        this in advance should check with ``min_image_limit``. Coordinates need
        not be wrapped into the box.
        """
        box = getattr(self, "cell_lengths", None)
        if box is None:
            return np.asarray(delta, dtype=float)
        d = np.asarray(delta, dtype=float)
        b = np.asarray(box, dtype=float)
        return d - b * np.round(d / b)

    @property
    def min_image_limit(self) -> float:
        """Half the smallest cell length, or infinity without periodic boundaries."""
        box = getattr(self, "cell_lengths", None)
        return math.inf if box is None else min(box) / 2.0

    @staticmethod
    def _remaining_content(handle) -> bool:
        """True when anything other than blank lines is left in the file."""
        for line in handle:
            if line.strip():
                return True
        return False

    # ------------------------------------------------------------------ #
    # Trajectory reader                                                   #
    # ------------------------------------------------------------------ #
    async def read_selected_atoms(self, rows, status=None):
        """Stream the trajectory, returning only the atoms at ``rows``.

        ``rows`` are 0-based positions within a frame. Returns a list with one
        entry per complete frame, each entry a list of ``(element, x, y, z)`` in
        the order given by ``rows`` -- or ``None`` if the read failed outright.

        Only the requested atoms are parsed: every other atom line is consumed
        with a bare ``readline()`` and never split or float-converted, which is
        what makes the read cheap on large systems.

        A blank line where a frame header belongs is treated as clean EOF (the
        trailing newline nearly every writer leaves). A truncated frame, an
        unreadable header, or a header whose atom count disagrees with
        ``self.num_atoms`` stops the read with a warning but **keeps** the
        frames already read, rather than discarding the whole run.
        """
        wanted = set(rows)
        total_frames = getattr(self, "total_frame_number", 0)
        frames = []

        try:
            with open(self.trajec, "r") as f:
                while True:
                    header = f.readline()
                    if not header.strip():
                        # EOF, or the trailing blank line at the end of the file.
                        # Say so only if real content follows: that means the
                        # frame framing has desynced.
                        if header and self._remaining_content(f):
                            await self.warning_function(
                                f"Unexpected blank line after frame {len(frames)}; "
                                f"the rest of the file was not read."
                            )
                        break

                    try:
                        declared = int(header.split()[0])
                    except (ValueError, IndexError):
                        await self.warning_function(
                            f"Frame {len(frames) + 1} has an unreadable header line; "
                            f"stopping after {len(frames)} frame(s)."
                        )
                        break
                    if declared != self.num_atoms:
                        await self.warning_function(
                            f"Frame {len(frames) + 1} declares an atom count of "
                            f"{declared}, but the file was opened with "
                            f"{self.num_atoms}; stopping after {len(frames)} frame(s)."
                        )
                        break

                    if not f.readline():  # comment line
                        await self.warning_function(
                            f"Frame {len(frames) + 1} is truncated; using the "
                            f"{len(frames)} complete frame(s) read so far."
                        )
                        break

                    selected = {}
                    truncated = False
                    for row in range(self.num_atoms):
                        line = f.readline()
                        if not line.strip():
                            truncated = True
                            break
                        if row not in wanted:
                            continue  # not needed: skip without parsing
                        tokens = line.split()
                        if len(tokens) < 4:
                            raise ValueError(
                                f"Invalid atom line format in frame {len(frames) + 1}."
                            )
                        selected[row] = (
                            tokens[0],
                            float(tokens[1]),
                            float(tokens[2]),
                            float(tokens[3]),
                        )

                    if truncated:
                        await self.warning_function(
                            f"Frame {len(frames) + 1} is truncated; using the "
                            f"{len(frames)} complete frame(s) read so far."
                        )
                        break

                    frames.append([selected[row] for row in rows])

                    if len(frames) % self.progress_stride == 0 and status is not None:
                        self._set_status(status(len(frames), total_frames))
                        await asyncio.sleep(0)  # let the window repaint

        except OSError as e:
            await self.warning_function(f"Could not open the trajectory file: {e}")
            return None
        except Exception as e:
            await self.warning_function(f"Error reading TRAJEC.xyz file: {e}")
            return None

        if not frames:
            await self.warning_function(
                "No frames were read from the trajectory file."
            )
            return None
        return frames

    # ------------------------------------------------------------------ #
    # Histogram, statistics and free energy                               #
    # ------------------------------------------------------------------ #
    def histogram_percent(self, values, bin_width: float, lower: float, upper: float):
        """Floor-bin ``values`` into ``[lower, upper)`` and return percentages.

        Returns ``(centers, percentages, total, out_of_range)``. Binning is
        ``floor((v - lower) / bin_width)``, matching the original scalar loops
        exactly. Values outside the range are excluded and counted rather than
        clamped into the end bin, so an out-of-range tail cannot silently
        inflate the last bar (the caller is expected to report the count).
        """
        values = np.asarray(values, dtype=float)
        num_bins = max(1, int((upper - lower) / bin_width))
        centers = [lower + bin_width * (i + 0.5) for i in range(num_bins)]

        indices = np.floor((values - lower) / bin_width).astype(int)
        in_range = (indices >= 0) & (indices < num_bins)
        counts = np.bincount(indices[in_range], minlength=num_bins)
        total = int(counts.sum())
        out_of_range = int((~in_range).sum())

        if total == 0:
            return centers, [0.0] * num_bins, 0, out_of_range
        percentages = [(c / total) * 100.0 for c in counts.tolist()]
        return centers, percentages, total, out_of_range

    @staticmethod
    def basic_stats(values) -> dict:
        """Mean / extremes / sample variance (ddof=1, matching statistics.variance)."""
        arr = np.asarray(values, dtype=float)
        return {
            "average": float(arr.mean()),
            "largest": float(arr.max()),
            "smallest": float(arr.min()),
            "variance": float(np.var(arr, ddof=1)) if arr.size > 1 else 0.0,
            "std_dev": float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0,
        }

    def _smoothing_window(self, n_points: int) -> int:
        """Choose a small odd-sized moving-average window for FE smoothing."""
        if n_points < 5:
            return 1
        window = max(5, min(11, n_points // 12))
        if window % 2 == 0:
            window += 1
        return min(window, n_points if n_points % 2 == 1 else n_points - 1)

    def smooth_series(self, values):
        """Gentle moving-average smoother with edge padding."""
        n_points = len(values)
        window = self._smoothing_window(n_points)
        if window <= 1:
            return list(values)
        pad = window // 2
        padded = np.pad(np.asarray(values, dtype=float), pad_width=pad, mode="edge")
        kernel = np.ones(window, dtype=float) / window
        return np.convolve(padded, kernel, mode="valid").tolist()

    def free_energy_from_histogram(self, centers, percentages, jacobian=None):
        """``-RT ln P`` per populated bin, optionally divided by a Jacobian factor.

        ``jacobian`` is a callable mapping a bin centre to its volume-element
        factor (``r**2`` for a distance, ``sin(theta)`` for an angle); ``None``
        means the coordinate has a uniform measure, as a dihedral does. Empty
        bins are skipped because ``ln 0`` is undefined. When ``self.smooth_fe``
        is set the probabilities are smoothed before the logarithm.

        Returns a list of ``[centre, free energy]`` pairs.
        """
        probabilities = []
        xs = []
        for x_i, pct in zip(centers, percentages):
            if pct <= 0.0:
                continue
            p = pct / 100.0
            if jacobian is not None:
                p = p / max(jacobian(x_i), 1e-12)
            xs.append(x_i)
            probabilities.append(max(p, 1e-12))

        if not probabilities:
            return []

        if getattr(self, "smooth_fe", False):
            probabilities = self.smooth_series(probabilities)

        rt = self.gas_constant * self.sim_temp
        return [
            [x_i, -rt * math.log(max(p, 1e-12))]
            for x_i, p in zip(xs, probabilities)
        ]

    # ------------------------------------------------------------------ #
    # Summary rendering                                                   #
    # ------------------------------------------------------------------ #
    def render_summary(self, pair_text: str, parameters, results) -> str:
        """Build the closing summary shared by all three tools.

        ``parameters`` and ``results`` are lists of ``(label, value)`` pairs.
        Keeping the layout here is what makes the three summary files and status
        boxes look the same.
        """
        # Width 22 keeps a visible gap after the longest label any of the three
        # tools uses ("Jacobian sin(theta):", "Wrap to [-180,180]:").
        col = 22
        lines = [
            f"{self.tool_title} - COMPLETED\n",
            "-" * (len(self.tool_title) + 12) + "\n",
            f"{'Trajectory:':<{col}}{os.path.basename(self.trajec)}\n",
            f"{'Selected atoms:':<{col}}{pair_text}\n",
            f"{'Frames processed:':<{col}}{self.frames_used}\n",
            f"{'Number of atoms:':<{col}}{self.num_atoms}\n",
            "\nPARAMETERS\n",
        ]
        lines += [f"{label + ':':<{col}}{value}\n" for label, value in parameters]
        lines.append("\nRESULTS\n")
        lines += [f"{label + ':':<{col}}{value}\n" for label, value in results]
        lines.append(f"\n{'Output directory:':<{col}}{self.output_dir}\n")
        return "".join(lines)
