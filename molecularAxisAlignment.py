from __future__ import annotations

import math
import os
from typing import List

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW


# Recenters the molecule, finds its principal spatial axes, and rewrites the XYZ
# coordinates so the smallest, middle, and largest extents align with x, y, and z.
class MolecularAxisAlignment:
    """Align an XYZ structure so its principal extents map to x, y, and z."""

    def __init__(self) -> None:
        self.input_file = ""
        self.output_file = ""

    @staticmethod
    def _normalize_dialog_result(dialog_result) -> str:
        if not dialog_result:
            return ""
        if isinstance(dialog_result, (list, tuple)):
            if not dialog_result:
                return ""
            dialog_result = dialog_result[0]
        path_value = getattr(dialog_result, "path", dialog_result)
        return str(path_value)

    @staticmethod
    def suggest_output_filename(input_file: str) -> str:
        base, ext = os.path.splitext(input_file)
        if not ext:
            ext = ".xyz"
        return f"{base}_aligned{ext}"

    def read_xyz(self, file_path: str) -> tuple[List[str], str, List[List[float]]]:
        with open(file_path, "r", encoding="utf-8") as file_handle:
            lines = [line.rstrip("\n") for line in file_handle]

        if len(lines) < 2:
            raise ValueError("The XYZ file must contain at least two header lines.")

        try:
            atom_count = int(lines[0].strip())
        except ValueError as exc:
            raise ValueError("The first line of the XYZ file must be an integer atom count.") from exc

        coordinate_lines = [line for line in lines[2:] if line.strip()]
        if len(coordinate_lines) < atom_count:
            raise ValueError("The XYZ file does not contain the expected number of atom records.")

        symbols: List[str] = []
        coordinates: List[List[float]] = []
        for index, line in enumerate(coordinate_lines[:atom_count], start=1):
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"Invalid XYZ atom line {index + 2}: '{line}'")
            symbol = parts[0]
            try:
                x_coord = float(parts[1])
                y_coord = float(parts[2])
                z_coord = float(parts[3])
            except ValueError as exc:
                raise ValueError(f"Invalid coordinates on line {index + 2}: '{line}'") from exc

            symbols.append(symbol)
            coordinates.append([x_coord, y_coord, z_coord])

        return symbols, lines[1], coordinates

    def write_xyz(
        self,
        file_path: str,
        symbols: List[str],
        comment: str,
        coordinates: List[List[float]],
    ) -> None:
        with open(file_path, "w", encoding="utf-8") as file_handle:
            file_handle.write(f"{len(symbols)}\n")
            file_handle.write(f"{comment}\n")
            for symbol, (x_coord, y_coord, z_coord) in zip(symbols, coordinates):
                file_handle.write(f"{symbol:<3} {x_coord: .8f} {y_coord: .8f} {z_coord: .8f}\n")

    @staticmethod
    def _centroid(coordinates: List[List[float]]) -> List[float]:
        atom_count = len(coordinates)
        if atom_count == 0:
            raise ValueError("At least one atom is required for alignment.")
        return [
            sum(coord[axis] for coord in coordinates) / atom_count
            for axis in range(3)
        ]

    @staticmethod
    def _subtract(point_a: List[float], point_b: List[float]) -> List[float]:
        return [point_a[index] - point_b[index] for index in range(3)]

    @staticmethod
    def _dot(vector_a: List[float], vector_b: List[float]) -> float:
        return sum(a * b for a, b in zip(vector_a, vector_b))

    @staticmethod
    def _cross(vector_a: List[float], vector_b: List[float]) -> List[float]:
        return [
            vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1],
            vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2],
            vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0],
        ]

    @staticmethod
    def _norm(vector: List[float]) -> float:
        return math.sqrt(sum(value * value for value in vector))

    def _normalize(self, vector: List[float]) -> List[float]:
        length = self._norm(vector)
        if length == 0.0:
            raise ValueError("Cannot normalize a zero-length vector.")
        return [value / length for value in vector]

    def _covariance_matrix(self, centered_coordinates: List[List[float]]) -> List[List[float]]:
        covariance = [[0.0, 0.0, 0.0] for _ in range(3)]
        for point in centered_coordinates:
            for row in range(3):
                for column in range(3):
                    covariance[row][column] += point[row] * point[column]

        scale = 1.0 / max(len(centered_coordinates), 1)
        for row in range(3):
            for column in range(3):
                covariance[row][column] *= scale
        return covariance

    def _jacobi_eigen_decomposition(
        self,
        matrix: List[List[float]],
        tolerance: float = 1.0e-12,
        max_iterations: int = 100,
    ) -> tuple[List[float], List[List[float]]]:
        eigenvectors = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
        working = [row[:] for row in matrix]

        for _ in range(max_iterations):
            largest_value = 0.0
            pivot_row = 0
            pivot_column = 1
            for row in range(3):
                for column in range(row + 1, 3):
                    candidate = abs(working[row][column])
                    if candidate > largest_value:
                        largest_value = candidate
                        pivot_row, pivot_column = row, column

            if largest_value < tolerance:
                break

            if abs(working[pivot_row][pivot_row] - working[pivot_column][pivot_column]) < tolerance:
                angle = math.pi / 4.0
            else:
                angle = 0.5 * math.atan2(
                    2.0 * working[pivot_row][pivot_column],
                    working[pivot_column][pivot_column] - working[pivot_row][pivot_row],
                )

            cosine = math.cos(angle)
            sine = math.sin(angle)

            for row in range(3):
                if row != pivot_row and row != pivot_column:
                    value_row = working[row][pivot_row]
                    value_column = working[row][pivot_column]
                    working[row][pivot_row] = cosine * value_row - sine * value_column
                    working[pivot_row][row] = working[row][pivot_row]
                    working[row][pivot_column] = sine * value_row + cosine * value_column
                    working[pivot_column][row] = working[row][pivot_column]

            diagonal_pp = working[pivot_row][pivot_row]
            diagonal_qq = working[pivot_column][pivot_column]
            off_diagonal = working[pivot_row][pivot_column]

            working[pivot_row][pivot_row] = (
                cosine * cosine * diagonal_pp
                - 2.0 * sine * cosine * off_diagonal
                + sine * sine * diagonal_qq
            )
            working[pivot_column][pivot_column] = (
                sine * sine * diagonal_pp
                + 2.0 * sine * cosine * off_diagonal
                + cosine * cosine * diagonal_qq
            )
            working[pivot_row][pivot_column] = 0.0
            working[pivot_column][pivot_row] = 0.0

            for row in range(3):
                eigenvector_row_p = eigenvectors[row][pivot_row]
                eigenvector_row_q = eigenvectors[row][pivot_column]
                eigenvectors[row][pivot_row] = cosine * eigenvector_row_p - sine * eigenvector_row_q
                eigenvectors[row][pivot_column] = sine * eigenvector_row_p + cosine * eigenvector_row_q

        eigenvalues = [working[index][index] for index in range(3)]
        eigenvector_columns = [
            self._normalize([eigenvectors[row][column] for row in range(3)])
            for column in range(3)
        ]
        return eigenvalues, eigenvector_columns

    def _deterministic_axis_sign(
        self,
        axis_vector: List[float],
        centered_coordinates: List[List[float]],
    ) -> List[float]:
        projections = [self._dot(point, axis_vector) for point in centered_coordinates]
        if not projections:
            return axis_vector
        max_index = max(range(len(projections)), key=lambda idx: abs(projections[idx]))
        if projections[max_index] < 0.0:
            return [-value for value in axis_vector]
        return axis_vector

    def _principal_axes(self, centered_coordinates: List[List[float]]) -> List[List[float]]:
        if len(centered_coordinates) == 1:
            return [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]

        covariance = self._covariance_matrix(centered_coordinates)
        eigenvalues, eigenvectors = self._jacobi_eigen_decomposition(covariance)

        axes_by_extent = sorted(
            zip(eigenvalues, eigenvectors),
            key=lambda item: item[0],
        )

        x_axis = self._deterministic_axis_sign(axes_by_extent[0][1], centered_coordinates)
        y_axis = self._deterministic_axis_sign(axes_by_extent[1][1], centered_coordinates)

        z_axis = self._cross(x_axis, y_axis)
        if self._norm(z_axis) < 1.0e-10:
            z_axis = axes_by_extent[2][1]
        z_axis = self._deterministic_axis_sign(self._normalize(z_axis), centered_coordinates)

        # Rebuild y from z and x to keep an orthonormal right-handed basis.
        y_candidate = self._cross(z_axis, x_axis)
        if self._norm(y_candidate) < 1.0e-10:
            y_axis = self._deterministic_axis_sign(axes_by_extent[1][1], centered_coordinates)
            y_axis = self._normalize(y_axis)
        else:
            y_axis = self._normalize(y_candidate)

        if self._dot(self._cross(x_axis, y_axis), z_axis) < 0.0:
            x_axis = [-value for value in x_axis]

        return [x_axis, y_axis, z_axis]

    def align_coordinates(self, coordinates: List[List[float]]) -> List[List[float]]:
        centroid = self._centroid(coordinates)
        centered_coordinates = [self._subtract(point, centroid) for point in coordinates]
        axes = self._principal_axes(centered_coordinates)

        aligned_coordinates: List[List[float]] = []
        for point in centered_coordinates:
            aligned_coordinates.append([
                self._dot(point, axes[0]),
                self._dot(point, axes[1]),
                self._dot(point, axes[2]),
            ])
        return aligned_coordinates

    def align_xyz_file(self, input_file: str, output_file: str | None = None) -> str:
        symbols, comment, coordinates = self.read_xyz(input_file)
        aligned_coordinates = self.align_coordinates(coordinates)

        final_output = output_file or self.suggest_output_filename(input_file)
        updated_comment = f"{comment} | aligned to principal spatial axes"
        self.write_xyz(final_output, symbols, updated_comment, aligned_coordinates)
        return final_output


class MolecularAxisAlignmentUI(MolecularAxisAlignment):
    """Toga window for XYZ principal-axis alignment."""

    def __init__(self, *args) -> None:
        del args
        super().__init__()
        self.main_window = toga.Window(title="Molecular Axis Alignment", size=(720, 260))
        self.input_file_input: toga.TextInput | None = None
        self.output_file_input: toga.TextInput | None = None
        self.status_label: toga.Label | None = None
        self.build_ui()
        self.main_window.show()

    def build_ui(self) -> None:
        label_style = Pack(width=140, margin_right=8)
        row_style = Pack(direction=ROW, margin_bottom=10)

        self.input_file_input = toga.TextInput(
            placeholder="Select an XYZ structure file",
            style=Pack(flex=1, margin_right=8),
        )
        self.output_file_input = toga.TextInput(
            placeholder="Suggested: molecule_aligned.xyz",
            style=Pack(flex=1, margin_right=8),
        )
        self.status_label = toga.Label(
            "Largest extent -> z, middle extent -> y, smallest extent -> x.",
            style=Pack(margin_top=10, color="#555555"),
        )

        content = toga.Box(
            style=Pack(direction=COLUMN, margin=18),
            children=[
                toga.Label(
                    "Align XYZ Structure to Principal Spatial Axes",
                    style=Pack(font_size=16, font_weight="bold", margin_bottom=12),
                ),
                toga.Box(
                    style=row_style,
                    children=[
                        toga.Label("Input file", style=label_style),
                        self.input_file_input,
                        toga.Button("Browse", on_press=self.select_input_file),
                    ],
                ),
                toga.Box(
                    style=row_style,
                    children=[
                        toga.Label("Output file", style=label_style),
                        self.output_file_input,
                        toga.Button("Save As", on_press=self.select_output_file),
                    ],
                ),
                toga.Button(
                    "Align Molecule",
                    on_press=self.run_alignment,
                    style=Pack(width=160, margin_top=8),
                ),
                self.status_label,
            ],
        )

        self.main_window.content = content

    async def select_input_file(self, widget) -> None:
        del widget
        try:
            selected = await self.main_window.dialog(
                toga.OpenFileDialog(title="Select XYZ File")
            )
            input_file = self._normalize_dialog_result(selected)
            if not input_file:
                return
            assert self.input_file_input is not None
            assert self.output_file_input is not None
            self.input_file = input_file
            self.input_file_input.value = input_file
            self.output_file = self.suggest_output_filename(input_file)
            self.output_file_input.value = self.output_file
        except Exception as exc:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to open file: {exc}"))

    async def select_output_file(self, widget) -> None:
        del widget
        try:
            assert self.output_file_input is not None
            suggested_name = (
                os.path.basename(self.output_file_input.value.strip())
                if self.output_file_input.value.strip()
                else "molecule_aligned.xyz"
            )
            selected = await self.main_window.dialog(
                toga.SaveFileDialog(
                    title="Save Aligned XYZ File",
                    suggested_filename=suggested_name,
                )
            )
            output_file = self._normalize_dialog_result(selected)
            if output_file:
                self.output_file = output_file
                self.output_file_input.value = output_file
        except Exception as exc:
            await self.main_window.dialog(toga.ErrorDialog("Error", f"Failed to select output file: {exc}"))

    async def run_alignment(self, widget) -> None:
        del widget
        assert self.input_file_input is not None
        assert self.output_file_input is not None
        assert self.status_label is not None

        input_file = self.input_file_input.value.strip()
        output_file = self.output_file_input.value.strip()

        if not input_file:
            await self.main_window.dialog(toga.ErrorDialog("Input Error", "Please choose an XYZ input file."))
            return

        if not output_file:
            output_file = self.suggest_output_filename(input_file)
            self.output_file_input.value = output_file

        try:
            final_output = self.align_xyz_file(input_file, output_file)
            self.status_label.text = f"Aligned structure saved as: {final_output}"
            self.status_label.style.color = "#1e6b34"
        except Exception as exc:
            self.status_label.text = f"Alignment failed: {exc}"
            self.status_label.style.color = "#9f1239"
