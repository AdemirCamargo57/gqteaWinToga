import csv
import json
import math
import os
import tempfile
from typing import List

import matplotlib.pyplot as plt
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from displayPlots import launch_plot_viewer, series_style_kwargs
from help import HelpGqteaWin

class PlotterBase:
    CPMD_PLOT_TYPE = "CPMD energy file"
    GQTEAMD_PLOT_TYPE = "gqteaMD energy file"
    JSON_PLOT_TYPE = "JSON plot file"
    TABLE_PLOT_TYPE = "Multi-column data file"

    def __init__(self):
        # Initialize x_axis to avoid attribute error
        self.x_axis = []
        self.data = []
        self.gqtea_headers = []
        self.gqtea_rows = []
        self.gqtea_y_switches = []
        self.json_figures = []
        self.json_file = ""
        self.table_data = None

    # ------------------------------------------------------------------ #
    # JSON plot manifests                                                  #
    #                                                                      #
    # The accepted format is the plot manifest every gQTEA analysis tool    #
    # already writes for the interactive viewer (displayPlots.save_plots),  #
    # so any figure produced elsewhere in the suite can be reopened here:   #
    #                                                                      #
    #   [{"x": [...], "y": [...], "xlabel": "...", "ylabel": "...",         #
    #     "title": "...", "xlim": [lo, hi], "ylim": [lo, hi]},              #
    #    {"series": [{"x": [...], "y": [...], "label": "..."}, ...], ...},  #
    #    {"type": "bars", "categories": [...],                              #
    #     "groups": [{"label": "...", "values": [...], "errors": [...]}],   #
    #     "line": {"label": "...", "values": [...], "ylabel": "..."}}]      #
    #                                                                      #
    # xlabel/ylabel/title/xlim/ylim are optional; so are a group's 'errors' #
    # and a bars figure's 'line'.                                           #
    # ------------------------------------------------------------------ #
    @staticmethod
    def _check_numeric_sequence(values, where: str) -> int:
        """Validate one x or y array and return its length."""
        if not isinstance(values, list):
            raise ValueError(f"{where} must be a list of numbers.")
        if not values:
            raise ValueError(f"{where} contains no data points.")
        for value in values:
            # bool is an int subclass, and NaN/inf cannot be plotted meaningfully
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{where} must contain only numbers.")
            if value != value or value in (float("inf"), float("-inf")):
                raise ValueError(f"{where} must contain only finite numbers.")
        return len(values)

    @classmethod
    def _check_limits(cls, figure: dict, key: str, where: str) -> None:
        if key not in figure:
            return
        limits = figure[key]
        if (not isinstance(limits, list) or len(limits) != 2
                or any(isinstance(v, bool) or not isinstance(v, (int, float))
                       for v in limits)):
            raise ValueError(f"{where}: '{key}' must be two numbers [low, high].")

    @classmethod
    def validate_plot_manifest(cls, payload) -> List[dict]:
        """Check a decoded JSON plot manifest and return its figures.

        Raises ValueError with a message naming the offending figure, so the
        user is told *which* entry is wrong rather than just that the file is
        bad.
        """
        if not isinstance(payload, list):
            raise ValueError(
                "The JSON file must contain a list of figures (a JSON array), "
                "for example: [{\"x\": [...], \"y\": [...]}]."
            )
        if not payload:
            raise ValueError("The JSON file contains no figures to plot.")

        for index, figure in enumerate(payload, start=1):
            where = f"Figure {index}"
            if not isinstance(figure, dict):
                raise ValueError(f"{where} is not a JSON object.")

            if figure.get("type") == "bars":
                categories = figure.get("categories")
                if not isinstance(categories, list) or not categories:
                    raise ValueError(
                        f"{where}: 'categories' must be a non-empty list."
                    )
                n_categories = len(categories)
                groups = figure.get("groups")
                if not isinstance(groups, list) or not groups:
                    raise ValueError(
                        f"{where}: 'groups' must be a non-empty list of bar groups."
                    )
                for g_index, group in enumerate(groups, start=1):
                    tag = f"{where}, group {g_index}"
                    if not isinstance(group, dict):
                        raise ValueError(f"{tag} is not a JSON object.")
                    if "values" not in group:
                        raise ValueError(f"{tag} needs a 'values' list.")
                    n_values = cls._check_numeric_sequence(group["values"], f"{tag}: 'values'")
                    if n_values != n_categories:
                        raise ValueError(
                            f"{tag}: 'values' must have the same length as "
                            f"'categories' ({n_values} vs {n_categories})."
                        )
                    if group.get("errors") is not None:
                        n_errors = cls._check_numeric_sequence(group["errors"], f"{tag}: 'errors'")
                        if n_errors != n_categories:
                            raise ValueError(
                                f"{tag}: 'errors' must have the same length as "
                                f"'categories' ({n_errors} vs {n_categories})."
                            )
                line = figure.get("line")
                if line is not None:
                    tag = f"{where}, line"
                    if not isinstance(line, dict) or "values" not in line:
                        raise ValueError(f"{tag} needs a 'values' list.")
                    n_line = cls._check_numeric_sequence(line["values"], f"{tag}: 'values'")
                    if n_line != n_categories:
                        raise ValueError(
                            f"{tag}: 'values' must have the same length as "
                            f"'categories' ({n_line} vs {n_categories})."
                        )
            elif "series" in figure:
                series = figure["series"]
                if not isinstance(series, list) or not series:
                    raise ValueError(
                        f"{where}: 'series' must be a non-empty list of curves."
                    )
                for s_index, curve in enumerate(series, start=1):
                    tag = f"{where}, series {s_index}"
                    if not isinstance(curve, dict):
                        raise ValueError(f"{tag} is not a JSON object.")
                    if "x" not in curve or "y" not in curve:
                        raise ValueError(f"{tag} needs both 'x' and 'y'.")
                    n_x = cls._check_numeric_sequence(curve["x"], f"{tag}: 'x'")
                    n_y = cls._check_numeric_sequence(curve["y"], f"{tag}: 'y'")
                    if n_x != n_y:
                        raise ValueError(
                            f"{tag}: 'x' and 'y' must have the same length "
                            f"({n_x} vs {n_y})."
                        )
                    # Optional styling; an unknown value would crash the viewer.
                    if "linestyle" in curve and curve["linestyle"] not in cls.SERIES_LINESTYLES:
                        raise ValueError(
                            f"{tag}: 'linestyle' must be one of "
                            f"{', '.join(cls.SERIES_LINESTYLES)}."
                        )
                    if "marker" in curve and curve["marker"] not in cls.SERIES_MARKERS:
                        raise ValueError(
                            f"{tag}: 'marker' must be one of "
                            f"{' '.join(cls.SERIES_MARKERS)}."
                        )
            elif "x" in figure and "y" in figure:
                n_x = cls._check_numeric_sequence(figure["x"], f"{where}: 'x'")
                n_y = cls._check_numeric_sequence(figure["y"], f"{where}: 'y'")
                if n_x != n_y:
                    raise ValueError(
                        f"{where}: 'x' and 'y' must have the same length "
                        f"({n_x} vs {n_y})."
                    )
            else:
                raise ValueError(
                    f"{where} must provide either 'x' and 'y', or 'series'."
                )

            cls._check_limits(figure, "xlim", where)
            cls._check_limits(figure, "ylim", where)

        return payload

    @classmethod
    def load_json_plot_file(cls, filepath: str) -> List[dict]:
        """Read and validate a JSON plot manifest from disk."""
        if not filepath or not os.path.isfile(filepath):
            raise ValueError(f"JSON file not found: {filepath}")

        try:
            with open(filepath, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"The file is not valid JSON (line {exc.lineno}, "
                f"column {exc.colno}): {exc.msg}."
            ) from exc
        except OSError as exc:
            raise ValueError(f"Could not read the JSON file: {exc}") from exc

        return cls.validate_plot_manifest(payload)

    @staticmethod
    def describe_json_figures(figures: List[dict]) -> str:
        """One readable line per figure, for the message panel."""
        lines = [f"Loaded {len(figures)} figure(s) from the JSON file:", ""]
        for index, figure in enumerate(figures, start=1):
            title = figure.get("title") or "(untitled)"
            if figure.get("type") == "bars":
                n_groups = len(figure.get("groups", []))
                n_categories = len(figure.get("categories", []))
                detail = f"{n_groups} groups, {n_categories} categories"
            elif "series" in figure:
                n_points = len(figure["series"][0].get("x", []))
                detail = f"{len(figure['series'])} curves, {n_points} points each"
            else:
                detail = f"{len(figure.get('x', []))} points"
            xlabel = figure.get("xlabel", "")
            ylabel = figure.get("ylabel", "")
            axes = f"  [{xlabel} vs {ylabel}]" if (xlabel or ylabel) else ""
            lines.append(f"  {index}. {title} - {detail}{axes}")
        lines.append("")
        lines.append("Press Plot to open them in the interactive viewer.")
        return "\n".join(lines)

    def is_json_plot_type(self):
        return self.plot_type_selection.value == self.JSON_PLOT_TYPE

    # ------------------------------------------------------------------ #
    # Multi-column tables                                                  #
    #                                                                      #
    # Column 1 is x; every further column is one y series plotted against  #
    # it. Blank lines and comment lines (TABLE_COMMENT_PREFIXES, plus any  #
    # trailing "# ..." on a line) are skipped wherever they appear. The    #
    # last non-comment, non-numeric line before the first data row is the  #
    # column-name header; earlier ones (e.g. a title) are ignored. Once    #
    # data has started, a non-numeric line is an error, never a header.    #
    # ------------------------------------------------------------------ #
    TABLE_COMMENT_PREFIXES = ("#", "!", "%", "@", "//")

    # Line-style choices offered per y column: label -> (matplotlib
    # linestyle, marker or None). Written into each manifest series as the
    # optional 'linestyle'/'marker' keys, which plotViewer passes to plot().
    TABLE_LINE_STYLES = {
        "Solid": ("-", None),
        "Dashed": ("--", None),
        "Dotted": (":", None),
        "Dash-dot": ("-.", None),
        "Line + markers": ("-", "o"),
        "Markers only": ("none", "o"),
    }
    DEFAULT_TABLE_LINE_STYLE = "Solid"
    # Values accepted for a series' optional 'linestyle'/'marker' keys.
    SERIES_LINESTYLES = ("-", "--", ":", "-.", "none")
    SERIES_MARKERS = ("o", "s", "^", "v", "D", "x", "+", "*", ".")

    @staticmethod
    def _split_table_line(line: str) -> List[str]:
        """Split on commas (CSV), else tabs (TSV), else any whitespace."""
        if "," in line:
            tokens = [token.strip() for token in line.split(",")]
        elif "\t" in line:
            tokens = [token.strip() for token in line.split("\t")]
        else:
            return line.split()
        # A trailing delimiter ("1,2,") leaves an empty last field.
        while tokens and not tokens[-1]:
            tokens.pop()
        return tokens

    @staticmethod
    def _parse_numeric_tokens(tokens: List[str]):
        """Return the tokens as floats, or None if any is not a number."""
        try:
            return [float(token) for token in tokens]
        except ValueError:
            return None

    @classmethod
    def parse_table_text(cls, text: str) -> dict:
        """Parse a multi-column numeric table.

        Returns {"x", "ys", "xlabel", "labels", "has_header", "n_rows"};
        raises ValueError naming the offending line when the table is unusable.
        """
        header = None
        header_line_no = None
        x_values: List[float] = []
        ys: List[List[float]] = []
        n_columns = None

        for line_no, raw in enumerate(text.splitlines(), start=1):
            line = raw.strip()
            if not line or line.startswith(cls.TABLE_COMMENT_PREFIXES):
                continue
            line = line.split("#", 1)[0].strip()  # trailing inline comment
            if not line:
                continue

            tokens = cls._split_table_line(line)
            values = cls._parse_numeric_tokens(tokens)

            if values is None:
                if n_columns is None:
                    header, header_line_no = tokens, line_no
                    continue
                raise ValueError(
                    f"Line {line_no} is not numeric: '{raw.strip()}'. Text is only "
                    "allowed before the data, as comments or a column-name header."
                )

            if any(math.isnan(v) or math.isinf(v) for v in values):
                raise ValueError(
                    f"Line {line_no} must contain only finite numbers: '{raw.strip()}'."
                )

            if n_columns is None:
                n_columns = len(values)
                if n_columns < 2:
                    raise ValueError(
                        f"The table needs at least two columns (x and one y); "
                        f"line {line_no} has {n_columns}."
                    )
                ys = [[] for _ in range(n_columns - 1)]
            elif len(values) != n_columns:
                raise ValueError(
                    f"Line {line_no} has {len(values)} values, but the table has "
                    f"{n_columns} columns."
                )

            x_values.append(values[0])
            for column, value in zip(ys, values[1:]):
                column.append(value)

        if n_columns is None:
            raise ValueError("The file contains no numeric data rows.")

        if header is not None:
            if len(header) != n_columns:
                raise ValueError(
                    f"The column-name header (line {header_line_no}) has "
                    f"{len(header)} names, but the data has {n_columns} columns."
                )
            xlabel, labels = header[0], header[1:]
        else:
            xlabel = "Column 1"
            labels = [f"Column {i}" for i in range(2, n_columns + 1)]

        return {
            "x": x_values,
            "ys": ys,
            "xlabel": xlabel,
            "labels": labels,
            "has_header": header is not None,
            "n_rows": len(x_values),
        }

    @classmethod
    def load_table_file(cls, filepath: str) -> dict:
        """Read and parse a multi-column table from disk."""
        if not filepath or not os.path.isfile(filepath):
            raise ValueError(f"Data file not found: {filepath}")
        try:
            with open(filepath, "r", encoding="utf-8", errors="replace") as handle:
                text = handle.read()
        except OSError as exc:
            raise ValueError(f"Could not read the data file: {exc}") from exc
        return cls.parse_table_text(text)

    @classmethod
    def table_to_figure(cls, table: dict, title: str = "", xlabel: str = "",
                        ylabel: str = "", styles: List[str] = None,
                        include: List[bool] = None) -> dict:
        """Express a parsed table as one plotViewer 'series' figure.

        Blank xlabel/ylabel fall back to the header's first name and "Value".
        ``styles`` holds one TABLE_LINE_STYLES key per y column (default Solid);
        ``include`` one flag per y column saying whether to plot it (default
        all). Both are indexed by column, so a style stays with its column
        whichever columns are left out.
        """
        n_series = len(table["ys"])
        if styles is None:
            styles = [cls.DEFAULT_TABLE_LINE_STYLE] * n_series
        if include is None:
            include = [True] * n_series
        if len(styles) != n_series:
            raise ValueError(
                f"The table has {n_series} y columns but {len(styles)} line "
                f"styles were given."
            )
        if len(include) != n_series:
            raise ValueError(
                f"The table has {n_series} y columns but {len(include)} plot "
                f"selections were given."
            )
        if not any(include):
            raise ValueError("Select at least one y column to plot.")

        series = []
        for y, label, style, wanted in zip(table["ys"], table["labels"], styles, include):
            if not wanted:
                continue
            if style not in cls.TABLE_LINE_STYLES:
                raise ValueError(f"Unknown line style '{style}' for '{label}'.")
            linestyle, marker = cls.TABLE_LINE_STYLES[style]
            curve = {"x": table["x"], "y": y, "label": label, "linestyle": linestyle}
            if marker:
                curve["marker"] = marker
            series.append(curve)

        return {
            "series": series,
            "xlabel": (xlabel or "").strip() or table["xlabel"],
            "ylabel": (ylabel or "").strip() or "Value",
            "title": title,
        }

    @staticmethod
    def describe_table(table: dict) -> str:
        """Readable summary of a loaded table, for the message panel."""
        source = "column-name header" if table["has_header"] else "no header (generic names)"
        lines = [
            f"Loaded {table['n_rows']} data rows, {len(table['ys'])} y-axis series.",
            f"Column names: {source}",
            "",
            f"  x-axis: {table['xlabel']}",
        ]
        lines += [f"  y{i}: {label}" for i, label in enumerate(table["labels"], start=1)]
        lines += ["", "Press Plot to open the figure in the interactive viewer."]
        return "\n".join(lines)

    def is_table_plot_type(self):
        return self.plot_type_selection.value == self.TABLE_PLOT_TYPE

    async def open_file_dialog(self, widget):
        try:
            is_gqtea = self.is_gqtea_plot_type()
            is_json = self.is_json_plot_type()
            is_table = self.is_table_plot_type()
            if is_json:
                dialog_title = "Open JSON plot file"
                file_types = ["json"]
            elif is_table:
                dialog_title = "Open multi-column data file"
                file_types = ["dat", "txt", "csv", "xvg"]
            elif is_gqtea:
                dialog_title = "Open gqteaMD energy file"
                file_types = ["*.csv", "*.dat", "*.log", "*.txt", "*.*"]
            else:
                dialog_title = "Open CPMD ENERGY File"
                file_types = ["*.*", "*.txt"]

            self.energy_file = await self.main_window.dialog(
                toga.OpenFileDialog(
                    title=dialog_title,
                    multiple_select=False,
                    file_types=file_types,
                )
            )

            if not self.energy_file:
                await self.main_window.info_dialog("Warning", "No file was selected!")
                return

            self.energy_file = str(self.energy_file)
            self.text_input_file.value = self.energy_file
            self.output_dir = os.path.dirname(self.energy_file)
            self.data = []
            if is_json:
                await self.parse_json_plot_file()
            elif is_table:
                await self.parse_table_file()
            elif is_gqtea:
                await self.parse_gqtea_energy_file()
            else:
                await self.parse_energy_file()

        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Failed to open file: {e}")
            )

    async def parse_energy_file(self):
        try:
            with open(self.energy_file, "r") as f:
                for idx, line in enumerate(f, start=1):
                    line_data = line.strip().split()
                    if len(line_data) != 8:
                        self.multi_line_text.value = f"ERROR: Line {idx} has invalid format: {line}\n"
                        await self.main_window.dialog(toga.ErrorDialog("ERROR",f'Line {idx} has invalid format: {line}'))
                        return

                    self.data.append(line_data)

            self.multi_line_text.value = f"Number of steps: {len(self.data)}\n"
            self.multi_line_text.value += "The ENERGY file has a valid format to be displayed.\n"

            time_step_str = self.text_input_time_step.value
            if not time_step_str:
                await self.main_window.info_dialog("Error", "Simulation time step is required.")
                return

            try:
                time_step = float(time_step_str)
            except ValueError:
                await self.main_window.info_dialog("Error", "Invalid format for simulation time step.")
                return

            self.compute_x_axis(time_step)

        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Failed to read file: {e}")
            )

    async def parse_gqtea_energy_file(self):
        try:
            delimiter = await self.detect_gqtea_delimiter()
            if delimiter is None:
                return

            if delimiter == "whitespace":
                headers, rows = await self.read_gqtea_whitespace_energy_file()
            else:
                headers, rows = await self.read_gqtea_csv_energy_file()

            if headers is None:
                return

            if not rows:
                await self.main_window.dialog(
                    toga.ErrorDialog("ERROR", "The gqteaMD energy file contains no data rows.")
                )
                return

            self.gqtea_headers = headers
            self.gqtea_rows = rows
            self.data = rows
            self.update_gqtea_column_controls(headers)
            self.multi_line_text.value = f"Number of steps: {len(rows)}\n"
            delimiter_label = "whitespace-separated" if delimiter == "whitespace" else "CSV"
            self.multi_line_text.value += (
                f"The gqteaMD energy file has a valid {delimiter_label} format to be displayed.\n"
            )

        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Failed to read file: {e}")
            )

    async def detect_gqtea_delimiter(self):
        with open(self.energy_file, "r") as f:
            for line in f:
                stripped_line = line.strip()
                if stripped_line:
                    return "csv" if "," in stripped_line else "whitespace"

        await self.main_window.dialog(
            toga.ErrorDialog("ERROR", "The gqteaMD energy file is empty.")
        )
        return None

    async def read_gqtea_csv_energy_file(self):
        with open(self.energy_file, "r", newline="") as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames:
                await self.main_window.dialog(
                    toga.ErrorDialog("ERROR", "The gqteaMD energy file has no header row.")
                )
                return None, None

            header_pairs = [(header, header.strip()) for header in reader.fieldnames]
            headers = [clean_header for _, clean_header in header_pairs]
            if not await self.validate_gqtea_headers(headers):
                return None, None

            rows = []
            for idx, row in enumerate(reader, start=2):
                if row.get(None):
                    await self.main_window.dialog(
                        toga.ErrorDialog(
                            "ERROR",
                            f"Line {idx} has too many values for the CSV header.",
                        )
                    )
                    return None, None

                cleaned_row = {}
                for original_header, header in header_pairs:
                    value = row.get(original_header, "")
                    try:
                        cleaned_row[header] = float(value)
                    except (TypeError, ValueError):
                        await self.main_window.dialog(
                            toga.ErrorDialog(
                                "ERROR",
                                f"Line {idx}, column '{header}' is not numeric: {value}",
                            )
                        )
                        return None, None
                rows.append(cleaned_row)

        return headers, rows

    async def read_gqtea_whitespace_energy_file(self):
        with open(self.energy_file, "r") as f:
            lines = [(idx, line.strip()) for idx, line in enumerate(f, start=1) if line.strip()]

        if not lines:
            await self.main_window.dialog(
                toga.ErrorDialog("ERROR", "The gqteaMD energy file is empty.")
            )
            return None, None

        _, header_line = lines[0]
        headers = header_line.split()
        if not await self.validate_gqtea_headers(headers):
            return None, None

        rows = []
        for idx, line in lines[1:]:
            values = line.split()
            if len(values) != len(headers):
                await self.main_window.dialog(
                    toga.ErrorDialog(
                        "ERROR",
                        f"Line {idx} has {len(values)} values, but the header line has {len(headers)} columns.",
                    )
                )
                return None, None

            cleaned_row = {}
            for header, value in zip(headers, values):
                try:
                    cleaned_row[header] = float(value)
                except ValueError:
                    await self.main_window.dialog(
                        toga.ErrorDialog(
                            "ERROR",
                            f"Line {idx}, column '{header}' is not numeric: {value}",
                        )
                    )
                    return None, None
            rows.append(cleaned_row)

        return headers, rows

    async def validate_gqtea_headers(self, headers):
        if any(not header for header in headers):
            await self.main_window.dialog(
                toga.ErrorDialog("ERROR", "The gqteaMD energy file contains an empty column name.")
            )
            return False
        if len(set(headers)) != len(headers):
            await self.main_window.dialog(
                toga.ErrorDialog("ERROR", "The gqteaMD energy file contains duplicate column names.")
            )
            return False
        return True

    def is_gqtea_plot_type(self):
        return self.plot_type_selection.value == self.GQTEAMD_PLOT_TYPE

    async def parse_json_plot_file(self):
        """Load the selected JSON manifest, or report exactly what is wrong."""
        self.json_figures = []
        self.json_file = ""
        try:
            figures = self.load_json_plot_file(self.energy_file)
        except ValueError as exc:
            self.multi_line_text.value = (
                f"Could not load the JSON plot file.\n\n{exc}"
            )
            await self.main_window.dialog(
                toga.ErrorDialog("Invalid JSON plot file", str(exc))
            )
            return

        self.json_figures = figures
        self.json_file = self.energy_file
        self.multi_line_text.value = self.describe_json_figures(figures)

    async def json_plot(self):
        """Show the loaded manifest, preferring the interactive viewer."""
        if not self.json_figures:
            await self.main_window.dialog(
                toga.InfoDialog(
                    "No data",
                    "Load a JSON plot file first, using the Browse button.",
                )
            )
            return

        # The manifest on disk is already exactly what plotViewer consumes, so
        # it can be handed over untouched (zoom/pan/save toolbar, and it works
        # in frozen builds too).
        if launch_plot_viewer(self.json_file):
            return

        # Fallback: render each figure to a PNG and show it in a Toga window,
        # matching how the other plot types in this module display figures.
        await self.json_plot_static()

    async def parse_table_file(self):
        """Load the selected multi-column table, or report what is wrong."""
        self.table_data = None
        self.update_table_style_controls([])
        try:
            table = self.load_table_file(self.energy_file)
        except ValueError as exc:
            self.multi_line_text.value = f"Could not load the data file.\n\n{exc}"
            await self.main_window.dialog(
                toga.ErrorDialog("Invalid data file", str(exc))
            )
            return

        self.table_data = table
        self.update_table_style_controls(table["labels"])
        self.multi_line_text.value = self.describe_table(table)

    async def table_plot(self):
        """Plot every y column against column 1, preferring the interactive viewer."""
        if not self.table_data:
            await self.main_window.dialog(
                toga.InfoDialog(
                    "No data",
                    "Load a multi-column data file first, using the Browse button.",
                )
            )
            return

        try:
            figure = self.table_to_figure(
                self.table_data,
                title=os.path.basename(self.energy_file),
                xlabel=self.table_xlabel_input.value,
                ylabel=self.table_ylabel_input.value,
                styles=self.selected_table_styles(),
                include=self.selected_table_columns(),
            )
        except ValueError as exc:
            await self.main_window.dialog(toga.ErrorDialog("Plot Error", str(exc)))
            return
        # plotViewer reads a manifest from disk; keep it in the temp dir so
        # nothing is left next to the user's data.
        with tempfile.NamedTemporaryFile(
            "w", delete=False, suffix=".json", prefix="gqtea_table_plot_",
            encoding="utf-8",
        ) as handle:
            json.dump([figure], handle)
            manifest_path = handle.name

        if launch_plot_viewer(manifest_path):
            return
        await self.render_figures_static([figure])

    async def json_plot_static(self):
        await self.render_figures_static(self.json_figures)

    async def render_figures_static(self, figures: List[dict]):
        """Fallback: draw manifest figures to PNGs shown in Toga windows."""
        try:
            for figure in figures:
                temp_filename = tempfile.NamedTemporaryFile(
                    delete=False, suffix=".png", dir=self.output_dir
                ).name
                figure_obj = plt.figure(figsize=(8, 6))
                if figure.get("type") == "bars":
                    # Deferred import: displayPlots must not be pulled in at
                    # module load time by the --plot-viewer child process (see
                    # displayPlots.py / plotViewer.py); this method only runs in
                    # the Toga (main) process, where the import is harmless.
                    from displayPlots import draw_bar_comparison
                    draw_bar_comparison(figure_obj.gca(), figure)
                else:
                    if "series" in figure:
                        for curve in figure["series"]:
                            plt.plot(curve.get("x", []), curve.get("y", []),
                                     label=curve.get("label", ""), antialiased=True,
                                     **series_style_kwargs(curve))
                        plt.legend()
                    else:
                        plt.plot(figure.get("x", []), figure.get("y", []),
                                 antialiased=True)
                    plt.xlabel(figure.get("xlabel", ""))
                    plt.ylabel(figure.get("ylabel", ""))
                    plt.title(figure.get("title", ""))
                    if "xlim" in figure:
                        plt.xlim(*figure["xlim"])
                    if "ylim" in figure:
                        plt.ylim(*figure["ylim"])
                plt.tight_layout()
                plt.savefig(temp_filename)
                plt.close()
                self.show_plot(temp_filename)
        except Exception as exc:
            await self.main_window.dialog(
                toga.ErrorDialog("Plot Error", f"Could not draw the figures: {exc}")
            )

    def compute_x_axis(self, time_step: float):
        x_values = [float(row[0]) for row in self.data]

        unit = self.unit_selection.value
        if unit == "steps":
            self.x_axis = x_values
        elif unit == "fs":
            self.x_axis = [x * time_step * 0.024188 for x in x_values]
        elif unit == "ps":
            self.x_axis = [(x * time_step * 0.024188) / 1000.0 for x in x_values]
        else:
            self.main_window.info_dialog("Error", "Unsupported x-axis unit selected.")
            self.x_axis = []

    def display_plot(self, plot_title: str, x_label: str, y_label: str, lines: List[dict]):

        # Save the temporary plot file
        temp_filename = tempfile.NamedTemporaryFile(
            delete=False, suffix=".png", dir=self.output_dir
        ).name

        if not self.x_axis:
            self.main_window.info_dialog("Error", "X-axis data is not available.")
            return

        plt.figure(figsize=(8, 6))
        for line in lines:
            plt.plot(self.x_axis, line["data"], label=line["label"], antialiased=True)

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(plot_title)
        plt.legend()
        plt.tight_layout()

        # Save plot to in-memory buffer
        plt.savefig(temp_filename)
        plt.close()

        # Display the plot using Toga's ImageView
        self.show_plot(temp_filename)

    def show_plot(self, temp_filename):
        # Load the image using Toga's Image class
        plot_image = toga.Image(temp_filename)
        # Create a new window for the plot
        plot_window = toga.Window(
            title="Plot",
            size=(900, 600),
        )
        # Create a box to hold the image
        plot_box = toga.Box(style=Pack(direction=COLUMN,flex=1))
        plot_window.content = plot_box
        # Create an ImageView to display the image
        plot_imageview = toga.ImageView(plot_image, style=Pack(flex=1, margin=10))
        plot_box.add(plot_imageview)
        # Show the window
        plot_window.show()


    def fictitious_and_ionic(self, widget=None):
        time_step = float(self.text_input_time_step.value)
        
        self.compute_x_axis(time_step)
        
        y_fictitious = [float(row[1]) for row in self.data]
        y_ionic = [float(row[4]) - float(row[3]) for row in self.data]

        lines = [
            {"data": y_fictitious, "label": "Fictitious Energy"},
            {"data": y_ionic, "label": "Ionic Kinetic Energy"}
        ]

        self.display_plot(
            plot_title="Fictitious and Ionic Kinetic Energy",
            x_label=self.get_x_label(),
            y_label="Energy (Ha)",
            lines=lines
        )

    def temperature(self, widget=None):
        time_step = float(self.text_input_time_step.value)       
        self.compute_x_axis(time_step)
        
        y_temperature = [float(row[2]) for row in self.data]

        lines = [
            {"data": y_temperature, "label": "Temperature (K)"}
        ]

        self.display_plot(
            plot_title="Simulation Temperature",
            x_label=self.get_x_label(),
            y_label="Temperature (K)",
            lines=lines
        )

    def khon_sham_energy(self, widget=None):
        time_step = float(self.text_input_time_step.value)       
        self.compute_x_axis(time_step)
        
        y_ksh_energy = [float(row[3]) for row in self.data]

        lines = [
            {"data": y_ksh_energy, "label": "Kohn-Sham Energy (Ha)"}
        ]

        self.display_plot(
            plot_title="Kohn-Sham Potential Energy",
            x_label=self.get_x_label(),
            y_label="Energy (Ha)",
            lines=lines
        )

    def kohn_sham_and_ionic(self, widget=None):
        time_step = float(self.text_input_time_step.value)       
        self.compute_x_axis(time_step)
        
        y_ionic_kinetic = [float(row[4]) for row in self.data]

        lines = [
            {"data": y_ionic_kinetic, "label": "Ionic Kinetic Energy (Ha)"}
        ]

        self.display_plot(
            plot_title="Kohn-Sham and Ionic Kinetic Energy",
            x_label=self.get_x_label(),
            y_label="Energy (Ha)",
            lines=lines
        )

    def total_energy(self, widget=None):
        time_step = float(self.text_input_time_step.value)       
        self.compute_x_axis(time_step)
        
        y_total_energy = [float(row[5]) for row in self.data]

        lines = [
            {"data": y_total_energy, "label": "Total Energy (Ha)"}
        ]

        self.display_plot(
            plot_title="Total Energy",
            x_label=self.get_x_label(),
            y_label="Energy (Ha)",
            lines=lines
        )

    def cpu_time(self, widget=None):
        time_step = float(self.text_input_time_step.value)       
        self.compute_x_axis(time_step)
        
        y_cpu_time = [float(row[7]) for row in self.data]

        lines = [
            {"data": y_cpu_time, "label": "CPU Time (s)"}
        ]

        self.display_plot(
            plot_title="CPU Time by Step",
            x_label=self.get_x_label(),
            y_label="Time (s)",
            lines=lines
        )

    def get_x_label(self) -> str:
        unit = self.unit_selection.value
        return "Steps" if unit == "steps" else f"Time ({unit})"

    async def gqtea_energy_plot(self):
        if not self.gqtea_rows:
            await self.main_window.info_dialog("Error", "No gqteaMD data available. Please load a valid data file.")
            return

        x_column = self.gqtea_x_axis_selection.value
        y_columns = [
            switch.text
            for switch in self.gqtea_y_switches
            if switch.value and switch.enabled
        ]

        if not x_column:
            await self.main_window.info_dialog("Info", "Please select an x-axis column.")
            return

        if not y_columns:
            await self.main_window.info_dialog("Info", "Please select at least one y-axis column.")
            return

        self.x_axis = [row[x_column] for row in self.gqtea_rows]
        lines = [
            {"data": [row[column] for row in self.gqtea_rows], "label": column}
            for column in y_columns
        ]

        self.display_plot(
            plot_title="gqteaMD Energy",
            x_label=x_column,
            y_label="Value",
            lines=lines,
        )


class PlotterUI(PlotterBase):
    def __init__(self,*args):
        super().__init__()
        self.layout_main_window(*args)

    def layout_main_window(self,widget):
        # Create the main window
        self.main_window = toga.Window(
            title="Energy File Plot",
            size=(760, 620),
        )

        # Define common styles
        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        button_style = Pack(margin=5, width=100)

        # Main container
        main_box = toga.Box(style=Pack(direction=COLUMN, margin=20))

        # Title
        title_label = toga.Label("Plot", style=heading_style)
        main_box.add(title_label)

        # Plot type selection
        plot_type_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0,5,10,5)))
        plot_type_label = toga.Label(
            "Plot type:",
            style=Pack(margin=(0,5,0,5), text_align=LEFT, width=200),
        )
        self.plot_type_selection = toga.Selection(
            items=[
                self.CPMD_PLOT_TYPE,
                self.GQTEAMD_PLOT_TYPE,
                self.JSON_PLOT_TYPE,
                self.TABLE_PLOT_TYPE,
            ],
            on_change=self.on_plot_type_change,
            style=Pack(flex=1, margin=(0,5,0,5)),
        )
        plot_type_box.add(plot_type_label)
        plot_type_box.add(self.plot_type_selection)
        main_box.add(plot_type_box)

        # Switches for plot selection
        self.cpmd_switch_box = toga.Box(
            style=Pack(direction=COLUMN, align_items=CENTER, margin=(0, 0, 5, 5))
        )

        self.switch_fictitious_ionic = toga.Switch(
            "Fictitious and Ionic Kinetic Energy Plot",
            style=Pack(text_align=LEFT, margin=(0,5,0,5),)
        )

        self.switch_temperature = toga.Switch(
            "Temperature Plot",
            style=Pack(text_align=LEFT, margin=(0,5,0,5),)
        )

        self.switch_ksh_energy = toga.Switch(
            "Potential Kohn-Sham Energy Plot",
            style=Pack(text_align=LEFT, margin=(0,5,0,5)),
        )

        self.switch_ksh_ionic = toga.Switch(
            "Kohn-Sham Plus Ionic Kinetic Energy",
            style=Pack(text_align=LEFT, margin=(0,5,0,5)),
        )

        self.switch_total_energy = toga.Switch(
            "Total Energy",
            style=Pack(text_align=LEFT, margin=(0,5,0,5)),
        )

        self.switch_cpu_time = toga.Switch(
            "CPU Time Plot",
            style=Pack(text_align=LEFT, margin=(0,5,0,5)),
        )

        switches = [
            self.switch_fictitious_ionic,
            self.switch_temperature,
            self.switch_ksh_energy,
            self.switch_ksh_ionic,
            self.switch_total_energy,
            self.switch_cpu_time
        ]

        for switch in switches:
            self.cpmd_switch_box.add(switch)

        self._section_slots = {}
        self._add_section(main_box, "cpmd_switch_box")

        # File selection section
        file_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0,5,0,5)))

        self.file_label = toga.Label(
            "Select cpmd ENERGY file:",
            style=Pack(margin=(0,5,0,5), text_align=LEFT, width=200),
        )

        self.text_input_file = toga.TextInput(
            placeholder="Click Browse to select CPMD ENERGY file",
            style=Pack(flex=1, margin=(0,5,0,5), color="blue"),
        )

        browse_button = toga.Button(
            "Browse", on_press=self.open_file_dialog, style=button_style
        )

        file_box.add(self.file_label)
        file_box.add(self.text_input_file)
        file_box.add(browse_button)
        main_box.add(file_box)

        # Simulation time step input
        self.time_step_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0,5,0,5)))

        self.time_step_label = toga.Label(
            "Simulation Time Step:",
            style=Pack(margin=(0,5,0,5), text_align=LEFT, width=200),
        )

        self.text_input_time_step = toga.TextInput(
            placeholder=" ",
            style=Pack(flex=1, margin=(5,5,0,5)),
        )
        self.text_input_time_step.value = 5.0
        
        self.time_step_box.add(self.time_step_label)
        self.time_step_box.add(self.text_input_time_step)
        

        # X-axis unit selection
        self.units_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(5,5,10,5)))

        self.units_label = toga.Label(
            "Select x-axis Unit: ",
            style=Pack(text_align=LEFT, width=200),
        )

        self.unit_selection = toga.Selection(
            items=["steps", "fs", "ps"],
            style=Pack(margin=(0,5,5,13)),
        )

        self.units_box.add(self.units_label)
        self.units_box.add(self.unit_selection)
        
        self._add_section(main_box, "time_step_box")
        self._add_section(main_box, "units_box")

        # gqteaMD column selection
        self.gqtea_column_box = toga.Box(
            style=Pack(direction=COLUMN, margin=(0,5,10,5))
        )

        self.gqtea_x_axis_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin=(0,5,0,5))
        )
        self.gqtea_x_axis_label = toga.Label(
            "Data x-axis column:",
            style=Pack(margin=(0,5,0,5), text_align=LEFT, width=200),
        )
        self.gqtea_x_axis_control_box = toga.Box(style=Pack(direction=ROW, flex=1))
        self.gqtea_x_axis_selection = toga.Selection(
            items=["Load a data file"],
            on_change=self.on_gqtea_x_axis_change,
            style=Pack(flex=1, margin=(0,5,0,5)),
        )
        self.gqtea_x_axis_control_box.add(self.gqtea_x_axis_selection)
        self.gqtea_x_axis_box.add(self.gqtea_x_axis_label)
        self.gqtea_x_axis_box.add(self.gqtea_x_axis_control_box)

        self.gqtea_y_columns_label = toga.Label(
            "Data y-axis columns:",
            style=Pack(margin=(0,5,0,5), text_align=LEFT),
        )
        self.gqtea_y_columns_box = toga.Box(style=Pack(direction=COLUMN, margin=(0,5,0,205)))

        self.gqtea_column_box.add(self.gqtea_x_axis_box)
        self.gqtea_column_box.add(self.gqtea_y_columns_label)
        # One switch per file column: scroll, so a wide file cannot grow the window.
        self.gqtea_column_box.add(toga.ScrollContainer(
            content=self.gqtea_y_columns_box,
            horizontal=False,
            style=Pack(height=110, margin=(0,5,0,5)),
        ))
        self._add_section(main_box, "gqtea_column_box")

        # Multi-column table: axis labels and one line style per y column
        self.table_options_box = toga.Box(
            style=Pack(direction=COLUMN, margin=(0,5,10,5))
        )
        self.table_xlabel_input = toga.TextInput(
            placeholder="blank = name of the first column",
            style=Pack(flex=1, margin=(0,5,0,5)),
        )
        self.table_ylabel_input = toga.TextInput(
            placeholder="blank = Value",
            style=Pack(flex=1, margin=(0,5,0,5)),
        )
        for text, field in (("Table x-axis label:", self.table_xlabel_input),
                            ("Table y-axis label:", self.table_ylabel_input)):
            row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0,5,0,5)))
            row.add(toga.Label(text, style=Pack(margin=(0,5,0,5), text_align=LEFT, width=200)))
            row.add(field)
            self.table_options_box.add(row)

        self.table_styles_label = toga.Label(
            "Y columns to plot and their line styles:",
            style=Pack(margin=(5,5,0,10), text_align=LEFT),
        )
        self.table_styles_box = toga.Box(style=Pack(direction=COLUMN))
        self.table_style_selections = []
        self.table_column_switches = []
        # Scrolls so a table with many columns cannot push the buttons off-screen.
        table_styles_scroll = toga.ScrollContainer(
            content=self.table_styles_box,
            horizontal=False,
            style=Pack(height=110, margin=(0,5,0,5)),
        )
        self.table_options_box.add(self.table_styles_label)
        self.table_options_box.add(table_styles_scroll)
        self._add_section(main_box, "table_options_box")

        # Multi-line text for messages
        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(0,5,0,5), font_size=12),
            readonly=True,
        )
        self.multi_line_text.value = HelpGqteaWin.Plotting_Options
        main_box.add(self.multi_line_text)

        # Action buttons
        button_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin_top=10)
        )

        self.btn_help = toga.Button(
            "Help",style=button_style, on_press=self.open_window_help)

        self.btn_execute = toga.Button(
            "Plot", style=button_style, on_press=self.workflow
        )
        self.btn_close = toga.Button(
            "Close", style=button_style, on_press=self.close_window
        )
        button_box.add(self.btn_execute)
        button_box.add(self.btn_help)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        # Trim to the active type's sections *before* setting the content: with
        # every section present the window would first grow to fit them all.
        self.on_plot_type_change(self.plot_type_selection)
        self.main_window.content = main_box
        self.main_window.show()

    async def workflow(self, widget):
        # JSON mode carries its own data and none of the switches below apply.
        if self.is_json_plot_type():
            await self.json_plot()
            return

        if self.is_table_plot_type():
            await self.table_plot()
            return

        if not self.data:
            await self.main_window.info_dialog("Error", "No data available. Please load a valid ENERGY file.")
            return

        if self.is_gqtea_plot_type():
            await self.gqtea_energy_plot()
            return

        if not any([
            self.switch_fictitious_ionic.value,
            self.switch_temperature.value,
            self.switch_ksh_energy.value,
            self.switch_ksh_ionic.value,
            self.switch_total_energy.value,
            self.switch_cpu_time.value
        ]):
            await self.main_window.info_dialog("Info", "Please select at least one plot option.")
            return

        if self.switch_fictitious_ionic.value:
            self.fictitious_and_ionic()

        if self.switch_temperature.value:
            self.temperature()

        if self.switch_ksh_energy.value:
            self.khon_sham_energy()

        if self.switch_ksh_ionic.value:
            self.kohn_sham_and_ionic()

        if self.switch_total_energy.value:
            self.total_energy()

        if self.switch_cpu_time.value:
            self.cpu_time()

    def on_plot_type_change(self, widget):
        is_gqtea = self.is_gqtea_plot_type()
        is_json = self.is_json_plot_type()
        is_table = self.is_table_plot_type()

        # In JSON and table modes the figure is fully described by the file,
        # so every other plot-specific control is switched off.
        is_cpmd = not is_gqtea and not is_json and not is_table

        if is_json:
            self.file_label.text = "Select JSON plot file:"
            self.text_input_file.placeholder = "Click Browse to select a JSON plot file"
        elif is_table:
            self.file_label.text = "Select multi-column data file:"
            self.text_input_file.placeholder = "Click Browse to select a data table (x, y1, y2, ...)"
        elif is_gqtea:
            self.file_label.text = "Select gqteaMD data file:"
            self.text_input_file.placeholder = "Click Browse to select gqteaMD energy file"
        else:
            self.file_label.text = "Select cpmd ENERGY file:"
            self.text_input_file.placeholder = "Click Browse to select CPMD ENERGY file"

        for switch in self.get_cpmd_switches():
            if not is_cpmd:
                switch.value = False
            switch.enabled = is_cpmd

        self.text_input_time_step.enabled = is_cpmd
        self.unit_selection.enabled = is_cpmd
        self.time_step_label.enabled = is_cpmd
        self.units_label.enabled = is_cpmd

        self.gqtea_x_axis_label.enabled = is_gqtea
        self.gqtea_x_axis_selection.enabled = is_gqtea and bool(self.gqtea_headers)
        self.gqtea_y_columns_label.enabled = is_gqtea
        for switch in self.gqtea_y_switches:
            switch.enabled = is_gqtea and switch.text != self.gqtea_x_axis_selection.value

        self.data = []
        self.x_axis = []
        self.json_figures = []
        self.json_file = ""
        self.table_data = None
        self.text_input_file.value = ""
        if is_json:
            self.multi_line_text.value = HelpGqteaWin.Json_Plot_Options
        elif is_table:
            self.multi_line_text.value = HelpGqteaWin.Table_Plot_Options
        else:
            self.multi_line_text.value = HelpGqteaWin.Plotting_Options

        # Table customisation belongs to the loaded file, so it is reset too.
        self.table_xlabel_input.value = ""
        self.table_ylabel_input.value = ""
        self.update_table_style_controls([])
        self.table_xlabel_input.enabled = is_table
        self.table_ylabel_input.enabled = is_table
        self.table_styles_label.enabled = is_table

        # Only the active type's section takes up room: stacking all of them
        # (disabled) made the window taller than a laptop screen.
        self.show_sections({
            "cpmd_switch_box": is_cpmd,
            "time_step_box": is_cpmd,
            "units_box": is_cpmd,
            "gqtea_column_box": is_gqtea,
            "table_options_box": is_table,
        })

    def _add_section(self, main_box, name):
        """Place section ``name`` in its own slot box, so it can later be
        removed and restored at the same position in the window."""
        slot = toga.Box(style=Pack(direction=COLUMN))
        slot.add(getattr(self, name))
        main_box.add(slot)
        self._section_slots[name] = slot

    def show_sections(self, shown_by_name):
        # Toga 0.5's Pack layout ignores display="none" (a hidden widget still
        # takes its full height), so a hidden section is taken out of its
        # slot; an empty slot has zero height. Remove before adding: otherwise
        # both types' sections coexist for a moment and the window grows to
        # fit them, and Toga never shrinks it back.
        for name, shown in shown_by_name.items():
            section, slot = getattr(self, name), self._section_slots[name]
            if not shown and section in slot.children:
                slot.remove(section)
        for name, shown in shown_by_name.items():
            section, slot = getattr(self, name), self._section_slots[name]
            if shown and section not in slot.children:
                slot.add(section)

    def update_table_style_controls(self, labels):
        """Rebuild one row per y column (in column order): a Switch choosing
        whether to plot it and a line-style Selection."""
        for child in list(self.table_styles_box.children):
            self.table_styles_box.remove(child)
        self.table_style_selections = []
        self.table_column_switches = []
        for label in labels:
            row = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0,5,0,5)))
            selection = toga.Selection(
                items=list(self.TABLE_LINE_STYLES),
                style=Pack(flex=1, margin=(0,5,0,5)),
            )
            selection.value = self.DEFAULT_TABLE_LINE_STYLE
            switch = toga.Switch(
                label,
                value=True,
                on_change=lambda widget, sel=selection: setattr(sel, "enabled", widget.value),
                style=Pack(margin=(0,5,0,5), text_align=LEFT, width=190),
            )
            row.add(switch)
            row.add(selection)
            self.table_styles_box.add(row)
            self.table_column_switches.append(switch)
            self.table_style_selections.append(selection)

    def selected_table_styles(self):
        return [selection.value for selection in self.table_style_selections]

    def selected_table_columns(self):
        return [bool(switch.value) for switch in self.table_column_switches]

    def get_cpmd_switches(self):
        return [
            self.switch_fictitious_ionic,
            self.switch_temperature,
            self.switch_ksh_energy,
            self.switch_ksh_ionic,
            self.switch_total_energy,
            self.switch_cpu_time,
        ]

    def update_gqtea_column_controls(self, headers):
        for child in list(self.gqtea_x_axis_control_box.children):
            self.gqtea_x_axis_control_box.remove(child)

        self.gqtea_x_axis_selection = toga.Selection(
            items=headers,
            on_change=self.on_gqtea_x_axis_change,
            style=Pack(flex=1, margin=(0,5,0,5)),
        )
        self.gqtea_x_axis_control_box.add(self.gqtea_x_axis_selection)

        for child in list(self.gqtea_y_columns_box.children):
            self.gqtea_y_columns_box.remove(child)

        self.gqtea_y_switches = []
        for header in headers:
            switch = toga.Switch(
                header,
                style=Pack(text_align=LEFT, margin=(0,5,0,5)),
            )
            self.gqtea_y_switches.append(switch)
            self.gqtea_y_columns_box.add(switch)

        self.gqtea_x_axis_selection.enabled = self.is_gqtea_plot_type()
        self.on_gqtea_x_axis_change(self.gqtea_x_axis_selection)

    def on_gqtea_x_axis_change(self, widget):
        x_column = self.gqtea_x_axis_selection.value
        is_gqtea = self.is_gqtea_plot_type()
        for switch in self.gqtea_y_switches:
            if switch.text == x_column:
                switch.value = False
                switch.enabled = False
            else:
                switch.enabled = is_gqtea
            

    def open_window_help(self, widget):

        window = toga.Window(title=f"Instructions to use plot mudule",
                             size = (700, 600),)
        
        help_box = toga.Box(style=Pack(direction=COLUMN, flex=1))
        multi_line_text = toga.MultilineTextInput(
            style=Pack(font_size=11, margin=(5, 5), flex=1)
        )
        multi_line_text.value = HelpGqteaWin.help_plots

        help_box.add(multi_line_text)

        window.content = help_box

        window.show()        
              
    def close_window(self, widget):
        self.main_window.close()

