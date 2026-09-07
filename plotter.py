import csv
import json
import os
import tempfile
from typing import List

import matplotlib.pyplot as plt
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from displayPlots import launch_plot_viewer
from help import HelpGqteaWin

class PlotterBase:
    CPMD_PLOT_TYPE = "CPMD energy file"
    GQTEAMD_PLOT_TYPE = "gqteaMD energy file"
    JSON_PLOT_TYPE = "JSON plot file"

    def __init__(self):
        # Initialize x_axis to avoid attribute error
        self.x_axis = []
        self.data = []
        self.gqtea_headers = []
        self.gqtea_rows = []
        self.gqtea_y_switches = []
        self.json_figures = []
        self.json_file = ""

    # ------------------------------------------------------------------ #
    # JSON plot manifests                                                  #
    #                                                                      #
    # The accepted format is the plot manifest every gQTEA analysis tool    #
    # already writes for the interactive viewer (displayPlots.save_plots),  #
    # so any figure produced elsewhere in the suite can be reopened here:   #
    #                                                                      #
    #   [{"x": [...], "y": [...], "xlabel": "...", "ylabel": "...",         #
    #     "title": "...", "xlim": [lo, hi], "ylim": [lo, hi]},              #
    #    {"series": [{"x": [...], "y": [...], "label": "..."}, ...], ...}]  #
    #                                                                      #
    # xlabel/ylabel/title/xlim/ylim are optional.                           #
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

            if "series" in figure:
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
            if "series" in figure:
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

    async def open_file_dialog(self, widget):
        try:
            is_gqtea = self.is_gqtea_plot_type()
            is_json = self.is_json_plot_type()
            if is_json:
                dialog_title = "Open JSON plot file"
                file_types = ["json"]
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

    async def json_plot_static(self):
        try:
            for figure in self.json_figures:
                temp_filename = tempfile.NamedTemporaryFile(
                    delete=False, suffix=".png", dir=self.output_dir
                ).name
                plt.figure(figsize=(8, 6))
                if "series" in figure:
                    for curve in figure["series"]:
                        plt.plot(curve.get("x", []), curve.get("y", []),
                                 label=curve.get("label", ""), antialiased=True)
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
            size=(760, 680),
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

        main_box.add(self.cpmd_switch_box)

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
        
        main_box.add(self.time_step_box)
        main_box.add(self.units_box)

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
        self.gqtea_column_box.add(self.gqtea_y_columns_box)
        main_box.add(self.gqtea_column_box)

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

        # Set the content of the main window
        self.main_window.content = main_box
        self.on_plot_type_change(self.plot_type_selection)
        self.main_window.show()

    async def workflow(self, widget):
        # JSON mode carries its own data and none of the switches below apply.
        if self.is_json_plot_type():
            await self.json_plot()
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

        # In JSON mode the figures are fully described by the file, so every
        # other plot-specific control is switched off.
        is_cpmd = not is_gqtea and not is_json

        if is_json:
            self.file_label.text = "Select JSON plot file:"
            self.text_input_file.placeholder = "Click Browse to select a JSON plot file"
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
        self.text_input_file.value = ""
        self.multi_line_text.value = (
            HelpGqteaWin.Json_Plot_Options if is_json else HelpGqteaWin.Plotting_Options
        )

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

