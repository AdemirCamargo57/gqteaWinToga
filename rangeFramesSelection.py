import os
import asyncio
import numpy as np
import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from help import HelpGqteaWin

class RangeFramesSelection:
    def __init__(self):
        self.trajec = None
        self.num_atoms = 0
        self.total_frame_number = 0
        self.r_interval = None
        self.min_r = None
        self.max_r = None
        self.atom_labels = None
        self.output_dir = None
        self.fout = None
        self.main_window = None

    def _set_progress(self, percent):
        self.progress_bar.value = max(0.0, min(100.0, percent))

    @staticmethod
    def _format_xyz_frame(num_atoms, frame_number, atom_data, comment_suffix=""):
        comment_line = f"Frame {frame_number}"
        if comment_suffix:
            comment_line = f"{comment_line} | {comment_suffix}"

        frame_lines = [f"{num_atoms}\n", f"{comment_line}\n"]
        for atom in atom_data:
            frame_lines.append(
                f"{atom[0]} {atom[1]:>14.6f} {atom[2]:>14.6f} {atom[3]:>14.6f}\n"
            )
        return "".join(frame_lines)

    @staticmethod
    def _write_summary_file(
        summary_path,
        average_distance,
        standard_deviation,
        total_frames,
        selected_frames,
    ):
        with open(summary_path, "w", encoding="utf-8") as summary_file:
            summary_file.write("Interatomic distance range selection summary\n")
            if average_distance is None:
                summary_file.write("Average distance: N/A\n")
            else:
                summary_file.write(f"Average distance: {average_distance:.6f} A\n")
            if standard_deviation is None:
                summary_file.write("Standard deviation: N/A\n")
            else:
                summary_file.write(f"Standard deviation: {standard_deviation:.6f} A\n")
            summary_file.write(f"Original total number of frames: {total_frames}\n")
            summary_file.write(f"Number of selected frames: {selected_frames}\n")

    async def warning_function(self, msg):
        await self.main_window.dialog(
            toga.InfoDialog("Error", str(msg))
        )

    async def read_params(self, widget):
        async def read_input(text_input, field_name, expected_type):
            value = text_input.value.strip()
            if not value:
                await self.warning_function(f"Please input a valid value for {field_name}.")
                return None
            try:
                if expected_type == list:
                    # Expects two space-separated atom labels
                    labels = [int(label) for label in value.split()]
                    if len(labels) != 2:
                        raise ValueError("Please input exactly two atom labels.")
                    return labels
                else:
                    return expected_type(value)
            except ValueError as e:
                await self.warning_function(f"Invalid format for {field_name}: {e}")
                return None

        async def read_interval(text_input, field_name):
            raw_value = text_input.value.strip()
            if not raw_value:
                await self.warning_function(f"Please input a valid value for {field_name}.")
                return None

            normalized_value = raw_value
            for separator in ("to", ",", ";"):
                normalized_value = normalized_value.replace(separator, " ")
            normalized_value = normalized_value.replace("-", " ")

            parts = [part for part in normalized_value.split() if part]
            if len(parts) != 2:
                await self.warning_function(
                    f"Invalid format for {field_name}. Use two values such as '1.5 2.3' or '1.5-2.3'."
                )
                return None

            try:
                min_r, max_r = map(float, parts)
            except ValueError as exc:
                await self.warning_function(f"Invalid format for {field_name}: {exc}")
                return None

            if min_r > max_r:
                await self.warning_function("The interval range must be provided in ascending order.")
                return None

            return min_r, max_r

        # Validate inputs and set attributes
        self.r_interval = await read_interval(self.textInput_r_interval, "r interval range")
        if self.r_interval is None:
            return
        self.min_r, self.max_r = self.r_interval
        self.atom_labels = await read_input(self.textInput_atom_labels, "atom labels", list)
        if self.atom_labels is None:
            return

    async def select_trajec(self, widget):
        """Let the user select a trajectory file and set related attributes."""
        trajec_file = await self.main_window.dialog(
            toga.OpenFileDialog(title="Select TRAJEC.xyz file")
        )
        if trajec_file:
            self.trajec = trajec_file
            self.textInput_file.value = trajec_file
            self.multi_line_text.value = f"Selected file: {self.trajec}\n"
            self.output_dir = os.path.dirname(trajec_file)
            # File will be opened only during bond_length to avoid leaving open files
            # self.fout = open(f"{self.output_dir}/min_free_energ_frames.xyz", "w")

    async def bond_length(self, widget):
        """
        Calculates the bond length between two given atom indices for each frame and
        saves qualifying frames to an output file (min_free_energ_frames.xyz).
        """
        idx1, idx2 = map(int, self.atom_labels)
        output_path = os.path.join(self.output_dir, "selected_distance_range_frames.xyz")
        average_frame_path = os.path.join(self.output_dir, "frame_closest_to_average_distance.xyz")
        summary_path = os.path.join(self.output_dir, "selected_distance_range_summary.txt")
        frame_number = 0
        selected_frame = 0
        selected_distances = []
        selected_frame_records = []
        file_size = max(os.path.getsize(self.trajec), 1)
        try:
            with open(self.trajec, "r") as f, open(output_path, "w") as fout:
                self._set_progress(0)
                while True:
                    title_line = f.readline()
                    if not title_line:
                        break
                    self.num_atoms = int(title_line.strip())
                    if self.num_atoms <= 0:
                        raise ValueError("Invalid number of atoms in the file")
                    comment_line = f.readline()
                    if not comment_line:
                        break

                    atom_data = []
                    for _ in range(self.num_atoms):
                        line = f.readline()
                        if not line:
                            raise ValueError("Unexpected end of file")
                        tokens = line.strip().split()
                        if len(tokens) < 4:
                            raise ValueError("Invalid atom line format")
                        element = tokens[0]
                        x, y, z = map(float, tokens[1:4])
                        atom_data.append((element, x, y, z))

                    # Extract positions
                    try:
                        elmt1, x1, y1, z1 = atom_data[idx1 - 1]
                        elmt2, x2, y2, z2 = atom_data[idx2 - 1]
                    except IndexError:
                        await self.warning_function(
                            f"Atom index out of range in frame {frame_number+1}."
                        )
                        return

                    bond_length = np.linalg.norm([x2 - x1, y2 - y1, z2 - z1])

                    # Write this frame if bond length is within the specified range
                    if self.min_r <= bond_length <= self.max_r:
                        selected_frame += 1
                        selected_distances.append(float(bond_length))
                        frame_content = self._format_xyz_frame(
                            self.num_atoms,
                            frame_number,
                            atom_data,
                            comment_suffix=f"distance = {bond_length:.6f} A",
                        )
                        selected_frame_records.append(
                            {
                                "frame_number": frame_number,
                                "distance": float(bond_length),
                                "content": frame_content,
                            }
                        )
                        fout.write(frame_content)

                    frame_number += 1
                    if frame_number % 100 == 0:
                        self._set_progress((f.tell() / file_size) * 100.0)
                        await asyncio.sleep(0)
                self._set_progress(100)
        except Exception as e:
            self._set_progress(0)
            await self.warning_function(f"Error reading TRAJEC.xyz file: {e}")

        average_distance = None
        standard_deviation = None
        closest_frame_number = None

        if selected_distances:
            average_distance = float(np.mean(selected_distances))
            standard_deviation = float(np.std(selected_distances))
            closest_record = min(
                selected_frame_records,
                key=lambda record: abs(record["distance"] - average_distance),
            )
            closest_frame_number = closest_record["frame_number"]
            with open(average_frame_path, "w", encoding="utf-8") as average_frame_file:
                average_frame_file.write(closest_record["content"])
            self._write_summary_file(
                summary_path,
                average_distance,
                standard_deviation,
                frame_number,
                selected_frame,
            )
        else:
            self._write_summary_file(
                summary_path,
                None,
                None,
                frame_number,
                selected_frame,
            )

        self.multi_line_text.value = f"\n\n  Number of frames analyzed: {frame_number}\n"
        self.multi_line_text.value += f"  Number of frames extracted: {selected_frame}\n"
        self.multi_line_text.value += f"  Extracted frames file: {output_path}\n"

        if average_distance is not None and standard_deviation is not None:
            self.multi_line_text.value += f"  Average distance: {average_distance:.6f} A\n"
            self.multi_line_text.value += f"  Standard deviation: {standard_deviation:.6f} A\n"
            self.multi_line_text.value += f"  Closest-to-average frame: {closest_frame_number}\n"
            self.multi_line_text.value += f"  Closest frame file: {average_frame_path}\n"
            self.multi_line_text.value += f"  Summary file: {summary_path}\n"
        else:
            self.multi_line_text.value += "  No frames were selected, so no closest-to-average frame file was created.\n"
            self.multi_line_text.value += f"  Summary file: {summary_path}\n"

class RangeFramesSelectionUI(RangeFramesSelection):
    def __init__(self, *args):
        super().__init__()
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        # Window
        self.main_window = toga.Window(
            title="Frame selection by interatomic distance range",
            size=(700, 600),
        )
        # Styles
        heading_style = Pack(font_size=18, font_weight="bold", margin_bottom=10)
        label_style = Pack(margin=5, text_align=LEFT, width=200)
        input_style = Pack(flex=1, margin=5)
        button_style = Pack(margin=5, width=100)
        box_style = Pack(direction=ROW, align_items=CENTER, margin_bottom=5)

        main_box = toga.Box(style=Pack(direction=COLUMN, margin=20))

        # Heading and progress
        box_1 = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin_bottom=10))
        main_box.add(box_1)
        box_1a = toga.Box(style=Pack(flex=1))
        box_1.add(box_1a)

        title_label = toga.Label("Frame selection by interatomic distance range", style=heading_style)
        box_1a.add(title_label)

        # Input fields
        input_fields = [
            ("Distance Range:", "Enter range such as 1.5 2.3 or 1.5-2.3", "textInput_r_interval"),
            ("Atom Labels:", "Enter two atom labels, separated by a space (e.g., 1 2)", "textInput_atom_labels"),
        ]
        for label_text, placeholder, attr_name in input_fields:
            box = toga.Box(style=box_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            box.add(label)
            box.add(text_input)
            main_box.add(box)

        # File selection
        file_box = toga.Box(style=box_style)
        file_label = toga.Label("Select Trajectory File:", style=label_style)
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select TRAJEC.xyz file", style=input_style
        )
        browse_button = toga.Button("Browse", on_press=self.select_trajec, style=button_style)
        file_box.add(file_label)
        file_box.add(self.textInput_file)
        file_box.add(browse_button)
        main_box.add(file_box)

        # Progress bar
        progress_box = toga.Box(style=Pack(direction=COLUMN))
        self.progress_bar = toga.ProgressBar(max=100)
        progress_box.add(self.progress_bar)
        main_box.add(progress_box)

        # Help text / Output
        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin_top=10, font_size=12)
        )
        self.multi_line_text.value = HelpGqteaWin.minFframesSelection
        main_box.add(self.multi_line_text)

        # Bottom buttons
        button_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin_top=10)
        )
        self.btn_execute = toga.Button(
            "Extract Frames", style=button_style, on_press=self.workflow
        )
        self.btn_close = toga.Button(
            "Close", style=button_style, on_press=self.closeTopLevel
        )
        button_box.add(self.btn_execute)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        self.main_window.content = main_box
        self.main_window.show()

    async def workflow(self, widget):
        await self.read_params(widget)
        if not self.trajec:
            await self.warning_function("No trajectory file selected.")
            return
        await self.bond_length(widget)

    def closeTopLevel(self, widget):
        self.main_window.close()
