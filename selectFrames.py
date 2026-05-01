import asyncio
import os

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT

from help import HelpGqteaWin
from framesCounter import FramesCounter


class SelectFrames(FramesCounter):
    progress_update_interval = 200

    def _set_progress(self, percent: float) -> None:
        self.progress_bar.value = max(0.0, min(100.0, percent))

    def _set_gaussian_fields_enabled(self, enabled: bool) -> None:
        gaussian_fields = [
            self.textInput_chk_filename,
            self.textInput_amount_mem,
            self.textInput_nproc,
            self.textInput_route,
            self.textInput_charge,
            self.textInput_multiplicity,
        ]
        for field in gaussian_fields:
            field.enabled = enabled

    def toggle_gaussian_fields(self, widget) -> None:
        self._set_gaussian_fields_enabled(bool(self.switch_generate_gaussian.value))

    def _append_output(self, text: str) -> None:
        current = self.multi_line_text.value.rstrip()
        if current:
            self.multi_line_text.value = f"{current}\n{text}"
        else:
            self.multi_line_text.value = text

    async def frames_counter(self, widget):
        """Count frames in TRAJEC.xyz while updating the shared progress bar."""
        await self.open_file_dialog(widget)
        if not hasattr(self, "trajec") or not self.trajec:
            return

        self.multi_line_text.value = "\nReading TRAJEC.xyz file ........\n"
        self._set_progress(0)

        frame_count = 0
        total_lines = 0
        file_size = max(os.path.getsize(self.trajec), 1)

        try:
            with open(self.trajec, "r") as f:
                while True:
                    title_line = f.readline()
                    if not title_line:
                        break

                    comment_line = f.readline()
                    if not comment_line:
                        break
                    total_lines += 2

                    for _ in range(self.num_atoms):
                        atom_line = f.readline().split()
                        total_lines += 1
                        if not atom_line or len(atom_line) != 4:
                            raise ValueError(
                                f"Error reading atom line {total_lines} in frame {frame_count}"
                            )

                    frame_count += 1
                    if frame_count % self.progress_update_interval == 0:
                        self._set_progress((f.tell() / file_size) * 100.0)
                        await asyncio.sleep(0)

                self._set_progress(100)

            number_of_atoms = f"Number of atoms  -->  {self.num_atoms}\n"
            number_of_lines = f"Number of lines  -->  {total_lines}\n"
            number_of_frames = f"Number of frames -->  {frame_count}\n"
            first_frame = "  TRAJEC.xyz FIRST FRAME\n"

            with open(self.trajec, "r") as f:
                frame_lines = [f.readline() for _ in range(self.num_atoms + 2)]

            update_text = (
                number_of_atoms
                + number_of_lines
                + number_of_frames
                + first_frame
                + "".join(frame_lines)
            )
            self.multi_line_text.value = update_text
            self.total_frame_number = frame_count
            self.last_frame_index = max(0, frame_count - 1)
            self.textInput_stop_frame.value = str(self.last_frame_index)
        except Exception as e:
            self._set_progress(0)
            await self.main_window.dialog(
                toga.InfoDialog("Error", f"Error reading TRAJEC.xyz file: {e}")
            )

    async def read_params(self, widget) -> bool:
        """Read the input parameters."""

        async def read_input(text_input, field_name, expected_type):
            value = text_input.value.strip()
            if not value:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Please input a valid value for {field_name}.")
                )
                return None
            try:
                return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        self.start_frame = await read_input(self.textInput_start_frame, "starting frame", int)
        if self.start_frame is None:
            return False

        self.skipped_frames = await read_input(
            self.textInput_skipped_frames, "number of frame to be skipped", int
        )
        if self.skipped_frames is None:
            return False

        self.stop_frame = await read_input(
            self.textInput_stop_frame, "stop frame to be collected", int
        )
        if self.stop_frame is None:
            return False

        if self.start_frame < 0:
            await self.main_window.dialog(
                toga.InfoDialog("Error", "Starting frame must be zero or greater.")
            )
            return False

        if self.skipped_frames < 0:
            await self.main_window.dialog(
                toga.InfoDialog("Error", "Number of skipped frames must be zero or greater.")
            )
            return False

        if self.stop_frame < self.start_frame:
            await self.main_window.dialog(
                toga.InfoDialog("Error", "Stop frame must be greater than or equal to the starting frame.")
            )
            return False

        if getattr(self, "total_frame_number", 0) and self.stop_frame > getattr(self, "last_frame_index", self.total_frame_number - 1):
            await self.main_window.dialog(
                toga.InfoDialog(
                    "Error",
                    f"Stop frame must be smaller than or equal to {getattr(self, 'last_frame_index', self.total_frame_number - 1)}.",
                )
            )
            return False

        self.chk_filename = ""
        self.amount_mem = None
        self.nproc = None
        self.route = ""
        self.charge = None
        self.multiplicity = None

        if self.switch_generate_gaussian.value:
            self.chk_filename = await read_input(
                self.textInput_chk_filename, "gaussian chk file name", str
            )
            if self.chk_filename is None:
                return False

            self.amount_mem = await read_input(
                self.textInput_amount_mem,
                "amount of gaussian memory allocation in GB",
                int,
            )
            if self.amount_mem is None:
                return False

            self.nproc = await read_input(self.textInput_nproc, "number of process", int)
            if self.nproc is None:
                return False

            self.route = await read_input(self.textInput_route, "gaussian section route", str)
            if self.route is None:
                return False

            self.charge = await read_input(
                self.textInput_charge, "total charge on the system", int
            )
            if self.charge is None:
                return False

            self.multiplicity = await read_input(self.textInput_multiplicity, "multiplicity", int)
            if self.multiplicity is None:
                return False

        update_lines = [
            f"{'Start frame':.<30} {self.start_frame:>20}",
            f"{'Number of frames skipped':.<30}{self.skipped_frames:>20}",
            f"{'Stop frame to be collected':.<30}{self.stop_frame:>20}",
            f"{'Generate Gaussian inputs':.<30}{str(bool(self.switch_generate_gaussian.value)):>20}",
        ]

        if self.switch_generate_gaussian.value:
            update_lines.extend(
                [
                    f"{'Gaussian chk file name':.<30}{self.chk_filename:>20}",
                    f"{'Amount of memory':.<30}{self.amount_mem:>20}",
                    f"{'Number of process':.<30}{self.nproc:>20}",
                    f"{'Gaussian section route':.<30}{self.route:>20}",
                    f"{'Total charge on the system':.<30}{self.charge:>20}",
                    f"{'Multiplicity':.<30}{self.multiplicity:>20}",
                ]
            )

        self.multi_line_text.value = "\n".join(update_lines) + "\n"
        return True

    async def select_frames(self) -> bool:
        try:
            start_frame = self.start_frame
            skip_frames = self.skipped_frames
            stop_frame = self.stop_frame
            step = skip_frames + 1

            self._set_progress(0)
            self._append_output("Selecting frames...")

            with open(self.trajec, "r") as f, open(
                f"{self.output_dir}/selected_frames.xyz", "w"
            ) as f2:
                frame = 0
                selected_count = 0
                while frame <= stop_frame:
                    should_collect = frame >= start_frame and ((frame - start_frame) % step == 0)

                    if should_collect:
                        for _ in range(self.num_atoms + 2):
                            line = f.readline()
                            if not line:
                                raise ValueError("Unexpected end of file while selecting frames.")
                            f2.write(line)
                        selected_count += 1
                    else:
                        for _ in range(self.num_atoms + 2):
                            line = f.readline()
                            if not line:
                                raise ValueError("Unexpected end of file while skipping frames.")

                    frame += 1
                    if stop_frame > 0 and (
                        frame % self.progress_update_interval == 0 or frame > stop_frame
                    ):
                        self._set_progress((frame / (stop_frame + 1)) * 100.0)
                        await asyncio.sleep(0)

            self._set_progress(100)
            self._append_output(f"Frames have been successfully selected: {selected_count}")
            return True
        except Exception as e:
            self._set_progress(0)
            self._append_output(f"An error occurred: {str(e)}")
            return False

    def gaussian_input(self):
        chk = self.chk_filename
        mem = self.amount_mem
        nproc = self.nproc
        key_words = self.route
        charge = self.charge
        multiplicity = self.multiplicity

        try:
            with open(f"{self.output_dir}/selected_frames.xyz", "r") as f:
                n = 1
                while True:
                    line = f.readline()
                    line = f.readline()
                    if line == "\n" or not line:
                        break

                    with open(f"{self.output_dir}/g16_input_{n}.gjf", "w") as g16_inputs:
                        g16_inputs.write(f"%chk={chk}\n%Mem={mem}\n%nproc={nproc}\n# {key_words}\n\n")
                        g16_inputs.write("Gaussian input created by gqteaWin\n\n")
                        g16_inputs.write(f"{charge} {multiplicity}\n")

                        for _ in range(self.num_atoms):
                            line = f.readline()
                            g16_inputs.write(line)

                    n += 1

            self._append_output(f"Gaussian inputs successfully created: {n - 1}")
        except Exception as e:
            self._append_output(f"An error occurred: {str(e)}")


class SelectFramesUI(SelectFrames):
    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        del widget
        self.main_window = toga.Window(
            title="CPMD Input Files",
            size=(760, 680),
        )

        heading_style = Pack(font_size=18, font_weight="bold", text_align=LEFT, margin=(0, 0, 10, 0))
        label_style = Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=230)
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=110)
        row_style = Pack(direction=ROW, margin=(0, 0, 5, 0), align_items="center")

        main_box = toga.Box(style=Pack(direction=COLUMN, margin=20))

        box_title = toga.Box(style=Pack(direction=ROW, margin=(0, 0, 10, 0)))
        sub_box_title = toga.Box(style=Pack(flex=1))
        sub_box_progress = toga.Box(style=Pack(width=220))
        box_title.add(sub_box_title)
        box_title.add(sub_box_progress)

        title_label = toga.Label("Select Frames From TRAJEC.xyz file", style=heading_style)
        self.progress_bar = toga.ProgressBar(max=100, style=Pack(width=220))
        sub_box_title.add(title_label)
        sub_box_progress.add(self.progress_bar)
        main_box.add(box_title)

        input_fields = [
            ("Starting frame:", "Starting frame to start extraction", "textInput_start_frame", ""),
            (
                "Number of frames to be skipped",
                "Enter number of frames to be skipped between collected frames",
                "textInput_skipped_frames",
                "0",
            ),
            ("Stop frame to be collected:", "Enter the stop frame to be collected", "textInput_stop_frame", ""),
        ]

        for label_text, placeholder, attr_name, default_value in input_fields:
            box = toga.Box(style=row_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(
                value=default_value,
                placeholder=placeholder,
                style=input_style,
            )
            setattr(self, attr_name, text_input)
            box.add(label)
            box.add(text_input)
            main_box.add(box)

        switch_box = toga.Box(style=Pack(direction=COLUMN, margin=(0, 0, 8, 0)))
        self.switch_generate_gaussian = toga.Switch(
            "Check box for gaussian input file generation",
            value=False,
            on_change=self.toggle_gaussian_fields,
            style=Pack(text_align=LEFT),
        )
        switch_box.add(self.switch_generate_gaussian)
        main_box.add(switch_box)

        gaussian_fields = [
            ("Gaussian chk file name:", "Enter the gaussian chk file name", "textInput_chk_filename"),
            (
                "Gaussian amount memory:",
                "Enter the amount memory, e.g. Mem=8GB",
                "textInput_amount_mem",
            ),
            ("Number of process:", "Enter the number of process, e.g., nproc=10", "textInput_nproc"),
            ("Gaussian section route:", "Enter the gaussian section route", "textInput_route"),
        ]

        for label_text, placeholder, attr_name in gaussian_fields:
            box = toga.Box(style=row_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            box.add(label)
            box.add(text_input)
            main_box.add(box)

        charge_multi_box = toga.Box(style=row_style)
        charge_multi_label = toga.Label("Total charge and multiplicity:", style=label_style)
        self.textInput_charge = toga.TextInput(
            placeholder="Enter the charge on the system",
            style=Pack(flex=1, margin=(5, 5)),
        )
        self.textInput_multiplicity = toga.TextInput(
            placeholder="Enter the multiplicity",
            style=Pack(flex=1, margin=(5, 5)),
        )
        charge_multi_box.add(charge_multi_label)
        charge_multi_box.add(self.textInput_charge)
        charge_multi_box.add(self.textInput_multiplicity)
        main_box.add(charge_multi_box)

        file_box = toga.Box(style=row_style)
        file_label = toga.Label("Select input file:", style=label_style)
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select TRAJEC.xyz file",
            style=input_style,
        )
        browse_button = toga.Button("Browse", on_press=self.frames_counter, style=button_style)
        file_box.add(file_label)
        file_box.add(self.textInput_file)
        file_box.add(browse_button)
        main_box.add(file_box)

        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(5, 0), font_size=11)
        )
        self.multi_line_text.value = "\nSelect frames from TRAJEC.xyz file and build gaussian input files"
        main_box.add(self.multi_line_text)

        button_box = toga.Box(style=Pack(direction=ROW, margin_top=5, align_items="center"))
        self.btn_execute = toga.Button("Select Frames", style=button_style, on_press=self.workflow)
        self.btn_help = toga.Button("Help", style=button_style, on_press=self.open_window_help)
        self.btn_close = toga.Button("Close", style=button_style, on_press=self.closeTopLevel)
        button_box.add(self.btn_execute)
        button_box.add(self.btn_help)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        self.main_window.content = main_box
        self.main_window.show()
        self.toggle_gaussian_fields(None)

    def open_window_help(self, widget):
        del widget
        window = toga.Window(title="Instruction to use selection frames module ")
        main_box = toga.Box(style=Pack(direction=COLUMN, flex=1))
        multi_line_text = toga.MultilineTextInput(
            style=Pack(font_size=11, margin=5, flex=1)
        )
        multi_line_text.value = HelpGqteaWin.help_selectFrames
        main_box.add(multi_line_text)
        window.content = main_box
        window.show()

    async def workflow(self, widget):
        if not await self.read_params(widget):
            return

        if not hasattr(self, "trajec") or not self.trajec:
            await self.main_window.dialog(
                toga.InfoDialog("Error", "No trajectory file selected.")
            )
            return

        if not await self.select_frames():
            return

        if self.switch_generate_gaussian.value:
            self.gaussian_input()

    def closeTopLevel(self, widget):
        self.main_window.close()
