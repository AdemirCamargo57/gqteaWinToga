import os, tempfile, toga
import matplotlib.pyplot as plt
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from help import HelpGqteaWin


class SHGeomAnalyzer:
    """Non-GUI logic for SH trajectory analysis."""

    def __init__(self):
        self.nStates = None
        self.root_path = None

    async def warning_function(self, msg):
        await self.main_window.dialog(
            toga.InfoDialog("Error", str(msg))
        )

    async def read_params(self):
        """Reads and validates parameters from UI text inputs."""
        try:
            self.nStates = int(self.textInput_nStates.value.strip())
        except (ValueError, AttributeError):
            await self.warning_function("Please input a valid integer for number of states.")
            return False
        return True

    async def select_dir(self, widget):
        self.root_path = await self.main_window.select_folder_dialog(
            title="Select SH root directory"
        )
        if self.root_path:
            self.textInput_file.value = self.root_path
            self.multi_line_text.value = f"\nSelected directory: {self.root_path}\n"

    async def check_sh_state_integrity(self, directory):
        file_path = os.path.join(directory, "SH_STATE.dat")
        errors = []

        try:
            with open(file_path, "r") as f:
                for line_num, line in enumerate(f, start=1):
                    tokens = line.strip().split()
                    if len(tokens) != 2:
                        errors.append(f"{file_path} --> Line {line_num}: Expected 2 items, found {len(tokens)}.")
                        continue
                    try:
                        _ = int(tokens[0])
                        state = int(tokens[1])
                    except ValueError:
                        errors.append(f"{file_path} --> Line {line_num}: Non-integer values detected.")
                        continue

                    if not (0 <= state <= self.nStates):
                        errors.append(f"{file_path} --> Line {line_num}: Invalid state {state}. Must be >=0 and <={self.nStates}.")
        except FileNotFoundError:
            errors.append(f"{file_path} not found.")

        if errors:
            await self.warning_function("\n".join(errors))
            return False
        return True

    async def check_trajec_integrity(self, directory):
        """Check TRAJEC.xyz integrity with detailed line numbers."""
        file_path = os.path.join(directory, "TRAJEC.xyz")
        errors = []

        try:
            with open(file_path, "r") as f:
                frame_num = 0
                current_line_number = 0

                while True:
                    num_atoms_line = f.readline()
                    current_line_number += 1
                    if not num_atoms_line:
                        break  # EOF

                    frame_num += 1

                    # Parse number of atoms
                    try:
                        num_atoms = int(num_atoms_line.strip())
                    except ValueError:
                        errors.append(
                            f"{file_path} --> Line {current_line_number}: Invalid number of atoms '{num_atoms_line.strip()}'"
                        )
                        break

                    # Read comment line
                    comment_line = f.readline()
                    current_line_number += 1
                    if not comment_line:
                        errors.append(
                            f"{file_path} --> Line {current_line_number}: Missing comment line."
                        )
                        break

                    # Read atom lines
                    for atom_idx in range(1, num_atoms + 1):
                        atom_line = f.readline()
                        current_line_number += 1

                        if not atom_line:
                            errors.append(
                                f"{file_path} --> Line {current_line_number}: Incomplete atom lines. Expected {num_atoms} atoms, found {atom_idx - 1}."
                            )
                            break

                        tokens = atom_line.strip().split()

                        if len(tokens) != 4:
                            errors.append(
                                f"{file_path} --> Line {current_line_number}: Expected 4 items (label + 3 coordinates), found {len(tokens)}."
                            )
                            continue

                        try:
                            _ = str(tokens[0])
                            x, y, z = map(float, tokens[1:])
                        except ValueError:
                            errors.append(
                                f"{file_path} --> Line {current_line_number}: Non-numeric coordinate detected."
                            )

        except FileNotFoundError:
            errors.append(f"{file_path} not found.")

        if errors:
            await self.warning_function("\n".join(errors))
            return False
        return True


    async def analyze(self):
        if not self.root_path or self.nStates is None:
            await self.warning_function("Please select a directory and input the number of states before running.")
            return

        directories = [os.path.join(self.root_path, d) for d in os.listdir(self.root_path) if os.path.isdir(os.path.join(self.root_path, d))]
        self.multi_line_text.value = '\nStarting analysis...\n\n'

        avg_state = {f"state{state}": 1 for state in range(self.nStates + 1)}
        nRunSH = 0

        for directory in directories:
            self.multi_line_text.value = f"\nProcessing directory: {directory}\n"

            # Integrity checks
            valid_state = await self.check_sh_state_integrity(directory)
            valid_trajec = await self.check_trajec_integrity(directory)

            if not (valid_state and valid_trajec):
                self.multi_line_text.value = f" Integrity check failed in {directory}. Skipping.\n"
                nRunSH += 1
                continue
            else:
                self.multi_line_text.value = f" Integrity check passed in {directory}\n"

            # Parse SH_STATE.dat
            state_per_frame = []
            with open(os.path.join(directory, "SH_STATE.dat"), "r") as f:
                for line in f:
                    tokens = line.strip().split()
                    state_per_frame.append(int(tokens[1]))

            # Count states
            state_counts = {}
            for state in state_per_frame:
                key = f"state{state}"
                state_counts[key] = state_counts.get(key, 0) + 1

            for key, value in state_counts.items():
                avg_state[key] += value

            # Process trajectory
            with open(os.path.join(directory, "TRAJEC.xyz"), 'r') as traj_file:
                frame_idx = 0
                while True:
                    num_atoms_line = traj_file.readline()
                    if not num_atoms_line:
                        break

                    comment_line = traj_file.readline()
                    natoms = int(num_atoms_line.strip())
                    atom_lines = [traj_file.readline() for _ in range(natoms)]
                    if not all(atom_lines):
                        self.multi_line_text.value = f"\nIncomplete frame in {directory}\n"
                        break

                    if frame_idx < len(state_per_frame):
                        state = state_per_frame[frame_idx]
                    else:
                        state = "unknown"

                    out_path = os.path.join(directory, f"state{state}.xyz")
                    with open(out_path, 'a') as out:
                        out.write(num_atoms_line)
                        out.write(comment_line)
                        out.writelines(atom_lines)

                    frame_idx += 1

        # Combine state files
        state_files = [f"state{n}.xyz" for n in range(1, self.nStates + 1)]
        for state in state_files:
            all_contents = []
            for subdir, dirs, files in os.walk(self.root_path):
                for filename in files:
                    if filename == state:
                        file_path = os.path.join(subdir, filename)
                        with open(file_path, 'r') as file:
                            all_contents.append(file.read())
            with open(os.path.join(self.root_path, state), 'w') as outfile:
                for item in all_contents:
                    outfile.write(item)

        # Remove blank lines
        for state in state_files:
            file_path = os.path.join(self.root_path, state)
            with open(file_path, 'r') as f:
                lines = f.readlines()
            non_blank_lines = [line for line in lines if line.strip()]
            with open(file_path, 'w') as f:
                f.writelines(non_blank_lines)

        # Calculate averages
        nRunSH_eff = len(directories) - nRunSH if (len(directories) - nRunSH) > 0 else 1
        for key in avg_state:
            avg_state[key] = avg_state[key] / nRunSH_eff

        # Write average state occupancy (excited states only)
        output_file = os.path.join(self.root_path, "sh_avg_perc.dat")
        try:
            with open(output_file, "w") as f:
                total = sum(avg_state.get(f"state{n}", 0) for n in range(1, self.nStates + 1))
                for n in range(1, self.nStates + 1):
                    key = f"state{n}"
                    avg = avg_state.get(key, 0)
                    prct = (avg * 100.0) / total if total else 0
                    f.write(f"{n:>4}    {avg:>10.3f} {prct:>10.3f}\n")
        except IOError as e:
            await self.warning_function(f"An error occurred while writing to {output_file}: {e}")
            return
        self.show_plot()

    def show_plot(self):
        """Display a plot in a new window."""
        avg_bar_file = tempfile.NamedTemporaryFile(delete=False, suffix=".png", dir=self.root_path).name
        x = []
        y = []
        file_path = os.path.join(self.root_path, "sh_avg_perc.dat")
        if not os.path.exists(file_path):
            self.multi_line_text.value = "\nsh_avg_perc.dat file not found.\n"
            return

        with open(file_path, 'r') as f:
            for linha in f:
                tokens = linha.strip().split()
                if len(tokens) >= 3:
                    x.append(float(tokens[0]))
                    y.append(float(tokens[2]))

        if not x or not y:
            self.multi_line_text.value = "\nNo data to plot.\n"
            return

        # Bar plot is more appropriate for state occupancy
        plt.figure(figsize=(8, 5))
        bars = plt.bar(x, y,color='red')
        plt.title("Average State Occupancy")
        plt.xlabel("State Number")
        plt.ylabel("Average Occupancy (%)")
        # Add value labels to each bar
        for bar in bars:
            height = bar.get_height()
            plt.annotate(f'{height:.1f}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom')
        plt.tight_layout()
        plt.savefig(avg_bar_file)
        plt.close()  # Correct usage

        # Load the image using Toga's Image class
        plot_image = toga.Image(avg_bar_file)

        # Create a new window for the plot
        plot_window = toga.Window(
            title="Average State Occupancy",
            size=(700, 600),
        )

        # Create a box to hold the image
        plot_box = toga.Box(style=Pack(flex=1))
        plot_window.content = plot_box

        # Create an ImageView to display the image
        plot_imageview = toga.ImageView(plot_image, style=Pack(flex=1))
        plot_box.add(plot_imageview)

        # Show the window
        plot_window.show()
        

class SHGeomAnalyzerUI(SHGeomAnalyzer):
    """UI for extracting and grouping trajectory frames based on states."""

    def __init__(self, *args):
        super().__init__()
        self.mainWindow()

    def mainWindow(self):
        self.main_window = toga.Window(
            title="Extract and group trajectory frames based on states",
            size=(700, 600)
        )

        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style = Pack(margin=(5, 5), text_align=LEFT, width=200)
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=100)
        box_style = Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0))

        main_box = toga.Box(style=Pack(direction=COLUMN, margin=20))

        # Title
        box_1 = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 10, 0)))
        main_box.add(box_1)
        box_1a, box_1b, box_1c = (toga.Box(style=Pack(width=w)) for w in (400, 150, 100))
        box_1.add(box_1a)
        box_1.add(box_1b)
        box_1.add(box_1c)
        title_label = toga.Label("Separate Frames by State", style=heading_style)
        empty_label = toga.Label(" ", style=heading_style)
        self.progress_label = toga.Label(" ", style=heading_style)
        box_1a.add(title_label)
        box_1b.add(empty_label)
        box_1c.add(self.progress_label)

        # Input fields
        input_fields = [
            ("Number of States:", "Enter the total number of states simulated", "textInput_nStates"),
        ]

        for label_text, placeholder, attr_name in input_fields:
            box = toga.Box(style=box_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            box.add(label)
            box.add(text_input)
            main_box.add(box)

        # Directory selection
        file_box = toga.Box(style=box_style)
        file_label = toga.Label("Select SH Root Directory:", style=label_style)
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select SH root directory", style=input_style
        )
        browse_button = toga.Button("Browse", on_press=self.select_dir, style=button_style)
        file_box.add(file_label)
        file_box.add(self.textInput_file)
        file_box.add(browse_button)
        main_box.add(file_box)

        # Output / help text
        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(10, 0), font_size=12)
        )
        self.multi_line_text.value = HelpGqteaWin.sh_geom_analyzer
        main_box.add(self.multi_line_text)

        # Bottom buttons
        button_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin_top=10))
        self.btn_execute = toga.Button(
            "Extract Frames", style=button_style, on_press=self.on_extract_frames
        )
        self.btn_close = toga.Button(
            "Close", style=button_style, on_press=self.close_main_window
        )
        button_box.add(self.btn_execute)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        self.main_window.content = main_box
        self.main_window.show()

    async def on_extract_frames(self, widget):
        valid = await self.read_params()
        if not valid:
            return
        await self.analyze()

    def close_main_window(self, widget):
        self.main_window.close()
