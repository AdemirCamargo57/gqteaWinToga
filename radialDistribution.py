import toga, tempfile
import numpy as np
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from help import HelpGqteaWin
from framesCounter import FramesCounter
from matplotlib import pyplot as plt


class RadialAnalyser(FramesCounter):
    def read_params(self, widget):
        """Read and validate the input parameters."""

        def read_input(text_input, field_name, expected_type):
            value = text_input.value.strip()
            if not value:
                self.main_window.error_dialog(
                    "Error", f"Please input a valid value for {field_name}."
                )
                return None
            try:
                if expected_type == "int_list":
                    labels = [int(label) for label in value.split()]
                    return labels
                elif expected_type == "float_list":
                    labels = [float(label) for label in value.split()]
                    return labels
                else:
                    return expected_type(value)
            except ValueError as e:
                self.main_window.error_dialog(
                    "Error", f"Invalid format for {field_name}: {e}"
                )
                return None

        # Read and validate inputs
        self.radius = read_input(
            self.textInput_radius, "radius for radial distribution function calculation", float
        )
        if self.radius is None:
            return

        self.bin_width = read_input(
            self.textInput_bin_width, "bin width for histogram distribution function", float
        )
        if self.bin_width is None:
            return

        self.delete_atom_list = read_input(
            self.textInput_atom_list, "list of atom labels to be excluded", "int_list"
        )
        if self.delete_atom_list is None:
            self.delete_atom_list = [0]
            
        if self.delete_atom_list == [0]:
            self.delete_atom_list = []

        self.shell_center = read_input(
            self.textInput_shell_center, "atom label at the center of the shell", int
        )
        if self.shell_center is None:
            return

        self.rdf_atom_symbol = read_input(
            self.textInput_atom_symbol, "atomic symbol for g(r) calculation", str
        )
        if self.rdf_atom_symbol is None:
            return

        self.cell_lattices = read_input(
            self.textInput_cell_lattices, "cell lattices a, b, and c", "float_list"
        )
        if self.cell_lattices is None:
            return        
        self.axis = read_input(
            self.textInput_axis, "x-axis for coordination number", "float_list")
        if self.axis is None:
            return

        # Display parameters summary
        update_text = (
            f"{'Maximum radius for RDF':.<30} {self.radius:>20}\n"
            f"{'Bin width for histogram distribution':.<30}{self.bin_width:>20}\n"
            f"{'Atom labels to be excluded':.<30}{str(self.delete_atom_list):>20}\n"
            f"{'Atom label at the center of the shell':.<30}{self.shell_center:>20}\n"
            f"{'Atomic symbol for g(r) calculation':.<30}{self.rdf_atom_symbol:>20}\n"
            f"{'Cell lattices a, b, and c':.<30}{str(self.cell_lattices):>20}\n"
        )
        self.multi_line_text.value = update_text

        self.ideal_density()

    def ideal_density(self):
        # Calculate volume and density based on first frame
        count_atom = 0
        with open(self.trajec, "r") as f:
            num_atoms = int(f.readline().strip())
            f.readline()
            for _ in range(num_atoms):
                parts = f.readline().split()
                if parts[0] == self.rdf_atom_symbol:
                    count_atom += 1
        a, b, c = self.cell_lattices
        self.volume = a * b * c
        self.rho = count_atom / self.volume


    def calcRDF(self):
        """
        Calculates the Radial Distribution Function using NumPy for performance.

        This version explicitly generates 27 periodic images of the system
        as requested, but uses vectorized operations to make the process
        faster than the original pure Python implementation.
        """
        # Setup parameters and data structures
        symbol = self.rdf_atom_symbol
        center_idx = self.shell_center - 1  # Use 0-based indexing
        bin_width = self.bin_width
        radius = self.radius
        num_bins = int(radius / bin_width)

        bins = np.linspace(0, radius, num_bins + 1)
        r = bins[:-1] # Radius values are the left edge of the bins

        total_histo = np.zeros(num_bins, dtype=np.float64)
        frame_count = 0

        # Define the 27 translation vectors for periodic images
        a, b, c = self.cell_lattices
        operations = [
            (dx, dy, dz)
            for dx in [0, a, -a]
            for dy in [0, b, -b]
            for dz in [0, c, -c]
        ]
        translation_vectors = np.array(operations)

        self.progress_bar.start()

        # 3. Process trajectory frame by frame
        with open(self.trajec, "r") as traj_file:
            while True:
                line_atoms = traj_file.readline()
                if not line_atoms:
                    break # End of file

                try:
                    num_atoms = int(line_atoms.strip())
                    traj_file.readline() # Skip comment line

                    frame_data = [traj_file.readline().split() for _ in range(num_atoms)]
                    atom_data = np.array(
                        [(parts[0], float(parts[1]), float(parts[2]), float(parts[3])) for parts in frame_data],
                        dtype=[('symbol', 'U10'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8')]
                    )
                except (ValueError, IndexError):
                    break # Stop on malformed frame

                # 4. Prepare coordinates for augmentation
                center_coords = np.array([atom_data['x'][center_idx], atom_data['y'][center_idx], atom_data['z'][center_idx]])

                # Create a mask to exclude specified atoms
                include_mask = np.ones(num_atoms, dtype=bool)
                if self.delete_atom_list:
                    delete_indices = np.array(self.delete_atom_list) - 1
                    include_mask[delete_indices] = False
                
                # The center atom should not be part of the augmented set for distance calculation
                include_mask[center_idx] = False

                base_coords = np.vstack([atom_data['x'][include_mask], atom_data['y'][include_mask], atom_data['z'][include_mask]]).T
                base_symbols = atom_data['symbol'][include_mask]
                
                if base_coords.shape[0] == 0:
                    frame_count += 1
                    continue

                # 5. Generate 27 periodic images using NumPy broadcasting
                # This is faster than a Python loop but uses significant memory.
                # It creates a new array of shape (n_atoms_included, 27, 3)
                # and then reshapes it to (n_atoms_included * 27, 3).
                augmented_coords = (base_coords[:, np.newaxis, :] + translation_vectors[np.newaxis, :, :]).reshape(-1, 3)
                
                # Repeat the symbol array to match the augmented coordinates
                augmented_symbols = np.repeat(base_symbols, 27)

                # 6. Filter for the target symbol and calculate distances
                symbol_mask = (augmented_symbols == symbol)
                final_target_coords = augmented_coords[symbol_mask]
                
                if final_target_coords.shape[0] == 0:
                    frame_count += 1
                    continue
                    
                # Calculate distances from the center atom to all atoms in the augmented set
                distances = np.sqrt(np.sum((final_target_coords - center_coords)**2, axis=1))

                # Filter distances to be within the desired radius
                distances_in_range = distances[distances <= radius]
                
                # 7. Update histogram using NumPy's optimized function
                histo, _ = np.histogram(distances_in_range, bins=bins)
                total_histo += histo

                frame_count += 1
                if frame_count % 100 == 0:
                    self.progress_bar.value = (frame_count / self.total_frame_number) * 100

        self.progress_bar.value = 100
        self.progress_bar.stop()
        if frame_count == 0:
            self.main_window.error_dialog("Error", "No frames were processed.")
            return

        # 8. Normalize and compute g(r)
        avg_histo = total_histo / frame_count
        shell_volumes = (4.0/3.0) * np.pi * (bins[1:]**3 - bins[:-1]**3)
        shell_volumes[shell_volumes < 1e-9] = 1.0 # Avoid division by zero
        
        real_density = avg_histo / shell_volumes
        g_r = real_density / self.rho
        coordination_num = np.cumsum(avg_histo)

        # 9. Save data and plot (This part remains the same)
        with open(f"{self.output_dir}/histo.dat", "w") as tmpFile:
            tmpFile.write("      r              g(r)             Integral\n\n")
            for i in range(num_bins):
                tmpFile.write(f"{r[i]:>10.5f}{g_r[i]:>20.7f}{coordination_num[i]:>20.7f}\n")

        # Plot g(r)
        rdf_file = tempfile.NamedTemporaryFile(delete=False, suffix="_rdf.png", dir=self.output_dir).name
        fig1, ax1 = plt.subplots()
        ax1.plot(r, g_r)
        ax1.set_title("Radial Distribution Function g(r)")
        ax1.set_xlabel("r (Å)")
        center_symbol = atom_data['symbol'][center_idx] if 'atom_data' in locals() else self.rdf_atom_symbol
        ax1.set_ylabel(f"g(r) ({center_symbol}{self.shell_center}---{symbol})")
        plt.savefig(rdf_file)
        plt.close(fig1)
        self.show_plot(rdf_file)

        # Plot coordination number
        coord_file = tempfile.NamedTemporaryFile(delete=False, suffix="_coordination.png", dir=self.output_dir).name
        fig2, ax2 = plt.subplots()
        ax2.plot(r, coordination_num)
        ax2.set_title("Coordination Number")
        ax2.set_xlabel("r (Å)")
        ax2.set_ylabel("Coordination number")
        ax2.set_xlim(0, self.axis[0])
        ax2.set_ylim(0, self.axis[1])
        plt.savefig(coord_file)
        plt.close(fig2)
        self.show_plot(coord_file)

    def show_plot(self, temp_filename):
        plot_image = toga.Image(temp_filename)
        plot_window = toga.Window(
            title="Plot",
            size=(900, 600),
        )
        plot_box = toga.Box(style=Pack(direction=COLUMN, flex=1))
        plot_window.content = plot_box
        plot_imageview = toga.ImageView(plot_image, style=Pack(flex=1, padding=10))
        plot_box.add(plot_imageview)
        plot_window.show()


class RadialFunctionUI(RadialAnalyser):
    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        # Create the main window
        self.main_window = toga.Window(
            title="Radial Distribution Function",
            size=(700, 600),
        )

        # Common styles
        heading_style = Pack(font_size=18, font_weight="bold", text_align=LEFT, padding=(0, 0, 10, 0))
        label_style = Pack(padding=(0, 0, 5, 5), text_align=LEFT, width=200)
        input_style = Pack(flex=1, padding=(5, 5))
        button_style = Pack(padding=5, width=100)
        box_style = Pack(direction=ROW, alignment=CENTER, padding=(0, 0, 5, 0))

        main_box = toga.Box(style=Pack(direction=COLUMN, padding=20))

        # Title and progress
        box_title = toga.Box(style=Pack(direction=ROW, alignment=CENTER, padding=(0, 0, 10, 0)))
        sub_title = toga.Box(style=Pack(width=410))
        sub_empty = toga.Box(style=Pack(width=140))
        sub_progress = toga.Box(style=Pack(width=100))
        title_lbl = toga.Label("Radial Distribution Function g(r)", style=heading_style)
        empty_lbl = toga.Label(" ", style=heading_style)
        self.progress_label = toga.Label(" ", style=heading_style)
        sub_title.add(title_lbl)
        sub_empty.add(empty_lbl)
        sub_progress.add(self.progress_label)
        box_title.add(sub_title)
        box_title.add(sub_empty)
        box_title.add(sub_progress)
        main_box.add(box_title)

        # Input fields
        input_fields = [
            ("Maximum radius for RDF:", "Enter the maximum radius for radial distribution function", "textInput_radius"),
            ("Bin width for histogram:", "Enter bin width for histogram distribution", "textInput_bin_width"),
            ("Atom labels to be excluded:", "Enter labels to exclude (0 for none)", "textInput_atom_list"),
            ("Shell center atom label:", "Enter atom label at shell center", "textInput_shell_center"),
            ("Atomic symbol for g(r):", "Enter atomic symbol for g(r) calculation", "textInput_atom_symbol"),
            ("Cell lattices (a b c):", "Enter cell lattices a, b, c in Å", "textInput_cell_lattices"),
            ("x and y coordination number plot:", "Enter the x-limit and y-limit separated by white space", "textInput_axis"),
        ]
        for label_text, placeholder, attr in input_fields:
            box = toga.Box(style=box_style)
            lbl = toga.Label(label_text, style=label_style)
            txt = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr, txt)
            box.add(lbl)
            box.add(txt)
            main_box.add(box)

        # File selector and progress bar
        file_box = toga.Box(style=box_style)
        file_lbl = toga.Label("Select input file:", style=label_style)
        self.textInput_file = toga.TextInput(placeholder="Click Browse to select TRAJEC.xyz file", style=input_style)
        browse_btn = toga.Button("Browse", on_press=self.frames_counter, style=button_style)
        file_box.add(file_lbl)
        file_box.add(self.textInput_file)
        file_box.add(browse_btn)
        main_box.add(file_box)
        self.progress_bar = toga.ProgressBar(max=100)
        main_box.add(self.progress_bar)

        # Parameters display
        self.multi_line_text = toga.MultilineTextInput(style=Pack(flex=1, padding=(5, 0), font_size=11))
        self.multi_line_text.value = "\nRadial Distribution Function calculation based on the TRAJEC.xyz file"
        main_box.add(self.multi_line_text)

        # Action buttons
        btn_box = toga.Box(style=Pack(direction=ROW, alignment=CENTER, padding_top=5))
        self.btn_input_params = toga.Button("Read Params", style=button_style, on_press=self.read_params)
        self.btn_execute = toga.Button("RDF calculation", style=button_style, on_press=self.workflow)
        self.btn_help = toga.Button("Help", style=button_style, on_press=self.open_window_help)
        self.btn_close = toga.Button("Close", style=button_style, on_press=self.closeTopLevel)
        for btn in [self.btn_input_params, self.btn_execute, self.btn_help, self.btn_close]:
            btn_box.add(btn)
        main_box.add(btn_box)

        self.main_window.content = main_box
        self.main_window.show()

    def workflow(self, widget):
        self.read_params(widget)
        self.calcRDF()

    def open_window_help(self, widget):
        window = toga.Window(title="Instructions to carry out Radial Distribution Function")
        box = toga.Box(style=Pack(direction=COLUMN, flex=1))
        ml = toga.MultilineTextInput(style=Pack(font_size=11, padding=(5,5), flex=1))
        ml.value = HelpGqteaWin.help_RDF
        box.add(ml)
        window.content = box
        window.show()

    def closeTopLevel(self, widget):
        self.main_window.close()



    

