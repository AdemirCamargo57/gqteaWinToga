import toga, os, asyncio
import numpy as np
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from help import HelpGqteaWin
from framesCounter import FramesCounter
from displayPlots import DisplayPlots


class RadialAnalyser(FramesCounter, DisplayPlots):
    async def read_params(self, widget) -> bool:
        """Read and validate the input parameters.

        Returns True when every parameter is present and valid; otherwise shows
        an error dialog and returns False. Callers must check the return value
        before running the calculation.
        """

        async def read_input(text_input, field_name, expected_type):
            value = text_input.value.strip()
            if not value:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Please input a valid value for {field_name}.")
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
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        async def read_optional_float(text_input, field_name, default):
            """Like ``read_input`` for a single float, but a blank field falls
            back to ``default`` instead of being treated as an error."""
            value = text_input.value.strip()
            if not value:
                return default
            try:
                return float(value)
            except ValueError as e:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        async def fail(message):
            await self.main_window.dialog(toga.InfoDialog("Error", message))
            return False

        # A trajectory must be loaded first: validation needs the atom count and
        # the calculation needs the frame count.
        if (
            not getattr(self, "trajec", None)
            or not hasattr(self, "num_atoms")
            or not hasattr(self, "total_frame_number")
        ):
            return await fail(
                "Please select an input TRAJEC.xyz file with the Browse button "
                "before reading parameters."
            )

        # Read and validate inputs. The cell lattices are read first because the
        # minimum-image scheme caps the usable radius at half the smallest box
        # length -- which is also the default applied when the radius is blank.
        self.cell_lattices = await read_input(
            self.textInput_cell_lattices, "cell lattices a, b, and c", "float_list"
        )
        if self.cell_lattices is None:
            return False
        if len(self.cell_lattices) != 3:
            return await fail("Enter exactly three cell lattices: a b c.")
        if any(length <= 0 for length in self.cell_lattices):
            return await fail("Cell lattices must be positive.")

        # The periodic-image scheme is only valid up to half the smallest box length.
        min_half = min(self.cell_lattices) / 2.0

        # Maximum radius for the RDF. A blank field defaults to the maximum
        # possible value, min(a, b, c) / 2, imposed by the minimum-image scheme.
        self.radius = await read_optional_float(
            self.textInput_radius,
            "radius for radial distribution function calculation",
            min_half,
        )
        if self.radius is None:
            return False
        if self.radius <= 0:
            return await fail("Maximum radius must be greater than zero.")
        if self.radius > min_half:
            return await fail(
                f"Maximum radius ({self.radius}) must not exceed half the smallest "
                f"cell lattice ({min_half:.4f} Angstrom) for the periodic-image "
                f"method to be valid."
            )

        self.bin_width = await read_input(
            self.textInput_bin_width, "bin width for histogram distribution function", float
        )
        if self.bin_width is None:
            return False
        if self.bin_width <= 0:
            return await fail("Bin width must be greater than zero.")
        if self.bin_width > self.radius:
            return await fail("Bin width must not exceed the maximum radius.")

        delete_list = await read_input(
            self.textInput_atom_list,
            "list of atom labels to be excluded (0 for none)",
            "int_list",
        )
        if delete_list is None:
            return False
        if delete_list == [0]:
            self.delete_atom_list = []
        elif any(label < 1 or label > self.num_atoms for label in delete_list):
            return await fail(
                f"Atom labels to exclude must be between 1 and {self.num_atoms} "
                f"(or 0 for none)."
            )
        else:
            self.delete_atom_list = delete_list

        self.shell_center = await read_input(
            self.textInput_shell_center, "atom label at the center of the shell", int
        )
        if self.shell_center is None:
            return False
        if not 1 <= self.shell_center <= self.num_atoms:
            return await fail(
                f"Shell-center atom label must be between 1 and {self.num_atoms}."
            )

        self.rdf_atom_symbol = await read_input(
            self.textInput_atom_symbol, "atomic symbol for g(r) calculation", str
        )
        if self.rdf_atom_symbol is None:
            return False

        # X/Y limits for the coordination-number plot. Both fields are optional:
        # a blank field falls back to its default (x = lattice a / 2, y = 10) so
        # the plot always has sensible bounds.
        x_limit = await read_optional_float(
            self.textInput_xlim, "x-axis range limit", self.cell_lattices[0] / 2.0
        )
        if x_limit is None:
            return False
        y_limit = await read_optional_float(
            self.textInput_ylim, "y-axis range limit", 10.0
        )
        if y_limit is None:
            return False
        if x_limit <= 0 or y_limit <= 0:
            return await fail("Coordination-number plot limits must be positive.")
        self.axis = [x_limit, y_limit]

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
        if self.rho <= 0:
            return await fail(
                f"No atoms with symbol '{self.rdf_atom_symbol}' were found in the "
                f"first frame; cannot compute g(r)."
            )

        return True

    def ideal_density(self) -> None:
        """Compute the ideal number density ``self.rho`` (atoms of the target
        symbol per unit cell volume) from the first frame, used to normalise g(r)."""
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

    async def frames_counter(self, widget):
        """Load the trajectory and count its frames, reporting progress as text
        in the status box.

        This overrides the shared ``FramesCounter.frames_counter`` so this tool
        can report status in the ``MultilineTextInput`` instead of the removed
        progress label (leaving the mixin unchanged for every other tool).
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
                    if not title_line:
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

    async def calc_rdf(self):
        """Calculate the Radial Distribution Function g(r).

        Periodic boundaries are handled with the vectorized minimum-image
        convention, which is valid for radius <= min(a, b, c) / 2 (enforced in
        read_params) and also tolerates unwrapped input coordinates.
        """
        # Setup parameters and data structures
        symbol = self.rdf_atom_symbol
        center_idx = self.shell_center - 1  # Use 0-based indexing
        bin_width = self.bin_width
        radius = self.radius
        num_bins = int(radius / bin_width)

        bins = np.linspace(0, radius, num_bins + 1)
        r = 0.5 * (bins[:-1] + bins[1:])  # Radius values are the bin centers

        total_histo = np.zeros(num_bins, dtype=np.float64)
        frame_count = 0

        # Orthorhombic box lengths for the minimum-image convention.
        box = np.array(self.cell_lattices, dtype=np.float64)

        # Running status log shown in the text box (this tool has no progress bar).
        center_symbol = None
        center_label = f"atom {self.shell_center}"
        status_header = (
            f"Reading trajectory: {os.path.basename(self.trajec)}\n"
            f"  atoms: {self.num_atoms}   frames: {self.total_frame_number}\n\n"
            f"Validating parameters... OK\n"
        )

        def calc_status(done):
            return (
                status_header
                + f"Calculating g(r)  [{center_label} --- {symbol}] ...\n"
                + f"  processed {done} / {self.total_frame_number} frames"
            )

        self.multi_line_text.value = calc_status(0)

        # 3. Process trajectory frame by frame
        malformed = False
        with open(self.trajec, "r") as traj_file:
            while True:
                line_atoms = traj_file.readline()
                if not line_atoms:
                    break  # Clean end of file
                if not line_atoms.strip():
                    break  # Trailing blank line(s): treat as end of data

                try:
                    num_atoms = int(line_atoms.strip())
                except ValueError:
                    malformed = True
                    break  # Header line is not a valid atom count

                traj_file.readline()  # Skip comment line
                frame_lines = [traj_file.readline() for _ in range(num_atoms)]
                if any(not line for line in frame_lines):
                    malformed = True
                    break  # Truncated final frame (fewer atom lines than declared)

                try:
                    atom_data = np.array(
                        [(parts[0], float(parts[1]), float(parts[2]), float(parts[3]))
                         for parts in (line.split() for line in frame_lines)],
                        dtype=[('symbol', 'U10'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8')]
                    )
                except (ValueError, IndexError):
                    malformed = True
                    break  # Malformed atom line within the frame

                # 4. Coordinates of the shell-center atom.
                center_coords = np.array([atom_data['x'][center_idx], atom_data['y'][center_idx], atom_data['z'][center_idx]])

                if center_symbol is None:  # remember the center element for the status/summary
                    center_symbol = atom_data['symbol'][center_idx]
                    center_label = f"{center_symbol}{self.shell_center}"

                # Build a mask selecting the target-symbol neighbours, excluding
                # the user-specified atoms and the center atom itself.
                include_mask = np.ones(num_atoms, dtype=bool)
                if self.delete_atom_list:
                    delete_indices = np.array(self.delete_atom_list) - 1
                    include_mask[delete_indices] = False
                include_mask[center_idx] = False
                target_mask = include_mask & (atom_data['symbol'] == symbol)

                target_coords = np.vstack([
                    atom_data['x'][target_mask],
                    atom_data['y'][target_mask],
                    atom_data['z'][target_mask],
                ]).T

                if target_coords.shape[0] == 0:
                    frame_count += 1
                    continue

                # 5. Minimum-image displacement to the center atom. Wrapping each
                # component into [-L/2, L/2] applies periodic boundaries correctly
                # (and tolerates unwrapped coordinates) for radius <= L/2.
                delta = target_coords - center_coords
                delta -= box * np.round(delta / box)
                distances = np.sqrt(np.sum(delta**2, axis=1))

                # 6. Bin the distances that fall within the requested radius.
                distances_in_range = distances[distances <= radius]
                histo, _ = np.histogram(distances_in_range, bins=bins)
                total_histo += histo

                frame_count += 1
                if frame_count % 100 == 0:
                    self.multi_line_text.value = calc_status(frame_count)
                    await asyncio.sleep(0)  # Yield so the UI can repaint

        if frame_count == 0:
            await self.main_window.dialog(
                toga.InfoDialog("Error", "No frames were processed.")
            )
            return
        if malformed:
            await self.main_window.dialog(
                toga.InfoDialog(
                    "Warning",
                    f"A malformed or truncated frame was encountered; "
                    f"processed {frame_count} of {self.total_frame_number} frames.",
                )
            )

        # 8. Normalize and compute g(r)
        avg_histo = total_histo / frame_count
        shell_volumes = (4.0/3.0) * np.pi * (bins[1:]**3 - bins[:-1]**3)
        shell_volumes[shell_volumes < 1e-9] = 1.0 # Avoid division by zero
        
        real_density = avg_histo / shell_volumes
        g_r = real_density / self.rho
        coordination_num = np.cumsum(avg_histo)

        # 9. Save the tabulated data next to the trajectory.
        output_file = os.path.join(
            self.output_dir, f"RDF_{center_symbol}{self.shell_center}_{symbol}.dat"
        )
        with open(output_file, "w") as out_file:
            out_file.write("      r              g(r)             Integral\n\n")
            for i in range(num_bins):
                out_file.write(f"{r[i]:>10.5f}{g_r[i]:>20.7f}{coordination_num[i]:>20.7f}\n")

        # 10. Hand both curves to the shared interactive plot viewer. save_png is
        # False so no static PNG image files are left next to the trajectory; the
        # figures are shown through the interactive viewer only.
        self.save_plots(
            1, r, g_r,
            "r (Å)",
            f"g(r) ({center_symbol}{self.shell_center}---{symbol})",
            "Radial Distribution Function g(r)",
            save_png=False,
        )
        self.save_plots(
            2, r, coordination_num,
            "r (Å)",
            "Coordination number",
            "Coordination Number",
            xlim=(0, self.axis[0]),
            ylim=(0, self.axis[1]),
            save_png=False,
        )
        self.display_plots()

        # 11. Write the completion summary to the status box.
        excluded = " ".join(str(a) for a in self.delete_atom_list) if self.delete_atom_list else "None"
        a, b, c = self.cell_lattices
        peak_idx = int(np.argmax(g_r))
        self.multi_line_text.value = (
            "RADIAL DISTRIBUTION FUNCTION - COMPLETED\n"
            "----------------------------------------\n"
            f"Center atom:        {center_symbol}{self.shell_center}\n"
            f"Target symbol:      {symbol}\n"
            f"Frames processed:   {frame_count} / {self.total_frame_number}\n"
            f"Max radius / bin:   {self.radius} / {self.bin_width} Angstrom\n"
            f"Excluded atoms:     {excluded}\n"
            f"Cell (a b c):       {a} {b} {c} Angstrom\n"
            f"Cell volume:        {self.volume:.2f} Angstrom^3\n"
            f"Density (rho):      {self.rho:.4f} Angstrom^-3\n"
            f"First g(r) peak:    r = {r[peak_idx]:.3f}, g = {g_r[peak_idx]:.3f}\n"
            f"Coordination({self.radius}): {coordination_num[-1]:.3f}\n"
            f"Output: {os.path.basename(output_file)}\n"
        )


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
        box_style = Pack(direction=ROW, align_items=CENTER, padding=(0, 0, 5, 0))

        main_box = toga.Box(style=Pack(direction=COLUMN, padding=20))

        # Title
        box_title = toga.Box(style=Pack(direction=ROW, align_items=CENTER, padding=(0, 0, 10, 0)))
        title_lbl = toga.Label("Radial Distribution Function g(r)", style=heading_style)
        box_title.add(title_lbl)
        main_box.add(box_title)

        # Input fields: (label, placeholder, attribute, default value)
        input_fields = [
            ("Maximum radius for RDF:", "Maximum radius for RDF (leave blank for default: half the smallest cell lattice)", "textInput_radius", ""),
            ("Bin width for histogram:", "Enter bin width for histogram distribution", "textInput_bin_width", "0.01"),
            ("Atom labels to be excluded:", "Enter labels to exclude (0 for none)", "textInput_atom_list", ""),
            ("Shell center atom label:", "Enter atom label at shell center", "textInput_shell_center", ""),
            ("Atomic symbol for g(r):", "Enter atomic symbol for g(r) calculation", "textInput_atom_symbol", ""),
            ("Cell lattices (a b c):", "Enter cell lattices a, b, c in Å", "textInput_cell_lattices", ""),
            ("X-axis range limit for RDF plot:", "Enter the x-max for the g(r) plot (leave blank for default: lattice_a / 2)", "textInput_xlim", ""),
            ("Y-axis range limit for RDF plot:", "Enter the y-max for the g(r) plot (leave blank for default: 10)", "textInput_ylim", ""),
        ]
        for label_text, placeholder, attr, default in input_fields:
            box = toga.Box(style=box_style)
            lbl = toga.Label(label_text, style=label_style)
            txt = toga.TextInput(value=default, placeholder=placeholder, style=input_style)
            setattr(self, attr, txt)
            box.add(lbl)
            box.add(txt)
            main_box.add(box)

        # File selector
        file_box = toga.Box(style=box_style)
        file_lbl = toga.Label("Select input file:", style=label_style)
        self.textInput_file = toga.TextInput(placeholder="Click Browse to select TRAJEC.xyz file", style=input_style)
        browse_btn = toga.Button("Browse", on_press=self.frames_counter, style=button_style)
        file_box.add(file_lbl)
        file_box.add(self.textInput_file)
        file_box.add(browse_btn)
        main_box.add(file_box)

        # Parameters display
        self.multi_line_text = toga.MultilineTextInput(style=Pack(flex=1, padding=(5, 0), font_size=11))
        self.multi_line_text.value = "\nRadial Distribution Function calculation based on the TRAJEC.xyz file"
        main_box.add(self.multi_line_text)

        # Action buttons
        btn_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, padding_top=5))
        self.btn_input_params = toga.Button("Read Params", style=button_style, on_press=self.read_params)
        self.btn_execute = toga.Button("RDF calculation", style=button_style, on_press=self.workflow)
        self.btn_help = toga.Button("Help", style=button_style, on_press=self.open_window_help)
        self.btn_close = toga.Button("Close", style=button_style, on_press=self.closeTopLevel)
        for btn in [self.btn_input_params, self.btn_execute, self.btn_help, self.btn_close]:
            btn_box.add(btn)
        main_box.add(btn_box)

        self.main_window.content = main_box
        self.main_window.show()

    async def workflow(self, widget) -> None:
        """Read and validate parameters, then run the g(r) calculation."""
        if not await self.read_params(widget):
            return
        await self.calc_rdf()

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



    

