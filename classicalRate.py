import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW,CENTER, LEFT
from displayPlots import DisplayPlots
import numpy as np

class ClassicalRate(DisplayPlots):
    def __init__(self):
        # Physical constants
        self.Boltzmann = 1.380649e-23       # Boltzmann constant in J/K
        self.gas_constant = 8.3144621       # Gas constant in J/(mol·K)
        self.planck = 6.62607015e-34        # Planck constant in J·s

    async def read_params(self, widget):
        """Read and validate the input parameters."""

        async def read_input(textInput, field_name, expected_type):
            value = textInput.value.strip()
            if not value:
                await self.main_window.dialog(toga.ErrorDialog(
                    "Error", f"Please input a valid value for {field_name}.")
                )
                return None
            try:
                if expected_type == list:
                    # Assuming atom labels are space-separated integers
                    labels = [float(label) for label in value.split()]
                    if len(labels) != 2:
                        raise ValueError("Please input exactly two atom labels.")
                    return labels
                else:
                    return expected_type(value)

            except ValueError as e:
                await self.main_window.dialog(toga.ErrorDialog(
                    "Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        # Read and validate inputs
        self.temperature = await read_input(
            self.textInput_temperature, "Temperature (K)", float
        )
        if self.temperature is None:
            return False

        self.activation_energy = await read_input(
            self.textInput_activation_energy, "Activation Energy Ea (kJ/mol)", float
        )
        if self.activation_energy is None:
            return False
        # Convert activation energy from kJ/mol to J/mol
        self.activation_energy *= 1000  # kJ/mol to J/mol

        self.entropy = await read_input(
            self.textInput_entropy, "Transition-State Entropy ΔS‡ (J/(mol·K))", float
        )
        if self.entropy is None:
            return False
        
        self.plot_range = await read_input(
            self.textInput_plot_range, "Transition-State Entropy ΔS‡ (J/(mol·K))", list
        )
        if self.plot_range is None:
            return False

        # Temperature range for plotting
        self.t_min, self.t_max = self.plot_range

        update_text = (
            f"{'Temperature for k(T)':<35} {self.temperature:>20} K\n"
            f"{'Activation Energy Ea':<35}{self.activation_energy / 1000:>20} kJ/mol\n"
            f"{'Transition-State Entropy ΔS‡':<35}{self.entropy:>20} J/(mol·K)\n"
        )
        self.multi_line_text.value = update_text
        return True

    def classical_rate_const(self, T, Ea, S):
        """Calculate the classical rate constant using Transition State Theory."""
        kB = self.Boltzmann
        h = self.planck
        R = self.gas_constant

        # Calculate rate constant
        rate_constant = (kB * T / h) * np.exp(S / R) * np.exp(-Ea / (R * T))
        return rate_constant

    def plot_rate_const(self, t_min, t_max, num_points=100):
        """Plot rate constant vs. temperature."""
        plot_range = np.linspace(t_min, t_max, num_points)
        Ea = self.activation_energy
        S = self.entropy

        rate_constants = [self.classical_rate_const(T, Ea, S) for T in plot_range]

        # Assuming save_plots is a method in DisplayPlots
        plot_number = 1
        x_axis = plot_range
        y_axis = rate_constants
        plot_xlabel = "Temperature (K)"
        plot_ylabel = "Rate Constant k(T) (s⁻¹)"
        plot_title = "Classical Rate Constant vs. Temperature"

        self.save_plots(plot_number, x_axis, y_axis, plot_xlabel, plot_ylabel, plot_title)

        output_file = f"{self.output_dir}/rate_constant_vs_temperature.dat"
        self.save_file(plot_range, rate_constants, output_file)
        
        # Plot lnK vs. 1/T and save to file
        self.plot_lnK_inv_T(plot_range, rate_constants)

    def plot_lnK_inv_T(self, temprt, rate_const):
        """Plot lnK vs. inverse of temperature and save to file."""

        # Calculate lnK
        lnK = [np.log(k) for k in rate_const]
        inv_T = [1 / T for T in temprt]

        # Assuming save_plots is a method in DisplayPlots
        plot_number = 2
        x_axis = inv_T
        y_axis = lnK
        plot_xlabel = "1/T (K⁻¹)"
        plot_ylabel = "ln(K) (s⁻¹)"
        plot_title = "Lnk vs 1/T"

        self.save_plots(plot_number, x_axis, y_axis, plot_xlabel, plot_ylabel, plot_title)

        output_file = f"{self.output_dir}/LnK_vs_invT.dat"
        self.save_file(inv_T, lnK, output_file)

    def save_file(self, plot_range, rate_constants, file_name):
        """Save the output to a file."""
        with open(file_name, "w") as f:
            for T in plot_range:
                for rate in rate_constants:
                    f.write(f"{T:.5} {rate:.5}\n")

    async def open_save_dir(self,widget):
        self.output_dir = await self.main_window.dialog(toga.SelectFolderDialog(
            "Select folder to save output files")
        )
        self.textInput_save_dir.value = self.output_dir

class ClassicalRateUI(ClassicalRate):
    def __init__(self,*args):
        super().__init__()
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        # Create the main window
        self.main_window = toga.Window(
            title="Classical Rate Constant Calculation",
            size=(700, 600),
        )

        # Define common styles
        heading_style = Pack(
            font_size=18, font_weight="bold", text_align=LEFT, padding=(0, 0, 10, 0)
        )
        label_style = Pack(
            padding=(0, 0, 5, 5), text_align=LEFT, width=200
        )
        input_style = Pack(flex=1, padding=(5, 5))
        button_style = Pack(padding=5, width=100)
        box_style = Pack(direction=ROW, align_items=CENTER, padding=(0, 0, 5, 0))

        # Main container
        main_box = toga.Box(style=Pack(direction=COLUMN, padding=20))

        # Title
        title_label = toga.Label("Classical Rate Constant Calculation", style=heading_style)
        main_box.add(title_label)

        # Input fields with labels
        input_fields = [
            (
                "Temperature (K):",
                "Enter the temperature in Kelvin for k(T) calculation",
                "textInput_temperature",
            ),
            (
                "Activation Energy Ea (kJ/mol):",
                "Enter the activation energy Ea in kJ/mol",
                "textInput_activation_energy",
            ),
            (
                "Transition-State Entropy ΔS‡ (J/(mol·K)):",
                "Enter the transition-state entropy ΔS‡ in J/(mol·K)",
                "textInput_entropy",
            ),
            (
                "Temperature range for plot:",
                "Enter the temperature range separated by space for plotting",
                "textInput_plot_range",
            ),
        ]

        for label_text, placeholder, attr_name in input_fields:
            box = toga.Box(style=box_style)
            label = toga.Label(label_text, style=label_style)
            textInput = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, textInput)
            box.add(label)
            box.add(textInput)
            main_box.add(box)

        box_save_dir = toga.Box(style=box_style)
        label_save_dir = toga.Label("Select directory:", style=label_style)
        self.textInput_save_dir = toga.TextInput(
            placeholder="Click Browse to select directory to save output files",
            style=input_style,
        )
        browse_save_dir = toga.Button(
            "Browse", on_press=self.open_save_dir, style=button_style
        )
        box_save_dir.add(label_save_dir)
        box_save_dir.add(self.textInput_save_dir)
        box_save_dir.add(browse_save_dir)
        main_box.add(box_save_dir)

        # Multi-line text input for output
        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, padding=(5, 0), font_size=11), readonly=True
        )
        self.multi_line_text.value = (
            "\n  CLASSICAL RATE CONSTANT CALCULATION\n"
        )
        main_box.add(self.multi_line_text)

        # Buttons at the bottom
        button_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, padding_top=5)
        )
        self.btn_read_params = toga.Button(
            "Read params", style=button_style, on_press=self.read_params
        )
        self.btn_calculate = toga.Button(
            "Calculate", style=button_style, on_press=self.workflow
        )
        self.btn_close = toga.Button(
            "Close", style=button_style, on_press=self.close_main_window
        )
        button_box.add(self.btn_read_params)
        button_box.add(self.btn_calculate)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        # Set the content of the main window
        self.main_window.content = main_box
        self.main_window.show()

    async def workflow(self, widget):
        success = await self.read_params(widget)  # Read and validate input parameters
        if not success:
            return

        rate_constant = self.classical_rate_const(
            self.temperature, self.activation_energy, self.entropy
        )

        self.multi_line_text.value += (
            f"\n K({self.temperature}K) = {rate_constant:.5e} s⁻¹\n"
        )

        # Plot rate constant vs. temperature
        self.plot_rate_const(self.t_min, self.t_max)
        self.display_plots()

    def close_main_window(self, widget):
        self.main_window.close()
