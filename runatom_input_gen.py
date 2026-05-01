import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, CENTER, BOLD
import subprocess
import platform
from pathlib import Path

# --- Help text extracted from INPUT_AE.txt and INPUT_GEN.txt ---

HELP_TEXT_AE = """
--- All-Electron (AE) Input Fields ---

[ ifae,ifpsp,ifprt,ifplw,ilogd ]
Controls the main program flow.
- ifae=1: Perform an all-electron calculation from scratch.
- ifpsp=0: Stop after the all-electron part.
- ifprt: Controls the level of print output (-3 to 1).
- ifplw: Generate data for wavefunction plots.
- ilogd: Number of l-values for which log-derivatives are calculated.

[ rlogd,emin,emax,nnt ]
Parameters for the log-derivative calculation grid.
- rlogd: Radius at which log-derivatives are calculated.
- emin, emax: Energy range for the calculation.
- nnt: Number of intervals in the energy mesh.

[ thresh,tol,damp,maxit ]
Convergence parameters for the self-consistent iteration.
- thresh: Convergence threshold for energy eigenvalues.
- tol: Tolerance for self-consistency.
- damp: Damping parameter for mixing.
- maxit: Maximum number of iterations.

[ z,xion,exfact ]
- z: Nuclear charge of the atom.
- xion: Net charge of the ion (0 for a neutral atom).
- exfact: Specifies the exchange-correlation functional (e.g., 0=LDA, 5=PBE).

[ rmax,aasf,bbsf ]
Parameters that define the radial mesh for the calculation.

[ ncspvs,irel ]
- ncspvs: Total number of states (core + valence).
- irel: Relativistic treatment (0=non-rel, 2=scalar-relativistic).

[ nnlz,wwnl,ee ] (One line for each state defined by ncspvs)
- nnlz: Quantum numbers n and l (e.g., 100 for 1s).
- wwnl: Electron occupation for this shell.
- ee: An initial guess for the energy eigenvalue.
"""

HELP_TEXT_GEN = """
--- Pseudopotential (GEN) Input Fields ---

NOTE: This assumes an all-electron file has already been generated.

[ ncores,nvales,nang ]
- ncores: Number of states to be treated as core.
- nvales: Number of states to be treated as valence.
- nang: l_max + 1 for the pseudization projectors.

[ keyps,ifpcor,rinner ]
- keyps: Type of pseudopotential (3 for Vanderbilt ultrasoft).
- ifpcor: Use partial core correction (1=yes, 0=no).
- rinner: Inner radius for pseudizing Q-functions.

[ nbeta,rcloc ]
- nbeta: Total number of beta projectors across all channels.
- rcloc: Cutoff radius for pseudizing the local potential.

[ rc ]
A list of cutoff radii for each angular momentum channel (s, p, d...).

[ lll,keyee,eeread,iptype ] (One line for each beta projector)
- lll: Angular momentum (l) of the projector.
- keyee: Source for the reference energy (e.g., >0 to use an AE eigenvalue).
- eeread: The reference energy value if keyee=0.
- iptype: Pseudization method (e.g., 2=soft, 4=norm-conserving).

[ npf,ptryc ] (Only if any iptype=2)
Parameters for the smoothness optimization of the pseudo-wavefunction.

[ lloc,keyee,eloc,iploctype ]
Defines how the local potential is constructed.
- lloc=-1: Simple polynomial matching.
- lloc>=0: Ensure correct scattering for channel lloc.

[ ifqopt,nqf,qtryc ]
Defines how Q-functions are pseudized.
- ifqopt: Optimization method (2 or 3 are common).
- nqf: Number of components in the polynomial expansion.
- qtryc: q-space cutoff for smoothness optimization.
"""


class RunatomGenerator:
    def __init__(self, *args):
        del args
        self.build_main_window()

    def build_main_window(self):
        """
        Constructs the main window, core widgets, and commands.
        """
        self.main_window = toga.Window(title="Runatom.x input generator", size=(800, 600))

        # --- Commands ---
        # Define a Help group if you want a custom one (optional)
        help_group = toga.Group("Help", order=1)

        self.help_command = toga.Command(
            self.show_help_window,
            text="Input Field Descriptions",
            tooltip="Show description for each input field",
            group=help_group,  # or group=toga.Group.HELP for the system Help menu
            order=1,
        )


        # --- Main Layout ---
        main_box = toga.Box(style=Pack(direction=COLUMN, margin=10))
        title_label = toga.Label(
            "Runatom.x Input Generator",
            style=Pack(text_align=CENTER, font_size=20, font_weight=BOLD, margin_bottom=10)
        )
        main_box.add(title_label)

        # --- Mode Selection (AE vs GEN) ---
        mode_box = toga.Box(style=Pack(direction=ROW, margin_bottom=10))
        mode_label = toga.Label("Select input file type:", style=Pack(margin_right=10))
        self.mode_selection = toga.Selection(
            items=["All-Electron (AE)", "Pseudopotential (GEN)"],
            on_change=self.on_mode_select,
            style=Pack(flex=1)
        )
        mode_box.add(mode_label)
        mode_box.add(self.mode_selection)
        main_box.add(mode_box)

        # --- Scrollable Form Area ---
        self.form_scroll_container = toga.ScrollContainer(horizontal=False)
        self.form_scroll_container.style.flex = 1
        self.form_box = toga.Box(style=Pack(direction=COLUMN, margin=10))
        self.form_scroll_container.content = self.form_box
        main_box.add(self.form_scroll_container)

        # --- Options Box (for Makefile switch) ---
        options_box = toga.Box(style=Pack(margin_top=10))
        self.makefile_switch = toga.Switch('Also generate Makefile', value=True)
        options_box.add(self.makefile_switch)
        main_box.add(options_box)

        # --- Styles ---
        button_style = Pack(margin=5, width=200)
        box_style = Pack(direction=ROW, margin=(0, 0, 5, 30))

        # --- Buttons Box (Row Layout) ---
        self.box_buttons = toga.Box(style=box_style)

        # --- Spacer to push buttons to the right ---
        spacer = toga.Box(style=Pack(flex=1))

        # --- Generate Button ---
        self.generate_button = toga.Button(
            "Generate",  # Give a meaningful label
            on_press=self.generate_file,
            style=button_style
        )

        # --- Help Button ---
        self.help_button = toga.Button(
            "Input Field Descriptions",
            on_press=self.show_help_window,
            style=button_style
        )

        # --- Add Widgets ---
        self.box_buttons.add(spacer)                # Spacer first
        self.box_buttons.add(self.generate_button)  # Then buttons
        self.box_buttons.add(self.help_button)

        # --- Add to main box ---
        main_box.add(self.box_buttons)


        self.main_window.content = main_box
        self.main_window.show()

        # Initialize with the AE form
        self.on_mode_select(self.mode_selection)

    async def show_help_window(self, widget):
        """
        Creates and shows a new window with help text corresponding to the
        currently selected mode.
        """
        mode = self.mode_selection.value
        if mode == "All-Electron (AE)":
            title = "All-Electron (AE) Help"
            content = HELP_TEXT_AE
        elif mode == "Pseudopotential (GEN)":
            title = "Pseudopotential (GEN) Help"
            content = HELP_TEXT_GEN
        else:
            return

        help_window = toga.Window(title=title, size=(700, 500))

        help_text_area = toga.MultilineTextInput(readonly=True, value=content, style=Pack(font_size=13,flex=1))
        help_box = toga.Box(style=Pack(direction=COLUMN, margin=10))
        help_box.add(help_text_area)
        
        help_window.content = help_box
        help_window.show()

    def on_mode_select(self, widget):
        """
        Clears the form and rebuilds it based on the selected mode (AE or GEN).
        """
        self.form_box.clear()
        if widget.value == "All-Electron (AE)":
            self.build_ae_form()
            self.generate_button.text = "Generate AE Input File"
        elif widget.value == "Pseudopotential (GEN)":
            self.build_gen_form()
            self.generate_button.text = "Generate GEN Input File"

    def build_ae_form(self):
        """
        Creates the input widgets for the All-Electron (AE) form.
        """
        self.ae_widgets = {}

        def create_ae_input(name, label_text, placeholder, default_value=""):
            box = toga.Box(style=Pack(direction=ROW, margin_bottom=5, align_items=CENTER))
            label = toga.Label(label_text, style=Pack(width=250, text_align="right", margin_right=10))
            text_input = toga.TextInput(value=default_value, placeholder=placeholder, style=Pack(flex=1))
            box.add(label)
            box.add(text_input)
            self.form_box.add(box)
            self.ae_widgets[name] = text_input
            return text_input

        # --- Build the form ---
        create_ae_input("line1", "ifae,ifpsp,ifprt,ifplw,ilogd:", "e.g., 1 0 0 0 3 (5 integers)", "1    0    0    0    3")
        create_ae_input("line2", "rlogd,emin,emax,nnt:", "e.g., 1.8 -2.4 1.6 40 (3f, 1i)", "1.80     -2.4       1.6       40")
        create_ae_input("line3", "thresh,tol,damp,maxit:", "e.g., 1.0d-10 1.0d-09 0.4 0 (2e, 1f, 1i)", "1.0d-10   1.0d-09   0.4        0")
        create_ae_input("title", "Title (Atom Symbol):", "e.g., h, pb", "h")
        create_ae_input("cfg", "Configuration (CFG):", "e.g., s1", "s1")
        create_ae_input("line5", "z,xion,exfact:", "e.g., 82. 0.0 5.0 (3 floats)", "1.0    0.0       5.0")
        create_ae_input("line6", "rmax,aasf,bbsf:", "e.g., 150.0 13.0 40.0 (3 floats)", "80.0       6.0      59.0")

        ncspvs_box = toga.Box(style=Pack(direction=ROW, margin_bottom=5, align_items=CENTER))
        ncspvs_label = toga.Label("ncspvs,irel (N States, Relativistic):", style=Pack(width=250, text_align="right", margin_right=10))
        self.ae_widgets['ncspvs_irel'] = toga.TextInput(value="1    2", placeholder="e.g., 15 2 (2 integers)", style=Pack(flex=1))
        self.ae_widgets['ncspvs_irel'].on_change = self.update_ae_states_inputs
        ncspvs_box.add(ncspvs_label)
        ncspvs_box.add(self.ae_widgets['ncspvs_irel'])
        self.form_box.add(ncspvs_box)
        
        self.ae_states_box = toga.Box(style=Pack(direction=COLUMN, margin_left=20))
        self.form_box.add(self.ae_states_box)
        self.update_ae_states_inputs(self.ae_widgets['ncspvs_irel'])

    def update_ae_states_inputs(self, widget):
        """
        Dynamically adds or removes input fields for atomic states based on `ncspvs`.
        """
        self.ae_states_box.clear()
        self.ae_widgets['states'] = []
        try:
            num_states = int(widget.value.strip().split()[0])
        except (ValueError, IndexError):
            num_states = 0
        
        for i in range(num_states):
            state_box = toga.Box(style=Pack(direction=ROW, margin_bottom=5, align_items=CENTER))
            state_label = toga.Label(f"State {i+1} (nnlz,wwnl,ee):", style=Pack(width=230, text_align="right", margin_right=10))
            default_val = "100  1.        -1.0" if i == 0 else ""
            state_input = toga.TextInput(value=default_val, placeholder="e.g., 100 2. -9000.0 (i, f, f)", style=Pack(flex=1))
            state_box.add(state_label)
            state_box.add(state_input)
            self.ae_states_box.add(state_box)
            self.ae_widgets['states'].append(state_input)
    
    def build_gen_form(self):
        """
        Creates the input widgets for the Pseudopotential Generation (GEN) form.
        """
        self.gen_widgets = {}

        def create_gen_input(name, label_text, placeholder, default_value=""):
            box = toga.Box(style=Pack(direction=ROW, margin_bottom=5, align_items=CENTER))
            label = toga.Label(label_text, style=Pack(width=250, text_align="right", margin_right=10))
            text_input = toga.TextInput(value=default_value, placeholder=placeholder, style=Pack(flex=1))
            box.add(label)
            box.add(text_input)
            self.form_box.add(box)
            self.gen_widgets[name] = text_input
            return text_input

        # Build the form using h_ps.adat.txt as a template
        self.form_box.add(toga.Label("General Settings", style=Pack(font_weight=BOLD, margin_top=10)))
        create_gen_input("line1", "ifae,ifpsp,ifprt,ifplw,ilogd:", "e.g., 0 2 1 1 3 (5 integers)", "0    2    1    1    3")
        create_gen_input("line2", "rlogd,emin,emax,nnt:", "e.g., 1.2 -2.4 1.6 80 (3f, 1i)", "1.2      -2.4       1.6       80")
        create_gen_input("line3", "thresh,tol,damp,maxit:", "e.g., 1.0d-11 1.0d-09 0.5 0", "1.0d-11   1.0d-09   0.5        0")
        create_gen_input("title", "Title (Atom Symbol):", "Must match AE run", "h")
        create_gen_input("cfg", "Configuration (CFG):", "e.g., s1", "s1")

        self.form_box.add(toga.Label("Pseudopotential Config", style=Pack(font_weight=BOLD, margin_top=10)))
        create_gen_input("states", "ncores,nvales,nang:", "e.g., 0 1 1 (3 integers)", "0    1    1")
        create_gen_input("bessel", "besrmax,besemin,besemax,besde:", "e.g., 10.0 20.0 40.0 10.0", "10.0      20.0      40.0      10.0")
        create_gen_input("pseudo_type", "keyps,ifpcor,rinner:", "e.g., 3 0 0.8 (2i, 1f)", "3    0   0.8")

        # Simplified dynamic sections for now
        self.form_box.add(toga.Label("Beta Projectors", style=Pack(font_weight=BOLD, margin_top=10)))
        create_gen_input("nbeta", "nbeta,rcloc:", "e.g., 1 0.8 (1i, 1f)", "1   0.8")
        create_gen_input("rc", "rc (for s, p, d...):", "e.g., 0.8 0.0 0.0", "0.8       0.0       0.0")
        self.form_box.add(toga.Label("Beta 1 (lll,keyee,eeread,iptype):", style=Pack(margin_left=20)))
        create_gen_input("beta1", "", "e.g., 0 1 0.0 2 (2i, 1f, 1i)", "0    1   0.0        2")

        self.form_box.add(toga.Label("Optimal Smoothness (if any iptype=2)", style=Pack(font_weight=BOLD, margin_top=10)))
        create_gen_input("smoothness", "npf,ptryc:", "e.g., 6 10.0 (1i, 1f)", "6  10.0")

        self.form_box.add(toga.Label("Local Potential", style=Pack(font_weight=BOLD, margin_top=10)))
        create_gen_input("local_pot", "lloc,keyee,eloc,iploctype:", "e.g., 1 0 0.0 1", "1    0   0.0        1")

        self.form_box.add(toga.Label("Q-Function Pseudization", style=Pack(font_weight=BOLD, margin_top=10)))
        create_gen_input("q_func", "ifqopt,nqf,qtryc:", "e.g., 3 8 10.0 (2i, 1f)", "3    8  10.0")

    async def generate_file(self, widget):
        """
        Gathers data, formats content, saves file(s), and asks to open the directory.
        """
        mode = self.mode_selection.value
        gen_makefile = self.makefile_switch.value
        
        try:
            if mode == "All-Electron (AE)":
                adat_content = self.format_ae_file()
                atom = self.ae_widgets['title'].value.strip()
                cfg = self.ae_widgets['cfg'].value.strip()
                filename = f"{atom}_ae_{cfg}.adat"
            elif mode == "Pseudopotential (GEN)":
                adat_content = self.format_gen_file()
                atom = self.gen_widgets['title'].value.strip()
                cfg = self.gen_widgets['cfg'].value.strip()
                filename = f"{atom}_ps.adat"
            else:
                return

            # Show save dialog to get the desired directory and filename
            file_path_obj = await self.main_window.dialog(
                toga.SaveFileDialog(
                    f"Save {mode} Input File",
                    suggested_filename=filename,
                    file_types=["adat", "txt"],
                )
            )

            if not file_path_obj:
                return  # User cancelled

            file_path = Path(file_path_obj)
            
            # Write the .adat file
            with open(file_path, 'w') as f:
                f.write(adat_content)

            saved_files_msg = [f"Successfully saved:\n- {file_path.name}"]

            # Write the Makefile if requested
            if gen_makefile:
                makefile_path = file_path.parent / "Makefile"
                makefile_content = self.format_makefile(atom, cfg)
                with open(makefile_path, 'w') as f:
                    f.write(makefile_content)
                saved_files_msg.append(f"- {makefile_path.name}")

            await self.main_window.dialog(
                toga.InfoDialog("Success", "\n".join(saved_files_msg))
            )

            # Ask to open the containing folder
            if await self.main_window.dialog(
                toga.QuestionDialog("Open Folder", "Do you want to open the save location?")
            ):
                await self.open_directory(file_path.parent)

        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Could not generate file: {e}")
            )
    
    async def open_directory(self, path):
        """Opens the specified directory in the system's file explorer."""
        system = platform.system()
        try:
            if system == "Windows":
                subprocess.Popen(["explorer", str(path)])
            elif system == "Darwin": # macOS
                subprocess.Popen(["open", str(path)])
            else: # Linux and other UNIX-like
                subprocess.Popen(["xdg-open", str(path)])
        except Exception as e:
            await self.main_window.dialog(
                toga.ErrorDialog("Error", f"Could not open directory: {e}")
            )

    def format_makefile(self, atom, cfg):
        """Generates the content for the Makefile."""
        return f"""
#################################################################
# This is the only section that usually needs to be modified
#
# Set parameters:
#   ATOM = chemical symbol of atom
#   CFG  = atomic or ionic configuration used for generation

ATOM= {atom}
CFG= {cfg}

#################################################################
#
# Standard operation:
# 
#   1.   make          Generate pseudopotential
#   2.                 Optionally, inspect potential by viewing *ps*out
#   3.   make install  Install pseudopotential file in target directory
#   4.   make clean    Clean up
#
# For other options, see below
#
#################################################################

# This should be root directory of a7.3.3 package:
ROOT= ../../..

PROG= ${{ROOT}}/Bin/runatom.x
POT=  ${{ROOT}}/Pot

AEG= ${{ATOM}}_ae_${{CFG}}
PSP= ${{ATOM}}_ps

default: pseudo

#---------------------------------------------------------------
# run all-electron calculations
#---------------------------------------------------------------

ae: ae_ground

ae_ground: ${{AEG}}.ae

${{AEG}}.ae: ${{AEG}}.adat
	${{PROG}} ${{AEG}}.adat ${{AEG}}.out ${{AEG}}.ae ${{AEG}}.atwf ${{AEG}}.logd dummy

#---------------------------------------------------------------
# generate pseudopotential
#---------------------------------------------------------------

pseudo: ${{PSP}}.uspp

${{PSP}}.uspp: ${{PSP}}.adat ${{AEG}}.ae
	${{PROG}} ${{PSP}}.adat ${{PSP}}.out ${{AEG}}.ae ${{PSP}}.atwf ${{PSP}}.logd ${{PSP}}.uspp

#---------------------------------------------------------------
# install pseudopotential in target directory
#---------------------------------------------------------------

install:
	cp ${{PSP}}.uspp ${{POT}}/`basename \\`pwd\\``.uspp
	- cp README ${{POT}}/`basename \\`pwd\\``.readme
	echo Installed to ${{POT}}/`basename \\`pwd\\``.uspp

# The above should have the effect of storing the potential
# as ${{POT}}/NAME.uspp where NAME is the last segment of the
# pathname of the current directory.

#---------------------------------------------------------------
# clean up
#---------------------------------------------------------------

clean:
	- rm *.out *.ae *.atwf *.logd *.uspp

# to clean up only pseudo outputs, but preserve results of AE runs:
psclean:
	- rm *_ps.out *_ps.atwf *_ps.logd *.uspp
"""

    def format_ae_file(self):
        """
        Gathers data from the AE form widgets and formats it into the .adat file format.
        """
        w = self.ae_widgets
        title_with_cfg = f"hydrogen {w['cfg'].value.strip()}" # Example title format
        lines = [
            f"{w['line1'].value:<30} ifae,ifpsp,ifprt,ifplw,ilogd (5i5)",
            f"{w['line2'].value:<35} rlogd,emin,emax,nnt (3f10.5,i5)",
            f"{w['line3'].value:<35} thresh,tol,damp,maxit (2e10.1,f10.5,i5)",
            f"{title_with_cfg:<30} title (a20)",
            f"{w['line5'].value:<30} z,xion,exfact (f7.2,2f10.5)",
            f"{w['line6'].value:<30} rmax,aasf,bbsf (3f10.5)",
            f"{w['ncspvs_irel'].value:<30} ncspvs,irel (2i5)"
        ]
        
        for state_input in w['states']:
            lines.append(f"{state_input.value:<30} nnlz,wwnl,ee (i4,f7.3,f14.6)")
        
        return "\n".join(lines)

    def format_gen_file(self):
        """
        Gathers data from the GEN form widgets and formats it. (Simplified)
        """
        w = self.gen_widgets
        title_with_cfg = f"hydrogen {w['cfg'].value.strip()}"
        lines = [
            f"{w['line1'].value:<30} ifae,ifpsp,ifprt,ifplw,ilogd (5i5)",
            f"{w['line2'].value:<35} rlogd,emin,emax,nnt (3f10.5,i5)",
            f"{w['line3'].value:<35} thresh,tol,damp,maxit (2e10.1,f10.5,i5)",
            f"{title_with_cfg:<30} title (a20)",
            f"{w['states'].value:<30} ncores,nvales,nang (3i5)",
            f"{w['bessel'].value:<40} besrmax,besemin,besemax,besde (4f10.5)",
            f"{w['pseudo_type'].value:<30} keyps,ifpcor,rinner (2i5,f10.5)",
            f"{w['nbeta'].value:<30} nbeta,rcloc (i5,f10.5)",
            f"{w['rc'].value:<30} rc (3f10.5)",
            f"{w['beta1'].value:<30} lll,keyee,eeread,iptype (2i5,f10.5,i5)",
            f"{w['smoothness'].value:<30} npf,ptryc (i5,f10.5)",
            f"{w['local_pot'].value:<30} lloc,keyee,eloc,iploctype (2i5,f10.5,i5)",
            f"{w['q_func'].value:<30} ifqopt,nqf,qtryc (2i5,f10.5)"
        ]

        return "\n".join(lines)



