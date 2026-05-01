# This module contains classes and methods to analyze hydrogen bonds from molecular dynamics simulations.
import os, toga
import numpy as np
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER
from help import HelpGqteaWin
import statistics as st
from framesCounter import FramesCounter
from displayPlots import DisplayPlots


class HbondAnalyzer(FramesCounter,DisplayPlots):

    async def warning_function(self, msg):
        await self.main_window.dialog(toga.InfoDialog("Error", f"{msg}"))

    async def read_params(self, widget):
        #self.atufs = 0.02418884326505  # Atomic time unit in femtoseconds.

        async def read_input(text_input, field_name, expected_type):
            value = text_input.value.strip()
            if not value:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Please input a valid value for {field_name}.")
                )
                return None
            try:
                if expected_type == list:
                    # Assuming atom labels are space-separated integers
                    labels = [int(label) for label in value.split()]
                    if len(labels) != 2:
                        raise ValueError("Please input exactly two atom labels.")
                    return labels
                else:
                    return expected_type(value)
            except ValueError as e:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid format for {field_name}: {e}")
                )
                return None

        # Read and validate inputs
        self.max_hbond = await read_input(self.textInput_max_hbond, "Maximum hbond length", float)
        if self.max_hbond is None:
            return

        self.dt = await read_input(self.textInput_time_step, "Simulation time step", float)
        if self.dt is None:
            return

        self.sampling = await read_input(self.textInput_sampling_interval, "Interval to collect frames", int)
        if self.sampling is None:
            return

        self.sim_temp = await read_input(self.textInput_temperature, "simulation temperature", float)
        if self.sim_temp is None:
            return

        self.atom_donor = await read_input(self.textInput_atom_donor, "Atom donor label", int)
        if self.atom_donor is None:
            return        
        
        self.atom_H = await read_input(self.textInput_atom_H, "Atom H label", int)
        if self.atom_H is None:
            return        
        
        self.atom_acceptor = await read_input(self.textInput_atom_acceptor, "Atom acceptor label", int)
        if self.atom_acceptor is None:
            return

        self.hist_bin_hbonds = await read_input(self.textInput_bin_hbonds, "bin width", float)
        if self.hist_bin_hbonds is None:
            return        
        
        self.hist_bin_angleBonds = await read_input(self.textInput_bin_angleBonds, "bin width", float)
        if self.hist_bin_angleBonds is None:
            return        

        update_text = (
            f"Maximum hbond length --> {self.max_hbond}\n"
            f"Simulation time step --> {self.dt}\n"
            f"Interval sampling --> {self.sampling}\n"
            f"Simulation temperature --> {self.sim_temp}\n"
            f"Donor atom label --> {self.atom_donor}\n"
            f"Hydrogen atom label --> {self.atom_H}\n"
            f"Acceptor atom label --> {self.atom_acceptor}\n"
        )
        self.multi_line_text.value = update_text
   

    async def calcHBond(self,widget):
        """_summary_
        This function calculates the D-H and H--A bond lengths, where D-H stands for
        donor bond lengths and H--A stands for acptor bond lengths. D is the atom that
        donate H and A is the atom that acept H.
        """

        atufs=0.02418884326505
        time_step = self.dt
        sampling  = self.sampling
        # frame_number = 0
        

        D = self.atom_donor
        H = self.atom_H
        A = self.atom_acceptor
        self.hBond = [] #List to store simulation time, D-H bonds,H--A bonds, and D-H--A angles
        
        time_btf = time_step*sampling*atufs/1000.0  #Time between two frames collected

        elm1,elm2,elm3 = 0, 0, 0
        coords1, coords2, coords3 = None, None, None

        with open(self.trajec,"r") as traj:
            self.progress_bar.start()
            line = traj.readline() #Read number of atoms in trajectory
            t = 1.0
            frame_count = 1
            while True:
                line = traj.readline() #Read comment line of trajec
                if not line:
                    break

                for n in range(self.num_atoms):
                    line = traj.readline().split()
                    
                    line[1] = float(line[1])       # Positions  x in angstroms
                    line[2] = float(line[2])       # Positions  y in angstroms
                    line[3] = float(line[3])       # Positions  z in angstroms
                    if D == n + 1:
                        elm1 = line[0] 
                        coords1 = [line[1],line[2],line[3]]
                    elif H == n + 1:
                        elm2 = line[0]
                        coords2 = [line[1],line[2],line[3]]
                    elif A == n + 1:
                        elm3 = line[0]
                        coords3 = [line[1],line[2],line[3]]

                # Ensure that all coordinates have been assigned
                if coords1 is None or coords2 is None or coords3 is None:
                    msg = f"Could not find coordinates for one or more atoms: D, H, A."
                    await self.warning_function(msg)
                    
                a1 = np.array(coords1) #Stands for donor   (D) coordinates
                a2 = np.array(coords2) #Stands for Hydrogen (H) coordinates
                a3 = np.array(coords3) #Stands for Acceptor (A) coordinates
                vec1 = a1 - a2
                vec2 = a3 - a2

                #Calculate the vector norm
                vec1_norm = np.linalg.norm(vec1) #Distance between D and H
                vec2_norm = np.linalg.norm(vec2) #Distance between H and A

                # Normalize the vectors
                vec1_normalized = vec1 / vec1_norm
                vec2_normalized = vec2 / vec2_norm

                # Calculate the cosine of the angle
                cos_angle = np.dot(vec1_normalized, vec2_normalized)

                # Calculate the angle in radians
                angle_rad = np.arccos(cos_angle)

                # Convert the angle to degrees
                angle_deg = np.degrees(angle_rad)            
            
                simTime = t*time_btf

                self.hBond.append([simTime,vec1_norm,vec2_norm,angle_deg])

                t += 1.0
                line = traj.readline() #Read the number of atoms
                if not line:
                    break
                
                frame_count += 1

                if (frame_count % 400) == 0:
                    progress_bar_increment = (frame_count / self.total_frame_number)*100
                    self.progress_bar.value = progress_bar_increment

                if frame_count == self.total_frame_number:
                    self.progress_bar.value = 100
                
            self.progress_bar.stop()

        self.file0 = f"{self.output_dir}/hBond_{elm1}{D}_{elm2}{H}_{elm3}{A}.dat"  
        self.file1 = f"{self.output_dir}/hBond_FE_{elm1}{D}-{elm2}{H}.dat"
        self.file2 = f"{self.output_dir}/hBond_dist_func_{elm1}{D}_{elm2}{H}_{elm3}{A}.dat"
        self.file3 = f"{self.output_dir}/hBond_summary_{elm1}{D}_{elm2}{H}_{elm3}{A}.txt"
        self.file4 = f"{self.output_dir}/hBond_angle_FE_{elm1}{D}_{elm2}{H}_{elm3}{A}.dat"
        self.file5 = f"{self.output_dir}/hBondAngle_FE{elm1}{D}_{elm2}{H}_{elm3}{A}.dat"
        self.file6 = f"{self.output_dir}/hBond_FE_{elm2}{H}--{elm3}{A}.dat"

        with open(self.file3, 'w') as file3:
            file3.write(f"Number of frames:........................{self.total_frame_number}\n")
            file3.write(f"Number of atoms in each frame............{self.num_atoms}\n")            
            file3.write(f"Selected atoms:..........................{elm1}{D}-{elm2}{H}-{elm3}{A}\n")
        
        with open(self.file0, 'w') as file0:
            file0.write(f"time (ps)        {elm1}{D}-{elm2}{H}    {elm2}{H}--{elm3}{A}     angle (°)\n")
            for k in range(len(self.hBond)):
                p0 = self.hBond[k][0]
                p1 = self.hBond[k][1]
                p2 = self.hBond[k][2]
                p3 = self.hBond[k][3]
                file0.write(f"{p0:>10.4f}{p1:>14.5f}{p2:>14.5f}{p3:>14.5f}\n")
   
        self.multi_line_text.value += f"{elm1}{D}-{elm2}{H}--{elm3}{A}\n\n"

        # Save the bond length plot file and plot
        x  = [sublist[0] for sublist in self.hBond]
        y1 = [sublist[1] for sublist in self.hBond]
        y2 = [sublist[2] for sublist in self.hBond]
        y3 = [sublist[3] for sublist in self.hBond]

        plot_xlabel = "Simulation time (ps)"
        plot_ylabel = "Bond length D-H (Å)"
        plot_title = f"{elm1}{D}-{elm2}{H}"
        plot_number = 1
        self.save_plots(plot_number, x ,y1, plot_xlabel, plot_ylabel, plot_title)  

        plot_xlabel = "Simulation time (ps)"
        plot_ylabel = "Hydrogen Bond length H---A (Å)"
        plot_title = f"{elm2}{H}---{elm3}{A}"
        plot_number = 2
        self.save_plots(plot_number, x ,y2, plot_xlabel, plot_ylabel, plot_title) 

        plot_xlabel = "Simulation time (ps)"
        plot_ylabel = "Hydrogen Bond Angle D-H---A (°)"
        plot_title = f"{elm1}{D}-{elm2}{H}---{elm3}{A}"
        plot_number = 3
        self.save_plots(plot_number, x ,y3, plot_xlabel, plot_ylabel, plot_title) 
       
        #To be used futher
        self.elm1 = elm1
        self.elm2 = elm2
        self.elm3 = elm3
        self.D = D
        self.H = H
        self.A = A
    
    def hBond_dist_function(self):
        bin_width_bond = self.hist_bin_hbonds
        bin_width_angle = self.hist_bin_angleBonds
        max_r = self.max_hbond

        DH_bond  = [sublist[1] for sublist in self.hBond]
        AH_bond  = [sublist[2] for sublist in self.hBond]
        DHA_angle= [sublist[3] for sublist in self.hBond]

        bins_bond = np.arange(0,max_r,bin_width_bond)
        bins_angle = np.arange(0,180.0,bin_width_angle)

        DH_histo,DH_edges = np.histogram(DH_bond, bins_bond, density=True)
        AH_histo,AH_edges = np.histogram(AH_bond, bins_bond, density=True)
        DHA_histo,DHA_edges = np.histogram(DHA_angle, bins_angle, density=True)

        self.DH = DH_histo
        self.AH = AH_histo
        self.DHA = DHA_histo

        dim_DH  = len(DH_histo)
        dim_DHA = len(DHA_histo)

        x_bond = [0.0 for i in range(dim_DH)]
        x_hAngle = [0.0 for i in range(dim_DHA)]

        for j in range(dim_DH):
            x_bond[j] = (DH_edges[j] + DH_edges[j+1])/2.0

        for j in range(dim_DHA):
            x_hAngle[j] = (DHA_edges[j] + DHA_edges[j+1])/2.0

        #To be used in the next function
        self.x_bond = x_bond
        self.x_hAngle = x_hAngle

        avg_DH = st.mean(DH_bond)
        avg_AH = st.mean(AH_bond)
        avg_DHA = st.mean(DHA_angle)

        largest_DH = max(DH_bond)
        largest_AH = max(AH_bond)
        largest_DHA = max(DHA_angle)

        smallest_DH = min(DH_bond)
        smallest_AH = min(AH_bond)
        smallest_DHA = min(DHA_angle)

        var_DH = st.variance(DH_bond)
        var_AH = st.variance(AH_bond)
        var_DHA = st.variance(DHA_angle)

        stdev_DH = st.stdev(DH_bond)
        stdev_AH = st.stdev(AH_bond)
        stdev_DHA = st.stdev(DHA_angle)

        with open(self.file3,'a') as file3:
            file3.write(f"The largest {self.elm1}{self.D}-{self.elm2}{self.H} bond:.............{largest_DH:.4f} Å\n")
            file3.write(f"The smallest {self.elm1}{self.D}-{self.elm2}{self.H} bond:............{smallest_DH:.4f} Å\n")
            file3.write(f"Average {self.elm1}{self.D}-{self.elm2}{self.H} bond:.................{avg_DH:.4f} Å\n")
            file3.write(f"Variance of {self.elm1}{self.D}-{self.elm2}{self.H} bond:.............{var_DH:.4f} Å\n")
            file3.write(f"Stdev of {self.elm1}{self.D}-{self.elm2}{self.H}:.....................{stdev_DH:.4f}\n")

            file3.write(f"The largest {self.elm2}{self.H}--{self.elm3}{self.A} bond:.............{largest_AH:.4f} Å\n")
            file3.write(f"The smallest {self.elm2}{self.H}-{self.elm3}{self.A} bond:.............{smallest_AH:.4f} Å\n")
            file3.write(f"Average {self.elm2}{self.H}-{self.elm3}{self.A} bond:..................{avg_AH:.4f} Å\n")
            file3.write(f"Variance of {self.elm2}{self.H}-{self.elm3}{self.A} bond:..............{var_AH:.4f} Å\n")
            file3.write(f"Stdev of {self.elm2}{self.H}-{self.elm3}{self.A}:......................{stdev_AH:.4f}\n")

            file3.write(f"The largest {self.elm1}{self.D}-{self.elm2}{self.H}--{self.elm3}{self.A} angle:........{largest_DHA:.4f} Å\n")
            file3.write(f"The smallest {self.elm1}{self.D}-{self.elm2}{self.H}--{self.elm3}{self.A} angle:.......{smallest_DHA:.4f} Å\n")
            file3.write(f"Average {self.elm1}{self.D}-{self.elm2}{self.H}--{self.elm3}{self.A} angle:............{avg_DHA:.4f} Å\n")
            file3.write(f"Variance of {self.elm1}{self.D}-{self.elm2}{self.H}--{self.elm3}{self.A} angle:........{var_DHA:.4f} Å\n")
            file3.write(f"Stdev of {self.elm1}{self.D}-{self.elm2}{self.H}--{self.elm3}{self.A} angle:...........{stdev_DHA:.4f}\n")

        plot_xlabel = "D-H Bond length (Å)"
        plot_ylabel = f"{self.elm1}{self.D}-{self.elm2}{self.H} distribution"
        plot_title = f"{self.elm1}{self.D}-{self.elm2}{self.H}"
        plot_number = 4
        self.save_plots(plot_number, x_bond ,DH_histo, plot_xlabel, plot_ylabel, plot_title) 

        plot_xlabel = "A--H hydrogen bond length (Å)"
        plot_ylabel = f"{self.elm2}{self.H}---{self.elm3}{self.A}"
        plot_title = f"{self.elm2}{self.H}---{self.elm3}{self.A}"
        plot_number = 5
        self.save_plots(plot_number, x_bond, AH_histo, plot_xlabel, plot_ylabel, plot_title) 

        plot_xlabel = "D-H--A Hydrogen bond angle (°)"
        plot_ylabel = "Hydrogen bond angle distribution"
        plot_title = f"{self.elm1}{self.D}-{self.elm2}{self.H}---{self.elm3}{self.A} Distribution"
        plot_number = 6
        self.save_plots(plot_number, x_hAngle, DHA_histo, plot_xlabel, plot_ylabel, plot_title) 

    def hBonds_free_energy(self):
        """
        Calculate the free energies for hBonds and hBond angles using
        the distribution function stored in the DH_histo,AH_histo and DHA_histo files
        """
        DH_FE = []
        AH_FE = []
        DHA_FE = []
        R = 0.001987204 #Ideal gas universal constant in kcal/mol.k
        T = self.sim_temp
        dim_DH = len(self.DH)
        dim_AH = len(self.AH)
        dim_DHA = len(self.DHA)
    
        for n in range(dim_DH):
            if self.DH[n] != 0.0:
                FE = -R*T*np.log(self.DH[n])
                DH_FE.append([self.x_bond[n],FE])

        for n in range(dim_AH):
            if self.AH[n] != 0.0:
                FE = -R*T*np.log(self.AH[n])
                AH_FE.append([self.x_bond[n],FE])

        for n in range(dim_DHA):
            if self.DHA[n] != 0.0:
                FE = -R*T*np.log(self.DHA[n])
                DHA_FE.append([self.x_hAngle[n],FE])

        with open(self.file1,"w") as file1:
            file1.write(f"    r (Å)      {self.elm1}{self.D}-{self.elm2}{self.H}\n")
            for n in range(len(DH_FE)):
                file1.write(f"{DH_FE[n][0]:>10.5f}{DH_FE[n][1]:>14.5f}\n")
        
        with open(self.file6,"w") as file6:
            file6.write(f"   r (Å)        {self.elm2}{self.H}--{self.elm3}{self.A}\n")
            for n in range(len(AH_FE)):
                file6.write(f"{AH_FE[n][0]:>10.5f}{AH_FE[n][1]:>14.5f}\n")
        
        with open(self.file4,"w") as file4:
            file4.write(f" Angle (°)     {self.elm1}{self.D}-{self.elm2}{self.H}--{self.elm3}{self.A}\n")
            for n in range(len(DHA_FE)):
                file4.write(f"{DHA_FE[n][0]:>10.5f}{DHA_FE[n][1]:14.5f}\n")

        
        
        #plot Free Energy
        x_DH = [sublist[0] for sublist in DH_FE]
        x_AH = [sublist[0] for sublist in AH_FE]
        x_DHA = [sublist[0] for sublist in DHA_FE]

        y_DH = [sublist[1] for sublist in DH_FE]
        y_AH = [sublist[1] for sublist in AH_FE]
        y_DHA = [sublist[1] for sublist in DHA_FE]

        with open(self.file3,"a") as file3:
            min_y_DH = min(y_DH)
            min_y_DH_idx = y_DH.index(min_y_DH)
            file3.write(f"The smallest {self.elm1}{self.D}-{self.elm2}{self.H} bond FE is {round(min_y_DH,3)} at {round(x_DH[min_y_DH_idx],3)} Å\n")
            
            min_y_AH = min(y_AH)
            min_y_AH_idx = y_AH.index(min_y_AH)
            file3.write(f"The smallest {self.elm2}{self.H}--{self.elm3}{self.A} bond FE is {round(min_y_AH,3)} at {round(x_AH[min_y_AH_idx],3)} Å\n")
            
            min_y_DHA = min(y_DHA)
            min_y_DHA_idx = y_DHA.index(min_y_DHA)
            file3.write(f"The smallest {self.elm1}{self.D}-{self.elm2}{self.H}--{self.elm3}{self.A} angle FE is {round(min_y_DHA,3)} at {round(x_DHA[min_y_DHA_idx],3)}°\n")

        plot_xlabel = f"{self.elm1}{self.D}-{self.elm2}{self.H} bond length (Å)"
        plot_ylabel = f"{self.elm1}{self.D}-{self.elm2}{self.H} free energy (kcal/mol)"
        plot_title = f"{self.elm1}{self.D}-{self.elm2}{self.H} bond free energy"
        plot_number = 7
        self.save_plots(plot_number, x_DH, y_DH, plot_xlabel, plot_ylabel, plot_title) 

        plot_xlabel = f"{self.elm2}{self.H}--{self.elm3}{self.A} bond length (Å)"
        plot_ylabel = f"{self.elm2}{self.H}--{self.elm3}{self.A} free energy (kcal/mol)"
        plot_title  = f"{self.elm2}{self.H}--{self.elm3}{self.A} bond free energy"
        plot_number = 8
        self.save_plots(plot_number, x_AH, y_AH, plot_xlabel, plot_ylabel, plot_title) 

        plot_xlabel = f"{self.elm1}{self.D}-{self.elm2}{self.H}---{self.elm3}{self.A} hBond Angles"
        plot_ylabel = f"{self.elm1}{self.D}-{self.elm2}{self.H}---{self.elm3}{self.A} free energy (kcal/mol)"
        plot_title  = f"{self.elm1}{self.D}-{self.elm2}{self.H}---{self.elm3}{self.A} hBond Angles Free Energy"
        plot_number = 9
        self.save_plots(plot_number, x_DHA, y_DHA, plot_xlabel, plot_ylabel, plot_title) 

    def summary_multi_line_text(self):
        with open(self.file3,"r") as f:
            lines = f.readlines()
            for line in lines:
                self.multi_line_text.value += line

class HBondUI(HbondAnalyzer):
    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        # Create the main window
        self.main_window = toga.Window(
            title="Hydrogen Bond Analysis from Molecular Dynamics Simulations",
            size=(700, 600),
        )

        # Define common styles
        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style = Pack(margin=(5, 5), text_align=LEFT, width=200)
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=100)
        box_style = Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0))

        # Main container
        main_box = toga.Box(style=Pack(direction=COLUMN, margin=20))

        box_1 = toga.Box(style=Pack(direction=ROW, 
                        align_items=CENTER, margin=(0,0,10,0))
        )

        main_box.add(box_1)
        box_1a = toga.Box(style=Pack(width=400))
        box_1b = toga.Box(style=Pack(width=150))
        box_1c = toga.Box(style=Pack(width=100))
        box_1.add(box_1a)
        box_1.add(box_1b)
        box_1.add(box_1c)

        # Title
        title_label = toga.Label("Hydrogen Bond Analysis", 
                                 style=heading_style
        )
        empty_label = toga.Label("",style=heading_style
        )
        self.progress_label = toga.Label(" ",style=heading_style
        )

        box_1a.add(title_label)
        box_1b.add(empty_label)
        box_1c.add(self.progress_label)

        # Input fields with labels
        input_fields = [
            (
                "Maximum hydrogen bond length:",
                "Enter the maximum hydrogen bond length value",
                "textInput_max_hbond",
            ),
            (
                "Simulation Time Step (a.u.):",
                "Enter the simulation time step used in the simulation", 
                "textInput_time_step"),
            (
                "Sampling Interval (frames):",
                "Enter the number of skipped frames to collect frames",
                "textInput_sampling_interval",
            ),
            (
                "Simulation Temperature (K):",
                "Enter the temperature used in the simulation run",
                "textInput_temperature",
            )

        ]

        for label_text, placeholder, attr_name in input_fields:
            box = toga.Box(style=box_style)
            label = toga.Label(label_text, style=label_style)
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            box.add(label)
            box.add(text_input)
            main_box.add(box)

        hist_box = toga.Box(style=box_style)
        hist_label = toga.Label("Histogram Bin Widths (Å):",style=label_style)
        self.textInput_bin_hbonds = toga.TextInput(placeholder="Bin width for hbond histogram", style= input_style)
        self.textInput_bin_angleBonds = toga.TextInput(placeholder="Bin width for hbond angle histogram", style= input_style)
        hist_box.add(hist_label)
        hist_box.add(self.textInput_bin_hbonds)
        hist_box.add(self.textInput_bin_angleBonds)
        main_box.add(hist_box)

        atom_box = toga.Box(style=Pack(direction=ROW))
        atom_labels = toga.Label("Atom labels:",style=label_style)
        self.textInput_atom_donor = toga.TextInput(
                                     placeholder="Hydrogen donor atom label",
                                     style=Pack(width=170,margin=(0,5,0,5))
        )
        self.textInput_atom_H = toga.TextInput(
                                placeholder="Hydrogen atom label",
                                style=Pack(width=130,margin=(0,5,0,5))
        )
        self.textInput_atom_acceptor = toga.TextInput(
                                       placeholder="Hydrogen acceptor atom label",
                                       style=Pack(width=170,margin=(0,5,0,5))
        )

        atom_box.add(atom_labels)
        atom_box.add(self.textInput_atom_donor)
        atom_box.add(self.textInput_atom_H)
        atom_box.add(self.textInput_atom_acceptor)

        main_box.add(atom_box)

        # File selection button
        file_box = toga.Box(style=box_style)
        file_label = toga.Label("Select Trajectory File:", style=label_style)
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select TRAJEC.xyz file", style=input_style
        )
        browse_button = toga.Button(
            "Browse", on_press=self.frames_counter, style=button_style
        )
        progress_box = toga.Box(style=Pack(direction=COLUMN,height=5,margin=(0,5,5,5)))
        self.progress_bar = toga.ProgressBar(max=100)
        progress_box.add(self.progress_bar)

        file_box.add(file_label)
        file_box.add(self.textInput_file)
        file_box.add(browse_button)
        main_box.add(file_box)
        main_box.add(progress_box)


        # Multi-line text input for help or output
        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(15,5,5,5), font_size=12)
        )
        self.multi_line_text.value = HelpGqteaWin.help_HB
        main_box.add(self.multi_line_text)

        # Buttons at the bottom
        button_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin=(5,5,5,5))
        )
        self.btn_execute = toga.Button(
            "Analyze", style=button_style, on_press=self.workflow
        )
        self.btn_close = toga.Button(
            "Close", style=button_style, on_press=self.closeTopLevel
        )
        button_box.add(self.btn_execute)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        # Set the content of the main window
        self.main_window.content = main_box
        self.main_window.show()

    async def workflow(self, widget):
        await self.read_params(widget)
        if not hasattr(self, "trajec") or not self.trajec:
            await self.warning_function("No trajectory file selected.")
            return
        
        await self.calcHBond(widget)
        self.hBond_dist_function()
        self.hBonds_free_energy()
        self.summary_multi_line_text()
        self.display_plots()

    def closeTopLevel(self, widget):
        self.main_window.close()

