
class AtomicData:
    atomic_masses = {
        "H": 1.00784,   # Hydrogen
        "He": 4.002602, # Helium
        "Li": 6.94,     # Lithium
        "Be": 9.0122,   # Beryllium
        "B": 10.81,     # Boron
        "C": 12.011,    # Carbon
        "N": 14.007,    # Nitrogen
        "O": 15.999,    # Oxygen
        "F": 18.998,    # Fluorine
        "Ne": 20.180,   # Neon
        "Na": 22.990,   # Sodium
        "Mg": 24.305,   # Magnesium
        "Al": 26.982,   # Aluminium
        "Si": 28.085,   # Silicon
        "P": 30.9738,   # Phosphorus
        "S": 32.06,     # Sulfur
        "Cl": 35.45,    # Chlorine
        "Ar": 39.95,    # Argon
        "K": 39.098,    # Potassium
        "Ca": 40.078,   # Calcium
        "Sc": 44.956,   # Scandium
        "Ti": 47.867,   # Titanium
        "V": 50.9415,   # Vanadium
        "Cr": 51.9961,  # Chromium
        "Mn": 54.938,   # Manganese
        "Fe": 55.845,   # Iron
        "Co": 58.933,   # Cobalt
        "Ni": 58.693,   # Nickel
        "Cu": 63.546,   # Copper
        "Zn": 65.38,    # Zinc
        "Ga": 69.723,   # Gallium
        "Ge": 72.63,    # Germanium
        "As": 74.9216,  # Arsenic
        "Se": 78.971,   # Selenium
        "Br": 79.904,   # Bromine
        "Kr": 83.798,   # Krypton
        "Rb": 85.468,   # Rubidium
        "Sr": 87.62,    # Strontium
        "Y": 88.906,    # Yttrium
        "Zr": 91.224,   # Zirconium
        "Nb": 92.906,   # Niobium
        "Mo": 95.95,    # Molybdenum
        "Tc": 98,       # Technetium
        "Ru": 101.07,   # Ruthenium
        "Rh": 102.91,   # Rhodium
        "Pd": 106.42,   # Palladium
        "Ag": 107.87,   # Silver
        "Cd": 112.41,   # Cadmium
        "In": 114.82,   # Indium
        "Sn": 118.71,   # Tin
        "Sb": 121.76,   # Antimony
        "Te": 127.6,    # Tellurium
        "I": 126.90,    # Iodine
        "Xe": 131.29,   # Xenon
        "Cs": 132.91,   # Cesium
        "Ba": 137.33,   # Barium
        "La": 138.91,   # Lanthanum
        "Ce": 140.12,   # Cerium
        "Pr": 140.91,   # Praseodymium
        "Nd": 144.24,   # Neodymium
        "Pm": 145,      # Promethium
        "Sm": 150.36,   # Samarium
        "Eu": 151.96,   # Europium
        "Gd": 157.25,   # Gadolinium
        "Tb": 158.93,   # Terbium
        "Dy": 162.50,   # Dysprosium
        "Ho": 164.93,   # Holmium
        "Er": 167.26,   # Erbium
        "Tm": 168.93,   # Thulium
        "Yb": 173.05,   # Ytterbium
        "Lu": 174.97,   # Lutetium
        "Hf": 178.49,   # Hafnium
        "Ta": 180.95,   # Tantalum
        "W": 183.84,    # Tungsten
        "Re": 186.21,   # Rhenium
        "Os": 190.23,   # Osmium
        "Ir": 192.22,   # Iridium
        "Pt": 195.08,   # Platinum
        "Au": 196.97,   # Gold
        "Hg": 200.59,   # Mercury
        "Tl": 204.38,   # Thallium
        "Pb": 207.2,    # Lead
        "Bi": 208.98,   # Bismuth
        "Po": 209,      # Polonium
        "At": 210,      # Astatine
        "Rn": 222,      # Radon
        "Fr": 223,      # Francium
        "Ra": 226,      # Radium
        "Ac": 227,      # Actinium
        "Th": 232.04,   # Thorium
        "Pa": 231.04,   # Protactinium
        "U": 238.03,    # Uranium
        "Np": 237,      # Neptunium
        "Pu": 244,      # Plutonium
        "Am": 243,      # Americium
        "Cm": 247,      # Curium
        "Bk": 247,      # Berkelium
        "Cf": 251,      # Californium
        "Es": 252,      # Einsteinium
        "Fm": 257,      # Fermium
        "Md": 258,      # Mendelevium
        "No": 259,      # Nobelium
        "Lr": 262,      # Lawerencium
        "Rf": 267,      # Rutherfordium
        "Db": 270,      # Dubnium
        "Sg": 271,      # Seaborgium
        "Bh": 274,      # Bohrium
        "Hs": 277,      # Hassium
        "Mt": 278,      # Meitnerium
        "Ds": 281,      # Darmstadtium
        "Rg": 282,      # Roentgenium
        "Cn": 285,      # Copernicium
        "Nh": 286,      # Nihonium
        "Fl": 289,      # Flerovium
        "Mc": 290,      # Moscovium
        "Lv": 293,      # Livermorium
        "Ts": 294,      # Tennessine
        "Og": 294,      # Oganesson
    }
    atomic_numbers = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
        "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
        "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
        "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
        "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
        "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
        "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
        "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
        "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74,
        "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
        "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
        "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98,
        "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105,
        "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112,
        "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
    }

    lang = {"C":"P","O":"P","H":"S","N":"P",
            "Ca":"D","P":"D","Fe":"D","F":"P","S":"P",
            "Cl":"D","Br":"D","I":"D"}

class Fonts:
    font_1 = {'color':  'darkred','weight': 'normal','size': 16}
    font_2 = {'color':  'darkred','weight': 'normal','size': 12}


class HelpGqteaWin:
    help_cpmd_input = f"""
STEPS TO FOLLOW TO GENERATE CPMD INPUT FILES:
1. Enter a 'Prefix file name' to be used for your generated input files (e.g., VitC).
2. Input the 'Total charge on the system' (e.g., 0).
3. Define the 'Periodic box parameters' by providing a, b, c, cosα, cosβ, and cosγ 
   separated by white space.
4. Input the 'Energy Cutoff (Ry)' (e.g., 25).
5. Input the 'Dual for ρ expansion' in planewaves (e.g., 8).
6. Select the appropriate 'Cell symmetry' from the dropdown menu 
   (e.g., 1 - CUBIC a=b=c α=β=γ=90°).
7. Click the 'Browse' button to select your starting molecular geometry file 
   in xyz format.

Here is an example input file in .xyz format:

5
Comment line or blank line
P     0.112   -0.092    0.064  
O     0.112    1.475   -0.490  
O     1.470   -0.876   -0.490  
O    -1.246   -0.876   -0.490  
O     0.112   -0.092    1.727

This file contains five atoms, with the first line indicating this number. 
The following lines specify the atom type and its position in Cartesian coordinates. 
Any comments or blank lines can be added, and they will be ignored by the program.
By filling out these fields and providing a valid .xyz file, you can successfully 
generate CPMD input files for your simulations.
"""
   
    help_collision = f"""
Steps to follow to set up initial velocity in collision molecular dynamics:

1 - Make sure you have the GEOMETRY.xyz file from the CPMD run.
2 - Load the GEOMETRY.xyz file into gqteaWin by pressing open file button. 
3 - Fill in all the entry boxes on the control panel.
4 - Press the "Exec" button to create the new GEOMETRY file with the 
    initial velocity of the attack molecule pointing in the direction 
    of the target position.

There is no direct way to study molecular collisions using the CPMD program. 
Therefore, some additional procedures are necessary to simulate the collisions. 
The following steps should be followed when studying collision simulations:

1 - Generate the system input using GaussView.
2 - Generate the CPMD input using the gqteaWin program.
3 - Run one step of the dynamics using CPMD.
4 - Download the GEOMETRY.xyz file generated by CPMD.
5 - Using the GEOMETRY.xyz file and the gqteaWin program, generate a new GEOMETRY 
    (Collision-GEOMETRY) file with an initial velocity chosen 
    for the attack molecule in the direction 
    of the chosen molecular target.
6 - Delete the old GEOMETRY file and rename Collision-GEOMETRY to GEOMETRY 
    and restart the CPMD calculation with the new GEOMETRY file by using 
    the GEOFILE keyword. This keyword is important because otherwise, the calculation would be 
    restarted from the old velocities stored in the RESTART.1 file.

    Restart cpmd.x with the following keywords:
    RESTART WAVEFUNCTION COORDINATES VELOCITIES GEOFILE LATEST

7 - Therefore, The attack molecule will have the kinetic energy previously stipulated 
    for the collision. 
8 - The trajectory can be visualized using the Visual Molecular Dynamics (VMD) software, 
    and the data analysis can be carried out using the gqteaWin program.
    
"""

    main_window = f"""
                        

                *********************************************
                * MAIN CONTRIBUTORS                         
                * Ademir J. CAMARGO      ajc@ueg.br             
                * Valter H. C. SILVA     fatioleg@ueg.br        
                * Solemar S. OLIVEIRA     solemar@ueg.br         
                * Hamilton B. NAPOLITANO   hamilton@ueg.br        
                * Luciano RIBEIRO         lribeiro@ueg.br        
                * Flávio O. Sanches      flavio.neto@ifg.edu.br 
                *********************************************
                
The gqteaWin provides tools to help analyze the results of the 
molecular dynamics run. The following features and applications 
are now available:

1. Geometry Analysis:
   . Bond length analysis
   . Bond angle analysis
   . Dihedral angle analysis
   . Hydrogen bond analysis
2. Structural
   . Radial distribution function
   . Mean Residence Time
   . Autocorrelation function
   . Solvent box
   . Solvent box with a single solute
   . Mixtures of two solvents and solute

3. CPMD:
    . cpmd inputs
        - Wave function optimization
        - Equilibration
        - Simulation for production
    . Select frames and create Gaussian inputs
    . Surface Hopping inputs
        - Spectra 
        - Excited-state molecular dynamics simulation
    . Molecular collision with initial velocity setup
    . Energy plots

4. Tools:
    . Remote access to a server
    . 3D visualization of the molecule
    . Coordinate transformation

"""
    
    minFframesSelection = f"""This application allows users to extract frames from molecular
    dynamics trajectory files in which the bond length between two selected atoms falls 
    within a user-defined interval. To use it, enter the minimum and maximum bond lengths, 
    specify the indices of the two atoms, and select the trajectory file in XYZ format. 
    After clicking “Extract Frames”, the program analyzes the trajectory and saves all 
    frames that satisfy the bond-length criterion to a new file in the same directory as the input file.

To select frames around the minimum Helmholtz free energy, follow the procedure below:

1. Perform a bond length analysis for the two atoms of interest by navigating to:
   Geometry → Bond Length Analysis
2. From the resulting Helmholtz free energy plot, identify the bond length corresponding to the minimum Helmholtz free energy.
3. Define a bond-length interval around this minimum value. For example, if the minimum 
   free energy occurs at a bond length of 1.5 Å, you may choose a range of 1.5±0.1 Å 
   (based on the standard deviation obtained from the bond length analysis). 
   In this case, the program will select frames with bond lengths between 1.4 Å and 1.6 Å.
4. Enter the indices of the two atoms involved in the bond.
5. Select the XYZ trajectory file containing the frames to be analyzed.
 Click “Extract Frames” to identify and extract all frames whose bond length lies within the interval defined in Step 3.
7. The selected frames will be written to a new file in the same folder as the original trajectory file.
8. In addition, a .txt file will be created in the same directory. This file will 
   contain the averagebond lengths and the stand deviation values associated 
   with each selected frame, and it may be used for further analysis or reference.
9. The extracted frames can be visualized using molecular visualization software using the 3D Molecular Viewer 
   available in gqteaWin or other software such as VMD, PyMOL, etc. 
   This allows you to examine the structural characteristics of the selected frames and gain 
   insights into the molecular behavior at the bond lengths corresponding to the minimum Helmholtz free energy.
10. The frame which interatomic distance closest to the average value is also saved in a separate 
    file named "frame_closest_to_average.xyz".
"""

    help_bond_analysis = f"""BOND LENGTH ANALYSIS MODULE

This python module calculates: 
 - bond lengths (angstroms); 
 - average bond length (angstroms); 
 - bond length distribution function;
 - free energy (kcal/mol) using the probability distribution function.
 - The free energy is calculated using the formula G = -R*T*ln(P(r)), where P(r) is 
   the probability distribution function of the bond lengths, i.e., P(r) is calculated 
   as the number of frames per bin divided by the total number of frames. 
   R is the Boltzmann constant, and T is the simulation temperature in kelvin. 

The analysis is based on the TRAJEC.xyz file from a cpmd run.

PURPOSE AND SIGNIFICANCE
The bond length analysis is used to study the behavior of atoms or molecules in 
a system over time. It provides information on the average bond length, bond length 
distribution function, and free energy, which are important parameters in understanding 
the behavior of the interatomic bonding.

STEPS TO FOLLOW:
1. Load the TRAJEC.xyz file from the cpmd run by pressing on browse button.
3. Fill in the entry boxes, namely:
   3.1 Select two atoms to calculate the bond length.
   3.2 Enter the simulation temperature in kelvin.
   3.3 Specify the number of bins for the bond length distribution function.
   3.4 Click the 'Exec' button to carry out the atomic bond analysis for the 
       two selected atoms

 References:
 - J Mol Model (2011) 17:2159–2168 DOI 10.1007/s00894-010-0939-6
 - Science 275, 817 (1997) DOI: 10.1126/science.275.5301.817
"""

    help_bond_angle = f"""
This python module calculates: 
 - atomic bond angles (degrees); 
 - average bond angle (degree); 
 - bond angle distribution function;
 - bond angle free energy (kcal/mol) using the probability 
   distribution function.
    
STEPS TO FOLLOW:
1. Load the TRAJEC.xyz file from the cpmd run by pressing open 
   file button.
2. Choose the destination folder for your gqteaWin analysis files 
   by clicking the 
   'Save Dir' button.
3. Fill in all the entry boxes on the control panel.
4. Click the 'Exec' button to carry out the atomic bond analysis 
   for the two selected atoms.

The bond angle free energy is calculated using the formula 
G = -R*T*ln(P(r)), where P(r) is the probability distribution 
function of the bond lengths p(r). The P(r) is calculated 
as the number of frames per bin divided by the total number of 
frames. R is the Boltzmann constant, and T is the simulation 
temperature in kelvin.

 References:

 - J Mol Model (2011) 17:2159–2168 DOI 10.1007/s00894-010-0939-6
 - Science 275, 817 (1997) DOI: 10.1126/science.275.5301.817
 """
    help_RDF = f"""
This python module calculates the radial pair distribution function 
(RDF) for the two selected atoms.
    
STEPS TO FOLLOW:
1. Load the TRAJEC.xyz file from the cpmd run by pressing open file 
   button.
2. Choose the destination folder for your gqteaWin analysis files by 
   clicking the 'Save Dir' button.
3. Fill in all the entry boxes on the control panel.
4. Click the 'Exec' button to carry out the atomic bond analysis for 
   the two selected atoms.

The radial pair distribution function (or simply radial distribution 
function - RDF) is a mathematical function that describes the probability 
density of finding a particle at a certain distance from another particle 
in a given system. It is commonly used in materials science, chemistry, 
and physics to analyze the spatial arrangement 
of particles in a system. The radial distribution function is defined as:

                 g(r) = 1 / (4πr²ρ) ∑ᵢ≠ⱼ δ(r - rᵢⱼ)

where r is the distance between two particles, ρ is the number density of 
particles in the system, and δ is the Dirac delta function. The sum is taken 
over all pairs of particles in the system, excluding pairs that include the 
same particle twice.The radial distribution function is often plotted as a 
function of r, and it provides information about the distribution of particles 
around a reference particle. If g(r) is close to unity, the particles are 
distributed randomly. If g(r) is greater than unity, the particles are more 
likely to be found at that distance from the reference particle. If g(r) is 
less than unity, the particles are less likely to be found at that distance. 
The radial distribution function is a powerful tool for characterizing the 
structural properties of a system, and it can provide insights into the nature 
of interactions between particles in the system.

 References:
Fischer-Colbrie, Bienenstock, Fuoss, Marcus. Phys. Rev. B (1988) 38, 12388
Jensen, K. M., Billinge, S. J. (2015). IUCrJ, 2(5), 481-489.
https://chem.libretexts.org/Bookshelves/Biological_Chemistry/Concepts_in_Biophysical_Chemistry_(Tokmakoff)/01%3A_Water_and_Aqueous_Solutions/01%3A_Fluids/1.02%3A_Radial_Distribution_Function#:~:text=g()%3D%CF%81(2,the%20radial%20pair%2Ddistribution%20function.

 """
    help_dihedral_angle = f"""
This python module calculates: 
 - dihedral angles in degrees; 
 - average dihedral angle in degrees; 
 - dihedral angle distribution function;
 - dihedral angle free energy (kcal/mol) using the probability distribution 
   function.
    
STEPS TO FOLLOW:
1. Load the TRAJEC.xyz file from the cpmd run by pressing open file button.
2. Choose the destination folder for your gqteaWin analysis files by clicking the 
   'Save Dir' button.
3. Fill in all the entry boxes on the control panel.
4. Click the 'Exec' button to carry out the atomic bond analysis for the two selected atoms.

The computation of the dihedral angles in degrees was carried out using the the following algorithm:
    a = np.array(coords1)
    b = np.array(coords2)
    c = np.array(coords3)
    d = np.array(coords4)

    ab = b - a
    bc = c - b
    cd = d - c

    normal1 = np.cross(ab, bc)
    normal2 = np.cross(bc, cd)

    x = np.dot(normal1, normal2)
    y = np.dot(np.cross(normal1, bc/np.linalg.norm(bc)), normal2)

    angle_rad = np.arctan2(y, x)
    angle_deg = np.degrees(angle_rad)

The dihedral angle free energy is calculated using the formula G = -R*T*ln(P(r)), where P(r) is 
the probability distribution function of the dihedral angles P(r). The P(r) is calculated 
as the number of frames per bin divided by the total number of frames. 
R is the Boltzmann constant, and T is the simulation temperature in kelvin.

 References:
 - https://en.wikipedia.org/wiki/Atan2
 - https://en.wikipedia.org/wiki/Dihedral_angle 
 - J Mol Model (2011) 17:2159–2168 DOI 10.1007/s00894-010-0939-6
 - Science 275, 817 (1997) DOI: 10.1126/science.275.5301.817
 """

    help_selectFrames = f"""
To select frames from a TRAJEC.xyz trajectory and generate a Gaussian input file,
follow these steps:

1. Make sure the TRAJEC.xyz file from the CPMD run is in the same folder as gqteaWin,
   or specify the full path to the file in the entry box.
2. Choose the destination folder for your gqteaWin analysis files by clicking the 
   'Save Dir' button.
3. Fill in all the entry boxes on the control panel.
4. To create a new trajectory with selected frames only, fill in the first three entries 
   on the left-hand side: start frame, number of frames to skip between two selected 
   frames (second entry), and final frame (third entry).
5. To generate the Gaussian input file, fill in the remaining entry boxes.
6. Select the check box.
7. Click the 'Load File' button to load the TRAJEC.xyz file.
8. Click the 'Exec' button.

A new trajectory, named newTRAJEC.xyz, will be created from the selected frames. 
The Gaussian input file will have the generic name g16_input_x.gjf.
"""    
    
    help_spectraSH = f"""

This tutorial explains how to use this Python module to extract frames from a CPMD run 
and prepare input files for spectra and surface hopping molecular dynamics simulations
using GEOMETRY.xyz and TRAJECTORY from cpmd program. 

STEPS TO FOLLOW:

1. Specify the full path to the GEOMETRY.xyz and TRAJECTORY files from the CPMD run
   by pressing the respectively buttons on the right side of the screen.
2. Fill in all the entry boxes in the control panel.
3. To create the input files for spectra and SH dynamics, fill in the six entries boxes
   on the control panel: start frame, number of frames to skip between two selected 
   frames (second entry), final frame (third entry), number of states, initial state for
   TDDFT dynamics, charge on the system, cell parameters, cutoff energy, and dual. 
   The place-holder in the entries are self explained. 
7. Click the 'Exec' button to generate the input files.
8. Submit the prefix-SPECTRA-xxx.inp file created by gqteaWin to be running on cpmd program:
   8.1 -> mpirun -np 20 cpmd.x prefix-SPECTRA-xxx.inp > prefix-SPECTRA-xxx.out &
   8.2 -> After the calculation has finished, replace the GEOMETRY file by the one created
          by gqteaWin (prefix-GEOMETRY-xxx). Do not forget to change its name: 
          prefix-GEOMETRY-xxx --> GEOMETRY
    8.3 -> Restart the cpmd calculation with with prefix-SHTDDFT-XXX.inp input.

Notes: This kind of calculation does not work with Car-Parrinello molecular dynamics
       or Vanderbilt ultrasoft pseudopotentials. As a result, all the calculation 
       carried out here is done using BOMD and Martin-Truller pseudopotentials.

The inputs are used by the CPMD program to simulate the behavior of molecules 
as they undergo "surface hopping" between different energy states using a version of the 
Tully's trajectory surface hopping algorithm that has been adapted to use time-dependent 
density functional theory (TDDFT) to calculate the forces on the excited state surfaces.
The simulation works by starting with an initial configuration of the molecule and then 
computing the forces on the different excited state surfaces using TDDFT. The program then 
selects which surface to use as the "running surface" based on the current energy and other
factors. If a surface hop occurs, the running surface changes and the simulation continues 
on the new surface.
To perform this simulation, you need to declare a sufficiently large number of excited 
states using the "STATES" keyword in the TDDFT section. You also need to specify the 
initial running surface using the "FORCE STATE" keyword. After the simulation, a series 
of output files are generated that contain information on the state amplitudes, coupling 
strength, and transition probabilities between different states. These files are needed to 
restart the simulation from a previous point.
To gather statistics, you need to run the simulation multiple times using different initial
coordinates and velocities.

References:
1. E. Tapavicza, I. Tavernelli, and U. Rothlisberger, Phys. Rev. Lett. 98, 023001 (2007).
2. I. Tavernelli, E. Tapavicza, and U. Rothlisberger, J. Mol. Struct. : THEOCHEM 914, 22 (2009).

"""

    help_HB = f"""
Hydrogen bonds analysis tutorial

The geometrical parameters that characterize a hydrogen bond are bond length and bond 
angle.

1. Bond Length: Hydrogen bond length refers to the distance between the hydrogen atom 
and the atom it is bonded to. This length is generally longer than a typical covalent 
bond but shorter than van der Waals interactions. The exact length can vary depending 
on the atoms involved and the surrounding environment, but it typically falls in the 
range of 1.5 to 2.2 Angstroms.

2. Bond Angle: The bond angle in a hydrogen bond is the angle between the donor atom, 
the hydrogen atom, and the acceptor atom. In an ideal hydrogen bond, this would be a 
straight line, i.e., 180 degrees, but in reality, there can be some deviation from this. 
The closer the hydrogen bond angle is to 180 degrees, the stronger the bond tends to be, 
because this arrangement allows for optimal overlap between the involved orbitals.

It's important to note that these parameters are not fixed; they can be influenced by 
many factors, including the specific atoms involved in the bond, the surrounding chemical 
environment, temperature, pressure, and so forth. Understanding these parameters can give 
us insight into the strength and properties of the hydrogen bond and how it might behave 
under different conditions.

STEPS TO FOLLOW:

1. Specify the full path to the GEOMETRY.xyz and TRAJECTORY files from the CPMD run
   by pressing the respectively buttons on the right side of the screen.
2. Fill in all the entry boxes in the control panel.
    2.1 Specify the bond bin width and angle bin width
    2.2 Specify the simulation time step in atomic units
    2.3 Specify the amount of frames sipet between two frames collected
    2.4 Chose a directory to save the calculation results by pressing the 'Save dir' buttom
    2.5 Specify the temperature that MD was run
    2.6 Specify the maximum length between H--A, that is, the max r to be considered a 
        hydrogen bond.
    2.7 Specify atom labels involved in the hydrogen bond, such as donnor, H, and Acceptor

8. Press 'Exec' button to carried out the hydrogen bond analysis

"""
    help_mrt = f"""
MRT - Mean Residence Time Analysis
The mean residence time in molecular dynamics simulation refers 
to the average duration that a molecule or atom spends within a 
particular region or state during the course of the simulation. 
It represents the typical amount of time that a molecule remains 
in a specific location or state before transitioning  to another 
region or  undergoing  a  particular event. This metric provides 
insights into the dynamics and kinetics of the system under study, 
allowing  researchers  to  analyze and understand the behavior of 
molecules and their interactions within a simulated environment.

STEPS TO FOLLOW TO CALCULATE MRT

1. Specify the full path to the TRAJEC.xyz files from a CPMD run
   by pressing the 'Open file' button on  the  right side on the 
   control panel.
2. Fill in all the entry boxes in the control panel.
    2.1 Specify the shell's inner radius in Å
    2.2 Specify the shell's outer radius in Å
    2.3 Specify the simulation time step in Å
    2.4 Specify  the  folder  to  save the mrt analysis results by 
        pressing the 'Save dir' buttom at the right panel
    2.5 Define the frame sampling rate, which refers to the  number 
        of frames omitted between two successive frames collected.
    2.6 Specify the labels (indexes) of the atoms to be excluded 
        from the mean residence time (MRT) analysis. Ensure that 
        the labels are separated by whitespace instead of commas.
    2.7 Specify atom labels (index) at the shell's center
    2.8 Specify the element symbol to investigate the mean residence 
        time within the shell.

8. Press 'Exec' button to carried out the MRT

The mean residence time and the coordination number are related in
thecontext  of  molecular  dynamics  simulations. The coordination 
number refers to the number of  neighboring  atoms  or  molecules 
that are in direct contact with a central atom or molecule. 
It provides information about the local environment and the extent 
of interactions surrounding a specific atom. In molecular dynamics 
simulations, the mean residence time can be influenced by the 
coordination number. Generally, a higher coordination number 
indicates a greater number of surrounding atoms or molecules, 
implying stronger interactions and a potentially longer mean 
residence time. This is because a higher coordination number implies
a more stable or constrained environment for the central atom or 
molecule, which can result in longer residence times within that 
specific region.However, it's  important  to  note  that  the 
relationship between mean residence time and coordination number can 
vary depending on the specific system being studied and the nature of 
the  interactions involved. Other factors, such as temperature, 
pressure, and  potential energy landscapes, can also impact the mean 
residence time independently of the coordination number. Therefore, 
it is  essential  to consider  multiple  factors  and analyze the 
system holistically when examining the  relationship between mean  
residence  time  and  coordination  number in  molecular  dynamics 
simulations.
"""

    help_vaf = f"""
TUTORIAL

VAF - Velocity Autocorrelation Function 

Velocity autocorrelation function (VAF), also known as 
the velocity autocorrelation coefficient, is a measure 
used in statistical mechanics and physics to quantify 
the correlation between the velocity of a particle at 
one time and its velocity at a later time. 
Mathematically, for a system of particles, the VAF 
for a single particle "i" can be expressed as:

C(t) = <v_i(0),v_i(t)>

Here, < > denotes the ensemble average, v_i(0) is the 
velocity of the particle "i" at time "0", and v_i(t) 
is the velocity at some later time "t". VAF represents 
the memory of the system or how much the system 
"remembers" its past states. A rapid decay of VAF means 
the system has little memory of past velocity states, 
while a slower decay indicates a longer memory.

Importance of VAF:

1. Understanding System Dynamics: The VAF gives us 
    insight into the dynamics of the system. By examining 
    how the VAF decays with time, we can get an idea of 
    how the velocity of a particle changes over time and 
    the nature of the interactions within the system.

2. Determining Diffusion Coefficients: In many-body 
    systems, the VAF is directly related to the 
    self-diffusion coefficient through the Green-Kubo 
    relations. The integral of the VAF over all time 
    gives the diffusion coefficient, which is a crucial 
    quantity in understanding transport phenomena.

3. Testing Molecular Dynamics Simulations: In molecular 
    dynamics simulations, the VAF is an important test for 
    the correctness of the simulations. The VAF can be 
    compared between experimental data and simulations to 
    assess the quality of the models and force fields used.

4. Spectral Properties: The Fourier transform of the 
    VAF gives the power spectral density of the velocity, 
    which can provide important information about the 
    frequency components of a particle's motion and the 
    underlying dynamical processes.

5. Study of Liquids and Gases: VAF is also crucial 
    in the study of the dynamics of liquids and gases, 
    particularly in understanding viscosity and thermal 
    conductivity. The Green-Kubo relations also link the 
    autocorrelation functions of the stress tensor in a fluid 
    to its shear viscosity and thermal conductivity.

STEPS TO FOLLOW TO CALCULATE MRT

1. Specify the full path to the TRAJECTORY files from 
    a CPMD run by pressing the 'Open' button on  
    the right side on the control panel.

2. Fill in all the entry boxes in the control panel.
    2.1 Specify the start frame
    2.2 Specify the amount of frames to use in the VAF 
        calculation.

3. Press 'Exec' button to carried out the MRT

"""
    help_solBox = f"""A solvation box, often used in molecular dynamics simulations, is a box that contains a solute (like a protein) and a solvent (like water). The purpose of the solvation box is to mimic a realistic biological environment for the solute.

There are two main things that a solvation box can do¹²:
1. Generate a box of solvent. This is done when you have a structure file with a box, but without atoms.
2. Solvate a solute configuration, like a protein, in a bath of solvent molecules. This is done when you have a solute and a solvent.

The box of solute is built by stacking the coordinates read from the coordinate file. These coordinates should be equilibrated in periodic boundary conditions to ensure a good alignment of molecules on the stacking interfaces¹².

Solvent molecules are removed from the box where the distance between any atom of the solute molecule and any atom of the solvent molecule is less than the sum of the scaled van der Waals radii of both atoms¹².

In some cases, a layer of water of a specified thickness can be placed around the solute¹². This is often done to simulate the hydration shell that surrounds proteins in a biological environment.

It's important to note that the usefulness of the van der Waals radii depends on the atom names, and thus varies widely with the force field¹²..

Source: Conversation with Bing, 19/12/2023
(1) gmx solvate — GROMACS 2018 documentation. https://manual.gromacs.org/documentation/2018/onlinehelp/gmx-solvate.html.
(2) gmx solvate - GROMACS 2023.3 documentation. https://manual.gromacs.org/current/onlinehelp/gmx-solvate.html.
(3) Solvating the Protein - University of Illinois Urbana-Champaign. https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-html/node7.html.
(4) en.wikipedia.org. https://en.wikipedia.org/wiki/Solvation.
"""
    mixture_solvent_box = f""" Mixture Solvent Box Builder: A Step-by-Step Guide
The Solvent Box Builder allows you to generate a solvent box with a mixture of two different solvents. 
This feature enables the simulation of complex environments where multiple solvent types coexist, such 
as aqueous-organic mixtures or binary solvent systems. Users can specify the composition of each solvent,
ensuring the correct ratio in the final solvent box. With this feature, you can:Select two different 
solvent molecules (Solvent A and Solvent B) in XYZ format. Specify the composition of the mixture 
(e.g., 40% Solvent A and 60% Solvent B). Define minimum distances between solvent molecules and the 
solute at the center. Ensure a specified minimum distance between solvent molecules to avoid overlap.
By combining two solvents in precise proportions, this feature expands the application's 
utility for more realistic molecular simulations.
"""

    shared_wall_box = f""" Shared-Wall Double Solvent Box Builder: A Step-by-Step Guide
The Shared-Wall Double Solvent Box Builder creates two adjacent solvent regions stacked along
the z axis. Box 1 and Box 2 share one common face, but no wall atoms or boundary markers are
written to the final XYZ file. Each region can use its own pure solvent, binary solvent
mixture, or binary solvent mixture with one centered solute molecule. The two regions share
the same a and b lattice dimensions and use separate c heights. Periodic minimum-image clash
detection is applied using the full combined box.
"""
    
    single_solute_solvent_box = f""" Mixture Solvent Box Builder: A Step-by-Step Guide
The Solvent Box Builder allows you to generate a solvent box with just one solute molecule simulating
the infinit dilution. This feature enables the simulation of complex environments where the solute
is solvated by the solevent molecules. The users can specify any kind of solvent and solute. The user
should define a minimum distances between solvent molecules and the solute and solvent distance, avoiding
overlap between the solute and solvent molecules.
"""

    solvent_box = f""" Solvent Box Builder: A Step-by-Step Guide
The solvent Box Builder allows you to generate a solvent box with a single solvent. 
This feature enables the simulation of complex environments, such as aqueous or organic solvent systems. 
With this feature, you can select a solvent molecule in XYZ format. Define minimum distances between 
solvent molecules. Ensure a specified minimum distance between solvent molecules to avoid overlap.
"""

    help_plots = """
CPMD Energy File Plotting ApplicationWelcome to the CPMD Energy File Plotting 
Application. This application allows you to read CPMD energy files, extract 
energy data, and visualize it through various types of plots, such as fictitious
and ionic kinetic energy, temperature, Kohn-Sham energy, total energy, and CPU time.

 "CPMD ENERGY FILE PLOT".

1. User Interface Overview

The application is divided into several sections:Plot Selection: A series of switches 
where you can choose the type of plot to generate. File Selection: A section for 
selecting the CPMD energy file. Time Step and Unit Selection: Inputs for specifying 
the simulation time step and selecting the x-axis unit. Action Buttons: Buttons 
to execute the plot generation or to close the application

2. Steps to Generate Plots
Step 1: Select a CPMD Energy File

Click the "Browse" Button:
This opens a file dialog that allows you to locate and select the CPMD energy file 
on your computer. The file format should be .txt or similar, containing data 
with 8 columns.

View the Selected File:

Once the file is selected, the file path will be shown in the "Select cpmd 
ENERGY file" input box.Verify File Content: If the file is in a valid format, 
the number of steps will be displayed in the multi-line text area, along 
with a confirmation that the ENERGY file has a valid format.

Step 2: Set the Simulation Time Step
Enter the Time Step:

In the "Simulation Time Step" input box, enter the time step for your 
simulation. This value is used for plotting time-related data, like 
when the x-axis unit is in femtoseconds (fs) or picoseconds (ps).

Select the X-Axis Unit:

Choose the unit for the x-axis from the drop-down menu. Available options 
are: Steps: Plots the x-axis based on simulation steps.fs (Femtoseconds): 
Plots time in femtoseconds, based on the simulation time step.
ps (Picoseconds): Plots time in picoseconds.

Step 3: Choose the Plot Type

Select the Type of Plot:

Turn on the switch for the plot type you want to generate. The following options 
are available: Fictitious and Ionic Kinetic Energy Plot: Plots fictitious 
and ionic kinetic energies from the CPMD simulation.Temperature Plot: 
Plots the simulation temperature in Kelvin. Kohn-Sham Energy Plot: 
Plots the Kohn-Sham potential energy in Hartrees (Ha).Kohn-Sham and Ionic 
Kinetic Energy Plot: Plots both Kohn-Sham and ionic kinetic energies.
Total Energy Plot: Plots the total energy (fictitious, Kohn-Sham, and 
ionic kinetic energies). CPU Time Plot: Plots CPU time spent per simulation step.
You Can Select Multiple Plots:

You can choose one or more plot types by enabling the respective switches. 
The application will generate plots for all the enabled options.

Step 4: Plot the Data
Click the "Plot" Button:

Once the file is selected, time step is set, and plot types are chosen, 
click the "Plot" button. The application will process the data and 
generate the requested plots.

View the Plot:

A new window will open displaying the plot. You can resize or close 
the window as needed.
Step 5: Close the Application
Click the "Close" Button:
To exit the application, simply click the "Close" button.
4. Plot Types Explained
Here is a description of the various plot types you can generate:

Fictitious and Ionic Kinetic Energy:

Displays the fictitious energy of the simulation (energy in fictitious 
particles), and the difference between ionic kinetic energy and 
the total energy.

Temperature Plot:

Displays the temperature of the simulation in Kelvin, which is derived 
from the kinetic energy of the particles.
Kohn-Sham Energy: Shows the Kohn-Sham potential energy, an essential 
quantity in Density Functional Theory (DFT) used to describe the 
electron interactions.
Kohn-Sham and Ionic Kinetic Energy:

Displays both the Kohn-Sham energy and the ionic kinetic energy in one 
plot, giving insight into how these energies 
evolve during the simulation.
Total Energy:

Combines fictitious, Kohn-Sham, and ionic kinetic energies into a single 
plot, providing an overview of the total energy during the simulation.

CPU Time:

Displays the time spent by the CPU per simulation step, helping to evaluate 
computational performance over time.

5. Handling Errors

File Format Issues:

If the selected file does not conform to the expected format (8 columns), 
an error will be displayed in the message area. Ensure that the file has 
valid data before attempting to plot.

Time Step Errors:

If the simulation time step is missing or invalid, an error message will 
appear. Double-check your input to ensure it’s a valid numerical value.

No Data Errors:

If no data is available or if the data is not loaded correctly, the 
application will display an error, prompting you to check your file 
and settings.

6. Example Workflow

Select File: Click "Browse" to select your energy.txt file.Set Time Step: 
Enter the time step, e.g., 0.5. Select X-Axis Unit: Choose "fs" for 
femtoseconds.

Select Plot Types: Turn on "Fictitious and Ionic

Kinetic Energy" and "Temperature Plot."

Generate Plot: Click "Plot." Two plot windows will open, each displaying 
the selected plot.

7. Additional Features

Help Section:

In the text area at the bottom, some help information from the HelpGqteaWin 
class is displayed. This can guide you through specific functions or 
additional tips related to CPMD energy analysis.

8. Troubleshooting

Plot Doesn’t Appear: Ensure that the energy file is in the correct format a
nd that you've entered a valid time step and selected a plot type.CPU Time 
Appears Unusual: The CPU time plot reflects the time taken per step and 
can vary based on the complexity of the simulation and hardware capabilities. 
By following this manual, you should be able to use the CPMD Energy File 
Plotting Application to visualize data from CPMD energy files and gain
valuable insights into the energy dynamics of your molecular simulations.

"""

    Plotting_Options = """
Choose from multiple plot types, including:
    - Fictitious and Ionic Kinetic Energy
    - Temperature
    - Potential Kohn-Sham Energy
    - Kohn-Sham Plus Ionic Kinetic Energy
    - Total Energy
    - CPU Time by Step
    - Flexible X-Axis Units: Select the unit for the x-axis from steps, 
      femtoseconds (fs), or picoseconds (ps).
Interactive GUI: User-friendly interface with options to browse files, select plots, and view results.
For more information on each plot type and how to use the application, 
refer to the Help section in the application manual.

"""

    sh_geom_analyzer = """
   This program is designed to process surface hopping (SH) simulation data by separating 
trajectory frames according to electronic states. It features a graphical user interface built 
with Python and the Toga framework. When you open the program, you will be prompted to input 
the total number of states present in your simulation. This number should reflect all possible 
electronic states, typically starting from zero.
    After entering the number of states, you should click on the "Browse" button to select the 
root directory that contains all the simulation subfolders. Each subfolder must include two 
essential files: SH_STATE.dat, which contains the time evolution of electronic states, and 
TRAJEC.xyz, which contains the atomic trajectories for each frame. The program scans all 
subfolders inside the selected root directory.
    Once the directory is selected and the number of states is entered, you can proceed by 
clicking the "Extract Frames" button. The program will start reading each SH_STATE.dat file to 
determine the state associated with each frame. Simultaneously, it will read the TRAJEC.xyz file 
to extract the atomic configurations for each frame. For each frame, the program checks the 
associated state and writes that frame into a separate file named stateX.xyz, where X is the 
state number. These files are generated within each subfolder and also consolidated in the root 
folder.
    At the end of the process, the program calculates the average occupation of each electronic 
state across all simulations. It generates an output file called sh_avg_perc.dat, which 
summarizes how frequently each state appears as a percentage. This file is saved in the root 
directory. The content of sh_avg_perc.dat lists each state number, the average count of 
frames in that state, and the corresponding percentage.
   If a required file is missing in a subfolder, the program will display a warning message in 
the output window and skip that subfolder. The user is notified in case the SH_STATE.dat or 
TRAJEC.xyz file is not found. The text window in the interface displays all the progress logs, 
including which directories are being processed, whether files are missing, and whether the 
process was successful. Once the processing is complete, the user can simply close the 
application by clicking the Close button.
    The program creates stateX.xyz files inside each subfolder because it processes each 
trajectory file individually. As it reads the TRAJEC.xyz file from a subfolder, it immediately 
classifies each frame based on its associated electronic state and writes the frames into 
separate files named according to the state, such as state0.xyz, state1.xyz, and so on. 
These output files are saved directly in the same subfolder where the input files are located.
    This approach is useful because it allows the user to check how each individual simulation 
contributes to the distribution of frames across the electronic states. It makes it easier 
to verify and debug the results of each specific simulation run before looking at the 
combined data.
    However, this step is not strictly necessary if the user is only interested in the overall 
combined trajectories separated by states. The program already performs a second step where it 
consolidates all the stateX.xyz files from every subfolder into a single file for each state 
in the root directory. These consolidated files contain the full trajectory information grouped 
by state across all subfolders.
"""
