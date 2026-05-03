# gqteaWinToga

gqteaWinToga is a desktop molecular analysis toolkit developed by the gQTEA group
for preparing, inspecting, converting, and analyzing molecular dynamics and
electronic-structure simulation data. The application provides a graphical
interface built with [Toga](https://toga.readthedocs.io/) and collects several
specialized tools for CPMD, Quantum ESPRESSO `cp.x`, ORCA, surface hopping, and
trajectory analysis workflows.

The current application title identifies the toolkit as **gQTEA-0.3.1 Molecular
Analysis Toolkit**.

## Capabilities

The main window organizes the tools into five tabs.

### Geometry

- Bond length analysis
- Bond angle analysis
- Dihedral angle analysis
- Hydrogen bond analysis

These modules work with molecular trajectory files, commonly `TRAJEC.xyz`, and
can calculate averages, distributions, and free-energy profiles from probability
distributions.

### Inputs

- CPMD input generation
- Surface hopping input builder
- Molecular collision setup
- Quantum ESPRESSO `cp.x` input builder
- ORCA input builder
- gqteaMD input builder
- Vanderbilt `runatom.x` input builder

### Structural Analysis

- Radial distribution function
- Mean residence time
- Legacy mean residence time workflow
- Velocity autocorrelation function
- Single-solute solvent box builder
- Mixed-solvent box builder

### Thermodynamics

- Classical rate constant calculation

### General Tools

- 3D molecular viewer
- CPMD energy plotting
- Molecular axis alignment
- Frame selection
- Frame selection by interatomic distance range
- CPMD input to XYZ conversion
- Surface hopping geometry analyzer
- `*.pos` to `trajec.xyz` conversion
- Force conversion from Quantum ESPRESSO `cp.x` `*.for` files
- Coordinate conversion
- XYZ Angstrom-to-Bohr conversion
- Unit conversion

## Repository Layout

This repository is currently organized as a flat Python source tree. The main
launcher is:

```text
gqteaWinToga.py
```

Each major tool is implemented in its own module, for example:

- `bond.py`, `bondAngle.py`, `dihedralAngle.py`, and `hbond.py` for geometry
  analysis
- `cpmdInput.py`, `cpx_input_builder.py`, `orca_input_builder.py`, and
  `gqteaMDinptBuilder.py` for input generation
- `molecularViewer.py` for OpenGL-based molecule visualization
- `plotter.py` and `displayPlots.py` for plotting
- `help.py` for application help text and shared atomic data

## Requirements

- Python 3.10 or newer is recommended.
- Windows is the primary target when using the Toga WinForms backend.
- A working graphical desktop session is required.

Python package dependencies used by the application include:

- `toga-winforms`
- `numpy`
- `scipy`
- `matplotlib`
- `glfw`
- `PyOpenGL`
- `PyOpenGL-accelerate`
- `PyMuPDF`

The code also uses Python standard-library modules such as `asyncio`, `csv`,
`dataclasses`, `decimal`, `html`, `http.server`, `math`, `os`, `pathlib`,
`random`, `re`, `statistics`, `subprocess`, `tempfile`, `threading`, and
`typing`.

## Installation

Clone the repository and create a virtual environment:

```powershell
git clone https://github.com/AdemirCamargo57/gqteaWinToga.git
cd gqteaWinToga
python -m venv venv
.\venv\Scripts\Activate.ps1
```

Install the runtime dependencies:

```powershell
python -m pip install --upgrade pip
python -m pip install toga-winforms numpy scipy matplotlib glfw PyOpenGL PyOpenGL-accelerate PyMuPDF
```

If you are already working inside this repository's `venv/src` directory, run the
install command from the active virtual environment.

## Running the Application

From the source directory, start the GUI with:

```powershell
python gqteaWinToga.py
```

The application opens a desktop window titled:

```text
gQTEA - grupo de Quimica Teorica e Estrutural de Anapolis
```

Select a tab, choose the desired tool, fill in the required file paths and
parameters, and run the calculation or conversion from the tool window.

## Usage Examples

### Analyze a Bond Length Distribution

1. Start the application with `python gqteaWinToga.py`.
2. Open **Geometry**.
3. Select **Bond length analysis**.
4. Load a CPMD `TRAJEC.xyz` trajectory file.
5. Enter the two atom indices, simulation temperature, and histogram bin count.
6. Run the analysis to generate bond-length statistics, distributions, and
   free-energy data.

### Plot CPMD Energy Data

1. Open **Tools**.
2. Select **Energy plots**.
3. Browse to a CPMD energy file.
4. Enter the simulation time step.
5. Select one or more plot types, such as temperature, Kohn-Sham energy, total
   energy, or CPU time.
6. Generate the plots.

### Build a CPMD Input

1. Open **Inputs**.
2. Select **CPMD Inputs**.
3. Provide a prefix, total charge, periodic box parameters, energy cutoff, dual
   for rho expansion, cell symmetry, and an XYZ starting geometry.
4. Generate the CPMD input files for wave-function optimization, equilibration,
   or production simulation workflows.

### Extract Surface Hopping Frames by Electronic State

1. Open **Tools**.
2. Select **SH Geometry Analyzer**.
3. Enter the number of electronic states.
4. Select a root directory containing simulation subfolders with `SH_STATE.dat`
   and `TRAJEC.xyz`.
5. Run the extractor to create state-separated `stateX.xyz` trajectory files and
   the summary file `sh_avg_perc.dat`.

## Supported Input and Output Files

Common files used by the toolkit include:

- `.xyz` molecular geometry and trajectory files
- `TRAJEC.xyz` CPMD trajectory files
- CPMD `GEOMETRY.xyz` files
- CPMD energy files
- Surface hopping `SH_STATE.dat` files
- Quantum ESPRESSO `cp.x` input and force files
- ORCA input files
- Vanderbilt `runatom.x` inputs

Generated outputs depend on the selected tool and may include analysis tables,
selected trajectory frames, converted geometry files, generated input decks, and
plots.

## Packaging for Windows

The project can be packaged as a Windows executable with
[`auto-py-to-exe`](https://pypi.org/project/auto-py-to-exe/) and installed with
Inno Setup.

Basic packaging workflow:

```powershell
cd <project-venv>
python -m pip install auto-py-to-exe
auto-py-to-exe
```

When packaging the application, ensure the required runtime packages are included:

- `toga-winforms`
- `matplotlib`
- `numpy`
- `scipy`
- `glfw`
- `PyOpenGL`
- `PyOpenGL-accelerate`
- `PyMuPDF`

`auto-py-to-exe` may not automatically find the GLFW DLL. If packaging fails or
the molecular viewer cannot start, manually include the DLL from a path similar
to:

```text
venv\Lib\site-packages\glfw\glfw3.dll
```

## Development Notes

- The application is launched from `gqteaWinToga.py`.
- Most user-facing modules expose a `*UI` class imported by the launcher.
- Keep new tools self-contained in their own module when possible, then register
  the UI class in the appropriate tab in `gqteaWinToga.py`.
- Prefer adding clear help text in `help.py` for new analysis workflows.

## Contributors

gqteaWinToga core development team:

- Ademir J. Camargo - ajc@ueg.br
- Valter H. C. Silva - fatioleg@ueg.br
- Solemar S. Oliveira - solemar@ueg.br
- Hamilton B. Napolitano - hamilton@ueg.br
- Luciano Ribeiro - lribeiro@ueg.br
- Flavio O. Sanches - flavio.neto@ifg.edu.br

## License

No license file is currently included in this repository. Add a `LICENSE` file
before distributing the project publicly so users know how they may use, modify,
and redistribute the code.
