# gqteaWinToga User Manual

## Index

- [Overview](#overview)
- [Installation and Requirements](#installation-and-requirements)
  - [Requirements](#requirements)
    - [Python packages](#python-packages)
    - [System requirements that pip cannot install](#system-requirements-that-pip-cannot-install)
    - [Optional packages](#optional-packages)
    - [External simulation programs](#external-simulation-programs)
  - [Install from Source](#install-from-source)
  - [Start the Program](#start-the-program)
- [Input File Basics](#input-file-basics)
  - [XYZ Files](#xyz-files)
  - [Atom Labels](#atom-labels)
  - [Common Units](#common-units)
- [Basic Workflows](#basic-workflows)
  - [Analyze a Bond Length](#analyze-a-bond-length)
  - [Run a QTAIM Topology Analysis](#run-a-qtaim-topology-analysis)
  - [Generate Natural Transition Orbitals](#generate-natural-transition-orbitals)
  - [Build CPMD Input Files](#build-cpmd-input-files)
  - [Select Frames and Generate Gaussian Inputs](#select-frames-and-generate-gaussian-inputs)
  - [Plot CPMD Energy Data](#plot-cpmd-energy-data)
- [Features and Commands](#features-and-commands)
  - [Main Window](#main-window)
  - [Geometry Tools](#geometry-tools)
    - [Bond Length Analysis](#bond-length-analysis)
    - [All Bond Distance Analysis](#all-bond-distance-analysis)
    - [Bond Angle Analysis](#bond-angle-analysis)
    - [All Bond Angle Analysis](#all-bond-angle-analysis)
    - [Dihedral Angle Analysis](#dihedral-angle-analysis)
    - [All Dihedral Angle Analysis](#all-dihedral-angle-analysis)
    - [Hydrogen Bond Analysis](#hydrogen-bond-analysis)
    - [Comparison of Molecular Geometric Parameters](#comparison-of-molecular-geometric-parameters)
  - [Input Builders](#input-builders)
    - [CPMD Inputs](#cpmd-inputs)
    - [Surface Hopping Input Builder](#surface-hopping-input-builder)
    - [Collision Input](#collision-input)
    - [Quantum ESPRESSO cp.x Input Builder](#quantum-espresso-cpx-input-builder)
    - [ORCA Input Builder](#orca-input-builder)
    - [gqteaMD Input Builder](#gqteamd-input-builder)
    - [Gaussian Fukui Input Builder](#gaussian-fukui-input-builder)
    - [Vanderbilt runatom.x Input Builder](#vanderbilt-runatomx-input-builder)
  - [Structural Tools](#structural-tools)
    - [Radial Distribution Function](#radial-distribution-function)
    - [Mean Residence Time](#mean-residence-time)
    - [Legacy Mean Residence Time](#legacy-mean-residence-time)
    - [Autocorrelation Function](#autocorrelation-function)
    - [Single Solute Solvent Box](#single-solute-solvent-box)
      - [System properties](#system-properties)
      - [Pre-optimization](#pre-optimization)
    - [Mixture of Two Solvent Box](#mixture-of-two-solvent-box)
    - [Shared-Wall Double Solvent Box](#shared-wall-double-solvent-box)
  - [Thermo Tools](#thermo-tools)
    - [Classical Rate Constant](#classical-rate-constant)
  - [General Tools](#general-tools)
    - [3D Molecular Viewer](#3d-molecular-viewer)
      - [Cube tab](#cube-tab)
    - [Electrostatic Surface Potential (ESP)](#electrostatic-surface-potential-esp)
    - [QTAIM Topology Analysis](#qtaim-topology-analysis)
    - [Natural Transition Orbitals (NTO)](#natural-transition-orbitals-nto)
      - [UV-Vis spectrum (Plot Spectrum)](#uv-vis-spectrum-plot-spectrum)
    - [Analysis of Fukui Functions](#analysis-of-fukui-functions)
      - [Molecular Design](#molecular-design)
        - [Opening and using the two windows](#opening-and-using-the-two-windows)
        - [Creating and modifying bonds](#creating-and-modifying-bonds)
        - [Running the pre-optimization](#running-the-pre-optimization)
        - [Importing and saving structure files](#importing-and-saving-structure-files)
        - [Selecting and moving one molecule](#selecting-and-moving-one-molecule)
        - [Measuring bond lengths, angles and dihedrals](#measuring-bond-lengths-angles-and-dihedrals)
        - [Setting a bond length, angle or dihedral](#setting-a-bond-length-angle-or-dihedral)
        - [Freezing parameters during pre-optimization](#freezing-parameters-during-pre-optimization)
        - [Adding hydrogens automatically](#adding-hydrogens-automatically)
        - [Replacing an atom's element](#replacing-an-atoms-element)
        - [Showing atom numbers](#showing-atom-numbers)
        - [Canvas background](#canvas-background)
        - [Bringing a molecule back from the 3D viewer](#bringing-a-molecule-back-from-the-3d-viewer)
        - [Mouse and keyboard reference](#mouse-and-keyboard-reference)
        - [Validation messages](#validation-messages)
    - [Energy Plots](#energy-plots)
    - [Molecular Axis Alignment](#molecular-axis-alignment)
    - [Select Frames](#select-frames)
    - [Frame Selection by Interatomic Distance Range](#frame-selection-by-interatomic-distance-range)
    - [CPMD Input to XYZ Converter](#cpmd-input-to-xyz-converter)
    - [SH Geometry Analyzer](#sh-geometry-analyzer)
    - [Convert cp.x .pos File to trajec.xyz](#convert-cpx-pos-file-to-trajecxyz)
    - [Compute Forces from cp.x .for File](#compute-forces-from-cpx-for-file)
    - [Coordinate Converter](#coordinate-converter)
    - [Unit Converter](#unit-converter)
- [Configuration](#configuration)
- [Practical Examples](#practical-examples)
  - [Example 1: Find a Bond-Length Free-Energy Minimum](#example-1-find-a-bond-length-free-energy-minimum)
  - [Example 2: Prepare a CPMD Collision Restart](#example-2-prepare-a-cpmd-collision-restart)
  - [Example 3: Group Surface Hopping Frames by State](#example-3-group-surface-hopping-frames-by-state)
  - [Example 4: Convert a cp.x Trajectory](#example-4-convert-a-cpx-trajectory)
  - [Example 5: Analyze All Solute Bond Distances](#example-5-analyze-all-solute-bond-distances)
  - [Example 6: Analyze All Solute Bond Angles](#example-6-analyze-all-solute-bond-angles)
  - [Example 7: Analyze All Solute Dihedral Angles](#example-7-analyze-all-solute-dihedral-angles)
  - [Example 8: Read a Fukui Function Figure (phenol, f⁺)](#example-8-read-a-fukui-function-figure-phenol-f)
  - [Example 9: QTAIM Analysis of a Hydrogen Bond (water dimer)](#example-9-qtaim-analysis-of-a-hydrogen-bond-water-dimer)
  - [Example 10: Natural Transition Orbitals of Uracil S3](#example-10-natural-transition-orbitals-of-uracil-s3)
- [Formulas and Methods (All Bond / Angle / Dihedral tools)](#formulas-and-methods-all-bond--angle--dihedral-tools)
- [Formulas and Methods (Bond Length / Bond Angle / Dihedral Angle tools)](#formulas-and-methods-bond-length--bond-angle--dihedral-angle-tools)
- [Formulas and Methods (Mean Residence Time tool)](#formulas-and-methods-mean-residence-time-tool)
- [Formulas and Methods (Electrostatic Surface Potential tool)](#formulas-and-methods-electrostatic-surface-potential-tool)
- [Formulas and Methods (Fukui Function tool)](#formulas-and-methods-fukui-function-tool)
- [Theoretical Background: the Quantum Theory of Atoms in Molecules](#theoretical-background-the-quantum-theory-of-atoms-in-molecules)
- [Formulas and Methods (QTAIM tab)](#formulas-and-methods-qtaim-tab)
- [Formulas and Methods (NTO tab)](#formulas-and-methods-nto-tab)
- [Formulas and Methods (Solvation box builders)](#formulas-and-methods-solvation-box-builders)
- [Troubleshooting](#troubleshooting)
  - ["No file was selected!"](#no-file-was-selected)
  - ["Please input a valid value for ..."](#please-input-a-valid-value-for-)
  - ["Invalid format for ... Please input exactly two atom labels."](#invalid-format-for--please-input-exactly-two-atom-labels)
  - ["Failed to open file" or "Failed to read file"](#failed-to-open-file-or-failed-to-read-file)
  - ["The file is empty!"](#the-file-is-empty)
  - ["Invalid line format in TRAJECTORY file"](#invalid-line-format-in-trajectory-file)
  - ["The stop frame must be less than the total number of frames"](#the-stop-frame-must-be-less-than-the-total-number-of-frames)
  - [ENERGY file format errors](#energy-file-format-errors)
  - [Solvent box insertion warnings](#solvent-box-insertion-warnings)
  - ["Pre-optimization not performed"](#pre-optimization-not-performed)
  - ["Spin multiplicity" or "Check the spin multiplicity"](#spin-multiplicity-or-check-the-spin-multiplicity)
  - [Missing `SH_STATE.dat` or `TRAJEC.xyz`](#missing-sh_statedat-or-trajecxyz)
  - [Plot windows do not appear](#plot-windows-do-not-appear)
  - [Molecular viewer does not open or shows no 3D view](#molecular-viewer-does-not-open-or-shows-no-3d-view)
  - [ESP tool messages](#esp-tool-messages)
  - [NTO tab messages](#nto-tab-messages)
  - [Fukui input builder messages](#fukui-input-builder-messages)
  - [Packaged Windows executable cannot find GLFW](#packaged-windows-executable-cannot-find-glfw)
- [Version, Contributors, and License](#version-contributors-and-license)

## Overview

gqteaWinToga is a desktop molecular analysis toolkit developed by the gQTEA group. It helps users prepare, inspect, convert, visualize, and analyze molecular simulation data from CPMD, Quantum ESPRESSO `cp.x`, ORCA, surface hopping workflows, and related molecular dynamics tools.

It also analyzes Gaussian electronic-structure results without running any external program. From a formatted checkpoint (`.fchk`) it maps the electrostatic surface potential, computes Fukui functions and QTAIM topology, and, with the output of a CIS/TDA/TD job, builds natural transition orbitals (NTOs) and the UV-Vis spectrum of the excited states.

The program is a graphical application built with Python and Toga. Most tools follow the same pattern: open a tool window from the main tabbed interface, browse to one or more input files, enter calculation parameters, then click the tool's action button to generate outputs.

The current application title identifies the toolkit as:

```text
gQTEA-0.9.5 Molecular Analysis Toolkit
```

## Installation and Requirements

### Requirements

Recommended environment:

- Python 3.10 or newer. The toolkit is developed and verified on Python 3.13.
- Windows with a working graphical desktop session.
- A Python virtual environment.

#### Python packages

All Python dependencies are listed in `installation_requirements.txt`, in the
source directory. Install them in one step:

```powershell
python -m pip install -r installation_requirements.txt
```

The file installs these packages:

| Package | Used for |
|---------|----------|
| `toga-winforms` | The GUI framework and its Windows backend. Also installs `toga-core`, `pythonnet` and `pillow`; Pillow draws the live construction view of the solvent-box builders and saves 3D-viewer images. |
| `numpy` | Trajectory parsing, geometry and statistics, used by most analysis tools. |
| `scipy` | The autocorrelation FFT, the molecular pre-optimizer, the solvent-box builders and their periodic pre-optimization (L-BFGS-B, periodic neighbour search, random rotations), and the ESP tool (Boys function, surface-mesh adjacency). |
| `matplotlib` | All plotting, in both the main window and the separate interactive plot viewer. |
| `glfw` | Windowing and mouse input for the 3D molecular viewer and the ESP surface viewer. |
| `PyOpenGL` | OpenGL bindings for the 3D molecular viewer and the ESP surface viewer. |
| `PyOpenGL-accelerate` | Optional speedups for `PyOpenGL`. The toolkit runs correctly without it, only more slowly in the 3D viewer. |
| `PyMuPDF` | Renders `cpx_input_description.pdf` inside the `cp.x` input builder's help viewer. Without it, only the in-app PDF help is lost. |

Verified working on Python 3.13.12 with `toga-winforms` 0.5.3, `numpy` 2.4.2,
`scipy` 1.17.0, `matplotlib` 3.10.8, `glfw` 2.10.0, `PyOpenGL` 3.1.10,
`PyOpenGL-accelerate` 3.1.10 and `PyMuPDF` 1.27.2.2.

Two native libraries that the 3D viewer needs — `glfw3.dll` and `freeglut` —
are bundled inside the `glfw` and `PyOpenGL` wheels, so a normal pip install
requires no separate download. A packaged `.exe` build is different; see
[Packaged Windows executable cannot find GLFW](#packaged-windows-executable-cannot-find-glfw).

#### System requirements that pip cannot install

Three things are not Python packages and cannot be installed from
`installation_requirements.txt`. They are listed as comments at the end of that
file as well.

**1. tkinter (Tcl/Tk)**

A standard-library module that is not distributed on PyPI. The interactive plot
viewer uses matplotlib's `TkAgg` backend, which requires it.

On Windows, tkinter is included with the official python.org installer as long
as the **tcl/tk and IDLE** component is selected. To add it to an existing
installation, re-run the installer, choose **Modify**, and enable that
component. Check whether it is present with:

```powershell
python -c "import tkinter; print(tkinter.TkVersion)"
```

If tkinter is missing, plot windows fail to open **silently**, with no error
message. See [Plot windows do not appear](#plot-windows-do-not-appear).

**2. A graphical desktop session with working OpenGL drivers**

The 3D molecular viewer opens a native OpenGL window. It cannot run over a
plain SSH session, in a headless container, or on a machine with no GPU driver
installed. Install the vendor driver for your graphics hardware (Intel, NVIDIA
or AMD); the Microsoft Basic Display Adapter driver that Windows falls back to
does not provide an OpenGL runtime. Remote Desktop sessions may work through a
software renderer, but slowly.

**3. Microsoft .NET runtime**

`toga-winforms` builds the interface on WinForms through `pythonnet`. Windows 10
and 11 ship a suitable .NET Framework, so no action is normally needed. A
stripped-down or Server Core installation may not include it.

#### Optional packages

None of these is needed to run the toolkit.

`numba` speeds up three tools with compiled kernels (identical results). Since
version 0.9.0 it is part of `installation_requirements.txt`, so the install
command above already includes it. It stays optional: if it cannot be installed
on your Python version, comment out the `numba` line in that file and the tools
use their NumPy paths. To add it to an existing installation:

```powershell
python -m pip install numba
```

- The ESP calculation of the
  [Electrostatic Surface Potential (ESP)](#electrostatic-surface-potential-esp)
  tool becomes about five times faster.
- Electron densities on grids become about four times faster. This speeds up
  the density grid of the ESP tool, and the grids and Hirshfeld charges of the
  [Analysis of Fukui Functions](#analysis-of-fukui-functions) tool.
- Orbital values on grids use the same compiled kernel, so the hole and
  particle isosurfaces and cube files of the
  [NTO tab](#natural-transition-orbitals-nto) are computed about ten times
  faster (uracil, one hole/particle pair at 0.2 bohr: 0.2 s instead of 2 s).
- Solvent packing in the three solvation-box builders
  ([Single Solute](#single-solute-solvent-box),
  [Mixture](#mixture-of-two-solvent-box),
  [Shared-Wall](#shared-wall-double-solvent-box)) becomes about 20–50 times
  faster than the original insertion loop, and 3–5 times faster than their
  NumPy engine; the box is identical for the same random seed.

It is picked up automatically when installed (the solvent builders' *Automatic*
packing engine uses it) and verified with `numba` 0.68.0 / `llvmlite` 0.50.0,
which accept `numpy` 2.4.2 without downgrading it. The first run compiles the
kernels for a few seconds; later runs reuse the compiled cache.

To run the automated test suite:

```powershell
python -m pip install pytest
python -m pytest tests/ -q
```

To build the standalone Windows executable:

```powershell
python -m pip install auto-py-to-exe
```

All three are listed as commented entries in `installation_requirements.txt`.

#### External simulation programs

Some workflows also require files produced by external chemistry programs, such
as CPMD, Quantum ESPRESSO `cp.x`, ORCA, Gaussian, or Vanderbilt `runatom.x`.
gqteaWinToga prepares and analyzes files for these programs, but it does not
replace the external simulation engines.

The Electrostatic Surface Potential tool (the ESP tab of the 3D Molecular Viewer) reads Gaussian formatted checkpoints
(`.fchk`/`.fch`), made with Gaussian's `formchk` utility. It does not need
Gaussian itself; if Gaussian's `cubegen` is installed it can optionally be used
as the ESP engine.

The NTO tab (natural transition orbitals and the UV-Vis spectrum) reads a
Gaussian `.fchk` together with the `.log`/`.out` of the same CIS, TDA or TD
job. Gaussian is needed only to produce those files; for accurate NTOs, run the
job with `IOp(9/40=4)` in the route (see
[Natural Transition Orbitals (NTO)](#natural-transition-orbitals-nto)).

### Install from Source

From PowerShell:

```powershell
git clone https://github.com/AdemirCamargo57/gqteaWinToga.git
cd gqteaWinToga
python -m venv venv
.\venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r installation_requirements.txt
```

If you already have the source tree and an active virtual environment, run the
last command in that environment.

To confirm the installation before starting the program:

```powershell
python -c "import toga, numpy, scipy, matplotlib, glfw, OpenGL, fitz; print('all packages import correctly')"
python -c "import tkinter; print('tkinter', tkinter.TkVersion)"
```

### Start the Program

Run the launcher from the source directory:

```powershell
python gqteaWinToga.py
```

The main window opens with five tabs:

- **Geometry**
- **Inputs**
- **Structural**
- **Thermo**
- **Tools**

Choose a tab, click a tool, complete the fields in the tool window, and run the calculation or conversion.

## Input File Basics

### XYZ Files

Many tools read XYZ geometry or trajectory files. A standard single-frame XYZ file has:

```text
5
Comment line
P  0.112  -0.092   0.064
O  0.112   1.475  -0.490
O  1.470  -0.876  -0.490
O -1.246  -0.876  -0.490
O  0.112  -0.092   1.727
```

The first line is the number of atoms. The second line is a comment. Each remaining line contains an atomic symbol and Cartesian coordinates.

Trajectory XYZ files repeat this block for each frame. CPMD trajectory analysis tools commonly expect a file named `TRAJEC.xyz`.

### Atom Labels

Most analysis tools use 1-based atom labels. For example, `1 2` means atom 1 and atom 2 in the XYZ frame, not Python-style zero-based indices.

### Common Units

- Distances: Angstrom.
- Angles: degrees.
- Temperature: kelvin.
- CPMD time step: atomic units unless a tool states otherwise.
- Some collision outputs distinguish between CPMD atomic-unit velocities and gqteaMD Angstrom/fs velocities.

## Basic Workflows

### Analyze a Bond Length

1. Start gqteaWinToga.
2. Open **Geometry**.
3. Click **Bond length analysis**.
4. Browse to `TRAJEC.xyz`.
5. Enter the maximum distance for the distribution, time step, sampling interval, temperature, two atom labels, and histogram bin width.
6. If the trajectory is periodic and unwrapped, enter **Cell lengths a b c**; leave it blank for an isolated system.
7. Keep **Show plots at the end** enabled if you want plot windows.
8. Keep **Save CSV outputs** enabled if you want CSV files in addition to `.dat` files.
9. Enable **Jacobian r² correction (PMF)** if you want a potential of mean force rather than `-RT ln P`.
10. Click **Analyze**.

Typical outputs include bond length time series, distribution data, free-energy data, summary text, and optional CSV files. Progress appears as text in the output box, and the run ends by showing a summary of the parameters used and the results obtained. Figures open in the interactive viewer; no PNG files are written.

### Run a QTAIM Topology Analysis

1. Run a Gaussian single point (or optimization) that keeps a checkpoint file, and convert it: `formchk name.chk name.fchk`.
2. Open **Tools > 3D Molecular Viewer** and choose the **QTAIM** tab.
3. Browse to `name.fchk`; optionally choose an output folder.
4. Leave the four starting-point switches on and the Advanced fields blank (Multiwfn's defaults).
5. Click **Run** and check that the summary says *Poincare-Hopf = 1 (satisfied)*.
6. Click **Display** to see the critical points and bond paths with the molecule; click a table row to locate a point.
7. Optionally set the **Canvas** (background, *Design style*, element symbols, atom numbers) and click **Export 3D image...**.
8. Choose a plane and click **Plot 2D map**; **Export 2D image...** saves it.
9. Read `name_QTAIM.txt` in the output folder.

The full procedure, with every option and what to do when points are missing, is in [QTAIM Topology Analysis](#qtaim-topology-analysis) (*Step-by-step procedure*); a worked example with real numbers is [Example 9](#example-9-qtaim-analysis-of-a-hydrogen-bond-water-dimer); the theory is in [Theoretical Background: the Quantum Theory of Atoms in Molecules](#theoretical-background-the-quantum-theory-of-atoms-in-molecules).

### Generate Natural Transition Orbitals

1. Run a Gaussian CIS, TDA or TD job that keeps a checkpoint, with `IOp(9/40=4)` in the route (for example `# PBE1PBE/6-31G* TD(NStates=5) IOp(9/40=4)`), and convert it: `formchk name.chk name.fchk`.
2. Open **Tools > 3D Molecular Viewer** and choose the **NTO** tab.
3. Browse to `name.fchk`; the Gaussian output `name.log` (or `name.out`) beside it is filled in. Click **Load**.
4. Choose an excited state and click **Run**; read the pair weights in the table.
5. Choose a pair, **Show** (hole, particle or both) and click **Display**.
6. **Export cubes** / **Export NTO .fchk** write to the output folder; **Export 3D image...** saves the view.
7. **Plot Spectrum** draws all the excited states as a GaussView-style UV-Vis spectrum and saves it as `name_UVVis.json`.

The details are in [Natural Transition Orbitals (NTO)](#natural-transition-orbitals-nto), a worked example in [Example 10](#example-10-natural-transition-orbitals-of-uracil-s3), and the method in [Formulas and Methods (NTO tab)](#formulas-and-methods-nto-tab).

### Build CPMD Input Files

1. Open **Inputs**.
2. Click **CPMD Inputs**.
3. Enter a prefix, total charge, periodic cell parameters, energy cutoff, and dual value.
4. Select a cell symmetry.
5. Browse to an XYZ starting geometry.
6. Click **Input Builder**.

The generated CPMD input files are based on the selected geometry and values entered in the form.

### Select Frames and Generate Gaussian Inputs

1. Open **Tools**.
2. Click **Select Frames**.
3. Enter the starting frame, number of skipped frames, and stop frame.
4. Browse to `TRAJEC.xyz`.
5. Optionally enable Gaussian input generation and fill in checkpoint name, memory, processor count, route section, charge, and multiplicity.
6. Click **Select Frames**.

The tool writes `selected_frames.xyz`. If Gaussian generation is enabled, it also writes files named like `g16_input_1.gjf`, `g16_input_2.gjf`, and so on.

### Plot CPMD Energy Data

1. Open **Tools**.
2. Click **Energy plots**.
3. Select the energy-file type.
4. Browse to the energy file.
5. Enter the simulation time step if needed.
6. Select the x-axis unit and one or more plot types.
7. Click **Plot**.

For CPMD ENERGY files, the program expects a valid numeric table with eight columns.

## Features and Commands

The program is GUI-based. The "commands" below refer to the buttons and controls in each tool window.

### Main Window

The main window organizes tools by category:

- **Geometry**: bond length, all bond distances, bond angle, all bond angles, dihedral angle, all dihedral angles, and hydrogen bond analysis.
- **Inputs**: CPMD, surface hopping, collision, Quantum ESPRESSO `cp.x`, ORCA, gqteaMD, and Vanderbilt `runatom.x` input builders.
- **Structural**: radial distribution function, mean residence time, autocorrelation function, and solvent box builders.
- **Thermo**: classical rate constant calculation.
- **Tools**: visualization, plotting, frame selection, converters, force processing, coordinate tools, and unit conversion.

### Geometry Tools

#### Bond Length Analysis

Use this tool to calculate bond lengths across a trajectory, average bond length, bond-length distribution, and Helmholtz free energy from the probability distribution.

Inputs:

- Maximum `r` for the distribution function.
- Simulation time step in atomic units.
- Sampling interval in frames.
- Simulation temperature in K.
- Two atom labels, for example `1 2`.
- Histogram bin width in Angstrom.
- Optional **Cell lengths a b c** in Angstrom.
- `TRAJEC.xyz` trajectory file.

Options:

- **Show plots at the end**: opens plot windows after the run.
- **Save CSV outputs**: writes CSV versions of the main data tables.
- **Jacobian r² correction (PMF)**: see below.

Periodic boundaries:

- If **Cell lengths a b c** is filled, the distance uses the orthorhombic **minimum-image convention**; leave it blank for an isolated or already-wrapped system.
- Minimum image is exact only up to half the smallest box length, so with cell lengths set, the maximum `r` must not exceed `min(a,b,c)/2`. The tool refuses to run otherwise.
- Coordinates do not need to be wrapped into the box.

Free energy:

- By default the tool reports `G(r) = -RT ln P(r)`, where `P(r)` is the fraction of frames falling in each bin.
- With **Jacobian r² correction (PMF)** enabled it reports the potential of mean force `W(r) = -RT ln[P(r)/r²]`, which removes the `4πr²` volume-element bias of sampling a separation in three dimensions. Both are defined only up to an additive constant, so compare shapes and well depths, not absolute values.

Bond lengths at or beyond the maximum `r` are excluded from the distribution and from the free energy. When that happens the tool reports how many frames were dropped, and records the count in the summary file.

Progress reporting:

- There is no progress bar and no frame-count label. Loading the trajectory and running the calculation both report as text in the large output box.
- When the run finishes, the box shows only a summary: a **PARAMETERS** block listing the values used, and a **RESULTS** block with the bond statistics, the lowest free energy, any excluded bonds, and the output directory. The same text is saved as `summary_<atom1>_<atom2>.txt`.

Typical outputs:

- `bond_<atom1>_<atom2>.dat`
- `bond_distribution_<atom1>_<atom2>.dat`
- `free_energy_<atom1>_<atom2>.dat`
- `summary_<atom1>_<atom2>.txt`
- Optional CSV files with matching names.

Plots open in the interactive viewer (zoom, pan, and save from its toolbar). No static PNG image files are written next to the trajectory; use the viewer's save button if you want an image.

#### All Bond Distance Analysis

Use **Geometry > All bond distance analysis** to calculate distance statistics for every connected atom pair in an XYZ trajectory. This tool is intended for complete solute geometry monitoring in a molecular dynamics run. For a solute in a solvent box, provide the solute atom indices so solvent-solvent and solute-solvent distances are ignored.

How connectivity is detected:

- Connectivity is re-evaluated **on every frame**: in each frame, two in-scope atoms are bonded when their distance is less than or equal to the **Maximum connection distance**.
- The reported pair set is the **union** of every pair bonded in at least one frame, so bonds that form or break mid-run are captured (not only those present in the first frame).
- Each pair's average distance, variance, and standard deviation are computed **only over the frames in which it is bonded**, and the tool also reports the **occurrence** (fraction of frames the pair is bonded).
- The default maximum connection distance is `1.7` Angstrom.
- If **Solute atom indices** is filled, only pairs where both atoms belong to the solute are considered.
- If **Cell lattices a b c** is filled, distances use the orthorhombic **minimum-image convention** (periodic boundaries); leave it blank for an isolated system.

Inputs:

- `trajectory.xyz` or `TRAJEC.xyz` trajectory file.
- **Maximum connection distance (A)** — **leave blank to use the default `1.7`**.
- **Cell lattices a b c (A)** — three positive box lengths for PBC; leave blank for no periodic boundaries. When used, the connection distance must not exceed half the smallest lattice.
- **Solute atom indices** — 1-based labels; **range syntax is supported**, e.g. `1-5 14-16 18 20` is equivalent to `1 2 3 4 5 14 15 16 18 20`. Whitespace and commas both separate tokens. Leave blank to use all atoms.
- Output TXT filename, default `all_bond_analysis_combined.txt`.

Step-by-step:

1. Open **Geometry > All bond distance analysis**.
2. Click **Browse** and select the trajectory XYZ file.
3. Enter the maximum connection distance, or leave it blank for the `1.7` default.
4. To apply periodic boundaries, enter the cell lattices `a b c`; otherwise leave blank.
5. If the system contains a solute in solvent, enter the solute atom labels in **Solute atom indices** (ranges allowed, e.g. `1-5 14-16 18 20`).
6. Enter the output filename or keep `all_bond_analysis_combined.txt`.
7. Click **Analyze**.
8. Wait until the final summary appears in the text area.
9. Open the statistics file and pair mapping file in the trajectory directory.

Output file:

The tool writes a **single self-describing file**, `all_bond_analysis_combined.txt`, that merges the pair identity, distance statistics, and run metadata. It begins with a commented header block:

```text
# gQTEA All Bond Distance Analysis
# frames_used <N>
# max_connection_distance <cutoff>
# atom_scope all_atoms            (or: # atom_scope solute_atoms + # solute_atom_indices ...)
# periodic_boundary none          (or: # cell_lengths a b c)
# row atom_i atom_j element_i element_j first_bonded_distance average variance standard_deviation occurrence_fraction frames_bonded
```

Each following row is one connected atom pair (atom labels 1-based):

```text
row  atom_i  atom_j  element_i  element_j  first_bonded_distance  average  variance  standard_deviation  occurrence_fraction  frames_bonded
```

Because the data rows include the two element-symbol columns, load the numeric columns with `np.loadtxt(path, usecols=(0,1,2,5,6,7,8,9,10))`, or use `np.genfromtxt(path, dtype=None, names=True, comments='#')` to keep the element strings.

#### Bond Angle Analysis

Use this tool to calculate bond angles across a trajectory, angle distributions, and free-energy profiles.

Inputs:

- Maximum angle for the distribution function.
- Simulation time step in atomic units.
- Sampling interval in frames.
- Simulation temperature in K.
- Three atom labels `i j k`, for example `1 2 3`. The vertex is the middle label. All three must be different.
- Histogram bin width in degrees.
- Optional **Cell lengths a b c** in Angstrom.
- `TRAJEC.xyz` trajectory file.

Options:

- **Show plots at the end**.
- **Save CSV outputs**.
- **Use sin(theta) Jacobian**: reports the potential of mean force `W(theta) = -RT ln[P(theta)/sin(theta)]` instead of `-RT ln P(theta)`, removing the volume-element bias of sampling an angle in three dimensions. Off by default.
- **Smooth free energy**: applies a small moving average to the probabilities before the logarithm, which tames noise at fine bin widths. **On by default** for this tool. Turn it off for the raw curve.

Periodic boundaries:

- If **Cell lengths a b c** is filled, both arms of the angle use the orthorhombic **minimum-image convention**; leave it blank for an isolated or already-wrapped system.
- Minimum image is exact only for arms shorter than `min(a,b,c)/2`. The tool warns when an arm exceeds that rather than failing silently.

Angles at or beyond the maximum angle are **excluded and reported**, not clamped into the last bin. The run ends showing only a summary: a **PARAMETERS** block and a **RESULTS** block, also saved as `summary_<i>_<j>_<k>.txt`. Progress appears as text in the output box; there is no progress bar or frame-count label. Figures open in the interactive viewer, and no PNG files are written.

#### All Bond Angle Analysis

Use **Geometry > All bond angle analysis** to calculate angle statistics for every connected `atom_i-atom_j-atom_k` angle in an XYZ trajectory. The central atom is `atom_j`. Connectivity is **re-evaluated on every frame**: a triplet exists in a frame when both of its bonds (`atom_j-atom_i` and `atom_j-atom_k`) are within the cutoff. The reported set is the union of every triplet that appears in at least one frame, each triplet's statistics are averaged only over the frames where it exists, and the tool reports the **occurrence** (fraction of frames the triplet is connected).

For a solute in a solvent box, provide the solute atom indices. The program then builds angles only from solute atoms, excluding any angle involving solvent atoms.

Inputs:

- `trajectory.xyz` or `TRAJEC.xyz` trajectory file.
- **Maximum connection distance (A)** — **leave blank to use the default `1.7`**.
- **Cell lattices a b c (A)** — three positive box lengths for PBC (orthorhombic minimum-image); leave blank for no periodic boundaries.
- **Solute atom indices** — 1-based labels; **range syntax is supported**, e.g. `1-5 14-16 18 20`. Leave blank to use all atoms.
- Output TXT filename, default `all_angle_analysis_combined.txt`.

Step-by-step:

1. Open **Geometry > All bond angle analysis**.
2. Click **Browse** and select the trajectory XYZ file.
3. Enter the maximum connection distance, or leave it blank for the `1.7` default.
4. To apply periodic boundaries, enter the cell lattices `a b c`; otherwise leave blank.
5. Enter solute atom labels if solvent atoms must be excluded (ranges allowed).
6. Enter the output filename or keep `all_angle_analysis_combined.txt`.
7. Click **Analyze**.
8. Wait until the text area reports that the analysis is completed.
9. Inspect the combined output file.

Output file:

A single self-describing file, `all_angle_analysis_combined.txt`, with a commented metadata header (`# frames_used`, `# max_connection_distance`, `# atom_scope`, and `# periodic_boundary none` / `# cell_lengths a b c`) followed by one row per triplet:

```text
# row atom_i atom_j atom_k element_i element_j element_k first_present_angle_degrees average_angle_degrees variance_degrees2 standard_deviation_degrees occurrence_fraction frames_present
```

Atom labels are 1-based. Because the rows include element-symbol columns, load numeric columns with `np.loadtxt(path, usecols=...)` or `np.genfromtxt(path, dtype=None, names=True, comments='#')`.

#### Dihedral Angle Analysis

Use this tool to calculate dihedral angles across a trajectory, dihedral distributions, and free-energy profiles.

Inputs:

- Maximum angle for the distribution.
- Simulation time step in atomic units.
- Sampling interval in frames.
- Simulation temperature in K.
- Four atom labels `i j k l`, for example `1 2 3 4`. All four must be different.
- Histogram bin width in degrees.
- Optional **Cell lengths a b c** in Angstrom.
- `TRAJEC.xyz` trajectory file.

Options:

- **Show plots at the end**.
- **Save CSV outputs**.
- **Wrap dihedral to [-180, 180]**: reports signed dihedral angles instead of a 0 to 360 degree range. With it on, the histogram spans `[-180, 180]` so negative angles are binned properly; with it off, angles map to `[0, 360)` and the histogram spans `[0, maximum angle]`.
- **Smooth free energy**: applies a small moving average to the probabilities before the logarithm. **On by default** for this tool.

There is deliberately **no Jacobian option** here. Unlike a bond length (`4*pi*r^2`) or a bond angle (`sin(theta)`), a dihedral has a uniform volume element, so there is nothing to divide out: `-RT ln P(phi)` is already the potential of mean force. The summary records this as `Jacobian: not applicable (uniform measure)`.

Periodic boundaries:

- If **Cell lengths a b c** is filled, all three connecting vectors use the orthorhombic **minimum-image convention**; leave it blank for an isolated or already-wrapped system.
- Minimum image is exact only for bonds shorter than `min(a,b,c)/2`. The tool warns when one exceeds that.

Angles outside the histogram range are **excluded and reported**, not clamped. The run ends showing only a summary, also saved as `summary_<i>_<j>_<k>_<l>.txt`. Progress appears as text in the output box; there is no progress bar or frame-count label. Figures open in the interactive viewer, and no PNG files are written.

#### All Dihedral Angle Analysis

Use **Geometry > All dihedral angle analysis** to calculate statistics for every connected `atom_i-atom_j-atom_k-atom_l` dihedral path in an XYZ trajectory. Connectivity is **re-evaluated on every frame**: a quadruplet exists in a frame when all three of its bonds (`atom_i-atom_j`, `atom_j-atom_k`, `atom_k-atom_l`) are within the cutoff. The reported set is the union of every quadruplet that appears in at least one frame, each quadruplet's statistics are averaged only over the frames where it exists, and the tool reports the **occurrence** (fraction of frames the quadruplet is connected).

For a solute in a solvent box, provide the solute atom indices. The program then builds only solute dihedrals and ignores all solvent atoms.

Dihedral convention:

- Dihedrals are signed angles in degrees.
- The range is `(-180, 180]`.
- This convention matches the single **Dihedral angle analysis** tool.

Inputs:

- `trajectory.xyz` or `TRAJEC.xyz` trajectory file.
- **Maximum connection distance (A)** — **leave blank to use the default `1.7`**.
- **Cell lattices a b c (A)** — three positive box lengths for PBC (orthorhombic minimum-image); leave blank for no periodic boundaries.
- **Solute atom indices** — 1-based labels; **range syntax is supported**, e.g. `1-5 14-16 18 20`. Leave blank to use all atoms.
- Output TXT filename, default `all_dihedral_analysis_combined.txt`.

Step-by-step:

1. Open **Geometry > All dihedral angle analysis**.
2. Click **Browse** and select the trajectory XYZ file.
3. Enter the maximum connection distance, or leave it blank for the `1.7` default.
4. To apply periodic boundaries, enter the cell lattices `a b c`; otherwise leave blank.
5. Enter solute atom labels if solvent atoms must be excluded (ranges allowed).
6. Enter the output filename or keep `all_dihedral_analysis_combined.txt`.
7. Click **Analyze**.
8. Wait until the text area reports that the analysis is completed.
9. Inspect the combined output file.

Output file:

A single self-describing file, `all_dihedral_analysis_combined.txt`, with a commented metadata header (`# frames_used`, `# max_connection_distance`, `# dihedral_convention signed_degrees_minus180_to_180`, `# atom_scope`, and `# periodic_boundary none` / `# cell_lengths a b c`) followed by one row per quadruplet:

```text
# row atom_i atom_j atom_k atom_l element_i element_j element_k element_l first_present_dihedral_degrees average_dihedral_degrees variance_degrees2 standard_deviation_degrees occurrence_fraction frames_present
```

Atom labels are 1-based. Because the rows include element-symbol columns, load numeric columns with `np.loadtxt(path, usecols=...)` or `np.genfromtxt(path, dtype=None, names=True, comments='#')`.

#### Hydrogen Bond Analysis

Use this tool to analyze hydrogen bonds from a trajectory.

Inputs:

- Maximum hydrogen bond length.
- Simulation time step in atomic units.
- Sampling interval in frames.
- Simulation temperature in K.
- Histogram bin widths for hydrogen bond length and angle.
- Donor atom label, hydrogen atom label, and acceptor atom label.
- `TRAJEC.xyz` trajectory file.

Click **Exec** or the main action button in the tool window to run the calculation.

#### Comparison of Molecular Geometric Parameters

Use this tool to see how a different environment changes a molecule's geometry. It compares the geometric parameters of the **same molecule** taken from two simulations — typically the isolated molecule and the same molecule in a water box, with the solvent already removed from the trajectory.

Inputs:

- **Parameter file 1 (reference)** and **Parameter file 2 (comparison)**: two combined parameter files written by the All bond distance, All bond angle or All dihedral angle analysis tools. Both must be of the same type; the type is read from each file's own title line.
- **Label for file 1** / **Label for file 2** (default `isolated` and `solvated`): these become part of the output column names, e.g. `average_isolated`.
- **Minimum occurrence fraction** (blank = 0, keep everything): ignores parameters present in fewer than this fraction of frames.
- **Relative difference columns** — two independent switches, either, both, or neither on:
  - **Symmetric percent difference** (on by default): `100*(avg_2-avg_1)/mean(avg_1,avg_2)`.
  - **Percent error vs file 1** (off by default): `100*(avg_2-avg_1)/avg_1`.

  Both are signed, file 2 minus file 1. Leaving both off writes the output file exactly as it was before these columns existed.
- **Plot** — **Open the comparison figure after comparing** (on by default) and **Parameters to plot** (blank = 25, the parameters that shifted most).
- **Output folder** (blank = next to parameter file 1) and **Output txt filename** (blank = `geometry_comparison_<type>.txt`).

Each file is read as soon as you select it, so a wrong or corrupt file is reported immediately, together with its type, frame count, atom scope and molecule size.

Click **Compare** to run, **Help** for the full in-program description.

##### How atoms are matched

Atoms are **not** matched by their raw index. A solvated run is normally analysed with solute atom indices, so its rows carry the atoms' indices in the **full box**: the oxygens of the reference molecule may be atoms 295, 296 and 297 there, while the same atoms are 29, 30 and 31 in the isolated run. The tool reads the `# solute_atom_indices` header and renumbers every atom back onto the molecule, so `C1-O295` in the box is recognised as the same bond as `C1-O29` isolated. The renumbering it applied is recorded in the output header.

Before comparing, the two files must agree on the parameter type, the molecule size, and the element of every atom. If any check fails the run is refused and the offending atom is named — a wrong atom mapping would produce results that look reasonable but are meaningless.

Parameters match regardless of atom order within a row: a bond is undirected, an angle keeps its vertex but may have its arms swapped, and a dihedral read backwards has the same signed value.

##### Why the occurrence filter matters

The all-\* tools list every pair that came within the connection distance in at least one frame, so a long trajectory collects contacts that are not real bonds — typically H–H pairs sitting just under the 1.7 Å cutoff in a fraction of a percent of the frames. Setting the fraction to e.g. `0.05` drops them. The filter is applied to **both** files before matching, so a rare contact cannot be reported as a parameter that one simulation is missing.

##### Output

A single text file with a metadata header (both source paths, labels, frame counts, atom scopes, the atom renumbering applied, the occurrence threshold and what it dropped, and the matched/unmatched counts) followed by three labelled sections:

```text
# section matched
# row atom_i atom_j element_i element_j average_isolated std_dev_isolated occurrence_isolated average_solvated std_dev_solvated occurrence_solvated difference_solvated_minus_isolated source_row_1 source_row_2
# section only_in_file_1
# row atom_i atom_j element_i element_j average std_dev occurrence source_row
# section only_in_file_2
```

Angle and dihedral files carry three and four atom/element columns instead of two. Atom labels are the **molecule-local** 1-based indices; `source_row_1` / `source_row_2` give the row number the values came from in each input file, so any row can be traced back. Distances are in ångströms and angles in degrees. Dihedral differences are wrapped into (−180, 180], so a shift from +179° to −179° is reported as +2°, not −358°.

Every comment line starts with `#`, so a section can be read back with `np.genfromtxt(path, dtype=None, names=True, comments='#')` once it is separated from the others.

After a successful run the text box reports the counts, the output path, and the ten parameters whose average changed most between the two simulations.

##### Relative differences and global metrics

With one or both **Relative difference columns** switches on, a percent column is appended after `difference_solvated_minus_isolated` for each switch that is on — appended, never inserted earlier, so a reader that indexes the columns above by position keeps working:

```text
# row ... difference_solvated_minus_isolated percent_difference percent_error_vs_isolated source_row_1 source_row_2
```

`percent_difference` is the symmetric form; `percent_error_vs_isolated` carries the reference file's label so two identically labelled comparisons stay distinguishable. Both are signed, file 2 minus file 1. A percentage whose denominator is too close to zero to mean anything — below `1e-6` Å for a bond, below `1.0°` for an angle or a dihedral — is written `nan`; for a dihedral the symmetric denominator is the circular mean of the two angles, not their arithmetic mean, so a shift between e.g. +179° and −179° does not divide by (almost) zero.

The header also carries the formula actually used, the denominator floor, the count of undefined rows, and global metrics over the whole matched set:

```text
# formula_percent_difference 100*(avg_2-avg_1)/((avg_1+avg_2)/2)
# formula_percent_error 100*(avg_2-avg_1)/avg_1
# percent_denominator_floor 1e-06 angstrom
# percent_undefined_rows 0
# global_metric_category bond distance
# global_metric_n 42
# mae_angstrom 0.01234567
# rmsd_angstrom 0.01567890
# mae_percent_difference 1.23456789
# rmsd_percent_difference 1.56789012
```

`mae_<unit>` / `rmsd_<unit>` are the mean absolute error and root-mean-square deviation over every matched parameter, in its own unit (ångströms or degrees); the `mae_percent_*` / `rmsd_percent_*` pair is added per selected switch, over the absolute value of that percentage, skipping undefined rows. With both switches off, none of the `formula_*`, `percent_*`, `global_metric_*`, `mae_*` or `rmsd_*` lines are written and the two extra columns are absent — the file is exactly what it always was.

##### What the message panel reports

After a successful comparison the panel summarises both files, the molecule size, the counts of matched and unmatched parameters, and the global metrics, then lists the parameters that shifted most. Each switch that is on adds its percentage to every row of that list:

```
Global metrics over 6 matched bond distances:
  MAE  0.01008 A     RMSD  0.01454 A
  MAE  0.8198 %     RMSD  1.1746 %   (percent_difference)
  MAE  0.8267 %     RMSD  1.1861 %   (percent_error)

Largest shifts (solvated - isolated):
  C2-O3                   1.41220 ->    1.44090     +0.02870 A    +2.01 %    +2.03 %
  O3-H4                   0.97120 ->    0.99060     +0.01940 A    +1.98 %    +2.00 %
  C1-C2                   1.52830 ->    1.53610     +0.00780 A    +0.51 %    +0.51 %
```

The percent columns appear in the same order as in the output file — symmetric difference first, then percent error — and a parameter whose percentage is undefined shows `nan` there without disturbing the rest of the row. With both switches off the list carries no percent columns at all.

##### The comparison figure

With **Open the comparison figure after comparing** on, a grouped bar chart shows both files' averages for the **Parameters to plot** parameters that shifted most (ranked by the size of the shift), each bar carrying a ±1 standard deviation error bar taken from its source file. If a percentage switch is on, that percentage is drawn as a line on a right-hand axis; with both switches on, the line always draws the symmetric percent difference, and the legend says so. The figure opens in the same interactive viewer (zoom/pan/save) as every other tool's plots, and no `.png` file is left next to the parameter files.

##### Reopening the figure later

The viewer writes a small `.json` figure manifest into the output folder when it opens. To see the chart again without re-running the comparison, open **Tools > Plots**, set **Plot type** to `JSON plot file`, and select that file — the Plotter understands the bar-chart manifest and reopens the figure exactly as it was, error bars, percentage line and all.

### Input Builders

#### CPMD Inputs

Use **Inputs > CPMD Inputs** to create CPMD input files from an XYZ geometry.

Fields:

- **Prefix file name**: prefix used for generated files, for example `VitC`.
- **Total charge on the system**.
- **Periodic box parameters**: `a b c cosA cosB cosC`, separated by spaces.
- **Energy Cutoff (Ry)**.
- **Dual for rho expansion**.
- **Cell symmetry**:
  - `1 - CUBIC a=b=c alpha=beta=gamma=90`
  - `6 - TETRAGONAL a=b!=c alpha=beta=gamma=90`
  - `8 - ORTHORHOMBIC a!=b!=c alpha=beta=gamma=90`
- XYZ input file.

Action:

- **Input Builder**: generates the CPMD input files.

#### Surface Hopping Input Builder

Use **Inputs > SH Input Builder** to prepare CPMD surface hopping and TDDFT input files.

Fields:

- Starting frame.
- Number of skipped frames.
- Final frame.
- Number of electronic states.
- Initial state for SHTDDFT.
- Prefix file name.
- Molecular dynamics time step.
- Total system charge.
- Periodic box parameters.
- Energy cutoff.
- Dual for rho expansion.
- Cell symmetry.
- Maximum number of steps.
- `GEOMETRY.xyz` input file.
- `TRAJECTORY` input file.

Actions:

- **Read SH parms**: reads and displays the entered parameters.
- **SH Input Builder**: generates surface hopping input files.
- **Help**: opens the tool help text.

#### Collision Input

Use **Inputs > Collision Input** to create new geometry/velocity files for collision molecular dynamics.

Fields:

- **MD engine**: select `cpmd` or `gqteaMD`.
- **Attacker atom indices**: one or more 1-based atom labels, for example `12 13 14`.
- **Initial velocity**: atomic units for CPMD; the label changes for gqteaMD.
- **Target xyz coordinates**: target point in Angstrom, for example `10.0 12.5 8.0`.
- `GEOMETRY` or `GEOMETRY.xyz` file.

Options:

- **Zero velocities of non-attacker atoms**.
- **Compute center-of-mass KE**.

Typical outputs:

- `summary.txt`
- `newGeometry.xyz`
- `with-vibration-GEOMETRY`
- `without-vibration-GEOMETRY`
- `gqteaMD-GEOMETRY`
- `gqteaMD-without-vibration-GEOMETRY`

For CPMD collision simulations, the help text recommends restarting CPMD with:

```text
RESTART WAVEFUNCTION COORDINATES VELOCITIES GEOFILE LATEST
```

#### Quantum ESPRESSO cp.x Input Builder

Use **Inputs > cp.x Input Builder** to build Quantum ESPRESSO `cp.x` input files.

Main actions:

- **Open XYZ**: loads atomic coordinates from an XYZ file.
- **Preview input**: opens a preview of the generated input.
- **Save cp.x input**: writes the input file.
- **Help**: opens the bundled `cpx_input_description.pdf`.

Important behavior:

- If the output name is blank, the default name is based on the XYZ file, such as `<xyzname>_cp.in`.
- `ATOMIC_SPECIES` is inferred from XYZ symbols.
- Default pseudopotential names use the pattern `<Element>.UPF`; edit them if your pseudopotential filenames differ.
- `ATOMIC_POSITIONS` can be written in `alat`, `bohr`, `crystal`, or `angstrom`.
- Optional sections include velocities, cell parameters, constraints, occupations, atomic forces, Wannier plotting, and autopilot text.

#### ORCA Input Builder

Use **Inputs > ORCA Input Builder** to prepare ORCA `.inp` files from XYZ coordinates.

Main actions:

- **Open XYZ**.
- **Preview ORCA input**.
- **Save ORCA input**.
- **Clear form**.
- **Research notes**.

Common settings:

- Output filename. If blank, the default is `<xyzname>_orca.inp`.
- Title/comment line.
- Method, basis, task, charge, multiplicity, and coordinate mode.
- Extra simple keywords.
- Parallel settings such as `%pal nprocs`.
- Memory setting such as `%maxcore`.
- SCF, geometry, CPCM, and advanced ORCA blocks.

#### gqteaMD Input Builder

Use **Inputs > gqteaMD Input Builder** to generate a gqteaMD TOML input file.

Main actions:

- **Open XYZ**.
- **Preview TOML**.
- **Save TOML**.
- **Clear form**.
- **Manual notes**.

Common settings:

- Output filename. If blank, the default is `<xyzname>_gqteaMD.toml`.
- Starting XYZ path.
- Force provider.
- Time step, number of steps, temperature, output intervals, restart behavior, and force-field sections.
- For xTB force provider, **OMP threads** writes `omp_num_threads` to the TOML so gqteaMD can set `OMP_NUM_THREADS` before xTB runs.
- Optional UFF or classical parameter blocks.

#### Gaussian Fukui Input Builder

Use **Inputs > Gaussian Fukui Input Builder** to prepare the three Gaussian calculations that [Tools > Analysis of Fukui Functions](#analysis-of-fukui-functions) needs. The input is written for **Gaussian 16 (Rev. C.01) under Linux**, either in WSL or on a computing cluster. It is **one** input file with three linked jobs, separated by `--Link1--`:

| Job | State | Calculation | Charge | Multiplicity | Checkpoints |
|---|---|---|---|---|---|
| 1 | N | `Opt` | Q | M | `%Chk=<name>_N.chk` |
| 2 | N+p | single point, `Geom=Check Guess=Read` | Q − p | M₊ | `%OldChk=<name>_N.chk`, `%Chk=<name>_Np1.chk` |
| 3 | N−q | single point, `Geom=Check Guess=Read` | Q + q | M₋ | `%OldChk=<name>_N.chk`, `%Chk=<name>_Nm1.chk` |

Job 1 optimizes the N state from the geometry you give. Jobs 2 and 3 read the **optimized** geometry (and, by default, the orbitals as the initial guess) from `<name>_N.chk`. The ions are **not** re-optimized, because the Fukui functions are derivatives at fixed nuclei. Gaussian copies the `%OldChk` file into each job's own `%Chk` before the job starts, so `<name>_N.chk` keeps the optimized N state. All jobs share the same method, basis set and extra keywords.

**What you must supply.** Some values cannot be safely defaulted:
- the starting geometry;
- the charge of the N state;
- **all three spin multiplicities** (N, N+p, N−q). The builder never assumes a spin state. A blank multiplicity stops Preview and Save with a message that names the state and lists the values its electron count allows;
- `%Mem` and `%NProcShared` for the machine that runs Gaussian;
- the output folder.

**Setup tab.**

| Field | Default | Notes |
|---|---|---|
| Job name | name of the loaded file (`molecule` otherwise) | Letters, digits, `.`, `_`, `-`; it becomes the file names. Loading `<name>_N.fchk` gives `<name>`. |
| Output folder | folder of the loaded file | Must exist. |
| `%Mem` | `2GB` | Integer plus a unit (`KB`, `MB`, `GB`, `TB`, `KW`, …). A bare number is refused: Gaussian would read it as 8-byte *words*. |
| `%NProcShared` | 4 | Cores of the WSL machine or the cluster node. |
| Print level | `#P` | `#P`, `#`, `#N`, `#T`. |
| Method | `B3LYP` | No R/U prefix is needed: Gaussian uses R for a closed-shell state and U for an open-shell one. An `R` prefix with an open-shell state warns. |
| Basis set | `6-31G(d)` | A diffuse basis (`+`, `aug-`) warns; see the analysis tool's caveats. |
| Extra keywords (all jobs) | `SCF=Tight NoSymm` | `NoSymm` keeps one orientation in all three jobs. |
| Optimization (job 1) | `Opt` | One `Opt` keyword towards a minimum, e.g. `Opt=Tight` or `Opt=(Tight,CalcFC)`. `TS`, `QST2`/`QST3`, `ModRedundant`, `ReadFC` and options that read input are refused. |

The **Route** lines show the route of job 1 and of jobs 2–3 as they will be written.

Some keywords are **refused** in *Extra keywords*:
- `Opt` (it belongs in the Optimization field, so that only the N state is optimized);
- `IRC`, `Scan`;
- `Geom=` and `Guess=` (the builder writes them for the ions);
- `Units=` (the coordinates are written in Å);
- `Gen`/`GenECP`, `ExtraBasis`, `Pseudo`/ECPs, and any option that reads an extra input section.

Others only produce a **warning**:
- `Freq` (it would run on all three jobs and is not needed);
- post-SCF methods and `TD` (the analysis reads only the `Total SCF Density`);
- diffuse basis sets.

**Molecule & states tab.**
- **Load geometry...** reads `.xyz`, `.gjf`/`.com`, `.fchk`/`.fch`, `.pdb`, `.mol`, `.mol2` and `.hin`. A `.gjf` or `.fchk` also sets the charge and multiplicity of the N state; the ion multiplicities are cleared. You can also paste or type `Element x y z` lines in Å; an element symbol, a numbered label (`C1`) or an atomic number is accepted. Z-matrices, freeze flags and fragment tags are refused, with the line number in the message. This is the **starting** geometry: job 1 optimizes it.
- **Charge of the N state**: blank = 0.
- **Multiplicity of the N state, of N+p and of N−q**: required. For ordinary organic molecules the lowest states are usually 1, 2, 2 for a closed-shell neutral and 2, 1, 1 for a neutral radical, but you must enter them. For transition-metal systems, compute several spin states of each ion and use the lowest in energy.
- **p / q**: 1 for ordinary f⁺/f⁻. Larger values are for degenerate frontier orbitals (see *Formulas and Methods (Fukui Function tool)*, subsection 5).
- The **states** lines show the charge, electron count, `.chk` name and multiplicity of each state as you type. A missing multiplicity shows as `?` with the values allowed; an impossible one is flagged.
- **Start the ions from the optimized N-state orbitals (Guess=Read)**: on by default. When off, the ions get Gaussian's default initial guess; `Geom=Check` is written either way.

**Buttons.**
- **Preview** shows the input on the *Preview* tab, with any warnings.
- **Save input** writes `<name>_fukui.gjf` and `<name>_formchk.sh`, asking before it replaces existing files. Both files use relative file names and Linux (LF) line endings. The *Preview* tab then lists the commands to run.

**Typical workflow.**
1. Load or paste the starting geometry, enter the three multiplicities, and check the Setup tab.
2. Press **Save input**.
3. In a Linux shell with Gaussian 16 set up, go to the folder that holds the input. Under WSL, `C:\jobs\phenol` is `/mnt/c/jobs/phenol`. On a cluster, copy the folder there first.
4. Run `g16 < <name>_fukui.gjf > <name>_fukui.log`, or submit the job through the cluster's queue.
5. Run `bash <name>_formchk.sh`. It runs `formchk <name>_N.chk <name>_N.fchk`, and the same for the other two, and reports any missing checkpoint.
6. Copy the three `.fchk` files back if they were made on a cluster. In *Tools > Analysis of Fukui Functions*, choose `<name>_N.fchk`. The other two files are filled in automatically, and the analysis checks that the three states share the same geometry, basis set and method.

If one linked job fails, Gaussian does not run the jobs after it. The script then reports which `.chk` is missing; the reason is in the `.log` file.

#### Vanderbilt runatom.x Input Builder

Use **Inputs > Vanderbilt runatom.x Input Builder** to generate Vanderbilt `runatom.x` input files.

Controls:

- Mode selection for all-electron or generation input.
- **Also generate Makefile** option.
- **Generate** button.
- **Help** button.

Depending on the selected mode, the form shows the relevant atomic, pseudopotential, and state fields. The generated input is saved through a Save dialog, and a Makefile can be written alongside it.

### Structural Tools

#### Radial Distribution Function

Use **Structural > Radial Distribution Function** to calculate the radial distribution function g(r) and the running coordination number around a chosen central atom, from a `TRAJEC.xyz` trajectory.

Inputs include:

- **Maximum radius for RDF** — the cutoff distance. It must not exceed half the smallest cell lattice; larger values are rejected (a periodic-boundary requirement). **Leave it blank to use the maximum allowed value automatically** (half the smallest cell lattice).
- **Bin width for histogram** — width of the distance bins (defaults to `0.01`).
- **Atom labels to be excluded** — space-separated 1-based labels to leave out, or `0` for none.
- **Shell center atom label** — the 1-based label of the central atom that g(r) is measured from.
- **Atomic symbol for g(r)** — the element counted as neighbours (e.g. `O`). It must be present in the trajectory.
- **Cell lattices (a b c)** — the three orthorhombic box lengths in Å.
- **X-axis range limit for RDF plot** — the x-max for the coordination-number plot. **Leave blank to default to half the a lattice** (`a / 2`).
- **Y-axis range limit for RDF plot** — the y-max for the coordination-number plot. **Leave blank to default to `10`**.
- **Select input file** — click **Browse** to load the `TRAJEC.xyz` trajectory (also counts the frames).

Coordinates do not need to be wrapped into the box; periodic images are handled internally by the minimum-image convention.

Actions:

- **Read Params**: reads and validates the settings and shows a summary. Invalid entries raise an explanatory dialog and stop.
- **RDF calculation**: runs the calculation. Live progress (loading, frame count) and a final completion summary — parameters used, cell volume and number density, the first g(r) peak, the coordination number, and the output file — are shown in the on-screen text area; malformed or truncated frames are reported.
- **Help**: shows usage instructions.
- **Close**: closes the window.

Progress is reported as text in the on-screen status area rather than a progress bar.

Typical outputs (written next to the trajectory):

- `RDF_<center><n>_<symbol>.dat` — three columns: r (bin center), g(r), and the integrated coordination number.
- Interactive plots of g(r) and the coordination number.

#### Mean Residence Time

Use **Structural > Mean Residence Time** for MRT analysis of an XYZ trajectory.

The trajectory field, the **Run MRT** / **Export Results** / **Clear Output** buttons, and the status line are always visible. The rest of the controls are grouped into four tabs.

##### Setup tab

- Time step (Δt between saved frames) and its unit.
- Cutoff radius. Take it from the first minimum of the corresponding g(r) — the [Radial Distribution](#radial-distribution-function) tool will give you that.
- Tolerance frames (`t*`) — bridges short excursions: a particle that leaves for no more than `t*` frames and returns is treated as never having left.
- **Cell lengths a b c** — optional orthorhombic periodic box. Leave blank for no periodic boundary; when given, distances use the minimum-image convention and the cutoff must not exceed half the shortest edge.
- **Scan t\* values** — optional list such as `0 1 2 5 10`. The summary then includes one row per value so you can see how sensitive the result is before quoting a number.
- **Output folder** — optional; blank writes next to the trajectory.

Under **Advanced** on the same tab:

- **Integration rule** — how `R(t)` is turned into a time constant. *Truncate at first zero crossing* (default) integrates up to the point where the curve first goes negative; *Fit a single exponential* fits `ln R = a − t/τ` and reports τ, which is what you want when `R(t)` has not decayed to zero inside the lag range; *Integrate the whole lag range* uses everything.
- **Maximum lag (frames)** — blank uses a tenth of the trajectory. Raise it if the correlation is still decaying at the end of the plotted range.
- **Censor events clipped by the start/end of the run** — on by default. Turn it off only to reproduce a number from an older version; leaving it off biases the mean residence time downwards.

The summary always states which rule was applied and how many lags it covered, so a run is reproducible from its own output. If it reports that all available lags were used, the curve never crossed zero and you should either raise the maximum lag or switch to the exponential fit.

##### Definitions tab

- **Reference mode** — what the shell is centred on: a single atom, or the geometric centre or centre of mass of a group.
- **Observed mode** — what is watched: a single atom, every atom in a list tracked separately, or groups followed by their centre.

Each definition box carries an inline hint showing the exact syntax and a worked example for the mode currently selected, because the syntax changes with the mode (`295` for one atom, `3-7,10,15-17` for a list, `1-3; 4-6` for groups). All atom indices are 1-based.

**Ranges.** Atom indices accept inclusive ranges as well as single values, so

```text
3-7,10,15-17
```

means exactly the same as

```text
3,4,5,6,7,10,15,16,17
```

Ranges and single indices can be mixed in any order, separated by commas or spaces, and they work inside group definitions too (`1-3; 4-6`). The order you type is preserved — useful because the tool labels observed objects in input order. This is the same syntax the [All Bond Distance Analysis](#all-bond-distance-analysis) tool uses for its solute selection — one shared, tested parser.

Ranges must ascend and start at 1. Malformed input is rejected with a message naming the offending token:

| Input | Result |
| --- | --- |
| `7-3` | `Range bounds must be ascending: '7-3'.` |
| `0-5` | `Atom indices must be positive integers starting at 1.` |
| `3-` | `Invalid index range: '3-'.` |
| `-5` | `Invalid index range: '-5'.` |
| `3-5-7` | `Invalid index range: '3-5-7'.` |
| `a-b` | `Invalid index range: 'a-b'.` |

A bad token rejects the whole entry rather than being silently skipped, so `1,2,7-3,9` is refused outright.

##### Results tab

The full summary appears here after a run — parameters, results, the `t*` scan table if you requested one, and notes on the definitions. **Show plots** opens three figures in the interactive viewer: the survival function `S(t)`, the intermittent `C(t)` and `R(t)` overlaid, and the event-duration distribution. No image files are left next to your data.

##### Help tab

The full reference for the tool, including what each reported quantity means.

##### Output files

`<name>_continuous_survival.dat`, `<name>_intermittent_correlation.dat`, `<name>_mrt_event_durations.dat`, and `<name>_mrt_summary.txt`.

##### The two continuous time constants

The summary reports two different quantities, and they are not interchangeable:

- **Mean residence time `<T>`** — the average duration of a complete residence event, with a block-averaged uncertainty. Events that were already in progress when the trajectory started, or had not ended when it stopped, are *censored* (excluded) and counted separately, because including them at their truncated length biases `<T>` downwards, worst for the longest-lived species.
- **IMM survival integral** — the integral of the survival function averaged over time origins, in the sense of Impey, Madden and McDonald (*J. Phys. Chem.* **87** (1983) 5071). Mathematically this is the mean *residual* time; it equals `<T>` only for exponential kinetics and is larger for broad duration distributions.

The intermittent time constant from `R(t)` is integrated only over a bounded lag range (a tenth of the trajectory by default, truncated at the first zero crossing), because the longest lags average over very few time origins and are noise. The summary states which rule was applied and how many lags were used.

The tolerance time `t*` has a strong effect on the result (Laage and Hynes, *J. Phys. Chem. B* **112** (2008) 7697), which is why the scan field exists — a residence time quoted from a single `t*` should be treated with caution.

##### Equations

Every equation the tool evaluates — with what it represents, where it is used, and its reference — is documented in [Formulas and Methods (Mean Residence Time tool)](#formulas-and-methods-mean-residence-time-tool).

##### A worked example

This example measures how long a water hydrogen stays hydrogen-bonded to a particular solute oxygen, in a periodic box of water.

**The system.** A `TRAJEC.xyz` from CPMD with 430 atoms per frame and 15121 frames: a C₁₄ solute in water, so the composition is C 14, H 280, O 136. The acceptor of interest is atom 295, an oxygen. The simulation cell is 15 × 15 × 22 Å.

**Choosing the cutoff.** Run **Structural > Radial Distribution Function** first, with atom 295 as the shell centre and `H` as the target. Its first peak is the hydrogen bond; the first minimum after that peak is where the first shell ends. For this system:

```text
first O295-H peak : r = 1.905 A
first minimum     : r = 2.405 A   <- use this as the cutoff
```

**What to enter.**

| Tab | Field | Value | Why |
| --- | --- | --- | --- |
| header | XYZ trajectory file | `C:\...\TRAJEC.xyz` | |
| Setup | Δt between saved frames | `1.0` | set your real value; `1.0` reports times in frames |
| Setup | unit | `fs` | |
| Setup | Cutoff radius (Å) | `2.405` | the first g(r) minimum above |
| Setup | Tolerance frames | `0` | start unbridged, then scan |
| Setup | Cell lengths a b c (Å) | `15 15 22` | enables minimum image |
| Setup | Scan t* values | `0 2 5 10 20` | reveals the `t*` sensitivity |
| Setup | Output folder | *(blank)* | writes next to the trajectory |
| Definitions | Reference mode | Single atom | the shell is centred on one acceptor |
| Definitions | Reference definition | `295` | 1-based index |
| Definitions | Observed mode | Each atom in a list, tracked separately | pool statistics over all candidate donors |
| Definitions | Observed definition | `15,16,17,...,294` | every H in the box |

Then press **Run MRT**, and read the summary on the **Results** tab.

**Reading the result.** For this system the tool reports:

```text
Average occupancy <h>: 0.00734480
Complete residence events: 272
Censored events (clipped by the start or end of the run): 4

Mean residence time <T>: 95.62500000 frames  +/- 64.68269375
IMM survival integral:   658.81120365 frames
Intermittent time from R(t): 1187.16839960 frames
  (integration rule: zero_crossing, 1512 of 1512 lags used)

TOLERANCE (t*) SENSITIVITY
------------------------------------------------------------------------
 t* (frames)   events              <T>     IMM integral
           0      272      95.62500000     658.81120365
           2      265      87.46415094     696.09164845
           5      255      91.04705882     696.88557485
          10      237      98.56962025     727.93365504
          20      188     116.59574468     724.58083047
```

Four things to take from it:

- $\langle h\rangle \times 280 \approx 2.06$, so O295 accepts about two hydrogen bonds on average — a sensible number for a carbonyl or ether oxygen, and a good sanity check that the cutoff and the reference atom are right.
- $\langle T\rangle = 95.6$ frames with an uncertainty of $\pm 64.7$. The uncertainty is 68 % of the value, so this is an order-of-magnitude statement, not a precise one.
- The `t*` scan moves $\langle T\rangle$ from 95.6 to 116.6 frames. Quote the scan, not a single number.
- **`1512 of 1512 lags used` is a warning.** By §9 that means $R(t)$ never crossed zero, so 1187 frames is a lower bound that grows if you raise the lag ceiling — not a converged time constant. For this trajectory the intermittent value should not be quoted at all; $\langle T\rangle$ and $\tau_{\text{IMM}}$ are the defensible numbers.

**A common mistake.** If the reference atom is a water oxygen, its own two covalently bonded hydrogens sit at about 0.96 Å and are permanently inside a 2.4 Å shell. They will dominate the statistics and produce enormous residence times. Either choose a reference that has no bonded hydrogens (as atom 295 here does not), or exclude the bonded hydrogens from the observed list.

#### Legacy Mean Residence Time

The legacy MRT workflow is retained for compatibility with older analyses.

Inputs include:

- Shell inner radius.
- Shell outer radius.
- Simulation time step, in **a.u.** (the default, as in the CPMD input) or **fs** — chosen beside the field.
- Sampling interval (MD steps between two saved frames).
- **Tolerance**, in **frames** or **fs** — chosen beside the field (default `0`; blank means `0`).
- Atom labels to exclude.
- Atom labels at the shell center.
- Element symbol to investigate.
- `TRAJEC.xyz` trajectory file.

**Tolerance.** The tolerance is the largest number of consecutive *sampled*
frames an atom may spend outside the shell without that excursion being
counted as an exit; a gap of up to that many frames between two in-shell
frames is treated as inside. With the default `0`, every exit from the shell
is counted. A tolerance in fs is converted to frames with

- time between sampled frames: Δt (fs) = simulation time step (fs) × sampling interval (steps)
  (an a.u. time step is first multiplied by 0.02418884 fs/a.u.);
- tolerance (frames) = tolerance (fs) ÷ Δt, **rounded down**.

Rounding down is not an approximation: an excursion of *n* frames lasts
*n*·Δt and is tolerated when *n*·Δt ≤ the time given, so a time between two
multiples of Δt acts as the lower one (2.5 Δt → 2 frames), and *n* frames and
*n*·Δt fs always give identical results. A value in frames must be a whole
number. The tolerance applied is shown, in frames and fs, when the run starts
and in `mrt_summary.dat`.

**Exits and MRT.** An exit is an inside → outside step whose excursion lasts
longer than the tolerance; an excursion still under way at the end of the
trajectory counts only if it already lasts longer than the tolerance. For each
atom, `mrt_summary.dat` gives the total residence time, the number of exits and
mrt = total residence time ÷ exits (or the total residence time, a lower bound,
when the atom was never seen to leave). Up to version 0.9.5 this file counted
entries into the shell (`exchange`) instead; the two differ by at most one per
atom.

**Speed and progress.** *Browse* counts the frames from the file's line count,
and *MRT calculation* reads the trajectory once, parsing only the shell center
and the atoms of the chosen symbol (2–3 s for a 351 MB, 15 000-frame
trajectory of 430 atoms; the earlier two-pass reader took about 20 s, with
identical output files). Both run in the background: a spinner and a status line beside it show
what is being done (counting frames, reading frame *n* of *N*, applying the
tolerance, writing files) and the window stays responsive. Frame 1 is checked
line by line; in later frames the atom-count line of every frame and the lines
of the atoms that are analysed are checked, and every frame must list the atoms
in the same order, with the same elements, as frame 1. Lines of atoms that are
not analysed are not read.

Typical outputs:

- `mrt.dat`
- `mrt_total.dat`
- `mrt_summary.dat`

#### Autocorrelation Function

Use **Structural > Autocorrelation function** to calculate the velocity
autocorrelation function (VAF), the position autocorrelation function (PAF), and
an optional power spectrum (PSD) from a CPMD `TRAJECTORY` file.

Inputs:

- **Starting frame for VAF and PAF** — first trajectory frame to analyze.
- **Stop frame for VAF and PAF** — last trajectory frame to analyze.
- **Number of frames for each ACF** — the maximum correlation lag (window length,
  in frames).
- **Simulation time step** (atomic units) and **Sampling interval** (frames
  between stored steps) — used to build the physical time/frequency axes.
- `GEOMETRY.xyz` file (defines the number of atoms and element symbols).
- `TRAJECTORY` file (positions and velocities, atomic units).

Default values (whole trajectory):

- When you load the `TRAJECTORY` file, any of the three frame fields left blank is
  filled automatically so that the analysis uses the **entire trajectory**:
  - *Starting frame* defaults to **1** (first frame).
  - *Stop frame* defaults to the **last available frame**.
  - *Number of frames for each ACF* defaults to about **10% of the trajectory
    length**, capped at **2000** frames (floor 2). Using ~10% keeps roughly 90% of
    the frames available as time origins for averaging, which gives a statistically
    reliable estimate while bounding the computational cost and keeping the power
    spectrum's frequency resolution sensible on long trajectories.
- Any value you type yourself is preserved; only blank fields are auto-filled.
  Load `GEOMETRY.xyz` first, then `TRAJECTORY`, so the atom count is known when the
  frames are counted.
- If the `TRAJECTORY` was produced by a **restarted/continued** CPMD run it may
  contain markers such as `<<<<<<  NEW DATA  >>>>>>` (and blank lines). These are
  detected and skipped automatically, and the output box reports how many were
  found. Frames on either side of a restart are treated as one continuous
  trajectory.

How the functions are defined:

- The VAF/PAF are averaged over **all available time origins** (overlapping
  windows), include **lag 0**, and are **normalized so that C(0) = 1**.
- Positions are **mean-subtracted per atom**, so the PAF measures fluctuations
  about each atom's average position rather than absolute coordinates.

Options:

- Save generated `newTRAJEC.xyz` (positions converted from Bohr to Angstrom).
- Compute the **power spectrum** (vibrational density of states) from the VAF.

Power-spectrum options (used when the PSD switch is enabled):

- **Mass-weight the spectrum** — when checked, each atom's velocity is weighted by
  its atomic mass, producing the true vibrational density of states (VDOS) and the
  file `VDOS.dat`; when unchecked, an equal-weight power spectrum is written to
  `PSD.dat`. Masses are taken from the elements in `GEOMETRY.xyz`.
- **PSD window** — apodization applied before the FFT to suppress spectral leakage:
  `Welch` (default), `Hann`, or `None`.
- **Zero-pad x** — zero-padding factor (default `5`); higher values interpolate a
  smoother spectrum without changing peak positions.
- **Frequency unit** — `cm^-1` (default), `THz`, or `Hz`.
- **Partial VDOS groups** — an optional field to overlay *partial* spectra
  (per element or per atom selection) on the total, which helps **assign
  vibrational modes** (i.e. see which atoms contribute at which frequency). Enter
  one or more groups separated by `;`. Each group is a comma/space list of element
  symbols and/or 1-based atom indices or ranges, for example:
  - `H,O` — one partial for hydrogens together with oxygens;
  - `H; C; O` — three separate partials, one per element;
  - `1-8; 15` — atoms 1–8 as one partial and atom 15 as another.
  Each partial is written to `PSD_<group>.dat` (or `VDOS_<group>.dat` when
  mass-weighting is on) and drawn as a labelled curve on the power-spectrum plot,
  together with the `Total`. Each curve is normalized to its own maximum so the
  peak positions of every group are visible. Leave the field blank for the total
  spectrum only.

The power spectrum is computed by the Wiener–Khinchin route: the VAF is
mean-removed, mirrored about *t* = 0, windowed, zero-padded, Fourier transformed,
and squared; the one-sided spectrum is normalized so its strongest peak is 1. The
frequency axis is built from the **Simulation time step** and **Sampling interval**
you provide, so set those correctly for the frequency scale to be physical.

The PAF, VAF and power-spectrum figures open as **interactive** matplotlib windows
(zoom, pan, live cursor coordinates, and a save button in the toolbar), which makes
it easy to zoom into spectral peaks and read off frequencies. If the program is run
as a packaged executable where an interactive window is unavailable, the figures are
shown as static images instead.

##### Power Spectrum of a Solute Only (e.g. an Organic Molecule in a Water Box)

A common case is a simulation of an organic molecule (solute) in a box of water
(solvent) where you want the power spectrum of **only the solute**, without the
solvent contribution. Use the **Partial VDOS groups** field, keeping these points
in mind:

- **Do not select by element.** The solute and water usually share elements
  (H and O), so an element token such as `O` would also pick up every water oxygen.
  To isolate the solute you must select it by **atom index**.
- **Use the solute's atom-index range/list.** The indices come from the atom order
  in your `GEOMETRY.xyz` file (the same order as the `TRAJECTORY`). In most solvated
  setups the solute is written **first**, followed by the water molecules. Open
  `GEOMETRY.xyz`, count the atoms belonging to the organic molecule, and enter that
  range. For a solute that is the first 15 atoms, type:

  ```text
  1-15
  ```

  You may also give a list or several ranges, e.g. `1-12,14,17` or `1-15,20-22`.
  If the solute atoms are not contiguous, list their actual indices.
- **Read the result from the partial curve/file.** The group produces a partial
  spectrum built from only the selected atoms, so water contributes nothing to it.
  It is saved as `PSD_1-15.dat` (written as `PSD_1_15.dat`) — or `VDOS_1_15.dat`
  when mass-weighting is enabled — and drawn as the `1-15` curve on the plot. That
  file/curve is the solute-only power spectrum.
- **Ignore the `Total` curve.** The tool still computes and overlays the `Total`
  spectrum for the whole system (solute + water); simply use the solute partial and
  ignore `Total`. Because each curve is normalized to its own maximum, the large
  water background does not distort the solute spectrum.
- **Tip:** enable **Mass-weight the spectrum** for a physically proper vibrational
  density of states of the molecule.

Typical outputs:

- `PAF.dat` — columns: lag (frames), normalized PAF.
- `VAF.dat` — columns: lag (frames), normalized VAF.
- `PSD.dat` (or `VDOS.dat` when mass-weighting is on) when the power spectrum is
  enabled — columns: frequency (in the chosen unit), normalized power.
- `newTRAJEC.xyz` when the trajectory export option is enabled.

#### Single Solute Solvent Box

Use **Structural > Single Solute solvent box** to create a box containing one solute and inserted solvent molecules. The window has three tabs — **Setup**, **Construction** and **Summary** — above a fixed run bar (**Build Solvent Box**, **Cancel**, **Close**), status line and progress bar.

Setup tab inputs and options:

- Box lattice vectors (three positive lengths in Å).
- Target density.
- Maximum insertion attempts.
- van der Waals scaling factor.
- Extra wall padding or minimum distance.
- Solvent XYZ file.
- Solute XYZ file.
- Randomly rotate solvent.
- Calculate density.
- Periodic minimum-image clash detection.
- Include centered solute.
- **Charge and spin**: solute charge, charge per solvent molecule (both integers, blank = 0) and the **total spin multiplicity** (2S+1, required).
- **Pre-optimization**: *Pre-optimize after building*, the mode (*Rigid molecules* or *Flexible molecules*), the **target contact scale** (blank = 0.90, or the vdW scale if that is larger) and the maximum number of iterations (blank = 500).
- **Advanced > Packing engine**: *Automatic* (uses numba when it is installed), *NumPy cell list*, *Numba kernel* or *Legacy* (the original loop). All engines produce exactly the same box for the same random seed; they differ only in speed.

The **Construction** tab shows the box being built: the cell edges, the solute (orange, larger) and the solvent molecules (blue), farther atoms paler. It is updated a few times per second while solvent is inserted and during the pre-optimization, and finally shows the finished box. **Turn left/right** and **Tilt up/down** rotate the view by 15°, and **Zoom +**, **Zoom −** and **Reset view** change the zoom or restore the default view. Boxes with more than 4000 atoms show one dot per solvent molecule.

##### System properties

After the build the summary reports the total number of atoms, the net charge, the total number of electrons and the spin multiplicity:

- Net charge = solute charge + (number of inserted solvent molecules × solvent molecule charge). The number of solvent molecules is only known after packing, which is why charges are entered per species.
- Electrons = Σ Z (atomic numbers of all atoms) − net charge.
- The multiplicity is **never inferred**: an even electron count allows 1, 3, 5, … and an odd one 2, 4, 6, …, so the electron count alone does not decide it. A blank field stops the build and asks for it. After the build it is checked: the value must have the opposite parity to the electron count and be at most N_e + 1. An impossible value is marked *INCONSISTENT* in the summary and a dialog lists the allowed values; correct it on the Setup tab and press **Update properties** on the Summary tab — the box is not rebuilt.

##### Pre-optimization

Random insertion only guarantees `d ≥ vdW scale × (r_i + r_j)`, so with the usual 0.8 many contacts sit right at that tight limit. The optional pre-optimization pushes them further apart. It minimizes an overlap penalty `(d² − σ²)²` for every pair of atoms in *different* molecules closer than `σ = target contact scale × (r_i + r_j)` (plus the extra solute clearance for solute–solvent pairs).

The target must be larger than the packing vdW scale to have any effect, since packing already satisfies its own scale. The default is 0.90, or the vdW scale if that is larger. A final penalty of zero means every contact now satisfies the target. The minimizer aims 1% beyond σ, so contacts settle just outside the cutoff rather than a hair inside it.

- **Periodic boundaries**: with minimum-image clash detection on, every distance uses the minimum image, so molecules interact across the faces of the cell. The contact range must be shorter than half the smallest cell length. If it is not, the packed box is still written and a dialog explains why the pre-optimization was skipped. With the switch off, molecules are kept inside the box by soft walls at the wall padding.
- **Rigid molecules** (default): each molecule moves only as a whole (translation and rotation), so bond lengths and angles are exactly those of the input files.
- **Flexible molecules**: atoms move individually. Bonds (detected from covalent radii) and angles are held near their input values, so molecules can adjust slightly to their neighbours. There are no torsion terms.
- **Molecule-aware wrapping**: afterwards every molecule is translated as a whole so that its centre of mass lies inside the cell. Molecules are never split across a face, so atoms of a molecule near a face can lie slightly outside the cell. That is the usual convention for periodic codes, and it keeps every molecule whole.
- The summary reports the penalty, the number of clashing pairs and the closest intermolecular contact before and after, and the largest change of an intramolecular distance. It also confirms that the bond connectivity of every molecule is unchanged.

Typical outputs:

- `single_solute_solvent.txt`
- `single_solute_box_cmass.xyz`
- `single_solute_solvent_box.xyz` — the packed box, as before.
- `single_solute_solvent_box.txt` — the summary, now including system properties and pre-optimization results.
- `single_solute_solvent_box_preopt.xyz` — the pre-optimized, molecule-wrapped box, only when pre-optimization ran. Its comment line records the cell lengths.

If the target number of solvent molecules cannot be inserted, the summary reports a warning and the number of successful insertions. **Cancel** stops the run; no output files are written for a cancelled run.

#### Mixture of Two Solvent Box

Use **Structural > Mixture of two solvent box** to build a mixed-solvent box, optionally with a centered solute.

Inputs and options:

- Solvent A XYZ file.
- Solvent B XYZ file.
- Optional solute XYZ file.
- Mixture composition.
- Box/density and insertion settings.
- Random solvent rotation.
- Density calculation.
- Periodic minimum-image clash detection.
- Insert centered solute.
- **Charge and spin**: solute charge, charge per solvent A molecule and per solvent B molecule (integers, blank = 0) and the required **total spin multiplicity**. If the multiplicity is missing or impossible, the message names the field and the Setup tab scrolls to it.
- **Pre-optimization** and **Advanced > Packing engine**, as described for the [Single Solute Solvent Box](#single-solute-solvent-box).

The window uses the same **Setup / Construction / Summary** tabs, **Cancel** button and live construction view as the single-solute builder; solvent A is drawn blue, solvent B green and the solute orange. The net charge is solute charge + (inserted A molecules × A charge) + (inserted B molecules × B charge), and the electron count and multiplicity check work exactly as for the single-solute box. All packing engines give the same box for the same seed.

Outputs:

- `mixture_solvent_box.xyz` — the packed box, as before.
- `mixture_solvent_box.txt` — the summary with composition, insertion statistics, system properties and pre-optimization results.
- `mixture_solvent_box_preopt.xyz` — the pre-optimized, molecule-wrapped box, only when pre-optimization ran.

#### Shared-Wall Double Solvent Box

Use **Structural > Shared-wall double solvent box** to build two adjacent solvent regions that share one face. The module is intended for interface models, bilayer-like solvent arrangements, or any setup where two different solvent environments should occupy neighboring slabs in the same simulation box.

The generated structure is a single combined XYZ file. Box 1 occupies the lower z region, from `z = 0` to `z = c1`, and Box 2 occupies the upper z region, from `z = c1` to `z = c1 + c2`. The two boxes share the same `a` and `b` lattice dimensions, but each box has its own z height. No wall atoms or boundary markers are written at the shared interface.

How it works:

- Reads one or two solvent XYZ templates for each box.
- Optionally reads one solute XYZ template for each box.
- Recenters each molecule template, using center of mass when atomic masses are available.
- Estimates the target number of solvent molecules from the requested density, slab volume, solvent composition, and optional solute mass.
- Randomly inserts solvent molecules into Box 1 and then Box 2.
- Optionally randomizes solvent orientation before insertion.
- Keeps atoms inside the assigned z slab and applies the requested wall padding at the outside boundaries.
- Rejects placements that clash using van der Waals radii scaled by the selected vdW factor.
- Uses the full combined box for periodic minimum-image clash detection when that option is enabled.
- Reports insertion statistics and warns if the target number of molecules cannot be reached.

Shared geometry inputs:

- **Shared a b**: two lattice dimensions used by both boxes.
- **Box 1 c**: z height of the lower region.
- **Box 2 c**: z height of the upper region.
- **Periodic minimum-image clash detection**: checks overlaps across periodic boundaries using the full `a x b x (c1 + c2)` box.

Per-box inputs and options:

- Wall padding.
- Target density in g/cm^3.
- Solvent A and Solvent B composition percentages.
- Maximum insertion attempts.
- vdW scale.
- Optional random seed.
- Optional solute gap.
- Rotate solvent molecules during insertion.
- Include a centered solute.
- Solvent A, Solvent B, and optional solute XYZ files.

Charge, spin and pre-optimization:

- **Total spin multiplicity** is at the top of the Setup tab, in **Shared geometry**, because it applies to the whole cell (required, never inferred).
- **Charges** are in each box's own section: **Solvent A charge**, **Solvent B charge** (per molecule) and **Solute charge**. They are integers; blank means 0. The net charge is the sum over both boxes.
- The electron count and the multiplicity check work as for the [Single Solute Solvent Box](#single-solute-solvent-box). If the multiplicity is missing or impossible, the message names the field, and the Setup tab scrolls to it with the cursor in it.
- **Pre-optimization** works on the combined `a × b × (c1 + c2)` cell, adapted to the two-region geometry:
  - a flat-bottom restraint keeps each molecule's centre of mass inside the region it was packed into (Box 1: `0 ≤ z ≤ c1`, Box 2: `c1 ≤ z ≤ c1 + c2`), so molecules relax across contacts but the interface composition is preserved. A note in the summary reports any molecule that ends outside its region;
  - with minimum-image detection off, soft walls act at each box's wall padding on the x/y faces and on the outer z faces only — never at the shared interface;
  - the target contact scale defaults to 0.90 or the larger of the two boxes' vdW scales; the extra solute clearance is the larger of the two boxes' solute gaps.
- **Advanced > Packing engine**: as for the single-solute builder; Box 2 is checked against everything already placed in Box 1, and all engines give the same structure for the same seeds.

The window has **Setup / Construction / Summary** tabs; the Construction view shows Box 1 being filled, then Box 2, then the pre-optimization.

Typical outputs:

- `shared_wall_box.xyz` — the packed structure, as before. It is still shifted so that its lowest atom sits at `z = 0`.
- `shared_wall_box.txt` — the summary, now with system properties and pre-optimization results.
- `shared_wall_box_preopt.xyz` — the pre-optimized structure, only when pre-optimization ran. It is **not** shifted: it is in cell coordinates (Box 1 from `z = 0`, Box 2 from `z = c1`) with every molecule wrapped whole into the cell, and the cell lengths are on its comment line.

Example workflow:

1. Open **Structural > Shared-wall double solvent box**.
2. Enter shared lattice values, for example `30.0 30.0`, in **Shared a b**.
3. Enter `25.0` for **Box 1 c** and `25.0` for **Box 2 c**. The combined box will have `c = 50.0`, with the shared interface at `z = 25.0`.
4. For Box 1, set **Density** to `0.95`, **Solvent A (%)** to `100`, **Solvent B (%)** to `0`, and browse for a Box 1 Solvent A XYZ file such as water.
5. For Box 2, set **Density** to `0.80`, **Solvent A (%)** to `100`, **Solvent B (%)** to `0`, and browse for a Box 2 Solvent A XYZ file such as an organic solvent.
6. Leave **Rotate Box 1** and **Rotate Box 2** enabled unless you need fixed solvent orientations.
7. Enable **Box 1 solute** or **Box 2 solute** only if that slab should contain a centered solute, then browse for the corresponding solute XYZ file and set a solute gap if needed.
8. Click **Build Shared-Wall Box**.
9. Inspect the progress and summary panel. If the summary warns that not all solvent molecules were inserted, try a larger slab, lower density, smaller vdW scale, or more insertion attempts.

### Thermo Tools

#### Classical Rate Constant

Use **Thermo > Classical rate constant** to calculate rate constants from classical transition-state-style inputs.

Inputs include numeric parameters such as temperature range and activation-energy-related fields. The tool validates each field and asks for an output directory.

Typical outputs:

- `rate_constant_vs_temperature.dat`
- `LnK_vs_invT.dat`

### General Tools

#### 3D Molecular Viewer

Use **Tools > 3D Molecular Viewer** to inspect XYZ structures and trajectories.

**Window layout.** The controls are arranged as a fixed header, nine tabs, and a status line:

- **Header** (always visible, whichever tab you are on): the XYZ file field with **Browse**, plus the two main actions — **Display Molecule/Trajectory**, which opens the 3D window, and **Save Current Frame XYZ**.
- **Display** tab: projection, atom style, atom and bond scale, maximum bond length, rotation, and atom labels. The atom styles are Line, CPK (ball and stick), vdW (space-filling) and **Design style**, which draws the molecule as the Molecular Design Canvas does: equal flat discs in the canvas colours with a thin outline, single bond lines, the element symbol centred on each disc and the atom number beside it, in colours that stay readable on the current background. The *Atom numbers* and *Atomic symbols* switches are independent (either, both or neither) and are the same settings as the QTAIM tab's *Atom labels*.
- **Dynamics** tab: frame navigation, trajectory playback, and per-frame bond recomputation.
- **Measure** tab: distances, angles, dihedrals, and atom coordinates.
- **Design** tab: draw, import and pre-optimize a molecule (see [Molecular Design](#molecular-design)).
- **Cube** tab: Gaussian cube files as isosurfaces and slices, grid arithmetic, cube and image export (see [Cube tab](#cube-tab)).
- **ESP** tab: the electrostatic potential mapped on electron-density isosurfaces of a Gaussian `.fchk`, its extrema, and ESP/density cube export (see [Electrostatic Surface Potential (ESP)](#electrostatic-surface-potential-esp)).
- **QTAIM** tab: Bader's topology analysis of a Gaussian `.fchk` density — critical points, their properties, bond and ring paths, a text report and a 2D plane map (see [QTAIM Topology Analysis](#qtaim-topology-analysis)).
- **NTO** tab: natural transition orbitals of a CIS/TDA/TD excited state from a Gaussian `.fchk` and its output — pair weights, hole and particle isosurfaces, cube and NTO `.fchk` export — and the GaussView-style UV-Vis spectrum of all the excited states (see [Natural Transition Orbitals (NTO)](#natural-transition-orbitals-nto)).
- **Box & Performance** tab: periodic box settings and playback performance options.
- **Status line** (bottom): reports the last action, such as trajectory loading progress or the value you just applied.

**Loading a trajectory.** After **Browse**, the file is read in the background: a spinner turns beside the status line, which shows the frames read so far and the fraction of the file, and the window stays responsive. **Browse** is unavailable until loading ends. A 351 MB, 15 000-frame trajectory of 430 atoms loads in about 3 s (about 7 s with the earlier line-by-line reader, which read identical frames). Reading follows the same rules as before: blank or non-numeric lines where a frame should start are skipped, atom lines may have extra columns after x, y and z, and a frame cut short by the end of the file, with a line of fewer than four fields or with a coordinate that is not a number ends the reading, keeping the frames before it. Line endings may be `\n` or `\r\n`.

A typical session is: **Browse** for the file, adjust anything you need on the **Display** tab, then click **Display Molecule/Trajectory**. The 3D window opens alongside the control window, so you can keep changing settings while it is open.

**Numeric fields are applied by pressing Enter.** Type a value into a field — maximum bond length, atom or bond scale, box edge lengths, playback delay, frame step, frame number, or a rotation angle — and press Enter to apply it. There are no separate "Set" buttons.

##### Display tab

- **Projection**: `Orthographic` or `Perspective`.
- **Atom style**: `Line style`, `CPK style`, or `vdW style`.
- **Atom scale** and **Bond scale**: relative size of the drawn atoms and bonds.
- **Max bond length**: upper distance limit, in Å, for two atoms to be drawn as bonded. Type a value and press Enter to apply it to every atom pair. Leave it blank (and press Enter to clear a previous value) to use the default, which depends on the two elements: the sum of their covalent radii plus 0.45 Å (for example 1.97 Å for C–C, 1.42 Å for O–H), or 1.5 Å for an element with no known radius. Two atoms farther apart than the applicable limit are never drawn as bonded.
- **Rotation**: rotate the displayed system about X, Y, and Z (described below).
- **Atom numbers** and **Atomic symbols**: label every atom in the 3D view. Both are **off by default**, so the view starts uncluttered — click individual atoms to read their numbers (see the **Measure** tab's **None** type), or turn **Atom numbers** on to label the whole structure at once.
- **Clear selection**: deselect every atom you have clicked, remove any atom numbers you revealed by clicking, and clear the **Measure** tab's field and result.

##### Dynamics tab

- `|<`, `<<`, `>>`, `>|` and the slider move through the trajectory. The field to the left of the slider jumps to a specific frame when you press Enter; the counter to the right shows `current / last`.
- **Play** / **Pause** animates the trajectory.
- **Loop mode**: `Loop` restarts at the beginning, `Once` stops at the end, `Rock` reverses direction at each end.
- **Frame step**: how many frames each step or playback advance skips. Use `<` and `>`, or type a value and press Enter.
- **Delay (s)**: seconds between frames during playback. Larger values play more slowly.
- **Auto-zoom to fit each frame**: rescales the view for every frame. Turn it off for a steadier view during playback.
- **Recompute bonds for every frame**: re-detects the bonds from each frame's own interatomic distances (described below). Off by default, in which case the bonds found in the first frame are kept for the whole trajectory.

##### Measure tab

- The measurement type starts on **None**, which measures nothing: left-clicking an atom simply shows its number on the canvas, and clicking it again hides the number. Use it to find the indices you want before switching to a real measurement.
- Choose **Bond length**, **Bond angle**, **Dihedral angle**, or **Atom coordinates**, then select the atoms **either** by typing their labels into the field **or** by clicking them in the 3D window. Atom indices are 1-based. The value appears in the **Result** field and on the 3D canvas.
- **Measure** computes the value for whatever is currently in the field — use it after typing. When you select atoms by clicking, the measurement runs on its own as soon as enough atoms are selected. Pressing it while the type is **None** asks you to pick a measurement type first.
- **Clear** empties the field, deselects every clicked atom, removes any atom numbers you revealed by clicking, and removes the measurement from the canvas.
- Switching the type clears the marks from the previous mode, so numbers revealed under **None** disappear when you start a measurement, and a half-finished selection disappears when you go back to **None**.
- The **Result** field gives the fully labelled measurement, for example `Measurement: d(1,2) = 1.5000 Å`. On the 3D canvas only the number is drawn (`1.5000`, `120.307`), so the value stays readable against the structure. Distances are in Å, angles and dihedrals in degrees.

##### Box & Performance tab

- **Box a, b, c**: orthorhombic periodic box edge lengths in Å. Fill in all three, then press Enter.
- **Box center**: where the displayed box is placed (described below).
- **Show simulation box**: draw or hide the box.
- **Reduce overlays during playback**: skips some overlay drawing to keep long trajectories smooth.

The saved current frame includes the atoms from the displayed frame and a comment noting the frame number.

##### Cube tab

The **Cube** tab shows Gaussian cube files (`.cube`/`.cub`) in the 3D window: orbitals, densities, density differences, potentials, Fukui functions, or any other volumetric data on a grid.

**Loading.** **Load cube(s)...** reads one or more files at once. What is accepted:

- cubes in bohr or in Å (a negative voxel count on lines 4–6 means Å);
- skewed cells and a different spacing along each axis, as written by periodic codes;
- orbital cubes (negative atom count). Each orbital becomes its own entry, for example `phenol MO 24`.

The list shows every cube with its grid size, smallest and largest value, and integral. The integral is the sum of the values times the voxel volume: the electron count, for a density. Cubes made by a grid operation are marked `(unsaved)`.

**Showing cubes.** Each row ends in a **Show** checkbox:

- Tick it to draw that cube on the canvas; clear it to hide the cube. Hiding keeps the cube loaded with all its settings.
- Tick two or more to draw those cubes together, superimposed on the same canvas (for example two orbitals, or a density and a Fukui function).
- New cubes start ticked.
- **Remove selected**, the **Isosurface** settings and **Save selected cube...** act on the cube chosen in the Isosurface section's **Cube** list, not on a clicked row. The **Visible** switch there and that cube's **Show** box are the same setting.

**Display** reopens the 3D window if you closed it (or brings it to the front) and draws every ticked cube again, with the slice if it is on. Nothing is lost when the window is closed.

**Which molecule is shown.** The first cube you load decides it:

- With nothing displayed, the cube's own atoms are shown at once.
- If the viewer already shows the cube's own molecule (for example the one the ESP tab put there), it is kept as is, with no question.
- If a different structure is displayed, a dialog asks whether to replace it with the cube's atoms. **No** keeps your structure, and the volume is drawn in the cube's own coordinates.
- Later cubes never change the structure. If their atoms differ from the displayed cube's, the status line says so.

Two switches control what is drawn:

- **Show molecule** hides or shows the atoms and bonds. Hidden atoms cannot be clicked.
- **Show volumetric data** hides or shows every surface and the slice.

**Background.** **Background colour...** opens the Windows colour picker, starting on the colour in use. You can also type a colour name (Black, White, Grey, ...) or `#rrggbb` in the field beside the swatch and press Enter. The colour applies at once and changes nothing else: the same cubes stay shown, nothing is recomputed. It is the 3D window's background, so the ESP tab's **Canvas background** sets the same colour; the last one you choose wins.

**Isosurface.** Choose the cube, type the **Isovalue** and press Enter.

- For signed data (orbitals, differences), the positive surface is drawn at +isovalue and the negative one at −isovalue, each in its own colour.
- **Show** draws both surfaces, `+ only`, or `− only`. Cubes with no negative values (densities) start at `+ only`.
- Starting isovalues are 0.02 for orbitals and signed data and 0.0004 for densities, as in GaussView. **Potential (ESP) cubes**, recognised from the words "electrostatic potential" that cubegen and the ESP tab write in the header, start at ±0.02 a.u. with both signs.
- **Why ESP isosurfaces need larger isovalues.** The potential decays slowly: at the faces of a grid sized for the density it is still about ±0.01 a.u. A surface at a smaller value lies partly outside the grid and is cut open there, which shows as flat walls and torn edges. The status line then says *the surface … reaches the edge of the grid and is cut open there; raise the isovalue*.
- **Positive colour** and **Negative colour** take a name (Blue, Red, Green, Yellow, Orange, Cyan, Magenta, Purple, White, Grey) or `#rrggbb`. Press Enter to apply.
- **Opacity** below 1 makes the surface translucent.
- **Visible** hides one cube without removing it (the same setting as its **Show** box in the list).
- If an isovalue cuts nothing, the status line says so and gives the cube's value range.

**ESP cubes load with their density cube.** Loading an ESP cube on its own and looking at its isosurfaces shows **surfaces of constant potential** (V = ±isovalue). That is a different picture from the ESP tab's, which is the **electron-density surface** coloured by V. So, when you load `<name>_ESP.cube`:

- if `<name>_density.cube` is beside it on the same grid, it is loaded too. The ESP tab's **Save ESP cube...** writes both when **Also save the electron density cube** is on. The density surface (ρ = 0.0004) is shown coloured by the ESP, which is the ESP tab's picture. The ESP cube's own surfaces start unticked; tick its **Show** box to add them. The status line says *`<name>_ESP` shown on the `<name>_density` surface; tick its Show box for its own isosurfaces.*
- if there is no density cube, the ESP cube's own surfaces are shown, and the status line names the file to load for the ESP on the molecular surface.
- a density cube on a different grid is not loaded; the status line says why.
- choosing both files yourself loads each once and pairs them the same way.

**Colour by another cube.** **Colour by** colours the chosen cube's surface with another loaded cube's values instead of the solid colours: the usual way to look at an ESP is the ESP cube coloured on the density surface. Loading an ESP cube sets this up by itself (above). To set it up by hand, for example with the `_ESP.cube` and `_density.cube` written by the ESP tab:

1. Load both cubes, and untick the ESP cube's **Show** box (its own isosurfaces are not wanted).
2. Under **Isosurface**, choose the density cube in **Cube** and the ESP cube in **Colour by**.

- The other cube is interpolated (trilinear) at every surface point, so it does not need the same grid. Points it does not cover are drawn grey, and the status line gives their share.
- Colours follow GaussView and the ESP tab: red = most negative, green = 0, blue = most positive, over ± the largest |value| on the surface. A colour bar with the range and a caption is drawn at the **bottom** of the 3D window, so it never covers the ESP tab's bar at the top.
- With several coloured surfaces, the bar belongs to the first one.
- For benzene, the range matches the ESP tab's exact calculation within 0.2 % (±2.658e-2 against ±2.664e-2 a.u.).
- **None (solid colours)** goes back to the positive and negative colours. Removing the colouring cube does the same.

**Slice.** A coloured plane through one cube's grid.

- **Grid axis** `a`, `b` or `c` is a lattice axis of the cube. For an orthogonal grid these are x, y and z; for a skewed cell the plane follows the cell.
- The slider chooses the plane.
- **Colours** `Symmetric ±max|V|` runs blue (−), white (0), red (+) over ± the largest |value| in the plane. `Min to max` runs red, green, blue.
- **Export slice data...** writes one line per grid point: x, y, z (Å) and the value.
- The slice is opaque, so it hides whatever lies behind it.

**Grid operations.**

- **A + B**, **A − B** and **c × A** create a new cube named after the operation, unless you type a **Result name**.
- **A + B** and **A − B** need identical grids: the same number of points, the same origin and the same axes. Otherwise they are refused, with the reason.
- If A and B list different atoms, the result keeps A's atoms and the status line says so.
- Typical uses: a density difference ρ(N+1) − ρ(N), or a sign change with c = −1.

**Export.**

- **Save selected cube...** writes the cube in bohr, in the layout cubegen and Multiwfn use (6 values per line, `%13.5E`). An orbital is written as a plain cube, with its MO number on the comment line.
- **Export image...** saves the 3D window exactly as shown (PNG, TIFF or BMP).

**Speed.** Isosurfaces are built by marching tetrahedra on the cube's own grid. When `numba` is installed, a compiled kernel scans the grid, about 4× faster on a 200³ grid and with far less memory; without numba the result is identical. A surface is built once per isovalue, so changing a colour, the opacity or the slice does not rebuild it.

**Limits.**

- Translucent surfaces are not sorted against each other.
- The volume stays fixed while a trajectory plays.
- Grids with more than about 250 points per axis are slow to mesh and are held fully in memory.
- Slice colours are linear: a density's slice shows mostly the regions near the nuclei.

##### Selecting atoms by clicking in the 3D view

What a left-click does depends on the **Measure** tab's type. A click always means a click *without* dragging — dragging still rotates the view.

With the type on **None** (the default), clicking an atom shows its 1-based number next to it in cyan, and clicking the same atom again removes it. Nothing is measured and the **Measure** field stays empty. Click as many atoms as you like; **Clear** removes all the numbers at once.

With any real measurement type selected, you can select atoms directly in the 3D view instead of typing their labels. The atom turns yellow, gets a cyan `#1`, `#2`, … tag showing the order in which you picked it, and its 1-based index is written into the **Measure** tab's field. Typing and clicking are always both available: whatever ends up in the field is what gets measured.

**Click order is the measurement order.** Selecting atoms 4, 7 and 9 in that order for a **Bond angle** gives the field `4,7,9`, so atom 7 is the vertex — exactly as if you had typed it. The `#1`/`#2`/`#3` tags on the canvas tell you which atom is which.

Once you have selected as many atoms as the measurement needs (two for a bond length, three for an angle, four for a dihedral, one for atom coordinates), the measurement runs automatically and the result appears in the **Result** field and on the canvas.

Selecting further atoms keeps working without any clearing step:

- Left-click a selected atom to deselect it.
- Left-click a new atom when the selection is already full and the **oldest** pick drops off. For a bond length, this lets you walk a measurement along a chain of atoms one click at a time, getting a new value at every click.
- Switching the measurement type to one that needs fewer atoms keeps your most recent picks.
- **Clear** (Measure tab) or **Clear selection** (Display tab) deselects everything.

Because a click always writes into the field, it will overwrite part of an entry you typed by hand. If you are working from typed labels, click **Measure** before clicking anywhere in the 3D window.

Selections stay attached to the same atoms as you step through or play the trajectory, because atom ordering is constant across frames. Note that **Reduce overlays during playback** (Box & Performance tab) hides all atom labels while it is on, including the index and pick-order labels.

The rotation controls on the **Display** tab let you rotate the displayed molecular system around the X, Y, and Z axes without changing the coordinates stored in the loaded structure or trajectory. Enter rotation increments in degrees, then either press Enter in one of the angle fields or click **Apply angles** to add those increments to the current molecular orientation. Click **Reset rotation** to return the displayed molecular orientation to zero rotation. The `<` and `>` buttons beside each axis start continuous rotation in the negative or positive direction; click the same arrow again to stop it, or click another arrow to switch to that axis and direction.

When the periodic box is shown, it rotates with the molecular system. The **Box center** selector on the **Box & Performance** tab controls where the displayed orthorhombic box is placed:

- **Geometric center**: centers the box on the molecular system's geometric center.
- **Bottom at z=0**: centers the box in the a-b plane while placing the bottom face of the box at `z = 0`.

**Recompute bonds for every frame** (Dynamics tab) controls whether connectivity follows the trajectory:

- **Off** (default): bonds are calculated once from the first frame and that same connectivity is drawn in every frame. This gives a steady picture during ordinary vibrations and rotations, but a bond that breaks later in the run stays drawn.
- **On**: for every frame, each interatomic distance is compared with the applicable limit — the **Max bond length** you typed on the Display tab, or the element-pair default when that field is blank — and a bond is drawn only while the distance is within it. Bonds therefore disappear when they break and appear when they form. Each frame's result is cached, so revisiting a frame is instant.

Example: in a trajectory where two atoms separate during a dissociation event, turn **Recompute bonds for every frame** on. As you step through or play the frames, the bond disappears in the first frame where the distance exceeds the limit. To make the bond break at a specific distance, type that distance in **Max bond length** and press Enter. Changing the limit takes effect immediately on the current frame.

##### Molecular Design

The **Design** tab of the 3D Molecular Viewer lets you draw a molecule with the
mouse — or import one from a structure file — pre-optimize its geometry, and
send the result straight into the 3D viewer, or save it as an XYZ, Gaussian
input, MOL, SDF, MOL2 or PDB file. It is meant for building a starting
structure quickly, without leaving the toolkit.

###### Opening and using the two windows

The Design tab itself is a launcher. It has two buttons:

- **Open periodic table** — opens the **Periodic Table** window, where you
  choose the element to draw with.
- **Open design canvas** — opens the **Molecular Design Canvas** window, where
  you draw the molecule and run the optimization.

Open both and place them side by side. They are ordinary windows: move, resize
or close either one at any time. Pressing a launcher button again brings the
matching window back to the front rather than opening a second copy, and the
canvas window has its own **Periodic table** button for the same purpose.

**Closing a window never loses your work.** The molecule is held by the tool,
not by the window, so you can close the canvas, reopen it, and find your
structure exactly as you left it. Closing the Molecular Viewer window closes
both design windows with it.

###### Selecting elements from the periodic table

Click any element in the periodic-table window to select it. The table shows
all 118 elements in the usual layout, tinted by block (s, p, d and f), with the
selected element highlighted in yellow.

The selection takes effect in the canvas window **immediately** — there is
nothing to confirm or apply. The current element is shown in three places at
once: under the periodic table, above the canvas, and on the Design tab. Carbon
is selected when you start.

###### Placing atoms on the canvas

Click an empty part of the canvas to place an atom of the selected element. The
atom is drawn as a coloured circle with its symbol inside, using the same
colours as the 3D viewer.

To build a molecule, select an element, place the atoms you need, then select
the next element and continue. Clicking an existing atom does **not** place a
second atom on top of it; it picks the atom for a measurement instead (see
[Measuring](#measuring-bond-lengths-angles-and-dihedrals)).

**Hydrogens are not added while you draw.** Draw them yourself, or press **Add
hydrogens** once the heavy atoms and bond types are in place (see
[Adding hydrogens automatically](#adding-hydrogens-automatically)).

###### Creating and modifying bonds

To create a bond, **press the left mouse button on one atom, drag to a second
atom, and release**. A dashed guide line follows the pointer while you drag. A
new bond is always a single bond.

To **grow** the molecule, press on an atom, drag to **empty space** and
release: a new atom of the selected element is added, already bonded to the
one you started from. It is placed in the direction you dragged, at the typical
single-bond length for the two elements (1.52 Å for C–C), so the drawing starts
out close to its relaxed geometry.

To change a bond, **click it**. Each click steps to the next bond type:

```
single  ->  double  ->  triple  ->  resonance  ->  single
```

Clicking a bond that is at the end of its cycle returns it to a single bond, so
the same click both raises and lowers the bond type and nothing is ever out of
reach. The status line names the new type after every click.

Which types are available depends on the two elements involved:

| Bond | Cycle |
|------|-------|
| C–C, C–N, N–N | single → double → triple → resonance → single |
| C–O, C–S | single → double → resonance → single |
| O–H, C–F, C–Cl | stays single (the status line explains why) |

The four types are drawn like this:

| Type | Appearance |
|------|------------|
| single | one solid line |
| double | two parallel solid lines |
| triple | three parallel solid lines |
| resonance | one solid line with a dashed line beside it |

**Resonance bonds** are for delocalized bonding: the aromatic bonds of benzene,
the two equivalent C–O bonds of a carboxylate, the C–N bond of an amide. A
resonance bond counts as one and a half bonds towards an atom's valence and
relaxes to a length between a single and a double bond. Draw benzene as six
carbons in a ring with six resonance bonds, and it optimizes to six equal
1.40 Å bonds in a flat ring.

To **delete** a bond, right-click it (without dragging). To delete an atom,
right-click the atom — every bond attached to it goes too. **Undo** steps back
through the last 50 edits, including **Clear**.

###### Rotating, panning and zooming: building in 3D

The canvas is three-dimensional: every atom has a depth, and you can turn and
move the view at any point while you build. The mouse works as in the 3D
Molecular Viewer's window (the one the QTAIM tab draws in): left-drag rotates,
right-drag moves, the wheel zooms.

| Action | Gesture |
|--------|---------|
| Rotate the molecule about its centre | Left-drag an empty part of the canvas, or a bond |
| Pan (move the view) | Right-drag anywhere on the canvas |
| Zoom in / out | Mouse wheel (zooms about the pointer), or **Zoom +** / **Zoom −** (about the canvas centre) |
| Back to the original view | **Reset view** |

A left-drag that starts **on an atom** still draws a bond or grows a new atom,
and a right-click **without dragging** still deletes, so the view gestures
never get in the way of editing.

Atoms further away are drawn smaller and paler, and bonds and atoms in front
are painted over those behind. Rotating, panning and zooming **only change how
you look at the molecule**: they are not edits, they are not on the Undo list,
and they do not discard an optimized geometry.

Where new atoms go in depth:

- A **clicked** atom is placed on the screen-facing plane through the centre
  of the molecule.
- A **grown** atom is placed at the same depth as the atom it grows from.

So, to build out of the drawing plane — a tetrahedral carbon, a substituent
above a ring — **rotate the view first**, then grow the new atom from the side
you can now see. When two atoms overlap on screen, a click picks the one in
front.

###### Running the pre-optimization

**Pre-optimize geometry** relaxes the drawing into a realistic three-dimensional
structure. The drawn positions are the starting point: a drawing that is still
flat is first lifted slightly off its plane, while one you have already built
in 3D starts exactly as drawn. The geometry is relaxed with a compact force
field: harmonic bond stretching with reference lengths from the covalent radii
adjusted for each bond type; angle bending at the VSEPR angle for each atom,
counting lone pairs, which is what makes water come out bent rather than
linear; a twofold torsion term on double and resonance bonds, which keeps
alkenes and aromatic rings flat; and a repulsive term that stops distant parts
of the molecule from passing through each other.

When the run finishes, a dialog reports the final energy, the gradient norm,
the iteration count and any warnings. A run that does not converge is reported
with the minimizer's own message, and its partially relaxed structure is still
available to send on.

**The canvas is updated with the optimized geometry.** The molecule relaxes in
place, without moving or turning, so you can rotate the view to inspect the
result right away. The update counts as one edit: **Undo** brings your
original drawing back. Further edits build on the relaxed structure, and
pre-optimizing again starts from it, which makes it easy to refine a molecule
in steps: draw, optimize, add more atoms, optimize again. If you change the
drawing while an optimization is still running, that result is discarded
rather than overwriting your edit, and the status line tells you to run it
again.

This is a **pre-optimizer**: it produces a reasonable starting geometry for a
real CPMD, cp.x or ORCA calculation, not a converged result. The energies it
reports are in relative units and mean nothing outside a single run. It is
tuned for main-group, mostly organic molecules; geometries with a steric number
of six (octahedral) come out distorted.

###### Sending the structure to the 3D viewer

**Send to 3D viewer** transfers the optimized structure and **opens the 3D
viewer window automatically**, bringing it to the front if it is already open.
The connectivity you drew is transferred with it, so the bonds you made are the
bonds you see. Once it is there, **Save Current Frame XYZ** writes the structure
to an XYZ file.

The 3D view draws every bond the same way whatever its type. Bond order shows
up in the geometry instead, where it belongs: a resonance C–C bond arrives
visibly shorter than a single one.

If the 3D window cannot be opened — a missing OpenGL driver, a blocked window,
a start that takes too long — a dialog reports the reason. The structure has
still been transferred, so **Display Molecule/Trajectory** and **Save Current
Frame XYZ** remain available.

Editing the drawing again discards the previous optimization, so the structure
you send is always the one you last optimized. Rotating, panning or zooming the
view does not.

###### Importing and saving structure files

The **File** row of the canvas window holds **Import structure...**, the
save-format list and **Save Frame**.

**Import structure...** replaces the drawing with the structure in a file, and
you can then edit it with every tool above. The import counts as one edit, so
**Undo** brings the previous drawing back. The molecule is centred on the
canvas, the view is reset, and it is zoomed out if that is needed to show all
of it. A dialog reports the number of atoms and bonds read, plus anything the
reader had to decide for you.

| Format | Extensions | Import | Save | Bonds kept in the file |
|--------|------------|--------|------|------------------------|
| XYZ | `.xyz` | yes | yes (default) | none |
| Gaussian input | `.gjf` `.com` `.gau` | yes | yes | bond orders, in a `Geom=Connectivity` section |
| Gaussian formatted checkpoint | `.fchk` `.fch` | yes | no | none |
| MDL Molfile V2000 | `.mol` | yes | yes | bond orders (resonance = aromatic, type 4) |
| MDL SD file | `.sdf` `.sd` | yes (first record) | yes | bond orders |
| Tripos MOL2 | `.mol2` | yes | yes | bond orders (resonance = `ar`) |
| PDB | `.pdb` `.ent` | yes | yes | connectivity (`CONECT`), no bond orders |

**Bonds on import.** When a file has no bonds (XYZ, fchk, a Gaussian input
without `Geom=Connectivity`, a MOL or MOL2 file with an empty bond block, a PDB
without `CONECT` records), they are detected from the interatomic distances
with the 3D viewer's rule (a pair is bonded within the sum of its covalent
radii plus 0.45 Å), and **every detected bond is single**. Before
pre-optimizing, click the double, triple and resonance bonds to set their type,
or they relax to single-bond lengths. A PDB file's `CONECT` records are used
where they exist; atoms they never mention, such as a protein's standard
residues, get distance-detected bonds.

**What the readers understand.** XYZ: element symbols, numbered labels (`C1`)
or atomic numbers, and extra columns after x y z. A comment line reading
`coordinates in bohr` (as the Coordinate Converter writes) is converted to
Ångström. Gaussian input: `Units=Bohr`, freeze codes (`C 0 x y z`), fragment and
ONIOM labels (`C(Fragment=1)`, `C-CA--0.1`) and comma separators. Dummy and
ghost centres (`X`, `Bq`, `Du`, `LP`) are skipped. Only the first structure is
read from an XYZ trajectory, a multi-record SD file, a multi-model PDB or a
linked (`--Link1--`) Gaussian job, and the import dialog says so.

**Not supported**, each refused with a message saying why: Gaussian Z-matrix
input (re-save it with Cartesian coordinates), inputs that read the geometry
from a checkpoint (`Geom=Check`; import the `.fchk` instead), MOL V3000,
Gaussian `.log` output, binary `.chk` files (convert them with `formchk`), CIF
and HyperChem `.hin`. Files of more than 1000 atoms are refused too; open large
systems in the 3D viewer with **Display Molecule/Trajectory**.

**Malformed files** are refused with the file name, the format and the line at
fault, for example `'water.xyz' could not be imported as XYZ: line 4: 'zero' is
not a number.` The drawing is left unchanged. Atoms closer than 0.4 Å are
refused as well, since that almost always means duplicated atoms or a file in
bohr.

**Save Frame** writes the structure on the canvas: the optimized geometry after
a pre-optimization, the drawing otherwise, in Ångström and centred on the
origin. The view (rotation, pan, zoom) is not applied. Pick the format in the
list beside the button; **XYZ** is the default. The save dialog suggests
`designed_molecule` with that format's extension. If you type a name ending in
a *different* supported extension (`frame.pdb`), that format is written
instead, so the content always matches the name. Any other name gets the
selected format.

The **Gaussian input** is a template to check before running. Its route is
`#P B3LYP/6-31G(d) Opt Geom=Connectivity`, its `%Chk` is named after the file,
and it assumes charge 0 and the lowest multiplicity the electron count allows
(1 for an even count, 2 for an odd one). The confirmation dialog repeats this.
A MOL or SD file holds at most 999 atoms and 999 bonds, and a PDB file at most
99 999 atoms; a structure that does not fit is reported, and no file is written.

The canvas window has a side panel to the right of the drawing with four
groups, **Measure**, **Frozen parameters**, **Molecule** and **Canvas**, used
by the features below.

###### Selecting and moving one molecule

When the canvas holds several separate molecules (fragments that share no
bond), you can move one of them without touching the others.

| Action | Gesture |
|--------|---------|
| Select a molecule | **Shift + click** any of its atoms or bonds |
| Clear the selection | **Shift + click** empty space, or **Deselect** in the Molecule panel |
| Move the selected molecule | **Ctrl + drag** (left button), anywhere on the canvas |
| Turn the selected molecule | **Ctrl + right-drag**, about the molecule's own centre |

The selected molecule is drawn with an **orange halo** around its atoms and
orange bonds, and the Molecule panel gives its size. Moving and turning are
rigid, so the molecule's bond lengths, angles and dihedrals never change, and
every other molecule stays exactly where it is. The molecule follows the
pointer on screen whatever the view rotation and zoom, and it turns the same
way a left-drag on empty space turns the whole view. Each drag is one **Undo** step.

A move is refused, with the reason in the status line, when:

- nothing is selected;
- the selection has since been bonded to another molecule, because moving only
  part of a molecule would stretch that bond. Shift + click it again to
  reselect the whole molecule;
- a frozen parameter links an atom of the selection to an atom of another
  molecule, because the move would break it. Remove the frozen parameter
  first.

Without Ctrl or Shift every gesture works exactly as before. **Ctrl +
right-click never deletes anything.**

###### Measuring bond lengths, angles and dihedrals

**Click atoms in order** (a click without dragging):

| Atoms picked | Measurement | Unit |
|--------------|-------------|------|
| 2 | bond length (any distance) | Å |
| 3 | bond angle, with the **second** atom as the vertex | ° |
| 4 | dihedral angle about the bond between the 2nd and 3rd atoms | ° |

Picked atoms get a magenta ring and their pick order (**#1**…**#4**), joined by
a dotted line in that order. The value is shown in the top-left corner of the
canvas and in the **Measure** panel, for example `Bond length 1-2 (C-C): 1.5200
Å`. The dihedral follows the IUPAC sign convention: looking along atom 2 → atom
3, it is positive when atom 1 must turn clockwise to eclipse atom 4. A
dihedral through three atoms in a line has no defined value and is shown as
undefined.

Clicking a picked atom again removes it from the pick list; picking a fifth
atom starts a new measurement with it; **Clear** removes all picks. Measuring
is not an edit and is not on the Undo list. A picked atom that is deleted
simply drops out of the measurement.

###### Setting a bond length, angle or dihedral

Pick the atoms (2, 3 or 4, as for measuring) and press **Edit value...**. The
dialog shows the parameter, the picked atoms, the current value with its unit,
and a field for the new value (Å for a length, degrees for an angle). Type the
value and press **Apply** or Enter.

One side of the molecule then moves **rigidly**, the way a model kit behaves:

| Parameter | What moves |
|-----------|------------|
| Bond length 1-2 | everything on atom 2's side of the bond, along the bond |
| Bond angle 1-2-3 | everything on atom 3's side, turned about atom 2 in the plane of the three atoms |
| Dihedral 1-2-3-4 | everything on atom 4's side of the 2-3 bond, turned about that bond |

All other bond lengths and angles stay as they were. The edit is one **Undo**
step. Tick **Freeze at the new value during pre-optimization** to freeze the
parameter as you set it. A parameter that is already frozen is re-frozen at the
new value.

The edit is refused, and the reason is shown in red inside the dialog (which
stays open so you can correct it), when:

- the value is not a number, or is out of range. A length must be 0.3–20 Å,
  and an angle greater than 0° and at most 180°. A dihedral may be any angle
  and is wrapped into (−180°, 180°];
- the bond to be cut is in a **ring**, so no side can move on its own. Ring
  geometry changes only through pre-optimization;
- the picked atoms are not bonded in a chain (1-2, 2-3, 3-4). The exception is
  a length between two atoms of **different molecules**, which moves the second
  molecule as a whole;
- a dihedral passes through three atoms in a line;
- the move would change **another frozen parameter**. Remove that one first.

If one of these applies before you type anything, for example a ring bond,
**Edit value...** reports it at once instead of opening the dialog.

###### Freezing parameters during pre-optimization

Press **Freeze** to hold the picked length, angle or dihedral at its current
value during **Pre-optimize geometry**. Alternatively, tick the freeze box in
the edit dialog to hold it at a typed value. Every frozen parameter is listed
under **Frozen parameters** (kind, atom numbers, value). Select a row and press
**Remove**, or press **Remove all**. Freezing and removing are on the Undo
list, and they do not discard a finished optimization, since no atom moves.

**How it works.** Frozen parameters are held **exactly at every step** of the
optimization, not merely pulled towards their values. Each trial geometry is
projected back onto the frozen values before its energy is evaluated, and the
minimizer only moves in directions that leave them unchanged. The result
dialog reports how many were held and the largest deviation (typically about
1e-10 Å or °). Everything that is not frozen relaxes as usual. Without frozen
parameters the pre-optimizer runs exactly as before.

**Limits, reported rather than ignored.** If any frozen parameter cannot be
held, the pre-optimization does not run, and a dialog names the parameter and
the reason:

- an angle frozen within 1° of 0° or 180°. The constraint has no defined
  direction there; the angle term already keeps a linear centre straight;
- a dihedral through three atoms in a line;
- frozen values that contradict each other, for example three lengths that no
  triangle can have. The dialog reports the largest violation; remove or
  change one of them.

Deleting an atom removes every frozen parameter that uses it.

###### Adding hydrogens automatically

**Add hydrogens** (Molecule panel) completes each atom's valence with hydrogen
atoms, according to its bonds. If a molecule is selected, only that molecule
is completed; otherwise every molecule is. The valences used are: C 4, N 3,
O 2, B 3, Si 4, halogens 1, S 2 (also 4 and 6), and P 3 (also 5). A resonance
bond counts 1.5, so an aromatic CH carbon gets one hydrogen.

Existing hydrogens are bonds like any other, so a complete atom gets nothing,
and pressing the button twice adds no duplicates. Atoms whose state the
drawing does not settle are **listed in a dialog and left unchanged**:

| Case | Example | Why nothing is added |
|------|---------|----------------------|
| More bonds than the valence | N with 4 bonds | probably a charge (ammonium), which the canvas does not record |
| Bonds between two valences | S with 3 bonds, P with 4 | ambiguous oxidation state |
| Half a bond over resonance | carboxylate O (1.5) | H or negative charge cannot be told apart |
| Aromatic N with two resonance bonds | pyridine / pyrrole | pyridine-like (no H) or pyrrole-like (N–H)? Add the H yourself for N–H |
| No valence rule | Fe, Na, noble gases | no reliable rule |

New hydrogens are placed at the VSEPR positions around each atom at the
pre-optimizer's reference X–H length (1.07 Å for C–H). Alkenes are kept
planar, and the first hydrogen on a tetrahedral centre is placed anti to a
neighbour, giving a zigzag. **Pre-optimize afterwards** for the final geometry.
The addition is one Undo step.

###### Replacing an atom's element

To change an existing atom into another element, for example a carbon into a
nitrogen to turn benzene into pyridine:

1. Select the new element in the **Periodic Table** window.
2. Tick **Click an atom to replace it with N** at the top right of the canvas
   window. The box always names the element currently selected.
3. Click the atom (or atoms) to change.

The atom keeps everything except its element:

- its **position**, so nothing moves;
- its **number**, so measurements, saved files and the numbering below stay
  the same;
- its **bonds and their types**;
- any **measurement pick**, **frozen parameter** or **molecule selection**
  that includes it.

Each replacement is one **Undo** step. Choose another element in the periodic
table at any time and keep clicking. Clicking an atom that already has the
selected element changes nothing, and the status line says so.

The bonds are **kept exactly as drawn**, even when the new element would not
normally form them. The status line then reports, without changing anything:

| Note | Example | What to do |
|------|---------|------------|
| `X now has n bonds where m is typical.` | a carbon with three bonds replaced by O | delete or change a bond, if you did not mean it |
| `The double bond to atom k is above the single bond X usually forms; click it to change it.` | a C=C carbon replaced by H | click the bond until it is single |
| `X has no covalent radius, so the structure cannot be pre-optimized.` | a very heavy element | choose another element before pre-optimizing |

While replace mode is on, **clicking an atom replaces it instead of picking it
for a measurement**. Every other gesture is unchanged: clicking empty space
places an atom, dragging bonds and grows atoms, Shift and Ctrl gestures select
and move molecules. Untick the box to pick atoms for measurements again.

###### Showing atom numbers

Tick **Atom numbers** in the View row (beside **Reset view**) to show every
atom's number beside it, in blue. These are the same 1-based numbers used by
measurements (`Bond length 1-2`), by the Frozen parameters list, in the
validation messages, and as the atom order of saved files. Deleting an atom
renumbers the ones after it, and the labels follow at once.

Each number is placed beside its own atom where it is readable: it does not
cover any atom, any bond, another number, a measurement tag (#1–#4), or the
measurement text in the corner. When several free spots exist, the one on the
side away from the rest of the molecule is preferred, so ring atoms are
numbered outside the ring and a chain end beyond the chain. If an atom is so
crowded that no spot is entirely free, the least obstructed one is used.

The numbers are laid out again whenever the drawing or the view changes. While
you drag (pan, rotate, or move a molecule), each number moves with its atom,
so the canvas stays smooth even for a large molecule, and the numbers are
rearranged when you release the mouse. Atom numbers are a view setting: they
are not an edit and do not appear in saved files.

###### Canvas background

The **Background** list in the Canvas panel switches the canvas colour between
presets (off-white by default, white, light grey, pale blue, cream, dark grey,
black). Alternatively, type a colour as `#rrggbb` in **Custom** and press
Enter. On a dark background, bonds, outlines and the measurement text are drawn
light so they stay visible. The background is a view setting: it is not an
edit, and it does not affect saved files.

###### Bringing a molecule back from the 3D viewer

Press **From 3D viewer** in the canvas window's File row, or **Bring the 3D
molecule to the design canvas** on the Design tab. Either copies the frame
currently shown in the 3D viewer to the canvas for further editing, keeping:

- the **elements** and **coordinates** (unchanged, apart from centring);
- the **bonds the 3D viewer draws** for that frame, judged by its current bond
  settings.

The 3D viewer stores no bond orders, so its bonds arrive as **single bonds**,
and the dialog reminds you to set double and resonance bonds before
pre-optimizing. The exception is a structure that was sent to the viewer from
this canvas and not changed since: its drawn bond orders are restored.

**Your drawing is never overwritten.** If the canvas already holds atoms, you
are asked first:

- **Yes** adds the 3D molecule *beside* the drawing, 2 Å to its right. Nothing
  is lost, and the view is reset to show both.
- **No** cancels, and nothing changes. To replace the drawing instead, press
  **Clear** on the canvas first (Undo brings it back), then bring the molecule
  again.

On an empty canvas the molecule is loaded straight away. Either way the load
is one Undo step, and the 1000-atom limit of the canvas applies.

###### Mouse and keyboard reference

| Gesture | Where | Result |
|---------|-------|--------|
| Click | empty space | place an atom of the selected element |
| Click | atom | pick it for a measurement (again: unpick) |
| Click | atom, with **Click an atom to replace it** ticked | give it the element selected in the periodic table |
| Click | bond | step its type: single → double → triple → resonance → single |
| Drag | atom → atom | bond them |
| Drag | atom → empty space | grow a bonded atom |
| Drag | empty space or bond | rotate the view |
| Right-drag | anywhere | pan (move) the view |
| Mouse wheel | anywhere | zoom in / out about the pointer |
| Right-click | atom or bond | delete it |
| **Shift + click** | atom or bond | select its molecule |
| **Shift + click** | empty space | clear the selection |
| **Ctrl + drag** | anywhere | move the selected molecule |
| **Ctrl + right-drag** | anywhere | turn the selected molecule |
| **Enter** | edit dialog / Custom colour | apply the value |

###### Validation messages

The structure is checked before every optimization. Some problems stop the run,
because there is no sensible geometry to produce; others are only reported.

**These stop the run** (an error dialog names the atom responsible):

| Message | What it means | What to do |
|---------|---------------|------------|
| `The sketch has no atoms.` | The canvas is empty. | Draw something first. |
| `Atom N (X) is not bonded to anything.` | A lone atom sits beside a molecule — almost always a stray click. | Bond it, or right-click it to delete it. |
| `No covalent radius is known for X …` | The element has no radius in the program's tables, so no bond length can be assigned. | Use a different element; the tables cover the main group, the transition metals, the lanthanides and the early actinides. |

**These are reported but do not stop the run:**

| Message | What it means | What to do |
|---------|---------------|------------|
| `Unusual valence on atom N (X): 5 bonds where 4 is typical.` | The atom has more bonds than the element normally forms. | Nothing, if you meant it — unusual valences are chemically real. Otherwise remove a bond. |
| `Isolated resonance bond between atom N (X) and atom M (Y) …` | A resonance bond with no other resonance bond at either end. Delocalization over one isolated bond has no meaning. | Give it a resonance neighbour, or click it round to a single or double bond. |
| `The sketch contains N disconnected fragments …` | More than one separate molecule on the canvas. | Nothing, if you meant it — they are optimized together in one structure. |

Other failures are reported in a dialog as they happen: an optimization that
hits its iteration limit, or a 3D viewer that cannot open, both name the
underlying cause and leave the structure available.

#### Electrostatic Surface Potential (ESP)

Use the **ESP** tab of the **3D Molecular Viewer** (Tools > 3D Molecular Viewer)
to compute the molecular electrostatic potential of a Gaussian wavefunction on one or more
electron-density isosurfaces, view it colour-coded in 3D the way GaussView 6
draws a mapped ESP, and list the maxima and minima of the potential on each
surface, following the quantitative molecular surface analysis of Multiwfn.
The equations and algorithms are documented in
[Formulas and Methods (Electrostatic Surface Potential tool)](#formulas-and-methods-electrostatic-surface-potential-tool).

**What it is for.** The ESP on a molecular surface shows where an approaching
charge feels the molecule as positive (blue: electron-poor regions such as acidic
H atoms or σ-holes) or negative (red: lone pairs, π clouds). Its surface maxima
and minima are the standard descriptors of hydrogen- and halogen-bond donor and
acceptor sites and of electrophilic/nucleophilic attack.

**Input.** A Gaussian formatted checkpoint, `.fchk` or `.fch`, made from the
binary checkpoint with Gaussian's `formchk` utility:

```powershell
formchk name.chk name.fchk
```

The file must contain the basis set and the `Total SCF Density` section, which
every Gaussian SCF job writes. That one density is used for both the surface
and the potential. Open-shell (unrestricted) files are accepted: the ESP depends
only on the total (α + β) density.

**Two different quantities.** Keep them apart when reading the results:

| Quantity | What it is | Unit |
|---|---|---|
| Electron-density **isovalue** ρ<sub>iso</sub> | Defines *where* the surface is: the surface is the set of points where ρ(**r**) = ρ<sub>iso</sub>. | a.u. = e/bohr³ |
| **ESP** V(**r**) | The property *mapped onto* that surface: the colours, the extrema and every value in the output. | a.u. = Hartree/e, or kcal/mol |

The tab, the results text, the colour bar and the output files all name the
isovalue as "rho … (electron density isovalue)" and the potential as "ESP".

##### Tab layout

The ESP tab is one scrolling page, top to bottom: **Input** (file and output
folder), the surface definition and the mapped property, the buttons **Run**,
**Show surface** and **Clear** with a status line and a progress bar,
**Results**, **ESP cube**, and the text box (the help before a run, the results
summary after it). The surface is drawn in the viewer's own 3D window, together
with the molecule, so the other tabs (Display, Measure, Cube, ...) work on it too.

**Input and setup**

| Control | Meaning | Default / range |
|---|---|---|
| Formatted checkpoint | The `.fchk`/`.fch` file (**Browse** filters on both extensions). | required |
| Output folder | Where the extrema file is written. | blank = the folder of the `.fchk` |
| Isovalue(s) ρ (a.u. = e/bohr³) | One or several density isovalues, separated by commas, semicolons or spaces, e.g. `0.0004, 0.001, 0.002`. Order is kept, duplicates are dropped. | blank = **0.000400** (GaussView's default); each value must lie in (0, 1) |
| Grid spacing (bohr) | Spacing of the cubic density grid the surfaces are built from. Smaller gives smoother surfaces and more surface points, and takes longer. | blank = **0.25** (Multiwfn's default for the ESP); 0.05–0.6 |
| Density | Read-only: always `Total SCF Density`. | — |
| ESP engine | **Built-in exact** (analytic integrals in this program) or **Gaussian cubegen** (Gaussian's own program). | Built-in exact |
| cubegen executable | Path of `cubegen.exe`; pre-filled when one is found in `%GAUSS_EXEDIR%`, on `PATH`, in `C:\G09W` or in `C:\G16W`. Used only by the cubegen engine. | auto-detected |

**Results**

| Control | Meaning |
|---|---|
| Surface shown | Which isovalue's surface the 3D window displays (filled after a run). |
| Colour-bar unit | `a.u.` or `kcal/mol` for the colour-bar labels. |
| Surface | `Opaque` or `Transparent (50 %)`; **Show extrema markers** on/off. |
| Canvas background | `Black` (default), `White`, `Dark grey` or `Light grey`, for the viewer's 3D window. On a light background the colour-bar labels and caption are drawn in black so they stay readable. |
| Export image... | Saves the 3D canvas as an image file (see below). |
| Text box | The help text before a run; the results summary after it. |

Once the surface is shown, changing any Results control updates the 3D window
immediately.

**Which molecule is shown.** **Show surface** puts the `.fchk` molecule in the
viewer. If the viewer shows nothing, or already shows that molecule, this
happens without a question. If it shows a different structure (a trajectory,
a designed molecule, a cube's atoms), you are asked once per run: **Yes**
replaces it with the `.fchk` molecule, **No** keeps it and draws the surface in
the `.fchk` coordinates.

**Exporting the image.** Open the 3D window with **Show surface**, arrange the
view, then press **Export image...**. The frame is captured the moment you press
the button, before the save dialog appears, so the file holds exactly what the
canvas shows: the same orientation, zoom, framing, surface colours, markers,
colour bar and background, at the window's own pixel size (resize the 3D window
for a larger image). Only the canvas is saved — no window title bar or border —
and the view is not reset or changed. The suggested name is
`<name>_ESP_rho<isovalue>.png`; PNG, TIFF and BMP are accepted (all lossless —
JPEG is refused because it would alter the colours). If the 3D window is
minimised, restore it and export again.

**ESP cube.** After a run, **Save ESP cube...** writes the electrostatic
potential V(**r**) on the run's density grid (the same box and spacing the
surfaces were built from) as a Gaussian cube file, in a.u. (Hartree/e), with
the run's engine. A save dialog chooses the folder and the name; it opens in
the output folder with the name `<name>_ESP.cube` (for `Eteno.fch`:
`Eteno_ESP.cube`). With **Also save the electron density cube (same grid)** on
(the default), the density ρ(**r**) used for the surfaces is written beside it
as `<name>_density.cube`: the pair that GaussView or VMD use to draw an
ESP-mapped density surface, and that the Cube tab can load. The potential is
computed at every grid point, so this takes a few seconds for a small molecule
(ethylene, 134 000 points, about 1 s with numba) and longer for a large one
(benzene, 256 000 points, about 14 s); the progress bar follows it. V is
infinite on a nucleus: a grid point that falls exactly on one is evaluated
0.001 bohr away and the status line says so, so the file never holds an
infinite value. If cubegen fails, the built-in engine is used and the status
line says so.

##### Steps

1. Open **Tools > 3D Molecular Viewer** and choose the **ESP** tab.
   **Browse** to the `.fchk`/`.fch` file. Optionally choose an output folder.
2. Type the isovalue(s) or leave the field blank for
   0.000400; leave the grid spacing blank unless you need a finer surface.
3. Keep **Built-in exact**, or choose **Gaussian cubegen** if Gaussian is
   installed. If cubegen is missing or fails, the built-in engine is used for
   the rest of the run and the summary says so in a `NOTE:` line.
4. Press **Run**. The work runs in the background; the progress bar and status
   line report each stage (reading, density grid, and for each isovalue: surface,
   ESP, extrema). When it finishes, the text box shows the summary and the
   extrema file is written automatically.
5. Press **Show surface** to draw the surface in the viewer's 3D window (it
   opens the window, or brings it to the front). Rotate with the left mouse
   button, zoom with the wheel. Pick the surface, unit, opacity, markers and
   canvas background under *Results*, and press **Export image...** to save the
   canvas as shown.
6. Optionally press **Save ESP cube...** to write the ESP (and density) cube.
7. **Clear** discards the results, removes the surface from the 3D window and
   restores the help text.

##### Results summary

For every surface the summary lists the number of surface points (vertices),
the surface area in Å², the minimum and maximum ESP on that surface in a.u. and
kcal/mol, and a table of every extremum: row number, minimum/maximum,
global/local, ESP in a.u. and kcal/mol, x y z in Å and the nearest atom. The
header repeats the file, the molecule (atoms, electrons, charge, multiplicity),
the density, the ESP engine, the grid spacing and the run time.

##### Output file

`<name>_ESP_extrema.txt` (for `Eteno.fch`: `Eteno_ESP_extrema.txt`), one file
for all isovalues of a run, overwritten by the next run. A commented header
records the source file, density, surface definition, mapped property, ESP
engine, grid spacing, vertex-merge distance, extremum criterion, the
Hartree→kcal/mol factor, the coordinate unit, the molecule, and one
`# surface …` line per isovalue (vertices, triangles, area in Å², minimum and
maximum ESP in a.u.). Then one row per extremum:

```text
# row rho_isovalue_au type scope esp_au esp_kcal_mol x_A y_A z_A nearest_atom
1 0.000400 minimum global -0.02096722 -13.1571 -2.245522 -0.062908 0.009717 C1
3 0.000400 maximum global 0.01166272 7.3185 -0.010724 -2.212172 -2.112729 H5
```

| Column | Meaning |
|---|---|
| `row` | Running number over the whole file; the same number labels the marker in the 3D window. |
| `rho_isovalue_au` | Density isovalue of the surface the extremum lies on (a.u.). |
| `type` | `minimum` or `maximum`. |
| `scope` | `global` (lowest/highest point of that surface) or `local`. Each surface has exactly one global minimum and one global maximum. |
| `esp_au`, `esp_kcal_mol` | ESP at the extremum, in Hartree/e and kcal/mol (1 a.u. = 627.5094740631 kcal/mol). |
| `x_A y_A z_A` | Position in Å, in the coordinate frame stored in the `.fchk` (`Current cartesian coordinates`). |
| `nearest_atom` | Element and 1-based number of the closest nucleus. |

Within a surface, minima are listed first (global, then local from lowest to
highest), then maxima (global, then local from highest to lowest). Load the
numeric columns with `numpy.loadtxt(..., comments="#", usecols=(0, 1, 4, 5, 6, 7, 8))`.

##### The 3D window

The colour bar runs across the top, as in GaussView: red is the most negative
ESP, green zero, blue the most positive, on a linear hue scale over the
symmetric range ±max|V| of the surface shown. The end labels use GaussView's
style (`-2.097e-2 … 2.097e-2`; `-1.316e1 … 1.316e1` in kcal/mol), and the caption
underneath states the unit and the density isovalue of the surface. The surface
is drawn opaque with a specular highlight on a black background. Cyan spheres
mark minima and orange spheres maxima (the larger sphere is the global one), each
labelled with its row number in the extrema file. This is the 3D Molecular
Viewer's own window, with its usual mouse controls.

##### Example: ethylene (B3LYP/6-31G)

| ρ<sub>iso</sub> (a.u.) | Surface points | Area (Å²) | Minima (π cloud, above and below the plane) | Maxima (one by each H) |
|---|---|---|---|---|
| 0.000400 | 6655 | 88.84 | −0.020967 a.u. (−13.16 kcal/mol) | +0.01166 a.u. (+7.32 kcal/mol) |
| 0.001000 | 5587 | 74.67 | −0.025772 a.u. (−16.17 kcal/mol) | +0.01755 a.u. (+11.02 kcal/mol) |

GaussView 6.1 shows a colour bar of ±2.072e-2 a.u. for the same file at 0.0004;
this tool shows ±2.097e-2. GaussView interpolates the ESP from a coarse cube
grid, while this tool evaluates it exactly at every surface point. The built-in
and cubegen engines agree to better than 2×10⁻⁷ a.u. at every point. Each run
takes about 1.5 s.

##### Implemented capabilities and limitations

Implemented:

- Gaussian `.fchk`/`.fch` with s, p, SP, d, f and g shells, Cartesian (6D/10F)
  or spherical (5D/7F/9G) — verified against Gaussian 09 `cubegen` for each kind.
- Closed- and open-shell SCF wavefunctions (HF or DFT); effective core
  potentials (see below).
- Several isovalues in one run, each with its own surface, extrema and range.
- Exact ESP (built-in) or Gaussian `cubegen`, with automatic fallback.
- Optional `numba` acceleration of the built-in engine (about 5× faster;
  identical results), used automatically when installed.

Limitations:

- **Density:** only `Total SCF Density` is used. For MP2, CI, CC or TD-DFT jobs
  the correlated or excited-state density is *not* read; the surface and ESP are
  those of the SCF reference.
- **Basis:** functions above g (l > 4) are rejected with a message.
- **Effective core potentials:** the stored density lacks the core electrons and
  the nuclear charges are the effective ones, exactly as in Multiwfn. This is a
  good approximation on a molecular surface, far outside the core, but not close
  to the heavy nucleus.
- **No periodic boundary conditions**; isolated molecules or clusters only.
- **Extremum positions** are surface-point positions, not interpolated or
  refined, so their precision is about one grid spacing (0.25 bohr ≈ 0.13 Å by
  default). Symmetry-equivalent extrema can differ by ~10⁻⁵ a.u. for the same
  reason (ethylene at 0.001: −0.025772 vs −0.025753 a.u.).
- **Colour range** is always the symmetric ±max|V| of the shown surface; a
  fixed or user-set range is not available, so colours are not directly
  comparable between molecules.
- **Marker labels** hidden behind the opaque surface are not drawn; rotate the
  molecule or switch to *Transparent* to read them.
- **Memory and time** grow with molecule size: the whole density grid is kept in
  memory, and vertex elimination is a pure-Python loop, the slowest step for
  very large surfaces.
- Not implemented (possible future additions): Multiwfn's surface statistics
  (positive/negative areas, average and variance of V, Π, ν), export of the
  coloured surface (e.g. `.pqr`/`.cube`), interpolation of extremum positions,
  a user-defined colour range, and h functions.

#### QTAIM Topology Analysis

Use the **QTAIM** tab of the **3D Molecular Viewer** (Tools > 3D Molecular
Viewer) for Bader's *Quantum Theory of Atoms in Molecules* analysis of the
electron density of a Gaussian wavefunction. It finds every critical point (CP)
of ρ(r), computes the topological and energy-density properties at each one,
traces the bond paths (and, in cages, the ring paths), writes a text report,
draws the result in the viewer's 3D window and plots a Multiwfn-style 2D map of
any plane. The algorithms and every default parameter are Multiwfn's (main
function 2, *Topology analysis*), and the results agree with Multiwfn 2026.7.15
to the digits it prints (see *Validation* below). No external program is
needed at run time. The theory is summarised in
[Formulas and Methods (QTAIM tab)](#formulas-and-methods-qtaim-tab).

##### Tab layout

- **Input**: the `.fchk`/`.fch` file and an optional output folder (blank = the
  folder of the `.fchk`). The density analysed is always `Total SCF Density`.
- **Critical point search – starting points**: one switch per Multiwfn search
  mode, all on by default: *nuclei* (mode 2), *atom pairs* (3, the midpoint of
  each pair), *triangles* (4) and *tetrahedra* (5, the centres of three and four
  atoms). Modes 2 and 3 find the nuclear and bond points; 4 and 5 add the ring
  and cage points.
- **Advanced (blank = Multiwfn default)**: every parameter of Multiwfn's
  *Set CP searching parameters* and *Set path generating parameters* menus; the
  field shows the default when blank.

  | Field | Default | Meaning |
  |---|---|---|
  | Maximum iterations | 120 | Newton iterations per starting point |
  | Step-size scale factor | 1 | multiplies every Newton step; 0.5 helps when a BCP is missed |
  | Gradient-norm convergence | 1e-6 | \|∇ρ\| required at a CP (waived for nuclear starts) |
  | Displacement convergence | 1e-7 | length of the last Newton step (bohr) |
  | Minimum distance between CPs | 0.03 bohr | a point closer to a known CP is a duplicate |
  | Skip atoms beyond vdW sum × | 1.5 | combinations of atoms farther apart are not used as starting points |
  | Singular Hessian if \|det\| below | 5e-22 | the iteration stops (far from the molecule) |
  | Trust radius | none | when set, longer Newton steps are shortened to it (bohr) |
  | Maximum points per path | 451 | a longer path is dropped |
  | Path step | 0.03 bohr | spacing of the path points |
  | Stop within this of a CP | 0.05 bohr | a path ends when this close to its target |

- **Run / Display / Clear**, a progress bar and a status line. *Run* writes
  `<name>_QTAIM.txt` automatically. *Display* shows the result in the 3D window
  (opening it if needed). The molecule of the `.fchk` is loaded silently when
  the viewer shows nothing or the same molecule; otherwise you are asked first.
- **Results**: the counts and the Poincaré–Hopf check, then a table of every CP
  (row, type, assignment, ρ, ∇²ρ, ellipticity; the columns are re-fitted to the
  table's width whenever it changes). Selecting a row highlights that
  CP in white and prints its position, ρ, ∇²ρ, G, V and H in the status line.
  Below the table, switches show or hide each CP type, the bond paths, the ring
  paths and the CP numbers.
- **Canvas (3D window)**: the look of the 3D window, with the Molecular Design
  Canvas's choices. *Background* offers the canvas presets (Off-white, White,
  Light grey, Pale blue, Cream, Dark grey, Black); *Custom background* takes any
  `#rrggbb` (press Enter). *Molecule style* is Line, CPK, vdW or **Design style**
  (the canvas look). While the style is CPK or Design, the same menu also lists
  **Increase atom size** and **Decrease atom size**: each choice makes the atoms
  1.2× larger or smaller (between 0.2× and 5× the normal size) and the menu
  returns to the style. Smaller atoms uncover critical points and paths near
  the nuclei; larger ones make the atom labels easier to read. *Atom labels* has
  independent *Element symbols* and *Atom numbers* switches. Style, atom size
  and labels are the viewer's own settings, so the Display tab (its *Atom
  scale* field) follows them and vice versa. In Design style the (3,−1), (3,+1)
  and (3,+3) points and their numbers are painted over the atom discs (a bond
  point of an X–H bond lies inside the H disc), a (3,−3) is covered by its own
  atom and its number is left out (it equals the atom number), and CP numbers
  sit at the lower right of each point while atom numbers sit at the upper
  right of each disc. CP numbers are cyan on a dark background and dark teal on
  a light one. **Export 3D image...** saves the 3D window exactly as shown
  (PNG, TIFF or BMP; the frame is taken when you press the button, before the
  file dialog opens).
- **2D plane map**: the plane (*Three atoms*, or the *XY*, *XZ* or *YZ* plane at
  an offset along its normal), the three plane atoms (blank = the first three
  atoms that are not collinear), the margin around the projected atoms (2 Å),
  the grid points per axis (200), and switches for the density contours, the
  gradient lines and the interbasin paths. **Plot 2D map** opens the map in the
  interactive plot viewer (zoom, pan, save) and saves its data as
  `<name>_QTAIM_map.json`, which **Plotter** (*JSON plot type*) reopens later.
  **Export 2D image...** saves the last plotted map, drawn exactly as in the
  plot viewer, as PNG, SVG, PDF or TIFF (200 dpi for the bitmap formats).
- The report text, below everything.

##### Step-by-step procedure

**A. Prepare the wavefunction (Gaussian)**

1. Choose the geometry. QTAIM describes the density of *one* structure, so use
   an optimized geometry (or the exact structure you want to study). For a
   complex, optimize the complex, not the separate molecules.
2. Choose the level of theory. Any HF or DFT method works; use a basis with
   polarization functions (6-31G(d) or better) — properties at bond points of
   polar bonds and hydrogen bonds change noticeably without them. All-electron
   basis sets are preferred: with an ECP the core density is missing (see the
   limitations below).
3. Keep the checkpoint file: put `%chk=name.chk` in the input. The density
   read is the SCF density (`Total SCF Density`); for MP2, CCSD or TD-DFT jobs
   the correlated density is not used.
4. After the job ends normally, convert the checkpoint:
   `formchk name.chk name.fchk` (Gaussian's `formchk` utility; `.fch` works
   too).

**B. Load and run**

5. Open **Tools > 3D Molecular Viewer** and choose the **QTAIM** tab.
6. **Formatted checkpoint**: Browse to `name.fchk`. **Output folder**: leave
   it blank to write next to the `.fchk`, or Browse to another folder.
7. **Critical point search – starting points**: leave all four switches on
   (nuclei, atom pairs, triangles, tetrahedra). For a large molecule with no
   rings you may turn off *triangles* and *tetrahedra* to save time.
8. **Advanced**: leave every field blank for Multiwfn's defaults.
9. Click **Run**. The progress bar follows the search, then the paths; small
   molecules take well under a second with `numba` installed.

**C. Check that the topology is complete**

10. Read the summary line under **Results**, e.g.
    `(3,-3) 6   (3,-1) 5   (3,+1) 0   (3,+3) 0     Poincare-Hopf = 1 (satisfied)`.
    - The number of (3,−3) points should equal the number of atoms (more means
      a non-nuclear attractor, shown as `NNA`).
    - Every expected bond should have a (3,−1) point, and every ring a (3,+1).
    - *NOT satisfied* means points are missing. Try, one at a time, and run
      again: **Step-size scale factor** 0.5; **Trust radius** 0.1 bohr; a
      larger **Maximum iterations** (e.g. 300); for points far from the atoms
      a smaller **Singular Hessian** threshold (e.g. 1e-25). A satisfied rule
      is necessary, not sufficient: also check the 3D view.
11. Read the `NOTE:` lines at the end of the report text (incomplete bond
    paths, ECP atoms, open shell).

**D. Read the results**

12. The **table** lists every point: row, type, assignment (`O1`, `O1-H2`,
    `ring`, `cage`), ρ, ∇²ρ and ellipticity. Click a row: the point is
    highlighted in white in the 3D view and its position, ρ, ∇²ρ, G, V and H
    appear in the status line.
13. The **report** `name_QTAIM.txt` (also shown at the bottom of the tab) has
    every property of every point and the bond-path lengths. How to interpret
    the numbers is explained in [Example 9](#example-9-qtaim-analysis-of-a-hydrogen-bond-water-dimer)
    and in [Theoretical Background](#theoretical-background-the-quantum-theory-of-atoms-in-molecules).

**E. Look at the result in 3D**

14. Click **Display**. The `.fchk` molecule is loaded into the viewer (you are
    asked first if the viewer shows another structure) and the 3D window opens.
15. Use the switches under the table to show or hide each point type, the bond
    paths, the ring paths and the CP numbers.
16. In **Canvas (3D window)** choose the **Background**, a **Molecule style**
    (*Design style* draws the molecule as the Molecular Design Canvas does) and
    the **Atom labels** (element symbols, atom numbers, both or neither). In CPK
    and Design style, choose **Increase atom size** / **Decrease atom size** in
    the same menu, as often as needed, until the atoms neither hide the
    critical points nor look too small.
17. Rotate and zoom with the mouse, then click **Export 3D image...** to save
    the window as shown (PNG, TIFF or BMP).

**F. Make a 2D plane map**

18. Choose the **Plane**:
    - *Three atoms*: type three atom numbers (e.g. `1, 2, 3`), or leave the
      field blank for the first three atoms that are not on one line. For a
      planar molecule any three of its atoms give its plane; for a hydrogen
      bond choose the donor, the H and the acceptor.
    - *XY*, *XZ* or *YZ*: a Cartesian plane, moved along its normal by
      **Offset** (Å).
19. Adjust **Margin** (Å around the atoms) and **Grid points** (more = smoother
    contours, slower), and the switches for density contours, gradient lines
    and interbasin paths.
20. Click **Plot 2D map**. The map opens in the interactive plot viewer (zoom,
    pan, save) and is saved as `name_QTAIM_map.json` (reopen it later with
    **Plotter**, *JSON plot type*). Only points and paths within 0.5 bohr of the
    plane are drawn.
21. **Export 2D image...** saves the map as PNG, SVG, PDF or TIFF.

**G. Start again**

22. **Clear** removes the results, the 3D points and paths and the map; change
    the file or the settings and click **Run** again. Each run rewrites
    `name_QTAIM.txt`.

##### Output file

`<name>_QTAIM.txt` follows the layout of the Fukui report: a commented header
(file, method/basis, density, electron count, search modes, every parameter,
engine, units), then these sections:

- `CRITICAL POINT SUMMARY` — counts per type and the Poincaré–Hopf sum.
- `CRITICAL POINTS` — row, type, x y z (Å) and assignment: the atom of a
  (3,−3) (`O1`; `NNA` for a non-nuclear attractor), the atom pair of a (3,−1)
  (`O1-H2`, from its bond path), `ring` or `cage`.
- `DENSITY PROPERTIES AT CRITICAL POINTS (a.u.)` — ρ, ∇²ρ, λ₁, λ₂, λ₃,
  ellipticity and η.
- `ENERGY DENSITIES AT CRITICAL POINTS (a.u.)` — G, K, V, H, |V|/G and H/ρ.
- `BOND PATHS` — for each BCP its two half-path lengths (to each atom), their
  total, the distance between the two (3,−3) points and the ratio total /
  distance (1 for a straight path, larger for a curved one).
- `RING PATHS` — (3,+1) → (3,+3) paths and their lengths (cages only).
- `NOTE:` lines (missing CPs, incomplete paths, ECP atoms, open shell) and
  `Written:` lines.

Positions and lengths are in Å; every property is in atomic units. CP rows
are numbered as in the table and the 3D view: (3,−3) by atom, then (3,−1),
(3,+1), (3,+3).

##### The 3D window and the 2D map

In 3D, the CPs are small spheres in Multiwfn's colours — magenta (3,−3),
orange (3,−1), yellow (3,+1), green (3,+3) — with cyan CP numbers; bond paths
are brown lines and ring paths green.

The 2D map uses Multiwfn's plane-map conventions: ρ contours at 0.001, 0.002,
0.004, 0.008, 0.01, … a.u. (Bader's 1, 2, 4, 8 × 10ⁿ series), grey gradient
lines leaving each nucleus, deep-blue interbasin paths (the traces of the
zero-flux surfaces) leaving each bond point, bold brown bond paths, and the
CPs as dots — brown (3,−3), blue (3,−1), orange (3,+1), green (3,+3). CPs and
path points more than 0.5 bohr from the plane are not drawn. Nuclear points are
labelled with their atom, the others with their CP number.

##### Validation

Every fixture below gives the same CPs, types, atom pairs and paths as
Multiwfn 2026.7.15 run with modes 2, 3, 4 and 5 and its default parameters:

| Molecule (B3LYP) | CPs (3,−3)/(3,−1)/(3,+1)/(3,+3) | max position difference | max property difference (relative) | max path-length difference |
|---|---|---|---|---|
| water 6-31G(d) | 3 / 2 / 0 / 0 | 1×10⁻¹⁰ bohr | 2×10⁻⁹ | 2×10⁻⁶ bohr |
| ethylene 6-31G | 6 / 5 / 0 / 0 | 6×10⁻¹⁰ bohr | 3×10⁻⁹ | 3×10⁻⁶ bohr |
| benzene 6-31G(d) | 12 / 12 / 1 / 0 | 1×10⁻⁹ bohr | 4×10⁻⁹ | 2×10⁻⁶ bohr |
| methyl radical UB3LYP/6-31G(d) | 4 / 3 / 0 / 0 | 1×10⁻⁹ bohr | 6×10⁻⁹ | 5×10⁻⁶ bohr |
| water dimer 6-31G(d) | 6 / 5 / 0 / 0 | 4×10⁻¹⁰ bohr | 3×10⁻⁹ | 4×10⁻⁶ bohr |
| bicyclo[1.1.1]pentane 6-31G(d) | 13 / 14 / 3 / 1 | 2×10⁻⁹ bohr | 7×10⁻⁹ | 5×10⁻⁶ bohr |

The property differences are at the level of Multiwfn's printed digits (ten
significant figures). With `numba` installed the run and the map are 4–46×
faster (benzene: run 0.06 s instead of 1.7 s, map 0.9 s instead of 14 s).

##### Implemented capabilities and limitations

Implemented: Gaussian `.fchk`/`.fch` with s–g shells (Cartesian or spherical),
closed and open shell; Multiwfn search modes 2–5 with all their parameters;
classification and Poincaré–Hopf check; ρ, ∇²ρ, Hessian eigenvalues,
ellipticity, η, G, K, V, H; bond paths and ring paths; text report; 3D view;
2D plane map with contours, gradient lines and interbasin paths; optional
`numba` acceleration.

Limitations:

- **Density**: only `Total SCF Density` (no correlated, excited-state or spin
  densities).
- **Effective core potentials**: the `.fchk` holds no core density (Multiwfn
  adds one from its EDF library; this tool does not), so the (3,−3) points of
  ECP atoms and properties close to them are unreliable. The run still
  completes and the report says so in a `NOTE:`.
- **Not implemented**: 3D interbasin surfaces, basin integration (AIM charges
  and volumes), user starting points (mode 1) and random-sphere search (mode 6),
  Shannon aromaticity and ring-normal curvature, ESP or other functions at the
  CPs, `.wfn`/`.wfx` input, periodic systems.

#### Natural Transition Orbitals (NTO)

Use the **NTO** tab of the **3D Molecular Viewer** (Tools > 3D Molecular Viewer)
to describe an electronic excitation by *natural transition orbitals*. A CIS,
TDA or TD excitation usually mixes several occupied → virtual MO pairs; the NTOs
(Martin, *J. Chem. Phys.* **2003**, 118, 4775) rotate the occupied orbitals among
themselves and the virtual orbitals among themselves so that one or a few
**hole → particle** pairs carry the whole excitation, each with a **weight**
(weight × 100 = percent of the excitation). The definition, conventions and
output layout are Multiwfn's (main function 18, subfunction 6), and the results
agree with Multiwfn 2026.7.15 to the digits it writes (see *Validation*). No
external program is needed at run time. The method is in
[Formulas and Methods (NTO tab)](#formulas-and-methods-nto-tab).

##### What you need

Two files **of the same Gaussian job**:

- the formatted checkpoint (`formchk name.chk name.fchk`): basis set and MOs;
- the Gaussian output (`.log`/`.out`) of the CIS, TDA or TD calculation, which
  holds the configuration coefficients. **Add `IOp(9/40=4)` to the route**:
  without it Gaussian prints only coefficients above 0.1, the coefficients no
  longer sum to 0.5, and the NTOs are approximate (the report says so).

Example route: `# PBE1PBE/6-31G* TD(NStates=5) IOp(9/40=4)`. An excited-state
optimisation output also works (its last excitation block, at the final
geometry, is used). Multiwfn's plain-text format (`Excited State  n  mult  eV`
headers followed by `i -> a  coefficient` lines; closed-shell coefficients
normalised to 0.5) is accepted in place of the Gaussian output.

**Load** checks that the two files belong together: same atoms and geometry
(1e-4 Å), same basis-function and electron counts, both restricted or both
unrestricted, every orbital number in range. A mismatch is refused with a
message naming it.

##### Tab layout

- **Input**: formatted checkpoint (choosing it fills in `name.log` or `name.out`
  beside it when the excitation field is empty), excitation file, output folder
  (blank = the folder of the `.fchk`), **Load**.
- **Excited state**: the list of states (number, label, eV, nm, f), the *Peak
  half-width (eV)* of the spectrum, **Run**, **Plot Spectrum**, **Clear**,
  status line and progress bar.
- **NTO pairs**: a one-line summary and the table of pairs with weight ≥ 0.0001
  (pair, spin, weight, %, cumulative %). Clicking a row selects that pair.
- **Display**: pair, **Show** (*Hole*, *Particle*, *Hole and particle*),
  isovalue (default 0.02 a.u.), grid spacing (default 0.2 bohr), the four
  isosurface colours (Hole + blue, Hole − red, Particle + green, Particle −
  yellow; a colour name or `#rrggbb`, Enter applies), opacity and **Display**.
  Colour, opacity and Show changes redraw at once without recomputing; a new
  pair, isovalue or spacing recomputes the orbitals (a fraction of a second for
  a molecule of uracil's size).
- **Export**: **Export cubes** (the chosen pair, or every listed pair with *All
  listed pairs*), **Export NTO .fchk**, **Send to Cube tab** (the displayed
  pair's two orbitals as unsaved Cube-tab entries, for slices, grid arithmetic
  or saving under another name) and **Export 3D image...** (the 3D window as
  shown; captured first, then the file name is asked).
- The report (the help text until the first run).

##### Step-by-step procedure

1. Run the Gaussian job with `IOp(9/40=4)` and a `%Chk`; convert the checkpoint
   with `formchk`.
2. In the NTO tab, browse to the `.fchk` (the output beside it is filled in) and
   click **Load**. The status line reports the number of states, whether the
   reference is restricted or unrestricted, and the number of MOs.
3. Choose a state and click **Run**. The report `name_state<n>_NTO.txt` is
   written to the output folder and shown at the bottom of the tab.
4. Check the **normalisation** in the report: the deviation from 0.5 (1 for an
   unrestricted reference) should be ≈ 0.00001. A larger value means missing
   coefficients (rerun with `IOp(9/40=4)`).
5. Read the **weights**. One pair near 1 means the excitation is that hole →
   particle transition. Several pairs of similar weight mean it is not a
   single-pair excitation (the report notes it when no pair exceeds 0.5).
6. Choose a pair and click **Display**. Compare hole and particle: where the
   electron density is taken from and where it goes (n → π\*, π → π\*, charge
   transfer between fragments, ...).
7. Export what you need. The cube files open in GaussView, VMD, Multiwfn or the
   Cube tab. The NTO `.fchk` opens in GaussView or Multiwfn, where the orbital
   "energies" are the weights. The dominant hole sits at the HOMO's place and
   the dominant particle at the LUMO's: for uracil, MOs 29 and 30.

##### Output files

| File | Content |
|---|---|
| `name_state<n>_NTO.txt` | header, Multiwfn's normalisation printout, the listed pairs with the hole/particle MO numbers they take in the NTO `.fchk`, the composition of each NTO in the canonical MOs (contributions ≥ 5 %), notes |
| `name_state<n>_NTO<k>_hole.cube`, `..._particle.cube` | orbital amplitude (a.u.) on the display grid, bohr; open shell: `name_state<n>_NTO<k>_alpha_hole.cube` etc. |
| `name_state<n>_NTO.fchk` | the input `.fchk` with the MO coefficients replaced by the NTOs and the orbital energies by their weights; everything else (density, basis) unchanged |
| `name_UVVis.json` | **Plot Spectrum**: every excited state (number, label, energy in eV, wavelength in nm, oscillator strength), the broadening settings and the sampled ε curve, as a plot manifest that **Plotter** (*JSON plot type*) reopens |

##### UV-Vis spectrum (Plot Spectrum)

**Plot Spectrum** draws *all* the states of the excitation file, whichever
state is selected, in the layout of GaussView's *Electronic Spectra* window:

- wavelength λ (nm) on the horizontal axis;
- one vertical stick per state, of height f (the oscillator strength; right
  axis), labelled with the state's number and label (e.g. `6: Singlet-B2`);
  dark states (f = 0, such as triplets) are drawn at zero height with their
  labels;
- the broadened absorption spectrum ε (L mol⁻¹ cm⁻¹, left axis), the sum of one
  Gaussian band per state (see
  [Formulas and Methods (NTO tab)](#formulas-and-methods-nto-tab), §7).

*Peak half-width (eV)* is the half width at half height of each band (blank =
0.333 eV, GaussView's default). The spectrum opens in the interactive plot
viewer (zoom, pan, save as image) and is saved at once as `name_UVVis.json`
in the output folder, so it can be reopened with **Plotter** and its numbers
reused elsewhere. It does not change the NTO results or the 3D view. The
oscillator strengths come from the Gaussian output; a plain-text excitation
file has none, so the button refuses it.

##### Interpretation and cautions

- **Weights.** For CIS and TDA they lie between 0 and 1 and sum to 1. For full TD
  (TDHF/TD-DFT) the de-excitation part is left out, as in Multiwfn and
  Gaussian, so a weight or the sum can slightly exceed 1. Much more than 1
  means the TD result itself is doubtful; TDA is then the better choice.
- **Open shell.** Alpha and beta NTOs are separate; both are listed, sorted by
  weight, and the largest tells which spin carries the excitation (for the
  methyl radical's first state, a beta pair with 0.98).
- **Signs and degeneracy.** The overall sign of a pair is arbitrary (here: the
  hole's largest MO coefficient is positive; hole and particle flip together).
  Pairs of equal weight (benzene's first state: two pairs of 0.4996) may come
  out as any rotation of one another. Only the space they span is defined, in
  this program as in every other.
- **When NTOs do not help.** If no pair dominates, the excitation genuinely
  involves several hole → particle pairs; discuss them together.

##### Validation

Five fixtures in `tests/data/nto/` are tested against Multiwfn 2026.7.15
(12 states): Multiwfn's own uracil example (PBE0/6-31G*, TD), formaldehyde TD
singlets and triplets and CIS/cc-pVDZ (pure d shells), the methyl radical
(UB3LYP, open shell) and benzene (degenerate pairs).

| Quantity | Agreement with Multiwfn |
|---|---|
| NTO weights | ≤ 3.3 × 10⁻⁹ (Multiwfn writes 8 significant digits) |
| hole and particle orbitals (sign- and rotation-free, 1500 points) | ≤ 6.4 × 10⁻⁹ of the orbital's maximum |
| uracil S3 (Multiwfn manual §4.18.6) | 0.865529, 0.134025, 0.000582, 0.000121, 0.000063; sum of the ten largest 1.000387 |

The exported NTO `.fchk` is read by Gaussian's cubegen (to cubegen's own
accuracy near the nuclei, ~5 × 10⁻⁵ of the maximum). The orbital values
themselves are checked by an exact identity: the squared occupied orbitals sum
to the electron density.

##### Limitations

- Gaussian files only (no ORCA, CP2K, GAMESS/Firefly, BDF).
- Single-reference methods only (CIS, TDA, TDHF, TD-DFT); not EOM-CC, CASSCF or
  other multi-reference excited states. No spin-flip TD and no restricted
  open-shell (ROHF/ROKS) references.
- The transition density is built from the excitation coefficients only (no
  X+Y option); NTOs between two excited states and the hole–electron analysis
  are not available.
- Accuracy is limited by the coefficients Gaussian prints (5 decimals, > 1e-4
  with `IOp(9/40=4)`): NTOs agree with Multiwfn, which reads the same numbers,
  to round-off, but with Gaussian's own `Pop=NTO` only to about 1e-4.
- Basis functions up to g (as the other `.fchk` tools); no periodic systems.

#### Analysis of Fukui Functions

**Tools > Analysis of Fukui Functions** shows where a molecule is most susceptible to nucleophilic (f⁺), electrophilic (f⁻) and radical (f⁰) attack, and the dual descriptor Δf = f⁺ − f⁻ (positive: nucleophilic attack; negative: electrophilic attack). The results indicate *likely* reactive regions. They are not predictions of reaction outcome: sterics, solvent and barriers are not included.

**Before you start (Gaussian).**
1. Optimize the molecule in the state in which it reacts (N electrons) and check that it is a minimum.
2. At that same geometry, run three single points with the same method and basis set: N, N+1 and N−1. Do not re-optimize the ions.

   | Molecule | N | N+1 | N−1 |
   |---|---|---|---|
   | neutral, closed shell | 0 1 | −1 2 | 1 2 |
   | neutral radical | 0 2 | −1 1 | 1 1 |

3. Run `formchk` on each `.chk`.

Steps 1–3 can be prepared by [Inputs > Gaussian Fukui Input Builder](#gaussian-fukui-input-builder): one Gaussian 16 linked job (optimization of N, then `Geom=Check` single points on N+1 and N−1) plus a formchk script.

**Setup tab.**
- The three `.fchk` files and an optional output folder (blank = next to the N file). Choosing a `<name>_N.fchk` with **Browse** fills in `<name>_Np1.fchk` and `<name>_Nm1.fchk` when those files exist and the fields are empty (the names the input builder writes).
- Fukui isovalue: blank = 0.005 a.u.
- ρ(N) mapping surface: blank = 0.01 a.u.
- Grid spacing: blank = 0.25 bohr.
- Atom-condensed indices: on by default. They are computed inside the tool from Hirshfeld charges of the three states (built-in free-atom densities, elements H–Ar); no external program is needed. Switch them off to save a few seconds. For other elements, or files that use effective core potentials, the run still finishes and a NOTE explains why.
- Write cube files: on by default.
- Cube grid spacing: blank = 0.1 bohr. The cubes are computed exactly at every point of their own grid, which covers the same box as the 0.25 bohr analysis grid but has about 15 times as many points, so GaussView and the 3D viewer's Cube tab draw smooth isosurfaces. Enter 0.25 to write the analysis grid itself (files about 15 times smaller). Acetone with 6-31++G(d,p) at 0.1 bohr: 161 × 156 × 154 points, about 52 MB per cube, 6 s for all four with `numba` installed.

**Run** checks the three files (same atoms, geometry, basis set and method; N+1 has one electron more, N−1 one fewer), then computes everything. The Results tab shows the report, which is also saved as `<name>_Fukui.txt`. It contains:
- the global indices (IP, EA, electronegativity, hardness = I − A, softness, electrophilicity, nucleophilicity);
- grid integrals: about 1 for f⁺/f⁻/f⁰ and 0 for Δf, a quick check that the files are not swapped;
- the condensed table q(N), q(N+1), q(N−1), f⁻, f⁺, f⁰, CDD, s⁻, s⁺, s⁰;
- a site ranking for each kind of attack;
- the surface extrema of each function on ρ(N) = 0.01 a.u., with the nearest atom;
- NOTES, e.g. degenerate frontier orbitals, diffuse basis functions, negative condensed values, unusual multiplicities.

Cubes: `<name>_fplus.cub`, `_fminus.cub`, `_fzero.cub`, `_dual.cub` (a.u. on a bohr grid; open them in GaussView, VMD, Multiwfn or the 3D viewer's Cube tab). The report header records the cube grid (`# cube_spacing_bohr`).

**Show 3D.** Pick a function and a display:
- *Isosurfaces (+/- isovalue)*: blue positive, red negative.
- *Positive isosurface only* (+isovalue, blue) or *Negative isosurface only* (−isovalue, red): one sign on its own, e.g. the nucleophilic (Δf > 0) or electrophilic (Δf < 0) regions of the dual descriptor. If the function never reaches that value, the view is cleared and the status line says so. Exported images get `_positive`/`_negative` in their suggested name.
- *Mapped on ρ(N) surface*: colour = value, with extrema markers (cyan minimum, orange maximum).

**Atoms** (Results tab): *Line*, *Ball and stick (CPK)* or *Space-filling (vdW)*, and *Show atom numbers* (1-based, the numbering used in the report). Both apply at once, with or without results. Space-filling spheres are larger than the Fukui lobes and hide them: use that style to judge how exposed a site is, then switch back.

Canvas background and **Export image...** work as in the ESP tool.

**Reading the results.**
- Compare only atoms that can actually react.
- For electrophilic attack, prefer large f⁻ and the most negative Δf. For nucleophilic attack, prefer the most positive Δf: f⁺ alone is the least reliable of the descriptors.
- If the descriptors disagree, compute the transition states.
- For degenerate HOMO/LUMO (benzene, C₆₀), use N+p/N−q states or Multiwfn's orbital-weighted Fukui functions (see the NOTE).

Theory, algorithms, validation and limitations: [Formulas and Methods (Fukui Function tool)](#formulas-and-methods-fukui-function-tool), and `docs/superpowers/specs/2026-10-16-fukui-function-design.md` for the verified Multiwfn workflow.

#### Energy Plots

Use **Tools > Energy plots** to visualize CPMD or gqteaMD energy files.

For CPMD ENERGY files, plot options include:

- Fictitious and ionic kinetic energy.
- Temperature.
- Kohn-Sham potential energy.
- Kohn-Sham plus ionic kinetic energy.
- Total energy.
- CPU time by step.

X-axis units:

- Steps.
- Femtoseconds.
- Picoseconds.

For gqteaMD energy files, select the appropriate file type and x-axis option in the tool.

##### JSON plot file

The **Plot type** dropdown offers a third option, **JSON plot file**, which redraws figures a gQTEA analysis tool has already produced — straight from the file, without recomputing anything.

Every analysis tool that shows interactive figures writes such a file next to your data. Select **JSON plot file**, press **Browse**, pick the file, and the message panel lists what it contains:

```text
Loaded 3 figure(s) from the JSON file:

  1. Continuous survival S(t), origin-averaged - 1512 points  [time (fs) vs S(t)]
  2. Intermittent correlation C(t) and rescaled R(t) - 2 curves, 1512 points each  [time (fs) vs correlation]
  3. Residence event duration distribution - 44 points  [event duration (fs) vs count]
```

Press **Plot** to open them in the interactive viewer (zoom, pan, save).

Because the file describes the figures completely, **all the other plot options are disabled** while this type is selected — the six CPMD plot switches, the simulation time step, the x-axis unit, and the gqteaMD column selectors. Switching back to either energy-file type restores them.

**Expected format.** A JSON array of figures, where each figure is either a single curve

```text
{"x": [...], "y": [...], "xlabel": "...", "ylabel": "...",
 "title": "...", "xlim": [lo, hi], "ylim": [lo, hi]}
```

or several labelled curves drawn together

```text
{"series": [{"x": [...], "y": [...], "label": "..."}, ...],
 "xlabel": "...", "ylabel": "...", "title": "..."}
```

`xlabel`, `ylabel`, `title`, `xlim` and `ylim` are optional. Within one curve, `x` and `y` must be lists of finite numbers of equal length.

Three tools write figures of a special kind, marked by a `"type"` key, and the Plotter redraws each exactly as the tool showed it:

| `"type"` | Written by | Drawn as |
| --- | --- | --- |
| `"bars"` | [Comparison of Molecular Geometric Parameters](#comparison-of-molecular-geometric-parameters) | grouped bars with error bars and an optional percentage line |
| `"plane_map"` | the QTAIM tab's *Plot 2D map* (`name_QTAIM_map.json`) | density contours, gradient and interbasin lines, critical points |
| `"uv_vis"` | the NTO tab's *Plot Spectrum* (`name_UVVis.json`) | the GaussView-style UV-Vis spectrum: broadened ε curve and labelled oscillator-strength sticks |

A UV-Vis figure has this form (every key is written by the NTO tab; the Plotter needs `states` and `curve`):

```text
{"type": "uv_vis", "title": "UV-Vis Spectrum",
 "states": [{"number": 1, "label": "Singlet-A\"", "energy_ev": 4.7968,
             "wavelength_nm": 258.47, "oscillator": 0.0001}, ...],
 "curve": {"wavelength_nm": [...], "epsilon": [...]},
 "broadening": {"shape": "gaussian", "hwhm_ev": 0.333, "sigma_ev": 0.39997, ...},
 "xlim": [lo, hi], "source": {"fchk": "...", "excitation_file": "..."}}
```

The message panel describes it as, for example:

```text
  1. UV-Vis Spectrum - 3 excited states  [Wavelength (nm) vs ε (L mol⁻¹ cm⁻¹)]
```

**Error handling.** A file that cannot be used is reported with a message identifying the figure at fault, so you can fix it without reading the JSON by hand:

| Problem | Message |
| --- | --- |
| File does not exist | `JSON file not found: <path>` |
| Not JSON, or empty | `The file is not valid JSON (line 1, column 1): Expecting value.` |
| Top level is an object | `The JSON file must contain a list of figures (a JSON array)...` |
| Empty array | `The JSON file contains no figures to plot.` |
| `x` and `y` differ in length | `Figure 1: 'x' and 'y' must have the same length (3 vs 2).` |
| Neither `x`/`y` nor `series` | `Figure 1 must provide either 'x' and 'y', or 'series'.` |
| A value is not a number | `Figure 1: 'y' must contain only numbers.` |
| UV-Vis figure without states | `Figure 1: 'states' must be a non-empty list of excited states.` |
| UV-Vis state without wavelength or strength | `Figure 1, state 2 needs 'wavelength_nm' and 'oscillator'.` |
| UV-Vis curve arrays of different length | `Figure 1: the curve's 'wavelength_nm' and 'epsilon' must have the same length (2000 vs 1999).` |

##### Multi-column data file

The fourth **Plot type**, **Multi-column data file**, plots any numeric table, whatever its number of columns. **Column 1 is the x axis; every other column becomes its own curve**, all drawn against the same x values in one figure in the interactive viewer. As in JSON mode, every other plot option is disabled while this type is selected.

```text
# RDF of O around Na+            <- comment, ignored
r      g(r)     integral         <- optional column-name header
0.05   0.000    0.000
0.10   0.012    0.001   # inline comments are ignored too
```

- Columns can be separated by spaces, tabs or commas.
- Blank lines and lines starting with `#`, `!`, `%`, `@` or `//` are comments and are skipped wherever they appear (so GROMACS `.xvg` files load as is).
- If the table has a **column-name header** (the last text line before the first row of numbers), the y-column names label the curves in the legend and the first name labels the x axis. Without a header the curves are called `Column 2`, `Column 3`, …. Names written inside a comment line (`# r g(r) integral`) are ignored like any other comment.
- The header is never read as data. Once the numbers start, every line must hold the same number of finite numbers.

**Customising the figure.** Two fields and a per-column list below the file row apply only to this plot type:

- **Table x-axis label** and **Table y-axis label** set the axis titles. Left blank, the x axis uses the first column's name (or `Column 1`) and the y axis reads `Value`.
- **Y columns to plot and their line styles** shows one row per y column once a file is loaded. The **switch** (labelled with the column name) chooses whether that column is plotted: all are on by default, and turning one off leaves it out of the figure and greys out its style. At least one must stay on. The **dropdown** sets that curve's style: **Solid** (default), **Dashed**, **Dotted**, **Dash-dot**, **Line + markers** or **Markers only**, independently of the others. The list scrolls when the table has many columns.

Labels and styles are reset when you load another file or change the plot type.

After **Browse**, the message panel lists the row count and the column names it found. Problems are reported with the line at fault, for example `Line 7 has 2 values, but the table has 3 columns.`, `Line 12 is not numeric: ...`, or `The column-name header (line 2) has 2 names, but the data has 3 columns.`

##### Bar chart data file

The fifth **Plot type**, **Bar chart data file**, draws a table with any number of columns as a **grouped bar chart**. **Column 1 gives one bar group per row** and may hold text labels (atom symbols such as `C1`, `O7`) or numbers; either way the values are shown on the x axis exactly as written. **Every other column is one bar series.** With several columns selected, each row's bars are drawn side by side and a legend names each series. As with the other table types, every other plot option is disabled while this type is selected.

A typical input is the condensed-index table written by **Analysis of Fukui Functions**:

```text
CONDENSED INDICES  (source: Hirshfeld ...)        <- title line, ignored
  atom      q(N)   q(N+1)   q(N-1)       f-       f+  ...   <- column-name header
   C1     -0.0390  -0.1675   0.0268   0.0658   0.1285  ...
   O7     -0.1978  -0.2444  -0.0518   0.1460   0.0465  ...
```

The layout rules are those of the multi-column data file (spaces, tabs or commas; comment and blank lines skipped; the header is the last text line before the first data row). The one difference is that a data row is a label **followed by numbers**: once the data starts, columns 2 onward must hold finite numbers, and a row with a different number of columns is reported with its line number (`Line 9: columns 2 onward must be numbers: ...`). A title line before the header must contain some text after its first word. A line such as `Step 5` reads as a label followed by a number, that is, as a data row, so put a `#` in front of it.

**Customising the figure.** The same **Table x-axis label** and **Table y-axis label** fields apply (blank: the first column's name and `Value`). Below them, **Columns to plot as bar series** shows one switch per data column once a file is loaded. All are on by default, at least one must stay on, and the series appear in the file's column order. Bars have no line style, so these rows have no style dropdown. For the Fukui example, switching on only `f-`, `f+` and `f0` compares the three Fukui functions atom by atom, and `CDD` on its own shows the dual descriptor's sign on each atom.

The figure opens in the interactive viewer (zoom, pan, save). As for the multi-column type, no file is written next to your data.

#### Molecular Axis Alignment

Use **Tools > Molecular Axis Alignment** to align an XYZ molecule along a chosen axis and save a new XYZ file.

Actions:

- **Browse**: select input XYZ.
- **Save As**: choose output path.
- Run the alignment action from the tool window.

#### Select Frames

Use **Tools > Select Frames** to extract a frame range from `TRAJEC.xyz`.

Inputs:

- Starting frame.
- Number of frames to skip between collected frames.
- Stop frame.
- Optional Gaussian settings.
- Charge and multiplicity when Gaussian input generation is enabled.

Outputs:

- `selected_frames.xyz`
- Optional Gaussian `.gjf` files.

#### Frame Selection by Interatomic Distance Range

Use **Tools > Frame selection by interatomic distance range** to extract frames where the distance between two atoms falls in a selected interval.

Inputs:

- Minimum distance.
- Maximum distance.
- Two atom labels.
- `TRAJEC.xyz`.

Typical outputs:

- A selected-frames XYZ file in the input directory.
- A TXT file with selected-frame distances and statistics.
- `frame_closest_to_average.xyz`.

This tool is useful after bond-length free-energy analysis, for example to select frames near the minimum Helmholtz free energy.

#### CPMD Input to XYZ Converter

Use **Tools > CPMD Input to XYZ Converter** to extract coordinates from a CPMD input file and write an XYZ file.

Actions:

- Browse for a CPMD input file.
- Enter an output filename.
- Click **Convert**.

The converter reads `CELL` values from the `&SYSTEM` block and atom coordinates from the `&ATOMS` block.

#### SH Geometry Analyzer

Use **Tools > SH Geometry Analyzer** to group surface hopping trajectory frames by electronic state.

Inputs:

- Number of states.
- Root directory containing simulation subfolders.

Each simulation subfolder should contain:

- `SH_STATE.dat`
- `TRAJEC.xyz`

Outputs:

- `stateX.xyz` files inside each subfolder.
- Consolidated `stateX.xyz` files in the root directory.
- `sh_avg_perc.dat` with average state occupancy percentages.

If a subfolder is missing `SH_STATE.dat` or `TRAJEC.xyz`, the tool reports a warning and skips that subfolder.

#### Convert cp.x .pos File to trajec.xyz

Use **Tools > Convert *.pos file to trajec.xyz** to convert a Quantum ESPRESSO `cp.x` position trajectory to XYZ.

Inputs:

- `cp.x` input file, used to read atom labels and `nat`.
- `*.pos` trajectory file.

Output:

- XYZ trajectory file with coordinates converted from Bohr to Angstrom.

#### Compute Forces from cp.x .for File

Use **Tools > Compute forces from cp.x *.for file** to convert a Quantum ESPRESSO force trajectory into a trajectory force file.

Inputs:

- `cp.x` input file.
- `*.for` force trajectory file.

Click **Convert** to generate the converted force output.

#### Coordinate Converter

Use **Tools > Coordinate converter** to convert between molecular coordinate file formats and, in the same step, between Cartesian coordinates in **Angstrom (Å)** and **bohr** (atomic units of length), in either direction.

Controls:

- Input format selection, with its **Unit** (Angstrom or Bohr).
- Output format selection, with its **Unit** (Angstrom or Bohr).
- Browse for input file.
- Choose output file.
- **Convert**.

Supported formats: XYZ, PDB, ENT, MOL, MOL2, Gaussian (`.gjf`/`.com`/`.gau`) and HIN. The converter writes a new file rather than modifying the original, and refuses an output path equal to the input path.

Angstrom ↔ bohr conversion:

- The factor is the CODATA 2022 Bohr radius, **1 bohr = 0.529177210544 Å** (1 Å ≈ 1.8897261259 bohr). Each coordinate is divided (Å → bohr) or multiplied (bohr → Å) by this value in double precision.
- The output format can equal the input format, so **XYZ (Angstrom) → XYZ (Bohr)** is a pure unit conversion (this replaces the former *Convert coordinates from Å to Bohr* tool). Choosing the same format *and* the same unit is rejected, since there would be nothing to convert.
- XYZ files do not record their unit, so set the input **Unit** to match your file. Converted XYZ files state it on the comment line, e.g. `Converted file (coordinates in bohr)`.
- Gaussian inputs declare their unit in the route section. When the output unit is Bohr the converter writes `Units=Bohr` on the route line, and when a Gaussian input containing `Units=AU` or `Units=Bohr` is loaded, the input unit is set to Bohr automatically.
- PDB, ENT, MOL, MOL2 and HIN are defined in Angstrom by their specifications, so their **Unit** selector is locked to Angstrom. To get one of these from a bohr file, select XYZ (or Gaussian) with Unit = Bohr as input and the Angstrom format as output.
- Precision: XYZ, Gaussian and HIN coordinates are written with 10 decimal places, so a file with up to 8 decimals converts Å → bohr → Å (or bohr → Å → bohr) back to its original digits. PDB keeps 3 decimals and MOL/MOL2 keep 4, as their fixed-column formats require.
- When the units change, the default output name gains the unit, e.g. `2criegee.xyz` → `2criegee_converted_bohr.xyz`.

Example (`2criegee.xyz`, first atom): `C 1.99447596 0.39215795 -0.00041257` (Å) becomes `C 3.7690133291 0.7410711236 -0.0007796443` (bohr).

#### Unit Converter

Use **Tools > Unit Converter** for scalar unit conversion.

Controls:

- Category.
- Value.
- From unit.
- To unit.
- **Convert**.

If the value is blank or a unit is unsupported, the converter shows an error explaining what to fix.

## Configuration

gqteaWinToga does not currently use a persistent user configuration file. Settings are entered directly in each tool window for each run.

Customization points:

- Choose output directories using **Browse**, **Save Dir**, **Save As**, or Save dialogs when available.
- Leave output-name fields blank in some builders to use automatic names such as `<xyzname>_cp.in`, `<xyzname>_orca.inp`, or `<xyzname>_gqteaMD.toml`.
- Toggle optional outputs such as CSV files, plot display, Gaussian input generation, random solvent rotation, density calculation, and minimum-image clash detection.
- Edit generated input files before running external simulation engines, especially pseudopotential names, basis settings, memory, processor counts, and advanced blocks.

Generated files are usually written in the selected output directory or beside the selected input file. Check the status text in each tool window after running a calculation.

## Practical Examples

### Example 1: Find a Bond-Length Free-Energy Minimum

1. Open **Geometry > Bond length analysis**.
2. Load `TRAJEC.xyz`.
3. Enter:

```text
Maximum r: 3.0
Simulation time step: 5
Sampling interval: 1
Temperature: 300
Atom labels: 1 2
Histogram bin width: 0.02
```

4. Click **Analyze**.
5. Read the minimum straight from the **RESULTS** block in the output box (`Lowest free energy: ... at ... Angstrom`), or open the generated `free_energy_*.dat` or its plot.
6. Use **Tools > Frame selection by interatomic distance range** to extract frames near that distance.

### Example 2: Prepare a CPMD Collision Restart

1. Run CPMD for one step and obtain `GEOMETRY.xyz`.
2. Open **Inputs > Collision Input**.
3. Select `cpmd` as the MD engine.
4. Enter attacker atom labels, initial velocity, and target coordinates.
5. Browse to `GEOMETRY.xyz`.
6. Click **Input Builder**.
7. Rename the generated collision geometry as needed for your CPMD restart.
8. Restart CPMD with the `GEOFILE` keyword so the new velocities are used.

### Example 3: Group Surface Hopping Frames by State

1. Organize simulation folders under one root directory.
2. Ensure each subfolder contains `SH_STATE.dat` and `TRAJEC.xyz`.
3. Open **Tools > SH Geometry Analyzer**.
4. Enter the number of states.
5. Browse to the root directory.
6. Click **Extract Frames**.
7. Review `stateX.xyz` files and `sh_avg_perc.dat`.

### Example 4: Convert a cp.x Trajectory

1. Open **Tools > Convert *.pos file to trajec.xyz**.
2. Browse to the `cp.x` input file.
3. Browse to the `*.pos` file.
4. Click **Convert**.
5. Use the generated XYZ trajectory in the viewer or analysis tools.

### Example 5: Analyze All Solute Bond Distances

Suppose atoms `1 2 3 4` are the solute and the remaining atoms are solvent. A small trajectory might begin like this:

```text
5
frame 1
C  0.0  1.0  0.0
C  0.0  0.0  0.0
C  1.0  0.0  0.0
H  1.0  0.0  1.0
O  0.0  0.5  0.0
5
frame 2
C  0.0  1.1  0.0
C  0.0  0.0  0.0
C  1.1  0.0  0.0
H  1.0  0.0  1.1
O  0.0  0.5  0.0
```

Atom `5` is a solvent atom close to the solute, but it must not contribute to the solute bond analysis.

1. Open **Geometry > All bond distance analysis**.
2. Browse to the trajectory file.
3. Set **Maximum connection distance (A)** to `1.5` (or leave blank for the `1.7` default).
4. Set **Solute atom indices** to (range syntax also works, e.g. `1-4`):

```text
1 2 3 4
```

5. Keep the output filename as:

```text
all_bond_analysis_combined.txt
```

6. Click **Analyze**.

The single combined file carries the run metadata and one row per pair (identity + first-bonded distance + statistics + occurrence):

```text
# gQTEA All Bond Distance Analysis
# frames_used 2
# max_connection_distance 1.5
# atom_scope solute_atoms
# solute_atom_indices 1 2 3 4
# periodic_boundary none
# row atom_i atom_j element_i element_j first_bonded_distance average variance standard_deviation occurrence_fraction frames_bonded
     1        1        2        C        C       1.00000000       1.05000000       0.00250000       0.05000000       1.00000000          2
```

Only solute-solute pairs are reported. To apply periodic boundaries, also set **Cell lattices a b c**; the header then records `# cell_lengths a b c` instead of `# periodic_boundary none`.

### Example 6: Analyze All Solute Bond Angles

Use the same idea for all connected solute angles. For a water-like solute with atoms `1 2 3` and a nearby solvent atom `4`:

```text
4
frame 1
H  0.0  0.0  0.0
O  1.0  0.0  0.0
H  1.0  1.0  0.0
H  1.0  0.5  0.0
4
frame 2
H  0.0  0.0  0.0
O  1.0  0.0  0.0
H  0.5  0.8660254  0.0
H  1.0  0.5  0.0
```

Atom `4` is solvent and should be excluded.

1. Open **Geometry > All bond angle analysis**.
2. Browse to the trajectory file.
3. Set **Maximum connection distance (A)** to `1.1`.
4. Set **Solute atom indices** to:

```text
1 2 3
```

5. Keep the output filename as:

```text
all_angle_analysis_combined.txt
```

6. Click **Analyze**.

The combined file carries the run metadata and one row per triplet (identity + first-present angle + statistics + occurrence):

```text
# gQTEA All Bond Angle Analysis
# frames_used 2
# max_connection_distance 1.1
# atom_scope solute_atoms
# solute_atom_indices 1 2 3
# periodic_boundary none
# row atom_i atom_j atom_k element_i element_j element_k first_present_angle_degrees average_angle_degrees variance_degrees2 standard_deviation_degrees occurrence_fraction frames_present
     1        1        2        3        H        O        H      90.00000000      74.99999995     225.00000163      15.00000005       1.00000000          2
```

The angle is interpreted as `atom_i-atom_j-atom_k`, so atom `2` is the central atom in this example. To apply periodic boundaries, also set **Cell lattices a b c**; the header then records `# cell_lengths a b c`.

### Example 7: Analyze All Solute Dihedral Angles

For a solute chain `1-2-3-4` in a solvent environment, use **All dihedral angle analysis**. Example trajectory:

```text
5
frame 1
C  0.0  1.0  0.0
C  0.0  0.0  0.0
C  1.0  0.0  0.0
H  1.0  0.0  1.0
O  0.0  0.5  0.0
5
frame 2
C  0.0  1.0  0.0
C  0.0  0.0  0.0
C  1.0  0.0  0.0
H  1.0  0.0 -1.0
O  0.0  0.5  0.0
```

Atom `5` is solvent and should not be used to build dihedral paths.

1. Open **Geometry > All dihedral angle analysis**.
2. Browse to the trajectory file.
3. Set **Maximum connection distance (A)** to `1.1`.
4. Set **Solute atom indices** to:

```text
1 2 3 4
```

5. Keep the output filename as:

```text
all_dihedral_analysis_combined.txt
```

6. Click **Analyze**.

The combined file carries the run metadata and one row per quadruplet (identity + first-present dihedral + statistics + occurrence):

```text
# gQTEA All Dihedral Angle Analysis
# frames_used 2
# max_connection_distance 1.1
# dihedral_convention signed_degrees_minus180_to_180
# atom_scope solute_atoms
# solute_atom_indices 1 2 3 4
# periodic_boundary none
# row atom_i atom_j atom_k atom_l element_i element_j element_k element_l first_present_dihedral_degrees average_dihedral_degrees variance_degrees2 standard_deviation_degrees occurrence_fraction frames_present
     1        1        2        3        4        C        C        C        H     -90.00000000       0.00000000    8100.00000000      90.00000000       1.00000000          2
```

Dihedrals are reported as signed degrees using the range `(-180, 180]`. To apply periodic boundaries, also set **Cell lattices a b c**; the header then records `# cell_lengths a b c`.

### Example 8: Read a Fukui Function Figure (phenol, f⁺)

This example explains one figure from **Tools > Analysis of Fukui Functions**: the f⁺ isosurfaces of phenol (C₆H₅OH) at B3LYP/6-31G(d). The three `.fchk` files came from one [Gaussian Fukui Input Builder](#gaussian-fukui-input-builder) job:
- the neutral molecule, optimized: charge 0, multiplicity 1, 50 electrons;
- the anion: −1 2, 51 electrons;
- the cation: +1 2, 49 electrons.

Because the builder names the files `<name>_N.fchk`, `<name>_Np1.fchk` and `<name>_Nm1.fchk`, choosing the N file filled in the other two automatically.

Settings on the *Results* tab:
- **Function**: `f+`;
- **Display**: *Isosurfaces (+/- isovalue)*;
- **Atoms**: *Line*;
- **Canvas background**: *Black*;
- Fukui isovalue: 0.005 a.u. (the default).

#### What the figure shows

f⁺(r) = ρ(N+1) − ρ(N) is the electron density of the anion minus that of the neutral molecule, both at the same geometry. It shows where an **added electron goes**, so its large positive regions mark the sites most susceptible to **nucleophilic attack**.

| Element of the figure | Meaning |
|---|---|
| **Blue** lobes | f⁺ ≥ +0.005 a.u.: the density *increases* when the electron is added. |
| **Red** lobes | f⁺ ≤ −0.005 a.u.: the density *decreases* (the other electrons relax away). |
| Colour bar `−5.000e-3 … 5.000e-3` | Only the ±isovalue thresholds. In this mode each lobe has a single colour, so the intermediate colours of the bar are not used. |
| Axis triad | x to the right, y up, z out of the screen. Phenol lies in the screen (xy) plane, so the ring is seen face-on. |
| Thin lines | The bonds (*Line* atom style). Here the OH group is at the upper right. |

#### How to read the lobes

1. **Four ring carbons carry large blue lobes: the two *ortho* and the two *meta* carbons.**
   - The carbon bearing OH (*ipso*) and the one opposite it (*para*) have almost no f⁺; only a small red dot sits near the para carbon.
   - This is the shape of phenol's lowest unoccupied orbital (LUMO). It is a π* orbital whose nodal plane passes through the ipso and para carbons, and the added electron mostly enters it.
2. **The lobes are π lobes, above and below the ring plane.**
   - Seen face-on, the lobe above the plane and the one below overlap on top of each carbon, so each looks like a single blob.
   - Rotate the view to look along the ring plane (edge-on) to see both halves.
3. **The red regions are relaxation, not a second site.**
   - The added π electron repels the other electrons, so the density falls slightly beside the ring carbons and around the oxygen.
   - These changes are smaller than the blue lobes.
4. **The small blue spheres on the hydrogens** are a minor spill-over of the added density onto the H atoms.

The condensed (per-atom) indices in the report put numbers on this picture; scroll the *Results* tab to the condensed table. As a reference, the verification run of the same calculation in Gaussian 16 gave:

| | ipso C (bears OH) | ortho C | meta C | para C | O |
|---|---|---|---|---|---|
| f⁺ | lower | 0.128, 0.123 | 0.128, 0.131 | lower | lower |
| f⁻ | 0.100 | 0.092 (one ortho C) | lower | 0.144 | **0.146** |

"lower" means the atom is not among the four largest values of that function. That run used the atom order of `tests/data/fukui/inputs/phenol_N.gjf`: C1 para, C2/C6 meta, C3/C5 ortho, C4 ipso, O12. Check the atom numbers of your own file (*Show atom numbers*) before comparing.

#### The global indices in the report

| Quantity | Value | Meaning |
|---|---|---|
| E(N), E(N+1), E(N−1) | −307.464861, −307.383596, −307.163412 Hartree | total energies of the three states |
| Vertical IP = E(N−1) − E(N) | 0.301449 Ha = 8.20 eV | energy to remove an electron |
| Vertical EA = E(N) − E(N+1) | **−0.081265 Ha = −2.21 eV** | **negative**: the anion lies *above* the neutral molecule |
| Electronegativity χ = (I + A)/2; chemical potential μ = −χ | 2.996 eV; −2.996 eV | |
| Hardness η = I − A | 0.382714 Ha = 10.41 eV | Multiwfn convention, without the factor ½ |

#### Cautions

- **The negative electron affinity limits how literally f⁺ can be read.**
  - Phenol does not bind an extra electron in the gas phase. In this calculation the anion holds together only because 6-31G(d) has no diffuse functions that would let the electron leave.
  - f⁺ is still a fair picture of the π* orbital that accepts an electron, and that is how it is normally used, but treat it qualitatively.
  - Adding diffuse functions is not a simple fix: with B3LYP-type functionals the N+1 density then tends to become over-diffuse. The analysis tool writes a NOTE when the basis has diffuse functions; see *Formulas and Methods (Fukui Function tool)*, subsection 13.
- **For phenol, f⁻ and the dual descriptor are usually more informative than f⁺.**
  - Phenol is an electron-rich ring that reacts mainly with electrophiles. Set **Function** to `f-` to see where electrophiles attack: the condensed f⁻ above is largest on O and on the para carbon, then on the ipso and ortho carbons, in line with OH directing substitution to the ortho and para positions.
  - `DD` (the dual descriptor, f⁺ − f⁻) shows both in one figure: positive where nucleophiles attack, negative where electrophiles attack.
- **Show extrema markers** has no effect in the *Isosurfaces* display. The markers appear with *Mapped on ρ(N) surface*.
- The mapping surface used in this run was ρ(N) = 0.001 a.u. (shown as `surface` in the report header). The default and recommended value is 0.01 a.u., because the Fukui functions are very small on the 0.001 surface. Use 0.01 when you switch to the mapped display.

### Example 9: QTAIM Analysis of a Hydrogen Bond (water dimer)

This example follows the [step-by-step procedure](#qtaim-topology-analysis) for the water dimer at B3LYP/6-31G(d), the file `tests/data/qtaim/water_dimer_631gd.fchk` shipped with the program (input: `tests/data/qtaim/inputs/water_dimer_631gd.gjf`, the S22 geometry, O···O 2.90 Å). Atoms: O1, H2, H3 form the donor molecule (H3 points at the acceptor), O4, H5, H6 the acceptor.

#### 1. Run

QTAIM tab → Browse to the file → all four starting-point switches on, Advanced blank → **Run**. Summary:

```
(3,-3) 6   (3,-1) 5   (3,+1) 0   (3,+3) 0     Poincare-Hopf = 1 (satisfied)
```

Six atoms give six (3,−3) points; four O–H bonds plus one hydrogen bond give five (3,−1) points; 6 − 5 = 1, so the topology is complete.

#### 2. The bond critical points (from the report)

| Row | Bond | ρ | ∇²ρ | ε | G | V | H | \|V\|/G | H/ρ |
|---|---|---|---|---|---|---|---|---|---|
| 7 | O1–H2 (covalent) | 0.3576 | −1.8730 | 0.022 | 0.0699 | −0.6081 | −0.5382 | 8.70 | −1.505 |
| 8 | O1–H3 (donor O–H) | 0.3470 | −1.8813 | 0.022 | 0.0648 | −0.5999 | −0.5351 | 9.26 | −1.542 |
| 9 | **H3···O4 (hydrogen bond)** | **0.0280** | **+0.0799** | 0.037 | 0.0221 | −0.0242 | **−0.0021** | **1.10** | −0.076 |
| 10, 11 | O4–H5, O4–H6 | 0.3540 | −1.9001 | 0.022 | 0.0643 | −0.6037 | −0.5393 | 9.39 | −1.524 |

All values in atomic units.

#### 3. How to read them

1. **There is a bond path from H3 to O4.** The (3,−1) point 9 lies on it. The path joins the (3,−3) points of H3 and O4, which are 1.98 Å apart, and is practically straight (path length / distance = 1.0003).
2. **Covalent O–H bonds** (rows 7, 8, 10, 11): large ρ (about 0.35 a.u.), strongly negative ∇²ρ (charge concentrated between the nuclei), strongly negative H, and |V|/G ≈ 9 (far above 2). This is the *shared-shell* signature.
3. **The hydrogen bond** (row 9) is a *closed-shell* interaction: ρ is small (0.028 a.u.) and ∇²ρ is positive (charge depleted at the bond point). Both values fall inside the Koch–Popelier ranges for hydrogen bonds (ρ 0.002–0.034 a.u., ∇²ρ 0.024–0.139 a.u.).
4. **But it is not purely electrostatic**: H is slightly negative (−0.0021 a.u.) and |V|/G = 1.10 lies between 1 and 2, the *intermediate* region of Espinosa et al. — a hydrogen bond with a small covalent contribution, as expected for the water dimer.
5. **The donor O–H bond is weakened by the hydrogen bond**: ρ at O1–H3 (0.3470) is lower than at the free O1–H2 (0.3576) of the same molecule.
6. **An estimate of the hydrogen-bond energy** (not computed by the program; by hand): the relation of Emamian, Lu *et al.* (2019) for neutral hydrogen bonds, BE ≈ −223.08 ρ<sub>BCP</sub> + 0.7423 kcal/mol, gives −223.08 × 0.0280 + 0.74 = **−5.5 kcal/mol**, close to the accepted binding energy of the water dimer (about −5 kcal/mol). Espinosa's older BE ≈ V/2 = −0.0121 Hartree = −7.6 kcal/mol overestimates it, as Emamian *et al.* found in general.

#### 4. Figures

- **3D**: **Display**; then Canvas → Background *White*, Molecule style *Design style*, Element symbols and Atom numbers on → **Export 3D image...**. The (3,−1) points of the O–H bonds lie inside the H discs, so they are painted over the discs; atom numbers sit at the upper right of each atom, CP numbers at the lower right of each point.
- **2D**: Plane *Three atoms*, atoms `1, 2, 3` (the donor molecule's plane, which also contains O4) → **Plot 2D map**. The map shows the contours of ρ, the gradient lines leaving each nucleus, the bond paths O1–H2, O1–H3 and H3···O4, and the deep-blue interbasin paths: the one through point 9 separates the two molecules. **Export 2D image...** saves it.

#### 5. A ring and a cage (bicyclo[1.1.1]pentane)

The second fixture, `tests/data/qtaim/bicyclopentane_631gd.fchk` (B3LYP/6-31G(d), optimized in D₃h symmetry and checked by a frequency calculation: no imaginary frequency; C–C 1.556 Å), shows the other two point types:

```
(3,-3) 13   (3,-1) 14   (3,+1) 3   (3,+3) 1     Poincare-Hopf = 1 (satisfied)
```

- 13 atoms; 14 bonds (six C–C, eight C–H); three four-membered rings (three (3,+1) points); one cage (one (3,+3) point): 13 − 14 + 3 − 1 = 1. Ring paths join each (3,+1) to the (3,+3) point.
- There is **no** bond path between the two bridgehead carbons, although they are only 1.88 Å apart: the density between them has no saddle point.
- The strained C–C bonds are **curved** ("banana bonds"): their path is 1.0027 × the internuclear distance (a straight bond gives 1.0000), bowed outwards from the cage. ρ<sub>BCP</sub> = 0.233 a.u. is a single-bond value (benzene's C–C bonds: 0.311 a.u.).
- The three ring points (ρ = 0.0994 a.u.) lie only 0.25 Å from the cage point (0.0969 a.u.), with almost the same density: in this small cage the density is very flat, and the topology is close to a *catastrophe point*, where a small change of geometry can make ring and cage points merge and disappear. Such flat regions are where the search is most sensitive to its parameters.

### Example 10: Natural Transition Orbitals of Uracil S3

This is the example of the Multiwfn manual (§4.18.6). The files are in `tests/data/nto/`: `uracil.fch` and `uracil.out`, from Gaussian 09, `# PBE1PBE/6-31G* TD IOp(9/40=4)`, three singlet states.

#### 1. Why NTOs

Gaussian lists the third excited state (6.018 eV, 206 nm, f = 0) as a mixture of four main MO pairs:

```
 Excited State   3:      Singlet-A"     6.0180 eV  206.02 nm  f=0.0000
      26 -> 30         0.54135
      26 -> 31        -0.20634
      28 -> 30        -0.15424
      28 -> 31         0.36715
```

The largest pair carries only 2 × 0.541² = 58.5 %, so no single MO pair describes the excitation.

#### 2. Run

NTO tab: browse to `uracil.fch` (`uracil.out` is filled in), **Load**, choose `3: Singlet-A" 6.0180 eV 206.02 nm f=0.0000`, **Run**. The report `uracil_state3_NTO.txt` gives:

```
NORMALIZATION (as Multiwfn)
  Sum of square of excitation coefficients                        0.500203
  Negative of the sum of square of de-excitation coefficients    -0.000196
  Sum of above two values                                         0.500007
  Deviation to expected normalization value (0.5)                 0.000007
  Sum of all NTO weights                                          1.000406

  pair  spin    weight         %   cum. %   hole  particle
   1   alpha   0.865529   86.55   86.55     29     30
   2   alpha   0.134025   13.40   99.96     28     31

  pair  1 alpha hole:     MO 26 74.9%, MO 28 25.1%
  pair  1 alpha particle: MO 30 68.9%, MO 31 30.4%
```

The weights are those of the Multiwfn manual (0.865529, 0.134025, ...). The deviation of 0.000007 shows that `IOp(9/40=4)` captured essentially all of the coefficients.

#### 3. How to read it

- One pair carries 86.6 % of the excitation: the NTO transformation has turned four comparable MO pairs into one dominant hole → particle pair plus a 13.4 % companion.
- **Display** pair 1 with *Hole and particle*. The hole is an in-plane lone pair on the carbonyl oxygen. The particle is a π\* orbital over the ring and the C=O groups. S3 is therefore an n → π\* excitation, consistent with its A″ symmetry and zero oscillator strength.
- Each NTO is a mixture of canonical MOs (hole 1 = 75 % MO 26 + 25 % MO 28). That is why the canonical MOs alone gave no clear picture.
- **Export NTO .fchk** writes `uracil_state3_NTO.fchk`. In GaussView or Multiwfn, MO 29 is hole 1 and MO 30 is particle 1, with "energies" 0.865529.

#### 4. Cautions

- The weights sum to 1.0004, slightly above 1, because this is full TD-DFT and the de-excitations (Σ = 0.000196) are not part of the NTOs. The report's note says so.
- The overall sign of the hole/particle pair is arbitrary. Only the relative phases within each orbital carry meaning.

#### 5. The UV-Vis spectrum

Press **Plot Spectrum** (the *Peak half-width* field left blank, i.e. 0.333 eV). All three states are plotted, whichever one is selected, and `uracil_UVVis.json` is written to the output folder:

| State | Label | λ (nm) | f |
| --- | --- | --- | --- |
| 1 | Singlet-A″ | 258.47 | 0.0001 |
| 2 | Singlet-A′ | 228.67 | 0.1317 |
| 3 | Singlet-A″ | 206.02 | 0.0000 |

- The spectrum is one band, peaking at 228.7 nm with ε = 5333 L mol⁻¹ cm⁻¹: the bright S2 alone gives 1.3062974 × 10⁸ × 0.1317 / 3226 cm⁻¹ = 5333.
- S1 and S3, both A″, are essentially dark: their sticks are at, or almost at, zero height, so the n → π\* state analysed above (S3) does not appear in the absorption at all (run its NTOs to see what S1 is). A spectrum alone would never reveal it; the NTOs and the f values together tell which excitations are bright and what they are.
- Reopen the file later with **Plotter** (*JSON plot file*), or change the half-width and press **Plot Spectrum** again to compare band widths.

## Formulas and Methods (All Bond / Angle / Dihedral tools)

This section documents every equation implemented in the three connectivity-based
tools — [allBondAnalysis.py](allBondAnalysis.py), [allAnglesAnalysis.py](allAnglesAnalysis.py),
and [allDihedralAnalysis.py](allDihedralAnalysis.py). The first four equations are
shared by all three; the last two are specific to the angle and dihedral tools.
Notation: $\mathbf{r}_a$ is the Cartesian position of atom $a$ in a given frame,
$r_\text{cut}$ is the *maximum connection distance*, and $\mathbf{L}=(a,b,c)$ are
the orthorhombic cell lattices.

### 1. Interatomic vector and minimum-image convention (PBC)

The raw vector between two atoms and, when a cell is supplied, its
minimum-image replacement (applied component-wise):

$$
\Delta\mathbf{r}_{ab} = \mathbf{r}_a - \mathbf{r}_b,
\qquad
\Delta\mathbf{r}_{ab} \;\leftarrow\; \Delta\mathbf{r}_{ab}
- \mathbf{L}\,\operatorname{round}\!\left(\frac{\Delta\mathbf{r}_{ab}}{\mathbf{L}}\right)
$$

- **Represents:** the shortest vector between two atoms under periodic boundary
  conditions; it lets bonds that span a box face be measured correctly and
  tolerates unwrapped coordinates.
- **Where/how used:** in the bond tool it is applied inside `_pair_distances`; in
  the angle and dihedral tools the same operation is the helper `_minimum_image`,
  reused by `_bonded_adjacency`, `_calculate_angles`, and `_calculate_dihedrals`.
  When the **Cell lattices** field is left blank the subtraction is skipped
  (`cell_lengths is None`) and raw vectors are used. Validity requires
  $r_\text{cut} \le \tfrac{1}{2}\min(a,b,c)$, which the tools enforce in
  `_prepare_run`.
- **References:** M. P. Allen and D. J. Tildesley, *Computer Simulation of
  Liquids*, 2nd ed., Oxford University Press (2017), §1.6.3; D. Frenkel and
  B. Smit, *Understanding Molecular Simulation*, 2nd ed., Academic Press (2002).

### 2. Interatomic distance and the connection criterion

$$
d_{ab} = \lVert \Delta\mathbf{r}_{ab} \rVert
       = \sqrt{\Delta x_{ab}^{2} + \Delta y_{ab}^{2} + \Delta z_{ab}^{2}},
\qquad
\text{atoms } a,b \text{ are bonded} \iff d_{ab} \le r_\text{cut}.
$$

- **Represents:** the Euclidean interatomic distance (using the minimum-image
  vector above) and the per-frame test that decides whether a pair is "bonded".
- **Where/how used:** the distance is the quantity averaged by the **bond** tool;
  in all three tools the criterion $d_{ab}\le r_\text{cut}$ is re-evaluated **every
  frame** to build that frame's connectivity graph (`_pair_distances` /
  `_bonded_adjacency`). An i-j-k angle is present in a frame only when both
  $d_{ji}\le r_\text{cut}$ and $d_{jk}\le r_\text{cut}$; an i-j-k-l dihedral only
  when all three of $d_{ij},d_{jk},d_{kl}\le r_\text{cut}$.
- **References:** Allen and Tildesley (2017), Ch. 1; standard Euclidean norm.

### 3. Running mean, population variance, and standard deviation (Welford)

For a quantity $x$ (a distance, angle, or dihedral) accumulated over the
$N$ frames in which the item is present, with $x_n$ the value in the $n$-th such
frame:

$$
\mu_n = \mu_{n-1} + \frac{x_n - \mu_{n-1}}{n},
\qquad
M_{2,n} = M_{2,n-1} + \left(x_n - \mu_{n-1}\right)\left(x_n - \mu_n\right),
$$

$$
\sigma^{2} = \frac{M_{2,N}}{N},
\qquad
\sigma = \sqrt{\sigma^{2}}.
$$

- **Represents:** the average $\mu\equiv\bar{x}$, the **population** variance
  $\sigma^2$ (divided by $N$, not $N-1$), and the standard deviation, computed in a
  single numerically stable streaming pass.
- **Where/how used:** implemented as the masked Welford update `_accumulate` in the
  bond tool and inline within `_accumulate_frame` in the angle/dihedral tools;
  finalized (variance $=M_2/N$, $\sigma=\sqrt{\sigma^2}$) in `_finalize`. Each
  connected pair/triplet/quadruplet keeps its own $(N,\mu,M_2)$ accumulator that
  advances **only in the frames where the item is present**.
- **References:** B. P. Welford, "Note on a Method for Calculating Corrected Sums
  of Squares and Products", *Technometrics* **4**(3), 419–420 (1962); D. E. Knuth,
  *The Art of Computer Programming*, Vol. 2, 3rd ed., §4.2.2.

### 4. Occurrence (fraction of frames present)

$$
f = \frac{N_\text{present}}{N_\text{frames}}
$$

- **Represents:** how persistent a connection is — the fraction of the analysed
  trajectory in which the pair/triplet/quadruplet is within the cutoff. $f=1$ means
  present in every frame; a smaller value flags a bond that forms or breaks.
- **Where/how used:** reported as `occurrence_fraction` (with the raw count
  `frames_bonded`/`frames_present`) in each combined output file and in the on-screen
  summary; $N_\text{present}$ is the Welford count $N$ and $N_\text{frames}$ is
  `self.num_frames`.
- **References:** standard occupancy/occurrence definition used in hydrogen-bond and
  contact analysis, e.g. Allen and Tildesley (2017), Ch. 2 (time averages).

### 5. Bond angle (angle tool)

For the triplet i-j-k with central atom $j$, using the minimum-image bond vectors
$\mathbf{v}_{ji}=\mathbf{r}_i-\mathbf{r}_j$ and $\mathbf{v}_{jk}=\mathbf{r}_k-\mathbf{r}_j$:

$$
\theta_{ijk} = \frac{180}{\pi}\,
\arccos\!\left(
\operatorname{clip}\!\left(
\frac{\mathbf{v}_{ji}\cdot\mathbf{v}_{jk}}
     {\lVert\mathbf{v}_{ji}\rVert\,\lVert\mathbf{v}_{jk}\rVert},\,
-1,\,1\right)\right)
$$

- **Represents:** the bond angle in degrees, in the range $[0,180]$. The clip to
  $[-1,1]$ guards against tiny floating-point overshoot of the cosine before
  `arccos`.
- **Where/how used:** the vectorized helper `_calculate_angles` computes this for
  every i-j-k built from the current frame's connectivity; a degenerate case
  (two atoms coincident, i.e. a zero-length vector) yields `NaN` and that frame is
  skipped for the affected triplet rather than aborting the run.
- **References:** dot-product angle definition; Allen and Tildesley (2017), §1.3.
  Numerical clipping of the argument to `arccos` is standard practice.

### 6. Dihedral (torsion) angle (dihedral tool)

For the quadruplet i-j-k-l, using the minimum-image bond vectors
$\mathbf{b}_1=\mathbf{r}_j-\mathbf{r}_i$, $\mathbf{b}_2=\mathbf{r}_k-\mathbf{r}_j$,
$\mathbf{b}_3=\mathbf{r}_l-\mathbf{r}_k$, and the plane normals
$\mathbf{n}_1=\mathbf{b}_1\times\mathbf{b}_2$, $\mathbf{n}_2=\mathbf{b}_2\times\mathbf{b}_3$:

$$
\phi_{ijkl} = \frac{180}{\pi}\,
\operatorname{atan2}\!\Big(
\big(\hat{\mathbf{n}}_1\times\hat{\mathbf{b}}_2\big)\cdot\hat{\mathbf{n}}_2,\;
\hat{\mathbf{n}}_1\cdot\hat{\mathbf{n}}_2
\Big),
\qquad
\hat{\mathbf{u}}\equiv\frac{\mathbf{u}}{\lVert\mathbf{u}\rVert}
$$

- **Represents:** the **signed** torsion angle in the range $(-180, 180]$. The
  `atan2` form gives the correct sign (chirality) directly, unlike a bare
  `arccos` of the normals' dot product.
- **Where/how used:** the vectorized helper `_calculate_dihedrals` computes this
  for every i-j-k-l built around each central bond $j$–$k$ in the current frame.
  A degenerate geometry (overlapping central atoms $\lVert\mathbf{b}_2\rVert=0$, or a
  collinear triple giving $\lVert\mathbf{n}_1\rVert=0$ or $\lVert\mathbf{n}_2\rVert=0$)
  yields `NaN` and is skipped for that frame. The convention matches the single
  **Dihedral angle analysis** tool and is recorded in the output header as
  `# dihedral_convention signed_degrees_minus180_to_180`.
- **References:** A. Blondel and M. Karplus, "New formulation for derivatives of
  torsion angles and improper torsion angles in molecular mechanics: Elimination of
  singularities", *J. Comput. Chem.* **17**(9), 1132–1141 (1996); the `atan2`
  construction is the standard singularity-free torsion definition (also widely known
  as the "Praxeolitic" formula).

## Formulas and Methods (Bond Length / Bond Angle / Dihedral Angle tools)

This section documents every equation implemented in the three single-item
geometry tools — [bond.py](bond.py) (one atom pair), [bondAngle.py](bondAngle.py)
(one i-j-k triplet), and [dihedralAngle.py](dihedralAngle.py) (one i-j-k-l
quadruplet). These tools track **one** selected item across the whole trajectory,
build its distribution, and derive a free-energy profile from it. Everything they
share lives in [analysisCommon.py](analysisCommon.py), so equations 1, 2, 6, 7, 8,
10, and 11 below are literally the same code in all three tools.

Notation: $\mathbf{r}_a$ is the Cartesian position of atom $a$ in a given frame,
$\mathbf{L}=(a,b,c)$ are the orthorhombic cell lengths, $x$ stands for whichever
coordinate the tool measures ($r$, $\theta$, or $\phi$), and $N$ is the number of
frames actually used.

> **Note on the sibling section.** The [All Bond / Angle / Dihedral
> tools](#formulas-and-methods-all-bond--angle--dihedral-tools) above share the
> geometric definitions (equations 3–5 here match equations 2, 5, and 6 there) but
> differ in two ways: they use a **running Welford population** variance
> ($\div N$), whereas these tools use the **sample** variance ($\div (N-1)$, see
> equation 7); and they have no free-energy step at all.

### 1. Physical constants and the simulation time axis

$$
\tau_\text{a.u.} = 0.02418884326505\ \text{fs},
\qquad
R = 1.987204\times10^{-3}\ \text{kcal}\,\text{mol}^{-1}\text{K}^{-1}
$$

$$
\Delta t_\text{ps}
= \frac{\Delta t_\text{a.u.}\; N_\text{samp}\; \tau_\text{a.u.}}{1000},
\qquad
t_n = n\,\Delta t_\text{ps},\quad n = 0,1,\dots,N-1
$$

- **Represents:** the conversion from CPMD's internal units to picoseconds.
  $\Delta t_\text{a.u.}$ is the **Simulation Time Step** field and
  $N_\text{samp}$ the **Sampling Interval (frames)** field — the number of MD
  steps between two *stored* frames — so their product is the real time between
  consecutive frames of `TRAJEC.xyz`. Note that the sampling interval only
  rescales the time axis; it does **not** subsample the trajectory.
- **Where/how used:** `GeometryAnalysisBase._time_increment_ps` and the class
  attributes `atufs` and `gas_constant`. The time column of every `.dat`/`.csv`
  time series is $t_n$; $R$ appears in equations 8 and 9.
- **References:** CODATA recommended values (atomic unit of time,
  $2.4188843265\times10^{-17}$ s; molar gas constant
  $8.314462618\ \text{J}\,\text{mol}^{-1}\text{K}^{-1}$, converted with
  $1\ \text{cal}=4.184\ \text{J}$).

### 2. Interatomic vector and minimum-image convention (optional PBC)

$$
\Delta\mathbf{r}_{ab} = \mathbf{r}_a - \mathbf{r}_b,
\qquad
\Delta\mathbf{r}_{ab} \;\leftarrow\; \Delta\mathbf{r}_{ab}
- \mathbf{L}\,\operatorname{round}\!\left(\frac{\Delta\mathbf{r}_{ab}}{\mathbf{L}}\right)
$$

- **Represents:** the shortest vector between two atoms under periodic boundary
  conditions, applied component-wise. It lets an item that spans a box face be
  measured correctly and tolerates unwrapped coordinates.
- **Where/how used:** `GeometryAnalysisBase.minimum_image`. Applied to the single
  pair vector in the bond tool, to **both arms** $\mathbf{v}_{ji},\mathbf{v}_{jk}$
  in the angle tool, and to **all three** bond vectors
  $\mathbf{b}_1,\mathbf{b}_2,\mathbf{b}_3$ in the dihedral tool. When **Cell
  lengths a b c** is left blank the subtraction is skipped
  (`cell_lengths is None`) and raw vectors are used.
- **Validity:** the convention is exact only for true separations up to
  $\tfrac{1}{2}\min(a,b,c)$. The bond tool *enforces* this by requiring
  $r_\text{max}\le\tfrac{1}{2}\min(a,b,c)$ in `read_params`. The angle and
  dihedral tools have no equivalent range parameter to bound, so they instead
  **warn** when a computed arm or bond exceeds that limit.
- **References:** M. P. Allen and D. J. Tildesley, *Computer Simulation of
  Liquids*, 2nd ed., Oxford University Press (2017), §1.6.3; D. Frenkel and
  B. Smit, *Understanding Molecular Simulation*, 2nd ed., Academic Press (2002).

### 3. Bond length (bond tool)

$$
r_n = \lVert \Delta\mathbf{r}_{ab} \rVert
    = \sqrt{\Delta x_{ab}^{2} + \Delta y_{ab}^{2} + \Delta z_{ab}^{2}}
$$

- **Represents:** the Euclidean distance in ångström between the two selected
  atoms in frame $n$, using the minimum-image vector of equation 2.
- **Where/how used:** `BondAnalyser.bond_length`, computed for all frames at once
  as an $(N,2,3)$ array. Written to `bond_<pair>.dat` as $(t_n, r_n)$.
- **References:** standard Euclidean norm; Allen and Tildesley (2017), Ch. 1.

### 4. Bond angle (angle tool)

For the triplet i-j-k with vertex $j$, using the minimum-image arms
$\mathbf{v}_{ji}=\mathbf{r}_i-\mathbf{r}_j$ and
$\mathbf{v}_{jk}=\mathbf{r}_k-\mathbf{r}_j$:

$$
\theta_n = \frac{180}{\pi}\,
\arccos\!\left(
\operatorname{clip}\!\left(
\frac{\mathbf{v}_{ji}\cdot\mathbf{v}_{jk}}
     {\lVert\mathbf{v}_{ji}\rVert\,\lVert\mathbf{v}_{jk}\rVert},\,
-1,\,1\right)\right)
\;\in\;[0,180]
$$

- **Represents:** the bond angle in degrees at the central atom. The clip to
  $[-1,1]$ guards against tiny floating-point overshoot of the cosine before
  `arccos`.
- **Where/how used:** `BondAngleAnalyser.bond_angle`, vectorized over all frames.
  A **degenerate** frame — either arm of zero length, i.e. an atom sitting on the
  vertex — has no defined angle; it is dropped, counted in `skipped_degenerate`,
  and reported in the summary. The time axis is taken from the original frame
  index, so dropping a frame does not shift later times.
- **References:** dot-product angle definition; Allen and Tildesley (2017), §1.3.
  Numerical clipping of the `arccos` argument is standard practice.

### 5. Dihedral (torsion) angle and range convention (dihedral tool)

For the quadruplet i-j-k-l, with minimum-image bond vectors
$\mathbf{b}_1=\mathbf{r}_j-\mathbf{r}_i$, $\mathbf{b}_2=\mathbf{r}_k-\mathbf{r}_j$,
$\mathbf{b}_3=\mathbf{r}_l-\mathbf{r}_k$ and plane normals
$\mathbf{n}_1=\mathbf{b}_1\times\mathbf{b}_2$,
$\mathbf{n}_2=\mathbf{b}_2\times\mathbf{b}_3$:

$$
\phi_n = \frac{180}{\pi}\,
\operatorname{atan2}\!\Big(
\big(\hat{\mathbf{n}}_1\times\hat{\mathbf{b}}_2\big)\cdot\hat{\mathbf{n}}_2,\;
\hat{\mathbf{n}}_1\cdot\hat{\mathbf{n}}_2
\Big)\;\in\;(-180,180],
\qquad
\hat{\mathbf{u}}\equiv\frac{\mathbf{u}}{\lVert\mathbf{u}\rVert}
$$

The reported value then depends on the **Wrap dihedral to [-180, 180]** switch:

$$
\phi_n' =
\begin{cases}
\phi_n, & \text{wrap on}\;\Rightarrow\;\phi' \in (-180,180], \\[4pt]
\phi_n + 360\;\text{if}\;\phi_n<0,\;\text{else}\;\phi_n,
  & \text{wrap off}\;\Rightarrow\;\phi' \in [0,360).
\end{cases}
$$

- **Represents:** the **signed** torsion angle. The `atan2` form gives the correct
  sign (chirality) directly, unlike a bare `arccos` of the normals' dot product,
  and is free of the singularity that afflicts the naive formulation.
- **Where/how used:** `DihedralAngleAnalyser.dihedral_angle`, vectorized over all
  frames. A degenerate geometry — $\lVert\mathbf{b}_2\rVert=0$, or a collinear
  triple giving $\lVert\mathbf{n}_1\rVert=0$ or $\lVert\mathbf{n}_2\rVert=0$ — has
  no dihedral plane; those frames are dropped and counted. **The wrap switch also
  sets the histogram range of equation 6**: $[-180,180]$ when on, $[0,x_\text{max}]$
  when off. This coupling matters — a $[0,x_\text{max}]$ histogram applied to
  signed angles would give every negative angle a negative bin index and silently
  discard it.
- **References:** A. Blondel and M. Karplus, "New formulation for derivatives of
  torsion angles and improper torsion angles in molecular mechanics: Elimination
  of singularities", *J. Comput. Chem.* **17**(9), 1132–1141 (1996); the `atan2`
  construction is the standard singularity-free torsion definition.

### 6. Histogram binning and the probability distribution

With bin width $\Delta x$ over the range $[x_\text{min}, x_\text{max})$:

$$
N_b = \left\lfloor \frac{x_\text{max}-x_\text{min}}{\Delta x} \right\rfloor,
\qquad
k_n = \left\lfloor \frac{x_n - x_\text{min}}{\Delta x} \right\rfloor,
\qquad
x_k = x_\text{min} + \Delta x\left(k + \tfrac{1}{2}\right)
$$

A sample is retained only when $0 \le k_n < N_b$; the rest are counted as
out-of-range. With $n_k$ the count in bin $k$ and
$N_\text{in} = \sum_{k} n_k$ the number retained:

$$
P_k = \frac{n_k}{N_\text{in}},
\qquad
P_k^{\%} = 100\,P_k,
\qquad
N_\text{out} = N - N_\text{in}
$$

- **Represents:** the normalized probability of finding the coordinate in bin $k$,
  reported as a percentage. $x_k$ is the **bin centre**, which is the abscissa
  written to the distribution and free-energy files.
- **Where/how used:** `GeometryAnalysisBase.histogram_percent`. The ranges are
  $[0, r_\text{max})$ for the bond tool, $[0, \theta_\text{max})$ for the angle
  tool, and either $[-180,180)$ or $[0,\phi_\text{max})$ for the dihedral tool
  (equation 5).
- **Important:** $P_k$ is normalized over the **retained** samples only, so
  out-of-range values bias both the distribution and the free energy derived from
  it. All three tools therefore report $N_\text{out}$ in a dialog and in the
  summary. Out-of-range values are **excluded, never clamped** into the end bin —
  clamping would silently inflate the last bar.
- **References:** standard frequency-histogram estimator of a probability density;
  Frenkel and Smit (2002), Ch. 4.

### 7. Sample mean, variance, and standard deviation

$$
\bar{x} = \frac{1}{N}\sum_{n=1}^{N} x_n,
\qquad
s^{2} = \frac{1}{N-1}\sum_{n=1}^{N}\left(x_n-\bar{x}\right)^{2},
\qquad
s = \sqrt{s^{2}}
$$

- **Represents:** the average value of the coordinate over the trajectory and its
  spread. The $N-1$ denominator is **Bessel's correction**, giving the unbiased
  *sample* variance.
- **Where/how used:** `GeometryAnalysisBase.basic_stats`, implemented as
  `numpy.var(..., ddof=1)` / `numpy.std(..., ddof=1)`. Reported in the
  `RESULTS` block of the summary and in `summary_*.txt`. For $N=1$ both are
  defined as $0$.
- **Note:** this differs deliberately from the **population** variance ($\div N$)
  used by the Welford accumulator in the *All* tools — those average over a
  varying set of frames per item, these over one fixed series.
- **References:** Bessel's correction, standard in any statistics text; the
  convention matches Python's `statistics.variance` and `statistics.stdev`, which
  these tools previously used.

### 8. Free energy from the probability distribution

$$
G(x_k) = -RT\,\ln P_k
$$

- **Represents:** the free-energy profile along the chosen coordinate, in
  kcal/mol, obtained by inverting the Boltzmann relation
  $P \propto \exp(-G/RT)$. Only **populated** bins ($P_k > 0$) contribute, since
  $\ln 0$ is undefined; empty bins leave gaps in the curve rather than spikes.
- **Where/how used:** `GeometryAnalysisBase.free_energy_from_histogram`, called by
  all three tools. $T$ is the **Simulation Temperature (K)** field. The minimum of
  the curve and the coordinate at which it occurs are reported in the summary.
- **Important:** $G$ is defined only up to an **additive constant**, and its
  absolute value depends on the bin width through $P_k$. Compare *shapes* and
  *barrier heights* between runs, not absolute numbers — and only between runs
  that used the same bin width.
- **Numerical guard:** probabilities are floored at $10^{-12}$ before the
  logarithm to avoid $-\infty$ from underflow.
- **References:** D. Chandler, *Introduction to Modern Statistical Mechanics*,
  Oxford University Press (1987), Ch. 7; Frenkel and Smit (2002), Ch. 7;
  M. E. Tuckerman, *Statistical Mechanics: Theory and Molecular Simulation*,
  Oxford University Press (2010), Ch. 8. See also J. Mol. Model. **17**, 2159–2168
  (2011), DOI 10.1007/s00894-010-0939-6, and *Science* **275**, 817 (1997),
  DOI 10.1126/science.275.5301.817, the references cited in the tools' own help
  text.

### 9. Jacobian (volume-element) correction and the potential of mean force

The three coordinates do **not** share the same measure. Their volume elements are

$$
dV = 4\pi r^{2}\,dr
\quad\text{(distance)},
\qquad
d\Omega \propto \sin\theta\,d\theta
\quad\text{(angle)},
\qquad
d\phi
\quad\text{(dihedral, uniform)}.
$$

Dividing out that geometric weight turns the raw profile of equation 8 into the
potential of mean force $W$:

$$
W(r) = -RT\,\ln\frac{P(r)}{r^{2}},
\qquad
W(\theta) = -RT\,\ln\frac{P(\theta)}{\sin\theta},
\qquad
W(\phi) = -RT\,\ln P(\phi).
$$

- **Represents:** the free energy with the *entropic* contribution of the
  coordinate's own volume element removed. Without it, a distance profile shows an
  apparent bias toward large $r$ and an angle profile toward $90°$, purely because
  there is more phase space there — not because of any interaction.
- **Where/how used:** the `jacobian` argument of
  `GeometryAnalysisBase.free_energy_from_histogram`.
  - **bond.py** — the **Jacobian r² correction (PMF)** switch, off by default.
  - **bondAngle.py** — the **Use sin(θ) Jacobian** switch, off by default.
  - **dihedralAngle.py** — **no switch, deliberately.** A dihedral has a uniform
    measure, so there is nothing to divide out and $-RT\ln P(\phi)$ is *already*
    the potential of mean force. The summary records this as
    `Jacobian: not applicable (uniform measure)`.
- **Numerical guard:** the Jacobian factor is floored at $10^{-12}$, which matters
  for $\sin\theta$ as $\theta\to0°$ or $180°$.
- **Relation to $g(r)$:** for a distance this is the same construction that links
  the radial distribution function to the PMF,
  $W(r) = -RT\ln g(r)$ with $P(r)\propto 4\pi r^{2} g(r)$ — see the
  **Radial Distribution Function** tool.
- **References:** Chandler (1987), §7.3; Tuckerman (2010), Ch. 8;
  D. Trzesniak, A.-P. E. Kunz and W. F. van Gunsteren, "A comparison of methods to
  compute the potential of mean force", *ChemPhysChem* **8**(1), 162–169 (2007).

### 10. Moving-average smoothing of the probabilities (optional)

With an odd window $w$ and edge padding (the first and last values are repeated
so the smoothed series keeps its length):

$$
\tilde{P}_k = \frac{1}{w}\sum_{m=-(w-1)/2}^{(w-1)/2} P_{k+m},
\qquad
w = \max\!\left(5,\;\min\!\left(11,\;\left\lfloor \frac{N_p}{12}\right\rfloor\right)\right)
$$

where $N_p$ is the number of populated bins; $w$ is incremented by one if even,
capped at $N_p$, and set to $1$ (no smoothing) when $N_p < 5$.

- **Represents:** a gentle low-pass filter applied to the probabilities **before**
  the logarithm of equation 8, which tames the noise that sparse bins produce at
  fine bin widths.
- **Where/how used:** `GeometryAnalysisBase.smooth_series`, enabled by the
  **Smooth free energy** switch. Default **on** for bondAngle.py and
  dihedralAngle.py and **off** for bond.py, which preserves each tool's historical
  output. The summary records the setting as `Smoothed FE`.
- **Caution:** smoothing is cosmetic, not physical. It broadens sharp minima and
  slightly shifts the reported free-energy minimum, so turn it off when quoting
  numbers and prefer a coarser bin width to a heavier filter.
- **References:** moving-average (boxcar) filter, standard signal processing; for
  the trade-off between binning and smoothing in free-energy estimates see
  Frenkel and Smit (2002), Ch. 7.

### 11. Quantities reported in the summary

Collecting the above, each run ends with:

$$
\bar{x},\quad s^{2},\quad s,\quad x_\text{min},\quad x_\text{max},
\qquad
G_\text{min} = \min_k G(x_k),
\qquad
x^{*} = \arg\min_{x_k} G(x_k)
$$

- **Represents:** the average, spread, and extremes of the coordinate, plus the
  depth and position of the deepest free-energy well.
- **Where/how used:** `GeometryAnalysisBase.render_summary` lays out the
  `PARAMETERS` and `RESULTS` blocks identically for all three tools; the same text
  is shown in the output box and written to `summary_*.txt`. When the Jacobian is
  active the label reads `Lowest PMF` instead of `Lowest free energy`.
- **References:** as for equations 7–9.

## Formulas and Methods (Mean Residence Time tool)

This section documents every equation implemented in
[meanResidenceTime.py](meanResidenceTime.py), in the order the tool evaluates
them. The logic lives in `MeanResidenceTimeCalculator`; the class is pure
NumPy and is exercised headlessly by
[tests/test_meanResidenceTime.py](tests/test_meanResidenceTime.py), whose
expected values are the analytic results quoted below.

> **Note on the two continuous time constants.** Equations 5 and 6 define two
> *different* quantities — the mean residence time $\langle T\rangle$ and the
> mean residual time $\tau_\text{IMM}$. They coincide only for exponential
> kinetics. Reporting one as though it were the other is the most common error
> in this analysis, so the tool prints both, separately labelled.

**Symbols used throughout.**

| Symbol | Meaning |
| --- | --- |
| $N$ | number of frames in the trajectory |
| $M$ | number of observed objects (atoms or groups) |
| $\Delta t$ | time between saved frames |
| $r_\text{cut}$ | cutoff radius defining the shell |
| $\mathbf{L}=(a,b,c)$ | orthorhombic cell edges (optional) |
| $\mathbf{R}(t)$ | position of the reference point in frame $t$ |
| $\mathbf{x}_i(t)$ | position of observed object $i$ in frame $t$ |
| $h_i(t)$ | shell indicator (equation 3) |
| $t^{*}$ | tolerance time, in frames |
| $n_j$ | length, in frames, of residence event $j$ |

### 1. Reference and observed positions

$$
\mathbf{X}_\text{geom}(t) = \frac{1}{n}\sum_{k=1}^{n}\mathbf{r}_k(t),
\qquad
\mathbf{X}_\text{com}(t) =
  \frac{\sum_{k=1}^{n} m_k\,\mathbf{r}_k(t)}{\sum_{k=1}^{n} m_k}
$$

- **Represents:** the collapse of a group of $n$ atoms to a single tracked point
  per frame — either the unweighted centroid or the mass-weighted centre of mass.
  For the `single_atom` and `atom_list_individual` modes the atom's own
  coordinates are used directly and neither formula is evaluated.
- **Where/how used:** `_geometric_center_series` and `_center_of_mass_series`,
  both vectorized over all frames at once (measured 73× faster than the
  per-frame loop they replaced, with identical results). Masses come from
  `AtomicData` in [help.py](help.py) and are looked up **only** for the two
  centre-of-mass modes, so a purely geometric run never fails on an element
  missing from the table.
- **References:** standard definitions; see M. P. Allen and D. J. Tildesley,
  *Computer Simulation of Liquids*, 2nd ed., Oxford University Press (2017), §2.
  Atomic masses follow the IUPAC 2021 standard atomic weights as tabulated in
  `AtomicData`.

### 2. Separation and the minimum-image convention (optional PBC)

$$
\mathbf{d}_i(t) = \mathbf{x}_i(t) - \mathbf{R}(t),
\qquad
\mathbf{d}_i(t) \;\leftarrow\; \mathbf{d}_i(t)
  - \mathbf{L}\,\operatorname{round}\!\left(\frac{\mathbf{d}_i(t)}{\mathbf{L}}\right)
$$

$$
D_i(t) = \lVert \mathbf{d}_i(t) \rVert
       = \sqrt{d_x^{2} + d_y^{2} + d_z^{2}}
$$

- **Represents:** the separation between the observed object and the reference
  point. The second expression is the orthorhombic minimum-image convention,
  applied component-wise: it replaces each separation by the shortest one among
  all periodic images. It is applied only when cell lengths are supplied.
- **Where/how used:** `build_occupancy`. The convention is exact only while
  $r_\text{cut} \le \tfrac{1}{2}\min(a,b,c)$, and the tool *enforces* that bound,
  refusing larger cutoffs in both `build_occupancy` and `_collect_inputs`.
  Coordinates need **not** be wrapped into the box — minimum image tolerates
  unwrapped trajectories, which is what CPMD normally writes.
- **References:** M. P. Allen and D. J. Tildesley, *Computer Simulation of
  Liquids*, 2nd ed., Oxford University Press (2017), §1.6.3.

### 3. Occupancy (the shell indicator)

$$
h_i(t) =
\begin{cases}
1, & D_i(t) \le r_\text{cut} \\[4pt]
0, & \text{otherwise}
\end{cases}
$$

- **Represents:** the binary "object $i$ is inside the shell in frame $t$"
  function. Every quantity below is computed from this $M \times N$ matrix and
  nothing else, so the choice of $r_\text{cut}$ propagates into all of them.
- **Where/how used:** `build_occupancy`, returned alongside the distances.
  Take $r_\text{cut}$ from the first minimum of the corresponding $g(r)$; the
  [Radial Distribution Function](#radial-distribution-function) tool supplies it.
- **References:** R. W. Impey, P. A. Madden and I. R. McDonald, "Hydration and
  mobility of ions in solution", *J. Phys. Chem.* **87**(25), 5071–5083 (1983).

### 4. Tolerance time $t^{*}$ (recrossing correction)

$$
\ldots\,1\;1\;
\underbrace{0\;\cdots\;0}_{g\ \text{frames}}
\;1\;1\,\ldots
\;\longrightarrow\;
\ldots\,1\;1\;
\underbrace{1\;\cdots\;1}_{g\ \text{frames}}
\;1\;1\,\ldots
\qquad \text{if } g \le t^{*}
$$

- **Represents:** bridging short excursions. A gap of $g$ zero-frames flanked by
  ones on both sides is overwritten with ones when $g \le t^{*}$, so a molecule
  that momentarily jitters back and forth across the cutoff is not counted as
  having left and returned.
- **Where/how used:** `apply_tolerance_to_occupancy`, applied to the occupancy
  matrix before any statistic is taken. Gaps at the very **start or end** of the
  trajectory are deliberately never bridged, because there is no evidence the
  object was inside before or after the run. Because $t^{*}$ can only merge
  events, never split them, $\langle T\rangle$ is non-decreasing in $t^{*}$.
- **References:** Impey, Madden and McDonald (1983), as above, introduced $t^{*}$.
  D. Laage and J. T. Hynes, "On the residence time for water in a solute
  hydration shell", *J. Phys. Chem. B* **112**(26), 7697–7701 (2008), showed the
  result is *extremely* sensitive to $t^{*}$ and that the conventional
  $t^{*}=2$ ps can seriously overestimate residence times for low-barrier
  exchange — which is why the tool provides a $t^{*}$ scan rather than a single
  value.

### 5. Residence events, censoring, and the mean residence time

$$
T_j = n_j\,\Delta t,
\qquad
\langle T\rangle = \frac{1}{N_\text{ev}}\sum_{j=1}^{N_\text{ev}} T_j
$$

- **Represents:** each maximal run of ones in $h_i(t)$ is one residence event of
  duration $T_j$; $\langle T\rangle$ is their mean. This is the direct
  residence-time estimator, free of the assumptions built into the
  correlation-function routes.
- **Where/how used:** `_runs`, `extract_event_durations` and
  `count_boundary_events`. Runs touching frame $0$ or frame $N-1$ are
  **censored** — excluded from $\{T_j\}$ and reported separately — because they
  were still in progress at the boundary and their true duration is unknown.
  Counting them at their truncated length biases $\langle T\rangle$ downwards,
  and does so worst for the longest-lived species. If the censored count is a
  large fraction of the events, the trajectory is too short for the process
  being measured.
- **References:** H. R. Sánchez, "Residence Times from Molecular Dynamics
  Simulations", *J. Phys. Chem. B* **126**(20), 3843–3857 (2022).

### 6. Continuous survival function and the IMM time constant

$$
S(l) = \frac{\sum_{j}\max\!\left(0,\; n_j - l\right)}{\sum_{j} n_j},
\qquad t = l\,\Delta t
$$

$$
\tau_\text{IMM} = \Delta t \sum_{l=0}^{l_\text{max}-1} S(l)
\;=\; \frac{\langle T^{2}\rangle}{2\,\langle T\rangle}
$$

- **Represents:** $S(l)$ is the probability that an object found inside the shell
  at an arbitrary **occupied frame** is still *continuously* inside $l$ frames
  later. The average runs over time **origins**, not over events, so a long event
  contributes as many origins as it has frames. This is the $n_\text{ion}(t)$ of
  Impey, Madden and McDonald. Its integral is mathematically the mean *residual*
  time $\langle T^{2}\rangle / (2\langle T\rangle)$, which equals
  $\langle T\rangle$ **only** for exponential kinetics and exceeds it whenever the
  duration distribution is broad — the usual case for hydration shells.
- **Where/how used:** `compute_origin_averaged_survival`. The closed form above,
  summing over run lengths $n_j$, is exactly equivalent to scanning every time
  origin explicitly (a unit test checks it against a brute-force origin scan) but
  costs $\mathcal{O}(\text{events})$ per lag instead of $\mathcal{O}(N^{2})$.
  For a set of uniform runs of $n$ frames the integral is exactly
  $\tfrac{1}{2}(n+1)\,\Delta t$, which is the analytic value the tests assert.
- **References:** Impey, Madden and McDonald (1983); A. E. García and L. Stiller,
  "Computation of the mean residence time of water in the hydration shells of
  biomolecules", *J. Comput. Chem.* **14**(11), 1396–1406 (1993); B. Halle and
  F. Persson, "Analysis of Protein Dynamics Simulations by a Stochastic Point
  Process Approach", *J. Chem. Theory Comput.* **9**(6), 2838–2848 (2013), for
  computing the survival function from the residence-time sequence rather than by
  averaging over initial times; S. Mukherjee, S. Mondal and B. Bagchi,
  *J. Chem. Phys.* **147**, 024901 (2017), for the broad, log-normal duration
  distributions that make the two time constants diverge.

### 7. Intermittent correlation function

$$
C(l) = \frac{\sum_{i=1}^{M}\sum_{t_0=0}^{N-l-1} h_i(t_0)\,h_i(t_0+l)}
            {\sum_{i=1}^{M}\sum_{t_0=0}^{N-l-1} h_i(t_0)},
\qquad
\langle h\rangle = \frac{1}{MN}\sum_{i=1}^{M}\sum_{t=0}^{N-1} h_i(t)
$$

$$
R(l) = \frac{C(l) - \langle h\rangle}{1 - \langle h\rangle}
$$

- **Represents:** $C(l)$ is the probability that an object inside the shell at
  $t_0$ is inside again at $t_0+l$, **whether or not it left in between**. It
  therefore probes exchange rather than uninterrupted residence, and is the
  natural partner to equation 6. Since $C(l)\to\langle h\rangle$ at long lag,
  $R(l)$ is the rescaling running from $R(0)=1$ to $R(\infty)=0$, which is what
  can be integrated to a time constant.
- **Where/how used:** `compute_intermittent_functions`. Note the asymmetry:
  $C(l)$ is normalized over the origins available *at that lag*, while
  $\langle h\rangle$ is the global mean occupancy over all frames, so the two
  denominators differ slightly. When $\langle h\rangle = 1$ the object never
  leaves and $R$ is defined as $1$ for all lags. $C$ is clipped to $[0,1]$ and
  $R$ to $[-1,1]$ to absorb floating-point overshoot.
- **References:** Impey, Madden and McDonald (1983). The same
  intermittent/continuous pairing underlies the hydrogen-bond kinetics of
  A. Luzar and D. Chandler, "Hydrogen-bond kinetics in liquid water", *Nature*
  **379**, 55–57 (1996).

### 8. Lag ceiling

$$
l_\text{max} =
\max\!\left(1,\;
  \min\!\left(N,\; \left\lfloor \tfrac{N}{10} \right\rfloor \right)
\right)
$$

- **Represents:** the largest lag at which the correlation functions are
  evaluated, when the **Maximum lag** field is left blank. At lag $l$ only
  $N-l$ time origins are available, so the longest lags average over a handful
  of samples and are noise; integrating them into a residence time is
  meaningless.
- **Where/how used:** `default_max_lag`, applied to both equation 6 and
  equation 7. It can be overridden from the **Advanced** group on the Setup tab,
  which is how you check whether a reported time constant has converged.
- **References:** Halle and Persson (2013), for the growth of statistical error
  with lag in this class of correlation function.

### 9. Turning $R(t)$ into a time constant

$$
\text{zero crossing:}\quad
\tau = \Delta t \sum_{l=0}^{l_0-1} R(l),
\qquad l_0 = \min\{\, l : R(l) < 0 \,\}
$$

$$
\text{exponential:}\quad
\ln R(l) = A - \frac{l\,\Delta t}{\tau}
\quad\text{(least squares over } \{\, l : R(l) > 0 \,\}\text{)},
\qquad \tau = -\frac{1}{\text{slope}}
$$

$$
\text{full:}\quad
\tau = \Delta t \sum_{l=0}^{l_\text{max}-1} R(l)
$$

- **Represents:** three ways to reduce the decaying correlation to a single
  time. Zero crossing integrates up to the point the curve first goes negative;
  the exponential rule assumes single-exponential decay and reports the fitted
  $\tau$; full integrates the whole evaluated range.
- **Where/how used:** `integrate_correlation`, selected by the **Integration
  rule** control. The exponential fit requires at least three lags with
  $R(l)>0$ and falls back to the zero-crossing rule if the fitted slope is not
  negative. The summary always records which rule was applied and how many lags
  it covered, so a run is reproducible from its own output.
- **Convergence check:** if the number of lags used **equals** $l_\text{max}$,
  then $R(t)$ never crossed zero, and $\tau$ is a lower bound that grows with
  $l_\text{max}$ rather than a converged time constant. Raise the lag ceiling and
  see whether the value is stable. If it keeps growing, the trajectory is too
  short to resolve the process and the intermittent time constant should not be
  quoted at all — report $\langle T\rangle$ and $\tau_\text{IMM}$ instead.
- **References:** the exponential form is the standard single-relaxation-time
  reading of the intermittent correlation function used by Impey, Madden and
  McDonald (1983) and, for hydrogen-bond kinetics, by A. Luzar and D. Chandler,
  *Nature* **379**, 55–57 (1996). H. R. Sánchez, *J. Phys. Chem. B* **126**,
  3843 (2022), analyses the assumptions this step introduces relative to the
  direct estimator of equation 5.

### 10. Uncertainty on the mean residence time

$$
\sigma\!\left(\langle T\rangle\right)
= \frac{1}{\sqrt{B}}
  \sqrt{\frac{1}{B-1}\sum_{b=1}^{B}\left(\bar{T}_b - \overline{\bar{T}}\right)^{2}},
\qquad B = 5
$$

- **Represents:** the standard error of $\langle T\rangle$, estimated from the
  means $\bar{T}_b$ of $B$ contiguous blocks of the event-duration series.
- **Where/how used:** `block_average_uncertainty`. Block averaging is used rather
  than a plain standard error because successive residence events are **not**
  independent samples; a naive $s/\sqrt{n}$ would understate the error. Returns
  NaN when there is not enough data to form two blocks. Only $\langle T\rangle$
  carries an uncertainty — $\tau_\text{IMM}$ and the intermittent time constant
  do not.
- **References:** H. R. Sánchez et al., "Evaluating the uncertainty in mean
  residual times: Estimators based on residence times from discrete time
  processes", *Commun. Nonlinear Sci. Numer. Simul.* (2024).

## Formulas and Methods (Electrostatic Surface Potential tool)

This section documents the theory and every algorithm implemented by the
[Electrostatic Surface Potential (ESP)](#electrostatic-surface-potential-esp)
tool, in the order the tool evaluates them. The code is split by
responsibility: [gaussianWavefunction.py](gaussianWavefunction.py) (file,
basis set, density), [espIntegrals.py](espIntegrals.py) (the potential;
[espNumba.py](espNumba.py) is its optional compiled twin),
[molecularSurface.py](molecularSurface.py) (surface and extrema) and
[espViewer.py](espViewer.py) (colours). None of the first four uses the GUI, and
all are exercised headlessly by the `tests/test_gaussianWavefunction.py`,
`test_espIntegrals.py`, `test_molecularSurface.py`, `test_espViewer.py` and
`test_electrostaticSurfacePotential.py` suites against Gaussian 09 `cubegen`
reference values.

**Units.** All quantities are computed in Hartree atomic units
($\hbar = m_e = e = 4\pi\varepsilon_0 = 1$): lengths in bohr, charges in units of
the elementary charge, the electron density in e/bohr³ and the potential in
Hartree/e. Conversion to Å and kcal/mol happens only when results are reported
(equation 12).

**Symbols used throughout.**

| Symbol | Meaning |
| --- | --- |
| $\mathbf{r}$ | a point in space |
| $\mathbf{R}_A,\ Z_A$ | position and charge of nucleus $A$ (the fchk `Nuclear charges`) |
| $\rho(\mathbf{r})$ | total electron density (α + β) |
| $\chi_\mu(\mathbf{r})$ | basis function $\mu$ |
| $P_{\mu\nu}$ | total density matrix (fchk `Total SCF Density`) |
| $V(\mathbf{r})$ | molecular electrostatic potential |
| $\rho_\text{iso}$ | density isovalue that defines the surface |
| $h$ | grid spacing (bohr) |
| $p,\ \mathbf{P}$ | exponent and centre of a Gaussian product (equation 4) |
| $F_n(T)$ | Boys function of order $n$ (equation 6) |

### 1. The molecular electrostatic potential

$$
V(\mathbf{r}) \;=\; \sum_{A} \frac{Z_A}{\lvert \mathbf{r}-\mathbf{R}_A\rvert}
\;-\; \int \frac{\rho(\mathbf{r}')}{\lvert \mathbf{r}-\mathbf{r}'\rvert}\, d^3 r'
$$

- **Represents:** the electrostatic potential created at $\mathbf{r}$ by the
  molecule's nuclei (positive term) and electrons (negative term). Equivalently,
  $V(\mathbf{r})$ is the interaction energy of a unit positive point charge
  placed at $\mathbf{r}$ with the molecule, so in kcal/mol it reads as kcal/mol
  per elementary charge. $V>0$ where the nuclei dominate (regions attractive to
  nucleophiles and to the negative end of a partner), $V<0$ where the electrons
  dominate (lone pairs, π clouds; attractive to electrophiles).
- **Assumptions:** (i) Born–Oppenheimer clamped nuclei, treated as point
  charges; (ii) the potential of the *unperturbed* molecule — the test charge
  does not polarise it, which is what makes $V$ a property of the molecule
  alone; (iii) $\rho$ is the SCF (Hartree–Fock or Kohn–Sham) total density
  stored in the file; (iv) the potential is evaluated in vacuum: if the
  Gaussian job used an implicit solvent (SCRF/PCM), the stored density is the
  solvent-polarised one, but the reaction-field potential of the continuum is
  **not** added; (v) with an effective core potential, $Z_A$ is the effective
  nuclear charge (atomic number minus core electrons) and $\rho$ lacks the core
  electrons — a frozen-core point-charge model that is accurate on the
  molecular surface, far from the core, and is also Multiwfn's choice.
- **Where/how used:** `espIntegrals.electrostatic_potential` returns
  `nuclear_potential − electronic_potential`. The nuclear term is the explicit
  sum; the electronic term is evaluated analytically (equations 2–6).
- **References:** J. S. Murray and P. Politzer, "The electrostatic potential:
  an overview", *WIREs Comput. Mol. Sci.* **1**, 153–163 (2011);
  P. Politzer and D. G. Truhlar (eds.), *Chemical Applications of Atomic and
  Molecular Electrostatic Potentials*, Plenum (1981).

### 2. Electron density from the density matrix

$$
\rho(\mathbf{r}) = \sum_{\mu\nu} P_{\mu\nu}\,\chi_\mu(\mathbf{r})\,\chi_\nu(\mathbf{r}),
\qquad
P_{\mu\nu} = P^{\alpha}_{\mu\nu} + P^{\beta}_{\mu\nu}
$$

so that the electronic term of equation 1 becomes

$$
\int \frac{\rho(\mathbf{r}')}{\lvert \mathbf{r}-\mathbf{r}'\rvert}\, d^3 r'
= \sum_{\mu\nu} P_{\mu\nu}
  \int \frac{\chi_\mu(\mathbf{r}')\,\chi_\nu(\mathbf{r}')}{\lvert \mathbf{r}-\mathbf{r}'\rvert}\, d^3 r' .
$$

- **Represents:** the density as a quadratic form in the basis functions.
  Using the total density matrix makes closed- and open-shell files identical
  to handle; the spin density is not needed for the ESP.
- **Where/how used:** `gaussianWavefunction.read_fchk` reads
  `Total SCF Density` (lower triangle, row by row) and symmetrises it;
  `electron_density` evaluates $\rho$ on points in chunks of 20 000 as
  $\sum_i (\boldsymbol{\chi}^\mathsf{T}\mathbf{P})_i\,\chi_i$. Check:
  $\operatorname{tr}(\mathbf{P}\mathbf{S}) = N_\text{electrons}$ with the
  analytic overlap matrix $\mathbf{S}$ (`overlap_matrix`), satisfied to 10⁻⁸
  for every test file.
- **References:** A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry*,
  Dover (1996), chapter 3.

### 3. Gaussian basis functions

A contracted Cartesian Gaussian centred on $\mathbf{A}$, with
$\mathbf{r}_A = \mathbf{r}-\mathbf{A} = (x_A, y_A, z_A)$:

$$
\chi(\mathbf{r}) = \sum_{k} d_k\, N(\alpha_k; i,j,l)\;
x_A^{\,i}\, y_A^{\,j}\, z_A^{\,l}\; e^{-\alpha_k r_A^2},
\qquad
N(\alpha; i,j,l) = \left(\frac{2\alpha}{\pi}\right)^{3/4}
\frac{(4\alpha)^{L/2}}{\sqrt{(2i-1)!!\,(2j-1)!!\,(2l-1)!!}}
$$

with $L = i+j+l$, contraction coefficients $d_k$ and exponents $\alpha_k$
taken from the fchk. A spherical (pure) function is a fixed combination of the
Cartesian components of the same shell,

$$
\chi_{L m} = \sum_{c} T_{m c}\,\chi_c ,
\qquad
\mathbf{P}^\text{Cart} = \mathbf{T}^\mathsf{T}\,\mathbf{P}\,\mathbf{T},
$$

so every integral is done once, in the Cartesian basis.

- **Represents:** the basis exactly as Gaussian defines it. Each Cartesian
  primitive is normalised on its own and the fchk contraction coefficients are
  applied as written (the contracted functions then come out normalised, which
  `test_trace_ps_is_the_electron_count` confirms to 10⁻⁷). An `SP` shell is an s
  and a p shell sharing exponents, the p part using the
  `P(S=P) Contraction coefficients`.
- **Where/how used:** `primitive_norm`, `component_angular_norm`,
  `pure_to_cartesian` (real solid harmonics, rows renormalised numerically;
  Gaussian's pure order $0, +1, -1, +2, -2, \dots$) and `basis_transform`. The
  Cartesian component order follows the fchk (`FCHK_CARTESIAN_ORDER`); for f
  shells it is XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ, which differs
  from the `.wfn` order. s to g shells (l ≤ 4) are supported.
- **References:** T. Helgaker, P. Jørgensen and J. Olsen, *Molecular
  Electronic-Structure Theory*, Wiley (2000), chapter 6 (real solid
  harmonics); T. Lu and F. Chen, Multiwfn source `fileIO.f90` (fch reader).

### 4. Gaussian product theorem and Hermite expansion

$$
e^{-\alpha r_A^2}\, e^{-\beta r_B^2} = K_{AB}\, e^{-p\, r_P^2},
\qquad
p = \alpha+\beta,\quad
\mathbf{P} = \frac{\alpha\mathbf{A}+\beta\mathbf{B}}{p},\quad
K_{AB} = e^{-\frac{\alpha\beta}{p}\lvert\mathbf{A}-\mathbf{B}\rvert^2}
$$

The Cartesian factors are expanded in Hermite Gaussians, one direction at a
time, with the coefficients $E^{ij}_t$ obtained by recursion
($X_{PA} = P_x - A_x$, $X_{PB} = P_x - B_x$, $E^{00}_0 = 1$):

$$
E^{i+1,j}_{t} = \frac{1}{2p}E^{ij}_{t-1} + X_{PA}E^{ij}_{t} + (t+1)E^{ij}_{t+1},
\qquad
E^{i,j+1}_{t} = \frac{1}{2p}E^{ij}_{t-1} + X_{PB}E^{ij}_{t} + (t+1)E^{ij}_{t+1}
$$

- **Represents:** the product of two primitives as a single Gaussian at
  $\mathbf{P}$ times Hermite polynomials — the step that makes the
  three-dimensional integral of equation 2 analytic.
- **Where/how used:** `hermite_coefficients` (vectorised over all primitive
  pairs of a shell pair); `_primitive_pairs` forms $p$, $\mathbf{P}$ and
  $K_{AB}$. The same coefficients give the overlap matrix
  ($S = K_{AB}(\pi/p)^{3/2}E^{x}_0E^{y}_0E^{z}_0$) used as a check.
- **References:** L. E. McMurchie and E. R. Davidson, "One- and two-electron
  integrals over Cartesian Gaussian functions", *J. Comput. Phys.* **26**,
  218–231 (1978); Helgaker, Jørgensen and Olsen (2000), chapter 9.

### 5. Electronic potential by the McMurchie–Davidson scheme

$$
\int \frac{G_a(\mathbf{r}')\,G_b(\mathbf{r}')}{\lvert\mathbf{r}'-\mathbf{C}\rvert}\,d^3r'
= \frac{2\pi}{p}\,K_{AB}\sum_{tuv} E^{ij}_t\,E^{kl}_u\,E^{mn}_v\;R_{tuv}(p,\mathbf{P}-\mathbf{C})
$$

$$
R^{n}_{000} = (-2p)^n F_n\!\left(p\,\lvert\mathbf{P}-\mathbf{C}\rvert^2\right),
\qquad
R^{n}_{t+1,u,v} = t\,R^{n+1}_{t-1,u,v} + X_{PC}\,R^{n+1}_{t,u,v}
$$

(and likewise in $u$ with $Y_{PC}$, in $v$ with $Z_{PC}$); $R_{tuv} = R^{0}_{tuv}$.
Contracting the density matrix into the Hermite coefficients **once**,

$$
D^{(ab)}_{tuv} = \frac{2\pi}{p}\,K_{AB}\,c_a c_b
\sum_{\text{components}} P^\text{Cart}_{ab}\,E^{ij}_t E^{kl}_u E^{mn}_v,
\qquad
\int \frac{\rho(\mathbf{r}')}{\lvert\mathbf{r}'-\mathbf{C}\rvert}\,d^3r'
= \sum_{\text{pairs}}\ \sum_{tuv} D^{(ab)}_{tuv}\,R_{tuv}(p,\mathbf{P}-\mathbf{C})
$$

- **Represents:** the exact (analytic) electronic potential. $D_{tuv}$ does
  not depend on the evaluation point $\mathbf{C}$, so for $N$ surface points
  the work is $N \times$ (number of primitive pairs) $\times$ (number of
  $tuv$), with no grid or interpolation error.
- **Where/how used:** `build_hermite_batches` computes $D_{tuv}$ for the shell
  pairs $a \ge b$ (off-diagonal pairs doubled), drops a primitive pair when
  $\lvert \text{prefactor}\rvert \cdot \max\lvert P\rvert < 10^{-14}$
  (`PAIR_SCREEN`), and groups pairs by $L_a+L_b$; `hermite_integrals` runs the
  $R$ recursion keeping two levels of $n$ alive; `electronic_potential`
  evaluates blocks of at most 2×10⁶ pair–point elements; `electrostatic_potential`
  spreads point chunks over up to 12 threads. When `numba` is installed,
  `espNumba.electronic_potential` performs the same sum in one compiled loop
  (agreement 10⁻¹⁰ a.u.; about 5× faster: 20 000 points on benzene/6-31G(d)
  take 0.9 s instead of 4.7 s).
- **Validation:** against Gaussian 09 `cubegen` at 300 points 0.5–8 bohr from
  the nuclei, the largest difference is below 4×10⁻⁷ a.u. for 5D and 6D d,
  pure f and g, Cartesian f, SP shells and an open-shell file. Far from a
  neutral molecule $V$ tends to $\boldsymbol{\mu}\cdot\mathbf{r}/r^3$ with the
  fchk dipole, which is also tested.
- **References:** McMurchie and Davidson (1978); Helgaker, Jørgensen and Olsen
  (2000), chapter 9 (Coulomb integrals).

### 6. Boys function

$$
F_n(T) = \int_0^1 u^{2n}\, e^{-T u^2}\, du
= \frac{\gamma\!\left(n+\tfrac12,\,T\right)}{2\,T^{\,n+1/2}},
\qquad
F_{n-1}(T) = \frac{2T\,F_n(T) + e^{-T}}{2n-1}
$$

- **Represents:** the one-dimensional integral left after the angular
  integration of the Coulomb operator; $\gamma$ is the lower incomplete gamma
  function. The downward recursion is numerically stable.
- **Where/how used:** `espIntegrals.boys` computes the highest order from
  `scipy.special.gammainc`/`gamma` (with the series $1/(2n+1) - T/(2n+3)$ for
  $T < 10^{-10}$) and recurses downward. `espNumba._boys_into` uses the power
  series $F_n(T) = e^{-T}\sum_k (2T)^k / [(2n+1)(2n+3)\cdots(2n+2k+1)]$ for
  $T < 60$ and the asymptotic form $F_n(T) \approx \frac{(2n-1)!!}{2^{n+1}}\sqrt{\pi/T^{2n+1}}$
  above (at $T = 30$ the asymptotic $F_8$ is still 10⁻⁶ off, hence the limit of
  60). Both agree with numerical quadrature to 10⁻¹⁰ relative and with each
  other to 10⁻¹⁴.
- **References:** S. F. Boys, *Proc. R. Soc. Lond. A* **200**, 542–554 (1950);
  Helgaker, Jørgensen and Olsen (2000), chapter 9.

### 7. The molecular surface and its grid

$$
S(\rho_\text{iso}) = \{\, \mathbf{r} : \rho(\mathbf{r}) = \rho_\text{iso} \,\}
$$

- **Represents:** the molecular surface as an electron-density contour. The
  isovalue is a convention, not a physical constant: 0.001 a.u. is the
  classic van der Waals surface of Bader et al., 0.002 a.u. is also common,
  and 0.0004 a.u. (this tool's default) is GaussView's. Lower isovalues give
  larger surfaces further from the nuclei, where $\lvert V\rvert$ is smaller
  (ethylene: minimum −0.0258 a.u. at 0.001, −0.0210 a.u. at 0.0004). Compare
  ESP values only between surfaces of the same isovalue.
- **Where/how used:** $\rho$ is evaluated on a cubic grid of spacing $h$
  (`grid_points`). `grid_box` starts from the atoms ± 1.7 van der Waals
  radii (Multiwfn's choice; Bondi/Mantina radii in `help.AtomicData.vdw_radii`,
  2.0 Å when missing) and pushes any face whose largest density is still
  $\ge \rho_\text{iso}/10$ out by 1 bohr, up to 30 rounds — necessary because
  Multiwfn's 1.7 factor is sized for 0.001 a.u. and is too small for 0.0004.
  Several isovalues share one grid, sized for the smallest.
- **References:** R. F. W. Bader, M. T. Carroll, J. R. Cheeseman and C. Chang,
  "Properties of atoms in molecules: atomic volumes", *J. Am. Chem. Soc.*
  **109**, 7968–7979 (1987); A. Bondi, *J. Phys. Chem.* **68**, 441 (1964);
  M. Mantina et al., *J. Phys. Chem. A* **113**, 5806 (2009).

### 8. Isosurface triangulation (marching tetrahedra)

Each grid cube is split into six tetrahedra sharing the main diagonal; a grid
point is *inside* when $\rho \ge \rho_\text{iso}$. In a tetrahedron with one or
three inside corners the surface is one triangle, with two it is two
triangles. On every cut edge from an inside point $\mathbf{r}_a$ to an outside
point $\mathbf{r}_b$ the bracket is first halved three times using the
**exact** density, then the vertex is placed by linear interpolation:

$$
\mathbf{r}_\text{v} = \mathbf{r}_a + \frac{\rho_a - \rho_\text{iso}}{\rho_a - \rho_b}\,(\mathbf{r}_b - \mathbf{r}_a)
$$

- **Represents:** a closed triangulated approximation of $S(\rho_\text{iso})$.
  After the bisections the vertices lie on the isosurface to
  $\lvert\rho - \rho_\text{iso}\rvert \approx 10^{-6}$ a.u. (ethylene).
- **Where/how used:** `marching_tetrahedra`. Vertices are identified by the
  grid edge they lie on, so neighbouring triangles share vertices and the mesh
  is connected — which the extremum test needs. Triangles are oriented so their
  normals point outward, along $-\nabla\rho$ (central differences, step
  10⁻³ bohr; `density_outward`, `orient_triangles`).
- **References:** A. Doi and A. Koide, "An efficient method of triangulating
  equi-valued surfaces by using tetrahedral cells", *IEICE Trans.* **E74**,
  214–224 (1991); T. Lu and F. Chen, *J. Mol. Graph. Model.* **38**, 314–323
  (2012).

### 9. Vertex merging and elimination; surface topology

Two edge-connected vertices closer than half a grid spacing are merged into
their midpoint,

$$
\lvert \mathbf{r}_i - \mathbf{r}_j \rvert < 0.5\,h
\;\Rightarrow\;
\mathbf{r}_i \leftarrow \tfrac12(\mathbf{r}_i + \mathbf{r}_j),
$$

after which vertices with exactly two neighbours, then those with exactly
three (their hole is re-capped by one triangle), are removed. Each closed
surface component must then satisfy Euler's relation

$$
\chi = V - E + F = 2 .
$$

- **Represents:** Multiwfn's clean-up of the uneven marching-tetrahedra mesh:
  it removes the very short edges and sparse vertices that would otherwise
  create spurious surface extrema. On ethylene it removes 71 % of the
  vertices, leaving a mean edge length of 0.248 bohr at $h = 0.25$ bohr.
- **Where/how used:** `eliminate_redundant_vertices` (sequential, in vertex
  order, exactly as Multiwfn; no minimum spacing is guaranteed afterwards,
  because each merge moves the survivor). `SurfaceMesh.component_euler_characteristics`
  checks every component; a value other than 2 adds a `NOTE:` to the results
  suggesting a smaller grid spacing. The surface area is
  $A = \sum_\triangle \tfrac12\lvert(\mathbf{b}-\mathbf{a})\times(\mathbf{c}-\mathbf{a})\rvert$
  (`SurfaceMesh.area`); vertex normals for shading are the normalised sum of the
  incident triangle normals (`vertex_normals`).
- **References:** Lu and Chen (2012), and Multiwfn `surfana.f90`.

### 10. Surface extrema

Let $\mathcal{N}_1(i)$ be the vertices joined to vertex $i$ by an edge and
$\mathcal{N}_2(i) = \bigl(\bigcup_{j\in\mathcal{N}_1(i)}\mathcal{N}_1(j)\bigr)\setminus\{i\}$
its second ring. With $V_i = V(\mathbf{r}_i)$:

$$
i \text{ is a local minimum} \iff V_i < V_k \;\; \forall\, k \in \mathcal{N}_1(i)\cup\mathcal{N}_2(i),
\qquad
i \text{ is a local maximum} \iff V_i > V_k \;\; \forall\, k \in \mathcal{N}_1(i)\cup\mathcal{N}_2(i)
$$

$$
i_\text{global min} = \arg\min_i V_i, \qquad i_\text{global max} = \arg\max_i V_i
$$

- **Represents:** Multiwfn's definition of surface extrema, copied without
  change. The inequality is **strict**, so equal values (a plateau) give no
  extremum, and requiring the second ring as well as the first suppresses
  shallow ripples one edge wide. Two extrema of the same kind are therefore
  always at least three edges apart (about 0.4 Å at the default spacing,
  where the mean edge is 0.248 bohr). There
  is **no** value tolerance and **no** automatic merging of nearby extrema; the
  position reported is the vertex position, not an interpolated one. The
  global extremum is always reported, as `global`, even when a tie keeps it out
  of the strict local test.
- **Where/how used:** `molecularSurface.find_surface_extrema` builds the
  $\mathcal{N}_1\cup\mathcal{N}_2$ pattern as a sparse matrix (`second_ring`) and
  compares each vertex with the row minimum/maximum. The nearest atom is the
  nucleus closest to the extremum.
- **References:** Lu and Chen (2012); T. Lu and F. Chen, "Multiwfn: A
  multifunctional wavefunction analyzer", *J. Comput. Chem.* **33**, 580–592
  (2012), and Multiwfn `surfana.f90`.

### 11. Colour mapping

$$
R = \max_i \lvert V_i\rvert,
\qquad
f_i = \operatorname{clip}\!\left(\frac{V_i + R}{2R},\,0,\,1\right),
\qquad
\text{hue}_i = 240^\circ \times f_i
$$

- **Represents:** GaussView's default appearance: a symmetric range about
  zero, red ($V=-R$) → yellow → green ($V=0$) → cyan → blue ($V=+R$), at full
  saturation and brightness. Because $R$ comes from the surface shown, the same
  colour means different potentials on different molecules or isovalues; read
  the bar labels.
- **Where/how used:** `espViewer.symmetric_range`, `esp_colours`, `hue_colours`;
  labels by `gaussview_label` (`-2.097e-2` style), converted to kcal/mol on
  request.

### 12. Units and conversion factors

| Quantity | Internal unit | Reported as | Factor |
|---|---|---|---|
| Length | bohr ($a_0$) | Å | $1\,a_0 = 0.529177210544$ Å (CODATA 2022; `coordinateConverter.BOHR_RADIUS_ANGSTROM`) |
| Electron density | e/bohr³ | e/bohr³ (a.u.) | $1$ e/bohr³ $= 6.7483$ e/Å³; 0.0004 a.u. = 0.0027 e/Å³ |
| Potential | Hartree/e | a.u. and kcal/mol | $1\,E_h = 627.5094740631$ kcal/mol (CODATA 2018 $E_h\,N_A/4184$ J kcal⁻¹; `HARTREE_TO_KCAL_MOL`) |
| Area | bohr² | Å² | $a_0^2 = 0.280028520$ Å² |

Multiwfn prints kcal/mol with the rounded factor 627.51 (relative difference
8×10⁻⁷, invisible at the printed precision). 1 a.u. of potential is also
27.211 V, which is not reported.

## Formulas and Methods (Fukui Function tool)

This section documents the theory and every algorithm behind the
[Analysis of Fukui Functions](#analysis-of-fukui-functions) tool, in the order
the tool evaluates them. The code is split by responsibility:
- [fukuiFunctionAnalyze.py](fukuiFunctionAnalyze.py): the three-state checks, the formulas, the grids and surfaces, the report, and the window;
- [gaussianWavefunction.py](gaussianWavefunction.py): the `.fchk` reader and the electron density;
- [densityNumba.py](densityNumba.py): the optional compiled kernels;
- [hirshfeldCharges.py](hirshfeldCharges.py) with [proatomDensities.py](proatomDensities.py): the Hirshfeld charges;
- [molecularSurface.py](molecularSurface.py): grids, isosurfaces and extrema, shared with the ESP tool;
- [cubeFile.py](cubeFile.py): cube output;
- [espViewer.py](espViewer.py): the 3D view.

Everything except the window and the viewer runs headless. It is tested by:
- `tests/test_fukuiFunctionAnalyze.py`, `test_hirshfeldCharges.py`, `test_gaussianWavefunction.py` and `test_cubeFile.py`;
- against reference values produced by Multiwfn 2026.7.15 and Gaussian 09.

No external program is needed at run time.

**Status of the results.** Fukui functions and the quantities derived from
them indicate *where* a molecule is most susceptible to a given kind of
attack. They do not predict reaction outcomes, rates or selectivities: steric
hindrance, solvent, distortion of the reactants and transition-state energies
are outside the model (subsection 13).

**Units.** Hartree atomic units throughout: lengths in bohr, energies in
Hartree, densities and Fukui functions in e/bohr³ (written "a.u."), charges in
units of the elementary charge e. Å and eV appear only in the report
($1\ \text{bohr} = 0.529177210544$ Å, $1\ E_h = 27.211386245988$ eV, CODATA).

**Symbols used throughout.**

| Symbol | Meaning |
| --- | --- |
| $N$ | number of electrons of the reference (N) state |
| $v(\mathbf{r})$ | external potential (the nuclei, at fixed positions) |
| $\rho_N(\mathbf{r})$ | electron density of the state with $N$ electrons (fchk `Total SCF Density`, α + β) |
| $E(N)$ | total energy of the state with $N$ electrons (fchk `Total Energy`) |
| $f^+, f^-, f^0$ | Fukui functions for nucleophilic, electrophilic and radical attack |
| $\Delta f$ | dual descriptor (DD, or CDD when condensed) |
| $p, q$ | electrons added (N+p state) and removed (N−q state); 1 normally |
| $q_A,\ Z_A,\ \mathbf{R}_A$ | charge, atomic number and position of atom $A$ |
| $\rho_A^0(r)$ | spherically averaged density of the free neutral atom $A$ |
| $w_A(\mathbf{r})$ | Hirshfeld weight of atom $A$ at $\mathbf{r}$ |
| $I,\ A$ | vertical ionization potential and electron affinity |
| $\mu,\ \eta,\ S,\ \omega$ | chemical potential, hardness, softness, electrophilicity |
| $h$ | spacing of the uniform grid (bohr) |

### 1. Conceptual DFT: response of the density to the electron number

Conceptual density functional theory treats the energy as a functional of
the electron number and the external potential, $E = E[N, v]$. Its first
derivatives are the chemical potential and the density:

$$
\mu = \left(\frac{\partial E}{\partial N}\right)_{v},
\qquad
\rho(\mathbf{r}) = \left(\frac{\delta E}{\delta v(\mathbf{r})}\right)_{N}.
$$

The mixed second derivative is the **Fukui function** (Parr and Yang, 1984).
By the Maxwell relation it is both the response of the density to a change in
$N$ at fixed nuclei, and the response of $\mu$ to a local change in the
potential:

$$
f(\mathbf{r})
= \left(\frac{\partial \rho(\mathbf{r})}{\partial N}\right)_{v}
= \left(\frac{\delta \mu}{\delta v(\mathbf{r})}\right)_{N},
\qquad
\int f(\mathbf{r})\, d^3r = 1 .
$$

The normalization follows from $\int \rho\, d^3r = N$.

A region where $f$ is large is where electron density is gained or lost most
readily when the molecule accepts or donates charge. In the frontier-orbital
picture, it is where an incoming reagent can interact most strongly.
"Fixed $v$" means fixed nuclei, which is why all three states must share one
geometry (subsection 4).

### 2. Finite differences: f⁺, f⁻ and f⁰

$E(N)$ and $\rho_N$ have a derivative discontinuity at integer $N$, so the
derivative is replaced by one-sided finite differences between states that
differ by whole electrons:

$$
f^{+}(\mathbf{r}) = \rho_{N+1}(\mathbf{r}) - \rho_{N}(\mathbf{r}),
\qquad
f^{-}(\mathbf{r}) = \rho_{N}(\mathbf{r}) - \rho_{N-1}(\mathbf{r}),
\qquad
f^{0}(\mathbf{r}) = \tfrac12\left[f^{+}(\mathbf{r}) + f^{-}(\mathbf{r})\right]
= \tfrac12\left[\rho_{N+1}(\mathbf{r}) - \rho_{N-1}(\mathbf{r})\right].
$$

| Function | Change of the molecule | Large where the site is susceptible to |
| --- | --- | --- |
| $f^{+}$ | accepts an electron | **nucleophilic** attack (an electron donor attacks here) |
| $f^{-}$ | donates an electron | **electrophilic** attack (an electron acceptor attacks here) |
| $f^{0}$ | average of both | **radical** attack (little net charge transfer) |

The names describe the **attacking reagent**. A site "susceptible to
electrophilic attack" is a nucleophilic site.

**Frozen orbitals.** If the orbitals did not relax when an electron is added
or removed, $f^+ = |\varphi_\text{LUMO}|^2$ and $f^- = |\varphi_\text{HOMO}|^2$
(frontier-orbital theory). The tool always uses the relaxed densities of
separate calculations, which include orbital relaxation.

**Reliability.** Comparison studies by the Multiwfn author find $f^-$
reasonably reliable for electrophilic-attack sites, and $f^+$ the weakest
descriptor for nucleophilic-attack sites. The dual descriptor (subsection 3)
and Hirshfeld charges do better there.

### 3. The dual descriptor

The second derivative of the density with respect to $N$, also the response
of the hardness to the potential, is the **dual descriptor** (Morell, Grand and
Toro-Labbé, 2005). Its finite-difference form is

$$
\Delta f(\mathbf{r})
= \left(\frac{\partial f(\mathbf{r})}{\partial N}\right)_{v}
\;\approx\; f^{+}(\mathbf{r}) - f^{-}(\mathbf{r})
= \rho_{N+1}(\mathbf{r}) - 2\rho_{N}(\mathbf{r}) + \rho_{N-1}(\mathbf{r}),
\qquad \int \Delta f\, d^3r = 0 .
$$

$\Delta f > 0$ marks regions favouring **nucleophilic** attack and
$\Delta f < 0$ regions favouring **electrophilic** attack. One function thus
shows both kinds of site. When ranking competing sites, compare how positive
or how negative $\Delta f$ is, not only its sign.

### 4. The three states and how they are checked

The tool reads three Gaussian formatted checkpoints. For each it takes:
- the atoms and coordinates;
- the basis set, which gives the density through equation 2 of the ESP section;
- `Total SCF Density`, `Total Energy` and the electron counts;
- the alpha orbital energies;
- the method and basis names on line 2 of the fchk.

The finite difference is only meaningful when the three files differ in
nothing but the electron count. A **ValueError**, shown as a dialog, is raised
when:

| Check | Tolerance / rule |
| --- | --- |
| same atoms, in the same order | atomic numbers equal |
| same geometry | every coordinate within $10^{-5}$ bohr; message: run single points, do not re-optimize the ions |
| same basis set | shells, angular momenta, centres and exponents equal (relative $10^{-10}$) |
| same method | method names equal after removing the R / U / RO prefix (`RB3LYP` = `UB3LYP`) |
| electron counts | $N_{+} > N > N_{-}$; otherwise the files are probably swapped |

With $p = N_{+} - N$ and $q = N - N_{-}$ (subsection 5), a **NOTE**, which is not an
error, is written when:
- the multiplicity changes by anything other than $p$ (or $q$);
- the file also contains a correlated density (e.g. `Total MP2 Density`). The SCF density is used.
- the basis set has diffuse functions (subsection 13);
- the N state has degenerate frontier orbitals that the files do not account for (subsection 5).

**Charge and multiplicity of the ions** for ordinary organic molecules:

| Reference state | N | N+1 | N−1 |
| --- | --- | --- | --- |
| closed-shell neutral | 0 1 | −1 2 | +1 2 |
| neutral radical | 0 2 | −1 1 | +1 1 |
| closed-shell cation | +1 1 | 0 2 | +2 2 |
| closed-shell anion | −1 1 | −2 2 | 0 2 |

For transition-metal systems, compute the lowest few spin states of each ion
and use the lowest in energy.

### 5. Degenerate frontier orbitals: N+p and N−q states

If the LUMO is $p$-fold and the HOMO $q$-fold degenerate, as in benzene or
C₆₀, adding or removing one electron breaks the molecular symmetry, and $f^\pm$
come out misleading. The remedy (Multiwfn manual, section 3.25.4) spreads the
change over the whole degenerate set:

$$
f^{+}_{Q} = \frac{\rho_{N+p} - \rho_{N}}{p},
\qquad
f^{-}_{Q} = \frac{\rho_{N} - \rho_{N-q}}{q},
\qquad
f^{0}_{Q} = \tfrac12\left(f^{+}_{Q} + f^{-}_{Q}\right),
\qquad
\Delta f_{Q} = f^{+}_{Q} - f^{-}_{Q}.
$$

The multiplicities are $p+1$ and $q+1$ for a closed-shell N state. The tool
infers $p$ and $q$ from the electron counts of the files, so no option needs
setting, and it applies the same division to the condensed indices
(subsection 7).

For a closed-shell N state it also measures the degeneracy itself, by
Multiwfn's rule: the number of occupied (virtual) orbitals within
**0.01 eV** of the HOMO (LUMO), counting the HOMO (LUMO) itself. If this
disagrees with the $p, q$ of the files, a NOTE suggests N+p/N−q states or
Multiwfn's orbital-weighted Fukui functions.

The global indices of subsection 6 need $E(N\pm1)$, so they are computed only
when $p = q = 1$.

### 6. Global indices

Fit $E$ by a parabola through $E(N-1)$, $E(N)$ and $E(N+1)$, with the vertical
energies $I = E(N-1) - E(N)$ and $A = E(N) - E(N+1)$. Its first derivative is
the chemical potential and its second derivative the hardness:

$$
\mu \approx -\frac{I + A}{2} = -\chi,
\qquad
\eta = \frac{\partial^2 E}{\partial N^2} \approx I - A,
\qquad
S = \frac{1}{\eta},
\qquad
\omega = \frac{\mu^2}{2\eta},
\qquad
N_\text{Nu} = E_\text{HOMO}(N) - E_\text{HOMO}(\text{TCE}).
$$

- $\chi$ is Mulliken's electronegativity.
- **Hardness convention.** Many papers write $\eta = (I-A)/2$. Like Multiwfn, this tool uses $\eta = I - A$, the finite-difference second derivative itself (the fundamental gap). $S$ and $\omega$ follow from it, so compare only with values computed by the same convention.
- $\omega$ is Parr's electrophilicity index.
- $N_\text{Nu}$ is the nucleophilicity index of Domingo and co-workers, referenced to tetracyanoethylene. The tool uses Multiwfn's $E_\text{HOMO}(\text{TCE}) = -0.335198\ E_h$, computed at B3LYP/6-31G\*. At any other level, compute TCE yourself and subtract.
- $N_\text{Nu}$ is reported only for a closed-shell N state, where the HOMO is unambiguous.
- If $I - A \le 0$, $S$ and $\omega$ are left undefined and a NOTE asks you to check the charges and multiplicities.

**Local versions.** The Fukui functions compare sites **within one molecule**:
$\int f = 1$, so $f$ thins out as molecules grow. To compare sites of
different molecules, scale $f$ by a global index. The local softness is
$s^{\pm,0}(\mathbf{r}) = S\, f^{\pm,0}(\mathbf{r})$; the report gives its
condensed form, $s_A = S\, f_A$.

### 7. Atom-condensed indices

Integrating a Fukui function over the region that belongs to atom $A$ turns it
into a difference of atomic populations $p_A = Z_A - q_A$ (Yang and Mortier,
1986). The nuclear charges cancel, leaving charge differences:

$$
f^{+}_{A} = \frac{q_A(N) - q_A(N+p)}{p},
\qquad
f^{-}_{A} = \frac{q_A(N-q) - q_A(N)}{q},
\qquad
f^{0}_{A} = \tfrac12\left(f^{+}_{A} + f^{-}_{A}\right),
\qquad
\text{CDD}_A = f^{+}_{A} - f^{-}_{A}.
$$

The values depend on how space is divided among atoms. **Hirshfeld charges**
are the recommended partition for this purpose: in a comparison of 15 charge
schemes they predicted regioselectivity best (Wang, Rong, Chattaraj and Liu,
2019), and Multiwfn uses them in its conceptual-DFT module. Mulliken charges
are discouraged. Hirshfeld condensed values are usually non-negative;
negative values can still occur (orbital relaxation, the partition). The tool
lists them in a NOTE, because they should not be read as reactivity.

The **site ranking** in the report lists, for each kind of attack, the five
atoms with the largest $f^-$, $f^+$ or $f^0$, and the most negative and most
positive CDD. Consider only atoms that can actually react: the ipso carbon of
phenol, for example, scores high but does not take part in substitution.

### 8. Hirshfeld (stockholder) partition

Hirshfeld (1977) shares the molecular density among the atoms in proportion to
the density each free atom would have at that point. The free atoms are placed
at the molecular geometry; their sum is the "promolecule":

$$
w_A(\mathbf{r}) = \frac{\rho_A^0\!\left(\lvert \mathbf{r}-\mathbf{R}_A\rvert\right)}
                       {\sum_B \rho_B^0\!\left(\lvert \mathbf{r}-\mathbf{R}_B\rvert\right)},
\qquad
\sum_A w_A(\mathbf{r}) = 1,
\qquad
p_A = \int w_A(\mathbf{r})\, \rho(\mathbf{r})\, d^3r,
\qquad
q_A = Z_A - p_A .
$$

**Free-atom densities** ([proatomDensities.py](proatomDensities.py)) are a
generated data module that must not be edited by hand.
- **Level.** Neutral ground-state atoms H–Ar at UB3LYP/cc-pVQZ: the level of Multiwfn's built-in densities for $Z \le 18$ (Multiwfn manual, Appendix 3). Ground-state multiplicities are H, Li, B, F, Na, Al, Cl 2; C, O, Si, S 3; N, P 4; the noble gases and Be, Mg 1.
- **Sphericalization.** Open-shell atoms, such as C with two p electrons, have non-spherical densities, so each density is averaged over a 110-point Lebedev sphere (exact for the angular dependence of an s–g density):
  $\rho_A^0(r) = \sum_k \omega_k\, \rho_A(r\,\hat{\mathbf{n}}_k)$, with $\sum_k \omega_k = 1$.
- **Table.** 600 radii, log-spaced from $10^{-5}$ to 40 bohr, evaluated without the shell screening of subsection 11, so that the tails stay positive.
- **Interpolation.** $\ln\rho_A^0$ is interpolated linearly in $r$. Inside $10^{-5}$ bohr the first value is used; beyond 40 bohr the density is zero.
- **Normalization check.** Each table integrates to $Z_A$ within $3\times10^{-4}$ relative (`test_each_table_holds_its_electrons`).
- **Rebuilding.** `tools/make_proatom_densities.py` regenerates the tables; it needs Gaussian 09.

**Final normalization.** The integrated populations are scaled so that they
sum exactly to the wavefunction's electron count:
$p_A \leftarrow p_A\, N / \sum_B p_B$. Then $\sum_A q_A$ equals the molecular
charge to rounding, as in Multiwfn's "Final atomic charges".

**Restrictions** (ValueError; the tool turns it into a NOTE and keeps all
real-space results):
- **Elements beyond Ar.** They have no table.
- **Files with effective core potentials** (nuclear charges ≠ atomic numbers). Their density lacks the core electrons, so comparing it with all-electron free atoms would give wrong charges.

### 9. Becke multicentre integration

The integral of equation (8) runs over a sum of atom-centred grids. A smooth
partition of unity removes the cusps and overlaps between them (Becke, 1988).

**Fuzzy cells.** For atoms $A \ne B$ at distance $R_{AB}$, and a point at
distances $r_A$, $r_B$ from them,

$$
\mu_{AB} = \frac{r_A - r_B}{R_{AB}},
\qquad
g(\mu) = \tfrac32\mu - \tfrac12\mu^3,
\qquad
s(\mu) = \tfrac12\left[1 - g\big(g(g(\mu))\big)\right],
$$

$$
P_A(\mathbf{r}) = \prod_{B \ne A} s(\mu_{AB}),
\qquad
W_A(\mathbf{r}) = \frac{P_A(\mathbf{r})}{\sum_C P_C(\mathbf{r})}.
$$

The three-fold iteration of $g$ makes $s$ a smooth step: 1 near atom $A$, 0
near atom $B$. Becke's atomic-size adjustment is not applied; the
normalization of subsection 8 absorbs the small integration error.

**Atom-centred quadrature.** Around each atom, the radial coordinate uses
Gauss–Chebyshev quadrature of the second kind, mapped to $[0, \infty)$:

$$
x_i = \cos\theta_i,\ \ \theta_i = \frac{i\pi}{n+1},\ \ i = 1..n,
\qquad
r_i = r_m \frac{1 + x_i}{1 - x_i},
\qquad
w_i = \frac{\pi}{n+1}\,\sin\theta_i\;\frac{2 r_m}{(1 - x_i)^2}\; r_i^{2}.
$$

$r_m$ is the atom's covalent radius (Cordero et al., 2008). The angular part
is a 434-point Lebedev sphere, exact for spherical harmonics up to degree 35.
With $n = 75$ these are Multiwfn's default 75 × 434 = 32 550 points per atom:

$$
\int F(\mathbf{r})\, d^3r \;\approx\; \sum_{A}\sum_{i,k}
w_i\, \omega_k\, W_A(\mathbf{r}_{A,ik})\, F(\mathbf{r}_{A,ik}),
\qquad
\mathbf{r}_{A,ik} = \mathbf{R}_A + r_i \hat{\mathbf{n}}_k .
$$

For phenol (13 atoms, 423 150 points) the grid integrates $\rho_N$ to the
electron count within $10^{-5}$. **One grid serves all three states**, because
they share the geometry.

### 10. The density on grids and surfaces

The density at any point is the exact expression of the ESP section,
equation 2: $\rho(\mathbf{r}) = \sum_{\mu\nu} P_{\mu\nu}\chi_\mu(\mathbf{r})\chi_\nu(\mathbf{r})$.
Every value used by the tool comes from this formula, not from interpolating
a grid: the grid nodes, the isosurface refinement, the surface vertices and
the Hirshfeld points.

**Uniform grid (for the isosurfaces and the cube files).** One cubic grid of
spacing $h$ (default 0.25 bohr) holds $\rho_N$, $\rho_{N+p}$ and $\rho_{N-q}$,
and from them $f^+$, $f^-$, $f^0$ and $\Delta f$.
- **Starting box.** The atoms ± 1.7 van der Waals radii, as in the ESP tool.
- **Growth.** Each face moves out by 1 bohr until it satisfies, at every face node,
  $$
  \max\!\left(\frac{\rho_N}{\rho_\text{s}},\ \frac{|f^+|}{f_\text{iso}},\ \frac{|f^-|}{f_\text{iso}}\right) < 0.1 .
  $$
  Here $\rho_\text{s}$ is the mapping-surface isovalue (default 0.01 a.u.) and $f_\text{iso}$ the Fukui isovalue (default 0.005 a.u.). Since $|f^0| \le \max|f^\pm|$ and $|\Delta f| \le |f^+| + |f^-|$, no isosurface of any of the four functions is cut by the box.

**Diagnostic integrals.** The report prints $h^3 \sum f$ over the grid. It
should be close to 1 for $f^{\pm,0}$ (about 0.999 for phenol) and to 0 for
$\Delta f$. A value near −1 means the N+1 and N−1 files were swapped.

**Isosurfaces ("lobes").** For each function and each sign $s = \pm1$, the
surface $s f(\mathbf{r}) = f_\text{iso}$ is built by the ESP tool's algorithm:
- marching tetrahedra on the grid of $s f$, with vertices shared through edge identities;
- three bisections of the **exact** $s f$ along each cut edge;
- vertex merging (closer than $0.5h$) and elimination;
- triangles oriented along $-\nabla(s f)$.

The positive and negative surfaces are built as separate meshes and drawn
together, or one at a time with *Positive/Negative isosurface only*.

**Mapped surface and extrema.** The four functions are also evaluated exactly
at the vertices of the isosurface $\rho_N = \rho_\text{s}$. The default
0.01 a.u. follows Multiwfn's manual, section 4.12.4: on the usual 0.001 a.u.
van der Waals surface the Fukui functions are too small to discriminate.
Minima and maxima on that surface use the ESP tool's rule: a vertex is a
minimum (maximum) if it is strictly below (above) all its first- and
second-ring neighbours, and the global extremum is always reported. Each
extremum is labelled with its nearest atom.

**Cube files.** `<name>_fplus.cub`, `_fminus.cub`, `_fzero.cub` and
`_dual.cub` contain the four functions in the standard Gaussian cube layout:
bohr, the atoms' numbers and positions, values `%13.5E`, six per line, z
fastest. VMD, Multiwfn and GaussView read them.

The cubes have their own grid, spacing $h_\text{cube}$ (default 0.1 bohr). It
starts at the analysis grid's origin and covers at least its box, and every
value is computed exactly at its grid point, never interpolated from the
analysis grid. Because $\rho = \chi^\mathsf{T} P\, \chi$ is linear in the
density matrix,

$$
f^+(\mathbf{r}) = \chi(\mathbf{r})^\mathsf{T}\,\frac{P_{N+p} - P_N}{p}\,\chi(\mathbf{r}),
\qquad
f^-(\mathbf{r}) = \chi(\mathbf{r})^\mathsf{T}\,\frac{P_N - P_{N-q}}{q}\,\chi(\mathbf{r}),
$$

so the basis functions are evaluated once per point and contracted with two
difference matrices instead of three densities. It is also more accurate:
near a nucleus the three densities are of order 10²–10³ a.u. while $f$ is
of order 10⁻², and subtracting them loses digits. For acetone the two routes
differ by at most $9\times10^{-14}$ a.u. A cube grid at the analysis spacing
writes the analysis grid itself. Grids beyond 25 million points (about
330 MB per file) are refused with the smallest spacing that fits.

### 11. Compiled kernels (Numba)

The grid work is dominated by evaluating the Cartesian basis functions

$$
\chi_j(\mathbf{r}) = N_j\, x_c^{a}\, y_c^{b}\, z_c^{c}
\sum_k d_k\, N(\alpha_k)\, e^{-\alpha_k r_c^2}
$$

at many points ($x_c, y_c, z_c$ relative to the shell centre). With `numba`
installed, [densityNumba.py](densityNumba.py) evaluates them in one compiled,
multi-threaded loop over points.
- **Screening.** A shell is set to zero at a point when even its most diffuse primitive satisfies $\alpha_\text{min} r_c^2 > 50$, i.e. it is below $e^{-50} \approx 2\times10^{-22}$ of its coefficient.
- **Contraction.** The step $\rho = \sum_{\mu\nu}P_{\mu\nu}\chi_\mu\chi_\nu$ stays a matrix product (BLAS) on blocks of 20 000 points.
- **Becke weights.** The same module computes the Becke weights of subsection 9.
- **Cube text.** Writing the values as `%13.5E` text in Python takes about 0.6 µs per value. [volumeNumba.py](volumeNumba.py) writes them in a compiled, parallel loop and produces the same bytes. A value whose sixth significant digit lies within $10^{-7}$ of a rounding tie, a non-finite value or an exponent beyond ±300 is left to Python, about 2 in 10⁷ values in practice. Measured: 3.75 M values, 2.3 s in Python vs 0.11 s compiled.

Without `numba`, the numpy path runs and gives the same numbers:

| Quantity | Agreement, compiled vs numpy |
| --- | --- |
| density, every basis type (pure/Cartesian d, f, g; SP shells; open shell) | ≤ $10^{-12}$ of the largest density |
| Becke weights | ≤ $10^{-12}$ of the largest weight |
| Hirshfeld populations | ≤ $10^{-9}$ e |

The Becke weights cannot agree to the last digit where they are tiny. Near a
cell boundary $\mu \to 1$, and the step $\tfrac12(1 - \mu)$ amplifies
differences in rounding order; those points carry negligible weight.

**Timing** (phenol, 0.25 bohr grid, 222 300 points):
- density: 0.83 s with numpy, 0.20 s compiled;
- whole run, including the condensed indices: about 8 s.

**Timing of the cube export** (acetone, 6-31++G(d,p), 0.1 bohr, 3.87 M points):
- f⁺ and f⁻ on the cube grid: 8.4 s as three densities, 3.9 s from the difference matrices;
- four cube files: 9.2 s of Python formatting, under 1 s compiled;
- whole export: about 6 s, against about 18 s by the old route.

The first run after installing `numba` spends a few seconds compiling; later
runs load the cache.

### 12. The 3D view

Surfaces use the ESP tool's colour scale, red → yellow → green (0) →
cyan → blue.
- **Isosurfaces.** The bar spans $\pm f_\text{iso}$, so negative lobes are red and positive lobes blue.
- **Mapped surface.** The bar spans $\pm\max|f|$ on that surface. A function that is positive everywhere, like $f^-$ on phenol, therefore uses only the green-to-blue half.

Extremum markers are cyan (minimum) and orange (maximum), larger for the
global ones, and numbered as in the report.

**Atoms row** (Results tab), applied at once with or without results:

| Choice | Atom radius (Å) | Bonds |
| --- | --- | --- |
| Line | 0.18 × visual radius, at least 0.06 | lines |
| Ball and stick (CPK) | visual radius (e.g. H 0.25, C 0.42, O 0.40) | cylinders |
| Space-filling (vdW) | van der Waals radius (e.g. H 1.20, C 1.70, O 1.52) | none needed |

"CPK" means ball-and-stick here, as in the 3D Molecular Viewer. **Show atom
numbers** labels every atom with its 1-based index, the numbering of the
report and the condensed table.

The surfaces are opaque. The Fukui lobes lie inside the van der Waals radii,
so space-filling atoms hide them; use that style to judge how exposed a site
is, then switch back. With atom numbers on, some labels may sit behind a lobe.

### 13. Limitations

- **Level of theory.** Removing or adding a whole electron is a large perturbation, and the result depends on how well the method describes the ions.
  - B3LYP-like functionals over-delocalize the extra electron of the anion (self-interaction error).
  - With **diffuse basis functions** the N+1 density can become Rydberg-like, so $f^+$ and $\Delta f$ lose chemical meaning (Multiwfn manual, section 4.22.3.1). The tool writes a NOTE when the basis set has diffuse functions.
  - Diffuse functions are, however, what makes $E(N+1)$, and so $A$, $\chi$ and $\omega$, accurate.
  - Multiwfn's manual recommends a range-separated functional (ωB97XD) with a triple-ζ basis without diffuse functions for $f^+$ and $\Delta f$, and calls B3LYP/6-31G\* the lowest acceptable level.
- **Unbound anions.** A negative vertical EA (phenol: −2.21 eV at B3LYP/6-31G\*) means the extra electron is unbound at that level. Read $f^+$ and $\Delta f$ with care.
- **Density.** Only the SCF (HF/DFT) density is used. Open-shell states contribute their total (α + β) density; spin-resolved Fukui functions are not computed.
- **Partition.** Condensed values depend on the partition (subsection 7). The real-space functions and the surface extrema do not.
- **Elements.** Hirshfeld charges cover H–Ar. Heavier elements and ECP files give a NOTE, and the real-space functions are still produced.
- **Not provided** (use Multiwfn main function 22 or 12):
  - Fukui potential and dual-descriptor potential;
  - orbital-weighted Fukui functions;
  - bond dual descriptor;
  - $\omega_\text{cubic}$ and the electrophilic descriptor;
  - averages over local atomic surfaces;
  - orbital contributions to $f$.
- **What Fukui functions cannot say.** Sterics, distortion of the reactants, non-covalent steering, solvent, concerted mechanisms and kinetics. When descriptors disagree, or when the answer matters, compare transition-state free-energy barriers.

### 14. Validation

| Quantity | Reference | Agreement | Test |
| --- | --- | --- | --- |
| $f^+, f^-, f^0, \Delta f$ at 300 grid nodes (phenol) | Multiwfn 2026.7.15 cubes, main function 22 option 3 | ≤ $2\times10^{-5}$ a.u. (measured $8\times10^{-6}$) | `test_grid_matches_multiwfn` |
| Global indices (phenol) | Multiwfn `CDFT.txt` | every printed digit | `test_global_indices_reproduce_multiwfn` |
| Hirshfeld charges, 8 molecules (H, Li, B–F, Na–Cl) | Multiwfn main function 7, option 1 | ≤ $3\times10^{-4}$ e (measured $1.4\times10^{-4}$) | `TestAgainstMultiwfn` |
| Hirshfeld charges, phenol N, N+1, N−1 | Multiwfn `CDFT.txt` | ≤ $2\times10^{-4}$ e | `test_phenol_three_states` |
| Condensed $q$, $f^\pm$, $f^0$, CDD (phenol) | Multiwfn `CDFT.txt` | ≤ $3\times10^{-4}$ (4 printed decimals) | `test_condensed_indices_reproduce_multiwfn` |
| Density, compiled vs numpy | the numpy path | ≤ $10^{-12}$ relative | `TestNumbaDensity` |

Multiwfn is not needed by the tool or the tests. The two dev-only scripts
`tests/data/fukui/make_reference.py` and `tests/data/hirshfeld/make_reference.py`
regenerate the reference files.

### 15. References

- R. G. Parr, W. Yang, *J. Am. Chem. Soc.* **106**, 4049 (1984): the Fukui function.
- W. Yang, W. J. Mortier, *J. Am. Chem. Soc.* **108**, 5708 (1986): condensed Fukui functions.
- C. Morell, A. Grand, A. Toro-Labbé, *J. Phys. Chem. A* **109**, 205 (2005): the dual descriptor.
- R. G. Parr, R. G. Pearson, *J. Am. Chem. Soc.* **105**, 7512 (1983): absolute hardness.
- W. Yang, R. G. Parr, *Proc. Natl. Acad. Sci. USA* **82**, 6723 (1985): softness and local softness.
- R. G. Parr, L. v. Szentpály, S. Liu, *J. Am. Chem. Soc.* **121**, 1922 (1999): the electrophilicity index.
- L. R. Domingo et al., *J. Org. Chem.* **73**, 4615 (2008): the nucleophilicity index.
- F. L. Hirshfeld, *Theor. Chim. Acta* **44**, 129 (1977): stockholder partitioning.
- A. D. Becke, *J. Chem. Phys.* **88**, 2547 (1988): multicentre integration.
- V. I. Lebedev, D. N. Laikov, *Dokl. Math.* **59**, 477 (1999): Lebedev quadrature.
- B. Cordero et al., *Dalton Trans.* 2832 (2008): covalent radii.
- B. Wang, C. Rong, P. K. Chattaraj, S. Liu, *Theor. Chem. Acc.* **138**, 124 (2019): Hirshfeld charges for regioselectivity.
- P. Fuentealba, P. Pérez, R. Contreras, *J. Chem. Phys.* **113**, 2544 (2000): condensed Fukui functions and their sign.
- P. Geerlings, F. De Proft, W. Langenaeker, *Chem. Rev.* **103**, 1793 (2003): conceptual DFT review.
- T. Lu, F. Chen, *J. Comput. Chem.* **33**, 580 (2012); T. Lu, *J. Chem. Phys.* **161**, 082503 (2024): Multiwfn.
- T. Lu, Q. Chen, in *Conceptual Density Functional Theory*, Wiley-VCH (2022), pp. 631–647: conceptual DFT in Multiwfn.
- T. Lu, *Multiwfn Manual* 2026.7.10: sections 3.25, 4.5.4, 4.7.3, 4.12.4, 4.22 and Appendix 3.

The full methodological discussion, with the verified Multiwfn workflow, is in
`docs/superpowers/specs/2026-10-16-fukui-function-design.md`.

## Theoretical Background: the Quantum Theory of Atoms in Molecules

This section explains the ideas behind the [QTAIM tab](#qtaim-topology-analysis). The exact formulas and algorithms the program uses are in [Formulas and Methods (QTAIM tab)](#formulas-and-methods-qtaim-tab). Atomic units are used throughout (ρ in e/bohr³, energy densities in Hartree/bohr³).

### 1. The idea

The *Quantum Theory of Atoms in Molecules* (QTAIM), developed by R. F. W. Bader and co-workers from the 1970s, defines atoms, bonds, rings and cages from a single physical quantity: the electron density

$$
\rho(\mathbf{r}) = N \int |\Psi(\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_N)|^2 \, ds_1\, d\mathbf{x}_2 \cdots d\mathbf{x}_N .
$$

ρ(**r**) is an observable (it is measured by X-ray diffraction) and does not depend on how the wavefunction is written — orbitals, basis functions or localization schemes. QTAIM reads the chemical structure from the *topology* of ρ: its maxima, minima and saddle points, and the lines of steepest ascent that connect them.

### 2. The gradient field and atomic basins

At every point, the gradient ∇ρ points in the direction in which the density increases fastest. Following ∇ρ from a point traces a **gradient path**. Almost every gradient path ends at a nucleus, because ρ has a sharp maximum at each nucleus: the nuclei are the **attractors** of the gradient field.

The set of points whose gradient paths end at the same nucleus is that nucleus's **atomic basin** Ω. A QTAIM *atom* is a nucleus together with its basin. Two neighbouring basins meet at an **interatomic surface**, which no gradient path crosses:

$$
\nabla\rho(\mathbf{r}) \cdot \mathbf{n}(\mathbf{r}) = 0 \quad \text{for every point } \mathbf{r} \text{ of the surface}
$$

(**n** is the unit normal of the surface). This *zero-flux condition* is what makes the QTAIM atom a proper quantum subsystem: within a zero-flux basin, quantum mechanics (Schwinger's principle of stationary action) defines the atom's energy and other properties, and they add up to the molecular values. Because the flux of ∇ρ through the surface of a zero-flux basin is zero,

$$
\int_\Omega \nabla^2\rho \, d\mathbf{r} = \oint_{S(\Omega)} \nabla\rho \cdot \mathbf{n} \, dS = 0 ,
$$

which is why kinetic energy densities that differ by a multiple of ∇²ρ (section 6) give the same atomic kinetic energy, and why the virial theorem holds atom by atom. Integrating ρ over a basin gives the atom's electron population N(Ω) and its charge q(Ω) = Z − N(Ω). *(The QTAIM tab finds the topology; it does not integrate the basins. Atomic charges and volumes are Multiwfn's basin analysis, main function 17.)*

### 3. Critical points

A **critical point** (CP) is a point where ∇ρ = 0. Its character comes from the **Hessian** of ρ, the 3 × 3 matrix of second derivatives ∂²ρ/∂x<sub>i</sub>∂x<sub>j</sub>. Its eigenvalues λ₁ ≤ λ₂ ≤ λ₃ are the curvatures of ρ along three perpendicular directions (the eigenvectors). A CP is labelled (ω, σ): the **rank** ω is the number of non-zero eigenvalues and the **signature** σ is the sum of their signs. In molecules practically all CPs have rank 3:

| Type | Curvatures | What it is | Name |
|---|---|---|---|
| (3,−3) | three negative | a maximum of ρ | nuclear CP (NCP): a nucleus, or rarely a non-nuclear attractor (NNA) |
| (3,−1) | two negative, one positive | a saddle: a maximum in a plane, a minimum along the third direction | bond CP (BCP) |
| (3,+1) | one negative, two positive | a saddle: a minimum in a plane, a maximum along the third direction | ring CP (RCP) |
| (3,+3) | three positive | a minimum of ρ | cage CP (CCP) |

- At a **BCP**, λ₃ > 0 points along the bond: ρ is a *minimum* along the bond path. λ₁ and λ₂ < 0 span the interatomic surface: ρ is a *maximum* within it.
- An **RCP** lies inside a ring of bonded atoms, at the lowest density of the ring surface; a **CCP** lies inside a cage enclosed by rings.
- **Non-nuclear attractors** are maxima away from nuclei. They occur in some metal clusters and in Li₂, and with poor basis sets.
- A CP with a zero eigenvalue (rank < 3) is **degenerate** and unstable: a tiny change of geometry turns it into ordinary CPs or makes it vanish.

The numbers of CPs of an isolated molecule obey the **Poincaré–Hopf relation**

$$
n_\text{NCP} - n_\text{BCP} + n_\text{RCP} - n_\text{CCP} = 1 .
$$

It is a check of completeness: if it fails, CPs are missing. If it holds, CPs may still be missing in pairs (one BCP and one RCP, for example), so the relation is necessary, not sufficient.

### 4. Bond paths and the molecular graph

From a BCP, two gradient paths climb, one to each of two nuclei (they start along the eigenvector of λ₃). Together they form the **bond path**, the line of maximum density linking the two nuclei. The set of all bond paths is the **molecular graph**, the QTAIM picture of the structure: in ordinary molecules it reproduces the bonds of the Lewis structure.

- A bond path can be **curved**. Its length is then longer than the internuclear distance; strained rings give "banana bonds" (in bicyclo[1.1.1]pentane the C–C paths are 0.27 % longer than the C–C distances, [Example 9](#example-9-qtaim-analysis-of-a-hydrogen-bond-water-dimer)).
- A BCP is generally not at the midpoint: it is displaced towards the less electronegative atom (or the smaller basin). In O–H bonds it lies much closer to H.
- Bond paths also appear for weak interactions: hydrogen bonds, halogen bonds, van der Waals contacts, metal–ligand bonds.
- **Interpretation.** Bader first described a bond path as the necessary and sufficient condition for two atoms to be bonded. Later, prompted by paths between atoms that are not bonded in any chemical sense (for example the H···H paths in planar biphenyl), he described it instead as a *privileged exchange channel* — a line of accumulated density shared by two atoms — and stressed that a bond path is not the same as a chemical bond (Bader, 2009). Read a bond path between weakly interacting atoms as evidence of an interaction, then judge its nature from the properties below.
- Topologies can change abruptly with geometry: a BCP and an RCP that approach each other merge and disappear (a *catastrophe*). Rings whose ρ at the RCP is close to ρ at a BCP, or ring and cage points with almost equal density, are close to such a point and fragile.

### 5. The Laplacian of the density

The Laplacian ∇²ρ = λ₁ + λ₂ + λ₃ measures whether charge is locally **concentrated** (∇²ρ < 0: ρ is larger than its average in the immediate neighbourhood) or **depleted** (∇²ρ > 0). Around an atom, −∇²ρ shows the electron shells, and in the valence shell its local maxima correspond to the bonded and non-bonded (lone) pairs of the Lewis and VSEPR models.

At a BCP, the sign of ∇²ρ tells which curvatures dominate:

- **Shared-shell (covalent) interactions**: charge is concentrated between the nuclei. ρ<sub>BCP</sub> is large (typically > 0.1 a.u.) and ∇²ρ < 0, because |λ₁| and |λ₂| (contraction towards the bond path) exceed λ₃; η = |λ₁|/λ₃ > 1.
- **Closed-shell interactions** (ionic bonds, hydrogen bonds, van der Waals contacts): charge is depleted at the BCP and concentrated in each atom separately. ρ<sub>BCP</sub> is small (typically < 0.05 a.u.), ∇²ρ > 0, and η < 1.

Koch and Popelier proposed that a hydrogen bond has a BCP with ρ between 0.002 and 0.034 a.u. and ∇²ρ between 0.024 and 0.139 a.u.

For a given pair of elements, ρ<sub>BCP</sub> grows with bond order. In the program's test molecules (B3LYP) the C–C value is 0.233 a.u. in bicyclo[1.1.1]pentane (single bonds), 0.311 in benzene (aromatic) and 0.326 in ethylene (double bond, 6-31G).

### 6. Energy densities

Two kinetic energy densities are used, because the kinetic energy density is not unique:

- the **Lagrangian** G(**r**) = ½ Σ<sub>i</sub> n<sub>i</sub> |∇φ<sub>i</sub>|², which is never negative;
- the **Hamiltonian** K(**r**) = −¼ Σ<sub>i</sub> n<sub>i</sub> (φ<sub>i</sub>*∇²φ<sub>i</sub> + φ<sub>i</sub>∇²φ<sub>i</sub>*), which can be negative.

They differ by K = G − ¼∇²ρ, which integrates to zero over a zero-flux basin, so both give the same atomic kinetic energy. The **local virial theorem** relates G and the potential energy density V(**r**) to the Laplacian:

$$
\tfrac14 \nabla^2\rho(\mathbf{r}) = 2G(\mathbf{r}) + V(\mathbf{r}) .
$$

G > 0 and V < 0 everywhere. A negative Laplacian therefore means that the potential energy dominates locally, a positive one that the kinetic energy does. The **total energy density**

$$
H(\mathbf{r}) = G(\mathbf{r}) + V(\mathbf{r}) = -K(\mathbf{r})
$$

is negative where the potential energy dominates. At a BCP:

- **Cremer and Kraka** criterion: H < 0 indicates a covalent contribution, even when ∇²ρ > 0. This is the case for strong hydrogen bonds and many metal–ligand bonds.
- **Espinosa et al.** use the ratio |V|/G:
  - below 1, a pure closed-shell interaction (∇²ρ > 0 and H > 0);
  - between 1 and 2, intermediate (∇²ρ > 0 but H < 0), partly covalent;
  - above 2, a shared-shell (covalent) interaction (∇²ρ < 0).
- H/ρ (the "bond degree") and G/ρ are also used. G/ρ is well below 1 for covalent bonds (0.2 for the O–H bonds of the water dimer) and near or above 1 for closed-shell interactions (0.8 for its hydrogen bond).

In the water dimer ([Example 9](#example-9-qtaim-analysis-of-a-hydrogen-bond-water-dimer)) the O–H bonds have |V|/G ≈ 9, and the hydrogen bond 1.10 with a slightly negative H: intermediate.

The **hydrogen-bond energy** has been estimated from BCP properties. Espinosa, Molins and Lecomte proposed BE ≈ V/2. Emamian, Lu *et al.* (2019) showed that a linear relation in ρ<sub>BCP</sub> is more reliable: BE ≈ −223.08 ρ + 0.7423 kcal/mol for neutral hydrogen bonds and BE ≈ −332.34 ρ − 1.0661 kcal/mol for charged ones. Both are estimates; the program reports ρ and V, not BE.

### 7. Ellipticity

At a BCP, the **ellipticity**

$$
\varepsilon = \frac{\lambda_1}{\lambda_2} - 1 \ge 0
$$

measures how much the density is drawn out in one direction perpendicular to the bond (λ₁ is the more negative curvature).

- ε = 0 for a cylindrically symmetric bond: single σ bonds and triple bonds.
- ε grows with π character in one plane. In the test molecules (B3LYP) ε is about 0.01–0.02 for C–H bonds, 0.209 for the C–C bonds of benzene and 0.252 for the C=C bond of ethylene (6-31G).
- A large ε at a weak BCP signals a structurally unstable bond, close to the catastrophe described in section 4.

The direction of the λ₂ eigenvector (the softer curvature) gives the orientation of the π system.

### 8. Ring and cage points

ρ at an RCP grows with the electron delocalization in the ring, and has been used as an aromaticity indicator (Howard and Krygowski, 1997). They found that the curvature of ρ perpendicular to the ring plane correlates even better. Benzene at B3LYP/6-31G(d) has ρ<sub>RCP</sub> = 0.0200 a.u. with ∇²ρ = +0.162 a.u. A CCP shows that rings enclose a cage, such as the three four-membered rings of bicyclo[1.1.1]pentane.

### 9. Practical considerations

- **Method and basis.** ρ at the BCPs of polar bonds and weak interactions depends on the basis set and, less, on the method. Compare values computed at the same level. Use polarization functions; add diffuse functions for anions and weak complexes.
- **Gaussian basis sets** have no cusp at the nucleus, so the (3,−3) point of a hydrogen atom lies a few hundredths of a bohr from the nucleus (about 0.06 bohr in water at 6-31G(d)). The program — like Multiwfn — reports the CP where it is, and measures bond-path ratios between the CPs.
- **Effective core potentials** remove the core electrons from the density. The (3,−3) point of such an atom and properties close to it are not meaningful; bond points farther out are usually still usable.
- **One structure.** The topology belongs to one geometry. For flexible or fluxional systems, analyse the structures that matter, such as minima, transition states or snapshots.
- **Correlated densities** (MP2, CCSD) can be analysed in principle, but this program reads the SCF density only.

### 10. References

- R. F. W. Bader, *Atoms in Molecules: A Quantum Theory*, Oxford University Press (1990).
- R. F. W. Bader, "A quantum theory of molecular structure and its applications", *Chem. Rev.* 91, 893 (1991).
- P. L. A. Popelier, *Atoms in Molecules: An Introduction*, Prentice Hall (2000).
- C. F. Matta, R. J. Boyd (eds.), *The Quantum Theory of Atoms in Molecules*, Wiley-VCH (2007).
- R. F. W. Bader, "Bond paths are not chemical bonds", *J. Phys. Chem. A* 113, 10391 (2009).
- D. Cremer, E. Kraka, *Angew. Chem. Int. Ed. Engl.* 23, 627 (1984) (energy density criterion).
- U. Koch, P. L. A. Popelier, *J. Phys. Chem.* 99, 9747 (1995) (hydrogen-bond criteria).
- E. Espinosa, E. Molins, C. Lecomte, *Chem. Phys. Lett.* 285, 170 (1998) (BE ≈ V/2).
- E. Espinosa, I. Alkorta, J. Elguero, E. Molins, *J. Chem. Phys.* 117, 5529 (2002) (|V|/G regions).
- S. Emamian, T. Lu, H. Kruse, H. Emamian, *J. Comput. Chem.* 40, 2868 (2019) (hydrogen-bond energies from ρ<sub>BCP</sub>).
- S. T. Howard, T. M. Krygowski, *Can. J. Chem.* 75, 1174 (1997) (ring critical points and aromaticity).
- T. Lu, F. Chen, *J. Comput. Chem.* 33, 580 (2012); T. Lu, *J. Chem. Phys.* 161, 082503 (2024) (Multiwfn, whose topology algorithms the program follows).

## Formulas and Methods (QTAIM tab)

### 1. Electron density and its derivatives

From the Cartesian density matrix P of `Total SCF Density` and the
normalised Cartesian Gaussians χ (the same basis conventions as the ESP and
Fukui tools):

$$
\rho = \sum_{ab} P_{ab}\,\chi_a\chi_b, \qquad
\nabla\rho = 2\sum_{ab} P_{ab}\,\chi_a\nabla\chi_b, \qquad
\partial_i\partial_j\rho = 2\sum_{ab} P_{ab}\left(\partial_i\chi_a\,\partial_j\chi_b + \chi_a\,\partial_i\partial_j\chi_b\right)
$$

The AO derivatives are analytic. With `numba` installed, one compiled kernel
evaluates them and sums P only over the AOs that are not screened out at each
point; it agrees with the numpy path to 10⁻¹⁵ relative.

### 2. Critical points

A critical point is a point where ∇ρ = 0. Its type (rank, signature) comes from
the eigenvalues λ₁ ≤ λ₂ ≤ λ₃ of the Hessian: three negative → (3,−3), a maximum
(nuclear, NCP); two → (3,−1) (bond, BCP); one → (3,+1) (ring, RCP); none →
(3,+3) (cage, CCP). A converged point with a zero eigenvalue (degenerate) is
counted in a note and not reported. For an isolated molecule all CPs satisfy
the Poincaré–Hopf relation

$$
n_{(3,-3)} - n_{(3,-1)} + n_{(3,+1)} - n_{(3,+3)} = 1 .
$$

It is a necessary check, not a sufficient one.

### 3. Search (Multiwfn's Newton method)

Starting points: the nuclei (mode 2), the midpoints of atom pairs (3), and the
centres of three (4) and four (5) atoms. A combination is skipped if any two of
its atoms are farther apart than 1.5 × the sum of their vdW radii (Bondi /
Mantina, `help.py`). From each start the Newton step

$$
\Delta\mathbf{r} = -s\,\mathbf{H}^{-1}\nabla\rho
$$

(scale s = 1, optionally capped at a trust radius) is repeated until
|Δr| < 10⁻⁷ bohr and |∇ρ| < 10⁻⁶ a.u. (the gradient test is waived for
nuclear starts, whose peaks are too sharp), for at most 120 iterations. The
start is abandoned when |det H| < 5×10⁻²² or ρ < 10⁻⁸ a.u. A converged point
closer than 0.03 bohr to one found before (in mode order) is a duplicate. A
(3,−3) more than 0.3 bohr from every nucleus is reported as a non-nuclear
attractor (NNA). With Gaussian basis sets an H nucleus' (3,−3) usually sits
~0.06 bohr off the nucleus, in Multiwfn too.

### 4. Properties at a critical point (Multiwfn's definitions)

$$
G = \tfrac12\sum_{ab} P_{ab}\,\nabla\chi_a\cdot\nabla\chi_b, \qquad
K = G - \tfrac14\nabla^2\rho, \qquad
V = -G - K = \tfrac14\nabla^2\rho - 2G, \qquad
H = G + V = -K
$$

G is the Lagrangian (positive-definite) kinetic energy density, K the
Hamiltonian one, V the potential energy density from the local virial theorem
and H the total energy density. Also reported: the ellipticity
ε = λ₁/λ₂ − 1, η = |λ₁|/λ₃ (Multiwfn's *eta index*, negative at a (3,−3)),
|V|/G and H/ρ. As a usual reading at a BCP: ∇²ρ < 0 and H < 0 indicate a shared
(covalent) interaction; ∇²ρ > 0 with |V|/G < 1 a closed-shell one (ionic,
hydrogen bond, van der Waals); 1 < |V|/G < 2 the intermediate region.

### 5. Gradient paths

Paths are followed with classical fourth-order Runge–Kutta on the unit vector
±∇ρ/|∇ρ| with a fixed step h = 0.03 bohr, so their points are evenly spaced,
as in Multiwfn. A path ends when it comes within 0.05 bohr of its target CP (the
CP is then its last point) and is dropped after 451 points.

- **Bond paths**: from each (3,−1), displaced by ±h along the eigenvector of
  λ₃, uphill to a (3,−3). The two (3,−3) reached give the BCP its atom pair.
  Each half's length is the sum of its segments.
- **Ring paths** (only when the molecule has a (3,+3)): from each (3,+1),
  displaced by ±h along the eigenvector of λ₁, downhill. A half that reaches a
  (3,+3) is kept; the other leaves to infinity, which is normal.

### 6. The 2D plane map

The plane is that of three atoms (origin at their centroid, u along atom 1 →
atom 2, v in the plane, orthonormal) or XY/XZ/YZ at an offset. The ρ grid
covers the projected atoms plus a margin. Gradient lines start on a 0.1-bohr
circle around each nucleus lying within 0.5 bohr of the plane and descend
along the in-plane projection of −∇ρ. They end at an in-plane (3,−1)/(3,+1),
at the map border, or after 2000 points. Interbasin paths descend the same
way from each in-plane (3,−1), along the in-plane projection of its
negative-curvature eigenvector (λ₁ or λ₂, whichever lies more in the plane).
CPs and bond- or ring-path points within 0.5 bohr of the plane are projected
onto it (Multiwfn's default display distance).

### 7. References

R. F. W. Bader, *Atoms in Molecules: A Quantum Theory* (Oxford, 1990);
T. Lu, F. Chen, *J. Comput. Chem.* 33, 580 (2012) and T. Lu, *J. Chem. Phys.*
161, 082503 (2024) (Multiwfn); Multiwfn manual §3.14 (topology analysis) and
§3.5.5 (CPs and paths on plane maps); E. Espinosa, I. Alkorta, J. Elguero,
E. Molins, *J. Chem. Phys.* 117, 5529 (2002) (the |V|/G regions of section 4).

## Formulas and Methods (NTO tab)

### 1. Excited states of single-reference methods

In CIS and TDA the excited state is a linear combination of singly excited
determinants, each moving one electron from an occupied MO *i* to a virtual MO
*a*, with configuration coefficients $w_i^a$:

$$
\Psi_{\text{exc}} = \sum_{i}^{\text{occ}}\sum_{a}^{\text{vir}} w_i^a\,\Phi_i^a .
$$

TDHF and TD-DFT add de-excitation coefficients $w_i'^a$ (Gaussian's `<-` lines),
with $\sum (w_i^a)^2 - \sum (w_i'^a)^2 = 1$. For a closed-shell reference
Gaussian prints one set of coefficients normalised to 0.5 instead of 1 (the α and
β halves are identical). An unrestricted reference has separate α and β
coefficients (`A`/`B` labels) normalised together to 1.

### 2. Transition matrix

For each spin, the transition matrix is the occupied × virtual array of the
excitation coefficients (Multiwfn's definition; de-excitations are left out, as
Gaussian's `Pop=NTO` also does):

$$
T_{ia} = w_i^a .
$$

The report repeats Multiwfn's check:
$\sum (w_i^a)^2$, $-\sum (w_i'^a)^2$, their sum and its deviation from 0.5
(restricted) or 1 (unrestricted).

### 3. Singular value decomposition

$$
T = U\,\Sigma\,V^{\mathsf T}, \qquad
\Sigma = \operatorname{diag}(\sigma_1 \ge \sigma_2 \ge \dots \ge 0),
$$

with $U$ (occupied × occupied) and $V$ (virtual × virtual) orthogonal. The
squared singular values are the eigenvalues of $TT^{\mathsf T}$ and of
$T^{\mathsf T}T$, which are the two matrices Multiwfn diagonalises separately.
The SVD gives the same values and pairs every particle with its hole by
construction: $v_k = T^{\mathsf T}u_k/\sigma_k$.

### 4. Holes, particles and weights

With $C_{\text{occ}}$ and $C_{\text{vir}}$ the coefficient matrices of the
occupied and virtual MOs (basis functions × MOs):

$$
\phi^{\text{hole}}_k = C_{\text{occ}}\,u_k, \qquad
\phi^{\text{particle}}_k = C_{\text{vir}}\,v_k, \qquad
\lambda_k = f\,\sigma_k^2,
$$

with $f = 2$ for a restricted reference (undoing Gaussian's 0.5
normalisation) and $f = 1$ for each spin of an unrestricted one. $100\,\lambda_k$
is the percentage of the excitation described by pair *k*. The NTOs are
orthonormal because $U$ and $V$ are orthogonal and the MOs are orthonormal.
The orbitals are rotated only among the occupied and among the virtual MOs,
so the ground-state density is unchanged. That is why the NTO `.fchk` keeps
its `Total SCF Density` section as it is. The composition of NTO *k* in
the canonical MOs is $u_{ik}^2$ (hole) and $v_{ak}^2$ (particle).

For CIS and TDA, $0 \le \lambda_k \le 1$ and $\sum_k\lambda_k = 1$. For TD the
omitted de-excitations let $\sum_k \lambda_k$ exceed 1 by $f\sum(w_i'^a)^2$.

### 5. Signs and degeneracy

Each pair can change sign as a whole, $(u_k, v_k) \to (-u_k, -v_k)$. This
program makes the largest-magnitude element of $u_k$ positive. When several
pairs share one weight, any orthogonal rotation of them is an equally valid
set; only the space they span is defined. The tests compare with Multiwfn
accordingly: sign-free for single pairs, and by subspace (an orthogonal
least-squares map) for degenerate groups.

### 6. Orbitals on a grid

The orbital value at a point is
$\phi(\mathbf r) = \sum_\mu c_\mu \chi_\mu(\mathbf r)$, evaluated exactly at
every grid point. The basis functions χ are the same normalised Cartesian
Gaussians as in the ESP, Fukui and QTAIM tools, and pure shells are mapped onto
Cartesians. The AOs are computed once per block of points (by the compiled
kernel when `numba` is installed) and multiplied by the hole and particle
coefficients together.

The grid starts from the atoms ± 1.7 van der Waals radii (Multiwfn's box). Any
face on which |hole| or |particle| still reaches isovalue/10 is pushed out by
1 bohr until none does. A pair on uracil at 0.2 bohr (≈ 0.4 million points) takes
about 0.3 s with `numba`. The isosurfaces at ±isovalue are triangulated as in
the Cube tab.

### 7. UV-Vis spectrum

*Plot Spectrum* uses the band shape of GaussView and of Gaussian's note
*UV/Visible Plots in Gaussian*: a Gaussian in wavenumber for each state *i*
with excitation energy $\tilde\nu_i$ and oscillator strength $f_i$,

$$
\varepsilon(\tilde\nu) = \sum_i 1.3062974\times10^{8}\,\frac{f_i}{\sigma}
\exp\!\left[-\left(\frac{\tilde\nu-\tilde\nu_i}{\sigma}\right)^2\right],
$$

with $\tilde\nu$, $\tilde\nu_i$ and σ in cm⁻¹ and ε in L mol⁻¹ cm⁻¹. The half
width at half height is $\sigma\sqrt{\ln 2}$, so the default 0.333 eV is
σ = 0.400 eV (3226 cm⁻¹). A band of strength f therefore peaks at
ε = 1.3062974 × 10⁸ f/σ. Wavelength, energy and wavenumber are converted with
λ (nm) = 1239.841984 / E (eV) and $\tilde\nu$ (cm⁻¹) = 10⁷ / λ. The curve is
sampled at 2000 wavelengths over a range that covers every band with f > 0
± 4σ and every dark state's stick ± 1σ. Sticks are drawn at the wavelengths
Gaussian prints.

### 8. References

- R. L. Martin, *J. Chem. Phys.* **2003**, 118, 4775 (natural transition orbitals).
- T. Lu, F. Chen, *J. Comput. Chem.* **2012**, 33, 580, and the Multiwfn manual,
  §3.21 and §4.18.6.

## Formulas and Methods (Solvation box builders)

This section documents the equations used by the three solvation-box builders
— [single_solute_solvent_box.py](single_solute_solvent_box.py),
[mixtureSolventBox.py](mixtureSolventBox.py) and
[sharedWallBox.py](sharedWallBox.py) — and by the modules they share:
[solventPacking.py](solventPacking.py) (packing),
[periodicPreOptimizer.py](periodicPreOptimizer.py) (pre-optimization) and
[systemProperties.py](systemProperties.py) (electrons and multiplicity).
They are exercised headlessly by the `tests/test_*` files of the same names.

**Symbols used throughout.**

| Symbol | Meaning |
| --- | --- |
| $\mathbf{L}=(a,b,c)$ | orthorhombic cell edges, Å (shared wall: $c = c_1 + c_2$) |
| $\rho$ | target density, g/cm³ |
| $m_u$ | atomic mass unit, $1.66053906660\times10^{-24}$ g |
| $M_X$ | molar mass of species $X$ (sum of atomic masses) |
| $r_i$ | van der Waals radius of atom $i$ |
| $s$ | packing vdW scale (the *vdW scale* field) |
| $s_c$ | target contact scale of the pre-optimization |
| $\delta$ | extra solute clearance (*Extra solute clearance* / *Solute gap*) |
| $\mathbf{d}_{ij}$ | minimum-image vector between atoms $i$ and $j$ |

### 1. Target number of solvent molecules

$$
N_\text{single} = \left\lfloor \frac{\rho V / m_u - M_\text{solute}}{M_\text{solvent}} \right\rfloor,
\qquad
N_\text{mixture} = \operatorname{round}\!\left( \frac{\rho V / m_u - M_\text{solute}}{x_A M_A + x_B M_B} \right),
\qquad V = abc \times 10^{-24}\ \text{cm}^3
$$

- **Represents:** how many solvent molecules bring the cell to the target density.
- **Where/how used:** `estimate_target_num_solvent_molecules` (single, rounds
  down) and `estimate_target_molecule_count` (mixture and each shared-wall
  region with its own slab volume $a\,b\,c_k$; rounds to nearest). The two
  rounding rules are historical and are kept so that seeded boxes do not
  change. $x_A, x_B$ are the composition fractions.

### 2. Minimum image

$$
\mathbf{d}_{ij} = \Delta\mathbf{r} - \mathbf{L} \odot \operatorname{round}(\Delta\mathbf{r} \oslash \mathbf{L}),
\qquad \Delta\mathbf{r} = \mathbf{r}_i - \mathbf{r}_j
$$

- **Represents:** the shortest separation between two atoms over all periodic
  images, component by component; exact for separations below half the
  smallest cell edge. `round` rounds half to even.
- **Where/how used:** every clash test while packing (when *Periodic
  minimum-image clash detection* is on) and every intermolecular distance in the
  pre-optimization. With the switch off, $\mathbf{d}_{ij} = \Delta\mathbf{r}$.

### 3. Packing acceptance

A randomly placed (and optionally randomly rotated) candidate is accepted only
if every atom pair with the solute and with molecules already accepted
satisfies

$$
|\mathbf{d}_{ij}| \ge s\,(r_i + r_j) \;(+\,\delta \text{ for solvent–solute pairs}),
$$

and every candidate atom lies inside the box (minus the wall padding; in the
shared-wall box, inside its own slab, with padding on the outer faces only).
Molecules whose centres are further apart than
$R_\text{cand} + R_m + 2.5$ Å are not compared atom by atom, where $R$ is the
molecule's bounding radius (single: including the atoms' vdW radii; mixture
and shared wall: geometric only).

- **Where/how used:** the *Legacy* engine evaluates it with the original Python
  loops; the *NumPy* and *Numba* engines evaluate the **same** inequality with a
  cell list and vectorised or compiled kernels, in the same floating-point
  order, so all engines accept and reject exactly the same candidates and a
  fixed seed gives the same box.
- **Mixture composition:** the next molecule is solvent A with probability
  $p_A = n_A / (n_A + n_B)$, $n_X = \max(0,\, x_X - f_X)$, where $f_X$ is the
  fraction of $X$ placed so far. The first molecule, and any molecule placed
  while both needs are zero, is A with probability $x_A$. The shared-wall
  builder skips the random draw entirely when a region's composition is 0 %
  or 100 %.

### 4. Final density

$$
\rho_\text{final} = \frac{m_u \sum_{i} m_i}{V}
$$

over all placed atoms (per region for the shared-wall box).

### 5. Pre-optimization: overlap penalty

$$
E_\text{overlap} = K \sum_{\substack{i<j \\ \text{different molecules} \\ |\mathbf{d}_{ij}| < \sigma_{ij}}}
\left( |\mathbf{d}_{ij}|^2 - \sigma_{ij}^2 \right)^2,
\qquad
\sigma_{ij} = s_c\,(r_i + r_j) \;(+\,\delta \text{ for solute–solvent pairs}),
\qquad K = 1
$$

$$
\frac{\partial E_\text{overlap}}{\partial \mathbf{r}_i} =
\sum_j 4K \left( |\mathbf{d}_{ij}|^2 - \sigma_{ij}^2 \right) \mathbf{d}_{ij}
$$

- **Represents:** a penalty that is exactly zero once no intermolecular pair is
  closer than $\sigma$ — the same form of criterion as packing, with the larger
  target contact scale $s_c$ (default $\max(0.90, s)$; packing already
  guarantees $s$, so relaxing to $s$ itself would do nothing).
- **Where/how used:** `OverlapPenalty`. Pairs come from a periodic k-d tree
  (`scipy.spatial.cKDTree` with `boxsize = L`), which requires
  $\max \sigma < \tfrac12 \min(a,b,c)$. Atoms of the same molecule never
  interact here. The minimizer works with $\sigma(1 + 0.01)$ so that contacts
  settle just outside $\sigma$; the reported penalty, clash count
  ($|\mathbf{d}_{ij}| < \sigma_{ij}$) and closest contact use $\sigma$ itself.

### 6. Pre-optimization: walls and regions

Without periodicity a harmonic wall keeps atoms inside the box:

$$
E_\text{wall} = K_w \sum_{i,\alpha} \left[ \min(x_{i\alpha} - \ell_{i\alpha},\, 0)^2 + \max(x_{i\alpha} - h_{i\alpha},\, 0)^2 \right],
\qquad K_w = 10
$$

with $\ell = $ padding and $h = L - $ padding (shared wall: per box on the x/y
faces, and in z only on the outer faces, never at the interface).

In the shared-wall box a flat-bottom restraint keeps the centre of mass $Z_m$ of
each molecule in the region it was packed into, $[z^\text{lo}_m, z^\text{hi}_m]$
= $[0, c_1]$ or $[c_1, c_1 + c_2]$:

$$
E_\text{slab} = K_w \sum_m \left[ \min(Z_m - z^\text{lo}_m,\, 0)^2 + \max(Z_m - z^\text{hi}_m,\, 0)^2 \right],
\qquad Z_m = \frac{\sum_{i \in m} m_i z_i}{\sum_{i \in m} m_i}
$$

### 7. Pre-optimization: rigid molecules

Each molecule $m$ has a translation $\mathbf{t}_m$ and a rotation vector
$\boldsymbol{\omega}_m$; its atoms are
$\mathbf{x}_i = \mathbf{R}(\boldsymbol{\omega}_m)\,\mathbf{b}_i + \mathbf{t}_m$,
with $\mathbf{b}_i$ the atom's offset from the molecule's centroid in the packed
box. With $\mathbf{g}_i = \partial E / \partial \mathbf{x}_i$:

$$
\frac{\partial E}{\partial \mathbf{t}_m} = \sum_{i \in m} \mathbf{g}_i,
\qquad
\frac{\partial E}{\partial \boldsymbol{\omega}_m} =
\boldsymbol{\tau} - A\,\boldsymbol{\omega} \times \boldsymbol{\tau}
+ B\,\boldsymbol{\omega} \times (\boldsymbol{\omega} \times \boldsymbol{\tau}),
\qquad
\boldsymbol{\tau} = \sum_{i \in m} (\mathbf{x}_i - \mathbf{t}_m) \times \mathbf{g}_i
$$

$$
A = \frac{1 - \cos\theta}{\theta^2}, \qquad B = \frac{\theta - \sin\theta}{\theta^3}, \qquad \theta = |\boldsymbol{\omega}|
$$

(Taylor series $A \approx \tfrac12 - \theta^2/24$, $B \approx \tfrac16 - \theta^2/120$
below $\theta = 10^{-3}$.) This is the transposed left Jacobian of SO(3) applied
to the torque. Bond lengths and angles are preserved to round-off.

### 8. Pre-optimization: flexible molecules

Every atom moves. Within each molecule, bonds are detected once on the input
geometry ($r \le r^\text{cov}_i + r^\text{cov}_j + 0.45$ Å), and the input
geometry itself is the minimum:

$$
E_\text{bond} = 300 \sum (r - r_0)^2, \qquad
E_\text{angle} = 150 \sum (\cos\theta - \cos\theta_0)^2 \;\; \text{or} \;\; 80 \sum (1 + \cos\theta) \text{ if } \theta_0 \ge 175^\circ,
$$

$$
E_\text{1-5} = 0.1 \sum \left( \frac{\sigma_{ij}}{r_{ij}} \right)^{12},
\qquad \sigma_{ij} = \min\!\left(1.8\,(r^\text{cov}_i + r^\text{cov}_j),\; 0.8\, r^\text{input}_{ij}\right)
$$

for pairs at least four bonds apart; $r_0$, $\theta_0$ and $r^\text{input}$
are taken from the input structure. There are no torsion terms. These are the
terms of the molecular-design pre-optimizer, built by
`ForceField.from_geometry` instead of from bond orders.

### 9. Minimization and wrapping

$E = E_\text{overlap} + E_\text{wall} + E_\text{slab}$ (+ the flexible terms)
is minimized with L-BFGS-B (`scipy.optimize.minimize`) up to the *Max
iterations* limit; reaching the limit is reported, not treated as an error.
Afterwards every molecule is translated as a whole:

$$
\mathbf{x}_i \leftarrow \mathbf{x}_i - \mathbf{L} \odot \left\lfloor \mathbf{X}_m \oslash \mathbf{L} \right\rfloor
\quad (i \in m),
$$

where $\mathbf{X}_m$ is the molecule's centre of mass, so that every centre of
mass lies in $[0, \mathbf{L})$ and no molecule is split across a face.

### 10. Electrons and spin multiplicity

$$
N_e = \sum_i Z_i - Q, \qquad
Q = q_\text{solute} + N_A\,q_A + N_B\,q_B
$$

(single solute: one solvent species; shared wall: summed over both regions),
where $q$ are the charges per molecule and $N_A, N_B$ the numbers of molecules
actually inserted. A multiplicity $M = 2S + 1$ is accepted when

$$
N_e \ge 0, \qquad (N_e + M) \text{ is odd}, \qquad 1 \le M \le N_e + 1 .
$$

$M$ is never derived from $N_e$: an even $N_e$ is compatible with
$M = 1, 3, 5, \ldots$ and an odd one with $M = 2, 4, 6, \ldots$, so the user's
value is checked and an impossible one is reported, not replaced.

## Troubleshooting

### "No file was selected!"

You closed a file dialog without choosing a file. Click **Browse** again and select the required input.

### "Please input a valid value for ..."

One of the numeric fields is blank or contains text that cannot be converted to a number. Check the field named in the message and enter a valid value.

### "Invalid format for ... Please input exactly two atom labels."

The atom-label field has the wrong number of labels. Bond analysis needs two labels, bond angle needs three, and dihedral angle needs four. Separate labels with spaces, not commas.

### "Failed to open file" or "Failed to read file"

The selected file could not be opened or parsed. Confirm that:

- The file still exists.
- You have permission to read it.
- The file format matches the selected tool.
- XYZ files begin with the atom count and contain valid coordinate lines.

### "The file is empty!"

The selected file has no readable content. Choose a valid geometry, trajectory, or input file.

### "Invalid line format in TRAJECTORY file"

The trajectory file does not match the expected numeric layout. Check that the file came from the expected simulation program and has not been truncated or edited incorrectly.

### "The stop frame must be less than the total number of frames"

The selected stop frame is outside the trajectory. Use the frame count displayed after browsing for the trajectory, then choose a smaller stop frame.

### ENERGY file format errors

The CPMD energy plotter expects a valid ENERGY file with eight numeric columns. If plotting fails, verify that the file is complete and that you selected the correct plot type.

### Solvent box insertion warnings

If the solvent box builder cannot insert all requested molecules, the target density or minimum-distance settings may be too restrictive. Try one or more of the following:

- Increase the box dimensions.
- Lower the target density.
- Reduce the minimum distance or van der Waals scaling factor.
- Increase the maximum number of insertion attempts.

Random insertion becomes slow as the box fills, because most late attempts are
rejected. With `numba` installed the *Automatic* packing engine is 20–50 times
faster than the original loop, so a larger number of attempts is affordable;
every engine produces the same box for the same random seed. The
pre-optimization cannot add molecules — it only relaxes the ones that were
placed.

### "Pre-optimization not performed"

The packed box was written, but the optional pre-optimization was skipped. The
dialog states the reason. The usual one is:

*The contact range (… A) must be shorter than half the smallest cell length.*
Under periodic boundaries every distance uses the minimum image, which is only
valid for contacts shorter than half the smallest cell edge. The contact range
is the largest `target contact scale × (r_i + r_j)` (plus the solute clearance).
Use a larger box, or a smaller target contact scale.

If the summary instead reports clashing pairs left after the pre-optimization,
the box is too dense for the requested contact scale: lower the density or the
target contact scale, or switch from *Rigid* to *Flexible* molecules.

### "Spin multiplicity" or "Check the spin multiplicity"

The total spin multiplicity (2S + 1) is required and is never guessed from the
electron count.

- **Missing or not a positive integer** (before the build): the Setup tab opens
  with the cursor in **Total spin multiplicity**. In the single-solute and
  mixture builders it is in the **Charge and spin** section; in the
  shared-wall builder it is at the top, in **Shared geometry**.
- **Impossible for this box** (after the build): an even electron count allows
  1, 3, 5, …; an odd one allows 2, 4, 6, …; and the value cannot exceed the
  number of electrons + 1. The message lists the allowed values. Correct the
  field (or the charges) and press **Update properties** on the Summary tab —
  the summary file is rewritten and the box is not rebuilt.

### Missing `SH_STATE.dat` or `TRAJEC.xyz`

The SH Geometry Analyzer skips subfolders that do not contain both files. Add the missing files or remove incomplete subfolders from the selected root directory.

### Plot windows do not appear

A tool reports that it finished, but no plot window opens and no error is
shown. The interactive plot viewer runs in a separate process and pins
matplotlib's `TkAgg` backend, so this almost always means **tkinter is
missing**. It fails silently rather than raising.

Check for it:

```powershell
python -c "import tkinter; print(tkinter.TkVersion)"
```

If that command fails, re-run the python.org installer, choose **Modify**, and
enable the **tcl/tk and IDLE** component. tkinter is a standard-library module
and cannot be installed with pip. See
[System requirements that pip cannot install](#system-requirements-that-pip-cannot-install).

In a packaged `.exe` build, the same symptom means the `TkAgg` backend was not
bundled. See the packaging notes in `readme.md`.

### Molecular viewer does not open or shows no 3D view

The viewer depends on OpenGL and GLFW. Work through these in order:

1. Confirm the packages are installed:

   ```powershell
   python -c "import glfw, OpenGL; print('ok')"
   ```

   If this fails, reinstall with
   `python -m pip install -r installation_requirements.txt`. The required
   `glfw3.dll` and `freeglut` libraries are bundled inside those wheels, so
   they do not need to be downloaded separately.

2. Confirm your graphics driver provides OpenGL in the current desktop
   session. Install the vendor driver for your hardware (Intel, NVIDIA or
   AMD); the Microsoft Basic Display Adapter driver does not include an
   OpenGL runtime.

3. The viewer cannot run over a plain SSH session or in a headless container.
   A Remote Desktop session may work through a software renderer, but slowly.

If the window still does not open, the failure reason is reported in a dialog
and in the status line at the bottom of the Molecular Viewer window.

### ESP tool messages

Errors from the Electrostatic Surface Potential tool appear in a dialog; notes
appear as `NOTE:` lines at the end of the results summary.

| Message | Cause | What to do |
|---|---|---|
| `… is missing fchk section(s): …` | The file is not a complete formatted checkpoint (or is a `.log`/`.chk`). | Re-create it with `formchk name.chk name.fchk`. |
| `Section '…' declares N values but holds M; the file looks truncated.` | The `.fchk` was cut short (copy or disk problem). | Re-run `formchk`. |
| `basis functions with l = 5 are not supported (highest supported is g, l = 4).` | The basis contains h (or higher) functions. | Use a basis set without h functions. |
| `Isovalue … is outside (0, 1) a.u. It is an electron-density value …` | An ESP-like or negative number was typed as isovalue. | Type a density value, e.g. `0.0004` or `0.001`. |
| `No isosurface at rho = … a.u. inside the grid.` | The isovalue is larger than the density anywhere on the grid. | Use a smaller isovalue. |
| `Could not enclose the isosurface: density stays above … on the grid boundary.` | The isovalue is so small that the surface keeps growing (very diffuse basis, anions). | Use a larger isovalue. |
| `NOTE: rho = …: surface components with V+F-E = [...] (2 expected …)` | A surface component is not a closed sphere-like surface — usually a grid too coarse for a thin neck between fragments. | Repeat with a smaller grid spacing (e.g. 0.15). |
| `NOTE: cubegen was not found; …` / `NOTE: cubegen failed: …` | The cubegen engine was chosen but cubegen is missing, or it failed on this file (for example a file written by a newer Gaussian than the cubegen). | Nothing required: the built-in engine was used. Fix the cubegen path to use cubegen. |
| `Could not write …` | The output folder is read-only or missing. | Choose another output folder. |
| `The 3D viewer could not be opened` | OpenGL/GLFW problem, as for the 3D molecular viewer. | See [Molecular viewer does not open or shows no 3D view](#molecular-viewer-does-not-open-or-shows-no-3d-view). |

### NTO tab messages

Errors from the NTO tab appear in a dialog; notes appear under `NOTES` in the
report `name_state<n>_NTO.txt`.

| Message | Cause | What to do |
| --- | --- | --- |
| `Press Load to read these files and list their excited states.` | The file fields changed since the last **Load**. | Press **Load** (or **Run**, which loads first). |
| `… and … describe different atoms` / `… have different geometries (atom n differs by …). Use the .fchk of the same job.` | The `.fchk` and the output belong to different molecules or geometries (e.g. a reoptimized `.fchk`). | Use the `.fchk` converted from the checkpoint of the same TD job. |
| `… has N basis functions but … has M` / `… alpha electrons …` | Same molecule, different basis set, charge or multiplicity. | Use the two files of one job. |
| `… is unrestricted (A/B labels) but … is restricted` (or the reverse) | The output and the `.fchk` come from a closed-shell and an open-shell calculation. | Use the two files of one job. |
| `… is a restricted open-shell (RO…) reference` | ROHF/ROKS-based TD is not supported. | Use an unrestricted (UHF/UKS) reference. |
| `… is a spin-flip pair; spin-flip excitations are not supported.` | A spin-flip TD output. | Not supported. |
| `… is a Gaussian output without excited states.` | The `.log` is of a ground-state job. | Run a CIS, TDA or TD job and use its output. |
| `… has no MO coefficients; use the .fchk of the CIS/TD job.` | The `.fchk` was written without orbitals. | Re-create it with `formchk` from the job's checkpoint. |
| `The configuration coefficients sum to … instead of 0.5 … Rerun Gaussian with IOp(9/40=4)` | Gaussian printed only the coefficients above 0.1 (no `IOp(9/40=4)`). | Rerun the TD job with `IOp(9/40=4)`; until then the NTOs are approximate. |
| `TDHF/TD-DFT state: the de-excitation coefficients … are not part of the NTOs …` | Normal for full TD: the weights can sum slightly above 1. | Nothing; well above 1 suggests using TDA. |
| `No NTO pair dominates …` | The excitation genuinely mixes several hole → particle pairs. | Discuss the listed pairs together. |
| `A grid at … bohr would have … million points. Use a grid spacing of at least … bohr.` | The grid spacing is too fine for the molecule. | Use the suggested spacing (or larger). |
| `No +0.02 surface for the hole.` (status line) | The orbital never reaches the isovalue with that sign. | Lower the isovalue. |
| `The excitation file gives no oscillator strengths …` (Plot Spectrum) | A plain-text excitation file was loaded; it has no f values. | Load the Gaussian output instead. |
| `The peak half-width must be a positive …` (Plot Spectrum) | The half-width field holds zero, a negative number or text. | Type a positive value in eV, or leave it blank for 0.333 eV. |
| `The plot viewer could not be started. The spectrum was saved as …` | The interactive viewer could not open (see [Plot windows do not appear](#plot-windows-do-not-appear)). | The spectrum is saved; reopen it with **Plotter** (*JSON plot file*). |

### Fukui input builder messages

| Message | Cause | What to do |
|---|---|---|
| `%Mem = '2000': a number without a unit means 8-byte words` | A memory value without a unit. | Type e.g. `2GB` or `500MB`. |
| `Keyword Opt: Optimization keywords go in the 'Optimization (job 1)' field …` | `Opt` was typed in *Extra keywords*, which apply to all three jobs. | Put it (with its options) in *Optimization (job 1)*. |
| `The N+1 state: enter its spin multiplicity; it is not assumed (51 electrons allow 2, 4, 6, ...)` | A multiplicity field is blank. | Enter the multiplicity of the lowest-energy spin state of that state. |
| `The N+1 state: multiplicity 1 is impossible (…)` | The multiplicity does not fit the electron count (even count: 1, 3, …; odd: 2, 4, …). | Correct it. |
| `<name>_Np1.chk not found: did that job run?` (from `<name>_formchk.sh`) | Gaussian stopped before that linked job (an earlier job failed, e.g. the optimization did not converge). | Read the end of the `.log` file, fix the cause and rerun. |
| `formchk not found: set up Gaussian 16 first` | The shell has no Gaussian environment. | `module load gaussian` on a cluster, or source `$g16root/g16/bsd/g16.profile`. |
| `Error parsing percent card` (in the Gaussian `.log`) | The input was run with Gaussian 09 Rev. A.02, which does not know `%OldChk`. | Run it with Gaussian 16 (or a later Gaussian 09 revision). |
| `$'\r': command not found` when running the `.sh` | The script was edited and saved with Windows line endings. | Save the input again from the builder, or run `dos2unix <name>_formchk.sh`. |

### Packaged Windows executable cannot find GLFW

When packaging with `auto-py-to-exe`, include the GLFW DLL manually if needed. It is usually located at a path similar to:

```text
venv\Lib\site-packages\glfw\glfw3.dll
```

## Version, Contributors, and License

Version shown by the application:

```text
gQTEA-0.9.5 Molecular Analysis Toolkit
```

Version 0.9.5 (October 2026) adds the [QTAIM Topology Analysis](#qtaim-topology-analysis) tab to the 3D Molecular Viewer:

- Bader's topology of the electron density of a Gaussian `.fchk`: all critical points (nuclear, bond, ring and cage) with ρ, ∇²ρ, the Hessian eigenvalues, ellipticity, η and the energy densities G, K, V and H; bond and ring paths; the Poincaré–Hopf check; and the `<name>_QTAIM.txt` report;
- Multiwfn's algorithms and default parameters, every one adjustable, validated against Multiwfn 2026.7.15 on six molecules;
- the critical points and paths drawn in the 3D window, with a canvas background, a new **Design style** (the Molecular Design Canvas look, also on the Display tab) and independent element-symbol and atom-number labels;
- a Multiwfn-style 2D plane map (density contours, gradient lines, interbasin and bond paths) in the interactive plot viewer;
- image export of the 3D view and of the 2D map;
- `numba` kernels for the density derivatives (4–46× faster analyses);
- a step-by-step procedure, [Example 9](#example-9-qtaim-analysis-of-a-hydrogen-bond-water-dimer) and a [theoretical background](#theoretical-background-the-quantum-theory-of-atoms-in-molecules) in this manual.

Version 0.9.0 (October 2026) turns the [Molecular Design](#molecular-design) canvas into a structure editor:

- [import](#importing-and-saving-structure-files) XYZ, Gaussian input and `.fchk`, MOL/SDF, MOL2 and PDB files, and save in all of them but `.fchk` (bond orders kept where the format allows);
- [select and move one molecule](#selecting-and-moving-one-molecule) (Shift+click, Ctrl+drag, Ctrl+right-drag);
- [measure](#measuring-bond-lengths-angles-and-dihedrals) and [set](#setting-a-bond-length-angle-or-dihedral) bond lengths, angles and dihedrals;
- [freeze](#freezing-parameters-during-pre-optimization) them during pre-optimization, held exactly at every step;
- [add hydrogens](#adding-hydrogens-automatically), [replace an atom's element](#replacing-an-atoms-element), [show atom numbers](#showing-atom-numbers) and [change the background](#canvas-background);
- [bring a molecule back from the 3D viewer](#bringing-a-molecule-back-from-the-3d-viewer) without overwriting a drawing.

Dragging, rotating and panning on the canvas are also about three times faster (five times for larger molecules). `numba` is now installed by default with `installation_requirements.txt`; it stays optional at run time.

Version 0.8.5 (October 2026) added the [Gaussian Fukui Input Builder](#gaussian-fukui-input-builder) on the Inputs tab. It writes the Gaussian 16 job that the Fukui analysis needs:

- one input file with three linked jobs: an optimization of the N state, then single points on N+1 and N−1 that read the optimized geometry from its checkpoint (`%OldChk` + `Geom=Check`), each into its own `.chk`;
- charges set from the electron count; all three spin multiplicities must be entered, none is assumed;
- a `formchk` script that converts the checkpoints to the `.fchk` files the analysis tool reads.

The files use relative names and Linux line endings, for Gaussian 16 under WSL or on a cluster. In *Analysis of Fukui Functions*, choosing the builder's `<name>_N.fchk` now fills in the N+1 and N−1 files. The manual gains [Example 8](#example-8-read-a-fukui-function-figure-phenol-f), which explains how to read an f⁺ figure of phenol.

Version 0.8.0 (October 2026) added the [Analysis of Fukui Functions](#analysis-of-fukui-functions) tool. From the N, N+1 and N−1 `.fchk` files it computes f⁺, f⁻, f⁰ and the dual descriptor, the global indices and built-in Hirshfeld condensed indices, with 3D views, atom styles and numbering, and a compiled density kernel when `numba` is installed.

Version 0.7.5 (October 2026) added the Electrostatic Surface Potential (ESP) tool and the solvation-box improvements: faster packing, optional periodic pre-optimization, system properties (atoms, electrons, spin multiplicity) and a live construction view.

Core development team:

- Ademir J. Camargo - ajc@ueg.br
- Valter H. C. Silva - fatioleg@ueg.br
- Solemar S. Oliveira - solemar@ueg.br
- Hamilton B. Napolitano - hamilton@ueg.br
- Luciano Ribeiro - lribeiro@ueg.br
- Flavio O. Sanches - flavio.neto@ifg.edu.br

No license file is currently included in the source tree. Before distributing the program publicly, add a `LICENSE` file so users know how they may use, modify, and redistribute the code.
