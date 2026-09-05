# gqteaWinToga User Manual

## Index

- [Overview](#overview)
- [Installation and Requirements](#installation-and-requirements)
  - [Requirements](#requirements)
  - [Install from Source](#install-from-source)
  - [Start the Program](#start-the-program)
- [Input File Basics](#input-file-basics)
  - [XYZ Files](#xyz-files)
  - [Atom Labels](#atom-labels)
  - [Common Units](#common-units)
- [Basic Workflows](#basic-workflows)
  - [Analyze a Bond Length](#analyze-a-bond-length)
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
  - [Input Builders](#input-builders)
    - [CPMD Inputs](#cpmd-inputs)
    - [Surface Hopping Input Builder](#surface-hopping-input-builder)
    - [Collision Input](#collision-input)
    - [Quantum ESPRESSO cp.x Input Builder](#quantum-espresso-cpx-input-builder)
    - [ORCA Input Builder](#orca-input-builder)
    - [gqteaMD Input Builder](#gqteamd-input-builder)
    - [Vanderbilt runatom.x Input Builder](#vanderbilt-runatomx-input-builder)
  - [Structural Tools](#structural-tools)
    - [Radial Distribution Function](#radial-distribution-function)
    - [Mean Residence Time](#mean-residence-time)
    - [Legacy Mean Residence Time](#legacy-mean-residence-time)
    - [Autocorrelation Function](#autocorrelation-function)
    - [Single Solute Solvent Box](#single-solute-solvent-box)
    - [Mixture of Two Solvent Box](#mixture-of-two-solvent-box)
    - [Shared-Wall Double Solvent Box](#shared-wall-double-solvent-box)
  - [Thermo Tools](#thermo-tools)
    - [Classical Rate Constant](#classical-rate-constant)
  - [General Tools](#general-tools)
    - [3D Molecular Viewer](#3d-molecular-viewer)
    - [Energy Plots](#energy-plots)
    - [Molecular Axis Alignment](#molecular-axis-alignment)
    - [Select Frames](#select-frames)
    - [Frame Selection by Interatomic Distance Range](#frame-selection-by-interatomic-distance-range)
    - [CPMD Input to XYZ Converter](#cpmd-input-to-xyz-converter)
    - [SH Geometry Analyzer](#sh-geometry-analyzer)
    - [Convert cp.x .pos File to trajec.xyz](#convert-cpx-pos-file-to-trajecxyz)
    - [Compute Forces from cp.x .for File](#compute-forces-from-cpx-for-file)
    - [Coordinate Converter](#coordinate-converter)
    - [Convert Coordinates from Angstrom to Bohr](#convert-coordinates-from-angstrom-to-bohr)
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
- [Formulas and Methods (All Bond / Angle / Dihedral tools)](#formulas-and-methods-all-bond--angle--dihedral-tools)
- [Formulas and Methods (Bond Length / Bond Angle / Dihedral Angle tools)](#formulas-and-methods-bond-length--bond-angle--dihedral-angle-tools)
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
  - [Missing `SH_STATE.dat` or `TRAJEC.xyz`](#missing-sh_statedat-or-trajecxyz)
  - [Molecular viewer does not open or shows no 3D view](#molecular-viewer-does-not-open-or-shows-no-3d-view)
  - [Packaged Windows executable cannot find GLFW](#packaged-windows-executable-cannot-find-glfw)
- [Version, Contributors, and License](#version-contributors-and-license)

## Overview

gqteaWinToga is a desktop molecular analysis toolkit developed by the gQTEA group. It helps users prepare, inspect, convert, visualize, and analyze molecular simulation data from CPMD, Quantum ESPRESSO `cp.x`, ORCA, surface hopping workflows, and related molecular dynamics tools.

The program is a graphical application built with Python and Toga. Most tools follow the same pattern: open a tool window from the main tabbed interface, browse to one or more input files, enter calculation parameters, then click the tool's action button to generate outputs.

The current application title identifies the toolkit as:

```text
gQTEA-0.4.0 Molecular Analysis Toolkit
```

## Installation and Requirements

### Requirements

Recommended environment:

- Python 3.10 or newer.
- Windows with a working graphical desktop session.
- A Python virtual environment.

Python packages used by the application include:

- `toga-winforms`
- `numpy`
- `scipy`
- `matplotlib`
- `glfw`
- `PyOpenGL`
- `PyOpenGL-accelerate`
- `PyMuPDF`

Some workflows also require files produced by external chemistry programs, such as CPMD, Quantum ESPRESSO `cp.x`, ORCA, Gaussian, or Vanderbilt `runatom.x`. gqteaWinToga prepares and analyzes files for these programs, but it does not replace the external simulation engines.

### Install from Source

From PowerShell:

```powershell
git clone https://github.com/AdemirCamargo57/gqteaWinToga.git
cd gqteaWinToga
python -m venv venv
.\venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install toga-winforms numpy scipy matplotlib glfw PyOpenGL PyOpenGL-accelerate PyMuPDF
```

If you already have the source tree and an active virtual environment, install the packages in that environment.

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

The newer MRT tool supports:

- XYZ trajectory selection.
- Time step and time-step unit.
- Cutoff radius.
- Tolerance frames.
- Reference definition.
- Observed atom definition.
- Run, export, and clear actions.

Outputs commonly include survival/correlation data, event durations, and text summaries.

#### Legacy Mean Residence Time

The legacy MRT workflow is retained for compatibility with older analyses.

Inputs include:

- Shell inner radius.
- Shell outer radius.
- Simulation time step.
- Sampling interval.
- Tolerance frames.
- Atom labels to exclude.
- Atom labels at the shell center.
- Element symbol to investigate.
- `TRAJEC.xyz` trajectory file.

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

Use **Structural > Single Solute solvent box** to create a box containing one solute and inserted solvent molecules.

Inputs and options:

- Box lattice vectors.
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

Typical outputs:

- `single_solute_solvent.txt`
- `single_solute_box_cmass.xyz`
- A generated solvent-box XYZ file.

If the target number of solvent molecules cannot be inserted, the summary reports a warning and the number of successful insertions.

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

Outputs include a generated XYZ structure and a TXT summary with composition and insertion statistics.

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

Typical outputs:

- `shared_wall_box.xyz`
- `shared_wall_box.txt`

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

**Window layout.** The controls are arranged as a fixed header, four tabs, and a status line:

- **Header** (always visible, whichever tab you are on): the XYZ file field with **Browse**, plus the two main actions — **Display Molecule/Trajectory**, which opens the 3D window, and **Save Current Frame XYZ**.
- **Display** tab: projection, atom style, atom and bond scale, maximum bond length, rotation, and atom labels.
- **Frames** tab: frame navigation and trajectory playback.
- **Measure** tab: distances, angles, dihedrals, and atom coordinates.
- **Box & Performance** tab: periodic box settings and playback performance options.
- **Status line** (bottom): reports the last action, such as trajectory loading progress or the value you just applied.

A typical session is: **Browse** for the file, adjust anything you need on the **Display** tab, then click **Display Molecule/Trajectory**. The 3D window opens alongside the control window, so you can keep changing settings while it is open.

**Numeric fields are applied by pressing Enter.** Type a value into a field — maximum bond length, atom or bond scale, box edge lengths, playback delay, frame step, frame number, or a rotation angle — and press Enter to apply it. There are no separate "Set" buttons.

##### Display tab

- **Projection**: `Orthographic` or `Perspective`.
- **Atom style**: `Line style`, `CPK style`, or `vdW style`.
- **Atom scale** and **Bond scale**: relative size of the drawn atoms and bonds.
- **Max bond length**: upper distance limit, in Å, for two atoms to be drawn as bonded.
- **Rotation**: rotate the displayed system about X, Y, and Z (described below).
- **Atom numbers** and **Atomic symbols**: label every atom in the 3D view.
- **Clear clicked labels**: remove the labels you added by clicking individual atoms.

##### Frames tab

- `|<`, `<<`, `>>`, `>|` and the slider move through the trajectory. The field to the left of the slider jumps to a specific frame when you press Enter; the counter to the right shows `current / last`.
- **Play** / **Pause** animates the trajectory.
- **Loop mode**: `Loop` restarts at the beginning, `Once` stops at the end, `Rock` reverses direction at each end.
- **Frame step**: how many frames each step or playback advance skips. Use `<` and `>`, or type a value and press Enter.
- **Delay (s)**: seconds between frames during playback. Larger values play more slowly.
- **Auto-zoom to fit each frame**: rescales the view for every frame. Turn it off for a steadier view during playback.

##### Measure tab

- Choose **Bond length**, **Bond angle**, **Dihedral angle**, or **Atom coordinates**, enter the atom labels, and click **Measure**. Atom indices are 1-based. The value appears in the **Result** field and on the 3D canvas.

##### Box & Performance tab

- **Box a, b, c**: orthorhombic periodic box edge lengths in Å. Fill in all three, then press Enter.
- **Box center**: where the displayed box is placed (described below).
- **Show simulation box**: draw or hide the box.
- **Reduce overlays during playback**: skips some overlay drawing to keep long trajectories smooth.
- **Bond mode**: how connectivity is updated during playback (described below).

The saved current frame includes the atoms from the displayed frame and a comment noting the frame number.

You can identify individual atoms directly in the 3D view by clicking on them. Left-click an atom (a click without dragging — dragging still rotates the view) and its 1-based index appears next to the atom on the canvas. This is useful for picking out only the solute atoms in a solvated `TRAJEC.xyz` without labeling every solvent atom. Left-click the same atom again to remove its label. The labels stay attached to the same atoms as you step through or play the trajectory, because atom ordering is constant across frames. Click **Clear clicked labels** on the **Display** tab (next to the atom-number and atomic-symbol switches) to remove all clicked labels at once.

The rotation controls on the **Display** tab let you rotate the displayed molecular system around the X, Y, and Z axes without changing the coordinates stored in the loaded structure or trajectory. Enter rotation increments in degrees, then either press Enter in one of the angle fields or click **Apply angles** to add those increments to the current molecular orientation. Click **Reset rotation** to return the displayed molecular orientation to zero rotation. The `<` and `>` buttons beside each axis start continuous rotation in the negative or positive direction; click the same arrow again to stop it, or click another arrow to switch to that axis and direction.

When the periodic box is shown, it rotates with the molecular system. The **Box center** selector on the **Box & Performance** tab controls where the displayed orthorhombic box is placed:

- **Geometric center**: centers the box on the molecular system's geometric center.
- **Bottom at z=0**: centers the box in the a-b plane while placing the bottom face of the box at `z = 0`.

Bond mode controls how the viewer updates connectivity while stepping through or playing a trajectory:

- **Static first frame**: calculates bonds from frame 0 and reuses that same connectivity for all frames. This is useful when you want stable visual connectivity during normal vibrations or rotations.
- **Dynamic cached**: recalculates bonds for each frame and stores the result for faster revisiting of frames. This is useful for trajectories where bonds may form or break and you still want smooth playback.
- **Dynamic live**: recalculates bonds every time the current frame is rendered, without using cached bond lists. This is useful when you are actively changing bond-length settings or want the freshest possible connectivity during inspection.

To use this feature, open **Tools > 3D Molecular Viewer**, load an XYZ trajectory, then choose the desired option from the **Bond mode** selector on the **Box & Performance** tab. The viewer updates the displayed bonds using the selected mode.

Example: if a trajectory shows two atoms separating during a dissociation event, choose **Dynamic cached** or **Dynamic live**. As you step through the frames, the bond disappears when the atom distance exceeds the current upper bond-length limit. If you choose **Static first frame**, that bond remains visible throughout playback because the viewer keeps the first-frame connectivity.

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

Use **Tools > Coordinate converter** to convert between molecular coordinate file formats.

Controls:

- Input format selection.
- Output format selection.
- Browse for input file.
- Choose output file.
- **Convert**.

Supported formats are determined by the converter UI and parser implementation. The converter writes a new file rather than modifying the original.

#### Convert Coordinates from Angstrom to Bohr

Use **Tools > Convert coordinates from A to Bohr** to convert XYZ coordinates from Angstrom to Bohr.

Inputs:

- Input XYZ file.
- Output filename.

Action:

- **Convert**.

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

### Missing `SH_STATE.dat` or `TRAJEC.xyz`

The SH Geometry Analyzer skips subfolders that do not contain both files. Add the missing files or remove incomplete subfolders from the selected root directory.

### Molecular viewer does not open or shows no 3D view

The viewer depends on OpenGL and GLFW. Confirm that `glfw`, `PyOpenGL`, and `PyOpenGL-accelerate` are installed and that your graphics driver supports OpenGL in the current desktop session.

### Packaged Windows executable cannot find GLFW

When packaging with `auto-py-to-exe`, include the GLFW DLL manually if needed. It is usually located at a path similar to:

```text
venv\Lib\site-packages\glfw\glfw3.dll
```

## Version, Contributors, and License

Version shown by the application:

```text
gQTEA-0.4.0 Molecular Analysis Toolkit
```

The launcher source notes that the program was revised and updated in March 2026.

Core development team:

- Ademir J. Camargo - ajc@ueg.br
- Valter H. C. Silva - fatioleg@ueg.br
- Solemar S. Oliveira - solemar@ueg.br
- Hamilton B. Napolitano - hamilton@ueg.br
- Luciano Ribeiro - lribeiro@ueg.br
- Flavio O. Sanches - flavio.neto@ifg.edu.br

No license file is currently included in the source tree. Before distributing the program publicly, add a `LICENSE` file so users know how they may use, modify, and redistribute the code.
