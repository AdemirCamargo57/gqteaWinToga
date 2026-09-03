# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

**gqteaWinToga** (title: *gQTEA Molecular Analysis Toolkit*, current app version 0.4.0) is a Windows desktop GUI for preparing, converting, and analyzing molecular dynamics and electronic-structure simulation data (CPMD, Quantum ESPRESSO `cp.x`, ORCA, surface hopping, trajectory analysis). It is built with the [Toga](https://toga.readthedocs.io/) GUI framework on the WinForms backend.

The git repository lives in `venv/src/` (this directory), which is nested inside a Python `venv` created against Python 3.13. Run all commands from `venv/src/`.

## Running & Development

There is no build step or linter configured — this is a plain Python source tree launched directly.

```powershell
# From venv/src/ with the venv activated (..\Scripts\Activate.ps1)
python gqteaWinToga.py
```

There is a small `pytest` suite under [tests/](tests/) (currently covering the RDF tool). Run it headless from `venv/src/`:

```powershell
python -m pytest tests/ -q
```

Tests drive the analysis classes directly with stubbed Toga widgets/dialogs (no GUI is created); [tests/conftest.py](tests/conftest.py) puts the flat source dir on `sys.path` and pins the `Agg` matplotlib backend. New tests should follow the same pattern — instantiate the `*Analyser` logic class, not the `*UI` class.

Runtime dependencies (install into the venv if missing): `toga-winforms`, `numpy`, `scipy`, `matplotlib`, `glfw`, `PyOpenGL`, `PyOpenGL-accelerate`, `PyMuPDF`. A graphical desktop session is required. The 3D molecular viewer needs `glfw3.dll` from `venv/Lib/site-packages/glfw/`.

Packaging to a Windows `.exe` is done manually with `auto-py-to-exe` (bundle the deps above; the GLFW DLL often must be added by hand). See `readme.md` for the full packaging notes.

## Architecture

**Flat module tree, one tool per file.** [gqteaWinToga.py](gqteaWinToga.py) is the only entry point. Its `gqteaWin(toga.App).startup()` imports one `*UI` class per tool and registers them as buttons across five `OptionContainer` tabs: **Geometry**, **Inputs** (`cpmd_tools`), **Structural**, **Thermo**, **Tools** (`general_tools`). Each button's `on_press` is the `*UI` class itself — clicking instantiates the class, which opens its own `toga.Window`.

**Standard per-tool pattern** (see [bond.py](bond.py) as the canonical example):
- A logic/analysis class (e.g. `BondAnalyser`) holding the numeric work and file I/O.
- A `*UI` subclass (e.g. `BondUI`) whose `__init__(self, *args)` calls `layout_main_window(widget)` to build the Toga window, then defines an `async workflow(self, widget)` that orchestrates: read/validate params → compute → write outputs → optionally show plots. `closeTopLevel` closes the window.
- Analysis classes multiply-inherit shared mixins rather than composing them.

**Shared mixins / modules** (inherited or imported by most tools):
- [framesCounter.py](framesCounter.py) — `FramesCounter` mixin: file-open dialogs, reads an `.xyz`/`TRAJEC.xyz`, sets `self.trajec`, `self.num_atoms`, `self.total_frame_number`, and `self.output_dir` (defaults to the trajectory's directory).
- [displayPlots.py](displayPlots.py) — `DisplayPlots` mixin: `save_plots(k, x, y, xlabel, ylabel, title, xlim=None, ylim=None, save_png=True)` records each curve (JSON-safe x/y data for the viewer, plus a static PNG unless `save_png=False`); `display_plots()` launches [plotViewer.py](plotViewer.py) in a **separate process** to show interactive matplotlib windows (zoom/pan/save toolbar). This now works in **frozen/packaged builds too**: from source it spawns `[python, plotViewer.py, manifest]`; when `sys.frozen`, it re-launches the app's own exe as `[sys.executable, "--plot-viewer", manifest]`, which the [gqteaWinToga.py](gqteaWinToga.py) entry point dispatches to `plotViewer.main` *before* importing the Toga/OpenGL stack. `plotViewer.py` pins the `TkAgg` backend (must be bundled — see readme packaging notes); the static-PNG Toga-image-window path (`_display_static`) remains only as a fallback for when the viewer process cannot be spawned at all. The main process pins the `Agg` backend so it never conflicts with the viewer's GUI loop. Shared by all analysis tools, so plotting changes here are global.
- [plotter.py](plotter.py) — `PlotterBase`/`PlotterUI`, standalone energy-plot tool that parses both CPMD and gqteaMD energy files.
- [help.py](help.py) — `AtomicData` (element `atomic_masses`, `atomic_numbers`, etc.), `Fonts`, and `HelpGqteaWin` (long help strings like `help_bond_analysis` shown in each tool's `MultilineTextInput`). This is the shared reference-data + help-text hub.

**Conventions observed across tools:**
- Async Toga event handlers throughout; user errors surface via `toga.InfoDialog` (often a `warning_function` helper) rather than exceptions.
- `.xyz` format is parsed manually: line 1 = atom count, line 2 = comment, then `num_atoms` lines of `element x y z`. Atom labels in the UI are **1-based**.
- Outputs (`.dat`, `.csv`, summary `.txt`, plots) are written next to the input file in `self.output_dir`. Time is derived from the atomic time unit constant `atufs = 0.02418884326505` fs.
- Files may contain large commented-out "previous version" blocks kept for reference (e.g. the bottom half of `bond.py`); the active code is the uncommented top.
- Two distinct trajectory formats are in play: CPMD **`TRAJEC.xyz`** (standard `.xyz`, positions in Å) and CPMD **`TRAJECTORY`** (no per-frame header; each line is `step x y z vx vy vz` in atomic units, positions *and* velocities). [autocorrelationFunction.py](autocorrelationFunction.py) reads the latter and needs `GEOMETRY.xyz` loaded first for the atom count. Note the `TRAJECTORY` file can contain CPMD **restart markers** (`<<<<<<  NEW DATA  >>>>>>`) and blank lines mid-file when a run is continued; any code reading it must skip non-data lines to keep atom framing aligned (see `_next_atom_tokens`).

## Interactive figures without PNG clutter (reusable recipe)

The pattern first applied in [radialDistribution.py](radialDistribution.py): show plots **only** through the interactive viewer (zoom/pan/save toolbar) — working in both source and packaged/frozen builds — and stop leaving static PNG files next to the user's data. The plumbing is already global in [displayPlots.py](displayPlots.py), [plotViewer.py](plotViewer.py), and the [gqteaWinToga.py](gqteaWinToga.py) entry point, so **applying it to another tool is a one-line change per plot call**:

1. In the tool, pass `save_png=False` to every `save_plots(...)` call (and keep the single `self.display_plots()` after them). That is the whole change — the tool already benefits from the frozen-build viewer via the shared mixin.
2. No change needed to `display_plots()`: from source it spawns `[python, plotViewer.py, manifest]`; when `sys.frozen` it re-launches the app's own exe as `[sys.executable, "--plot-viewer", manifest]`. The [gqteaWinToga.py](gqteaWinToga.py) entry point intercepts `--plot-viewer <manifest>` **before** importing the Toga/OpenGL stack and dispatches to `plotViewer.main`.
3. Packaging requirement (once, global): the frozen build must bundle `tkinter` / the **`TkAgg`** backend that `plotViewer.py` pins, or the interactive window silently fails to open (see readme packaging notes). If a build shows *no* window rather than a PNG, `TkAgg` is missing from the bundle.
4. `save_multiseries_plot` (overlay figures) takes the same `save_png` flag — pass `save_png=False` there too for tools that overlay curves (e.g. total + partial VDOS in `autocorrelationFunction.py`).
5. `_display_static` (static PNG in Toga image windows) remains only as a fallback for when the viewer subprocess cannot be spawned at all; with `save_png=False` there is no PNG to fall back to, which is the intended trade-off.
6. Tests: assert `glob("*.png") == []` in the output dir after the workflow while the `.dat`/`.csv` output is still written, and (for the shared behavior) that `display_plots()` spawns the viewer under both `sys.frozen = False` and `True`. See [tests/test_radialDistribution.py](tests/test_radialDistribution.py) for the reference tests.

## Tool-specific notes

- **[autocorrelationFunction.py](autocorrelationFunction.py)** — VAF/PAF are averaged over all overlapping time origins, include lag 0, and are normalized to `C(0)=1`; positions are mean-subtracted per atom before the PAF. The correlation uses an **FFT (Wiener–Khinchin) estimator** (`_fft_autocorrelation`, via `scipy.fft` with `workers=-1`): O(N log N) instead of the old O(N·max_lag) lag loop, verified identical to the direct method to ~1e-15 (~22× faster on the compute). `_is_data_line` is a cheap `len(tokens) >= 7` test (restart markers have 4 tokens) — do not re-add float-parsing there; the values are parsed once when stored. Note: a `np.fromstring` bulk loader was tried and rejected (10× slower, ~1 GB) — keep the line-by-line loader. On `TRAJECTORY` load, blank frame fields are auto-filled for the whole trajectory (`_apply_frame_defaults`): start=1, stop=last frame, and *Number of frames for each ACF* = `_default_num_frame_acf()` ≈ 10% of frames capped at 2000. The power spectrum (`calculate_psd`) is the Wiener–Khinchin route on the VAF: mean-removal → mirror about t=0 → window (Welch/Hann/None) → zero-pad → `rfft` → `|·|²`, one-sided, peak-normalized, with a physical frequency axis (cm⁻¹/THz/Hz from the time-step × sampling-interval spacing). Mass-weighting (masses from [help.py](help.py) `AtomicData`) computes a separate mass-weighted VACF (`self.vaf_mw`) in the same correlation loop and emits `VDOS.dat`; the plain spectrum emits `PSD.dat`. PSD UI controls: `switch_mass_weight`, `selection_window`, `textInput_zero_pad`, `selection_freq_unit`, `textInput_partial`. **Partial VDOS** (`textInput_partial`, parsed by `_parse_partial_groups`): `;`-separated groups of element symbols and/or 1-based atom indices/ranges; each group's atoms give a partial VACF (`self.partial_vacfs`) whose spectrum is overlaid on the total (each normalized to its own max) and saved as `PSD_<group>.dat`/`VDOS_<group>.dat` — used for vibrational-mode assignment. Overlays use `DisplayPlots.save_multiseries_plot` + the `series` manifest format in [plotViewer.py](plotViewer.py). All plot calls here pass `save_png=False` (see the *Interactive figures without PNG clutter* recipe), so figures are shown only through the interactive viewer and **no static PNG files are written** next to the trajectory.
- **[radialDistribution.py](radialDistribution.py)** — g(r) around a single **shell-center** atom (`shell_center`, 1-based) to all atoms of a target symbol, from `TRAJEC.xyz`. Periodic boundaries use the vectorized **minimum-image convention** (`delta -= box * round(delta/box)`), which is exact only for `radius ≤ min(a,b,c)/2` — this bound is enforced in `read_params` (which returns `bool`; `workflow` aborts if it fails). Coordinates need **not** be wrapped into the box: min-image tolerates unwrapped input, verified by a box-shift-invariance test. Distances are histogrammed into bins of the requested width, normalized by shell volume and the first-frame ideal density (`self.rho`, from `ideal_density`), and the reported `r` is the **bin center**. `self.axis` = user x/y limits for the coordination-number plot only. Output `RDF_<center><n>_<symbol>.dat` (r, g(r), integral) is written next to the trajectory; both plots go through the shared `DisplayPlots` mixin (the coordination plot passes `xlim`/`ylim` to `save_plots`, an optional/back-compatible extension also honored by [plotViewer.py](plotViewer.py)). This tool passes `save_png=False`, so **no static PNG image files are written** next to the trajectory — the figures are shown only through the interactive viewer (which now also works in frozen builds). The frame loop distinguishes clean EOF/trailing-blank lines from malformed/truncated frames (warns and reports processed-vs-expected). Covered by [tests/test_radialDistribution.py](tests/test_radialDistribution.py) (analytic simple-cubic coordination numbers + validation matrix).
- **[molecularViewer.py](molecularViewer.py)** — GLFW + OpenGL 3D viewer for `.xyz`/`TRAJEC.xyz`, split into `MolecularViewer` (GL/render/trajectory logic) and `MolecularViewerUI` (Toga window). The GL window runs in its own daemon thread (`main_loop`); **GLFW input callbacks fire on that render thread** during `glfw.poll_events`, so picking touches GL state on the same thread as rendering (no cross-thread GL calls). Shared mutable state (frames, `picked_atoms`, measurement overlay) is guarded by `self._state_lock`. Text labels are drawn with GLUT bitmap fonts via `_draw_text_3d` (needs `_glut_ready`; `glutInit` is best-effort in `_ensure_glut_ready`). **Click-to-identify picking**: a left-click *without* dragging (drag threshold >3 px, tracked by `_left_press_pos`/`_left_dragged` so left-drag still rotates) toggles the atom's 1-based index label on the canvas, drawn in `pick_label_color` (cyan, distinct from the yellow measurement text). `main_loop` snapshots `(modelview, projection, viewport)` into `self._pick_view` **inside the molecule push matrix, right before `render_frame`**, so the captured transform matches the space atoms are drawn in (camera + view-rotation + per-molecule rotation). `_pick_atom_at` maps the GLFW top-left cursor to OpenGL's bottom-left window origin and picks the nearest projected atom within `pick_pixel_threshold` (15 px), tie-broken by smallest depth; the projection math lives in `_project_atom`, a **pure numpy `gluProject` equivalent** (matrices in glGet column-major layout, so `vec @ matrix`) that needs no GL context and is unit-tested headless in [tests/test_molecularViewer.py](tests/test_molecularViewer.py). Picked labels persist across frames (atom ordering is constant); the **Clear labels** button (`clear_atom_labels` → `clear_picked_atoms`) empties the set. The live GL capture, GLFW callback wiring, and Toga button require a display and are not covered by the headless tests.

## Adding a new tool

1. Create a self-contained `myTool.py` exposing a `MyToolUI` class following the `layout_main_window` / `workflow` pattern above, reusing the `FramesCounter`/`DisplayPlots` mixins where relevant.
2. Import it in [gqteaWinToga.py](gqteaWinToga.py) and add `("Button label", MyToolUI)` to the appropriate tab list.
3. Add help text to `HelpGqteaWin` in [help.py](help.py) and, if user-facing, document it in [USER_MANUAL.md](USER_MANUAL.md).

## Docs

- [readme.md](readme.md) — capabilities, install, packaging, contributors.
- [USER_MANUAL.md](USER_MANUAL.md) — end-user workflows for every tool.
- `cpx_input_description.pdf` — reference for the `cp.x` input builder.
