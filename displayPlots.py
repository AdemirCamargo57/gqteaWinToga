import tempfile, os, sys, json, subprocess
import matplotlib
matplotlib.use('Agg') # non-interactive backend in the main (Toga) process
import matplotlib.pyplot as plt


def launch_plot_viewer(manifest_path):
    """Spawn plotViewer.py on an existing manifest file. True if it started.

    From source, sys.executable is a Python interpreter and we hand it
    plotViewer.py. In a frozen/packaged build there is no separate interpreter:
    sys.executable is the app's own .exe, so we re-launch it with the
    --plot-viewer flag, which the entry point (gqteaWinToga.py) dispatches to
    plotViewer before loading the GUI. Either way the figures are fully
    interactive.

    This is the single implementation of that dispatch; callers that already
    hold a manifest on disk (e.g. plotter.py's JSON plot type) use it directly.
    """
    try:
        if getattr(sys, "frozen", False):
            cmd = [sys.executable, "--plot-viewer", manifest_path]
        else:
            viewer = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "plotViewer.py"
            )
            cmd = [sys.executable, viewer, manifest_path]
        subprocess.Popen(cmd)
        return True
    except Exception:
        return False


def series_style_kwargs(curve):
    """plt.plot keyword args for a manifest series' optional styling.

    'linestyle' and 'marker' are optional keys of a "series" curve (written by
    plotter.py's table plot); a curve without them draws exactly as before.
    """
    kwargs = {}
    if "linestyle" in curve:
        kwargs["linestyle"] = curve["linestyle"]
    if "marker" in curve:
        kwargs["marker"] = curve["marker"]
    return kwargs


def draw_bar_comparison(axes, figure_entry, font_style=None):
    """Draw one "bars" manifest entry onto `axes`; return the twin axis or None.

    Public (no leading underscore) because plotViewer.py imports it: it is part
    of this module's interface, not an internal detail. Kept out of the class so
    the standalone viewer can reuse the exact same
    rendering: a figure must not look different depending on which process drew
    it. The percentage line sits on a twin right-hand axis because it shares no
    scale with the absolute values in angstroms or degrees.
    """
    font_style = font_style or {'color': 'darkred', 'weight': 'normal', 'size': 14}
    categories = figure_entry.get("categories", [])
    groups = figure_entry.get("groups", [])
    positions = range(len(categories))
    count = max(len(groups), 1)
    width = 0.8 / count

    handles = []
    for index, group in enumerate(groups):
        offset = (index - (count - 1) / 2.0) * width
        bars = axes.bar(
            [p + offset for p in positions],
            group.get("values", []),
            width=width,
            yerr=group.get("errors"),
            capsize=4,
            label=group.get("label", ""),
        )
        handles.append(bars)

    axes.set_xticks(list(positions))
    axes.set_xticklabels(categories, rotation=60, ha="right", fontsize=8)
    axes.set_xlabel(figure_entry.get("xlabel", ""), fontdict=font_style)
    axes.set_ylabel(figure_entry.get("ylabel", ""), fontdict=font_style)
    axes.set_title(figure_entry.get("title", ""), fontdict=font_style)
    axes.grid(True, axis="y", alpha=0.3)

    twin = None
    line = figure_entry.get("line")
    if line is not None:
        twin = axes.twinx()
        drawn = twin.plot(
            list(positions), line.get("values", []),
            color="black", marker="o", markersize=4, linewidth=1.5,
            label=line.get("label", ""),
        )
        twin.axhline(0.0, color="grey", linestyle="--", linewidth=0.8)
        twin.set_ylabel(line.get("ylabel", ""), fontdict=font_style)
        handles.extend(drawn)

    # One legend for everything: two y axes would otherwise produce two.
    axes.legend(handles=handles, loc="best", fontsize=9)
    return twin


# Multiwfn's plane-map colours (manual 3.5.5): brown (3,-3), blue (3,-1),
# orange (3,+1), green (3,+3); bold dark-brown bond paths, deep-blue interbasin paths.
PLANE_MAP_POINT_COLOURS = {"(3,-3)": "#8B4513", "(3,-1)": "blue", "(3,+1)": "orange", "(3,+3)": "green"}
PLANE_MAP_LINE_STYLES = {
    "gradient": {"color": "0.55", "linewidth": 0.6},
    "interbasin": {"color": "#00008B", "linewidth": 1.4},
    "bond_path": {"color": "#5C2E0B", "linewidth": 2.4},
    "ring_path": {"color": "#2E8B57", "linewidth": 1.6},
}


def draw_plane_map(axes, figure_entry, font_style=None):
    """Draw one "plane_map" manifest entry (the QTAIM 2D map) onto `axes`.

    Public for the same reason as draw_bar_comparison: plotViewer.py and the
    static fallback must draw the identical figure. Contours of the density
    grid (levels as given), then the lines (gradient, interbasin, bond and
    ring paths), then the critical points with their labels; equal aspect.
    """
    font_style = font_style or {'color': 'darkred', 'weight': 'normal', 'size': 14}
    xmin, xmax, ymin, ymax = figure_entry.get("extent", (0.0, 1.0, 0.0, 1.0))
    grid = figure_entry.get("grid")
    levels = figure_entry.get("contours") or []
    if grid is not None and levels:
        rows, cols = len(grid), len(grid[0])
        xs = [xmin + (xmax - xmin) * i / (cols - 1) for i in range(cols)]
        ys = [ymin + (ymax - ymin) * j / (rows - 1) for j in range(rows)]
        axes.contour(xs, ys, grid, levels=sorted(levels), colors="black", linewidths=0.6)
    for line in figure_entry.get("lines", []):
        pts = line.get("points", [])
        if len(pts) >= 2:
            axes.plot([p[0] for p in pts], [p[1] for p in pts],
                      **PLANE_MAP_LINE_STYLES.get(line.get("kind"), {"color": "grey", "linewidth": 0.8}))
    for kind, colour in PLANE_MAP_POINT_COLOURS.items():
        chosen = [p for p in figure_entry.get("points", []) if p.get("type") == kind]
        if not chosen:
            continue
        axes.scatter([p["x"] for p in chosen], [p["y"] for p in chosen], s=36 if kind == "(3,-3)" else 22,
                     color=colour, edgecolors="black", linewidths=0.4, zorder=5, label=kind)
        for p in chosen:
            if p.get("label"):
                axes.annotate(p["label"], (p["x"], p["y"]), xytext=(4, 4), textcoords="offset points",
                              fontsize=8, color=colour, zorder=6,
                              bbox={"boxstyle": "round,pad=0.1", "fc": "white", "ec": "none", "alpha": 0.7})
    axes.set_xlim(xmin, xmax)
    axes.set_ylim(ymin, ymax)
    axes.set_aspect("equal")
    axes.set_xlabel(figure_entry.get("xlabel", ""), fontdict=font_style)
    axes.set_ylabel(figure_entry.get("ylabel", ""), fontdict=font_style)
    axes.set_title(figure_entry.get("title", ""), fontdict=font_style)
    if figure_entry.get("points"):
        axes.legend(loc="upper right", fontsize=8)


PLANE_MAP_IMAGE_SUFFIXES = (".png", ".svg", ".pdf", ".tif", ".tiff")


def save_plane_map_image(figure_entry, path, dpi=200):
    """Write a "plane_map" entry to an image file; returns the path written.

    Drawn by draw_plane_map on an off-screen Agg figure (no pyplot state), so
    the file matches the interactive viewer. No suffix means PNG; other
    suffixes than PLANE_MAP_IMAGE_SUFFIXES are refused with a ValueError.
    """
    import os
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    root, suffix = os.path.splitext(str(path))
    if not suffix:
        path, suffix = root + ".png", ".png"
    if suffix.lower() not in PLANE_MAP_IMAGE_SUFFIXES:
        raise ValueError(f"Save the map as {', '.join(s[1:] for s in PLANE_MAP_IMAGE_SUFFIXES)}, "
                         f"not '{suffix}'.")
    figure = Figure(figsize=(9, 8))
    FigureCanvasAgg(figure)
    draw_plane_map(figure.add_subplot(), figure_entry)
    figure.tight_layout()
    figure.savefig(str(path), dpi=dpi)
    return str(path)


# GaussView's Electronic Spectra window: ivory plot area, dotted grid, black
# curve and sticks, bold black titles, a dark-red frame round the figure.
UV_VIS_BACKGROUND = "#FFFFF0"
UV_VIS_FRAME = "#7F0000"
UV_VIS_FONT = {"color": "black", "weight": "bold", "size": 14}


def draw_uv_vis(axes, figure_entry, font_style=None):
    """Draw one "uv_vis" manifest entry (uvVisSpectrum.spectrum_entry); return the f axis.

    The broadened molar absorptivity on the left axis and one stick per excited
    state, of height f, on a twin right-hand axis, each stick labelled with the
    state's number and label - the layout of GaussView's UV-Vis spectrum.
    ``font_style`` is accepted for the common renderer signature; the GaussView
    look uses its own bold black titles.
    """
    curve = figure_entry.get("curve", {})
    states = figure_entry.get("states", [])
    axes.set_facecolor(UV_VIS_BACKGROUND)
    axes.figure.patch.set_edgecolor(UV_VIS_FRAME)
    axes.figure.patch.set_linewidth(3)
    axes.plot(curve.get("wavelength_nm", []), curve.get("epsilon", []), color="black", linewidth=1.2)
    top = max(curve.get("epsilon", []) or [0.0])
    axes.set_ylim(0.0, top * 1.1 if top > 0 else 1.0)
    if figure_entry.get("xlim"):
        axes.set_xlim(*figure_entry["xlim"])
    axes.grid(True, linestyle=":", color="0.6")
    axes.set_xlabel(figure_entry.get("xlabel", ""), fontdict=UV_VIS_FONT)
    axes.set_ylabel(figure_entry.get("ylabel", ""), fontdict=UV_VIS_FONT)
    axes.set_title(figure_entry.get("title", ""), fontdict=UV_VIS_FONT)

    twin = axes.twinx()
    wavelengths = [s["wavelength_nm"] for s in states]
    strengths = [s["oscillator"] for s in states]
    twin.vlines(wavelengths, 0.0, strengths, colors="black", linewidth=1.2)
    f_top = max(strengths or [0.0])
    twin.set_ylim(0.0, f_top * 1.35 if f_top > 0 else 0.1)
    for s in states:
        twin.annotate(f"{s['number']}: {s.get('label', '')}", (s["wavelength_nm"], s["oscillator"]),
                      xytext=(0, 3), textcoords="offset points", rotation=90, ha="center", va="bottom",
                      fontsize=8)
    twin.set_ylabel(figure_entry.get("ylabel_right", ""), fontdict=UV_VIS_FONT)
    return twin


HEAT_MAP_ROTATE_LABELS = 12   # more column labels than this are drawn vertically


def draw_heat_map(axes, figure_entry, font_style=None):
    """Draw one "heat_map" manifest entry (the Excitation tab's composition / CT matrices).

    ``matrix`` is a list of rows; ``row_labels`` / ``col_labels`` label them;
    ``vmin`` / ``vmax`` (optional, null = automatic) fix the colour scale; the
    colour bar carries ``colorbar``. Returns the image.
    """
    import numpy as np

    font_style = font_style or {}
    matrix = np.asarray(figure_entry.get("matrix", [[]]), dtype=float)
    image = axes.imshow(matrix, aspect="auto", cmap="viridis", interpolation="nearest",
                        vmin=figure_entry.get("vmin"), vmax=figure_entry.get("vmax"))
    rows, cols = figure_entry.get("row_labels", []), figure_entry.get("col_labels", [])
    axes.set_yticks(range(len(rows)))
    axes.set_yticklabels(rows)
    axes.set_xticks(range(len(cols)))
    axes.set_xticklabels(cols, rotation=90 if len(cols) > HEAT_MAP_ROTATE_LABELS else 0)
    axes.set_xlabel(figure_entry.get("xlabel", ""), fontdict=font_style)
    axes.set_ylabel(figure_entry.get("ylabel", ""), fontdict=font_style)
    axes.set_title(figure_entry.get("title", ""), fontdict=font_style)
    axes.figure.colorbar(image, ax=axes, label=figure_entry.get("colorbar", ""))
    return image


ENERGY_LEVEL_WIDTH = 0.6        # width of one column's levels (x units; columns are 1 apart)
ENERGY_LEVEL_DEGENERACY = 0.003  # eV; closer levels are drawn side by side (molecularOrbitals.DEGENERACY_EV)
ENERGY_LEVEL_COLOURS = {"occupied": "#1F3A93", "singly": "#D35400", "virtual": "#7F7F7F"}


def _degenerate_groups(levels, tolerance=ENERGY_LEVEL_DEGENERACY):
    """Levels sorted by energy, split where neighbours differ by more than ``tolerance`` eV."""
    groups = []
    for level in sorted(levels, key=lambda l: l["energy_ev"]):
        if groups and level["energy_ev"] - groups[-1][-1]["energy_ev"] <= tolerance:
            groups[-1].append(level)
        else:
            groups.append([level])
    return groups


def draw_energy_levels(axes, figure_entry, font_style=None):
    """Draw one "energy_levels" manifest entry (molecularOrbitals.energy_levels_entry).

    One column per MO set (alpha and beta side by side for an unrestricted
    file). Occupied levels are solid (singly occupied ones orange), virtual
    levels dashed grey; degenerate levels share their slot side by side, as in
    GaussView's MO editor. Each level carries its MO number; HOMO and LUMO are
    named and the gap is marked with an arrow. ``ylim`` opens the view on the
    frontier orbitals; the rest is reached by zooming or panning.
    """
    font_style = font_style or {}
    columns = figure_entry.get("columns", [])
    half = ENERGY_LEVEL_WIDTH / 2.0
    # Occupancy 1 is "singly occupied" only in a one-set (restricted-open)
    # diagram; in an alpha or beta column it is a full occupation.
    one_set = len(columns) == 1
    for c, column in enumerate(columns):
        for group in _degenerate_groups(column.get("levels", [])):
            width = ENERGY_LEVEL_WIDTH / len(group)
            for k, level in enumerate(group):
                x0 = c - half + k * width + 0.04 * width
                x1 = c - half + (k + 1) * width - 0.04 * width
                occupancy = level.get("occupancy", 0)
                if occupancy <= 0:
                    kind = "virtual"
                else:
                    kind = "singly" if one_set and occupancy < 2 else "occupied"
                axes.hlines(level["energy_ev"], x0, x1, colors=ENERGY_LEVEL_COLOURS[kind], linewidth=2.0,
                            linestyles="solid" if occupancy > 0 else "dashed")
            axes.text(c - half - 0.03, group[0]["energy_ev"], ",".join(str(l["number"]) for l in group),
                      ha="right", va="center", fontsize=7, clip_on=True)
        by_number = {l["number"]: l for l in column.get("levels", [])}
        homo, lumo = by_number.get(column.get("homo")), by_number.get(column.get("lumo"))
        for level, name in ((homo, "HOMO"), (lumo, "LUMO")):
            if level is not None:
                axes.text(c + half + 0.03, level["energy_ev"], f"{name} {level['energy_ev']:.3f} eV",
                          ha="left", va="center", fontsize=8, clip_on=True)
        if homo is not None and lumo is not None:
            x = c + half + 0.02
            axes.annotate("", xy=(x, lumo["energy_ev"]), xytext=(x, homo["energy_ev"]),
                          arrowprops={"arrowstyle": "<->", "color": "black", "linewidth": 0.8})
            axes.text(x + 0.03, 0.5 * (homo["energy_ev"] + lumo["energy_ev"]),
                      f"gap {lumo['energy_ev'] - homo['energy_ev']:.3f} eV", ha="left", va="center",
                      fontsize=8, clip_on=True)
    axes.set_xlim(-0.9, len(columns) - 1 + 1.2)
    labels = [column.get("label", "") for column in columns]
    axes.set_xticks(range(len(columns)) if any(labels) else [])
    if any(labels):
        axes.set_xticklabels(labels)
    if figure_entry.get("ylim"):
        axes.set_ylim(*figure_entry["ylim"])
    axes.grid(True, axis="y", linestyle=":", color="0.8")
    axes.set_ylabel(figure_entry.get("ylabel", ""), fontdict=font_style)
    axes.set_title(figure_entry.get("title", ""), fontdict=font_style)


class DisplayPlots():

    font_style = {'color':  'darkred','weight': 'normal','size': 14}
    saved_plot_files = []
    saved_plot_data = []

    def save_plots(self, k, x ,y, plot_xlabel, plot_ylabel, plot_title, xlim=None, ylim=None, save_png=True):
        # Optional xlim/ylim are (low, high) tuples; when given they are applied
        # here and forwarded to the viewer. Pass save_png=False to skip writing a
        # static PNG to disk (used by tools that rely solely on the interactive
        # viewer and do not want left-over image files next to their data); the
        # raw curve data recorded below is what the interactive viewer consumes.
        if save_png:
            temp_filename = tempfile.NamedTemporaryFile(delete=False, suffix=".png", dir=self.output_dir).name
            plt.figure(k)
            plt.plot(x, y,antialiased=True)
            plt.xlabel(plot_xlabel, fontdict = self.font_style)
            plt.ylabel(plot_ylabel, fontdict = self.font_style)
            plt.title(plot_title, fontdict = self.font_style)
            if xlim is not None:
                plt.xlim(*xlim)
            if ylim is not None:
                plt.ylim(*ylim)
            plt.savefig(temp_filename)
            plt.close()
            self.saved_plot_files.append(temp_filename)

        # Raw data for the interactive viewer (convert range/ndarray -> JSON-safe
        # lists of floats).
        entry = {
            "x": [float(v) for v in x],
            "y": [float(v) for v in y],
            "xlabel": plot_xlabel,
            "ylabel": plot_ylabel,
            "title": plot_title,
        }
        if xlim is not None:
            entry["xlim"] = [float(xlim[0]), float(xlim[1])]
        if ylim is not None:
            entry["ylim"] = [float(ylim[0]), float(ylim[1])]
        self.saved_plot_data.append(entry)

    def save_multiseries_plot(self, k, series, plot_xlabel, plot_ylabel, plot_title, save_png=True):
        """Record a figure holding several labelled curves overlaid with a legend.

        `series` is a list of (x, y, label) tuples. Used e.g. to overlay a total
        spectrum with its per-element partial spectra. Pass save_png=False to skip
        writing a static PNG to disk (see save_plots); the raw series data recorded
        below is what the interactive viewer consumes.
        """
        if save_png:
            temp_filename = tempfile.NamedTemporaryFile(delete=False, suffix=".png", dir=self.output_dir).name
            plt.figure(k)
            for x, y, lbl in series:
                plt.plot(x, y, antialiased=True, label=lbl)
            plt.xlabel(plot_xlabel, fontdict=self.font_style)
            plt.ylabel(plot_ylabel, fontdict=self.font_style)
            plt.title(plot_title, fontdict=self.font_style)
            plt.legend()
            plt.savefig(temp_filename)
            plt.close()
            self.saved_plot_files.append(temp_filename)

        self.saved_plot_data.append({
            "series": [
                {"x": [float(v) for v in x], "y": [float(v) for v in y], "label": lbl}
                for x, y, lbl in series
            ],
            "xlabel": plot_xlabel,
            "ylabel": plot_ylabel,
            "title": plot_title,
        })

    def save_bar_comparison_plot(self, k, categories, groups, plot_xlabel,
                                 plot_ylabel, plot_title, line=None, save_png=True):
        """Record a grouped bar chart with error bars and an optional overlaid line.

        `categories` names the bar positions; `groups` is a list of
        (label, values, errors) tuples, one bar per category per group, with
        `errors` either a sequence of symmetric error-bar half-lengths or None.
        `line` is an optional (label, values, y2label) tuple drawn on a twin
        right-hand axis -- used to put a relative measure (a percentage) beside
        absolute values that share no scale with it.

        Pass save_png=False to skip writing a static PNG (see save_plots); the
        raw data recorded below is what the interactive viewer consumes.
        """
        entry = {
            "type": "bars",
            "categories": [str(c) for c in categories],
            "groups": [],
            "xlabel": plot_xlabel,
            "ylabel": plot_ylabel,
            "title": plot_title,
        }
        for label, values, errors in groups:
            group = {"label": label, "values": [float(v) for v in values]}
            if errors is not None:
                group["errors"] = [float(e) for e in errors]
            entry["groups"].append(group)
        if line is not None:
            line_label, line_values, line_ylabel = line
            entry["line"] = {
                "label": line_label,
                "values": [float(v) for v in line_values],
                "ylabel": line_ylabel,
            }

        if save_png:
            temp_filename = tempfile.NamedTemporaryFile(delete=False, suffix=".png", dir=self.output_dir).name
            figure = plt.figure(k)
            axes = figure.gca()
            draw_bar_comparison(axes, entry)
            plt.tight_layout()
            plt.savefig(temp_filename)
            plt.close(figure)
            self.saved_plot_files.append(temp_filename)

        self.saved_plot_data.append(entry)

    def display_plots(self):
        # Preferred path: launch an interactive matplotlib viewer in a separate
        # process (zoom/pan/save toolbar + live cursor coordinates) so it does
        # not block the Toga event loop.
        #
        # From source, sys.executable is a Python interpreter and we hand it
        # plotViewer.py. In a frozen/packaged build there is no separate
        # interpreter: sys.executable is the app's own .exe, so we re-launch it
        # with the --plot-viewer flag, which the entry point (gqteaWinToga.py)
        # dispatches to plotViewer before loading the GUI. Either way the figures
        # are fully interactive; the static-PNG path below is only a fallback for
        # when the viewer process cannot be spawned at all.
        if self.saved_plot_data:
            try:
                manifest = tempfile.NamedTemporaryFile(
                    delete=False, suffix=".json", dir=self.output_dir, mode="w"
                )
                json.dump(self.saved_plot_data, manifest)
                manifest.close()
                if not launch_plot_viewer(manifest.name):
                    raise RuntimeError("viewer process could not be started")
                self.saved_plot_data = []
                self.saved_plot_files = []
                return
            except Exception:
                pass  # fall back to static image windows below

        self._display_static()

    def _display_static(self):
        """Fallback: show each saved PNG in a Toga image window (non-interactive)."""
        # Deferred here (not a module-level import): plotViewer.py, which runs
        # as the --plot-viewer child process, imports this module for
        # draw_bar_comparison and must never load the Toga/OpenGL stack --
        # this is the only method in the module that touches Toga.
        import toga
        from toga.style import Pack

        for plot_filename in self.saved_plot_files:

            # Load the image using Toga's Image class
            plot_image = toga.Image(plot_filename)

            # Create a new window for the plot
            plot_window = toga.Window(
                title=" ",
                size=(700, 600),
            )

            # Create a box to hold the image
            plot_box = toga.Box(style=Pack(flex=1))
            plot_window.content = plot_box

            # Create an ImageView to display the image
            plot_imageview = toga.ImageView(plot_image, style=Pack(flex=1))
            plot_box.add(plot_imageview)

            # Show the window
            plot_window.show()

        # Clear the lists of saved plots after displaying them
        self.saved_plot_files = []
        self.saved_plot_data = []
