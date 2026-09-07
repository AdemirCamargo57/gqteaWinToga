import toga, tempfile, os, sys, json, subprocess
from toga.style import Pack
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
