"""Standalone interactive plot viewer for gqteaWinToga.

Reads a JSON manifest describing one or more figures and shows them with
matplotlib's default *interactive* backend (zoom, pan, cursor coordinates and
save toolbar). It is launched by ``DisplayPlots.display_plots`` in a separate
process so the interactive Tk/Qt event loop never blocks the Toga GUI.

Manifest format (list of figures)::

    [{"x": [...], "y": [...], "xlabel": "...", "ylabel": "...", "title": "...",
      "xlim": [lo, hi], "ylim": [lo, hi]}, ...]

(``xlim``/``ylim`` are optional.)

Usage::

    python plotViewer.py <manifest.json>
"""
import sys
import json

import matplotlib
# Pin an interactive GUI backend explicitly. In a packaged (frozen) build
# matplotlib's autodetection can fall back to the non-interactive Agg backend,
# which would silently produce no window; TkAgg ships with CPython (tkinter) and
# is the interactive backend already used when running from source.
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

FONT_STYLE = {"color": "darkred", "weight": "normal", "size": 14}


def build_figures(figures):
    """Create one matplotlib figure per manifest entry (no blocking show).

    Each entry is either a single curve ("x"/"y") or several labelled curves
    ("series": [{"x", "y", "label"}, ...]) drawn overlaid with a legend.
    """
    for i, fig in enumerate(figures, 1):
        plt.figure(i)
        if "series" in fig:
            for s in fig["series"]:
                plt.plot(s.get("x", []), s.get("y", []), antialiased=True,
                         label=s.get("label", ""))
            plt.legend()
        else:
            plt.plot(fig.get("x", []), fig.get("y", []), antialiased=True)
        plt.xlabel(fig.get("xlabel", ""), fontdict=FONT_STYLE)
        plt.ylabel(fig.get("ylabel", ""), fontdict=FONT_STYLE)
        plt.title(fig.get("title", ""), fontdict=FONT_STYLE)
        if "xlim" in fig:
            plt.xlim(*fig["xlim"])
        if "ylim" in fig:
            plt.ylim(*fig["ylim"])
        plt.grid(True, alpha=0.3)
        plt.tight_layout()


def main(manifest_path):
    with open(manifest_path, "r") as fh:
        figures = json.load(fh)
    build_figures(figures)
    plt.show()   # blocks in this standalone process until all windows close


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        print("usage: python plotViewer.py <manifest.json>", file=sys.stderr)
