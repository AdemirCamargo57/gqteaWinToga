"""Tests for the shared plotting mixin and the standalone viewer.

Run from ``venv/src/``::

    python -m pytest tests/ -q

The bar-comparison figure is the third manifest shape (after the single curve
and the labelled-series overlay). These tests pin its manifest contract and,
just as importantly, that the two older shapes still render unchanged -- every
analysis tool shares this module.
"""
import glob
import os
import subprocess
import sys

import pytest

from displayPlots import DisplayPlots

# The --plot-viewer child process (gqteaWinToga.py) imports plotViewer, which
# imports displayPlots for draw_bar_comparison, *before* the Toga/OpenGL stack
# is loaded -- that ordering is the whole point of the --plot-viewer dispatch,
# so a frozen build can drop Toga from the viewer's half of the app. Both
# modules must therefore be importable without ever pulling in toga. Checked
# in a subprocess with a clean interpreter, since this pytest session may
# already have imported toga (e.g. via another test module).
SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _import_without_toga(module_name):
    result = subprocess.run(
        [sys.executable, "-c",
         f"import {module_name}, sys\n"
         "print('TOGA_ABSENT' if 'toga' not in sys.modules else 'TOGA_PRESENT')"],
        cwd=SRC_DIR,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def test_importing_displayplots_does_not_load_toga():
    assert _import_without_toga("displayPlots") == "TOGA_ABSENT"


def test_importing_plotviewer_does_not_load_toga():
    assert _import_without_toga("plotViewer") == "TOGA_ABSENT"


def fresh(tmp_path):
    """A DisplayPlots with per-instance buffers.

    saved_plot_files/saved_plot_data are CLASS attributes on the mixin, so a
    test that did not shadow them would see another test's figures.
    """
    plotter = DisplayPlots()
    plotter.output_dir = str(tmp_path)
    plotter.saved_plot_files = []
    plotter.saved_plot_data = []
    return plotter


CATEGORIES = ["C1-C2", "C2-O3"]
GROUPS = [
    ("isolated", [1.52, 1.43], [0.01, 0.02]),
    ("solvated", [1.54, 1.42], [0.02, 0.01]),
]


def test_bar_plot_records_a_typed_manifest_entry(tmp_path):
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES, GROUPS, "parameter", "bond distance (A)", "t",
        save_png=False,
    )

    entry = plotter.saved_plot_data[0]
    assert entry["type"] == "bars"
    assert entry["categories"] == CATEGORIES
    assert [group["label"] for group in entry["groups"]] == ["isolated", "solvated"]
    assert entry["groups"][0]["values"] == [1.52, 1.43]
    assert entry["groups"][0]["errors"] == [0.01, 0.02]
    assert "line" not in entry


def test_bar_plot_records_the_overlaid_line_when_given(tmp_path):
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES, GROUPS, "parameter", "bond distance (A)", "t",
        line=("%dr", [1.25, -0.84], "percent difference (%)"),
        save_png=False,
    )

    line = plotter.saved_plot_data[0]["line"]
    assert line["label"] == "%dr"
    assert line["values"] == [1.25, -0.84]
    assert line["ylabel"] == "percent difference (%)"


def test_bar_plot_omits_errors_that_were_not_supplied(tmp_path):
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES, [("isolated", [1.52, 1.43], None)],
        "parameter", "d (A)", "t", save_png=False,
    )

    assert "errors" not in plotter.saved_plot_data[0]["groups"][0]


def test_bar_plot_writes_no_png_when_asked_not_to(tmp_path):
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES, GROUPS, "parameter", "d (A)", "t", save_png=False,
    )

    assert plotter.saved_plot_files == []
    assert glob.glob(os.path.join(str(tmp_path), "*.png")) == []


def test_bar_plot_can_still_write_a_static_png(tmp_path):
    """The _display_static fallback needs a real image when one is asked for."""
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES, GROUPS, "parameter", "d (A)", "t",
        line=("%dr", [1.25, -0.84], "percent (%)"), save_png=True,
    )

    assert len(plotter.saved_plot_files) == 1
    assert len(glob.glob(os.path.join(str(tmp_path), "*.png"))) == 1


def test_values_are_json_safe_floats(tmp_path):
    """numpy scalars must not reach json.dump."""
    import json
    np = pytest.importorskip("numpy")
    plotter = fresh(tmp_path)

    plotter.save_bar_comparison_plot(
        1, CATEGORIES,
        [("isolated", np.array([1.52, 1.43]), np.array([0.01, 0.02]))],
        "parameter", "d (A)", "t",
        line=("%dr", np.array([1.25, -0.84]), "percent (%)"),
        save_png=False,
    )

    json.dumps(plotter.saved_plot_data)  # raises TypeError on numpy scalars


# --------------------------------------------------------------------------- #
# The viewer renders all three manifest shapes                                  #
# --------------------------------------------------------------------------- #
def test_viewer_renders_a_bar_figure_with_a_twin_axis():
    import matplotlib.pyplot as plt
    import plotViewer

    plt.close("all")
    plotViewer.build_figures([{
        "type": "bars",
        "categories": CATEGORIES,
        "groups": [
            {"label": "isolated", "values": [1.52, 1.43], "errors": [0.01, 0.02]},
            {"label": "solvated", "values": [1.54, 1.42], "errors": [0.02, 0.01]},
        ],
        "line": {"label": "%dr", "values": [1.25, -0.84],
                 "ylabel": "percent difference (%)"},
        "xlabel": "parameter", "ylabel": "bond distance (A)", "title": "t",
    }])

    figure = plt.figure(1)
    # Two axes: the bars' own, plus the twin carrying the percentage line.
    assert len(figure.axes) == 2
    assert figure.axes[0].get_xticklabels()[0].get_text() == "C1-C2"
    plt.close("all")


def test_viewer_renders_a_bar_figure_without_a_line():
    import matplotlib.pyplot as plt
    import plotViewer

    plt.close("all")
    plotViewer.build_figures([{
        "type": "bars",
        "categories": CATEGORIES,
        "groups": [{"label": "isolated", "values": [1.52, 1.43]}],
        "xlabel": "parameter", "ylabel": "d (A)", "title": "t",
    }])

    assert len(plt.figure(1).axes) == 1
    plt.close("all")


def test_viewer_still_renders_the_two_older_manifest_shapes():
    """Regression guard: every other tool writes these, and they carry no
    'type' key at all."""
    import matplotlib.pyplot as plt
    import plotViewer

    plt.close("all")
    plotViewer.build_figures([
        {"x": [0, 1], "y": [0, 1], "xlabel": "x", "ylabel": "y", "title": "curve"},
        {"series": [{"x": [0, 1], "y": [1, 0], "label": "a"}],
         "xlabel": "x", "ylabel": "y", "title": "overlay"},
    ])

    assert plt.figure(1).axes[0].get_title() == "curve"
    assert plt.figure(2).axes[0].get_title() == "overlay"
    assert len(plt.figure(2).axes[0].get_legend().get_texts()) == 1
    plt.close("all")
