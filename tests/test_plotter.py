"""Tests for the energy-plot tool (``plotter.py``), JSON plot type.

Run from ``venv/src/``::

    python -m pytest tests/ -q

The Toga GUI is never instantiated. ``validate_plot_manifest`` and
``load_json_plot_file`` are pure, so they are driven directly; the UI wiring is
checked by constructing ``PlotterUI`` with ``__new__`` and stub widgets.

The JSON accepted here is the *plot manifest* every gQTEA analysis tool already
writes for the interactive viewer (see ``displayPlots.save_plots``), so a file
produced by any tool in the suite can be re-opened and re-plotted later.
"""
import json

import pytest

from plotter import PlotterBase


# --------------------------------------------------------------------------- #
# Fixtures                                                                      #
# --------------------------------------------------------------------------- #
def single_curve(n=5):
    return {
        "x": [float(i) for i in range(n)],
        "y": [float(i * i) for i in range(n)],
        "xlabel": "time (fs)",
        "ylabel": "S(t)",
        "title": "Continuous survival",
    }


def multi_series(n=4):
    return {
        "series": [
            {"x": [float(i) for i in range(n)],
             "y": [float(i) for i in range(n)], "label": "C(t)"},
            {"x": [float(i) for i in range(n)],
             "y": [float(-i) for i in range(n)], "label": "R(t)"},
        ],
        "xlabel": "time (fs)",
        "ylabel": "correlation",
        "title": "Intermittent correlation",
    }


def bar_manifest():
    return {
        "type": "bars",
        "categories": ["C1-C2", "C2-O3"],
        "groups": [
            {"label": "isolated", "values": [1.52, 1.43], "errors": [0.01, 0.02]},
            {"label": "solvated", "values": [1.54, 1.42], "errors": [0.02, 0.01]},
        ],
        "line": {"label": "%dr", "values": [1.25, -0.84], "ylabel": "percent (%)"},
        "xlabel": "parameter", "ylabel": "bond distance (A)", "title": "t",
    }


def write(tmp_path, payload, name="manifest.json"):
    path = tmp_path / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return str(path)


# --------------------------------------------------------------------------- #
# Accepting well-formed manifests                                               #
# --------------------------------------------------------------------------- #
class TestValidManifests:
    def test_single_curve(self):
        figures = PlotterBase.validate_plot_manifest([single_curve()])
        assert len(figures) == 1

    def test_multi_series(self):
        figures = PlotterBase.validate_plot_manifest([multi_series()])
        assert len(figures) == 1

    def test_mixed_manifest(self):
        payload = [single_curve(), multi_series(), single_curve()]
        assert len(PlotterBase.validate_plot_manifest(payload)) == 3

    def test_optional_limits_are_accepted(self):
        fig = single_curve()
        fig["xlim"] = [0.0, 4.0]
        fig["ylim"] = [0, 16]
        PlotterBase.validate_plot_manifest([fig])

    def test_labels_are_optional(self):
        PlotterBase.validate_plot_manifest([{"x": [1, 2], "y": [3, 4]}])

    def test_integer_data_is_accepted(self):
        PlotterBase.validate_plot_manifest([{"x": [1, 2, 3], "y": [4, 5, 6]}])

    def test_a_real_file_round_trips(self, tmp_path):
        path = write(tmp_path, [single_curve(), multi_series()])
        figures = PlotterBase.load_json_plot_file(path)
        assert len(figures) == 2

    def test_bars_manifest_is_accepted(self):
        figures = PlotterBase.validate_plot_manifest([bar_manifest()])
        assert len(figures) == 1

    def test_bars_manifest_without_errors_or_line_is_accepted(self):
        fig = bar_manifest()
        del fig["line"]
        del fig["groups"][0]["errors"]
        del fig["groups"][1]["errors"]
        PlotterBase.validate_plot_manifest([fig])

    def test_mixed_manifest_with_bars_curve_and_series(self):
        payload = [bar_manifest(), single_curve(), multi_series()]
        assert len(PlotterBase.validate_plot_manifest(payload)) == 3


# --------------------------------------------------------------------------- #
# Rejecting malformed manifests                                                 #
# --------------------------------------------------------------------------- #
class TestInvalidManifests:
    def test_missing_file(self, tmp_path):
        with pytest.raises(ValueError, match="not found"):
            PlotterBase.load_json_plot_file(str(tmp_path / "nope.json"))

    def test_not_json_at_all(self, tmp_path):
        path = tmp_path / "bad.json"
        path.write_text("this is not JSON {{{", encoding="utf-8")
        with pytest.raises(ValueError, match="not valid JSON"):
            PlotterBase.load_json_plot_file(str(path))

    def test_empty_file(self, tmp_path):
        path = tmp_path / "empty.json"
        path.write_text("", encoding="utf-8")
        with pytest.raises(ValueError, match="not valid JSON"):
            PlotterBase.load_json_plot_file(str(path))

    def test_top_level_object_instead_of_array(self):
        with pytest.raises(ValueError, match="list of figures"):
            PlotterBase.validate_plot_manifest({"x": [1], "y": [2]})

    def test_empty_array(self):
        with pytest.raises(ValueError, match="no figures"):
            PlotterBase.validate_plot_manifest([])

    def test_entry_is_not_an_object(self):
        with pytest.raises(ValueError, match="Figure 2"):
            PlotterBase.validate_plot_manifest([single_curve(), "oops"])

    def test_entry_without_x_y_or_series(self):
        with pytest.raises(ValueError, match="x.*y.*series"):
            PlotterBase.validate_plot_manifest([{"title": "nothing to plot"}])

    def test_x_and_y_length_mismatch(self):
        with pytest.raises(ValueError, match="same length"):
            PlotterBase.validate_plot_manifest([{"x": [1, 2, 3], "y": [1, 2]}])

    def test_empty_curve(self):
        with pytest.raises(ValueError, match="no data"):
            PlotterBase.validate_plot_manifest([{"x": [], "y": []}])

    def test_non_numeric_values(self):
        with pytest.raises(ValueError, match="numbers"):
            PlotterBase.validate_plot_manifest([{"x": [1, 2], "y": [1, "two"]}])

    def test_x_is_not_a_list(self):
        with pytest.raises(ValueError, match="list"):
            PlotterBase.validate_plot_manifest([{"x": 1, "y": [1]}])

    def test_series_not_a_list(self):
        with pytest.raises(ValueError, match="series"):
            PlotterBase.validate_plot_manifest([{"series": {"x": [1], "y": [2]}}])

    def test_empty_series_list(self):
        with pytest.raises(ValueError, match="series"):
            PlotterBase.validate_plot_manifest([{"series": []}])

    def test_series_entry_missing_y(self):
        with pytest.raises(ValueError, match="series 1"):
            PlotterBase.validate_plot_manifest([{"series": [{"x": [1, 2]}]}])

    def test_series_entry_length_mismatch(self):
        payload = [{"series": [{"x": [1, 2, 3], "y": [1, 2]}]}]
        with pytest.raises(ValueError, match="same length"):
            PlotterBase.validate_plot_manifest(payload)

    def test_bad_xlim(self):
        fig = single_curve()
        fig["xlim"] = [1.0]
        with pytest.raises(ValueError, match="xlim"):
            PlotterBase.validate_plot_manifest([fig])

    def test_error_message_identifies_the_figure(self):
        payload = [single_curve(), single_curve(), {"x": [1], "y": [1, 2]}]
        with pytest.raises(ValueError, match="Figure 3"):
            PlotterBase.validate_plot_manifest(payload)

    def test_nan_is_rejected(self):
        """json.load turns NaN into float('nan'); a plot of NaN is not useful."""
        with pytest.raises(ValueError, match="numbers"):
            PlotterBase.validate_plot_manifest(
                [{"x": [1.0, 2.0], "y": [1.0, float("nan")]}]
            )

    def test_bars_group_values_length_mismatch(self):
        fig = bar_manifest()
        fig["groups"][0]["values"] = [1.52]  # one category short
        with pytest.raises(ValueError, match="Figure 1.*same length"):
            PlotterBase.validate_plot_manifest([fig])

    def test_bars_group_errors_length_mismatch(self):
        fig = bar_manifest()
        fig["groups"][0]["errors"] = [0.01]
        with pytest.raises(ValueError, match="Figure 1.*errors.*same length"):
            PlotterBase.validate_plot_manifest([fig])

    def test_bars_line_length_mismatch(self):
        fig = bar_manifest()
        fig["line"]["values"] = [1.25]
        with pytest.raises(ValueError, match="Figure 1.*line.*same length"):
            PlotterBase.validate_plot_manifest([fig])

    def test_bars_missing_categories(self):
        fig = bar_manifest()
        del fig["categories"]
        with pytest.raises(ValueError, match="Figure 1.*categories"):
            PlotterBase.validate_plot_manifest([fig])

    def test_bars_empty_categories(self):
        fig = bar_manifest()
        fig["categories"] = []
        with pytest.raises(ValueError, match="Figure 1.*categories"):
            PlotterBase.validate_plot_manifest([fig])

    def test_bars_missing_groups(self):
        fig = bar_manifest()
        del fig["groups"]
        with pytest.raises(ValueError, match="Figure 1.*groups"):
            PlotterBase.validate_plot_manifest([fig])

    def test_bars_group_missing_values(self):
        fig = bar_manifest()
        del fig["groups"][0]["values"]
        with pytest.raises(ValueError, match="Figure 1.*values"):
            PlotterBase.validate_plot_manifest([fig])

    def test_bars_error_names_the_offending_figure(self):
        fig = bar_manifest()
        fig["groups"][0]["values"] = [1.52]
        payload = [single_curve(), fig]
        with pytest.raises(ValueError, match="Figure 2"):
            PlotterBase.validate_plot_manifest(payload)


# --------------------------------------------------------------------------- #
# describe_json_figures                                                         #
# --------------------------------------------------------------------------- #
class TestDescribeJsonFigures:
    def test_single_curve_point_count(self):
        text = PlotterBase.describe_json_figures([single_curve(5)])
        assert "5 points" in text

    def test_series_curve_count(self):
        text = PlotterBase.describe_json_figures([multi_series(4)])
        assert "2 curves, 4 points each" in text

    def test_bars_group_and_category_count(self):
        text = PlotterBase.describe_json_figures([bar_manifest()])
        assert "2 groups, 2 categories" in text


# --------------------------------------------------------------------------- #
# Plot-type plumbing                                                            #
# --------------------------------------------------------------------------- #
class FakeWidget:
    def __init__(self, value=None, text=""):
        self.value = value
        self.text = text
        self.enabled = True
        self.placeholder = ""


class TestPlotTypeWiring:
    def test_the_three_plot_types_are_distinct(self):
        types = {PlotterBase.CPMD_PLOT_TYPE,
                 PlotterBase.GQTEAMD_PLOT_TYPE,
                 PlotterBase.JSON_PLOT_TYPE}
        assert len(types) == 3

    def test_json_type_is_detected(self):
        obj = PlotterBase.__new__(PlotterBase)
        obj.plot_type_selection = FakeWidget(value=PlotterBase.JSON_PLOT_TYPE)
        assert obj.is_json_plot_type()
        assert not obj.is_gqtea_plot_type()

    def test_existing_types_are_not_json(self):
        obj = PlotterBase.__new__(PlotterBase)
        for value in (PlotterBase.CPMD_PLOT_TYPE, PlotterBase.GQTEAMD_PLOT_TYPE):
            obj.plot_type_selection = FakeWidget(value=value)
            assert not obj.is_json_plot_type()

    def test_gqtea_detection_is_unchanged(self):
        obj = PlotterBase.__new__(PlotterBase)
        obj.plot_type_selection = FakeWidget(value=PlotterBase.GQTEAMD_PLOT_TYPE)
        assert obj.is_gqtea_plot_type()


class TestControlSwitching:
    @staticmethod
    def bare_ui():
        from plotter import PlotterUI
        ui = PlotterUI.__new__(PlotterUI)
        ui.gqtea_headers = []
        ui.gqtea_y_switches = []
        ui.json_figures = []
        ui.data = []
        ui.x_axis = []
        for name in ("switch_fictitious_ionic", "switch_temperature",
                     "switch_ksh_energy", "switch_ksh_ionic",
                     "switch_total_energy", "switch_cpu_time"):
            setattr(ui, name, FakeWidget(value=True))
        for name in ("text_input_time_step", "unit_selection", "time_step_label",
                     "units_label", "gqtea_x_axis_label", "gqtea_x_axis_selection",
                     "gqtea_y_columns_label", "file_label", "text_input_file",
                     "multi_line_text"):
            setattr(ui, name, FakeWidget(value=""))
        ui.plot_type_selection = FakeWidget(value=PlotterBase.JSON_PLOT_TYPE)
        return ui

    def test_json_mode_disables_every_other_option(self):
        ui = self.bare_ui()
        ui.on_plot_type_change(ui.plot_type_selection)

        for switch in ui.get_cpmd_switches():
            assert switch.enabled is False, "CPMD switches must be disabled"
            assert switch.value is False, "CPMD switches must be cleared"
        for widget in (ui.text_input_time_step, ui.unit_selection,
                       ui.time_step_label, ui.units_label,
                       ui.gqtea_x_axis_label, ui.gqtea_x_axis_selection,
                       ui.gqtea_y_columns_label):
            assert widget.enabled is False, "all plot-specific options must be off"

    def test_json_mode_relabels_the_file_row(self):
        ui = self.bare_ui()
        ui.on_plot_type_change(ui.plot_type_selection)
        assert "JSON" in ui.file_label.text
        assert "JSON" in ui.text_input_file.placeholder

    def test_switching_away_from_json_restores_cpmd_controls(self):
        ui = self.bare_ui()
        ui.on_plot_type_change(ui.plot_type_selection)
        ui.plot_type_selection.value = PlotterBase.CPMD_PLOT_TYPE
        ui.on_plot_type_change(ui.plot_type_selection)

        for switch in ui.get_cpmd_switches():
            assert switch.enabled is True
        assert ui.text_input_time_step.enabled is True
        assert ui.unit_selection.enabled is True
        assert ui.gqtea_x_axis_label.enabled is False

    def test_switching_away_from_json_restores_gqtea_controls(self):
        ui = self.bare_ui()
        ui.on_plot_type_change(ui.plot_type_selection)
        ui.plot_type_selection.value = PlotterBase.GQTEAMD_PLOT_TYPE
        ui.on_plot_type_change(ui.plot_type_selection)

        assert ui.gqtea_x_axis_label.enabled is True
        assert ui.gqtea_y_columns_label.enabled is True
        for switch in ui.get_cpmd_switches():
            assert switch.enabled is False

    def test_changing_type_clears_loaded_json(self):
        ui = self.bare_ui()
        ui.json_figures = [single_curve()]
        ui.plot_type_selection.value = PlotterBase.CPMD_PLOT_TYPE
        ui.on_plot_type_change(ui.plot_type_selection)
        assert ui.json_figures == [], "stale figures must not survive a mode change"
