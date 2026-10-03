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


class FakeBox:
    def __init__(self):
        self.children = []

    def add(self, child):
        self.children.append(child)

    def remove(self, child):
        self.children.remove(child)


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
                     "multi_line_text", "table_xlabel_input", "table_ylabel_input",
                     "table_styles_label"):
            setattr(ui, name, FakeWidget(value=""))
        ui.table_styles_box = FakeBox()
        ui.table_style_selections = []
        ui.table_column_switches = []
        ui._section_slots = {}
        for name in ("cpmd_switch_box", "time_step_box", "units_box",
                     "gqtea_column_box", "table_options_box"):
            section, slot = FakeWidget(), FakeBox()
            slot.add(section)
            setattr(ui, name, section)
            ui._section_slots[name] = slot
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


# --------------------------------------------------------------------------- #
# Multi-column table plot type                                                  #
# --------------------------------------------------------------------------- #
class TestParseTableText:
    def test_plain_numeric_columns_get_generic_labels(self):
        table = PlotterBase.parse_table_text("0 1 2\n1 3 4\n2 5 6\n")
        assert table["x"] == [0.0, 1.0, 2.0]
        assert table["ys"] == [[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]]
        assert table["xlabel"] == "Column 1"
        assert table["labels"] == ["Column 2", "Column 3"]
        assert table["has_header"] is False

    def test_header_names_become_legend_labels(self):
        text = "time E_kin E_pot E_tot\n0 1 2 3\n1 4 5 6\n"
        table = PlotterBase.parse_table_text(text)
        assert table["xlabel"] == "time"
        assert table["labels"] == ["E_kin", "E_pot", "E_tot"]
        assert table["has_header"] is True
        assert len(table["ys"]) == 3

    def test_arbitrary_number_of_y_columns(self):
        n_cols = 12
        rows = [" ".join(str(float(r * c)) for c in range(n_cols)) for r in range(5)]
        table = PlotterBase.parse_table_text("\n".join(rows))
        assert len(table["ys"]) == n_cols - 1
        assert all(len(y) == 5 for y in table["ys"])

    def test_single_y_column(self):
        table = PlotterBase.parse_table_text("0 10\n1 20\n")
        assert table["ys"] == [[10.0, 20.0]]

    def test_comment_lines_are_ignored(self):
        text = (
            "# produced by some tool\n"
            "! fortran-style comment\n"
            "% matlab-style comment\n"
            "@ xmgrace directive\n"
            "// c-style comment\n"
            "r g(r) integral\n"
            "# a comment between header and data\n"
            "0.1 0.0 0.0\n"
            "\n"
            "0.2 1.5 0.3   # trailing inline comment\n"
            "# trailing comment\n"
        )
        table = PlotterBase.parse_table_text(text)
        assert table["x"] == [0.1, 0.2]
        assert table["labels"] == ["g(r)", "integral"]
        assert table["ys"] == [[0.0, 1.5], [0.0, 0.3]]

    def test_comment_only_preamble_without_header(self):
        text = "# x y1 y2\n# another comment\n1 2 3\n4 5 6\n"
        table = PlotterBase.parse_table_text(text)
        assert table["has_header"] is False
        assert table["labels"] == ["Column 2", "Column 3"]
        assert table["x"] == [1.0, 4.0]

    def test_header_is_never_read_as_data(self):
        # A header that *looks* half-numeric must still not become a data row.
        text = "x 1 2\n0 5 6\n1 7 8\n"
        table = PlotterBase.parse_table_text(text)
        assert table["x"] == [0.0, 1.0]
        assert table["labels"] == ["1", "2"]

    def test_csv_with_header(self):
        text = "step, temperature (K), pressure\n0, 300, 1.0\n1, 301.5, 1.1\n"
        table = PlotterBase.parse_table_text(text)
        assert table["xlabel"] == "step"
        assert table["labels"] == ["temperature (K)", "pressure"]
        assert table["ys"][0] == [300.0, 301.5]

    def test_tab_separated(self):
        table = PlotterBase.parse_table_text("t\ta\tb\n0\t1\t2\n1\t3\t4\n")
        assert table["labels"] == ["a", "b"]

    def test_scientific_notation(self):
        table = PlotterBase.parse_table_text("1.0E-3 -2.5e+02\n2.0e-3 1e2\n")
        assert table["x"] == [1.0e-3, 2.0e-3]
        assert table["ys"] == [[-250.0, 100.0]]

    def test_last_text_line_before_data_is_the_header(self):
        text = "My simulation title\nx alpha beta\n0 1 2\n"
        table = PlotterBase.parse_table_text(text)
        assert table["labels"] == ["alpha", "beta"]


class TestParseTableTextErrors:
    def test_empty(self):
        with pytest.raises(ValueError, match="no numeric data"):
            PlotterBase.parse_table_text("")

    def test_only_comments_and_header(self):
        with pytest.raises(ValueError, match="no numeric data"):
            PlotterBase.parse_table_text("# c\nx y\n")

    def test_single_column(self):
        with pytest.raises(ValueError, match="at least two columns"):
            PlotterBase.parse_table_text("1\n2\n3\n")

    def test_ragged_row_names_the_line(self):
        with pytest.raises(ValueError, match="Line 3.*3 values.*2"):
            PlotterBase.parse_table_text("x y\n0 1\n1 2 3\n")

    def test_text_inside_data_names_the_line(self):
        with pytest.raises(ValueError, match="Line 3.*not numeric"):
            PlotterBase.parse_table_text("0 1\n1 2\nfoo bar\n")

    def test_non_finite_values_rejected(self):
        with pytest.raises(ValueError, match="Line 2.*finite"):
            PlotterBase.parse_table_text("0 1\n1 nan\n")

    def test_header_column_count_mismatch(self):
        with pytest.raises(ValueError, match="header.*2 names.*3 columns"):
            PlotterBase.parse_table_text("x y\n0 1 2\n")


class TestLoadTableFile:
    def test_reads_from_disk(self, tmp_path):
        path = tmp_path / "data.dat"
        path.write_text("# c\nx a b\n0 1 2\n1 3 4\n", encoding="utf-8")
        table = PlotterBase.load_table_file(str(path))
        assert table["labels"] == ["a", "b"]

    def test_missing_file(self, tmp_path):
        with pytest.raises(ValueError, match="not found"):
            PlotterBase.load_table_file(str(tmp_path / "nope.dat"))

    def test_table_to_figure_is_a_valid_manifest(self):
        table = PlotterBase.parse_table_text("x a b\n0 1 2\n1 3 4\n")
        figure = PlotterBase.table_to_figure(table, title="data.dat")
        assert [s["label"] for s in figure["series"]] == ["a", "b"]
        assert figure["xlabel"] == "x"
        assert figure["title"] == "data.dat"
        PlotterBase.validate_plot_manifest([figure])

    def test_describe_table(self):
        table = PlotterBase.parse_table_text("x a b\n0 1 2\n1 3 4\n")
        text = PlotterBase.describe_table(table)
        assert "2 data rows" in text
        assert "2 y-axis series" in text
        assert "a" in text and "b" in text


class TestTableModeWiring:
    def test_four_plot_types_are_distinct(self):
        types = {PlotterBase.CPMD_PLOT_TYPE, PlotterBase.GQTEAMD_PLOT_TYPE,
                 PlotterBase.JSON_PLOT_TYPE, PlotterBase.TABLE_PLOT_TYPE}
        assert len(types) == 4

    def test_table_type_is_detected(self):
        obj = PlotterBase.__new__(PlotterBase)
        obj.plot_type_selection = FakeWidget(value=PlotterBase.TABLE_PLOT_TYPE)
        assert obj.is_table_plot_type()
        assert not obj.is_json_plot_type()
        assert not obj.is_gqtea_plot_type()

    def test_table_mode_disables_every_other_option(self):
        ui = TestControlSwitching.bare_ui()
        ui.table_data = None
        ui.plot_type_selection.value = PlotterBase.TABLE_PLOT_TYPE
        ui.on_plot_type_change(ui.plot_type_selection)
        for switch in ui.get_cpmd_switches():
            assert switch.enabled is False
            assert switch.value is False
        for widget in (ui.text_input_time_step, ui.unit_selection,
                       ui.gqtea_x_axis_label, ui.gqtea_x_axis_selection,
                       ui.gqtea_y_columns_label):
            assert widget.enabled is False
        assert "column" in ui.file_label.text.lower()

    def test_changing_type_clears_loaded_table(self):
        ui = TestControlSwitching.bare_ui()
        ui.table_data = PlotterBase.parse_table_text("0 1\n1 2\n")
        ui.plot_type_selection.value = PlotterBase.CPMD_PLOT_TYPE
        ui.on_plot_type_change(ui.plot_type_selection)
        assert ui.table_data is None


# --------------------------------------------------------------------------- #
# Table plot customisation: axis labels and per-column line styles             #
# --------------------------------------------------------------------------- #
def _table():
    return PlotterBase.parse_table_text("t a b c\n0 1 2 3\n1 4 5 6\n")


class TestTableFigureCustomisation:
    def test_default_axis_labels(self):
        figure = PlotterBase.table_to_figure(_table())
        assert figure["xlabel"] == "t"
        assert figure["ylabel"] == "Value"

    def test_custom_axis_labels(self):
        figure = PlotterBase.table_to_figure(
            _table(), xlabel="Time (fs)", ylabel="Energy (Ha)")
        assert figure["xlabel"] == "Time (fs)"
        assert figure["ylabel"] == "Energy (Ha)"

    def test_blank_axis_labels_fall_back_to_defaults(self):
        figure = PlotterBase.table_to_figure(_table(), xlabel="   ", ylabel="")
        assert figure["xlabel"] == "t"
        assert figure["ylabel"] == "Value"

    def test_default_style_is_a_solid_line(self):
        figure = PlotterBase.table_to_figure(_table())
        for curve in figure["series"]:
            assert curve["linestyle"] == "-"
            assert "marker" not in curve

    def test_each_column_gets_its_own_style(self):
        figure = PlotterBase.table_to_figure(
            _table(), styles=["Dashed", "Markers only", "Line + markers"])
        a, b, c = figure["series"]
        assert a["linestyle"] == "--" and "marker" not in a
        assert b["linestyle"] == "none" and b["marker"] == "o"
        assert c["linestyle"] == "-" and c["marker"] == "o"

    def test_every_offered_style_is_a_valid_manifest(self):
        for style in PlotterBase.TABLE_LINE_STYLES:
            figure = PlotterBase.table_to_figure(_table(), styles=[style] * 3)
            PlotterBase.validate_plot_manifest([figure])

    def test_wrong_number_of_styles_rejected(self):
        with pytest.raises(ValueError, match="3 y columns.*2 line styles"):
            PlotterBase.table_to_figure(_table(), styles=["Solid", "Dashed"])

    def test_unknown_style_rejected(self):
        with pytest.raises(ValueError, match="Unknown line style"):
            PlotterBase.table_to_figure(_table(), styles=["Solid", "Wavy", "Solid"])

    def test_styled_figure_renders_with_matplotlib(self):
        import matplotlib.pyplot as plt
        from plotViewer import build_figures
        figure = PlotterBase.table_to_figure(
            _table(), styles=["Dotted", "Dash-dot", "Markers only"])
        build_figures([figure])
        lines = plt.gcf().gca().get_lines()
        assert [ln.get_linestyle() for ln in lines] == [":", "-.", "None"]
        assert lines[2].get_marker() == "o"
        plt.close("all")


class TestSeriesStyleValidation:
    def test_unknown_linestyle_rejected(self):
        fig = multi_series()
        fig["series"][0]["linestyle"] = "wavy"
        with pytest.raises(ValueError, match="Figure 1, series 1.*linestyle"):
            PlotterBase.validate_plot_manifest([fig])

    def test_unknown_marker_rejected(self):
        fig = multi_series()
        fig["series"][1]["marker"] = 3
        with pytest.raises(ValueError, match="Figure 1, series 2.*marker"):
            PlotterBase.validate_plot_manifest([fig])

    def test_unstyled_series_still_valid(self):
        PlotterBase.validate_plot_manifest([multi_series()])


class TestTableStyleControls:
    def test_mode_change_clears_style_rows_and_axis_labels(self):
        ui = TestControlSwitching.bare_ui()
        ui.table_data = None
        ui.table_styles_box.add(FakeWidget())
        ui.table_style_selections = [FakeWidget(value="Dashed")]
        ui.table_xlabel_input.value = "old x"
        ui.plot_type_selection.value = PlotterBase.CPMD_PLOT_TYPE
        ui.on_plot_type_change(ui.plot_type_selection)
        assert ui.table_styles_box.children == []
        assert ui.table_style_selections == []
        assert ui.table_xlabel_input.value == ""
        assert ui.table_xlabel_input.enabled is False
        assert ui.table_ylabel_input.enabled is False

    def test_table_mode_enables_axis_label_inputs(self):
        ui = TestControlSwitching.bare_ui()
        ui.table_data = None
        ui.plot_type_selection.value = PlotterBase.TABLE_PLOT_TYPE
        ui.on_plot_type_change(ui.plot_type_selection)
        assert ui.table_xlabel_input.enabled is True
        assert ui.table_ylabel_input.enabled is True

    def test_selected_styles_are_read_in_column_order(self):
        ui = TestControlSwitching.bare_ui()
        ui.table_style_selections = [FakeWidget(value="Dotted"),
                                     FakeWidget(value="Solid")]
        assert ui.selected_table_styles() == ["Dotted", "Solid"]


# --------------------------------------------------------------------------- #
# Table plot: choosing which y columns to plot                                  #
# --------------------------------------------------------------------------- #
class TestTableColumnSelection:
    def test_all_columns_plotted_by_default(self):
        figure = PlotterBase.table_to_figure(_table())
        assert [s["label"] for s in figure["series"]] == ["a", "b", "c"]

    def test_only_selected_columns_are_plotted(self):
        figure = PlotterBase.table_to_figure(_table(), include=[True, False, True])
        assert [s["label"] for s in figure["series"]] == ["a", "c"]
        assert figure["series"][1]["y"] == [3.0, 6.0]

    def test_styles_stay_attached_to_their_column(self):
        figure = PlotterBase.table_to_figure(
            _table(), styles=["Dashed", "Dotted", "Markers only"],
            include=[False, True, True])
        b, c = figure["series"]
        assert (b["label"], b["linestyle"]) == ("b", ":")
        assert (c["label"], c["marker"]) == ("c", "o")

    def test_no_column_selected_is_rejected(self):
        with pytest.raises(ValueError, match="at least one y column"):
            PlotterBase.table_to_figure(_table(), include=[False, False, False])

    def test_wrong_number_of_include_flags_rejected(self):
        with pytest.raises(ValueError, match="3 y columns.*2 plot selections"):
            PlotterBase.table_to_figure(_table(), include=[True, False])

    def test_partial_selection_is_a_valid_manifest(self):
        figure = PlotterBase.table_to_figure(_table(), include=[False, True, False])
        PlotterBase.validate_plot_manifest([figure])


class TestTableColumnSwitches:
    def test_selected_columns_are_read_in_column_order(self):
        ui = TestControlSwitching.bare_ui()
        ui.table_column_switches = [FakeWidget(value=True), FakeWidget(value=False),
                                    FakeWidget(value=True)]
        assert ui.selected_table_columns() == [True, False, True]

    def test_mode_change_clears_column_switches(self):
        ui = TestControlSwitching.bare_ui()
        ui.table_data = None
        ui.table_column_switches = [FakeWidget(value=True)]
        ui.plot_type_selection.value = PlotterBase.CPMD_PLOT_TYPE
        ui.on_plot_type_change(ui.plot_type_selection)
        assert ui.table_column_switches == []


# --------------------------------------------------------------------------- #
# Window height: only the active plot type's section takes up space            #
# --------------------------------------------------------------------------- #
class TestSectionVisibility:
    SECTIONS = {
        PlotterBase.CPMD_PLOT_TYPE: {"cpmd_switch_box", "time_step_box", "units_box"},
        PlotterBase.GQTEAMD_PLOT_TYPE: {"gqtea_column_box"},
        PlotterBase.JSON_PLOT_TYPE: set(),
        PlotterBase.TABLE_PLOT_TYPE: {"table_options_box"},
    }

    @pytest.mark.parametrize("plot_type", list(SECTIONS))
    def test_only_the_active_section_is_shown(self, plot_type):
        ui = TestControlSwitching.bare_ui()
        ui.table_data = None
        ui.plot_type_selection.value = plot_type
        ui.on_plot_type_change(ui.plot_type_selection)
        all_sections = set().union(*self.SECTIONS.values())
        shown = {name for name in all_sections
                 if getattr(ui, name) in ui._section_slots[name].children}
        assert shown == self.SECTIONS[plot_type]

    def test_switching_back_shows_the_section_again(self):
        ui = TestControlSwitching.bare_ui()
        ui.table_data = None
        ui.plot_type_selection.value = PlotterBase.TABLE_PLOT_TYPE
        ui.on_plot_type_change(ui.plot_type_selection)
        ui.plot_type_selection.value = PlotterBase.CPMD_PLOT_TYPE
        ui.on_plot_type_change(ui.plot_type_selection)
        assert ui._section_slots["cpmd_switch_box"].children == [ui.cpmd_switch_box]
        assert ui._section_slots["table_options_box"].children == []

    def test_repeated_changes_never_add_a_section_twice(self):
        ui = TestControlSwitching.bare_ui()
        ui.table_data = None
        for plot_type in (PlotterBase.CPMD_PLOT_TYPE, PlotterBase.CPMD_PLOT_TYPE,
                          PlotterBase.TABLE_PLOT_TYPE, PlotterBase.CPMD_PLOT_TYPE):
            ui.plot_type_selection.value = plot_type
            ui.on_plot_type_change(ui.plot_type_selection)
        assert all(len(slot.children) <= 1 for slot in ui._section_slots.values())
