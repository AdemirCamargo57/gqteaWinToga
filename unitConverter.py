from __future__ import annotations

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW

from unitConverterBackEnd import EnergyConverter


class UnitConvertUI:
    """Window-based Toga interface for the unit converter."""

    def __init__(self, *args) -> None:
        del args
        self.converter = EnergyConverter()
        self.category = "energy"
        self.category_selection: toga.Selection | None = None
        self.value_input: toga.TextInput | None = None
        self.from_selection: toga.Selection | None = None
        self.to_selection: toga.Selection | None = None
        self.result_label: toga.Label | None = None
        self.status_label: toga.Label | None = None
        self.category_label_to_key = {
            "Energy": "energy",
            "Length": "length",
            "Density": "density",
            "Time": "time",
        }
        self.default_units = {
            "energy": ("hartree", "ev"),
            "length": ("meter", "nanometer"),
            "density": ("amu/angstrom^3", "g/cm^3"),
            "time": ("second", "femtosecond"),
        }
        self.current_unit_labels: dict[str, str] = {}
        self.main_window = toga.Window(title="Unit Converter", size=(520, 340))
        self.build_ui()
        self.main_window.show()

    def build_ui(self) -> None:
        self.category_selection = toga.Selection(
            items=list(self.category_label_to_key.keys()),
            style=Pack(width=150, margin_bottom=10),
            on_change=self.handle_category_change,
        )
        self.value_input = toga.TextInput(
            placeholder="Enter a value (for example 1.23e-4)",
            style=Pack(flex=1, margin_bottom=10),
        )
        self.from_selection = toga.Selection(
            items=[],
            style=Pack(flex=1, margin_right=5),
        )
        self.to_selection = toga.Selection(
            items=[],
            style=Pack(flex=1, margin_left=5),
        )

        convert_button = toga.Button(
            "Convert",
            on_press=self.handle_convert,
            style=Pack(margin_top=10, width=140),
        )

        self.result_label = toga.Label(
            "Result will appear here.",
            style=Pack(margin_top=12, font_size=12),
        )
        self.status_label = toga.Label(
            "",
            style=Pack(margin_top=6, color="#555555"),
        )

        content = toga.Box(
            style=Pack(direction=COLUMN, margin=20),
            children=[
                toga.Label(
                    "Unit Converter",
                    style=Pack(font_size=18, font_weight="bold", margin_bottom=10),
                ),
                toga.Label("Category", style=Pack(margin_bottom=4)),
                self.category_selection,
                toga.Label("Value", style=Pack(margin_bottom=4)),
                self.value_input,
                toga.Box(
                    style=Pack(direction=ROW, margin_top=8),
                    children=[
                        toga.Box(
                            style=Pack(direction=COLUMN, flex=1, margin_right=6),
                            children=[
                                toga.Label("From", style=Pack(margin_bottom=4)),
                                self.from_selection,
                            ],
                        ),
                        toga.Box(
                            style=Pack(direction=COLUMN, flex=1, margin_left=6),
                            children=[
                                toga.Label("To", style=Pack(margin_bottom=4)),
                                self.to_selection,
                            ],
                        ),
                    ],
                ),
                convert_button,
                self.result_label,
                self.status_label,
            ],
        )

        self.main_window.content = content
        self.category_selection.value = "Energy"
        self.update_unit_selections()

    def update_unit_selections(self) -> None:
        assert self.from_selection is not None
        assert self.to_selection is not None
        assert self.status_label is not None

        self.current_unit_labels = self.converter.get_display_map(self.category)
        unit_items = list(self.current_unit_labels.values())
        self.from_selection.items = unit_items
        self.to_selection.items = unit_items

        default_from_key, default_to_key = self.default_units[self.category]
        self.from_selection.value = self.current_unit_labels[default_from_key]
        self.to_selection.value = self.current_unit_labels[default_to_key]

        if self.category == "energy":
            self.status_label.text = "Energy conversions are normalized through joules."
        elif self.category == "length":
            self.status_label.text = "Length conversions are normalized through meters."
        elif self.category == "density":
            self.status_label.text = "Density conversions are normalized through grams per cubic centimeter."
        else:
            self.status_label.text = "Time conversions are normalized through seconds. Month = year/12 and year = 365 days."
        self.status_label.style.color = "#555555"

    def handle_category_change(self, widget: toga.Selection) -> None:
        if widget.value is None:
            return
        self.category = self.category_label_to_key[widget.value]
        if self.result_label is None:
            return
        self.result_label.text = "Result will appear here."
        self.update_unit_selections()

    def handle_convert(self, widget: toga.Button) -> None:
        del widget
        assert self.category_selection is not None
        assert self.value_input is not None
        assert self.from_selection is not None
        assert self.to_selection is not None
        assert self.result_label is not None
        assert self.status_label is not None

        try:
            if self.from_selection.value is None or self.to_selection.value is None:
                raise ValueError("Please choose both units before converting.")

            from_unit = self._label_to_unit_key(self.from_selection.value)
            to_unit = self._label_to_unit_key(self.to_selection.value)
            result_text = self.converter.describe_result(
                self.value_input.value,
                from_unit,
                to_unit,
                category=self.category,
                use_labels=True,
            )
            self.result_label.text = result_text
            self.status_label.text = "Conversion completed successfully."
            self.status_label.style.color = "#1e6b34"
        except ValueError as exc:
            self.result_label.text = "Conversion failed."
            self.status_label.text = str(exc)
            self.status_label.style.color = "#9f1239"

    def _label_to_unit_key(self, label: str) -> str:
        for key, display_label in self.current_unit_labels.items():
            if display_label == label:
                return key
        raise ValueError(f"Unsupported unit '{label}'.")
