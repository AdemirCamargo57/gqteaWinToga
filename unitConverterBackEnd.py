from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation, localcontext
from typing import Dict, Iterable, Optional


AVOGADRO_CONSTANT = Decimal("6.02214076e23")
ELEMENTARY_CHARGE = Decimal("1.602176634e-19")
HARTREE_ENERGY = Decimal("4.3597447222060e-18")
CALORIE_TO_JOULE = Decimal("4.184")
ATOMIC_TIME_UNIT = Decimal("2.4188843265857e-17")


@dataclass(frozen=True)
class UnitDefinition:
    key: str
    label: str
    to_base_factor: Decimal
    aliases: tuple[str, ...] = field(default_factory=tuple)


class EnergyConverter:
    """Backend service for unit conversion.

    Each unit family is normalized through a category-specific base unit to keep
    the logic simple and make additional unit families straightforward to add.
    """

    def __init__(self) -> None:
        self._categories: Dict[str, Dict[str, UnitDefinition]] = {
            "energy": self._build_energy_units(),
            "length": self._build_length_units(),
            "density": self._build_density_units(),
            "time": self._build_time_units(),
        }

    def _build_energy_units(self) -> Dict[str, UnitDefinition]:
        calorie_per_mole = CALORIE_TO_JOULE / AVOGADRO_CONSTANT
        kilocalorie_per_mole = Decimal("1000") * calorie_per_mole

        units = (
            UnitDefinition("hartree", "Hartree", HARTREE_ENERGY, aliases=("eh",)),
            UnitDefinition("joule", "Joule", Decimal("1"), aliases=("j",)),
            UnitDefinition(
                "cal/mol",
                "cal/mol",
                calorie_per_mole,
                aliases=("calorie/mol",),
            ),
            UnitDefinition(
                "kcal/mol",
                "kcal/mol",
                kilocalorie_per_mole,
                aliases=("kilocalorie/mol",),
            ),
            UnitDefinition(
                "ev",
                "Electronvolt (eV)",
                ELEMENTARY_CHARGE,
                aliases=("eV", "electronvolt"),
            ),
        )
        return {unit.key: unit for unit in units}

    def _build_length_units(self) -> Dict[str, UnitDefinition]:
        units = (
            UnitDefinition("meter", "Meter (m)", Decimal("1"), aliases=("m", "metre")),
            UnitDefinition(
                "kilometer",
                "Kilometer (km)",
                Decimal("1000"),
                aliases=("km", "kilometre"),
            ),
            UnitDefinition(
                "centimeter",
                "Centimeter (cm)",
                Decimal("0.01"),
                aliases=("cm", "centimetre"),
            ),
            UnitDefinition(
                "millimeter",
                "Millimeter (mm)",
                Decimal("0.001"),
                aliases=("mm", "millimetre"),
            ),
            UnitDefinition(
                "micrometer",
                "Micrometer (µm)",
                Decimal("1e-6"),
                aliases=("um", "µm", "micron"),
            ),
            UnitDefinition(
                "nanometer",
                "Nanometer (nm)",
                Decimal("1e-9"),
                aliases=("nm", "nanometre"),
            ),
            UnitDefinition(
                "angstrom",
                "Angstrom (Å)",
                Decimal("1e-10"),
                aliases=("angstrom", "ångström", "ang", "å"),
            ),
            UnitDefinition(
                "bohr",
                "Bohr unit (a₀)",
                Decimal("5.29177210903e-11"),
                aliases=("a0", "a_0", "a₀", "bohr radius"),
            ),
            UnitDefinition("mile", "Mile (mi)", Decimal("1609.344"), aliases=("mi",)),
            UnitDefinition("yard", "Yard (yd)", Decimal("0.9144"), aliases=("yd",)),
            UnitDefinition("foot", "Foot (ft)", Decimal("0.3048"), aliases=("ft",)),
            UnitDefinition("inch", "Inch (in)", Decimal("0.0254"), aliases=("in",)),
            UnitDefinition(
                "light_year",
                "Light Year (ly)",
                Decimal("9.4607304725808e15"),
                aliases=("light year", "ly", "light-year"),
            ),
            UnitDefinition(
                "parsec",
                "Parsec (pc)",
                Decimal("3.0856775814913673e16"),
                aliases=("pc",),
            ),
        )
        return {unit.key: unit for unit in units}

    def _build_time_units(self) -> Dict[str, UnitDefinition]:
        second = Decimal("1")
        minute = Decimal("60") * second
        hour = Decimal("60") * minute
        day = Decimal("24") * hour
        week = Decimal("7") * day
        year = Decimal("365") * day
        month = year / Decimal("12")

        units = (
            UnitDefinition("second", "Second (s)", second, aliases=("s", "sec")),
            UnitDefinition("millisecond", "Millisecond (ms)", Decimal("1e-3"), aliases=("ms",)),
            UnitDefinition("minute", "Minute (min)", minute, aliases=("min",)),
            UnitDefinition("hour", "Hour (h)", hour, aliases=("h", "hr")),
            UnitDefinition("day", "Day (d)", day, aliases=("d",)),
            UnitDefinition("week", "Week", week, aliases=("wk", "w")),
            UnitDefinition("month", "Month", month, aliases=("mo",)),
            UnitDefinition("year", "Year (y)", year, aliases=("y", "yr")),
            UnitDefinition("decade", "Decade", Decimal("10") * year, aliases=("dec",)),
            UnitDefinition("century", "Century", Decimal("100") * year, aliases=("c",)),
            UnitDefinition("millennium", "Millennium", Decimal("1000") * year, aliases=("mil",)),
            UnitDefinition("microsecond", "Microsecond (µs)", Decimal("1e-6"), aliases=("us", "µs")),
            UnitDefinition("nanosecond", "Nanosecond (ns)", Decimal("1e-9"), aliases=("ns",)),
            UnitDefinition("picosecond", "Picosecond (ps)", Decimal("1e-12"), aliases=("ps",)),
            UnitDefinition("femtosecond", "Femtosecond (fs)", Decimal("1e-15"), aliases=("fs",)),
            UnitDefinition(
                "atomic_time_unit",
                "Atomic time unit (a.u.)",
                ATOMIC_TIME_UNIT,
                aliases=("au", "a.u.", "atomic unit of time"),
            ),
        )
        return {unit.key: unit for unit in units}

    def _build_density_units(self) -> Dict[str, UnitDefinition]:
        amu_per_cubic_angstrom = (Decimal("1e24") / AVOGADRO_CONSTANT)

        units = (
            UnitDefinition(
                "amu/angstrom^3",
                "atomic mass unit per cubic ångström (amu/Å³)",
                amu_per_cubic_angstrom,
                aliases=("amu/å^3", "amu/a^3", "amu/ang^3", "amu/å³", "amu/a3"),
            ),
            UnitDefinition(
                "g/cm^3",
                "grams per cubic centimeter (g/cm³)",
                Decimal("1"),
                aliases=("g/cm3", "gram per cubic centimeter"),
            ),
            UnitDefinition(
                "kg/l",
                "kilograms per liter (kg/L)",
                Decimal("1"),
                aliases=("kg/liter", "kg/L"),
            ),
        )
        return {unit.key: unit for unit in units}

    def list_categories(self) -> list[str]:
        return sorted(self._categories.keys())

    def list_units(self, category: str = "energy") -> list[str]:
        self._require_category(category)
        return list(self._categories[category].keys())

    def get_display_map(self, category: str = "energy") -> Dict[str, str]:
        self._require_category(category)
        return {
            unit.key: unit.label for unit in self._categories[category].values()
        }

    def parse_value(self, raw_value: object) -> Decimal:
        if raw_value is None:
            raise ValueError("Please enter a value to convert.")

        text = str(raw_value).strip()
        if not text:
            raise ValueError("Please enter a value to convert.")

        try:
            return Decimal(text)
        except InvalidOperation as exc:
            raise ValueError(
                "Enter a valid number. Scientific notation is supported."
            ) from exc

    def convert(
        self,
        value: Decimal | int | float | str,
        from_unit: str,
        to_unit: str,
        category: str = "energy",
    ) -> Decimal:
        self._require_category(category)
        source = self._get_unit(category, from_unit)
        target = self._get_unit(category, to_unit)
        numeric_value = self.parse_value(value)

        with localcontext() as ctx:
            ctx.prec = 28
            value_in_base_units = numeric_value * source.to_base_factor
            return value_in_base_units / target.to_base_factor

    def format_result(self, value: Decimal | int | float | str) -> str:
        numeric_value = self.parse_value(value)
        if numeric_value.is_zero():
            return "0"

        magnitude = abs(numeric_value)
        if magnitude >= Decimal("1e6") or magnitude < Decimal("1e-4"):
            return f"{numeric_value:.10E}"

        normalized = numeric_value.normalize()
        text = format(normalized, "f")
        return text.rstrip("0").rstrip(".") if "." in text else text

    def describe_result(
        self,
        value: Decimal | int | float | str,
        from_unit: str,
        to_unit: str,
        category: str = "energy",
        use_labels: bool = False,
    ) -> str:
        source = self._get_unit(category, from_unit)
        target = self._get_unit(category, to_unit)
        converted = self.convert(value, from_unit, to_unit, category=category)
        source_text = source.label if use_labels else source.key
        target_text = target.label if use_labels else target.key
        return (
            f"{self.format_result(self.parse_value(value))} {source_text} = "
            f"{self.format_result(converted)} {target_text}"
        )

    def _require_category(self, category: str) -> None:
        if category not in self._categories:
            raise ValueError(f"Unsupported category: {category}")

    def _get_unit(self, category: str, unit_name: str) -> UnitDefinition:
        if unit_name in self._categories[category]:
            return self._categories[category][unit_name]

        lowered = unit_name.lower()
        for unit in self._categories[category].values():
            if lowered == unit.key.lower() or lowered in {
                alias.lower() for alias in unit.aliases
            }:
                return unit

        available = ", ".join(self.list_units(category))
        raise ValueError(f"Unsupported unit '{unit_name}'. Available units: {available}")

    def iter_unit_labels(self, category: str = "energy") -> Iterable[tuple[str, str]]:
        self._require_category(category)
        for unit in self._categories[category].values():
            yield unit.key, unit.label
