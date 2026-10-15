"""Atom count, electron count and spin-multiplicity checks for a built box.

Shared by the solvation-box builders. Pure Python, no Toga: every problem is a
``ValueError`` (or an ``inconsistent`` result) with a message the UI shows as a
dialog.

The multiplicity is **never inferred**. An even electron count allows a singlet,
a triplet, a quintet...; an odd one a doublet, a quartet... -- the count alone
does not single one out, so the user always supplies it and this module only
checks that it is possible: ``(N_e + M)`` odd and ``1 <= M <= N_e + 1``.
"""
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence

from help import AtomicData


@dataclass
class SystemProperties:
    atom_count: int
    nuclear_charge: int          # sum of atomic numbers
    net_charge: int
    electron_count: int
    multiplicity: Optional[int]
    multiplicity_valid: bool
    message: str                 # empty when the multiplicity is valid

    def summary_lines(self) -> List[str]:
        if self.multiplicity is None:
            mult = "not set"
        elif self.multiplicity_valid:
            mult = str(self.multiplicity)
        else:
            mult = f"{self.multiplicity} (INCONSISTENT: {self.message})"
        return [
            "--- System properties ---\n",
            f"Total number of atoms: {self.atom_count}\n",
            f"Net charge: {self.net_charge:+d}\n",
            f"Total number of electrons: {self.electron_count}\n",
            f"Total spin multiplicity: {mult}\n",
        ]


def parse_integer(text: str, field_name: str, default: Optional[int] = None) -> Optional[int]:
    """Parse an integer field; blank gives ``default``. ``1.0`` is accepted, ``1.5`` is not."""
    value = (text or "").strip()
    if not value:
        return default
    try:
        number = float(value)
    except ValueError:
        raise ValueError(f"{field_name} must be an integer (got '{value}').")
    if not number.is_integer():
        raise ValueError(f"{field_name} must be an integer (got '{value}').")
    return int(number)


def parse_multiplicity(text: str) -> int:
    """The total spin multiplicity field: required, integer, >= 1."""
    value = parse_integer(text, "Total spin multiplicity")
    if value is None:
        raise ValueError(
            "Please enter the total spin multiplicity (2S+1). It is not inferred from the "
            "electron count: an even count allows 1, 3, 5, ... and an odd count 2, 4, 6, ..."
        )
    if value < 1:
        raise ValueError("Total spin multiplicity must be a positive integer (2S+1 >= 1).")
    return value


def nuclear_charge(elements: Iterable[str]) -> int:
    total = 0
    for element in elements:
        number = AtomicData.atomic_numbers.get(element)
        if number is None:
            raise ValueError(f"Unknown element '{element}': no atomic number is available.")
        total += int(number)
    return total


def allowed_multiplicities(electron_count: int, limit: int = 4) -> List[int]:
    """The first ``limit`` multiplicities compatible with ``electron_count``."""
    if electron_count < 0:
        return []
    first = 1 if electron_count % 2 == 0 else 2
    values = list(range(first, electron_count + 2, 2))
    return values[:limit]


def multiplicity_problem(electron_count: int, multiplicity: int) -> str:
    """Empty string when ``multiplicity`` is possible for ``electron_count``."""
    if electron_count < 0:
        return f"the net charge leaves {electron_count} electrons"
    if multiplicity < 1:
        return "the multiplicity must be at least 1"
    if multiplicity > electron_count + 1:
        return f"{electron_count} electrons allow at most multiplicity {electron_count + 1}"
    if (electron_count + multiplicity) % 2 == 0:
        parity = "even" if electron_count % 2 == 0 else "odd"
        allowed = ", ".join(str(m) for m in allowed_multiplicities(electron_count))
        return f"{electron_count} electrons ({parity}) allow multiplicity {allowed}, ..."
    return ""


def total_charge(species_counts: Dict[str, int], species_charges: Dict[str, int]) -> int:
    """Net charge = sum over species of (number of molecules x charge per molecule)."""
    return int(sum(int(count) * int(species_charges.get(name, 0)) for name, count in species_counts.items()))


def compute_properties(elements: Sequence[str], net_charge: int,
                       multiplicity: Optional[int]) -> SystemProperties:
    z_total = nuclear_charge(elements)
    electrons = z_total - int(net_charge)
    if multiplicity is None:
        valid, message = False, "no multiplicity was given"
    else:
        message = multiplicity_problem(electrons, int(multiplicity))
        valid = not message
    return SystemProperties(
        atom_count=len(elements),
        nuclear_charge=z_total,
        net_charge=int(net_charge),
        electron_count=electrons,
        multiplicity=multiplicity,
        multiplicity_valid=valid,
        message=message,
    )
