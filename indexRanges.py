"""The project's single parser of compact atom-index ranges ("1-5 14-16 18").

Toga-free, so numerical modules can use it; allBondAnalysis re-exports it.
"""

from typing import Iterable, List


def parse_solute_index_ranges(text: str, sort_unique: bool = True) -> List[int]:
    """Expand compact atom-index syntax into an explicit 1-based list.

    Accepts whitespace- and/or comma-separated tokens, where each token is
    either a single positive integer (``18``) or an inclusive range
    (``14-16``), so ``"3-7,10,15-17"`` means ``[3,4,5,6,7,10,15,16,17]``. An
    empty/blank string returns ``[]``. Raises ``ValueError`` on non-positive
    values, descending ranges, or malformed tokens.

    ``sort_unique=True`` (the default, used by the three all-* connectivity
    tools) returns a sorted, de-duplicated list. ``sort_unique=False``
    preserves the order the user typed and keeps duplicates -- meanResidenceTime
    needs that, because it labels observed objects in input order.
    """
    ordered: List[int] = []
    for token in text.replace(",", " ").split():
        if "-" in token:
            parts = token.split("-")
            if len(parts) != 2 or not parts[0] or not parts[1]:
                raise ValueError(f"Invalid index range: '{token}'.")
            try:
                low, high = int(parts[0]), int(parts[1])
            except ValueError as exc:
                raise ValueError(f"Invalid index range: '{token}'.") from exc
            if low <= 0 or high <= 0:
                raise ValueError("Atom indices must be positive integers starting at 1.")
            if high < low:
                raise ValueError(f"Range bounds must be ascending: '{token}'.")
            ordered.extend(range(low, high + 1))
        else:
            try:
                value = int(token)
            except ValueError as exc:
                raise ValueError(f"Invalid atom index: '{token}'.") from exc
            if value <= 0:
                raise ValueError("Atom indices must be positive integers starting at 1.")
            ordered.append(value)

    return sorted(set(ordered)) if sort_unique else ordered


def format_index_ranges(indices: Iterable[int], sep: str = ", ") -> str:
    """Inverse of parse_solute_index_ranges: [1, 2, 3, 7, 9, 10] -> "1-3, 7, 9-10".

    Sorted and de-duplicated; ``sep`` joins the parts (", " for a text field,
    "," for a compact XYZ comment line). Empty input gives "".
    """
    parts: List[str] = []
    start = previous = None
    for value in sorted(set(int(v) for v in indices)):
        if start is None:
            start = previous = value
        elif value == previous + 1:
            previous = value
        else:
            parts.append(str(start) if start == previous else f"{start}-{previous}")
            start = previous = value
    if start is not None:
        parts.append(str(start) if start == previous else f"{start}-{previous}")
    return sep.join(parts)
