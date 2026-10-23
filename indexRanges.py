"""The project's single parser of compact atom-index ranges ("1-5 14-16 18").

Toga-free, so numerical modules can use it; allBondAnalysis re-exports it.
"""

from typing import List


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
