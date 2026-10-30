"""The periodic table as data: layout, blocks, groups and periods.

No Toga here (a subprocess test pins it), so numerical modules such as
runatomDefaults can use the layout. periodicTable re-exports these names for
the picker window and its existing importers.
"""
from typing import List, Optional, Tuple

# The periodic table as it is drawn: 18 columns per row, ``None`` for a gap.
# The two f-block rows sit below a blank spacer row, the way a printed table
# lays them out. A test asserts this contains all 118 elements exactly once.
PERIODIC_TABLE_LAYOUT: List[List[Optional[str]]] = [
    ["H", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, "He"],
    ["Li", "Be", None, None, None, None, None, None, None, None, None, None, "B", "C", "N", "O", "F", "Ne"],
    ["Na", "Mg", None, None, None, None, None, None, None, None, None, None, "Al", "Si", "P", "S", "Cl", "Ar"],
    ["K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"],
    ["Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"],
    ["Cs", "Ba", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"],
    ["Fr", "Ra", "Ac", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"],
    [None] * 18,
    [None, None, None, "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", None],
    [None, None, None, "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", None],
]
F_BLOCK_FIRST_ROW = 8

# Block colours, so the table reads as a periodic table rather than a wall of
# identical squares.
BLOCK_COLORS = {
    "s": "#ffd9d9",
    "p": "#d9e8ff",
    "d": "#e2f5df",
    "f": "#f3e2ff",
}
SELECTED_ELEMENT_COLOR = "#ffe066"


def table_symbols() -> List[str]:
    """Every element in the layout, in reading order."""
    return [symbol for row in PERIODIC_TABLE_LAYOUT for symbol in row if symbol]


def element_block(symbol: str) -> str:
    """Which block of the periodic table an element belongs to."""
    for row_index, row in enumerate(PERIODIC_TABLE_LAYOUT):
        if symbol not in row:
            continue
        column = row.index(symbol)
        if row_index >= F_BLOCK_FIRST_ROW:
            return "f"
        if column <= 1 or symbol == "He":
            return "s"
        if column >= 12:
            return "p"
        return "d"
    return "p"


def group_and_period(symbol: str) -> Tuple[int, int]:
    """(IUPAC group 1-18, period); lanthanides and actinides count as group 3 (the La/Ac column)."""
    for row_index, row in enumerate(PERIODIC_TABLE_LAYOUT):
        if symbol in row:
            if row_index >= F_BLOCK_FIRST_ROW:
                return 3, row_index - F_BLOCK_FIRST_ROW + 6
            return row.index(symbol) + 1, row_index + 1
    raise KeyError(symbol)
