"""The Toga-free periodic-table data (layout, blocks, groups and periods)."""
import subprocess
import sys
from pathlib import Path

import pytest

import periodicTableData as ptd


def test_module_does_not_import_toga():
    code = "import sys, periodicTableData; assert 'toga' not in sys.modules, 'toga imported'"
    subprocess.run([sys.executable, "-c", code], check=True, cwd=Path(__file__).resolve().parents[1])


def test_the_window_module_still_exports_the_data():
    import periodicTable

    assert periodicTable.PERIODIC_TABLE_LAYOUT is ptd.PERIODIC_TABLE_LAYOUT
    assert periodicTable.BLOCK_COLORS is ptd.BLOCK_COLORS
    assert periodicTable.element_block is ptd.element_block
    assert periodicTable.table_symbols is ptd.table_symbols


@pytest.mark.parametrize("symbol, expected", [
    ("H", (1, 1)), ("He", (18, 1)), ("Mg", (2, 3)), ("Ti", (4, 4)), ("Sn", (14, 5)),
    ("La", (3, 6)), ("Ce", (3, 6)), ("Lu", (3, 6)), ("Ac", (3, 7)), ("U", (3, 7)), ("Rn", (18, 6)),
])
def test_group_and_period(symbol, expected):
    assert ptd.group_and_period(symbol) == expected


def test_every_element_has_a_group_and_period():
    for symbol in ptd.table_symbols():
        group, period = ptd.group_and_period(symbol)
        assert 1 <= group <= 18 and 1 <= period <= 7


def test_an_unknown_symbol_is_a_key_error():
    with pytest.raises(KeyError):
        ptd.group_and_period("Xx")
