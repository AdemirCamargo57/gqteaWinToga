"""Tests for the solute-index range parser used by the All Bond Distance tool.

``parse_solute_index_ranges`` expands compact range syntax like
``"1-5 14-16 18 20"`` into an explicit, sorted, de-duplicated list of 1-based
indices. It is a pure function (no Toga), so it is unit-tested directly.
"""
import pytest

from allBondAnalysis import parse_solute_index_ranges


def test_expands_ranges_and_singletons():
    assert parse_solute_index_ranges("1-5 14-16 18 20") == [1, 2, 3, 4, 5, 14, 15, 16, 18, 20]


def test_plain_list_is_unchanged():
    assert parse_solute_index_ranges("1 2 3 4") == [1, 2, 3, 4]


def test_commas_are_accepted_as_separators():
    assert parse_solute_index_ranges("1-3, 7, 9-10") == [1, 2, 3, 7, 9, 10]


def test_overlaps_and_duplicates_are_deduplicated_and_sorted():
    assert parse_solute_index_ranges("3 1-3 2") == [1, 2, 3]


def test_single_value_range_is_allowed():
    assert parse_solute_index_ranges("5-5") == [5]


def test_empty_string_returns_empty_list():
    assert parse_solute_index_ranges("   ") == []


@pytest.mark.parametrize("text", [
    "5-2",        # descending range
    "0",          # not a positive (1-based) index
    "-3",         # negative
    "1--3",       # malformed
    "a-b",        # non-numeric
    "2-",         # missing bound
    "1.5",        # not an integer
])
def test_invalid_input_raises(text):
    with pytest.raises(ValueError):
        parse_solute_index_ranges(text)
