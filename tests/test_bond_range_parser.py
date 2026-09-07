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


# --------------------------------------------------------------------------- #
# sort_unique=False: order-preserving mode used by meanResidenceTime.py         #
# --------------------------------------------------------------------------- #
def test_default_still_sorts_and_deduplicates():
    """The three all-* tools rely on this; the default must not change."""
    assert parse_solute_index_ranges("5 1-3 2") == [1, 2, 3, 5]


def test_order_preserving_mode_keeps_input_order():
    assert parse_solute_index_ranges("9,1,5", sort_unique=False) == [9, 1, 5]


def test_order_preserving_mode_keeps_duplicates():
    assert parse_solute_index_ranges("3,3,1", sort_unique=False) == [3, 3, 1]


def test_order_preserving_mode_expands_ranges_in_place():
    assert parse_solute_index_ranges("10,3-5,1", sort_unique=False) == [10, 3, 4, 5, 1]


def test_both_modes_agree_on_an_already_sorted_unique_input():
    text = "3-7,10,15-17"
    assert (parse_solute_index_ranges(text)
            == parse_solute_index_ranges(text, sort_unique=False))


@pytest.mark.parametrize("text", ["7-3", "0-5", "3-", "-5", "3-5-7", "a-b"])
def test_validation_is_identical_in_both_modes(text):
    with pytest.raises(ValueError):
        parse_solute_index_ranges(text, sort_unique=False)
