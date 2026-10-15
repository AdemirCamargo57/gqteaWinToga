import subprocess
import sys

import pytest

import systemProperties as sp


WATER = ["O", "H", "H"]


def test_water_electrons_and_singlet():
    props = sp.compute_properties(WATER * 10, net_charge=0, multiplicity=1)
    assert props.atom_count == 30
    assert props.electron_count == 100
    assert props.multiplicity_valid
    assert props.message == ""


def test_charge_removes_electrons():
    # Hydronium-like: 3 H + O with +1 -> 10 electrons.
    props = sp.compute_properties(["O", "H", "H", "H"], net_charge=1, multiplicity=1)
    assert props.electron_count == 10
    assert props.multiplicity_valid


def test_parity_mismatch_is_reported_not_corrected():
    props = sp.compute_properties(WATER, net_charge=0, multiplicity=2)
    assert not props.multiplicity_valid
    assert "1, 3, 5, 7" in props.message
    # The value the user gave is kept; nothing is silently replaced.
    assert props.multiplicity == 2


def test_odd_electron_count_allows_doublet_and_quartet():
    props = sp.compute_properties(["N", "O"], net_charge=0, multiplicity=2)  # NO radical, 15 e
    assert props.electron_count == 15
    assert props.multiplicity_valid
    assert sp.allowed_multiplicities(15) == [2, 4, 6, 8]


def test_triplet_is_allowed_for_even_count():
    # O2 is a triplet: the even count must not force a singlet.
    assert sp.compute_properties(["O", "O"], 0, 3).multiplicity_valid


def test_multiplicity_above_electron_bound():
    props = sp.compute_properties(["H"], net_charge=0, multiplicity=4)
    assert not props.multiplicity_valid
    assert "at most multiplicity 2" in props.message


def test_negative_electron_count():
    props = sp.compute_properties(["H"], net_charge=2, multiplicity=1)
    assert props.electron_count == -1
    assert not props.multiplicity_valid


def test_unknown_element_is_named():
    with pytest.raises(ValueError, match="Xx"):
        sp.compute_properties(["O", "Xx"], 0, 1)


def test_missing_multiplicity_is_not_inferred():
    props = sp.compute_properties(WATER, 0, None)
    assert not props.multiplicity_valid
    assert props.multiplicity is None
    with pytest.raises(ValueError, match="not inferred"):
        sp.parse_multiplicity("   ")


@pytest.mark.parametrize("text,expected", [("1", 1), (" 3 ", 3), ("2.0", 2)])
def test_parse_multiplicity_ok(text, expected):
    assert sp.parse_multiplicity(text) == expected


@pytest.mark.parametrize("text", ["0", "-1", "1.5", "abc"])
def test_parse_multiplicity_rejects(text):
    with pytest.raises(ValueError):
        sp.parse_multiplicity(text)


def test_parse_integer_default_and_errors():
    assert sp.parse_integer("", "Solute charge", default=0) == 0
    assert sp.parse_integer("-2", "Solute charge") == -2
    with pytest.raises(ValueError, match="Solute charge"):
        sp.parse_integer("0.5", "Solute charge")


def test_total_charge_per_species():
    assert sp.total_charge({"solute": 1, "solvent": 20}, {"solute": -1, "solvent": 0}) == -1
    assert sp.total_charge({"solute": 1, "A": 3, "B": 2}, {"solute": 2, "A": -1, "B": 1}) == 1


def test_summary_lines_mark_inconsistency():
    text = "".join(sp.compute_properties(WATER, 0, 2).summary_lines())
    assert "Total number of atoms: 3" in text
    assert "Total number of electrons: 10" in text
    assert "INCONSISTENT" in text


def test_module_does_not_import_toga():
    code = "import sys, systemProperties; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=str(__import__('pathlib').Path(__file__).resolve().parents[1]))
    assert out.stdout.strip() == "False", out.stderr
