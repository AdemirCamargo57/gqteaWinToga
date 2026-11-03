"""The data tables of the Gaussian 16 input builder (gaussianKeywords, no Toga)."""
import pytest

import gaussianKeywords as gk

LABEL_TABLES = {
    "OPT_TARGETS": gk.OPT_TARGETS, "FORCE_CONSTANTS": gk.FORCE_CONSTANTS,
    "IRC_DIRECTIONS": gk.IRC_DIRECTIONS, "NMR_METHODS": gk.NMR_METHODS, "STATE_SPINS": gk.STATE_SPINS,
    "SPIN_PREFIXES": gk.SPIN_PREFIXES, "DIFFUSE": gk.DIFFUSE, "POLARIZATION": gk.POLARIZATION,
    "DISPERSION": gk.DISPERSION, "GEOMETRY_SOURCES": gk.GEOMETRY_SOURCES, "GUESS_METHODS": gk.GUESS_METHODS,
    "POP_LEVELS": gk.POP_LEVELS, "SOLVATION_MODELS": gk.SOLVATION_MODELS,
    **{f"FAMILY_METHODS[{family}]": table for family, table in gk.FAMILY_METHODS.items()},
}


@pytest.mark.parametrize("name", sorted(LABEL_TABLES))
def test_route_values_are_unique_so_a_route_reads_back_to_one_label(name):
    values = [v.upper() for v in LABEL_TABLES[name].values() if v]
    assert len(values) == len(set(values)), name


@pytest.mark.parametrize("name", sorted(LABEL_TABLES))
def test_reverse_maps_every_route_value_to_its_label(name):
    table = LABEL_TABLES[name]
    back = gk.reverse(table)
    for label, value in table.items():
        if value:
            assert back[value.upper()] == label


def test_choice_lists_have_unique_entries():
    for values in (gk.JOB_TYPES, gk.STATES, gk.FAMILIES, gk.BASES, gk.SOLVENTS, gk.PRINT_LEVELS):
        assert len(values) == len(set(values))


def test_every_family_has_a_route_form():
    for family in gk.FAMILIES:
        assert (family in gk.FAMILY_KEYWORDS) != (family in gk.FAMILY_METHODS), family


def test_default_method_of_each_choice_family_is_one_of_its_methods():
    for family, method in gk.DEFAULT_METHOD.items():
        assert method in gk.FAMILY_METHODS[family]


def test_keywords_gaussian_spells_differently_from_the_label():
    assert gk.FAMILY_METHODS[gk.FAMILY_DFT]["M06-2X"] == "M062X"
    assert gk.FAMILY_METHODS[gk.FAMILY_DFT]["wB97X-D"] == "wB97XD"
    assert gk.FAMILY_METHODS[gk.FAMILY_DFT]["PBE0 (PBE1PBE)"] == "PBE1PBE"


def test_solvent_names_hold_no_comma_or_blank_because_they_become_one_route_option():
    for solvent in gk.SOLVENTS:
        assert "," not in solvent and " " not in solvent


def test_fukui_force_constants_exclude_reading_the_checkpoint():
    assert set(gk.FUKUI_FORCE_CONSTANTS) < set(gk.FORCE_CONSTANTS)
    assert not any(gk.FORCE_CONSTANTS[label] in gk.READ_FORCE_CONSTANTS for label in gk.FUKUI_FORCE_CONSTANTS)


def test_flag_fields_name_route_options():
    assert gk.FREQ_FLAG_FIELDS == {"freq_raman": "Raman", "freq_vcd": "VCD", "freq_anharmonic": "Anharmonic"}
    assert set(gk.GUESS_FLAG_FIELDS.values()) == {"Mix", "Always", "Only", "Save"}
    assert set(gk.POP_FLAG_FIELDS.values()) == {"MK", "Hirshfeld", "NBO", "NPA", "SaveNTO"}


def test_every_panel_keyword_has_an_owner():
    for key in ("OPT", "FREQ", "IRC", "STABLE", "NMR", "TD", "TDA", "EMPIRICALDISPERSION", "IOP", "SCF",
                "NOSYMM", "POLAR", "GEOM", "GUESS", "POP", "SCRF"):
        assert key in gk.PANEL_OF_KEYWORD
