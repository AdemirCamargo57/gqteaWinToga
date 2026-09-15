"""Tests for the Molecular Geometry Comparator tool (``MolGeomComparator.py``).

Run from ``venv/src/``::

    python -m pytest tests/ -q

Only the pure helpers and the logic class ``MolGeomCalculator`` are exercised
(no Toga GUI is created), following the same pattern as the other tool tests.

The tool compares the combined parameter files written by ``allBondAnalysis``,
``allAnglesAnalysis`` and ``allDihedralAnalysis`` for the *same* molecule
simulated twice (typically isolated vs. in a water box). The matching problem is
an index-remapping problem: in the solvated run the solute atoms carry their
global indices in the box (e.g. the oxygens at 295-297), so every row is
canonicalised to a *solute-local* index via the ``# solute_atom_indices`` header
before the two files are matched.
"""
import asyncio
import os

import pytest

from MolGeomComparator import (
    BOND_KIND,
    ANGLE_KIND,
    DIHEDRAL_KIND,
    PERCENT_DIFFERENCE_MODE,
    PERCENT_ERROR_MODE,
    PERCENT_MODES,
    MolGeomCalculator,
    MolGeomComparatorUI,
    angular_difference,
    canonical_identity,
    circular_mean,
    default_output_name,
    detect_parameter_kind,
    parse_optional_fraction,
    percent_difference,
    percent_error,
    read_parameter_file,
    resolve_output_dir,
    sanitize_label,
)


# --------------------------------------------------------------------------- #
# Fixture builders: write files byte-compatible with the three all-* writers    #
# --------------------------------------------------------------------------- #
BOND_COLUMNS = (
    "row atom_i atom_j element_i element_j first_bonded_distance "
    "average variance standard_deviation occurrence_fraction frames_bonded"
)
ANGLE_COLUMNS = (
    "row atom_i atom_j atom_k element_i element_j element_k "
    "first_present_angle_degrees average_angle_degrees variance_degrees2 "
    "standard_deviation_degrees occurrence_fraction frames_present"
)
DIHEDRAL_COLUMNS = (
    "row atom_i atom_j atom_k atom_l element_i element_j element_k element_l "
    "first_present_dihedral_degrees average_dihedral_degrees variance_degrees2 "
    "standard_deviation_degrees occurrence_fraction frames_present"
)


def write_parameter_file(
    path,
    kind_title,
    columns,
    rows,
    *,
    frames=1,
    max_distance=1.7,
    solute=None,
    cell=None,
    extra_header=(),
):
    """Write a combined parameter file in the exact layout of the all-* tools.

    ``rows`` holds ``(atoms, elements, average, std_dev, occurrence)`` tuples;
    the first-value and variance columns are filled with plausible values since
    the comparator never reads them.
    """
    with open(path, "w") as out:
        out.write(f"# {kind_title}\n")
        out.write(f"# frames_used {frames}\n")
        out.write(f"# max_connection_distance {max_distance:g}\n")
        for line in extra_header:
            out.write(f"# {line}\n")
        if solute is None:
            out.write("# atom_scope all_atoms\n")
        else:
            out.write("# atom_scope solute_atoms\n")
            out.write("# solute_atom_indices " + " ".join(str(i) for i in solute) + "\n")
        if cell is None:
            out.write("# periodic_boundary none\n")
        else:
            out.write("# cell_lengths {:g} {:g} {:g}\n".format(*cell))
        out.write(f"# {columns}\n")
        for row_index, (atoms, elements, average, std_dev, occurrence) in enumerate(rows, start=1):
            fields = [f"{row_index:>6d}"]
            fields += [f"{atom:>8d}" for atom in atoms]
            fields += [f"{element:>8s}" for element in elements]
            variance = std_dev * std_dev
            fields.append(f"{average:>16.8f}")          # first_* column
            fields.append(f"{average:>16.8f}")          # average column
            fields.append(f"{variance:>16.8f}")
            fields.append(f"{std_dev:>16.8f}")
            fields.append(f"{occurrence:>16.8f}")
            fields.append(f"{int(round(occurrence * frames)):>10d}")
            out.write(" ".join(fields) + "\n")
    return str(path)


def write_bond_file(path, rows, **kwargs):
    return write_parameter_file(
        path, "gQTEA All Bond Distance Analysis", BOND_COLUMNS, rows, **kwargs
    )


def write_angle_file(path, rows, **kwargs):
    return write_parameter_file(
        path, "gQTEA All Bond Angle Analysis", ANGLE_COLUMNS, rows, **kwargs
    )


def write_dihedral_file(path, rows, **kwargs):
    kwargs.setdefault("extra_header", ("dihedral_convention signed_degrees_minus180_to_180",))
    return write_parameter_file(
        path, "gQTEA All Dihedral Angle Analysis", DIHEDRAL_COLUMNS, rows, **kwargs
    )


def bond_row(i, j, ei, ej, average, std_dev=0.0, occurrence=1.0):
    return ((i, j), (ei, ej), average, std_dev, occurrence)


def angle_row(i, j, k, ei, ej, ek, average, std_dev=0.0, occurrence=1.0):
    return ((i, j, k), (ei, ej, ek), average, std_dev, occurrence)


def dihedral_row(i, j, k, l, ei, ej, ek, el, average, std_dev=0.0, occurrence=1.0):
    return ((i, j, k, l), (ei, ej, ek, el), average, std_dev, occurrence)


# --------------------------------------------------------------------------- #
# Parameter-kind detection                                                      #
# --------------------------------------------------------------------------- #
def test_detect_parameter_kind_recognises_the_three_titles():
    assert detect_parameter_kind("# gQTEA All Bond Distance Analysis") is BOND_KIND
    assert detect_parameter_kind("# gQTEA All Bond Angle Analysis") is ANGLE_KIND
    assert detect_parameter_kind("# gQTEA All Dihedral Angle Analysis") is DIHEDRAL_KIND


def test_detect_parameter_kind_rejects_an_unknown_title():
    with pytest.raises(ValueError) as excinfo:
        detect_parameter_kind("# gQTEA Radial Distribution Function")

    assert "Radial Distribution Function" in str(excinfo.value)


def test_bond_kind_uses_angstrom_and_the_angle_kinds_use_degrees():
    assert BOND_KIND.unit == "angstrom"
    assert ANGLE_KIND.unit == "degrees"
    assert DIHEDRAL_KIND.unit == "degrees"
    assert BOND_KIND.n_atoms == 2
    assert ANGLE_KIND.n_atoms == 3
    assert DIHEDRAL_KIND.n_atoms == 4


def test_only_the_dihedral_kind_is_periodic():
    assert not BOND_KIND.periodic
    assert not ANGLE_KIND.periodic
    assert DIHEDRAL_KIND.periodic


# --------------------------------------------------------------------------- #
# Reading a file                                                                #
# --------------------------------------------------------------------------- #
def test_read_parameter_file_reads_rows_metadata_and_kind(tmp_path):
    path = write_bond_file(
        tmp_path / "iso.txt",
        [bond_row(1, 2, "C", "C", 1.54834003), bond_row(1, 3, "C", "O", 1.20555185)],
        frames=1,
    )

    parsed = read_parameter_file(path)

    assert parsed.kind is BOND_KIND
    assert parsed.frames_used == 1
    assert parsed.atom_scope == "all_atoms"
    assert len(parsed.rows) == 2
    assert parsed.rows[0].global_atoms == (1, 2)
    assert parsed.rows[0].elements == ("C", "C")
    assert parsed.rows[0].average == pytest.approx(1.54834003)
    assert parsed.rows[0].source_row == 1
    assert parsed.rows[1].source_row == 2


def test_read_parameter_file_keeps_std_dev_and_occurrence(tmp_path):
    path = write_bond_file(
        tmp_path / "solv.txt",
        [bond_row(1, 2, "C", "C", 1.52968232, std_dev=0.03377609, occurrence=0.5)],
        frames=100,
    )

    row = read_parameter_file(path).rows[0]

    assert row.std_dev == pytest.approx(0.03377609)
    assert row.occurrence == pytest.approx(0.5)
    assert row.frames == 50


def test_all_atoms_scope_leaves_indices_unchanged(tmp_path):
    path = write_bond_file(tmp_path / "iso.txt", [bond_row(1, 29, "C", "O", 1.2)])

    parsed = read_parameter_file(path)

    assert parsed.rows[0].local_atoms == (1, 29)
    assert parsed.num_local_atoms == 29


def test_solute_scope_remaps_global_indices_to_solute_local_indices(tmp_path):
    """The real case: the solute's oxygens sit at 295-297 in the water box."""
    solute = list(range(1, 29)) + [295, 296, 297]
    path = write_bond_file(
        tmp_path / "solv.txt",
        [
            bond_row(1, 295, "C", "O", 1.23841101),
            bond_row(9, 296, "C", "O", 1.38216013),
            bond_row(28, 297, "H", "O", 1.03381182),
        ],
        solute=solute,
    )

    parsed = read_parameter_file(path)

    assert parsed.atom_scope == "solute_atoms"
    assert parsed.num_local_atoms == 31
    assert parsed.rows[0].local_atoms == (1, 29)
    assert parsed.rows[1].local_atoms == (9, 30)
    assert parsed.rows[2].local_atoms == (28, 31)
    assert parsed.local_to_global[29] == 295
    assert parsed.local_to_global[31] == 297


def test_solute_scope_remap_is_independent_of_the_listed_order(tmp_path):
    """The header list is sorted before ranking, as the all-* writers sort it."""
    path = write_bond_file(
        tmp_path / "solv.txt",
        [bond_row(295, 1, "O", "C", 1.2)],
        solute=[295, 1, 296],
    )

    assert read_parameter_file(path).rows[0].local_atoms == (2, 1)


def test_element_table_is_built_from_the_rows(tmp_path):
    path = write_bond_file(
        tmp_path / "iso.txt",
        [bond_row(1, 2, "C", "C", 1.5), bond_row(2, 3, "C", "O", 1.4)],
    )

    elements = read_parameter_file(path).elements_by_local

    assert elements == {1: "C", 2: "C", 3: "O"}


def test_angle_and_dihedral_files_read_three_and_four_atoms(tmp_path):
    angle_path = write_angle_file(
        tmp_path / "ang.txt", [angle_row(1, 2, 3, "H", "C", "H", 109.47)]
    )
    dihedral_path = write_dihedral_file(
        tmp_path / "dih.txt", [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", -59.5)]
    )

    angle = read_parameter_file(angle_path)
    dihedral = read_parameter_file(dihedral_path)

    assert angle.kind is ANGLE_KIND
    assert angle.rows[0].local_atoms == (1, 2, 3)
    assert angle.rows[0].average == pytest.approx(109.47)
    assert dihedral.kind is DIHEDRAL_KIND
    assert dihedral.rows[0].local_atoms == (1, 2, 3, 4)
    assert dihedral.rows[0].average == pytest.approx(-59.5)


def test_columns_are_located_by_name_not_by_position(tmp_path):
    """A future column reorder must not silently shift the parsed values."""
    path = tmp_path / "reordered.txt"
    with open(path, "w") as out:
        out.write("# gQTEA All Bond Distance Analysis\n")
        out.write("# frames_used 10\n")
        out.write("# atom_scope all_atoms\n")
        out.write(
            "# row occurrence_fraction average atom_i atom_j element_i element_j "
            "standard_deviation frames_bonded\n"
        )
        out.write("     1 1.00000000 1.50000000        1        2        C        C "
                  "0.01000000         10\n")

    row = read_parameter_file(str(path)).rows[0]

    assert row.global_atoms == (1, 2)
    assert row.average == pytest.approx(1.5)
    assert row.std_dev == pytest.approx(0.01)


# --------------------------------------------------------------------------- #
# Reading: rejection of malformed files                                         #
# --------------------------------------------------------------------------- #
def test_read_rejects_a_file_that_is_not_a_gqtea_parameter_file(tmp_path):
    path = tmp_path / "random.txt"
    path.write_text("31\nsome xyz comment\nC 0.0 0.0 0.0\n")

    with pytest.raises(ValueError) as excinfo:
        read_parameter_file(str(path))

    assert "random.txt" in str(excinfo.value)


def test_read_rejects_a_file_without_a_row_column_header(tmp_path):
    path = tmp_path / "noheader.txt"
    path.write_text(
        "# gQTEA All Bond Distance Analysis\n"
        "# frames_used 1\n"
        "# atom_scope all_atoms\n"
        "     1        1        2        C        C 1.5 1.5 0.0 0.0 1.0 1\n"
    )

    with pytest.raises(ValueError) as excinfo:
        read_parameter_file(str(path))

    assert "# row" in str(excinfo.value)


def test_read_rejects_a_missing_required_column(tmp_path):
    path = tmp_path / "nocolumn.txt"
    path.write_text(
        "# gQTEA All Bond Distance Analysis\n"
        "# frames_used 1\n"
        "# atom_scope all_atoms\n"
        "# row atom_i atom_j element_i element_j standard_deviation occurrence_fraction\n"
        "     1        1        2        C        C 0.0 1.0\n"
    )

    with pytest.raises(ValueError) as excinfo:
        read_parameter_file(str(path))

    assert "average" in str(excinfo.value)


def test_read_rejects_a_truncated_data_row_naming_its_line(tmp_path):
    path = tmp_path / "short.txt"
    good = write_bond_file(tmp_path / "good.txt", [bond_row(1, 2, "C", "C", 1.5)])
    lines = open(good).read().splitlines()
    lines.append("     2        1")
    path.write_text("\n".join(lines) + "\n")

    with pytest.raises(ValueError) as excinfo:
        read_parameter_file(str(path))

    assert "line" in str(excinfo.value).lower()


def test_read_rejects_a_row_whose_atom_is_outside_the_solute_scope(tmp_path):
    path = write_bond_file(
        tmp_path / "solv.txt",
        [bond_row(1, 400, "C", "O", 1.2)],
        solute=[1, 2, 295],
    )

    with pytest.raises(ValueError) as excinfo:
        read_parameter_file(path)

    assert "400" in str(excinfo.value)


def test_read_rejects_a_file_that_gives_one_atom_two_elements(tmp_path):
    path = write_bond_file(
        tmp_path / "corrupt.txt",
        [bond_row(1, 2, "C", "C", 1.5), bond_row(2, 3, "N", "O", 1.4)],
    )

    with pytest.raises(ValueError) as excinfo:
        read_parameter_file(path)

    message = str(excinfo.value)
    assert "2" in message and "C" in message and "N" in message


def test_read_rejects_a_file_with_no_data_rows(tmp_path):
    path = write_bond_file(tmp_path / "empty.txt", [])

    with pytest.raises(ValueError) as excinfo:
        read_parameter_file(path)

    assert "no parameter rows" in str(excinfo.value).lower()


def test_read_tolerates_blank_lines_and_a_trailing_newline(tmp_path):
    good = write_bond_file(tmp_path / "good.txt", [bond_row(1, 2, "C", "C", 1.5)])
    path = tmp_path / "blanks.txt"
    path.write_text(open(good).read() + "\n\n")

    assert len(read_parameter_file(str(path)).rows) == 1


def test_solute_scope_counts_every_listed_atom_even_if_unbonded(tmp_path):
    """A solute atom that never appears in a row still belongs to the molecule."""
    path = write_bond_file(
        tmp_path / "solv.txt",
        [bond_row(1, 295, "C", "O", 1.2)],
        solute=[1, 2, 3, 295],
    )

    assert read_parameter_file(path).num_local_atoms == 4


# --------------------------------------------------------------------------- #
# Pure helpers: canonical identity, circular difference, naming                 #
# --------------------------------------------------------------------------- #
def test_bond_identity_ignores_the_order_of_the_two_atoms():
    assert canonical_identity(BOND_KIND, (2, 1)) == canonical_identity(BOND_KIND, (1, 2))


def test_angle_identity_keeps_the_vertex_and_sorts_the_arms():
    assert canonical_identity(ANGLE_KIND, (3, 2, 1)) == canonical_identity(ANGLE_KIND, (1, 2, 3))


def test_angle_identity_distinguishes_a_different_vertex():
    assert canonical_identity(ANGLE_KIND, (1, 2, 3)) != canonical_identity(ANGLE_KIND, (2, 1, 3))


def test_dihedral_identity_ignores_reversal_of_the_chain():
    """Reversing i-j-k-l leaves the signed dihedral unchanged, so it is one parameter."""
    assert canonical_identity(DIHEDRAL_KIND, (4, 3, 2, 1)) == canonical_identity(
        DIHEDRAL_KIND, (1, 2, 3, 4)
    )


def test_dihedral_identity_distinguishes_a_different_chain():
    assert canonical_identity(DIHEDRAL_KIND, (1, 2, 3, 4)) != canonical_identity(
        DIHEDRAL_KIND, (1, 3, 2, 4)
    )


def test_angular_difference_wraps_across_the_plus_minus_180_seam():
    """+179 deg to -179 deg is a 2 deg shift, not -358 deg."""
    assert angular_difference(179.0, -179.0) == pytest.approx(2.0)
    assert angular_difference(-179.0, 179.0) == pytest.approx(-2.0)


def test_angular_difference_is_plain_subtraction_away_from_the_seam():
    assert angular_difference(10.0, 20.0) == pytest.approx(10.0)
    assert angular_difference(20.0, 10.0) == pytest.approx(-10.0)


def test_angular_difference_reports_a_half_turn_as_positive_180():
    assert angular_difference(0.0, 180.0) == pytest.approx(180.0)


# --------------------------------------------------------------------------- #
# Relative-difference formulas                                                  #
# --------------------------------------------------------------------------- #
def test_percent_difference_is_signed_and_relative_to_the_mean():
    # mean 1.5, difference +0.1 -> 100 * 0.1 / 1.5
    assert percent_difference(BOND_KIND, 1.45, 1.55) == pytest.approx(100 * 0.1 / 1.5)
    assert percent_difference(BOND_KIND, 1.55, 1.45) == pytest.approx(-100 * 0.1 / 1.5)


def test_percent_difference_is_symmetric_in_magnitude():
    """Swapping the two files must only flip the sign -- that is the point of
    the symmetric formula: neither simulation is privileged."""
    forward = percent_difference(BOND_KIND, 1.4321, 1.5987)
    backward = percent_difference(BOND_KIND, 1.5987, 1.4321)
    assert forward == pytest.approx(-backward)


def test_percent_error_uses_file_1_as_the_reference():
    # 100 * (1.55 - 1.45) / 1.45
    assert percent_error(BOND_KIND, 1.45, 1.55) == pytest.approx(100 * 0.1 / 1.45)


def test_percent_error_is_not_symmetric():
    """Unlike the symmetric formula, swapping the files changes the magnitude."""
    forward = percent_error(BOND_KIND, 1.45, 1.55)
    backward = percent_error(BOND_KIND, 1.55, 1.45)
    assert abs(forward) != pytest.approx(abs(backward))


def test_dihedral_percentages_use_the_wrapped_numerator():
    """+179 -> -179 is a +2 deg shift, not -358 deg."""
    value = percent_difference(DIHEDRAL_KIND, 179.0, -179.0)
    # numerator +2, circular-mean denominator 180
    assert value == pytest.approx(100 * 2.0 / 180.0)


def test_circular_mean_takes_the_midpoint_on_the_circle():
    assert circular_mean(179.0, angular_difference(179.0, -179.0)) == pytest.approx(180.0)
    assert circular_mean(10.0, angular_difference(10.0, 20.0)) == pytest.approx(15.0)


def test_percent_is_undefined_when_the_denominator_falls_below_the_floor():
    import math
    # A dihedral averaging a few tenths of a degree: floor is 1.0 deg.
    assert math.isnan(percent_difference(DIHEDRAL_KIND, 0.4, 0.8))
    assert math.isnan(percent_error(DIHEDRAL_KIND, 0.4, 0.8))
    # A bond never trips its 1e-6 A floor.
    assert not math.isnan(percent_difference(BOND_KIND, 1.4, 1.5))
    # A normal angle does not trip the 1.0 deg floor either.
    assert not math.isnan(percent_difference(ANGLE_KIND, 109.5, 111.0))


def test_each_kind_carries_its_own_percent_floor():
    assert BOND_KIND.percent_floor == pytest.approx(1e-6)
    assert ANGLE_KIND.percent_floor == pytest.approx(1.0)
    assert DIHEDRAL_KIND.percent_floor == pytest.approx(1.0)


def test_percent_modes_have_a_fixed_order():
    """Column order must not depend on which switch the user ticked first."""
    assert PERCENT_MODES == (PERCENT_DIFFERENCE_MODE, PERCENT_ERROR_MODE)


def test_sanitize_label_makes_a_column_safe_token():
    assert sanitize_label("isolated") == "isolated"
    assert sanitize_label("in water box") == "in_water_box"
    assert sanitize_label("gas-phase (300 K)") == "gas-phase_300_K"


def test_sanitize_label_falls_back_when_nothing_usable_remains():
    assert sanitize_label("   ", fallback="file_1") == "file_1"
    assert sanitize_label("!!!", fallback="file_2") == "file_2"


def test_default_output_name_follows_the_parameter_type():
    assert default_output_name(BOND_KIND) == "geometry_comparison_bond.txt"
    assert default_output_name(ANGLE_KIND) == "geometry_comparison_angle.txt"
    assert default_output_name(DIHEDRAL_KIND) == "geometry_comparison_dihedral.txt"


# --------------------------------------------------------------------------- #
# Comparison                                                                    #
# --------------------------------------------------------------------------- #
def make_calculator(tmp_path, rows_1, rows_2, *, writer=write_bond_file, **kwargs):
    """Build a calculator over two freshly written files; defaults for the rest."""
    file_1 = writer(tmp_path / "file_1.txt", rows_1, **kwargs.pop("kwargs_1", {}))
    file_2 = writer(tmp_path / "file_2.txt", rows_2, **kwargs.pop("kwargs_2", {}))
    return MolGeomCalculator(file_1, file_2, **kwargs)


def test_matched_parameters_carry_both_averages_and_their_difference(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.54834003)],
        [bond_row(1, 2, "C", "C", 1.52968232)],
    )

    calculator.run()

    assert len(calculator.matched) == 1
    match = calculator.matched[0]
    assert match.average_1 == pytest.approx(1.54834003)
    assert match.average_2 == pytest.approx(1.52968232)
    assert match.difference == pytest.approx(-0.01865771)


def test_matching_ignores_the_atom_order_within_a_row(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.50)],
        [bond_row(2, 1, "C", "C", 1.55)],
    )

    calculator.run()

    assert len(calculator.matched) == 1
    assert calculator.matched[0].difference == pytest.approx(0.05)


def test_matching_remaps_solvated_global_indices_onto_the_isolated_molecule(tmp_path):
    """The real case: C1-O295 in the box is the same bond as C1-O29 isolated."""
    solute = list(range(1, 29)) + [295, 296, 297]
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 29, "C", "O", 1.20555185), bond_row(28, 31, "H", "O", 0.97231769)],
        [bond_row(1, 295, "C", "O", 1.23841101), bond_row(28, 297, "H", "O", 1.03381182)],
        kwargs_2={"solute": solute, "frames": 15121},
    )

    calculator.run()

    assert len(calculator.matched) == 2
    assert calculator.only_in_file_1 == []
    assert calculator.only_in_file_2 == []
    assert calculator.matched[0].difference == pytest.approx(0.03285916)


def test_unmatched_parameters_are_partitioned_by_source_file(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(1, 3, "C", "O", 1.4)],
        [bond_row(1, 2, "C", "C", 1.5), bond_row(2, 3, "C", "O", 1.3)],
    )

    calculator.run()

    assert len(calculator.matched) == 1
    assert [row.local_atoms for row in calculator.only_in_file_1] == [(1, 3)]
    assert [row.local_atoms for row in calculator.only_in_file_2] == [(2, 3)]


def test_matched_parameters_are_sorted_by_atom_index(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(2, 3, "C", "O", 1.4), bond_row(1, 2, "C", "C", 1.5)],
        [bond_row(1, 2, "C", "C", 1.5), bond_row(2, 3, "C", "O", 1.4)],
    )

    calculator.run()

    assert [match.local_atoms for match in calculator.matched] == [(1, 2), (2, 3)]


def test_matched_parameters_record_the_row_number_in_each_source_file(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(2, 3, "C", "O", 1.4)],
        [bond_row(2, 3, "C", "O", 1.3), bond_row(1, 2, "C", "C", 1.5)],
    )

    calculator.run()

    by_atoms = {match.local_atoms: match for match in calculator.matched}
    assert by_atoms[(1, 2)].source_row_1 == 1
    assert by_atoms[(1, 2)].source_row_2 == 2
    assert by_atoms[(2, 3)].source_row_1 == 2
    assert by_atoms[(2, 3)].source_row_2 == 1


def test_dihedral_differences_are_wrapped_not_subtracted(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 179.0)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", -179.0)],
        writer=write_dihedral_file,
    )

    calculator.run()

    assert calculator.matched[0].difference == pytest.approx(2.0)


def test_angle_differences_are_plain_subtraction(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [angle_row(1, 2, 3, "H", "C", "H", 109.5)],
        [angle_row(1, 2, 3, "H", "C", "H", 111.0)],
        writer=write_angle_file,
    )

    calculator.run()

    assert calculator.matched[0].difference == pytest.approx(1.5)


def test_largest_shifts_ranks_by_absolute_difference(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.50), bond_row(2, 3, "C", "O", 1.40),
         bond_row(3, 4, "O", "H", 1.00)],
        [bond_row(1, 2, "C", "C", 1.51), bond_row(2, 3, "C", "O", 1.10),
         bond_row(3, 4, "O", "H", 1.05)],
    )

    calculator.run()

    assert [match.local_atoms for match in calculator.largest_shifts(2)] == [(2, 3), (3, 4)]


def test_matched_parameters_carry_both_percentages(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
    )

    calculator.run()

    match = calculator.matched[0]
    assert match.percent_difference == pytest.approx(100 * 0.1 / 1.5)
    assert match.percent_error == pytest.approx(100 * 0.1 / 1.45)


def test_percentages_are_computed_even_when_no_mode_is_selected(tmp_path):
    """The switches gate the output columns, not the arithmetic."""
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.45)],
        [bond_row(1, 2, "C", "C", 1.55)],
        percent_modes=(),
    )

    calculator.run()

    assert calculator.matched[0].percent_difference == pytest.approx(100 * 0.1 / 1.5)


def test_a_near_zero_dihedral_gets_an_undefined_percentage(tmp_path):
    import math
    calculator = make_calculator(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.40)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.80)],
        writer=write_dihedral_file,
    )

    calculator.run()

    match = calculator.matched[0]
    assert math.isnan(match.percent_difference)
    assert match.difference == pytest.approx(0.40)  # the absolute shift survives


def test_the_calculator_counts_undefined_percentages(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.40),
         dihedral_row(2, 3, 4, 5, "C", "C", "H", "H", 120.0)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 0.80),
         dihedral_row(2, 3, 4, 5, "C", "C", "H", "H", 122.0)],
        writer=write_dihedral_file,
        percent_modes=PERCENT_MODES,
    )

    calculator.run()

    assert calculator.percent_undefined == 1


# --------------------------------------------------------------------------- #
# Comparison: the occurrence filter                                             #
# --------------------------------------------------------------------------- #
def test_occurrence_filter_drops_rows_below_the_threshold(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(3, 4, "H", "H", 1.69, occurrence=0.0004)],
        [bond_row(1, 2, "C", "C", 1.5), bond_row(3, 4, "H", "H", 1.68, occurrence=0.0004)],
        min_occurrence=0.01,
    )

    calculator.run()

    assert [match.local_atoms for match in calculator.matched] == [(1, 2)]
    assert calculator.filtered_out_1 == 1
    assert calculator.filtered_out_2 == 1


def test_default_threshold_keeps_every_parameter(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(3, 4, "H", "H", 1.69, occurrence=0.0004)],
        [bond_row(1, 2, "C", "C", 1.5), bond_row(3, 4, "H", "H", 1.68, occurrence=0.0004)],
    )

    calculator.run()

    assert len(calculator.matched) == 2
    assert calculator.filtered_out_1 == 0


def test_filtering_happens_before_matching_so_no_phantom_unmatched_row(tmp_path):
    """A pair that survives the cutoff in one file only must not look unmatched."""
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(3, 4, "H", "H", 1.69, occurrence=0.5)],
        [bond_row(1, 2, "C", "C", 1.5), bond_row(3, 4, "H", "H", 1.68, occurrence=0.0004)],
        min_occurrence=0.01,
    )

    calculator.run()

    assert len(calculator.matched) == 1
    assert [row.local_atoms for row in calculator.only_in_file_1] == [(3, 4)]
    assert calculator.only_in_file_2 == []


# --------------------------------------------------------------------------- #
# Comparison: aborts                                                            #
# --------------------------------------------------------------------------- #
def test_comparing_two_different_parameter_types_aborts(tmp_path):
    file_1 = write_bond_file(tmp_path / "b.txt", [bond_row(1, 2, "C", "C", 1.5)])
    file_2 = write_angle_file(tmp_path / "a.txt", [angle_row(1, 2, 3, "H", "C", "H", 109.5)])

    with pytest.raises(ValueError) as excinfo:
        MolGeomCalculator(file_1, file_2).run()

    message = str(excinfo.value)
    assert "bond distance" in message and "bond angle" in message


def test_element_mismatch_aborts_naming_the_atom_and_both_elements(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(2, 3, "C", "O", 1.4)],
        [bond_row(1, 2, "C", "C", 1.5), bond_row(2, 3, "C", "N", 1.4)],
    )

    with pytest.raises(ValueError) as excinfo:
        calculator.run()

    message = str(excinfo.value)
    assert "atom 3" in message
    assert "'O'" in message and "'N'" in message
    assert "file_1.txt" in message and "file_2.txt" in message


def test_different_molecule_sizes_abort_with_both_counts(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5)],
        [bond_row(1, 2, "C", "C", 1.5), bond_row(2, 3, "C", "O", 1.4)],
    )

    with pytest.raises(ValueError) as excinfo:
        calculator.run()

    message = str(excinfo.value)
    assert "2" in message and "3" in message


def test_no_matched_parameters_at_all_aborts(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(3, 4, "C", "C", 1.5)],
        [bond_row(1, 3, "C", "C", 1.5), bond_row(2, 4, "C", "C", 1.5)],
    )

    with pytest.raises(ValueError) as excinfo:
        calculator.run()

    assert "no " in str(excinfo.value).lower()


def test_a_negative_occurrence_threshold_is_rejected(tmp_path):
    calculator = make_calculator(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5)],
        [bond_row(1, 2, "C", "C", 1.5)],
        min_occurrence=-0.1,
    )

    with pytest.raises(ValueError) as excinfo:
        calculator.run()

    assert "occurrence" in str(excinfo.value).lower()


# --------------------------------------------------------------------------- #
# Output file                                                                   #
# --------------------------------------------------------------------------- #
def header_values(text):
    """Parse the ``# key value`` metadata block into a dict."""
    values = {}
    for line in text.splitlines():
        if not line.startswith("#"):
            continue
        content = line.lstrip("#").strip()
        if not content:
            continue
        key, _, value = content.partition(" ")
        if key in ("row", "section"):
            continue
        values[key] = value.strip()
    return values


def section_of(text, name):
    """The lines of one ``# section <name>`` block: its '# row' header and its rows."""
    lines = text.splitlines()
    start = lines.index(f"# section {name}")
    block = []
    for line in lines[start + 1:]:
        if line.startswith("# section "):
            break
        if line.strip() and line.strip() != "#":
            block.append(line)
    return block


def section_columns(text, name):
    return section_of(text, name)[0].lstrip("#").strip().split()


def section_rows(text, name):
    return [line.split() for line in section_of(text, name) if not line.startswith("#")]


def run_and_write(tmp_path, rows_1, rows_2, *, out="out.txt", **kwargs):
    calculator = make_calculator(tmp_path, rows_1, rows_2, **kwargs)
    calculator.run()
    path = calculator.write_results(str(tmp_path / out))
    return calculator, path, open(path).read()


def test_write_results_creates_the_file_and_returns_its_absolute_path(tmp_path):
    _, path, _ = run_and_write(
        tmp_path, [bond_row(1, 2, "C", "C", 1.5)], [bond_row(1, 2, "C", "C", 1.6)]
    )

    assert os.path.isabs(path)
    assert os.path.isfile(path)


def test_output_header_records_both_sources_and_the_parameter_type(tmp_path):
    _, _, text = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5)],
        [bond_row(1, 2, "C", "C", 1.6)],
        kwargs_1={"frames": 1},
        kwargs_2={"frames": 15121},
    )

    header = header_values(text)
    assert text.splitlines()[0] == "# gQTEA Molecular Geometry Comparison"
    assert header["parameter_type"] == "bond"
    assert header["value_unit"] == "angstrom"
    assert header["file_1"].endswith("file_1.txt")
    assert header["file_2"].endswith("file_2.txt")
    assert header["file_1_frames"] == "1"
    assert header["file_2_frames"] == "15121"
    assert header["file_1_atom_scope"] == "all_atoms"


def test_output_header_counts_matched_and_unmatched_parameters(tmp_path):
    _, _, text = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(1, 3, "C", "O", 1.4)],
        [bond_row(1, 2, "C", "C", 1.6), bond_row(2, 3, "C", "O", 1.3)],
    )

    header = header_values(text)
    assert header["matched_parameters"] == "1"
    assert header["only_in_file_1"] == "1"
    assert header["only_in_file_2"] == "1"


def test_output_header_records_the_occurrence_threshold_and_what_it_dropped(tmp_path):
    _, _, text = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(3, 4, "H", "H", 1.69, occurrence=0.0004)],
        [bond_row(1, 2, "C", "C", 1.6), bond_row(3, 4, "H", "H", 1.68, occurrence=1.0)],
        min_occurrence=0.01,
    )

    header = header_values(text)
    assert float(header["minimum_occurrence_fraction"]) == pytest.approx(0.01)
    assert header["file_1_rows_below_threshold"] == "1"
    assert header["file_2_rows_below_threshold"] == "0"


def test_output_header_records_the_index_remap_of_a_solvated_file(tmp_path):
    """A six-atom stand-in for the real case, where the oxygens sit at 295-297."""
    _, _, text = run_and_write(
        tmp_path,
        [
            bond_row(1, 2, "C", "C", 1.50), bond_row(2, 3, "C", "C", 1.51),
            bond_row(3, 4, "C", "O", 1.40), bond_row(4, 5, "O", "O", 1.45),
            bond_row(5, 6, "O", "O", 1.46),
        ],
        [
            bond_row(1, 2, "C", "C", 1.52), bond_row(2, 3, "C", "C", 1.53),
            bond_row(3, 295, "C", "O", 1.42), bond_row(295, 296, "O", "O", 1.47),
            bond_row(296, 297, "O", "O", 1.48),
        ],
        kwargs_2={"solute": [1, 2, 3, 295, 296, 297]},
    )

    header = header_values(text)
    assert "4:295" in header["file_2_local_to_global"]
    assert "6:297" in header["file_2_local_to_global"]


def test_no_remap_line_is_written_when_the_numbering_is_unchanged(tmp_path):
    _, _, text = run_and_write(
        tmp_path, [bond_row(1, 2, "C", "C", 1.5)], [bond_row(1, 2, "C", "C", 1.6)]
    )

    assert "local_to_global" not in text


def test_matched_columns_are_named_after_the_two_labels(tmp_path):
    _, _, text = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5)],
        [bond_row(1, 2, "C", "C", 1.6)],
        label_1="isolated",
        label_2="in water box",
    )

    columns = section_columns(text, "matched")
    assert columns == [
        "row", "atom_i", "atom_j", "element_i", "element_j",
        "average_isolated", "std_dev_isolated", "occurrence_isolated",
        "average_in_water_box", "std_dev_in_water_box", "occurrence_in_water_box",
        "difference_in_water_box_minus_isolated",
        "source_row_1", "source_row_2",
    ]


def test_identical_labels_are_disambiguated_in_the_column_names(tmp_path):
    _, _, text = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5)],
        [bond_row(1, 2, "C", "C", 1.6)],
        label_1="run",
        label_2="run",
    )

    columns = section_columns(text, "matched")
    assert "average_run_1" in columns
    assert "average_run_2" in columns


def test_matched_rows_hold_the_values_and_the_traceback_columns(tmp_path):
    _, _, text = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.54834003, std_dev=0.0)],
        [bond_row(1, 2, "C", "C", 1.52968232, std_dev=0.03377609)],
    )

    columns = section_columns(text, "matched")
    row = dict(zip(columns, section_rows(text, "matched")[0]))
    assert row["atom_i"] == "1" and row["atom_j"] == "2"
    assert row["element_i"] == "C" and row["element_j"] == "C"
    assert float(row["average_isolated"]) == pytest.approx(1.54834003)
    assert float(row["std_dev_solvated"]) == pytest.approx(0.03377609)
    assert float(row["difference_solvated_minus_isolated"]) == pytest.approx(-0.01865771)
    assert row["source_row_1"] == "1" and row["source_row_2"] == "1"


def test_unmatched_sections_list_the_rows_of_each_file(tmp_path):
    _, _, text = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(1, 3, "C", "O", 1.4)],
        [bond_row(1, 2, "C", "C", 1.6), bond_row(2, 3, "C", "O", 1.3)],
    )

    only_1 = section_rows(text, "only_in_file_1")
    only_2 = section_rows(text, "only_in_file_2")
    assert len(only_1) == 1 and len(only_2) == 1
    columns = section_columns(text, "only_in_file_1")
    row = dict(zip(columns, only_1[0]))
    assert (row["atom_i"], row["atom_j"]) == ("1", "3")
    assert float(row["average"]) == pytest.approx(1.4)
    assert row["source_row"] == "2"


def test_empty_unmatched_sections_still_appear_with_their_header(tmp_path):
    _, _, text = run_and_write(
        tmp_path, [bond_row(1, 2, "C", "C", 1.5)], [bond_row(1, 2, "C", "C", 1.6)]
    )

    assert "# section only_in_file_1" in text
    assert section_rows(text, "only_in_file_1") == []


def test_angle_output_carries_three_atom_columns(tmp_path):
    _, _, text = run_and_write(
        tmp_path,
        [angle_row(1, 2, 3, "H", "C", "H", 109.5)],
        [angle_row(1, 2, 3, "H", "C", "H", 111.0)],
        writer=write_angle_file,
    )

    columns = section_columns(text, "matched")
    assert columns[1:7] == ["atom_i", "atom_j", "atom_k", "element_i", "element_j", "element_k"]
    assert header_values(text)["value_unit"] == "degrees"


def test_dihedral_output_carries_four_atom_columns(tmp_path):
    _, _, text = run_and_write(
        tmp_path,
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 59.0)],
        [dihedral_row(1, 2, 3, 4, "H", "C", "C", "H", 61.0)],
        writer=write_dihedral_file,
    )

    columns = section_columns(text, "matched")
    assert columns[1:9] == [
        "atom_i", "atom_j", "atom_k", "atom_l",
        "element_i", "element_j", "element_k", "element_l",
    ]
    assert header_values(text)["parameter_type"] == "dihedral"


def test_writing_before_running_the_comparison_is_refused(tmp_path):
    calculator = make_calculator(
        tmp_path, [bond_row(1, 2, "C", "C", 1.5)], [bond_row(1, 2, "C", "C", 1.6)]
    )

    with pytest.raises(ValueError) as excinfo:
        calculator.write_results(str(tmp_path / "out.txt"))

    assert "run" in str(excinfo.value).lower()


def test_output_file_can_be_reread_column_by_column(tmp_path):
    """The matched section must survive a name-based numeric read."""
    np = pytest.importorskip("numpy")
    _, path, text = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.50), bond_row(2, 3, "C", "O", 1.40)],
        [bond_row(1, 2, "C", "C", 1.55), bond_row(2, 3, "C", "O", 1.30)],
    )

    columns = section_columns(text, "matched")
    rows = section_rows(text, "matched")
    differences = np.array(
        [float(row[columns.index("difference_solvated_minus_isolated")]) for row in rows]
    )
    np.testing.assert_allclose(differences, [0.05, -0.10], atol=1e-8)


def test_summary_lines_report_the_counts_and_the_output_path(tmp_path):
    calculator, path, _ = run_and_write(
        tmp_path,
        [bond_row(1, 2, "C", "C", 1.5), bond_row(1, 3, "C", "O", 1.4)],
        [bond_row(1, 2, "C", "C", 1.6), bond_row(2, 3, "C", "O", 1.3)],
    )

    summary = calculator.summary_text(path)

    assert "bond distance" in summary
    assert "Matched parameters: 1" in summary
    assert path in summary
    assert "C1-C2" in summary


# --------------------------------------------------------------------------- #
# UI input helpers                                                              #
# --------------------------------------------------------------------------- #
def test_blank_occurrence_field_means_keep_everything():
    assert parse_optional_fraction("") == 0.0
    assert parse_optional_fraction("   ") == 0.0


def test_occurrence_field_parses_a_fraction():
    assert parse_optional_fraction("0.05") == pytest.approx(0.05)


def test_occurrence_field_rejects_text_that_is_not_a_number():
    with pytest.raises(ValueError) as excinfo:
        parse_optional_fraction("a bit")

    assert "number" in str(excinfo.value).lower()


def test_occurrence_field_rejects_a_value_outside_zero_to_one():
    with pytest.raises(ValueError):
        parse_optional_fraction("-0.1")
    with pytest.raises(ValueError):
        parse_optional_fraction("1.5")


def test_blank_output_folder_falls_back_to_the_first_file_folder(tmp_path):
    reference = tmp_path / "sub" / "file_1.txt"
    reference.parent.mkdir()
    reference.write_text("x")

    assert resolve_output_dir(str(reference), "") == str(tmp_path / "sub")


def test_a_chosen_output_folder_wins(tmp_path):
    reference = tmp_path / "file_1.txt"
    reference.write_text("x")
    chosen = tmp_path / "elsewhere"
    chosen.mkdir()

    assert resolve_output_dir(str(reference), str(chosen)) == str(chosen)


# --------------------------------------------------------------------------- #
# UI plumbing that can be checked without a display                             #
# --------------------------------------------------------------------------- #
class FakeWidget:
    def __init__(self, value=""):
        self.value = value
        self.placeholder = ""


class FakeWindow:
    """Records which Toga dialog type was raised.

    Toga dialogs keep their title/message inside the backend impl and expose
    neither, so the type name is what can be asserted portably.
    """

    def __init__(self):
        self.dialogs = []

    async def dialog(self, dlg):
        self.dialogs.append(type(dlg).__name__)
        return None


class TestUIPlumbing:
    @staticmethod
    def bare_ui():
        """A MolGeomComparatorUI with stub widgets and no Toga window."""
        ui = MolGeomComparatorUI.__new__(MolGeomComparatorUI)
        ui.file_1 = None
        ui.file_2 = None
        ui.parsed_1 = None
        ui.parsed_2 = None
        ui.output_dir = os.getcwd()
        ui.textInput_file_1 = FakeWidget()
        ui.textInput_file_2 = FakeWidget()
        ui.textInput_label_1 = FakeWidget("isolated")
        ui.textInput_label_2 = FakeWidget("solvated")
        ui.textInput_min_occurrence = FakeWidget()
        ui.textInput_output_dir = FakeWidget()
        ui.textInput_output = FakeWidget()
        ui.multi_line_text = FakeWidget()
        ui.main_window = FakeWindow()
        return ui

    def test_accepts_the_button_toga_passes(self):
        """gqteaWinToga registers the class itself as on_press, so __init__ is
        handed the Button. It must take *args like every other tool."""
        import inspect
        params = list(inspect.signature(MolGeomComparatorUI.__init__).parameters.values())
        assert any(p.kind is inspect.Parameter.VAR_POSITIONAL for p in params), (
            "MolGeomComparatorUI.__init__ must accept *args"
        )

    def test_read_params_refuses_to_run_without_both_files(self):
        ui = self.bare_ui()

        assert asyncio.run(ui.read_params()) is False
        assert ui.main_window.dialogs == ["InfoDialog"]

    def test_read_params_accepts_a_path_typed_into_the_field(self, tmp_path):
        """The file fields are editable, so a pasted path must work like Browse."""
        ui = self.bare_ui()
        ui.textInput_file_1.value = write_bond_file(
            tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.5)]
        )
        ui.textInput_file_2.value = write_bond_file(
            tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.6)]
        )

        assert asyncio.run(ui.read_params()) is True
        assert ui.file_1 == str(tmp_path / "file_1.txt")
        assert ui.file_2 == str(tmp_path / "file_2.txt")

    def test_read_params_rejects_a_typed_path_that_does_not_exist(self, tmp_path):
        ui = self.bare_ui()
        ui.textInput_file_1.value = str(tmp_path / "missing.txt")
        ui.textInput_file_2.value = write_bond_file(
            tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.6)]
        )

        assert asyncio.run(ui.read_params()) is False
        assert ui.main_window.dialogs == ["InfoDialog"]

    def test_read_params_refuses_the_same_file_twice(self, tmp_path):
        ui = self.bare_ui()
        path = write_bond_file(tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.5)])
        ui.file_1 = path
        ui.file_2 = path

        assert asyncio.run(ui.read_params()) is False
        assert ui.main_window.dialogs == ["InfoDialog"]

    def test_read_params_refuses_a_bad_occurrence_fraction(self, tmp_path):
        ui = self.bare_ui()
        ui.file_1 = write_bond_file(tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.5)])
        ui.file_2 = write_bond_file(tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.6)])
        ui.textInput_min_occurrence.value = "half"

        assert asyncio.run(ui.read_params()) is False
        assert ui.main_window.dialogs == ["InfoDialog"]

    def test_read_params_resolves_a_relative_output_name_against_the_folder(self, tmp_path):
        ui = self.bare_ui()
        ui.file_1 = write_bond_file(tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.5)])
        ui.file_2 = write_bond_file(tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.6)])
        ui.output_dir = str(tmp_path)
        ui.textInput_output.value = "result.txt"

        assert asyncio.run(ui.read_params()) is True
        assert ui.output_file == str(tmp_path / "result.txt")

    def test_read_params_keeps_an_absolute_output_name(self, tmp_path):
        ui = self.bare_ui()
        ui.file_1 = write_bond_file(tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.5)])
        ui.file_2 = write_bond_file(tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.6)])
        absolute = str(tmp_path / "explicit.txt")
        ui.textInput_output.value = absolute

        assert asyncio.run(ui.read_params()) is True
        assert ui.output_file == absolute

    def test_read_params_names_the_output_after_the_parameter_type_when_blank(self, tmp_path):
        ui = self.bare_ui()
        ui.file_1 = write_angle_file(tmp_path / "file_1.txt", [angle_row(1, 2, 3, "H", "C", "H", 109.5)])
        ui.file_2 = write_angle_file(tmp_path / "file_2.txt", [angle_row(1, 2, 3, "H", "C", "H", 111.0)])
        ui.parsed_1 = read_parameter_file(ui.file_1)
        ui.output_dir = str(tmp_path)

        assert asyncio.run(ui.read_params()) is True
        assert ui.output_file == str(tmp_path / "geometry_comparison_angle.txt")

    def test_workflow_reports_a_bad_file_as_a_dialog_not_a_traceback(self, tmp_path):
        """A parse failure must never escape the async handler."""
        ui = self.bare_ui()
        broken = tmp_path / "broken.txt"
        broken.write_text("not a gqtea file\n")
        ui.file_1 = str(broken)
        ui.file_2 = write_bond_file(tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.6)])
        ui.output_dir = str(tmp_path)

        asyncio.run(ui.workflow(None))

        assert ui.main_window.dialogs == ["InfoDialog"]

    def test_workflow_writes_the_file_and_reports_the_summary(self, tmp_path):
        ui = self.bare_ui()
        ui.file_1 = write_bond_file(tmp_path / "file_1.txt", [bond_row(1, 2, "C", "C", 1.50)])
        ui.file_2 = write_bond_file(tmp_path / "file_2.txt", [bond_row(1, 2, "C", "C", 1.55)])
        ui.output_dir = str(tmp_path)

        asyncio.run(ui.workflow(None))

        assert ui.main_window.dialogs == []
        assert os.path.isfile(tmp_path / "geometry_comparison_bond.txt")
        assert "Matched parameters: 1" in ui.multi_line_text.value
