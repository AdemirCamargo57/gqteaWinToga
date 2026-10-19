"""Tests for molecularFileFormats: structure-file import and export.

Every reader is checked on a well-formed file and on the malformed cases it
must refuse with a message that names the line; every writer is checked by
reading its own output back. Fixtures from tests/data are used where the
project already holds real files (GaussView's ethylene input with
Geom=Connectivity, Gaussian fchk files).
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import molecularFileFormats as mff
from molecularFileFormats import (
    FORMATS,
    StructureData,
    detect_bonds,
    export_notes,
    format_for_path,
    parse_structure_text,
    read_structure_file,
    structure_to_text,
    sybyl_atom_type,
)

DATA = Path(__file__).parent / "data"

WATER_XYZ = """3
water
O   0.000000   0.000000   0.117300
H   0.000000   0.757200  -0.469200
H   0.000000  -0.757200  -0.469200
"""


def _ethylene(orders=(2.0, 1.0, 1.0, 1.0, 1.0)) -> StructureData:
    coordinates = np.array([
        [0.000, 0.6676, 0.0], [0.000, -0.6676, 0.0],
        [0.9236, 1.2405, 0.0], [-0.9236, 1.2405, 0.0],
        [0.9236, -1.2405, 0.0], [-0.9236, -1.2405, 0.0],
    ])
    pairs = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]
    bonds = [(i, j, order) for (i, j), order in zip(pairs, orders)]
    return StructureData(["C", "C", "H", "H", "H", "H"], coordinates, bonds=bonds, title="ethylene")


def _benzene() -> StructureData:
    angles = np.deg2rad(np.arange(6) * 60.0)
    ring = np.c_[1.39 * np.cos(angles), 1.39 * np.sin(angles), np.zeros(6)]
    hydrogens = np.c_[2.47 * np.cos(angles), 2.47 * np.sin(angles), np.zeros(6)]
    bonds = [(k, (k + 1) % 6, 1.5) for k in range(6)] + [(k, k + 6, 1.0) for k in range(6)]
    bonds = [(min(i, j), max(i, j), o) for i, j, o in bonds]
    return StructureData(["C"] * 6 + ["H"] * 6, np.vstack([ring, hydrogens]), bonds=sorted(bonds))


# ------------------------------------------------------------------
# Module boundaries
# ------------------------------------------------------------------
def test_the_module_does_not_import_toga():
    result = subprocess.run(
        [sys.executable, "-c", "import molecularFileFormats, sys; print('toga' in sys.modules)"],
        capture_output=True, text=True,
    )
    assert result.stdout.strip() == "False", result.stderr


def test_the_bohr_radius_matches_the_coordinate_converter():
    from coordinateConverter import BOHR_RADIUS_ANGSTROM

    assert mff.BOHR_RADIUS_ANGSTROM == BOHR_RADIUS_ANGSTROM


# ------------------------------------------------------------------
# Formats and file names
# ------------------------------------------------------------------
@pytest.mark.parametrize("name, key", [
    ("a.xyz", "XYZ"), ("A.XYZ", "XYZ"), ("a.gjf", "Gaussian"), ("a.com", "Gaussian"),
    ("a.fchk", "FCHK"), ("a.fch", "FCHK"), ("a.mol", "MOL"), ("a.sdf", "SDF"),
    ("a.mol2", "MOL2"), ("a.pdb", "PDB"), ("a.ent", "PDB"),
])
def test_the_extension_selects_the_format(name, key):
    assert format_for_path(name).key == key


def test_an_unknown_extension_lists_the_supported_formats():
    with pytest.raises(ValueError) as error:
        format_for_path("notes.txt")
    message = str(error.value)
    assert "'.txt' is not a supported structure format" in message
    for fmt in FORMATS:
        assert fmt.label in message


def test_a_file_without_extension_is_refused():
    with pytest.raises(ValueError, match="no file extension"):
        format_for_path("GEOMETRY")


def test_a_binary_checkpoint_gets_a_formchk_hint():
    with pytest.raises(ValueError, match="formchk"):
        format_for_path("job.chk")


def test_fchk_can_be_imported_but_not_saved():
    assert format_for_path("job.fchk", "read").key == "FCHK"
    with pytest.raises(ValueError, match="can be imported but not saved"):
        format_for_path("job.fchk", "write")


def test_every_writable_format_has_a_writer():
    structure = _ethylene()
    for fmt in mff.WRITABLE_FORMATS:
        assert structure_to_text(structure, fmt.key)


# ------------------------------------------------------------------
# Shared checks
# ------------------------------------------------------------------
def test_bonds_are_detected_with_the_viewer_rule():
    structure = parse_structure_text(WATER_XYZ, "XYZ")
    assert structure.bonds == [(0, 1, 1.0), (0, 2, 1.0)]
    assert "no bonds" in structure.notes[0] and "2 bonds were detected" in structure.notes[0]


def test_an_element_without_a_radius_is_left_unbonded():
    assert detect_bonds(["C", "Og"], [[0, 0, 0], [1.5, 0, 0]]) == []


def test_overlapping_atoms_are_refused():
    text = "2\n\nC 0 0 0\nC 0.1 0 0\n"
    with pytest.raises(ValueError, match="atoms 1 and 2 are only 0.100 A apart"):
        parse_structure_text(text, "XYZ")


def test_a_missing_file_is_reported(tmp_path):
    with pytest.raises(ValueError, match="could not be opened"):
        read_structure_file(tmp_path / "missing.xyz")


def test_a_binary_file_is_refused(tmp_path):
    path = tmp_path / "fake.xyz"
    path.write_bytes(b"3\n\0\0\0binary")
    with pytest.raises(ValueError, match="binary file"):
        read_structure_file(path)


def test_errors_name_the_file_and_the_format(tmp_path):
    path = tmp_path / "broken.xyz"
    path.write_text("3\ncomment\nO 0 0 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"'broken.xyz' could not be imported as XYZ: line 1 announces 3"):
        read_structure_file(path)


def test_a_windows_code_page_title_is_accepted(tmp_path):
    path = tmp_path / "accent.xyz"
    path.write_bytes(WATER_XYZ.replace("water", "\xe1gua").encode("latin-1"))
    assert read_structure_file(path).title == "\xe1gua"


# ------------------------------------------------------------------
# XYZ
# ------------------------------------------------------------------
def test_xyz_atoms_and_coordinates_are_read():
    structure = parse_structure_text(WATER_XYZ, "XYZ")
    assert structure.elements == ["O", "H", "H"]
    assert structure.title == "water"
    assert np.allclose(structure.coordinates[1], [0.0, 0.7572, -0.4692])


def test_xyz_accepts_labels_atomic_numbers_and_extra_columns():
    text = "3\n\n8 0 0 0.1173 0.0 0.0 0.0\nH1 0 0.7572 -0.4692\nh2 0 -0.7572 -0.4692 extra\n"
    assert parse_structure_text(text, "XYZ").elements == ["O", "H", "H"]


def test_xyz_written_in_bohr_by_the_converter_is_converted():
    text = "2\nConverted file (coordinates in bohr)\nH 0 0 0\nH 1.4 0 0\n"
    structure = parse_structure_text(text, "XYZ")
    assert structure.coordinates[1, 0] == pytest.approx(1.4 * mff.BOHR_RADIUS_ANGSTROM)
    assert any("bohr" in note for note in structure.notes)


def test_only_the_first_frame_of_a_trajectory_is_imported():
    structure = parse_structure_text(WATER_XYZ * 3, "XYZ")
    assert structure.atom_count == 3
    assert any("3 frames" in note for note in structure.notes)


def test_dummy_atoms_are_skipped():
    text = WATER_XYZ.replace("3\n", "4\n", 1) + "X 0 0 3\n"
    structure = parse_structure_text(text, "XYZ")
    assert structure.atom_count == 3
    assert any("dummy" in note for note in structure.notes)


@pytest.mark.parametrize("text, message", [
    ("", "empty"),
    ("water\n", "line 1 should hold the number of atoms"),
    ("0\n\n", "atom count on line 1 is 0"),
    ("2\n\nO 0 0 0\nH 0 0\n", "line 4: expected 'element x y z'"),
    ("2\n\nO 0 0 0\nH 0 zero 1\n", "line 4: 'zero' is not a number"),
    ("2\n\nO 0 0 0\nQq 0 1 1\n", "line 4: 'Qq' is not an element"),
    ("2\n\nO 0 0 0\nH 0 nan 1\n", "not a finite number"),
])
def test_malformed_xyz_is_refused_with_the_line(text, message):
    with pytest.raises(ValueError, match=message):
        parse_structure_text(text, "XYZ")


# ------------------------------------------------------------------
# Gaussian input
# ------------------------------------------------------------------
def test_gaussview_connectivity_keeps_the_double_bond():
    structure = read_structure_file(DATA / "fchk" / "inputs" / "eteno_b3lyp_631g.gjf")
    assert structure.elements == ["C", "H", "H", "C", "H", "H"]
    assert structure.bonds == [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 2.0), (3, 4, 1.0), (3, 5, 1.0)]
    assert structure.notes == []


def test_a_gaussian_input_without_connectivity_gets_detected_bonds():
    structure = read_structure_file(DATA / "fchk" / "inputs" / "benzene_631gd.gjf")
    assert structure.atom_count == 12
    assert len(structure.bonds) == 12
    assert any("detected" in note for note in structure.notes)


GAUSSIAN_DECORATED = """%NProcShared=4
%Mem=2GB
#P B3LYP/6-31G(d) Opt Units=Bohr ! a comment

title line one
title line two

0,1
O(Fragment=1)   0   0.0  0.0  0.2217
H1             -1   0.0  1.4309 -0.8867
H-HO--0.4           0.0 -1.4309 -0.8867 L
C-Bq                3.0  3.0  3.0

--Link1--
%OldChk=x.chk
#P B3LYP/6-31G(d) Geom=Check Guess=Read

second job

0 1

"""


def test_gaussian_labels_freeze_codes_units_and_linked_jobs():
    structure = parse_structure_text(GAUSSIAN_DECORATED, "Gaussian")
    assert structure.elements == ["O", "H", "H"]
    assert structure.title == "title line one title line two"
    assert structure.coordinates[1, 1] == pytest.approx(1.4309 * mff.BOHR_RADIUS_ANGSTROM)
    notes = " ".join(structure.notes)
    assert "Link1" in notes and "Units=Bohr" in notes and "ghost" in notes


def test_a_gaussian_resonance_bond_order_is_kept():
    text = ("# HF/STO-3G Geom=Connectivity\n\nt\n\n0 1\nC 0 0 0\nC 1.40 0 0\nO 2.0 1.2 0\n\n"
            " 1 2 1.5\n 2 3 2.0\n 3\n\n")
    assert parse_structure_text(text, "Gaussian").bonds == [(0, 1, 1.5), (1, 2, 2.0)]


def test_an_unsupported_gaussian_bond_order_becomes_single():
    text = "# HF/STO-3G Geom=Connectivity\n\nt\n\n0 1\nC 0 0 0\nC 1.5 0 0\n\n 1 2 0.5\n 2\n\n"
    structure = parse_structure_text(text, "Gaussian")
    assert structure.bonds == [(0, 1, 1.0)]
    assert any("order other than" in note for note in structure.notes)


@pytest.mark.parametrize("text, message", [
    ("just text\n", "no route section"),
    ("# HF/STO-3G\n\nt\n\n", "ends before the charge and multiplicity"),
    ("# HF/STO-3G\n\nt\n\nzero one\nC 0 0 0\n\n", "line 5: expected the charge and spin multiplicity"),
    ("# HF/STO-3G\n\nt\n\n0 1\nO\nH 1 0.96\nH 1 0.96 2 104.5\n\n", "line 6: 'O' is not a Cartesian"),
    ("# HF/STO-3G\n\nt\n\n0 1\nO 0 0 0\nH 1 0.96 2 104.5\n\n", "line 7: .* Z-matrix"),
    ("# HF/STO-3G Geom=Check\n\nt\n\n0 1\n\n", "Geom=Check"),
    ("# HF/STO-3G Geom=AllCheck\n\n", "Geom=Check"),
])
def test_malformed_gaussian_input_is_refused(text, message):
    with pytest.raises(ValueError, match=message):
        parse_structure_text(text, "Gaussian")


def test_the_gaussian_writer_round_trips_with_connectivity():
    structure = _benzene()
    text = structure_to_text(structure, "Gaussian", name="my_benzene")
    lines = text.splitlines()
    assert lines[0] == "%Chk=my_benzene.chk"
    assert lines[1] == "#P B3LYP/6-31G(d) Opt Geom=Connectivity"
    assert "0 1" in lines
    back = parse_structure_text(text, "Gaussian")
    assert back.bonds == structure.bonds
    assert np.allclose(back.coordinates, structure.coordinates, atol=1e-8)


def test_the_gaussian_title_avoids_characters_gaussian_forbids():
    from fukuiInput import TITLE_FORBIDDEN

    title = structure_to_text(_ethylene(), "Gaussian").splitlines()[3]
    assert not set(title) & set(TITLE_FORBIDDEN)


def test_an_odd_electron_count_is_written_as_a_doublet():
    methyl = StructureData(["C", "H", "H", "H"], np.array(
        [[0, 0, 0], [1.08, 0, 0], [-0.54, 0.935, 0], [-0.54, -0.935, 0]], dtype=float))
    assert "0 2" in structure_to_text(methyl, "Gaussian").splitlines()
    assert "multiplicity 2" in export_notes(methyl, "Gaussian")[0]


# ------------------------------------------------------------------
# Gaussian formatted checkpoint
# ------------------------------------------------------------------
@pytest.mark.parametrize("name", ["eteno_b3lyp_631g.fchk", "water_631gd_6d.fchk", "benzene_631gd.fchk"])
def test_fchk_geometry_matches_the_wavefunction_reader(name):
    from gaussianWavefunction import read_fchk

    structure = read_structure_file(DATA / "fchk" / name)
    wavefunction = read_fchk(DATA / "fchk" / name)
    assert structure.elements == wavefunction.atom_symbols()
    assert np.allclose(structure.coordinates, wavefunction.coordinates * mff.BOHR_RADIUS_ANGSTROM)


def test_a_truncated_fchk_is_refused(tmp_path):
    lines = (DATA / "fchk" / "water_631gd_6d.fchk").read_text().splitlines()
    end = next(k for k, line in enumerate(lines) if line.startswith("Current cartesian coordinates"))
    path = tmp_path / "cut.fchk"
    path.write_text("\n".join(lines[:end + 2]) + "\n")
    with pytest.raises(ValueError, match="truncated"):
        read_structure_file(path)


def test_a_file_without_geometry_sections_is_not_an_fchk():
    with pytest.raises(ValueError, match="'Atomic numbers' section is missing"):
        parse_structure_text("title\nSP RB3LYP 6-31G\nNumber of atoms   I   3\n", "FCHK")


# ------------------------------------------------------------------
# MDL molfile and SD file
# ------------------------------------------------------------------
def test_the_molfile_round_trips_bond_orders_and_resonance():
    structure = _benzene()
    text = structure_to_text(structure, "MOL")
    assert text.splitlines()[3].endswith("999 V2000")
    back = parse_structure_text(text, "MOL")
    assert back.bonds == structure.bonds
    assert np.allclose(back.coordinates, structure.coordinates, atol=1e-4)


def test_a_molfile_without_a_bond_block_gets_detected_bonds():
    text = structure_to_text(StructureData(["O", "H", "H"], parse_structure_text(WATER_XYZ, "XYZ").coordinates), "MOL")
    back = parse_structure_text(text, "MOL")
    assert len(back.bonds) == 2
    assert any("detected" in note for note in back.notes)


def test_the_coordinate_converters_molfile_is_readable(tmp_path):
    from coordinateConverter import CoordinateConverter

    source = tmp_path / "water.xyz"
    source.write_text(WATER_XYZ)
    converter = CoordinateConverter()
    converter.input_file, converter.output_file = str(source), str(tmp_path / "water.mol")
    converter.write_atoms(converter.read_xyz(), "MOL")

    structure = read_structure_file(tmp_path / "water.mol")
    assert structure.elements == ["O", "H", "H"]
    assert len(structure.bonds) == 2


@pytest.mark.parametrize("text, message", [
    ("t\n\n\n", "too short"),
    ("t\n\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n", "V3000"),
    ("t\n\n\nabcdef\n", "line 4 should be the counts line"),
    ("t\n\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0\n",
     "declares 2 atoms and 0 bonds, but the file ends"),
    ("t\n\n\n  1  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0\n  1  2  1\n",
     "line 6: the bond refers to atom 2"),
])
def test_malformed_molfiles_are_refused(text, message):
    with pytest.raises(ValueError, match=message):
        parse_structure_text(text, "MOL")


def test_only_the_first_sd_record_is_imported():
    record = structure_to_text(_ethylene(), "SDF")
    assert record.endswith("$$$$\n")
    structure = parse_structure_text(record + record, "SDF")
    assert structure.atom_count == 6
    assert any("2 records" in note for note in structure.notes)


# ------------------------------------------------------------------
# Tripos MOL2
# ------------------------------------------------------------------
def test_mol2_round_trips_with_sybyl_types():
    structure = _benzene()
    text = structure_to_text(structure, "MOL2")
    assert " C.ar " in text and "     1     1     2 ar" in text
    back = parse_structure_text(text, "MOL2")
    assert back.elements == structure.elements
    assert back.bonds == structure.bonds


@pytest.mark.parametrize("element, orders, expected", [
    ("C", [1.0, 1.0, 1.0, 1.0], "C.3"), ("C", [2.0, 1.0, 1.0], "C.2"), ("C", [3.0, 1.0], "C.1"),
    ("C", [2.0, 2.0], "C.1"), ("C", [1.5, 1.5, 1.0], "C.ar"), ("N", [1.0, 1.0, 1.0], "N.3"),
    ("N", [1.0] * 4, "N.4"), ("N", [1.5, 1.5], "N.ar"), ("O", [2.0], "O.2"), ("O", [1.0, 1.0], "O.3"),
    ("Cl", [1.0], "Cl"), ("H", [1.0], "H"),
])
def test_sybyl_types_follow_the_bond_orders(element, orders, expected):
    assert sybyl_atom_type(element, orders) == expected


def test_mol2_generic_types_fall_back_on_the_atom_name():
    text = ("@<TRIPOS>MOLECULE\nm\n2 1\nSMALL\n\n@<TRIPOS>ATOM\n"
            "1 CL1 0 0 0 Hal\n2 C2 1.77 0 0 C.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 2 nc\n")
    structure = parse_structure_text(text, "MOL2")
    assert structure.elements == ["Cl", "C"]
    assert structure.bonds == [(0, 1, 1.0)]


@pytest.mark.parametrize("text, message", [
    ("@<TRIPOS>MOLECULE\nm\n", "no @<TRIPOS>ATOM section"),
    ("@<TRIPOS>ATOM\n1 C1 0 0\n", "line 2: expected 'id name x y z type"),
    ("@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 C2 1.5 0 0 C.3\n@<TRIPOS>BOND\n1 1 2 quadruple\n",
     "unknown bond type 'quadruple'"),
    ("@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n@<TRIPOS>BOND\n1 1 9 1\n", "atom 9, which is not in the ATOM section"),
])
def test_malformed_mol2_is_refused(text, message):
    with pytest.raises(ValueError, match=message):
        parse_structure_text(text, "MOL2")


# ------------------------------------------------------------------
# PDB
# ------------------------------------------------------------------
def test_pdb_round_trips_connectivity_without_orders():
    structure = _ethylene()
    text = structure_to_text(structure, "PDB")
    atom_line = text.splitlines()[1]
    assert atom_line.startswith("HETATM    1  C1  UNL A   1")
    assert atom_line[76:78] == " C"
    back = parse_structure_text(text, "PDB")
    assert [(i, j) for i, j, _ in back.bonds] == [(i, j) for i, j, _ in structure.bonds]
    assert all(order == 1.0 for _, _, order in back.bonds)
    assert np.allclose(back.coordinates, structure.coordinates, atol=1e-3)
    assert "no bond orders" in export_notes(structure, "PDB")[0]


def _pdb_atom(serial, name, x, element="", alt=" "):
    return (f"ATOM  {serial:5d} {name}{alt}ALA A   1    {x:8.3f}{0.0:8.3f}{0.0:8.3f}"
            f"  1.00  0.00          {element:>2s}")


def test_pdb_elements_come_from_the_name_when_the_column_is_blank():
    text = "\n".join([_pdb_atom(1, " CA ", 0.0), _pdb_atom(2, "CA  ", 5.0), _pdb_atom(3, " N  ", 1.47)])
    assert parse_structure_text(text, "PDB").elements == ["C", "Ca", "N"]


def test_pdb_keeps_the_first_model_and_location():
    text = "\n".join([
        "MODEL        1", _pdb_atom(1, " C  ", 0.0, "C"), _pdb_atom(2, " O  ", 1.2, "O"),
        _pdb_atom(3, " O  ", 1.3, "O", alt="B"), "ENDMDL",
        "MODEL        2", _pdb_atom(1, " C  ", 0.0, "C"), "ENDMDL",
    ])
    structure = parse_structure_text(text, "PDB")
    assert structure.elements == ["C", "O"]
    notes = " ".join(structure.notes)
    assert "several models" in notes and "alternate-location" in notes


def test_pdb_atoms_outside_conect_get_detected_bonds():
    text = "\n".join([
        _pdb_atom(1, " C  ", 0.0, "C"), _pdb_atom(2, " C  ", 1.52, "C"),
        _pdb_atom(3, " O  ", 10.0, "O"), _pdb_atom(4, " H  ", 10.96, "H"),
        "CONECT    1    2", "CONECT    2    1",
    ])
    structure = parse_structure_text(text, "PDB")
    assert structure.bonds == [(0, 1, 1.0), (2, 3, 1.0)]


def test_pdb_coordinates_that_are_not_numbers_are_refused():
    with pytest.raises(ValueError, match="line 1: the coordinates in columns 31-54"):
        parse_structure_text("ATOM      1  C   ALA A   1    abcdefgh   0.000   0.000           C", "PDB")


# ------------------------------------------------------------------
# Writers
# ------------------------------------------------------------------
@pytest.mark.parametrize("key, tolerance, keeps_orders", [
    ("XYZ", 1e-8, False), ("Gaussian", 1e-8, True), ("MOL", 1e-4, True),
    ("SDF", 1e-4, True), ("MOL2", 1e-4, True), ("PDB", 1e-3, False),
])
def test_every_writer_reads_back(key, tolerance, keeps_orders):
    structure = _ethylene()
    back = parse_structure_text(structure_to_text(structure, key), key)
    assert back.elements == structure.elements
    assert np.allclose(back.coordinates, structure.coordinates, atol=tolerance)
    assert [(i, j) for i, j, _ in back.bonds] == [(i, j) for i, j, _ in structure.bonds]
    if keeps_orders:
        assert back.bonds == structure.bonds


def test_size_limits_of_fixed_width_formats_are_reported():
    many = StructureData(["H"] * 1000, np.c_[np.arange(1000.0), np.zeros(1000), np.zeros(1000)])
    with pytest.raises(ValueError, match="at most 999 atoms"):
        structure_to_text(many, "MOL")
    far = StructureData(["H"], np.array([[20000.0, 0.0, 0.0]]))
    with pytest.raises(ValueError, match="cannot exceed"):
        structure_to_text(far, "PDB")


def test_an_unwritable_format_is_refused():
    with pytest.raises(ValueError, match="cannot be saved as FCHK"):
        structure_to_text(_ethylene(), "FCHK")
