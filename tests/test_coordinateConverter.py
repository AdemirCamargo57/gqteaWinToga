"""Angstrom <-> bohr support in coordinateConverter.py.

Drives the headless ``CoordinateConverter`` logic class directly; no Toga
window is created. ``data/2criegee.xyz`` is a 15-atom, 8-decimal Angstrom XYZ
(three stacked Criegee-intermediate geometries) with trailing blank lines.
"""

import shutil
from pathlib import Path

import pytest

from coordinateConverter import (
    BOHR_RADIUS_ANGSTROM,
    UNIT_ANGSTROM,
    UNIT_BOHR,
    CoordinateConverter,
    convert_atom_units,
    convert_length,
)

DATA = Path(__file__).parent / "data" / "2criegee.xyz"


def _read_xyz_rows(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    count = int(lines[0])
    rows = []
    for line in lines[2:2 + count]:
        symbol, x, y, z = line.split()[:4]
        rows.append((symbol, x, y, z))
    return lines[1], rows


def _converter(input_file, output_file, in_fmt="XYZ", out_fmt="XYZ",
               in_unit=UNIT_ANGSTROM, out_unit=UNIT_BOHR):
    conv = CoordinateConverter()
    conv.set_input_file(str(input_file))
    conv.input_format, conv.output_format = in_fmt, out_fmt
    conv.input_unit, conv.output_unit = in_unit, out_unit
    conv.output_file = str(output_file)
    return conv


@pytest.fixture
def criegee(tmp_path):
    target = tmp_path / "2criegee.xyz"
    shutil.copy(DATA, target)
    return target


# --------------------------------------------------------------------------
# Conversion factor and arithmetic
# --------------------------------------------------------------------------

def test_bohr_radius_is_codata_2022():
    assert BOHR_RADIUS_ANGSTROM == 0.529177210544


def test_convert_length_directions():
    assert convert_length(1.0, UNIT_BOHR, UNIT_ANGSTROM) == BOHR_RADIUS_ANGSTROM
    assert convert_length(BOHR_RADIUS_ANGSTROM, UNIT_ANGSTROM, UNIT_BOHR) == 1.0
    assert convert_length(1.0, UNIT_ANGSTROM, UNIT_BOHR) == pytest.approx(1.8897261259077822, rel=1e-15)


def test_same_unit_is_exact_identity():
    value = 3.68672253
    assert convert_length(value, UNIT_ANGSTROM, UNIT_ANGSTROM) == value
    assert convert_length(value, UNIT_BOHR, UNIT_BOHR) == value


def test_round_trip_within_one_ulp():
    rng_values = [0.0, -2.23096882, 1e-8, 123.456789012, -9876.54321]
    for value in rng_values:
        back = convert_length(convert_length(value, UNIT_ANGSTROM, UNIT_BOHR), UNIT_BOHR, UNIT_ANGSTROM)
        assert back == pytest.approx(value, rel=2.3e-16, abs=0.0)


def test_unknown_unit_rejected():
    with pytest.raises(ValueError, match="Unsupported coordinate unit"):
        convert_length(1.0, "nm", UNIT_BOHR)


def test_convert_atom_units_does_not_mutate_input():
    atoms = [{"symbol": "C", "x": 1.0, "y": 2.0, "z": 3.0}]
    converted = convert_atom_units(atoms, UNIT_ANGSTROM, UNIT_BOHR)
    assert atoms[0]["x"] == 1.0
    assert converted[0]["symbol"] == "C"
    assert converted[0]["z"] == 3.0 / BOHR_RADIUS_ANGSTROM


# --------------------------------------------------------------------------
# File round trip on 2criegee.xyz
# --------------------------------------------------------------------------

def test_xyz_angstrom_to_bohr_on_criegee(criegee, tmp_path):
    out = tmp_path / "2criegee_bohr.xyz"
    assert _converter(criegee, out).convert_file() == 15

    _, original = _read_xyz_rows(criegee)
    comment, converted = _read_xyz_rows(out)
    assert comment == "Converted file (coordinates in bohr)"
    assert [row[0] for row in converted] == [row[0] for row in original]
    for (_, *xyz_a), (_, *xyz_b) in zip(original, converted):
        for a, b in zip(xyz_a, xyz_b):
            assert float(b) == pytest.approx(float(a) / BOHR_RADIUS_ANGSTROM, abs=5e-11)
    # First carbon: 1.99447596 Angstrom = 3.7690133291 bohr
    assert converted[0][1] == "3.7690133291"


def test_criegee_round_trip_recovers_every_digit(criegee, tmp_path):
    bohr = tmp_path / "bohr.xyz"
    back = tmp_path / "back.xyz"
    _converter(criegee, bohr).convert_file()
    _converter(bohr, back, in_unit=UNIT_BOHR, out_unit=UNIT_ANGSTROM).convert_file()

    _, original = _read_xyz_rows(criegee)
    _, recovered = _read_xyz_rows(back)
    for (sym_a, *xyz_a), (sym_b, *xyz_b) in zip(original, recovered):
        assert sym_a == sym_b
        for a, b in zip(xyz_a, xyz_b):
            # 10 written decimals: rounding back to the input's 8 is exact.
            assert f"{float(b):.8f}" == a


def test_bohr_to_angstrom_round_trip_from_bohr_input(tmp_path):
    src = tmp_path / "in_bohr.xyz"
    src.write_text("2\nbohr input\nO 1.23456789 -2.00000001 0.00000003\nH 0.5 0.25 -0.125\n",
                   encoding="utf-8")
    ang = tmp_path / "ang.xyz"
    back = tmp_path / "back.xyz"
    _converter(src, ang, in_unit=UNIT_BOHR, out_unit=UNIT_ANGSTROM).convert_file()
    _converter(ang, back).convert_file()
    _, rows = _read_xyz_rows(back)
    assert [f"{float(v):.8f}" for v in rows[0][1:]] == ["1.23456789", "-2.00000001", "0.00000003"]


def test_angstrom_xyz_to_pdb_unchanged_numbers(criegee, tmp_path):
    out = tmp_path / "c.pdb"
    _converter(criegee, out, out_fmt="PDB", out_unit=UNIT_ANGSTROM).convert_file()
    first_atom = next(l for l in out.read_text().splitlines() if l.startswith("ATOM"))
    assert float(first_atom[30:38]) == pytest.approx(1.994, abs=1e-3)


def test_bohr_xyz_to_pdb_lands_in_angstrom(criegee, tmp_path):
    bohr = tmp_path / "bohr.xyz"
    pdb = tmp_path / "c.pdb"
    _converter(criegee, bohr).convert_file()
    _converter(bohr, pdb, out_fmt="PDB", in_unit=UNIT_BOHR, out_unit=UNIT_ANGSTROM).convert_file()
    first_atom = next(l for l in pdb.read_text().splitlines() if l.startswith("ATOM"))
    assert float(first_atom[30:38]) == pytest.approx(1.994, abs=1e-3)


# --------------------------------------------------------------------------
# Gaussian declares its unit in the route section
# --------------------------------------------------------------------------

def test_gaussian_bohr_output_declares_units(criegee, tmp_path):
    out = tmp_path / "c.gjf"
    _converter(criegee, out, out_fmt="Gaussian").convert_file()
    text = out.read_text(encoding="utf-8")
    assert "#p hf/3-21g Units=Bohr" in text

    reread = CoordinateConverter()
    reread.set_input_file(str(out))
    assert reread.input_format == "Gaussian"
    assert reread.input_unit == UNIT_BOHR
    atoms = reread.read_atoms("Gaussian")
    assert atoms[0]["x"] == pytest.approx(1.99447596 / BOHR_RADIUS_ANGSTROM, abs=5e-11)


def test_gaussian_angstrom_output_has_no_units_keyword(criegee, tmp_path):
    out = tmp_path / "c.gjf"
    _converter(criegee, out, out_fmt="Gaussian", out_unit=UNIT_ANGSTROM).convert_file()
    assert "Units" not in out.read_text(encoding="utf-8")
    assert CoordinateConverter.guess_gaussian_unit(str(out)) == UNIT_ANGSTROM


@pytest.mark.parametrize("route", ["#p b3lyp Units=AU", "#P units(au,deg) opt", "#p Units=(Bohr,Rad)"])
def test_gaussian_unit_guess_variants(tmp_path, route):
    path = tmp_path / "g.com"
    path.write_text(f"{route}\n\ntitle\n\n0 1\nH 0 0 0\n\n", encoding="utf-8")
    assert CoordinateConverter.guess_gaussian_unit(str(path)) == UNIT_BOHR


# --------------------------------------------------------------------------
# Validation and output naming
# --------------------------------------------------------------------------

def test_same_format_same_unit_rejected(criegee, tmp_path):
    conv = _converter(criegee, tmp_path / "o.xyz", out_unit=UNIT_ANGSTROM)
    with pytest.raises(ValueError, match="formats and units are the same"):
        conv.convert_file()


@pytest.mark.parametrize("fmt", ["PDB", "ENT", "MOL", "MOL2", "HIN"])
def test_bohr_rejected_for_angstrom_only_output(criegee, tmp_path, fmt):
    conv = _converter(criegee, tmp_path / "o.out", out_fmt=fmt, out_unit=UNIT_BOHR)
    with pytest.raises(ValueError, match="Angstrom by specification"):
        conv.convert_file()


def test_bohr_rejected_for_angstrom_only_input(criegee, tmp_path):
    conv = _converter(criegee, tmp_path / "o.xyz", in_fmt="PDB", in_unit=UNIT_BOHR, out_unit=UNIT_ANGSTROM)
    with pytest.raises(ValueError, match="Input format PDB"):
        conv.validate_conversion()


def test_refuses_to_overwrite_input(criegee):
    conv = _converter(criegee, criegee)
    with pytest.raises(ValueError, match="overwrite the input"):
        conv.convert_file()


def test_set_input_file_forces_angstrom_for_pdb(tmp_path):
    conv = CoordinateConverter()
    conv.input_unit = UNIT_BOHR
    conv.set_input_file(str(tmp_path / "x.pdb"))
    assert conv.input_unit == UNIT_ANGSTROM


class _StubWidget:
    def __init__(self, value=None):
        self.value = value
        self.enabled = True
        self.text = ""


def _stub_ui(input_format="XYZ", output_format="XYZ"):
    from coordinateConverter import CoordinateConverterUI

    ui = CoordinateConverterUI.__new__(CoordinateConverterUI)
    CoordinateConverter.__init__(ui)
    ui.input_format, ui.output_format = input_format, output_format
    ui.input_unit_selection = _StubWidget()
    ui.output_unit_selection = _StubWidget()
    ui.output_format_selection = _StubWidget()
    ui.output_file_label = _StubWidget()
    return ui


def test_ui_locks_angstrom_only_formats():
    ui = _stub_ui(output_format="PDB")
    ui.output_unit = UNIT_BOHR
    ui.sync_unit_controls()
    assert ui.output_unit == UNIT_ANGSTROM
    assert ui.output_unit_selection.value == UNIT_ANGSTROM
    assert ui.output_unit_selection.enabled is False
    assert ui.input_unit_selection.enabled is True


def test_ui_unit_change_updates_output_name(criegee):
    ui = _stub_ui()
    ui.input_file = str(criegee)
    ui.on_output_unit_change(_StubWidget(UNIT_BOHR))
    assert ui.output_unit == UNIT_BOHR
    assert ui.output_file.endswith("2criegee_converted_bohr.xyz")
    assert ui.output_file_label.text.endswith("2criegee_converted_bohr.xyz")


class _RecordingWindow:
    def __init__(self):
        self.dialogs = []

    async def dialog(self, dialog):
        self.dialogs.append(dialog)


@pytest.mark.filterwarnings("error::DeprecationWarning")
def test_conversion_dialogs_use_non_deprecated_api(criegee, tmp_path):
    import asyncio

    import toga

    ui = _stub_ui()
    ui.main_window = _RecordingWindow()
    ui.set_input_file(str(criegee))
    ui.output_unit = UNIT_BOHR
    ui.output_file = str(tmp_path / "out.xyz")
    asyncio.run(ui.perform_conversion(None))
    assert isinstance(ui.main_window.dialogs[-1], toga.InfoDialog)

    ui.output_unit = UNIT_ANGSTROM  # same format and unit -> error path
    asyncio.run(ui.perform_conversion(None))
    assert isinstance(ui.main_window.dialogs[-1], toga.ErrorDialog)


def test_ui_offers_same_format_as_output():
    ui = _stub_ui()
    assert "XYZ" in ui.get_output_formats()


def test_default_output_path_tags_unit_change(criegee):
    conv = CoordinateConverter()
    conv.set_input_file(str(criegee))
    conv.output_format = "XYZ"
    conv.output_unit = UNIT_BOHR
    assert conv.default_output_path().endswith("2criegee_converted_bohr.xyz")
    conv.output_unit = UNIT_ANGSTROM
    conv.output_format = "PDB"
    assert conv.default_output_path().endswith("2criegee_converted.pdb")
