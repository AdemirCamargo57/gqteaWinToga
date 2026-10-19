"""Structure files for the molecular design canvas: import and export.

The design canvas ([molecularDesign.py](molecularDesign.py)) imports a
structure from a file and saves the one it shows. This module is the file half
of that: it turns text into a :class:`StructureData` (elements, Cartesian
coordinates in Angstrom, bonds with orders) and back. It knows nothing about
the canvas, its pixel coordinates or Toga, so every reader and writer is
tested headless (tests/test_molecularFileFormats.py).

Supported formats (:data:`FORMATS`)::

    Format                         Read  Write  Bonds kept
    XYZ (.xyz)                     yes   yes    none -- detected from distances
    Gaussian input (.gjf .com .gau) yes  yes    orders, with Geom=Connectivity
    Gaussian fchk (.fchk .fch)     yes   no     none -- detected from distances
    MDL Molfile V2000 (.mol)       yes   yes    orders (aromatic = resonance)
    MDL SD file (.sdf .sd)         yes   yes    orders; first record only
    Tripos MOL2 (.mol2)            yes   yes    orders (ar = resonance)
    PDB (.pdb .ent)                yes   yes    connectivity (CONECT), no orders

Formats left out on purpose, because they cannot be read or written reliably
with what the project holds: Gaussian Z-matrix input (no internal-coordinate
builder; refused with a message), MOL V3000, Gaussian output logs (several
orientations per file), CIF (needs a unit cell the canvas does not have) and
HyperChem HIN (no reference file to validate against, and
coordinateConverter's HIN dialect is its own, not HyperChem's field layout).

When a file holds no connectivity, bonds are detected the way the 3D viewer
and the pre-optimizer do it -- a pair is bonded when its distance is at most
``r_i + r_j + 0.45`` Angstrom with covalent radii from help.py -- and every
detected bond is single. Every problem is a ``ValueError`` whose message names
the file, the format and, where there is one, the line.
"""
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

from fukuiInput import (
    DEFAULT_BASIS,
    DEFAULT_METHOD,
    MIN_ATOM_DISTANCE,
    element_symbol,
    format_geometry,
    keyword_name,
    route_tokens,
)
from help import AtomicData
from molecularPreOptimizer import GEOMETRY_BOND_TOLERANCE
from systemProperties import allowed_multiplicities, nuclear_charge

# coordinateConverter imports Toga, so the CODATA 2022 value is repeated here
# (as hirshfeldCharges does); a test keeps the two equal.
BOHR_RADIUS_ANGSTROM = 0.529177210544

# Bond orders as the canvas stores them (see molecularSketch).
SINGLE, RESONANCE, DOUBLE, TRIPLE = 1.0, 1.5, 2.0, 3.0
SUPPORTED_ORDERS = (SINGLE, RESONANCE, DOUBLE, TRIPLE)

# Labels that mark a dummy or ghost centre, not an atom: skipped on import.
DUMMY_LABELS = {"X", "XX", "BQ", "DU", "LP"}

# MDL bond types 1-3 are the orders; 4 is aromatic, the canvas's resonance.
MOL_BOND_ORDERS = {1: SINGLE, 2: DOUBLE, 3: TRIPLE, 4: RESONANCE}
MOL_BOND_TYPES = {SINGLE: 1, DOUBLE: 2, TRIPLE: 3, RESONANCE: 4}
# Tripos bond types; "am" (amide) is single by the SYBYL convention.
MOL2_BOND_ORDERS = {"1": SINGLE, "2": DOUBLE, "3": TRIPLE, "ar": RESONANCE,
                    "am": SINGLE, "du": SINGLE, "un": SINGLE}
MOL2_BOND_TYPES = {SINGLE: "1", DOUBLE: "2", TRIPLE: "3", RESONANCE: "ar"}

# Fixed-width field limits of the formats that have them.
MOL_MAX_COUNT = 999
MOL_MAX_COORDINATE = 9999.9999
PDB_MAX_ATOMS = 99999
PDB_MAX_COORDINATE = 9999.999


@dataclass(frozen=True)
class StructureFormat:
    """One file format: what it is called, its extensions, what it can do."""

    key: str
    label: str
    extensions: Tuple[str, ...]
    readable: bool
    writable: bool
    bonds: str  # what the format keeps of the bonds, for the documentation

    @property
    def display_name(self) -> str:
        return f"{self.label} ({', '.join(self.extensions)})"


FORMATS: Tuple[StructureFormat, ...] = (
    StructureFormat("XYZ", "XYZ", (".xyz",), True, True, "none"),
    StructureFormat("Gaussian", "Gaussian input", (".gjf", ".com", ".gau"), True, True,
                    "bond orders (Geom=Connectivity)"),
    StructureFormat("FCHK", "Gaussian formatted checkpoint", (".fchk", ".fch"), True, False, "none"),
    StructureFormat("MOL", "MDL Molfile V2000", (".mol",), True, True, "bond orders"),
    StructureFormat("SDF", "MDL SD file", (".sdf", ".sd"), True, True, "bond orders"),
    StructureFormat("MOL2", "Tripos MOL2", (".mol2",), True, True, "bond orders"),
    StructureFormat("PDB", "PDB", (".pdb", ".ent"), True, True, "connectivity, no bond orders"),
)
FORMATS_BY_KEY: Dict[str, StructureFormat] = {fmt.key: fmt for fmt in FORMATS}
READABLE_FORMATS = tuple(fmt for fmt in FORMATS if fmt.readable)
WRITABLE_FORMATS = tuple(fmt for fmt in FORMATS if fmt.writable)

# Extensions people try that are deliberately unsupported, with what to do.
UNSUPPORTED_HINTS = {
    ".chk": "Binary Gaussian checkpoints cannot be read. Convert it first: formchk name.chk name.fchk",
    ".log": "Gaussian output files are not read (they hold several orientations). "
            "Import the input file or the .fchk instead.",
    ".out": "Program output files are not read. Import the input file, or an .xyz of the geometry.",
    ".cif": "CIF crystal structures are not supported: the design canvas has no unit cell.",
    ".hin": "HyperChem .hin files are not supported.",
    ".zmat": "Z-matrix files are not supported; only Cartesian coordinates can be imported.",
}


@dataclass
class StructureData:
    """A molecule as files hold it.

    ``coordinates`` is ``(N, 3)`` in Angstrom; ``bonds`` are ``(i, j, order)``
    with 0-based atom positions and ``i < j``. ``notes`` collects what the
    reader wants the user to know (frames skipped, bonds detected, ...).
    """

    elements: List[str]
    coordinates: np.ndarray
    bonds: List[Tuple[int, int, float]] = field(default_factory=list)
    title: str = ""
    notes: List[str] = field(default_factory=list)

    @property
    def atom_count(self) -> int:
        return len(self.elements)

    def neighbour_orders(self) -> List[List[float]]:
        """Per atom, the orders of the bonds it takes part in."""
        orders: List[List[float]] = [[] for _ in self.elements]
        for i, j, order in self.bonds:
            orders[i].append(order)
            orders[j].append(order)
        return orders


# ----------------------------------------------------------------------
# Formats and file names
# ----------------------------------------------------------------------
def supported_formats_text(purpose: str = "read") -> str:
    formats = READABLE_FORMATS if purpose == "read" else WRITABLE_FORMATS
    return ", ".join(fmt.display_name for fmt in formats)


def format_for_path(path, purpose: str = "read") -> StructureFormat:
    """The format a file name's extension stands for, or a ``ValueError`` saying why not."""
    name = Path(path).name
    extension = Path(path).suffix.lower()
    action = "imported" if purpose == "read" else "saved"
    supported = f"Supported formats: {supported_formats_text(purpose)}."
    if not extension:
        raise ValueError(f"'{name}' has no file extension, so its format is unknown. {supported}")
    for fmt in FORMATS:
        if extension in fmt.extensions:
            if (fmt.readable if purpose == "read" else fmt.writable):
                return fmt
            raise ValueError(f"{fmt.label} files ({extension}) can be imported but not {action}. {supported}")
    hint = UNSUPPORTED_HINTS.get(extension, f"'{extension}' is not a supported structure format.")
    raise ValueError(f"'{name}' cannot be {action}. {hint} {supported}")


# ----------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------
def read_structure_file(path) -> StructureData:
    """Read the structure in ``path``; the format comes from the extension."""
    path = Path(path)
    fmt = format_for_path(path, "read")
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise ValueError(f"'{path.name}' could not be opened: {exc.strerror or exc}") from exc
    if b"\0" in raw:
        raise ValueError(f"'{path.name}' is a binary file, not a {fmt.label} text file.")
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")  # every byte decodes; titles may be in a Windows code page
    try:
        return parse_structure_text(text, fmt.key)
    except ValueError as exc:
        raise ValueError(f"'{path.name}' could not be imported as {fmt.label}: {exc}") from exc


def parse_structure_text(text: str, format_key: str) -> StructureData:
    """Parse ``text`` as the format ``format_key``, then check and complete it."""
    parsers = {
        "XYZ": _parse_xyz,
        "Gaussian": _parse_gaussian,
        "FCHK": _parse_fchk,
        "MOL": _parse_mol,
        "SDF": _parse_sdf,
        "MOL2": _parse_mol2,
        "PDB": _parse_pdb,
    }
    parser = parsers.get(format_key)
    if parser is None:
        raise ValueError(f"Unsupported structure format: {format_key}.")
    structure, bonds_known = parser(text)
    return _finish(structure, bonds_known, FORMATS_BY_KEY[format_key].label)


def _finish(structure: StructureData, bonds_known: bool, label: str) -> StructureData:
    """Checks every format shares, and bond detection where the file had none."""
    if not structure.elements:
        raise ValueError("the file contains no atoms.")
    coordinates = np.asarray(structure.coordinates, dtype=float).reshape(-1, 3)
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("some coordinates are not finite numbers.")
    structure.coordinates = coordinates

    if len(coordinates) > 1:
        close = sorted(cKDTree(coordinates).query_pairs(MIN_ATOM_DISTANCE))
        if close:
            i, j = close[0]
            distance = float(np.linalg.norm(coordinates[i] - coordinates[j]))
            raise ValueError(
                f"atoms {i + 1} and {j + 1} are only {distance:.3f} A apart. The file probably "
                "repeats atoms, or its coordinates are not in Angstrom."
            )

    if not bonds_known:
        structure.bonds = detect_bonds(structure.elements, coordinates)
        if len(coordinates) > 1:
            structure.notes.append(
                f"The {label} file holds no bonds: {len(structure.bonds)} bonds were detected from "
                "the interatomic distances, all as single bonds. Click a bond to set its type "
                "before pre-optimizing, or a double bond will relax to a single-bond length."
            )

    bonds, seen = [], set()
    for i, j, order in structure.bonds:
        if i == j:
            raise ValueError(f"a bond joins atom {i + 1} to itself.")
        pair = (min(i, j), max(i, j))
        if pair not in seen:  # a bond listed twice is kept once
            seen.add(pair)
            bonds.append((pair[0], pair[1], order))
    structure.bonds = sorted(bonds)
    return structure


def detect_bonds(elements: Sequence[str], coordinates,
                 tolerance: float = GEOMETRY_BOND_TOLERANCE) -> List[Tuple[int, int, float]]:
    """Single bonds wherever ``d <= r_i + r_j + tolerance`` (the 3D viewer's rule).

    Atoms whose element has no covalent radius are left unbonded.
    """
    coordinates = np.asarray(coordinates, dtype=float).reshape(-1, 3)
    radii = np.array([AtomicData.covalent_radii.get(e, np.nan) for e in elements], dtype=float)
    known = np.isfinite(radii)
    if known.sum() < 2:
        return []
    reach = 2.0 * np.nanmax(radii) + tolerance
    bonds = []
    for i, j in sorted(cKDTree(coordinates).query_pairs(reach)):
        if not (known[i] and known[j]):
            continue
        if np.linalg.norm(coordinates[i] - coordinates[j]) <= radii[i] + radii[j] + tolerance:
            bonds.append((i, j, SINGLE))
    return bonds


def _number(token: str, where: str) -> float:
    """A finite float; Fortran ``D`` exponents are accepted."""
    try:
        value = float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        raise ValueError(f"{where}: '{token}' is not a number.") from None
    if not math.isfinite(value):
        raise ValueError(f"{where}: '{token}' is not a finite number.")
    return value


def _element(token: str, where: str) -> Optional[str]:
    """Element symbol from a symbol, a numbered label or an atomic number.

    ``None`` for a dummy or ghost centre (X, Bq, Du, LP), which is skipped.
    """
    if token.upper() in DUMMY_LABELS:
        return None
    symbol = element_symbol(token)
    if symbol is None:
        raise ValueError(f"{where}: '{token}' is not an element symbol or atomic number.")
    return symbol


def _skipped_note(count: int) -> List[str]:
    if not count:
        return []
    return [f"{count} dummy or ghost centres (X, Bq, Du, LP) were skipped."]


# --- XYZ ---------------------------------------------------------------
_XYZ_UNIT_RE = re.compile(r"coordinates\s+in\s+bohr", re.IGNORECASE)


def _parse_xyz(text: str) -> Tuple[StructureData, bool]:
    lines = text.splitlines()
    while lines and not lines[0].strip():
        lines.pop(0)
    if not lines:
        raise ValueError("the file is empty.")
    first = lines[0].split()
    if not first or not re.fullmatch(r"\d+", first[0]):
        raise ValueError(
            f"line 1 should hold the number of atoms, but reads '{lines[0].strip()}'. "
            "An XYZ file starts with the atom count, then a comment line."
        )
    count = int(first[0])
    if count == 0:
        raise ValueError("the atom count on line 1 is 0.")
    if len(lines) < 2 + count:
        raise ValueError(
            f"line 1 announces {count} atoms, but the file ends after "
            f"{max(0, len(lines) - 2)} atom lines."
        )
    comment = lines[1].strip()

    elements, coordinates, skipped = [], [], 0
    for number, line in enumerate(lines[2:2 + count], start=3):
        tokens = line.split()
        where = f"line {number}"
        if len(tokens) < 4:
            raise ValueError(f"{where}: expected 'element x y z', found '{line.strip()}'.")
        symbol = _element(tokens[0], where)
        xyz = [_number(token, where) for token in tokens[1:4]]
        if symbol is None:
            skipped += 1
            continue
        elements.append(symbol)
        coordinates.append(xyz)

    notes = _skipped_note(skipped)
    coordinates = np.array(coordinates, dtype=float).reshape(-1, 3)
    if _XYZ_UNIT_RE.search(comment):
        coordinates = coordinates * BOHR_RADIUS_ANGSTROM
        notes.append("The comment line says the coordinates are in bohr; they were converted to Angstrom.")

    frames, position = 1, 2 + count
    while position < len(lines):
        header = lines[position].split()
        if not header or not re.fullmatch(r"\d+", header[0]):
            break
        frames += 1
        position += int(header[0]) + 2
    if frames > 1:
        notes.append(f"The file holds {frames} frames (a trajectory); the first frame was imported.")
    return StructureData(elements, coordinates, title=comment, notes=notes), False


# --- Gaussian input ----------------------------------------------------
# A Gaussian atom label: symbol or number, then optional "(...)" parameters
# and "-type-charge" MM decorations, e.g. C(Fragment=1) or C-CA--0.1.
_GAUSSIAN_ATOM_RE = re.compile(r"^\s*([^\s,(]+(?:\([^)]*\))?[^\s,]*)[\s,]*(.*)$")
_INTEGER_RE = re.compile(r"[+-]?\d+")


def _gaussian_label(label: str, where: str) -> Optional[str]:
    plain = label.split("(", 1)[0]
    parts = plain.split("-")
    if any(part.upper() == "BQ" for part in parts[1:]):
        return None  # C-Bq: a ghost atom
    return _element(parts[0], where)


def _route_options(token: str) -> List[str]:
    return [option for option in re.split(r"[=(),\s]+", token.upper())[1:] if option]


def _parse_gaussian(text: str) -> Tuple[StructureData, bool]:
    numbered = [(number, line.split("!", 1)[0].rstrip()) for number, line in
                enumerate(text.splitlines(), start=1)]
    notes: List[str] = []
    for index, (_, line) in enumerate(numbered):
        if line.strip().lower() == "--link1--":
            numbered = numbered[:index]
            notes.append("The file holds several linked jobs (--Link1--); the first job's geometry was imported.")
            break

    position = 0
    while position < len(numbered) and (not numbered[position][1].strip()
                                         or numbered[position][1].lstrip().startswith("%")):
        position += 1
    if position >= len(numbered) or not numbered[position][1].lstrip().startswith("#"):
        raise ValueError("no route section (a line starting with '#') was found, so this is not a "
                         "Gaussian input file.")
    route = []
    while position < len(numbered) and numbered[position][1].strip():
        route.append(numbered[position][1].strip())
        position += 1

    bohr = connectivity = False
    for token in route_tokens(" ".join(route).lstrip("#")):
        name, options = keyword_name(token), _route_options(token)
        if name in ("GEOM", "GEOMETRY"):
            if any(option in ("CHECK", "ALLCHECK", "STEP", "NGEOM") for option in options):
                raise ValueError("the route reads the geometry from a checkpoint (Geom=Check), so the "
                                 "file holds no coordinates. Import the matching .fchk instead.")
            connectivity = any(option in ("CONNECTIVITY", "CONNECT") for option in options)
        elif name == "UNITS":
            bohr = any(option in ("AU", "BOHR") for option in options)

    def skip_blank(position):
        while position < len(numbered) and not numbered[position][1].strip():
            position += 1
        return position

    position = skip_blank(position)
    title = []
    while position < len(numbered) and numbered[position][1].strip():
        title.append(numbered[position][1].strip())
        position += 1
    position = skip_blank(position)
    if position >= len(numbered):
        raise ValueError("the file ends before the charge and multiplicity line.")
    number, line = numbered[position]
    tokens = line.replace(",", " ").split()
    if not tokens or len(tokens) % 2 or not all(_INTEGER_RE.fullmatch(t) for t in tokens):
        raise ValueError(f"line {number}: expected the charge and spin multiplicity (e.g. '0 1'), "
                         f"found '{line.strip()}'.")
    position += 1

    elements, coordinates, file_index, skipped = [], [], {}, 0
    atom_number = 0
    while position < len(numbered) and numbered[position][1].strip():
        number, line = numbered[position]
        position += 1
        where = f"line {number}"
        match = _GAUSSIAN_ATOM_RE.match(line)
        values = re.split(r"[\s,]+", match.group(2).strip()) if match and match.group(2).strip() else []
        xyz_tokens = _gaussian_cartesian_tokens(values)
        if xyz_tokens is None:
            raise ValueError(
                f"{where}: '{line.strip()}' is not a Cartesian 'element x y z' line. Z-matrix "
                "(internal-coordinate) input cannot be imported; re-save the input with "
                "Cartesian coordinates."
            )
        atom_number += 1
        symbol = _gaussian_label(match.group(1), where)
        xyz = [_number(token, where) for token in xyz_tokens]
        if symbol is None:
            skipped += 1
            continue
        file_index[atom_number] = len(elements)
        elements.append(symbol)
        coordinates.append(xyz)
    notes.extend(_skipped_note(skipped))

    coordinates = np.array(coordinates, dtype=float).reshape(-1, 3)
    if bohr:
        coordinates = coordinates * BOHR_RADIUS_ANGSTROM
        notes.append("The route has Units=Bohr; the coordinates were converted to Angstrom.")

    bonds: List[Tuple[int, int, float]] = []
    bonds_known = False
    if connectivity:
        position = skip_blank(position)
        bonds, lines_read, odd_orders = _parse_gaussian_connectivity(numbered, position, file_index)
        bonds_known = lines_read > 0
        if odd_orders:
            notes.append(f"{odd_orders} bonds with an order other than 1, 1.5, 2 or 3 were imported as single bonds.")
        if not bonds_known:
            notes.append("The route asks for Geom=Connectivity but no connectivity section follows.")
    structure = StructureData(elements, coordinates, bonds=bonds, title=" ".join(title), notes=notes)
    return structure, bonds_known


def _is_number(token: str) -> bool:
    try:
        return math.isfinite(float(token.replace("D", "E").replace("d", "e")))
    except ValueError:
        return False


def _gaussian_cartesian_tokens(values: List[str]) -> Optional[List[str]]:
    """The x, y, z tokens of a Cartesian atom line, or ``None`` for anything else.

    ``values`` follow the atom label: ``x y z``, or ``flag x y z`` with a
    freeze code of 0 or below, either one possibly followed by ONIOM layer
    fields. A Z-matrix line references earlier atoms by positive integer
    (``H 1 0.96 2 104.5``) and is refused.
    """
    def is_reference(token):
        return re.fullmatch(r"\d+", token) is not None and int(token) > 0

    if len(values) >= 4 and _INTEGER_RE.fullmatch(values[0]) and int(values[0]) <= 0 \
            and all(_is_number(v) for v in values[1:4]):
        return values[1:4]  # freeze code: C 0 x y z / C -1 x y z
    if len(values) >= 4 and is_reference(values[0]) and is_reference(values[2]):
        return None
    if len(values) >= 3 and all(_is_number(v) for v in values[:3]):
        return values[:3]
    return None


def _parse_gaussian_connectivity(numbered, position, file_index):
    """The Geom=Connectivity block: 'atom  partner order  partner order ...'."""
    bonds, seen, lines_read, odd_orders = [], set(), 0, 0
    while position < len(numbered) and numbered[position][1].strip():
        number, line = numbered[position]
        position += 1
        tokens = line.split()
        where = f"line {number}"
        if not _INTEGER_RE.fullmatch(tokens[0]) or len(tokens) % 2 == 0:
            raise ValueError(f"{where}: a connectivity line is 'atom partner order ...', found '{line.strip()}'.")
        lines_read += 1
        atom = int(tokens[0])
        for partner_token, order_token in zip(tokens[1::2], tokens[2::2]):
            if not _INTEGER_RE.fullmatch(partner_token):
                raise ValueError(f"{where}: '{partner_token}' is not an atom number.")
            partner, order = int(partner_token), _number(order_token, where)
            if order <= 0.0 or atom not in file_index or partner not in file_index:
                continue
            pair = tuple(sorted((file_index[atom], file_index[partner])))
            if pair[0] == pair[1] or pair in seen:
                continue
            if order not in SUPPORTED_ORDERS:
                odd_orders += 1
                order = SINGLE
            seen.add(pair)
            bonds.append((pair[0], pair[1], order))
    return bonds, lines_read, odd_orders


# --- Gaussian formatted checkpoint -------------------------------------
def _fchk_array(lines: List[str], name: str, kind: str) -> np.ndarray:
    """One array section of an fchk, read without parsing the rest of the file.

    Same layout as gaussianWavefunction.read_fchk_sections (name in columns
    1-43, type in 44, ``N=`` count), which parses every section -- far too
    slow for a large checkpoint when only the geometry is wanted.
    """
    for index, line in enumerate(lines):
        if line[:43].strip() != name or len(line) < 45 or line[43] != kind:
            continue
        rest = line[44:].split()
        if not rest or not rest[0].startswith("N="):
            break
        count = int(rest[0][2:]) if len(rest[0]) > 2 else int(rest[1])
        tokens: List[str] = []
        for value_line in lines[index + 1:]:
            if value_line[:1] != " " or len(tokens) >= count:
                break
            tokens.extend(value_line.split())
        if len(tokens) != count:
            raise ValueError(f"section '{name}' declares {count} values but holds {len(tokens)}; "
                             "the file looks truncated.")
        try:
            return np.array([int(t) if kind == "I" else float(t.replace("D", "E")) for t in tokens])
        except ValueError:
            raise ValueError(f"section '{name}' contains a value that is not a number.") from None
    raise ValueError(f"the '{name}' section is missing; this is not a Gaussian formatted checkpoint.")


def _parse_fchk(text: str) -> Tuple[StructureData, bool]:
    lines = text.splitlines()
    if len(lines) < 3:
        raise ValueError("the file is too short to be a Gaussian formatted checkpoint.")
    numbers = _fchk_array(lines, "Atomic numbers", "I").astype(int)
    coordinates = _fchk_array(lines, "Current cartesian coordinates", "R")
    if coordinates.size != 3 * numbers.size:
        raise ValueError(f"{numbers.size} atomic numbers but {coordinates.size} coordinates.")
    symbols = {z: s for s, z in AtomicData.atomic_numbers.items()}
    elements, kept, skipped = [], [], 0
    for index, z in enumerate(numbers):
        if z == 0:
            skipped += 1  # ghost atom
            continue
        if z not in symbols:
            raise ValueError(f"atom {index + 1} has atomic number {z}, which is not an element.")
        elements.append(symbols[z])
        kept.append(index)
    coordinates = coordinates.reshape(-1, 3)[kept] * BOHR_RADIUS_ANGSTROM
    return StructureData(elements, coordinates, title=lines[0].strip(), notes=_skipped_note(skipped)), False


# --- MDL molfile and SD file ------------------------------------------
def _parse_mol(text: str, lines: Optional[List[str]] = None) -> Tuple[StructureData, bool]:
    lines = text.splitlines() if lines is None else lines
    if len(lines) < 4:
        raise ValueError("the file is too short: a molfile has three header lines and a counts line.")
    counts = lines[3]
    if "V3000" in counts.upper():
        raise ValueError("V3000 molfiles are not supported; save the structure as a V2000 molfile.")
    try:
        atom_count, bond_count = int(counts[0:3]), int(counts[3:6])
    except ValueError:
        raise ValueError(f"line 4 should be the counts line (atom and bond counts in columns 1-6), "
                         f"found '{counts.strip()}'.") from None
    if atom_count <= 0:
        raise ValueError("the counts line declares no atoms.")
    if len(lines) < 4 + atom_count + bond_count:
        raise ValueError(f"the counts line declares {atom_count} atoms and {bond_count} bonds, "
                         "but the file ends before them.")

    elements, coordinates, file_index, skipped = [], [], {}, 0
    for offset, line in enumerate(lines[4:4 + atom_count]):
        where = f"line {5 + offset}"
        try:
            xyz, symbol_token = [float(line[0:10]), float(line[10:20]), float(line[20:30])], line[31:34].strip()
        except ValueError:
            tokens = line.split()  # not in the fixed columns: try free format
            if len(tokens) < 4:
                raise ValueError(f"{where}: expected an atom line 'x y z symbol ...', found '{line.strip()}'.") from None
            xyz, symbol_token = [_number(t, where) for t in tokens[:3]], tokens[3]
        symbol = _element(symbol_token, where)
        if symbol is None:
            skipped += 1
            continue
        file_index[offset + 1] = len(elements)
        elements.append(symbol)
        coordinates.append(xyz)

    bonds, query_bonds = [], 0
    for offset, line in enumerate(lines[4 + atom_count:4 + atom_count + bond_count]):
        where = f"line {5 + atom_count + offset}"
        try:
            first, second, kind = int(line[0:3]), int(line[3:6]), int(line[6:9])
        except ValueError:
            raise ValueError(f"{where}: expected a bond line 'atom atom type', found '{line.strip()}'.") from None
        for atom in (first, second):
            if not 1 <= atom <= atom_count:
                raise ValueError(f"{where}: the bond refers to atom {atom}, but there are {atom_count} atoms.")
        if first not in file_index or second not in file_index:
            continue
        order = MOL_BOND_ORDERS.get(kind)
        if order is None:
            query_bonds += 1
            order = SINGLE
        bonds.append((file_index[first], file_index[second], order))

    notes = _skipped_note(skipped)
    if query_bonds:
        notes.append(f"{query_bonds} query bonds (types 5-8) were imported as single bonds.")
    structure = StructureData(elements, np.array(coordinates, dtype=float), bonds=bonds,
                              title=lines[0].strip(), notes=notes)
    # A molfile with atoms but no bond block (as some converters write) has no
    # connectivity at all, rather than a set of isolated atoms.
    return structure, bond_count > 0 or atom_count == 1


def _parse_sdf(text: str) -> Tuple[StructureData, bool]:
    lines = text.splitlines()
    records = [index for index, line in enumerate(lines) if line.strip() == "$$$$"]
    end = records[0] if records else len(lines)
    structure, bonds_known = _parse_mol(text, lines[:end])
    if records:
        trailing = any(line.strip() for line in lines[records[-1] + 1:])
        total = len(records) + (1 if trailing else 0)
        if total > 1:
            structure.notes.append(f"The file holds {total} records; the first one was imported.")
    return structure, bonds_known


# --- Tripos MOL2 --------------------------------------------------------
def _parse_mol2(text: str) -> Tuple[StructureData, bool]:
    sections: Dict[str, List[Tuple[int, str]]] = {}
    current = None
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if stripped.startswith("@<TRIPOS>"):
            current = stripped[len("@<TRIPOS>"):].upper()
            if current in sections:
                current = None  # a second molecule: only the first is read
            else:
                sections[current] = []
            continue
        if current is not None and stripped and not stripped.startswith("#"):
            sections[current].append((number, stripped))
    if "ATOM" not in sections:
        raise ValueError("there is no @<TRIPOS>ATOM section.")

    title = sections.get("MOLECULE", [(0, "")])[0][1]
    elements, coordinates, file_index, skipped = [], [], {}, 0
    for number, line in sections["ATOM"]:
        tokens = line.split()
        where = f"line {number}"
        if len(tokens) < 6:
            raise ValueError(f"{where}: expected 'id name x y z type ...', found '{line}'.")
        xyz = [_number(t, where) for t in tokens[2:5]]
        type_symbol = tokens[5].split(".")[0]
        if type_symbol.upper() in DUMMY_LABELS:
            symbol = None
        else:
            symbol = element_symbol(type_symbol)
            if symbol is None:  # a generic type ("Any", "Hal"): fall back on the atom name
                match = re.match(r"[A-Za-z]{1,2}", tokens[1])
                symbol = element_symbol(match.group(0)) if match else None
            if symbol is None:
                raise ValueError(f"{where}: neither the atom type '{tokens[5]}' nor the name "
                                 f"'{tokens[1]}' identifies an element.")
        if symbol is None:
            skipped += 1
            continue
        file_index[tokens[0]] = len(elements)
        elements.append(symbol)
        coordinates.append(xyz)

    atom_ids = {line.split()[0] for _, line in sections["ATOM"]}
    bonds = []
    for number, line in sections.get("BOND", []):
        tokens = line.split()
        if len(tokens) < 4:
            raise ValueError(f"line {number}: expected 'id atom atom type', found '{line}'.")
        kind = tokens[3].lower()
        if kind == "nc":
            continue  # "not connected"
        for atom in tokens[1:3]:
            if atom not in atom_ids:
                raise ValueError(f"line {number}: the bond refers to atom {atom}, which is not in the ATOM section.")
        if tokens[1] not in file_index or tokens[2] not in file_index:
            continue
        order = MOL2_BOND_ORDERS.get(kind)
        if order is None:
            raise ValueError(f"line {number}: unknown bond type '{tokens[3]}'.")
        bonds.append((file_index[tokens[1]], file_index[tokens[2]], order))

    structure = StructureData(elements, np.array(coordinates, dtype=float), bonds=bonds,
                              title=title, notes=_skipped_note(skipped))
    return structure, bool(sections.get("BOND")) or len(elements) == 1


# --- PDB ---------------------------------------------------------------
def _pdb_element(line: str, where: str) -> Optional[str]:
    """Columns 77-78; otherwise the atom name, whose first two columns hold
    the element right-justified (' CA ' is an alpha carbon, 'CA  ' calcium)."""
    token = line[76:78].strip()
    if token:
        return _element(token, where)
    name = line[12:16]
    if name[:1] in (" ", "") or name[:1].isdigit():
        letter = name[1:2].strip()
        if not letter:
            raise ValueError(f"{where}: the atom has neither an element (columns 77-78) nor a name.")
        return _element(letter, where)
    two = element_symbol(name[:2].strip())
    return two if two is not None else _element(name[:1], where)


def _parse_pdb(text: str) -> Tuple[StructureData, bool]:
    elements, coordinates, serial_index, skipped = [], [], {}, 0
    conect: List[Tuple[int, int]] = []
    notes: List[str] = []
    title = ""
    models_ended = alternates = 0
    for number, line in enumerate(text.splitlines(), start=1):
        record = line[:6].strip().upper()
        where = f"line {number}"
        if record in ("TITLE", "COMPND") and not title:
            title = line[10:].strip()
        elif record == "ENDMDL":
            models_ended += 1
        elif record in ("ATOM", "HETATM") and not models_ended:
            if line[16:17] not in (" ", "", "A"):
                alternates += 1
                continue
            try:
                xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            except ValueError:
                raise ValueError(f"{where}: the coordinates in columns 31-54 are not numbers.") from None
            symbol = _pdb_element(line, where)
            serial = line[6:11].strip()
            if symbol is None:
                skipped += 1
                continue
            if serial:
                serial_index[serial] = len(elements)
            elements.append(symbol)
            coordinates.append(xyz)
        elif record == "CONECT":
            fields = [line[start:start + 5].strip() for start in (6, 11, 16, 21, 26)]
            if not fields[0]:
                raise ValueError(f"{where}: a CONECT record without an atom serial number.")
            conect.extend((fields[0], partner) for partner in fields[1:] if partner)

    if models_ended:
        notes.append("The file holds several models; the first model was imported.")
    if alternates:
        notes.append(f"{alternates} alternate-location atoms were skipped (the first location was kept).")
    notes.extend(_skipped_note(skipped))
    coordinates = np.array(coordinates, dtype=float).reshape(-1, 3)

    pairs = set()
    for first, second in conect:
        if first in serial_index and second in serial_index and serial_index[first] != serial_index[second]:
            pairs.add(tuple(sorted((serial_index[first], serial_index[second]))))
    if not pairs:
        return StructureData(elements, coordinates, title=title, notes=notes), False

    # CONECT usually lists only the hetero groups: atoms it never mentions
    # (a protein's standard residues) get distance-detected bonds instead.
    bonds = [(i, j, SINGLE) for i, j in sorted(pairs)]
    listed = {atom for pair in pairs for atom in pair}
    detected = [bond for bond in detect_bonds(elements, coordinates)
                if bond[0] not in listed and bond[1] not in listed]
    notes.append("PDB files store connectivity but no bond orders: every bond was imported as single.")
    if detected:
        notes.append(f"{len(detected)} bonds between atoms without CONECT records were detected from distances.")
    return StructureData(elements, coordinates, bonds=bonds + detected, title=title, notes=notes), True


# ----------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------
def structure_to_text(structure: StructureData, format_key: str, name: str = "molecule") -> str:
    """The structure as the text of a ``format_key`` file; ``name`` is the file stem."""
    writers = {
        "XYZ": lambda: write_xyz_text(structure),
        "Gaussian": lambda: write_gaussian_text(structure, name),
        "MOL": lambda: write_mol_text(structure),
        "SDF": lambda: write_mol_text(structure) + "$$$$\n",
        "MOL2": lambda: write_mol2_text(structure),
        "PDB": lambda: write_pdb_text(structure),
    }
    writer = writers.get(format_key)
    if writer is None:
        raise ValueError(f"Structures cannot be saved as {format_key}. "
                         f"Supported formats: {supported_formats_text('write')}.")
    if not structure.elements:
        raise ValueError("There is no structure to save.")
    return writer()


def export_notes(structure: StructureData, format_key: str) -> List[str]:
    """What a saved file does not keep, for the confirmation dialog."""
    has_orders = any(order != SINGLE for _, _, order in structure.bonds)
    if format_key == "PDB" and has_orders:
        return ["PDB files store connectivity but no bond orders: every bond reads back as single."]
    if format_key == "Gaussian":
        charge, multiplicity = gaussian_charge_and_multiplicity(structure)
        return [f"The Gaussian input assumes charge {charge} and multiplicity {multiplicity}, and its route "
                f"({gaussian_route()}) is a template: check both before running it."]
    return []


def write_xyz_text(structure: StructureData) -> str:
    lines = [str(structure.atom_count), structure.title]
    for element, (x, y, z) in zip(structure.elements, structure.coordinates):
        lines.append(f"{element:<3s} {x: .8f} {y: .8f} {z: .8f}")
    return "\n".join(lines) + "\n"


def gaussian_route() -> str:
    return f"#P {DEFAULT_METHOD}/{DEFAULT_BASIS} Opt Geom=Connectivity"


def gaussian_charge_and_multiplicity(structure: StructureData) -> Tuple[int, int]:
    """Neutral, lowest spin the electron count allows (1 if even, 2 if odd)."""
    return 0, allowed_multiplicities(nuclear_charge(structure.elements))[0]


def write_gaussian_text(structure: StructureData, name: str = "molecule") -> str:
    charge, multiplicity = gaussian_charge_and_multiplicity(structure)
    stem = re.sub(r"[^A-Za-z0-9._-]+", "", name) or "molecule"
    lines = [f"%Chk={stem}.chk", gaussian_route(), "",
             f"Built with the gQTEA molecular design tool. Charge {charge} and multiplicity "
             f"{multiplicity} are assumed, check them before running.", "",
             f"{charge} {multiplicity}"]
    lines.extend(format_geometry(list(zip(structure.elements, *structure.coordinates.T))))
    lines.append("")
    partners: List[List[Tuple[int, float]]] = [[] for _ in structure.elements]
    for i, j, order in structure.bonds:
        partners[i].append((j, order))
    for atom, bonded in enumerate(partners):
        lines.append(f" {atom + 1}" + "".join(f" {j + 1} {order:.1f}" for j, order in sorted(bonded)))
    lines.append("")
    return "\n".join(lines) + "\n"


def write_mol_text(structure: StructureData) -> str:
    if structure.atom_count > MOL_MAX_COUNT or len(structure.bonds) > MOL_MAX_COUNT:
        raise ValueError(f"A V2000 molfile holds at most {MOL_MAX_COUNT} atoms and {MOL_MAX_COUNT} bonds; "
                         f"this structure has {structure.atom_count} atoms and {len(structure.bonds)} bonds.")
    if np.abs(structure.coordinates).max() > MOL_MAX_COORDINATE:
        raise ValueError(f"A molfile coordinate cannot exceed {MOL_MAX_COORDINATE} A in magnitude.")
    lines = [structure.title[:80], f"{'  gQTEA':<20s}3D", "",
             f"{structure.atom_count:3d}{len(structure.bonds):3d}" + "  0" * 8 + "999 V2000"]
    for element, (x, y, z) in zip(structure.elements, structure.coordinates):
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {element:<3s} 0" + "  0" * 11)
    for i, j, order in structure.bonds:
        lines.append(f"{i + 1:3d}{j + 1:3d}{MOL_BOND_TYPES[order]:3d}" + "  0" * 4)
    lines.append("M  END")
    return "\n".join(lines) + "\n"


def sybyl_atom_type(element: str, orders: Sequence[float]) -> str:
    """A SYBYL atom type from the element and its bond orders (C.ar, C.2, O.3 ...)."""
    if element == "C":
        if RESONANCE in orders:
            return "C.ar"
        if TRIPLE in orders or list(orders).count(DOUBLE) >= 2:
            return "C.1"
        return "C.2" if DOUBLE in orders else "C.3"
    if element == "N":
        if RESONANCE in orders:
            return "N.ar"
        if TRIPLE in orders:
            return "N.1"
        if DOUBLE in orders:
            return "N.2"
        return "N.4" if len(orders) >= 4 else "N.3"
    if element in ("O", "S"):
        return f"{element}.2" if DOUBLE in orders or RESONANCE in orders else f"{element}.3"
    if element == "P":
        return "P.3"
    return element


def write_mol2_text(structure: StructureData) -> str:
    lines = ["@<TRIPOS>MOLECULE", structure.title or "molecule",
             f"{structure.atom_count} {len(structure.bonds)} 1 0 0", "SMALL", "NO_CHARGES", "",
             "@<TRIPOS>ATOM"]
    orders = structure.neighbour_orders()
    for index, (element, (x, y, z)) in enumerate(zip(structure.elements, structure.coordinates), start=1):
        atom_type = sybyl_atom_type(element, orders[index - 1])
        lines.append(f"{index:7d} {element + str(index):<8s} {x:10.4f} {y:10.4f} {z:10.4f} "
                     f"{atom_type:<6s} 1 UNL1 0.0000")
    lines.append("@<TRIPOS>BOND")
    for index, (i, j, order) in enumerate(structure.bonds, start=1):
        lines.append(f"{index:6d} {i + 1:5d} {j + 1:5d} {MOL2_BOND_TYPES[order]}")
    return "\n".join(lines) + "\n"


def _pdb_atom_name(element: str, index: int) -> str:
    """Four columns, a one-letter element starting in column 14 as the standard asks."""
    label = f"{element}{index}"
    if len(element) == 1 and len(label) <= 3:
        return f" {label:<3s}"
    return f"{label[:4]:<4s}" if len(label) <= 4 else f"{element:<4s}"


def write_pdb_text(structure: StructureData) -> str:
    if structure.atom_count > PDB_MAX_ATOMS:
        raise ValueError(f"A PDB file holds at most {PDB_MAX_ATOMS} atoms; this structure has "
                         f"{structure.atom_count}.")
    if np.abs(structure.coordinates).max() > PDB_MAX_COORDINATE:
        raise ValueError(f"A PDB coordinate cannot exceed {PDB_MAX_COORDINATE} A in magnitude.")
    lines = [f"COMPND    {structure.title[:70]}".rstrip()]
    for index, (element, (x, y, z)) in enumerate(zip(structure.elements, structure.coordinates), start=1):
        lines.append(f"HETATM{index:5d} {_pdb_atom_name(element, index)} UNL A   1    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element.upper():>2s}  ")
    partners: List[List[int]] = [[] for _ in structure.elements]
    for i, j, _ in structure.bonds:
        partners[i].append(j)
        partners[j].append(i)
    for atom, bonded in enumerate(partners):
        bonded = sorted(bonded)
        for start in range(0, len(bonded), 4):
            lines.append(f"CONECT{atom + 1:5d}" + "".join(f"{j + 1:5d}" for j in bonded[start:start + 4]))
    lines.append("END")
    return "\n".join(lines) + "\n"
