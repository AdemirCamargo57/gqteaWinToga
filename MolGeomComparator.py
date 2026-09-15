"""Comparison of molecular geometric parameters between two simulations.

This tool answers one question: how does the environment shift the geometry of a
molecule? It takes the combined parameter files written by the three all-*
analysis tools ([allBondAnalysis.py], [allAnglesAnalysis.py],
[allDihedralAnalysis.py]) for the *same* molecule simulated twice -- typically
isolated in vacuum and solvated in a water box -- matches the parameters that
correspond to each other, and reports each one side by side with the difference.

The matching problem is an index-remapping problem. A solvated run is analysed
with ``solute_atom_indices``, so its rows carry the atoms' *global* indices in
the box (the oxygens of the reference molecule sit at 295-297, for instance),
while the isolated run numbers the very same atoms 29-31. Every row is therefore
canonicalised to a **solute-local** index -- the atom's rank in the sorted
solute list, which is the identity in the isolated file -- before matching. The
remap is verified against the element of every shared atom rather than assumed:
a disagreement aborts the run naming the offending atom, because a silently
wrong mapping produces plausible-looking nonsense.

The three input formats differ only in how many atoms each parameter names and
what its value column is called, so one table-driven reader serves all of them,
locating columns **by name** from the file's own ``# row`` header line.
"""
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import toga
from toga.style import Pack
from toga.style.pack import LEFT

from help import HelpGqteaWin


# --------------------------------------------------------------------------- #
# Parameter kinds: the three combined-file formats, described declaratively     #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ParameterKind:
    """One of the three combined-parameter file formats.

    ``periodic`` marks a coordinate whose values live on a circle, so that
    differences must be wrapped into (-180, 180]; only the signed dihedral is.
    """

    key: str
    title: str
    label: str
    n_atoms: int
    average_column: str
    std_column: str
    occurrence_column: str
    frames_column: str
    unit: str
    periodic: bool

    @property
    def atom_columns(self) -> Tuple[str, ...]:
        return ATOM_COLUMN_NAMES[: self.n_atoms]

    @property
    def element_columns(self) -> Tuple[str, ...]:
        return ELEMENT_COLUMN_NAMES[: self.n_atoms]

    @property
    def required_columns(self) -> Tuple[str, ...]:
        return self.atom_columns + self.element_columns + (
            self.average_column,
            self.std_column,
            self.occurrence_column,
        )

    @property
    def unit_symbol(self) -> str:
        return "A" if self.unit == "angstrom" else "deg"


ATOM_COLUMN_NAMES = ("atom_i", "atom_j", "atom_k", "atom_l")
ELEMENT_COLUMN_NAMES = ("element_i", "element_j", "element_k", "element_l")

BOND_KIND = ParameterKind(
    key="bond",
    title="gQTEA All Bond Distance Analysis",
    label="bond distance",
    n_atoms=2,
    average_column="average",
    std_column="standard_deviation",
    occurrence_column="occurrence_fraction",
    frames_column="frames_bonded",
    unit="angstrom",
    periodic=False,
)

ANGLE_KIND = ParameterKind(
    key="angle",
    title="gQTEA All Bond Angle Analysis",
    label="bond angle",
    n_atoms=3,
    average_column="average_angle_degrees",
    std_column="standard_deviation_degrees",
    occurrence_column="occurrence_fraction",
    frames_column="frames_present",
    unit="degrees",
    periodic=False,
)

DIHEDRAL_KIND = ParameterKind(
    key="dihedral",
    title="gQTEA All Dihedral Angle Analysis",
    label="dihedral angle",
    n_atoms=4,
    average_column="average_dihedral_degrees",
    std_column="standard_deviation_degrees",
    occurrence_column="occurrence_fraction",
    frames_column="frames_present",
    unit="degrees",
    periodic=True,
)

PARAMETER_KINDS = (BOND_KIND, ANGLE_KIND, DIHEDRAL_KIND)


def detect_parameter_kind(title_line: str) -> ParameterKind:
    """Map a file's leading ``# gQTEA ...`` title line onto its format."""
    title = title_line.lstrip("#").strip()
    for kind in PARAMETER_KINDS:
        if title == kind.title:
            return kind

    known = "\n".join(f"  - {kind.title}" for kind in PARAMETER_KINDS)
    raise ValueError(
        f"Unrecognised parameter file: its first line reads '{title}'.\n"
        f"Expected one of:\n{known}"
    )


# --------------------------------------------------------------------------- #
# Parsed representation                                                         #
# --------------------------------------------------------------------------- #
@dataclass
class ParameterRow:
    """One parameter (a bond, angle or dihedral) as read from one file.

    ``global_atoms`` are the 1-based indices printed in the file; ``local_atoms``
    are the same atoms renumbered within the analysed molecule, which is the
    frame the two files are matched in.
    """

    source_row: int
    global_atoms: Tuple[int, ...]
    local_atoms: Tuple[int, ...]
    elements: Tuple[str, ...]
    average: float
    std_dev: float
    occurrence: float
    frames: int


@dataclass
class GeometryParameterFile:
    """A parsed combined-parameter file."""

    path: str
    kind: ParameterKind
    metadata: Dict[str, str]
    rows: List[ParameterRow]
    local_to_global: Dict[int, int]
    elements_by_local: Dict[int, str]

    @property
    def name(self) -> str:
        return os.path.basename(self.path)

    @property
    def atom_scope(self) -> str:
        return self.metadata.get("atom_scope", "all_atoms")

    @property
    def frames_used(self) -> Optional[int]:
        try:
            return int(self.metadata["frames_used"])
        except (KeyError, ValueError):
            return None

    @property
    def num_local_atoms(self) -> int:
        """How many atoms the analysed molecule has.

        With a solute scope this is the length of the header's atom list, so an
        atom that happens to form no bond still counts; otherwise the largest
        index any row mentions is the best available estimate.
        """
        return max(self.local_to_global) if self.local_to_global else 0

    def remapped_atoms(self) -> Dict[int, int]:
        """Local -> global for the atoms whose numbering actually changed."""
        return {
            local: global_index
            for local, global_index in sorted(self.local_to_global.items())
            if local != global_index
        }


def _parse_solute_map(metadata: Dict[str, str], path: str) -> Optional[Dict[int, int]]:
    """Build local -> global from ``# solute_atom_indices``.

    The all-* writers print this list sorted, and a solute atom's rank in it is
    exactly its position in a file that analysed the molecule on its own, so the
    rank is the local index.
    """
    if metadata.get("atom_scope") != "solute_atoms":
        return None

    raw = metadata.get("solute_atom_indices", "").split()
    if not raw:
        raise ValueError(
            f"{os.path.basename(path)}: the header says 'atom_scope solute_atoms' "
            "but no '# solute_atom_indices' line lists them."
        )

    try:
        indices = sorted({int(value) for value in raw})
    except ValueError:
        raise ValueError(
            f"{os.path.basename(path)}: '# solute_atom_indices' must be a list of "
            "whole numbers."
        )

    return {local: global_index for local, global_index in enumerate(indices, start=1)}


def read_parameter_file(path: str) -> GeometryParameterFile:
    """Read one combined-parameter file into a `GeometryParameterFile`.

    Raises ``ValueError`` -- naming the file and, where it helps, the offending
    line -- for anything the comparator cannot trust.
    """
    path = os.path.abspath(path)
    name = os.path.basename(path)

    try:
        with open(path, "r") as handle:
            lines = handle.read().splitlines()
    except OSError as exc:
        raise ValueError(f"Could not read {name}: {exc}")

    title_line = next((line for line in lines if line.strip()), "")
    if not title_line.lstrip().startswith("#"):
        raise ValueError(
            f"{name} does not look like a gQTEA parameter file: it must start "
            "with a '# gQTEA All ... Analysis' title line."
        )

    try:
        kind = detect_parameter_kind(title_line)
    except ValueError as exc:
        raise ValueError(f"{name}: {exc}")

    metadata: Dict[str, str] = {}
    columns: Optional[List[str]] = None
    column_index: Dict[str, int] = {}
    rows: List[ParameterRow] = []
    solute_map: Optional[Dict[int, int]] = None
    local_to_global: Dict[int, int] = {}
    elements_by_local: Dict[int, str] = {}

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("#"):
            content = line.lstrip("#").strip()
            if not content:
                continue
            key, _, value = content.partition(" ")
            if key == "row":
                # The column-name header; everything after it is data. "row" is
                # itself the name of the leading row-number column.
                columns = ["row"] + value.split()
                column_index = {column: position for position, column in enumerate(columns)}
                missing = [c for c in kind.required_columns if c not in column_index]
                if missing:
                    raise ValueError(
                        f"{name} is missing the column(s) {', '.join(missing)} from its "
                        f"'# row' header, so it cannot be read as a {kind.label} file."
                    )
                solute_map = _parse_solute_map(metadata, path)
                if solute_map is not None:
                    # Seed the molecule from the header list, so an atom that
                    # forms no bond still counts towards its size.
                    local_to_global.update(solute_map)
            else:
                metadata[key] = value.strip()
            continue

        if columns is None:
            raise ValueError(
                f"{name}: found a data row on line {line_number} before the "
                "'# row' column header. The file looks truncated or edited."
            )

        fields = line.split()
        if len(fields) != len(columns):
            raise ValueError(
                f"{name}, line {line_number}: expected {len(columns)} columns "
                f"({' '.join(columns)}) but found {len(fields)}."
            )

        row = _build_row(
            kind=kind,
            fields=fields,
            column_index=column_index,
            solute_map=solute_map,
            name=name,
            line_number=line_number,
            source_row=len(rows) + 1,
        )
        rows.append(row)

        for local, global_index, element in zip(row.local_atoms, row.global_atoms, row.elements):
            local_to_global[local] = global_index
            known = elements_by_local.setdefault(local, element)
            if known != element:
                raise ValueError(
                    f"{name}, line {line_number}: atom {local} is given as "
                    f"'{known}' on an earlier row but as '{element}' here. "
                    "The file is inconsistent."
                )

    if columns is None:
        raise ValueError(
            f"{name} has no '# row' column header, so its columns cannot be "
            "identified. Re-run the analysis tool that produced it."
        )

    if not rows:
        raise ValueError(f"{name} contains no parameter rows.")

    return GeometryParameterFile(
        path=path,
        kind=kind,
        metadata=metadata,
        rows=rows,
        local_to_global=local_to_global,
        elements_by_local=elements_by_local,
    )


def _build_row(
    *,
    kind: ParameterKind,
    fields: Sequence[str],
    column_index: Dict[str, int],
    solute_map: Optional[Dict[int, int]],
    name: str,
    line_number: int,
    source_row: int,
) -> ParameterRow:
    """Turn one split data line into a `ParameterRow`, remapping its atoms."""

    def value_of(column: str) -> str:
        return fields[column_index[column]]

    try:
        global_atoms = tuple(int(value_of(column)) for column in kind.atom_columns)
    except ValueError:
        raise ValueError(
            f"{name}, line {line_number}: the atom index columns must be whole numbers."
        )

    try:
        average = float(value_of(kind.average_column))
        std_dev = float(value_of(kind.std_column))
        occurrence = float(value_of(kind.occurrence_column))
    except ValueError:
        raise ValueError(
            f"{name}, line {line_number}: the value columns must be numbers."
        )

    frames = 0
    if kind.frames_column in column_index:
        try:
            frames = int(value_of(kind.frames_column))
        except ValueError:
            frames = 0

    elements = tuple(value_of(column) for column in kind.element_columns)

    if solute_map is None:
        local_atoms = global_atoms
    else:
        global_to_local = {g: l for l, g in solute_map.items()}
        try:
            local_atoms = tuple(global_to_local[atom] for atom in global_atoms)
        except KeyError as exc:
            raise ValueError(
                f"{name}, line {line_number}: atom {exc.args[0]} is not listed in "
                "'# solute_atom_indices', so it cannot be placed within the "
                "molecule. The header and the rows disagree."
            )

    return ParameterRow(
        source_row=source_row,
        global_atoms=global_atoms,
        local_atoms=local_atoms,
        elements=elements,
        average=average,
        std_dev=std_dev,
        occurrence=occurrence,
        frames=frames,
    )


# --------------------------------------------------------------------------- #
# Pure helpers: identity, circular differences, naming                          #
# --------------------------------------------------------------------------- #
def _canonical_permutation(kind: ParameterKind, atoms: Sequence[int]) -> Tuple[int, ...]:
    """The ordering that makes equivalent atom tuples compare equal.

    A bond is undirected; an angle keeps its vertex but may swap its arms; a
    dihedral read backwards has the same signed value, so i-j-k-l and l-k-j-i
    are one parameter.
    """
    if kind.n_atoms == 2:
        return (0, 1) if atoms[0] <= atoms[1] else (1, 0)
    if kind.n_atoms == 3:
        return (0, 1, 2) if atoms[0] <= atoms[2] else (2, 1, 0)
    return (0, 1, 2, 3) if tuple(atoms) <= tuple(reversed(atoms)) else (3, 2, 1, 0)


def canonical_identity(kind: ParameterKind, atoms: Sequence[int]) -> Tuple[int, ...]:
    """The key two files are matched on: atom indices in canonical order."""
    return tuple(atoms[position] for position in _canonical_permutation(kind, atoms))


def canonical_elements(
    kind: ParameterKind, atoms: Sequence[int], elements: Sequence[str]
) -> Tuple[str, ...]:
    """The elements permuted to follow `canonical_identity`."""
    return tuple(elements[position] for position in _canonical_permutation(kind, atoms))


def angular_difference(first: float, second: float) -> float:
    """``second - first`` for a signed angle, wrapped into (-180, 180].

    Without the wrap a dihedral that moves from +179 deg to -179 deg -- a 2 deg
    shift -- would be reported as -358 deg.
    """
    difference = (second - first + 180.0) % 360.0 - 180.0
    if difference <= -180.0:
        difference += 360.0
    return difference


def sanitize_label(text: str, fallback: str = "file") -> str:
    """Turn a user-typed label into a token safe for a column name."""
    cleaned = re.sub(r"[^A-Za-z0-9-]+", "_", text.strip()).strip("_")
    return cleaned or fallback


def default_output_name(kind: ParameterKind) -> str:
    return f"geometry_comparison_{kind.key}.txt"


@dataclass
class MatchedParameter:
    """One geometric parameter found in both simulations."""

    local_atoms: Tuple[int, ...]
    elements: Tuple[str, ...]
    average_1: float
    std_1: float
    occurrence_1: float
    source_row_1: int
    average_2: float
    std_2: float
    occurrence_2: float
    source_row_2: int
    difference: float

    @property
    def atom_label(self) -> str:
        """A readable identity such as ``C1-C2`` or ``H16-C3-H17``."""
        return "-".join(
            f"{element}{atom}" for element, atom in zip(self.elements, self.local_atoms)
        )


# --------------------------------------------------------------------------- #
# The comparison                                                                #
# --------------------------------------------------------------------------- #
class MolGeomCalculator:
    """
    All functions required to compare two molecular structural parameters are provided here.
    """

    MAX_REPORTED_MISMATCHES = 5

    def __init__(
        self,
        file_1: str,
        file_2: str,
        label_1: str = "isolated",
        label_2: str = "solvated",
        min_occurrence: float = 0.0,
    ) -> None:
        self.file_1 = file_1
        self.file_2 = file_2
        self.label_1 = sanitize_label(label_1, fallback="file_1")
        self.label_2 = sanitize_label(label_2, fallback="file_2")
        self.min_occurrence = min_occurrence

        self.parsed_1: Optional[GeometryParameterFile] = None
        self.parsed_2: Optional[GeometryParameterFile] = None
        self.matched: List[MatchedParameter] = []
        self.only_in_file_1: List[ParameterRow] = []
        self.only_in_file_2: List[ParameterRow] = []
        self.filtered_out_1 = 0
        self.filtered_out_2 = 0

    # -- public API -------------------------------------------------------- #
    @property
    def kind(self) -> Optional[ParameterKind]:
        return self.parsed_1.kind if self.parsed_1 is not None else None

    def run(self) -> "MolGeomCalculator":
        """Read both files, check they describe the same molecule, and match.

        Raises ``ValueError`` with a specific message for anything that would
        make the comparison meaningless.
        """
        self._validate_threshold()
        self.parsed_1 = read_parameter_file(self.file_1)
        self.parsed_2 = read_parameter_file(self.file_2)
        self._validate_same_parameter_type()
        self._validate_same_molecule()
        self._compare()
        return self

    def largest_shifts(self, count: int = 10) -> List[MatchedParameter]:
        """The parameters the second environment changed most."""
        return sorted(self.matched, key=lambda match: -abs(match.difference))[:count]

    # -- validation -------------------------------------------------------- #
    def _validate_threshold(self) -> None:
        if self.min_occurrence < 0.0 or self.min_occurrence > 1.0:
            raise ValueError(
                "The minimum occurrence fraction must lie between 0 and 1 "
                f"(got {self.min_occurrence:g})."
            )

    def _validate_same_parameter_type(self) -> None:
        if self.parsed_1.kind is self.parsed_2.kind:
            return
        raise ValueError(
            f"The two files hold different parameter types: "
            f"{self.parsed_1.name} is a {self.parsed_1.kind.label} file and "
            f"{self.parsed_2.name} is a {self.parsed_2.kind.label} file. "
            "Compare two files produced by the same analysis tool."
        )

    def _validate_same_molecule(self) -> None:
        """Check the two files describe the same molecule, atom by atom.

        The index remap is only as good as this check: if it is wrong, the
        elements stop agreeing, and matching the wrong atoms together would
        produce plausible-looking nonsense rather than an obvious failure.
        """
        size_1 = self.parsed_1.num_local_atoms
        size_2 = self.parsed_2.num_local_atoms
        if size_1 != size_2:
            raise ValueError(
                f"The two files describe molecules of different size: "
                f"{self.parsed_1.name} has {size_1} atoms and "
                f"{self.parsed_2.name} has {size_2}. "
                "Both runs must analyse the same molecule -- check the solute "
                "atom indices used for the solvated run."
            )

        mismatches = [
            f"  atom {local}: '{element}' in {self.parsed_1.name} but "
            f"'{self.parsed_2.elements_by_local[local]}' in {self.parsed_2.name}"
            for local, element in sorted(self.parsed_1.elements_by_local.items())
            if local in self.parsed_2.elements_by_local
            and self.parsed_2.elements_by_local[local] != element
        ]
        if not mismatches:
            return

        shown = "\n".join(mismatches[: self.MAX_REPORTED_MISMATCHES])
        if len(mismatches) > self.MAX_REPORTED_MISMATCHES:
            shown += f"\n  ... and {len(mismatches) - self.MAX_REPORTED_MISMATCHES} more"
        raise ValueError(
            "The two files do not describe the same molecule -- these atoms "
            f"carry different elements:\n{shown}\n"
            "Atom numbering is compared within the molecule, so check the "
            "solute atom indices and the atom order of the two runs."
        )

    # -- matching ---------------------------------------------------------- #
    def _index_by_identity(
        self, rows: Sequence[ParameterRow], kind: ParameterKind
    ) -> Dict[Tuple[int, ...], ParameterRow]:
        indexed: Dict[Tuple[int, ...], ParameterRow] = {}
        for row in rows:
            indexed.setdefault(canonical_identity(kind, row.local_atoms), row)
        return indexed

    def _compare(self) -> None:
        kind = self.parsed_1.kind

        # Filter first: a cutoff artifact that survives in one file only would
        # otherwise be reported as a parameter the other simulation is missing.
        rows_1 = [row for row in self.parsed_1.rows if row.occurrence >= self.min_occurrence]
        rows_2 = [row for row in self.parsed_2.rows if row.occurrence >= self.min_occurrence]
        self.filtered_out_1 = len(self.parsed_1.rows) - len(rows_1)
        self.filtered_out_2 = len(self.parsed_2.rows) - len(rows_2)

        indexed_1 = self._index_by_identity(rows_1, kind)
        indexed_2 = self._index_by_identity(rows_2, kind)

        self.matched = []
        for identity in sorted(set(indexed_1) & set(indexed_2)):
            row_1 = indexed_1[identity]
            row_2 = indexed_2[identity]
            if kind.periodic:
                difference = angular_difference(row_1.average, row_2.average)
            else:
                difference = row_2.average - row_1.average
            self.matched.append(
                MatchedParameter(
                    local_atoms=identity,
                    elements=canonical_elements(kind, row_1.local_atoms, row_1.elements),
                    average_1=row_1.average,
                    std_1=row_1.std_dev,
                    occurrence_1=row_1.occurrence,
                    source_row_1=row_1.source_row,
                    average_2=row_2.average,
                    std_2=row_2.std_dev,
                    occurrence_2=row_2.occurrence,
                    source_row_2=row_2.source_row,
                    difference=difference,
                )
            )

        self.only_in_file_1 = [
            indexed_1[identity] for identity in sorted(set(indexed_1) - set(indexed_2))
        ]
        self.only_in_file_2 = [
            indexed_2[identity] for identity in sorted(set(indexed_2) - set(indexed_1))
        ]

        if not self.matched:
            raise ValueError(
                f"No {kind.label} parameter is present in both files, so there is "
                "nothing to compare. The two runs may describe different molecules, "
                "or the minimum occurrence fraction may be too high."
            )

    # -- output ------------------------------------------------------------ #
    @property
    def column_labels(self) -> Tuple[str, str]:
        """The tokens that name each file's columns, kept distinct."""
        if self.label_1 == self.label_2:
            return f"{self.label_1}_1", f"{self.label_2}_2"
        return self.label_1, self.label_2

    def matched_columns(self) -> List[str]:
        kind = self.parsed_1.kind
        first, second = self.column_labels
        return (
            ["row"]
            + list(kind.atom_columns)
            + list(kind.element_columns)
            + [f"average_{first}", f"std_dev_{first}", f"occurrence_{first}"]
            + [f"average_{second}", f"std_dev_{second}", f"occurrence_{second}"]
            + [f"difference_{second}_minus_{first}", "source_row_1", "source_row_2"]
        )

    def unmatched_columns(self) -> List[str]:
        kind = self.parsed_1.kind
        return (
            ["row"]
            + list(kind.atom_columns)
            + list(kind.element_columns)
            + ["average", "std_dev", "occurrence", "source_row"]
        )

    def write_results(self, output_file: str) -> str:
        """Write the comparison as one self-describing, three-section text file."""
        if self.parsed_1 is None or self.parsed_2 is None or not self.matched:
            raise ValueError(
                "There is no comparison to write. Call run() before write_results()."
            )

        output_file = os.path.abspath(output_file)
        os.makedirs(os.path.dirname(output_file) or os.getcwd(), exist_ok=True)

        kind = self.parsed_1.kind
        with open(output_file, "w") as out:
            out.write("# gQTEA Molecular Geometry Comparison\n")
            out.write(f"# parameter_type {kind.key}\n")
            out.write(f"# value_unit {kind.unit}\n")
            self._write_source_header(out, 1, self.parsed_1, self.label_1, self.filtered_out_1)
            self._write_source_header(out, 2, self.parsed_2, self.label_2, self.filtered_out_2)
            out.write(f"# minimum_occurrence_fraction {self.min_occurrence:g}\n")
            out.write(f"# matched_parameters {len(self.matched)}\n")
            out.write(f"# only_in_file_1 {len(self.only_in_file_1)}\n")
            out.write(f"# only_in_file_2 {len(self.only_in_file_2)}\n")

            matched_columns = self.matched_columns()
            out.write("#\n# section matched\n")
            out.write("# " + " ".join(matched_columns) + "\n")
            for row_index, match in enumerate(self.matched, start=1):
                out.write(self._format_matched_row(row_index, match))

            unmatched_columns = self.unmatched_columns()
            for section, rows in (
                ("only_in_file_1", self.only_in_file_1),
                ("only_in_file_2", self.only_in_file_2),
            ):
                out.write(f"#\n# section {section}\n")
                out.write("# " + " ".join(unmatched_columns) + "\n")
                for row_index, row in enumerate(rows, start=1):
                    out.write(self._format_unmatched_row(row_index, row, kind))

        return output_file

    def _write_source_header(
        self,
        out,
        position: int,
        parsed: GeometryParameterFile,
        label: str,
        filtered_out: int,
    ) -> None:
        frames = parsed.frames_used
        out.write(f"# file_{position} {parsed.path}\n")
        out.write(f"# file_{position}_label {label}\n")
        out.write(f"# file_{position}_frames {frames if frames is not None else 'unknown'}\n")
        out.write(f"# file_{position}_atom_scope {parsed.atom_scope}\n")
        remapped = parsed.remapped_atoms()
        if remapped:
            # Only the atoms whose numbering differs, so a row can be traced
            # back to the source file without widening every row.
            mapping = " ".join(f"{local}:{global_index}" for local, global_index in remapped.items())
            out.write(f"# file_{position}_local_to_global {mapping}\n")
        out.write(f"# file_{position}_rows_below_threshold {filtered_out}\n")

    @staticmethod
    def _format_matched_row(row_index: int, match: MatchedParameter) -> str:
        fields = [f"{row_index:>6d}"]
        fields += [f"{atom:>8d}" for atom in match.local_atoms]
        fields += [f"{element:>8s}" for element in match.elements]
        fields += [
            f"{match.average_1:>16.8f}",
            f"{match.std_1:>16.8f}",
            f"{match.occurrence_1:>16.8f}",
            f"{match.average_2:>16.8f}",
            f"{match.std_2:>16.8f}",
            f"{match.occurrence_2:>16.8f}",
            f"{match.difference:>16.8f}",
            f"{match.source_row_1:>12d}",
            f"{match.source_row_2:>12d}",
        ]
        return " ".join(fields) + "\n"

    @staticmethod
    def _format_unmatched_row(row_index: int, row: ParameterRow, kind: ParameterKind) -> str:
        atoms = canonical_identity(kind, row.local_atoms)
        elements = canonical_elements(kind, row.local_atoms, row.elements)
        fields = [f"{row_index:>6d}"]
        fields += [f"{atom:>8d}" for atom in atoms]
        fields += [f"{element:>8s}" for element in elements]
        fields += [
            f"{row.average:>16.8f}",
            f"{row.std_dev:>16.8f}",
            f"{row.occurrence:>16.8f}",
            f"{row.source_row:>12d}",
        ]
        return " ".join(fields) + "\n"

    # -- reporting --------------------------------------------------------- #
    def summary_text(self, output_file: str, shifts: int = 10) -> str:
        """The closing report shown in the tool's text box."""
        kind = self.parsed_1.kind
        first, second = self.column_labels
        unit = kind.unit_symbol

        lines = [
            "Comparison completed.",
            f"Parameter type: {kind.label}",
            "",
            f"File 1 ({first}): {self.parsed_1.path}",
            f"  frames {self.parsed_1.frames_used}, scope {self.parsed_1.atom_scope}, "
            f"{len(self.parsed_1.rows)} parameters",
            f"File 2 ({second}): {self.parsed_2.path}",
            f"  frames {self.parsed_2.frames_used}, scope {self.parsed_2.atom_scope}, "
            f"{len(self.parsed_2.rows)} parameters",
            "",
            f"Molecule size: {self.parsed_1.num_local_atoms} atoms "
            "(elements agree in both files)",
            f"Minimum occurrence fraction: {self.min_occurrence:g} "
            f"(dropped {self.filtered_out_1} and {self.filtered_out_2} rows)",
            f"Matched parameters: {len(self.matched)}",
            f"Only in file 1: {len(self.only_in_file_1)}",
            f"Only in file 2: {len(self.only_in_file_2)}",
            "",
            f"Output file:\n{output_file}",
            "",
            f"Largest shifts ({second} - {first}):",
        ]

        for match in self.largest_shifts(shifts):
            lines.append(
                f"  {match.atom_label:<20s} {match.average_1:>10.5f} -> "
                f"{match.average_2:>10.5f}   {match.difference:>+10.5f} {unit}"
            )

        remapped = self.parsed_2.remapped_atoms()
        if remapped:
            preview = " ".join(f"{l}->{g}" for l, g in list(remapped.items())[:8])
            if len(remapped) > 8:
                preview += " ..."
            lines += [
                "",
                f"Atom numbering remapped in file 2 (molecule -> box): {preview}",
            ]

        return "\n".join(lines)


# --------------------------------------------------------------------------- #
# UI input helpers (pure, so the validation rules can be tested headless)        #
# --------------------------------------------------------------------------- #
def parse_optional_fraction(text: str, default: float = 0.0) -> float:
    """Read an optional fraction field; a blank field means `default`."""
    stripped = (text or "").strip()
    if not stripped:
        return default

    try:
        value = float(stripped)
    except ValueError:
        raise ValueError(
            f"'{stripped}' is not a number. Enter a minimum occurrence fraction "
            "between 0 and 1, or leave the field blank to keep every parameter."
        )

    if not 0.0 <= value <= 1.0:
        raise ValueError(
            f"The minimum occurrence fraction must lie between 0 and 1 (got {value:g})."
        )
    return value


def resolve_output_dir(reference_file: str, chosen: str) -> str:
    """Where the results go: the chosen folder, or the first file's folder."""
    stripped = (chosen or "").strip()
    if stripped:
        return os.path.abspath(stripped)
    return os.path.dirname(os.path.abspath(reference_file)) or os.getcwd()


class MolGeomComparatorUI:
    """
    The graphical user interface is provided here.
    """

    def __init__(self, *args) -> None:
        self.file_1 = None
        self.file_2 = None
        self.parsed_1 = None
        self.parsed_2 = None
        self.output_dir = os.getcwd()
        self.layout_main_window(*args)

    async def warning_function(self, title: str, message: str) -> None:
        await self.main_window.dialog(toga.InfoDialog(title, message))

    # -- layout ------------------------------------------------------------ #
    def layout_main_window(self, widget) -> None:
        self.main_window = toga.Window(
            title="Comparison of Molecular Geometric Parameters",
            size=(860, 640),
        )

        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style = Pack(margin=(5, 5), text_align=LEFT, width=260)
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=120)
        row_style = Pack(direction="row", margin=(0, 0, 5, 0))

        main_box = toga.Box(style=Pack(direction="column", margin=20))

        title_row = toga.Box(style=Pack(direction="row", margin=(0, 0, 10, 0)))
        title_box = toga.Box(style=Pack(width=800))
        title_box.add(toga.Label("Comparison of Molecular Geometric Parameters", style=heading_style))
        title_row.add(title_box)
        main_box.add(title_row)

        def form_row(label_text, widget_to_add, browse=None):
            row = toga.Box(style=row_style)
            row.add(toga.Label(label_text, style=label_style))
            row.add(widget_to_add)
            if browse is not None:
                row.add(browse)
            main_box.add(row)

        self.textInput_file_1 = toga.TextInput(
            placeholder="Click Browse to select the first parameter file (e.g. the isolated molecule)",
            style=input_style,
        )
        form_row(
            "Parameter file 1 (reference):",
            self.textInput_file_1,
            toga.Button("Browse", on_press=self.open_file_1_dialog, style=button_style),
        )

        self.textInput_file_2 = toga.TextInput(
            placeholder="Click Browse to select the second parameter file (e.g. the solvated molecule)",
            style=input_style,
        )
        form_row(
            "Parameter file 2 (comparison):",
            self.textInput_file_2,
            toga.Button("Browse", on_press=self.open_file_2_dialog, style=button_style),
        )

        self.textInput_label_1 = toga.TextInput(
            value="isolated", placeholder="isolated", style=input_style
        )
        form_row("Label for file 1:", self.textInput_label_1)

        self.textInput_label_2 = toga.TextInput(
            value="solvated", placeholder="solvated", style=input_style
        )
        form_row("Label for file 2:", self.textInput_label_2)

        self.textInput_min_occurrence = toga.TextInput(
            placeholder="Default: 0 (keep every parameter); e.g. 0.05 drops rare contacts",
            style=input_style,
        )
        form_row("Minimum occurrence fraction:", self.textInput_min_occurrence)

        self.textInput_output_dir = toga.TextInput(
            placeholder="Leave blank to write next to parameter file 1",
            style=input_style,
        )
        form_row(
            "Output folder:",
            self.textInput_output_dir,
            toga.Button("Browse", on_press=self.choose_output_dir, style=button_style),
        )

        self.textInput_output = toga.TextInput(
            placeholder="Leave blank for geometry_comparison_<type>.txt",
            style=input_style,
        )
        form_row("Output txt filename:", self.textInput_output)

        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(10, 0), font_size=12)
        )
        self.multi_line_text.value = (
            "This module compares the geometric parameters of the same molecule taken from "
            "two simulations -- typically the isolated molecule and the same molecule in a "
            "water box whose solvent has already been removed.\n\n"
            "Select two combined parameter files written by the All bond distance, All bond "
            "angle or All dihedral angle analysis tools. Both files must be of the same type.\n\n"
            "Atoms are matched within the molecule, not by their raw index: when the solvated "
            "run was analysed with solute atom indices, the tool uses that header to renumber "
            "its atoms back onto the isolated molecule, and it checks the element of every atom "
            "before comparing anything."
        )
        main_box.add(self.multi_line_text)

        button_row = toga.Box(style=Pack(direction="row", margin=(10, 0, 0, 0)))
        self.btn_execute = toga.Button("Compare", style=button_style, on_press=self.workflow)
        self.btn_help = toga.Button("Help", style=button_style, on_press=self.open_window_help)
        self.btn_close = toga.Button("Close", style=button_style, on_press=self.closeTopLevel)
        for button in (self.btn_execute, self.btn_help, self.btn_close):
            button_row.add(button)
        main_box.add(button_row)

        self.main_window.content = main_box
        self.main_window.show()

    # -- file selection ---------------------------------------------------- #
    async def open_file_1_dialog(self, widget) -> None:
        await self._select_input_file(1)

    async def open_file_2_dialog(self, widget) -> None:
        await self._select_input_file(2)

    async def _select_input_file(self, position: int) -> None:
        """Pick a parameter file and read it at once, so a bad file is caught here."""
        try:
            selected_file = await self.main_window.dialog(
                toga.OpenFileDialog(f"Open parameter file {position}")
            )
        except ValueError:
            await self.warning_function("Error", "Open file was canceled.")
            return

        if selected_file is None:
            await self.warning_function("Warning", "No file was selected.")
            return

        path = str(selected_file)
        try:
            parsed = read_parameter_file(path)
        except ValueError as exc:
            await self.warning_function("Error", str(exc))
            return

        if position == 1:
            self.file_1 = path
            self.parsed_1 = parsed
            self.textInput_file_1.value = path
            self.output_dir = os.path.dirname(os.path.abspath(path)) or os.getcwd()
        else:
            self.file_2 = path
            self.parsed_2 = parsed
            self.textInput_file_2.value = path

        self._report_selection()

    def _report_selection(self) -> None:
        """Describe what has been loaded so far in the text box."""
        lines = []
        for position, parsed in ((1, self.parsed_1), (2, self.parsed_2)):
            if parsed is None:
                lines.append(f"Parameter file {position}: not selected yet.")
                continue
            remapped = parsed.remapped_atoms()
            lines.append(
                f"Parameter file {position}: {parsed.name}\n"
                f"  type {parsed.kind.label}, {len(parsed.rows)} parameters, "
                f"{parsed.frames_used} frames\n"
                f"  atom scope {parsed.atom_scope}, molecule size {parsed.num_local_atoms} atoms"
                + (f", {len(remapped)} atoms renumbered from the box" if remapped else "")
            )

        if self.parsed_1 is not None and self.parsed_2 is not None:
            if self.parsed_1.kind is self.parsed_2.kind:
                lines.append("\nBoth files are of the same type. Press Compare.")
            else:
                lines.append(
                    f"\nThe two files hold different parameter types "
                    f"({self.parsed_1.kind.label} and {self.parsed_2.kind.label}). "
                    "Select two files produced by the same analysis tool."
                )

        self.multi_line_text.value = "\n".join(lines)

    async def choose_output_dir(self, widget) -> None:
        try:
            selected_dir = await self.main_window.dialog(
                toga.SelectFolderDialog("Select the output folder")
            )
        except ValueError:
            await self.warning_function("Error", "Folder selection was canceled.")
            return

        if selected_dir is None:
            await self.warning_function("Warning", "No folder was selected.")
            return

        self.textInput_output_dir.value = str(selected_dir)

    # -- parameters -------------------------------------------------------- #
    async def read_params(self) -> bool:
        """Validate the form. Returns False (after a dialog) if anything is wrong."""
        # The path fields are editable, so a typed or pasted path counts as a
        # selection just as Browse does; the field wins when it disagrees.
        for position, field in ((1, self.textInput_file_1), (2, self.textInput_file_2)):
            typed = (field.value or "").strip()
            if typed and typed != getattr(self, f"file_{position}"):
                setattr(self, f"file_{position}", typed)
                # Whatever was read at selection time describes another file now.
                setattr(self, f"parsed_{position}", None)

        if not self.file_1 or not self.file_2:
            await self.warning_function(
                "Error", "Select both parameter files before comparing."
            )
            return False

        for position, path in ((1, self.file_1), (2, self.file_2)):
            if not os.path.isfile(path):
                await self.warning_function(
                    "Error", f"Parameter file {position} does not exist:\n{path}"
                )
                return False

        if os.path.abspath(self.file_1) == os.path.abspath(self.file_2):
            await self.warning_function(
                "Error",
                "The same file was selected twice. Choose two different parameter files.",
            )
            return False

        try:
            self.min_occurrence = parse_optional_fraction(self.textInput_min_occurrence.value)
        except ValueError as exc:
            await self.warning_function("Error", str(exc))
            return False

        self.label_1 = self.textInput_label_1.value.strip() or "isolated"
        self.label_2 = self.textInput_label_2.value.strip() or "solvated"

        self.output_dir = resolve_output_dir(self.file_1, self.textInput_output_dir.value)

        output_name = self.textInput_output.value.strip()
        if not output_name:
            # Before the first run the type is known only if file 1 parsed at
            # selection time; otherwise the workflow renames it once it knows.
            kind = self.parsed_1.kind if self.parsed_1 is not None else BOND_KIND
            output_name = default_output_name(kind)

        if os.path.isabs(output_name):
            self.output_file = output_name
        else:
            self.output_file = os.path.join(self.output_dir, output_name)

        return True

    # -- workflow ---------------------------------------------------------- #
    async def workflow(self, widget) -> None:
        if not await self.read_params():
            return

        calculator = MolGeomCalculator(
            self.file_1,
            self.file_2,
            label_1=self.label_1,
            label_2=self.label_2,
            min_occurrence=self.min_occurrence,
        )

        self.multi_line_text.value = (
            "Reading both parameter files, checking that they describe the same "
            "molecule, and matching their parameters..."
        )

        try:
            calculator.run()
        except ValueError as exc:
            self.multi_line_text.value = f"The comparison was not carried out.\n\n{exc}"
            await self.warning_function("Error", str(exc))
            return

        self.parsed_1 = calculator.parsed_1
        self.parsed_2 = calculator.parsed_2

        # The auto-named output follows the type, which is certain only now.
        if not self.textInput_output.value.strip():
            self.output_file = os.path.join(
                self.output_dir, default_output_name(calculator.kind)
            )

        try:
            output_file = calculator.write_results(self.output_file)
        except (ValueError, OSError) as exc:
            self.multi_line_text.value = f"The results could not be written.\n\n{exc}"
            await self.warning_function("Error", f"Could not write the output file: {exc}")
            return

        self.multi_line_text.value = calculator.summary_text(output_file)

    def open_window_help(self, widget) -> None:
        window = toga.Window(title="Instructions to compare molecular geometric parameters")
        box = toga.Box(style=Pack(direction="column", flex=1))
        text = toga.MultilineTextInput(style=Pack(font_size=11, margin=(5, 5), flex=1))
        text.value = HelpGqteaWin.help_mol_geom_comparator
        box.add(text)
        window.content = box
        window.show()

    def closeTopLevel(self, widget) -> None:
        self.main_window.close()
