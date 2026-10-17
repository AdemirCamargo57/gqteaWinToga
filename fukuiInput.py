"""
Gaussian 16 input for a finite-difference Fukui function calculation.

One ``<stem>_fukui.gjf`` holds three linked jobs (``--Link1--``) for Gaussian
16 (Rev. C.01) under Linux -- WSL or a cluster:

  1. N state: geometry optimization from the given Cartesian geometry,
     checkpoint ``<stem>_N.chk``;
  2. N+p state: single point, ``%OldChk=<stem>_N.chk`` + ``Geom=Check``
     (the optimized N geometry), its own ``<stem>_Np1.chk``;
  3. N-q state: the same, ``<stem>_Nm1.chk``.

``%OldChk`` is copied into the job's own ``%Chk`` before the job starts, so
the N checkpoint is never overwritten. ``<stem>_formchk.sh`` turns the three
checkpoints into ``<stem>_N.fchk`` ... for Tools > Analysis of Fukui
Functions. File names are relative and the files use LF line endings: the
jobs are run from the folder that holds the input.

No Toga here (``fukuiInputBuilder`` is the window).
"""

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from help import AtomicData
from systemProperties import allowed_multiplicities, multiplicity_problem, nuclear_charge

PRINT_LEVELS = ("#P", "#", "#N", "#T")
DEFAULT_PRINT_LEVEL = "#P"
DEFAULT_METHOD = "B3LYP"
DEFAULT_BASIS = "6-31G(d)"          # the level of the phenol test fixtures
DEFAULT_KEYWORDS = "SCF=Tight NoSymm"
DEFAULT_OPT_KEYWORD = "Opt"
DEFAULT_MEMORY = "2GB"
DEFAULT_NPROC = 4
DEFAULT_STEM = "molecule"

INPUT_SUFFIX = "_fukui.gjf"
SCRIPT_SUFFIX = "_formchk.sh"
FCHK_EXTENSION = ".fchk"
ROUTE_WIDTH = 72                    # route lines are wrapped below Gaussian's 80 columns
MIN_ATOM_DISTANCE = 0.4             # Angstrom; closer atoms are a typing error, not a molecule
TITLE_FORBIDDEN = "@#!-_\\"         # characters Gaussian does not allow in the title section

MEMORY_UNITS = ("KB", "MB", "GB", "TB", "KW", "MW", "GW", "TW")
_MEMORY_RE = re.compile(r"^(\d+)\s*([KMGT][BW])$", re.IGNORECASE)
_STEM_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_LABEL_RE = re.compile(r"^([A-Za-z]{1,2})\d*$")

# Route keywords that cannot be in the shared route (the builder writes Opt, Geom and Guess itself).
FORBIDDEN_KEYWORDS = {
    "OPT": "Optimization keywords go in the 'Optimization (job 1)' field: only the N state is optimized; "
           "the ions are single points on its geometry.",
    "IRC": "IRC is not part of a Fukui calculation.",
    "IRCMAX": "IRCMax is not part of a Fukui calculation.",
    "SCAN": "Scan is not part of a Fukui calculation.",
    "GEOM": "Geom= is managed by the builder (Geom=Check in the ion jobs).",
    "GUESS": "Guess= is managed by the builder (see 'Start the ions from the N-state orbitals').",
    "UNITS": "Units= is not allowed: the coordinates are written in Angstrom.",
    "PSEUDO": "Pseudo= (ECP) needs extra input, and the Hirshfeld indices do not support ECPs.",
    "EXTRABASIS": "ExtraBasis needs an extra input section, which the builder does not write.",
    "CHKBASIS": "ChkBasis is not needed: every job has the basis set in its route.",
}
# Opt options that need extra input or do not lead to a minimum.
FORBIDDEN_OPT_OPTIONS = ("TS", "QST2", "QST3", "SADDLE", "CONICAL", "MODREDUNDANT", "ADDREDUNDANT",
                         "ADDGIC", "GIC", "READFC", "RCFC", "RDOPT")
GENERAL_BASIS = ("GEN", "GENECP")
POST_SCF_METHODS = ("MP2", "MP3", "MP4", "MP5", "CCSD", "CCD", "QCISD", "CISD", "CID", "CIS", "BD",
                    "CASSCF", "SAC-CI", "EOM")
POST_SCF_KEYWORDS = ("TD", "TDA", "CIS", "EOMCCSD", "SAC-CI")

Atom = Tuple[str, float, float, float]

_SYMBOLS = {int(z): s for s, z in AtomicData.atomic_numbers.items()}
_SYMBOLS_UPPER = {s.upper(): s for s in AtomicData.atomic_numbers}


# ----------------------------------------------------------------------
# Field parsers
# ----------------------------------------------------------------------
def parse_memory(text: Optional[str]) -> str:
    """'2GB', '500mb', '16 GB' -> '2GB', '500MB', '16GB'. A bare number is refused."""
    value = (text or "").strip()
    if re.fullmatch(r"\d+", value):
        raise ValueError(f"%Mem = '{value}': a number without a unit means 8-byte words to Gaussian. "
                         "Add a unit, e.g. 2GB or 500MB.")
    match = _MEMORY_RE.match(value)
    if not match or int(match.group(1)) <= 0:
        raise ValueError(f"%Mem must be a positive integer followed by a unit "
                         f"({', '.join(MEMORY_UNITS)}), e.g. 2GB (got '{value}').")
    return f"{int(match.group(1))}{match.group(2).upper()}"


def parse_nproc(text: Optional[str]) -> int:
    value = (text or "").strip()
    if not re.fullmatch(r"\d+", value) or int(value) < 1:
        raise ValueError(f"%NProcShared must be a positive integer (got '{value}').")
    return int(value)


def validate_stem(text: Optional[str]) -> str:
    value = (text or "").strip()
    if not _STEM_RE.match(value) or not value.isascii():
        raise ValueError(f"Job name '{value}': use letters, digits, '.', '_' or '-' only, starting with a "
                         "letter or digit (it becomes the .gjf and .chk file names).")
    return value


def parse_opt_keyword(text: Optional[str]) -> str:
    """The job-1 optimization keyword: 'Opt' or 'Opt=...' / 'Opt(...)' towards a minimum."""
    value = (text or "").strip()
    tokens = route_tokens(value)
    if len(tokens) != 1 or keyword_name(tokens[0]) != "OPT":
        raise ValueError(f"Optimization (job 1) must be one Opt keyword, e.g. Opt or Opt=(Tight,CalcFC) "
                         f"(got '{value}').")
    options = re.split(r"[=(),\s]+", tokens[0].upper())[1:]
    bad = [o for o in options if o in FORBIDDEN_OPT_OPTIONS or "READ" in o]
    if bad:
        raise ValueError(f"Optimization (job 1): option {bad[0]} is not supported (it needs extra input "
                         "or does not optimize to a minimum).")
    return value


def element_symbol(token: str) -> Optional[str]:
    if token.isdigit():
        return _SYMBOLS.get(int(token))
    match = _LABEL_RE.match(token)
    return _SYMBOLS_UPPER.get(match.group(1).upper()) if match else None


def parse_geometry(text: Optional[str], check_distances: bool = True) -> List[Atom]:
    """Cartesian 'Element x y z' lines in Angstrom (symbol, numbered label or atomic number).

    A full .xyz file (count line + comment line) is accepted as pasted.
    ``check_distances=False`` skips the O(n^2) overlap test (live counts while typing).
    """
    lines = [(k, l.strip()) for k, l in enumerate((text or "").splitlines(), start=1)]
    lines = [(k, l) for k, l in lines if l]
    if len(lines) >= 2 and re.fullmatch(r"\d+", lines[0][1]) and int(lines[0][1]) == len(lines) - 2:
        lines = lines[2:]
    if not lines:
        raise ValueError("The molecular geometry is empty: paste Cartesian 'Element x y z' lines "
                         "(Angstrom) or use 'Load geometry...'.")
    atoms: List[Atom] = []
    for number, line in lines:
        tokens = line.split()
        if len(tokens) != 4:
            raise ValueError(f"Geometry line {number} ('{line}'): expected 'Element x y z' in Angstrom. "
                             "Only Cartesian coordinates are supported (no Z-matrix, freeze flags or "
                             "fragment tags).")
        symbol = element_symbol(tokens[0])
        if symbol is None:
            raise ValueError(f"Geometry line {number}: unknown element '{tokens[0]}'.")
        try:
            xyz = [float(t) for t in tokens[1:]]
        except ValueError:
            raise ValueError(f"Geometry line {number} ('{line}'): the coordinates are not numbers.")
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError(f"Geometry line {number} ('{line}'): the coordinates must be finite.")
        atoms.append((symbol, *xyz))
    for i in range(len(atoms) if check_distances else 0):
        for j in range(i + 1, len(atoms)):
            d = math.dist(atoms[i][1:], atoms[j][1:])
            if d < MIN_ATOM_DISTANCE:
                raise ValueError(f"Geometry: atoms {i + 1} and {j + 1} are only {d:.3f} A apart.")
    return atoms


def format_geometry(atoms: Sequence[Atom]) -> List[str]:
    """Fixed-width Angstrom lines, the layout Gaussian itself writes."""
    return [f"{s:<2s}{x:16.8f}{y:16.8f}{z:16.8f}" for s, x, y, z in atoms]


# ----------------------------------------------------------------------
# Route section
# ----------------------------------------------------------------------
def route_tokens(text: Optional[str]) -> List[str]:
    """Split keywords on blanks outside parentheses: 'SCF=(Tight, XQC) NoSymm' -> 2 tokens."""
    tokens, current, depth = [], "", 0
    for ch in (text or "").strip():
        depth += (ch == "(") - (ch == ")")
        if ch.isspace() and depth == 0:
            if current:
                tokens.append(current)
            current = ""
        else:
            current += ch
    if current:
        tokens.append(current)
    return tokens


def keyword_name(token: str) -> str:
    return re.split(r"[=(]", token, maxsplit=1)[0].upper()


def strip_prefix(method: str) -> Tuple[str, str]:
    """('R'|'U'|'RO'|'', rest). 'REV...' functionals are not an R prefix."""
    m = method.upper()
    if m.startswith("REV"):
        return "", m
    for prefix in ("RO", "R", "U"):
        if m.startswith(prefix) and len(m) > len(prefix):
            return prefix, m[len(prefix):]
    return "", m


def check_route(method: str, basis: str, keywords: str) -> Tuple[List[str], List[str]]:
    """(errors, warnings) of the method / basis / keyword fields shared by the three jobs."""
    errors: List[str] = []
    warnings: List[str] = []
    method, basis = (method or "").strip(), (basis or "").strip()
    if not method or len(method.split()) != 1:
        errors.append("Method must be one keyword, e.g. B3LYP.")
    elif "/" in method:
        errors.append("Method: enter the basis set in the Basis set field, not after a '/'.")
    if len(basis.split()) > 1:
        errors.append("Basis set must be one keyword, e.g. 6-31G(d).")
    if basis.upper() in GENERAL_BASIS:
        errors.append(f"Basis set {basis}: Gen/GenECP need a basis-set section that the builder does not write.")
    for token in route_tokens(keywords):
        name = keyword_name(token)
        if name in FORBIDDEN_KEYWORDS:
            errors.append(f"Keyword {token}: {FORBIDDEN_KEYWORDS[name]}")
        elif "READ" in token.upper():
            errors.append(f"Keyword {token}: options that read extra input sections are not supported.")
        elif name == "FREQ":
            warnings.append("Freq in the shared keywords runs a frequency job on all three states; it is not "
                            "needed for Fukui functions.")
        elif name in POST_SCF_KEYWORDS:
            warnings.append(f"Keyword {token}: the Fukui analysis reads only the Total SCF Density, "
                            "not excited-state or correlated densities.")
    _, bare = strip_prefix(method)
    if any(bare.startswith(m) for m in POST_SCF_METHODS):
        warnings.append(f"Method {method}: the Fukui analysis reads only the Total SCF Density "
                        "(the HF density here), not the correlated density.")
    b = basis.upper()
    if "+" in b or b.startswith("AUG-") or "D-AUG" in b:
        warnings.append(f"Basis set {basis} has diffuse functions: good for the electron affinity, but "
                        "with self-interaction-prone functionals the N+1 density can become over-diffuse "
                        "(the analysis tool repeats this note).")
    if not basis and not errors:
        warnings.append("No basis set given: the method keyword alone is written.")
    return errors, warnings


def build_route(print_level: str, method: str, basis: str, keywords: str, extra: Sequence[str] = ()) -> List[str]:
    """The route section, wrapped at ROUTE_WIDTH; ``extra`` (Opt, Geom=Check, ...) goes last."""
    model = f"{method.strip()}/{basis.strip()}" if basis.strip() else method.strip()
    words = [print_level, model, *route_tokens(keywords), *extra]
    lines, current = [], ""
    for word in words:
        if current and len(current) + 1 + len(word) > ROUTE_WIDTH:
            lines.append(current)
            current = word
        else:
            current = f"{current} {word}" if current else word
    lines.append(current)
    return lines


# ----------------------------------------------------------------------
# The job
# ----------------------------------------------------------------------
@dataclass
class FukuiJobSpec:
    """Everything the window collects; ``prepare`` validates it.

    A multiplicity of None means "not given": ``prepare`` asks for it rather than assuming one.
    """
    stem: str
    output_dir: str
    geometry: str
    charge: int
    multiplicity: Optional[int]
    plus_multiplicity: Optional[int]
    minus_multiplicity: Optional[int]
    memory: str = DEFAULT_MEMORY
    nproc: int = DEFAULT_NPROC
    print_level: str = DEFAULT_PRINT_LEVEL
    method: str = DEFAULT_METHOD
    basis: str = DEFAULT_BASIS
    keywords: str = DEFAULT_KEYWORDS
    opt_keyword: str = DEFAULT_OPT_KEYWORD
    p: int = 1
    q: int = 1
    guess_read: bool = True


@dataclass
class FukuiState:
    label: str            # "N", "N+1", "N-1"
    tag: str              # "N", "Np1", "Nm1" (file names)
    words: str            # "N", "N plus 1", "N minus 1" (title line)
    charge: int
    multiplicity: int
    n_electrons: int
    chk: str              # as written after %Chk= (relative to the input's folder)
    fchk: Path            # where <stem>_formchk.sh puts the formatted checkpoint


@dataclass
class FukuiJob:
    spec: FukuiJobSpec
    atoms: List[Atom]
    states: List[FukuiState]
    memory: str
    opt_keyword: str
    input_path: Path
    script_path: Path
    warnings: List[str] = field(default_factory=list)


def state_tags(p: int, q: int):
    """(label, tag, words, k, sign) of the N, N+p, N-q states."""
    return (("N", "N", "N", 0, 0),
            (f"N+{p}", f"Np{p}", f"N plus {p}", p, -1),
            (f"N-{q}", f"Nm{q}", f"N minus {q}", q, +1))


def _title_text(stem: str) -> str:
    cleaned = "".join(" " if ch in TITLE_FORBIDDEN else ch for ch in stem)
    return " ".join(cleaned.split()) or DEFAULT_STEM


def _allowed_text(n_electrons: int) -> str:
    return ", ".join(str(m) for m in allowed_multiplicities(n_electrons, 3)) + ", ..."


class FukuiInputBuilder:
    """Validation, the linked-job text, the formchk script and the files."""

    def prepare(self, spec: FukuiJobSpec) -> FukuiJob:
        """Validate every field; ValueError lists all the problems found."""
        errors: List[str] = []
        warnings: List[str] = []

        def attempt(fn, *args):
            try:
                return fn(*args)
            except ValueError as exc:
                errors.append(str(exc))
                return None

        stem = attempt(validate_stem, spec.stem)
        folder_text = (spec.output_dir or "").strip()
        folder = Path(folder_text).resolve() if folder_text else None
        if folder is None or not folder.is_dir():
            errors.append(f"Output folder '{folder_text}' does not exist; choose an existing folder."
                          if folder_text else "Output folder: choose where the input files are written.")
        memory = attempt(parse_memory, spec.memory)
        if not isinstance(spec.nproc, int) or spec.nproc < 1:
            errors.append(f"%NProcShared must be a positive integer (got {spec.nproc}).")
        if spec.print_level not in PRINT_LEVELS:
            errors.append(f"Print level must be one of {', '.join(PRINT_LEVELS)} (got '{spec.print_level}').")
        opt_keyword = attempt(parse_opt_keyword, spec.opt_keyword)
        route_errors, route_warnings = check_route(spec.method, spec.basis, spec.keywords)
        errors += route_errors
        warnings += route_warnings
        for name, k in (("p", spec.p), ("q", spec.q)):
            if not isinstance(k, int) or k < 1:
                errors.append(f"{name} (electrons added/removed) must be a positive integer (got {k}).")
        atoms = attempt(parse_geometry, spec.geometry)

        states: List[FukuiState] = []
        if atoms is not None and not errors:
            z_total = nuclear_charge(a[0] for a in atoms)
            prefix, _ = strip_prefix(spec.method)
            multiplicities = (spec.multiplicity, spec.plus_multiplicity, spec.minus_multiplicity)
            for (label, tag, words, k, sign), mult in zip(state_tags(spec.p, spec.q), multiplicities):
                charge = spec.charge + sign * k
                n_electrons = z_total - charge
                if n_electrons < 1:
                    errors.append(f"The {label} state (charge {charge}) would have no electrons left.")
                    continue
                if mult is None:
                    errors.append(f"The {label} state: enter its spin multiplicity; it is not assumed "
                                  f"({n_electrons} electrons allow {_allowed_text(n_electrons)}).")
                    continue
                problem = multiplicity_problem(n_electrons, mult)
                if problem:
                    errors.append(f"The {label} state: multiplicity {mult} is impossible ({problem}).")
                if k and spec.multiplicity is not None and abs(mult - spec.multiplicity) != k:
                    warnings.append(f"The {label} state has multiplicity {mult} (N state: {spec.multiplicity}); "
                                    f"adding/removing {k} electron(s) usually changes it by {k}. Confirm this "
                                    "is the lowest-energy spin state.")
                if prefix == "R" and mult > 1:
                    warnings.append(f"Method {spec.method} forces a closed-shell wavefunction, but the {label} "
                                    f"state has multiplicity {mult}; remove the R so Gaussian chooses R or U.")
                states.append(FukuiState(label, tag, words, charge, mult, n_electrons, f"{stem}_{tag}.chk",
                                         folder / f"{stem}_{tag}{FCHK_EXTENSION}"))
        if errors:
            raise ValueError("\n".join(f"- {e}" for e in errors))
        return FukuiJob(spec, atoms, states, memory, opt_keyword, folder / f"{stem}{INPUT_SUFFIX}",
                        folder / f"{stem}{SCRIPT_SUFFIX}", warnings)

    def input_text(self, job: FukuiJob) -> str:
        """Job 1 optimizes N; jobs 2-3 read its geometry (Geom=Check) into their own checkpoints."""
        spec = job.spec
        title = _title_text(spec.stem)
        reference = job.states[0]
        blocks = []
        for k, state in enumerate(job.states):
            lines = []
            if k == 0:
                extra, what = (job.opt_keyword,), " optimization"
            else:
                extra, what = ("Geom=Check", "Guess=Read") if spec.guess_read else ("Geom=Check",), ""
                lines.append(f"%OldChk={reference.chk}")
            lines += [f"%Chk={state.chk}", f"%Mem={job.memory}", f"%NProcShared={spec.nproc}"]
            lines += build_route(spec.print_level, spec.method, spec.basis, spec.keywords, extra)
            charge = f"minus {-state.charge}" if state.charge < 0 else str(state.charge)  # no '-' in a title
            lines += ["", f"{title} Fukui {state.words} state{what} charge {charge} "
                          f"multiplicity {state.multiplicity}", "",
                      f"{state.charge} {state.multiplicity}"]
            if k == 0:
                lines += format_geometry(job.atoms)
            lines.append("")
            blocks.append("\n".join(lines))
        return "\n--Link1--\n".join(blocks) + "\n"

    def script_text(self, job: FukuiJob) -> str:
        """A bash script that converts the three checkpoints, run where Gaussian ran."""
        lines = ["#!/bin/bash",
                 "# gQTEA Fukui input builder: convert the three checkpoints to formatted",
                 "# checkpoints for Tools > Analysis of Fukui Functions. Run it where Gaussian",
                 f"# ran (WSL or the cluster), after the jobs finish:  bash {job.script_path.name}",
                 'cd "$(dirname "$0")" || exit 1',
                 "if ! command -v formchk >/dev/null 2>&1; then",
                 "  echo \"formchk not found: set up Gaussian 16 first (e.g. 'module load gaussian', or source",
                 "  \\$g16root/g16/bsd/g16.profile).\" >&2",
                 "  exit 1",
                 "fi",
                 "status=0",
                 "convert() {",
                 '  if [ ! -f "$1" ]; then',
                 '    echo "$1 not found: did that job run? Gaussian stops the remaining linked jobs when one '
                 'fails; see the .log file." >&2',
                 "    status=1",
                 "    return",
                 "  fi",
                 '  formchk "$1" "$2" || { echo "formchk failed for $1" >&2; status=1; }',
                 "}"]
        lines += [f'convert "{s.chk}" "{s.fchk.name}"' for s in job.states]
        lines.append("exit $status")
        return "\n".join(lines) + "\n"

    def write_files(self, job: FukuiJob) -> List[Path]:
        """Write the .gjf and the .sh with LF line endings; ValueError on any OS error."""
        written = []
        for path, text in ((job.input_path, self.input_text(job)), (job.script_path, self.script_text(job))):
            try:
                with open(path, "w", encoding="ascii", errors="strict", newline="\n") as handle:
                    handle.write(text)
            except (OSError, UnicodeEncodeError) as exc:
                raise ValueError(f"Could not write {path}: {exc}") from exc
            written.append(path)
        return written


def sibling_state_files(n_path) -> Tuple[Optional[Path], Optional[Path]]:
    """For '<stem>_N.fchk' (or .fch): the existing '<stem>_Np<k>' and '<stem>_Nm<k>' files, else None."""
    n_path = Path(n_path)
    match = re.match(r"^(.*)_N$", n_path.stem)
    if not match or n_path.suffix.lower() not in (".fchk", ".fch"):
        return None, None
    stem = match.group(1)

    def find(letter: str) -> Optional[Path]:
        found = sorted(p for p in n_path.parent.glob(f"{stem}_N{letter}*")
                       if re.fullmatch(rf"{re.escape(stem)}_N{letter}\d+", p.stem)
                       and p.suffix.lower() in (".fchk", ".fch"))
        same_ext = [p for p in found if p.suffix.lower() == n_path.suffix.lower()]
        candidates = same_ext or found
        return candidates[0] if len(candidates) == 1 else None

    return find("p"), find("m")
