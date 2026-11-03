"""
Gaussian 16 input builder: the job model, validation, route, input text and reopening (no Toga).

``GaussianJobSpec`` (job-wide fields + one ``GaussianStepSpec`` per linked step) ->
``GaussianInputBuilder.prepare`` (every problem in one ValueError, one '- Step n: ...' line each) ->
``GaussianJob`` -> ``input_text`` / ``script_text`` / ``write_files``. ``read_input_text`` is the
inverse used by Open input. The Fukui Job Type is delegated to ``fukuiInput`` unchanged. Files are
ASCII with LF endings and relative checkpoint names: the job runs under Linux from the input's folder.
"""

import dataclasses
import math
import re
from dataclasses import dataclass, field
from pathlib import Path, PureWindowsPath
from typing import Dict, List, Optional, Sequence, Tuple

import fukuiInput as fi
import gaussianKeywords as gk
from systemProperties import allowed_multiplicities, multiplicity_problem, nuclear_charge

PARENTHESIS_KEYWORDS = ("IOP",)        # written as IOp(9/40=4), the form every Gaussian manual uses


# ----------------------------------------------------------------------
# Route keywords
# ----------------------------------------------------------------------
@dataclass
class RouteKeyword:
    """One route keyword and its options, e.g. Opt=(TS,CalcFC) -> ('Opt', ['TS', 'CalcFC'])."""
    name: str
    options: List[str] = field(default_factory=list)

    @property
    def key(self) -> str:
        return self.name.upper()


def option_key(option: str) -> str:
    """'MaxCycles=50' -> 'MAXCYCLES': two options with the same key set the same thing."""
    return option.split("=", 1)[0].strip().upper()


def _outer_parentheses(text: str) -> str:
    """'(a,b)' -> 'a,b' when the first '(' closes at the very end; anything else unchanged."""
    if not (text.startswith("(") and text.endswith(")")):
        return text
    depth = 0
    for position, ch in enumerate(text):
        depth += (ch == "(") - (ch == ")")
        if depth == 0 and position < len(text) - 1:
            return text
    return text[1:-1]


def _split_options(text: str) -> List[str]:
    """Split on commas and blanks outside parentheses."""
    options, current, depth = [], "", 0
    for ch in text:
        if depth == 0 and (ch == "," or ch.isspace()):
            if current:
                options.append(current)
            current = ""
            continue
        depth += (ch == "(") - (ch == ")")
        current += ch
    if current:
        options.append(current)
    return options


def keyword_text_name(token: str) -> str:
    """The keyword part of a route token as written: 'Opt=(TS)' -> 'Opt', 'IOp(9/40=4)' -> 'IOp'."""
    return re.split(r"[=(]", token, maxsplit=1)[0]


def parse_keyword(token: str) -> RouteKeyword:
    name = keyword_text_name(token)
    rest = token[len(name):]
    if rest.startswith("="):
        rest = rest[1:]
    return RouteKeyword(name, _split_options(_outer_parentheses(rest)))


def format_keyword(keyword: RouteKeyword) -> str:
    options = keyword.options
    if not options:
        return keyword.name
    if keyword.key in PARENTHESIS_KEYWORDS:
        return f"{keyword.name}({','.join(options)})"
    if len(options) == 1 and "=" not in options[0]:
        return f"{keyword.name}={options[0]}"
    return f"{keyword.name}=({','.join(options)})"


def merge_keywords(panel: Sequence[RouteKeyword], extra_text: str) -> Tuple[List[RouteKeyword], List[str]]:
    """The panels' keywords plus Additional Keywords (spec 5.1).

    A keyword already written by a panel takes the extra options; an option repeated with the same
    key is kept once, and a valued option (X=...) given another value is an error naming both places.
    """
    merged = [RouteKeyword(k.name, list(k.options)) for k in panel]
    panel_keys = {k.key for k in panel}
    errors: List[str] = []
    for token in fi.route_tokens(extra_text):
        if "/" in keyword_text_name(token):
            errors.append(f"Additional Keywords: '{token}' is a method/basis set; choose them in the "
                          "Method panel.")
            continue
        extra = parse_keyword(token)
        target = next((k for k in merged if k.key == extra.key), None)
        if target is None:
            merged.append(extra)
            continue
        for option in extra.options:
            same = [o for o in target.options if option_key(o) == option_key(option)]
            if not same:
                target.options.append(option)
            elif "=" in option and same[0].upper() != option.upper():
                owner = (f"the {gk.PANEL_OF_KEYWORD.get(extra.key, 'Job Type')} panel"
                         if extra.key in panel_keys else "an earlier Additional Keyword")
                errors.append(f"Additional Keywords: {token} sets {option}, but {owner} sets {same[0]}.")
    return merged, errors


# ----------------------------------------------------------------------
# One step's panels
# ----------------------------------------------------------------------
@dataclass
class GaussianStepSpec:
    """Every panel of one linked step, as the window shows it (labels, not keywords).

    A multiplicity of None means "not given": ``prepare`` asks for it rather than assuming one.
    """
    # Job Type
    job_type: str = gk.JOB_ENERGY
    opt_target: str = gk.OPT_MINIMUM
    force_constants: str = gk.FC_DEFAULT
    opt_tight: bool = False
    opt_max_cycles: Optional[int] = None
    freq_raman: bool = False
    freq_vcd: bool = False
    freq_anharmonic: bool = False
    irc_direction: str = "Both"
    irc_force_constants: str = gk.FC_CALCFC
    irc_max_points: Optional[int] = None
    irc_step_size: Optional[int] = None
    stable_opt: bool = False
    nmr_method: str = "GIAO"
    nmr_spinspin: bool = False
    fukui_p: int = 1
    fukui_q: int = 1
    fukui_plus_multiplicity: Optional[int] = None
    fukui_minus_multiplicity: Optional[int] = None
    fukui_guess_read: bool = True
    # Method
    state: str = gk.STATE_GROUND
    n_states: int = gk.DEFAULT_N_STATES
    root: int = 1
    state_spin: str = "Default"
    print_all_coefficients: bool = False
    spin: str = "Default"
    family: str = gk.FAMILY_DFT
    method_name: str = "B3LYP"
    method_other: str = ""
    basis: str = "6-31G"
    basis_diffuse: str = "None"
    basis_polarization: str = "(d)"
    basis_other: str = ""
    dispersion: str = "None"
    # Molecule
    geometry_source: str = gk.GEOM_COORDINATES
    charge: int = 0
    multiplicity: Optional[int] = None
    modredundant: str = ""
    # Title
    title: str = ""
    # Link 0
    memory: str = fi.DEFAULT_MEMORY
    nproc: int = fi.DEFAULT_NPROC
    write_chk: bool = True
    chk: str = ""
    old_chk: str = ""
    extra_link0: str = ""
    # General
    print_level: str = fi.DEFAULT_PRINT_LEVEL
    scf_tight: bool = False
    scf_qc: bool = False
    nosymm: bool = False
    polar: bool = False
    connectivity: bool = False
    # Guess
    guess: str = "Default"
    guess_mix: bool = False
    guess_always: bool = False
    guess_only: bool = False
    guess_save: bool = False
    # Pop
    pop_level: str = "Default"
    pop_mk: bool = False
    pop_hirshfeld: bool = False
    pop_nbo: bool = False
    pop_npa: bool = False
    pop_savento: bool = False
    # Solvation
    solvation: str = gk.SOLVATION_NONE
    solvent: str = "Water"
    solvent_other: str = ""
    # Add. Inp. and the Additional Keywords field
    additional_input: str = ""
    additional_keywords: str = ""


FIELD_LABELS: Dict[str, str] = {
    "job_type": "Job type", "opt_target": "Optimize to a", "force_constants": "Force constants",
    "opt_tight": "Tight convergence", "opt_max_cycles": "Maximum optimization cycles",
    "freq_raman": "Raman", "freq_vcd": "VCD", "freq_anharmonic": "Anharmonic",
    "irc_direction": "IRC direction", "irc_force_constants": "IRC force constants",
    "irc_max_points": "IRC MaxPoints", "irc_step_size": "IRC StepSize", "stable_opt": "Reoptimize",
    "nmr_method": "NMR method", "nmr_spinspin": "Spin-spin couplings",
    "fukui_p": "p (electrons added)", "fukui_q": "q (electrons removed)",
    "fukui_plus_multiplicity": "Multiplicity of N+p", "fukui_minus_multiplicity": "Multiplicity of N-q",
    "fukui_guess_read": "Ions start from the N-state orbitals",
    "state": "Electronic state", "n_states": "Number of states", "root": "State of interest (Root)",
    "state_spin": "States", "print_all_coefficients": "Print all excitation coefficients",
    "spin": "Spin", "family": "Method family", "method_name": "Method", "method_other": "Other method",
    "basis": "Basis set", "basis_diffuse": "Diffuse functions", "basis_polarization": "Polarization",
    "basis_other": "Other basis set", "dispersion": "Empirical dispersion",
    "geometry_source": "Geometry source", "charge": "Charge", "multiplicity": "Multiplicity",
    "modredundant": "Scan/freeze lines", "title": "Title", "memory": "%Mem", "nproc": "%NProcShared",
    "write_chk": "Write a checkpoint (%Chk)", "chk": "%Chk", "old_chk": "%OldChk",
    "extra_link0": "Extra Link 0 lines", "print_level": "Print level", "scf_tight": "SCF=Tight",
    "scf_qc": "SCF=QC", "nosymm": "NoSymm", "polar": "Polar", "connectivity": "Write connectivity",
    "guess": "Initial guess", "guess_mix": "Mix HOMO and LUMO", "guess_always": "Always generate",
    "guess_only": "Only do guess", "guess_save": "Save orbitals", "pop_level": "Population output",
    "pop_mk": "MK charges", "pop_hirshfeld": "Hirshfeld (CM5) charges", "pop_nbo": "NBO", "pop_npa": "NPA",
    "pop_savento": "Save NTOs", "solvation": "Solvation model", "solvent": "Solvent",
    "solvent_other": "Other solvent", "additional_input": "Additional input",
    "additional_keywords": "Additional Keywords",
}

# Fields chosen from a list; method_name depends on the family (method_choices).
CHOICES: Dict[str, Sequence[str]] = {
    "job_type": gk.JOB_TYPES, "opt_target": tuple(gk.OPT_TARGETS), "force_constants": tuple(gk.FORCE_CONSTANTS),
    "irc_direction": tuple(gk.IRC_DIRECTIONS), "irc_force_constants": tuple(gk.FORCE_CONSTANTS),
    "nmr_method": tuple(gk.NMR_METHODS), "state": gk.STATES, "state_spin": tuple(gk.STATE_SPINS),
    "spin": tuple(gk.SPIN_PREFIXES), "family": gk.FAMILIES, "basis": gk.BASES,
    "basis_diffuse": tuple(gk.DIFFUSE), "basis_polarization": tuple(gk.POLARIZATION),
    "dispersion": tuple(gk.DISPERSION), "geometry_source": tuple(gk.GEOMETRY_SOURCES),
    "print_level": gk.PRINT_LEVELS, "guess": tuple(gk.GUESS_METHODS), "pop_level": tuple(gk.POP_LEVELS),
    "solvation": tuple(gk.SOLVATION_MODELS), "solvent": gk.SOLVENTS,
}


def method_choices(family: str) -> Tuple[str, ...]:
    return tuple(gk.FAMILY_METHODS[family]) if family in gk.FAMILY_METHODS else (gk.NO_CHOICE,)


def choices(field_name: str, step: Optional[GaussianStepSpec] = None) -> Tuple[str, ...]:
    if field_name == "method_name":
        return method_choices((step or GaussianStepSpec()).family)
    return tuple(CHOICES[field_name])


# ----------------------------------------------------------------------
# Route of one step
# ----------------------------------------------------------------------
def no_basis(step: GaussianStepSpec) -> bool:
    """ZINDO, semi-empirical and compound methods take no basis set."""
    if step.state == gk.STATE_ZINDO:
        return True
    return step.state != gk.STATE_CIS and step.family in gk.NO_BASIS_FAMILIES


def state_options(step: GaussianStepSpec) -> List[str]:
    if step.state == gk.STATE_GROUND:
        return []
    options = [f"NStates={step.n_states}"]
    if step.root != 1:
        options.append(f"Root={step.root}")
    spin = gk.STATE_SPINS.get(step.state_spin)
    if spin:
        options.append(spin)
    return options


def method_keyword(step: GaussianStepSpec, with_prefix: bool = True) -> str:
    """'UB3LYP', 'CIS(NStates=3)', 'PM6', ... (the part before the '/')."""
    prefix = gk.SPIN_PREFIXES.get(step.spin, "") if with_prefix else ""
    if step.state in gk.MODEL_STATES:
        name, options = gk.STATE_KEYWORDS[step.state], state_options(step)
        return prefix + (f"{name}({','.join(options)})" if options else name)
    if step.family in gk.FAMILY_KEYWORDS:
        return prefix + gk.FAMILY_KEYWORDS[step.family]
    keyword = gk.FAMILY_METHODS.get(step.family, {}).get(step.method_name)
    return prefix + (keyword if keyword else step.method_other.strip())


def basis_keyword(step: GaussianStepSpec) -> str:
    if no_basis(step):
        return ""
    if step.basis == gk.OTHER:
        return step.basis_other.strip()
    if step.basis in gk.EXTENSIBLE_POPLE:
        diffuse = gk.DIFFUSE.get(step.basis_diffuse, "")
        polarization = gk.POLARIZATION.get(step.basis_polarization, "")
        return f"{step.basis[:-1]}{diffuse}G{polarization}"
    return step.basis


def model_token(step: GaussianStepSpec) -> str:
    method, basis = method_keyword(step), basis_keyword(step)
    return f"{method}/{basis}" if basis else method


def opt_options(step: GaussianStepSpec) -> List[str]:
    options = [o for o in (gk.OPT_TARGETS.get(step.opt_target), gk.FORCE_CONSTANTS.get(step.force_constants)) if o]
    if step.opt_tight:
        options.append("Tight")
    if step.opt_max_cycles:
        options.append(f"MaxCycles={step.opt_max_cycles}")
    if step.job_type == gk.JOB_SCAN or step.modredundant.strip():
        options.append("ModRedundant")
    return options


def _flags(step: GaussianStepSpec, table: Dict[str, str]) -> List[str]:
    return [option for name, option in table.items() if getattr(step, name)]


def job_keywords(step: GaussianStepSpec) -> List[RouteKeyword]:
    job, out = step.job_type, []
    if job in gk.OPTIMIZING_JOBS:
        out.append(RouteKeyword("Opt", opt_options(step)))
    if job in gk.FREQUENCY_JOBS:
        out.append(RouteKeyword("Freq", _flags(step, gk.FREQ_FLAG_FIELDS)))
    if job == gk.JOB_IRC:
        options = [o for o in (gk.IRC_DIRECTIONS.get(step.irc_direction),
                               gk.FORCE_CONSTANTS.get(step.irc_force_constants)) if o]
        if step.irc_max_points:
            options.append(f"MaxPoints={step.irc_max_points}")
        if step.irc_step_size:
            options.append(f"StepSize={step.irc_step_size}")
        out.append(RouteKeyword("IRC", options))
    if job == gk.JOB_STABILITY:
        out.append(RouteKeyword("Stable", ["Opt"] if step.stable_opt else []))
    if job == gk.JOB_NMR:
        method = gk.NMR_METHODS.get(step.nmr_method)
        out.append(RouteKeyword("NMR", ([method] if method else []) + (["SpinSpin"] if step.nmr_spinspin else [])))
    return out


def method_keywords(step: GaussianStepSpec) -> List[RouteKeyword]:
    """TD/TDA, dispersion and the coefficient IOp (CIS and ZINDO live in the model token)."""
    out = []
    if step.state in (gk.STATE_TD, gk.STATE_TDA):
        out.append(RouteKeyword(gk.STATE_KEYWORDS[step.state], state_options(step)))
    dispersion = gk.DISPERSION.get(step.dispersion)
    if dispersion:
        out.append(RouteKeyword("EmpiricalDispersion", [dispersion]))
    if step.print_all_coefficients:
        out.append(RouteKeyword("IOp", [gk.COEFFICIENTS_IOP]))
    return out


def general_keywords(step: GaussianStepSpec) -> List[RouteKeyword]:
    """General, Molecule (Geom), Guess, Pop and Solvation panels, in that order."""
    out = []
    scf = (["Tight"] if step.scf_tight else []) + (["QC"] if step.scf_qc else [])
    if scf:
        out.append(RouteKeyword("SCF", scf))
    if step.nosymm:
        out.append(RouteKeyword("NoSymm"))
    if step.polar:
        out.append(RouteKeyword("Polar"))
    source = gk.GEOMETRY_SOURCES.get(step.geometry_source)
    geom = ([source] if source else []) + (["Connectivity"] if step.connectivity else [])
    if geom:
        out.append(RouteKeyword("Geom", geom))
    guess = gk.GUESS_METHODS.get(step.guess)
    guess_options = ([guess] if guess else []) + _flags(step, gk.GUESS_FLAG_FIELDS)
    if guess_options:
        out.append(RouteKeyword("Guess", guess_options))
    level = gk.POP_LEVELS.get(step.pop_level)
    pop = ([level] if level else []) + _flags(step, gk.POP_FLAG_FIELDS)
    if pop:
        out.append(RouteKeyword("Pop", pop))
    model = gk.SOLVATION_MODELS.get(step.solvation)
    if model:
        solvent = step.solvent_other.strip() if step.solvent == gk.OTHER else step.solvent
        out.append(RouteKeyword("SCRF", [model, f"Solvent={solvent}"]))
    return out


def panel_keywords(step: GaussianStepSpec) -> List[RouteKeyword]:
    return job_keywords(step) + method_keywords(step) + general_keywords(step)


def route_keywords(step: GaussianStepSpec) -> Tuple[List[RouteKeyword], List[str]]:
    """The panels' keywords merged with Additional Keywords; second item = merge errors."""
    return merge_keywords(panel_keywords(step), step.additional_keywords)


def route_lines(step: GaussianStepSpec) -> List[str]:
    """The route section, wrapped below 80 columns (merge errors are ignored here; prepare reports them)."""
    keywords, _ = route_keywords(step)
    return fi.build_route(step.print_level, method_keyword(step), basis_keyword(step),
                          " ".join(format_keyword(k) for k in keywords))


# ----------------------------------------------------------------------
# Scan / freeze lines (Opt=ModRedundant)
# ----------------------------------------------------------------------
_MODRED_COUNTS = {letter: count for letter, count in gk.MODRED_TYPES.values()}


def parse_modredundant(text: str, n_atoms: Optional[int] = None) -> List[Tuple[str, Tuple[int, ...], str]]:
    """'B 1 2 S 10 0.05' / 'D 1 2 3 4 F' lines -> (letter, atoms, 'S'|'F'); ValueError names the line."""
    parsed = []
    for number, raw in enumerate((text or "").splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        where = f"Scan/freeze line {number} ('{line}')"
        tokens = line.split()
        letter = tokens[0].upper()
        if letter not in _MODRED_COUNTS:
            raise ValueError(f"{where}: start with B, A or D (bond, angle, dihedral).")
        count = _MODRED_COUNTS[letter]
        atom_tokens = tokens[1:1 + count]
        if len(atom_tokens) != count or not all(t.isdigit() for t in atom_tokens):
            raise ValueError(f"{where}: {letter} needs {count} atom numbers.")
        atoms = tuple(int(t) for t in atom_tokens)
        if len(set(atoms)) != count or min(atoms) < 1:
            raise ValueError(f"{where}: give {count} different atoms, numbered from 1.")
        if n_atoms is not None and max(atoms) > n_atoms:
            raise ValueError(f"{where}: atom {max(atoms)} does not exist (the molecule has {n_atoms} atoms).")
        rest = [t.upper() for t in tokens[1 + count:]]
        if rest == ["F"]:
            action = "F"
        elif len(rest) == 3 and rest[0] == "S":
            if not rest[1].isdigit() or int(rest[1]) < 1:
                raise ValueError(f"{where}: the number of steps must be a positive integer.")
            try:
                size = float(rest[2])
            except ValueError:
                size = 0.0
            if size == 0.0 or not math.isfinite(size):
                raise ValueError(f"{where}: the step size must be a non-zero number.")
            action = "S"
        else:
            raise ValueError(f"{where}: end with F (freeze) or S <steps> <step size> (scan).")
        parsed.append((letter, atoms, action))
    return parsed


def modredundant_line(kind: str, atoms_text: str, action: str, steps_text: str, size_text: str) -> str:
    """The Molecule panel's helper row -> one line, checked by the same rules as typed lines."""
    letter, _ = gk.MODRED_TYPES[kind]
    atoms = [t for t in re.split(r"[\s,]+", (atoms_text or "").strip()) if t]
    code = gk.MODRED_ACTIONS[action]
    words = [letter, *atoms, code] + ([steps_text.strip(), size_text.strip()] if code == "S" else [])
    line = " ".join(w for w in words if w)
    parse_modredundant(line)
    return line


# ----------------------------------------------------------------------
# The job
# ----------------------------------------------------------------------
INPUT_SUFFIX = ".gjf"
SCRIPT_HEADER = ("# gQTEA Gaussian 16 input builder: convert the checkpoints to formatted",
                 "# checkpoints for the 3D Molecular Viewer's analysis tabs. Run it where Gaussian")
ASCII_FIELDS = ("title", "extra_link0", "additional_input", "additional_keywords", "method_other",
                "basis_other", "solvent_other")
FORCE_CONSTANT_OPTIONS = {"CALCFC", "CALCALL", "READFC", "RCFC"}


@dataclass
class GaussianJobSpec:
    """The whole job: file names, the step-1 geometry and one GaussianStepSpec per linked step."""
    stem: str = fi.DEFAULT_STEM
    output_dir: str = ""
    geometry: str = ""
    steps: List[GaussianStepSpec] = field(default_factory=lambda: [GaussianStepSpec()])
    write_formchk: bool = True


@dataclass
class PreparedStep:
    """One validated step, ready to write."""
    link0: List[str]
    route: List[str]
    title: Optional[str]            # None with Geom=AllCheck (no title section)
    charge_line: Optional[str]      # None with Geom=AllCheck
    coordinates: bool
    connectivity: List[str]
    modredundant: List[str]
    additional_input: List[str]
    chk: Optional[str]


@dataclass
class GaussianJob:
    spec: GaussianJobSpec
    atoms: Optional[List[fi.Atom]]
    steps: List[PreparedStep]
    input_path: Path
    script_path: Optional[Path]
    chk_files: List[str]
    warnings: List[str] = field(default_factory=list)
    fukui: Optional[fi.FukuiJob] = None


def _relative_name(value: str, label: str) -> str:
    name = (value or "").strip()
    if not name or any(ch.isspace() for ch in name) or not name.isascii():
        raise ValueError(f"{label} '{name}': use a file name without blanks (ASCII letters, digits, . _ -).")
    if PureWindowsPath(name).drive or "/" in name or "\\" in name:
        raise ValueError(f"{label} '{name}': give a file name relative to the input's folder "
                         "(the job runs under Linux from there).")
    return name


def _one_keyword(text: str) -> bool:
    text = (text or "").strip()
    return bool(text) and len(text.split()) == 1 and "/" not in text and "," not in text


def _checkpoint_readers(keywords: Sequence[RouteKeyword]) -> List[str]:
    """What in the route reads the checkpoint: 'Geom=AllCheck', 'Guess=Read', 'Opt=ReadFC', ..."""
    readers = []
    for keyword in keywords:
        options = {o.upper() for o in keyword.options}
        if keyword.key == "GEOM":
            readers += [f"Geom={o}" for o in ("AllCheck", "Check") if o.upper() in options]
        elif keyword.key == "GUESS":
            readers += [f"Guess={o}" for o in gk.READ_GUESSES if o.upper() in options]
        elif keyword.key in ("OPT", "IRC"):
            readers += [f"{keyword.name}={o}" for o in gk.READ_FORCE_CONSTANTS if o.upper() in options]
    return readers


def connectivity_lines(atoms: Sequence[fi.Atom]) -> List[str]:
    """The Geom=Connectivity section, single bonds from the distances (the 3D viewer's rule)."""
    from molecularFileFormats import detect_bonds    # numpy + scipy, needed only for this option
    bonds = detect_bonds([a[0] for a in atoms], [a[1:] for a in atoms])
    partners: Dict[int, List[int]] = {i: [] for i in range(len(atoms))}
    for i, j, _ in bonds:
        partners[i].append(j)
    return [" ".join([f" {i + 1}"] + [f"{j + 1} 1.0" for j in sorted(partners[i])]) for i in range(len(atoms))]


def _trim_blank_lines(text: str) -> List[str]:
    lines = [line.rstrip() for line in (text or "").splitlines()]
    while lines and not lines[0]:
        lines.pop(0)
    while lines and not lines[-1]:
        lines.pop()
    return lines


def _prepare_step(step: GaussianStepSpec, stem: str, atoms, written: set,
                  write_formchk: bool) -> Tuple[List[str], List[str], Optional[PreparedStep]]:
    """(errors, warnings, prepared step or None) of one step; messages carry no step prefix."""
    errors: List[str] = []
    warnings: List[str] = []

    def attempt(fn, *args):
        try:
            return fn(*args)
        except ValueError as exc:
            errors.append(str(exc))
            return None

    for name, values in CHOICES.items():
        if getattr(step, name) not in values:
            errors.append(f"{FIELD_LABELS[name]}: '{getattr(step, name)}' is not one of the choices.")
    if step.family in gk.FAMILIES and step.method_name not in method_choices(step.family):
        errors.append(f"Method: '{step.method_name}' is not a {step.family} method.")
    for name in ASCII_FIELDS:
        if not getattr(step, name).isascii():
            errors.append(f"{FIELD_LABELS[name]}: use ASCII characters only (Gaussian input files are ASCII).")
    if errors:
        return errors, warnings, None

    # Link 0
    memory = attempt(fi.parse_memory, step.memory)
    if not isinstance(step.nproc, int) or step.nproc < 1:
        errors.append(f"%NProcShared must be a positive integer (got {step.nproc}).")
    chk = attempt(_relative_name, step.chk.strip() or f"{stem}.chk", "%Chk") if step.write_chk else None
    old_chk = attempt(_relative_name, step.old_chk, "%OldChk") if step.old_chk.strip() else None
    extra_link0 = [line.strip() for line in step.extra_link0.splitlines() if line.strip()]
    errors += [f"Extra Link 0 line '{line}' must start with %." for line in extra_link0 if not line.startswith("%")]

    # Method and basis
    excited = step.state != gk.STATE_GROUND
    if (step.state not in gk.MODEL_STATES and step.family == gk.FAMILY_DFT and step.method_name == gk.OTHER
            and not _one_keyword(step.method_other)):
        errors.append("Other method: enter one keyword, e.g. B97D3 (no blanks, no '/').")
    if not no_basis(step):
        if step.basis == gk.OTHER and not _one_keyword(step.basis_other):
            errors.append("Other basis set: enter one keyword, e.g. def2TZVPD (no blanks, no '/').")
        if step.basis in gk.GENERAL_BASES and not step.additional_input.strip():
            errors.append(f"Basis set {step.basis}: put the basis-set section in the Add. Inp. panel.")
    if step.solvation != gk.SOLVATION_NONE and step.solvent == gk.OTHER and not _one_keyword(step.solvent_other):
        errors.append("Other solvent: enter one Gaussian solvent name, e.g. Heptane.")
    if excited:
        if not isinstance(step.n_states, int) or step.n_states < 1:
            errors.append("Number of states must be at least 1.")
        elif not isinstance(step.root, int) or not 1 <= step.root <= step.n_states:
            errors.append(f"State of interest (Root) must be between 1 and {step.n_states}.")
        if step.state == gk.STATE_TD and step.family not in (gk.FAMILY_HF, gk.FAMILY_DFT):
            errors.append("TD-SCF needs a Hartree-Fock or DFT method.")
        if step.state == gk.STATE_TDA and step.family != gk.FAMILY_DFT:
            errors.append("TDA needs a DFT method.")
    elif step.print_all_coefficients:
        errors.append("Print all excitation coefficients needs an excited state (TD-SCF, TDA, CIS or ZINDO).")

    # Geometry, charge, multiplicity
    source = step.geometry_source
    prefix = gk.SPIN_PREFIXES[step.spin]
    mult = step.multiplicity
    if step.connectivity and source != gk.GEOM_COORDINATES:
        errors.append("Write connectivity needs coordinates; this step reads its geometry from the checkpoint.")
    if source != gk.GEOM_ALLCHECK:
        electrons = None if atoms is None else nuclear_charge(a[0] for a in atoms) - step.charge
        if mult is None:
            allowed = (f" ({electrons} electrons allow "
                       f"{', '.join(str(m) for m in allowed_multiplicities(electrons, 3))}, ...)"
                       if electrons is not None and electrons >= 1 else "")
            errors.append(f"Multiplicity: enter it; it is not assumed{allowed}.")
        elif electrons is not None:
            problem = multiplicity_problem(electrons, mult) if electrons >= 1 else "no electrons are left"
            if problem:
                errors.append(f"Multiplicity {mult} is impossible ({problem}).")
    if mult is not None and mult > 1 and prefix == "R":
        errors.append(f"The R (restricted) prefix needs multiplicity 1 (got {mult}); use RO or U.")
    if gk.STATE_SPINS.get(step.state_spin) and (prefix in ("U", "RO") or (mult is not None and mult != 1)):
        errors.append(f"States = {step.state_spin}: Singlets, Triplets and 50-50 need a closed-shell reference "
                      "(multiplicity 1, no U or RO prefix).")
    title = step.title.strip() or fi._title_text(stem)
    bad = sorted(set(title) & set(fi.TITLE_FORBIDDEN))
    if bad:
        errors.append(f"Title: Gaussian does not allow {' '.join(bad)} in the title.")

    # Scan / freeze lines
    parsed: List[Tuple[str, Tuple[int, ...], str]] = []
    modred_ok = True
    if step.modredundant.strip():
        if step.job_type not in gk.OPTIMIZING_JOBS:
            errors.append("Scan/freeze lines need the Optimization, Opt+Freq or Scan job type.")
            modred_ok = False
        else:
            parsed = attempt(parse_modredundant, step.modredundant, len(atoms) if atoms else None)
            modred_ok = parsed is not None
            parsed = parsed or []
            if any(p[2] == "S" for p in parsed) and step.job_type != gk.JOB_SCAN:
                errors.append("Scan lines (S) need the Scan job type; Optimization takes freeze lines (F) only.")
    if step.job_type == gk.JOB_SCAN and modred_ok and not any(p[2] == "S" for p in parsed):
        errors.append("The Scan job type needs at least one scan line (S); add one in the Molecule panel.")

    # Route and checkpoints
    keywords, merge_errors = route_keywords(step)
    errors += merge_errors
    readers = _checkpoint_readers(keywords)
    if readers and old_chk is None and (chk is None or chk.lower() not in written):
        errors.append(f"{', '.join(readers)} reads a checkpoint, but no earlier step writes "
                      f"{chk or 'one'}; set %OldChk to an existing checkpoint.")
    options = {k.key: {o.upper() for o in k.options} for k in keywords}
    if "TS" in options.get("OPT", set()) and not FORCE_CONSTANT_OPTIONS & options["OPT"]:
        warnings.append("Optimization to a TS without force constants (CalcFC, CalcAll or ReadFC) often fails; "
                        "CalcFC is recommended.")
    if "IRC" in options and not FORCE_CONSTANT_OPTIONS & options["IRC"]:
        warnings.append("IRC without force constants: give CalcFC, CalcAll or ReadFC.")
    if write_formchk and (excited or step.family in gk.POST_SCF_FAMILIES) and "DENSITY" not in options:
        warnings.append("The .fchk will hold the SCF density, not the excited-state or correlated one; add "
                        "Density=Current to the Additional Keywords if an analysis needs it.")
    if errors:
        return errors, warnings, None

    link0 = ([f"%OldChk={old_chk}"] if old_chk else []) + ([f"%Chk={chk}"] if chk else []) + \
        [f"%Mem={memory}", f"%NProcShared={step.nproc}"] + extra_link0
    route = fi.build_route(step.print_level, method_keyword(step), basis_keyword(step),
                           " ".join(format_keyword(k) for k in keywords))
    all_check = source == gk.GEOM_ALLCHECK
    return errors, warnings, PreparedStep(
        link0=link0, route=route, title=None if all_check else title,
        charge_line=None if all_check else f"{step.charge} {mult}",
        coordinates=source == gk.GEOM_COORDINATES,
        connectivity=connectivity_lines(atoms) if step.connectivity else [],
        modredundant=[line.strip() for line in step.modredundant.splitlines() if line.strip()],
        additional_input=_trim_blank_lines(step.additional_input), chk=chk)


def linked_step(step: GaussianStepSpec) -> GaussianStepSpec:
    """A new step continuing from ``step``: same method and Link 0, geometry and orbitals from its checkpoint.

    Charge and multiplicity are cleared: Geom=AllCheck reads them from the checkpoint, and a step
    switched to another source must be given them (they are never assumed).
    """
    return dataclasses.replace(step, job_type=gk.JOB_ENERGY, geometry_source=gk.GEOM_ALLCHECK,
                               guess="Read checkpoint", title="", modredundant="", additional_input="",
                               connectivity=False, old_chk="", charge=0, multiplicity=None)


def step_label(index: int, step: GaussianStepSpec) -> str:
    """'2: TD-SCF B3LYP/6-31G(d), Geom=AllCheck' for the step list (index is 1-based)."""
    what = step.state if step.state != gk.STATE_GROUND and step.job_type == gk.JOB_ENERGY else step.job_type
    source = gk.GEOMETRY_SOURCES.get(step.geometry_source)
    return f"{index}: {what} {model_token(step)}" + (f", Geom={source}" if source else "")


class GaussianInputBuilder:
    """Validation, the input text, the formchk script and the files."""

    def __init__(self):
        self._fukui = fi.FukuiInputBuilder()

    def prepare(self, spec: GaussianJobSpec) -> GaussianJob:
        """Validate everything; ValueError lists every problem, one '- ' line each."""
        if any(s.job_type == gk.JOB_FUKUI for s in spec.steps):
            return self._prepare_fukui(spec)
        errors: List[str] = []
        warnings: List[str] = []
        stem = None
        try:
            stem = fi.validate_stem(spec.stem)
        except ValueError as exc:
            errors.append(str(exc))
        folder_text = (spec.output_dir or "").strip()
        folder = Path(folder_text).resolve() if folder_text else None
        if folder is None or not folder.is_dir():
            errors.append(f"Output folder '{folder_text}' does not exist; choose an existing folder."
                          if folder_text else "Output folder: choose where the input files are written.")
        if not spec.steps:
            errors.append("The job has no steps.")
        atoms = None
        if any(s.geometry_source == gk.GEOM_COORDINATES for s in spec.steps):
            try:
                atoms = fi.parse_geometry(spec.geometry)
            except ValueError as exc:
                errors.append(str(exc))
        name = stem or fi.DEFAULT_STEM
        written: set = set()
        prepared: List[PreparedStep] = []
        for index, step in enumerate(spec.steps, start=1):
            step_errors, step_warnings, item = _prepare_step(step, name, atoms, written, spec.write_formchk)
            errors += [f"Step {index}: {e}" for e in step_errors]
            warnings += [f"Step {index}: {w}" for w in step_warnings]
            if item is not None:
                prepared.append(item)
            if step.write_chk:
                written.add((step.chk.strip() or f"{name}.chk").lower())
        if errors:
            raise ValueError("\n".join(f"- {e}" for e in errors))
        chk_files = list(dict.fromkeys(p.chk for p in prepared if p.chk))
        if len([p for p in prepared if p.chk]) > len(chk_files):
            warnings.append("Several steps write the same checkpoint: it keeps only the last step's results, "
                            "and that is what the formchk script converts.")
        script_path = None
        if spec.write_formchk and chk_files:
            script_path = folder / f"{stem}{fi.SCRIPT_SUFFIX}"
        elif spec.write_formchk:
            warnings.append("No step writes a checkpoint (%Chk), so no formchk script is written.")
        return GaussianJob(spec, atoms, prepared, folder / f"{stem}{INPUT_SUFFIX}", script_path, chk_files, warnings)

    def _prepare_fukui(self, spec: GaussianJobSpec) -> GaussianJob:
        if len(spec.steps) != 1:
            raise ValueError("- The Fukui functions job type is already three linked jobs (N, N+p, N-q): "
                             "it must be the only step.")
        step = spec.steps[0]
        reference = fukui_defaults(step)
        errors = [f"The Fukui functions job type sets {FIELD_LABELS[name]} itself; choose the job type again "
                  "to reset it." for name in FUKUI_FIXED if getattr(step, name) != getattr(reference, name)]
        keywords, merge_errors = fukui_keywords(step)
        errors += merge_errors
        if any(k.key == "OPT" for k in keywords):
            errors.append("The Fukui functions job type: set the N-state optimization with the Optimization "
                          "options of the Job Type panel, not in Additional Keywords.")
        if errors:
            raise ValueError("\n".join(f"- {e}" for e in errors))
        job = self._fukui.prepare(fukui_spec(spec))
        return GaussianJob(spec, job.atoms, [], job.input_path, job.script_path, [s.chk for s in job.states],
                           list(job.warnings), fukui=job)

    def input_text(self, job: GaussianJob) -> str:
        if job.fukui is not None:
            return self._fukui.input_text(job.fukui)
        blocks = []
        for step in job.steps:
            lines = list(step.link0) + list(step.route) + [""]
            if step.title is not None:
                lines += [step.title, "", step.charge_line]
                if step.coordinates:
                    lines += fi.format_geometry(job.atoms)
                lines.append("")
            for section in (step.connectivity, step.modredundant, step.additional_input):
                if section:
                    lines += list(section) + [""]
            blocks.append("\n".join(lines))
        return "\n--Link1--\n".join(blocks) + "\n"

    def script_text(self, job: GaussianJob) -> Optional[str]:
        if job.fukui is not None:
            return self._fukui.script_text(job.fukui)
        if job.script_path is None:
            return None
        return fi.formchk_script(job.script_path.name, SCRIPT_HEADER,
                                 [(chk, Path(chk).stem + fi.FCHK_EXTENSION) for chk in job.chk_files])

    def write_files(self, job: GaussianJob) -> List[Path]:
        """Write the .gjf (and the .sh) with LF line endings; ValueError on any OS error."""
        if job.fukui is not None:
            return self._fukui.write_files(job.fukui)
        written = []
        files = [(job.input_path, self.input_text(job))]
        if job.script_path is not None:
            files.append((job.script_path, self.script_text(job)))
        for path, text in files:
            try:
                data = text.encode("ascii")     # before open(): a failure must not truncate an existing file
            except UnicodeEncodeError as exc:
                raise ValueError(f"Could not write {path}: {exc}") from exc
            try:
                with open(path, "wb") as handle:
                    handle.write(data)
            except OSError as exc:
                raise ValueError(f"Could not write {path}: {exc}") from exc
            written.append(path)
        return written


# ----------------------------------------------------------------------
# The Fukui Job Type (delegated to fukuiInput, unchanged)
# ----------------------------------------------------------------------
# Fields the Fukui linked job writes itself; they must hold fukui_defaults' values.
FUKUI_FIXED = ("opt_target", "force_constants", "state", "print_all_coefficients", "spin", "geometry_source",
               "modredundant", "title", "write_chk", "chk", "old_chk", "extra_link0", "connectivity", "guess",
               "guess_mix", "guess_always", "guess_only", "guess_save", "additional_input")


def fukui_defaults(step: GaussianStepSpec) -> GaussianStepSpec:
    """The step as the Fukui Job Type uses it: the former builder's SCF=Tight NoSymm, unused controls reset."""
    force_constants = step.force_constants if step.force_constants in gk.FUKUI_FORCE_CONSTANTS else gk.FC_DEFAULT
    return dataclasses.replace(
        step, job_type=gk.JOB_FUKUI, opt_target=gk.OPT_MINIMUM, force_constants=force_constants,
        state=gk.STATE_GROUND, print_all_coefficients=False, spin="Default", geometry_source=gk.GEOM_COORDINATES,
        modredundant="", title="", write_chk=True, chk="", old_chk="", extra_link0="", connectivity=False,
        guess="Default", guess_mix=False, guess_always=False, guess_only=False, guess_save=False,
        additional_input="", scf_tight=True, nosymm=True)


def fukui_opt_options(step: GaussianStepSpec) -> List[str]:
    """Options of the N-state optimization (towards a minimum; no ReadFC, no ModRedundant)."""
    options = [o for o in (gk.FORCE_CONSTANTS.get(step.force_constants),) if o]
    if step.opt_tight:
        options.append("Tight")
    if step.opt_max_cycles:
        options.append(f"MaxCycles={step.opt_max_cycles}")
    return options


def fukui_keywords(step: GaussianStepSpec) -> Tuple[List[RouteKeyword], List[str]]:
    """Keywords shared by the three Fukui jobs (panels + Additional Keywords, merged as in spec 5.1)."""
    return merge_keywords(method_keywords(step) + general_keywords(step), step.additional_keywords)


def fukui_spec(spec: GaussianJobSpec) -> fi.FukuiJobSpec:
    """The panels of the single Fukui step as the FukuiJobSpec that fukuiInput validates and writes."""
    step = spec.steps[0]
    shared = [format_keyword(k) for k in fukui_keywords(step)[0]]
    return fi.FukuiJobSpec(
        stem=(spec.stem or "").strip(), output_dir=(spec.output_dir or "").strip(), geometry=spec.geometry,
        charge=step.charge, multiplicity=step.multiplicity,
        plus_multiplicity=step.fukui_plus_multiplicity, minus_multiplicity=step.fukui_minus_multiplicity,
        memory=step.memory, nproc=step.nproc, print_level=step.print_level,
        method=method_keyword(step, with_prefix=False), basis=basis_keyword(step),
        keywords=" ".join(shared),
        opt_keyword=format_keyword(RouteKeyword("Opt", fukui_opt_options(step))),
        p=step.fukui_p, q=step.fukui_q, guess_read=step.fukui_guess_read)


# ----------------------------------------------------------------------
# Live texts of the window (never raise)
# ----------------------------------------------------------------------
def route_preview(spec: GaussianJobSpec, index: int) -> str:
    """The route line(s) of step ``index`` (0-based), with merge problems on extra lines."""
    step = spec.steps[index]
    try:
        if step.job_type == gk.JOB_FUKUI:
            f = fukui_spec(spec)
            args = (f.print_level, f.method, f.basis, f.keywords)
            ion_extra = ("Geom=Check", "Guess=Read") if f.guess_read else ("Geom=Check",)
            errors = fukui_keywords(step)[1]
            return ("job 1:     " + " ".join(fi.build_route(*args, (f.opt_keyword,))) + "\n"
                    "jobs 2-3:  " + " ".join(fi.build_route(*args, ion_extra))
                    + "".join(f"\n{e}" for e in errors))
        keywords, errors = route_keywords(step)
        route = fi.build_route(step.print_level, method_keyword(step), basis_keyword(step),
                               " ".join(format_keyword(k) for k in keywords))
        return " ".join(route) + "".join(f"\n{e}" for e in errors)
    except (KeyError, ValueError) as exc:
        return f"Cannot build the route yet: {exc}"


def _spin_text(electrons: int, multiplicity: Optional[int]) -> str:
    if multiplicity is None:
        allowed = ", ".join(str(m) for m in allowed_multiplicities(electrons, 3))
        return f"multiplicity ? (enter it: {electrons} electrons allow {allowed}, ...)"
    problem = multiplicity_problem(electrons, multiplicity) if electrons >= 1 else "no electrons left"
    return f"multiplicity {multiplicity}" + (f"   <- impossible: {problem}" if problem else "")


def molecule_summary(spec: GaussianJobSpec, index: int) -> str:
    """Charge, electron count and multiplicity of step ``index`` (the Fukui job: of its three states)."""
    step = spec.steps[index]
    if step.geometry_source == gk.GEOM_ALLCHECK:
        return "Geometry, charge and multiplicity come from the checkpoint."
    if not (spec.geometry or "").strip():
        return "Enter or load a geometry to see the electron count."
    try:
        atoms = fi.parse_geometry(spec.geometry, check_distances=False)
    except ValueError as exc:
        return f"Cannot count electrons yet: {exc}"
    z_total = nuclear_charge(a[0] for a in atoms)
    if step.job_type != gk.JOB_FUKUI:
        electrons = z_total - step.charge
        return (f"{len(atoms)} atoms, charge {step.charge:+d}, {electrons} electrons, "
                f"{_spin_text(electrons, step.multiplicity)}").replace("charge +0", "charge 0")
    stem = (spec.stem or "").strip() or fi.DEFAULT_STEM
    lines = []
    multiplicities = (step.multiplicity, step.fukui_plus_multiplicity, step.fukui_minus_multiplicity)
    for (label, tag, _, k, sign), mult in zip(fi.state_tags(step.fukui_p, step.fukui_q), multiplicities):
        charge = step.charge + sign * k
        electrons = z_total - charge
        lines.append(f"{label:<5} charge {charge:+d}, {electrons} electrons -> {stem}_{tag}.chk, "
                     f"{_spin_text(electrons, mult)}")
    return "\n".join(lines).replace("charge +0", "charge 0")


# ----------------------------------------------------------------------
# Open input...: a Gaussian input back into the panels
# ----------------------------------------------------------------------
GAUSSIAN_INPUT_SUFFIXES = (".gjf", ".com", ".gau")
_POPLE_RE = re.compile(r"^(6-31|6-311)(\+{0,2})G(\*{1,2}|\([^)]*\))?$", re.IGNORECASE)
_STAR_POLARIZATION = {"*": "(d)", "**": "(d,p)"}
_INTEGER_RE = re.compile(r"[+-]?\d+")


def _int_option(option: str) -> Optional[int]:
    value = option.split("=", 1)[1].strip() if "=" in option else ""
    return int(value) if value.isdigit() and int(value) > 0 else None


def _bare_method(upper: str, original: Optional[str] = None) -> Optional[dict]:
    """Panel fields of a method keyword without prefix (upper case), or None if no panel knows it.

    ``original`` is the same text in the file's own case (CIS/ZINDO options keep their spelling)."""
    families = {v.upper(): k for k, v in gk.FAMILY_KEYWORDS.items()}
    families["MP4"] = gk.FAMILY_MP4
    if upper in families:
        return {"family": families[upper], "method_name": gk.NO_CHOICE}
    for family, table in gk.FAMILY_METHODS.items():
        back = gk.reverse(table)
        if upper in back:
            return {"family": family, "method_name": back[upper]}
    name = keyword_text_name(upper)
    states = {"CIS": gk.STATE_CIS, "ZINDO": gk.STATE_ZINDO}
    if name in states:
        return {"state": states[name], "family": gk.FAMILY_HF, "method_name": gk.NO_CHOICE,
                "_options": parse_keyword(original or upper).options}
    return None


def _method_fields(text: str) -> Optional[dict]:
    """'UB3LYP' -> {'spin': 'Unrestricted (U)', 'family': 'DFT', 'method_name': 'B3LYP'}; None if unknown."""
    upper = text.upper()
    fields = _bare_method(upper, text)
    if fields is not None:
        return {"spin": "Default", **fields}
    for prefix in ("RO", "R", "U"):
        if upper.startswith(prefix):
            fields = _bare_method(upper[len(prefix):], text[len(prefix):])
            if fields is not None:
                return {"spin": gk.reverse(gk.SPIN_PREFIXES)[prefix], **fields}
    return None


class _RouteReader:
    """Puts route tokens back into a GaussianStepSpec; what no panel models is kept as extra keywords."""

    def __init__(self, step: GaussianStepSpec, notes: List[str], where: str):
        self.step, self.notes, self.where = step, notes, where
        self.jobs: List[str] = []
        self.job_tokens: List[Tuple[str, str, Optional[int]]] = []   # (job keyword, token, index in extra)
        self.modredundant = False
        self.extra: List[str] = []
        self.handlers = {"OPT": self._opt, "FREQ": self._freq, "IRC": self._irc, "STABLE": self._stable,
                         "NMR": self._nmr, "TD": self._td, "TDA": self._td,
                         "EMPIRICALDISPERSION": self._dispersion, "IOP": self._iop, "SCF": self._scf,
                         "NOSYMM": self._switch("nosymm"), "POLAR": self._switch("polar"), "GEOM": self._geom,
                         "GUESS": self._guess, "POP": self._pop, "SCRF": self._scrf}

    def read(self, tokens: Sequence[str]) -> None:
        model_seen = False
        for token in tokens:
            if not model_seen and self._model(token):
                model_seen = True
                continue
            keyword = parse_keyword(token)
            handler = self.handlers.get(keyword.key)
            leftover = handler(keyword) if handler else None
            index = len(self.extra)
            if leftover is None:
                self.extra.append(token)
            elif leftover:
                self.extra.append(format_keyword(RouteKeyword(keyword.name, leftover)))
            if keyword.key in self.JOB_KEYS:
                self.job_tokens.append((keyword.key, token, index if len(self.extra) > index else None))
        if not model_seen:
            self.step.family, self.step.method_name, self.step.basis = gk.FAMILY_HF, gk.NO_CHOICE, "STO-3G"
            self.notes.append(f"{self.where}: the route names no method; Gaussian's default HF/STO-3G is used.")

    JOB_KEYS = ("OPT", "FREQ", "IRC", "STABLE", "NMR")

    def kept_extra(self) -> List[str]:
        return [t for t in self.extra if t is not None]

    def _keep_unused_jobs(self, used: Sequence[str]) -> None:
        """A job keyword the chosen job type does not write goes back as an Additional Keyword."""
        dropped = []
        for key, token, index in self.job_tokens:
            if key in used:
                continue
            if index is not None:
                self.extra[index] = None             # its leftover options travel inside the whole token
            self.extra.append(token)
            dropped.append(token)
        if dropped:
            self.notes.append(f"{self.where}: {', '.join(dropped)} cannot share the job type with the other job "
                              "keywords; kept as Additional Keywords.")

    def job_type(self, scan_lines: bool) -> str:
        jobs = self.jobs
        kinds = [j for j in ("OPT", "IRC", "STABLE", "NMR") if j in jobs]
        if "OPT" in jobs:
            used = ("OPT", "FREQ")
        elif "FREQ" in jobs:
            used = ("FREQ",)
        else:
            used = tuple(kinds[:1])
        self._keep_unused_jobs(used)
        if "OPT" in jobs:
            return gk.JOB_OPT_FREQ if "FREQ" in jobs else (gk.JOB_SCAN if scan_lines else gk.JOB_OPT)
        if "FREQ" in jobs:
            return gk.JOB_FREQ
        return {"IRC": gk.JOB_IRC, "STABLE": gk.JOB_STABILITY, "NMR": gk.JOB_NMR}.get(kinds[0] if kinds else "",
                                                                                    gk.JOB_ENERGY)

    # --- model chemistry ---
    def _model(self, token: str) -> bool:
        method_text, slash, basis_text = token.partition("/")
        if slash and method_text.count("(") != method_text.count(")"):
            return False                                     # IOp(9/40=4)
        fields = _method_fields(method_text)
        if fields is None:
            if not slash:
                return False
            fields = {"spin": "Default", "family": gk.FAMILY_DFT, "method_name": gk.OTHER,
                      "method_other": method_text}
        options = fields.pop("_options", [])
        for name, value in fields.items():
            setattr(self.step, name, value)
        dropped = self._state_options(options) if self.step.state in gk.MODEL_STATES else []
        if dropped:
            self.notes.append(f"{self.where}: option(s) {', '.join(dropped)} of {self.step.state} have no panel "
                              "and were dropped; add them back by hand if they matter.")
        if slash:
            self._basis(basis_text)
        elif not no_basis(self.step):
            self.step.basis = "STO-3G"
            self.notes.append(f"{self.where}: the route names no basis set; Gaussian's default STO-3G is used.")
        return True

    def _basis(self, text: str) -> None:
        step = self.step
        match = _POPLE_RE.match(text)
        if match:
            polarization = (match.group(3) or "")
            polarization = _STAR_POLARIZATION.get(polarization, polarization).lower()
            back = {v.lower(): k for k, v in gk.POLARIZATION.items()}
            if polarization in back:
                step.basis = match.group(1) + "G"
                step.basis_diffuse = {v: k for k, v in gk.DIFFUSE.items()}[match.group(2)]
                step.basis_polarization = back[polarization]
                return
        known = {b.upper(): b for b in gk.POPLE_BASES + gk.OTHER_BASES}
        if text.upper() in known:
            step.basis = known[text.upper()]
        else:
            step.basis, step.basis_other = gk.OTHER, text

    def _state_options(self, options: Sequence[str]) -> List[str]:
        step, rest = self.step, []
        spins = gk.reverse(gk.STATE_SPINS)
        step.n_states, step.root, step.state_spin = gk.DEFAULT_N_STATES, 1, "Default"
        for option in options:
            key = option_key(option)
            if key == "NSTATES" and _int_option(option):
                step.n_states = _int_option(option)
            elif key == "ROOT" and _int_option(option):
                step.root = _int_option(option)
            elif option.upper() in spins:
                step.state_spin = spins[option.upper()]
            else:
                rest.append(option)
        return rest

    # --- keywords (each returns the options it does not model; None = keep the token as it is) ---
    def _opt(self, keyword: RouteKeyword) -> List[str]:
        self.jobs.append("OPT")
        rest, fc = [], gk.reverse(gk.FORCE_CONSTANTS)
        for option in keyword.options:
            key = option_key(option)
            if key == "TS":
                self.step.opt_target = gk.OPT_TS
            elif key in fc:
                self.step.force_constants = fc[key]
            elif key == "TIGHT":
                self.step.opt_tight = True
            elif key == "MAXCYCLES" and _int_option(option):
                self.step.opt_max_cycles = _int_option(option)
            elif key == "MODREDUNDANT":
                self.modredundant = True
            else:
                rest.append(option)
        return rest

    def _flags(self, keyword: RouteKeyword, table: Dict[str, str]) -> List[str]:
        back = {option.upper(): name for name, option in table.items()}
        rest = []
        for option in keyword.options:
            if option.upper() in back:
                setattr(self.step, back[option.upper()], True)
            else:
                rest.append(option)
        return rest

    def _freq(self, keyword: RouteKeyword) -> List[str]:
        self.jobs.append("FREQ")
        return self._flags(keyword, gk.FREQ_FLAG_FIELDS)

    def _irc(self, keyword: RouteKeyword) -> List[str]:
        self.jobs.append("IRC")
        step, rest = self.step, []
        step.irc_force_constants = gk.FC_DEFAULT
        directions, fc = gk.reverse(gk.IRC_DIRECTIONS), gk.reverse(gk.FORCE_CONSTANTS)
        for option in keyword.options:
            key = option_key(option)
            if key in directions:
                step.irc_direction = directions[key]
            elif key in fc:
                step.irc_force_constants = fc[key]
            elif key == "MAXPOINTS" and _int_option(option):
                step.irc_max_points = _int_option(option)
            elif key == "STEPSIZE" and _int_option(option):
                step.irc_step_size = _int_option(option)
            else:
                rest.append(option)
        return rest

    def _stable(self, keyword: RouteKeyword) -> List[str]:
        self.jobs.append("STABLE")
        rest = [o for o in keyword.options if o.upper() != "OPT"]
        self.step.stable_opt = len(rest) < len(keyword.options)
        return rest

    def _nmr(self, keyword: RouteKeyword) -> List[str]:
        self.jobs.append("NMR")
        methods, rest = gk.reverse(gk.NMR_METHODS), []
        for option in keyword.options:
            key = option.upper()
            if key == "GIAO":
                self.step.nmr_method = "GIAO"
            elif key in methods:
                self.step.nmr_method = methods[key]
            elif key == "SPINSPIN":
                self.step.nmr_spinspin = True
            else:
                rest.append(option)
        return rest

    def _td(self, keyword: RouteKeyword) -> List[str]:
        self.step.state = gk.STATE_TD if keyword.key == "TD" else gk.STATE_TDA
        return self._state_options(keyword.options)

    def _dispersion(self, keyword: RouteKeyword) -> Optional[List[str]]:
        back = gk.reverse(gk.DISPERSION)
        if len(keyword.options) == 1 and keyword.options[0].upper() in back:
            self.step.dispersion = back[keyword.options[0].upper()]
            return []
        return None

    def _iop(self, keyword: RouteKeyword) -> List[str]:
        rest = [o for o in keyword.options if o.replace(" ", "") != gk.COEFFICIENTS_IOP]
        self.step.print_all_coefficients = len(rest) < len(keyword.options)
        return rest

    def _scf(self, keyword: RouteKeyword) -> List[str]:
        rest = []
        for option in keyword.options:
            if option.upper() == "TIGHT":
                self.step.scf_tight = True
            elif option.upper() == "QC":
                self.step.scf_qc = True
            else:
                rest.append(option)
        return rest

    def _switch(self, name: str):
        def handler(keyword: RouteKeyword) -> Optional[List[str]]:
            if keyword.options:
                return None
            setattr(self.step, name, True)
            return []
        return handler

    def _geom(self, keyword: RouteKeyword) -> List[str]:
        sources, rest = gk.reverse(gk.GEOMETRY_SOURCES), []
        for option in keyword.options:
            key = option.upper()
            if key in sources:
                self.step.geometry_source = sources[key]
            elif key in ("CONNECTIVITY", "CONNECT"):
                self.step.connectivity = True
            else:
                rest.append(option)
        return rest

    def _guess(self, keyword: RouteKeyword) -> List[str]:
        methods = gk.reverse(gk.GUESS_METHODS)
        flags = {option.upper(): name for name, option in gk.GUESS_FLAG_FIELDS.items()}
        rest = []
        for option in keyword.options:
            key = option.upper()
            if key in methods:
                self.step.guess = methods[key]
            elif key in flags:
                setattr(self.step, flags[key], True)
            else:
                rest.append(option)
        return rest

    def _pop(self, keyword: RouteKeyword) -> List[str]:
        levels = gk.reverse(gk.POP_LEVELS)
        flags = {option.upper(): name for name, option in gk.POP_FLAG_FIELDS.items()}
        rest = []
        for option in keyword.options:
            key = option.upper()
            if key in levels:
                self.step.pop_level = levels[key]
            elif key in flags:
                setattr(self.step, flags[key], True)
            else:
                rest.append(option)
        return rest

    def _scrf(self, keyword: RouteKeyword) -> List[str]:
        models = gk.reverse(gk.SOLVATION_MODELS)
        solvents = {s.upper(): s for s in gk.SOLVENTS if s != gk.OTHER}
        step, rest = self.step, []
        step.solvation, step.solvent = "IEFPCM", "Water"   # Gaussian's defaults for a bare SCRF
        for option in keyword.options:
            key = gk.SOLVATION_ALIASES.get(option.upper(), option.upper())
            if key in models:
                step.solvation = models[key]
            elif option_key(option) == "SOLVENT" and "=" in option:
                name = option.split("=", 1)[1].strip()
                if name.upper() in solvents:
                    step.solvent = solvents[name.upper()]
                else:
                    step.solvent, step.solvent_other = gk.OTHER, name
            else:
                rest.append(option)
        return rest


def _split_blocks(text: str) -> List[List[str]]:
    blocks, current = [], []
    for raw in (text or "").splitlines():
        if raw.strip().startswith("!"):
            continue                                  # a whole-line comment is not a blank line
        line = raw.split("!", 1)[0].rstrip()          # '!' starts a comment in Gaussian input
        if line.strip().lower() == "--link1--":
            blocks.append(current)
            current = []
        else:
            current.append(line)
    blocks.append(current)
    return [b for b in blocks if any(line.strip() for line in b)]


def _sections(lines: Sequence[str]) -> List[List[str]]:
    sections, current = [], []
    for line in lines:
        if line.strip():
            current.append(line)
        elif current:
            sections.append(current)
            current = []
    if current:
        sections.append(current)
    return sections


def _reopened_name(value: str, label: str, where: str, notes: List[str]) -> str:
    name = PureWindowsPath(value.strip()).name
    if name != value.strip():
        notes.append(f"{where}: {label}={value.strip()} becomes {name} (the job runs from the input's folder).")
    return name


def _read_charge(step: GaussianStepSpec, line: str, where: str) -> None:
    tokens = line.replace(",", " ").split()
    if len(tokens) != 2 or not all(_INTEGER_RE.fullmatch(t) for t in tokens):
        raise ValueError(f"{where}: '{line.strip()}' is not one charge and one multiplicity "
                         "(per-fragment charges are not supported).")
    step.charge, step.multiplicity = int(tokens[0]), int(tokens[1])


def _read_block(lines: Sequence[str], index: int, stem: str, notes: List[str]
                ) -> Tuple[GaussianStepSpec, Optional[str]]:
    """One linked step -> (its panels, its geometry text or None)."""
    where = f"Step {index}"
    sections = _sections(lines)
    head = sections[0]
    route_start = next((k for k, line in enumerate(head) if line.lstrip().startswith("#")), None)
    if route_start is None:
        raise ValueError(f"{where}: no route section (a line starting with '#') was found.")
    step = GaussianStepSpec(write_chk=False)

    extra_link0, seen = [], set()
    for line in head[:route_start]:
        text = line.strip()
        if not text.startswith("%"):
            raise ValueError(f"{where}: '{text}' before the route is not a Link 0 (%) line.")
        key, _, value = text[1:].partition("=")
        key, value = key.strip().lower(), value.strip()
        seen.add(key)
        if key == "mem":
            step.memory = value
        elif key in ("nprocshared", "nproc") and value.isdigit():
            step.nproc = int(value)
        elif key == "chk":
            step.write_chk = True
            name = _reopened_name(value, "%Chk", where, notes)
            step.chk = "" if name == f"{stem}.chk" else name
        elif key == "oldchk":
            step.old_chk = _reopened_name(value, "%OldChk", where, notes)
        else:
            extra_link0.append(text)
    step.extra_link0 = "\n".join(extra_link0)
    if "mem" not in seen:
        notes.append(f"{where}: no %Mem in the file; {fi.DEFAULT_MEMORY} is used.")

    route_text = " ".join(line.strip() for line in head[route_start:])
    match = re.match(r"#([PNTpnt]?)(?=\s|$)", route_text)
    step.print_level = "#" + (match.group(1).upper() if match else "")
    reader = _RouteReader(step, notes, where)
    reader.read(fi.route_tokens(route_text[match.end():] if match else route_text[1:]))

    rest = sections[1:]
    geometry = None
    title_lines: List[str] = []
    if step.geometry_source == gk.GEOM_CHECK:
        if len(rest) < 2:
            raise ValueError(f"{where}: Geom=Check needs a title and a charge/multiplicity line.")
        title_lines, charge_lines, rest = rest[0], rest[1], rest[2:]
        _read_charge(step, charge_lines[0], where)
    elif step.geometry_source == gk.GEOM_COORDINATES:
        if len(rest) < 2:
            raise ValueError(f"{where}: the file ends before the molecule specification.")
        title_lines, molecule, rest = rest[0], rest[1], rest[2:]
        _read_charge(step, molecule[0], where)
        try:
            atoms = fi.parse_geometry("\n".join(molecule[1:]), check_distances=False)
        except ValueError as exc:
            raise ValueError(f"{where}: {exc}") from exc
        geometry = "\n".join(fi.format_geometry(atoms))
    if title_lines:
        title = " ".join(line.strip() for line in title_lines)
        clean = " ".join("".join(" " if ch in fi.TITLE_FORBIDDEN else ch for ch in title).split())
        if clean != title:
            notes.append(f"{where}: the title had characters Gaussian does not allow; they became blanks.")
        step.title = "" if clean == fi._title_text(stem) else clean

    if step.connectivity:
        if rest and rest[0][0].split()[0].isdigit():
            rest = rest[1:]
            notes.append(f"{where}: the connectivity section is rebuilt from the distances when saving.")
        else:
            notes.append(f"{where}: Geom=Connectivity but no connectivity section follows.")
    scan_lines = False
    if reader.modredundant:
        if rest:
            step.modredundant = "\n".join(line.strip() for line in rest[0])
            rest = rest[1:]
            try:
                scan_lines = any(p[2] == "S" for p in parse_modredundant(step.modredundant))
            except ValueError:
                scan_lines = False      # prepare reports the line
        else:
            notes.append(f"{where}: Opt=ModRedundant but no coordinate lines follow.")
    step.job_type = reader.job_type(scan_lines)
    step.additional_input = "\n\n".join("\n".join(section) for section in rest)
    step.additional_keywords = " ".join(reader.kept_extra())
    return step, geometry


def read_input_text(text: str, stem: str = fi.DEFAULT_STEM, output_dir: str = ""
                    ) -> Tuple[GaussianJobSpec, List[str]]:
    """A Gaussian input -> (job spec, notes). ValueError names the step."""
    blocks = _split_blocks(text)
    if not blocks:
        raise ValueError("The input is empty.")
    notes: List[str] = []
    steps: List[GaussianStepSpec] = []
    geometry = ""
    for index, lines in enumerate(blocks, start=1):
        step, step_geometry = _read_block(lines, index, stem, notes)
        if step_geometry is not None:
            if not geometry:
                geometry = step_geometry
            elif step_geometry != geometry:
                raise ValueError(f"Step {index} has its own geometry, different from the first one; the builder "
                                 "writes one geometry per job.")
        steps.append(step)
    if len(steps) == 3 and " Fukui " in f" {steps[0].title} ":
        notes.append("This is a Fukui job written by the builder; it is shown as three linked steps. Choose the "
                     "Fukui functions job type to rebuild it from its panels.")
    return GaussianJobSpec(stem=stem, output_dir=output_dir, geometry=geometry, steps=steps), notes


def read_input_file(path) -> Tuple[GaussianJobSpec, List[str]]:
    path = Path(path)
    if path.suffix.lower() not in GAUSSIAN_INPUT_SUFFIXES:
        raise ValueError(f"'{path.name}' is not a Gaussian input ({', '.join(GAUSSIAN_INPUT_SUFFIXES)}).")
    try:
        text = path.read_bytes().decode("utf-8-sig", errors="replace")
    except OSError as exc:
        raise ValueError(f"'{path.name}' could not be opened: {exc.strerror or exc}") from exc
    stem = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in path.stem) or fi.DEFAULT_STEM
    try:
        return read_input_text(text, stem, str(path.resolve().parent))
    except ValueError as exc:
        raise ValueError(f"{path.name}: {exc}") from exc
