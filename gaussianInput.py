"""
Gaussian 16 input builder: the job model, validation, route, input text and reopening (no Toga).

``GaussianJobSpec`` (job-wide fields + one ``GaussianStepSpec`` per linked step) ->
``GaussianInputBuilder.prepare`` (every problem in one ValueError, one '- Step n: ...' line each) ->
``GaussianJob`` -> ``input_text`` / ``script_text`` / ``write_files``. ``read_input_text`` is the
inverse used by Open input. The Fukui Job Type is delegated to ``fukuiInput`` unchanged. Files are
ASCII with LF endings and relative checkpoint names: the job runs under Linux from the input's folder.
"""

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import fukuiInput as fi
import gaussianKeywords as gk

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
