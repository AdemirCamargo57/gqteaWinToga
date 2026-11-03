"""
Data tables of the Gaussian 16 input builder: what every control writes into the route.

Modelled on GaussView 6.1's Gaussian Calculation Setup (its help page gausscalcopts.htm lists each
control's keyword). Keys are the labels the window shows; values are Gaussian keywords or options
(None = write nothing). ``gaussianInput`` builds routes from these tables and reads routes back
through ``reverse``, so a label and its keyword are defined only here. No Toga.
"""

from typing import Dict, Mapping, Optional

from fukuiInput import DEFAULT_PRINT_LEVEL, PRINT_LEVELS

OTHER = "Other..."
NO_CHOICE = "(none)"          # the method list of a family that has no sub-choice (HF, MP2, ...)

# --- Job Type -----------------------------------------------------------------
JOB_ENERGY = "Energy"
JOB_OPT = "Optimization"
JOB_FREQ = "Frequency"
JOB_OPT_FREQ = "Opt+Freq"
JOB_IRC = "IRC"
JOB_SCAN = "Scan (relaxed)"
JOB_STABILITY = "Stability"
JOB_NMR = "NMR"
JOB_FUKUI = "Fukui functions (N, N+p, N-q)"
JOB_TYPES = (JOB_ENERGY, JOB_OPT, JOB_FREQ, JOB_OPT_FREQ, JOB_IRC, JOB_SCAN, JOB_STABILITY, JOB_NMR,
             JOB_FUKUI)
OPTIMIZING_JOBS = (JOB_OPT, JOB_OPT_FREQ, JOB_SCAN)
FREQUENCY_JOBS = (JOB_FREQ, JOB_OPT_FREQ)

OPT_MINIMUM = "Minimum"
OPT_TS = "TS (Berny)"
OPT_TARGETS: Dict[str, Optional[str]] = {OPT_MINIMUM: None, OPT_TS: "TS"}

FC_DEFAULT = "Default"
FC_CALCFC = "Calculate at first point (CalcFC)"
FORCE_CONSTANTS: Dict[str, Optional[str]] = {
    FC_DEFAULT: None,
    FC_CALCFC: "CalcFC",
    "Calculate at all points (CalcAll)": "CalcAll",
    "Read from checkpoint (ReadFC)": "ReadFC",
}
FUKUI_FORCE_CONSTANTS = (FC_DEFAULT, FC_CALCFC, "Calculate at all points (CalcAll)")
READ_FORCE_CONSTANTS = ("ReadFC", "RCFC")          # options that read the checkpoint

IRC_DIRECTIONS: Dict[str, Optional[str]] = {"Both": None, "Forward only": "Forward", "Reverse only": "Reverse"}
NMR_METHODS: Dict[str, Optional[str]] = {"GIAO": None, "CSGT": "CSGT", "IGAIM": "IGAIM", "All": "All"}
FREQ_FLAG_FIELDS = {"freq_raman": "Raman", "freq_vcd": "VCD", "freq_anharmonic": "Anharmonic"}

# --- Method -------------------------------------------------------------------
STATE_GROUND = "Ground state"
STATE_TD = "TD-SCF"
STATE_TDA = "TDA"
STATE_CIS = "CIS"
STATE_ZINDO = "ZINDO"
STATES = (STATE_GROUND, STATE_TD, STATE_TDA, STATE_CIS, STATE_ZINDO)
STATE_KEYWORDS = {STATE_TD: "TD", STATE_TDA: "TDA", STATE_CIS: "CIS", STATE_ZINDO: "ZINDO"}
MODEL_STATES = (STATE_CIS, STATE_ZINDO)            # written as the method of the model token: CIS(...)/basis
STATE_SPINS: Dict[str, Optional[str]] = {"Default": None, "Singlets": "Singlets", "Triplets": "Triplets",
                                         "50-50": "50-50"}
DEFAULT_N_STATES = 3                               # Gaussian's own default
COEFFICIENTS_IOP = "9/40=4"                        # print every CI coefficient (the NTO/Excitation tabs need it)

SPIN_PREFIXES: Dict[str, str] = {"Default": "", "Restricted (R)": "R", "Unrestricted (U)": "U",
                                 "Restricted-open (RO)": "RO"}

FAMILY_HF = "Hartree-Fock"
FAMILY_DFT = "DFT"
FAMILY_SEMI = "Semi-empirical"
FAMILY_MP2 = "MP2"
FAMILY_MP4 = "MP4(SDTQ)"
FAMILY_CCSD = "CCSD"
FAMILY_CCSDT = "CCSD(T)"
FAMILY_COMPOUND = "Compound"
FAMILIES = (FAMILY_HF, FAMILY_DFT, FAMILY_SEMI, FAMILY_MP2, FAMILY_MP4, FAMILY_CCSD, FAMILY_CCSDT,
            FAMILY_COMPOUND)
FAMILY_KEYWORDS = {FAMILY_HF: "HF", FAMILY_MP2: "MP2", FAMILY_MP4: "MP4(SDTQ)", FAMILY_CCSD: "CCSD",
                   FAMILY_CCSDT: "CCSD(T)"}       # families without a sub-choice
FAMILY_METHODS: Dict[str, Dict[str, Optional[str]]] = {
    FAMILY_DFT: {"B3LYP": "B3LYP", "CAM-B3LYP": "CAM-B3LYP", "wB97X-D": "wB97XD", "M06-2X": "M062X",
                 "M06": "M06", "PBE0 (PBE1PBE)": "PBE1PBE", "PBEPBE": "PBEPBE", "B3PW91": "B3PW91",
                 "BP86": "BP86", "TPSSTPSS": "TPSSTPSS", "APFD": "APFD", "HSEH1PBE": "HSEH1PBE",
                 "LC-wHPBE": "LC-wHPBE", OTHER: None},
    FAMILY_SEMI: {"PM7": "PM7", "PM6": "PM6", "AM1": "AM1", "PM3": "PM3", "PDDG": "PDDG"},
    FAMILY_COMPOUND: {"CBS-QB3": "CBS-QB3", "CBS-4M": "CBS-4M", "G4": "G4", "G4MP2": "G4MP2", "W1BD": "W1BD"},
}
NO_BASIS_FAMILIES = (FAMILY_SEMI, FAMILY_COMPOUND)
POST_SCF_FAMILIES = (FAMILY_MP2, FAMILY_MP4, FAMILY_CCSD, FAMILY_CCSDT)
DEFAULT_METHOD = {FAMILY_DFT: "B3LYP", FAMILY_SEMI: "PM7", FAMILY_COMPOUND: "CBS-QB3"}

POPLE_BASES = ("STO-3G", "3-21G", "6-31G", "6-311G")
EXTENSIBLE_POPLE = ("6-31G", "6-311G")             # the ones the diffuse/polarization lists apply to
OTHER_BASES = ("cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ",
               "Def2SVP", "Def2TZVP", "Def2TZVPP", "Def2QZVP", "LanL2DZ", "SDD", "DGDZVP", "Gen", "GenECP")
BASES = POPLE_BASES + OTHER_BASES + (OTHER,)
GENERAL_BASES = ("Gen", "GenECP")
DIFFUSE: Dict[str, str] = {"None": "", "+": "+", "++": "++"}
POLARIZATION: Dict[str, str] = {"None": "", "(d)": "(d)", "(d,p)": "(d,p)", "(2d,p)": "(2d,p)",
                                "(2df,2pd)": "(2df,2pd)", "(3df,3pd)": "(3df,3pd)"}
DISPERSION: Dict[str, Optional[str]] = {"None": None, "GD3": "GD3", "GD3BJ": "GD3BJ", "GD2": "GD2"}

# --- Molecule -----------------------------------------------------------------
GEOM_COORDINATES = "Coordinates (the step-1 geometry)"
GEOM_ALLCHECK = "Checkpoint: everything (Geom=AllCheck)"
GEOM_CHECK = "Checkpoint: new charge and multiplicity (Geom=Check)"
GEOMETRY_SOURCES: Dict[str, Optional[str]] = {GEOM_COORDINATES: None, GEOM_ALLCHECK: "AllCheck",
                                              GEOM_CHECK: "Check"}
MODRED_TYPES = {"Bond": ("B", 2), "Angle": ("A", 3), "Dihedral": ("D", 4)}
MODRED_ACTIONS = {"Scan": "S", "Freeze": "F"}

# --- Guess / Pop / Solvation --------------------------------------------------
GUESS_METHODS: Dict[str, Optional[str]] = {"Default": None, "Core Hamiltonian": "Core",
                                           "Extended Huckel": "Huckel", "Read checkpoint": "Read",
                                           "Read if possible (TCheck)": "TCheck"}
READ_GUESSES = ("Read",)                           # TCheck falls back to a new guess, so it needs no checkpoint
GUESS_FLAG_FIELDS = {"guess_mix": "Mix", "guess_always": "Always", "guess_only": "Only", "guess_save": "Save"}

POP_LEVELS: Dict[str, Optional[str]] = {"Default": None, "Full": "Full", "Regular": "Regular",
                                        "Minimal": "Minimal", "None": "None"}
POP_FLAG_FIELDS = {"pop_mk": "MK", "pop_hirshfeld": "Hirshfeld", "pop_nbo": "NBO", "pop_npa": "NPA",
                   "pop_savento": "SaveNTO"}

SOLVATION_NONE = "None (gas phase)"
SOLVATION_MODELS: Dict[str, Optional[str]] = {SOLVATION_NONE: None, "IEFPCM": "IEFPCM", "CPCM": "CPCM",
                                              "SMD": "SMD"}
SOLVATION_ALIASES = {"PCM": "IEFPCM"}              # read back: SCRF=PCM is IEFPCM in Gaussian 16
SOLVENTS = ("Water", "Acetonitrile", "Methanol", "Ethanol", "DiMethylSulfoxide", "Acetone", "Chloroform",
            "DiChloroMethane", "TetraHydroFuran", "Toluene", "Benzene", "n-Hexane", "CarbonTetraChloride",
            "DiethylEther", OTHER)

# Route keyword (upper case) -> the panel that writes it, for messages.
PANEL_OF_KEYWORD = {
    "OPT": "Job Type", "FREQ": "Job Type", "IRC": "Job Type", "STABLE": "Job Type", "NMR": "Job Type",
    "TD": "Method", "TDA": "Method", "EMPIRICALDISPERSION": "Method", "IOP": "Method",
    "SCF": "General", "NOSYMM": "General", "POLAR": "General", "GEOM": "Molecule/General",
    "GUESS": "Guess", "POP": "Pop", "SCRF": "Solvation",
}


def reverse(table: Mapping[str, Optional[str]]) -> Dict[str, str]:
    """Route value (upper case) -> label; entries that write nothing are skipped."""
    return {value.upper(): label for label, value in table.items() if value}
