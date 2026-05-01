from __future__ import annotations

import base64
import html
import http.server
import re
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW


# -----------------------------------------------------------------------------
# Metadata layer
# -----------------------------------------------------------------------------

@dataclass
class KeywordSpec:
    name: str
    kind: str = "str"  # str, int, float, bool, choice
    default: Any = ""
    required: bool = False
    choices: List[str] = field(default_factory=list)
    help_text: str = ""


@dataclass
class CardSpec:
    name: str
    unit_options: List[str] = field(default_factory=list)
    help_text: str = ""


class CPInputBuilder:
    """Logic class for cp.x input generation."""

    BOHR_PER_ANG = 1.8897261246257702

    ATOMIC_MASSES = {
        "H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182, "B": 10.811,
        "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984032, "Ne": 20.1797,
        "Na": 22.98976928, "Mg": 24.3050, "Al": 26.9815386, "Si": 28.0855,
        "P": 30.973762, "S": 32.065, "Cl": 35.453, "Ar": 39.948, "K": 39.0983,
        "Ca": 40.078, "Sc": 44.955912, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961,
        "Mn": 54.938045, "Fe": 55.845, "Co": 58.933195, "Ni": 58.6934, "Cu": 63.546,
        "Zn": 65.38, "Ga": 69.723, "Ge": 72.64, "As": 74.92160, "Se": 78.96,
        "Br": 79.904, "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62, "Y": 88.90585,
        "Zr": 91.224, "Nb": 92.90638, "Mo": 95.96, "Tc": 98.0, "Ru": 101.07,
        "Rh": 102.90550, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.411, "In": 114.818,
        "Sn": 118.710, "Sb": 121.760, "Te": 127.60, "I": 126.90447, "Xe": 131.293,
        "Cs": 132.9054519, "Ba": 137.327, "La": 138.90547, "Ce": 140.116,
        "Pr": 140.90765, "Nd": 144.242, "Sm": 150.36, "Eu": 151.964, "Gd": 157.25,
        "Tb": 158.92535, "Dy": 162.500, "Ho": 164.93032, "Er": 167.259,
        "Tm": 168.93421, "Yb": 173.054, "Lu": 174.9668, "Hf": 178.49, "Ta": 180.94788,
        "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.084,
        "Au": 196.966569, "Hg": 200.59, "Tl": 204.3833, "Pb": 207.2, "Bi": 208.98040,
    }

    SECTION_SPECS: Dict[str, List[KeywordSpec]] = {
        "CONTROL": [
            KeywordSpec("calculation", "choice", "cp", True, ["cp", "scf", "nscf", "relax", "vc-relax", "vc-cp", "cp-wf", "vc-cp-wf"]),
            KeywordSpec("title", "str", "MD Simulation "),
            KeywordSpec("verbosity", "choice", "low", False, ["debug", "high", "medium", "low", "default", "minimal"]),
            KeywordSpec("isave", "int", 100),
            KeywordSpec("restart_mode", "choice", "restart", False, ["from_scratch", "restart", "reset_counters"]),
            KeywordSpec("nstep", "int", 50),
            KeywordSpec("iprint", "int", 10),
            KeywordSpec("tstress", "bool", False),
            KeywordSpec("tprnfor", "bool", False),
            KeywordSpec("dt", "float", "1.D0"),
            KeywordSpec("outdir", "str", "./"),
            KeywordSpec("saverho", "bool", False),
            KeywordSpec("prefix", "str", "cp"),
            KeywordSpec("ndr", "int", 50),
            KeywordSpec("ndw", "int", 50),
            KeywordSpec("tabps", "bool", False),
            KeywordSpec("max_seconds", "float", "1.D+7"),
            KeywordSpec("etot_conv_thr", "float", "1.0D-4"),
            KeywordSpec("forc_conv_thr", "float", "1.0D-3"),
            KeywordSpec("ekin_conv_thr", "float", "1.0D-6"),
            KeywordSpec("disk_io", "choice", "default", False, ["default", "high"]),
            KeywordSpec("memory", "choice", "default", False, ["default", "small"]),
            KeywordSpec("pseudo_dir", "str", "$HOME/espresso/pseudo/"),
            KeywordSpec("tefield", "bool", False),
        ],
        "SYSTEM": [
            KeywordSpec("ibrav", "int", 0, True),
            KeywordSpec("celldm(1)", "float", ""),
            KeywordSpec("celldm(2)", "float", ""),
            KeywordSpec("celldm(3)", "float", ""),
            KeywordSpec("celldm(4)", "float", ""),
            KeywordSpec("celldm(5)", "float", ""),
            KeywordSpec("celldm(6)", "float", ""),
            KeywordSpec("A", "float", ""),
            KeywordSpec("B", "float", ""),
            KeywordSpec("C", "float", ""),
            KeywordSpec("cosAB", "float", ""),
            KeywordSpec("cosAC", "float", ""),
            KeywordSpec("cosBC", "float", ""),
            KeywordSpec("nat", "int", 0, True),
            KeywordSpec("ntyp", "int", 0, True),
            KeywordSpec("nbnd", "int", ""),
            KeywordSpec("tot_charge", "float", 0.0),
            KeywordSpec("tot_magnetization", "float", ""),
            KeywordSpec("ecutwfc", "float", "", True),
            KeywordSpec("ecutrho", "float", ""),
            KeywordSpec("nr1", "int", ""), KeywordSpec("nr2", "int", ""), KeywordSpec("nr3", "int", ""),
            KeywordSpec("nr1s", "int", ""), KeywordSpec("nr2s", "int", ""), KeywordSpec("nr3s", "int", ""),
            KeywordSpec("nr1b", "int", ""), KeywordSpec("nr2b", "int", ""), KeywordSpec("nr3b", "int", ""),
            KeywordSpec("occupations", "choice", "fixed", False, ["fixed", "ensemble"]),
            KeywordSpec("degauss", "float", "0.D0"),
            KeywordSpec("smearing", "choice", "gaussian", False, ["gaussian", "fermi-dirac", "hermite-delta", "gaussian-splines", "cold-smearing", "marzari-vanderbilt", "0", "-1"]),
            KeywordSpec("nspin", "choice", 1, False, ["1", "2"]),
            KeywordSpec("ecfixed", "float", 0.0),
            KeywordSpec("qcutz", "float", 0.0),
            KeywordSpec("q2sigma", "float", 0.1),
            KeywordSpec("input_dft", "str", ""),
            KeywordSpec("exx_fraction", "float", ""),
            KeywordSpec("lda_plus_u", "bool", False),
            KeywordSpec("vdw_corr", "choice", "none", False, ["none", "grimme-d2", "DFT-D", "TS", "ts-vdw", "tkatchenko-scheffler", "XDM", "xdm"]),
            KeywordSpec("london_s6", "float", 0.75),
            KeywordSpec("london_rcut", "float", 200),
            KeywordSpec("ts_vdw", "bool", False),
            KeywordSpec("ts_vdw_econv_thr", "float", "1.D-6"),
            KeywordSpec("ts_vdw_isolated", "bool", False),
            KeywordSpec("assume_isolated", "choice", "none", False, ["none", "makov-payne", "m-p", "mp"]),
        ],
        "ELECTRONS": [
            KeywordSpec("electron_maxstep", "int", 100),
            KeywordSpec("electron_dynamics", "choice", "none", False, ["none", "sd", "damp", "verlet", "cg"]),
            KeywordSpec("conv_thr", "float", "1.D-6"),
            KeywordSpec("niter_cg_restart", "int", 20),
            KeywordSpec("efield", "float", "0.D0"),
            KeywordSpec("epol", "choice", 3, False, ["1", "2", "3"]),
            KeywordSpec("emass", "float", "400.D0"),
            KeywordSpec("emass_cutoff", "float", "2.5D0"),
            KeywordSpec("orthogonalization", "choice", "ortho", False, ["ortho", "Gram-Schmidt"]),
            KeywordSpec("ortho_eps", "float", "1.D-8"),
            KeywordSpec("ortho_max", "int", 300),
            KeywordSpec("ortho_para", "int", 0),
            KeywordSpec("electron_damping", "float", "0.1D0"),
            KeywordSpec("electron_velocities", "choice", "default", False, ["zero", "default", "change_step"]),
            KeywordSpec("electron_temperature", "choice", "not_controlled", False, ["nose", "rescaling", "not_controlled"]),
            KeywordSpec("ekincw", "float", "0.001D0"),
            KeywordSpec("fnosee", "float", "1.D0"),
            KeywordSpec("startingwfc", "choice", "random", False, ["atomic", "random"]),
            KeywordSpec("tcg", "bool", False),
            KeywordSpec("maxiter", "int", 100),
            KeywordSpec("passop", "float", "0.3D0"),
            KeywordSpec("pre_state", "bool", False),
            KeywordSpec("n_inner", "int", 2),
            KeywordSpec("niter_cold_restart", "int", 1),
            KeywordSpec("lambda_cold", "float", "0.03D0"),
            KeywordSpec("grease", "float", "1.D0"),
            KeywordSpec("ampre", "float", "0.D0"),
        ],
        "IONS": [
            KeywordSpec("ion_dynamics", "choice", "none", False, ["none", "sd", "damp", "verlet", "bfgs"]),
            KeywordSpec("ion_positions", "choice", "default", False, ["default", "from_input"]),
            KeywordSpec("ion_velocities", "choice", "default", False, ["default", "zero", "random", "from_input", "change_step"]),
            KeywordSpec("ion_damping", "float", "0.2D0"),
            KeywordSpec("ion_radius", "float", "1.D0"),
            KeywordSpec("iesr", "int", 0),
            KeywordSpec("ion_nstepe", "int", 1),
            KeywordSpec("remove_rigid_rot", "bool", False),
            KeywordSpec("ion_temperature", "choice", "not_controlled", False, ["not_controlled", "nose", "rescaling"]),
            KeywordSpec("tempw", "float", 300.0),
            KeywordSpec("fnosep", "float", "1.D0"),
            KeywordSpec("tolp", "float", ""),
            KeywordSpec("nhpcl", "int", 0),
            KeywordSpec("nhptyp", "int", 0),
            KeywordSpec("nhgrp", "int", 0),
            KeywordSpec("fnhscl", "float", "1.D0"),
            KeywordSpec("ndega", "int", -1),
            KeywordSpec("tranp", "bool", False),
            KeywordSpec("amprp", "float", "0.D0"),
            KeywordSpec("greasp", "float", "1.D0"),
        ],
        "CELL": [
            KeywordSpec("cell_parameters", "choice", "default", False, ["default", "from_input"]),
            KeywordSpec("cell_dynamics", "choice", "none", False, ["none", "sd", "damp", "pr", "w"]),
            KeywordSpec("cell_velocities", "choice", "default", False, ["default", "zero", "from_input"]),
            KeywordSpec("cell_damping", "float", "2.D0"),
            KeywordSpec("press", "float", 0.0),
            KeywordSpec("wmass", "float", ""),
            KeywordSpec("cell_factor", "float", 1.2),
            KeywordSpec("cell_temperature", "choice", "not_controlled", False, ["not_controlled", "nose", "rescaling"]),
            KeywordSpec("temph", "float", 300.0),
            KeywordSpec("fnoseh", "float", "1.D0"),
            KeywordSpec("greash", "float", "1.D0"),
            KeywordSpec("cell_dofree", "choice", "all", False, ["all", "shape", "volume", "2Dxy", "2Dshape", "x", "y", "z", "xy", "xz", "yz"]),
        ],
        "WANNIER": [
            KeywordSpec("wf_efield", "int", 0),
            KeywordSpec("wf_switch", "bool", False),
            KeywordSpec("sw_len", "int", 1),
            KeywordSpec("efx0", "float", 0.0), KeywordSpec("efy0", "float", 0.0), KeywordSpec("efz0", "float", 0.0),
            KeywordSpec("efx1", "float", 0.0), KeywordSpec("efy1", "float", 0.0), KeywordSpec("efz1", "float", 0.0),
            KeywordSpec("wfsd", "int", 1),
            KeywordSpec("wfdt", "float", 1.0),
            KeywordSpec("maxwfdt", "float", 1.0),
            KeywordSpec("nit", "int", 10),
            KeywordSpec("nsd", "int", 10),
            KeywordSpec("wf_q", "float", 300.0),
            KeywordSpec("wf_friction", "float", 0.3),
            KeywordSpec("nsteps", "int", 100),
            KeywordSpec("tolw", "float", "1.D-6"),
            KeywordSpec("adapt", "bool", False),
            KeywordSpec("calwf", "int", 0),
            KeywordSpec("nwf", "int", 0),
            KeywordSpec("wffort", "int", 1),
            KeywordSpec("writev", "bool", False),
            KeywordSpec("exx_neigh", "int", 60),
            KeywordSpec("exx_dis_cutoff", "float", 8.0),
            KeywordSpec("exx_poisson_eps", "float", "1.0D-6"),
            KeywordSpec("exx_use_cube_domain", "bool", False),
            KeywordSpec("exx_ps_rcut_self", "float", 6.0),
            KeywordSpec("exx_ps_rcut_pair", "float", 5.0),
            KeywordSpec("exx_me_rcut_self", "float", 10.0),
            KeywordSpec("exx_me_rcut_pair", "float", 7.0),
        ],
    }

    CARD_SPECS: Dict[str, CardSpec] = {
        "ATOMIC_SPECIES": CardSpec("ATOMIC_SPECIES"),
        "ATOMIC_POSITIONS": CardSpec("ATOMIC_POSITIONS", ["alat", "bohr", "crystal", "angstrom"]),
        "ATOMIC_VELOCITIES": CardSpec("ATOMIC_VELOCITIES"),
        "CELL_PARAMETERS": CardSpec("CELL_PARAMETERS", ["alat", "bohr", "angstrom"]),
        "REF_CELL_PARAMETERS": CardSpec("REF_CELL_PARAMETERS", ["alat", "bohr", "angstrom"]),
        "CONSTRAINTS": CardSpec("CONSTRAINTS"),
        "OCCUPATIONS": CardSpec("OCCUPATIONS"),
        "ATOMIC_FORCES": CardSpec("ATOMIC_FORCES"),
        "PLOT_WANNIER": CardSpec("PLOT_WANNIER"),
        "AUTOPILOT": CardSpec("AUTOPILOT"),
    }

    def parse_xyz(self, filepath: str) -> List[Tuple[str, float, float, float]]:
        path = Path(filepath)
        lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
        if len(lines) < 2:
            raise ValueError("The XYZ file is too short.")
        try:
            nat = int(lines[0])
        except ValueError as exc:
            raise ValueError("The first line of the XYZ file must contain the number of atoms.") from exc
        atom_lines = lines[2:2 + nat]
        if len(atom_lines) != nat:
            raise ValueError("The XYZ file does not contain the expected number of atomic coordinate lines.")
        atoms: List[Tuple[str, float, float, float]] = []
        for idx, line in enumerate(atom_lines, start=1):
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"Invalid XYZ line {idx + 2}: '{line}'")
            symbol = self.normalize_symbol(parts[0])
            x, y, z = map(float, parts[1:4])
            atoms.append((symbol, x, y, z))
        return atoms

    @staticmethod
    def normalize_symbol(raw: str) -> str:
        raw = raw.strip()
        if not raw:
            raise ValueError("Empty atomic symbol in XYZ file.")
        raw = re.sub(r"[^A-Za-z]", "", raw)
        if not raw:
            raise ValueError("Invalid atomic label in XYZ file.")
        return raw[0].upper() + raw[1:].lower()

    def species_from_atoms(self, atoms: List[Tuple[str, float, float, float]]) -> List[str]:
        seen = []
        for symbol, *_ in atoms:
            if symbol not in seen:
                seen.append(symbol)
        return seen

    def default_species_rows(self, atoms: List[Tuple[str, float, float, float]]) -> List[Dict[str, str]]:
        rows = []
        for symbol in self.species_from_atoms(atoms):
            mass = self.ATOMIC_MASSES.get(symbol, 0.0)
            rows.append({
                "symbol": symbol,
                "mass": f"{mass:.6f}" if mass else "",
                "pseudo": f"{symbol}.UPF",
            })
        return rows

    def format_value(self, spec: KeywordSpec, value: Any) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, str):
            stripped = value.strip()
            if stripped == "":
                return None
            value = stripped
        if spec.kind == "bool":
            if isinstance(value, str):
                normalized = value.strip().lower()
                truthy = normalized in {"true", ".true.", "1", "yes", "y"}
            else:
                truthy = bool(value)
            return ".true." if truthy else ".false."
        if spec.kind in {"str", "choice"}:
            text = str(value)
            if text.startswith("'") and text.endswith("'"):
                return text
            # QE character values are usually quoted.
            if spec.kind == "choice" or spec.name not in {"pseudo_dir", "outdir"}:
                return f"'{text}'"
            return f"'{text}'"
        return str(value)

    def build_namelist(self, name: str, values: Dict[str, Any], include_defaults: bool = True) -> str:
        specs = self.SECTION_SPECS[name]
        lines = [f"&{name}"]
        for spec in specs:
            raw = values.get(spec.name, spec.default)
            if raw in (None, "") and spec.required:
                raise ValueError(f"Required keyword '{spec.name}' in section &{name} is missing.")
            if raw in (None, "") and not include_defaults:
                continue
            value = raw if raw not in (None, "") else spec.default
            if value in (None, ""):
                continue
            formatted = self.format_value(spec, value)
            if formatted is None:
                continue
            lines.append(f"   {spec.name} = {formatted},")
        lines.append("/")
        return "\n".join(lines)

    def build_atomic_species_card(self, rows: List[Dict[str, str]]) -> str:
        if not rows:
            raise ValueError("ATOMIC_SPECIES is empty.")
        lines = ["ATOMIC_SPECIES"]
        for row in rows:
            symbol = row.get("symbol", "").strip()
            mass = row.get("mass", "").strip()
            pseudo = row.get("pseudo", "").strip()
            if not symbol or not mass or not pseudo:
                raise ValueError("Each ATOMIC_SPECIES row must contain symbol, mass, and pseudopotential.")
            lines.append(f" {symbol} {mass} {pseudo}")
        return "\n".join(lines)

    def build_atomic_positions_card(self, atoms: List[Tuple[str, float, float, float]], unit: str = "angstrom", if_pos: Optional[List[Tuple[str, str, str]]] = None) -> str:
        lines = [f"ATOMIC_POSITIONS ({unit})"]
        if_pos = if_pos or []
        for i, (symbol, x, y, z) in enumerate(atoms):
            line = f" {symbol} {x: .10f} {y: .10f} {z: .10f}"
            if i < len(if_pos) and any(v.strip() for v in if_pos[i]):
                flags = " ".join(v.strip() or "1" for v in if_pos[i])
                line += f" {flags}"
            lines.append(line)
        return "\n".join(lines)

    def build_atomic_velocities_card(self, rows: List[Dict[str, str]]) -> str:
        if not rows:
            return ""
        lines = ["ATOMIC_VELOCITIES"]
        for row in rows:
            label = row.get("label", "").strip()
            vx = row.get("vx", "").strip()
            vy = row.get("vy", "").strip()
            vz = row.get("vz", "").strip()
            if not label:
                continue
            if not vx or not vy or not vz:
                raise ValueError("Each ATOMIC_VELOCITIES row must contain label, vx, vy, and vz.")
            lines.append(f" {label} {vx} {vy} {vz}")
        return "\n".join(lines)

    def build_vector_card(self, card_name: str, unit: str, rows: List[Dict[str, str]]) -> str:
        if not rows:
            return ""
        lines = [f"{card_name} ({unit})"]
        for row in rows[:3]:
            x = row.get("x", "").strip()
            y = row.get("y", "").strip()
            z = row.get("z", "").strip()
            if not x or not y or not z:
                raise ValueError(f"Each vector in {card_name} must contain x, y, and z.")
            lines.append(f" {x} {y} {z}")
        return "\n".join(lines)

    def build_constraints_card(self, text: str) -> str:
        text = text.strip()
        return f"CONSTRAINTS\n{text}" if text else ""

    def build_occupations_card(self, text: str) -> str:
        text = text.strip()
        return f"OCCUPATIONS\n{text}" if text else ""

    def build_atomic_forces_card(self, rows: List[Dict[str, str]]) -> str:
        if not rows:
            return ""
        lines = ["ATOMIC_FORCES"]
        for row in rows:
            label = row.get("label", "").strip()
            fx = row.get("fx", "").strip()
            fy = row.get("fy", "").strip()
            fz = row.get("fz", "").strip()
            if not label:
                continue
            if not fx or not fy or not fz:
                raise ValueError("Each ATOMIC_FORCES row must contain label, fx, fy, and fz.")
            lines.append(f" {label} {fx} {fy} {fz}")
        return "\n".join(lines)

    def build_plot_wannier_card(self, text: str) -> str:
        text = text.strip()
        return f"PLOT_WANNIER\n{text}" if text else ""

    def build_autopilot_card(self, text: str) -> str:
        text = text.strip()
        if not text:
            return ""
        return f"AUTOPILOT\n{text}\nENDRULES"

    def generate_input_text(self, payload: Dict[str, Any]) -> str:
        blocks = []
        for section in ["CONTROL", "SYSTEM", "ELECTRONS", "IONS", "CELL", "WANNIER"]:
            if section in payload["sections"]:
                blocks.append(self.build_namelist(section, payload["sections"][section], include_defaults=True))

        blocks.append(self.build_atomic_species_card(payload["atomic_species"]))
        blocks.append(self.build_atomic_positions_card(
            payload["atoms"],
            unit=payload.get("atomic_positions_unit", "angstrom"),
            if_pos=payload.get("if_pos_flags", []),
        ))

        optional_cards = [
            self.build_atomic_velocities_card(payload.get("atomic_velocities", [])),
            self.build_vector_card("CELL_PARAMETERS", payload.get("cell_parameters_unit", "angstrom"), payload.get("cell_parameters", [])),
            self.build_vector_card("REF_CELL_PARAMETERS", payload.get("ref_cell_parameters_unit", "angstrom"), payload.get("ref_cell_parameters", [])),
            self.build_constraints_card(payload.get("constraints_text", "")),
            self.build_occupations_card(payload.get("occupations_text", "")),
            self.build_atomic_forces_card(payload.get("atomic_forces", [])),
            self.build_plot_wannier_card(payload.get("plot_wannier_text", "")),
            self.build_autopilot_card(payload.get("autopilot_text", "")),
        ]
        blocks.extend([card for card in optional_cards if card.strip()])
        return "\n\n".join(blocks) + "\n"

    def save_input(self, xyz_path: str, input_text: str, filename: Optional[str] = None) -> Path:
        xyz = Path(xyz_path)
        output_name = filename.strip() if filename else f"{xyz.stem}_cp.in"
        out_path = xyz.with_name(output_name)
        out_path.write_text(input_text, encoding="utf-8", newline="\n")
        return out_path


class CPInputBuilderUI:
    """Window-based GUI class for cp.x input generation with Toga."""

    HELP_PDF_CANDIDATES = (
        Path(__file__).with_name("cpx_input_description.pdf"),
    )

    def __init__(self, *args, app: Optional[toga.App] = None):
        resolved_app = app if isinstance(app, toga.App) else None
        self.app = resolved_app
        self.formal_name = getattr(resolved_app, "formal_name", "CP x Input Builder")

        self.builder = CPInputBuilder()
        self.xyz_path: Optional[str] = None
        self.atoms: List[Tuple[str, float, float, float]] = []
        self.species_rows: List[Dict[str, toga.Widget]] = []
        self.if_pos_rows: List[Tuple[toga.TextInput, toga.TextInput, toga.TextInput]] = []
        self.velocity_rows: List[Dict[str, toga.TextInput]] = []
        self.force_rows: List[Dict[str, toga.TextInput]] = []
        self.cell_rows: List[Tuple[toga.TextInput, toga.TextInput, toga.TextInput]] = []
        self.ref_cell_rows: List[Tuple[toga.TextInput, toga.TextInput, toga.TextInput]] = []
        self.preview_window: Optional[toga.Window] = None
        self.help_window: Optional[toga.Window] = None
        self.help_pdf_path: Optional[Path] = None
        self.help_webview: Optional[toga.WebView] = None
        self.help_zoom: float = 1.35
        self.help_zoom_label: Optional[toga.Label] = None
        self.help_temp_dir: Optional[tempfile.TemporaryDirectory[str]] = None
        self.help_httpd: Optional[http.server.ThreadingHTTPServer] = None
        self.help_server_thread: Optional[threading.Thread] = None
        self.help_server_port: Optional[int] = None

        if self.app is None:
            self.main_window = toga.Window(title=self.formal_name, size=(750, 700))
        else:
            self.main_window = toga.MainWindow(title=self.formal_name)

        self.section_widgets: Dict[str, Dict[str, toga.Widget]] = {}
        self.section_boxes: Dict[str, toga.Box] = {}
        self.action_buttons: List[toga.Button] = []

        self.xyz_label = toga.Label("No XYZ file loaded.", style=Pack(margin_bottom=8))
        self.output_name_input = toga.TextInput(value="", placeholder="Leave blank to use <xyzname>_cp.in", style=Pack(width=260, margin_left=8))

        top_actions = toga.Box(style=Pack(direction=ROW, margin_bottom=10))
        load_btn = toga.Button("Open XYZ", on_press=self.open_xyz_file, style=Pack(margin_right=8))
        save_btn = toga.Button("Save cp.x input", on_press=self.save_cp_input, style=Pack(margin_right=8))
        preview_btn = toga.Button("Preview input", on_press=self.preview_input, style=Pack(margin_right=8))
        help_btn = toga.Button("Help", on_press=self.open_help_pdf, style=Pack())
        self.action_buttons = [load_btn, save_btn, preview_btn, help_btn]
        top_actions.add(load_btn)
        top_actions.add(save_btn)
        top_actions.add(preview_btn)
        top_actions.add(help_btn)

        file_row = toga.Box(style=Pack(direction=ROW, margin_bottom=10))
        file_row.add(toga.Label("Output filename:", style=Pack(margin_top=6, margin_right=8)))
        file_row.add(self.output_name_input)

        content = toga.Box(style=Pack(direction=COLUMN, margin=10))
        content.add(top_actions)
        content.add(self.xyz_label)
        content.add(file_row)

        self.scroll = toga.ScrollContainer(style=Pack(flex=1))
        scroll_content = toga.Box(style=Pack(direction=COLUMN, margin=5))

        for section in ["CONTROL", "SYSTEM", "ELECTRONS", "IONS", "CELL", "WANNIER"]:
            sec_box = self.make_section_box(section)
            self.section_boxes[section] = sec_box
            scroll_content.add(sec_box)

        scroll_content.add(self.make_atomic_species_box())
        scroll_content.add(self.make_positions_options_box())
        scroll_content.add(self.make_optional_cards_box())

        self.scroll.content = scroll_content
        content.add(self.scroll)

        self.main_window.content = content
        self.main_window.show()
    def make_section_box(self, section_name: str) -> toga.Box:
        outer = toga.Box(style=Pack(direction=COLUMN, margin_bottom=14))
        outer.add(toga.Label(f"&{section_name}", style=Pack(font_weight="bold", margin_bottom=6)))
        widget_map: Dict[str, toga.Widget] = {}

        for spec in self.builder.SECTION_SPECS[section_name]:
            row = toga.Box(style=Pack(direction=ROW, margin_bottom=4))
            row.add(toga.Label(spec.name, style=Pack(width=180, margin_top=6, margin_right=8)))
            widget: toga.Widget
            if spec.kind == "choice":
                items = [str(v) for v in spec.choices]
                widget = toga.Selection(items=items, value=str(spec.default), style=Pack(width=220))
            elif spec.kind == "bool":
                widget = toga.Selection(items=[".false.", ".true."], value=".true." if spec.default else ".false.", style=Pack(width=120))
            else:
                widget = toga.TextInput(value="" if spec.default == "" else str(spec.default), style=Pack(width=220))
                if spec.default == "":
                    widget.placeholder = spec.name
            widget_map[spec.name] = widget
            row.add(widget)
            if spec.required:
                row.add(toga.Label("required", style=Pack(margin_top=6, margin_left=8)))
            outer.add(row)

        self.section_widgets[section_name] = widget_map
        return outer

    def make_atomic_species_box(self) -> toga.Box:
        outer = toga.Box(style=Pack(direction=COLUMN, margin_bottom=14))
        outer.add(toga.Label("ATOMIC_SPECIES", style=Pack(font_weight="bold", margin_bottom=6)))
        self.species_container = toga.Box(style=Pack(direction=COLUMN))
        outer.add(self.species_container)
        return outer

    def rebuild_species_box(self, rows: List[Dict[str, str]]):
        self.species_container.clear()
        self.species_rows.clear()
        header = toga.Box(style=Pack(direction=ROW, margin_bottom=4))
        for text, width in [("Symbol", 80), ("Mass", 140), ("Pseudopotential", 220)]:
            header.add(toga.Label(text, style=Pack(width=width, font_weight="bold")))
        self.species_container.add(header)
        for row_data in rows:
            row = toga.Box(style=Pack(direction=ROW, margin_bottom=4))
            sym = toga.TextInput(value=row_data["symbol"], style=Pack(width=80, margin_right=6))
            mass = toga.TextInput(value=row_data["mass"], style=Pack(width=140, margin_right=6))
            pseudo = toga.TextInput(value=row_data["pseudo"], style=Pack(width=220, margin_right=6))
            row.add(sym)
            row.add(mass)
            row.add(pseudo)
            self.species_rows.append({"symbol": sym, "mass": mass, "pseudo": pseudo})
            self.species_container.add(row)

    def make_positions_options_box(self) -> toga.Box:
        outer = toga.Box(style=Pack(direction=COLUMN, margin_bottom=14))
        outer.add(toga.Label("ATOMIC_POSITIONS", style=Pack(font_weight="bold", margin_bottom=6)))

        unit_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        unit_row.add(toga.Label("Coordinate unit", style=Pack(width=180, margin_top=6, margin_right=8)))
        self.atomic_positions_unit = toga.Selection(items=["alat", "bohr", "crystal", "angstrom"], value="angstrom", style=Pack(width=160))
        unit_row.add(self.atomic_positions_unit)
        outer.add(unit_row)

        outer.add(toga.Label("Coordinates loaded from XYZ", style=Pack(margin_bottom=4)))
        self.atomic_positions_text = toga.MultilineTextInput(style=Pack(height=140, margin_bottom=8))
        self.atomic_positions_text.placeholder = "Load an XYZ file to display the ATOMIC_POSITIONS coordinates here"
        outer.add(self.atomic_positions_text)

        outer.add(toga.Label("Optional if_pos flags (blank = omitted)", style=Pack(margin_bottom=6)))
        self.if_pos_container = toga.Box(style=Pack(direction=COLUMN))
        outer.add(self.if_pos_container)
        return outer

    def rebuild_if_pos_box(self):
        self.if_pos_container.clear()
        self.if_pos_rows.clear()
        if not self.atoms:
            self.if_pos_container.add(toga.Label("Load an XYZ file to populate atomic rows."))
            return
        header = toga.Box(style=Pack(direction=ROW, margin_bottom=4))
        for text, width in [("Atom", 70), ("if_pos(1)", 80), ("if_pos(2)", 80), ("if_pos(3)", 80)]:
            header.add(toga.Label(text, style=Pack(width=width, font_weight="bold")))
        self.if_pos_container.add(header)
        for symbol, x, y, z in self.atoms:
            row = toga.Box(style=Pack(direction=ROW, margin_bottom=4))
            row.add(toga.Label(symbol, style=Pack(width=70, margin_top=6)))
            a = toga.TextInput(value="", style=Pack(width=80, margin_right=4))
            b = toga.TextInput(value="", style=Pack(width=80, margin_right=4))
            c = toga.TextInput(value="", style=Pack(width=80))
            row.add(a)
            row.add(b)
            row.add(c)
            self.if_pos_rows.append((a, b, c))
            self.if_pos_container.add(row)

    def render_atomic_positions_text(self):
        if not self.atoms:
            self.atomic_positions_text.value = ""
            return
        unit = self.atomic_positions_unit.value or "angstrom"
        self.atomic_positions_text.value = self.builder.build_atomic_positions_card(self.atoms, unit=unit)

    def set_action_buttons_enabled(self, enabled: bool):
        for button in self.action_buttons:
            button.enabled = enabled

    @classmethod
    def resolve_help_pdf_path(cls) -> Optional[Path]:
        for candidate in cls.HELP_PDF_CANDIDATES:
            if candidate.exists():
                return candidate
        return None

    @staticmethod
    def build_help_anchor_id(page_number: int, y_value: float, zoom: float) -> str:
        return f"dest-p{page_number + 1}-y{int(round(y_value * zoom))}"

    def build_help_html(self, pdf_path: Path) -> str:
        try:
            import fitz  # PyMuPDF
        except ImportError as exc:
            raise RuntimeError(
                "In-app PDF viewing requires PyMuPDF. Install it with: pip install pymupdf"
            ) from exc

        zoom = self.help_zoom
        matrix = fitz.Matrix(zoom, zoom)
        sections: List[str] = []

        with fitz.open(pdf_path) as doc:
            destination_anchors: Dict[Tuple[int, int], str] = {}
            page_anchor_html: Dict[int, List[str]] = {page_number: [] for page_number in range(len(doc))}

            for page_number, page in enumerate(doc):
                for link in page.get_links():
                    target_page = link.get("page")
                    target_point = link.get("to")
                    if target_page is None or target_point is None:
                        continue
                    key = (target_page, int(round(target_point.y * zoom)))
                    if key in destination_anchors:
                        continue
                    anchor_id = self.build_help_anchor_id(target_page, target_point.y, zoom)
                    destination_anchors[key] = anchor_id
                    page_anchor_html[target_page].append(
                        f'<div id="{anchor_id}" class="dest-anchor" style="top:{target_point.y * zoom:.1f}px;"></div>'
                    )

            for page_number, page in enumerate(doc):
                pixmap = page.get_pixmap(matrix=matrix, alpha=False)
                encoded = base64.b64encode(pixmap.tobytes("png")).decode("ascii")
                page_width = pixmap.width
                page_height = pixmap.height
                page_title = f"Page {page_number + 1}"

                overlays: List[str] = []
                for link in page.get_links():
                    rect = link.get("from")
                    if rect is None:
                        continue

                    href: Optional[str] = None
                    if link.get("kind") == fitz.LINK_URI and link.get("uri"):
                        href = html.escape(link["uri"], quote=True)
                    elif link.get("page") is not None and link.get("to") is not None:
                        target_point = link["to"]
                        key = (link["page"], int(round(target_point.y * zoom)))
                        anchor_id = destination_anchors.get(key)
                        if anchor_id is not None:
                            href = f"#{anchor_id}"

                    if not href:
                        continue

                    overlays.append(
                        (
                            '<a class="page-link" href="{href}" title="Open link" '
                            'style="left:{left:.1f}px; top:{top:.1f}px; width:{width:.1f}px; height:{height:.1f}px;"></a>'
                        ).format(
                            href=href,
                            left=rect.x0 * zoom,
                            top=rect.y0 * zoom,
                            width=max((rect.x1 - rect.x0) * zoom, 8.0),
                            height=max((rect.y1 - rect.y0) * zoom, 8.0),
                        )
                    )

                sections.append(
                    """
                    <section class="page-section">
                      <div class="page-title">{page_title}</div>
                      <div class="page-frame" id="page-{page_index}" style="width:{page_width}px; height:{page_height}px;">
                        {anchors}
                        <img class="page-image" src="data:image/png;base64,{encoded}" alt="{page_title}" width="{page_width}" height="{page_height}">
                        {overlays}
                      </div>
                    </section>
                    """.format(
                        page_title=html.escape(page_title),
                        page_index=page_number + 1,
                        page_width=page_width,
                        page_height=page_height,
                        encoded=encoded,
                        anchors="".join(page_anchor_html.get(page_number, [])),
                        overlays="".join(overlays),
                    )
                )

        zoom_percent = int(round(zoom * 100))
        return f"""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>cp.x input description</title>
  <style>
    :root {{
      color-scheme: light;
      --page-shadow: 0 12px 28px rgba(17, 24, 39, 0.16);
      --page-border: #d8dde6;
      --page-bg: #ffffff;
      --canvas-bg: #eef2f7;
      --ink-soft: #475569;
      --link-hover: rgba(37, 99, 235, 0.18);
      --link-outline: rgba(37, 99, 235, 0.55);
    }}
    * {{
      box-sizing: border-box;
    }}
    html {{
      scroll-behavior: smooth;
    }}
    body {{
      margin: 0;
      padding: 18px 0 28px;
      background: linear-gradient(180deg, #f8fafc 0%, var(--canvas-bg) 100%);
      font-family: "Segoe UI", Tahoma, sans-serif;
      color: #0f172a;
    }}
    .doc {{
      width: max-content;
      margin: 0 auto;
    }}
    .meta {{
      margin: 0 0 16px;
      text-align: center;
      color: var(--ink-soft);
      font-size: 13px;
    }}
    .page-section {{
      margin: 0 auto 22px;
    }}
    .page-title {{
      margin: 0 0 8px 8px;
      font-size: 13px;
      font-weight: 600;
      color: var(--ink-soft);
      letter-spacing: 0.02em;
      text-transform: uppercase;
    }}
    .page-frame {{
      position: relative;
      background: var(--page-bg);
      border: 1px solid var(--page-border);
      border-radius: 10px;
      overflow: hidden;
      box-shadow: var(--page-shadow);
    }}
    .page-image {{
      display: block;
      user-select: none;
      -webkit-user-drag: none;
    }}
    .page-link {{
      position: absolute;
      display: block;
      border-radius: 3px;
      text-decoration: none;
      background: transparent;
      outline: 1px solid transparent;
      transition: background-color 120ms ease, outline-color 120ms ease;
    }}
    .page-link:hover,
    .page-link:focus {{
      background: var(--link-hover);
      outline-color: var(--link-outline);
    }}
    .dest-anchor {{
      position: absolute;
      left: 0;
      width: 1px;
      height: 1px;
    }}
  </style>
</head>
<body>
  <main class="doc">
    <p class="meta">Zoom: {zoom_percent}% | Hover highlighted areas to follow PDF links.</p>
    {''.join(sections)}
  </main>
</body>
</html>
"""

    async def refresh_help_view(self):
        if self.help_webview is None or self.help_pdf_path is None:
            return
        url = self.write_help_html_file()
        await self.help_webview.load_url(url)
        if self.help_zoom_label is not None:
            self.help_zoom_label.text = f"Zoom: {int(round(self.help_zoom * 100))}%"

    async def zoom_help_in(self, widget):
        if self.help_pdf_path is None:
            return
        self.help_zoom = min(3.0, self.help_zoom + 0.15)
        await self.refresh_help_view()

    async def zoom_help_out(self, widget):
        if self.help_pdf_path is None:
            return
        self.help_zoom = max(0.6, self.help_zoom - 0.15)
        await self.refresh_help_view()

    async def zoom_help_reset(self, widget):
        if self.help_pdf_path is None:
            return
        self.help_zoom = 1.35
        await self.refresh_help_view()

    def handle_help_window_close(self, widget, **kwargs):
        self.help_window = None
        self.help_webview = None
        self.help_zoom_label = None
        return True

    def cleanup_help_resources(self):
        if self.help_window is not None:
            try:
                if not getattr(self.help_window, "closed", False):
                    self.help_window.close()
            except Exception:
                pass

        self.help_window = None
        self.help_webview = None
        self.help_zoom_label = None

        if self.help_httpd is not None:
            try:
                self.help_httpd.shutdown()
                self.help_httpd.server_close()
            except Exception:
                pass
            self.help_httpd = None

        self.help_server_thread = None
        self.help_server_port = None

        if self.help_temp_dir is not None:
            try:
                self.help_temp_dir.cleanup()
            except Exception:
                pass
            self.help_temp_dir = None

    def ensure_help_server(self):
        if self.help_httpd is not None and self.help_temp_dir is not None and self.help_server_port is not None:
            return

        self.help_temp_dir = tempfile.TemporaryDirectory(prefix="cpx_help_")
        directory = self.help_temp_dir.name

        class HelpRequestHandler(http.server.SimpleHTTPRequestHandler):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, directory=directory, **kwargs)

            def log_message(self, format, *args):
                return

        httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), HelpRequestHandler)
        server_thread = threading.Thread(target=httpd.serve_forever, daemon=True)
        server_thread.start()

        self.help_httpd = httpd
        self.help_server_thread = server_thread
        self.help_server_port = httpd.server_address[1]

    def write_help_html_file(self) -> str:
        self.ensure_help_server()
        if self.help_temp_dir is None or self.help_server_port is None or self.help_pdf_path is None:
            raise RuntimeError("Help viewer server is not initialized.")

        html_path = Path(self.help_temp_dir.name) / "help_viewer.html"
        html_path.write_text(self.build_help_html(self.help_pdf_path), encoding="utf-8")
        cache_buster = int(time.time() * 1000)
        return f"http://127.0.0.1:{self.help_server_port}/help_viewer.html?ts={cache_buster}"

    def make_optional_cards_box(self) -> toga.Box:
        outer = toga.Box(style=Pack(direction=COLUMN, margin_bottom=20))
        outer.add(toga.Label("Optional cards", style=Pack(font_weight="bold", margin_bottom=6)))

        # ATOMIC_VELOCITIES
        outer.add(toga.Label("ATOMIC_VELOCITIES (one row per atom label if you want to include it)", style=Pack(margin_bottom=4)))
        self.atomic_velocities_box = toga.MultilineTextInput(style=Pack(height=80, margin_bottom=10))
        self.atomic_velocities_box.placeholder = "Example:\nH 0.0000 0.0000 0.0000\nO 0.0000 0.0000 0.0000"
        outer.add(self.atomic_velocities_box)

        # CELL_PARAMETERS
        cell_title = toga.Box(style=Pack(direction=ROW, margin_bottom=4))
        cell_title.add(toga.Label("CELL_PARAMETERS", style=Pack(width=150, font_weight="bold")))
        self.cell_parameters_unit = toga.Selection(items=["alat", "bohr", "angstrom"], value="angstrom", style=Pack(width=160))
        cell_title.add(self.cell_parameters_unit)
        outer.add(cell_title)
        self.cell_parameters_text = toga.MultilineTextInput(style=Pack(height=70, margin_bottom=10))
        self.cell_parameters_text.placeholder = "Three lines: v1x v1y v1z\nv2x v2y v2z\nv3x v3y v3z"
        outer.add(self.cell_parameters_text)

        ref_title = toga.Box(style=Pack(direction=ROW, margin_bottom=4))
        ref_title.add(toga.Label("REF_CELL_PARAMETERS", style=Pack(width=150, font_weight="bold")))
        self.ref_cell_parameters_unit = toga.Selection(items=["alat", "bohr", "angstrom"], value="angstrom", style=Pack(width=160))
        ref_title.add(self.ref_cell_parameters_unit)
        outer.add(ref_title)
        self.ref_cell_parameters_text = toga.MultilineTextInput(style=Pack(height=70, margin_bottom=10))
        self.ref_cell_parameters_text.placeholder = "Optional, same three-line format"
        outer.add(self.ref_cell_parameters_text)

        outer.add(toga.Label("CONSTRAINTS", style=Pack(font_weight="bold", margin_bottom=4)))
        self.constraints_text = toga.MultilineTextInput(style=Pack(height=70, margin_bottom=10))
        self.constraints_text.placeholder = "Paste the CONSTRAINTS content here"
        outer.add(self.constraints_text)

        outer.add(toga.Label("OCCUPATIONS", style=Pack(font_weight="bold", margin_bottom=4)))
        self.occupations_text = toga.MultilineTextInput(style=Pack(height=70, margin_bottom=10))
        self.occupations_text.placeholder = "Paste the OCCUPATIONS content here"
        outer.add(self.occupations_text)

        outer.add(toga.Label("ATOMIC_FORCES", style=Pack(font_weight="bold", margin_bottom=4)))
        self.atomic_forces_text = toga.MultilineTextInput(style=Pack(height=80, margin_bottom=10))
        self.atomic_forces_text.placeholder = "Example:\nH 0.0 0.0 0.0\nO 0.0 0.0 0.0"
        outer.add(self.atomic_forces_text)

        outer.add(toga.Label("PLOT_WANNIER", style=Pack(font_weight="bold", margin_bottom=4)))
        self.plot_wannier_text = toga.MultilineTextInput(style=Pack(height=70, margin_bottom=10))
        self.plot_wannier_text.placeholder = "One Wannier index per line"
        outer.add(self.plot_wannier_text)

        outer.add(toga.Label("AUTOPILOT", style=Pack(font_weight="bold", margin_bottom=4)))
        self.autopilot_text = toga.MultilineTextInput(style=Pack(height=120, margin_bottom=10))
        self.autopilot_text.placeholder = "Example:\non_step = 31 : dt = 5.0\non_step = 91 : iprint = 100"
        outer.add(self.autopilot_text)

        note = (
            "Notes: ATOMIC_SPECIES is inferred from the XYZ symbols, with default atomic masses and a placeholder pseudopotential name '<Element>.UPF'. "
            "For cp.x, required information such as pseudopotentials, ecutwfc, and usually cell data must be supplied manually."
        )
        outer.add(toga.Label(note, style=Pack(margin_top=8)))
        return outer

    async def open_xyz_file(self, widget):
        self.set_action_buttons_enabled(False)
        try:
            file = await self.main_window.dialog(toga.OpenFileDialog(title="Select XYZ coordinates file", file_types=["xyz"]))
            if not file:
                return
            self.xyz_path = str(file)
            self.atoms = self.builder.parse_xyz(self.xyz_path)
            self.xyz_label.text = f"Loaded XYZ: {self.xyz_path} | Atoms: {len(self.atoms)}"

            # Update SYSTEM defaults dependent on XYZ
            system_widgets = self.section_widgets["SYSTEM"]
            system_widgets["nat"].value = str(len(self.atoms))
            system_widgets["ntyp"].value = str(len(self.builder.species_from_atoms(self.atoms)))

            self.rebuild_species_box(self.builder.default_species_rows(self.atoms))
            self.render_atomic_positions_text()
            self.rebuild_if_pos_box()
        except Exception as exc:
            await self.main_window.dialog(toga.ErrorDialog(title="XYZ read error", message=str(exc)))
        finally:
            self.set_action_buttons_enabled(True)

    async def open_help_pdf(self, widget):
        try:
            pdf_path = self.resolve_help_pdf_path()
            if pdf_path is None:
                raise FileNotFoundError(
                    "Could not find cpx_input_description.pdf. Place it next to cpx_input_builder.py or update the help PDF path."
                )
            self.help_pdf_path = pdf_path

            if self.help_window is not None and getattr(self.help_window, "closed", False):
                self.handle_help_window_close(self.help_window)

            if self.help_window is None or self.help_webview is None:
                help_window = toga.Window(
                    title="cp.x input description",
                    size=(980, 760),
                    on_close=self.handle_help_window_close,
                )
                help_box = toga.Box(style=Pack(direction=COLUMN, flex=1, margin=10))

                toolbar = toga.Box(style=Pack(direction=ROW, margin_bottom=10))
                toolbar.add(toga.Label("Document viewer", style=Pack(margin_top=6, margin_right=12, font_weight="bold")))
                toolbar.add(toga.Button("Zoom -", on_press=self.zoom_help_out, style=Pack(margin_right=8)))
                toolbar.add(toga.Button("Zoom +", on_press=self.zoom_help_in, style=Pack(margin_right=8)))
                toolbar.add(toga.Button("Reset zoom", on_press=self.zoom_help_reset, style=Pack(margin_right=12)))
                self.help_zoom_label = toga.Label("", style=Pack(margin_top=6))
                toolbar.add(self.help_zoom_label)

                self.help_webview = toga.WebView(style=Pack(flex=1))
                help_box.add(toolbar)
                help_box.add(self.help_webview)
                help_window.content = help_box
                self.help_window = help_window

            self.help_window.show()
            await self.refresh_help_view()
        except Exception as exc:
            await self.main_window.dialog(toga.ErrorDialog(title="Help file error", message=str(exc)))

    def collect_section_values(self, section_name: str) -> Dict[str, Any]:
        values: Dict[str, Any] = {}
        for spec in self.builder.SECTION_SPECS[section_name]:
            widget = self.section_widgets[section_name][spec.name]
            value = getattr(widget, "value", None)
            if spec.kind == "bool":
                values[spec.name] = str(value).strip().lower() in {".true.", "true", "1", "yes"}
            else:
                values[spec.name] = value
        return values

    @staticmethod
    def parse_simple_rows(text: str, ncols: int, keys: List[str]) -> List[Dict[str, str]]:
        rows: List[Dict[str, str]] = []
        for raw in text.splitlines():
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("!"):
                continue
            parts = line.split()
            if len(parts) < ncols:
                raise ValueError(f"Invalid row: '{line}'")
            rows.append({key: parts[i] for i, key in enumerate(keys)})
        return rows

    @staticmethod
    def parse_vector_rows(text: str) -> List[Dict[str, str]]:
        rows: List[Dict[str, str]] = []
        for raw in text.splitlines():
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("!"):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Invalid vector row: '{line}'")
            rows.append({"x": parts[0], "y": parts[1], "z": parts[2]})
        return rows

    def collect_payload(self) -> Dict[str, Any]:
        if not self.xyz_path or not self.atoms:
            raise ValueError("Load an XYZ file first.")

        sections = {name: self.collect_section_values(name) for name in self.builder.SECTION_SPECS}

        atomic_species = [
            {k: widget.value for k, widget in row.items()} for row in self.species_rows
        ]
        if_pos_flags = [(a.value, b.value, c.value) for a, b, c in self.if_pos_rows]

        payload = {
            "sections": sections,
            "atoms": self.atoms,
            "atomic_species": atomic_species,
            "atomic_positions_unit": self.atomic_positions_unit.value or "angstrom",
            "if_pos_flags": if_pos_flags,
            "atomic_velocities": self.parse_simple_rows(self.atomic_velocities_box.value or "", 4, ["label", "vx", "vy", "vz"]),
            "cell_parameters_unit": self.cell_parameters_unit.value or "angstrom",
            "cell_parameters": self.parse_vector_rows(self.cell_parameters_text.value or ""),
            "ref_cell_parameters_unit": self.ref_cell_parameters_unit.value or "angstrom",
            "ref_cell_parameters": self.parse_vector_rows(self.ref_cell_parameters_text.value or ""),
            "constraints_text": self.constraints_text.value or "",
            "occupations_text": self.occupations_text.value or "",
            "atomic_forces": self.parse_simple_rows(self.atomic_forces_text.value or "", 4, ["label", "fx", "fy", "fz"]),
            "plot_wannier_text": self.plot_wannier_text.value or "",
            "autopilot_text": self.autopilot_text.value or "",
        }
        return payload

    async def preview_input(self, widget):
        try:
            payload = self.collect_payload()
            text = self.builder.generate_input_text(payload)
            preview_window = toga.Window(title="cp.x input preview", size=(900, 700))
            preview_box = toga.Box(style=Pack(direction=COLUMN, flex=1, margin=10))
            preview_scroll = toga.ScrollContainer(style=Pack(flex=1))
            preview_text = toga.MultilineTextInput(value=text, readonly=True, style=Pack(flex=1))
            preview_scroll.content = preview_text
            preview_box.add(preview_scroll)
            preview_window.content = preview_box
            preview_window.show()
            self.preview_window = preview_window
        except Exception as exc:
            await self.main_window.dialog(toga.ErrorDialog(title="Preview error", message=str(exc)))

    async def save_cp_input(self, widget):
        try:
            payload = self.collect_payload()
            text = self.builder.generate_input_text(payload)
            saved_path = self.builder.save_input(self.xyz_path, text, self.output_name_input.value or None)
            await self.main_window.dialog(toga.InfoDialog(title="Input saved", message=f"cp.x input saved to:\n{saved_path}"))
        except Exception as exc:
            await self.main_window.dialog(toga.ErrorDialog(title="Save error", message=str(exc)))


class CPInputBuilderApp(toga.App):
    """Standalone Toga app wrapper for the cp.x input builder."""

    def startup(self):
        self.ui = CPInputBuilderUI(app=self)
        self.main_window = self.ui.main_window

    def on_exit(self):
        self.ui.cleanup_help_resources()
        return True
