from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Tuple

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW


@dataclass
class ChoiceSpec:
    name: str
    choices: List[str] = field(default_factory=list)
    default: str = ""
    help_text: str = ""


class ORCAInputBuilder:
    """Backend logic for building common ORCA input files."""

    SIMPLE_SPECS: Dict[str, ChoiceSpec] = {
        "method": ChoiceSpec(
            "Method",
            ["B3LYP", "PBE0", "BP86", "M06-2X", "HF", "MP2", "DLPNO-CCSD(T)", "XTB2", "Custom"],
            "B3LYP",
        ),
        "basis": ChoiceSpec(
            "Basis set",
            ["def2-SVP", "def2-TZVP", "def2-TZVPP", "ma-def2-TZVP", "cc-pVDZ", "cc-pVTZ", "x2c-TZVPall", "Custom"],
            "def2-SVP",
        ),
        "job": ChoiceSpec(
            "Run type",
            ["SP", "Engrad", "Opt", "TightOpt", "VeryTightOpt", "Freq", "NumFreq", "Opt Freq", "TightOpt Freq", "Opt NumFreq"],
            "SP",
        ),
        "dispersion": ChoiceSpec(
            "Dispersion",
            ["None", "D3BJ", "D4"],
            "None",
        ),
        "scf_keyword": ChoiceSpec(
            "SCF keyword",
            ["Default", "NormalSCF", "TightSCF", "VeryTightSCF", "SlowConv"],
            "Default",
        ),
        "ri_keyword": ChoiceSpec(
            "RI keyword",
            ["Default", "RI", "RIJCOSX", "NoRI"],
            "Default",
        ),
        "print_keyword": ChoiceSpec(
            "Print keyword",
            ["Default", "MiniPrint", "SmallPrint", "NormalPrint", "LargePrint", "PrintBasis", "PrintMOs"],
            "Default",
        ),
        "solvation_model": ChoiceSpec(
            "Solvation",
            ["None", "CPCM", "SMD"],
            "None",
        ),
        "coordinates_mode": ChoiceSpec(
            "Coordinates mode",
            ["Inline XYZ", "XYZFILE reference"],
            "Inline XYZ",
        ),
        "scf_convergence": ChoiceSpec(
            "SCF convergence block",
            ["Default", "Loose", "Medium", "Strong", "Tight", "VeryTight", "Extreme"],
            "Default",
        ),
    }

    RESEARCH_SUMMARY = (
        "ORCA 6.1 input essentials used in this builder:\n"
        "- General input structure: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html\n"
        "- Coordinates and xyzfile syntax: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html\n"
        "- Basic run types: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basics.html\n"
        "- Parallel %pal usage: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html\n"
        "- SCF convergence keywords: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html\n"
        "- Basis-set usage: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html\n"
        "- Implicit solvation (CPCM / SMD): https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html\n"
        "- Geometry optimization keywords: https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html\n\n"
        "Implemented in this first version:\n"
        "- One main ! keyword line with common method, basis, run type, SCF, RI, dispersion, print, and solvent options.\n"
        "- Optional %pal, %maxcore, %scf, %geom, and %cpcm blocks.\n"
        "- Inline '* xyz charge mult' or external '* xyzfile charge mult filename' coordinates.\n"
        "- Free-form advanced blocks for uncommon ORCA options."
    )

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

    @staticmethod
    def split_keywords(raw: str) -> List[str]:
        tokens: List[str] = []
        for line in raw.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens.extend(stripped.split())
        return tokens

    @staticmethod
    def parse_positive_int(raw: str, label: str) -> Optional[int]:
        value = (raw or "").strip()
        if not value:
            return None
        parsed = int(value)
        if parsed < 1:
            raise ValueError(f"{label} must be a positive integer.")
        return parsed

    @staticmethod
    def parse_nonnegative_int(raw: str, label: str) -> Optional[int]:
        value = (raw or "").strip()
        if not value:
            return None
        parsed = int(value)
        if parsed < 0:
            raise ValueError(f"{label} must be a non-negative integer.")
        return parsed

    def build_simple_keywords(self, payload: Dict[str, Any]) -> List[str]:
        tokens: List[str] = []

        method = payload["method_custom"].strip() if payload["method"] == "Custom" else payload["method"]
        basis = payload["basis_custom"].strip() if payload["basis"] == "Custom" else payload["basis"]
        if not method:
            raise ValueError("An ORCA method keyword is required.")
        if not basis:
            raise ValueError("A basis-set keyword is required.")
        tokens.append(method)
        tokens.append(basis)

        job_value = payload["job"]
        if job_value and job_value != "SP":
            tokens.extend(job_value.split())

        for key in ["dispersion", "scf_keyword", "ri_keyword", "print_keyword"]:
            value = payload[key]
            if value and value != "Default" and value != "None":
                tokens.append(value)

        solvation_model = payload["solvation_model"]
        solvent = payload["solvent"].strip()
        if solvation_model != "None":
            if not solvent:
                raise ValueError("Provide a solvent name when using CPCM or SMD.")
            tokens.append(f"{solvation_model}({solvent})")

        tokens.extend(self.split_keywords(payload.get("extra_simple_keywords", "")))
        return tokens

    @staticmethod
    def build_block(name: str, lines: List[str]) -> str:
        useful = [line.rstrip() for line in lines if line and line.strip()]
        if not useful:
            return ""
        block_lines = [f"%{name}"]
        block_lines.extend(f"  {line}" for line in useful)
        block_lines.append("end")
        return "\n".join(block_lines)

    def build_scf_block(self, payload: Dict[str, Any]) -> str:
        lines: List[str] = []
        convergence = payload["scf_convergence"]
        max_iter = self.parse_positive_int(payload.get("scf_maxiter", ""), "SCF MaxIter")
        if convergence != "Default":
            lines.append(f"Convergence {convergence}")
        if max_iter is not None:
            lines.append(f"MaxIter {max_iter}")
        raw_text = payload.get("scf_block_text", "").strip()
        if raw_text:
            lines.extend(raw_text.splitlines())
        return self.build_block("scf", lines)

    def build_geom_block(self, payload: Dict[str, Any]) -> str:
        lines: List[str] = []
        max_iter = self.parse_positive_int(payload.get("geom_maxiter", ""), "Geometry MaxIter")
        if max_iter is not None:
            lines.append(f"MaxIter {max_iter}")
        raw_text = payload.get("geom_block_text", "").strip()
        if raw_text:
            lines.extend(raw_text.splitlines())
        return self.build_block("geom", lines)

    def build_cpcm_block(self, payload: Dict[str, Any]) -> str:
        raw_text = payload.get("cpcm_block_text", "").strip()
        if not raw_text:
            return ""
        return self.build_block("cpcm", raw_text.splitlines())

    def build_parallel_block(self, payload: Dict[str, Any]) -> str:
        nprocs = self.parse_positive_int(payload.get("nprocs", ""), "Number of processes")
        if nprocs is None or nprocs <= 1:
            return ""
        return self.build_block("pal", [f"nprocs {nprocs}"])

    def build_maxcore_block(self, payload: Dict[str, Any]) -> str:
        maxcore = self.parse_positive_int(payload.get("maxcore", ""), "MaxCore")
        if maxcore is None:
            return ""
        return f"%maxcore {maxcore}"

    def build_coordinates_block(
        self,
        atoms: List[Tuple[str, float, float, float]],
        xyz_path: str,
        charge: int,
        multiplicity: int,
        mode: str,
    ) -> str:
        if mode == "XYZFILE reference":
            filename = Path(xyz_path).name
            return f"* xyzfile {charge} {multiplicity} {filename}"

        lines = [f"* xyz {charge} {multiplicity}"]
        for symbol, x, y, z in atoms:
            lines.append(f"{symbol:<3} {x: .10f} {y: .10f} {z: .10f}")
        lines.append("*")
        return "\n".join(lines)

    def validate_payload(self, payload: Dict[str, Any]) -> None:
        charge = int((payload.get("charge", "0") or "0").strip())
        multiplicity = int((payload.get("multiplicity", "1") or "1").strip())
        if multiplicity < 1:
            raise ValueError("Multiplicity must be at least 1.")
        payload["charge"] = charge
        payload["multiplicity"] = multiplicity

    def generate_input_text(self, payload: Dict[str, Any]) -> str:
        self.validate_payload(payload)
        atoms = payload.get("atoms", [])
        xyz_path = payload.get("xyz_path", "")
        if not atoms or not xyz_path:
            raise ValueError("Load an XYZ file first.")

        title = payload.get("title", "").strip()
        notes = payload.get("notes", "").strip()
        pieces: List[str] = []

        if title:
            pieces.append(f"# {title}")
        if notes:
            for line in notes.splitlines():
                pieces.append(f"# {line}".rstrip())

        simple_keywords = self.build_simple_keywords(payload)
        pieces.append("! " + " ".join(simple_keywords))

        optional_blocks = [
            self.build_parallel_block(payload),
            self.build_maxcore_block(payload),
            self.build_scf_block(payload),
            self.build_geom_block(payload),
            self.build_cpcm_block(payload),
        ]
        pieces.extend(block for block in optional_blocks if block)

        advanced_blocks = payload.get("advanced_blocks_text", "").strip()
        if advanced_blocks:
            pieces.append(advanced_blocks)

        pieces.append(
            self.build_coordinates_block(
                atoms=atoms,
                xyz_path=xyz_path,
                charge=payload["charge"],
                multiplicity=payload["multiplicity"],
                mode=payload["coordinates_mode"],
            )
        )

        return "\n\n".join(piece for piece in pieces if piece.strip()) + "\n"

    def save_input(self, xyz_path: str, input_text: str, filename: Optional[str] = None) -> Path:
        xyz = Path(xyz_path)
        output_name = filename.strip() if filename else f"{xyz.stem}_orca.inp"
        out_path = xyz.with_name(output_name)
        out_path.write_text(input_text, encoding="utf-8", newline="\n")
        return out_path


class ORCAInputBuilderUI:
    """Toga GUI for building ORCA input files."""

    def __init__(self, *args, app: Optional[toga.App] = None):
        del args
        resolved_app = app if isinstance(app, toga.App) else None
        self.app = resolved_app
        self.formal_name = getattr(resolved_app, "formal_name", "ORCA Input Builder")
        self.builder = ORCAInputBuilder()

        self.xyz_path: Optional[str] = None
        self.atoms: List[Tuple[str, float, float, float]] = []
        self.preview_window: Optional[toga.Window] = None
        self.help_window: Optional[toga.Window] = None

        if self.app is None:
            self.main_window = toga.Window(title=self.formal_name, size=(650, 700))
        else:
            self.main_window = toga.MainWindow(title=self.formal_name)

        self.action_buttons: List[toga.Button] = []
        self._build_ui()
        self.main_window.show()

    def _build_ui(self):
        content = toga.Box(style=Pack(direction=COLUMN, margin=10))

        top_actions = toga.Box(style=Pack(direction=ROW, margin_bottom=10))
        load_btn = toga.Button("Open XYZ", on_press=self.open_xyz_file, style=Pack(margin_right=8))
        preview_btn = toga.Button("Preview ORCA input", on_press=self.preview_input, style=Pack(margin_right=8))
        save_btn = toga.Button("Save ORCA input", on_press=self.save_input, style=Pack(margin_right=8))
        clear_btn = toga.Button("Clear form", on_press=self.clear_form, style=Pack(margin_right=8))
        help_btn = toga.Button("Research notes", on_press=self.show_research_notes, style=Pack())
        self.action_buttons = [load_btn, preview_btn, save_btn, clear_btn, help_btn]
        for button in self.action_buttons:
            top_actions.add(button)

        self.xyz_label = toga.Label("No XYZ file loaded.", style=Pack(margin_bottom=8))

        self.output_name_input = toga.TextInput(
            value="",
            placeholder="Leave blank to use <xyzname>_orca.inp",
            style=Pack(width=280, margin_left=8),
        )

        output_row = toga.Box(style=Pack(direction=ROW, margin_bottom=10))
        output_row.add(toga.Label("Output filename:", style=Pack(width=130, margin_top=6)))
        output_row.add(self.output_name_input)

        self.scroll = toga.ScrollContainer(style=Pack(flex=1))
        scroll_content = toga.Box(style=Pack(direction=COLUMN, margin=5))

        scroll_content.add(self.make_text_metadata_box())
        scroll_content.add(self.make_coordinates_box())
        scroll_content.add(self.make_simple_keywords_box())
        scroll_content.add(self.make_blocks_box())

        self.scroll.content = scroll_content

        content.add(top_actions)
        content.add(self.xyz_label)
        content.add(output_row)
        content.add(self.scroll)

        self.main_window.content = content
        self.update_custom_field_state()

    def make_text_metadata_box(self) -> toga.Box:
        outer = toga.Box(style=Pack(direction=COLUMN, margin_bottom=14))
        outer.add(toga.Label("Job metadata", style=Pack(font_weight="bold", margin_bottom=6)))

        title_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        title_row.add(toga.Label("Title/comment", style=Pack(width=130, margin_top=6)))
        self.title_input = toga.TextInput(placeholder="Optional comment line written as '# ...'", style=Pack(flex=1))
        title_row.add(self.title_input)
        outer.add(title_row)

        notes_row = toga.Box(style=Pack(direction=COLUMN, margin_bottom=4))
        notes_row.add(toga.Label("Extra comment lines", style=Pack(margin_bottom=4)))
        self.notes_input = toga.MultilineTextInput(style=Pack(height=65))
        self.notes_input.placeholder = "Optional notes written as ORCA comment lines"
        notes_row.add(self.notes_input)
        outer.add(notes_row)
        return outer

    def make_coordinates_box(self) -> toga.Box:
        outer = toga.Box(style=Pack(direction=COLUMN, margin_bottom=14))
        outer.add(toga.Label("Coordinates", style=Pack(font_weight="bold", margin_bottom=6)))

        mode_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        mode_row.add(toga.Label("Coordinate mode", style=Pack(width=130, margin_top=6)))
        self.coordinates_mode = toga.Selection(
            items=self.builder.SIMPLE_SPECS["coordinates_mode"].choices,
            value=self.builder.SIMPLE_SPECS["coordinates_mode"].default,
            style=Pack(width=180),
        )
        mode_row.add(self.coordinates_mode)
        outer.add(mode_row)

        charge_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        charge_row.add(toga.Label("Charge", style=Pack(width=130, margin_top=6)))
        self.charge_input = toga.TextInput(value="0", style=Pack(width=100, margin_right=20))
        charge_row.add(self.charge_input)
        charge_row.add(toga.Label("Multiplicity", style=Pack(width=90, margin_top=6)))
        self.multiplicity_input = toga.TextInput(value="1", style=Pack(width=100))
        charge_row.add(self.multiplicity_input)
        outer.add(charge_row)

        outer.add(toga.Label("Loaded coordinates preview", style=Pack(margin_bottom=4)))
        self.coordinates_preview = toga.MultilineTextInput(readonly=True, style=Pack(height=140))
        self.coordinates_preview.placeholder = "Load an XYZ file to preview the ORCA coordinate section."
        outer.add(self.coordinates_preview)
        return outer

    def make_simple_keywords_box(self) -> toga.Box:
        outer = toga.Box(style=Pack(direction=COLUMN, margin_bottom=14))
        outer.add(toga.Label("Simple ! keywords", style=Pack(font_weight="bold", margin_bottom=6)))

        self.method_selection = self.make_choice_row(outer, "Method", "method", on_change=self.on_method_change)
        self.method_custom = self.make_text_row(outer, "Custom method", "Only used when Method = Custom")
        self.basis_selection = self.make_choice_row(outer, "Basis set", "basis", on_change=self.on_basis_change)
        self.basis_custom = self.make_text_row(outer, "Custom basis", "Only used when Basis set = Custom")
        self.job_selection = self.make_choice_row(outer, "Run type", "job")
        self.dispersion_selection = self.make_choice_row(outer, "Dispersion", "dispersion")
        self.scf_keyword_selection = self.make_choice_row(outer, "SCF keyword", "scf_keyword")
        self.ri_keyword_selection = self.make_choice_row(outer, "RI keyword", "ri_keyword")
        self.print_keyword_selection = self.make_choice_row(outer, "Print keyword", "print_keyword")
        self.solvation_selection = self.make_choice_row(outer, "Solvation", "solvation_model")
        self.solvent_input = self.make_text_row(outer, "Solvent", "Examples: Water, Acetonitrile, Methanol")

        extra_row = toga.Box(style=Pack(direction=COLUMN, margin_bottom=6))
        extra_row.add(toga.Label("Extra simple keywords", style=Pack(margin_bottom=4)))
        self.extra_simple_keywords = toga.MultilineTextInput(style=Pack(height=70))
        self.extra_simple_keywords.placeholder = "Examples:\nUKS\nGrid5 FinalGrid6\nUseSym"
        extra_row.add(self.extra_simple_keywords)
        outer.add(extra_row)
        return outer

    def make_blocks_box(self) -> toga.Box:
        outer = toga.Box(style=Pack(direction=COLUMN, margin_bottom=14))
        outer.add(toga.Label("Optional ORCA blocks", style=Pack(font_weight="bold", margin_bottom=6)))

        parallel_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        parallel_row.add(toga.Label("nprocs", style=Pack(width=130, margin_top=6)))
        self.nprocs_input = toga.TextInput(value="", placeholder="e.g. 4", style=Pack(width=100, margin_right=20))
        parallel_row.add(self.nprocs_input)
        parallel_row.add(toga.Label("MaxCore (MB/core)", style=Pack(width=120, margin_top=6)))
        self.maxcore_input = toga.TextInput(value="", placeholder="e.g. 2000", style=Pack(width=120))
        parallel_row.add(self.maxcore_input)
        outer.add(parallel_row)

        scf_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        scf_row.add(toga.Label("SCF Convergence", style=Pack(width=130, margin_top=6)))
        self.scf_convergence_selection = toga.Selection(
            items=self.builder.SIMPLE_SPECS["scf_convergence"].choices,
            value=self.builder.SIMPLE_SPECS["scf_convergence"].default,
            style=Pack(width=150, margin_right=20),
        )
        scf_row.add(self.scf_convergence_selection)
        scf_row.add(toga.Label("SCF MaxIter", style=Pack(width=90, margin_top=6)))
        self.scf_maxiter_input = toga.TextInput(value="", placeholder="e.g. 200", style=Pack(width=120))
        scf_row.add(self.scf_maxiter_input)
        outer.add(scf_row)

        geom_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        geom_row.add(toga.Label("Geom MaxIter", style=Pack(width=130, margin_top=6)))
        self.geom_maxiter_input = toga.TextInput(value="", placeholder="e.g. 150", style=Pack(width=120))
        geom_row.add(self.geom_maxiter_input)
        outer.add(geom_row)

        outer.add(toga.Label("%scf extra lines", style=Pack(margin_bottom=4)))
        self.scf_block_text = toga.MultilineTextInput(style=Pack(height=70, margin_bottom=8))
        self.scf_block_text.placeholder = "Examples:\nSOSCFStart 0.00033\nCNVDamp true"
        outer.add(self.scf_block_text)

        outer.add(toga.Label("%geom extra lines", style=Pack(margin_bottom=4)))
        self.geom_block_text = toga.MultilineTextInput(style=Pack(height=70, margin_bottom=8))
        self.geom_block_text.placeholder = "Examples:\nCalc_Hess true\nRecalc_Hess 5"
        outer.add(self.geom_block_text)

        outer.add(toga.Label("%cpcm extra lines", style=Pack(margin_bottom=4)))
        self.cpcm_block_text = toga.MultilineTextInput(style=Pack(height=70, margin_bottom=8))
        self.cpcm_block_text.placeholder = "Examples:\ncds_cpcm 2\nDRACO true"
        outer.add(self.cpcm_block_text)

        outer.add(toga.Label("Advanced raw blocks", style=Pack(margin_bottom=4)))
        self.advanced_blocks_text = toga.MultilineTextInput(style=Pack(height=120))
        self.advanced_blocks_text.placeholder = "Paste complete ORCA blocks here, for example:\n%basis\n  AuxJ \"def2/J\"\nend"
        outer.add(self.advanced_blocks_text)
        return outer

    def make_choice_row(
        self,
        parent: toga.Box,
        label_text: str,
        spec_key: str,
        on_change=None,
    ) -> toga.Selection:
        row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        row.add(toga.Label(label_text, style=Pack(width=130, margin_top=6)))
        spec = self.builder.SIMPLE_SPECS[spec_key]
        widget = toga.Selection(items=spec.choices, value=spec.default, on_change=on_change, style=Pack(width=180))
        row.add(widget)
        parent.add(row)
        return widget

    def make_text_row(self, parent: toga.Box, label_text: str, placeholder: str) -> toga.TextInput:
        row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        row.add(toga.Label(label_text, style=Pack(width=130, margin_top=6)))
        widget = toga.TextInput(placeholder=placeholder, style=Pack(flex=1))
        row.add(widget)
        parent.add(row)
        return widget

    def set_action_buttons_enabled(self, enabled: bool):
        for button in self.action_buttons:
            button.enabled = enabled

    def update_custom_field_state(self):
        method_is_custom = (self.method_selection.value or "") == "Custom"
        basis_is_custom = (self.basis_selection.value or "") == "Custom"
        self.method_custom.enabled = method_is_custom
        self.basis_custom.enabled = basis_is_custom
        if not method_is_custom:
            self.method_custom.value = ""
        if not basis_is_custom:
            self.basis_custom.value = ""

    def on_method_change(self, widget):
        self.update_custom_field_state()

    def on_basis_change(self, widget):
        self.update_custom_field_state()

    async def clear_form(self, widget):
        del widget
        self.xyz_path = None
        self.atoms = []
        self.xyz_label.text = "No XYZ file loaded."
        self.output_name_input.value = ""
        self.title_input.value = ""
        self.notes_input.value = ""
        self.coordinates_mode.value = self.builder.SIMPLE_SPECS["coordinates_mode"].default
        self.charge_input.value = "0"
        self.multiplicity_input.value = "1"
        self.method_selection.value = self.builder.SIMPLE_SPECS["method"].default
        self.basis_selection.value = self.builder.SIMPLE_SPECS["basis"].default
        self.job_selection.value = self.builder.SIMPLE_SPECS["job"].default
        self.dispersion_selection.value = self.builder.SIMPLE_SPECS["dispersion"].default
        self.scf_keyword_selection.value = self.builder.SIMPLE_SPECS["scf_keyword"].default
        self.ri_keyword_selection.value = self.builder.SIMPLE_SPECS["ri_keyword"].default
        self.print_keyword_selection.value = self.builder.SIMPLE_SPECS["print_keyword"].default
        self.solvation_selection.value = self.builder.SIMPLE_SPECS["solvation_model"].default
        self.solvent_input.value = ""
        self.extra_simple_keywords.value = ""
        self.nprocs_input.value = ""
        self.maxcore_input.value = ""
        self.scf_convergence_selection.value = self.builder.SIMPLE_SPECS["scf_convergence"].default
        self.scf_maxiter_input.value = ""
        self.geom_maxiter_input.value = ""
        self.scf_block_text.value = ""
        self.geom_block_text.value = ""
        self.cpcm_block_text.value = ""
        self.advanced_blocks_text.value = ""
        self.update_custom_field_state()
        self.render_coordinates_preview()

    def render_coordinates_preview(self):
        if not self.atoms or not self.xyz_path:
            self.coordinates_preview.value = ""
            return
        try:
            charge = int((self.charge_input.value or "0").strip())
            multiplicity = int((self.multiplicity_input.value or "1").strip())
        except ValueError:
            self.coordinates_preview.value = "Charge and multiplicity must be integers."
            return

        self.coordinates_preview.value = self.builder.build_coordinates_block(
            atoms=self.atoms,
            xyz_path=self.xyz_path,
            charge=charge,
            multiplicity=multiplicity,
            mode=self.coordinates_mode.value or "Inline XYZ",
        )

    async def open_xyz_file(self, widget):
        self.set_action_buttons_enabled(False)
        try:
            file = await self.main_window.dialog(
                toga.OpenFileDialog(title="Select XYZ coordinates file", file_types=["xyz"])
            )
            if not file:
                return
            self.xyz_path = str(file)
            self.atoms = self.builder.parse_xyz(self.xyz_path)
            self.xyz_label.text = f"Loaded XYZ: {self.xyz_path} | Atoms: {len(self.atoms)}"
            self.render_coordinates_preview()
        except Exception as exc:
            await self.main_window.dialog(toga.ErrorDialog(title="XYZ read error", message=str(exc)))
        finally:
            self.set_action_buttons_enabled(True)

    def collect_payload(self) -> Dict[str, Any]:
        if not self.xyz_path or not self.atoms:
            raise ValueError("Load an XYZ file first.")

        payload = {
            "xyz_path": self.xyz_path,
            "atoms": self.atoms,
            "title": self.title_input.value or "",
            "notes": self.notes_input.value or "",
            "charge": self.charge_input.value or "0",
            "multiplicity": self.multiplicity_input.value or "1",
            "coordinates_mode": self.coordinates_mode.value or "Inline XYZ",
            "method": self.method_selection.value or "B3LYP",
            "method_custom": self.method_custom.value or "",
            "basis": self.basis_selection.value or "def2-SVP",
            "basis_custom": self.basis_custom.value or "",
            "job": self.job_selection.value or "SP",
            "dispersion": self.dispersion_selection.value or "None",
            "scf_keyword": self.scf_keyword_selection.value or "Default",
            "ri_keyword": self.ri_keyword_selection.value or "Default",
            "print_keyword": self.print_keyword_selection.value or "Default",
            "solvation_model": self.solvation_selection.value or "None",
            "solvent": self.solvent_input.value or "",
            "extra_simple_keywords": self.extra_simple_keywords.value or "",
            "nprocs": self.nprocs_input.value or "",
            "maxcore": self.maxcore_input.value or "",
            "scf_convergence": self.scf_convergence_selection.value or "Default",
            "scf_maxiter": self.scf_maxiter_input.value or "",
            "geom_maxiter": self.geom_maxiter_input.value or "",
            "scf_block_text": self.scf_block_text.value or "",
            "geom_block_text": self.geom_block_text.value or "",
            "cpcm_block_text": self.cpcm_block_text.value or "",
            "advanced_blocks_text": self.advanced_blocks_text.value or "",
        }
        return payload

    async def preview_input(self, widget):
        try:
            self.render_coordinates_preview()
            payload = self.collect_payload()
            text = self.builder.generate_input_text(payload)
            preview_window = toga.Window(title="ORCA input preview", size=(900, 720))
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

    async def save_input(self, widget):
        try:
            self.render_coordinates_preview()
            payload = self.collect_payload()
            text = self.builder.generate_input_text(payload)
            saved_path = self.builder.save_input(self.xyz_path, text, self.output_name_input.value or None)
            await self.main_window.dialog(
                toga.InfoDialog(title="Input saved", message=f"ORCA input saved to:\n{saved_path}")
            )
        except Exception as exc:
            await self.main_window.dialog(toga.ErrorDialog(title="Save error", message=str(exc)))

    async def show_research_notes(self, widget):
        if self.help_window is not None and not getattr(self.help_window, "closed", False):
            self.help_window.show()
            return

        self.help_window = toga.Window(title="ORCA input research notes", size=(860, 620))
        help_box = toga.Box(style=Pack(direction=COLUMN, margin=10, flex=1))
        help_text = toga.MultilineTextInput(
            value=self.builder.RESEARCH_SUMMARY,
            readonly=True,
            style=Pack(flex=1),
        )
        help_box.add(help_text)
        self.help_window.content = help_box
        self.help_window.show()


class ORCAInputBuilderApp(toga.App):
    """Standalone wrapper app for the ORCA input builder."""

    def startup(self):
        self.ui = ORCAInputBuilderUI(app=self)
        self.main_window = self.ui.main_window
