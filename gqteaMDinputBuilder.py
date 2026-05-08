from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW


class GqteaMDInputBuilder:
    """Backend logic for building gqteaMD TOML input files."""

    FORCE_PROVIDER_CHOICES = ["harmonic", "uff", "xtb", "gaussian", "classical"]
    BOOL_CHOICES = ["false", "true"]
    UFF_ELECTROSTATICS_CHOICES = ["auto", "true", "false"]
    UFF_EXCLUSION_CHOICES = ["exclude_12_13", "exclude_12", "none"]
    UFF_LJ_CUTOFF_CHOICES = ["plain", "shift"]
    XTB_METHOD_CHOICES = ["GFN2-xTB", "GFN1-xTB", "GFN0-xTB", "GFN-FF", "Custom"]
    GAUSSIAN_COMMAND_CHOICES = ["g16", "g16.exe", "g09", "g09.exe", "Custom"]
    ROUTE_CHOICES = [
        "# B3LYP/6-31G(d) Force NoSymm SCF=Tight",
        "# HF/3-21G Force NoSymm SCF=Tight",
        "# AM1 Force NoSymm SCF=(MaxCycle=1000)",
        "Custom",
    ]

    MANUAL_SUMMARY = (
        "gqteaMD TOML sections used by this builder:\n"
        "- [input]: starting XYZ geometry.\n"
        "- [cell]: orthorhombic a, b, c box lengths in angstrom.\n"
        "- [dynamics]: velocity Verlet timestep_fs and steps.\n"
        "- [force_provider]: harmonic, classical, uff, xtb, or gaussian.\n"
        "- [output]: trajectory and log settings. GEOMETRY is written automatically.\n"
        "- [restart]: optional restart writing and resume behavior.\n\n"
        "Gaussian notes:\n"
        "- gqteaMD supports command, route, charge, multiplicity, nproc, memory, chk, and workdir.\n"
        "- gqteaMD automatically adds Force and NoSymm to the route if missing.\n"
        "- gqteaMD currently creates per-step checkpoint files automatically.\n"
        "- The default trajectory is TRAJEC.xyz and the default log is <xyzname>_gqteaMD.log.\n"
        "- gqteaMD rewrites GEOMETRY at every calculation step with symbol, positions, velocities, and forces.\n"
        "- resume_from_RESTART resumes directly from RESTART.\n"
        "- resume_from_GEOMETRY reads positions, velocities, and forces from GEOMETRY and energies from RESTART.\n"
        "- The chk and memory fields are optional Gaussian provider settings; blank fields are omitted.\n\n"
        "UFF notes:\n"
        "- gqteaMD detects bonds from covalent radii, assigns simple UFF atom types, and builds "
        "a topology with bonds, angles, torsions, inversions, exclusions, 1-4 pairs, and optional charges.\n"
        "- Current UFF forces use harmonic bonds, harmonic angles, and Lennard-Jones interactions.\n"
        "- UFF now uses bond-order/electronegativity corrected bond lengths, analytic angle forces, "
        "and finite-difference torsion/inversion forces.\n"
        "- UFF charges activate fixed-charge Coulomb electrostatics unless electrostatics is set false.\n"
        "- UFF nonbonded settings include 1-2/1-3 exclusions, 1-4 LJ/Coulomb scaling, and plain or shifted LJ cutoffs."
        "\n- UFF topology can be controlled with explicit bonds, angles, torsions, and inversions; omitted lists are generated."
        "\n- UFF can use a Verlet neighbor list for cutoff nonbonded interactions; the skin controls rebuild frequency."
        "\n- Validation examples include uff_water_charged.toml and uff_ethene_explicit_topology.toml.\n\n"
        "xTB notes:\n"
        "- gqteaMD can use xTB results for single-point energies and Cartesian forces through the xtb-python ASE calculator.\n"
        "- On Windows, set command to the full xtb.exe path, for example C:/xTB/xtb-6.7.1/bin/xtb.exe.\n"
        "- xTB energies are used in eV and forces in eV/angstrom, matching gqteaMD internal units.\n"
        "- Typical settings include method, charge, multiplicity, accuracy, electronic_temperature, max_iterations, and solvent.\n"
        "- omp_num_threads controls OMP_NUM_THREADS for xTB calculations when supported by the gqteaMD runtime.\n"
        "- Install the optional gqteaMD xTB dependencies, or install ASE and xtb-python in the active environment.\n"
        "- use_unwrapped_positions should usually remain true to avoid passing broken molecules across periodic boundaries."
    )

    @staticmethod
    def toml_string(value: str) -> str:
        return json.dumps((value or "").strip())

    @staticmethod
    def parse_float(raw: str, label: str, *, required: bool = True) -> Optional[float]:
        value = (raw or "").strip()
        if not value:
            if required:
                raise ValueError(f"{label} is required.")
            return None
        return float(value)

    @staticmethod
    def parse_positive_int(raw: str, label: str, *, required: bool = True) -> Optional[int]:
        value = (raw or "").strip()
        if not value:
            if required:
                raise ValueError(f"{label} is required.")
            return None
        parsed = int(value)
        if parsed < 1:
            raise ValueError(f"{label} must be a positive integer.")
        return parsed

    @staticmethod
    def parse_int(raw: str, label: str, *, required: bool = True) -> Optional[int]:
        value = (raw or "").strip()
        if not value:
            if required:
                raise ValueError(f"{label} is required.")
            return None
        return int(value)

    @staticmethod
    def parse_bool(raw: str) -> bool:
        return (raw or "false").strip().lower() == "true"

    @staticmethod
    def default_log_name(xyz_path: str) -> str:
        stem = Path(xyz_path).stem or "gqteaMD"
        return f"{stem}_gqteaMD.log"

    @staticmethod
    def parse_atom_types(raw: str) -> list[str]:
        return [item.strip() for item in (raw or "").replace("\n", ",").split(",") if item.strip()]

    @staticmethod
    def parse_float_list(raw: str, label: str) -> list[float]:
        values = [item.strip() for item in (raw or "").replace("\n", ",").split(",") if item.strip()]
        try:
            return [float(item) for item in values]
        except ValueError as exc:
            raise ValueError(f"{label} must contain comma-separated numbers.") from exc

    @staticmethod
    def inspect_xyz(filepath: str) -> tuple[int, str]:
        path = Path(filepath)
        lines = path.read_text(encoding="utf-8").splitlines()
        if len(lines) < 2:
            raise ValueError("The XYZ file must contain at least two lines.")
        try:
            atom_count = int(lines[0].strip())
        except ValueError as exc:
            raise ValueError("The first line of the XYZ file must contain the atom count.") from exc
        if len(lines[2:]) < atom_count:
            raise ValueError("The XYZ file has fewer atom lines than declared.")
        return atom_count, lines[1].strip()

    def generate_toml(self, payload: Dict[str, Any]) -> str:
        xyz_path = (payload["xyz"] or "").strip()
        if not xyz_path:
            raise ValueError("Choose or type an input XYZ file.")

        force_type = (payload["force_type"] or "harmonic").lower()
        if force_type not in self.FORCE_PROVIDER_CHOICES:
            raise ValueError(f"Unsupported force provider: {force_type}")

        lines: list[str] = [
            "# gqteaMD input generated by gqteaWinToga",
            "",
            "[input]",
            f"xyz = {self.toml_string(xyz_path)}",
            "",
            "[cell]",
            f"a = {self.parse_float(payload['cell_a'], 'Cell a')}",
            f"b = {self.parse_float(payload['cell_b'], 'Cell b')}",
            f"c = {self.parse_float(payload['cell_c'], 'Cell c')}",
            "",
            "[dynamics]",
            f"timestep_fs = {self.parse_float(payload['timestep_fs'], 'Timestep')}",
            f"steps = {self.parse_positive_int(payload['steps'], 'Steps')}",
            "",
            "[force_provider]",
            f"type = {self.toml_string(force_type)}",
        ]

        force_tail: list[str] = []
        if force_type == "harmonic":
            lines.append(
                f"k_ev_per_angstrom2 = {self.parse_float(payload['harmonic_k'], 'Harmonic force constant')}"
            )
        elif force_type == "uff":
            lines.append(
                f"bond_detection_scale = {self.parse_float(payload['uff_bond_detection_scale'], 'UFF bond detection scale')}"
            )
            cutoff = self.parse_float(payload.get("uff_cutoff_angstrom", ""), "UFF cutoff", required=False)
            if cutoff is not None:
                lines.append(f"cutoff_angstrom = {cutoff}")
            atom_types = self.parse_atom_types(payload.get("uff_atom_types", ""))
            if atom_types:
                lines.append("atom_types = [" + ", ".join(self.toml_string(item) for item in atom_types) + "]")
            charges = self.parse_float_list(payload.get("uff_charges", ""), "UFF charges")
            if charges:
                lines.append("charges = [" + ", ".join(str(charge) for charge in charges) + "]")
            electrostatics = (payload.get("uff_electrostatics", "auto") or "auto").strip().lower()
            if electrostatics != "auto":
                lines.append(f"electrostatics = {str(self.parse_bool(electrostatics)).lower()}")
            lines.append(f"nonbonded_exclusions = {self.toml_string(payload.get('uff_nonbonded_exclusions', 'exclude_12_13'))}")
            lines.append(f"lj_14_scale = {self.parse_float(payload.get('uff_lj_14_scale', '1.0'), 'UFF LJ 1-4 scale')}")
            lines.append(
                f"electrostatic_14_scale = {self.parse_float(payload.get('uff_electrostatic_14_scale', '1.0'), 'UFF electrostatic 1-4 scale')}"
            )
            lines.append(f"lj_cutoff_mode = {self.toml_string(payload.get('uff_lj_cutoff_mode', 'plain'))}")
            lines.append(f"use_neighbor_list = {str(self.parse_bool(payload.get('uff_use_neighbor_list', 'true'))).lower()}")
            lines.append(
                f"neighbor_skin_angstrom = {self.parse_float(payload.get('uff_neighbor_skin_angstrom', '2.0'), 'UFF neighbor skin')}"
            )
            bond_orders_text = (payload.get("uff_bond_orders_text", "") or "").strip()
            if bond_orders_text:
                lines.append("")
                lines.extend(bond_orders_text.splitlines())
            topology_text = (payload.get("uff_topology_text", "") or "").strip()
            if topology_text:
                lines.append("")
                lines.extend(topology_text.splitlines())
        elif force_type == "xtb":
            method = payload["xtb_method_custom"].strip() if payload["xtb_method"] == "Custom" else payload["xtb_method"]
            if not method:
                raise ValueError("xTB method is required.")
            command = (payload.get("xtb_command", "") or "").strip()
            if command:
                lines.append(f"command = {self.toml_string(command)}")
            lines.extend(
                [
                    f"method = {self.toml_string(method)}",
                    f"charge = {self.parse_float(payload['xtb_charge'], 'xTB charge')}",
                    f"multiplicity = {self.parse_positive_int(payload['xtb_multiplicity'], 'xTB multiplicity')}",
                    f"accuracy = {self.parse_float(payload['xtb_accuracy'], 'xTB accuracy')}",
                    f"electronic_temperature = {self.parse_float(payload['xtb_electronic_temperature'], 'xTB electronic temperature')}",
                    f"max_iterations = {self.parse_positive_int(payload['xtb_max_iterations'], 'xTB max iterations')}",
                    f"omp_num_threads = {self.parse_positive_int(payload.get('xtb_omp_num_threads', '1'), 'xTB OMP threads')}",
                    f"solvent = {self.toml_string(payload.get('xtb_solvent', 'none'))}",
                    f"cache_api = {str(self.parse_bool(payload.get('xtb_cache_api', 'true'))).lower()}",
                    f"use_unwrapped_positions = {str(self.parse_bool(payload.get('xtb_use_unwrapped_positions', 'true'))).lower()}",
                ]
            )
        elif force_type == "gaussian":
            command = payload["gaussian_command_custom"].strip() if payload["gaussian_command"] == "Custom" else payload["gaussian_command"]
            if not command:
                raise ValueError("Gaussian command is required.")
            route = payload["gaussian_route_custom"].strip() if payload["gaussian_route"] == "Custom" else payload["gaussian_route"]
            if not route:
                raise ValueError("Gaussian route section is required.")
            if not route.lstrip().startswith("#"):
                route = "# " + route.strip()

            lines.extend(
                [
                    f"command = {self.toml_string(command)}",
                    f"route = {self.toml_string(route)}",
                    f"charge = {self.parse_int(payload['gaussian_charge'], 'Gaussian charge')}",
                    f"multiplicity = {self.parse_positive_int(payload['gaussian_multiplicity'], 'Gaussian multiplicity')}",
                ]
            )
            nproc = self.parse_positive_int(payload.get("gaussian_nproc", ""), "Gaussian nproc", required=False)
            if nproc is not None:
                lines.append(f"nproc = {nproc}")
            memory = (payload.get("gaussian_memory", "") or "").strip()
            if memory:
                lines.append(f"memory = {self.toml_string(memory)}")
            chk_file = (payload.get("gaussian_chk", "") or "").strip()
            if chk_file:
                lines.append(f"chk = {self.toml_string(chk_file)}")
            workdir = (payload.get("gaussian_workdir", "") or "").strip() or "gaussian_steps"
            lines.append(f"workdir = {self.toml_string(workdir)}")
        elif force_type == "classical":
            cutoff = self.parse_float(payload.get("classical_cutoff_angstrom", ""), "Classical cutoff", required=False)
            if cutoff is not None:
                lines.append(f"cutoff_angstrom = {cutoff}")
            lines.append(f"exclude_bonded = {str(self.parse_bool(payload.get('classical_exclude_bonded', 'true'))).lower()}")
            atom_types = self.parse_atom_types(payload.get("classical_atom_types", ""))
            if atom_types:
                lines.append("atom_types = [" + ", ".join(self.toml_string(item) for item in atom_types) + "]")
            bonds_text = (payload.get("classical_bonds_text", "") or "").strip()
            if bonds_text:
                lines.append("")
                lines.extend(bonds_text.splitlines())
            lj_text = (payload.get("classical_lj_text", "") or "").strip()
            if lj_text:
                force_tail.extend(["", lj_text])

        trajectory = (payload.get("trajectory", "") or "").strip() or "TRAJEC.xyz"
        log = (payload.get("log", "") or "").strip() or self.default_log_name(xyz_path)

        lines.extend(
            [
                "",
                *force_tail,
                "",
                "[output]",
                f"trajectory = {self.toml_string(trajectory)}",
                f"log = {self.toml_string(log)}",
                f"log_interval = {self.parse_positive_int(payload['log_interval'], 'Log interval')}",
            ]
        )

        if self.parse_bool(payload.get("include_restart", "true")):
            lines.extend(
                [
                    "",
                    "[restart]",
                    f"path = {self.toml_string(payload['restart_path'])}",
                    f"interval = {self.parse_positive_int(payload['restart_interval'], 'Restart interval')}",
                    f"resume_from_RESTART = {str(self.parse_bool(payload['restart_resume_from_RESTART'])).lower()}",
                    f"resume_from_GEOMETRY = {str(self.parse_bool(payload['restart_resume_from_GEOMETRY'])).lower()}",
                ]
            )

        return "\n".join(lines).replace("\n\n\n", "\n\n").rstrip() + "\n"

    def save_toml(self, payload: Dict[str, Any], toml_text: str, filename: Optional[str] = None) -> Path:
        xyz_path = Path(payload["xyz"])
        output_name = (filename or "").strip() or f"{xyz_path.stem}_gqteaMD.toml"
        out_path = xyz_path.with_name(output_name)
        out_path.write_text(toml_text, encoding="utf-8", newline="\n")
        return out_path


class GqteaMDInputBuilderUI:
    """Toga GUI for building gqteaMD TOML input files."""

    def __init__(self, *args, app: Optional[toga.App] = None):
        del args
        resolved_app = app if isinstance(app, toga.App) else None
        self.app = resolved_app
        self.formal_name = getattr(resolved_app, "formal_name", "gqteaMD Input Builder")
        self.builder = GqteaMDInputBuilder()
        self.preview_window: Optional[toga.Window] = None
        self.help_window: Optional[toga.Window] = None

        if self.app is None:
            self.main_window = toga.Window(title=self.formal_name, size=(600, 700))
        else:
            self.main_window = toga.MainWindow(title=self.formal_name)

        self.action_buttons: list[toga.Button] = []
        self._build_ui()
        self.main_window.show()

    def _build_ui(self):
        content = toga.Box(style=Pack(direction=COLUMN, margin=10))

        top_actions = toga.Box(style=Pack(direction=ROW, margin_bottom=10))
        load_btn = toga.Button("Open XYZ", on_press=self.open_xyz_file, style=Pack(margin_right=8))
        preview_btn = toga.Button("Preview TOML", on_press=self.preview_toml, style=Pack(margin_right=8))
        save_btn = toga.Button("Save TOML", on_press=self.save_toml, style=Pack(margin_right=8))
        clear_btn = toga.Button("Clear form", on_press=self.clear_form, style=Pack(margin_right=8))
        help_btn = toga.Button("Manual notes", on_press=self.show_manual_notes)
        self.action_buttons = [load_btn, preview_btn, save_btn, clear_btn, help_btn]
        for button in self.action_buttons:
            top_actions.add(button)

        self.xyz_status_label = toga.Label("No XYZ file loaded.", style=Pack(margin_bottom=8))
        self.output_name_input = toga.TextInput(
            value="",
            placeholder="Leave blank to use <xyzname>_gqteaMD.toml",
            style=Pack(width=300, margin_left=8),
        )
        output_row = toga.Box(style=Pack(direction=ROW, margin_bottom=10))
        output_row.add(toga.Label("Output filename:", style=Pack(width=130, margin_top=6)))
        output_row.add(self.output_name_input)

        self.scroll = toga.ScrollContainer(style=Pack(flex=1))
        scroll_content = toga.Box(style=Pack(direction=COLUMN, margin=5))
        scroll_content.add(self.make_input_box())
        scroll_content.add(self.make_cell_dynamics_box())
        scroll_content.add(self.make_force_provider_box())
        scroll_content.add(self.make_output_box())
        scroll_content.add(self.make_restart_box())
        self.scroll.content = scroll_content

        content.add(top_actions)
        content.add(self.xyz_status_label)
        content.add(output_row)
        content.add(self.scroll)
        self.main_window.content = content
        self.rebuild_force_options()

    def make_input_box(self) -> toga.Box:
        outer = self.make_section("Input geometry")
        row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        row.add(toga.Label("XYZ path", style=Pack(width=130, margin_top=6)))
        self.xyz_path_input = toga.TextInput(placeholder="initial.xyz", style=Pack(flex=1))
        row.add(self.xyz_path_input)
        outer.add(row)
        return outer

    def make_cell_dynamics_box(self) -> toga.Box:
        outer = self.make_section("Cell and dynamics")
        cell_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        cell_row.add(toga.Label("Cell a, b, c", style=Pack(width=130, margin_top=6)))
        self.cell_a_input = toga.TextInput(value="20.0", style=Pack(width=90, margin_right=8))
        self.cell_b_input = toga.TextInput(value="20.0", style=Pack(width=90, margin_right=8))
        self.cell_c_input = toga.TextInput(value="20.0", style=Pack(width=90))
        cell_row.add(self.cell_a_input)
        cell_row.add(self.cell_b_input)
        cell_row.add(self.cell_c_input)
        outer.add(cell_row)

        dyn_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        dyn_row.add(toga.Label("Timestep fs", style=Pack(width=130, margin_top=6)))
        self.timestep_input = toga.TextInput(value="0.5", style=Pack(width=100, margin_right=18))
        dyn_row.add(self.timestep_input)
        dyn_row.add(toga.Label("Steps", style=Pack(width=60, margin_top=6)))
        self.steps_input = toga.TextInput(value="100", style=Pack(width=100))
        dyn_row.add(self.steps_input)
        outer.add(dyn_row)
        return outer

    def make_force_provider_box(self) -> toga.Box:
        outer = self.make_section("Force provider")
        row = toga.Box(style=Pack(direction=ROW, margin_bottom=8))
        row.add(toga.Label("Provider", style=Pack(width=130, margin_top=6)))
        self.force_provider_selection = toga.Selection(
            items=self.builder.FORCE_PROVIDER_CHOICES,
            value="harmonic",
            on_change=self.on_force_provider_change,
            style=Pack(width=180),
        )
        row.add(self.force_provider_selection)
        outer.add(row)

        self.force_options_container = toga.Box(style=Pack(direction=COLUMN))
        outer.add(self.force_options_container)
        return outer

    def make_output_box(self) -> toga.Box:
        outer = self.make_section("Output")
        self.trajectory_input = self.make_text_row(outer, "Trajectory", "TRAJEC.xyz", "TRAJEC.xyz")
        self.log_input = self.make_text_row(outer, "Log", "<xyzname>_gqteaMD.log", "")
        outer.add(
            toga.Label(
                "GEOMETRY is written automatically at every calculation step with positions, velocities, and forces.",
                style=Pack(margin_bottom=6),
            )
        )
        self.log_interval_input = self.make_text_row(outer, "Log interval", placeholder="Write LOG every N calculation steps. Default value: 1",)
        return outer

    def make_restart_box(self) -> toga.Box:
        outer = self.make_section("Restart")
        row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        row.add(toga.Label("Include restart", style=Pack(width=130, margin_top=6)))
        self.include_restart_selection = toga.Selection(items=self.builder.BOOL_CHOICES, value="true", style=Pack(width=80))
        row.add(self.include_restart_selection)
        outer.add(row)
        self.restart_path_input = self.make_text_row(outer, "Restart path", "RESTART", "RESTART")
        self.restart_interval_input = self.make_text_row(
            outer,
            "Interval",
            placeholder="Write RESTART every N calculation steps. Default value: 5",
        )

        restart_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        restart_row.add(toga.Label("resume_from_RESTART", style=Pack(width=170, margin_top=6)))
        self.restart_resume_selection = toga.Selection(items=self.builder.BOOL_CHOICES, value="false", style=Pack(width=80))
        restart_row.add(self.restart_resume_selection)
        outer.add(restart_row)

        geometry_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        geometry_row.add(toga.Label("resume_from_GEOMETRY", style=Pack(width=170, margin_top=6)))
        self.restart_outputs_selection = toga.Selection(items=self.builder.BOOL_CHOICES, value="false", style=Pack(width=80))
        geometry_row.add(self.restart_outputs_selection)
        outer.add(geometry_row)
        return outer

    def rebuild_force_options(self):
        self.force_options_container.clear()
        force_type = self.force_provider_selection.value or "harmonic"
        if force_type == "harmonic":
            self.harmonic_k_input = self.make_text_row(
                self.force_options_container,
                "k eV/angstrom2",
                "0.1",
                "0.1",
            )
        elif force_type == "uff":
            self.uff_bond_detection_scale_input = self.make_text_row(
                self.force_options_container,
                "Bond scale",
                "1.2",
                "1.2",
            )
            self.uff_cutoff_input = self.make_text_row(
                self.force_options_container,
                "Cutoff",
                "Optional, e.g. 10.0",
                "10.0",
            )
            self.uff_atom_types_input = self.make_text_row(
                self.force_options_container,
                "Atom types",
                "Optional comma-separated explicit UFF types",
                "",
            )
            self.uff_charges_input = self.make_text_row(
                self.force_options_container,
                "Charges",
                "Optional comma-separated charges, one per atom",
                "",
            )
            self.uff_electrostatics_selection = self.make_selection_row(
                self.force_options_container,
                "Electrostatics",
                self.builder.UFF_ELECTROSTATICS_CHOICES,
                "auto",
            )
            self.uff_nonbonded_exclusions_selection = self.make_selection_row(
                self.force_options_container,
                "Exclusions",
                self.builder.UFF_EXCLUSION_CHOICES,
                "exclude_12_13",
            )
            scale_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
            scale_row.add(toga.Label("1-4 scales", style=Pack(width=130, margin_top=6)))
            self.uff_lj_14_scale_input = toga.TextInput(value="1.0", placeholder="LJ", style=Pack(width=90, margin_right=8))
            self.uff_electrostatic_14_scale_input = toga.TextInput(
                value="1.0",
                placeholder="Coulomb",
                style=Pack(width=90, margin_right=18),
            )
            scale_row.add(self.uff_lj_14_scale_input)
            scale_row.add(self.uff_electrostatic_14_scale_input)
            scale_row.add(toga.Label("LJ cutoff", style=Pack(width=70, margin_top=6)))
            self.uff_lj_cutoff_mode_selection = toga.Selection(
                items=self.builder.UFF_LJ_CUTOFF_CHOICES,
                value="plain",
                style=Pack(width=100),
            )
            scale_row.add(self.uff_lj_cutoff_mode_selection)
            self.force_options_container.add(scale_row)
            neighbor_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
            neighbor_row.add(toga.Label("Neighbor list", style=Pack(width=130, margin_top=6)))
            self.uff_use_neighbor_list_selection = toga.Selection(
                items=self.builder.BOOL_CHOICES,
                value="true",
                style=Pack(width=90, margin_right=8),
            )
            self.uff_neighbor_skin_input = toga.TextInput(
                value="2.0",
                placeholder="Skin A",
                style=Pack(width=100),
            )
            neighbor_row.add(self.uff_use_neighbor_list_selection)
            neighbor_row.add(self.uff_neighbor_skin_input)
            self.force_options_container.add(neighbor_row)
            bond_orders_label = toga.Label("Bond-order TOML lines", style=Pack(margin_bottom=4))
            self.force_options_container.add(bond_orders_label)
            self.uff_bond_orders_text = toga.MultilineTextInput(style=Pack(height=80, margin_bottom=8))
            self.uff_bond_orders_text.placeholder = (
                "Optional example:\n"
                "bond_orders = [\n"
                "  { atoms = [0, 1], order = 2.0 },\n"
                "]"
            )
            self.force_options_container.add(self.uff_bond_orders_text)
            topology_label = toga.Label("Topology TOML lines", style=Pack(margin_bottom=4))
            self.force_options_container.add(topology_label)
            self.uff_topology_text = toga.MultilineTextInput(style=Pack(height=120, margin_bottom=8))
            self.uff_topology_text.placeholder = (
                "Optional examples:\n"
                "bonds = [{ atoms = [0, 1], order = 1.0 }]\n"
                "angles = [{ atoms = [1, 0, 2] }]\n"
                "torsions = []\n"
                "inversions = []"
            )
            self.force_options_container.add(self.uff_topology_text)
        elif force_type == "gaussian":
            self.gaussian_command_selection = self.make_selection_row(
                self.force_options_container,
                "Command",
                self.builder.GAUSSIAN_COMMAND_CHOICES,
                "g16",
                self.on_gaussian_choice_change,
            )
            self.gaussian_command_custom_input = self.make_text_row(
                self.force_options_container,
                "Custom command",
                "Full path or executable name",
                "",
            )
            self.gaussian_route_selection = self.make_selection_row(
                self.force_options_container,
                "Route",
                self.builder.ROUTE_CHOICES,
                "# B3LYP/6-31G(d) Force NoSymm SCF=Tight",
                self.on_gaussian_choice_change,
            )
            self.gaussian_route_custom_input = self.make_text_row(
                self.force_options_container,
                "Custom route",
                "# method/basis Force NoSymm SCF=Tight",
                "",
            )
            charge_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
            charge_row.add(toga.Label("Charge/mult", style=Pack(width=130, margin_top=6)))
            self.gaussian_charge_input = toga.TextInput(value="0", style=Pack(width=80, margin_right=8))
            self.gaussian_multiplicity_input = toga.TextInput(value="1", style=Pack(width=80, margin_right=18))
            charge_row.add(self.gaussian_charge_input)
            charge_row.add(self.gaussian_multiplicity_input)
            charge_row.add(toga.Label("nproc", style=Pack(width=50, margin_top=6)))
            self.gaussian_nproc_input = toga.TextInput(value="4", style=Pack(width=80))
            charge_row.add(self.gaussian_nproc_input)
            self.force_options_container.add(charge_row)
            self.gaussian_workdir_input = self.make_text_row(
                self.force_options_container,
                "Workdir",
                "gaussian_steps",
                "gaussian_steps",
            )
            self.gaussian_chk_input = self.make_text_row(
                self.force_options_container,
                "chk file",
                "Optional note, e.g. step.chk",
                "",
            )
            self.gaussian_memory_input = self.make_text_row(
                self.force_options_container,
                "Memory",
                "Optional, e.g. 1500MB for 32-bit Gaussian or 4GB for 64-bit Gaussian",
                "",
            )
            self.on_gaussian_choice_change(None)
        elif force_type == "xtb":
            self.xtb_method_selection = self.make_selection_row(
                self.force_options_container,
                "Method",
                self.builder.XTB_METHOD_CHOICES,
                "GFN2-xTB",
                self.on_xtb_choice_change,
            )
            self.xtb_method_custom_input = self.make_text_row(
                self.force_options_container,
                "Custom method",
                "Custom xTB method name",
                "",
            )
            self.xtb_command_input = self.make_text_row(
                self.force_options_container,
                "Command",
                "xtb executable name or full path",
                "xtb",
            )
            xtb_charge_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
            xtb_charge_row.add(toga.Label("Charge/mult", style=Pack(width=130, margin_top=6)))
            self.xtb_charge_input = toga.TextInput(value="0", style=Pack(width=80, margin_right=8))
            self.xtb_multiplicity_input = toga.TextInput(value="1", style=Pack(width=80, margin_right=18))
            xtb_charge_row.add(self.xtb_charge_input)
            xtb_charge_row.add(self.xtb_multiplicity_input)
            xtb_charge_row.add(toga.Label("max iter", style=Pack(width=70, margin_top=6)))
            self.xtb_max_iterations_input = toga.TextInput(value="250", style=Pack(width=80))
            xtb_charge_row.add(self.xtb_max_iterations_input)
            self.force_options_container.add(xtb_charge_row)
            xtb_numeric_row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
            xtb_numeric_row.add(toga.Label("Acc/temp K", style=Pack(width=130, margin_top=6)))
            self.xtb_accuracy_input = toga.TextInput(value="1.0", style=Pack(width=90, margin_right=8))
            self.xtb_electronic_temperature_input = toga.TextInput(value="300.0", style=Pack(width=100, margin_right=18))
            xtb_numeric_row.add(self.xtb_accuracy_input)
            xtb_numeric_row.add(self.xtb_electronic_temperature_input)
            xtb_numeric_row.add(toga.Label("OMP threads", style=Pack(width=90, margin_top=6)))
            self.xtb_omp_num_threads_input = toga.TextInput(value="1", style=Pack(width=70))
            xtb_numeric_row.add(self.xtb_omp_num_threads_input)
            self.force_options_container.add(xtb_numeric_row)
            self.xtb_solvent_input = self.make_text_row(
                self.force_options_container,
                "Solvent",
                "none, water, methanol, ...",
                "none",
            )
            self.xtb_cache_api_selection = self.make_selection_row(
                self.force_options_container,
                "Cache API",
                self.builder.BOOL_CHOICES,
                "true",
            )
            self.xtb_use_unwrapped_positions_selection = self.make_selection_row(
                self.force_options_container,
                "Unwrapped pos",
                self.builder.BOOL_CHOICES,
                "true",
            )
            self.force_options_container.add(
                toga.Label(
                    "xTB provides potential energy and atomic forces for gqteaMD propagation.",
                    style=Pack(margin_bottom=6),
                )
            )
            self.on_xtb_choice_change(None)
        elif force_type == "classical":
            self.classical_cutoff_input = self.make_text_row(
                self.force_options_container,
                "Cutoff",
                "Optional, e.g. 10.0",
                "10.0",
            )
            self.classical_exclude_bonded_selection = self.make_selection_row(
                self.force_options_container,
                "Exclude bonded",
                self.builder.BOOL_CHOICES,
                "true",
            )
            self.classical_atom_types_input = self.make_text_row(
                self.force_options_container,
                "Atom types",
                "Optional comma-separated atom types",
                "",
            )
            bonds_label = toga.Label("Bonds TOML lines", style=Pack(margin_bottom=4))
            self.force_options_container.add(bonds_label)
            self.classical_bonds_text = toga.MultilineTextInput(style=Pack(height=95, margin_bottom=8))
            self.classical_bonds_text.placeholder = (
                "Example:\n"
                "bonds = [\n"
                "  { atoms = [0, 1], k_ev_per_angstrom2 = 45.0, r0_angstrom = 0.9572 },\n"
                "]"
            )
            self.force_options_container.add(self.classical_bonds_text)
            lj_label = toga.Label("Lennard-Jones tables", style=Pack(margin_bottom=4))
            self.force_options_container.add(lj_label)
            self.classical_lj_text = toga.MultilineTextInput(style=Pack(height=120, margin_bottom=8))
            self.classical_lj_text.placeholder = (
                "[force_provider.lennard_jones.O]\n"
                "epsilon_ev = 0.0067\n"
                "sigma_angstrom = 3.1507"
            )
            self.force_options_container.add(self.classical_lj_text)

    def make_section(self, title: str) -> toga.Box:
        outer = toga.Box(style=Pack(direction=COLUMN, margin_bottom=14))
        outer.add(toga.Label(title, style=Pack(font_weight="bold", margin_bottom=6)))
        return outer

    def make_text_row(self, parent: toga.Box, label_text: str, placeholder: str, value: str = "") -> toga.TextInput:
        row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        row.add(toga.Label(label_text, style=Pack(width=130, margin_top=6)))
        widget = toga.TextInput(value=value, placeholder=placeholder, style=Pack(flex=1))
        row.add(widget)
        parent.add(row)
        return widget

    def make_selection_row(
        self,
        parent: toga.Box,
        label_text: str,
        items: list[str],
        value: str,
        on_change=None,
    ) -> toga.Selection:
        row = toga.Box(style=Pack(direction=ROW, margin_bottom=6))
        row.add(toga.Label(label_text, style=Pack(width=130, margin_top=6)))
        widget = toga.Selection(items=items, value=value, on_change=on_change, style=Pack(width=280))
        row.add(widget)
        parent.add(row)
        return widget

    def set_action_buttons_enabled(self, enabled: bool):
        for button in self.action_buttons:
            button.enabled = enabled

    def on_force_provider_change(self, widget):
        del widget
        self.rebuild_force_options()

    def on_gaussian_choice_change(self, widget):
        del widget
        command_is_custom = (self.gaussian_command_selection.value or "") == "Custom"
        route_is_custom = (self.gaussian_route_selection.value or "") == "Custom"
        self.gaussian_command_custom_input.enabled = command_is_custom
        self.gaussian_route_custom_input.enabled = route_is_custom
        if not command_is_custom:
            self.gaussian_command_custom_input.value = ""
        if not route_is_custom:
            self.gaussian_route_custom_input.value = ""

    def on_xtb_choice_change(self, widget):
        del widget
        method_is_custom = (self.xtb_method_selection.value or "") == "Custom"
        self.xtb_method_custom_input.enabled = method_is_custom
        if not method_is_custom:
            self.xtb_method_custom_input.value = ""

    async def open_xyz_file(self, widget):
        del widget
        self.set_action_buttons_enabled(False)
        try:
            file = await self.main_window.dialog(
                toga.OpenFileDialog(title="Select gqteaMD starting XYZ file", file_types=["xyz"])
            )
            if not file:
                return
            self.xyz_path_input.value = str(file)
            self.log_input.value = self.builder.default_log_name(str(file))
            atom_count, comment = self.builder.inspect_xyz(str(file))
            self.xyz_status_label.text = f"Loaded XYZ: {file} | Atoms: {atom_count} | Comment: {comment}"
        except Exception as exc:
            await self.main_window.dialog(toga.ErrorDialog(title="XYZ read error", message=str(exc)))
        finally:
            self.set_action_buttons_enabled(True)

    async def clear_form(self, widget):
        del widget
        self.xyz_path_input.value = ""
        self.xyz_status_label.text = "No XYZ file loaded."
        self.output_name_input.value = ""
        self.cell_a_input.value = "20.0"
        self.cell_b_input.value = "20.0"
        self.cell_c_input.value = "20.0"
        self.timestep_input.value = "0.5"
        self.steps_input.value = "100"
        self.force_provider_selection.value = "harmonic"
        self.rebuild_force_options()
        self.trajectory_input.value = "TRAJEC.xyz"
        self.log_input.value = ""
        self.log_interval_input.value = "1"
        self.include_restart_selection.value = "true"
        self.restart_path_input.value = "RESTART"
        self.restart_interval_input.value = ""
        self.restart_resume_selection.value = "false"
        self.restart_outputs_selection.value = "false"

    def collect_payload(self) -> Dict[str, Any]:
        force_type = self.force_provider_selection.value or "harmonic"
        payload: Dict[str, Any] = {
            "xyz": self.xyz_path_input.value or "",
            "cell_a": self.cell_a_input.value or "",
            "cell_b": self.cell_b_input.value or "",
            "cell_c": self.cell_c_input.value or "",
            "timestep_fs": self.timestep_input.value or "",
            "steps": self.steps_input.value or "",
            "force_type": force_type,
            "trajectory": self.trajectory_input.value or "TRAJEC.xyz",
            "log": self.log_input.value or "",
            "log_interval": self.log_interval_input.value or "1",
            "include_restart": self.include_restart_selection.value or "false",
            "restart_path": self.restart_path_input.value or "RESTART",
            "restart_interval": self.restart_interval_input.value or "5",
            "restart_resume_from_RESTART": self.restart_resume_selection.value or "false",
            "restart_resume_from_GEOMETRY": self.restart_outputs_selection.value or "false",
        }
        if force_type == "harmonic":
            payload["harmonic_k"] = self.harmonic_k_input.value or "0.1"
        elif force_type == "uff":
            payload["uff_bond_detection_scale"] = self.uff_bond_detection_scale_input.value or "1.2"
            payload["uff_cutoff_angstrom"] = self.uff_cutoff_input.value or ""
            payload["uff_atom_types"] = self.uff_atom_types_input.value or ""
            payload["uff_charges"] = self.uff_charges_input.value or ""
            payload["uff_electrostatics"] = self.uff_electrostatics_selection.value or "auto"
            payload["uff_nonbonded_exclusions"] = self.uff_nonbonded_exclusions_selection.value or "exclude_12_13"
            payload["uff_lj_14_scale"] = self.uff_lj_14_scale_input.value or "1.0"
            payload["uff_electrostatic_14_scale"] = self.uff_electrostatic_14_scale_input.value or "1.0"
            payload["uff_lj_cutoff_mode"] = self.uff_lj_cutoff_mode_selection.value or "plain"
            payload["uff_use_neighbor_list"] = self.uff_use_neighbor_list_selection.value or "true"
            payload["uff_neighbor_skin_angstrom"] = self.uff_neighbor_skin_input.value or "2.0"
            payload["uff_bond_orders_text"] = self.uff_bond_orders_text.value or ""
            payload["uff_topology_text"] = self.uff_topology_text.value or ""
        elif force_type == "xtb":
            payload["xtb_method"] = self.xtb_method_selection.value or "GFN2-xTB"
            payload["xtb_method_custom"] = self.xtb_method_custom_input.value or ""
            payload["xtb_command"] = self.xtb_command_input.value or ""
            payload["xtb_charge"] = self.xtb_charge_input.value or "0"
            payload["xtb_multiplicity"] = self.xtb_multiplicity_input.value or "1"
            payload["xtb_accuracy"] = self.xtb_accuracy_input.value or "1.0"
            payload["xtb_electronic_temperature"] = self.xtb_electronic_temperature_input.value or "300.0"
            payload["xtb_max_iterations"] = self.xtb_max_iterations_input.value or "250"
            payload["xtb_omp_num_threads"] = self.xtb_omp_num_threads_input.value or "1"
            payload["xtb_solvent"] = self.xtb_solvent_input.value or "none"
            payload["xtb_cache_api"] = self.xtb_cache_api_selection.value or "true"
            payload["xtb_use_unwrapped_positions"] = self.xtb_use_unwrapped_positions_selection.value or "true"
        elif force_type == "gaussian":
            payload["gaussian_command"] = self.gaussian_command_selection.value or "g16"
            payload["gaussian_command_custom"] = self.gaussian_command_custom_input.value or ""
            payload["gaussian_route"] = self.gaussian_route_selection.value or "# B3LYP/6-31G(d) Force NoSymm SCF=Tight"
            payload["gaussian_route_custom"] = self.gaussian_route_custom_input.value or ""
            payload["gaussian_charge"] = self.gaussian_charge_input.value or "0"
            payload["gaussian_multiplicity"] = self.gaussian_multiplicity_input.value or "1"
            payload["gaussian_nproc"] = self.gaussian_nproc_input.value or ""
            payload["gaussian_workdir"] = self.gaussian_workdir_input.value or "gaussian_steps"
            payload["gaussian_chk"] = self.gaussian_chk_input.value or ""
            payload["gaussian_memory"] = self.gaussian_memory_input.value or ""
        elif force_type == "classical":
            payload["classical_cutoff_angstrom"] = self.classical_cutoff_input.value or ""
            payload["classical_exclude_bonded"] = self.classical_exclude_bonded_selection.value or "true"
            payload["classical_atom_types"] = self.classical_atom_types_input.value or ""
            payload["classical_bonds_text"] = self.classical_bonds_text.value or ""
            payload["classical_lj_text"] = self.classical_lj_text.value or ""
        return payload

    async def preview_toml(self, widget):
        del widget
        try:
            payload = self.collect_payload()
            text = self.builder.generate_toml(payload)
            preview_window = toga.Window(title="gqteaMD TOML preview", size=(900, 720))
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

    async def save_toml(self, widget):
        del widget
        try:
            payload = self.collect_payload()
            text = self.builder.generate_toml(payload)
            saved_path = self.builder.save_toml(payload, text, self.output_name_input.value or None)
            await self.main_window.dialog(
                toga.InfoDialog(title="TOML saved", message=f"gqteaMD input saved to:\n{saved_path}")
            )
        except Exception as exc:
            await self.main_window.dialog(toga.ErrorDialog(title="Save error", message=str(exc)))

    async def show_manual_notes(self, widget):
        del widget
        if self.help_window is not None and not getattr(self.help_window, "closed", False):
            self.help_window.show()
            return
        self.help_window = toga.Window(title="gqteaMD input builder notes", size=(820, 560))
        help_box = toga.Box(style=Pack(direction=COLUMN, margin=10, flex=1))
        help_text = toga.MultilineTextInput(
            value=self.builder.MANUAL_SUMMARY,
            readonly=True,
            style=Pack(flex=1),
        )
        help_box.add(help_text)
        self.help_window.content = help_box
        self.help_window.show()


class GqteaMDInputBuilderApp(toga.App):
    """Standalone wrapper app for the gqteaMD input builder."""

    def startup(self):
        self.ui = GqteaMDInputBuilderUI(app=self)
        self.main_window = self.ui.main_window


def main():
    return GqteaMDInputBuilderApp("gqteaMD Input Builder", "br.ueg.gqtea.gqteamdinput")


if __name__ == "__main__":
    app = main()
    app.main_loop()
