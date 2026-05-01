# Read CPMD GEOMETRY.xyz and build collision input files
# Module refactored and improved by A. J. Camargo, 11/09/2025
import pathlib
import datetime
from math import sqrt

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT, CENTER

from help import AtomicData, HelpGqteaWin
from framesCounter import FramesCounter

BOHR_PER_ANGSTROM = 1.0 / 0.52917720859  # Å -> bohr
AU_VELOCITY_TO_M_PER_S = 2.18769126379e6  # a.u. -> m/s
AU_VELOCITY_TO_ANGSTROM_PER_FS = AU_VELOCITY_TO_M_PER_S * 1.0e-5
MD_ENGINE_CHOICES = ["cpmd", "gqteaMD"]


class CollisionSetUp(FramesCounter):
    """
    Build CPMD restart GEOMETRY/velocity for collision MD.
    Reads a CPMD one-step GEOMETRY.xyz (Å + a.u. vel) and creates:
      - output_dir/newGeometry.xyz (echo of parsed data, Å + a.u.)
      - output_dir/with-vibration-GEOMETRY (bohr + a.u., CPMD GEOFILE body)
      - output_dir/without-vibration-GEOMETRY (optional, zeros non-attackers)
      - output_dir/summary.txt (run summary; per-atom & optional COM energies)
    """

    # -------- File I/O helpers --------
    def _ensure_output_dir(self):
        if not hasattr(self, "output_dir") or not self.output_dir:
            base = pathlib.Path(getattr(self, "trajec", ".")).resolve().parent
            self.output_dir = str(base / "collision-output")
        pathlib.Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    # -------- Data loading --------
    async def readGeometry(self, widget):
        """Read CPMD GEOMETRY.xyz or gqteaMD GEOMETRY data."""
        try:
            await self.open_file_dialog(widget)  # expected to set self.trajec
        except Exception as e:
            await self.main_window.dialog(
                toga.InfoDialog("File chooser error", str(e))
            )
            return

        if not getattr(self, "trajec", None):
            await self.main_window.dialog(
                toga.InfoDialog("No file", "No file selected.")
            )
            return

        if hasattr(self, "textInput_file"):
            self.textInput_file.value = self.trajec

        self.coords = []
        try:
            with open(self.trajec, "r", encoding="utf-8") as f:
                lines = f.readlines()
        except FileNotFoundError:
            await self.main_window.dialog(toga.InfoDialog("Error", "File not found."))
            return
        except Exception as e:
            await self.main_window.dialog(
                toga.InfoDialog("Error", f"Could not read file: {e}")
            )
            return

        if not lines:
            await self.main_window.dialog(
                toga.InfoDialog("Warning", "The selected file is empty.")
            )
            return

        self.multi_line_text.value = "".join(lines)

        coords = []
        try:
            self.num_atoms = int(lines[0].strip())
            data_start = 2
            data_lines = lines[data_start:data_start + self.num_atoms]
            if len(data_lines) < self.num_atoms:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", "File does not contain enough atom data.")
                )
                return
        except ValueError:
            data_start = 0
            data_lines = [line for line in lines if line.split() and not line.lstrip().startswith("#")]
            self.num_atoms = len(data_lines)

        for offset, line in enumerate(data_lines):
            i = data_start + offset
            parts = line.split()
            if len(parts) < 7:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Line {i+1} is incomplete.")
                )
                return
            atom = parts[0]
            try:
                x, y, z, vx, vy, vz = map(float, parts[1:7])
                if len(parts) >= 10:
                    fx, fy, fz = map(float, parts[7:10])
                else:
                    fx, fy, fz = 0.0, 0.0, 0.0
            except ValueError:
                await self.main_window.dialog(
                    toga.InfoDialog("Error", f"Invalid numeric data at line {i+1}.")
                )
                return
            coords.append([atom, x, y, z, vx, vy, vz, fx, fy, fz])

        self.raw_spv = coords
        self.spv = [row[:] for row in coords]
        self.num_atoms = len(coords)
        self._ensure_output_dir()

    # -------- Pipeline --------
    async def execute(self, widget):
        """Run the full pipeline after user inputs are set."""
        if not hasattr(self, "raw_spv"):
            await self.main_window.dialog(
                toga.InfoDialog("Missing file", "Load a GEOMETRY.xyz first (Browse).")
            )
            return

        try:
            idx_str = self.textInput_atomic_labels.value.strip()
            if not idx_str:
                raise ValueError("Attacker atom indices are required (1-based).")
            attacker_idx = [int(i) for i in idx_str.split()]
            for i in attacker_idx:
                if i < 1 or i > self.num_atoms:
                    raise ValueError(
                        f"Attacker atom index {i} out of range [1, {self.num_atoms}]."
                    )
            self.attacker_mol_idx = attacker_idx

            target_pos_list = self.textInput_coords_site.value.strip().split()
            if len(target_pos_list) != 3:
                raise ValueError("Target position must have exactly 3 numbers (x y z).")
            self.target_pos = [float(i) for i in target_pos_list]

            self.md_engine = self.get_selected_md_engine()

            init_vel_str = self.textInput_initial_velocity.value.strip()
            if not init_vel_str:
                raise ValueError(f"Initial attacker velocity is required ({self.initial_velocity_unit_label()}).")
            self.init_attacker_vel_input = float(init_vel_str)
            if self.md_engine == "gqteaMD":
                self.init_attacker_vel = self.init_attacker_vel_input / AU_VELOCITY_TO_ANGSTROM_PER_FS
            else:
                self.init_attacker_vel = self.init_attacker_vel_input

        except ValueError as ve:
            await self.main_window.dialog(toga.InfoDialog("Input error", str(ve)))
            return
        except Exception as e:
            await self.main_window.dialog(
                toga.InfoDialog("Error", f"An error occurred: {e}")
            )
            return

        self.multi_line_text.value = ""
        self.multi_line_text.value += (
            f"Attacker atom indices (1-based):    {self.attacker_mol_idx}\n"
        )
        self.multi_line_text.value += f"Target position (Å):                {self.target_pos}\n"
        if self.md_engine == "gqteaMD":
            self.multi_line_text.value += (
                f"Initial velocity (Å/fs):            {self.init_attacker_vel_input}\n\n"
            )
        else:
            self.multi_line_text.value += (
                f"Initial velocity (a.u.):            {self.init_attacker_vel}\n\n"
            )
        self.multi_line_text.value += f"MD engine:                          {self.md_engine}\n\n"

        # Compute and write artifacts
        self.prepare_geometry_for_engine()
        self.setInitialVelocity()
        if self.md_engine == "cpmd":
            self.write_atomic_data()
        self.write_geometry_for_engine()
        self.kineticEnergy()     # now includes per-atom energies (+ optional COM)
        self.set_zero_velocity() # prints the with/without vibration content to UI when requested

        # Show summary
        try:
            with open(f"{self.output_dir}/summary.txt", "r", encoding="utf-8") as f:
                lines = f.readlines()
            self.multi_line_text.value += "\n>>> Summary <<<\n"
            self.multi_line_text.value += "".join(lines)
        except Exception as e:
            await self.main_window.dialog(
                toga.InfoDialog("Error", f"Could not read summary.txt: {e}")
            )

    # -------- Writers / calculators --------
    def write_atomic_data(self):
        """Write parsed/updated atomic data back out (Å + a.u.)."""
        self._ensure_output_dir()
        try:
            with open(f"{self.output_dir}/newGeometry.xyz", "w", encoding="utf-8") as f:
                f.write(f"{self.num_atoms}\n\n")
                for sym, x, y, z, vx, vy, vz, *_ in self.spv:
                    f.write(
                        "{:<2}{:>20.12f}{:>20.12f}{:>20.12f}{:>20.12f}{:>20.12f}{:>20.12f}\n"
                        .format(sym, x, y, z, vx, vy, vz)
                    )
        except Exception as e:
            self.multi_line_text.value += f"\nError writing atomic data: {e}\n"

    def get_selected_md_engine(self):
        engine = getattr(getattr(self, "md_engine_selection", None), "value", "cpmd")
        if engine not in MD_ENGINE_CHOICES:
            raise ValueError(f"Unsupported MD engine: {engine}")
        return engine

    def prepare_geometry_for_engine(self):
        """
        Reset the working geometry from the loaded file and normalize velocity
        units to a.u. for internal calculations.
        """
        self.spv = [row[:] for row in self.raw_spv]
        if self.md_engine == "gqteaMD":
            for row in self.spv:
                row[4] /= AU_VELOCITY_TO_ANGSTROM_PER_FS
                row[5] /= AU_VELOCITY_TO_ANGSTROM_PER_FS
                row[6] /= AU_VELOCITY_TO_ANGSTROM_PER_FS

    def initial_velocity_unit_label(self):
        if getattr(self, "md_engine", None) == "gqteaMD" or self.get_selected_md_engine() == "gqteaMD":
            return "Å/fs"
        return "a.u."

    def write_geometry_for_engine(self):
        if self.md_engine == "cpmd":
            self.GEOMETRY()
        elif self.md_engine == "gqteaMD":
            self.write_gqteamd_geometry()
        else:
            self.multi_line_text.value += f"\nUnsupported MD engine: {self.md_engine}\n"

    def GEOMETRY(self):
        """
        Create the CPMD GEOFILE body (bohr + a.u.) as 'with-vibration-GEOMETRY'.
        Does NOT mutate self.spv and does NOT print to UI here.
        """
        self._ensure_output_dir()
        try:
            vib_path = pathlib.Path(self.output_dir) / "with-vibration-GEOMETRY"
            with open(vib_path, "w", encoding="utf-8") as f_vib:
                for _, x, y, z, vx, vy, vz, *_ in self.spv:
                    xb, yb, zb = x * BOHR_PER_ANGSTROM, y * BOHR_PER_ANGSTROM, z * BOHR_PER_ANGSTROM
                    f_vib.write(
                        "{:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f}\n"
                        .format(xb, yb, zb, vx, vy, vz)
                    )
        except Exception as e:
            self.multi_line_text.value += f"\nError writing GEOMETRY file: {e}\n"

    def write_gqteamd_geometry(self):
        """
        Create a gqteaMD GEOMETRY file with symbols, Angstrom positions,
        Angstrom/fs velocities, and force columns.
        """
        self._ensure_output_dir()
        try:
            geometry_path = pathlib.Path(self.output_dir) / "gqteaMD-GEOMETRY"
            with open(geometry_path, "w", encoding="utf-8") as f_geom:
                f_geom.write(f"{self.num_atoms}\n")
                f_geom.write(
                    "Generated by gqteaWinToga collision setup; cartesian Angstrom positions, Angstrom/fs velocities\n"
                )
                for sym, x, y, z, vx, vy, vz, fx, fy, fz in self.spv:
                    vx_afs = vx * AU_VELOCITY_TO_ANGSTROM_PER_FS
                    vy_afs = vy * AU_VELOCITY_TO_ANGSTROM_PER_FS
                    vz_afs = vz * AU_VELOCITY_TO_ANGSTROM_PER_FS
                    f_geom.write(
                        "{:<2} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f}\n"
                        .format(sym, x, y, z, vx_afs, vy_afs, vz_afs, fx, fy, fz)
                    )
            self.multi_line_text.value += (
                "\n>>> gqteaMD GEOMETRY written to gqteaMD-GEOMETRY <<<\n"
                "Rename this file to GEOMETRY and set resume_from_GEOMETRY = true.\n"
                "Velocities are written in Angstrom/fs for gqteaMD.\n"
            )
        except Exception as e:
            self.multi_line_text.value += f"\nError writing gqteaMD GEOMETRY file: {e}\n"

    def kineticEnergy(self):
        """
        Compute kinetic energies:
          - Per-atom energies (kJ/mol and kcal/mol) for attacker atoms
          - Molar mass of attacker set
          - Optional center-of-mass kinetic energy for the attacker set (toggle)
        Uses molar masses (kg/mol) and velocities in m/s:
          E_molar [J/mol] = 1/2 * M [kg/mol] * v^2
        """
        self._ensure_output_dir()
        try:
            # Per-atom energies for attacker atoms
            lines = []
            total_M_kg_per_mol = 0.0

            for idx1 in self.attacker_mol_idx:
                i = idx1 - 1
                sym, _, _, _, vx_au, vy_au, vz_au, *_ = self.spv[i]
                try:
                    m_g_per_mol = AtomicData.atomic_masses[sym]
                except KeyError:
                    raise ValueError(f"Unknown atomic mass for symbol '{sym}'")
                m_kg_per_mol = m_g_per_mol / 1000.0
                total_M_kg_per_mol += m_kg_per_mol

                v_ms = sqrt(vx_au*vx_au + vy_au*vy_au + vz_au*vz_au) * AU_VELOCITY_TO_M_PER_S
                E_J_per_mol = 0.5 * m_kg_per_mol * v_ms * v_ms
                E_kJ_per_mol = E_J_per_mol / 1000.0
                E_kcal_per_mol = E_kJ_per_mol * 0.239005736

                lines.append(
                    f"KE(atom {idx1:>3d} {sym:>2s}) ...............: "
                    f"{E_kJ_per_mol:.6E} kJ/mol | {E_kcal_per_mol:.6E} kcal/mol"
                )

            # Optional COM kinetic energy (mass-weighted)
            com_enabled = getattr(self, "switch_com", None) and self.switch_com.value
            if com_enabled and len(self.attacker_mol_idx) >= 1:
                # v_COM = sum(m_i v_i) / sum(m_i), with m_i in kg/mol and v_i in m/s
                sum_m_vx = sum_m_vy = sum_m_vz = 0.0
                for idx1 in self.attacker_mol_idx:
                    i = idx1 - 1
                    sym, _, _, _, vx_au, vy_au, vz_au, *_ = self.spv[i]
                    m_kg_per_mol = AtomicData.atomic_masses[sym] / 1000.0
                    sum_m_vx += m_kg_per_mol * (vx_au * AU_VELOCITY_TO_M_PER_S)
                    sum_m_vy += m_kg_per_mol * (vy_au * AU_VELOCITY_TO_M_PER_S)
                    sum_m_vz += m_kg_per_mol * (vz_au * AU_VELOCITY_TO_M_PER_S)

                if total_M_kg_per_mol > 0.0:
                    vx_com = sum_m_vx / total_M_kg_per_mol
                    vy_com = sum_m_vy / total_M_kg_per_mol
                    vz_com = sum_m_vz / total_M_kg_per_mol
                    v_com = sqrt(vx_com*vx_com + vy_com*vy_com + vz_com*vz_com)

                    E_com_J_per_mol = 0.5 * total_M_kg_per_mol * v_com * v_com
                    E_com_kJ_per_mol = E_com_J_per_mol / 1000.0
                    E_com_kcal_per_mol = E_com_kJ_per_mol * 0.239005736

                    lines.append(
                        f"COM KE (attackers) ...............: "
                        f"{E_com_kJ_per_mol:.6E} kJ/mol | {E_com_kcal_per_mol:.6E} kcal/mol"
                    )

            # Write results
            with open(f"{self.output_dir}/summary.txt", "a", encoding="utf-8") as f:
                f.write(f"Molar mass of attackers ..........: {total_M_kg_per_mol:.6E} kg/mol\n")
                for L in lines:
                    f.write(L + "\n")
                f.write("--------------------------------------------\n\n")

        except Exception as e:
            self.multi_line_text.value += f"\nError calculating kinetic energy: {e}\n"

    def setInitialVelocity(self):
        """
        Set attacker COM velocity to point towards the target.
        CPMD preserves attacker internal velocities relative to the original COM.
        gqteaMD uses a pure translational attacker velocity toward the target.
        """
        self._ensure_output_dir()
        try:
            with open(f"{self.output_dir}/summary.txt", "w", encoding="utf-8") as f:
                now = datetime.datetime.now()
                f.write(f"File created at {now:%H:%M} on {now:%d}/{now:%m}/{now:%Y}\n\n")
                f.write(f"Total number of atoms .............: {self.num_atoms}\n\n")
                f.write(f"MD engine option active ...........: {self.md_engine}\n")
                if self.md_engine == "gqteaMD":
                    f.write("Velocity unit in this summary .....: Angstrom/fs\n")
                    f.write("gqteaMD GEOMETRY velocities are written in Angstrom/fs, not a.u.\n\n")
                else:
                    f.write("Velocity unit in this summary .....: a.u.\n\n")
                total_mass = 0.0
                com_x = com_y = com_z = 0.0
                com_vx = com_vy = com_vz = 0.0

                for idx1 in self.attacker_mol_idx:
                    i = idx1 - 1
                    sym, x, y, z, vx, vy, vz, *_ = self.spv[i]
                    try:
                        mass = AtomicData.atomic_masses[sym]
                    except KeyError:
                        raise ValueError(f"Unknown atomic mass for symbol '{sym}'")
                    total_mass += mass
                    com_x += mass * x
                    com_y += mass * y
                    com_z += mass * z
                    com_vx += mass * vx
                    com_vy += mass * vy
                    com_vz += mass * vz

                if total_mass <= 0.0:
                    raise ValueError("Attacker molecule has zero total mass.")

                com_x /= total_mass
                com_y /= total_mass
                com_z /= total_mass
                com_vx /= total_mass
                com_vy /= total_mass
                com_vz /= total_mass

                dx = self.target_pos[0] - com_x
                dy = self.target_pos[1] - com_y
                dz = self.target_pos[2] - com_z
                dist = sqrt(dx*dx + dy*dy + dz*dz)
                if dist == 0.0:
                    raise ValueError("Target position coincides with attacker center of mass.")

                ux, uy, uz = dx / dist, dy / dist, dz / dist
                collision_vx = ux * self.init_attacker_vel
                collision_vy = uy * self.init_attacker_vel
                collision_vz = uz * self.init_attacker_vel

                self.collision_velocity = [collision_vx, collision_vy, collision_vz]

                f.write(
                    f"Attacker COM position (Å) ..........: {com_x:.6E} {com_y:.6E} {com_z:.6E}\n"
                )
                if self.md_engine == "gqteaMD":
                    collision_display = [
                        collision_vx * AU_VELOCITY_TO_ANGSTROM_PER_FS,
                        collision_vy * AU_VELOCITY_TO_ANGSTROM_PER_FS,
                        collision_vz * AU_VELOCITY_TO_ANGSTROM_PER_FS,
                    ]
                    f.write(
                        "Collision velocity vector (Å/fs) ..: "
                        f"{collision_display[0]:.6E} {collision_display[1]:.6E} {collision_display[2]:.6E}\n\n"
                    )
                else:
                    f.write(
                        f"Collision velocity vector (a.u.) ...: {collision_vx:.6E} {collision_vy:.6E} {collision_vz:.6E}\n\n"
                    )

                for idx1 in self.attacker_mol_idx:
                    i = idx1 - 1
                    if self.md_engine == "gqteaMD":
                        self.spv[i][4] = collision_vx
                        self.spv[i][5] = collision_vy
                        self.spv[i][6] = collision_vz
                    else:
                        internal_vx = self.spv[i][4] - com_vx
                        internal_vy = self.spv[i][5] - com_vy
                        internal_vz = self.spv[i][6] - com_vz
                        self.spv[i][4] = internal_vx + collision_vx
                        self.spv[i][5] = internal_vy + collision_vy
                        self.spv[i][6] = internal_vz + collision_vz

                    vel_au = sqrt(
                        self.spv[i][4] ** 2 + self.spv[i][5] ** 2 + self.spv[i][6] ** 2
                    )
                    vel_ms = vel_au * AU_VELOCITY_TO_M_PER_S
                    if self.md_engine == "gqteaMD":
                        vel_afs = vel_au * AU_VELOCITY_TO_ANGSTROM_PER_FS
                        f.write(
                            f"Velocity of atom {idx1} ...............: {vel_afs:.6E} Å/fs or {vel_ms:.6E} m/s\n"
                        )
                    else:
                        f.write(
                            f"Velocity of atom {idx1} ...............: {vel_au:.6E} a.u. or {vel_ms:.6E} m/s\n"
                        )
        except Exception as e:
            self.multi_line_text.value += f"\nError setting initial velocities: {e}\n"

    def set_zero_velocity(self):
        """
        If switch is on, zero velocities of non-attacker atoms using the freshly
        written 'with-vibration-GEOMETRY' (bohr + a.u.), and print that file to UI.
        """
        if not getattr(self, "switch", None) or not self.switch.value:
            return

        self._ensure_output_dir()
        if self.md_engine == "gqteaMD":
            self.write_gqteamd_without_vibration_geometry()
            return

        src = pathlib.Path(self.output_dir) / "with-vibration-GEOMETRY"
        dst = pathlib.Path(self.output_dir) / "without-vibration-GEOMETRY"

        try:
            with open(src, "r", encoding="utf-8") as f:
                lines = f.readlines()

            if len(lines) < self.num_atoms:
                raise ValueError(
                    "with-vibration-GEOMETRY does not contain enough lines."
                )

            self.multi_line_text.value += (
                "\n>>> The following geometry sets velocities of non-attackers to zero <<<\n\n"
            )
            collision_vx, collision_vy, collision_vz = getattr(
                self, "collision_velocity", [0.0, 0.0, 0.0]
            )
            with open(dst, "w", encoding="utf-8") as f2:
                for atom_idx in range(self.num_atoms):
                    parts = lines[atom_idx].split()
                    if len(parts) < 6:
                        raise ValueError(f"Incomplete data at atom line {atom_idx+1}.")

                    x, y, z = map(float, parts[0:3])
                    if (atom_idx + 1) in self.attacker_mol_idx:
                        vx, vy, vz = collision_vx, collision_vy, collision_vz
                    else:
                        vx, vy, vz = 0.0, 0.0, 0.0

                    f2.write(
                        "{:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f}\n"
                        .format(x, y, z, vx, vy, vz)
                    )
                    self.multi_line_text.value += (
                        "{:<5d}{:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f}\n"
                        .format(atom_idx + 1, x, y, z, vx, vy, vz)
                    )
        except Exception as e:
            self.multi_line_text.value += f"\nError setting zero velocities: {e}\n"

    def write_gqteamd_without_vibration_geometry(self):
        dst = pathlib.Path(self.output_dir) / "gqteaMD-without-vibration-GEOMETRY"
        try:
            collision_vx, collision_vy, collision_vz = getattr(
                self, "collision_velocity", [0.0, 0.0, 0.0]
            )
            self.multi_line_text.value += (
                "\n>>> The following gqteaMD geometry sets non-attackers to zero and attackers to the collision COM velocity <<<\n\n"
            )
            with open(dst, "w", encoding="utf-8") as f_geom:
                f_geom.write(f"{self.num_atoms}\n")
                f_geom.write(
                    "Generated by gqteaWinToga collision setup; cartesian Angstrom positions, Angstrom/fs velocities\n"
                )
                for atom_idx, (sym, x, y, z, vx, vy, vz, fx, fy, fz) in enumerate(self.spv, start=1):
                    if atom_idx in self.attacker_mol_idx:
                        out_vx, out_vy, out_vz = collision_vx, collision_vy, collision_vz
                    else:
                        out_vx, out_vy, out_vz = 0.0, 0.0, 0.0
                    out_vx *= AU_VELOCITY_TO_ANGSTROM_PER_FS
                    out_vy *= AU_VELOCITY_TO_ANGSTROM_PER_FS
                    out_vz *= AU_VELOCITY_TO_ANGSTROM_PER_FS

                    line = (
                        "{:<2} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f} {:>20.12f}\n"
                        .format(sym, x, y, z, out_vx, out_vy, out_vz, fx, fy, fz)
                    )
                    f_geom.write(line)
                    self.multi_line_text.value += line
        except Exception as e:
            self.multi_line_text.value += f"\nError setting gqteaMD zero velocities: {e}\n"


class CollisionUI(CollisionSetUp):
    def __init__(self, *args):
        self.layout_main_window(*args)

    def layout_main_window(self, widget):
        """Create & arrange the GUI."""
        self.main_window = toga.Window(
            title="CPMD Input Files to Simulate Collision Molecular Dynamics",
            size=(720, 620),
        )

        heading_style = Pack(font_size=18, font_weight="bold", margin=(0, 0, 10, 0))
        label_style = Pack(margin=(5, 5), text_align=LEFT, width=240)
        input_style = Pack(flex=1, margin=(5, 5))
        button_style = Pack(margin=5, width=130)
        box_style = Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0))

        main_box = toga.Box(style=Pack(direction=COLUMN, margin=20))

        # Title
        box_title = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 10, 0)))
        title_label = toga.Label("Collision Input Files", style=heading_style)
        box_title.add(title_label)
        main_box.add(box_title)

        # Inputs
        engine_box = toga.Box(style=box_style)
        engine_label = toga.Label("MD engine:", style=label_style)
        self.md_engine_selection = toga.Selection(
            items=MD_ENGINE_CHOICES,
            value="cpmd",
            on_change=self.update_initial_velocity_label,
            style=Pack(width=160, margin=(5, 5)),
        )
        engine_box.add(engine_label)
        engine_box.add(self.md_engine_selection)
        main_box.add(engine_box)

        input_fields = [
            ("Attacker atom indices (1-based):", "e.g. 12 13 14", "textInput_atomic_labels"),
            ("Initial velocity (a.u.):", "e.g. 0.010", "textInput_initial_velocity"),
            ("Target xyz coordinates (Å):", "e.g. 10.0 12.5 8.0", "textInput_coords_site"),
        ]
        for label_text, placeholder, attr_name in input_fields:
            box = toga.Box(style=box_style)
            label = toga.Label(label_text, style=label_style)
            if attr_name == "textInput_initial_velocity":
                self.initial_velocity_label = label
            text_input = toga.TextInput(placeholder=placeholder, style=input_style)
            setattr(self, attr_name, text_input)
            box.add(label)
            box.add(text_input)
            main_box.add(box)

        # Switches row
        switches = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 5)))
        self.switch = toga.Switch("Zero velocities of non-attacker atoms")
        self.switch_com = toga.Switch("Compute center-of-mass KE")
        switches.add(self.switch)
        switches.add(self.switch_com)

        # File selection
        file_box = toga.Box(style=box_style)
        file_label = toga.Label("Select GEOMETRY file:", style=label_style)
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select GEOMETRY or GEOMETRY.xyz",
            style=Pack(flex=1, margin=(5, 5), color="blue"),
            readonly=True,
        )
        browse_button = toga.Button("Browse", on_press=self.readGeometry, style=button_style)
        file_box.add(file_label)
        file_box.add(self.textInput_file)
        file_box.add(browse_button)

        main_box.add(file_box)
        main_box.add(switches)

        # Output area
        self.multi_line_text = toga.MultilineTextInput(
            style=Pack(flex=1, margin=(10, 0), font_size=12)
        )
        self.multi_line_text.value = HelpGqteaWin.help_collision
        main_box.add(self.multi_line_text)

        # Buttons
        button_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin_top=10))
        self.btn_execute = toga.Button("Input Builder", style=button_style, on_press=self.execute)
        self.btn_close = toga.Button("Close", style=button_style, on_press=self.closeTopLevel)
        button_box.add(self.btn_execute)
        button_box.add(self.btn_close)
        main_box.add(button_box)

        self.main_window.content = main_box
        self.main_window.show()

    def update_initial_velocity_label(self, widget=None):
        if not hasattr(self, "initial_velocity_label"):
            return
        if self.get_selected_md_engine() == "gqteaMD":
            self.initial_velocity_label.text = "Initial velocity (Å/fs):"
        else:
            self.initial_velocity_label.text = "Initial velocity (a.u.):"

    def closeTopLevel(self, widget):
        self.main_window.close()
