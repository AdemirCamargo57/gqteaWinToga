"""Molecular orbitals of a Gaussian .fchk for the 3D Molecular Viewer's Orbitals tab.

Numerical only: no Toga (a subprocess test keeps it out). The tab
(``orbitalTab.OrbitalTabController``) lists the orbitals of
``load_orbitals`` and draws the grids of ``orbital_volumes``.

* **Catalogue.** Restricted (one MO set, occupancy 2/0), unrestricted (alpha
  and beta sets, 1/0 each) and restricted-open (one set, 2/1/0 from
  n_alpha/n_beta) SCF wavefunctions. Labels count from each spin's HOMO.
* **Orbitals on a grid.** ``gaussianWavefunction.orbital_values`` (one AO pass
  for any number of MOs, the compiled AO kernel when numba is installed) on
  one cubic grid whose faces are below isovalue / 10 for every requested MO
  (``molecularSurface.grid_box``, as the NTO tab). Signs are the fchk
  coefficients as written, so the phases are GaussView's and cubegen's.
* **Files.** Orbital cubes, the MO energy table and the energy-level diagram
  (a plotViewer manifest of type "energy_levels", drawn by
  ``displayPlots.draw_energy_levels``).

Lengths are bohr; energies Hartree unless a name says eV.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

import cubeFile as cf
import gaussianWavefunction as gw
import molecularSurface as ms
from help import AtomicData
from indexRanges import parse_solute_index_ranges
from ntoAnalysis import check_grid_size, resolve_output_dir  # noqa: F401  (resolve_output_dir: re-exported)

RESTRICTED, UNRESTRICTED, RESTRICTED_OPEN = "restricted", "unrestricted", "restricted-open"
ALPHA, BETA = "alpha", "beta"
HARTREE_TO_EV = 27.211386245988   # CODATA 2018, as the Fukui and Excitation tools
BOHR_RADIUS_ANGSTROM = cf.BOHR_RADIUS_ANGSTROM

ORBITAL_ISOVALUE = 0.02           # GaussView's MO default (a.u.)
DENSITY_ISOVALUE = 0.0004         # |phi|^2; = 0.02^2, the same surface as the default phases
GRID_PRESETS = {"Coarse": 3, "Medium": 6, "Fine": 12}   # points per bohr (cubegen's -2/-3/-4)
DEFAULT_GRID = "Medium"
SLAB_POINTS = 2_000_000           # grid points evaluated at a time
VALUES_PER_BATCH = 50_000_000     # grid values held at once when writing several cubes (~400 MB)
MAX_EXPORT_ORBITALS = 50          # cubes written by one range export
LEVEL_WINDOW = 10                 # the diagram opens on HOMO-10 ... LUMO+10,
LEVEL_SPAN_EV = 20.0              # ... but never more than 20 eV below the HOMO / above the LUMO (core levels)
DEGENERACY_EV = 0.003             # levels closer than this (~1e-4 Ha) are drawn side by side
CUBE_COMMENT = "Made by the gQTEA 3D Molecular Viewer (Orbitals tab)"


def preset_spacing(name: str) -> float:
    """Grid spacing (bohr) of a Coarse/Medium/Fine preset."""
    if name not in GRID_PRESETS:
        raise ValueError(f"Unknown grid quality '{name}': choose {', '.join(GRID_PRESETS)}.")
    return 1.0 / GRID_PRESETS[name]


# ----------------------------------------------------------------------
# Catalogue
# ----------------------------------------------------------------------
def orbital_label(number: int, homo: Optional[int]) -> str:
    """HOMO, HOMO-1, ..., LUMO, LUMO+1, ... counted from ``homo`` (None: no occupied MO)."""
    homo = homo or 0
    if number <= homo:
        return "HOMO" if number == homo else f"HOMO-{homo - number}"
    return "LUMO" if number == homo + 1 else f"LUMO+{number - homo - 1}"


@dataclass(frozen=True)
class OrbitalInfo:
    number: int          # 1-based, as Gaussian and GaussView number them
    spin: str            # ALPHA or BETA (restricted files: ALPHA)
    energy: float        # Hartree
    occupancy: float     # 2/1/0 restricted(-open), 1/0 unrestricted
    label: str           # relative to the spin's HOMO

    @property
    def energy_ev(self) -> float:
        return self.energy * HARTREE_TO_EV


@dataclass
class OrbitalCatalogue:
    wavefunction: gw.FchkWavefunction = field(repr=False)
    kind: str
    orbitals: Dict[str, List[OrbitalInfo]]          # ALPHA always; BETA for unrestricted files
    homo: Dict[str, Optional[int]]

    @property
    def path(self) -> str:
        return self.wavefunction.path

    @property
    def stem(self) -> str:
        return Path(self.wavefunction.path).stem

    @property
    def n_mo(self) -> int:
        return len(self.orbitals[ALPHA])

    @property
    def spins(self) -> Tuple[str, ...]:
        return tuple(self.orbitals)

    @property
    def unrestricted(self) -> bool:
        return self.kind == UNRESTRICTED

    def orbital(self, spin: str, number: int) -> OrbitalInfo:
        if not 1 <= number <= self.n_mo:
            raise ValueError(f"MO {number} does not exist: {Path(self.path).name} has {self.n_mo} MOs.")
        return self.orbitals[spin][number - 1]

    def lumo(self, spin: str) -> Optional[int]:
        homo = self.homo[spin] or 0
        return homo + 1 if homo < self.n_mo else None

    def gap(self, spin: str) -> Optional[float]:
        """LUMO - HOMO of one spin (Hartree), None when either is missing."""
        homo, lumo = self.homo[spin], self.lumo(spin)
        if not homo or lumo is None:
            return None
        return self.orbital(spin, lumo).energy - self.orbital(spin, homo).energy

    def spin_text(self, spin: str) -> str:
        """'' for one-set files, ' alpha' / ' beta' for unrestricted ones."""
        return f" {spin}" if self.unrestricted else ""

    def summary(self) -> str:
        wfn = self.wavefunction
        parts = [f"{wfn.method}/{wfn.basis_set}, {wfn.n_electrons} electrons, {self.n_mo} MOs ({self.kind})."]
        for spin in self.spins:
            homo, lumo, gap = self.homo[spin], self.lumo(spin), self.gap(spin)
            name = f"{spin.capitalize()} " if self.unrestricted else ""
            text = []
            if homo:
                text.append(f"{name}HOMO {homo}: {self.orbital(spin, homo).energy_ev:.3f} eV")
            if lumo is not None:
                text.append(f"LUMO {lumo}: {self.orbital(spin, lumo).energy_ev:.3f} eV")
            if gap is not None:
                text.append(f"gap {gap * HARTREE_TO_EV:.3f} eV")
            parts.append(", ".join(text) + ".")
        return " ".join(parts)


def _occupancies(kind: str, spin: str, n_mo: int, n_alpha: int, n_beta: int) -> np.ndarray:
    numbers = np.arange(1, n_mo + 1)
    if kind == RESTRICTED:
        return np.where(numbers <= n_alpha, 2.0, 0.0)
    if kind == RESTRICTED_OPEN:
        return np.where(numbers <= n_beta, 2.0, np.where(numbers <= n_alpha, 1.0, 0.0))
    return np.where(numbers <= (n_alpha if spin == ALPHA else n_beta), 1.0, 0.0)


def orbital_catalogue(wfn: gw.FchkWavefunction) -> OrbitalCatalogue:
    """The MOs of an SCF wavefunction. Every problem is a ValueError naming it."""
    name = Path(wfn.path).name
    if wfn.alpha_mo is None:
        raise ValueError(f"{name} has no molecular orbitals (no 'Alpha MO coefficients' section).")
    if wfn.alpha_orbital_energies is None:
        raise ValueError(f"{name} has no 'Alpha Orbital Energies' section.")
    n_mo = wfn.alpha_mo.shape[1]
    if wfn.is_unrestricted:
        kind = UNRESTRICTED
        if wfn.beta_orbital_energies is None:
            raise ValueError(f"{name} has beta MO coefficients but no 'Beta Orbital Energies' section.")
    else:
        kind = RESTRICTED if wfn.n_alpha == wfn.n_beta else RESTRICTED_OPEN
    energies = {ALPHA: wfn.alpha_orbital_energies}
    if kind == UNRESTRICTED:
        energies[BETA] = wfn.beta_orbital_energies
    for spin, values in energies.items():
        if len(values) != n_mo:
            raise ValueError(f"{name}: {len(values)} {spin} orbital energies for {n_mo} MOs.")
    if not 0 <= wfn.n_beta <= wfn.n_alpha <= n_mo:
        raise ValueError(f"{name}: {wfn.n_alpha} alpha and {wfn.n_beta} beta electrons do not fit "
                         f"{n_mo} MOs.")

    orbitals, homos = {}, {}
    for spin, values in energies.items():
        occupancy = _occupancies(kind, spin, n_mo, wfn.n_alpha, wfn.n_beta)
        occupied = np.nonzero(occupancy > 0)[0]
        homo = int(occupied[-1]) + 1 if len(occupied) else None
        homos[spin] = homo
        orbitals[spin] = [OrbitalInfo(k + 1, spin, float(values[k]), float(occupancy[k]),
                                      orbital_label(k + 1, homo)) for k in range(n_mo)]
    return OrbitalCatalogue(wfn, kind, orbitals, homos)


def load_orbitals(path) -> OrbitalCatalogue:
    """Read a .fchk/.fch and list its MOs."""
    return orbital_catalogue(gw.read_fchk(path))


# ----------------------------------------------------------------------
# Orbital selections typed by the user
# ----------------------------------------------------------------------
_TOKEN = re.compile(r"^(?:(HOMO|LUMO)([+-]\d+)?|(\d+))$", re.IGNORECASE)
_NUMBER_RANGE = re.compile(r"^\d+(?:-\d+)?$")


def _resolve(token: str, catalogue: OrbitalCatalogue, spin: str) -> int:
    match = _TOKEN.match(token)
    if match is None:
        raise ValueError(f"'{token}' is not an orbital: use a number, HOMO, LUMO, HOMO-n or LUMO+n.")
    if match.group(3) is not None:
        return int(match.group(3))
    homo = catalogue.homo[spin] or 0
    base = homo if match.group(1).upper() == "HOMO" else homo + 1
    return base + int(match.group(2) or 0)


def parse_orbital_spec(text: str, catalogue: OrbitalCatalogue, spin: str = ALPHA) -> List[int]:
    """MO numbers from text such as ``20-25, 30`` or ``HOMO-2:LUMO+2`` (sorted, unique, 1-based).

    Number ranges use '-' (the project's range parser); ranges with HOMO/LUMO
    use ':' because '-' is their offset. Separators: commas or spaces.
    """
    parts = [p for p in re.split(r"[,\s]+", (text or "").strip()) if p]
    if not parts:
        raise ValueError("Enter the orbitals, e.g. 'HOMO-2:LUMO+2' or '20-25, 30'.")
    numbers = set()
    for part in parts:
        if ":" in part:
            ends = part.split(":")
            if len(ends) != 2 or not all(ends):
                raise ValueError(f"'{part}': a range is two orbitals joined by ':', e.g. HOMO-2:LUMO+2.")
            a, b = (_resolve(e, catalogue, spin) for e in ends)
            numbers.update(range(min(a, b), max(a, b) + 1))
        elif _NUMBER_RANGE.match(part):
            try:
                numbers.update(parse_solute_index_ranges(part))
            except ValueError:  # its messages speak of atom indices
                raise ValueError(f"'{part}' is not a valid MO range: write it from low to high with "
                                 "numbers from 1, e.g. 20-25.") from None
        else:
            numbers.add(_resolve(part, catalogue, spin))
    for number in sorted(numbers):
        if not 1 <= number <= catalogue.n_mo:
            raise ValueError(f"MO {number} does not exist: {Path(catalogue.path).name} has "
                             f"{catalogue.n_mo} MOs.")
    return sorted(numbers)


# ----------------------------------------------------------------------
# Grids
# ----------------------------------------------------------------------
def _vdw_radii_bohr(wfn) -> np.ndarray:
    return np.array([AtomicData.vdw_radii.get(s, AtomicData.DEFAULT_VDW_RADIUS)
                     for s in wfn.atom_symbols()]) / BOHR_RADIUS_ANGSTROM


def _check_grid_request(spacing: float, isovalue: float) -> None:
    if not (spacing > 0 and np.isfinite(spacing)):
        raise ValueError(f"The grid spacing must be a positive number of bohr, not {spacing}.")
    if not (isovalue > 0 and np.isfinite(isovalue)):
        raise ValueError(f"The isovalue must be a positive number, not {isovalue}.")


def orbital_grid(catalogue: OrbitalCatalogue, spin: str, numbers: Sequence[int], spacing: float,
                 isovalue: float, use_numba: Optional[bool] = None) -> Tuple[np.ndarray, Tuple[int, int, int]]:
    """Origin and shape of one grid for all ``numbers``: atoms +/- 1.7 r_vdW, each face grown until
    every MO is below isovalue / 10 on it. The 25-million-point limit is checked here."""
    _check_grid_request(spacing, isovalue)
    wfn = catalogue.wavefunction
    for number in numbers:
        catalogue.orbital(spin, number)  # range check
    coefficients = gw.mo_coefficients(wfn, spin)[:, [n - 1 for n in numbers]]

    def envelope(points):
        return np.abs(gw.orbital_values(wfn, coefficients, points, use_numba=use_numba)).max(axis=1)

    origin, shape = ms.grid_box(wfn.coordinates, _vdw_radii_bohr(wfn), envelope, isovalue, spacing)
    check_grid_size(shape, spacing)
    return np.asarray(origin, dtype=float), tuple(int(n) for n in shape)


def evaluate_orbitals(catalogue: OrbitalCatalogue, spin: str, numbers: Sequence[int], origin: np.ndarray,
                      spacing: float, shape: Tuple[int, int, int],
                      progress: Optional[Callable[[float], None]] = None,
                      use_numba: Optional[bool] = None) -> np.ndarray:
    """Every value exact: (len(numbers), nx, ny, nz), evaluated in slabs of x planes."""
    wfn = catalogue.wavefunction
    coefficients = gw.mo_coefficients(wfn, spin)[:, [n - 1 for n in numbers]]
    nx, ny, nz = shape
    values = np.empty((len(numbers), nx, ny, nz))
    planes = max(1, SLAB_POINTS // (ny * nz))
    for x0 in range(0, nx, planes):
        x1 = min(nx, x0 + planes)
        points = ms.grid_points(origin + np.array([x0 * spacing, 0.0, 0.0]), spacing, (x1 - x0, ny, nz))
        slab = gw.orbital_values(wfn, coefficients, points, use_numba=use_numba)
        values[:, x0:x1] = slab.T.reshape(len(numbers), x1 - x0, ny, nz)
        if progress is not None:
            try:
                progress(x1 / nx)
            except Exception:
                pass  # reporting must never abort a calculation
    return values


def orbital_title(catalogue: OrbitalCatalogue, info: OrbitalInfo, density: bool = False) -> str:
    what = "|MO|^2" if density else "MO"
    return (f"{what} {info.number} ({info.label}){catalogue.spin_text(info.spin)}, occupancy {info.occupancy:g}, "
            f"energy {info.energy:.6f} Ha ({info.energy_ev:.4f} eV)")


def _volume(catalogue: OrbitalCatalogue, info: OrbitalInfo, origin, spacing, values,
            density: bool = False) -> cf.CubeVolume:
    wfn = catalogue.wavefunction
    return cf.CubeVolume(orbital_title(catalogue, info, density), CUBE_COMMENT,
                         np.asarray(wfn.atomic_numbers), np.asarray(wfn.nuclear_charges, dtype=float),
                         wfn.coordinates.copy(), np.asarray(origin, dtype=float).copy(), spacing * np.eye(3),
                         np.ascontiguousarray(values), orbital=None if density else info.number)


def orbital_volumes(catalogue: OrbitalCatalogue, spin: str, numbers: Sequence[int],
                    spacing: float = 1.0 / GRID_PRESETS[DEFAULT_GRID], isovalue: float = ORBITAL_ISOVALUE,
                    progress: Optional[Callable[[float], None]] = None,
                    use_numba: Optional[bool] = None) -> List[cf.CubeVolume]:
    """The orbitals ``numbers`` of ``spin`` (amplitude, a.u.) on one shared grid."""
    numbers = list(numbers)
    if not numbers:
        raise ValueError("Choose at least one orbital.")
    origin, shape = orbital_grid(catalogue, spin, numbers, spacing, isovalue, use_numba)
    values = evaluate_orbitals(catalogue, spin, numbers, origin, spacing, shape, progress, use_numba)
    return [_volume(catalogue, catalogue.orbital(spin, n), origin, spacing, values[k])
            for k, n in enumerate(numbers)]


def density_volume(catalogue: OrbitalCatalogue, volume: cf.CubeVolume, spin: str, number: int) -> cf.CubeVolume:
    """|phi|^2 of an amplitude volume, on the same grid (no new evaluation)."""
    info = catalogue.orbital(spin, number)
    return _volume(catalogue, info, volume.origin, float(volume.axes[0, 0]), np.square(volume.values),
                   density=True)


# ----------------------------------------------------------------------
# Files
# ----------------------------------------------------------------------
def cube_filename(stem: str, number: int, spin: Optional[str] = None, density: bool = False) -> str:
    """<stem>_MO<n>[_alpha|_beta][_density].cube (spin only for unrestricted files)."""
    return f"{stem}_MO{number}{'_' + spin if spin else ''}{'_density' if density else ''}.cube"


def write_cubes(catalogue: OrbitalCatalogue, folder: str, spin: str, numbers: Sequence[int], spacing: float,
                isovalue: float, density: bool = False, progress: Optional[Callable[[float], None]] = None,
                use_numba: Optional[bool] = None) -> List[str]:
    """One cube per MO, all on one grid (so the Cube tab can add or subtract them). Returns the paths.

    The orbitals are evaluated in batches of at most VALUES_PER_BATCH grid
    values, so a long range does not hold every grid in memory.
    """
    numbers = list(numbers)
    if not numbers:
        raise ValueError("Choose at least one orbital.")
    if len(numbers) > MAX_EXPORT_ORBITALS:
        raise ValueError(f"{len(numbers)} orbitals requested; export at most {MAX_EXPORT_ORBITALS} cubes at a time.")
    origin, shape = orbital_grid(catalogue, spin, numbers, spacing, isovalue, use_numba)
    batch = max(1, VALUES_PER_BATCH // int(np.prod(shape)))
    spin_part = spin if catalogue.unrestricted else None
    written: List[str] = []
    try:
        Path(folder).mkdir(parents=True, exist_ok=True)
        for start in range(0, len(numbers), batch):
            chunk = numbers[start:start + batch]

            def step(fraction, done=start, size=len(chunk)):
                if progress is not None:
                    progress((done + fraction * size) / len(numbers))

            values = evaluate_orbitals(catalogue, spin, chunk, origin, spacing, shape, step, use_numba)
            for k, number in enumerate(chunk):
                grid = np.square(values[k]) if density else values[k]
                volume = _volume(catalogue, catalogue.orbital(spin, number), origin, spacing, grid, density)
                path = Path(folder) / cube_filename(catalogue.stem, number, spin_part, density)
                written.append(cf.write_cube_volume(path, volume, use_numba))
    except OSError as exc:
        raise ValueError(f"Cannot write the cube files in {folder}: {exc}") from exc
    return written


def write_volume(volume: cf.CubeVolume, path, use_numba: Optional[bool] = None) -> str:
    """Write one already computed volume; OSError becomes a ValueError naming the file."""
    try:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        return cf.write_cube_volume(path, volume, use_numba)
    except OSError as exc:
        raise ValueError(f"Cannot write {path}: {exc}") from exc


def energy_table_filename(stem: str) -> str:
    return f"{stem}_MO_energies.txt"


def energy_table_text(catalogue: OrbitalCatalogue) -> str:
    wfn = catalogue.wavefunction
    lines = [
        "# gQTEA Molecular Orbital Energies",
        f"# file  {Path(catalogue.path).name}",
        f"# method  {wfn.method}",
        f"# basis  {wfn.basis_set}",
        f"# wavefunction  {catalogue.kind}",
        f"# electrons  {wfn.n_electrons} (alpha {wfn.n_alpha}, beta {wfn.n_beta})",
        f"# charge  {wfn.charge}",
        f"# multiplicity  {wfn.multiplicity}",
        f"# basis_functions  {wfn.n_basis}",
        f"# molecular_orbitals  {catalogue.n_mo}",
        f"# hartree_to_ev  {HARTREE_TO_EV}",
    ]
    for spin in catalogue.spins:
        suffix = f"_{spin}" if catalogue.unrestricted else ""
        for name, number in (("homo", catalogue.homo[spin]), ("lumo", catalogue.lumo(spin))):
            if number:
                info = catalogue.orbital(spin, number)
                lines.append(f"# {name}{suffix}  {number}  {info.energy:.6f} Ha  {info.energy_ev:.4f} eV")
        gap = catalogue.gap(spin)
        if gap is not None:
            lines.append(f"# gap{suffix}  {gap:.6f} Ha  {gap * HARTREE_TO_EV:.4f} eV")
    lines.append("# mo spin occupancy label energy_hartree energy_ev")
    for spin in catalogue.spins:
        for info in catalogue.orbitals[spin]:
            lines.append(f"{info.number:5d} {spin:5s} {info.occupancy:g} {info.label:9s} "
                         f"{info.energy:14.6f} {info.energy_ev:12.4f}")
    return "\n".join(lines) + "\n"


def write_energy_table(catalogue: OrbitalCatalogue, folder: str) -> str:
    path = Path(folder) / energy_table_filename(catalogue.stem)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(energy_table_text(catalogue), encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"Cannot write {path}: {exc}") from exc
    return str(path)


def levels_filename(stem: str) -> str:
    return f"{stem}_MO_levels.json"


def energy_levels_entry(catalogue: OrbitalCatalogue) -> dict:
    """The "energy_levels" plotViewer figure: one column per spin set, opened on HOMO-10 ... LUMO+10
    (at most LEVEL_SPAN_EV beyond the frontier, so a small molecule's core levels do not squash it)."""
    columns, window = [], []
    for spin in catalogue.spins:
        homo, lumo = catalogue.homo[spin], catalogue.lumo(spin)
        levels = [{"number": o.number, "energy_ev": o.energy_ev, "occupancy": o.occupancy, "label": o.label}
                  for o in catalogue.orbitals[spin]]
        columns.append({"label": spin.capitalize() if catalogue.unrestricted else "",
                        "levels": levels, "homo": homo, "lumo": lumo})
        first = max(1, (homo or 1) - LEVEL_WINDOW)
        last = min(catalogue.n_mo, (lumo or catalogue.n_mo) + LEVEL_WINDOW)
        bottom = levels[(homo or 1) - 1]["energy_ev"] - LEVEL_SPAN_EV
        top = levels[(lumo or catalogue.n_mo) - 1]["energy_ev"] + LEVEL_SPAN_EV
        window += [max(levels[first - 1]["energy_ev"], bottom), min(levels[last - 1]["energy_ev"], top)]
    lo, hi = min(window), max(window)
    pad = max(1.0, 0.05 * (hi - lo))
    wfn = catalogue.wavefunction
    return {
        "type": "energy_levels",
        "title": f"Molecular orbital energies: {Path(catalogue.path).name}",
        "ylabel": "Energy (eV)",
        "columns": columns,
        "ylim": [lo - pad, hi + pad],
        "source": {"fchk": Path(catalogue.path).name, "method": wfn.method, "basis": wfn.basis_set,
                   "wavefunction": catalogue.kind},
    }


def write_levels(entry: dict, path) -> str:
    """A one-figure manifest that plotViewer and the Plotter's JSON type both open."""
    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps([entry]), encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"Cannot write {path}: {exc}") from exc
    return str(path)
