"""Natural transition orbitals (NTOs) of a CIS / TDA / TD excited state.

Pure numpy, no Toga: the 3D Molecular Viewer's NTO tab (ntoTab) builds on
it and the tests drive it headless. Lengths are bohr.

Multiwfn's definition (manual section 3.21.6; excittrans.f90, subroutine NTO):
the transition matrix of one spin is T[i, a] = w(i -> a) over the occupied
(i) and virtual (a) MOs, from the excitation coefficients only - the
de-excitation ("<-") coefficients of TDHF/TD-DFT are ignored, which is also
what Gaussian's Pop=NTO does. Multiwfn diagonalises T T^T and T^T T; here
T = U S V^T is decomposed once (SVD), which gives the same eigenvalues
(S^2) and pairs each particle with its hole (v_k = T^T u_k / s_k):

    weight_k   = f s_k^2      (f = 2 for a restricted reference, whose
                               coefficients Gaussian normalises to 0.5;
                               1 per spin for an unrestricted one)
    hole_k     = C_occ u_k
    particle_k = C_vir v_k

weight_k x 100 is the percentage of the excitation carried by pair k.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

import cubeFile as cf
import excitationData as ed
import gaussianWavefunction as gw
import molecularSurface as ms
from help import AtomicData

LISTED_WEIGHT = 1e-4           # pairs listed in the table and the report (at least the largest)
NORMALIZATION_TOLERANCE = 0.01  # |deviation| above which the coefficients look incomplete
DOMINANT_WEIGHT = 0.5          # below this the largest pair does not "dominate"
COMPOSITION_THRESHOLD = 0.05   # MO contributions listed for each NTO
SPIN_ORDER = (ed.SPIN_ALPHA, ed.SPIN_BETA)
KINDS = ("hole", "particle")

DEFAULT_SPACING_BOHR = 0.2
ORBITAL_ISOVALUE = 0.02        # GaussView's / the Cube tab's orbital default
MAX_GRID_POINTS = 25_000_000   # as the Fukui tool's cubes (~13 bytes per value in a cube file)
SLAB_POINTS = 2_000_000        # points evaluated at a time
BOHR_RADIUS_ANGSTROM = cf.BOHR_RADIUS_ANGSTROM
CUBE_COMMENT = "Made by the gQTEA 3D Molecular Viewer (NTO tab); orbital amplitude (a.u.)"
FCHK_VALUES_PER_LINE = 5


@dataclass
class NTOSet:
    """The NTOs of one spin (the only set of a restricted reference)."""
    spin: str
    n_occ: int
    weights: np.ndarray       # (min(n_occ, n_vir),), descending
    hole_mo: np.ndarray       # U (n_occ, n_occ): hole k = column k-1, over the occupied MOs
    particle_mo: np.ndarray   # V (n_vir, n_vir): particle k = column k-1, over the virtual MOs
    holes: np.ndarray         # (n_basis, n_occ) in the fchk basis
    particles: np.ndarray     # (n_basis, n_vir)

    def composition(self, kind: str, pair: int, threshold: float = COMPOSITION_THRESHOLD
                    ) -> List[Tuple[int, float]]:
        """(MO number, fraction) of the MOs making up hole/particle ``pair`` (1-based), largest first."""
        if kind == "hole":
            column, offset = self.hole_mo[:, pair - 1], 1
        elif kind == "particle":
            column, offset = self.particle_mo[:, pair - 1], self.n_occ + 1
        else:
            raise ValueError(f"kind must be 'hole' or 'particle', not '{kind}'.")
        fractions = column ** 2
        order = np.argsort(-fractions, kind="stable")
        return [(int(j) + offset, float(fractions[j])) for j in order if fractions[j] >= threshold]


@dataclass(frozen=True)
class NTOPair:
    spin: str
    pair: int      # 1-based, in this spin's weight order
    weight: float


@dataclass
class NTOResults:
    fchk_path: str
    excitation_path: str
    wavefunction: gw.FchkWavefunction = field(repr=False)
    excitations: ed.ExcitationFile = field(repr=False)
    state: ed.ExcitedState = field(repr=False)
    restricted: bool
    factor: float
    sets: Dict[str, NTOSet]
    sum_excitation: float       # sum of w^2 over "->" lines
    sum_deexcitation: float     # sum of w'^2 over "<-" lines
    expected_norm: float        # 0.5 restricted, 1 unrestricted
    notes: List[str]
    elapsed: float = 0.0

    @property
    def deviation(self) -> float:
        """Multiwfn's "Deviation to expected normalization value"."""
        return self.sum_excitation - self.sum_deexcitation - self.expected_norm

    def weight_sum(self, spin: Optional[str] = None) -> float:
        sets = [self.sets[spin]] if spin else self.sets.values()
        return float(sum(s.weights.sum() for s in sets))

    def listed_pairs(self, minimum: float = LISTED_WEIGHT) -> List[NTOPair]:
        """Pairs of every spin with weight >= ``minimum``, largest first (at least the largest pair)."""
        pairs = [NTOPair(spin, k + 1, float(w)) for spin in SPIN_ORDER if spin in self.sets
                 for k, w in enumerate(self.sets[spin].weights)]
        pairs.sort(key=lambda p: (-p.weight, SPIN_ORDER.index(p.spin), p.pair))
        listed = [p for p in pairs if p.weight >= minimum]
        return listed or pairs[:1]

    def nto_set(self, spin: str) -> NTOSet:
        if spin not in self.sets:
            raise ValueError(f"There are no {spin} NTOs (the reference is "
                             f"{'restricted' if self.restricted else 'unrestricted'}).")
        return self.sets[spin]

    def pair_orbitals(self, spin: str, pair: int) -> Tuple[np.ndarray, np.ndarray, float]:
        """(hole coefficients, particle coefficients, weight) of ``pair`` (1-based), fchk basis."""
        s = self.nto_set(spin)
        if not 1 <= int(pair) <= len(s.weights):
            raise ValueError(f"There is no {spin} NTO pair {pair}; the pairs are 1-{len(s.weights)}.")
        k = int(pair) - 1
        return s.holes[:, k], s.particles[:, k], float(s.weights[k])

    @property
    def stem(self) -> str:
        return Path(self.fchk_path).stem


def transition_matrix(state: ed.ExcitedState, spin: str, n_occ: int, n_mo: int) -> np.ndarray:
    """T[i-1, a-n_occ-1] = w(i -> a) over ``spin``'s excitation lines; de-excitations are ignored."""
    t = np.zeros((n_occ, n_mo - n_occ))
    for c in state.configurations:
        if c.spin != spin or c.deexcitation:
            continue
        if not (1 <= c.occupied <= n_occ < c.virtual <= n_mo):
            raise ValueError(f"Line {c.line}: {c.occupied} -> {c.virtual} is not an occupied -> virtual pair "
                             f"for {spin} ({n_occ} occupied of {n_mo} orbitals).")
        t[c.occupied - 1, c.virtual - n_occ - 1] = c.coefficient
    return t


def natural_transition_orbitals(t: np.ndarray, c_mo: np.ndarray, n_occ: int, factor: float,
                                spin: str) -> NTOSet:
    """SVD of ``t``; the joint sign of each pair makes its hole's largest MO coefficient positive."""
    u, s, vt = np.linalg.svd(t, full_matrices=True)
    for k in range(len(s)):
        j = int(np.argmax(np.abs(u[:, k])))
        if u[j, k] < 0:
            u[:, k] *= -1.0
            vt[k] *= -1.0
    v = vt.T
    return NTOSet(spin, n_occ, factor * s ** 2, u, v, c_mo[:, :n_occ] @ u, c_mo[:, n_occ:] @ v)


# ----------------------------------------------------------------------
# Grids and cubes
# ----------------------------------------------------------------------
def resolve_output_dir(fchk_path: str, chosen: Optional[str]) -> str:
    """The chosen folder, or the folder of the .fchk file when blank."""
    chosen = (chosen or "").strip()
    return chosen if chosen else str(Path(fchk_path).resolve().parent)


def check_grid_size(shape: Tuple[int, int, int], spacing: float) -> None:
    points = int(np.prod(shape))
    if points > MAX_GRID_POINTS:
        smallest = np.ceil(100.0 * spacing * (points / MAX_GRID_POINTS) ** (1.0 / 3.0)) / 100.0
        raise ValueError(f"A grid at {spacing:g} bohr would have {points / 1e6:.1f} million points. "
                         f"Use a grid spacing of at least {smallest:g} bohr.")


def _vdw_radii_bohr(wfn) -> np.ndarray:
    return np.array([AtomicData.vdw_radii.get(s, AtomicData.DEFAULT_VDW_RADIUS)
                     for s in wfn.atom_symbols()]) / BOHR_RADIUS_ANGSTROM


def _orbital_title(results: NTOResults, spin: str, pair: int, kind: str, weight: float) -> str:
    spin_text = "" if results.restricted else f" ({spin})"
    return f"NTO {pair} {kind}{spin_text} of excited state {results.state.number}, weight {weight:.6f}"


def pair_volumes(results: NTOResults, spin: str, pair: int, spacing: float = DEFAULT_SPACING_BOHR,
                 isovalue: float = ORBITAL_ISOVALUE, progress: Optional[Callable[[float], None]] = None,
                 use_numba: Optional[bool] = None) -> Tuple[cf.CubeVolume, cf.CubeVolume]:
    """Hole and particle of one pair on one grid, every value computed exactly.

    The box is Multiwfn's (atoms +/- 1.7 r_vdW), grown until both orbitals are
    below isovalue / 10 on every face (molecularSurface.grid_box).
    """
    if not (spacing > 0 and np.isfinite(spacing)):
        raise ValueError(f"The grid spacing must be a positive number of bohr, not {spacing}.")
    if not (isovalue > 0 and np.isfinite(isovalue)):
        raise ValueError(f"The isovalue must be a positive number, not {isovalue}.")
    hole, particle, weight = results.pair_orbitals(spin, pair)
    wfn = results.wavefunction
    both = np.column_stack([hole, particle])

    def envelope(points):
        return np.abs(gw.orbital_values(wfn, both, points, use_numba=use_numba)).max(axis=1)

    origin, shape = ms.grid_box(wfn.coordinates, _vdw_radii_bohr(wfn), envelope, isovalue, spacing)
    check_grid_size(shape, spacing)
    nx, ny, nz = shape
    values = np.empty((2, nx, ny, nz))
    planes = max(1, SLAB_POINTS // (ny * nz))
    for x0 in range(0, nx, planes):
        x1 = min(nx, x0 + planes)
        points = ms.grid_points(origin + np.array([x0 * spacing, 0.0, 0.0]), spacing, (x1 - x0, ny, nz))
        slab = gw.orbital_values(wfn, both, points, use_numba=use_numba)
        values[:, x0:x1] = slab.T.reshape(2, x1 - x0, ny, nz)
        if progress is not None:
            progress(x1 / nx)
    return tuple(cf.CubeVolume(_orbital_title(results, spin, pair, kind, weight), CUBE_COMMENT,
                               np.asarray(wfn.atomic_numbers), np.asarray(wfn.nuclear_charges, dtype=float),
                               wfn.coordinates.copy(), np.asarray(origin, dtype=float).copy(),
                               spacing * np.eye(3), np.ascontiguousarray(values[k]))
                 for k, kind in enumerate(KINDS))


def cube_filename(stem: str, state_number: int, pair: int, kind: str, spin: Optional[str] = None) -> str:
    """<stem>_state<n>_NTO<k>_hole.cube; an unrestricted run inserts the spin before the kind."""
    spin_part = f"_{spin}" if spin else ""
    return f"{stem}_state{state_number}_NTO{pair}{spin_part}_{kind}.cube"


def nto_fchk_filename(stem: str, state_number: int) -> str:
    return f"{stem}_state{state_number}_NTO.fchk"


def report_filename(stem: str, state_number: int) -> str:
    return f"{stem}_state{state_number}_NTO.txt"


# ----------------------------------------------------------------------
# NTO formatted checkpoint (Multiwfn's layout)
# ----------------------------------------------------------------------
def _padded(values: np.ndarray, n: int) -> np.ndarray:
    out = np.zeros(n)
    out[:min(n, len(values))] = values[:n]
    return out


def nto_orbital_sections(results: NTOResults) -> Dict[str, np.ndarray]:
    """Replacement values of the orbital-energy and MO-coefficient sections.

    Multiwfn's order: the occupied NTOs by ascending weight (the dominant hole
    takes the HOMO's place), the virtual NTOs by descending weight (the
    dominant particle takes the LUMO's place); the "energies" are the weights.
    """
    sections = {}
    for spin, s in results.sets.items():
        n_vir = s.particles.shape[1]
        coefficients = np.hstack([s.holes[:, ::-1], s.particles])            # (n_basis, n_mo)
        energies = np.concatenate([_padded(s.weights, s.n_occ)[::-1], _padded(s.weights, n_vir)])
        name = "Alpha" if spin == ed.SPIN_ALPHA else "Beta"
        sections[f"{name} Orbital Energies"] = energies
        sections[f"{name} MO coefficients"] = coefficients.T.ravel()         # stored MO by MO
    return sections


def _fchk_value_lines(values: np.ndarray) -> List[str]:
    return ["".join(f"{v:16.8E}" for v in values[i:i + FCHK_VALUES_PER_LINE])
            for i in range(0, len(values), FCHK_VALUES_PER_LINE)]


def nto_fchk_text(results: NTOResults) -> str:
    """The input .fchk with only its orbital energies and MO coefficients replaced by the NTOs.

    Rotations among occupied and among virtual orbitals leave the SCF density
    unchanged, so every other section stays valid (GaussView, cubegen, Multiwfn).
    """
    try:
        lines = Path(results.fchk_path).read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise ValueError(f"Cannot read {results.fchk_path}: {exc}") from exc
    replacements = nto_orbital_sections(results)
    out: List[str] = []
    done = set()
    n = 0
    while n < len(lines):
        line = lines[n]
        out.append(line)
        n += 1
        name = line[:43].strip() if line[:1] != " " else ""
        if name in replacements:
            values = replacements[name]
            declared = int(line.split("N=")[1])
            if declared != values.size:
                raise ValueError(f"{Path(results.fchk_path).name}: '{name}' holds {declared} values, "
                                 f"but the NTOs need {values.size}.")
            while n < len(lines) and lines[n][:1] == " ":
                n += 1  # drop the old values
            out.extend(_fchk_value_lines(values))
            done.add(name)
    missing = sorted(set(replacements) - done)
    if missing:
        raise ValueError(f"{Path(results.fchk_path).name} lacks the section(s) {', '.join(missing)}.")
    return "\n".join(out) + "\n"


# ----------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------
def _composition_text(s: NTOSet, kind: str, pair: int) -> str:
    parts = [f"MO {mo} {100 * f:.1f}%" for mo, f in s.composition(kind, pair)]
    return ", ".join(parts) if parts else "no MO above 5%"


def report_text(results: NTOResults, written: Sequence[str] = ()) -> str:
    wfn, state = results.wavefunction, results.state
    listed = results.listed_pairs()
    total = results.weight_sum()
    lines = [
        "# gQTEA Natural Transition Orbitals",
        "# NTOs of one excited state from the excitation coefficients (Multiwfn's definition; de-excitations "
        "ignored).",
        f"# formatted checkpoint   {Path(results.fchk_path).name}",
        f"# excitation file        {Path(results.excitation_path).name} ({results.excitations.source})",
        f"# method / basis         {wfn.method} / {wfn.basis_set}",
        f"# reference              {'restricted' if results.restricted else 'unrestricted'}; "
        f"electrons {wfn.n_electrons} (alpha {wfn.n_alpha}, beta {wfn.n_beta}); "
        f"{wfn.n_basis} basis functions, {wfn.n_mo} MOs",
        f"# excited state          {state.title()}  multiplicity {state.multiplicity or 'undefined'}",
        f"# weight factor          {results.factor:g} (weight = factor x singular value^2)",
    ]
    norm = [("Sum of square of excitation coefficients", results.sum_excitation),
            ("Negative of the sum of square of de-excitation coefficients", -results.sum_deexcitation),
            ("Sum of above two values", results.sum_excitation - results.sum_deexcitation),
            (f"Deviation to expected normalization value ({results.expected_norm:g})", results.deviation),
            ("Sum of all NTO weights", total)]
    if not results.restricted:
        norm += [(f"Sum of {spin} NTO weights", results.weight_sum(spin)) for spin in SPIN_ORDER]
    lines += ["", "NORMALIZATION (as Multiwfn)"] + [f"  {label:<62}{value:10.6f}" for label, value in norm]
    lines += [
        "",
        f"NTO PAIRS (weight >= {LISTED_WEIGHT:g}; hole and particle MO numbers in the exported NTO .fchk)",
        "  pair  spin    weight         %   cum. %   hole  particle",
    ]
    cumulative = 0.0
    for p in listed:
        s = results.sets[p.spin]
        cumulative += p.weight
        lines.append(f"  {p.pair:2d}   {p.spin:<5}   {p.weight:.6f}  {100 * p.weight:6.2f}  {100 * cumulative:6.2f}"
                     f"  {s.n_occ + 1 - p.pair:5d}  {s.n_occ + p.pair:5d}")
    lines += ["", "COMPOSITION IN THE CANONICAL MOs (contributions of 5% or more)"]
    for p in listed:
        s = results.sets[p.spin]
        for kind in KINDS:
            lines.append(f"  pair {p.pair:2d} {p.spin:<5} {kind + ':':<9} {_composition_text(s, kind, p.pair)}")
    lines += ["", "NOTES"]
    lines += [f"  - {n}" for n in results.notes] or ["  none"]
    if written:
        lines += ["", "FILES"] + [f"  {w}" for w in written]
    return "\n".join(lines) + "\n"


def _write(path: Path, text: str) -> str:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"Could not write {path}: {exc}") from exc
    return str(path)


def _file_key(path) -> Tuple[str, int, int]:
    stat = Path(path).stat()
    return str(Path(path).resolve()), stat.st_mtime_ns, stat.st_size


class NTOAnalyser:
    """Reads an .fchk + excitation file pair (cached) and builds the NTOs of one state."""

    def __init__(self):
        self._loaded = None   # (key, (wfn, data, notes))

    def load(self, fchk_path, excitation_path):
        """(wavefunction, excitation file, notes); ValueError when the files do not belong together."""
        try:
            key = (_file_key(fchk_path), _file_key(excitation_path))
        except OSError as exc:
            raise ValueError(f"Cannot read {exc.filename}: {exc.strerror}") from exc
        if self._loaded is not None and self._loaded[0] == key:
            return self._loaded[1]
        wfn = gw.read_fchk(fchk_path)
        data = ed.read_excitation_file(excitation_path)
        notes = ed.check_against_wavefunction(data, wfn)
        self._loaded = (key, (wfn, data, notes))
        return self._loaded[1]

    def run(self, fchk_path, excitation_path, state_number: int,
            progress: Optional[Callable[[float, str], None]] = None) -> NTOResults:
        start = time.perf_counter()

        def report(fraction, message):
            if progress is not None:
                try:
                    progress(fraction, message)
                except Exception:
                    pass  # reporting must never abort the analysis

        report(0.05, "Reading the files...")
        wfn, data, load_notes = self.load(fchk_path, excitation_path)
        state = data.state(int(state_number))
        restricted = not wfn.is_unrestricted
        factor = 2.0 if restricted else 1.0
        spins = (ed.SPIN_ALPHA,) if restricted else SPIN_ORDER
        sets = {}
        for n, spin in enumerate(spins):
            report(0.3 + 0.6 * n / len(spins), f"Building the {spin} NTOs...")
            n_occ = wfn.n_alpha if spin == ed.SPIN_ALPHA else wfn.n_beta
            t = transition_matrix(state, spin, n_occ, wfn.n_mo)
            sets[spin] = natural_transition_orbitals(t, gw.mo_coefficients(wfn, spin), n_occ, factor, spin)
        x2, y2 = state.norm_terms()
        results = NTOResults(
            fchk_path=str(fchk_path), excitation_path=str(excitation_path), wavefunction=wfn,
            excitations=data, state=state, restricted=restricted, factor=factor, sets=sets,
            sum_excitation=x2, sum_deexcitation=y2, expected_norm=0.5 if restricted else 1.0,
            notes=list(load_notes))
        results.notes += self._notes(results)
        results.elapsed = time.perf_counter() - start
        report(1.0, "Done.")
        return results

    @staticmethod
    def _notes(r: NTOResults) -> List[str]:
        notes = []
        if abs(r.deviation) > NORMALIZATION_TOLERANCE:
            notes.append(f"The configuration coefficients sum to {r.sum_excitation - r.sum_deexcitation:.6f} "
                         f"instead of {r.expected_norm:g} (deviation {r.deviation:+.6f}): small coefficients are "
                         "probably missing. Rerun Gaussian with IOp(9/40=4); these NTOs are approximate.")
        if r.sum_deexcitation > 0:
            notes.append(f"TDHF/TD-DFT state: the de-excitation coefficients (sum of squares "
                         f"{r.sum_deexcitation:.6f}) are not part of the NTOs, as in Multiwfn and Gaussian, so "
                         f"the weights sum to {r.weight_sum():.6f}; a value well above 1 suggests using TDA.")
        top = r.listed_pairs()[0].weight
        if top < DOMINANT_WEIGHT:
            notes.append(f"No NTO pair dominates (largest weight {top:.4f}): this transition is not one hole -> "
                         "particle pair; discuss the listed pairs together.")
        return notes

    # --- outputs ------------------------------------------------------
    @staticmethod
    def write_report(results: NTOResults, output_dir: str) -> str:
        """<stem>_state<n>_NTO.txt in ``output_dir``; OSError becomes ValueError."""
        path = Path(output_dir) / report_filename(results.stem, results.state.number)
        return _write(path, report_text(results))

    @staticmethod
    def write_cubes(results: NTOResults, output_dir: str, pairs: Sequence[NTOPair],
                    spacing: float = DEFAULT_SPACING_BOHR, isovalue: float = ORBITAL_ISOVALUE,
                    progress: Optional[Callable[[float], None]] = None,
                    use_numba: Optional[bool] = None) -> List[str]:
        """Hole and particle cubes of each pair; returns the written paths in order."""
        written = []
        for n, p in enumerate(pairs):
            def step(fraction, n=n):
                if progress is not None:
                    progress((n + fraction) / len(pairs))
            volumes = pair_volumes(results, p.spin, p.pair, spacing, isovalue, step, use_numba)
            for kind, volume in zip(KINDS, volumes):
                name = cube_filename(results.stem, results.state.number, p.pair, kind,
                                     None if results.restricted else p.spin)
                path = Path(output_dir) / name
                try:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    written.append(cf.write_cube_volume(path, volume))
                except OSError as exc:
                    raise ValueError(f"Could not write {path}: {exc}") from exc
        return written

    @staticmethod
    def write_nto_fchk(results: NTOResults, path) -> str:
        return _write(Path(path), nto_fchk_text(results))
