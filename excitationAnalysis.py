"""Electron excitation analysis of a CIS / TDA / TD excited state (Multiwfn's main function 18).

Pure numpy/scipy, no Toga: the 3D Molecular Viewer's Excitation tab
(excitationTab) builds on it and the tests drive it headless. Lengths are bohr.

Multiwfn's definitions (manual section 3.21; excittrans.f90, read in the
2026.10.18 source), for configuration coefficients w (i -> a, matrix X) and
w' (i <- a, matrix Y) of one spin, f = 2 for a restricted reference (Gaussian
normalises its coefficients to 0.5) and 1 per spin otherwise:

    P_hole = f C_occ (X X^T - Y Y^T) C_occ^T      (hole density matrix)
    P_ele  = f C_vir (X^T X - Y^T Y) C_vir^T      (electron density matrix)
    P_T    = f C_occ (X + Y) C_vir^T              (transition density matrix)

Every quantity linear in these densities (integrals, centroids, RMSD, D,
H, t, transition dipole) is evaluated exactly from dipole and second-moment
integrals (espIntegrals.moment_matrices); the non-linear ones (Sm, Sr, HDI,
EDI, Lambda, Hirshfeld contributions) on a Becke molecular grid
(hirshfeldCharges.molecular_grid). Multiwfn evaluates all of them on a
uniform grid; it also drops configurations with |w| < 0.01 from the cross
terms, which ``cross_term_threshold`` reproduces for the validation tests.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

import cubeFile as cf
import espIntegrals as ei
import excitationData as ed
import gaussianWavefunction as gw
import hirshfeldCharges as hc
from indexRanges import parse_solute_index_ranges

HARTREE_TO_EV = 27.211386245988   # CODATA 2018
BOHR_RADIUS_ANGSTROM = cf.BOHR_RADIUS_ANGSTROM
SPIN_ORDER = (ed.SPIN_ALPHA, ed.SPIN_BETA)
MAJOR_TRANSITION_THRESHOLD = 0.05  # Multiwfn: 10 x compthres (0.5 %) = 5 %


# ----------------------------------------------------------------------
# Coefficients and density matrices
# ----------------------------------------------------------------------
@dataclass
class SpinCoefficients:
    spin: str
    n_occ: int
    x: np.ndarray   # (n_occ, n_vir) excitation coefficients, i -> a
    y: np.ndarray   # (n_occ, n_vir) de-excitation coefficients, i <- a


@dataclass
class DensityMatrices:
    """Hole, electron and transition density matrices in the fchk basis."""
    hole: np.ndarray
    electron: np.ndarray
    transition: np.ndarray   # not symmetric; real-space densities use (P + P^T) / 2

    @property
    def transition_symmetric(self) -> np.ndarray:
        return 0.5 * (self.transition + self.transition.T)


def occupied_count(wfn, spin: str) -> int:
    return wfn.n_alpha if spin == ed.SPIN_ALPHA else wfn.n_beta


def coefficient_matrices(state: ed.ExcitedState, wfn) -> Dict[str, SpinCoefficients]:
    """X and Y of every spin present in the reference (alpha only when restricted)."""
    spins = SPIN_ORDER if wfn.is_unrestricted else (ed.SPIN_ALPHA,)
    out = {}
    for spin in spins:
        n_occ = occupied_count(wfn, spin)
        n_vir = wfn.n_mo - n_occ
        out[spin] = SpinCoefficients(spin, n_occ, np.zeros((n_occ, n_vir)), np.zeros((n_occ, n_vir)))
    for c in state.configurations:
        if c.spin not in out:
            raise ValueError(f"Line {c.line}: a {c.spin} pair in a restricted reference.")
        s = out[c.spin]
        if not (1 <= c.occupied <= s.n_occ < c.virtual <= wfn.n_mo):
            raise ValueError(f"Line {c.line}: {c.occupied} -> {c.virtual} is not an occupied -> virtual pair "
                             f"for {c.spin} ({s.n_occ} occupied of {wfn.n_mo} orbitals).")
        target = s.y if c.deexcitation else s.x
        target[c.occupied - 1, c.virtual - s.n_occ - 1] = c.coefficient
    return out


def _gram(m: np.ndarray, kind: str, threshold: float) -> np.ndarray:
    """m m^T (hole, over occupied) or m^T m (electron, over virtual).

    With ``threshold`` > 0 the off-diagonal (cross) terms use only the entries
    with |w| >= threshold, the diagonal (local) terms all of them: Multiwfn's
    ``skippair`` screening in hole_electron / calc_holeele.
    """
    full = m @ m.T if kind == "hole" else m.T @ m
    if threshold <= 0.0:
        return full
    kept = np.where(np.abs(m) >= threshold, m, 0.0)
    cross = kept @ kept.T if kind == "hole" else kept.T @ kept
    return cross - np.diag(np.diag(cross)) + np.diag(np.diag(full))


def density_matrices(coeffs: Dict[str, SpinCoefficients], wfn, factor: float,
                     cross_term_threshold: float = 0.0) -> DensityMatrices:
    n = wfn.n_basis
    hole, electron, transition = np.zeros((n, n)), np.zeros((n, n)), np.zeros((n, n))
    for s in coeffs.values():
        c = gw.mo_coefficients(wfn, s.spin)
        co, cv = c[:, :s.n_occ], c[:, s.n_occ:]
        hole += factor * co @ (_gram(s.x, "hole", cross_term_threshold)
                               - _gram(s.y, "hole", cross_term_threshold)) @ co.T
        electron += factor * cv @ (_gram(s.x, "ele", cross_term_threshold)
                                   - _gram(s.y, "ele", cross_term_threshold)) @ cv.T
        transition += factor * co @ (s.x + s.y) @ cv.T
    return DensityMatrices(hole, electron, transition)


# ----------------------------------------------------------------------
# Analytic (moment) indices
# ----------------------------------------------------------------------
def basis_moments(wfn) -> Dict[str, np.ndarray]:
    """Overlap, dipole and diagonal second-moment integrals in the fchk basis (origin 0)."""
    t = wfn.transform
    return {key: t @ value @ t.T for key, value in ei.moment_matrices(wfn).items()}


def _trace(p: np.ndarray, m: np.ndarray) -> float:
    return float(np.sum(p * m))  # tr(P M) for symmetric M


@dataclass
class MomentIndices:
    """Multiwfn's hole-electron indices (manual 3.21.1.1, Theory 3), exact. Lengths bohr."""
    hole_integral: float
    electron_integral: float
    transition_integral: float
    hole_centroid: np.ndarray
    electron_centroid: np.ndarray
    sigma_hole: np.ndarray
    sigma_electron: np.ndarray
    transition_dipole: np.ndarray   # a.u., -integral of r T(r)

    @property
    def separation(self) -> np.ndarray:
        return self.electron_centroid - self.hole_centroid

    @property
    def d_components(self) -> np.ndarray:
        return np.abs(self.separation)

    @property
    def d_index(self) -> float:
        return float(np.linalg.norm(self.separation))

    @property
    def dipole_change(self) -> np.ndarray:
        """Change of the (unrelaxed) dipole moment, a.u.: -(R_ele - R_hole) x mean integral."""
        return -self.separation * 0.5 * (self.hole_integral + self.electron_integral)

    @property
    def sigma_hole_norm(self) -> float:
        return float(np.linalg.norm(self.sigma_hole))

    @property
    def sigma_electron_norm(self) -> float:
        return float(np.linalg.norm(self.sigma_electron))

    @property
    def delta_sigma(self) -> np.ndarray:
        return self.sigma_electron - self.sigma_hole

    @property
    def delta_sigma_index(self) -> float:
        return self.sigma_electron_norm - self.sigma_hole_norm

    @property
    def h(self) -> np.ndarray:
        return 0.5 * (self.sigma_electron + self.sigma_hole)

    @property
    def h_ct(self) -> float:
        """|(H_x D_x, H_y D_y, H_z D_z)| / D, as Multiwfn computes it (the manual prints |H . u_CT|)."""
        d = self.d_index
        return float(np.linalg.norm(self.h * self.separation) / d) if d > 0 else 0.0

    @property
    def h_index(self) -> float:
        return 0.5 * (self.sigma_electron_norm + self.sigma_hole_norm)

    @property
    def t_index(self) -> float:
        return self.d_index - self.h_ct

    def oscillator(self, energy_ev: float) -> float:
        mu = self.transition_dipole
        return float(2.0 / 3.0 * energy_ev / HARTREE_TO_EV * (mu @ mu))


def _centroid_and_sigma(p: np.ndarray, m: Dict[str, np.ndarray]) -> Tuple[float, np.ndarray, np.ndarray]:
    n = _trace(p, m["s"])
    centroid = np.array([_trace(p, m[k]) for k in "xyz"]) / n
    second = np.array([_trace(p, m[k + k]) for k in "xyz"]) / n
    return n, centroid, np.sqrt(np.maximum(second - centroid ** 2, 0.0))


def moment_indices(dm: DensityMatrices, moments: Dict[str, np.ndarray]) -> MomentIndices:
    n_h, r_h, s_h = _centroid_and_sigma(dm.hole, moments)
    n_e, r_e, s_e = _centroid_and_sigma(dm.electron, moments)
    t = dm.transition_symmetric
    return MomentIndices(n_h, n_e, _trace(t, moments["s"]), r_h, r_e, s_h, s_e,
                         -np.array([_trace(t, moments[k]) for k in "xyz"]))


# ----------------------------------------------------------------------
# Orbital-based quantities
# ----------------------------------------------------------------------
@dataclass
class MOContribution:
    spin: str
    mo: int          # 1-based, as in the fchk / Gaussian output
    hole: float      # fraction of the hole carried by this occupied MO
    electron: float  # fraction of the electron carried by this virtual MO


def mo_contributions(coeffs: Dict[str, SpinCoefficients], factor: float) -> List[MOContribution]:
    """Theta_i = f sum_a (w^2 - w'^2), Theta_a = f sum_i (w^2 - w'^2) (manual Theory 2); nonzero only."""
    out = []
    for spin in SPIN_ORDER:
        if spin not in coeffs:
            continue
        s = coeffs[spin]
        w2 = s.x ** 2 - s.y ** 2
        hole, electron = factor * w2.sum(axis=1), factor * w2.sum(axis=0)
        out += [MOContribution(spin, i + 1, float(v), 0.0) for i, v in enumerate(hole) if v != 0.0]
        out += [MOContribution(spin, s.n_occ + a + 1, 0.0, float(v)) for a, v in enumerate(electron) if v != 0.0]
    return out


def _orbital_label(kind: str, offset: int) -> str:
    if kind == "H":
        return "H" if offset == 0 else f"H{offset}"
    return "L" if offset == 0 else f"L+{offset}"


def major_transitions(state: ed.ExcitedState, wfn, threshold: float = MAJOR_TRANSITION_THRESHOLD) -> str:
    """Multiwfn's 'H-4 -> L 81.9%, H-4 -> L+2 12.1%' line (excittrans.f90, majorMOtrans)."""
    restricted = not wfn.is_unrestricted
    order = sorted(range(len(state.configurations)), key=lambda k: -abs(state.configurations[k].coefficient))
    parts = []
    for k in order:
        c = state.configurations[k]
        weight = c.coefficient ** 2 * (2 if restricted else 1)
        if not weight > threshold:
            continue
        homo = occupied_count(wfn, c.spin)
        spin = "" if restricted else ("a" if c.spin == ed.SPIN_ALPHA else "b")
        left = "H" + spin + (str(c.occupied - homo) if c.occupied != homo else "")
        right = "L" + spin + (f"+{c.virtual - homo - 1}" if c.virtual != homo + 1 else "")
        percent = f"{100 * weight:5.1f}".strip()
        arrow, sign = (" <- ", "-") if c.deexcitation else (" -> ", "")
        parts.append(f"{left}{arrow}{right} {sign}{percent}%")
    return ", ".join(parts)


@dataclass
class GhostHunter:
    value: float   # Hartree; an excitation energy below it suggests a ghost CT state
    term1: float   # sum w^2 (e_a - e_i) / sum w^2 over excitations
    term2: float   # 1 / D


def orbital_energies(wfn, spin: str) -> np.ndarray:
    energies = wfn.alpha_orbital_energies if spin == ed.SPIN_ALPHA else wfn.beta_orbital_energies
    if energies is None:
        raise ValueError(f"{wfn.path} has no {spin} orbital energies.")
    return energies


def ghost_hunter(state: ed.ExcitedState, wfn, d_bohr: float) -> GhostHunter:
    """M_AC of Campetella et al. (JCC 38, 2151) with the corrected form (JCP 154, 204102), as Multiwfn."""
    lines = [c for c in state.configurations if not c.deexcitation]
    w2 = np.array([c.coefficient ** 2 for c in lines])
    gaps = np.array([orbital_energies(wfn, c.spin)[c.virtual - 1] - orbital_energies(wfn, c.spin)[c.occupied - 1]
                     for c in lines])
    term1 = float(w2 @ gaps / w2.sum())
    term2 = 1.0 / d_bohr if d_bohr > 0 else float("inf")
    return GhostHunter(term1 - term2, term1, term2)


@dataclass
class OrbitalIndex:
    """Delta_r (bohr) or Lambda with the contribution of each excitation line (combined coefficient K)."""
    value: float
    pairs: List[Tuple[ed.Configuration, float, float]]   # (line, K, contribution)


def combined_coefficients(state: ed.ExcitedState) -> Tuple[List[Tuple[ed.Configuration, float]], float]:
    """Multiwfn's K = w + w' of the same pair and its normaliser N (delta_r / lambda_excit).

    N sums K^2 over the '->' lines and w'^2 over the '<-' lines (each '<-' line
    keeps its own coefficient in Multiwfn's list), as excittrans.f90 does.
    """
    back = {(c.occupied, c.virtual, c.spin): c.coefficient for c in state.configurations if c.deexcitation}
    combined, norm = [], 0.0
    for c in state.configurations:
        if c.deexcitation:
            norm += c.coefficient ** 2
            continue
        k = c.coefficient + back.get((c.occupied, c.virtual, c.spin), 0.0)
        combined.append((c, k))
        norm += k * k
    return combined, norm


def lambda_index(state: ed.ExcitedState, wfn, grid, use_numba: Optional[bool] = None) -> OrbitalIndex:
    """Lambda of Peach et al. (JCP 128, 044118): K^2-weighted overlap of |phi_i| and |phi_a| (Becke grid)."""
    combined, norm = combined_coefficients(state)
    overlaps = {}
    for spin in SPIN_ORDER:
        lines = [c for c, _ in combined if c.spin == spin]
        if not lines:
            continue
        occ = sorted({c.occupied for c in lines})
        vir = sorted({c.virtual for c in lines})
        c_mo = gw.mo_coefficients(wfn, spin)
        columns = c_mo[:, [m - 1 for m in occ + vir]]
        ovlp = np.zeros((len(occ), len(vir)))
        for k in range(0, len(grid.points), hc.POINT_CHUNK):
            values = np.abs(gw.orbital_values(wfn, columns, grid.points[k:k + hc.POINT_CHUNK], use_numba=use_numba))
            weighted = values[:, :len(occ)] * grid.weights[k:k + hc.POINT_CHUNK, None]
            ovlp += weighted.T @ values[:, len(occ):]
        overlaps[spin] = (dict(zip(occ, range(len(occ)))), dict(zip(vir, range(len(vir)))), ovlp)
    pairs, total = [], 0.0
    for c, k in combined:
        io, iv, ovlp = overlaps[c.spin]
        contribution = k * k / norm * float(ovlp[io[c.occupied], iv[c.virtual]])
        pairs.append((c, k, contribution))
        total += contribution
    return OrbitalIndex(total, pairs)


# ----------------------------------------------------------------------
# Becke-grid indices
# ----------------------------------------------------------------------
@dataclass
class GridIndices:
    """Non-linear hole-electron indices on the Becke grid (densities clipped at 0, as Multiwfn)."""
    sm: float               # integral of min(rho_hole, rho_ele)
    sr: float               # integral of sqrt(rho_hole rho_ele)
    hdi: float              # 100 sqrt(integral of rho_hole^2)
    edi: float              # 100 sqrt(integral of rho_ele^2)
    hole_integral: float    # quadrature checks of the analytic integrals
    electron_integral: float


def molecular_grid(wfn, use_numba: Optional[bool] = None):
    return hc.molecular_grid(wfn.coordinates, wfn.atom_symbols(), use_numba=use_numba)


def point_densities(wfn, dm: DensityMatrices, points: np.ndarray, use_numba: Optional[bool] = None
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """rho_hole, rho_ele and the transition density T at ``points`` (bohr), one AO pass."""
    t = wfn.transform
    cart = [t.T @ m @ t for m in (dm.hole, dm.electron, dm.transition_symmetric)]
    hole, electron, transition = gw.matrix_densities(wfn, cart, points, use_numba=use_numba)
    return hole, electron, transition


def grid_indices(wfn, dm: DensityMatrices, grid, use_numba: Optional[bool] = None) -> GridIndices:
    sums = np.zeros(6)
    for k in range(0, len(grid.points), hc.POINT_CHUNK):
        w = grid.weights[k:k + hc.POINT_CHUNK]
        hole, electron, _ = point_densities(wfn, dm, grid.points[k:k + hc.POINT_CHUNK], use_numba)
        h, e = np.maximum(hole, 0.0), np.maximum(electron, 0.0)
        sums += [w @ np.minimum(h, e), w @ np.sqrt(h * e), w @ (h * h), w @ (e * e), w @ hole, w @ electron]
    return GridIndices(float(sums[0]), float(sums[1]), 100.0 * float(np.sqrt(sums[2])),
                       100.0 * float(np.sqrt(sums[3])), float(sums[4]), float(sums[5]))


# ----------------------------------------------------------------------
# Atom and fragment contributions, IFCT
# ----------------------------------------------------------------------
PARTITION_MULLIKEN, PARTITION_HIRSHFELD = "mulliken", "hirshfeld"
PARTITIONS = (PARTITION_MULLIKEN, PARTITION_HIRSHFELD)


@dataclass
class AtomContributions:
    """Fractions of hole and electron on each atom (manual Theory 2), normalised to 1 as Multiwfn does."""
    method: str
    hole: np.ndarray
    electron: np.ndarray
    negative: bool          # a raw (Mulliken) contribution was negative

    @property
    def overlap(self) -> np.ndarray:
        """sqrt(hole x electron) where both are positive, else 0 (Multiwfn)."""
        both = (self.hole > 0) & (self.electron > 0)
        return np.where(both, np.sqrt(np.where(both, self.hole * self.electron, 0.0)), 0.0)

    @property
    def difference(self) -> np.ndarray:
        """Electron minus hole: the change of the atom's electron population."""
        return self.electron - self.hole


def basis_atoms(wfn) -> np.ndarray:
    """0-based atom of every fchk basis function."""
    return np.repeat([s.atom for s in wfn.shells], [s.n_functions for s in wfn.shells])


def _screened_gram(m: np.ndarray, kind: str, screen: float) -> np.ndarray:
    """sum over coefficient pairs sharing an orbital, keeping only |w_i w_j| >= screen (Multiwfn's 0.001)."""
    if kind == "ele":
        m = m.T  # electron: pairs share the occupied orbital, i.e. a row of X
    out = np.zeros((m.shape[0], m.shape[0]))
    for column in m.T:
        nz = np.flatnonzero(column)
        if nz.size == 0:
            continue
        outer = np.outer(column[nz], column[nz])
        out[np.ix_(nz, nz)] += np.where(np.abs(outer) >= screen, outer, 0.0)
    return out


def _normalised(raw_hole: np.ndarray, raw_electron: np.ndarray, method: str) -> AtomContributions:
    negative = bool((raw_hole < 0).any() or (raw_electron < 0).any())
    return AtomContributions(method, raw_hole / raw_hole.sum(), raw_electron / raw_electron.sum(), negative)


def mulliken_contributions(dm: DensityMatrices, overlap: np.ndarray, wfn, pair_screen: float = 0.0,
                           coeffs: Optional[Dict[str, SpinCoefficients]] = None,
                           factor: float = 2.0) -> AtomContributions:
    """Mulliken-like partition: Theta_A = sum_{mu on A} (P S)_mu,mu, then renormalised.

    ``pair_screen`` > 0 rebuilds P from ``coeffs`` dropping coefficient pairs with
    |w_i w_j| below it, as Multiwfn's atmcontri_holeele_Mulliken (0.001) does.
    """
    hole, electron = dm.hole, dm.electron
    if pair_screen > 0.0:
        if coeffs is None:
            raise ValueError("pair_screen needs the coefficient matrices.")
        n = wfn.n_basis
        hole, electron = np.zeros((n, n)), np.zeros((n, n))
        for s in coeffs.values():
            c = gw.mo_coefficients(wfn, s.spin)
            co, cv = c[:, :s.n_occ], c[:, s.n_occ:]
            hole += factor * co @ (_screened_gram(s.x, "hole", pair_screen)
                                   - _screened_gram(s.y, "hole", pair_screen)) @ co.T
            electron += factor * cv @ (_screened_gram(s.x, "ele", pair_screen)
                                       - _screened_gram(s.y, "ele", pair_screen)) @ cv.T
    atoms = basis_atoms(wfn)
    raw_hole = np.bincount(atoms, weights=np.einsum("ij,ji->i", hole, overlap), minlength=wfn.n_atoms)
    raw_electron = np.bincount(atoms, weights=np.einsum("ij,ji->i", electron, overlap), minlength=wfn.n_atoms)
    return _normalised(raw_hole, raw_electron, PARTITION_MULLIKEN)


def hirshfeld_contributions(wfn, dm: DensityMatrices, grid, use_numba: Optional[bool] = None) -> AtomContributions:
    """Hirshfeld partition: Theta_A = integral of w_A(r) rho(r), then renormalised. H-Ar, no ECP."""
    hc.check_hirshfeld_support(wfn)
    symbols = wfn.atom_symbols()
    raw = np.zeros((2, wfn.n_atoms))
    for k in range(0, len(grid.points), hc.POINT_CHUNK):
        pts = grid.points[k:k + hc.POINT_CHUNK]
        share = hc.hirshfeld_shares(pts, wfn.coordinates, symbols)
        hole, electron, _ = point_densities(wfn, dm, pts, use_numba)
        w = grid.weights[k:k + hc.POINT_CHUNK]
        raw[0] += (w * hole) @ share
        raw[1] += (w * electron) @ share
    return _normalised(raw[0], raw[1], PARTITION_HIRSHFELD)


@dataclass(frozen=True)
class Fragment:
    name: str
    atoms: Tuple[int, ...]   # 0-based
    automatic: bool = False  # the 'Rest' fragment added for atoms left out


def parse_fragments(text: str, n_atoms: int) -> Tuple[List[Fragment], List[str]]:
    """'1-10; 11-13' -> fragments (';' separates them, ranges as allBondAnalysis); atoms left out form 'Rest'."""
    text = (text or "").strip()
    if not text:
        return [], []
    fragments, owner = [], {}
    for n, part in enumerate([p for p in text.split(";") if p.strip()], 1):
        atoms = parse_solute_index_ranges(part)
        for a in atoms:
            if a > n_atoms:
                raise ValueError(f"Fragment {n}: atom {a} does not exist (the molecule has {n_atoms} atoms).")
            if a in owner:
                raise ValueError(f"Atom {a} is in fragments {owner[a]} and {n}.")
            owner[a] = n
        fragments.append(Fragment(f"Fragment {n}", tuple(a - 1 for a in atoms)))
    notes = []
    rest = tuple(a for a in range(n_atoms) if a + 1 not in owner)
    if rest:
        fragments.append(Fragment("Rest", rest, automatic=True))
        notes.append(f"Atoms not in any fragment ({len(rest)}) form the fragment 'Rest'.")
    return fragments, notes


@dataclass
class IFCTResult:
    """Interfragment charge transfer (Multiwfn manual 3.21.8)."""
    fragments: List[Fragment]
    hole: np.ndarray                # fragment fractions, negatives set to 0
    electron: np.ndarray
    matrix: np.ndarray              # Q[R, S] = hole_R electron_S: electron moved from R to S
    clamped: bool                   # a negative fragment contribution was set to 0

    @property
    def population_change(self) -> np.ndarray:
        """Delta p_R = sum_S Q[S, R] - sum_S Q[R, S]."""
        return self.matrix.sum(axis=0) - self.matrix.sum(axis=1)

    @property
    def intrinsic_ct(self) -> float:
        return 100.0 * float(self.matrix.sum() - np.trace(self.matrix))

    @property
    def intrinsic_le(self) -> float:
        return 100.0 * float(np.trace(self.matrix))

    @property
    def apparent_ct(self) -> Optional[float]:
        """|Delta p| of either fragment, two fragments only."""
        return 100.0 * abs(float(self.population_change[0])) if len(self.fragments) == 2 else None

    def net_transfers(self) -> List[Tuple[int, int, float]]:
        """(R, S, Q[R,S] - Q[S,R]) for R < S, largest net transfer first (0-based fragment indices)."""
        out = []
        for r in range(len(self.fragments)):
            for s in range(r + 1, len(self.fragments)):
                net = float(self.matrix[r, s] - self.matrix[s, r])
                out.append((r, s, net) if net >= 0 else (s, r, -net))
        return sorted(out, key=lambda item: -item[2])


def ifct(atoms: AtomContributions, fragments: List[Fragment]) -> IFCTResult:
    hole = np.array([atoms.hole[list(f.atoms)].sum() for f in fragments])
    electron = np.array([atoms.electron[list(f.atoms)].sum() for f in fragments])
    clamped = bool((hole < 0).any() or (electron < 0).any())
    hole, electron = np.maximum(hole, 0.0), np.maximum(electron, 0.0)
    return IFCTResult(list(fragments), hole, electron, np.outer(hole, electron), clamped)


def atom_ct_matrix(atoms: AtomContributions) -> np.ndarray:
    """Atom-atom charge-transfer matrix hole_A electron_B, negatives set to 0 (Multiwfn's atmCTmat.txt)."""
    return np.maximum(np.outer(atoms.hole, atoms.electron), 0.0)


def _mo_centroids(wfn, moments, spin: str, mos: List[int]) -> Dict[int, np.ndarray]:
    c = gw.mo_coefficients(wfn, spin)
    out = {}
    for mo in mos:
        v = c[:, mo - 1]
        out[mo] = np.array([v @ moments[k] @ v for k in "xyz"])
    return out


def delta_r(state: ed.ExcitedState, wfn, moments: Dict[str, np.ndarray]) -> OrbitalIndex:
    """Delta_r of Guido et al. (JCTC 9, 3118): K^2-weighted distance of the MO centroids (bohr)."""
    combined, norm = combined_coefficients(state)
    centroids = {}
    for spin in SPIN_ORDER:
        mos = sorted({m for c, _ in combined if c.spin == spin for m in (c.occupied, c.virtual)})
        if mos:
            centroids[spin] = _mo_centroids(wfn, moments, spin, mos)
    pairs, total = [], 0.0
    for c, k in combined:
        cen = centroids[c.spin]
        contribution = k * k / norm * float(np.linalg.norm(cen[c.virtual] - cen[c.occupied]))
        pairs.append((c, k, contribution))
        total += contribution
    return OrbitalIndex(total, pairs)
