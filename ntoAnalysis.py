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
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

import excitationData as ed
import gaussianWavefunction as gw

LISTED_WEIGHT = 1e-4           # pairs listed in the table and the report (at least the largest)
NORMALIZATION_TOLERANCE = 0.01  # |deviation| above which the coefficients look incomplete
DOMINANT_WEIGHT = 0.5          # below this the largest pair does not "dominate"
COMPOSITION_THRESHOLD = 0.05   # MO contributions listed for each NTO
SPIN_ORDER = (ed.SPIN_ALPHA, ed.SPIN_BETA)


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
