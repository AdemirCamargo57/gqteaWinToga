"""QTAIM topology analysis of a Gaussian .fchk: run, report and 2D plane map.

``QTAIMAnalyser`` (headless, no Toga) reads the Total SCF Density of a
.fchk/.fch, finds its critical points and bond/ring paths with
qtaimTopology (Multiwfn's algorithms and defaults), and writes the
``<stem>_QTAIM.txt`` report in the layout of the Fukui tool's report. The
QTAIM tab of the 3D Molecular Viewer (qtaimTab) drives it.

Lengths are bohr internally; the report gives positions and lengths in
Angstrom and every property in atomic units.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

import gaussianWavefunction as gw
import qtaimTopology as qt
from qtaimTopology import BOHR_RADIUS_ANGSTROM, BCP, CCP, CP_TYPES, NCP, RCP

REPORT_SUFFIX = "_QTAIM.txt"
MAP_SUFFIX = "_QTAIM_map.json"
DEFAULT_MODES = (qt.MODE_NUCLEI, qt.MODE_PAIRS, qt.MODE_TRIANGLES, qt.MODE_TETRAHEDRA)

ProgressCallback = Callable[[float, str], None]


def resolve_output_dir(fchk_path: str, chosen: Optional[str]) -> str:
    """The chosen folder, or the folder of the .fchk file when blank."""
    chosen = (chosen or "").strip()
    return chosen if chosen else str(Path(fchk_path).resolve().parent)


@dataclass
class QTAIMResults:
    fchk_path: str
    wavefunction: gw.FchkWavefunction
    modes: Tuple[int, ...]
    search: qt.CPSearchParameters
    paths: qt.PathParameters
    critical_points: List[qt.CriticalPoint]
    bond_paths: List[qt.TopologyPath]
    ring_paths: List[qt.TopologyPath]
    notes: List[str] = field(default_factory=list)
    engine: str = "numpy"
    elapsed: float = 0.0

    def counts(self) -> Dict[str, int]:
        return qt.cp_counts(self.critical_points)

    def poincare_hopf(self) -> Tuple[int, bool]:
        return qt.poincare_hopf(self.counts())

    def cp(self, index: int) -> qt.CriticalPoint:
        return self.critical_points[index - 1]

    def paths_of(self, index: int) -> List[qt.TopologyPath]:
        return [p for p in self.bond_paths + self.ring_paths if p.start_cp == index]


def atom_label(symbols: Sequence[str], atom: int) -> str:
    return "NNA" if atom < 0 else f"{symbols[atom]}{atom + 1}"


def assignment_label(cp: qt.CriticalPoint, symbols: Sequence[str]) -> str:
    """O1 / NNA for a (3,-3), O1-H2 for a (3,-1) with a bond path, ring / cage, '-' otherwise."""
    a = cp.assignment
    if cp.kind == NCP:
        return atom_label(symbols, a[1]) if a and a[0] == "atom" else "NNA"
    if cp.kind == BCP:
        return "-".join(atom_label(symbols, i) for i in a[1:]) if a and a[0] == "pair" else "-"
    return "ring" if cp.kind == RCP else "cage"


class QTAIMAnalyser:
    """Headless QTAIM topology analysis of the Total SCF Density."""

    def run(self, fchk_path: str, modes: Sequence[int] = DEFAULT_MODES,
            search: Optional[qt.CPSearchParameters] = None, paths: Optional[qt.PathParameters] = None,
            use_numba: Optional[bool] = None, progress: Optional[ProgressCallback] = None) -> QTAIMResults:
        """Critical points, their properties, bond paths and (with cage points) ring paths."""
        def report(fraction: float, message: str):
            if progress is not None:
                try:
                    progress(fraction, message)
                except Exception:
                    pass  # reporting must never abort an analysis

        start = time.perf_counter()
        modes = tuple(sorted(set(int(m) for m in modes)))
        if not modes:
            raise ValueError("Select at least one CP search mode.")
        unknown = [m for m in modes if m not in qt.MODE_LABELS]
        if unknown:
            raise ValueError(f"Unknown CP search mode(s) {unknown}; use 2, 3, 4 or 5.")
        search = (search or qt.CPSearchParameters()).validate()
        paths = (paths or qt.PathParameters()).validate()
        report(0.0, f"Reading {Path(fchk_path).name}...")
        wfn = gw.read_fchk(fchk_path)
        engine = "numba (compiled)" if gw._resolve_numba(use_numba) else "numpy"

        report(0.02, "Searching critical points...")
        cps, notes = qt.find_critical_points(
            wfn, modes, search, use_numba=use_numba,
            progress=lambda f: report(0.02 + 0.73 * f, f"Searching critical points... {f:.0%}"))
        report(0.75, "Tracing bond paths...")
        bonds, path_notes = qt.bond_paths(wfn, cps, paths, use_numba=use_numba)
        notes += path_notes
        report(0.9, "Tracing ring paths...")
        rings = qt.ring_paths(wfn, cps, paths, use_numba=use_numba)

        results = QTAIMResults(fchk_path=str(fchk_path), wavefunction=wfn, modes=modes, search=search,
                               paths=paths, critical_points=cps, bond_paths=bonds, ring_paths=rings,
                               notes=notes, engine=engine)
        total, satisfied = results.poincare_hopf()
        if not satisfied:
            notes.append(f"The Poincare-Hopf sum is {total}, not 1: some critical points are missing. "
                         "Try more search modes, a smaller step-size scale or a trust radius.")
        ecp = [f"{s}{i + 1}" for i, (s, z, q) in enumerate(zip(wfn.atom_symbols(), wfn.atomic_numbers,
                                                                  wfn.nuclear_charges)) if abs(z - q) > 1e-6]
        if ecp:
            notes.append(f"ECP atoms ({', '.join(ecp)}): the .fchk holds no core density, so their (3,-3) "
                         "points and the properties near them are unreliable (Multiwfn adds a core density "
                         "from its EDF library; this tool does not).")
        if wfn.n_alpha != wfn.n_beta:
            notes.append("Open-shell wavefunction: the topology is that of the total density (alpha + beta).")
        results.elapsed = time.perf_counter() - start
        report(1.0, "Done.")
        return results

    # --- report ----------------------------------------------------------
    @staticmethod
    def _e(x: float) -> str:
        return f"{x:17.10E}"

    def report_text(self, results: QTAIMResults, written: Sequence[str] = ()) -> str:
        wfn = results.wavefunction
        symbols = wfn.atom_symbols()
        s, p = results.search, results.paths
        counts = results.counts()
        total, satisfied = results.poincare_hopf()
        trust = "none" if s.trust_radius is None else f"{s.trust_radius:g} bohr"
        lines = [
            "# gQTEA QTAIM Topology Analysis",
            "# Critical points and gradient paths of the electron density (Bader's QTAIM; Multiwfn's algorithms).",
            f"# file                   {Path(results.fchk_path).name}",
            f"# method / basis         {wfn.method} / {wfn.basis_set}",
            f"# density                {gw.DENSITY_SECTION}",
            f"# electrons              {wfn.n_electrons}   charge {wfn.charge}   multiplicity {wfn.multiplicity}",
            f"# atoms                  {wfn.n_atoms}",
            f"# search_modes           {', '.join(qt.MODE_LABELS[m] for m in results.modes)}",
            f"# cp_search              max_iterations {s.max_iterations}; step_scale {s.step_scale:g}; "
            f"gradient {s.gradient_tolerance:g}; displacement {s.displacement_tolerance:g};",
            f"#                        min_cp_distance {s.min_cp_distance:g} bohr; vdw_cutoff_scale "
            f"{s.vdw_cutoff_scale:g}; singular {s.singular_determinant:g}; trust_radius {trust}",
            f"# path_generation        max_points {p.max_points}; step {p.step:g} bohr; "
            f"stop_distance {p.stop_distance:g} bohr",
            f"# engine                 {results.engine}",
            "# units                  positions and lengths in Angstrom; properties in a.u.",
            "",
            "CRITICAL POINT SUMMARY",
            "  " + "   ".join(f"{t} {counts[t]}" for t in CP_TYPES),
            f"  Poincare-Hopf  {counts[NCP]} - {counts[BCP]} + {counts[RCP]} - {counts[CCP]} = {total}   "
            + ("satisfied" if satisfied else "NOT satisfied"),
            "",
            "CRITICAL POINTS",
            "   row  type          x (A)        y (A)        z (A)   assignment",
        ]
        for cp in results.critical_points:
            x, y, z = cp.position * BOHR_RADIUS_ANGSTROM
            lines.append(f"  {cp.index:4d}  {cp.kind}  {x:11.6f}  {y:11.6f}  {z:11.6f}   "
                         f"{assignment_label(cp, symbols)}")
        lines += ["", "DENSITY PROPERTIES AT CRITICAL POINTS (a.u.)",
                  "   row  type                  rho         laplacian           lambda1           lambda2"
                  "           lambda3   ellipticity        eta"]
        for cp in results.critical_points:
            q = cp.properties
            lines.append(f"  {cp.index:4d}  {cp.kind} " + " ".join(self._e(v) for v in (q.rho, q.laplacian,
                                                                                    *q.eigenvalues))
                         + f"  {q.ellipticity:12.6f} {q.eta:10.6f}")
        lines += ["", "ENERGY DENSITIES AT CRITICAL POINTS (a.u.; K = G - lap/4, V = -G - K, H = G + V = -K)",
                  "   row  type                    G                 K                 V                 H"
                  "        |V|/G         H/rho"]
        for cp in results.critical_points:
            q = cp.properties
            lines.append(f"  {cp.index:4d}  {cp.kind} " + " ".join(self._e(v) for v in (q.g, q.k, q.v, q.h))
                         + f"  {q.v_over_g:11.6f}  {q.h_over_rho:12.6f}")
        lines += ["", "BOND PATHS (each (3,-1) to its two (3,-3); lengths along the path; ratio = total / ncp_distance)",
                  "  path   BCP  atoms            length_1 (A)  length_2 (A)  total (A)  ncp_distance (A)  ratio"]
        number = 0
        for cp in results.critical_points:
            halves = sorted((b for b in results.bond_paths if b.start_cp == cp.index), key=lambda b: b.end_cp)
            if cp.kind != BCP or len(halves) != 2:
                continue
            number += 1
            l1, l2 = (h.length * BOHR_RADIUS_ANGSTROM for h in halves)
            ends = [results.cp(h.end_cp).position for h in halves]
            distance = float(np.linalg.norm(ends[0] - ends[1])) * BOHR_RADIUS_ANGSTROM
            ratio = (l1 + l2) / distance if distance > 0 else float("nan")
            lines.append(f"  {number:4d}  {cp.index:4d}  {assignment_label(cp, symbols):15s} {l1:12.6f}  "
                         f"{l2:12.6f} {l1 + l2:10.6f}  {distance:16.6f}  {ratio:6.4f}")
        if results.ring_paths:
            lines += ["", "RING PATHS ((3,+1) to (3,+3))", "  path   RCP   CCP   length (A)"]
            for k, r in enumerate(results.ring_paths, start=1):
                lines.append(f"  {k:4d}  {r.start_cp:4d}  {r.end_cp:4d}  {r.length * BOHR_RADIUS_ANGSTROM:11.6f}")
        lines.append("")
        for note in results.notes:
            lines.append(f"NOTE: {note}")
        for path in written:
            lines.append(f"Written: {path}")
        return "\n".join(lines) + "\n"

    def write_outputs(self, results: QTAIMResults, output_dir: str) -> List[str]:
        """Write <stem>_QTAIM.txt into ``output_dir``; returns the written paths."""
        report = Path(output_dir) / (Path(results.fchk_path).stem + REPORT_SUFFIX)
        try:
            report.parent.mkdir(parents=True, exist_ok=True)
            report.write_text(self.report_text(results, [str(report)]), encoding="utf-8")
        except OSError as exc:
            raise ValueError(f"Could not write {report}: {exc}") from exc
        return [str(report)]
