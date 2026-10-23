"""The text report of the Excitation tab: descriptors with definitions and units, the classification with its
evidence, limitations and references. No Toga; imported lazily by excitationAnalysis.write_report."""

from __future__ import annotations

from pathlib import Path
from typing import List, Sequence

import numpy as np

import excitationAnalysis as xa
import excitationClassification as xc
import excitationData as ed

A = xa.BOHR_RADIUS_ANGSTROM
EV = xa.HARTREE_TO_EV
MO_CONTRIBUTION_THRESHOLD = 0.01

LIMITATIONS = [
    "Hole, electron and transition densities are the unrelaxed ones built from the configuration coefficients "
    "(Multiwfn's definition); relaxed densities (Z-vector, Gaussian density=) are not used.",
    "Gaussian prints the coefficients with 5 decimals and, without IOp(9/40=4), only those above 0.1; the "
    "normalisation check above shows whether the printed set is complete.",
    "Mulliken-like contributions are unreliable with diffuse basis functions and can be negative; Hirshfeld "
    "contributions need free-atom densities (H-Ar here) and all-electron basis sets.",
    "The classification summarises the descriptors with the thresholds quoted in its evidence; D, Sr and t follow "
    "published criteria, the Delta sigma and diffuse-electron thresholds are heuristics of this program. It is not a "
    "substitute for inspecting the hole and electron isosurfaces (Multiwfn manual 4.18.4).",
    "Fragment-based statements (IFCT) depend on how the fragments are defined; with the whole molecule as one "
    "fragment the CT percentage is zero by construction.",
    "Lambda and the Becke-grid indices converge slowly where orbitals have nodes; they are evaluated on 75 x 434 "
    "points per atom (Multiwfn's automatic 15 x 170 differs by up to ~1e-2 in Lambda).",
]


def _vec(v, digits: int = 4) -> str:
    return "  ".join(f"{x:{digits + 6}.{digits}f}" for x in v)


def _spin_tag(results, spin: str) -> str:
    return "" if results.restricted else ("A" if spin == ed.SPIN_ALPHA else "B")


def _references(results) -> List[str]:
    keys = ["multiwfn", "hole_electron", "le_bahers", "delta_r", "lambda", "ghost", "ghost_fix", "becke",
            "exciton", "mulliken" if results.partition == xa.PARTITION_MULLIKEN else "hirshfeld"]
    if results.ifct is not None:
        keys.append("ifct")
    return xc.references_for(keys)


def report_text(results, written: Sequence[str] = ()) -> str:
    r, wfn, state, m, g = results, results.wavefunction, results.state, results.moments, results.grid
    c = r.classification
    lines = [
        "# gQTEA Electron Excitation Analysis",
        "# Hole-electron analysis and related descriptors of one excited state (Multiwfn's definitions, main "
        "function 18), computed without Multiwfn.",
        f"# formatted checkpoint   {Path(r.fchk_path).name}",
        f"# excitation file        {Path(r.excitation_path).name} ({r.excitations.source})",
        f"# method / basis         {wfn.method} / {wfn.basis_set}",
        f"# reference              {'restricted (closed shell)' if r.restricted else 'unrestricted'}; electrons "
        f"{wfn.n_electrons} (alpha {wfn.n_alpha}, beta {wfn.n_beta}); {wfn.n_basis} basis functions, {wfn.n_mo} MOs"
        f"{'; diffuse functions present' if r.has_diffuse else ''}",
        f"# excited state          {state.title()}  multiplicity {state.multiplicity or 'undefined'}",
        "# evaluation             exact moment integrals (centroids, RMSD, D, H, t, transition dipole); Becke grid "
        "75 x 434 points per atom (Sm, Sr, HDI, EDI, Lambda, Hirshfeld)",
        f"# atom partition         {'Mulliken-like' if r.partition == xa.PARTITION_MULLIKEN else 'Hirshfeld'}",
    ]
    norm = [("Sum of square of excitation coefficients", r.sum_excitation),
            ("Negative of the sum of square of de-excitation coefficients", -r.sum_deexcitation),
            ("Sum of above two values", r.sum_excitation - r.sum_deexcitation),
            (f"Deviation to expected normalization value ({r.expected_norm:g})", r.deviation)]
    lines += ["", "NORMALIZATION (as Multiwfn)"] + [f"  {label:<62}{value:10.6f}" for label, value in norm]

    lines += ["", "MAJOR MO TRANSITIONS (contribution > 5 %; H = HOMO, L = LUMO)", f"  {r.major or 'none above 5 %'}"]

    lines += ["", f"MO CONTRIBUTIONS TO HOLE AND ELECTRON (>= {100 * MO_CONTRIBUTION_THRESHOLD:g} %)",
              "  MO         hole %   electron %"]
    total_h = sum(x.hole for x in r.mo_contributions)
    total_e = sum(x.electron for x in r.mo_contributions)
    for x in r.mo_contributions:
        if max(abs(x.hole), abs(x.electron)) >= MO_CONTRIBUTION_THRESHOLD:
            lines.append(f"  {x.mo:5d}{_spin_tag(r, x.spin):<2}{100 * x.hole:10.3f}{100 * x.electron:12.3f}")
    lines.append(f"  sum of all MOs  hole {100 * total_h:.3f} %   electron {100 * total_e:.3f} %")

    rows = [
        ("Integral of hole", f"{m.hole_integral:.6f}", "", "ideal 1"),
        ("Integral of electron", f"{m.electron_integral:.6f}", "", "ideal 1"),
        ("Integral of transition density", f"{m.transition_integral:.6f}", "", "ideal 0"),
        ("Sm index", f"{g.sm:.5f}", "", "integral of min(rho_hole, rho_ele): hole-electron overlap"),
        ("Sr index", f"{g.sr:.5f}", "", "integral of sqrt(rho_hole rho_ele): hole-electron overlap, 0 to 1"),
        ("Centroid of hole (X Y Z)", _vec(m.hole_centroid * A), "A", "integral of r rho_hole / integral of rho_hole"),
        ("Centroid of electron (X Y Z)", _vec(m.electron_centroid * A), "A", "same for the electron"),
        ("D_x D_y D_z", _vec(m.d_components * A, 3), "A", "|R_ele - R_hole| per axis"),
        ("D index", f"{m.d_index * A:.3f}", "A", "distance between the hole and electron centroids (CT length)"),
        ("Dipole change (X Y Z)", _vec(m.dipole_change), "a.u.", "-(R_ele - R_hole), unrelaxed"),
        ("RMSD of hole (X Y Z, norm)", _vec(np.append(m.sigma_hole, m.sigma_hole_norm) * A, 3), "A",
         "spatial extent of the hole"),
        ("RMSD of electron (X Y Z, norm)", _vec(np.append(m.sigma_electron, m.sigma_electron_norm) * A, 3), "A",
         "spatial extent of the electron"),
        ("Delta sigma (X Y Z, overall)", _vec(np.append(m.delta_sigma, m.delta_sigma_index) * A, 3), "A",
         "RMSD of electron minus RMSD of hole"),
        ("H_x H_y H_z", _vec(m.h * A, 3), "A", "mean RMSD of hole and electron per axis"),
        ("H_CT", f"{m.h_ct * A:.3f}", "A", "mean spread along the CT direction"),
        ("H index", f"{m.h_index * A:.3f}", "A", "(|sigma_hole| + |sigma_ele|) / 2"),
        ("t index", f"{m.t_index * A:.3f}", "A", "D - H_CT: > 0 when hole and electron are separated"),
        ("HDI", f"{g.hdi:.2f}", "", "hole delocalization index, 100 sqrt(int rho_hole^2): smaller = more delocalized"),
        ("EDI", f"{g.edi:.2f}", "", "electron delocalization index, same for the electron"),
    ]
    width = max(len(n) for n, *_ in rows)
    lines += ["", "HOLE-ELECTRON INDICES (Multiwfn manual 3.21.1.1)"]
    lines += [f"  {name:<{width}}  {value:>34} {unit:<5} {meaning}" for name, value, unit, meaning in rows]

    gh = r.ghost
    lines += ["", "ORBITAL-BASED INDICES",
              f"  Delta_r                {r.delta_r.value * A:10.4f} A    K^2-weighted distance between the "
              "centroids of the MOs of each pair (Guido et al.)",
              f"  Lambda                 {r.lambda_.value:10.4f}      K^2-weighted overlap of |phi_i| and |phi_a| "
              "(Peach et al.)",
              f"  Ghost-hunter index     {gh.value * EV:10.3f} eV   1st term {gh.term1 * EV:.3f} eV, 2nd term "
              f"{gh.term2 * EV:.3f} eV (1/D); excitation energy {state.energy_ev:.3f} eV",
              f"  Electron on diffuse basis functions {100 * r.diffuse_electron:.1f} %" if r.has_diffuse else
              "  The basis set has no diffuse functions."]

    mu = m.transition_dipole
    lines += ["", "TRANSITION ELECTRIC DIPOLE MOMENT (ground -> this state)",
              f"  calculated (a.u.)      {_vec(mu)}   |mu| {np.linalg.norm(mu):.4f}",
              f"  oscillator strength    calculated {m.oscillator(state.energy_ev):.4f}" +
              (f", printed {state.oscillator:.4f}" if state.oscillator is not None else "")]
    if r.log_dipole is not None:
        lines.append(f"  Gaussian log (a.u.)    {_vec(r.log_dipole)}   |mu| {np.linalg.norm(r.log_dipole):.4f}")

    atoms = r.atoms
    method = "Mulliken-like" if atoms.method == xa.PARTITION_MULLIKEN else "Hirshfeld"
    lines += ["", f"ATOM CONTRIBUTIONS ({method}; percent; overlap = sqrt(hole x electron); diff = electron - hole)",
              "  atom          hole %  electron %   overlap %     diff %"]
    for k, symbol in enumerate(wfn.atom_symbols()):
        lines.append(f"  {k + 1:4d}({symbol:<2})  {100 * atoms.hole[k]:9.2f}{100 * atoms.electron[k]:12.2f}"
                     f"{100 * atoms.overlap[k]:12.2f}{100 * atoms.difference[k]:11.2f}")

    if r.ifct is not None:
        f = r.ifct
        names = [x.name for x in f.fragments]
        lines += ["", "FRAGMENTS AND INTERFRAGMENT CHARGE TRANSFER (IFCT, Multiwfn manual 3.21.8)"]
        for x, h, e, d in zip(f.fragments, f.hole, f.electron, f.population_change):
            atoms_text = _ranges([a + 1 for a in x.atoms])
            lines.append(f"  {x.name:<12} atoms {atoms_text:<24} hole {100 * h:6.2f} %  electron {100 * e:6.2f} %  "
                         f"population change {d:+.5f}")
        lines.append("  Q[R, S] (electrons moved from R to S; diagonal = redistribution within R):")
        lines.append("  " + " " * 14 + "".join(f"{n:>14}" for n in names))
        for name, row in zip(names, f.matrix):
            lines.append(f"  {name:<14}" + "".join(f"{q:14.5f}" for q in row))
        for a, b, net in f.net_transfers():
            lines.append(f"  net {names[a]} -> {names[b]}: {net:.5f} e")
        lines.append(f"  intrinsic CT {f.intrinsic_ct:.3f} %, intrinsic LE {f.intrinsic_le:.3f} %"
                     + (f"; apparent CT {f.apparent_ct:.3f} %, apparent LE {100 - f.apparent_ct:.3f} %"
                        if f.apparent_ct is not None else ""))

    co = r.coulomb
    lines += ["", "HOLE-ELECTRON COULOMB ATTRACTION",
              f"  between C_hole and C_ele (Multiwfn's option 18)   {co.gaussian:.6f} a.u. ({co.gaussian * EV:.4f} eV)"]
    if co.exact is not None:
        lines.append(f"  between rho_hole and rho_ele (manual definition)  {co.exact:.6f} a.u. ({co.exact * EV:.4f} eV)"
                     f"; {xa.COULOMB_MODE_TEXT[co.mode]}")
    else:
        lines.append("  between rho_hole and rho_ele (manual definition)  not computed (tick 'Exact Coulomb "
                     "attraction - fast' or '- precise')")

    lines += ["", "CLASSIFICATION (descriptor-based)", f"  Verdict: {c.verdict}", f"  {c.rule}", "  Evidence:"]
    lines += [f"    - {e}" for e in c.evidence]
    if c.checks:
        lines += ["  Independent orbital-based checks:"] + [f"    - {x}" for x in c.checks]
    lines.append(f"  Confidence: {c.confidence}")
    lines.append(f"  Fragment-based assignment: {c.fragment_assignment or 'none (no fragments defined)'}")
    lines.append(f"  Orbital character: {c.orbital_note}")
    if c.warnings:
        lines += ["  Warnings:"] + [f"    - {w}" for w in c.warnings]

    lines += ["", "LIMITATIONS"] + [f"  - {x}" for x in LIMITATIONS]
    lines += ["", "REFERENCES"] + [f"  [{k}] {x}" for k, x in enumerate(_references(r), 1)]
    lines += ["", "NOTES"] + ([f"  - {n}" for n in r.notes] or ["  none"])
    if written:
        lines += ["", "FILES"] + [f"  {w}" for w in written]
    return "\n".join(lines) + "\n"


def _ranges(atoms: List[int]) -> str:
    """[1, 2, 3, 5] -> '1-3,5'."""
    parts, start = [], None
    for k, a in enumerate(atoms):
        if start is None:
            start = a
        if k + 1 == len(atoms) or atoms[k + 1] != a + 1:
            parts.append(str(start) if start == a else f"{start}-{a}")
            start = None
    return ",".join(parts)
