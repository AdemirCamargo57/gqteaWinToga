"""Gaussian formatted checkpoint (.fchk / .fch) reader, basis set and density.

Pure numpy, no Toga: the electrostatic surface potential tool builds on it
and the tests drive it headless. Every length is in bohr.

Basis conventions follow Multiwfn's fch reader (fileIO.f90 readfch), and
were verified against Gaussian 09 cubegen on water 5D/6D, cc-pVTZ (pure f),
cc-pVQZ (pure g), formaldehyde 10F (Cartesian f), ethylene (SP shells) and
the methyl radical (open shell):

* each Cartesian primitive is normalised on its own, x^i y^j z^k
  exp(-a r^2) with N = (2a/pi)^(3/4) (4a)^(L/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!),
  and multiplied by the fchk contraction coefficient as is (the contracted
  functions then come out normalised: tr(P S) = N_electrons to 1e-8);
* an SP shell (type -1) is an s shell plus a p shell sharing exponents, the
  p part taking ``P(S=P) Contraction coefficients``;
* a pure shell (type <= -2) is mapped onto its Cartesian components by
  ``pure_to_cartesian``, so all integrals run in one Cartesian code path;
* ``Total SCF Density`` (lower triangle, row-wise) is the density used.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import comb, pi, sqrt
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from help import AtomicData

# Cartesian component order written to the fchk for each angular momentum
# (Multiwfn define.f90: s, p, d as Gaussian pop=full; f and g in the fch
# order, which is not the .wfn order). Entries are powers (i, j, k) of x, y, z.
FCHK_CARTESIAN_ORDER: Dict[int, Tuple[Tuple[int, int, int], ...]] = {
    0: ((0, 0, 0),),
    1: ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    2: ((2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)),
    3: ((3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0),
        (2, 0, 1), (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)),
    4: ((0, 0, 4), (0, 1, 3), (0, 2, 2), (0, 3, 1), (0, 4, 0),
        (1, 0, 3), (1, 1, 2), (1, 2, 1), (1, 3, 0), (2, 0, 2),
        (2, 1, 1), (2, 2, 0), (3, 0, 1), (3, 1, 0), (4, 0, 0)),
}
MAX_ANGULAR_MOMENTUM = max(FCHK_CARTESIAN_ORDER)
DENSITY_SECTION = "Total SCF Density"
REQUIRED_SECTIONS = (
    "Atomic numbers", "Nuclear charges", "Current cartesian coordinates",
    "Shell types", "Number of primitives per shell", "Shell to atom map",
    "Primitive exponents", "Contraction coefficients", DENSITY_SECTION,
)
SYMBOL_BY_ATOMIC_NUMBER = {z: s for s, z in AtomicData.atomic_numbers.items()}


def cartesian_components(l: int) -> Tuple[Tuple[int, int, int], ...]:
    return FCHK_CARTESIAN_ORDER[l]


def double_factorial(n: int) -> int:
    """n!! with (-1)!! = 0!! = 1."""
    result = 1
    while n > 1:
        result *= n
        n -= 2
    return result


def component_angular_norm(powers: Tuple[int, int, int]) -> float:
    """1/sqrt((2i-1)!!(2j-1)!!(2k-1)!!), the component part of a Cartesian norm."""
    i, j, k = powers
    return 1.0 / sqrt(double_factorial(2 * i - 1) * double_factorial(2 * j - 1)
                      * double_factorial(2 * k - 1))


def primitive_norm(alpha: np.ndarray, l: int) -> np.ndarray:
    """(2a/pi)^(3/4) (4a)^(l/2), the exponent part of a Cartesian norm."""
    return (2.0 * alpha / pi) ** 0.75 * (4.0 * alpha) ** (l / 2.0)


def _same_centre_overlap(p1, p2) -> float:
    """Overlap of two normalised Cartesian components with one centre and exponent."""
    value = 1.0
    for a, b in zip(p1, p2):
        if (a + b) % 2:
            return 0.0
        value *= double_factorial(a + b - 1)
    return value * component_angular_norm(p1) * component_angular_norm(p2)


def _solid_harmonic_monomials(l: int, m: int) -> Dict[Tuple[int, int, int], float]:
    """Real solid harmonic S_lm as {(i, j, k): coefficient} over x^i y^j z^k.

    Helgaker, Jorgensen & Olsen, Molecular Electronic-Structure Theory,
    eq. 6.4.48; the overall scale is irrelevant (rows are renormalised).
    """
    am = abs(m)
    vm = 0.0 if m >= 0 else 0.5
    terms: Dict[Tuple[int, int, int], float] = {}
    v_values = [vm + n for n in range(int(np.floor(am / 2.0 - vm)) + 1)]
    for t in range((l - am) // 2 + 1):
        for u in range(t + 1):
            for v in v_values:
                two_v = int(round(2 * v))
                coef = ((-1) ** int(round(t + v - vm)) * (0.25 ** t) * comb(l, t)
                        * comb(l - t, am + t) * comb(t, u) * comb(am, two_v))
                key = (2 * t + am - 2 * u - two_v, 2 * u + two_v, l - 2 * t - am)
                terms[key] = terms.get(key, 0.0) + coef
    return terms


def pure_to_cartesian(l: int) -> np.ndarray:
    """T with pure_p = sum_c T[p, c] cart_c (both normalised).

    Rows follow Gaussian's pure order 0, +1, -1, +2, -2, ...; columns the
    fchk Cartesian order of ``cartesian_components(l)``.
    """
    comps = cartesian_components(l)
    index = {c: n for n, c in enumerate(comps)}
    order = [0]
    for m in range(1, l + 1):
        order += [m, -m]
    overlap = np.array([[_same_centre_overlap(a, b) for b in comps] for a in comps])
    rows = []
    for m in order:
        row = np.zeros(len(comps))
        for powers, coef in _solid_harmonic_monomials(l, m).items():
            row[index[powers]] += coef / component_angular_norm(powers)
        row /= sqrt(row @ overlap @ row)
        rows.append(row)
    return np.array(rows)


def read_fchk_sections(text: str) -> Dict[str, object]:
    """Every integer/real scalar and array section of an fchk, keyed by name.

    A header line has the name in columns 1-43, the type letter in column 44
    and either a scalar or ``N=  count``; array values follow on lines that
    start with a blank. Character and logical sections are skipped.
    """
    lines = text.splitlines()
    if len(lines) < 3:
        raise ValueError("The file is too short to be a Gaussian formatted checkpoint.")
    sections: Dict[str, object] = {}
    n = 2  # line 1 = title, line 2 = job type / method / basis
    while n < len(lines):
        line = lines[n]
        n += 1
        if not line.strip() or line[0] == " " or len(line) < 45:
            continue
        name, kind, rest = line[:43].strip(), line[43], line[44:].split()
        if not rest:
            continue
        if rest[0].startswith("N="):
            count = int(rest[0][2:]) if len(rest[0]) > 2 else int(rest[1])
            tokens: List[str] = []
            while n < len(lines) and lines[n][:1] == " ":
                tokens.extend(lines[n].split())
                n += 1
            try:
                if kind == "I":
                    values = np.array([int(t) for t in tokens], dtype=np.int64)
                elif kind == "R":
                    values = np.array([float(t.replace("D", "E")) for t in tokens])
                else:
                    continue
            except ValueError as exc:
                raise ValueError(f"Section '{name}' contains a value that is not a number.") from exc
            if values.size != count:
                raise ValueError(
                    f"Section '{name}' declares {count} values but holds {values.size}; "
                    "the file looks truncated."
                )
            sections[name] = values
        elif kind == "I":
            sections[name] = int(rest[0])
        elif kind == "R":
            sections[name] = float(rest[0].replace("D", "E"))
    return sections


@dataclass
class Shell:
    l: int
    pure: bool
    center: np.ndarray       # bohr
    atom: int                # 0-based atom index
    exponents: np.ndarray
    coefficients: np.ndarray  # fchk contraction coefficients, as written

    @property
    def n_cart(self) -> int:
        return (self.l + 1) * (self.l + 2) // 2

    @property
    def n_functions(self) -> int:
        return 2 * self.l + 1 if self.pure else self.n_cart


@dataclass
class FchkWavefunction:
    path: str
    title: str
    atomic_numbers: np.ndarray
    nuclear_charges: np.ndarray   # effective charges (differ from Z with an ECP)
    coordinates: np.ndarray       # (n_atoms, 3) bohr
    shells: List[Shell]
    density: np.ndarray           # symmetric, in the fchk basis
    n_electrons: int
    charge: int
    multiplicity: int
    dipole: Optional[np.ndarray] = None
    transform: np.ndarray = field(default=None, repr=False)     # fchk -> Cartesian
    cart_density: np.ndarray = field(default=None, repr=False)  # T^T P T

    @property
    def n_atoms(self) -> int:
        return len(self.atomic_numbers)

    def atom_symbols(self) -> List[str]:
        return [SYMBOL_BY_ATOMIC_NUMBER.get(int(z), "X") for z in self.atomic_numbers]


def read_fchk(path) -> FchkWavefunction:
    """Parse a .fchk/.fch file. Every problem is a ValueError naming it."""
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ValueError(f"Cannot read {path}: {exc}") from exc
    sections = read_fchk_sections(text)
    missing = [name for name in REQUIRED_SECTIONS if name not in sections]
    if missing:
        raise ValueError(f"{path.name} is missing fchk section(s): {', '.join(missing)}.")

    coords = np.asarray(sections["Current cartesian coordinates"], dtype=float).reshape(-1, 3)
    z = np.asarray(sections["Atomic numbers"])
    if len(z) != len(coords):
        raise ValueError(f"{path.name}: {len(z)} atomic numbers but {len(coords)} coordinate triples.")
    shell_types = sections["Shell types"]
    nprim = sections["Number of primitives per shell"]
    shell_atom = sections["Shell to atom map"]
    exps = sections["Primitive exponents"]
    coefs = sections["Contraction coefficients"]
    sp_coefs = sections.get("P(S=P) Contraction coefficients")

    shells: List[Shell] = []
    start = 0
    for stype, count, atom in zip(shell_types.tolist(), nprim.tolist(), shell_atom.tolist()):
        sl = slice(start, start + count)
        start += count
        centre = coords[atom - 1]
        if stype == -1:
            if sp_coefs is None:
                raise ValueError(f"{path.name}: SP shell without 'P(S=P) Contraction coefficients'.")
            shells.append(Shell(0, False, centre, atom - 1, exps[sl].copy(), coefs[sl].copy()))
            shells.append(Shell(1, False, centre, atom - 1, exps[sl].copy(), sp_coefs[sl].copy()))
            continue
        l = abs(stype)
        if l > MAX_ANGULAR_MOMENTUM:
            raise ValueError(
                f"{path.name}: basis functions with l = {l} are not supported "
                f"(highest supported is g, l = {MAX_ANGULAR_MOMENTUM})."
            )
        shells.append(Shell(l, stype <= -2, centre, atom - 1, exps[sl].copy(), coefs[sl].copy()))
    if start != exps.size:
        raise ValueError(f"{path.name}: shell primitive counts do not match the exponent list.")

    nbf = sum(s.n_functions for s in shells)
    tri = sections[DENSITY_SECTION]
    if tri.size != nbf * (nbf + 1) // 2:
        raise ValueError(
            f"{path.name}: '{DENSITY_SECTION}' has {tri.size} elements, "
            f"expected {nbf * (nbf + 1) // 2} for {nbf} basis functions."
        )
    density = np.zeros((nbf, nbf))
    density[np.tril_indices(nbf)] = tri
    density = density + density.T - np.diag(np.diag(density))

    lines = text.splitlines()
    wfn = FchkWavefunction(
        path=str(path),
        title=lines[0].strip() if lines else "",
        atomic_numbers=z,
        nuclear_charges=np.asarray(sections["Nuclear charges"], dtype=float),
        coordinates=coords,
        shells=shells,
        density=density,
        n_electrons=int(sections.get("Number of electrons", 0)),
        charge=int(sections.get("Charge", 0)),
        multiplicity=int(sections.get("Multiplicity", 1)),
        dipole=np.asarray(sections["Dipole Moment"]) if "Dipole Moment" in sections else None,
    )
    wfn.transform = basis_transform(shells)
    wfn.cart_density = wfn.transform.T @ density @ wfn.transform
    return wfn


def basis_transform(shells: List[Shell]) -> np.ndarray:
    """Block-diagonal (n_functions, n_cartesian) map onto normalised Cartesians."""
    nbf = sum(s.n_functions for s in shells)
    ncart = sum(s.n_cart for s in shells)
    t = np.zeros((nbf, ncart))
    row = col = 0
    for s in shells:
        block = pure_to_cartesian(s.l) if s.pure else np.eye(s.n_cart)
        t[row:row + s.n_functions, col:col + s.n_cart] = block
        row += s.n_functions
        col += s.n_cart
    return t


def evaluate_cartesian_aos(wfn: FchkWavefunction, points: np.ndarray) -> np.ndarray:
    """(n_points, n_cartesian) values of every normalised Cartesian component."""
    points = np.asarray(points, dtype=float)
    out = np.empty((points.shape[0], sum(s.n_cart for s in wfn.shells)))
    col = 0
    for s in wfn.shells:
        d = points - s.center
        r2 = np.einsum("ij,ij->i", d, d)
        radial = np.exp(-np.outer(r2, s.exponents)) @ (s.coefficients * primitive_norm(s.exponents, s.l))
        for c in cartesian_components(s.l):
            out[:, col] = (radial * component_angular_norm(c)
                           * d[:, 0] ** c[0] * d[:, 1] ** c[1] * d[:, 2] ** c[2])
            col += 1
    return out


def electron_density(wfn: FchkWavefunction, points: np.ndarray, chunk: int = 20000) -> np.ndarray:
    """rho(r) = sum_ab P_ab chi_a(r) chi_b(r) at every point (a.u.)."""
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    rho = np.empty(points.shape[0])
    for a in range(0, points.shape[0], chunk):
        ao = evaluate_cartesian_aos(wfn, points[a:a + chunk])
        rho[a:a + chunk] = np.einsum("ij,ij->i", ao @ wfn.cart_density, ao)
    return rho
