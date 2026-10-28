"""Results of the Vanderbilt USPP program runatom: the .out report and the .uspp file.

No Toga here (a subprocess test pins it). read_uspp follows Utils/read_uspp.f
(files of version 7.3.2 or later): gfortran sequential unformatted records, each
framed by its 4-byte little-endian length before and after; integers are 4 bytes.
"""

import re
import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

EIGENVALUE_FLAG_RY = 1e-3      # Doc/TUTORIAL: a test is "probably OK" below 1 mRy
GHOST_TOLERANCE_RY = 0.1       # above the highest AE eigenvalue of a channel

_ROW = re.compile(r"^\s*(\d{3})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
_LROW = re.compile(r"^\s*l=\s*(\d+)\s+(.*)$")
_TOTAL = re.compile(r"total energy =\s*(-?\d+\.\d+)")
_ZV = re.compile(r"\|\s*z =\s*[-\d.]+\s+zv =\s*(-?\d+\.\d+)")   # the boxed pseudopotential report


@dataclass(frozen=True)
class EigenRow:
    nlm: int
    all_electron: float
    pseudo: float
    diff: float

    @property
    def l(self) -> int:
        return (self.nlm // 10) % 10


@dataclass
class OutReport:
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    comparison: List[EigenRow] = field(default_factory=list)   # last table in the file
    bessel: Dict[int, List[float]] = field(default_factory=dict)  # last US-Hamiltonian block
    total_energy: Optional[float] = None
    zv: Optional[float] = None


def _last_index(lines: Sequence[str], marker: str) -> Optional[int]:
    hits = [i for i, line in enumerate(lines) if marker in line]
    return hits[-1] if hits else None


def parse_out(text: str) -> OutReport:
    lines = text.splitlines()
    report = OutReport()
    for i, line in enumerate(lines):
        low = line.lower()
        if "***error" in low or "*** error" in low:
            report.errors.append(" ".join(l.strip() for l in lines[i:i + 3] if l.strip()))
        elif "warning" in low and line.strip() not in report.warnings:
            report.warnings.append(line.strip())
    start = _last_index(lines, "comparison of all-electron and pseudo eigenvalues")
    if start is not None:
        for line in lines[start + 1:]:
            match = _ROW.match(line)
            if match:
                report.comparison.append(EigenRow(int(match.group(1)), *map(float, match.group(2, 3, 4))))
            elif report.comparison and line.strip():
                break
    start = _last_index(lines, "eigenvalues of US hamiltonian")
    if start is not None:
        for line in lines[start + 1:]:
            match = _LROW.match(line)
            if not match:
                break
            report.bessel[int(match.group(1))] = [float(v) for v in match.group(2).split()]
    totals = _TOTAL.findall(text)
    report.total_energy = float(totals[-1]) if totals else None
    zv = _ZV.search(text)
    report.zv = float(zv.group(1)) if zv else None
    return report


def ghost_candidates(report: OutReport, tolerance: float = GHOST_TOLERANCE_RY) -> List[int]:
    """Channels with more bound Bessel-basis states than all-electron valence states (Doc/TUTORIAL)."""
    ghosts = []
    for l in sorted({row.l for row in report.comparison}):
        expected = [row.all_electron for row in report.comparison if row.l == l]
        found = report.bessel.get(l, [])
        if sum(1 for e in found if e < max(expected) + tolerance) > len(expected):
            ghosts.append(l)
    return ghosts


def uspp_records(data: bytes) -> List[bytes]:
    records, position = [], 0
    while position < len(data):
        index = len(records) + 1
        if position + 4 > len(data):
            raise ValueError(f"record {index}: truncated length marker")
        (length,) = struct.unpack_from("<i", data, position)
        end = position + 4 + length
        if length < 0 or end + 4 > len(data):
            raise ValueError(f"record {index}: length {length} runs past the end of the file")
        (trailer,) = struct.unpack_from("<i", data, end)
        if trailer != length:
            raise ValueError(f"record {index}: leading and trailing lengths differ ({length} vs {trailer}); "
                             "not a gfortran unformatted file")
        records.append(data[position + 4:end])
        position = end + 4
    return records


class _Record:
    def __init__(self, payload: bytes, index: int):
        self.payload, self.index, self.position = payload, index, 0

    def _take(self, dtype: str, count: int) -> np.ndarray:
        size = np.dtype(dtype).itemsize * count
        if self.position + size > len(self.payload):
            raise ValueError(f"record {self.index} is shorter than the USPP 7.3.6 layout expects")
        values = np.frombuffer(self.payload, dtype=dtype, count=count, offset=self.position).copy()
        self.position += size
        return values

    def ints(self, count: int) -> np.ndarray:
        return self._take("<i4", count)

    def reals(self, count: int) -> np.ndarray:
        return self._take("<f8", count)

    def int(self) -> int:
        return int(self.ints(1)[0])

    def real(self) -> float:
        return float(self.reals(1)[0])

    def text(self, count: int) -> str:
        if self.position + count > len(self.payload):
            raise ValueError(f"record {self.index} is shorter than the USPP 7.3.6 layout expects")
        value = self.payload[self.position:self.position + count].decode("latin-1")
        self.position += count
        return value

    def end(self):
        if self.position != len(self.payload):
            raise ValueError(f"record {self.index} is longer than the USPP 7.3.6 layout expects")


@dataclass
class UsppData:
    version: Tuple[int, int, int]
    date: Tuple[int, int, int]
    title: str
    z: float
    zv: float
    exfact: float
    valence: List[Tuple[int, float, float]]      # (nnlz, occupation, eigenvalue Ry)
    mesh: int
    etot: float
    keyps: int
    ifpcor: int
    nang: int
    lloc: int
    eloc: float
    ifqopt: int
    nqf: int
    qtryc: float
    rinner: np.ndarray
    irel: int
    rc: np.ndarray
    kkbeta: int
    l: np.ndarray
    eee: np.ndarray
    beta: np.ndarray            # (nbeta, kkbeta)
    ddd0: np.ndarray            # (nbeta, nbeta), symmetric
    ddd: np.ndarray
    qqq: np.ndarray
    qfunc: np.ndarray           # (nbeta, nbeta, kkbeta)
    qfcoef: np.ndarray          # (nbeta, nbeta, 2*nang-1, nqf)
    iptype: np.ndarray
    npf: int
    ptryc: float
    rcloc: float
    vloc0: np.ndarray
    rpcor: Optional[float]
    rspsco: Optional[np.ndarray]
    vloc: np.ndarray
    rsatom: np.ndarray
    r: np.ndarray
    rab: np.ndarray
    wf: np.ndarray              # (nchi, mesh)

    @property
    def nbeta(self) -> int:
        return len(self.l)


def read_uspp(path) -> UsppData:
    records = [_Record(p, i) for i, p in enumerate(uspp_records(Path(path).read_bytes()), 1)]
    stream = iter(records)

    def record() -> _Record:
        try:
            return next(stream)
        except StopIteration:
            raise ValueError("the file ends before the USPP 7.3.6 layout does") from None

    r = record(); version = tuple(int(v) for v in r.ints(3)); date = tuple(int(v) for v in r.ints(3)); r.end()
    if version < (7, 3, 2):
        raise ValueError(f"file version {version} predates the 7.3.2 layout this reader follows")
    r = record(); title = r.text(20).strip(); z, zv, exfact = r.real(), r.real(), r.real()
    nval, mesh = r.int(), r.int(); etot = r.real(); r.end()
    r = record(); valence = [(r.int(), r.real(), r.real()) for _ in range(nval)]; r.end()
    r = record(); keyps, ifpcor = r.int(), r.int(); r.real(); r.end()
    r = record(); nang, lloc, eloc, ifqopt, nqf, qtryc = r.int(), r.int(), r.real(), r.int(), r.int(), r.real(); r.end()
    n_l = 2 * nang - 1
    r = record(); rinner = r.reals(n_l); r.end()
    r = record(); irel = r.int(); r.end()
    r = record(); rc = r.reals(nang); r.end()
    r = record(); nbeta, kkbeta = r.int(), r.int(); r.end()
    l, eee = np.zeros(nbeta, int), np.zeros(nbeta)
    beta = np.zeros((nbeta, kkbeta))
    ddd0, ddd, qqq = (np.zeros((nbeta, nbeta)) for _ in range(3))
    qfunc = np.zeros((nbeta, nbeta, kkbeta))
    qfcoef = np.zeros((nbeta, nbeta, n_l, nqf))
    for j in range(nbeta):
        r = record(); l[j], eee[j] = r.int(), r.real(); beta[j] = r.reals(kkbeta); r.end()
        for k in range(j, nbeta):
            r = record()
            ddd0[j, k] = ddd0[k, j] = r.real()
            ddd[j, k] = ddd[k, j] = r.real()
            qqq[j, k] = qqq[k, j] = r.real()
            qfunc[j, k] = qfunc[k, j] = r.reals(kkbeta)
            qfcoef[j, k] = qfcoef[k, j] = r.reals(nqf * n_l).reshape(n_l, nqf)
            r.end()
    r = record(); iptype = r.ints(nbeta); npf, ptryc = r.int(), r.real(); r.end()
    r = record(); rcloc = r.real(); vloc0 = r.reals(mesh); r.end()
    rpcor = rspsco = None
    if ifpcor > 0:
        r = record(); rpcor = r.real(); r.end()
        r = record(); rspsco = r.reals(mesh); r.end()
    arrays = []
    for _ in range(4):                       # vloc, rsatom, r, rab
        r = record(); arrays.append(r.reals(mesh)); r.end()
    r = record(); nchi = r.int(); r.end()
    r = record(); wf = r.reals(mesh * nchi).reshape(nchi, mesh); r.end()
    if next(stream, None) is not None:
        raise ValueError("the file has more records than the USPP 7.3.6 layout")
    return UsppData(version, date, title, z, zv, exfact, valence, mesh, etot, keyps, ifpcor, nang,
                    lloc, eloc, ifqopt, nqf, qtryc, rinner, irel, rc, kkbeta, l, eee, beta, ddd0,
                    ddd, qqq, qfunc, qfcoef, iptype, npf, ptryc, rcloc, vloc0, rpcor, rspsco,
                    *arrays, wf)


def job_summary(stages: List[Tuple[str, float, OutReport]], uspp: Optional[UsppData] = None) -> str:
    """The Results tab: per stage time, warnings, eigenvalue table and ghost check."""
    lines = []
    for label, seconds, report in stages:
        lines.append(f"{label}: finished in {seconds:.2f} s")
        lines += [f"  WARNING: {w}" for w in report.warnings]
        if report.comparison:
            lines.append("    nlm     all-elec       pseudo         diff   (Ry)")
            for row in report.comparison:
                flag = "  <-- more than 1 mRy" if abs(row.diff) > EIGENVALUE_FLAG_RY else ""
                lines.append(f"    {row.nlm:3d} {row.all_electron:12.6f} {row.pseudo:12.6f} {row.diff:12.6f}{flag}")
        if report.comparison and report.bessel:
            ghosts = ghost_candidates(report)
            lines.append("  ghost check: " + (f"possible ghost state for l = {', '.join(map(str, ghosts))} "
                                              "(see 'eigenvalues of US hamiltonian' in the .out file)"
                                              if ghosts else "no extra bound state below the expected eigenvalues"))
        if report.total_energy is not None:
            lines.append(f"  total energy = {report.total_energy:.7f} Ry")
        lines.append("")
    if uspp is not None:
        lines.append(f"Pseudopotential file: {uspp.title}, zv = {uspp.zv:g}, nbeta = {uspp.nbeta}, "
                     f"mesh = {uspp.mesh}, kkbeta = {uspp.kkbeta}, "
                     f"rc = {', '.join(f'{v:g}' for v in uspp.rc)} bohr, "
                     f"core correction = {'yes' if uspp.ifpcor else 'no'}")
    return "\n".join(lines)
