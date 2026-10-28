"""Input files (.adat) of the Vanderbilt USPP program runatom (USPP 7.3.6).

No Toga here (a subprocess test pins it). runatom reads every input line with a
fixed-column Fortran format (rwsubs.f, subroutine readin). gfortran ignores
blanks inside a numeric field and scales an F/E field that has no decimal point
by 10**-d: "1 0 0" in an i5 field is 100 and "2" in an f10.5 field is 2e-5.
Text typed by a user is therefore never copied into a file: values are parsed,
checked and written column-exactly here.
"""

import math
import re
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

LABEL_COLUMN = 38           # labels start here, or 2 columns after the data
MAX_FORTRAN_NAME = 40       # runatom's file-name arguments are character*40
MAX_REFORM_NAME = 36        # reform.f reads its file names with format a36
SAFE_NAME = re.compile(r"^[A-Za-z0-9._-]+$")

_FORMAT_ITEM = re.compile(r"^(\d*)([iafe])(\d+)(?:\.(\d+))?$")
_FORTRAN_REAL = re.compile(r"^([+-]?)(\d*)(\.?)(\d*)(?:[ED]([+-]?\d+)|([+-]\d+))?$")


def parse_format(fmt: str) -> List[Tuple[str, int, int]]:
    """'2i5,f10.5' -> [('i', 5, 0), ('i', 5, 0), ('f', 10, 5)]."""
    fields = []
    for item in fmt.split(","):
        match = _FORMAT_ITEM.match(item.strip())
        if match is None:
            raise ValueError(f"unsupported Fortran format item {item!r}")
        repeat = int(match.group(1) or 1)
        fields.extend([(match.group(2), int(match.group(3)), int(match.group(4) or 0))] * repeat)
    return fields


def read_field(kind: str, width: int, decimals: int, raw: str):
    """One input field as gfortran reads it (BLANK='NULL': blanks ignored)."""
    if kind == "a":
        return raw
    text = raw.replace(" ", "")
    if kind == "i":
        if not text:
            return 0
        if not re.fullmatch(r"[+-]?\d+", text):
            raise ValueError(f"Bad integer {raw.strip()!r}")
        return int(text)
    if not text:
        return 0.0
    match = _FORTRAN_REAL.match(text.upper())
    if match is None or not (match.group(2) or match.group(4)):
        raise ValueError(f"Bad value during floating point read: {raw.strip()!r}")
    sign, whole, point, frac, exp_e, exp_sign = match.groups()
    exponent = int(exp_e or exp_sign or 0)
    if point:
        return float(f"{sign}{whole or '0'}.{frac or '0'}e{exponent}")
    return float(f"{sign}{whole}e{exponent - decimals}")   # implied decimal point


def read_record(line: str, fmt: str, count: Optional[int] = None) -> list:
    """The first `count` fields of one record (a short record is padded with blanks)."""
    fields = parse_format(fmt)
    if count is not None:
        fields = fields[:count]     # a READ stops when its list is exhausted
    line = line.rstrip("\r\n")
    values, position = [], 0
    for kind, width, decimals in fields:
        values.append(read_field(kind, width, decimals, line[position:position + width].ljust(width)))
        position += width
    return values


def _float_texts(value: float):
    """Spellings of a float, most precise first; every one has a decimal point."""
    text = repr(value)
    if "e" in text:
        mantissa, exponent = text.split("e")
        yield f"{mantissa if '.' in mantissa else mantissa + '.0'}E{int(exponent)}"
    else:
        yield text
    if value == 0.0 or 1e-3 <= abs(value) < 1e15:
        for decimals in range(12, 0, -1):
            yield f"{value:.{decimals}f}"
    for digits in range(12, 0, -1):
        yield f"{value:.{digits}E}"


def write_field(kind: str, width: int, decimals: int, value) -> str:
    """One output field; floats keep one leading blank unless only the exact value fits."""
    if kind == "a":
        text = str(value)
        if len(text) > width:
            raise ValueError(f"{text!r} is longer than {width} characters")
        return text.ljust(width)
    if kind == "i":
        if isinstance(value, float) and not value.is_integer():
            raise ValueError(f"{value} is not an integer")
        text = str(int(value))
        if len(text) > width:
            raise ValueError(f"{value} does not fit in {width} columns")
        return text.rjust(width)
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{value!r} is not a finite number")
    for text in _float_texts(number):
        if len(text) < width or (len(text) == width and float(text) == number):
            return text.rjust(width)
    raise ValueError(f"{value} does not fit in {width} columns")


def format_line(fmt: str, values: Sequence, label: str) -> List[str]:
    """Data written column-exactly, then '<label> (<fmt>)' after the columns read.

    More values than the format has fields continue on a new line, as Fortran
    format reversion does when it reads them.
    """
    fields = parse_format(fmt)
    lines = []
    for start in range(0, max(len(values), 1), len(fields)):
        chunk = values[start:start + len(fields)]
        lines.append("".join(write_field(k, w, d, v) for (k, w, d), v in zip(fields, chunk)))
    full_width = sum(width for _, width, _ in fields)   # every column one line of the format reads
    column = max(LABEL_COLUMN, full_width + 2, len(lines[0]) + 2)
    lines[0] = lines[0].ljust(column) + f"{label} ({fmt})"
    return lines


LINE_FORMATS = {
    "controls": ("6i5", "ifae,ifpsp,ifprt,ifplw,ilogd"),
    "logd": ("3f10.5,i5", "rlogd,emin,emax,nnt"),
    "convergence": ("2e10.1,f10.5,i5", "thresh,tol,damp,maxit"),
    "title": ("a20", "title"),
    "atom": ("f5.0,2f10.5", "z,xion,exfact"),
    "mesh": ("3f10.5", "rmax,aasf,bbsf"),
    "nstates": ("2i5", "ncspvs,irel"),
    "shell": ("i4,f7.3,f14.6", "nnlz,wwnl,ee"),
    "partition": ("3i5", "ncores,nvales,nang"),
    "bessel": ("4f10.5", "besrmax,besemin,besemax,besde"),
    "keyps": ("2i5,2f10.5", "keyps,ifpcor,rinner,rpcor"),
    "rinner": ("5f10.5", "rinner for L=1..2*nang-1"),
    "nbeta": ("i5,f10.5", "nbeta,rcloc"),
    "rc": ("4f10.5", "rc for l=0..nang-1"),
    "beta": ("2i5,f10.5,i5", "lll,keyee,eeread,iptype"),
    "npf": ("i5,f10.5", "npf,ptryc"),
    "local": ("2i5,f10.5,i5", "lloc,keyee,eloc,iploctype"),
    "qfunc": ("2i5,f10.5,i5", "ifqopt,nqf,qtryc,nfix"),
    "qfix": ("i5,2f10.5", "ibfix,rfix,qfix"),
}


@dataclass
class Controls:
    """The four lines every .adat file starts with."""
    ifae: int
    ifpsp: int
    ifprt: int
    ifplw: int
    ilogd: int
    rlogd: float
    emin: float
    emax: float
    nnt: int
    thresh: float
    tol: float
    damp: float
    maxit: int
    title: str


@dataclass
class Shell:
    nnlz: int       # n, l, 0 as three digits (320 = 3d)
    wwnl: float     # electrons (negative: orbital skipped)
    ee: float       # starting guess for the eigenvalue (Ry)


@dataclass
class AEInput:
    controls: Controls
    z: float
    xion: float
    exfact: float
    rmax: float
    aasf: float
    bbsf: float
    irel: int
    shells: List[Shell]


@dataclass
class Beta:
    l: int
    keyee: int
    eeread: float
    iptype: int


@dataclass
class QFix:
    ibeta: int
    rfix: float
    qfix: float


@dataclass
class PSInput:
    controls: Controls
    ncores: int
    nvales: int
    nang: int
    besrmax: float
    besemin: float
    besemax: float
    besde: float
    keyps: int
    ifpcor: int
    rinner1: float
    rpcor: float
    rinner: List[float]      # 2*nang-1 values, written only when rinner1 == 0
    rcloc: float
    rc: List[float]          # nang values
    betas: List[Beta]
    npf: int                 # npf, ptryc written only when some beta has iptype 2
    ptryc: float
    lloc: int
    keyee_loc: int
    eloc: float
    iploctype: int
    ifqopt: int
    nqf: int
    qtryc: float
    fixes: List[QFix]


@dataclass
class TestInput:
    controls: Controls
    ncores: int
    nvales: int
    nang: int
    besrmax: float
    besemin: float
    besemax: float
    besde: float


@dataclass
class TestConfig:
    """A transferability test: the generation configuration with other occupations,
    or (library cases) a test configuration with its own shell list."""
    name: str
    xion: float
    occupations: Dict[int, float]                # nnlz -> electrons; others unchanged
    test_input: Optional[TestInput] = None       # None: derived from the generation input
    ae_input: Optional[AEInput] = None           # its own all-electron input (other shells); wins over occupations


@dataclass
class CaseInputs:
    symbol: str
    generation_name: str
    ae: AEInput
    ps: PSInput
    tests: List[TestConfig] = field(default_factory=list)


class AdatError(ValueError):
    def __init__(self, line_number: int, message: str):
        super().__init__(f"line {line_number}: {message}")
        self.line_number = line_number


def _line(key: str, values: Sequence) -> List[str]:
    fmt, label = LINE_FORMATS[key]
    try:
        return format_line(fmt, list(values), label)
    except ValueError as exc:
        raise ValueError(f"{label}: {exc}") from None


def _control_lines(c: Controls) -> List[str]:
    return (_line("controls", [c.ifae, c.ifpsp, c.ifprt, c.ifplw, c.ilogd])
            + _line("logd", [c.rlogd, c.emin, c.emax, c.nnt])
            + _line("convergence", [c.thresh, c.tol, c.damp, c.maxit])
            + _line("title", [c.title]))


def format_adat(data: Union[AEInput, PSInput, TestInput]) -> str:
    lines = _control_lines(data.controls)
    if isinstance(data, AEInput):
        lines += _line("atom", [data.z, data.xion, data.exfact])
        lines += _line("mesh", [data.rmax, data.aasf, data.bbsf])
        lines += _line("nstates", [len(data.shells), data.irel])
        for s in data.shells:
            lines += _line("shell", [s.nnlz, s.wwnl, s.ee])
    else:
        lines += _line("partition", [data.ncores, data.nvales, data.nang])
        lines += _line("bessel", [data.besrmax, data.besemin, data.besemax, data.besde])
    if isinstance(data, PSInput):
        lines += _line("keyps", [data.keyps, data.ifpcor, data.rinner1, data.rpcor])
        if data.rinner1 == 0.0:
            lines += _line("rinner", data.rinner)
        lines += _line("nbeta", [len(data.betas), data.rcloc])
        lines += _line("rc", data.rc)
        for b in data.betas:
            lines += _line("beta", [b.l, b.keyee, b.eeread, b.iptype])
        if any(b.iptype == 2 for b in data.betas):
            lines += _line("npf", [data.npf, data.ptryc])
        lines += _line("local", [data.lloc, data.keyee_loc, data.eloc, data.iploctype])
        lines += _line("qfunc", [data.ifqopt, data.nqf, data.qtryc, len(data.fixes)])
        for f in data.fixes:
            lines += _line("qfix", [f.ibeta, f.rfix, f.qfix])
    return "\n".join(lines) + "\n"


class _Reader:
    def __init__(self, text: str):
        self.lines = text.replace("\r\n", "\n").split("\n")
        if self.lines and self.lines[-1] == "":
            self.lines.pop()            # the newline that ends the last record
        self.index = 0

    def read(self, key: str, count: Optional[int] = None) -> list:
        fmt, label = LINE_FORMATS[key]
        per_line = len(parse_format(fmt))
        wanted = per_line if count is None else count
        values: list = []
        while len(values) < wanted or not values and wanted == 0:
            if self.index >= len(self.lines):
                raise AdatError(self.index + 1, f"end of file while reading {label}")
            line_number, line = self.index + 1, self.lines[self.index]
            self.index += 1
            try:
                values += read_record(line, fmt, min(per_line, wanted - len(values)))
            except ValueError as exc:
                raise AdatError(line_number, f"{label}: {exc}") from None
            if wanted == 0:
                break
        return values


def parse_adat(text: str) -> Union[AEInput, PSInput, TestInput]:
    """An .adat file read exactly as rwsubs.f reads it (extra lines are ignored, as there)."""
    r = _Reader(text)
    ifae, ifpsp, ifprt, ifplw, ilogd = r.read("controls", 5)
    rlogd, emin, emax, nnt = r.read("logd")
    thresh, tol, damp, maxit = r.read("convergence")
    (title,) = r.read("title")
    controls = Controls(ifae, ifpsp, ifprt, ifplw, ilogd, rlogd, emin, emax, nnt,
                        thresh, tol, damp, maxit, title.rstrip())
    if ifae == 1:
        if ifpsp != 0:
            raise AdatError(1, "an all-electron run and a pseudopotential run in one file is not "
                               "supported here; use separate files as the Work/ library does")
        z, xion, exfact = r.read("atom")
        rmax, aasf, bbsf = r.read("mesh")
        ncspvs, irel = r.read("nstates")
        shells = [Shell(*r.read("shell")) for _ in range(max(ncspvs, 0))]
        return AEInput(controls, z, xion, exfact, rmax, aasf, bbsf, irel, shells)
    if ifae != 0 or ifpsp not in (1, 2):
        raise AdatError(1, f"ifae = {ifae}, ifpsp = {ifpsp}: not an all-electron, "
                           "generation or test input")
    ncores, nvales, nang = r.read("partition")
    if not 1 <= nang <= 4:
        raise AdatError(r.index, f"nang = {nang}; runatom accepts 1 to 4")
    bessel = r.read("bessel")
    if ifpsp == 1:
        return TestInput(controls, ncores, nvales, nang, *bessel)
    keyps, ifpcor, rinner1, rpcor = r.read("keyps")
    rinner = r.read("rinner", 2 * nang - 1) if rinner1 == 0.0 else []
    nbeta, rcloc = r.read("nbeta")
    rc = r.read("rc", nang)
    betas = [Beta(*r.read("beta")) for _ in range(max(nbeta, 0))]
    npf, ptryc = r.read("npf") if any(b.iptype == 2 for b in betas) else (0, 0.0)
    lloc, keyee_loc, eloc, iploctype = r.read("local")
    ifqopt, nqf, qtryc, nfix = r.read("qfunc")
    fixes = [QFix(*r.read("qfix")) for _ in range(max(nfix, 0))]
    return PSInput(controls, ncores, nvales, nang, *bessel, keyps, ifpcor, rinner1, rpcor,
                   rinner, rcloc, rc, betas, npf, ptryc, lloc, keyee_loc, eloc, iploctype,
                   ifqopt, nqf, qtryc, fixes)


def file_stems(case: CaseInputs) -> Dict[str, object]:
    """File stems in the library's naming (Work/000-Ex Makefiles)."""
    s = case.symbol
    return {"ae": f"{s}_ae_{case.generation_name}", "ps": f"{s}_ps",
            "tests": [(f"{s}_ae_{t.name}", f"{s}_test_{t.name}") for t in case.tests]}


def test_ae_input(case: CaseInputs, test: TestConfig) -> AEInput:
    if test.ae_input is not None:
        return test.ae_input
    shells = [replace(s, wwnl=test.occupations.get(s.nnlz, s.wwnl)) for s in case.ae.shells]
    return replace(case.ae, xion=test.xion, shells=shells)


def test_input(case: CaseInputs, test: TestConfig) -> TestInput:
    if test.test_input is not None:
        return test.test_input
    ps = case.ps
    controls = replace(ps.controls, ifae=0, ifpsp=1, ifplw=0)
    return TestInput(controls, ps.ncores, ps.nvales, ps.nang,
                     ps.besrmax, ps.besemin, ps.besemax, ps.besde)


def case_files(case: CaseInputs) -> Dict[str, str]:
    stems = file_stems(case)
    files = {f"{stems['ae']}.adat": format_adat(case.ae), f"{stems['ps']}.adat": format_adat(case.ps)}
    for test, (ae_stem, test_stem) in zip(case.tests, stems["tests"]):
        files[f"{ae_stem}.adat"] = format_adat(test_ae_input(case, test))
        files[f"{test_stem}.adat"] = format_adat(test_input(case, test))
    return files


def write_case(case: CaseInputs, folder) -> List[Path]:
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    for name, text in case_files(case).items():
        path = folder / name
        path.write_text(text, encoding="ascii", newline="\n")
        written.append(path)
    return written


def _makefile_vars(text: str) -> Dict[str, str]:
    return {m.group(1): m.group(2).strip()
            for m in re.finditer(r"^([A-Z][A-Z0-9_]*)\s*=\s*(.*?)\s*$", text, re.M)}


def read_case_folder(folder) -> CaseInputs:
    """A Work/ library case (Makefile + .adat files) as CaseInputs."""
    folder = Path(folder)
    variables = _makefile_vars((folder / "Makefile").read_text(encoding="latin-1"))
    symbol, generation = variables["ATOM"], variables["CFG"]
    # test configurations: CFT, CFT1.. (most cases) or CF1, CF2 (000-Ex-ca-sp-tests)
    test_keys = sorted((k for k in variables if re.fullmatch(r"CF(T\d*|\d+)", k)),
                       key=lambda k: int(re.sub(r"\D", "", k) or 0))

    def load(name: str):
        try:
            return parse_adat((folder / name).read_text(encoding="latin-1"))
        except AdatError as exc:
            raise AdatError(exc.line_number, f"{name}: {exc}") from None

    ae, ps = load(f"{symbol}_ae_{generation}.adat"), load(f"{symbol}_ps.adat")
    reference = {s.nnlz: s.wwnl for s in ae.shells}
    tests = []
    for key in test_keys:
        name = variables[key]
        if not (folder / f"{symbol}_test_{name}.adat").exists():
            continue
        test_ae = load(f"{symbol}_ae_{name}.adat")
        test_file = load(f"{symbol}_test_{name}.adat")
        if [s.nnlz for s in test_ae.shells] != [s.nnlz for s in ae.shells]:
            # 12 library cases test configurations with another shell list: keep their own input
            tests.append(TestConfig(name, test_ae.xion, {}, test_file, test_ae))
            continue
        changed = {s.nnlz: s.wwnl for s in test_ae.shells if s.wwnl != reference[s.nnlz]}
        tests.append(TestConfig(name, test_ae.xion, changed, test_file))
    return CaseInputs(symbol, generation, ae, ps, tests)


def makefile_text(case: CaseInputs) -> str:
    """A Makefile in the Work/000-Ex template's style, for running the case by hand."""
    tests = case.tests
    head = [f"ATOM= {case.symbol}", f"CFG= {case.generation_name}"]
    head += [f"CFT{k}= {t.name}" for k, t in enumerate(tests, 1)]
    names = ["", "ROOT= ../../..", "PROG= ${ROOT}/Bin/runatom.x", "POT=  ${ROOT}/Pot", "",
             "AEG= ${ATOM}_ae_${CFG}", "PSP= ${ATOM}_ps"]
    names += [f"AET{k}= ${{ATOM}}_ae_${{CFT{k}}}" for k in range(1, len(tests) + 1)]
    names += [f"TST{k}= ${{ATOM}}_test_${{CFT{k}}}" for k in range(1, len(tests) + 1)]
    rules = ["", "default: pseudo", "",
             "${AEG}.ae: ${AEG}.adat",
             "\t${PROG} ${AEG}.adat ${AEG}.out ${AEG}.ae ${AEG}.atwf ${AEG}.logd dummy", "",
             "pseudo: ${PSP}.uspp", "${PSP}.uspp: ${PSP}.adat ${AEG}.ae",
             "\t${PROG} ${PSP}.adat ${PSP}.out ${AEG}.ae ${PSP}.atwf ${PSP}.logd ${PSP}.uspp", ""]
    rules.append("test: " + " ".join(f"${{TST{k}}}.out" for k in range(1, len(tests) + 1)))
    for k in range(1, len(tests) + 1):
        rules += ["", f"${{AET{k}}}.ae: ${{AET{k}}}.adat",
                  f"\t${{PROG}} ${{AET{k}}}.adat ${{AET{k}}}.out ${{AET{k}}}.ae ${{AET{k}}}.atwf ${{AET{k}}}.logd dummy",
                  f"${{TST{k}}}.out: ${{TST{k}}}.adat ${{AET{k}}}.ae ${{PSP}}.uspp",
                  f"\t${{PROG}} ${{TST{k}}}.adat ${{TST{k}}}.out ${{AET{k}}}.ae ${{TST{k}}}.atwf ${{TST{k}}}.logd ${{PSP}}.uspp"]
    rules += ["", "clean:", "\t- rm *.out *.ae *.atwf *.logd *.uspp", ""]
    header = ["# Written by gqteaWinToga (Vanderbilt USPP Generator) for USPP 7.3.6.",
              "# make: generate the pseudopotential; make test: run the transferability tests."]
    return "\n".join(header + head + names + rules)


# Capacities compiled into runatom (runatom.f parameters; initsubs.f for the mesh)
MAX_MESH = 1000
MAX_SHELLS = 26
MAX_BETAS = 10
MAX_VALENCE = 5
MAX_NANG = 4
MAX_BETAS_PER_L = 4
MAX_POLY_TERMS = 20
TITLE_WIDTH = 20
OCCUPATION_TOLERANCE = 1e-3
XC_FUNCTIONALS = {0.0: "Ceperley-Alder (LDA)", -1.0: "Wigner (LDA)", -2.0: "Hedin-Lundqvist (LDA)",
                  -3.0: "Gunnarsson-Lundqvist (LDA)", 1.0: "CA + Becke88 + LYP (GGA)",
                  2.0: "CA + Becke88 (GGA)", 3.0: "CA + Becke88 + Perdew86 (GGA)",
                  4.0: "Perdew-Wang 91 (GGA)", 5.0: "PBE (GGA)"}

FORM_FIELDS = {
    "symbol": ("Element symbol (file names)", "s", False, "e.g. ti gives ti_ae_<configuration>.adat and ti_ps.adat"),
    "generation": ("Generation configuration", "s", False, "a short name such as s2d1"),
    "title": ("Title (at most 20 characters)", "s", False, ""),
    "ae_print": ("ifprt  ifplw  ilogd", "iii", False, "printout level, wavefunction plot file (0/1), number of l for log derivatives"),
    "ae_logd": ("rlogd  emin  emax  nnt", "fffi", False, "log-derivative radius (bohr) and energy mesh (Ry)"),
    "ae_conv": ("thresh  tol  damp  maxit", "fffi", False, "0 = program default (1e-6, 1e-5, 0.5, 250)"),
    "ae_atom": ("z  xion  exfact", "fff", False, "exfact: 0 CA, -1 Wigner, -2 HL, -3 GL, 1 BLYP, 2 B88, 3 BP86, 4 PW91, 5 PBE"),
    "ae_mesh": ("rmax  aasf  bbsf  irel", "fffi", False, "radial mesh (0 = 80, 6, 40); irel 0 non-relativistic, 2 scalar-relativistic"),
    "ae_shells": ("Shells: nnlz  wwnl  ee", "iff", True, "one shell per line, core first, valence last; e.g. 300 2.0 -3.0"),
    "ps_print": ("ifprt  ifplw  ilogd", "iii", False, ""),
    "ps_logd": ("rlogd  emin  emax  nnt", "fffi", False, "rinner must not exceed rlogd"),
    "ps_conv": ("thresh  tol  damp  maxit", "fffi", False, ""),
    "ps_partition": ("ncores  nvales  nang", "iii", False, "ncores + nvales = number of shells; nang = lmax + 1 of the projectors"),
    "ps_bessel": ("besrmax  besemin  besemax  besde", "ffff", False, "ghost check in a Bessel basis: box radius (bohr), cutoffs (Ry)"),
    "ps_keyps": ("keyps  ifpcor  rinner  rpcor", "iiff", False, "keyps 3; ifpcor 1 = nonlinear core correction; rinner 0 = per-L values below"),
    "ps_rinner": ("rinner per L (if rinner = 0)", "f*", False, "2*nang - 1 values; ignored when rinner is not 0"),
    "ps_rcloc": ("rcloc", "f", False, "radius for the local potential"),
    "ps_rc": ("rc per l", "f*", False, "nang values (s, p, d, f)"),
    "ps_betas": ("Projectors: l  keyee  eeread  iptype", "iifi", True, "one per line, l non-decreasing; keyee > 0 takes that valence eigenvalue"),
    "ps_npf": ("npf  ptryc", "if", False, "written only when a projector has iptype 2"),
    "ps_local": ("lloc  keyee  eloc  iploctype", "iifi", False, "lloc -1 = polynomial inside rcloc"),
    "ps_qfunc": ("ifqopt  nqf  qtryc", "iif", False, "GGA functionals need ifqopt 3"),
    "ps_fixes": ("Q constraints: ibfix  rfix  qfix", "iff", True, "optional; ifqopt 2 or 3 only"),
    "tests": ("Tests: name  xion  nnlz=occupation ...", "t", True, "e.g. s2d2 0.0 320=2  (unlisted shells keep the generation occupation)"),
}


def parse_line_values(text: str, kinds: str, name: str) -> list:
    """Space- or comma-separated values typed by a user, typed by `kinds`."""
    tokens = [t for t in re.split(r"[\s,]+", (text or "").strip()) if t]
    if kinds.endswith("*"):
        if not tokens:
            raise ValueError(f"{name}: enter at least one value")
        kinds = kinds[0] * len(tokens)
    if len(tokens) != len(kinds):
        raise ValueError(f"{name}: expected {len(kinds)} values, got {len(tokens)}")
    values = []
    for token, kind in zip(tokens, kinds):
        try:
            values.append(int(token) if kind == "i" else float(token.lower().replace("d", "e")))
        except ValueError:
            what = "an integer" if kind == "i" else "a number"
            raise ValueError(f"{name}: {token!r} is not {what}") from None
    return values


def parse_rows(text: str, kinds: str, name: str) -> List[list]:
    rows = [line for line in (text or "").splitlines() if line.strip()]
    return [parse_line_values(line, kinds, f"{name}, row {k}") for k, line in enumerate(rows, 1)]


def _number(value: float) -> str:
    text = repr(float(value))
    return text[:-2] if text.endswith(".0") else text


def _text(values) -> str:
    return "  ".join(str(v) if isinstance(v, int) else _number(v) for v in values)


def _controls_fields(prefix: str, c: Controls) -> Dict[str, str]:
    return {f"{prefix}_print": _text([c.ifprt, c.ifplw, c.ilogd]),
            f"{prefix}_logd": _text([c.rlogd, c.emin, c.emax, c.nnt]),
            f"{prefix}_conv": _text([c.thresh, c.tol, c.damp, c.maxit])}


def ae_fields(ae: AEInput) -> Dict[str, str]:
    fields = _controls_fields("ae", ae.controls)
    fields.update(ae_atom=_text([ae.z, ae.xion, ae.exfact]),
                  ae_mesh=_text([ae.rmax, ae.aasf, ae.bbsf, ae.irel]),
                  ae_shells="\n".join(_text([s.nnlz, s.wwnl, s.ee]) for s in ae.shells))
    return fields


def ps_fields(ps: PSInput) -> Dict[str, str]:
    fields = _controls_fields("ps", ps.controls)
    fields.update(
        ps_partition=_text([ps.ncores, ps.nvales, ps.nang]),
        ps_bessel=_text([ps.besrmax, ps.besemin, ps.besemax, ps.besde]),
        ps_keyps=_text([ps.keyps, ps.ifpcor, ps.rinner1, ps.rpcor]),
        ps_rinner=_text(ps.rinner), ps_rcloc=_text([ps.rcloc]), ps_rc=_text(ps.rc),
        ps_betas="\n".join(_text([b.l, b.keyee, b.eeread, b.iptype]) for b in ps.betas),
        ps_npf=_text([ps.npf, ps.ptryc]) if any(b.iptype == 2 for b in ps.betas) else "",
        ps_local=_text([ps.lloc, ps.keyee_loc, ps.eloc, ps.iploctype]),
        ps_qfunc=_text([ps.ifqopt, ps.nqf, ps.qtryc]),
        ps_fixes="\n".join(_text([f.ibeta, f.rfix, f.qfix]) for f in ps.fixes))
    return fields


def case_to_fields(case: CaseInputs) -> Dict[str, str]:
    fields = {"symbol": case.symbol, "generation": case.generation_name, "title": case.ae.controls.title}
    fields.update(ae_fields(case.ae))
    fields.update(ps_fields(case.ps))
    fields["tests"] = "\n".join(          # only occupation changes fit this field
        "  ".join([t.name, _number(t.xion)] + [f"{n}={_number(w)}" for n, w in t.occupations.items()])
        for t in case.tests if t.ae_input is None)
    return fields


def unrepresentable_tests(case: CaseInputs) -> List[str]:
    """Tests whose all-electron input has its own shell list: the Tests field cannot hold them."""
    return [t.name for t in case.tests if t.ae_input is not None]


def _values(fields: Dict[str, str], key: str):
    label, kinds, multiline, _ = FORM_FIELDS[key]
    text = fields.get(key, "")
    return parse_rows(text, kinds, label) if multiline else parse_line_values(text, kinds, label)


def _controls_from(fields: Dict[str, str], prefix: str, ifae: int, ifpsp: int, title: str) -> Controls:
    ifprt, ifplw, ilogd = _values(fields, f"{prefix}_print")
    rlogd, emin, emax, nnt = _values(fields, f"{prefix}_logd")
    thresh, tol, damp, maxit = _values(fields, f"{prefix}_conv")
    return Controls(ifae, ifpsp, ifprt, ifplw, ilogd, rlogd, emin, emax, nnt, thresh, tol, damp, maxit, title)


def _parse_tests(text: str) -> List[TestConfig]:
    tests = []
    for k, line in enumerate((l for l in (text or "").splitlines() if l.strip()), 1):
        tokens = line.split()
        if len(tokens) < 2:
            raise ValueError(f"Tests, row {k}: write a name and xion, then nnlz=occupation pairs")
        (xion,) = parse_line_values(tokens[1], "f", f"Tests, row {k}, xion")
        occupations = {}
        for pair in tokens[2:]:
            nnlz, _, value = pair.partition("=")
            if not value:
                raise ValueError(f"Tests, row {k}: {pair!r} is not nnlz=occupation")
            (n,) = parse_line_values(nnlz, "i", f"Tests, row {k}, {pair}")
            (w,) = parse_line_values(value, "f", f"Tests, row {k}, {pair}")
            occupations[n] = w
        tests.append(TestConfig(tokens[0], xion, occupations))
    return tests


def case_from_fields(fields: Dict[str, str]) -> CaseInputs:
    """The window's text fields as CaseInputs; ValueError names the field at fault."""
    title = (fields.get("title") or "").rstrip()      # as parse_adat keeps it (a20, trailing blanks)
    z, xion, exfact = _values(fields, "ae_atom")
    rmax, aasf, bbsf, irel = _values(fields, "ae_mesh")
    shells = [Shell(*row) for row in _values(fields, "ae_shells")]
    ae = AEInput(_controls_from(fields, "ae", 1, 0, title), z, xion, exfact, rmax, aasf, bbsf, irel, shells)
    ncores, nvales, nang = _values(fields, "ps_partition")
    keyps, ifpcor, rinner1, rpcor = _values(fields, "ps_keyps")
    betas = [Beta(*row) for row in _values(fields, "ps_betas")]
    npf, ptryc = _values(fields, "ps_npf") if any(b.iptype == 2 for b in betas) else (0, 0.0)
    lloc, keyee_loc, eloc, iploctype = _values(fields, "ps_local")
    ifqopt, nqf, qtryc = _values(fields, "ps_qfunc")
    ps = PSInput(_controls_from(fields, "ps", 0, 2, title), ncores, nvales, nang, *_values(fields, "ps_bessel"),
                 keyps, ifpcor, rinner1, rpcor,
                 _values(fields, "ps_rinner") if rinner1 == 0.0 else [],
                 _values(fields, "ps_rcloc")[0], _values(fields, "ps_rc"), betas, npf, ptryc,
                 lloc, keyee_loc, eloc, iploctype, ifqopt, nqf, qtryc,
                 [QFix(*row) for row in _values(fields, "ps_fixes")])
    return CaseInputs((fields.get("symbol") or "").strip(), (fields.get("generation") or "").strip(),
                      ae, ps, _parse_tests(fields.get("tests", "")))


def radial_mesh_size(z: float, rmax: float, aasf: float, bbsf: float) -> int:
    """Number of radial points runatom will use (initsubs.f, subroutine rinit)."""
    rmax, aasf, bbsf = rmax or 80.0, aasf or 6.0, bbsf or 40.0
    a = math.exp(-aasf) / z
    mesh = 2 + int(math.log(rmax / a + 1.0) / (1.0 / bbsf))
    return mesh + 1 if mesh % 2 == 0 else mesh


def check_ae(ae: AEInput) -> List[str]:
    """runatom's own all-electron checks (rwsubs.f readin, initsubs.f rinit) plus this tool's limits."""
    problems = []
    if len(ae.controls.title) > TITLE_WIDTH or not ae.controls.title.isascii():
        problems.append(f"The title must be at most {TITLE_WIDTH} ASCII characters.")
    if ae.z <= 0:
        problems.append("z must be positive.")
    if ae.exfact not in XC_FUNCTIONALS:
        problems.append(f"exfact = {ae.exfact:g} is not one of {sorted(XC_FUNCTIONALS)}.")
    if ae.irel not in (0, 2):
        problems.append("irel must be 0 or 2: runatom generates and tests pseudopotentials "
                        "only without the Dirac equation (irel = 1).")
    if not 1 <= len(ae.shells) <= MAX_SHELLS:
        problems.append(f"Give between 1 and {MAX_SHELLS} shells (runatom's limit).")
    for s in ae.shells:
        n, l = s.nnlz // 100, (s.nnlz // 10) % 10
        if not (100 <= s.nnlz <= 999 and s.nnlz % 10 == 0 and l < n and l <= 3):
            problems.append(f"Shell {s.nnlz}: write n, l and 0 as three digits (320 = 3d), with l < n and l <= 3.")
    electrons = sum(max(s.wwnl, 0.0) for s in ae.shells)
    if abs(ae.z - electrons - ae.xion) > OCCUPATION_TOLERANCE:
        problems.append(f"The occupations add up to {electrons:g} electrons, but z - xion = {ae.z - ae.xion:g}.")
    if ae.z > 0:
        mesh = radial_mesh_size(ae.z, ae.rmax, ae.aasf, ae.bbsf)
        if mesh > MAX_MESH:
            problems.append(f"The radial mesh would have {mesh} points; runatom allows {MAX_MESH}. Lower rmax or bbsf.")
    c = ae.controls
    if min(c.thresh, c.tol, c.damp) < 0 or c.maxit < 0:
        problems.append("thresh, tol, damp and maxit cannot be negative.")
    return problems


def check_ps(ps: PSInput, ae: AEInput) -> List[str]:
    """runatom's generation checks (rwsubs.f readin, lines 300-673) plus this tool's limits."""
    problems = []
    if ps.keyps != 3:
        problems.append("keyps must be 3 (the only scheme runatom implements).")
    if ps.ncores + ps.nvales != len(ae.shells):
        problems.append(f"ncores + nvales = {ps.ncores + ps.nvales}, but the all-electron input has {len(ae.shells)} shells.")
    if not 1 <= ps.nvales <= MAX_VALENCE:
        problems.append(f"nvales must be between 1 and {MAX_VALENCE}.")
    if not 1 <= ps.nang <= MAX_NANG:
        return problems + [f"nang must be between 1 and {MAX_NANG}."]
    n_l = 2 * ps.nang - 1
    rinner = ps.rinner if ps.rinner1 == 0.0 else [ps.rinner1] * n_l
    if ps.rinner1 == 0.0 and len(ps.rinner) != n_l:
        problems.append(f"rinner = 0 needs {n_l} per-L values.")
    if rinner and max(rinner) > ps.controls.rlogd:
        problems.append(f"rinner ({max(rinner):g}) must not exceed rlogd ({ps.controls.rlogd:g}).")
    if len(ps.rc) != ps.nang:
        problems.append(f"Give {ps.nang} rc values, one per l.")
    if not 1 <= len(ps.betas) <= MAX_BETAS:
        problems.append(f"Give between 1 and {MAX_BETAS} projectors.")
    ls = [b.l for b in ps.betas]
    if ls != sorted(ls):
        problems.append("List the projectors with l non-decreasing (s, then p, then d).")
    for l in sorted(set(ls)):
        if ls.count(l) > MAX_BETAS_PER_L:
            problems.append(f"At most {MAX_BETAS_PER_L} projectors per l (l = {l} has {ls.count(l)}).")
    if ls and max(ls) + 1 != ps.nang:
        problems.append(f"nang = {ps.nang}, but the highest projector l is {max(ls)}: nang must be lmax + 1.")
    for k, b in enumerate(ps.betas, 1):
        if b.keyee < 0:
            problems.append(f"Projector {k}: keyee < 0 (wavefunctions from another reference "
                            "configuration) is not supported by this tool.")
        elif b.keyee > ps.nvales:
            problems.append(f"Projector {k}: keyee = {b.keyee} exceeds nvales = {ps.nvales}.")
        if not 0 <= b.iptype <= 4:
            problems.append(f"Projector {k}: iptype must be 0 to 4.")
    smooth = any(b.iptype == 2 for b in ps.betas)
    if smooth and not 4 < ps.npf <= MAX_POLY_TERMS:
        problems.append(f"npf must be between 5 and {MAX_POLY_TERMS} when a projector has iptype 2.")
    if ps.keyee_loc > ps.nvales:
        problems.append(f"Local potential: keyee = {ps.keyee_loc} exceeds nvales = {ps.nvales}.")
    if ps.keyee_loc > 0 and ps.iploctype == 2 and not smooth:
        problems.append("Local potential: iploctype 2 needs npf and ptryc, i.e. a projector with iptype 2.")
    if ps.ifqopt > 3:
        problems.append("ifqopt must be -1 (no pseudization), 0, 1, 2 or 3.")
    if 1 <= ps.ifqopt <= 3 and not 3 < ps.nqf <= MAX_POLY_TERMS:
        problems.append(f"ifqopt {ps.ifqopt} needs 3 < nqf <= {MAX_POLY_TERMS}.")
    if ps.fixes and ps.ifqopt < 2:
        problems.append("Q-function constraints need ifqopt 2 or 3.")
    if len(ps.fixes) > len(ps.betas):
        problems.append("There are more Q-function constraints than projectors.")
    previous = -1
    for k, f in enumerate(ps.fixes, 1):
        if f.ibeta <= previous or f.ibeta > len(ps.betas):
            problems.append(f"Constraint {k}: ibfix must increase and be at most nbeta.")
        previous = f.ibeta
        if rinner and f.rfix > min(rinner):
            problems.append(f"Constraint {k}: rfix must not exceed rinner.")
    if ae.exfact >= 1 and ps.ifqopt != 3:
        problems.append("GGA functionals (exfact >= 1) need ifqopt = 3: runatom refuses -1, 0, 1 and 2.")
    return problems


def check_case(case: CaseInputs) -> List[str]:
    """Every problem runatom would stop on, found before it runs; [] when ready."""
    problems = []
    names = [case.symbol, case.generation_name] + [t.name for t in case.tests]
    for name in names:
        if not SAFE_NAME.match(name or ""):
            problems.append(f"Name {name!r}: use only letters, digits, '.', '_' and '-'.")
    test_names = [t.name for t in case.tests]
    if len(set(test_names)) != len(test_names) or case.generation_name in test_names:
        problems.append("Test names must be different from each other and from the generation configuration.")
    if problems:
        return problems
    stems = file_stems(case)
    longest = max([stems["ae"], stems["ps"]] + [s for pair in stems["tests"] for s in pair], key=len)
    if len(longest) + len(".atwf") > MAX_FORTRAN_NAME:
        return [f"The file name {longest}.atwf is longer than {MAX_FORTRAN_NAME} characters, which "
                "runatom cannot read; shorten the symbol or the configuration names."]
    try:      # check the values runatom will read, after the column-exact round trip
        ae = parse_adat(format_adat(case.ae))
        ps = parse_adat(format_adat(case.ps))
    except ValueError as exc:
        return [str(exc)]
    problems += [f"All-electron: {p}" for p in check_ae(ae)]
    problems += [f"Pseudopotential: {p}" for p in check_ps(ps, ae)]
    known = {s.nnlz for s in ae.shells}
    for test in case.tests:
        unknown = sorted(set(test.occupations) - known) if test.ae_input is None else []
        if unknown:
            problems.append(f"Test {test.name}: shells {unknown} are not in the all-electron shell list.")
            continue
        try:
            test_ae = parse_adat(format_adat(test_ae_input(case, test)))
        except ValueError as exc:
            problems.append(f"Test {test.name}: {exc}")
            continue
        problems += [f"Test {test.name}: {p}" for p in check_ae(test_ae)]
    return problems
