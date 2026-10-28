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
