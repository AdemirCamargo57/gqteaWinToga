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
from typing import List, Optional, Sequence, Tuple

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
    lines[0] = lines[0].ljust(max(LABEL_COLUMN, len(lines[0]) + 2)) + f"{label} ({fmt})"
    return lines
