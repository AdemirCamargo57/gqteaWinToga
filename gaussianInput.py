"""
Gaussian 16 input builder: the job model, validation, route, input text and reopening (no Toga).

``GaussianJobSpec`` (job-wide fields + one ``GaussianStepSpec`` per linked step) ->
``GaussianInputBuilder.prepare`` (every problem in one ValueError, one '- Step n: ...' line each) ->
``GaussianJob`` -> ``input_text`` / ``script_text`` / ``write_files``. ``read_input_text`` is the
inverse used by Open input. The Fukui Job Type is delegated to ``fukuiInput`` unchanged. Files are
ASCII with LF endings and relative checkpoint names: the job runs under Linux from the input's folder.
"""

import re
from dataclasses import dataclass, field
from typing import List, Sequence, Tuple

import fukuiInput as fi
import gaussianKeywords as gk

PARENTHESIS_KEYWORDS = ("IOP",)        # written as IOp(9/40=4), the form every Gaussian manual uses


# ----------------------------------------------------------------------
# Route keywords
# ----------------------------------------------------------------------
@dataclass
class RouteKeyword:
    """One route keyword and its options, e.g. Opt=(TS,CalcFC) -> ('Opt', ['TS', 'CalcFC'])."""
    name: str
    options: List[str] = field(default_factory=list)

    @property
    def key(self) -> str:
        return self.name.upper()


def option_key(option: str) -> str:
    """'MaxCycles=50' -> 'MAXCYCLES': two options with the same key set the same thing."""
    return option.split("=", 1)[0].strip().upper()


def _outer_parentheses(text: str) -> str:
    """'(a,b)' -> 'a,b' when the first '(' closes at the very end; anything else unchanged."""
    if not (text.startswith("(") and text.endswith(")")):
        return text
    depth = 0
    for position, ch in enumerate(text):
        depth += (ch == "(") - (ch == ")")
        if depth == 0 and position < len(text) - 1:
            return text
    return text[1:-1]


def _split_options(text: str) -> List[str]:
    """Split on commas and blanks outside parentheses."""
    options, current, depth = [], "", 0
    for ch in text:
        if depth == 0 and (ch == "," or ch.isspace()):
            if current:
                options.append(current)
            current = ""
            continue
        depth += (ch == "(") - (ch == ")")
        current += ch
    if current:
        options.append(current)
    return options


def keyword_text_name(token: str) -> str:
    """The keyword part of a route token as written: 'Opt=(TS)' -> 'Opt', 'IOp(9/40=4)' -> 'IOp'."""
    return re.split(r"[=(]", token, maxsplit=1)[0]


def parse_keyword(token: str) -> RouteKeyword:
    name = keyword_text_name(token)
    rest = token[len(name):]
    if rest.startswith("="):
        rest = rest[1:]
    return RouteKeyword(name, _split_options(_outer_parentheses(rest)))


def format_keyword(keyword: RouteKeyword) -> str:
    options = keyword.options
    if not options:
        return keyword.name
    if keyword.key in PARENTHESIS_KEYWORDS:
        return f"{keyword.name}({','.join(options)})"
    if len(options) == 1 and "=" not in options[0]:
        return f"{keyword.name}={options[0]}"
    return f"{keyword.name}=({','.join(options)})"


def merge_keywords(panel: Sequence[RouteKeyword], extra_text: str) -> Tuple[List[RouteKeyword], List[str]]:
    """The panels' keywords plus Additional Keywords (spec 5.1).

    A keyword already written by a panel takes the extra options; an option repeated with the same
    key is kept once, and a valued option (X=...) given another value is an error naming both places.
    """
    merged = [RouteKeyword(k.name, list(k.options)) for k in panel]
    panel_keys = {k.key for k in panel}
    errors: List[str] = []
    for token in fi.route_tokens(extra_text):
        if "/" in keyword_text_name(token):
            errors.append(f"Additional Keywords: '{token}' is a method/basis set; choose them in the "
                          "Method panel.")
            continue
        extra = parse_keyword(token)
        target = next((k for k in merged if k.key == extra.key), None)
        if target is None:
            merged.append(extra)
            continue
        for option in extra.options:
            same = [o for o in target.options if option_key(o) == option_key(option)]
            if not same:
                target.options.append(option)
            elif "=" in option and same[0].upper() != option.upper():
                owner = (f"the {gk.PANEL_OF_KEYWORD.get(extra.key, 'Job Type')} panel"
                         if extra.key in panel_keys else "an earlier Additional Keyword")
                errors.append(f"Additional Keywords: {token} sets {option}, but {owner} sets {same[0]}.")
    return merged, errors
