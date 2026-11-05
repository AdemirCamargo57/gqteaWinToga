"""
The final optimized geometry of a Gaussian output (.log/.out), for the Gaussian 16 Input Builder (no Toga).

A Gaussian output is a sequence of jobs: each --Link1-- step of the input and each job Gaussian chains
itself (the Freq part of Opt Freq). A job ends with "Normal termination"; the next one starts with
"Initial command:" or "Link1:  Proceeding to internal job step number". A job optimizes when its route
has an Opt keyword. The route decides because the Freq job of Opt Freq also prints
"Optimization completed" (verified on G16 C.01). The geometry imported is the last orientation block
before the last "Optimization completed" of the last optimizing job that has one (Input orientation
preferred: it keeps the input's frame). An optimization that never completed is refused, never imported.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from fukuiInput import element_symbol

GAUSSIAN_BANNERS = ("Entering Gaussian System", "Gaussian, Inc.")
NORMAL_TERMINATION = "Normal termination of Gaussian"
ERROR_TERMINATION = "Error termination"
COMPLETED = "Optimization completed"
STOPPED = "Optimization stopped"
ORIENTATIONS = ("Input orientation", "Standard orientation", "Z-Matrix orientation")   # preference order
CHARGE_RE = re.compile(r"^\s*Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
PRINT_LEVEL_RE = re.compile(r"^#[PNTpnt]?(?=\s|$)")
OPT_RE = re.compile(r"(?i)(?:^|\s)opt(?=$|[\s=(])")

Atom = Tuple[str, float, float, float]


@dataclass
class OptimizedGeometry:
    atoms: List[Atom]                 # element, x, y, z in Angstrom, in the file's atom order
    charge: Optional[int]
    multiplicity: Optional[int]
    job_number: int                   # 1-based job that holds the completed optimization
    job_count: int
    orientation: str                  # one of ORIENTATIONS
    notes: List[str] = field(default_factory=list)


@dataclass
class _Job:
    number: int
    lines: List[str]
    route: str
    completions: List[int]            # indexes (in lines) of "Optimization completed"
    normal: bool                      # ended with Normal termination

    @property
    def optimizes(self) -> bool:
        return bool(OPT_RE.search(PRINT_LEVEL_RE.sub("", self.route.strip(), count=1)))

    def failure(self) -> str:
        """Why this job's optimization did not complete, as far as the file says."""
        for k, line in enumerate(self.lines):
            if STOPPED in line:
                detail = self.lines[k + 1].strip().lstrip("-").strip() if k + 1 < len(self.lines) else ""
                return f"the optimization stopped ({' '.join(detail.split())})" if detail else \
                    "the optimization stopped"
        if any(ERROR_TERMINATION in line for line in self.lines):
            return "Gaussian ended with an error termination during the optimization"
        return "the file ends before the optimization completed (the job is still running or was cut off)"


def _route(lines: List[str]) -> str:
    """The route echoed between dashed lines: from the first line starting with '#' to the next dashes."""
    start = next((k for k, line in enumerate(lines) if line.strip().startswith("#")), None)
    if start is None:
        return ""
    words = []
    for line in lines[start:]:
        if line.strip().startswith("---"):
            break
        words.append(line.strip())
    return " ".join(words)


def _jobs(lines: List[str]) -> List[_Job]:
    bounds, start = [], 0
    for k, line in enumerate(lines):
        if NORMAL_TERMINATION in line:
            bounds.append((start, k + 1, True))
            start = k + 1
    if any(line.strip() for line in lines[start:]):
        bounds.append((start, len(lines), False))
    jobs = []
    for number, (a, b, normal) in enumerate(bounds, start=1):
        part = lines[a:b]
        jobs.append(_Job(number, part, _route(part), [k for k, line in enumerate(part) if COMPLETED in line],
                         normal))
    return jobs


def _block(lines: List[str], start: int, name: str) -> List[Atom]:
    """The atoms of the orientation block whose title is at ``start``."""
    atoms, k = [], start + 5          # title, dashes, two heading lines, dashes
    while k < len(lines) and not lines[k].strip().startswith("---"):
        parts = lines[k].split()
        try:
            z, xyz = int(parts[1]), tuple(float(v) for v in parts[-3:])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{name}, line {k + 1}: cannot read the orientation row '{lines[k].strip()}'.") from exc
        if z <= 0:
            raise ValueError(f"{name}: atom {parts[0]} of the optimized geometry is a ghost or dummy atom "
                             f"(atomic number {z}); the builder takes real atoms only.")
        symbol = element_symbol(str(z))
        if symbol is None:
            raise ValueError(f"{name}: atomic number {z} is not an element.")
        atoms.append((symbol, *xyz))
        k += 1
    if not atoms:
        raise ValueError(f"{name}, line {start + 1}: the orientation block holds no atoms.")
    return atoms


def _converged_block(job: _Job, name: str) -> Tuple[List[Atom], str]:
    """The orientation block of the geometry step the last completion belongs to (Input orientation first).

    Each geometry step prints each kind of orientation at most once, so the step is everything after the
    previous block of the same kind as the last one.
    """
    lines, done = job.lines, job.completions[-1]
    titles: Dict[int, str] = {k: t for k, line in enumerate(lines[:done]) for t in ORIENTATIONS if f"{t}:" in line}
    if not titles:
        raise ValueError(f"{name}: job {job.number} completed its optimization but printed no geometry.")
    last = max(titles)
    earlier = [k for k in titles if k < last and titles[k] == titles[last]]
    floor = max(earlier) if earlier else -1
    step = {t: k for k, t in sorted(titles.items()) if k > floor}
    title = next(t for t in ORIENTATIONS if t in step)
    return _block(lines, step[title], name), title


def parse_optimized_geometry(text: str, name: str) -> OptimizedGeometry:
    """The last completed optimization of a Gaussian output; ValueError (naming ``name``) otherwise."""
    if not any(banner in text for banner in GAUSSIAN_BANNERS):
        raise ValueError(f"{name} is not a Gaussian output (no Gaussian banner); load a Gaussian .log/.out, "
                         "an input file or an .xyz.")
    jobs = _jobs(text.splitlines())
    optimizing = [j for j in jobs if j.optimizes]
    if not optimizing:
        raise ValueError(f"{name} contains no geometry optimization (no route has Opt); load its input file "
                         "or an .xyz instead.")
    done = [j for j in optimizing if j.completions]
    if not done:
        last = optimizing[-1]
        raise ValueError(f"{name}: no successfully optimized geometry was found: in job {last.number}, "
                         f"{last.failure()}. The last geometry is not imported because it is not optimized.")
    job = done[-1]
    atoms, orientation = _converged_block(job, name)
    charges = [m for m in (CHARGE_RE.match(line) for line in job.lines[:job.completions[-1]]) if m]
    charge, mult = (int(charges[-1].group(1)), int(charges[-1].group(2))) if charges else (None, None)
    notes = [f"job {j.number}'s optimization did not complete ({j.failure()}); the geometry of job {job.number}'s "
             "completed optimization was imported." for j in optimizing if j.number > job.number]
    if len(job.completions) > 1:
        notes.append(f"job {job.number} completed {len(job.completions)} optimizations (a relaxed scan); the "
                     "last one's geometry was imported.")
    if not job.normal:
        notes.append(f"job {job.number} did not end with Normal termination; its optimization had completed, "
                     "so that geometry was imported.")
    return OptimizedGeometry(atoms, charge, mult, job.number, len(jobs), orientation, notes)


def read_optimized_geometry(path) -> OptimizedGeometry:
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ValueError(f"'{path.name}' could not be opened: {exc.strerror or exc}") from exc
    return parse_optimized_geometry(text, path.name)
