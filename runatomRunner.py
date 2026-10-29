"""Running the Vanderbilt USPP program runatom as a child process.

No Toga here (a subprocess test pins it). runatom stops on errors with
`call exit(1)` (134 places) or `stop` (exit status 0) and takes its six file
names as character*40 command-line arguments. It therefore runs as a separate
process, with its working directory set to the job folder and short relative
names. A stage succeeds when the exit status is 0, its completion marker is on
stderr and its output files exist.
"""

import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional, Sequence

from runatomInput import MAX_FORTRAN_NAME, MAX_REFORM_NAME, CaseInputs, file_stems
from runatomOutput import OutReport, parse_out

DEFAULT_STAGE_TIMEOUT = 120.0     # s; the slowest measured stage takes 0.7 s
POLL_INTERVAL = 0.05              # s between cancel/timeout checks
ENV_VAR = "GQTEA_RUNATOM"
FROZEN_SUBDIR = "uspp"
VENDOR_DIR = Path(__file__).resolve().parent / "third_party" / "uspp-736"
VENDOR_EXE = VENDOR_DIR / "bin" / "runatom.exe"
COMPLETION_MARKERS = {"ae": "all electron calculation completed",
                      "pseudo": "writing the pseudopotential to file",
                      "test": "pseudopotential test converged"}
_SCF = re.compile(r"completed self-consistent cycle\s*(\d+)\s+delta\s*=\s*([-+0-9.EeDd]+)")
_NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)


class RunatomError(RuntimeError):
    def __init__(self, kind: str, stage_label: str, message: str):
        super().__init__(f"{stage_label}: {message}")
        self.kind, self.stage_label, self.message = kind, stage_label, message


class RunatomCancelled(RuntimeError):
    pass


def _bundle_dir() -> Optional[Path]:
    if getattr(sys, "frozen", False):
        return Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent)) / FROZEN_SUBDIR
    return None


def _find_program(explicit: Optional[str], env_var: str, exe_name: str, vendor_path: Path,
                  path_names: Sequence[str]) -> Optional[str]:
    """Explicit path, the environment variable, the frozen bundle, third_party, PATH; None if absent."""
    candidates = []
    if explicit and explicit.strip():
        candidates.append(explicit.strip())
    if os.environ.get(env_var):
        candidates.append(os.environ[env_var])
    bundle = _bundle_dir()
    if bundle is not None:
        candidates.append(str(bundle / exe_name))
    candidates.append(str(vendor_path))
    for name in path_names:
        found = shutil.which(name)
        if found:
            candidates.append(found)
    return next((c for c in candidates if os.path.isfile(c)), None)


def find_runatom(explicit: Optional[str] = None) -> Optional[str]:
    return _find_program(explicit, ENV_VAR, "runatom.exe", VENDOR_EXE, ("runatom.exe", "runatom.x"))


def library_dirs() -> List[Path]:
    """Where the USPP input library (Work/) may be, most preferred first.

    Frozen (PyInstaller onedir): <install dir>/third_party/uspp-736/Work, the
    same under the bundle's data folder (_internal/ since PyInstaller 6), then
    <bundle>/uspp/Work (the layout readme.md documents). From source: the
    repository copy. Never relative to the current working directory."""
    if not getattr(sys, "frozen", False):
        return [VENDOR_DIR / "Work"]
    install = Path(sys.executable).resolve().parent
    bundle = Path(getattr(sys, "_MEIPASS", install))
    relative = Path("third_party") / VENDOR_DIR.name / "Work"
    candidates = [install / relative, bundle / relative, bundle / FROZEN_SUBDIR / "Work"]
    return list(dict.fromkeys(candidates))


def library_dir() -> Path:
    """The first existing library folder, or the preferred one when none exists."""
    candidates = library_dirs()
    return next((c for c in candidates if c.is_dir()), candidates[0])


def runatom_command(exe: str) -> List[str]:
    return [exe]


@dataclass(frozen=True)
class Stage:
    kind: str           # "ae", "pseudo" or "test"
    label: str
    stem: str           # <stem>.adat, .out, .atwf, .logd
    ae_file: str        # .ae written (ae) or read (pseudo, test)
    uspp_file: str      # "dummy" for an all-electron stage

    def arguments(self) -> List[str]:
        return [f"{self.stem}.adat", f"{self.stem}.out", self.ae_file,
                f"{self.stem}.atwf", f"{self.stem}.logd", self.uspp_file]

    def required(self) -> List[str]:
        files = [f"{self.stem}.out"]
        files += [self.ae_file] if self.kind == "ae" else []
        files += [self.uspp_file] if self.kind == "pseudo" else []
        return files

    def products(self) -> List[str]:
        return self.required() + [f"{self.stem}.atwf", f"{self.stem}.logd"]


@dataclass
class StageResult:
    stage: Stage
    seconds: float
    report: OutReport
    returncode: int


def plan_stages(case: CaseInputs) -> List[Stage]:
    """The Work/000-Ex Makefile's chain: AE, generation, then per test AE + test."""
    stems = file_stems(case)
    ae_ref, uspp = f"{stems['ae']}.ae", f"{stems['ps']}.uspp"
    stages = [Stage("ae", f"All-electron ({case.generation_name})", stems["ae"], ae_ref, "dummy"),
              Stage("pseudo", "Pseudopotential generation", stems["ps"], ae_ref, uspp)]
    for test, (ae_stem, test_stem) in zip(case.tests, stems["tests"]):
        stages.append(Stage("ae", f"All-electron ({test.name})", ae_stem, f"{ae_stem}.ae", "dummy"))
        stages.append(Stage("test", f"Transferability test ({test.name})", test_stem, f"{ae_stem}.ae", uspp))
    return stages


def progress_text(info: dict) -> str:
    text = f"Stage {info['stage_index'] + 1}/{info['stage_count']}: {info['stage']}"
    if info.get("iteration") is not None:
        return f"{text} - SCF cycle {info['iteration']}, delta = {info['delta']:.2e}"
    return f"{text} - {info['message']}"


def _tail(stderr_text: str) -> str:
    lines = [l.strip() for l in stderr_text.splitlines() if l.strip() and "self-consistent cycle" not in l]
    return " | ".join(lines[-3:])


def classify(stage: Stage, returncode: int, stderr_text: str, out_text: str, workdir: Path) -> Optional[RunatomError]:
    def error(kind, message):
        return RunatomError(kind, stage.label, message)

    if "Fortran runtime error" in stderr_text:
        detail = " ".join(l.strip() for l in stderr_text.splitlines()
                          if l.startswith("At line") or "Fortran runtime error" in l)
        if "unit = 10" in stderr_text:
            return error("input_format", f"runatom could not read {stage.stem}.adat ({detail}); "
                                         "the file does not follow runatom's fixed columns.")
        return error("runtime_error", detail)
    if "negative densities" in stderr_text:
        return error("negative_density", "the pseudo charge density became negative while descreening; "
                                         "reduce rinner or add Q-function constraints (Doc/TUTORIAL).")
    lines = out_text.splitlines()
    for i, line in enumerate(lines):
        if "***error" in line or "*** error" in line:
            detail = " ".join(l.strip() for l in lines[i:i + 3] if l.strip())
            return error("not_converged" if "could not be converged" in detail else "program_error", detail)
    if returncode != 0:
        return error("runtime_error", f"runatom stopped with exit status {returncode}: {_tail(stderr_text)}")
    if COMPLETION_MARKERS[stage.kind] not in stderr_text:
        tail = _tail(stderr_text)
        return error("missing_output", "runatom stopped before the end of this stage" + (f": {tail}" if tail else "."))
    missing = [name for name in stage.required() if not (workdir / name).is_file()]
    if missing:
        return error("missing_output", "expected output not written: " + ", ".join(missing))
    return None


def _pump(stream, lines: "queue.Queue"):
    for raw in iter(stream.readline, b""):
        lines.put(raw.decode("latin-1"))
    lines.put(None)


def _stop(process: subprocess.Popen):
    if process.poll() is None:
        process.kill()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass


def _remove(workdir: Path, names: Sequence[str]):
    for name in names:
        try:
            (workdir / name).unlink()
        except FileNotFoundError:
            pass


def run_stage(command: Sequence[str], stage: Stage, workdir, progress: Optional[Callable[[dict], None]] = None,
              cancel: Optional[threading.Event] = None, timeout: float = DEFAULT_STAGE_TIMEOUT,
              stage_index: int = 0, stage_count: int = 1) -> StageResult:
    workdir = Path(workdir)
    arguments = stage.arguments()
    too_long = [a for a in arguments if len(a) > MAX_FORTRAN_NAME]
    if too_long:
        raise ValueError(f"file names longer than {MAX_FORTRAN_NAME} characters: {', '.join(too_long)}")

    def report(message: str, iteration: Optional[int] = None, delta: Optional[float] = None):
        if progress is None:
            return
        try:
            progress({"stage_index": stage_index, "stage_count": stage_count, "stage": stage.label,
                      "message": message, "iteration": iteration, "delta": delta})
        except Exception:      # reporting must never abort a run
            pass

    _remove(workdir, stage.products())
    report("starting")
    started = time.monotonic()
    process = subprocess.Popen(list(command) + arguments, cwd=str(workdir), stdin=subprocess.DEVNULL,
                               stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, creationflags=_NO_WINDOW)
    lines: "queue.Queue" = queue.Queue()
    pump = threading.Thread(target=_pump, args=(process.stderr, lines), daemon=True)
    pump.start()
    collected, finished_reading = [], False
    try:
        while not (finished_reading and process.poll() is not None):
            if cancel is not None and cancel.is_set():
                raise RunatomCancelled(f"{stage.label}: cancelled")
            if time.monotonic() - started > timeout:
                raise RunatomError("timeout", stage.label, f"no result after {timeout:.0f} s; the run was stopped.")
            try:
                line = lines.get(timeout=POLL_INTERVAL)
            except queue.Empty:
                continue
            if line is None:
                finished_reading = True
                continue
            collected.append(line)
            match = _SCF.search(line)
            if match:
                report(f"SCF cycle {match.group(1)}", int(match.group(1)),
                       float(match.group(2).upper().replace("D", "E")))
            elif line.strip():
                report(line.strip())
    except (RunatomCancelled, RunatomError):
        _stop(process)
        _remove(workdir, stage.products())
        raise
    finally:
        _stop(process)
        pump.join(timeout=5)       # the killed process closes the pipe; never close it under the reader
        process.stderr.close()
    seconds = time.monotonic() - started
    out_path = workdir / f"{stage.stem}.out"
    out_text = out_path.read_text(encoding="latin-1") if out_path.is_file() else ""
    problem = classify(stage, process.returncode, "".join(collected), out_text, workdir)
    if problem is not None:
        _remove(workdir, [n for n in stage.products() if n != f"{stage.stem}.out"])  # keep the report
        raise problem
    return StageResult(stage, seconds, parse_out(out_text), process.returncode)


def run_job(command: Sequence[str], stages: Sequence[Stage], workdir,
            progress: Optional[Callable[[dict], None]] = None, cancel: Optional[threading.Event] = None,
            timeout: float = DEFAULT_STAGE_TIMEOUT) -> List[StageResult]:
    """All stages in order; the first failure raises RunatomError (or RunatomCancelled)."""
    return [run_stage(command, stage, workdir, progress, cancel, timeout, index, len(stages))
            for index, stage in enumerate(stages)]


REFORM_ENV_VAR = "GQTEA_REFORM"
VENDOR_REFORM = VENDOR_DIR / "bin" / "reform.exe"
REFORM_DONE = "pseudopotential has been written out"
FORMATTED_LABEL = "Formatted file for CPMD"
REFORM_MISSING_NOTE = ("reform.exe was not found, so no formatted file for CPMD was written "
                       "(the binary .uspp is complete).")


def find_reform(explicit: Optional[str] = None) -> Optional[str]:
    return _find_program(explicit, REFORM_ENV_VAR, "reform.exe", VENDOR_REFORM, ("reform.exe", "reform.x"))


def reform_command(exe: str) -> List[str]:
    return [exe]


def formatted_name(case: CaseInputs) -> str:
    return f"{file_stems(case)['ps']}_formatted.uspp"


def cpmd_atom_line(name: str) -> str:
    """The &ATOMS species line that makes CPMD read the file as a formatted Vanderbilt potential."""
    return f"*{name} FORMATTED"


def write_formatted(command: Sequence[str], workdir, uspp_name: str, formatted: str,
                    timeout: float = DEFAULT_STAGE_TIMEOUT) -> Path:
    """The CPMD formatted copy of a binary .uspp, written by USPP's reform (itype 1), with LF endings."""
    workdir = Path(workdir)
    for name in (uspp_name, formatted):
        if len(name) > MAX_REFORM_NAME:
            raise ValueError(f"reform reads file names of at most {MAX_REFORM_NAME} characters: {name}")
    target = workdir / formatted
    _remove(workdir, [formatted])                  # reform opens its output with status='new'
    try:
        done = subprocess.run(list(command), cwd=str(workdir), input=f"1\n{uspp_name}\n{formatted}\n".encode("ascii"),
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout,
                              creationflags=_NO_WINDOW)
    except subprocess.TimeoutExpired:
        _remove(workdir, [formatted])
        raise RunatomError("timeout", FORMATTED_LABEL, f"no result after {timeout:.0f} s; the conversion was stopped.")
    stdout, stderr = done.stdout.decode("latin-1"), done.stderr.decode("latin-1")
    if done.returncode != 0 or REFORM_DONE not in stdout or not target.is_file():
        _remove(workdir, [formatted])
        detail = _tail(stderr)
        raise RunatomError("formatted", FORMATTED_LABEL,
                           f"reform could not convert {uspp_name}" + (f": {detail}" if detail else "."))
    data = target.read_bytes()
    if b"\r\n" in data:            # MinGW gfortran writes CRLF; CPMD runs on Linux
        target.write_bytes(data.replace(b"\r\n", b"\n"))
    return target


@dataclass
class JobOutcome:
    results: List[StageResult]           # completed stages, in order
    error: Optional[RunatomError]        # the stage error that stopped the job, or None
    failed_stage: Optional[Stage]        # the stage that raised it, or None
    formatted: Optional[Path]            # CPMD formatted file written after the generation stage, or None
    formatted_note: str                  # "" or why no formatted file was written


def clean_job_outputs(stages: Sequence[Stage], workdir, formatted_name: str) -> None:
    """Delete every product of every planned stage and the formatted file (missing files are fine)."""
    names = [name for stage in stages for name in stage.products()] + [formatted_name]
    _remove(Path(workdir), list(dict.fromkeys(names)))


def run_case_job(command: Sequence[str], reform_command: Optional[Sequence[str]], stages: Sequence[Stage],
                 workdir, formatted_name: str, progress: Optional[Callable[[dict], None]] = None,
                 cancel: Optional[threading.Event] = None, timeout: float = DEFAULT_STAGE_TIMEOUT) -> JobOutcome:
    """The window's job: clean, run every stage, convert for CPMD right after the generation.

    A stage error ends the job and is returned with the stages done so far, so a
    transferability test that fails never hides the pseudopotential already
    generated. RunatomCancelled is raised.
    """
    workdir = Path(workdir)
    clean_job_outputs(stages, workdir, formatted_name)
    outcome = JobOutcome([], None, None, None, "")
    for index, stage in enumerate(stages):
        try:
            result = run_stage(command, stage, workdir, progress, cancel, timeout, index, len(stages))
        except RunatomError as exc:
            outcome.error, outcome.failed_stage = exc, stage
            return outcome
        outcome.results.append(result)
        if stage.kind != "pseudo":
            continue
        if cancel is not None and cancel.is_set():
            raise RunatomCancelled(f"{stage.label}: cancelled")
        if reform_command is None:
            outcome.formatted_note = REFORM_MISSING_NOTE
            continue
        try:
            outcome.formatted = write_formatted(reform_command, workdir, stage.uspp_file, formatted_name, timeout)
        except (RunatomError, ValueError, OSError) as exc:
            outcome.formatted_note = f"The formatted file for CPMD was not written: {exc}"
    return outcome
