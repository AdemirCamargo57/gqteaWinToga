"""The main menu's tools, and loading each tool's module only when needed.

The main window used to import every tool module before it could appear:
2.5 s with a warm disk cache, 11 s cold, most of it numba, OpenGL and scipy
pulled in by the 3D Molecular Viewer. A menu button now names a module and a
class (``ToolSpec``); ``ToolLoader`` imports the module on one background
thread, either when the button is pressed or earlier, through a preload that
the app starts once its window is on screen. A pressed tool jumps the queue.

This module imports neither Toga nor any tool at module level (a subprocess
test checks it), so importing it costs nothing at startup.
"""

import collections
import importlib
import threading
from concurrent.futures import Future
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, Optional, Tuple


@dataclass(frozen=True)
class ToolSpec:
    """One menu button: its label and the class it opens (``module.attribute``)."""

    label: str
    module: str
    attribute: str


TOOL_TABS: Tuple[Tuple[str, Tuple[ToolSpec, ...]], ...] = (
    ("Geometry", (
        ToolSpec("Bond length analysis", "bond", "BondUI"),
        ToolSpec("All bond distance analysis", "allBondAnalysis", "allBondAnalysisUI"),
        ToolSpec("Bond angle analysis", "bondAngle", "BondAngleUI"),
        ToolSpec("All bond angle analysis", "allAnglesAnalysis", "allAnglesAnalysisUI"),
        ToolSpec("Dihedral angle analysis", "dihedralAngle", "DihedralUI"),
        ToolSpec("All dihedral angle analysis", "allDihedralAnalysis", "allDihedralAnalysisUI"),
        ToolSpec("Hydrogen bond analysis", "hbond", "HBondUI"),
        ToolSpec("Comparison of Molecular Geometric Parameters", "MolGeomComparator", "MolGeomComparatorUI"),
    )),
    ("Inputs", (
        ToolSpec("CPMD Inputs", "cpmdInput", "CpmdInputUI"),
        ToolSpec("SH Input Builder", "spectraSH", "SurfaceHoppingUI"),
        ToolSpec("Collision Input", "collision", "CollisionUI"),
        ToolSpec("cp.x Input Builder", "cpx_input_builder", "CPInputBuilderUI"),
        ToolSpec("ORCA Input Builder", "orca_input_builder", "ORCAInputBuilderUI"),
        ToolSpec("gqteaMD Input Builder", "gqteaMDinputBuilder", "GqteaMDInputBuilderUI"),
        ToolSpec("Gaussian Fukui Input Builder", "fukuiInputBuilder", "FukuiInputBuilderUI"),
        ToolSpec("Vanderbilt runatom.x Input Builder", "runatom_input_gen", "RunatomGenerator"),
    )),
    ("Structural", (
        ToolSpec("Radial Distribution Function", "radialDistribution", "RadialFunctionUI"),
        ToolSpec("Mean Residence Time Statistical", "meanResidenceTime_statistical", "MeanResidenceTimeUI"),
        ToolSpec("Mean Residence Time Classical", "meanResidenceTime_classical", "ClassicalMeanResidenceTimeUI"),
        ToolSpec("Autocorrelation function", "autocorrelationFunction", "AutoCorrelationFunctionUI"),
        ToolSpec("Single Solute solvent box", "single_solute_solvent_box", "SingleSoluteSolventBoxUI"),
        ToolSpec("Mixture of two solvent box", "mixtureSolventBox", "MixtureSolventBoxUI"),
        ToolSpec("Shared-wall double solvent box", "sharedWallBox", "SharedWallBoxUI"),
    )),
    ("Thermo", (
        ToolSpec("Classical rate constant", "classicalRate", "ClassicalRateUI"),
    )),
    ("Tools", (
        ToolSpec("3D Molecular Viewer", "molecularViewer", "MolecularViewerUI"),
        ToolSpec("Plots", "plotter", "PlotterUI"),
        ToolSpec("Analysis of Fukui Functions", "fukuiFunctionAnalyze", "FukuiFunctionAnalyzeUI"),
        ToolSpec("Molecular Axis Alignment", "molecularAxisAlignment", "MolecularAxisAlignmentUI"),
        ToolSpec("Select Frames", "selectFrames", "SelectFramesUI"),
        ToolSpec("Frame selection by interatomic distance range", "rangeFramesSelection", "RangeFramesSelectionUI"),
        ToolSpec("CPMD Input to XYZ Converter", "cpmdInputToXYZ", "CpmdInputToXYZUI"),
        ToolSpec("SH Geometry Analyzer", "sh_geom_analyzer", "SHGeomAnalyzerUI"),
        ToolSpec("Convert *.pos file to trajec.xyz", "cpTrajec2xyz", "CPtraj2xyzUI"),
        ToolSpec("Compute forces from cp.x *.for file", "cpxForces", "CPforcesUI"),
        ToolSpec("Coordinate converter", "coordinateConverter", "CoordinateConverterUI"),
        ToolSpec("Unit Converter", "unitConverter", "UnitConvertUI"),
    )),
)

# Imported before any tool. displayPlots pins matplotlib's Agg backend for the
# whole (Toga) process; when every tool was imported at startup that pin was
# always in place, and a tool that imports pyplot without displayPlots
# (sh_geom_analyzer) relied on it.
PREREQUISITE_MODULES: Tuple[str, ...] = ("displayPlots",)


def all_tools() -> Tuple[ToolSpec, ...]:
    """Every menu entry, in menu order."""
    return tuple(spec for _, specs in TOOL_TABS for spec in specs)


def tool_modules() -> Tuple[str, ...]:
    """Every tool module, in menu order, without repeats."""
    return tuple(dict.fromkeys(spec.module for spec in all_tools()))


def _frozen_build_imports():  # pragma: no cover - never called
    """Literal imports so a PyInstaller build (auto-py-to-exe) bundles every tool.

    ToolLoader imports the tools by name, which PyInstaller's static analysis
    cannot see; it does see these statements. This function is never called.
    tests/test_toolLoader.py keeps the list equal to ``TOOL_TABS``.
    """
    import displayPlots
    import bond, allBondAnalysis, bondAngle, allAnglesAnalysis, dihedralAngle
    import allDihedralAnalysis, hbond, MolGeomComparator
    import cpmdInput, spectraSH, collision, cpx_input_builder, orca_input_builder
    import gqteaMDinputBuilder, fukuiInputBuilder, runatom_input_gen
    import radialDistribution, meanResidenceTime_statistical, meanResidenceTime_classical
    import autocorrelationFunction, single_solute_solvent_box, mixtureSolventBox, sharedWallBox
    import classicalRate
    import molecularViewer, plotter, fukuiFunctionAnalyze, molecularAxisAlignment
    import selectFrames, rangeFramesSelection, cpmdInputToXYZ, sh_geom_analyzer
    import cpTrajec2xyz, cpxForces, coordinateConverter, unitConverter


class ToolLoader:
    """Imports tool modules on one background thread.

    ``request(name)`` returns a ``concurrent.futures.Future`` for the module and
    puts it at the front of the queue; ``start_preload()`` queues every tool
    module behind it. One thread does all the importing, so two imports never
    race on the same packages. A failed import is reported once and forgotten,
    so the next request tries again.

    ``on_loaded(name, error)`` is called **on the loader thread** after each
    tool module attempt (``error`` is None on success); a GUI must marshal it.
    """

    def __init__(
        self,
        modules: Iterable[str],
        prerequisites: Iterable[str] = PREREQUISITE_MODULES,
        on_loaded: Optional[Callable[[str, Optional[BaseException]], None]] = None,
        import_module: Callable[[str], object] = importlib.import_module,
    ):
        self.modules = tuple(dict.fromkeys(modules))
        self.prerequisites = tuple(prerequisites)
        self.on_loaded = on_loaded
        self._import_module = import_module
        self._lock = threading.Lock()
        self._queue = collections.deque()
        self._futures: Dict[str, Future] = {}
        self._loaded = set()
        self._failed = set()
        self._worker: Optional[threading.Thread] = None
        self._stopped = False

    # -- queries -------------------------------------------------------
    def is_loaded(self, name: str) -> bool:
        with self._lock:
            return name in self._loaded

    def loaded_count(self) -> int:
        with self._lock:
            return len(self._loaded.intersection(self.modules))

    def failed_count(self) -> int:
        with self._lock:
            return len(self._failed)

    def total_count(self) -> int:
        return len(self.modules)

    def is_busy(self) -> bool:
        with self._lock:
            return bool(self._queue) or self._worker is not None

    # -- requests ------------------------------------------------------
    def request(self, name: str) -> Future:
        """A future for module ``name``, imported next if it is not yet loaded."""
        with self._lock:
            future = self._futures.get(name)
            if future is None:
                future = self._futures[name] = Future()
            if not future.done():
                if name in self._queue:
                    self._queue.remove(name)
                self._queue.appendleft(name)
                self._stopped = False
                self._start_worker()
            return future

    def start_preload(self) -> None:
        """Queue every tool module not yet loaded, behind any pending request."""
        with self._lock:
            self._stopped = False
            for name in self.modules:
                if name not in self._loaded and name not in self._queue:
                    self._futures.setdefault(name, Future())
                    self._queue.append(name)
            self._start_worker()

    def stop(self) -> None:
        """Drop the queue; the import in progress (if any) still finishes."""
        with self._lock:
            self._stopped = True
            self._queue.clear()

    # -- worker --------------------------------------------------------
    def _start_worker(self) -> None:
        # Called with the lock held.
        if self._worker is None and self._queue:
            self._worker = threading.Thread(target=self._run, name="gqtea-tool-loader", daemon=True)
            self._worker.start()

    def _run(self) -> None:
        while True:
            with self._lock:
                if self._stopped or not self._queue:
                    self._worker = None
                    return
                name = self._queue.popleft()
                future = self._futures[name]
            try:
                for prerequisite in self.prerequisites:
                    self._import_module(prerequisite)
                module = self._import_module(name)
            except BaseException as error:  # an ImportError, a SyntaxError, anything
                with self._lock:
                    self._futures.pop(name, None)
                    self._failed.add(name)
                future.set_exception(error)
                self._report(name, error)
            else:
                with self._lock:
                    self._loaded.add(name)
                    self._failed.discard(name)
                future.set_result(module)
                self._report(name, None)

    def _report(self, name: str, error: Optional[BaseException]) -> None:
        if self.on_loaded is None or name not in self.modules:
            return
        try:
            self.on_loaded(name, error)
        except Exception:
            pass  # reporting must never stop the loader
