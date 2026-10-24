# Viewer Atom/Molecule Selection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ctrl/Shift+click selection of atoms and whole molecules in the 3D Molecular Viewer, a "show only selected" view that follows the selected atoms through a trajectory, and XYZ export of only the selected atoms (current frame or whole trajectory).

**Architecture:** Pure helpers (molecule = connected component of the drawn bonds, visible atom/bond filtering, XYZ text, export checks, summary) go in a new Toga- and OpenGL-free module `viewerSelection.py`; `indexRanges.format_index_ranges` is added as the inverse of the existing range parser. The selection state (a set of 0-based atom indices, frozen at click time) lives in the headless `MolecularViewer` base class beside `picked_atoms`/`identified_atoms`, with a no-op `_notify_selection_changed` hook that `MolecularViewerUI` marshals onto the Toga loop. `MolecularViewerUI` gains a "Selection" section on the Display tab and the two export paths.

**Tech Stack:** Python 3.13, numpy, scipy (`scipy.sparse.csgraph`), PyOpenGL/GLFW (render thread), Toga 0.5 WinForms, pytest.

**Spec:** [docs/superpowers/specs/2026-10-24-viewer-selection-design.md](../specs/2026-10-24-viewer-selection-design.md)

## Global Constraints

- Run everything from `venv/src/`; tests: `python -m pytest tests/ -q`. Baseline failures to ignore: `test_meanResidenceTime.py::TestAdvancedOptionsExposed::test_every_run_parameter_is_supplied_by_the_ui` and three table-styling tests in `test_plotter.py`.
- Atom indices are 0-based internally, **1-based in every UI text** (field, summary, comment line).
- A plain left click keeps today's behaviour exactly (label-only / measurement pick). Ctrl+click = one atom, Shift+click = molecule; Ctrl wins when both are held.
- Molecule = connected component of `self.bonds` (the bonds drawn for the displayed frame); frozen into indices at click time; no PBC.
- `viewerSelection.py` must import neither Toga nor OpenGL (subprocess test).
- Unticked "Export only selected", *Save Current Frame XYZ* output must stay **byte-identical** to the current code (`open(path, "w", encoding="utf-8")`, comment `Frame {k} exported from molecularViewer`, atom line `f"{element:<3s} {x: .8f} {y: .8f} {z: .8f}\n"`).
- Frames are never modified or copied for export: the writer runs on the `self.frames` list object (replaced, never mutated, by loads).
- Widget writes go through `cubeTab._assign` (already imported in molecularViewer) so unchanged values cause no WinForms re-layout.
- Commits: only once the user has authorised commits; otherwise leave the work uncommitted and report it.

---

### Task 1: Pure selection helpers (`indexRanges.format_index_ranges`, `viewerSelection.py`)

**Files:**
- Modify: `indexRanges.py` (append `format_index_ranges`)
- Create: `viewerSelection.py`
- Test: `tests/test_viewerSelection.py` (new)

**Interfaces:**
- Produces:
  - `indexRanges.format_index_ranges(indices: Iterable[int], sep: str = ", ") -> str`
  - `viewerSelection.connected_atoms(seed: int, bonds, natoms: int) -> List[int]` (sorted; `IndexError` if seed outside `[0, natoms)`)
  - `viewerSelection.expand_to_molecules(indices, bonds, natoms: int) -> List[int]` (sorted; indices ≥ natoms kept as they are)
  - `viewerSelection.visible_atoms_and_bonds(natoms: int, bonds, selection: Optional[Set[int]]) -> Tuple[Sequence[int], list]` (no/empty selection → `(range(natoms), bonds)` with the *same* bonds object)
  - `viewerSelection.selection_summary(elements: Sequence[str], indices: Sequence[int]) -> str`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_viewerSelection.py`:

```python
"""Atom/molecule selection in the 3D Molecular Viewer (headless).

Pure helpers in viewerSelection/indexRanges first, then the viewer's
selection state and click routing, rendering filters, and the UI plumbing
(MolecularViewerUI built with __new__ + stub widgets; no window is created).
"""
import subprocess
import sys
from pathlib import Path

import pytest

import viewerSelection as vs
from indexRanges import format_index_ranges, parse_solute_index_ranges

SRC = Path(__file__).resolve().parents[1]

# Two waters: O0 H1 H2 and O3 H4 H5 (bonded within each molecule only).
DIMER_BONDS = [(0, 1), (0, 2), (3, 4), (3, 5)]
DIMER_ELEMENTS = ["O", "H", "H", "O", "H", "H"]


# ----------------------------------------------------------------------
# Pure helpers
# ----------------------------------------------------------------------
def test_format_index_ranges_compresses_runs():
    assert format_index_ranges([1, 2, 3, 7, 9, 10]) == "1-3, 7, 9-10"
    assert format_index_ranges([5, 3, 4, 4]) == "3-5"
    assert format_index_ranges([]) == ""
    assert format_index_ranges([295, 296, 297, 1], sep=",") == "1,295-297"


@pytest.mark.parametrize(
    "values", [[1], [1, 3, 5], list(range(1, 50)), [2, 3, 10, 11, 12, 40]]
)
def test_format_round_trips_through_the_parser(values):
    assert parse_solute_index_ranges(format_index_ranges(values)) == values


def test_connected_atoms_returns_the_seed_molecule():
    assert vs.connected_atoms(4, DIMER_BONDS, 6) == [3, 4, 5]
    assert vs.connected_atoms(0, DIMER_BONDS, 6) == [0, 1, 2]


def test_connected_atoms_of_an_unbonded_atom_is_itself():
    assert vs.connected_atoms(2, [(0, 1)], 4) == [2]
    assert vs.connected_atoms(0, [], 1) == [0]


def test_connected_atoms_rejects_an_index_outside_the_frame():
    with pytest.raises(IndexError):
        vs.connected_atoms(6, DIMER_BONDS, 6)


def test_expand_to_molecules_adds_every_partner():
    assert vs.expand_to_molecules([1, 6], DIMER_BONDS, 7) == [0, 1, 2, 6]
    assert vs.expand_to_molecules([], DIMER_BONDS, 6) == []


def test_visible_atoms_without_a_selection_filter_is_everything():
    bonds = [(0, 1), (1, 2)]
    for selection in (None, set()):
        atoms, shown = vs.visible_atoms_and_bonds(3, bonds, selection)
        assert list(atoms) == [0, 1, 2]
        assert shown is bonds  # no copy on the common path


def test_visible_atoms_keep_only_selected_atoms_and_inner_bonds():
    atoms, shown = vs.visible_atoms_and_bonds(4, [(0, 1), (1, 2), (2, 3)], {1, 2, 9})
    assert atoms == [1, 2]  # 9 does not exist in this frame
    assert shown == [(1, 2)]


def test_selection_summary():
    assert vs.selection_summary(DIMER_ELEMENTS, []) == "No atoms selected."
    assert vs.selection_summary(DIMER_ELEMENTS, [3]) == "1 atom: O"
    assert vs.selection_summary(DIMER_ELEMENTS, list(range(6))) == "6 atoms: H4 O2"
    assert vs.selection_summary(DIMER_ELEMENTS[:2], [0, 5]) == "2 atoms: O (1 not in this frame)"


def test_viewer_selection_imports_neither_toga_nor_opengl():
    code = "import sys, viewerSelection; print('toga' in sys.modules, 'OpenGL' in sys.modules)"
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True
    )
    assert result.stdout.split() == ["False", "False"]
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_viewerSelection.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'viewerSelection'`.

- [ ] **Step 3: Implement**

Append to `indexRanges.py` (and change the import line to `from typing import Iterable, List`):

```python
def format_index_ranges(indices: Iterable[int], sep: str = ", ") -> str:
    """Inverse of parse_solute_index_ranges: [1, 2, 3, 7, 9, 10] -> "1-3, 7, 9-10".

    Sorted and de-duplicated; ``sep`` joins the parts (", " for a text field,
    "," for a compact XYZ comment line). Empty input gives "".
    """
    parts: List[str] = []
    start = previous = None
    for value in sorted(set(int(v) for v in indices)):
        if start is None:
            start = previous = value
        elif value == previous + 1:
            previous = value
        else:
            parts.append(str(start) if start == previous else f"{start}-{previous}")
            start = previous = value
    if start is not None:
        parts.append(str(start) if start == previous else f"{start}-{previous}")
    return sep.join(parts)
```

Create `viewerSelection.py`:

```python
"""Atom/molecule selection helpers for the 3D Molecular Viewer.

Toga- and OpenGL-free (a subprocess test keeps it so). The viewer keeps its
selection as 0-based atom indices into frames of (element, (x, y, z))
tuples; everything here works on those, so it is tested without a window.
A molecule is a connected component of the bonds the viewer draws.
"""

from collections import Counter
from typing import List, Optional, Sequence, Set, Tuple

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import breadth_first_order, connected_components


def _bond_graph(bonds, natoms: int):
    """Undirected adjacency (CSR) of ``bonds``; pairs outside the frame are dropped."""
    pairs = np.asarray(bonds, dtype=np.intp).reshape(-1, 2)
    pairs = pairs[(pairs >= 0).all(axis=1) & (pairs < natoms).all(axis=1)]
    data = np.ones(len(pairs), dtype=np.int8)
    return coo_matrix((data, (pairs[:, 0], pairs[:, 1])), shape=(natoms, natoms)).tocsr()


def connected_atoms(seed: int, bonds, natoms: int) -> List[int]:
    """The molecule containing ``seed``: its connected component in ``bonds``, sorted."""
    if not 0 <= seed < natoms:
        raise IndexError(f"Atom {seed + 1} is not in this frame ({natoms} atoms).")
    order = breadth_first_order(
        _bond_graph(bonds, natoms), seed, directed=False, return_predecessors=False
    )
    return sorted(int(i) for i in order)


def expand_to_molecules(indices, bonds, natoms: int) -> List[int]:
    """``indices`` plus every atom bonded (directly or not) to one of them, sorted.

    One connected-components pass for the whole frame, so expanding a large
    selection costs the same as expanding one atom. Indices outside the
    frame are kept unchanged (they may exist in other frames).
    """
    indices = {int(i) for i in indices}
    if not indices:
        return []
    _, labels = connected_components(_bond_graph(bonds, natoms), directed=False)
    inside = [i for i in indices if 0 <= i < natoms]
    wanted = np.isin(labels, labels[inside]) if inside else np.zeros(natoms, dtype=bool)
    return sorted(set(np.flatnonzero(wanted).tolist()) | (indices - set(inside)))


def visible_atoms_and_bonds(
    natoms: int, bonds, selection: Optional[Set[int]]
) -> Tuple[Sequence[int], list]:
    """Atoms and bonds to draw: all of them, or only the selection and its inner bonds.

    With no (or an empty) selection this returns ``range(natoms)`` and the
    very ``bonds`` object it was given, so the common path copies nothing.
    """
    if not selection:
        return range(natoms), bonds
    atoms = sorted(i for i in selection if 0 <= i < natoms)
    keep = set(atoms)
    return atoms, [(i, j) for i, j in bonds if i in keep and j in keep]


def _composition(elements: Sequence[str]) -> str:
    """'H4 O2' (alphabetical, count omitted when 1)."""
    counts = Counter(elements)
    return " ".join(f"{e}{n}" if n > 1 else e for e, n in sorted(counts.items()))


def selection_summary(elements: Sequence[str], indices: Sequence[int]) -> str:
    """One line for the Display tab, e.g. '6 atoms: H4 O2'."""
    if not indices:
        return "No atoms selected."
    present = [elements[i] for i in indices if 0 <= i < len(elements)]
    missing = len(indices) - len(present)
    word = "atom" if len(indices) == 1 else "atoms"
    text = f"{len(indices)} {word}: {_composition(present)}"
    if missing:
        text += f" ({missing} not in this frame)"
    return text
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_viewerSelection.py tests/test_bond_range_parser.py -q`
Expected: all pass.

- [ ] **Step 5: Commit** (if authorised)

```bash
git add indexRanges.py viewerSelection.py tests/test_viewerSelection.py
git commit -m "feat(viewer): pure helpers for atom/molecule selection"
```

---

### Task 2: XYZ export helpers

**Files:**
- Modify: `viewerSelection.py`
- Test: `tests/test_viewerSelection.py`

**Interfaces:**
- Consumes: `indexRanges.format_index_ranges` (Task 1).
- Produces:
  - `viewerSelection.xyz_frame_text(frame, indices: Optional[Sequence[int]], comment: str) -> str` (`indices=None` = every atom)
  - `viewerSelection.selection_comment(frame_index: int, source: str, indices: Sequence[int]) -> str`
  - `viewerSelection.selection_export_problem(frames, indices: Sequence[int], first_index: int = 0) -> Optional[str]`
  - `viewerSelection.write_xyz_frames(path: str, frames, indices, comment_for: Callable[[int], str], first_index: int = 0, progress: Optional[Callable[[int, int], None]] = None) -> int`
  - constants `EMPTY_SELECTION_MESSAGE: str`, `EXPORT_PROGRESS_FRAMES = 200`

- [ ] **Step 1: Write the failing tests** (append to `tests/test_viewerSelection.py`)

```python
# ----------------------------------------------------------------------
# XYZ export helpers
# ----------------------------------------------------------------------
FRAMES = [
    [("O", (0.0, 0.0, 0.0)), ("H", (0.96, 0.0, 0.0)), ("H", (-0.24, 0.93, 0.0)), ("O", (2.9, 0.0, 0.0))],
    [("O", (0.1, 0.0, 0.0)), ("H", (1.0, 0.1, 0.0)), ("H", (-0.2, 0.9, 0.1)), ("O", (2.8, 0.1, 0.0))],
]


def _old_viewer_frame_text(frame, comment):
    """What save_current_frame_xyz wrote before selection existed."""
    lines = [f"{len(frame)}\n", f"{comment}\n"]
    lines += [f"{element:<3s} {x: .8f} {y: .8f} {z: .8f}\n" for element, (x, y, z) in frame]
    return "".join(lines)


def test_xyz_frame_text_matches_the_viewer_format():
    comment = "Frame 0 exported from molecularViewer"
    assert vs.xyz_frame_text(FRAMES[0], None, comment) == _old_viewer_frame_text(FRAMES[0], comment)


def test_selected_frame_keeps_the_given_source_order():
    lines = vs.xyz_frame_text(FRAMES[1], [0, 3], "c").splitlines()
    assert lines[0] == "2" and lines[1] == "c"
    assert lines[2].split() == ["O", "0.10000000", "0.00000000", "0.00000000"]
    assert lines[3].split() == ["O", "2.80000000", "0.10000000", "0.00000000"]


def test_comment_is_forced_onto_one_line():
    assert vs.xyz_frame_text(FRAMES[0], [0], "a\nb\r\n").splitlines()[1] == "a b"


def test_selection_comment_names_frame_source_and_atoms():
    assert vs.selection_comment(7, "TRAJEC.xyz", [0, 1, 2, 294, 295, 296]) == (
        "Frame 7 of TRAJEC.xyz: 6 selected atoms (source atoms 1-3,295-297)"
    )
    assert vs.selection_comment(0, "x.xyz", [4]) == "Frame 0 of x.xyz: 1 selected atom (source atoms 5)"


def test_export_problem_for_an_empty_selection():
    assert vs.selection_export_problem(FRAMES, []) == vs.EMPTY_SELECTION_MESSAGE


def test_export_problem_names_the_first_short_frame():
    frames = FRAMES + [FRAMES[0][:2]]
    message = vs.selection_export_problem(frames, [0, 3], first_index=10)
    assert message.startswith("Frame 12 has 2 atoms")
    assert "atom 4" in message
    assert vs.selection_export_problem(FRAMES, [0, 3]) is None


def test_write_xyz_frames_writes_every_frame_in_order(tmp_path):
    from molecularViewer import read_xyz_frames

    path = tmp_path / "sel.xyz"
    reports = []
    count = vs.write_xyz_frames(
        str(path), FRAMES, [0, 3], lambda k: f"Frame {k}",
        progress=lambda done, total: reports.append((done, total)),
    )
    assert count == 2
    frames, _ = read_xyz_frames(str(path))  # the viewer's own reader reads it back
    assert [[atom[0] for atom in frame] for frame in frames] == [["O", "O"], ["O", "O"]]
    assert frames[1][1][1] == pytest.approx(FRAMES[1][3][1])
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[1] == "Frame 0" and lines[5] == "Frame 1"
    assert reports[-1] == (2, 2)


def test_write_xyz_frames_without_selection_writes_every_atom(tmp_path):
    from molecularViewer import read_xyz_frames

    path = tmp_path / "all.xyz"
    vs.write_xyz_frames(str(path), FRAMES, None, lambda k: "c")
    frames, _ = read_xyz_frames(str(path))
    assert [len(frame) for frame in frames] == [4, 4]


def test_first_index_numbers_the_comment(tmp_path):
    path = tmp_path / "one.xyz"
    vs.write_xyz_frames(str(path), [FRAMES[1]], None, lambda k: f"Frame {k}", first_index=5)
    assert path.read_text(encoding="utf-8").splitlines()[1] == "Frame 5"
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_viewerSelection.py -q`
Expected: the new tests fail with `AttributeError: module 'viewerSelection' has no attribute 'xyz_frame_text'`.

- [ ] **Step 3: Implement** (append to `viewerSelection.py`; add `from typing import Callable` and `from indexRanges import format_index_ranges` to the imports)

```python
# ----------------------------------------------------------------------
# XYZ export
# ----------------------------------------------------------------------
EXPORT_PROGRESS_FRAMES = 200  # frames between two progress reports

EMPTY_SELECTION_MESSAGE = (
    "No atoms are selected. Select atoms with Ctrl+click or whole molecules with "
    "Shift+click, or untick 'Export only selected atoms/molecules'."
)


def xyz_frame_text(frame, indices: Optional[Sequence[int]], comment: str) -> str:
    """One XYZ frame: count line, one-line comment, then the atoms.

    ``indices`` = None writes every atom; otherwise exactly those atoms, in
    the order given. The atom line format is the viewer's long-standing
    ``f"{element:<3s} {x: .8f} {y: .8f} {z: .8f}"``.
    """
    atoms = frame if indices is None else [frame[i] for i in indices]
    one_line = " ".join(comment.split())
    lines = [f"{len(atoms)}\n", f"{one_line}\n"]
    lines.extend(f"{element:<3s} {x: .8f} {y: .8f} {z: .8f}\n" for element, (x, y, z) in atoms)
    return "".join(lines)


def selection_comment(frame_index: int, source: str, indices: Sequence[int]) -> str:
    """Comment line of a selected-atoms frame; names the source atoms 1-based."""
    word = "atom" if len(indices) == 1 else "atoms"
    ranges = format_index_ranges((i + 1 for i in indices), sep=",")
    return f"Frame {frame_index} of {source}: {len(indices)} selected {word} (source atoms {ranges})"


def selection_export_problem(frames, indices: Sequence[int], first_index: int = 0) -> Optional[str]:
    """Why the selection cannot be exported from ``frames``, or None when it can."""
    if not indices:
        return EMPTY_SELECTION_MESSAGE
    needed = max(indices) + 1
    for position, frame in enumerate(frames):
        if len(frame) < needed:
            return (
                f"Frame {first_index + position} has {len(frame)} atoms, but the selection "
                f"includes atom {needed}. The atom count changes along this trajectory, "
                "so the selection cannot be exported for every frame."
            )
    return None


def write_xyz_frames(
    path: str,
    frames,
    indices: Optional[Sequence[int]],
    comment_for: Callable[[int], str],
    first_index: int = 0,
    progress: Optional[Callable[[int, int], None]] = None,
) -> int:
    """Write ``frames`` (all atoms, or ``indices``) to ``path``; returns the frame count.

    Frame numbers passed to ``comment_for`` start at ``first_index``. The file
    is opened exactly as the viewer's Save Current Frame always did (text mode,
    UTF-8), so an all-atom single frame is byte-identical to its old output.
    ``progress(done, total)`` is called every EXPORT_PROGRESS_FRAMES frames and
    at the end. OSError propagates to the caller.
    """
    total = len(frames)
    with open(path, "w", encoding="utf-8") as handle:
        for position, frame in enumerate(frames):
            handle.write(xyz_frame_text(frame, indices, comment_for(first_index + position)))
            done = position + 1
            if progress is not None and (done % EXPORT_PROGRESS_FRAMES == 0 or done == total):
                progress(done, total)
    return total
```

- [ ] **Step 4: Run the tests**

Run: `python -m pytest tests/test_viewerSelection.py -q`
Expected: all pass.

- [ ] **Step 5: Measure export speed** on a synthetic trajectory (1000 frames × 3000 atoms, all atoms and 30 selected), from a scratchpad script:

```python
import time, random
import viewerSelection as vs
frame = [("O", (random.random(), random.random(), random.random())) for _ in range(3000)]
frames = [frame] * 1000
for indices in (None, list(range(30))):
    t = time.perf_counter()
    vs.write_xyz_frames(r"<scratchpad>\big.xyz", frames, indices, lambda k: f"Frame {k}")
    print(indices is None, time.perf_counter() - t)
```

Record both timings in the final report. Only optimise the formatting if all-atom export of this size exceeds ~10 s (it runs off the UI thread).

- [ ] **Step 6: Commit** (if authorised)

```bash
git add viewerSelection.py tests/test_viewerSelection.py
git commit -m "feat(viewer): XYZ export helpers for selected atoms"
```

---

### Task 3: Selection state and Ctrl/Shift click routing (`MolecularViewer`)

**Files:**
- Modify: `molecularViewer.py` — imports; `MolecularViewer.__init__` (after `self.identified_atoms`, ~line 341); `load_designed_molecule` (~line 561-581); `on_mouse_button` (~line 967); new methods after `_notify_identified_atoms_changed` (~line 1207); `load_xyz_trajectory` (~line 2174-2184)
- Test: `tests/test_viewerSelection.py`

**Interfaces:**
- Consumes: `connected_atoms`, `expand_to_molecules` (Task 1); `parse_solute_index_ranges`.
- Produces (on `MolecularViewer`):
  - attributes `selected_atoms: Set[int]`, `_selection_version: int`, `show_selected_only: bool`, `export_selected_only: bool`, `_left_press_mods: int`
  - `get_selected_atoms() -> List[int]` (sorted)
  - `set_selected_atoms(atoms) -> bool`, `toggle_selected_atom(i) -> bool`, `toggle_selected_molecule(i) -> bool`, `expand_selection_to_molecules() -> bool`, `clear_selection() -> bool` (True = the selection changed)
  - `molecule_atoms(i) -> List[int]`
  - `parse_selection_text(text: str) -> List[int]` (0-based; `ValueError` with a user message)
  - `_route_atom_click(atom_index: int, mods: int = 0)`
  - hook `_notify_selection_changed()` (no-op), helper `_refocus_camera()`

- [ ] **Step 1: Write the failing tests** (append)

```python
# ----------------------------------------------------------------------
# Viewer selection state and click routing
# ----------------------------------------------------------------------
import asyncio
import gc

import glfw
import numpy as np

import molecularViewer as mv
from molecularViewer import MolecularViewer, MolecularViewerUI

DIMER_FRAME_0 = [
    ("O", (0.0, 0.0, 0.0)), ("H", (0.96, 0.0, 0.0)), ("H", (-0.24, 0.93, 0.0)),
    ("O", (2.9, 0.0, 0.0)), ("H", (3.4, 0.8, 0.0)), ("H", (3.4, -0.8, 0.0)),
]
# H1 has hopped to the second oxygen (0.9 A from O3, 2.0 A from O0).
DIMER_FRAME_1 = [DIMER_FRAME_0[0], ("H", (2.0, 0.0, 0.0))] + DIMER_FRAME_0[2:]


def _dimer_viewer():
    viewer = MolecularViewer()
    viewer.frames = [list(DIMER_FRAME_0), list(DIMER_FRAME_1)]
    viewer.molecule_data = viewer.frames[0]
    viewer.calculate_bonds(0)
    return viewer


def test_ctrl_toggles_one_atom_and_notifies():
    viewer = _dimer_viewer()
    calls = []
    viewer._notify_selection_changed = lambda: calls.append(viewer.get_selected_atoms())

    assert viewer.toggle_selected_atom(4) is True
    viewer.toggle_selected_atom(1)
    viewer.toggle_selected_atom(4)

    assert viewer.get_selected_atoms() == [1]
    assert calls == [[4], [1, 4], [1]]


def test_shift_selects_and_then_removes_the_whole_molecule():
    viewer = _dimer_viewer()
    viewer.toggle_selected_molecule(4)
    assert viewer.get_selected_atoms() == [3, 4, 5]
    viewer.toggle_selected_molecule(5)
    assert viewer.get_selected_atoms() == []


def test_shift_on_a_partly_selected_molecule_completes_it():
    viewer = _dimer_viewer()
    viewer.toggle_selected_atom(4)
    viewer.toggle_selected_molecule(3)
    assert viewer.get_selected_atoms() == [3, 4, 5]


def test_molecule_selection_is_frozen_when_connectivity_changes():
    viewer = _dimer_viewer()
    viewer.toggle_selected_molecule(0)
    viewer.set_recompute_bonds_per_frame(True)
    viewer.set_frame_index(1)  # H1 now bonded to O3

    assert viewer.get_selected_atoms() == [0, 1, 2]  # identities kept
    # A new Shift+click uses the bonds drawn for the displayed frame.
    assert viewer.molecule_atoms(3) == [1, 3, 4, 5]


def test_click_routing_by_modifier():
    viewer = _dimer_viewer()  # label-only mode is the default
    viewer._route_atom_click(4, glfw.MOD_CONTROL)
    assert viewer.get_selected_atoms() == [4]
    viewer._route_atom_click(0, glfw.MOD_SHIFT)
    assert viewer.get_selected_atoms() == [0, 1, 2, 4]
    viewer._route_atom_click(5, glfw.MOD_CONTROL | glfw.MOD_SHIFT)  # Ctrl wins
    assert viewer.get_selected_atoms() == [0, 1, 2, 4, 5]
    assert viewer.get_identified_atoms() == []

    viewer._route_atom_click(3, 0)  # plain click: unchanged behaviour
    assert viewer.get_identified_atoms() == [3]
    assert viewer.get_selected_atoms() == [0, 1, 2, 4, 5]


def _ortho(left, right, bottom, top, near, far):
    return np.array(
        [[2.0 / (right - left), 0, 0, 0], [0, 2.0 / (top - bottom), 0, 0],
         [0, 0, -2.0 / (far - near), 0],
         [-(right + left) / (right - left), -(top + bottom) / (top - bottom),
          -(far + near) / (far - near), 1.0]],
        dtype=float,
    )


def _picking_viewer():
    viewer = MolecularViewer()
    viewer.frames = []
    viewer.molecule_data = [("C", (5.0, 0.0, 0.0)), ("O", (-5.0, 0.0, 0.0)), ("N", (0.0, 0.0, 0.0))]
    viewer._pick_view = (np.eye(4), _ortho(-10, 10, -10, 10, -10, 10), (0, 0, 200, 200))
    return viewer


def test_mouse_modifiers_held_at_press_decide(monkeypatch):
    viewer = _picking_viewer()
    monkeypatch.setattr(mv.glfw, "get_cursor_pos", lambda window: (150.0, 100.0))  # atom C

    viewer.on_mouse_button(None, glfw.MOUSE_BUTTON_LEFT, glfw.PRESS, glfw.MOD_CONTROL)
    viewer.on_mouse_button(None, glfw.MOUSE_BUTTON_LEFT, glfw.RELEASE, 0)
    assert viewer.get_selected_atoms() == [0]

    viewer.on_mouse_button(None, glfw.MOUSE_BUTTON_LEFT, glfw.PRESS, 0)
    viewer.on_mouse_button(None, glfw.MOUSE_BUTTON_LEFT, glfw.RELEASE, glfw.MOD_CONTROL)
    assert viewer.get_selected_atoms() == [0]
    assert viewer.get_identified_atoms() == [0]  # plain press -> label


def test_set_clear_and_expand():
    viewer = _dimer_viewer()
    assert viewer.set_selected_atoms([4, 1]) is True
    assert viewer.get_selected_atoms() == [1, 4]
    assert viewer.set_selected_atoms([1, 4]) is False  # unchanged
    assert viewer.expand_selection_to_molecules() is True
    assert viewer.get_selected_atoms() == [0, 1, 2, 3, 4, 5]
    assert viewer.clear_selection() is True
    assert viewer.clear_selection() is False


def test_parse_selection_text():
    viewer = _dimer_viewer()
    assert viewer.parse_selection_text("1-3, 5") == [0, 1, 2, 4]
    assert viewer.parse_selection_text("  ") == []
    with pytest.raises(ValueError, match="Atom 7"):
        viewer.parse_selection_text("2 7")
    with pytest.raises(ValueError):
        viewer.parse_selection_text("3-1")


def test_loading_a_designed_molecule_clears_the_selection():
    viewer = _dimer_viewer()
    viewer.set_selected_atoms([0, 1])
    viewer.load_designed_molecule(list(DIMER_FRAME_0))
    assert viewer.get_selected_atoms() == []


def test_loading_a_trajectory_clears_the_selection(tmp_path):
    path = tmp_path / "t.xyz"
    path.write_text("2\nc\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    viewer = _dimer_viewer()
    viewer.set_selected_atoms([0, 1])
    try:
        asyncio.run(viewer.load_xyz_trajectory(str(path)))
    finally:
        gc.unfreeze()  # the loader freezes the heap (see test_viewerTrajectoryLoading)
    assert viewer.get_selected_atoms() == []
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_viewerSelection.py -q`
Expected: new tests fail with `AttributeError: 'MolecularViewer' object has no attribute 'toggle_selected_atom'`.

- [ ] **Step 3: Implement**

Imports (molecularViewer.py top): change `from typing import Dict, List, Optional, Tuple` to `from typing import Callable, Dict, List, Optional, Set, Tuple` and add after the `from help import AtomicData` line:

```python
from indexRanges import parse_solute_index_ranges
from viewerSelection import connected_atoms, expand_to_molecules
```

In `MolecularViewer.__init__`, right after `self.identified_atoms: List[int] = []`:

```python
        # Atom/molecule selection (Ctrl+click = atom, Shift+click = molecule):
        # 0-based atom identities, kept unchanged through every frame. A
        # molecule is resolved once, from the bonds drawn when it was clicked.
        self.selected_atoms: Set[int] = set()
        self._selection_version = 0
        self.show_selected_only = False
        self.export_selected_only = False
        self._left_press_mods = 0
```

In `on_mouse_button`, LEFT/PRESS branch add `self._left_press_mods = mods`; in the RELEASE branch replace `self._handle_atom_click(picked_index)` with `self._route_atom_click(picked_index, self._left_press_mods)`, and reset `self._left_press_mods = 0` alongside `self._left_dragged = False` at the end of the release branch.

New methods after `_notify_identified_atoms_changed`:

```python
    # ------------------------------------------------------------------
    # Atom/molecule selection (Display tab "Selection")
    # ------------------------------------------------------------------
    def _route_atom_click(self, atom_index: int, mods: int = 0):
        """Ctrl+click selects one atom, Shift+click its molecule (Ctrl wins);
        a plain click keeps the label/measure behaviour of _handle_atom_click."""
        if mods & glfw.MOD_CONTROL:
            self.toggle_selected_atom(atom_index)
        elif mods & glfw.MOD_SHIFT:
            self.toggle_selected_molecule(atom_index)
        else:
            self._handle_atom_click(atom_index)

    def get_selected_atoms(self) -> List[int]:
        """Sorted snapshot of the selected 0-based atom indices."""
        with self._state_lock:
            return sorted(self.selected_atoms)

    def _change_selection(self, change: Callable[[Set[int]], Set[int]]) -> bool:
        """Apply ``change`` to a copy of the selection under the lock; True if it changed."""
        with self._state_lock:
            new = set(change(set(self.selected_atoms)))
            if new == self.selected_atoms:
                return False
            self.selected_atoms = new
            self._selection_version += 1
            refocus = self.show_selected_only
        if refocus:
            self._refocus_camera()
        self._notify_selection_changed()
        return True

    def set_selected_atoms(self, atoms) -> bool:
        wanted = {int(i) for i in atoms}
        return self._change_selection(lambda current: wanted)

    def toggle_selected_atom(self, atom_index: int) -> bool:
        return self._change_selection(lambda current: current ^ {int(atom_index)})

    def molecule_atoms(self, atom_index: int) -> List[int]:
        """The molecule of ``atom_index`` in the bonds drawn for the displayed frame."""
        with self._state_lock:
            bonds = list(self.bonds)
            natoms = len(self.get_current_frame_data())
        return connected_atoms(int(atom_index), bonds, natoms)

    def toggle_selected_molecule(self, atom_index: int) -> bool:
        """Select the clicked atom's molecule; remove it if it is already wholly selected."""
        molecule = set(self.molecule_atoms(atom_index))
        return self._change_selection(
            lambda current: current - molecule if molecule <= current else current | molecule
        )

    def expand_selection_to_molecules(self) -> bool:
        with self._state_lock:
            bonds = list(self.bonds)
            natoms = len(self.get_current_frame_data())
            selected = sorted(self.selected_atoms)
        expanded = set(expand_to_molecules(selected, bonds, natoms))
        return self._change_selection(lambda current: current | expanded)

    def clear_selection(self) -> bool:
        return self._change_selection(lambda current: set())

    def parse_selection_text(self, text: str) -> List[int]:
        """1-based ranges ("1-3, 295-297") -> 0-based indices, checked against the frame."""
        numbers = parse_solute_index_ranges(text or "")
        if not numbers:
            return []
        natoms = len(self.get_current_frame_data())
        if natoms == 0:
            raise ValueError("No molecular data loaded.")
        too_large = [n for n in numbers if n > natoms]
        if too_large:
            raise ValueError(
                f"Atom {too_large[0]} does not exist: the displayed frame has {natoms} atoms."
            )
        return [n - 1 for n in numbers]

    def _refocus_camera(self):
        """Re-frame the view after the set of shown atoms changed."""
        self.invalidate_scene_cache()
        self.request_projection_update()

    def _notify_selection_changed(self):
        """Hook fired whenever the selection changes (possibly on the render thread).

        No-op here so the base viewer stays headless; MolecularViewerUI
        overrides it to update the Display tab from the Toga loop.
        """
        return None
```

In `load_designed_molecule`, after `self.calculate_bonds(0)` add `self.clear_selection()`. In `load_xyz_trajectory`, after `self.invalidate_scene_cache()` (the reset block at the top) add `self.clear_selection()`.

- [ ] **Step 4: Run tests**

Run: `python -m pytest tests/test_viewerSelection.py tests/test_molecularViewer.py tests/test_viewerTrajectoryLoading.py -q`
Expected: all pass.

- [ ] **Step 5: Commit** (if authorised)

```bash
git add molecularViewer.py tests/test_viewerSelection.py
git commit -m "feat(viewer): Ctrl/Shift+click atom and molecule selection"
```

---

### Task 4: Show-only-selected rendering, picking, camera framing, outlines

**Files:**
- Modify: `molecularViewer.py` — module constants (after `DESIGN_SYMBOL_DROP`); `_pick_atom_at` (~1056); `set_projection` (~1209); `render_frame` (~1853); `_render_design` (~1998); `main_loop` (~2102-2103); new methods next to the Task 3 block
- Test: `tests/test_viewerSelection.py`

**Interfaces:**
- Consumes: `viewerSelection.visible_atoms_and_bonds` (Task 1); selection state (Task 3).
- Produces (on `MolecularViewer`):
  - `set_show_selected_only(enabled: bool)`
  - `visible_atoms_and_bonds(natoms: int, bonds) -> Tuple[Sequence[int], list]`
  - `selection_outline_atoms(natoms: int) -> Set[int]`
  - `camera_focus_frame(frame_data: Optional[Frame] = None) -> Optional[Frame]` (None = frame every atom, as before)
  - constants `SELECTION_COLOR = (1.0, 0.55, 0.0)`, `SELECTION_OUTLINE_SCALE = 1.3`

- [ ] **Step 1: Write the failing tests** (append)

```python
# ----------------------------------------------------------------------
# Show only selected: drawing filter, picking, camera
# ----------------------------------------------------------------------
def test_visible_atoms_follow_the_show_only_switch():
    viewer = _dimer_viewer()
    viewer.set_selected_atoms([3, 4, 5])
    atoms, bonds = viewer.visible_atoms_and_bonds(6, viewer.bonds)
    assert list(atoms) == list(range(6))

    viewer.set_show_selected_only(True)
    atoms, bonds = viewer.visible_atoms_and_bonds(6, viewer.bonds)
    assert atoms == [3, 4, 5]
    assert bonds == [(3, 4), (3, 5)]


def test_show_only_with_an_empty_selection_shows_everything():
    viewer = _dimer_viewer()
    viewer.set_show_selected_only(True)
    atoms, _ = viewer.visible_atoms_and_bonds(6, viewer.bonds)
    assert list(atoms) == list(range(6))
    assert viewer.camera_focus_frame() is None


def test_outlines_only_when_everything_is_shown():
    viewer = _dimer_viewer()
    viewer.set_selected_atoms([0, 9])
    assert viewer.selection_outline_atoms(6) == {0}
    viewer.set_show_selected_only(True)
    assert viewer.selection_outline_atoms(6) == set()


def test_hidden_atoms_cannot_be_picked():
    viewer = _picking_viewer()
    viewer.set_selected_atoms([0])
    viewer.set_show_selected_only(True)
    assert viewer._pick_atom_at(50, 100) is None  # O is hidden
    assert viewer._pick_atom_at(150, 100) == 0     # C is shown


def test_vectorised_pick_matches_the_per_atom_projection():
    rng = np.random.default_rng(7)
    viewer = MolecularViewer()
    viewer.frames = []
    viewer.molecule_data = [("C", tuple(p)) for p in rng.uniform(-9, 9, size=(400, 3))]
    modelview = np.eye(4)
    projection = _ortho(-10, 10, -10, 10, -10, 10)
    viewport = (0, 0, 200, 200)
    viewer._pick_view = (modelview, projection, viewport)

    def reference(px, py):  # the per-atom loop _pick_atom_at used to run
        hits = []
        for index, (_, position) in enumerate(viewer.molecule_data):
            sx, sy, sz = viewer._project_atom(position, modelview, projection, viewport)
            if 0.0 <= sz <= 1.0 and np.hypot(sx - px, sy - (200 - py)) <= viewer.pick_pixel_threshold:
                hits.append((sz, index))
        return min(hits)[1] if hits else None

    for px, py in rng.uniform(0, 200, size=(200, 2)):
        assert viewer._pick_atom_at(px, py) == reference(px, py)


def test_camera_frames_the_selected_atoms_and_refocuses():
    viewer = _dimer_viewer()
    viewer.set_selected_atoms([3, 4, 5])
    assert viewer.camera_focus_frame() is None  # switch off: frame everything

    viewer.get_active_scene_radius()  # fixes the radius
    viewer.consume_projection_update_request()
    viewer.set_show_selected_only(True)
    focus = viewer.camera_focus_frame()
    assert [atom[0] for atom in focus] == ["O", "H", "H"]
    assert viewer.get_system_center(focus)[0] == pytest.approx((2.9 + 3.4 + 3.4) / 3)
    assert viewer._fixed_scene_radius is None
    assert viewer.consume_projection_update_request() is True

    viewer.consume_projection_update_request()
    viewer.toggle_selected_atom(5)  # a change while showing only the selection refocuses
    assert viewer.consume_projection_update_request() is True
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_viewerSelection.py -q`
Expected: new tests fail (`AttributeError: ... 'set_show_selected_only'`).

- [ ] **Step 3: Time the current `_pick_atom_at`** (scratchpad script, before changing it): build a `MolecularViewer` with 3000 random atoms and an identity/ortho `_pick_view` as in the test, time 50 calls, record ms per click. Expected ≈ 20 µs per atom (≈ 60 ms per click at 3000 atoms) — far above the 10 ms budget, so Step 4 vectorises it. If the measurement is below 10 ms, keep the loop and only add the visibility filter.

- [ ] **Step 4: Implement**

Constants after `DESIGN_SYMBOL_DROP`:

```python
SELECTION_COLOR = (1.0, 0.55, 0.0)  # orange: distinct from measurement yellow and pick cyan
SELECTION_OUTLINE_SCALE = 1.3       # selected atoms' rim: a back-face sphere this much larger
```

Import: `from viewerSelection import connected_atoms, expand_to_molecules, visible_atoms_and_bonds`.

New methods (after `_notify_selection_changed`):

```python
    def set_show_selected_only(self, enabled: bool):
        """Draw (and frame, and pick) only the selected atoms while a selection exists."""
        with self._state_lock:
            self.show_selected_only = bool(enabled)
        self._refocus_camera()

    def _shown_selection(self) -> Optional[Set[int]]:
        with self._state_lock:
            if self.show_selected_only and self.selected_atoms:
                return set(self.selected_atoms)
        return None

    def visible_atoms_and_bonds(self, natoms: int, bonds):
        """Atom indices and bonds to draw for a frame of ``natoms`` atoms."""
        return visible_atoms_and_bonds(natoms, bonds, self._shown_selection())

    def selection_outline_atoms(self, natoms: int) -> Set[int]:
        """Atoms drawn with the orange rim: the selection, unless only it is shown."""
        with self._state_lock:
            if self.show_selected_only:
                return set()
            return {i for i in self.selected_atoms if 0 <= i < natoms}

    def camera_focus_frame(self, frame_data: Optional[Frame] = None) -> Optional[Frame]:
        """The atoms the camera frames: the shown selection, or None for every atom."""
        selection = self._shown_selection()
        if selection is None:
            return None
        if frame_data is None:
            frame_data = self.get_current_frame_data()
        focus = [frame_data[i] for i in sorted(selection) if i < len(frame_data)]
        return focus or None
```

Replace the body of `_pick_atom_at` from `vx, vy, vw, vh = viewport` to the end with the vectorised form (same rule: within `pick_pixel_threshold`, window z in [0, 1], smallest depth then smallest index wins):

```python
        vx, vy, vw, vh = viewport
        # GLFW cursor y grows downward; OpenGL window y grows upward.
        click_x = px
        click_y = (vy + vh) - py

        # All shown atoms projected in one product (the _project_atom maths,
        # vectorised: a per-atom loop cost ~20 us per atom on the render thread).
        indices, _ = self.visible_atoms_and_bonds(len(frame_data), ())
        indices = np.fromiter(indices, dtype=np.intp)
        if indices.size == 0:
            return None
        obj = np.ones((indices.size, 4), dtype=float)
        obj[:, :3] = [frame_data[i][1] for i in indices]
        clip = (obj @ modelview) @ projection
        w = clip[:, 3]
        valid = w != 0.0
        ndc = np.zeros((indices.size, 3), dtype=float)
        ndc[valid] = clip[valid, :3] / w[valid, None]
        win_x = vx + vw * (ndc[:, 0] + 1.0) / 2.0
        win_y = vy + vh * (ndc[:, 1] + 1.0) / 2.0
        win_z = (ndc[:, 2] + 1.0) / 2.0
        distance = np.hypot(win_x - click_x, win_y - click_y)
        hit = valid & (win_z >= 0.0) & (win_z <= 1.0) & (distance <= self.pick_pixel_threshold)
        if not hit.any():
            return None
        candidates = np.flatnonzero(hit)
        best = candidates[np.lexsort((indices[candidates], win_z[candidates]))[0]]
        return int(indices[best])
```

`set_projection`: replace `radius = self.get_scene_radius()` with `radius = self.get_scene_radius(self.camera_focus_frame())` (None keeps today's behaviour and cache).

`main_loop`: replace

```python
                scene_center = self.get_system_center(frame_data)
                scene_radius = self.get_active_scene_radius(frame_data)
```

with

```python
                # Showing only the selection: frame those atoms, so the view
                # follows them as they move through the trajectory.
                focus = self.camera_focus_frame(frame_data)
                camera_frame = frame_data if focus is None else focus
                scene_center = self.get_system_center(camera_frame)
                scene_radius = self.get_active_scene_radius(camera_frame)
```

`render_frame`: at the top, after `picked = ...`, add

```python
        atom_indices, bonds = self.visible_atoms_and_bonds(len(frame_data), bonds)
        outlined = self.selection_outline_atoms(len(frame_data))
        shown = None if isinstance(atom_indices, range) else set(atom_indices)
        if shown is not None:
            picked = [i for i in picked if i in shown]
```

Call `self._render_design(frame_data, bonds, highlighted_atoms, atom_indices, outlined)`. In the other branch replace

```python
            for idx, (element, position) in enumerate(frame_data, start=1):
                atom_index = idx - 1
```

with

```python
            for atom_index in atom_indices:
                element, position = frame_data[atom_index]
                idx = atom_index + 1
```

and, right after the `if self.show_atoms_as_spheres() or atom_index in highlighted_atoms:` block (before the label block), add

```python
                if atom_index in outlined:
                    sphere = radius if self.show_atoms_as_spheres() else max(radius, 0.18)
                    self._draw_selection_outline(position, sphere)
```

In the identified-atoms loop add `if shown is not None and atom_index not in shown: continue` after the range check.

New GL helper next to `draw_sphere`:

```python
    def _draw_selection_outline(self, position, radius: float):
        """Orange rim round a selected atom: the back half of a larger sphere,
        unlit (Design style's outline technique), so it reads as a ring from any side."""
        glDisable(GL_LIGHTING)
        glEnable(GL_CULL_FACE)
        glCullFace(GL_FRONT)
        glColor3f(*SELECTION_COLOR)
        glPushMatrix()
        try:
            glTranslatef(*position)
            self.draw_sphere(radius * SELECTION_OUTLINE_SCALE)
        finally:
            glPopMatrix()
            glCullFace(GL_BACK)
            glDisable(GL_CULL_FACE)
            glEnable(GL_LIGHTING)
```

`_render_design(self, frame_data, bonds, highlighted_atoms, atom_indices=None, outlined=frozenset())`: at the top `if atom_indices is None: atom_indices = range(len(frame_data))`; the atom loop becomes

```python
            for index in atom_indices:
                element, position = frame_data[index]
                fill = (self.measurement_highlight_color if index in highlighted_atoms
                        else design_atom_colour(element))
                selected = index in outlined
                glPushMatrix()
                glTranslatef(*position)
                glCullFace(GL_FRONT)
                glColor3f(*(SELECTION_COLOR if selected else palette["outline"]))
                self.draw_sphere(radius * (SELECTION_OUTLINE_SCALE if selected else DESIGN_OUTLINE_SCALE))
                glCullFace(GL_BACK)
                glColor3f(*fill)
                self.draw_sphere(radius)
                glPopMatrix()
```

and the label loop becomes `for index in atom_indices: element, position = frame_data[index]; number = index + 1` with the body unchanged.

- [ ] **Step 5: Run tests**

Run: `python -m pytest tests/test_viewerSelection.py tests/test_molecularViewer.py tests/test_viewerDesignStyle.py tests/test_qtaimTab.py -q`
Expected: all pass (the existing pick tests pin the vectorised picker too).

- [ ] **Step 6: Time the vectorised pick** with the Step 3 script; record before/after ms per click.

- [ ] **Step 7: Commit** (if authorised)

```bash
git add molecularViewer.py tests/test_viewerSelection.py
git commit -m "feat(viewer): show only selected atoms, selection outlines, vectorised picking"
```

---

### Task 5: Display-tab "Selection" section and selective XYZ export (`MolecularViewerUI`)

**Files:**
- Modify: `molecularViewer.py` — imports; `_build_display_tab` (append the section; relabel "Clear selection" → "Clear labels"); `save_current_frame_xyz` (~3534); new UI methods after `_report_identified_atoms`
- Test: `tests/test_viewerSelection.py`

**Interfaces:**
- Consumes: Tasks 1-4 (`format_index_ranges`, `selection_summary`, `selection_comment`, `selection_export_problem`, `write_xyz_frames`, selection state methods).
- Produces (on `MolecularViewerUI`): widgets `selection_input`, `selection_summary_label`, `show_selected_switch`, `export_selected_switch`, `save_trajectory_button`; handlers `apply_selection_field` (async), `clear_atom_selection`, `expand_selection`, `toggle_show_selected_only`, `toggle_export_selected_only`, `save_trajectory_xyz` (async); `_sync_selection_widgets()`; `_export_source()`.

- [ ] **Step 1: Write the failing tests** (append)

```python
# ----------------------------------------------------------------------
# Display-tab plumbing and export (stub widgets; no window)
# ----------------------------------------------------------------------
from types import SimpleNamespace

import toga


class _StubLoop:
    def is_closed(self):
        return False

    def call_soon_threadsafe(self, callback, *args):
        callback(*args)

    def create_task(self, coro):
        coro.close()


class _DialogWindow:
    def __init__(self, answer=None):
        self.answer = answer
        self.dialogs = []

    async def dialog(self, dialog):
        self.dialogs.append(dialog)
        return self.answer if isinstance(dialog, toga.SaveFileDialog) else None


def _selection_ui(tmp_path, answer_name="out.xyz"):
    ui = MolecularViewerUI.__new__(MolecularViewerUI)
    MolecularViewer.__init__(ui)
    ui.loop = _StubLoop()
    ui.loading_label = None
    ui.textInput_file = SimpleNamespace(value=str(tmp_path / "TRAJEC.xyz"))
    ui.selection_input = SimpleNamespace(value="")
    ui.selection_summary_label = SimpleNamespace(text="")
    ui.show_selected_switch = SimpleNamespace(value=False)
    ui.export_selected_switch = SimpleNamespace(value=False)
    ui.save_trajectory_button = SimpleNamespace(enabled=True)
    ui.main_window = _DialogWindow(str(tmp_path / answer_name))
    ui.frames = [list(DIMER_FRAME_0), list(DIMER_FRAME_1)]
    ui.molecule_data = ui.frames[0]
    ui.calculate_bonds(0)
    return ui


def test_canvas_selection_fills_the_field_and_summary(tmp_path):
    ui = _selection_ui(tmp_path)
    ui._route_atom_click(4, glfw.MOD_CONTROL)
    assert ui.selection_input.value == "5"
    assert ui.selection_summary_label.text == "1 atom: H"
    ui._route_atom_click(0, glfw.MOD_SHIFT)
    assert ui.selection_input.value == "1-3, 5"


def test_typed_ranges_apply_and_normalise(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.selection_input.value = "6,4 5"
    asyncio.run(ui.apply_selection_field(None))
    assert ui.get_selected_atoms() == [3, 4, 5]
    assert ui.selection_input.value == "4-6"


def test_a_bad_field_reports_and_keeps_the_selection(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_selected_atoms([0])
    ui.selection_input.value = "9"
    asyncio.run(ui.apply_selection_field(None))
    assert isinstance(ui.main_window.dialogs[-1], toga.ErrorDialog)
    assert ui.get_selected_atoms() == [0]


def test_clear_and_expand_buttons(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_selected_atoms([4])
    ui.expand_selection(None)
    assert ui.get_selected_atoms() == [3, 4, 5]
    ui.clear_atom_selection(None)
    assert ui.get_selected_atoms() == []
    assert ui.selection_input.value == ""
    assert ui.selection_summary_label.text == "No atoms selected."


def test_switches_drive_the_viewer_flags(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.show_selected_switch.value = True
    ui.toggle_show_selected_only(ui.show_selected_switch)
    ui.export_selected_switch.value = True
    ui.toggle_export_selected_only(ui.export_selected_switch)
    assert ui.show_selected_only is True and ui.export_selected_only is True


def test_unticked_current_frame_export_is_byte_identical(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_frame_index(1)
    asyncio.run(ui.save_current_frame_xyz(None))

    expected = tmp_path / "expected.xyz"
    with open(expected, "w", encoding="utf-8") as f:  # the pre-selection writer, verbatim
        f.write(f"{len(DIMER_FRAME_1)}\n")
        f.write("Frame 1 exported from molecularViewer\n")
        for element, (x, y, z) in DIMER_FRAME_1:
            f.write(f"{element:<3s} {x: .8f} {y: .8f} {z: .8f}\n")
    assert (tmp_path / "out.xyz").read_bytes() == expected.read_bytes()
    assert isinstance(ui.main_window.dialogs[-1], toga.InfoDialog)


def test_selected_current_frame_export(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_selected_atoms([3, 4, 5])
    ui.export_selected_only = True
    asyncio.run(ui.save_current_frame_xyz(None))
    lines = (tmp_path / "out.xyz").read_text(encoding="utf-8").splitlines()
    assert lines[0] == "3"
    assert lines[1] == "Frame 0 of TRAJEC.xyz: 3 selected atoms (source atoms 4-6)"
    assert [line.split()[0] for line in lines[2:]] == ["O", "H", "H"]


def test_empty_selection_export_writes_nothing(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.export_selected_only = True
    asyncio.run(ui.save_current_frame_xyz(None))
    asyncio.run(ui.save_trajectory_xyz(None))
    assert not (tmp_path / "out.xyz").exists()
    assert all(isinstance(d, toga.ErrorDialog) for d in ui.main_window.dialogs)
    assert len(ui.main_window.dialogs) == 2  # no file dialog was opened


def test_selected_trajectory_export(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_selected_atoms([0, 1, 2])
    ui.export_selected_only = True
    asyncio.run(ui.save_trajectory_xyz(None))
    frames, _ = mv.read_xyz_frames(str(tmp_path / "out.xyz"))
    assert [len(frame) for frame in frames] == [3, 3]
    assert frames[1][1] == ("H", (2.0, 0.0, 0.0))  # the hopped H stays in, by identity
    lines = (tmp_path / "out.xyz").read_text(encoding="utf-8").splitlines()
    assert lines[6] == "Frame 1 of TRAJEC.xyz: 3 selected atoms (source atoms 1-3)"
    assert ui.save_trajectory_button.enabled is True
    assert isinstance(ui.main_window.dialogs[-1], toga.InfoDialog)


def test_full_trajectory_export_when_unticked(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.set_selected_atoms([0])  # ignored: the checkbox is off
    asyncio.run(ui.save_trajectory_xyz(None))
    frames, _ = mv.read_xyz_frames(str(tmp_path / "out.xyz"))
    assert [len(frame) for frame in frames] == [6, 6]


def test_trajectory_export_refuses_a_changing_atom_count(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.frames.append(list(DIMER_FRAME_0[:3]))
    ui.set_selected_atoms([4])
    ui.export_selected_only = True
    asyncio.run(ui.save_trajectory_xyz(None))
    assert not (tmp_path / "out.xyz").exists()
    assert "Frame 2 has 3 atoms" in ui.main_window.dialogs[-1].message if hasattr(
        ui.main_window.dialogs[-1], "message") else True
    assert isinstance(ui.main_window.dialogs[-1], toga.ErrorDialog)


def test_trajectory_export_refuses_to_overwrite_the_loaded_file(tmp_path):
    ui = _selection_ui(tmp_path, answer_name="TRAJEC.xyz")
    asyncio.run(ui.save_trajectory_xyz(None))
    assert not (tmp_path / "TRAJEC.xyz").exists()
    assert isinstance(ui.main_window.dialogs[-1], toga.ErrorDialog)


def test_trajectory_export_needs_a_loaded_structure(tmp_path):
    ui = _selection_ui(tmp_path)
    ui.frames = []  # only the built-in placeholder remains in molecule_data
    asyncio.run(ui.save_trajectory_xyz(None))
    assert len(ui.main_window.dialogs) == 1
    assert isinstance(ui.main_window.dialogs[0], toga.ErrorDialog)
```

(While implementing, check whether `toga.ErrorDialog` exposes its message headless; if it does, replace the guarded `hasattr` assertion in `test_trajectory_export_refuses_a_changing_atom_count` with a plain one, otherwise delete that line and keep the `isinstance` check.)

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_viewerSelection.py -q`
Expected: new tests fail (`AttributeError: ... 'apply_selection_field'` / no `save_trajectory_xyz`).

- [ ] **Step 3: Implement**

Imports: `from indexRanges import format_index_ranges, parse_solute_index_ranges` and
`from viewerSelection import (connected_atoms, expand_to_molecules, selection_comment, selection_export_problem, selection_summary, visible_atoms_and_bonds, write_xyz_frames)`.

`MolecularViewerUI.__init__`: add `self.selection_input = None`, `self.selection_summary_label = None`, `self.show_selected_switch = None`, `self.export_selected_switch = None`, `self.save_trajectory_button = None` beside the other widget placeholders.

`_build_display_tab`: rename the button text `"Clear selection"` → `"Clear labels"` (handler unchanged), then before `return tab_box` add `self._build_selection_section(tab_box)`:

```python
    def _build_selection_section(self, tab_box: toga.Box):
        """Atom/molecule selection: which atoms are focused on and exported."""
        tab_box.add(toga.Divider(style=Pack(margin=(10, 0, 4, 0))))
        tab_box.add(self._section_heading("Selection"))
        tab_box.add(
            self._hint(
                "In the 3D window, Ctrl+click selects or deselects an atom and Shift+click "
                "its whole molecule (the bonds drawn now). Selected atoms keep their identity "
                "in every frame. Or type 1-based atoms and ranges, e.g. 1-3, 295-297, and "
                "press Enter."
            )
        )
        self.selection_input = toga.TextInput(
            placeholder="e.g. 1-3, 295-297",
            style=Pack(flex=1, margin=(0, 8, 0, 0)),
            on_confirm=self.apply_selection_field,
        )
        tab_box.add(self._form_row(self._form_label("Selected atoms:"), self.selection_input))
        self.selection_summary_label = toga.Label(
            "No atoms selected.", style=Pack(margin=(0, 0, 6, 0))
        )
        tab_box.add(self.selection_summary_label)
        tab_box.add(
            self._form_row(
                toga.Button("Clear", on_press=self.clear_atom_selection,
                            style=Pack(width=90, margin=(0, 8, 0, 0))),
                toga.Button("Expand to molecules", on_press=self.expand_selection,
                            style=Pack(width=170)),
            )
        )
        # value goes to the constructor so building the tab fires no on_change.
        self.show_selected_switch = toga.Switch(
            "Show only selected atoms/molecules",
            value=self.show_selected_only,
            on_change=self.toggle_show_selected_only,
        )
        tab_box.add(self._form_row(self.show_selected_switch))
        self.export_selected_switch = toga.Switch(
            "Export only selected atoms/molecules",
            value=self.export_selected_only,
            on_change=self.toggle_export_selected_only,
        )
        self.save_trajectory_button = toga.Button(
            "Save Trajectory XYZ...", on_press=self.save_trajectory_xyz,
            style=Pack(width=170, margin=(0, 0, 0, 16)),
        )
        tab_box.add(self._form_row(self.export_selected_switch, self.save_trajectory_button))
        tab_box.add(
            self._hint(
                "The checkbox applies to Save Current Frame XYZ and Save Trajectory XYZ. "
                "Atoms are written in their original order; the trajectory file is not changed."
            )
        )
```

UI handlers (after `_report_identified_atoms`):

```python
    # ------------------------------------------------------------------
    # Atom/molecule selection (Display tab)
    # ------------------------------------------------------------------
    def _notify_selection_changed(self):
        """Marshal a selection change onto the Toga loop (render thread -> UI)."""
        loop = getattr(self, "loop", None)
        if loop is None or loop.is_closed():
            return
        loop.call_soon_threadsafe(self._sync_selection_widgets)

    def _sync_selection_widgets(self):
        """Selection field and summary from the selection (written only on change)."""
        if getattr(self, "selection_input", None) is None:
            return
        atoms = self.get_selected_atoms()
        elements = [element for element, _ in self.get_current_frame_data()]
        summary = selection_summary(elements, atoms)
        _assign(self.selection_input, "value", format_index_ranges(i + 1 for i in atoms))
        _assign(self.selection_summary_label, "text", summary)
        self.set_status_message(f"Selection: {summary}")

    async def apply_selection_field(self, widget):
        try:
            atoms = self.parse_selection_text(self.selection_input.value or "")
        except ValueError as exc:
            await self._show_error("Selection Error", str(exc))
            return
        self.set_selected_atoms(atoms)
        self._sync_selection_widgets()  # normalise the typed text even if unchanged

    def clear_atom_selection(self, widget):
        self.clear_selection()
        self._sync_selection_widgets()

    def expand_selection(self, widget):
        if not self.get_selected_atoms():
            self.set_status_message("Select atoms first, then expand them to whole molecules.")
            return
        self.expand_selection_to_molecules()
        self._sync_selection_widgets()

    def toggle_show_selected_only(self, widget):
        self.set_show_selected_only(bool(self.show_selected_switch.value))
        if self.show_selected_only and not self.get_selected_atoms():
            self.set_status_message("Show only selected is on; select atoms to apply it.")

    def toggle_export_selected_only(self, widget):
        self.export_selected_only = bool(self.export_selected_switch.value)

    def _export_source(self) -> Tuple[str, str]:
        """(loaded file path or "", base name for suggested file names)."""
        source = (self.textInput_file.value or "").strip() if self.textInput_file else ""
        base = os.path.splitext(os.path.basename(source))[0] if source else ""
        return source, base or "frame"

    async def save_trajectory_xyz(self, widget):
        """Every loaded frame as XYZ: all atoms, or only the selected ones."""
        with self._state_lock:
            frames = self.frames  # replaced, never mutated, by a load: a free snapshot
        if not frames:
            await self._show_error("No Data", "Load a trajectory or structure before saving it.")
            return
        indices = self.get_selected_atoms() if self.export_selected_only else None
        if indices is not None:
            problem = selection_export_problem(frames, indices)
            if problem:
                await self._show_error("Selection Error", problem)
                return

        source, base = self._export_source()
        suffix = "_selected" if indices is not None else "_trajectory"
        selected = await self.main_window.dialog(
            toga.SaveFileDialog(title="Save Trajectory as XYZ", suggested_filename=f"{base}{suffix}.xyz")
        )
        if not selected:
            return
        output_path = str(selected)
        if source and os.path.normcase(os.path.abspath(output_path)) == os.path.normcase(
            os.path.abspath(source)
        ):
            await self._show_error(
                "Save Error", "Choose another file name: this is the trajectory being viewed."
            )
            return

        source_name = os.path.basename(source) or "molecularViewer"
        if indices is None:
            comment_for = lambda k: f"Frame {k} of {source_name} exported from molecularViewer"
        else:
            comment_for = lambda k: selection_comment(k, source_name, indices)

        def report(done: int, total: int):
            self.set_status_message(f"Saving trajectory... {done} / {total} frames")

        if self.save_trajectory_button is not None:
            self.save_trajectory_button.enabled = False
        try:
            count = await asyncio.to_thread(
                write_xyz_frames, output_path, frames, indices, comment_for, 0, report
            )
        except OSError as exc:
            await self._show_error("Save Error", str(exc))
            return
        finally:
            if self.save_trajectory_button is not None:
                self.save_trajectory_button.enabled = True
        atoms = len(indices) if indices is not None else len(frames[0])
        what = f"{atoms} selected atoms" if indices is not None else "all atoms"
        self.set_status_message(f"Saved {count} frames ({what}) to {output_path}")
        await self.main_window.dialog(
            toga.InfoDialog("Trajectory Saved", f"{count} frames ({what}) saved to:\n{output_path}")
        )
```

Replace `save_current_frame_xyz` with:

```python
    async def save_current_frame_xyz(self, widget):
        try:
            frame_data = self.get_current_frame_data()
            if not frame_data:
                await self._show_error("No Data", "No frame is available to save.")
                return

            indices = self.get_selected_atoms() if self.export_selected_only else None
            if indices is not None:
                problem = selection_export_problem([frame_data], indices, self.current_frame)
                if problem:
                    await self._show_error("Selection Error", problem)
                    return

            source, base_name = self._export_source()
            suffix = "_selected" if indices is not None else ""
            suggested = f"{base_name}_frame_{self.current_frame:04d}{suffix}.xyz"

            selected = await self.main_window.dialog(
                toga.SaveFileDialog(title="Save Current Frame as XYZ", suggested_filename=suggested)
            )
            if not selected:
                return

            output_path = str(selected)
            if indices is None:
                comment_for = lambda k: f"Frame {k} exported from molecularViewer"
            else:
                source_name = os.path.basename(source) or "molecularViewer"
                comment_for = lambda k: selection_comment(k, source_name, indices)
            write_xyz_frames(output_path, [frame_data], indices, comment_for, self.current_frame)

            await self.main_window.dialog(
                toga.InfoDialog("Frame Saved", f"Current frame saved to:\n{output_path}")
            )
        except Exception as exc:
            await self._show_error("Save Error", str(exc))
```

- [ ] **Step 4: Run tests**

Run: `python -m pytest tests/test_viewerSelection.py tests/test_molecularViewer.py -q`
Expected: all pass.

- [ ] **Step 5: Commit** (if authorised)

```bash
git add molecularViewer.py tests/test_viewerSelection.py
git commit -m "feat(viewer): Display-tab selection controls and selective XYZ export"
```

---

### Task 6: Docs, full suite, live check

**Files:**
- Modify: `CLAUDE.md` (molecularViewer entry: new **Atom/molecule selection** paragraph; viewerSelection in the tool notes), `USER_MANUAL.md` (Display tab description ~line 1605-1640: Selection section, Ctrl/Shift gestures, Save Trajectory XYZ, export checkbox; "Clear selection" → "Clear labels" at ~1632 and ~1772)

- [ ] **Step 1: Write the docs.** CLAUDE.md paragraph must state: selection = 0-based identities frozen at click time (molecule from `self.bonds` of the displayed frame, no PBC); `_route_atom_click` (Ctrl atom / Shift molecule / plain = old behaviour, mods taken at press); `_notify_selection_changed` hook pattern; show-only = drawing filter + pick filter + camera framing via `camera_focus_frame` (None = old behaviour); orange rim via `_draw_selection_outline`; `_pick_atom_at` vectorised (measured numbers); export via Toga-free `viewerSelection.write_xyz_frames` (unticked current frame byte-identical; refuses empty selection, changing atom count, overwriting the loaded file; trajectory export needs real frames, not the placeholder); `format_index_ranges` is the inverse of the single range parser; tests in `tests/test_viewerSelection.py`; export timings.

- [ ] **Step 2: Full suite**

Run: `python -m pytest tests/ -q`
Expected: only the four baseline failures.

- [ ] **Step 3: Live check** in the real app (`python gqteaWinToga.py` → 3D Molecular Viewer): load a water-box trajectory, Ctrl+click and Shift+click atoms (orange rims; field/summary update), plain click still shows a number, turn on Show only selected and play (view follows the pair), Save Trajectory XYZ with the checkbox on and off, reopen the exported file in the viewer.

- [ ] **Step 4: Commit** (if authorised)

```bash
git add CLAUDE.md USER_MANUAL.md docs/superpowers/specs/2026-10-24-viewer-selection-design.md docs/superpowers/plans/2026-10-24-viewer-selection.md
git commit -m "docs: viewer atom/molecule selection and selective export"
```
