"""Atom/molecule selection helpers for the 3D Molecular Viewer.

Toga- and OpenGL-free (a subprocess test keeps it so). The viewer keeps its
selection as 0-based atom indices into frames of (element, (x, y, z))
tuples; everything here works on those, so it is tested without a window.
A molecule is a connected component of the bonds the viewer draws.
"""

from collections import Counter
from typing import Callable, List, Optional, Sequence, Set, Tuple

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import breadth_first_order, connected_components

from indexRanges import format_index_ranges


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
    inside = [i for i in indices if 0 <= i < natoms]
    expanded = set(indices) - set(inside)
    if inside:
        _, labels = connected_components(_bond_graph(bonds, natoms), directed=False)
        expanded.update(np.flatnonzero(np.isin(labels, labels[inside])).tolist())
    return sorted(expanded)


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
