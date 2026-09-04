"""Tests for the All Bond Distance Analysis tool (``allBondAnalysis.py``).

Run from ``venv/src/``::

    python -m pytest tests/ -q

Only the logic class ``AllBondAnalysis`` is exercised (no Toga GUI is created),
following the same pattern as the RDF tests. Statistics are checked against
hand-computed analytic references rather than a frozen golden file.

Bond connectivity is re-evaluated **every frame**: a pair contributes to its
mean/variance only in the frames where it is within the cutoff, and the tool
reports how often each pair is bonded (occurrence). Periodic boundaries are
optional (orthorhombic minimum-image) when a cell is supplied.
"""
import os
import math
import asyncio

import numpy as np
import pytest

from allBondAnalysis import AllBondAnalysis, BondPair


# --------------------------------------------------------------------------- #
# Fixture builders                                                              #
# --------------------------------------------------------------------------- #
def write_xyz(path, frames, symbols, comment="frame"):
    """Write a standard multi-frame XYZ file."""
    with open(path, "w") as f:
        for coords in frames:
            f.write(f"{len(symbols)}\n{comment}\n")
            for s, (x, y, z) in zip(symbols, coords):
                f.write(f"{s} {x:.6f} {y:.6f} {z:.6f}\n")


def two_atom_frames(separations):
    """Two atoms on the x-axis at the given per-frame separations."""
    return [np.array([[0.0, 0.0, 0.0], [d, 0.0, 0.0]]) for d in separations]


# --------------------------------------------------------------------------- #
# Frame reading & first-frame metadata (unchanged behavior)                     #
# --------------------------------------------------------------------------- #
def test_read_first_frame_sets_metadata(tmp_path):
    traj = tmp_path / "t.xyz"
    write_xyz(traj, two_atom_frames([1.0]), ["O", "O"])

    analyzer = AllBondAnalysis(str(traj))
    elements, coords = analyzer.read_first_frame()

    assert elements == ["O", "O"]
    assert analyzer.num_atoms == 2
    assert coords.shape == (2, 3)
    np.testing.assert_allclose(coords[1], [1.0, 0.0, 0.0])


def test_read_frames_counts_all_frames(tmp_path):
    traj = tmp_path / "t.xyz"
    write_xyz(traj, two_atom_frames([1.0, 1.1, 1.2]), ["O", "O"])
    analyzer = AllBondAnalysis(str(traj))
    assert len(list(analyzer._read_frames())) == 3


# --------------------------------------------------------------------------- #
# Default maximum connection distance                                           #
# --------------------------------------------------------------------------- #
def test_default_connection_distance_is_1_7():
    assert AllBondAnalysis.DEFAULT_CONNECTION_DISTANCE == pytest.approx(1.7)
    assert AllBondAnalysis("x.xyz").max_connection_distance == pytest.approx(1.7)


# --------------------------------------------------------------------------- #
# Per-frame statistics, occurrence, and the frozen-frame fix                    #
# --------------------------------------------------------------------------- #
def test_statistics_over_all_bonded_frames(tmp_path):
    """A pair within cutoff every frame: stats over all frames, occurrence 1.0."""
    traj = tmp_path / "pair.xyz"
    write_xyz(traj, two_atom_frames([1.0, 1.1, 1.2]), ["O", "O"])

    analyzer = AllBondAnalysis(str(traj), max_connection_distance=1.5)
    stats = analyzer.analyze()

    assert len(stats) == 1
    mean, variance, std = stats[0]
    assert mean == pytest.approx(1.1)
    assert variance == pytest.approx(0.02 / 3.0)
    assert std == pytest.approx(math.sqrt(0.02 / 3.0))
    pair = analyzer.connected_atom_pairs[0]
    assert pair.frames_bonded == 3
    assert analyzer.num_frames == 3


def test_statistics_only_over_bonded_frames(tmp_path):
    """A pair that leaves the cutoff mid-run contributes only its bonded frames."""
    traj = tmp_path / "pair.xyz"
    write_xyz(traj, two_atom_frames([1.0, 2.0, 1.2]), ["O", "O"])  # frame 2 not bonded

    analyzer = AllBondAnalysis(str(traj), max_connection_distance=1.5)
    stats = analyzer.analyze()

    mean, variance, std = stats[0]
    assert mean == pytest.approx(1.1)           # (1.0 + 1.2) / 2, frame 2 excluded
    assert variance == pytest.approx(0.01)      # ((-0.1)^2 + 0.1^2) / 2
    assert analyzer.connected_atom_pairs[0].frames_bonded == 2
    assert analyzer.num_frames == 3


def test_bond_forming_later_is_captured(tmp_path):
    """Regression for the frozen-frame bug: a pair NOT bonded in frame 1 is still
    analyzed once it comes within the cutoff in a later frame."""
    traj = tmp_path / "pair.xyz"
    write_xyz(traj, two_atom_frames([3.0, 1.0, 1.1]), ["O", "O"])  # unbonded frame 1

    analyzer = AllBondAnalysis(str(traj), max_connection_distance=1.5)
    stats = analyzer.analyze()

    assert len(analyzer.connected_atom_pairs) == 1
    assert stats[0][0] == pytest.approx(1.05)   # mean of frames 2 and 3
    assert analyzer.connected_atom_pairs[0].frames_bonded == 2


def test_first_bonded_distance_recorded(tmp_path):
    """The reference distance is the first frame in which the pair is bonded."""
    traj = tmp_path / "pair.xyz"
    write_xyz(traj, two_atom_frames([3.0, 1.0, 1.1]), ["O", "O"])
    analyzer = AllBondAnalysis(str(traj), max_connection_distance=1.5)
    analyzer.analyze()
    assert analyzer.connected_atom_pairs[0].first_bonded_distance == pytest.approx(1.0)


def test_solute_scope_restricts_pairs(tmp_path):
    """With solute indices {0,1}, only the A-B pair is considered (C excluded)."""
    traj = tmp_path / "chain.xyz"
    coords = np.array([[0.0, 0, 0], [1.0, 0, 0], [1.4, 0, 0]])  # A-B, B-C both <=1.5
    write_xyz(traj, [coords], ["A", "B", "C"])

    analyzer = AllBondAnalysis(str(traj), max_connection_distance=1.5,
                               solute_atom_indices=[0, 1])
    analyzer.analyze()
    pairs = [(p.atom_i, p.atom_j) for p in analyzer.connected_atom_pairs]
    assert pairs == [(0, 1)]


def test_no_bonded_pairs_raises(tmp_path):
    """If no pair is ever within the cutoff, statistics cannot be produced."""
    traj = tmp_path / "pair.xyz"
    write_xyz(traj, two_atom_frames([3.0, 3.0]), ["O", "O"])
    analyzer = AllBondAnalysis(str(traj), max_connection_distance=1.5)
    with pytest.raises(ValueError, match="No connected atom pairs"):
        analyzer.analyze()


# --------------------------------------------------------------------------- #
# Periodic boundary conditions (optional orthorhombic minimum-image)            #
# --------------------------------------------------------------------------- #
def test_pbc_detects_bond_across_boundary(tmp_path):
    """Atoms near opposite faces are bonded under the minimum-image convention."""
    traj = tmp_path / "pbc.xyz"
    coords = np.array([[0.5, 0, 0], [9.7, 0, 0]])  # raw 9.2 A, min-image 0.8 A
    write_xyz(traj, [coords], ["O", "O"])

    with_pbc = AllBondAnalysis(str(traj), max_connection_distance=1.5,
                               cell_lengths=[10.0, 10.0, 10.0])
    stats = with_pbc.analyze()
    assert with_pbc.connected_atom_pairs[0].first_bonded_distance == pytest.approx(0.8)
    assert stats[0][0] == pytest.approx(0.8)

    # Without a cell the same pair is 9.2 A apart -> not bonded.
    without = AllBondAnalysis(str(traj), max_connection_distance=1.5)
    with pytest.raises(ValueError, match="No connected atom pairs"):
        without.analyze()


def test_pbc_box_shift_invariance(tmp_path):
    """Translating every atom by one full box must not change the statistics."""
    base = tmp_path / "base.xyz"
    shifted = tmp_path / "shifted.xyz"
    coords = np.array([[0.5, 0, 0], [9.7, 0, 0]])
    write_xyz(base, [coords, coords], ["O", "O"])
    write_xyz(shifted, [coords + np.array([10.0, 0, 0]), coords + np.array([10.0, 0, 0])], ["O", "O"])

    box = [10.0, 10.0, 10.0]
    a = AllBondAnalysis(str(base), max_connection_distance=1.5, cell_lengths=box)
    b = AllBondAnalysis(str(shifted), max_connection_distance=1.5, cell_lengths=box)
    assert a.analyze()[0] == pytest.approx(b.analyze()[0])


def test_pbc_rejects_cutoff_above_half_box(tmp_path):
    traj = tmp_path / "pbc.xyz"
    write_xyz(traj, two_atom_frames([1.0]), ["O", "O"])
    analyzer = AllBondAnalysis(str(traj), max_connection_distance=2.0,
                               cell_lengths=[3.0, 3.0, 3.0])  # half box = 1.5 < 2.0
    with pytest.raises(ValueError, match="half"):
        analyzer.analyze()


# --------------------------------------------------------------------------- #
# Async path: identical result + progress callback                              #
# --------------------------------------------------------------------------- #
def test_async_matches_sync_and_reports_progress(tmp_path):
    traj = tmp_path / "pair.xyz"
    write_xyz(traj, two_atom_frames([1.0, 1.1, 1.2]), ["O", "O"])

    sync_stats = AllBondAnalysis(str(traj), max_connection_distance=1.5).analyze()

    calls = []
    a = AllBondAnalysis(str(traj), max_connection_distance=1.5)
    async_stats = asyncio.run(
        a.analyze_async(progress_callback=lambda done, total: calls.append((done, total)),
                        ui_update_interval=1)
    )
    assert async_stats == pytest.approx(sync_stats)
    assert calls[-1] == (3, 3)


# --------------------------------------------------------------------------- #
# Output files: statistics + row-to-pair mapping (now with occurrence)          #
# --------------------------------------------------------------------------- #
def test_write_results_single_combined_file(tmp_path):
    """Both files are merged into one self-describing file: metadata header plus
    one row per pair carrying identity, first-bonded distance, stats, occurrence."""
    traj = tmp_path / "chain.xyz"
    coords = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]])  # A-B, B-C in; A-C (2.0) out
    write_xyz(traj, [coords, coords], ["A", "B", "C"])

    analyzer = AllBondAnalysis(str(traj), max_connection_distance=1.5)
    analyzer.analyze()
    output_file = analyzer.write_results(str(tmp_path / "combined.txt"))

    # A single path is returned and no separate *_pairs.txt is created.
    assert isinstance(output_file, str)
    assert os.path.exists(output_file)
    assert not os.path.exists(str(tmp_path / "combined_pairs.txt"))

    text = open(output_file).read()
    assert "# frames_used 2" in text
    assert "# max_connection_distance 1.5" in text
    assert "# atom_scope all_atoms" in text
    assert "# periodic_boundary none" in text

    column_header = [l for l in text.splitlines() if l.startswith("# row")][0]
    for col in ("atom_i", "first_bonded_distance", "average",
                "occurrence_fraction", "frames_bonded"):
        assert col in column_header

    rows = [l for l in text.splitlines() if l and not l.startswith("#")]
    assert len(rows) == 2                                   # A-B and B-C
    first = rows[0].split()
    assert first[1:5] == ["1", "2", "A", "B"]               # identity, 1-based labels
    assert float(first[5]) == pytest.approx(1.0)            # first_bonded_distance (A-B)
    assert float(first[6]) == pytest.approx(1.0)            # average
    assert float(first[9]) == pytest.approx(1.0)            # occurrence_fraction
    assert int(first[10]) == 2                              # frames_bonded


def test_write_results_records_cell_and_solute_metadata(tmp_path):
    traj = tmp_path / "chain.xyz"
    coords = np.array([[0.0, 0, 0], [1.0, 0, 0], [1.4, 0, 0]])
    write_xyz(traj, [coords], ["A", "B", "C"])

    analyzer = AllBondAnalysis(str(traj), max_connection_distance=1.5,
                               solute_atom_indices=[0, 1], cell_lengths=[10.0, 10.0, 10.0])
    analyzer.analyze()
    output_file = analyzer.write_results(str(tmp_path / "s.txt"))

    text = open(output_file).read()
    assert "# atom_scope solute_atoms" in text
    assert "# solute_atom_indices 1 2" in text
    assert "# cell_lengths 10" in text


def test_write_results_without_statistics_raises(tmp_path):
    analyzer = AllBondAnalysis(str(tmp_path / "none.xyz"))
    with pytest.raises(ValueError, match="No statistics"):
        analyzer.write_results(str(tmp_path / "out.txt"))


# --------------------------------------------------------------------------- #
# Input validation & malformed trajectories                                     #
# --------------------------------------------------------------------------- #
def test_missing_trajectory_file_raises():
    analyzer = AllBondAnalysis(trajectory_file=None)
    with pytest.raises(ValueError, match="No trajectory file"):
        list(analyzer._read_frames())


def test_empty_trajectory_raises(tmp_path):
    traj = tmp_path / "empty.xyz"
    traj.write_text("")
    analyzer = AllBondAnalysis(str(traj))
    with pytest.raises(ValueError, match="does not contain any frames"):
        analyzer.read_first_frame()


def test_inconsistent_atom_count_raises(tmp_path):
    traj = tmp_path / "bad.xyz"
    with open(traj, "w") as f:
        f.write("2\nc\nA 0 0 0\nB 1 0 0\n")
        f.write("3\nc\nA 0 0 0\nB 1 0 0\nC 2 0 0\n")
    analyzer = AllBondAnalysis(str(traj), max_connection_distance=1.5)
    with pytest.raises(ValueError, match="Inconsistent number of atoms"):
        analyzer.analyze()


def test_truncated_frame_raises(tmp_path):
    traj = tmp_path / "trunc.xyz"
    with open(traj, "w") as f:
        f.write("2\nc\nA 0 0 0\n")
    analyzer = AllBondAnalysis(str(traj))
    with pytest.raises(ValueError, match="Unexpected end of file"):
        analyzer.read_first_frame()


def test_invalid_coordinates_raise(tmp_path):
    traj = tmp_path / "coords.xyz"
    with open(traj, "w") as f:
        f.write("2\nc\nA 0 0 0\nB x 0 0\n")
    analyzer = AllBondAnalysis(str(traj))
    with pytest.raises(ValueError, match="Invalid coordinates"):
        analyzer.read_first_frame()


@pytest.mark.parametrize("cell", [
    [10.0, 10.0],            # not three lengths
    [10.0, 0.0, 10.0],       # non-positive
])
def test_invalid_cell_raises(tmp_path, cell):
    traj = tmp_path / "pair.xyz"
    write_xyz(traj, two_atom_frames([1.0]), ["O", "O"])
    with pytest.raises(ValueError):
        AllBondAnalysis(str(traj), cell_lengths=cell).analyze()


@pytest.mark.parametrize("indices,message", [
    ([], "cannot be empty"),
    ([0, 0], "duplicates"),
    ([5], "between"),
    ([-1], "between"),
])
def test_solute_index_validation(tmp_path, indices, message):
    traj = tmp_path / "pair.xyz"
    write_xyz(traj, two_atom_frames([1.0, 1.1]), ["O", "O"])
    analyzer = AllBondAnalysis(str(traj), solute_atom_indices=indices)
    with pytest.raises(ValueError, match=message):
        analyzer.read_first_frame()
