"""Tests for the All Bond Angle Analysis tool (``allAnglesAnalysis.py``).

Run from ``venv/src/``. Only the logic class ``AllAnglesAnalysis`` is exercised
(no Toga GUI). Mirrors the all-bond tool: connectivity is re-evaluated every
frame, each i-j-k angle is averaged only over the frames where both of its
bonds (j-i and j-k) are within the cutoff, occurrence is reported, PBC is
optional (orthorhombic minimum-image), the connection distance defaults to 1.7,
solute indices accept range syntax, and results go to one combined file.
"""
import os
import math
import asyncio

import numpy as np
import pytest

from allAnglesAnalysis import AllAnglesAnalysis, AngleTriplet


def write_xyz(path, frames, symbols, comment="frame"):
    with open(path, "w") as f:
        for coords in frames:
            f.write(f"{len(symbols)}\n{comment}\n")
            for s, (x, y, z) in zip(symbols, coords):
                f.write(f"{s} {x:.6f} {y:.6f} {z:.6f}\n")


def bent_frames(d2_per_frame, d1=1.0):
    """H(0)-O(1)-H(2) with a fixed 90 deg angle; O-H2 distance varies per frame.

    O-H1 = d1 (constant), O-H2 = d2 (per frame). i-k stays > 1.5 A so no spurious
    bond closes the triangle.
    """
    return [
        np.array([[d1, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, d2, 0.0]])
        for d2 in d2_per_frame
    ], ["H", "O", "H"]


# --------------------------------------------------------------------------- #
# Metadata + default                                                            #
# --------------------------------------------------------------------------- #
def test_read_first_frame_sets_metadata(tmp_path):
    traj = tmp_path / "t.xyz"
    frames, symbols = bent_frames([1.4])
    write_xyz(traj, frames, symbols)
    a = AllAnglesAnalysis(str(traj))
    elements, coords = a.read_first_frame()
    assert elements == ["H", "O", "H"]
    assert a.num_atoms == 3
    assert coords.shape == (3, 3)


def test_default_connection_distance_is_1_7():
    assert AllAnglesAnalysis.DEFAULT_CONNECTION_DISTANCE == pytest.approx(1.7)
    assert AllAnglesAnalysis("x.xyz").max_connection_distance == pytest.approx(1.7)


# --------------------------------------------------------------------------- #
# Per-frame occurrence and the frozen-frame fix                                 #
# --------------------------------------------------------------------------- #
def test_angle_present_every_frame(tmp_path):
    traj = tmp_path / "a.xyz"
    frames, symbols = bent_frames([1.4, 1.4])
    write_xyz(traj, frames, symbols)
    a = AllAnglesAnalysis(str(traj), max_connection_distance=1.5)
    stats = a.analyze()
    assert len(stats) == 1
    mean, variance, std = stats[0]
    assert mean == pytest.approx(90.0)
    assert variance == pytest.approx(0.0)
    triplet = a.angle_triplets[0]
    assert (triplet.atom_i, triplet.atom_j, triplet.atom_k) == (0, 1, 2)
    assert triplet.frames_present == 2
    assert a.num_frames == 2


def test_angle_only_over_present_frames(tmp_path):
    """Frame 2 breaks the O-H2 bond, so that frame does not contribute."""
    traj = tmp_path / "a.xyz"
    frames, symbols = bent_frames([1.4, 2.0, 1.4])  # frame 2: O-H2 = 2.0 > 1.5
    write_xyz(traj, frames, symbols)
    a = AllAnglesAnalysis(str(traj), max_connection_distance=1.5)
    a.analyze()
    assert a.angle_triplets[0].frames_present == 2
    assert a.num_frames == 3


def test_angle_forming_later_is_captured(tmp_path):
    """Regression for the frozen-frame bug: a triplet absent in frame 1 is
    still analyzed once both its bonds appear in a later frame."""
    traj = tmp_path / "a.xyz"
    frames, symbols = bent_frames([2.0, 1.4, 1.4])  # frame 1: no triplet
    write_xyz(traj, frames, symbols)
    a = AllAnglesAnalysis(str(traj), max_connection_distance=1.5)
    stats = a.analyze()
    assert len(a.angle_triplets) == 1
    assert a.angle_triplets[0].frames_present == 2
    assert stats[0][0] == pytest.approx(90.0)


def test_no_angle_triplets_raises(tmp_path):
    traj = tmp_path / "a.xyz"
    frames, symbols = bent_frames([2.0])  # O-H2 too long -> only one bond
    write_xyz(traj, frames, symbols)
    a = AllAnglesAnalysis(str(traj), max_connection_distance=1.5)
    with pytest.raises(ValueError, match="No connected angle triplets"):
        a.analyze()


def test_solute_scope_restricts_triplets(tmp_path):
    """A nearby solvent atom bonded to the center is excluded by solute scope."""
    traj = tmp_path / "a.xyz"
    coords = np.array([[1.0, 0, 0], [0.0, 0, 0], [0.0, 1.4, 0], [0.0, 0.0, 1.0]])
    write_xyz(traj, [coords], ["H", "O", "H", "X"])  # atom 3 (X) also bonds O
    a = AllAnglesAnalysis(str(traj), max_connection_distance=1.5,
                          solute_atom_indices=[0, 1, 2])
    a.analyze()
    triplets = [(t.atom_i, t.atom_j, t.atom_k) for t in a.angle_triplets]
    assert triplets == [(0, 1, 2)]


# --------------------------------------------------------------------------- #
# PBC (optional orthorhombic minimum-image)                                     #
# --------------------------------------------------------------------------- #
def test_pbc_forms_angle_across_boundary(tmp_path):
    traj = tmp_path / "pbc.xyz"
    # O near x=0 face; H_i across the boundary (min-image O-H_i = 1.4 A); H_k in +y.
    coords = np.array([[8.8, 0, 0], [0.2, 0, 0], [0.2, 1.4, 0]])
    write_xyz(traj, [coords], ["H", "O", "H"])

    with_pbc = AllAnglesAnalysis(str(traj), max_connection_distance=1.5,
                                 cell_lengths=[10.0, 10.0, 10.0])
    stats = with_pbc.analyze()
    assert stats[0][0] == pytest.approx(90.0)
    assert with_pbc.angle_triplets[0].frames_present == 1

    without = AllAnglesAnalysis(str(traj), max_connection_distance=1.5)
    with pytest.raises(ValueError, match="No connected angle triplets"):
        without.analyze()


def test_pbc_box_shift_invariance(tmp_path):
    base = tmp_path / "base.xyz"
    shifted = tmp_path / "shifted.xyz"
    frames, symbols = bent_frames([1.4, 1.4])
    write_xyz(base, frames, symbols)
    write_xyz(shifted, [c + np.array([10.0, 0, 0]) for c in frames], symbols)
    box = [10.0, 10.0, 10.0]
    a = AllAnglesAnalysis(str(base), max_connection_distance=1.5, cell_lengths=box)
    b = AllAnglesAnalysis(str(shifted), max_connection_distance=1.5, cell_lengths=box)
    assert a.analyze()[0] == pytest.approx(b.analyze()[0])


def test_pbc_rejects_cutoff_above_half_box(tmp_path):
    traj = tmp_path / "pbc.xyz"
    frames, symbols = bent_frames([1.4])
    write_xyz(traj, frames, symbols)
    a = AllAnglesAnalysis(str(traj), max_connection_distance=2.0,
                          cell_lengths=[3.0, 3.0, 3.0])
    with pytest.raises(ValueError, match="half"):
        a.analyze()


# --------------------------------------------------------------------------- #
# Async parity + combined output                                                #
# --------------------------------------------------------------------------- #
def test_async_matches_sync(tmp_path):
    traj = tmp_path / "a.xyz"
    frames, symbols = bent_frames([1.4, 1.4, 1.4])
    write_xyz(traj, frames, symbols)
    sync = AllAnglesAnalysis(str(traj), max_connection_distance=1.5).analyze()
    a = AllAnglesAnalysis(str(traj), max_connection_distance=1.5)
    async_stats = asyncio.run(a.analyze_async())
    assert async_stats == pytest.approx(sync)


def test_write_results_single_combined_file(tmp_path):
    traj = tmp_path / "a.xyz"
    frames, symbols = bent_frames([1.4, 1.4])
    write_xyz(traj, frames, symbols)
    a = AllAnglesAnalysis(str(traj), max_connection_distance=1.5,
                          cell_lengths=[10.0, 10.0, 10.0])
    a.analyze()
    output_file = a.write_results(str(tmp_path / "angles.txt"))

    assert isinstance(output_file, str)
    assert not os.path.exists(str(tmp_path / "angles_triplets.txt"))
    text = open(output_file).read()
    assert "# gQTEA All Bond Angle Analysis" in text
    assert "# frames_used 2" in text
    assert "# max_connection_distance 1.5" in text
    assert "# cell_lengths 10 10 10" in text
    header = [l for l in text.splitlines() if l.startswith("# row")][0]
    for col in ("atom_i", "first_present_angle_degrees", "average_angle_degrees",
                "occurrence_fraction", "frames_present"):
        assert col in header
    rows = [l for l in text.splitlines() if l and not l.startswith("#")]
    assert len(rows) == 1
    fields = rows[0].split()
    assert fields[1:7] == ["1", "2", "3", "H", "O", "H"]
    assert float(fields[8]) == pytest.approx(90.0)   # average angle
    assert float(fields[11]) == pytest.approx(1.0)   # occurrence
    assert int(fields[12]) == 2                       # frames_present


# --------------------------------------------------------------------------- #
# Validation                                                                    #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("cell", [[10.0, 10.0], [10.0, 0.0, 10.0]])
def test_invalid_cell_raises(tmp_path, cell):
    traj = tmp_path / "a.xyz"
    frames, symbols = bent_frames([1.4])
    write_xyz(traj, frames, symbols)
    with pytest.raises(ValueError):
        AllAnglesAnalysis(str(traj), cell_lengths=cell).analyze()


def test_ranges_supported_via_shared_parser(tmp_path):
    """Solute range parsing is the same shared helper used by the bond tool."""
    from allAnglesAnalysis import parse_solute_index_ranges
    assert parse_solute_index_ranges("1-3 5") == [1, 2, 3, 5]
