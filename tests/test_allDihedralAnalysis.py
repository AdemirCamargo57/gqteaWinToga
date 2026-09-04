"""Tests for the All Dihedral Angle Analysis tool (``allDihedralAnalysis.py``).

Run from ``venv/src/``. Only the logic class ``AllDihedralAnalysis`` is exercised
(no Toga GUI). Mirrors the all-bond / all-angle tools: connectivity is
re-evaluated every frame, each i-j-k-l dihedral is averaged only over the frames
where all three of its bonds (i-j, j-k, k-l) are within the cutoff, occurrence is
reported, PBC is optional (orthorhombic minimum-image), the connection distance
defaults to 1.7, solute indices accept ranges, and results go to one combined file.
"""
import os
import asyncio

import numpy as np
import pytest

from allDihedralAnalysis import AllDihedralAnalysis, DihedralQuadruplet


def write_xyz(path, frames, symbols, comment="frame"):
    with open(path, "w") as f:
        for coords in frames:
            f.write(f"{len(symbols)}\n{comment}\n")
            for s, (x, y, z) in zip(symbols, coords):
                f.write(f"{s} {x:.6f} {y:.6f} {z:.6f}\n")


def planar_chain(kl=1.2):
    """Planar i-j-k-l chain (trans, dihedral = 180 deg). Non-adjacent atoms stay
    > 1.5 A apart so no spurious bond closes a ring. ``kl`` sets the k-l distance
    (set > 1.5 to break the terminal bond)."""
    i = [0.0, 0.0, 0.0]
    j = [1.2, 0.0, 0.0]
    k = [1.2, 1.2, 0.0]
    l = [1.2 + kl, 1.2, 0.0]
    return np.array([i, j, k, l]), ["C", "C", "C", "C"]


# --------------------------------------------------------------------------- #
# Metadata + default                                                            #
# --------------------------------------------------------------------------- #
def test_read_first_frame_sets_metadata(tmp_path):
    traj = tmp_path / "t.xyz"
    coords, symbols = planar_chain()
    write_xyz(traj, [coords], symbols)
    a = AllDihedralAnalysis(str(traj))
    elements, xyz = a.read_first_frame()
    assert elements == ["C", "C", "C", "C"]
    assert a.num_atoms == 4


def test_default_connection_distance_is_1_7():
    assert AllDihedralAnalysis.DEFAULT_CONNECTION_DISTANCE == pytest.approx(1.7)
    assert AllDihedralAnalysis("x.xyz").max_connection_distance == pytest.approx(1.7)


# --------------------------------------------------------------------------- #
# Per-frame occurrence and the frozen-frame fix                                 #
# --------------------------------------------------------------------------- #
def test_dihedral_present_every_frame(tmp_path):
    traj = tmp_path / "d.xyz"
    coords, symbols = planar_chain()
    write_xyz(traj, [coords, coords], symbols)
    a = AllDihedralAnalysis(str(traj), max_connection_distance=1.5)
    stats = a.analyze()
    assert len(stats) == 1
    mean, variance, std = stats[0]
    assert abs(mean) == pytest.approx(180.0)
    assert variance == pytest.approx(0.0)
    quad = a.dihedral_quadruplets[0]
    assert (quad.atom_i, quad.atom_j, quad.atom_k, quad.atom_l) == (0, 1, 2, 3)
    assert quad.frames_present == 2
    assert a.num_frames == 2


def test_dihedral_only_over_present_frames(tmp_path):
    """Frame 2 breaks the k-l bond, so that frame does not contribute."""
    traj = tmp_path / "d.xyz"
    present, symbols = planar_chain(kl=1.2)
    broken, _ = planar_chain(kl=2.0)   # k-l = 2.0 > 1.5 -> quadruplet absent
    write_xyz(traj, [present, broken, present], symbols)
    a = AllDihedralAnalysis(str(traj), max_connection_distance=1.5)
    a.analyze()
    assert a.dihedral_quadruplets[0].frames_present == 2
    assert a.num_frames == 3


def test_dihedral_forming_later_is_captured(tmp_path):
    """Regression for the frozen-frame bug: a quadruplet absent in frame 1 is
    still analyzed once all its bonds appear in a later frame."""
    traj = tmp_path / "d.xyz"
    present, symbols = planar_chain(kl=1.2)
    broken, _ = planar_chain(kl=2.0)
    write_xyz(traj, [broken, present, present], symbols)
    a = AllDihedralAnalysis(str(traj), max_connection_distance=1.5)
    stats = a.analyze()
    assert len(a.dihedral_quadruplets) == 1
    assert a.dihedral_quadruplets[0].frames_present == 2
    assert abs(stats[0][0]) == pytest.approx(180.0)


def test_no_dihedral_quadruplets_raises(tmp_path):
    traj = tmp_path / "d.xyz"
    broken, symbols = planar_chain(kl=2.0)
    write_xyz(traj, [broken], symbols)
    a = AllDihedralAnalysis(str(traj), max_connection_distance=1.5)
    with pytest.raises(ValueError, match="No connected dihedral quadruplets"):
        a.analyze()


def test_solute_scope_restricts_quadruplets(tmp_path):
    """A nearby solvent atom is excluded by the solute scope."""
    traj = tmp_path / "d.xyz"
    coords, symbols = planar_chain(kl=1.2)
    coords = np.vstack([coords, [1.2, 0.0, 1.0]])   # solvent atom bonded to j
    write_xyz(traj, [coords], symbols + ["X"])
    a = AllDihedralAnalysis(str(traj), max_connection_distance=1.5,
                            solute_atom_indices=[0, 1, 2, 3])
    a.analyze()
    quads = [(q.atom_i, q.atom_j, q.atom_k, q.atom_l) for q in a.dihedral_quadruplets]
    assert quads == [(0, 1, 2, 3)]


# --------------------------------------------------------------------------- #
# PBC (optional orthorhombic minimum-image)                                     #
# --------------------------------------------------------------------------- #
def test_pbc_box_shift_invariance(tmp_path):
    base = tmp_path / "base.xyz"
    shifted = tmp_path / "shifted.xyz"
    coords, symbols = planar_chain()
    write_xyz(base, [coords, coords], symbols)
    write_xyz(shifted, [coords + np.array([10.0, 0, 0]), coords + np.array([10.0, 0, 0])], symbols)
    box = [10.0, 10.0, 10.0]
    a = AllDihedralAnalysis(str(base), max_connection_distance=1.5, cell_lengths=box)
    b = AllDihedralAnalysis(str(shifted), max_connection_distance=1.5, cell_lengths=box)
    assert a.analyze()[0] == pytest.approx(b.analyze()[0])


def test_pbc_reconnects_chain_across_boundary(tmp_path):
    """A chain whose i-j bond spans the periodic boundary is only fully connected
    (and thus yields a dihedral) when the cell is supplied."""
    traj = tmp_path / "pbc.xyz"
    # Shift atom i one full box in +x so the i-j bond wraps around.
    coords, symbols = planar_chain(kl=1.2)
    coords[0] = coords[0] + np.array([10.0, 0.0, 0.0])  # i now at x=10.0
    write_xyz(traj, [coords], symbols)

    with_pbc = AllDihedralAnalysis(str(traj), max_connection_distance=1.5,
                                   cell_lengths=[10.0, 10.0, 10.0])
    stats = with_pbc.analyze()
    assert abs(stats[0][0]) == pytest.approx(180.0)

    without = AllDihedralAnalysis(str(traj), max_connection_distance=1.5)
    with pytest.raises(ValueError, match="No connected dihedral quadruplets"):
        without.analyze()


def test_pbc_rejects_cutoff_above_half_box(tmp_path):
    traj = tmp_path / "pbc.xyz"
    coords, symbols = planar_chain()
    write_xyz(traj, [coords], symbols)
    a = AllDihedralAnalysis(str(traj), max_connection_distance=2.0,
                            cell_lengths=[3.0, 3.0, 3.0])
    with pytest.raises(ValueError, match="half"):
        a.analyze()


# --------------------------------------------------------------------------- #
# Async parity + combined output                                                #
# --------------------------------------------------------------------------- #
def test_async_matches_sync(tmp_path):
    traj = tmp_path / "d.xyz"
    coords, symbols = planar_chain()
    write_xyz(traj, [coords, coords, coords], symbols)
    sync = AllDihedralAnalysis(str(traj), max_connection_distance=1.5).analyze()
    a = AllDihedralAnalysis(str(traj), max_connection_distance=1.5)
    async_stats = asyncio.run(a.analyze_async())
    assert async_stats == pytest.approx(sync)


def test_write_results_single_combined_file(tmp_path):
    traj = tmp_path / "d.xyz"
    coords, symbols = planar_chain()
    write_xyz(traj, [coords, coords], symbols)
    a = AllDihedralAnalysis(str(traj), max_connection_distance=1.5,
                            cell_lengths=[10.0, 10.0, 10.0])
    a.analyze()
    output_file = a.write_results(str(tmp_path / "dih.txt"))

    assert isinstance(output_file, str)
    assert not os.path.exists(str(tmp_path / "dih_quadruplets.txt"))
    text = open(output_file).read()
    assert "# gQTEA All Dihedral Angle Analysis" in text
    assert "# frames_used 2" in text
    assert "# max_connection_distance 1.5" in text
    assert "# dihedral_convention signed_degrees_minus180_to_180" in text
    assert "# cell_lengths 10 10 10" in text
    header = [l for l in text.splitlines() if l.startswith("# row")][0]
    for col in ("atom_i", "first_present_dihedral_degrees", "average_dihedral_degrees",
                "occurrence_fraction", "frames_present"):
        assert col in header
    rows = [l for l in text.splitlines() if l and not l.startswith("#")]
    assert len(rows) == 1
    fields = rows[0].split()
    assert fields[1:9] == ["1", "2", "3", "4", "C", "C", "C", "C"]  # atoms + elements
    assert abs(float(fields[9])) == pytest.approx(180.0)   # first_present_dihedral
    assert abs(float(fields[10])) == pytest.approx(180.0)  # average dihedral
    assert float(fields[13]) == pytest.approx(1.0)         # occurrence_fraction
    assert int(fields[14]) == 2                             # frames_present


# --------------------------------------------------------------------------- #
# Validation                                                                    #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("cell", [[10.0, 10.0], [10.0, 0.0, 10.0]])
def test_invalid_cell_raises(tmp_path, cell):
    traj = tmp_path / "d.xyz"
    coords, symbols = planar_chain()
    write_xyz(traj, [coords], symbols)
    with pytest.raises(ValueError):
        AllDihedralAnalysis(str(traj), cell_lengths=cell).analyze()


def test_ranges_supported_via_shared_parser():
    from allDihedralAnalysis import parse_solute_index_ranges
    assert parse_solute_index_ranges("1-4 6") == [1, 2, 3, 4, 6]
