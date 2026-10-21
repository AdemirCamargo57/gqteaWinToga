"""3D viewer: loading an .xyz trajectory (read_xyz_frames + load_xyz_trajectory).

The bulk reader must give exactly what the line-by-line reader it replaced
gave -- the same frames (lists of (element, (x, y, z)) tuples), the same byte
offsets and the same stopping point on a malformed file -- so that reader is
kept here as the reference and every case is compared with it. The window
half (spinner, Browse disabled) is driven headless through the hooks.
"""
import asyncio
import gc
import threading
from pathlib import Path

import pytest

import molecularViewer as mv
from molecularViewer import MolecularViewer, read_xyz_frames

FIXTURE = Path(__file__).parent / "data" / "mrt_classical" / "na_water.xyz"


@pytest.fixture(autouse=True)
def nothing_stays_frozen():
    """A load freezes the collector's heap; never let that leak into other tests."""
    yield
    gc.unfreeze()


def reference_reader(path):
    """load_xyz_trajectory's loop before the bulk reader (frames, offsets)."""
    frames, offsets = [], []
    with open(path, "r", encoding="utf-8") as f:
        while True:
            frame_offset = f.tell()
            header = f.readline()
            if not header:
                break
            header = header.strip()
            if not header:
                continue
            try:
                num_atoms = int(header)
            except ValueError:
                continue
            comment = f.readline()
            if comment == "":
                break
            frame_data, malformed = [], False
            for _ in range(num_atoms):
                atom_line = f.readline()
                if atom_line == "":
                    malformed = True
                    break
                parts = atom_line.split()
                if len(parts) < 4:
                    malformed = True
                    break
                try:
                    x, y, z = map(float, parts[1:4])
                except ValueError:
                    malformed = True
                    break
                frame_data.append((parts[0], (x, y, z)))
            if malformed:
                break
            if frame_data:
                frames.append(frame_data)
                offsets.append(frame_offset)
    return frames, offsets


def frame(atoms, comment="c"):
    return [str(len(atoms)), comment] + atoms


WATER = ["O 0.0 0.0 0.1173", "H 0.0 0.7572 -0.4692", "H 0.0 -0.7572 -0.4692"]
MOVED = ["O 0.1 0.0 0.1173", "H 0.1 0.7572 -0.4692", "H 0.1 -0.7572 -0.4692"]

CASES = {
    "plain": frame(WATER) + frame(MOVED) + frame(WATER),
    "blank and junk lines between frames": ["", "junk"] + frame(WATER) + ["", "  "] + frame(MOVED) + ["???"],
    "extra columns": frame([a + " 0.5 1" for a in WATER]) + frame(MOVED),
    "atom count changes": frame(WATER) + frame(WATER[:2]) + frame(["C 1 2 3"] + MOVED),
    "element order changes": frame(WATER) + frame(["H 0 0 0", "O 1 1 1", "H 2 2 2"]) + frame(WATER),
    "zero and negative atom counts": frame(WATER) + ["0", "empty"] + ["-2", "negative"] + frame(MOVED),
    "truncated last frame": frame(WATER) + frame(MOVED)[:-1],
    "missing comment at the end": frame(WATER) + ["3"],
    "short line mid-file": frame(WATER) + frame(["O 0 0", "H 1 1 1", "H 2 2 2"]) + frame(WATER),
    "bad number mid-file": frame(WATER) + frame(["O 0 0 x", "H 1 1 1", "H 2 2 2"]) + frame(WATER),
    # 4 x atoms fields in all, but one line short and one long: must not be read in bulk
    "short and long line in one frame": frame(WATER) + frame(["O 0 0", "H 1 1 1 1", "H 2 2 2"]) + frame(WATER),
    "scientific notation and signs": frame(["O 1e-3 -2.5E+00 +3", "H .5 -.25 1.", "H 0 0 0"]) * 2,
    "numeric element labels": frame(["8 0 0 0.1173", "1 0 0.7572 -0.4692", "1 0 -0.7572 -0.4692"]) * 3,
    "non-ascii comment": frame(WATER, "água — café") + frame(MOVED, "ångström"),
}


def write(path, lines, end="\n", final=True):
    path.write_bytes((end.join(lines) + (end if final else "")).encode("utf-8"))
    return path


@pytest.mark.parametrize("case", sorted(CASES))
@pytest.mark.parametrize("end", ["\n", "\r\n"])
def test_same_frames_and_offsets_as_the_reference_reader(tmp_path, case, end):
    path = write(tmp_path / "t.xyz", CASES[case], end)
    assert read_xyz_frames(str(path)) == reference_reader(path)


@pytest.mark.parametrize("chunk", [1, 2, 3, 4, 5, 7, 11])
def test_block_boundaries_do_not_change_the_result(tmp_path, monkeypatch, chunk):
    monkeypatch.setattr(mv, "XYZ_READ_CHUNK_LINES", chunk)
    for case, lines in CASES.items():
        path = write(tmp_path / "t.xyz", lines)
        assert read_xyz_frames(str(path)) == reference_reader(path), case


def test_no_final_newline(tmp_path):
    path = write(tmp_path / "t.xyz", CASES["plain"], final=False)
    frames, offsets = read_xyz_frames(str(path))
    assert len(frames) == 3 and (frames, offsets) == reference_reader(path)


def test_the_31_atom_300_frame_fixture():
    frames, offsets = read_xyz_frames(str(FIXTURE))
    assert len(frames) == 300 and all(len(f) == 31 for f in frames)
    assert (frames, offsets) == reference_reader(FIXTURE)
    element, (x, y, z) = frames[0][0]
    assert element == "Na" and all(type(v) is float for v in (x, y, z))


def test_empty_file(tmp_path):
    path = write(tmp_path / "t.xyz", [], final=False)
    assert read_xyz_frames(str(path)) == ([], [])


def test_progress_reaches_the_end(tmp_path, monkeypatch):
    monkeypatch.setattr(mv, "XYZ_READ_CHUNK_LINES", 4)
    seen = []
    read_xyz_frames(str(write(tmp_path / "t.xyz", CASES["plain"])), lambda n, f: seen.append((n, f)))
    assert len(seen) > 2 and seen[-1] == (3, 1.0)
    assert [f for _, f in seen] == sorted(f for _, f in seen)


def test_invalid_utf8_is_an_error_as_before(tmp_path):
    path = tmp_path / "t.xyz"
    path.write_bytes(b"1\nc\nO 0 0 \xff\n")
    with pytest.raises(UnicodeDecodeError):
        reference_reader(path)
    with pytest.raises(UnicodeDecodeError):
        read_xyz_frames(str(path))


@pytest.mark.parametrize("enabled", [True, False])
def test_the_garbage_collector_is_restored(tmp_path, enabled):
    path = write(tmp_path / "t.xyz", CASES["plain"])
    states = []
    try:
        gc.enable() if enabled else gc.disable()
        read_xyz_frames(str(path), lambda n, f: states.append(gc.isenabled()))
        assert gc.isenabled() is enabled
        assert len(states) >= 2 and not any(states[:-1])    # paused while reading
        assert states[-1] is enabled                         # the closing report comes after
    finally:
        gc.enable()


def test_freeze_moves_the_frames_out_of_the_collector(tmp_path):
    path = write(tmp_path / "t.xyz", CASES["plain"])
    gc.unfreeze()
    read_xyz_frames(str(path))
    assert gc.get_freeze_count() == 0                  # only on request
    frames, _ = read_xyz_frames(str(path), freeze=True)
    assert gc.get_freeze_count() > 0 and gc.isenabled()
    assert gc.is_tracked(frames[0][0])                 # frozen, not untracked: freed as usual


# ------------------------------------------------------------------
# load_xyz_trajectory
# ------------------------------------------------------------------
class RecordingViewer(MolecularViewer):
    def __init__(self):
        super().__init__()
        self.events, self.threads = [], []

    def _loading_started(self):
        self.events.append("started")

    def _loading_finished(self):
        self.events.append("finished")


def test_load_installs_the_frames(tmp_path):
    viewer = RecordingViewer()
    path = write(tmp_path / "t.xyz", CASES["plain"])
    asyncio.run(viewer.load_xyz_trajectory(str(path)))
    frames, offsets = reference_reader(path)
    assert viewer.frames == frames and viewer.frame_offsets == offsets
    assert viewer.molecule_data is viewer.frames[0] and viewer.current_frame == 0
    assert viewer.frame_elements == ["O", "H", "H"]
    assert viewer.bonds                             # frame 0's bonds were detected
    assert viewer.events == ["started", "finished"]
    assert viewer._status_message == "3 frames loaded."


def test_reading_runs_in_a_worker_thread(tmp_path, monkeypatch):
    viewer = RecordingViewer()
    threads = []
    original = mv.read_xyz_frames

    def reader(*args):
        threads.append(threading.get_ident())
        return original(*args)
    monkeypatch.setattr(mv, "read_xyz_frames", reader)
    asyncio.run(viewer.load_xyz_trajectory(str(write(tmp_path / "t.xyz", CASES["plain"]))))
    assert threads and threads[0] != threading.get_ident()


def test_loading_status_is_reported(tmp_path, monkeypatch):
    monkeypatch.setattr(MolecularViewer, "LOADING_STATUS_INTERVAL", 0.0)
    monkeypatch.setattr(mv, "XYZ_READ_CHUNK_LINES", 5)
    viewer = RecordingViewer()
    messages = []
    viewer.set_status_message = messages.append
    asyncio.run(viewer.load_xyz_trajectory(str(write(tmp_path / "t.xyz", CASES["plain"]))))
    assert messages[0] == "Loading trajectory..."
    assert any(m.startswith("Loading trajectory... ") and "frames (" in m for m in messages)
    assert messages[-1] == "3 frames loaded."


def test_no_frames(tmp_path):
    viewer = RecordingViewer()
    asyncio.run(viewer.load_xyz_trajectory(str(write(tmp_path / "t.xyz", ["junk", "more junk"]))))
    assert viewer.frames == [] and viewer.molecule_data == [] and viewer.frame_elements == []
    assert viewer._status_message == "No valid frames found in XYZ file."
    assert viewer.events == ["started", "finished"]


def test_a_read_error_still_ends_the_loading_state(tmp_path):
    viewer = RecordingViewer()
    with pytest.raises(OSError):
        asyncio.run(viewer.load_xyz_trajectory(str(tmp_path / "missing.xyz")))
    assert viewer.events == ["started", "finished"]
    assert viewer._status_message.startswith("I/O error while reading file")


def test_a_load_freezes_and_the_next_load_unfreezes(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(mv.gc, "unfreeze", lambda: calls.append("unfreeze"))
    original = mv.read_xyz_frames

    def reader(path, progress=None, freeze=False):
        calls.append(("read", freeze))
        return original(path, progress)               # (no real freeze inside the test)
    monkeypatch.setattr(mv, "read_xyz_frames", reader)
    viewer = RecordingViewer()
    for name in ("a.xyz", "b.xyz"):
        asyncio.run(viewer.load_xyz_trajectory(str(write(tmp_path / name, CASES["plain"]))))
    assert calls == ["unfreeze", ("read", True), "unfreeze", ("read", True)]


def test_closing_the_window_unfreezes(monkeypatch):
    calls = []
    monkeypatch.setattr(mv.gc, "unfreeze", lambda: calls.append("unfreeze"))
    ui = mv.MolecularViewerUI.__new__(mv.MolecularViewerUI)
    ui.design_controller = None
    assert ui._on_main_window_close(None) is True
    assert calls == ["unfreeze"]


def test_a_new_load_replaces_the_previous_trajectory(tmp_path):
    viewer = RecordingViewer()
    asyncio.run(viewer.load_xyz_trajectory(str(write(tmp_path / "a.xyz", CASES["plain"]))))
    asyncio.run(viewer.load_xyz_trajectory(str(write(tmp_path / "b.xyz", frame(["C 0 0 0"])))))
    assert viewer.frames == [[("C", (0.0, 0.0, 0.0))]] and viewer.frame_elements == ["C"]


# ------------------------------------------------------------------
# The window's hooks: spinner and Browse
# ------------------------------------------------------------------
class StubSpinner:
    def __init__(self):
        self.events = []

    def start(self):
        self.events.append("start")

    def stop(self):
        self.events.append("stop")


def test_the_window_spins_and_blocks_browse_while_loading(tmp_path):
    from types import SimpleNamespace
    ui = mv.MolecularViewerUI.__new__(mv.MolecularViewerUI)
    MolecularViewer.__init__(ui)
    ui.loading_indicator, ui.browse_button = StubSpinner(), SimpleNamespace(enabled=True)
    seen = []
    original = mv.read_xyz_frames

    def reader(*args):
        seen.append((list(ui.loading_indicator.events), ui.browse_button.enabled))
        return original(*args)
    mv_reader, mv.read_xyz_frames = mv.read_xyz_frames, reader
    try:
        ui.update_frame_label = lambda: None
        asyncio.run(ui.load_xyz_trajectory(str(write(tmp_path / "t.xyz", CASES["plain"]))))
    finally:
        mv.read_xyz_frames = mv_reader
    assert seen == [(["start"], False)]
    assert ui.loading_indicator.events == ["start", "stop"] and ui.browse_button.enabled
