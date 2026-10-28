"""runatom.exe against the Linux reference build (spec §7). Skipped without runatom.exe.

Rule: a disagreement is investigated in the build (flags, compiler), never fixed
by loosening a tolerance.
"""

import shutil
import time

import numpy as np
import pytest

import runatomInput as ri
import runatomOutput as ro
import runatomRunner as rr
from runatom_fixtures import LIBRARY_DIR, REFERENCE_CASES, case_info, library_case_dir, reference_bytes, reference_names, reference_text

EXE = rr.find_runatom()
pytestmark = pytest.mark.skipif(EXE is None, reason="runatom.exe is not available")

EIGEN_TOL = 2e-6          # Ry: one unit in the last printed place
BESSEL_TOL = 1e-5         # Ry
MESH_RTOL = 1e-12
POT_RTOL, POT_ATOL = 1e-7, 1e-10
PROJ_RTOL = 1e-6
# Source/BENIGN_BUG: the 2x2 s-channel rotation (ttt) is supersensitive. The Pb s channel is the
# documented case; the other s channels were measured to be the same effect: projectors differ by
# more than a sign while the eigenvalues of the D and Q blocks agree to < 1e-9 relative.
UNSTABLE = {"082-Pb-ca-d-vgrp": {0}, "000-Ex-ca-sp-test": {0}, "008-O-ca--vgrp": {0},
            "029-Cu-gpbe-n-bm": {0}, "056-Ba-gpw-nsp-campos": {0}, "079-Au-ca--ts": {0}}


def _npf_is_defined(iptype) -> bool:
    """npf and ptryc are read only for iptype = 2 (rwsubs.f); otherwise runatom writes them uninitialised."""
    return bool(np.any(np.asarray(iptype) == 2))


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    """Every reference case run once with runatom.exe on the original library inputs."""
    done = {}
    for name in REFERENCE_CASES:
        work = tmp_path_factory.mktemp(name)
        shutil.copytree(library_case_dir(name), work, dirs_exist_ok=True)
        case = ri.read_case_folder(work)
        started = time.monotonic()
        rr.run_job(rr.runatom_command(EXE), rr.plan_stages(case), work)
        done[name] = (work, case, time.monotonic() - started)
    return done


@pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
def test_eigenvalue_tables_match_the_reference(runs, name):
    work = runs[name][0]
    for filename in (n for n in reference_names(name) if n.endswith(".out")):
        ours = ro.parse_out((work / filename).read_text(encoding="latin-1"))
        ref = ro.parse_out(reference_text(name, filename))
        assert [r.nlm for r in ours.comparison] == [r.nlm for r in ref.comparison], filename
        for a, b in zip(ours.comparison, ref.comparison):
            assert abs(a.all_electron - b.all_electron) <= EIGEN_TOL and abs(a.pseudo - b.pseudo) <= EIGEN_TOL, filename
        assert sorted(ours.bessel) == sorted(ref.bessel), filename
        for l in ref.bessel:
            assert np.allclose(ours.bessel[l], ref.bessel[l], rtol=0, atol=BESSEL_TOL), (filename, l)
        if ref.total_energy is not None:
            assert abs(ours.total_energy - ref.total_energy) <= EIGEN_TOL, filename


def _signs(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.where(np.einsum("ij,ij->i", a, b) < 0, -1.0, 1.0)


def _block_eigenvalues(matrix: np.ndarray, l: np.ndarray, channel: int) -> np.ndarray:
    index = np.flatnonzero(l == channel)
    return np.linalg.eigvalsh(matrix[np.ix_(index, index)])


@pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
def test_pseudopotential_file_matches_the_reference(runs, name, tmp_path):
    atom = case_info(name)["atom"]
    ours = ro.read_uspp(runs[name][0] / f"{atom}_ps.uspp")
    (tmp_path / "ref.uspp").write_bytes(reference_bytes(name, f"{atom}_ps.uspp"))
    ref = ro.read_uspp(tmp_path / "ref.uspp")
    for attr in ("version", "title", "z", "zv", "exfact", "mesh", "keyps", "ifpcor", "nang", "lloc",
                 "ifqopt", "nqf", "irel", "kkbeta") + (("npf",) if _npf_is_defined(ref.iptype) else ()):
        assert getattr(ours, attr) == getattr(ref, attr), attr
    assert np.array_equal(ours.l, ref.l) and np.array_equal(ours.iptype, ref.iptype)
    assert np.allclose(ours.r, ref.r, rtol=MESH_RTOL, atol=0) and np.allclose(ours.rab, ref.rab, rtol=MESH_RTOL, atol=0)
    for attr in ("vloc", "vloc0", "rsatom"):
        assert np.allclose(getattr(ours, attr), getattr(ref, attr), rtol=POT_RTOL, atol=POT_ATOL), attr
    unstable = UNSTABLE.get(name, set())
    stable = np.array([l not in unstable for l in ref.l])
    s = _signs(ours.beta, ref.beta)
    scale = np.abs(ref.beta).max()
    assert np.allclose((s[:, None] * ours.beta)[stable], ref.beta[stable], rtol=PROJ_RTOL, atol=PROJ_RTOL * scale)
    pair = np.outer(stable, stable)
    for attr in ("ddd", "qqq"):
        mine = np.outer(s, s) * getattr(ours, attr)
        assert np.allclose(mine[pair], getattr(ref, attr)[pair], rtol=PROJ_RTOL, atol=1e-9), attr
    for channel in unstable:      # representation-independent check (BENIGN_BUG)
        for attr in ("ddd", "qqq"):
            assert np.allclose(_block_eigenvalues(getattr(ours, attr), ours.l, channel),
                               _block_eigenvalues(getattr(ref, attr), ref.l, channel), rtol=1e-6, atol=1e-9)


@pytest.mark.parametrize("name", sorted(REFERENCE_CASES))
def test_reformatted_inputs_give_an_identical_pseudopotential(runs, name, tmp_path):
    work, case, _ = runs[name]
    ri.write_case(ri.replace(case, tests=[]), tmp_path)        # our writer, same values
    rr.run_job(rr.runatom_command(EXE), rr.plan_stages(ri.replace(case, tests=[])), tmp_path)
    stem = ri.file_stems(case)["ps"]
    original = ro.uspp_records((work / f"{stem}.uspp").read_bytes())
    rewritten = ro.uspp_records((tmp_path / f"{stem}.uspp").read_bytes())
    assert len(original) == len(rewritten)
    data = ro.read_uspp(work / f"{stem}.uspp")
    # The iptype record (iptype(1:nbeta), npf, ptryc) comes right before the (rcloc, vloc0) record,
    # the one of 8*(mesh+1) bytes. npf and ptryc are uninitialised unless a projector has iptype 2.
    masked = None
    if not _npf_is_defined(data.iptype):
        masked = next(k for k, record in enumerate(original) if len(record) == 8 * (data.mesh + 1)) - 1
    for index, (a, b) in enumerate(zip(original, rewritten)):
        if index == 0:
            continue                              # record 1 holds the run date
        if index == masked:
            a, b = a[:-12], b[:-12]
        assert a == b, f"record {index + 1}"


def test_sample_output_of_2006_is_reproduced(runs):
    work = runs["000-Ex-ca-sp-test"][0]
    sample = LIBRARY_DIR / "000-Ex" / "Sample_output"
    for filename in ("ti_ps.out", "ti_test_s2d2.out"):
        ours = ro.parse_out((work / filename).read_text(encoding="latin-1")).comparison
        old = ro.parse_out((sample / filename).read_text(encoding="latin-1")).comparison
        assert [r.nlm for r in ours] == [r.nlm for r in old]
        for a, b in zip(ours, old):
            assert abs(a.all_electron - b.all_electron) <= EIGEN_TOL and abs(a.pseudo - b.pseudo) <= EIGEN_TOL


def test_unreadable_library_input_is_reported_as_an_input_problem(tmp_path):
    shutil.copytree(LIBRARY_DIR / "023-V" / "023-V-gpw-nsp-campos", tmp_path, dirs_exist_ok=True)
    stage = rr.Stage("pseudo", "V generation", "V_ps", "V_ae_ref.ae", "V_ps.uspp")
    rr.run_stage(rr.runatom_command(EXE), rr.Stage("ae", "V AE", "V_ae_ref", "V_ae_ref.ae", "dummy"), tmp_path)
    with pytest.raises(rr.RunatomError) as info:
        rr.run_stage(rr.runatom_command(EXE), stage, tmp_path)
    assert info.value.kind == "input_format"


@pytest.mark.parametrize("folder", ["020-Ca/020-Ca-gpw-nsp-campos", "078-Pt/078-Pt-gpw-n-campos",
                                    "029-Cu/029-Cu-gpw-n-campos"])
def test_imported_library_case_runs_from_the_window_fields(folder, tmp_path):
    # These cases tune damp, maxit, ee or the test input per test; the window's fields hold only
    # name, xion and occupations, so the imported inputs must be attached back to run at all.
    imported = ri.read_case_folder(LIBRARY_DIR / folder)
    case = ri.attach_imported_tests(ri.case_from_fields(ri.case_to_fields(imported)), imported)
    assert ri.check_case(case) == []
    ri.write_case(case, tmp_path)
    stages = rr.plan_stages(case)
    results = rr.run_job(rr.runatom_command(EXE), stages, tmp_path)
    assert len(results) == len(stages) and stages[-1].kind == "test"


def test_every_reference_case_is_fast(runs):
    assert max(seconds for _, _, seconds in runs.values()) < 10.0
