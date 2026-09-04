"""Pytest bootstrap for the gqteaWinToga source tree.

The tools use flat, bare-name imports (``from framesCounter import ...``), so the
source directory (the parent of this ``tests/`` folder) must be on ``sys.path``.
We also pin matplotlib to the non-interactive Agg backend before anything imports
pyplot, so tests never try to open a GUI window.
"""
import os
import sys

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)


@pytest.fixture(autouse=True)
def _force_agg_backend():
    """Keep every test on the non-interactive Agg backend.

    ``plotViewer`` calls ``matplotlib.use("TkAgg")`` at import time (it is the
    interactive viewer, meant to run in its own subprocess). When a test imports
    it in-process, that switches the whole pytest session to TkAgg, so a later
    test rendering a PNG tries to spin up Tk inside the running interpreter and
    fails. Re-pinning Agg before each test isolates that contamination.
    """
    import matplotlib
    matplotlib.use("Agg", force=True)
    yield
