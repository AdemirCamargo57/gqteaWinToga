"""Pytest bootstrap for the gqteaWinToga source tree.

The tools use flat, bare-name imports (``from framesCounter import ...``), so the
source directory (the parent of this ``tests/`` folder) must be on ``sys.path``.
We also pin matplotlib to the non-interactive Agg backend before anything imports
pyplot, so tests never try to open a GUI window.
"""
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)
