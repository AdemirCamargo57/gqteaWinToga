"""A MolecularViewer owned by a tool window that follows its atom selection.

The 3D Molecular Viewer's own window (MolecularViewerUI) handles selection
changes itself; a viewer opened from another tool window (the Fukui
analysis, the All bond / angle / dihedral analyses) reports them through
``selection_listener``, called on ``loop`` (the Toga loop) because a
Ctrl/Shift-click changes the selection on the render thread.
"""
from typing import Callable, Optional

import glfw

from molecularViewer import MolecularViewer


class SelectionViewer(MolecularViewer):
    """MolecularViewer with a selection listener and a thread-safe close request."""

    def __init__(self):
        # Before super(): the base constructor may already touch the selection hooks.
        self.selection_listener: Optional[Callable[[], None]] = None
        super().__init__()

    def _notify_selection_changed(self):
        listener, loop = self.selection_listener, getattr(self, "loop", None)
        if listener is None or loop is None or loop.is_closed():
            return
        loop.call_soon_threadsafe(listener)

    def request_close(self):
        """Ask GLFW to close the window (glfwSetWindowShouldClose is thread-safe)."""
        if self.window is not None:
            try:
                glfw.set_window_should_close(self.window, True)
            except Exception:
                pass
