"""Write-only-on-change for Toga widgets (no Toga import; any window can use it).

Split out of cubeTab, which re-exports it as ``_assign``, so a window that
needs only this helper does not import cubeTab's numba/OpenGL chain.
"""
import math


def assign_if_changed(widget, attr: str, value) -> None:
    """Write a widget property only when it actually changes.

    On WinForms every Toga widget write re-lays out the whole viewer window
    (~70 ms measured, seven tabs of controls). Loading one cube used to write
    ~25 mostly unchanged values in a row, which showed as the window
    flickering for seconds. Numbers are compared with a tolerance because a
    slider reads back a rounded value.
    """
    current = getattr(widget, attr, None)
    if isinstance(value, float) and isinstance(current, (int, float)) and not isinstance(current, bool):
        if math.isclose(current, value, rel_tol=1e-9, abs_tol=1e-9):
            return
    elif current == value:
        return
    setattr(widget, attr, value)
