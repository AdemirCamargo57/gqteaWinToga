"""
Locating and running Gaussian's stand-alone utilities (cubegen, formchk).

No Toga here: the ESP tool, the Fukui input builder and their tests share it.
"""

import os
import shutil
from typing import Optional

# Where G09W / G16W install themselves on Windows by default.
DEFAULT_GAUSSIAN_FOLDERS = (r"C:\G09W", r"C:\G16W")


def find_gaussian_utility(name: str, explicit: Optional[str] = None) -> Optional[str]:
    """Path of a Gaussian utility: explicit, $GAUSS_EXEDIR, PATH, usual folders; None if absent."""
    exe = name if name.lower().endswith(".exe") else name + ".exe"
    candidates = []
    if explicit and explicit.strip():
        candidates.append(explicit.strip())
    exedir = os.environ.get("GAUSS_EXEDIR")
    if exedir:
        candidates += [os.path.join(d, exe) for d in exedir.split(os.pathsep) if d]
    which = shutil.which(name)
    if which:
        candidates.append(which)
    candidates += [os.path.join(d, exe) for d in DEFAULT_GAUSSIAN_FOLDERS]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None
