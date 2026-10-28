"""Stand-in for USPP's reform (itype 1) in tests; never the real program.

Reads itype, input and output names from stdin like reform.f. FAKE_REFORM_MODE:
ok (writes the output with CRLF line endings, as MinGW gfortran does), fail
(missing input), silent (exit 0 without the completion message).
"""

import os
import sys

mode = os.environ.get("FAKE_REFORM_MODE", "ok")
itype, source, target = (sys.stdin.readline().strip() for _ in range(3))
print(" selected itype =", itype)
if os.path.exists(target):        # reform.f opens its output with status='new'
    sys.stderr.write(f"Fortran runtime error: Cannot open file '{target}': File exists\n")
    sys.exit(2)
if mode == "fail" or not os.path.exists(source):
    sys.stderr.write(f"Fortran runtime error: Cannot open file '{source}': No such file or directory\n")
    sys.exit(2)
with open(target, "w", newline="\r\n") as f:
    f.write("    7    3    6\nformatted\n")
if mode != "silent":
    print("1pseudopotential has been written out")
