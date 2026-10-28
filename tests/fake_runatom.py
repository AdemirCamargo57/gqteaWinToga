"""Stand-in for runatom in tests (never the real program).

Called like runatom: input out ae atwf logd uspp. FAKE_RUNATOM_MODE selects:
ok, error, stop0, runtime, negative, hang. Every call is logged to fake_calls.jsonl.
"""

import json
import os
import sys
import time

SAMPLE_OUT = """ comparison of all-electron and pseudo eigenvalues (Ry)

   nlm    all-elec      pseudo        diff

   200   -1.000000   -1.000001   -0.000001
   210   -0.500000   -0.500000    0.000000

 -----------------------------------------------
   eigenvalues of US hamiltonian
   l= 0      -1.000000   0.300000
   l= 1      -0.500000   0.700000
"""

mode = os.environ.get("FAKE_RUNATOM_MODE", "ok")
inp, out, ae, atwf, logd, uspp = sys.argv[1:7]
with open("fake_calls.jsonl", "a") as log:
    log.write(json.dumps({"argv": sys.argv[1:], "cwd": os.getcwd()}) + "\n")
kind = "ae" if uspp == "dummy" else ("test" if "_test_" in inp else "pseudo")
err = sys.stderr
err.write("beginning execution pseudopotential program version  7.3.6\n")
err.flush()
if mode == "runtime":
    err.write(f"At line 165 of file rwsubs.f (unit = 10, file = '{inp}')\n"
              "Fortran runtime error: Bad value during floating point read\n")
    sys.exit(2)
if mode == "negative":
    err.write(" *** stopping in vhxc: negative densities ***\n")
    sys.exit(1)
if mode == "stop0":
    sys.exit(0)
for cycle in range(1, 4):
    err.write(f"completed self-consistent cycle{cycle:4d} delta ={10.0 ** -cycle:15.7E}\n".replace("E", "D"))
    err.flush()
if mode == "hang":
    for name in (out, ae if kind == "ae" else uspp):
        open(name, "w").close()
    while True:
        time.sleep(0.05)
if mode == "error":
    with open(out, "w") as f:
        f.write(" ***error in subroutine readin\n nang =  5 but idim5 =  4\n")
    sys.exit(1)
err.write({"ae": " all electron calculation completed\n",
           "pseudo": " writing the pseudopotential to file\n",
           "test": " pseudopotential test converged\n"}[kind])
with open(out, "w") as f:
    f.write(" value of mesh generated in rinit is   851\n" if kind == "ae" else SAMPLE_OUT)
if kind == "ae":
    open(ae, "wb").write(b"ae")
if kind == "pseudo":
    open(uspp, "wb").write(b"uspp")
