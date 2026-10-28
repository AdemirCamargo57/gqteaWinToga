Reference outputs of the original USPP 7.3.6 runatom, and of reform (itype 1,
the formatted file CPMD reads), from a Linux gfortran build (see MANIFEST.txt),
for 8 library cases chosen to exercise every input branch:
LDA/PBE/PW91, semicore, nonlinear core correction, per-L rinner, iptype 2/3/4,
ifqopt 0/2/3, nfix constraints and transferability tests.
Regenerate only with make_reference.sh (WSL with gfortran); never by hand.
The .uspp files carry the run date in record 1; comparisons mask it.
