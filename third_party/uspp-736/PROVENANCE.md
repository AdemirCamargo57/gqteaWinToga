# USPP 7.3.6 in gqteaWinToga

Upstream: Vanderbilt ultrasoft pseudopotential package uspp-736 (D. Vanderbilt,
Rutgers University, December 2006). Licence: GNU GPL v2 (`Source/LICENSE`),
with the additional condition in `AA-README` (no profit-making use, no inclusion
in commercial code without arrangement). Download URL and archive checksum:
**not recorded yet — add them when known.**

The files outside `patches/`, `bin/`, `build_runatom.ps1`, `PROVENANCE.md`,
`SOURCE_SHA256.txt` and `.gitattributes` are unmodified upstream files
(`SOURCE_SHA256.txt`, checked by tests/test_runatomVendor.py).

## Change made to build runatom.exe

`patches/0001-gfortran-intrinsic-idate.patch` comments out `external idate` in
`Source/runatom.f`. With gfortran the declaration hides the IDATE intrinsic and
the link fails (`undefined reference to 'idate_'`). The change only affects the
date printed in the .out header and stored in the .uspp file. It is applied to a
temporary copy by `build_runatom.ps1`; the file in `Source/` stays unmodified.

## Known defects of the upstream input library (left as shipped)

- `Work/023-V/023-V-gpw-nsp-campos/V_ps.adat` line 2: the minus sign of emin
  falls in the rlogd field; runatom stops with "Bad value during floating point read".
- `Work/041-Nb/041-Nb-gpw-nsp-campos/041-Nb-gpw-nsp-campos_us_pw91_13elec` (and
  the copy in the .tar next to it): 10240 bytes of zeros, not a pseudopotential.
- `006-C-gpw-n-campos` and `024-Cr-gpw-nsp-campos` have ifpcor = 0 although their
  names and READMEs claim a nonlinear core correction; `048-Cd-gpw-campos` has
  ifpcor = 1 but no flags field in its name.
- `039-Y-ca-sp-vgrp` is a retracted entry (README only).
- `Doc/TUTORIAL` suggests starting with ifqopt = -1, which runatom refuses for
  every GGA functional (rwsubs.f: ifqopt < 3 with exfact >= 1 is an error).
