#!/bin/bash
# Dev-time only, never run by the test suite. Builds the reference runatom.x
# with Linux gfortran from third_party/uspp-736/Source plus the idate patch,
# runs 8 library cases and stores their .out and .uspp files (gzipped) here.
# Run from venv/src inside WSL:  bash tests/data/runatom/make_reference.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$(cd "$HERE/../../../third_party/uspp-736" && pwd)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
FLAGS="-O1 -std=legacy -fallow-argument-mismatch"

mkdir "$TMP/build"; cp "$SRC"/Source/*.f "$TMP/build/"
sed -i 's/^      external idate\(\r\?\)$/c     external idate\1/' "$TMP/build/runatom.f"
[ "$(grep -c '^c     external idate' "$TMP/build/runatom.f")" = 1 ]
(cd "$TMP/build" && gfortran $FLAGS -o runatom.x *.f)
X="$TMP/build/runatom.x"
mkdir "$TMP/reform"; cp "$SRC/Utils/reform.f" "$TMP/reform/"          # unmodified
(cd "$TMP/reform" && gfortran $FLAGS -o reform.x reform.f)
R="$TMP/reform/reform.x"

CASES="000-Ex/000-Ex-ca-sp-test 001-H/001-H-gpbe--bm 008-O/008-O-ca--vgrp
       026-Fe/026-Fe-gpw-n-campos 029-Cu/029-Cu-gpbe-n-bm 079-Au/079-Au-ca--ts
       082-Pb/082-Pb-ca-d-vgrp 056-Ba/056-Ba-gpw-nsp-campos"
{ echo "compiler: $(gfortran --version | head -1)"; echo "flags: $FLAGS";
  echo "source: third_party/uspp-736/Source + patches/0001-gfortran-intrinsic-idate.patch"; } > "$HERE/MANIFEST.txt"

for c in $CASES; do
  name=$(basename "$c"); dst="$HERE/$name"; rm -rf "$dst"; mkdir -p "$dst"
  run="$TMP/$name"; cp -r "$SRC/Work/$c" "$run"; cd "$run"
  atom=$(sed -n 's/^ATOM= *//p' Makefile | tr -d '\r')
  cfg=$(sed -n 's/^CFG= *//p' Makefile | tr -d '\r')
  tests=$(sed -n 's/^CFT[0-9]*= *//p' Makefile | tr -d '\r')
  "$X" "${atom}_ae_${cfg}.adat" "${atom}_ae_${cfg}.out" "${atom}_ae_${cfg}.ae" "${atom}_ae_${cfg}.atwf" "${atom}_ae_${cfg}.logd" dummy 2>/dev/null
  "$X" "${atom}_ps.adat" "${atom}_ps.out" "${atom}_ae_${cfg}.ae" "${atom}_ps.atwf" "${atom}_ps.logd" "${atom}_ps.uspp" 2>/dev/null
  done_tests=""
  for t in $tests; do
    [ -f "${atom}_test_${t}.adat" ] || continue
    "$X" "${atom}_ae_${t}.adat" "${atom}_ae_${t}.out" "${atom}_ae_${t}.ae" "${atom}_ae_${t}.atwf" "${atom}_ae_${t}.logd" dummy 2>/dev/null
    "$X" "${atom}_test_${t}.adat" "${atom}_test_${t}.out" "${atom}_ae_${t}.ae" "${atom}_test_${t}.atwf" "${atom}_test_${t}.logd" "${atom}_ps.uspp" 2>/dev/null
    done_tests="$done_tests $t"
  done
  printf '1\n%s\n%s\n' "${atom}_ps.uspp" "${atom}_ps_formatted.uspp" | "$R" > /dev/null
  for f in *.out "${atom}_ps.uspp" "${atom}_ps_formatted.uspp"; do gzip -9 -n -c "$f" > "$dst/$f.gz"; done
  { echo "case: $c"; echo "atom: $atom"; echo "cfg: $cfg"; echo "tests:$done_tests"; } > "$dst/CASE.txt"
  cd "$HERE"
done
echo "reference written to $HERE"
