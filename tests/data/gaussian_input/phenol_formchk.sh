#!/bin/bash
# gQTEA Fukui input builder: convert the three checkpoints to formatted
# checkpoints for Tools > Analysis of Fukui Functions. Run it where Gaussian
# ran (WSL or the cluster), after the jobs finish:  bash phenol_formchk.sh
cd "$(dirname "$0")" || exit 1
if ! command -v formchk >/dev/null 2>&1; then
  echo "formchk not found: set up Gaussian 16 first (e.g. 'module load gaussian', or source
  \$g16root/g16/bsd/g16.profile)." >&2
  exit 1
fi
status=0
convert() {
  if [ ! -f "$1" ]; then
    echo "$1 not found: did that job run? Gaussian stops the remaining linked jobs when one fails; see the .log file." >&2
    status=1
    return
  fi
  formchk "$1" "$2" || { echo "formchk failed for $1" >&2; status=1; }
}
convert "phenol_N.chk" "phenol_N.fchk"
convert "phenol_Np1.chk" "phenol_Np1.fchk"
convert "phenol_Nm1.chk" "phenol_Nm1.fchk"
exit $status
