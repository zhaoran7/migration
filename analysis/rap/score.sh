#!/usr/bin/env bash
set -euo pipefail
trait=${1:?Usage: score.sh score chromosome}
chr=${2:?}
[[ $trait =~ ^prs_[a-z0-9_]+$ && $chr =~ ^([1-9]|1[0-9]|2[0-2])$ ]]
out="${MOVE_RESULTS:-$MOVE_ROOT/results}/scores/$trait"
mkdir -p "$out"
plink2 --pfile "${IMP_PREFIX:?}$chr" --keep "${KEEP:?}" --threads "${MOVE_CORES:-1}" \
  --score "${WEIGHTS:?}/$trait.score" 1 2 3 header-read ignore-dup-ids cols=scoresums \
  --out "$out/chr$chr"
