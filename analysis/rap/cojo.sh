#!/usr/bin/env bash
set -euo pipefail
trait=${1:?Usage: cojo.sh trait chromosome}
chr=${2:?}
[[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ ]]
out="${MOVE_RESULTS:-$MOVE_ROOT/results}/cojo/$trait"
gcta64 --bfile "${LD_PREFIX:?}$chr" --maf .01 --diff-freq .2 \
  --cojo-file "$out/chr$chr.ma" --cojo-slct --out "$out/chr$chr"
