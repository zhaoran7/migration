#!/usr/bin/env bash
set -euo pipefail
trait=${1:?Usage: weights.sh score chromosome discovery_N}
chr=${2:?}
n=${3:?}
[[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ ]]
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python "${PRSCS:?}/PRScs.py" --ref_dir="${PRS_REF:?}" --bim_prefix="${PRS_BIM:?}" \
  --sst_file="${PRS_SUM:?}/$trait.txt" --n_gwas="$n" \
  --out_dir="${WEIGHTS:?}/$trait" --chrom="$chr"
