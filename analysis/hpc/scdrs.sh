#!/usr/bin/env bash
set -euo pipefail
atlas=${1:?Usage: scdrs.sh atlas h5ad group raw_counts}
h5ad=${2:?}; group=${3:?}; raw=${4:?True or False}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
out="${SCDRS_OUT:?}"
mkdir -p "$out/group/$atlas"
scdrs compute-score --h5ad-file "$h5ad" --h5ad-species human \
  --gs-file "${SCDRS_GENES:?}" --gs-species human --ctrl-match-opt mean_var --weight-opt vs \
  --flag-filter-data True --flag-raw-count "$raw" --n-ctrl 1000 \
  --flag-return-ctrl-raw-score False --flag-return-ctrl-norm-score True \
  --out-folder "$out/group/$atlas"
scdrs perform-downstream --h5ad-file "$h5ad" \
  --score-file "$out/group/$atlas/@.full_score.gz" --out-folder "$out/group/$atlas" \
  --group-analysis "$group" --flag-filter-data True --flag-raw-count "$raw" \
  --knn-n-neighbors 15 --knn-n-pcs 20
python "$(dirname "$0")/umap.py" "$h5ad" "${group%%,*}" "$out/coord/$atlas.umap.tsv.gz"
