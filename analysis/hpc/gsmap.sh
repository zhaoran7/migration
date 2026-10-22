#!/usr/bin/env bash
set -euo pipefail
trait=${1:?Usage: gsmap.sh trait sample h5ad annotation layer}
sample=${2:?}; h5ad=${3:?}; annotation=${4:?}; layer=${5:-counts}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# This is the quick_mode interface used for the published analysis, not the newer quick-mode API.
args=(--workdir "${GSMAP_OUT:?}" --sample_name "$sample" --gsMap_resource_dir "${GSMAP_REF:?}"
  --hdf5_path "$h5ad" --annotation "$annotation" --data_layer "$layer"
  --trait_name "$trait" --sumstats_file "${GSMAP_SUM:?}/$trait.sumstats.gz"
  --max_processes "${MOVE_CORES:-1}")
[[ -z ${HOMOLOG:-} ]] || args+=(--homolog_file "$HOMOLOG")
gsmap quick_mode "${args[@]}"
