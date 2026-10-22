#!/usr/bin/env bash
set -euo pipefail

# 1. Explicit stages; submit step 2 only after step 1 has finished.
stage=${1:?Usage: gwas.sh 1|2|sum migration_distance|delta_tdi_ukb_cond [chromosome]}
trait=${2:?}
case "$trait" in migration_distance|delta_tdi_ukb_cond) ;; *) exit 2 ;; esac
root=${MOVE_ROOT:?}
pheno="$root/data/phe/$trait.pheno"
name=$trait
[[ $trait != migration_distance ]] || name=migration_distance_pc40
out="${MOVE_RESULTS:-$root/results}/gwas/$name"
mkdir -p "$out"
cd "$out"
cov=age,age2,sex
for pc in {1..40}; do cov+=",PC$pc"; done
[[ $trait != delta_tdi_ukb_cond ]] || cov+=,tdi_birth
common=(--phenoFile "$pheno" --phenoCol "$trait" --qt --apply-rint
  --covarFile "$pheno" --covarColList "$cov" --catCovarList center,gen_batch
  --maxCatLevels 120 --threads "${MOVE_CORES:-1}" --bsize 1000)
case "$stage" in
  1)
    regenie --step 1 --bed "${TYP_BED:?}" --lowmem --lowmem-prefix tmp --out step1 "${common[@]}"
    ;;
  2)
    chr=${3:?chromosome required}
    [[ $chr =~ ^([1-9]|1[0-9]|2[0-2])$ ]]
    test -s step1_pred.list
    regenie --step 2 --pgen "${IMP_PREFIX:?}$chr" --extract "${EXTRACT_PREFIX:?}$chr.extract" --write-samples \
      --pred step1_pred.list --out "chr$chr" "${common[@]}"
    ;;
  sum)
    bash "$MOVE_CODE/run.sh" sumstats regenie "$out" "$trait" "$name"
    ;;
  *) exit 2 ;;
esac
