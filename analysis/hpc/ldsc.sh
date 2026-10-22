#!/usr/bin/env bash
set -euo pipefail
stage=${1:?Usage: ldsc.sh munge|h2|rg trait [second trait]}
trait=${2:?}
app=${LDSC:?path to the LDSC checkout}
ref=${LDSC_REF:?directory with eur_w_ld_chr and w_hm3.snplist}
out=${MOVE_LDSC:-$PWD/ldsc}
mkdir -p "$out"/{sumstats,h2,rg}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
case "$stage" in
  munge)
    bash "$(dirname "$0")/../run.sh" ldsc prepare "${GWAS_CLEAN:?}/$trait.gz" \
      "${LDSC_MAP:?HapMap3 GRCh38 chromosome-position-rsID map}" "$out/sumstats/$trait.input.gz"
    python "$app/munge_sumstats.py" --sumstats "$out/sumstats/$trait.input.gz" \
      --snp SNP --a1 EA --a2 NEA --N-col N --signed-sumstats BETA,0 --p P \
      --merge-alleles "$ref/w_hm3.snplist" --chunksize 10000 --out "$out/sumstats/$trait"
    ;;
  h2)
    python "$app/ldsc.py" --h2 "$out/sumstats/$trait.sumstats.gz" \
      --ref-ld-chr "$ref/eur_w_ld_chr/" --w-ld-chr "$ref/eur_w_ld_chr/" --out "$out/h2/$trait"
    ;;
  rg)
    other=${3:?second trait required}
    python "$app/ldsc.py" --rg "$out/sumstats/$trait.sumstats.gz,$out/sumstats/$other.sumstats.gz" \
      --ref-ld-chr "$ref/eur_w_ld_chr/" --w-ld-chr "$ref/eur_w_ld_chr/" --out "$out/rg/$trait.$other"
    ;;
  *) exit 2 ;;
esac
