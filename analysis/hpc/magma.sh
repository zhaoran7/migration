#!/usr/bin/env bash
set -euo pipefail
stage=${1:?Usage: magma.sh annotate|genes trait}
trait=${2:?}
out="${MAGMA_OUT:?}"
mkdir -p "$out"
case "$stage" in
  annotate)
    magma --annotate window=10,10 --snp-loc "${MAGMA_SNPS:?}" --gene-loc "${MAGMA_GENES:?}" --out "$out/annotation"
    ;;
  genes)
    magma --bfile "${MAGMA_LD:?}" --pval "${MAGMA_SUM:?}/$trait.tsv" use=SNP,P ncol=N \
      --gene-annot "$out/annotation.genes.annot" --out "$out/$trait"
    ;;
  *) exit 2 ;;
esac
