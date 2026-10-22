#!/usr/bin/env bash
set -euo pipefail
export MOVE_CODE="$(cd "$(dirname "$0")" && pwd)"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ ${1:--h} == -h || ${1:-} == --help ]]; then
  printf '%s\n' 'Usage: bash analysis/run.sh module [arguments]'
  awk -F '\t' 'NR>1 && $5!="internal" {
    if ($2!=group) {group=$2; printf "\n%s\n",group}
    printf "  %-46s %s [%s]\n",$5,$4,$6
  }' "$MOVE_CODE/modules.tsv"
  printf '\n%s\n' 'RAP genotype stages: gwas 1|2|sum trait [chr]; score trait chr'
  exit 0
fi
module=$1
shift
if [[ $module == gwas || $module == score ]]; then
  exec bash "$MOVE_CODE/rap/$module.sh" "$@"
fi
path=$(awk -F '\t' -v module="$module" 'NR>1 && $1==module && $5!="internal" {print $3}' "$MOVE_CODE/modules.tsv")
[[ -n $path ]] || { printf 'Unknown module: %s\n' "$module" >&2; exit 2; }
exec "${RSCRIPT:-Rscript}" "$MOVE_CODE/$path" "$@"
