#!/usr/bin/env bash
set -euo pipefail
# Run within RAP Swiss Army Knife / a RAP workstation, not on a local UKB mirror.
module=${1:?Usage: rap/run.sh module [arguments]}
shift
: "${MOVE_ROOT:?set a writable RAP working directory}"
: "${MOVE_CODE:?set the analysis directory inside the RAP job}"
export MOVE_PHE=${MOVE_PHE:-$MOVE_ROOT/data/phe}
export MOVE_UKB=${MOVE_UKB:-$MOVE_ROOT/data/ukb}
export MOVE_RESULTS=${MOVE_RESULTS:-$MOVE_ROOT/results}
export MOVE_SOURCE=${MOVE_SOURCE:-$MOVE_ROOT/source}
export MOVE_REGION=${MOVE_REGION:-$MOVE_RESULTS/region}
bash "$MOVE_CODE/run.sh" "$module" "$@"
case "$module:${1:-}" in ukb:*|geo:*|pheno:*|import:*|scores:*|gwas:sum)
  mkdir -p /tmp/prepared/data
  cp -a "$MOVE_PHE" /tmp/prepared/data/phe
  ;;
esac
