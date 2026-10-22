#!/usr/bin/env bash
set -euo pipefail
# One independent job per call. No dependency jobs or automatic polling.
name=${1:?Usage: submit.sh name command [arguments]}
shift
log=${MOVE_LOG:-$PWD/log}
mkdir -p "$log"
printf -v command '%q ' "$@"
bsub -q "${QUEUE:-ser}" -n "${MOVE_CORES:-1}" -R 'span[hosts=1]' \
  -J "$name" -o "$log/$name.%J.log" -e "$log/$name.%J.err" "$command"
