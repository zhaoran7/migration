#!/usr/bin/env bash
set -euo pipefail
# The code and data directories below must already exist in the RAP project mount.
module=${1:?Usage: rap/submit.sh module [arguments]}
shift
: "${RAP_IMAGE:?DNAnexus file ID of the built Docker image archive}"
: "${RAP_CODE:?absolute /mnt/project/... path to the uploaded analysis directory}"
: "${RAP_INPUT:?absolute /mnt/project/... path to the prepared data bundle directory}"
: "${RAP_DEST:?DNAnexus project:/folder destination}"
[[ $module =~ ^[a-z]+$ ]]
printf -v args '%q ' "$@"
printf -v code '%q' "$RAP_CODE"
printf -v data '%q' "$RAP_INPUT/data"
printf -v previous '%q' "$RAP_INPUT/results"
printf -v settings '%q' "$RAP_INPUT/paths.sh"
command="set -euo pipefail; work=\$(mktemp -d /home/dnanexus/move.XXXXXX); \
cp -a $data \"\$work/data\"; \
chmod -R u+rwX \"\$work/data\"; \
if [ -d $previous ]; then cp -a $previous /tmp/results; chmod -R u+rwX /tmp/results; fi; \
if [ -f $settings ]; then set -a; source $settings; set +a; fi; \
export MOVE_ROOT=\"\$work\" MOVE_CODE=$code MOVE_RESULTS=/tmp/results MOVE_SOURCE=/tmp/source; \
bash $code/rap/run.sh $module $args"
dx run swiss-army-knife -iimage_file="$RAP_IMAGE" -imount_inputs=true \
  -icmd="$command" --instance-type "${RAP_INSTANCE:-mem2_ssd1_v2_x8}" \
  --destination "$RAP_DEST" --name "move.$module" --brief -y
