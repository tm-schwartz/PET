#!/usr/bin/env bash
# Stage pet-fdg cohort inputs onto Merlin's shared /home-nfs.
# Run from the NFS-mounted workspace host (which sees /nfs/masi), NOT on merlin.
#   container/merlin_stage.sh [N_SUBJECTS] [MERLIN_LOGIN]
# Selection logic mirrors batch_run.pick_series: path-sorted glob, first
# whole-body AC series per subject, first N subjects.
set -euo pipefail
N=${1:-50}
MERLIN=${2:-merlin-01}
BIDS=/nfs/masi/schwat1/bids
REMOTE=/home-nfs/schwat1/pet-fdg

all=$(
  find "$BIDS" -mindepth 3 -maxdepth 3 -path '*/PT/*' \
       -name '*AC*3DWB*PT.nii.gz' | LC_ALL=C sort |
  awk '{ mrn=$0; sub(/.*sub-/,"",mrn); sub(/_.*/,"",mrn);
         if (!(mrn in seen)) { seen[mrn]=1; print } }'
)
mapfile -t all_arr <<< "$all"
series=$(printf '%s\n' "${all_arr[@]:0:$N}")
echo "$series" | wc -l | xargs echo "series to stage:"

mkdir -p /tmp/pet_stage
echo "$series" | sed "s|$BIDS/||" > /tmp/pet_stage/series.txt
ssh "$MERLIN" "mkdir -p $REMOTE/bids $REMOTE/out $REMOTE/templates"
rsync -a --info=progress2 --files-from=/tmp/pet_stage/series.txt \
  "$BIDS/" "$MERLIN:$REMOTE/bids/"

rsync -a --partial --info=progress2 \
  /home-local/schwat1/pet-fdg.sif "$MERLIN:$REMOTE/pet-fdg.sif"
rsync -a --info=progress2 \
  /home-local/schwat1/pet-fdg-test/templates/ "$MERLIN:$REMOTE/templates/"

ssh "$MERLIN" "mkdir -p $REMOTE/out && cp /dev/null $REMOTE/out/.ready && \
  chmod +x $REMOTE/pet-fdg.sif 2>/dev/null || true"
echo "staged to $MERLIN:$REMOTE"
