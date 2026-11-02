#!/bin/sh
# Stream the ACT DR6 lensing simulations one at a time into data/act_dr6/sims (160 MB each; 64 GB in all, against 28 GB of
# disk): fetch sim i when fewer than KEEP unprocessed files wait, where a consumer (kappa_templates.py sims) deletes each file
# once used and leaves done_NNNN. usage: sh act_dr6_simstream.sh <variant> <first> <last> [keep]; writes ALLDONE to
# data/act_dr6/simstream_<variant>.log at the end. (R-122's simstream.sh, lost with the gitignored data/, rewritten.)
set -u
variant="$1"; first="$2"; last="$3"; keep="${4:-3}"
here=$(cd "$(dirname "$0")" && pwd)
dir="$here/data/act_dr6/sims"; mkdir -p "$dir"
log="$here/data/act_dr6/simstream_$variant.log"
i=$first
while [ "$i" -le "$last" ]; do
  n=$(printf '%04d' "$i")
  f="$dir/kappa_alm_sim_act_dr6_lensing_v1_${variant}_$n.fits"
  if [ -e "$dir/done_${variant}_$n" ] || [ -e "$f" ]; then i=$((i + 1)); continue; fi
  while [ "$(ls "$dir" | grep -c "kappa_alm_sim_act_dr6_lensing_v1_${variant}_.*\.fits$")" -ge "$keep" ]; do sleep 5; done
  if sh "$here/act_dr6_fetch.sh" "maps/$variant/simulations/kappa_alm_sim_act_dr6_lensing_v1_${variant}_$n.fits" "$f"; then
    echo "$(date -u +%H:%M:%S) fetched $n" >> "$log"
  else
    echo "$(date -u +%H:%M:%S) FAILED $n" >> "$log"
  fi
  i=$((i + 1))
done
echo "ALLDONE" >> "$log"
