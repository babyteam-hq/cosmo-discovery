#!/bin/sh
# Stream the ACT DR6 lensing simulation inputs (input_kappa_alm_sim_NNNN.fits, 149 MB each) into data/act_dr6/inputs for the
# transfer function, as act_dr6_simstream.sh does for the reconstructions (consumer: z5_lens.py simstore --map input).
# usage: sh act_dr6_inputstream.sh <first> <last> [keep]; writes ALLDONE to data/act_dr6/simstream_input.log at the end.
set -u
first="$1"; last="$2"; keep="${3:-3}"
here=$(cd "$(dirname "$0")" && pwd)
dir="$here/data/act_dr6/inputs"; mkdir -p "$dir"
log="$here/data/act_dr6/simstream_input.log"
i=$first
while [ "$i" -le "$last" ]; do
  n=$(printf '%04d' "$i")
  f="$dir/input_kappa_alm_sim_$n.fits"
  if [ -e "$dir/done_input_$n" ] || [ -e "$f" ]; then i=$((i + 1)); continue; fi
  while [ "$(ls "$dir" | grep -c "input_kappa_alm_sim_.*\.fits$")" -ge "$keep" ]; do sleep 5; done
  if sh "$here/act_dr6_fetch.sh" "sim_inputs/kappa_alm/input_kappa_alm_sim_$n.fits" "$f"; then
    echo "$(date -u +%H:%M:%S) fetched $n" >> "$log"
  else
    echo "$(date -u +%H:%M:%S) FAILED $n" >> "$log"
  fi
  i=$((i + 1))
done
echo "ALLDONE" >> "$log"
