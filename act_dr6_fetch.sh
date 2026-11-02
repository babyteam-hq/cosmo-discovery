#!/bin/sh
# Fetch one file of the ACT DR6 lensing v1 release from NERSC's portal (R-122 kept this in the gitignored data/; R-129 moved
# it into the repo; R-146 added the stall abort: < 1 MB/s for 20 s retries). usage: sh act_dr6_fetch.sh <path under dr6_lensing_v1/> <destination>; retries 5 times.
set -eu
src="$1"; dst="$2"
base="https://portal.nersc.gov/project/act/dr6_lensing_v1"
i=0
while [ $i -lt 5 ]; do
  if curl -sS -f -m 900 -Y 1000000 -y 20 -o "$dst.part" "$base/$src"; then mv "$dst.part" "$dst"; exit 0; fi
  i=$((i + 1)); sleep $((i * 5))
done
rm -f "$dst.part"; echo "act_dr6_fetch: failed $src" >&2; exit 1
