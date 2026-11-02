#!/bin/sh
# R-153: accumulate the detector pattern over one image list with N workers.  Usage: sh alp_prun.sh <name> [n]
name=$1; nw=${2:-4}
cd "$(dirname "$0")"
mkdir -p data/alp/logs
i=0
while [ $i -lt $nw ]; do
  nohup python -u alp_pattern.py data/alp/list_$name.csv $name $i $nw > data/alp/logs/pat_$name.$i.log 2>&1 &
  i=$((i+1))
done
