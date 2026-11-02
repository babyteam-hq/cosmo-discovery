#!/bin/sh
# R-153: reduce one dwarf's (or fake field's) image list with N workers in parallel, subtracting the fixed detector
# pattern built without that dwarf's own images (alp_pattern.py; fake fields use every dwarf's).
# Usage: sh alp_run.sh <name> <ra> <dec> [n_workers]   -> data/alp/redp/<name>/
name=$1; ra=$2; dec=$3; nw=${4:-4}
cd "$(dirname "$0")"
mkdir -p data/alp/redp/$name data/alp/logs
i=0
while [ $i -lt $nw ]; do
  nohup python -u alp_reduce.py --pattern data/alp/pattern $name mean data/alp/list_$name.csv data/alp/redp/$name $ra $dec $i $nw > data/alp/logs/$name.$i.log 2>&1 &
  i=$((i+1))
done
