#!/bin/sh
# R-143: Gaia DR3 G<16 over the R-138 field box via IRSA TAP: 10x10 deg RA/Dec tiles (RA 144-254, Dec -12..48), 4 at a time.
# Great-circle polygon edges overlap neighbours slightly; the loader keeps each row only inside its own tile's ranges.
# Usage: sh spherex_gaia_fetch.sh data/gaia
cd "$1" || exit 1
for ra in $(seq 144 10 244); do for dec in -12 -2 8 18 28 38; do
  f=g_${ra}_${dec}.csv; [ -s $f ] && continue
  r2=$((ra+10)); d2=$((dec+10))
  ( curl -sS --max-time 900 "https://irsa.ipac.caltech.edu/TAP/sync" --data-urlencode "FORMAT=csv" --data-urlencode "QUERY=SELECT ra, dec, phot_g_mean_mag FROM gaia_dr3_source WHERE CONTAINS(POINT('ICRS',ra,dec), POLYGON('ICRS', $ra,$dec, $r2,$dec, $r2,$d2, $ra,$d2))=1 AND phot_g_mean_mag < 16" > $f.tmp && head -c 6 $f.tmp | grep -q "^ra,dec" && mv $f.tmp $f ) &
  while [ $(jobs -p | wc -l) -ge 4 ]; do sleep 2; done
done; done; wait; echo DONE
