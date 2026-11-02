#!/bin/sh
# R-127: every delay fit, three chains side by side (4 cores). Run from research/cosmo-discovery.
D="python3 ztf_lens_delays.py"
L=logs
(
$D j1004 C A --band zr --tmin 600 --tmax 1050 --step 2 --mlknot 365 --nmock 30 --plot > $L/delay_j1004_CA_zr.log 2>&1
$D j1004 C A --band zg --tmin 600 --tmax 1050 --step 2 --mlknot 365 --nmock 30 > $L/delay_j1004_CA_zg.log 2>&1
$D j1004 C B --band zr --tmin 600 --tmax 1050 --step 2 --mlknot 365 --nmock 30 --plot > $L/delay_j1004_CB_zr.log 2>&1
$D j1004 C B --band zg --tmin 600 --tmax 1050 --step 2 --mlknot 365 --nmock 30 > $L/delay_j1004_CB_zg.log 2>&1
$D j1004 B A --band zr --tmin -100 --tmax 150 --step 1 --mlknot 365 --nmock 30 --plot > $L/delay_j1004_BA_zr.log 2>&1
$D j1004 B A --band zg --tmin -100 --tmax 150 --step 1 --mlknot 365 --nmock 30 > $L/delay_j1004_BA_zg.log 2>&1
$D j1004 C D --band zr --tmin 1900 --tmax 2720 --step 4 --mlknot 365 --nmock 30 --min-overlap 120 --plot > $L/delay_j1004_CD_zr.log 2>&1
$D j1029 A B --band zr --tmin 450 --tmax 1050 --step 2 --mlknot 365 --nmock 30 --plot > $L/delay_j1029_AB_zr.log 2>&1
$D j1029 A B --band zg --tmin 450 --tmax 1050 --step 2 --mlknot 365 --nmock 30 > $L/delay_j1029_AB_zg.log 2>&1
echo chain1 done > $L/chain1.done
) &
(
$D j2222 A B --band zr --tmin -120 --tmax 200 --step 1 --mlknot 365 --bin 5 --nmock 30 --plot > $L/delay_j2222_AB_zr.log 2>&1
$D j2222 A B --band zg --tmin -120 --tmax 200 --step 1 --mlknot 365 --bin 5 --nmock 30 > $L/delay_j2222_AB_zg.log 2>&1
$D j2222 C A --band zr --tmin 450 --tmax 1000 --step 2 --mlknot 365 --bin 5 --nmock 30 --plot > $L/delay_j2222_CA_zr.log 2>&1
$D j2222 C A --band zg --tmin 450 --tmax 1000 --step 2 --mlknot 365 --bin 5 --nmock 30 > $L/delay_j2222_CA_zg.log 2>&1
$D j0335 A B --band zr --tmin 100 --tmax 900 --step 2 --mlknot 365 --bin 5 --nmock 30 --plot > $L/delay_j0335_AB_zr.log 2>&1
$D j0335 A C --band zr --tmin -500 --tmax 250 --step 2 --mlknot 365 --bin 5 --nmock 30 --plot > $L/delay_j0335_AC_zr.log 2>&1
$D j0335 A B --band zg --tmin 100 --tmax 900 --step 2 --mlknot 365 --bin 5 --nmock 30 > $L/delay_j0335_AB_zg.log 2>&1
$D j0335 A C --band zg --tmin -500 --tmax 250 --step 2 --mlknot 365 --bin 5 --nmock 30 > $L/delay_j0335_AC_zg.log 2>&1
echo chain2 done > $L/chain2.done
) &
(
$D j1153 A1 A2 --band zr --tmin 700 --tmax 1800 --step 4 --mlknot 365 --bin 5 --nmock 30 --plot > $L/delay_j1153_A1A2_zr.log 2>&1
$D j1153 A1 A3 --band zr --tmin 400 --tmax 1200 --step 4 --mlknot 365 --bin 5 --nmock 30 --plot > $L/delay_j1153_A1A3_zr.log 2>&1
$D j1153 A1 A2 --band zg --tmin 700 --tmax 1800 --step 4 --mlknot 365 --bin 5 --nmock 30 > $L/delay_j1153_A1A2_zg.log 2>&1
$D j1153 A1 A3 --band zg --tmin 400 --tmax 1200 --step 4 --mlknot 365 --bin 5 --nmock 30 > $L/delay_j1153_A1A3_zg.log 2>&1
$D j1326 A B --band zr --tmin -2500 --tmax 2500 --step 4 --mlknot 365 --bin 5 --nmock 30 --plot > $L/delay_j1326_AB_zr.log 2>&1
$D j1326 A B --band zg --tmin -2500 --tmax 2500 --step 4 --mlknot 365 --bin 5 --nmock 30 > $L/delay_j1326_AB_zg.log 2>&1
until grep -q done $L/fetch_j0542.log; do sleep 20; done
python3 ztf_lens_phot.py j0542 > $L/phot_j0542.log 2>&1
$D j0542 AB C --band zr --tmin -2500 --tmax 2500 --step 4 --mlknot 365 --bin 5 --nmock 30 --plot > $L/delay_j0542_ABC_zr.log 2>&1
$D j0542 AB C --band zg --tmin -2500 --tmax 2500 --step 4 --mlknot 365 --bin 5 --nmock 30 > $L/delay_j0542_ABC_zg.log 2>&1
echo chain3 done > $L/chain3.done
) &
wait
