# Results (R-127, cosmo-deep-7, 2026-11-01)

ZTF time delays of the cluster- and group-scale lensed quasars from scene-model photometry on the public difference
images, validated on the published delays of SDSS J1004+4112 and J1029+2623, and each new delay tested against the lens
model's prediction published before any measurement. This file is the worker's own. Scripts: `ztf_lens_targets.py`
(the fixed target list, positions, published and predicted delays), `ztf_lens_fetch.py` (IRSA IBE cutouts),
`ztf_lens_phot.py` (the scene model), `ztf_lens_delays.py` (curve shifting with mock errors). Data under
`data/ztf_lens/<target>/` (gitignored): `epochs.csv`, `cut/` (one difference-image cutout and PSF stamp per epoch),
`ref/` (reference-image cutout and PSF-fit catalogue per field/ccd/quadrant/filter), `phot.csv`, `delay_*.json`.

## Trajectory check (run in the dispatching session, 2026-11-01, context-free Opus; not rerun here)

"VERDICT: partly novel: measuring wide-separation lensed-quasar time delays from public ZTF forced-photometry light curves
is already done (Nunez-Pizarro et al. 2026, A&A 706 L12, GraL J1651-0417, ~10" quad; arXiv id not found), and PSF-forward-
model deblending for delays is established on dedicated monitoring (TDCOSMO XVII, arXiv:2504.02932); what remains open is
ZTF delays for cluster/group-scale systems (J1004, J1029, J2222, J1326, COOL J0542 etc.), scene-model forced photometry on
ZTF difference images, and a sample-wide comparison with pre-published cluster-model predictions". Prior work the brief
names, cited not rerun: arXiv:2206.08597 (J1004 delays), arXiv:1505.06187 (J2222 delays), Fohlmeister et al. 2013
(arXiv:1207.5776, J1029), arXiv:1609.08848 (J2222 D-F predictions), arXiv:2209.03972 (COOL J0542), arXiv:2305.14317
(COOL J0335 predictions), arXiv:2605.13514 (COOL J1153 predictions), arXiv:2608.07470 (J1651 ZTF delays vs HST models).

## Reachability from the cloud

IRSA answers: the IBE image search (`/ibe/search/ztf/products/sci?POS=`), the science-, difference- and reference-image
cutouts (`?center=&size=` on the `.fits` and `.fits.fz` products) and the reference PSF-fit catalogues, at ~3-5 files/s
over 6-10 connections. The ZTF forced-photometry service (ztfweb.ipac.caltech.edu) needs an account (401), so the
scene model on the difference images is the product, not a fallback. arXiv e-prints answer 200 (all positions, delays
and predictions below come from the papers' LaTeX sources); VizieR via IDIA answers (Gaia DR3, DR16Q); CDS is cut.

## Targets and pairs, fixed before any light curve was extracted

ZTF public epochs in g+r (i not used), JD-2450000 of first and last epoch, and the baseline; "primary set" is the
(field, ccd, quadrant) holding most epochs, the one reference the light curves are built on.

| system | z_s / z_l | images in ZTF | g+r epochs | first-last | baseline (d) | primary set |
|---|---|---|---|---|---|---|
| SDSS J1004+4112 | 1.734 / 0.68 | A, B, C, D (r 18.7-20.5; A-B 3.8") | 2295 | 8203-11077 | 2874 | 711/10/4 (94 %) |
| SDSS J1029+2623 | 2.197 / 0.584 | A, B(+C, 1.8" off, blended) (r 18.3-18.7) | 1498 | 8203-11078 | 2875 | 621/11/3 (89 %) |
| SDSS J2222+2745 | 2.805 / 0.49 | A, B, C (r 20.7-21.1); D, E, F at r 23-25 are below ZTF's per-epoch depth | 3468 | 8235-11135 | 2900 | 645/12/1 (99 %) |
| SDSS J1326+4806 | 2.078 / 0.396 | S, N (G 20.7, 20.9; 21.06" apart, 6.1" SE and 15.0" NW of the BCG) | 2120 | 8199-11077 | 2878 | 757/12/4 (96 %) |
| COOL J0542-2125 | 1.84 / 0.61 | A+B (2.4" apart, blended), C (g 20.7 brightest) | 696 | 8377-11103 | 2726 | 256/13/4 (83 %) |
| COOL J0335-1927 | 3.27 / 0.4178 | A, B, C | 1022 | 8332-11284 | 2952 | 301/1/2 (98 %) |
| COOL J1153+0755 A | 1.524 / 0.4301 | A1, A2, A3 (not point sources in DECaLS; faint) | 968 | 8203-11220 | 3017 | 472/15/1 (87 %) |

Not in the list: SDSS J0909+4449 (3 images, 14" ; Shu et al. 2018) was not in the brief and has no pre-published delay
prediction found; it is a loose end. The GraL compilation's delay table (arXiv:2603.12810) holds galaxy-scale systems
only besides J2222 and J1651 (both commented out of the table in the source); no further wide system without a ZTF
delay came from it.

Pairs, with the published delay (to validate) or the pre-published prediction (to test), and whether the 2018-2026
baseline can reach it (overlap = baseline minus delay; a pair needs at least a season, ~150 d, of overlapping variable
data; said here before any fit). Convention: dt_XY = t_Y - t_X, positive when Y lags X.

| system | pair | published delay | prediction (pre-measurement) | overlap (d) | reachable? |
|---|---|---|---|---|---|
| J1004 | B->A | 43.0 ± 0.3 (Munoz+2022; B leads A) | | 2830 | yes (validation) |
| J1004 | C->B | 782.2 ± 0.4 | | 2090 | yes (validation) |
| J1004 | C->A | 825.2 ± 0.5 | | 2050 | yes (validation) |
| J1004 | C->D | 2458.5 ± 1.0 | | 415 | marginal: D's first 1.1 yr vs C's last 1.1 yr, D at r 20.5 |
| J1029 | A->B+C | 744 ± 10 (Fohlmeister+2013, 90 %; A leads) | | 2130 | yes (validation) |
| J2222 | A->B | 47.7 ± 6.0 (Dahle+2015, 95 %) | | 2850 | yes if the r ~ 21 images vary enough per season |
| J2222 | C->A | 722 ± 24 (Dahle+2015, 95 %) | | 2180 | same |
| J2222 | A->D,E,F | | 502 ± 68, 611 ± 75, 415 ± 72 (Sharon+2017) | | no: D, E, F are 2-4 mag below ZTF's per-epoch limit |
| J1326 | S<->N | none | none published (SIS: see below) | 2878 - |dt| | yes for |dt| < ~2500 |
| J0542 | A+B<->C | none | none numerical (arXiv:2209.03972 says "multiple years") | 2726 - |dt| | measurable, nothing to test |
| J0335 | A->B | none | 499 (+141/-146), syst. 240-700 (Napier+2023) | 2450 | yes |
| J0335 | A->C | none | -127 (+83/-17), syst. -300 to -30 | 2820 | yes |
| J1153 | A1->A2 | none | 1253 (+42/-50), space-based update 1146 (+107/-54) (Solhaug+2026) | 1760 | yes if the images are detected |
| J1153 | A1->A3 | none | 789 (+55/-59), update 788 (+48/-55) | 2230 | same |

The test, as fixed here: a prediction is "confirmed" when the measured delay's 1-sigma interval overlaps the prediction's
own interval (the paper's statistical interval; the systematic range where the paper gives one is reported beside it),
"rejected" when the two intervals are more than their quadrature-summed widths apart, and "untested" when ZTF cannot
measure the pair. The sample-wide statement is the count.

## Photometry (ztf_lens_phot.py) and its check

Each epoch's difference-image cutout (60", the pipeline `scimrefdiffimg`) is fitted with the pipeline's own
difference-image PSF stamp (`diffimgpsf`, 25 x 25) placed at the fixed image positions, one free flux per image plus a
constant background; a global astrometric shift is fitted when any image reaches S/N 8 in the epoch (it never did at
r > 19, so the ZTF WCS is taken as is). Difference flux to AB microJansky with the epoch's MAGZP; total flux = the
reference PSF-fit catalogue flux (the reference zero point) + difference flux; epochs with INFOBITS != 0, seeing > 4" or
reduced chi2 > 3 dropped; same-night epochs averaged; magnitudes with a 0.01 mag floor. Only the primary (field, ccd,
quadrant) per band is used, so each image has one reference.

Check against the public ZTF light-curve service (PSF-fit catalogue photometry on the science images, DR23, catflags 0),
J1004 r band, nightly means: image A (3.8" from B) service minus mine median -0.017 mag, rms 0.104, correlation 0.81 over
493 nights; image C (isolated) -0.008, rms 0.132, correlation 0.88 over 476 nights. The service's own single-epoch errors
at r 19.5-20 are 0.12-0.15 mag, so the two agree at their noise. Night-to-night rms of my nightly bins is 1.4-1.5 times
the formal error (e.g. J1004 r: A 0.18 vs 0.12 expected, C 0.21 vs 0.16), i.e. ~0.1 mag of per-night noise the
difference-image formal errors miss (resampling correlations, zero-point rms 0.03); the delay mocks carry it as an excess
term (below), so the quoted delay errors include it.

## Delay estimator (ztf_lens_delays.py) and the validation on J1004 and J1029

For a pair (X, Y) and a trial delay tau, Y's epochs are shifted back by tau and both curves are fitted jointly by one
cubic B-spline (knots every 30 d, a second-difference ridge so season gaps stay tame) for the intrinsic variability,
plus a slow spline in Y's own time (knots every 365 d) for microlensing and the magnitude offset. chi2(tau) on a grid
(1 d for short delays, 2-4 d for long ones); the delay is the parabolic minimum. Errors from 30 mock pairs per fit: the
best model at a true delay drawn within a few steps of the minimum, plus the formal noise and the measured excess
residual rms (0.10-0.13 mag) as exponentially correlated noise (10 d), refitted on a window around the minimum; the
error is the quadrature sum of the mocks' mean bias and scatter (the PyCS convention). This is deliberately
conservative: the chi2 profile's own delta-chi2 = 1 interval is 1-3 d wide for the J1004 pairs, the mock error 11-28 d.

| system | pair | band | ZTF delay (d) | published (d) | difference |
|---|---|---|---|---|---|
| J1004 | C->A (14.6") | r | 825.8 ± 28.2 | 825.2 ± 0.5 (Munoz+2022) | +0.6 (0.0 sigma) |
| J1004 | C->A | g | 848.1 ± 19.1 | 825.2 | +23 (1.2 sigma) |
| J1004 | C->B (14.6") | r | 782.6 ± 11.2 | 782.2 ± 0.4 | +0.4 (0.0 sigma) |
| J1004 | C->B | g | 785.8 ± 18.7 | 782.2 | +3.6 (0.2 sigma) |
| J1004 | B->A (3.8", blended) | r | 53.7 ± 13.0 | 43.0 ± 0.3 | +10.7 (0.8 sigma) |
| J1004 | B->A | g | 61.8 ± 13.6 | 43.0 | +18.8 (1.4 sigma) |
| J1004 | C->D (12") | r | not recovered: the minimum runs to the grid's end (2703, overlap 170 d); secondary minima 2376, 2396; the published 2458.5 ± 1.0 leaves 415 d of overlap, D's first season at r 20.5 against C's last, noisy one | 2458.5 ± 1.0 | untested by this baseline |

The three resolved pairs come back within 1.2 sigma, the two r-band long delays within 1 d of the published values (the
2.5-mag rise of 2018-2021 in A and B, seen in C from 2020, is what fixes them: `delay_CA_zr.png`). The blended 3.8"
pair is 11-19 d high in both bands, with every setting tried (knot 15-50 d, polynomial or spline microlensing: 54-77 d);
a bias of the scene model on pairs closer than about two seeing widths, where PSF errors move flux between the images,
or short-timescale microlensing in A (known in this system) that the 365-d spline cannot follow. Pairs closer than 4"
are therefore not claimed below (J1004 A-B, J1029 B-C, COOL J0542 A-B).

## Results (the overlap-scored estimator, 30 mocks; every number also in the note in `finding.md`)

The first runs scored chi2 over all points; a long trial delay then won by leaving little to compare (J1326 ran to
-2275 d with a flat profile). The statistic was changed to chi2 per point over the overlap of the shifted curves with a
25 % overlap floor; the J1004 and J1029 validation numbers below are from that version. Faint images (r > 21) were
binned over 5 d in flux before the magnitude conversion, after a nightly S/N cut was found to truncate their noise
(J2222 C's first fits).

Setting dependence (r band; knot 20 / 50 d, microlensing knot 200 / 600 d, bin 3 / 10 d, ridge 0.3 / 3; `logs/robust/`):

| pair | default | knot20 | knot50 | ml200 | ml600 | bin3 | bin10 | ridge0.3 | ridge3 |
|---|---|---|---|---|---|---|---|---|---|
| J1029 A->B | 706.5 | 724.9 | 710.5 | 731.3 | 725.5 | 708.3 | 707.6 | 708.9 | 709.2 |
| J2222 C->A | 884.7 | 880.2 | 698.2 | 884.7 | 884.8 | 794.0 | 462.1 | 884.5 | 884.8 |
| J0335 A->B | 900 (edge) | 553.1 | 486.4 | 900 (edge) | 560.4 | 561.3 | 900 (edge) | 561.5 | 900 (edge) |
| J0335 A->C | 225.0 | 226.8 | -462.8 | 225.4 | 225.5 | 221.6 | -441.2 | 227.1 | -442.8 |
| J1153 A1->A2 | 1769 | 1800 (edge) | 1789 | 1800 (edge) | 1734 | 1800 (edge) | 1582 | 1800 (edge) | 1769 |
| J1153 A1->A3 | 914.7 | 905.5 | 913.0 | 914.8 | 914.2 | 943.4 | 918.5 | 914.2 | 916.6 |
| J1326 A->B (12 % floor) | -1928.6 | -1924.0 | -1949.2 | -1943.3 | -1942.7 | -1974.9 | -1921.4 | poly2: -1974.7 | |

J1326 in g with the 12 % floor: +2487.7 (the other end of the scan), other minima -2316 to -1980; both bands have 0.2 mag
of excess scatter and the central |tau| < 1500 d region is worse by delta-chi2 > 100 in both, which is the signature of
noise-dominated curves meeting at their quiet ends, not of a delay. COOL J0542 image C has no source in either ZTF
reference catalogue (g 22.1 by the scene fit on the reference cutout; 13 usable g bins, none in r): no measurement.

## The table

| system | pair | ZTF delay (d) | settings spread | published / predicted (d) | status |
|---|---|---|---|---|---|
| SDSS J1004+4112 | C->A | r 825.7 ± 25.7; g 847.0 ± 21.6 | | 825.2 ± 0.5 measured | recovered |
| SDSS J1004+4112 | C->B | r 781.6 ± 21.0; g 786.0 ± 24.2 | | 782.2 ± 0.4 measured | recovered |
| SDSS J1004+4112 | B->A (3.8") | r 52.2 ± 14.6; g 60.7 ± 13.2 | 54-77 | 43.0 ± 0.3 measured | 9-18 d high; blended |
| SDSS J1004+4112 | C->D | overlap 14 % of the baseline | | 2458.5 ± 1.0 measured | untested |
| SDSS J1029+2623 | A->B(+C) | r 706.5 ± 3.7; g 739.9 ± 15.8 | r 707-731 | 744 ± 10 (90 %) measured | g recovered; r 1-2 sigma low |
| SDSS J2222+2745 | A->B | r 42.6 ± 5.8; g unusable | r 33-50 | 47.7 ± 6.0 (95 %) measured | recovered |
| SDSS J2222+2745 | C->A | r 884.7 ± 17.5; g 461 ± 21 | r 462-885 | 722 ± 24 (95 %) measured | not recovered (C at r 21.2) |
| SDSS J2222+2745 | A->D, E, F | images at r 23-25 | | 502 ± 68, 611 ± 75, 415 ± 72 predicted | untestable |
| SDSS J1326+4806 | A->B | none (r -1930, g +2490: opposite scan ends) | | none | no measurement |
| COOL J0542-2125 | AB->C | C not in the ZTF reference | | none numerical | no measurement |
| COOL J0335-1927 | A->B | r 486-561 (5/8) or 900 (3/8); g 636 ± 153 | | 499 (+141/-146), syst 240-700 predicted | consistent, unstable |
| COOL J0335-1927 | A->C | r +225 ± 20 (5/8) or -445 ± 25 (3/8); g +220 ± 36 | | -127 (+83/-17), syst -300 to -30 predicted | tension (sign) |
| COOL J1153+0755 | A1->A2 | runs to the baseline's end (1770-1800) | | 1253 (+42/-50) predicted | unconstrained |
| COOL J1153+0755 | A1->A3 | r 914.7 ± 10.8 | 905-943 | 789 (+55/-59); 788 (+48/-55) predicted | 2.1 sigma high |


## Finding check (`blind-review.md`, context-free Opus, 2026-11-01, the note alone in a file)

VERDICT: novel but not sound: the single new delay (J1153 A1→A3) and the claimed tensions sit in the faint regime where the method misplaced J2222 C→A by roughly 160–190 d, and the error model only measures scatter around the chosen minimum, so neither the ±10.8 d error nor the 2.1σ tension is supported.

Strongest objection, in one sentence: the one new delay and both tensions sit in the brightness regime where the pipeline's only test against a published delay failed by 7-10 of its quoted sigma, and mocks built around the chosen minimum cannot say whether the mode itself is an alias of the seasonal sampling; the referee asked for the same scan on faint mock curves with the real sampling and noise and a rate of catastrophic errors.

Injection-recovery, run after the verdict (the test the referee named; `ztf_lens_delays.py --inject`, 40 mocks per
case from the model fitted at the injected delay, the data's own sampling, formal noise and excess as 10-d correlated
noise, scanned over the full grid; `logs/inject/`, `data/ztf_lens/<target>/inject_*.json`):

| pair (r) | injected (d) | found within 40 d | where the misses land | the data's minimum |
|---|---|---|---|---|
| J1004 C->A (bright control) | 825 | 35/40 (825 ± 13) | 784, 873-894, 999 | 825.7 |
| J2222 C->A | 722 (published) | 25/40 (723 ± 10) | 7 of 15 at 876-920, the data's mode; others 482-635, 941-1000 | 884.7 |
| J1153 A1->A3 | 789 (predicted) | 30/40 (788 ± 5) | 5 of 10 at 880-916, the data's mode; 430, 957, 1152-1173 | 914.7 |
| J1153 A1->A3 | 915 (measured) | 40/40 (914 ± 8) | none | 914.7 |
| J0335 A->C | -127 (predicted) | 9/40 (-144 ± 17) | 18 of 31 at +116 to +234, the data's mode; 11 at -178 to -499 | +225.0 |
| J0335 A->C | +225 (measured) | 20/40 (224 ± 14) | -434 to -113, +132 to +182 | +225.0 |

So the catastrophic-error rate of the estimator in the faint regime is 12 % for J1153 at the predicted delay (with
half of those errors landing on the measured 915), 38 % for J2222 at its published delay (half of them on the data's
885), and 78 % for J0335 A->C at its predicted delay, with the measured +225 the single most likely wrong answer. The
J2222 C->A miss and the J0335 A->C "sign tension" are what this estimator does to a true 722 and a true -127 with these
data; they are not evidence against the published or predicted values. J1153 A1->A3 = 915 is recovered every time when
true and produced by a true 789 about one mock in eight: a candidate delay, not a measurement that tests the prediction.
Parked at the referee's verdict (a physics point, not a presentation one); the next pass, if anyone takes it, is the
note rewritten around this table with no tension claimed and the 915 d candidate carrying its 12 % alias probability,
then a second referee pass. Attribution corrected after the referee: the COOL J0335 paper (arXiv:2305.14317) is Napier et
al. 2023 and the COOL J0542 paper (arXiv:2209.03972) Martinez et al. 2023, swapped above and in the note; and Dyrland
2019's reanalysis of J2222 (dt_AB = 42.44 +1.36/-1.44, dt_AC = -696.65 +2.0/-2.1) is the later published value, which the
ZTF A->B 42.6 ± 5.8 matches.
