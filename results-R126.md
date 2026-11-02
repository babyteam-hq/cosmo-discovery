# Results (R-126, cosmo-deep-6, 2026-11-01)

Capture-recapture across the gravitational-wave search catalogues: an injection-free count of the detectable mergers in
O1-O4a from the overlap pattern of the independent searches on the same strain. This file is the worker's own. Scripts:
`gw_lists.py` (builds the per-run candidate lists from the catalogue papers' tables), `gw_capture_recapture.py` (the
estimators, the mocks, the injection-based VT ratios); data under `data/gw/` (gitignored: the arXiv LaTeX sources and
PDFs of the catalogue papers, the 4-OGC search output, the assembled `lists_<run>.csv`, the run logs `cr_data.log` and
`cr_mock.log`). Both scripts run in `data/gw/venv` (numpy, scipy, h5py, pandas, statsmodels), a few minutes in all.

## Trajectory check (run in the dispatching session before this worker started; verbatim)

"VERDICT: novel: no prior work answers this question". Closest prior work named there: the unified p_astro paper (PRD 108,
083043, 2023; injection-based), Malz et al. arXiv:2607.07272 (ML ranking of sub-threshold candidates by multi-pipeline
consensus, GWTC-5.0), arXiv:2401.08709 and arXiv:2312.06631 (catalogue overlaps described qualitatively), arXiv:2303.17628
(union list), the GWTC-5.0 methods paper arXiv:2605.27224 and the selection-function emulator arXiv:2408.16828
(injection-based). The referee's method caveat (shared strain and templates make the lists positively dependent through
SNR; closed-population estimators are biased low; interaction terms cannot identify a dependence through one shared
variable) is answered by construction below: the headline estimator conditions on SNR and chirp mass, residual dependence
is measured, and the estimate is presented as a lower bound where it is one.

## Reachability from the cloud (2026-11-01)

gwosc.org fails TLS at the proxy ("unable to get local issuer certificate"; www.gw-openscience.org resets), so the GWOSC
event lists with their per-pipeline columns and the O3/O4a injection sets (GWOSC and Zenodo, which 403s) were not used.
Everything came from the papers themselves: arXiv e-prints (LaTeX sources with the tables' macro files for GWTC-2,
GWTC-2.1, GWTC-3, 4-OGC, IAS O1, O2, O3a, O3b and IAS-HM; PDFs for GWTC-1 and the GWTC-4.0 papers, read with pdftotext)
and raw.githubusercontent.com (the 4-OGC search output `4-OGC_small.hdf`, 61 MB, every candidate with IFAR > 1e-5 yr).
dcc.ligo.org and gracedb.ligo.org (public superevents, O4 low-latency triggers per pipeline) answer 200 and were not needed.

## The lists

One row per candidate; a candidate is "in" a list when that search reports it above the threshold. Two thresholds are
carried throughout: FAR < 1/yr in that search (the LVK's "significant" line; for IAS the per-bank IFAR as published,
which carries no trials factor, and for 4-OGC the catalogue IFAR), and p_astro > 0.5 in that search. The per-pipeline
columns exist only for catalogued candidates, so a signal below every list's threshold is unobserved by construction;
that is what the estimators infer.

| run | lists (independent code bases in bold) | source of the per-list numbers | candidates in the union (FAR < 1/yr / p_astro > 0.5) |
|---|---|---|---|
| O1-O2 | **PyCBC**, **GstLAL**, **cWB** (LVK), 4-OGC (PyCBC-based, Nitz et al.), **IAS** | GWTC-1 Tables I-II (PDF); 4-OGC HDF; IAS arXiv:1902.10341 (O1), 1904.07214 (O2) | 17 HL BBH / (GWTC-1 gives no per-pipeline p_astro) |
| O3a, O3b | cWB, **GstLAL**, **MBTA**, **PyCBC** (broad), PyCBC-BBH (LVK), 4-OGC, **IAS** (the (2,2) O3a and O3b searches and IAS-HM, one list, best IFAR / p_astro) | GWTC-3 and GWTC-2.1 event macros (`event_macros.tex`, `marginal_macros.tex`, the PE macros for chirp masses); GWTC-2 for cWB in O3a; 4-OGC HDF; IAS arXiv:2201.02252, 2311.06061, 2312.06631 (new-event and cross-catalogue tables) | 62 / 97 HL BBH; 69 / 82 all candidates, LVK lists |
| O4a | cWB-BBH, GstLAL, MBTA, PyCBC (LVK) | GWTC-4.0 Tables 1, 2 (per-pipeline FAR, SNR, p_astro for the 129 candidates with p_astro > 0.5 in at least one pipeline) and 3 (chirp masses) | 87 / 128 |

Shared templates and code: PyCBC-BBH is PyCBC on a BBH-focused bank; 4-OGC is the PyCBC code base run by an external
group with its own bank, ranking statistic and data selection; GstLAL, MBTA and PyCBC share the LVK template-bank
construction methods but not the banks or ranking; cWB is template-free; IAS is an independent code base, bank and
statistic (Hanford-Livingston coincident data only, BBH banks only). For the runs with external lists the population is
therefore restricted to Hanford-Livingston coincident BBH candidates (both detectors observing, chirp mass > 5 Msun or
p_BBH > 0.5), the population every list searched; the LVK-only analyses use every candidate.

The reported SNR of a candidate is the maximum network SNR among the lists that report it (the IAS network SNR where it
is the only list). Chirp masses are detector-frame: the LVK PE medians (GWTC-2.1 and GWTC-3 PE macros), GWTC-4.0 Table 3's
source-frame value times (1+z), the 4-OGC template, or the IAS PE value times (1+z), in that order of preference.

## The mock check (`gw_capture_recapture.py --mock`; `data/gw/cr_mock.log`, `cr_mock_oracle.log`)

Populations of 20,000 (or 3,000, configuration F) signals with true network SNR distributed as SNR^-4 above 5 (a
Euclidean volume), a shared noise term on top (sigma 0.6, the common strain realisation; zero in configuration B) giving
each signal a base SNR, and for each of 5 or 7 lists its own noise (sigma 0.3-0.5), a threshold between 8.0 and 9.3, and
a random drop (5 % or 20 %: vetoes, gaps). Configuration D gives one list a chirp-mass-dependent threshold (cWB-like).
The estimators see the union and each candidate's reported SNR (the maximum over the lists that found it, as in the
data). Truth per bin is the number of signals, found or not, whose base SNR lies in the bin.

| configuration (3 realisations each; 5 for F) | union | Chao lower bound | log-linear independence / pairwise | HT given reported SNR and Mc | MH odds ratio | reported SNR 8-8.5: seen / HT / truth | 8.5-9.5 | 9.5-10.5 | > 10.5 |
|---|---|---|---|---|---|---|---|---|---|
| A: 5 lists, shared noise 0.6, own 0.3, thresholds 8.0-9.0, 5 % drops | 5217-5261 | 5714-5760 ± 52 | 5218-5262 / 5445-5488 | 5227-5272 (± 1) | 1.8-2.0 | 549-582 / 554-587 / 850-899 | 1347-1370 / 1351-1375 / 1188-1250 | 887-944 / same / 754-835 | 2371-2433 / same / 2121-2215 |
| B: as A, no shared noise (own 0.5) | 5348-5479 | 5946-6155 | 5352-5483 / 6070-6275 | 5381-5521 | 1.4 | 481-507 / 500-526 / 798-816 | 1404-1441 / 1423-1462 / 1180-1193 | 934-997 / 936-1000 / 733-757 | 2493-2588 / same / 2082-2187 |
| C: as A, 20 % drops | 5062-5121 | 5370-5503 | 5077-5138 / 5164-5246 | 5086-5148 | 1.0 | 477-566 / 482-572 / 811-863 | 1236-1286 / 1247-1297 / 1132-1208 | 878-897 / 883-902 / 763-802 | 2379-2469 / 2383-2473 / 2185-2230 |
| D: as A, one list needs 2 more SNR per e-fold of lower chirp mass | 5544-5604 | 6287-6409 | 5549-5608 / 6299-6431 | 5581-5631 | 1.3-1.5 | 583-603 / 592-616 / 874-894 | 1370-1406 / 1378-1414 / 1173-1227 | 844-893 / same / 752-811 | 2353-2443 / same / 2153-2219 |
| E: 7 lists, thresholds 8.0-9.3 (the O3 shape) | 5302-5398 | 5756-5837 | 5302-5398 / 5450-5551 | 5308-5403 | 1.7-1.9 | 571-646 / 574-650 / 861-882 | 1278-1307 / 1280-1309 / 1188-1209 | 895-958 / same / 831-853 | 2483-2563 / same / 2165-2246 |
| F: as A, 3,000 signals (catalogue-sized) | 766-822 | 826-915 ± 17-24 | 766-822 / 806-846 | 767-824 | 1.6-1.8 | 81-98 / 82-99 / 120-143 | 183-212 / 184-213 / 167-199 | 125-160 / same / 114-137 | 330-389 / same / 295-350 |

What the mocks show. (1) Every estimator returns the union to within 0-15 % (HT given the reported SNR: +0.1 to +0.8 %;
the independence log-linear model: +0.1 %; the pairwise log-linear model and Chao's bound: +4 to +15 %), and the
Mantel-Haenszel odds ratio stays above 1 (1.4-2.0) in every configuration with shared noise and even in B without it
(1.4): the selection on the maximum of noisy statistics is itself a positive dependence given the reported SNR. (2) The
truth in the lowest reported-SNR bin is 1.5 times the union (899 vs 549 in A; 136 vs 117 in F), and no estimator sees it:
the signals missed by every list at base SNR 8-8.5 are the ones whose noise fell low in every list, and the
reported SNR of the found ones is biased high by the same selection. Above a reported SNR of 9.5 the union exceeds the
truth per bin (upward scatter of the reported maximum), so the HT per bin is right in the sense that nothing is missing
there. (3) The oracle (`cr_mock_oracle.log`): the same estimator conditioned on the true base SNR recovers the truth in
every bin from 8.5 up (A: 1189 vs 1188, 754 vs 754, 2215 vs 2215) and 80-90 % of it in the 8-8.5 bin (780 of 899; F: 118
of 136). The method is sound when the conditioning variable is the shared variable itself, and returns the union when
the conditioning variable is the selected maximum. The catalogues give only the latter (the PE network SNR of the
confident events is a less biased measurement, but the marginal candidates, which carry the information, have none).
(4) The nested null: a list that is a random 60 % subsample of another returns the larger list exactly (Lincoln-Petersen
4902 ± 0 for 4902), as it must.

Consequence for the data: the capture-recapture totals below are lower bounds on the number of detectable signals, equal
to the union to within a few per cent, and the near-threshold incompleteness that the injection campaigns measure is
invisible to the overlap pattern of lists that share the strain. The quantity the overlap does measure without bias is the
relative catchability of the lists on the same candidates, which is what the injection-based relative sensitive volumes
are compared with below.

## Estimators (`gw_capture_recapture.py`; full output in `data/gw/cr_data.log`, 200 bootstrap replicates)

On the capture-history matrix of each sample (one row per candidate in the union, one column per list): Lincoln-Petersen
(Chapman) per pair; Chao's (1987) lower bound from the candidates seen by exactly one and exactly two lists; the
sample-coverage estimate (Chao & Tsay 1998, independence form); Poisson log-linear models on the 2^K − 1 observed
cells with main effects and with all pairwise interactions (the unobserved cell extrapolated; AIC); and the headline, a
Huggins-Alho conditional-likelihood model in which each list's catchability is logistic in the reported network SNR and
log chirp mass, the lists conditionally independent given both, and the total the Horvitz-Thompson sum over the observed
candidates of 1 / P(found by at least one list) ("HT"), bootstrapped over candidates (68 % intervals). "N_w" is the same
sum with each candidate weighted by its probability of being real (the maximum p_astro any list assigns it), so that the
noise candidates a p_astro > 0.5 list admits do not count as missed signals; "expected real" is the sum of those weights.
Residual dependence given SNR: the Mantel-Haenszel odds ratio of each pair of lists pooled over five reported-SNR strata
(1 = conditionally independent; > 1 = positive dependence left, so the HT total is a lower bound).

| sample | threshold | union (expected real) | LVK union (real) | per list | Chao lower bound | log-linear indep. / pairwise | HT given SNR, Mc | N_w (p_astro-weighted) | MH odds ratio (median) | LVK share of N_w |
|---|---|---|---|---|---|---|---|---|---|---|
| O1-O2, HL BBH, 5 lists | FAR < 1/yr | 17 (16.0) | 10 (10.0) | PyCBC 8, GstLAL 10, cWB 7, 4-OGC 12, IAS 15 | 21.2 ± 4.9 | 17.1 / 17.0 | 17.3 (17.0-18.6) | 16.3 (15.9-17.3) | 2.4 (1 pair usable) | 0.61 |
| O1-O2, all, 3 LVK lists | FAR < 1/yr | 11 | 11 | PyCBC 9, GstLAL 11, cWB 7 | 11.1 ± 0.4 | 11.0 / 11.0 | 11.0 (11.0-11.0) | 11.0 | – | 1.00 |
| O3, HL BBH, 7 lists | FAR < 1/yr | 62 (60.7) | 51 (50.4) | cWB 24, GstLAL 44, MBTA 32, PyCBC 34, PyCBC-BBH 44, 4-OGC 44, IAS 59 | 67.6 ± 4.6 | 62.0 / 64.5 ± 4.2 | 62.4 (62.1-63.0) | **61.1 (60.7-61.7)** | 1.24 | **0.83** |
| O3a, HL BBH, 7 lists | FAR < 1/yr | 37 (36.4) | 31 (30.7) | cWB 15, GstLAL 26, MBTA 21, PyCBC 21, PyCBC-BBH 27, 4-OGC 25, IAS 36 | 40.6 ± 3.9 | – | 37.3 (37.0-37.6) | 36.7 (36.4-37.1) | 2.4 | 0.84 |
| O3b, HL BBH, 7 lists | FAR < 1/yr | 25 (24.3) | 20 (19.7) | cWB 9, GstLAL 18, MBTA 11, PyCBC 13, PyCBC-BBH 17, 4-OGC 19, IAS 23 | 27.0 ± 2.6 | – | 25.4 (25.0-28.4) | 24.6 (24.3-25.9) | 1.0 | 0.80 |
| O3, HL BBH, 7 lists | p_astro > 0.5 | 97 (84.2) | 65 (61.9) | cWB 25, GstLAL 44, MBTA 46, PyCBC 37, PyCBC-BBH 58, 4-OGC 61, IAS 80 | 497 ± 310 (40 singletons, 2 doubletons: not usable) | 97.3 / diverges | 105.5 (102.2-115.5) | **88.5 (87.3-94.0)** | 3.15 | **0.70** |
| O3a, HL BBH, 7 lists | p_astro > 0.5 | 52 (45.9) | 36 (34.9) | cWB 15, GstLAL 27, MBTA 28, PyCBC 22, PyCBC-BBH 33, 4-OGC 34, IAS 48 | 214 ± 180 | – | 54.2 (52.9-57.1) | 47.0 (46.2-48.6) | 4.3 | 0.74 |
| O3b, HL BBH, 7 lists | p_astro > 0.5 | 45 (38.3) | 29 (27.0) | cWB 10, GstLAL 17, MBTA 18, PyCBC 15, PyCBC-BBH 25, 4-OGC 27, IAS 32 | 287 ± 264 | – | 51.3 (49.5-62.9) | 41.6 (40.7-50.7) | 1.9 | 0.65 |
| O3, all candidates, 5 LVK lists | FAR < 1/yr | 69 (64.1) | 69 | cWB 26, GstLAL 58, MBTA 36, PyCBC 37, PyCBC-BBH 47 | 147 ± 51 | 69.5 / 256 ± 173 | 70.7 (70.0-72.7) | 65.1 (63.9-67.2) | 3.3 | 0.98 |
| O3, all candidates, 5 LVK lists | p_astro > 0.5 | 82 (74.7) | 82 | cWB 28, GstLAL 55, MBTA 49, PyCBC 42, PyCBC-BBH 62 | 129 ± 24 | 82.9 / 288 ± 138 | 84.8 (83.8-88.5) | 76.3 (75.3-77.7) | 3.3 | 0.98 |
| O4a, all candidates, 4 LVK lists | FAR < 1/yr | 87 (85.5) | 87 | cWB 44, GstLAL 77, MBTA 50, PyCBC 70 | 97.8 ± 6.7 | 87.4 / diverges | 88.4 (87.7-90.5) | **86.8 (86.3-89.0)** | 2.3 | 0.99 |
| O4a, all candidates, 4 LVK lists | p_astro > 0.5 | 128 (115.2) | 128 | cWB 54, GstLAL 87, MBTA 75, PyCBC 102 | 179 ± 21 | 130.3 / 817 ± 521 | 130.9 (129.9-133.9) | **117.0 (115.7-119.6)** | 1.8 | 0.98 |
| O4a, BBH, 4 LVK lists | p_astro > 0.5 | 83 (81.7) | 83 | cWB 47, GstLAL 74, MBTA 62, PyCBC 74 | 86.7 ± 3.3 | – | 83.2 (83.1-83.6) | 81.8 (81.7-82.1) | 3.2 | 1.00 |

By reported SNR (O3 HL BBH, p_astro > 0.5): SNR < 8.5: 25 seen (18.2 expected real), 7 of them LVK, N_w 21.5
(17.5-27.8); 8.5-9.5: 27 (21.7), 14 LVK, N_w 22.5 (19.1-26.7); 9.5-10.5: 10, all LVK, N_w 9.7; above 10.5: 35, 34 LVK,
N_w 34.7. At FAR < 1/yr the same sample: SNR < 8.5: 7 seen, 1 LVK, N_w 6.7 (4.8-9.4); 8.5-9.5: 12, 7 LVK, N_w 11.4
(8.3-14.5); above 9.5: 43, all LVK. O4a (FAR < 1/yr): SNR < 8.5: 2; 8.5-9.5: 19, N_w 19.5 (16.1-24.7); 9.5-10.5: 20, N_w
19.6; every bin's N_w within 3 % of the seen.

By chirp mass (detector frame; O3 HL BBH, p_astro > 0.5; seen (expected real) / LVK (real) / N_w (68 %) / Chao):
Mc 5-15: 22 (18.9) / 14 (13.6) / 19.2 (15.6-23.3) / 58 ± 26; Mc 15-30: 10 (9.0) / 9 (8.4) / 9.2 (6.5-11.5) / 13.0 ± 4.5;
Mc 30-60: 39 (36.0) / 31 (30.1) / 37.8 (33.8-43.9) / 80 ± 49; Mc > 60: 26 (20.4) / 11 (9.8) / 22.3 (18.5-27.7) / 207 ± 199.
The LVK catalogue's share of the expected real signals is therefore 0.71, 0.91, 0.80 and 0.44 in the four bins; the
heaviest bin is where the IAS higher-harmonic search (arXiv:2312.06631) reports 21 of the 26 candidates, with p_astro
0.5-0.9 from its own rate model. At FAR < 1/yr: Mc 5-15: 14 (13.6) / 12 / N_w 13.7; 15-30: 7 / 6 / 7.0; 30-60: 32 (31.3) /
26 / 31.6 (27.5-34.7); > 60: 9 / 7 / 8.8. O4a (FAR < 1/yr): Mc 5-15: 16 / N_w 16.0; 15-30: 13 / 13.0; 30-60: 35 / 35.0
(30.8-41.0); > 60: 21 / 20.7.

Catchability curves (the fitted logistic p_i at chirp mass 20, FAR < 1/yr). O3 LVK lists, SNR 8 / 9 / 10 / 12: cWB
0.08 / 0.11 / 0.14 / 0.21; GstLAL 0.41 / 0.54 / 0.67 / 0.86; MBTA 0.23 / 0.28 / 0.35 / 0.48; PyCBC 0.16 / 0.23 / 0.31 /
0.50; PyCBC-BBH 0.51 / 0.55 / 0.58 / 0.64; any 0.83 / 0.90 / 0.95 / 0.99. O4a: cWB 0.15 / 0.17 / 0.20 / 0.27; GstLAL
0.13 / 0.49 / 0.87 / 1.00; MBTA 0.46 / 0.50 / 0.55 / 0.64; PyCBC 0.36 / 0.57 / 0.75 / 0.94; any 0.74 / 0.91 / 0.99 / 1.00.
In the seven-list O3 sample the IAS list's catchability is 0.88-0.95 at every SNR (its per-bank IFAR carries no trials
factor, so at a nominal 1/yr it is the most permissive list), which is why that union is complete to 0.99 on its own terms.

## Relative catchability vs the injection-based relative VT (`gw_capture_recapture.py --compare`)

The quantity the overlap measures without the threshold bias is each pipeline's share of the candidates the union holds.
Here, for the five LVK O3 pipelines at p_astro > 0.5 (the GWTC-3 VT table's threshold), the share of the expected real
signals (sum of max p_astro) per chirp-mass bin, with a binomial error, against GWTC-3's all-O3 sensitive hypervolume of
that pipeline divided by the "any pipeline" hypervolume at the nearest injected chirp mass (Table in `detection.tex`,
printed by `--vt`; 35+35 Msun: cWB 0.49, GstLAL 0.79, MBTA 0.66, PyCBC 0.62, PyCBC-BBH 0.82; 35+20: 0.42, 0.77, 0.61,
0.61, 0.83; 20+20: 0.33, 0.79, 0.68, 0.67, 0.84; 10+5: 0.12, 0.81, 0.74, 0.87, 0.84; 5+5: 0.07, 0.82, 0.62, 0.87, 0.68).

| Mc (det.) | candidates (expected real) | cWB | GstLAL | MBTA | PyCBC | PyCBC-BBH |
|---|---|---|---|---|---|---|
| 5-15 (VT at 10+5) | 17 (15.8) | 0.06 ± 0.06 (VT 0.12, −0.9σ) | 0.69 ± 0.12 (0.81, −1.0σ) | 0.74 ± 0.11 (0.74, 0.0σ) | 0.85 ± 0.09 (0.87, −0.2σ) | 0.89 ± 0.08 (0.84, +0.7σ) |
| 15-30 (VT at 35+20) | 13 (12.1) | 0.41 ± 0.14 (0.42, 0.0σ) | 0.74 ± 0.13 (0.77, −0.3σ) | 0.61 ± 0.14 (0.61, 0.0σ) | 0.50 ± 0.14 (0.61, −0.8σ) | 0.72 ± 0.13 (0.83, −0.8σ) |
| 30-60 (VT at 35+35) | 34 (31.2) | 0.41 ± 0.09 (0.49, −0.9σ) | 0.79 ± 0.07 (0.79, 0.0σ) | 0.70 ± 0.08 (0.66, +0.5σ) | 0.48 ± 0.09 (0.62, −1.6σ) | 0.81 ± 0.07 (0.82, −0.2σ) |
| > 60 (VT at 35+35, extrapolated) | 11 (9.7) | 0.61 ± 0.16 (0.49, +0.8σ) | 0.61 ± 0.16 (0.79, −1.1σ) | 0.51 ± 0.16 (0.66, −0.9σ) | 0.51 ± 0.16 (0.62, −0.7σ) | 1.00 ± 0.09 (0.82, +2.1σ) |
| all BBH (5-1000) | 75 (68.8) | 0.36 ± 0.06 | 0.73 ± 0.05 | 0.67 ± 0.06 | 0.57 ± 0.06 | 0.84 ± 0.04 |

Twenty comparisons in the four BBH bins: 18 within 1.1σ, PyCBC (broad bank) at Mc 30-60 at −1.6σ and PyCBC-BBH above
Mc 60 at +2.1σ (11 of 11 candidates, against an injection point at a lower mass than the bin). Chi-square 14.9 for 20
(the injection-table values taken without their own 3-10 % errors). The injection-based relative sensitive volumes of the
five pipelines are consistent with their relative catchability on the real catalogue at the 0.05-0.15 level per bin. The
same shares at FAR < 1/yr (not the VT table's threshold; `--compare` prints them) move GstLAL up to 0.86 ± 0.05 and MBTA
down to 0.58 ± 0.06 over all BBH, as expected from the two thresholds' different meaning per pipeline. GWTC-2.1's O3a BBH
population VTs (FAR < 2/yr; `--vt`): GstLAL 0.85, MBTA 0.61, PyCBC 0.63, PyCBC-BBH 0.83 of "any", against the O3a
seven-list sample's LVK shares at FAR < 1/yr of 26/31, 21/31, 21/31, 27/31 = 0.84, 0.68, 0.68, 0.87 (± 0.07-0.08).

## Finding check (the eighth section of `finding.md`, alone in a file, to a context-free Opus referee)

Verbatim: "VERDICT: novel but not sound: the injection-free total equals the union by construction (it conditions on a capture-defined SNR), and its claimed lower-bound status is undercut by noise-inflated max-p_astro weighting across incompatible thresholds, leaving only a low-power and not like-for-like relative-sensitivity check."

MISS on a method point; parked, no revision pass (the gate rule: one pass on a presentation point, park on a physics or
method point). The referee grants the novelty (no prior capture-recapture across gravitational-wave searches; closest
Laycock arXiv:1701.03801 on X-ray transients, the pipeline tallies arXiv:2303.17628 and arXiv:2401.08709, the
injection-free sensitivity work of Essick arXiv:2307.02765, 2606.14045, and the emulator arXiv:2408.16828) and objects
that (1) conditioning on the reported SNR, which exists only for captured candidates and is the selected maximum, makes
the total the union by construction, as the mocks here show; (2) the max-p_astro weights come from incompatible rate
models and thresholds, so the "lower bound" is not guaranteed for real signals (the Chao 497 ± 310 from 40 singletons is
the noise contamination showing); (3) the ratio check's 20 comparisons share candidates and the "any" denominator (the
chi-square is not for 20 degrees of freedom), compare population bins with fixed-mass injection points, test only
ratios to "any", and pass over the two largest deviations. Path to pass: a covariate that is not capture-defined (the PE
network SNR for all candidates, marginal ones included, which the catalogues do not publish) or lists with independent
noise (different detectors or data); a population-weighted VT per bin and a candidate bootstrap for the joint
chi-square. Headline of the request: no gate pass; the method's limit (the overlap of searches on one strain cannot
see the near-threshold selection function) is established by the mocks and is the useful product; the relative
catchability vs VT agreement stands as a ±10 % consistency statement.
