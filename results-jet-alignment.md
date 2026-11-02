# Results (R-160, cosmo-deep-28, 2026-11-02)

Do the four quasar jet alignment regions (QJARs) of Mandarakas et al. 2021 (arXiv:2106.06546, A&A 653, A123) stay aligned
on VLBI sources that are not in the claim's sample? Proposal: `open-problems/round7/t1-jet-alignment.md`. Trajectory check
(2026-11-02, in the brief): "VERDICT: novel: no prior work answers this question".

**Outcome: STOPPED AT THE STEP 0 FEASIBILITY GATE** (about 35 min after start). Too few independent sources fall inside
the frozen regions. Even with every generous choice, the test reaches only 1.2-2.4 sigma at half the published alignment
strength. At the full published strength it reaches 2.4-4.8 sigma, and 2.4-3.5 sigma on the primary (spectroscopic-z)
sample. In every bracket the independent fraction is f = n'/n = 0.04-0.40, far below the proposal's own stop line of
f >= 0.35 summed over the regions. No jet position angle was measured, and no alignment statistic was computed. Step 1
(the method check) was not run, because the gate comes first.

Scripts: `jet_common.py` (frozen regions, parallel transport, axial statistics), `jet_sample.py` (source table:
positions, first-image dates, redshifts), `jet_feas.py` (forecast). The data are in `data/jet/` (gitignored).

## F0. Access (fetch times; `data/` started empty)

| file | source | size | time |
|---|---|---|---|
| 2106.06546, 2203.13750, 2002.00982 (PDF) | arxiv.org | 3 files | seconds |
| rfc_2026b_cat.txt (22,839 sources) | astrogeo.org/sol/rfc/rfc_2026b | 5 MB | 1 s |
| Astrogeo image index and per-source file listings (22,381 directories) | astrogeo.org/images/&lt;name&gt;/, 24-32 parallel requests | 22,381 small pages | 7.6 min (a few DNS failures, re-fetched; none left empty) |
| OCARS, revision 27-Oct-2026 (22,095 sources) | http://www.gaoran.ru/english/as/ac_vlbi/ocars.txt (the gao.spb.ru host in the proposal redirects to an https host that fails the TLS handshake through the proxy) | 4 MB | seconds |
| Plavin et al. 2022 Tables 1 and 2 (J/ApJS/260/4) | vizier.cfa.harvard.edu asu-tsv, with redirects followed (`curl -L`); the CDS hosts time out | 9,263 + 16,552 rows | 3 s |

Fetched but not used for anything here: the Plavin tables. No Mandarakas / Blinov PA table exists on VizieR (both
"not found"), so the claim's per-source PAs are not public.

## F1. The frozen template (from 2106.06546 Table 1, before looking)

| QJAR | RA | Dec | radius | z range | S.L. | mean PA (transported to centre) |
|---|---|---|---|---|---|---|
| 1 | 36.18 | +25.95 | 14.82 deg | 0.02-1.53 | 5.54 sigma | -16.9 deg |
| 2 | 168.75 | +47.54 | 18.11 deg | 0.36-1.89 | 6.74 sigma | -17.7 deg |
| 3 | 180.15 | +23.56 | 16.07 deg | 0.29-0.65 | 5.72 sigma | -25.6 deg |
| 4 | 201.83 | +50.97 | 18.27 deg | 0.07-0.35 | 5.38 sigma | +37.6 deg |

The claim used Astrogeo as it stood on 2020-06-10 (15,235 sources, 7,290 "good" visually confirmed jets, so a 48 %
good-jet rate) and OCARS redshifts. The per-region source counts n are not published.

## F2. The independent sample (`jet_sample.py`)

- **NEW** means every Astrogeo map of the source was observed after 2020-06-10, so the source cannot be in the claim
  sample. This is conservative: a source observed earlier but added to the database later is counted as old. Of the
  21,343 image directories with an RFC/OCARS position, 5,581 are new and 15,762 old (the claim had 15,235). The new
  sources' first epochs fall in 2020 (H2): 1,635; 2021: 1,556; 2022: 1,515; 2026: 868; 2023-24: 7. A further 1,038
  directories have no RFC or OCARS position (non-RFC names), and so no redshift either.
- **Redshift rule** (fixed before any count). The primary rule takes OCARS's fiducial z with flag 's', excluding
  origins b (Quaia, which is Gaia spectrophotometry) and c (Makarov), and falls back to OCARS's SDSS column.
  Extended (B) adds OCARS photometric, Quaia and Makarov z. Most generous (C) also adds OLD sources whose OCARS z comes
  from a post-2020 compilation (origin d DESI DR1, e, f, g): those had PAs in 2020 but no redshift, so they were outside
  the claim's 3D regions.

| QJAR | sources in cap (new) | new, spec z in range (A) | + photo-z (B) | + old with post-2020 z (C) | old with NED/SDSS z in range |
|---|---|---|---|---|---|
| 1 | 453 (84) | 7 | 20 | 32 | 97 |
| 2 | 721 (162) | 34 | 39 | 46 | 256 |
| 3 | 552 (121) | 4 | 7 | 7 | 59 |
| 4 | 761 (167) | 12 | 13 | 14 | 60 |
| all | | **57** | 79 | 99 | 472 |

These counts are before any jet detection. QJAR3 (0.29 < z < 0.65) and QJAR4 (0.07 < z < 0.35) are narrow,
low-redshift slices. The new VLBI sources are fainter blazar candidates, and few of them have a spectroscopic redshift
there.

## F3. Forecast (`jet_feas.py`; positions, dates and redshifts only)

The S-test D of 2106.06546 behaves like the mean resultant length R of the doubled angles. Under random PAs,
<D> = 0.886/sqrt(n) and sigma_D = 0.463/sqrt(n), so the published S.L. Z gives R = (0.886 + 0.463 Z)/sqrt(n). n is
bracketed by 0.40-0.55 x the old sources now in each region with a NED/SDSS z (the claim's 48 % good-jet rate; some of
today's NED redshifts postdate 2020). On n' new good PAs, the V-test against the published axis has
E[V] = sqrt(2 n') R c, with c = 0.94 for ~10 deg PA errors. The four regions are combined as Stouffer sum V_k / 2. The
good-jet rate on the new, fainter sources is bracketed by 0.30-0.48.

| sample | good-jet rate | claim n / old | f_k (QJAR1-4) | Z, published strength | Z, half strength |
|---|---|---|---|---|---|
| A | 0.30 | 0.40 | 0.05, 0.10, 0.05, 0.15 | 2.79 | 1.40 |
| A | 0.30 | 0.55 | 0.04, 0.07, 0.04, 0.11 | 2.38 | 1.19 |
| A | 0.48 | 0.40 | 0.09, 0.16, 0.08, 0.24 | 3.53 | 1.77 |
| A | 0.48 | 0.55 | 0.06, 0.12, 0.06, 0.17 | 3.01 | 1.51 |
| B | 0.30 | 0.40 | 0.15, 0.12, 0.09, 0.16 | 3.43 | 1.71 |
| B | 0.30 | 0.55 | 0.11, 0.08, 0.07, 0.12 | 2.92 | 1.46 |
| B | 0.48 | 0.40 | 0.25, 0.19, 0.15, 0.26 | 4.34 | 2.17 |
| B | 0.48 | 0.55 | 0.18, 0.14, 0.11, 0.19 | 3.70 | 1.85 |
| C | 0.30 | 0.40 | 0.25, 0.14, 0.09, 0.18 | 3.78 | 1.89 |
| C | 0.30 | 0.55 | 0.18, 0.10, 0.07, 0.13 | 3.22 | 1.61 |
| C | 0.48 | 0.40 | 0.40, 0.22, 0.15, 0.28 | 4.78 | 2.39 |
| C | 0.48 | 0.55 | 0.29, 0.16, 0.11, 0.20 | 4.08 | 2.04 |

**Gate reading.** The published S.L. comes from a grid search over sky position, five scales and 50 redshift cuts,
followed by a particle-swarm optimisation of centre, radius and z range, with no look-elsewhere correction. So the
published R is an upper bound (winner's curse), and the half-strength column is the realistic one. It is below 3 sigma
in all twelve brackets (1.2-2.4). The full-strength column straddles 3 sigma on the primary sample (2.4-3.5). A null at
these numbers could not exclude a real alignment at half the published strength, and a positive result could not
reach 3 sigma unless the claimed amplitude were real at full, optimised strength. **Stop**, as the brief and the
proposal's own line (f >= 0.35) require.

## What stands, and what would reopen it

- **What stands:** 57 new VLBI sources with a spectroscopic redshift fall inside the four frozen (RA, Dec, z) regions,
  out of 5,581 sources first imaged after the claim's cut-off. That caps any independent re-test of the QJARs at
  ~1.2-2.4 sigma for a half-strength effect. The proposal's estimate (n' of about 1,800-2,500 new PAs and f = 0.25-0.6)
  counted new sources over the whole sky, not inside the narrow redshift slices.
- **To reopen:** (1) Redshifts for the 534 new sources inside the caps (QJAR1-4: 84 / 162 / 121 / 167), e.g. from DESI
  DR2 or targeted spectroscopy. At a ~20 % in-slice yield that would roughly double sample A. (2) The claim authors'
  per-region n and source lists, which would remove the n bracket. (3) A non-independent but cheap check, not run here:
  the QJAR S-test recomputed on the claim-era sources with Plavin et al. 2022's uv-plane jet directions (9,220 AGN).
  That tests whether the alignment depends on the PA method; it does not test new sources.
