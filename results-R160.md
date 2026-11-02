# Results (R-159 / brief R-160, cosmo-deep-27, 2026-11-02)

Do the Faraday-derotated polarisation angles of radio quasars align at 150-500 Mpc comoving, re-tested with SPICE-RACS
DR2 (arXiv:2605.16917) against the JVAS/CLASS 8.4 GHz claim (Tiwari & Jain, arXiv:1201.5180; Pelgrims & Hutsemekers,
arXiv:1503.03482)? Proposal: `open-problems/round6/s3-pol-align.md`. Trajectory check (2026-11-02, in the brief): "partly
novel ... still open is a test at 150-500 Mpc with an independent instrument and broadband Faraday-derotated angles
(SPICE-RACS), plus a check of the same sources against JVAS".

**Outcome: STOPPED AT THE STEP 0 FEASIBILITY GATE** (about 25 min after start). The flat-spectrum quasar subset cannot
test the published amplitude at 3 sigma at any scale: 0.1-0.8 sigma under the neutral transfer reading, 0.0-0.1 sigma under
the conservative one. The deciding number is the same-source transfer: SPICE-RACS derotated angles carry rho = 0.40 of
the JVAS angle. No alignment statistic was computed on real angles. Step 1 (the method check) was not run, because the gate
comes first.

Scripts: `r160_extract.py` (streams the 5 GB gzip and keeps the polarised rows), `r160_qso.py` (quasar list),
`r160_jvas.py` (anchor), `r160_sample.py` (de-duplication, matches, spectral index), `r160_common.py` (parallel transport,
statistics), `r160_feas.py` (forecast). Data in `data/r160/` (gitignored). Log: `logs/r160_feas.log`.

## F0. Access (fetch times; `data/` started empty)

| file | source | size | time |
|---|---|---|---|
| spice-racs.dr2.fits.gz (9,294,225 rows x 125 cols) | data.csiro.au DAP collection 76100 (not on CASDA TAP, which has DR1 only) | 4.97 GB | 9 min 44 s |
| DR16Q_v4.fits | data.sdss.org/sas/dr16/eboss/qso/DR16Q | 1.7 GB | ~4.5 min |
| QSO_cat_iron_cumulative_v0.fits (DESI DR1 QSO) | data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/QSO/iron | 816 MB | ~2 min |
| milliquas.fits.zip (v8) | quasars.org | 42 MB | seconds |
| JVAS/CLASS polarisation, J/MNRAS/376/371 table2 (12,743 rows) | vizier.cfa.harvard.edu (the CDS hosts time out from here) | 1 MB | seconds |
| NVSS (S > 3 mJy) | vizier.cfa.harvard.edu VIII/65 | 46 MB | ~2 min |
| RACS-mid, RACS-high sources (S > 3 mJy) | CASDA TAP sync, 5-degree Dec strips (async results 404) | 1.6 M rows | ~3.5 min |
| RACS-low2 components (S > 2 mJy) | CASDA TAP sync, strips | 3.45 M rows | ~2.5 min |

All the inputs were in hand within 17 min. The 9.6 GB uncompressed table was parsed straight from the gzip stream in about 90 s.

## F1. The anchor, reproduced from the JVAS catalogue (`r160_jvas.py`)

Set 1 (polarised flux > 1 mJy): **4,400 sources**, the same as Tiwari & Jain. 2D neighbours, S_D maximised with parallel transport,
shuffled-angle null (their Procedure I), 1,000 shuffles:

| n_v | S_D | null | sigma | pairwise C (shuffle null) | median radius |
|---|---|---|---|---|---|
| 6 | 0.3801 | 0.3678 +- 0.0052 | +2.4 | 0.0230 (0.0053 +- 0.0065) | 2.4 deg |
| **10** | **0.3012** | **0.2840 +- 0.0048** | **+3.6** | **0.0180 (0.0034 +- 0.0049)** | **3.3 deg** |
| 15 | 0.2409 | 0.2317 +- 0.0050 | +1.8 | 0.0110 (0.0026 +- 0.0040) | 4.2 deg |

Tiwari & Jain quote a null of 0.2836 +- 0.00495 at n_v = 10, which this reproduces. **Published amplitude, as a pairwise
amplitude: dC_JVAS = 0.0146 +- 0.0049 for pairs within ~3.3 deg.** Their "150 Mpc" is this angle, with every source placed
at z = 1. It is an angular claim: no published 3D test reaches below ~700 Mpc (1503.03482 Fig. 2; 1809.01270: 719 Mpc).

## F2. Samples (`r160_qso.py`, `r160_sample.py`)

- SPICE-RACS goodRM (the catalogue's own cut, polarised S/N >= 8): 333,173 rows. **De-duplicated (rows within 10" are one
  source; the row nearest its tile centre is kept): 188,217.** Duplicates matter here: a duplicated source pair is a
  zero-separation pair with identical angles.
- Spectroscopic quasars: DESI DR1 (ZWARN = 0) 1,589,836, plus DR16Q 413,345 and Milliquas v8 Q/A 74,999 not already
  listed, z >= 0.1. That is 2,078,180 in all.
- **Polarised quasars (3" match): 4,028** (DESI 2,038, DR16Q 660, Milliquas 1,330), median z 1.45.
- Spectral index: the least-squares slope over RACS-low2 (887.5 MHz, components summed by source), RACS-mid (1367.5),
  NVSS (1400) and RACS-high (1655.5), matched within 10". SPICE's own `total_I_flux` is a component flux: against
  RACS-mid it gives a median alpha of +0.25 with a 3-unit spread, so it is not used. Flat means alpha >= -0.5:
  7,377 of 104,735 polarised sources (7 %), and 1,145 of the 1,812 JVAS matches (63 %), as the JVAS selection predicts.
- **Flat-spectrum polarised quasars (the headline sample): 1,547.** Flat and compact (RACS-mid single Gaussian,
  deconvolved major axis < 10"): 237.
- JVAS set-1 sources with a goodRM SPICE-RACS counterpart (5"): **1,812** (the proposal estimated ~1,000).

## F3. Same-source transfer (a calibration; it uses real angles but no alignment statistic; disclosed)

On the 1,812 JVAS sources: **rho = <cos 2(chi0_SPICE - PA_JVAS)> = 0.401 +- 0.015.** Checks on the conventions:
- With the sign of chi0 flipped, -0.004. With chi0 + 90 deg, -0.401. So the conventions agree.
- In the observed frame, before derotation, it is 0.016. The derotation is what recovers the 8.4 GHz angle.

How rho varies across the sample:
- **It rises steeply with fractional polarisation at 0.9 GHz:** 0.15 (p < 2 %, N = 630), 0.48 (2-5 %, N = 880), 0.70
  (p > 5 %, N = 302). Depolarised sources have scrambled 1 GHz angles.
- It hardly depends on |RM|: 0.43 / 0.41 / 0.40 / 0.32 for |RM| < 10, 10-30, 30-100 and > 100 rad m^-2.
- It falls with the chi0 error: 0.48 / 0.41 / 0.31 for errors < 5, 5-10 and 10-20 deg.

## F4. Forecast (`r160_feas.py`; positions, redshifts and uniform random angles only)

The pairwise amplitude C is taken over pairs within R comoving (3D) or within 3.3 deg (2D). Pairs closer than 1' are
dropped. sigma(C) comes from 300 random-angle realisations on the real positions.

Each source's angle correlates with its JVAS angle by rho = rho_J rho_S. Which instrument tracks an underlying aligned
axis better is not measurable with two tracers. So the published amplitude reaches SPICE-RACS angles somewhere between
dC rho^2 = 0.0024 (JVAS the clean tracer), dC = 0.0146 (both equally good: the neutral reading) and dC / rho^2 = 0.091
(SPICE the clean tracer).

| sample | scale | pairs | sigma(C) | S/N, dC rho^2 | **S/N, dC (neutral)** | S/N, dC/rho^2 | 3-sigma amplitude |
|---|---|---|---|---|---|---|---|
| **flat quasars (1,547)** | 150 Mpc | 39 | 0.109 | 0.0 | **0.1** | 0.8 | 0.33 |
| | 300 Mpc | 265 | 0.042 | 0.1 | **0.4** | 2.2 | 0.13 |
| | 500 Mpc | 1,110 | 0.022 | 0.1 | **0.7** | 4.2 | 0.065 |
| | 2D, 3.3 deg | 1,899 | 0.017 | 0.1 | **0.8** | 5.2 | 0.052 |
| all quasars (4,028) | 150 Mpc | 341 | 0.037 | 0.1 | 0.4 | 2.4 | 0.11 |
| | 300 Mpc | 2,406 | 0.015 | 0.2 | 1.0 | 6.2 | 0.044 |
| | 500 Mpc | 10,188 | 0.0071 | 0.3 | 2.1 | 12.8 | 0.021 |
| | 2D, 3.3 deg | 16,316 | 0.0057 | 0.4 | 2.6 | 15.9 | 0.017 |
| flat & compact (237) | 500 Mpc | 25 | 0.14 | 0.0 | 0.1 | 0.6 | 0.42 |
| all flat SPICE sources, 2D, no z (7,377; outside the brief's sample) | 3.3 deg | 26,837 | 0.0044 | 0.5 | 3.3 | 20.6 | 0.013 |

Two readings change the numbers but not the verdict:
- **Pair weighting** by rho_i rho_j (from the p-bin rho above) does not rescue the conservative reading: 0.0-0.2 sigma for
  flat quasars, 0.6 sigma for all quasars in 2D.
- **"3D-only" reading.** In this reading the JVAS 2D excess is made only by pairs within R in 3D, so the 3D amplitude is
  dC rho^2 / f(R). Here f(R) is the fraction of the 3.3-deg pairs that lie within R: 0.019 / 0.089 / 0.18 at 150 / 300 / 500 Mpc.
  At rho^2 this gives flat-subset forecasts of 1.1 / 0.6 / 0.6 sigma.

**Verdict of the gate.** Under the neutral reading the headline subset reaches 0.1-0.8 sigma; under the conservative
reading, at most 0.1. Only the most favourable transfer reading (SPICE angles a perfect tracer and JVAS's 0.4-correlated
ones the noisy copy) gives 2.2-5.2 sigma at 300-500 Mpc and in 2D, and still only 0.8 at 150 Mpc. That reading cannot be
established from these data, and under it a null would not refute the JVAS claim. The brief's stop condition holds:
**stop and report.**

## What stands, and what would reopen it

1. **A measured transfer (new):** the broadband-derotated 0.9 GHz angle recovers the 8.4 GHz JVAS angle to rho = 0.40
   (1,812 sources; 0.70 above 5 % fractional polarisation). Without derotation it recovers nothing (0.016). This is the
   ceiling on any SPICE-RACS re-test of a JVAS-angle claim: the published amplitude is diluted by up to rho^2 = 0.16.
2. **The JVAS anchor is an angular claim at 3 sigma** (0.0146 +- 0.0049 pairwise, 3.3 deg), reproduced to the
   published precision. Its "150 Mpc" assumes z = 1 for every source. So a 3D test at 150 Mpc tests a different
   statement, and the sample sizes now available leave it at 39 flat-quasar pairs.
3. **To reopen:**
   - (a) Calibrate the transfer with a third tracer of the same sources (a triad: e.g. WSRT/VLA 5 GHz angles or a
     second-epoch 8 GHz survey gives rho_J and rho_S separately). That settles which reading holds. Under the neutral
     reading, the redshift-free 2D test on all 7,377 flat SPICE-RACS sources reaches 3.3 sigma.
   - (b) Wait for POSSUM (a rough scaling, not computed). It has several times SPICE-RACS's RM density. If the
     flat-quasar count scaled by k, pairs would scale by k^2 and the neutral-reading S/N by about k. For k ~ 5 that is
     ~3.5 sigma at 500 Mpc and ~4 sigma in 2D, with 150 Mpc still below 1 sigma. Fainter polarised sources match
     spectroscopic quasars less often, so k is an upper bound.
