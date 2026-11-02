# cosmo-deep-18: the Lyman-alpha forest x Compton-y cross-correlation at z = 2-3.5 (R-148)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md` (passed
sections: samples and statistic fixed before the fit, the claim its interval, errors from simulations), `open-problems/round2/p6-lya-y.md`
(the proposal as the referee saw it). This job uses Compton-y maps: read the y-map lessons first, `results-R131.md`, `results-R132.md`
and the `### cosmo-deep-11` and `### cosmo-deep-12` notes (CIB and radio residuals sank two stacks; the CIB-deprojected ACT map and the
frequency dependence are this job's check, decided before any number is seen). Reuse the ACT DR6 y-map and Websky access from R-136
and R-140. Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min); land after every milestone. Rows in `results-R148.md`. Long jobs with `nohup python -u` and an `until`
watcher; foreground waits cap at 10 min, bare `sleep` is blocked, match processes with `ps -eo pid,args`. `data/` starts empty:
report fetch times.

## Trajectory check: already run (2026-11-01, context-free Opus, 5 web searches and 8 fetches), verbatim:
"VERDICT: novel: no prior work answers this question". What it found: the only forest x SZ measurement is Croft, Banday & Hernquist
2005 (astro-ph/0512380: 3,488 SDSS DR3 forests x WMAP-1, a weak upper limit); quasar host tSZ stacks (1502.01723, Crichton et al.
2016, 1907.11731) measure host halos, not forest-weighted pressure, but bear on the main systematic; <bP_e> tomography with galaxies
stops near z ~ 1 (2502.08850, 2108.01600). Read 0512380 in full. Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour; stop and report if it fails)
Reach DESI DR1 Lya forest data (delta_F, or the spectra and continuum fits) and the ACT DR6 y maps (standard and CIB-deprojected)
and Planck PR4 NILC y; count forest pixels in the overlap per redshift bin; forecast the cross-correlation significance for the
hydrodynamic-simulation prediction of <bP_e> at z = 2-3.5 (from a stated published model) using the y-map noise and the Websky
y + CIB skies. Also forecast the CIB leakage into the standard map. If the forecast for the prediction is below 2.5 sigma on the
CIB-deprojected map, the result is an upper limit: say so in step 0 and continue only if the forecast limit is below the predictions
by a factor that would be informative (state the factor); otherwise stop and report.

## Method
1. Fix before any cross-correlation, in `results-R148.md`: the forest sample (DESI DR1 quasars, the rest-frame Lya window, masks
   for DLAs and BALs as DESI flags them), the redshift bins (2.0-2.5, 2.5-3.0, 3.0-3.6), the y maps (CIB-deprojected ACT as headline,
   standard ACT and Planck as checks), the estimator (pixel-weighted cross-correlation in angular bins 1'-60', or the angular
   cross-power, with the weights stated), the masks (point sources, clusters above a stated y, the y-map masks), and the model
   (b_F x <bP_e>(z) times the projection, with b_F from DESI's own Lya auto-correlation; the quasar-host contribution at the forest's
   own sightline modelled or excluded by an angular cut, stated).
2. Errors from rotated and flipped y maps and forests shuffled between sightlines; nulls: ACT minus Planck difference maps, the
   frequency dependence (the standard vs CIB-deprojected difference predicted from Websky), the forest red-side (no Lya) control.
3. The claim is the table: per redshift bin, the cross-correlation amplitude and <bP_e> with its interval (or upper limit), on each
   map, and against the published predictions. No sentence about feedback beyond what the numbers say.
4. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above and anyone who measured high-redshift pressure or forest cross-correlations with CMB
   secondaries). Verdict verbatim in `results-R148.md` and under the section. One revision pass on a presentation point; park on a
   physics point. Nothing on the site.
5. End: `### cosmo-deep-18` in `LOG.md`. Final report: headline, feasibility numbers, the table, nulls,
   friction. Budget 6 h.
