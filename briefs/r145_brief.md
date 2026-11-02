# cosmo-deep-16: do early- and late-type galaxies share one weak-lensing radial acceleration relation? Settled with spectroscopic isolation on independent shear (R-145)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md` (the
passed sections: samples and statistic fixed before the fit, the claim its interval), `open-problems/round2/p2-wl-rar.md` (the
proposal as the referee saw it) and `open-problems.md`. Follow the agents' working rules: worktree on staging; keep-alive (a background loop running
the session keep-alive every 20 min); land after every milestone with `git pull --rebase origin staging
&& git push origin HEAD:staging`. Rows in `results-R145.md`. Long jobs with `nohup python -u` and an `until` watcher; foreground
waits cap at 10 min, bare `sleep` is blocked, match processes with `ps -eo pid,args`; OMP_NUM_THREADS=1; disk ~28 GB, memory ~15 GB.
`research/cosmo-discovery/data/` starts empty: report fetch times.

## Trajectory check: already run (2026-11-01, context-free Opus, 5 web searches and 5 fetches), verbatim:
"VERDICT: partly novel: whether ETG and LTG lensing RARs differ is already measured on KiDS-1000 with contradictory answers
(arXiv:2106.11677 >=6 sigma difference; arXiv:2310.15248 / 2406.09685 same relation under stricter photo-z isolation); what remains
open is an independent-shear (DES Y3/HSC) measurement with spectroscopic (DESI BGS) isolation and mock-calibrated leakage that
decides between them". Read both papers' methods in full before fixing anything: the dispute turns on lens isolation and on the
baryonic-mass estimates, so those are the two steps this job must do differently and state. Also cite Milgrom 2013 (arXiv:1305.3516,
red/blue lens split, earlier). Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour; stop and report if it fails)
Reach DESI DR1 BGS (data.desi.lbl.gov) and the DES Y3 shape catalogue (DES data release site; or HSC PDR/Y3 shear where public);
count isolated DESI BGS lenses (the isolation criterion fixed below) inside the shear footprint, split by type; forecast the
significance of a difference of the size Brouwer et al. report (and of zero) from the shear noise. If the forecast cannot separate
"same relation" from Brouwer's difference at 3 sigma, stop and report: that is the answer, at an hour's cost.

## Method
1. Fix before any shear is measured, in `results-R145.md`: the lens sample (DESI DR1 BGS, stellar masses from a stated SED fit or
   the DESI value-added catalogue), the isolation criterion (no neighbour above a stated fraction of the lens stellar mass within a
   stated projected radius and line-of-sight velocity, using spectroscopic redshifts, with the photometric-redshift fill-in for
   fainter neighbours stated), the type split (Sersic index and colour, both, as Brouwer did), the baryonic masses (stellar plus a
   stated gas prescription, the same for both types; and Mistele's prescription as the alternative), the radial bins, the conversion
   from excess surface density to g_obs (the deprojection both papers used), and the statistic (the ETG/LTG ratio of g_obs at fixed
   g_bar, averaged over a fixed g_bar range).
2. Measure excess surface density around each lens class with the public shear, with the standard corrections (boost, random
   subtraction, source redshift cuts). Covariance from jackknife and from random-point catalogues; nulls: cross shear, random lenses.
3. Calibrate isolation leakage on public mock catalogues with the DESI BGS selection (or a stated substitute), as the referee asks.
4. The claim is the ratio with its interval, under both baryonic prescriptions, and which of the two published answers it favours.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the papers above and anyone who measured type-split lensing RARs or halo masses by colour). Verdict
   verbatim in `results-R145.md` and under the section. One revision pass on a presentation point; park on a physics point.
   Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-16` in `LOG.md`. Final report: headline, feasibility numbers, the ratio table, nulls,
   friction. Budget 6 h.
