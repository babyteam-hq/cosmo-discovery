# cosmo-deep-14: the first tomography of the near-infrared background from real SPHEREx data, per spectral channel against DESI spectroscopic slices (R-138)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md`
(the passes were written as their intervals on samples and statistics fixed before the fit), `screen-2/q2-spherex-tomography.md`
(the proposal as the referee saw it) and `screen-2/verdicts.md`. This is a new data neighbourhood for the project: no code to reuse
beyond the general patterns. Follow the agents' working rules: worktree on staging; keep-alive (a background loop running
the session keep-alive every 20 min); land after every milestone with
`git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R138.md` (own file). Long jobs with
`nohup python -u` into a log file and an `until` watcher; OMP_NUM_THREADS=1 in workers; disk ~28 GB, memory ~15 GB: stream images,
never hold the survey in memory.

## Trajectory check: already run (2026-11-01, context-free Opus, in the dispatching session), verbatim:
"VERDICT: partly novel: the SPHEREx×galaxy-slice b dI/dz tomography observable and method are forecast in arXiv:2109.10914 (method
demonstrated on GALEX×SDSS in arXiv:1810.00885), and a measured broadband NIR b dI/dz at 1.1/1.8 µm for z<1 exists from CIBER×DESI
Legacy/HSC (arXiv:2608.12116); what remains open is the first measurement with real SPHEREx data, per spectral channel against
DESI spectroscopic slices, separating H-alpha/[OIII]/Pa-alpha line emission from the continuum". Read in full before computing:
arXiv:2109.10914 (the estimator and forecast; use its notation), arXiv:1810.00885 (the clustering-redshift intensity estimator and
its nulls), arXiv:2608.12116 (the CIBER measurement: compare at 1.1 and 1.8 µm where the channels overlap). The referee warns the
SPHEREx team is likely preparing this; work fast and land early. Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour; stop and report if it fails)
1. Reach the SPHEREx Quick Release spectral images and the zodiacal model at IRSA from the cloud (IRSA's SPHEREx overview_qr page;
   IRSA worked for earlier jobs); find the deepest coverage (the ecliptic-pole deep fields) and how many images overlap DESI DR1.
2. Reach DESI DR1 spectroscopic catalogues (BGS, LRG, ELG, QSO) and their random catalogues (data.desi.lbl.gov, public DR1).
3. Estimate the signal-to-noise per channel and slice from the overlap area and the forecast's scaling. If the overlap with DESI
   DR1 is too small for any channel to reach 3 sigma, or the images cannot be fetched, stop: write the row and report. That is a
   useful answer and costs an hour.

## The question
The bias-weighted near-infrared background intensity per unit redshift, b dI/dz, in each SPHEREx channel, from cross-correlating
masked, zodi-subtracted SPHEREx images with DESI spectroscopic galaxy overdensity in redshift slices. Lines (H-alpha, [OIII],
Pa-alpha) should appear where their observed wavelength crosses a channel; the claim is b dI/dz per channel and slice with
intervals, the line amplitude against the continuum where detected, and the comparison with CIBER at 1.1 and 1.8 µm.

## Method
1. Fix before any cross-correlation, in `results-R138.md`: the fields, channels (or channel groups if per-channel S/N is too low),
   redshift slices (dz = 0.05-0.1 over the DESI range), the masks (Gaia stars and Legacy Survey sources above a fixed flux, with a
   fixed radius law), the angular scales (a fixed annulus in comoving Mpc per slice), the estimator (that of arXiv:1810.00885, with
   the galaxy bias from each slice's own clustering), and the line model (lines at their observed wavelengths, continuum smooth in z).
2. Nulls: redshift-shuffled galaxies; cross-correlation with Gaia stars; the zodiacal model alone; half-survey difference maps (and
   detector or time halves if the release provides them). Errors: jackknife over fields or sky patches, checked against the
   difference maps.
3. The claim is the table: per channel and slice, b dI/dz with its interval; the line detections with significance; the CIBER
   comparison. A null is reported as limits per channel. No sentence about galaxy evolution beyond what the table says.
4. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus
   agent, section alone; hints: the three works above and anyone who measured SPHEREx or near-infrared background
   cross-correlations). Verdict verbatim in `results-R138.md` and under the section. One revision pass on a presentation point;
   park on a physics point. Nothing on the site; a pass is the headline of your report.
5. End: `### cosmo-deep-14` in `LOG.md`. Final report: headline, the feasibility numbers, the table, the null
   checks, friction. Budget 6 h; stop at the budget with everything landed.
