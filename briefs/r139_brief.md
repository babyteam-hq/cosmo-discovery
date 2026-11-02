# cosmo-deep-15: the clustering bias of quasars at 3.5 < z < 5 from CMB lensing (R-139)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md`
(the eRASS1 x ACT lensing section is the pass this job copies in shape: samples and statistic fixed before the fit, the claim its
interval, covariance from released simulations), `open-problems/round1/g-qso-z4-lensing.md` (the proposal as the referee saw it)
and `open-problems.md` (how it was chosen). Reuse `erass1_act.py` (ACT DR6 lensing release handling, masks, the 400 simulations,
`simstream.sh`). `research/cosmo-discovery/data/` is gitignored and empty in a new container: refetch (NERSC for ACT; PLA or NERSC
for Planck PR4 lensing and its simulations; SDSS DR16Q from SDSS; DESI DR1 quasars from data.desi.lbl.gov) and report the refetch
time. Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min); land after every milestone with `git pull --rebase origin staging && git push origin HEAD:staging`. Rows in
`results-R139.md` (own file). Long jobs with `nohup python -u` into a log file and an `until` watcher; OMP_NUM_THREADS=1 in
workers; disk ~28 GB, memory ~15 GB.

## Trajectory check: already run (2026-11-01, context-free Opus, 3 web searches and 4 arXiv API queries), verbatim:
"VERDICT: novel: no prior work answers this question". What it found: quasar x CMB-lensing bias measurements stop at z <= 3.5
(DESI DR1 x Planck PR4, arXiv:2506.22416 and 2512.17865, top bin centred near z ~ 2.75; Quaia x ACT/Planck, arXiv:2507.08798 and
2402.05761; eBOSS x Planck, arXiv:1809.04196); at z ~ 4 the bias is known only from the quasar auto-correlation (Shen et al. 2007,
astro-ph/0702214, b ~ 13-14 with a dependence on the treatment of negative points) and quasar-LBG cross-correlation (He et al.
arXiv:1704.08461, b ~ 3-6, and arXiv:1507.05292), which disagree. The referee's cautions: the 3-5 sigma forecast looks optimistic
(Planck noise, partial ACT overlap, ~10^4 quasars), magnification bias at z ~ 4 is a leading systematic through the luminosity-
function slope, and check whether the DESI DR1 x PR4 top bin leaks information above z = 3.5. Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour; stop and report if it fails)
Count the quasars at 3.5 < z < 5 in SDSS DR16Q and DESI DR1 (with the selection masks), the overlap with ACT DR6 and Planck PR4,
and forecast the cross-power significance for b = 3, 6 and 13 (Gaussian covariance with the released lensing noise; the
magnification term included). If even b = 13 forecasts below 2.5 sigma, stop: write the row and report. If b = 13 is detectable
but b = 3 is not, the result can still separate the two published answers: say so and continue.

## The question
The linear bias of quasars at 3.5 < z < 5 from the cross-power of their overdensity with CMB lensing convergence, and with it the
discrimination between the auto-correlation (b ~ 13) and quasar-galaxy cross-correlation (b ~ 3-6) answers. The claim is b with
its interval, per survey and combined, and the probability ratio between the two published values.

## Method
1. Fix before any cross-spectrum, in `results-R139.md`: the samples (DR16Q and DESI DR1 at 3.5 < z < 5, the redshift-quality
   cuts, the magnitude limits), the masks (survey selection, lensing masks, |b| > 20 deg), the multipole range, the estimator
   (pseudo-C_ell, mode-coupling corrected), the model (linear bias with the spectroscopic dN/dz, plus magnification bias with the
   slope s measured from the samples' own luminosity function at the flux limit, its uncertainty propagated), and the fit.
2. Covariance from the 400 ACT simulations and the Planck PR4 simulations; nulls: the lensing curl mode, the quasar randoms, an
   ACT-Planck consistency test on the overlap, and the same pipeline at 2.2 < z < 3.5 reproducing the published DESI x PR4 bias
   as the end-to-end check.
3. The claim is the table: per survey and lensing map, N quasars, the cross-power significance, b with its interval, and the
   combined value; the probability ratio for b = 13 versus b = 4. No sentence about duty cycles beyond what b says.
4. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus
   agent, section alone; hints: the works above and anyone who measured high-redshift quasar bias with lensing). Verdict verbatim
   in `results-R139.md` and under the section. One revision pass on a presentation point; park on a physics point. Nothing on
   the site; a pass is the headline of your report.
5. End: `### cosmo-deep-15` in `LOG.md`. Final report: headline, the feasibility numbers, the table, the null
   checks, the refetch time, friction. Budget 6 h; stop at the budget with everything landed.
