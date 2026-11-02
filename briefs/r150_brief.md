# cosmo-deep-19: the circumgalactic dust reddening profile around galaxies at z ~ 0.8-1.6 (R-150)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md` (passed
sections: samples and statistic fixed before the fit, the claim its interval), `open-problems/round2/p7-dust-z1.md` (the proposal as
the referee saw it). Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min); land after every milestone. Rows in `results-R150.md`. Long jobs with `nohup python -u` and an `until`
watcher; foreground waits cap at 10 min, bare `sleep` is blocked, match processes with `ps -eo pid,args`. `data/` starts empty:
report fetch times.

## Trajectory check: already run (2026-11-01, context-free Opus, 4 web searches and 6 fetches), verbatim:
"VERDICT: partly novel: absorber-selected reddening (Mg II, arXiv:0706.0898, and the dust budget in Mg II halo gas, arXiv:1204.1978)
and individual quasar sightlines inside 20 kpc (GOTOQs, arXiv:2601.02480) already constrain dust around z ~ 0.5–1 galaxies, and
galaxy-selected stacked profiles exist only at z ≲ 0.5 (0902.4240, 2503.04098, 2606.24159); a stacked reddening profile against
impact parameter around galaxies at z ~ 1–1.5, and how it has changed since z ~ 0.3, remains open". It qualifies as posed (the LSST
Science Book, arXiv:0912.0201, names halo-dust evolution as a goal) and never measured. Read Ménard et al. 2010 (0902.4240) in full
first: its estimator, its nulls and its magnification treatment are the template. The referee's warning: DESI quasars are
colour-selected, so reddened quasars drop out of the background sample. Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour; stop and report if it fails)
Count DESI DR1 foreground ELGs (0.8 < z < 1.6) and LRGs, background quasars at z > 2 (and z > foreground + 0.5), and pairs within
20 kpc-3 Mpc; measure the per-quasar reddening-estimator noise from the data; forecast the profile significance for Ménard's
z ~ 0.3 amplitude evolved to z ~ 1.2 under (a) no evolution and (b) the dust-to-gas scaling of a stated model. Quantify the
selection loss: the fraction of quasars that DESI's targeting would lose at the reddening expected, from the published targeting
selection applied to reddened composites. If the forecast is below 3 sigma, or the selection bias exceeds the signal and cannot be
modelled, stop and report.

## Method
1. Fix before any stack, in `results-R150.md`: foreground samples and redshift bins, background quasar cuts, the reddening estimators
   (continuum slope against the DESI composite at each quasar's redshift; g-r, r-z, z-W1 colour residuals), the impact-parameter
   bins, the Mg II flagging at the foreground redshift (reported with and without), the magnification term (fitted from the same
   pairs, as Ménard did), the selection-bias correction (from step 0), and the fit (A_V(r_perp) as a power law, compared with
   Ménard's z ~ 0.3 profile as a ratio with its interval).
2. Nulls: foregrounds behind the quasar, random foreground positions, stars as backgrounds, Galactic-dust (SFD/Planck) splits;
   errors by jackknife.
3. The claim is the table: per redshift bin, A_V at fixed impact parameters and the power-law amplitude and slope with intervals,
   and the ratio to z ~ 0.3. No sentence about outflows beyond what the ratio says.
4. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above and anyone who measured halo dust reddening at z > 0.5). Verdict verbatim in
   `results-R150.md` and under the section. One revision pass on a presentation point; park on a physics point. Nothing on the site.
5. End: `### cosmo-deep-19` in `LOG.md`. Final report: headline, feasibility numbers, the table, nulls,
   friction. Budget 6 h.
