# cosmo-deep-24: first measurement of Doppler magnification, the galaxy size-number-count dipole, at z = 0.15-0.45 (R-156)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate, including the method check
1b), `finding.md` (passed sections: samples and statistic fixed before the fit, the claim its interval), and
`open-problems/round4/r2-doppler-mag.md` with its verdict in `open-problems/round4/verdicts.md`. Follow the agents' working rules: worktree on
staging; keep-alive (a background loop running the session keep-alive every 20 min); land after every
milestone. Rows in `results-R156.md`. Long jobs with `nohup python -u` and an `until` watcher; foreground waits cap at 10 min, bare
`sleep` is blocked, match processes with `ps -eo pid,args`. `data/` starts empty: report fetch times.

## Trajectory check: already run (2026-11-01, context-free Opus, 2 web searches and 3 fetches), verbatim:
"VERDICT: partly novel: the velocity-induced size/distance–density correlation (the same physics as Doppler convergence, via
Fundamental-Plane/Tully–Fisher sizes) is already measured in DESI DR1 BGS at z ~ 0.07 and is consistent with GR (DESI DR1 PV Survey,
arXiv:2512.03230, with 2512.03229/2512.03231); no one has measured the Doppler-magnification size–number-count dipole with raw galaxy
sizes, nor anything at z ≈ 0.15–0.45 (only forecasts such as arXiv:2511.17858, 1610.05946)". It qualifies as posed (forecast) and
never measured; the claim is restricted to z = 0.15-0.45 (and the raw-size estimator), where the peculiar-velocity surveys cannot
reach. Read 1610.05946 (the estimator), 2511.17858 (the DESI forecast) and 2512.03230 (what the PV survey already measured) in full
first. Do not rerun the trajectory check.

## Step 0: feasibility gate (first hour; stop and report if it fails)
Count DESI DR1 BGS galaxies at 0.15 < z < 0.45 with Legacy Surveys DR9/DR10 Tractor sizes (shape_r and its error) in the overlap;
measure the size scatter at fixed magnitude and colour (the noise of the convergence estimator) from the data; scale the 2511.17858
forecast to that noise and number. Name what could erase or fake the dipole: PSF-size correlations with density, blending in dense
regions, size-magnitude selection, the lensing convergence (even, so it should cancel in the dipole, but check the wide-angle leak),
and estimate each. If the forecast against zero is below 3 sigma per the combined range, stop and report.

## Step 1: method check (blind-review.md 1b), before measuring
The fixed plan (sample and size standardisation, the estimator (the odd dipole of the size-density cross-correlation, Bonvin et al.),
the wide-angle and evolution terms, the PSF and blending systematics with their nulls, the A_D fit and its theory template, the
covariance (jackknife and mocks if available), the redshift bins) alone in a file; the method-check prompt verbatim in a fresh
context-free Opus agent. Close every point it names before measuring, or stop. Verdict verbatim in `results-R156.md`.

## Method
1. Standardised size residual per galaxy (log size against its magnitude-colour-redshift relation, fitted on the same sample).
2. The dipole of the residual-density cross-correlation in bins of separation and redshift; fit A_D against the GR template.
3. Nulls: the even (monopole, quadrupole) parts where the dipole template predicts none; PSF size and seeing maps cross-correlated
   with density; the dipole with randomised redshifts; the estimator on stars.
4. The claim is A_D per bin and combined with its interval, against zero and against GR. No sentence about modified gravity beyond
   what A_D says.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above). Verdict verbatim in `results-R156.md` and under the section. One revision pass on a
   presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-24` in `LOG.md`. Final report: headline, feasibility numbers, method-check verdict, the
   table, nulls, friction. Budget 6 h.
