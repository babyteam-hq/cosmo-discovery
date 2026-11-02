# cosmo-deep-10: archival strongly lensed supernovae by forced photometry at known lens positions, below the alert threshold (R-130)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md`
(what passed and what missed; the eRASS1 x ACT section passed by being written as its interval on samples and statistics fixed
before the fit; R-127's ZTF scene photometry section records what ZTF difference imaging can and cannot do: resolved sources
brighter than r ~ 20.5 are fine, r ~ 21 is not; reuse its IRSA access code). Follow the agents' working rules: worktree on staging; keep-alive (a
background loop running the session keep-alive every 20 min); land after every milestone with
`git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R130.md` (own file). Long jobs with
`nohup python -u` into a log file and an `until` watcher; logs to files.

## Trajectory check: already run (2026-11-01, context-free Opus, in the dispatching session), verbatim:
"VERDICT: partly novel: archival per-lens lensed-SN searches of the same known-lens catalogues exist (Sheu et al. arXiv:2301.03578
on 5,807 Legacy Survey lenses in DECam epochs; Magee et al. arXiv:2303.15439 cross-matching ZTF transients to known lenses, null;
Townsend et al. arXiv:2405.18589 ZTF archive candidates); what remains open is ZTF/ATLAS forced photometry at every lens position
below the alert threshold, with injection-recovery efficiency and a per-lens rate upper limit compared against the
Goldstein/Wojtak/Arendse forecasts". Other prior work it named: Shu et al. 2018 arXiv:1803.07569 (rates in known lenses), the
LCOLSS pilot arXiv:2111.01680 (114 lenses monitored, nothing), arXiv:2406.00052 (ZTF lensed-SN detectability), the review
arXiv:2406.13519 (about 7 lensed SNe per year expected in ZTF's footprint from MasterLens). Read Sheu et al. in full first: if its
body contains an efficiency-calibrated rate limit, the open part narrows to the ZTF/ATLAS sub-threshold search and the forecast
comparison; say so in the first row. Do not rerun the trajectory check.

## The question
Has any strongly lensed supernova or lensed transient occurred behind the known galaxy-scale strong lenses during 2018-2025 that
the alert-based searches missed (they trigger only above the ZTF alert threshold, about r ~ 20.5), and what is the lensed-SN rate
per lens class, with an injection-calibrated efficiency, compared with the forecasts that set Rubin's expectations? Positive: a
new lensed SN or lensed transient (an object, defensible on its own). Null: a rate upper limit per lens class at the forecast
level, with the efficiency measured rather than assumed.

## Data (public; verify reachability from the cloud first; zenodo.org and GitHub release assets 403; IRSA worked for R-127)
Lens catalogues: the DESI Legacy Survey strong-lens candidates (Huang et al. 2020/2021, Storfer et al. 2022; grades A and B),
SLACS, BELLS, SL2S, the Lemon et al. compilation, MasterLens; ZTF forced photometry: the IRSA forced-photometry service (batch
requests by position; it is slow at thousands of positions, so prioritise by Einstein radius and lens-galaxy brightness, and
state the subset) or R-127's difference-image cutout code; ATLAS forced photometry (the public server, rate-limited; a subset).

## Method
1. Fix the lens list (by catalogue and grade), the time window, the detection statistic and the lens classes before any light
   curve is pulled, in `results-R130.md`. Fix the injection set: lensed Ia and core-collapse light curves (magnified by the
   lens model's typical mu, at the source redshift distribution of the forecasts) injected into the real cadence and noise at the
   same positions, for the efficiency per lens class and per peak magnitude.
2. Pull the forced photometry at each lens position (and at each image position where known); detect transients in the
   difference light curves with the fixed statistic; for each candidate, the standard vetting (host, variability class, the
   lens model's delay range, colour); report every candidate with its vetting, not only survivors.
3. The rate: the number of lensed-SN candidates surviving vetting against the injection-calibrated expected number per lens
   class from the Goldstein et al. 2019, Wojtak et al. 2019 and Arendse et al. 2024 forecasts scaled to the lens list and the
   window; the upper limit (Poisson, 95 %) where none survive.
4. The claim is the table: per lens class, N lenses, the efficiency, the candidates, the forecast expectation, the ratio or limit
   with its interval. No sentence about lens-population physics beyond what the ratio says.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus
   agent, section alone; hints: the prior works above and anyone who did forced-photometry lensed-transient searches). Verdict
   verbatim in `results-R130.md` and under the section. One revision pass on a presentation point; park on a physics point.
   Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-10` in `LOG.md`. Final report: headline, the table, the injection efficiencies,
   every candidate, what Sheu et al. turned out to contain, friction. Budget 8 h; stop at the budget with everything landed.
