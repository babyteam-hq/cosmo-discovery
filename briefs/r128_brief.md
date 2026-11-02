# cosmo-deep-8: lensing magnification as a natural instrument for the JWST z > 8 bright-end and mass excess (R-128)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md`
(what passed and what missed; the one pass of this kind, the eRASS1 x ACT section, was written as its interval on samples fixed
before the fit, and that is what passed). Follow the agents' working rules: worktree on staging; keep-alive (a background loop running
the session keep-alive every 20 min); land after every milestone with
`git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R128.md` (own file). Long jobs with
`nohup python -u` into a log file and an `until` watcher; logs to files. Python 3.13 container: pin packages that need it.

## Trajectory check: already run (2026-11-01, context-free Opus, in the dispatching session), verbatim:
"VERDICT: partly novel: magnification bias as a falsification tool for photo-z contamination (arXiv:2507.09142, HST z~4-5) and
lens-model-spread propagation in lensed vs blank UV LFs (arXiv:2205.11526, 1610.00283; UNCOVER z>9 arXiv:2312.05030) are done; a
regression of de-magnified M_UV, stellar mass and size on log mu for JWST z>8 cluster samples with lens-model spread, testing
whether the z>9 bright/massive excess depends on mu, remains open". Other prior work it named: Bouwens et al. 2022 (ApJ 927, 81,
sizes vs lens models), Atek et al. 2018, Chemerynska et al. UNCOVER LF (arXiv:2312.05030: bright sources at low mu ~ 1.3, two lens
models, no mu-dependence test), UNCOVER z > 9 candidates (arXiv:2305.01793), GLIMPSE (arXiv:2509.24881, 2511.07542), a z = 6-8
source-plane reassessment behind Abell 2744 (arXiv:2512.00730; read its full text first: it may propagate a lens-model ensemble),
a nine-field z ~ 7-10 LF paper (arXiv:2610.32307; read it: check for a lensed-vs-blank split), Mason et al. 2015. Do not rerun the
trajectory check; read the two flagged papers in full before computing and record what they do and do not do.

## The question
Foreground magnification mu is set by the cluster, not by the background source, so after de-magnification the intrinsic
properties of z > 8 galaxies (M_UV, stellar mass, size at fixed luminosity) must not depend on mu beyond what the selection
through magnification predicts. Any residual dependence is lens-model error, selection error or photo-z contamination feeding
the claimed excess of bright and massive z > 9 galaxies. The result is an interval either way: a slope consistent with the
selection model (the lensed-field excess is robust to the lensing step, which strengthens the case that it is physical) or a
slope that is not (the excess partly comes from the lensing step).

## Data (public; verify reachability from the cloud first; zenodo.org and GitHub release assets 403; MAST, the DJA (DAWN JWST
Archive) and the surveys' own catalogue hosts are to be tried)
UNCOVER / MegaScience Abell 2744 catalogues and the public lens models (several families per field: the UNCOVER release, Furtak
et al., and others), GLASS-JWST, CANUCS, RELICS-JWST, with their released magnification maps; the blank fields (CEERS, JADES,
PRIMER, COSMOS-Web) catalogues as the mu ~ 1 arm; the high-z candidate lists (photometric z > 8, and spectroscopic where the DJA
has them).

## Method
1. Fix the sample and the statistic before any fit, in `results-R128.md`: the z > 8 photometric selection (one consistent cut
   across fields, the catalogues' own photo-z posteriors), the magnification range, the lens-model families used and how their
   spread enters (per-source mu error from the model ensemble), the regressions (de-magnified M_UV on log mu; stellar mass on
   log mu at fixed M_UV; size residual on log mu), the forward model of the selection (number counts through mu with the survey
   depth and the LF, which predicts the expected slope for a correct lens model and no contamination).
2. Compute the regressions per field and combined, with the mu error propagated; the expected slope from the forward model;
   bootstrap and jackknife-by-field errors; the test of whether the bright (M_UV < -20) and massive (log M* > 9.5) z > 9 sources
   sit at higher or lower mu than the parent sample at the same observed magnitude.
3. Lensed-vs-blank: the de-magnified bright-end counts in the cluster fields against CEERS/JADES at matched depth, with the
   lens-model spread and cosmic variance as errors, as a formal comparison.
4. The claim is the slopes with their intervals and the expected slopes, and the lensed-vs-blank ratio with its interval. No
   sentence about LCDM or early dark energy beyond what the intervals say.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus
   agent, section alone; hints: the prior works above and anyone who tested mu-dependence of de-magnified high-z properties).
   Verdict verbatim in `results-R128.md` and under the section. One revision pass on a presentation point; park on a physics
   point. Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-8` in `LOG.md`. Final report: headline, the table (field, N, slope, interval,
   expected), the lensed-vs-blank ratio, what the two flagged papers turned out to contain, friction. Budget 6 h.
