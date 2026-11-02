# cosmo-deep-6: capture-recapture across gravitational-wave search pipelines: an injection-free count of detectable mergers (R-126)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md`
(what passed and what missed; the misses were claims stronger than their intervals). Follow the agents' working rules: worktree on staging;
keep-alive (a background loop running the session keep-alive every 20 min); land after every
milestone with `git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R126.md` (own file). Long
jobs with `nohup python -u` into a log file and an `until` watcher; logs to files.

## Trajectory check: already run (2026-11-01, context-free Opus, in the dispatching session), verbatim:
"VERDICT: novel: no prior work answers this question". Closest prior work it found: the unified p_astro paper (PRD 108, 083043,
2023; injection-based), Malz et al. arXiv:2607.07272 (ML ranking of subthreshold candidates, GWTC-5.0), arXiv:2401.08709 and
arXiv:2312.06631 (catalogue overlaps described qualitatively), arXiv:2303.17628 (union list), GWTC-5.0 methods arXiv:2605.27224
and the selection-function emulator arXiv:2408.16828 (injection-based). Do not rerun the trajectory check; cite these.
The referee's method caveat, which the analysis must answer by construction: the pipelines share the same strain and similar
templates, so their detections are strongly positively dependent through SNR; closed-population capture-recapture is biased low
under such dependence and log-linear interaction terms cannot identify dependence that runs through one shared variable. So: bin
by network SNR (and chirp mass), model catchability as a function of SNR within bins, use estimators built for dependent lists
(Chao's lower bound; sample-coverage; latent-class / Rasch models), and present the result as a lower bound where that is what the
estimator gives.

## The question
Treat the independent search catalogues on the same public strain (the LVK pipelines GstLAL, PyCBC, MBTA, cWB as reported per
event in GWTC-1 to GWTC-4.0 on GWOSC, plus the external 4-OGC (Nitz et al.) and IAS (Venumadhav, Olsen, Mehta et al.) catalogues)
as overlapping imperfect registries of the same population, and estimate the total number of real mergers above a stated SNR
or FAR threshold in O1-O4a from the overlap pattern, per observing run and per chirp-mass bin, without injections. Compare with
the number the LVK injection-based sensitive volume (VT) implies for the same thresholds (from the public injection sets and the
population papers). Positive: the two disagree beyond their errors in some mass range, which is a statement about the selection
function behind every rate, mass-spectrum and dark-siren result. Null: an injection-free consistency check of VT, the first.

## Data (public; verify reachability from the cloud first; zenodo.org and GitHub release assets 403; GWOSC's own site and the
catalogue papers' tables are to be tried; the IAS and 4-OGC catalogues are in their papers' tables and GitHub source trees)
GWOSC event lists with per-pipeline FAR/SNR/p_astro (GWTC-2.1, GWTC-3, GWTC-4.0 include per-pipeline columns), 4-OGC, IAS O1-O3
catalogues; the LVK O3 and O4a injection sets if reachable (GWOSC hosts some; else the sensitivity curves in the population papers).

## Method
1. Build the lists: one row per candidate, with flags for which of the 4-6 lists found it above a common threshold (define the
   threshold per list carefully: FAR < 1/yr, or p_astro > 0.5; state both). Mind that LVK pipelines' per-event columns exist only
   for the catalogued events: a candidate below every pipeline's threshold is unobserved by construction, which is what the
   estimator infers. State clearly which "lists" are independent runs and which share templates.
2. Estimators: Lincoln-Petersen per pair, Chao lower bound, log-linear models with pairwise interactions (fit, AIC), sample-coverage,
   and a latent-class model with catchability depending on SNR and chirp mass. Bootstrap errors. Report per run (O1+O2, O3a, O3b,
   O4a) and per chirp-mass bin.
3. Compare with the VT-implied counts at the same thresholds (the population papers' sensitive volume and rate posteriors; or the
   injection sets run through the thresholds). The claim is the table: estimated totals with intervals (lower bounds where that is
   what the method gives), the VT-implied totals, their ratio per bin.
4. Nulls and checks: simulate the dependence (a mock population with a shared SNR and per-pipeline noise) to show what the
   estimators return when the truth is known; a null where one list is a random subsample of another.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus
   agent, section alone; hints: the prior works above and anyone who estimated the GW selection function without injections).
   Verdict verbatim in `results-R126.md` and under the section. One revision pass on a presentation point; park on a physics
   point. Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-6` in `LOG.md`. Final report: headline, the table, the mock check, friction.
   Budget 6 h; stop at the budget with everything landed.
