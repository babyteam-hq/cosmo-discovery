# cosmo-deep-3: which dataset pays for the DESI DR2 dynamical-dark-energy preference (Shapley attribution and a specification curve) (R-123)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate; run the
trajectory check before any compute, the finding check after writing), `finding.md` (what passed and what missed, verdicts verbatim:
the misses were claims stronger than their intervals). Follow the agents' working rules: worktree on staging; keep-alive (a background loop running
the session keep-alive every 20 min); land after every milestone with
`git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R123.md` (own file). Long jobs with
`nohup python -u` into a log file and an `until` watcher; logs to files.

## The question
DESI DR2 BAO + CMB + supernovae prefer w0wa over LCDM at a significance that depends on the supernova sample (about 2.8 sigma with
Pantheon+, 3.8 with Union3, 4.2 with DES-SN5YR, as DESI reports). Collaborations report leave-one-out tests, which are not
additive and depend on order. Nobody has published the unique fair attribution: the Shapley value of each dataset's contribution to
delta-chi2(w0wa - LCDM) over all 2^n coalitions, nor a preregistered specification curve of the significance across every defensible
analysis choice with its variance decomposed by choice. Both are interval-valued results whichever way they come out: either one
dataset carries the preference (a sharp statement about where the systematics must be) or it is spread (an argument for physics).

## Data (public; check each is reachable from the cloud before planning around it; GitHub release assets and zenodo.org 403)
DESI DR2 BAO likelihood (data.desi.lbl.gov; the Cobaya/desilike package), Planck PR4 compressed and full likelihoods (CamSpec or
HiLLiPoP/LoLLiPoP via cobaya), ACT DR6 lensing (LAMBDA may 403; NERSC's ACT portal https://portal.nersc.gov/project/act/ is
reachable), Pantheon+, Union3 and DES-SN5YR distance likelihoods (their GitHub repositories' source trees are reachable even when
release assets are not), DES Y3 3x2pt if cheap. The compressed CMB (theta*, omega_b, omega_bc, or the Planck distance priors) makes
the 64-coalition grid cheap; the full CMB likelihood is the check on the headline coalition.

## Method
1. Trajectory check first (blind-review.md's before-computing prompt, verbatim, one fresh context-free Opus agent): the question,
   the data, the method, the expected precision; hints: anyone who applied Shapley or game-theoretic attribution to cosmological
   dataset combinations, or published a factorial/specification-curve sweep of the DESI w0wa significance. "Partly novel" with a
   thin open part is a no (entry rule); then do the specification curve alone if that is what is open, or stop and report.
2. Reproduce first: DESI DR2's published delta-chi2 and significance for the three SN samples with your pipeline. If those do not
   come back within 0.3 sigma, stop and report why.
3. Shapley: for n = 6 datasets (DESI DR2 BAO; CMB; Pantheon+; Union3; DES-SN5YR; DES Y3 or ACT lensing), the 64 coalition
   minimisations of delta-chi2(w0wa - LCDM) with a profile likelihood (minimiser, several starts), exact Shapley values, and each
   dataset's chi2 sacrifice at the joint best fit. Repeat for sum m_nu < 0 ("effective negative neutrino mass") if time allows.
4. Specification curve: enumerate and freeze (in the file, before computing) the defensible choices: CMB likelihood, lensing in/out,
   SN sample and calibration, BAO release (DR1/DR2), parametrisation (CPL; a thawing or binned form), priors, neutrino-mass
   prior, metric (delta-chi2 profile; parameter-shift; evidence where cheap). Run the grid with emulators (CosmoPower/CONNECT)
   validated against the published points, or compressed likelihoods; rerun with the full code where the emulator is outside its
   training hull. Plot the ranked curve; decompose the variance in significance by choice (Shapley or ANOVA).
5. The claim is the table: Shapley shares with bootstrap errors, the curve's range, the variance attribution. No sentence about
   "new physics" or "systematics" beyond what the numbers state.
6. Write the result as a `finding.md` section (eighth), then the finding check (blind-review.md's after-writing prompt, verbatim,
   fresh context-free Opus agent, the section alone in a file; hints as in step 1). Verdict verbatim in `results-R123.md` and under
   the section. One revision pass if the objection is a presentation point; park on a physics point. Nothing goes on the site; if
   it passes, say so as the headline and the main session writes the paper.
7. End: `### cosmo-deep-3` in `LOG.md`. Final report: headline (passed or not, verdict verbatim), the
   Shapley table, the curve's range and top attribution, the reproduction check, friction. Budget 6 h; stop at the budget with
   everything landed.
