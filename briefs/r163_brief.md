# cosmo-deep-25b: the spin-asymmetry claims tested on Euclid Q1, restated on the pre-registered raw result (one revision)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Your Request number is the R-N in your prompt; use it for the keep-alive and commit messages.
Read first: `LOG.md` (handoff, `### cosmo-deep-25`), `research/cosmo-discovery/blind-review.md`, `results-R158.md`
(all rows, the method-check verdict and the finding-check verdict at line ~203), the R-158 section of `finding.md`, and R-158's
scripts. Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min); land after every milestone. Rows continue in `results-R158.md`. No new data or classification is needed:
work from R-158's landed outputs (refetch only if a number must be recomputed, and say so).

## Why this job
R-158's gate, verbatim: "VERDICT: novel but not sound: the claim that the JADES excess is confined to the published selection rests on
a post-hoc full-field test with an uncalibrated HST threshold, no NIRCam PSF or eps calibration, and a comparison against galaxies the
published pipeline rejected, so the gap between the two subsets (~2.7σ, about the original result's own significance) cannot locate the
excess, while the Euclid corrected values hinge on a model-based correction for a strongly chiral VIS PSF and a misclassification rate
calibrated only on simulations." Novelty held. Both objections are to claims the result does not need: the post-hoc JADES localisation,
and corrected amplitudes. The pre-registered raw Euclid measurement excludes the published amplitudes (JADES +20 %, +9.3 % at matched
z; the 24.8 % / 11.4 % Galactic-pole dipole; a 10.5 % HST-axis dipole) by margins far larger than any plausible PSF or
misclassification correction. This is the one revision, on presentation.

## Method
1. The headline becomes the raw, model-free Euclid Q1 result (A = +0.44 +- 0.19 % on 291,652 discs; the dipoles on the pre-registered
   axes, raw) and the exclusions of the published amplitudes, stated with the worst-case correction bound: the largest bias that the
   chiral VIS PSF and the misclassification rate could produce, bounded from the data (focal-plane and size/PSF splits, the mirrored
   null) rather than from the simulation model, and shown to be far below the excluded amplitudes. Say plainly that a ~1 % raw excess
   cannot be interpreted without the correction, and claim nothing about it.
2. Remove the JADES localisation claim (the post-hoc full-field test) from the findings; keep only the pre-registered statement that the
   classifier reproduces the published labels (181/182) and amplitude on the published selection, as a check of the classifier.
3. COSMOS stays a "not constraining" row.
4. Revise the section, then the finding check again (after-writing prompt verbatim, fresh context-free Opus agent, section alone; hints
   as R-158's). Last pass: park on any further physics point. Nothing on the site; a pass is the headline of your report.
5. End: append to `### cosmo-deep-25`. Final report: headline, the bound, the table, friction. Budget 3 h.
