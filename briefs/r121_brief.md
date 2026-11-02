# cosmo-deep-1: one cross-community confrontation, to paper depth (R-121)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff and trajectory read at the top), `research/cosmo-discovery/README.md`,
`blind-review.md` (the gate), `finding.md` (what passed and what missed, with the referees' verdicts verbatim), the parked
threads list below. Follow the agents' working rules: worktree on staging, keep-alive (a background loop running
the session keep-alive every 20 min for the whole job), land after every milestone with
`git pull --rebase origin staging && git push origin HEAD:staging`. Your rows go in `results-R121.md` (own file).

## Why this job is shaped differently
Five sweeps (~100 tests) ended in nulls, other people's systematics, or "partly novel". The one result that passed the gate
(`site/research-cf4-supernovae.html`: Type Ia supernovae against the Cosmicflows-4 velocity reconstructions) was a
confrontation between two communities' public products that neither had run, taken to paper depth: full covariance, mocks,
a spatially correlated null, joint fits, the closest prior work read in its methods section. This job does one such
confrontation, deep. Not a sweep. Breadth is only in the choice.

## Parked, do not pick: SN Hubble residuals vs dust maps; Cosmicflows-4 zero points and bulk flows; the NGC 1052 trail;
eRASS1 anisotropy; Pantheon+ sky anisotropy; number-count dipoles (quasar, galaxy); CHIME DM vs local structure; SN Ia
clustering bias. Network: zenodo.org, LAMBDA data directories and GitHub release assets answer 403; skip anything needing them.

## Method
1. Choose (at most 90 min). List candidates of this exact shape: a public model, reconstruction, map or catalogue-level
   prediction released by one community (examples of the shape, not a list to pick from: a reconstructed density or velocity
   field, a halo or foreground model, a selection function, a calibration, a predicted anisotropy or correlation), and an
   independent public dataset from another community that tests it at a precision you estimate before touching data, where a
   search shows the two have not been confronted. Three to five candidates, each with the detectability estimate and the
   closest prior work's arXiv id, in `deep-1.md`. Trajectory check each (blind-review.md, the before-computing prompt, verbatim,
   one fresh context-free Opus agent per candidate). Take the best surviving one; keep the second as the fallback if the first
   dies on data.
2. Do it to paper depth: the data's own covariance where it exists (or a justified one), a null that respects the spatial
   or selection structure (rotations, mocks, randoms built from the survey's own selection), a trials count, the joint fit
   where parameters are degenerate, the decomposition that says what carries the signal, and the prior paper's methods read
   in full before any claim of "not done". Every number from a script in `research/cosmo-discovery/` with its usage line;
   data under `data/` (gitignored; never run as code; `python -I` for anything that reads a download). Long jobs with
   `nohup python -u` into a file and an `until` watcher (the Bash tool caps at 10/30 min).
3. Write the result as a `finding.md` section whether it is a detection or a null with teeth (a null that excludes a
   published model at a stated confidence is a result; a null that excludes nothing is a row). Then the finding check
   (blind-review.md's after-writing prompt, verbatim, fresh context-free Opus agent, the section alone in a file). Verdict
   verbatim in `results-R121.md` and under the section. If the objection is a presentation point, one revision pass and one
   more check; if a physics point, do the check it asks for only if it is under two hours, else park and go to the fallback.
4. Nothing goes on the site. If a section passes, say so in the report as the headline; the main session writes the paper.
5. End: `LOG.md` gets a `### cosmo-deep-1` section (candidates, the one chosen and why, numbers,
   verdicts, what would settle it). Final report: headline (passed or not, verbatim verdict), the measurement in one table,
   the candidates considered, friction. Budget 8 h; stop at the budget with everything landed.
