# cosmo-deep-19b: the large-scale dust reddening around z ~ 1 emission-line galaxies, as its upper limit (R-152)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff, `### cosmo-deep-19`), `research/cosmo-discovery/blind-review.md`,
`results-R150.md` (all rows), the R-150 section of `finding.md`, and R-150's scripts. Follow the agents' working rules: worktree on staging; keep-alive
(a background loop running the session keep-alive every 20 min); land after every milestone. Rows continue
in `results-R150.md`. Long jobs with `nohup python -u` and an `until` watcher; foreground waits cap at 10 min, bare `sleep` is blocked,
match processes with `ps -eo pid,args`. R-150 fetched its data (2.9 GB) in 88 s.

## Why this job
R-150's gate, verbatim: "VERDICT: novel but not sound: the LRG detection and the ELG–LRG contrast are not separated from an LRG-tied
systematic (stellar null 1.8 ± 0.7, strongest in z−W1) and a magnification leak of comparable size, and the 0.25–3 Mpc scales cannot
show halo dust." Novelty held, and the referee called the ELG limit the robust part. The fix is the claim, not the data: this is the
one revision.

## Method
1. The headline becomes the ELG result alone: around DESI DR1 ELGs at 0.8 < z < 1.6, the reddening of z > 2 quasars at 0.25-3 Mpc is
   rho = 0.04 +- 0.29 times Ménard et al.'s z ~ 0.36 profile at the same comoving scales (< 0.61 at 95 %), stated as large-scale
   (two-halo-dominated) reddening, not halo dust; say what scales would be needed for halo dust and why the data cannot reach them
   (~10^3 pairs inside 130 kpc).
2. Check that the ELG limit is free of the referee's two problems, with numbers: the stellar null and the magnification leak computed
   for the ELG sample alone (not the LRGs), and the limit's change when each is propagated. If either moves the ELG limit by more than
   its error, say so and state the limit with it included.
3. The LRG rows move to a "not separated from systematics" row with the stellar null beside them, and no ELG-LRG contrast is claimed.
4. Compare the limit, at matched comoving scales, with Ménard's profile and with the absorber-selected dust budget at z ~ 1
   (arXiv:1204.1978), as a ratio with its interval.
5. Revise the section, then the finding check again (after-writing prompt verbatim, fresh context-free Opus agent, section alone;
   hints as R-150's). Last pass: park on any further physics point. Nothing on the site; a pass is the headline of your report.
6. End: append to `### cosmo-deep-19`. Final report: headline, the ELG systematics numbers, the table, friction. Budget 3 h.
