# cosmo-deep-17: CMB lensing of z ~ 5 Lyman-break galaxies, and the amplitude of structure at z ~ 5 (R-146)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md` (the passed
sections: samples and statistic fixed before the fit, the claim its interval, covariance from released simulations),
`open-problems/round2/p1-z5-lensing.md` (the proposal as the referee saw it), and `results-R139.md` (the quasar x CMB-lensing job:
reuse its pseudo-C_ell, magnification and ACT/Planck handling; and its lesson: inconsistent sky legs averaged without inflating the
error failed the gate). Reuse `erass1_act.py` for the ACT DR6 release and its 400 simulations. Follow the agents' working rules: worktree on staging;
keep-alive (a background loop running the session keep-alive every 20 min); land after every milestone.
Rows in `results-R146.md`. Long jobs with `nohup python -u` and an `until` watcher; foreground waits cap at 10 min, bare `sleep` is
blocked, match processes with `ps -eo pid,args`. `data/` starts empty: report fetch times (R-139: 9 min catalogues and maps, 85 min
for the 400 ACT simulations).

## Trajectory check: already run (2026-11-01, context-free Opus, 5 web searches and 8 fetches), verbatim:
"VERDICT: novel: no prior work answers this question". What it found: the only z >= 4 measurement is Miyatake et al. 2022
(arXiv:2103.15862: ~1.47M HSC g-dropouts at z ~ 4 x Planck PR3 kappa, 5.1 sigma, halo mass ~2.9e11 h^-1 M_sun, a joint
(Omega_m, sigma_8) at z ~ 4 somewhat low against Planck); forecasts only otherwise (arXiv:1904.13378, 2401.13166, 2511.22243).
Read 2103.15862 in full before fixing anything (whether it used r-dropouts at all; its interloper and magnification treatment). The
referee's cautions: the r-dropout count is unverified; at z ~ 5 low-z interlopers and the CIB matter more. Do not rerun the
trajectory check.

## Step 0: feasibility gate (first hour; stop and report if it fails)
Reach the HSC public r-dropout catalogues (HSC-SSP PDR3 and the GOLDRUSH dropout samples, with their selection completeness and
interloper estimates), count the r-dropouts with masks, the overlap with ACT DR6 and Planck PR4/PR3 lensing, and forecast the
cross-power significance at the bias expected for z ~ 5 LBGs (b ~ 4-6). Reproduce the z ~ 4 g-dropout x Planck result first as the
end-to-end check (Miyatake's 5.1 sigma within its error). If the z ~ 5 forecast is below 3 sigma, stop and report.

## Method
1. Fix before any cross-spectrum, in `results-R146.md`: the samples (r-dropouts by the published colour cuts, magnitude limits,
   masks), the multipole range, the estimator (pseudo-C_ell, mode coupling corrected), the model (linear bias with the published
   dropout dN/dz including the interloper fraction, plus magnification bias from the measured number-count slope, CIB leakage
   bounded with the CIB-deprojected / profile-hardened lensing variants), and the fit: b at z ~ 5 from kappa x g, and, with the
   dropouts' auto-spectrum, the amplitude sigma_8(z ~ 5) (or the growth amplitude relative to Planck LCDM) as Miyatake did at z ~ 4.
2. Covariance from the 400 ACT simulations and the Planck simulations; nulls: lensing curl, random catalogues, ACT-Planck agreement
   on the overlap (and if the legs disagree, inflate the error by their scatter, as R-139's gate demanded), the z ~ 4 reproduction.
3. The claim is the table: N, the cross-power significance, b(z ~ 5) and the growth amplitude with intervals, against Planck LCDM,
   with the z ~ 4 reproduction beside it.
4. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus agent,
   section alone; hints: the works above and anyone who measured structure growth at z >= 4). Verdict verbatim in `results-R146.md`
   and under the section. One revision pass on a presentation point; park on a physics point. Nothing on the site.
5. End: `### cosmo-deep-17` in `LOG.md`. Final report: headline, feasibility numbers, the table, nulls,
   friction. Budget 6 h.
