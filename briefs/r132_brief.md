# cosmo-deep-12: the one-halo thermal SZ signal of group catalogues at 10^12.5 to 10^13.5 solar masses on the ACT DR6 Compton-y map (R-132)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md`
(the eRASS1 x ACT lensing section is the one pass of this kind: samples and statistic fixed before the fit, the claim its interval;
the misses since were claims that outran their intervals or methods the referee could name a flaw in). R-131 (running now) is
stacking the same ACT DR6 and Planck y maps on the eRASS1 clusters with `erass1_ysz.py`: reuse its y-map access, UPP templates,
footprint and rotation nulls; do not duplicate its fetches (check `research/cosmo-discovery/data/` on staging's notes). Follow
the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive every
20 min); land after every milestone with `git pull --rebase origin staging && git push origin HEAD:staging`. Rows in
`results-R132.md` (own file). Long jobs with `nohup python -u` into a log file and an `until` watcher; OMP_NUM_THREADS=1 in workers.

## Trajectory check: already run (2026-11-01, context-free Opus, in the dispatching session), verbatim:
"VERDICT: partly novel: Planck-resolution stacks on spectroscopic and DESI DR9 groups already measured the group Y-M down to about
10^12.5 M_sun (Lim et al. arXiv:1710.06856; arXiv:2201.12591), and ACT DR6 stacks exist only for galaxies or LRGs (arXiv:2607.08721,
2609.08939, 2608.00164) with no group-catalogue Y-M; the ACT-resolution one-halo Y-M amplitude and slope on Tully 2015/2MRS/DESI
DR9 groups below 10^13.5 M_sun remains open". Other prior work it named: Greco et al. arXiv:1608.04160 (Planck y x SDSS Yang
groups, six mass bins), arXiv:1706.03753 (the two-halo term dominates low-mass Planck stacks), arXiv:2510.12553 (ACT kSZ on the
DESI DR9 Yang groups), arXiv:2507.04476 (ACT DR5 tSZ clusters). Read in full before computing: arXiv:2607.08721 (Battaglia et al.
2026: ACT DR6 y on Planck's Locally Bright Galaxies; a factor-of-two excess in galaxies hosting radio sources), arXiv:2609.08939
(Moore et al. 2026: ACT DR6+Planck y on DESI LRGs, BGS and eROMaPPer clusters; could not fit the two-halo term; lowest-mass bins go
negative with radio sources in) and arXiv:2608.00164 (Siegel et al. 2026: ACT on DESI LRGs with WL masses at 10^13-10^14, pressure a
factor two below FLAMINGO). Their radio-source handling and mass floors are the two things this job must do better. Do not rerun the
trajectory check.

## The question
The group-scale Y-M relation (10^12.5 to 10^13.5 M_sun) is where feedback breaks self-similarity, and the published measurements are
at Planck's 10' resolution, where the one-halo term is blended with the two-halo term and neighbours. ACT DR6's y map resolves the
one-halo term. Stack it on the group catalogues with halo masses (the DESI DR9 Yang et al. 2021 groups, where their masses are
reliable, above 10^13 by the referee's warning; the Tully 2015 2MRS groups at low z; the SDSS Yang groups) in halo-mass bins,
fit the one-halo Y-M amplitude and slope with the two-halo term fitted from the catalogue's own clustering, and compare with the
Planck-resolution relations and with the self-similar slope. The claim is the relation: amplitude and slope with intervals per
catalogue, and the ratio to the Planck-era values with its interval.

## Data (public; verify reachability; zenodo.org, LAMBDA data directories and GitHub release assets 403; NERSC and VizieR work)
ACT DR6 Compton-y maps, masks and (if released) simulations or noise realisations (NERSC ACT portal; R-131's notes say where);
the Planck 2015 MILCA y map for the comparison; the DESI DR9 group catalogue (Yang et al. 2021, public), the Tully 2015 group
catalogue (VizieR J/AJ/149/171), the SDSS DR7 Yang groups; a radio-source catalogue reachable from the cloud (NVSS/FIRST via VizieR
or the Legacy Survey's radio flags) for the contamination mask.

## Method
1. Fix before any stack, in `results-R132.md`: the catalogues and their mass floors, the mass-redshift bins, the aperture (a UPP
   matched filter at each group's R500 from its halo mass, convolved with the map beam; a fixed aperture as check), the radio-source
   treatment (mask groups with a radio source above a stated flux within the aperture; report with and without), the two-halo
   model (from the catalogue's measured group-group correlation times the UPP, or a halo-model two-halo with free amplitude), the
   statistic (one-halo Y_500 per bin after two-halo subtraction; A and B of Y ∝ A (M/M_pivot)^B E(z)^(2/3) D_A^-2).
2. Stack with the per-group mass uncertainty propagated (the catalogues' quoted scatter); covariance from rotated-map stacks and
   random positions (and y-map sims where available); the stack on the groups' own neighbours as the two-halo check.
3. Fit A and B per catalogue; the same on the Planck map at Planck's resolution as the like-for-like comparison with Lim et al.
   and arXiv:2201.12591; the CIB-deprojected ACT variant as the foreground check; the radio-masked and unmasked results side by side.
4. The claim is the table: per bin Y_500 with its interval; A, B with intervals per catalogue and map; the Planck comparison. No
   sentence about feedback models beyond what the slope and amplitude say.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus
   agent, section alone; hints: the prior works above). Verdict verbatim in `results-R132.md` and under the section. One revision
   pass on a presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-12` in `LOG.md`. Final report: headline, the table, the radio check, friction.
   Budget 6 h; stop at the budget with everything landed.
