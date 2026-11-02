# cosmo-deep-9: a contamination budget for the ACT DR6 and Planck PR4 CMB-lensing maps from templates that cannot carry signal (R-129)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate), `finding.md`
(what passed and what missed; the eRASS1 x ACT section passed by being written as its interval on samples and statistics fixed
before the fit). R-122's `erass1_act.py` already fetches and handles the ACT DR6 release, its masks and its 400 simulations: reuse
it. Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min); land after every milestone with `git pull --rebase origin staging && git push origin HEAD:staging`. Rows in
`results-R129.md` (own file). Long jobs with `nohup python -u` into a log file and an `until` watcher. Disk is about 28 GB and
memory 15 GB: stream the simulations one at a time as R-122 did (`simstream.sh`).

## Trajectory check: already run (2026-11-01, context-free Opus, in the dispatching session), verbatim:
"VERDICT: partly novel: κ×stellar-density/extinction (and PSF) templates propagated into the κ×galaxy bias via w_κS·w_gS/w_SS was
already done for DES Y1/Y3 with SPT+Planck lensing (Omori et al. arXiv:1810.02342; Chang et al. arXiv:2203.12440); what remains
open is the same direct κ×template test on ACT DR6 / Planck PR4 maps with Gaia DR3 magnitude/colour bins and Legacy/KiDS
templates, propagated into the unWISE, DESI LRG and Quaia amplitudes, whose papers rely on galaxy-side weights and mask variation
instead". Other prior work it named: Omori et al. arXiv:1810.02441 (shear x lensing, PSF templates), Giannantonio et al.
arXiv:1507.05551, Krolewski et al. arXiv:1909.07412 (unWISE x Planck; Gaia only galaxy-side), Kim et al. arXiv:2407.04606 and
Sailer et al. arXiv:2407.04607 (ACT DR6 x DESI LRG; mask variation, foreground-only sims), Farren et al. arXiv:2309.05659 (ACT DR6
x unWISE), Alonso et al. arXiv:2306.17748 (Quaia x Planck). Use the DES method (Omori 2019 Sec. V; Chang 2022 App. B1) as the
statistic so the result is like-for-like, and cite it as the method's origin. Do not rerun the trajectory check.

## The question
Milky Way stars, survey depth, PSF, airmass and the scan pattern cannot be physically correlated with CMB lensing, so the
cross-power of the ACT DR6 (and Planck PR4) convergence with those templates measures contamination of the lensing map, and its
propagation through each galaxy sample's own template dependence (w_κS w_gS / w_SS) is the bias on that sample's published
κ×galaxy amplitude. The result is an interval either way: a detected contamination with a sign and size (a correctable
systematic on published S8-from-cross-correlation values) or an upper limit the cross-correlation papers can cite.

## Data (public; verify reachability from the cloud first; zenodo.org and GitHub release assets 403)
ACT DR6 lensing v1 (NERSC portal: baseline, cibdeproj, tonly, polonly, galcut040; analysis masks; the 400 simulations), Planck
PR3 lensing (PLA; PR4 if reachable from the Planck Legacy Archive or NERSC); Gaia DR3 source density maps per G bin and colour
(Gaia archive TAP or the HEALPix count products; chunk queries as R-120 learned); DESI Legacy Survey DR9/DR10 random-catalogue
systematics maps (depth, PSF size, airmass, stellar density, EBV); KiDS-1000 and DES Y3 public survey-property maps where
reachable; Planck dust and CIB maps; the unWISE, DESI LRG and Quaia selection templates where their papers release them.

## Method
1. Fix the template list, the band powers (L ranges), the latitude bands and the pass criterion in `results-R129.md` before any
   cross-spectrum; the null is the same template crossed with each of the 400 ACT simulations (and rotated maps as a check), so
   every band power's null distribution is empirical.
2. Pseudo-C_ell (NaMaster or a healpy equivalent with mode-coupling corrected) of κ with each template on the shared footprint;
   per template, per latitude band, for the baseline and each map variant; the PTE of each against the simulations.
3. Propagation: for each published κ×galaxy sample whose template dependence is public (or measurable from its public catalogue
   and randoms: DESI LRG from the LSS catalogues, Quaia from its released selection function, unWISE from its released maps), the
   DES statistic X_S = <κS><Sg>/<SS> per template, summed in quadrature across templates, as a fraction of the published
   amplitude and of its error.
4. The claim is the table: per template and map, the cross-power's significance and the implied bias on each published
   amplitude with its interval, or the upper limit. No sentence about S8 beyond what the table says.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus
   agent, section alone; hints: the prior works above). Verdict verbatim in `results-R129.md` and under the section. One revision
   pass on a presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-9` in `LOG.md`. Final report: headline, the table, the simulation-null check,
   friction. Budget 6 h.
