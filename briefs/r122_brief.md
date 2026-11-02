# cosmo-deep-2: the eRASS1 mass scale with ACT DR6 CMB lensing (R-122)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top, and `### cosmo-deep-1`), `research/cosmo-discovery/blind-review.md`,
the sixth section of `finding.md` (the Planck version of this measurement, with its three referee verdicts verbatim and their
objections), `results-R121.md`, and `erass1_cmblens.py` (R-121's script: reuse its catalogue handling, model, filters and fits).
Follow the agents' working rules: worktree on staging; keep-alive (a background loop running the session keep-alive
every 20 min for the whole job); land after every milestone with `git pull --rebase origin staging && git push origin HEAD:staging`.
Rows in `results-R122.md` (own file). Long jobs with `nohup python -u` into a log file and an `until` watcher; logs to files.

## What this is
R-121 stacked Planck PR3 CMB lensing on the eRASS1 clusters: 10σ, one-halo mass scale 0.8-1.0 ± 0.2 of eROSITA's weak-lensing
calibration. Three blind referees called it new and all three said "not sound", the last one because the mass-scale claim was
stronger than its intervals (every 95 % b_M interval includes 0.5) and the firmer exclusion fixed the two-halo term to the eRASS1
cosmology under test. Planck cannot do better. ACT DR6 can: about three times the signal over the eRASS1 footprint it covers.

## Data (reachable from the cloud; LAMBDA is not, NERSC is)
ACT DR6 lensing v1, public: https://portal.nersc.gov/project/act/dr6_lensing_v1/maps/<variant>/ with variants `baseline`,
`cibdeproj`, `tonly`, `polonly`, `galcut040`; files `kappa_alm_data_act_dr6_lensing_v1_<variant>.fits` (92 MB),
`mask_act_dr6_lensing_v1_healpix_nside_4096_<variant>.fits`, `N_L_kk_...txt`, `kappa_filter_...txt`, and
`simulations/kappa_alm_sim_act_dr6_lensing_v1_<variant>_<i>.fits` (400 simulated reconstructions). Read the release's README
in that tree first (conventions: the alm are already mean-field subtracted and normalised; the filter file is the band). Cite
Qu et al. 2024 (arXiv:2304.05202), Madhavacheril et al. 2024 (arXiv:2304.05203), the maps paper. The eRASS1 catalogues from the
eROSITA site as R-121 did. Everything under `data/` (gitignored); never run a download as code; `python -I` for readers.

## Method (paper depth; every referee objection from the Planck note answered by construction)
1. Baseline map, with the release mask; the same footprint selection as R-121 (clean within 30', ≥ 80 % clean within 2°), and
   report how many eRASS1 clusters survive, by mass bin, and how many of those the mask's holes remove at fixed mass (the referee's
   point: point-source holes correlate with clusters).
2. Model as R-121's, with two changes the referees required: the predicted κ per cluster averaged over the cluster's mass posterior
   (the catalogue's M500_L/H as a lognormal, and the M500_PDF column as a check, both reported), and the two-halo term fitted free
   with a halo-exclusion form (or at least its amplitude free and its σ8 varied 0.76-0.88 as a systematic), never fixed to the eRASS1
   cosmology.
3. Covariance and nulls from the 400 simulations (the same stack on each simulated map, so the covariance carries the map's own
   noise and mask), plus mask-matched randoms (random positions passing the same footprint selection) for the null-mean
   correction, plus the rotated-catalogue nulls R-121 used, all three reported; Hartlap on the sims.
4. The statistic: the one-halo amplitude with the two-halo free, per sample (primary in footprint; M500 < 4e14; the cosmology
   sample; cosmology below 4e14), and b_M with its 68 % and 95 % intervals. Checks: cibdeproj (the CIB-leakage objection), tonly,
   polonly, galcut040; concentration ×0.7/×1.3; the Planck PR3 stack on the same clusters for comparison.
5. The claim is the interval, nothing beyond it: if the 95 % b_M interval excludes 0.5, say so and only that; if it excludes 0.8,
   say that. No morphology split, no S8 sentence stronger than what b_M allows. Trials: fix the samples and statistic before the
   first fit, record every amplitude computed.
6. Write the result as a `finding.md` section (seventh), then the finding check (blind-review.md's after-writing prompt, verbatim,
   fresh context-free Opus agent, the section alone in a file; hints for closest prior work: anyone who stacked CMB lensing on
   eROSITA or X-ray-selected clusters, or tested the eRASS1 mass scale with an independent observable). Verdict verbatim in
   `results-R122.md` and under the section. One revision pass if the objection is a presentation point; park on a physics point.
   Nothing goes on the site; if it passes, the headline of your report says so and the main session writes the paper.
7. End: `### cosmo-deep-2` in `LOG.md` (numbers, verdict, what would settle it). Final report: headline
   (passed or not, verdict verbatim), the measurement in one table (sample, n, A_1h, A_2h, b_M 68/95, detection), the sim-vs-rotation
   null comparison, friction. Budget 6 h; stop at the budget with everything landed.
