# cosmo-deep-5: FRB dispersion-measure excess against the Compton-y maps: the temperature of the baryons the FRBs find (R-125)

Goal (the only one): a novel, sound cosmological result on public data, reproducible from `research/cosmo-discovery/`.
Read first: `LOG.md` (handoff at the top), `research/cosmo-discovery/blind-review.md` (the gate; trajectory
check before any compute, finding check after writing), `finding.md` (what passed and what missed; the misses were claims stronger
than their intervals, and R-124 was killed as "not novel" after a "partly novel" trajectory verdict: if the trajectory check says
partly novel, the open part must stand alone as a result or you stop). Follow the agents' working rules: worktree on staging; keep-alive (a
background loop running the session keep-alive every 20 min); land after every milestone with
`git pull --rebase origin staging && git push origin HEAD:staging`. Rows in `results-R125.md` (own file). Long jobs with
`nohup python -u` into a log file and an `until` watcher; logs to files.

## The question
A fast radio burst's extragalactic dispersion measure is the electron column along its sightline; the thermal Sunyaev-Zeldovich
Compton-y is the electron pressure column. Their ratio, measured on the same sightlines, is a density-weighted electron temperature
of the gas the FRBs pass through. FRB-baryon work so far correlates DM with galaxy or halo counts; tSZ work stacks on galaxies and
clusters; a catalogue large enough to regress one on the other at the sightline level only came with CHIME/FRB Catalog 2 (about
3,600 bursts, 2026). Positive: a slope of DM excess on y (or of y on DM excess, with the selection in mind) gives that temperature and
says whether the gas the FRBs find is the 10^6-10^7 K warm-hot medium. Null with a detected DM excess: most of those electrons are
cool (below about 10^5.5 K) diffuse gas, outside tSZ's reach. Either is a measurement with an interval.

## Data (public; verify reachability from the cloud first; zenodo.org, LAMBDA data directories and GitHub release assets 403)
CHIME/FRB Catalog 2 (assets.chime-frb.ca; R-120's note: the catalogue lives at catalogs/catalog2/releases/<date>/data/, the
bucket listing is paginated), one burst per repeating source (R-119's lesson: repeaters inflate everything), Galactic DM from
NE2001 and YMW16 (pygedm) plus a halo term, |b| > 20 or so; Planck 2015 MILCA and NILC y maps (PLA, as R-121 used; files over ~1 GB
are cut by the proxy, resume with curl -C -), ACT DR6 Compton-y (NERSC's ACT portal https://portal.nersc.gov/project/act/ is
reachable: look for the DR6 ILC y maps and their CIB-deprojected variants), the Planck CIB and dust maps for nulls.

## Method
1. Trajectory check first (blind-review.md's before-computing prompt, verbatim, one fresh context-free Opus agent): hints: anyone
   who cross-correlated or regressed FRB dispersion measures against tSZ or Compton-y maps, or measured the DM-y relation or the
   temperature of FRB-traced gas (search "FRB dispersion measure tSZ", "FRB Compton-y cross-correlation", "DM y relation", 2022-26;
   include CHIME Catalog 1 work and ASKAP/CRAFT localised-FRB stacks). "Partly novel" with a thin open part is a no: stop and report.
2. The estimator: for each burst, DM_exc = DM - DM_MW(NE2001 or YMW16) - DM_halo (state the halo model; carry both NE2001 and YMW16);
   y sampled from the map at the burst position and in annuli (CHIME localisations are arcminutes to a degree: use the catalogue's
   per-burst localisation ellipse and average y over it, and show the result is stable to the aperture); regress DM_exc on y with
   the catalogue's selection (S/N, declination, exposure) as covariates, and y on DM_exc. Expected scaling: for gas at temperature T,
   y / DM ~ (k_B T / m_e c^2) x (sigma_T) x units; write the conversion out and check it on a cluster sightline.
3. Nulls that respect the structure: y re-sampled at the burst positions rotated in RA (declination and CHIME exposure preserved),
   the same with the CIB-deprojected y variants, the regression against the Planck CIB and dust maps (a positive slope there is
   contamination), and the regression restricted to bursts with DM_exc above and below the median.
4. The claim is the slope with its interval and the implied temperature interval, or the upper limit; and separately the DM
   excess detection itself. No sentence about feedback or the missing baryons beyond what the interval says.
5. Write the result as a `finding.md` section, then the finding check (after-writing prompt verbatim, fresh context-free Opus
   agent, section alone; hints as in step 1). Verdict verbatim in `results-R125.md` and under the section. One revision pass on a
   presentation point; park on a physics point. Nothing on the site; a pass is the headline of your report.
6. End: `### cosmo-deep-5` in `LOG.md`. Final report: headline, the table (map, N bursts, slope, interval,
   implied T, nulls), friction. Budget 6 h; stop at the budget with everything landed.
