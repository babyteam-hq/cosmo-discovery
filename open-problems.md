# Open problems: how candidates are found and chosen (living; refine when progress stalls)

Will, 2026-11-01: explore the open-problem space in an unbiased way, extract and validate the "most open" ones we can legitimately
contribute to, and refine this process in an ongoing way when progress stalls.

## Stages

1. **Map (unbiased).** A fixed taxonomy of twelve areas covering cosmology, chosen before looking at what we have tried (below).
   One harvester per area reads primary sources (2024-2026 reviews' open-question sections, tension papers, "future work"
   statements, decadal and Snowmass lists, forecasts never carried out, anomalies reported and not chased) and lists the open
   problems with the arXiv id that states each one open. Our past attempts are given only so overlap is flagged, never as a filter.
2. **Extract.** For each problem, the harvester names the measurement on public, reachable data that could move it, and why it
   looks unpublished (what it searched). A problem with no measurement we could make is kept in the map, marked "not ours".
3. **Quick kill.** The harvester's own search; anything clearly done is marked and dropped there.
4. **Validate.** The top candidates across areas (ranked by openness, data readiness, payoff) go to the blind trajectory check in
   `blind-review.md`, prompt unchanged. Only a verbatim "novel" is dispatched.
5. **Execute.** One worker per "novel", with a feasibility gate in its first hour; the finding check at the end.

## Yield log (refine the stage that stalls; one line per round)

| Round | Date | Problems mapped | Measurements extracted | Survived quick kill | Referee "novel" | Dispatched | Passed gate | Change made after |
|---|---|---|---|---|---|---|---|---|
| screen-1 | 2026-11-01 | (ad hoc, 7) | 7 | 7 | 1 | 1 (R-137) | **pass** (after one revision, R-140) | ideation anchored on past passes |
| screen-2 | 2026-11-01 | (dataset scout, 3) | 3 | 3 | 0 | 1 exception (R-138) | miss, revising (R-143) | dataset-first yields "same question, new data"; switch to problem-first (this file) |
| map-1 batch 1 | 2026-11-01 | 12 areas, ~170 rows | ~150 | 24 top picks | 0 of 4 (A ZTF fs8: answered at 6 %; B CHIME lensing FAR: calibrated before; C LRD variability: answered at N<=9; D dark-force dipole: never measured, but DESI too weak) | 0 | - | Stated open problems are posed by definition, so they come back 'partly novel'; rule refined: posed-but-never-measured with discriminating precision qualifies. Subagents reported a 'shared session' web-search limit; checked: the parent's search works, so the cap is per agent. Fix: each agent spends its searches on targeted queries and uses the arXiv listing/API (export.arxiv.org) for the rest; batches run in parallel. |
| map-1 batch 2 | 2026-11-01 | (same map) | | | 1 of 4 (G quasar bias z~4, novel; E not novel; F answered elsewhere; H posed, precision unforecast) | 1 (R-139) | miss, parked | Full searches (arXiv API past the cap) gave sharper verdicts than batch 1's degraded ones. |
| map-1 batch 3 | 2026-11-01 | (same map) | | | 0 of 4 (I, J, K, L all answered on other data) | 0 | - | Stalled stage: the quick kill (harvesters' searches capped, 'unpublished' rested on one search). Round 2: strict quick kill on the arXiv API before listing; target regime extensions (a known measurement never taken to a redshift/mass/scale/sample public data now reaches) and disagreements between published answers, the shape of the one novel (G). Map-1 totals: 12 refereed, 1 novel. |
| map-2 batch 1 | 2026-11-01 | 6 harvesters x 2 areas | ~40 | 10 (strict: >=3 arXiv queries each) | 2 of 4 (z~5 lensing novel; WL RAR settles a disagreement) | 2 (R-145 stopped at its feasibility gate, ~190 deg2 overlap; R-146) | pending | The strict quick kill and the extension/disagreement target lifted the yield from 1 in 12 to 2 in 4. |
| map-2 batch 2 | 2026-11-01 | (same map) + R-145 by-product | | | p5 isolation leakage (settles a dispute), p6 Lya x y (novel), p7 dust z~1 (posed, never measured), p8 streams (epicyclic test only), p9 quasar tSZ (held); R-145 stopped at its gate | 2 (R-147, R-148) + 3 queued | pending | Feasibility gates stop bad bets in ~1 h (R-145); a footprint check now goes in each proposal. |
| execution stall | 2026-11-01 | | | | | R-143, R-146, R-147, R-150 missed 'novel but not sound'; R-145, R-148 stopped at gates | | Stalled stage moved from novelty to soundness. Change: a blind method check on each pre-registration before computing (blind-review.md 1b; in every brief from R-153). |
| map-3 | 2026-11-01 | 6 harvesters (overlap + feasibility checks) | ~45 screened | 6 (one harvester returned 0 of 15) | 1 of 4 qualifies (q4 ALP SPHEREx); q1 SZ lensing held on precision; q2, q3 answered | 1 (R-153) | pending | The feasibility checks moved kills from the worker to the harvest; the referee now kills mostly on 'answered with other data'. |
| map-4 (+ Giant Arc from map-3) | 2026-11-01 | 3 harvesters, 3-sigma floor anchored on published significance, no local code | ~25 screened | 4 (probes B 0) | Giant Arc novel; Doppler magnification qualifies; rotational kSZ, relativistic dipole no | 2 (R-155, R-156) | **R-155 pass** (first check, with the method check); R-156 pending | Survivors are thinning: most of what public data allows is done or below 3 sigma until Euclid DR1 / eROSITA DR2 photons. |
| map-5 (Euclid day-one) | 2026-11-02 | 2 harvesters | | 3 | spin asymmetry (R-158), LRD lensing (held for DR1); LRD census no | 1 + 1 held | pending | Euclid DR1-Foundation slips into November. |
| map-6 (less-mined data types) | 2026-11-02 | 2 harvesters | ~60 | 4 | polarisation alignment (R-159); eRASS1 x shear held; RM pairs no | 1 | pending | |
| map-7 (unretested anomaly claims) | 2026-11-02 | 2 harvesters | ~22 | 1 | jet alignment regions novel | 1 | pending | The Giant Arc shape; survey teams re-test most anomalies within ~2 years, so the window is narrow. |
| map-8 (Gaia DR4 day-one; claims C) | 2026-11-02 | 2 harvesters | ~15 | 1 | tidal tails (held for DR4) | 0 now | - | Public-data candidates are exhausted for now: rounds 5-8 put 4 workers through, all stopped at gates or parked. No round 9 until Euclid DR1 or Gaia DR4 is public; two briefs wait for them. |

## Taxonomy (fixed for round 1)

1. Expansion history, distance ladder, H0
2. Dark energy and the DESI w0-wa signal
3. CMB: primary, polarisation, reionisation optical depth, large-angle anomalies
4. Large-scale structure growth and S8
5. Clusters, groups and baryons: feedback, missing baryons, the circumgalactic medium
6. Dark matter on small scales: Milky Way, Gaia, dwarfs, streams, lensing substructure
7. Early universe: inflation, primordial non-Gaussianity, the PTA background, primordial black holes
8. Neutrinos and light relics: mass sum, N_eff, the negative-mass preference
9. High-redshift galaxies and reionisation: JWST, 21 cm, quasars
10. Transients as probes: GW sirens, FRBs, supernovae, lensed transients, time delays
11. Black holes across cosmic time: SMBH growth, little red dots, AGN demographics
12. Foundations and anomalies: isotropy, constants, gravity tests, unexplained catalogue features

## Tried by us (for overlap flags only)

Sky dipoles and anisotropy in counts; SN residuals vs dust; Cosmicflows-4 flows and SN vs CF4 fields (passed); eRASS1 clusters x
ACT/Planck lensing (passed) and Compton-y stacks (parked); group tSZ on ACT (parked); void lensing; FRB DM x y; GW host bias;
lensed SNe in ZTF (running); CMB-lensing template systematics; Quaia x lensing; JWST massive-galaxy tallies; supercluster mass
outside members (running); SPHEREx NIR tomography (running).
