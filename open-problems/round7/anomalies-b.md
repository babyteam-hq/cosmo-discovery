# Round 7, harvester B: unretested anomaly claims (tracers, high z, transients, small-scale DM, fundamental physics)

2026-11-02. Search and reasoning only; no code run, no data downloaded. Target: a claim at >= 3 sigma (2018-2026) that
nobody has re-tested on independent data, where public data now allow a re-test that separates "real" from "chance" at
>= 3 sigma (the Giant Arc shape, R-155). Skipped as taken: Giant Arc / Big Ring / Huge-LQG, CMB cooling, spin asymmetry,
quasar polarisation alignment, sky dipoles (incl. alpha dipole, Migkas, bulk flows).

## Survivors: none

Every candidate below already has an independent-data re-test, never reached 3 sigma, or cannot be re-tested at 3 sigma
on public data. Searches: 10 extended web searches and 7 arXiv API queries (listed per row).

| Candidate (claim, arXiv, significance) | Why it is out | Evidence |
|---|---|---|
| Hercules-Corona Borealis Great Wall in GRBs (2008.03679; 1510.01933) | Weak (p = 0.017 in 2312.10050), disputed by isotropic Monte Carlo (2006.00141), and already tested on independent tracers: SDSS volume-limited quasars at 1.6 < z <= 2.1 (7,847 QSOs) find homogeneity ("Large-scale homogeneity in the distribution of quasars in the HCB region", 2022). A GRB overdensity also doesn't fix a quasar amplitude, so no frozen template. | web (extended) + arXiv API `all:"Hercules-Corona Borealis"` |
| Galaxy rotation coherent with neighbour motion to 6 Mpc (Lee et al. 1908.10972, CALIFA) | 2.8 sigma only; re-tested on SAMI (2207.03752): 1.7 sigma, attributed to scatter | web (extended) |
| CHIME FRB x galaxy correlation and high-host-DM subpopulation (2106.04354); FRB DM excess through foreground halos (Connor & Ravi 2107.13692) | Re-tested on Catalog 2: 2506.08932 (5.1 sigma DM x galaxies), 2604.22105 (3,455 sightlines, ten tracers; overdense sightlines carry larger DM), 2601.21336 | web (extended) + map row 10-transients #6 |
| Quasar Hubble diagram > 4 sigma off LCDM (Risaliti & Lusso 1811.02590) | eRASS1 (2508.00732, 519 QSOs) too shallow in z to test it; SN-calibrated re-analysis (2408.15547) finds no tension; the dispute is the non-evolution assumption, not the data | web (extended) |
| MOND external field effect in rotation curves, > 4 sigma (Chae 2009.11525, 2109.04745) | Re-tested on independent LITTLE THINGS curves (2605.27217: weak evidence against); 2112.00026 shows LCDM produces an EFE-like signal, so a re-test cannot separate the readings | arXiv API `"external field effect" AND rotation` |
| Screened fifth force from HI-optical centroid offsets, 6.6 sigma (Desmond 1807.01482) | Withdrawn by the same authors after upgrades (2009.08743: no evidence on 0.3-8 Mpc); significance depended on the noise model; FAST/WALLABY beams too coarse for the predicted offsets | web (extended) + arXiv API `ti:"fifth force"` |
| Rotation-measure excess behind Mg II absorbers, 3.7 sigma (1710.10396) | Done on independent data: SPICE-RACS DR2 x DESI DR1, 4.5 sigma (2605.16924) | web (extended) + arXiv API `abs:"rotation measure" AND Mg II` |
| Excess strong Mg II toward blazars/FSRQs | Null on new samples (2507.17866) | arXiv API |
| Vast Polar Structure orbital alignment (Pawlowski) | New dwarfs already tested against the plane (2609.02245, 2209.12422); 2205.02860 argues consistency with LCDM | arXiv API `"vast polar structure"` |
| BBH mergers with AGN flares (GW190521/ZTF19abanrhr, Graham 2023) | Re-tested on GWTC-4 x ZTF DR23 (Zhu et al.; 7 +24/-5 % associated) and ALeRCE O4 search (2603.04342) | web (extended) |
| Black-hole cosmological coupling k = 3 (2302.07878) | Many independent tests (Gaia BH binaries, JWST AGN, LIGO) | reasoning (known literature) |
| Fornax dwarfs lack dark matter (Asencio 2208.02265) | Depends on subjective "disturbed" classification; no frozen statistic transfers to another cluster; same-data follow-up 2502.04091 | web (extended) |
| Overmassive high-z BHs (Pacucci 2308.12331) | Re-measured on larger JWST AGN samples; selection-bias rebuttals; no clean frozen template | reasoning |
| ISW excess of supervoids A_ISW = 5.2 +/- 1.6 (1811.07812) | ~2.6 sigma; CMB area (harvester A) | map 12-foundations |
| Andromeda satellite lopsidedness; Cen A co-rotating plane | ~3 sigma on 30-40 objects; too few new satellites for a 3 sigma re-test | reasoning |
| EDGES 21 cm, 3.5 keV line, GRB LIV, wide binaries, Keplerian MW decline, neutrino-TDE | Re-tested (SARAS 3, LHAASO, Gaia DR3 rebuttals, many MW analyses) or outside 2018-2026 | reasoning |

## Note for the next round

In these five areas, the >= 3 sigma claims since 2018 with positional templates have mostly been re-tested within two
years by the survey teams (CHIME Catalog 2, SPICE-RACS x DESI, GWTC-4 x ZTF). What's left is single objects (no
independent data) or model-mediated claims (no frozen statistic). The next independent data that reopen this are
Gaia DR4 (wide binaries, VPOS, December 2026) and Euclid DR1 (cluster-dwarf morphology across clusters, LRD pairs).
