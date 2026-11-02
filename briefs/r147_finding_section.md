# Photometric-redshift isolation leaks unevenly by galaxy type: on DESI DR1 BGS spectroscopic truth, the KiDS-bright isolation criterion leaves red galaxies with a real > 10 % neighbour within 0.3 Mpc 1.5 times as often as blue ones (30 % against 20 %), mostly as hidden group satellites, and that predicts a spurious early/late offset in the lensing radial acceleration relation of 0.04 dex at g_bar > 10⁻¹³ m s⁻², 0.21-0.26 of the published 0.14-0.19 dex, rising to 0.41-0.51 of it over the full range

## Claim

Brouwer et al. 2021 (arXiv:2106.11677) found isolated KiDS-1000 early- and late-type galaxies ≥ 6σ apart in the weak-lensing
radial acceleration relation (RAR). Their measure is log10⟨g_obs,E/g_obs,L⟩: 0.14 dex (Sérsic split) and 0.19 dex (u−r split) at
g_bar > 10⁻¹³ m s⁻² (R < 0.3 Mpc/h70), and 0.17 and 0.27 dex over all 15 g_bar bins. They assumed photo-z errors degrade the
isolation of both types equally. Mistele et al. 2024 (arXiv:2310.15248) find one relation when they use 4 Mpc/h70 isolation and
their own ETG masses.

Here the KiDS-bright isolation (no neighbour with > 10 % of the lens stellar mass within 3 or 4 Mpc in 3D, from photo-zs) is
emulated on the DESI DR1 BGS group catalogue. The emulation uses the published KiDS-bright photo-z error distribution, five seeds
and about 780,000 candidate lenses per seed. Each photo-z-isolated lens is then checked against the spectroscopic redshifts.
Each real > 10 % neighbour within R_isol and 500 km/s has its group halo (NFW, catalogue halo mass) added to the lens's ESD in
Brouwer's g_bar bins. The lens's own profile is taken as the same for both types (the McGaugh relation), so any ETG−LTG offset
that comes out is spurious by construction.

- **Contamination (3 Mpc, all masses, mass-matched as Brouwer).** Within 0.3 / 0.5 / 1 / 3 Mpc, 0.30 / 0.43 / 0.63 / 0.92 of
  photo-z-isolated red (rest-frame g−r > 0.75) lenses have a real > 10 % neighbour, against 0.20 / 0.31 / 0.53 / 0.89 of blue
  ones. The gap is 0.10-0.19 at ≤ 0.5 Mpc in every Brouwer stellar-mass bin from 10^8.5 to 10^11 M☉ and shrinks at 3 Mpc.
  - The 4 Mpc criterion barely changes these numbers (0.27 / 0.39 against 0.18 / 0.28).
  - The cause: 28 % of photo-z-isolated red lenses are satellites in a host group, against 12 % of blue ones.
- **Predicted spurious offset log10⟨g_obs,E/g_obs,L⟩.** Median with 68 % and 95 % ranges, over seeds × bootstrap with halo-mass
  draws:

| Split (Brouwer's published value) | Criterion | g_bar > 10⁻¹³ | fraction | all 15 bins | fraction |
|---|---|---|---|---|---|
| colour, rest g−r (u−r: 0.19 / 0.27) | 3 Mpc | 0.041 (0.038-0.046; 0.036-0.048) | 0.22 (0.20-0.24) | 0.121 (0.111-0.136; 0.109-0.140) | 0.45 (0.41-0.50) |
| colour | 4 Mpc | 0.040 (0.036-0.045; 0.034-0.048) | 0.21 (0.19-0.24) | 0.112 (0.104-0.126; 0.100-0.131) | 0.41 (0.39-0.47) |
| spectral QUIESCENT (u−r) | 3 Mpc | 0.048 (0.044-0.054; 0.043-0.057) | 0.26 (0.23-0.29) | 0.130 (0.120-0.143; 0.117-0.148) | 0.48 (0.44-0.53) |
| spectral QUIESCENT | 4 Mpc | 0.047 (0.042-0.051; 0.039-0.056) | 0.25 (0.22-0.27) | 0.122 (0.111-0.135; 0.107-0.140) | 0.45 (0.41-0.50) |
| morphology DEV vs EXP/REX (Sérsic: 0.14 / 0.17) | 3 Mpc | 0.031 (0.029-0.035; 0.027-0.039) | 0.22 (0.20-0.25) | 0.086 (0.080-0.091; 0.076-0.096) | 0.51 (0.47-0.54) |
| morphology | 4 Mpc | 0.032 (0.028-0.037; 0.024-0.040) | 0.23 (0.20-0.27) | 0.083 (0.077-0.091; 0.072-0.096) | 0.49 (0.45-0.53) |

- **Radius dependence (colour, 3 Mpc).** The predicted ETG−LTG offset in log g_obs is 0.009 dex at R = 0.04 Mpc, 0.048 at 0.12,
  0.088 at 0.21, 0.14 at 0.37 and 0.19-0.22 at 0.7-2 Mpc. The predicted ETG bias itself runs from 0.02 dex to 0.56 dex.

The claim is this table and the radius dependence. It says type-dependent leakage accounts for about a fifth to a quarter of
the published offset where Brouwer judged isolation reliable (a third for the depth-complete subsample), and for up to half of it
over the full range. It does not say what the remainder is.

## Data and method

- **Data.** DESI DR1 BGS group VAC (bgs-groups-dr1-v1.0; r < 19.5, 2.23 M galaxies). It gives fastspecfit stellar masses,
  rest-frame g−r, the spectral QUIESCENT flag, group membership, central/satellite and an abundance-matched group halo mass.
  Legacy Survey MORPHTYPE comes from the DESI DR1 LSS BGS_BRIGHT catalogue. H0 = 70, Ωm = 0.3, physical distances.
- **Lenses.** Spectroscopic, 0.05 < z < 0.3. In the emulated photo-z view, Brouwer's cuts apply: 0.1 < z_ph < 0.5 and
  log M*_ph < 11, in his bins [8.5, 10.3, 10.6, 10.8, 11.0]. All galaxies with a redshift serve as neighbours.
- **Photo-z emulation.** Bilicki et al. 2021's (arXiv:2101.06010) fit to KiDS-bright − GAMA, a generalized Lorentzian with
  a = 2.613 and s = 0.0149, drawn per galaxy. Its wings are the outliers. Masses rescale by D_L² with 0.12 dex scatter, as in
  Brouwer's MICE emulation.
- **Isolation.** Applied in the photo-z view in 3D (projected separation plus photo-z line-of-sight distance), at 3 Mpc
  (Brouwer) and 4 Mpc (Mistele), with 10 % of the stellar mass.
- **Truth.** Catalogue redshifts, a projected radius and |Δv| < 500 km/s.
- **ESD model.** Each contaminating neighbour's group halo is an untruncated NFW (Duffy et al. 2008 c(M, z)) at its catalogue
  mass, taken as M_200m. It is placed at the group central's projected offset and azimuthally averaged about the lens.
  - A lens that is a satellite gets its host halo.
  - A lens that is a central with a > 10 % satellite gets its group halo minus the median halo of catalogue centrals without
    satellites at its stellar mass.
- **RAR conversion.** g_obs = 4GΔΣ, Brouwer's SIS conversion, in his 15 log g_bar bins 10⁻¹⁵-5 × 10⁻¹². Each lens enters at
  R = (G M_gal/g_bar)^½ with Brouwer's cold-gas M_gal. Weights per lens and bin are (R/D_A)² Σ_crit⁻² for a KiDS-like source
  n(z).
- **Statistic.** Types are mass-matched by reweighting to the overlap of their M*_ph histograms. The statistic is
  δ = log10(mean over bins of (1 + b_E)/(1 + b_L)), Brouwer's statistic, with b the contamination ESD over the shared own
  profile.
- **Intervals.** 5 seeds, each with a fresh photo-z draw and fresh halo masses (per halo N(0, 0.2) dex plus a global N(0, 0.15)),
  each with 100 Poisson bootstraps over lenses.
- **Fixed before any number was computed** (results-R147.md, timestamped): sample, splits, error model, criteria, statistic,
  ESD model, intervals, checks, and what each outcome would mean. Nothing was changed after. One choice was made at
  implementation: the morphology catalogue has no Sérsic index, so its SER class (65 % of rows) is dropped as the
  pre-registration said.

## Checks

| Check (colour, 3 Mpc) | δ at g_bar > 10⁻¹³ / all bins | Fraction |
|---|---|---|
| Spectroscopically isolated lenses | 0 (asserted for every run) | — |
| Gaussian σ_z = 0.02(1+z) | 0.043 / 0.128 | 0.23 / 0.47 |
| Spectroscopic neighbours only | 0.040 / 0.119 | 0.21 / 0.44 |
| No scatter (the 3D criterion's own redshift-space leak) | 0.012 / 0.044 | 0.06 / 0.16 |
| Lenses r < 17 (complete to 10 % neighbours in BGS) | 0.059 (0.041-0.073) / 0.099 | 0.31 / 0.37 |
| Same, spectral split | 0.066 (0.050-0.084) / 0.106 | 0.35 / 0.39 |
| Per Brouwer mass bin | 0.038-0.045 / 0.11-0.18 | — |

## What it means

- Brouwer's equal-contamination assumption fails. At fixed stellar mass, photometric isolation leaves red, quiescent or
  bulge-dominated galaxies with a real massive neighbour 1.2-1.6 times as often as blue ones inside 0.5 Mpc, mainly because
  they are more often unrecognised satellites.
- At R < 0.3 Mpc, where both KiDS papers trust their isolation, this bias reaches the published offset only in part: 0.03-0.05
  dex of 0.14-0.19 dex. Most of the small-radius split is not leakage of > 10 % neighbours. That is consistent with Mistele et
  al.'s attribution of much of it to the ETG mass-to-light ratio and hot gas, but it does not test that attribution.
- At larger radii, leakage predicts an offset of order 0.2 dex, comparable to the published one. So the full-range 0.27 dex (u−r)
  is contaminated to roughly half.
- The predicted leakage barely changes between 3 and 4 Mpc, because the photo-z line-of-sight error (~100 Mpc) dwarfs both
  radii. Leakage of > 10 % neighbours therefore does not reproduce the R_isol dependence of the ETG RAR that Mistele et al.
  report.

## Weakest points

- **The own profile.** The fractional bias is taken relative to a type-independent McGaugh profile, without the environment
  beyond R_isol.
  - Over all lenses at R > 0.3 Mpc, the contamination adds 1.05 times that profile. Brouwer's MICE test found 0.30 relative to
    the truly isolated total ESD.
  - The two agree only if the environment beyond R_isol is about 2.5 times the own profile there.
  - If it is less, the large-radius biases, and the full-range fractions, are overestimates. At R < 0.3 Mpc the denominator
    question is minor.
- **Halo masses.** The abundance-matched group masses, the satellite assignments of the group finder, untruncated NFW halos and
  Duffy concentrations carry the ESD. The ±0.2 dex per-halo and ±0.15 dex global draws are the only modelling uncertainty
  propagated.
- **Sample differences.**
  - DESI BGS (r < 19.5, z_spec < 0.3, median z 0.20) is shallower and closer than KiDS-bright (r < 20, z_ph < 0.5, mean 0.25).
    The depth-complete r < 17 subsample raises the fraction to 0.31-0.35. Higher redshift lengthens the photo-z line-of-sight
    error, so both differences probably make the main numbers low.
  - Rest-frame g−r and QUIESCENT stand in for observer-frame u−r. DEV vs EXP/REX on 40 % of lenses stands in for Sérsic n ≷ 2.
- **Omitted terms.** Neighbours under 10 % of the lens mass and fibre-collision redshift assignments are not modelled as
  contaminants.

## To settle it

- A lightcone mock with a KiDS-bright selection to r < 20 and z < 0.5 (e.g. AbacusSummit or MICE) run through the same
  pipeline. Its truly isolated total ESD would fix the denominator, by type.
- A KiDS-bright × GAMA (or DESI DR2 BGS in KiDS-N) spectroscopic subsample. The isolated lensing RAR measured there with the
  two isolation selections would measure the bias directly.
