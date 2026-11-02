# Results (R-145, cosmo-deep-16, 2026-11-02)

Question: do isolated early- and late-type galaxies share one weak-lensing radial acceleration relation (RAR)? This run
uses DESI DR1 BGS lenses with spectroscopic isolation and DES Y3 shear. The two published answers disagree on the same
KiDS-1000 data:
- Brouwer et al. 2021 (arXiv:2106.11677) find a difference of at least 6σ. log⟨g_obs,E/g_obs,L⟩ is 0.17 dex with a
  Sérsic split and 0.27 dex with a u−r split, over the full g_bar range. At g_bar > 10⁻¹³ m s⁻² it is 0.14 dex and 0.19 dex.
- Mistele et al. 2024 (arXiv:2310.15248) find one relation. They use 4 Mpc isolation and their own masses: Q_ETG = 1.4, hot
  gas for ETGs and cold gas for LTGs, which makes ETG baryonic masses about 1.7 times Brouwer's.
- Milgrom 2013 (arXiv:1305.3516) did an earlier red/blue split of CFHTLS lenses (Brimioulle+13). Red lenses need a baryonic
  M/L about 0.6 dex higher at fixed luminosity.

Trajectory verdict on file (brief): "partly novel": the independent-shear, spectroscopic-isolation measurement is what is
open. Scripts: `r145_isolate.py`, `r145_forecast.py`, `r145_photoz_leak.py`. Data are under `data/r145/` (gitignored).

## Step 0: feasibility gate. FAILED; stopped as the brief says

**Fetch times (data/ started empty).**

| Product | Size | Time | Notes |
|---|---|---|---|
| DESI DR1 BGS group catalogue (`vac/dr1/bgs-groups/v1.0/bgs-groups-dr1-v1.0.fits`) | 290 MB | 21 s | Williams et al. 2026. BGS BRIGHT r < 19.5, 3-pass footprint, 2,229,741 galaxies with fastspecfit stellar masses, rest-frame g−r and a QUIESCENT flag. 79 % DESI spectra, 1 % SDSS, 12 % nearest-neighbour z, 1 % photo-z. |
| DES Y3 GLIMPSE mass-map mask (`y3a2_files/massmaps/glimpse_mask.fits`) | 100 MB | 15 s | Used as the footprint: nside 1024, 4,742 deg². |
| DES Y3 KS map | 100 MB | 10 s | Unmasked, so not used. |
| DES Y3 metacal shapes (`y3kp_cats/DESY3_metacal_v03-004.h5`) | 312 GB | not fetched | Range requests are served (206). First byte takes 0.6-110 s at a cold offset, with several timeouts in 10 tries. The HDF5 metadata walk got 3 reads in ~25 min and never reached the column layout. The whole file is 11× the disk. A column-wise read is possible in principle, as large sequential ranges after the layout is known, but was not needed once the forecast failed. |

**Lens counts in the DES Y3 × DESI DR1 BGS overlap.** Candidates have a spectroscopic redshift, 0.05 < z < 0.30,
9.5 < log M*/M☉ < 11.0 (h70, from LOGMSTAR + 0.31) and r < 19.5: there are 46,015.

A neighbour counts against a lens if it has more than f times the lens stellar mass and lies within R projected. For a
neighbour with a spectroscopic redshift there is also a window: |Δv| < 500 km/s. A neighbour without a spectrum is placed
at the lens redshift, with its mass rescaled by D_L², so the fill-in is conservative.

| Isolation | Isolated lenses | ETG (quiescent) / LTG |
|---|---|---|
| R = 3 Mpc, f = 0.1 (Brouwer's radius and fraction) | **7** | 3 / 4 |
| R = 3 Mpc, f = 0.1, Δv < 1500 km/s | 7 | |
| R = 3 Mpc, f = 1.0 (no more massive neighbour) | 236 | |
| R = 2 Mpc, f = 0.1 | 79 | |
| R = 1 Mpc, f = 0.1 | 2,399 | 727 / 1,672 |
| R = 1 Mpc, f = 1.0 | 6,760 | 2,646 / 4,114 |
| R = 0.5 Mpc, f = 0.1 | 12,343 | 3,899 / 8,444 |
| none | 46,015 | 16,450 / 29,565 |
| none, wide (0.05 < z < 0.5, 9.0 < log M* < 11.6) | 84,371 | 39,397 / 44,974 |
| R = 1 Mpc, f = 0.1, wide | 8,187 | 3,978 / 4,209 |

BGS-wide, outside the DES cut, the 3 Mpc / 10 % criterion leaves 41 of 906,484 candidates. The conservative fill-in removes
most lenses that spectroscopic neighbours alone would pass. With spectroscopic neighbours only, about 9 % pass at 3 Mpc,
which would be about 4,000 lenses in the overlap. That is still below what the gate needs (next table).

**Forecast (`r145_forecast.py`).**
- **Noise:** shape noise only, so it is an optimistic bound. DES Y3 has n_eff = 5.59 arcmin⁻², σ_e = 0.261, and a smooth
  n(z) with mean 0.63. Sources are taken at z_s > z_l + 0.1.
- **Signal and binning:** Brouwer's 15 g_bar bins, with the statistic over the 7 bins at g_bar > 10⁻¹³. g_obs = 4GΔΣ, with
  the expected signal from the McGaugh RAR.
- **Mass and split:** Brouwer's baryonic mass. The type split is the QUIESCENT flag, standing in for the colour split.
- **Statistic:** Δ = the mean over those bins of log(g_obs,ETG/g_obs,LTG).

| Lens set | σ(Δ), dex | Brouwer Sérsic 0.14 dex | Brouwer colour 0.19 dex |
|---|---|---|---|
| R = 3 Mpc, f = 0.1 (7 lenses) | 4.0 | 0.0σ | 0.0σ |
| R = 1 Mpc, f = 0.1 | 0.30 | 0.5σ | 0.6σ |
| R = 1 Mpc, f = 1.0 | 0.14 | 1.0σ | 1.4σ |
| R = 0.5 Mpc, f = 0.1 | 0.13 | 1.1σ | 1.5σ |
| none (46,015) | 0.063 | 2.2σ | 3.0σ |
| none, wide (84,371) | 0.038 | 3.7σ | 5.0σ |

**Gate verdict: fail.** No isolated sample in the overlap separates "same relation" from Brouwer's difference at 3σ, even
from shape noise alone.
- Only lenses with no isolation reach 3σ. Those are the samples where 2-halo terms dominate the comparison, and that is not
  the question.
- At the isolation radius both papers use (3-4 Mpc, 10 %), spectroscopic redshifts leave 7 lenses.
- At 1 Mpc, σ(Δ) is 0.14-0.30 dex. That is as large as the effect.
- The constraints: the overlap of DESI DR1's 3-pass BGS footprint (NTILE_ALT ≥ 3) with DES Y3 is about **190 deg²**
  (occupied nside-256 pixels; 98 % of lenses at −7° < dec < +4°). The proposal assumed 700-1,400 deg². With the BGS flux
  limit on top, that gives ~20 times fewer lenses than KiDS-bright's ~10⁶.
- HSC Y3 shapes need a registered account. KiDS-N overlaps DESI, but its shear is not independent.

## By-product: how isolated are photo-z-isolated lenses? (`r145_photoz_leak.py`, no shear, no gate)

This turns the dispute's own question into a measurement with the spectroscopic redshifts already in hand. It is not the
mock calibration the referee asked for.

**Method.** Brouwer's isolation is applied to BGS with every redshift scattered by 0.02(1+z), the KiDS-bright ANNz2 scatter.
- Masses are rescaled by D_L².
- The criterion is no neighbour with more than 10 % of the lens mass within 3 Mpc in 3D.
- Lens cuts: 0.1 < z_phot < 0.5, 10 < log M* < 11, z_spec < 0.3.
- For each photo-z-isolated lens, the test is whether a real neighbour with more than 10 % of its mass lies within R
  projected and 500 km/s.
- Only spectroscopic neighbours are used in either view, 60,000 lenses per seed.

**Results** (fraction of photo-z-isolated lenses that are not isolated spectroscopically; seeds 1 and 2):

| log M* | Type | N (seed 1) | within 0.5 Mpc | within 1 Mpc | within 3 Mpc |
|---|---|---|---|---|---|
| 10.0-10.5 | ETG (quiescent) | 1,996 | 0.48 / 0.48 | 0.69 / 0.70 | 0.94 / 0.93 |
| 10.0-10.5 | LTG | 8,341 | 0.30 / 0.29 | 0.53 / 0.53 | 0.88 / 0.89 |
| 10.5-11.0 | ETG | 10,205 | 0.38 / 0.39 | 0.60 / 0.60 | 0.90 / 0.91 |
| 10.5-11.0 | LTG | 11,986 | 0.27 / 0.27 | 0.49 / 0.48 | 0.87 / 0.86 |

- Statistical errors are 0.003-0.011.
- 54 % of lenses pass photo-z isolation. 9 % pass spectroscopic isolation at 3 Mpc and 40 % at 1 Mpc.
- At fixed stellar mass, photo-z-isolated ETGs have a real neighbour with more than 10 % of their mass within 0.5 Mpc
  11-18 points more often than LTGs (0.38-0.48 vs 0.27-0.30).
- That is the direction that raises the ETG ΔΣ at the R ~ 0.1-0.5 Mpc of g_bar ~ 10⁻¹³, which is where Brouwer's offset
  sits. It also fits Mistele's finding that ETGs, not LTGs, move as R_isol grows.

**Caveats.**
- BGS (r < 19.5) is shallower than KiDS-bright (r < 20).
- QUIESCENT is not KiDS u−r.
- Gaussian scatter has no catastrophic outliers.
- The size of the ΔΣ bias this implies needs the neighbours' halo masses (or a mock). It is not computed here.

## What would make the measurement feasible

- DESI DR2/DR3 BGS: more passes and a wider overlap with DES, including dec < −10°.
- Public HSC Y3 shapes.
- A softer, stated isolation (no more massive neighbour within 1 Mpc, f = 1) on a lens sample about 5 times larger.
- Or the question can be reframed without shear: the type-dependent isolation leakage above, turned into a ΔΣ bias with a
  mock lightcone (AbacusSummit BGS) or the BGS group catalogue's halo masses. That would be a shear-free test of whether
  Brouwer's 0.14-0.19 dex is leakage. It needs its own trajectory check.
