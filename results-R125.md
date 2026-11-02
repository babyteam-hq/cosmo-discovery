# Results (R-125, cosmo-deep-5, 2026-11-01)

The CHIME/FRB Catalog 2 dispersion-measure excess against the Compton-y maps: the density-weighted electron temperature of
the gas FRBs find. This file is the worker's own. **Killed at the trajectory check ("not novel"); nothing was downloaded or
computed**, per `blind-review.md` ("a 'not novel' verdict kills the row before any download") and the brief ("'partly novel'
with a thin open part is a no: stop and report").

## Trajectory check (before any compute; `blind-review.md`, context-free Opus, 2026-11-01)

| proposal | verdict (verbatim) |
|---|---|
| DM_exc = DM − DM_MW (NE2001 and YMW16) − DM_halo per burst, one burst per repeating source, \|b\| > 20; Planck 2015 MILCA/NILC and ACT DR6 y (and CIB-deprojected variants) averaged over the per-burst localisation ellipse; regressions of y on DM_exc and DM_exc on y with S/N, declination and exposure as covariates, the slope converted to a density-weighted T_e; RA-rotation, CIB/dust-map and DM-split nulls; the mean-y excess on FRB sightlines; expected 3σ at T ≳ 10^6 K with N ≈ 2,000-2,500 | VERDICT: not novel: Takahashi et al. arXiv:2511.02155 already detected the FRB DM × Planck/ACT Compton-y correlation and converted it to an electron temperature (~2×10⁷ K), and Sharma et al. arXiv:2604.22105 and Zhai et al. arXiv:2608.09014 already measured CHIME Catalog 2 DM × Planck y (3.8σ; WHIM T ≈ 2.4×10⁶ K) |

The proposal as the referee saw it is the brief's method written as a sweep row (question, data, method, expected precision).

Verified afterwards on arXiv (abstract pages fetched; all six exist and say what the referee says):

- Takahashi, Ioka, Shirasaki & Osato 2025 (arXiv:2511.02155, PASJ): first detection of the angular cross-correlation between
  cosmological DM and Compton-y, 133 localised FRBs over 1'-1000' against the Planck and ACT y maps; A = 2.01 ± 0.50 (4.0σ,
  Planck) and 1.23 ± 0.82 (1.5σ, ACT) relative to the Planck 2018 ΛCDM halo-model prediction (HMx); isothermal reading
  ⟨T_e⟩ ≈ 2 × 10^7 K. This is the brief's headline quantity (y/DM as a temperature) on both of the brief's map families.
- Sharma, Krause, Ravi, Anbajagane, Connor et al. 2026 (arXiv:2604.22105, ApJ): 3,455 unique CHIME/FRB Catalog 2 sources
  (~15' localisations) cross-correlated with ten tracers, among them tSZ (3.8σ), the CIB (4.0σ), soft X-ray background
  (4.1σ), X-ray clusters (5.0σ); tSZ × DM and SXRB × DM amplitudes match a moderate-feedback model at 0.5σ, weak feedback
  disfavoured at 3.5σ by SXRB × DM. The same catalogue and sightlines the brief would have regressed, with the CIB null
  included as a tracer.
- Zhai, Wang, Dai & Cen 2026 (arXiv:2608.09014): CHIME/FRB Catalog 2 DMs × the Planck tSZ y map with virialised clusters
  masked; positive correlation at > 99.77 %; joint inference f_WHIM = 0.48 (68 % 0.27-0.61) at a mean WHIM temperature of
  2.4 × 10^6 K. Answers "is it the warm-hot medium" directly.
- Sharma, Ravi et al. 2026 (arXiv:2605.01994): 3,455 CHIME sightlines stacked on 1,288 BOSS voids, a DM deficit at 3.2σ;
  with the existing tSZ void stacks, T_e ≲ (1.1 ± 0.7) × 10^6 K for void gas: the DM-plus-y temperature logic applied to voids.
- Muñoz & Loeb 2018 (arXiv:1809.04074): the original forecast of DM × tSZ as the WHIM's temperature probe (per-burst DM and y,
  ~10^4 CHIME bursts).
- Medlock & Nagai 2026 (arXiv:2608.06455): the Takahashi and Sharma detections re-read with feedback and non-thermal pressure
  separately parameterised.

What the referee named as not yet published, and why it does not stand alone: the per-burst regression of y on DM_exc (both
directions, errors in variables) is a re-parameterisation at zero lag of the covariance those papers measure; the pairing of
CHIME Catalog 2 with ACT DR6 and the CIB-deprojected y variants may be unpublished but have partial precedents (Takahashi on
ACT, Sharma on the CIB); the RA-rotation and dust-map nulls are robustness checks, not a question. Decision (agent, under the
brief's rule): stop; no download, no compute. R-124's lesson holds: a brief written from an ideation list, not from a
literature harvest, met a 2025-26 literature that had already done it; the FRB-baryon field is now moving at several papers a
quarter on this exact catalogue (R-115 had already closed the FRB-LSS rows on arXiv:2506.08932 and 2604.22105).

## Reachability notes (checked before the verdict as the brief asked; nothing downloaded)

- assets.chime-frb.ca: Catalog 2 at `catalogs/catalog2/releases/20251231/data/` (chimefrbcat2.csv 4.1 MB, .fits 7.6 MB,
  .json 11 MB, .npy 24 MB), listing answers 200.
- Planck Legacy Archive (`pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=`): 200 for
  `COM_CompMap_Compton-SZMap-milca-ymaps_2048_R2.00.fits`, `COM_CompMap_Compton-SZMap-nilc-ymaps_2048_R2.00.fits`,
  `COM_CompMap_CIB-GNILC-F545_2048_R2.00.fits`, `COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits`; 404 for the
  `COM_CompMap_YSZ_R2.02` bundle names.
- NERSC ACT portal: `https://portal.nersc.gov/project/act/dr6.02/` answers 200 (the DR6 release directory; the y maps would
  be under it); LAMBDA's ACT product table answers 404 from this network.
- Python: the container has numpy and pandas only; a venv with scipy, astropy, healpy, statsmodels installs cleanly, but
  `pygedm` fails to build a wheel on Python 3.13 (would need the Fortran/C extension toolchain or a Python ≤ 3.12 venv).

Budget used: about 25 minutes of 6 hours.
