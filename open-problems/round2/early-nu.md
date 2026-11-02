# Round 2: early universe, neutrinos and light relics

Quick kill run 2026-11-01. One candidate survives. The arXiv API rate-limited
(HTTP 429) after six queries, so the later checks were done with WebSearch
(extended).

## 1. Primordial helium (Y_p, hence N_eff) from DESI DR1 metal-poor star-forming galaxies

**Question.** What is Y_p when it comes from a homogeneous sample of thousands of
low-metallicity star-forming galaxies, instead of the few dozen hand-picked HII
regions behind every current value? Does it side with the low EMPRESS value or
the standard LBT value?

**Shape.** A disagreement, with a regime extension in sample size.
- Low: EMPRESS XV, arXiv:2506.24050: Y_p = 0.2402 +/- 0.0040, N_eff = 2.54
  (+0.20/-0.25), a mild (~2 sigma) tension with the Standard Model. The earlier
  EMPRESS value was 0.2370 +0.0034/-0.0033 (arXiv:2203.09617).
- Standard: LBT Y_p Project IV, arXiv:2601.22238: Y_p = 0.2458 +/- 0.0013
  (regression gives 0.2451 +/- 0.0022), matching BBN with Planck's baryon
  density (0.2467).
- The tension is ~1.3 sigma against EMPRESS XV and ~2.5 sigma against EMPRESS
  2022. It drives the low-N_eff and lepton-asymmetry papers (arXiv:2610.29069,
  2610.39579, 2603.13226).
- Closest prior survey work: Kurichin et al. 2021 (MNRAS 502, 3045), with about
  100 SDSS HII regions plus HeBCD, got 0.2462 +/- 0.0022. DESI DR2 direct-method
  abundances for 49,959 galaxies (arXiv:2601.02463) give N, O, Ne, S and Ar but
  no helium. DESI EDR found 95 confirmed extremely metal-poor galaxies
  (arXiv:2312.00300, 2408.04019).

**Public data.**
- DESI DR1 spectra and redshift catalogues: https://data.desi.lbl.gov/public/dr1/
  (spectro/redux/iron; zcatalog).
- The DR1 FastSpecFit emission-line VAC, under the DR1 vac/ tree, for
  preselection. Check which He I lines it carries. Plan to refit He I 3889,
  4026, 4471, 5876, 6678 and 7065 ourselves.
- He I emissivities (Porter et al.) and the Aver et al. MCMC helium method
  (arXiv:2010.04180), both published.

**Method.**
1. Select DR1 star-forming galaxies with an [O III]4363 detection, 12+log(O/H) < 8.0,
   EW(Hbeta) > 50 A and S/N(He I 5876) > 20. Get direct-method T_e and O/H.
2. Fit He+/H+ per galaxy with the Aver-style MCMC (n_e, T_e, reddening, underlying
   absorption, collisional and fluorescent terms), using He I 3889 and 7065 to
   constrain density in place of He I 10830, which falls outside DESI's range.
3. Regress Y on O/H (and on N/H and S/H) to get Y_p. Get errors from a bootstrap
   over galaxies plus Monte Carlo over line fluxes.
4. Nulls and injections: split by redshift, fibre aperture fraction and
   spectrograph/petal, and inject synthetic He I lines of known Y into real DESI
   continua to measure the bias of the line fits and absorption corrections.

**Expected precision; can it decide?**
- Statistical error, for a few hundred qualifying galaxies: about 0.0015-0.002.
- Without He I 10830, the density-temperature degeneracy brings back systematics
  of about 0.003-0.005, the pre-2014 Izotov-era level.
- So it probably separates the two answers (a difference of 0.0056) at only
  about 1-1.5 sigma. It would decide clearly only if the injection tests hold
  the systematics below about 0.002.
- Its real worth is as an independent, homogeneous and selection-controlled
  answer, with the slope dY/d(O/H) from many more points. Weak to moderate.

**Quick-kill queries (7; none found Y_p from DESI spectra).**
1. WebSearch "primordial helium abundance DESI dwarf galaxies spectra Yp measurement":
   found LBT IV and SDSS (Kurichin 2021) only.
2. arXiv all:"primordial helium" AND all:DESI: 6 hits, all BAO+BBN combinations or forecasts.
3. arXiv abs:helium AND abs:DESI AND abs:metal-poor: 1 hit (dark-radiation fit, no spectra).
4. arXiv abs:"helium abundance" AND abs:"Dark Energy Spectroscopic Instrument":
   2 hits, both BAO plus literature Y_p.
5. arXiv abs:"metal-poor" AND abs:DESI AND abs:galaxies: 19 hits; DESI EMPG and
   direct-abundance papers, none measuring helium.
6. WebSearch "EMPRESS primordial helium ... LBT Yp ... discrepancy": set out the
   tension (2506.24050 vs 2601.22238). There is no third, survey-based determination.
7. WebSearch "DESI He I 5876 helium abundance star-forming galaxies Y_P": no DESI
   result. Found a statistical reanalysis of published Y_p values
   (arXiv:2606.03192), which uses no new data.

## Dropped in the quick kill

- **PTA background anisotropy x LSS (NANOGrav 15 yr):** already done with 2MPZ
  templates (arXiv:2603.02110) and with DESI, which gave an upper bound. Moving to
  the MPTA southern sky (anisotropy maps in arXiv:2412.01214; PPTA anisotropy in
  arXiv:2602.11529) would give a null far above the expected signal. 3 queries.
- **MPTA 7 nHz hotspot (p = 0.015 before trial factors) tested in other PTAs:**
  not a two-answer disagreement. PPTA DR3 found its own unrelated ~2 sigma
  hotspot, and individual-binary searches already limit those sky positions.
  2 queries, 1 abstract read.
- **PBH fraction from SN Ia magnification (extension to DES 5 yr):** done
  (arXiv:2410.07956, alpha < 0.12). 1 query.
- **Lyman-alpha P1D running/slope vs Planck (eBOSS disagreement) with DESI DR1:**
  done (arXiv:2601.21432; alpha_s = 0.0014 +/- 0.0041, N_eff = 3.02 +/- 0.10).
  Thermal relics from DR1 Lyman-alpha are also done (arXiv:2609.22581). 1 query.
- **Neutrino BAO phase shift in DESI DR1 (beta_phi = 2.7 +/- 1.7, arXiv:2412.05990)
  vs Standard Model 1, extended to the Lyman-alpha forest:** an arXiv check was
  blocked by rate limits. A Lyman-alpha-only sigma(beta_phi) would be several
  units, too weak to decide, so it was dropped on precision. 1 query.
- **f_NL from LoTSS DR2 / RACS scale-dependent bias:** no published value found,
  but sigma(f_NL) of order 50 or more cannot compete with DESI or Planck. The large
  scales it needs overlap our dipole and count-anisotropy work, so it was
  dropped. 1 query.
