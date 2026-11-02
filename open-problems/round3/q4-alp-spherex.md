# Proposed test: eV-mass axion-like-particle (ALP) dark-matter decay line in SPHEREx spectral images of dwarf spheroidals

**Question.** Does SPHEREx, an instrument independent of JWST, see a two-photon decay line from ALP dark matter
(m_a = 0.5-3.3 eV, line at lambda = 2.48 um / (m_a/eV)) in the dwarf spheroidals? If not, what is the upper limit
on g_agamma, and does it reach the coupling band where the two JWST analyses disagree?

**Shape.** Mainly a disagreement, with a regime extension: the first decay-line search on SPHEREx data, and the
first all-sky 0.75-5 um spectral survey used for this.
- **Strong answer:** Pinetti, arXiv:2503.11753, *First constraints on QCD axion dark matter using JWST observations*,
  from the public JWST blank-sky archive. It claims limits that exclude the QCD axion down to about 0.1 eV.
  Janish & Pinetti, arXiv:2310.15395 (PRL 134, 071002), is the earlier analysis of the GN-z11 sky spectrum.
- **Weak answer:** arXiv:2609.10674, *Re-examining the sensitivity of JWST to decaying axion dark matter*. It does
  its own reduction with a Gaussian-process continuum, data-driven error rescaling and a spurious-signal correction,
  and gets NIRSpec limits of "few x 1e-11 GeV^-1". It says the strongest claim is about 10x stronger in coupling
  (about 100x in flux) and "exceeds [the sensitivity] attainable by any analysis of the datasets from which [it was]
  derived". Saha et al., arXiv:2503.14582, agrees with it to within 10% before its correction.
- **The disputed step:** JWST instrument modelling and pipeline background handling. 2609.10674 shows that its own
  sensitivity saturates at about 25 observations, because the systematic mismodelling grows with the data. More
  JWST data or another JWST reanalysis cannot settle it. An instrument with a different calibration chain, the
  dSph's own D-factor contrast, and an ON-OFF within single exposures removes that step.
  - Caveat: our mapping of 2609.10674's reference numbers to arXiv ids comes from the order of its reference list
    ([49] 2310.15395, [63] 2503.11753, [68] 2311.04987, [69] 2503.14582). Confirm it before citing.
- **Prior SPHEREx work:** forecasts only, arXiv:2412.12286 (dSphs, LMC, MW halo in the deep fields; full 2-yr
  mission; no flat-field or stellar-light systematics) and 2102.00580 (cross-correlation with LSST). No analysis
  of SPHEREx data exists.

**Public data.**
- SPHEREx Quick Release spectral images (QR1-QR3; QR3 has been at IRSA since 2026-09-16; weekly releases since
  July 2025). Their units are MJy/sr. The zodiacal model is supplied in a separate ZODI extension and **not**
  subtracted, and the per-pixel wavelength is in the CWAVE/CBAND images:
  https://irsa.ipac.caltech.edu/Missions/spherex.html , https://registry.opendata.aws/spherex-qr/ ,
  docs https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/ . The IRSA Mosaic Tool
  (2026-05-07) builds spectral cubes.
- dSph decay D-factors: Song, Hiroshima & Murase 2024, JCAP 05, 087 (the set 2412.12286 uses). Use Pace & Strigari
  2019 as the cross-check.
- Star masks: the SPHEREx source catalogue, plus Gaia DR3 / Legacy Surveys DR10 for resolved RGB stars.
- Our SPHEREx tooling from the tomography work (`spherex_*.py`) already reads QR images.

**Method.**
1. For Draco, Ursa Minor, Sculptor, Fornax, Sextans, Carina, Leo I/II and Bootes I, cut every QR exposure
   touching r < 1.5 deg. Mask stars. Keep each pixel's wavelength (CWAVE) and bandwidth (CBAND), and difference ON
   annuli against OFF pixels **in the same exposure and LVF row** (the same wavelength), which cancels the
   flat-field and zodi terms.
2. In each wavelength channel, fit the ON-OFF surface brightness against radius with three profiles: the
   predicted D(r) (NFW or cored), the dSph's stellar profile (Plummer, from published structural parameters; it is
   more compact than the dark matter) times a stellar-population template spectrum, and a planar zodi gradient.
3. Profile-likelihood amplitude of the line (R = 40 Gaussian at fixed lambda) on a grid of m_a. Look-elsewhere
   over about 100 channels. Stack the dwarfs with their D-factors (+/- their D-factor errors) as nuisances.
4. Errors and nulls:
   - errors from the same fit on 200+ blank fields at matched ecliptic latitude (fake dSphs) and on rotated or
     scrambled annuli;
   - injection of a line at known g into real QR pixels, to measure the bias of the full chain;
   - the result is g_agamma(m_a) with a 95% interval, or a detection.

**Overlap (from footprints).**
- SPHEREx is all-sky, and by QR3 it has made at least 2 full-sky passes (one per 6 months since May 2025), so all
  nine dSphs are covered.
- Draco (about 9 deg from the NEP) and Ursa Minor (about 17 deg) get above-average coverage. They dominate the
  forecast; each alone is about 1.4x weaker than the combination (2412.12286).
- Each dSph's ON region (r < 0.5 deg) is about 2.6e5 SPHEREx 6.2" pixels per channel, before masking.

**Feasibility (signal, what could erase it, S/N).**

*The signal.* The line brightness is I = Gamma D_Omega / 4pi, with Gamma = g^2 m^3 / 64pi.
- Draco: D(0.5 deg) of about 6e18 GeV cm^-2 (an assumed order of magnitude; take the exact value from Song+2024),
  so D_Omega is about 2.5e22 GeV cm^-2 sr^-1.
- At m_a = 2 eV and g = 1e-12 GeV^-1, I = 1.9e-4 nW m^-2 sr^-1. Spread over one R = 40 channel, that is
  vI_v of about 7.6e-3 nW m^-2 sr^-1.
- The line scales as g^2: 0.19 nW at g = 5e-12, and 6.8 nW at 3e-11.

*The noise, and what reaches the disputed band.*
- Per-pixel, per-channel noise is about 27 nW m^-2 sr^-1 at 0.75 um after the full 2 yr (2412.12286). With the
  roughly 1.3 yr now public it is about 33 nW.
- Stacking 2.6e5 pixels gives 0.065 nW, so Draco alone reaches a 2-sigma g of about 4e-12 at m_a = 2 eV.
  Combining the dSphs gives about 3e-12.
- Both are below 2609.10674's few x 1e-11, and at the level of 2503.11753's claim.
- Because g scales as m^-3/2, the reach weakens to about 1.1e-11 at m_a = 1 eV, still inside the disputed band.
  At 0.5 eV it is about 3e-11: no discrimination there.

*What could erase it.*
1. **Sky-derived flat fields or background subtraction.** In a linear-variable filter (LVF) instrument, a
   sky-uniform line at fixed lambda is degenerate with a per-pixel flat. The same step removed the line in some JWST
   pipelines. The ON-OFF within one exposure and LVF row, plus the D(r) radial shape, is immune to this. An
   all-sky Milky Way halo analysis would not be: it is statistically 10x deeper (5e-13 in the NEP deep field) but
   systematics-bound, so it is kept only as a secondary result.
2. **dSph starlight.** Draco's centre is about 24 AB mag/arcsec^2 in the near-IR, about 20% of the zodi there and a
   few % averaged over 0.5 deg. It must be modelled to about 0.5% of the stellar flux at R = 40 (CO 2.3 um and H2O
   bands). The stellar profile is more compact than D(r), and the line test sits in the 0.3-1 deg annulus, where
   starlight is below about 1% of the zodi. This is the main risk; the injections measure it.
3. **ISM lines** (Pa-alpha, Br-alpha, PAH 3.3 um) at the dSphs' |b| > 30 deg. Mask the m_a channels they occupy.
4. **D-factor uncertainty** (about +/- 0.2 dex for the classical dSphs) moves g by about 25%. It cannot close a 10x
   gap.

**Expected precision; does it decide?**
- At m_a of about 1.5-3 eV it does decide. A null excludes g of about 3-5e-12 with data independent of JWST, which
  makes 2503.11753's claimed exclusion physically secure whatever its JWST analysis did. A line at 5e-12 to 3e-11
  would mean the strong JWST claim was wrong and the conservative reading of 2609.10674 was right.
- Below about 1 eV it gives only an independent few x 1e-11 limit.
- First-year QR calibration maturity is the unforecast risk. The blank-field nulls show whether the floor sits at
  the photon noise or above it.

**Quick-kill queries (no SPHEREx-data ALP or decaying-DM line search found).**
1. WebSearch: "SPHEREx axion-like particle decay line search eV dark matter Milky Way halo spectral images
   constraint". Found only forecasts: 2412.12286, 2102.00580, 2012.00771, and the EPJC LIM x LSS forecast.
2. arXiv API: `abs:SPHEREx AND (abs:axion OR abs:"dark matter" OR abs:decay)`, 12 hits. None uses survey data for
   DM; the data papers are commissioning, persistence and thermospheric He.
3. arXiv site search, "SPHEREx dark matter decay": 3 hits (2502.03328 LORRI COB, 2412.12286, 2102.00580), all
   forecasts or other data.
4. WebSearch: "JWST blank-sky spectra axion decay line limit ... NIRSpec". This set out the JWST disagreement
   (2503.11753 / 2310.15395 vs 2609.10674, with 2503.14582 and 2311.04987) and found no SPHEREx result.
5. WebSearch: "SPHEREx public data release 2026 ...". QR3 is out (2026-09-16), the all-sky cubes are not yet
   released, and the only early analysis is the spectral map-making preprint.
