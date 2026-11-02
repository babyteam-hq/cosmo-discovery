# Proposed test: x1. The X-ray x cosmic-shear cross-correlation with eRASS1 instead of ROSAT, AGN masked to the eRASS1 limit

**Question.** Two published ROSAT x DES Y3 shear analyses disagree on the halo mass where feedback has expelled half the gas:
log10(M_c/Msun) = 13.64 +0.08/-0.12 (arXiv:2309.11129, 25 sigma detection) against 14.83 +0.16/-0.23 (arXiv:2412.12081,
X-ray + Planck tSZ x shear). The second paper says modelling of the X-ray AGN contamination (and non-thermal pressure) must
be fixed to resolve the tension. FLAMINGO (arXiv:2602.02484, v5) finds that the strong-feedback models favoured by kSZ and
eROSITA gas-fraction stacks (arXiv:2512.02954, 2410.19905) are incompatible with the ROSAT x DES signal, but only after it
modelled the attenuated, unresolved AGN. The 25-group eRASS1 gas budget (arXiv:2608.17735) also disfavours strong feedback.
The disputed step is the unresolved AGN in a ROSAT map. Does an eRASS1 map with point sources masked to about 3-4x deeper
than ROSAT give M_c near 10^13.6 or near 10^14.8, and does the cross-correlation with the eRASS1 clusters masked (group
regime, 1e13-1e14 Msun, where strong and fiducial feedback differ most) still allow strong feedback?

**Shape.** Disagreement (13.64 vs 14.83, about 5 sigma apart with the published errors; the two fits share the DES shear,
so the gap is in the X-ray modelling) plus a regime extension in mass (FLAMINGO: the ROSAT signal is dominated by
1e14-1e15 Msun halos). Prior work: 2309.11129, 2412.12081, 2602.02484, 2410.22397 (eRASS1 diffuse auto-spectrum,
S8 = 0.80 +0.02/-0.01, which names cross-correlations with lensing as future work), 2608.17735, 2512.02954, 1911.08494
(eROSITA x shear forecast, simulations only).

**Data.**
- eRASS1 (eROSITA-DE DR1) per-tile 0.2-2.3 keV images, exposure maps, and the point-source and cluster catalogues:
  https://erosita.mpe.mpg.de/dr1/ (images and catalogues are public; DR2 is catalogue-only but gives deeper masks).
- DES Y3 shape catalogue (metacalibration) and n(z): https://des.ncsa.illinois.edu/releases/y3a2/Y3key-catalogs
- ROSAT RASS maps for the like-for-like re-run: https://heasarc.gsfc.nasa.gov/docs/rosat/rass.html
- Optional check: ACT DR6 lensing convergence (https://lambda.gsfc.nasa.gov/product/act/actadv_prod_table.html); no X-ray
  map x CMB-lensing measurement exists either (query 6).

**Method (four lines).**
1. Count-rate map (HEALPix nside 1024) from eRASS1 0.5-2 keV images over exposure in the DES overlap; point sources masked
   at three flux limits (none, ROSAT-equivalent, full eRASS1), plus a version with eRASS1 clusters masked to 2 R500.
2. Pseudo-C_ell (NaMaster) against the four DES Y3 tomographic shear bins, 30 < ell < 2000, as in 2309.11129.
3. Errors from 1,000+ shear-rotation nulls (random spin per galaxy: keeps shape noise and mask), Gaussian covariance checked
   against them; B-modes as a null; the ROSAT re-run must reproduce the published 2309.11129 numbers (known answer).
4. Fit the same baryonification model (M_c, Gamma, plus an AGN term with free amplitude) per masking level; the dispute is
   decided by how M_c and the AGN amplitude move from the ROSAT-depth mask to the full eRASS1 mask.

**Overlap and feasibility.**
- eRASS1-DE (Galactic l > 180 deg) x DES Y3: 4,060 deg2 (eRASS1 cosmology papers, arXiv:2402.08458), nearly all of DES Y3
  (4,143 deg2). The ROSAT measurement used the same shear on the same area, so the shape noise is the same.
- eRASS1 resolves point sources to about 5e-14 erg/s/cm2 (0.5-2 keV) at mean exposure (deeper towards the south ecliptic
  pole), against about 1-2e-13 for RASS, so more of the AGN in the cross-correlation are removed rather than modelled.
  FLAMINGO puts the eROSITA-like gas-only cross signal at about 14x the ROSAT one in count units.
- S/N anchor: 25 sigma for ROSAT x DES Y3 (2309.11129). The shear is the same, so eRASS1 should reach at least 25 sigma on
  the full map. With clusters masked, my estimate is 8-12 sigma (clusters carry a half to two-thirds of the signal at
  ell < 1000). That is an estimate, not a forecast. The gate should forecast it from a FLAMINGO-like halo model first.

**Expected separation and anchor.** M_c to about +/-0.1-0.15 dex (the 2309.11129 precision with the same shear), against
a 1.2 dex gap between the two published values: more than 5 sigma once the deeper mask fixes the AGN term. In the
cluster-masked group regime, strong-feedback FLAMINGO groups have about half the fiducial gas fraction (2608.17735 rules out
the strong variants at 2.5-8 sigma on 25 groups), so their X-ray emission (density squared) is about 3-4x lower. A
group-regime detection at 8 sigma under fiducial feedback would put strong feedback at least 4 sigma away.

**What could erase it.**
- Unresolved AGN below the eRASS1 limit still correlate with shear. This is the disputed step, reduced, not removed; the
  test is the trend of M_c and the AGN amplitude with masking depth.
- Galactic foregrounds (local hot bubble, the eROSITA bubbles at the DES edges) and the eRASS1 exposure pattern don't
  correlate with shear. They only add noise, which the null spread checks.
- Solar-wind charge exchange is small in eRASS1 (2410.23345).
- The eRASS1 tile edges and the ~30 arcsec PSF limit ell to below ~3000.
- Scoop risk is real: the ROSAT x DES authors and the eROSITA-DE consortium both name this as the next step, and the
  proprietary eRASS:4/5 data would let the consortium do it deeper.

**Overlap with our work.** We did eRASS1 clusters x Planck lensing (passed, 10 sigma) and parked the Compton-y stacks.
This is a different observable: the diffuse map against galaxy shear, with clusters masked. The round-4 feedback row
(probes-b #4) was dropped because the gas-fraction stacks came down to a modelling choice. Here the disputed step is the
ROSAT AGN term, which the eRASS1 masks change directly.

**Quick-kill queries and results.**
1. WebSearch (extended) "eRASS1 diffuse X-ray background cross-correlation CMB lensing ACT DR6 hot gas": no X-ray map x
   lensing measurement; the eRASS1 APS paper 2410.22397 is auto-only and names cross-correlations as future work.
2. WebSearch (extended) "eROSITA X-ray map cross-correlation weak lensing shear baryonic feedback detection 2025 2026":
   only ROSAT x DES (2309.11129), forecasts, and stacked-profile work (2512.02954).
3. WebSearch (extended) "eRASS1 X-ray map cosmic shear cross-correlation DES Y3 KiDS eROSITA diffuse emission baryons
   Ferreira": no eRASS1 x shear map measurement; FLAMINGO 2602.02484 compares only with ROSAT x DES.
4. arXiv API abs:"cosmic shear" AND (abs:X-ray OR abs:eROSITA), all 31 entries: every measured X-ray x shear
   cross-correlation is ROSAT x DES (2309.11129, 2412.12081; compared in 2602.02484).
5. arXiv API abs:eRASS1 AND (abs:shear OR abs:"cross-correlation" OR abs:"cross-power"): 8 entries, all cluster catalogue
   x shear (mass calibration), none with the diffuse map.
6. arXiv API abs:"X-ray" AND abs:"CMB lensing" AND abs:cross: 5 entries, none a measurement.
