# Proposed test: Do early- and late-type galaxies share one weak-lensing radial acceleration relation? (disagreement)

**Question.** At fixed baryonic mass, do isolated early-type galaxies (ETGs) and late-type galaxies (LTGs) lie on the same
lensing RAR (g_obs vs g_bar, 1e-15 to 1e-12 m/s^2)?

**The two answers (same KiDS-1000 data):**
- Brouwer et al. 2021, arXiv:2106.11677: the ETG and LTG RARs differ by **at least 6 sigma** (split by Sersic index and u-r colour).
- Mistele et al. 2023/24, arXiv:2310.15248 (and Lelli et al. arXiv:2411.02499): with a **stricter isolation criterion** and
  consistent stellar and gas masses, ETGs and LTGs lie **on the same** joint RAR. They also say their lenses are probably less
  isolated than the formal criterion implies (photo-z isolation), and that the lensing RAR is robust only down to g_bar ~ 1e-14.

The disagreement turns on isolation, which both teams judged from photometric redshifts. A spectroscopic 3D isolation on
independent shear decides it.

**Public data.**
- DESI DR1 BGS spectroscopic redshifts (lenses and isolation): https://data.desi.lbl.gov/doc/releases/dr1/
- Legacy Surveys DR10 photometry and Tractor Sersic fits (stellar masses, morphology): https://www.legacysurvey.org/dr10/
- DES Y3 shape catalogue (metacalibration) and n(z): https://des.ncsa.illinois.edu/releases/y3a2/Y3key-catalogs
- HSC-SSP Y3 shear catalogue (confirm public access first): https://hsc-release.mtk.nao.ac.jp/
- KiDS-1000 shapes, for a same-shear cross-check: https://kids.strw.leidenuniv.nl/DR4/
- Nulls and calibration: random lens points; a public lightcone mock (AbacusSummit, MICE) with the BGS selection for isolation leakage.

**Method.**
1. Select BGS lenses with no neighbour above 10% of their stellar mass within 3 Mpc projected and ±Δv (spectroscopic); split ETG/LTG by Sersic n and colour exactly as 2106.11677 does, and also as 2310.15248 does.
2. Stack DES Y3 (and HSC Y3) tangential shear around each class in bins of baryonic mass, and convert ΔΣ to g_obs with the Brouwer/Mistele deprojection.
3. Fit the ETG minus LTG offset Δlog g_obs at fixed g_bar with a jackknife-plus-random-point covariance; repeat with photo-z isolation to reproduce both papers' cuts.
4. Calibrate isolation leakage (the 2-halo term) on the mock lightcone; report Δ with a 68% interval.

**Expected precision.** The DES Y3 x DESI DR1 overlap (roughly 700-1400 deg^2, to be measured) plus HSC Y3 (~416 deg^2,
deeper) gives lensing S/N comparable to KiDS-1000's 1006 deg^2. A 6-sigma ETG/LTG offset would show at several sigma;
zero offset (Mistele) would be bounded well below the Brouwer amplitude. **Can tell the answers apart: yes**, if the
overlap is near the upper end. Step one is measuring the overlap.

**Quick kill (6 queries).**
- WebSearch "weak lensing radial acceleration relation isolated galaxies DES Y3 OR HSC OR DESI 2025 2026 ...": KiDS work only (2106.11677, 2310.15248), a Zenodo KiDS reproduction pipeline, and 2607.29314 (a theory fit to the KiDS data).
- WebSearch "radial acceleration relation weak lensing early-type late-type difference DES/HSC/UNIONS/DESI ...": no lensing RAR from anything but KiDS.
- arXiv `"radial acceleration relation" lensing` (13 results) and `"radial acceleration relation" weak lensing`: the lensing RARs are all KiDS-based.
- arXiv "isolated galaxies weak lensing early-type late-type acceleration": only 2106.11677 and 2310.15248.
- arXiv "galaxy-galaxy lensing MOND UNIONS/HSC/DES isolated lenses": 2409.13329 (DECaLS lensing around SDSS centrals for MOND gas masses; no ETG/LTG RAR split); nothing else.
- arXiv "flat rotation curves weak lensing isolated" and "weak lensing isolated galaxies rotation curves": theory only (2504.12662, 2411.18221).

---
