# Round 2: transients and foundations, regime extensions and disagreements

Mapped 2026-11-01. Quick kill: WebSearch (extended) plus arXiv search. The arXiv API (`export.arxiv.org`) returned
HTTP 429 for this whole session (other agents share it), so the arXiv queries went through arXiv's HTML search
(`arxiv.org/search/?searchtype=all&order=-announced_date_first`), which covers all years, 2018-2026 included.
Each query and what it found is listed per candidate. "No hit" means none in these queries.

**Result: 3 survivors, all in foundations (tests of gravity, catalogue anomalies). None survived in transients.**

---

## 1. Do early- and late-type galaxies share one weak-lensing radial acceleration relation? (disagreement)

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

## 2. Gravitational redshift of galaxies in clusters with DESI DR1 (regime extension)

**Question.** The stacked member-minus-centre velocity offset Δv(R) in clusters, and its radial profile, with errors small
enough to separate GR from f(R) and to see the relativistic terms of arXiv:2503.11585.

**Prior work it extends.** Wojtak et al. 2011 (arXiv:1109.6571, SDSS DR7, about -10 km/s); Sadeh et al. 2015
(arXiv:1410.5262, SDSS+BOSS); Jimeno et al. 2015 (arXiv:1410.6050); Mpetha et al. 2021 (arXiv:2102.11156, SPIDERS,
~2,500 clusters, -13.5 ± 4.7 km/s, which says ~500,000 galaxies are needed to separate GR and f(R)); Rosselli et al. 2022
(arXiv:2206.05313, SDSS, 3,058 clusters at z < 0.5, -11.4 ± 3.3 km/s, marginal tension with strong-field f(R)).
**All are SDSS depth (main sample r < 17.77, or BOSS), all at z < 0.5, all ±3-5 km/s.** DESI DR1 BGS (r < 19.5) reaches
about 2 mag fainter members, and DR1 LRGs carry members to z ~ 0.8. The DESI DR1 PV paper (arXiv:2512.03230) says
gravitational redshift isn't yet accounted for; arXiv:2610.37782 forecasts infall (F_G), not gravitational redshift.

**Public data.**
- DESI DR1 redshifts (BGS, LRG) and DR1 mocks: https://data.desi.lbl.gov/doc/releases/dr1/
- Cluster centres: DESI Legacy Imaging cluster catalogues (Zou et al.; Wen & Han), redMaPPer (http://risa.stanford.edu/redmapper/), eRASS1 clusters for X-ray centres (https://erosita.mpe.mpg.de/dr1/).
- Spectroscopic BCG redshifts: SDSS DR17 (https://www.sdss4.org/dr17/) plus DESI DR1.

**Method.**
1. Match DESI DR1 members to clusters with spectroscopic BCGs; Δv = c(z_gal - z_BCG)/(1+z_BCG) within 0-4 Mpc, binned in R and in mass (richness, L_X).
2. Fit the Δv distribution as a cluster Gaussian plus interloper model (as in 2206.05313) for the mean shift per bin.
3. Model GR with the transverse-Doppler, light-cone and surface-brightness terms (2503.11585; Kaiser arXiv:1303.3663), and f(R) and DGP profiles.
4. Errors: bootstrap over clusters. Nulls: shuffled BCG redshifts and random centres. Mis-centring check with X-ray vs BCG centres (arXiv:2603.25339 warns mis-centring can mimic modified gravity).

**Expected precision.** Fainter BGS members give several times the SDSS member count per cluster, so the ~500k-member target
of 2102.11156 is within reach at z < 0.4 if DR1 fibre completeness in cluster cores holds up. DR1 tiles are incomplete in
dense regions, and that's the main risk. My estimate is ±1.5-2 km/s overall, about twice as tight as Rosselli, which
separates GR from the strong-field f(R) profile that 2206.05313 already found marginal. A z = 0.5-0.8 LRG-member bin,
a regime not reached before, would be about ±4-6 km/s. **Can tell apart: GR vs strong-field f(R), probably yes; GR vs DGP, no.**
Scoop risk: a DESI group gravitational-redshift paper was "in prep" in 2024 (Lahav, Geneva talk); none found published.

**Quick kill (5 queries).**
- arXiv "gravitational redshift galaxy clusters DESI" (29 results): no DESI measurement; forecasts only (2610.37782).
- arXiv `"gravitational redshift" galaxy clusters` (~50 results, all years): measurements are SDSS/BOSS/SPIDERS only (1109.6571, 1410.5262, 1410.6050, 2102.11156, 2206.05313); the rest are theory or forecasts (2603.25339, 2601.19861, 2512.13221, 2503.11585).
- arXiv `"gravitational redshift" clusters galaxies stacked`: theory only since 2014.
- WebSearch "measurement gravitational redshift galaxy clusters DESI DR1 BGS LRG stacked velocity offset ...": no DESI measurement; SDSS/SPIDERS only.
- Abstracts read: 2102.11156, 2206.05313, 2603.25339.

---

## 3. CMB "cooling" around nearby galaxies, filaments and voids: an independent test (disagreement)

**Question.** Is there an achromatic CMB temperature decrement within the halos of nearby large spirals (z < 0.02), and an
excess in nearby voids, beyond ΛCDM plus known foregrounds?

**The two answers.**
- Claimed by one group: Luparello et al. 2022 (arXiv:2206.14217), Hansen et al. 2023 (arXiv:2305.00268), Garcia Lambas et al.
  2023 (arXiv:2310.13755, which explains the Cold Spot as this foreground), Hansen et al. 2024 (arXiv:2411.15307: 5.7 sigma,
  3-4 sigma after look-elsewhere), Dominguez Feldman et al. 2025 (arXiv:2506.08833, filaments, published in PRD) and Hansen
  et al. 2025 (arXiv:2506.08832: hot voids at z < 0.03, ~6.5 sigma, accepted in A&A). All of these use Planck PR3
  SMICA/SEVEM with 2MRS tracers and 100-10,000 simulations.
- ΛCDM with known foregrounds: no such signal. The tSZ of these halos is ~1-2 uK and chromatic (2506.08833 itself puts its
  100-217 GHz difference down to tSZ). Nobody outside the group has published a test, confirmation or rebuttal.

Strictly, only one side is a published measurement. This is a claimed >5 sigma detection set against the standard model,
never checked outside the group.

**Public data.**
- Planck PR4 (NPIPE) maps and their 600 end-to-end simulations: https://portal.nersc.gov/project/cmb/planck2020/ ; PR3 component maps and FFP10 sims: https://pla.esac.esa.int/
- ACT DR6 maps (f090/f150/f220), for sub-degree halo profiles and a second instrument: https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/
- Tracers: 2MRS (http://tdc-www.harvard.edu/2mrs/) to reproduce the claim; new, independent tracers from the Siena Galaxy Atlas 2020 (sizes of ~380k large galaxies, https://www.legacysurvey.org/sga/sga2020/) with DESI DR1 BGS / 6dF redshifts; DESIVAST DR1 voids (https://data.desi.lbl.gov/doc/releases/dr1/vac/desivast/).

**Method.**
1. Freeze the claimed profile and cuts of 2411.15307 (spirals with R > 8.5 kpc, z < 0.015, pixels inside the halo) as the hypothesis, before looking at the data.
2. Apply them to galaxies **not in 2MRS** (SGA with DESI/6dF redshifts) and to DESIVAST voids, so the sky positions are new.
3. Stack ΔT profiles on Planck PR4 component-separated and single-frequency maps and on ACT DR6 f090/f150/f220; fit an amplitude A (claimed profile = 1) and the SED (achromatic vs tSZ vs dust).
4. Errors from the 600 PR4 sims (plus many Gaussian ΛCDM realisations), with the same masks; null: rotated tracer positions.

**Expected precision.** The CMB sky is shared, so a second instrument tests only processing artefacts; the independent test
is the new tracer positions. If SGA roughly triples the large-spiral sample at z < 0.02 outside 2MRS (still to be checked),
a real signal gives A = 1 at about 3 sigma, and chance or selection gives A < 0.3 at 95%. Frequency differences cancel the
primary CMB, so the SED test is limited by noise and foreground residuals, not by cosmic variance: well under 1 uK.
**Can tell apart: yes** for the halo claim; the void claim is weaker, because the local voids are few and large.

**Quick kill (6 queries, plus 3 papers read).**
- arXiv "CMB cooling galaxies foreground spiral": only the group's 2411.15307.
- arXiv "CMB temperature nearby galaxies 2MRS": only the group's 2310.13755 and 2206.14217.
- arXiv "CMB temperature decrement nearby galaxies Hansen Luparello": only 2310.13755.
- arXiv "local voids CMB temperature excess sign change ISW": only 2506.08832 itself.
- WebSearch "CMB temperature decrement around nearby large spiral galaxies Luparello Hansen independent test rebuttal foreground ACT SPT": group papers only; no rebuttal, no ACT/SPT test.
- WebSearch "independent test CMB cooling nearby galaxies filaments ... ACT DR6 OR Planck PR4 OR no evidence 2025 2026": no outside test.
- Read: the 2506.08833 HTML (PR3 SMICA, SEVEM frequency check, 100-1000 sims) and the 2411.15307 and 2506.08832 abstracts.

---

## Dropped in the quick kill

| candidate | shape | why dropped | queries |
|---|---|---|---|
| T_CMB(z) from SZ clusters to z ~ 2 with ACT DR6 | extension | ACT DR5 already used 370 clusters to z = 1.4 (arXiv:2106.12467, α = 0.017 +0.029/-0.032), limited by inter-band calibration; the joint α is -0.001 ± 0.012. DR6 adds little, and there's no disagreement to settle | 1 WebSearch (arXiv API blocked); dropped on precision |
| Δα/α from Mg II/Si IV alkali doublets in DESI absorbers | extension | Low-resolution centroids plus the 2-3e-5 DESI wavelength-calibration floor (arXiv:2404.03123) can't reach the 1e-6 of high-resolution work or test the Webb dipole | 1 WebSearch |
| Δα/α from DESI DR1 [O III], with a sky test | extension | The pilot 2404.03123 already tested spatial variation; JWST covers 2.5 < z < 9.5 (2405.08977) and X-Shooter 1.1 < z < 3.7 (2310.17947); calibration-limited | 3 arXiv |
| Supervoid ISW excess with DESI DR1 LRG supervoids | disagreement (Kovács 1811.07812, A ~ 5, vs Hang 2105.11936, A = -0.10 ± 0.69) | eBOSS quasar supervoids already took it above z = 0.8 (Kovács 2021, arXiv:2107.13038); DESI shares its sky and CMB with BOSS/Legacy, so it isn't an independent third method; overlaps our void work | 3 arXiv + 1 WebSearch |
| Galaxy spin-direction asymmetry (Shamir, JADES 2502.18781, vs the Iye and Patel & Desmond nulls) | disagreement | Crowded back-and-forth on many data sets (SDSS, HSC 2410.15269/2411.08723, JWST 2403.17271); its dipole form is on our skip list; the answer turns on classifier choices, not on data | 2 WebSearch |
| Time dilation in quasar variability (Hawkins null vs Lewis & Brewer 2306.04053) | disagreement | Answered (2306.04053, 2501.04171, Pan-STARRS 2106.05399) | 1 arXiv |
| FRB width time dilation vs z | new | Only a single-repeater correction exists (2402.02360); with ~1 dex intrinsic width scatter and z < 1 the exponent can't be pinned (σ ~ 2) | 1 arXiv |
| SN Ia light-curve time dilation at z > 1.5 (HST/JWST) | extension | DES already pins b to ±0.005 out to z ~ 1.2 (MNRAS 533, 3365); ~25 HST SNe would give σ_b ~ 0.1, and no live alternative is left to test | 1 WebSearch |

Not searched (gaps, not conclusions): lensed-quasar delays from public ZTF light curves, GW strong lensing, and
kilonova and superluminous-SN candles.
