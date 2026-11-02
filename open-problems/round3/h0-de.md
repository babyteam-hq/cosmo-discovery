# Round 3: expansion history, H0 and dark energy, quick-kill survivors

Checked 2026-11-01. Searches: WebSearch (extended), the arXiv API
(export.arxiv.org, which answered this session) and arXiv site search.
Overlap and feasibility numbers come from the public catalogues where we
could query them (NOIRLab Data Lab TAP on `desi_dr1.zpix`; the DESI DR1
public mock directory), otherwise from the papers' published tables.

**One survivor, and it is weak.** H0 and the w0wa signal are mined out on
public data. Almost every measurement that could separate 67 from 73 km/s/Mpc,
or move the low-z SN anchor behind w0wa, is done, belongs to a team's own
programme, needs data that isn't public yet (Euclid DR1 slipped to a staged
release starting November 2026; Gaia DR4 is due December 2026), or carries
too little signal. The survivor passes the overlap and physics checks but
decides no disagreement.

## Killed or dropped

- **DESI DR1 redshifts of SN Ia hosts vs Pantheon+ redshifts** (disagreement:
  Steinhardt et al. 2020 on redshift errors biasing cosmology vs Carr et al.
  2112.01471; WiFeS 2402.13484 found a 4.3e-5 mean offset, possibly a missing
  geocentric correction). Unpublished on DESI (4 queries; MOST Hosts
  2405.03857 checks host identity, not redshift offsets). **We ran the
  overlap:** of 667 Pantheon+ SNe with z_HEL < 0.15, 131 have a ZWARN=0 DESI
  DR1 galaxy within 15 arcsec of the SN (Pantheon+ gives no host coordinates
  for most, so the match is by SN position); 118 agree within 0.002. Mean
  offset (-0.7 +/- 3.3) x 1e-5 (SEM), NMAD 1.4e-4, so ~1.3e-5 with a robust
  estimator. **Dropped on physics:** a 4e-5 offset is 0.003 mag at z = 0.03
  and 0.1 km/s/Mpc in H0. The disputes it touches (the 0.04 mag low-z offset
  behind w0wa, the 6 km/s/Mpc tension) are 10-50x larger, so even a clean
  detection decides nothing. (The thirteen >0.002 mismatches are mostly
  SN-position matches that land on a neighbouring galaxy.)
- **Gaia XP cross-calibration of the low-z SN surveys (CfA, CSP, Foundation)**
  against the 0.04 mag offset: done inside Dovekie (2506.05471 sec. 2.5 and
  3.2.2 integrate Gaia XP spectra through each survey's filters; per-survey
  offsets <= 0.027 mag, Table 3).
- **DESI DR1 Lya full-shape Alcock-Paczynski at z = 2.33**: done (2509.15308,
  1.6% on AP; joint with DR2 BAO in 2602.18761).
- **Fine-structure-constant drift from DESI [OIII] (coupled quintessence)**:
  done on DESI ELGs (2404.03123), and limited by a 0.002-0.003 A wavelength
  distortion that mimics the signal. Equivalence-principle bounds on the
  coupling make the cosmological signal smaller still.
- **PNLF distances to SN Ia hosts from the public MUSE archive**: the
  Jacoby/Ciardullo programme (2024 paper II: H0 = 74.2 +/- 7.2 +/- 3.7; only
  four SNe Ia in three galaxies have usable archival MUSE). Can't separate
  67 from 73.
- **Low-z deceleration (q0) from DESI DR1 FP distances, SN-free** (to test the
  low-z anchor): passive FP zero-point evolution over 0.01 < z < 0.1
  (~0.02-0.03 mag in r) and z-dependent Malmquist bias are as large as the
  ~0.02 mag q0 signal. Dropped on physics.
- **Dark sirens with GWTC-4 and the DESI DR1 spectroscopic catalogue**:
  ~10-20% on H0 at best; can't tell the answers apart.
- **Angular BAO from DESI DR1 BGS at z < 0.4** (extending 2610.10215, which
  starts at z = 0.4): angular-only BAO in thin low-z shells is far weaker than
  DESI's own 3D BGS BAO, so it adds nothing.
- **Anything on Euclid DR1**: not public yet (staged release, first part
  November 2026).

---

## 1. BAO from cosmic voids in DESI DR1 (void-galaxy and void-void correlations)

- **Question:** What are D_V/r_d (and D_M/r_d, D_H/r_d where the S/N allows)
  from the void-galaxy and void-void correlation functions of DESI DR1, per
  tracer and redshift bin, alone and combined with the galaxy BAO? In
  particular, the first galaxy-traced void BAO at 1.1 < z < 1.6 (ELG2).
- **Shape:** regime extension (sample and redshift). Void BAO exists on BOSS
  (Kitaura et al. 2016, 1511.04405, first detection; Zhao et al. 2020,
  1802.03990: +18% at 0.2 < z < 0.5) and eBOSS (Zhao et al. 2022, 2110.03824,
  LRG+ELG, 5-15% tighter D_V; eBOSS quasar voids at z = 1.48, 2208.06238, no
  gain on data). DESI DR1 has void catalogues (2610.18348: size function,
  profiles, bias; DESIVAST for BGS) but no void BAO.
- **Public data:**
  - DESI DR1 LSS clustering catalogues and randoms:
    https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/
  - DESI DR1 EZmocks, 1000 realisations, dark (LRG, ELG, QSO) and bright
    (BGS), with fibre-assignment emulation (`*_ffa_*_clustering`), checked:
    https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/EZmock/dark/v1/ ;
    AbacusSummit DR1 mocks for the template:
    https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/AbacusSummit/
  - Void finder: DIVE Delaunay voids as in Zhao et al.
    (https://github.com/cheng-zhao/DIVE) or Revolver
    (https://github.com/seshnadathur/Revolver); reconstruction with pyrecon
    (https://github.com/cosmodesi/pyrecon).
- **Method:**
  1. Reconstruct each tracer and bin as DESI DR1 did, then find Delaunay
     voids on the reconstructed galaxies and on each EZmock the same way.
  2. Measure the void-galaxy and void-void correlation functions (radius cut
     chosen on mocks to maximise the BAO S/N) and the galaxy auto.
  3. Fit alpha_iso (alpha_AP where the S/N allows) with a template from
     AbacusSummit voids and covariance from the 1000 EZmocks, blinded until
     the mock tests pass.
  4. Report void-only alpha and the joint galaxy+void D_V/r_d per bin, with
     the EZmock-derived correlation, and the gain over DESI DR1's galaxy-only
     BAO.
- **Overlap and feasibility:** a single survey, so the overlap is the DR1 LSS
  footprint: BGS 7,473 deg^2 (300,017 redshifts), LRG 5,840 deg^2 (2.14 M),
  ELG 5,914 deg^2 (2.43 M), QSO 7,249 deg^2 (856,652) (2404.03000 Table 2).
  Physics: nothing in the estimator removes the BAO mode. The voids trace the
  same density field, so their information partly repeats the reconstructed
  galaxy BAO; that redundancy is the main limit, and the EZmock
  cross-covariance captures it. Fibre assignment changes the void population,
  but the ffa EZmocks carry the same incompleteness, so it enters the
  covariance and template rather than as an unmodelled bias. On eBOSS QSOs
  the gain vanished on data (2208.06238), so the ELG2 bin may show none.
- **Expected precision:** void-only alpha_iso ~1.5-3% per bin (scaling the
  BOSS and eBOSS void results by DR1's ~2.5-3x effective volume); joint
  galaxy+void D_V/r_d 5-15% tighter than DESI DR1's galaxy-only values (e.g.
  ~0.86% -> ~0.75-0.8% at z = 0.93). That moves the w0wa significance by well
  under 0.5 sigma. **It decides no disagreement:** it gives a clean, partly
  independent BAO cross-check on a new sample, and nothing more. It is also
  likely on DESI's own list (Zhao and Kitaura are in DESI), so it only pays
  if done quickly.
- **Quick-kill queries (5):** arXiv API all:void AND all:BAO AND all:DESI (4
  hits: 2607.13362, 2606.02712 and 2604.05466 are not void BAO; 2208.06238 is
  eBOSS); arXiv search "void BAO DESI", newest first (the same 4); WebSearch
  extended "void BAO DESI DR1 voids baryon acoustic oscillation measurement"
  (angular BAO 2610.10215 and methodology 2107.02950 only); WebSearch
  extended "void BAO galaxy-void cross-correlation DESI emission line
  galaxies ... multi-tracer voids reconstruction" (2610.18348 void structure,
  2605.12482 f(R) multipoles, Roman/Uchuu simulation work; no BAO); arXiv
  API id_list 2110.03824, 2208.06238, 1802.03990 for the prior precisions.
  No DESI void BAO measurement was found.

## Quick-kill queries for the dropped candidates

- Host redshifts: WebSearch extended "DESI redshifts of type Ia supernova host
  galaxies compared Pantheon+ redshift systematic offset" (WiFeS 2402.13484
  only); arXiv API DESI AND supernova AND host AND redshifts (15 hits, none
  compares redshifts); 2405.03857 abstract and full text (MOST Hosts: host
  identity vs Pantheon+, no z offset); arXiv API redshift AND errors AND
  supernova AND cosmological (2402.13484 is the only independent-spectra
  check). Overlap computed by box queries on `desi_dr1.zpix` through
  https://datalab.noirlab.edu/tap/sync (scratch script, not committed).
- Gaia XP: WebSearch extended "Gaia XP synthetic photometry cross-calibration
  low-redshift supernova surveys CSP CfA Foundation zero point offsets";
  2506.05471 abstract and full text.
- Lya AP: WebSearch extended "DESI DR1 Lyman-alpha forest full-shape
  Alcock-Paczynski measurement Cuceu" (2509.15308).
- Alpha drift: WebSearch extended "fine-structure constant variation DESI
  quasar [OIII] doublet DR1 measurement" (2404.03123).
- PNLF: WebSearch extended "planetary nebula luminosity function Hubble
  constant MUSE type Ia supernova host galaxies 2025".
- Ladder disputes: WebSearch extended "2026 arXiv disagreement CCHP SH0ES JWST
  TRGB JAGB distances same galaxies discrepancy cause unresolved" (2601.00650,
  2606.26831: the differences trace to sample selection and zero points that
  the teams already contest on their own data).
- Euclid: WebSearch "Euclid DR1 public data release date 2026".
