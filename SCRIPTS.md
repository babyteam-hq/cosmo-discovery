# Cosmo discovery

Hunting a novel cosmological result on public data. Every question gets a literature check before any compute
(what's been done is listed in `LOG.md`), and every test and every finding goes through the blind
review in `blind-review.md`: a context-free Opus referee, its own literature search, a fixed verdict line, recorded verbatim.

Data (gitignored, `data/`): Pantheon+SH0ES distances and STAT+SYS covariance from
github.com/PantheonPlusSH0ES/DataRelease, `4_DISTANCES_AND_COVAR/`.

- `frames.py`: Hubble residuals (flat LCDM, Om=0.334, calibrators out, z>0.01) vs dipole and axisymmetric
  quadrupole in equatorial, galactic, ecliptic and CMB-dipole frames, full-covariance GLS, for zHD/zCMB/zHEL.
  Result 2026-10-31: zHEL shows the CMB dipole (sanity, 6 sigma); zCMB the known residual bulk flow; zHD nothing
  above 2 sigma except the CMB-axis quadrupole (p=0.03) Macpherson 2023 already reports. Ecliptic: null.
- `fit_dipole.py`: bulk-flow (1/z) vs constant dipole, redshift shells. Not run: already published
  (Tsagas 2024, arXiv:2608.07209, 2605.03004).
- `vpec_scale.py`: scale factor beta on Pantheon+'s 2M++ peculiar-velocity correction (beta=1 means the flow model's
  amplitude is right), overall and by distance shell. Result 2026-10-31: beta=0.87±0.14; no shell past 2 sigma;
  shells vs one beta p=0.20. Null: the corrections' amplitude is right.
- `cf4_vs_2mpp.py`, `cf4_decompose.py`, `cf4_variants.py`, `cf4_2023.py`, `cf4_robust.py`: Pantheon+ residuals vs the
  Cosmicflows-4 reconstructed velocity fields (2023 FITS grids and the 2025 CF4++ npz, from
  projets.ip2i.in2p3.fr/cosmicflows, into `data/`). Result 2026-10-31, see `finding.md`: the fields' bulk flow is
  supported (beta ~1), their non-dipole structure rejected (beta 0.10±0.07 for CF4++, 13 sigma from 1).
- `bulkflow_free.py`: the bulk flow the Pantheon+ SNe measure on their own (free 3-vector, GLS), vs the CF4++ and 2M++
  mean vectors over the same SNe. 312±44 km/s toward (l,b)=(311,18); see the addendum in `finding.md`.

R-115 (cosmo-sweep, 2026-10-31): `sweep.md` (31 literature-checked tests, ranked), `results.md` (a row per test run),
`finding.md` second section (the dust-map result). Scripts: `dust_maps.py` (residuals vs SFD/CSFD/GNILC/MWEBV, scatter
anisotropy), `dust_robust.py`, `dust_resolution.py`, `host_in_sfd.py`, `dust_emissionfree.py` (Edenhofer 2023 stellar
map), `ztf_tests.py` (ZTF SN Ia DR2 replication, dust and CF4++), `cf4_notin.py`, `cf4_des.py` (DES Dovekie low-z),
`cf4_axes.py` (grid convention scan), `bulk_flow.py` (free bulk-flow vectors per sample), `env_density.py`,
`frb_dm_dipole.py` (CHIME Catalog 2). Data (gitignored, `data/`): `dustmaps/` (SFD from github.com/kbarbary/sfddata,
CSFD and Planck via the dustmaps package, Edenhofer 2023 from Zenodo), `ztf/` (ztfcosmo.in2p3.fr/download/data),
`des5yr_repo/` (sparse clone of des-science/DES-SN5YR), `cf4_table2.csv` (VizieR J/ApJ/944/94 via vizier.idia.ac.za),
`chimefrbcat2.csv` (assets.chime-frb.ca).

R-116 (cosmo-sweep-2, 2026-10-31): `dust_settle.py` (settling the dust candidate: Bayestar 2019 and Lenz 2017 emission-free
maps, point and 1°-averaged, with the slope each would show if GNILC were right; the GNILC-for-SFD swap propagated through
the SALT2 response; the Galactic-axis terms before and after; GNILC dust temperature), `dust_salt_response.py` (the SALT2
response to a Milky-Way E(B-V) error, sncosmo). Data (gitignored, `data/`): `dustmaps/bayestar/bayestar2019.h5` and
`dustmaps/lenz2017/ebv_lhd.hpx.fits` from the Harvard Dataverse API (`/api/access/datafile/3424724` and `/3036047`; send a
browser user agent, the dustmaps fetcher's Python one gets a 403), `dustmaps/planck/COM_CompMap_Dust-GNILC-Model-
Temperature_2048_R2.01.fits` from the Planck Legacy Archive, Edenhofer 2023 `mean_and_std_healpix.fits` from Zenodo 8187943
(3.25 GB; the proxy cuts long transfers, resume with `curl -C -`). Verdict in `finding.md`: partial. `sweep-2.md`,
`results.md` (R-116 section) for the second sweep.
Sweep 2 scripts: `trail_ngc1052.py` (NGC 1052 trail kinematics, membership, distances, chance alignment, DESI spectra),
`cf4_dust.py` (CF4 distances vs dust maps and hemisphere zero points), `erass1_nz.py`, `dust_stellar2.py` (Zhou 2024 and
Mudur 2023 stellar maps at the SNe), `dust_mapcorr.py`, `quaia_dust.py`, `quaia_zerr.py`, `desi_nz_dust.py`,
`desi_lrg_dipole.py`, `sn_kappa.py`, `dust_colour.py` and `dust_colour2.py` (the SN colour test). Data (gitignored): `cf4/`
(VizieR J/ApJ/944/94 table2 via vizier.idia.ac.za), `desi_trail/` (DESI DR1 redrock healpix 36785/36787), `desi_lss/` (DR1
LSScats v1.5 LRG and QSO clustering catalogues and one LRG random per cap), `quaia/` (Zenodo 8060755), `erass1/` (VizieR
J/A+A/685/A106 via IDIA), `zhou/` (data.desi.lbl.gov/public/papers/mws/desi_dust/y2/v1/maps), `mudur/` (Zenodo 7411344,
15' version), `dustmaps/csfd/` (Zenodo 8207175), `pr4/` (github.com/carronj/planck_PR4_lensing release "Data").
- `trail_roman.py`: the DF2-DF4 strip in Roman+2021's independent LSBG catalogue (42 objects, transcribed from ar5iv): 6 on
  the strip vs 1.9 expected (p=0.01, radius-preserving null); richest strip over all pair-lines p=0.12 at 5.2'. See results.md.
- `dm_content.py`: M_dyn/M_* inside r_1/2 for DF2, DF4, DF9 at D = 13..21.7 Mpc for every published dispersion; the
  dark-matter deficit holds at every cell (1-3 vs 11-35 for a normal halo). Table 2 of the trail note.
- `erass1_lxt.py`: eRASS1 cluster L_X-T anisotropy toward the Migkas axis and its exposure dependence (data: eRASS1 primary
  cluster catalogue from erosita.mpe.mpg.de into `data/erass1/`). The deficit (4-7 sigma raw, the published ~9 % H0 amplitude)
  vanishes with a log-exposure term or at fixed exposure. Third finding; see finding.md and the site note.

R-124 (cosmo-deep-4, 2026-11-01): `periodic_cw_limits.py`, the NANOGrav 15-year continuous-wave limits against the
periodic-quasar binary candidates as a population (`results-R124.md`). Trajectory check "partly novel": NANOGrav's own
targeted 15-year search on the 111 CRTS candidates (Agarwal et al. 2025, arXiv:2508.16534) covers the brief's `enterprise`
pass; the open part (the population tally with mass and mass-ratio nuisances, the PTF/PS1/DES/ZTF catalogues, the
period-doubling convention) is the cheap pass: 0 of 192 candidates excluded at their claimed parameters (0-1.9 %), the
closest a factor 6 below its limit. Data (gitignored, `data/`): `ng15/` (Zenodo 8067506, the CW analysis release with
`15yr_cw_3d_limits_v4.npz`, the sky-and-frequency 95 % limit map; 8092346 noise spectra; 16051178 the 15-year timing data,
unused), `periodic/` (arXiv e-print sources of the candidate papers; the targeted paper's PDF and its `pdftotext`). Zenodo
answers from the cloud; every VizieR/CDS host is cut by the proxy.

R-125 (cosmo-deep-5, 2026-11-01): `results-R125.md` only. CHIME/FRB Catalog 2 DM excess vs the Planck/ACT Compton-y maps
(the temperature of FRB-traced gas): trajectory check "not novel" (arXiv:2511.02155, 2604.22105, 2608.09014, verified on
arXiv), killed before any download; reachability of the CHIME, PLA and NERSC ACT DR6 sources recorded there.
