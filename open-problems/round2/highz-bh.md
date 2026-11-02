# Round 2: high-z galaxies, reionisation and black holes, quick-kill survivors

Checked 2026-11-01. The arXiv API (export.arxiv.org) returned HTTP 429 on every
call this session, so the arXiv checks used arXiv's own search
(arxiv.org/search, all fields, newest first), which covers the same index.
Query counts include those plus WebSearch (extended).

**Result: no new survivors.** Two candidates passed the quick kill, but parallel
round-2 agents had already written both up, so they're cross-referenced here
rather than duplicated:

- **CMB lensing x z~5 HSC r-dropouts on ACT DR6** (regime extension of
  Miyatake et al. 2022, arXiv:2103.15862, z~4 with Planck PR3 only): written up
  as `cmb-lss.md` #1. Our own 10 queries agree: nothing at z>=4.5. The
  only later LBG x lensing measurement is UNIONS at z~2.5 x Planck
  (2607.09943). The rest are forecasts (1904.13378, 2511.22243, 2401.13166).
- **tSZ around quasar hosts, feedback or gravity** (the disagreement is Ruan et al.
  1502.01723 and Crichton et al. 1510.05656, which call it feedback, against
  Cen & Safarzadeh 1507.07934 and Soergel et al. 1703.00920, which put it down to the host
  halo. Hall et al. 1907.11731 got 3.8 sigma at z>1.91 and couldn't decide). It's written up as
  `cmb-lss.md` #2 (DESI DR1 x ACT DR6, with the halo mass fixed by DESI DR2
  clustering, 2606.27115) and `baryons-dm.md` #2 (the z = 3.5-5 extension). Our 10
  queries found nothing newer than 2310.09445 (an ALMA stack of 10 quasars at z 2.2-2.7).

## Killed in the quick kill (question already answered on some data)

| candidate | shape | queries | killed by |
|---|---|---|---|
| Halo mass of radio-loud vs radio-quiet quasars beyond z~2.3 (DESI/eBOSS x LoTSS DR2) | extension of Retana-Montenegro & Rottgering 1611.08630 (z 0.3-2.3) and Shen et al. 0810.4144 | 8 (6 arXiv, 1 WebSearch, 1 abstract) | Yue et al. 2601.20576: 225k eBOSS quasars in LoTSS DR2, clustering vs jet fraction; jet-dominated hosts sit in 10^13-10^14 h^-1 Msun halos |
| Obscured vs unobscured quasar halo mass (third method) | disagreement: DiPompeo/Petter (1411.0527, 1705.05306, 2302.00690) vs X-ray results (1409.7693) | 5 | Petter et al. 2302.00690 already found it evolves with redshift, and HSC Type I/II HODs exist (2510.11780); no clean third method on public data |
| Quasar clustering vs BH mass at fixed luminosity, z 2-3.5 | disagreement (Shen 0810.4144, Krolewski 1501.03898) | 4 | DESI DR1 PAC (2607.21749), DESI DR2 clustering (2606.27115), HOD vs M_BH (2508.08851); also `11-black-holes.md` #10 |
| LRD clustering / halo mass from public catalogues | disagreement (2411.11534 vs 2506.04004) | 3 | Crowded: 2606.09721, 2605.27903, 2610.23806; already `09-high-z.md` / `11-black-holes.md` #7 |
| UV-background Gamma_HI(z) from the statistical line-of-sight proximity effect in DESI | extension of Dall'Aglio 0906.1484 (z 2.3-4.6) and Calverley 1011.5850 | 8 | DESI Y1 proximity zones measured (2512.18281); method for DESI published (2504.03848) |
| Mean free path of ionising photons at z 3-5 from DESI stacks | extension of Worseck 2014 / Zhu 2023 | 1 (decisive) | DESI Y1, 3.2<=z<=4.6 (2411.15838) |
| L_X-L_UV relation evolution from eRASS1 stacks of quasars | disagreement (Lusso & Risaliti vs Khadka & Ratra) | 3 | 2508.00732, 2512.05807 (eROSITA x SDSS) |
| Quasar bias at z~4 from quasar x HSC g-dropout cross-correlation | disagreement (Shen 2007 high vs He et al. 1704.08461 / Uchiyama 1704.06050 lower; Garcia-Vergara 1701.01114 high) | 4 | Same question as the skip-list item "quasar bias at z>3.5 from CMB lensing" (`round1/g-qso-z4-lensing.md`); not re-proposed. If that lensing attempt lacks S/N, DESI DR1 z 3.5-4.5 quasars x GOLDRUSH g-dropouts is the natural second method |
| Reionisation clustering boost of LAEs at z=6.6 via LAE x LBG cross | extension of Ouchi et al. 2018 LAE auto-clustering | 3 | LAE clustering at 5.7/6.6/7.0 already constrains it (SILVERRUSH, Umeda et al.); HSC NB fields (~25 deg^2) too small for the cross to beat the auto by much |

Data for the two cross-referenced survivors (GOLDRUSH, ACT DR6 lensing and y
maps, DESI DR1) are listed in `cmb-lss.md` and `baryons-dm.md`.
