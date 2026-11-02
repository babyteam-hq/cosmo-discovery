"""The target list of R-127 (cosmo-deep-7): cluster- and group-scale lensed quasars with image separations >= 3" in the
ZTF footprint (Dec > -31), their image positions from the discovery papers, the published delays and the lens-model
predictions published before any measurement. Fixed before any light curve was extracted (2026-11-01).

Positions: J1004 Inada et al. 2003 (astro-ph/0312427, Table 1, SDSS astrometry); J1029 Gaia DR3 (A, B and C; Inada et
al. 2006 astro-ph/0611275 for the system); J2222 Dahle et al. 2013 (arXiv:1211.1091, Table 2, SDSS astrometry; images
D, E, F at r = 23-25 are beyond ZTF and are listed only for the record); J1326 Shu et al. 2019 (arXiv:1909.02010; the BCG
at 201.50006, 48.11208; A and B are each about 10.5" from it, located on the ZTF reference image); COOL J0542 Martinez et
al. 2023 (arXiv:2209.03972, Table of lens-model constraints, images 1.1-1.3); COOL J0335 Napier et al. 2023
(arXiv:2305.14317, Table 1); COOL J1153 Solhaug et al. 2026 (arXiv:2605.13514, Table 2, Magellan slit positions).

Delays: "published" are measured delays to validate against; "predicted" are model predictions published before any
measurement of that pair, with the interval the paper gives. Sign convention: dt_XY > 0 means image Y lags image X
(X leads), i.e. t_Y - t_X, as in Munoz et al. 2022 for J1004 (dt_AC = 825 means C leads A by 825 d) and Dahle et al.
2015 for J2222 (tau_AB = 47.7, A leads B; tau_CA = 722, C leads A).
"""

TARGETS = {
    "j1004": {
        "name": "SDSS J1004+4112", "z_s": 1.734, "z_l": 0.68, "center": (151.1446, 41.2112),
        "images": {  # Inada et al. 2003 Table 1 (J2000), r mags 18.70, 19.10, 19.73, 20.51
            "A": (151.144975, 41.210914), "B": (151.145458, 41.211886),
            "C": (151.140929, 41.209672), "D": (151.141900, 41.213597),
        },
        "published": {  # Munoz et al. 2022 (arXiv:2206.08597), joint BIC, days; order of arrival C, B (+782), A (+825), D (+2458)
            ("B", "A"): (43.01, 0.27, "Munoz+2022 dt_AB (B leads A)"),
            ("C", "A"): (825.23, 0.46, "Munoz+2022 dt_AC"),
            ("C", "B"): (782.20, 0.43, "Munoz+2022 dt_BC"),
            ("C", "D"): (2458.47, 1.02, "Munoz+2022 dt_DC"),
        },
        "predicted": {},
    },
    "j1029": {
        "name": "SDSS J1029+2623", "z_s": 2.197, "z_l": 0.584, "center": (157.3085, 26.3915),
        "images": {  # Gaia DR3; C is 1.8" from B and blends with it in ZTF (B+C fitted as one source at B)
            "A": (157.308112, 26.388320), "B": (157.309369, 26.394478),
        },
        "published": {("A", "B"): (744.0, 10.0, "Fohlmeister+2013 dt_AB (A leads B and C by 744 d; 90 % interval)")},
        "predicted": {},
    },
    "j2222": {
        "name": "SDSS J2222+2745", "z_s": 2.805, "z_l": 0.49, "center": (335.5367, 27.7605),
        "images": {  # Dahle et al. 2013 Table 2 (J2000); r = 20.67, 21.03, 21.14
            "A": (335.537667, 27.760528), "B": (335.536667, 27.761111), "C": (335.532917, 27.760500),
        },
        "published": {
            ("A", "B"): (47.7, 6.0, "Dahle+2015 tau_AB (A leads B; 95 % interval)"),
            ("C", "A"): (722.0, 24.0, "Dahle+2015 tau_CA (C leads A; 95 % interval)"),
        },
        "predicted": {  # Sharon et al. 2017 (arXiv:1609.08848) for D, E, F: r = 23-25, beyond ZTF; recorded, not testable
            ("A", "D"): (502, 68, "Sharon+2017 tau_AD"), ("A", "E"): (611, 75, "Sharon+2017 tau_AE"),
            ("A", "F"): (415, 72, "Sharon+2017 tau_AF"),
        },
    },
    "j1326": {
        "name": "SDSS J1326+4806", "z_s": 2.078, "z_l": 0.396, "center": (201.50006, 48.11208),
        "images": {  # Gaia DR3: the two G ~ 20.8 sources 21.06" apart on either side of the BCG (6.1" SE and 15.0" NW of it);
            "A": (201.500737, 48.110416),  # = the BOSS quasar SDSS J132600.18+480637.4 (DR16Q, z 2.081), Shu's image A
            "B": (201.496635, 48.115585),  # Shu's image B (z 2.078 from WHT)
        },
        "published": {}, "predicted": {},
    },
    "j0542": {
        "name": "COOL J0542-2125", "z_s": 1.84, "z_l": 0.61, "center": (85.7395, -21.4290),
        "images": {  # Martinez et al. 2023 Table (images 1.1, 1.2, 1.3); 1.1 and 1.2 (A, B) are 2.4" apart and blend in
            # ZTF, so they are fitted as one source AB at their mean position; C is 1.3
            "AB": (85.738883, -21.426067), "C": (85.74191404, -21.43245443),
        },
        "published": {},
        "predicted": {},  # arXiv:2209.03972 gives no numerical prediction ("could be multiple years"); none to test
    },
    "j0335": {
        "name": "COOL J0335-1927", "z_s": 3.27, "z_l": 0.4178, "center": (53.7671, -19.4649),
        "images": {  # Napier et al. 2023 Table 1
            "A": (53.769891, -19.462959), "B": (53.766807, -19.464290), "C": (53.764431, -19.466858),
        },
        "published": {},
        "predicted": {  # best-fit model, 1-sigma statistical; systematic ranges 240 < dt_AB < 700, -300 < dt_AC < -30
            ("A", "B"): (499, (146, 141), "Napier+2023 dt_AB (B lags A); syst. range 240-700"),
            ("A", "C"): (-127, (17, 83), "Napier+2023 dt_AC (C leads A); syst. range -300 to -30"),
        },
    },
    "j1153": {
        "name": "COOL J1153+0755 (quasar A)", "z_s": 1.524, "z_l": 0.4301, "center": (178.3318, 7.9356),
        "images": {  # Solhaug et al. 2026 Table 2 (J2000)
            "A1": (178.335243, 7.933551), "A2": (178.331806, 7.936514), "A3": (178.328974, 7.936987),
        },
        "published": {},
        "predicted": {  # ground-based model (the paper's first), 1-sigma; space-based update in the same paper: 1146 (+107/-54), 788 (+48/-55)
            ("A1", "A2"): (1253, (50, 42), "Solhaug+2026 dt_A2,A1 (A1 leads)"),
            ("A1", "A3"): (789, (59, 55), "Solhaug+2026 dt_A3,A1 (A1 leads)"),
        },
    },
}
