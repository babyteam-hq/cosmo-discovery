"""R-158: Euclid target lists (cuts from r158_plan.md), seeded random order. Fields by position; LDN1641 dropped."""
import os, numpy as np
from astropy.table import Table, unique
from r158_common import DATA
g = Table.read(os.path.join(DATA, "euclid_gal.fits"))
g = unique(g, keys="object_id")
for c in g.colnames:
    if g[c].dtype.kind in "OUS": g[c] = np.array([float(v) if str(v) not in ("", "--", "None") else np.nan for v in g[c]])
vis = -2.5 * np.log10(np.array(g["flux_detection_total"], float)) + 23.9
ra, dec = np.array(g["ra"]), np.array(g["dec"])
field = np.where(dec > 50, "EDF-N", np.where((dec > -33) & (dec < -23), "EDF-F", np.where(dec < -40, "EDF-S", "other")))
g["field"] = field; g["vis"] = vis
pl = np.nan_to_num(np.array(g["point_like_prob"], float)); sp = np.nan_to_num(np.array(g["spurious_prob"], float))
z = np.array(g["phz_median"], float); re = np.array(g["sersic_sersic_vis_radius"], float)
q = np.array(g["sersic_sersic_vis_axis_ratio"], float); n = np.array(g["sersic_sersic_vis_index"], float)
base = (field != "other") & (pl < 0.5) & (sp < 0.5) & (z > 0.3) & (z < 1.2) & (vis < 23)
disc = base & (re > 0.5) & (re < 3) & (q > 0.4) & (n < 2.5)
ell = base & (re > 0.5) & (re < 3) & (q > 0.6) & (n > 3.5) & (vis < 22)
print({f: int((disc & (field == f)).sum()) for f in ("EDF-N", "EDF-F", "EDF-S")}, "disc", disc.sum(), "ell", ell.sum(),
      "other-field", int((field == "other").sum()))
rng = np.random.default_rng(158)
for tag, m, cap in (("disc", disc, None), ("ell", ell, 4000)):
    t = g[m]; t = t[rng.permutation(len(t))]
    if cap: t = t[:cap]
    t["size_arcsec"] = np.clip(8 * np.array(t["sersic_sersic_vis_radius"]), 6, 20)
    t.write(os.path.join(DATA, f"euclid_targets_{tag}.fits"), overwrite=True)
s = Table.read(os.path.join(DATA, "euclid_star.fits")); s = unique(s, keys="object_id")
s["size_arcsec"] = 6.0; s = s[rng.permutation(len(s))]
s.write(os.path.join(DATA, "euclid_targets_stars.fits"), overwrite=True); print("stars", len(s))
