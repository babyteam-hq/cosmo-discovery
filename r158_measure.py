"""R-158 measurement: classify the prepped stamps with both classifiers at the simulation thresholds, then A per
field and redshift bin (binomial + tile bootstrap), dilution-corrected with eps from the chiral injections and the
armless bound, the mirror null, and the per-galaxy table. Usage: r158_measure.py cosmos|euclid_disc|euclid_ell|euclid_stars|jades|cosmos_ell|cosmos_star
Thresholds/eps: the dataset's own simulations (euclid_* use 'euclid_disc' sims, cosmos_* use 'cosmos')."""
import os, sys, json, numpy as np, pandas as pd, torch
from r158_common import DATA, HERE, cosmos_classes
import r158_classify as C
torch.set_num_threads(4)
ds = sys.argv[1]
simds = "cosmos" if ds.startswith("cosmos") else ("jades" if ds.startswith("jades") else "euclid_disc")
th = json.load(open(os.path.join(DATA, "thresholds_sims.json")))
if simds not in th: simds = "cosmos"          # JADES: no star-based PSF; COSMOS (HST, 0.03"/px) thresholds stand in
T = th[simds]
p = dict(np.load(os.path.join(DATA, f"prep_{ds.split('_')[0] if ds.startswith('cosmos') else ds}.npz")))
if ds.startswith("cosmos"):
    cls, z = cosmos_classes(); want = {"cosmos": "disc", "cosmos_ell": "ell", "cosmos_star": "star"}[ds]
    m = cls[p["ids"]] == want; p = {k: v[m] for k, v in p.items()}
    meta = pd.DataFrame({"z": z.lp_zbest.values[p["ids"]], "tile": z.acstile.values[p["ids"]], "ra": z.ra.values[p["ids"]],
                         "dec": z.dec.values[p["ids"]], "field": "COSMOS", "fwhm_pix": 0.095 / 0.03})
elif ds.startswith("euclid"):
    from astropy.table import Table
    tag = ds.split("_")[1]
    t = Table.read(os.path.join(DATA, {"disc": "euclid_targets_disc.fits", "ell": "euclid_targets_ell2000.fits",
                                       "stars": "euclid_targets_stars800.fits"}[tag]))
    t = t.to_pandas().set_index("object_id").loc[p["ids"]]
    meta = pd.DataFrame({"z": t.get("phz_median", pd.Series(np.nan, index=t.index)).values, "tile": t.tileid.values,
                         "ra": t.ra.values, "dec": t.dec.values,
                         "field": t["field"].str.decode("ascii").values if "field" in t and t["field"].dtype == object and isinstance(t["field"].iloc[0], bytes) else t.get("field", pd.Series("EDF", index=t.index)).values,
                         "fwhm_pix": 1.7})
elif ds == "jadesfield":
    from astropy.table import Table
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = pd.read_csv(os.path.join(DATA, "jadesfield_cat.csv")).set_index("id").loc[p["ids"]]
    g = Table.read(os.path.join(DATA, "euclid_gal.fits")); g = g[(g["ra"] > 52.9) & (g["ra"] < 53.4) & (g["dec"] > -28.05) & (g["dec"] < -27.6)]
    i, d2, _ = SkyCoord(c.ra.values, c.dec.values, unit="deg").match_to_catalog_sky(SkyCoord(g["ra"], g["dec"], unit="deg"))
    zz = np.where(d2 < 1 * u.arcsec, np.array(g["phz_median"])[i], np.nan)
    meta = pd.DataFrame({"z": zz, "tile": (c.x.values // 2048) * 100 + (c.y.values // 2048), "ra": c.ra.values, "dec": c.dec.values,
                         "field": "GOODS-S", "fwhm_pix": 0.07 / 0.03})
else:
    s = pd.read_csv(os.path.join(DATA, "jades_shamir_z.csv")).set_index("id").loc[p["ids"]]
    meta = pd.DataFrame({"z": s.z_euclid.values, "tile": 0, "ra": s.ra.values, "dec": s.dec.values, "field": "GOODS-S",
                         "fwhm_pix": 0.07 / 0.03, "shamir": s.shamir.values})
net = C.make_net(); net.load_state_dict(torch.load(os.path.join(DATA, "net.pt"))); net.eval()
X = torch.tensor(p["X"].astype(np.float32) if p["X"].dtype != np.float32 and len(p["X"]) < 50000 else p["X"])[:, None]
with torch.no_grad():
    f = torch.cat([net(X[b:b + 512].float()) for b in range(0, len(X), 512)]).numpy()
    fm = torch.cat([net(torch.flip(X[b:b + 512].float(), (3,))) for b in range(0, len(X), 512)]).numpy()
S = -C.fft_handedness_batch(p["X"])        # sign fixed in the plan (B = -S)
Sm = -C.fft_handedness_batch(p["X"][:, :, ::-1])
# image position angle from second moments (for the PSF-orientation test)
g = np.arange(64) - 31.5; yy, xx = np.meshgrid(g, g, indexing="ij")
qxx, qyy, qxy = [np.empty(len(p["X"])) for _ in range(3)]
for b in range(0, len(p["X"]), 20000):
    w = np.clip(p["X"][b:b + 20000].astype(np.float32), 0, None); w = w / (w.sum((1, 2))[:, None, None] + 1e-12)
    qxx[b:b + 20000], qyy[b:b + 20000], qxy[b:b + 20000] = (w * xx * xx).sum((1, 2)), (w * yy * yy).sum((1, 2)), (w * xx * yy).sum((1, 2))
pa = (0.5 * np.degrees(np.arctan2(2 * qxy, qxx - qyy))) % 180
out = meta.assign(id=p["ids"], f=f, f_mirror=fm, S=S, S_mirror=Sm, r50=p["r50"], snr=p["snr"], pa=pa,
                  r50_fwhm=p["r50"] / meta.fwhm_pix.values,
                  clsA=np.where(np.abs(f) > T["f"]["t"], np.sign(f), 0), clsB=np.where(np.abs(S) > T["S"]["t"], np.sign(S), 0))
out.to_csv(os.path.join(DATA, f"classified_{ds}.csv"), index=False)
print(ds, "N", len(out), "classified A", int((out.clsA != 0).sum()), "B", int((out.clsB != 0).sum()),
      "mirror exact:", float(np.max(np.abs(f + fm))), float(np.max(np.abs(S + Sm))))
