"""R-158: Euclid Q1 catalogue from IRSA's merged-objects HATS parquet on S3 (MER + morphology + PHZ), reading only the
needed column chunks of each of the 114 files (12 threads). A pyarrow.dataset scan with a filter ran 22 min without
finishing; this takes ~2 min. Galaxies: VIS < 23 and Sersic VIS radius > 0.4"; stars: point_like_prob > 0.95, 19<VIS<23."""
import os, time, numpy as np, fsspec, pyarrow.parquet as pq, pyarrow as pa
from concurrent.futures import ThreadPoolExecutor
from astropy.table import Table
from r158_common import DATA
fs = fsspec.filesystem("s3", anon=True)
root = "nasa-irsa-euclid-q1/contributed/q1/merged_objects/hats/euclid_q1_merged_objects-hats/dataset"
cols = ["tileid", "object_id", "ra", "dec", "mer_flux_detection_total", "mer_point_like_prob", "mer_spurious_prob",
        "morph_sersic_sersic_vis_radius", "morph_sersic_sersic_vis_axis_ratio", "morph_sersic_sersic_vis_index", "phz_phz_median"]
files = [x for x in fs.find(root) if x.endswith(".parquet")]
def one(f):
    for k in range(4):
        try:
            T = pq.ParquetFile(fs.open(f, block_size=2 ** 20, cache_type="none"), pre_buffer=True).read(columns=cols)
            d = {c: T.column(c).to_numpy(zero_copy_only=False) for c in cols}
            fl = np.nan_to_num(d["mer_flux_detection_total"].astype(float)); r = np.nan_to_num(d["morph_sersic_sersic_vis_radius"].astype(float))
            pl = np.nan_to_num(d["mer_point_like_prob"].astype(float)); sp = np.nan_to_num(d["mer_spurious_prob"].astype(float), nan=1)
            g = (fl > 2.29) & (r > 0.4); s = (pl > 0.95) & (fl > 2.29) & (fl < 91.2) & (sp < 0.1)
            return {c: v[g] for c, v in d.items()}, {c: v[s] for c, v in d.items()}
        except Exception as e:
            err = e; time.sleep(2 ** k)
    raise err
t0 = time.time(); G, S = [], []
with ThreadPoolExecutor(12) as ex:
    for j, (g, s) in enumerate(ex.map(one, files)):
        G.append(g); S.append(s)
        if j % 10 == 0: print(j, len(files), round(time.time() - t0), flush=True)
ren = lambda n: n.replace("mer_", "").replace("morph_", "").replace("phz_phz_", "phz_")
for L, name in ((G, "gal"), (S, "star")):
    T = Table({ren(c): np.concatenate([x[c] for x in L]) for c in cols})
    T.write(os.path.join(DATA, f"euclid_{name}.fits"), overwrite=True); print(name, len(T))
print("done", round(time.time() - t0))
