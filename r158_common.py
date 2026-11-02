"""R-158 shared helpers: IRSA TAP with retries, paths, and Shamir's pre-registered anchors."""
import io, os, time, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "r158")
os.makedirs(DATA, exist_ok=True)

def fetch(url, out=None, tries=6, timeout=600):
    for k in range(tries):
        r = subprocess.run(["curl", "-sS", "-L", "--fail", "--max-time", str(timeout), url] + (["-o", out] if out else []),
                           capture_output=True)
        if r.returncode == 0:
            return out if out else r.stdout
        time.sleep(2 ** k)
    raise RuntimeError(f"fetch failed {url}: {r.stderr.decode()[:300]}")

def tap(adql, fmt="csv", base="https://irsa.ipac.caltech.edu/TAP/sync"):
    import urllib.parse
    url = base + "?" + urllib.parse.urlencode({"QUERY": adql, "FORMAT": fmt})
    return fetch(url)

def tap_table(adql):
    from astropy.table import Table
    return Table.read(io.BytesIO(tap(adql, fmt="votable")), format="votable")

def tap_async(adql, out, base="https://irsa.ipac.caltech.edu/TAP/async", fmt="votable", poll=15, maxwait=3600):
    """Submit an async TAP job, wait, download the result to `out`."""
    import re, urllib.parse
    for k in range(6):
        r = subprocess.run(["curl", "-sS", "-i", "--max-time", "120", "-X", "POST", base,
                            "--data-urlencode", f"QUERY={adql}", "-d", "LANG=ADQL", "-d", f"FORMAT={fmt}", "-d", "PHASE=RUN"],
                           capture_output=True, text=True)
        m = re.search(r"^[Ll]ocation:\s*(\S+)", r.stdout, re.M)
        if m: break
        time.sleep(2 ** k)
    job = m.group(1)
    t0 = time.time()
    while time.time() - t0 < maxwait:
        try:
            ph = fetch(job + "/phase", timeout=60).decode().strip()
        except Exception:
            ph = "?"
        if ph == "COMPLETED": break
        if ph in ("ERROR", "ABORTED"):
            raise RuntimeError(fetch(job + "/error", timeout=60).decode()[:2000])
        time.sleep(poll)
    fetch(job + "/results/result", out=out)
    return out

def read_vot(path, adql_cols):
    """IRSA's VOTable comes back with col_N names; rename from the SELECT list (comma-separated, aliases stripped)."""
    from astropy.table import Table
    t = Table.read(path)
    names = [c.strip().split(".")[-1].split()[-1] for c in adql_cols.split(",")]
    t.rename_columns(t.colnames[:len(names)], names)
    return t

def cosmos_classes():
    """Per-row class of the COSMOS ZEST x COSMOS2020 table: 'disc' (the pre-registered sample, 0.3<z<1.2), 'disc_edge'
    (0.25-0.3 or 1.2-1.3, fetched as margin, not used), 'ell', 'star', ''."""
    from astropy.table import Table
    z = Table.read(os.path.join(DATA, "cosmos_zest_z.fits")).to_pandas()
    zz = z.lp_zbest
    base = (z.type == 2) & (z.r50 > 9.5) & (z.acs_elongation < 2.5) & (z.junkflag == 0)
    cls = np.where(base & (zz > 0.3) & (zz < 1.2), "disc", np.where(base & (zz > 0.25) & (zz < 1.3), "disc_edge",
          np.where((z.type == 1) & (z.r50 > 9.5) & (z.junkflag == 0), "ell",
          np.where((z.acs_class_star > 0.95) & (z.acs_mag_auto > 19) & (z.acs_mag_auto < 23), "star", ""))))
    return cls, z

def prep_subset(ds):
    """prep arrays restricted to the science sample (COSMOS 'disc'; Euclid 'disc' dir; JADES all)."""
    p = dict(np.load(os.path.join(DATA, f"prep_{ds}.npz")))
    if ds == "cosmos":
        cls, _ = cosmos_classes(); m = cls[p["ids"]] == "disc"
        p = {k: v[m] for k, v in p.items()}
    return p
