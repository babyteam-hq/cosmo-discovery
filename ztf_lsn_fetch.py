"""ZTF r-band difference-image cutouts at the known-lens positions (R-130, cosmo-deep-10).

For each target of lsn_targets.csv in rank order: list every public ZTF science image covering it (IRSA IBE search),
then fetch a 30" cutout of the pipeline difference image (`scimrefdiffimg.fits.fz`, sci minus reference) for every
r-band epoch. g-band is fetched only for candidates (--filters zg,zr --only <id>). Everything lands under
data/ztf_lsn/<id>/ (gitignored): epochs.csv and cut/<stem>_diff.fits; a `done` marker when the target's files are in.
The IRSA access code is R-127's (ztf_lens_fetch.py), with keep-alive sessions per thread; IRSA answers 503 above
~30 connections, so 24 workers share one queue across targets.

Usage (from research/cosmo-discovery):
    python ztf_lsn_fetch.py [--start 0] [--stop 440] [--workers 24] [--size 30] [--filters zr] [--only ID,ID]
"""
from __future__ import annotations

import argparse
import csv
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "ztf_lsn"
IBE = "https://irsa.ipac.caltech.edu/ibe"
_local = threading.local()


def session() -> requests.Session:
    if not hasattr(_local, "s"):
        _local.s = requests.Session()
    return _local.s


def get(url: str, dest: Path, tries: int = 5) -> str:
    if dest.exists() and dest.stat().st_size > 0:
        return "cached"
    for i in range(tries):
        try:
            r = session().get(url, timeout=120)
        except requests.RequestException:
            time.sleep(2 * (i + 1))
            continue
        if r.status_code == 404:
            dest.with_suffix(dest.suffix + ".404").touch()
            return "404"
        if r.status_code == 200 and r.content[:6] == b"SIMPLE":
            dest.write_bytes(r.content)
            return "ok"
        time.sleep(2 * (i + 1) + (3 if r.status_code == 503 else 0))
    return "fail"


def epoch_list(tid: str, ra: float, dec: float) -> list[dict]:
    out = DATA / tid / "epochs.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if not out.exists() or out.stat().st_size == 0:
        for i in range(4):
            try:
                r = session().get(f"{IBE}/search/ztf/products/sci", params={"POS": f"{ra},{dec}", "ct": "csv"}, timeout=600)
                if r.status_code == 200 and r.text.startswith("in_row_id"):
                    out.write_text(r.text)
                    break
                # anything else (IRSA's VOTable error page, e.g. a transient name-resolution fault) is retried
            except requests.RequestException:
                pass
            time.sleep(5 * (i + 1))
        else:
            return []
    rows = list(csv.DictReader(out.open()))
    return rows if rows and "filtercode" in rows[0] else []


def sci_base(e: dict) -> str:
    ffd = e["filefracday"]
    return (f"{IBE}/data/ztf/products/sci/{ffd[:4]}/{ffd[4:8]}/{ffd[8:]}/"
            f"ztf_{ffd}_{int(e['field']):06d}_{e['filtercode']}_c{int(e['ccdid']):02d}_o_q{e['qid']}_")


def jobs_for(tid: str, ra: float, dec: float, filters: list[str], size: float) -> list[tuple[str, Path]]:
    epochs = [e for e in epoch_list(tid, ra, dec) if e["filtercode"] in filters and e["imgtypecode"] == "o"]
    cut = DATA / tid / "cut"
    cut.mkdir(exist_ok=True)
    center = f"center={ra},{dec}&size={size}arcsec&gzip=false"
    jobs = []
    for e in epochs:
        stem = f"{e['filefracday']}_{int(e['field']):06d}_{e['filtercode']}_c{int(e['ccdid']):02d}_q{e['qid']}"
        dest = cut / f"{stem}_diff.fits"
        if dest.with_suffix(dest.suffix + ".404").exists() or (dest.exists() and dest.stat().st_size > 0):
            continue
        jobs.append((sci_base(e) + "scimrefdiffimg.fits.fz?" + center, dest))
    return jobs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int, default=10**6)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--size", type=float, default=30.0)
    ap.add_argument("--filters", default="zr")
    ap.add_argument("--only", default="")
    ap.add_argument("--targets", default="lsn_targets.csv", help="target list (lsn_targets.csv, or val_targets.csv for the validation SNe)")
    a = ap.parse_args()
    t = pd.read_csv(HERE / a.targets)  # itertuples names the "class" column _3
    t = t.iloc[a.start:a.stop]
    if a.only:
        t = t[t.id.isin(a.only.split(","))]
    filters = a.filters.split(",")
    t0 = time.time()
    total = 0
    # the epoch lists three targets ahead, so the cutout queue never waits on a search
    with ThreadPoolExecutor(a.workers) as ex, ThreadPoolExecutor(3) as lister:
        rows = list(t.itertuples())
        listed = {r.id: lister.submit(jobs_for, r.id, r.ra, r.dec, filters, a.size) for r in rows[:3]}
        for i, r in enumerate(rows):
            for nxt in rows[i + 1:i + 4]:
                if nxt.id not in listed:
                    listed[nxt.id] = lister.submit(jobs_for, nxt.id, nxt.ra, nxt.dec, filters, a.size)
            jobs = listed.pop(r.id).result()
            counts: dict[str, int] = {}
            if jobs:
                futs = [ex.submit(get, u, d) for u, d in jobs]
                for f in as_completed(futs):
                    s = f.result()
                    counts[s] = counts.get(s, 0) + 1
            total += len(jobs)
            marker = DATA / r.id / ("done_" + "_".join(filters))
            if counts.get("fail", 0) == 0:
                marker.touch()
            el = time.time() - t0
            print(f"[{r.rank:3d} {r.id} {r._3}] {len(jobs)} files {counts} | "
                  f"{total} files in {el:.0f}s ({total / max(el, 1):.1f}/s)", flush=True)


if __name__ == "__main__":
    sys.exit(main())
