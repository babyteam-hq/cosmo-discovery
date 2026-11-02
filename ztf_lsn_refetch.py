"""Refetch of the R-135 sample under a throttled IRSA (R-141, cosmo-deep-10b).

The same files as ztf_lsn_fetch.py (30" r-band difference cutouts) and ztf_lsn_psfcal.py (difference-image PSF stamps),
for the targets of a contiguous rank range, with three changes for speed only:
  - the pre-registered epoch cut (INFOBITS == 0, SEEING <= 4", MAGLIM >= 19; ztf_lsn_detect.clean) is applied to the
    IBE search rows before fetching, so the epochs the analysis drops are never downloaded;
  - one queue across all targets in rank order (the finished prefix grows contiguously) and a short jittered retry on
    IRSA's 503s (which come back at once at any concurrency above ~4 here, 2026-11-02) instead of a long back-off;
  - the PSF stamp of every 10th clean r epoch (in time order) is queued with the cutouts; ztf_lsn_psfcal.py --use-cached
    then measures every stamp present.
Writes the same tree as ztf_lsn_fetch.py, and the marker done_zr_clean per finished target.

Usage (from research/cosmo-discovery):
    python ztf_lsn_refetch.py [--start 0] [--stop 185] [--workers 20] [--skip ID,ID]
"""
from __future__ import annotations

import argparse
import random
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

from ztf_lsn_fetch import DATA, HERE, epoch_list, sci_base, session
from ztf_lsn_psfcal import psf_url


ARGS_G = False


def get(url: str, dest: Path, tries: int = 40) -> str:
    if dest.exists() and dest.stat().st_size > 0:
        return "cached"
    for _ in range(tries):
        try:
            r = session().get(url, timeout=90)
        except requests.RequestException:
            time.sleep(random.uniform(1, 4))
            continue
        if r.status_code == 404:
            dest.with_suffix(dest.suffix + ".404").touch()
            return "404"
        if r.status_code == 200 and r.content[:6] == b"SIMPLE":
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(r.content)
            tmp.rename(dest)
            return "ok"
        time.sleep(random.uniform(0.3, 2.0))
    return "fail"


def jobs_for(tid: str, ra: float, dec: float, size: float = 30.0, gwin: tuple[float, float] | None = None
             ) -> list[tuple[str, Path]]:
    ep = [e for e in epoch_list(tid, ra, dec) if e["filtercode"] == "zr" and e["imgtypecode"] == "o"]
    if gwin is not None:  # g band (candidates' colour vetting): clean g epochs inside the JD window, no PSF stamps
        ep = [e for e in epoch_list(tid, ra, dec) if e["filtercode"] == "zg" and e["imgtypecode"] == "o"
              and gwin[0] <= float(e["obsjd"]) <= gwin[1]]
    ep = [e for e in ep if int(e["infobits"]) == 0 and float(e["seeing"]) <= 4.0 and float(e["maglimit"]) >= 19.0]
    ep.sort(key=lambda e: float(e["obsjd"]))
    cut, pdir = DATA / tid / "cut", DATA / tid / "psf"
    cut.mkdir(parents=True, exist_ok=True)
    pdir.mkdir(exist_ok=True)
    center = f"center={ra},{dec}&size={size}arcsec&gzip=false"
    jobs = []
    for i, e in enumerate(ep):
        stem = f"{e['filefracday']}_{int(e['field']):06d}_{e['filtercode']}_c{int(e['ccdid']):02d}_q{e['qid']}"
        jobs.append((sci_base(e) + "scimrefdiffimg.fits.fz?" + center, cut / f"{stem}_diff.fits"))
        if i % 10 == 0 and gwin is None:
            jobs.append((psf_url(stem), pdir / f"{stem}_psf.fits"))
    return [(u, d) for u, d in jobs if not d.with_suffix(d.suffix + ".404").exists()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int, default=185)
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--skip", default="")
    ap.add_argument("--only", default="")
    ap.add_argument("--gwin", default="", help="JD0:JD1: fetch the clean g epochs in this window instead (with --only)")
    a = ap.parse_args()
    t = pd.read_csv(HERE / "lsn_targets.csv")
    t = t[t.id.isin(a.only.split(","))] if a.only else t.iloc[a.start:a.stop]
    t = t[~t.id.isin(a.skip.split(","))]
    gwin = tuple(float(x) for x in a.gwin.split(":")) if a.gwin else None
    global ARGS_G
    ARGS_G = gwin is not None
    rows = list(t.itertuples())
    t0 = time.time()
    lock = threading.Lock()
    stats = {"ok": 0, "cached": 0, "404": 0, "fail": 0}
    with ThreadPoolExecutor(4) as lister, ThreadPoolExecutor(a.workers) as ex:
        listed = [lister.submit(jobs_for, r.id, r.ra, r.dec, 30.0, gwin) for r in rows]
        pending = []
        for r, lf in zip(rows, listed):
            jobs = lf.result()
            futs = [ex.submit(get, u, d) for u, d in jobs]
            pending.append((r, futs))
            # report finished targets in rank order without blocking the submission of later ones for long
            while pending and all(f.done() for f in pending[0][1]):
                report(*pending.pop(0), stats, lock, t0)
            while len(pending) > 6:  # keep a bounded look-ahead so the prefix stays contiguous
                pr, pf = pending.pop(0)
                for f in pf:
                    f.result()
                report(pr, pf, stats, lock, t0)
        for pr, pf in pending:
            for f in pf:
                f.result()
            report(pr, pf, stats, lock, t0)


def report(r, futs, stats, lock, t0) -> None:
    c: dict[str, int] = {}
    for f in futs:
        s = f.result()
        c[s] = c.get(s, 0) + 1
    with lock:
        for k, v in c.items():
            stats[k] = stats.get(k, 0) + v
        if c.get("fail", 0) == 0 and not ARGS_G:
            (DATA / r.id / "done_zr_clean").touch()
        el = time.time() - t0
        print(f"[{r.rank:3d} {r.id} {r._3}] {len(futs)} files {c} | new {stats['ok']} in {el:.0f}s "
              f"({stats['ok'] / max(el, 1):.1f}/s)", flush=True)


if __name__ == "__main__":
    sys.exit(main())
