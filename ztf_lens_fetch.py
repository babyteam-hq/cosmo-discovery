"""ZTF difference-image cutouts for the cluster-scale lensed quasars (R-127, cosmo-deep-7).

For one target: list every public ZTF science image covering it (IRSA IBE search), then fetch for each g and r epoch a
cutout of the pipeline difference image (sci minus the matched reference, `scimrefdiffimg.fits.fz`) and the difference
image's PSF stamp (`diffimgpsf.fits`), plus one reference-image cutout and the reference PSF-fit catalogue per
(field, ccd, quadrant, filter). Everything lands under data/ztf_lens/<target>/ (gitignored).

Usage (from research/cosmo-discovery):
    python ztf_lens_fetch.py <target> [--size 60] [--workers 8] [--filters zg,zr]
Targets are the keys of TARGETS in ztf_lens_targets.py.
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import threading

import requests

from ztf_lens_targets import TARGETS

_local = threading.local()


def session() -> requests.Session:
    """One keep-alive session per thread: the proxy opens a tunnel per connection and tolerates reuse better."""
    if not hasattr(_local, "s"):
        _local.s = requests.Session()
    return _local.s

IBE = "https://irsa.ipac.caltech.edu/ibe"
DATA = Path(__file__).resolve().parent / "data" / "ztf_lens"


def get(url: str, dest: Path, tries: int = 5) -> str:
    """Download url to dest unless present. Returns 'ok', 'cached', '404' or 'fail'."""
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
        time.sleep(2 * (i + 1))
    return "fail"


def epoch_list(target: str, ra: float, dec: float) -> list[dict]:
    out = DATA / target / "epochs.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if not out.exists():
        r = requests.get(f"{IBE}/search/ztf/products/sci", params={"POS": f"{ra},{dec}", "ct": "csv"}, timeout=600)
        r.raise_for_status()
        out.write_text(r.text)
    return list(csv.DictReader(out.open()))


def sci_base(e: dict) -> str:
    ffd = e["filefracday"]
    return (f"{IBE}/data/ztf/products/sci/{ffd[:4]}/{ffd[4:8]}/{ffd[8:]}/"
            f"ztf_{ffd}_{int(e['field']):06d}_{e['filtercode']}_c{int(e['ccdid']):02d}_o_q{e['qid']}_")


def ref_base(e: dict) -> str:
    f = f"{int(e['field']):06d}"
    return (f"{IBE}/data/ztf/products/ref/{f[:3]}/field{f}/{e['filtercode']}/ccd{int(e['ccdid']):02d}/q{e['qid']}/"
            f"ztf_{f}_{e['filtercode']}_c{int(e['ccdid']):02d}_q{e['qid']}_")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("target")
    ap.add_argument("--size", type=float, default=60.0, help="cutout side in arcsec")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--filters", default="zg,zr")
    a = ap.parse_args()
    t = TARGETS[a.target]
    ra, dec = t["center"]
    epochs = [e for e in epoch_list(a.target, ra, dec)
              if e["filtercode"] in a.filters.split(",") and e["imgtypecode"] == "o"]
    cut = DATA / a.target / "cut"
    cut.mkdir(exist_ok=True)
    center = f"center={ra},{dec}&size={a.size}arcsec&gzip=false"
    jobs = []
    for e in epochs:
        stem = f"{e['filefracday']}_{int(e['field']):06d}_{e['filtercode']}_c{int(e['ccdid']):02d}_q{e['qid']}"
        b = sci_base(e)
        jobs.append((b + "scimrefdiffimg.fits.fz?" + center, cut / f"{stem}_diff.fits"))
        jobs.append((b + "diffimgpsf.fits", cut / f"{stem}_psf.fits"))
    refs = {(e["field"], e["filtercode"], e["ccdid"], e["qid"]): e for e in epochs}
    refdir = DATA / a.target / "ref"
    refdir.mkdir(exist_ok=True)
    for k, e in refs.items():
        stem = f"{int(e['field']):06d}_{e['filtercode']}_c{int(e['ccdid']):02d}_q{e['qid']}"
        b = ref_base(e)
        jobs.append((b + "refimg.fits?" + center, refdir / f"{stem}_refimg.fits"))
        jobs.append((b + "refpsfcat.fits", refdir / f"{stem}_refpsfcat.fits"))
    jobs = [j for j in jobs if not j[1].with_suffix(j[1].suffix + ".404").exists()]
    print(f"{a.target}: {len(epochs)} epochs in {a.filters}, {len(refs)} reference sets, {len(jobs)} files", flush=True)
    counts: dict[str, int] = {}
    t0 = time.time()
    with ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(get, u, d): d for u, d in jobs}
        for n, f in enumerate(as_completed(futs), 1):
            s = f.result()
            counts[s] = counts.get(s, 0) + 1
            if n % 200 == 0:
                print(f"  {n}/{len(jobs)} {counts} {time.time() - t0:.0f}s", flush=True)
    print(f"{a.target}: done {counts} in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    sys.exit(main())
