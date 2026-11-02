"""R-156: stream the DESI DR1 LS-DR9 tractorphot VAC (46 GB, 104 files) and keep only the rows we need.

Kept: every BGS_BRIGHT clustering galaxy (by TARGETID), plus a star sample (Tractor TYPE PSF with a Gaia
parallax at >5 sigma or a proper motion at >5 sigma), with the size, PSF, blending and photometry columns.
Output: data/r156/tractor/hpNNN.npz, one per file, so a crashed run resumes where it stopped.

usage: python -u r156_fetch_tractor.py [nproc]
"""
import os, sys, re, time, urllib.request
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data', 'r156')
OUT = os.path.join(D, 'tractor')
URL = 'https://data.desi.lbl.gov/public/dr1/vac/dr1/lsdr9-photometry/iron/v1.1/observed-targets/tractorphot/'
KEEP = ['TARGETID', 'RELEASE', 'TYPE', 'RA', 'DEC', 'MASKBITS', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_IVAR_R',
        'FIBERFLUX_R', 'FIBERTOTFLUX_R', 'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z',
        'NOBS_R', 'FRACFLUX_R', 'FRACMASKED_R', 'FRACIN_R', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z',
        'PSFDEPTH_R', 'GALDEPTH_R', 'NEA_R', 'BLOB_NEA_R', 'RCHISQ_R', 'SERSIC', 'SHAPE_R', 'SHAPE_R_IVAR',
        'SHAPE_E1', 'SHAPE_E2', 'PARALLAX', 'PARALLAX_IVAR', 'PMRA', 'PMRA_IVAR', 'PMDEC', 'PMDEC_IVAR',
        'GAIA_PHOT_G_MEAN_MAG']
CODES = {'L': 'i1', 'B': 'u1', 'I': 'i2', 'J': 'i4', 'K': 'i8', 'E': 'f4', 'D': 'f8'}


def parse_header(raw):
    cards = [raw[i:i + 80].decode('latin1') for i in range(0, len(raw), 80)]
    return cards


def read_exactly(r, n):
    buf = bytearray()
    while len(buf) < n:
        chunk = r.read(n - len(buf))
        if not chunk:
            break
        buf += chunk
    return bytes(buf)


def dtype_from(cards):
    names, forms = {}, {}
    for c in cards:
        m = re.match(r'(TTYPE|TFORM)(\d+)\s*=\s*\'([^\']*)\'', c)
        if m:
            (names if m.group(1) == 'TTYPE' else forms)[int(m.group(2))] = m.group(3).strip()
    fields = []
    for k in sorted(names):
        f = forms[k]
        m = re.match(r'(\d*)([A-Z])', f)
        rep = int(m.group(1) or 1)
        code = m.group(2)
        if code == 'A':
            fields.append((names[k], 'S%d' % rep))
        elif rep == 1:
            fields.append((names[k], '>' + CODES[code]))
        else:
            fields.append((names[k], '>' + CODES[code], (rep,)))
    return np.dtype(fields)


def work(name):
    out = os.path.join(OUT, name.replace('.fits', '.npz'))
    if os.path.exists(out):
        return name, 'cached', 0
    ids = np.load(os.path.join(D, 'bgs_targetids.npy'))
    t0 = time.time()
    for attempt in range(5):
        try:
            r = urllib.request.urlopen(URL + name, timeout=120)
            # primary HDU, then the table header
            hdr = b''
            nhdr = 0
            cards = []
            while True:
                block = read_exactly(r, 2880)
                hdr += block
                if any(block[i:i + 80].startswith(b'END ') for i in range(0, 2880, 80)):
                    nhdr += 1
                    if nhdr == 2:
                        cards = parse_header(hdr[hdr.find(b'XTENSION'):])
                        break
            dt = dtype_from(cards)
            nrow = int([c for c in cards if c.startswith('NAXIS2')][0].split('=')[1].split('/')[0])
            assert dt.itemsize == int([c for c in cards if c.startswith('NAXIS1')][0].split('=')[1].split('/')[0])
            keep = []
            step = 20000
            done = 0
            while done < nrow:
                n = min(step, nrow - done)
                buf = read_exactly(r, n * dt.itemsize)
                a = np.frombuffer(buf, dtype=dt, count=n)
                tid = a['TARGETID']
                isbgs = np.isin(tid, ids)
                star = (a['TYPE'] == b'PSF') & (
                    (a['PARALLAX'] * np.sqrt(a['PARALLAX_IVAR']) > 5) |
                    (np.hypot(a['PMRA'] * np.sqrt(a['PMRA_IVAR']), a['PMDEC'] * np.sqrt(a['PMDEC_IVAR'])) > 5))
                sel = isbgs | star
                if sel.any():
                    s = a[sel]
                    rec = {k: s[k].astype(s[k].dtype.newbyteorder('=')) if s[k].dtype.kind != 'S' else s[k] for k in KEEP}
                    rec['ISBGS'] = isbgs[sel]
                    keep.append(rec)
                done += n
            r.close()
            res = {k: np.concatenate([x[k] for x in keep]) if keep else np.array([]) for k in KEEP + ['ISBGS']}
            np.savez(out + '.tmp.npz', **res)
            os.replace(out + '.tmp.npz', out)
            return name, len(res['TARGETID']), time.time() - t0
        except Exception as e:  # network hiccup: retry with backoff
            time.sleep(2 ** attempt)
            err = e
    return name, 'FAILED %r' % err, time.time() - t0


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    html = urllib.request.urlopen(URL, timeout=60).read().decode()
    names = sorted(set(re.findall(r'(tractorphot-nside4-hp\d+-iron\.fits)"', html)))
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    t0 = time.time()
    with Pool(nproc) as p:
        for name, n, dt in p.imap_unordered(work, names):
            print('%s %s %.0fs (elapsed %.0fs)' % (name, n, dt, time.time() - t0), flush=True)
    print('DONE %.0fs' % (time.time() - t0), flush=True)
