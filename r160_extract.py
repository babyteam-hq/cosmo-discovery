"""R-160: stream SPICE-RACS DR2 (spice-racs.dr2.fits.gz, 9.3M rows x 1030 B) and keep the polarised rows.

The gzip is never decompressed to disk: the FITS binary table is parsed in chunks from the gzip stream and only the
columns below, for rows with snr_polint >= 6 or goodRM_flag, are written to data/r160/spice_pol.fits.

Usage: python -I r160_extract.py data/r160/spice-racs.dr2.fits.gz data/r160/spice_pol.fits
"""
import gzip
import sys

import numpy as np
from astropy.table import Table

KEEP = [
    "cat_id", "source_id", "ra", "dec", "total_I_flux", "total_I_flux_err", "peak_I_flux", "maj_axis", "min_axis",
    "dc_maj_axis", "s_code", "N_Gaus", "tile_id", "sbid", "beamdist", "separation_tile_centre", "l_tile_centre",
    "m_tile_centre", "rm", "rm_err", "polint", "polint_err", "stokesQ", "stokesU", "stokesQ_err", "stokesU_err",
    "polangle", "polangle_err", "derot_polangle", "derot_polangle_err", "fracpol", "refwave_sq_pol", "reffreq_pol",
    "snr_polint", "leakage", "snr_flag", "leakage_flag", "channel_flag", "stokesI_fit_flag", "complex_sigma_add_flag",
    "complex_M2_CC_flag", "complex_flag", "local_rm_flag", "is_blended_flag", "goodI_flag", "goodRM_flag",
    "spectral_index", "spectral_index_err", "stokesI_model_order", "stokesI", "l", "b", "epoch", "start_time",
    "type", "notes", "Ncomp", "sigma_add", "rmsf_fwhm", "nn_rm_med", "nn_rm_mad_std",
]

FMT = {"D": ">f8", "E": ">f4", "K": ">i8", "J": ">i4", "L": "S1"}


def parse_header(f):
    blocks = []
    while True:
        b = f.read(2880).decode("ascii")
        blocks.append(b)
        if any(b[i:i + 80].startswith("END ") or b[i:i + 80].rstrip() == "END" for i in range(0, 2880, 80)):
            return "".join(blocks)


def cards(h):
    out = {}
    for i in range(0, len(h), 80):
        c = h[i:i + 80]
        if "=" in c[:10]:
            k = c[:8].strip()
            v = c[10:].split("/")[0].strip().strip("'").strip()
            out[k] = v
    return out


def main(src, dst):
    f = gzip.open(src, "rb")
    parse_header(f)  # primary
    c = cards(parse_header(f))
    nrow, rowlen, nf = int(c["NAXIS2"]), int(c["NAXIS1"]), int(c["TFIELDS"])
    names, fmts = [], []
    for i in range(1, nf + 1):
        t = c[f"TFORM{i}"]
        code = t[-1]
        rep = int(t[:-1]) if t[:-1] else 1
        names.append(c[f"TTYPE{i}"])
        fmts.append(f"S{rep}" if code == "A" else FMT[code])
    dt = np.dtype({"names": names, "formats": fmts})
    assert dt.itemsize == rowlen, (dt.itemsize, rowlen)
    chunk = 200_000
    parts, done = [], 0
    while done < nrow:
        n = min(chunk, nrow - done)
        buf = f.read(n * rowlen)
        a = np.frombuffer(buf, dtype=dt, count=n)
        good = a["goodRM_flag"] == b"T"
        sel = (a["snr_polint"] >= 6) | good
        sub = a[sel]
        parts.append(np.array(sub[KEEP]))
        done += n
        if done % 2_000_000 < chunk:
            print(f"{done}/{nrow} rows, kept {sum(len(p) for p in parts)}", flush=True)
    allk = np.concatenate(parts)
    t = Table(allk)
    for k in KEEP:
        col = t[k]
        if col.dtype.kind == "S" and col.dtype.itemsize == 1 and k.endswith("flag"):
            t[k] = col == b"T"
        elif col.dtype.kind == "S":
            t[k] = np.char.strip(np.char.decode(col, "ascii"))
        elif col.dtype.byteorder == ">":
            t[k] = col.astype(col.dtype.newbyteorder("="))
    t.write(dst, overwrite=True)
    print(f"wrote {len(t)} rows to {dst}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
