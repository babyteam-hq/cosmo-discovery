"""cosmo-deep-28 (R-160): common pieces for the quasar jet alignment region (QJAR) re-test.

Frozen template from Mandarakas et al. 2021 (arXiv:2106.06546, A&A 653, A123), Table 1. Nothing here is tuned on
the new sample.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data', 'jet')

# id, RA, Dec, radius (deg), zmin, zmax, published S.L. (sigma), published mean PA (deg, parallel-transported to the centre)
QJARS = [
    ('QJAR1', 36.18, 25.95, 14.82, 0.02, 1.53, 5.54, -16.9),
    ('QJAR2', 168.75, 47.54, 18.11, 0.36, 1.89, 6.74, -17.7),
    ('QJAR3', 180.15, 23.56, 16.07, 0.29, 0.65, 5.72, -25.6),
    ('QJAR4', 201.83, 50.97, 18.27, 0.07, 0.35, 5.38, 37.6),
]
CLAIM_CUTOFF = '2020_06_10'  # Astrogeo database version used by the claim


def unitvec(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)


def angsep(ra1, dec1, ra2, dec2):
    a, b = unitvec(ra1, dec1), unitvec(ra2, dec2)
    return np.degrees(np.arccos(np.clip((a * b).sum(-1), -1, 1)))


def basis(ra, dec):
    """Local north and east unit vectors."""
    ra, dec = np.radians(ra), np.radians(dec)
    north = np.stack([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)], -1)
    east = np.stack([-np.sin(ra), np.cos(ra), np.zeros_like(ra)], -1)
    return north, east


def transport(pa, ra, dec, ra0, dec0):
    """Parallel-transport position angles (deg, N through E) at (ra, dec) along great circles to (ra0, dec0).

    Returns the PA measured at the centre. Implemented as the rotation about the axis A x B that takes A to B
    (equivalent to Eq. 1 of 2106.06546, xi' = xi + chi_b - chi_a)."""
    pa = np.radians(np.asarray(pa, float))
    ra = np.asarray(ra, float); dec = np.asarray(dec, float)
    A = unitvec(ra, dec); B = unitvec(np.full_like(ra, ra0), np.full_like(dec, dec0))
    nA, eA = basis(ra, dec)
    v = np.cos(pa)[:, None] * nA + np.sin(pa)[:, None] * eA
    axis = np.cross(A, B); s = np.linalg.norm(axis, axis=-1); c = (A * B).sum(-1)
    k = np.where(s[:, None] > 1e-12, axis / np.maximum(s, 1e-12)[:, None], 0.0)
    # Rodrigues rotation by angle theta with sin=s, cos=c
    vr = v * c[:, None] + np.cross(k, v) * s[:, None] + k * (k * v).sum(-1)[:, None] * (1 - c)[:, None]
    nB, eB = basis(np.array([ra0]), np.array([dec0]))
    return np.degrees(np.arctan2((vr * eB).sum(-1), (vr * nB).sum(-1)))


def axial_mean(pa_deg, w=None):
    """Mean axial direction (deg, in (-90, 90]) and resultant length D of doubled angles."""
    t = np.radians(2 * np.asarray(pa_deg, float))
    w = np.ones_like(t) if w is None else np.asarray(w, float)
    C = (w * np.cos(t)).sum() / w.sum(); S = (w * np.sin(t)).sum() / w.sum()
    return np.degrees(0.5 * np.arctan2(S, C)), np.hypot(C, S)


def axdiff(a, b):
    """Axial difference a - b wrapped to [-90, 90)."""
    return (np.asarray(a, float) - np.asarray(b, float) + 90.0) % 180.0 - 90.0


def in_region(q, ra, dec, z):
    _, ra0, dec0, rad, zmin, zmax = q[:6]
    return (angsep(ra, dec, ra0, dec0) <= rad) & (z >= zmin) & (z <= zmax)
