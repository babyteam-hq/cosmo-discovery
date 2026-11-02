"""Limber theory for quasar overdensity x CMB lensing convergence at 3.5 < z < 5 (R-139, cosmo-deep-15).

C_ell^{kappa g} = b * C^{kappa m}[n(z)] + (5 s - 2) * C^{kappa mu}[n(z)], Planck 2018 LCDM, CAMB halofit matter power.
The clustering kernel is the normalised spectroscopic n(z); the magnification kernel is the lensing efficiency of the same n(z).
Also C^{gg} (b^2, for the Gaussian covariance and the forecast) and C^{kappa kappa}.
"""
import numpy as np

H0, OMBH2, OMCH2, NS, AS = 67.36, 0.02237, 0.1200, 0.9649, 2.1e-9
C_KMS = 299792.458


class Theory:
    def __init__(self, lmax=2000, zmax=6.0):
        import camb
        p = camb.set_params(H0=H0, ombh2=OMBH2, omch2=OMCH2, ns=NS, As=AS, mnu=0.06, NonLinear=camb.model.NonLinear_both,
                            WantTransfer=True, kmax=20, redshifts=[0.0])
        self.res = camb.get_background(p)
        self.pk = camb.get_matter_power_interpolator(p, nonlinear=True, hubble_units=False, k_hunit=False, kmax=20,
                                                     zmax=1100 if False else 12, var1='delta_tot', var2='delta_tot')
        self.om = (OMBH2 + OMCH2 + 0.06 / 93.14) / (H0 / 100) ** 2
        self.chistar = self.res.conformal_time(0) - self.res.tau_maxvis
        self.z = np.linspace(0.005, zmax + 2.0, 1200)
        self.chi = self.res.comoving_radial_distance(self.z)
        self.Hz = self.res.hubble_parameter(self.z) / C_KMS          # 1/Mpc
        self.lmax = lmax
        # CMB lensing kernel W_k(chi) = 1.5 Om H0^2 (1+z) chi (chi*-chi)/chi*
        h0 = H0 / C_KMS
        self.Wk = 1.5 * self.om * h0 ** 2 * (1 + self.z) * self.chi * (self.chistar - self.chi) / self.chistar
        # kappa from z > z grid max is missing; C^{kk} computed separately from CAMB when needed

    def kernels(self, zc, nz):
        """the clustering kernel W_g(chi) = n(z) H(z) (normalised in z) and the magnification kernel per unit (5s-2)."""
        n = np.interp(self.z, zc, nz, left=0, right=0)
        n /= np.trapezoid(n, self.z)
        Wg = n * self.Hz
        h0 = H0 / C_KMS
        # lensing efficiency of the source distribution: g(chi) = int_{chi}^{inf} dz' n(z') (chi'-chi)/chi'
        g = np.zeros_like(self.z)
        for i in range(len(self.z)):
            sel = self.z > self.z[i]
            if sel.any():
                g[i] = np.trapezoid(n[sel] * (self.chi[sel] - self.chi[i]) / self.chi[sel], self.z[sel])
        Wmu = 1.5 * self.om * h0 ** 2 * (1 + self.z) * self.chi * g
        return Wg, Wmu

    def limber(self, W1, W2, ells=None):
        ells = np.arange(self.lmax + 1) if ells is None else ells
        out = np.zeros(len(ells))
        dchi = np.gradient(self.chi)
        for j, l in enumerate(ells):
            if l < 2:
                continue
            k = (l + 0.5) / self.chi
            ok = k < 19.9
            p = self.pk.P(self.z[ok], k[ok], grid=False)
            out[j] = np.sum(W1[ok] * W2[ok] / self.chi[ok] ** 2 * p * dchi[ok])
        return out

    def spectra(self, zc, nz):
        """unit-bias kappa-g, unit-(5s-2) kappa-mu, unit-bias gg, g-mu and mu-mu (for C^{gg} with magnification)."""
        Wg, Wmu = self.kernels(zc, nz)
        return dict(km=self.limber(self.Wk, Wg), kmu=self.limber(self.Wk, Wmu), gg=self.limber(Wg, Wg),
                    gmu=self.limber(Wg, Wmu), mumu=self.limber(Wmu, Wmu), kk=self.limber(self.Wk, self.Wk))


def forecast(sp, b, s, nbar_sr, fsky, nkk, lmin, lmax):
    """Gaussian S/N of the kappa-g cross power (amplitude fixed at the model), lmin..lmax."""
    l = np.arange(lmin, lmax + 1)
    ckg = b * sp['km'][l] + (5 * s - 2) * sp['kmu'][l]
    cgg = b ** 2 * sp['gg'][l] + 2 * b * (5 * s - 2) * sp['gmu'][l] + (5 * s - 2) ** 2 * sp['mumu'][l] + 1.0 / nbar_sr
    var = (ckg ** 2 + (sp['kk'][l] + nkk[l]) * cgg) / ((2 * l + 1) * fsky)
    return np.sqrt(np.sum(ckg ** 2 / var))
