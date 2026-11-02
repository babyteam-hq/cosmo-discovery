"""Pantheon+ Hubble residuals vs dipole/quadrupole in several frames, full covariance GLS.
usage: python frames.py data/Pantheon+SH0ES.dat data/Pantheon+SH0ES_STAT+SYS.cov"""
import sys, numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.stats import chi2 as chi2d
C=299792.458; OM=0.334
def dl(z):
    e=lambda x:1/np.sqrt(OM*(1+x)**3+1-OM); return (1+z)*C/70*np.array([quad(e,0,zi)[0] for zi in z])
def unit(ra,dec):
    ra,dec=np.radians(ra),np.radians(dec); return np.stack([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)],-1)
R_gal=np.array([[-0.0548755604,-0.8734370902,-0.4838350155],[0.4941094279,-0.4448296300,0.7469822445],[-0.8676661490,-0.1980763734,0.4559837762]])
eps=np.radians(23.4392911); R_ecl=np.array([[1,0,0],[0,np.cos(eps),np.sin(eps)],[0,-np.sin(eps),np.cos(eps)]])
cmb_dip=unit(167.942,-6.944)  # RA,Dec of CMB dipole apex (l=264.02,b=48.25)
def main(dat,covf):
    d=pd.read_csv(dat,sep=r'\s+')
    with open(covf) as f: n=int(f.readline()); cov=np.loadtxt(f).reshape(n,n)
    for zcol in ['zHD','zCMB','zHEL']:
        m=(d.IS_CALIBRATOR==0)&(d[zcol]>0.01); idx=np.where(m)[0]
        z=d[zcol].values[idx]; y=d.m_b_corr.values[idx]-(5*np.log10(dl(z))+25)
        n_eq=unit(d.RA.values[idx],d.DEC.values[idx])
        L=np.linalg.cholesky(cov[np.ix_(idx,idx)]); w=lambda X:np.linalg.solve(L,X)
        def fit(cols):
            X=np.column_stack([np.ones(len(z))]+cols); Xw,yw=w(X),w(y)
            b,*_=np.linalg.lstsq(Xw,yw,rcond=None); r=yw-Xw@b
            return b,np.sqrt(np.diag(np.linalg.inv(Xw.T@Xw))),r@r
        _,_,c0=fit([])
        print(f"\n=== {zcol} N={len(z)} chi2_0={c0:.1f}")
        for fname,R in [('equatorial',np.eye(3)),('galactic',R_gal),('ecliptic',R_ecl@np.eye(3))]:
            nn=n_eq@R.T
            # latitude-only (axisymmetric) terms: sin(lat) [dipole along pole], P2(sin lat) [quadrupole along pole]
            s=nn[:,2]
            for tname,cols,dof in [('pole dipole sin(lat)',[s],1),('pole quad P2(sin lat)',[1.5*s*s-0.5],1),
                                   ('full dipole',[nn[:,0],nn[:,1],nn[:,2]],3)]:
                b,e,c=fit(cols); dc=c0-c
                print(f"  {fname:10s} {tname:22s} dchi2={dc:6.2f} dof={dof} p={chi2d.sf(dc,dof):.3f}  coef={np.round(b[1:],4)} ±{np.round(e[1:],4)}")
        s=n_eq@cmb_dip
        for tname,cols in [('CMB-dipole axis dipole',[s]),('CMB-dipole axis quad',[1.5*s*s-0.5])]:
            b,e,c=fit(cols); dc=c0-c; print(f"  {'cmb-axis':10s} {tname:22s} dchi2={dc:6.2f} dof=1 p={chi2d.sf(dc,1):.3f}  coef={np.round(b[1:],4)} ±{np.round(e[1:],4)}")
main(*sys.argv[1:])
