# Proposed test: the bias and redshift distribution of LoTSS-DR3 radio sources from CMB lensing and DESI clustering redshifts

Question: what is the linear bias of the LOFAR Two-metre Sky Survey DR3 radio-source population (13.7 million sources, 88 % of the northern sky, released February 2026), measured by its cross-correlation with CMB lensing, and does it agree with the bias implied by its own angular auto-correlation once its redshift distribution is measured from clustering redshifts with DESI DR1 galaxies? The claim: b (or its amplitude relative to a fixed b(z) model) with its interval, per flux cut, and the ratio of lensing to auto-correlation amplitudes.

Data: LoTSS-DR3 source catalogue and its noise/depth maps (SURF DOI 10.25606/SURF.LoTSS-DR3); Planck PR4 lensing convergence with its released simulations; ACT DR6 lensing with its 400 simulations on the overlap; DESI DR1 galaxies and quasars for clustering redshifts.

Method: fixed before computing: flux cuts (1.5 mJy and 2 mJy, a 5 mJy check), the mask (depth above a threshold, |b| > 20 deg), the multipole range. Pseudo-C_ell (NaMaster) of the radio overdensity with kappa; the auto-spectrum with shot noise subtracted; dN/dz from clustering redshifts in DESI tomographic bins; fit b with the dN/dz uncertainty propagated. Covariance from the lensing simulations; nulls from the lensing curl mode, random catalogues following the depth map, ACT-vs-Planck agreement, flux-cut stability.

Expected precision: about 40 sigma on the cross-spectrum; b to about 3 %.
