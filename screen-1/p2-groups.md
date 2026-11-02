# Proposed test: the mass scale of DESI Legacy Survey galaxy groups from stacked ACT DR6 CMB lensing

Question: are the halo masses assigned to galaxy groups in the DESI Legacy Imaging Survey DR9 group catalogue (Yang et al. 2021, abundance-matching masses) correct in scale? The ratio of the CMB-lensing mass to the catalogue mass, per mass and redshift bin, with its interval.

Data: Yang et al. 2021 DESI DR9 group catalogue (halo masses from abundance matching, photo-z and spec-z), ACT DR6 lensing convergence with its 400 released simulations, Planck PR4 lensing as a check.

Method: groups with catalogue mass 10^13.5-10^14.5 Msun, 4 mass bins x 2 redshift bins (0.1-0.4, 0.4-0.8), fixed before the fit. Stack kappa at the group centres with a matched NFW filter at each group's catalogue mass convolved with the beam; fit the mass-scale ratio with the two-halo term from the groups' own bias, and the miscentring distribution as a nuisance with a stated prior. Covariance from the 400 simulations; random-position and rotated-map nulls.

Expected precision: about 10^5 groups above 10^13.5 in the ACT footprint; mass-scale ratio to about 5-8 % per bin, 3 % overall.
