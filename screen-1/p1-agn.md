# Proposed test: the host halo mass of eROSITA X-ray AGN from stacked ACT DR6 CMB lensing

Question: what is the mean host dark-matter halo mass of X-ray-selected AGN in the eRASS1 (eROSITA DR1) point-source catalogue, measured directly by CMB lensing, and does it agree with the halo mass inferred from the same sample's clustering?

Data: eRASS1 Main point-source catalogue with its Legacy Survey DR10 counterparts and photometric redshifts (Salvato et al. 2024); ACT DR6 CMB lensing convergence map, mask and its 400 released simulations; Planck PR4 lensing as a check.

Method: select AGN with secure counterparts, split into 3 X-ray luminosity bins and 2 redshift bins (0.2-0.8, 0.8-1.6), fixed before any stack. Stack the ACT kappa map at their positions (and compute the kappa x AGN-overdensity cross-power as the like-for-like statistic); fit a halo-model profile (one-halo NFW plus two-halo with the sample's own bias free) to obtain the effective halo mass per bin; errors from the 400 simulations; nulls from random positions and rotated maps. Compare with the clustering-based halo masses of eROSITA AGN (eFEDS and eRASS1 clustering papers) as a ratio with its interval.

Expected precision: about 10^5 AGN over the ACT footprint; kappa x AGN cross-power detection of order 10-20 sigma; effective halo mass to about 0.1-0.15 dex per luminosity bin.
