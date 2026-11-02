# Proposed test: tomography of the near-infrared extragalactic background from SPHEREx quick-release images crossed with DESI galaxies

Question: what is the bias-weighted near-infrared background intensity per unit redshift, b dI/dz, in each SPHEREx spectral channel, measured by cross-correlating SPHEREx images with DESI galaxies in redshift slices? The redshift-channel structure should show the H-alpha (and [OIII], Pa-alpha) line emission moving through the channels; the claim is b dI/dz per channel and redshift bin with intervals, and the line amplitude relative to the continuum.

Data: SPHEREx Quick Release (QR2) spectral images and their zodiacal-light model at IRSA (weekly since July 2025, reprocessed October 2025); DESI DR1 spectroscopic galaxies (BGS, LRG, ELG) and the Legacy Survey for masking; Gaia for stars.

Method: fixed before computing: the fields (the deepest SPHEREx coverage overlapping DESI DR1), the channels, the redshift slices, the masks (stars and galaxies brighter than a fixed limit), the angular scales. Statistic: the angular cross-correlation of each channel's masked, zodi-subtracted image with each galaxy slice's overdensity, deprojected against the zodiacal model, converted to b dI/dz with the slice's galaxy bias from its own clustering (the clustering-redshift estimator). Nulls: redshift-shuffled galaxies, cross-correlation with stars, the zodiacal model alone, half-survey difference maps; errors from jackknife over fields and from the half-difference maps.

Expected precision: a few percent per channel and slice at the forecast level; in practice set by map-making systematics.
