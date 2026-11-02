# Proposed test: synchrotron emission from the cosmic web between DESI LRG pairs, stacked in LOFAR LoTSS-DR2 144 MHz images

Question: is there diffuse radio emission (and so magnetised plasma) in the filaments between galaxy pairs? Earlier stacks claimed a detection of synchrotron from filaments and others disputed it. Stack LoTSS-DR2 images between physically close DESI DR1 luminous-red-galaxy pairs and measure the mean filament surface brightness, or an upper limit, with an interval; convert to a limit on the filament magnetic field under stated assumptions.

Data: LoTSS-DR2 mosaics (and the low-resolution 20" and 60" images), source catalogues for masking; DESI DR1 LRG spectroscopic catalogue and randoms.

Method: fixed before stacking: pair selection (projected 4-10 Mpc, line-of-sight separation < 6 Mpc/h for physical pairs; unphysical pairs at the same projected separation but separated in redshift as the control), the compact-source subtraction or masking (catalogue sources above a stated flux, with a fixed mask radius law), the frame (pair-aligned, rescaled), the statistic (mean surface brightness in the inter-pair rectangle minus the control stack, minus the two galaxies' own profiles). Errors from rotated stacks and pair bootstrap; nulls from unphysical and rotated pairs.

Expected precision: of order 10^5-10^6 pairs in the LoTSS-DR2 footprint; surface-brightness sensitivity a few tenths of a microjansky per beam, enough to test the claimed signal level.
