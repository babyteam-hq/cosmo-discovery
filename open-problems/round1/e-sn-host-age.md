# Proposed test: the slope of SN Ia Hubble residuals with host stellar age, from DESI DR1 host spectra

Question: do standardised Type Ia supernova luminosities depend on the age of their host's stellar population? A strong progenitor-age bias has been claimed to change the dark-energy inference from supernovae (and to bring it into line with DESI BAO). Measure the Hubble-residual slope with host age, in mag per Gyr with an interval, from ages derived from the hosts' own DESI spectra.

Data: Pantheon+ and ZTF SN Ia DR2 (public light-curve fits and host coordinates); DESI DR1 spectra of the hosts (matched by position; public); stellar-population fits (e.g. FastSpecFit value-added catalogue, or our own full-spectrum fits with stated templates).

Method: fixed before fitting: the SN sample and cuts, the host match radius and the fibre-aperture correction (the DESI fibre samples the host's centre, not the SN's site: the local-versus-global distinction is stated), the age estimator (light-weighted age from full-spectrum fitting, with its posterior), the model (Hubble residual = a + b × age + c × mass step, with the mass step fitted jointly so the age slope is not a proxy for mass). Interval from the fit with age uncertainties propagated (errors-in-variables). Nulls: ages shuffled among hosts; residuals against the age of a random DESI galaxy at the same redshift.

Expected precision: a few hundred SNe with DESI host spectra; the slope to about +-0.01 mag/Gyr, enough to test the claimed ~0.03 mag/Gyr.
