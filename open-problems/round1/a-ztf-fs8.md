# Proposed test: the growth rate f sigma_8 at z ~ 0.03 from ZTF SN Ia DR2 peculiar velocities

Question: what is f sigma_8 in the local universe measured from the peculiar velocities of the Zwicky Transient Facility's second Type Ia supernova data release (about 1,000 SNe Ia at z < 0.06 with homogeneous selection over the northern sky)?

Data: ZTF SN Ia DR2 light curves, SALT2/SALT3 parameters and host data (public); the 2M++ and Cosmicflows-4 density fields for the velocity-density comparison; public N-body light-cone mocks or our own lognormal mocks matching the ZTF selection.

Method: fixed before fitting: the sample cuts (DR2 quality cuts, z < 0.06, |b| > 10 deg), the standardisation, the velocity estimator (log-distance ratio), the likelihood (the maximum-likelihood velocity covariance method, and the velocity-density comparison beta = f/b as the second estimator), and the treatment of the DR2 calibration caveat (fit with zero-point offsets per filter and field as nuisance parameters). Errors and bias from mocks with the ZTF selection and cadence; nulls from shuffled positions.

Expected precision: f sigma_8 to about 15-20 % from ~1,000 SNe, the forecast in the ZTF peculiar-velocity forecast paper.
