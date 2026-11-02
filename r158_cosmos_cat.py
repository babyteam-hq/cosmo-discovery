"""R-158: COSMOS sample. ZEST (Scarlata+07, ACS F814W, I<24) x COSMOS2020 Classic LePhare photo-z, matched at 0.5"."""
import numpy as np
from astropy.table import Table, join
from astropy.coordinates import SkyCoord
import astropy.units as u
from r158_common import *
fz, fc = os.path.join(DATA, "zest.vot"), os.path.join(DATA, "c2020.vot")
t0 = time.time()
if not os.path.exists(fz):
    tap_async("SELECT ra,dec,acs_mag_auto,acs_class_star,acs_clean,acs_mask,acs_elongation,acs_theta_image,r50,r80,rpet,sersic_n_gim2d,ell_gim2d,type,bulg,stellarity,junkflag,acstile FROM cosmos_morph_zurich_1", fz)
if not os.path.exists(fc):
    tap_async("SELECT id,ra,dec,lp_zbest,lp_zpdf_l68,lp_zpdf_u68,lp_type,flag_combined FROM cosmos2020c WHERE lp_zbest > 0.2 AND lp_zbest < 1.4", fc)
print("fetch s", round(time.time() - t0))
z = read_vot(fz, "ra,dec,acs_mag_auto,acs_class_star,acs_clean,acs_mask,acs_elongation,acs_theta_image,r50,r80,rpet,sersic_n_gim2d,ell_gim2d,type,bulg,stellarity,junkflag,acstile")
c = read_vot(fc, "id,ra,dec,lp_zbest,lp_zpdf_l68,lp_zpdf_u68,lp_type,flag_combined")
print(len(z), len(c))
cz = SkyCoord(z["ra"], z["dec"], unit="deg"); cc = SkyCoord(c["ra"], c["dec"], unit="deg")
i, d2, _ = cz.match_to_catalog_sky(cc)
ok = d2 < 0.5 * u.arcsec
for k in ["lp_zbest", "lp_zpdf_l68", "lp_zpdf_u68", "lp_type", "flag_combined"]:
    z[k] = np.where(ok, np.array(c[k])[i], np.nan)
z["acstile"] = np.array(z["acstile"]).astype(str)
z.write(os.path.join(DATA, "cosmos_zest_z.fits"), overwrite=True)
print("matched", ok.sum())
