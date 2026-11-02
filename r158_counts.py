"""R-158 step 0: Euclid Q1 counts of resolved disc galaxies at 0.3<z<1.2 (MER x morphology x PHZ, IRSA)."""
import sys, numpy as np
from astropy.table import Table
from r158_common import *
out = os.path.join(DATA, "euclid_q1_discs.vot")
if not os.path.exists(out):
    t0 = time.time()
    tap_async("""SELECT m.object_id, m.ra, m.dec, m.tileid, m.flux_detection_total, m.segmentation_area, m.ellipticity,
      m.semimajor_axis, m.position_angle, m.point_like_prob, m.spurious_prob, m.det_quality_flag,
      mo.sersic_sersic_vis_radius, mo.sersic_sersic_vis_axis_ratio, mo.sersic_sersic_vis_index, mo.has_spiral_arms_yes, mo.has_spiral_arms_no,
      mo.smooth_or_featured_featured_or_disk, p.phz_median, p.phz_70_int1, p.phz_70_int2
      FROM euclid_q1_mer_catalogue m JOIN euclid_q1_mer_morphology mo ON m.object_id=mo.object_id
      JOIN euclid_q1_phz_photo_z p ON m.object_id=p.object_id
      WHERE m.flux_detection_total > 2.29 AND mo.sersic_sersic_vis_radius > 0.4 AND p.phz_median > 0.25 AND p.phz_median < 1.3""", out)
    print("fetch s", round(time.time() - t0))
t = Table.read(out)
print(len(t), t.colnames)
