# using the chi2 minimization algorithm to find the best fit solar flare location

import math
import random
import numpy as np
from stix.flare_pipeline import skylut
from stix.spice import solo

#open area ratio obtained by  the simulator
#detectorIndex: ratio

INF = 1e10  #math.inf is not supported by javascript
SAMPLE_POINTS = 200 * 200


def get_measured_cfl_pixel_open_area(counts, counts_error, det_fluence,
                                     det_fluence_error):
    """
        Estimate the open area for each pixel 
        Arguments:
         counts: CFL counts
         counts_error: CFL counts error
         det_fluence: detector fluence, in units of counts /mm2
    """
    if det_fluence <= 0:
        return np.zeros_like(counts), np.zeros_like(counts_error) + INF

    ce = np.array(counts)
    ce_err = np.array(counts_error)

    if np.sum(ce_err) == 0:
        ce_err = np.sqrt(ce)

    cfl_area = ce / det_fluence
    cfl_area_err = ce_err / det_fluence

    big_pixel_area = 9.6685
    max_area = np.max(cfl_area)
    if max_area > big_pixel_area:
        #this can happen due to statistics fluctuation
        cfl_area, cfl_area_err = cfl_area * big_pixel_area / max_area, cfl_area_err * big_pixel_area / max_area

    return cfl_area, cfl_area_err


def fit_location(counts,
                 count_errors,
                 mean_fluence,
                 mean_fluence_error,
                 flare_utc,
                 use_small_pixels=True,
                 use_det_fluence=True):
    #counts: CFL counts
    #mean_fluence:  detector summed  fluence if no grid counts/mm2
    """
      determine flare location using the least-mean squares method
      counts: 1x12 ndarray
            counts of CFL
      mean_fluence:
           detected averaged mean fluence  (counts/mm2) assuming there is no front and rear grid
      flare_utc:  str
            flare UTC
      
      Returns
        dictionary
        include flare location, chi2 map, etc.
    """

    max_pixel_index = 12 if use_small_pixels else 8

    obs_area, obs_area_err = get_measured_cfl_pixel_open_area(
        counts, count_errors, mean_fluence, mean_fluence_error)

    solo_location = solo.SoloEphemeris.get_solo_ephemeris(flare_utc, flare_utc, 1)
    sun_angular_diameter = solo_location.get('sun_angular_diameter', 0) * 60
    try:
        sun_angular_radius = sun_angular_diameter * 0.5
    except TypeError:
        sun_angular_radius = sun_angular_diameter[0] * 0.5
    
    X = skylut.lut[:, 0]
    Y = skylut.lut[:, 1]
    abs_X = np.abs(X)
    abs_Y = np.abs(Y)
    margin = 1.3

    pat = skylut.lut[:, 2:]
    #chi2=np.apply_along_axis(get_chi2, 1, pat, use_small_pixels, obs_area, obs_area_err)
    #replaced with the line below. it is x10 faster
    chi2 = np.sum(np.divide((obs_area - pat)**2,
                            obs_area_err**2,
                            where=obs_area_err != 0),
                  axis=1)

    min_idx = np.argmin(chi2)
    best_fit_area = pat[min_idx]
    min_chi2 = chi2[min_idx]
    delta_chi2 = chi2 - min_chi2

    sel = delta_chi2 < 50
    #only sends a part of the chi2map to client
    _delta_chi2 = delta_chi2[sel]
    chi2map_xy = skylut.lut[sel, 0:2]
    chi2map_x = chi2map_xy[:, 0]
    chi2map_y = chi2map_xy[:, 1]
    min_pos = skylut.lut[min_idx, 0:2]
    #nn_xy = cfl_nn.predict(counts, mean_fluence)
    res = {
        'min_pos': min_pos.tolist(),
        'min_chi2': min_chi2,
        'cfl_open_area_best_fit': best_fit_area,  #before normalized
        'cfl_open_area_obs': obs_area.tolist(),
        'cfl_open_area_obs_err': obs_area_err.tolist(),
        'sun_angular_diameter': sun_angular_diameter,
        'angle_unit': 'arcsec',
        #'nn_xy': nn_xy,
        'grid': {
            'x': skylut.X.tolist(),
            'y': skylut.Y.tolist()
        },
        'delta_chi2': {
            'x': chi2map_x.tolist(),
            'y': chi2map_y.tolist(),
            'z': _delta_chi2.tolist()
        }
    }
    return delta_chi2, res
