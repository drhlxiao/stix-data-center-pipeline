"""
 get flare ephemeris for stix flares
"""

import sys
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from datetime import datetime
import numpy as np
import spiceypy as spice
from sunpy.coordinates.frames import HeliocentricEarthEcliptic, HeliographicStonyhurst
from stix.spice import stix_datetime
from stix.spice import spice_manager
from astropy import constants as const


def get_ang_vectors(v1, v2):
    vector1 = np.array([v1.x.value, v1.y.value, v1.z.value])
    vector2 = np.array([v2.x.value, v2.y.value, v2.z.value])
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    angle = np.arccos(dot_product)  #angle in radian
    return np.degrees(angle)


def get_flare_spice(sun_x, sun_y, flare_utc, observer='earth'):
    #sun_x, sun_y in unites of arcsec
    date_flare = stix_datetime.utc2datetime(flare_utc)
    et_stix = spice.datetime2et(date_flare)
    flare_coord = [sun_x, sun_y]
    # Obtain the coordinates of Solar Orbiter
    earth_hee_spice, ltime_earth_sun = spice.spkpos(
        'EARTH',
        et_stix,
        'SOLO_HEE_NASA',  #  Reference frame of the output position vector of the object 
        'NONE',
        'SUN')

    earth_hee_spice = earth_hee_spice * u.km
    # Convert the coordinates to HEE
    earth_hee = HeliocentricEarthEcliptic(earth_hee_spice,
                                          obstime=Time(date_flare).isot,
                                          representation_type='cartesian')
    solo_hee_spice, ltime_sun_solo = spice.spkpos('SOLO', et_stix,
                                                  'SOLO_HEE_NASA', 'NONE',
                                                  'SUN')
    solo_hee_spice = solo_hee_spice * u.km
    # Convert the coordinates to HEE
    solo_hee = HeliocentricEarthEcliptic(solo_hee_spice,
                                         obstime=Time(date_flare).isot,
                                         representation_type='cartesian')
    if observer == 'earth':
        obs_coord = earth_hee.transform_to(
            HeliographicStonyhurst(obstime=date_flare))
    else:
        obs_coord = solo_hee.transform_to(
            HeliographicStonyhurst(obstime=date_flare))

    flare_skycoord = SkyCoord(flare_coord[0] * u.arcsec,
                              flare_coord[1] * u.arcsec,
                              obstime=date_flare,
                              observer=obs_coord,
                              frame='helioprojective')

    flare_hee = flare_skycoord.transform_to(
        HeliocentricEarthEcliptic(obstime=date_flare))
    v_flare_earth = earth_hee.cartesian - flare_hee.cartesian
    flare_hee_cart = flare_hee.cartesian
    flare_earth_r = np.sqrt(v_flare_earth.x**2 + v_flare_earth.y**2 +
                            v_flare_earth.z**2)

    v_flare_solo = -flare_hee.cartesian + solo_hee.cartesian
    flare_solo_r = np.sqrt(v_flare_solo.x**2 + v_flare_solo.y**2 +
                           v_flare_solo.z**2)
    ves = solo_hee.cartesian - earth_hee.cartesian
    v_earth_solo = np.sqrt(ves.x**2 + ves.y**2 + ves.z**2)

    flare_earth_t = flare_earth_r.to(u.m) / const.c
    flare_solo_t = flare_solo_r.to(u.m) / const.c
    owlt = v_earth_solo.to(u.m) / const.c

    flare_theta_stix = get_ang_vectors(flare_hee_cart, v_flare_solo)
    flare_theta_earth = get_ang_vectors(flare_hee_cart, v_flare_earth)

    return {
        'flare_earth_lt': flare_earth_t.value,
        'flare_solo_lt': flare_solo_t.value,
        'owlt': owlt.value,
        'dt': (flare_earth_t - flare_solo_t).value,
        'flare_solo_r': flare_solo_r.to(u.au).value,
        'dt_solar_center': ltime_earth_sun - ltime_sun_solo,
        'flare_utc': flare_utc,
        'theta_flare_norm_earth_deg': flare_theta_earth,
        'theta_flare_norm_solo_deg': flare_theta_stix,
        'earth_sun_ltime': ltime_earth_sun,
        'sun_solo_ltime': ltime_sun_solo,
        'observer': observer
    }
