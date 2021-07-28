import os
import math
import numpy as np
from datetime import datetime, timedelta
import spiceypy as spice
import astropy.units as u
from astropy import constants as const
from stix.spice import helio as hsp
from stix.spice import stix_datetime
from stix.utils import bson


solo_spice_min_unix= stix_datetime.utc2unix('2020-02-10T05:00:00Z')


def compute_earth_sun_so_angle(solo_positions):
    size = solo_positions.x.value.size
    vec_earth_sun = np.array([-1, 0, 0] * size).reshape(size, 3)

    vec_solo_sun = np.array([
        solo_positions.x.value, solo_positions.y.value, solo_positions.z.value
    ]).T
    fun = lambda v: math.sqrt(v.dot(v))
    mag = np.apply_along_axis(fun, 1, vec_solo_sun)
    product = np.array(
        [np.dot(v1, v2) for v1, v2 in zip(vec_earth_sun, vec_solo_sun)])
    return np.degrees(np.arccos(product / mag))

def get_sun_earth_light_time(obstime):
    dt=stix_datetime.utc2datetime(obstime)
    et=spice.datetime2et(dt)
    [state, ltime] = spice.spkezr( 'Earth', et,      'J2000', 'LT+S',   'Sun')
    return ltime

def get_solo_ephemeris(start_utc,
                       end_utc,
                       observer='SUN',
                       frame='GSE',
                       num_steps=200, return_json=False):
    '''
      calculate solo orbiter orbit using spice kernel data
      Args:
        start_utc: start_utc string
        end_utc:   end utc string
        frame:    coordinate frame supported by SPICE
        num_steps:  number of data points. 
      Returns:
        orbit data which is a python dictionary
    '''
    orbiter = hsp.Trajectory('Solar Orbiter')
    #starttime = stix_datetime.utc2datetime(start_utc)
    start_unix=stix_datetime.utc2unix(start_utc)
    end_unix=stix_datetime.utc2unix(end_utc)
    if start_unix< solo_spice_min_unix:
        start_unix= solo_spice_min_unix

    if end_unix<start_unix:
        end_unix=start_unix

    ut_space=np.linspace(start_unix, end_unix, num_steps)
    times = []
    utc_times = []
    for t in ut_space:
        dt=stix_datetime.unix2datetime(t)
        times.append(dt)
        utc_times.append(dt.strftime("%Y-%m-%dT%H:%M:%SZ"))
    result = {}
    try:
        orbiter.generate_positions(times, observer, frame)
        orbiter.change_units(u.au)
        dist_to_earth = np.sqrt((orbiter.x.value + 1)**2 + orbiter.y.value**2 +
                                orbiter.z.value**2)  #distance to earth
        light_time = (dist_to_earth * u.au).to(u.m) / const.c
        #time_to_earth = (1 * u.au).to(u.m) / const.c
        #time_to_earth = get_sun_earth_light_time(start_utc)
        time_to_solo = orbiter.r.to(u.m) / const.c
        sun_open_angle=const.R_sun.to(u.m)/orbiter.r.to(u.m)

        
        sun_angular_diameter_arcmin=np.degrees(np.arctan(sun_open_angle.value))*60.*2

        sun_earth_ltime= [get_sun_earth_light_time(utc)  for utc in utc_times]
        lt_diff = np.array(sun_earth_ltime) - time_to_solo.value

        earth_sun_solo_angles = compute_earth_sun_so_angle(orbiter)
        elevation= np.degrees(np.arctan2(orbiter.z.value, orbiter.r.value))

        result = {
            'ref_frame':'GSE',
            'observer':'Sun',
            'aunit': 'deg',
            'lunit': 'au',
            'vunit':'km/s',
            'tunit':'s',
            'utc': utc_times,
            'x': orbiter.x.value,
            'y': orbiter.y.value,
            'z': orbiter.z.value,
            'sun_solo_r': orbiter.r.value,
            'earth_solo_r': dist_to_earth,
            'speed': orbiter.speed.value,
            'owlt': light_time.value,
            'sun_earth_ltime':sun_earth_ltime,
            'light_time_diff': lt_diff,
            'earth_sun_solo_angle': earth_sun_solo_angles,
            'sun_angular_diameter':sun_angular_diameter_arcmin,
            'elevation': elevation,
        }
    except Exception as e:
        raise
        result = {'error': str(e)}
    if return_json:
        return bson.dict_to_json(result)
    return result


if __name__ == '__main__':
    from pprint import pprint
    result = get_solo_ephemeris('2020-04-10T00:00:00', '2020-04-11T00:00:00')
    pprint(result)
