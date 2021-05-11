import os
import math
import numpy as np
from datetime import datetime, timedelta
import spiceypy as spice
import astropy.units as u
from astropy import constants as const
from stix.spice import helio as hsp
from stix.spice import stix_datetime

from spice.spice_manager import SpiceManager
spice_manager = SpiceManager.get_instance()
solo_spice_min_datetime = stix_datetime.utc2datetime('2020-02-10T05:00:00Z')


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


def get_solo_ephemeris(start_utc,
                       end_utc,
                       observer='SUN',
                       frame='GSE',
                       num_steps=200):
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
    starttime = stix_datetime.utc2datetime(start_utc)

    if starttime < solo_spice_min_datetime:
        starttime = solo_spice_min_datetime
    endtime = stix_datetime.utc2datetime(end_utc)
    if endtime < starttime:
        endtime = starttime

    times = []
    utc_times = []
    span = endtime - starttime

    step_seconds = span.total_seconds() / num_steps
    step = timedelta(minutes=1)
    m, s = divmod(step_seconds, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    if d > 0 or h > 0 or m > 0:
        step = timedelta(days=d, hours=h, minutes=m)

    while starttime <= endtime:
        times.append(starttime)
        utc_times.append(starttime.strftime("%Y-%m-%dT%H:%M:%SZ"))
        starttime += timedelta(days=1)
    result = {}
    try:
        orbiter.generate_positions(times, observer, frame)
        orbiter.change_units(u.au)
        dist_to_earth = np.sqrt((orbiter.x.value + 1)**2 + orbiter.y.value**2 +
                                orbiter.z.value**2)  #distance to earth
        light_time = (dist_to_earth * u.au).to(u.m) / const.c
        time_to_earth = (1 * u.au).to(u.m) / const.c
        time_to_solo = orbiter.r.to(u.m) / const.c
        sun_open_angle=const.R_sun.to(u.m)/orbiter.r.to(u.m)

        
        sun_angular_diameter_arcmin=np.degrees(np.arctan(sun_open_angle.value))*60.*2

        lt_diff = time_to_earth.value - time_to_solo.value

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
            'light_time_diff': lt_diff,
            'earth_sun_solo_angle': earth_sun_solo_angles,
            'sun_angular_diameter':sun_angular_diameter_arcmin,
            'elevation': elevation,
        }
    except Exception as e:
        print(e)
        result = {'error': str(e)}
    return result


if __name__ == '__main__':
    from pprint import pprint
    result = get_solo_ephemeris('2020-04-10T00:00:00', '2020-04-11T00:00:00')
    pprint(result)
