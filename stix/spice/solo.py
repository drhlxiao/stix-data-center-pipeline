import os
import math
import numpy as np
from datetime import datetime, timedelta
import spiceypy as spice
import astropy.units as u
from astropy import constants as const
from stix.spice import helio as hsp
from stix.spice import stix_datetime

solo_spice_min_unix= stix_datetime.utc2unix('2020-02-10T05:00:00Z')


def compute_earth_sun_so_angle(solo_sun):

    size = solo_sun.x.value.size
    vec_earth_sun = np.array([1, 0, 0] * size).reshape(size, 3)
    vec_solo_sun = np.array([
        solo_sun.x.value, solo_sun.y.value, solo_sun.z.value
    ]).T
    fun = lambda v: math.sqrt(v.dot(v))
    mag = np.apply_along_axis(fun, 1, vec_solo_sun)
    product = np.array(
        [np.dot(v1, v2) for v1, v2 in zip(vec_earth_sun, vec_solo_sun)])
    return np.degrees(np.arccos(product / mag))

def get_earth_spice_HEE(datetimes):
    ets=[spice.datetime2et(t) for t in datetimes]
    pos_vel, light_times= spice.spkezr('Earth', ets,      'SOLO_HEE_NASA', 'LT+S',   'Sun')
    positions = np.array(pos_vel)[:, :3] * u.km
    return {'light_times':np.array(light_times),'positions':positions.to(u.au)}

def get_solo_ephemeris(start_utc,
                       end_utc,
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
    observer='Earth'
    frame='SOLO_HEE_NASA'
    target='SOLO'
    orbiter_earth= hsp.Trajectory(target)
    orbiter_sun= hsp.Trajectory(target)
    earth_hee= hsp.Trajectory('Earth')
    #starttime = stix_datetime.utc2datetime(start_utc)
    start_unix=stix_datetime.utc2unix(start_utc)
    end_unix=stix_datetime.utc2unix(end_utc)
    if start_unix< solo_spice_min_unix:
        start_unix= solo_spice_min_unix

    if end_unix<start_unix:
        end_unix=start_unix

    step=(end_unix-start_unix)/num_steps
    if step>12*3600:
        num_steps=int((end_unix-start_unix)/(12*3600))
    

    ut_space=np.linspace(start_unix, end_unix, num_steps)
    times = []
    utc_times = []
    for t in ut_space:
        dt=stix_datetime.unix2datetime(t)
        times.append(dt)
        utc_times.append(dt.strftime("%Y-%m-%dT%H:%M:%SZ"))
    result = {}
    try:
        orbiter_earth.generate_positions(times, 'Earth', frame)
        orbiter_sun.generate_positions(times, 'SUN', frame)
        earth_hee.generate_positions(times, 'SUN', frame)

        orbiter_earth.change_units(u.au)
        orbiter_sun.change_units(u.au)


        solo_dist_to_earth = orbiter_earth.r.value



        sun_open_angle=const.R_sun.to(u.m)/orbiter_sun.r.to(u.m)

        sun_angular_diameter_arcmin=np.degrees(np.arctan(sun_open_angle.value))*60.*2

        lt_diff = earth_hee.light_times - orbiter_sun.light_times 

        earth_sun_solo_angles = compute_earth_sun_so_angle(orbiter_sun)
        elevation= np.degrees(np.arctan2(orbiter_sun.z.value, orbiter_sun.r.value))

        result = {
            'ref_frame':frame,
            'observer':observer,
            'aunit': 'deg',
            'lunit': 'au',
            'vunit':'km/s',
            'tunit':'s',
            'utc': utc_times,
            'x': -orbiter_sun.x.value,
            'y': -orbiter_sun.y.value,
            'z': -orbiter_sun.z.value,
            'sun_solo_r': orbiter_sun.r.value,
            'earth_solo_r': orbiter_earth.r.value,
            'speed': orbiter_sun.speed.value,
            'owlt': orbiter_earth.light_times,
            #'sun_earth_ltime':sun_earth_ltime,
            'light_time_diff': lt_diff,
            'earth_sun_solo_angle': earth_sun_solo_angles,
            'sun_angular_diameter':sun_angular_diameter_arcmin,
            'elevation': elevation,
        }
    except Exception as e:
        raise
        print(e)
        result = {'error': str(e)}
    return result


if __name__ == '__main__':
    from pprint import pprint
    result = get_solo_ephemeris('2020-04-10T00:00:00', '2020-04-11T00:00:00')
    pprint(result)
