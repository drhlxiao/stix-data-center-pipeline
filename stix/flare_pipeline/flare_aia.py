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
from stix.spice import time_utils as sdt
from stix.spice import spice_manager
from stix.spice.solo import SoloEphemeris
from astropy import constants as const
import drms
from sunpy.net import Fido
from sunpy.net import attrs as a

def get_ang_vectors(v1, v2):
    vector1 = np.array([v1.x.value, v1.y.value, v1.z.value])
    vector2 = np.array([v2.x.value, v2.y.value, v2.z.value])
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    angle = np.arccos(dot_product)  #angle in radian
    return np.degrees(angle)


def get_flare_spice(sun_x, sun_y, flare_utc, observer='earth'):
    """
        get flare auxiliary data from SPICE kernel

    """
    #sun_x, sun_y in unites of arcsec
    date_flare = sdt.utc2datetime(flare_utc)
    et_stix = spice.utc2et(date_flare)
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
    """
    if observer == 'earth':
        obs_coord = earth_hee#.transform_to(
            #HeliographicStonyhurst(obstime=date_flare))
    else:
        obs_coord = solo_hee#.transform_to(
            #HeliographicStonyhurst(obstime=date_flare))
            """
    obs_coord=earth_hee if observer=='earth' else solo_hee

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

    earth_flare_solo_deg= get_ang_vectors(v_flare_earth, v_flare_solo)


    res={
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
        'earth_flare_solo_deg': earth_flare_solo_deg,
        'sun_solo_ltime': ltime_sun_solo,
        'observer': observer
    }
    return res
    
def get_AIA_observer(date_obs):
    """Get reference coordinate for SDO/AIA coordinate frame. If no AIA frame information available, use Earth reference frame (can result in differences up to 10" at the limb)"""
    date_obs_str = date_obs.strftime("%Y.%b.%d_%H:%M:%S")
    client = drms.Client()
    kstr = 'HGLN_OBS,HGLT_OBS,DSUN_OBS,RSUN_OBS,DATE-OBS,RSUN_REF,EXPTIME,INSTRUME,WAVELNTH,WAVEUNIT,TELESCOP,LVL_NUM,CROTA2'
    qstr = f'aia.lev1_uv_24s[{date_obs}/1m@30s]' #if no UV data available can also try EUV: qstr=f'aia.lev1_euv_12s[{start_utc}/1m@30s]'
    try:
        df = client.query(qstr, key=kstr)
        meta = df.iloc[0].to_dict()
        obstime = meta['DATE-OBS']
        sdo_observer = SkyCoord(meta['HGLN_OBS']*u.deg,meta['HGLT_OBS']*u.deg,meta['DSUN_OBS']*u.m,frame=HeliographicStonyhurst,obstime=obstime)
    except Exception: #timeout - use Earth position as approximation then
        eph = SoloEphemeris()
        positions = eph.get_earth_spice_HEE([pd.to_datetime(stix_bp_map.meta['date_obs'])])['positions'][0]
        obstime = stix_bp_map.meta['date_obs']
        sdo_observer = SkyCoord(*positions,frame=HeliocentricEarthEcliptic,obstime=obstime,representation_type='cartesian')
    sdo_refcoord = SkyCoord(0*u.arcsec,0*u.arcsec,observer=sdo_observer,frame='helioprojective',obstime=obstime) #helioprojective instead of stonyhurst
    return sdo_refcoord

def find_cutout_coords(stix_bp_map,padding=10*u.arcsec):
    """Find the coordinates of AIA cutout that upon reprojection
    corresponds to the area of the STIX back projection map. Returns None if the STIX area is not visible from AIA
    
    Inputs:
    stix_bp_map: sunpy.map.Map
    (can also write this to use CFL flare position or any other input that contains obstime and central coordinate)
    padding: Quantity
    how much to pad the coordinates
    
    Outputs:
    bottom_left_coord: SkyCoord or None
    bottom left coord of cutout to request (AIA POV at obstime stix_bp_map.meta['date-obs'])
    
    top_right_coord: SkyCoord or None
    top right coord of cutout to request (AIA POV at obstime stix_bp_map.meta['date-obs'])
    
    For usage with:
    cutout = a.jsoc.Cutout(
    bottom_left_coord,
    top_right=top_right_coord,
    tracking=True
    )
    """
    # Get AIA pointing info
    sdo_refcoord = get_AIA_observer(pd.to_datetime(stix_bp_map.meta['date-obs'])) #warning no light-time correction
    
    # Get STIX map extent
    stix_bottom_left = stix_bp_map.bottom_left_coord
    stix_top_right = stix_bp_map.top_right_coord
    
    # STIX coords to AIA coords
    bottom_left_transformed = stix_bottom_left.transform_to(sdo_refcoord)
    top_right_transformed = stix_top_right.transform_to(sdo_refcoord)
    
    # Check for NaNs and deal with them
    bottom_left_xy, top_right_xy = [] ,[]
    for bl,tr in [(bottom_left_transformed.Tx, top_right_transformed.Tx), (bottom_left_transformed.Ty, top_right_transformed.Ty)]:
        if np.isnan(bl) and np.isnan(tr):
            # All off-disk or not jointly observed.
            return None, None
        elif np.isnan(bl):  # Replace with +-1228"
            bottom_left_xy.append(-1228*u.arcsec+padding)
        elif np.isnan(tr):
            top_right_xy.append(1228*u.arcsec-padding)
        else:
            bottom_left_xy.append(bl)
            top_right_xy.append(tr)
    
    # Pad
    bottom_left_coord = SkyCoord(bottom_left_xy[0] + np.sign(bottom_left_xy[0])*padding,
                                 bottom_left_xy[1] + np.sign(bottom_left_xy[1])*padding,
                                 frame=bottom_left_transformed.frame)
    top_right_coord = SkyCoord(top_right_xy[0] + np.sign(top_right_xy[0])*padding,
                                 top_right_xy[1] + np.sign(top_right_xy[1])*padding,
                                 frame=top_right_transformed.frame)

    return bottom_left_coord, top_right_coord
    
def download_AIA_cutout(time_int, wave = 1600, series='aia_lev1_uv_24s', cutout_coords=False, single_result=True, aia_file_path=None):
    """Download AIA cutout corresponding to STIX flare."""
    if type(time_int[0]) == str:
        time_int[0] = dt.strptime(time_int[0],'%Y-%m-%dT%H:%M:%S')
        time_int[1] = dt.strptime(time_int[1],'%Y-%m-%dT%H:%M:%S')
    
    wave = a.Wavelength(wave*u.angstrom)
    #instr= sn.attrs.Instrument(instrument)
    time = a.Time(time_int[0],time_int[1])
    qs = [time,wave, series]

    if cutout_coords != False:
        cutout = a.jsoc.Cutout(cutout_coords[0],top_right=cutout_coords[1])
    
    qs.append(cutout)
    #qs.append(a.jsoc.Notify(jsoc_email)) # Set this in config

#    if sample: #Quantity
#        sample = a.Sample(sample)
#        qs.append(sample)

    res = Fido.search(*qs)
    if single_result: # Take the first result
        res = res[0]

    if not aia_file_path: 
        files = Fido.fetch(res,path='./') #set path in config - destination for downloaded files
    else: 
        files = Fido.fetch(res,path=f"{aia_file_path}/")
    return res

if __name__=='__main__':
    
    download_AIA_cutout(['2022-01-01T00:00:00', '2022-01-01T01:00:00']) 
