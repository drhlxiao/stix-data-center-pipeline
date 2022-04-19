import astropy.units as u
import numpy as np

from astropy.coordinates import SkyCoord
from datetime import datetime as dt
from sunpy.coordinates import Helioprojective
from sunpy.coordinates.frames import HeliocentricEarthEcliptic, HeliographicStonyhurst
from sunpy.map.maputils import _verify_coordinate_helioprojective, coordinate_is_on_solar_disk
from sunpy.map import Map, make_fitswcs_header
import drms
from astropy.wcs import WCS
import spiceypy
from astropy.time import Time

        
def jointly_visible(date_obs: dt,flare_loc_obs: tuple, obs_in: str = 'SOLO', obs_out: str = 'Earth') -> SkyCoord:
    '''Check if a solar flare whose location is given in one spacecraft's observation fram is visible from another spacecraft
    
    Inputs:
    
    date_obs: datetime
        Datetime (at earth) for which observation took place
    flare_loc_obs: tuple
        Tuple (x,y) in arcsec of flare location as seen by observing spacecraft
    obs_in: str or HeliocentricEarthEcliptic, default='SOLO'
        Name (SPICE identifier) of observing spacecraft, or HEE frame with spacecraft trajectory informatino
    obs_out: str or HeliocentricEarthEcliptic, default='Earth'
        Name of spacecraft for which joint visibility should be calculated, or HEE frame with spacecraft trajectory informatino. For flares close to the AIA limb, using Earth as the observer rather than AIA can result in an error of up to 10" in x
    
    Returns:
    
    flare_coord_out: SkyCoord
        HPC coordinate of flare as seen by obs_out spacecraft
    '''
    
    if isinstance(obs_in, HeliocentricEarthEcliptic):
        hee_in=obs_in
    else:
        hee_in = coordinates_body(date_obs, obs_in)
        
    if isinstance(obs_out, HeliocentricEarthEcliptic):
        hee_out=obs_out
    else:
        hee_out = coordinates_body(date_obs, obs_out)

    flare_coord_in = SkyCoord(*flare_loc_obs*u.arcsec,obstime=date_solo,observer=hee_in.transform_to(HeliographicStonyhurst(obstime=date_obs)),frame='helioprojective')

    flare_refcoord= SkyCoord(0*u.arcsec,0*u.arcsec,obstime=date_obs,observer=hee_out.transform_to(HeliographicStonyhurst(obstime=date_obs)),
                               frame='helioprojective')
                               
    flare_coord_out= flare_coord_in.transform_to(flare_refcoord.frame)

    #is flare_coord_out on-disk?
    if coordinate_is_on_solar_disk(flare_coord_out):
        return flare_coord_out
    else:
        raise ValueError(f"The input coordinate {flare_coord_in} is behind the limb from the view of {obs_out}.")

def is_visible_from_earth(date_solo: dt,flare_loc_solo: tuple) -> bool:
    '''wrapper returns bool '''
    try:
        _ = jointly_visible(date_solo,flare_loc_solo)
        return True
    except ValueError:
        return False

def is_visible_from_SO(date_earth: dt,flare_loc_earth: tuple,obs_out='SOLO') -> bool:
    '''wrapper returns bool '''
    try:
        _ = visible_from_SO(date_earth,flare_loc_earth,obs_in='Earth',obs_out=obs_out)
        return True
    except ValueError:
        return False
                
def coordinates_body(date_body,body_name,light_time=False):
    """
    Derive the coordinates of the given celestial body and then return them in
    Heliocentric Earth Ecliptic (HEE) coordinates. Optionally also return light travel time between Earth and body
    """

    # Observing time
    obstime = spiceypy.datetime2et(date_body)

    # Obtain the coordinates of Solar Orbiter
    if body_name == 'SOLO' or body_name == 'SO' or body_name == 'Solar Orbiter':
        ref_frame = 'SOLO_HEE_NASA'
    else:
        ref_frame = 'HEE'
    hee_spice, lighttimes = spiceypy.spkpos(body_name, obstime,
                                     ref_frame, #  Reference frame of the output position vector of the object
                                     'NONE', 'SUN')
    hee_spice = hee_spice * u.km

    # Convert the coordinates to HEE
    body_hee = HeliocentricEarthEcliptic(hee_spice,
                                          obstime=Time(date_body).isot,
                                          representation_type='cartesian')
    if not light_time:
        # Return the HEE coordinates of the body
        return body_hee
    else:
        return body_hee,lighttimes
                
def get_observer(date_in,obs='Earth',wcs=True,sc=False,wlen=1600,out_shape=(4096,4096),scale=None,rsun=False,hee=False):
    '''Get observer information. Get WCS object if requested. Return Observer as SkyCoord at (0",0") if requested.'''
    if obs in ['AIA','SDO']:
        observer,wcs_out=get_AIA_observer(date_in,wcs=True,sc=sc,wlen=wlen)
    else:
        if not hee:
            hee = coordinates_body(date_in, obs)
        observer=hee.transform_to(HeliographicStonyhurst(obstime=date_in))
        if sc:
            if rsun: #explicity set rsun
                observer = SkyCoord(0*u.arcsec, 0*u.arcsec,obstime=date_in,observer=observer,rsun=rsun,frame='helioprojective')
            else:
                observer = SkyCoord(0*u.arcsec, 0*u.arcsec,obstime=date_in,observer=observer,frame='helioprojective')

    if wcs == False:
        return observer
    else:
        try:
            return observer,out_wcs
        except UnboundLocalError:
            out_wcs=get_wcs(date_in,observer,out_shape=out_shape,scale=scale)
            return observer,out_wcs


def get_AIA_observer(date_obs,sc=False,wcs=True,wlen=1600):
    '''downloading required keywords from FITS headers. Resulting coordinates can differ from using Earth as observer by as much as 3" (on the limb) but usally will not be more than .5" otherwise.'''
    date_obs_str=dt.strftime(date_obs,"%Y.%b.%d_%H:%M:%S")
    client = drms.Client()
    kstr='CUNIT1,CUNIT2,CRVAL1,CRVAL2,CDELT1,CDELT2,CRPIX1,CRPIX2,CTYPE1,CTYPE2,HGLN_OBS,HGLT_OBS,DSUN_OBS,RSUN_OBS,DATE-OBS,RSUN_REF,CRLN_OBS,CRLT_OBS,EXPTIME,INSTRUME,WAVELNTH,WAVEUNIT,TELESCOP,LVL_NUM,CROTA2'
    if wlen in [1600,1700]:
        qstr=f'aia.lev1_uv_24s[{date_obs_str}/1m@30s]'
    else:
        qstr=f'aia.lev1_euv_12s[{date_obs_str}/1m@30s]'
    df = client.query(qstr, key=kstr)#'aia.lev1_euv_12s[2018.01.01_TAI/1d@12h] #'aia.lev1_euv_12s[2018.01.01_05:00:20/1m@30s]'
    if df.empty: #event not in database yet or other error
        observer=get_observer(date_obs,obs='Earth',sc=sc)
        import warnings
        warnings.warn(f"FITS headers for {date_obs} not available, using Earth observer" )
    else:
        try:
            meta=df.where(df.WAVELNTH==wlen).dropna().iloc[0].to_dict()
        except IndexError:
            meta=df.iloc[0].to_dict()
        if np.isnan(meta['CRPIX1']):
            meta['CRPIX1']=0.
        if np.isnan(meta['CRPIX2']):
            meta['CRPIX2']=0.
        fake_map=Map(np.zeros((10,10)),meta) #could probably do a faster way but this is convenient
        observer=fake_map.coordinate_frame.observer
    
    if sc:
        observer = SkyCoord(0*u.arcsec, 0*u.arcsec,obstime=date_obs,observer=observer,frame='helioprojective')

    if wcs:
        #while we're here and it's convenient...
        if df.empty:
            #wcs=get_wcs(date_obs, 'Earth')
            wcs=observer[1]
            observer=observer[0]
        else:
            wcs=WCS(meta)
        return observer,wcs
    return observer
    
def get_wcs(date_obs,obs_body,out_shape=(4096,4096),scale=None):
    '''generalized way to get WCS'''
    if isinstance(obs_body,str):
        obs_body=get_observer(date_obs, obs=obs_body)
    elif isinstance(obs_body,SkyCoord):
        ref_coord=obs_body
    else:
        ref_coord = SkyCoord(0*u.arcsec,0*u.arcsec,obstime=date_obs,
                              observer=obs_body,frame='helioprojective')
                              
    if scale:
        scale=(1., 1.)*ref_coord.observer.radius/u.AU*u.arcsec/u.pixel
    out_header = make_fitswcs_header(out_shape,ref_coord,scale=scale)
    return WCS(out_header)



