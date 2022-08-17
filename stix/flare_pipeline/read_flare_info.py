"""
 imaging inputs gen
"""
import os
import json
import numpy as np
from datetime import datetime, timedelta

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits

from sunpy import map
from sunpy.map import make_fitswcs_header
from sunpy.coordinates.frames import HeliocentricEarthEcliptic,HeliographicStonyhurst

from stix.spice import solo
#TO install stixdcpy: pip install git+https://github.com/drhlxiao/stixdcpy.git
##### Constants


def generate_idl_inputs(flare_id, 
    energy_range_science_channel_upper_limit: int,
    energy_range_science_channel_lower_limit: int,
    vis_fwdfit_map_filename: str,
    bp_map_filename: str,
    stx_quicklook_fwdfit_path: str):
    """
    Arguments:
        background_L1A_fits_filename: str,
            path of the L1A background file
        flaring_data_L1A_fits_filename: str,
            path of the L1A science file
        flare_start_UTC: str,
            start of the time interval to consider for image reconstruction
        flare_end_UTC: str,
            end of the time interval to consider for image reconstruction
        energy_range_science_channel_upper_limit: int,
            lower bound of the energy range to consider for image reconstruction
        energy_range_science_channel_upper_limit: int,
            lower bound of the energy range to consider for image reconstruction
        vis_fwdfit_map_filename: str,
            path and name of the VIS_FWDFIT_PSO map reconstructed from STIX data
        bp_map_filename: str,
            path and name of the Back Projection map reconstructed from STIX data
        ssw_home: str,
            path of the SSW folder
        idl_home: str,
            path of the IDL installation folder
        stx_quicklook_fwdfit_path: str,
            path of the folder containing 'stx_quicklook_fwdfit.pro'
    Returns:
        Fits file of the FWDFIT_PSO map of the flaring event
    """
    
    this_time = Time(flare_start_UTC).to_datetime()
    #solo_hee = coordinates_body(this_time, observer)
    stix_aux=solo.SoloEphemeris.get_solar_limb_stix_fov(flare_start_UTC)
    try:
        apparent_radius_sun = stix_aux['sun_angular_diameter']*0.5*60 #in units of arcsec
        roll_angle_solo = stix_aux['roll']  #roll angle in degrees
        hee_spice=stix_aux['solo_hee']
    except KeyError:
        raise ValueError('No auxiliary data available for the requested time')

    hee_spice = np.array(hee_spice[0]) * u.km
    # Convert the coordinates to HEE
    solo_hee = HeliocentricEarthEcliptic(hee_spice,
                                          obstime=Time(this_time).isot,
                                          representation_type='cartesian')

    solo_hgs = solo_hee.transform_to(HeliographicStonyhurst(obstime=this_time))

    B0 = solo_hgs.lat.deg # Heliographic latitude (B0 angle)
    L0 = solo_hgs.lon.deg # Heliographic longitude (L0 angle)
    #ROLL ANGLE Solar Orbiter in degree
    
    # RSUN
    
    
    # Define dictionary containing the inputs for the IDL code
    inputs = {'background_filename': background_L1A_fits_filename,
              'signal_filename': flaring_data_L1A_fits_filename,
              'signal_time_range':[flare_start_UTC, flare_end_UTC],
              'flare_end_UTC': flare_end_UTC,
              'energy_range_sci_channel': [ energy_range_science_channel_lower_limit, energy_range_science_channel_upper_limit],
              'vis_fwdfit_map_filename': vis_fwdfit_map_filename,
              'bp_map_filename': bp_map_filename, 
              'L0': L0,
              'B0': B0,
              'apparent_radius_sun': apparent_radius_sun,
              'roll_angle_solo': roll_angle_solo
              }
            
    with open('inputs.json','w') as f:
        json.dump(inputs,f)
        print('inputs have been written to inputs.json')
    

