#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Flare processing
   - attach ephemeris data to flares
   - determine goes class
   - store results in the flare list database
"""
import numpy as np
from stix.core import mongo_db as db
from stix.spice import time_utils
from stix.utils import bson
from stix.core import config, logger
from stix.spice import solo
from datetime import datetime

from astropy import units as u
from sunpy.coordinates.frames import HeliocentricEarthEcliptic
from hek_event_handler import HEKEventHandler
from stix_goes_fit import STIX_GOES_fit

threshold = 0

logger = logger.get_logger()
mdb = db.MongoDB()
flare_db = mdb.get_collection('flares')

GOES_STIX_COEFFS = [-6.90, 0.07505,
                    0.07064]  #must used with bkg subtracted counts
GOES_STIX_ERROR_LUT = {'bins': [(0, 10)], 'errors': [0.3]}
"""
Erica's data
GOES_STIX_COEFFS=[-8.0780854936,0.6641698797] 
GOES_STIX_ERROR_LUT={'bins':[
(2.0,2.5),(2.5,3.0),(3.,3.5),(3.5,4.0),(4,4.5),(4.5 5.0)
],
'errors':[0.2404348065,0.2850726823,0.1928913381,0.277,0.139,0.273,0.217]
}
"""

def save_config(config:dict):
    """
    save GOES STIX coefficients in mongodb
    Parameters
    ---
    config: dict
        it has a structure like the one below
        {'coeff': [-6.9, 0.07, 0.07],
        'error_bar_lut':{'bins': [(0, 10)], 'errors': [0.3]}
        }
    """
    config['type']='GOES_STIX_COEFF'
    config['date']=datetime.now()
    sw_config_db=mdb.get_collection('sw_config')
    sw_config_db.save(config)
    
    
def read_lastest_config():
    """
    Read the latest GOES STIX coefficients from Mongodb
    Parameters
    ----
    dt: datetime

    Returns
    ---
    None if not found or a dictionary 
        
    """
    sw_config_db=mdb.get_collection('sw_config')
    docs=list(sw_config_db.find({'type':'GOES_STIX_COEFF'}).sort('date',-1).limit(1))
    if docs:
        return docs[0]
    return None


def get_first_element(obj):
    if not isinstance(obj, dict):
        return obj
    new_obj = {}
    for key, val in obj.items():
        if isinstance(val, list):
            new_obj[key] = val[0]
        else:
            new_obj[key] = val
    return new_obj


def reconstruct_hee(eph_list,obstime):
    '''Reconstruct SkyCoord from ephemeris list of coordinates, assuming units are km
    eph_list:list
        list of ephemeris coordinates [x,y,z]
    obstime:  (tuple, list, str, pandas.Timestamp, pandas.Series, pandas.DatetimeIndex, datetime.datetime, datetime.date, numpy.datetime64, numpy.ndarray, astropy.time.Time)
        The time of the observation. This is used to determine the position of solar-system bodies (e.g., the Sun and the Earth) as needed to define the origin and orientation of the frame.
    Returns:
        eph_hee: HeliocentricEarthEcliptic
    '''
    return HeliocentricEarthEcliptic(*clist*u.km,obstime=obstime,representation_type='cartesian')

def goes_flux_to_class(x, frac=True):
    x = float(f'{x:.1e}')
    if x == 0:
        return 'NA'
    elif x < 1e-7:
        return 'A'
    elif x < 1e-6:
        return f'B{x/1e-7:.1f}' if frac else f'B{x/1e-7:.0f}'
    elif x < 1e-5:
        return f'C{x/1e-6:.1f}' if frac else f'C{x/1e-6:.0f}'
    elif x < 1e-4:
        return f'M{x/1e-5:.1f}' if frac else f'M{x/1e-5:.0f}'
    else:
        return f'X{x/1e-4:.1f}' if frac else f'X{x/1e-4:.0f}'


def get_goes_info(start, end,peak_time=None,hek=False,eph_hee=False):
    """
    start: float
        Unix start time
    end: float
        Unix end time
    peak_time: datetime, default=None
        Time of peak STIX counts
    hek: bool, default False
        Use HEK to get GOES flux corresponding to AIA flares visible from SOLO
    eph_hee: HeliocentricEarthEcliptic, default=False
        Observer frame of SOLO, for calculating if AIA-derived flare locations are visible from SOLO
    """
    data = mdb.get_goes_fluxes(start, end, satellite=16)
    last_time = 0
    start_times = []
    low_name = '0.1-0.8nm'
    high_name = '0.05-0.4nm'
    peak_time_low = 0
    peak_flux_low = 0
    peak_flux_high = 0
    peak_time_high = 0
    bkg_low = 0
    goes_source='GOES 16'
    
    if hek:
        if not peak_time:
            hek_time=unix2datetime(start)
        else:
            hek_time=peak_time
        hek_event=hek.HEKEventHandler(hek_time,obs_out='SOLO',hee_out=eph_hee, single_result=True).df #searches HEK for any AIA flares within 5 minutes of the (lighttime corrected) STIX peak time
        if hek.visible_from_SOLO[0]: #if AIA flare is also seen by STIX, use GOES flux at AIA flare peak time
            goes_source=f"HEK {hek_event.frm_name[0]}"
            event_peaktime=utc2unix(hek_event.event_peaktime[0])
            peak_time_low, peak_flux_low, peak_time_high, peak_flux_high=event_peaktime, goes_flux, event_peaktime, 0.
            if hek_event.fl_goescls[0]:
                goes_class=hek_event.goes_class[0]
            else:
                goes_class=goes_flux_to_class(hek_event.fl_peakflux[0])
            try:
                bkg_low=[d.get('background', 0) for d in data if d['energy'] == low_name and np.isclose(d['unix_time'],event_peaktime,atol=60)][0] #within 60 seconds of AIA peak time (with what time resolution are goes fluxes stored?)
            except IndexError:
                pass #keep bkg_low=0
        else:
            use_goes=True
    
    if not hek or use_goes:
        for d in data:
            unix = d['unix_time']
            flux = d['flux']
            if d['energy'] == low_name and flux > peak_flux_low:
                peak_flux_low = flux
                peak_time_low = unix
                bkg_low = d.get('background', 0)

            if d['energy'] == high_name and flux > peak_flux_high:
                peak_flux_high = flux
                peak_time_high = unix

        goes_class = goes_flux_to_class(peak_flux_low)
    return (peak_time_low, peak_flux_low, peak_time_high, peak_flux_high,
            goes_class, bkg_low,goes_source)


def find_goes_class_flares_in_file(file_id,hek=False):
    logger.info(f'processing flares in file {file_id}')
    fdb = mdb.get_collection('flares')
    flares = flare_db.find({'run_id': file_id, 'hidden': False})
    if not flares:
        logger.error(f'Flare not found in file {file_id} !')
        return
    goes_class_list = []
    for doc in flares:
        peak_utc, goes_class, goes_estimated = get_flare_goes_class(doc,hek=hek) #or hek event hnadler....
        goes_class_list.append((peak_utc, goes_class, goes_estimated))

    return goes_class_list

def update_fit(qdict={"visible_from_SOLO",True}):
    """Update STIX-GOES fit using all flares in the database that meet certain conditions"""
    fdb = mdb.get_collection('flares')
    #get flares that meet certain requirements - look angle, visible_from_SOLO=True, goes_source='HEK SSW Latest Events' etc
    #TO DO
    selected_flares= flare_db.find(qdict) #flare_db.find... I think it will return a list
    df=pd.DataFrame(selected_flares) #assuming list of dicts
    #any renaming necessary. stix_goes_fit needs:
    #peak_utc
    #_id
    #peak_counts_corrected (ideally background-subtracted peak counts, corrected to 1 AU)
    #GOES_flux (W m^2 s^1)
    
    fit=STIX_GOES_fit(df).do_fit().bin_residuals()
    config_dict={}
    config_dict['coeff']= [fit.slope,fit.intercept],
    config_dict['error_bar_lut']={}
    config_dict['error_bar_lut']['bins']=[(bin,bins[i+1]) for i,bin in enumerate(bins[:-1])]
    config_dict['error_bar_lut']['errors']=[fit[f'bin_{bin}_sigma'] for bin in fit.bins]
    config_dict['number_flares']=fit.nflares
    config_dict['r_value']=fit.rvalue
    save_config(config_dict)

def estimate_goes_class(counts: float,
                        dsun: float,
                        coeffs: list,
                        error_lut: dict,
                        default_error=0.3):
    """
    estimate goes class
       see https://pub023.cs.technik.fhnw.ch/wiki/index.php?title=GOES_Flux_vs_STIX_counts

    Parameters:
        peak_counts:  float
            STIX background subtracted peak counts
        dsun: float
            distance between solar orbiter and the sun in units of au
        coeff: list
            polynomial function coefficients
        error_lut: dict
            a look-up table contains errors  in log goes flux
            for example, errors_lut={'bins':[[0,0.5],[0.5,2]],'errors':[0.25,0.25]}
            bins contains bin edges and errors the corresponding error of the bin
        default_error: float
            default error in log goes flux when it can not be found in the lut
    Returns: dict
        predicted GOES class and limits
    """

    cnts = counts * dsun**2
    if cnts <= 0:
        return {'min': None, 'max': None, 'max': None}
    x = np.log10(cnts)

    g = lambda y: sum([coeffs[i] * y**i for i in range(len(coeffs))])
    f = lambda y: goes_flux_to_class(10**g(y), frac=False)

    error = default_error  #default error
    try:
        bin_range = error_lut['bins']
        errors = error_lut['errors']
        for b, e in zip(bin_range, errors):
            if b[0] <= x <= b[1]:
                error = e
                #print(error)
                break
    except (KeyError, IndexError, TypeError):
        pass

    result = {
        'min': f(x - error),
        'center': f(x),
        'max': f(x + error),
        'parameters': {
            'coeff': coeffs,
            'error': error
        }
    }
    return result


def get_flare_goes_class(doc,hek=False):
    """
    estimate flare GOES class, and store data into MongoDB
    Parameters
    ---
    doc: dict
        flare document in MongoDB
    hek: bool, default False
        Use HEK database to determine if flare is visible from both STIX and AIA. If so, return GOES class measured from AIA event peak time closest to STIX (lighttime corrected) event time. If not, return goes_class=None
    Returns
    ----
        peak_time_low: datetime
            UTC time
        goes_class: str
            GOES class ... either as seen at AIA event peak, or None
        estimated_class: str
            estimated GOES class for flare, using ...

    """
    start_unix = doc['start_unix']
    end_unix = doc['end_unix']
    start_utc = time_utils.unix2utc(start_unix)
    end_utc = time_utils.unix2utc(end_unix)
    peak_utc = doc['peak_utc']

    peak_counts = doc['LC_statistics']['lc0']['signal_max']

    if peak_counts <= threshold:
        logger.info(f"Ignored flares {doc['_id']}, peak counts < {threshold}")
        return

    ephs = solo.SoloEphemeris.get_solo_ephemeris(peak_utc, peak_utc, 1)
    eph = get_first_element(ephs)
    dsun = eph.get('sun_solo_r', 1)
    try:
        delta_lt = eph['light_time_diff']
    except (KeyError, IndexError):
        delta_lt = 0
    peak_utc_corrected=unix2datetime(utc2utix(peak_utc)+delta_lt)
    peak_time_low, peak_flux_low, peak_time_high, peak_flux_high, goes_class, bkg_low,goes_source = get_goes_info(
        start_unix + delta_lt, end_unix + delta_lt,peak_time=peak_utc_corrected,hek=hek,eph_hee=reconstruct_hee(eph.get('solo_hee',1)),peak_utc) #if hek pass in ephemeris too
    bkg_subtracted_counts = peak_counts - doc['LC_statistics']['lc0'][
        'bkg_median']
    estimated_class = estimate_goes_class(bkg_subtracted_counts, dsun,
                                          GOES_STIX_COEFFS,
                                          GOES_STIX_ERROR_LUT)

    is_ok = goes_class > estimated_class[
        'min'] and goes_class < estimated_class['max']
    updated_data = {
        'goes': {
            'low': {
                'unix_time': peak_time_low,
                'utc': time_utils.unix2utc(peak_time_low),
                'flux': peak_flux_low,
                'background': bkg_low
            },
            'high': {
                'unix_time': peak_time_high,
                'utc': time_utils.unix2utc(peak_time_high),
                'flux': peak_flux_high
            },
            'class': goes_class,
            'estimated_class': estimated_class,
            'goes_source':goes_source,
        },
        'ephemeris': eph,
    }
    flare_db.update_one({'_id': doc['_id']}, {'$set': updated_data})
    return time_utils.unix2utc(peak_time_low), goes_class, estimated_class


if __name__ == '__main__':
    import sys
    file_ids = []
    if len(sys.argv) < 2:
        print('flare_processing file_number')
    elif len(sys.argv) == 2:
        file_ids.append(int(sys.argv[1]))
    else:
        file_ids = [x for x in range(int(sys.argv[1]), int(sys.argv[2]) + 1)]
    try:
        hek=bool(sys.argv[3]) #use 3rd arg for hek=True or hek=False (use HEK database to get AIA flare locations and GOES flux/class values)
        if hek.lower() == 'true':
            hek=True
    except IndexError:
        hek=False
    for _id in file_ids:
        find_goes_class_flares_in_file(_id,hek=hek)
