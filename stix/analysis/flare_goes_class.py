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
from stix.core import config
from stix.spice import solo
from stix.core import logger
from datetime import datetime

threshold = 0

logger = logger.get_logger()
mdb = db.MongoDB()
flare_db = mdb.get_collection('flares')

GOES_STIX_COEFFS = [-6.90, 0.07505,
                    0.07064]  # must used with bkg subtracted counts
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


def save_config(config: dict):
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
    config['type'] = 'GOES_STIX_COEFF'
    config['date'] = datetime.now()
    sw_config_db = mdb.get_collection('sw_config')
    sw_config_db.save(config)


def read_latest_config():
    """
    Read the latest GOES STIX coefficients from Mongodb
    Parameters
    ----
    dt: datetime

    Returns
    ---
    None if not found or a dictionary 

    """
    sw_config_db = mdb.get_collection('sw_config')
    docs = list(
        sw_config_db.find({
            'type': 'GOES_STIX_COEFF'
        }).sort('date', -1).limit(1))
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


def get_goes_info(start, end):
    """
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
            goes_class, bkg_low)


def find_goes_class_flares_in_file(file_id):
    logger.info(f'processing flares in file {file_id}')
    fdb = mdb.get_collection('flares')
    flares = flare_db.find({'run_id': file_id, 'hidden': False})
    if not flares:
        logger.error(f'Flare not found in file {file_id} !')
        return
    goes_class_list = []
    for doc in flares:
        peak_utc, goes_class, goes_estimated = get_flare_goes_class(doc)
        if peak_utc is None:
            continue

        goes_class_list.append((peak_utc, goes_class, goes_estimated))

    return goes_class_list


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

    invalid_result={'min': 'None', 'center': 'None', 'max': 'None'}
    cnts = counts * dsun**2
    if cnts <= 0:
        return invalid_result, invalid_result

    x = np.log10(cnts)

    def g(y): return sum([coeffs[i] * y**i for i in range(len(coeffs))])
    def f(y): return goes_flux_to_class(10**g(y), frac=False)

    error = default_error  # default error
    try:
        bin_range = error_lut['bins']
        errors = error_lut['errors']
        for b, e in zip(bin_range, errors):
            if b[0] <= x <= b[1]:
                error = e
                # print(error)
                break
    except (KeyError, IndexError, TypeError):
        pass

    goes_class = {
        'min': f(x - error),
        'center': f(x),
        'max': f(x + error),
        'parameters': {
            'coeff': coeffs,
            'error': error
        }
    }
    goes_flux = {
        'min': g(x - error),
        'center': g(x),
        'max': g(x + error),
    }
    return goes_class, goes_flux

def compute_flare_goes_for_flares(id_start, id_end):
    fdb = mdb.get_collection('flares')
    flares = flare_db.find({'_id': {"$gte":id_start, "$lte":id_end}, 'hidden': False})
    if not flares:
        logger.error(f'Flare not found !')
        return
    goes_class_list = []
    for doc in flares:
        peak_utc, goes_class, goes_estimated = get_flare_goes_class(doc)
        if peak_utc is None:
            continue
        goes_class_list.append((peak_utc, goes_class, goes_estimated))
    return goes_class_list

def get_flare_goes_class(doc):
    """
    estimate flare GOES class, and store data into MongoDB
    Parameters
    ---
    doc: dict
        flare document in MongoDB
    Returns
    ----
        None

    """
    start_unix = doc['start_unix']
    end_unix = doc['end_unix']
    start_utc = time_utils.unix2utc(start_unix)
    end_utc = time_utils.unix2utc(end_unix)
    peak_utc = doc['peak_utc']
    try:
        peak_counts = doc['LC_statistics']['lc0']['signal_max']
    except KeyError:
        logger.warning(
            f'Failed to calculate flare goes class for flare {doc["_id"]}, due to insufficient information'
        )
        return None, None, None

    if peak_counts <= threshold:
        logger.info(f"Ignored flares {doc['_id']}, peak counts < {threshold}")
        return None, None, None

    ephs = solo.SoloEphemeris.get_solo_ephemeris(peak_utc, peak_utc, 1)
    eph = get_first_element(ephs)
    dsun = eph.get('sun_solo_r', 1)

    try:
        delta_lt = eph['light_time_diff']
    except (KeyError, IndexError):
        delta_lt = 0
    peak_time_low, peak_flux_low, peak_time_high, peak_flux_high, goes_class, bkg_low = get_goes_info(
        start_unix + delta_lt, end_unix + delta_lt)
    bkg_subtracted_counts = peak_counts - doc['LC_statistics']['lc0'][
        'bkg_median']
    estimated_class, estimated_flux = estimate_goes_class(bkg_subtracted_counts, dsun,
                                          GOES_STIX_COEFFS,
                                          GOES_STIX_ERROR_LUT)

    # is_ok = goes_class > estimated_class[
    #    'min'] and goes_class < estimated_class['max']
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
            'estimated_flux': estimated_flux,
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

    if sys.argv[1]=='flareid' and len(sys.argv)==4:
        flare_id_start=int(sys.argv[2])
        flare_id_end=int(sys.argv[3])
        compute_flare_goes_for_flares(flare_id_start, flare_id_end)



    elif len(sys.argv) == 2:
        file_ids.append(int(sys.argv[1]))
    else:
        file_ids = [x for x in range(int(sys.argv[1]), int(sys.argv[2]) + 1)]
    for _id in file_ids:
        find_goes_class_flares_in_file(_id)
