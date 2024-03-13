"""
 create flare images for science data
 Procedure to create image for science data
 Step 1:
    prepare inputs for imaging, store the data into flare_images database
"""
import os
import sys
import json
import numpy as np
import subprocess
import random
import time
import uuid
from datetime import datetime, timedelta

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits

from sunpy import map
from stix.core import config
from stix.core import mongo_db as db
from stix.analysis.science_l1 import ScienceL1

from stix.core import logger
from stix.spice import solo
from stix.core.energy_bins import StixEnergyBins
from stix.utils import bson
from stix.spice import time_utils as stu
from stix.spice import spice_manager as spm
from pathlib import Path

spm.spice.load_kernels(reload_all=True)

##### Constants
mdb = db.MongoDB()
logger = logger.get_logger()
quicklook_path = config.get_config('pipeline.daemon.flare_images')
bsd_db = mdb.get_collection('bsd')
flare_db = mdb.get_collection('flares')
fits_db = mdb.get_collection('fits')
flare_images_db = mdb.get_collection('flare_images')
req_db = mdb.get_collection('data_requests')
MIN_COUNTS = 10000
DATA_LEVEL_FOR_IMAGING='L1'
BKG_FILE_MAX_DAY_OFFSET=180

ASPECT_TIME_TOR = 300  #time tolerance for looking for aspect solution

last_bkg_fits_doc = None
TIME_MARGIN=4



def get_sun_center_from_aspect(y_srf, z_srf):
    """
    Frederic message on slack  Jun 22, 2022
    So, as it is done in the IDL imaging pipeline, the images have to be shifted by:
    3:06
    (Y_SRF + offset_X, -Z_SRF + offset_Y)
    3:07
    the (offset_X, offset_Y) is the residual systematic error that I've measured on a small sample of events already a year ago or so... my estimate is: +60'', +8''
    3:07
    but these values will be re-evaluated at some point, better store them in a file or a variable that you can easily update when need be
    """
    offset_x, offset_y = y_srf + 60, -z_srf + 8
    #here it is sun center, so we should add the minus sign
    return (offset_x, offset_y)



def attach_aspect_solutions(start_unix, end_unix, config):
    """
    attach aspect solutions
    """
    #return
    mean_time=(start_unix+end_unix)/2.
    docs=mdb.get_aspect_solutions(start_unix, end_unix)
    first = None
    second = None
    for doc in docs:
        if doc['unix_time'] < mean_time:
            first = doc
        if doc['unix_time'] >= mean_time:
            second = doc
            break
    def interp(doc_a, doc_b, key, mean_time, elem=None):
        """
            linear interpolation between two consecutive  data plots
        """
        t0, t1 =doc_a['unix_time'], doc_b['unix_time']
        if t0==t1:
            if elem is None:
                return doc_a[key]
            else:
                return doc_a[key][elem]
        if elem is None:
            y0=doc_a[key]
            y1=doc_b[key]
        else:
            y0=doc_a[key][elem]
            y1=doc_b[key][elem]
        return (y0*(t1-mean_time)+y1*(mean_time-t0))/(t1-t0)

    if first is None or second is None:
        config['aux']['error']='Aspect solution not found'
    else:
        config['aux']['L0']=interp(first, second, 'solo_loc_carrington_lonlat', mean_time, elem=0)
        config['aux']['B0']=interp(first, second, 'solo_loc_carrington_lonlat', mean_time, elem=1)
        config['aux']['sun_center']=get_sun_center_from_aspect(
                interp(first, second, 'y_srf', mean_time),
                interp(first, second, 'z_srf', mean_time))
        config['aux']['rsun']=interp(first, second, 'spice_disc_size', mean_time)
        config['aux']['roll']=interp(first, second, 'roll_angle_rpy', mean_time, elem=0)
        config['aux']['dsun']=interp(first, second, 'solo_loc_carrington_dist', mean_time, elem=0)*1000
        # distance to sun, in units of km
        config['aux']['data_source_file']=doc['filename']
        config['aux']['data_source_type']='Aspect'
        config['aux']['utc']=stu.unix2utc(mean_time)
        config['aux']['comments']=f'Interpolation of  two measurements at {stu.unix2utc(first["unix_time"])} and {stu.unix2utc(second["unix_time"])}' 



def create_imaging_tasks_for_BSD(bsd_id, 
                        min_counts=MIN_COUNTS,
                        duration=60,
                        energy_bands=[[4, 10], [16, 28]],
                        bkg_max_day_off=BKG_FILE_MAX_DAY_OFFSET, overwritten=False):
    """
    This function only register imaging tasks
    Parameters
    ----------
        min_counts: float 
            minimal counts per time bin
        min_duration: float
            minimal time per bin
        duration: float
            integration time in units of seconds
    """
    logger.info(f'\n\n\nCreating imaging tasks for BSD # {bsd_id}')
    exclude_imaging_existing= not overwritten

    docs,msg =mdb.find_bsd_for_imaging(bsd_id, exclude_imaging_existing)
    if not docs :
        logger.warning(f'No BSD doc or FITS found for BSD # {bsd_id}, due to   {msg} ')
        return
    for doc in docs:
        _create_imaging_tasks_for_BSD(doc, 
                        min_counts,
                        duration,
                        energy_bands,
                        bkg_max_day_off, overwritten)

def _create_imaging_tasks_for_BSD(doc, 
                        min_counts=MIN_COUNTS,
                        duration=60,
                        energy_bands=[[4, 10], [16, 28]],
                        bkg_max_day_off=BKG_FILE_MAX_DAY_OFFSET, overwritten=True):


    fits_doc=doc['bsd_fits']['fits']
    bsd_id=doc['bsd_fits']['_id']


    try:
        fname = os.path.join(fits_doc['path'], fits_doc['filename'])
    except (IndexError, KeyError):
        logger.warning(
            f'No signal data found for bsd {bsd_id} ')
        return
    if not os.path.isfile(fname):
        logger.warning(
            f'No signal data found for BSD {bsd_id}')
        return

    flare_doc=doc['flare']

    flare_id=flare_doc['flare_id']
    uid=doc['bsd_fits']['unique_id']

    flare_image_id = mdb.get_next_flare_image_id()

    
    flare_start_unix=max(flare_doc['start_unix'], fits_doc['data_start_unix']+1)
    flare_end_unix=min(flare_doc['end_unix'], fits_doc['data_end_unix']-1)

    logger.info(f"Flare time: {stu.unix2utc(flare_start_unix)}, duration {flare_end_unix-flare_start_unix}")
    


    try:

        bkg_meta= find_background(flare_start_unix, bkg_max_day_off=30)
    except IndexError:
        bkg_meta = None
    if not bkg_meta:
        logger.warning(
            f'No background data found for BSD {bsd_id}')
        return

    
    if flare_doc['peak_counts']> 3000:
        start, end= flare_doc['peak_unix_time']-30, flare_doc['peak_unix_time']+30
    else:
        start, end= fits_doc['data_start_unix']+TIME_MARGIN, fits_doc['data_end_unix']-TIME_MARGIN

    energies =[[4,10]] 

    try:
        if flare_doc['LC_statistics']['upper_ilc'] >1:
            energies.append([16,28])
    except (KeyError, TypeError):
        pass



    start=max(start, fits_doc['data_start_unix']+TIME_MARGIN)
    end=min(end, fits_doc['data_end_unix']-TIME_MARGIN)


    try:
        signal_utc = stu.unix2utc( (start+end)*0.5 )
        B0, L0, roll, rsun, dsun, solo_hee, solo_sun_r, sun_center = solo.SoloEphemeris.get_ephemeris_for_imaging(
            signal_utc)
    except Exception as e:
        logger.warning(f'No ephemeris data found for {bsd_id} (uid {uid}), error: {str(e)}')
        return


    boxes=[{'total_counts':  0,
                    'pixel_counts':  0,
                    'energy_range_keV': erange,
                    'unix_time_range': [start,  end],
                    'utc_range': [stu.unix2utc(start),  stu.unix2utc(end)]}
                    for erange in energies ]

    num_images = 0
    flare_image_ids = []
    for box in boxes:
        energy_range_str = f'{box["energy_range_keV"][0]}-{box["energy_range_keV"][1]}'
        logger.info(f'box: {box["energy_range_keV"]}')
        start_utc_str = stu.utc2filename(box['utc_range'][0])
        fits_prefix = f'stix_sci_preview_{bsd_id}_{flare_id}_uid_{uid}_{energy_range_str}keV_{start_utc_str}_{flare_image_id}'
        folder = os.path.join(quicklook_path, str(uid))
        if not os.path.exists(folder):
            os.makedirs(folder)
        num_images += 1
        task_id = uuid.uuid4().hex[0:10]
        spectral_model = 'auto'

        ospex_config = {
            'model': spectral_model
        }
        config = {
            'filename': fname,
            'bsd_id': bsd_id,
            'flare_id':flare_id,
            'unique_id': uid,
            'task_id': task_id,
            'author': 'bot',
            'hidden': False,
            'num_idl_calls': 0,
            'run_type': 'auto',
            'idl_status': '',
            'idl_message':'',#placeholder
            'aux': {
                'B0': B0.to(u.deg).value,
                'L0': L0.to(u.deg).value,
                'roll': roll.to(u.deg).value,
                'rsun': rsun.to(u.arcsec).value,
                'dsun': dsun.to(u.m).value,
                'sun_center': sun_center,
                'data_source_type': 'SPICE',
                'solo_sun_r': solo_sun_r.to(u.au).value,
                'solo_hee': solo_hee.to(u.km).value
            },
            'background': bkg_meta,
            'start_unix': box['unix_time_range'][0],
            'end_unix': box['unix_time_range'][1],
            'duration': box['unix_time_range'][1] - box['unix_time_range'][0],
            'pixel_counts': box['pixel_counts'],
            'energy_range': box["energy_range_keV"],  #energy in time 
            'utc_range': box['utc_range'],
            'total_counts': box['total_counts'],
            'idl_config': {
                'folder': folder,
                'prefix': fits_prefix,
                'fwdfit_shape':
                'ellipse' if box["energy_range_keV"][1] < 16 else 'multi',
                'ospex': ospex_config
            },
            'fits': {}, #placeholder
            'figs': {}, #placeholder
        }
        
        attach_aspect_solutions(box['unix_time_range'][0] - ASPECT_TIME_TOR,
                                box['unix_time_range'][1] + ASPECT_TIME_TOR,
                                config)

        imaging_inputs = bson.dict_to_json(config)
        imaging_inputs.update({
            'creation_time': datetime.now(), 
            '_id':  flare_image_id,
            'signal_data_type': 'PixelData'
            })
        logger.info(f"Inserting data into db for Flare #{flare_id} ")
        flare_image_ids.append(flare_image_id)
        mdb.insert_flare_image(imaging_inputs)
        flare_image_id += 1
    logger.info(f'{num_images} images will be created')
    if num_images > 0:
        if bsd_id >= 0:
            bsd_db.update_one({'_id': bsd_id},
                          {'$set': {
                              'image_ids': flare_image_ids
                          }},
                          upsert=False)
        flare_db.update_one({'_id': flare_doc['_id']},
                          {'$set': {
                              'image_ids': flare_image_ids
                          }})
        #push

def create_imaging_tasks_for_BSDs_between(start_id, end_id, overwritten=False):
    for i in range(start_id, end_id+1):
        create_imaging_tasks_for_BSD(i)

def create_imaging_tasks_for_bsd(start_id, end_id):
    for i in range(start_id, end_id+1):
        create_imaging_tasks_for_bsd(i)



def register_imaging_tasks_for_last_bsd(max_age= 18000 ):
    """
    called by L1 data processing pipeline
    imaging requests are submitted to database  after parsing the raw TM
    create recent flares
    """
    min_unix=stu.get_now()-max_age * 86400
    logger.info(f"registering imaging tasks for flares detected since{stu.unix2utc(min_unix)}")
    for _id in [d['_id'] for d in  bsd_db.find({'start_unix_time': {'$gte':min_unix},
        'image_ids.0':{'$exists':False} }, {'_id':1}).sort('_id', 1)]:
        create_imaging_tasks_for_BSD(_id)


def find_background(peak_time, bkg_max_day_off=15):
    bkg_fits_doc= mdb.find_best_background_fits(peak_time,
                           bkg_max_day_off,
                           level=DATA_LEVEL_FOR_IMAGING)
    if not bkg_fits_doc:
        print("Not found for :",stu.unix2utc(peak_time))
        return None

    bkg_fname = os.path.join(bkg_fits_doc['merg'][0]['path'],
                         bkg_fits_doc['merg'][0]['filename'])
    if not os.path.isfile(bkg_fname):
        logger.warning(
            f'{bkg_fname} not exists!')
        return None
    background = {
                'filename': bkg_fname,
                'req_form_id': bkg_fits_doc['_id'],
                'fits_id': bkg_fits_doc['merg'][0]['_id'],
                'unique_id': bkg_fits_doc['merg'][0]['request_id'],
            }

    return background





def update_auxiliary_data():
    db_flare_images = mdb.get_collection('flare_images')
    for doc in db_flare_images.find():
        attach_aspect_solutions(doc['start_unix'] - ASPECT_TIME_TOR,
                                doc['end_unix'] + ASPECT_TIME_TOR, doc)
        doc['num_idl_calls'] = 0
        #reset
        logger.info(f'updating {doc["_id"]}...')

        db_flare_images.update_one({'_id': doc['_id']}, {'$set': doc})


def main():
    if len(sys.argv) == 1:
        print('Usage:')
        print('imaging_take_manager ')
        register_imaging_tasks_for_last_bsd()

    elif len(sys.argv)==2:
        create_imaging_tasks_for_BSD(int(sys.argv[1]) )
    elif len(sys.argv)==3:
        create_imaging_tasks_for_BSDs_between(int(sys.argv[1]), 
                int(sys.argv[2]) )
    #else:
    #update_ephmeris()
if __name__ == '__main__':
    main()
