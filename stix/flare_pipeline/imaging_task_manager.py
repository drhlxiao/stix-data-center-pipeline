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
from stix.utils import energy_bins as eb
from stix.utils import bson
from stix.spice import time_utils as stu
from pathlib import Path

##### Constants
mdb = db.MongoDB()
logger = logger.get_logger()
quicklook_path = config.get_config('pipeline.daemon.flare_images')
bsd_db = mdb.get_collection('bsd')
flare_images_db= mdb.get_collection('flare_images')
req_db = mdb.get_collection('data_requests')
MIN_COUNTS=10000

ASPECT_TIME_TOR=300 #time tolerance for looking for aspect solution 

last_bkg_fits_doc={}
def get_sun_center_from_aspect(y_srf,  z_srf):
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
    offset_x, offset_y=y_srf + 60, -z_srf +8
    #here it is sun center, so we should add the minus sign
    return (-offset_x, -offset_y)


def attach_aspect_solutions(start_unix, end_unix, config):
    """
    attach aspect solutions
    """
    #return
    mean_time=(start_unix+end_unix)/2.
    docs=mdb.get_aspect_solutions(start_unix, end_unix)
    min_dt=np.inf
    for doc in docs:
        try:
            if np.abs(doc['unix_time']-mean_time)<min_dt:
                config['aux']['L0'], config['aux']['B0']=doc['solo_loc_carrington_lonlat']
                #config['aux']['sun_center']=(-doc['y_srf'],doc['z_srf'])
                config['aux']['sun_center']=get_sun_center_from_aspect(doc['y_srf'], doc['z_srf'])
                #Not sure about the sign, need to be confirmed
                config['aux']['rsun']=doc['spice_disc_size']
                config['aux']['roll']=doc['roll_angle_rpy'][0]
                config['aux']['dsun']=doc['solo_loc_carrington_dist'][0]*1000 # distance to sun, in units of km
                config['aux']['data_source_file']=doc['filename']
                config['aux']['data_source_type']='Aspect'
                config['aux']['utc']=stu.unix2utc(doc['unix_time'])
                min_dt=doc['unix_time']
        except (KeyError, IndexError):
            pass 
def register_imaging_task_for_science_data(bsd_ids=[]):
    if not bsd_ids:
        #only processing for data which are not background
        docs = bsd_db.find({
            'name': 'L1',
            'synopsis.is_background': False,
            'flare_image_ids.0': {
                '$exists': False
                }
            }).sort('_id', -1)
    else:
        docs = bsd_db.find({
            'name': 'L1',
            'synopsis.is_background': False,
            '_id': {
                '$in': bsd_ids
                }
            }).sort('_id', -1)

    for doc in docs:
        logger.info(f'Adding imaging tasks for BSD {doc["_id"]}')
        queue_imaging_tasks(doc)

def queue_imaging_tasks(doc,
        min_counts=MIN_COUNTS,
        duration=60,
        energy_bands=[[4, 10], [16, 28]],
        bkg_max_day_off=60):
    """
    This function only pushes parameters for imaging to  the database
    Parameters
    ----------
        min_counts: float 
            minimal counts per time bin
        min_duration: float
            minimal time per bin
        duration: float
            integration time in units of seconds
    """
    if doc['name']!='L1':
        logger.info(f"{doc['_id']} is not L1 science data")
        return
    uid = int(doc.get('unique_id', -1))
    bsd_id = doc['_id']

    flare_image_id=mdb.get_next_flare_image_id()

    bsd_start_unix = doc['start_unix']
    bsd_end_unix = doc['end_unix']

    fits_doc = mdb.get_bsd_fits_info_by_request_id(uid)
    if not fits_doc:
        logger.warning(f'No signal Fits file found for {bsd_id} (uid {uid})')
        return
    boxes_energy_low = np.min(np.array(energy_bands))
    boxes_energy_high = np.max(np.array(energy_bands))

    box_emin_sci, box_emax_sci = eb.keV2sci(boxes_energy_low,
            boxes_energy_high)
    #energy range

    bkg_fits_docs = list(
            mdb.find_L1_background(bsd_start_unix,
                bkg_max_day_off,
                emin=box_emin_sci,
                emax=box_emax_sci))
    #find background data acquired within bkg_max_day_off days
    if not bkg_fits_docs:
        logger.warning(
                f'No background Fits file found for {bsd_id} (uid {uid}), trying to use last one')
        bkg_fits_docs=last_bkg_fits_doc
    try:
        bkg_fits = bkg_fits_docs[0]  # select the most recent one
        fname = os.path.join(fits_doc[0]['path'], fits_doc[0]['filename'])
    except (IndexError,KeyError):
        logger.warning(f'No background data found for BSD {bsd_id} (uid {uid})')
        return
    if not os.path.isfile(fname):
        logger.warning(f'No background data found for BSD {bsd_id} (uid {uid})')
        return

    last_bkg_fits_doc[0]=bkg_fits
    try:
        signal_utc = stu.unix2utc((bsd_start_unix + bsd_end_unix) / 2.)
        B0, L0, roll, rsun, dsun, solo_hee, solo_sun_r, sun_center= solo.SoloEphemeris.get_ephemeris_for_imaging(
                signal_utc)
    except ValueError:
        logger.warning(f'No ephemeris data found for {bsd_id} (uid {uid})')
        return
    #sun_center=ephem['sun_center']

    #find flares in the time frame
    try:
        l1 = ScienceL1.from_fits(fname)
    except IndexError:
        logger.warning(f'Failed to read fits file for BSD # {bsd_id} (uid {uid})')
        return
    bkg_fname = os.path.join(bkg_fits['merg'][0]['path'],
            bkg_fits['merg'][0]['filename'])


    bsd_flare_time_ranges =[ {
        'start': max(x['start_unix'],  bsd_start_unix),  
        'peak':  x['peak_unix_time']  if bsd_start_unix <= x['peak_unix_time'] <= bsd_end_unix else None,
        'end':  min(x['end_unix'],  bsd_end_unix)
            } for x in mdb.find_flares_by_time_range(bsd_start_unix, bsd_end_unix)]

    boxes = l1.get_time_ranges_for_imaging(bsd_flare_time_ranges, energy_bands, min_counts, duration )
    #two energies at the peak

    if boxes is None:
        logger.warning(f'No boxes selected for BSD {bsd_id} (uid {uid})')
        return
    num_images = 0
    flare_image_ids=[]

    for box in boxes:
        energy_range_str=f'{box["energy_range_keV"][0]}-{box["energy_range_keV"][1]}'
        start_utc_str=stu.utc2filename(box['utc_range'][0])
        fits_prefix = f'stix_imaging_spectroscopy_sci_{bsd_id}_uid_{uid}_{energy_range_str}keV_{start_utc_str}_{flare_image_id}'
        folder=os.path.join(quicklook_path, str(uid))
        if not os.path.exists(folder):
            os.makedirs(folder)
        num_images += 1
        task_id= uuid.uuid4().hex[0:10]
        
        config={
                    'filename': fname,
                    'bsd_id': doc['_id'],
                    'raw_file_id':doc['run_id'],
                    'unique_id': uid,
                    'task_id':task_id,
                    'author':'bot',
                    'hidden':False,
                    'num_idl_calls':0,
                    'run_type':'auto',
                    'idl_status':'',
                    'aux': {
                        'B0': B0.to(u.deg).value,
                        'L0': L0.to(u.deg).value,
                        'roll': roll.to(u.deg).value,
                        'rsun': rsun.to(u.arcsec).value,
                        'dsun': dsun.to(u.m).value,
                        'sun_center':sun_center,
                        'data_source_type':'SPICE',
                        'solo_sun_r':solo_sun_r.to(u.au).value,
                        'solo_hee':solo_hee.to(u.km).value
                        },
                    'background': {
                        'filename': bkg_fname,
                        'req_form_id': bkg_fits['_id'],
                        'fits_id': bkg_fits['merg'][0]['_id'],
                        'unique_id': bkg_fits['merg'][0]['request_id'],
                        },
                    'start_unix': box['unix_time_range'][0],
                    'end_unix': box['unix_time_range'][1],
                    'duration':box['unix_time_range'][1]-box['unix_time_range'][0],
                    'pixel_counts': box['pixel_counts'],
                    'energy_range': box["energy_range_keV"], #energy in time 
                    'utc_range':box['utc_range'],
                    'total_counts':box['total_counts'],
                    'idl_config':{'folder':folder,  'prefix':fits_prefix, 'fwdfit_shape':'ellipse' if box["energy_range_keV"][1] <16 else 'multi'},
                    'fits':{},
                    'figs':{}
                }

        attach_aspect_solutions(box['unix_time_range'][0]-ASPECT_TIME_TOR, box['unix_time_range'][1]+ASPECT_TIME_TOR, config)

        imaging_inputs = bson.dict_to_json(config)
        imaging_inputs['creation_time']=datetime.now()
        logger.info(f"Inserting data into db for bsd #{bsd_id} ")
        imaging_inputs['_id']=flare_image_id
        flare_image_ids.append(flare_image_id)
        mdb.insert_flare_image(imaging_inputs)
        flare_image_id+=1
    logger.info(f'{num_images} images will be created')
    if num_images > 0:
        bsd_db.update_one({'_id':doc['_id']},{'$set':{'flare_image_ids': flare_image_ids}}, upsert=False)
        #push 

def register_imaging_tasks_for_file(file_id):
    """
    called by L1 data processing pipeline
    imaging requests are submitted to database  after parsing the raw TM
    """
    docs=bsd_db.find({'run_id':file_id, 'name':'L1'})
    ids=[int(x['_id']) for x in docs]
    register_imaging_task_for_science_data(ids)

def update_auxiliary_data():
    db_flare_images=mdb.get_collection('flare_images')
    for doc in db_flare_images.find():
        attach_aspect_solutions(doc['start_unix']-ASPECT_TIME_TOR, doc['end_unix']+ASPECT_TIME_TOR, doc)
        doc['num_idl_calls']=0
        logger.info(f'updating {doc["_id"]}...')
        db_flare_images.update_one({'_id':doc['_id']}, {'$set':doc})





if __name__ == '__main__':
    update_auxiliary_data()
    """
    if len(sys.argv) == 1:
        print('Usage:')
        print('imaging_take_manager regtask [ids]')
    elif len(sys.argv)==2:
        register_imaging_task_for_science_data()
    else:
        ids = [int(i) for i in sys.argv[2:]]
        register_imaging_task_for_science_data(ids)
    """


