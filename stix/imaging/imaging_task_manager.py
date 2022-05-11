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
from stix.flare_pipeline import image_viewer as imv

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

last_bkg_fits_doc={}

def register_imaging_task_for_science_data(bsd_ids=[]):
    if not bsd_ids:
        #only processing for data which are not background
        docs = bsd_db.find({
            'name': 'L1',
            'synopsis.is_background': False,
            'flare_image_ids': {
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
        interval=600,
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
        interval: float
            time interval between two selected time periods
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
    sun_center=ephem['sun_center']

    #find flares in the time frame
    bsd_flare_time_ranges = [[
        x['start_unix'] if x['start_unix'] > bsd_start_unix else
        bsd_start_unix, x['peak_unix_time']
        if bsd_start_unix <= x['peak_unix_time'] <= bsd_end_unix else None,
        x['end_unix'] if x['end_unix'] < bsd_end_unix else bsd_end_unix
        ] for x in mdb.find_flares_by_time_range(bsd_start_unix, bsd_end_unix)]

    if not bsd_flare_time_ranges:
        logger.warning(f'No flares found for {bsd_id} (uid {uid})')
        return
    try:
        l1 = ScienceL1.from_fits(fname)
    except IndexError:
        logger.warning(f'Failed to read fits file for BSD # {bsd_id} (uid {uid})')
        return

    bkg_fname = os.path.join(bkg_fits['merg'][0]['path'],
            bkg_fits['merg'][0]['filename'])
    if not bsd_flare_time_ranges:
        logger.warning(f'No flares found for {bsd_id} (uid {uid})')
        return
    #(bsd_flare_time_ranges, energy_bands)
    boxes = l1.get_time_ranges_for_imaging(energy_bands,
            bsd_flare_time_ranges, min_counts, duration, interval)
    if boxes is None:
        logger.warning(f'No boxes selected for BSD {bsd_id} (uid {uid})')
        return
    num_images = 0
    flare_image_ids=[]

    for box in boxes:
        for ie, energy in enumerate(energy_bands):
            energy_range_str='-'.join([str(x) for x in energy])
            start_utc_str=stu.utc2filename(box['utc_range'][0])
            fits_prefix = f'stix_ql_image_sci_{bsd_id}_uid_{uid}_{energy_range_str}keV_{start_utc_str}_{flare_image_id}'

            folder=os.path.join(quicklook_path, str(uid))
            if not os.path.exists(folder):
                os.makedirs(folder)
            if box['counts_enough'][ie]:
                task_id= uuid.uuid4().hex[0:10]
                num_images += 1
                config={
                    'filename': fname,
                    'bsd_id': doc['_id'],
                    'unique_id': uid,
                    'task_id':task_id,
                    'author':'bot',
                    'hidden':False,
                    'idl_args':idl_args, 
                    'num_idl_calls':0,
                    'run_type':'auto',
                    'idl_status':'',
                    'aux': {
                        'B0': B0.to(u.deg).value,
                        'L0': L0.to(u.deg).value,
                        'roll': roll.to(u.deg).value,
                        'rsun': rsun.to(u.arcsec).value,
                        'dsun': dsun.to(u.m).value,
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
                    'energy_range': energy, #energy in time 
                    'utc_range':box['utc_range'],
                    'total_counts':box['total_counts'][ie],
                    'idl_config':{'folder':folder,  'prefix':fits_prefix, 'fwdfit_shape':'ellipse' if energy[1] <15 else 'multi'},
                    'fits':{},
                    'figs':{}
                    }

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




if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Usage:')
        print('flare_image_creator <regtask | runidl | tosvg> [ids]')
    elif len(sys.argv)==2:
        register_imaging_task_for_science_data()
    else:
        ids = [int(i) for i in sys.argv[2:]]
        register_imaging_task_for_science_data(ids)


