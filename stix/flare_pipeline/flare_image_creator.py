"""
 create flare images for science data
 Procedure to create image for science data
 Step 1:
    prepare inputs for imaging, store the data into flare_images database
 step 2:
    imaging daemon get inputs from the database can create inputs 
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
SSW_HOME = '/data2/ssw'
IDL_HOME = '/opt/idl88/idl88'
IDL_SCRIPT_PATH = '/data/scripts/imaging'

from pprint import pprint
HOST = 'https://datacenter.stix.i4ds.net/'


def execute_script(shell_filename, idl_script, verbose=False):
    """
    execute script
    """
    logger.info(f"Executing script {shell_filename} {idl_script}...")
    subprocess.call(['chmod', 'u+x', shell_filename])
    cmd_output = subprocess.check_call([shell_filename, idl_script])
    #        shell=True,
    #                            stderr=subprocess.PIPE,
    #                            stdout=subprocess.PIPE)
    #print('done, checking output')
    #check_for_errors(cmd_output, verbose)


def check_for_errors(output, verbose):
    """
    Check IDL output to try and decide if an error has occurred
    """
    stdout = output.stdout.decode('utf-8')
    stderr = output.stderr.decode('utf-8')
    # NOTE: For some reason, not only errors are output to stderr so we
    # have to check it for certain keywords to see if an error occurred
    if 'execution halted' in stderr.lower():
        raise RuntimeError(stderr)
    if 'failed to acquire license' in stderr.lower():
        raise RuntimeError(stderr)
    if verbose:
        logger.info(f'{stderr}\n{stdout}')


def register_imaging_task_for_science_data(bsd_ids=[]):
    if not bsd_ids:
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
        queue_imaging_tasks(doc)

def queue_imaging_tasks(doc,
        min_counts=2000,
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
                f'No background Fits file found for {bsd_id} (uid {uid})')
        return

    try:
        signal_utc = stu.unix2utc((bsd_start_unix + bsd_end_unix) / 2.)
        B0, L0, roll, rsun, solo_hee, solo_sun_r = solo.SoloEphemeris.get_ephemeris_for_imaging(
                signal_utc)
    except ValueError:
        logger.warning(f'No ephemeris data found for {bsd_id} (uid {uid})')
        return

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
        bkg_fits = bkg_fits_docs[0]  # select the most recent one
        fname = os.path.join(fits_doc[0]['path'], fits_doc[0]['filename'])
    except (IndexError,KeyError):
        logger.warning(f'No background data for {bsd_id} (uid {uid})')
        return
    try:
        l1 = ScienceL1.from_fits(fname)
    except IndexError:
        logger.warning(f'Failed to read fits file for BSD # {bsd_id} (uid {uid})')
        return

    bkg_fname = os.path.join(bkg_fits['merg'][0]['path'],
            bkg_fits['merg'][0]['filename'])
    ibox = 0
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
            fits_prefix = f'sci_{bsd_id}_uid_{uid}_{ibox}_{ie}_{num_images}_{box["utc_range"][0]}'
            output_filenames = [
                    os.path.join(quicklook_path, fits_prefix + ext)
                    for ext in ['_fwfit.fits', '_bp.fits', '.png']
                    ]
            if box['counts_enough'][ie]:
                idl_script_uid=uuid.uuid4().hex[0:10]
                idl_args=[[
                    bkg_fname, fname, box['utc_range'][0], box['utc_range'][1],
                    box['energy_range_keV'][ie][0],
                    box['energy_range_keV'][ie][1], 
                    output_filenames[0],
                    output_filenames[1],
                    round(L0.to(u.deg).value, 4),
                    round(B0.to(u.deg).value, 4),
                    round(rsun.to(u.arcsec).value, 4),
                    round(roll.to(u.deg).value, 4)
                    ], bkg_fname, fname, idl_script_uid]
                num_images += 1
                imaging_inputs = bson.dict_to_json({
                    'filename': fname,
                    'bsd_id': doc['_id'],
                    'unique_id': uid,
                    'idl_args':idl_args, 
                    'num_idl_calls':0,
                    'run_type':'auto',
                    'idl_status':'',
                    'aux': {
                        'B0': B0.to(u.deg).value,
                        'L0': L0.to(u.deg).value,
                        'roll': roll.to(u.deg).value,
                        'rsun': rsun.to(u.arcsec).value,
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
                    'energy_range': energy,
                    'utc_range':box['utc_range'],
                    'total_counts':box['total_counts'][ie],
                    'fits':output_filenames[0:2],
                    'figs':[]
                    })
                logger.info(f"Inserting data into db for bsd #{bsd_id}")
                imaging_inputs['_id']=flare_image_id
                flare_image_ids.append(flare_image_id)
                mdb.insert_flare_image(imaging_inputs)
                flare_image_id+=1
        if num_images > 0:
            bsd_db.update_one({'_id':doc['_id']},{'$set':{'flare_image_ids': flare_image_ids}}, upsert=False)
        #push 

        

def call_idl(inputs, bkg_fits, sig_fits, uuid='test'):
    """
    """
    parameters = ','.join(
            [f"'{x}'" if isinstance(x, str) else str(x) for x in inputs])
    bkg_filename = os.path.basename(bkg_fits)
    sig_filename = os.path.basename(sig_fits)
    script_lines = [
            "!PATH=!PATH",
            f'; The two fits files can be downloaded via the links below:',
            f';{HOST}/download/fits/filename/{bkg_filename}',
            f';{HOST}/download/fits/filename/{sig_filename}',
            f'.run {IDL_SCRIPT_PATH}/stix_image_reconstruction.pro',
            f'stx_image_reconstruct, {parameters}', 'exit'
            ]
    sc_fname = os.path.join(IDL_SCRIPT_PATH, f'top_{uuid}.pro')
    f = open(sc_fname, 'w')
    for l in script_lines:
        f.write(l + '\n')
    f.close()
    try:
        execute_script(os.path.join(IDL_SCRIPT_PATH, 'stix_imaging.sh'), sc_fname)
    except RuntimeError:
        logger.error('IDL runtime error')
        return False

    return True


def create_images_in_queue():
    cursor=flare_images_db.find({'num_idl_calls':0, 'idl_status':''}).sort('_id',-1)
    create_images_for_bsd_docs(cursor)

def create_images_for_science_data(bsd_id):
    cursor=flare_images_db.find({'bsd_id':bsd_id})
    create_images_for_bsd_docs(cursor)

def create_qk_figures(doc, update_db=False, output_folder=None):
    """
    write images to the same folder if output_folder is None
    """
    try:
        bp_fname=doc['fits'][0]
        solo_hee=np.array(doc['aux']['solo_hee'])*u.km
        fnames=[doc['fits'][0], doc['fits'][1]]
        start_utc, end_utc=doc['utc_range']
        energy_range=doc['energy_range']
        figs=imv.images_to_graph(fnames[0], fnames[1], solo_hee, start_utc, end_utc, energy_range, output_folder)
    except Exception as e:
        logger.error(e)
        return None
    if update_db:
        if isinstance(figs,str):
            figs=[figs]
        updates={'$set':{'figs':figs}}
        flare_images_db.update_one({'_id':doc['_id']}, updates)
    return figs

def create_figures_ids_between(start_id, end_id):
    for i in range(start_id, end_id):
        doc=flare_images_db.find_one({'_id':i})
        if not doc:
            logger.warning("Failed to create figures for DocID:{i}")
            continue
        logger.info("Creating images for BSD#{doc['bsd_id']}, DocID:{i}")
        create_qk_figures(doc, update_db=True, output_folder=None)

def create_images_for_bsd_docs(cursor):
    for doc in cursor:
        args=doc.get('idl_args',None)
        if args is not None:
            logger.info(f"Processing {doc['_id']}")
            logger.info(f": prameters {str(args)}")
            success = call_idl(args[0], args[1], args[2], args[3])
            logger.info(f"End of processing {doc['_id']}, status: {success}")
            updates={'num_idl_calls':doc['num_idl_calls']+1, 'idl_status':success}
            if success:
                figs=create_qk_figures(doc, update_db=None, output_folder=None)
                if figs:
                    updates['figs']=figs

            updates={'$set':updates}
            flare_images_db.update_one({'_id':doc['_id']}, updates)
            

def register_imaging_tasks_for_file(file_id):
    docs=bsd_db.find({'run_id':file_id, 'name':'L1'})
    ids=[int(x['_id']) for x in docs]
    register_imaging_task_for_science_data(ids)




if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Usage:')
        print('flare_image_creator <regtask | runidl | tosvg> [ids]')
    elif sys.argv[1] == 'regtask':
        if len(sys.argv)==2:
            register_imaging_task_for_science_data()
        else:
            ids = [int(i) for i in sys.argv[2:]]
            register_imaging_task_for_science_data(ids)


