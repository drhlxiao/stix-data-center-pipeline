"""
 Generate inputs for imaging software
"""
import os
import sys
import json
import numpy as np
import subprocess
from datetime import datetime, timedelta

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits

from sunpy import map
from sunpy.map import make_fitswcs_header
from sunpy.coordinates.frames import HeliocentricEarthEcliptic,HeliographicStonyhurst
from stix.core import config
from stix.core import mongo_db as db
from stix.analysis.science_l1 import  ScienceL1 

from stix.core import logger
from stix.spice import solo
from stix.utils import energy_bins as eb
from stix.utils import bson 
from stix.spice import time_utils as stu
from pathlib import Path


#TO install stixdcpy: pip install git+https://github.com/drhlxiao/stixdcpy.git
##### Constants
mdb = db.MongoDB()
logger = logger.get_logger()
quicklook_path = config.get_config(
    'pipeline.daemon.flare_images')
bsd_db=mdb.get_collection('bsd')

req_db=mdb.get_collection('data_requests')
SSW_HOME='/data2/ssw'
IDL_HOME='/opt/idl88/idl88'
IDL_SCRIPT_PATH=Path(__file__).resolve().parent.parent/ 'idl'


from pprint import pprint
HOST='https://datacenter.stix.i4ds.net/'

def execute_script(shell_filename, verbose=False):
    """
    execute script
    """
    subprocess.call(['chmod', 'u+x', shell_filename])
    cmd_output = subprocess.run([shell_filename], shell=True, stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE)
    check_for_errors(cmd_output, verbose)


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
        print(f'{stderr}\n{stdout}')

def generate_inputs_for_science_data(bsd_ids=[]):
    if not bsd_ids:
        docs=bsd_db.find({'name':'L1', 'synopsis.is_background':False,'imaging':{'$exists':False}}).sort('_id', -1)
    else:
        docs=bsd_db.find({'name':'L1', 'synopsis.is_background':False,'_id':{'$in':bsd_ids}}).sort('_id', -1)

    for doc in docs:
        generate_imaging_inputs(doc)
        
def call_idl(inputs, bkg_fits, sig_fits):
    parameters=','.join([f"'{x}'" if isinstance(x, str) else str(x) for x in inputs])
    bkg_filename=os.path.basename(bkg_fits)
    sig_filename=os.path.basename(sig_fits)
    script_lines=["!PATH=!PATH",
            f'; The two fits files can be downloaded via the links below:',
            f';{HOST}/download/fits/filename/{bkg_filename}',
            f';{HOST}/download/fits/filename/{sig_filename}',
                f'.run {IDL_SCRIPT_PATH}/stix_image_reconstruction.pro',
                f'stx_image_reconstruct, {parameters}',
                'exit']
    sc_fname=os.path.join(IDL_SCRIPT_PATH, 'top.pro')
    f=open(sc_fname, 'w')
    for l in script_lines:
        f.write(l+'\n')
    print('updated', sc_fname)
    f.close()
    print("executing script")
    try:
        execute_script(IDL_SCRIPT_PATH/'stix_imaging.sh')
    except RuntimeError:
        print("Run time error!")
        return False
    
    print("success!")
    return True


        


def generate_imaging_inputs(doc, min_counts=5000,
        min_duration=30,
                          imaging_energies=[[4,10],[16,28]], bkg_max_day_off=30,
                          overlap_time=0.5 ):
    """
    min_counts: minimal counts per time bin
    min_duration: minimal time per bin
    """
    uid=doc['unique_id']
    bsd_id=doc['_id']
    #print(uid)
    uid=int(uid)

    #find flare times

    bsd_start_unix=doc['start_unix']
    bsd_end_unix=doc['end_unix']

    #signal utc
    fits_doc = mdb.get_bsd_fits_info_by_request_id(uid)           
    if not fits_doc:
        logger.warning(f'No signal Fits file found for {bsd_id} (uid {uid})')
        return
    boxes_energy_low=np.min(np.array(imaging_energies))
    boxes_energy_high=np.max(np.array(imaging_energies))
    box_emax_sci, box_emax_sci=eb.keV2sci(boxes_energy_low, boxes_energy_high)
    #energy range

    bkg_fits_docs = list(mdb.find_L1_background(bsd_start_unix, bkg_max_day_off, 
                                                emin=box_emax_sci, emax=box_emax_sci))
    #find background data acquired within bkg_max_day_off days
    if not bkg_fits_docs:
        logger.warning(f'No background Fits file found for {bsd_id} (uid {uid})')
        return

    try:    
        signal_utc=stu.unix2utc((bsd_start_unix+bsd_end_unix)/2.)
        B0,L0, roll, rsun=solo.SoloEphemeris.get_B0_L0_roll_radius(signal_utc)
    except ValueError:
        logger.warning(f'No ephemeris data found for {bsd_id} (uid {uid})')
        return
    
    bsd_flare_time_ranges=[[x['start_unix'] if x['start_unix'] > bsd_start_unix else bsd_start_unix, 
        x['peak_unix_time'] if bsd_start_unix <= x['peak_unix_time']<=bsd_end_unix else None ,
        x['end_unix'] if x['end_unix'] < bsd_end_unix else bsd_end_unix ] for x in  mdb.find_flares_by_time_range(bsd_start_unix, bsd_end_unix)]
    if not bsd_flare_time_ranges:
        logger.warning(f'No flares found for {bsd_id} (uid {uid})')
        return
    bkg_fits=bkg_fits_docs[0] # select the most recent one
    fname=os.path.join(fits_doc[0]['path'], fits_doc[0]['filename'])
    #signal filename 

    l1=ScienceL1.from_fits(fname)
    bkg_fname=os.path.join(bkg_fits['merg'][0]['path'], 
            bkg_fits['merg'][0]['filename'])
    box_energies=[]
    box_sci_range=[]
    box_time_ranges=[]
    out_filename=[]
    num_boxes=0
    ibox=0

    #determine time ranges

        

    if not bsd_flare_time_ranges:
        logger.warning(f'No flares found for {bsd_id} (uid {uid})')
        return 
    print(bsd_flare_time_ranges, imaging_energies)

    boxes = l1.get_time_ranges_for_imaging(imaging_energies, bsd_flare_time_ranges, min_counts, min_duration)
    print(boxes)
    if boxes is None:
        logger.warning(f'No time bins found for {bsd_id} (uid {uid})')
        return 
    num_images=0
    for tb in boxes:
        tb['fits']=[[]]*len(imaging_energies)
        for ie, energy in enumerate(imaging_energies):
            fits_prefix=f'sci_{bsd_id}_uid_{uid}_{ibox}_{ie}'
            out_fits_fnames=[os.path.join(quicklook_path, fits_prefix+ext) for ext in ['_fwfit.fits', '_bp.fits']]
            num_images+=1
            if tb['counts_enough'][ie]:
                success=call_idl([bkg_fname, fname, tb['utc_range'][0],tb['utc_range'][1],
                    tb['energy_range_sci'][ie][0],
                    tb['energy_range_sci'][ie][1],
                    out_fits_fnames[0],
                    out_fits_fnames[1],
                    round(L0,4), round(B0,4), round(rsun,4),round(roll,4)], bkg_fname, fname)
                if not success:
                    continue
                print("success, output:", out_fits_fnames)
                
                tb['fits'][ie]=out_fits_fnames

    if num_images==0:
        logger.warning(f'No image created for {bsd_id} (uid {uid}), counts too low')
        return

    imaging_inputs=bson.dict_to_json({
        'signal':{
            'filename':fname,
            'bsd_id':doc['_id'],
            'unique_id':uid,
            },
        'config':{
            'min_counts_per_bin':min_counts, 
            'path':quicklook_path,
      
        'fame_overlap_time':overlap_time,
            'min_duration':min_duration,
        },
        'aux':{
            'B0':B0,
            'L0':L0,
            'roll':roll,
            'rsun':rsun,
        },
     'background':{
            'filename':bkg_fname,
            'req_form_id':bkg_fits['_id'],
            'fits_id':bkg_fits['merg'][0]['_id'],
            'unique_id':bkg_fits['merg'][0]['request_id'],
            },
        'boxes': boxes
 
    })
    logger.info(f"Inserting data into db for bsd #{bsd_id}")
    bsd_db.update_one({'_id':bsd_id}, {'$set':{'imaging': imaging_inputs}}, upsert=False)
if __name__=='__main__':
    if len(sys.argv)==1:
        generate_inputs_for_science_data()   
    else:
        ids=[int(i) for i in sys.argv[1:]]
        generate_inputs_for_science_data(ids)   



# In[ ]:




