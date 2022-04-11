"""
 Generate inputs for imaging software
"""
import os
import sys
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
from stix.core import config
from stix.core import mongo_db as db
from stix.analysis.science_l1 import  ScienceL1 

from stix.core import logger
from stix.spice import solo
from stix.utils import energy_bins as eb
from stix.utils import bson 
from stix.spice import time_utils as stu
#TO install stixdcpy: pip install git+https://github.com/drhlxiao/stixdcpy.git
##### Constants
mdb = db.MongoDB()
logger = logger.get_logger()
quicklook_path = config.get_config(
    'pipeline.daemon.flare_images')
bsd_db=mdb.get_collection('bsd')

req_db=mdb.get_collection('data_requests')


from pprint import pprint


def generate_inputs_for_science_data(bsd_ids=[]):
    if not bsd_ids:
        docs=bsd_db.find({'name':'L1', 'synopsis.is_background':False,'imaging':{'$exists':False}}).sort('_id', -1)
    else:
        docs=bsd_db.find({'name':'L1', 'synopsis.is_background':False,'_id':{'$in':bsd_ids}}).sort('_id', -1)

    for doc in docs:
        create_images_for_bsd(doc)
        
        
def create_images_for_bsd(doc, min_counts=5000,
        min_duration=30,
                          imaging_energies=[[4,10],[16,28]], bkg_max_day_off=30,
                          overlap_time=0.5 ):
    """
    bulk science data document
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

    bsd_flare_time_ranges=[[x['start_unix'], x['end_unix']] for x in  mdb.find_flares_by_time_range(bsd_start_unix, bsd_end_unix)]
    if not bsd_flare_time_ranges:
        logger.warning(f'No flares found for {bsd_id} (uid {uid})')
        return
    bsd_flare_time_ranges=np.array(bsd_flare_time_ranges)
    signal_time_range=[max(np.min(bsd_flare_time_ranges), bsd_start_unix),
                        min(np.min(bsd_flare_time_ranges), bsd_end_unix)]

    #signal utc
    fits_doc = mdb.get_bsd_fits_info_by_request_id(uid)           
    if not fits_doc:
        logger.warning(f'No signal Fits file found for {bsd_id} (uid {uid})')
        return
    try:    
        signal_utc=stu.unix2utc(bsd_flare_time_ranges[0][0])
        B0,L0, roll, rsun=solo.SoloEphemeris.get_B0_L0_roll_radius(signal_utc)
    except ValueError:
        raise
        #print('No aux data for ',uid)
        return
    
    logger.info(f'finding fits file for request {uid}')
    fname=os.path.join(fits_doc[0]['path'], fits_doc[0]['filename'])
    l1=ScienceL1.from_fits(fname)
    #tbins,cnts=l1.get_time_bins_for_imaging(bsd_flare_times)

    boxes_energy_low=np.min(np.array(imaging_energies))
    boxes_energy_high=np.max(np.array(imaging_energies))
    box_emax_sci, box_emax_sci=eb.keV2sci(boxes_energy_low, boxes_energy_high)
    
    bkg_fits_docs = list(mdb.find_L1_background(bsd_start_unix, bkg_max_day_off, 
                                                emin=box_emax_sci, emax=box_emax_sci))
    
    if not bkg_fits_docs:
        #print(f'No background file find for {doc["_id"]}')
        logger.warning(f'No background Fits file found for {bsd_id} (uid {uid})')
        return

    bkg_fits=bkg_fits_docs[0] # select the most recent one
    #print("BKG:", bkg_fits_docs)
    #pprint(bkg_fits)
    bkg_fname=os.path.join(bkg_fits['merg'][0]['path'], bkg_fits['merg'][0]['filename'])
    l1=ScienceL1.from_fits(fname)
    box_energies=[]
    box_sci_range=[]
    box_time_ranges=[]
    out_filename=[]
    num_boxes=0
    ibox=0
    
    
    for e in imaging_energies:
        tbins,cnts=l1.get_time_bins_for_imaging(e[0], e[1], min_counts, min_duration , signal_time_range)
        if tbins is None:
            logger.error('Failed to compute time bins')
            continue
        for tb in tbins.tolist():
        
            box_energies.append(e)
            box_sci_range.append(eb.keV2sci(float(e[0]),float(e[1])))
            box_time_ranges.append([stu.unix2utc(tb[0]-overlap_time), stu.unix2utc(tb[1]+overlap_time)] )
            fits_prefix=f'sci_{bsd_id}_uid_{uid}_{ibox}'
            out_fits_fnames=[fits_prefix+ext for ext in ['_fwfit.fits', '_bp.fits']]
            ibox+=1
            
                         
    
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
        'boxes':{
                      'time_range_unixs':tbins,
                'time_range_utcs':box_time_ranges,
                'counts':cnts,
               'energies_keV':box_energies,
                'box_energies_sci':box_sci_range,
        },
 
    })
    logger.info(f"Inserting data into db for bsd #{bsd_id}")
    bsd_db.update_one({'_id':bsd_id}, {'$set':{'imaging': imaging_inputs, 'start_utc':stu.unix2utc(start_unix)}}, upsert=False)
if __name__=='__main__':
    if len(sys.argv)==1:
        generate_inputs_for_science_data()   
    else:
        ids=[int(i) for i in sys.argv[1:]]
        generate_inputs_for_science_data(ids)   



# In[ ]:




