#!/usr/bin/python3
'''
L1 flare data analysis pipeline
  what it does here includes
  - do background subtraction for flare data
  - find flare location solution
  - pre-requirements:
     flare data should be processed by flare_goes_class.py and flare_detection.py and sci_packets_analyzer.py 
    1) flare_goes_class prepare ephemeris for flare data
    2) sci_packets_analyzer prepares spectra and  check if it is background
    3) flare detection tells flare start time, end time, flare peak counts
Author: Hualin Xiao, Aug. 2021
email: hualin.xiao@fhnw.ch
'''
import sys
import os
import json
import numpy as np
from datetime import datetime
from stix.core import stix_datatypes as sdt
from stix.spice import stix_datetime
from stix.core import mongo_db as db
from stix.core import stix_logger
from stix.core import config
from stix.flare_pipeline import flare_location_fitter as flf
from pprint import pprint

mdb = db.MongoDB()
logger = stix_logger.get_logger()
DET_ID_START=20
DET_ID_END=32
DET_SENSITVE_AREA=80.96
GRID_OPEN_AREA_RATIO=np.array([0.292312,
		0.272487,
		0.252282,
		0.257189,
		0.281502,
		0.27254,
		0.292103,
		0.260207,
		0.335213,
		0.305518,
		0.335412,
		0.254877,
		0.264532,
		0.255886,
		0.304781,
		0.335199,
		0.30505,
		0.258,
		0.257187,
		0.252324,
		0.281503,
		0.260178,
		0.280854,
		0.257033,
		0.264449,
		0.261743,
		0.292299,
		0.272491,
		0.264514,
		0.254682])

BKG_MIN_DURATION=180
class SciL1Processor(object):
    def __init__(self):
        self.dt = []
        self.bsd_db= mdb.get_collection_bsd()

    def get_background_bsd(self, signal_unix_time, emin, emax, background_bsd_id=None):
        #find/load the most recent L1 data for background subtraction at a given signal time 
        if background_bsd_id is not None:
            return self.bsd_db.find_one({'_id':background_bsd_id})
        bsd_ealier=list(self.bsd_db.find({'synopsis.is_background':True,
            'synopsis.duration':{'$gt':BKG_MIN_DURATION},
            'synopsis.emin':{'$lte':emin},
            'synopsis.emax':{'$gte':emax},
            'start_unix':{'$lt':signal_unix_time}}).sort('start_unix',-1).limit(1))
        bsd_later=list(self.bsd_db.find({'synopsis.is_background':True, 
            'synopsis.duration':{'$gt':BKG_MIN_DURATION},
            'synopsis.emin':{'$lte':emin},
            'synopsis.emax':{'$gte':emax},
            'start_unix':{'$gt':signal_unix_time}}).sort('start_unix',1).limit(1))
        try:
            t0=np.abs(bsd_ealier[0]['start_unix']-signal_unix_time)

        except (KeyError, IndexError):
            t0=np.inf
        try:
            t1=np.abs(bsd_later[0]['start_unix']-signal_unix_time)
        except (KeyError, IndexError):
            t1=np.inf
        closet= np.min(t1,t0)
        
        if closet==np.inf:
            return None
        elif closet==t0:
            return bsd_ealier[0]
        else:
            return bsd_later[0]
        
    
    def get_flare_spectra_and_location(self, bsd_doc):
        results={}
        _id=bsd_doc['_id']
        try:
            flare_ids=bsd_doc['synopsis']['flares']['flare_ids']
        except KeyError:
            print(f'no flare information found in BSD #{_id}!')
            return None

        if len(flare_ids)>1:
            print(f'BSD {_id} contains several flares,but only the first will be processed!')
        flare_id=flare_ids[0]
        signal_counts=np.array(bsd_doc['synopsis']['flares']['flare_pixel_counts']['counts']) #already energy integrated 1d 384 vector
        signal_emin=bsd_doc['synopsis']['flares']['flare_pixel_counts']['emin']
        signal_emax=bsd_doc['synopsis']['flares']['flare_pixel_counts']['emax']
        signal_duration=bsd_doc['synopsis']['flares']['integration_times'][0]
        bkg_doc=get_background_bsd(bsd_doc['start_unix'])
        if not bkg_doc:
            print(f'No background found for {_id}')
            return


        bkg_rates_full_energy=np.array(bkg_doc['synopsis']['pixel_rate_spec']) #rate are not energy integrated, dimemsions 384x32
        bkg_duration=bkg_doc['synopsis']['duration']

        bkg_rates=np.sum(bkg_rates_full_energy[:,emin:emax], axis=0) #only select bins in the energy range, energy integrated counts
        bkg_rate_error=np.sqrt(bkg_rates*bkg_duration) #only statistic error 

        bkg_cfl_rate=bkg_rates[8*12:9*12]
        bkg_cfl_rate_error=bkg_rate_error[8*12:9*12]
        cfl_pixel_counts=signal_counts[8*12:9*12] #already energy bin integrated

        cfl_pixel_bkgsub_counts=cfl_pixel_counts - signal_duration * bkg_cfl_rate
        cfl_pixel_bkgsub_count_error=np.sqrt(cfl_pixel_counts  + (signal_duration*bkg_cfl_rate_error) **2)


        sig_bkgsub_counts=signal_counts - signal_duration*bkg_rates #background subtracted counts
        sig_bkgsub_counts_error=np.sqrt(signal_counts+(signal_duration*bkg_rate_error)**2)

        det_bkgsub_fluence=np.array([np.sum(sig_bkgsub_counts[idet*12:(idet*12+12)])/(GRID_OPEN_AREA_RATIO[idet]*DET_SENSITVE_AREA)
            for idet in range(32)]) 
        #background subtracted and energy integrated counts, normalized by open area ratio

        det_bkgsub_fluence_error=np.array([np.sqrt(np.sum(sig_bkgsub_counts[idet*12:(idet*12+12)]))/(GRID_OPEN_AREA_RATIO[idet]*DET_SENSITVE_AREA)
            for idet in range(32)]) 
        #detector counts errors, only consider statistical uncertainties

        mean_fluence=np.mean(det_bkgsub_fluence[DET_ID_START: DET_ID_END])
        mean_fluence_error=np.sqrt(np.sum([det_bkgsub_fluence_error[i]**2 
            for i in range(DET_ID_START, DET_ID_END)]))
        flare_utc=sdt.unix2utc(bsd_doc['start_unix'])

        res=flf.fit_location(sig_bkgsub_counts, sig_bkgsub_counts_error, 
                mean_fluence, mean_fluence_error, flare_utc,  use_small_pixels=True, use_det_fluence=True)
        data={'flare_location':
                {
                    'solution':res,
                    'inputs':{
            
                    'detector_mean_fluence_units':'counts / mm2',
                    'detector_mean_fluence':mean_fluence, 
                    'detector_mean_fluence_error':mean_fluence_error,
                    'sig_bkgsub_counts':sig_bkgsub_counts,
                    'sig_bkgsub_counts_error':sig_bkgsub_counts_error,
                    'background_bsd_id':bkg_doc['_id'],
                    }
                }
                }
        pprint(data)


    def process_L1_BSD_in_file(self, file_id):
        collection = mdb.get_collection_bsd()
        bsd_cursor = collection.find({'run_id': file_id, 'SPID':54115}).sort('_id', 1)
        for doc in bsd_cursor:
            self.get_flare_spectra_and_location(doc)
            break





if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('flare_location_solver run_id')
        print('flare_location_solver run_id_start id_end')
    elif len(sys.argv)==2:
        process_L1_BSD_in_file(int(sys.argv[1]))
    else:
        for i in range(int(sys.argv[1]),int(sys.argv[2])+1):
            process_L1_BSD_in_file(i)

