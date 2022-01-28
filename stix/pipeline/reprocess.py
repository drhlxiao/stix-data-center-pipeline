#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : reprocess_file.py
# @description  : reprocess file using the latest spice kernel data
# @author       : Hualin Xiao
# @date         : Feb. 11, 2020
#

import os
import sys
import glob
import time
from datetime import datetime
from stix.core import config
from stix.spice import datetime
from stix.core import mongo_db 
from stix.core import logger, idb, parser
from stix.fits import fits_creator
from  stix.analysis import calibration
from  stix.analysis import background_estimation as bkg
from  stix.analysis import flare_detection
from  stix.analysis import sci_packets_analyzer
from stix.core import spice_manager as spm
logger = logger.get_logger()

ENABLE_FITS_CREATION= True
DO_BULK_SCIENCE_DATA_MERGING= True
DO_FLARE_SEARCH=True
SCI_PACKET_SPIDS= ['54114', '54115', '54116', '54117', '54143', '54125']

daemon_config=config.get_config('pipeline.daemon')
MDB=mongo_db.MongoDB()
MAX_ALLOWED_DIFF=0.5
MAX_DELTA_TIME=24*3600

QL_SPIDS = [
     54118,
     54119,
     54120,
     54121,
     54122,
     54124,
     54125
]
SCI_SPIDS=[
    54114, #: 'L0',
    54115,#: 'L1',
    54116,#: 'L2',
    54117,#: 'L3',
    54143,#: 'L3',
    ]


def purge_fits(file_ids):
    collection_fits=MDB.get_collection('fits')
    for i in file_ids:
        print('Purging fits file created for file:',i)
        #MDB.get_collection('fits').delete_many(
        fits_docs=collection_fits.find({'file_id':i})
        for doc in fits_docs:
            filename=os.path.join(doc['path'],doc['filename'])
            try:
                print('removing: ',filename)
                os.remove(filename)
            except OSError as e:
                pass
        collection_fits.delete_many({'file_id': i})
def find_and_fix_raw_files(start_run, end_run=-1,  min_delta_time=0):
    if end_run==-1:
        end_run=start_run
    raw_db=MDB.get_collection('raw_files')
    docs=raw_db.find({'_id':{'$gte':start_run, '$lte':end_run},'hidden':False}).sort('_id',1)
    spice_fname=spm.spice.get_kernel_filename()
    print('Last spice kernel:',spice_fname)
    file_ids=[]
    for doc in docs:
        _id=doc['_id']
        file_ids.append(_id)
        start_unix=doc['data_start_unix_time']
        end_unix=doc['data_stop_unix_time']
        start_scet=doc.get('data_start_scet',0)
        end_scet=doc.get('data_end_scet',0)
        #if start_scet==0 or end_scet ==0 :
        #    print(f'Ignored: File # {_id}, no SCET found')
        #    continue
        if start_unix==0 or end_unix==0:
            print(f'Ignored: File # {_id}, Invalid unix')
            continue


        print(start_unix, start_scet)
        new_start_unix=datetime.scet2unix(start_scet)
        new_end_unix=datetime.scet2unix(end_scet)
        delta_start=abs(start_unix-new_start_unix)
        delta_end=abs(end_unix-new_end_unix)
        print('_id',', Delta Time:', delta_start, ', ' , delta_end)
        if delta_start < min_delta_time and  delta_end< min_delta_time:
            continue
        if delta_start > MAX_DELTA_TIME or delta_end >= MAX_DELTA_TIME:
            continue

        print(f'Updating file info for file: {_id}')

        updates={'$set':{
            'data_start_unix_time': new_start_unix,
            'data_stop_unix_time': new_end_unix,
            'spice_kernel': spice_fname,
            'update_time': datetime.get_now(),
            'version': doc.get('version',0)+1
            }
            }
        raw_db.update_one({'_id': _id}, updates, upsert=False)
    return file_ids

def fix_packets(run_id):
    print('Fixing packets for file:',run_id)
    pkt_db=MDB.get_collection('packets')
    packets=pkt_db.find({'run_id':run_id}).sort('_id',1)
    num=0
    errors=0
    for pkt in packets:
        header=pkt['header']
        _id=pkt['_id']
        try:
            updates={}
            if header.get('SCET', 0)>0:
                coarse=header['coarse_time']
                fine=header['fine_time']
                scet=header['SCET']
                utc=datetime.scet2utc(coarse,fine)
                unix=datetime.scet2unix(scet)
                updates['header.UTC']=utc
                updates['header.obt_utc']=utc
                updates['header.unix_time']=unix
            if header.get('SPID',0) in  QL_SPIDS:
                if pkt['parameters'][1][0]=='NIX00445':
                    scet_t0=pkt['parameters'][1][1]
                    utc=datetime.scet2utc(scet_t0)
                    updates['parameters.1.2']=utc
            if header.get('SPID',0) in  SCI_SPIDS:
                if pkt['parameters'][12][0]=='NIX00402':
                    scet_t0=pkt['parameters'][12][1]
                    utc=datetime.scet2utc(scet_t0)
                    updates['parameters.12.2']=utc
            if updates:
                pkt_db.update_one({'_id':_id},{'$set':updates}, upsert=False)
            num+=1
        except Exception as e:
            print(str(e))
            errors+=1


    print(f'{num} packets updated for File #{run_id}, Errors: {errors}')

def main(start_file, end_file, min_dt, force=False):
    spm.spice.refresh_kernels()
    if force:
        file_ids=range(start_file, end_file+1)
    else:
        file_ids=find_and_fix_raw_files(start_file, end_file, min_dt)
    for _id in file_ids:
        fix_packets(_id)
        try:
            print('Creating fits files for ',_id)
            fits_creator.create_fits(_id, daemon_config['fits_path'])
        except Exception as e:
            logger.error(str(e))
        try:
            print('creating L1 json file')
            sci_packets_merger.process(_id)
        except Exception as e:
            print(str(e))

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        start_id=int(sys.argv[1])
        end_id=start_id if len(sys.argv)==2 else int(sys.argv[2])
        min_dt=0 if len(sys.argv)<=3 else int(sys.argv[3])
        main(start_id, end_id, min_dt)
    else:
        print(' process start_id, [end_end], [min_dt]')

