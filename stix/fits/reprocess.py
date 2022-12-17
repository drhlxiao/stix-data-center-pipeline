'''
reprocess all fits files
'''
import sys
import os
import argparse
from collections import defaultdict
from datetime import datetime, timedelta
from itertools import chain
from pathlib import Path

from stix.core import datatypes as sdt
from stix.fits.io import hk_fits_writer as hkw
from stix.spice import time_utils as st
from stix.fits import fits_creator as fc
from stix.core import mongo_db, logger
logger = logger.get_logger()
import pandas as pd
db = mongo_db.MongoDB()
db = mongo_db.MongoDB()
#logger = get_logger(__name__)
bsd_db = db.get_collection('bsd')
dry_run=False
PATH='/data/fits'
fits_db= db.collection_fits
def remove_fits_except_science(start_unix, end_unix):
    curs=fits_db.find({'data_start_unix':{'$lte':end_unix}, 'data_end_unix':{'$gte' :start_unix}, 'product_group':{'$ne': 'science'} , 'level':'L1A' })
    for cur in curs:
        fits_filename = os.path.join(cur['path'], cur['filename'])
        logger.info(f'Removing older version of FITS file : {fits_filename}')
        if not dry_run:
            os.unlink(fits_filename)
            fits_db.delete_one({'_id':cur['_id']})

def remove_science_fits(request_id, start_unix, tpad = 300 ):
    curs=fits_db.find({'data_start_unix':{'$lte':start_unix-tpad}, 
        'data_start_unix':{'$gte' :start_unix + tpad}, 'product_group': 'science', 'request_id':request_id , 'level':'L1A'})
    for cur in curs:
        fits_filename = os.path.join(cur['path'], cur['filename'])
        logger.info(f'Removing older version of FITS file : {fits_filename}')
        if not dry_run:
            os.unlink(fits_filename)
            fits_db.delete_one({'_id':cur['_id']})



def reprocess_low_latency(start_date, end_date):
    daterange = pd.date_range(start_date, end_date)
    tpad=60
    for d in daterange:
        start_unix = d.timestamp()
        end_unix = start_unix +86400
        date_str=d.date().strftime("%Y-%m-%d")
        try:
            #remove_fits_except_science(start_unix-60, end_unix+60)
            fc.create_daily_low_latency_fits(date_str, tpad=60, version=3)

        except Exception as e:
            logger.error(str(e))

def process_science(bsd_id_start, bsd_id_end):
    for _id in range(bsd_id_start, bsd_id_end):
        bsd = bsd_db.find_one({
            '_id': _id
        })
        if not bsd:
            continue

        spid=bsd['SPID']
        if spid == 54143:
            continue
        
        uid = bsd.get('unique_id',None)
        if not uid:
            #don't recreate those aspect files
            continue
    
        pids = bsd['packet_ids']
        pids.sort()
        pkts = db.get_packets_by_ids(pids)
        logger.info(
            f'Creating fits file for bsd #{bsd["_id"]}, number of packets: {len(pids)}'
        )
        file_id = bsd['run_id']

        product = sdt.SPID_MAP[spid]
        logger.info(f'{spid},{product}')
        start_unix_time=bsd['start_unix_time']
        remove_science_fits(uid, start_unix_time)
        fc.create_fits_for_packets(file_id,
                                pkts,
                                spid,
                                product,
                                True,
                                base_path_name=PATH,
                                overwrite=True,
                                version=3,
                                remove_duplicates=True)

#reprocess_low_latency('2020-04-12T00:00:00', '2022-10-26T00:00:00')
#process_science(1, 14161)
