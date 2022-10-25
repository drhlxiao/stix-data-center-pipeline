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
dry_run=True
PATH='/data/fits'
fits_db= db.collection_fits
def remove_fits_except_science(start_unix, end_unix):
    curs=fits_db.find({'data_start_unix':{'$lte':end_unix}, 'data_end_unix':{'$gte' :start_unix}, 'product_group':{'$ne': 'science'}  })
    for cur in curs:
        fits_filename = os.path.join(cur['path'], cur['filename'])
        logger.info(f'Removing older version of FITS file : {fits_filename}')
        if not dry_run:
            os.unlink(fits_filename)
            fits_db.delete_one({'_id':cur['_id']})

def remove_science_fits(request_id, start_unix, tpad = 300 ):
    curs=fits_db.find({'data_start_unix':{'$lte':start_unix-tpad}, 'data_start_unix':{'$gte' :start_unix + tpad}, 'product_group': 'science', 'request_id':request_id })
    for cur in curs:
        fits_filename = os.path.join(cur['path'], cur['filename'])
        logger.info(f'Removing older version of FITS file : {fits_filename}')
        if not dry_run:
            os.unlink(fits_filename)
            fits_db.delete_one({'_id':cur['_id']})



def reprocess_low_latency(start_date, end_date, remove=False):
    daterange = pd.date_range(start_date, end_date)
    for d in daterange:
        start_unix = d.timestamp()
        end_unix = 86400 + start_unix
        logger.info(f'Creating daily fits file for {d}')
        try:
            remove_fits(start_unix,end_unix)
            fc.create_continous_low_latency_fits(start_unix,
                                      end_unix,
                                      output_path=PATH, 
                                      overwrite=True,
                                      version=1,
                                      run_type='daily')

        except Exception as e:
            logger.error(str(e))

def preprocess_science(bsd_id_start, bsd_id_end):
    for _id in range(bsd_id_start, bsd_id_end):
        bsd = bsd_db.find({
            '_id': _id
        }).max_time_ms(300 * 1000)
        pids = bsd['packet_ids']
        pids.sort()
        pkts = db.get_packets_by_ids(pids)
        logger.info(
            f'Creating fits file for bsd #{bsd["_id"]}, number of packets: {len(pids)}'
        )
        file_id = bsd['run_id']
        uid = bsd['unique_id']
        spid = bsd['SPID']
        product = sdt.SPID_MAP[spid]
        logger.info(f'{spid},{product}')
        start_unix_time=bsd['start_unix_time']
        remove_science_fits(uid, start_unix_time)
        create_fits_for_packets(file_id,
                                pkts,
                                spid,
                                product,
                                True,
                                output_path=PATH,
                                overwrite=True,
                                version=3,
                                remove_duplicates=True)

reprocess_low_latency('2020-04-14T00:00:00', '2020-10-25:00:00')
