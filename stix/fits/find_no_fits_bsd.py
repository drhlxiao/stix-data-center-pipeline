#find bulk science data that don't have fits file created    
# hualin.xiao created on Feb. 11, 2022
import sys
import os
import json
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime
from matplotlib import pyplot as plt
from stix.core import mongo_db as db
from stix.spice import time_utils as sdt
from stix.fits import fits_creator
import matplotlib.ticker as mticker
from pathlib import Path
import pymongo
mdb = db.MongoDB()
connect = pymongo.MongoClient()
db = connect["stix"]
col_bsd=db['bsd']
col_fits=db['fits']


def create_fits(bsd_id_start):
    bsd_id_end=bsd_id_start
    print('Creating fits for :',bsd_id_start)
    try:
        fits_creator.create_fits_for_bulk_science(bsd_id_start, bsd_id_end,output_path='/data/fits')
    except Exception as e:
        raise

def main(num_docs=1000):
    cursor = col_bsd.find().sort('_id',-1).limit(num_docs) if not None else col_bsd.find().sort('_id',-1)
    for doc in  cursor:
        request_id=doc.get('unique_id',None)
        if request_id is None:
            create_fits(doc['_id'])
            continue
        fdoc=col_fits.find_one({'request_id':request_id})
        if not fdoc:
            create_fits(doc['_id'])
            continue

        if not Path(fdoc['path'],fdoc['filename']).is_file():
            create_fits(doc['_id'])
if __name__=='__main__':
    print('A script to create FITS files for bulk science data that do not have  FITS files created') 
    main()
