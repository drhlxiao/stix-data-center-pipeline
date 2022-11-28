#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : science data pipeline
# @description  : create fits file, create science data snapshot, and submit imaging task
# @author       : Hualin Xiao
# @date         : Nov. 28, 2022
#

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
from stix.flare_pipeline import imaging_task_manager as itm
from stix.pipeline import task_manager 
from stix.analysis import sci_packets_analyzer as scia
from stix.core import logger
from pathlib import Path
mdb = db.MongoDB()

logger = logger.get_logger()

def create_fits(bsd_id_start):
    bsd_id_end=bsd_id_start
    print('Creating fits for :',bsd_id_start)
    try:
        fits_creator.create_fits_for_bulk_science(bsd_id_start, bsd_id_end,output_path='/data/fits')
    except Exception as e:
        raise

def create_fits_recent(last_num=1000, reproduce_all=False):
    bsd_query={}
    col_bsd=mdb.get_collection('bsd')
    col_fits=mdb.get_collection('fits')
    for doc in col_bsd.find(bsd_query).sort('_id',-1).limit(last_num):
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
            logger.info("Generation cache for science data ")
            scia.process_science_request_doc(doc)
            itm.register_imaging_task_for_science_data([doc['_id']])
