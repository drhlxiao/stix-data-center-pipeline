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
from stix.analysis import sci_packets_analyzer as scip
import matplotlib.ticker as mticker
from pathlib import Path
import pymongo
mdb = db.MongoDB()
connect = pymongo.MongoClient()
db = connect["stix"]
col_bsd=db['bsd']

#start_id=0
#end_id=100000
#bsd_query={'_id':{'$gte':start_bsd_id, 'lte':end_bsd_id}}
bsd_query={"SPID":{'$exists':True}, 'unique_id':{'$exists':True}}

def create_snapshot(doc):
    try:
        print("processing:", doc['_id'])
        scip.process_science_request_doc(doc)
    except Exception as e:
        raise

for doc in col_bsd.find(bsd_query).sort('_id',-1):
    request_id=doc.get('unique_id',None)
    spid=doc.get('SPID',None)
    if request_id is None or spid is None:
        continue
    snapshot=doc.get('level1', None)
    if snapshot is None:
        create_snapshot(doc)
    elif  not Path(snapshot).is_file():
        create_snapshot(doc)
