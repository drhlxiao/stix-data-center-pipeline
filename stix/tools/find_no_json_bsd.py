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
from stix.spice import datetime as sdt
from stix.fits import fits_creator
from stix.analysis import sci_packets_analyzer as spa 
import matplotlib.ticker as mticker
from pathlib import Path
import pymongo
mdb = db.MongoDB()
connect = pymongo.MongoClient()
db = connect["stix"]
col_bsd=db['bsd']
col_fits=db['fits']

#start_id=0
#end_id=100000
#bsd_query={'_id':{'$gte':start_bsd_id, 'lte':end_bsd_id}}
bsd_query={'name':'L1'}

for doc in col_bsd.find(bsd_query).sort('_id',-1):
    level1=doc.get('level1',None)
    to_do=False
    if level1 is None:
        to_do=True
    elif not Path(level1).is_file():
        to_do=True
    if to_do:
        print("creating json :", doc['_id'])
        spa.process_science_request(doc['_id'])
