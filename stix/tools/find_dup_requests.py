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
#bsd_query={'name':'xray-cpd'}
bsd_query={}
total=0
last=0
last_id=0
last_type=''
fnames=[]
for doc in col_bsd.find(bsd_query).sort('start_unix_time',1):
    start=doc.get('start_unix_time',None)
    header_time =doc.get('header_unix_time',None)
    start_utc=sdt.unix2utc(start)
    SPID=doc.get('name',None)
    start=int(start)
    request_id=doc.get('unique_id',None)
    if request_id is None:
        continue
    _id=doc['_id']

    fname=f'{start_utc}_{request_id}_{SPID}'
    fnames.append(fname)


fnames.sort()
last_name=''

for fname in fnames:
    if last_name==fname:
        total+=1
        print("Duplicates:", fname)
    last_name=fname
print("Total duplicates:", total)


    


