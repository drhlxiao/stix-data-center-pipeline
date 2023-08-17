#find bulk science data that don't have fits file created    
# hualin.xiao created on Feb. 11, 2023
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
connect = pymongo.MongoClient(port=9123)
db = connect["stix"]
col_bsd=db['bsd']
db_req=db['data_requests']


#start_id=0
#end_id=100000
bsd_query={'_id':{'$gte':23662, '$lte':23699}}
#bsd_query={'name':'xray-cpd'}
total=0
last=0
last_id=0
last_type=''
fnames=[]
for doc in col_bsd.find(bsd_query).sort('_id',1):
    start=doc.get('start_unix',None)
    end=doc.get('end_unix',None)
    #print(doc.get('unique_id',None))
    doc_req = list(db_req.find({'start_unix': {'$gte':start-60, '$lte':start+60},
        'end_unix': {'$gte':end-60, '$lte':end+60}, 
        }))
    if len(doc_req)==1:
        db_req.update_one({'_id': doc_req[0]['_id']},{'$push':{'unique_ids':doc['unique_id']}})
        #print({'_id': doc_req[0]['_id']},{'$push':{'unique_ids':doc['unique_id']}})
    elif doc_req:
        print('multiple found: ')
        print(doc_req)
        print('end')
    else:
        print(f'Not found for {doc["_id"]}')




