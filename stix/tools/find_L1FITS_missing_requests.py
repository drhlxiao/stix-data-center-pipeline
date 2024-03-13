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

# Execute the aggregation pipeline
docs = col_bsd.find({})
f=open('L1_missing.csv','w')
f.write("Unique_ID,BSD_DB_Entry_ID, Data Start UTC\n")
num=0
for doc in docs:
    if 'unique_id' not in doc:
        continue
    fdoc=col_fits.find_one({'request_id':doc['unique_id'],'level':'L1'},{'_id':1})
    if fdoc:
        continue
    line=f'{doc["unique_id"]},{doc["_id"]},{sdt.unix2utc(doc["start_unix_time"])}\n'
    print(line)
    num+=1
    f.write(line)

print("Total files:",num)

