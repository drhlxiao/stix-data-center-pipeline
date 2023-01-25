# move fits files to different folders instead of a single file
# hualin.xiao created on Feb. 11, 2022
import sys
import os
import json
import numpy as np
import shutil

from datetime import datetime
from stix.core import mongo_db as db
from stix.spice import time_utils as sdt
from pathlib import Path
import pymongo
mdb = db.MongoDB()
connect = pymongo.MongoClient()
db = connect["stix"]
db_fits=db['fits']
for doc in db_fits.find():
    start=doc.get('data_start_unix')
    dt= sdt.unix2datetime(start)
    sub_folder=dt.strftime("%Y%m")
    old_path=doc['path']
    if sub_folder in old_path:
        print("sub folder exists:", old_path)
        continue
    new_path =os.path.join(doc['path'], sub_folder)
    os.makedirs(new_path, exist_ok=True)
    old_filename= os.path.join(doc['path'], doc['filename']) 
    new_filename=os.path.join(new_path, doc['filename']) 
    try:
        shutil.move( old_filename, new_filename)        
        print('mv ', doc['_id'],  old_filename, new_filename)
        db_fits.update_one({'_id':doc['_id']}, {'$set':{
            'path':new_path
            }}, upsert=False)
    except Exception as e:
        print(e)


