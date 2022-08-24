import os
import sys

import numpy as np
from stix.core import mongo_db
from stix.spice import time_utils as sdt
mdb= mongo_db.MongoDB()
db=mdb.get_db()
cal_db=db['calibration_runs']
flare_db=db['flares']

def is_sun_quiet(start_unix, duration, max_counts):
    fdocs=flare_db.find({'peak_unix_time':{'$gt': start_unix,'$lte':start_unix+duration}})
    return fdocs

for doc in cal_db.find():
    start_unix=doc['start_unix_time']
    duration=doc['duration']
    try:
        counts=np.array(doc['counts'][0])
    except Exception as e:
        continue
    cfl_counts=counts[8*12:8*12+8]
    min_cnts, max_cnts=np.min(cfl_counts), np.max(cfl_counts)


    if max_cnts>min_cnts+10*np.sqrt(max_cnts) and is_sun_quiet(start_unix, duration, 1e4):
        print('Occulted flares calibration run: ', 2*(max_cnts-min_cnts)/(min_cnts+max_cnts), doc['_id'])









