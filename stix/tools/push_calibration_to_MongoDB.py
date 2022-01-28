import os
import sys

import numpy as np
from stix.core import mongo_db
from stix.spice import datetime as sdt
mdb= mongo_db.MongoDB()
db=mdb.get_db()
cal_db=db['calibration_runs']
flare_db=db['flares']

def is_sun_quiet(start_unix, duration, max_counts):
    fdocs=flare_db.find({'peak_unix_time':{'$gt': start_unix,'$lte':start_unix+duration},'total_counts':{'$lt':max_counts}})
    return fdoc

for doc in cal_db.find():
    start_unix=doc['start_unix_time']
    duration=doc['duration']
    counts=np.array(doc['counts'][1])
    cfl_counts=counts[8*12:8*12+8]
    min_cnts, max_cnts=np.min(cfl_counts), np.max(cfl_counts)


    if max_cnts>2*min_cnts and is_sun_quiet(start_unix, duration, 2e4):
        print('Occulted flares calibration run', doc['_id'])









