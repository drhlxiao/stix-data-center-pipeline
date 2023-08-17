import sys
import numpy as np
from stix.core import datatypes as sdt
from stix.spice import time_utils
from stix.core import mongo_db as mdb
import pymongo
from stix.core import metadata as m
from datetime import datetime

m_db=mdb.MongoDB()
lcdb=m_db.get_collection('ql_lightcurves')

sp=m.StixQuickLookPacketAnalyzer(m_db)

def process_file(fid):
    docs=lcdb.find({}).sort('_id',-1)
    for lc in docs:
        time=lc['time']
        rcr=lc['rcr']
        print(time_utils.unix2utc(time[0]))
        sp.record_rcr_status_change_events(time,rcr)


for i in range(3305):
    print(i)
    process_file(i)
    

