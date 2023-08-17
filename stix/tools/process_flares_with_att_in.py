""""
Here is how the flares with ATT in is processed:
    1) ql_analyzer, goes through all QL data, and find periods when ATT was in
    2) Period info is recorded in the event database
    3) this script does count correction for the periods with ATT inserted
    4) corrected counts are recorded in another database "qlc_att
    5) this scripts delete flares closed to the periods with ATT inserted
    6) load longer LCs and reconstruct flares
    7) estimate GOES class for the new inserted flares
    8) Web page updated
    9) You need to flag "data_requested" for flares in the database, otherwise,new data requests will be created
    
"""

import os
import sys
import math
import matplotlib
import numpy as np
from scipy import signal
from pprint import pprint
from datetime import datetime, timedelta

from stix.spice import time_utils as st
from stix.utils import bson as bs
from stix.core import mongo_db as db

from stix.analysis import correct_qlc_with_att as cqlc
from stix.analysis import flare_detection as fld

mdb = db.MongoDB()

start_utc, end_utc='2020-07-15T00:00:00','2023-07-17T00:00:00'
start_unix, end_unix = st.utc2unix(start_utc), st.utc2unix(end_utc)

#print processing new flares
print("correcting counts for ATT...")
cqlc.process_new(start_utc,end_utc)

print(f"Find ATT in time range in:  {start_utc} - {end_utc}")
time_ranges = mdb.find_att_in_time_ranges(start_unix, end_unix)
for tr in time_ranges:
    tgap=80
    start_utc = st.unix2utc(tr[0] - tgap)
    end_utc = st.unix2utc(tr[1] + tgap)
    print("reprocess flares:",start_utc," - ", end_utc)
    cqlc.remove_flares_with_attenuator_in(tr[0], tr[1])
    #remove reported flares -80 sec and +80 sec ATT inserted
    lc_margin = 4*3600
    fld.find_flares_in_time_range(tr[0]- lc_margin ,tr[1]+lc_margin)
    # margin is required for flare identification 




