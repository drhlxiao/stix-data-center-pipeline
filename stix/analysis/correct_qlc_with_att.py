"""
a script to correct light curves for ATT inserted
@Author: hualin.xiao@fhnw.ch
@Date: 2023-08-10

"""

import os
import sys
sys.path.append('.')
from scipy import signal
from pprint import pprint
import numpy as np
import math
import matplotlib
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from stix.core import datatypes as sdt
from stix.core import mongo_db as db
from stix.utils import bson as bs
from stix.spice import time_utils as st

from stix.analysis import ql_analyzer as qla
from stix.core import logger
logger = logger.get_logger()

mdb = db.MongoDB()
#qla.set_db(mdb)
qlc_att_db=mdb.get_collection('qlc_att_in')

def correct_ql_counts_in_time_range(start_unix, end_unix, tgap=20):
    """
    correct counts for periods for ATT is in,
    we need to specify for a gap to get the max counts before ATT inserted for correction
    """
    start_utc, end_utc=st.unix2utc(start_unix), st.unix2utc(end_unix)
    logger.info(f"Find ATT in time range in:  {start_utc} - {end_utc}")

    time_ranges=mdb.find_att_in_time_ranges(start_unix, end_unix)
    for tr in time_ranges:
        start_utc=st.unix2utc(tr[0]-tgap) 
        end_utc=st.unix2utc(tr[1]+tgap) 
        logger.info(f"Correcting LCs for {start_utc} - {end_utc}")
        try:
            res=qla.LightCurveMerger.from_database(start_utc, end_utc)
            doc=res.correct_att_in_counts(tgap=tgap)
        except Exception as e:
            print(e)
            continue
         
        doc.update({
                'start_unix':tr[0],
                'end_unix':tr[1],
                'start_utc':start_utc,
                'end_utc':end_utc,
                })
        
        doc=bs.dict_to_json(doc)
        qlc_att_db.update_one({'start_unix':tr[0],'end_unix':tr[1]},{'$set':doc},upsert=True)
def process_new(start_unix=None, end_unix=None):
    """
    make sure that data already processed is not processed again
    """
    if start_unix is None or end_unix is None:
        docs=qlc_att_db.find().sort('end_unix',-1).limit(1)
        docs=list(docs)
        end_unix =st.get_now('unix')
        try:
            start_unix = docs[0]['end_unix']
        except Exception as e:
            logger.error(str(e))
            start_unix=None
            end_unix=None

    if start_unix is not None and end_unix is not None:
        correct_ql_counts_in_time_range(start_unix, end_unix)
        return


    logger.warn("Failed to determine start time or end time. You did not specified or could not find the information in the database")


    


if __name__=='__main__':
    start, end = None, None
    if len(sys.argv)==3:
        start=st.utc2unix(sys.argv[1])
        end=st.utc2unix(sys.argv[2])
    process_new(start, end)





