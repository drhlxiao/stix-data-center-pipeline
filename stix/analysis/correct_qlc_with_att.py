"""
a script to correct light curves for ATT inserted
@Author: hualin.xiao@fhnw.ch
@Date: 2023-08-10

"""
import os
import sys
import math
import numpy as np
import matplotlib
from scipy import signal
from pprint import pprint
from datetime import datetime, timedelta
from matplotlib import pyplot as plt

from stix.core import logger
from stix.analysis import ql_analyzer as qla
from scipy import interpolate
from stix.spice import time_utils as st
from stix.utils import bson as bs
from stix.core import mongo_db as db
from stix.core import datatypes as sdt

sys.path.append('.')


logger = logger.get_logger()

mdb = db.MongoDB()
qlc_att_db = mdb.get_collection('qlc_att_in')
flare_db = mdb.get_collection('flares')


def remove_flares_with_attenuator_in(start_unix, end_unix, tgap=80):
    """
    delete solar flares reported when att was in
    """
    # light curves are smoothed
    logger.info(
        f'changing status of flare detected from {st.unix2utc(start_unix)} to {st.unix2utc(end_unix)}'
    )


    flare_db.update_many(
        {
            'peak_unix_time': {
                '$gte': start_unix - tgap,
                '$lte': end_unix + tgap
            }
        }, {
            '$set': {
                'hidden': True,
                'comment': 'removed because ATT was inserted'
            }
        })
    # we don't really delete them





def correct_ql_counts_in_time_range(start_unix, end_unix, tgap=20):
    """
    correct counts for periods for ATT is in,
    we need to specify for a gap to get the max counts before ATT inserted for correction

    this must be executed after flare detection in the pipeline

    """
    start_utc, end_utc = st.unix2utc(start_unix), st.unix2utc(end_unix)
    logger.info(f"Find ATT in time range in:  {start_utc} - {end_utc}")
    time_ranges = mdb.find_att_in_time_ranges(start_unix, end_unix)
    for tr in time_ranges:
        start_utc = st.unix2utc(tr[0] - tgap)
        end_utc = st.unix2utc(tr[1] + tgap)
        remove_flares_with_attenuator_in(tr[0], tr[1])

        logger.info(f"Correcting LCs for {start_utc} - {end_utc}")
        try:
            res = qla.LightCurveMerger.from_database(start_utc, end_utc)
            doc = res.correct_att_in_counts(tgap=tgap)
        except Exception as e:
            print(e)
            raise
            continue
        doc.update({
            'start_unix': tr[0],
            'end_unix': tr[1],
            'start_utc': start_utc,
            'end_utc': end_utc,
        })

        doc = bs.dict_to_json(doc)
        qlc_att_db.update_one({
            'start_unix': tr[0],
            'end_unix': tr[1]
        }, {'$set': doc},
            upsert=True)


def process_new(start=None, end=None):
    """
    make sure that data already processed is not processed again
    """
    

    if start is None or end is None:
        docs = qlc_att_db.find().sort('end_unix', -1).limit(1)
        docs = list(docs)
        end_unix = st.get_now('unix')
        try:
            start_unix = docs[0]['end_unix']
        except Exception as e:
            logger.error(str(e))
            start_unix = None
            end_unix = None
    else:
        try:
            start_unix, end_unix=float(start), float(end)
        except ValueError:
            start_unix = st.utc2unix(start)
            end_unix = st.utc2unix(end)

    if start_unix is not None and end_unix is not None:
        correct_ql_counts_in_time_range(start_unix, end_unix)
        return

    logger.warn(
        "Failed to determine start time or end time. You did not specified or could not find the information in the database"
    )


if __name__ == '__main__':
    start, end = None, None
    if len(sys.argv) == 3:
        start = st.utc2unix(sys.argv[1])
        end = st.utc2unix(sys.argv[2])
    process_new(start, end)
