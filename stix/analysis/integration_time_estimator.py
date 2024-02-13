#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Compute integration times
    Steps:
    1) Find the rotation buffer integration times
    2) estimate integration time  

"""

import os
import sys
from scipy import signal
import numpy as np
import math
from matplotlib import pyplot as plt
from stix.core import mongo_db as db
from stix.spice import time_utils
from stix.core import logger
from stix.spice import solo
from stix.analysis import ql_analyzer as qla

logger = logger.get_logger()

mdb = db.MongoDB()

default_rot_tmin=0.5
default_rot_tmax=20
default_rot_threshold=1800
#rotation buffer configuration


def get_rotating_buffer_config(max_time, min_time=None):
    rot_config={312: [], 313: [], 314:[]}
    config_db= mdb.get_collection('config')

    if not min_time:
        min_time=max_time



    for _id in rot_config.keys():
        # First Query
        query1 = {
            'type': 'asw',
            'parameter': _id,
            'execution_unix': { '$lte': min_time }
        }

        result1 = config_db.find(query1).sort('execution_unix', -1).limit(1)

        # Second Query
        query2 = {
            'type': 'asw',
            'parameter': _id,
            'execution_unix': { '$gte': min_time, '$lte': max_time }
        }

        result2 = config_db.find(query2)

        # Combine Results
        rows= list(result1) + list(result2)
        for row in rows:
            rot_config[_id].append( {'time': row['execution_unix'], 'value': row['value']})
    return rot_config

def get_rot_config(configs, unix_time):
    tmin=default_rot_tmin
    tmax=default_rot_tmax
    thr=default_rot_threshold
    

    for row in configs[312]:
        if unix_time > row['time']:
            tmin=0.1*(int(row['value'])+1)
    for row in configs[313]:
        if unix_time > row['time']:
            tmax=0.1*(int(row['value'])+1)
    for row in configs[314]:
        if unix_time > row['time']:
            thr=row['value']

    return tmin, tmax, thr


def process_file(file_id):
    print("Processing file:", file_id)
    packets = mdb.select_packets_by_run(file_id, SPIDs=54118)
    if not packets:
        return
    data = qla.LightCurveAnalyzer.parse(packets)
    if data is None:
        return
    unix_time = data['time']  # set to the center of a bin
    configs = get_rotating_buffer_config(unix_time[-1], unix_time[0])

    lightcurve = data['lcs'][0] + data['lcs'][1]  # 4-15 keV
    res = {'t': [], 'tbin': [], 'counts': [],  'num_tbins': []}
    s = 0
    last_time = 0

    

    for t, c in zip(unix_time, lightcurve):
        t = int(t)
        c = int(c)
        if last_time == 0:
            last_time = t
            continue

        min_tbin, max_tbin, thr=get_rot_config(configs, t)
        print("Final settings:", min_tbin,max_tbin,thr)


        if s > 0 and (s >= thr or t - last_time >= max_tbin
                or c >= thr):  # close the opening time bin anyway
            # close last time bin
            tbin = round(min([max_tbin, t - last_time]), 1)
            res['tbin'].append(tbin)
            res['num_tbins'].append(1)
            res['t'].append([last_time, t])
            res['counts'].append(s)
            # print(last_time, t,t-last_time,  c, tbin, s)
            s = 0
            last_time = t

        if c < thr:
            s += c
        else:
            num_tbins = math.ceil(c / thr)
            max_bins = math.ceil(4. / min_tbin)
            num_tbins = 1 if num_tbins < 1 else num_tbins
            num_tbins = max_bins if num_tbins > max_bins else num_tbins
            tbin = round(4. / num_tbins, 1)  # time step 0.1s
            res['tbin'].append(tbin)
            res['num_tbins'].append(num_tbins)
            res['t'].append([t - 4, t])
            res['counts'].append(c / num_tbins)
            # print(last_time, t, t-last_time, c, tbin, c/num_tbins)
            s = 0
            last_time = t

    res['start_unix'] = res['t'][0][0]
    res['end_unix'] = res['t'][-1][1]
    res['start_utc'] = time_utils.unix2utc(res['t'][0][0])
    res['end_utc'] = time_utils.unix2utc(res['t'][-1][1])
    res['file_id'] = file_id
    mdb.insert_time_bins(res, delete_existing=True)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('estimate_integration_times file_id_start, file_id_end')
    elif len(sys.argv) >= 2:
        start_id = int(sys.argv[1])
        end_id = start_id
        if len(sys.argv) == 3:
            end_id = int(sys.argv[2])
        for i in range(start_id, end_id + 1):
            process_file(i)
