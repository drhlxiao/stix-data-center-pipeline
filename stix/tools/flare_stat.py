#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Flare detection algorithm
    Procedure:
    1) smoothing the lowest energy bin lightcurve using a butterworth low-pass filter
    2) searching for local maxima in the smoothed lightcurve
    This routine is called after a raw telemetry is being processed

"""

import os
import sys
from scipy import signal
import numpy as np
import math
import matplotlib
from matplotlib import pyplot as plt
from stix.core import stix_datatypes as sdt
from stix.core import mongo_db as db
from stix.spice import stix_datetime
from stix.core import stix_logger
#from stix.spice import solo
from stix.analysis import ql_analyzer as qla
logger = stix_logger.get_logger()
matplotlib.use('Agg')

mdb = db.MongoDB()

PEAK_MIN_NUM_POINTS = 7  #  peak duration must be greater than 28 seconds, used to determine energy range upper limit,

SPID = 54118

terminal = False


def info(msg):
    if terminal:
        print(msg)
        return
    logger.info(msg)



def update_flare_stat(_id): 
    print('updating flare ',_id)
    fdb=mdb.get_collection('flares')
    doc=fdb.find_one({'_id':_id})
    if not doc:
        print('can not find flare ', _id)
        return
    start_unix=doc['start_unix']
    span=doc['end_unix']-start_unix
    packets = mdb.get_quicklook_packets('lc',start_unix, span, 'header.unix_time')
    if not packets:
        info(f'No QL LC packets found for run {file_id}')
        return None
    data=qla.LightCurveAnalyzer.parse(packets)

    
    peak_min_width=15,
    peak_min_distance=150,
    rel_height=0.9,
    snapshot_path='.'
    unix_time = data['time']
    lightcurve = data['lcs'][0]

    med = np.median(lightcurve)
    prominence = 2 * np.sqrt(med)
    noise_rms = np.sqrt(med)
    baseline = med

    height = med + 3 * np.sqrt(med)
    stat = mdb.get_nearest_qllc_statistics(unix_time[0], max_limit=500)
    result = {}
    LC_statistics = []
    if not stat:
        print('can not find background for ',_id)
        return

    if stat['std'][0] < 2 * math.sqrt(
            stat['median'][0]):  #valid background
        height = stat['median'][0] + 2 * stat['std'][0]
        baseline = stat['median'][0]
        prominence = 2 * stat['std'][0]
        noise_rms = stat['std'][0]

    flare_stats = {'bkg_time_range':  (stat['start_utc'], stat['end_utc'])}

    upper_bin = 0
    for ilc in range(0, 5):  
        flare_stat = {}
        flare_stat['bkg_median'] = stat['median'][ilc]
        flare_stat['bkg_sigma'] = stat['std'][ilc]
        lc_cnts = data['lcs'][ilc]
        flare_stat['signal_median'] = int(np.median(lc_cnts))
        flare_stat['signal_max'] = int(np.max(lc_cnts))
        flare_stat['signal_min'] = int(np.min(lc_cnts))
        num_2sigma = int(
            np.sum(
                lc_cnts > stat['median'][ilc] + 2 * stat['std'][ilc]))
        flare_stat['num_points_above_2sigma'] = num_2sigma
        if num_2sigma > PEAK_MIN_NUM_POINTS:  #longer than half minutes
            upper_bin = ilc
        flare_stat['num_points_above_3sigma'] = int(
            np.sum(
                lc_cnts > stat['median'][ilc] + 3 * stat['std'][ilc]))
        flare_stats['lc' + str(ilc)] = flare_stat
    flare_stats['upper_ilc'] = upper_bin

    print(flare_stats)
        
    fdb.update_one({'_id': _id},{'$set':{'LC_statistics':flare_stats}})



if __name__ == '__main__':
    import sys
    terminal = True
    if len(sys.argv) < 2:
        print('flare_detection flare id')
    elif len(sys.argv) == 2:
        res = update_flare_stat(int(sys.argv[1]))
    else:
        for i in range(int(sys.argv[1]),
                int(sys.argv[2])):
            update_flare_stat(i)
