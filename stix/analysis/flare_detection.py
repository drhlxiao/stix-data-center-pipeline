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
from stix.core import datatypes as sdt
from stix.core import mongo_db as db
from stix.spice import time_utils as st
from stix.core import logger
#from stix.spice import solo
from stix.analysis import ql_analyzer as qla
logger = logger.get_logger()
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


def merge_intervals(x):
    #sort the intervals by its first value
    x.sort(key=lambda x: x[0])
    m = []
    m.append(x[0])
    overlaps = lambda a, b: b[0] > a[0] and b[0] < a[1]
    for i in range(1, len(x)):
        pm = m.pop()
        if overlaps(pm, x[i]):
            m.append((pm[0], max(pm[1], x[i][1])))
        else:
            m.append(pm)
            m.append(x[i])
    return m


def smooth(y, N=15):

    y_padded = np.pad(y, (N // 2, N - 1 - N // 2), mode='edge')
    y_smooth = np.convolve(y_padded, np.ones((N, )) / N, mode='valid')
    return y_smooth


def get_lightcurve_data(file_id):
    packets = mdb.select_packets_by_run(file_id, SPID)
    if not packets:
        info(f'No QL LC packets found for run {file_id}')
        return None
    return qla.LightCurveAnalyzer.parse(packets)

def make_lightcurve_snapshot(data, docs, snapshot_path):
    '''
                '_id': first_id + i,
                'run_id': result['run_id'],
                'hidden': hidden,
                'peak_counts': result['peak_counts'][i],
                'peak_utc': result['peak_utc'][i],
                'peak_unix_time': result['peak_unix_time'][i],
    '''
    #print(flare_list)
    inserted_ids = docs['inserted_ids']
    #print(inserted_ids)

    for i, inserted_id in enumerate(inserted_ids):
        if inserted_id == None:
            continue
        _id = inserted_id
        #print('_id', inserted_id)
        min_height = docs['conditions']['min_height']

        start_unix = docs['start_unix'][i]
        end_unix = docs['end_unix'][i]
        flare_id = docs['flare_id'][i]
        margin = 300
        where = np.where((data['time'] > start_unix - margin)
                         & (data['time'] < end_unix + margin))
        unix_ts = data['time'][where]
        t_since_t0 = unix_ts - docs['peak_unix_time'][i]
        lc = data['lcs'][0][where]
        peak_counts = docs['peak_counts'][i]

        fig = plt.figure(figsize=(6, 2))
        plt.plot(t_since_t0, lc, label="4-10 keV LC")
        plt.plot(t_since_t0,
                 data['lc_smoothed'][where],
                 label='1-min moving mean')
        T0 = st.unix2utc(docs['peak_unix_time'][i])
        xmin = docs['start_unix'][i] - docs['peak_unix_time'][i]
        xmax = docs['end_unix'][i] - docs['peak_unix_time'][i]
        t70=[docs['PH70_unix'][i][0]-docs['peak_unix_time'][i],
                docs['PH70_unix'][i][1]-docs['peak_unix_time'][i]]

        plt.plot([0], [peak_counts], marker='+', color='cyan', markersize=15)

        ylow = peak_counts - 1.1 * docs['properties']['prominences'][i]

        plt.axvline(xmin, linestyle='dashed', color='C3')
        plt.axvline(xmax, linestyle='dashed', color='C3')

        plt.vlines(t70, ymin= 0.4 * peak_counts, ymax=0.9*peak_counts, linestyle='dashed', color='b')

        plt.text(xmin, 0.9 * peak_counts, 'Start', color='C3')
        plt.text(0.8 * xmax, 0.9 * peak_counts, 'End', color='C3')
        baseline = docs['baseline']
        plt.hlines(docs['properties']['width_heights'][i],
                   xmin=xmin,
                   xmax=xmax,
                   linewidth=2,
                   color='C2')

        plt.xlabel(f'Seconds since {T0}')
        plt.ylabel('Counts in 4 s')
        plt.title(f'Flare #{flare_id}')
        filename = os.path.join(snapshot_path,
                                f'flare_lc_{_id}_{flare_id}.png')
        fig.tight_layout()
        #print(filename)
        plt.savefig(filename, dpi=300)
        mdb.set_flare_lc_filename(_id, filename)
        plt.close()
        plt.clf()


def major_peaks(lefts, rights):
    #remove small peaks
    num = lefts.size
    major = [True] * num
    for i in range(num):
        a = (lefts[i], rights[i])
        for j in range(num):
            if i == j:
                continue
            b = (lefts[j], rights[j])

            if a[0] >= b[0] and a[1] <= b[1]:  #a in b
                major[i] = False
    return major


def find_flares_in_one_file(run_id,
           peak_min_width=15,
           peak_min_distance=150,
           rel_height=0.9,
           snapshot_path='.'):
    data = get_lightcurve_data(run_id)
    print(f'Deleting flares in file #{run_id}')
    mdb.delete_flares_of_file(run_id)
    if not data:
        info(f'No QL LC packets found for file {run_id}')
        return 0
    auxilary={'run_id':run_id}
    return find_flares_in_data(data, peak_min_width, peak_min_distance,rel_height, snapshot_path, auxilary)

def find_flares_in_data(data, 
           peak_min_width=15,
           peak_min_distance=150,
           rel_height=0.9,
           snapshot_path='.', auxilary=None):

    unix_time = data['time']
    lightcurve = data['lcs'][0]
    med = np.median(lightcurve)
    prominence = 2 * np.sqrt(med)
    noise_rms = np.sqrt(med)
    baseline = med

    height = med + 3 * np.sqrt(med)
    stat = mdb.get_nearest_lc_stats(unix_time[0], max_limit=500)
    if stat:
        #only use the lowest lightcurve
        if stat['std'][0] < 2 * math.sqrt(
                stat['median'][0]):  #valid background
            height = stat['median'][0] + 2 * stat['std'][0]
            baseline = stat['median'][0]
            prominence = 2 * stat['std'][0]
            noise_rms = stat['std'][0]

    lc_smoothed = smooth(lightcurve)
    result = {}
    xpeaks, properties = signal.find_peaks(
        lc_smoothed,
        height=height,
        prominence=prominence,
        width=peak_min_width,
        rel_height=rel_height,  #T90 calculation
        distance=peak_min_distance)
    if xpeaks.size == 0:
        info(f'No peaks found for file {run_id}')
        return 0
    info('Number of peaks:{}'.format(xpeaks.size))
    conditions = {
        'peak_min_width': peak_min_width,
        'peak_min_distance': peak_min_distance,
        'rel_height': rel_height,
        'prominence': prominence,
        'bkg_statistics': stat,
        'min_height': height,
    }
    doc = {}


    peak_values = properties['peak_heights']
    peak_unix_times = unix_time[xpeaks]
    peaks_utc = [st.unix2utc(x) for x in peak_unix_times]
    flare_ids = [
        st.unix2datetime(x).strftime("%y%m%d%H%M")
        for x in peak_unix_times
    ]

    majors = major_peaks(properties['left_ips'], properties['right_ips'])
    range_indexs = [(int(properties['left_ips'][i]),
                     int(properties['right_ips'][i]))
                    for i in range(xpeaks.size)]
    total_counts = [int(np.sum(lightcurve[r[0]:r[1]])) for r in range_indexs]
    LC_statistics = []
    if stat:#calculate statistics for all light curves, used for data requests
        for ipeak in range(xpeaks.size):
            flare_stats = {'bkg_time_range':  (stat['start_utc'], stat['end_utc'])}
            upper_bin = 0
            for ilc in range(0, 5):  
                flare_stat = {}
                flare_stat['bkg_median'] = stat['median'][ilc]
                flare_stat['bkg_sigma'] = stat['std'][ilc]
                r = range_indexs[ipeak]
                lc_cnts = data['lcs'][ilc][r[0]:r[1]]
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
            LC_statistics.append(flare_stats)

    seconds_per_bin = 4
    durations = np.array([(r[1] - r[0]) * seconds_per_bin
                          for r in range_indexs])
    cps = total_counts / durations

    total_signal_counts = []
    for i, (cnts, dur) in enumerate(zip(total_counts, durations)):
        sig_cnts = cnts - dur * baseline / seconds_per_bin
        if sig_cnts < 0:
            sig_cnts = cnts - dur * properties['width_heights'][
                i] / seconds_per_bin
        total_signal_counts.append(sig_cnts)

    flare_start_unix=unix_time[properties['left_ips'].astype(int)].tolist()
    flare_end_unix=unix_time[properties['right_ips'].astype(int)].tolist()

    #peak_width_T90=np.vstack((unix_time[properties['left_ips'].astype(int)], unix_time[properties['right_ips'].astype(int)])).T

    H70_res= signal.peak_widths(lc_smoothed, xpeaks, rel_height=0.7) #H70 calculation
    peak_width_H70= np.vstack((unix_time[H70_res[2].astype(int)], unix_time[H70_res[3].astype(int)])).T

    H80_res= signal.peak_widths(lc_smoothed, xpeaks, rel_height=0.7) #H70 calculation
    peak_width_H80= np.vstack((unix_time[H80_res[2].astype(int)], unix_time[H80_res[3].astype(int)])).T

    H50_res= signal.peak_widths(lc_smoothed, xpeaks, rel_height=0.5) #H50 calculation
    peak_width_H50= np.vstack((unix_time[H50_res[2].astype(int)], unix_time[H50_res[3].astype(int)])).T
    #calculate peak width at 50% of the maximum count

    doc = {
        'num_peaks': xpeaks.size,
        'peak_unix_time': peak_unix_times.tolist(),
        'peak_counts': peak_values.tolist(),
        'peak_utc': peaks_utc,
        'flare_id': flare_ids,
        'baseline': baseline,
        'noise_rms': noise_rms,
        'total_counts': total_counts,
        'total_signal_counts': total_signal_counts,
        'duration': durations.tolist(),
        'mean_cps': cps.tolist(),
        'conditions': conditions,
        'peak_width_bins': properties['widths'].tolist(),  #width of height
        'width_height': properties['width_heights'].tolist(
        ),  # height of the width, background level
        'peak_prominence': properties['prominences'].tolist(),
        #'peak_T90_widths': peak_width_T90.tolist(),
        'PH70_unix': peak_width_H70.tolist(),
        'PH80_unix': peak_width_H80.tolist(),
        'PH50_unix': peak_width_H50.tolist(),
        'start_unix':flare_start_unix,
        'end_unix': flare_end_unix,
        'is_major': majors,
        'LC_statistics': LC_statistics,
    }
    if isinstance(auxilary, dict):
        doc.update(auxilary)

    mdb.save_flare_info(doc)
    doc['properties'] = properties
    data['lc_smoothed'] = lc_smoothed
    make_lightcurve_snapshot(data, doc, snapshot_path)
    return xpeaks.size

def find_flares_in_files(fid_start, fid_end, img_path='/data/flare_lc'):
    for i in range(fid_start, fid_end + 1):
        print(f'deleting flares of Files {i}')
        mdb.delete_flares_of_file(i)
    for i in range(fid_start, fid_end + 1):
        search_flares(i, snapshot_path=img_path)


if __name__ == '__main__':
    import sys
    terminal = True
    if len(sys.argv) < 2:
        print('flare_detection file_number')
    elif len(sys.argv) == 2:
        res = search_flares(int(sys.argv[1]), snapshot_path='/data/flare_lc')
        print('Number of peaks:', res)
    else:
        find_flares_in_files(int(sys.argv[1]),
                       int(sys.argv[2]),
                       img_path='/data/flare_lc')
