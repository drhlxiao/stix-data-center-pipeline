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
sys.path.append('.')
from scipy import signal
import numpy as np
import math
import matplotlib
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from stix.core import datatypes as sdt
from stix.core import mongo_db as db
from stix.spice import time_utils as st

from stix.analysis import ql_analyzer as qla
from stix.core import logger
logger = logger.get_logger()

mdb = db.MongoDB()
debug = False

PEAK_MIN_NUM_POINTS = 7  #  peak duration must be greater than 28 seconds, used to determine energy range upper limit,
niter=900
#matplotlib.use('qtagg' if debug else 'agg')
try:
    import ROOT
    ROOT_EXISTS=True
except ImportError:
    ROOT_EXISTS=False


def get_lightcurve_baseline(counts, niter):
    """
    background estimation based on the "Sensitive Nonlinear Iterative Peak (SNIP) clipping algorithm"
    """
    source = np.copy(counts)
    s= ROOT.TSpectrum()
    nbins=counts.size
    niter = min(int((nbins-1)/2), niter) 
    #niter must be > nbins-1/2 according to the manual
    if niter < 150:
        #this algorithm dosen't work 
        return None

    res=s.Background(source,nbins, niter, #number of iterations
            1, #direction, decreasing 1, increasing 0
            2, #filterOrder
            1, # boolean,  smoothing
            3, #smooth window
            0) #compton
    return source




def find_flare_time_ranges(lc_times, lc_counts, peaks, props, lc_baseline,  noise_rms):
    """
     calculate peak width at height 
     arguments
     lc_times: list
        light curve time series
    lc_count: list
        light curve counts
    """
    #print("Threshold:", threshold)

    num_peaks = peaks.size
    imax = len(lc_counts)
    boundaries = []
    left_ips, right_ips, start_times, end_times = [], [], [], []

    for k in range(num_peaks):
        peak = peaks[k]
        i = peak
        while i >= 0 and lc_baseline[i] +noise_rms <= lc_counts[i]:
            i -= 1
        left_ip = i
        i = peak
        while i < imax and lc_baseline[i] + noise_rms <= lc_counts[i]:
            i += 1
        right_ip = i

        if left_ip==right_ip:
            left_ip-=7
            right_ip+=8
            #at least last one minutes

        boundaries.append(True if left_ip <= 0 or right_ip >= imax else False)

        right_ip = min(right_ip, imax - 1)
        left_ip = max(0, left_ip)  #falls at the edge






        #print(i_min, i_max, left_ip, right_ip)
        left_ips.append(left_ip)
        right_ips.append(right_ip)
        start_times.append(lc_times[left_ip])
        end_times.append(lc_times[right_ip])
    return start_times, end_times, left_ips, right_ips, boundaries


def smooth(y, win=15):
    """
    moving averaging algorithm
    Parameters
    ----
    y: np.array
        data points
    win: int
        window

    """
    y_padded = np.pad(y, (win // 2, win - 1 - win // 2), mode='edge')
    y_smoothed = np.convolve(y_padded, np.ones((win, )) / win, mode='valid')
    return y_smoothed


def create_quicklook_plot(data, docs, lc_output_dir,lc_baseline=None, same_plot=False):
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
    tstart=data['time'][0]
    if same_plot:
        #fig,(ax, ax2)= plt.subplots(2,1)
        fig,ax= plt.subplots(1,1)
        lc = data['lcs'][0]
        dt=data['time']-tstart
        ax.plot(dt, lc, label="STIX 4-10 keV LC")
        ax.plot(dt, data['lc_smoothed'], label="Smoothed LC")
        if lc_baseline is not None:
            ax.plot(dt, lc_baseline, label="Baseline")
            #ax2.plot(dt, data['lc_smoothed']-lc_baseline, label="Baseline")
            #ax2.set_xlabel(f'Start at {st.unix2utc(tstart)}')
            #ax2.set_ylabel('Counts')
            #ax2.set_title('Baseline subtracted counts')
            #ax2.set_yscale('log')
            
    
    num_ids=len(inserted_ids)

    for i, inserted_id in enumerate(inserted_ids):
        if inserted_id == None:
            continue
        _id = inserted_id
        #print('_id', inserted_id)
        min_height = docs['conditions']['min_height']
        is_major = docs['is_major']

        start_unix = docs['start_unix'][i]
        end_unix = docs['end_unix'][i]
        flare_id = docs['flare_id'][i]
        margin = 300
        where = np.where((data['time'] > start_unix - margin)
                & (data['time'] < end_unix + margin))
        unix_ts = data['time'][where]
        T0 = docs['peak_unix_time'][i] if not same_plot else tstart
        t_since_t0 = unix_ts - T0
        lc = data['lcs'][0][where]
        peak_counts = docs['peak_counts'][i]

        if not same_plot:
            fig,ax= plt.subplots(1,1)
            ax.plot(t_since_t0, lc, label="4-10 keV LC")
            ax.plot(t_since_t0,
                data['lc_smoothed'][where],
                label='1-min moving mean')
        xmin = docs['start_unix'][i] - T0
        xmax = docs['end_unix'][i] - T0

        ax.plot([docs['peak_unix_time'][i] -T0], [peak_counts], marker='+', color='red')
        
        ax.axvspan(xmin, xmax, color='cyan', ec='b', alpha=0.5)

        filename = os.path.join(lc_output_dir,
                f'flare_lc_{_id}_{flare_id}.png')
        if not same_plot or i == num_ids-1:
            ax.set_xlabel(f'T [s] - Start at {st.unix2utc(T0)}')
            ax.set_ylabel('Counts')
            filename = os.path.join(lc_output_dir,
                f'flare_lc_{_id}.png')
            ax.set_title(f'Flare #{flare_id} (major: {is_major[i]})')
            ax.set_yscale('log')
            ax.legend()
            plt.tight_layout()
            plt.savefig(filename, dpi=300)
            logger.info(f'Creating file:{filename}')
            if debug:
                plt.show()
            plt.close()
            plt.clf()
        mdb.set_flare_lc_filename(_id, filename)


def find_major_peaks(lefts, rights, peak_values):
    """
        merge overlapped time ranges
    """
    num = len(lefts)
    #print("NUM:", num)
    major = [True] * num
    for i in range(num):
        a = (lefts[i], rights[i])
        for j in range(num):
            if i == j:
                continue
            b = (lefts[j], rights[j])
            if a[0] < b[1] and a[1] > b[0] and peak_values[i] < peak_values[j]:
                major[i] = False
    return major


def find_flares_in_one_file(run_id,
        peak_min_width=15,
        peak_min_distance=75,
        rel_height=0.9,
        lc_output_dir='.'):
    packets = mdb.select_packets_by_run(run_id, 54118)
    if not packets:
        logger.info(f'No QL LC packets found for Run {run_id}')
        return 0
    data = qla.LightCurveAnalyzer.parse(packets)
    if not data:
        logger.info(f'No QL LC packets found for Run {run_id}')
        return 0
    auxilary = {'run_id': run_id}
    return find_flares_in_data(data, peak_min_width, peak_min_distance,
            rel_height, lc_output_dir, auxilary)


def find_flares_in_recent_LC(start_date_off=-3, end_date_off=0, lc_path='.'):
    """
       it is possible that flares are not detected if a flare 
    """
    now_unix = datetime.now().timestamp()
    start_unix = now_unix + start_date_off * 86400
    end_unix = now_unix + end_date_off * 86400
    return find_flares_in_time_range(start_unix,
            end_unix,
            lc_output_dir=lc_path)


def find_flares_in_time_range(start_unix,
        end_unix,
        peak_min_width=15,
        peak_min_distance=75,
        rel_height=0.9,
        lc_output_dir='.'):
    packets = mdb.get_LC_pkt_by_tw(start_unix, end_unix - start_unix)
    if not packets:
        logger.info(f'No QL LC packets found ')
        return 0
    data = qla.LightCurveAnalyzer.parse(packets)
    auxilary = {'run_id': -1}
    if not data:
        logger.info(f'No QL LC packets found ')
        return 0
    return find_flares_in_data(data, peak_min_width, peak_min_distance,
            rel_height, lc_output_dir, auxilary)


def find_flares_in_data(data,
        peak_min_width=15,
        peak_min_distance=150,
        rel_height=0.9,
        lc_output_dir='.',
        auxilary=None):
    """
    find flare peak times in data
    Parameters
    -----------
    peak_min_width: int
        the minimal width (number of data points) of a flare, not considered as a flare if the duration is less that limit
    peak_min_distance: int
        min. interval (number of data points) between two peaks 

    """
    if auxilary is None:
        auxilary={'run_id':-1} 

    unix_time = data['time']
    lightcurve = data['lcs'][0]
    lc_smoothed = smooth(lightcurve)
    lc_baseline=get_lightcurve_baseline(lc_smoothed, niter)
    stat = mdb.get_nearest_lc_stats(unix_time[0], max_limit=500)




    conf_set=False
    if lc_baseline is None:
        #only use the lowest lightcurve for flare identification 
        logger.warning('Failed to estimated the baseline, probably due to a too short duration of, trying to use quiet sun period estimation')
        if stat['std'][0] < 2 * math.sqrt(
                stat['median'][0]):  #valid background
            height = stat['median'][0] + 2 * stat['std'][0]
            baseline = stat['median'][0]
            prominence = 2 * stat['std'][0]
            noise_rms = stat['std'][0]
            conf_set=True
    else:
        #med = np.median(lightcurve)
        #these default values are still under testing
        med=np.median(lc_baseline)
        baseline=np.min(lc_baseline)
        height = baseline+ 3 * np.sqrt(baseline) #statistic 
        prominence = 2 * np.sqrt(baseline)
        noise_rms = np.sqrt(baseline)
        baseline = baseline
        conf_set=True
        #height=300

    if not conf_set :
        logger.error('Flare detection failed, no sufficient background data!')
        return

        


    result = {}

    #2 * np.sqrt(np.mean((lc_smoothed-lightcurve)**2))
    print(f'Peak finding parameters: {noise_rms=}, {prominence=}, {peak_min_distance=}, {height=}')

    xpeaks, properties = signal.find_peaks(
            lc_smoothed,
            height=height,
            prominence=prominence,
            width=peak_min_width,
            rel_height=rel_height,  #T90 calculation
            distance=peak_min_distance)
    if xpeaks.size == 0:
        logger.info(f'No peaks found for file {run_id}')
        return 0
    logger.info('Number of peaks:{}'.format(xpeaks.size))
    conditions = {
            'peak_min_width': peak_min_width,
            'peak_min_distance': peak_min_distance,
            'rel_height': rel_height,
            'prominence': prominence,
            'bkg_statistics': stat,
            'min_height': height,
            }
    doc = {}

    #threshold=height

    peak_values = properties['peak_heights']
    peak_unix_times = unix_time[xpeaks]
    peaks_utc = [st.unix2utc(x) for x in peak_unix_times]
    flare_ids = [
            st.unix2datetime(x).strftime("%y%m%d%H%M") for x in peak_unix_times
            ]
    #print('\n'.join(flare_ids))
    #

    threshold = baseline + noise_rms
    flare_start_times, flare_end_times, left_ips, right_ips, boundaries = find_flare_time_ranges(
            unix_time, lc_smoothed, xpeaks, properties, lc_baseline, noise_rms)
    #find flare time ranges

    is_major_flags = find_major_peaks(left_ips, right_ips, peak_values)
    total_counts = [
            int(np.sum(lightcurve[r0:r1])) for r0, r1 in zip(left_ips, right_ips)
            ]
    LC_statistics = []
    if stat:  
        for ipeak in range(xpeaks.size):
            flare_stats = {
                    'bkg_time_range': (stat['start_utc'], stat['end_utc'])
                    }
            upper_bin = 0
            for ilc in range(0, 5):
                flare_stat = {}
                lc_cnts = data['lcs'][ilc][left_ips[ipeak]:right_ips[ipeak]]
                if not np.any(lc_cnts):
                    continue
                num_2sigma = int(
                        np.sum(
                            lc_cnts > stat['median'][ilc] + 2 * stat['std'][ilc]))
                num_3sigma=int(
                        np.sum(
                            lc_cnts > stat['median'][ilc] + 3 * stat['std'][ilc]))
                if num_2sigma > PEAK_MIN_NUM_POINTS:  #longer than half minutes
                    upper_bin = ilc
                flare_stat={
                        'bkg_median': stat['median'][ilc],
                        'bkg_sigma': stat['std'][ilc],
                            'signal_median': int(np.median(lc_cnts)),
                            'signal_max': int(np.max(lc_cnts)),
                            'signal_min': int(np.min(lc_cnts)),
                            'num_points_above_2sigma': num_2sigma, 
                            'num_points_above_3sigma':  num_3sigma} 

                flare_stats['lc' + str(ilc)] = flare_stat

            flare_stats['upper_ilc'] = upper_bin
            LC_statistics.append(flare_stats)

    seconds_per_bin = 4

    durations = np.array([(r[1] - r[0]) * seconds_per_bin
        for r in zip(left_ips, right_ips)])
    cps = total_counts / durations

    total_signal_counts = []
    for i, (cnts, dur) in enumerate(zip(total_counts, durations)):
        sig_cnts = cnts - dur * baseline / seconds_per_bin
        if sig_cnts < 0:
            sig_cnts = cnts - dur * properties['width_heights'][
                    i] / seconds_per_bin
        total_signal_counts.append(sig_cnts)

    peak_width_T90 = np.vstack(
            (unix_time[properties['left_ips'].astype(int)],
                unix_time[properties['right_ips'].astype(int)])).T

    unix_time[properties['left_ips'].astype(int)].tolist()
    flare_end_unix = unix_time[properties['right_ips'].astype(int)].tolist()

    H10_res = signal.peak_widths(lc_smoothed, xpeaks,
            rel_height=0.1)  #H10 calculation
    peak_width_H10 = np.vstack((unix_time[H10_res[2].astype(int)],
        unix_time[H10_res[3].astype(int)])).T

    H70_res = signal.peak_widths(lc_smoothed, xpeaks,
            rel_height=0.7)  #H70 calculation
    peak_width_H70 = np.vstack((unix_time[H70_res[2].astype(int)],
        unix_time[H70_res[3].astype(int)])).T

    H80_res = signal.peak_widths(lc_smoothed, xpeaks,
            rel_height=0.7)  #H70 calculation
    peak_width_H80 = np.vstack((unix_time[H80_res[2].astype(int)],
        unix_time[H80_res[3].astype(int)])).T

    H50_res = signal.peak_widths(lc_smoothed, xpeaks,
            rel_height=0.5)  #H50 calculation
    peak_width_H50 = np.vstack((unix_time[H50_res[2].astype(int)],
        unix_time[H50_res[3].astype(int)])).T
    #calculate peak width at 50% of the maximum count

    doc = {
            'num_peaks': xpeaks.size,
            'peak_unix_time': peak_unix_times.tolist(),
            'peak_counts': peak_values.tolist(),
            'peak_utc': peaks_utc,
            'flare_id': flare_ids,
            'baseline': baseline,
            'threshold': threshold,
            'noise_rms': noise_rms,
            'total_counts': total_counts,
            'total_signal_counts': total_signal_counts,
            'duration': durations.tolist(),
            'mean_cps': cps.tolist(),
            'conditions': conditions,
            'is_boundary': boundaries,
            'peak_width_bins': properties['widths'].tolist(),  #width of height
            'width_height': properties['width_heights'].tolist(
                ),  # height of the width, background level
            'peak_prominence': properties['prominences'].tolist(),
            'time_ranges': {
                'PH70': peak_width_H70.tolist(),
                'PH80': peak_width_H80.tolist(),
                'PH50': peak_width_H50.tolist(),
                'PH10': peak_width_H10.tolist()
                },
            'start_unix': flare_start_times,
            'end_unix': flare_end_times,
            'is_major': is_major_flags,
            'LC_statistics': LC_statistics,
            }
    doc.update(auxilary)

    mdb.save_flare_info(doc)

    #if a flare is not major it will be ignored
    doc['properties'] = properties
    data['lc_smoothed'] = lc_smoothed
    same_plot=True if debug else False
    create_quicklook_plot(data, doc, lc_output_dir, lc_baseline,  same_plot)
    return xpeaks.size


def find_flares_in_files(fid_start, fid_end, img_path='/data/flare_lc'):
    for i in range(fid_start, fid_end + 1):
        print(f'deleting flares of Files {i}')
        mdb.delete_flares_of_file(i)
    for i in range(fid_start, fid_end + 1):
        find_flares_in_one_file(i, lc_output_dir=img_path)


if __name__ == '__main__':
    import sys
    #debug = True
    if len(sys.argv) < 2:
        print('flare_detection file_number')
    elif len(sys.argv) == 2:
        res = find_flares_in_one_file(int(sys.argv[1]),
                lc_output_dir='/data/flare_lc')
        print('Number of peaks:', res)
    else:
        find_flares_in_files(int(sys.argv[1]),
                int(sys.argv[2]),
                img_path='/data/flare_lc')
