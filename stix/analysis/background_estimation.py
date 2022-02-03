#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Construct STIX QL background level history
    Procedure:
       requesting QL data hour by hour from Mongo database
       Check if STIX is quiet 
       calculate mean, median, rms and save them to database 

"""

import os
import sys
import math
from scipy import signal
import numpy as np
from stix.core import datatypes as sdt
from stix.core import mongo_db as db
from stix.spice import time_utils as st

mdb = db.MongoDB()

MARGIN = 0
FRAME_SPAN = 3600

SPID = 54118
terminal = False


def is_quiet(median, mean, _max, _std):
    if _std[0] < 3 * math.sqrt(median[0]):
        return True
    return False


def process_file(file_id):

    ql_db = mdb.get_collection('quick_look')
    query = {'run_id': file_id, 'SPID': 54118}
    docs = list(ql_db.find(query).sort('start_unix_time', 1))
    if not docs:
        print(f"File {file_id} has no QL LC report")
        return
    start_unix = docs[0]['start_unix_time']
    end_unix = docs[-1]['stop_unix_time']
    start = start_unix
    end = start_unix
    while end < end_unix:
        span = FRAME_SPAN
        end = start + span
        get_QL_background(start, end, file_id)
        start = end


def get_QL_background(start, end, file_id):
    span = end - start
    bkg_db = mdb.get_collection('lc_stats')
    if bkg_db.find({
            'start_unix': {
                '$lt': (start + end) / 2
            },
            'end_unix': {
                '$gt': (start + end) / 2
            }
    }).count() > 0:
        #background exists
        return

    unix_time = []
    packets = mdb.get_LC_pkt_by_tw(start, span)
    last_unix = 0
    start_unix = 0
    lightcurves = {}
    for pkt in packets:
        packet = sdt.Packet(pkt)
        try:
            if not packet.isa(SPID):
                continue
        except:
            continue
        #fig = None
        if packet['header']['unix_time'] < last_unix:
            continue

        last_unix = packet['header']['unix_time']
        if start_unix == 0:
            start_unix = last_unix

        scet_coarse = packet[1].raw
        scet_fine = packet[2].raw
        start_scet = scet_coarse + scet_fine / 65536.
        int_duration = (packet[3].raw + 1) * 0.1
        num_lc_points = packet.get('NIX00270/NIX00271')[0]
        lc = packet.get('NIX00270/NIX00271/*.eng')[0]
        rcr = packet.get('NIX00275/*.raw')
        UTC = packet['header']['UTC']
        num_lcs = len(lc)
        for i in range(len(lc[0])):
            t = st.scet2unix(start_scet + i * int_duration)
            if t > end:
                break
            if t < start:
                continue
            for j in range(num_lcs):
                if j not in lightcurves:
                    lightcurves[j] = []
                lightcurves[j].append(lc[j][i])
            unix_time.append(t)

    if not lightcurves:
        return None
    median = []
    mean = []
    std = []
    amax = []
    amin = []
    num_points = len(unix_time)
    for _id in lightcurves:
        lc = lightcurves[_id]
        if not lc:
            continue
        lca = np.array(lc)
        median.append(np.median(lca))
        mean.append(np.mean(lca))
        std.append(np.std(lca))
        amax.append(max(lc))
        amin.append(min(lc))

    doc = {
        'file_id': file_id,
        'is_quiet': is_quiet(median, mean, max, std),
        'num_points': num_points,
        'start_unix': unix_time[0],
        'end_unix': unix_time[-1],
        'start_utc': st.unix2datetime(unix_time[0]),
        'end_utc': st.unix2datetime(unix_time[-1]),
        'median': median,
        'mean': mean,
        'std': std,
        'max': amax,
        'min': amin,
    }
    mdb.insert_lc_stats(doc)


if __name__ == '__main__':
    import sys
    terminal = True
    if len(sys.argv) == 1:
        print('process file file_number')
    elif len(sys.argv) == 2:
        process_file(int(sys.argv[1]))
    else:
        for i in range(265, 423):
            process_file(i)
