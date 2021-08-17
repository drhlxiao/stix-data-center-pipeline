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
from stix.core import stix_datatypes as sdt
from stix.core import mongo_db as db
from stix.spice import stix_datetime
from stix.core import stix_logger
from stix.spice import solo
from stix.analysis import ql_analyzer as qla
logger = stix_logger.get_logger()

mdb = db.MongoDB()

def get_rotating_buffer_config(ut):
    conf={
            'min_tbin':1,
            'max_tbin':20,
            'threshold':1800
            }
    db_pkt=mdb.get_collection('packets')
    conf_pkt=list(db_pkt.find({'header.name':'ZIX37013', 'header.unix_time':{'$lt': ut}}).sort('header.unix_time',-1))
    if not conf_pkt:
        print('Using default configuration because no rotation buffer configuration TC for time: ', ut, ' was found!')
    else:
        parameters=conf_pkt[0]['parameters']
        conf['min_tbin']=(int(parameters[0][1])+1)*0.1
        conf['max_tbin']=(int(parameters[1][1])+1)*0.1
        conf['threshold']=parameters[2][1]
    return conf


def process_file(file_id):
    print("Processing file:",file_id)
    packets = mdb.select_packets_by_run(file_id, SPIDs=54118)
    if not packets:
        return
    data=qla.LightCurveAnalyzer.parse(packets)
    if data is None:
        return 
    unix_time = data['time'] #set to the center of a bin
    lightcurve = data['lcs'][0]+data['lcs'][1] # 4-15 keV
    conf=get_rotating_buffer_config(unix_time[0])
    res={'t':[],'tbin':[], 'counts':[], 'config':conf, 'num_tbins':[]}
    s=0
    last_time = 0
    thr=conf['threshold']
    max_tbin=conf['max_tbin']
    min_tbin=conf['min_tbin']


    for t, c in zip(unix_time,lightcurve):
        t=int(t)
        c=int(c)
        if last_time==0:
            last_time=t
            continue
        if s>0 and (s>=thr or t-last_time>=max_tbin or c>=thr):  #close the opening time bin anyway
            #close last time bin
            tbin=round(min([max_tbin, t-last_time]),1)
            res['tbin'].append(tbin)
            res['num_tbins'].append(1)
            res['t'].append([last_time,t])
            res['counts'].append(s)
            #print(last_time, t,t-last_time,  c, tbin, s)
            s=0
            last_time=t

        if c<thr:
            s+=c
        else:
            num_tbins=math.ceil(c/thr)
            max_bins=math.ceil(4./min_tbin)
            num_tbins=1 if num_tbins<1 else num_tbins
            num_tbins=max_bins if num_tbins>max_bins else num_tbins
            tbin=round(4./num_tbins,1) #time step 0.1s
            res['tbin'].append(tbin)
            res['num_tbins'].append(num_tbins)
            res['t'].append([t-4, t])
            res['counts'].append(c/num_tbins)
            #print(last_time, t, t-last_time, c, tbin, c/num_tbins)
            s=0
            last_time=t


    res['start_unix']=res['t'][0][0]
    res['end_unix']=res['t'][-1][1]
    res['start_utc']=stix_datetime.unix2utc(res['t'][0][0])
    res['end_utc']=stix_datetime.unix2utc(res['t'][-1][1])
    res['file_id']=file_id
    mdb.insert_time_bins(res, delete_existing=True)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('estimate_integration_times file_id_start, file_id_end')
    elif len(sys.argv) >= 2:
        start_id=int(sys.argv[1])
        end_id=start_id
        if len(sys.argv)==3:
            end_id=int(sys.argv[2])
        for i in range(start_id, end_id+1):
            process_file(i)


