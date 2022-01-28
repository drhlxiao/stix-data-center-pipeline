#!/usr/bin/python3
# export request from db to json
#@Author: Hualin Xiao
#@Date:   June. 5, 2020

import sys
import os
import math
import json
from spice import datetime
from core import datatypes as sdt
from datetime import timedelta
from core import mongodb_api
import numpy as np
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
import config
bsd_lc_path= config.get_path('bsd_lcs')
STIX_MDB = mongodb_api.MongoDB()
collection=STIX_MDB.get_collection('bsd_req_forms')

data_levels={'L0':0, 'L1':1,'L2':2,'L3':3,'Spectrogram':4, 'Aspect':5}
MAX_DURATION_PER_REQUEST=6550.
def enable_data_transfer(mask=31399711):
    return {"name":
                "AIXF210A",
                'actionTime':'00:00:00',
                "parameters": [[
                    "XF210A01",
                    mask,
                ]  ]
            }


'''
#def get_uid(unix, duration, request_type, version=0):
    #LEVEL = {'L0': 0, 'L1': 1, 'L2': 2, 'L3': 3, 'Spectrogram': 4, 'Aspect': 5}
    #level = 0
    #try:
    #    level = LEVEL[request_type]
    except KeyError:
    #    level=6

    unix=int(unix)
    start_datetime=datetime.unix2datetime(unix)
    #year = start_datetime.year
    #month = start_datetime.month
    #day = start_datetime.day
    #hour = start_datetime.hour
    #minute = start_datetime.minute
    #if level==5:
    #    return int(f'{start_datetime:%y%m%d%H%M%S}{duration:05}')

    seq_of_day=STIX_MDB.select_user_data_requests_by_date(start_datetime).count()+1
    
    #request_id = (year & 0xF) << 28
    #request_id |= (month & 0xF) << 24
    #request_id |= (day & 0x1F) << 19
    #request_id |= (hour & 0x1F) << 14
    #request_id |= (minute & 0x3F) << 8
    #request_id |= (level & 0xF) << 4
    #request_id += version
    #while not STIX_MDB.select_user_data_request_forms_by_uid(request_id):
    #    request_id+=1
    #    print(request_id)
    request_id=f'{start_datetime.strftime("%y%m%d")}{seq_of_day:04}'
    return request_id

    '''


def form_aspect_request(start_unix, stop_unix,summing):
    start_obt = int(datetime.unix2scet(start_unix))
    stop_obt = int(datetime.unix2scet(stop_unix))
    duration=stop_obt-start_obt
    tm_packets=math.ceil(8*(64/summing)*duration/4096.)

    volume=8*(64/summing)*duration+26*tm_packets
    return {"name":
                "AIXF422A",
                'actionTime':'00:10:00',
                'data_volume':volume,
                "parameters": [[
                    "XF422A01",
                    start_obt,
                ], [
                    "XF422A02",
                    stop_obt,
                ], ["XF422A03", summing]]
            }


def form_bsd_request_sequence(uid, start_unix,
                 level,
                 detector_mask,
                 tmin,
                 tmax,
                 tunit,
                 emin,
                 emax,
                 eunit,
                 pixel_mask=0xfff,
                 action_time='00:00:10'):

    start_obt = int(datetime.unix2scet(start_unix))
    num_detectors = sum([(detector_mask & (1 << i)) >> i
                         for i in range(0, 32)])
    num_pixels = sum([((pixel_mask & (1 << i)) >> i) for i in range(0, 12)])
    #print(num_detectors, num_pixels)
    T = tmax / tunit
    M = num_detectors
    P = num_pixels
    E = (emax - emin + 1)/(eunit+1)
    data_volume=0
    if level == 1:
        data_volume = 1.1 * T * (M * P * (E + 4) + 16)
    elif level == 4:
        data_volume = 1.1 * T * (E + 4)
        

    parameters = [
        ['XF417A01', uid],
        ['XF417A02', level],
        ['XF417A03', start_obt],
        ['XF417A04', 0],  #subseconds
        ['XF417A05', "0x%X" % detector_mask],
        #['PIX00071', 1],  #number of structure
        ['XF417A06', tmin],
        ['XF417A07', tmax],
        ['XF417A08', tunit],
        ['XF417A09', emin],
        ['XF417A10', emax],
        ['XF417A11', eunit]
    ]
    return {
        'name': 'AIXF417A',
        'actionTime': action_time,
        'uid': uid,
        'data_volume': data_volume,
        'parameters': parameters
    }


def form_mask_config_telecommands(detector_mask, pixel_mask, level):
    """  #index, 
    #see https://docs.google.com/spreadsheets/d/1xRrmUTjNuLFie8NlmUF9fU7Ncomk5KwBJQwoKgQUgNY/edit?usp=sharing
    323	DETECTOR_MASK_L0	PARAM_INT32_MAX	0	PARAM_INT32_MAX
324	DETECTOR_MASK_L1	PARAM_INT32_MAX	0	PARAM_INT32_MAX
325	DETECTOR_MASK_L2	PARAM_INT32_MAX	0	PARAM_INT32_MAX
326	DETECTOR_MASK_L3	PARAM_INT32_MAX	0	PARAM_INT32_MAX
327	DETECTOR_MASK_L4	PARAM_INT32_MAX	0	PARAM_INT32_MAX328	PIXEL_MASK_L0	0xFFF	0	0xFFF
329	PIXEL_MASK_L1R0	0xFFF	0	0xFFF
330	PIXEL_MASK_L1R1	0xFFF	0	0xFFF
331	PIXEL_MASK_L1R2	0xFFF	0	0xFFF
332	PIXEL_MASK_L1R3	0xFFF	0	0xFFF
333	PIXEL_MASK_L1R4	0xFFF	0	0xFFF
334	PIXEL_MASK_L1R5	0xF	0	0xFFF
335	PIXEL_MASK_L1R6	0xF	0	0xFFF
336	PIXEL_MASK_L1R7	0xF	0	0xFFF
337	PIXEL_MASK_L2R0P1	0x888	0	0xFFF
338	PIXEL_MASK_L2R0P2	0x444	0	0xFFF
339	PIXEL_MASK_L2R0P3	0x222	0	0xFFF
340	PIXEL_MASK_L2R0P4	0x111	0	0xFFF
341	PIXEL_MASK_L2R0P5	0	0	0xFFF
342	PIXEL_MASK_L2R1P1	0x888	0	0xFFF
343	PIXEL_MASK_L2R1P2	0x444	0	0xFFF
344	PIXEL_MASK_L2R1P3	0x222	0	0xFFF
345	PIXEL_MASK_L2R1P4	0x111	0	0xFFF
346	PIXEL_MASK_L2R1P5	0	0	0xFFF
347	PIXEL_MASK_L2R2P1	0x888	0	0xFFF
348	PIXEL_MASK_L2R2P2	0x444	0	0xFFF
349	PIXEL_MASK_L2R2P3	0x222	0	0xFFF
350	PIXEL_MASK_L2R2P4	0x111	0	0xFFF
351	PIXEL_MASK_L2R2P5	0	0	0xFFF
352	PIXEL_MASK_L2R3P1	0x888	0	0xFFF
353	PIXEL_MASK_L2R3P2	0x444	0	0xFFF
354	PIXEL_MASK_L2R3P3	0x222	0	0xFFF
355	PIXEL_MASK_L2R3P4	0x111	0	0xFFF
356	PIXEL_MASK_L2R3P5	0	0	0xFFF
357	PIXEL_MASK_L2R4P1	0x8	0	0xFFF
358	PIXEL_MASK_L2R4P2	0x4	0	0xFFF
359	PIXEL_MASK_L2R4P3	0x2	0	0xFFF
360	PIXEL_MASK_L2R4P4	0x1	0	0xFFF
361	PIXEL_MASK_L2R4P5	0	0	0xFFF
362	PIXEL_MASK_L2R5P1	0x8	0	0xFFF
363	PIXEL_MASK_L2R5P2	0x4	0	0xFFF
364	PIXEL_MASK_L2R5P3	0x2	0	0xFFF
365	PIXEL_MASK_L2R5P4	0x1	0	0xFFF
366	PIXEL_MASK_L2R5P5	0	0	0xFFF
367	PIXEL_MASK_L2R6P1	0x8	0	0xFFF
368	PIXEL_MASK_L2R6P2	0x4	0	0xFFF
369	PIXEL_MASK_L2R6P3	0x2	0	0xFFF
370	PIXEL_MASK_L2R6P4	0x1	0	0xFFF
371	PIXEL_MASK_L2R6P5	0	0	0xFFF
372	PIXEL_MASK_L2R7P1	0x8	0	0xFFF
373	PIXEL_MASK_L2R7P2	0x4	0	0xFFF
374	PIXEL_MASK_L2R7P3	0x2	0	0xFFF
375	PIXEL_MASK_L2R7P4	0x1	0	0xFFF
376	PIXEL_MASK_L2R7P5	0	0	0xFFF
377	PIXEL_MASK_L3R0P1	0x888	0	0xFFF
378	PIXEL_MASK_L3R0P2	0x444	0	0xFFF
379	PIXEL_MASK_L3R0P3	0x222	0	0xFFF
380	PIXEL_MASK_L3R0P4	0x111	0	0xFFF
381	PIXEL_MASK_L3R0P5	0xFFF	0	0xFFF
382	PIXEL_MASK_L3R1P1	0x888	0	0xFFF
383	PIXEL_MASK_L3R1P2	0x444	0	0xFFF
384	PIXEL_MASK_L3R1P3	0x222	0	0xFFF
385	PIXEL_MASK_L3R1P4	0x111	0	0xFFF
386	PIXEL_MASK_L3R1P5	0xFFF	0	0xFFF
387	PIXEL_MASK_L3R2P1	0x888	0	0xFFF
388	PIXEL_MASK_L3R2P2	0x444	0	0xFFF
389	PIXEL_MASK_L3R2P3	0x222	0	0xFFF
390	PIXEL_MASK_L3R2P4	0x111	0	0xFFF
391	PIXEL_MASK_L3R2P5	0xFFF	0	0xFFF
392	PIXEL_MASK_L3R3P1	0x888	0	0xFFF
393	PIXEL_MASK_L3R3P2	0x444	0	0xFFF
394	PIXEL_MASK_L3R3P3	0x222	0	0xFFF
395	PIXEL_MASK_L3R3P4	0x111	0	0xFFF
396	PIXEL_MASK_L3R3P5	0xFFF	0	0xFFF
397	PIXEL_MASK_L3R4P1	0x888	0	0xFFF
398	PIXEL_MASK_L3R4P2	0x444	0	0xFFF
399	PIXEL_MASK_L3R4P3	0x222	0	0xFFF
400	PIXEL_MASK_L3R4P4	0x111	0	0xFFF
401	PIXEL_MASK_L3R4P5	0xFFF	0	0xFFF
402	PIXEL_MASK_L3R5P1	0x8	0	0xFFF
403	PIXEL_MASK_L3R5P2	0x4	0	0xFFF
404	PIXEL_MASK_L3R5P3	0x2	0	0xFFF
405	PIXEL_MASK_L3R5P4	0x1	0	0xFFF
406	PIXEL_MASK_L3R5P5	0xF	0	0xFFF
407	PIXEL_MASK_L3R6P1	0x8	0	0xFFF
408	PIXEL_MASK_L3R6P2	0x4	0	0xFFF
409	PIXEL_MASK_L3R6P3	0x2	0	0xFFF
410	PIXEL_MASK_L3R6P4	0x1	0	0xFFF
411	PIXEL_MASK_L3R6P5	0xF	0	0xFFF
412	PIXEL_MASK_L3R7P1	0x8	0	0xFFF
413	PIXEL_MASK_L3R7P2	0x4	0	0xFFF
414	PIXEL_MASK_L3R7P3	0x2	0	0xFFF
415	PIXEL_MASK_L3R7P4	0x1	0	0xFFF
416	PIXEL_MASK_L3R7P5	0xF	0	0xFFF
"""
    dmask_map={4: '424', 1: '323', 2:'324', 3:'236',0:'322'}  
    pmask_map={4: '425', 1: '328', 2:'336', 3:'377',0:'327'}  
    #see https://docs.google.com/spreadsheets/d/1xRrmUTjNuLFie8NlmUF9fU7Ncomk5KwBJQwoKgQUgNY/edit?usp=sharing
    #note this has to be updated according to RCR level
    return [{
                "name":
                "AIXF414A",
                'actionTime':'08:00:00',
                "parameters": [[
                    "XF414A01",
                    dmask_map[level],
                ], [
                    "XF414A02",
                    "0",
                ], ["XF414A03", "0x%X" % detector_mask]]
            },
            {
                "name":
                "AIXF414A",
                'actionTime':'00:00:01',
                "parameters": [[
                    "XF414A01",
                    pmask_map[level],
                ], [
                    "XF414A02",
                    "0",
                ], ["XF414A03", "0x%X" % pixel_mask]]
            },
        ]

def parse_int(text):
    if '0x' in text:
        return int(text,0)
    return int(text)


def attach_TC_aux_info(TC, form):
    TC.update({
        'author':form['author'],
        'subject':form['subject'],
        'db_id':form['_id'],
        })

def create_occurrences(collection, _ids, forms):
    last_detector_mask=0
    last_pixel_mask=0
    total_volume=0
    last_level=-1;

    requests = { 'occurrences':  []  }
    TC_enabled_transfer=enable_data_transfer()
    requests['occurrences'].append(TC_enabled_transfer)

    #seq_of_day=STIX_MDB.select_user_data_requests_by_date(start_datetime).count()+1
    next_uids={}
    for form in  forms:
        print('_id:',form['_id'])
        start_utc = form['start_utc']
        start_unix= datetime.utc2unix(start_utc)
        level=data_levels[form['request_type']]
        dt = int(form['duration'])
        unique_ids=[]
        #lc_filename=make_lightcurve(form['_id'], start_unix, dt)
        #if lc_filename:
        #    form['lc_filename']=lc_filename

        start_date=datetime.utc2datetime(start_utc)
        start_date_str=start_date.strftime('%y%m%d')
        if start_date_str not in next_uids:
            next_uids[start_date_str]=STIX_MDB.get_user_data_request_next_unique_id(start_date, _ids)


        if level==5:
            #aspect data
            TC=form_aspect_request(start_unix, start_unix+dt, int(form['averaging']))
            TC['uid']=next_uids[start_date_str]
            attach_TC_aux_info(TC, form)
            unique_ids=[next_uids[start_date_str]]
            next_uids[start_date_str]+=1
            requests['occurrences'].append(TC)
            total_volume += TC['data_volume']
        else:
            #science data
            mask_TCs=[]
            detector_mask=parse_int(form['detector_mask'])
            pixel_mask= parse_int(form['pixel_mask'])
            tbin=int(form['time_bin'])
            emin=int(form['emin'])
            emax=int(form['emax'])
            eunit=int(form['eunit'])
            if detector_mask!=last_detector_mask or pixel_mask!=last_pixel_mask or level!=last_level:
                mask_TCs=form_mask_config_telecommands(detector_mask, pixel_mask, level )
                #rcr 
                requests['occurrences'].extend(mask_TCs)
            num_TCs=math.ceil(dt/MAX_DURATION_PER_REQUEST)
            #create several TCs for long requests
            if num_TCs<1:
                num_TCs=1
            last_id=0
            for i in range(0,num_TCs):
                T0=start_unix+i*MAX_DURATION_PER_REQUEST
                deltaT=dt-i*MAX_DURATION_PER_REQUEST
                if deltaT>MAX_DURATION_PER_REQUEST:
                    deltaT=MAX_DURATION_PER_REQUEST
                uid=next_uids[start_date_str]
                next_uids[start_date_str]+=1
                #get_uid(T0, deltaT, form['request_type'], 0)
                TC = form_bsd_request_sequence(uid, T0, level, detector_mask, 0, deltaT * 10, tbin * 10, emin,
                                  emax, eunit-1, pixel_mask)
                unique_ids.append(uid)
                attach_TC_aux_info(TC, form)
                requests['occurrences'].append(TC)
                total_volume += TC['data_volume']
            last_pixel_mask=pixel_mask
            last_detector_mask=detector_mask
            last_level=level
        form['unique_ids']=unique_ids
        form['export_time']=datetime.get_now()
        collection.save(form)


    requests['total_volume']=total_volume
    return requests

if __name__=='__main__':
    id_start = 830
    id_end = 1450
    flares=[
                '2105211927',
                '2105130858',
                '2105122233',
                '2105110325',
                '2105091355',
                '2105090434',
                '2105081840',
                '2105071900',
                '2105070036',
                '2105052229',
                '2105040617',
                '2105040506',
                '2105021215',
                '2105021150',
                '2105021125'
     ]
    subjects=[f'Flare # {x}' for x in flares]
    ids=[830,831,832, 1472, 1471, 1473,1474]
    aspects=range(1030,1043)
    ids.extend(aspects)

    prio_ids=[]
    all_forms=[]
    forms1= collection.find({'$or':[{'_id': {'$in':ids}}, {'subject': {'$in':subjects}}]}).sort([('request_type',-1),('detector_mask',1),('pixel_mask',-1)])
    for f in forms1:
        if f['hidden']==True:
            continue
        #print(f['_id'])
        prio_ids.append(f['_id'])
        #all_forms.append(f)
    #candiate_ids=[]
    forms2= collection.find({'_id': {'$gte':1192}}).sort([('request_type',-1),('detector_mask',1),('pixel_mask',-1)])
    for f in forms2:
        if f['hidden']:
            continue
        if f['_id'] not in prio_ids:
            prio_ids.append(f['_id'])
            if len(prio_ids)>140:
                break
        #prio_ids.append(f['_id'])
        #all_forms.append(f)

    forms3= collection.find({'_id': {'$in':prio_ids}}).sort([('request_type',-1),('detector_mask',1),('pixel_mask',-1)])
    requests=create_occurrences(collection,prio_ids, forms3)
    print(requests)
    print('number of requests:', len(requests['occurrences']) )
    print('volume:', requests['total_volume']/(1024*1024.) )
    with open('STP151_data_request.json', 'w') as txtfile:
        json.dump(requests, txtfile)
    

