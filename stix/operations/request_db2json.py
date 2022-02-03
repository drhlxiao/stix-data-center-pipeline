#!/usr/bin/python3
# export request from db to json
#@Author: Hualin Xiao
#@Date:   Feb. 11, 2020

import sys
from pprint import pprint
from stix.spice import datetime
from datetime import datetime
import pymongo
import json
import math

data_levels = {
    'L0': 0,
    'L1': 1,
    'L2': 2,
    'L3': 3,
    'Spectrogram': 4,
    'Aspect': 5
}
MAX_DURATION_PER_REQUEST = 6550.


def get_uid(start_unix, level):
    start_datetime = datetime.unix2datetime(start_unix)
    year = start_datetime.year
    month = start_datetime.month
    day = start_datetime.day
    hour = start_datetime.hour
    minute = start_datetime.minute
    request_id = (year & 0xF) << 28
    request_id |= (month & 0xF) << 24
    request_id |= (day & 0x1F) << 19
    request_id |= (hour & 0x1F) << 14
    request_id |= (minute & 0x3F) << 8
    request_id |= (level & 0xF) << 4
    return request_id


def form_aspect_request(start_unix, stop_unix, summing):
    start_obt = int(datetime.unix2scet(start_unix))
    stop_obt = int(datetime.unix2scet(stop_unix))
    duration = stop_obt - start_obt
    tm_packets = math.ceil(8 * (64 / summing) * duration / 4096.)

    volume = 8 * (64 / summing) * duration + 26 * tm_packets
    return {
        "name":
        "AIXF422A",
        'actionTime':
        '00:02:00',
        'data_volume':
        volume,
        "parameters": [[
            "XF422A01",
            start_obt,
        ], [
            "XF422A02",
            stop_obt,
        ], ["XF422A03", summing]]
    }


def form_bsd_request(uid,
                     start_unix,
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
    #uid = get_uid(start_unix, level)
    start_obt = int(datetime.unix2scet(start_unix))
    num_detectors = sum([(detector_mask & (1 << i)) >> i
                         for i in range(0, 32)])
    num_pixels = sum([((pixel_mask & (1 << i)) >> i) for i in range(0, 12)])
    #print(num_detectors, num_pixels)
    T = tmax / tunit
    M = num_detectors
    P = num_pixels
    E = (emax - emin + 1) / (eunit + 1)
    data_volume = 0
    if level == 1:
        data_volume = 1.1 * T * (M * P * (E + 4) + 16)
    elif level == 4:
        data_volume = 1.1 * T * (E + 4)

    parameters = [
        ['PIX00076', uid],
        ['PIX00070', level],
        ['PIX00072', start_obt],
        ['PIX00073', 0],  #subseconds
        ['PIX00248', "0x%X" % detector_mask],
        ['PIX00071', 1],  #number of structure
        ['PIX00077', tmin],
        ['PIX00078', tmax],
        ['PIXX0079', tunit],
        ['PIX00200', emin],
        ['PIX00201', emax],
        ['PIX00202', eunit]
    ]
    return {
        'name': 'ZIXX3801',
        'actionTime': action_time,
        'uid': uid,
        'data_volume': data_volume,
        'parameters': parameters
    }


def form_bsd_request_sequence(uid,
                              start_unix,
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
    #uid = get_uid(start_unix, level)
    start_obt = int(datetime.unix2scet(start_unix))
    num_detectors = sum([(detector_mask & (1 << i)) >> i
                         for i in range(0, 32)])
    num_pixels = sum([((pixel_mask & (1 << i)) >> i) for i in range(0, 12)])
    #print(num_detectors, num_pixels)
    T = tmax / tunit
    M = num_detectors
    P = num_pixels
    E = (emax - emin + 1) / (eunit + 1)
    data_volume = 0
    if level == 1:
        data_volume = 1.1 * T * (M * P * (E + 4) + 16)
    elif level == 4:
        data_volume = 1.1 * T * (E + 4)

    parameters = [
        ['XF417B01', uid],
        ['XF417B02', level],
        ['XF417B03', start_obt],
        ['XF417B04', 0],  #subseconds
        ['XF417B05', "0x%X" % detector_mask],
        #['PIX00071', 1],  #number of structure
        ['XF417B06', tmin],
        ['XF417B07', tmax],
        ['XF417B08', tunit],
        ['XF417B09', emin],
        ['XF417B10', emax],
        ['XF417B11', eunit]
    ]
    return {
        'name': 'AIXF417B',
        'actionTime': action_time,
        'uid': uid,
        'data_volume': data_volume,
        'parameters': parameters
    }


def form_mask_config_telecommands(detector_mask, pixel_mask):
    return [
        {
            "name":
            "AIXF414A",
            "parameters": [[
                "XF414A01",
                "424",
            ], [
                "XF414A02",
                "0",
            ], ["XF414A03", "0x%X" % detector_mask]]
        },
        {
            "name":
            "AIXF414A",
            "parameters": [[
                "XF414A01",
                "425",
            ], [
                "XF414A02",
                "0",
            ], ["XF414A03", "0x%X" % pixel_mask]]
        },
    ]


def parse_int(text):
    if '0x' in text:
        return int(text, 0)
    return int(text)


def main(_ids, json_filename, server='localhost', port=27017):
    try:
        connect = pymongo.MongoClient(server, port)
        db = connect["stix"]
        collection = db['data_requests']
    except Exception as e:
        print('request')
        print(str(e))

    print('request2')
    requst_forms = list(collection.find({'_id': {'$in': _ids}}))

    last_detector_mask = 0
    last_pixel_mask = 0
    total_volume = 0

    msg = ['Author, subject, request_id, data level/type, data volume (kB) ']
    with open(json_filename, 'w') as json_file:
        requests = {'occurrences': []}
        for form in requst_forms:
            start_utc = form['start_utc']
            start_unix = datetime.utc2unix(start_utc)
            level = data_levels[form['request_type']]
            dt = int(form['duration'])
            subject = form['subject']
            author = form['author']
            if level == 5:
                TC = form_aspect_request(start_unix, start_unix + dt,
                                         int(form['averaging']))
                TC['author'] = author
                TC['subject'] = subject

                TC['uid'] = None
                requests['occurrences'].append(TC)
                total_volume += TC['data_volume']

            else:
                mask_TCs = []
                detector_mask = parse_int(form['detector_mask'])
                pixel_mask = parse_int(form['pixel_mask'])
                tbin = int(form['time_bin'])
                emin = int(form['emin'])
                emax = int(form['emax'])
                eunit = int(form['eunit'])
                if detector_mask != last_detector_mask or pixel_mask != last_pixel_mask:
                    mask_TCs = form_mask_config_telecommands(
                        detector_mask, pixel_mask)
                    requests['occurrences'].extend(mask_TCs)

                num_TCs = math.ceil(dt / MAX_DURATION_PER_REQUEST)

                unique_ids = []

                for i in range(0, num_TCs):
                    T0 = start_unix + i * MAX_DURATION_PER_REQUEST
                    deltaT = dt - i * MAX_DURATION_PER_REQUEST
                    if deltaT > MAX_DURATION_PER_REQUEST:
                        deltaT = MAX_DURATION_PER_REQUEST
                    uid = form['unique_id']
                    TC = form_bsd_request_sequence(uid, T0, level,
                                                   detector_mask, 0,
                                                   deltaT * 10, tbin * 10,
                                                   emin, emax, eunit - 1,
                                                   pixel_mask)
                    unique_ids.append(TC['uid'])
                    TC['author'] = author
                    TC['subject'] = subject

                    requests['occurrences'].append(TC)
                    total_volume += TC['data_volume']
                    #form['unique_id']=' '.join(unique_ids)
                last_pixel_mask = pixel_mask
                last_detector_mask = detector_mask

            msg.append('{},{},{}, {},{}'.format(author, subject, level,
                                                TC['uid'],
                                                TC['data_volume'] / 1024.))

        json_file.write(json.dumps(requests, indent=4))

        print('Operation requests have been written to ', json_filename)
        print('data volume: {} KB ({} MB)'.format(
            total_volume / 1024., total_volume / (1024 * 1024)))
        print('\n'.join(msg))


if __name__ == '__main__':
    port = 9000
    if len(sys.argv) >= 4:
        ids = [x for x in range(int(sys.argv[1]), int(sys.argv[2]))]
        if len(sys.argv) == 5:
            port = int(sys.argv[4])

        print(ids)
        main(ids, sys.argv[3], port=port)
    else:
        print(' run start_id, stop_id, filename.json, port')
