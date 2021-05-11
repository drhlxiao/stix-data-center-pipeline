#!/usr/bin/python3
# request to json
#@Author: Hualin Xiao
#@Date:   Feb. 11, 2020

import sys
from pprint import pprint
from stix.spice import stix_datetime
from datetime import datetime
import json

Request = {'occurrences': []}


def get_uid(start_unix, level):
    start_datetime = stix_datetime.unix2datetime(start_unix)
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


def form_request(start_unix,
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
    uid = get_uid(start_unix, level)
    start_obt = int(stix_datetime.unix2scet(start_unix))
    num_detectors = sum([(detector_mask & (1 << i)) >> i
                         for i in range(0, 32)])
    num_pixels = sum([((pixel_mask & (1 << i)) >> i) for i in range(0, 12)])
    #print(num_detectors, num_pixels)
    T = tmax / tunit
    M = num_detectors
    P = num_pixels
    E = emax - emin + 1
    data_volume = 1.1 * T * (M * P * (E + 4) + 16)
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


def create_request(start_utc,  duration, detector_mask,pixel_mask,  json_filename):
    #start utc : utc str
    unix_t0 = stix_datetime.utc2unix(start_utc)
    json_file = open(json_filename, 'w')
    detector_mask = 0xffffffff - 0x300
    pixel_mask = 0xff  #big pixels
    requests = {
        'occurrences': [
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
    }
    total_volume = 0

    for i in range(0, 16):
        dt = 5760
        unix = unix_t0 + i * dt
        tbin = 20
        emin = 3
        emax = 17
        eunit = 0
        #print(stix_datetime.unix2utc(unix))
        TC = form_request(unix, 1, detector_mask, 0, dt * 10, tbin * 10, emin,
                          emax, eunit, pixel_mask)
        total_volume += TC['data_volume']
        print('{} , {},   {}, {}, {}, {}, {}, {}, {},{} '.format(
            stix_datetime.unix2utc(unix), stix_datetime.unix2utc(unix + dt),
            tbin, detector_mask, pixel_mask, emin, emax, eunit,
            TC['data_volume'], TC['uid']))
        requests['occurrences'].append(TC)
    json_file.write(json.dumps(requests, indent=4))
    print('Operation requests have been written to crab.json!')
    print('data volume: {} KB ({} MB)'.format(total_volume / 1024.,
                                              total_volume / (1024 * 1024)))


main()
