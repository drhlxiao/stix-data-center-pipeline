#!/usr/bin/python
'''
 Extract stix configuration parameters  from telecommand packets and save the information to the database
 Author: Hualin Xiao
 Date: Aug. 31, 2021
'''
import sys

import os
from stix.spice import time_utils
import numpy as np

STIX_CONFIG_TC_NAMES = [
    'ZIX37703',
    'ZIX39019',
    'ZIX39004',
    'ZIX36605',  #ASIC register write, asic latency, and HV
     'ZIX37013',
     'ZIX36605',
    'ZIX37018'
]
'''
TC_GROUPS = {
    'elut': ['ZIX37703', 'ZIX37008'],  #STIX upload ELUT, APPLY LUT, 
    'detector': ['ZIX39019', 'ZIX39004',
                 'ZIX36605'],  #ASIC register write, asic latency, and HV
    'hv': ['ZIX36605'],
     'rot': ['ZIX37013'],
    'asw': ['ZIX37018']
}

'''

ELUT_ENERGY_BINS = [
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36,
    40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150
]


#4:15 keV
def attach_timestamp(header, result):
    utc = header['UTC']
    dt = time_utils.utc2datetime(utc)
    unix = time_utils.utc2unix(utc)
    date = dt.strftime('%Y%m%d')
    result.update({
        'execution_utc': utc,
        'execution_unix': unix,
        'execution_date': date
    })


def parse_elut_upload_telecommand(packet):
    '''extract elut from Telemcommand packets
    '''
    eluts = []
    header = packet['header']
    parameters = packet['parameters']
    if header['TMTC'] != 'TC':
        return []
    if header['name'] != 'ZIX37703':
        return []
    num = parameters[1][1]
    children = parameters[1][3]
    for i in range(0, num):
        i0 = i * 36
        det = children[i0 + 0][1]
        pix = children[i0 + 1][1]
        pixel_id = (det - 1) * 12 + pix
        #if channel_id in ELUT:
        #    continue
        ebins = [children[j + i0 + 2][1] for j in range(0, 32)]
        adcs = np.array(ebins[1:len(ELUT_ENERGY_BINS) + 1]) / 4.
        a, b = np.polyfit(ELUT_ENERGY_BINS, adcs, 1)
        slope = round(a, 5)
        offset = round(b, 5)
        utc = header['UTC']
        elut = {
            'type': 'elut',
            'parameter': 'elut',
            'pixel_id': pixel_id,
            'ebins': ebins,
            'adc': adcs.tolist(),
            'state': header['state'],
            'slope': slope,
            'sub_id':i,
            'offset': offset,
        }
        print(pixel_id, utc)
        attach_timestamp(header, elut)
        eluts.append(elut)

    return eluts


def parse_config_telecommand(packet):
    result = {}
    results=[]
    header = packet['header']
    parameters = packet['parameters']
    if header['name'] == 'ZIX37018':
        result = {
            'type': 'asw',
            'parameter': parameters[0][1],
            'sub_id': parameters[1][1],
            'value': parameters[2][1],
        }
        attach_timestamp(header, result)
        results.append(result)
    elif header['name']=='ZIX37013':
        result = {
            'type': 'asw',
            'parameter': 312,#rotating buffer
            'value':int(parameters[0][1]),
            'state': header['state'],
            'sub_id': 0,
            'category':'rotbuf:min_time_bin',
            'comment':'Updated via ZIX37013'
        }
        #used by the data request command to estimated data
        attach_timestamp(header, result)
        results.append(result)

        result = {
            'type': 'asw',
            'parameter': 313,#rotating buffer
            'value':int(parameters[1][1]) ,
            'sub_id': 0,
            'state': header['state'],
            'category':'rotbuf:max_time_bin',
            'comment':'Updated via ZIX37013'
        }
        #used by the data request command to estimated data
        attach_timestamp(header, result)
        results.append(result)

        result = {
            'type': 'asw',
            'parameter': 314,#rotating buffer
            'value':int(parameters[2][1]) ,
            'state': header['state'],
            'sub_id': 0,
            'category':'rotbuf:threshold',
            'comment':'Updated via ZIX37013'
        }
        #used by the data request command to estimated data
        attach_timestamp(header, result)
        results.append(result)


    elif header['name'] == 'ZIX39004':
        result = {
            'type': 'asic',
            'parameter': 'latency',
            'value': parameters[0][1],
            'sub_id': 0,
            'category':'asic:latency',
            'state': header['state']
        }
        attach_timestamp(header, result)
        results.append(result)
    elif header['name'] == 'ZIX36605':
        #if parameters[1][2] in hv_config:
        #only the last one
        result = {
            'type': 'hvps',
            'parameter': parameters[1][2],
            'value': parameters[2][1],
            'sub_id': 0,
            'category':'hvps:',
            'state': header['state']
        }
        attach_timestamp(header, result)
        results.append(result)

    elif header['name'] == 'ZIX39019':
        num_struct = parameters[0][1]
        children = parameters[0][3]
        for i in range(0, num_struct):
            i0 = i * 48
            detector = children[i0][1]
            #if detector in asic_config:
            #    #only the last one
            #    continue
            reg_offsets = [
                4, 5, 6, 38, 40, 41, 42, 43, 44, 45, 46, 47
            ]  #register 3 threshold, its value is meaningless here
            thr_offsets = [32, 21, 14, 7, 35, 24, 11, 6, 36, 27, 17, 9]  #pixels
            result = {
                'type': 'asic',
                'parameter': 'threshold',
                'detector': detector,
                'register_mask': children[i0 + 2][1],
                'registers': [children[i0 + k][1] for k in reg_offsets],
                'sub_id': i,
                'thresholds': [children[i0 + k][1] for k in thr_offsets],
                'state': header['state'],
            }
            attach_timestamp(header, result)
            results.append(result)

    return results


def parse_telecommand_packet(run_id, packet_id, packet, db_collection=None):
    if not packet:
        return
    header = packet['header']
    if header['TMTC'] != 'TC' or 'state' not in header:
        return
        #skip those TC without confirmation
    if header['name'] not in STIX_CONFIG_TC_NAMES:
        return
    results = []
    if header['name'] == 'ZIX37703':
        results = parse_elut_upload_telecommand(packet)
    else:
        results = parse_config_telecommand(packet)

    if db_collection is not None:
        for doc in results:
            if doc:
                sub_id=doc['sub_id']
                doc.update({'run_id': run_id, 'packet_id': packet_id })
                db_collection.update_one({'packet_id': packet_id, 'sub_id':sub_id},
                                         {'$set': doc},
                                         upsert=True)
                #update or insert if not exist


def process_file(run_id):
    from stix.core import mongo_db as db
    mdb = db.MongoDB()
    print('Processing file:', run_id)
    config_db = mdb.get_collection('config')

    config_db.delete_many({'run_id':run_id})
    #delete existing entries

    packets = mdb.get_telecommands_in_file(run_id)

    for packet in packets:
        run_id = packet['run_id']
        packet_id = packet['_id']
        parse_telecommand_packet(run_id, packet_id, packet, config_db)


if __name__ == '__main__':
    import sys
    terminal = True
    if len(sys.argv) < 2:
        print('usage: instrumment_config  <file_id> [<file_id_end>]')
    elif len(sys.argv) == 2:
        process_file(int(sys.argv[1]))
    else:
        for i in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
            process_file(i)
