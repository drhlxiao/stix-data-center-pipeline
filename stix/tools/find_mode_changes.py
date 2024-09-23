import pymongo
import sys
sys.path.append('.')
from stix.core import datatypes as sdt
from stix.spice import time_utils as dt
connect = pymongo.MongoClient('localhost', 27017)
stix = connect['stix']
db = stix['packets']
NIX00012={
            "0": "RESET",
            "1": "BOOT",
            "2": "SAFE",
            "3": "MAINT",
            "4": "CONFIG",
            "5": "NOMINAL"
        }

def request_hk_packets(start_utc, end_utc):
    names = HK.keys()
    start_unix = datetime.utc2unix(start_utc)
    end_unix = datetime.utc2unix(end_utc)
    cur = db.find({
        'header.SPID': 54102,
        'header.unix_time': {
            '$gt': start_unix,
            '$lt': end_unix
        }
    }).sort('header.unix_time', 1)
    print(f'Number of packets:{cur.count()}')
    last_time = 0
    csv_filename = f'hk_{start_utc}_{end_utc}.csv'
    csv_filename = csv_filename.replace(':', '')
    fcsv = open(csv_filename, 'w')
    print('request data...')
    i = 0
    indexes = []
    last_mode = 'OFF'


    for pkt in cur:
        header = pkt['header']
        parameters = pkt['parameters']


        if header['unix_time'] <= last_time:
            continue
        this_time = header['unix_time']



        raw_mode = NIX00012[str(parameters[3][1])]

        if this_time - last_time > 300:
            fcsv.write(f'{dt.unix2utc(last_time)}, {dt.unix2utc(this_time)},"NO_HK","{this_time - last_time}"\n')
            print(f'{dt.unix2utc(last_time)}, {dt.unix2utc(this_time)},"NO_HK","{this_time - last_time}"\n')

        last_time = this_time

        if raw_mode == last_mode:
            continue
        else:
            fcsv.write(f'{dt.unix2utc(last_time)}, {dt.unix2utc(this_time)}, {last_mode}, {raw_mode}\n')
            print(f'{dt.unix2utc(last_time)}, {dt.unix2utc(this_time)}, {last_mode}, {raw_mode}\n')


        last_mode = raw_mode

        i+=1

    print('Done.')


request_hk_packets('2020-01-01T00:00:00', '2024-10-04T00:00:00')
