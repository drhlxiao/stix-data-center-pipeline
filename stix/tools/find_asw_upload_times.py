import pymongo
import sys
from stix.core import datatypes as sdt
from stix.spice import time_utils as stu
connect = pymongo.MongoClient('localhost', 27017)
stix = connect['stix']
db = stix['packets']


#'NIXD0001': 'SW Version Number',
#    'NIXD0004': 'IDPU identifier',
params=['NIXD0001','NIXD0004']

def request_hk_packets(start_utc, end_utc):
    start_unix = stu.utc2unix(start_utc)
    end_unix = stu.utc2unix(end_utc)
    cur = db.find({
        'header.SPID': 54102,
        'header.unix_time': {
            '$gt': start_unix,
            '$lt': end_unix
        }
    }).sort('header.unix_time', 1)
    print(f'Number of packets:{cur.count()}')
    last_time = 0
    fmain= open('main.csv', 'w')
    fred= open('red.csv', 'w')
    i = 0
    irow=0
    indexes = []
    last_idpu=0
    last_asw=0
    last_asw_red=0
    for pkt in cur:
        header = pkt['header']
        parameters = pkt['parameters']
        if header['unix_time'] <= last_time:
            continue
        irow+=1
        if parameters[47][0]!=params[0] and parameters[53][0]!=params[1]:
            raise TypeError("packet format may changed")
            #version, ipdu
        this_asw=parameters[47][1]
        this_idpu=parameters[53][1]
        if this_asw!=last_asw or this_idpu!=last_idpu:
            last_asw=this_asw
            last_idpu=this_idpu
            line=f'{header.utc}, {this_asw}, {this_idpu}\n'
            print(line)
            fcsv.write(line)



    print('Done.')


#request_hk_packets('2020-05-20T00:00:00', '2021-01-01T00:00:00')
request_hk_packets('2020-04-15T00:00:00', '2022-08-04T00:00:00')

