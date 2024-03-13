import pymongo
import sys
from stix.core import datatypes as sdt
from stix.spice import time_utils as stu
connect = pymongo.MongoClient('localhost', 9123)
stix = connect['stix']
db = stix['packets']


#'NIXD0001': 'SW Version Number',
#    'NIXD0004': 'IDPU identifier',
params=['NIX00150','NIX00152']
param_list={}

def request_asw_dump_packets(start_utc, end_utc):
    start_unix = stu.utc2unix(start_utc)
    end_unix = stu.utc2unix(end_utc)
    cur = db.find({
        'header.SPID': 54129,
        'header.unix_time': {
            '$gt': start_unix,
            '$lt': end_unix
        }
    }).sort('header.unix_time', 1)
    #print(f'Number of packets:{cur.count()}')
    for pkt in cur:
        header = pkt['header']
        parameters = pkt['parameters']
        if parameters[0][0]!=params[0]:
            raise "Invalid packet selected!"

        _id=parameters[0][1]
        val=parameters[2][1]

        param_list[_id]=val
        
    f=open("param_list.csv","w")
    for key, val in param_list.items():
        hexs=hex(val)
        f.write(f'{key},{val},{hexs}\n')
        print(key)




#request_hk_packets('2020-05-20T00:00:00', '2021-01-01T00:00:00')
request_asw_dump_packets('2023-12-01T00:00:00', '2023-12-10T00:00:00')

