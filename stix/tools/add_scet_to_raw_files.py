import sys
import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
raw_db= db['raw_files']
packets_db= db['packets']

from stix.core import decompressor
from stix.spice import datetime
from stix.core import datatypes as std

def get_packet(file_id):
    start_scet = 0
    end_scet = 0
    for pkt in packets_db.find({'run_id':file_id}).sort('_id',1):
        start_scet=pkt['header'].get('SCET',0)
        if start_scet>0:
            break
    for pkt in packets_db.find({'run_id':file_id}).sort('_id',-1).limit(1):
        end_scet=pkt['header'].get('SCET',0)
        if end_scet>0:
            break
    return start_scet, end_scet

        
def main(start_id, end_id):
    cursor=raw_db.find({'hidden':False, '_id':{'$gte': start_id, '$lte':end_id}})
    for raw in cursor:
        file_id=int(raw['_id'])
        start_scet, end_scet=get_packet(file_id)
        raw['data_start_scet']=start_scet
        raw['data_end_scet']=end_scet
        print(file_id, start_scet, end_scet)
        raw_db.update({'_id':raw['_id']}, raw)
                


    

    
    
