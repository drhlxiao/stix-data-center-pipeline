import sys
import os
sys.path.append('.')
from stix.core import mongo_db as db

def analyze(mdb, packets):
    """
        analyze service 5 packets
        mdb: mongodb instance
        packets:  service 5 packets

    """
    for packet in packets:
        header, params=packet['header'], packet['parameters']
        if header['SPID'] == 54352:
            unique_id=int(params[1][1])
            mdb.set_sci_request_failed(unique_id, packet['_id'])

def process_all():
    mdb = db.MongoDB()

    packets=mdb.collection_packets.find({'header.SPID':54352})
    analyze(mdb, packets)


if __name__=='__main__':
    process_all()






