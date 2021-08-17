# author: Hualin Xiao
# pre-process science data, merge bulk science data packets and write merged data to json files
# so that web client side could load the data quickly
import sys
import os
import json
import numpy as np
from datetime import datetime
from stix.core import stix_datatypes as sdt
from stix.spice import stix_datetime
from stix.core import mongo_db as db
from stix.core import stix_logger
from stix.core import config
mdb = db.MongoDB()
logger = stix_logger.get_logger()
level1_products_path = config.get_config(
    'pipeline.daemon.level1_products_path')

DATA_REQUEST_REPORT_SPIDS = [54114, 54115, 54116, 54117, 54143, 54125]
DATA_REQUEST_REPORT_NAME = {
    54114: 'L0',
    54115: 'L1',
    54116: 'L2',
    54117: 'L3',
    54143: 'L4',
    54125: 'ASP'
}

PROCESS_METHODS = {
    'normal': [54114, 54117, 54143, 54125],
    'yield': [54115, 54116]
}




def merge(file_id, remove_existing=True):
    collection = mdb.get_collection_bsd()
    bsd_cursor = collection.find({'run_id': file_id}).sort('_id', 1)
    for doc in bsd_cursor:
        spid = int(doc['SPID'])
        logger.info(f'processing bsd id: {doc["_id"]}, spid:{spid}')
        uid=doc.get('unique_id',None)
        if uid is not None:
            form=mdb.get_bsd_req_form_by_uid(int(uid))
            collection.update_one({'_id': doc['_id']}, {'$set':{'request_form': form}})


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('sci_packet_merger  run_id')
        print('sci_packet_merger  run_id_start id_end')
    elif len(sys.argv)==2:
        merge(int(sys.argv[1]))
    else:
        for i in range(int(sys.argv[1]),int(sys.argv[2])+1):
            print('process:',i)
            merge(i)

