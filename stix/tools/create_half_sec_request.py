# author: Hualin Xiao
# pre-process science data, merge bulk science data packets and write merged data to json files
# so that web client side could load the data quickly
import sys
import os
import json
import numpy as np
from datetime import datetime
from stix.core import datatypes as sdt
from stix.spice import datetime
from stix.core import mongo_db as db
from stix.core import logger
from stix.core import config
mdb = db.MongoDB()

def run():
    start_id=4201
    db= mdb.get_collection('data_requests')
    docs= db.find({'start_utc':{'$gt': '2021-09-06T13:11:08'}, 'request_type':'Spectrogram','time_bin':'1'})
    for doc in docs:
        print('find ', doc['_id'])
        doc['_id']=start_id
        doc['data_volume']=0
        doc['subject']='Half sec. request, '+doc['subject']
        doc['data_volume_upper_limit']=0
        doc['time_bin']=0.5
        doc['description']+=f'. Similar to Request # {doc["unique_ids"]} except time bin changed from 1 sec to 0.5 sec'
        doc['unique_ids']=[]
        print(doc)
        start_id+=1
        db.insert_one(doc)




if __name__ == '__main__':
    run()
