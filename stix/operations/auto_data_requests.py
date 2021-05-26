import sys
sys.path.append('.')
import os
import math
from datetime import datetime
from stix.core import stix_datatypes as sdt
from stix.spice import stix_datetime
from stix.core import mongo_db as db
from stix.core import stix_logger
from stix.core import config
mdb = db.MongoDB()
logger = stix_logger.get_logger()
db=mdb.get_db()
bsd_form=db['bsd_req_forms']
flare_collection=db['flare_tbc']
time_margin=120
def create_request(flare_entry_id, run_id, flare_id,start_utc, duration):
    if list(bsd_form.find().sort(
        '_id', -1)):
        print(f'data request for Flare {flare_id} already exists.')
        return

    try:
        current_id= bsd_form.find().sort(
            '_id', -1).limit(1)[0]['_id'] + 1
    except IndexError:
        current_id=0
    data=create_template(current_id, flare_id, start_utc, duration)
    data['_id']=current_id
    data['flare_entry_id']=flare_entry_id
    data['creation_time']=datetime.now()
    data['run_id']=run_id
    data['status']=0
    data['hidden']=False
    bsd_form.insert_one(data)
    print('Created a request for Flare {flare_id}, ID:{current_id}')

def create_template(flare_id,start_utc, duration, left_margin=0, right_margin=0):

    T = float(duration) 
    M = 32
    P = 12
    emax=21
    emin=0
    E = 22 #22 energy bins
    data_volume=0
    if level == 1:
        data_volume = 1.1 * T * (M * P * (E + 4) + 16)
    if left_margin!=0:
        start_utc=stix_datetime.unix2utc(stix_datetime.utc2unix(start_utc)+left_margin)

    duration=int(duration-left_margin+right_margin)

    return {
    "data_volume" : str(math.floor(data_volume)),
    "execution_date" : '',
    "author" : "Hualin",
    "email" : "stix-obs@fhnw.ch",
    "subject" : str(flare_id),
    "purpose" : "Flare",
    "request_type" : "L1",
    "start_utc" : str(start_utc),
    "duration" : str(duration),
    "time_bin" : "1",
    "detector_mask" : "0xFFFFFFFF",
    "pixel_mask" : "0xFFF",
    "emin" : "0",
    "emax" : "21",
    "eunit" : "1",
    "description" : f"flare {flare_id}, Created by robot",
    "volume" : str(int(data_volume)),
    "unique_ids" : []
        }
def create_data_request_templates(file_ids):
    flares=flare_collection.find({'run_id':{'$in': file_ids}, 'hidden':False})
    for flare in flares:
        flare_id=flare['flare_id']
        start_unix=flare['start_unix']
        start_utc=stix_datetime.unix2utc(start_unix-time_margin)
        duration=flare['duration']+2*time_margin
        create_request(flare_entry_id, file_id, flare_id,start_utc, duration)





    
