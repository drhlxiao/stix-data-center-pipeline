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
DATA_LEVEL_NAMES={0: 'L0',1:'L1',2: 'L2',3: 'L3',4: 'Spectrogram',5: 'Aspect'}
EMAX_MAP=[8, 13, 17, 23, 28] # emax if there are prominence the corresponding lightcurve

logger = stix_logger.get_logger()
db=mdb.get_db()
bsd_form=db['bsd_req_forms']
flare_collection=db['flares_tbc']
L4_TIME_MARGIN=120
L1_TIME_MARGIN=300

L4_REQUEST_GROUP_MAX_TIME_GAP=3600

def create_requests(flare_docs):
    flares=list(flare_docs)
    if len(flares)==0:
        print('No flares selected')
        return
    for doc in flares:
        create_l1_request(doc)
    create_l4_groups(flares)

def create_l4_groups(flare_docs):
    group=[]
    groups=[]
    last_time=flare_docs[0]['peak_unix_time']
    for doc in flare_docs:
        this_time=doc['peak_unix_time']
        if this_time-last_time>L4_REQUEST_GROUP_MAX_TIME_GAP:
            groups.append(group)
            group=[]
        last_time=this_time
        group.append(doc)
    print('number of L4 groups', len(groups))
    

    for gp in groups:
        start_unix=gp[0]['start_unix']
        end_unix=gp[-1]['end_unix']
        flare_ids=[d['flare_id'] for d in gp]
        start_utc=stix_datetime.unix2utc(start_unix)
        duration=end_unix-start_unix
        emax=13
        flare_entry_ids=[x['_id'] for x in gp]
        run_ids=[x['run_id'] for x in gp]
        try:
            ilc=max([d['LC_statistics']['upper_ilc']  for d in gp])
            emax=EMAX_MAP[ilc]
            print(f'Max energy bin changed to {emax}')
        except:
            pass
        create_template(flare_ids, flare_entry_ids, run_ids, start_utc, 
                duration, emax, left_margin=-L4_TIME_MARGIN, right_margin=L4_TIME_MARGIN, level=4)

        




def create_l1_request(doc):
    flare_id=doc['flare_id']
    flare_entry_id=doc['_id']
    start_utc=stix_datetime.unix2utc(doc['start_unix'])
    duration=int(doc['duration'])
    flare_entry_ids=doc['_id']
    run_ids=doc['run_id']
    emax=13
    try:
        ilc=doc['LC_statistics']['upper_ilc']
        emax=EMAX_MAP[ilc]
    except:
        pass
    create_template(flare_id, flare_entry_ids, run_ids, start_utc, duration, emax, -L1_TIME_MARGIN, L1_TIME_MARGIN, level=1)

def create_template(flare_ids, flare_entry_ids, run_ids, start_utc, duration, emax=13, left_margin=0, right_margin=0,
        level=1):
    level_name=DATA_LEVEL_NAMES[level]
    if list(bsd_form.find({'flare_id':flare_ids, 'request_type':level_name}).sort(
        '_id', -1)):
        print(f'data request for Flare {flare_ids} already exists.')
        return
    try:
        current_id= bsd_form.find().sort(
            '_id', -1).limit(1)[0]['_id'] + 1
    except IndexError:
        current_id=0

    if level not in [1,4]:
        print('Not supported data level')
        return

    if left_margin!=0:
        start_utc=stix_datetime.unix2utc(stix_datetime.utc2unix(start_utc)+left_margin)

    if isinstance(flare_ids, list):
        if len(flare_ids)==1:
            flare_ids=flare_ids[0]

    duration=int(duration-left_margin+right_margin)
    detector_mask='0xFFFFFFFF' if level==1 else '0xFFFFFCFF'
    num_detectors = 32 if level==1 else 30
    emin=0
    eunit=1
    num_pixels = 12
    tunit=1
    tmax=duration
    T = tmax / tunit
    M = num_detectors
    P = num_pixels
    E = (emax - emin + 1)/(eunit+1)
    data_volume=0
    if level == 1:
        data_volume = 1.1 * T * (M * P * (E + 4) + 16)
    elif level == 4:
        data_volume = 1.1 * T * (E + 4)

    form={
    "data_volume" : str(math.floor(data_volume)),
    "execution_date" : '',
    "author" : "Hualin",
    "email" : "stix-obs@fhnw.ch",
    "subject" :  f"Flare # {flare_ids}",
    "purpose" : "Solar Flare",
    "request_type" : level_name,
    "start_utc" : str(start_utc),
    "duration" : str(duration),
    "time_bin" : "1",
    "flare_id" : flare_ids,
    'flare_entry_ids':flare_entry_ids,
    "detector_mask" : detector_mask, 
    "creation_time": datetime.utcnow(),
    "pixel_mask" : "0xFFF",
    "emin" : str(emin),
    "emax" : str(emax),
    'hidden':False,
    'run_id': run_ids,
    'status':0,
    "eunit" : str(eunit),
    '_id':current_id, 
    "description" : f"{level_name} request for Flare {flare_ids}",
    "volume" : str(int(data_volume)),
    "unique_ids" : []
        }
    print(f'Inserting request {form["_id"]}, type: {level_name}')
    bsd_form.insert_one(form)
#def create_data_request_for_files(file_ids, cnts_threshold=400):
#    flares=flare_collection.find({'run_id':{'$in': file_ids}, 'peak_counts':{'$gt':cnts_threshold}, 'hidden':False}).sort('peak_unix_time',1)
#    create_requests(flares)

def create_data_request_for_flares(start_flare_id, end_flare_id, cnts_threshold=400):
    query={'_id':{'$gte': start_flare_id, '$lte':end_flare_id},  'peak_counts':{'$gt':cnts_threshold}, 'hidden':False}
    flares=flare_collection.find(query).sort('peak_unix_time',1)
    print('counts',flares.count())
    create_requests(flares)



if __name__=='__main__':
    if len(sys.argv)<3:
        print('auto_req <begin_flare_id> <end_flare_id> [threshold:400]')
    else:
        threshold=400
        if len(sys.argv)==4:
            threshold=int(sys.argv[3])
        start_id=int(sys.argv[1])
        end_id=int(sys.argv[2])
        create_data_request_for_flares(start_id, end_id)


    
