import pymongo
import sys
import json
from stix.spice import time_utils as sdt
from pprint import pprint
connect = pymongo.MongoClient(port=9123)
db = connect['stix']['events']

#ltp_id='LTP07_v1'


def create_report(_id, subject, start, end, descr, ltp_id):
    report={
    "_id" : _id,
    "author" : "Hualin",
    "subject" :subject,
    "event_type" : "norminal",
    "start_utc" : start,
    "end_utc" : end,
    "description" : descr,
    "submit" : "",
    "creation_time" : sdt.now(),
    "ltp_id":ltp_id,
    "status" : 0,
    "hidden" : False
    }
    pprint(report)
    db.insert_one(report)

name_map={
        'STIX_BASIC': 'STIX in nominal observation mode',
        'STIX_ANALYSIS':'Data request window',
        }


def import_timeline(filename, ltp_id ):
    f=open(filename)
    data = json.load(f)
    observations=data['observations']
    try:
        next_id = db.find({}).sort(
            '_id', -1).limit(1)[0]['_id'] + 1
    except IndexError:
        next_id=0

    for obs in observations:
        name=obs['name']
        start=obs['startDate']
        end=obs['endDate']
        if name=='STIX_BASIC':
            continue
        vol=obs['numberParameters']['SCI_REQUEST_VOLUME']
        data_vol=f"Max science data volume: {vol['value']} {vol['unit']}"
        subject=name_map.get(name, name)

        create_report(next_id, subject, start, end, descr=data_vol, ltp_id=ltp_id)
        next_id+=1



if __name__=='__main__':
    if len(sys.argv)==3:
        filename=sys.argv[1]
        ltp_id=sys.argv[2]
        import_timeline(filename, ltp_id)
    else:
        print('./import  <filename> <LTP_ID>')
