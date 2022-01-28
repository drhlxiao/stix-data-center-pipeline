import csv
import pymongo
from dateutil import parser as dtparser
connect = pymongo.MongoClient()
mdb = connect["stix"]
bsd_db= mdb['bsd']
req_db= mdb['bsd_req_forms']

bsd_docs=bsd_db.find({})
for d in bsd_docs:
    unique_id=d.get('unique_id',None)
    print(unique_id)
    if unique_id is not None:
        query={'unique_ids': int(unique_id), 'hidden':False}
        req_form = req_db.find_one(query)
        if req_form:
            d['request_form'] = req_form
            print(f"updating {d['_id']}")
            bsd_db.update_one({'_id':d['_id']},{'$set': {'request_form': req_form}})
        else:
            print('Not found')

