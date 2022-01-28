import sys
from pprint import pprint
from datetime import datetime
from stix.spice import datetime as sdt

import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
bsd_req= db['bsd_req_forms']
docs=bsd_req.find({}).sort('_id',1)
for  doc in docs:
    doc['start_unix']=sdt.utc2unix(doc['start_utc'])
    doc['end_unix']=sdt.utc2unix(doc['start_utc'])+float(doc['duration'])
    _id=doc['_id']
    print(_id)
    bsd_req.replace_one({'_id': _id}, doc)



