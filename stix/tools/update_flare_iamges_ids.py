import sys
import os
import json
import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
col_flare=db['flare_images']
flare_id=0
for doc in col_flare.find({}).sort('_id',-1):
    _id=doc['_id']
    print('_id',id)
    doc['_id']=flare_id
    col_flare.save(doc)
    col_flare.delete_one({'_id':_id})
    flare_id+=1
