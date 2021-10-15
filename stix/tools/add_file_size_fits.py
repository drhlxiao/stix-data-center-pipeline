import sys
import pymongo
import os
connect = pymongo.MongoClient()
db = connect["stix"]
fits_db= db['fits']

cursor=fits_db.find({})
for raw in cursor:
    fn=os.path.join(raw['path'],raw['filename'])
    sz=os.path.getsize(fn)
    #sz=100
    fits_db.update_one({'_id':raw['_id']}, {'$set':{'file_size':sz}})


            






