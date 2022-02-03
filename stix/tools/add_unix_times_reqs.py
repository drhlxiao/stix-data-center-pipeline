import pymongo
from stix.spice import datetime
connect = pymongo.MongoClient()
db = connect["stix"]
packet_db= db['data_requests']
cursor=packet_db.find()
for doc in cursor:
    doc['start_unix']=datetime.utc2unix(doc['start_utc'])
    doc['end_unix']=doc['start_unix']+float(doc['duration'])
    doc['end_utc']=datetime.unix2utc(doc['end_unix'])

    packet_db.replace_one({'_id':doc['_id']}, doc)
    print('fixing:',doc['_id'])
            


    

    
    
