import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
packet_db= db['data_requests']
cursor=packet_db.find()
for doc in cursor:
    doc['unique_ids']=[doc['unique_id']]
    packet_db.replace_one({'_id':doc['_id']}, doc)
    print('fixing:')
            


    

    
    
