import sys
import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
qldb= db['quick_look']
def main():
    cursor=qldb.find({'run_id':{'$gte': 555}})
    for raw in cursor:
        stop=raw['start_unix_time']+  raw['duration']
        print(raw['_id'])
        qldb.update({'_id':raw['_id']}, {'$set':{'stop_unix_time': stop}})

                
    
main()

    
    
