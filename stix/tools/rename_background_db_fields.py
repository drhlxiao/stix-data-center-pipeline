import sys
import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
db= self.db['qllc_statistics']
def main():
    cursor=db.find({'run_id':{'$gte': 555}})
    for raw in cursor:
        stop=raw['start_unix_time']+  raw['duration']
        print(raw['_id'])
        db.update({'_id':raw['_id']}, {'$set':{'stop_unix_time': stop}})

                
    
main()

    
    
