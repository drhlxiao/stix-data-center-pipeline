import sys
from pprint import pprint
from stix.spice import datetime

import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
col_iors= db['iors']

#names=['ZIX36005','ZIX36004','AIXF060A','AIXF061A']
names=['AIXF414A'] #load prameters
start_utc='2021-12-09:00:00'
end_utc='2022-11-15T00:00:00'

start_unix=datetime.utc2unix(start_utc)
end_unix=datetime.utc2unix(end_utc)
query_string={'startUnix': { '$gte': start_unix,    '$lt':end_unix },
            'status':{'$gt':0},
            'occurrences':{'$elemMatch':{ 'name': {'$in':names}}},
        }
print(query_string)
iors=col_iors.find(query_string).sort('_id',-1)

results=[]
for  ior in iors:
    occurrences=ior['occurrences']
    for tc in occurrences:
        if tc['name'] in names:
            results.append(tc)

pprint(results)

