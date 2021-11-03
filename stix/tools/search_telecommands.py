import sys
from pprint import pprint
from stix.spice import stix_datetime

import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
col_iors= db['iors']

names=['ZIX36005','ZIX36004','AIXF060A','AIXF061A']
start_utc='2020-04-15T00:00:00'
end_utc='2021-04-15T00:00:00'

start_unix=stix_datetime.utc2unix(start_utc)
end_unix=stix_datetime.utc2unix(end_utc)
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
        results.append(tc)

pprint(results)

