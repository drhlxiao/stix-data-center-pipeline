import sys
import os
import pprint
from stix.core import mongo_db
print('connecting server')
mdb= mongo_db.MongoDB()

db=mdb.get_db()['packets']

ids=[1244, 1355]

for _id in ids:
    pkts=list(db.find({'run_id':_id}).sort('_id',1))
    import pickle
    with open(f'tm_{_id}.pkl', 'wb') as f:
        pickle.dump(pkts, f)


