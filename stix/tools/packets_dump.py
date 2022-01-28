import sys
import os
import pprint
from stix.core import mongo_db
mdb= mongo_db.MongoDB()

db=mdb.get_db()['packets']

pkts=list(db.find({'run_id':1355}).sort('_id',1))
import pickle
with open('tm_1355.pkl', 'wb') as f:
    pickle.dump(pkts, f)


