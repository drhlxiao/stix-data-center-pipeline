import sys
import os
from stix.core import mongo_db 
from stix.spice import datetime 
from stix.pipeline import parser_pipeline as pp
mdb = mongo_db.MongoDB()
MAX_ALLOWED_DIFF=0.5

def test(start_run, end_run=-1):
    runs=[]
    if end_run==-1:
        runs=[start_run]
    else:
        runs=range(start_run,end_run)
    raw_db=mdb.get_collection('raw_files')
    filenames=[]
    docs=raw_db.find({'_id':{'$gte':start_run, '$lte':end_run},'hidden':False}).sort('_id',1)
    for doc in docs:
        start_unix=doc['data_start_unix_time']
        end_unix=doc['data_stop_unix_time']
        start_scet=doc['data_start_scet']
        end_scet=doc['data_end_scet']
        if start_scet==0 or end_scet ==0 :
            continue
        new_start_unix=datetime.scet2unix(start_scet)
        new_end_unix=datetime.scet2unix(end_scet)
        delta_start=abs(start_unix-new_start_unix)
        delta_end=abs(end_unix-new_end_unix)
        if delta_start>MAX_ALLOWED_DIFF or delta_end>MAX_ALLOWED_DIFF:
            print(doc['_id'],' DTstart:', delta_start, ', DTend (s): ' , delta_end)




if __name__=='__main__':
    if len(sys.argv)==1:
        print("usage: ./test <run_begin> [run_end]")
    elif len(sys.argv)==2:
        test(int(sys.argv[1]))
    elif len(sys.argv)==3:
        test(int(sys.argv[1]), int(sys.argv[2]))
