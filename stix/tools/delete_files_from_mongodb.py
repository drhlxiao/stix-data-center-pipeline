import sys
import os
from stix.core import mongo_db 
def delete_one(start_run, end_run=-1):
    runs=[]
    if end_run==-1:
        runs=[start_run]
    else:
        runs=range(start_run,end_run)
    print('runs to delete:{}'.format(str(runs)))
    ret=input('Are you sure to delete them ? Y/N: ')
    if ret=='Y':
        ans=input('Keep raw file and delete the reset? Y/N: ')
        keep_raw=False if ans=='N' else True
        mdb = mongo_db.MongoDB()
        mdb.delete_runs(runs, keep_raw)
        print('deleted')


if __name__=='__main__':
    if len(sys.argv)==1:
        print("usage: ./delete_files_from_mongodb <run_begin> [run_end]")
    elif len(sys.argv)==2:
        delete_one(int(sys.argv[1]))
    elif len(sys.argv)==3:
        delete_one(int(sys.argv[1]), int(sys.argv[2]))
