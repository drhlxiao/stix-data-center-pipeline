# find requests that executed multiple times
# hualin.xiao created on Feb. 11, 2022
import sys
import os
import pymongo

connect = pymongo.MongoClient()
db = connect["stix"]
col_bsd=db['bsd']
col_pkt=db['packets']

bsd_query={}
f=open('bsd_multiple.log','w')


for doc in col_bsd.find(bsd_query).sort('_id',-1):
    request_id=doc.get('unique_id',None)
    if request_id is None:
        continue
    packet_ids=doc['packet_ids']
    print(doc["_id"])
    
    segs=[0]*4
    for pkt in col_pkt.find({'_id':{'$in':packet_ids}}).sort('_id',1):
        seg=pkt['header']['seg_flag']
        segs[seg]+=1
        this_utc=pkt['header']['UTC']

        if segs[1] >1 or segs[2]>1:
            msg=f">>Multiple requests:{doc['_id']}, {doc['unique_id']}, {this_utc}/{last_utc}" 
            print(msg)
            f.write(msg)
            f.flush()
            break
        last_utc=this_utc
f.close()

