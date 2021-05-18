import sys
sys.path.append('.')
from pprint import pprint
from datetime import datetime
from stix.spice import stix_datetime as sdt

import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
bsd_req= db['bsd_req_forms']



def update_request(doc, new_id):
    T=float(doc['duration'])
    emax=float(doc['emax'])
    emin=float(doc['emin'])
    eunit=float(doc['eunit'])
    E = (emax - emin + 1)/(eunit+1)
    data_volume = 1.1 * T * (E + 4)
    doc.update({
    "_id":new_id,
    "data_volume" : str(int(data_volume)),
    "request_type" : 'Spectrogram',
    "author" : 'Hualin Xiao',
    "email" : 'hualin.xiao@fhnw.ch',
    "time_bin" : "1",
    "detector_mask" : "0xFFFFFCFF",
    "description" :f"{doc['description']}  L4 request for L1 Req # {doc['_id']}",
    "volume" : str(int(data_volume)),
    "creation_time":datetime.now(),
    "unique_ids" : []
        })
    print(doc)
    #bsd_req.save(doc)
    return  data_volume


start_utc='2021-01-01T00:00:00'
end_utc='2021-06-15T00:00:00'
start_unix=sdt.utc2unix(start_utc)
end_unix=sdt.utc2unix(end_utc)
max_id=[x for x in bsd_req.find().sort('_id',-1).limit(1)][0]['_id']
print('max id:', max_id)

query_string={'request_type':'L1','hidden':False, '_id':{'$lte':828}}
docs=bsd_req.find(query_string).sort('_id',1)
new_id=max_id+1
total_volume=0
num_req=0
for  doc in docs:
    try:
        ut=sdt.utc2unix(doc['start_utc'])
    except KeyError:
        print(doc)
        continue
    if  ut>start_unix and ut<end_unix:
        total_volume+=update_request(doc, new_id)
        new_id+=1
        num_req+=1


print('Total requests:',num_req)
print('Total volume:',total_volume/(1024.*1024))
