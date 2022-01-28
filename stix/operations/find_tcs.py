import pymongo
connect = pymongo.MongoClient('localhost', 9000)
stix = connect['stix']
collection = stix['iors']
dreq_collection = stix['bsd_req_forms']


def search_occurrences(search_strings):
    results = []
    iors = collection.find().sort('_id', -1)
    print(iors.count())
    for ior in iors:
        if 'occurrences' not in ior:
            continue
        if len(ior['occurrences']) == 0:
            continue
        for occ in ior['occurrences']:
            for _str in search_strings:
                if str(_str) in str(occ):
                    results.append(occ)
                    break
    return results


# In[21]:

requests_ids = set([
    1169781776, 1169843728, 1169870352, 1169895952, 1169922576, 1169949200,
    1169974800, 1170001424, 1170027024, 1170053648, 1170211344, 1170236944,
    1170263568, 1170289168, 1170315792, 1170327568, 1170342416, 1170367760
])

occurrences = search_occurrences(requests_ids)

# In[22]:

len(occurrences)

# In[25]:

print(occurrences[-1])

# In[50]:

import sys
from stix.spice import datetime
import pymongo
connect = pymongo.MongoClient('localhost', 9000)
stix = connect['stix']
collection = stix['operation_requests']
dreq_collection = stix['bsd_req_forms']
inserted = []


def create_request(_id, occ):
    param = occ['parameters']
    uid = int(param[0][1])
    if uid in inserted:
        return
    inserted.append(uid)

    print(int(param[2][1]))
    print(str(datetime.scet2utc(int(param[2][1]))))
    doc = {
        "_id": _id,
        "author": "Robot",
        "subject": f"Crab Re-request {uid} ",
        "purpose": "4-5 keV crab observation data",
        "execution_date": "IX5 Day 2",
        "unique_id": uid + 1,
        "request_type": "L" + param[1][1],
        "start_utc": str(datetime.scet2utc(int(param[2][1]))),
        "duration": str(int(int(param[7][1]) / 10)),
        "averaging": "",
        "time_bin": str(int(int(param[8][1]) / 10)),
        "detector_mask": param[4][1],
        "pixel_mask": "0xFFF",
        "emin": "1",
        "emax": "2",
        "eunit": "1",
        "description": "crab data request, energy range 4-5 keV",
        "submit": "",
        "volume": "643738",
        "creation_time":
        datetime.utc2datetime("2020-10-23T19:05:46.280Z"),
        "status": '0',
        "hidden": False,
        "data_volume": "643738"
    }
    print(f'writting {_id}\n')
    dreq_collection.save(doc)


# In[51]:

_id = 214
for occ in occurrences:
    create_request(_id, occ)
    _id += 1

# In[ ]:
