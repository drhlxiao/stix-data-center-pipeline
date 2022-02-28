
import sys
sys.path.append('.')
import pymongo
connect = pymongo.MongoClient('localhost', 9000)
stix = connect['stix']
collection = stix['iors']
dreq_collection = stix['data_requests']


iors = collection.find().sort('_id', -1)
print('start')
for ior in iors:
    print(ior['_id'])
    if 'occurrences' not in ior:
        continue
    if len(ior['occurrences']) == 0:
        continue
    for occ in ior['occurrences']:
        if occ['name']=='AIXF414A' :
            if int(occ['parameters'][0][1]) == 477:
                print(ior['_id'], occ)


