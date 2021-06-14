import pymongo
connect = pymongo.MongoClient()
#('localhost',9000)
stix = connect['stix']
collection = stix['iors']
dreq_collection = stix['bsd_req_forms']


def search_occurrences(start_parameter_id, end_parameter_id):
    results = {}
    iors = collection.find().sort('_id', -1)
    print(iors.count())
    for ior in iors:
        if 'occurrences' not in ior:
            continue
        if len(ior['occurrences']) == 0:
            continue
        for occ in ior['occurrences']:
            if occ['name'] == 'AIXF414A':
                parameter_id = int(occ['parameters'][0][1])
                if parameter_id >= start_parameter_id and parameter_id <= end_parameter_id:
                    if parameter_id not in results:
                        results[parameter_id] = []

                    occ['_id'] = ior['_id']
                    results[parameter_id].append(occ)
    return results


# In[21]:
parameter_map = {
    '273': 'RCR_ENABLED',
    '274': 'RCR_L0',
    '275': 'RCR_L1',
    '276': 'RCR_L2',
    '277': 'RCR_L3',
    '278': 'RCR_B0',
    '279': 'RCR_GROUP_MASK',
    '280': 'RCR_PIXEL_0',
    '281': 'RCR_PIXEL_1',
    '282': 'RCR_PIXEL_2N',
    '283': 'RCR_PIXEL_2S',
    '284': 'RCR_PIXEL_3N1',
    '285': 'RCR_PIXEL_3N2',
    '286': 'RCR_PIXEL_3S1',
    '287': 'RCR_PIXEL_3S2',
    '288': 'RCR_PIXEL_4N1',
    '289': 'RCR_PIXEL_4N2',
    '290': 'RCR_PIXEL_4N3',
    '291': 'RCR_PIXEL_4N4',
    '292': 'RCR_PIXEL_4S1',
    '293': 'RCR_PIXEL_4S2',
    '294': 'RCR_PIXEL_4S3',
    '295': 'RCR_PIXEL_4S4',
    '296': 'RCR_PIXEL_5',
    '297': 'RCR_PIXEL_6_1',
    '298': 'RCR_PIXEL_6_2',
    '299': 'RCR_PIXEL_7_1',
    '300': 'RCR_PIXEL_7_2',
    '301': 'RCR_PIXEL_7_3',
    '302': 'RCR_PIXEL_7_4',
    '303': 'RCR_PIXEL_B0',
    '304': 'RCR_PIXEL_B1',
    '305': 'RCR_PIXEL_B2',
    '306': 'RCR_PIXEL_B3',
    '307': 'RCR_PIXEL_B4',
    '308': 'RCR_PIXEL_B5',
    '309': 'RCR_PIXEL_B6',
    '310': 'RCR_PIXEL_B7',
    '311': 'REGIM',
}

occurrences = search_occurrences(274, 311)

len(occurrences)

for i, occs in occurrences.items():
    print('-------------------------')
    print(i)
    print(parameter_map[f'{i}'])
    #print(occs)
    values = [x['parameters'][2][1] for x in occs]
    ids = [x['_id'] for x in occs]
    print(values)
    print(ids)
