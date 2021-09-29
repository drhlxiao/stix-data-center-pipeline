import json
import sys
sys.path.append('/home/xiaohl/web')
from core import mongodb_api
from core import user_request_export as ure
mdb = mongodb_api.MongoDB()
data_request_database = mdb.get_collection('bsd_req_forms')
MiB = 1024 * 1024

STP = 169

MAX_NUM_TCs = 340
SIZE_LIMITS = [10 * MiB] * 7
NUM_IORS = 7
start_id = 2200
end_id = 3500
#data request entry id, can be found on the web page

aspect_requests_ids = []
rows = []
selected_ids = []


def write_json(file_id, data):
    fname = f'STP_{STP}_{file_id}.json'
    print("Writing ", fname)
    f = open(fname, 'w')
    j = json.dumps(data, indent=4)
    f.write(j)
    f.close()


def add_a_request(form):
    if form['_id'] not in selected_ids:
        rows.append(form)
        selected_ids.append(form['_id'])


forms = data_request_database.find({
    '_id': {
        '$gte': start_id,
        '$lte': end_id
    },
    'hidden': False
}).sort('_id', 1)

for form in forms:
    uids = form.get('unique_ids', [])
    if uids:
        ior_id = mdb.get_ior_ids_by_request_uids(uids)
        tms = mdb.get_data_request_tm_entries(uids)
        print(form['_id'], ior_id, tms, form['request_type'])
        if ior_id == -1 and tms == 0 and form['request_type'] != 'Aspect':
            add_a_request(form)
    else:
        #uid has not been assigned
        if form['request_type'] != 'Aspect':

            add_a_request(form)

print('Number of requests:', len(rows))

for i in range(NUM_IORS):
    total_volume = 0
    total_selected = 0
    ids = []

    for row in rows:
        size = float(row['data_volume_upper_limit'])
        if len(ids) < MAX_NUM_TCs and total_volume < SIZE_LIMITS[i]:
            selected = True
            #here defines the strategy of selections
            if size > 0.3 * SIZE_LIMITS[
                    i] and total_volume > 0.25 * SIZE_LIMITS[i]:
                selected = False
            if size > 0.1 * SIZE_LIMITS[
                    i] and total_volume > 0.5 * SIZE_LIMITS[i]:
                selected = False
            if size > 0.05 * SIZE_LIMITS[
                    i] and total_volume > 0.7 * SIZE_LIMITS[i]:
                selected = False

            if selected:
                total_volume += size
                ids.append(row['_id'])
                total_selected += 1
                rows.remove(row)

    if i == 0:
        ids.extend(aspect_requests_ids)
    print(f"TCs:{len(ids)}, Total volume:{total_volume/MiB:.3f}")
    print(f'Left: {len(rows)}')
    #data = ure.create_occurrences(data_request_database, ids)
    write_json(i, data)

print("total number of selected requests:", total_selected)
