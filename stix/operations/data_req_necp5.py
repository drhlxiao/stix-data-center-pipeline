#!/usr/bin/python
"""
   Export data requests from stix data center
   Author: Hualin Xiao
   Date:   Sept. 21, 2021

"""

import json
import requests
HOST='https://pub023.cs.technik.fhnw.ch/'
URLS={'export':f'{HOST}/api/request/user-data-requests/ior/export',
    'query':f'{HOST}/api/request/user-data-requests/pending',
    }

MiB = 1024 * 1024
token='9b3fcb0a60153715665aac305c9954e5'

# don't modify above lines 


STP = 169
MAX_NUM_TCs = 340
#max number of telecommands in each IOR

NUM_IORS = 7  
#number  of IORs to be created

SIZE_LIMITS = [10 * MiB] * NUM_IORS
# data volume limits 


start_id = 3490  
end_id = 3520
# data request entry ids can be found on the web page
# processed or executed  data requests will be ignored 

aspect_requests_ids = []
#Put the IDs of aspect data requests here





def fetch_data_request_from_server(start_id, end_id):
    url=URLS['query']
    form={'start_id':start_id, 'end_id':end_id, 'token':token}
    response = requests.post(url, data=form)
    data = response.json()
    return data

def download_requests(file_number, selected_ids):
    """
        Exporting  data requests to json file
        
    """
    if not selected_ids:
        print(f"Warning: the number of data requests is zero for IOR {file_number}")
        return
    url=URLS['export']
    ids=','.join([str(i) for i in selected_ids])
    form={'ids':ids,'token':token}
    fname = f'STP_{STP}_{file_number}.json'
    response = requests.post(url, data=form)
    data = response.json()
    print("Writing ", fname)
    f = open(fname, 'w')
    j = json.dumps(data, indent=4)
    f.write(j)
    f.close()



def main():
    req_forms=fetch_data_request_from_server(start_id, end_id)
    if 'error' in req_forms:
        print(req_forms)
        return
    data_requests = [ form  for form in req_forms  if form['request_type'] != 'Aspect']
    #exclude aspect requests, one needs to include aspect data request manually 

    print('Number of requests:', len(data_requests))
    for i in range(NUM_IORS):
        total_volume = 0
        total_selected = 0
        ids = []

        for req in data_requests:
            size = float(req['worst_case_data_volume'])
            if len(ids) < MAX_NUM_TCs and total_volume < SIZE_LIMITS[i]:
                selected = True
                #here defines the strategy of  selections
                #trying to include  as many small size requests as possible

                #don't include big volume requests when almost reaching the limit
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
                    ids.append(req['_id'])
                    total_selected += 1
                    data_requests.remove(req)

        if i == 0 and aspect_requests_ids:
            #manually add aspect data requests in the first IOR
            ids.extend(aspect_requests_ids)
        print(f"TCs:{len(ids)}, Total volume:{total_volume/MiB:.3f}")
        print(f'Nb. of remaining selected data requests: {len(data_requests)}')
        download_requests(i, ids)

    print("total number of selected requests:", total_selected)

if __name__=='__main__':
    main()
