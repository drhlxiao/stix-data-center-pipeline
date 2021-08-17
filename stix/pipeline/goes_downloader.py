#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @author       : Hualin Xiao
# @date         : May. 11, 2021

""" goes x-ray flux grabber
    workflow:
    1) grab GOES x-ray flux data 
    2) extract start time and stop time information
    3) write the indexing information to mongo db and write the data to disk

    Author: Hualin Xiao
    Date: 2020-06-18
"""
import sys
import os
import json
import time
import requests
import pymongo
from datetime import datetime
from dateutil import parser as dtparser
from stix.core import config

GEOS_DATA_DIRECTORY=config.get_config('pipeline.daemon.goes_lc_path')

def utc2unix(utc):
    if isinstance(utc, str):
        if not utc.endswith('Z'):
            utc += 'Z'
        try:
            return dtparser.parse(utc).timestamp()
        except:
            return 0
    elif isinstance(utc, int) or isinstance(utc, float):
        return utc
    else:
        return 0
def save_geos(data, filename=None):
    
    utc_now= datetime.utcnow()
    if not filename:
        filename=os.path.join(GEOS_DATA_DIRECTORY, utc_now.strftime('%Y%m%d%H%M.json'))
        with open(filename,'w') as f:
            f.write(json.dumps(data))
    
    start_unix=utc2unix(data[0]['time_tag'])
    end_unix=utc2unix(data[-1]['time_tag'])
    print(start_unix,end_unix,filename)
    connect = pymongo.MongoClient()
    db = connect["stix"]
    collection_goes= db['goes']
    base_name= os.path.basename(filename)
    abspath = os.path.abspath(filename)
    path = os.path.dirname(abspath)
    doc={'creation_time':utc_now,
        'filename':base_name,
        'path':path,
        'start_unix':start_unix,
        'stop_unix':end_unix
        }
    print(doc)
    _id=collection_goes.save(doc)
    #print('inserted {}.'.format(_id))
    connect.close()


def main():
    url='http://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json'
    print('downloading goes light curve')
    r = requests.get(url)
    data=r.json()
    print('save goes')
    save_geos(data, None)

def import_data(filename):
    with open(filename) as f:
        data=json.loads(f.read())
        save_geos(data, filename)


if __name__=='__main__':
    if len(sys.argv)==1:
        main()
    else:
        import_data(sys.argv[1])



