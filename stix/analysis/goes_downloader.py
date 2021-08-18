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
from stix.spice import stix_datetime as sdt 
from stix.core import mongo_db as db
import math

mdb = db.MongoDB()

class GOES(object):
    def __init__(self):
        self.db= mdb.get_collection('goes_fluxes')
        try:
            self.last_unix_time = self.db.find().sort('unix_time',-1).limit(1)[0]['unix_time']
        except:
            self.last_unix_time=0
        print('Last unix time ', self.last_unix_time)
    def get_max_time(self):
        return self.last_unix_time
    def save_geos_fluxes(self, data):
        num = 0
        for d in data:
            d['unix_time']=sdt.utc2unix(d['time_tag'])
            if d['unix_time']<self.last_unix_time:
                #don't insert again
                continue
            ret=self.db.update_one(
                    {'unix_time':d['unix_time'], 'energy':d['energy']},
                    {'$set': d},
                    upsert=True)
            num+=1
            self.last_unix_time=d['unix_time']
        print(f'{num} entries inserted ')
    def download(self, max_unix=math.inf):
        if max_unix<self.last_unix_time:
            return
            

        url='http://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json'
        print('downloading GOES light curve')
        r = requests.get(url)
        data=r.json()
        self.save_geos_fluxes(data)
    def import_data(self, filename):
        with open(filename) as f:
            data=json.loads(f.read())
            self.save_geos_fluxes(data)
    def loop(self, wait=48*3600):
        while True:
            try:
                self.download()
            except Exception as e:
                print(e)
            print('Waiting ...')
            time.sleep(wait)





if __name__=='__main__':
    p=GOES()
    if len(sys.argv)==1:
        p.download()
    else:
        p.import_data(sys.argv[1])



