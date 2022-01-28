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
import ROOT
import numpy as np
from datetime import datetime, timedelta
from dateutil import parser as dtparser
from stix.core import config
from stix.spice import datetime as sdt 
from stix.core import mongo_db as db
from matplotlib import pyplot as plt
import math

mdb = db.MongoDB()




class GOES(object):
    def __init__(self):
        self.db= mdb.get_collection('goes_fluxes')
        try:
            self.last_unix_time = self.db.find().sort('unix_time',-1).limit(1)[0]['unix_time']
        except:
            self.last_unix_time=0
        print('Last UTC ', sdt.unix2utc(self.last_unix_time))

    def estimate_goes_background(self, start_unix_time,
            duration=7*86400, 
            w=900, 
            create_ql_plot=False):
        energy="0.1-0.8nm"
        rows = self.db.find({'unix_time':{'$gt': start_unix_time, '$lt':start_unix_time+duration},
            'energy':energy
            }).sort('unix_time',1)
        
        counts=[]
        keys=[]
        last_time=0
        for row in rows:
            if row['unix_time'] <=last_time:
                continue
            counts.append(row['flux'] )
            keys.append({'unix_time':row['unix_time'], 'energy':energy})
            last_time=row['unix_time']

        nbins = len(counts)
        if nbins==0:
            return
        source = np.copy(counts)
        s= ROOT.TSpectrum()
        for i,c in enumerate(counts):
            source[i]=c
        s.Background(source,nbins,w,1,2,0,3,0)
        for   key, bkg in zip(keys, source):
            
            if key['unix_time'] < start_unix_time+0.2*duration :
                #to avoid bad estimates  at edges
                continue
            if key['unix_time'] > start_unix_time+0.2*duration:
                break

            self.db.update_many(
                    key,
                    {'$set': {'background':bkg}
                        },
                    upsert=False)
        if create_ql_plot:
            fig=plt.figure()
            plt.plot(counts)
            plt.plot(source)

            utc=sdt.unix2utc(start_unix_time)
            plt.xlabel('Time (min) since {utc}')
            plt.title('GOES-flux T0: '+utc)
            plt.yscale('log')
            goes_lc_path=config.get_config('pipeline.daemon.goes_lc_path')
            fname=os.path.join(goes_lc_path, f'goes_flux_with_baseline_{utc}.png')
            print(fname)
            plt.savefig(fname)
                
            

    

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
        print(f'Last timestamp: {sdt.unix2utc(self.last_unix_time)}')

        duration=7*86400
        try:
            self.estimate_goes_background(self.last_unix_time - duration, duration, w=900 )
        except Exception as e:
            print(e)

    def download(self, max_unix=math.inf):
        if max_unix<self.last_unix_time:
            return
        url='http://services.swpc.noaa.gov/json/goes/primary/xrays-3-day.json'
        print('Downloading GOES x-ray flux')
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
        for w in range(200):
            start_dt=datetime.now()-timedelta(days=w*2.9)
            print(start_dt)
            start_unix=start_dt.timestamp()
            p.estimate_goes_background(start_unix,  w=900, create_ql_plot= True)


