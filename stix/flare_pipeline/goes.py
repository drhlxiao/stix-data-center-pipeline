#!/usr/bin/python3
''' create goes light curves for STIX flares
   author: Hualin Xiao, 2021-08-20
'''
import sys
import os
import json
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime
from matplotlib import pyplot as plt
from stix.core import mongo_db as db
from stix.spice import stix_datetime as sdt
import matplotlib.ticker as mticker

mdb = db.MongoDB()
db=mdb.get_collection('flares')

energy_map= {
            '0.1-0.8nm': 'GOES low', #1.5 keV - 12 keV
            '0.05-0.4nm': 'GOES high' # 3 - 25
}
def plot_flare_goes(folder,_id, overwrite=False):
    #_id is the flare db entry number
    flare_doc=mdb.find_one({'_id':_id})
    if not flare_doc:
        print(f'Flare {_id} does not exist!')
        return
    flare_id=flare_doc['flare_id']
    key='goes'
    if  mdb.get_flare_pipeline_products(_id, key) and overwrite == False:
        print(f'GOES LC for Flare #{_id} has been created!')
        return 
    fname=os.path.join(folder, f'flare_{_id}_{flare_id}_goes.png')
    #start=flare_doc['start_unix']
    #end=flare_doc['end_unix']
    light_time_diff=flare_doc['ephemeris']['light_time_diff']
    start=flare_doc['start_unix']+light_time_diff-120
    end=flare_doc['end_unix'] +light_time_diff + 120

    start_utc=sdt.unix2utc(start)
    end_utc=sdt.unix2utc(end)

    if plot_goes(start_utc, end_utc):
        mdb.update_flare_pipeline_products(_id, 'goes', [fname])

def plot_goes(start_utc, end_utc, fig_filename=None): #start unix time and end unix time
    flux = {}
    start_times={}
    num=0
    start=sdt.utc2unix(start_utc)
    end=sdt.utc2unix(end_utc)
    print(start,end, start_utc, end_utc)
    try:
        data = mdb.get_goes_fluxes(start, end)
        last_time=0
        for d in data:
            unix = d['unix_time']
            if unix<last_time:
                continue
            num+=1
            last_time=unix
            if d['energy'] not in flux:
                flux[d['energy']] = [[], []]
            flux[d['energy']][1].append(d['flux'])
            flux[d['energy']][0].append(d['unix_time'])
            if d['energy'] not in start_times:
                start_times[d['energy']]=d['time_tag']
    except Exception as e:
        print("ERROR:", e)
        return False
    if num==0:
        print('GOES LC not available')
        return None
    fig, ax = plt.subplots()
    t0_utc=None
    for key,val in flux.items():
        if not t0_utc:
            t0_utc=start_times[key]
        t=[datetime.fromtimestamp(x) for x in val[0]]
        ax.plot(t,val[1],label=energy_map.get(key,'unknown'))
    ax.set_ylabel(r'Watts m$^{-2}$')
    ax.set_yscale('log')
    labels = ['A', 'B', 'C', 'M', 'X']
    ax.set_title('GOES x-ray flux')
    ax.set_ylim(1e-9, 1e-2)
    ax2 = ax.twinx()
    ax2.set_yscale("log")
    ax2.set_ylim(1e-9, 1e-2)
    labels = ['A', 'B', 'C', 'M', 'X']
    centers = np.logspace(-7.5, -3.5, len(labels))
    ax2.yaxis.set_minor_locator(mticker.FixedLocator(centers))
    ax2.set_yticklabels(labels, minor=True)
    ax2.set_yticklabels([])
    ax.yaxis.grid(True, 'major')
    ax.xaxis.grid(False, 'major')
    ax.legend(loc='upper right')
    # beautify the x-labels
    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    fig.tight_layout()
    if fig_filename:
        plt.savefig(fig_filename, dpi=100)
    else:
        plt.show()

    return True
if __name__=='__main__':
    plot_goes('2021-08-17T00:00:00','2021-08-17T10:00:00')
