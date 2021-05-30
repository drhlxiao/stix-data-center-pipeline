
import sys
sys.path.append('.')
import os
import json
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime
from matplotlib import pyplot as plt
from stix.core import mongo_db as db
from stix.spice import stix_datetime
import matplotlib.ticker as mticker
mdb = db.MongoDB()
energy_map= {
            '0.1-0.8nm': 'GOES low', #1.5 keV - 12 keV
            '0.05-0.4nm': 'GOES high' # 3 - 25
}
def get_class(x):
    if x<1e-7:
        return 'A'
    elif x<1e-6:
        return f'B{x/1e-7:.1f}'
    elif x<1e-5:
        return f'C{x/1e-6:.1f}'
    elif x<1e-4:
        return f'M{x/1e-5:.1f}'
    else:
        return f'X{x/1e-4:.1ff}'




def plot_goes(folder,_id, flare_id, start_utc, end_utc, peak_utc=None, overwrite=False):
    key='goes'
    if  mdb.get_flare_joint_obs(_id, key) and overwrite == False:
        print(f'GOES LC for Flare {flare_id} was not created!')
        return 

    flux = {}
    start_times={}
    num=0
    try:
        start = stix_datetime.utc2unix(start_utc)
        end = stix_datetime.utc2unix(end_utc)
        files= list(mdb.get_goes_x_ray_flux_file_list(start, end))
        last_time=0
        for item in files:
            filename = os.path.join(item['path'], item['filename'])
            json_file = open(filename, 'r')
            data = json.load(json_file)
            for d in data:
                unix = stix_datetime.utc2unix(d['time_tag'])
                if unix < start:
                    continue
                if unix > end:
                    break
                if unix<last_time:
                    continue
                num+=1

                last_time=unix

                if d['energy'] not in flux:
                    flux[d['energy']] = [[], []]
                    #time_tags[d['energy']]=[]
                flux[d['energy']][1].append(d['flux'])
                flux[d['energy']][0].append(stix_datetime.utc2unix(d['time_tag']))
                if d['energy'] not in start_times:
                    start_times[d['energy']]=d['time_tag']

    except Exception as e:
        print("ERROR:", e)
        return False
    if num==0:
        return
    fig, ax = plt.subplots()
    max_flux=0
    t0_utc=None
    for key,val in flux.items():
        if not t0_utc:
            t0_utc=start_times[key]
        #t0_unix=stix_datetime.utc2unix(t0_utc)
        #t=np.array(val[0])-t0_unix
        t=[datetime.fromtimestamp(x) for x in val[0]]
        if key=='0.1-0.8nm':
            max_flux=np.max(val[1])
            imax=np.array(val[1]).argmax()
            goes_peak_unix_time=val[0][imax]
        ax.plot(t,val[1],label=energy_map.get(key,'unknown'))

    goes_peak_utc=stix_datetime.unix2utc(goes_peak_unix_time)
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

    fname=os.path.join(folder, f'flare_{_id}_{flare_id}_goes.png')
    print(fname)
    fig.tight_layout()
    plt.savefig(fname, dpi=100)
    mdb.update_flare_joint_obs(_id, 'goes', [fname])
    mdb.update_flare_field(_id, 'goes', {'flux': max_flux, 'class':get_class(max_flux), 'goes_peak_unix_time': goes_peak_unix_time, 'goes_peak_utc':goes_peak_utc})
    print(_id, {'flux': max_flux, 'class':get_class(max_flux)})
    return True
if __name__=='__main__':
    plot_goes('.',0,2,'2021-05-05T00:00:00','2021-05-05T01:10:00', '2021-05-05T02:05:00', True)
