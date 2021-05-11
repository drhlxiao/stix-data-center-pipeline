
import sys
sys.path.append('.')
import os
import json
import numpy as np
from matplotlib import pyplot as plt
from stix.core import mongo_db as db
from stix.spice import stix_datetime
mdb = db.MongoDB()
energy_map= {
            '0.1-0.8nm': 'GOES low',
            '0.05-0.4nm': 'GOES high'
}

def plot_goes(folder,_id, flare_id, start_utc, end_utc, overwrite=False):
    key='goes'
    if  mdb.get_flare_joint_obs(_id, key) and overwrite == False:
        print(f'GOES LC for Flare {flare_id} was not created!')
        return 

    flux = {}
    start_times={}
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

                last_time=unix

                if d['energy'] not in flux:
                    flux[d['energy']] = [[], []]
                    #time_tags[d['energy']]=[]
                flux[d['energy']][1].append(d['flux'])
                flux[d['energy']][0].append(stix_datetime.utc2unix(d['time_tag']))
                if d['energy'] not in start_times:
                    start_times[d['energy']]=d['time_tag']

    except Exception as e:
        print(e)
        return False
    fig=plt.figure()
    t0_utc=''
    for key,val in flux.items():
        t0_utc=start_times[key]
        t0_unix=stix_datetime.utc2unix(t0_utc)
        t=np.array(val[0])-t0_unix
        plt.plot(t,val[1],label=energy_map.get(key,'unknown'))
    plt.xlabel(f'Start time: {t0_utc}' )
    plt.ylabel(r'Flux (W/m$^2$)')
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.title('GOES x-ray flux')
    fname=os.path.join(folder, f'flare_{_id}_{flare_id}_goes.png')
    print(fname)
    plt.savefig(fname, dpi=100)
    mdb.update_flare_joint_obs(_id, key, [fname])
    return True
if __name__=='__main__':
    plot_goes('.',0,'2021-05-05T00:00:00','2021-05-05T01:10:00')
