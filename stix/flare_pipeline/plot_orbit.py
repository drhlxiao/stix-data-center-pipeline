#Hualin Xiao (hualin.xiao@fhnw.ch)
import sys
import os
import matplotlib.patches as patches
from matplotlib import pyplot as plt
from stix.spice import solo
from stix.core import mongo_db as db


from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
mdb = db.MongoDB()

def plot_solo_location(folder,_id, flare_id ,peak_utc,  overwrite=False):
    key='loc'
    data=solo.get_solo_ephemeris(peak_utc, peak_utc)
    if  mdb.get_flare_joint_obs(_id, key) and overwrite == False:
        print(f'Location for Flare {flare_id} was not created !')
        return data
    if not data or 'error' in data:
        print(f'data not available. Location for Flare {flare_id} was not created !')
        return None
    fig, ax = plt.subplots()
    sun= patches.Circle((0.0, 0.0), 0.25, alpha=0.8, fc='yellow')
    plt.text(0,0,'Sun')
    ax.add_patch(sun)

    earth= patches.Circle((-1,0), 0.12, alpha=0.8, fc='green')
    ax.add_patch(earth)
    ax.text(-1,0,'Earth')

    ax.plot(data['x'], data['y'], 'x') 
    ax.text(data['x'][-1],data['y'][-1],'SOLO')
    ax.set_xlabel('X (au)')
    ax.set_ylabel('Y (au)')
    ax.set_xlim(-2,2)
    ax.set_ylim(-2,2)
    ax.set_aspect('equal')
    ax.grid()
    
    plt.title(f'SOLO Location at {peak_utc}')

    fname=os.path.join(folder, f'solo_loc_flare_{_id}_{flare_id}.png')
    plt.savefig(fname)
    mdb.update_flare_joint_obs(_id, key, [fname])
    return data

if __name__=='__main__':
    plot_solo_location('.',301,301, '2021-04-01T00:00:00',True)
