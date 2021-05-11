#Hualin Xiao (hualin.xiao@fhnw.ch)
import sys
sys.path.append('.')
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
    if  mdb.get_flare_joint_obs(_id, key) and overwrite == False:
        print(f'Location for Flare {flare_id} was not created !')
        return 
    data=solo.get_solo_ephemeris(peak_utc, peak_utc)
    if not data or 'error' in data:
        print(f'data not available. Location for Flare {flare_id} was not created !')
        return
    fig, ax = plt.subplots(subplot_kw=dict(aspect="equal"))
    sun= patches.Circle((0.0, 0.0), 0.15, alpha=0.8, fc='yellow')
    plt.text(0,0,'Sun')
    ax.add_patch(sun)

    earth= patches.Circle((0.0, -1), 0.05, alpha=0.8, fc='blue')
    ax.add_patch(earth)
    ax.text(0,-1,'Earth')

    ax.plot(data['x'], data['y'], 'x') 
    ax.text(data['x'][-1],data['y'][-1],'SOLO')
    ax.set_xlabel('X (au)')
    ax.set_ylabel('Y (au)')
    plt.title(f'SOLO Location at {peak_utc}')

    fname=os.path.join(folder, f'solo_loc_flare_{_id}_{flare_id}.png')
    plt.savefig(fname)
    mdb.update_flare_joint_obs(_id, key, [filename])

if __name__=='__main__':
    plot_solo_location('.',0,0, '2021-04-01T00:00:00')
