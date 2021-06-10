#Author Hualin Xiao (hualin.xiao@fhnw.ch)
#plot AIA image for STIX flares
import sys
sys.path.append('.')
import os
import io
import astropy.units as u
from matplotlib import pyplot as plt
from sunpy.net import Fido, attrs as a
from sunpy.map import Map
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from stix.spice import stix_datetime
import tempfile

from stix.core import mongo_db as db
mdb = db.MongoDB()

def plot_aia(folder,_id, flare_id ,utc_start,  wavelen=131, overwrite=False):
    key='aia131'
    if  mdb.get_flare_joint_obs(_id, key) and overwrite == False:
        print(f'GOES LC for Flare {flare_id} is not created!')
        return 
    unix_start= stix_datetime.utc2unix(utc_start)
    utc_end = stix_datetime.unix2datetime(unix_start +
                                          60).strftime('%Y-%m-%dT%H:%M:%S')
    sdo_query = Fido.search(a.Time(utc_start, utc_end), a.Instrument('AIA'),
                            a.Wavelength(wavelen* u.angstrom))
    temp_folder=tempfile.gettempdir()
    sdo_res = Fido.fetch(sdo_query[0], progress=False, path=temp_folder)
    if not sdo_res:
        print('data not available')
        return ''
    sdo = Map(sdo_res[0])
    fig = plt.figure(figsize=(6, 6), dpi=100, facecolor='white')
    ax = fig.add_subplot(projection=sdo)
    sdo.plot(clip_interval=[1, 100] * u.percent, axes=ax)
    plt.colorbar()
    sdo.draw_limb()
    fname=os.path.join(folder, f'AIA_{wavelen}_{_id}_{flare_id}.png')
    plt.savefig(fname, dpi=100)
    mdb.update_flare_joint_obs(_id, key, [fname])

    return fname

if __name__=='__main__':
    plot_aia('.',0,0,'2021-04-01T00:00:00', True)
