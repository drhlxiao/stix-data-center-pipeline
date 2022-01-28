# Importing sunpy_soar registers the client with sunpy
import tempfile
import os

import sunpy_soar
from sunpy.net import Fido
from sunpy.net.attrs import Instrument, Level, Time
from sunpy_soar.attrs import Identifier
from sunpy.map import Map
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
import astropy.units as u
from dateutil import parser as dtparser
import sunpy.visualization.colormaps
from stix.core import mongo_db as db
mdb = db.MongoDB()
def utc2datetime(utc):
    if not utc.endswith('Z'):
        utc += 'Z'
    return dtparser.parse(utc)
# Create search attributes
def plot_eui(folder,_id, flare_id, start_utc, end_utc, sun_angular_diameter_arcmin, overwrite=False):
    key='eui'
    if  mdb.get_flare_joint_obs(_id, key) and overwrite == False:
        print(f'EUI image for Flare {flare_id} is not created!')
        return 
    instrument = Instrument('EUI')
    time = Time(start_utc, end_utc)
    level = Level(1)
    identifier = Identifier('EUI-FSI174-IMAGE')

    result = Fido.search(instrument & time & level & identifier)
    if not result:
        return None
    num=len(result)
    row=result[int(num/2.)]
    temp_folder=tempfile.gettempdir()
    files= Fido.fetch(row, progress=False, path=temp_folder)
    sdo = Map(files[0])
    cmap = plt.get_cmap('solar orbiterfsi174')
    fig = plt.figure(figsize=(6, 6), dpi=100)
    ax = fig.add_subplot(projection=sdo)
    sdo.plot(clip_interval=[1, 100] * u.percent, cmap=cmap, axes=ax)
    plt.colorbar()
    sdo.draw_limb()
    fname=os.path.join(folder, f'EUI_fsi174_flare_{_id}_{flare_id}.png')
    print(fname)
    plt.savefig(fname, dpi=100)
    mdb.update_flare_joint_obs(_id, key, [fname])
if __name__=='__main__':
    plot_eui('.', 0, 0, '2021-02-01T00:00:00','2021-02-02T00:10:00')
