# Importing sunpy_soar registers the client with sunpy
import tempfile
import os

import sunpy_soar
from sunpy.net import Fido
from sunpy.net.attrs import Instrument, Level, Time
from sunpy_soar.attrs import Product
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
    instrument = Instrument('EUI')
    time = Time(start_utc, end_utc)
    level = Level(2)
    identifier = Product('EUI-FSI174-IMAGE')
    result = Fido.search(instrument & time & level & identifier)
    if not result:
        return None
    num=len(result)
    row=result[int(num/2.)]
    temp_folder=tempfile.gettempdir()
    files= Fido.fetch(row, progress=False, path=temp_folder)
    eui_map = Map(files[0])
    cmap = plt.get_cmap('solar orbiterfsi174')
    fig = plt.figure(figsize=(6, 6), dpi=100)
    ax = fig.add_subplot(projection=eui_map)
    eui_map.plot(clip_interval=[1, 100] * u.percent, cmap=cmap, axes=ax)
    plt.colorbar()
    eui_map.draw_limb()
    fname=os.path.join(folder, f'EUI_fsi174_flare_{_id}_{flare_id}.png')
    print(fname)
    plt.savefig(fname, dpi=100)
    plt.show()
if __name__=='__main__':
    plot_eui('.', 0, 0, '2022-03-01T00:00:00','2022-03-02T00:10:00',60)
