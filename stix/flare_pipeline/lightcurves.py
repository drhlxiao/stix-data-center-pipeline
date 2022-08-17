#plot STIX lightcurve
#Hualin Xiao (hualin.xiao@fhnw.ch)
# March 01, 2021
import sys
import os
import numpy as np
import math
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from stix.core import mongo_db as db
from stix.core import datatypes as sdt
from stix.spice import time_utils as ut
from stix.analysis import ql_analyzer as qla


mdb = db.MongoDB()
EBINS = [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150, math.inf]
def get_energy_bins(emask):
    ebins=[]
    for i in range(32):
        if emask & (1 << i) != 0:
            ebins.append(i);
    names=[]
    for i in range(len(ebins) - 1):
        begin = ebins[i]
        end = ebins[i + 1]
        if end == 32:
            names.append(f'{EBINS[begin]} keV –⁠ Emax')
        elif end < 32:
            names.append(f'{EBINS[begin]}  – {EBINS[end]} keV')
        else:
            names.append('')
    return names

def get_lightcurve_data(start_utc, end_utc, sort_field='header.unix_time'):

    """Request 

    Args:
        start_unix ([type]): [description]
        duration ([type]): data selection duration
        packet_type ([type]): 
            {'lc':54118,'bkg':54119, 'qlspec':54120, 'var':54121, 'flare':54122}

        sort_field (str, optional): [description]. Defaults to 'header.unix_time'.
    """
    start_unix=ut.utc2unix(start_utc)
    duration=ut.utc2unix(end_utc)-start_unix
    packets=mdb.get_LC_pkt_by_tw(
            start_unix,
            duration)
    return qla.LightCurveAnalyzer.parse(packets)



def plot_lc_and_save(folder, _id, event_name, start_utc, end_utc,  overwrite=True, hvline=None, light_time=0, event_type="Flare #"):
    key='stixlc'
    if  mdb.get_flare_pipeline_products(_id, key) and overwrite == False:
        print(f'{event_type}{event_name} STIX LCs were not created!')
        return 
    fig,ax=plot_lc(start_utc, end_utc,hvline, light_time, event_type, event_name)
    #fig.tight_layout()
    if not fig:
        return None

    filename = os.path.join(folder, f'ql_lc_{_id}_{event_name}.png')
    fig.savefig(filename, dpi=100)
    print(filename)
    mdb.update_flare_pipeline_products(_id, key, [filename])
    return filename

def plot_lc(start_utc, end_utc, fill_between_times=[], light_time=0, event_type='', event_name='', ax=None):
    print("Plotting LC for time range:", start_utc, end_utc)
    data=get_lightcurve_data(start_utc, end_utc)
    if data['num']==0:

        print('No LC data')
        return None, None
    earth_dt = [ut.unix2datetime(x+light_time) for x in  data['time']]
    if ax is None:
        fig, ax = plt.subplots(figsize=(12,6))
    else:
        fig=plt.gcf()

    for i,lc in data['lcs'].items():
        ax.plot(earth_dt, lc, label=data['energy_bins']['names'][i])

    if fill_between_times:
        dts=[ut.utc2datetime(x) for x in fill_between_times]
        if len(dts):
            ax.axvspan(dts[0], dts[1], alpha=0.5, color='cyan')

    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)

    xlabel=f'UTC + {light_time:.02f} s (4 sec time bins)' if light_time!=0 else 'UTC (4 sec time bins)'
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Counts')
    title=f'STIX QL Light Curves'
    ax.set_title(title)
    ax.set_yscale('log')
    ax.legend(loc='upper right', fontsize='x-small', framealpha=0.0)
    return fig, ax
def plot_QL_lc_for_bsd(bsd_id, fill_between_times=[], ax=None):
    bsd_db=mdb.get_collection('bsd')
    doc=bsd_db.find_one({'_id':bsd_id}, {'start_unix':1, 'end_unix':1})
    try:
        plot_lc(ut.unix2utc(doc['start_unix']), ut.unix2utc(doc['end_unix']), fill_between_times, ax=ax)
    except (TypeError, KeyError, ValueError):

        print(f'Failed to plot light curve for {bsd_id}')



if __name__=='__main__':
    plot_lc_and_save('.',0,0, '2021-03-01T08:50:00','2021-03-01T08:52:00')
