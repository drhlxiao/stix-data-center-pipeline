#reduce lc packet data volume before sending to browsers
#Hualin Xiao (hualin.xiao@fhnw.ch)
# March 01, 2021
import sys
sys.path.append('.')
import os
import numpy as np
import math
from matplotlib import pyplot as plt
from pprint import pprint
from stix.core import mongo_db as db
from stix.core import stix_datatypes as sdt
from stix.spice import stix_datetime
from stix.analysis import ql_analyzer as qla


mdb = db.MongoDB()
EBINS = [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150, 'Emax']
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

class LCMinifier(object):
    def __init__(self):
        pass

    @staticmethod
    def request_by_tw(start_utc, end_utc, sort_field='header.unix_time'):

        """Request 

        Args:
            start_unix ([type]): [description]
            duration ([type]): data selection duration
            packet_type ([type]): 
                {'lc':54118,'bkg':54119, 'qlspec':54120, 'var':54121, 'flare':54122}

            sort_field (str, optional): [description]. Defaults to 'header.unix_time'.
        """
        start_unix=stix_datetime.utc2unix(start_utc)
        duration=stix_datetime.utc2unix(end_utc)-start_unix
        packets=mdb.get_LC_pkt_by_tw(
                start_unix,
                duration)
        return qla.LightCurveAnalyzer.parse(packets)



def plot_stix_lc(folder, _id, flare_id, start_utc, end_utc,  overwrite=True, T0_utc=None, light_time=0):
    key='stixlc'
    if  mdb.get_flare_joint_obs(_id, key) and overwrite == False:
        print(f'Flare {flare_id} STIX LCs were not created!')
        return 
    data=LCMinifier.request_by_tw(start_utc, end_utc)
    if data['num']==0:
        print('No LC data')
        return
    if not T0_utc:
        t_since_t0 = data['time']- data['time'][0]
        start_utc=stix_datetime.unix2utc(data['time'][0])
    else:
        start_utc=T0_utc
        t0=stix_datetime.utc2unix(T0_utc)
        t_since_t0 = data['time']- t0

    if light_time!=0:
        start_utc= f'{start_utc} + {light_time:.02f} s'

    fig = plt.figure()
    for i,lc in data['lcs'].items():
        plt.plot(t_since_t0, lc, label=data['energy_bins']['names'][i])

    plt.xlabel(f'Start: {start_utc}')
    plt.ylabel('Counts / 4 s')
    title=f'STIX QL LCs (Flare #{flare_id})'
    plt.title(title)
    filename = os.path.join(folder, f'stix_flare_ql_lc_{_id}_{flare_id}.png')
    fig.tight_layout()
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.savefig(filename, dpi=100)
    print(filename)
    mdb.update_flare_joint_obs(_id, key, [filename])
    #print(filename)
    #plt.plot(t_since_t0,data['lc_smoothed'][where], label='1-min moving mean')
if __name__=='__main__':
    plot_stix_lc('.',0,0, '2021-03-01T08:50:00','2021-03-01T08:52:00')
