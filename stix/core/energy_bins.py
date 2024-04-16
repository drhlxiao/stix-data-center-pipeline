#A python class to manage STIX energy bins
#it needs to be maintained 
# History:
#    Created on 2023-08-02 by  Hualin Xiao 
import time

import numpy as np
import pandas as pd
from astropy.time import Time
from stix.core import config

#MISSION_END_TIME = 2524608000  # 2050-01-01:00:00

class StixEnergyBins:
    _default_ebins = config.instrument_config['default_ebins']

    _history = config.instrument_config['ebins_history']

    @staticmethod
    def get_ebin(ibin, unix_time=0):
        return StixEnergyBins.get_energy_bins(unix_time)[ibin]

    @staticmethod
    def get_energy_bins( unix_time=0):
        if unix_time == 0:
            return StixEnergyBins._default_ebins

        for entry in StixEnergyBins._history:
            if entry['time_range'][0] <= unix_time < entry['time_range'][1]:
                return entry['energy_bins']

        return StixEnergyBins._default_ebins

    @staticmethod
    def keV2sci(elow_keV, ehigh_keV, unix_time):
        ebins=StixEnergyBins.get_energy_bins(unix_time)
        ebins_low=np.array(ebins)[:-1]
        ebins_high=np.array(ebins)[1:]
        sel=np.where( (elow_keV<=ebins_low) &  (ebins_high<=ehigh_keV))
        return [np.min(sel), np.max(sel)]
    @staticmethod
    def sci2keV(elow_sci, ehigh_sci, unix_time):
        ebins=StixEnergyBins.get_energy_bins(unix_time)
        ebins_low=np.array(ebins)[:-1]
        ebins_high=np.array(ebins)[1:]
        return ebins_low[elow_sci], ebins_high[ehigh_sci]

    @staticmethod
    def get_elut_edge(offset, gain, energy):
        """
        compute elut edge 
        """
        return round(offset+ energy/gain)


    
    @staticmethod
    def get_emask_energy_bins(emask,unix_time):
        ebins=[]
        for i in range(32):
            if emask & (1 << i) != 0:
                ebins.append(i);
        names=[]
        energy_bins=StixEnergyBins.get_energy_bins(unix_time)
        sci_edges=[]
        for i in range(len(ebins) - 1):
            begin = ebins[i]
            end = ebins[i + 1]
            sci_edges.append([begin,end])
            if end == 32:
                names.append(f'{energy_bins[begin]} keV –⁠ inf')
            elif end < 32:
                names.append(f'{energy_bins[begin]}  – {energy_bins[end]} keV')
            else:
                names.append('')
        return {'names':names, 'sci_bin_edges':sci_edges}


    @staticmethod
    def get_energy_bins_dict(obs_time):
        if isinstance(obs_time, Time):
            t=obs_time.to_datetime()
        else:
            t=pd.to_datetime(obs_time,utc=True)
        unix_time=t.timestamp()
        ebins=StixEnergyBins.get_energy_bins(unix_time)
        res={}
        for i in range(32):
            res[i]={ 'channel_edge': i,
                    'energy_edge': ebins[i],
                    'e_lower': ebins[i],
                    'e_upper': ebins[i+1],
                    'bin_width': ebins[i+1]-ebins[i],
                    }
        
        return res

if __name__=='__main__':
    from pprint import pprint
    a=StixEnergyBins.get_energy_bins_dict('2023-01-01T00:00:00')
    b=StixEnergyBins.get_energy_bins_dict('2023-08-01T00:00:00')
    pprint(a)
    pprint(b)
    

