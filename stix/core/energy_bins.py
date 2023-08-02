#A python class to manage STIX energy bins
#it needs to be maintained 
# History:
#    Created on 2023-08-02 by  Hualin Xiao 
import time

import numpy as np
import pandas as pd
from astropy.time import Time

MISSION_END_TIME = 2524608000  # 2050-01-01:00:00

class StixEnergyBins:
    _default_ebins = [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70,
                      76, 84, 100, 120, 150, np.inf]

    _history = [
        {
            'time_range': [1690761600, MISSION_END_TIME],  # 2023-07-31T00:00:00 to sometime in the future
            'energy_bins': [0.0, 4.0, 5.5, 6.3, 6.56, 6.8, 7.3, 8.0, 9.0,
                            10.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0,
                            28.0, 32.0, 36.0, 40.0, 45.0, 50.0, 56.0, 63.0, 70.0, 76.0, 84.0, 100.0, 120.0, 150.0, np.inf]
        },
        {
            'time_range': [1675036800, 1675353600],  # 2023-01-30 to 2023-02-02T16:00
            'energy_bins': [0, 4, 4.45, 4.95, 5.45, 5.95, 6.45, 6.95, 7.35, 7.75, 8.25, 8.75, 9.25, 10,
                            10.5, 11, 12, 13, 15, 18, 21, 25, 28, 32, 36, 43, 50, 59, 70, 84, 110, 150, np.inf]
        }
    ]

    def get_ebin(ibin, unix_time=0):
        return StixEnergyBins.get_energy_bins(unix_time)[ibin]

    def get_energy_bins( unix_time=0):
        if unix_time == 0:
            return StixEnergyBins._default_ebins

        for entry in StixEnergyBins._history:
            if entry['time_range'][0] <= unix_time < entry['time_range'][1]:
                return entry['energy_bins']

        return StixEnergyBins._default_ebins

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
    

