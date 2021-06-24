#!/usr/bin/python3

import re
import os
import glob
from datetime import datetime
from dateutil import parser as dtparser
from astropy.time import Time

import spiceypy
from stix.core import stix_logger
from stix.core import config
from stix.core import mongo_db 

logger = stix_logger.get_logger()

MDB=mongo_db.MongoDB()

loaded_kernels=[]


# SOLAR ORBITER naif identifier
class SpiceManager:
    """taken from https://issues.cosmos.esa.int/solarorbiterwiki/display/SOSP/Translate+from+OBT+to+UTC+and+back
    """
    def __init__(self):
        self.loaded_kernels=[]
        self.last_sclk_file=None
        self.refresh_kernels()
        
    def refresh_kernels(self):
        for pattern in config.get_spice():
            for fname in glob.glob(pattern):
                if fname not in self.loaded_kernels:
                    self.loaded_kernels.append(fname)
                    MDB.insert_spice_kernel(fname)
        latest_kernels=MDB.get_latest_spice_kernels()
        for kernel_type,kernels in latest_kernels.items():
            for kernel in kernels:
                fname=os.path.join(kernel['path'],kernel['filename'])
                try:
                    spiceypy.furnsh(fname)
                    if 'sclk' in fname:
                        self.last_sclk_file=fname
                except spiceypy.utils.exceptions.SpiceNOSUCHFILE:
                    #print(f'Failed to load {fname}')
                    pass

        

    def get_last_sclk_filename(self):
        return self.last_sclk_file



    def obt2utc(self, obt_string):
        # Obt to Ephemeris time (seconds past J2000)
        ephemeris_time = spiceypy.scs2e(-144, obt_string)
        # Ephemeris time to Utc
        # Format of output epoch: ISOC (ISO Calendar format, UTC)
        # Digits of precision in fractional seconds: 3
        return spiceypy.et2utc(ephemeris_time, "ISOC", 3)

    def utc2obt(self, utc_string):
        # Utc to Ephemeris time (seconds past J2000)
        ephemeris_time = spiceypy.utc2et(utc_string)
        # Ephemeris time to Obt
        #return ephemeris_time
        obt_string = spiceypy.sce2s(-144, ephemeris_time)
        time_fields = re.search('\/(.*?):(\d*)', obt_string)
        group = time_fields.groups()
        try:
            return int(group[0]) + int(group[1]) / 65536.
        except Exception as e:
            logger.warning(str(e))
            return 0

    def scet2utc(self, coarse, fine=0):
        obt_string = '{}:{}'.format(coarse, fine)
        #print(obt_string)
        return self.obt2utc(obt_string)

    def utc2scet(self, utc):
        # Utc to Ephemeris time (seconds past J2000)
        ephemeris_time = spiceypy.utc2et(utc)
        # Ephemeris time to Obt
        return spiceypy.sce2s(-144, ephemeris_time)


spice= SpiceManager()


