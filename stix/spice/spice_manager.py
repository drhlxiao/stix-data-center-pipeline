#!/usr/bin/python3

import re
import os
import sys
sys.path.append('./')
import glob
import spiceypy
from pathlib import Path
from datetime import datetime
from dateutil import parser as dtparser
from astropy.time import Time

from stix.core import stix_logger
from stix.core import config
from stix.core import mongo_db 

NUM_KERNEL_FILES_LIMIT=10

logger = stix_logger.get_logger()

mdb=mongo_db.MongoDB()
# SOLAR ORBITER naif identifier
class SpiceManager:
    """taken from https://issues.cosmos.esa.int/solarorbiterwiki/display/SOSP/Translate+from+OBT+to+UTC+and+back
    """
    def __init__(self):

        self.version_date=datetime.strptime('19700101', "%Y%m%d")

        self.loaded_kernel_filename=None
        self.latest_mk=None
        self.load_kernels()

    def get_kernel_filename(self):
        return self.loaded_kernel_filename
    def load_kernels(self):
        spice_folder=config.get_config('spice')
        mk_folder=os.path.join(spice_folder,'mk')
        os.chdir(mk_folder)
        #pred_mk='solo_ANC_soc-pred-mk.tm'
        #print(mk_folder)
        self.latest_mk=None
        for filename in glob.glob(f'{mk_folder}/solo_ANC_soc-flown-mk*.tm'):
            date_str=re.findall(r"\d{4}\d{2}\d{2}", filename)
            if date_str:
                fdt=datetime.strptime(date_str[0], "%Y%m%d")
                if fdt>self.version_date:
                    self.latest_mk=os.path.basename(filename)
                    self.version_date=fdt
        #if latest_mk!=None and utc<self.version_date:
        
        if self.loaded_kernel_filename !=  self.latest_mk and self.latest_mk !=None:
            spiceypy.furnsh(self.latest_mk)
            self.loaded_kernel_filename=self.latest_mk
        else:
            print(f'Skipped! {self.latest_mk} has been loaded!')


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


