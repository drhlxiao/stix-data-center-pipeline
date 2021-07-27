#!/usr/bin/python3

import re
import os
import glob
from pathlib import Path
from datetime import datetime
from dateutil import parser as dtparser
from astropy.time import Time

import spiceypy
from stix.core import stix_logger
from stix.core import config
from stix.core import mongo_db 

NUM_KERNEL_FILES_LIMIT=10

logger = stix_logger.get_logger()

mdb=mongo_db.MongoDB()
collection_spice=mdb.get_collection('spice')
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
                    self.insert_spice_kernel(fname)
        latest_kernels=self.get_latest_spice_kernels()
        for kernel_type,kernels in latest_kernels.items():
            for kernel in kernels:
                fname=os.path.join(kernel['path'],kernel['filename'])
                print(f'Loading type {kernel_type}: {fname}')
                try:
                    spiceypy.furnsh(fname)
                    if 'sclk' in fname:
                        self.last_sclk_file=fname
                except spiceypy.utils.exceptions.SpiceNOSUCHFILE:
                    print(f'Failed to load {fname}')
                    pass

    def get_latest_spice_kernels(self):
        types=list(collection_spice.distinct('type'))
        results={}
        for key in types:
            num=1 if key in ['sclk','spk'] else NUM_KERNEL_FILES_LIMIT
            sort_field='file_date' if key != 'spk' else 'counter'
            rows=list(collection_spice.find({'type':key}).sort(sort_field, -1).limit(num))

            if rows:
                rows.reverse()
                results[key]=rows
        return results 
    def get_spice_kernels(self):
        return collection_spice.find({}).sort('file_date',1)

    def insert_spice_kernel(self, filename):
        next_id=0
        doc={}
        pfilename=Path(filename)
        fdate=re.findall(r"\d{4}\d{2}\d{2}", filename)
        date_str=fdate[0] if fdate else '19700101'
        if 'orbit' in filename:
            ver=re.findall(r"_V(\d)_", filename)
            ver2=re.findall(r"_V(\d{2})_", filename)
            ltp=re.findall(r"_L(\d{3})_", filename)
            counter=re.findall(r"_(\d{5})_", filename)
            try:
                doc['LTP']=int(ltp[0])
            except:
               pass
            try:
                doc['version']=int(ver[0])
            except:
                pass
            try:
                doc['version2']=int(ver2[0])
            except:
                pass
            try:
                doc['counter']=int(counter[0])
            except:
                pass
        doc['file_date']=datetime.strptime(date_str, "%Y%m%d")
        doc['path']=pfilename.parent.as_posix()
        doc['filename']=pfilename.name
        doc['type']=pfilename.parent.name
        if collection_spice:
            is_found=collection_spice.find_one({'filename':doc['filename']})
            if is_found:
                return
        try:
            next_id=collection_spice.find({}).sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            pass

        doc['_id']=next_id
        doc['entry_time']=datetime.now()
        collection_spice.save(doc)
        

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


