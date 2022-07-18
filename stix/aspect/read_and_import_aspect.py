"""
    Read stix auxiliary fits files and import the data to NoSQL
    Hualin Xiao July 18, 2022
"""
import os
import numpy as np
from astropy.io import fits
from datetime import datetime


from stix.spice import time_utils as sdt
from stix.core import logger
logger = logger.get_logger()

mdb = db.MongoDB()

def read_fits_to_dict(fname):
    """
    read fits file and convert data to dict
    """
    hdul = fits.open(fname)
    data=hdul['DATA'].data
    res=[]
    now=datetime.now()
    keys=data.dtype.names
    for i in range(data.size):
        row={}
        for key in  keys:
            val=data[i][key]
            if isinstance(val, np.ndarray):
                val=val.tolist()
            row[key]=val
        row['filename']=os.path.basename(fname)
        row['unix_time']=sdt.utc2unix(row['time_utc'])
        row['import_date']=now
        res.append(row)
    return res

def import_auxiliary(fname):
    db=mdb.get_collection('aspect')
    if db.find_one({'filename':os.path.basename(fname)}):
        logger.warn(f'{fname} ignored because it has been imported!')
        return
    rows=read_fits_to_dict(fname)
    num=db.insert_many(rows)
    logger.info(f"Inserted entries {len(row)}, File: {fname}")

if __name__=='__main__':
    if len(sys.argv)!==1:
        logger.info('read_and_import_aspect <filename')
    else:
        import_auxiliary(sys.argv[1])


