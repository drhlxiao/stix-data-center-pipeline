from astropy.io import fits
from stix.spice import time_utils as sdt
from stix.spice import spice_manager as spm
from astropy.time import Time
import glob
def update_header(fname):
    hdul=fits.open(fname) 
    try:
        rsun=hdul['PRIMARY'].header['RSUN_ARC']
        print('Ignore:'+fname)
    except KeyError as e:
        try:
            mid=hdul['PRIMARY'].header['DATE_AVG']
            start=hdul['PRIMARY'].header['DATE_BEG']
        except KeyError:
            print('[ ERROR ]NO DATE info: '+fname)
            return None
    
        start_dt=Time(start)
        center_dt=Time(mid)
        try:
            eph_header= spm.spice.get_fits_headers(start_time=start_dt,
                                                average_time=center_dt)
        
        except Exception as e:
            print('[ ERROR ]No ephme data'+fname)
            return None

        print('[ INFO1 ]updating '+fname)
        hdul['PRIMARY'].header.update(eph_header)
        hdul.writeto(fname, overwrite=True, checksum=True)
        #hdul.writeto(full_path, overwrite=True, checksum=True)

for fname in glob.glob('/data/fits/*.fits'):

    update_header(fname)
