#this is script is to the fits files, packets and QL data
# hualin.xiao created on May. 24, 2022
# The issue is explained here:  https://www.cs.technik.fhnw.ch/elog/Operations/2

import os
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.table.table import QTable
from datetime import datetime

from stix.core import mongo_db as db
from stix.core import datatypes as sdt
from stix.spice import time_utils as st

from pathlib import Path

mdb = db.MongoDB()
qb_bsd = mdb.get_collection("bsd")
qb_pkts = mdb.get_collection("packets")
db_ql = mdb.get_collection("packets")
db_fits = mdb.get_collection("fits")

start_utc = '2023-01-30T00:00:00'
end_utc = '2023-02-02T16:00:00'
start_unix = st.utc2unix(start_utc)
end_unix = st.utc2unix(end_utc)
elow = [
    0, 4, 4.45, 4.95, 5.45, 5.95, 6.45, 6.95, 7.35, 7.75, 8.25, 8.75, 9.25, 10,
    10.5, 11, 12, 13, 15, 18, 21, 25, 28, 32, 36, 43, 50, 59, 70, 84, 110, 150
]
ehigh = [
    4, 4.45, 4.95, 5.45, 5.95, 6.45, 6.95, 7.35, 7.75, 8.25, 8.75, 9.25, 10,
    10.5, 11, 12, 13, 15, 18, 21, 25, 28, 32, 36, 43, 50, 59, 70, 84, 110, 150,
    np.inf
]
#elow=[  0.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,  11.,  12.,  13.,
#        14.,  15.,  16.,  18.,  20.,  22.,  25.,  28.,  32.,  36.,  40.,
#        45.,  50.,  56.,  63.,  70.,  76.,  84., 100., 120., 150.]
#ehigh=[  4.,   5.,   6.,   7.,   8.,   9.,  10.,  11.,  12.,  13.,  14.,
#        15.,  16.,  18.,  20.,  22.,  25.,  28.,  32.,  36.,  40.,  45.,
#        50.,  56.,  63.,  70.,  76.,  84., 100., 120., 150.,  np.inf]

def correct_energy_bins(fname):
    hdul = fits.open(fname)
    energies = QTable()
    energies['channel'] = range(len(elow))
    energies['e_low'] = elow * u.keV
    energies['e_high'] = ehigh * u.keV
    energy_hdu = fits.BinTableHDU(energies)
    energy_hdu.name = 'ENERGIES'
    hdul['PRIMARY'].header['remarks'] = 'Fine energy binning test'
    hdul['ENERGIES'] = energy_hdu
    print("Correcting file:", fname)
    hdul.writeto(fname, overwrite=True, checksum=True)


def correct_FITS():
    query = {
        "data_start_unix": {
            '$gt': start_unix,
            '$lte': end_unix
        },
        'data_end_unix': {
            '$gt': start_unix,
            '$lte': end_unix
        },
        'product_group': 'science',
        'packet_spid': {
            '$in': [54115, 54143]
        }
    }
    fits_docs = db_fits.find(query)
    for doc in fits_docs:
        fname = os.path.join(doc['path'], doc['filename'])
        correct_energy_bins(fname)


correct_FITS()
