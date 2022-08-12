#!/usr/bin/python
"""
    Import data from STIX data archive
    Hualin Xiao July 18, 2022
"""
import os
import re
import sys

sys.path.append('/opt/stix/parser/')
import re
import glob
import numpy as np
from astropy.io import fits
from datetime import datetime

from stix.core import mongo_db as db
from stix.utils import bson

from stix.spice import time_utils as sdt
from stix.core import logger

logger = logger.get_logger()

mdb = db.MongoDB()

DEFAULT_ASP_PATH_PATTEN = "/data/pub099/fits/L2/*/*/*/*/*aux*.fits"
DATA_ARCHIVE_FITS_PATH = "/data/pub099/fits/**/*.fits"
processed_list = []

DATA_ARCHIVE_FILE_INFO = {
    'L1_stix-ql-background': ('ql_background', 'L1', 'quicklook', 54119),
    'L1_stix-ql-flareflag': ('flareflag_location', 'L1', 'quicklook', 54122),
    'L1_stix-ql-lightcurve': ('ql_light_curves', 'L1', 'quicklook', 54118),
    'L1_stix-ql-spectra': ('ql_spectrogram', 'L1', 'quicklook', 54120),
    'L1_stix-ql-variance': ('ql_variance', 'L1', 'quicklook', 54121),
    'L1_stix-hk-maxi': ('hk_maxi', 'L1', 'housekeeping', 54102),
    'L1_stix-hk-mini': ('hk_mini', 'L1', 'housekeeping', 54101),
    'L1_stix-cal-energy': ('calibration_spectrum', 'L1', 'quicklook', 54124),
    'L1_stix-sci-xray-cpd': ('xray_l1_user', 'L1', 'science', 54115),
    'L1_stix-sci-xray-spec': ('spectrogram_user', 'L1', 'science', 54143),
    'L1_stix-sci-aspect-burst': ('aspect', 'L1', 'auxiliary', 54125),
    'L2_stix-aux-auxiliary': ('auxiliary', 'L2', 'auxiliary', 54102),
    'L2_stix-hk-maxi': ('hk_maxi', 'L2', 'housekeeping', 54102)
}


def read_data_archive_fits_meta(filename, file_info):
    meta = {}
    try:
        hdu = fits.open(filename)
    except (OSError, FileNotFoundError):
        logger.warn(
            f'Could not open FITS file:{filename}, it may be corrupted!')
        return None
    try:
        date_range = (hdu[0].header['DATE-BEG'], hdu[0].header['DATE-END'])
    except KeyError:
        logger.warn(
            f'Could not read start time and end from FITS file:{filename}')
        return None
    creation_time = hdu[0].header.get('DATE', None)
    creation_time = sdt.utc2datetime(
        creation_time) if creation_time is not None else datetime.now()
    meta = {
        '_id': mdb.get_next_fits_id(),
        'filename': os.path.basename(filename),
        'path': os.path.dirname(filename),
        'complete': True,
        'level': file_info[1],
        'packet_spid': file_info[3],
        'file_id': '',
        'product_type': file_info[0],
        'product_group': file_info[2],
        'file_size': os.path.getsize(filename),
        'creation_time': creation_time,
        'data_start_unix': sdt.utc2unix(date_range[0]),
        'data_end_unix': sdt.utc2unix(date_range[1]),
    }

    if 'L1_stix-sci' in filename:
        #extract request id from filename
        request_id = [
            int(x) for x in re.findall(r"V\d\d_(\d+)-\d+.fits", filename)
        ]
        if request_id:
            meta['request_id'] = request_id[0]
    return meta


def import_data_archive_products(path=DATA_ARCHIVE_FITS_PATH):
    logger.info(f'Checking folder:{path}...')
    fits_db = mdb.get_collection('fits')
    file_types = DATA_ARCHIVE_FILE_INFO.keys()
    for fname in glob.glob(path, recursive=True):
        print("FILENAME:", fname)
        file_type = [t for t in file_types if t in fname]
        basename = os.path.basename(fname)
        if not file_type:
            logger.info(f'{basename} not supported type. Skipped!')
            continue
        if fits_db.find_one({'filename': basename}):
            logger.info(f'{basename} already inserted in the database')
            continue

        file_type = file_type[0]
        logger.info(f'Add {basename} to fits file database')
        meta = read_data_archive_fits_meta(fname,
                                           DATA_ARCHIVE_FILE_INFO[file_type])
        if meta:
            logger.info(f'inserting metadata for {basename} to fits_db..')
            fits_db.insert_one(meta)
            #update if
        if 'L2_stix-aux-auxiliary' in fname:
            import_auxiliary(fname)


def import_all_aspect_solutions(path_patten=DEFAULT_ASP_PATH_PATTEN):
    """
        Import all aspect solutions to database
    """

    for fname in glob.iglob(path_patten):
        if fname not in processed_list:
            import_auxiliary(fname)
            processed_list.append(fname)
            #reduce checking of database


def read_fits_to_dict(fname):
    """
    read fits file and convert data to dict
    """
    hdul = fits.open(fname)
    data = hdul['DATA'].data
    res = []
    now = datetime.now()
    keys = data.dtype.names
    for i in range(data.size):
        row = {key: data[i][key] for key in keys}
        row['filename'] = os.path.basename(fname)
        try:
            row['unix_time'] = sdt.utc2unix(row['time_utc'])
        except KeyError:
            logger.warn(f'key time_utc  not found in aux file:{fname} !')
            continue
        res.append(row)
    return bson.dict_to_json(res)


def import_auxiliary(fname):
    logger.info(f'Importing aspect solutions from {basename} ...')
    asp_db = mdb.get_collection('aspect')
    if asp_db.find_one({'filename': os.path.basename(fname)}):
        logger.warn(f'{fname} ignored because it has been imported!')
        return
    rows = read_fits_to_dict(fname)
    if rows:
        num = asp_db.insert_many(rows)
        logger.info(f"Inserted entries {len(rows)}, File: {fname}")
    else:
        logger.warn(f'{fname} contains empty valid entries!')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        logger.info('read_and_import_aspect <filename')
        import_data_archive_products()
    else:
        import_auxiliary(sys.argv[1])
