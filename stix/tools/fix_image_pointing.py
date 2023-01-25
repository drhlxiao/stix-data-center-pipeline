#!/usr/bin/python
"""
Image viewer, create image from fits files which are created by idl imaging software 

April 27, 2022
"""
import os
import sys
import matplotlib
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.io import fits
from matplotlib import pyplot as plt
from dateutil.parser import parse as dtparse
from stix.core import logger
from stix.spice import time_utils as ut

from datetime import datetime

import sunpy
from sunpy.map import make_fitswcs_header

from sunpy.coordinates.frames import HeliocentricEarthEcliptic,HeliographicStonyhurst

from stix.flare_pipeline import lightcurves
from stix.core import mongo_db as db

logger = logger.get_logger()
mdb = db.MongoDB()
flare_images_db= mdb.get_collection('flare_images')


def fix_clean_map_fits_header(doc, image_filename):
    """
        to add dsun to clean maps
        idl script dosen't write dsun to clean fits files, we do a fix here 
        
    """
    hduls=fits.open(image_filename) 
    dsun=(doc['aux']['dsun'], 'S/C distance to Sun (meters)')
    for hdu in hduls:
        hdu.header['DSUN_OBS']=dsun
        hdu.header['DSUN']=dsun
    hduls.writeto(image_filename, overwrite=True, checksum=True)
    logger.info(f"Adding more keywords to {image_filename}")


