#!/usr/bin/env python3  
# -*- coding: utf-8 -*- 
#----------------------------------------------------------------------------
# Created By  : Hualin Xiao (hualin.xiao)
# Created Date: Sept 5, 2022
# version =1.0
"""
    Various tools to manipulate stix images
    python version must be greater than  3.9
"""
import os
import sys
import sunpy
import sunpy.map

import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy import units as u
import matplotlib.pyplot as plt
from photutils.aperture import CircularAperture
from photutils.detection import DAOStarFinder
from shapely.geometry import Polygon, Point

from stix.core import mongo_db as db
from stix.spice import time_utils as sdt
from stix.spice import solo 
from stix.core import logger
from stix.utils import bson

mdb = db.MongoDB()
flare_image_db = mdb.get_collection('flare_images')
logger = logger.get_logger()
FIG_PATH='/data/quicklook/images/'

def export_images(doc, fig_path=FIG_PATH)
    """
        extract meta data from an image
    """
    try:
        filename=doc['fits']['image_clean']
        clean= sunpy.map.Map(filename)
    except (TypeError, KeyError, IOError):
        return
    try:
        filename=doc['fits']['image_em']
        em= sunpy.map.Map(filename)
    except (TypeError, KeyError, IOError):
        return
    erange=f"{doc['energy_range'][0]} – {doc['energy_range'][1]} keV"
    utc=f"{doc['utc_range'][0]}"
    fig=plt.figure(figsize=(9,4))
    ax = fig.add_subplot(121, projection=clean)
    clean.plot(cmap='std_gamma_2', axes=ax)
    ax.set_aspect('equal')
    ax1 = fig.add_subplot(122, projection=em)
    clean.plot(cmap='std_gamma_2', axes=ax1)
    plt.suptitle(f'{doc["_id"]} {erange} {utc}')
    fname=os.path.join(fig_path, f'image_{doc["_id"]}.png')
    print(fname)
    fig.savefig(fname)

    


if __name__=='__main__':
    ignore_existing=False
    query={'fits.image_clean':{'$exists':True}}
    for doc in flare_image_db.find(query).sort('_id', 1):
        try:
            export_images(doc)
        except Exception as e:
            logger.error(str(e))


