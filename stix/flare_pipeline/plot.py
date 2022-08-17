#!/usr/bin/python
"""
Image viewer, create image from fits files which are created by idl imaging software 

April 27, 2022
"""
import os
import sys
import matplotlib
from matplotlib import cm
import numpy as np
from stix.imaging.images import plot

logger = logger.get_logger()
mdb = db.MongoDB()
flare_images_db= mdb.get_collection('flare_images')

SMALL_SIZE = 9
PDF_FIGURE_SIZE=(17, 7)


matplotlib.rc('font', size=SMALL_SIZE)
matplotlib.rc('axes', titlesize=SMALL_SIZE)
#matplotlib.rcParams['axes.titlepad']=20
#plt.rcParams['figure.constrained_layout.use'] = True
def plot(_id):
    doc=flare_images_db.find_one({'_id':_id})
    if doc:
        logger.info(f"Creating images for doc: {_id}...")
        try:
            plot_imaging_results(doc)
        except Exception as e:
            logger.error(e)
    else:
        logger.info(f"doc {_id} does not exist")
        
