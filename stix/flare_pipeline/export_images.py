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
from astropy import units as u
import matplotlib.pyplot as plt
from stix.core import mongo_db as db
from stix.spice import time_utils as sdt
from stix.spice import solo 
from stix.core import logger
from astropy import units as u

mdb = db.MongoDB()
flare_image_db = mdb.get_collection('flare_images')
logger = logger.get_logger()
thermal_fig_path='/data/quicklook/images/thermal'
nonthermal_fig_path='/data/quicklook/images/nonthermal'

grid_spacing=5* u.deg
color='w'
def export_images(_id):
    """
        extract meta data from an image
    """
    print("processing:", _id)
    doc=flare_image_db.find_one({'_id':_id})
    try:
        filename=doc['fits']['image_clean']
        clean= sunpy.map.Map(filename)
    except Exception as e:
        return
    try:
        filename=doc['fits']['image_em']
        em= sunpy.map.Map(filename)
    except Exception as e:
        return
    erange=f"{doc['energy_range'][0]} – {doc['energy_range'][1]} keV"
    utc=f"{doc['utc_range'][0]}"
    fig=plt.figure(figsize=(9,4))
    ax = fig.add_subplot(121, projection=clean)
    clean.plot(cmap='std_gamma_2', axes=ax)
    clean.draw_grid(color=color, ls='--', grid_spacing=grid_spacing)
    clean.draw_limb(axes=ax, color=color, alpha=0.5)
    ax.set_aspect('equal')
    ax1 = fig.add_subplot(122, projection=em)
    em.plot(cmap='std_gamma_2', axes=ax1)

    em.draw_grid(color=color, ls='--', grid_spacing=grid_spacing)
    em.draw_limb(axes=ax1, color=color, alpha=0.5)
    if doc['energy_range'][1]<15:
        fname=os.path.join(thermal_fig_path, f'image_{doc["_id"]}.png')
    else:
        fname=os.path.join(nonthermal_fig_path, f'image_{doc["_id"]}.png')
    
    plt.suptitle(f'Image #{doc["_id"]} {erange} {utc}')
    print(fname)
    fig.savefig(fname)

    


if __name__=='__main__':
    if len(sys.argv[1])==3:
        print('export <start_id> <end_id>')
    else:
        for i in range(int(sys.argv[1]), int(sys.argv[2])):
            try:
                export_images(i)
            except Exception as e:
                print(e)


