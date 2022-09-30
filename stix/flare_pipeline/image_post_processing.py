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

def extract_image_meta(flare_image_doc, which_image='image_em',
        contour_level=50,  savefig=True, fig_path=None, ignore_existing=True):
    """
        extract meta data from an image
    """
    if  ignore_existing:
        image_meta=flare_image_doc.get('image_meta',None)
        if image_meta:
            logger.info(f'Imaging run  {flare_image_doc["_id"]} processed, ignore')
            return
    
    try:
        filename=flare_image_doc['fits'][which_image]
        imap = sunpy.map.Map(filename)
    except (TypeError, KeyError, IOError):
        #not pixel data, data not ready
        logger.error(f'Imaging run  {flare_image_doc["_id"]} processed, ignore')
        return

    mean, median, std = sigma_clipped_stats(imap.data, sigma=10)
    daofind = DAOStarFinder(fwhm=7, threshold=3. * std) # 3 sigma above standard deviation
    sources = daofind(imap.data - median)
    if not sources:
        logger.info(f'No source found for #: {flare_image_doc["_id"]}')
        return


    peak_pix_coords = np.transpose((sources['xcentroid'], sources['ycentroid']))


    peak_wcs_coords=[]
    flare_utc=sdt.unix2utc(flare_image_doc['start_unix'])
    flare_aux=None
    for pos in peak_pix_coords: #iterate over peaks
        peak_coord=imap.pixel_to_world(pos[0]*u.pix, pos[1]*u.pix)
        peak_coord_wcs=(peak_coord.Tx.value, peak_coord.Ty.value)
        if flare_aux is None:
            flare_aux =solo.SoloEphemeris.get_flare_spice(peak_coord_wcs[0], peak_coord_wcs[1], flare_utc, observer='earth')
        peak_wcs_coords.append(peak_coord_wcs)

    clevels = np.array([contour_level])*u.percent

    fig=plt.figure(figsize=(9,4))

    ax = fig.add_subplot(121, projection=imap)
    imap.plot(cmap='std_gamma_2', axes=ax, title="")
    imap.draw_grid(color='w', ls='--', grid_spacing=10 * u.deg)
    imap.draw_limb(axes=ax, color='w', alpha=0.5)

    ax.set_aspect('equal')
    ax1 = fig.add_subplot(122, projection=imap)

    cs=imap.draw_contours(clevels, axes=ax1)
    apertures = CircularAperture(peak_pix_coords, r=20.)
    norm = ImageNormalize(stretch=SqrtStretch())
    ax1.imshow(imap.data, cmap='Greys', origin='lower', norm=norm,
               interpolation='nearest')
    apertures.plot(color='blue', lw=1.5, alpha=1, axes=ax1)


    meta={'sources':[]}

    for i in range(len(clevels)):
        contour = cs.collections[i]
        if not contour:
            logger.warn(f'No contour found for #: {flare_image_doc["_id"]}')
            continue
        paths=contour.get_paths()
        for i, p in enumerate(paths):
            #number of levels
            vs_trans = p.vertices.T
            
            wcs_vtx = imap.pixel_to_world(vs_trans[0]*u.pix, vs_trans[1]*u.pix)
            vtx=np.vstack((wcs_vtx.Tx.value, wcs_vtx.Ty.value)).T

            if vtx.size<3:
                logger.warn(f'No contour found for#: {flare_image_doc["_id"]}')
                continue


            pgon = Polygon(vtx)
     
            meta['sources'].append({'contour':{
                'coord':vtx,  
                'image': which_image, 
                'level':contour_level,
                'rsun':flare_image_doc['aux']['rsun'],
                'area':pgon.area
                },
                'peaks':[p_wcs for p_wcs in peak_wcs_coords if pgon.contains(Point(p_wcs))]
                    })

            
            ax1.plot(vs_trans[0],vs_trans[1], label=f'area:{pgon.area:.2f}')
    ax1.legend()
    ax1.set_aspect('equal')
    if fig_path is None:
        fig_path = flare_image_doc["idl_config"]["folder"]
    fname=os.path.join(fig_path, f'image_sources_{flare_image_doc["_id"]}.png')
    plt.suptitle(f'{flare_image_doc["_id"]}')
    fig.savefig(fname)
    logger.info(f'Writing image filename:{fname}')
    meta['flare_aux']=flare_aux
    meta_json = bson.dict_to_json(meta)
    flare_image_db.update_one({'_id':flare_image_doc['_id']},{'$set':{'image_meta': meta_json}})

    


if __name__=='__main__':
    ignore_existing=False
    query={'figs.0':{'$exists':True}, 'image_meta':{'$exists':False}} if ignore_existing else {'figs.0':{'$exists':True}}
    for doc in flare_image_db.find(query).sort('_id', -1):
        try:
            extract_image_meta(doc, ignore_existing=ignore_existing)
        except Exception as e:
            logger.error(str(e))


