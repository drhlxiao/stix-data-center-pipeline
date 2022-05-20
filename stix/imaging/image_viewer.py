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

SMALL_SIZE = 9
matplotlib.rc('font', size=SMALL_SIZE)
matplotlib.rc('axes', titlesize=SMALL_SIZE)
#matplotlib.rcParams['axes.titlepad']=20
#plt.rcParams['figure.constrained_layout.use'] = True

CMAP='std_gamma_2' #color map
GRID_COLOR='w'
def create_images(_id):
    doc=flare_images_db.find_one({'_id':_id})
    if doc:
        logger.info(f"Creating images for doc: {_id}...")
        try:
            plot_stix_images(doc)
        except Exception as e:
            logger.error(e)
    else:
        logger.info(f"doc {_id} does not exist")
def rotate_map(m, recenter=False):
    #further checks are required
    #rotate map
    return m.rotate(angle=(m.meta['crota2'])*u.deg, recenter=recenter)

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


def create_images_in_queue(num=10000):
    cursor=flare_images_db.find({'fits':{'$gt':{}}, 'figs.0':{'$exists':False}}).sort('_id',-1).limit(num)
    for doc in cursor:
        logger.info("Creating images for doc: {doc['_id']}...")
        try:
            plot_stix_images(doc)
        except Exception as e:
            logger.error(e)


def create_figures_ids_between(start_id, end_id):
    for i in range(start_id, end_id):
        doc=flare_images_db.find_one({'_id':i})
        if not doc:
            logger.warning(f"Failed to create figures for DocID:{i}")
            continue
        logger.info(f"Creating images for BSD#{doc['bsd_id']}, DocID:{i}")
        plot_stix_images(doc )

            

def plot_stix_images(doc ):
    """
    Plot STIX images and save them to image files
    Parameters

    """
    if not doc.get('fits',None):
        logger.error('Fits file not found')
        return
    # key should be like "image_xxx"

    name_map={'clean':'CLEAN', 
            'full':'Back-projection (full-disk)', 
            'fwd':'Forward-fit',
            'em':'EM', 
            'bp':'Back-projection'}

    ar=u.def_unit("arcsecs",1*u.arcsec)
    am=u.def_unit("meters",1*u.m)
    ad=u.def_unit("degrees",1*u.deg)
    u.add_enabled_units([ar, am, ad])
    #define unit arcsecs is the same arcsec

    maps, map_names=[],[]
    doc_fits=doc['fits']
    find_key=lambda x: [k for k in keys if x in k]

    fig = plt.figure(figsize=(12,7), dpi=100, facecolor='white', constrained_layout=True)
    ax_lc= fig.add_subplot(231)
    bsd_id=doc['bsd_id']
    energy_range=doc['energy_range']
    start_utc, end_utc=doc['utc_range']
    bsd_uid=doc['unique_id']
    bkg_uid=doc['background']['unique_id']

    lightcurves.plot_QL_lc_for_bsd(bsd_id, fill_between_times=[start_utc, end_utc], ax=ax_lc)

    text_xy=[0.02,0.95]


    mbp_full=sunpy.map.Map(doc_fits['image_full_disk'])

    mbp_full= rotate_map(mbp_full, recenter=True)

    ax_mbp_full= fig.add_subplot(232, projection=mbp_full)

    mbp_full.plot(vmin=mbp_full.max()*0.2, 
                    cmap=CMAP,
                    axes=ax_mbp_full, title="")#  title="Back-projection (full)",)

    mbp_full.draw_grid(color=GRID_COLOR, ls='--', grid_spacing=10*u.deg)
    mbp_full.draw_limb(axes=ax_mbp_full, color=GRID_COLOR,alpha=0.5)

    ax_mbp_full.text(text_xy[0], text_xy[1],'Back-projection (full-disk)',
     horizontalalignment='left',   verticalalignment='center',    transform = ax_mbp_full.transAxes, color=GRID_COLOR)





    mbp=sunpy.map.Map(doc_fits['image_bp'])

    mbp= rotate_map(mbp)

    ax_mbp= fig.add_subplot(233, projection=mbp)
    mbp.plot(  cmap=CMAP,
                    axes=ax_mbp,  title="")

    mbp.draw_grid(color=GRID_COLOR, ls='--', grid_spacing=10*u.deg)
    mbp.draw_limb(axes=ax_mbp, color=GRID_COLOR,alpha=0.5)

    ax_mbp.text(text_xy[0], text_xy[1],'Back-projection',
     horizontalalignment='left',   verticalalignment='center',    transform = ax_mbp.transAxes, color=GRID_COLOR)

    
    levels=np.array([0.5])
    cs=mbp.draw_contours(levels*100*u.percent)
    clevels =levels*mbp.max()
    plt.clabel(cs, inline=1,  fmt={x: f'{levels[i]*100:.0f} %'    for i, x in enumerate(clevels) })


    fix_clean_map_fits_header(doc, doc_fits['image_clean'])
    #add dsun to the header


    mclean=sunpy.map.Map(doc_fits['image_clean'])


    mclean=mclean[4]
    mclean= rotate_map(mclean)

    ax_clean= fig.add_subplot(234, projection=mclean)
    mclean.plot( cmap=CMAP,
                    axes=ax_clean, title="")

    mclean.draw_grid(color=GRID_COLOR, ls='--', grid_spacing=10*u.deg)
    mclean.draw_limb(axes=ax_clean, color=GRID_COLOR,alpha=0.5)

    ax_clean.text(text_xy[0], text_xy[1],'CLEAN',
     horizontalalignment='left',   verticalalignment='center',    transform = ax_clean.transAxes, color=GRID_COLOR)


    mem=sunpy.map.Map(doc_fits['image_em'])

    mem= rotate_map(mem)

    ax_em= fig.add_subplot(235, projection=mem)
    mem.plot( cmap=CMAP, axes=ax_em, title="")


    mem.draw_grid(color=GRID_COLOR, ls='--', grid_spacing=10*u.deg)
    mem.draw_limb(axes=ax_em, color=GRID_COLOR,alpha=0.5)
    ax_em.text(text_xy[0], text_xy[1],'EM',
     horizontalalignment='left',   verticalalignment='center',    transform = ax_em.transAxes, color=GRID_COLOR)


    try:
        fwdshape=f"({doc['idl_config']['fwdfit_shape']})"
    except (KeyError, TypeError):
        fwdshape=''
    mfwd=sunpy.map.Map(doc_fits['image_fwdfit'])

    mfwd= rotate_map(mfwd)

    ax_fwd= fig.add_subplot(236, projection=mfwd)
    mfwd.plot( cmap=CMAP, axes=ax_fwd, title="")

    mfwd.draw_grid(color=GRID_COLOR, ls='--', grid_spacing=10*u.deg)
    mfwd.draw_limb(axes=ax_fwd, color=GRID_COLOR,alpha=0.5)

    ax_fwd.text(text_xy[0], text_xy[1],f'Forward-fit {fwdshape}',
     horizontalalignment='left',   verticalalignment='center',    transform = ax_fwd.transAxes, color=GRID_COLOR)


    img_fname=os.path.join(doc['idl_config']['folder'],f'{doc["idl_config"]["prefix"]}.png')
    pdf_img_fname=os.path.join(doc['idl_config']['folder'],f'{doc["idl_config"]["prefix"]}.pdf')
    image_id_str = f'(#{doc["_id"]})'

    duration=ut.utc2unix(end_utc)-ut.utc2unix(start_utc)


    plt.suptitle(f'Start UTC {start_utc}; Exp. time: {duration:.1f} sec; Energy range: {energy_range[0]} – {energy_range[1]} keV ', fontsize=12)
    plt.subplots_adjust(
            top=0.85,
                    wspace=0.2, 
                    hspace=0.4)
    plt.savefig(img_fname, format='png')

    plt.savefig(pdf_img_fname, format='pdf')
    logger.info(f"Images have been written to file:{img_fname}")

    updates={'$set':{'figs':[{'quicklook':img_fname}], 'processing_date':datetime.now()}}
    flare_images_db.update_one({'_id':doc['_id']}, updates)



def test_image_viewer():
    create_figures_ids_between(0,1)
if __name__=='__main__':
    if len(sys.argv)==1:
        create_images_in_queue()
    elif len(sys.argv)==2:
        create_images(int(sys.argv[1]))
    elif len(sys.argv)==3:
        for i in range(int(sys.argv[1]), int(sys.argv[2])):
            print("Creating images for Entry #",i)
            create_images(i)
