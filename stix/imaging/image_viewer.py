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
from stix.analysis.science_l1 import ScienceL1

from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime
from scipy import ndimage

import sunpy
from sunpy.map import make_fitswcs_header

from sunpy.coordinates.frames import HeliocentricEarthEcliptic,HeliographicStonyhurst
from stix.flare_pipeline import lightcurves
from stix.core import mongo_db as db

logger = logger.get_logger()
mdb = db.MongoDB()
flare_images_db= mdb.get_collection('flare_images')

SMALL_SIZE = 9
PDF_FIGURE_SIZE=(17, 6)


matplotlib.rc('font', size=SMALL_SIZE)
matplotlib.rc('axes', titlesize=SMALL_SIZE)
#matplotlib.rcParams['axes.titlepad']=20
#plt.rcParams['figure.constrained_layout.use'] = True

CMAP='std_gamma_2' #color map
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

def zoom(m, ax, scale=4):
    xc, yc=ndimage.center_of_mass(m.data)
    y0,y1=ax.get_ylim()
    x0,x1=ax.get_xlim()
    if x0<xc<x1:
        ax.set_xlim(  xc- (xc-x0)/scale, xc + (x1 - xc)/scale) 
    if y0<yc<y1:
        ax.set_ylim(  yc- (yc-y0)/scale, yc + (y1 - yc)/scale) 

def plot_flare_image(imap, fig, panel_grid=111, title='', descr='', draw_image=True, contour_levels=[], zoom_scale=1, 
        cmap=CMAP, color='w',   grid_spacing=10*u.deg, 
        text_xy=[0.02,0.95], desc_xy=[0.95, 0.98], vmin=None):
    """
    img_map: sunpy map
    panel_grid: subplot 
    """
    ax = fig.add_subplot(panel_grid, projection=imap)
    if draw_image:
        if vmin is None:
            imap.plot(cmap=cmap, axes=ax, title="")
        else:
            imap.plot(cmap=cmap, axes=ax, title="", vmin=vmin*imap.max())
    imap.draw_grid(color=color, ls='--', grid_spacing=10*u.deg)
    imap.draw_limb(axes=ax, color=color,alpha=0.5)
    if title:
        ax.text(text_xy[0], text_xy[1],title,
             horizontalalignment='left',   verticalalignment='center',    transform = ax.transAxes, color=color)
    if descr:
        ax.text(desc_xy[0], desc_xy[1],descr, 
                horizontalalignment='right', verticalalignment='top',    transform = ax.transAxes, color=color)
    if contour_levels:
        print(contour_levels)
        clevels =np.array(contour_levels)*imap.max()
        cs=imap.draw_contours(clevels)
        plt.clabel(cs, inline=1,  fmt={x: f'{contour_levels[i]*100:.0f} %'    for i, x in enumerate(clevels) })
    if zoom_scale !=1:
        zoom(imap, ax, zoom_scale)    
    ax.set_aspect('equal')
    return  ax
    




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

def create_title_page(pdf, img_id=0, obs='',expt='',sig_id='',bkg_id='',erange='', fig=None, ax=None):
    if fig is None or ax is None:
        fig, ax=plt.subplots( figsize=PDF_FIGURE_SIZE)
    title = 'STIX preview images'

    descr=f'Observation begin UTC:{obs}\nExposure time: {round(expt,2)}\nSignal data unique ID: {sig_id}\nBackground data unique ID: {bkg_id} \nEnergy range:{erange} (keV)'
    remarks=f'Created on {datetime.today().strftime("%B %d, %Y")} \n  by STIX Data Center preview image creator'
    ax.axis("off") 
    ax.text(0.1,0.9,title, transform=fig.transFigure, size=24, ha="left")
    ax.text(0.1,0.8,descr, transform=fig.transFigure, size=14, ha="left", va='top')
    ax.text(0.1,0.5,remarks, transform=fig.transFigure, size=14, ha="left", va='top')
    pdf.savefig()
    plt.close()


def plot_stix_images(doc ):
    """
    Plot STIX images and save them to image files
    Parameters
    """
    if not doc.get('fits',None):
        logger.error('Fits file not found')
        return
    # key should be like "image_xxx"

    ar=u.def_unit("arcsecs",1*u.arcsec)
    am=u.def_unit("meters",1*u.m)
    ad=u.def_unit("degrees",1*u.deg)
    u.add_enabled_units([ar, am, ad])
    #define unit arcsecs is the same arcsec

    img_fname=os.path.join(doc['idl_config']['folder'],f'{doc["idl_config"]["prefix"]}.png')
    pdf_img_fname=os.path.join(doc['idl_config']['folder'],f'{doc["idl_config"]["prefix"]}.pdf')


    maps, map_names=[],[]
    doc_fits=doc['fits']
    find_key=lambda x: [k for k in keys if x in k]

    bsd_id=doc['bsd_id']
    energy_range=doc['energy_range']
    energy_range_str=f'{energy_range[0]} - {energy_range[1]} keV'
    start_utc, end_utc=doc['utc_range']
    bsd_uid=doc['unique_id']
    bkg_uid=doc['background']['unique_id']
    duration=ut.utc2unix(end_utc)-ut.utc2unix(start_utc)
    duration=round(duration,2)


    text_xy=[0.02,0.95]


    
    fig = plt.figure(figsize=(12,7), dpi=100, facecolor='white', constrained_layout=True)

    ax_lc= fig.add_subplot(231)
    lightcurves.plot_QL_lc_for_bsd(bsd_id, fill_between_times=[start_utc, end_utc], ax=ax_lc)




    #Full-disk image


    mbp_full=sunpy.map.Map(doc_fits['image_full_disk'])
    mbp_full= rotate_map(mbp_full, recenter=True)
    #bp
    mbp=sunpy.map.Map(doc_fits['image_bp'])
    mbp= rotate_map(mbp)

    # CLEAN map
    fix_clean_map_fits_header(doc, doc_fits['image_clean'])
    #add dsun to the header
    mclean=sunpy.map.Map(doc_fits['image_clean'])
    mclean=mclean[4]
    mclean= rotate_map(mclean)
    #MEM
    mem=sunpy.map.Map(doc_fits['image_em'])
    mem= rotate_map(mem)
    # FWD-fit
    mfwd=sunpy.map.Map(doc_fits['image_fwdfit'])
    mfwd= rotate_map(mfwd)

    descr_full=f'OBS_BEG: {start_utc} \nExposure time: {duration} (s) \nEnergy range 4 - 10 keV'
    descr=f'OBS_BEG: {start_utc} \nExposure time: {duration} (s) \nEnergy range :{energy_range_str} '

    maps=[mbp_full, mbp, mclean, mem, mfwd]
    try:
        fwdshape=f"({doc['idl_config']['fwdfit_shape']})"
    except (KeyError, TypeError):
        fwdshape=''

    titles=['Back-projection (full disk)',
            'Back-projection',
            'CLEAN',
            'MEM',
            f'VIS_FWDFIT ({fwdshape})']
    panel_ids=[232, 233, 234, 235, 236]
    for i, imap in enumerate(maps):
        vmin=0.4 if i==0 else None
        plot_flare_image(imap, fig, panel_grid=panel_ids[i], title=titles[i], 
                descr='', draw_image=True, contour_levels=[], zoom_scale=1, vmin=vmin)

    image_id_str = f'(#{doc["_id"]})'
    plt.suptitle(f'Start UTC {start_utc}; Exposure time: {duration:.1f} sec; Energy range: {energy_range[0]} – {energy_range[1]} keV ', fontsize=12)
    plt.subplots_adjust(
            top=0.85,
                    wspace=0.2, 
                    hspace=0.4)

    plt.savefig(img_fname, format='png')
    logger.info(f"Images have been written to file:{img_fname}")
    updates={'$set':{'figs':[{'quicklook':img_fname}], 'processing_date':datetime.now()}}
    flare_images_db.update_one({'_id':doc['_id']}, updates)


    ## Print plots to pdf
    logger.info('Creating PDF...')
    with PdfPages(pdf_img_fname) as pdf:


        create_title_page(pdf, doc['_id'], obs=start_utc,expt=duration,sig_id=bsd_uid,bkg_id=bkg_uid,erange=energy_range)


        pfig, (ax_lc_pdf, ax_spec)=plt.subplots(1,2,  figsize=PDF_FIGURE_SIZE)
        lightcurves.plot_QL_lc_for_bsd(bsd_id, fill_between_times=[start_utc, end_utc], ax=ax_lc_pdf)

        l1=ScienceL1.from_fits(doc['filename'])
        selection_box={'trange':[start_utc, end_utc], 'erange':energy_range}
        l1.plot_spectrogram(ax_spec, selection_box)



        pdf.savefig(pfig)
        pfig=plt.figure(figsize=PDF_FIGURE_SIZE)
        vmin=0.3
        plot_flare_image(maps[0], pfig, panel_grid=111, title=titles[0], 
                descr='', draw_image=True, contour_levels=[], zoom_scale=1, vmin=vmin)
        plt.suptitle(descr_full, fontsize=12)
        pdf.savefig(pfig)



        levels=np.array([0.3,  0.5, 0.7, 0.9])
        for i, imap in enumerate(maps):
            if i==0:
                continue
            pfig=plt.figure(figsize=PDF_FIGURE_SIZE, constrained_layout=True)
            plot_flare_image(imap, pfig, panel_grid=131, title=titles[i], 
                    descr='', draw_image=True, contour_levels=[], zoom_scale=1 )

            plot_flare_image(imap, pfig, panel_grid=132, title=titles[i] +' (3x)', 
                    descr='', draw_image=True, contour_levels=[], zoom_scale=3)
            ax=plot_flare_image(imap, pfig, panel_grid=133, title=titles[i], 
                    descr='', draw_image=False, contour_levels=[0.3, 0.5, 0.7], zoom_scale=3, color='k')
            ax.set_xlabel('solar_x [arcsec]')
            ax.set_ylabel('solar_y [arcsec]')


            plt.suptitle(descr, fontsize=12)
            plt.subplots_adjust(
                top=0.95,
                        wspace=0.2, 
                        hspace=1.5)
            pdf.savefig(pfig)
            plt.close()







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
