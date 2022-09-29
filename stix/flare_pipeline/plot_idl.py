#!/usr/bin/python
"""
    filename: plot_idl.py

        Plot imaging and ospex spectral fitting results
    History:
        April 27, 2022, first version
        Aug 17, 2022, included code to plot ospex spectral fitting result
    Author(s): Hualin Xiao (hualin.xiao@fhnw.ch)

"""
import os
import sys
import matplotlib
import numpy as np
from scipy import ndimage

from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates

from datetime import datetime, timedelta as td

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.io import fits
from dateutil.parser import parse as dtparse

import sunpy
from sunpy.map import make_fitswcs_header
from sunpy.coordinates.frames import HeliocentricEarthEcliptic, HeliographicStonyhurst

from stix.flare_pipeline import lightcurves
from stix.flare_pipeline import ospex
from stix.flare_pipeline import sdo_aia
from stix.flare_pipeline import imaging_task_manager as itm
from stix.core import mongo_db as db
from stix.core import logger
from stix.utils import bson
from stix.spice import time_utils as ut
from stix.spice.solo import SoloEphemeris
from stix.analysis.science_l1 import ScienceL1

ar = u.def_unit("arcsecs", 1 * u.arcsec)
am = u.def_unit("meters", 1 * u.m)
ad = u.def_unit("degrees", 1 * u.deg)
u.add_enabled_units([ar, am, ad])
#define unit arcsecs is the same arcsec
matplotlib.use('Agg')
logger = logger.get_logger()
mdb = db.MongoDB()
flare_image_db = mdb.get_collection('flare_images')

SMALL_SIZE = 9

DEFAULT_PLOT_DPI = 300

matplotlib.rc('font', size=SMALL_SIZE)
matplotlib.rc('axes', titlesize=SMALL_SIZE)
#matplotlib.rcParams['axes.titlepad']=20
#plt.rcParams['figure.constrained_layout.use'] = True

CMAP = 'std_gamma_2'  #color map



def plot_idl(doc_id: int, create_aia=False):
    doc = flare_image_db.find_one({'_id': doc_id})
    if not doc:
        logger.error(f'Doc {doc_id} not found in the database')
        return

    try:
        plot_imaging_and_ospex_results(doc)
    except Exception as e:
        #don't raise any exception
        logger.error(e)
    if not create_aia:
        return
    doc = flare_image_db.find_one({'_id': doc_id})
    try:
        logger.info(f'Creating aia images for {doc["_id"]}..')
        plot_aia(doc)
    except Exception as e:
        raise
        logger.error(e)
        #don't raise any exception


def rotate_map(m, recenter=False):
    """ rotate a map
    
    """
    if m.meta['crota2'] == 0:
        #do nothing
        return m
    return m.rotate(angle=(m.meta['crota2']) * u.deg, recenter=recenter)

def plot_aspect_data(filename, start_unix, end_unix, flare_sun_x=0, flare_sun_y=0):
    asp_docs=mdb.get_aspect_solutions(start_unix, end_unix)
    timestamps=[]
    fig, axes=plt.subplots(1,2, figsize=(7,3))
    flare_unix_time=(start_unix+end_unix)/2.

    sun_x=[]
    sun_y=[]
    for doc in asp_docs:
        timestamps.append(ut.unix2datetime(doc['unix_time']))
        x, y = itm.get_sun_center_from_aspect(doc['y_srf'], doc['z_srf'])
        sun_x.append(x)
        sun_y.append(y)
    if timestamps:
        axes[0].plot(timestamps, sun_x, label="sun_x")
        axes[1].plot(timestamps, sun_y, label="sun_y")
    else:
        plt.suptitle('Aspect solution not available')

    flare_dt=[ut.unix2datetime(flare_unix_time)]
    axes[0].plot(flare_dt, [flare_sun_x],'x',  label="offset_x used")
    axes[1].plot(flare_dt, [flare_sun_y],'x',  label="offset_y used")
    axes[0].set_xlabel('Time')
    axes[0].set_ylabel('X (arcsec)')
    axes[0].legend()

    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    axes[0].xaxis.set_major_locator(locator)
    axes[0].xaxis.set_major_formatter(formatter)

    axes[1].set_xlabel('Time')
    axes[1].set_ylabel('Y (arcsec)')
    axes[1].legend()
    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    axes[1].xaxis.set_major_locator(locator)
    axes[1].xaxis.set_major_formatter(formatter)

    fig.suptitle('STIX aspect solutions')
    #fig.autofmt_xdate()
    fig.tight_layout()
    logger.info(f'Writing aspect solution to :{filename}')
    fig.savefig(filename, dpi=DEFAULT_PLOT_DPI)
    



def zoom(m, ax, ratio=2):
    """Zoom in/out image
    Parameters:
    -----------
    m: sunpy.map.Map
        sunpy map
    ax: matplotlib axes
    ratio:  float
        zoom ratio
    Returns
    -----------
    None
    """
    max_coord= np.where(np.max(m.data)==m.data)
    xc,yc=max_coord[0][0], max_coord[1][0]
    y0, y1 = ax.get_ylim()
    x0, x1 = ax.get_xlim()
    if x0 < xc < x1:
        ax.set_xlim(xc - (xc - x0) / ratio, xc + (x1 - xc) / ratio)
    if y0 < yc < y1:
        ax.set_ylim(yc - (yc - y0) / ratio, yc + (y1 - yc) / ratio)


def plot_flare_image(imap,
                     fig,
                     panel_grid=111,
                     title='',
                     descr='',
                     draw_image=True,
                     contour_levels=[],
                     zoom_ratio=1,
                     cmap=CMAP,
                     color='w',
                     grid_spacing=10 * u.deg,
                     text_xy=[0.02, 0.95],
                     desc_xy=[0.95, 0.98],
                     vmin=None):
    """
    img_map: sunpy map
    panel_grid: subplot 
    """
    ax = fig.add_subplot(panel_grid, projection=imap)
    if draw_image:
        if vmin is None:
            imap.plot(cmap=cmap, axes=ax, title="")
        else:
            imap.plot(cmap=cmap, axes=ax, title="", vmin=vmin * imap.max())
    imap.draw_grid(color=color, ls='--', grid_spacing=10 * u.deg)
    imap.draw_limb(axes=ax, color=color, alpha=0.5)
    if title:
        ax.text(text_xy[0],
                text_xy[1],
                title,
                horizontalalignment='left',
                verticalalignment='center',
                transform=ax.transAxes,
                color=color)
    if descr:
        ax.text(desc_xy[0],
                desc_xy[1],
                descr,
                horizontalalignment='right',
                verticalalignment='top',
                transform=ax.transAxes,
                color=color)
    if contour_levels:
        #print(contour_levels)
        clevels = np.array(contour_levels) * imap.max()
        cs = imap.draw_contours(clevels)
        proxy = [
            plt.Rectangle((1, 1), 2, 1, fc=pc.get_edgecolor()[0])
            for pc in cs.collections
        ]
        legends = [
            f'{contour_levels[i]*100:.0f} %' for i, x in enumerate(clevels)
        ]
        plt.legend(proxy, legends)
    if zoom_ratio != 1:
        zoom(imap, ax, zoom_ratio)
    ax.set_aspect('equal')
    return ax


def fix_map_fits_header(image_filename):
    """
       fits headers are not compatible with sunpy we are fixing them here
    """
    try:
        hduls = fits.open(image_filename)
        for hdu in hduls:
            hdu.header['DATE-OBS']=ut.utc2isoformat(hdu.header['DATE_OBS'])
        hduls.writeto(image_filename, overwrite=True, checksum=True)
    except Exception as e:
        logger.error(str(e))
    logger.info(f"Adding more keywords to {image_filename}")


def create_images_in_queue(num=10000):
    cursor = flare_image_db.find({
        'fits': {
            '$gt': {}
        },
        'figs.0': {
            '$exists': False
        }
    }).sort('_id', -1).limit(num)
    for doc in cursor:
        logger.info("Creating images for doc: {doc['_id']}...")
        try:
            plot_imaging_and_ospex_results(doc)
        except Exception as e:
            #don't raise any exception
            logger.error(e)


def create_images_for_ids_between(start_id, end_id):
    for i in range(start_id, end_id):
        doc = flare_image_db.find_one({'_id': i})
        if not doc:
            logger.warning(f"Failed to create figures for DocID:{i}")
            continue
        logger.info(f"Creating images for BSD#{doc['bsd_id']}, DocID:{i}")
        plot_imaging_and_ospex_results(doc)


def plot_imaging_and_ospex_results(doc):
    """
    Plot STIX images and ospex result and save them to pdf and pdf files 
    
    """
    new_values = {'processing_date': datetime.now()}
    figs=[]

    try:
        ospex_fig_obj = plot_ospex_results(doc)
    except Exception as e:
        logger.error("Error occurred when creating spectral fitting result..")

    img_fname, report=plot_images(doc, ospex_fig_obj)
    if img_fname is not None:
        figs.append({'images': img_fname})

    if report is not None:
        if ospex_fig_obj:
            if 'output' in ospex_fig_obj:
                report['ospex']={'filename':ospex_fig_obj['output'], 
                    'title':'Spectral fitting result'}
        new_values['report']=report



    if ospex_fig_obj is not None:
        new_values['ospex_meta'] = bson.dict_to_json(ospex_fig_obj['meta'])
        figs.append({'spectrum': ospex_fig_obj['output']})

    new_values['figs'] = figs
    updates = {'$set': new_values}
    flare_image_db.update_one({'_id': doc['_id']}, updates)

def plot_ospex_results(task_doc, dpi=DEFAULT_PLOT_DPI):
    """
    plot ospex results
    """

    try:
        fits_filename = task_doc['fits']
    except (KeyError, TypeError):
        logger.error('No fits information found in the database')
        return None
    if not fits_filename:
        logger.error('No fits file information in the database')
        return None
    # key should be like "image_xxx"

    ospex_fig_obj = None
    ospex_fname = fits_filename.get('ospex_results', None)
    #ospex result fits filename 
    out_folder=task_doc['idl_config']['folder']
    fout_prefix=task_doc["idl_config"]["prefix"]
    fname= os.path.join(
        out_folder,
        f'{fout_prefix}_ospex.png')
    #png filename

    if ospex_fname:
        logger.info("Creating OSPEX spectral fitting result...")
        ospex_fig = plt.figure(figsize=(6, 7), dpi=dpi, facecolor='white')
        try:
            ospex_fig_obj = ospex.plot_ospex(ospex_fname, ospex_fig)
            ospex_fig.savefig(fname, dpi=dpi)
            logger.info(f'Writing ospex results to {fname}')
        except Exception as e:
            ospex_fig_obj = None
            #don't raise any exception
            logger.warning(str(e))
            logger.warning('Failed to create spectral fitting result plot')
        ospex_fig_obj['output']=fname
        ospex_fig_obj['figure']=ospex_fig

    return ospex_fig_obj



def plot_images(task_doc,  ospex_fig_obj=None, dpi=DEFAULT_PLOT_DPI, create_report=True):
    if task_doc.get('signal_data_type', None) != 'PixelData':
        logger.warning('Images can not only created for PixelData!')
        return None
    try:
        fits_filename = task_doc['fits']
    except (KeyError, TypeError):
        logger.error('No fits information found in the database')
        return None
    if not fits_filename:
        logger.error('No fits file information in the database')
        return None
    # key should be like "image_xxx" 
    for _, fname in fits_filename.items():
        fix_map_fits_header(fname)
    out_folder=task_doc['idl_config']['folder']
    fout_prefix=task_doc["idl_config"]["prefix"]
    img_fname = os.path.join(out_folder,
                             f'{fout_prefix}.png')

    fig = plt.figure(figsize=(12, 7), dpi=dpi, facecolor='white')

    #create images

    maps, map_names = [], []
    find_key = lambda x: [k for k in keys if x in k]

    bsd_id = task_doc['bsd_id']


    energy_range = task_doc['energy_range']
    energy_range_str = f'{energy_range[0]} - {energy_range[1]} keV'
    start_utc, end_utc = task_doc['utc_range']
    bsd_uid = task_doc['unique_id']
    bkg_uid = task_doc['background']['unique_id']
    duration = ut.utc2unix(end_utc) - ut.utc2unix(start_utc)
    duration = round(duration, 2)

    text_xy = [0.02, 0.95]

    ax_lc = fig.add_subplot(231)
    try: 
        lightcurves.plot_QL_lc_for_bsd(bsd_id,
                                   fill_between_times=[start_utc, end_utc],
                                   ax=ax_lc)
    except Exception as e:
        logger.error(str(e))
    #Full-disk image


    logger.info("Creating full disk images...")

    mbp_full = sunpy.map.Map(fits_filename['image_full_disk'])
    mbp_full = rotate_map(mbp_full, recenter=True)
    #bp
    logger.info("Creating BP images...")
    mbp = sunpy.map.Map(fits_filename['image_bp'])
    mbp = rotate_map(mbp)

    mclean = sunpy.map.Map(fits_filename['image_clean'])
    #mclean=mclean[4]
    mclean = rotate_map(mclean)
    #MEM
    mem = sunpy.map.Map(fits_filename['image_em'])
    mem = rotate_map(mem)
    # FWD-fit
    mfwd = sunpy.map.Map(fits_filename['image_fwdfit'])
    mfwd = rotate_map(mfwd)
    aspect_info = 'Aspect corrected' if task_doc['aux'].get(
        'data_source_type', None) == 'Aspect' else 'Aspect not corrected'

    descr_full = f'{start_utc} – {end_utc} \n4 - 10 keV  {aspect_info}'
    descr = f'{start_utc} – {end_utc}\n{energy_range_str} {aspect_info} '

    maps = [mbp_full, mbp, mclean, mem, mfwd]
    map_names=['Full-disk', 'Back-projection', 'CLEAN', 'MEM', 'VIS_FWDFIT']
    try:
        fwdshape = f"({task_doc['idl_config']['fwdfit_shape']})"
    except (KeyError, TypeError):
        fwdshape = ''

    titles = [
        'Back-projection (Full disk)',
        'Back-projection',  # might want to edit this for composite maps
        'CLEAN',
        'MEM',
        f'VIS_FWDFIT {fwdshape}'
    ]
    panel_ids = [232, 233, 234, 235, 236]
    for i, imap in enumerate(maps):
        vmin = 0.4 if i == 0 else None
        plot_flare_image(imap,
                         fig,
                         panel_grid=panel_ids[i],
                         title=titles[i],
                         descr='',
                         draw_image=True,
                         contour_levels=[],
                         zoom_ratio=1,
                         vmin=vmin)

    image_id_str = f'(#{task_doc["_id"]})'

    fig.suptitle(descr, fontsize=10)
    fig.subplots_adjust(top=0.85, wspace=0.4, hspace=0.4)

    fig.savefig(img_fname, format='png', dpi=dpi)

    logger.info(f"Images have been written to file:{img_fname}")


    logger.info('Creating plots for detailed report ...')
    report = task_doc.get('report',{})
    if create_report:
        #create high resolution plots for analysis report
        pfig, (ax_lc_pdf,
               ax_spec) = plt.subplots(1, 2, figsize=(12, 5))

        #plot light curves and spectrogram
        lightcurves.plot_QL_lc_for_bsd(bsd_id,
                                       fill_between_times=[start_utc, end_utc],
                                       ax=ax_lc_pdf)

        l1 = ScienceL1.from_fits(task_doc['filename'])
        selection_box = {
            'trange': [start_utc, end_utc],
            'erange': energy_range
        }
        l1.plot_spectrogram(ax_spec, selection_box)

        #plt.subplots_adjust(top=0.95, wspace=0.2, hspace=2)
        pfig.tight_layout()

        lc_and_spec_fname = os.path.join(out_folder,
                                 f'{fout_prefix}_lc_and_spec.png')

        plt.savefig(lc_and_spec_fname, dpi=DEFAULT_PLOT_DPI)
        report['A_lc_and_spec']=      {'filename':lc_and_spec_fname, 
                    'title':'Light curves and spectrogram'
                    }
        #the prefix is used for sorting 

        levels = np.array([0.3, 0.5, 0.7, 0.9])

        for i, imap in enumerate(maps):
            if i == 0:
                continue
            pfig = plt.figure(figsize=(13,7))
            

            plot_flare_image(imap,
                             pfig,
                             panel_grid=121,
                             title=titles[i],
                             descr='',
                             draw_image=True,
                             contour_levels=[],
                             zoom_ratio=1)
            ax = plot_flare_image(imap,
                                  pfig,
                                  panel_grid=122,
                                  title=titles[i],
                                  descr='',
                                  draw_image=False,
                                  contour_levels=[0.3, 0.5, 0.7],
                                  zoom_ratio=1,
                                  color='k')
            ax.set_xlabel('solar_x [arcsec]')
            ax.set_ylabel('solar_y [arcsec]')
            pfig.suptitle(descr, fontsize=10)
            pfig.subplots_adjust(top=0.95, wspace=0.2, hspace=4)
            img_filename = os.path.join(out_folder,
                             f'{fout_prefix}_{map_names[i]}.png')
            try:
                logger.info(img_filename)
                pfig.savefig(img_filename)
                report[f'B_image_{map_names[i]}']={'title':f'{map_names[i]} map',
                        'filename':img_filename} 
            except IndexError:
                logger.error("Index error when saving figures ")
            plt.close()

        asp_filename = os.path.join(out_folder,
                         f'{fout_prefix}_aspect.png')

        logger.info("Creating aspect solution plot ...")
        try:
            plot_aspect_data(asp_filename, task_doc['start_unix'], task_doc['end_unix'], task_doc['aux']['sun_center'][0], task_doc['aux']['sun_center'][1])
            report['C_aspect']={'filename':asp_filename, 
                    'title':'STIX pointing information'}
        except Exception as e:
            logger.error(str(e))
    return img_fname, report

def plot_aia(doc, wavelen=1600):
    try:
        image_em=doc['fits']['image_em']
    except (KeyError, TypeError):
        logger.info('Could not create AIA image, can not read info from STIX EM image !')
        return
   
    def exists(objs, ts):
        for o in objs:
            if o['type'] == ts:
                return True
        return False

    report=doc.get('report',{})
    aia_meta=doc.get('aia',{})
    out_folder=doc['idl_config']['folder']
    fout_prefix=doc["idl_config"]["prefix"]


    aia_rep_map, stix_bp_map, aia_map =sdo_aia.get_projected_aia_map(image_em, wavelen)
    if aia_map is None:
        logger.info('AIA image is None!')
        return


    orbit_fname= os.path.join(out_folder,  f'{fout_prefix}_orbit.png')
    fig_orbit=sdo_aia.plot_orbit(aia_map, stix_bp_map)
    logger.info(f'creating {orbit_fname}')
    fig_orbit.savefig(orbit_fname)
    report['D_location']={'filename':orbit_fname,
                        'title':'SolO location'}


    if aia_rep_map is None or aia_map is None:
        logger.info('Could not create AIA image, AIA map is none!')
    else:
        aia_fits= os.path.join(out_folder,  f'{fout_prefix}_aia_{wavelen}.fits')
        aia_rep_fits= os.path.join(out_folder,  f'{fout_prefix}_aia_{wavelen}_reprojected.fits')

        aia_map.save(aia_fits, overwrite=True)
        aia_rep_map.save(aia_rep_fits, overwrite=True)

        erange=f"{doc['energy_range'][0]} – {doc['energy_range'][1]} keV"
        fig_aia=sdo_aia.plot_map_reproj(aia_map, aia_rep_map, stix_bp_map, stix_descr=f'STIX EM {erange} ')


        aia_fname= os.path.join(out_folder,  f'{fout_prefix}_aia_{wavelen}.png')
        logger.info(f'creating {aia_fname}')
        fig_aia.savefig(aia_fname)
        aia_meta[f'{wavelen}']={'map': aia_fits,  'rep_map': aia_rep_fits}
        report[f'F_aia-{wavelen}']={'filename':aia_fname,
                        'title':f'AIA {wavelen} image and STIX image'}
    updates={'report':report}
    if aia_meta:
        updates['aia']=aia_meta 

    flare_image_db.update_one({'_id': doc['_id']},{'$set':updates})




    

def create_all_for_all():
    docs = flare_image_db.find({'signal_data_type':'PixelData', 'aia':{'$exists':False}}).sort('_id',-1)
    for doc in docs:
        try:
            logger.info(f'Creating aia images for {doc["_id"]}..')
            plot_idl(doc['_id'], True)
        except Exception as e:
            logger.error(e)
            #don't raise any exception





def test():
    create_images_for_ids_between(0, 1)
    create_images_in_queue()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage plot_idl <doc_id>')
    else:
        plot_idl(int(sys.argv[1]))
