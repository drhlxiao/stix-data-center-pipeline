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

logger = logger.get_logger()
mdb = db.MongoDB()
flare_image_db = mdb.get_collection('flare_images')

SMALL_SIZE = 9
PDF_FIGURE_SIZE_NORMAL = (16, 7)
PDF_FIGURE_SIZE_SHORT = (16, 5)
PDF_FIGURE_SIZE_LONG = (16, 10)

DEFAULT_PLOT_DPI = 300

matplotlib.rc('font', size=SMALL_SIZE)
matplotlib.rc('axes', titlesize=SMALL_SIZE)
#matplotlib.rcParams['axes.titlepad']=20
#plt.rcParams['figure.constrained_layout.use'] = True

CMAP = 'std_gamma_2'  #color map


def plot_idl(doc_id: int):
    doc = flare_image_db.find_one({'_id': doc_id})
    if not doc:
        logger.error(f'Doc {doc_id} not found in the database')
        return

    try:
        plot_imaging_and_ospex_results(doc)
    except Exception as e:
        #don't raise any exception
        logger.error(e)


def rotate_map(m, recenter=False):
    """ rotate a map
    
    """
    if m.meta['crota2'] == 0:
        #do nothing
        return m
    return m.rotate(angle=(m.meta['crota2']) * u.deg, recenter=recenter)



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
            print("#############", ut.utc2isoformat(hdu.header['DATE_OBS']))
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


def create_pdf_title_page(pdf,
                          img_id=0,
                          asp_source=None,
                          start_utc='',
                          end_utc='',
                          expt='',
                          sig_id='',
                          bkg_id='',
                          erange='',
                          fig=None,
                          ax=None):
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=PDF_FIGURE_SIZE_SHORT)
    title = f'STIX imaging and spectroscopy report #{img_id}'
    descr = f'''Time range:{start_utc} – {end_utc}\nIntegration time: {round(expt,2)}\nSignal file UID: {sig_id}\nBackground file UID: {bkg_id}\nEnergy range: {erange} (keV)'''
    remark_lines = [
        f'Absolute flare location accuracy ~ 1 arcmin'
        if asp_source is None else f'Auxiliary data file: {asp_source}',
        f'\nCreator: STIX Data Center Web Imaging Tool\nCreation time: {datetime.today().strftime("%B %d, %Y")}'
    ]
    remarks = ''.join(remark_lines)
    ax.axis("off")
    ax.text(0.1, 0.9, title, transform=fig.transFigure, size=24, ha="left")
    ax.text(0.1,
            0.8,
            descr,
            transform=fig.transFigure,
            size=14,
            ha="left",
            va='top')
    ax.text(0.1,
            0.5,
            remarks,
            transform=fig.transFigure,
            size=12,
            ha="left",
            va='top')
    pdf.savefig()
    plt.close()


def reverse_colormap(palette_name):
    """reverses matplotlib colormap"""
    if not palette_name.endswith('_r'):
        new_cdata = cm.revcmap(plt.get_cmap(palette_name)._segmentdata)
        new_cmap = matplotlib.colors.LinearSegmentedColormap(
            f'{palette_name}_r', new_cdata)
        return new_cmap
    else:
        return None


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

    img_fname=plot_images(doc, ospex_fig_obj)
    if img_fname is not None:
        figs.append({'images': img_fname})

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

    fname= os.path.join(
        task_doc['idl_config']['folder'],
        f'{task_doc["idl_config"]["prefix"]}_ospex.png')
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

def plot_images(task_doc,  ospex_fig_obj=None, dpi=DEFAULT_PLOT_DPI):
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

    pdf_fname = os.path.join(task_doc['idl_config']['folder'],
                             f'{task_doc["idl_config"]["prefix"]}.pdf')
    img_fname = os.path.join(task_doc['idl_config']['folder'],
                             f'{task_doc["idl_config"]["prefix"]}.png')

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
    lightcurves.plot_QL_lc_for_bsd(bsd_id,
                                   fill_between_times=[start_utc, end_utc],
                                   ax=ax_lc)
    #Full-disk image



    mbp_full = sunpy.map.Map(fits_filename['image_full_disk'])
    mbp_full = rotate_map(mbp_full, recenter=True)
    #bp
    mbp = sunpy.map.Map(fits_filename['image_bp'])
    mbp = rotate_map(mbp)

    if fits_filename.get('image_aia', None) is not None:
        comp_map = sunpy.map.map(mbp,
                                 sunpy.map.Map(fits_filename['image_aia']).submap(
                                     mbp.bottom_left_coord,
                                     top_right=mbp.top_right_coord),
                                 composite=True)
        levels = [70, 80, 90]
        comp_map.set_levels(index=1, levels=levels, percent=True)

        aia_plot_settings = comp_map.get_plot_settings(0)
        aia_plot_settings['cmap'] = reverse_colormap(aia_plot_settings['cmap'])
        comp_map.set_plot_settings(0, aia_plot_settings)
        mbp = comp_map

    # CLEAN map
    #fix_clean_map_fits_header(task_doc, fits_filename['image_clean'])
    #add dsun to the header
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
    fig.subplots_adjust(top=0.85, wspace=0.2, hspace=0.4)

    fig.savefig(img_fname, format='png', dpi=dpi)

    logger.info(f"Images have been written to file:{img_fname}")


    ## Print plots to pdf
    logger.info('Creating PDF...')
    with PdfPages(pdf_fname) as pdf:
        aspect_source_file = task_doc['aux'].get('data_source_file', None)
        create_pdf_title_page(pdf,
                              task_doc['_id'],
                              aspect_source_file,
                              start_utc=start_utc,
                              end_utc=end_utc,
                              expt=duration,
                              sig_id=bsd_uid,
                              bkg_id=bkg_uid,
                              erange=energy_range)
        pfig, (ax_lc_pdf,
               ax_spec) = plt.subplots(1, 2, figsize=PDF_FIGURE_SIZE_SHORT)

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

        plt.subplots_adjust(top=0.95, wspace=0.2, hspace=2)
        pdf.savefig(pfig)


        levels = np.array([0.3, 0.5, 0.7, 0.9])

        for i, imap in enumerate(maps):
            if i == 0:
                continue
            pfig = plt.figure(figsize=PDF_FIGURE_SIZE_NORMAL)
            

            plot_flare_image(imap,
                             pfig,
                             panel_grid=131,
                             title=titles[i],
                             descr='',
                             draw_image=True,
                             contour_levels=[],
                             zoom_ratio=1)
            plot_flare_image(imap,
                             pfig,
                             panel_grid=132,
                             title=titles[i] + ' (2x)',
                             descr='',
                             draw_image=True,
                             contour_levels=[],
                             zoom_ratio=3)
            ax = plot_flare_image(imap,
                                  pfig,
                                  panel_grid=133,
                                  title=titles[i],
                                  descr='',
                                  draw_image=False,
                                  contour_levels=[0.3, 0.5, 0.7],
                                  zoom_ratio=2,
                                  color='k')
            ax.set_xlabel('solar_x [arcsec]')
            ax.set_ylabel('solar_y [arcsec]')
            plt.suptitle(descr, fontsize=10)
            plt.subplots_adjust(top=0.95, wspace=0.2, hspace=2)
            pdf.savefig(pfig)
            plt.close()

        if ospex_fig_obj:
            ospex_fig=ospex_fig_obj['figure']
            ospex_fig.set_figheight(6)
            ospex_fig.set_figwidth(7)
            pdf.savefig(ospex_fig)
        d = pdf.infodict()
        d['Title'] = 'STIX preview images and spectral fitting results'
        d['Author'] = 'STIX web imaging tool v1.0 (contact: hualin.xiao@fhnw.ch)'
        d['CreationDate'] = datetime.now()

    return img_fname



def test():
    create_images_for_ids_between(0, 1)
    create_images_in_queue()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage plot_idl <doc_id>')
    else:
        plot_idl(int(sys.argv[1]))
