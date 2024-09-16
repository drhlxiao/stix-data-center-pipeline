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
import time
import numpy as np
from scipy import ndimage

from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter
from datetime import datetime, timedelta as td
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
#from stix.flare_pipeline import sdo_aia
from stix.flare_pipeline import solo_eui 
from stix.flare_pipeline import imaging_task_manager as itm
from stix.core import mongo_db as db
from stix.core import logger
from stix.utils import bson
from stix.spice import time_utils as ut
from stix.spice import solo 
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



def plot_idl(anydoc, create_stix_images=True, create_eui=True, create_aia=False):
    if isinstance(anydoc, int):
        doc = flare_image_db.find_one({'_id': anydoc})
        doc_id=anydoc
    else:
        doc=anydoc
        doc_id=doc['_id']
    logger.info(f"Processing {doc_id}")
    flare_image_db.update_one({'_id': doc_id}, {'$set': {'py_calls':1, 'plotter':{'start': ut.now()  } }})
    #will not be called again

    if not doc:
        logger.error(f'Doc {doc_id} not found in the database')
        return
    if create_stix_images:
        try:
            plot_imaging_and_ospex_results(doc)
            flare_image_db.update_one({'_id': doc_id}, {'$set': {'plotter':{'end': ut.now()} }})
        except Exception as e:
            flare_image_db.update_one({'_id': doc_id}, {'$set': {'plotter':{'end': ut.now(), 'error':str(e)  } }})
            logger.error(e)
    return


def rotate_map(m, recenter=False):
    """ rotate a map
    
    """
    if isinstance(m, list):
        print("[WARN] MAP  is a list, the last one is taken")
        m=m[-1]


    if m.meta['crota2'] == 0:
        #do nothing
        return m
    return m.rotate(angle=(m.meta['crota2']) * u.deg, recenter=recenter)

def plot_aspect_data(filename, start_unix, end_unix, flare_sun_x=0, flare_sun_y=0, padding=600):
    end_unix=end_unix+padding
    start_unix=start_unix-padding
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
        axes[0].plot(timestamps, sun_x, '--o', label="Measured", alpha=0.5)
        axes[1].plot(timestamps, sun_y,'--o', label="Measured", alpha=0.5)
    else:
        plt.suptitle('Aspect solution not available')

    flare_dt=[ut.unix2datetime(flare_unix_time)]
    axes[0].plot(flare_dt, [flare_sun_x],'+',  label="Applied")
    axes[1].plot(flare_dt, [flare_sun_y],'+',  label="Applied")
    #axes[0].set_xlabel('Time')
    axes[0].set_ylabel('X (arcsec)')
    axes[0].legend()

    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    axes[0].xaxis.set_major_locator(locator)
    axes[0].xaxis.set_major_formatter(formatter)

    #axes[1].set_xlabel('Time')
    axes[1].set_ylabel('Y (arcsec)')
    axes[1].legend()
    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    axes[1].xaxis.set_major_locator(locator)
    axes[1].xaxis.set_major_formatter(formatter)
    axes[0].set_title("X of the Sun Center")
    axes[1].set_title("Y of the Sun Center")

    #fig.suptitle('STIX aspect solutions \nand the solution used for pointing correction')
    #fig.autofmt_xdate()
    fig.tight_layout()
    logger.info(f'Writing aspect solution to :{filename}')
    fig.savefig(filename, dpi=DEFAULT_PLOT_DPI)
    
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
                     grid_spacing=5* u.deg,
                     text_xy=[0.02, 0.95],
                     desc_xy=[0.95, 0.98],
                     vmin=None):
    try:
        _plot_flare_image(imap,
                     fig,
                     panel_grid,
                     title,
                     descr,
                     draw_image,
                     contour_levels,
                     zoom_ratio,
                     cmap,
                     color,
                     grid_spacing,
                     text_xy,
                     desc_xy,
                     vmin)
    except Exception as e:

        print(str(e))

def _plot_flare_image(imap,
                     fig,
                     panel_grid=111,
                     title='',
                     descr='',
                     draw_image=True,
                     contour_levels=[],
                     zoom_ratio=1,
                     cmap=CMAP,
                     color='w',
                     grid_spacing=5* u.deg,
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
    imap.draw_grid(color=color, ls='--', grid_spacing=grid_spacing)
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
        clevels = contour_levels * u.percent
        cs = imap.draw_contours(clevels)
        h,_ = cs.legend_elements()
        labels=[f'{x}%' for x in contour_levels]
        ax.legend(h, labels,  framealpha=0)
    ax.set_aspect('equal')
    return ax


def fix_map_fits_header(image_filename):
    """
       fits headers are not compatible with sunpy we are fixing them here
    """
    try:
        hduls = fits.open(image_filename)
        for hdu in hduls:
            try:
                hdu.header['DATE-OBS']=ut.utc2isoformat(hdu.header['DATE_OBS'])
            except KeyError:
                try:
                    hdu.header['DATE-OBS']=ut.utc2isoformat(hdu.header['DATE_AVG'])
                except KeyError:
                    hdu.header['DATE-OBS']='1970-01-01T00:00:00'
                    

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
    logger.info("Start to create STIX images and spectral fitting plot...")

    try:
        ospex_fig_obj = plot_ospex_results(doc)
    except Exception as e:
        logger.error("Error occurred when creating spectral fitting result..")

    report={}
    if ospex_fig_obj:
        if 'output' in ospex_fig_obj:
            report['ospex']={'filename':ospex_fig_obj['output'], 
                'title':'Spectral fitting result'}

    img_fname, im_report, ephemeris=plot_images(doc, ospex_fig_obj)
    if img_fname is not None:
        figs.append({'images': img_fname})
    if isinstance(im_report, dict):
        report.update(im_report)

    new_values['report']=report
    if ephemeris:
        new_values['flare_ephemeris']=ephemeris



    if ospex_fig_obj is not None:
        new_values['ospex_meta'] = bson.dict_to_json(ospex_fig_obj['meta'])
        figs.append({'spectrum': ospex_fig_obj['output']})

    new_values['figs'] = figs
    updates = {'$set': new_values}
    plt.close('all')
    flare_image_db.update_one({'_id': doc['_id']}, updates)

def plot_ospex_results(task_doc, dpi=DEFAULT_PLOT_DPI):
    """
    plot ospex results
    """

    try:
        fits_filename = task_doc['fits']
    except (KeyError, TypeError):
        logger.error('No imaging or spectral fitting fits file found in the database')
        return None
    if not fits_filename:
        logger.error('No imaging or spectral fitting fits file found in the database')
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

def update_ephemeris_all():
    docs = flare_image_db.find()
    for doc in docs:
        try:
            fits_filename = doc['fits']
            em_file=fits_filename['image_em']
        except Exception as e:
            print('EM file not found')
            continue
        
        mem = sunpy.map.Map(em_file)
        mem = rotate_map(mem)
        
        try:
            flare_ephm=get_flare_ephemeris(mem)
        except Exception as e:
            print(e)
            flare_ephm=None

        if flare_ephm:
            flare_image_db.update_one({'_id': doc['_id']},{'$set':{'flare_ephemeris':flare_ephm}})



def plot_images(task_doc,  ospex_fig_obj=None, dpi=DEFAULT_PLOT_DPI, create_report=True):
    if task_doc.get('signal_data_type', None) != 'PixelData':
        logger.warning('Images can not only created for PixelData!')
        return None, None, None
    try:
        fits_filename = task_doc['fits']
    except (KeyError, TypeError):
        logger.error('No fits information found in the database')
        return None,None, None
    if not fits_filename:
        logger.error('No fits file information in the database')
        return None,None, None
    # key should be like "image_xxx" 
    for _, fname in fits_filename.items():
        fix_map_fits_header(fname)
    out_folder=task_doc['idl_config']['folder']
    fout_prefix=task_doc["idl_config"]["prefix"]
    img_fname = os.path.join(out_folder,
                             f'{fout_prefix}.png')

    fig = plt.figure(figsize=(12, 7), dpi=dpi, facecolor='white')

    #create images

    #maps, map_names = [], []
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


    image_names={
    'image_full_disk': 'Full-disk',
    'image_bp': 'Back-projection',
    'image_clean': 'CLEAN',
    'image_em': 'MEM',
    'image_fwdfit': 'VIS_FWDFIT'
    }  
    try:
        fwdshape = f"({task_doc['idl_config']['fwdfit_shape']})"
    except (KeyError, TypeError):
        fwdshape = ''

    tiles={
            'image_full_disk': 'Back-projection (Full disk)',
            'image_bp': 'Back-projection',  # might want to edit this for composite maps
            'image_clean':'CLEAN',
            'image_em':'MEM GE',
            'image_fwdfit':f'VIS_FWDFIT {fwdshape}'
    }

    valid_maps={}


    for key,val in image_names.items():
        if key =='image_full_disk':
            continue
            #don't make full disk images
        
        recenter=True if key == 'image_full_disk' else False
        if fname:=fits_filename.get(key, None):
            print("Loading map:",val)
            mb= sunpy.map.Map(fname)
            mb= rotate_map(mb, recenter=recenter)
            valid_maps[key]=mb
 
    try:
        flare_ephm=get_flare_ephemeris(valid_maps['image_em'])
    except Exception as e:
        flare_ephm={'error':str(e)}



    # FWD-fit


    aspect_info = 'Aspect corrected' if task_doc['aux'].get(
        'data_source_type', None) == 'Aspect' else 'Aspect not corrected'

    descr_full = f'STIX 4 – 10 keV image ({aspect_info})\n {start_utc} – {end_utc} UT '
    descr= f'STIX {energy_range_str} image ({aspect_info})\n {start_utc} – {end_utc} UT'

    logger.info("Plotting images...")
    panel_ids = [232, 233, 234, 235, 236]
    i=0
    for key, imap in valid_maps.items():
        vmin = 0.4 if i == 0 else None
        logger.info(f"Creating image # {i}")
        title=tiles[key]
        plot_flare_image(imap,
                         fig,
                         panel_grid=panel_ids[i],
                         title=title,
                         descr='',
                         draw_image=True,
                         contour_levels=[],
                         zoom_ratio=1,
                         vmin=vmin)
        i+=1

    image_id_str = f'(#{task_doc["_id"]})'

    try:
        fig.suptitle(descr, fontsize=10)
        fig.subplots_adjust(top=0.85, wspace=0.3, hspace=0.4)
        fig.savefig(img_fname, format='png', dpi=dpi)
    except Exception as e:
        logger.error(f"Failed to create image {str(e)}")
        pass

    logger.info(f"Images have been written to file:{img_fname}")


    logger.info('Creating plots for detailed report ...')
    report = task_doc.get('report',{})
    if create_report:
        #create high resolution plots for analysis report
        pfig, (ax_lc_pdf,
               ax_spec) = plt.subplots(1, 2, figsize=(12, 5))

        #plot light curves and spectrogram
        logger.info('Creating light curves ...')
        lightcurves.plot_QL_lc_for_bsd(bsd_id,
                                       fill_between_times=[start_utc, end_utc],
                                       ax=ax_lc_pdf)

        logger.info('Creating spectrograms ...')

        l1 = ScienceL1.from_fits(task_doc['filename'])
        selection_box = {
            'trange': [start_utc, end_utc],
            'erange': energy_range
        }
        l1.plot_spectrogram(ax_spec, selection_box)

        plt.subplots_adjust(top=0.95, wspace=0.3, hspace=0.2)


        try:
            lc_and_spec_fname = os.path.join(out_folder,
                                 f'{fout_prefix}_lc_and_spec.png')
            plt.savefig(lc_and_spec_fname, dpi=DEFAULT_PLOT_DPI)
            report['A_lc_and_spec']=      {'filename':lc_and_spec_fname, 
                        'title':'Light curves and spectrogram'
                    }
        except ValueError:
            #raised when the fits file only contains one time bin
            pass


        #the prefix is used for sorting 

        levels = np.array([0.3, 0.5, 0.7, 0.9])

        i=0

        for map_key, imap in valid_maps.items():
            #if i == 0:
            #    i+=1
            #    continue
            pfig = plt.figure(figsize=(13,7))
            
            logger.info(f'Creating flare image # {map_key}')
            title=tiles.get(map_key,' ')
            

            plot_flare_image(imap,
                             pfig,
                             panel_grid=121,
                             title=title,
                             descr='',
                             draw_image=True,
                             contour_levels=[],
                             zoom_ratio=1)
            ax = plot_flare_image(imap,
                                  pfig,
                                  panel_grid=122,
                                  title=title,
                                  descr='',
                                  draw_image=False,
                                  contour_levels=[30,40, 50,60,70, 80,90],
                                  zoom_ratio=1,
                                  color='k')
            if ax:
                ax.set_xlabel('solar_x [arcsec]')
                ax.set_ylabel('solar_y [arcsec]')

            pfig.suptitle(descr, fontsize=12)
            pfig.subplots_adjust(top=0.95, wspace=0.3)
            #pfig.tight_layout()
            image_name=image_names.get(map_key,'')
            img_filename = os.path.join(out_folder,
                             f'{fout_prefix}_{image_name}.png')
            try:
                logger.info(img_filename)
                pfig.savefig(img_filename)
                report[f'B_image_{image_name}']={'title':f'{image_name} map',
                        'filename':img_filename} 
            except IndexError:
                logger.error("Index error when saving figures ")
            plt.close()
            i+=1

        asp_filename = os.path.join(out_folder,
                         f'{fout_prefix}_aspect.png')

        logger.info("Creating aspect solution plot ...")
        try:
            plot_aspect_data(asp_filename, task_doc['start_unix'], 
                    task_doc['end_unix'], task_doc['aux']['sun_center'][0], task_doc['aux']['sun_center'][1])
            report['C_aspect']={'filename':asp_filename, 
                    'title':'STIX pointing information', 'subtitle':
                    "x and y of the sun's center in the STIX coord. frame calibrated by the STIX aspect system as well as x and y applied for correction"
                    }
        except Exception as e:
            logger.error(str(e))
        logger.info("done ...")
        
    return img_fname, report, flare_ephm

def plot_aia(_id, wavelen=1600):
    doc=flare_image_db.find_one({'_id':_id})
    try:
        image_clean=doc['fits']['image_clean']
    except (KeyError, TypeError):
        logger.info('Could not create AIA image, failed to find stix clean map !')
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


    aia_rep_map, stix_bp_map, aia_map =sdo_aia.get_projected_aia_map(image_clean, wavelen)
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
        descr= f'AIA {wavelen} reprojected {aia_map.date}  \n STIX {erange} {stix_bp_map.date}'
        fig_aia=sdo_aia.plot_map_reproj(aia_map, aia_rep_map, stix_bp_map, descr=descr)


        aia_fname= os.path.join(out_folder,  f'{fout_prefix}_aia_{wavelen}.png')
        logger.info(f'creating {aia_fname}')
        fig_aia.savefig(aia_fname)
        aia_meta[f'{wavelen}']={'map': aia_fits,  'rep_map': aia_rep_fits}
        report[f'F_aia-{wavelen}']={'filename':aia_fname,
                        'title':f'AIA {wavelen} reprojected image and STIX image'}
    updates={'report':report}
    if aia_meta:
        updates['aia']=aia_meta 

    flare_image_db.update_one({'_id': doc['_id']},{'$set':updates})




def create_stix_for_all(start_id, end_id=8700):
    for _id in range(start_id, end_id):
        try:
            plot_idl(_id, True, True)
        except Exception as e:
            logger.error(e)
            #don't raise any exception

def get_flare_ephemeris(flare_map):
    data = flare_map.data
    peak_index = np.where(data==data.max())
    peak_coord = flare_map.pixel_to_world(peak_index[0]*u.pix, peak_index[1]*u.pix)
    flare_utc=peak_coord.obstime.value
    sun_x, sun_y=peak_coord.Tx.value, peak_coord.Ty.value
    sun_x,sun_y=sun_x[0],sun_y[0]
    meta=solo.SoloEphemeris.get_flare_spice(sun_x,sun_y,flare_utc, observer='earth')

    meta.update({
        'flare_x':sun_x,
        'flare_y': sun_y,
        #'flare_x_units':'arcsec',
        #'flare_y_units':'arcsec',
        })

    if 'observer' in meta:
        meta.pop('observer')
    
    return meta
        
    
    



def create_all_for_all():
    docs = flare_image_db.find({'signal_data_type':'PixelData', 'aia':{'$exists':False}}).sort('_id',-1)
    for doc in docs:
        try:
            logger.info(f'Creating aia images for {doc["_id"]}..')
            plot_idl(doc['_id'], True)
        except Exception as e:
            logger.error(e)
            #don't raise any exception

def process_latest():
    docs = flare_image_db.find( {'py_calls':0 }).sort('_id',-1).limit(1)
    for doc in docs:
        plot_idl(doc['_id'], True)
        
def process_all():
    while True:
        try:
            process_latest()
        except:
            pass
        time.sleep(10)





def test():
    create_images_for_ids_between(0, 1)
    create_images_in_queue()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage plot_idl <doc_id> or plot_idl to start daemon')
        #process_all()
        update_ephemeris_all()
    else:
        plot_idl(int(sys.argv[1]))
