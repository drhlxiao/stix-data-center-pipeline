"""
reproject AIA image to STIX coordinate frame

Originally prepared by Erica
Adapted to STIX data center by Hualin Xiao

"""


import os
import sys
sys.path.append('.')
import uuid
import numpy as np
import sunpy
from datetime import datetime
import astropy.units as u
import pandas as pd
from matplotlib import cm
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord

from astropy.wcs import WCS
from astropy.time import Time

from sunpy.coordinates.frames import HeliocentricEarthEcliptic, HeliographicStonyhurst
from astropy import constants as const
from datetime import timedelta 
from sunpy.net import attrs as a
import drms
from sunpy.map import Map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames
from sunpy.coordinates import get_body_heliographic_stonyhurst
from stix.spice.solo import SoloEphemeris

from stix.spice import time_utils as sdt
from stix.spice import spice_manager
from stix.core import logger

ar=u.def_unit("arcsecs",1*u.arcsec)
u.add_enabled_units([ar])


log = logger.get_logger()
from aiapy.calibrate import normalize_exposure, register, update_pointing

AIA_DOWNLOAD_PATH = '/data/tmp/'
def get_AIA_observer(date_obs):
    """Get reference coordinate for SDO/AIA coordinate frame. 
    If no AIA frame information available, use Earth reference frame (can result in differences up to 10" at the limb)"""
    date_obs_str = date_obs.strftime("%Y.%b.%d_%H:%M:%S")
    client = drms.Client()
    kstr = 'HGLN_OBS,HGLT_OBS,DSUN_OBS,RSUN_OBS,DATE-OBS,RSUN_REF,EXPTIME,INSTRUME,WAVELNTH,WAVEUNIT,TELESCOP,LVL_NUM,CROTA2'
    # if no UV data available can also try EUV: qstr=f'aia.lev1_euv_12s[{start_utc}/1m@30s]'

    qstr = f'aia.lev1_uv_24s[{date_obs_str}/1m@30s]'

    try:
        df = client.query(qstr, key=kstr)
        meta = df.iloc[0].to_dict()
        obstime = meta['DATE-OBS']
        sdo_observer = SkyCoord(meta['HGLN_OBS'] * u.deg,
                                meta['HGLT_OBS'] * u.deg,
                                meta['DSUN_OBS'] * u.m,
                                frame=HeliographicStonyhurst,
                                obstime=obstime)
    except Exception:  # timeout - use Earth position as approximation then
        eph = SoloEphemeris()
        positions = eph.get_earth_spice_HEE(
            [sdt.utc2datetime(stix_bp_map.meta['date_obs'])])['positions'][0]
        obstime = stix_bp_map.meta['date_obs']
        sdo_observer = SkyCoord(*positions,
                                frame=HeliocentricEarthEcliptic,
                                obstime=obstime,
                                representation_type='cartesian')
    sdo_refcoord = SkyCoord(
        0 * u.arcsec,
        0 * u.arcsec,
        observer=sdo_observer,
        frame='helioprojective',
        obstime=obstime)  # helioprojective instead of stonyhurst
    return sdo_refcoord

def find_cutout_coords(stix_bp_map, earth_time, padding=10 * u.arcsec):
    """Find the coordinates of AIA cutout that upon reprojection
    corresponds to the area of the STIX back projection map. Returns None if the STIX area is not visible from AIA

    Inputs:
    stix_bp_map: sunpy.map.Map
    (can also write this to use CFL flare position or any other input that contains obstime and central coordinate)
    padding: Quantity
    how much to pad the coordinates

    Outputs:
    bottom_left_coord: SkyCoord or None
    bottom left coord of cutout to request (AIA POV at obstime stix_bp_map.meta['date-obs'])

    top_right_coord: SkyCoord or None
    top right coord of cutout to request (AIA POV at obstime stix_bp_map.meta['date-obs'])

    For usage with:
    cutout = a.jsoc.Cutout(
    bottom_left_coord,
    top_right=top_right_coord,
    tracking=True
    )
    """
    
    
    sdo_refcoord = get_AIA_observer(earth_time)  # warning no light-time correction

    # Get STIX map extent
    stix_bottom_left = stix_bp_map.bottom_left_coord
    stix_top_right = stix_bp_map.top_right_coord

    # STIX coords to AIA coords
    bottom_left_transformed = stix_bottom_left.transform_to(sdo_refcoord)
    top_right_transformed = stix_top_right.transform_to(sdo_refcoord)

    # Check for NaNs and deal with them
    bottom_left_xy, top_right_xy = [], []
    for bl, tr in [(bottom_left_transformed.Tx, top_right_transformed.Tx),
                   (bottom_left_transformed.Ty, top_right_transformed.Ty)]:
        if np.isnan(bl) and np.isnan(tr):
            # All off-disk or not jointly observed.
            return None, None
        if np.isnan(bl): 
            bottom_left_xy.append(-1228 * u.arcsec + padding)
        else:
            bottom_left_xy.append(bl)

        if np.isnan(bl):  # Replace with +-1228"
            top_right_xy.append(1228 * u.arcsec - padding)
        else:
            top_right_xy.append(tr)


    # Pad
    bottom_left_coord = SkyCoord(
        bottom_left_xy[0] + np.sign(bottom_left_xy[0]) * padding,
        bottom_left_xy[1] + np.sign(bottom_left_xy[1]) * padding,
        frame=bottom_left_transformed.frame)
    top_right_coord = SkyCoord(
        top_right_xy[0] + np.sign(top_right_xy[0]) * padding,
        top_right_xy[1] + np.sign(top_right_xy[1]) * padding,
        frame=top_right_transformed.frame)

    return bottom_left_coord, top_right_coord


def download_AIA(time_int,
                        wavelen=1600,
                        cutout_coords=False,
                        single_result=True,
                        aia_file_path=AIA_DOWNLOAD_PATH):
    """Download AIA cutout corresponding to STIX flare."""
    if type(time_int[0]) == str:
        time_int[0] = dt.strptime(time_int[0], '%Y-%m-%dT%H:%M:%S')
        time_int[1] = dt.strptime(time_int[1], '%Y-%m-%dT%H:%M:%S')

    series=a.jsoc.Series('aia.lev1_uv_24s') if wavelen==1600 else a.jsoc.Series('aia.lev1_euv_12s')

    #series=a.Sample(24*u.second) if wavelen==1600 else a.Sample(12*u.second)



    wavelen = a.Wavelength(wavelen * u.angstrom)
    time = a.Time(time_int[0], time_int[1])

    qs = [time, wavelen, series]
    if cutout_coords != False:
        cutout = a.jsoc.Cutout(cutout_coords[0], top_right=cutout_coords[1])

    qs.append(cutout)
    jsoc_email='stix.datacenter@gmail.com'#'test@gmail.com'
    qs.append(a.jsoc.Notify(jsoc_email))  # Set this in config

    res = Fido.search(*qs)
    print("Fido result:", res)
    if single_result:  # Take the first result
        res = res[0]
    files = Fido.fetch(res, path=aia_file_path, overwrite=True)
    if isinstance(files, list):
        files=files[0]

    return files
    """
    aia_map = sunpy.map.Map(files)
    if isinstance(aia_map, list):
        aia_map=aia_map[0]
    #if cutout_coords != False:
    #    aia_map = aia_map.submap(cutout_coords[0], top_right=cutout_coords[1])
    uid=str(uuid.uuid1())
    out_fname =os.path.join(aia_file_path,  f"{uid}.fits")
    aia_map.save(out_fname)
    return out_fname
    """





def aia_preprocessing(fname, expnorm=True, path=None):
    '''Convert from level 1 to level 1.5 data (aia_prep.pro using aiapy instead of IDL)
    from https://aiapy.readthedocs.io/en/latest/generated/gallery/prepping_level_1_data.html'''
    log.info(f"Opening filename:{fname}")
    m = Map(fname)
    try:
        m_registered = register(update_pointing(m))
    except ValueError:  # not full-disk image
        m_registered = m

    if expnorm:
        m_out = normalize_exposure(m_registered)
    else:
        m_out = m_registered
    uid=str(uuid.uuid1())
    out_fname =os.path.join(AIA_DOWNLOAD_PATH,  f"{uid}.fits")
    m_out.save(out_fname)
    return out_fname

def find_reprojected_extent(smap,
                        sobs,
                            swcs=None,
                            full_disk=False,
                            nan_threshold=0.5):
    """Find extent in NAXIS1 and NAXIS2 of reprojected map
    Inputs:
    smap: sunpy.map.Map
    sobs: SkyCoord reference coordinate
    swcs: astropy.wcs
    """
    if full_disk:
        rsun_arcsec = int(
            sunpy.map.solar_angular_radius(sobs).value)  #assuming 1":1px scale
        return -rsun_arcsec, rsun_arcsec, -rsun_arcsec, rsun_arcsec

    if not swcs:  #make the WCS from observer information
        swcs = WCS(sunpy.map.make_fitswcs_header((1, 1), sobs))

    # Obtain the pixel locations of the edges of the reprojected map
    edges_pix = np.concatenate(
        sunpy.map.map_edges(smap)
    )  #what about if off disk? currently doesn't work correctly in that case
    edges_coord = smap.pixel_to_world(edges_pix[:, 0], edges_pix[:, 1])
    new_edges_coord = edges_coord.transform_to(sobs)
    new_edges_xpix, new_edges_ypix = swcs.world_to_pixel(new_edges_coord)

    #check for NaNs...
    if new_edges_xpix[np.isnan(new_edges_xpix)].size > 0 or new_edges_ypix[
            np.isnan(new_edges_ypix)].size > 0:
        cc = sunpy.map.all_coordinates_from_map(smap).transform_to(sobs)
        on_disk = sunpy.map.coordinate_is_on_solar_disk(cc)
        on_disk_coordinates = cc[on_disk]
        nan_percent = 1. - len(on_disk_coordinates) / len(cc.flatten())
        if nan_percent > nan_threshold:
            raise ValueError(
                f"Warning - {nan_percent*100:.1f}% of pixels in reprojected map are NaN!"
            )

    # Determine the extent needed - use of nanmax/nanmin means only on-disk coords are considered
    left, right = np.nanmin(new_edges_xpix), np.nanmax(new_edges_xpix)
    bottom, top = np.nanmin(new_edges_ypix), np.nanmax(new_edges_ypix)
    return left, right, bottom, top



def modify_header(smap, sobs, swcs=None, full_disk=False):
    """Modify the header information to be passed to sunpy.map.Map.reproject_to()
    to contain WCS keywords corresponding to the new location and extent of the reprojected map

    Inputs:
    smap: sunpy.map.Map
    sobs: SkyCoord reference coordinate
    swcs: astropy.wcs
    """
    left, right, bottom, top = find_reprojected_extent(smap,
                                                       sobs,
                                                       swcs,
                                                       full_disk=full_disk)
    # Adjust the CRPIX and NAXIS values
    modified_header = sunpy.map.make_fitswcs_header((1, 1), sobs)
    modified_header['crpix1'] -= left
    modified_header['crpix2'] -= bottom
    modified_header['naxis1'] = int(np.ceil(right - left))
    modified_header['naxis2'] = int(np.ceil(top - bottom))
    return modified_header


def reproject_map(m):
    """Reprojects AIA map into STIX frame. If reprojection does not contain a reasonable amount of non-NaN pixels, return None (no AIA map will be underplotted)
    
    Inputs:
    m : sunpy.map.Map
    AIA (cutout) map from FITS file"""
    # Get reference coordinate
    eph = SoloEphemeris().get_solo_ephemeris(m.meta['date-obs'],
                                             m.meta['date-obs'],
                                             num_steps=1)
    solo_obs = SkyCoord(eph['solo_hee'][0][0] * u.km,
                        eph['solo_hee'][0][1] * u.km,
                        eph['solo_hee'][0][2] * u.km,
                        representation_type='cartesian',
                        frame=HeliocentricEarthEcliptic,
                        obstime=Time(m.meta['date-obs']))
    solo_refcoord = SkyCoord(0 * u.arcsec,
                             0 * u.arcsec,
                             frame='helioprojective',
                             observer=solo_obs,
                             obstime=Time(m.meta['date-obs']))
    # Reproject map
    try:
        reprojected_aia_header = modify_header(m, solo_refcoord)
        reprojected_aia = m.reproject_to(reprojected_aia_header)
        return reprojected_aia
    except Exception as e:
        raise
        log.error(f'Failed to rotate AIA image {str(e)}')
        return None


def reverse_colormap(palette_name):
    """reverses matplotlib colormap"""
    if not palette_name.endswith('_r'):
        new_cdata = cm.revcmap(plt.get_cmap(palette_name)._segmentdata)
        new_cmap = matplotlib.colors.LinearSegmentedColormap(
            f'{palette_name}_r', new_cdata)
        return new_cmap
    else:
        return None


def  get_projected_aia_map(stix_image_fname, wavelen=171):

    stix_bp_map = sunpy.map.Map(stix_image_fname)

    stix_dt_obs = sdt.utc2datetime(stix_bp_map.meta['date_obs'])
    earth_time=sdt.utc2datetime(stix_bp_map.meta['date_ear']) + timedelta(seconds=stix_bp_map.meta['exptime']/2.) 
    #obs start time  at earth

    t = SoloEphemeris() # STIX location
    sdo_refcoord = get_AIA_observer(earth_time) # AIA location
    bl,tr = find_cutout_coords(stix_bp_map, earth_time) 
    # bottom-left and top-right coordinates of STIX back projection map, transformed to the AIA-observer coordinate frame
    # download AIA cutout
    time_int = [earth_time-timedelta(seconds=15), earth_time + timedelta(seconds = 15)]

    log.info('Downloading AIA full image ...')
    filename = download_AIA(time_int, wavelen=wavelen, 
            cutout_coords=(bl, tr), 
            aia_file_path=AIA_DOWNLOAD_PATH)
    log.info('Calibrating cutout...')
    aia_prepped = aia_preprocessing(filename) # calibrate to level 1.5
    if not aia_prepped:
        log.error('Failed to calibrate the image...')
        return None,stix_bp_map, aia_map
    aia_map = sunpy.map.Map(aia_prepped)
    log.info('Re-projecting map...')
    reprojected_aia = reproject_map(aia_map) # a sunpy map
    return reprojected_aia, stix_bp_map,  aia_map
def plot_map_reproj(aia_map, reprojected_map, stix_map, levels = np.array([0.3, 0.5, 0.7, 0.9]), stix_descr=''):
    """
    Plot the original map, reprojected map and observer locations
    Parameters
    ----------
    map : `sunpy.map.Map`
        The input map to be reprojected
    reprojected_map :  `sunpy.map.Map`
        The reprojected map
    Returns
    -------
    `matplotlib.Figure`
        Figure showing the original map, reprojected map, and the observer locations
    """
    color='w'
    

    fig = plt.figure(figsize=(12, 4))
    ax1 = fig.add_subplot(1, 2, 1, projection=aia_map)

    aia_map.plot(axes=ax1)#, title=f'{aia_map.detector} map {aia_map.date}')

    #aia_map.draw_grid(color='w', ls='--', grid_spacing=10 * u.deg)
    #aia_map.draw_limb(axes=ax1, color='w', alpha=0.5)

    ax2 = fig.add_subplot(1, 2, 2, projection=reprojected_map)
    title=f'AIA reprojected {aia_map.date}'
    reprojected_map.plot(axes=ax2, title=title)
    try:
        clevels = np.array(levels) * stix_map.max()
        cs = stix_map.draw_contours(clevels, axes=ax2)
        proxy = [plt.Rectangle((0, 0), 2, 0.1, fc=pc.get_edgecolor()[0])
                for pc in cs.collections]
        legends = [f'{levels[i]*100:.0f} %' for i, x in enumerate(clevels)]
        ax2.legend(proxy, legends, title=stix_descr,frameon=False, bbox_to_anchor=(1.01, 1), loc='upper right', fontsize=6)
    except Exception as e:
        pass
    return fig
def plot_orbit(aia_map, reprojected_map):
    fig = plt.figure(figsize=(6, 4))
    new_observer = reprojected_map.observer_coordinate
    original_observer = aia_map.observer_coordinate
    sun_coords = get_body_heliographic_stonyhurst("sun", aia_map.date)
    ax = fig.add_subplot(1, 1, 1, projection="polar")
    ax.plot(new_observer.lon.to(u.rad),
             new_observer.radius.to(u.AU),
             '+', ms=10,
             label="SolO")
    ax.plot(sun_coords.lon.to(u.rad),
             sun_coords.radius.to(u.AU),
             'o', ms=10,
             label="Sun")
    ax.plot(original_observer.lon.to(u.rad),
             original_observer.radius.to(u.AU),
             'o', ms=7,
                 label="Earth")
    ax.legend(loc='upper right', bbox_to_anchor=(1.30, 1) )
    ax.set_title(f'SolO location  {aia_map.date}')
    return fig

def test():
    wavelen=171
    aia_rep_map, stix_bp_map, aia_map =get_projected_aia_map('/home/xiaohl/Downloads/stix_image.fits',wavelen)
    plot_orbit(aia_map, aia_rep_map)
    plt.show()
    erange=''
    plot_map_reproj(aia_map, aia_rep_map, stix_bp_map, stix_descr=f'STIX map {erange}')
    plt.show()



if __name__=='__main__':
    test()
