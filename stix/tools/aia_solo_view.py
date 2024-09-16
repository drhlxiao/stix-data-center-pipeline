'''
Created by Andrea Francesco Battaglia
andrea-battaglia@ethz.ch

Last modification: 27-May-2021
'''

import sys
############################################################
##### Imports
import astropy.units as u
import glob
import numpy as np
import os
import spiceypy as spice
import sunpy
import warnings
warnings.filterwarnings('ignore')

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from datetime import datetime, timedelta
from reproject import reproject_adaptive, reproject_interp, reproject_exact
from sunpy.coordinates import Helioprojective
from sunpy.coordinates.frames import HeliocentricEarthEcliptic, HeliographicStonyhurst
from sunpy.map import Map, make_fitswcs_header
from sunpy.map.maputils import all_coordinates_from_map
import astropy.units as u
import astropy.io as fits
import matplotlib.pyplot as plt

from aiapy.calibrate import normalize_exposure, register, update_pointing
from aiapy.calibrate import register, update_pointing
from astropy.io import fits
from astropy.time import Time, TimeDelta
from reproject import reproject_exact
from sunpy.coordinates import RotatedSunFrame, transform_with_sun_center
from sunpy.net import attrs as a
from sunpy.net import Fido
from stix.spice import time_utils as sdt

import tempfile
from stix.core import mongo_db as db

mdb = db.MongoDB()

DEFAULT_DOWNLOAD_PATH= '/data/tmp/'
from stix.core import mongo_db, logger
logger = logger.get_logger()


############################################################
##### Main function
def convert_aia_to_solo_view(aia_map,  date_solo=None, center=None, fov=None): 
    """
    Return the input map as seen by Solar Orbiter.

    The map has to be a full disk map. If `center` and `fov` 
    are specified, the returned maps include also the cutout
    of the full disk map, centered on `center` and with the
    dimensions given by `fov`. 

    Parameters
    ------------
    map : `sunpy.map.Map`
        Input map to convert as seen by Solar Orbiter
    path_kernel : string
        String containing the path of the folder in which the
        SPICE kernels are stored (for the moment works only
        with the SPICE kernel I think)
    date : optional `astropy.time.Time` or `datetime.datetime`
        Date used to obtain the Solar Orbiter position. 
        DEFAULT: the time of the input map
    center : optional array_like (float) or array_like (int)
        Center of the region of interest as seen from Earth, 
        given in arcsec. This region of interest will be the
        center of the returned map, but with coordinates as 
        seen by Solar Orbiter.
        DEFAULT: full disk map, without a region of interest
    fov : optional 'float' or optional array_like (float) or 
                    array_like (int)
        Field of view of the region of interest as seen from
        Earth, given in arcmin.
        If `center` is specified, fov has also to be
        specified, otherwise a default value of 5 arcmin in
        both x and y direction is given.
        DEFAULT: full disk map, without specify the fov
    """

    # I don't know if this is a good way of converting
    # string format to datetime.
    if date_solo == None:
        date_solo = datetime.strptime(aia_map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f')

    # If center is specified and not the FOV, set default of 5 arcmin
    if center != None and fov == None:
        fov = [5,5]

    # Set the real radius of the Sun (and not the solar disk)
    aia_map.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value('m')

    # Get the HEE coordinates of Solar Orbiter
    solo_hee = get_solo_coord(date_solo)

    # Mask the off disk data (array type must be float!)
    hpc_coords = all_coordinates_from_map(aia_map)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / aia_map.rsun_obs
    aia_map = check_float(aia_map) # maybe there is a better way of doing this
    aia_map.data[r > 1.0] = np.nan
    
    # Solar Orbiter reference coordinates
    solo_ref_coord = SkyCoord(0*u.arcsec, 
                              0*u.arcsec,
                              obstime=date_solo,
                              observer=solo_hee.transform_to(HeliographicStonyhurst(obstime=date_solo)),
                              frame='helioprojective')

    # Dimension of the output FOV (in pixels)
    out_shape = (3000, 3000) # (4096, 4096)
    dsun_solo = float(np.sqrt(solo_hee.x**2+solo_hee.y**2+solo_hee.z**2)/u.km*1e3)
    dsun_earth = float(aia_map.dsun/u.m)
    factor = dsun_earth / dsun_solo
    pixel_size = [float(aia_map.scale[0]/u.arcsec*u.pixel) * factor,
                      float(aia_map.scale[1]/u.arcsec*u.pixel) * factor]
    # Create a FITS-WCS header from a coordinate object
    out_header = make_fitswcs_header(
                                     out_shape,
                                     solo_ref_coord,
                                     scale=pixel_size*u.arcsec/u.pixel,
                                     instrument="SOLO-AIA",
                                     observatory="SOLO",
                                     wavelength=aia_map.wavelength)

    # Standard World Coordinate System (WCS) transformation
    out_wcs = WCS(out_header)

    # Transform to HGS coordinates
    out_wcs.heliographic_observer = solo_hee.transform_to(HeliographicStonyhurst(obstime=date_solo))

    # Image reprojection
    #output, _ = reproject_adaptive(map, out_wcs, out_shape) # Can give memory problems
    output, _ = reproject_interp(aia_map, out_wcs, out_shape) # The fastest algorithm
    #output, _ = reproject_exact(map, out_wcs, out_shape) # The slowest algorithm

    # 2D map as seen from Solar Orbiter
    solo_map = Map((output, out_header))
    
    # If center == None, then return only the full map, otherwise also the sub-map
    if center == None:
        return solo_map    
    else:
        center_solo_hpc = roi_hpc_SOLO(aia_map, center, date_solo, solo_hee)

        # To convert FOV as seen from Solar Orbiter
        ratio_dist = solo_map.dsun/aia_map.dsun

        # Coordinates of the sub-frame
        bl_solo = SkyCoord(center_solo_hpc.Tx - (fov[0]*60/(ratio_dist*2))*u.arcsec, 
                        center_solo_hpc.Ty - (fov[1]*60/(ratio_dist*2))*u.arcsec, 
                        frame=solo_map.coordinate_frame)
        tr_solo = SkyCoord(center_solo_hpc.Tx + (fov[0]*60/(ratio_dist*2))*u.arcsec, 
                        center_solo_hpc.Ty + (fov[1]*60/(ratio_dist*2))*u.arcsec, 
                        frame=solo_map.coordinate_frame)

        # Extract the sub-map
        solo_submap = solo_map.submap(bottom_left=bl_solo, top_right=tr_solo)

        return solo_map, solo_submap

def get_solo_coord(date_solo):
    """
    Load the kernel needed in order to derive the
    coordinates of Solar Orbiter and then return them in
    Heliocentric Earth Ecliptic (HEE) coordinates.

    ****************************************************
    Maybe this function already exists or maybe it can
    be improved in order to add to STIXCore, since it 
    can be used in many different context I think. 
    ****************************************************
    """



    # Observing time (to get the SOLO coordinates)
    date_solo_dt=Time(date_solo).to_datetime()
    et_solo = spice.datetime2et(date_solo_dt)

    # Obtain the coordinates of Solar Orbiter
    solo_hee_spice, _ = spice.spkpos('SOLO', et_solo, 'SOLO_HEE_NASA', 'NONE', 'SUN')
    solo_hee_spice = solo_hee_spice * u.km

    # Convert the coordinates to HEE
    solo_hee = HeliocentricEarthEcliptic(solo_hee_spice, 
                                         obstime=Time(date_solo).isot, 
                                         representation_type='cartesian')
    
    # Return the HEE coordinates of Solar Orbiter
    return solo_hee

def check_float(imap):
    """
    Check if the data contained in map are float. If it is
    not the case, change the format.
    """

    if imap.data.dtype.kind == 'i':
        imap = Map((imap.data.astype('float'), imap.meta))
    
    return imap

def roi_hpc_SOLO(aia_map, coord, date_solo, solo_hee):
    """
    Takes the coordinates of the region of interest (ROI) as
    seen from Earth (on the surface of the Sun) and
    transform them to hpc coordinates as seen (always on the
    Sun) from Solar Orbiter
    """

    # SkyCoord of the ROI as seen from Earth
    roi_earth_hpc = SkyCoord(coord[0]*u.arcsec, 
                             coord[1]*u.arcsec, 
                             frame=aia_map.coordinate_frame)
    
    # Assume ROI to be on the surface of the Sun
    roi_inter = roi_earth_hpc.transform_to(HeliocentricEarthEcliptic)
    third_dim = 1*u.Rsun

    # ROI location in HEE
    roi_hee = SkyCoord(roi_inter.lon, roi_inter.lat, third_dim, 
                       frame=HeliocentricEarthEcliptic(obstime=date_solo))

    # Since now we have the full 3D coordinate of the ROI position
    # given in HEE, we can now transform that coordinated as seen
    # from Solar Orbiter and give them in Helioprojective coordinates
    roi_solo_hpc = roi_hee.transform_to(Helioprojective(obstime=date_solo, 
                                                        observer=solo_hee.transform_to(HeliographicStonyhurst(obstime=date_solo))))
    
    return roi_solo_hpc
############################################################
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

    out_fname =os.path.join(path,  f"{uid}.fits")
    m_out.save(out_fname)
    return out_fname
def find_cutout_coords(stix_bp_map, earth_time,  padding=10 * u.arcsec):
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


def create_aia_image_as_seen_by_STIX( stix_map,  wavelength=1600):
    '''
    Usually we rotate the 1600 map, since has most likely
    less projection effects compared to the others AIA
    bands
    '''
    earth_time=sdt.utc2datetime(stix_bp_map.meta['date_ear']) + timedelta(seconds=stix_bp_map.meta['exptime']/2.) 
    logger.info("Finding STIX image coord range...")

    bottom_left, top_right=find_cutout_coords(stix_map, earth_time)



    time_start = Time(earth_time)-15*u.s
    time_end = Time(earth_time)+15*u.s
    # Download and calibrate AIA

    #series=24*u.second if wavelength = 1600 else 12 *u.second


        
    logger.info("Downloading AIA...")
    # Search for the data
    query = Fido.search(
                    a.Time(time_start, time_end), a.Instrument.aia,
                    #a.Physobs.intensity,
                    #a.Sample(series),
                    a.Wavelength(wavelength*u.angstrom)

    )
    # Download the data
    logger.info("Downloading FITS files ...")
    result=Fido.fetch(query[0], progress=True, path=DEFAULT_DOWNLOAD_PATH)

    if not result:
        print("AIA 1600 data not available")
        return None
    if isinstance(result, list):
        result=result[0]

    aia_prepped = aia_preprocessing(result, path=DEFAULT_DOWNLOAD_PATH) # calibrate to level 1.5
    aia_map = sunpy.map.Map(aia_prepped)
    # Calibrate the map
    if isinstance(aia_map, list):
        aia_map=aia_map[0]


    aia_map= representation_type
    except Exception as e:
        logger.error("Failed to calibrate AIA")

    if aia_map is None:
        return None
    aia_map_proj= convert_aia_to_solo_view(aia_map)
    aia_submap_proj = aia_map_proj.submap(bottom_left, top_right=top_right)
    return aia_submap





def test():
    create_aia_image_as_seen_by_STIX('/tmp', '2020-07-01 00:00:00')

if __name__=='__main__':
    test()




