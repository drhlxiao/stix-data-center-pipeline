'''
Created by Andrea Francesco Battaglia
andrea-battaglia@ethz.ch

Last modification: 27-May-2021
'''

import sys
sys.path.append('.')
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
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from datetime import datetime
from reproject import reproject_adaptive, reproject_interp, reproject_exact
from sunpy.coordinates import Helioprojective
from sunpy.coordinates.frames import HeliocentricEarthEcliptic, HeliographicStonyhurst
from sunpy.map import Map, make_fitswcs_header
from sunpy.map.maputils import all_coordinates_from_map
import astropy.units as u
import astropy.io as fits
import matplotlib.pyplot as plt

from aiapy.calibrate import register, update_pointing
from astropy.io import fits
from astropy.time import Time, TimeDelta
from reproject import reproject_exact
from sunpy.coordinates import RotatedSunFrame, transform_with_sun_center
from sunpy.net import attrs as a
from sunpy.net import Fido

import tempfile
from stix.spice import spice_manager
from stix.core import mongo_db as db
mdb = db.MongoDB()

############################################################
##### Main function
def download_calibrate_aia(time_start, time_end, wavelength): 
    """
    Download AIA maps found in the interval [time_start, 
    time_end]. Then, the maps are processed to level 1.5
    and the pointing information in the header is updated. 

    Parameters
    ------------
    time_start : `tuple`, `list`, `str`, `pandas.Timestamp`, 
                 `pandas.Series`, `pandas.DatetimeIndex`, 
                 `datetime.datetime`, `datetime.date`, 
                 `numpy.datetime64`, `numpy.ndarray`, 
                 `astropy.time.Time`
        The start time of the range
    time_end : `tuple`, `list`, `str`, `pandas.Timestamp`, 
               `pandas.Series`, `pandas.DatetimeIndex`, 
               `datetime.datetime`, `datetime.date`, 
               `numpy.datetime64`, `numpy.ndarray`, 
               `astropy.time.Time`
        The end time of the range
    wavelength : `int`
        Wavelength of the maps to calibrate
    """
    
    # Search for the data
    query = Fido.search(a.Instrument.aia,
                    a.Physobs.intensity,
                    a.Wavelength(wavelength*u.angstrom),
                    a.Time(time_start, time_end),
    )
    # Download the data
    temp_folder=tempfile.gettempdir()
    result=Fido.fetch(query[0], progress=True, path=temp_folder)
    if not result:
        print("AIA 1600 data not available")
        return None
    aia_map = sunpy.map.Map(result)
    # Calibrate the map
    aia_map = register(update_pointing(aia_map))
    # Return the calibrated map
    return aia_map

############################################################
##### Main function
def as_seen_by_SOLO(map,  date_solo=None, center=None, fov=None): 
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
        date_solo = datetime.strptime(map.meta['date-obs'], '%Y-%m-%dT%H:%M:%S.%f')

    # If center is specified and not the FOV, set default of 5 arcmin
    if center != None and fov == None:
        fov = [5,5]

    # Set the real radius of the Sun (and not the solar disk)
    map.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value('m')

    # Get the HEE coordinates of Solar Orbiter
    solo_hee = coordinates_SOLO(date_solo)

    # Mask the off disk data (array type must be float!)
    hpc_coords = all_coordinates_from_map(map)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / map.rsun_obs
    map = check_float(map) # maybe there is a better way of doing this
    map.data[r > 1.0] = np.nan
    
    # Solar Orbiter reference coordinates
    solo_ref_coord = SkyCoord(0*u.arcsec, 
                              0*u.arcsec,
                              obstime=date_solo,
                              observer=solo_hee.transform_to(HeliographicStonyhurst(obstime=date_solo)),
                              frame='helioprojective')

    # Dimension of the output FOV (in pixels)
    out_shape = (3000, 3000) # (4096, 4096)

    # Create a FITS-WCS header from a coordinate object
    out_header = make_fitswcs_header(
                                     out_shape,
                                     solo_ref_coord,
                                     scale=(1.2, 1.2)*u.arcsec/u.pixel,
                                     instrument="SOLO-AIA",
                                     observatory="SOLO",
                                     wavelength=map.wavelength)

    # Standard World Coordinate System (WCS) transformation
    out_wcs = WCS(out_header)

    # Transform to HGS coordinates
    out_wcs.heliographic_observer = solo_hee.transform_to(HeliographicStonyhurst(obstime=date_solo))

    # Image reprojection
    #output, _ = reproject_adaptive(map, out_wcs, out_shape) # Can give memory problems
    output, _ = reproject_interp(map, out_wcs, out_shape) # The fastest algorithm
    #output, _ = reproject_exact(map, out_wcs, out_shape) # The slowest algorithm

    # 2D map as seen from Solar Orbiter
    solo_map = Map((output, out_header))
    
    # If center == None, then return only the full map, otherwise also the sub-map
    if center == None:
        return solo_map    
    else:
        center_solo_hpc = roi_hpc_SOLO(map, center, date_solo, solo_hee)

        # To convert FOV as seen from Solar Orbiter
        ratio_dist = solo_map.dsun/map.dsun

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
############################################################


############################################################
##### Support functions
def coordinates_SOLO(date_solo):
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
    et_solo = spice.datetime2et(date_solo)

    # Obtain the coordinates of Solar Orbiter
    solo_hee_spice, _ = spice.spkpos('SOLO', et_solo, 'SOLO_HEE_NASA', 'NONE', 'SUN')
    solo_hee_spice = solo_hee_spice * u.km

    # Convert the coordinates to HEE
    solo_hee = HeliocentricEarthEcliptic(solo_hee_spice, 
                                         obstime=Time(date_solo).isot, 
                                         representation_type='cartesian')
    
    # Return the HEE coordinates of Solar Orbiter
    return solo_hee

def check_float(map):
    """
    Check if the data contained in map are float. If it is
    not the case, change the format.
    """

    if map.data.dtype.kind == 'i':
        map = Map((map.data.astype('float'), map.meta))
    
    return map

def roi_hpc_SOLO(map, coord, date_solo, solo_hee):
    """
    Takes the coordinates of the region of interest (ROI) as
    seen from Earth (on the surface of the Sun) and
    transform them to hpc coordinates as seen (always on the
    Sun) from Solar Orbiter
    """

    # SkyCoord of the ROI as seen from Earth
    roi_earth_hpc = SkyCoord(coord[0]*u.arcsec, 
                             coord[1]*u.arcsec, 
                             frame=map.coordinate_frame)
    
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
def plot(folder,_id, flare_id ,time_at_peak, wavelen=1600,  overwrite=False):
    '''
    Usually we rotate the 1600 map, since has most likely
    less projection effects compared to the others AIA
    bands
    '''
    wlth = wavelen

    '''
    Since I don't know how you deal with the SPICE kernel,
    here I put the variable in which the path to the MK folder
    is needed. One could also modify this directly in
    SolO.as_seen_by_SOLO
    '''
    key='aia1600'
    if  mdb.get_flare_joint_obs(_id, key) and overwrite == False:
        print(f'GOES LC for Flare {flare_id} is not created!')
        return 

    time_start = Time(time_at_peak)-6*u.s
    time_end = Time(time_at_peak)+6*u.s
    # Download and calibrate AIA
    aia_map = download_calibrate_aia(time_start, time_end, wlth)
    if aia_map is None:
        return ''
    solo_map = as_seen_by_SOLO(aia_map)
    str_title = 'AIA '+str(solo_map.meta['wavelnth'])+             r' $\AA$'+' as seen by Solar Orbiter \n  '+solo_map.meta['date-obs']
    plt.figure(figsize=(8, 7), dpi=100, facecolor='white')
    fig = plt.figure(figsize=(7, 7), dpi=100, facecolor='white')
    ax = fig.add_subplot(111, projection=solo_map)
    plt.suptitle(str_title)
    solo_map.plot(vmin=0, vmax=np.nanmax(solo_map.data)/10., 
                  cmap='sdoaia'+str(np.int(solo_map.meta['wavelnth'])), 
                  axes=ax, title=False)
    #plt.colorbar()
    solo_map.draw_grid(color='k', ls='--')
    fname=os.path.join(folder, f'AIA_{wlth}_{_id}_{flare_id}.png')
    #lon, lat = ax.coords
    plt.savefig(fname, dpi=100)
    mdb.update_flare_joint_obs(_id, key, [fname])
    #plt.show()
    return fname





