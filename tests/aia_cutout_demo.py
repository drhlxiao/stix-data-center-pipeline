## demonstrate how to use reproject_aia branch of drhlxiao/pystix 

import pystix
from pystix.flare_pipeline import flare_aia
from pystix.imaging import image_viewer
from pystix.spice.solo import SoloEphemeris
from datetime import timedelta as td
import sunpy
import pandas as pd
from sunpy.net import attrs as a

ar=u.def_unit("arcsecs",1*u.arcsec)
u.add_enabled_units([ar])

# Given a STIX back projection map, find the corresponding AIA 1600 cutout to download
bp_file = '/Users/wheatley/Documents/Solar/STIX/tutorials/stix_tutorials/sample_data/stix_2110280107_4-10keV_20211028T152451_bp_cutout.fits'
stix_bp_map = sunpy.map.Map(bp_file)
stix_dt_obs = pd.to_datetime(stix_bp_map.meta['date_obs'])

t = SoloEphemeris() # STIX location
sdo_refcoord = flare_aia.get_AIA_observer(stix_dt_obs) # AIA location
bl,tr = flare_aia.find_cutout_coords(stix_bp_map) # bottom-left and top-right coordinates of STIX back projection map, transformed to the AIA-observer coordinate frame

# download AIA cutout
time_int = [stix_dt_obs, stix_dt_obs + td(seconds = 30)]# time intervaln to use for cutout. AIA UV observations are ~ 24s apart 
result = flare_aia.download_AIA_cutout(time_int, series = a.jsoc.Series.aia_lev1_uv_24s, cutout_coords=(bl, tr), aia_file_path='/Users/wheatley/Documents/Solar/STIX/tutorials/stix_tutorials/sample_data/', jsoc_email = 'erica.lastufka@fhnw.ch')
print(result)
# load AIA cutout
aia_file = '/Users/wheatley/Documents/Solar/STIX/tutorials/stix_tutorials/sample_data/aia.lev1_uv_24s.2021-10-28T152459Z.1600.image.fits' #level 1 data
aia_prepped = flare_aia.aia_prep_py(aia_file) # calibrate to level 1.5
aia_map = sunpy.map.Map(aia_prepped)

# Reproject AIA map
reprojected_aia = image_viewer.reproject_map(aia_map) # a sunpy map
# save as FITS file, put in database, make a request that has this filename in doc_fits['image_aia']
# to be done by Hualin
# call plot_stix_images()
# the reprojected AIA image should be plotted underneath the back projection contours
