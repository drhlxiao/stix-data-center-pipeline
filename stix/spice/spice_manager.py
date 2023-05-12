#!/usr/bin/python3

import re
import os
import sys
import glob
import spiceypy
import numpy as np
from pathlib import Path
from datetime import datetime
from dateutil import parser as dtparser
from astropy.time import Time
from astropy import units as u
from stix.core import logger
from stix.core import config
from stix.core import mongo_db
from spiceypy.utils.exceptions import SpiceBADPARTNUMBER, SpiceINVALIDSCLKSTRING

NUM_KERNEL_FILES_LIMIT = 10

logger = logger.get_logger()

mdb = mongo_db.MongoDB()
# SOLAR ORBITER naif identifier

IGNORE_ATT = True
#loading att files takes a lot time


class _SpiceManager(object):
    """taken from https://issues.cosmos.esa.int/solarorbiterwiki/display/SOSP/Translate+from+OBT+to+UTC+and+back
    """
    __instance = None

    @staticmethod
    def get_instance():
        if not _SpiceManager.__instance:
            _SpiceManager()
        return _SpiceManager.__instance

    #singleton
    #make sure only one instance is created
    def __init__(self):
        if _SpiceManager.__instance:
            logger.warning('SPICE kernel manager already initialized')
            pass
        else:
            _SpiceManager.__instance = self

        self.version_date = datetime.strptime('19700101', "%Y%m%d")

        self.loaded_kernel_filename = None
        self.latest_mk = None
        self.load_kernels()

    def get_kernel_files_from_mk(self, mk_filename):
        fnames = []
        with open(mk_filename) as f:
            content = f.read()
            kernels = re.findall(r"\$KERNELS\/(.*)'", content)
            for k in kernels:
                fnames.append(k)
        return fnames

    def get_kernel_filename(self):
        return self.loaded_kernel_filename

    def load_kernels(self):
        cwd = os.getcwd()
        spice_folder = config.get_config('spice')
        mk_folder = os.path.join(spice_folder, 'mk')
        self.latest_mk = None
        for filename in glob.glob(f'{mk_folder}/solo_ANC_soc-flown-mk*.tm'):
            date_str = re.findall(r"\d{4}\d{2}\d{2}", filename)
            if date_str:
                fdt = datetime.strptime(date_str[0], "%Y%m%d")
                if fdt > self.version_date:
                    self.latest_mk = os.path.basename(filename)
                    self.version_date = fdt

        #if latest_mk!=None and utc<self.version_date:

        if self.loaded_kernel_filename == self.latest_mk or self.latest_mk is None:
            logger.info(
                f'SPICE kernel loaded already: {self.loaded_kernel_filename}.')
            return

        os.chdir(mk_folder)
        #try:
        #    logger.info(f'Loading kernel from mk file:{self.latest_mk}')
        #    spiceypy.furnsh(os.path.join(mk_folder, self.latest_mk))
        #except spiceypy.utils.exceptions.SpiceNOSUCHFILE:
        #    logger.warning('Failed to load mk file, try loading one by one')

        fnames = self.get_kernel_files_from_mk(self.latest_mk)
        num_ignored, num_loaded = 0, 0
        logger.info(f'Loading kernel files...')
        for fname in fnames:
            if 'flown-att' in fname and IGNORE_ATT:
                num_ignored += 1
                continue
            try:
                num_loaded += 1
                spiceypy.furnsh(os.path.join(spice_folder, fname))
            except spiceypy.utils.exceptions.SpiceNOSUCHFILE:
                logger.warning(f'Failed to load: {fname}')
        if num_ignored > 0:
            logger.warning(f'Number of ignored ATT files:{num_ignored}')
        logger.info(f'Number of loaded kernel files : {num_loaded}')
        self.loaded_kernel_filename = self.latest_mk
        os.chdir(cwd)
        #change back to folder

    def obt2utc(self, obt_string):
        """
         Obt to Ephemeris time (seconds past J2000)
         Ephemeris time to Utc
         Format of output epoch: ISOC (ISO Calendar format, UTC)
         Digits of precision in fractional seconds: 3
         """
        try:
            ephemeris_time = spiceypy.scs2e(-144, obt_string)
            return spiceypy.et2utc(ephemeris_time, "ISOC", 3)
        except Exception as e:
            raise e

    def utc2obt(self, utc_string):
        # Utc to Ephemeris time (seconds past J2000)
        ephemeris_time = spiceypy.utc2et(utc_string)
        # Ephemeris time to Obt
        #return ephemeris_time
        obt_string = spiceypy.sce2s(-144, ephemeris_time)
        time_fields = re.search('\/(.*?):(\d*)', obt_string)
        group = time_fields.groups()
        try:
            return int(group[0]) + int(group[1]) / 65536.
        except Exception as e:
            logger.warning(str(e))
            return 0

    def scet2utc(self, coarse, fine=0):
        obt_string = '{}:{}'.format(coarse, fine)
        #print(obt_string)
        return self.obt2utc(obt_string)

    def utc2scet(self, utc):
        # Utc to Ephemeris time (seconds past J2000)
        ephemeris_time = spiceypy.utc2et(utc)
        # Ephemeris time to Obt
        return spiceypy.sce2s(-144, ephemeris_time)

    def get_fits_headers(self, *, start_time, average_time):
        try:
            et = spiceypy.scs2e(-144, str(average_time))
        except (SpiceBADPARTNUMBER, SpiceINVALIDSCLKSTRING):
            et = spiceypy.utc2et(average_time.isot)
        #et=self.utc2et(average_time.isoformat())

        # HeliographicStonyhurst
        solo_sun_hg, sun_solo_lt = spiceypy.spkezr('SOLO', et,
                                                   'SUN_EARTH_CEQU', 'None',
                                                   'Sun')

        # Convert to spherical and add units
        hg_rad, hg_lon, hg_lat = spiceypy.reclat(solo_sun_hg[:3])
        hg_rad = hg_rad * u.km
        hg_lat, hg_lon = (hg_lat * u.rad).to('deg'), (hg_lon * u.rad).to('deg')
        # Calculate radial velocity add units
        rad_vel, *_ = spiceypy.reclat(solo_sun_hg[3:])
        rad_vel = rad_vel * (u.km / u.s)

        rsun_arc = np.arcsin((1 * u.R_sun) / hg_rad).decompose().to('arcsec')

        solo_sun_hee, _ = spiceypy.spkezr('SOLO', et, 'SOLO_HEE', 'None',
                                          'Sun')
        solo_sun_hci, _ = spiceypy.spkezr('SOLO', et, 'SOLO_HCI', 'None',
                                          'Sun')
        solo_sun_hae, _ = spiceypy.spkezr('SOLO', et, 'SUN_ARIES_ECL', 'None',
                                          'Sun')
        solo_sun_heeq, _ = spiceypy.spkezr('SOLO', et, 'SOLO_HEEQ', 'None',
                                           'Sun')
        solo_sun_gse, earth_solo_lt = spiceypy.spkezr('SOLO', et,
                                                      'EARTH_SUN_ECL', 'None',
                                                      'Earth')
        sun_earth_hee, sun_earth_lt = spiceypy.spkezr('Earth', et, 'SOLO_HEE',
                                                      'None', 'Sun')

        precision = 2
        headers = (
            ('SPICE_MK', self.loaded_kernel_filename,
             'SPICE meta kernel file'),
            ('RSUN_ARC', rsun_arc.to_value('arcsec'),
             '[arcsec] Apparent photospheric solar radius'),
            # ('CAR_ROT', ,), Doesn't make sense as we don't have a crpix
            ('HGLT_OBS', np.around(hg_lat.to_value('deg'), precision),
             '[deg] s/c heliographic latitude (B0 angle)'),
            ('HGLN_OBS', np.around(hg_lon.to_value('deg'), precision),
             '[deg] s/c heliographic longitude'),
            # Not mistake same values know by different terms
            ('CRLT_OBS', np.around(hg_lat.to_value('deg'), precision),
             '[deg] s/c Carrington latitude (B0 angle)'),
            ('CRLN_OBS', np.around(hg_lon.to_value('deg'), precision),
             '[deg] s/c Carrington longitude (L0 angle)'),
            ('DSUN_OBS', np.around(hg_rad.to_value('m'),
                                   precision), '[m] s/c distance from Sun'),
            ('HEEX_OBS',
             np.around((solo_sun_hee[0] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Earth Ecliptic X'),
            ('HEEY_OBS',
             np.around((solo_sun_hee[1] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Earth Ecliptic Y'),
            ('HEEZ_OBS',
             np.around((solo_sun_hee[2] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Earth Ecliptic Z'),
            ('HCIX_OBS',
             np.around((solo_sun_hci[0] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Inertial X'),
            ('HCIY_OBS',
             np.around((solo_sun_hci[1] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Inertial Y'),
            ('HCIZ_OBS',
             np.around((solo_sun_hci[2] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Inertial Z'),
            ('HCIX_VOB',
             np.around(
                 (solo_sun_hci[3] * (u.km / u.s)).to_value('m/s'),
                 precision), '[m/s] s/c Heliocentric Inertial X Velocity'),
            ('HCIY_VOB',
             np.around(
                 (solo_sun_hci[4] * (u.km / u.s)).to_value('m/s'),
                 precision), '[m/s] s/c Heliocentric Inertial Y Velocity'),
            ('HCIZ_VOB',
             np.around(
                 (solo_sun_hci[5] * (u.km / u.s)).to_value('m/s'),
                 precision), '[m/s] s/c Heliocentric Inertial Z Velocity'),
            ('HAEX_OBS',
             np.around((solo_sun_hae[0] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Aries Ecliptic X'),
            ('HAEY_OBS',
             np.around((solo_sun_hae[1] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Aries Ecliptic Y'),
            ('HAEZ_OBS',
             np.around((solo_sun_hae[0] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Aries Ecliptic Z'),
            ('HEQX_OBS',
             np.around((solo_sun_heeq[0] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Earth Equatorial X'),
            ('HEQY_OBS',
             np.around((solo_sun_heeq[1] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Earth Equatorial Y'),
            ('HEQZ_OBS',
             np.around((solo_sun_heeq[2] * u.km).to_value('m'),
                       precision), '[m] s/c Heliocentric Earth Equatorial Z'),
            ('GSEX_OBS',
             np.around((solo_sun_gse[0] * u.km).to_value('m'),
                       precision), '[m] s/c Geocentric Solar Ecliptic X'),
            ('GSEY_OBS',
             np.around((solo_sun_gse[1] * u.km).to_value('m'),
                       precision), '[m] s/c Geocentric Solar Ecliptic Y'),
            ('GSEZ_OBS',
             np.around((solo_sun_gse[2] * u.km).to_value('m'),
                       precision), '[m] s/c Geocentric Solar Ecliptic Y'),
            ('OBS_VR', np.around(rad_vel.to_value('m/s'), precision),
             '[m/s] Radial velocity of spacecraft relative to Sun'),
            ('EAR_TDEL', np.around(sun_earth_lt - sun_solo_lt, precision),
             '[s] Time(Sun to Earth) - Time(Sun to S/C)'),
            ('SUN_TIME', np.around(sun_solo_lt,
                                   precision), '[s] Time(Sun to s/c)'),
            ('DATE_EAR', (start_time + np.around(
                (sun_earth_lt - sun_solo_lt), precision) * u.s).fits,
             'Start time of observation, corrected to Earth'),
            ('DATE_SUN',
             (start_time - np.around(sun_solo_lt, precision) * u.s).fits,
             'Start time of observation, corrected to Su'),
        )

        return headers


spice = _SpiceManager()
