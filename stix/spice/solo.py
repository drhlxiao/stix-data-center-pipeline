import sys
sys.path.append('.')
import os
import math
import numpy as np
from datetime import datetime, timedelta
import spiceypy as sp
import astropy.units as u
from astropy import constants as const

from spiceypy.utils import support_types as spiceytypes
import astropy.time as time
from scipy.spatial.transform import Rotation 
import astropy.coordinates as astrocoords
import sunpy
import sunpy.coordinates as suncoords
from scipy import linalg
from astropy.coordinates import SkyCoord
from astropy.time import Time


from sunpy.coordinates.frames import HeliocentricEarthEcliptic, HeliographicStonyhurst

from stix.spice import spice_manager
from stix.core import logger
from stix.utils import bson
from stix.spice import time_utils as stix_datetime

logger = logger.get_logger()

MAX_STEPS=10000

# Mapping from SPICE frame name to (frame, frame kwargs)
spice_astropy_frame_mapping = {
    'IAU_SUN': (suncoords.HeliographicCarrington,
                {'observer': suncoords.HeliographicStonyhurst(
                    0 * u.deg, 0 * u.deg, sunpy.sun.constants.radius)}),
}
SOLAR_ORBITER_ID = -144
SOLAR_ORBITER_SRF_FRAME_ID = -144000
SOLAR_ORBITER_STIX_ILS_FRAME_ID = -144851
SOLAR_ORBITER_STIX_OPT_FRAME_D = -144852

MIN_UNIX_TIM_LIMIT= stix_datetime.utc2unix('2020-02-10T05:00:00Z')



def linspace_datetime(start_utc,
                      end_utc,
                      num_steps=100,
                      max_num_steps=None,
                      min_time=None):
    ""
    ""
    start_unix = utc2unix(start_utc)
    end_unix = utc2unix(end_utc)
    if min_time:
        start_unix = max(start_unix, min_time)
        end_unix = max(end_unix, min_time)
    if max_num_steps:
        num_steps = min(max_num_steps, num_steps)
    ut = np.linspace(start_unix, end_unix, num_steps)
    return [datetime.fromtimestamp(u) for u in ut]


def vsep(v1,v2):
    """
        Find the separation angle in radians 
        between two double precision, 3-dimensional vectors. This angle is defined as zero if either vector is zero.
    """
    vector1=np.array([v1.x.value,v1.y.value,v1.z.value])
    vector2=np.array([v2.x.value,v2.y.value,v2.z.value])
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    angle = np.arccos(dot_product) #angle in radian
    return np.degrees(angle)



class Body:
    """
    A class for a single body.

    Parameters
    ----------
    body : `int` or `str`
        Either the body ID code or the body name.
    """
    def __init__(self, body):
        if isinstance(body, int):
            self.id = body
        elif isinstance(body, str):
            self.name = body
        else:
            raise ValueError('body must be an int or str')

    def __repr__(self):
        return f'{super().__repr__()}, name={self.name}, id={self.id}'

    def __eq__(self, other):
        return isinstance(other, Body) and other.id == self.id

    @property
    def id(self):
        """Body id code."""
        return self._id

    @id.setter
    def id(self, id):
        self._id = id
        try:
            self._name = sp.bodc2n(id)
        except spiceytypes.SpiceyError:
            raise ValueError(f'id "{id}" not known by SPICE')

    @property
    def name(self):
        """Body name."""
        return self._name

    @name.setter
    def name(self, name):
        self._name = name
        try:
            self._id = sp.bodn2c(name)
        except spiceytypes.SpiceyError:
            raise ValueError(f'Body name "{name}" not known by SPICE')


class Trajectory:
    """
    A class for the trajectory of a single body.

    Objects are initially created using only the body. To perform
    the actual trajectory calculation run :meth:`generate_positions`.
    The generated positions are then available via. the attributes
    :attr:`times`, :attr:`x`, :attr:`y`, and :attr:`z`.

    Parameters
    ----------
    target : str
        Name of the target. The name must be present in the loaded kernels.

    Notes
    -----
    When an instance of this class is created, a leapseconds kernel and a
    planets kernel are both automatically loaded.

    See also
    --------
    furnish : for loading in local spice kernels.
    """
    def __init__(self, target):
        self._target = Body(target)
        self._generated = False
        self.light_times=[]

        self.positions=[]
    def generate_positions(self, times, observing_body, frame,
                           abcorr=None):
        """
        Generate positions from a spice kernel.

        Parameters
        ----------
        times : time like
            An object that can be parsed by `~astropy.time.Time`.
        observing_body : str or int
            The observing body. Output position vectors are given relative to
            the position of this body. See
            https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
            for a list of bodies.
        frame : str
            The coordinate system to return the positions in. See
            https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html
            for a list of frames.
        abcorr : str, optional
            By default no aberration correciton is performed.
            See the documentaiton at
            https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html
            for allowable values and their effects.
        """
        times = time.Time(times)
        # Spice needs a funny set of times
        fmt = '%Y %b %d, %H:%M:%S'
        spice_times = [sp.str2et(time.strftime(fmt)) for time in times]
        self._abcorr = str(abcorr)

        # Do the calculation
        pos_vel, lightTimes = sp.spkezr(
            self.target.name, spice_times, frame, self._abcorr,
            observing_body)

        positions = np.array(pos_vel)[:, :3] * u.km
        velocities = np.array(pos_vel)[:, 3:] * u.km / u.s
        self.light_times=np.array(lightTimes)
        self.positions=np.array(pos_vel)[:,:3]

        self._frame = frame
        self._times = time.Time(times)
        self._velocities = velocities
        self._x = positions[:, 0]
        self._y = positions[:, 1]
        self._z = positions[:, 2]
        self._vx = velocities[:, 0]
        self._vy = velocities[:, 1]
        self._vz = velocities[:, 2]
        self._generated = True
        self._observing_body = Body(observing_body)

    @property
    def observing_body(self):
        '''
        Observing `Body`. The position vectors are all specified relative to
        this body.
        '''
        return self._observing_body

    @property
    def spice_frame(self):
        """
        The coordinate frame used by SPICE.
        """
        return self._spice_frame

    @property
    def times(self):
        '''
        A :class:`~astropy.time.Time` object containing the times sampled.
        '''
        return self._times

    @property
    def x(self):
        '''
        x coordinates of position.
        '''
        return self._x

    @property
    def y(self):
        '''
        y coordinates of position.
        '''
        return self._y

    @property
    def z(self):
        '''
        z coordinates of position.
        '''
        return self._z

    @property
    def r(self):
        '''
        Magnitude of position vectors.
        '''
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)

    @property
    def coords(self):
        """
        A :class:`~astropy.coordinates.SkyCoord` object.
        """
        if self._frame not in spice_astropy_frame_mapping:
            raise ValueError(f'Current frame "{self._frame}" not in list of '
                             f'known coordinate frames implemented in astropy '
                             f'or sunpy ({spice_astropy_frame_mapping})')
        if self._abcorr.lower() != 'none':
            raise NotImplementedError(
                'Can only convert to astropy coordinates if the aberration '
                'correction is set to "none" '
                f'(currently set to {self._abcorr})')

        frame = spice_astropy_frame_mapping[self._frame][0]

        # Override kwargs for sunpy < 2
        if (frame == suncoords.HeliographicCarrington and
                int(sunpy.__version__[0]) < 2):
            kwargs = {}

        kwargs = spice_astropy_frame_mapping[self._frame][1]
        coords = astrocoords.SkyCoord(
            self.x, self.y, self.z,
            frame=frame, representation_type='cartesian',
            obstime=self.times,
            **kwargs)
        coords.representation_type = frame().default_representation
        return coords

    @property
    def vx(self):
        """
        x component of velocity.
        """
        return self._vx

    @property
    def vy(self):
        """
        y component of velocity.
        """
        return self._vy

    @property
    def vz(self):
        """
        z component of velocity.
        """
        return self._vz

    @property
    def velocity(self):
        """
        Velocity.

        Returned as a shape ``(n, 3)`` array, where the ``n`` axis
        is the time axis.
        """
        return self._velocities

    @property
    def speed(self):
        '''
        Speed (magnitude of velocity vectors).
        '''
        return np.sqrt(self.vx**2 + self.vy**2 + self.vz**2)

    @property
    def generated(self):
        '''
        ``True`` if positions have been generated, ``False`` otherwise.
        '''
        return self._generated

    @property
    def target(self):
        '''
        The `Body` whose coordinates are being calculated.
        '''
        return self._target

    def change_units(self, unit):
        """
        Convert the positions to different units.

        Parameters
        ----------
        unit : astropy.units.Quantity
            Must be a unit of length (e.g. km, m, AU).
        """
        self._x = self._x.to(unit)
        self._y = self._y.to(unit)
        self._z = self._z.to(unit)




class SoloEphemeris(object):
    @staticmethod
    def compute_earth_sun_so_angle(solo_sun):

        size = solo_sun.x.value.size
        vec_earth_sun = np.array([1, 0, 0] * size).reshape(size, 3)
        vec_solo_sun = np.array([
            solo_sun.x.value, solo_sun.y.value, solo_sun.z.value
        ]).T
        fun = lambda v: math.sqrt(v.dot(v))
        mag = np.apply_along_axis(fun, 1, vec_solo_sun)
        product = np.array(
            [np.dot(v1, v2) for v1, v2 in zip(vec_earth_sun, vec_solo_sun)])
        return np.degrees(np.arccos(product / mag))

    @staticmethod
    def get_earth_spice_HEE(datetimes):
        ets=[sp.datetime2et(t) for t in datetimes]
        pos_vel, light_times= sp.spkezr('Earth', ets,      'SOLO_HEE_NASA', 'LT+S',   'Sun')
        positions = np.array(pos_vel)[:, :3] * u.km
        return {'light_times':np.array(light_times),'positions':positions.to(u.au)}

    @staticmethod
    def get_solo_ephemeris(start_utc,
                           end_utc,
                           num_steps=200):
        '''
          calculate solo orbiter orbit using spice kernel data
          Parameters:
            start_utc: start_utc string
            end_utc:   end utc string
            frame:    coordinate frame supported by SPICE
            num_steps:  number of data points. 
          Returns:
            orbit data which is a python dictionary
        '''
        observer='Earth'
        frame='SOLO_HEE_NASA'
        target='SOLO'
        orbiter_earth= Trajectory(target)
        orbiter_sun= Trajectory(target)
        earth_hee= Trajectory('Earth')
        #starttime = stix_datetime.utc2datetime(start_utc)
        start_unix=stix_datetime.utc2unix(start_utc)
        end_unix=stix_datetime.utc2unix(end_utc)

        start_unix=max(MIN_UNIX_TIM_LIMIT,start_unix)
        end_unix=max(start_unix, end_unix)
        num_steps=max(int((end_unix-start_unix)/(12*3600)), num_steps)
        ut_space=np.linspace(start_unix, end_unix, num_steps)

        times = []
        utc_times = []
        for t in ut_space:
            dt=stix_datetime.unix2datetime(t)
            times.append(dt)
            utc_times.append(dt.strftime("%Y-%m-%dT%H:%M:%SZ"))
        result = {}

        try:
            orbiter_earth.generate_positions(times, 'Earth', frame)
            orbiter_sun.generate_positions(times, 'SUN', frame)
            earth_hee.generate_positions(times, 'SUN', frame)

            orbiter_earth.change_units(u.au)
            orbiter_sun.change_units(u.au)


            solo_dist_to_earth = orbiter_earth.r.value

            sun_open_angle=const.R_sun.to(u.m)/orbiter_sun.r.to(u.m)
            sun_angular_diameter_arcmin=np.degrees(np.arctan(sun_open_angle.value))*60.*2
            lt_diff = earth_hee.light_times - orbiter_sun.light_times 
            earth_sun_solo_angles = SoloEphemeris.compute_earth_sun_so_angle(orbiter_sun)
            elevation= np.degrees(np.arctan2(orbiter_sun.z.value, orbiter_sun.r.value))
            #orientations=SoloEphemeris.get_solo_orientations(times)
            hee_pos=orbiter_sun.positions.tolist()
            #km

            result = {
                'ref_frame':frame,
                'observer':observer,
                'aunit': 'deg',
                'lunit': 'au',
                'vunit':'km/s',
                'tunit':'s',
                'utc': utc_times,
                'solo_hee':hee_pos,
                'x': -orbiter_sun.x.value,
                'y': -orbiter_sun.y.value,
                'z': -orbiter_sun.z.value,
                'sun_solo_r': orbiter_sun.r.value,
                'earth_solo_r': orbiter_earth.r.value,
                'speed': orbiter_sun.speed.value,
                #'orientation':orientations,
                'owlt': orbiter_earth.light_times,
                #'sun_earth_ltime':sun_earth_ltime,
                'light_time_diff': lt_diff,
                'earth_sun_solo_angle': earth_sun_solo_angles,
                'sun_angular_diameter':sun_angular_diameter_arcmin,
                'elevation': elevation,
        }
        except Exception as e:
            result = {'error': str(e)}

        return bson.dict_to_json(result)



    @staticmethod
    def get_flare_spice(sun_x,sun_y,flare_utc, observer='earth'):
        #sun_x, sun_y in unites of arcsec
        date_flare = stix_datetime.utc2datetime(flare_utc)
        et_stix = sp.datetime2et(date_flare)
        flare_coord= [sun_x, sun_y]
        # Obtain the coordinates of Solar Orbiter
        earth_hee_spice, ltime_earth_sun = sp.spkpos('EARTH', et_stix, 
                                         'SOLO_HEE_NASA', #  Reference frame of the output position vector of the object 
                                         'NONE', 'SUN')
        
        earth_hee_spice = earth_hee_spice * u.km
        # Convert the coordinates to HEE
        earth_hee = HeliocentricEarthEcliptic(earth_hee_spice, 
                                              obstime=Time(date_flare).isot, 
                                              representation_type='cartesian')
        solo_hee_spice,ltime_sun_solo = sp.spkpos('SOLO', et_stix, 'SOLO_HEE_NASA', 'NONE', 'SUN')
        solo_hee_spice = solo_hee_spice * u.km
        # Convert the coordinates to HEE
        solo_hee = HeliocentricEarthEcliptic(solo_hee_spice, 
                                         obstime=Time(date_flare).isot, 
                                         representation_type='cartesian')
        if observer=='earth':
            obs_coord=earth_hee#.transform_to(HeliographicStonyhurst(obstime=date_flare))
        else:
            obs_coord=solo_hee#.transform_to(HeliographicStonyhurst(obstime=date_flare) )
                                         
        flare_skycoord= SkyCoord(flare_coord[0]*u.arcsec, 
                              flare_coord[1]*u.arcsec,
                              obstime=date_flare,
                              observer=obs_coord,
                          frame='helioprojective')
        
        flare_hee = flare_skycoord.transform_to(HeliocentricEarthEcliptic(obstime=date_flare))
        v_flare_earth=earth_hee.cartesian-flare_hee.cartesian
        flare_hee_cart=flare_hee.cartesian
        flare_earth_r=np.sqrt(v_flare_earth.x**2+v_flare_earth.y**2+v_flare_earth.z**2)
        
        v_flare_solo=-flare_hee.cartesian+solo_hee.cartesian
        flare_solo_r=np.sqrt(v_flare_solo.x**2+v_flare_solo.y**2+v_flare_solo.z**2)
        ves=solo_hee.cartesian-earth_hee.cartesian
        v_earth_solo=np.sqrt(ves.x**2+ves.y**2+ves.z**2)
        
        flare_earth_t=flare_earth_r.to(u.m) / const.c
        flare_solo_t=flare_solo_r.to(u.m) / const.c
        owlt=v_earth_solo.to(u.m) / const.c

        flare_theta_stix=vsep(flare_hee_cart, v_flare_solo)
        flare_theta_earth=vsep(flare_hee_cart, v_flare_earth)

        return {'flare_earth_lt':flare_earth_t.value,
                'flare_solo_lt':flare_solo_t.value,
                'owlt':owlt.value,
                'dt': (flare_earth_t-flare_solo_t).value,
                'flare_solo_r':flare_solo_r.to(u.au).value,
                'dt_solar_center':ltime_earth_sun-ltime_sun_solo,
                'flare_utc':flare_utc,
                'theta_flare_norm_earth_deg':flare_theta_earth,
                'theta_flare_norm_solo_deg':flare_theta_stix,
                'earth_sun_ltime': ltime_earth_sun,
                'sun_solo_ltime': ltime_sun_solo,
                'observer':observer
    }
    @staticmethod
    def get_orientation(dt, frame1='SOLO_SRF', ref_frame='SOLO_SUN_RTN'):
        """
        Get the orientation or roll, pith and yaw of STIX (ILS or OPT).
        Taken from https://github.com/i4Ds/STIXCore/blob/9a765a33f2e924ead669b9b99afc1e41a4d2d8e8/stixcore
        /ephemeris/tests/test_position.py#L28-L40
        Parameters
        ----------
        date : `datetime.datetime`
            Date at which to obtain orientation information
        frame : `str`
            Name of the coordinate frame
        Returns
        -------
        tuple
            Roll, pitch and yaw of the spacecraft

        
   SOLO mission specific generic frames:

      SOLO_SUN_RTN                Sun Solar Orbiter Radial-Tangential-Normal
      SOLO_SOLAR_MHP              S/C-centred mirror helioprojective
      SOLO_IAU_SUN_2009           Sun Body-Fixed based on IAU 2009 report
      SOLO_IAU_SUN_2003           Sun Body-Fixed based on IAU 2003 report
      SOLO_GAE                    Geocentric Aries Ecliptic at J2000 (GAE)
      SOLO_GSE                    Geocentric Solar Ecliptic at J2000 (GSE)
      SOLO_HEE                    Heliocentric Earth Ecliptic at J2000 (HEE)
      SOLO_VSO                    Venus-centric Solar Orbital (VSO)

   Heliospheric Coordinate Frames developed for the NASA STEREO mission:

      SOLO_ECLIPDATE              Mean Ecliptic of Date Frame
      SOLO_HCI                    Heliocentric Inertial Frame
      SOLO_HEE_NASA               Heliocentric Earth Ecliptic Frame
      SOLO_HEEQ                   Heliocentric Earth Equatorial Frame
      SOLO_GEORTN                 Geocentric Radial Tangential Normal Frame

   Heliocentric Generic Frames(*):

      SUN_ARIES_ECL               Heliocentric Aries Ecliptic   (HAE)
      SUN_EARTH_CEQU              Heliocentric Earth Equatorial (HEEQ)
      SUN_EARTH_ECL               Heliocentric Earth Ecliptic   (HEE)
      SUN_INERTIAL                Heliocentric Inertial         (HCI)

   Geocentric Generic Frames:

      EARTH_SUN_ECL   (*)         Geocentric Solar Ecliptic     (GSE)
      EARTH_MECL_MEQX (*)         Earth Mean Ecliptic and Equinox of date
                                  frame (Auxiliary frame for EARTH_SUN_ECL)
      EARTH_MECL_MEQX_J2000       Earth Mean Ecliptic and Equinox at J2000
                                  frame (Auxiliary frame for SOLO_GSE and
                                  SOLO_HEE)


        """
        

        et = sp.datetime2et(dt)
        sc = sp.sce2c(SOLAR_ORBITER_ID, et) #convert to clock ticks
        #tol=sp.sctiks(int(sc), "1:000")
        tol= 1.0
        #cmat, sc= sp.ckgp(SOLAR_ORBITER_SRF_FRAME_ID, sc, tol, 'SOLO_ECLIP_NORM')
        #cmat, sc= sp.ckgp(SOLAR_ORBITER_SRF_FRAME_ID, sc, tol, 'SOLO_EQUAT_NORM')
        frame_id=SOLAR_ORBITER_SRF_FRAME_ID if frame1=='SOLO_SRF' else SOLAR_ORBITER_STIX_ILS_FRAME_ID
        cmat, sc= sp.ckgp(frame_id, sc, tol, ref_frame)
        #get c-matrix orientiation information
        roll, pitch, yaw = sp.m2eul(cmat, 1, 2, 3)
        #matrix to Euler angles
        # Following lines taken from from get_sunspice_roll.pro
        #https://www.heliodocs.com/xdoc/xdoc_list.php?dir=$SSW/packages/sunspice/idl
        roll, pitch, yaw = np.rad2deg([roll, pitch, yaw])#convert to deg
        #take from Frederic idl code
        #if yaw <0:
        #    yaw+=360

        """
        The lines below are taken from the IDL script at
        https://www.heliodocs.com/xdoc/xdoc_list.php?dir=$SSW/packages/sunspice/idl
        """
        pitch = -pitch
        if abs(roll) > 360: 
            roll=roll - math.copysign(360,roll)
        if abs(pitch) > 90:
            pitch= math.copysign(360,pitch)-pitch
            yaw=yaw-math.copysign(360, yaw)
            roll=roll-math.copysign(360, roll)

        #yaw = yaw %360
        #yaw=yaw - 180 #I don't know if it is right to do so


        return [roll, pitch, yaw]
    @staticmethod
    def get_solo_orientations(start_utc, end_utc, frame1='SOLO_SRF', ref_frame='SOLO_SUN_RTN', num_steps=200):
        datetimes=linspace_datetime(start_utc, end_utc, num_steps, MAX_STEPS)
        valid_att=[]
        valid_t=[]
        for dt in datetimes:
            try:
                res=SoloEphemeris.get_orientation(dt, frame1, ref_frame)
                valid_t.append(dt.timestamp()) 
                valid_att.append(res)
            except Exception as e: 
                logger.error(e)
        if not valid_att:
            return {'error': 'Data not available!'}
        euler_angles=np.array(valid_att).T 
        euler_angles=euler_angles.tolist()
        res={'unix_time':valid_t,'roll': euler_angles[0],'pitch': euler_angles[1],'yaw': euler_angles[2],'error':'',
                'frame1':frame1, 'ref_frame':ref_frame, 'units':'deg'}
        return res
    @staticmethod
    def to_stix_frame(rtn_coord, cmat_inv):
        """
        Sun center to STIX coordinates 
        """
        new_coord=np.dot(cmat_inv,rtn_coord)
        solo_r=np.sqrt(np.sum(rtn_coord**2))
        x_arcsec=np.arctan(new_coord[1]/solo_r )*180*3600/np.pi
        y_arcsec=np.arctan(new_coord[2]/solo_r )*180*3600/np.pi
        return [x_arcsec, y_arcsec]

    @staticmethod
    def aspect_correction(utc:str, stix_x:float,stix_y:float,
        ref_frame='SOLO_SUN_RTN'
            ):
        """
        coord_stix_frame: 
         [x0,y0]
         do aspect correct
        """
        dts=[stix_datetime.utc2datetime(utc)]
        times = time.Time(dts)
        fmt = '%Y %b %d, %H:%M:%S'
        spice_times = [sp.str2et(time.strftime(fmt)) for time in times]
        pos_vel, lightTimes = sp.spkezr('SUN', spice_times, ref_frame, 'None',
                    'SOLO')
        positions = np.array(pos_vel)[:, :3] * u.km
        et = sp.datetime2et(dts)
        sc = sp.sce2c(SOLAR_ORBITER_ID, et) #convert to clock ticks
        tol= 1.0
        frame_id= SOLAR_ORBITER_STIX_ILS_FRAME_ID
        cmat, _ = sp.ckgp(frame_id, sc, tol, ref_frame)
        sun_loc=np.sqrt(np.sum(positions**2))
        #rsun=const.R_sun.to(u.m).value
        coord_rtn=np.array([sun_loc,stix_x,stix_y])
        return SoloEphemeris.to_stix_frame(coord_rtn, cmat)

    @staticmethod
    def get_ephemeris_for_imaging(obs_utc):
        """
        calculate B0, L0, roll and radius for imaging software
        Parameters
        ----
        obs_utc: str
            observation time
        Returns
        ---
         B0, L0, roll, radius: float
            the first three in units of degrees and the last in arcsec
            
        """
        this_time = stix_datetime.utc2datetime(obs_utc)
        #solo_hee = coordinates_body(this_time, observer)
        stix_aux=SoloEphemeris.get_solar_limb_stix_fov(obs_utc)
        try:
            rsun_deg= np.degrees(np.arctan(stix_aux['rsun']/stix_aux['solo_sun_r']))* u.deg#in units of arcmin
            roll = stix_aux['roll']  #roll angle in degrees
            hee_spice=stix_aux['solo_hee']
            sun_center=stix_aux['sun_center']
        except KeyError:
            raise ValueError('No auxiliary data available for the requested time')

        hee_spice = np.array(hee_spice[0]) * u.km
        # Convert the coordinates to HEE
        solo_hee = HeliocentricEarthEcliptic(hee_spice,
                                              obstime=Time(this_time).isot,
                                              representation_type='cartesian')

        solo_hgs = solo_hee.transform_to(HeliographicStonyhurst(obstime=this_time))
        B0 = solo_hgs.lat.deg # Heliographic latitude (B0 angle)
        L0 = solo_hgs.lon.deg # Heliographic longitude (L0 angle)
        #ROLL ANGLE Solar Orbiter in degree
        return B0*u.deg, L0*u.deg, roll*u.deg,rsun_deg ,stix_aux['solo_hee']*u.km, stix_aux['solo_sun_r']*u.m,sun_center

    @staticmethod
    def get_solar_limb_stix_fov(utc, ref_frame='SOLO_SUN_RTN'):
        dts=[stix_datetime.utc2datetime(utc)]
        times = time.Time(dts)
        fmt = '%Y %b %d, %H:%M:%S'
        spice_times = [sp.str2et(time.strftime(fmt)) for time in times]
        pos_vel, lightTimes = sp.spkezr('SUN', spice_times, ref_frame, 'None',
                    'SOLO')
        positions = np.array(pos_vel)[:, :3] * u.km

        solo_hee_spice,ltime_sun_solo = sp.spkpos('SOLO',spice_times, 'SOLO_HEE_NASA', 'NONE', 'SUN')
        solo_hee_spice = solo_hee_spice.tolist()

        et = sp.datetime2et(dts)
        sc = sp.sce2c(SOLAR_ORBITER_ID, et) #convert to clock ticks
        tol= 1.0
        frame_id= SOLAR_ORBITER_STIX_ILS_FRAME_ID
        cmat, _ = sp.ckgp(frame_id, sc, tol, ref_frame)
        cmat_inv=linalg.inv(cmat)
        roll, pitch, yaw = sp.m2eul(cmat, 1, 2, 3)
        #matrix to Euler angles
        # Following lines taken from from get_sunspice_roll.pro
        #https://www.heliodocs.com/xdoc/xdoc_list.php?dir=$SSW/packages/sunspice/idl
        roll, pitch, yaw = np.rad2deg([roll, pitch, yaw])

        sun_loc=np.sqrt(np.sum(positions**2))
        rsun=const.R_sun.to(u.m).value
        sun_center=np.array([sun_loc.to(u.m).value, 0, 0])
        nsew_coords=[
                [sun_loc.to(u.m).value, 0, rsun],
                [sun_loc.to(u.m).value, 0, -rsun],
                [sun_loc.to(u.m).value,  -rsun,0],
                [sun_loc.to(u.m).value,  rsun,0],
                       ]

        limbs=[[sun_loc.to(u.m).value, rsun*np.cos(theta), rsun*np.sin(theta)] for theta in np.linspace(0,2*np.pi, 50)
                ]
        sun_center_stix_frame=SoloEphemeris.to_stix_frame(sun_center, cmat_inv)

        nsew_stix_frame=[SoloEphemeris.to_stix_frame(np.array(loc), cmat_inv) for loc in nsew_coords]

        limbs_stix_frame=np.array([SoloEphemeris.to_stix_frame(np.array(loc), cmat_inv) for loc in limbs])

        sun_outside_stix_fov=False if np.max(np.abs(limbs_stix_frame))<3600 else True
        pitch = -pitch
        if abs(roll) > 360: 
            roll=roll - math.copysign(360,roll)
        if abs(pitch) > 90:
            pitch= math.copysign(360,pitch)-pitch
            yaw=yaw-math.copysign(360, yaw)
            roll=roll-math.copysign(360, roll)

        return {
            'sun_center':sun_center_stix_frame, 
            'frame':'STIX_ILS',
            'ref_frame':ref_frame,
            'cmat':cmat.tolist(),
            'cmat_inv':cmat_inv.tolist(),
            'rsun':rsun, #units of m
            'solo_sun_r':sun_loc.to(u.m).value, #unit of m
            'time':utc,
            'fov':{'x':[-3600, 3600, 3600, -3600, -3600], 'y':[3600, 3600, -3600, -3600, 3600]},
            'limb':{'x':limbs_stix_frame[:,0].tolist(), 'y':limbs_stix_frame[:,1].tolist()},
            'nsew':nsew_stix_frame,
            'aunit':'arcsec',
            'sun_outside_stix_fov':sun_outside_stix_fov, 
            'roll':roll,
            'pitch':pitch,
            'yaw':yaw,
            'lunit':'m',
            'solo_hee':solo_hee_spice, #units of km
            'aunit':['arcsec','deg']
            }    

if __name__ == '__main__':
    from pprint import pprint
    #result = SoloEphemeris.get_solo_ephemeris('2020-04-10T00:00:00', '2020-04-11T00:00:00')
    #pprint(result)
    spm=spice_manager.spice
    dt=stix_datetime.utc2datetime('2021-07-15T00:00:00')
    spm.utc2scet('2021-10-15T00:00:00')
    ort=SoloEphemeris.get_orientation(dt)
    print(ort)
