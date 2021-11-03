""""
Quick Look fits file definitions
"""
from datetime import timedelta, datetime
import astropy.units as u
import numpy as np
from astropy.io import fits
# from stix_parser.io.fits.tm_to_fits import generate_filename
from astropy.table.operations import unique
from astropy.table.table import QTable

from stix.spice.stix_datetime import datetime_to_scet

from stix.core.stix_logger import get_logger
logger = get_logger(__name__)

sec_in_day = 24 * 60 * 60

'''
#commented by Hualin 
class FitsL0Processor:
    def __init__(self, archive_path):
        self.archive_path = archive_path
        

    def write_fits(self, product):
        # if product.type == '21-6' and product.control['ssid'][0] in [20]:
        #     for prod in product.to_requests():
        #         pass

        for prod in product.to_days():
            filename = self.generate_filename(prod, version=1)
            if getattr(prod, 'ssid', False):
                parts = [prod.level, prod.control['service_type'][0],
                     prod.control['service_subtype'][0], prod.ssid]
            else:
                parts = [prod.level, prod.control['service_type'][0],
                     prod.control['service_subtype'][0]]
            path = self.archive_path.joinpath(*[str(x) for x in parts])
            path.mkdir(parents=True, exist_ok=True)
            fitspath = path / filename
            if fitspath.exists():
                logger.info('Fits file %s exists appending data', fitspath.name)
                existing = prod.from_fits(fitspath)
                if np.abs([((len(existing.data['data'][i])/2) - (existing.control['data_len'][i]+7))
                           for i in range(len(existing.data))]).sum() > 0:
                    raise ValueError()
                logger.debug('Existing %s \n New %s', existing, prod)
                prod = prod + existing
                logger.debug('Combined %s', prod)

            # control = unique(prod.control, ['scet_coarse', 'scet_fine', 'seq_count'])
            # data = prod.data[np.isin(prod.data['control_index'], control['index'])]

            control = prod.control
            data = prod.data

            if np.abs([((len(data['data'][i]) / 2) - (control['data_len'][i] + 7))
                       for i in range(len(data))]).sum() > 0:
                raise ValueError()

            primary_header = self.generate_primary_header(prod, filename)
            primary_hdu = fits.PrimaryHDU()
            primary_hdu.header.update(primary_header)
            primary_hdu.header.update({'HISTORY': 'Processed by STIX'})

            control_hdu = fits.BinTableHDU(control)
            control_hdu.name = 'CONTROL'
            data_hdu = fits.BinTableHDU(data)
            data_hdu.name = 'DATA'
            hdul = fits.HDUList([primary_hdu, control_hdu, data_hdu])

            logger.info(f'Writing fits file to {path / filename}')
            hdul.writeto(path / filename, overwrite=True, checksum=True)

    def generate_primary_header(self, product, filename):
        """
        Generate primary header cards.

        Parameters
        ----------
        filename : str
            Filename

        Returns
        -------
        tuple
            List of header cards as tuples (name, value, comment)
        """
        headers = (
            # Name, Value, Comment
            ('TELESCOP', 'SOLO/STIX', 'Telescope/Sensor name'),
            ('INSTRUME', 'STIX', 'Instrument name'),
            ('OBSRVTRY', 'Solar Orbiter', 'Satellite name'),
            ('FILENAME', filename, 'FITS filename'),
            ('DATE', datetime.now().isoformat(timespec='milliseconds'),
             'FITS file creation date in UTC'),
            ('OBT_BEG', f'{product.control["scet_coarse"][0]}:{product.control["scet_fine"][0]}'),
            ('OBT_END', f'{product.control["scet_coarse"][-1]}:{product.control["scet_fine"][-1]}'),
            ('TIMESYS', 'OBT', 'System used for time keywords'),
            ('LEVEL', 'L)', 'Processing level of the data'),
            ('ORIGIN', 'STIX Team, FHNW', 'Location where file has been generated'),
            ('CREATOR', 'STIX-SWF', 'FITS creation software'),
            ('VERSION', 1, 'Version of data product'),
            ('OBS_MODE', 'Nominal '),
            ('VERS_SW', 1, 'Software version')
        )
        return headers

    def generate_filename(self, product, version, status=''):
        """
        Generate fits file name with SOLO conventions.

        Parameters
        ----------
        product : stix_parser.product.BaseProduct
            Product
        version : int
            Version of this product
        status : str
            Status of the packets

        Returns
        -------
        str
            The filename
        """
        user_req = getattr(product.control, 'request_id', '')
        if user_req:
            user_req = f'_{user_req}'

        sec_in_day = 24* 60 * 60
        scet_obs = (product.control["scet_coarse"][0] // sec_in_day) * sec_in_day
        ssid = getattr(product, 'ssid', '')
        if ssid:
            ssid = f"-{ssid}"
        name = f"{product.control['service_type'][0]}-{product.control['service_subtype'][0]}" \
               f"{ssid}"
        return f'solo_{product.level}_stix-{name}' \
               f'_{scet_obs}_V{version:02d}{status}.fits'

'''

class FitsL1Processor:
    def __init__(self, archive_path, db_id, ver, run_type):
        self.archive_path = archive_path

        if run_type=='file': #to be consistent with the old convention
            run_type=''
        else:
            run_type='_'+run_type

        self.run_type=run_type
        self.fits_db_id=db_id
        self.meta=[]
        self.version=ver
    def get_meta_data(self):
        return self.meta


    def write_fits(self, product):
        self.meta=[]
        if callable(getattr(product, 'to_days', None)):
            #print("# to_days")
            products = product.to_days()
        else:
            #print("# to_request")
            products = product.to_requests()
        for prod in products:
            #create an fits for each request
            filename = self.generate_filename(product=prod)
            start_date = prod.obs_beg.to_datetime()
            if prod.type == 'ql':
                start_date = prod.obs_avg.to_datetime()

            #path = self.archive_path.joinpath(*[prod.level, format(start_date.year, '04d'),
            #                                    format(start_date.month, '02d'),
            #                                    format(start_date.day, '02d'),
            #                                    prod.type.upper()])
            path=self.archive_path
            #path.mkdir(parents=True, exist_ok=True)
            fitspath = path / filename
            logger.info(f'Fits filename:{fitspath.as_posix()}')
            #if fitspath.exists():
            #    logger.info(f'Fits file %s exists appending data {fitspath.name}')
            #    existing = prod.from_fits(fitspath)
            #    # logger.debug('Existing %s \n New %s', existing, prod)
            #    prod = prod + existing
            #    # logger.debug('Combined %s', prod)

            control = prod.control
            data = prod.data

            elow, ehigh = prod.get_energies()

            energies = QTable()
            energies['channel'] = range(len(elow))
            energies['e_low'] = elow * u.keV
            energies['e_high'] = ehigh * u.keV

            # Convert time to be relative to start date
            data['time'] = (data['time'] - prod.obs_beg).to(u.s)

            primary_header = self.generate_primary_header(filename, prod)
            primary_hdu = fits.PrimaryHDU()
            primary_hdu.header.update(primary_header)
            primary_hdu.header.update({'HISTORY': 'Processed by STIX'})

            control_hdu = fits.BinTableHDU(control)
            control_hdu.name = 'CONTROL'
            data_hdu = fits.BinTableHDU(data)
            data_hdu.name = 'DATA'
            energy_hdu = fits.BinTableHDU(energies)
            energy_hdu.name = 'ENERGIES'

            hdul = fits.HDUList([primary_hdu, control_hdu, data_hdu, energy_hdu])

            #logger.debug(f'Writing fits file to {path / filename}')
            full_path=path/filename
            if full_path.is_file():
                print("Removing existing fits:", str(full_path))
                full_path.unlink()

            hdul.writeto(full_path, overwrite=True, checksum=True)
            _meta={'data_start_unix': prod.obs_beg.to_datetime().timestamp(),
                    'data_end_unix': prod.obs_end.to_datetime().timestamp(),
                    '_id': self.fits_db_id,
                    'filename':filename
                    }
            if 'request_id' in control.colnames:
                _meta['request_id']=int(control['request_id'][0])

            self.fits_db_id+=1
            self.meta.append(_meta)

    def generate_filename(self, product=None, status=''):
        """
        Generate fits file name with SOLO conventions.

        Parameters
        ----------
        product : stix_parser.product.BaseProduct
            Product
        version : int
            Version of this product
        status : str
            Status of the packets

        Returns
        -------
        str
            The filename
        """
        status = ''
        if status:
            status = f'_{status}'

        user_req = ''
        if 'request_id' in product.control.colnames:
            user_req = f"-{product.control['request_id'][0]}"

        tc_control = ''
        if 'tc_packet_seq_control' in product.control.colnames and user_req != '':
            tc_control = f'_{product.control["tc_packet_seq_control"][0]}'

        if product.type == 'ql':
            date_range = product.obs_avg.to_datetime().strftime("%Y%m%d")
        else:
            start_obs = product.obs_beg.to_datetime().strftime("%Y%m%dT%H%M%S")
            end_obs = product.obs_end.to_datetime().strftime("%Y%m%dT%H%M%S")
            date_range = f'{start_obs}-{end_obs}'

        fits_filename=f'solo_{product.level}_stix-{product.type}-' \
               f'{product.name.replace("_", "-")}{user_req}{self.run_type}' \
               f'_{date_range}_{self.fits_db_id:06d}_V{self.version:02d}.fits'
        return fits_filename
               #f'_{date_range}_V{version:02d}{status}{tc_control}_{self.fits_db_id}.fits'

    def generate_primary_header(self, filename, product):
        """
        Generate primary header cards.

        Parameters
        ----------
        filename : str
            Filename
        scet_coarse : numpy.ndarray
            SCET coarse time
        scet_fine : numpy.ndarray
            SCET fine time
        obs_beg : datetime.datetime
            Begging of observation
        obs_avg : datetime.datetime
           Averagea of observation
        obs_end : datetime.datetime
            End of observation

        Returns
        -------
        tuple
            List of header cards as tuples (name, value, comment)
        """
        headers = (
            # Name, Value, Comment
            ('TELESCOP', 'SOLO/STIX', 'Telescope/Sensor name'),
            ('INSTRUME', 'STIX', 'Instrument name'),
            ('OBSRVTRY', 'Solar Orbiter', 'Satellite name'),
            ('FILENAME', filename, 'FITS filename'),
            ('DATE', datetime.now().isoformat(timespec='milliseconds'),
             'FITS file creation date in UTC'),
            ('OBT_BEG', datetime_to_scet(product.obs_beg)),
            ('OBT_END', datetime_to_scet(product.obs_end)),
            ('TIMESYS', 'UTC', 'System used for time keywords'),
            ('LEVEL', 'L1A', 'Processing level of the data'),
            ('ORIGIN', 'STIX Team, FHNW', 'Location where file has been generated'),
            ('CREATOR', 'STIX-SWF', 'FITS creation software'),
            ('VERSION', self.version, 'Version of data product'),
            ('OBS_MODE', 'Nominal '),
            ('VERS_SW', 1, 'Software version'),
            ('DATE_OBS', product.obs_beg.fits,
             'Start of acquisition time in UT'),
            ('DATE_BEG', product.obs_beg.fits),
            ('DATE_AVG', product.obs_avg.fits),
            ('DATE_END', product.obs_end.fits),
            ('MJDREF', product.obs_beg.mjd),
            ('DATEREF', product.obs_beg.fits),
            ('OBS_TYPE', 'LC'),
            # TODO figure out where this info will come from
            ('SOOP_TYP', 'SOOP'),
            ('OBS_ID', 'obs_id'),
            ('TARGET', 'Sun')
        )
        return headers

# def energy_bands(num_energies):
#     """
#     Generate the energy band fits extension.
#
#     Parameters
#     ----------
#     num_energies : int
#         The number of energies
#
#     Returns
#     -------
#     astropy.io.fits.BinTableHDU
#         The Energy Extension
#
#     """
#     # Extension for energy channel data
#     energy_columns = (
#         fits.Column(name='CHANNEL', format='B', array=np.zeros(num_energies)),
#         fits.Column(name='E_MIN', format='E', unit='KeV', array=np.zeros(num_energies)),
#         fits.Column(name='E_MAX', format='E', unit='KeV', array=np.zeros(num_energies))
#     )
#
#     energy_coldefs = fits.ColDefs(energy_columns)
#     enenrgy_hdu = fits.BinTableHDU.from_columns(energy_coldefs)
#     enenrgy_hdu.name = 'ENEBAND'
#
#     return enenrgy_hdu
