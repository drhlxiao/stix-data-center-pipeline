"""
Low Latency Date Pipeline processing
"""

import logging
import shutil
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from os import getenv
from pathlib import Path

import numpy as np
from astropy.wcs import wcs

from parser.core import datatypes as sdt
from parser.core.parser import StixTCTMParser
from parser.datetime import datetime_to_scet
from parser.io.fits.processors import get_products, SPID_MAP
from parser.products.quicklook import LightCurve, FlareFlagAndLocation
from parser.utils import logger

Y_M_D_H_M = "%Y%m%d%H%M"
DIR_ENVNAMES = [('requests', 'instr_input_requests'), ('output', 'instr_output')]
REQUEST_GLOB = 'request_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9]*'
LLDP_VERSION = "00.07.00"


def get_paths():
    """
    Return paths for lldp processing.

    Returns
    -------
    dict
        Dictionary of paths
    """
    paths = {}
    for (name, envname) in DIR_ENVNAMES:
        path = getenv(envname)
        if path is None:
            raise ValueError('Environment variable with name "%s" was not found.' % (envname,))
        paths[name] = Path(path)

    paths['products'] = paths['output'].joinpath('products')
    paths['failed'] = paths['output'].joinpath('failed')
    paths['temporary'] = paths['output'].joinpath('temporary')
    paths['logs'] = paths['output'].joinpath('logs')
    return paths


PATHS = get_paths()

logger = logger.get_logger(__name__)
logger.propagate = False
file_handler = logging.FileHandler(PATHS['logs'] / 'lldp_processing.log')
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s - %(name)s',
                              datefmt='%Y-%m-%dT%H:%M:%SZ')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


class RequestException(Exception):
    """
    Exception during processing of request.
    """
    pass


def generate_lldp_filename(name, obs_beg, obs_end, curtime, status):
    """
    Generate fits filename using LLDP conventions.

    e.g. 'solo_LL01_stix-lightcurve_0628187840-0628189995_V201912191105C.fits'

    Parameters
    ----------
    name : str
        Name of data e.g lighcuve, variance, spectra
    obs_beg : datetime.datetime
        Start of observation
    obs_end : datetime.datetime
        End of observation
    curtime : datetime.datetime
        Current time
    status : str
        Status of data complete,  incomplete, or unknown (C, I U)
    Returns
    -------
    str
        Filename
    """
    obt_start = datetime_to_scet(obs_beg)
    obt_end = datetime_to_scet(obs_end)

    # Extract coarse time
    obt_start = obt_start.split(':')[0]
    obt_end = obt_end.split(':')[0]

    timestamp = curtime.strftime(Y_M_D_H_M)
    return f'solo_LL01_stix-{name}_{obt_start}-{obt_end}_V{timestamp}{status}.fits'


def generate_lldp_primary_header(filename, scet_coarse, scet_fine, obs_end, curtime, status):
    """
    Generate LLDP fits file primary header.

    Parameters
    ----------
    filename
    scet_coarse
    scet_fine
    obs_end
    curtime
    status

    Returns
    -------

    """

    headers = (
        # Name, Value, Comment
        ('TELESCOP', 'SOLO/STIX', 'Telescope/Sensor name'),
        ('INSTRUME', 'STIX', 'Instrument name'),
        ('OBSRVTRY', 'Solar Orbiter', 'Satellite name'),
        ('FILENAME', filename, 'FITS filename'),
        ('DATE', curtime.now().isoformat(timespec='milliseconds'),
         'FITS file creation date in UTC'),
        ('OBT_BEG', f'{scet_coarse}:{scet_fine}'),
        ('OBT_END', datetime_to_scet(obs_end)),
        ('TIMESYS', 'OBT', 'System used for time keywords'),
        ('LEVEL', 'LL01', 'Processing level of the data'),
        ('ORIGIN', 'Solar Orbiter SOC, ESAC', 'Location where file has been generated'),
        ('CREATOR', 'STIX-LLDP', 'FITS creation software'),
        ('VERSION', curtime.strftime(Y_M_D_H_M), 'Version of data product'),
        ('OBS_MODE', 'Nominal'),
        ('VERS_SW', LLDP_VERSION, 'Software version'),
        ('COMPLETE', status)
    )
    return headers


def process_request(request, outputdir):
    """
    Process at LLDP request.

    Parameters
    ----------
    request : `pathlib.Path`
        Path to directory containing request
    outputdir :
        Path to directory to store outputs

    Raises
    ------
    RequestException
        There was an error processing the request
    """
    errors = []
    parser = StixTCTMParser()
    logger.info('Processing %s', request.name)
    tmtc_files = list(request.joinpath('telemetry').glob('*.xml'))
    if len(tmtc_files) != 1:
        raise RequestException('Expected one tmtc file found %s.', len(tmtc_files))
    tree = ET.parse(tmtc_files[0])
    root = tree.getroot()
    packets = []
    for packet_node in root.findall('.//PktRawResponseElement'):
        # packet_id = packet_node.attrib.get('packetID')
        packet_hex = list(packet_node.getiterator('Packet'))[0].text.strip()
        packet = parser.parse_hex(packet_hex)[0]
        packet = sdt.Packet(packet)
        packets.append(packet)

    packets.sort(key=lambda x: (x.seq_count, x.obt_utc))

    # Only process light curve (54118) and flare flag and location (54122)
    complete_products, incomplete_products = get_products(packets, spids=[54122, 54118])

    for product, data in complete_products.items():
        product_name = [v for k, v in SPID_MAP.items() if k == product][0]
        for packets in data:
            parsed_packets = sdt.Packet.merge(packets, product)
            if parsed_packets['SPID'][0] == 54122:
                stx_prod = FlareFlagAndLocation(parsed_packets)
            elif parsed_packets['SPID'][0] == 54118:
                stx_prod = LightCurve(parsed_packets)
            else:
                continue

            cur_time = datetime.now()
            status = 'C'
            filename = generate_lldp_filename(product_name, stx_prod.obs_beg, stx_prod.obs_end,
                                              cur_time, status)
            primary_header = generate_lldp_primary_header(filename, stx_prod.scet_coarse,
                                                          stx_prod.scet_fine, stx_prod.obs_end,
                                                          cur_time, status)
            try:
                hdul = stx_prod.to_hdul()
            except ValueError as ve:
                logger.error(ve, exc_info=True)
                errors.append(ve)
                continue

            wcs = [('CRPIX1', 0.0, 'Pixel coordinate of reference point'),
                   ('CRPIX2', 0.0, 'Pixel coordinate of reference point'),
                   ('CDELT1', 2/3.0, '[deg] Coordinate increment at reference point'),
                   ('CDELT2', 2/3.0, '[deg] Coordinate increment a reference point'),
                   ('CTYPE1', 'YLIF-TAN', 'Coordinate type code'),
                   ('CTYPE2', 'ZLIF-TAN', 'Coordinate type code'),
                   ('CRVAL1', 0.0, '[deg] Coordinate value at reference point'),
                   ('CRVAL2', 0.0, '[deg] Coordinate value at reference point'),
                   ('CUNIT1', 'deg', 'Units of coordinate increment and value'),
                   ('CUNIT2', 'deg', 'Units of coordinate increment and value')]

            hdul[0].data = np.zeros((3, 3), dtype=np.uint8)
            hdul[0].header.update(primary_header)
            hdul[0].header.update(wcs)
            hdul[0].header.update({'HISTORY': 'Processed by STIX LLDP VM'})

            hdul.writeto(outputdir / filename, checksum=True)

    if errors:
        raise RequestException(errors)


if __name__ == "__main__":
    logger.root.setLevel(logging.DEBUG)
    logger.info('LLDP pipeline starting')
    while True:
        requests = list(PATHS['requests'].glob(REQUEST_GLOB))
        logger.info('Found %s requests to process.', len(requests))
        for request in PATHS['requests'].glob(REQUEST_GLOB):
            if PATHS['failed'].joinpath(request.name).exists():
                logger.info('%s found in failed (failed).', request.name)
            elif PATHS['temporary'].joinpath(request.name).exists():
                logger.info('%s found in temporary (processing or failed)', request.name)
            elif PATHS['products'].joinpath(request.name).exists():
                logger.info('%s found in products (completed)', request.name)
            else:
                logger.info('%s to be praocessed', request.name)
                temp_path = PATHS['temporary'] / request.name
                temp_path.mkdir(exist_ok=True, parents=True)
                try:
                    process_request(request, temp_path)
                except Exception as e:
                    logger.error('%s error while processing', request.name)
                    logger.error(e, exc_info=True)
                    failed_path = PATHS['failed']
                    shutil.move(temp_path.as_posix(), failed_path.as_posix())
                    continue

                shutil.move(temp_path.as_posix(), PATHS['products'].as_posix())
                logger.info('Finished Processing %s', request.name)

        time.sleep(5)
