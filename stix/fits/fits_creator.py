import sys
import os
import argparse
from collections import defaultdict
from datetime import datetime, timedelta
from itertools import chain
from pathlib import Path

from stix.core import datatypes as sdt

from stix.fits.io.processors import FitsL1Processor
from stix.fits.io import hk_fits_writer as hkw
from stix.spice import time_utils as st
from stix.fits.products.housekeeping import MiniReport, MaxiReport
from stix.fits.products.quicklook import LightCurve, Background, Spectra, Variance, \
    FlareFlagAndLocation, CalibrationSpectra, TMManagementAndFlareList
from stix.fits.products.science import XrayL0, Aspect, XrayL1, XrayL2, XrayL3, Spectrogram
#from stix.utils.logger import get_logger

from stix.core import mongo_db, logger
from stix.utils import checksum
logger = logger.get_logger()

db = mongo_db.MongoDB()
#logger = get_logger(__name__)
bsd_db = db.get_collection('bsd')
fits_db = db.get_collection('fits')
packets_db = db.get_collection('packets')

FITS_PATH = '/data/fits/'

DATA_LEVEL = 'L1A'


def create_fits_for_packets(file_id,
                            packets,
                            spid,
                            product,
                            is_complete,
                            base_path_name = FITS_PATH,
                            overwrite = True,
                            version = 3,
                            remove_duplicates=True,
                            run_type='file'):
    
    try:
        _create_fits_for_packets(file_id, packets, spid, product, is_complete,
                                 base_path_name, overwrite, version,
                                 remove_duplicates, run_type)
    except Exception as e:
        #raise e
        logger.error(e)


def _create_fits_for_packets(file_id,
                             packets,
                             spid,
                             product,
                             is_complete,
                             base_path_name=FITS_PATH,
                             overwrite=True,
                             version=1,
                             remove_duplicates=True,
                             run_type='file'):
    """
    Process a sequence containing one or more packets for a given product.

    Parameters
    ----------
    packets: list
        Packets
    spid : int
        SPID
    product : basestring
        Product name
    basepath : pathlib.Path
        Path
    pacekt_type:
        complete packets or incomplete packets
    overwrite : bool (optional)
        False (default) will raise error if fits file exits, True overwrite existing file
    """

    # For HK merge all stand alone packets in request
    if isinstance(file_id, str):
        file_id = int(file_id)
    if not packets:
        logger.warning('No packets found!')
        return
    logger.info('Merging packets...')

    parsed_packets = sdt.Packet.merge(packets,
                                      spid,
                                      value_type='raw',
                                      remove_duplicates=remove_duplicates)
    #eng_packets = sdt.Packet.merge(packets, spid, value_type='eng',remove_duplicates=remove_duplicates)
    eng_packets = None
    prod = None
    logger.info('Converting merged packets to fits...')

    try:
        if product == 'hk-mini':
            prod = MiniReport(parsed_packets)
            product_type = 'housekeeping'
        elif product == 'hk-maxi':
            prod = MaxiReport(parsed_packets)
            product_type = 'housekeeping'
        elif product == 'ql-lc':
            prod = LightCurve.from_packets(parsed_packets, eng_packets)
            product_type = 'quicklook'
        elif product == 'ql-bkg':
            prod = Background.from_packets(parsed_packets, eng_packets)
            product_type = 'quicklook'
        elif product == 'ql-spec':
            prod = Spectra.from_packets(parsed_packets, eng_packets)
            product_type = 'quicklook'
        elif product == 'ql-var':
            prod = Variance.from_packets(parsed_packets, eng_packets)
            product_type = 'quicklook'
        elif product == 'ql-ffl':
            prod = FlareFlagAndLocation.from_packets(parsed_packets)
            product_type = 'quicklook'
        elif product == 'ql-cal':
            prod = CalibrationSpectra.from_packets(parsed_packets, eng_packets)
            product_type = 'quicklook'
        elif product == 'ql-tm':
            prod = TMManagementAndFlareList.from_packets(parsed_packets)
            product_type = 'quicklook'
        elif product == 'xray-rpd':
            prod = XrayL0.from_packets(parsed_packets, eng_packets)
            product_type = 'science'
        elif product == 'xray-cpd':
            prod = XrayL1.from_packets(parsed_packets, eng_packets)
            product_type = 'science'
        elif product == 'xray-scpd':
            prod = XrayL2.from_packets(parsed_packets, eng_packets)
            product_type = 'science'
        elif product == 'xray-vis':
            prod = XrayL3.from_packets(parsed_packets, eng_packets)
            product_type = 'science'
        elif product == 'xray-spec':
            prod = Spectrogram.from_packets(parsed_packets, eng_packets)
            product_type = 'science'
        elif product == 'aspect':
            prod = Aspect.from_packets(parsed_packets, eng_packets)
            product_type = 'science'
        else:
            logger.warning(
                f'Not FITS writer for product {product}, SPID {spid}. Ignored')
            return
        #prod are a , containing a control and data in the dictionary type

        logger.info("Preparing data for writing fits file...")

        try:
            obs_time = parsed_packets['unix_time'][0]
            dt = st.unix2datetime(obs_time)
            sub_folder=dt.strftime("%Y%m")
        except Exception as e :
            logger.error(e)
            sub_folder='000000'

            

        base_path = Path(base_path_name)/sub_folder
        base_path.mkdir(parents=True, exist_ok=True)

        fits_uid= db.get_next_fits_id()
        metadata_entries = []
        #write extracted information to fits files
        if product_type == 'housekeeping':
            metadata = hkw.write_fits(base_path, fits_uid, prod, product,
                                      overwrite, version, run_type)
            metadata_entries = [metadata]
        else:
            logger.info(f"Creating fits structure for file {fits_uid}")
            fits_processor = FitsL1Processor(base_path, fits_uid, version,
                                             run_type)
            logger.info(f"Writing fits file {fits_uid}")
            fits_processor.write_fits(prod)
            metadata_entries = fits_processor.get_metadata()
        for metadata in metadata_entries:
            try:
                abs_path = base_path / metadata['filename']
                file_size = abs_path.stat().st_size
            except:
                file_size = 0

            md5=checksum.get_file_md5(abs_path)

            doc = {
                #'_id':fits_uid,
                'md5':md5,
                'packet_id_start': parsed_packets['min_id'],
                'packet_id_end': parsed_packets['max_id'],
                'packet_spid': spid,
                'num_packets': parsed_packets['num_packets'],
                'file_id': file_id,
                'product_type': product,
                'product_group': product_type,
                #'data_start_unix':meta['data_start_unix'],
                #'data_end_unix':meta['data_end_unix'],
                #'filename': meta['filename'],
                'complete': is_complete,
                'run_type': run_type,
                'version': version,
                'level': DATA_LEVEL,
                'creation_time': datetime.utcnow(),
                'path': base_path.as_posix(),
                'file_size': file_size

            }
            doc.update(metadata)
            print(doc)
            db.write_fits_index_info(doc)
            logger.info(f'Created  fits file:  {metadata["filename"]}')
    except Exception as e:
        #raise e
        logger.error(str(e))
        #raise e


def purge_fits_for_science_request(unique_id: int):
    cursor = fits_db.find({'request_id': unique_id})
    current_version = 1
    for cur in cursor:
        current_version = cur.get('version', 0)

        try:
            fits_filename = os.path.join(cur['path'], cur['filename'])
            logger.info(f'Removing older version of FITS file : {fits_filename}')
            os.unlink(fits_filename)
        except Exception as e:
            logger.warning(
                f'Failed to remove fits file:{fits_filename} due to: {str(e)}')
    logger.info(
        f'deleting fits collections for request with unique_id: {unique_id}')
    cursor = fits_db.delete_many({'request_id': unique_id})
    return current_version


def purge_fits_for_raw_file(file_id):
    print(f'Removing existing fits files for {file_id}')
    cursor = fits_db.find({'file_id': int(file_id)})
    for cur in cursor:
        try:
            fits_filename = os.path.join(cur['path'], cur['filename'])
            logger.info(f'Removing file: {fits_filename}')
            os.unlink(fits_filename)
        except Exception as e:
            logger.warning(
                f'Failed to remove fits file:{fits_filename} due to: {str(e)}'
            )
    logger.info(f'deleting fits collections for file: {file_id}')
    cursor = fits_db.delete_many({'file_id': int(file_id)})


def create_continous_low_latency_fits(start_unix,
                                      end_unix,
                                      output_path=FITS_PATH,
                                      overwrite=True,
                                      version=1,
                                      run_type='daily'):
    file_id = -1
    logger.info(f"data time span:{end_unix-start_unix}")
    for spid, product in sdt.LOW_LATENCY_TYPES.items():
        #print(spid, product, start_unix, end_unix)
        #if spid!=54118:
        #    print(" we only create lc this time >>>")
        #    continue

        if spid in sdt.QL_SPID_MAP.keys():
            packets = db.get_quicklook_packets(sdt.QL_SPID_MAP[spid],
                                               start_unix,
                                               end_unix - start_unix,
                                               sort_field='header.unix_time')
        elif spid == 54102 or spid == 54101:
            packets = packets_db.find({
                'header.unix_time': {
                    '$gte': start_unix,
                    '$lte': end_unix
                },
                'header.SPID': spid
            }).sort('header.unix_time', 1).max_time_ms(300 * 1000)
        create_fits_for_packets(file_id,
                                packets,
                                spid,
                                product,
                                is_complete=True,
                                base_path_name=output_path,
                                overwrite=overwrite,
                                version=version,
                                remove_duplicates=True,
                                run_type=run_type)


def create_daily_low_latency_fits(date, path=FITS_PATH, tpad=0, version=1):
    start_datetime = f'{date}T00:00:00'
    print("creating daily fits file for data on ", start_datetime)
    start_unix = st.utc2unix(start_datetime)
    end_unix = 86400 + start_unix
    create_continous_low_latency_fits(start_unix-tpad,
                                      end_unix+tpad,
                                      output_path=path,
                                      overwrite=True,
                                      version=version,
                                      run_type='daily')


def create_low_latency_fits_between_dates(date_start,
                                          date_end,
                                          path=FITS_PATH):
    start_datetime = f'{date_start}T00:00:00'
    start_unix = st.utc2unix(start_datetime)
    end_datetime = f'{date_end}T00:00:00'
    end_unix = st.utc2unix(end_datetime)
    while start_unix <= end_unix:
        create_continous_low_latency_fits(start_unix,
                                          start_unix + 86400,
                                          output_path=path,
                                          overwrite=True,
                                          version=1,
                                          run_type='daily')
        start_unix += 86400


def create_low_latency_fits_relative_days(relative_start,
                                          relative_end,
                                          path=FITS_PATH):
    if relative_start > 0:
        relative_start = -relative_start
    if relative_end > 0:
        relative_end = -relative_end

    start = datetime.utcnow() + timedelta(days=relative_start)
    end = datetime.utcnow() + timedelta(days=relative_end)
    date_start = start.strftime('%Y-%m-%d')
    date_end = end.strftime('%Y-%m-%d')
    create_low_latency_fits_between_dates(date_start, date_end, path)


def create_fits_for_bulk_science(bsd_id_start,
                                 bsd_id_end,
                                 output_path=FITS_PATH,
                                 overwrite=True,
                                 version=1):
    """
        create fits file for bulk science data 
        Parameters:
        bsd_id_start:  bulk science data entry id
        bsd_id_end:  bulk science data end id
    """
    bsd_docs = bsd_db.find({
        '_id': {
            '$gte': bsd_id_start,
            '$lte': bsd_id_end
        }
    }).max_time_ms(300 * 1000)
    for bsd in bsd_docs:
        pids = bsd['packet_ids']
        pkts = db.get_packets_by_ids(pids)
        logger.info(
            f'Creating fits file for bsd #{bsd["_id"]}, number of packets: {len(pids)}'
        )
        file_id = bsd['run_id']
        spid = bsd['SPID']
        product = sdt.SPID_MAP[spid]
        logger.info(f'{spid},{product}')
        create_fits_for_packets(file_id,
                                pkts,
                                spid,
                                product,
                                True,
                                output_path,
                                overwrite,
                                version,
                                remove_duplicates=True)


def create_daily_low_latency_fits_for_all(path=FITS_PATH, end_date_offset=7):
    last_date = db.get_last_daily_fits_file_date()
    if last_date:
        start_date_noon = last_date + timedelta(hours=12)
        end_date = datetime.now() - timedelta(hours=end_date_offset * 24)
        start_date_str = start_date_noon.strftime('%Y-%m-%d')
        logger.info(f'Creating daily fits file from {start_date_str}')
        end_date_str = end_date.strftime('%Y-%m-%d')
        create_low_latency_fits_between_dates(start_date_str, end_date_str,
                                              path)


def recreate_fits_for_science_data(start_date_offset: int,
                                   end_date_offset: int,
                                   path=FITS_PATH,
                                   purge_old=True,
                                   version=2):
    """
    recreate fits files for science data which were processed in the given time range
    """
    start_dt = datetime.now() - timedelta(hours=start_date_offset * 24)
    end_dt = datetime.now() - timedelta(hours=end_date_offset * 24)
    bsd_docs = bsd_db.find({
        'process_time': {
            '$gte': start_dt,
            '$lte': end_dt
        }
    }).max_time_ms(300 * 1000)
    for doc in bsd_docs:
        _id = doc['_id']
        unique_id = doc['unique_id']
        if purge_old:
            version = purge_fits_for_science_request(unique_id)
            version = version + 1
        create_fits_for_bulk_science(_id,
                                     _id,
                                     path,
                                     overwrite=True,
                                     version=version)


def create_fits(file_id, output_path, overwrite=True, version=1):
    if overwrite:
        purge_fits_for_raw_file(file_id)

    file_id = int(file_id)
    spid_packets = db.get_file_spids(file_id)
    #print(spid_packets)
    if not spid_packets:
        logger.warning(f'File {file_id} has no packet!')
        return
    for spid in spid_packets:
        spid = int(spid)
        if spid not in sdt.SPID_MAP.keys():
            logger.warning(f'Not supported spid : {spid}')
            continue
        logger.info(f'Requesting packets of file {file_id} from MongoDB')
        sort_field = 'header.unix_time'
        #if spid==54118:
        logger.info(f'Querying packets for SPID: {spid}')
        product = sdt.SPID_MAP[spid]
        #process science packets
        if spid in sdt.SCI_REPORT_SPIDS:
            bsd_docs = db.get_bsd_docs_by_run_id(file_id, spid)
            for bsd in bsd_docs:
                pids = bsd['packet_ids']
                pkts = db.get_packets_by_ids(pids)
                logger.info(f'Creating fits file for bsd #{bsd["_id"]}')
                create_fits_for_packets(file_id,
                                        pkts,
                                        spid,
                                        product,
                                        True,
                                        output_path,
                                        overwrite,
                                        version,
                                        remove_duplicates=True)

            continue

        #process low latency packets
        cursor = db.select_packets_by_run(file_id, [spid],
                                          sort_field=sort_field)

        is_complete = True
        #if housekeeping, process all of them once
        if spid in [54101, 54102]:
            create_fits_for_packets(file_id, cursor, spid, product,
                                    is_complete, output_path, overwrite,
                                    version)
            continue
        packets = []
        received_first = False
        received_last = False

        #iterate over packets of the same type
        hashes = []
        for i, pkt in enumerate(cursor):
            #process packets except HK and science data
            if pkt['hash'] in hashes:
                #remove duplicated packets
                continue
            hashes.append(pkt['hash'])

            seg_flag = int(pkt['header']['seg_flag'])

            if seg_flag == 3:  #'stand-alone packet':
                packets = [pkt]
                received_first = True
                received_last = True

            elif seg_flag == 1:  # 'first packet':
                packets = [pkt]

                received_first = True
                received_last = False

            elif seg_flag == 0:  #'continuation packet':
                if packets:
                    packets.append(pkt)
                else:
                    packets = [pkt]
            elif seg_flag == 2:  # 'last packet':
                received_last = True

                if packets:
                    packets.append(pkt)
                else:
                    packets = [pkt]
            if received_first and received_last:
                is_complete = True
                #print('complete packets')
                create_fits_for_packets(file_id, packets, spid, product,
                                        is_complete, output_path, overwrite,
                                        version)
                packets = []  #clean container after processing packets
                received_first = False
                received_last = False
        if packets:  #unprocessed packets
            #incomplete packets, QL packets are always incomplete
            is_complete = False
            create_fits_for_packets(file_id, packets, spid, product,
                                    is_complete, output_path, overwrite,
                                    version)
            logger.warning(
                f'Incomplete report {spid} (packets[0]["_id"]) found but fits file still created'
            )


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    required = ap.add_argument_group('Required arguments')
    optional = ap.add_argument_group('Optional arguments')

    optional.add_argument("-p",
                          dest='path',
                          default=FITS_PATH,
                          required=False,
                          help="Output fits path. ")

    optional.add_argument("-f",
                          dest='file_id',
                          default=None,
                          required=False,
                          help="File ID")
    optional.add_argument("-fs",
                          dest='file_id_range',
                          nargs=2,
                          default=None,
                          required=False,
                          help="File ID range")

    optional.add_argument("-d",
                          dest='date',
                          default=None,
                          required=False,
                          help="Start time to select LL data")
    optional.add_argument("-ds",
                          dest='date_range',
                          nargs=2,
                          default=None,
                          required=False,
                          help="End time to select LL data")
    optional.add_argument(
        "-da",
        dest='daily_all',
        nargs=2,
        default=None,
        required=False,
        help=
        "Generate daily FITS files for low latency data starting from the end date of last processing run "
    )

    optional.add_argument("-b",
                          dest='bsd_id',
                          default=None,
                          required=False,
                          help="bsd id")
    optional.add_argument("-bs",
                          dest='bsd_id_range',
                          nargs=2,
                          default=None,
                          required=False,
                          help="Bulk science data id range")
    optional.add_argument(
        "-br",
        dest='recreate_date_offsets',
        nargs=2,
        default=None,
        required=False,
        help=
        "Date offset range of fit files of science data requests to be recreated"
    )

    args = vars(ap.parse_args())
    path = args['path']

    if args['bsd_id']:
        bsd_id_start = int(args['bsd_id'])
        bsd_id_end = bsd_id_start
        logger.info(f'Creating fits for bsd #{bsd_id_start}')
        create_fits_for_bulk_science(bsd_id_start,
                                     bsd_id_end,
                                     path,
                                     overwrite=True,
                                     version=1)
    if args['bsd_id_range']:
        bsd_id_start = int(args['bsd_id_range'][0])
        bsd_id_end = int(args['bsd_id_range'][1])
        create_fits_for_bulk_science(bsd_id_start,
                                     bsd_id_end,
                                     path,
                                     overwrite=True,
                                     version=1)
    if args['file_id']:
        create_fits(int(args['file_id']), path, overwrite=True, version=1)
    if args['file_id_range']:
        for i in range(int(args['file_id_range'][0]),
                       int(args['file_id_range'][1]) + 1):
            create_fits(i, path, overwrite=True, version=1)
    if args['daily_all']:
        create_daily_low_latency_fits_for_all()

    if args['date_range']:
        date_start = args['date_range'][0]
        date_end = args['date_range'][1]
        create_low_latency_fits_between_dates(date_start, date_end, path)
    if args['date']:
        date = args['date']
        create_daily_low_latency_fits(date, path)
    if args['recreate_date_offsets']:
        start_off = int(args['recreate_date_offsets'][0])
        end_off = int(args['recreate_date_offsets'][1])
        recreate_fits_for_science_data(start_off,
                                       end_off,
                                       path=FITS_PATH,
                                       purge_old=True,
                                       version=2)
