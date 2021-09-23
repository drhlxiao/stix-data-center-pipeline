#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : parser_pipeline.py
# @description  : STIX packet parser pipeline. It detects new files in the specified folder,
#                 parses the packets and stores the decoded packets in the MongoDB
# @author       : Hualin Xiao
# @date         : Feb. 11, 2020
#

import os
import sys
import glob
from datetime import datetime
from stix.core import config
from stix.spice import stix_datetime
from stix.core import mongo_db
from stix.core import stix_logger
from stix.core import stix_parser
from stix.fits import fits_creator
from stix.analysis import calibration_chisquare
from stix.analysis import background_estimation as bkg
from stix.analysis import flare_detection
from stix.analysis import sci_packets_analyzer
from stix.analysis import integration_time_estimator
from stix.analysis import flare_goes_class as fgc
from stix.analysis import goes_downloader as gdl
from stix.spice import spice_manager as spm
from stix.flare_pipeline import flare_L1_analyzer as fla



S20_EXCLUDED = True
DO_CALIBRATIONS = True
ENABLE_FITS_CREATION = True
DO_BULK_SCIENCE_DATA_MERGING = True
FIND_FLARES = True
RUN_L1_FLARE_ANALYZER=True

SCI_PACKET_SPIDS = ['54114', '54115', '54116', '54117', '54143', '54125']
DO_BACKGROUND_ESTIMATION = True
ESTIMATE_ROTATING_BUFFER_TIME_BINS = True
daemon_config = config.get_config('pipeline.daemon')
noti_config = config.get_config('pipeline.daemon.notification')
mongodb_config = config.get_config('pipeline.mongodb')

MDB = mongo_db.MongoDB(mongodb_config['host'], mongodb_config['port'],
                       mongodb_config['user'], mongodb_config['password'])

HOST = config.HTTP_PREFIX

logger = stix_logger.get_logger()

goes=gdl.GOES()

def get_now():
    return datetime.now().isoformat()


def create_notification(raw_filename, service_5_headers, summary, num_flares, goes_class_list):
    file_id = summary['_id']
    start = stix_datetime.unix2utc(summary['data_start_unix_time'])
    end = stix_datetime.unix2utc(summary['data_stop_unix_time'])
    content = f'New file: {raw_filename}\nObservation time: {start} - {end} \nRaw packets: {HOST}/view/packet/file/{file_id}\n'
    try:
        if '54102' in summary['summary']['spid'] or '54101' in summary[
                'summary']['spid']:
            content += f'\nHousekeeping data: {HOST}/view/plot/housekeeping/file/{file_id}\n'
        if '54118' in summary['summary']['spid']:
            content += f'\nLight curves: {HOST}/view/plot/lightcurves?run={file_id}\n'
        content += f'\nL1A FITS files: {HOST}/view/list/fits/file/{file_id}\n'
        if summary['calibration_run_ids']:
            content += f'\nCalibration runs: {HOST}/view/plot/calibration/file/{file_id}\n'
        if [x for x in summary['summary']['spid'] if x in SCI_PACKET_SPIDS]:
            content += f'\nScience data: {HOST}/view/list/bsd/file/{file_id}\n'
    except Exception as e:
        logger.error(e)
    if service_5_headers:
        content += '\nSTIX Service 5 packets:\n'
        for header in service_5_headers:
            content += '\tAt {}, TM({},{}) {}\n'.format(
                header['UTC'], header['service_type'],
                header['service_subtype'], header['descr'])
    else:
        content += 'No Service 5 packet found in the file.\n'

    if num_flares > 0:
        content += '''\n{} solar flare(s) identified in the file\n \n'''.format(
            num_flares)
    else:
        content += '\n No solar flare detected.\n'
    if goes_class_list:
        content+='Peak UTC *  GOES class\n'
        for fl in goes_class_list:
            content+=f'{fl[0]}  -  {fl[1]} \n'


    doc = {
        'title': 'STIX operational message',
        'group': 'operations',
        'content': content,
        'time': stix_datetime.get_now(),
        'released': False,
        'is_sent': False,
        'file': file_id
    }
    _id = MDB.insert_notification(doc)
    return _id


def clear_ngnix_cache():
    '''
        remove ngnix cache if ngnix cache folder is defined in the configuration file
    '''
    files = glob.glob(daemon_config['ngnix_cache'])
    logger.info('Removing nginx cache..')
    for fname in files:
        try:
            os.remove(fname)
        except OSError as e:
            logger.error(str(e))
    logger.info('Nginx cache removed')

def process_one(filename):
    file_id = MDB.get_file_id(filename)
    notification_ids=[]
    if file_id == -2:
        summary=process('FM', filename, True, debugging=True)
        try:
            notification_ids.append(summary['notification_id'])
        except (TypeError, KeyError):
            pass
    MDB.release_notifications(notification_ids)


def process(instrument, filename, notification_enabled=True, debugging=False):
    spm.spice.load_kernels()
    #always load the latest kernel files
    base = os.path.basename(filename)
    name = os.path.splitext(base)[0]
    num_flares = 0
    log_path = daemon_config['log_path']
    log_filename = os.path.join(log_path, name + '.log')
    logger.set_logger(log_filename, level=3)
    if debugging:
        logger.enable_debugging()
    parser = stix_parser.StixTCTMParser()
    parser.set_MongoDB_writer(mongodb_config['host'], mongodb_config['port'],
                              mongodb_config['user'],
                              mongodb_config['password'], '', filename,
                              instrument)
    logger.info('{}, processing {} ...'.format(get_now(), filename))
    if S20_EXCLUDED:
        parser.exclude_S20()
    #parser.set_store_binary_enabled(False)
    parser.set_packet_buffer_enabled(False)
    service_5_headers = None
    goes_class_list=None

    try:
        parser.parse_file(filename)
        service_5_headers = parser.get_stix_alerts()
    except Exception as e:
        logger.error(str(e))
        return  None
    summary = parser.done()
    if not summary:
        return None

    file_id = summary['_id']

    if DO_BACKGROUND_ESTIMATION:
        logger.info('Background estimation..')
        try:
            bkg.process_file(file_id)
        except Exception as e:
            logger.error(str(e))

    if FIND_FLARES:
        logger.info('Searching for flares..')
        try:
            num_flares = flare_detection.find_flares(
                file_id, snapshot_path=daemon_config['flare_lc_snapshot_path'])
            if num_flares>0:
                goes_class_list=fgc.find_goes_class_flares_in_file(file_id)

            summary['num_flares']=num_flares
        except Exception as e:
            logger.error(str(e))
    if ESTIMATE_ROTATING_BUFFER_TIME_BINS:
        try:
            integration_time_estimator.process_file(file_id)
        except Exception as e:
            logger.error(str(e))



    if ENABLE_FITS_CREATION:
        logger.info('Creating fits files...')
        try:
            fits_creator.create_fits(file_id, daemon_config['fits_path'])
        except Exception as e:
            logger.error(str(e))

    if DO_BULK_SCIENCE_DATA_MERGING:
        logger.info(
            'merging bulk science data and preparing bsd json files...')
        try:
            sci_packets_analyzer.process_packets_in_file(file_id)
        except Exception as e:
            logger.error(str(e))
    if RUN_L1_FLARE_ANALYZER:
        flp=fla.FlareDataAnalyzer()
        try:
            flp.process_L1_BSD_in_file(file_id)
        except Exception as e:
            logger.error(str(e))

    if notification_enabled:
        logger.info('Creating notification...')
        try:
            notif_id = create_notification(base, service_5_headers, summary,
                                           num_flares, goes_class_list)
            summary['notification_id']=notif_id
        except Exception as e:
            logger.info(str(e))

    if DO_CALIBRATIONS:
        logger.info('Starting calibration spectrum analysis...')
        try:
            calibration_run_ids = summary['calibration_run_ids']
            report_path = daemon_config['calibration_report_path']
            cal=Calibration()
            for run_id in calibration_run_ids:
                cal.process_one_run(run_id,create_pdf=True, pdf_path=report_path)
        except Exception as e:
            logger.error(str(e))


    clear_ngnix_cache()
    return summary


def main():
    filelist = {}
    print('checking new files ...')
    num_processed = 0
    notification_ids=[]
    for instrument, selectors in daemon_config['data_source'].items():
        for pattern in selectors:
            filenames = glob.glob(pattern)
            for filename in filenames:
                if os.path.getsize(filename) == 0:
                    continue
                file_id = MDB.get_file_id(filename)
                if file_id == -2:
                    if instrument not in filelist:
                        filelist[instrument] = []
                    filelist[instrument].append(filename)
    for instrument, files in filelist.items():
        goes.download()
        for filename in files:
            print('Processing file:', filename)
            summary=process(instrument, filename , True)
            try:
                notification_ids.append(summary['notification_id'])
            except Exception as e: 
                print(e)
                pass
            num_processed += 1
    MDB.release_notifications(notification_ids)
    return num_processed


if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        process_one(sys.argv[1])
