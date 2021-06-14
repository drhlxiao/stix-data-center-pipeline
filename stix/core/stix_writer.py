#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @title        : stix_writer.py
# @description  : write decoded packets to mongodb, binary data files or python pickle files
# @author       : Hualin Xiao
# @date         : Feb. 27, 2019
import datetime
import gzip
import os
import pickle
import time
import math

import pymongo

from stix.core import stix_global
from stix.core import config
from stix.spice import stix_datetime
from stix.core import stix_logger
from stix.core import stix_metadata as meta
from stix.spice import spice_manager as spm

logger = stix_logger.get_logger()
MONGODB_CONFIG = config.get_config()['pipeline']['mongodb']


class StixPacketWriter(object):
    def write_all(self, packets):
        pass

    def write_one(self, packet):
        pass

    def set_filename(self, fname):
        pass

    def set_summary(self, summary):
        pass

    def close(self):
        pass

    def is_processed(self, filename):
        return False


class StixPickleWriter(StixPacketWriter):
    def __init__(self, filename):
        super(StixPickleWriter, self).__init__()
        self.filename = filename
        self.packet_counter = 0
        self.fout = None
        self.packets = []
        self.run = None
        if filename.endswith('.pklz'):
            self.fout = gzip.open(filename, 'wb')
        else:
            self.fout = open(filename, 'wb')

    def register_run(self,
                     in_filename,
                     filesize=0,
                     comment='',
                     idb_version=''):
        self.run = {
            'Input': in_filename,
            'Output': self.filename,
            'filsize': filesize,
            'comment': comment,
            'idb_version': idb_version,
            'Date': datetime.datetime.now().isoformat()
        }

    def write_all(self, packets):
        if self.fout:
            data = {'run': self.run, 'packet': packets}
            pickle.dump(data, self.fout)
            self.fout.close()

    def write_one(self, packet):
        self.packets.append(packet)

    def close(self):
        self.write_all(self.packets)


class StixBinaryWriter(StixPacketWriter):
    def __init__(self, filename):
        super(StixBinaryWriter, self).__init__()
        self.filename = filename
        self.packet_counter = 0
        self.fout = None
        self.packets = []
        self.num_success = 0
        try:
            self.fout = open(self.filename, 'wb')
        except IOError:
            logger.error('IO error. Can not create file:{}'.format(filename))

    def register_run(self,
                     in_filename,
                     filesize=0,
                     comment='',
                     idb_version=''):
        pass
        #not write them to binary file
    def get_num_sucess(self):
        return self.num_success

    def write_one(self, packet):
        if self.fout:
            try:
                raw = packet['bin']
                self.fout.write(raw)
                self.num_success += 1
            except KeyError:
                logger.warning('binary data not available')

    def write_all(self, packets):
        if self.fout:
            for packet in packets:
                self.write_one(packet)

    def close(self):
        if self.fout:
            self.fout.close()


class StixMongoDBWriter(StixPacketWriter):
    """write data to   MongoDB"""
    def __init__(self,
                 server=MONGODB_CONFIG['host'],
                 port=MONGODB_CONFIG['port'],
                 username=MONGODB_CONFIG['user'],
                 password=MONGODB_CONFIG['password']):
        super(StixMongoDBWriter, self).__init__()

        self.ipacket = 0
        self.packets = []
        self.start_unix = math.inf
        self.end_unix = 0
        self.start_scet = math.inf
        self.end_scet = 0
        self.db = None
        self.filename = ''
        self.path = ''
        self.summary = ''
        self.collection_packets = None
        self.current_packet_id = 0
        self.collection_raw_files = None
        self.current_run_id = 0
        self.start = -1
        self.end = -1
        self.run_info = None
        try:
            self.connect = pymongo.MongoClient(server,
                                               port,
                                               username=username,
                                               password=password)
            self.db = self.connect["stix"]
            self.collection_packets = self.db['packets']
            self.collection_raw_files = self.db['raw_files']
        except Exception as e:
            logger.error(str(e))

        self.science_report_analyzer = meta.StixScienceReportAnalyzer(self.db)

    def is_processed(self, in_filename):
        filename = os.path.basename(in_filename)
        abspath = os.path.abspath(in_filename)
        path = os.path.dirname(abspath)
        try:
            run = self.collection_raw_files.find_one({
                'path': path,
                'filename': filename
            })
            for x in run:
                return True
        except Exception as e:
            logger.error(str(e))

        return False

    def register_run(self,
                     in_filename,
                     filesize=0,
                     comment='',
                     idb_version='',
                     instrument=''):
        try:
            self.current_run_id = self.collection_raw_files.find().sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            self.current_run_id = 0
            # first entry

        try:
            self.current_packet_id = self.collection_packets.find().sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            self.current_packet_id = 0

        log_filename = logger.get_log_filename()
        self.filename = os.path.basename(in_filename)
        abspath = os.path.abspath(in_filename)
        self.path = os.path.dirname(abspath)

        self.run_info = {
            'filename': self.filename,
            'path': self.path,
            'comment': comment,
            'log': log_filename,
            'hidden': False,
            'creation_time': datetime.datetime.now(),
            'run_start_unix_time': time.time(),
            'data_start_scet': 0,
            'data_end_scet': 0,
            'run_stop_unix_time': 0,
            'data_start_unix_time': 0,
            'data_stop_unix_time': 0,
            '_id': self.current_run_id,
            'status': stix_global.UNKNOWN,
            'summary': '',
            'spice_sclk': spm.spice.get_last_sclk_filename(),
            'filesize': filesize,
            'instrument': instrument,
            'idb_version': idb_version
        }
        #print(self.run_info)

        self.inserted_run_id = self.collection_raw_files.insert_one(
            self.run_info).inserted_id

    def set_filename(self, fname):
        self.filename = os.path.basename(fname)
        abspath = os.path.abspath(fname)
        self.path = os.path.dirname(abspath)

    def set_summary(self, summary):
        self.summary = summary

    def write_all(self, packets):
        for packet in packets:
            self.write_one(packets)

    def write_one(self, packet):
        header_unix = packet['header']['unix_time']
        scet = packet['header'].get('SCET', 0)

        if stix_datetime.is_scet_valid(scet):
            if scet < self.start_scet:
                self.start_scet = scet
            if scet > self.end_scet:
                self.end_scet = scet
        if stix_datetime.is_unix_time_valid(header_unix):
            if header_unix < self.start_unix:
                self.start_unix = header_unix
            if header_unix > self.end_unix:
                self.end_unix = header_unix

        #insert header

        packet['run_id'] = self.current_run_id
        packet['_id'] = self.current_packet_id

        self.science_report_analyzer.start(self.current_run_id,
                                           self.current_packet_id, packet)

        try:
            self.collection_packets.insert_one(packet)
        except Exception as e:
            logger.error('Error occurred when inserting packet to MongoDB')
            logger.error(str(e))
            logger.info('Packet:' + str(packet['header']))
            raise
            return

        self.current_packet_id += 1
        self.ipacket += 1

    def close(self):
        if not self.collection_raw_files:
            logger.warning('MongoDB is not initialized ')
            return None
        logger.info('{} packets have been inserted into MongoDB'.format(
            self.ipacket))
        run = self.collection_raw_files.find_one({'_id': self.inserted_run_id})
        if run:
            if self.start_unix == math.inf:
                self.start_unix = 0
            if self.start_scet == math.inf:
                self.start_scet = 0
            run['data_start_unix_time'] = self.start_unix
            run['data_stop_unix_time'] = self.end_unix

            run['data_start_scet'] = self.start_scet
            run['data_end_scet'] = self.end_scet

            run['run_stop_unix_time'] = time.time()
            run['filename'] = self.filename
            run['path'] = self.path
            run['status'] = stix_global.OK
            #status ==1 if success  0
            run['summary'] = self.summary
            run['calibration_run_ids'] = self.science_report_analyzer.get_calibration_run_ids(
            )
            self.collection_raw_files.save(run)
            logger.info('File info updated successfully.')
            logger.info('File ID:{}'.format(run['_id']))
            return run
        logger.error('File info not updated.')
        return None
