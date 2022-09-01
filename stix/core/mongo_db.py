#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @title        : MongoDB.py
# @description  : Mongodb reader
# @author       : Hualin Xiao
# @date         : May. 12, 2019
#import json
import os
import re
import pprint
import uuid
import bson
import pymongo
from pathlib import Path
from datetime import datetime,timedelta
from stix.core import logger
logger = logger.get_logger()
NUM_KERNEL_FILES_LIMIT=10

NUM_MAX_PACKETS = 20000
MAX_REQUEST_LC_TIME_SPAN_DAYS = 3
MAX_NUM_RETURN_RECORDS = 20000
QL_SPIDS = {
    'lc': 54118,
    'bkg': 54119,
    'qlspec': 54120,
    'var': 54121,
    'flare': 54122
}


class MongoDB(object):
    def __init__(self, server='localhost', port=27017, user='', pwd=''):
        self.filename = None
        self.packets = []
        self.db = None
        self.collection_packets = None
        self.collection_raw_files = None
        self.collection_calibration = None
        self.collection_ql = None
        try:
            if server == 'localhost' and user == '' and pwd == '':
                self.connect = pymongo.MongoClient(server, port)
            else:
                self.connect = pymongo.MongoClient(server,
                                                   port,
                                                   username=user,
                                                   password=pwd,
                                                   authSource='stix')
            self.db = self.connect["stix"]
            self.collection_packets = self.db['packets']
            self.collection_raw_files = self.db['raw_files']
            self.collection_calibration = self.db['calibration_runs']
            self.collection_ql = self.db['quick_look']
            self.collection_science_data = self.db['bsd']
            self.collection_data_requests = self.db['data_requests']
            self.collection_fits = self.db['fits']
            self.collection_events = self.db['events']
            self.col_goes = self.db['goes_fluxes']
            self.collection_flares = self.db['flares']
            self.collection_lc_stats = self.db['lc_stats']
            self.collection_notifications = self.db['notifications']
            self.collection_time_bins= self.db['time_bins']
            #self.collection_aspect_solutions= self.db['aspect_solutions']
            self.collection_qlspectra=self.db['ql_spectra']
            self.collection_qllc=self.db['ql_lightcurves']
            self.collection_qloc=self.db['ql_flare_locs']
            self.collection_sw_config=self.db['sw_config']
            self.collection_flare_images= self.db['flare_images']
            self.collection_aspect= self.db['aspect']

            self.col_user_groups=self.db['user_groups']
        except Exception as e:
            raise IndexError('Error occurred while initializing mongodb: {}'.format(
                str(e)))

    def get_db(self):
        return self.db

    def get_group_users(self,group):
        return list(self.col_user_groups.find({'group': group}))



    def get_bsd_req_form_by_uid(self, uid):
        return self.collection_data_requests.find_one({'unique_ids': int(uid)})

    def get_collection(self, colname):
        try:
            return self.db[colname]
        except:
            return None

    def get_aspect_solutions(self, start_unix, end_unix):
        return self.collection_aspect.find({'unix_time':{'$gte':start_unix,'$lte':end_unix}})



    def get_collection_calibration(self):
        return self.collection_calibration

    def get_collection_processing(self):
        return self.collection_raw_files

    def get_collection_packets(self):
        return self.collection_packets

    def get_collection_bsd(self):
        return self.collection_science_data


    def is_connected(self):
        if self.db:
            return True
        else:
            return False

    def get_raw_info(self, run_id):
        return self.collection_raw_files.find_one({'_id': int(run_id)})

    def get_goes_fluxes(self, start_unix, stop_unix, satellite=16):
        res = []
        query_string = {
                'unix_time': {
                    '$lt': stop_unix, '$gt': start_unix
                },
                'satellite':satellite
            }
        return  self.col_goes.find(query_string).sort('unix_time', 1)

    def get_houskeeping_packets_of_file(self,file_id):
        return self.collection_packets.find({'_id':file_id, 'header.SPID':54102})

    def get_LC_pkt_by_tw(self, start_unix_time, span):
        if not self.collection_ql:
            yield []
        span = float(span)
        start_unix_time = float(start_unix_time)
        max_duration = 3600 * 24 * MAX_REQUEST_LC_TIME_SPAN_DAYS
        span = span if span <= max_duration else max_duration
        stop_unix_time = start_unix_time + span
        SPID = 54118
        query_string = {
            'stop_unix_time': {
                '$gt': start_unix_time
            },
            'start_unix_time': {
                '$lt': stop_unix_time
            },
            'SPID': SPID
        }
        ret = self.collection_ql.find(query_string, {
            'packet_id': 1
        }).sort('_id', 1)
        packet_ids = [x['packet_id'] for x in ret]
        if not packet_ids:
            yield []
        if packet_ids:
            for _id in packet_ids:
                #query_string = {'_id': {'$in': packet_ids}}
                yield self.collection_packets.find_one({'_id': _id})

    def get_bsd_docs_by_run_id(self, run_id, spid, complete_only=True):
        query={'run_id':int(run_id), 'SPID':int(spid)}
        if complete_only:
            query.update({'first_pkt':{'$gte':0}, 'last_pkt':{'$gte':0 }})
        docs=self.collection_science_data.find(query)
        return docs

    def get_packets_of_bsd_request(self, record_id, header_only=True):
        packets = []
        requests = list(
            self.collection_science_data.find({
                '_id': record_id
            }).limit(1))
        if not requests:
            return []
        request = requests[0]
        SPID = request['SPID']
        if request['SPID'] not in [54114, 54115, 54116, 54117, 54143, 54125]:
            return []
        packet_ids = request['packet_ids']
        if header_only:
            return self.collection_packets.find({'_id': {
                '$in': packet_ids
            }}, {'header': 1})
        #return self.collection_packets.find({'_id': {'$in': packet_ids}}).sort('header.unix_time',1)
        #removed the line below to improve performance
        packet_ids.sort()
        for pid  in packet_ids:
            yield self.collection_packets.find_one({'_id':  pid})


    def get_filename_of_run(self, run_id):
        cursor = self.collection_raw_files.find({'_id': int(run_id)})
        for x in cursor:
            return x['filename']
        return ''

    def get_raw_file_info(self, _id):
        cursor = self.collection_raw_files.find({'_id': int(run_id)})
        for x in cursor:
            return x
        return None

    def get_file_id(self, filename):
        if not self.collection_raw_files:
            return -1
        basename = os.path.basename(filename)
        abspath = os.path.abspath(filename)
        path = os.path.dirname(abspath)
        cursor = self.collection_raw_files.find({'filename': basename})
        for x in cursor:
            return x['_id']

        return -2

    def select_packets_by_id(self, pid):
        return self.collection_packets.find({'_id': int(pid)})

    def delete_one_run(self, run_id, keep_raw=False):
        run_id = int(run_id)
        cursor = self.collection_packets.delete_many({'run_id': run_id})

        if keep_raw:
            cursor = self.collection_raw_files.update_one({'_id': run_id},{'$set':{'hidden':True}})
        else:
            cursor = self.collection_raw_files.delete_many({'_id': run_id})

        cursor = self.collection_calibration.delete_many(
            {'run_id': run_id})

        self.collection_lc_stats.delete_many({'file_id': run_id})

        cursor = self.collection_ql.delete_many({'run_id': run_id})

        cursor = self.collection_science_data.delete_many(
                {'run_id': run_id})
        cursor = self.collection_flares.delete_many({'run_id': run_id})
        fits_docs = self.collection_fits.find({'file_id': run_id})
        for doc in fits_docs:
            filename = os.path.join(doc['path'], doc['filename'])
            try:
                os.remove(filename)
            except OSError as e:
                pass
        cursor = self.collection_fits.delete_many({'file_id': int(run_id)})

        self.collection_time_bins.delete_many({'file_id': int(run_id)})
        self.collection_qlspectra.delete_many({'run_id': int(run_id)})
        self.collection_qllc.delete_many({'run_id': int(run_id)})

    def delete_runs(self, runs, keep_raw=False):
        for run in runs:
            self.delete_one_run(run, keep_raw)

    def select_packets_by_run(self,
                              run_id,
                              SPIDs=[],
                              sort_field='_id',
                              order=1):
        if not isinstance(SPIDs, list):
            SPIDs = [SPIDs]
        pkts = []
        query_string = {'run_id': int(run_id)}
        if SPIDs:
            query_string = {
                'run_id': int(run_id),
                'header.SPID': {
                    '$in': SPIDs
                }
            }
        pkts = self.collection_packets.find(query_string).sort(
            sort_field, order)
        return pkts

    def close(self):
        if self.connect:
            self.connect.close()

    def select_all_runs(self, order=-1):
        runs = self.collection_raw_files.find().sort('_id', order)
        return runs

    def set_run_ql_pdf(self, _id, pdf_filename):
        run = self.collection_raw_files.find_one({'_id': _id})
        run['quicklook_pdf'] = pdf_filename
        self.collection_raw_files.save(run)


    def get_run_ql_pdf(self, _id):
        run = self.collection_raw_files.find_one({'_id': _id})
        if 'quicklook_pdf' in run:
            return run['quicklook_pdf']
        return None

    def get_calibration_run_data(self, calibration_id):
        _id = int(calibration_id)
        return self.collection_calibration.find({'_id': _id})

    def update_calibration_analysis_report(self, calibration_id, data):
        _id = int(calibration_id)
        self.collection_calibration.update_one({'_id': _id}, {'$set':{'analysis_report':data}})

    def get_calibration_runs_for_processing(self):
        #search for calibration runs which have  not been yet processed
        docs = []
        try:
            docs = list(self.collection_calibration.find())
        except Exception as e:
            pass
        if docs:
            return [
                doc['_id'] for doc in docs if ('spectra' in doc) and (
                    'sbspec_formats' in doc) and 'analysis_report' not in doc
            ]
        return []

    def write_fits_index_info(self, doc):
        self.collection_fits.save(doc)

    def get_next_fits_id(self):
        try:
            return self.collection_fits.find({}).sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            return 0
    def get_next_flare_image_id(self):
        try:
            return self.collection_flare_images.find({}).sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            return 0
    def insert_flare_image(self,doc):
        self.collection_flare_images.save(doc)


    def get_next_flare_id(self):
        try:
            return self.collection_flares.find().sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            return 0

    def save_flare_info(self, result):
        """
            write flare info into database
        """
        if 'run_id' not in result:
            return None
        if result['num_peaks'] == 0:
            return None

        try:
            cursor = self.collection_flares.delete_many(
                {'run_id': int(result['run_id'])})
            #delete
        except Exception as e:
            pass

        first_id = self.get_next_flare_id()

        new_inserted_flares = []

        num_inserted = 0
        inserted_ids = [None] * result['num_peaks']
        for i in range(result['num_peaks']):
            if not result['is_major'][i]:
                #not major peak
                continue

            peak_unix = float(result['peak_unix_time'][i])
            time_window = 90
            hidden = False
            exists = self.collection_flares.find_one({
                'peak_unix_time': {
                    '$gt': peak_unix - time_window,
                    '$lt': peak_unix + time_window
                }
            })
            #if two peaks 
            if exists:
                continue
            doc = {
                '_id': first_id + num_inserted,
                'hidden': hidden,
                'published': False
            }
            for key in result:
                if isinstance(result[key], list):
                    if len(result[key]) == result['num_peaks']:
                        doc[key] = result[key][i]
                else:
                    doc[key] = result[key]

            num_inserted += 1
            doc['peak_index'] = i
            new_inserted_flares.append(doc)

            inserted_ids[i] = doc['_id']
            logger.info(f'Inserting flare: {doc["_id"]}')
            self.collection_flares.insert_one(doc)
        #return new_inserted_flares
        result['inserted_ids'] = inserted_ids

    def set_flare_lc_filename(self, _id, lc_filename):
        doc = self.collection_flares.find_one({'_id': _id})
        if doc:
            doc['lc_path'] = os.path.dirname(lc_filename)
            doc['lc_filename'] = os.path.basename(lc_filename)
            self.collection_flares.replace_one({'_id': _id}, doc)

    def get_flare_pipeline_products(self, _id, key):
        doc = self.collection_flares.find_one({'_id': _id})
        if doc:
            if key == 'wiki':
                return doc.get('wiki', None)
            if 'pipeline' not in doc:
                return None
            return doc['pipeline'].get(key, None)
        return None

    def update_flare_info(self, _id, key, value):
        doc = self.collection_flares.find_one({'_id': _id})
        if doc:
            doc[key] = value
            self.collection_flares.replace_one({'_id': _id}, doc)

    def update_flare_pipeline_products(self, _id, key, value):
        doc = self.collection_flares.find_one({'_id': _id})
        if doc:
            if 'sim_obs' not in doc:
                doc['sim_obs'][key] = value
            self.collection_flares.replace_one({'_id': _id}, doc)

    def get_quicklook_packets(self,
                              packet_type,
                              start_unix_time,
                              span,
                              sort_field='_id'):
        span = float(span)
        start_unix_time = float(start_unix_time)
        if span > 3600 * 24 * MAX_REQUEST_LC_TIME_SPAN_DAYS:  #max 3 days
            return []
        stop_unix_time = start_unix_time + span
        SPID = QL_SPIDS[packet_type]
        collection = self.collection_ql
        if not collection:
            return []
        query_string ={
                'stop_unix_time': {
                    '$gt': start_unix_time
                },
                'start_unix_time': {
                    '$lt': stop_unix_time
                },
                'SPID': SPID
            }
        ret = collection.find(query_string, {'packet_id': 1}).sort('_id', 1)
        packet_ids = [x['packet_id'] for x in ret]

        if packet_ids:
            query_string = {'_id': {'$in': packet_ids}}
            cursor = self.collection_packets.find(query_string).sort(
                sort_field, 1)
            return cursor
        return []


    def get_telecommands_in_file(self,file_id):
        return self.collection_packets.find({'header.TMTC':'TC', 'run_id': int(file_id)})


    def get_file_spids(self, file_id):
        return [
            int(x) for x in self.collection_packets.distinct(
                'header.SPID', {'run_id': int(file_id)}) if x
        ]

    def delete_flares_of_file(self, run_id):
        self.collection_flares.delete_many({'run_id': int(run_id)})


    def search_flares_by_tw(self,
                            start_unix,
                            duration,
                            num=MAX_NUM_RETURN_RECORDS,
                            threshold=0):
        #flares to be confirmed
        results = []
        query_string = {
            'peak_unix_time': {
                '$gte': start_unix,
                '$lt': start_unix + duration
            },
            'peak_counts': {
                '$gte': threshold
            }
        }
        results = self.collection_flares.find(query_string).sort(
                'peak_unix_time', -1).limit(num)
        return results

    def insert_lc_stats(self, doc):
        next_id = 0
        try:
            next_id = self.collection_lc_stats.find({}).sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            pass
        doc['_id'] = next_id
        self.collection_lc_stats.save(doc)

 
    def get_nearest_lc_stats(self, start_unix, max_limit=500):
        try:
            right_closest = list(
                self.collection_lc_stats.find({
                    'start_unix': {
                        '$gte': start_unix
                    },
                    'is_quiet': True,
                    'max.0': {
                        '$lt': max_limit
                    }
                }).sort('start_unix', 1).limit(1))
            left_closest = list(
                self.collection_lc_stats.find({
                    'start_unix': {
                        '$lte': start_unix
                    },
                    'is_quiet': True,
                    'max.0': {
                        '$lt': max_limit
                    }
                }).sort('start_unix', -1).limit(1))
            if right_closest and not left_closest:
                return left_closest[0]
            if not right_closest and left_closest:
                return left_closest[0]
            if not right_closest and not left_closest:
                return None
            left_unix = left_closest[0]['start_unix']
            right_unix = right_closest[0]['start_unix']
            if start_unix - left_unix > right_unix - start_unix:
                return right_closest[0]
            return left_closest[0]
        except:
            return None


    def get_last_daily_fits_file_date(self):
        flist=list(self.collection_fits.find({'run_type':'daily'}).sort('data_end_unix',-1).limit(1))
        if flist:
            return datetime.fromtimestamp(flist[0]['data_end_unix'])
        return None

    def insert_time_bins(self,doc, delete_existing):
        file_id=doc['file_id']
        if delete_existing:
            self.collection_time_bins.delete_many({'file_id':file_id})
        self.collection_time_bins.save(doc)
    def get_bsd_fits_info_by_request_id(self, request_id, max_num_fits=1000):
        return self.collection_fits.find({
                'request_id': request_id
            }).sort('_id', -1).limit(max_num_fits)
    def find_L1_background(self, unix_time, max_days_off=14, emin=3, emax=17):
        start_unix=unix_time-max_days_off*86400
        query=[{
                '$match': {'request_type':'L1','purpose':'Background','detector_mask':'0xFFFFFFFF', 
                'emax':{'$gte':str(emax)},
                'emin':{'$lte':str(emin)},
                'start_unix':{'$gt':start_unix, '$lt':unix_time}
            }},
            {
                '$lookup': {
                    'from': 'fits',
                    'localField': 'unique_ids.0',
                    'foreignField': 'request_id',
                    'as': 'merg'
                }
            },  {
                '$match': {
                    'merg.request_id': {
                        '$exists': True
                    }
                }
                },
              {'$sort':{'start_unix':-1}
            }
            ]
        return self.collection_data_requests.aggregate(query)
       
    def get_packets_by_ids(self,id_list):
        for _id in id_list:
            yield self.collection_packets.find_one({'_id': _id})

    def find_flares_by_time_range(self, start_unix, end_unix):
        return self.collection_flares.find({'start_unix':{'$lt':end_unix},
            'end_unix':{'$gt':start_unix}} ).sort('start_unix',1)
        

    def get_quicklook_packets_of_run(self, packet_type, run):
        collection = None
        SPID = QL_SPIDS[packet_type]
        query_string = {'run_id': run, 'header.SPID': SPID}
        cursor = self.collection_packets.find(query_string).sort('_id', 1)
        return cursor
        return []


#if __name__ == '__main__':
#    mdb = MongoDB()
#    #print(mdb.get_packet_for_header(318))
#    mdb.set_quicklook_pdf(0, '/data/a.pdf')
