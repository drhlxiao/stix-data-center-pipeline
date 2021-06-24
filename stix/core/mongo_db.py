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
from datetime import datetime,timedelta
import uuid
import bson
import pymongo
from pathlib import Path
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
            self.collection_data_requests = self.db['bsd']
            self.col_bsd_forms = self.db['bsd_req_forms']
            self.collection_fits = self.db['fits']
            self.collection_events = self.db['events']
            self.col_goes = self.db['goes']
            #self.col_flare_tbc= self.db['flare_tbc']
            self.collection_flares = self.db['flares']
            self.collection_qllc_statistics = self.db['qllc_statistics']
            self.collection_notifications = self.db['notifications']
            self.collection_spice = self.db['spice']
            self.collection_time_bins= self.db['time_bins']
            self.collection_aspect_solutions= self.db['aspect_solutions']

        except Exception as e:
            print('Error occurred while initializing mongodb: {}'.format(
                str(e)))

    def get_db(self):
        return self.db

    def get_bsd_req_form_by_uid(self, uid):
        return self.col_bsd_forms.find_one({'unique_ids': int(uid)})

    def get_collection(self, colname):
        try:
            return self.db[colname]
        except:
            return None

    def get_collection_calibration(self):
        return self.collection_calibration

    def get_collection_processing(self):
        return self.collection_raw_files

    def get_collection_packets(self):
        return self.collection_packets

    def get_collection_bsd(self):
        return self.collection_data_requests

    def is_connected(self):
        if self.db:
            return True
        else:
            return False

    def get_raw_info(self, run_id):
        if self.collection_raw_files:
            cursor = self.collection_raw_files.find_one({'_id': int(run_id)})
            return cursor
        return None

    def get_goes_x_ray_flux_file_list(self, start_unix, stop_unix):
        runs = []
        if self.col_goes:
            #query_string = {'$or':[{'start_unix':{'$gte': start_unix,  '$lt': stop_unix}},
            #    {'stop_unix': {  '$gte': start_unix, '$lt': stop_unix  }  }]}
            query_string = {
                'start_unix': {
                    '$lt': stop_unix
                },
                'stop_unix': {
                    '$gt': start_unix
                }
            }
            runs = self.col_goes.find(query_string).sort('_id', 1)
        return runs

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

    def get_packets_of_bsd_request(self, record_id, header_only=True):
        packets = []
        requests = list(
            self.collection_data_requests.find({
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
        return self.collection_packets.find({'_id': {'$in': packet_ids}})

    def get_filename_of_run(self, run_id):
        if self.collection_raw_files:
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
        #,'path':path})
        for x in cursor:
            return x['_id']

        return -2

    def select_packets_by_id(self, pid):
        if self.collection_packets:
            cursor = self.collection_packets.find({'_id': int(pid)})
            return cursor
        return []

    def delete_one_run(self, run_id):
        run_id = int(run_id)
        if self.collection_packets:
            cursor = self.collection_packets.delete_many({'run_id': run_id})

        if self.collection_raw_files:
            cursor = self.collection_raw_files.delete_many({'_id': run_id})

        if self.collection_calibration:
            cursor = self.collection_calibration.delete_many(
                {'run_id': run_id})

        if self.collection_qllc_statistics:
            self.collection_qllc_statistics.delete_many({'file_id': run_id})

        if self.collection_ql:
            cursor = self.collection_ql.delete_many({'run_id': run_id})

        if self.collection_data_requests:
            cursor = self.collection_data_requests.delete_many(
                {'run_id': run_id})
        if self.collection_flares:
            cursor = self.collection_flares.delete_many({'run_id': run_id})
        if self.collection_fits:
            fits_docs = self.collection_fits.find({'file_id': run_id})
            for doc in fits_docs:
                filename = os.path.join(doc['path'], doc['filename'])
                try:
                    os.remove(filename)
                except OSError as e:
                    pass
            cursor = self.collection_fits.delete_many({'file_id': int(run_id)})

        if self.collection_time_bins:
            self.collection_time_bins.delete_many({'file_id': int(run_id)})

    def delete_runs(self, runs):
        for run in runs:
            self.delete_one_run(run)

    def select_packets_by_run(self,
                              run_id,
                              SPIDs=[],
                              sort_field='_id',
                              order=1):
        if not isinstance(SPIDs, list):
            SPIDs = [SPIDs]
        pkts = []
        if self.collection_packets:
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
        if self.collection_raw_files:
            runs = self.collection_raw_files.find().sort('_id', order)
            return runs
        return None

    def set_run_ql_pdf(self, _id, pdf_filename):
        if self.collection_raw_files:
            run = self.collection_raw_files.find_one({'_id': _id})
            run['quicklook_pdf'] = pdf_filename
            self.collection_raw_files.save(run)


    def get_run_ql_pdf(self, _id):
        if self.collection_raw_files:
            run = self.collection_raw_files.find_one({'_id': _id})
            if 'quicklook_pdf' in run:
                return run['quicklook_pdf']
        return None

    def get_calibration_run_data(self, calibration_id):
        _id = int(calibration_id)
        rows = []
        if self.collection_calibration:
            rows = self.collection_calibration.find({'_id': _id})
        return rows

    def update_calibration_analysis_report(self, calibration_id, data):
        #save calibration data analysis results to mongodb
        _id = int(calibration_id)
        if self.collection_calibration:
            doc = self.collection_calibration.find_one({'_id': _id})
            doc['analysis_report'] = data
            self.collection_calibration.save(doc)

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
        if self.collection_fits:
            self.collection_fits.save(doc)

    def get_next_fits_id(self):
        try:
            return self.collection_fits.find({}).sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            return 0

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

        if self.collection_flares:
            first_id = self.get_next_flare_id()
        new_inserted_flares = []

        num_inserted = 0
        inserted_ids = [None] * result['num_peaks']
        for i in range(result['num_peaks']):
            if not result['is_major'][i]:
                #not major peak
                continue

            peak_unix = float(result['peak_unix_time'][i])
            time_window = 300
            hidden = False
            exists = self.collection_flares.find_one({
                'peak_unix_time': {
                    '$gt': peak_unix - time_window,
                    '$lt': peak_unix + time_window
                }
            })
            #if it is exists in db
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
            self.collection_flares.save(doc)
        #return new_inserted_flares
        result['inserted_ids'] = inserted_ids

    def set_flare_lc_filename(self, _id, lc_filename):
        doc = self.collection_flares.find_one({'_id': _id})
        if doc:
            doc['lc_path'] = os.path.dirname(lc_filename)
            doc['lc_filename'] = os.path.basename(lc_filename)
            self.collection_flares.replace_one({'_id': _id}, doc)

    def get_flare_joint_obs(self, _id, key):
        doc = self.collection_flares.find_one({'_id': _id})
        if doc:
            if key == 'wiki':
                return doc.get('wiki', None)
            if 'joint_obs' not in doc:
                return None
            return doc['joint_obs'].get(key, None)
        return None

    def update_flare_field(self, _id, key, value):
        doc = self.collection_flares.find_one({'_id': _id})
        if doc:
            doc[key] = value
            self.collection_flares.replace_one({'_id': _id}, doc)

    def update_flare_joint_obs(self, _id, key, value):
        doc = self.collection_flares.find_one({'_id': _id})
        if doc:
            if key != 'wiki':
                if 'joint_obs' not in doc:
                    doc['joint_obs'] = {}
                doc['joint_obs'][key] = value
            else:
                doc['wiki'] = value
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
        query_string = {
            "$and": [{
                'stop_unix_time': {
                    '$gt': start_unix_time
                }
            }, {
                'start_unix_time': {
                    '$lt': stop_unix_time
                }
            }, {
                'SPID': SPID
            }]
        }
        ret = collection.find(query_string, {'packet_id': 1}).sort('_id', 1)
        packet_ids = [x['packet_id'] for x in ret]

        if packet_ids:
            query_string = {'_id': {'$in': packet_ids}}
            cursor = self.collection_packets.find(query_string).sort(
                sort_field, 1)
            return cursor
        return []

    def get_file_spids(self, file_id):
        return [
            int(x) for x in self.collection_packets.distinct(
                'header.SPID', {'run_id': int(file_id)}) if x
        ]

    def delete_flares_of_file(self, run_id):
        if self.collection_flares:
            self.collection_flares.delete_many({'run_id': int(run_id)})

    def search_flares_by_tw(self,
                            start_unix,
                            duration,
                            num=MAX_NUM_RETURN_RECORDS,
                            threshold=0):
        #flares to be confirmed
        results = []
        if self.collection_flares:
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

    def insert_qllc_statistics(self, doc):
        next_id = 0
        try:
            next_id = self.collection_qllc_statistics.find({}).sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            pass
        doc['_id'] = next_id
        self.collection_qllc_statistics.save(doc)

    def get_latest_spice_kernels(self):
        types=list(self.collection_spice.distinct('type'))
        results={}
        for key in types:
            rows=list(self.collection_spice.find({'type':key}).sort('file_date', -1).limit(NUM_KERNEL_FILES_LIMIT))
            if rows:
                rows.reverse()
                results[key]=rows
        return results


    def insert_spice_kernel(self, filename):
        next_id = 0
        doc = {}

        pfilename = Path(filename)

        fdate = re.findall(r"\d{4}\d{2}\d{2}", filename)
        date_str = fdate[0] if fdate else '19700101'
        doc['file_date'] = datetime.strptime(date_str, "%Y%m%d")

        doc['path'] = pfilename.parent.as_posix()
        doc['filename'] = pfilename.name
        doc['type'] = pfilename.parent.name
        if self.collection_spice:
            is_found = self.collection_spice.find_one(
                {'filename': doc['filename']})
            if is_found:
                return
        try:
            next_id = self.collection_spice.find({}).sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            pass

        doc['_id'] = next_id
        doc['entry_time'] = datetime.now()
        self.collection_spice.save(doc)

    def get_nearest_qllc_statistics(self, start_unix, max_limit=500):
        try:
            right_closest = list(
                self.collection_qllc_statistics.find({
                    'start_unix': {
                        '$gte': start_unix
                    },
                    'is_quiet': True,
                    'max.0': {
                        '$lt': max_limit
                    }
                }).sort('start_unix', 1).limit(1))
            left_closest = list(
                self.collection_qllc_statistics.find({
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

    def insert_notification(self, doc, is_sent=False):
        next_id = 0
        try:
            next_id = self.collection_notifications.find({}).sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            pass
        doc['_id'] = next_id
        doc['is_sent'] = is_sent
        self.collection_notifications.save(doc)
    def insert_time_bins(self,doc, delete_existing):
        file_id=doc['file_id']
        if delete_existing:
            self.collection_time_bins.delete_many({'file_id':file_id})
        self.collection_time_bins.save(doc)

    

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
