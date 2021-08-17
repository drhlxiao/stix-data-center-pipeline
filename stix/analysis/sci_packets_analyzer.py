#!/usr/bin/python3
# author: Hualin Xiao
# pre-process science data, merge bulk science data packets and write merged data to json files
# so that web client side could load the data quickly
import sys
import os
import json
import numpy as np
from datetime import datetime
from stix.core import stix_datatypes as sdt
from stix.spice import stix_datetime
from stix.core import mongo_db as db
from stix.core import stix_logger
from stix.core import config
mdb = db.MongoDB()
logger = stix_logger.get_logger()
level1_products_path = config.get_config(
    'pipeline.daemon.level1_products_path')

DATA_REQUEST_REPORT_SPIDS = [54114, 54115, 54116, 54117, 54143, 54125]
DATA_REQUEST_REPORT_NAME = {
    54114: 'L0',
    54115: 'L1',
    54116: 'L2',
    54117: 'L3',
    54143: 'L4',
    54125: 'ASP'
}

PROCESS_METHODS = {
    'normal': [54114, 54117, 54143, 54125],
    'yield': [54115, 54116]
}


def get_process_method(spid):
    for key, val in PROCESS_METHODS.items():
        if spid in val:
            return key
    return 'unknown'


class StixBulkL0Analyzer(object):
    def __init__(self):
        self.request_id = -1
        self.packet_unix = 0
        self.start_unix=None
        self.end_unix=None

        self.auxiliary = {
            'pixel_mask': [],
            'detector_mask': [],
            'rcr': [],
            'dt': [],
            'triggers': [],
            'time': []
        }
        self.boxes = {
            'detector': [],
            'pixel': [],
            'energy': [],
            'time': [],
            'counts': [],
        }

    def format_report(self):
        report = {
            'packet_unix': self.packet_unix,
            'boxes': self.boxes,
            'auxiliary': self.auxiliary,
            'request_id':self.request_id,
            'start_unix':self.start_unix,
            'end_unix':self.end_unix
        }
        return report

    def merge(self, cursor):

        for pkt in cursor:
            packet = sdt.Packet(pkt)
            self.request_id = packet[3].raw
            # print(self.request_id)
            self.packet_unix = packet['unix_time']
            T0 = stix_datetime.scet2unix(packet[12].raw)

            num_structures = packet[13].raw
            eacc_SKM = (packet[5].raw, packet[6].raw, packet[7].raw)
            trig_SKM = (packet[9].raw, packet[10].raw, packet[11].raw)

            trig_idx = 1
            # uncompressed counts
            if sum(trig_SKM) > 0:
                trig_idx = 2

            counts_idx = 1  # raw value
            if sum(eacc_SKM) > 0:
                counts_idx = 2
            children = packet[13].children

            for i in range(0, num_structures):
                offset = i * 23
                unix_time = children[offset][1] * 0.1 + T0
                if self.start_unix is None:
                    self.start_unix=unix_time
                self.end_unix=unix_time

                self.auxiliary['time'].append(unix_time)
                self.auxiliary['rcr'].append(children[offset + 1][1])
                dt = children[offset + 2][1] * 0.1
                self.auxiliary['dt'].append(dt)
                self.auxiliary['pixel_mask'].append(children[offset + 4][1])
                self.auxiliary['detector_mask'].append(children[offset + 5][1])
                self.auxiliary['triggers'].append(
                    [children[k + 6][trig_idx] for k in range(0, 16)])
                # id 6-22, 16 trigger accumulators
                num_samples = children[22][1]
                samples = children[22][3]
                for j in range(0, num_samples):
                    k = j * 5
                    pixel = samples[k + 1][1]
                    detector = samples[k + 2][1]
                    energy_bin = samples[k + 3][1]
                    num_bits = samples[k + 4][1]
                    continous_counts = samples[k + 4][3]
                    counts = 1
                    if num_bits == 1:
                        counts = continous_counts[0][counts_idx]
                    if num_bits == 2:
                        counts = continous_counts[0][counts_idx] << 8
                        counts += continous_counts[1][counts_idx]
                    self.boxes['detector'].append(detector)
                    self.boxes['pixel'].append(pixel)
                    self.boxes['energy'].append(energy_bin)
                    self.boxes['time'].append(unix_time)
                    self.boxes['counts'].append(counts)

        return self.format_report()


class StixBulkL1L2Analyzer(object):
    def __init__(self):
        self.time = []
        self.rcr = []
        self.dt = []
        self.pixel_mask = []
        self.detector_mask = []
        self.triggers = []
        self.hits = []
        #self.hits_time_int = np.zeros((33, 12, 32))  # total hits
        #self.hits_time_ene_int = np.zeros(
        #    (33, 12))  # energy integrated total hits
        # self.T0=[]
        self.eacc_SKM = ()
        self.trig_SKM = ()
        self.request_id = -1
        self.packet_unix = 0
        self.groups = []
        self.pixel_total_counts = [0] * 384
        self.num_time_bins = []
        self.start_unix=None
        self.end_unix=None


    def format_report(self):
        report = {
            'request_id': self.request_id,
            'packet_unix': self.packet_unix,
            'groups': self.groups,
            'trig_skm': self.trig_SKM,
            'eacc_skm': self.eacc_SKM,
            'start_unix': self.start_unix,
            'end_unix': self.end_unix,
            'pixel_total_counts': self.pixel_total_counts,
        }
        return report

    def get_spectrum_pixel_indexes(self, detector_mask, pixel_mask_array):
        detectors = [0] * 32
        pixel_mask = 0
        for p in pixel_mask_array:
            pixel_mask |= p
        pixels = [0] * 12
        for i in range(32):
            if (detector_mask & (1 << i)) != 0:
                detectors[i] = 1
        for j in range(12):
            if (pixel_mask & (1 << j)) != 0:
                pixels[j] = 1

        pixel_indexes = []
        for i in range(32):
            for j in range(12):
                if detectors[i] == 0 or pixels[j] == 0:
                    continue
                pixel_indexes.append(i * 12 + j)
        return pixel_indexes

    def merge(self, cursor):
        for pkt in cursor:
            packet = sdt.Packet(pkt)
            self.request_id = packet[3].raw
            self.packet_unix = packet['unix_time']
            T0 = stix_datetime.scet2unix(packet[12].raw)
            num_structures = packet[13].raw
            self.eacc_SKM = (packet[5].raw, packet[6].raw, packet[7].raw)
            self.trig_SKM = (packet[9].raw, packet[10].raw, packet[11].raw)

            trig_idx = 1
            # uncompressed counts
            if sum(self.trig_SKM) > 0:
                trig_idx = 2
            counts_idx = 1  # raw value
            if sum(self.eacc_SKM) > 0:
                counts_idx = 2
            children = packet[13].children
            group = {}
            last_time_bin = 0
            self.num_time_bins = 0
            for i in range(0, num_structures):
                offset = i * 22
                time = children[offset][1] * 0.1 + T0
                if self.start_unix is None:
                    self.start_unix=time

                self.end_unix=time
                if time != last_time_bin:
                    self.num_time_bins += 1
                    last_time_bin = time

                rcr = children[offset + 1][1]
                pixel_mask = [
                    e[1] for e in children[offset + 2][3] if 'NIXG' not in e[0]
                ]
                # exclude NIXG parameters

                detector_mask = children[offset + 3][1]
                integrations = children[offset + 4][1]
                triggers = [children[k + 5][trig_idx] for k in range(0, 16)]
                # id 6-22, 16 trigger accumulators
                num_samples = children[21 + offset][1]
                samples = children[21 + offset][3]
                energies = []
                pixel_indexes = self.get_spectrum_pixel_indexes(
                    detector_mask, pixel_mask)
                for j in range(0, num_samples):
                    k = j * 4
                    E1_low = samples[k + 1][1]
                    E2_high = samples[k + 2][1]
                    pixel_counts = []
                    for idx, e in enumerate(samples[k + 3][3]):
                        pixel_counts.append(e[counts_idx])
                        self.pixel_total_counts[
                            pixel_indexes[idx]] += e[counts_idx]
                    energies.append(
                        [E1_low, E2_high,
                         sum(pixel_counts), pixel_counts])
                group = {
                    'time': time,
                    'rcr': rcr,
                    'pixel_mask': pixel_mask,
                    'detector_mask': detector_mask,
                    'triggers': triggers,
                    'integrations': integrations,
                    'pixel_indexes': pixel_indexes,
                    'energies': energies,
                }
                self.groups.append(group)
        return self.format_report()


class StixBulkL3Analyzer(object):
    def __init__(self):
        self.time = []
        self.rcr = []
        self.dt = []
        self.pixel_mask = []
        self.detector_mask = []
        self.triggers = []
        # self.T0=[]
        self.eacc_SKM = ()
        self.trig_SKM = ()
        self.request_id = -1
        self.packet_unix = 0
        self.groups = []
        self.start_unix=None
        self.end_unix=None

    def format_report(self):
        report = {
            'request_id': self.request_id,
            'packet_unix': self.packet_unix,
            'groups': self.groups,
            'trig_skm': self.trig_SKM,
            'eacc_skm': self.eacc_SKM,
            'start_unix':self.start_unix,
            'end_unix':self.end_unix
        }
        return report

    def merge(self, cursor):
        for pkt in cursor:
            packet = sdt.Packet(pkt)
            self.request_id = packet[3].raw
            self.packet_unix = packet['unix_time']
            T0 = stix_datetime.scet2unix(packet[12].raw)
            num_structures = packet[13].raw
            self.eacc_SKM = (packet[5].raw, packet[6].raw, packet[7].raw)
            self.trig_SKM = (packet[9].raw, packet[10].raw, packet[11].raw)

            trig_idx = 1
            # uncompressed counts
            if sum(self.trig_SKM) > 0:
                trig_idx = 2
            counts_idx = 1  # raw value
            if sum(self.eacc_SKM) > 0:
                counts_idx = 2
            children = packet[13].children
            group = {}
            for i in range(0, num_structures):
                offset = i * 31
                time = children[offset][1] * 0.1 + T0
                if self.start_unix is None:
                    self.start_unix=time
                self.end_unix=time
                rcr = children[offset + 1][1]
                integrations = children[offset + 2][1]
                pixel_mask = [
                    children[offset + 4][1], children[offset + 6][1],
                    children[offset + 8][1]
                ]
                detector_mask = children[offset + 13][1]
                triggers = [children[k + 14][trig_idx] for k in range(0, 16)]
                # id 6-22, 16 trigger accumulators
                num_samples = children[30 + offset][1]
                samples = children[30 + offset][3]
                subgroups = []
                for j in range(0, num_samples):
                    k = j * 5
                    E1_low = samples[k + 1][1]
                    E2_high = samples[k + 2][1]
                    flux = samples[k + 3][1]
                    num_vis = samples[k + 4][1]
                    visiblity_root = samples[k + 4][3]
                    visibilities = [
                        (visiblity_root[m][1], visiblity_root[m + 1][1],
                         visiblity_root[m + 2][1]) for m in range(0, num_vis)
                    ]
                    subgroups.append([E1_low, E2_high, flux, visibilities])
                group = {
                    'time': time,
                    'rcr': rcr,
                    'pixel_mask': pixel_mask,
                    'detector_mask': detector_mask,
                    'triggers': triggers,
                    'integrations': integrations,
                    'subgroups': subgroups,
                }

                self.groups.append(group)
        return self.format_report()


class StixBulkL4Analyzer(object):
    # spectrogram
    def __init__(self):
        self.time = []
        self.rcr = []
        self.dt = []
        self.pixel_mask = []
        self.detector_mask = []
        self.triggers = []
        self.start_time = 0
        self.eacc_SKM = ()
        self.trig_SKM = ()
        self.request_id = -1
        self.packet_unix = 0
        self.groups = []
        self.lightcurves = []
        self.groups = []
        self.num_time_bins = 0
        self.start_unix=None
        self.end_unix=None

    def format_report(self):
        report = {
            'request_id': self.request_id,
            'packet_unix': self.packet_unix,
            'groups': self.groups,
            'num_time_bins': self.num_time_bins,
            'start_unix':self.start_unix,
            'end_unix':self.end_unix
        }
        return report

    def merge(self, cursor):
        last_timestamp=None
        for pkt in cursor:
            packet = sdt.Packet(pkt)
            self.request_id = packet[3].raw
            self.packet_unix = packet['unix_time']
            T0 = stix_datetime.scet2unix(packet[12].raw)
            num_structures = packet[13].raw
            self.eacc_SKM = (packet[5].raw, packet[6].raw, packet[7].raw)
            self.trig_SKM = (packet[9].raw, packet[10].raw, packet[11].raw)

            trig_idx = 1
            if sum(self.trig_SKM) > 0:
                trig_idx = 2
            counts_idx = 1  # raw value
            if sum(self.eacc_SKM) > 0:
                counts_idx = 2
            children = packet[13].children
            group = {}
            if last_timestamp is None:
                last_timestamp = T0
            self.num_time_bins = 0


            for i in range(0, num_structures):
                offset = i * 10
                pixel_mask = children[offset + 1][1]
                detector_mask = children[offset + 2][1]
                rcr = children[offset + 3][1]
                E1 = children[offset + 5][1]
                E2 = children[offset + 6][1]
                Eunit = children[offset + 7][1]

                num_samples = children[8 + offset][1]
                samples = children[8 + offset][3]
                subgroups = []

                for j in range(0, num_samples):
                    k = j * 3
                    timestamp = samples[k + 0][1] * 0.1 + T0
                    if self.start_unix is None:
                        self.start_unix=timestamp
                    self.end_unix=timestamp
                    dT = timestamp - last_timestamp

                    if dT > 0:
                        self.num_time_bins += 1

                    last_timestamp = timestamp
                    triggers = samples[k + 1][trig_idx]
                    num_energies = samples[k + 2][1]
                    energy_children = samples[k + 2][3]
                    lcs = [m[counts_idx] for m in energy_children]
                    subgroups.append((timestamp, triggers, lcs, dT))

                group = {
                    'rcr': rcr,
                    'pixel_mask': pixel_mask,
                    'detector_mask': detector_mask,
                    'E1': E1,
                    'E2': E2,
                    'Eunit': Eunit,
                    'subgroups': subgroups,
                }
                self.groups.append(group)
        return self.format_report()


class StixBulkAspectAnalyzer(object):
    def __init__(self):
        self.start_unix=None
        self.end_unix=None
        pass

    def merge(self, cursor):
        readouts = [[], [], [], []]
        read_time = []
        packet_utc = ''
        start_time = 0
        for pkt in cursor:
            packet = sdt.Packet(pkt)
            packet_utc = packet['UTC']
            T0 = stix_datetime.scet2unix(
                packet[1].raw) + packet[2].raw / 65536.
            if self.start_unix is None:
                self.start_unix=T0

            if start_time == 0:
                start_time = T0
            dt = packet[4].raw / 64 * (1 / 64.)  # has to be fixed
            children = packet[5].children
            for i, param in enumerate(children):
                readouts[i % 4].append(param[1])
                if i % 4 == 0:
                    ut=dt * int(i / 4) + T0
                    read_time.append(ut)
                    self.end_unix=ut

        return {
            'packet_utc': packet_utc,
            'readouts': readouts,
            'read_time': read_time,
            'start_time': start_time,
            'start_unix':self.start_unix,
            'end_unix':self.end_unix

        }


def process(file_id):
    try:
        merge(file_id)
    except Exception as e:
        print(str(e))
        logger.error(str(e))


def merge(file_id, remove_existing=True):
    collection = mdb.get_collection_bsd()
    bsd_cursor = collection.find({'run_id': file_id}).sort('_id', 1)
    for doc in bsd_cursor:
        spid = int(doc['SPID'])
        logger.info(f'processing bsd id: {doc["_id"]}, spid:{spid}')
        cursor = mdb.get_packets_of_bsd_request(doc['_id'], header_only=False)
        if not cursor:
            continue
        result = None
        if spid == 54125:
            analyzer = StixBulkAspectAnalyzer()
            result = analyzer.merge(cursor)
        elif spid == 54114:
            analyzer = StixBulkL0Analyzer()
            result = analyzer.merge(cursor)
        elif spid in [54115, 54116]:
            analyzer = StixBulkL1L2Analyzer()
            result = analyzer.merge(cursor)

        elif spid == 54117:
            analyzer = StixBulkL3Analyzer()
            result = analyzer.merge(cursor)
        elif spid == 54143:
            analyzer = StixBulkL4Analyzer()
            result = analyzer.merge(cursor)
        if result:
            date_str=datetime.now().strftime("%y%m%d%H")
            existing_fname=doc.get('level1','')
            if existing_fname:
                old_file= os.path.join(level1_products_path,existing_fname)
                try:
                    os.remove(old_file)
                except Exception:
                    pass
            json_filename = os.path.join(level1_products_path,
                                         f'L1_{doc["_id"]}_{date_str}.json')
            start_unix=result.get('start_unix',0)
            end_unix=result.get('end_unix',0)
            with open(json_filename, 'w') as outfile:
                json.dump(result, outfile)
            collection.update_one({'_id': doc['_id']}, {'$set':{'level1': json_filename, 
                'start_unix':start_unix, 'end_unix':end_unix}})


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('sci_packet_merger  run_id')
        print('sci_packet_merger  run_id_start id_end')
    elif len(sys.argv)==2:
        process(int(sys.argv[1]))
    else:
        for i in range(int(sys.argv[1]),int(sys.argv[2])+1):
            print('process:',i)
            process(i)

