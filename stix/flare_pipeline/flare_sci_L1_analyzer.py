#!/usr/bin/python3
# author: Hualin Xiao
# pre-process science data, merge bulk science data packets and write merged data to json files
# so that web client side could load the data quickly
import sys
sys.path.append('.')
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

class SciL1Analyzer(object):
    def __init__(self):
        self.dt = []
        self.detector_mask = []
        self.triggers = []
        self.hits = []
        self.eacc_SKM = ()
        self.trig_SKM = ()
        self.start_unix=None
        self.end_unix=None
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
    def compute(self, cursor, min_unix_time=None, max_unix_time=None): 
        #calculate spectrogram
        detector_counts=np.zeros((32,32)) #32 energy spectrum and 32 energy bins
        energy_map= []
        for pkt in cursor:
            packet = sdt.Packet(pkt)
            self.request_id = packet[3].raw
            self.packet_unix = packet['unix_time']
            T0 = stix_datetime.scet2unix(packet[12].raw)
            num_structures = packet[13].raw
            self.eacc_SKM = (packet[5].raw, packet[6].raw, packet[7].raw)
            self.trig_SKM = (packet[9].raw, packet[10].raw, packet[11].raw)
            counts_idx = 1  # raw value
            if sum(self.eacc_SKM) > 0:
                counts_idx = 2
            children = packet[13].children
            group = {}
            last_time_bin = 0
            self.num_time_bins = 0
            emin=32
            emax=0
            for i in range(0, num_structures):
                offset = i * 22
                time = children[offset][1] * 0.1 + T0

                if min_unix_time is not None and max_unix_time  is not None:
                    if time>max_unix_time or time<min_unix_time:
                        continue 
                if self.start_unix is None:
                    self.start_unix=time
                self.end_unix=time
                rcr = children[offset + 1][1]
                pixel_mask = [
                    e[1] for e in children[offset + 2][3] if 'NIXG' not in e[0]
                ]
                detector_mask = children[offset + 3][1]
                integrations = children[offset + 4][1]
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
                    ebin_name=(E1_low, E2_high)
                    if ebin_name not in energy_map:
                        energy_map.append(ebin_name)
                    i_energy=energy_map.index(ebin_name)

                    for idx, e in enumerate(samples[k + 3][3]):
                        pixel_counts.append(e[counts_idx])
                        detector_id=int(pixel_indexes[idx]/12)
                        detector_counts[detector_id][i_energy]+= e[counts_idx]
        duration=self.end_unix-self.start_unix
        detector_count_rates=detector_counts/duration
        
        int(np.array(energy_map).max())
        return {
                'rates':detector_count_rates.tolist(),
                'total_counts':detector_counts.tolist(),
                'description':' 32 detectors x 32 spectra, rates (cnts/s)',
                'emin':int(np.array(energy_map).min()),
                'emax':int(np.array(energy_map).max()),
                'energy_map':energy_map,
                'start_unix':self.start_unix,
                'end_unix':self.end_unix,
                'duration':self.end_unix-self.start_unix,
                }


def process_L1_BSD_in_file(file_id):
    collection = mdb.get_collection_bsd()
    bsd_cursor = collection.find({'run_id': file_id, 'SPID':54115}).sort('_id', 1)
    start_unix=None
    end_unix=None
    flare_location=None
    for doc in bsd_cursor:
        db_id=doc['_id']
        #print(doc)
        try:
            bsd_start_unix=doc['request_form']['start_unix']
            bsd_end_unix=doc['request_form']['end_unix']
        except KeyError or TypeError:
            print('Skipped, failed to find data request form for BSD #',db_id)
            print('run tools/attach_data_req*.py')
            continue
        bsd_duration=bsd_end_unix-bsd_start_unix
        flares=list(mdb.search_flares_by_tw(bsd_start_unix, bsd_duration, threshold=0))
        packet_cursor= mdb.get_packets_of_bsd_request(db_id, header_only=False)
        analyzer = SciL1Analyzer()
        result = analyzer.compute(packet_cursor, bsd_start_unix, bsd_end_unix)
        start_unix=result['start_unix']
        end_unix=result['end_unix']
        #data real start time and end time
        duration=end_unix-start_unix
        is_background=True if flares else False
        flare_info=None
        is_background=True

        if flares:
            is_background=False
            flare_doc=flares[0]
            flare_info={
                    'flare_id': flares[0]['flare_id'],
                    'flare_entry_id':flares[0]['_id'],
                    'peak_utc':flares[0]['peak_utc']
                    }
            bkg_cursor= collection.find({'db_id': {'$lt':db_id}, 'SPID':54115}).sort('_id', -1).limit(1)
            flare_level1_result=analyzer_flare_level1(result,flare_doc)

        print('Processing BSD #', db_id)


        collection.update_one({'_id': doc['_id']}, {'$set':
            {'is_background': is_background, 
            'detector_statistics':result, 
            'flare_info':flare_info,
            'start_unix':start_unix,
            'end_unix':end_unix,
            'flare_entry_id':flare_entry_id,
            'flare_level1_result': flare_level1_result,
                }})


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('flare_location_solver run_id')
        print('flare_location_solver run_id_start id_end')
    elif len(sys.argv)==2:
        process_L1_BSD_in_file(int(sys.argv[1]))
    else:
        for i in range(int(sys.argv[1]),int(sys.argv[2])+1):
            process_L1_BSD_in_file(i)

