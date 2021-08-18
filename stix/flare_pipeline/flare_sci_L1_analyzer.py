#!/usr/bin/python3
# author: Hualin Xiao
# Compute flare location, do background subtraction
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
DET_ID_START=20
DET_ID_END=32
DET_SENSITVE_AREA=80.96
GRID_OPEN_AREA_RATIO=np.array([0.292312,
		0.272487,
		0.252282,
		0.257189,
		0.281502,
		0.27254,
		0.292103,
		0.260207,
		0.335213,
		0.305518,
		0.335412,
		0.254877,
		0.264532,
		0.255886,
		0.304781,
		0.335199,
		0.30505,
		0.258,
		0.257187,
		0.252324,
		0.281503,
		0.260178,
		0.280854,
		0.257033,
		0.264449,
		0.261743,
		0.292299,
		0.272491,
		0.264514,
		0.254682])

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
        self.time_integrated_counts=None
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
    def get_flare_spectra_and_location(self,flare_bsd_doc,flare_doc):
        #bsd_doc, bulk science doc for this run
        results={}
        bsd_id=flare_bsd_doc['_id']

        if not self.time_integrated_counts:
            print('run compute_detector_counts first')
            return None
        bkg_cursor= list(collection.find({'_id': {'$lt':bsd_id},
            'request_form.emin':flare_bsd_doc['request_form']['emin'],
            'request_form.emax':{'$gte': flare_bsd_doc['request_form']['emax']},
            'request_form.detector_mask':flare_bsd_doc['request_form']['detector_mask'],
            'request_form.eunit':flare_bsd_doc['request_form']['eunit'],
            'is_background':True,
            'SPID':54115}).sort('_id', -1).limit(1))
        #background should have the same configuration
        if not bkg_cursor:
            print("Can not find background bsd data for bulk science data #",bsd_id)
            return None
        bkg_doc=bkg_cursor[0]
        bkg_stat=bkg_doc['time_integrated_counts']
        sig_counts=np.array(self.time_integrated_counts['detector_total_counts'])


        sig_duration=self.time_integrated_counts['duration']
        sig_energy_map=self.time_integrated_counts['energy_map']
        min_bin=0
        max_bin=8
        for ibin, ebin in enumerate(sig_energy_map):
            if ebin[0]>=1: #4 keV
                min_bin=ibin
                break
        for ibin, ebin in enumerate(sig_energy_map):
            if ebin[1]>=8: #12 keV
                max_bin=ibin
                break
        results['bkg_bsd_id']=bkg_doc['_id']

        bkg_rates=np.array(bkg_stat['detector_rates'])
        bkg_rate_error=np.array(bkg_stat['detector_rate_error'])

        bkg_cfl_rates=np.array(bkg_stat['pixel_rates'][8*12:9*12,:])
        bkg_cfl_rate_error=np.array(bkg_stat['pixel_rate_error'][8*12:9*12,:])
        cfl_pixel_counts=np.array(self.time_integrated_counts['pixel_total_counts'][8*12:9*12,:])

        cfl_pixel_bkgsub_counts=cfl_pixel_counts - sig_duration * bkg_cfl_rates
        cfl_pixel_bkgsub_count_error=np.sqrt(cfl_pixel_counts  + sig_duration*sig_duration*bkg_cfl_rate_error **2)

        bkg_rates=bkg_rates[:,0:len(sig_energy_map)]  #make the same size
        bkg_rate_error=bkg_rate_error[:,0:len(sig_energy_map)]  #make the same size
        sig_bkgsub_counts=sig_counts - sig_duration*bkg_rates
        sig_bkgsub_count_error=np.sqrt(sig_counts+sig_duration*sig_duration*bkg_rate_error**2)

        results['detector_spectra_subtracted_bkg']=sig_bkgsub_counts
        results['detector_spectra_subtracted_bkg_error']=sig_bkgsub_count_error
        

        det_mean_counts=np.mean(np.sum(sig_bkgsub_counts[DET_ID_START:DET_ID_END, min_bin:max_bin], axis=1)/GRID_OPEN_AREA_RATIO[DET_ID_END:DET_SENSITVE_AREA])
        det_mean_counts_error=np.sum(sig_bkgsub_count_error[DET_ID_START:DET_ID_END, min_bin:max_bin],axis=1)**2/GRID_OPEN_AREA_RATIO[DET_ID_END:DET_SENSITVE_AREA])

        
        results['cfl_pixel_spectra_subtracted_bkg_4_12_keV']=cfl_pixel_bkgsub_counts[:, min_bin:max_bin]
        results['cfl_pixel_spectra_subtracted_bkg_error_4_12_keV']=cfl_pixel_bkgsub_count_error[:, min_bin:max_bin]
        detector_4_12_keV_counts=np.sum(sig_bkgsub_counts[:, min_bin:max_bin],axis=1)


        return results
        
        
    def compute_detector_counts(self, cursor, min_unix_time=None, max_unix_time=None): 
        #calculate spectrogram
        detector_counts=np.zeros((32,32)) #32 energy spectrum and 32 energy bins
        pixel_counts=np.zeros((32*12,32)) #32 energy spectrum and 32 energy bins
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
                        pixel_id=pixel_indexes[idx]
                        detector_counts[detector_id][i_energy]+= e[counts_idx]
                        pixel_counts[pixel_id][i_energy]+= e[counts_idx]
        duration=self.end_unix-self.start_unix
        num_energy_bin=len(energy_map)
        detector_counts=detector_counts[:,0:num_energy_bin]
        detector_count_rates=detector_counts/duration
        detector_count_rate_errors=np.sqrt(detector_counts)/duration
        
        pixel_counts=pixel_counts[:,0:num_energy_bin]
        pixel_count_rates=pixel_counts/duration
        pixel_count_rate_errors=np.sqrt(pixel_counts)/duration

        self.time_integrated_counts={
                'detector_rates':detector_count_rates.tolist(),
                'detector_rate_error':detector_count_rate_errors.tolist(),
                'detector_total_counts':detector_counts.tolist(),
                'detector_description':f'32 detectors x {num_energy_bin} spectra, rates (cnts/s)',

                'pixel_rates':pixel_count_rates.tolist(),
                'pixel_rate_error':pixel_count_rate_errors.tolist(),
                'pixel_total_counts':pixel_counts.tolist(),
                'pixel_description':f'384 pixels x {num_energy_bin} spectra, rates (cnts/s)',

                'emin':int(np.array(energy_map).min()),
                'emax':int(np.array(energy_map).max()),
                'energy_map':energy_map,
                'start_unix':self.start_unix, #flare start time if flare data
                'end_unix':self.end_unix,
                'duration':self.end_unix-self.start_unix,
                }
        return self.time_integrated_counts


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
            bsd_start_unix=doc['start_unix']
            bsd_end_unix=doc['end_unix']
        except KeyError or TypeError:
            print('Skipped, failed to find data request form for BSD #',db_id)
            print('run tools/attach_data_req*.py')
            continue
        bsd_duration=bsd_end_unix-bsd_start_unix
        flares=list(mdb.search_flares_by_tw(bsd_start_unix, bsd_duration, threshold=0))
        packet_cursor= mdb.get_packets_of_bsd_request(db_id, header_only=False)
        analyzer = SciL1Analyzer()
        result = analyzer.compute_detector_counts(packet_cursor, bsd_start_unix, bsd_end_unix)
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
            flare_level1_result=analyzer.get_flare_spectra_and_location(bsd_doc, flare_doc)

        print('Processing BSD #', db_id)
        collection.update_one({'_id': doc['_id']}, {'$set':
            {'is_background': is_background, 
            'time_integrated_counts':result, 
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

