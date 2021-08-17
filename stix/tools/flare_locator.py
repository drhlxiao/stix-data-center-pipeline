# fit flare location manually, for an given data requests

import sys
import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
bsd_db= db['bsd']
packet_db= db['packets']
def get_packets(unique_id):
    unique_id=int(unique_id)
    cursor=bsd_db.find({'unique_id': unique_id})
    for row in cursor:
        packet_ids=row['packet_ids']
        print('Number of packets:', len(packet_ids))
        return packet_db.find({'_id':{'$in':packet_ids}}).sort('_id',1)
    return []

                
ebins= [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150, 1e3];  

ebin_map={
        '8-12':[4,9],
        '8-10':[4,7]
        }


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


    def format_report(self):
        report = {
            'request_id': self.request_id,
            'packet_unix': self.packet_unix,
            'groups': self.groups,
            'trig_skm': self.trig_SKM,
            'eacc_skm': self.eacc_SKM,
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
    def process(self,unique_id, start_utc, end_utc,range_name):
        emin,emax=ebin_map[range_name]
        start_unix=stix_datetime.utc2unix(start_utc)
        end_unix=stix_datetime.utc2unix(end_utc)
        packets=get_packets(unique_id,start_unix, end_unix, emin,emax)
        return merge(packets,start_unix, end_unix, emin, emax) 
    def merge(self, cursor, start_unix, end_unix, emin, emax):
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
            last_time_bin = 0
            self.num_time_bins = 0
            for i in range(0, num_structures):
                offset = i * 22
                time = children[offset][1] * 0.1 + T0

                if time<start_unix or time>end_unix:
                    continue


                if time != last_time_bin:
                    self.num_time_bins += 1
                    last_time_bin = time

                rcr = children[offset + 1][1]
                pixel_mask = [
                    e[1] for e in children[offset + 2][3] if 'NIXG' not in e[0]
                ]
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
                    if E1_low<emin or E2_high>emax:
                        continue


                    pixel_counts = []
                    for idx, e in enumerate(samples[k + 3][3]):
                        pixel_counts.append(e[counts_idx])
                        self.pixel_total_counts[
                            pixel_indexes[idx]] += e[counts_idx]
                #    energies.append(
                #        [E1_low, E2_high,
                #         sum(pixel_counts), pixel_counts])
                #group = {
                #    'time': time,
                #    'rcr': rcr,
                #    'pixel_mask': pixel_mask,
                #    'detector_mask': detector_mask,
                #    'triggers': triggers,
                #    'integrations': integrations,
                #    'pixel_indexes': pixel_indexes,
                #    'energies': energies,
                #}
                #self.groups.append(group)
        return {'emin':emin,
                'emax':emax,
                'start': stix_datetime.unix2utc(start_unix),
                'end': stix_datetime.unix2utc(end_unix),
                'counts':self.pixel_total_counts,
                'cfl':self.pixel_total_counts[8*12, 8*12+12]
                }

        #return self.format_report()




p=StixBulkL1L2Analyzer()
counts=p.process(unique_id, start_utc, end_utc, e_range)
