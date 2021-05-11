import sys
import os
sys.path.append(os.path.abspath('../../'))
from stix.core import stix_datatypes as sdt
from stix.spice import stix_datetime
from stix.core import stix_logger
from pprint import pprint
logger = stix_logger.get_logger()
import numpy as np


class StixCalibrationReportAnalyzer(object):
    """
      Capture calibration reports and fill calibration information into a MongoDB collection

    """

    def __init__(self):
        self.report = None

        self.total_counts = np.zeros((8,12*32), dtype=np.int32)
        self.packet_indexing= np.zeros((8,12*32), dtype=np.int32)

        #self.db_collection = None
        self.current_calibration_run_id = 0
        #self.db_collection = collection
        #self.background_spectra=[[] for x in range(0,12*32)]
        #try:
        #    self.current_calibration_run_id = self.db_collection.find().sort(
        #        '_id', -1).limit(1)[0]['_id'] + 1
        #except IndexError:
        #    self.current_calibration_run_id = 0

    def capture(self, run_id, packet_id, pkt):
        #if not self.db_collection:
        #    return
        packet = sdt.Packet(pkt)
        if packet.SPID != 54124:
            return

        #detector_mask=packet[10].raw
        #pixel_mask=packet[12].raw
        sbspec_mask=packet[13].raw
        #for i in range(0,8):
        #    if sbspec_mask & (1<<i) !=0:
        #        sbspec_status[i]=True
        #print('number of subspec:', sbspec_status)


        subspec_formats=[(packet[i*4+15].raw, packet[i*4+16].raw, packet[i*4+17].raw) for i in range(0,8)]

        detector_ids = packet.get('NIX00159/NIXD0155')[0]
        pixel_ids = packet.get('NIX00159/NIXD0156')[0]
        subspc_ids=packet.get('NIX00159/NIXD0157')[0]
        sbspec_spectra = packet.get('NIX00159/NIX00146/*.eng')[0]
        #print('dector ids')
        #print(len(detector_ids))
        #print('subspectra length')
        #print(len(sbspec_spectra))
        #return

        #return
        for i in range(0, len(detector_ids)):
            sbspec_id=subspc_ids[i]
            detector_id=detector_ids[i]
            pixel_id=pixel_ids[i]
            sbspec_spectrum=sbspec_spectra[i]
            pprint(sbspec_spectrum)
            return
            num_counts=sum(sbspec_spectrum)
            self.total_counts[sbspec_id][detector_id*12+pixel_id] = num_counts
            #self.total_counts.append([
            #    packet_id, detector_ids[ispec], pixel_ids[ispec], sbspec_id, total_counts])
            #self.background_sbspec_spectrum[detector_ids[ispec]*12+pixel_ids[ispec]]=spectrum

        if packet['seg_flag'] in [1, 3]:
            #first or standalone packet
            self.report = {
                'run_id': run_id,
                'packet_ids': [packet_id],
                'sbspec_mask':sbspec_mask,
                'sbspec_formats':subspec_formats,
                'header_unix_time': packet['unix_time'],
                '_id': self.current_calibration_run_id,
            }
        else:
            if not self.report:
                logger.warning('The first calibration report is missing!')
            else:
                self.report['packet_ids'].append(packet_id)

        if packet['seg_flag'] in [2, 3]:
            #last packet or standalone 
            if not self.report:
                logger.warning(
                        'A calibration run (last packet ID:{}) is not recorded due to missing the first packet!'.format(packet_id)
                )

            param_dict = packet.children_as_dict()
            self.report['duration'] = param_dict['NIX00122'][0].raw
            scet = param_dict['NIX00445'][0].raw
            self.report['SCET'] = scet
            self.report['start_unix_time'] = stix_datetime.scet2unix(scet)
            self.report['auxiliary'] = packet['parameters'][0:14]
            #Don't copy repeaters 
            self.report['total_counts'] = self.total_counts
            #self.report['spectra'] = self.background_spectra
            #spectra
            #self.db_collection.insert_one(self.report)
            pprint(self.report)
            return

            self.current_calibration_run_id += 1
            self.report = None
            self.total_counts = np.zeros((8,12*32))
            #self.background_spectra=[[] for x in range(0,12*32)]



class Plugin:
    """ don't modify here """

    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        print("Plugin  loaded ...")

    def run(self):
        # your code goes here
        print('current row')
        print(self.current_row)
        #if len(self.packets) > 1:
        #    pprint.pprint(self.packets[self.current_row])
        cal=StixCalibrationReportAnalyzer()
        for i,packet in enumerate(self.packets):
            cal.capture(0, i, packet)


