#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    Create metadata from science reports. The information from science packets is written into mongodb
    @Author: Hualin Xiao
    @Date: Nov. 2019
"""
import sys
import numpy as np
from stix.core import stix_datatypes as sdt
from stix.spice import stix_datetime
from stix.core import stix_logger
logger = stix_logger.get_logger()
DATA_REQUEST_REPORT_SPIDS = [54114, 54115, 54116, 54117, 54143, 54125]
DATA_REQUEST_REPORT_NAME = {
    54114: 'L0',
    54115: 'L1',
    54116: 'L2',
    54117: 'L3',
    54143: 'L4',
    54125: 'ASP'
}
QL_REPORT_SPIDS = [54118, 54119, 54121, 54120, 54122]
EBINS = [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150, 'Emax']


class StixScienceReportAnalyzer(object):
    def __init__(self, db):
        self.calibration_analyzer = StixCalibrationReportAnalyzer(
            db['calibration_runs'])
        self.ql_analyzer = StixQuickLookReportAnalyzer(db)
        self.user_request_analyzer = StixUserDataRequestReportAnalyzer(db)

    def start(self, run_id, packet_id, packet):
        self.calibration_analyzer.capture(run_id, packet_id, packet)
        self.ql_analyzer.capture(run_id, packet_id, packet)
        self.user_request_analyzer.capture(run_id, packet_id, packet)

    def get_calibration_run_ids(self):
        return self.calibration_analyzer.get_calibration_run_ids()


class StixCalibrationReportAnalyzer(object):
    """
      Capture calibration reports and fill calibration information into MongoDB 

    """
    def __init__(self, collection):
        self.report = None
        self.sbspec_counts = np.zeros((8, 12 * 32), dtype=np.int32)
        self.packet_index = np.zeros((8, 12 * 32), dtype=np.int32)
        self.db_collection = None
        self.current_calibration_run_id = 0
        self.db_collection = collection
        self.calibration_run_ids = []
        self.spectra = []
        #self.background_spectra=[[] for x in range(0,12*32)]
        try:
            self.current_calibration_run_id = self.db_collection.find().sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            self.current_calibration_run_id = 0

    def get_calibration_run_ids(self):
        return self.calibration_run_ids

    def capture(self, run_id, packet_id, pkt):
        if not self.db_collection:
            return
        packet = sdt.Packet(pkt)
        if not packet['parameters']:
            return
        if packet.SPID != 54124:
            return

        detector_mask = packet.get_one('NIX00407')[1]  #raw value
        pixel_mask = packet.get_one('NIXD0407')[1]  #raw value

        detector_ids = packet.get('NIX00159/NIXD0155')[0]
        pixel_ids = packet.get('NIX00159/NIXD0156')[0]
        subspc_ids = packet.get('NIX00159/NIXD0157')[0]
        sbspec_spectra = packet.get('NIX00159/NIX00146/*.eng')[0]

        start_index = packet.index('NIXD0129')

        sbspec_formats = [(packet[i * 4 + start_index].raw,
                           packet[i * 4 + start_index + 1].raw,
                           packet[i * 4 + start_index + 2].raw)
                          for i in range(0, 8)]

        for i in range(0, len(detector_ids)):
            sbspec_id = subspc_ids[i]
            detector_id = detector_ids[i]
            pixel_id = pixel_ids[i]
            sbspec_spectrum = sbspec_spectra[i]
            try:
                self.spectra.append(
                    (detector_id, pixel_id, sbspec_id,
                     sbspec_formats[sbspec_id][2],
                     int(sbspec_formats[sbspec_id][1]) + 1, sbspec_spectrum))
            except Exception as e:
                logger.warning(
                    'Error occurred when formatting spectra: {}'.format(
                        str(e)))

            num_counts = 0
            try:
                num_counts = sum(sbspec_spectrum)
            except TypeError:
                logger.warning(
                    "Counts not decompressed. File id: {}, Packet id:{}".
                    format(run_id, packet_id))
            self.sbspec_counts[sbspec_id][detector_id * 12 +
                                          pixel_id] = num_counts
            self.packet_index[sbspec_id][detector_id * 12 +
                                         pixel_id] = packet_id

        if packet['seg_flag'] in [1, 3]:
            sbspec_mask = packet[13].raw
            sbspec_status = [False for i in range(0, 8)]
            for i in range(0, 8):
                if sbspec_mask & (1 << i) != 0:
                    sbspec_status[i] = True

            self.report = {
                'run_id': run_id,
                'packet_ids': [packet_id],
                'sbspec_formats': sbspec_formats,
                'sbspec_status': sbspec_status,
                'sbspec_mask': sbspec_mask,
                'header_unix_time': packet['unix_time'],
                '_id': self.current_calibration_run_id,
            }
        else:
            #continuation packet
            if not self.report:
                logger.warning('The first calibration report is missing!')
            else:
                self.report['packet_ids'].append(packet_id)

        if packet['seg_flag'] in [2, 3]:
            if not self.report:
                logger.warning(
                    'A calibration run (last packet ID:{}) is'
                    ' not recorded due to missing the first packet!'.format(
                        packet_id))
                return

            param_dict = packet.children_as_dict()
            self.report['duration'] = param_dict['NIX00122'][0].raw
            scet = param_dict['NIX00445'][0].raw
            self.report['SCET'] = scet
            self.report['start_unix_time'] = stix_datetime.scet2unix(scet)
            spectra_index = packet.index('NIX00159')
            self.report['auxiliary'] = packet['parameters'][0:spectra_index]
            #Don't copy repeaters
            self.report['sbspec_counts_sum'] = np.sum(self.sbspec_counts,
                                                      axis=1).tolist()
            self.report['counts'] = self.sbspec_counts.tolist()
            self.report['spectra'] = self.spectra
            self.report['packet_index'] = self.packet_index.tolist()
            #packet id, used by web applications
            #spectra
            self.db_collection.insert_one(self.report)
            self.calibration_run_ids.append(self.current_calibration_run_id)

            #reset report
            self.current_calibration_run_id += 1
            self.report = None
            self.spectra = []
            self.counts = np.zeros((8, 12 * 32), dtype=np.int32)
            self.packet_index = np.zeros((8, 12 * 32), dtype=np.int32)


class StixQuickLookReportAnalyzer(object):
    """
    capture quicklook reports and fill packet information into a MongoDB collection

    """
    def __init__(self, db):
        self.stix_db=db
        self.qlspec_db=db.get_collection('ql_spectra')
        self.qllc_db=db.get_collection('ql_lightcurves')
        self.qlloc_db=db.get_collection('ql_flarelocations')
        self.ql_db_collection = db.get_collection('quick_look')
        self.report = None
        try:
            self.current_report_id = self.ql_db_collection.find().sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            self.current_report_id = 0
    def get_qllc_energy_bins(self, emask):
        ebins=[]
        for i in range(32):
            if emask & (1 << i) != 0:
                ebins.append(i);
        names=[]
        sci_edges=[]
        for i in range(len(ebins) - 1):
            begin = ebins[i]
            end = ebins[i + 1]
            sci_edges.append([begin,end])
            if end == 32:
                names.append(f'{EBINS[begin]} keV –⁠ Emax')
            elif end < 32:
                names.append(f'{EBINS[begin]}  – {EBINS[end]} keV')
            else:
                names.append('')
        return {'names':names, 'sci_bin_edges':sci_edges}

    def write_qllc_to_db(self, packet, run_id, packet_id):
        lightcurves = {}
        unix_time = []
        energy_bins={}
        last_time=0
        if not isinstance(packet, sdt.Packet):
            packet = sdt.Packet(packet)
        if not packet.isa(54118) or not packet.is_valid():
            return

        scet_coarse = packet[1].raw
        scet_fine = packet[2].raw
        start_scet = scet_coarse + scet_fine / 65536.
        int_duration = (packet[3].raw + 1) * 0.1
        detector_mask = packet[4].raw
        pixel_mask = packet[6].raw
        num_lc = packet[17].raw
        compression_s = packet[8].raw
        compression_k = packet[9].raw
        compression_m = packet[10].raw
        if not energy_bins:
            energy_bin_mask= packet[16].raw
            energy_bins=self.get_qllc_energy_bins(energy_bin_mask)
        num_lc_points = packet.get('NIX00270/NIX00271')[0]
        lc = packet.get('NIX00270/NIX00271/*.eng')[0]
        rcr = packet.get('NIX00275/*.raw')
        trig= packet.get('NIX00273/*.raw')
        UTC = packet['header']['UTC']
        for i in range(len(lc)):
            if i not in lightcurves:
                lightcurves[str(i)] = []
            lightcurves[str(i)].extend(lc[i])
        unix_time.extend([
            stix_datetime.scet2unix(start_scet + x * int_duration)
            for x in range(num_lc_points[0])
        ])

        if not lightcurves:
            return 
        doc={'run_id':run_id, 'packet_id':packet_id, 'time': unix_time, 'lightcurves': lightcurves,'rcr':rcr,'trig':trig,
                'energy_bins':energy_bins,'num':len(unix_time),'start_unix': unix_time[0],'end_unix':unix_time[-1]}
        self.qllc_db.save(doc)

    def write_ql_spec_to_db(self, packet, run_id, packet_id):
        #write ql spectra to ql_spectra database
        if not isinstance(packet, sdt.Packet):
            packet=sdt.Packet(packet)
        if not packet.isa(54120) or not packet.is_valid():
            return

        scet_coarse = packet[1].raw
        scet_fine = packet[2].raw
        start_unix=stix_datetime.scet2unix(scet_coarse, scet_fine)
        tbin= (packet[3].raw + 1)*0.1
        num_samples=packet[14].raw
        samples=packet[14].children
        for i in range(num_samples):
            start_i=35*i
            detector=samples[start_i][1]
            spectra=[samples[1+j+start_i][2] for j in range(32)]
            trig=samples[33+start_i][2]
            t=start_unix+samples[34+start_i][1]*tbin#number of integrations after the first one
            doc={'detector':detector,
                 'spectra':spectra,
                 'triggers':trig,
                 'run_id':run_id,
                 'obs_time': t,
                 'tbin':tbin,
                 'packet_id':packet_id
                }
            self.qlspec_db.save(doc)
    def write_ql_flare_loc_to_db(self, packet, run_id, packet_id):
        #write ql spectra to ql_spectra database
        if not isinstance(packet, sdt.Packet):
            packet=sdt.Packet(packet)
        if not packet.isa(54122) or not packet.is_valid():
            return
        scet_coarse = packet[1].raw
        scet_fine = packet[2].raw
        start_unix=stix_datetime.scet2unix(scet_coarse, scet_fine)
        tbin= (packet[3].raw + 1)*0.1
        num_samples=packet[4].raw
        samples=packet[4].children
        #parameters=['NIX00283','NIX00284','NIXD0060','NIXD0061']
        for i in range(num_samples):
            offset=i*7
            locz=samples[offset+5][1]
            locy=samples[offset+6][1]
            if locz==0 and locy==0:
                continue
            non_thermal_index=samples[offset+3][1]
            thermal_index=samples[offset+4][1]
            
            doc={
                 'locz':locz,
                 'locy':locy,
                 'nonthermal_index':non_thermal_index,
                 'thermal_index':thermal_index,
                 'run_id':run_id,
                 'obs_time': start_unix+i*tbin,
                 'packet_id':packet_id
                }
            self.qlloc_db.save(doc)

    def capture(self, run_id, packet_id, pkt):
        if not self.ql_db_collection:
            return
        packet = sdt.Packet(pkt)

        if packet.SPID not in QL_REPORT_SPIDS:
            return
        if not packet.is_valid():
            return

        start_coarse_time = 0
        start_fine_time = 0
        integrations = 0
        detector_mask = 0
        pixel_mask = 0
        start_unix_time = 0
        duration = 0
        start_coarse_time = packet[1].raw
        start_fine_time = packet[2].raw
        integrations = packet[3].raw
        points = 0

        if packet.SPID == 54118:
            #QL LC
            #detector_mask = packet[4].raw
            #pixel_mask = packet[6].raw
            points = packet.get('NIX00270/NIX00271')[0][0]
            try:
                self.write_qllc_to_db(packet, run_id,packet_id)
            except Exception as e: #in case a broken packet
                raise
                print(str(e))
        elif packet.SPID == 54119:
            points = packet[15].raw
        elif packet.SPID == 54121:
            points = packet[13].raw
        elif packet.SPID == 54120:
            points = packet[14].raw
            try:
                self.write_ql_spec_to_db(packet, run_id,packet_id)
            except TypeError or IndexError: #in case a broken packet
                pass
        elif packet.SPID == 54122:
            points = packet[4].raw
            self.write_ql_flare_loc_to_db(packet, run_id,packet_id)
            #flare flag report

        duration = points * 0.1 * (integrations + 1)
        start_unix_time = stix_datetime.scet2unix(start_coarse_time,
                                                  start_fine_time)
        start_scet = start_coarse_time + start_fine_time / 65536.
        report = {
            '_id': self.current_report_id,
            'run_id': run_id,
            'SPID': packet.SPID,
            'packet_id': packet_id,
            'start_unix_time': start_unix_time,
            'stop_unix_time': start_unix_time + duration,
            'start_scet': start_scet,
            'packet_header_time': packet['unix_time'],
            'integrations': integrations,
            'duration': duration,
            'stop_scet': start_scet + duration,
            #'detector_mask': detector_mask,
            #'pixel_mask': pixel_mask
        }
        self.ql_db_collection.insert_one(report)
        self.current_report_id += 1


class StixUserDataRequestReportAnalyzer(object):
    def __init__(self, stix_db):
        self.db = stix_db['bsd']
        self.db_bsd_forms = stix_db['bsd_req_forms']
        self.last_unique_id = -1
        self.last_request_spid = -1
        self.packet_ids = []
        self.has_first_packet=False
        self.start_time = 0
        self.stop_time = 0
        try:
            self.current_id = self.db.find().sort('_id',
                                                  -1).limit(1)[0]['_id'] + 1
        except IndexError:
            self.current_id = 0

    def capture(self, run_id, packet_id, pkt):
        if not self.db:
            return
        packet = sdt.Packet(pkt)
        if packet.SPID not in DATA_REQUEST_REPORT_SPIDS:
            return
        if not packet['parameters']:
            return
        if packet.SPID == 54125:
            #aspect data
            try:
                self.process_aspect(run_id, packet_id, packet)
            except Exception as e:
                logger.warning(
                    'Can not parse aspect packet #{} in file #{}'.format(
                        packet_id, run_id))
                #print(str(e))
                #print(packet)
        else:
            self.process_bulk_science(run_id, packet_id, packet)

    def process_bulk_science(self, run_id, packet_id, packet):
        try:
            unique_id = packet[3].raw
            start = packet[12].raw
        except Exception as e:
            logger.warning(str(e))
            return
        if packet['seg_flag'] in [1, 3]:  #first or standalone
            self.packet_ids = []
            self.has_first_packet=True
        self.packet_ids.append(packet_id)

        if packet['seg_flag'] in [2, 3]:
            #if self.last_unique_id!=unique_id or self.last_request_spid != packet['SPID']:
            #the last or standalone
            if self.db:
                report = {
                    'start_unix_time': stix_datetime.scet2unix(start),
                    'start_scet': start,
                    'unique_id': unique_id,
                    'is_complete':self.has_first_packet,
                    'packet_ids': self.packet_ids,
                    'run_id': run_id,
                    'SPID': packet['SPID'],
                    'name': DATA_REQUEST_REPORT_NAME[packet['SPID']],
                    'header_unix_time': packet['unix_time'],
                    'header_scet': packet['SCET'],
                    '_id': self.current_id,
                }
                try:
                    req_form = self.db_bsd_forms.find_one(
                        {'unique_ids': int(unique_id)})
                    if req_form:
                        report['request_form'] = req_form
                except Exception as e:
                    print(e)
                    pass

                self.db.insert_one(report)
                self.current_id += 1
                self.has_first_packet=False
                self.packet_ids = []
                


        self.last_unique_id = unique_id
        self.last_request_spid = packet['SPID']

    def process_aspect(self, run_id, packet_id, packet):
        start_obt = packet[1].raw + packet[2].raw / 65536
        summing = packet[3].raw
        samples = packet[4].raw
        duration = samples / 64. * summing

        if packet['seg_flag'] in [1, 3]:
            self.packet_ids = []
            self.start_obt_time = start_obt
        self.packet_ids.append(packet_id)

        end_time = start_obt + duration

        if packet['seg_flag'] in [2, 3]:
            #the last or standalone
            if self.db:
                report = {
                    'start_unix_time':
                    stix_datetime.scet2unix(self.start_obt_time),
                    'end_unix_time': stix_datetime.scet2unix(end_time),
                    'start_scet': self.start_obt_time,
                    'end_scet': end_time,
                    'packet_ids': self.packet_ids,
                    'SPID': packet['SPID'],
                    'run_id': run_id,
                    'name': DATA_REQUEST_REPORT_NAME[packet['SPID']],
                    'header_unix_time': packet['unix_time'],
                    'header_scet': packet.get('SCET', 0),
                    '_id': self.current_id,
                }
                self.db.insert_one(report)
                self.current_id += 1
