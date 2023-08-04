import os
import sys
import numpy as np
from stix.spice import time_utils
from stix.core import datatypes as sdt
from stix.core.energy_bins import StixEnergyBins

QLLC_SPID = 54118
QLBKG_SPID= 54120


class QLAnalyzer(object):
    def __init__(self):
        pass
    @classmethod
    def parse(cls, packets):
        pass
class LightCurveAnalyzer(QLAnalyzer):

    def __init__(self):
        self.data=None
    
    @staticmethod
    def compute_ql_counts_with_BKG(self):
        """
        compute QL curves using BKG light curves
        """
        pass





    
    @classmethod
    def parse(cls, packets, exclude_duplicated=True):
        if not packets:
            return None
        lightcurves = {}
        unix_time = []
        energy_bins={}
        rcrs=[]
        last_time=0
        for pkt in packets:
            packet = sdt.Packet(pkt)
            if not packet.isa(QLLC_SPID):
                continue
            #fig = None

            scet_coarse = packet[1].raw
            scet_fine = packet[2].raw
            start_scet = scet_coarse + scet_fine / 65536.

            if start_scet<=last_time and exclude_duplicated:
                continue
            last_time=start_scet
            int_duration = (packet[3].raw + 1) * 0.1

            detector_mask = packet[4].raw
            pixel_mask = packet[6].raw

            num_lc = packet[17].raw

            compression_s = packet[8].raw
            compression_k = packet[9].raw
            compression_m = packet[10].raw

            num_lc_points = packet.get('NIX00270/NIX00271')[0]
            lc = packet.get('NIX00270/NIX00271/*.eng')[0]
            rcr = packet.get('NIX00275/*.raw')
            UTC = packet['header']['UTC']
            rcrs.extend(rcr[0])
            for i in range(len(lc)):
                if i not in lightcurves:
                    lightcurves[i] = []
                lightcurves[i].extend(lc[i])
            unix_time.extend([
                time_utils.scet2unix(start_scet + x * int_duration)
                for x in range(num_lc_points[0])
            ])

            if not energy_bins:
                energy_bin_mask= packet[16].raw
                tstart=time_utils.scet2unix(start_scet)
                energy_bins=StixEnergyBins.get_emask_energy_bins(energy_bin_mask, tstart)

        if not lightcurves:
            return None
        return {'time': np.array(unix_time), 'lcs': {x: np.array(lightcurves[x]) for x in lightcurves}, 'rcr':np.array(rcrs),
                'energy_bins':energy_bins,'num':len(unix_time),'start_unix': unix_time[0],'end_unix':unix_time[-1]}

class BackgroundReportAnalyzer(QLAnalyzer):
    @classmethod
    def parse(cls, packets):
        if not packets:
            return None
        results=[]
        for pkt in packets:
            packet = sdt.Packet(pkt)
            if not packet.isa(QLBKG_SPID):
                continue
            scet_coarse = packet[1].raw
            scet_fine = packet[2].raw
            start_unix=time_utils.scet2unix(scet_coarse, scet_fine)
            tbin= (packet[3].raw + 1)*0.1
            num_samples=packet[14].raw
            samples=packet[14].children
            for i in range(num_samples):
                detector=samples[0+35*i][1]
                spectra=[samples[j+35*i][2] for j in range(32)]
                trig=samples[33+35*i][2]
                t=start_unix+samples[34+35*i][1]*tbin#number of integrations after the first one
                struct={
                        'detector':detector,
                        'spectra':spectra,
                        'triggers':trig,
                        'obs_time': t,
                        'tbin':tbin,
                        'packet_id':pkt['_id']
                        }
                results.push(struct)
        return results


