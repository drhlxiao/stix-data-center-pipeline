import os
import sys

import numpy as np
from datetime import datetime
from stix.core import stix_datatypes as sdt

SPID = 54118


class LightCurveAnalyzer(object):
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.iql = 0
        self.h2counter = []

    def run(self, pdf):
        print('Number of packets : {}'.format(len(self.packets)))
        #with PdfPages(filename) as pdf:
        figsize = (12, 8)
        isub = 0
        spectra_container = []
        T0 = 0
        for pkt in self.packets:
            packet = sdt.Packet(pkt)
            if not packet.isa(SPID):
                continue
            seg_flag = packet['seg_flag']
            if seg_flag in (1, 3):
                #first packet
                self.h2counter.clear()

            fig = None

            scet_coarse = packet[1].raw
            scet_fine = packet[2].raw
            int_duration = packet[3].raw + 1

            detector_mask = packet[4].raw
            pixel_mask = packet[6].raw

            num_lc = packet[17].raw

            compression_s = packet[8].raw
            compression_k = packet[9].raw
            compression_m = packet[10].raw

            num_lc_points = packet.get('NIX00270/NIX00271')[0]

            light_curve = packet.get('NIX00270/NIX00271/*.eng')[0]
            triggers = packet.get('NIX00273/*.eng')
            rcr = packet.get('NIX00275/*.raw')

            UTC = packet['header']['UTC']

            fig = plt.figure(figsize=figsize)
            plt.axis('off')
            title = ' QL sum LC: # {} \n Packets received at: {} \n T0: {} \
            \n SCET: {} \n Comp_S: {} \n Comp_K: {} \n Comp_M:{}'.format(
                self.iql, UTC, int_duration * 0.1,
                scet_coarse + scet_fine / 65536., compression_s, compression_k,
                compression_m)
            plt.text(0.5, 0.5, title, ha='center', va='center')
            pdf.savefig()
            plt.close()
            #fig.clf()

            fig = plt.figure(figsize=figsize)
            ax = plt.subplot(111)
            for ilc, lc in enumerate(light_curve):
                print("Length of light curve :")
                print(len(lc))
                ax.plot(lc, label='LC {}'.format(ilc))

            ax.legend()
            ax.set_title('QL lightcurves')
            ax.set_xlabel('Time / {} (s)'.format(int_duration * 0.1))
            ax.set_ylabel('Counts in {} s '.format(int_duration * 0.1))
            pdf.savefig()
            plt.close()
            fig = plt.figure(figsize=figsize)
            ax = plt.subplot(111)
            for lt in triggers:
                print("Length of trigger counter accumulator:")
                print(len(lt))
                ax.plot(lt, label='triggers')
                ax.set_title('Triggers')
                ax.set_xlabel('Time / {} (s)'.format(int_duration * 0.1))
                ax.set_ylabel('Counts in {} s '.format(int_duration * 0.1))
            pdf.savefig()
            plt.close()

            fig = plt.figure(figsize=figsize)
            ax = plt.subplot(111)
            for lt in rcr:
                print("Length of RCR:")
                print(len(lt))
                ax.plot(lt, label='RCR')
                ax.set_title('RCR')
                ax.set_xlabel('Time / {} (s)'.format(int_duration * 0.1))
                ax.set_ylabel('RCR ')
            pdf.savefig()
            plt.close()
            self.iql += 1
