# plugin to draw current
import os

import sys
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
from stix.core import datatypes as sdt

SPID = 54118


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.iql = 0

    def run(self):
        figsize = (12, 8)
        isub = 0
        spectra_container = []
        T0 = 0
        pkt = self.packets[self.current_row]
        packet = sdt.Packet(pkt)
        if not packet.isa(SPID):
            print('Not a background packet')
            return
        seg_flag = packet['seg_flag']
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

        light_curve = packet.get('NIX00270/NIX00271/*.eng')[0]
        triggers = packet.get('NIX00273/*.eng')
        rcr = packet.get('NIX00275/*.raw')

        fig = plt.figure(figsize=figsize)

        plt.subplot(311)
        for ilc, lc in enumerate(light_curve):
            print("Length of light curve :")
            print(len(lc))
            plt.plot(lc, label='LC {}'.format(ilc))
            #plt.plot(lc, label='LC {}'.format(energy_bins[ilc]))

        plt.title('QL lightcurve ( {} )'.format(packet['header']['UTC']))
        plt.xlabel('Time / {} (s)'.format(int_duration * 0.1))
        plt.ylabel('Counts in {} s '.format(int_duration * 0.1))
        plt.legend()

        plt.subplot(312)
        for lt in triggers:
            print("Length of trigger counter accumulator:")
            print(len(lt))
            plt.plot(lt, label='triggers')
        plt.title('Triggers')
        plt.xlabel('Time / {} (s)'.format(int_duration * 0.1))
        plt.ylabel('Counts in {} s '.format(int_duration * 0.1))
        plt.subplot(313)
        for lt in rcr:
            plt.plot(lt, label='RCR')
        plt.title('RCR')
        plt.xlabel('Time / {} (s)'.format(int_duration * 0.1))
        plt.ylabel('RCR ')
        plt.tight_layout()
        plt.show()
