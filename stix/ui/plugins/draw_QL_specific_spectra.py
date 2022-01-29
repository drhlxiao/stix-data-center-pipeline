# -*- encoding: utf-8 -*-
""" A Qt GUI to plot STIX quicklook specific spectra

In order to use it, you need to select a background report in the 
GUI and execute this plugin in the plugin manager
"""
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
from stix.core import datatypes as sdt

SPID = 54120


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
        group_data=packet[14]['children']
        num_samples=len(group_data)
        print('num samples:')
        num_groups=int(num_samples/35)
        print(num_samples, num_groups)
        #if int(num_samples)%35 !=0:
        #    print('Unrecognized format')
        #    return
        spectra={}
        detector_channels=[]
        triggers=[]
        for i in range(0, num_groups):
            detector=group_data[i*35][1]
            spectrum_parameters=[x[2] for x in group_data[i*35+1:i*35+32]] #eng value
            spectra[detector]=spectrum_parameters
            detector_channels.append(detector)
            triggers.append(group_data[i*35+33][2])




        fig = plt.figure(figsize=figsize)
        plt.subplot(211)
        for detector, spectrum in spectra.items():
            plt.plot(spectrum, label='detector {}'.format(detector))
            #plt.plot(lc, label='LC {}'.format(energy_bins[ilc]))
        plt.title('QL specific spectrum ( {} )'.format(packet['header']['UTC']))
        plt.xlabel('Energy channel')
        plt.ylabel('Counts in {} s '.format(int_duration * 0.1))
        plt.legend()

        plt.subplot(212)
        plt.plot(detector_channels,triggers)
        plt.xlabel('Detector #')
        plt.ylabel('Counts in {} s '.format(int_duration * 0.1))
        plt.tight_layout()
        plt.show()
