# -*- encoding: utf-8 -*-
""" A Qt GUI to plot STIX quicklook background reports

In order to use it, you need to select a quicklook background report in the 
GUI and execute this plugin in the plugin manager
"""
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
from stix.core import stix_datatypes as sdt

SPID = 54119


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

        compression_s = packet[5].raw
        compression_k = packet[6].raw
        compression_m = packet[7].raw

        light_curve = packet.get('NIX00270/NIX00277/*.eng')[0]
        triggers = packet.get('NIX00273/*.eng')

        fig = plt.figure(figsize=figsize)

        plt.subplot(211)

        #energy_bins=['0 - 10 keV', '10 - 15 keV' , '15 - 25 keV','25 - 50 keV' , '50 - 150 keV']
        for ilc, lc in enumerate(light_curve):
            print("Length of light curve :")
            print(len(lc))
            plt.plot(lc, label='LC {}'.format(ilc))
            #plt.plot(lc, label='LC {}'.format(energy_bins[ilc]))

        plt.title('QL background ( {} )'.format(packet['header']['UTC']))
        plt.xlabel('Time / {} (s)'.format(int_duration * 0.1))
        plt.ylabel('Counts in {} s '.format(int_duration * 0.1))

        plt.subplot(212)
        for lt in triggers:
            plt.plot(lt, label='triggers')
        plt.title('Triggers')
        plt.xlabel('Time / {} (s)'.format(int_duration * 0.1))
        plt.ylabel('Counts in {} s '.format(int_duration * 0.1))
        plt.tight_layout()
        plt.show()
