# -*- encoding: utf-8 -*-
""" A Qt GuI plugin to plot calibration spectra 

In order to use it, you need to select a calibration report in the 
GUI and execute this plugin in the plugin manager
"""
import os
import sys
import numpy as np

from datetime import datetime
from stix.core import datatypes as sdt
SPID = 54124
import detector_geometry as dg


def plot(h2counter):
    data=np.zeros((32,12))
    for e in h2counter:
        data[e[0]][e[1]]=e[2]
    print(data)
    dg.create_STIX_svg(data.tolist(),'./calibration.svg')


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.ical = 0
        self.h2counter = []
        self.spectra_container = {}

    def run(self):
        print('Number of packets : {}'.format(len(self.packets)))
        num_packets = len(self.packets)
        num_read = 0
        while self.current_row < num_packets:
            pkt = self.packets[self.current_row]
            packet = sdt.Packet(pkt)
            self.current_row += 1

            if not packet.isa(SPID):
                continue
            seg_flag = packet['seg_flag']

            detector_ids = packet.get('NIX00159/NIXD0155')[0]
            pixels_ids = packet.get('NIX00159/NIXD0156')[0]
            spectra = packet.get('NIX00159/NIX00146/*')[0]

            live_time = packet[4].raw
            quiet_time = packet[3].raw
            num_read += 1

            for i, spec in enumerate(spectra):
                if sum(spec) > 0:
                    det = detector_ids[i]
                    pixel = pixels_ids[i]
                    self.spectra_container[(det, pixel)] = spec
                    self.h2counter.append((det, pixel, sum(spec)))

            if seg_flag in [2, 3]:
                break
        num = len(self.spectra_container)
        if num == 0:
            print('spectra empty')
            return
        print('Number of packets read:', num_read)
        plot(self.h2counter)
