# -*- encoding: utf-8 -*-
""" A Qt GuI plugin to plot calibration spectra 

In order to use it, you need to select a calibration report in the 
GUI and execute this plugin in the plugin manager
"""
import os
import sys
import numpy as np

from datetime import datetime
from tkinter import filedialog
from tkinter import *
from ROOT import TTree, TFile, TH1F, TH2F

from stix.core import datatypes as sdt
SPID = 54124




def get_filename():
    root = Tk()
    root.filename = filedialog.asksaveasfilename(
        initialdir=".",
        title="Set output filename",
        filetypes=(("ROOT", "*.root"), ("all files", "*.*")))
    filename = root.filename
    root.withdraw()
    return filename
class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.ical = 0
        self.h2counter = []
        self.root_file=None

    def run(self):
        print('Number of packets : {}'.format(len(self.packets)))
        num_packets = len(self.packets)
        num_read = 0

        root_filename=get_filename()
        if not root_filename:
            print("ROOT filename not set")
            return

        file_root=TFile(root_filename,'recreate')

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
                    print(' {}-{} has {} counts'.format(det,pixel, sum(spec)))
                    hspc=TH1F("hcal_{}_{}".format(det,pixel), 
                            "Calibration spectrum of det {}, pixel {}; Energy channel; Counts".format(det,pixel),
                            len(spec), 0, len(spec))
                    for k, count in enumerate(spec):
                        hspc.SetBinContent(k+1,count)
                    self.h2counter.append((det, pixel, sum(spec)))
                    hspc.Write('hcal_{}_{}'.format(det,pixel))
            if seg_flag in [2, 3]:
                break

        file_root.cd()
        h2=TH2F('h2',"hits; Detector; Pixel; Counts",32, 0,32,12, 0,12)
        for e in self.h2counter:
            h2.Fill(e[0]+0.5, e[1]+0.5, e[2])
        h2.Write()
        file_root.Close()
        print('Output:',root_filename)
        
