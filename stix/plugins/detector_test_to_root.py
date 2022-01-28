# -*- encoding: utf-8 -*-
""" A  Qt parser GUI plugin used to plot ASIC on-demand readouts, threshold readouts and 
    dark current readouts. 
In order to use it, you need to select a background report in the 
GUI and execute this plugin in the plugin manager
"""

import os
import sys
from ROOT import TTree, TFile

import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
from stix.core import datatypes as sdt
from tkinter import filedialog
from tkinter import *

SPIDs = [54132, 54133, 54134, 54130]
#On demand readout or baseline

DESCR = {
    54132: 'On-demand readout',
    54133: 'Baseline readouts',
    54134: 'Dark current readouts',
    54130: 'threshold scan'
}
PIXEL_TO_ASIC = [26, 15, 8, 1, 29, 18, 5, 0, 30, 21, 11, 3]


def get_channel(ch):
    try:
        return PIXEL_TO_ASIC[ch]
    except IndexError:
        return -1


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
        self.iql = 0

    def run(self):
        isub = 0
        fig = None
        fname = get_filename()

        f = TFile(fname, 'recreate')
        t = TTree('stix', 'detector test')

        coarse_time = np.empty((1), dtype=np.uint32)
        packet_id= np.empty((1), dtype=np.uint32)
        spid = np.empty((1), dtype='int32')
        detector = np.empty((1), dtype='int32')
        channel = np.empty((1), dtype='int32')
        ADC_mean = np.empty((1), dtype='int32')
        ADC_stdev = np.empty((1), dtype='int32')
        desc = np.empty((32), dtype='byte')

        #des= array('c', 32*[''])
        run = np.empty((1), dtype='int32')

        t.Branch('spid', spid, 'spid/I')
        t.Branch('detector', detector, 'detector/I')
        t.Branch('channel', channel, 'channel/I')
        t.Branch('mean', ADC_mean, 'mean/I')
        t.Branch('std', ADC_stdev, 'std/I')
        t.Branch('run', run, 'run/I')
        t.Branch('desc', desc, 'desc[32]/C')
        t.Branch('coarse_time', coarse_time, 'coarse_time/I')
        t.Branch('packet_id', packet_id, 'packet_id/I')

        is_thr = False
        num_read = 0
        current_run = 0
        packet_id[0]=0
        for pkt in self.packets:
            packet = sdt.Packet(pkt)
            if packet.SPID not in SPIDs and num_read == 0:
                continue
            num_read += 1
            if packet.SPID not in SPIDs:
                print('Not ADC readout report')
                continue
            if packet.SPID == 54130:
                is_thr = True
            desc = DESCR[packet.SPID]

            print('Analyzing packet # {}'.format(self.current_row))
            coarse_time[0]= packet.coarse_time
            packet_id[0]+=1
            det = packet.get('NIX00104/NIX00100')[0]
            chan = packet.get('NIX00104/NIX00105')[0]
            adc = []
            std = []
            selector = 'NIX00104/NIX00106'
            if is_thr:
                selector = 'NIX00104/NIX00108'
            else:
                std = packet.get('NIX00104/NIX00107')[0]
            adc = packet.get(selector)[0]
            for k, d in enumerate(det):
                spid[0] = packet.SPID
                detector[0] = d
                channel[0] = chan[k]
                if get_channel(chan[k]) == -1:
                    continue

                ADC_mean[0] = adc[k]
                if not is_thr:
                    ADC_stdev[0] = std[k]
                run[0] = current_run
                t.Fill()

            if packet.seg_flag in [2, 3]:
                current_run += 1
            self.current_row += 1

        print('Number of packets:', num_read)
        t.Write()
        f.Close()
