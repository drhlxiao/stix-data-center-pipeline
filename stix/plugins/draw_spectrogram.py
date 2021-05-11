# plugin to draw current
import os

import sys
sys.path.append(os.path.abspath('../../'))
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt
from stix.core import stix_datatypes as sdt
from stix.spice import stix_datetime

SPID = 54143


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
        delta_time=[]
        spectra=[]
        num_spectra=[]

        for pkt in self.packets:
            packet = sdt.Packet(pkt)
            if not packet.isa(SPID):
                continue
            seg_flag = packet['seg_flag']
            start_obt= packet[12].raw
            t0= packet.get('NIX00403/NIX00089/NIX00441')
            spec= packet.get('NIX00403/NIX00089/NIX00270/*')[0][0]
            num= packet.get('NIX00403/NIX00089/NIX00270')[0][0]
            t=[ 0.1* x+start_obt for x in t0[0][0]]
            delta_time.extend(t)
            spectra.extend(spec)
            num_spectra.extend(num)
        

        fig = plt.figure(figsize=figsize)
        lc={}
        print(num_spectra[0])
        print(spectra[0])

        for k in range(0,num_spectra[0]):

            lc[k]=[ x[k] for x in spectra]
            plt.plot(delta_time, lc[k], label='LC -  E{}'.format(k))
        plt.title('STIX spectrogram')
        plt.xlabel('SCET (s)')
        plt.ylabel('Counts in 20 seconds')
        plt.legend()
        plt.show()
        plt.show()
            
