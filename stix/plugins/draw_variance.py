# plugin to draw current
import os

import sys
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt
from stix.core import stix_datatypes as sdt
from stix.spice import stix_datetime

SPID = 54121


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.iql = 0

    def run(self):
        figsize = (12, 8)
        isub = 0
        T0 = 0
        variances=[]
        obts=[]
        sample_per_variance=0

        for pkt in self.packets:
            packet = sdt.Packet(pkt)
            if not packet.isa(SPID):
                continue
            seg_flag = packet['seg_flag']
            start_obt= packet[1].raw+packet[2].raw/65535
            sample_per_variance=packet[4].raw
            integration_duration=(packet[3].raw+1)*0.1

            var= packet.get('NIX00280/*.eng')[0]
            t=[ i*integration_duration+start_obt for i in range(0,len(var))]
            obts.extend(t)
            variances.extend(var)
        

        fig = plt.figure(figsize=figsize)
        plt.plot(obts,variances)
        plt.title('STIX Variance')
        plt.xlabel('SCET (s)')
        plt.ylabel('Variance of {} samples'.format(sample_per_variance))
        plt.legend()
        plt.show()

