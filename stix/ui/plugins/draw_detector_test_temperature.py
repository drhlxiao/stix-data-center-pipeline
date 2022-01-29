# plugin to draw current
import os

import sys
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
from stix.core import datatypes as sdt

SPID = 54131


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.iql = 0

    def run(self):
        pkt = self.packets[self.current_row]
        packet = sdt.Packet(pkt)
        if not packet.isa(SPID):
            print('Not a temperature test report')
            return
        fig = None

        UTC = packet['header']['UTC']

        detector_number = packet.get('NIX00103/NIX00100')
        temp_mean = packet.get('NIX00103/NIX00101.eng')
        temp_stdev = packet.get('NIX00103/NIX00102.eng')

        plt.errorbar(detector_number[0], temp_mean[0], yerr=temp_stdev[0])
        plt.title('ASIC temperature at {}'.format(UTC))
        plt.xlabel('Detector #')
        plt.ylabel('Temperature (deg)')
        plt.show()
