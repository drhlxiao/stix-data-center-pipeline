import os

import sys

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from datetime import datetime
from core import packet_analyzer as sta

analyzer = sta.analyzer()

SPIDs = [54131]


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.ical = 0
        self.h2counter = []

    def run(self, filename):
        print('Number of packets : {}'.format(len(self.packets)))
        with PdfPages(filename) as pdf:
            figsize = (12, 8)

            isub = 0
            spectra_container = []
            T0 = 0
            for packet in self.packets:
                try:
                    if int(packet['header']['SPID']) in SPIDs:
                        continue
                except ValueError:
                    continue
                header = packet['header']
                fig = None
                analyzer.load(packet)
                detector_ids = analyzer.to_array('NIX00103/NIX00100')[0]
                asic_tmean = analyzer.to_array('NIX00103/NIX00101')[0]
                asic_tstd = analyzer.to_array('NIX00103/NIX00102')[0]
